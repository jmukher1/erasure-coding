function [X, D, iters, relres, errorIterations, checkpointTime, recoveryTime, e2eTime, stats] = ...
    inverse_power_shifted_checkpoint_restart(A, X0, k, possibleErrorIterations, rowDeletedStart, varargin)
% The same shifted block inverse power iteration as inverse_power_shifted, protected by checkpoint/restart
% instead of erasure coding (baseline for the comparison).
%
% A fault at executed iteration it loses A, its factorization and the
% iterate. Recovery reloads A from stable storage and either restarts from
% X0 ("restart") or reloads the last checkpoint and refactors ("static",
% checkpoint every static_interval iterations). "none" ignores faults
% (unprotected reference run). Solver options as in inverse_power_shifted ('coding' and
% 'codingScheme' are accepted and ignored); every other name/value pair goes
% to cr_parse_options.
%
% Timing (wall clock, seconds, tic/toc)
%   setupTime       initial factorization of A - sigma*I
%   checkpointTime  writing the checkpoints (plus writing A if charge_save_A)
%   recoveryTime    reload of A (only if charge_load_A) + reload of the
%                   checkpoint + refactorization after each fault
%   reworkTime      execution time of the iterations that were rolled back
%   e2eTime         whole call (minus the reload of A when it is not charged)

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'tol',          1e-8);
    addParameter(p, 'maxiter',      1000);
    addParameter(p, 'sigma',        0);
    addParameter(p, 'blockFactor',  2);
    addParameter(p, 'shiftMode',    'auto');
    addParameter(p, 'shiftTol',     1e-2);
    addParameter(p, 'shiftGap',     1e-8);
    addParameter(p, 'rrMode',       'gram');
    addParameter(p, 'padSeed',      42);
    addParameter(p, 'shiftSafety',  1.1);
    addParameter(p, 'shiftFloor',   1e-8);
    addParameter(p, 'refactorGain', 0.9);
    addParameter(p, 'maxRefactor',  30);
    addParameter(p, 'stagnWindow',  20);
    addParameter(p, 'stagnFactor',  0.9);
    for ign = {'coding', 'codingScheme'}
        addParameter(p, ign{1}, []);
    end
    parse(p, varargin{:});
    o  = p.Results;
    um = p.Unmatched;  crArgs = reshape([fieldnames(um) struct2cell(um)]', 1, []);
    cr = cr_parse_options(crArgs{:});
    if ~ismember(o.shiftMode, {'auto', 'fixed'})
        error('inverse_power_shifted_checkpoint_restart:badShiftMode', 'shiftMode must be ''auto'' or ''fixed''.');
    end
    if ~ismember(o.rrMode, {'gram', 'qr'})
        error('inverse_power_shifted_checkpoint_restart:badRrMode', 'rrMode must be ''gram'' or ''qr''.');
    end

    tol = o.tol;  maxiter = o.maxiter;
    [n, s] = size(X0);
    m = min(round(o.blockFactor * s), n);
    if size(X0, 2) < m
        sc = norm(X0, 'fro') / max(sqrt(numel(X0)), 1);
        if ~isfinite(sc) || sc <= 0, sc = 1; end
        X0 = [X0, sc * randn(RandStream('mt19937ar', 'Seed', o.padSeed), n, m - size(X0, 2))];
    elseif size(X0, 2) > m
        X0 = X0(:, 1:m);
    end
    normA  = norm(A, 'fro');
    normA1 = norm(A, 1);
    faults = sort(unique(possibleErrorIterations(possibleErrorIterations > 0)));
    faults = faults(:).';
    if cr.interval_mode == "none" && ~isempty(faults)
        warning('inverse_power_shifted_checkpoint_restart:unprotected', ...
                'interval_mode = "none": the %d scheduled fault(s) are ignored.', numel(faults));
        faults = [];
    end

    ctx = [];  saveATime = 0;
    if cr.interval_mode ~= "none"
        ctx = cr_store('init', A, cr);
        saveATime = ctx.saveATime;
    end
    cleanStore = onCleanup(@() cr_store('close', ctx));

    tAll = tic;
    tS = tic;
    [X, ~] = qr(X0, 0);
    sigma0 = o.sigma;
    [solveShift, okShift, sigma] = shiftSolver(A, sigma0, n, []);
    if ~okShift
        error('inverse_power_shifted_checkpoint_restart:indefiniteShift', ...
              'A - %g*I is not positive definite; supply a smaller sigma.', sigma0);
    end
    setupTime = toc(tS);

    r0 = 0;  D = [];  nRefactor = 0;  BXprev = [];  nGramFallback = 0;  stagnatedFlag = false;
    relres = [];  absres = [];  iters = [];  errorIterations = [];  sigHist = [];  iterWall = [];
    relresLog = [];
    checkpointTime = 0;  recoveryTime = 0;  nCkpt = 0;  ckptIters = [];  ckptCost = [];
    loadATime = 0;  excludedLoadA = 0;  loadZTime = 0;  rebuildTime = 0;
    execLogical = [];  execDur = [];  wasted = false(0);
    firedFault = false(size(faults));
    converged = false;

    T = Inf;
    if cr.interval_mode == "static"
        T = min(cr.static_interval, max(maxiter, 1));
        [ctx, tC] = cr_store('write', ctx, packCkpt(0));
        checkpointTime = checkpointTime + tC;  nCkpt = 1;  ckptIters = 0;  ckptCost = tC;
        if cr.cr_verbose, fprintf('[CR] checkpoint cost C = %.4f s (iteration-0 checkpoint)\n', tC); end
    end

    fprintf('inverse power (shifted) with %s fault tolerance: m = %d, s = %d\n', cr.interval_mode, m, s);
    fprintf('iter\tResidual Norm\tRelative Residual Norm\n');
    tLoop0 = tic;
    j = 0;  it = 0;
    itMax = maxiter * (numel(faults) + 1) + numel(faults);
    while it < itMax
        it = it + 1;

        % fault: reload A, restart or roll back and refactor
        %      (timed as recoveryTime)
        hit = find(faults == it & ~firedFault, 1);
        if ~isempty(hit)
            firedFault(hit) = true;
            errorIterations(end+1) = it;
            fprintf('\n[CR] fault before executed iteration %d (node rows %d..%d): A, factor and iterate lost\n', ...
                    it, rowDeletedStart, rowDeletedStart + k - 1);
            tr = tic;
            A = [];  X = [];  solveShift = [];  BXprev = [];
            [A, tLA] = cr_store('readA', ctx);  loadATime = loadATime + tLA;
            if cr.interval_mode == "restart"
                tb = tic;
                [X, ~] = qr(X0, 0);
                [solveShift, ~, sigma] = shiftSolver(A, sigma0, n, []);
                rebuildTime = rebuildTime + toc(tb);
                c = 0;  r0 = 0;  D = [];  nRefactor = 0;  relresLog = [];
            else
                [P, tLZ] = cr_store('readZ', ctx);  loadZTime = loadZTime + tLZ;
                tb = tic;
                [solveShift, okR, sigma] = shiftSolver(A, P.sigma, n, []);
                if ~okR
                    warning('inverse_power_shifted_checkpoint_restart:refactor', ...
                            'Refactor at sigma = %g failed; using the initial shift.', P.sigma);
                    [solveShift, ~, sigma] = shiftSolver(A, sigma0, n, []);
                end
                rebuildTime = rebuildTime + toc(tb);
                c = P.j;  X = P.X;  r0 = P.r0;  D = P.D;  nRefactor = P.nRefactor;  relresLog = P.relresLog;
            end
            recoveryTime = recoveryTime + toc(tr);
            if ~cr.charge_load_A
                recoveryTime  = recoveryTime - tLA;
                excludedLoadA = excludedLoadA + tLA;
            end
            wasted(execLogical > c) = true;
            fprintf('[CR] recovered to logical iteration %d (%d iterations to redo) in %.4f s\n', c, j - c, toc(tr));
            j = c;
        end

        % one shifted inverse power step + Rayleigh-Ritz
        j  = j + 1;
        ti = tic;
        if isempty(BXprev), BX = X; else, BX = BXprev; end
        Z = solveShift(BX);
        gramOK = false;
        if strcmp(o.rrMode, 'gram')
            Q = Z ./ max(vecnorm(Z), realmin);
            AQ_mat = A * Q;  BQ_mat = Q;
            AQ = Q' * AQ_mat;  AQ = (AQ + AQ') / 2;
            BQ = Q' * BQ_mat;  BQ = (BQ + BQ') / 2;
            [~, pG] = chol(BQ);
            if pG == 0 && all(isfinite(AQ(:)))
                [U, Dmat] = eig(AQ, BQ, 'chol');
                gramOK = all(isfinite(U(:)));
            end
            if ~gramOK, nGramFallback = nGramFallback + 1; end
        end
        if ~gramOK
            [Q, ~] = qr(Z, 0);
            AQ_mat = A * Q;  BQ_mat = Q;
            AQ = Q' * AQ_mat;  BQ = Q' * BQ_mat;
            AQ = (AQ + AQ') / 2;  BQ = (BQ + BQ') / 2;
            [U, Dmat] = eig(AQ, BQ);
        end
        [D, idx] = sort(real(diag(Dmat)), 'ascend');
        U = U(:, idx);
        X  = Q * U;
        AX = AQ_mat * U;
        BX = BQ_mat * U;
        BXprev = BX;
        if j == 1
            r0 = norm(AX(:, 1:s) - BX(:, 1:s) .* D(1:s)', 'fro');
        end
        [resNorm, relresNorm] = compute_residual(AX(:, 1:s), BX(:, 1:s), diag(D(1:s)), r0, normA);
        relres(end+1) = relresNorm;  absres(end+1) = resNorm;  iters(end+1) = it;
        relresLog(end+1) = relresNorm;
        sigHist(end+1) = sigma;

        stop = false;
        sw = o.stagnWindow;
        if isfinite(sw) && numel(relresLog) > sw && j > sw && relresNorm >= tol
            recent = relresLog(end - sw);
            if isfinite(recent) && recent > 0 && relresNorm > o.stagnFactor * recent
                fprintf('\nInverse power stopped at logical iteration %d: residual stalled at %.3e\n', j, relresNorm);
                stagnatedFlag = true;  stop = true;
            end
        end
        fprintf('%d\t%e\t%e\t%e\tsigma=%.6g\n', it, resNorm, relresNorm, tol, sigma);

        if ~stop && strcmp(o.shiftMode, 'auto') && nRefactor < o.maxRefactor && relresNorm >= tol
            th1 = D(1);
            r1  = norm(AX(:, 1) - th1 * BX(:, 1));
            if th1 > 0 && r1 / th1 < o.shiftTol
                cand = min([th1 - o.shiftSafety * r1, (1 - o.shiftGap) * th1, th1 - 100 * eps * normA1]);
                if cand > sigma && (th1 - cand) < o.refactorGain * (th1 - sigma)
                    [solveNew, okNew, sigNew] = shiftSolver(A, cand, n, sigma);
                    if okNew && sigNew > sigma
                        solveShift = solveNew;  sigma = sigNew;  nRefactor = nRefactor + 1;
                        fprintf('  [shift] sigma -> %.10g (theta_1 = %.10g, refactor #%d)\n', sigma, th1, nRefactor);
                    end
                    solveNew = [];
                end
            end
        end
        execLogical(end+1) = j;  execDur(end+1) = toc(ti);  wasted(end+1) = false;
        iterWall(end+1) = toc(tLoop0);
        if stop, break; end
        if relresNorm < tol
            converged = true;
            fprintf('Converged in %d logical iterations (%d executed, rel. residual = %e)\n', j, it, relresNorm);
            break;
        end
        if j >= maxiter, break; end

        % periodic checkpoint (timed as checkpointTime)
        if isfinite(T) && mod(j, T) == 0
            [ctx, tC] = cr_store('write', ctx, packCkpt(j));
            checkpointTime = checkpointTime + tC;
            nCkpt = nCkpt + 1;  ckptIters(end+1) = j;  ckptCost(end+1) = tC;
        end
    end

    sKeep = min(s, size(X, 2));
    X = X(:, 1:sKeep);
    D = D(1:sKeep);
    e2eTime = toc(tAll) - excludedLoadA;
    if cr.charge_save_A
        e2eTime = e2eTime + saveATime;  checkpointTime = checkpointTime + saveATime;
    end
    reworkTime = sum(execDur(wasted));

    stats = struct();
    stats.method          = 'inverse_power_shifted_cr';
    stats.sigma           = 'smallestabs';
    stats.rrMode          = o.rrMode;
    stats.nGramFallback   = nGramFallback;
    stats.iters           = iters;
    stats.res             = relres;
    stats.absres          = absres;
    stats.eigenvalues     = D;
    stats.converged       = converged;
    stats.stagnated       = stagnatedFlag;
    stats.n_iters         = numel(iters);
    stats.logicalIters    = j;
    stats.faultIterations = errorIterations;
    stats.r0              = r0;
    stats.blockWidth      = m;
    stats.nRefactor       = nRefactor;
    stats.shiftMode       = o.shiftMode;
    stats.sigmaHist       = sigHist;
    stats.setupTime       = setupTime;
    stats.iterWall        = iterWall;
    stats.cr = struct('mode', cr.interval_mode, 'interval', T, ...
        'nCkpt', nCkpt, 'ckptIters', ckptIters, 'ckptCost', ckptCost, ...
        'checkpointTime', checkpointTime, 'recoveryTime', recoveryTime, ...
        'loadATime', loadATime, 'chargeLoadA', cr.charge_load_A, 'excludedLoadATime', excludedLoadA, 'loadZTime', loadZTime, 'rebuildTime', rebuildTime, ...
        'reworkTime', reworkTime, 'lostIters', nnz(wasted), 'nFaults', numel(errorIterations), ...
        'saveATime', saveATime, 'chargeSaveA', cr.charge_save_A, ...
        'bytesA', getf(ctx, 'bytesA'), 'bytesZ', getf(ctx, 'bytesZ'), ...
        'execLogical', execLogical, 'execDur', execDur, 'wasted', wasted);
    stats.timing = struct('setupTime', setupTime, 'checkpointTime', checkpointTime, ...
        'recoveryTime', recoveryTime, 'reworkTime', reworkTime, ...
        'solvTime', max(e2eTime - setupTime - checkpointTime - recoveryTime, 0), 'e2eTime', e2eTime);

    fprintf('\n[CR-SUMMARY] mode %s | interval %g | checkpoints %d | faults %d | lost iters %d\n', ...
            cr.interval_mode, T, nCkpt, numel(errorIterations), nnz(wasted));
    fprintf('e2eTime(s)\tcheckpointTime(s)\trecoveryTime(s)\treworkTime(s)\n');
    fprintf('%f\t%f\t%f\t%f\n', e2eTime, checkpointTime, recoveryTime, reworkTime);

    function P = packCkpt(jj)
        P = struct('j', jj, 'X', X, 'sigma', sigma, 'r0', r0, 'D', D, ...
                   'nRefactor', nRefactor, 'relresLog', relresLog);
    end
end

function [solve, ok, sigmaUsed] = shiftSolver(A, sigma, n, sigmaLo)
    if isempty(sigmaLo), maxTry = 1; sigmaLo = sigma; else, maxTry = 20; end
    cand = sigma;
    for attempt = 1:maxTry
        [solve, ok] = tryChol(A, cand, n);
        if ok, sigmaUsed = cand; return; end
        cand = (sigmaLo + cand) / 2;
        if abs(cand - sigmaLo) <= eps(max(abs(sigmaLo), 1)), break; end
    end
    sigmaUsed = sigma;
end

function [solve, ok] = tryChol(A, sigma, n)
    Ashift = A - sigma * speye(n);
    Ashift = (Ashift + Ashift') / 2;
    solve = [];
    if issparse(Ashift)
        [Rc, p, Sc] = chol(Ashift);
    else
        [Rc, p] = chol(Ashift);  Sc = speye(n);
    end
    ok = (p == 0);
    if ~ok, return; end
    solve = @(V) Sc * (Rc \ (Rc' \ (Sc' * V)));
end

function v = getf(S, f)
    if isstruct(S) && isfield(S, f), v = S.(f); else, v = 0; end
end
