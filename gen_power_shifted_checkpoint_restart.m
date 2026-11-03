function [X, D, iters, relres, errorIterations, checkpointTime, recoveryTime, e2eTime, stats] = ...
    gen_power_shifted_checkpoint_restart(A, X0, k, possibleErrorIterations, rowDeletedStart, varargin)
% The same shifted block power iteration as gen_erasure_power_shifted using checkpoint/restart
% instead of erasure coding (baseline for the comparison). 

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'tol',         1e-8);
    addParameter(p, 'maxiter',     2000);
    addParameter(p, 'blockFactor', 2);
    addParameter(p, 'shiftMode',   'auto');
    addParameter(p, 'shiftFrac',   0.5);
    addParameter(p, 'padSeed',     42);
    for ign = {'coding', 'codingScheme'}
        addParameter(p, ign{1}, []);
    end
    parse(p, varargin{:});
    o  = p.Results;
    um = p.Unmatched;  crArgs = reshape([fieldnames(um) struct2cell(um)]', 1, []);
    cr = cr_parse_options(crArgs{:});
    if ~ismember(o.shiftMode, {'auto', 'none'})
        error('gen_power_shifted_checkpoint_restart:badShiftMode', 'shiftMode must be ''auto'' or ''none''.');
    end

    [n, s] = size(X0);
    tol = o.tol;  maxiter = o.maxiter;
    m = min(round(o.blockFactor * s), n);
    if size(X0, 2) < m
        sc = norm(X0, 'fro') / max(sqrt(numel(X0)), 1);
        if ~isfinite(sc) || sc <= 0, sc = 1; end
        X0 = [X0, sc * randn(RandStream('mt19937ar', 'Seed', o.padSeed), n, m - size(X0, 2))];
    elseif size(X0, 2) > m
        X0 = X0(:, 1:m);
    end
    normA  = norm(A, 'fro');
    faults = sort(unique(possibleErrorIterations(possibleErrorIterations > 0)));
    faults = faults(:).';
    if cr.interval_mode == "none" && ~isempty(faults)
        warning('gen_power_shifted_checkpoint_restart:unprotected', ...
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
    X = X0;  sigma = 0;  r0 = 1;  D = [];  AXprev = [];
    relres = [];  absres = [];  iters = [];  errorIterations = [];  sigHist = [];  iterWall = [];
    checkpointTime = 0;  recoveryTime = 0;  nCkpt = 0;  ckptIters = [];  ckptCost = [];
    loadATime = 0;  excludedLoadA = 0;  loadZTime = 0;
    execLogical = [];  execDur = [];  wasted = false(0);
    firedFault = false(size(faults));
    converged = false;

    % checkpoint interval T: static_interval in "static" mode 
    T = Inf;
    if cr.interval_mode == "static"
        T = min(cr.static_interval, max(maxiter, 1));
        [ctx, tC] = cr_store('write', ctx, packCkpt(0));
        checkpointTime = checkpointTime + tC;  nCkpt = 1;  ckptIters = 0;  ckptCost = tC;
        if cr.cr_verbose, fprintf('[CR] checkpoint cost C = %.4f s (iteration-0 checkpoint)\n', tC); end
    end

    fprintf('block power (shifted) with %s fault tolerance: m = %d, s = %d\n', cr.interval_mode, m, s);
    fprintf('iter\tResidual Norm\tRelative Residual Norm\n');
    tLoop0 = tic;
    j = 0;  it = 0;
    itMax = maxiter * (numel(faults) + 1) + numel(faults);
    while it < itMax
        it = it + 1;

        % fault: reload A, restart or roll back (timed as recoveryTime)
        hit = find(faults == it & ~firedFault, 1);
        if ~isempty(hit)
            firedFault(hit) = true;
            errorIterations(end+1) = it;
            fprintf('\n[CR] fault before executed iteration %d (node rows %d..%d): A and iterate lost\n', ...
                    it, rowDeletedStart, rowDeletedStart + k - 1);
            tr = tic;
            A = [];  X = [];  AXprev = [];
            [A, tLA] = cr_store('readA', ctx);  loadATime = loadATime + tLA;
            if cr.interval_mode == "restart"
                c = 0;  X = X0;  sigma = 0;  r0 = 1;  D = [];
            else
                [P, tLZ] = cr_store('readZ', ctx);  loadZTime = loadZTime + tLZ;
                c = P.j;  X = P.X;  sigma = P.sigma;  r0 = P.r0;  D = P.D;
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

        % one shifted power step with Rayleigh-Ritz
        j  = j + 1;
        ti = tic;
        if isempty(AXprev), Y = A * X; else, Y = AXprev; end
        Z = Y;
        if sigma ~= 0, Z = Z - sigma * X; end
        sigHist(end+1) = sigma;
        [Q, ~] = qr(Z, 0);
        AQ_mat = A * Q;
        AQ = Q' * AQ_mat;   BQ = Q' * Q;
        AQ = (AQ + AQ') / 2;  BQ = (BQ + BQ') / 2;
        [U, Dmat] = eig(AQ, BQ, 'chol');
        [D, idx] = sort(real(diag(Dmat)), 'descend');
        U = U(:, idx);
        X = Q * U;
        AX = AQ_mat * U;
        if strcmp(o.shiftMode, 'auto') && D(end) > 0
            sigma = o.shiftFrac * D(end);
        else
            sigma = 0;
        end
        BX = Q * U;
        AXprev = AX;
        if j == 1
            r0 = norm(AX(:, 1:s) - BX(:, 1:s) .* D(1:s)', 'fro');
        end
        [resNorm, relresNorm] = compute_residual(AX(:, 1:s), BX(:, 1:s), diag(D(1:s)), r0, normA);
        relres(end+1) = relresNorm;  absres(end+1) = resNorm;  iters(end+1) = it;
        execLogical(end+1) = j;  execDur(end+1) = toc(ti);  wasted(end+1) = false;
        iterWall(end+1) = toc(tLoop0);
        fprintf('%d\t%e\t%e\t%e\n', it, resNorm, relresNorm, tol);
        if relresNorm < tol
            converged = true;
            fprintf('Converged in %d logical iterations (%d executed, rel. residual = %e)\n', j, it, relresNorm);
            break;
        end
        if j >= maxiter, break; end

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
    stats.method          = 'block_power_shifted_cr';
    stats.sigma           = 'largestabs';
    stats.shiftMode       = o.shiftMode;
    stats.sigmaHist       = sigHist;
    stats.iters           = iters;
    stats.res             = relres;
    stats.absres          = absres;
    stats.eigenvalues     = D;
    stats.converged       = converged;
    stats.n_iters         = numel(iters);
    stats.logicalIters    = j;
    stats.faultIterations = errorIterations;
    stats.r0              = r0;
    stats.blockWidth      = m;
    stats.iterWall        = iterWall;
    stats.setupTime       = 0;
    stats.cr = struct('mode', cr.interval_mode, 'interval', T, ...
        'nCkpt', nCkpt, 'ckptIters', ckptIters, 'ckptCost', ckptCost, ...
        'checkpointTime', checkpointTime, 'recoveryTime', recoveryTime, ...
        'loadATime', loadATime, 'chargeLoadA', cr.charge_load_A, 'excludedLoadATime', excludedLoadA, 'loadZTime', loadZTime, 'rebuildTime', 0, ...
        'reworkTime', reworkTime, 'lostIters', nnz(wasted), 'nFaults', numel(errorIterations), ...
        'saveATime', saveATime, 'chargeSaveA', cr.charge_save_A, ...
        'bytesA', getf(ctx, 'bytesA'), 'bytesZ', getf(ctx, 'bytesZ'), ...
        'execLogical', execLogical, 'execDur', execDur, 'wasted', wasted);
    stats.timing = struct('setupTime', 0, 'checkpointTime', checkpointTime, ...
        'recoveryTime', recoveryTime, 'reworkTime', reworkTime, ...
        'solvTime', max(e2eTime - checkpointTime - recoveryTime, 0), 'e2eTime', e2eTime);

    fprintf('\n[CR-SUMMARY] mode %s | interval %g | checkpoints %d | faults %d | lost iters %d\n', ...
            cr.interval_mode, T, nCkpt, numel(errorIterations), nnz(wasted));
    fprintf('e2eTime(s)\tcheckpointTime(s)\trecoveryTime(s)\treworkTime(s)\n');
    fprintf('%f\t%f\t%f\t%f\n', e2eTime, checkpointTime, recoveryTime, reworkTime);

    function P = packCkpt(jj)
        P = struct('j', jj, 'X', X, 'sigma', sigma, 'r0', r0, 'D', D);
    end
end

function v = getf(S, f)
    if isstruct(S) && isfield(S, f), v = S.(f); else, v = 0; end
end
