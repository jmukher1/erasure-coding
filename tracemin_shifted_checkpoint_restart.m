function [X, Thi, iters, res, relres, errorIterations, checkpointTime, recoveryTime, e2eTime, stats] = ...
    tracemin_shifted_checkpoint_restart(A, X0, s, possibleErrorIterations, SIGMA, TOL, maxiter, opts, sopts, varargin)
% Same shifted TraceMin iteration as tracemin_shifted, using checkpoint/restart instead of erasure
% coding to be used as baseline for the comparison.
%
 
if nargin < 7 || isempty(maxiter), maxiter = 1000; end
if nargin < 8, opts  = []; end
if nargin < 9 || isempty(sopts), sopts = struct(); end
cr = cr_parse_options(varargin{:});

innerSolve = 'minres';
if isstruct(opts) && isfield(opts, 'innerSolve')
    innerSolve = opts.innerSolve;
    opts = rmfield(opts, 'innerSolve');
end
if isempty(opts), opts = tracemin_opts(); else, opts = tracemin_opts(opts); end
opts.innerSolve = innerSolve;
if ~isfield(sopts, 'shiftHold'), sopts.shiftHold = 1;      end

[n, ~] = size(A);
if isempty(X0), X0 = orth(randn(n, s)); end
B      = speye(n);
if ~strcmpi(SIGMA, 'smallestabs')
    error('tracemin_shifted_checkpoint_restart:sigma', 'only SIGMA = ''smallestabs'' is supported.');
end
normA0 = norm(A, 'fro');
faults = sort(unique(possibleErrorIterations(possibleErrorIterations > 0)));
faults = faults(:).';
if cr.interval_mode == "none" && ~isempty(faults)
    warning('tracemin_shifted_checkpoint_restart:unprotected', ...
            'interval_mode = "none": the %d scheduled fault(s) are ignored.', numel(faults));
    faults = [];
end
ITER_MAX = min(maxiter, n^3);

ctx = [];  saveATime = 0;
if cr.interval_mode ~= "none"
    ctx = cr_store('init', A, cr);
    saveATime = ctx.saveATime;
end
cleanStore = onCleanup(@() cr_store('close', ctx));

tAll = tic;
tS = tic;
[Aw, Bw, Z, ops, sigma] = tmSetup(A, B, X0, s, opts);
setupTime = toc(tS);
A = [];

r0 = 0;  state = [];  nRefactor = 0;  lastFault = 0;  resLog = [];
res = [];  relres = [];  iters = [];  errorIterations = [];  iterWall = [];
sigHist = [];  stagnated = false;  bStop = 0;
checkpointTime = 0;  recoveryTime = 0;
nCkpt = 0;  ckptIters = [];  ckptCost = [];
loadATime = 0;  excludedLoadA = 0;  loadZTime = 0;  rebuildTime = 0;
execLogical = [];  execDur = [];  wasted = false(0);
firedFault = false(size(faults));
% checkpoint interval T: static_interval in "static" mode, none otherwise;
% iteration-0 checkpoint
T = Inf;
if cr.interval_mode == "static"
    T = min(cr.static_interval, max(ITER_MAX, 1));
    [ctx, tC] = cr_store('write', ctx, packCkpt(0));
    checkpointTime = checkpointTime + tC;
    nCkpt = nCkpt + 1;  ckptIters(end+1) = 0;  ckptCost(end+1) = tC;
    if cr.cr_verbose, fprintf('[CR] checkpoint cost C = %.4f s (iteration-0 checkpoint)\n', tC); end
end

fprintf('\nTraceMin (shifted) with %s fault tolerance, %s, tol %g\n', cr.interval_mode, SIGMA, TOL);
fprintf('\nIter \t, ResNorm \t, RelResNorm \t, AbsResNorm \t, Tolerance \t, sigma\n');

tLoop0 = tic;
j  = 0;
it = 0;
itMax = ITER_MAX * (numel(faults) + 1) + numel(faults);
while it < itMax
    it = it + 1;

    % fault: reload A, then restart from X0 or roll back to the last
    %      checkpoint (timed as recoveryTime)
    hit = find(faults == it & ~firedFault, 1);
    if ~isempty(hit)
        firedFault(hit) = true;
        errorIterations(end+1) = it;
        fprintf('\n[CR] fault before executed iteration %d (logical %d): in-memory A, factor and iterate lost\n', it, j + 1);
        tr = tic;
        Aw = [];  Bw = [];  ops = [];  Z = [];
        [Ad, tLA] = cr_store('readA', ctx);   loadATime = loadATime + tLA;
        if cr.interval_mode == "restart"
            tb = tic;
            [Aw, Bw, Z, ops, sigma] = tmSetup(Ad, B, X0, s, opts);
            rebuildTime = rebuildTime + toc(tb);
            c = 0;  r0 = 0;  state = [];  nRefactor = 0;  resLog = [];
        else
            [P, tLZ] = cr_store('readZ', ctx);  loadZTime = loadZTime + tLZ;
            tb = tic;
            Aw = Ad;  Bw = B;
            sigma = P.sigma;
            Fn = spdFactor(Aw - sigma * Bw);
            if ~Fn.ok
                warning('tracemin_shifted_checkpoint_restart:refactor', 'Refactor at sigma = %g failed; using sigma = 0.', sigma);
                sigma = 0;  Fn = spdFactor(Aw);
            end
            ops = buildOps(Aw, Bw, Fn, sigma);  Fn = [];
            rebuildTime = rebuildTime + toc(tb);
            c = P.j;  Z = P.Z;  r0 = P.r0;  state = P.state;  nRefactor = P.nRefactor;
            resLog = P.resLog;  lastFault = P.lastFault;
        end
        Ad = [];
        recoveryTime = recoveryTime + toc(tr);
        if ~cr.charge_load_A
            recoveryTime  = recoveryTime - tLA;
            excludedLoadA = excludedLoadA + tLA;
        end
        wasted(execLogical > c) = true;
        fprintf('[CR] recovered to logical iteration %d (%d iterations to redo) in %.4f s\n', c, j - c, toc(tr));
        j = c;
    end

    j  = j + 1;
    ti = tic;
    [Z, bStop, Thi, X, r0, resNorm, relresNorm, ~, state] = ...
        tracemin_single_iter_shifted(Aw, [], Bw, Z, j, s, r0, TOL, [], struct(), ...
                                     normA0, opts, ops, state);
    res(end+1)    = resNorm;
    relres(end+1) = relresNorm;
    iters(end+1)  = it;
    sigHist(end+1) = sigma;
    resLog(end+1) = resNorm;

    if bStop == 1
        execLogical(end+1) = j;  execDur(end+1) = toc(ti);  wasted(end+1) = false;
        iterWall(end+1) = toc(tLoop0);
        fprintf('\nTraceMin (shifted, %s, %s) converged at logical iteration %d (executed %d)\n', ...
                SIGMA, cr.interval_mode, j, it);
        break;
    end

    if ~strcmp(opts.shiftMode, 'none') && nRefactor < opts.maxRefactor && ...
            (lastFault == 0 || j - lastFault >= sopts.shiftHold)
        sigNew = state.sigmaProposed;
        th1    = state.theta1;
        if sigNew > sigma && (th1 - sigNew) < opts.refactorGain * (th1 - sigma)
            for backoff = 1:4
                Fn = spdFactor(Aw - sigNew * Bw);
                if Fn.ok
                    ops = buildOps(Aw, Bw, Fn, sigNew);  Fn = [];
                    sigma = sigNew;
                    nRefactor = nRefactor + 1;
                    if opts.verbose
                        fprintf('    [shift] sigma -> %.10e (theta_1 = %.10e, refactor #%d)\n', sigma, th1, nRefactor);
                    end
                    break;
                end
                Fn = [];
                sigNew = (sigma + sigNew) / 2;
            end
        end
    end
    execLogical(end+1) = j;  execDur(end+1) = toc(ti);  wasted(end+1) = false;
    iterWall(end+1) = toc(tLoop0);

    w = opts.stagnWindow;
    w0 = max(1, lastFault);
    if isfinite(w) && j > 2 * w && (j - lastFault) > 2 * w && ...
            min(resLog(end-w+1:end)) > opts.stagnFactor * min(resLog(w0:end-w))
        stagnated = true;
        fprintf('\nStagnated: no %.0f%% improvement in %d iterations, stopping at %d\n', ...
                100 * (1 - opts.stagnFactor), w, j);
        break;
    end
    if j >= ITER_MAX, break; end

    % periodic checkpoint (timed as checkpointTime)
    if isfinite(T) && mod(j, T) == 0
        [ctx, tC] = cr_store('write', ctx, packCkpt(j));
        checkpointTime = checkpointTime + tC;
        nCkpt = nCkpt + 1;  ckptIters(end+1) = j;  ckptCost(end+1) = tC;
    end
end

sk  = min(s, size(Thi, 1));
Thi = Thi(1:sk, 1:sk);
X   = X(:, 1:sk);
e2eTime = toc(tAll) - excludedLoadA;               % reload of A excluded unless charged
if cr.charge_save_A
    e2eTime = e2eTime + saveATime;
    checkpointTime = checkpointTime + saveATime;
end
reworkTime = sum(execDur(wasted));

stats = struct();
stats.method          = 'tracemin_shifted_cr';
stats.sigma           = SIGMA;
stats.iters           = iters;
stats.res             = relres;
stats.relres          = relres;
stats.absres          = res;
stats.eigenvalues     = diag(Thi);
stats.converged       = bStop == 1;
stats.stagnated       = stagnated;
stats.n_iters         = numel(iters);
stats.logicalIters    = j;
stats.faultIterations = errorIterations;
stats.r0              = r0;
stats.sigmaHist       = sigHist;
stats.sigmaFinal      = sigma;
stats.nRefactor       = nRefactor;
stats.blockWidth      = size(Z, 2);
stats.innerSolve      = opts.innerSolve;
stats.setupTime       = setupTime;
stats.iterWall        = iterWall;
stats.Xtrue           = X;
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
    'solvTime', max(e2eTime - setupTime - checkpointTime - recoveryTime, 0), ...
    'e2eTime', e2eTime);

fprintf('\n[CR-SUMMARY] mode %s | interval %g | checkpoints %d | faults %d | lost iters %d\n', ...
        cr.interval_mode, T, nCkpt, numel(errorIterations), nnz(wasted));
fprintf('e2eTime(s)\tcheckpointTime(s)\trecoveryTime(s)\treworkTime(s)\n');
fprintf('%f\t%f\t%f\t%f\n', e2eTime, checkpointTime, recoveryTime, reworkTime);

    function P = packCkpt(jj)
        st = state;
        if isstruct(st)
            if isfield(st, 'Xtrue'),     st.Xtrue = [];          end
            if isfield(st, 'lastInner'), st.lastInner = struct(); end
        end
        P = struct('j', jj, 'Z', Z, 'sigma', sigma, 'r0', r0, 'state', {st}, ...
                   'nRefactor', nRefactor, 'resLog', resLog, 'lastFault', lastFault);
    end
end

% Setup: block padding and factorization of A.
function [Aw, Bw, Z, ops, sigma] = tmSetup(A, B, X0, s, opts)
    n = size(A, 1);
    if ~isempty(opts.seed), rng(opts.seed); end
    m = min(max(size(X0, 2), opts.blockFactor * s), n);
    Z = X0;
    if size(Z, 2) < m
        padSeed = 42; if ~isempty(opts.seed), padSeed = opts.seed; end
        Z = [Z, randn(RandStream('mt19937ar', 'Seed', padSeed), n, m - size(Z, 2))];
    end
    Aw = A;  Bw = B;
    F0 = spdFactor(Aw);
    if ~F0.ok, error('tracemin_shifted_checkpoint_restart:notSPD', 'A is not SPD.'); end
    sigma = 0;
    ops   = buildOps(Aw, Bw, F0, sigma);
end

function v = getf(S, f)
    if isstruct(S) && isfield(S, f), v = S.(f); else, v = 0; end
end

function ops = buildOps(A, B, F, sigma)
    ops = struct('applyA',    @(U) A * U, ...
                 'applyB',    @(U) B * U, ...
                 'applyPrec', F.solve, ...
                 'sigma',     sigma);
    ops.precExact = true;
    ops.fac       = F;
    ops.diagA     = full(diag(A));
    ops.diagB     = full(diag(B));
    ops.trueRes   = [];
end

function F = spdFactor(M)
    if issparse(M), M = (M + M') / 2; end
    if isdiag(M)
        dM = full(diag(M));
        F.ok = all(dM > 0);
        F.solve = @(U) U ./ dM;
        F.kind = 'diag';
        return;
    end
    if issparse(M)
        [R, p, q] = chol(M, 'vector');
        F.kind = 'sparse';
    else
        [R, p] = chol(M);
        q = 1:size(M, 1);
        F.kind = 'dense';
    end
    F.ok = (p == 0);
    if F.ok
        F.solve = @(U) cholSolve(R, q, U);
        F.R = R;  F.q = q;
    else
        F.solve = [];
    end
end

function Xs = cholSolve(R, q, U)
    Xs = zeros(size(U));
    Xs(q, :) = R \ (R' \ U(q, :));
end
