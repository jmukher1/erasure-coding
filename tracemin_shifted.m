function [X, Thi, iters, res, relres, errorIterations, eCodingGenTime, decomp2Time, restructTime, stats, Xtrue] = ...
    tracemin_shifted(A, B, X0, normA0, Za, Ra, Qb, Sb, k, s, ...
    possibleErrorIterations, rowDeletedStart, TOL, maxiter, opts, E, SIGMA, sopts)
% tracemin_shifted - TraceMin iteration with erasure coding and an optional
% shift sigma (preconditioner = exact factor of A - sigma*B).
%
% Computes the s smallest eigenpairs of (A, B); only SIGMA = 'smallestabs'
% is supported.
%
% At a fault iteration rows/cols F = rowDeletedStart..+k-1 of A and B are
% erased. The iteration then continues on the coded pencil
% A' = P'AP, B' = P'BP, P = [I(:,C), E], applied through the coding blocks
% Za = A*E, Ra = E'AE, Qb = B*E, Sb = E'BE (never formed explicitly). The
% shifted solve on the coded pencil uses a Cholesky factor of the surviving
% block K_CC and a k x k Schur complement. Rows F of the iterate are refilled
% per opts.faultReset. Xtrue = P*X is the eigenvector estimate of the
% original problem.
%
% Timing returned to tracemin_erasure_shifted:
%   restructTime  fault handling (coded operators + refactorization),
%                 minus simTime and diagTime
%   decomp2Time   the refactorization part of restructTime
%   stats.timing.simTime   zeroing the erased rows/cols (fault simulation)
%   stats.timing.diagTime  eigenvector recovery Xtrue = P*X (diagnostic)
% simTime and diagTime are excluded from e2eTime by the caller.

innerSolve = 'minres';
if nargin >= 15 && isstruct(opts) && isfield(opts, 'innerSolve')
    innerSolve = opts.innerSolve;
    opts = rmfield(opts, 'innerSolve');
end
if nargin < 15 || isempty(opts)
    opts = tracemin_opts();
else
    opts = tracemin_opts(opts);
end
opts.innerSolve = innerSolve;
if nargin < 16, E = []; end
if nargin < 17 || isempty(SIGMA), SIGMA = 'smallestabs'; end
if nargin < 18 || isempty(sopts), sopts = struct(); end
if ~isfield(sopts, 'postFault'), sopts.postFault = 'coded'; end
if ~isfield(sopts, 'shiftHold'), sopts.shiftHold = 1;       end

tFun0   = tic;
if ~strcmpi(SIGMA, 'smallestabs')
    error('tracemin_shifted:sigma', 'only SIGMA = ''smallestabs'' is supported.');
end
[n, ~]  = size(A);

% block width m = blockFactor*s; extra columns are seeded random vectors
if ~isempty(opts.seed), rng(opts.seed); end
m = min(max(size(X0, 2), opts.blockFactor * s), n);
Z = X0;
if size(Z, 2) < m
    padSeed = 42; if ~isempty(opts.seed), padSeed = opts.seed; end
    Z = [Z, randn(RandStream('mt19937ar', 'Seed', padSeed), n, m - size(Z, 2))];
end

F0 = spdFactor(A);
if ~F0.ok
    error('tracemin_shifted:notSPD', 'A is not SPD.');
end

r0 = 0; state = []; Fset = []; op = struct(); Bit = B; cod = [];
relres = []; Xtrue = []; res = []; iters = [];
errorIterations = []; eCodingGenTime = 0; restructTime = 0; decomp2Time = 0;
sigHist = []; codedHist = []; nRefactor = 0; stagnated = false; lastFault = 0;
diagTime = 0;
simTime  = 0;
nFullRefactor = 0;

ITER_MAX = min(maxiter, n^3);
fprintf('\nInside tracemin_shifted (%s) with tol %d', SIGMA, TOL);

sigma = 0;
ops   = buildOps(A, B, F0, sigma);
F0    = [];

fprintf('\nTraceMin Method (shifted)\n')
fprintf('\nIter \t, ResNorm \t, RelResNorm \t, AbsResNorm \t, Tolerance \t, sigma\n');

setupTime = toc(tFun0);
tLoop0    = tic;
iterWall  = [];

for i1 = 1:ITER_MAX
    % fault injection (timed as restructTime)
    if any(possibleErrorIterations == i1)
        errorIterations(end+1) = i1;
        lastFault = i1;
        fprintf("\nNumber of rows/columns erased = %d starting from row: %d in iteration %d\n", k, rowDeletedStart, i1);
        t2 = tic;
        diagAtFault = diagTime;
        simAtFault  = simTime;
        Ffault   = (rowDeletedStart : (rowDeletedStart + k - 1))';
        Chealthy = setdiff((1:n)', Ffault);

        if isempty(E)
            warning('tracemin_shifted:noCoding', 'Fault at iteration %d but no coding matrix E: erasure ignored.', i1);
        elseif ~isempty(Fset)
            % later fault on the already active coded pencil: the operator is
            % unchanged, only the lost rows of the iterate are refilled
            if ~isequal(Ffault, Fset)
                warning('tracemin_shifted:faultWindow', ...
                        ['Fault at iteration %d hits rows %d..%d, but the coded pencil is ' ...
                         'built on rows %d..%d; only Z(F,:) is reset.'], ...
                        i1, Ffault(1), Ffault(end), Fset(1), Fset(end));
            end
            Z = resetFaultRows(Z, Ffault, Chealthy, opts.faultReset, ops);
            fprintf('    [fault] coded pencil already active; Z(F,:) reset (%s) on rows %d..%d\n', ...
                    opts.faultReset, Ffault(1), Ffault(end));
        else
            if ~strcmp(sopts.postFault, 'coded')
                warning('tracemin_shifted:postFault', ...
                        'postFault = ''%s'' is not implemented; using the coded pencil.', sopts.postFault);
            end
            % first fault: switch to the coded pencil
            Ek = E(:, 1:k);
            Fset = Ffault;
            cod  = struct('E', Ek, 'F', Ffault, 'C', Chealthy, 'n', n, ...
                          'Za', Za(:, 1:k), 'Ra', Ra(1:k, 1:k), 'Qb', Qb(:, 1:k), 'Sb', Sb(1:k, 1:k));
            op   = cod;

            % simulate the erasure (timed as simTime)
            tsim = tic;
            A(Ffault, :) = 0;  A(:, Ffault) = 0;
            B(Ffault, :) = 0;  B(:, Ffault) = 0;
            simTime = simTime + toc(tsim);

            Bit  = codedOperator(makeCodedApply(B, cod.F, cod.Qb, cod.Sb), n);

            % factor the coded shifted operator (timed as decomp2Time)
            tfac = tic;
            trueRes = struct('E', Ek, 'F', Ffault);
            [opsN, okN] = codedOps(A, B, cod, sigma, trueRes);
            if ~okN && sigma ~= 0
                fprintf('    [fault] coded A - sigma*B not SPD at sigma = %.6e, dropping the shift\n', sigma);
                sigma = 0;  [opsN, okN] = codedOps(A, B, cod, sigma, trueRes);
            end
            if ~okN
                error('tracemin_shifted:codedNotSPD', ...
                      'Coded pencil is not SPD (cond(E_F) too large?).');
            end
            ops = opsN;  opsN = [];
            decomp2Time = decomp2Time + toc(tfac);
            nFullRefactor = nFullRefactor + 1;

            Z = resetFaultRows(Z, Ffault, Chealthy, opts.faultReset, ops);

            fprintf('    [fault] coded pencil active on rows %d..%d, Z(F,:) reset (%s)\n', ...
                    Ffault(1), Ffault(end), opts.faultReset);
        end

        restructTime = restructTime + toc(t2) - (diagTime - diagAtFault) - (simTime - simAtFault);
    end

    % one TraceMin iteration
    [Z, bStop, Thi, X, r0, resNorm, relresNorm, ~, state] = ...
        tracemin_single_iter_shifted(A, [], Bit, Z, i1, s, r0, TOL, Fset, op, ...
                                           normA0, opts, ops, state);

    res(end+1)     = resNorm;
    relres(end+1)  = relresNorm;
    iters(end+1)   = i1;
    iterWall(end+1) = toc(tLoop0);
    sigHist(end+1) = sigma;
    codedHist(end+1) = state.codedResNorm;
    if isfield(state, 'Xtrue') && ~isempty(state.Xtrue), Xtrue = state.Xtrue; end

    if bStop == 1
        fprintf('\nTracemin (shifted, %s) converged at iteration %d\n', SIGMA, i1);
        break;
    end

    % accept the proposed shift if it gains enough; refactor A - sigma*B
    %      (backs off if the shifted matrix is not SPD)
    if ~strcmp(opts.shiftMode, 'none') && nRefactor < opts.maxRefactor && ...
            (lastFault == 0 || i1 - lastFault >= sopts.shiftHold)
        sigNew = state.sigmaProposed;
        th1    = state.theta1;
        if sigNew > sigma && (th1 - sigNew) < opts.refactorGain * (th1 - sigma)
            for backoff = 1:4
                if isempty(Fset)
                    Fn = spdFactor(A - sigNew * B);
                    okN = Fn.ok;
                    if okN, opsN = buildOps(A, B, Fn, sigNew, ops.trueRes); end
                    Fn = [];
                else
                    [opsN, okN] = codedOps(A, B, cod, sigNew, ops.trueRes);
                end
                if okN
                    sigma = sigNew;
                    ops   = opsN;  opsN = [];
                    nRefactor = nRefactor + 1;
                    if opts.verbose
                        fprintf('    [shift] sigma -> %.10e (theta_1 = %.10e, refactor #%d)\n', ...
                                sigma, th1, nRefactor);
                    end
                    break;
                end
                sigNew = (sigma + sigNew) / 2;
            end
        end
    end

    % stop if the residual has stopped improving since the last fault
    w = opts.stagnWindow;
    w0 = max(1, lastFault);
    if isfinite(w) && i1 > 2 * w && (i1 - lastFault) > 2 * w && ...
            min(res(end-w+1:end)) > opts.stagnFactor * min(res(w0:end-w))
        stagnated = true;
        fprintf('\nStagnated: no %.0f%% improvement in %d iterations, stopping at %d\n', ...
                100 * (1 - opts.stagnFactor), w, i1);
        break;
    end
end

% eigenvector estimate of the original problem, Xtrue = P*X (timed as diagTime)
if exist('trueRes', 'var') == 1 && isstruct(trueRes) && ~isempty(trueRes.F)
    tvr = tic;
    Xtrue = recoverTrueX(X, trueRes.E, trueRes.F);
    diagTime = diagTime + toc(tvr);
end

sk  = min(s, size(Thi, 1));
Thi = Thi(1:sk, 1:sk);
X   = X(:, 1:sk);
if ~isempty(Xtrue), Xtrue = Xtrue(:, 1:min(s, size(Xtrue, 2))); end

stats = struct();
stats.iters           = iters;
stats.res             = res;
stats.relres          = relres;
stats.eigenvalues     = diag(Thi);
stats.converged       = (exist('bStop', 'var') == 1) && bStop == 1;
stats.stagnated       = stagnated;
stats.n_iters         = numel(iters);
stats.faultIterations = errorIterations;
stats.r0              = r0;
stats.normM = [];  stats.normMinv = [];
if exist('Ek', 'var') == 1 && exist('Ffault', 'var') == 1 && ~isempty(Ek)
    [stats.normM, stats.normMinv] = pencilCongruenceNorms(Ek, Ffault);
end
stats.sigma           = SIGMA;
stats.sigmaHist       = sigHist;
stats.sigmaFinal      = sigma;
stats.codedres        = codedHist;
stats.nRefactor       = nRefactor;
stats.blockWidth      = m;
stats.innerSolve      = opts.innerSolve;
stats.setupTime       = setupTime;
stats.iterWall        = iterWall;
stats.postFault       = sopts.postFault;
if ~isempty(state) && isfield(state, 'innerIters')
    stats.innerIters      = state.innerIters;
    stats.innerItersTotal = sum(state.innerIters);
    stats.nShifted        = state.nShifted;
    stats.nLocked         = state.nLocked;
    stats.nDeflated       = state.nDeflated;
end
stats.timing = struct('eCodingGenTime', eCodingGenTime, ...
                      'restructTime',   restructTime, ...
                      'refactorTime',   decomp2Time, ...
                      'diagTime',       diagTime, ...
                      'simTime',        simTime, ...
                      'nFullRefactor',  nFullRefactor);
end

% Refill the erased rows F of the iterate block Z.
function Z = resetFaultRows(Z, F, C, mode, ops)
    switch mode
        case 'random'
            rmsC = sqrt(mean(Z(C, :).^2, 1));
            Z(F, :) = randn(numel(F), size(Z, 2)) .* rmsC;
        case 'zero'
            Z(F, :) = 0;
        case 'jacobi'
            Z(F, :) = 0;
            AZ  = ops.applyA(Z);
            dAF = ops.diagA(F);
            Zf  = -AZ(F, :) ./ dAF;
            if all(isfinite(Zf(:))) && norm(Zf, 'fro') <= norm(Z(C, :), 'fro')
                Z(F, :) = Zf;
            end
    end
end

% Handle U -> (P'MP)*U for M with rows/cols F zeroed, using ME = M_orig*E
% and Smm = E'*M_orig*E.
function h = makeCodedApply(M, F, ME, Smm)
    MEt = ME.';
    D   = full(Smm - ME(F, :) - ME(F, :).');
    h   = @(U) codedApplyPre(M, F, ME, MEt, D, U);
end

function Y = codedApplyPre(M, F, ME, MEt, D, U)
    UF = U(F, :);
    Y  = M * U + ME * UF;
    Y(F, :) = Y(F, :) + MEt * U + D * UF;
end

% Operator bundle for the coded pencil at shift sigma. K' = A' - sigma*B' is
% solved by block elimination: Cholesky of K_CC = (A - sigma*B)(C,C) and the
% k x k Schur complement SC = K_FF - K_CF' * K_CC^{-1} * K_CF.
function [ops, ok] = codedOps(A, B, cod, sigma, trueRes)
    ops = [];  C = cod.C;  F = cod.F;  n = cod.n;  E = cod.E;  k = numel(F);
    FC = spdFactor(A(C, C) - sigma * B(C, C));
    ok = FC.ok;
    if ~ok, return; end
    if sigma ~= 0
        KCFm = cod.Za(C, :) - sigma * cod.Qb(C, :);
    else
        KCFm = cod.Za(C, :);
    end
    keepW = numel(C) * k <= 8e7;
    if keepW, W = zeros(numel(C), k); else, W = []; end
    SC   = cod.Ra - sigma * cod.Sb;
    for j0 = 1:64:k
        J  = j0:min(j0 + 63, k);
        WJ = FC.solve(full(KCFm(:, J)));
        SC(:, J) = SC(:, J) - KCFm' * WJ;
        if keepW, W(:, J) = WJ; end
    end
    SC   = (SC + SC') / 2;
    [L, p] = chol(SC, 'lower');
    ok = (p == 0);
    if ~ok, return; end
    dA = full(diag(A));  dA(F) = diag(cod.Ra);
    dB = full(diag(B));  dB(F) = diag(cod.Sb);
    ops = struct('applyA',    makeCodedApply(A, F, cod.Za, cod.Ra), ...
                 'applyB',    makeCodedApply(B, F, cod.Qb, cod.Sb), ...
                 'applyPrec', @(X) codedSolve(FC, KCFm, W, L, C, F, n, X), ...
                 'sigma',     sigma);
    ops.precExact = true;
    ops.fac       = FC;
    ops.diagA     = dA;
    ops.diagB     = dB;
    ops.trueRes   = trueRes;
end

% Z = (A' - sigma*B') \ X, with W = K_CC^{-1} K_CF when it was kept.
function Z = codedSolve(FC, KCFm, W, L, C, F, n, X)
    YC = FC.solve(X(C, :));
    ZF = L' \ (L \ (X(F, :) - KCFm' * YC));
    Z  = zeros(n, size(X, 2));
    if ~isempty(W)
        Z(C, :) = YC - W * ZF;
    else
        Z(C, :) = FC.solve(X(C, :) - KCFm * ZF);
    end
    Z(F, :) = ZF;
end

% Operator bundle for the original pencil; F is the factor of A - sigma*B.
function ops = buildOps(A, B, F, sigma, trueRes)
    if nargin < 5, trueRes = []; end
    ops = struct('applyA',    @(U) A * U, ...
                 'applyB',    @(U) B * U, ...
                 'applyPrec', F.solve, ...
                 'sigma',     sigma);
    ops.precExact = true;
    ops.fac       = F;
    ops.diagA   = full(diag(A));
    ops.diagB   = full(diag(B));
    ops.trueRes = trueRes;
end

% Cholesky (fill-reducing ordering for sparse M); F.ok = false if not SPD.
function F = spdFactor(M)
    if issparse(M)
        M = (M + M') / 2;
    end
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
