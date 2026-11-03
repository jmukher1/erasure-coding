function [Z, bStop, Thi, X, r0, resNorm, relresNorm, absresNorm, state] = ...
        tracemin_single_iter_shifted(A, dA, B, Z, iter, s, r0, TOL, Fset, op, normA0, opts, ops, state)
% tracemin_single_iter_shifted - one TraceMin iteration on the pencil (A, B).
%
%   1. B-orthonormalize the block Z -> V
%   2. Rayleigh-Ritz: H = V'*A*V, Ritz pairs (Thi, X = V*Y)
%   3. residual of the first s Ritz pairs, relative to r0 (r0 is set at iter 1)
%   4. TraceMin update Z = X - Delta with the exact (Schur) correction
%        Y = K^{-1} B X,  Z = Y (X'B K^{-1} B X)^{-1},   K = A - sigma*B
%
% A, B and the preconditioner are applied through the operator bundle ops
% (applyA, applyPrec, sigma), built by tracemin_shifted for either the
% original or the coded pencil. B may be a codedOperator after a fault.
% state carries per-iteration diagnostics and the proposed next shift.

    [opts, innerSolve] = splitInnerSolve(opts);
    if nargin < 14 || isempty(state)
        state = struct('nDeflated', 0, 'innerIters', [], 'nShifted', [], ...
                       'nLocked', [], 'perColRes', [], ...
                       'sigmaProposed', 0, 'thetaS', Inf, 'lastInner', struct(), ...
                       'theta1', Inf, 'r1rel', Inf, 'thetaM', Inf, ...
                       'Xtrue', [], 'codedResNorm', []);
    end
    if ~isfield(ops, 'sigma') || isempty(ops.sigma)
        ops.sigma = 0;
    end

    % B-orthonormalization and Rayleigh-Ritz
    [V, nDef, BV] = BOrthonormalizeStable(B, Z, opts.dropTol);
    state.nDeflated = state.nDeflated + nDef;

    m    = size(V, 2);
    sEff = min(s, m);

    W = ops.applyA(V);
    H = V' * W;
    H = (H + H') / 2;
    [Y, Thi] = eig(H);
    [ssort, idx] = sort(real(diag(Thi)));
    Thi = diag(ssort);
    Y   = Y(:, idx);

    X     = V * Y;
    BX    = BV * Y;
    AX    = W * Y;
    theta = ssort;

    % per-column relative residuals (used for locking)
    Rcols  = AX - BX .* theta';
    perCol = sqrt(sum(Rcols .* Rcols, 1))';
    denom  = max(abs(theta) .* sqrt(sum(X .* X, 1))', realmin);
    perCol = perCol ./ denom;

    % proposed shift for the next iteration (accepted or not by the caller)
    r1 = norm(AX(:, 1) - theta(1) * BX(:, 1));
    r1rel = r1 / max(abs(theta(1)), realmin);
    if r1rel < opts.shiftTol
        sigCand = theta(1) - max(opts.shiftSafety * r1, opts.shiftFloor * abs(theta(1)));
        state.sigmaProposed = max(sigCand, ops.sigma);
    else
        state.sigmaProposed = ops.sigma;
    end
    state.thetaS = theta(min(sEff, m));
    state.theta1 = theta(1);
    state.thetaM = theta(min(m, numel(theta)));
    state.r1rel  = r1rel;

    % residual of the first s Ritz pairs
    Res     = AX(:, 1:sEff) - BX(:, 1:sEff) * Thi(1:sEff, 1:sEff);
    resNorm = norm(Res, 'fro');
    hasTrueRes = isfield(ops, 'trueRes') && ~isempty(ops.trueRes) && ~isempty(ops.trueRes.F);
    if hasTrueRes
        state.Xtrue = [];
    else
        state.Xtrue = X;
    end
    state.codedResNorm = resNorm;

    if iter == 1
        r0 = resNorm;
    end
    absresNorm = resNorm / normA0;
    relresNorm = resNorm / r0;

    bStop = 0;
    if relresNorm < TOL
        bStop = 1;
        state.perColRes(:, end+1) = perCol;
        if opts.verbose
            fprintf('%d \t, %e \t, %e \t, %e \t, %e (converged)\n', ...
                    iter, resNorm, relresNorm, absresNorm, TOL);
        end
        return;
    end

    locked = perCol < opts.lockTol;
    switch opts.shiftMode
        case 'none'
            shiftOn = false(m, 1);
        case 'common'
            shiftOn = repmat(ops.sigma ~= 0, m, 1);
        otherwise
            error('tracemin:badShiftMode', 'shiftMode must be ''common'' or ''none''.');
    end
    if ~strcmp(innerSolve, 'schur') || ~isfield(ops, 'precExact') || ~isequal(ops.precExact, true)
        error('tracemin:innerSolve', 'only the exact Schur inner solve (innerSolve = ''schur'') is supported.');
    end

    % TraceMin update with the exact correction
    Y = ops.applyPrec(BX);
    S = BX' * Y;
    S = (S + S') / 2;
    [Rs, pS] = chol(S);
    if pS == 0
        Z = (Y / Rs) / Rs';
    else
        Z = Y / S;
    end

    state.innerIters(end+1)   = m;
    state.nShifted(end+1)     = sum(shiftOn);
    state.nLocked(end+1)      = sum(locked);
    state.perColRes(:, end+1) = perCol;
    state.lastInner           = struct('totalIters', m, 'method', 'schur');

    if opts.verbose
        fprintf('%d \t, %e \t, %e \t, %e \t, %e \t, sigma = %.6e\n', ...
                iter, resNorm, relresNorm, absresNorm, TOL, ops.sigma);
    end
end

function [opts, innerSolve] = splitInnerSolve(opts)
    innerSolve = 'minres';
    if nargin < 1 || isempty(opts)
        opts = tracemin_opts();
    else
        if isfield(opts, 'innerSolve')
            innerSolve = opts.innerSolve;
            opts = rmfield(opts, 'innerSolve');
        end
        if ~isfield(opts, 'blockFactor')
            opts = tracemin_opts(opts);
        end
    end
    opts.innerSolve = innerSolve;
end
