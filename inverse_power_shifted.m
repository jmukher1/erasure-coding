function [X, D, Xr, F, iters, relres, errorIterations, eCodingGenTime, augTime, solvTime, vectorRecoveryTime, restructTime, e2eTime, stats] = ...
    inverse_power_shifted(A, X0, k, possibleErrorIterations, rowDeletedStart, varargin)
% inverse_power_shifted - erasure-coded block inverse power method
% (smallest eigenvalues) with shift-and-invert, Z = (A' - sigma*B')^{-1} B'X.
%
    t = datetime('now');

    p = inputParser;
    addParameter(p, 'tol',          1e-8);
    addParameter(p, 'maxiter',      1000);
    addParameter(p, 'sigma',        0);
    addParameter(p, 'blockFactor',  2);
    addParameter(p, 'shiftMode',    'auto');
    addParameter(p, 'shiftTol',     1e-2);
    addParameter(p, 'shiftGap',     1e-8);
    addParameter(p, 'rrMode',       'gram');
    addParameter(p, 'reportSpike',  true);
    addParameter(p, 'faultReset',   'random');
    addParameter(p, 'xRecovery',    'local');
    addParameter(p, 'padSeed',      42);
    addParameter(p, 'shiftSafety',  1.1);
    addParameter(p, 'shiftFloor',   1e-8);
    addParameter(p, 'refactorGain', 0.9);
    addParameter(p, 'maxRefactor',  30);
    addParameter(p, 'stagnWindow',  20);
    addParameter(p, 'stagnFactor',  0.9);
    addParameter(p, 'coding',       'interleaved');
    addParameter(p, 'codingScheme', 'Random');
    parse(p, varargin{:});

    tol          = p.Results.tol;
    maxiter      = p.Results.maxiter;
    sigma        = p.Results.sigma;
    blockFactor  = p.Results.blockFactor;
    shiftMode    = p.Results.shiftMode;
    stagnWindow  = p.Results.stagnWindow;
    stagnFactor  = p.Results.stagnFactor;
    shiftSafety  = p.Results.shiftSafety;
    refactorGain = p.Results.refactorGain;
    maxRefactor  = p.Results.maxRefactor;
    shiftTol     = p.Results.shiftTol;
    shiftGap     = p.Results.shiftGap;
    rrMode       = p.Results.rrMode;
    if ~ismember(rrMode, {'gram', 'qr'})
        error('inverse_power_shifted:badRrMode', 'rrMode must be ''gram'' or ''qr''.');
    end
    if ~ismember(shiftMode, {'auto', 'fixed'})
        error('inverse_power_shifted:badShiftMode', 'shiftMode must be ''auto'' or ''fixed''.');
    end
    nGramFallback = 0;

    [n, s] = size(X0);
    if blockFactor < 1
        error('inverse_power_shifted:badBlockFactor', 'blockFactor must be >= 1.');
    end

    % block width m = blockFactor*s; extra columns are seeded random vectors
    m = min(round(blockFactor * s), n);
    if size(X0, 2) < m
        sc = norm(X0, 'fro') / max(sqrt(numel(X0)), 1);
        if ~isfinite(sc) || sc <= 0
            sc = 1;
        end
        padStream = RandStream('mt19937ar', 'Seed', p.Results.padSeed);
        X0 = [X0, sc * randn(padStream, n, m - size(X0, 2))];
    elseif size(X0, 2) > m
        X0 = X0(:, 1:m);
    end
    if m > s
        fprintf('inverse power: oversampling to m = %d columns (= %g * s, s = %d)\n', ...
                m, blockFactor, s);
    end

    eCodingGenTime = 0; augTime = 0; restructTime = 0;
    solvTime = 0; e2eTime = 0; vectorRecoveryTime = 0;
    simTime  = 0;
    diagTime = 0;
    nRecovered = 0;

    %  coding matrix E (n x k)
    t01 = datetime('now');
    E0 = makeErasureCodingMatrix(A, n, k, p.Results.coding, p.Results.codingScheme);
    t1  = datetime('now');
    eCodingGenTime = seconds(t1 - t01);

    %  coding blocks R = A*E, S = E'*A*E, T = E'*E
    R0 = A * E0;
    S0 = E0' * R0;
    T0 = E0' * E0;
    t3 = datetime('now');
    augTime = seconds(t3 - t1);

    normA  = norm(A, 'fro');
    normA1 = norm(A, 1);
    r0 = 0;
    stagnatedFlag = false;

    [X, ~] = qr(X0, 0);

    F = [];   C = (1:n)';
    E = [];   R = [];   S = [];   T = [];
    nRefactor = 0;

    [solveShift, okShift, sigma] = ipmShiftedSolver(A, sigma, n, E0, R0, S0, T0, [], (1:n)', []);
    if ~okShift
        error('inverse_power_shifted:indefiniteShift', ...
              'A - %g*I is not positive definite; supply a smaller sigma.', sigma);
    end

    relres = [];
    absres = [];
    iters  = [];
    errorIterations = [];
    BXprev = [];
    spikeIter = []; spikeAbs = []; spikeRel = [];

    fprintf('iter\tResidual Norm\tRelative Residual Norm\n');

    for iter = 1:maxiter

        %  fault injection
        Fnew = faultfunc_matlab2(iter, possibleErrorIterations, ...
                                 rowDeletedStart, rowDeletedStart + k - 1);

        if ~isempty(Fnew)
            fprintf('\n new faults detected at iter %d \n', iter);
            t11 = datetime('now');
            simBefore  = simTime;
            diagBefore = diagTime;

            Fprev = F;
            F  = sort(union(F, Fnew));
            nF = length(F);
            assert(nF <= k, 'inverse_power_shifted:tooManyFaults', ...
                'Accumulated faults (%d) exceed tolerated k (%d).', nF, k);
            errorIterations(end+1) = iter;

            % simulate the fault: reset the lost iterate rows, zero rows/cols F of A
            tsim = datetime('now');
            survivingRows = setdiff(1:n, Fnew);
            resetScale = norm(X(survivingRows, :), 'fro') / sqrt(numel(survivingRows) * size(X, 2));
            switch p.Results.faultReset
                case 'zero'
                    X(Fnew, :) = 0;
                otherwise
                    X(Fnew, :) = resetScale * randn(length(Fnew), size(X, 2));
            end
            A(Fnew, :) = 0;
            A(:, Fnew) = 0;
            simTime = simTime + seconds(datetime('now') - tsim);

            E = E0(:, 1:nF);
            R = R0(:, 1:nF);
            S = S0(1:nF, 1:nF);
            T = T0(1:nF, 1:nF);
            C = setdiff(1:n, F)';
            assert(isequal(F(:)', F(1):F(1)+nF-1), ...
                'inverse_power_shifted:nonContiguousF', 'the fault block must be contiguous.');

            % residual of the damaged iterate on (A', B')
            if p.Results.reportSpike && exist('D', 'var') && iter > 1
                td = tic;
                Xd = X(:, 1:s) ./ max(vecnorm(X(:, 1:s)), realmin);
                AXd = Aop(A, R, S, Xd, F);
                BXd = Bop(E, T, Xd, F);
                [spikeAbs, spikeRel] = compute_residual(AXd, BXd, diag(D(1:s)), r0, normA);
                spikeIter = iter;
                diagTime = diagTime + toc(td);
            end

            if strcmp(p.Results.xRecovery, 'local') && exist('D', 'var') && numel(D) >= size(X, 2)
                Cr  = C;
                ACC = A(Cr, Cr);  RC = R0(Cr, 1:nF);  EC = E0(Cr, 1:nF);
                Sk  = S0(1:nF, 1:nF);  Tk = T0(1:nF, 1:nF);
                YC  = X(Cr, :);   th = D(1:size(X, 2));
                AY  = ACC * YC;   RtY = RC' * YC;   EtY = EC' * YC;
                first = isempty(Fprev);
                if first
                    G1  = RC - ACC * EC;
                    rz  = find(any(G1, 2));
                    G1d = full(G1(rz, :));
                    G2a = full(Sk - RC' * EC);  G2b = full(Tk - EC' * EC);
                else
                    rz  = zeros(0, 1);
                    G1d = zeros(0, nF);
                    G2a = full(Sk);  G2b = full(Tk);
                end
                nRec = 0;  YF = zeros(nF, size(X, 2));  okc = false(1, size(X, 2));
                for j = 1:size(X, 2)
                    Gs  = [G1d; G2a - th(j) * G2b];
                    b1  = th(j) * YC(:, j) - AY(:, j);
                    rhs = [b1(rz); th(j) * EtY(:, j) - RtY(:, j)];
                    [Qs, Rs] = qr(Gs, 0);
                    dR = abs(diag(Rs));
                    if ~isempty(dR) && min(dR) > 1e-12 * max(dR)
                        YF(:, j) = Rs \ (Qs' * rhs);  okc(j) = true;  nRec = nRec + 1;
                    end
                end
                if first
                    X(Cr, okc) = YC(:, okc) - EC * YF(:, okc);
                end
                X(F, okc) = YF(:, okc);
                nRecovered = nRecovered + nRec;
                fprintf('  [xRecovery] rebuilt erased iterate rows of %d/%d columns from coded data\n', nRec, size(X, 2));
            end

            [solveShift, okShift, sigma] = ipmShiftedSolver(A, sigma, n, E0, R0, S0, T0, F, C, 0);
            if ~okShift
                error('inverse_power_shifted:codedNotSPD', ...
                      'Coded pencil A'' - sigma*B'' is not SPD even at sigma = 0.');
            end
            t12 = datetime('now');
            restructTime = restructTime + seconds(t12 - t11) - (simTime - simBefore) - (diagTime - diagBefore);
            BXprev = [];
        end

        if isempty(BXprev)
            BX = Bop(E, T, X, F);
        else
            BX = BXprev;
        end

        t_s = datetime('now');
        Z = solveShift(BX);
        solvTime = solvTime + seconds(datetime('now') - t_s);

        % Rayleigh-Ritz on (A', B')
        gramOK = false;
        if strcmp(rrMode, 'gram')
            Q = Z ./ max(vecnorm(Z), realmin);
            AQ_mat = Aop(A, R, S, Q, F);
            BQ_mat = Bop(E, T, Q, F);
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
            AQ_mat = Aop(A, R, S, Q, F);
            BQ_mat = Bop(E, T, Q, F);
            AQ = Q' * AQ_mat;  AQ = (AQ + AQ') / 2;
            BQ = Q' * BQ_mat;  BQ = (BQ + BQ') / 2;
            [U, Dmat] = eig(AQ, BQ);
        end
        [D, idx] = sort(real(diag(Dmat)), 'ascend');
        U = U(:, idx);
        X = Q * U;
        AX = AQ_mat * U;
        BX = BQ_mat * U;
        BXprev = BX;

        % residual of the current pencil, relative to r0
        AXres = AX(:, 1:s);
        BXres = BX(:, 1:s);
        if iter == 1
            r0 = norm(AXres - BXres .* D(1:s)', 'fro');
        end
        [resNorm, relresNorm] = compute_residual(AXres, BXres, diag(D(1:s)), r0, normA);
        if ~isempty(spikeIter) && spikeIter == iter
            fprintf('  [reportSpike] fault iteration reports the damaged-iterate residual %.3e (recovered: %.3e)\n', ...
                    spikeRel, relresNorm);
            resNorm = spikeAbs;  relresNorm = spikeRel;
            spikeIter = [];
        end
        relres(end+1) = relresNorm;
        absres(end+1) = resNorm;
        iters(end+1)  = iter;

        sinceFault = iter;
        if ~isempty(errorIterations)
            sinceFault = iter - max(errorIterations);
        end
        if isfinite(stagnWindow) && numel(relres) > stagnWindow && ...
           sinceFault > stagnWindow && relresNorm >= tol
            recent = relres(end - stagnWindow);
            if isfinite(recent) && recent > 0 && relresNorm > stagnFactor * recent
                fprintf(['\nInverse power stopped at iteration %d: residual stalled at ' ...
                         '%.3e (abs %.3e), no %.0f%% improvement over %d iterations (sigma = %.10g).\n'], ...
                        iter, relresNorm, resNorm, 100*(1-stagnFactor), stagnWindow, sigma);
                stagnatedFlag = true;
                break;
            end
        end

        fprintf('%d\t%e\t%e\t%e\tsigma=%.6g\n', iter, resNorm, relresNorm, tol, sigma);

        if strcmp(shiftMode, 'auto') && nRefactor < maxRefactor && relresNorm >= tol
            th1 = D(1);
            r1  = norm(AX(:, 1) - th1 * BX(:, 1));
            if th1 > 0 && r1 / th1 < shiftTol
                cand = min([th1 - shiftSafety * r1, (1 - shiftGap) * th1, ...
                            th1 - 100 * eps * normA1]);
                if cand > sigma && (th1 - cand) < refactorGain * (th1 - sigma)
                    [solveNew, okNew, sigNew] = ...
                        ipmShiftedSolver(A, cand, n, E0, R0, S0, T0, F, C, sigma);
                    if okNew && sigNew > sigma
                        solveShift = solveNew;
                        sigma      = sigNew;
                        nRefactor  = nRefactor + 1;
                        fprintf('  [shift] sigma -> %.10g (theta_1 = %.10g, refactor #%d)\n', ...
                                sigma, th1, nRefactor);
                    end
                end
            end
        end

        if relresNorm < tol
            fprintf('Converged in %d iterations (rel. residual = %e)\n', iter, relresNorm);
            break;
        end
    end

    t_rec = datetime('now');
    Xr = X;
    if ~isempty(F)
        nF   = length(F);
        E_r  = E0(:, 1:nF);
        EF_r = E0(F, 1:nF);
        XF   = X(F, :);
        C    = setdiff(1:n, F);
        Xr(C, :) = X(C, :) + E_r(C, :) * XF;
        Xr(F, :) = EF_r * XF;
    end
    vectorRecoveryTime = seconds(datetime('now') - t_rec);

    sKeep = min(s, size(Xr, 2));
    Xr    = Xr(:, 1:sKeep);
    D     = D(1:sKeep);
    X     = Xr;

    e2eTime = seconds(datetime('now') - t) - vectorRecoveryTime - simTime - diagTime;

    fprintf('\nsolvTime \t e2eTime \t eCodingGenTime \t augTime \t vecRecovTime \t restructTime\n');
    fprintf('%f \t %f \t %f \t %f \t %f \t %f\n', ...
        solvTime, e2eTime, eCodingGenTime, augTime, vectorRecoveryTime, restructTime);

    stats = struct();
    stats.method          = 'inverse_power_shifted';
    stats.rrMode          = rrMode;
    stats.xRecovery       = p.Results.xRecovery;
    stats.nRecovered      = nRecovered;
    stats.nGramFallback   = nGramFallback;
    stats.sigma           = 'smallestabs';
    stats.iters           = iters;
    stats.res             = relres;
    stats.eigenvalues     = D;
    stats.converged       = ~isempty(relres) && relres(end) < tol;
    stats.stagnated       = stagnatedFlag;
    stats.n_iters         = numel(iters);
    stats.faultIterations = errorIterations;
    stats.faultRows       = F;
    stats.normM = [];  stats.normMinv = [];
    if ~isempty(F)
        [stats.normM, stats.normMinv] = pencilCongruenceNorms(E0(:, 1:numel(F)), F);
    end
    stats.r0         = r0;
    stats.blockWidth = m;
    stats.nRefactor  = nRefactor;
    stats.shiftMode  = shiftMode;
    stats.timing = struct( ...
        'eCodingGenTime',     eCodingGenTime, ...
        'augTime',            augTime, ...
        'restructTime',       restructTime, ...
        'solvTime',           solvTime, ...
        'vectorRecoveryTime', vectorRecoveryTime, ...
        'simTime',            simTime, ...
        'diagTime',           diagTime, ...
        'e2eTime',            e2eTime);
end


function [solve, ok, sigmaUsed] = ipmShiftedSolver(A, sigma, n, E0, R0, S0, T0, F, C, sigmaLo)
    solve = []; ok = false; sigmaUsed = sigma;
    if nargin < 10 || isempty(sigmaLo)
        maxTry  = 1;
        sigmaLo = sigma;
    else
        maxTry = 20;
    end
    cand = sigma;
    for attempt = 1:maxTry
        [solve, ok] = tryBuildShift(A, cand, n, E0, R0, S0, T0, F, C);
        if ok
            sigmaUsed = cand;
            return;
        end
        cand = (sigmaLo + cand) / 2;
        if abs(cand - sigmaLo) <= eps(max(abs(sigmaLo), 1))
            break;
        end
    end
    sigmaUsed = sigma;
end

function [R, p, S] = cholAuto(M)
    if issparse(M)
        [R, p, S] = chol(M);
    else
        [R, p] = chol(M);
        S = speye(size(M, 1));
    end
    if p ~= 0
        R = [];
        S = speye(size(M, 1));
    end
end

function [solve, ok] = tryBuildShift(A, sigma, n, E0, R0, S0, T0, F, C)
    solve = [];
    Ashift = A - sigma * speye(n);
    Ashift = (Ashift + Ashift') / 2;

    if isempty(F)
        [Rc, p, Sc] = cholAuto(Ashift);
        ok = (p == 0);
        if ~ok, return; end
        solve = @(V) Sc * (Rc \ (Rc' \ (Sc' * V)));
        return;
    end

    nF      = numel(F);
    E       = E0(:, 1:nF);
    R_shift = R0(:, 1:nF) - sigma * E;
    S_shift = S0(1:nF, 1:nF) - sigma * T0(1:nF, 1:nF);

    [Rc, p, Sc] = cholAuto(Ashift(C, C));
    ok = (p == 0);
    if ~ok, return; end
    solveCC = @(V) Sc * (Rc \ (Rc' \ (Sc' * V)));

    R_CF     = full(R_shift(C, :));
    AinvR_CF = solveCC(R_CF);
    SC       = full((S_shift + S_shift') / 2 - R_CF' * AinvR_CF);
    SC       = (SC + SC') / 2;
    if ~all(isfinite(SC(:)))
        ok = false; return;
    end
    [L_SC, pS] = chol(SC, 'lower');
    ok = (pS == 0);
    if ~ok, return; end
    solve = @(V) ipmSchurSolve(solveCC, AinvR_CF, L_SC, R_CF, V, F, C, n);
end

function Z = ipmSchurSolve(solveCC, AinvR_CF, L_SC, R_CF, X, F, C, n_full)
    AinvX_C = solveCC(X(C, :));
    RHS_F   = X(F, :) - R_CF' * AinvX_C;
    Z_F     = L_SC' \ (L_SC \ RHS_F);
    Z = zeros(n_full, size(X, 2));
    Z(C, :) = AinvX_C - AinvR_CF * Z_F;
    Z(F, :) = Z_F;
end

function Y = Aop(A, R, S, X, F)
    Y = A * X;
    if ~isempty(F)
        XF     = X(F, :);
        Y      = Y      + R * XF;
        Y(F,:) = Y(F,:) + R' * X;
        Y(F,:) = Y(F,:) + (S - R(F,:) - R(F,:)') * XF;
    end
end

function Y = Bop(E, T, X, F)
    Y = X;
    if ~isempty(F)
        XF     = X(F, :);
        Y(F,:) = 0;
        Y      = Y      + E * XF;
        Y(F,:) = Y(F,:) + E' * X;
        Y(F,:) = Y(F,:) + (T - E(F,:) - E(F,:)') * XF;
    end
end
