function [X, D, Xr, F, iters, relres, errorIterations, eCodingGenTime, augTime, solvTime, vectorRecoveryTime, restructTime, e2eTime, stats] = ...
    gen_erasure_power_shifted(A, X0, k, possibleErrorIterations, rowDeletedStart, varargin)
% erasure-coded block power method (largest eigenvalues) with an optional shift
%   A                       n x n symmetric matrix
%   X0                      n x s starting block (padded to m = blockFactor*s)
%   k                       number of erased rows per fault (= columns of E)
%   possibleErrorIterations iterations at which a fault is injected ([] = none)
%   rowDeletedStart         first erased row (rows rowDeletedStart..+k-1)
%   name/value: tol, maxiter, blockFactor, coding, codingScheme, shiftMode
%               ('auto' | 'none'), shiftFrac, faultReset ('random' | 'zero'),
%               padSeed, reportSpike, postFault ('reconstruct' | 'coded'),
%               reconDrop, xRecovery ('local' | 'none')
%
% postFault 'reconstruct': after a fault the erased rows/cols of A are rebuilt
% from R = A*E and the iteration continues on (A, I); with xRecovery 'local'
% the erased rows of the iterate are rebuilt from the eigen-equation.
% postFault 'coded': the iteration continues on the coded pencil (A', B').
%
% Timing (wall clock, seconds)
%   eCodingGenTime      build E
%   augTime             build R = A*E, S = E'*A*E, T = E'*E
%   restructTime        fault handling (zeroing excluded, see simTime)
%   vectorRecoveryTime  Xr = P*X after the loop
%   solvTime            iteration loop minus restructTime
%   e2eTime             whole call
% simTime (zeroing the erased rows/cols, i.e. simulating the fault)  

    t = datetime("now");

    p = inputParser;
    addParameter(p, 'tol',          1e-8);
    addParameter(p, 'maxiter',      2000);
    addParameter(p, 'blockFactor',  2);
    addParameter(p, 'coding',       'interleaved');
    addParameter(p, 'codingScheme', 'Random');
    addParameter(p, 'shiftMode',    'auto');
    addParameter(p, 'shiftFrac',    0.5);
    addParameter(p, 'faultReset',   'random');
    addParameter(p, 'padSeed',      42);
    addParameter(p, 'reportSpike',  true);
    addParameter(p, 'postFault',    'reconstruct');
    addParameter(p, 'reconDrop',    100);
    addParameter(p, 'xRecovery',    'local');
    parse(p, varargin{:});
    postFault = p.Results.postFault;
    if ~ismember(postFault, {'reconstruct', 'coded'})
        error('gen_erasure_power_shifted:badPostFault', 'postFault must be ''reconstruct'' or ''coded''.');
    end
    if ~ismember(p.Results.xRecovery, {'local', 'none'})
        error('gen_erasure_power_shifted:badXRecovery', 'xRecovery must be ''local'' or ''none''.');
    end

    tol         = p.Results.tol;
    maxiter     = p.Results.maxiter;
    blockFactor = p.Results.blockFactor;
    shiftMode   = p.Results.shiftMode;
    shiftFrac   = p.Results.shiftFrac;
    if ~ismember(shiftMode, {'auto', 'none'})
        error('gen_erasure_power_shifted:badShiftMode', 'shiftMode must be ''auto'' or ''none''.');
    end
    if blockFactor < 1
        error('gen_erasure_power_shifted:badBlockFactor', 'blockFactor must be >= 1.');
    end

    sigma    = 0;
    sigHist  = [];
    nReconstruct = 0;
    diagTime = 0;
    simTime  = 0;
    [n, s] = size(X0);

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
        fprintf('block power: oversampling to m = %d columns (= %g * s, s = %d)\n', ...
                m, blockFactor, s);
    end

    X = X0;
    eCodingGenTime = 0; augTime = 0; restructTime = 0;
    solvTime = 0; e2eTime = 0; vectorRecoveryTime = 0; r0 = 1;

    % coding matrix E (n x k)
    t01 = datetime("now");
    E0 = makeErasureCodingMatrix(A, n, k, p.Results.coding, p.Results.codingScheme);
    t1  = datetime("now");
    eCodingGenTime = seconds(t1 - t01);

    % coding blocks R = A*E, S = E'*A*E, T = E'*E
    R0 = A * E0;
    S0 = E0' * R0;
    T0 = E0' * E0;
    t3 = datetime("now");
    augTime = seconds(t3 - t1);

    F               = [];
    EF              = [];
    G               = [];
    relres          = [];
    iters           = [];
    errorIterations = [];
    normA = norm(A, 'fro');
    if strcmp(postFault, 'reconstruct'), A_orig = A; end
    E = []; R = []; S = []; T = [];
    AXprev = [];

    fprintf('iter\tResidual Norm\tRelative Residual Norm\n');
    t4 = datetime("now");                       % start of the iteration loop

    for iter = 1:maxiter

        % fault injection
        Fnew = faultfunc_matlab2(iter, possibleErrorIterations, ...
                                 rowDeletedStart, rowDeletedStart + k - 1);

        if ~isempty(Fnew)
            fprintf('\nnew faults detected\n');
            fprintf('Set %d as errorIterations(end+1) \n', iter);
            errorIterations(end+1) = iter;
            t11 = datetime("now");
            diagAtFault = diagTime;
            simAtFault  = simTime;

            Fprospective = sort(union(F, Fnew));
            assert(length(Fprospective) <= k, ...
                'gen_erasure_power_shifted:tooManyFaults', ...
                'Accumulated faults (%d) exceed tolerated k (%d).', length(Fprospective), k);

            % simulate the fault: reset the lost iterate rows, zero rows/cols F of A
            tsim = tic;
            XFold = X(Fnew, :);
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
            simTime = simTime + toc(tsim);

            F  = Fprospective;
            nF = length(F);
            EF = E0(F, :);

            % coded pencil: factor of E_F'*E_F, used by invBop
            if strcmp(postFault, 'coded')
                AXprev = [];
                EF_nF  = EF(:, 1:nF);
                GramEF = EF_nF' * EF_nF;
                condEst = condest(GramEF);
                regTol  = 1e8;
                if condEst > regTol
                    lambda = regTol * eps * norm(GramEF, 'fro');
                    GramEF = GramEF + lambda * speye(size(GramEF, 1));
                    fprintf('  [invBop regularization] nF=%d, cond=%.3e -> ridge lambda=%.3e\n', ...
                        nF, condEst, lambda);
                end
                try
                    G = chol(GramEF, 'upper');
                    G_isChol = true;
                catch
                    G = decomposition((GramEF + GramEF') / 2, 'ldl');
                    G_isChol = false;
                end
            end

            % residual of the damaged iterate, reported for this iteration
            % (plot only; timed as diagTime)
            spikeAbs = []; spikeRel = [];
            if p.Results.reportSpike && exist('D', 'var') && iter > 1
                td = tic;
                Xd = X(:, 1:s) ./ max(vecnorm(X(:, 1:s)), realmin);
                if strcmp(postFault, 'reconstruct')
                    [spikeAbs, spikeRel] = compute_residual(A_orig * Xd, Xd, diag(D(1:s)), r0, normA);
                else
                    AXd = Aop(A, R0(:, 1:nF), S0(1:nF, 1:nF), Xd, F);
                    BXd = Bop(E0(:, 1:nF), T0(1:nF, 1:nF), Xd, F);
                    [spikeAbs, spikeRel] = compute_residual(AXd, BXd, diag(D(1:s)), r0, normA);
                end
                diagTime = diagTime + toc(td);
            end

            % rebuild rows/cols F of A from R = A*E, then the erased rows of
            % the iterate from the eigen-equation; continue on (A, I)
            if strcmp(postFault, 'reconstruct')
                Cr   = setdiff((1:n)', F);
                EFb  = full(E0(F, 1:nF));
                ECb  = E0(Cr, 1:nF);
                RCb  = R0(Cr, 1:nF);
                RFb  = R0(F, 1:nF);
                if issparse(A)
                    iEF  = inv(EFb);
                    ACC  = A(Cr, Cr);
                    ACF  = full(RCb - ACC * ECb) / EFb;
                    nzC  = eps * (full(abs(ACC) * abs(ECb) + abs(RCb)) * abs(iEF));
                    ACF(abs(ACF) < p.Results.reconDrop * nzC) = 0;
                    AFF  = full(RFb - ACF' * ECb) / EFb;
                    nzF  = eps * (full(abs(ACF') * abs(ECb) + abs(RFb)) * abs(iEF));
                    nzF  = max(nzF, nzF');
                    AFF  = (AFF + AFF') / 2;
                    AFF(abs(AFF) < p.Results.reconDrop * nzF) = 0;
                    A(Cr, F) = sparse(ACF);  A(F, Cr) = sparse(ACF');  A(F, F) = sparse(AFF);
                else
                    AE   = A * E0(:, 1:nF);
                    ACF  = (RCb - AE(Cr, :)) / EFb;
                    AFF  = (RFb - ACF' * ECb) / EFb;
                    AFF  = (AFF + AFF') / 2;
                    A(Cr, F) = ACF;          A(F, Cr) = ACF';          A(F, F) = AFF;
                end
                nReconstruct = nReconstruct + 1;
                fprintf('  [reconstruct] A rebuilt from R = A*E (rows %d..%d)\n', F(1), F(end));

                nRec = 0;
                if strcmp(p.Results.xRecovery, 'local') && exist('D', 'var') && numel(D) >= size(X, 2)
                    AFFf = full(AFF);  nRec = 0;
                    for j = 1:size(X, 2)
                        Mj = D(j) * eye(nF) - AFFf;
                        if rcond(Mj) > 1e-12
                            xF = Mj \ (ACF' * X(Cr, j));
                            if all(isfinite(xF)), X(F, j) = xF; nRec = nRec + 1; end
                        end
                    end
                    fprintf('  [xRecovery] rebuilt erased rows of %d/%d iterate columns\n', nRec, size(X, 2));
                end

                if ~isempty(AXprev) && nRec == size(X, 2)
                    AXprev = AXprev + A(:, F) * (X(F, :) - XFold);
                else
                    AXprev = [];
                end
                F = [];  EF = [];  G = [];
            end

            restructTime = restructTime + seconds(datetime("now") - t11) ...
                           - (diagTime - diagAtFault) - (simTime - simAtFault);
        end

        if ~isempty(F)
            nF = length(F);
            E  = E0(:, 1:nF);
            R  = R0(:, 1:nF);
            S  = S0(1:nF, 1:nF);
            T  = T0(1:nF, 1:nF);
            EF = EF(:, 1:nF);
        else
            E = []; R = []; S = []; T = []; EF = [];
        end

        % shifted power step
        if isempty(AXprev)
            Y = Aop(A, R, S, X, F);
        else
            Y = AXprev;
        end
        if isempty(F)
            Z = Y;
        else
            Z = invBop(G, G_isChol, E, EF, Y, F);
        end
        if sigma ~= 0
            Z = Z - sigma * X;
        end
        sigHist(end+1) = sigma;

        % Rayleigh-Ritz on (A', B')
        [Q, ~] = qr(Z, 0);
        AQ_mat = Aop(A, R, S, Q, F);
        BQ_mat = Bop(E, T, Q, F);
        AQ = Q' * AQ_mat;   AQ = (AQ + AQ') / 2;
        BQ = Q' * BQ_mat;   BQ = (BQ + BQ') / 2;
        [U, Dmat] = eig(AQ, BQ, 'chol');
        [D, idx] = sort(real(diag(Dmat)), 'descend');
        U = U(:, idx);
        X = Q * U;
        AX = AQ_mat * U;
        BX = BQ_mat * U;
        AXprev = AX;

        % shift for the next step
        if strcmp(shiftMode, 'auto') && D(end) > 0
            sigma = shiftFrac * D(end);
        else
            sigma = 0;
        end

        % residual of the current pencil, relative to r0
        AXres = AX(:, 1:s);
        BXres = BX(:, 1:s);
        if iter == 1
            r0 = norm(AXres - BXres .* D(1:s)', 'fro');
        end
        [resNorm, relresNorm] = compute_residual(AXres, BXres, diag(D(1:s)), r0, normA);
        if exist('spikeRel', 'var') && ~isempty(spikeRel) && ~isempty(errorIterations) && errorIterations(end) == iter
            fprintf('  [reportSpike] fault iteration reports the damaged-iterate residual %.3e (recovered: %.3e)\n', ...
                    spikeRel, relresNorm);
            resNorm = spikeAbs;  relresNorm = spikeRel;
            spikeRel = [];
        end
        relres(end+1) = relresNorm;
        iters(end+1)  = iter;

        fprintf('%d\t%e\t%e\t%e\n', iter, resNorm, relresNorm, tol);
        if relresNorm < tol
            fprintf('Converged in %d iterations (computable rel. residual = %e)\n', iter, relresNorm);
            break;
        end
    end

    % eigenvectors of the original problem 
    t_rec = datetime("now");
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
    vectorRecoveryTime = seconds(datetime("now") - t_rec);

    sKeep = min(s, size(Xr, 2));
    Xr    = Xr(:, 1:sKeep);
    X     = X(:,  1:sKeep);
    D     = D(1:sKeep);

    t5 = datetime("now");
    e2eTime  = seconds(t5 - t) - diagTime - simTime - vectorRecoveryTime;
    solvTime = seconds(t5 - t4) - restructTime - vectorRecoveryTime - diagTime - simTime;

    fprintf('\nsolvTime \t e2eTime \t eCodingGenTime\t augTime \t vecRecovTime \t restructTime\n');
    fprintf('%f \t %f \t %f \t %f \t %f \t %f\n', ...
        solvTime, e2eTime, eCodingGenTime, augTime, vectorRecoveryTime, restructTime);

    stats = struct();
    stats.method          = 'block_power_shifted';
    stats.shiftMode       = shiftMode;
    stats.sigmaHist       = sigHist;
    stats.postFault       = postFault;
    stats.nReconstruct    = nReconstruct;
    stats.sigma           = 'largestabs';
    stats.iters           = iters;
    stats.res             = relres;
    stats.eigenvalues     = D;
    stats.converged       = ~isempty(relres) && relres(end) < tol;
    stats.stagnated       = ~isempty(relres) && ~stats.converged && numel(relres) < maxiter;
    stats.n_iters         = numel(iters);
    stats.faultIterations = errorIterations;
    stats.faultRows       = F;
    % ||M||_2, ||M^{-1}||_2 of the congruence used (computed after timing stops)
    stats.normM = [];  stats.normMinv = [];
    if ~isempty(F)
        [stats.normM, stats.normMinv] = pencilCongruenceNorms(E0(:, 1:numel(F)), F);
    end
    stats.r0 = r0;
    stats.timing = struct( ...
        'eCodingGenTime',     eCodingGenTime, ...
        'augTime',            augTime, ...
        'restructTime',       restructTime, ...
        'solvTime',           solvTime, ...
        'vectorRecoveryTime', vectorRecoveryTime, ...
        'diagTime',           diagTime, ...
        'simTime',            simTime, ...
        'e2eTime',            e2eTime);
end


% A' * X for the coded pencil (rows/cols F of A are zero)
function Y = Aop(A, R, S, X, F)
    Y = A * X;
    if ~isempty(F)
        XF     = X(F, :);
        Y      = Y      + R * XF;
        Y(F,:) = Y(F,:) + R' * X;
        Y(F,:) = Y(F,:) + (S - R(F,:) - R(F,:)') * XF;
    end
end

% B' * X, B' = P'P
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

function Y = invBop(G, G_isChol, E, EF, X, F)
    Y      = X;
    Y(F,:) = 0;
    W = E' * X;
    V = X(F,:) + EF' * X(F,:);
    if G_isChol
        GinvW   = G \ (G' \ W);
        GinvV   = G \ (G' \ V);
        GinvVmW = G \ (G' \ (V - W));
    else
        GinvW   = G \ W;
        GinvV   = G \ V;
        GinvVmW = G \ (V - W);
    end
    Y      = Y + E * GinvW - E * GinvV;
    I_F    = eye(length(F));
    Y(F,:) = Y(F,:) + (I_F + EF) * GinvVmW;
end
