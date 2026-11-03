function [V, nDeflated, BV] = BOrthonormalizeStable(B, Z, dropTol)
% B-orthonormalize the columns of Z (V'*B*V = I). 

    if nargin < 3 || isempty(dropTol)
        dropTol = 1e-12;
    end

    [n, m] = size(Z);
    nDeflated = 0;

    cn = sqrt(sum(Z .* Z, 1));
    cn(~isfinite(cn) | cn == 0) = 1;
    V = Z ./ cn;

    for pass = 1:2
        BVp = B * V;
        G = V' * BVp;
        G = (G + G') / 2;

        [Q, D] = eig(G);
        d = real(diag(D));
        dmax = max(d);

        if ~isfinite(dmax) || dmax <= 0
            error('BOrthonormalizeStable:rankLoss', ...
                  'Gram matrix has no positive eigenvalues; the block has collapsed entirely.');
        end

        keep = d > dropTol * dmax;
        if pass == 1
            nDeflated = sum(~keep);
        end

        C  = Q(:, keep) ./ sqrt(d(keep))';
        V  = V * C;
        BV = BVp * C;
    end

    guard = 0;
    while size(V, 2) < m && guard < 5
        guard = guard + 1;
        W = randn(n, m - size(V, 2));

        W = W - V * (V' * (B * W));
        W = W - V * (V' * (B * W));

        Gw = W' * (B * W);
        Gw = (Gw + Gw') / 2;
        [Qw, Dw] = eig(Gw);
        dw = real(diag(Dw));
        if isempty(dw) || max(dw) <= 0
            break;
        end
        kw = dw > dropTol * max(dw);
        if ~any(kw)
            break;
        end
        V = [V, W * (Qw(:, kw) ./ sqrt(dw(kw))')];
    end
    if guard > 0
        BV = B * V;
    end
end
