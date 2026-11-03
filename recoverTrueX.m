function Xr = recoverTrueX(X, E, F)
% recoverTrueX - map a block X of the coded pencil back to the original
% problem, Xr = P*X with P = I except P(:,F) = E:
%   Xr(F,:) = E(F,:)*X(F,:),  Xr(C,:) = X(C,:) + E(C,:)*X(F,:)

    if isempty(F) || isempty(E)
        Xr = X;
        return;
    end

    nF = numel(F);
    if size(E, 2) > nF
        E = E(:, 1:nF);
    elseif size(E, 2) < nF
        error('recoverTrueX:tooFewCodingColumns', ...
              'E has %d columns but |F| = %d.', size(E, 2), nF);
    end

    Xr       = X;
    Xr(F, :) = 0;
    Xr       = Xr + E * X(F, :);
end
