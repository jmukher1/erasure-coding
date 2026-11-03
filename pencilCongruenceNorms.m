function [normM, normMinv] = pencilCongruenceNorms(E, F)
% pencilCongruenceNorms - ||M||_2 and ||M^{-1}||_2 of the congruence
% P = [I(:,C), E] (C = complement of F) that maps the coded pencil to the
% original one, computed from the k x k matrix E'*E. Not part of any timing.
    F  = F(:);
    k  = numel(F);
    E  = full(E(:, 1:k));
    n  = size(E, 1);
    EF = E(F, :);
    T  = E' * E;
    Gc = T - EF' * EF;
    Gc = (Gc + Gc') / 2;
    [V, Dg] = eig(Gc);
    R  = diag(sqrt(max(diag(Dg), 0))) * V';
    K  = [eye(k), R; R', (T + T') / 2];
    lam = eig((K + K') / 2);
    if n - k > k
        lam = [lam; 1];
    end
    normM    = sqrt(max(lam));
    normMinv = 1 / sqrt(max(min(lam), realmin));
end
