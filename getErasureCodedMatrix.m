function E = getErasureCodedMatrix(n, k)
% getErasureCodedMatrix - sparse staggered erasure-coding matrix E (n x k).
    p = min(k, max(6, ceil(2 * log(max(k, 2)))));
    E_seed = randn(n, p) * sqrt(k / (n * p));
    E = formStaggeredMatrix(E_seed, k);
end
