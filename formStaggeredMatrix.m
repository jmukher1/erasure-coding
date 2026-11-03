function E = formStaggeredMatrix(ED, k)
% formStaggeredMatrix - place the p entries of row i of ED in columns
% mod(i-1, k)+1, ..., mod(i+p-2, k)+1 of the n x k matrix E.
    [n, p] = size(ED);
    E = zeros(n, k);
    for i = 1 : n
        js = 1 + mod(i-1, k);
        for jr = 1 : p
            j = 1 + mod(js + jr - 2, k);
            E(i, j) = ED(i, jr);
        end
    end
    E = sparse(E);
end
