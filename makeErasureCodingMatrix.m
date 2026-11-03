function E = makeErasureCodingMatrix(A, n, k, coding, codingScheme)
% makeErasureCodingMatrix - n x k erasure-coding matrix E.
%
% 'default': staggered sparse matrix (getErasureCodedMatrix).
% 'interleaved': row i has the single entry sqrt(k/n) in column mod(i-1,k)+1,
% so every contiguous block of k rows of E is a scaled permutation matrix
% (E_F is perfectly conditioned wherever the erased window lands) and
% ||E||_2 is about 1.
    if nargin < 4 || isempty(coding),       coding = 'interleaved'; end
    if nargin < 5 || isempty(codingScheme), codingScheme = 'Random'; end
    if ~strcmpi(string(codingScheme), "Random")
        error('makeErasureCodingMatrix:badScheme', 'only codingScheme ''Random'' is supported.');
    end
    switch lower(coding)
        case 'default'
            E = getErasureCodedMatrix(n, k);
        case 'interleaved'
            r = (1:n)';
            E = sparse(r, mod(r-1, k) + 1, sqrt(k / n), n, k);
        otherwise
            error('makeErasureCodingMatrix:badCoding', ...
                  'coding must be ''default'' or ''interleaved'' (got ''%s'').', coding);
    end
end
