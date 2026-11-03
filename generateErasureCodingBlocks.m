function [Z, R, Q, S] = generateErasureCodingBlocks(A, B, E)
% generateErasureCodingBlocks - coding blocks of the pencil (A, B):
%   Z = A*E,  R = E'*A*E,  Q = B*E,  S = E'*B*E

    Z = A * E;
    R = E'* Z;
    Q = B * E;
    S = E' * Q;

end
