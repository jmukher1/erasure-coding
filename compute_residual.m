function [resNorm, relresNorm, absresNorm] = compute_residual(AX, BX, Lambda, refNorm, A0norm)
% compute_residual - resNorm    = ||A*X - B*X*Lambda||_F
%                     relresNorm = resNorm / refNorm   (refNorm = initial residual r0)
%                     absresNorm = resNorm / A0norm 

    R = AX - BX * Lambda;
    resNorm = norm(R, 'fro');

    relresNorm = resNorm / max(refNorm, realmin);
    absresNorm = resNorm / max(A0norm, realmin);
end
