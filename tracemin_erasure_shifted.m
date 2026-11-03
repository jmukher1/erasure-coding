function [X, Thi, iters, res, relres, errorIterations, eCodingGenTime, augTime, solvTime, restructTime, ...
    e2eTime, stats] = tracemin_erasure_shifted(A, X0, ...
    frac, s, possibleErrorIterations, rowDeletedStart, codingScheme, SIGMA, ...
    perturbationScheme, perturbationParameter, TOL, maxiter, opts, sopts)
% tracemin_erasure_shifted - erasure-coded TraceMin for s eigenpairs of
% (A, I) at the smallest end (SIGMA = 'smallestabs'), with timing.
%
 
if nargin < 14, sopts = struct(); end
if nargin < 12 || nargin > 14
    error('tracemin_erasure_shifted:nargin', 'expected 12 to 14 input arguments.');
end

innerSolve = 'minres';
if nargin >= 13 && isstruct(opts) && isfield(opts, 'innerSolve')
    innerSolve = opts.innerSolve;
    opts = rmfield(opts, 'innerSolve');
end
if nargin < 13 || isempty(opts)
    opts = tracemin_opts();
else
    opts = tracemin_opts(opts);
end
opts.innerSolve = innerSolve;
if isempty(maxiter)
    maxiter = 1000;
end

t = datetime("now");                              % start of e2eTime
[n,~] = size(A);
B = speye(n);
k = ceil(frac*n);

if perturbationScheme == PerturbationScheme.alpha
    A = A + perturbationParameter .* A;
    B = B + perturbationParameter .* A;
elseif perturbationScheme == PerturbationScheme.mu
    A = A + perturbationParameter .* B;
    B = B + perturbationParameter .* B;
end

% coding matrix E (n x k)
t01 = datetime("now");
E = makeErasureCodingMatrix(A, n, k, opts.coding, codingScheme);
t1 = datetime("now");
eCodingGenTime = seconds(t1 - t01);

% coding blocks A*E, E'*A*E, B*E, E'*B*E
[Z, R, Q, S] = generateErasureCodingBlocks(A, B, E);
t3 = datetime("now");
augTime = seconds(t3 - t1);

% TraceMin iterations (fault handling is timed inside)
normA = norm(A, 'fro');
t4 = datetime("now");
[X, Thi, iters, res, relres, errorIterations, ~, decomp2Time, restructTime, tmStats, Xtrue] = tracemin_shifted(A, B, X0, normA, Z, R, Q, S, k, s, ...
            possibleErrorIterations, rowDeletedStart, TOL, maxiter, opts, E, SIGMA, sopts);
t4s = datetime("now");

diagTime = 0; simTime = 0;
if isfield(tmStats, 'timing') && isfield(tmStats.timing, 'diagTime')
    diagTime = tmStats.timing.diagTime;
end
if isfield(tmStats, 'timing') && isfield(tmStats.timing, 'simTime')
    simTime = tmStats.timing.simTime;
end
solvTime = max(seconds(t4s - t4) - restructTime - diagTime - simTime, 0);

if perturbationScheme == PerturbationScheme.alpha
    Thi = Thi ./ (1 + perturbationParameter);
elseif perturbationScheme == PerturbationScheme.mu
    Thi = Thi - (perturbationParameter .* eye(s));
end

t5 = datetime("now");
e2eTime = seconds(t5 - t) - diagTime - simTime;
fprintf('\nsolvTime (s), e2eTime (s), eCodingGenTime (s), augTime (s), restructTime (s)\n');
fprintf('\n%f, %f, %f, %f, %f\n', solvTime, e2eTime, eCodingGenTime, augTime, restructTime);

stats = struct();
stats.method          = 'tracemin';
stats.sigma           = SIGMA;
stats.iters           = iters;
stats.res             = relres;
stats.relres          = relres;
stats.absres          = res;
stats.eigenvalues     = diag(Thi);
stats.converged       = tmStats.converged;
stats.n_iters         = numel(iters);
stats.faultIterations = errorIterations;
if isfield(tmStats, 'r0')
    stats.r0 = tmStats.r0;
end
for fwd = {'blockWidth','innerIters','innerItersTotal','nShifted','nLocked', ...
           'nDeflated','perColRes','opts','nRefactor','sigmaFinal','stagnated', ...
           'codedres','sigmaHist','innerSolve', ...
           'setupTime','iterWall','postFault','normM','normMinv'}
    if isfield(tmStats, fwd{1})
        stats.(fwd{1}) = tmStats.(fwd{1});
    end
end
stats.Xtrue = Xtrue;
if isfield(tmStats, 'timing') && isfield(tmStats.timing, 'nFullRefactor')
    stats.nFullRefactor = tmStats.timing.nFullRefactor;
end
stats.timing = struct( ...
    'eCodingGenTime', eCodingGenTime, ...
    'augTime',        augTime, ...
    'restructTime',   restructTime, ...
    'refactorTime',   decomp2Time, ...           % part of restructTime
    'diagTime',       diagTime, ...
    'simTime',        simTime, ...
    'solvTime',       solvTime, ...
    'e2eTime',        e2eTime);
end
