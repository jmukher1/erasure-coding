# erasure-coding

MATLAB code for "Fault Oblivious Eigenvalue Solver".

## Setup

Put the matrices in `data/` (`<name>.mtx` sparse, `<name>.mat` dense), then from
the `erasure-coding` folder in MATLAB:

```matlab
addpath(pwd);  addpath(fullfile(pwd, 'driver'));

A = loadSparse('Pres_Poisson');  A = (A + A') / 2;
n = size(A, 1);  s = 10;              % s eigenpairs
frac = 0.01;  k = ceil(frac * n);     % erase 1% of the rows
rng(0);  X0 = orth(randn(n, s));      % starting block
fault = 10;                           % iteration of the fault
row   = 100;                          % erased rows: row .. row+k-1
```

Without erasure, pass `[]` as the fault iteration and `-1` as the first row.

## TraceMin (smallest eigenvalues)

```matlab
opts = struct('innerSolve', 'schur');

% no erasure
[X, Thi, iters, ~, relres, ~, ~, ~, ~, ~, ~, stats] = tracemin_erasure_shifted( ...
    A, X0, frac, s, [], -1, 'Random', 'smallestabs', 'None', 'None', 1e-10, 1500, opts);

% erasure coding, one fault
[X, Thi, iters, ~, relres, ~, ~, ~, ~, ~, ~, stats] = tracemin_erasure_shifted( ...
    A, X0, frac, s, fault, row, 'Random', 'smallestabs', 'None', 'None', 1e-10, 1500, opts);

lambda = diag(Thi);
```

## Inverse power (smallest eigenvalues)

```matlab
% no erasure
[~, D, Xr, ~, iters, relres, ~, ~, ~, ~, ~, ~, ~, stats] = inverse_power_shifted( ...
    A, X0, k, [], -1, 'tol', 1e-10, 'maxiter', 1500);

% erasure coding, one fault
[~, D, Xr, ~, iters, relres, ~, ~, ~, ~, ~, ~, ~, stats] = inverse_power_shifted( ...
    A, X0, k, fault, row, 'tol', 1e-10, 'maxiter', 1500);
```

## Block power (largest eigenvalues)

```matlab
% no erasure
[~, D, Xr, ~, iters, relres, ~, ~, ~, ~, ~, ~, ~, stats] = gen_erasure_power_shifted( ...
    A, X0, k, [], -1, 'tol', 1e-10, 'maxiter', 1500);

% erasure coding, one fault
[~, D, Xr, ~, iters, relres, ~, ~, ~, ~, ~, ~, ~, stats] = gen_erasure_power_shifted( ...
    A, X0, k, fault, row, 'tol', 1e-10, 'maxiter', 1500);
```

`D` holds the eigenvalues and `Xr` the eigenvectors. `relres` is the residual
history relative to the first iteration.

## Timing

Every solver returns `stats.timing` (seconds):

| Field | Measures |
|---|---|
| `eCodingGenTime` | building the coding matrix E |
| `augTime` | building the coding blocks A*E, E'*A*E |
| `restructTime` | recovering from the fault |
| `solvTime` | the iterations |
| `e2eTime` | the whole call |

 