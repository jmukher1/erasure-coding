function TOL = selectTolerance(SIGMA, dataset)
% selectTolerance - TraceMin stopping tolerance (smallest eigenvalues) per
% dataset. TOL is a reduction factor on the initial residual r0 (see compute_residual).
    switch dataset
        case {'bcsstk18', 'finan512', 'wathen120'}
            TOL = 1e-12;
        case 'bodyy6'
            TOL = 1e-9;
        case {'bodyy4', 'bodyy5', 'bundle1', 'Pres_Poisson', 'wathen100', ...
              'CIFAR10_SPD_TrainImages', 'COCO_SPD_TrainImages', 'MNIST_SPD_TrainImages'}
            TOL = 1e-10;
        otherwise
            TOL = 1e-14;
    end
end
