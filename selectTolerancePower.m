function TOL = selectTolerancePower(SIGMA, dataset)
% stopping tolerance of the inverse power ('smallestabs') and block power ('largestabs')  
    if strcmp(SIGMA, 'largestabs')
        switch dataset
            case 'bcsstk18',  TOL = 7e-10;
            otherwise,        TOL = 1e-10;
        end
    else
        switch dataset
            case 'bcsstk18'
                TOL = 1e-11;
            case {'finan512', 'wathen120'}
                TOL = 1e-12;
            case 'MNIST_SPD_TrainImages'
                TOL = 8e-10;
            case {'bodyy4', 'bodyy5', 'bodyy6', 'bundle1', 'Pres_Poisson', 'wathen100', ...
                  'CIFAR10_SPD_TrainImages', 'COCO_SPD_TrainImages'}
                TOL = 1e-10;
            otherwise
                TOL = 1e-14;
        end
    end
end
