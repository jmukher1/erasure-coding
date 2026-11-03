function A = loadDense(dataset)
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
    M = load(fullfile(dataDir, [dataset '.mat']));
    if contains(dataset, "Test")
        A = M.TestCov;
    elseif contains(dataset, "Train")
        A = M.TrainCov;
    elseif contains(dataset, "Validation")
        A = M.ValCov;
    end
end
