function A = loadSparse(dataset)
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
    [A, ~, ~, ~, ~, ~, ~] = mmread(fullfile(dataDir, [dataset '.mtx']));
end
