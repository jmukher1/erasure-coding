function varargout = cr_store(action, varargin)
% stable-storage layer of the checkpoint/restart baselines. 

    switch lower(action)
        case 'init',  varargout{1} = doInit(varargin{:});
        case 'write', [varargout{1}, varargout{2}] = doWrite(varargin{:});
        case 'reada', [varargout{1}, varargout{2}] = doReadA(varargin{:});
        case 'readz', [varargout{1}, varargout{2}] = doReadZ(varargin{:});
        case 'close', doClose(varargin{:});
        case 'savea'
            [A, file, ver] = varargin{:};
            t0 = tic;  saveStruct(file, struct('A', A), ver);
            varargout{1} = file;  varargout{2} = toc(t0);
        otherwise, error('cr_store:badAction', 'Unknown action %s.', action);
    end
end

function ctx = doInit(A, cr)
    parent = cr.ckpt_dir;
    if isempty(parent), parent = tempdir; end
    if ~exist(parent, 'dir'), mkdir(parent); end
    [~, leaf] = fileparts(tempname);
    ctx.dir   = fullfile(parent, ['cr_' leaf]);
    mkdir(ctx.dir);
    ctx.fileZ = fullfile(ctx.dir, 'ckpt.mat');
    ctx.ver   = cr.mat_version;
    ctx.keep  = cr.keep_files;
    if isfield(cr, 'a_file') && ~isempty(cr.a_file)
        ctx.fileA = cr.a_file;
        ctx.saveATime = cr.a_save_time;
    else
        ctx.fileA = fullfile(ctx.dir, 'A.mat');
        t0 = tic;
        saveStruct(ctx.fileA, struct('A', A), ctx.ver);
        ctx.saveATime = toc(t0);
    end
    d = dir(ctx.fileA);  ctx.bytesA = d.bytes;
    ctx.bytesZ = 0;
    ctx.nWrite = 0;
end

function [ctx, tC] = doWrite(ctx, P)
    tmp = [ctx.fileZ(1:end-4) '_tmp.mat'];
    t0 = tic;
    saveStruct(tmp, P, ctx.ver);
    movefile(tmp, ctx.fileZ, 'f');
    tC = toc(t0);
    d = dir(ctx.fileZ);  ctx.bytesZ = d.bytes;
    ctx.nWrite = ctx.nWrite + 1;
end

function [A, tL] = doReadA(ctx)
    t0 = tic;
    S  = load(ctx.fileA, 'A');
    A  = S.A;
    tL = toc(t0);
end

function [P, tL] = doReadZ(ctx)
    t0 = tic;
    P  = load(ctx.fileZ);
    tL = toc(t0);
end

function doClose(ctx)
    if isempty(ctx) || ~isstruct(ctx) || ctx.keep, return; end
    if exist(ctx.dir, 'dir'), rmdir(ctx.dir, 's'); end
end

function saveStruct(file, S, ver)
    w = whos('S');
    if w.bytes > 2^31 - 1 && ~strcmp(ver, '-v7.3')
        ver = '-v7.3';
    end
    save(file, '-struct', 'S', ver);
end
