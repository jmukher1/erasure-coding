function cr = cr_parse_options(varargin)
% name/value options of the checkpoint/restart baselines 

    p = inputParser;
    p.FunctionName = 'cr_parse_options';
    addParameter(p, 'interval_mode',   "static");
    addParameter(p, 'static_interval', 5);
    addParameter(p, 'ckpt_dir',        []);
    addParameter(p, 'mat_version',     '-v7');
    addParameter(p, 'charge_save_A',   false);
    addParameter(p, 'keep_files',      false);
    addParameter(p, 'cr_verbose',      1);
    addParameter(p, 'a_file',          []);
    addParameter(p, 'a_save_time',     0);
    addParameter(p, 'charge_load_A',   []);
    parse(p, varargin{:});
    cr = p.Results;

    cr.interval_mode = lower(string(cr.interval_mode));
    valid = ["none", "restart", "static"];
    if ~ismember(cr.interval_mode, valid)
        error('cr_parse_options:badMode', 'interval_mode must be one of: %s.', strjoin(valid, ', '));
    end
    if ~ismember(cr.mat_version, {'-v6', '-v7', '-v7.3'})
        error('cr_parse_options:badVersion', 'mat_version must be ''-v6'', ''-v7'' or ''-v7.3''.');
    end
    if cr.interval_mode == "static" && ~(isscalar(cr.static_interval) && cr.static_interval >= 1)
        error('cr_parse_options:badInterval', 'static_interval must be a positive integer.');
    end
    cr.static_interval = round(cr.static_interval);
    if isempty(cr.charge_load_A)
        cr.charge_load_A = (cr.interval_mode ~= "restart");
    end
    cr.charge_load_A = logical(cr.charge_load_A);
    if ~isempty(cr.a_file)
        cr.a_file = char(cr.a_file);
        if ~exist(cr.a_file, 'file')
            error('cr_parse_options:noAFile', 'a_file %s does not exist.', cr.a_file);
        end
    end
end
