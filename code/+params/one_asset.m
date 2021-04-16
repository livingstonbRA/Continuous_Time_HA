function [outparams, n] = one_asset(param_opts)
	import params.scf2019struct
    
    scf = scf2019struct();

	shocks = [-1, -500, -5000, 1, 500, 5000];

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / (scf.quarterly_earnings * 4);
    shared_params.numeraire_in_dollars = (scf.quarterly_earnings * 4);
    shared_params.no_transitory_incrisk = false;
    shared_params.nb = 150;
    shared_params.nb_KFE = 150;
    shared_params.bgrid_term1_weight = 0;
    shared_params.bgrid_term1_curv = 1;
    shared_params.b_gcurv_pos = 0.2;
    shared_params.OneAsset = true;
    shared_params.Bequests = false;
    shared_params.r_b = 0.01 / 4;

    shared_params.bmax = 500;
    shared_params.rho = 0;

    calibration = shared_params;
    calibration.calibration_vars = {'rho'};
    calibration.calibration_bounds = {[-0.01, 0.02]};
    calibration.calibration_stats = {'totw'};
    calibration.calibration_targets = [scf.mean_totw];
    calibration.calibration_scales = [1];
    
    incomedirs = {...
        'continuous_a/no_measurement_error',...
        'continuous_a/measurement_error_20pc',...
        'continuous_a/measurement_error_33pc',...
        'continuous_a/measurement_error_50pc'};

    IncomeDescriptions = {...
        'cont_a, no meas err',...
        'cont_a, meas err 20pc',...
        'cont_a, meas err 33pc',...
        'cont_a, meas err 50pc'};

    params = {};

    for ii = 1:4
        params{ii} = calibration;
        params{ii}.name = IncomeDescriptions{ii}; 
        params{ii}.income_dir = incomedirs{ii};
        params{ii}.IncomeDescr = IncomeDescriptions{ii};
    end

    %% --------------------------------------------------------------------
    % HOUSEKEEPING, DO NOT CHANGE
    % ---------------------------------------------------------------------
    n = numel(params);
    outparams = params{param_opts.param_index};
end