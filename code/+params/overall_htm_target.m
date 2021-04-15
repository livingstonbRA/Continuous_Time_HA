function [outparams, n] = overall_htm_target(param_opts)
	import params.scf2019struct
    
    scf = scf2019struct();

	shocks = [-1, -500, -5000, 1, 500, 5000];

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / (scf.quarterly_earnings * 4);
    shared_params.numeraire_in_dollars = (scf.quarterly_earnings * 4);
    shared_params.no_transitory_incrisk = false;
    shared_params.Bequests = false;
    shared_params.r_b = 0.01 / 4;

    anninc = shared_params.numeraire_in_dollars;
    shared_params.a_lb = 500 / anninc;
    
    incomedirs = {'continuous_a/no_measurement_error',...
        'continuous_a/measurement_error_20pc',...
        'continuous_a/measurement_error_33pc',...
        'continuous_a/measurement_error_50pc'};

    IncomeDescriptions = {'cont_a, no meas err',...
        'cont_a, meas err 20pc',...
        'cont_a, meas err 33pc',...
        'cont_a, meas err 50pc'};

    experiment = false;
    if experiment
        iy = 1;
        params = shared_params;
        params.calibration_vars = {'rho', 'r_b'};
        params.calibration_stats = {'median_totw', 'median_liqw'};
        params.calibration_targets = [1.54, 0.05];
        params.calibration_scales = [1, 100];
        params.income_dir = incomedirs{iy};
        params.IncomeDescr = IncomeDescriptions{iy};
        param_opts.param_index = 1;

        params.kappa1 = 2;
        params.kappa2 = 0.5;
        params.rho = 0.012;
        params.r_b = 0.005;
        params.r_a = 0.013;

        params.sd_r = 0.01;
        params.SDU = true;
        params.invies = 1 / 1.5;
        params.riskaver = 2;
        params.kappa1 = 0.1;
        params.r_a = 0.02;

        rho_bds = [0.003, 0.02];
        r_b_bds = [0.004, 0.0085];
        params.KFE_maxiters = 1e6;

        % Set calibrator
        params.calibration_bounds = {rho_bds, r_b_bds};
        params.calibration_backup_x0 = {};

        params = {params};
    else
        params = {};

        %% TARGET MEDIAN TOTAL WEALTH AND MEDIAN LIQUID WEALTH
        % Iterate over r_a, rho
        median_calibration = shared_params;
        median_calibration.calibration_vars = {'rho', 'r_a'};

        kappa_1s = [0.2:0.2:1, 1.5:0.5:5];
        kappa_2s = [0.5, 1.0, 1.5];

        calibrations = {median_calibration};

        ii = 1;
        group_num = 0;
        for icalibration = [1]
            for kappa2 = kappa_2s
                for iy = 1:3
                    group_num = group_num + 1;
                    for kappa1 = kappa_1s
                        params = [params {calibrations{icalibration}}];
                        params{ii}.name = sprintf('iy=%d, kappa2=%g', iy, kappa2);
                        params{ii}.kappa2 = kappa2;
                        params{ii}.income_dir = incomedirs{iy};
                        params{ii}.IncomeDescr = IncomeDescriptions{iy};
                        params{ii}.group_num = group_num;
                        
                        if params{ii}.no_transitory_incrisk
                            % params{ii}.rho = 0.001;
                            % params{ii}.r_a = 0.0052;
                            % params{ii}.calibration_bounds = {[0.0008, 0.003],...
                            % [shared_params.r_b + 0.0003, 0.009]};
                            % params{ii}.calibration_backup_x0 = {};
                        else
                            params{ii}.kappa1 = kappa1;
                            params{ii}.rho = 0.005;
                            params{ii}.r_a = 0.007;

                            rho_bds = [0.0005, 0.03];
                            % r_a_bds = [0.008, 0.02];
                            r_a_bds = [0.003, 0.04];
                            params{ii}.KFE_maxiters = 1e6;
                            % params{ii}.a_lb = 0.3;

                            % params{ii}.rho = mean(rho_bds);
                            % params{ii}.r_a = mean(r_a_bds);

                            % Set calibrator
                            params{ii}.calibration_bounds = {rho_bds, r_a_bds};
                            params{ii}.calibration_backup_x0 = {};
                        end
                        % params{ii}.calibration_stats = {'diff_median', 'median_liqw'};
                        % params{ii}.calibration_targets = [1.49, 0.05];
                        % params{ii}.calibration_scales = [1, 10];

                        % params{ii}.calibration_stats = {'diff_mean', 'liqw'};
                        % params{ii}.calibration_targets = [4.1-0.56, 0.56];
                        params{ii}.calibration_stats = {'totw', 'median_liqw'};
                        params{ii}.calibration_targets = [scf.mean_totw, scf.median_liqw];
                        params{ii}.calibration_scales = [1, 10];

                        ii = ii + 1;
                    end
                end
            end
        end
    end
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end