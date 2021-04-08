classdef ParamsDefaults < handle
	properties (SetAccess=protected)
        % Run options
        ComputeMPCS = false;
        ComputeMPCS_illiquid = false;
        SimulateMPCS = false;
        ComputeMPCS_news = false;
        SimulateMPCS_news = false;
        Bequests = true;
        OneAsset = false;
        DealWithSpecialCase;
        SolveNoRisk = true;
        fast = false;
        saveGrids = false;
        makePlots = false;
        IncomeDescr = '';
        
        % Name/nums for the parameterization
        name = 'unnamed';
        group_num = NaN;
        label;

        % Numeric identifier for the parameterization
        param_index;

        % Path to the income variables
        income_dir;

        % Turn off transitory income risk
        no_transitory_incrisk = false;

        %% -------------------------------------------
    	% Liquid Asset Grid Parameters
    	% --------------------------------------------
        % Min value of liquid assets
    	bmin = 0;

    	% Minimum value of liquid assets before
    	% borrowing wedge kicks in
    	b_soft_constraint = 0;

    	% Max value for liquid assets
    	bmax = 20;

        bgrid_term1_weight = 0;

        bgrid_term1_curv = 0.8;
        
    	% Curvature of positive section of liquid asset,
    	% lower value implies more curvature
    	b_gcurv_pos = 0.35;

    	% Curvature of negative section of liquid asset,
    	% lower value implies more curvature
    	b_gcurv_neg = 0.4;

    	% Number of points on liquid asset grid
        nb = 45;
        
        % Number of points on positive section of liquid
        % asset grid, defaults to nb
    	nb_pos;

    	% Number of points on negative section of liquid
    	% asset grid, defaults to 0
        nb_neg;

            % Number of points on liquid asset grid for the KFE
    	nb_KFE = 45;

    	% Number of points on positive section of liquid
    	% asset grid for the KFE, defaults to nb_KFE
    	nb_pos_KFE;

    	% Number of points on negative section of liquid
    	% asset grid for the KFE, defaults to 0
        nb_neg_KFE;

            %% -------------------------------------------
    	% Illiquid Asset Grid Parameters
    	% --------------------------------------------
    	% Min value of illiquid assets
    	amin = 0;

    	% Max value of illiquid assets
        amax = 500;

        % Number of points on the illiquid asset grid
        na = 45;

        agrid_term1_weight = 0;

        agrid_term1_curv = 0.6;
        
        % Curvature of positive section of illiquid asset,
		% lower value implies more curvature
        a_gcurv = 0.3;

        % Number of points on illiquid asset grid for
        % the KFE
        na_KFE = 45;

        %% -------------------------------------------
    	% Other Parameters
    	% --------------------------------------------
        % Number of income states
        ny;

        % Number of states in z-dimension
        nz;
        
	    %% -------------------------------------------
    	% Returns
    	% --------------------------------------------
    	% Liquid returns, at a quarterly rate
    	r_b = 0.02 / 4;

    	% Illiquid returns, at a quarterly rate
        r_a = 0.06 / 4;

        % Annuities, 0 turns off annuities
		perfectannuities = 0;

        % Wedge between borrowing rate and r_b on the
        % liquid asset, only relevant with bmin < 0
        borrwedge = 0.08/4;

        % Liquid asset borrowing rate, set according
        % to borrwedge
        r_b_borr;

        % Standard deviation of returns risk
        sigma_r = 0;

        % Boolean indicator, set to true to include
        % returns risk in the KFE
        retrisk_KFE = false;

    	%% -------------------------------------------
    	% Preferences
    	% --------------------------------------------

    	% Stochastic differential utility
        SDU = false;

        % Coefficient of risk aversion, can be set
        % to a row vector to accomodate RA heterogeneity
		riskaver = 1;

        % Inverse of the IES, only set this for SDU
        % preferences. Otherwise it is automaticaly set
        % equal to 'riskaver'.
        invies = 1;

        % Array of risk aversion coefficients, set
        % automatically
    	riskaver_fulldim;

    	% Transition into death, at a quarterly rate
    	deathrate = 1 / 200;

    	% Time discount factor, at a quarterly rate
    	rho = 0.015;

    	% Grid for rho values if there is rho
    	% heterogeneity, automatically set to
    	% rho + rho_grid
    	rhos;

    	% Grid to accomodate rho heterogeneity
    	rho_grid = [0];

        %% -------------------------------------------
        % CALIBRATION
        % --------------------------------------------

    	% ------------ targets ----------------------------
        calibrate;
        calibrator;
    	calibration_vars;
        calibration_stats;
        calibration_targets;
        calibration_bounds;
        calibration_backup_x0;
        calibration_maxiters = 200;
        calibration_crit = 1e-5;
        calibration_scales = [];

        % ------------ taxes ------------------------------
        transfer = 0; % transfer to households 
        wagetax = 0; % tax rate on wage income

    	% ------------ approximation parameters -----------
    	% HJB loop
    	hjb_options;
    	HJB_maxiters = 3000; % maximal allowable number of HJB iterations
    	HJB_tol = 1e-8; % critical value
    	HJB_delta = 1e6; % step size
        HJB_implicit = false;

    	% Howard improvement step in HJB loop
    	HIS_maxiters = 10; % total number of Howard improvement steps
    	HIS_start = 2; % when in HJB loop do Howard improvement steps begin?
    	HIS_tol= 1e-5; % critical value

    	% KFE loop
    	kfe_options;
    	KFE_maxiters = 1e5; % maximal allowable number of KFE iterations
    	KFE_tol = 1e-8; % critical value
    	KFE_delta = 1e6; %1e6; % step size
        KFE_iterative = true;

        % Options for MPC computations
        MPC_delta = 0.005;
        MPC_interp_method = "linear";
        mpc_options;
        mpc_options_illiquid;

        % Options for MPCs out of news
        MPCS_News_delta_terminal = 1e-3;
        MPCS_News_compute_mpcs;
        mpcs_news_options;

        % Options for MPC simulations
        MPCSim_T = 200;
        MPCSim_n = 5e4;
        MPCSim_interp_method = "makima";
        mpcsim_options;

        % Endogenous labor parameters
        HOURS_maxiters = 1e3;
        endogenous_labor = false;
        labor_disutility = 0;
        frisch = 0.5;
        MPL = 1;
        MPK = 1;
        
        % ----------- statistics variables -------------------------------
        epsilon_HtM = [0 0.005 0.01 0.02 0.05 0.1 0.15]; % for looking at fraction HtM
        dollars_HtM = [500, 1000, 2000, 5000, 10000];
        wpercentiles = [10 25 50 90 95 98 99 99.9];
        mpc_shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        numeraire_in_dollars;
        decomp_thresholds = [0 0.01 0.05];

        % ----------- adjustment costs -----------------------------------
        kappa0 = 0;
        kappa1 = 5;
        kappa2 = 1;

        a_lb = 0.25;

        % deposit share     
        directdeposit = 0;
    end
end