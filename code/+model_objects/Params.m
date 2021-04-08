classdef Params < model_objects.ParamsDefaults
    % This class stores the parameters of the model and contains methods to
    % adjust them for frequency and other factors.
    %
    % See the base class, ParamsDefaults, for default values.

    methods
        function obj = Params(varargin)
            obj.set_from_structure(varargin{:});
            
            obj.kfe_options = struct(...
				'maxiters', obj.KFE_maxiters,...
				'tol', obj.KFE_tol,...
				'delta', obj.KFE_delta,...
				'iterative', obj.KFE_iterative);

            obj.hjb_options = struct(...
				'delta', obj.HJB_delta,...
				'implicit', obj.HJB_implicit,...
				'HIS_maxiters', obj.HIS_maxiters,...
				'HIS_tol', obj.HIS_tol,...
				'HIS_start', obj.HIS_start);

            obj.mpc_options = struct(...
            	'delta', obj.MPC_delta,...
            	'interp_method', obj.MPC_interp_method,...
                'liquid_mpc', true);

            obj.mpc_options_illiquid = obj.mpc_options;
            obj.mpc_options_illiquid.liquid_mpc = false;

            obj.mpcs_news_options = struct(...
				'delta', obj.MPC_delta,...
				'delta_terminal', obj.MPCS_News_delta_terminal,...
				'save_policies', obj.SimulateMPCS_news,...
				'compute_mpcs', obj.ComputeMPCS_news);

            obj.mpcsim_options = struct(...
            	'T', obj.MPCSim_T,...
            	'n', obj.MPCSim_n,...
            	'interp_method', obj.MPCSim_interp_method);

            % Check for other heterogeneity
            if (numel(obj.rho_grid)>1) && (numel(obj.riskaver)>1)
            	error('Cannot have both rho and riskaver heterogeneity')
            else
                obj.nz = max(numel(obj.rho_grid), numel(obj.riskaver));
            end

            obj.rhos = obj.rho + obj.rho_grid;
            obj.riskaver_fulldim = reshape(obj.riskaver, [1 1 numel(obj.riskaver) 1]);

            if ~obj.SDU
                obj.invies = obj.riskaver;
            end

            % adjust interest rates
            obj.r_b_borr = obj.r_b + obj.borrwedge;

            % adjust if oneasset is selected
            if obj.OneAsset
                obj.na = 2;
                obj.na_KFE = 2;
                obj.kappa0 = 1e8;
                obj.r_a = obj.r_b;
                obj.ComputeMPCS_illiquid = false;
            end

            if obj.fast
                obj.nb = 16;
                obj.nb_pos = 16;
                obj.nb_neg = 0;
                obj.nb_KFE = 15;
                obj.nb_pos_KFE = 15;
                obj.nb_neg_KFE = 0;

                % illiquid grid parameters
                if obj.OneAsset == 0
                    obj.na = 14;
                    obj.na_KFE = 13;
                end
                
                obj.MPCSim_T = 100;
                obj.MPCSim_n = 1e3;
            end

            % Set default grid sizes
            if isempty(obj.nb_pos)
                obj.nb_pos = obj.nb;
            end
            if isempty(obj.nb_pos_KFE)
                obj.nb_pos_KFE = obj.nb_KFE;
            end
            
            obj.nb_neg = obj.nb - obj.nb_pos;
            obj.nb_neg_KFE = obj.nb_KFE - obj.nb_pos_KFE;

            if obj.calibrate
                obj.set_calibrator();
            end
        end

        function set_from_structure(obj, varargin)
    		parser = inputParser;

    		import model_objects.ParamsDefaults

    		defaults = ParamsDefaults();
    		fields = properties(ParamsDefaults);
            struct_options = varargin{1};

    		for k = 1:numel(fields)
    			field = fields{k};
                if ~ismember(field, {'hjb_options', 'kfe_options'})
                    addParameter(parser, field, defaults.(field));
                end
                
                if isfield(struct_options, field)
                    if ~isempty(struct_options.(field))
                        cleaned.(field) = struct_options.(field);
                    end
                end
            end
            
    		parse(parser, cleaned);

    		for k = 1:numel(fields)
    			field = fields{k};
                if ~ismember(field, {'hjb_options', 'kfe_options'})
                    obj.(field) = parser.Results.(field);
                end
    		end
    	end
        
        function obj = set(obj, field, new_val, quiet)
        	% Sets the value of a parameter.
        	%
        	% Inputs
        	% ------
        	%
        	% field : A string containing the parameter
        	%	name.
        	%
        	% new_val : The desired value of the parameter.
        	%
        	% quiet : An optional argument that, when it
        	%	evaluates to true, suppresses printing
        	%	to the screen.

        	field = char(field);

        	KFE_option_passed = false;
        	if numel(field) > 4
        		if strcmp(field(1:3), 'KFE')
        			KFE_option_passed = true;
        		end
        	end

        	if strcmp(field, 'rho')
        		obj.rhos = new_val + obj.rho_grid;
        	elseif KFE_option_passed
        		assert(isprop(obj.kfe_options, field(5:end)), "Invalid KFE option");
        		obj.kfe_options.set(field(5:end), new_val);
    		elseif ~isprop(obj, field)
    			error("Requested field is not an attribute of Params.");
        	end

            obj.(field) = new_val;

            if ~exist('quiet', 'var')
            	quiet = false;
            end

            if ~quiet
            	disp(strcat(field, sprintf(" has been reset to %.9f", new_val)));
            end
        end

        function set_calibrator(obj)
            import model_objects.HACTCalibrator

            calibrator = HACTCalibrator(obj, obj.calibration_vars,...
                obj.calibration_stats, obj.calibration_targets);

            if ~isempty(obj.calibration_scales)
                calibrator.set_fscale(obj.calibration_scales);
            end
            
            if ~isempty(obj.calibration_backup_x0)
                calibrator.add_backup_x0(obj.calibration_backup_x0{:});
            end

            calibrator.set_param_bounds(obj.calibration_bounds{:});
            calibrator.set_handle(obj);
            obj.calibrator = calibrator;
        end

        function print(obj)
            fprintf('\n\nSelected parameterization %d:\n',obj.param_index) 
            fprintf('%s\n\n',obj.name)

            fprintf('Chosen parameters were...\n\n')

            fprintf('\tb_soft_constraint = %f\n',obj.b_soft_constraint)
            fprintf('\tbmin = %f\n',obj.bmin)
            fprintf('\tbmax = %f\n',obj.bmax)
            fprintf('\tb_gcurv_pos = %f\n',obj.b_gcurv_pos)
            fprintf('\tb_gcurv_neg = %f\n',obj.b_gcurv_neg)
            fprintf('\tnb = %i\n',obj.nb)
            fprintf('\tnb_pos = %i\n',obj.nb_pos)
            fprintf('\tnb_KFE = %i\n',obj.nb_KFE)
            fprintf('\tnb_pos_KFE = %i\n',obj.nb_pos_KFE)
            fprintf('\n')

            if obj.OneAsset == 0
                fprintf('\tna = %i\n',obj.na)
                fprintf('\tna_KFE = %i\n',obj.na_KFE)
                fprintf('\tamin = %f\n',obj.amin)
                fprintf('\tamax = %f\n',obj.amax)
                fprintf('\ta_gcurv = %f\n',obj.a_gcurv)
                fprintf('\n')

                fprintf('\ta_lb = %f\n',obj.a_lb)
                fprintf('\n')
            end

            fprintf('\tBequests = %i\n',obj.Bequests)
            fprintf('\n')

            fprintf('\tr_b = %f\n',obj.r_b)
            fprintf('\tr_b_borr = %f\n',obj.r_b_borr)
            fprintf('\tr_a = %f\n',obj.r_a)
            fprintf('\tperfectannuities = %f\n',obj.perfectannuities)
            fprintf('\triskaver = %f\n',obj.riskaver)
            fprintf('\tdeathrate = %f\n',obj.deathrate)
            fprintf('\n\n')
        end

        function label_out = quantity2label(obj, val)
            if ~isempty(obj.numeraire_in_dollars)
                dollars = abs(val) * obj.numeraire_in_dollars;
                if val < 0
                    pref = '-$';
                else
                    pref = '$';
                end

                label_out = sprintf('%s%g', pref, dollars);
            else
                pref = '';
                label_out = sprintf('%g', val);
            end
        end
    end
end