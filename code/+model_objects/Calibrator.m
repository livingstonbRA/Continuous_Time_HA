classdef Calibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		options;
		variables;
		target_names;
		target_values;
		target_result;
		lbounds = [];
		ubounds = [];
		fscale = [];
		n;
		x0;
		iter = 1;
		dnorm = 1e5;

		ix_curr = 1;

		main_handle;

		results_handle;

		solver_handle;
	end

	methods
		function obj = Calibrator(params, variables,...
			target_names, target_values)
			obj.construct_options_struct(params);

			obj.variables = variables;
			obj.n = numel(variables);

			obj.target_names = target_names;
			obj.target_values = target_values;

			if obj.n ~= numel(obj.target_names)
				error("Number of instruments and targets don't match")
			elseif numel(obj.target_values) ~= numel(obj.target_names)
				error("Too many/few values provided for targets")
			end

			obj.fscale = ones(obj.n, 1);

			x0_1 = zeros(1, obj.n);
			for i_var = 1:obj.n
				x0_1(i_var) = params.(obj.variables{i_var});
			end
			
			obj.x0 = {x0_1};
		end

		function construct_options_struct(obj, params)
			obj.options = struct();
		end

		function add_backup_x0(obj, varargin)
			for ii = 1:numel(varargin)
				obj.x0 = [obj.x0, varargin{ii}];
			end
		end

		function set_fscale(obj, fscale_in)
			assert(numel(fscale_in) == obj.n, 'Invalid number of scaling factors');
			obj.fscale = fscale_in(:);
		end

		function set_param_bounds(obj, varargin)
			nv = numel(varargin);
			assert(nv==obj.n, "Too many or too few bounds provided");

			for ii = 1:nv
				var_bounds = varargin{ii};
				obj.lbounds(ii) = var_bounds(1);
				obj.ubounds(ii) = var_bounds(2);
			end

			for ix0 = 1:numel(obj.x0)
				x0 = obj.x0{ix0};
				for ii = 1:obj.n
					if x0(ii) < obj.lbounds(ii)
						x0(ii) = (5*obj.lbounds(ii) + obj.ubounds(ii))/6;
					elseif x0(ii) > obj.ubounds(ii)
						x0(ii) = (obj.lbounds(ii) + 5*obj.ubounds(ii))/6;
                    end
                    obj.x0{ix0} = x0;
				end
			end
		end

		function set_handle(obj, p)
			if isempty(obj.main_handle)
				error('Derived class has not initialized main_handle')
			end

			obj.solver_handle = @(x) obj.fn_handle(x, p);
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);
			quiet = true;

			for i_var = 1:obj.n
				current_params.set(obj.variables{i_var}, x(i_var), quiet);
			end
			results = obj.main_handle(current_params);

			fprintf('  -- function evaluation %d --\n    evaluated at: ', obj.iter)
			v = zeros(obj.n, 1);
			for i_var = 1:obj.n
				v(i_var) = obj.get_results_value(results, obj.target_names{i_var});
				fprintf('%s = %g', obj.variables{i_var}, x(i_var))
				if i_var < obj.n
					fprintf(", ")
				end
			end
			dv = v(:) - obj.target_values(:);
			obj.target_result = v(:);
			obj.dnorm = norm(dv);
			dv = dv .* obj.fscale;

			fprintf('\n    target variables: ')
			for i_var = 1:obj.n
				fprintf('%s = %g', obj.target_names{i_var}, v(i_var))
				if i_var < obj.n
					fprintf(", ")
				end
			end
			fprintf('\n    norm: %f\n', norm(dv))

			dv = obj.adjust_dv(results, current_params, dv);

			if obj.n > 1
				% Guide solver to reduce r_a if illiq asset very high
				dv = [dv; max(v(1) - 200, 0) * current_params.r_a];

				% Guide solver to increase r_a if illiq asset very low
				dv = [dv; 10 * min(v(1) - 1, 0) * (0.1 - current_params.r_a)];
			end

			obj.reset_param_options(current_params);
			obj.iter = obj.iter + 1;
		end

		function value = get_results_value(obj, results, variable_name)
			value = nan;
		end

		function dv = adjust_dv(obj, results, current_params, dv)
			% Do nothing unless overriden
			dv = dv;
		end

		function solver_args = get_args(obj)
			x0 = obj.get_next_x0();
			if ~isempty(obj.ubounds) && ~isempty(x0)
				solver_args = {x0, obj.lbounds, obj.ubounds};
			elseif ~isempty(x0)
				solver_args = {x0};
			else
				solver_args = {};
			end
		end

		function turn_off_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, 0, quiet);
			end
		end

		function reset_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, obj.options.(props{ip}), quiet);
			end
		end

		function x0curr = get_next_x0(obj)
            if (obj.ix_curr <= numel(obj.x0))
                x0curr = obj.x0{obj.ix_curr};
                obj.ix_curr = obj.ix_curr + 1;
            else
                x0curr = [];
            end
        end
	end
end
