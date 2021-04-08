classdef HJBSolver < handle

	properties (Constant)
		% Default options
		defaults = struct(...
				'implicit', false,...
				'delta', 1e5,...
				'HIS_maxiters', 0,...
				'HIS_tol', 1e-5,...
				'HIS_start', 2 ...
				);
 
	end

	properties (SetAccess=protected)
		p;

		income;

		% Total number of states.
		n_states;

		% Number of states per income level.
		states_per_income;

		% Number of states along each dim, a vector.
		shape;

		% Sparse matrix of discount factors.
		rho_mat;
        
        % An HJBOptions object.
        options;

        % Current iteration, needed for the HIS.
        current_iteration = 0;

        % Used only for SDU
        sdu_k_adj = 0;
        risk_adj = 0;
	end

	methods
		function obj = HJBSolver(p, income, varargin)
			% See the properties block above or view this class'
			% documentation to see the requirements of the
			% input parameters.

			obj.p = p;
			obj.income = income;

			obj.options = obj.parse_options(varargin{:});

			obj.n_states = p.nb * p.na * p.nz * income.ny;
			obj.states_per_income = p.nb * p.na * p.nz;
			obj.shape = [p.nb p.na p.nz income.ny];

			% ---------------------------------------------------------
			% Discount Factor Matrix
			% ---------------------------------------------------------
			obj.create_rho_matrix();
		end

		function V_update = solve(obj, grd, A, u, V, varargin)
			% Updates the value function (grd only used for SDU derived class)
			if obj.options.implicit
				V_update = obj.solve_implicit(A, u, V);
			else
				V_update = obj.solve_implicit_explicit(A, u, V);
			end
		end
	end

	methods (Access=protected)
		%%--------------------------------------------------------------
	    % Implicit Updating
	    % --------------------------------------------------------------
		function Vn1 = solve_implicit(obj, A, u, V)
			assert(obj.p.deathrate == 0,...
				"Fully implicit updating does not support death.")

			A_income = obj.income.full_income_transition_matrix(obj.p, V);

			% Note: risk_adj = 0 unless using SDU w/risky asset
		    RHS = obj.options.delta * (u(:) + obj.risk_adj(:)) + Vn(:);
	        
	        B = (obj.rho_mat - A - A_income) * obj.options.delta + speye(obj.n_states);
	        Vn1 = B \ RHS;
	        Vn1 = reshape(Vn1, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		%%--------------------------------------------------------------
	    % Implicit-Explicit Updating
	    % --------------------------------------------------------------
		function [Vn1, Bk_inv] = solve_implicit_explicit(obj, A, u, V)
			import computation.hjb_divisor
			obj.current_iteration = obj.current_iteration + 1;

			% Update value function
			u_k = reshape(u, [], obj.income.ny);
			Vn_k = reshape(V, [], obj.income.ny);
			Vn1 = zeros(obj.states_per_income, obj.income.ny);
			for k = 1:obj.income.ny
				[inctrans, inctrans_k] = obj.get_implicit_explicit_inctrans(k);
				
				% Solve HJB
				Bk = hjb_divisor(obj.options.delta, obj.p.deathrate, k,...
					A, inctrans, obj.rho_mat);
	        	Vn1(:,k) = obj.update_Vk_implicit_explicit(...
	        		Vn_k, u_k, k, Bk, inctrans_k);
	        end

	        Vn1 = reshape(Vn1, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		function Vn1_k = update_Vk_implicit_explicit(...
			obj, V_k, u_k, k, Bk, inctrans_k)
        	indx_k = ~ismember(1:obj.income.ny, k);

            offdiag_inc_term = sum(...
            	squeeze(inctrans_k(:,indx_k)) .* V_k(:,indx_k), 2);

            RHSk = obj.options.delta * (u_k(:,k) + offdiag_inc_term) + V_k(:,k);
            
            % Note: sdu_k_adj = 0 unless using SDU w/risky asset
            Vn1_k = Bk \ (RHSk + obj.sdu_k_adj);
        end

		function options = parse_options(obj, varargin)
			import computation.HJBSolver
			import aux.parse_keyvalue_pairs

			defaults = HJBSolver.defaults;
			options = parse_keyvalue_pairs(defaults, varargin{:});
		end

		function obj = create_rho_matrix(obj)
			import aux.sparse_diags
			if obj.options.implicit
				% discount factor values
		        if numel(obj.p.rhos) > 1
		            rhocol = repmat(kron(obj.p.rhos(:), ones(obj.p.nb*obj.p.na, 1)), obj.income.ny, 1);
		            obj.rho_mat = sparse_diags(rhocol, 0);
		        else
		            obj.rho_mat = obj.p.rho * speye(obj.n_states);
		        end
		    else
		    	if numel(obj.p.rhos) > 1
			        rhocol = kron(obj.p.rhos(:), ones(obj.p.nb*obj.p.na,1));
			        obj.rho_mat = sparse_diags(rhocol, 0);
			    else
			        obj.rho_mat = obj.p.rho * speye(obj.states_per_income);
			    end
	    	end
		end

		function Vn2_k = howard_improvement_step(obj, Vn1_k, u_k, Bik_all)
		    % Technique to speed up convergence.

		    for jj = 1:obj.options.HIS_maxiters
		        Vn2_k = zeros(obj.states_per_income, obj.income.ny);
		        for kk = 1:obj.income.ny
		            indx_k = ~ismember(1:obj.income.ny, kk);
		            
		            Vkp_stacked = sum(...
		                	repmat(obj.income.ytrans(kk, indx_k), obj.states_per_income,1)...
		                	.* Vn1_k(:,indx_k), 2);
		            qk = obj.options.delta * (u_k(:,kk) + Vkp_stacked) + Vn1_k(:,kk);
		            Vn2_k(:,kk) = Bik_all{kk} * qk;
		        end

		        dst = max(abs(Vn2_k(:) - Vn1_k(:)));
		        Vn1_k = Vn2_k;
		        if dst < obj.options.HIS_tol
		            break
		        end
    		end
		end

		function [inctrans, inctrans_k] = get_implicit_explicit_inctrans(obj, k)
			inctrans = obj.income.ytrans(k,k) * speye(obj.states_per_income);
			inctrans_k = repmat(obj.income.ytrans(k,:), obj.states_per_income, 1);
		end
	end
end