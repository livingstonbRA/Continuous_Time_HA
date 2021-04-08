classdef Model < handle
	properties
		% Parameters of the model.
		p;

		% Asset grids for the HJB.
		grids_HJB;

		% Asset grids for the KFE.
		grids_KFE;

		% The income process.
		income;

		% Instance of HJBSolver class.
		hjb_solver;

		% Instance of KFESolver class.
		kfe_solver;

		% Risk adjustment for SDU and risky asset
		risk_adj;

		% Instance of TransitionMatrixConstructor class.
		A_constructor_HJB;

		% Options
		options = struct();
	end

	methods
		function obj = Model(params, gridsHJB, gridsKFE, inc, varargin)
			obj.p = params;
			obj.grids_HJB = gridsHJB;
			obj.grids_KFE = gridsKFE;
			obj.income = inc;

			parser = inputParser;
		    addOptional(parser, 'quiet', false);
		    parse(parser, varargin{:});
		    obj.options.quiet = parser.Results.quiet;
		end

		function initialize(obj)
			import computation.TransitionMatrixConstructor
			import computation.HJBSolver
			import computation.HJBSolverSDU
			import computation.KFESolver

			if obj.p.SDU
				obj.hjb_solver = HJBSolverSDU(obj.p, obj.income, obj.p.hjb_options);
			else
				obj.hjb_solver = HJBSolver(obj.p, obj.income, obj.p.hjb_options);
			end

			returns_risk = obj.p.sigma_r > 0;
			obj.A_constructor_HJB = TransitionMatrixConstructor(obj.p, obj.income,...
				obj.grids_HJB, returns_risk);

			obj.kfe_solver = KFESolver(obj.p, obj.income,...
				obj.grids_KFE, obj.p.kfe_options, 'quiet', obj.options.quiet);
		end

		function [HJB, KFE, Au] = solve(obj, varargin)
			if obj.p.endogenous_labor
				% Find hours policy function at bmin over wage values.
				import model_objects.CRRA
				import model_objects.Frisch

				inc_mat = reshape(obj.income.y.vec, [1, 1, 1, obj.income.ny]);
				hours_bc_HJB = 0.5 * ones(obj.p.nb, 1, 1, obj.income.ny);
				hours_bc_KFE = 0.5 * ones(obj.p.nb_KFE, 1, 1, obj.income.ny);

				% hours_bc_last = hours_bc;
				for ih = 1:obj.p.HOURS_maxiters
					con_HJB = obj.income.nety_HJB_liq_hourly(hours_bc_HJB);
					con_KFE = obj.income.nety_KFE_liq_hourly(hours_bc_KFE);

					u1_HJB = CRRA.marginal_utility(con_HJB, obj.p.invies);
					u1_KFE = CRRA.marginal_utility(con_KFE, obj.p.invies);

					v1_HJB = (1-obj.p.directdeposit) * (1-obj.p.wagetax)...
						* hours_bc_HJB .* inc_mat .* u1_HJB;
					v1_KFE = (1-obj.p.directdeposit) * (1-obj.p.wagetax)...
						* hours_bc_KFE .* inc_mat .* u1_KFE;

					hours_bc_HJB = Frisch.inv_marginal_disutility(v1_HJB,...
						obj.p.labor_disutility, obj.p.frisch);
					hours_bc_KFE = Frisch.inv_marginal_disutility(v1_KFE,...
						obj.p.labor_disutility, obj.p.frisch);
					hours_bc_HJB = min(hours_bc_HJB, 1);
					hours_bc_KFE = min(hours_bc_KFE, 1);
					% max(abs(hours_bc(:)-hours_bc_last(:)))
					% hours_bc_last = hours_bc;
				end
			else
				hours_bc_HJB = 1;
				hours_bc_KFE = 1;
			end

			%% --------------------------------------------------------------------
			% HJB
			% ---------------------------------------------------------------------
			HJB = obj.solve_HJB(hours_bc_HJB);

			% Interpolate value function from HJB grids into KFE grids, if necessary
			if isequal(obj.grids_HJB.b.vec, obj.grids_KFE.b.vec) ...
				&& isequal(obj.grids_HJB.a.vec, obj.grids_KFE.a.vec)
				KFE = HJB;
			else
				V_KFE = zeros([obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny]);
				for iz = 1:obj.p.nz
					for iy = 1:obj.income.ny
						Vinterp = griddedInterpolant({obj.grids_HJB.b.vec, obj.grids_HJB.a.vec},...
							HJB.Vn(:,:,iz,iy), 'spline', 'nearest');
						V_KFE(:,:,iz,iy) = Vinterp({obj.grids_KFE.b.vec, obj.grids_KFE.a.vec});
					end
				end

				KFE = computation.find_policies(obj.p, obj.income,...
		    		obj.grids_KFE, V_KFE, hours_bc_KFE);
				KFE.Vn = V_KFE;
			end

			%% --------------------------------------------------------------------
			% KFE
			% ---------------------------------------------------------------------
			import computation.TransitionMatrixConstructor

			% True if returns should be treated as risky in the KFE
			returns_risk = (obj.p.sigma_r > 0) && (obj.p.retrisk_KFE == 1);
		    A_constructor_kfe = TransitionMatrixConstructor(obj.p,...
		    	obj.income, obj.grids_KFE, returns_risk);
		    Au = A_constructor_kfe.construct(KFE, KFE.Vn);   

			if obj.income.norisk
				wealth = NaN;
			else
				KFE.g = obj.kfe_solver.solve(Au);
				[mu_b, mu_a] = compute_wealth(KFE.g, obj.grids_KFE);
				if ~obj.options.quiet
					fprintf('    --- E[b] = %f ---\n', mu_b)
					fprintf('    --- E[a] = %f ---\n\n', mu_a)
				end
			end
		end
	end

	methods (Access=protected)
		function HJB = solve_HJB(obj, hours_bc)
			import computation.make_initial_guess

			[Vn, gguess] = make_initial_guess(obj.p, obj.grids_HJB,...
					obj.grids_KFE, obj.income);

			if (obj.income.norisk & ~obj.options.quiet)
				fprintf('    --- Iterating over HJB (no inc risk) ---\n')
			elseif ~obj.options.quiet
				fprintf('    --- Iterating over HJB ---\n')
			end

		    dst = 1e5;
			for nn	= 1:obj.p.HJB_maxiters
				[HJB, V_deriv_risky_asset_nodrift] = computation.find_policies(...
					obj.p, obj.income, obj.grids_HJB, Vn, hours_bc);

			    % Construct transition matrix 
		        [A, stationary] = obj.A_constructor_HJB.construct(HJB, Vn);

		        % Update value function
			    Vn1 = obj.hjb_solver.solve(obj.grids_HJB, A, HJB.u, Vn,...
			    	V_deriv_risky_asset_nodrift, stationary); % ignored unless SDU

				% check for convergence
			    Vdiff = Vn1 - Vn;
			    Vn = Vn1;
			    dst = max(abs(Vdiff(:)));
		        if ((nn==1) || (mod(nn,25)==0)) & ~obj.options.quiet
			    	fprintf('\tHJB iteration = %i, distance = %e\n', nn, dst);
		        end

		        if dst < obj.p.HJB_tol
		        	if ~obj.options.quiet
			        	fprintf('\tHJB converged after %i iterations\n', nn);
			        end
				 	break
				end
				check_if_not_converging(dst, nn);
		    end

		    if (dst > obj.p.HJB_tol)
		        error("HJB didn't converge");
		    end
		    
		    % Store value function and policies on both grids
		    HJB.Vn = Vn;
		end
	end
end

function check_if_not_converging(dst, nn)
	if dst>10 && nn>500
	 	% Not going to converge, throw exception
	 	msgID = 'HACT_Model:NotConverging';
	    msg = 'The HJB does not appear to be converging';
	    HJBException = MException(msgID,msg);
	    throw(HJBException)
	end
end

function [mu_b, mu_a] = compute_wealth(g, grdKFE)
	pmf = g .* grdKFE.trapezoidal.matrix;

	vals_b = pmf .* grdKFE.b.vec;
	mu_b = sum(vals_b(:));

	vals_a = pmf .* grdKFE.a.wide;
	mu_a = sum(vals_a(:));
end