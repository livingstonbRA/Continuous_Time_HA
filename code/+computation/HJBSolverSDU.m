classdef HJBSolverSDU < computation.HJBSolver

	properties
		sdu_adj;
	end


	methods
		function obj = HJBSolverSDU(p, income, varargin)
			obj = obj@computation.HJBSolver(p, income, varargin{:});
		end

		function V_update = solve(obj, grd, A, u, V, varargin)
			% Updates the value function.
			obj.check_inputs(A, u, V);

			obj.risk_adj = obj.compute_risk_adjustment_for_nodrift_case(...
	    		grd, V, varargin{:});

			if obj.options.implicit
				V_update = obj.solve_implicit(A, u, V);
			else
				obj.sdu_adj = obj.income.income_transitions_SDU(obj.p, V);

				V_update = obj.solve_implicit_explicit(A, u, V);
			end
		end
	end

	methods (Access=protected)
		function risk_adj = compute_risk_adjustment_for_nodrift_case(...
			obj, grd, Vn, V_deriv_risky_asset_nodrift, stationary)
			% computes the adjustment term when returns are risky
			% and there is no drift in the risky assets for some
			% asset > 0 cases
			
			% Parameters
			% ---------
			% p : a Params object
			%
			% grd : a Grid object
			%
			% V_deriv_risky_asset_nodrift : approximation of the
		    %	first derivative of V for the case with no drift
		    %	in the risky asset
		    %
		    % stationary : boolean mask to indicate states where
		    %	there is no drift in risky asset
		    %
		    % Vn : the value function over states
		    %
		    % Returns
		    % -------
		    % risk_adj : an array of shape (nb, na, nz, ny) with
		    % 	zeros everywhere except for states where risky
		    %	asset drift is zero, where the array contains
		    %	the term with Va^2 

			if ~isempty(stationary) & (obj.p.sigma_r > 0)
				% there are states with neither backward nor forward drift,
				% need to compute additional term for risk
				if obj.p.invies == 1
					risk_adj = (1-p.riskaver) * V_deriv_risky_asset_nodrift .^ 2;
				else
					risk_adj = V_deriv_risky_asset_nodrift .^ 2 ./ Vn ...
						* (obj.p.invies - obj.p.riskaver) / (1-obj.p.invies);
				end

				risk_adj = risk_adj .* (grd.a.wide * obj.p.sigma_r) .^ 2 / 2;
				risk_adj(~stationary) = 0;
				risk_adj = reshape(risk_adj, [], obj.p.ny);
			else
				risk_adj = 0;
			end
		end

		function [inctrans, inctrans_k] = get_implicit_explicit_inctrans(obj, k)
			inctrans = aux.sparse_diags(obj.sdu_adj(:,k,k), 0);
			inctrans_k = squeeze(obj.sdu_adj(:,k,:));
		end

		function Vn1_k = update_Vk_implicit_explicit(...
			obj, V_k, u_k, k, Bk, inctrans_k)   

			if obj.p.sigma_r > 0
	           	obj.sdu_k_adj = obj.options.delta * obj.risk_adj(:,k);
	        end

			Vn1_k = update_Vk_implicit_explicit@computation.HJBSolver(...
				obj, V_k, u_k, k, Bk, inctrans_k);
		end
	end
end