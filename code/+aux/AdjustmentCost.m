classdef AdjustmentCost < handle
	properties
		a_lb;
		kappa0;
		kappa1;
		kappa2;
	end

	methods
		function obj = AdjustmentCost(a_lb, kappa0, kappa1, kappa2)
			if nargin > 0
				obj.set_kappas(a_lb, kappa0, kappa1, kappa2);
			end
		end

		function set_from_params(obj, p)
			obj.set_kappas(p.a_lb, p.kappa0, p.kappa1, p.kappa2);
		end

		function set_kappas(obj, a_lb, kappa0, kappa1, kappa2)
			obj.a_lb = a_lb;
			obj.kappa0 = kappa0;
			obj.kappa1 = kappa1;
			obj.kappa2 = kappa2;
		end

		function cost_out = compute_cost(obj, d, a_grid)
			a_scaled = max(a_grid, obj.a_lb);

			cost_linear = obj.kappa0 * abs(d);
			cost_concave = obj.kappa1 / (1 + obj.kappa2) ...
				* abs(d) .^ (1 + obj.kappa2)  ./ (a_scaled .^ obj.kappa2);

			cost_out = cost_linear + cost_concave;
		end

		function cost_deriv = compute_deriv(obj, d, a_grid)
			a_scaled = max(a_grid, obj.a_lb);
			d_scaled = abs(d) ./ a_scaled;

			tmp = obj.kappa0 + obj.kappa1 .* d_scaled .^ obj.kappa2;
			cost_deriv = sign(d) .* tmp;
		end

		function d_opt = opt_deposits(obj, Vb, Va, a)
			dpos_term = max(Va ./ Vb - 1 - obj.kappa0, 0);
			dneg_term = min(Va ./ Vb - 1 + obj.kappa0, 0);

			combined = (dpos_term + dneg_term) ./ obj.kappa1;

			a_scaled = max(a, obj.a_lb);
			d_opt = sign(combined) .* a_scaled ...
				.* abs(combined) .^ (1 / obj.kappa2);
		end
	end

	methods(Static)
		function adj_cost = cost(d, a_grid, p)
			% adjustment cost function chi(d)

			% Parameters
			% ----------
			% d : deposit rate
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% adj_cost : the adjustment cost, same shape as d

			d_scaled = d./max(a_grid,p.a_lb);
    		adj_cost = max(a_grid,p.a_lb) .* (p.chi0 * abs(d_scaled) + 1/(1+p.chi2) * (abs(d_scaled).^(1+p.chi2) * p.chi1^(-p.chi2)));
		end

		function chi_prime = derivative(d, a_grid, p)
			% derivative of the adjustment cost function, chi'(d)

			% Parameters
			% ----------
			% d : deposit rate
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% chi_prime : the derivative, same shape as d

			d_scaled = d ./ max(a_grid, p.a_lb);

			linear_term = sign(d) * p.chi0;
			power_term = sign(d) .* (abs(d_scaled) / p.chi1) .^ p.chi2;
			chi_prime = linear_term + power_term;
		end

		function d = derivative_inverse(y, a_grid, p)
		    % inverse of the derivative of the adjustment cost function,
		    % i.e. (chi'(d))^{-1}

		    % Parameters
			% ----------
			% y : input values at which to evaluate the functikon
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% d : the inverse of the derivative of the adjustment cost function,
			%	same shape as y

			d = sign(y) .* max(a_grid,p.a_lb) .* p.chi1 .* (abs(y) - p.chi0).^(1/p.chi2);
		end
	end
end