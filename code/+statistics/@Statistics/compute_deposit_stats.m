function compute_deposit_stats(obj)
	%% CREATE LABELS
	obj.adjcosts = struct();
	obj.adjcosts.kappa0 = obj.empty_stat(...
		'kappa0, adj cost coeff on first (linear) term', 2);
	obj.adjcosts.kappa1 = obj.empty_stat(...
		'kappa1, adj cost coeff on second (power) term', 2);
	obj.adjcosts.kappa2 = obj.empty_stat(...
		'kappa2, power on second term', 2);
	obj.adjcosts.kappa_var = obj.empty_stat(...
		'kappa1 ^(-1/kappa2), low values cause mass at high a', 2);
	obj.adjcosts.a_lb = obj.empty_stat(...
		'a_lb, parameter s.t. max(a, a_lb) used for adj cost', 2);
	obj.adjcosts.adj_cost_fn = obj.empty_stat(...
		'Adjustment cost function', 2);
	obj.adjcosts.mean_cost = obj.empty_stat(...
        'Mean adjustment cost, E[chi(d, a)]', 2);
	obj.adjcosts.mean_d_div_a = obj.empty_stat(...
        'Mean ratio of deposits to assets, E[|d| / max(a, a_lb)]', 2);
	obj.adjcosts.mean_chi_div_d = obj.empty_stat(...
		'E[chi(d,a)/abs(d) | d != 0]', 2);

	npct = numel(obj.p.wpercentiles);
	obj.adjcosts.chi_div_d_pctiles = cell(1, npct);
	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);
		obj.adjcosts.chi_div_d_pctiles{ip} = obj.empty_stat(...
			sprintf('chi/abs(d), %gth pctile condl on d != 0', pct_at), 2);
	end

	%% COMPUTATIONS
	if obj.p.OneAsset
		return
	end

	obj.adjcosts.kappa0.value = obj.p.kappa0;
	obj.adjcosts.kappa1.value = obj.p.kappa1;
	obj.adjcosts.kappa2.value = obj.p.kappa2;
	obj.adjcosts.kappa_var.value = obj.p.kappa1 ^ (-1/obj.p.kappa2);
	obj.adjcosts.a_lb.value = obj.p.a_lb;

	term1 = sprintf("%g |d|", obj.p.kappa0);
	term2 = sprintf("(%g / (1 + %g)) |d / max(a,%g)| ^ (1 + %g) * max(a,%g)",...
		obj.p.kappa1, obj.p.kappa2, obj.p.a_lb, obj.p.kappa2, obj.p.a_lb);
	fn_form = strcat("cost(d,a)", " = ", term1, " + ", term2);
	obj.adjcosts.adj_cost_fn.value = fn_form;

	% Adj cost statistics
	adj_cost_obj = aux.AdjustmentCost();
	adj_cost_obj.set_from_params(obj.p);
	chii = adj_cost_obj.compute_cost(obj.model.d,...
    	shiftdim(obj.agrid, -1));
	obj.adjcosts.mean_cost.value = obj.expectation(chii);

	% Mean abs(d) / a
    d_div_a = abs(obj.model.d) ./ max(...
    	shiftdim(obj.agrid, -1), obj.p.a_lb);
    obj.adjcosts.mean_d_div_a.value = obj.expectation(d_div_a);

    % Conditional distribution of chi / |d|
    chii = adj_cost_obj.compute_cost(obj.model.d,...
    	shiftdim(obj.agrid, -1));
    chi_div_d = chii ./ abs(obj.model.d);
    nonzdeposits = abs(obj.model.d) > 1.0e-7;
    chi_div_d = chi_div_d(:);
    chi_div_d = chi_div_d(nonzdeposits(:));
    condpmf = obj.pmf(:);
    condpmf = condpmf(nonzdeposits(:));
    condpmf = condpmf / sum(condpmf);

    interpolant = aux.pctile_interpolant(chi_div_d, condpmf);
    obj.adjcosts.mean_chi_div_d.value = condpmf' * chi_div_d;
	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);
		obj.adjcosts.chi_div_d_pctiles{ip}.value = interpolant(pct_at / 100);
	end
end