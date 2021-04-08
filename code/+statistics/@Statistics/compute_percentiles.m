function compute_percentiles(obj)
	%% CREATE LABELS
	npct = numel(obj.p.wpercentiles);
    obj.lwpercentiles = cell(1, npct);
    obj.iwpercentiles = cell(1, npct);
    obj.wpercentiles = cell(1, npct);
	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);
		obj.lwpercentiles{ip} = obj.empty_stat(...
			sprintf('b, %gth pctile', pct_at));
		obj.iwpercentiles{ip} = obj.empty_stat(...
			sprintf('a, %gth pctile', pct_at), 2);
		obj.wpercentiles{ip} = obj.empty_stat(...
			sprintf('w, %gth pctile', pct_at), 2);
	end

	obj.median_liqw = obj.empty_stat('b, median');
	obj.median_illiqw = obj.empty_stat('a, median', 2);
	obj.median_totw = obj.empty_stat('w, median', 2);

	%% COMPUTATIONS
	lw_pctile_interp = aux.pctile_interpolant(obj.bgrid, obj.pmf_b);
	iw_pctile_interp = aux.pctile_interpolant(obj.agrid, obj.pmf_a);
	w_pctile_interp = aux.pctile_interpolant(obj.wealthmat, obj.pmf_w);

	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);
		obj.lwpercentiles{ip}.value = lw_pctile_interp(pct_at/100);
		obj.iwpercentiles{ip}.value = iw_pctile_interp(pct_at/100);
		obj.wpercentiles{ip}.value = w_pctile_interp(pct_at/100);
	end

	obj.median_liqw.value = lw_pctile_interp(0.5);
	obj.median_illiqw.value = iw_pctile_interp(0.5);
	obj.median_totw.value = w_pctile_interp(0.5);
	obj.diff_median = struct('value', obj.median_totw.value - obj.median_liqw.value);
end