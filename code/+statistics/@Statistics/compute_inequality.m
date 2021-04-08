function compute_inequality(obj)
	import aux.interp_integral_alt
	import aux.direct_gini
	import aux.multi_sum

	%% CREATE LABELS
	obj.w_top10share = obj.empty_stat('w, Top 10% share', 2);
	obj.w_top1share = obj.empty_stat('w, Top 1% share', 2);
	obj.lw_top10share = obj.empty_stat('b, Top 10% share');
	obj.lw_top1share = obj.empty_stat('b, Top 1% share');
	obj.iw_top10share = obj.empty_stat('a, Top 10% share', 2);
	obj.iw_top1share = obj.empty_stat('a, Top 1% share', 2);
	obj.wgini = obj.empty_stat('Gini coefficient, wealth');

	%% COMPUTATIONS
	% Top wealth shares
	tmp = sortrows([obj.wealthmat(:), obj.pmf_w(:)]);
	values_w = cumsum(tmp(:,1) .* tmp(:,2) / obj.totw.value);
	[cdf_u, iu] = unique(cumsum(tmp(:,2)), 'last');
	wcumshare_interp = griddedInterpolant(cdf_u, values_w(iu), 'pchip', 'nearest');
	obj.w_top10share.value = 1 - wcumshare_interp(0.9);
	obj.w_top1share.value = 1 - wcumshare_interp(0.99);

	% Top liquid wealth shares
	values_b = cumsum(obj.bgrid .* obj.pmf_b / obj.liqw.value);
	[cdf_b_u, iu_b] = unique(cumsum(obj.pmf_b), 'last');
	bcumshare_interp = griddedInterpolant(cdf_b_u, values_b(iu_b), 'pchip', 'nearest');
	obj.lw_top10share.value = 1 - bcumshare_interp(0.9);
	obj.lw_top1share.value = 1 - bcumshare_interp(0.99);

	% Top illiquid wealth shares
	if ~obj.p.OneAsset
		values_a = cumsum(obj.agrid .* obj.pmf_a(:) / obj.illiqw.value);
		[cdf_a_u, iu_a] = unique(cumsum(obj.pmf_a), 'last');
		iwshare_interp = griddedInterpolant(cdf_a_u, values_a(iu_a), 'pchip', 'nearest');
		obj.iw_top10share.value = 1 - iwshare_interp(0.9);
		obj.iw_top1share.value = 1 - iwshare_interp(0.99);
	end

	% Gini coefficient
	obj.wgini.value = direct_gini(obj.wealthmat, obj.pmf_w);
end