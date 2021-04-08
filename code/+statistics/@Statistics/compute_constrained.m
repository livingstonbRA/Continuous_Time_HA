function compute_constrained(obj)
	import aux.multi_sum

	%% CREATE LABELS
	neps = numel(obj.p.epsilon_HtM);
    obj.constrained_illiq = cell(1, neps);
    obj.constrained_illiq_pct = cell(1, neps);
    obj.constrained_liq = cell(1, neps);
    obj.constrained_liq_pct = cell(1, neps);
    obj.constrained = cell(1, neps);
    obj.constrained_pct = cell(1, neps);
    for ip = 1:neps
		htm = obj.p.epsilon_HtM(ip);
		obj.constrained_liq{ip} = obj.empty_stat(...
			sprintf('b <= %g', htm));
		obj.constrained_liq_pct{ip} = obj.empty_stat(...
			sprintf('b <= %g%% mean ann inc', htm * 100));
		obj.constrained_illiq{ip} = obj.empty_stat(...
			sprintf('a <= %g', htm), 2);
		obj.constrained_illiq_pct{ip} = obj.empty_stat(...
			sprintf('a <= %g%% mean ann inc', htm * 100), 2);
		obj.constrained{ip} = obj.empty_stat(...
			sprintf('w <= %g', htm), 2);
		obj.constrained_pct{ip} = obj.empty_stat(...
			sprintf('w <= %g%% mean ann inc', htm * 100), 2);
	end

	ndollars = numel(obj.p.dollars_HtM);
	obj.constrained_illiq_dollars = cell(1, ndollars);
	obj.constrained_liq_dollars = cell(1, ndollars);
	obj.constrained_dollars = cell(1, ndollars);
	for ip = 1:ndollars
		obj.constrained_illiq_dollars{ip} = obj.empty_stat(...
			sprintf('a <= $%g', obj.p.dollars_HtM(ip)), 2);
		obj.constrained_liq_dollars{ip} = obj.empty_stat(...
			sprintf('b <= $%g', obj.p.dollars_HtM(ip)));
		obj.constrained_dollars{ip} = obj.empty_stat(...
			sprintf('w <= $%g', obj.p.dollars_HtM(ip)), 2);
	end

	obj.liqw_lt_ysixth = obj.empty_stat('b_i <= y_i / 6');
   	obj.liqw_lt_ytwelfth = obj.empty_stat('b_i <= y_i / 12');
   	obj.w_lt_ysixth = obj.empty_stat('w_i <= y_i / 6', 2);
   	obj.w_lt_ytwelfth = obj.empty_stat('w_i <= y_i / 12', 2);

   	obj.WHtM_over_HtM_biweekly = obj.empty_stat(...
		'P(WHtM) / P(HtM), HtM in terms of y/6', 2);
	obj.WHtM_over_HtM_weekly = obj.empty_stat(...
		'P(WHtM) / P(HtM), HtM in terms of y/12', 2);

	%% COMPUTATIONS
	lw_interp = griddedInterpolant(obj.bgrid, cumsum(obj.pmf_b), 'pchip', 'nearest');
	iw_interp = griddedInterpolant(obj.agrid, cumsum(obj.pmf_a), 'pchip', 'nearest');

	sorted_mat = sortrows([obj.wealthmat(:), obj.pmf_w(:)]);
	[w_u, iu] = unique(sorted_mat(:,1), 'last');
	cdf_w = cumsum(sorted_mat(:,2));
	w_interp = griddedInterpolant(w_u, cdf_w(iu), 'pchip', 'nearest');

	% Share with wealth less than value (units of mean ann inc)
    for ip = 1:numel(obj.p.epsilon_HtM)
		htm = obj.p.epsilon_HtM(ip);
		obj.constrained_liq{ip}.value = lw_interp(htm);
		obj.constrained_liq_pct{ip}.value = lw_interp(htm);
		obj.constrained_illiq{ip}.value = iw_interp(htm);
		obj.constrained_illiq_pct{ip}.value = iw_interp(htm);
		obj.constrained{ip}.value = w_interp(htm);
		obj.constrained_pct{ip}.value = w_interp(htm);
	end

	% Share with wealth less than value (units of dollars)
	if ~isempty(obj.p.numeraire_in_dollars)
		for ip = 1:numel(obj.p.dollars_HtM)
			htm = obj.p.dollars_HtM(ip) / obj.p.numeraire_in_dollars;
			obj.constrained_illiq_dollars{ip}.value = iw_interp(htm);
			obj.constrained_liq_dollars{ip}.value = lw_interp(htm);
			obj.constrained_dollars{ip}.value = w_interp(htm);
		end
	end

	% Share with liquid wealth / (quarterly earnings) < value
	kernel_options = struct();
	kernel_options.ktype = 'triweight';
    kernel_options.h = 0.3;
    kernel_options.rescale_and_log = true;

	by_ratio = obj.bgrid ./ obj.income.y.wide;
	pmf_by = multi_sum(obj.pmf, [2, 3]);
	tmp = sortrows([by_ratio(:), pmf_by(:)]);
	by_interp = computation.KernelSmoother(...
		tmp(:,1), tmp(:,2), 0.3, kernel_options);
	
	obj.liqw_lt_ysixth.value = by_interp.keval(1/6);
   	obj.liqw_lt_ytwelfth.value = by_interp.keval(1/12);

	% Wealth / (quarterly earnings) < epsilon
	wy_ratio = obj.wealthmat ./ obj.income.y.wide;
	pmf_wy = sum(obj.pmf, 3);
	tmp = sortrows([wy_ratio(:), pmf_wy(:)]);
	wy_interp = computation.KernelSmoother(...
		tmp(:,1), tmp(:,2), 0.3, kernel_options);
	
	obj.w_lt_ysixth.value = wy_interp.keval(1/6);
	obj.w_lt_ytwelfth.value = wy_interp.keval(1/12);

	% HtM Ratios
	obj.WHtM_over_HtM_biweekly.value = ...
		1 - obj.w_lt_ysixth.value / obj.liqw_lt_ysixth.value;

	obj.WHtM_over_HtM_weekly.value = ...
		1 - obj.w_lt_ytwelfth.value / obj.liqw_lt_ytwelfth.value;
end