function add_mpcs_news(obj, mpc_obj)
	empty_stat = obj.sfill([], []);

	% Shock in one quarter
	empty_mpc_struct = struct(...
		'shock_normalized', empty_stat,...
		'shock', empty_stat,...
		'quarterly', empty_stat...
	);

	mpcs_stats = empty_mpc_struct;
	nshocks = numel(obj.p.mpc_shocks);
	for ishock = 1:nshocks
		shock = obj.p.mpc_shocks(ishock);
		shock_label = obj.p.quantity2label(shock);

		mpcs_stats(ishock) = empty_mpc_struct;

		mpcs_stats(ishock).shock_normalized = obj.sfill(...
			shock * 100, 'Size of shock next quarter, (% of mean ann inc)');

		mpcs_stats(ishock).shock = obj.sfill(shock_label,...
			'Size of shock next quarter');

		tmp = 100 * mpc_obj.mpcs(ishock).avg_1_quarterly;
		label = sprintf(...
			'Quarterly MPC (%%), out of %s next quarter',...
			shock_label);
		mpcs_stats(ishock).quarterly = obj.sfill(tmp, label);
	end
	obj.mpcs_news_one_quarter = mpcs_stats;

	% Shock in one year
	empty_mpc_struct = struct(...
		'shock_normalized', empty_stat,...
		'shock', empty_stat,...
		'quarterly', empty_stat,...
		'annual', empty_stat...
	);

	mpcs_stats = empty_mpc_struct;
	nshocks = numel(obj.p.mpc_shocks);
	for ishock = 1:nshocks
		shock = obj.p.mpc_shocks(ishock);
		shock_label = obj.p.quantity2label(shock);

		mpcs_stats(ishock) = empty_mpc_struct;

		mpcs_stats(ishock).shock_normalized = obj.sfill(...
			shock * 100, 'Size of shock next year, (% of mean ann inc)');

		mpcs_stats(ishock).shock = obj.sfill(shock_label,...
			'Size of shock next year');

		tmp = 100 * mpc_obj.mpcs(ishock).avg_4_quarterly;
		label = sprintf(...
			'Quarterly MPC (%%), out of %s next year',...
			shock_label);
		mpcs_stats(ishock).quarterly = obj.sfill(tmp, label);

		tmp = 100 * mpc_obj.mpcs(ishock).avg_4_annual;
		label = sprintf(...
			'Annual MPC (%%), out of %s next year',...
			shock_label);
		mpcs_stats(ishock).annual = obj.sfill(tmp, label);
	end
	obj.mpcs_news_one_year = mpcs_stats;
end