classdef StatsTable < tables.BaseTable
	methods
		function output_table = create(obj, params, stats)
			obj.output = table();
			for ip = 1:obj.n_cols
				p_ip = params(ip);
				stats_ip = stats{ip};

				obj.current_column = table();
				obj.intro_stats_table(p_ip, stats_ip);

				obj.income_stats_table(stats_ip);
				obj.wealth_stats_table(stats_ip);
				obj.mpc_size_table(stats_ip);
				obj.mpc_sign_table(stats_ip);
				obj.mpc_htm_table(stats_ip);

				obj.decomp_norisk_table(p_ip, stats_ip);

				shock = stats_ip.mpcs(5).shock.value;
				obj.mpc_comparison(stats_ip, shock);
				obj.mpc_comparison_pct(stats_ip, shock);

				if obj.one_asset_only
					assets = {'b'};
				else
					assets = {'w', 'b', 'a'};
				end
				obj.percentiles_tables(stats_ip, assets);
                obj.illiquid_mpcs_table(stats_ip);
                obj.adj_costs_table(stats_ip);
                obj.other_params_table(stats_ip);
                obj.other_stats_table(stats_ip);

                obj.add_column(ip)
			end

			output_table = obj.output;
		end

		function intro_stats_table(obj, p, stats)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_entries = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
				stats.beta_A
			};

			obj.update_current_column(out, new_entries);
		end

		function income_stats_table(obj, stats)
			panel_name = 'Income Statistics';
			out = obj.new_table_with_header(panel_name);

			% mean_ann_earnings.value = 1.0;
			% mean_ann_earnings.label = 'Mean gross annual earnings';
			% mean_ann_earnings.indicator = 0;

			% stdev_log_gross_earnings.value = 0.710;
			% stdev_log_gross_earnings.label = 'Stdev log ann gross earnings';
			% stdev_log_gross_earnings.indicator = 0;

			% stdev_log_net_earnings.value = 0.710;
			% stdev_log_net_earnings.label = 'Stdev log ann net earnings';
			% stdev_log_net_earnings.indicator = 0;

			new_entries = {
				stats.mean_gross_y_annual
				stats.std_log_gross_y_annual
				stats.std_log_net_y_annual
			};

			obj.update_current_column(out, new_entries);
		end

		function wealth_stats_table(obj, stats)
			panel_name = 'Wealth Statistics';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.totw
				stats.liqw
				stats.median_totw
				stats.median_liqw
				stats.sav0
				stats.constrained_liq_pct{1}
				stats.constrained_liq_pct{2}
				stats.constrained_liq_pct{3}
				stats.constrained_liq_pct{4}
				stats.constrained_liq_pct{5}
				stats.constrained_liq_pct{6}
				stats.constrained_liq_pct{7}
				stats.constrained_liq_dollars{1}
				stats.constrained_liq_dollars{2}
				stats.constrained_liq_dollars{3}
				stats.constrained_liq_dollars{4}
				stats.constrained_liq_dollars{5}
				stats.liqw_lt_ysixth
				stats.liqw_lt_ytwelfth
				stats.WHtM_over_HtM_biweekly
				stats.WHtM_over_HtM_weekly
				stats.w_top10share
				stats.w_top1share
				stats.wgini
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_size_table(obj, stats)
			panel_name = 'MPC size effects';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(4).quarterly
				stats.mpcs(6).quarterly
				stats.mpcs(4).annual
				stats.mpcs(6).annual
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_sign_table(obj, stats)
			panel_name = 'MPC sign effects';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(1).quarterly
				stats.mpcs(2).quarterly
				stats.mpcs(3).quarterly
				stats.mpcs(1).annual
				stats.mpcs(2).annual
				stats.mpcs(3).annual
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_htm_table(obj, stats)
			panel_name = 'MPC for HtM households';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(5).quarterly_htm
				stats.mpcs(5).quarterly_whtm
				stats.mpcs(5).quarterly_phtm
				stats.mpcs(5).annual_htm
				stats.mpcs(5).annual_whtm
				stats.mpcs(5).annual_phtm
			};

			obj.update_current_column(out, new_entries);
		end

		function decomp_norisk_table(obj, p, stats)
			if (obj.two_asset_only) || (~obj.decomp_norisk_present)
				return
			end

			panel_name = 'Decomps of E[MPC] wrt RA and no inc risk, $500 shock';
			out = obj.new_table_with_header(panel_name);

			tmp = stats.mpcs(5).quarterly;
			tmp.label = 'Quarterly MPC (%)'
			new_entries = {
				tmp
				stats.decomp_norisk(1).term1_pct
			};
			obj.update_current_column(out, new_entries);

			for ithresh = 1:numel(p.decomp_thresholds)
				threshold = p.decomp_thresholds(ithresh);
				panel_name = sprintf('For HtM threshold #%d', ithresh);
				out = obj.new_table_with_header(panel_name);

				new_entries = {
					stats.decomp_norisk(ithresh).term2;
					stats.decomp_norisk(ithresh).term3;
					stats.decomp_norisk(ithresh).term4;
				};

				obj.update_current_column(out, new_entries);
			end
		end

		function mpc_comparison(obj, stats, shock)
			if ~obj.decomp_baseline_present
				return
			end

			panel_name = sprintf(...
				'Decomposition of E[MPC] - E[MPC_b] out of %s',...
				shock);
			panel_name = strcat(panel_name,...
				', relative to baseline');
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.decomp_baseline.mean_mpc_diff
				stats.decomp_baseline.term1
				stats.decomp_baseline.term2
				stats.decomp_baseline.term2a(3)
				stats.decomp_baseline.term2b(3)
				stats.decomp_baseline.term2c(3)
				stats.decomp_baseline.term3
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_comparison_pct(obj, stats, shock)
			if ~obj.decomp_baseline_present
				return
			end

			panel_name = ...
				'Decomposition as % of E[MPC] - E[MPC_b]';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.decomp_baseline.term1_pct
				stats.decomp_baseline.term2_pct
				stats.decomp_baseline.term2a_pct(3)
				stats.decomp_baseline.term2b_pct(3)
				stats.decomp_baseline.term2c_pct(3)
				stats.decomp_baseline.term3_pct
			};

			obj.update_current_column(out, new_entries);
		end

		function percentiles_tables(obj, stats, assets)
			for ia = 1:numel(assets)
				switch assets{ia}
					case 'w'
						panel_name = 'Wealth percentiles';
						new_entries = stats.wpercentiles(:);
					case 'b'
						panel_name = 'Liquid wealth percentiles';
						new_entries = stats.lwpercentiles(:);
					case 'a'
						panel_name = 'Illiquid wealth percentiles';
						new_entries = stats.iwpercentiles(:);
				end
				out = obj.new_table_with_header(panel_name);
				obj.update_current_column(out, new_entries);
			end
		end

		function illiquid_mpcs_table(obj, stats)
			if obj.one_asset_only
				return
			end

			panel_name = 'Illiquid mpcs';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.illiquid_mpcs(1).quarterly
				stats.illiquid_mpcs(2).quarterly
				stats.illiquid_mpcs(3).quarterly
				stats.illiquid_mpcs(4).quarterly
				stats.illiquid_mpcs(5).quarterly
				stats.illiquid_mpcs(6).quarterly
				stats.illiquid_mpcs(4).annual
				stats.illiquid_mpcs(5).annual
				stats.illiquid_mpcs(6).annual
			};

			obj.update_current_column(out, new_entries);
		end

		function adj_costs_table(obj, stats)
			if obj.one_asset_only
				return
			end

			cost_lhs = 'cost(d,a)';
			cost_rhs = 'k0 * |d| + k1 * |d| ^ (1 + k2) / (1 + k2)';
			cost_fn = strcat(cost_lhs, ' = ', cost_rhs);

			panel_name = 'Adjustment cost statistics';
			panel_header = strcat(panel_name, ', ', cost_fn);
			out = obj.new_table_with_header(panel_header);

			new_entries = {
				stats.adjcosts.kappa0
				stats.adjcosts.kappa1
				stats.adjcosts.kappa2
				stats.adjcosts.kappa_var
				stats.adjcosts.a_lb
				stats.adjcosts.mean_cost
				stats.adjcosts.mean_d_div_a
				stats.adjcosts.mean_chi_div_d
				stats.adjcosts.chi_div_d_pctiles{1}
				stats.adjcosts.chi_div_d_pctiles{2}
				stats.adjcosts.chi_div_d_pctiles{3}
				stats.adjcosts.chi_div_d_pctiles{4}
				stats.adjcosts.chi_div_d_pctiles{5}
				stats.adjcosts.chi_div_d_pctiles{7}
			};

			obj.update_current_column(out, new_entries);
		end

		function other_params_table(obj, stats)
			if obj.one_asset_only
				return
			end

			panel_name = 'Other parameters';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.params.group_num
				stats.params.r_b
				stats.params.r_a
				stats.params.deathrate
				stats.params.bequests
				stats.params.bmax
				stats.params.amax
				stats.params.borrowlim
				stats.params.riskaver
				stats.params.numeraire
				stats.params.income_descr
			};

			obj.update_current_column(out, new_entries);
		end

		function other_stats_table(obj, stats)
			panel_name = 'Other statistics';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.constrained_illiq_pct{1}
				stats.constrained_illiq_pct{2}
				stats.constrained_illiq_pct{3}
				stats.constrained_illiq_pct{4}
				stats.constrained_illiq_pct{5}
				stats.constrained_illiq_pct{6}
				stats.constrained_illiq_dollars{1}
				stats.constrained_illiq_dollars{2}
				stats.constrained_illiq_dollars{3}
				stats.constrained_illiq_dollars{4}
				stats.constrained_illiq_dollars{5}
				stats.constrained_pct{1}
				stats.constrained_pct{2}
				stats.constrained_pct{3}
				stats.constrained_pct{4}
				stats.constrained_pct{5}
				stats.constrained_pct{6}
				stats.constrained_dollars{1}
				stats.constrained_dollars{2}
				stats.constrained_dollars{3}
				stats.constrained_dollars{4}
				stats.constrained_dollars{5}
				stats.w_lt_ysixth
				stats.w_lt_ytwelfth
			};

			obj.update_current_column(out, new_entries);
		end
	end
end