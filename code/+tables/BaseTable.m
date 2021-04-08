classdef BaseTable < handle
	properties (SetAccess = protected)
		mpcs_present = false;
		illiquid_mpcs_present = false;
		mpcs_news_present = false;
		one_asset_only = false;
		two_asset_only = false;
		decomp_baseline_present = false;
		decomp_norisk_present = false;
		selected_cases;

		n_cols;

		current_column;
	end

	properties
		output;
		decomp_incrisk_alt;
	end

	methods
		function obj = BaseTable(params, stats)
			obj.set_options(params, stats);
			% obj.filter_experiments(params);
		end

		function set_options(obj, params, stats)
			obj.mpcs_present = any([params.ComputeMPCS]);
			obj.illiquid_mpcs_present = any([params.ComputeMPCS_illiquid]);
			obj.mpcs_news_present = any([params.ComputeMPCS_news]);
			obj.decomp_norisk_present = any(...
				cellfun(@(x) x.decomp_norisk_completed, stats));

			for ii = 1:numel(stats)
				if stats{ii}.decomp_baseline_present
					obj.decomp_baseline_present = true;
					break
				end
			end
			obj.one_asset_only = all([params.OneAsset]);
			obj.two_asset_only = all(~[params.OneAsset]);

			obj.n_cols = numel(params);
		end

		% function filter_experiments(obj, params, use_all)
		% 	if nargin < 3
		% 		use_all = true;
		% 	end

		% 	if use_all || isempty(obj.included_groups)
		% 		obj.selected_cases = 1:numel(params);
		% 	else
		% 		all_names = {params.group};
		% 		inames = [];
		% 		for ii = 1:numel(all_names)
		% 			if ismember(obj.included_groups, all_names{ii})
		% 				inames = [inames, ii];
		% 			end
		% 		end

		% 		obj.selected_cases = unique(inames);
		% 	end

		% 	obj.n_cols = numel(obj.selected_cases);
		% end

		function update_current_column(obj, table_in, stats_in)
			if nargin < 3
				tmp = table_in;
			else
				tmp = obj.construct_from_stats(table_in, stats_in);
			end
			obj.current_column = [obj.current_column; tmp];
		end

		function table_out = construct_from_stats(obj, table_in, stats_in)
			vals = {};
			labels = {};
			jj = 1;
			for ii = 1:numel(stats_in)
				if obj.one_asset_only && (stats_in{ii}.indicator == 2)
					continue
				elseif obj.two_asset_only && (stats_in{ii}.indicator == 1)
					continue
				end

				vals{jj} = stats_in{ii}.value;
				labels{jj} = stats_in{ii}.label;
				jj = jj + 1;
			end

			table_to_append = table(vals(:),...
				'VariableNames', {'results'},...
				'RowNames', labels(:));

			table_out = [table_in; table_to_append];
		end

		function add_column(obj, ip)
			column_label = sprintf('Specification%d', ip);
			obj.current_column.Properties.VariableNames = {column_label};
			obj.output = [obj.output, obj.current_column];
		end

		function new_table = new_table_with_header(obj, header_name)
			header_formatted = strcat('____', header_name);
			new_table = table({NaN},...
				'VariableNames', {'results'},...
				'RowNames', {header_formatted});
		end
	end
end