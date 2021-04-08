classdef ComparisonDecomp < handle
	properties
		p0;
		p1;

		stats0;
		stats1;

		income;
		nthresholds;

		grids;

		income;
		results;
	end

	methods
		function obj = ComparisonDecomp(s0, s1)
			obj.p0 = s0.p;
			obj.p1 = s1.p;
			
			obj.income = s0.income;

			if obj.p0.OneAsset
				obj.grids = {s0.bgrid};
			else
				obj.grids = {s0.bgrid, s0.agrid};
			end

			obj.nthresholds = numel(p0.decomp_thresholds);

			obj.stats0 = s0.stats;
			obj.stats1 = s1.stats;

			obj.initialize();

			if isequaln(s0, s1)
		        return
		    end

		    obj.perform_decomp();
		end

		function initialize(obj)
			nan_vec = NaN(obj.nthresholds, 1);
			obj.results = struct('Em1_less_Em0', NaN,...
				'term1', NaN, 'term2', NaN,...
				'term3', NaN, 'term2a', nan_vec,..
				'term2b', nan_vec, 'term2c', nan_vec);
		end

		function perform_decomp(obj)
			reshape_dims = [obj.p0.nb_KFE, obj.p0.na_KFE,...
				obj.p0.nz*obj.income0.ny];

			m0 = reshape(obj.stats0.mpcs_over_ss{5}, reshape_dims);
		    pmf0 = obj.stats0.pmf;
		    [m0_x, pmf0_x] = aux.collapse_mpcs(m0, pmf0);
		    Em0 = dot(m0(:), pmf0(:));

		    reshape_dims = [obj.p1.nb_KFE, obj.p1.na_KFE,...
				obj.p1.nz*obj.income0.ny];
		    m1 = reshape(stats1.mpcs_over_ss{5}, reshape_dims);
		    pmf1 = obj.stats1.pmf;
		    [m1_x, pmf1_x] = aux.collapse_mpcs(m1, pmf1);
		    Em1 = dot(m1(:), pmf1(:));
		end
	end
end