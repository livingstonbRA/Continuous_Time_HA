classdef Decomp < handle
	properties
		p;
		income;
		bgrid;
		nb;
		nthresholds;

		pmf;
		pmf_b;

		stats;
		m_ra;
		Empc;

		mpcs_b;

		cdf_b_interp;
		mpc_integral;

		results_norisk;
		results_RA;
	end

	methods
		function obj = Decomp(p, bgrid, stats, income)
			obj.p = p;
			obj.bgrid = bgrid;
			obj.nb = numel(bgrid);
			obj.income = income;

			obj.nthresholds = numel(p.decomp_thresholds);

			obj.pmf = reshape(stats.pmf, obj.nb, []);
			obj.pmf_b = sum(obj.pmf, 2);

			obj.stats = stats;

			obj.initialize();
		end

		function initialize(obj)
			nan_vec = NaN(1, obj.nthresholds);
			obj.results_norisk = struct('term1', nan_vec,...
				'term2', nan_vec, 'term3', nan_vec,...
				'term4', nan_vec, 'completed', false);

			obj.results_RA = struct('RAmpc', NaN,...
				'Em1_less_mRA', NaN, 'term1', NaN,...
				'term2', NaN, 'term3', NaN);
		end

		function compute(obj)
			obj.make_initial_computations();

			if any(isnan(obj.mpcs_b(:)))
				return;
			end

			obj.decomp_norisk();
			obj.decomp_RA();
		end

		function make_initial_computations(obj)
			import aux.interp_integral_alt

			r_b_adj = obj.p.r_b + obj.p.deathrate * obj.p.perfectannuities;
			obj.m_ra = (obj.p.rho + obj.p.deathrate - r_b_adj)...
				/ obj.p.riskaver + r_b_adj;

			

			% Compute mpc(b) = E[mpc(b,yP,ib) | b]
			mpcs_states = obj.stats.mpcs_over_ss{5};
			obj.mpcs_b = obj.collapse_mpcs(mpcs_states, obj.pmf);
			obj.Empc = dot(obj.mpcs_b, obj.pmf_b);

			% Interpolant for cdf(b)
			cdf_b = cumsum(obj.pmf_b);

		    obj.cdf_b_interp = griddedInterpolant(...
		    	obj.bgrid, cdf_b, 'pchip', 'nearest');

		    % Interpolant for mpc(b) * g(b)
		    obj.mpc_integral = interp_integral_alt(...
		    	{obj.bgrid}, obj.mpcs_b, obj.pmf_b);
		end

		function decomp_norisk(obj)
			import aux.interp_integral_alt

			% Compute mpc_norisk(b) = E[mpc_norisk(b,ib) | b]
			mpcs_states_nr = obj.stats.other.mpcs_nr(5).mpcs(:,1);
			mpcs_states_nr = reshape(mpcs_states_nr, obj.nb, []);
			mpcs_b_nr = mpcs_states_nr(:,1);

			if any(isnan(mpcs_b_nr(:)))
				return;
			end

			% Interpolant for mpc_norisk(b) * g(b)
		    norisk_integ_interp = interp_integral_alt(...
		    	{obj.bgrid}, mpcs_b_nr, obj.pmf_b);

		    for ia = 1:obj.nthresholds
		    	abar = obj.p.decomp_thresholds(ia);
		        
		        pbelow = obj.cdf_b_interp(abar);
		    	pabove = 1 - pbelow;
		        obj.results_norisk.term1(ia) = obj.m_ra;
		        obj.results_norisk.term2(ia) = obj.mpc_integral(abar) - obj.m_ra * pbelow ;
		        obj.results_norisk.term3(ia) = (mpcs_b_nr(:)' * obj.pmf_b(:) - norisk_integ_interp(abar))...
		        	- obj.m_ra * pabove;
		        obj.results_norisk.term4(ia) = (obj.Empc - obj.mpc_integral(abar))...
		        	- (mpcs_b_nr(:)' * obj.pmf_b(:) - norisk_integ_interp(abar));
		    end
		    obj.results_norisk.completed = true;
		end

		function decomp_RA(obj)
			% Find E[mpc|a=3.5]
			psmall = obj.pmf_b < 1e-9;
			winterp = griddedInterpolant(obj.bgrid(~psmall),...
				obj.mpcs_b(~psmall), 'pchip', 'nearest');
			mpc_atmean = winterp(obj.stats.liqw.value);

			obj.results_RA.RAmpc = obj.m_ra;
		    obj.results_RA.Em1_less_mRA = obj.Empc - obj.m_ra;
		    obj.results_RA.term1 = mpc_atmean - obj.m_ra;
		    obj.results_RA.term2 = 0;
		    obj.results_RA.term3 = (obj.Empc - obj.m_ra) - (mpc_atmean - obj.m_ra);
		end

		function mpcs_b = collapse_mpcs(obj, mpcs_states, pmf)
			mpcs_states = reshape(mpcs_states, obj.nb, []);
			mpcs_b = sum(mpcs_states .* pmf, 2)...
				./ obj.pmf_b;

			pmf_b_small = obj.pmf_b < 1e-8;
			mpcs_b(pmf_b_small) = mean(mpcs_states(pmf_b_small,:), 2);
		end
	end
end