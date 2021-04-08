classdef MPCs < handle
	% This class contains methods that compute MPCs
	% over certain time periods. Baseline expected
	% consumption is found, and then this quantity--via
	% interpolation--is used to find expected consumption 
	% given an asset shock. Then MPCs are computed.
	%
	% Once the class is instantiated, the 'mpcs' results structure
	% is populated with NaN's. Calling solve() will compute
	% the MPCs.

	properties (Constant)
		% Default option values.
		defaults = struct(...
					'delta', 0.01,...
					'interp_method', 'linear',...
					'liquid_mpc', true,...
					'no_inc_risk', false);
    end
    
	properties (SetAccess=protected)

		% An object with at least the following attributes:
		%
		%	deathrate
		%	- The Poisson rate of death.
		%
		%	nb_KFE, na_KFE, nz
		%	- Number of points on each grid.
		p;

		% An object with at least the following attributes:
		%
		%	ytrans
		%	- The square income transition matrix, rows should
		%	  sum to zero.
		%
		%	ny
		%	- The number of income states.
		income;

		% A Grid object defined on the KFE grids.
		grids;

		% Options used internally by this class.
		options;

		% Feynman-Kac divisor matrix
		FKmat; 

		% cumulative consumption
		cumcon; % current state
		cum_con_baseline; % baseline
		cum_con_shock = cell(1,6); % shocked

		% income transitions w/o diagonal
		ytrans_offdiag;

		% results structure
		mpcs = struct();

		solved = false;
	end

	methods
		function obj = MPCs(p, income, grids, varargin)
			% class constructor

			% Required Inputs
			% ---------------
			% p : An object that satisfies the requirements
			%	laid out by the class properties.
			%	
			% income : An object that satisfies the requirements
			%	laid out by the class properties.
			%
			% grids : A Grid object on the KFE grid.
			%
			% Optional Key-Value Inputs
			% -------------------------
			% delta : Step size used while iterating over the
			%	Feynman-Kac equation, default = 0.02.
			%
			% interp_method : Interpolation method passed to
			%	griddedInterpolant, default = 'linear'. See
			%	the MATLAB documentation for griddedInterpolant.
			%
			% Returns
			% -------
			% obj : an MPCs object	

			obj.options = parse_options(varargin{:});

			obj.p = p;
			obj.income = income;
			obj.grids = grids;

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			for ii = 1:6
				obj.mpcs(ii).mpcs = NaN;
				obj.mpcs(ii).quarterly = NaN(4,1);
                obj.mpcs(ii).annual = NaN;

                obj.mpcs(ii).quarterly_htm = NaN;
                obj.mpcs(ii).annual_htm = NaN;
                obj.mpcs(ii).quarterly_whtm = NaN;
                obj.mpcs(ii).annual_whtm = NaN;
                obj.mpcs(ii).quarterly_phtm = NaN;
                obj.mpcs(ii).annual_phtm = NaN;
			end
		end

		function solve(obj, KFE, A, varargin)
			% computes the MPCs using Feynman-Kac by calling
			% a sequence of class methods

			% Parameters
			% ----------
			% KFE : a structure containing the policy functions
			%	on the KFE grids
			%
			% pmf : the equilibrium probability mass function,
			%	of shape (nb_KFE, na_KFE, nz, ny)
			%
			% A : the transition matrix for the KFE
			%
			% Directly Modifies
			% -----------------
			% obj.mpcs : the computed MPC statistics
			%
			% Note
			% ----
			% This method also modifies other class properties
			% via other methods.

			if obj.solved
				error('Already solved, create another instance')
            end

			obj.FKmat = computation.FeynmanKac.divisor(...
                obj.p, obj.income, obj.options.delta, A, true);
			obj.iterate_backward(KFE);

			if obj.options.no_inc_risk
				ishocks = 5;
			else
				ishocks = 1:6;
			end

			for ishock = ishocks
				obj.cumulative_consumption_with_shock(ishock);
				obj.computeMPCs(ishock, varargin{:});
			end

			obj.solved = true;
		end
	end

	methods (Access=private)

		function iterate_backward(obj, KFE)
			% iterates backward four quarters on the Feynman-Kac equation
			% to compute cum_con_baseline
			%
			% each quarter is split into 1/delta subperiods, which
			% are then iterated over

			% Parameters
			% ----------
			% KFE : a structure containing the policy functions and the
			%	value function, on the KFE grid
			%
			% Directly Modifies
			% -----------------
			% obj.cum_con_baseline : cumulative consumption at a given
			%	quarter for the baseline case, where the second dimension
			%	is the quarter; of shape (nb_KFE*na_KFE*nz*ny, num_quarters)
			%
			% Note
			% ----
			% This method also modifies other class properties
			% via other methods.

			import computation.FeynmanKac

			dim = obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz*obj.income.ny;
			obj.cumcon = zeros(dim,4);
			for it = 4:-obj.options.delta:obj.options.delta
				if mod(it*4,1) == 0
					fprintf('\tUpdating baseline cumulative consumption, quarter=%0.2f\n',it)
				end

				for period = ceil(it):4
					% When 'it' falls to 'period', start updating
					% that 'period'
					obj.cumcon(:,period) = FeynmanKac.update(obj.p, obj.grids,...
						obj.income, obj.cumcon(:,period), obj.FKmat,...
						KFE.c, obj.options.delta);
				end
			end

			obj.cum_con_baseline = zeros(dim,4);
			obj.cum_con_baseline(:,1) = obj.cumcon(:,1);
			for period = 2:4
				obj.cum_con_baseline(:,period) = obj.cumcon(:,period)...
					- obj.cumcon(:,period-1);
			end

		    obj.cum_con_baseline = reshape(obj.cum_con_baseline, [], 4);
		end

		function cumulative_consumption_with_shock(obj, ishock)
			% use baseline cumulative consumption to approximate
			% cumulative consumption for households presented with
			% an income shock in the first period

			% Parameters
			% ----------
			% ishock : the index of the shock, in reference to the
			%	shocks vector contained in the Params object used to
			%	instantiate this class
			%
			% Modifies
			% --------
			% obj.cum_con_shock : a cell array, indexed by shock, containing
			%	cumulative consumption over states for a given period; each
			%	cell contains an array of shape (nb_KFE*na_KFE*nz*ny, num_periods)

			shock = obj.p.mpc_shocks(ishock);

			if obj.options.liquid_mpc
				bgrid_mpc_vec = obj.grids.b.vec + shock;
				agrid_mpc_vec = obj.grids.a.vec;
			else
				bgrid_mpc_vec = obj.grids.b.vec;
				agrid_mpc_vec = obj.grids.a.vec + shock;
			end

			if (shock < 0) && obj.options.liquid_mpc
	            below_bgrid = bgrid_mpc_vec < obj.grids.b.vec(1);
	            bgrid_mpc_vec(below_bgrid) = obj.grids.b.vec(1);
	            some_below = any(below_bgrid);
	        elseif shock < 0
	        	below_agrid = agrid_mpc_vec < obj.grids.a.vec(1);
	            agrid_mpc_vec(below_agrid) = obj.grids.a.vec(1);
	            some_below = any(below_agrid);
	        end

	        % grids for interpolation
	        interp_grids = {obj.grids.b.vec, obj.grids.a.vec,...
	        	obj.grids.z.vec, obj.income.y.vec};
	       	value_grids = {bgrid_mpc_vec, agrid_mpc_vec,...
	       		obj.grids.z.vec, obj.income.y.vec};

        	if (obj.income.ny > 1) && (obj.p.nz > 1)
                inds = 1:4;
            elseif (obj.income.ny==1) && (obj.p.nz > 1)
                inds = 1:3;
            elseif obj.income.ny > 1
                inds = [1, 2, 4];
            else
                inds = [1, 2];
            end

            interp_grids = interp_grids(inds);
            value_grids = value_grids(inds);

			reshape_vec = [obj.p.nb_KFE obj.p.na_KFE obj.p.nz obj.income.ny];
			for period = 1:4
				% cumulative consumption in 'period'
	            con_period = reshape(obj.cum_con_baseline(:,period), reshape_vec);
	            mpcinterp = griddedInterpolant(interp_grids, squeeze(con_period), 'linear');

	            obj.cum_con_shock{ishock}(:,period) = reshape(mpcinterp(value_grids), [], 1);

	            if (shock < 0) && some_below && (period==1)
	                temp = reshape(obj.cum_con_shock{ishock}(:,period), reshape_vec);

	                reshape_bottom = reshape_vec;
	                if obj.options.liquid_mpc
	                	reshape_bottom(1) = 1;
	                	con_bottom = reshape(con_period(1,:,:,:), reshape_bottom);
		                temp(below_bgrid,:,:,:) = con_bottom + shock...
		                	+ obj.grids.b.vec(below_bgrid) - obj.grids.b.vec(1);
		            else
		            	reshape_bottom(2) = 1;
	                	con_bottom = reshape(con_period(:,1,:,:), reshape_bottom);
	                	a_bottom = obj.grids.a.vec(below_agrid);
	                	a_bottom = reshape(a_bottom, 1, []);
		            	temp(:,below_agrid,:,:) = con_bottom + shock...
		                	+ a_bottom - obj.grids.a.vec(1);
	            	end
	                obj.cum_con_shock{ishock}(:,period) = temp(:);                      
	            end
	        end
		end

		function computeMPCs(obj, ishock, pmf)
			% compute MPCs using the cumulative consumption arrays
			% found previously

			% Parameters
			% ----------
			% pmf : the equilibrium probability mass function of the
			%	baseline, of shape (nb_KFE, na_KFE, nz, ny)
			%
			% ishock : the shock index, in reference to the shock vector
			%	in the Params object used to instantiate this class
			%
			% Modifies
			% --------
			% obj.mpcs : the final MPC statistics computed from this class,
			%	a structure array of size nshocks
			
			shock = obj.p.mpc_shocks(ishock);

			% MPCs out of a shock at beginning of quarter 0
			mpcs = (obj.cum_con_shock{ishock} - obj.cum_con_baseline) / shock;
			if ishock == 5
				obj.mpcs(5).mpcs = mpcs;
			end

			if ~obj.options.no_inc_risk
				% Unconditional mean MPC
				obj.mpcs(ishock).quarterly = dot(mpcs(:,1), pmf(:));
				obj.mpcs(ishock).annual = dot(sum(mpcs, 2), pmf(:));

				% Conditional on HtM
				assets = obj.grids.b.vec;
				biweekly_y = shiftdim(obj.income.y.vec, -3) / 6;
				htm = assets <= biweekly_y;
				htm = repmat(htm, [1, obj.p.na_KFE, obj.p.nz, 1]);
				pmf_htm = pmf;
				pmf_htm(~htm) = 0;
				pmf_htm = pmf_htm ./ sum(pmf_htm(:));
				obj.mpcs(ishock).quarterly_htm = dot(mpcs(:,1), pmf_htm(:));
				obj.mpcs(ishock).annual_htm = dot(sum(mpcs, 2), pmf_htm(:));

				% Conditional on PHtM
				assets = obj.grids.b.vec + shiftdim(obj.grids.a.vec, -1);
				biweekly_y = shiftdim(obj.income.y.vec, -3) / 6;
				phtm = assets <= biweekly_y;
				phtm = repmat(phtm, [1, 1, obj.p.nz, 1]);
				pmf_phtm = pmf;
				pmf_phtm(~phtm) = 0;
				pmf_phtm = pmf_phtm ./ sum(pmf_phtm(:));
				obj.mpcs(ishock).quarterly_phtm = dot(mpcs(:,1), pmf_phtm(:));
				obj.mpcs(ishock).annual_phtm = dot(sum(mpcs, 2), pmf_phtm(:));

				% Conditional on WHtM
				whtm = htm & (~phtm);
				pmf_whtm = pmf;
				pmf_whtm(~whtm) = 0;
				pmf_whtm = pmf_whtm ./ sum(pmf_whtm(:));
				obj.mpcs(ishock).quarterly_whtm = dot(mpcs(:,1), pmf_whtm(:));
				obj.mpcs(ishock).annual_whtm = dot(sum(mpcs, 2), pmf_whtm(:));
			end
		end
	end
end

function options = parse_options(varargin)
	import computation.MPCs
	import aux.parse_keyvalue_pairs

	defaults = MPCs.defaults;
	options = parse_keyvalue_pairs(defaults, varargin{:});
end