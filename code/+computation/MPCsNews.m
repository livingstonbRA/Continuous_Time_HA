classdef MPCsNews < handle
	% This class provides properties and methods 
	% for solving for policy functions given a future shock
	% by iterating backward on the dynamic HJB.
	%
	% MPCs out of news can be computed with this class
	% by iterating over the Feynman-Kac equation while
	% iterating over the dynamic HJB.
	%
	% Allows for saving the policy functions 's', 'c', and
	% 'd' at given time steps for simulating MPCs out of news.

	properties (Constant)
		% Default option values.
		defaults = struct(...
					'delta', 0.025,...
					'delta_terminal', 1e-3,...
					'save_policies', false,...
					'compute_mpcs', true...
					);
	end

	properties (SetAccess = protected)

		% parameters, grids of the model
		p;
		income;
		grids;

		% Value function at given instant in time.
		V;

		% Transition matrix constructor.
		A_constructor_HJB;

		% Transition matrix constructor, for A
		% without returns risk (only used if returns risk
		% is on but turned off for the KFE).
		A_constructor_FK;

		% Transition matrix at given instant in time.
		A_HJB;

		% Transition matrix at given instant in time,
		% without returns risk (only used if returns risk
		% is on but turned off for the KFE).
		A_FK;

		% Structure containing the policy functions at
		% given instance in time.
		KFEint;
        
        % Row vector containing the shock indices to be
        % used.
        shocks;

        % Vector containing the "time until shock" values
        % at which point the policy functions are to be saved.
        % Only used if simulating MPCs out of news.
        savedTimesUntilShock;

        % Cumulative consumption over states.
		cumcon;

		% Income transitions with the diagonal removed.
		ytrans_offdiag;

		% Cumulative consumption for the q1, q4 shock.
		cum_con_q1 = cell(1,6);
		cum_con_q4 = cell(1,6);

		% Results structure.
		mpcs = struct();

		% A sparse diagonal matrix with rho values on
		% the diagonal.
		rho_diag;

		% Total number of states.
		n_states;

		% States in each income block.
		states_per_income;

		% Row vector with dimensions of state space.
		state_dims;

		% Options used by the solver.
		options;

	end

	methods
		function obj = MPCsNews(p, income, grids, shocks, varargin)
			% class constructor

			% Parameters
			% ----------
			% p : a Params object
			%
			% income : an Income object
			%
			% grids : a Grid object
			%
			% shocks : a vector of shock indices, in reference to the mpc shocks
			%	in 'params'
			%
			% Returns
			% -------
			% obj : a MPCsNews object

			% store references to important model objects
			obj.p = p;
			obj.income = income;
			obj.grids = grids;
            obj.shocks = shocks;

            obj.options = parse_options(varargin{:});

            obj.states_per_income = p.nb_KFE * p.na_KFE * p.nz;
            obj.n_states = obj.states_per_income * income.ny;
            obj.state_dims = [p.nb_KFE p.na_KFE p.nz income.ny];

			obj.ytrans_offdiag = income.ytrans - diag(diag(income.ytrans));

			% Initialize the results
			obj.mpcs = struct();
			for ishock = 1:6
				obj.mpcs(ishock).avg_1_quarterly = NaN;
				obj.mpcs(ishock).avg_4_quarterly = NaN(4,1);
				obj.mpcs(ishock).avg_4_annual = NaN;
            end

            import computation.TransitionMatrixConstructor

            % Initialize constructors of the A matrix
            returns_risk = (obj.p.sigma_r > 0);
		    obj.A_constructor_HJB = TransitionMatrixConstructor(...
                	obj.p, obj.income, obj.grids, returns_risk);

		    % If returns are risky and returns risk is used in the KFE,
		    % this A constructor will be used instead:
		    if returns_risk && (~obj.p.retrisk_KFE)
		    	obj.A_constructor_FK = TransitionMatrixConstructor(...
		    		obj.p, obj.income, obj.grids, false);
		    end

		    if numel(obj.p.rhos) > 1
		        rhocol = kron(obj.p.rhos', ones(obj.p.nb_KFE*obj.p.na_KFE, 1));
		        obj.rho_diag = aux.sparse_diags(rhocol, 0);
		    else
		        obj.rho_diag = obj.p.rho * speye(obj.states_per_income);
			end
		end

		function solve(obj, KFE, pmf, cum_con_baseline)
			% calls various class methods to solve the model backward
			% and compute MPCs

			% Parameters
			% ----------
			% KFE : a structure containing baseline policy functions and
			%	the value function on the KFE grid
			%
			% pmf : the equilibrium distribution over the KFE grid, of
			%	shape (nb_KFE, na_KFE, nz, ny)
			%
			% cum_con_baseline : cumulative consumption for each period
			%	in the baseline, used to compute MPCs
			%
			% Directly Modifies
			% -----------------
			% obj.savedTimesUntilShock : a vector indicating the time until
			%	the shock at a given iteration, used to find the appropriate
			%	policy function given the current subperiod; only modified/used
			%	if SimulateMPCS_news is turned on
			%
			% Note
			% ----
			% This method calls various other class methods and so indirectly
			% modifies other class properties.
            
            if obj.options.save_policies
				savedTimesUntilShock = [4:-obj.options.delta:obj.options.delta];
				round_factor = 1 / obj.options.delta;
				obj.savedTimesUntilShock = round(savedTimesUntilShock*round_factor)/round_factor;
            end
            
            % loop over shocks
			for ishock = obj.shocks
                fprintf('    --- Shock = %f ---\n', obj.p.mpc_shocks(ishock))
				success = obj.getTerminalCondition(KFE,ishock);
				if ~success
					return
				end
				obj.iterateBackwards(ishock);

				if obj.options.compute_mpcs
					obj.computeMPCs(pmf,ishock,cum_con_baseline);
				end
            end
            fprintf('\n')
		end

		function success = getTerminalCondition(obj, KFE, ishock)
			% this method finds the value function the instant before
			% the shock is applied

			% Parameters
			% ----------
			% KFE : structure containing the policy functions and the value
			%	function on the KFE grid
			%
			% ishock : the index of the shock, in reference to the vector
			%	of shock sizes in the Params object used to instantiate this
			%	class
			%
			% Returns
			% -------
			% success : a flag taking the value of true if the terminal
			%	condition was found, false otherwise
			%
			% Modifies
			% --------
			% obj.V : the value function given a future shock; this method
			%	initializes the value function with the terminal value
			%	function, of shape (nb_KFE, na_KFE, nz, ny)
			%
			% obj.KFEint : a structure containing the policy functions on
			%	the KFE grids; this method initializes these policy functions
			%	with the terminal policy functions
			%
			% obj.A : the transition matrix over the KFE grids; this method
			%	initializes A with the terminal transition matrix

			import aux.sparse_diags

            shock = obj.p.mpc_shocks(ishock);

            success = false;

			% Get the guess of terminal value function
			% V_{T+1}as V(b+shock,a,y)

			Vg_terminal = obj.get_value_at_instant_of_shock(KFE, shock);

			reshape_vec = [obj.p.nb_KFE,obj.p.na_KFE,obj.p.nz,obj.p.ny];
			Vg_terminal = reshape(Vg_terminal,reshape_vec);
			Vg_terminal_k = reshape(Vg_terminal,[],obj.p.ny);

			% iterate with implicit-explicit scheme to get V_terminal
		    V_terminal = Vg_terminal;
		    V_terminal_k = reshape(V_terminal, [], obj.p.ny);

		    hours = 1; % endog labor not coded for MPCs out of news
            
		    for ii = 1:5000
		    	KFE_terminal = computation.find_policies(...
		    		obj.p, obj.income, obj.grids, V_terminal, hours);
		    	A_terminal = obj.A_constructor_HJB.construct(KFE_terminal);

		    	u_k = reshape(KFE_terminal.u, [], obj.income.ny);

		    	V1_terminal_k = zeros(obj.states_per_income, obj.income.ny);
                obj.options.delta_terminal = 1e-3;
		    	for k = 1:obj.p.ny
		    		ind1 = 1+obj.states_per_income * (k-1);
			    	ind2 = obj.states_per_income * k;
		    		Ak = A_terminal(ind1:ind2, ind1:ind2);
                    
					indx_k = ~ismember(1:obj.p.ny,k);
		    		if ~obj.p.SDU
			    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k), obj.states_per_income, 1) ...
	                            .* V_terminal_k(:,indx_k),2);
			    		inctrans = obj.income.ytrans(k,k) * speye(obj.states_per_income);
			    	else
			    		ez_adj = obj.income.income_transitions_SDU(obj.p, V_terminal);
			    		Vk_stacked = sum(squeeze(ez_adj(:,k,indx_k)) .* V_terminal_k(:,indx_k), 2);
			    		inctrans = sparse_diags(ez_adj(:,k,k), 0);
			    	end

			    	RHS = obj.options.delta_terminal* (u_k(:,k) + Vk_stacked)...
			           		+ obj.options.delta_terminal*Vg_terminal_k(:,k)/obj.options.delta + V_terminal_k(:,k);
			        LHS = obj.rho_diag + (1/obj.options.delta + 1/obj.options.delta_terminal + ...
			        			obj.p.deathrate) * speye(obj.states_per_income) - Ak - inctrans;
			       	LHS = obj.options.delta_terminal * LHS;
		        	V1_terminal_k(:,k) = LHS \ RHS;
		    	end

    			dst = max(abs(V1_terminal_k(:)-V_terminal(:)));
		        if (ii==1) || (mod(ii,10) == 0)
		            fprintf('\tFinding terminal value fn for mpc out of news, iter = %d, diff = %d\n',ii,dst)
		        end
		        
		        V_terminal_k = V1_terminal_k;
                V_terminal = reshape(V_terminal_k, obj.state_dims);
		        
		        if dst < 1e-7
		        	fprintf('\tFound terminal value function after %i iterations\n',ii);
		        	success = true;
		            break
		        end
		    end

		    obj.V = V_terminal;
		    obj.KFEint = KFE_terminal;
            obj.A_HJB = A_terminal;

            if (obj.p.sigma_r > 0) && (~obj.p.retrisk_KFE)
				obj.A_FK = obj.A_constructor_FK.construct(obj.KFEint);
			end
		end

		function V_at_shock = get_value_at_instant_of_shock(obj, KFE, shock)
			% V_at_shock is needed to find the value function an instant
			% before the shock

			interp_grids = {obj.grids.b.vec, obj.grids.a.vec,...
				obj.grids.z.vec, obj.income.y.vec};
			if (obj.p.ny > 1) && (obj.p.nz > 1)
				inds = 1:4;
            elseif obj.p.ny > 1
                inds = [1, 2, 4];
            elseif obj.p.nz > 1
                inds = 1:3;
			else
				inds = 1:2;
			end

			Vinterp = griddedInterpolant(interp_grids(inds),squeeze(KFE.Vn), 'linear');
			points = {obj.grids.b.vec(:)+shock, obj.grids.a.vec(:),...
					obj.grids.z.vec(:), obj.income.y.vec(:)};
			V_at_shock = Vinterp(points(inds));
		end

		function iterateBackwards(obj, ishock)
			% iterate over the dynamic HJB, solver for policy functions, etc...

			% Parameters
			% ----------
			% ishock : the shock index, in reference to the shock sizes in the
			%	Params object used to instantiate this class
			%
			% Directly Modifies
			% -----------------
			% obj.cum_con_q1 : a cell array indexed by shock, where each entry
			%	contains the cumulative consumption over states given a shock in
			%	period 1 (where initial period is taken to be 0)
			%
			% obj.cum_con_q4 : a cell array indexed by shock, where each entry
			%	contains the cumulative consumption over states given a shock in
			%	period 4 (where initial period is taken to be 0)
			%
			% Note
			% ----
			% This method calls various other class methods and so indirectly
			% modifies other class properties.

			import computation.hjb_divisor
			import aux.sparse_diags
            import computation.FeynmanKac
			
            shock = obj.p.mpc_shocks(ishock);
			obj.cumcon = zeros(obj.n_states, 4);

			for it = 4:-obj.options.delta:obj.options.delta
				timeUntilShock = obj.options.delta + 4 - it;
				round_factor = 1 / obj.options.delta;
                timeUntilShock = round(timeUntilShock * round_factor) / round_factor;
				if mod(it,0.5) == 0
		            fprintf('\tUpdating policy functions given news, quarter=%0.2f\n',it)
                end

                if obj.p.SDU
                	ez_adj = obj.income.income_transitions_SDU(obj.p, obj.V);
			    end
                
                u_k = reshape(obj.KFEint.u,[],obj.p.ny);
                V_k = reshape(obj.V,[],obj.p.ny);
                V_k1 = zeros(size(V_k));
                for k = 1:obj.income.ny
                    ind1 = 1+obj.states_per_income*(k-1);
			    	ind2 = obj.states_per_income*k;
                    Ak = obj.A_HJB(ind1:ind2,ind1:ind2);
                    
                    indx_k = ~ismember(1:obj.income.ny, k);

                    if ~obj.p.SDU
                    	inctrans = obj.income.ytrans(k,k) * speye(obj.states_per_income);
			    		Vk_stacked 	= sum(repmat(obj.income.ytrans(k,indx_k), [obj.states_per_income 1]) ...
	                            		.* V_k(:,indx_k), 2);
			    	else
			    		inctrans = sparse_diags(ez_adj(:,k,k), 0);
			    		Vk_stacked 	= sum(squeeze(ez_adj(:,k,indx_k)) .* V_k(:,indx_k),2);
			    	end
			    	Bk = hjb_divisor(obj.options.delta, obj.p.deathrate, k, obj.A_HJB, inctrans, obj.rho_diag);
                    V_k1(:,k) = Bk \ (obj.options.delta * (u_k(:,k) + Vk_stacked) + V_k(:,k));
                end
                obj.V = reshape(V_k1, obj.state_dims);

		        % Find policies a fraction of a period back
		        obj.update_policies();
		        
		        % Find A matrix a fraction of a period back
		        obj.update_A_matrix();

			    if obj.options.compute_mpcs
				    if (obj.p.sigma_r > 0) && (~obj.p.retrisk_KFE)
	                    FKmats = FeynmanKac.divisor(obj.p, obj.income,...
	                    			obj.options.delta, obj.A_FK, true);
	                else
	                	FKmats = FeynmanKac.divisor(obj.p, obj.income,...
	                				obj.options.delta, obj.A_HJB, true);
	                end

			        for period = ceil(it):4
			        	obj.cumcon(:,period) = FeynmanKac.update(obj.p, obj.grids,...
							obj.income, obj.cumcon(:,period), FKmats, obj.KFEint.c,...
							obj.options.delta);
	                end

			        if (it == 3 + obj.options.delta)
			        	obj.cum_con_q1{ishock} = obj.cumcon(:,4);
			        end
                end

			    if obj.options.save_policies && ismember(timeUntilShock, obj.savedTimesUntilShock)
			    	% save policy function
			    	index = find(obj.savedTimesUntilShock==timeUntilShock);
                    
                    obj.savePolicies(index, ishock);
                end
            end

            if obj.options.compute_mpcs
	            obj.cum_con_q4{ishock}(:,1) = obj.cumcon(:,1);
	            for period = 2:4
	                obj.cum_con_q4{ishock}(:,period) =...
	                	obj.cumcon(:,period) - obj.cumcon(:,period-1);
	            end
	        end
		end

		function update_policies(obj)
			% updates the policy functions during backward iteration

			% Modifies
			% --------
			% obj.KFEint : the intermediate policy functions over the
			%	KFE grids

			hours = 1;
			obj.KFEint = computation.find_policies(...
				obj.p, obj.income, obj.grids, obj.V, hours);
		end

		function update_A_matrix(obj)
			% updates the transition matrix during backward iteration

			% Modifies
			% --------
			% obj.A_HJB : the transition matrix used for the policy functions,
			%	and when returns risk is ommitted from the KFE, is also used to
			%	solve the Feynman-Kac equations
			%
			% obj.A_FK : 
			obj.A_HJB = obj.A_constructor_HJB.construct(obj.KFEint);

			if (obj.p.sigma_r > 0) && (~obj.p.retrisk_KFE)
				obj.A_FK = obj.A_constructor_FK.construct(obj.KFEint);
			end
		end

		function computeMPCs(obj, pmf, ishock, cum_con_baseline)
			% performs MPC computations

            shock = obj.p.mpc_shocks(ishock);
            
            mpcs_1_quarterly = (obj.cum_con_q1{ishock} - cum_con_baseline(:,1)) / shock;
            obj.mpcs(ishock).avg_1_quarterly = mpcs_1_quarterly(:)' * pmf(:);

            mpcs_4_quarterly = (obj.cum_con_q4{ishock} - cum_con_baseline) / shock;
            obj.mpcs(ishock).avg_4_quarterly = mpcs_4_quarterly' * pmf(:);
            obj.mpcs(ishock).avg_4_annual = sum(obj.mpcs(ishock).avg_4_quarterly);
		end

		function savePolicies(obj, index, ishock)
			% saves the policy functions for MPC simulation
			
        	c = obj.KFEint.c;
        	s = obj.KFEint.s;
        	d = obj.KFEint.d;

	    	name = sprintf('policy%ishock%i.mat', index, ishock);
	    	spath = fullfile('temp', name);
	    	save(spath, 'c', 's', 'd')
        end
	end

end

function options = parse_options(varargin)
	import computation.MPCsNews
	import aux.parse_keyvalue_pairs

	defaults = MPCsNews.defaults;
	options = parse_keyvalue_pairs(defaults, varargin{:});
end