classdef FeynmanKac
	methods (Static)
		function cumcon_update = update(p, grids, income, cumcon_t, FKmats, c, stepsize)
			% Computes inflows for death
			%
			% Parameters
			% ----------
			% cumcon_t_k : cumulative consumption for period t,
			%	of shape (nb_KFE*na_KFE*nz, ny)
			
			cumcon_update = zeros(size(cumcon_t));
			cumcon_t_k = reshape(cumcon_t, [], income.ny);
			for k = 1:income.ny
				indx_k = ~ismember(1:income.ny, k);
				ytrans_cc_k = sum(income.ytrans(k,indx_k) .* cumcon_t_k(:,indx_k), 2);

		        if p.Bequests
		        	deathin_cc_k = p.deathrate * cumcon_t_k(:,k);
		        else
		        	cumcon_t_z_k = reshape(cumcon_t_k, [p.nb_KFE*p.na_KFE p.nz income.ny]);
		            deathin_cc_k = p.deathrate * cumcon_t_z_k(grids.loc0b0a,:,k)';
		            deathin_cc_k = kron(deathin_cc_k, ones(p.nb_KFE*p.na_KFE, 1));
		        end

		        ind1 = 1 + p.na_KFE * p.nb_KFE * p.nz * (k-1);
		        ind2 = p.na_KFE * p.nb_KFE * p.nz * k;
		        RHS = reshape(c(:,:,:,k), [], 1) + ytrans_cc_k + deathin_cc_k ...
		                    + cumcon_t_k(:,k) / stepsize;
		        cumcon_update(ind1:ind2) = FKmats{k} * RHS;
			end
		end

		function B = divisor(p, income, stepsize, A, invert)
			% Computes the (k,k')-income blocks down the main diagonal for
			% left-hand-size of the Feynman-Kac equation.
			%
			% Required Inputs
			% ---------------
			% p : An object with the following attributes:
			%	nb_KFE, na_KFE, nz
			%	- Grid sizes.
			%
			%	deathrate
			%	- The Poisson death rate.
			%
			% income : An object with the following attributes:
			%	ytrans
			%	- The square income transition matrix. Should
			%	  have row sums of zero.
			%
			%	ny
			%	- The number of income states.
			%
			%	stepsize
			%	- The time step for iteration. Smaller values increase
			%	  the accuracy of the results but will increase
			%	  computational time.
			%
			% Optional Inputs
			% ---------------
			% invert : Boolean indicator for the return type. If
			%	true, then an object B is returned s.t. 
			%	B{k} * A is roughly equivalent to the backslash
			%	operator. If false, B is the array on the
			%	left-hand-side of the Feynman-Kac equation.

			if nargin == 4
				invert = false;
			end

			B = cell(income.ny, 1);
			states_per_income = p.nb_KFE * p.na_KFE * p.nz;
			for k = 1:income.ny
				ind1 = 1 + states_per_income * (k-1);
		    	ind2 = states_per_income * k;

		        B{k} = speye(states_per_income) * (...
		        			1 / stepsize + p.deathrate - income.ytrans(k,k)...
		                	) - A(ind1:ind2,ind1:ind2);

		        if invert
		        	B{k} = inverse(B{k});
		        end
			end
		end
	end
end