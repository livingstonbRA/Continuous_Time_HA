classdef Grid < handle
    % This class creates and stores the asset grids.
   
	properties (SetAccess=protected)
		% Model parameters
        p;
        
        % The liquid asset grids.
		b = struct();

		% The illiquid asset grids.
		a = struct();

		% Total wealth.
		w = struct();

		% Dimension for preference/other heterogeneity.
		z = struct();

		% Number of points on the illiquid asset grid.
		na;

		% Number of points on the liquid asset grid.
        nb;

        % Number of points on the non-negative section of
        % the liquid asset grid.
        nb_pos;

        % Number of points on the negative section of
        % the liquid asset grid.
        nb_neg;

        % Linear index for b = 0 along the second dimension.
        loc0b;

		% Linear index for a = 0 along the second dimension.
        loc0a;

        % Linear index for b = 0 and a = 0 in a stacked b/a
        % vector.
        % i.e. repmat(b, na) == 0 and kron(agrid, ones(nb, 1)) == 0
        loc0b0a;

        % Grid deltas for trapezoidal integration.
        da_tilde;
        db_tilde;
        trapezoidal = struct();

        % String indicating which grid is being created, 'HJB'
        % or 'KFE'.
        gtype;

        % Number of points on the income grid.
        ny;

        % Number of states in the extra dimension of
		% heterogeneity.
        nz;

        % Grid for the values of the discount factor, rho,
        % before adding the rho.
        rho_grid = 0; % e.g. [-0.001,0,0.001]

        % A grid of rho values. In particular, the sum of
        % rho and rho_grid, a vector.
        rhos;
	end

	methods
	    function obj = Grid(p, ny, gtype)
	    	obj.p = p;

	    	% Pass gtype = 'HJB' or 'KFE'
	    	obj.gtype = gtype;
	    	obj.ny = ny;
            obj.nz = p.nz;

            % Set grid sizes from params
	    	if strcmp(gtype,'KFE')
	    		% Use different-sized KFE grid
                obj.nb = p.nb_KFE;
				obj.nb_neg = p.nb_neg_KFE;
                obj.nb_pos = p.nb_pos_KFE;
                obj.na = p.na_KFE;
	    	else
                obj.nb = p.nb;
				obj.nb_neg = p.nb_neg;
                obj.nb_pos = p.nb_pos;
                obj.na = p.na;
            end

            obj.create_agrids();
            obj.create_bgrids();
            obj.generate_variables();
        end

        function obj = generate_variables(obj)
        	obj.finite_diff('a');
		    obj.finite_diff('b');

			obj.loc0b = find(obj.b.vec==0);
			obj.loc0a = find(obj.a.vec==0);

			% location of b = 0 and a = 0 in (b,a)-space
			bLong = repmat(obj.b.vec, obj.na, 1);
			aLong = kron(obj.a.vec, ones(obj.nb, 1));
			obj.loc0b0a = find((bLong==0) & (aLong==0));
            
            obj.construct_trapezoidal_grid();
            obj.create_zgrids();
        end

	    function create_agrids(obj)
	    	% Creates the illiquid asset grids.
	    	%
	    	% Modifies
	    	% -------
	    	% obj.a.vec :  A column vector grid.
	    	%
	    	% obj.a.wide : An array grid of shape (1, na, 1, 1).
	    	%
	    	% obj.a.matrix : An array the full size of the model, (nb, na, nz, ny)

			grid_vec = model_objects.GridConstruction.create_positive_asset_grid(...
				obj.p.amax, obj.p.a_gcurv, obj.na, 'term1_curv', obj.p.agrid_term1_curv,...
				'term1_wt', obj.p.agrid_term1_weight);
            obj.a.vec = grid_vec(:);
            obj.a.wide = shiftdim(grid_vec, -1);
			assert(all(diff(obj.a.vec)>0), 'agrid not strictly increasing')
	    end

	    function create_bgrids(obj, bgrid_vec)
	    	% Creates the liquid asset grids.
	    	%
	    	% Modifies
	    	% --------
	    	% obj.b.vec : A column vector for the liquid asset grid.
	    	%
	    	% obj.b.matrix : An array of the liquid asset grid,
	    	%	shape (nb, na, nz, ny).

	    	msg = ['Number of non-negative and negative liquid grid points',...
	    			'do not sum to total number of liquid grid points'];
	    	assert(obj.nb == obj.nb_pos + obj.nb_neg, msg)
	    
			bgridpos = model_objects.GridConstruction.create_positive_asset_grid(...
				obj.p.bmax, obj.p.b_gcurv_pos, obj.nb_pos, 'term1_curv', obj.p.bgrid_term1_curv,...
				'term1_wt', obj.p.bgrid_term1_weight);

			if obj.nb_neg > 0
				bgridneg = model_objects.GridConstruction.create_negative_asset_grid(...
					obj.p.bmin, obj.p.b_gcurv_neg, obj.nb_neg);

			    obj.b.vec = [bgridneg; bgridpos];
			else
			    obj.b.vec = bgridpos;
			end

		    assert(all(diff(obj.b.vec)>0), 'bgrid not strictly increasing')
	    end

	    function obj = set(obj, varname, value)
	    	obj.(varname) = value;
	    end

	    function clean(obj)
	    	% Clears some variables to reduce object size
        	obj.trapezoidal = [];
        	obj.z = [];
        	obj.da_tilde = [];
        	obj.db_tilde = [];
        end
	end

	methods (Access=private)
	    function obj = finite_diff(obj, variable)
			% This function finds forward and backward differences.
			%
			% Parameters
			% ----------
			% variable : either the liquid asset, 'b', or illiquid, 'a'
			%
			% Modifies
			% -------
			% obj.<variable>.dF : the forward first difference of the selected grid
			%
			% obj.<variable>.dB : the backward first difference of the selected grid
			%
			% obj.d<variable>_tilde : the mean of dF and dB, adjusted at the boundaries

			if strcmp(variable,'a')
				input_grid = obj.a.vec;
				repvec = [1 obj.nb];
			elseif strcmp(variable,'b')
				input_grid = obj.b.vec;
				repvec = [1 obj.na];
			else
				error('invalid variable choice')
			end

			% Forward diff
			dF = NaN(size(input_grid));
			dF(1:end-1) = input_grid(2:end) - input_grid(1:end-1);
			dF(end) = dF(end-1);

			% Backward diff
			dB = NaN(size(input_grid));
			dB(2:end) = input_grid(2:end) - input_grid(1:end-1);
			dB(1) = dB(2);

			% Trapezoidal rule
		    d_tilde = (dF + dB)/2;
			d_tilde(1) = dF(1)/2;
			d_tilde(end) = dB(end)/2;
		    
		    dF = repmat(dF, repvec);
		    dB = repmat(dB, repvec);
		    
		    if strcmp(variable,'a')
		    	obj.da_tilde = d_tilde;
		        obj.a.dF = permute(dF,[2 1]);
		        obj.a.dB = permute(dB,[2 1]);
		    else
		    	obj.db_tilde = d_tilde;
		    	obj.b.dF = dF;
		    	obj.b.dB = dB;
		    end
		end

		function construct_trapezoidal_grid(obj)
			% Constructs the grid of db * da.
			import aux.sparse_diags

			obj.trapezoidal.vec = kron(obj.da_tilde, obj.db_tilde);

			tmp = obj.db_tilde * shiftdim(obj.da_tilde, -1);
		    obj.trapezoidal.matrix = repmat(tmp,...
		    	[1, 1, obj.nz, obj.ny]);
			obj.trapezoidal.diagm = sparse_diags(...
				obj.trapezoidal.matrix(:), 0);
		end


		function create_zgrids(obj)
			% Creates an integer index grid for the z-dimension.
			obj.z.vec = (1:obj.nz)';
			obj.z.wide = shiftdim(obj.z.vec, -2);
        end
	end
end