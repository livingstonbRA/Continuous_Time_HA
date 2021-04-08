function [V, gg] = make_initial_guess(p, grids, gridsKFE, income)
    % Makes initial guesses for the value function and distribution
    %
    % Parameters
    % ----------
    % p : a Params object
    %
    % grids : a Grid object which holds the HJB grids
    %
    % gridsKFE : a Grid object which holds the KFE grids
    %
    % income : an Income object
    %
    % Returns
    % -------
    % V : value function guess, of shape (nb, na, nz, ny)
    %
    % gg : distribution guess, of shape (nb_KFE, na_KFE, nz, ny)

    import aux.repmat_auto
    import aux.sparse_diags

	nb = p.nb; na = p.na; 
	nz = p.nz; ny = income.ny;
	dim = nb * na * nz * ny;

	ss_dims = [nb, na, nz, ny];
	repmat_to_ss = @(arr) repmat_auto(arr, ss_dims);
	reshape_to_ss = @(arr) reshape(arr, ss_dims);

    if numel(p.rhos) > 1
        rho_mat = p.deathrate + shiftdim(p.rhos, -2);
        rho_mat = repmat_to_ss(rho_mat);
        rho_mat = sparse_diags(rho_mat(:), 0);
    else
        rho_mat = (p.deathrate + p.rho) * speye(dim);
    end

    %% --------------------------------------------------------------------
    % GUESS FOR VALUE FUNCTION
    % ---------------------------------------------------------------------
	% Liquid returns grid
    r_b_mat = p.r_b .* (grids.b.vec >=0 ) +  p.r_b_borr .* (grids.b.vec < 0);

	% Ensure a numerically feasible guess
	r_b_mat = max(r_b_mat, 0.001);
	r_a_adj = max(p.r_a, 0.001);

	% Consumption guess
	c_0 = (1 - p.directdeposit) * (1 - p.wagetax) * income.y.wide + p.transfer ...
        + (r_a_adj + p.deathrate * p.perfectannuities) * grids.a.wide...
        + (r_b_mat + p.deathrate * p.perfectannuities) .* grids.b.vec;

    c_0 = repmat_to_ss(c_0);

    prefs = model_objects.Preferences();
    if p.SDU
        rho_bc = shiftdim(p.rhos, -2);
        prefs.set_SDU(p.invies, rho_bc + p.deathrate);
    else
        prefs.set_crra(p.invies);
    end
    u = prefs.u(c_0);

    inctrans = income.full_income_transition_matrix(p, u);
    
    if p.sigma_r > 0
        % Vaa term
        deltas = grids.a.dB + grids.a.dF;
        deltas(:, 1) = 2 * grids.a.dF(:, 1);
        deltas(:, na) = 2 * grids.a.dB(:, na);

        updiag = repmat(1 ./ grids.a.dF, [1 1 nz ny]);
        updiag(:,na,:,:) = repmat(1 ./ grids.a.dB(:,na), [1 1 nz ny]);

        centdiag = - repmat(1 ./ grids.a.dF + 1 ./ grids.a.dB, [1 1 nz ny]);
        centdiag(:,1,:,:) = repmat(1 ./ grids.a.dF(:,1), [1 1 nz ny]);
        centdiag(:,na,:,:) = -repmat(1 ./ grids.a.dB(:,na), [1 1 nz ny]);

        lowdiag = repmat(1 ./ grids.a.dB, [1 1 nz ny]);
        lowdiag(:,1,:,:) = -repmat(1 ./ grids.a.dF(:,1), [1 1 nz ny]);
        
        risk_adj = (grids.a.wide .* p.sigma_r) .^ 2;

        updiag = risk_adj .* updiag ./ deltas;
        centdiag = risk_adj .* centdiag ./ deltas;
        lowdiag = risk_adj .* lowdiag ./ deltas;

        Arisk = aux.sparse_diags(...
            [lowdiag(:), centdiag(:), updiag(:)],...
            [-nb, 0, nb]);
    else
        Arisk = sparse(dim, dim);
    end

    V = (rho_mat - inctrans - Arisk) \ u(:);
    V = reshape_to_ss(V);

    %% --------------------------------------------------------------------
    % GUESS FOR EQUILIBRIUM DISTRIBUTION
    % ---------------------------------------------------------------------
    gg0 = ones(p.nb_KFE,p.na_KFE,p.nz,income.ny);
    gg0 = gg0 .* permute( repmat(...
    	income.ydist, [1 p.nb_KFE p.na_KFE p.nz]), [2 3 4 1]);

    if p.OneAsset
        gg0(:,gridsKFE.a.vec>0,:,:) = 0;
    end
    gg0 = gg0 / sum(gg0(:));
    gg0 = gg0 ./ gridsKFE.trapezoidal.matrix;
    gg = gg0;
end
