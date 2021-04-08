function decomp = decomp_nobc_norisk(baseline, no_bc)
    % Decomposition w.r.t. the no-risk model and the no-bc model

    import aux.interpolate_integral

    grdKFE = baseline.grdKFE;
    p = baseline.p;
    income = baseline.income;

    assets = grdKFE.b.vec + grdKFE.a.wide;
    assets = assets(:);

    % Initialize to NaN
    for ia = 1:numel(p.decomp_thresholds)
        decomp.term1(ia) = NaN;
        decomp.term2(ia) = NaN;
        decomp.term3(ia) = NaN;
        decomp.term4(ia) = NaN;
        decomp.term5(ia) = NaN;
    end

    % Check if required MPCs are available
    mpcs_available = (p.ComputeMPCS == 1) && (p.OneAsset == 1) && (p.SolveNoRisk == 1);
    if ~mpcs_available
        return
    end
    
    %% --------------------------------------------------------------------
    % RA model WITHOUT borrowing constraint
    % ---------------------------------------------------------------------
    RA = struct();
    tmp = p.rho + p.deathrate - p.deathrate*p.perfectannuities - p.r_b;
    RA.mpc = p.r_b + tmp / p.riskaver;

    %% --------------------------------------------------------------------
    % RA model WITH borrowing constraint
    % ---------------------------------------------------------------------
    RA_with_BC = struct();

    mpcs = baseline.stats.mpcs_nr(5).mpcs(:,1);
    mpcs = reshape(mpcs, [], p.nz);
    RA_with_BC.mpcs = repmat(mpcs, 1, income.ny);

    %% --------------------------------------------------------------------
    % HA model WITH borrowing constraint (baseline)
    % ---------------------------------------------------------------------
    HA_with_BC = struct();

    HA_with_BC.Empc = baseline.stats.mpcs(5).avg_0_quarterly(1);
    HA_with_BC.mpcs = reshape(...
        baseline.stats.mpcs(5).mpcs(:,1), [], income.ny*p.nz);
    HA_with_BC.pmf = reshape(baseline.stats.pmf, [], income.ny*p.nz);

    % Interpolant for cdf
    HA_with_BC.cdf_interp = cell(income.ny*p.nz, 1);
    for ii = 1:income.ny*p.nz
        sortedAssetDistribution = sortrows([assets(:) HA_with_BC.pmf(:,ii)]);
        [assetVals, uniqueInds] = unique(sortedAssetDistribution(:,1), 'last');
        cumg = cumsum(sortedAssetDistribution(:,2));
        HA_with_BC.cdf_interp{ii} = griddedInterpolant(assetVals, cumg(uniqueInds), 'linear');
    end

    %% --------------------------------------------------------------------
    % HA model WITHOUT borrowing constraint (loose borrowing constraint)
    % ---------------------------------------------------------------------
    HA = struct();

    mpcs = reshape(...
        no_bc.stats.mpcs(5).mpcs(:,1), no_bc.p.nb_KFE, no_bc.p.na_KFE, income.ny*p.nz);

    non_negative_pts = no_bc.p.nb_neg_KFE+1:no_bc.p.nb_KFE;
    mpcs = mpcs(non_negative_pts, :, :);
    HA.mpcs = reshape(mpcs, [], income.ny*p.nz);

    %% --------------------------------------------------------------------
    % Prepare interpolants for integrals over MPCs from zero to eps
    % ---------------------------------------------------------------------
    % For integral over mpcs for HA model with budget constraint
    HA_with_BC.mpc_integral_interp = cell(income.ny*p.nz, 1);
    for ii = 1:income.ny*p.nz
        HA_with_BC.mpc_integral_interp{ii} = interpolate_integral(...
            assets, HA_with_BC.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    % For integral over mpcs for HA model without budget constraint
    HA.mpc_integral_interp = cell(income.ny*p.nz, 1);
    for ii = 1:income.ny*p.nz
        HA.mpc_integral_interp{ii} = interpolate_integral(...
            assets, HA.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    % For integral over mpcs for RA model with budget constraint
    RA_with_BC.mpc_integral_interp = cell(income.ny*p.nz, 1);
    for ii = 1:income.ny*p.nz
        RA_with_BC.mpc_integral_interp{ii} = interpolate_integral(...
            assets, RA_with_BC.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    %% --------------------------------------------------------------------
    % Compute expectations taken wrt stationary distribution of baseline
    % ---------------------------------------------------------------------
    RA_with_BC.Empc =  HA_with_BC.pmf(:)' * RA_with_BC.mpcs(:);
    HA.Empc = HA_with_BC.pmf(:)' * HA.mpcs(:);

    %% --------------------------------------------------------------------
    % Decomposition
    % ---------------------------------------------------------------------
    for ia = 1:numel(p.decomp_thresholds)
        threshold = p.decomp_thresholds(ia);

        % Precompute P(b<threshold)
        p_htm = 0;
        for ii = 1:income.ny*p.nz
            p_htm = p_htm + HA_with_BC.cdf_interp{ii}(threshold);
        end

        % Term 1: RA MPC
        decomp.term1(ia) = RA.mpc;

        % Term 2: HtM effect
        mpc_integral = 0;
        for ii = 1:income.ny*p.nz
            mpc_integral = mpc_integral...
                + HA_with_BC.mpc_integral_interp{ii}(threshold);
        end
        decomp.term2(ia) = mpc_integral - RA.mpc * p_htm;

        % Term 3: Borrowing constraint effect
        mpc_integral = RA_with_BC.Empc;
        for ii = 1:income.ny*p.nz
            mpc_integral = mpc_integral...
                - RA_with_BC.mpc_integral_interp{ii}(threshold);
        end
        decomp.term3(ia) = mpc_integral - RA.mpc * (1-p_htm);

        % Term 4: Income risk effect
        mpc_integral = HA.Empc;
        for ii = 1:income.ny*p.nz
            mpc_integral = mpc_integral...
                - HA.mpc_integral_interp{ii}(threshold);
        end
        decomp.term4(ia) = mpc_integral - RA.mpc * (1-p_htm);

        % Term 5: Interaction
        decomp.term5(ia) = HA_with_BC.Empc - decomp.term1(ia)...
            - decomp.term2(ia) - decomp.term3(ia) - decomp.term4(ia);
    end
end

function interpolant = interpolate_integral(grid_values, integrand_values, pmf, varargin)
    % Creates an interpolant that approximates the value of the integral
    % int_0^{epsilon} values(a)g(a)da for a given epsilon.
    %
    % Parameters
    % ----------
    % gridValues : Values at which the integrand is evaluated.
    %
    % integrand_values : Values of the integrand.
    %
    % pmf : The probability mass function over states.
    %
    % Results
    % -------
    % interpolant : A griddedInterpolant object such that interpolant(x)
    %   is the approximated value of the integral from 0 to x.

    is_sorted = validate_inputs(...
        grid_values, integrand_values, pmf, varargin{:});

    if ~is_sorted
        sorted_inputs = sortrows([grid_values(:) integrand_values(:) pmf(:)]);
        grid_values = sorted_inputs(:,1);
        integrand_values = sorted_inputs(:,2);
        pmf = sorted_inputs(:,3);
    end

    integral_values = cumsum(integrand_values .* pmf);
    [grid_unique, unique_inds] = unique(grid_values, 'last');
    integral_unique = integral_values(unique_inds);

    interpolant = griddedInterpolant(...
        grid_unique, integral_unique, 'pchip', 'nearest');
end

function is_sorted = validate_inputs(grid_values, integrand_values, pmf , varargin)
    assert(numel(grid_values(:)) == numel(integrand_values(:)),...
        "Inputs have inconsistent shapes");
    assert(numel(grid_values(:)) == numel(pmf(:)),...
        "Inputs have inconsistent shapes");
    assert(~(isempty(grid_values) || isempty(integrand_values) || isempty(pmf)),...
        "One or more imputs is empty");

    is_sorted = false;
    if ~isempty(varargin)
        if islogical(varargin{1})
            is_sorted = varargin{1};
        end
    end
end