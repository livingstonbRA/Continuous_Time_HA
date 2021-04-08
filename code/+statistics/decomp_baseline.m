function decomp = decomp_baseline(s0, s1)
    % Decomposition of E[mpc1] - E[mpc0]
    % 0 is baseline
    % 1 is the experiment
    
    p0 = s0.p;
    stats0 = s0.stats;
    p1 = s1.p;
    stats1 = s1.stats;
    
    decomp.Em1_less_Em0 = NaN;
    decomp.term1 = NaN;
    decomp.term2 = NaN;
    decomp.term3 = NaN;
    decomp.term2a = NaN(numel(p0.decomp_thresholds),1);
    decomp.term2b = NaN(numel(p0.decomp_thresholds),1);
    decomp.term2c = NaN(numel(p0.decomp_thresholds),1);
    
    if isequaln(s0, s1)
        return
    end

%     bgrid = stats0.grdKFE.b.vec;
%     agrid = stats0.grdKFE.a.vec;
    bgrid = stats0.bgrid;
    agrid = stats0.agrid;
    
    ny = stats0.ny;

    reshape_dims = [p0.nb_KFE, p0.na_KFE, p0.nz*ny];

    m0 = reshape(stats0.mpcs_over_ss{5}, reshape_dims);
    pmf0 = stats0.pmf;
    [m0_x, pmf0_x] = aux.collapse_mpcs(m0, pmf0);
    Em0 = dot(m0(:), pmf0(:));

    reshape_dims = [p1.nb_KFE, p1.na_KFE, p1.nz*ny];
    m1 = reshape(stats1.mpcs_over_ss{5}, reshape_dims);
    pmf1 = stats1.pmf;
    [m1_x, pmf1_x] = aux.collapse_mpcs(m1, pmf1);
    Em1 = dot(m1(:), pmf1(:));

    if p0.OneAsset
        grids = {bgrid};
    else
        grids = {bgrid, agrid};
    end

    import aux.interp_integral_alt
    m0g0interp = interp_integral_alt(grids, m0_x, pmf0_x);
    m1g0interp = interp_integral_alt(grids, m1_x, pmf0_x);
    m0g1interp = interp_integral_alt(grids, m0_x, pmf1_x);
    m1g1interp = interp_integral_alt(grids, m1_x, pmf1_x);

    % Main decomposition
    decomp.Em1_less_Em0 = Em1 - Em0;
    decomp.term1 = dot(m1_x(:) - m0_x(:), pmf0_x(:)); 
    decomp.term2 = dot(m0_x(:), pmf1_x(:) - pmf0_x(:));
    decomp.term3 = dot(m1_x(:) - m0_x(:), pmf1_x(:) - pmf0_x(:));

    for ia = 1:numel(p0.decomp_thresholds)
        x = p0.decomp_thresholds(ia);

        if p0.OneAsset
            b0_a0 = x;
            b0_amax = x;
            bmax = inf;
        elseif ~p0.OneAsset
            b0_a0 = [x, x];
            b0_amax = [x, inf];
            bmax_amax = [inf, inf];
            bmax_a0 = [inf, x];
        end

        decomp.term2a(ia) = m0g1interp(b0_a0) - m0g0interp(b0_a0);

        if p0.OneAsset
            decomp.term2b(ia) = (m0g1interp(inf) - m0g1interp(x)) ...
                - (m0g0interp(inf) - m0g0interp(x));
        else
            decomp.term2b(ia) = (m0g1interp(b0_amax) - m0g1interp(b0_a0)) ...
                - (m0g0interp(b0_amax) - m0g0interp(b0_a0));
            decomp.term2c(ia) = (m0g1interp(bmax_amax) - m0g1interp(b0_amax)) ...
                -(m0g0interp(bmax_amax) - m0g0interp(b0_amax));
        end
    end
end