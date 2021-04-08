function decomp = decomp_twoasset_oneasset(s0,s1)
    % Decomposition of E[mpc1] - E[mpc0]
    % 0 is the one-asset model
    % 1 is the two-asset model
    
    % Since it must hold that b <= n, I use NaN to represent
    % points with b > n. Sums must use 'omitnan' option.
    
    % values around which to take decompositions, not incl 0
    epsilons = [0.01,0.05];
    n_eps = numel(epsilons) + 1;
    
    % check that a one-asset baseline was passed
    if isempty(s0)
        decomp = struct();
        
        decomp.term1_mpcfn = NaN;
        decomp.term2_networth = NaN;
        decomp.term3_interaction = NaN;
        
        decomp.term1a = NaN(n_eps,1);
        decomp.term1b = NaN(n_eps,1);
        decomp.term1c = NaN(n_eps,1);
        decomp.term2a = NaN(n_eps,1);
        decomp.term2b = NaN(n_eps,1);
        
        decomp.mean_mpc_twoasset = NaN;
        decomp.mean_mpc_oneasset = NaN;
        decomp.Em1_minus_Em0 = NaN;
        return
    end
    
    % unpack objects
    p0 = s0.p;
    stats0 = s0.stats;
    p1 = s1.p;
    stats1 = s1.stats;
    grdKFE0 = s0.grdKFE;
    grdKFE1 = s1.grdKFE;
    income0 = s0.income;
    
    %% --------------------------------------------------------------------
	% CREATE NEW GRIDS
	% ---------------------------------------------------------------------
    n_networth = 500;
    n_b = p1.nb_KFE;
    
    % create net worth grid
    nmin = 0;
    nmax = max(p0.bmax,p1.bmax+p1.amax);
    
    curv = 0.2;
    ngrid = linspace(0,1,n_networth)';
    ngrid = ngrid .^(1/curv);
    ngrid = nmin + (nmax - nmin) * ngrid;
    for ia = 1:n_networth-1
        if ngrid(ia+1) - ngrid(ia) < p0.min_grid_spacing
            ngrid(ia+1) = ngrid(ia) + p0.min_grid_spacing;
        else
            break
        end
	end
    
    ngrid_long = repmat(ngrid,n_b,1);
    
    dnB = zeros(n_networth,1);
    dnB(2:end) = diff(ngrid);
    dnB(1) = dnB(2);
    
    dnF = zeros(n_networth,1);
    dnF(1:end-1) = diff(ngrid);
    dnF(end) = dnF(end-1);
    
    % trapezoidal rule with endpts
    dn_trapz = (dnB + dnF) / 2;
    dn_trapz(1) = dnF(1) / 2;
    dn_trapz(end) = dnB(end) / 2;
    
    % b grid
    bgrid = grdKFE1.b.vec;
    bgrid_long = kron(bgrid,ones(n_networth,1));
    
    dbB = zeros(p1.nb_KFE,1);
    dbB(2:end) = diff(bgrid);
    dbB(1) = dbB(2);
    
    dbF = zeros(p1.nb_KFE,1);
    dbF(1:end-1) = diff(bgrid);
    dbF(end) = dbF(end-1);
    
    % trapezoidal rule 
    db_trapz = (dbB + dbF) / 2;
    db_trapz(1) = dbF(1) / 2;
    db_trapz(end) = dbB(end) / 2;
    
    dnb = dn_trapz * db_trapz';

    % adjust dnb to ignore area below n = b line
    for ib = 2:n_b
        [~,imax] = max(ngrid-bgrid(ib)>=0,[],1);
        dnb(imax,ib) = dnb(imax,ib) ...
            * (1/2+(ngrid(imax)-bgrid(ib)) / (dn_trapz(imax)/2)/2);
    end
    
    %% --------------------------------------------------------------------
	% TWO ASSET CASE
	% ---------------------------------------------------------------------
    pmf1 = stats1.pmf;
    
    % integrate out y, z to get g(b,a)
    pmf1_b_a = sum(sum(pmf1,4),3);
    volumes = reshape(grdKFE1.trapezoidal.vec,p1.nb_KFE,p1.na_KFE);
    pdf1_b_a = pmf1_b_a ./ volumes;
    
    % interpolate to get distr over (n,b)
    grids = {grdKFE1.a.vec,grdKFE1.b.vec};
    pdf1_a_b = permute(pdf1_b_a,[2 1]);
    distr_interp = griddedInterpolant(grids,pdf1_a_b,'linear','none');
    pdf1_n_b = distr_interp(ngrid_long-bgrid_long,bgrid_long);
    pdf1_n_b = reshape(pdf1_n_b,[n_networth,n_b]);
    pmf1_n_b = pdf1_n_b .* dnb;
    pmf1_n_b = pmf1_n_b ./ sum(pmf1_n_b(:),'omitnan');
    
    % marginal distr over n
    pmf1_n = sum(pmf1_n_b,2,'omitnan');
    
    % condl distr b|n
    pmf_b_condl_n = pmf1_n_b ./ pmf1_n;
    pmf_b_condl_n = permute(pmf_b_condl_n,[2 1]);

    % get mpcs over (b,a)
    m1 = stats1.mpcs(5).mpcs(:,1);
    pmf1 = stats1.pmf;

    P1ab = sum(reshape(pmf1,[],income0.ny*p1.nz),2);
    m1 = reshape(m1,[],income0.ny*p1.nz) .* reshape(pmf1,[],income0.ny*p1.nz);
    m1 = sum(m1,2) ./ P1ab;
    m1 = reshape(m1,[p1.nb_KFE p1.na_KFE]);
    pmf1 = reshape(P1ab,[p1.nb_KFE,p1.na_KFE]);
    
    %% --------------------------------------------------------------------
	% ONE ASSET CASE
	% ---------------------------------------------------------------------

    m0 = stats0.mpcs(5).mpcs(:,1);
    pmf0 = stats0.pmf;
    
    P0ab = sum(reshape(pmf0,[],income0.ny),2);
    m0 = reshape(m0,[],income0.ny*p0.nz) .* reshape(pmf0,[],income0.ny*p0.nz);
    m0 = sum(m0,2) ./ P0ab;
    m0 = reshape(m0,[p0.nb_KFE p0.na_KFE]);
    m0 = m0(:,1);
    P0ab = reshape(P0ab,[p0.nb_KFE p0.na_KFE]);
    pmf0 = P0ab(:,1);
    
    b_tilde = (grdKFE0.b.dF(:,1) + grdKFE0.b.dB(:,1)) / 2;
    b_tilde(1) = grdKFE0.b.dF(1,1) / 2;
    b_tilde(end) = grdKFE0.b.dB(end,1) / 2;
    
    pdf0_b = pmf0 ./ b_tilde;
    pdf0_interp = griddedInterpolant(grdKFE0.b.vec,pdf0_b,'linear','none');
    pdf0_n = pdf0_interp(ngrid);
    pmf0_n = pdf0_n .* dn_trapz;
    pmf0_n = pmf0_n / sum(pmf0_n,'omitnan');
    
    %% --------------------------------------------------------------------
	% INTERPOLATE MPCs ONTO NEW GRIDS
	% ---------------------------------------------------------------------
    
    % interpolate two-asset mpcs onto (n,a) grids
    grids = {grdKFE1.a.vec,grdKFE1.b.vec};
    m1_switched = permute(m1,[2 1]);
    mpc_interp_twoasset = griddedInterpolant(grids,m1_switched,'linear','none');
    
    mpc_twoasset = mpc_interp_twoasset(ngrid_long-bgrid_long,bgrid_long);
    mpc_twoasset = reshape(mpc_twoasset,[n_networth,n_b]);

    % two asset mpcs on net worth grid
    mpc_twoasset_nw = sum(mpc_twoasset .* pmf1_n_b,2,'omitnan');
    mpc_twoasset_nw = mpc_twoasset_nw ./ sum(pmf1_n_b,2,'omitnan');
    
    % interpolate one-asset mpcs onto n grid
    mpc_interp = griddedInterpolant(grdKFE0.b.vec,m0,'linear');
    mpc_oneasset =  mpc_interp(ngrid);

    %% --------------------------------------------------------------------
    % SETUP FOR DECOMPOSITION
    % ---------------------------------------------------------------------
    decomp = struct();

    m1g1 = sum(mpc_twoasset_nw .* pmf1_n,1,'omitnan');
    m0g0 = sum(mpc_oneasset .* pmf0_n,1,'omitnan');

    decomp.mean_mpc_twoasset = m1g1;
    decomp.mean_mpc_oneasset = m0g0;
    decomp.Em1_minus_Em0 = m1g1 - m0g0;

%     fprintf('m1g1 = %.5f\n',decomp.mean_mpc_twoasset)
%     fprintf('actual mean mpc = %.5f\n',stats1.mpcs(5).avg_0_t(1))
%     
%     fprintf('\n\nm0g0 = %.5f\n',decomp.mean_mpc_oneasset)
%     fprintf('actual mean mpc = %.5f\n',stats0.mpcs(5).avg_0_t(1))

    % important integrals
    m0g1 = sum(mpc_oneasset .* pmf1_n,1,'omitnan');
    m1g0 = sum(mpc_twoasset_nw .* pmf0_n,1,'omitnan');

    %% --------------------------------------------------------------------
    % EFFECT OF MPC FUNCTION
    % ---------------------------------------------------------------------
    decomp.term1_mpcfn = m1g0 - m0g0;
    
    % decomposition around epsilon = 0
    term1a_part1 = mpc_twoasset(1,1) * pmf_b_condl_n(1,1) * pmf0_n(1);
    term1a_part2 = mpc_oneasset(1) * pmf_b_condl_n(1,1) * pmf0_n(1);
    decomp.term1a(1) = sum(term1a_part1,1,'omitnan')...
        - sum(term1a_part2,1,'omitnan');
    term1b_part1 = mpc_twoasset(2:end,1).* pmf_b_condl_n(1,2:end)' .* pmf0_n(2:end);
    term1b_part2 = mpc_oneasset(2:end).* pmf_b_condl_n(1,2:end)' .* pmf0_n(2:end);
    decomp.term1b(1) = sum(term1b_part1,1,'omitnan') - sum(term1b_part2,1,'omitnan');

    decomp.term1c = 0;
    for in = 2:n_networth
        increment_part1 = sum(mpc_twoasset(in,2:end)' .* pmf_b_condl_n(2:end,in) * pmf0_n(in),1,'omitnan');
        increment_part2 = -sum(mpc_oneasset(in) .* pmf_b_condl_n(2:end,in) * pmf0_n(in),1,'omitnan');
        
        decomp.term1c(1) = decomp.term1c + increment_part1 + increment_part2;
    end
    
    % decomposition around epsilon > 0
    % first term (integral from n = 0 to n = n*)
    integral_0_nstar = zeros(n_networth,1);
    for in = 1:n_networth
        integrand = sum(mpc_twoasset(in,:)' .* pmf_b_condl_n(:,in),'omitnan');
        if in == 1
            integral_0_nstar(1) = integrand * pmf0_n(in);
        elseif isfinite(pmf0_n(in))
            integral_0_nstar(in) = integrand * pmf0_n(in)...
                + integral_0_nstar(in-1);
        end
    end
    
    integral_interp = griddedInterpolant(ngrid,integral_0_nstar,'linear','none');
    i_eps = 2;
    for epsilon = epsilons
        decomp.term1a(i_eps) = integral_interp(epsilon);
        i_eps = i_eps + 1;
    end
    
    % second term
    i_eps = 2;
    for epsilon = epsilons
        integral_over_n = 0;
        for in = 1:n_networth
            ibmax = find(bgrid<=epsilon,1,'last');
            integrand = sum(mpc_twoasset(in,1:ibmax)' .* pmf_b_condl_n(1:ibmax,in),'omitnan');
            
            if isfinite(integrand) && isfinite(pmf0_n(in))
                integral_over_n = integral_over_n + integrand * pmf0_n(in);
            end
        end
        
        decomp.term1b(i_eps) = integral_over_n - decomp.term1a(i_eps);
        decomp.term1c(i_eps) = decomp.term1_mpcfn - decomp.term1a(i_eps)...
            - decomp.term1b(i_eps);
        i_eps = i_eps + 1;
    end

    %% --------------------------------------------------------------------
    % EFFECT OF THE DISTRIBUTION OF NET WORTH
    % ---------------------------------------------------------------------
    decomp.term2_networth = m0g1 - m0g0;
    
    % decomposition of net worth term around epsilon = 0
    decomp.term2a(1) = mpc_oneasset(1) * (pmf1_n(1) - pmf0_n(1));
    term2b_part1 = sum(mpc_oneasset(2:end) .* pmf1_n(2:end),1,'omitnan');
    term2b_part2 = sum(mpc_oneasset(2:end) .* pmf0_n(2:end),1,'omitnan');
    decomp.term2b(1) = term2b_part1 - term2b_part2;
    
    % decomposition of net worth term around epsilon > 0
    integral_0_nstar = zeros(n_networth,1);
    for in = 1:n_networth
        % integral from 0 to ngrid(in) of m0(n) * [g1(n) - g0(n)]
        term1 = sum(mpc_oneasset(1:in) .* pmf1_n(1:in),1,'omitnan');
        term2 = - sum(mpc_oneasset(1:in) .* pmf0_n(1:in),1,'omitnan');
        integral_0_nstar(in) = term1 + term2;
    end
    integral_interp = griddedInterpolant(ngrid,integral_0_nstar,'linear','none');
    
    i_eps = 2;
    for epsilon = epsilons
        decomp.term2a(i_eps) = integral_interp(epsilon);
        decomp.term2b(i_eps) = decomp.term2_networth - decomp.term2a(i_eps);
        i_eps = i_eps + 1;
    end
    
    %% --------------------------------------------------------------------
    % INTERACTION
    % ---------------------------------------------------------------------
    mdiff = mpc_twoasset_nw - mpc_oneasset;
    mdiffg0 = sum(mdiff .* pmf0_n,1,'omitnan');
    mdiffg1 = sum(mdiff .* pmf1_n,1,'omitnan');
    decomp.term3_interaction = mdiffg1 - mdiffg0;

    decomp