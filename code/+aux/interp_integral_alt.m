function interpolant = interp_integral_alt(grids, integrand_values, pmf)
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
	%	is the approximated value of the integral from 0 to x.

	nassets = numel(grids);
	weighted = integrand_values .* pmf;

	integral_values = zeros(size(pmf));
	nb = size(pmf, 1);
	na = size(pmf, 2);

	bvec = reshape(1:nb, [], 1);
	avec = reshape(1:na, 1, []);

	for ii = 1:nb
		for jj = 1:na
			mask = (bvec <= ii) & (avec <= jj);
			integral_values(ii,jj) = sum(weighted(mask));
        end
    end

	interpolant = griddedInterpolant(...
		grids, integral_values, 'spline', 'nearest');
end
