function [edges, counts, hist_out] = smoothed_histogram(...
	agrid, pmf, nbins, amax, rescale_dist)
	if nargin < 5
		rescale_dist = false;
	end

	if rescale_dist
		too_large = agrid > amax;
		agrid = agrid(~too_large);
		pmf = pmf(~too_large);
		pmf = pmf / sum(pmf(:));
	end

	a_cdf = cumsum(pmf);
	[agrid_u, iu] = unique(agrid, 'last');
	a_cdf = a_cdf(iu);

	cdf_interp = griddedInterpolant(agrid_u, a_cdf, 'pchip', 'nearest');

	amin = agrid(1);
	spacing = (amax - amin) / nbins;
	edges = 1:nbins;
	edges = amin + edges * spacing;

	counts = zeros(nbins, 1);
	for ibin = 1:nbins
		bin_start = edges(ibin) - spacing;
		bin_end = edges(ibin);

		P_lt_start = cdf_interp(bin_start);

		if ibin < nbins
			P_lt_end = cdf_interp(bin_end);
		else
			P_lt_end = 1;
		end

		counts(ibin) = P_lt_end - P_lt_start;
	end

	edges = [amin, edges];

	if nargout > 2
		fig = figure();
		hist_out = histogram(...
			'Parent', fig, 'BinEdges', edges,...
			'BinCounts', counts);
		ylabel("Probability density")
	end
end