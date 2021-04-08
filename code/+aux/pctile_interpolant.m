function interpolant = pctile_interpolant(values, pmf)
	sorted_mat = sortrows([values(:), pmf(:)]);
	[cdf_u, iu] = unique(cumsum(sorted_mat(:,2)), 'last');

	if numel(iu) > 2
		values_u = sorted_mat(iu,1);
		interpolant = griddedInterpolant(cdf_u, values_u, 'pchip', 'nearest');
	else
		interpolant = @(x) NaN;
	end
end