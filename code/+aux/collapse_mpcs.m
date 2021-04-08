function [mpcs_x, pmf_x] = collapse_mpcs(mpcs, pmf)

	nb = size(mpcs, 1);
	na = size(mpcs, 2);

	one_asset = (na <= 2);
	if one_asset
		mpcs = reshape(mpcs, nb, []);
		pmf = reshape(pmf, nb, []);
		last_dim = 2;
	else
		mpcs = reshape(mpcs, nb, na, []);
		pmf = reshape(pmf, nb, na, []);
		last_dim = 3;
	end

	pmf_x = sum(pmf, last_dim);

	mpcs_x = size(pmf_x);
	psmall = pmf_x < 1e-9;
	mpcs_unw_means = mean(mpcs, last_dim);
	mpcs_x(psmall) = mpcs_unw_means(psmall);

	mpcs_w_means = sum(mpcs .* pmf ./ pmf_x, last_dim);
	mpcs_x(~psmall) = mpcs_w_means(~psmall);

	if one_asset
		mpcs_x = reshape(mpcs_x, nb, []);
	else
		mpcs_x = reshape(mpcs_x, nb, na, []);
	end
end