function gini = direct_gini(level, distr)
	% computes the approximate Gini index
	
	% Parameters
	% ----------
	%	level : levels of the desired variable
	%		    (e.g. wealth)
	%
	%	distr : probability mass associated with levels
	%
	% Returns
	% -------
	%	gini : the Gini index, a scalar

    % sort distribution and levels by levels
    sort1 = sortrows([level(:),distr(:)]);
    level_sort = sort1(:,1);
    dist_sort  = sort1(:,2);
    S = [0;cumsum(dist_sort .* level_sort)];
    gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
end