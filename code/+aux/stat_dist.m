function stat_dist = stat_dist(lambda)
	% Computes the stationary distribution for
	% the Markov transition matrix, lambda.

	validate_input(lambda);

	% Get the eigenvalues
	eigvl = eig(lambda); 

	% Get the eigenvectors
	[eigvc,~] = eig(lambda);
	[~,pos] = min(abs(eigvl));

	stat_dist = eigvc(:,pos); 
	stat_dist = stat_dist./sum(stat_dist);
end

function validate_input(lambda)
	assert(ismatrix(lambda),...
		"Input must be a square matrix");
	assert(size(lambda, 1) == size(lambda, 2),...
		"Input must be a square matrix");
end