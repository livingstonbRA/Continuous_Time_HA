function M = sparse_diags(diags, offsets)
	% A wrapper for spdiags() which adds zeros to the
	% beginning/end of the input vectors depending on
	% the offsets, prior to constructing the matrix.
	% The size of the output matrix is inferred from the
	% number of rows in 'diags'.

	validate_inputs(diags, offsets);

	adj = diags;
	offsets = reshape(offsets, 1, []);

	for col = 1:numel(offsets)
		k = offsets(col);
		if k < 0
			adj(:, col) = [diags(-k+1:end, col); zeros(-k, 1)];
		elseif k > 0
			adj(:, col) = [zeros(k, 1); diags(1:end-k, col)];
		end
	end

	n = size(diags, 1);
	M = spdiags(adj, offsets, n, n);
end

function validate_inputs(diags, offsets)
	valid = ismatrix(diags) && ismatrix(offsets);
	if ~valid
		error('HACTLib:sparse_diags:NotMatrix',...
			'At least one of the inputs is not a matrix')
	elseif size(diags, 2) ~= numel(offsets)
		error('HACTLib:sparse_diags:IncompatibleShapes',...
			'Number of diagonals not equal to number of offsets')
	elseif max(abs(offsets)) > numel(diags)
		error('HACTLib:sparse_diags:IncompatibleShapes',...
			'One or more offsets is larger than the diagonal')
	elseif numel(offsets) ~= numel(unique(offsets))
		error('HACTLib:sparse_diags:Duplicates',...
			'Diagonal offsets must be unique')
	end
end