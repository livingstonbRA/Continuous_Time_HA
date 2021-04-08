function arr_out = repmat_auto(arr, shape_out)
	shape_in = size(arr);
	ndims_out = numel(shape_out);
	ndims_in = numel(shape_in);

	shape_out_matched = shape_out(1:ndims_in);
	incorrect_dims = (shape_out_matched ~= shape_in);

	if any(shape_in(incorrect_dims) ~= 1)
		error("Input shape and output shape do not match")
	end

	rep_vec = ones(1, ndims_in);
	rep_vec(incorrect_dims) = shape_out_matched(incorrect_dims);

	if ndims_in < ndims_out
		shape_out_unmatched = shape_out(ndims_in+1:end);
		rep_vec = [rep_vec, shape_out_unmatched];
	end

	arr_out = repmat(arr, rep_vec);
end