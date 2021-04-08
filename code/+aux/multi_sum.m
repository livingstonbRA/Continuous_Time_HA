function arr_out = multi_sum(arr, dims, flatten)
	if nargin < 3
		flatten = false;
	end

	arr_out = arr;
	for ddim = dims
		arr_out = sum(arr_out, ddim);
	end

	if flatten
		arr_out = arr_out(:);
	end
end