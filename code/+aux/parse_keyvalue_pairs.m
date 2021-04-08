function options = parse_keyvalue_pairs(defaults, varargin)
	% This function takes key-value arguments and returns
	% a structure. If a field in the 'defaults' structure
	% wasn't passed as an argument, it's value is set from
	% 'defaults'. If a key is passed that is not a field
	% of 'defaults', an error is thrown by the parser.

	parser = inputParser;
    if isstruct(defaults)
        attributes = fieldnames(defaults)';
    else
        attributes = properties(defaults)';
    end
	for field = attributes
		addParameter(parser, char(field), defaults.(char(field)));
	end
	parse(parser, varargin{:});

	options = parser.Results;
end