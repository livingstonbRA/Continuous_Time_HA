classdef Frisch

	methods (Static)
		function d = labor_disutility(h, weight, frisch)
			d = weight * (h .^ (1 + 1/frisch)) ./ (1 + 1/frisch);
		end

		function md = marginal_disutility(h, weight, frisch)
			md = weight * (h .^ (1/frisch));
		end

		function inv_md = inv_marginal_disutility(v, weight, frisch)
			inv_md = (v./weight) .^ frisch;
		end
	end
end