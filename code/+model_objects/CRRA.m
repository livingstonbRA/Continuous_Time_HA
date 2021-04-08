classdef CRRA
	% Provides functions for CRRA utility, with or
	% without IES heterogeneity.

	methods (Static)
		function u = utility(c, invies)
			u_log = log(c);
			u_other = c .^ (1 - invies) ./ (1 - invies);

			u_log(~isfinite(u_log)) = 0;
			u_other(~isfinite(u_other)) = 0;

			mask = (invies == 1);
			u = mask .* u_log + (~mask) .* u_other;
		end

		function muc = marginal_utility(c, invies)
			muc = c .^ (-invies); 
		end

		function c = u1inv(v, invies, zdim)
			c = v .^ (-1 ./ invies);
		end
	end
end