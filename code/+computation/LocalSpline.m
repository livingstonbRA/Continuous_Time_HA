classdef LocalSpline < handle
	properties
		sp;
		yy;
		x;

		lb_v;
		ub_v;
	end

	methods
		function obj = LocalSpline()
		end

		function yy = fit_spline(obj, v, cdf_v, x)
			% First find unsmoothed percentile of query value
			cdf_interp = griddedInterpolant(v, cdf_v,...
				'pchip', 'nearest');
			pct_x = cdf_interp(x);

			if pct_x == 0
				yy = v(1);
				return
			elseif pct_x >= 1
				yy = v(end);
				return
			end

			% Set region for interpolation
			curv = @(x) (x + x .^ 2 + x .^ 3) / 3;
			lb = curv(pct_x);
			ub = 1 - curv(1 - pct_x);

			[cdf_v_u, iu] = unique(cdf_v);
			icdf_interp = griddedInterpolant(cdf_v_u, v(iu),...
				'pchip', 'nearest');
			obj.lb_v = icdf_interp(lb);
			obj.ub_v = icdf_interp(ub);

			v_use = (v >= obj.lb_v) & (v <= obj.ub_v);

			obj.sp = polyfit(v(v_use), cdf_v(v_use), 3);
			obj.yy = polyval(obj.sp, x);
			obj.x = x;

			yy = obj.yy;
		end

		function show_fit(obj, v, cdf_v)
			plot(v, cdf_v);

			hold on

			fplot(@(x) polyval(obj.sp, x), [obj.lb_v, obj.ub_v]);

			plot(obj.x, obj.yy, '-o')

			legend('Raw data', 'Fitted spline', 'Chosen point')
			xlim([obj.lb_v, obj.ub_v])
			ylim('auto')
		end
	end
end

