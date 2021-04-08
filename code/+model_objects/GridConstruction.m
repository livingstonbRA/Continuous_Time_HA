classdef GridConstruction
	methods (Static)
		function vgrid = create_curved_grid(vmin, vmax, curv, npts, varargin)
			parser = inputParser;
			addOptional(parser, 'term1_curv', 1);
			addOptional(parser, 'term1_wt', 0);
			addOptional(parser, 'flip', false);
			parse(parser, varargin{:});

			term1_curv = parser.Results.term1_curv;
			term1_wt = parser.Results.term1_wt;
			flip = parser.Results.flip;

			if flip
				vgrid = linspace(1, 0, npts)';
			else
				vgrid = linspace(0, 1, npts)';
			end

			vgrid = term1_wt * vgrid .^ (1 / term1_curv) ...
				+ vgrid .^ (1 / curv);

			if flip
				vgrid = 1 - vgrid;
			end
				
			vgrid = vgrid / max(abs(vgrid));
			vgrid = vmin + (vmax - vmin) * vgrid;
		end

		function gridpos = create_positive_asset_grid(vmax, curv, npts, varargin)
			gridpos = model_objects.GridConstruction.create_curved_grid(0, vmax, curv, npts, varargin{:});
		end

		function gridneg = create_negative_asset_grid(vmin, curv, npts, varargin)
			npts_neg1 = ceil(npts / 2) + 1;
			npts_neg2 = npts - npts_neg1 + 2;
			mid_neg = vmin / 2;

			% part of grid close to borrowing limit
			gridneg1 = model_objects.GridConstruction.create_curved_grid(vmin, mid_neg, curv, npts_neg1, varargin{:});

			% part of grid close to zero
			gridneg2 = model_objects.GridConstruction.create_curved_grid(mid_neg, 0, curv, npts_neg2, varargin{:}, 'flip', true);
		   	gridneg2(1) = []; % remove midpoint
		    gridneg2(end) = []; % remove 0

	        gridneg = [gridneg1; gridneg2];
		end
	end
end