classdef TransitionMatrixConstructor < handle
    % This class constructs the A matrix for the HJB or KFE
    % First, instantiate the class with the required arguments,
    % then call the construct method with the policy functions
    % and value function as arguments.

    properties (SetAccess = protected)
    	p;

    	income;

    	% A Grid object.
        grids;

        nb;
        na;
        ny;
        nz;
        n_states;

        % Shape of the full state space, a row vector.
        shape;

        % Function handle to auto-repmat to full state space dims.
        repmat_to_ss;

        % returns risk is on or off
        returns_risk;

        % Grid deltas
        asset_dB;
        asset_dF;
        asset_dSum; % (deltaB + deltaF)

        % Boolean mask for top of the b (1-asset) or a (2-asset) grid
        top;

        % Pre-computed term when returns are risky
        risk_term;

        % Offsets from the diagonal for risky returns.
        offsets_for_rr;

        % 'HJB' or 'KFE' grid
        gridtype;
    end

    methods
        function obj = TransitionMatrixConstructor(p, income, grids, returns_risk)
            % Class constructor.
            %
            % Parameters
            % ----------
            % p : An object with the required attributes outlined in the
            %	  Properties block.
            %
            % income : An object with the required attributes outlined in
            %	  the Properties block.
            %
            % grids : A Grid object.
            %
            % returns_risk : A boolean indicator for whether or not returns
            %	  risk is to be included in the transition matrix.

            import aux.repmat_auto

            obj.nb = numel(grids.b.vec);
            obj.na = numel(grids.a.vec);
            obj.nz = grids.nz;
            obj.ny = numel(income.y.vec);
            obj.n_states = obj.nb * obj.na * obj.nz * obj.ny;
            obj.shape = [obj.nb obj.na obj.nz obj.ny];
            obj.repmat_to_ss = @(arr) repmat_auto(arr, obj.shape);

            obj.gridtype = grids.gtype;
            obj.grids = grids;
            obj.p = p;
            obj.income = income;

            obj.returns_risk = returns_risk;

            if p.sigma_r > 0
                obj.perform_returns_risk_computations();
            end
        end

        %% --------------------------------------------------------------------
        % Construct the A Matrix
        % ---------------------------------------------------------------------
        function [A, stationary, drifts] = construct(obj, model, V)
            % Constructs the transition matrix.
            %
            % Parameters
            % ----------
            % model : A structure containing the policy functions, saving (s)
            %   and deposits (d) are used here.
            %
            % V : The value function, of shape (nb, na, nz, ny).
            %
            % Returns
            % -------
            % A : The sparse transition matrix of shape (nb*na*nz*ny, nb*na*nz*ny).
            %
            % stationary : A boolean mask indicating states where drift in the
            %   risky asset was neither backward nor forward. Risk computations
            %   for these states will be included outside the A matrix.
            
            drifts = obj.compute_asset_drifts(model);

            A = obj.compute_liquid_transitions(drifts.b_B, drifts.b_F);
            A = A + obj.compute_illiquid_transitions(drifts.a_B, drifts.a_F);

            stationary = [];
            if obj.returns_risk
            	A = A + obj.compute_Vaa_terms();
            end  

            if obj.p.SDU && obj.returns_risk
            	[Arisk_Va, stationary] = obj.compute_Va_terms(drifts, V);
            	A = A + Arisk_Va;
            end
        end
    end

    methods (Access=protected)
    	%% --------------------------------------------------------------------
        % Compute Asset Drifts
        % ---------------------------------------------------------------------
        function drifts = compute_asset_drifts(obj, model)
        	% Approximates the drifts of the liquid and illiquid assets.
        	% The input variable 'model' must contain the policy
        	% functions 's' and 'd'.

            if strcmp(obj.gridtype,'KFE')
            	adrift = model.d + obj.income.nety_KFE_illiq_hourly(model.h);
                drifts.a_B = min(adrift, 0);
                drifts.a_F = max(adrift, 0);

                drifts.b_B = min(model.s_c + model.s_d, 0);
                drifts.b_F = max(model.s_c + model.s_d, 0);
            elseif strcmp(obj.gridtype, 'HJB')
            	adrift = obj.income.nety_HJB_illiq_hourly(model.h);
                drifts.a_B = min(model.d, 0) + min(adrift, 0);
                drifts.a_F = max(model.d, 0) + max(adrift, 0);

                drifts.b_B = min(model.s_d, 0) + min(model.s_c, 0);
                drifts.b_F = max(model.s_d, 0) + max(model.s_c, 0);    
            end
        end

    	%% --------------------------------------------------------------------
        % Liquid Asset Transitions
        % ---------------------------------------------------------------------
    	function A_liquid = compute_liquid_transitions(obj, bdriftB, bdriftF)
    		% Computes the sparse, square matrix for liquid asset transitions,
    		% adjusted for grid deltas.

    		lowdiag = - bdriftB ./ obj.grids.b.dB; 
            lowdiag(1,:,:,:) = 0;

            centerdiag = bdriftB ./ obj.grids.b.dB - bdriftF ./ obj.grids.b.dF; 
            centerdiag(1,:,:,:) = - bdriftF(1,:,:,:) ./ obj.grids.b.dF(1,:); 
            centerdiag(obj.nb,:,:,:) = bdriftB(obj.nb,:,:,:) ./ obj.grids.b.dB(obj.nb,:);

            updiag = bdriftF ./ obj.grids.b.dF; 
            updiag(obj.nb,:,:,:) = 0;

            A_liquid = aux.sparse_diags(...
                [lowdiag(:), centerdiag(:), updiag(:)],...
                [-1, 0, 1]);
        end

    	%% --------------------------------------------------------------------
        % Illiquid Asset Transitions
        % ---------------------------------------------------------------------
        function A_illiquid = compute_illiquid_transitions(obj, adriftB, adriftF)
        	% Computes the sparse, square matrix for illiquid asset transitions,
    		% adjusted for grid deltas.

    		lowdiag = -adriftB ./ obj.grids.a.dB; 
            lowdiag(:,1,:,:) = 0;

            centerdiag = adriftB ./ obj.grids.a.dB - adriftF ./ obj.grids.a.dF; 
            centerdiag(:,1,:,:) = -adriftF(:,1,:,:) ./ obj.grids.a.dF(:,1); 
            centerdiag(:,obj.na,:,:) = adriftB(:,obj.na,:,:) ./ obj.grids.a.dB(:,obj.na);

            updiag = adriftF ./ obj.grids.a.dF; 
            updiag(:,obj.na,:,:) = zeros(obj.nb,1,obj.nz,obj.ny);

            A_illiquid = aux.sparse_diags(...
                [lowdiag(:), centerdiag(:), updiag(:)],...
                [-obj.nb, 0, obj.nb]);
        end

        %% --------------------------------------------------------------------
        % Computations for Returns Risk
        % ---------------------------------------------------------------------
        function perform_returns_risk_computations(obj)
            % Pre-computes terms needed for risky returns.
            %
            % Modified
            % --------
            % obj.asset_dB : deltaBackward for the grid of the risky asset
            %
            % obj.asset_dF : deltaForward for the grid of the risky asset
            %
            % obj.top : indicator for top of risky asset grid
            %
            % obj.risk_term : (asset * sigma_r) ^ 2
            %
            % obj.offsets_for_rr : vector of offsets from the diagonal
            %   for the risky asset
            %
            % obj.asset_dSum : deltaBackward + deltaForward

            obj.top = false(obj.nb, obj.na, obj.nz, obj.ny);
            if obj.p.OneAsset == 1
                obj.asset_dB = obj.repmat_to_ss(obj.grids.b.dB);
                obj.asset_dF = obj.repmat_to_ss(obj.grids.b.dF);

                obj.top(obj.nb, :, :, :) = true;
                obj.risk_term = (obj.grids.b.vec * obj.p.sigma_r) .^ 2;
                obj.offsets_for_rr = [-1, 0, 1];
            else
                obj.asset_dB = obj.repmat_to_ss(obj.grids.a.dB);
                obj.asset_dF = obj.repmat_to_ss(obj.grids.a.dF);

                obj.top(:, obj.na, :, :) = true;
                obj.risk_term = (obj.grids.a.wide * obj.p.sigma_r) .^ 2;
                obj.offsets_for_rr = [-obj.nb, 0, obj.nb];
            end

            obj.asset_dF(obj.top) = obj.asset_dB(obj.top);
            obj.asset_dSum = obj.asset_dB + obj.asset_dF;
        end

        function Arisk_Vaa = compute_Vaa_terms(obj)
            % This function computes the (1/2) * (a * sigma_r) ^2 * Vaa component
            % of the A matrix.
            %
            % Returns
            % -------
            % Arisk_Vaa : the component of the A matrix s.t. Arisk_Vaa * V produces
            %   a second-order finite difference of V, after multiplying by the
            %   (1/2) * (a * sigma_r) ^ 2

            % (1/2) * (a * sigma_r) ^ 2 * V_{i+1} * (dF*(dB+dF)/2)
            V_i_plus_1_term = obj.risk_term ./ (obj.asset_dF .* obj.asset_dSum);

            % - (1/2) * (a * sigma_r) ^ 2 * V_i * (1/dB + 1/dF) / ((dB+dF)/2)
            V_i_term = - obj.risk_term .* (1./obj.asset_dB + 1./obj.asset_dF) ./ obj.asset_dSum;
            
            % (1/2) * (a * sigma_r) ^ 2 * V_{i-1} * (dB*(dB+dF)/2)
            V_i_minus_1_term = obj.risk_term ./ (obj.asset_dB .* obj.asset_dSum);

            % impose V' = 0 at the top of the grid
            V_i_term(obj.top) = - obj.risk_term(obj.top) ./ (obj.asset_dB(obj.top) .* obj.asset_dSum(obj.top));
            V_i_plus_1_term(obj.top) = 0;

            arr = [V_i_minus_1_term(:), V_i_term(:), V_i_plus_1_term(:)];
            Arisk_Vaa = aux.sparse_diags(arr, obj.offsets_for_rr);
        end

        function [Arisk_Va, stationary] = compute_Va_terms(obj, drifts, V)
            % Computes the (1/2) * (a * sigma_r) ^ 2 * Va ^2 ... term when
            % returns are risky. Only applies to SDU case.
            %
            % Parameters
            % ----------
            %
            % V : the value function, of shape (nb, nz, nz, ny)
            %
            % driftB : the backward drift of the risky asset, where negative
            %
            % driftF : the forward drift of the risky asset, where positive
            %
            % Returns
            % -------
            % Arisk_Va : after multiplying by V, this is the component of the A matrix that
            %   produces the Va ^ 2 term
            %
            % stationary : a boolean mask indicating states where drift in the
            %   risky asset was neither backward nor forward. risk computations
            %   for these states will be included outside the A matrix

            assert(obj.p.SDU, "This term only applies to SDU")

            import computation.fd_firstorder
            if obj.p.OneAsset
            	V1 = fd_firstorder(V, obj.grids.b.dB, obj.grids.b.dF, 1);
            	driftB = drifts.b_B;
            	driftF = drifts.b_F;
            else
            	V1 = fd_firstorder(V, obj.grids.a.dB, obj.grids.a.dF, 2);
            	driftB = drifts.a_B;
            	driftF = drifts.a_F;
            end

            % [V1B, V1F] = obj.first_diffs(V);
            V1 = (driftB < 0) .* V1.B + (driftF > 0) .* V1.F;

            % Now add zeta(a) * Va^2 term.
            if obj.p.invies == 1
                SDU_adj = obj.risk_term .* V1 * (1 - obj.p.riskaver);
            else
                SDU_adj = obj.risk_term .* V1 ./...
                	V * (obj.p.invies - obj.p.riskaver) / (1 - obj.p.invies);
            end

            updiag = zeros(obj.shape);
            lowdiag = zeros(obj.shape);

            lowdiag(driftB < 0) = - SDU_adj(driftB < 0) ./ obj.asset_dB(driftB < 0);
            updiag(driftF > 0) = SDU_adj(driftF > 0) ./ obj.asset_dF(driftF > 0);

            if obj.p.OneAsset
                lowdiag(1, :, :, :) = 0;
                updiag(obj.nb, :, :, :) = 0;
            else
                lowdiag(:, 1, :, :) = 0;
                updiag(:, obj.na, :, :) = 0;
            end

            centdiag = - lowdiag - updiag;
            Arisk_Va = aux.sparse_diags(...
                [lowdiag(:), centdiag(:), updiag(:)],...
                obj.offsets_for_rr);

            % Some states may have no drift, we will need to deal with them
            % separately by adding directly to the right_hand_side of the HJB.
            stationary = (driftB >= 0) & (driftF <= 0);
        end
    end
end

