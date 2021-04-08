
function [d, Ic_special] = upwind_deposits(Vb, Va, adjcost, opt_d)
    % Performs upwinding for the deposit rate.
    %
    % Inputs
    % ------
    % Vb : A structure containing the approximate derivatives
    %   of V w.r.t. the liquid asset, from backward/forward
    %   first differences.
    %
    % Va : A structure containing the approximate derivatives
    %   of V w.r.t. the illiquid asset, from backward/forward
    %   first differences.
    %
    % adjcost : A function handle s.t. adjcost(d) returns
    %   an array of adjustment costs of shape size(d).
    %
    % opt_d : A function handle s.t. opt_d(Vb, Va) returns
    %   an array of the optimal deposit rate of shape
    %   size(d).
    %
    % Returns
    % -------
    % d : The deposit rate, of shape size(Vb).

    nb = size(Vb.B, 1);
	na = size(Vb.B, 2);

	% Forward difference in a, backward in b
    dFB = opt_d(Vb.B, Va.F);
    dFB(:,na,:,:) = 0;
    dFB(1,:,:,:) = 0;
    HdFB = Va.F .* dFB - Vb.B .* (dFB + adjcost(dFB));
    HdFB(:,na,:,:) = -1.0e12;
    HdFB(1,:,:,:) = -1.0e12;
    validFB = (dFB > 0) & (HdFB > 0);

    % Backward difference in a, forward in b
    dBF = opt_d(Vb.F, Va.B);
    dBF(:,1,:,:) = 0;
    dBF(nb,:,:,:) = 0;
    HdBF = Va.B .* dBF - Vb.F .* (dBF + adjcost(dBF));
    HdBF(:,1,:,:) = -1.0e12;
    HdBF(nb,:,:,:) = -1.0e12;
    validBF = (dBF <= -adjcost(dBF)) & (HdBF > 0);

    % Backward differences in a and b
    dBB = opt_d(Vb.B_adj, Va.B);
    dBB(:,1,:,:) = 0;
    HdBB = Va.B .* dBB - Vb.B_adj .* (dBB + adjcost(dBB));
    HdBB(:,1,:,:) = -1.0e12;
    validBB = (dBB > - adjcost(dBB)) & (dBB <= 0) & (HdBB > 0);

    % TODO: Code special case
    Ic_special = false(size(Vb.B));
    d_special = zeros(size(Vb.B));

    % Choice indicators
    IcFB = validFB & (~validBF | (HdFB >= HdBF)) ...
        & (~validBB | (HdFB >= HdBB)) & (~Ic_special);
    IcBF = validBF & (~validFB | (HdBF >= HdFB)) ...
        & (~validBB | (HdBF >= HdBB)) & (~Ic_special);
    IcBB = (~validFB | (HdBB >= HdFB)) & (~validBF | (HdBB >= HdBF)) ...
        & validBB & (~Ic_special);
    Ic00  = (~validFB) & (~validBF) & (~validBB) & (~Ic_special);

    Isum = Ic_special + IcFB + IcBF + IcBB + Ic00;
    if ~isequal(Isum, ones(size(Isum)))
        error('HACT:upwind_deposits:LogicError',...
            'Logicals do not sum to unity')
    end

    d   = IcFB .* dFB + IcBF .* dBF + IcBB .* dBB...
            + Ic_special .* d_special;
end