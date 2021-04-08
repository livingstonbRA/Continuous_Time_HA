function [H_special,c_special,d_special] = deal_with_special_case(p,income,grids,r_b_mat,VaB)
	% this function computes the policy functions associated
	% with the special case that b = bmin, a > 0 and households
	% withdraw only enough to consume, so that a does not accumulate

	if strcmp(grids.gtype,'HJB')
		nb = p.nb;
		na = p.na;
	else
		nb = p.nb_KFE;
		na = p.na_KFE;
	end
	
	lambda = zeros(nb,na,p.nz,income.ny);
	rb = r_b_mat(1,1,1,1) * p.bmin;
	for ia = 2:p.na
	for iz = 1:p.nz
	for k = 1:p.ny
		l_fn = @(x) aux.lambda_function(x,VaB(1,ia,iz,k),rb,grids.a.vec(ia),income.y.vec(k),p,iz);

		options = optimset('TolX',1e-8,'TolFun',1e-11);
        [lambda(1,ia,iz,k),fval] = fminbnd(@(x) l_fn(x)^2,0,1e9,options);

        if fval > 1e-10
            error('no convergence for special case')
        end
	end
	end
	end

	chi1inv_arg = VaB ./ lambda - 1;
    d_special = aux.AdjustmentCost.derivative_inverse(chi1inv_arg,grids.a.matrix,p);
    d_special(~isfinite(d_special)) = 0;


    c_special = zeros(nb,na,p.nz,p.ny);
    if na == 2
    	riskaver_mat = p.riskaver;
    else
    	riskaver_mat = reshape(p.riskaver,[1 numel(p.riskaver)]);
    end
    c_special(1,2:na,:,:) = lambda(1,2:na,:,:) .^(-1./riskaver_mat);

    H_special = - 1e12 * ones(nb,na,p.nz,p.ny);
    H_special(1,2:na,:,:) = VaB(1,2:na,:,:) .* d_special(1,2:na,:,:);
end