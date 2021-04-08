function deriv = lambda_function_deriv(lambda,Va,a_grid,p)
    chi1inv_arg = Va ./ lambda - 1;
	
    term1 = - (1/p.riskaver) * lambda .^(-1/p.riskaver-1);
    
    term2a = max(a_grid,p.a_lb) .* p.chi1/p.chi2...
        .* (-p.chi0 - abs(chi1inv_arg)) .^(1/p.chi2-1);
    term2 = - Va .* lambda .^(-2) .* term2a;
    
    d_temp = AdjustmentCost.derivative_inverse(chi1inv_arg,a_grid,p);
    chi1_chi1inv_neg = - p.chi0 - (-d_temp./(a_grid*p.chi1)).^(p.chi2);
    chi1_chi1inv_pos = p.chi0 + (d_temp./(a_grid*p.chi1)).^(p.chi2);
    chi1_chi1inv = (d_temp<0) .* chi1_chi1inv_neg + (d_temp>0) .* chi1_chi1inv_pos;
    term3 = chi1_chi1inv .* term2;
    
    deriv = term1 + term2 + term3;
    

end