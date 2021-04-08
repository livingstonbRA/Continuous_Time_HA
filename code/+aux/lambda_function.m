function fval = lambda_function(lambda,Va,rb,a_grid,income,p,iz)

    term1 = - rb;

	term2 = lambda .^(-1/p.riskaver(iz));

	chi1_inv_arg = Va ./ lambda - 1;
	term3 = AdjustmentCost.derivative_inverse(chi1_inv_arg,a_grid,p);

	term4 = AdjustmentCost.cost(term3,a_grid,p);

	term5 = - income;

	fval = term1 + term2 + term3 + term4 + term5;

end