function CI = stoye_CI(lower_estim, upper_estim, varcov, signif_level, optimopts)

    % Stoye (2009) confidence interval for partially identified parameter

    % Inputs to routine
    Delta = upper_estim-lower_estim; % Point estimate of length of identif. set
    sigma_lower = sqrt(varcov(1,1)); % Std. dev. of lower bound estimate
    sigma_upper = sqrt(varcov(2,2)); % Std. dev. of upper bound estimate
    rho = varcov(1,2)/(sigma_lower*sigma_upper); % Correlation of lower and upper bound estimates
    
    % Numerically minimize CI length subject to coverage constraints
    c_opt = fmincon(@(c) [sigma_lower sigma_upper]*c', ...
                 [lower_estim-norminv(1-signif_level/2)*sigma_lower upper_estim+norminv(1-signif_level/2)*sigma_upper], ...
                 [], [], [], [], [], [], ...
                 @(c) deal(0, 1-signif_level - [stoye_bound(c, rho, Delta/sigma_upper); stoye_bound(fliplr(c), rho, Delta/sigma_lower)]), ...
                 optimopts);
    
	% Confidence interval
    CI = [lower_estim-sigma_lower*c_opt(1) upper_estim+sigma_upper*c_opt(2)];
                
end