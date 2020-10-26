function integr = stoye_bound(c, rho, Delta_over_sigma)

    % Integral used in Stoye (2009) confidence interval routine
    % See Appendix B of Stoye (2009)
    
    integr = integral(@(z) normcdf((rho*z + c(2) + Delta_over_sigma)/sqrt(1-rho^2)).*normpdf(z), ...
                      -10, ...
                      c(1));
end