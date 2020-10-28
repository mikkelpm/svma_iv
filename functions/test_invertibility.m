function [wald_stat,df,pval] = test_invertibility(VAR)

    % Granger causality pre-test for invertibility
    
    % Dimensions
    [k,n_x] = size(VAR.VAR_coeff);
    p = (k-1)/n_x;
    
    % Var-cov matrix of vec(VAR_coeff)
    varcov = kron(VAR.Sigma_u,inv(VAR.X'*VAR.X));
    
    % Selection matrix picking out lags of z in each y equation
    the_eye = eye(n_x);
    sel = [repmat(the_eye(:,n_x),p,1); 0]==1;
    
    % Test in each y equation separately
    wald_stat.eqns = nan(1,n_x-1);
    df.eqns = p*ones(1,n_x-1);
    pval.eqns = nan(1,n_x-1);
    for i=1:n_x-1
        the_coef = VAR.VAR_coeff(:,i);
        the_varcov = varcov((i-1)*k+1:i*k,(i-1)*k+1:i*k);
        wald_stat.eqns(i) = the_coef(sel)'*(the_varcov(sel,sel)\the_coef(sel));
        pval.eqns(i) = 1-chi2cdf(wald_stat.eqns(i), df.eqns(i));
    end
    
    % Overall test for all y equations jointly
    sel2 = [repmat(sel,n_x-1,1); zeros(k,1)]==1;
    coef_vec = VAR.VAR_coeff(:);
    wald_stat.all = coef_vec(sel2)'*(varcov(sel2,sel2)\coef_vec(sel2));
    df.all = p*(n_x-1);
    pval.all = 1-chi2cdf(wald_stat.all, df.all);

end