function lag = selectlag_IC(data,maxlag,penalty);

% preliminaries
n_x = size(data,2);

% compute information criterion for each possible lag length
IC  = NaN(maxlag,1);
for i = 1:maxlag
    
    % set lag length
    p      = i;
    
    % estimate VAR
    data_est   = data(maxlag-p+1:end,:);
    T          = size(data_est,1);
    T_VAR      = T - p;
    X          = lagmatrix(data_est,1:p);
    X          = X(p+1:end,:);
    Y          = data_est(p+1:end,:);
    VAR_coeff  = [X ones(length(X),1)]\Y;
    VAR_res    = Y-[X ones(length(X),1)]*VAR_coeff;
%     Sigma_u    = (VAR_res'*VAR_res)/(T_VAR-n_x*p-1);
    Sigma_u    = (VAR_res'*VAR_res)/T_VAR;
    
    % compute information criterion
    T       = size(data,1)-maxlag;
%     IC(i)  = log(det(Sigma_u)) + i * n_x^2 * 2/T;
    IC(i)  = log(det(Sigma_u)) + i * n_x^2 * penalty(T);
    
end

lag = find(IC==min(IC));