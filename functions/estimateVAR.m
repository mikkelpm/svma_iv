function VAR_sim = estimateVAR(data,settings)

% Least-squares VAR estimation

% preliminaries
n_x = size(data,2);
T   = size(data,1);

% select lag length
if settings.select_VAR_simlaglength == 1
    p = selectlag_IC(data,settings.max_simlaglength,settings.penalty);
elseif settings.select_VAR_simlaglength == 0
    p = settings.VAR_simlaglength;
end
T_VAR  = T - p;

% estimate VAR
X          = lagmatrix(data,1:p);
X          = X(p+1:end,:);
X_expand   = [X ones(length(X),1)];
Y          = data(p+1:end,:);
VAR_coeff  = X_expand\Y;
VAR_res    = Y-X_expand*VAR_coeff;
% Sigma_u    = (VAR_res'*VAR_res)/(T_VAR-n_x*p-1);
Sigma_u    = (VAR_res'*VAR_res)/T_VAR;

% collect results
VAR_sim.T_VAR     = T_VAR;
VAR_sim.VAR_coeff = VAR_coeff;
VAR_sim.Sigma_u   = Sigma_u;
VAR_sim.VAR_res   = VAR_res;
VAR_sim.X         = X_expand;
VAR_sim.laglength = p;

end