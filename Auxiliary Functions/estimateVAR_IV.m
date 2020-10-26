function VAR_sim = estimateVAR_IV(data_y,data_z,settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_y = size(data_y,2);
T   = size(data_y,1);

%----------------------------------------------------------------
% VAR in Macro Aggregates
%----------------------------------------------------------------

% select lag length
if settings.select_VAR_simlaglength == 1
    p = selectlag_IC(data_y,settings.max_simlaglength,settings.penalty);
elseif settings.select_VAR_simlaglength == 0
    p = settings.VAR_simlaglength;
end
T_VAR  = T - p;

% estimate VAR
X          = lagmatrix(data_y,1:p);
X          = X(p+1:end,:);
Y          = data_y(p+1:end,:);
VAR_coeff  = [X ones(length(X),1)]\Y;
VAR_res    = Y-[X ones(length(X),1)]*VAR_coeff;
Sigma_u    = (VAR_res'*VAR_res)/(T_VAR-n_y*p-1);

% collect results
VAR_sim.T_VAR       = T_VAR;
VAR_sim.VAR_coeff_y = VAR_coeff;
VAR_sim.Sigma_u_y   = Sigma_u;
VAR_sim.VAR_res_y   = VAR_res;
VAR_sim.laglength   = p;

%----------------------------------------------------------------
% Instrument Covariance
%----------------------------------------------------------------

Xlag = lagmatrix([data_y data_z], 1:p);
Xlag = [Xlag(p+1:end,:) ones(VAR_sim.T_VAR,1)];
VAR_sim.IV_beta = Xlag\data_z(p+1:end);
VAR_sim.IV_res = data_z(p+1:end)-Xlag*VAR_sim.IV_beta; % Residuals

VAR_sim.cov_uz = VAR_sim.VAR_res_y'*VAR_sim.IV_res/VAR_sim.T_VAR; % Covariance of VAR reduced-form residuals and IV residual
VAR_sim.gamma = VAR_sim.cov_uz/sqrt(VAR_sim.cov_uz'*(VAR_sim.Sigma_u_y\VAR_sim.cov_uz)); % Impact impulse responses to identified shock

end