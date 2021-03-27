function VAR_IV_boot = bootstrapVAR_IV(VAR,model,data,settings)

% Homoskedastic recursive residual VAR bootstrap, with external IV

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_y    = model.n_y;
n_x    = model.n_x;
T      = settings.T;
T_VAR  = VAR.T_VAR;
p      = VAR.laglength;
n_boot = settings.n_boot;

VAR_res_y   = VAR.VAR_res_y;
IV_res      = VAR.IV_res;
IV_beta     = VAR.IV_beta;
VAR_coeff_y = VAR.VAR_coeff_y;

%----------------------------------------------------------------
% VAR i.i.d. residual bootstrap
%----------------------------------------------------------------

VAR_coeff_y_boot = NaN(p*n_y+1,n_y,n_boot);
Sigma_u_y_boot = NaN(n_y,n_y,n_boot);
impact_irs_boot = NaN(n_y,n_boot);

% disp('Bootstrapping VAR...');

for i_boot = 1:n_boot
    res_inds = ceil(size(VAR_res_y,1)*rand(T_VAR,1)); % Random indices of bootstrap residuals
    u_boot = VAR_res_y(res_inds,:); % reshuffle VAR residuals
    v_boot = IV_res(res_inds); % Reshuffle IV residual
    
    % create new artificial data
    Y_boot = NaN(T,n_y);
    Z_boot = NaN(T,1);
    
    % Random block of initial observations
    init_inds = randi(T_VAR)+(1:p);
    Y_boot(1:p,:) = data.y(init_inds,:);
    Z_boot(1:p) = data.z(init_inds);
    Xlag = lagmatrix(Y_boot,1:p);
    Xlag = Xlag(p+1,:);
    Xlag_IV = lagmatrix([Y_boot Z_boot],1:p);
    Xlag_IV = Xlag_IV(p+1,:);
    
    for t = 1:T_VAR
        Y_boot(p+t,:) = [Xlag,1] * VAR_coeff_y + u_boot(t,:);
        Z_boot(p+t) = [Xlag_IV,1] * IV_beta + v_boot(t,:);
        Xlag = [Y_boot(p+t,:),Xlag(1:end-n_y)];
        Xlag_IV = [Y_boot(p+t,:) Z_boot(p+t) Xlag_IV(1:end-n_y-1)];
    end
    
    settings.select_VAR_simlaglength = 0;
    settings.VAR_simlaglength = p;
    VAR_boot = estimateVAR_IV(Y_boot,Z_boot,settings);
    VAR_coeff_y_boot(:,:,i_boot) = VAR_boot.VAR_coeff_y;
    Sigma_u_y_boot(:,:,i_boot) = VAR_boot.Sigma_u_y;
    impact_irs_boot(:,i_boot) = VAR_boot.impact_irs;
    
end

% collect results
VAR_IV_boot.VAR_coeff_y = VAR_coeff_y_boot;
VAR_IV_boot.Sigma_u_y   = Sigma_u_y_boot;
VAR_IV_boot.impact_irs  = impact_irs_boot;

end