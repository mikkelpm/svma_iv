function VAR_boot = bootstrapVAR(VAR,model,data,settings)

%----------------------------------------------------------------
% Estimate VAR
%----------------------------------------------------------------

% preliminaries
n_x    = model.n_x;
data   = data.x;
T      = settings.T;
p      = VAR.laglength;
n_boot = settings.n_boot;

% sample size
T_VAR = T - p;

% estimate VAR
X          = lagmatrix(data,1:p);
X          = X(p+1:end,:);
Y          = data(p+1:end,:);
VAR_coeff  = [X ones(length(X),1)]\Y;
VAR_res    = Y-[X ones(length(X),1)]*VAR_coeff;

%----------------------------------------------------------------
% Bootstrap
%----------------------------------------------------------------

VAR_coeff_boot = NaN(p*n_x+1,n_x,n_boot);
Sigma_u_boot = NaN(n_x,n_x,n_boot);
data_start = data(1:p,:);

for i = 1:n_boot
    if T >= 1000
        if mod(i,10) == 0
            disp(['I am at iteration ' num2str(i)])
        end
    end
    u_boot = VAR_res(ceil(size(VAR_res,1)*rand(T_VAR,1)),:); % reshuffle residuals
    
    % create new artificial data
    data_boot = NaN(T,n_x);
    data_boot(1:p,:) = data_start;
    Xlag = X(1,:);
    for t = 1:T_VAR
        data_boot(p+t,:) = [Xlag,1] * VAR_coeff + u_boot(t,:);
        Xlag = [data_boot(p+t,:),Xlag(1:end-n_x)];
    end

    % estimate VAR on artificial data
    X_boot                 = lagmatrix(data_boot,1:p);
    X_boot                 = X_boot(p+1:end,:);
    Y_boot                 = data_boot(p+1:end,:);
    VAR_coeff_boot(:,:,i)  = [X_boot ones(length(X_boot),1)]\Y_boot;
    VAR_res_boot           = Y_boot-[X_boot ones(length(X_boot),1)]*VAR_coeff_boot(:,:,i);
    Sigma_u_boot(:,:,i)    = (VAR_res_boot'*VAR_res_boot)/(T_VAR-n_x*p-1);
end

% collect results
VAR_boot.VAR_coeff = VAR_coeff_boot;
VAR_boot.Sigma_u   = Sigma_u_boot;

end