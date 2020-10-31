function [R2_LB,R2_UB,R2_true] = R2_IS(yzt_aux,model,settings,R2_hor,alpha);

% Identified set for degree of invertibility out to time t+l

%----------------------------------------------------------------
% Get Inputs
%----------------------------------------------------------------

Var_zt        = yzt_aux.Var_zt;
Sigma_yzt     = yzt_aux.Sigma_yzt;
Sigma_y_big   = yzt_aux.Sigma_y_big;
VMA_hor       = settings.VMA_hor;
n_y           = model.n_y;
alpha_true    = alpha.alpha_true;
alpha_LB      = alpha.alpha_LB;
alpha_UB      = alpha.alpha_UB;

%----------------------------------------------------------------
% Auxiliary Computations
%----------------------------------------------------------------

Cov_zty_lagged = zeros(1,n_y*VMA_hor);
for i = 1:R2_hor
    Cov_zty_lagged(1,1+(i-1)*n_y:i*n_y) = Sigma_yzt(R2_hor-i+1,:);
end

Var_zt_yhor = Var_zt - Cov_zty_lagged * Sigma_y_big^(-1) * Cov_zty_lagged';

Rtilde2 = 1 - Var_zt_yhor/Var_zt;

%----------------------------------------------------------------
% Truth
%----------------------------------------------------------------

if isempty(alpha.alpha_true) == 0

R2_true = Var_zt/alpha_true^2 * Rtilde2;

else
    
R2_true = [];
    
end

%----------------------------------------------------------------
% Upper Bound
%----------------------------------------------------------------

R2_UB = Var_zt/alpha_LB^2 * Rtilde2;

%----------------------------------------------------------------
% Lower Bound
%----------------------------------------------------------------

R2_LB = Var_zt/alpha_UB^2 * Rtilde2;

end