function [alpha_LB,alpha_UB,alpha_true,alpha_plot] = alpha_IS(yzt_aux,model,settings);

%----------------------------------------------------------------
% Get Inputs
%----------------------------------------------------------------

Var_zt    = yzt_aux.Var_zt;
s_y       = yzt_aux.s_y;
s_yzt     = yzt_aux.s_yzt;
n_y       = model.n_y;
ngrid     = settings.alpha_ngrid;
bnd_recov = settings.bnd_recov;
hor_pred  = round(settings.VMA_hor/2)-1; % largest you can go with a two-sided prediction

%----------------------------------------------------------------
% Truth
%----------------------------------------------------------------

if isfield(model,'alpha') == 1

alpha_true = model.alpha;

else
    
alpha_true = [];

end

%----------------------------------------------------------------
% Upper Bound
%----------------------------------------------------------------

alpha_UB = sqrt(Var_zt);

%----------------------------------------------------------------
% Lower Bound
%----------------------------------------------------------------

if bnd_recov == 0

alpha_LB_fun = @(omega) sqrt(real(ctranspose(s_yzt(omega)) * s_y(omega)^(-1) * s_yzt(omega)));
omega_grid = linspace(0,2*pi,ngrid)';
alpha_LB_vals = NaN(ngrid,1);
for i = 1:ngrid
    alpha_LB_vals(i) = alpha_LB_fun(omega_grid(i));
end
alpha_LB = max(alpha_LB_vals);

alpha_plot.alpha_LB_vals = alpha_LB_vals;
alpha_plot.omega_grid    = omega_grid;

elseif bnd_recov == 1

Var_y = NaN((2*hor_pred+1) * n_y, (2*hor_pred+1) * n_y);
for i = 1:2*hor_pred+1
    for j = 1:2*hor_pred+1
        if i > j
            Var_y(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = yzt_aux.Cov_y{1+abs(i-j)};
        else
            Var_y(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = yzt_aux.Cov_y{1+abs(i-j)}';
        end
    end
end
Cov_zy = zeros(1,(2*hor_pred+1)*n_y);
for i = 1:hor_pred+1
    Cov_zy(1,1+(i-1)*n_y:i*n_y) = yzt_aux.Sigma_yzt(hor_pred+1-i+1,:);
end

alpha_LB = sqrt(Cov_zy * Var_y^(-1) * Cov_zy'); 

alpha_plot.alpha_LB_vals = [];
alpha_plot.omega_grid    = [];

end

end