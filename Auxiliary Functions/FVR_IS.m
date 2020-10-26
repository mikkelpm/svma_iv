function [FVR_LB,FVR_UB,FVR_true] = FVR_IS(yzt_aux,model,settings,alpha)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

% imports

n_y      = model.n_y;
FVR_hor  = settings.FVR_hor;
FVR_true = NaN(size(FVR_hor,2),n_y);
FVR_UB   = NaN(size(FVR_hor,2),n_y);
FVR_LB   = NaN(size(FVR_hor,2),n_y);

VMA_hor       = settings.VMA_hor;
Sigma_y_yhor  = yzt_aux.Sigma_y_big;
Cov_y         = yzt_aux.Cov_y;

% get forecasting variance

yzt_aux.Var_y_yhor = NaN(n_y,n_y,FVR_hor(end));

if settings.use_KF == 0

for hor = 1:FVR_hor(end)
    if hor >= 2
        for j = 2:VMA_hor-(hor-1)
            Sigma_y_yhor(1:n_y,1+(j-1)*n_y:j*n_y) = Cov_y{1+abs(1-(j+(hor-1)))}';
        end
        Sigma_y_yhor(1:n_y,1+(VMA_hor-hor+2-1)*n_y:end) = 0;
        for i = 2:VMA_hor-(hor-1)
            Sigma_y_yhor(1+(i-1)*n_y:i*n_y,1:n_y) = Cov_y{1+abs(i+(hor-1)-1)};
        end
        Sigma_y_yhor(1+(VMA_hor-hor+2-1)*n_y:end,1:n_y) = 0;
    end
    yzt_aux.Var_y_yhor(:,:,hor) = Sigma_y_yhor(1:n_y,1:n_y) ...
        - Sigma_y_yhor(1:n_y,n_y+1:end) * Sigma_y_yhor(n_y+1:end,n_y+1:end)^(-1) * Sigma_y_yhor(n_y+1:end,1:n_y);
end

else
    
A_KF          = yzt_aux.A_KF;
B_KF          = yzt_aux.B_KF;
C_KF          = yzt_aux.C_KF;
cvar_states_1 = yzt_aux.cvar_states_1;

for hor = 1:FVR_hor(end)
    if hor == 1
        cvar_states = cvar_states_1;
    else
        cvar_states = A_KF * cvar_states * A_KF' + B_KF * B_KF';
    end
    yzt_aux.Var_y_yhor(:,:,hor) = C_KF * cvar_states * C_KF';
end

end

%----------------------------------------------------------------
% Truth
%----------------------------------------------------------------

if isempty(alpha.alpha_true) == 0

for hor_indx = 1:size(FVR_hor,2)
    hor = FVR_hor(hor_indx);
    for j = 1:n_y
        FVR_true(hor_indx,j) = FVR_fun(j,hor,yzt_aux,alpha.alpha_true);
    end
end

else
    
FVR_true = [];
    
end

%----------------------------------------------------------------
% Upper Bound
%----------------------------------------------------------------

for hor_indx = 1:size(FVR_hor,2)
    hor = FVR_hor(hor_indx);
    for j = 1:n_y
        FVR_UB(hor_indx,j) = FVR_fun(j,hor,yzt_aux,alpha.alpha_LB);
    end
end

%----------------------------------------------------------------
% Lower Bound
%----------------------------------------------------------------

for hor_indx = 1:size(FVR_hor,2)
    hor = FVR_hor(hor_indx);
    for j = 1:n_y
        FVR_LB(hor_indx,j) = FVR_fun(j,hor,yzt_aux,alpha.alpha_UB);
    end
end

end