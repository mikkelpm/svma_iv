function yzt_aux = get2ndmoments_VAR(VAR,model,settings)

% Compute autocovariances implied by reduced-form VAR

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

VAR_coeff = VAR.VAR_coeff;
Sigma_u   = VAR.Sigma_u;
n_x       = model.n_x;
n_y       = model.n_y;
p         = VAR.laglength;
VMA_hor   = settings.VMA_hor;

%----------------------------------------------------------------
% Wold Representation
%----------------------------------------------------------------

Theta_Wold                = NaN(VMA_hor,n_x,n_x);
VAR_coeff                 = [VAR_coeff(1:p*n_x,:);zeros(1+VMA_hor*n_x-p*n_x,n_x)];

for i = 1:VMA_hor
    if i == 1
        Theta_Wold(i,:,:) = Sigma_u^(0.5);
    else
        Theta_Wold(i,:,:) = zeros(n_x,n_x);
        for j = 1:i-1
            Theta_Wold(i,:,:) = squeeze(Theta_Wold(i,:,:)) + VAR_coeff(1+(j-1)*n_x:j*n_x,:)' * squeeze(Theta_Wold(i-j,:,:)); % just compute Cov(x_t, eps_{t-h})
        end
    end
end

Theta_Wold_y      = Theta_Wold(:,1:n_y,:);
Theta_Wold_z      = Theta_Wold(:,n_y+1,:);
Theta_Wold_zt     = squeeze(Theta_Wold(1,n_y+1,:));

%----------------------------------------------------------------
% Variances and Covariances
%----------------------------------------------------------------

Cov_x      = cell(VMA_hor,1);
for i = 1:VMA_hor
    Cov_x{i} = 0;
    for j = 1:(VMA_hor-i+1)
        Cov_x{i} = Cov_x{i} + squeeze(Theta_Wold(j,:,:)) * squeeze(Theta_Wold(j+i-1,:,:))';
    end
end

Cov_y      = cell(VMA_hor,1);
for i = 1:VMA_hor
    Cov_y{i} = Cov_x{i}(1:n_y,1:n_y);
end

Var_zt = Theta_Wold_zt' * Theta_Wold_zt;

Sigma_yzt = NaN(VMA_hor,n_y);
for i = 1:VMA_hor
    Sigma_yzt(i,:) = squeeze(Theta_Wold(i,1:end-1,:)) * squeeze(Theta_Wold(1,end,:));
end

Sigma_y_big = NaN(n_y*VMA_hor,n_y*VMA_hor);
for i = 1:VMA_hor
    for j = 1:VMA_hor
        if i > j
            Sigma_y_big(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_y{1+abs(i-j)};
        else
            Sigma_y_big(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_y{1+abs(i-j)}';
        end
    end
end

if settings.use_KF == 1

A_KF = zeros(n_x*p,n_x*p);
for i = 1:p
    A_KF(1:n_x,1+(i-1)*n_x:i*n_x) = VAR_coeff(1+(i-1)*n_x:i*n_x,:)';
end
A_KF(n_x+1:n_x*p,1:n_x*(p-1)) = eye(n_x*(p-1));

B_KF = zeros(n_x*p,n_x);
B_KF(1:n_x,:) = Sigma_u^(0.5);

C_KF = zeros(n_y,n_x*p);
C_KF(:,1:n_y) = eye(n_y);

T_KF = 1e2;
y_KF = zeros(T_KF,n_y);
[~,cvar_states_1] = Kalman_filter(A_KF,B_KF,C_KF,y_KF);

else
    
A_KF          = [];
B_KF          = [];
C_KF          = [];
cvar_states_1 = [];

end

%----------------------------------------------------------------
% Spectral Densities
%----------------------------------------------------------------

s_y   = @(omega) sd_fun(omega,Theta_Wold_y);
s_yzt = @(omega) cross_sd_fun(omega,Sigma_yzt);

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

yzt_aux.Theta_Wold        = Theta_Wold;
yzt_aux.Theta_Wold_y      = Theta_Wold_y;
yzt_aux.Theta_Wold_z      = Theta_Wold_z;
yzt_aux.Theta_Wold_zt     = Theta_Wold_zt;
yzt_aux.Cov_x             = Cov_x;
yzt_aux.Cov_y             = Cov_y;
yzt_aux.Var_zt            = Var_zt;
yzt_aux.Sigma_yzt         = Sigma_yzt;
yzt_aux.Sigma_u           = Sigma_u;
yzt_aux.s_y               = s_y;
yzt_aux.s_yzt             = s_yzt;
yzt_aux.Sigma_y_big       = Sigma_y_big;
yzt_aux.cvar_states_1     = cvar_states_1;
yzt_aux.A_KF              = A_KF;
yzt_aux.B_KF              = B_KF;
yzt_aux.C_KF              = C_KF;

end