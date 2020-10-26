function [IRF,FVD,M,tot_weights] = pop_analysis(model,settings);

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_s   = model.n_s;
n_x   = model.n_x;
n_eps = model.n_eps;
n_nu  = model.n_nu;
n_xi  = model.n_xi;
n_y   = model.n_y;
A_x   = model.ABCD.A_x;
B_x   = model.ABCD.B_x;
C_x   = model.ABCD.C_x;
D_x   = model.ABCD.D_x;
A_y   = model.ABCD.A_y;
B_y   = model.ABCD.B_y;
C_y   = model.ABCD.C_y;
D_y   = model.ABCD.D_y;
VMA_hor = settings.VMA_hor;
FVD_hor = settings.FVD_hor;

%----------------------------------------------------------------
% IRFs
%----------------------------------------------------------------

IRF = NaN(VMA_hor,n_y,n_xi);

for shock = 1:n_xi
    for i = 1:VMA_hor
        if i == 1
            resp = D_y';
            IRF(1,:,shock) = resp(shock,:);
        elseif i == 2
            resp = (C_y* B_y)';
            IRF(2,:,shock) = resp(shock,:);
        else
            resp = (C_y * A_y^(i-2) * B_y)';
            IRF(i,:,shock) = resp(shock,:);
        end
    end
end

%----------------------------------------------------------------
% FVDs
%----------------------------------------------------------------

FVD = NaN(size(FVD_hor,2),n_y);
for hor_indx = 1:size(FVD_hor,2)
    hor = FVD_hor(hor_indx);
    for var = 1:n_y
        num = sum(squeeze(IRF(1:hor,var,1)).^2);
        denom = sum(sum(squeeze(IRF(1:hor,var,:)).^2));
        FVD(hor_indx,var) = num/denom;
    end
end

%----------------------------------------------------------------
% Link Reduced-Form Errors and Structural Shocks
%----------------------------------------------------------------

% derive Sigma and K

Sigma_states = eye(n_s);
K_states = zeros(n_s,n_y);

dist = 1;
tol = 10^(-10);
relax = 0.9;

while dist >= tol
    Sigma_upd = (A_y - K_states*C_y)*Sigma_states*(A_y-K_states*C_y)' + B_y*B_y' + K_states*(D_y*D_y')*K_states' - B_y * D_y'*K_states' - K_states * D_y * B_y';
    K_upd = (A_y * Sigma_upd * C_y' + B_y * D_y') * (C_y * Sigma_upd * C_y' + D_y * D_y')^(-1);
    
    dist1 = max(max(abs(Sigma_upd - Sigma_states)));
    dist2 = max(max(abs(K_upd - K_states)));
    
    Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
    K_states = relax * K_states + (1-relax) * K_upd;
    
    dist = max(dist1, dist2);
end

% get the M matrices

aux = [A_y, zeros(n_s,n_s); K_states * C_y, A_y - K_states * C_y];
aux_lag = NaN(n_y,n_xi,VMA_hor);

for i = 1:VMA_hor
    aux_lag(:,:,i) = [C_y, -C_y] * aux^(i-1) * [B_y; K_states*D_y];
end

M = NaN(n_y,n_xi,VMA_hor+1);

M(:,:,1) = D_y;

for i = 2:VMA_hor+1
    M(:,:,i) = aux_lag(:,:,i-1); % the columns of this are the loadings
end

%----------------------------------------------------------------
% Total Available Weights
%----------------------------------------------------------------

tot_weights = NaN(n_xi,VMA_hor+1);

tot_weights(:,1) = diag(D_y' * (C_y * Sigma_states * C_y' + D_y * D_y')^(-1) * D_y);
for i = 2:VMA_hor+1
    tot_weights(:,i) = diag(([C_y, -C_y] * aux^(i-2) * [B_y; K_states*D_y])' ...
        * (C_y * Sigma_states * C_y' + D_y * D_y')^(-1) * ([C_y, -C_y] * aux^(i-2) * [B_y; K_states*D_y]));
end

end