function VAR_pop = popVAR(model,settings)

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

VAR_pop.laglength = settings.VAR_poplaglength;
VAR_laglength     = VAR_pop.laglength;

%----------------------------------------------------------------
% VAR in Macro Aggregates + Instrument
%----------------------------------------------------------------

% derive Sigma and K

Sigma_states = eye(n_s);
K_states = zeros(n_s,n_x);

dist = 1;
tol = 10^(-10);
relax = 0.9;

while dist >= tol
    Sigma_upd = (A_x - K_states*C_x)*Sigma_states*(A_x-K_states*C_x)' + B_x*B_x' + K_states*(D_x*D_x')*K_states' - B_x * D_x'*K_states' - K_states * D_x * B_x';
    K_upd = (A_x * Sigma_upd * C_x' + B_x * D_x') * (C_x * Sigma_upd * C_x' + D_x * D_x')^(-1);
    
    dist1 = max(max(abs(Sigma_upd - Sigma_states)));
    dist2 = max(max(abs(K_upd - K_states)));
    
    Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
    K_states = relax * K_states + (1-relax) * K_upd;
    
    dist = max(dist1, dist2);
end

% get VAR

VAR_pop.VAR_coeff = NaN(VAR_laglength*n_x+1,n_x);

for i = 1:VAR_laglength
    VAR_pop.VAR_coeff(1+(i-1)*n_x:i*n_x,:) = (C_x * (A_x - K_states * C_x)^(i-1) * K_states)';
end

VAR_pop.Sigma_u = C_x * Sigma_states * C_x' + D_x * D_x';

%----------------------------------------------------------------
% VAR in Macro Aggregates
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

% get VAR

VAR_pop.VAR_coeff_y = NaN(VAR_laglength*n_y+1,n_y);

for i = 1:VAR_laglength
    VAR_pop.VAR_coeff_y(1+(i-1)*n_y:i*n_y,:) = (C_y * (A_y - K_states * C_y)^(i-1) * K_states)';
end

VAR_pop.Sigma_u_y = C_y * Sigma_states * C_y' + D_y * D_y';

% get covariance

VAR_pop.cov_uz = model.alpha * model.ABCD.D_y(:,1);
VAR_pop.gamma  = VAR_pop.cov_uz ./ sqrt(VAR_pop.cov_uz' * VAR_pop.Sigma_u_y^(-1) * VAR_pop.cov_uz);

end