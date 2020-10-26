function [data] = simulate_data(model,settings)

T      = settings.T;
n_eps  = model.n_eps;
n_nu   = model.n_nu;
n_xi   = model.n_xi;
n_s    = model.n_s;
n_x    = model.n_x;
A_x    = model.ABCD.A_x;
B_x    = model.ABCD.B_x;
C_x    = model.ABCD.C_x;
D_x    = model.ABCD.D_x;

shocks = normrnd(0,1,T,n_xi);

data_states = NaN(T+1,n_s);
data_obs    = NaN(T+1,n_x);

data_states(1,:) = zeros(1,n_s);
data_obs(1,:)    = zeros(1,n_x);

for t = 2:T+1
    data_states(t,:) = A_x * data_states(t-1,:)' + B_x * shocks(t-1,:)';
    data_obs(t,:)    = C_x * data_states(t-1,:)' + D_x * shocks(t-1,:)';
end

data_obs = data_obs(2:end,:);

data.y = data_obs(:,1:end-1);
data.z = data_obs(:,end); 
data.x = data_obs;