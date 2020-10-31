function ABCD = ABCD_fun(model);

% ABCD representation of Dynare model

% prepare inputs

decision = model.decision;
n_s      = model.n_s;
n_xi     = model.n_xi;
obs_y    = model.obs_y;
obs_x    = model.obs_x;

if ~strcmp(model.shock,'fg')

% ABCD representation in macro aggregates y

A_y = decision(1:n_s,1:n_s); A_y = A_y';
B_y = decision(n_s+1:n_s+n_xi,1:n_s); B_y = B_y';
C_y = decision(1:n_s,obs_y); C_y = C_y';
D_y = decision(n_s+1:n_s+n_xi,obs_y); D_y = D_y';

% ABCD representation in all variables x = (y,z)

A_x = decision(1:n_s,1:n_s); A_x = A_x';
B_x = decision(n_s+1:n_s+n_xi,1:n_s); B_x = B_x';
C_x = decision(1:n_s,obs_x); C_x = C_x';
D_x = decision(n_s+1:n_s+n_xi,obs_x); D_x = D_x';

elseif strcmp(model.shock,'fg')
    
% ABCD representation in macro aggregates y

A1_y = decision(1:n_s,1:14); A1_y = A1_y';
A2_y = zeros(2,n_s); A2_y(2,15) = 1;
A3_y = decision(1:n_s,15:n_s-2); A3_y = A3_y';
A_y = [A1_y; A2_y; A3_y];

B1_y = decision(n_s+1:n_s+n_xi,1:14); B1_y = B1_y';
B2_y = zeros(2,n_xi); B2_y(1,1) = 1;
B3_y = decision(n_s+1:n_s+n_xi,15:n_s-2); B3_y = B3_y';
B_y = [B1_y; B2_y; B3_y];

C_y = decision(1:n_s,obs_y); C_y = C_y';

D_y = decision(n_s+1:n_s+n_xi,obs_y); D_y = D_y';

% ABCD representation in all variables x = (y,z)

A1_x = decision(1:n_s,1:14); A1_x = A1_x';
A2_x = zeros(2,n_s); A2_x(2,15) = 1;
A3_x = decision(1:n_s,15:n_s-2); A3_x = A3_x';
A_x = [A1_x; A2_x; A3_x];

B1_x = decision(n_s+1:n_s+n_xi,1:14); B1_x = B1_x';
B2_x = zeros(2,n_xi); B2_x(1,1) = 1;
B3_x = decision(n_s+1:n_s+n_xi,15:n_s-2); B3_x = B3_x';
B_x = [B1_x; B2_x; B3_x];

C_x = decision(1:n_s,obs_x); C_x = C_x';

D_x = decision(n_s+1:n_s+n_xi,obs_x); D_x = D_x';

end

% collect results

ABCD.A_y = A_y;
ABCD.B_y = B_y;
ABCD.C_y = C_y;
ABCD.D_y = D_y;
ABCD.A_x = A_x;
ABCD.B_x = B_x;
ABCD.C_x = C_x;
ABCD.D_x = D_x;

end