clear all;
addpath('../functions');

% Monte Carlo study of SVMA-IV and SVAR-IV procedures


%% Settings

% DGP
model.dgp           = 0;        % DGP number 0-8 (0=baseline), see below

% Parameters of interest
settings.FVR_hor    = [1 4];    % Horizons for FVR/FVD

% MCMC settings
mcmc.n_rep          = 5e3;      % Number of Monte Carlo repetitions
rng(2018, 'twister');           % Seed random number generator
rng(2018+model.dgp, 'twister'); % Seed random number generator

% Inference settings
settings.ic         = 'aic';    % Lag length information criterion
settings.n_boot     = 1e3;      % No. of bootstrap draws
settings.signif     = 0.1;      % Significance level

% Numerical settings
settings.VMA_hor    = 50;       % Maximal horizon in Wold/structural VMA representation
settings.use_kalman = true;     % Use Kalman filter for FVR computations?

% File name
mcmc.save_folder    = 'results'; % Folder in which to store files
mcmc.save_name      = sprintf('%s%d', 'dgp', model.dgp); % File name


%% Define true DGP

% Baseline parameters
model.rho_y    = 0.5;
model.rho_z    = 0;
model.rho_zy   = 0;
model.zeta     = 0;
model.sigma_nu = 1;
model.T        = 250;
model.p        = 1; % Can be at most 4 in this code

% Experiment-specific parameters
switch model.dgp
    case 1
        model.rho_y = 0.9;
    case 2
        model.rho_z = 0.8;
        model.rho_zy = 0.3;
    case 3
        model.zeta = 0.5;
    case 4
        model.zeta = 2;
    case 5
        model.sigma_nu = 2;
    case 6
        model.T = 100;
    case 7
        model.T = 500;
    case 8
        model.p = 4;
end

disp('True parameters:');
fields = fieldnames(model);
for i_f=1:length(fields)
    disp([fields{i_f}, ' = ', num2str(model.(fields{i_f}))]);
end
clearvars fields i_f;


%% ABCD representation

% Model parameters
model.Xi_1      = [model.rho_y 0; 0.5 0.5];
for i_p=2:4
    model.(sprintf('%s%d','Xi_',i_p)) = (i_p<=model.p)/i_p^2*model.Xi_1;
end
clearvars i_p;
model.Theta_0   = chol([1 0.8; 0.8 1],'lower');
model.Theta_1   = model.zeta * model.Theta_0;
model.Psi_1     = model.rho_z;
model.Lambda_1  = model.rho_zy * ones(1,size(model.Theta_0,1));
model.alpha     = 1;

% Dimensions
model.n_y   = size(model.Xi_1,1);
model.n_z   = 1;
model.n_x   = model.n_y + model.n_z;
model.n_eps = size(model.Theta_0,1);
model.n_nu  = 1;
model.n_xi  = model.n_eps + model.n_nu;
model.n_s   = 4*model.n_y+model.n_z+model.n_eps;

% ABCD representation
% See Fernandez-Villaverde, Rubio-Ramirez, Sargent & Watson (AER 2007)
model.ABCD.A_x = [model.Xi_1 model.Xi_2 model.Xi_3 model.Xi_4 zeros(model.n_y,model.n_z) model.Theta_1; ...
                     eye(model.n_y) zeros(model.n_y,3*model.n_y+model.n_z+model.n_eps); ...
                     zeros(model.n_y,model.n_y) eye(model.n_y) zeros(model.n_y,2*model.n_y+model.n_z+model.n_eps); ...
                     zeros(model.n_y,2*model.n_y) eye(model.n_y) zeros(model.n_y,model.n_y+model.n_z+model.n_eps); ...
                     model.Lambda_1 zeros(model.n_z,3*model.n_y) model.Psi_1 zeros(model.n_z,model.n_eps); ...
                     zeros(model.n_eps,4*model.n_y+model.n_z+model.n_eps)];
model.ABCD.B_x = [model.Theta_0 zeros(model.n_y,1); ...
                     zeros(3*model.n_y,model.n_eps+1); ...
                     model.alpha zeros(model.n_z,model.n_eps-1) model.sigma_nu; ...
                     eye(model.n_eps) zeros(model.n_eps,1)];
model.ABCD.C_x = [model.Xi_1 model.Xi_2 model.Xi_3 model.Xi_4 zeros(model.n_y,model.n_z) model.Theta_1; ...
                     model.Lambda_1 zeros(model.n_z,3*model.n_y) model.Psi_1 zeros(model.n_z,model.n_eps)];
model.ABCD.D_x = [model.Theta_0 zeros(model.n_y,1); ...
                     model.alpha zeros(model.n_z,model.n_eps-1) model.sigma_nu];

model.ABCD.A_y = model.ABCD.A_x;
model.ABCD.B_y = model.ABCD.B_x;
model.ABCD.C_y = [model.Xi_1 model.Xi_2 model.Xi_3 model.Xi_4 zeros(model.n_y,model.n_z) model.Theta_1];
model.ABCD.D_y = [model.Theta_0 zeros(model.n_y,1)];

settings.T = model.T;


%% Compute population IRFs + FVDs

settings.CI_for_R2_inv   = 1; % Construct CI for R2_inv?
settings.CI_for_R2_recov = 0; % Construct CI for R2_recov?
settings.CI_for_FVR      = 1; % Construct CI for FVR?
settings.CI_for_FVD      = 0; % Construct CI for FVD?
settings.VAR_poplaglength = settings.VMA_hor;
settings.FVD_hor = settings.FVR_hor;
settings.use_KF = settings.use_kalman;
settings.alpha_ngrid = 1;
settings.bnd_recov = 1;       % DGP is always recoverable

disp('Doing the population analysis...');

[model.IRF,model.FVD] = pop_analysis(model,settings); % True IRF/FVD

VAR_pop         = popVAR(model,settings);
yzt_aux         = get2ndmoments_VAR(VAR_pop,model,settings);
bounds_pop      = get_IS(yzt_aux,model,settings); % True bounds
model.R2_inv    = bounds_pop.R2_inv_UB * bounds_pop.alpha_LB^2/model.alpha^2; % True R2_inv
model.FVR       = bounds_pop.FVR_UB * bounds_pop.alpha_LB^2/model.alpha^2; % True FVR

disp('...done!');


%% Monte Carlo loop

% Parallel computing object
delete(gcp('nocreate'));
num_workers = str2num(getenv('SLURM_CPUS_PER_TASK'));
if ~isempty(num_workers)
    poolobj = parpool('local', num_workers);
else
    poolobj = parpool('local');
end
clearvars num_workers;

% Results arrays
laglengths = nan(1,mcmc.n_rep);
svma_alpha_cis = nan(2,mcmc.n_rep);
svma_R2_inv_cis = nan(2,mcmc.n_rep);
svma_FVR_cis = nan(length(settings.FVR_hor),model.n_y,2,mcmc.n_rep);
svar_FVR_cis = nan(length(settings.FVR_hor),model.n_y,2,mcmc.n_rep);

rngs = randi(2^32-1,mcmc.n_rep,1); % Random seeds for each repetition

disp('Running the MC loop...');
timer = tic;

parfor i_MC = 1:mcmc.n_rep
% for i_MC = 1:mcmc.n_rep
    
    % Seed RNG
    rng(rngs(i_MC), 'twister');
    
    % Simulate data
    the_data = simulate_data(model,settings);

    % SVMA-IV inference
    [the_bounds,~,~,~,the_VAR_OLS] ...
        = SVMAIV_estim(the_data.y, the_data.z, ...
                      'ic', settings.ic, ...
                      'ic_max', sqrt(settings.T), ... % Max lag length as function of sample size
                      'n_boot', settings.n_boot, ...
                      'signif', settings.signif, ...
                      'compute_R2_recov', false, ...
                      'compute_FVD', false, ...
                      'horiz', settings.FVR_hor, ...
                      'ci_param', false, ...
                      'optim_opts', [], ...
                      'use_kalman', settings.use_kalman, ...
                      'VMA_hor', settings.VMA_hor, ...
                      'verbose', false);
    
    % SVAR-IV inference
    [~,the_svar_FVD] = SVARIV_estim(the_data.y, the_data.z, ...
                      'p', the_VAR_OLS.laglength, ... % Same lag length as SVMA procedure
                      'n_boot', settings.n_boot, ...
                      'signif', settings.signif, ...
                      'horiz', settings.FVR_hor, ...
                      'verbose', false);
    
    % Store results
    laglengths(i_MC)            = the_VAR_OLS.laglength;
    svma_alpha_cis(:,i_MC)      = [the_bounds.ci.lower.alpha, the_bounds.ci.upper.alpha];
    svma_R2_inv_cis(:,i_MC)     = [the_bounds.ci.lower.R2_inv, the_bounds.ci.upper.R2_inv];
    svma_FVR_cis(:,:,:,i_MC)    = cat(3,the_bounds.ci.lower.FVR, the_bounds.ci.upper.FVR);
    svar_FVR_cis(:,:,:,i_MC)    = cat(3,the_svar_FVD.ci.lower, the_svar_FVD.ci.upper);
    
    % Print progress
    if mod(i_MC,ceil(mcmc.n_rep/50))==0
        fract = i_MC/mcmc.n_rep;
        offs = floor(50*fract);
        fprintf(['%' num2str(offs+3) 'd%s\n'], round(100*fract), '%');
    end

end

elapsed_time = toc(timer);
disp('...done!');

disp('Elapsed minutes:');
disp(elapsed_time/60);

clearvars timer;


%% Save results

mkdir(mcmc.save_folder);
save(fullfile(mcmc.save_folder, mcmc.save_name));

delete(gcp('nocreate'));

