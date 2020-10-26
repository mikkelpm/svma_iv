%% MONTE CARLO STUDY
% Mikkel Plagborg-Moller and Christian Wolf
% This version: 05/21/2018

%% HOUSEKEEPING

clc
clear all
close all

addpath('../Auxiliary Functions')

rng('shuffle')

%% SPECIFY EXPERIMENTS

disp('I am doing a LP-IV analysis.')

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

% set benchmark

rho_y_benchmark    = 0.5;
rho_z_benchmark    = 0;
rho_zy_benchmark   = 0;
theta_benchmark    = 0;
sigma_nu_benchmark = 1;
T_benchmark        = 250;

% set selection vectors

n_exp     = 8;

rho_y_all    = rho_y_benchmark * ones(n_exp,1); rho_y_all(2) = 0.9;
rho_z_all    = rho_z_benchmark * ones(n_exp,1); rho_z_all(3) = 0.8;
rho_zy_all   = rho_zy_benchmark * ones(n_exp,1); rho_zy_all(3) = 0.3;
theta_all    = theta_benchmark * ones(n_exp,1); theta_all(4) = 0.5; theta_all(5) = 2;
sigma_nu_all = sigma_nu_benchmark * ones(n_exp,1); sigma_nu_all(6) = 2;
T_all        = T_benchmark * ones(n_exp,1); T_all(7) = 100; T_all(8) = 500;

% placeholders

coverage_para_all.alpha    = cell(n_exp,1);
coverage_IS.alpha          = cell(n_exp,1);
coverage_para_all.R2_inv   = cell(n_exp,1);
coverage_IS.R2_inv         = cell(n_exp,1);
coverage_para_all.R2_recov = cell(n_exp,1);
coverage_IS_all.R2_recov   = cell(n_exp,1);
coverage_para_all.FVR      = cell(n_exp,1);
coverage_IS_all.FVR        = cell(n_exp,1);
coverage_IS_all.FVD        = cell(n_exp,1);

%% SETTINGS

%----------------------------------------------------------------
% Simulation and VAR Estimation
%----------------------------------------------------------------

settings.VAR_poplaglength        = 50; % population VAR lag length
settings.select_VAR_simlaglength = 1;  % estimated VAR: pre-select or choose lag length?
settings.use_KF                  = 1; % use Kalman filter for FVR computations?

%----------------------------------------------------------------
% Monte Carlo
%----------------------------------------------------------------

settings.n_MC = 5000;

%----------------------------------------------------------------
% Bootstrap
%----------------------------------------------------------------

settings.n_boot          = 1000; % bootstrap draws
settings.signif_level    = 0.1; % significance level
settings.optimopts       = optimoptions('fmincon', 'Display', 'notify'); % options for Stoye CI construction
settings.CI_for_R2_inv   = 1; % construct CI for R2_inv?
settings.CI_for_R2_recov = 0; % construct CI for R2_recov?
settings.CI_for_FVR      = 1; % construct CI for FVR?
settings.CI_for_FVD      = 0; % construct CI for FVD?

settings.fields          = {'alpha_LB', 'alpha_UB', 'R2_inv_LB', 'R2_inv_UB', 'R2_recov_LB', 'R2_recov_UB', 'FVR_LB', 'FVR_UB', 'FVD_LB'};
settings.fields_param    = {'alpha', 'R2_inv', 'R2_recov', 'FVR'};

%----------------------------------------------------------------
% Identified Set Characterization
%----------------------------------------------------------------

settings.VMA_hor        = 50; % maximal horizon in Wold/structural VMA representation
settings.alpha_ngrid    = 1000; % grid points for lower bound on alpha
settings.bnd_recov      = 1; % naive recoverability-based lower bound on alpha only?
settings.FVR_hor        = [1 4]; % horizons for FVR analysis
settings.FVR_var        = 2;
settings.FVD_hor        = [1 4]; % horizons for FVD analysis
settings.FVD_var        = 2;

%% MAIN LOOP

for i_exp = 1:n_exp
    
disp(['I am doing experiment ' num2str(i_exp)])

%% SETTINGS

%----------------------------------------------------------------
% Simulation and VAR Estimation
%----------------------------------------------------------------

settings.T                       = T_all(i_exp); % sample size for simulation
settings.VAR_simlaglength        = round(sqrt(settings.T)/2); % estimated VAR lag length (if pre-set)
settings.max_simlaglength        = round(sqrt(settings.T)); % maximal lag length for estimated VAR (when chosen)
settings.penalty                 = @(T) 2/T; % penalty term for lag length selection

%% VAR REPRESENTATION

%----------------------------------------------------------------
% Model Specification
%----------------------------------------------------------------

% raw model parameters

MC_model.rho_y   = rho_y_all(i_exp);
MC_model.Xi_1    = [MC_model.rho_y 0; 0.5 0.5];
MC_model.Theta_0 = chol([1 0.8; 0.8 1],'lower');
MC_model.theta   = theta_all(i_exp);
MC_model.Theta_1 = MC_model.theta * MC_model.Theta_0;

MC_model.Psi_1    = rho_z_all(i_exp);
MC_model.Lambda_1 = rho_zy_all(i_exp) * ones(1,size(MC_model.Theta_0,1));
MC_model.alpha    = 1;
MC_model.sigma_nu = sigma_nu_all(i_exp);

% model size

MC_model.n_y   = size(MC_model.Xi_1,1);
MC_model.n_z   = 1;
MC_model.n_x   = MC_model.n_y + MC_model.n_z;
MC_model.n_eps = size(MC_model.Theta_0,1);
MC_model.n_nu  = 1;
MC_model.n_xi  = MC_model.n_eps + MC_model.n_nu;
MC_model.n_s   = MC_model.n_y + MC_model.n_z + MC_model.n_eps;

% ABCD representation

MC_model.ABCD.A_x = [MC_model.Xi_1 zeros(MC_model.n_y,MC_model.n_z) MC_model.Theta_1; ...
                     MC_model.Lambda_1 MC_model.Psi_1 zeros(MC_model.n_z,MC_model.n_eps); ...
                     zeros(MC_model.n_eps,MC_model.n_s)];
MC_model.ABCD.B_x = [MC_model.Theta_0 zeros(MC_model.n_y,1); ...
                     MC_model.alpha zeros(MC_model.n_z,MC_model.n_eps-1) MC_model.sigma_nu; ...
                     eye(MC_model.n_eps) zeros(MC_model.n_eps,1)];
MC_model.ABCD.C_x = [MC_model.Xi_1 zeros(MC_model.n_y,MC_model.n_z) MC_model.Theta_1; ...
                     MC_model.Lambda_1 MC_model.Psi_1 zeros(MC_model.n_z,MC_model.n_eps)];
MC_model.ABCD.D_x = [MC_model.Theta_0 zeros(MC_model.n_y,1); ...
                     MC_model.alpha zeros(MC_model.n_z,MC_model.n_eps-1) MC_model.sigma_nu];

MC_model.ABCD.A_y = MC_model.ABCD.A_x;
MC_model.ABCD.B_y = MC_model.ABCD.B_x;
MC_model.ABCD.C_y = [MC_model.Xi_1 zeros(MC_model.n_y,MC_model.n_z) MC_model.Theta_1];
MC_model.ABCD.D_y = [MC_model.Theta_0 zeros(MC_model.n_y,1)];

%----------------------------------------------------------------
% Describe Experiment
%----------------------------------------------------------------

disp('The current experiment sets:')
disp(['rho_y = ' num2str(MC_model.rho_y)])
disp(['rho_z = ' num2str(rho_z_all(i_exp))])
disp(['rho_zy = ' num2str(rho_zy_all(i_exp))])
disp(['theta = ' num2str(MC_model.theta)])
disp(['sigma_nu = ' num2str(MC_model.sigma_nu)])
disp(['T = ' num2str(settings.T)])

%----------------------------------------------------------------
% Get Population IRFs + FVDs
%----------------------------------------------------------------

[MC_model.IRF,MC_model.FVD] = pop_analysis(MC_model,settings);

%----------------------------------------------------------------
% Get Population Estimands
%----------------------------------------------------------------

VAR_pop    = popVAR(MC_model,settings);
yzt_aux    = get2ndmoments_VAR(VAR_pop,MC_model,settings);

bounds_pop          = get_IS(yzt_aux,MC_model,settings);
bounds_pop.R2_inv   = bounds_pop.R2_inv_UB * bounds_pop.alpha_LB^2/MC_model.alpha^2;
bounds_pop.R2_recov = bounds_pop.R2_recov_UB * bounds_pop.alpha_LB^2/MC_model.alpha^2;

MC_model.FVR = bounds_pop.FVR_UB * bounds_pop.alpha_LB^2/MC_model.alpha^2;

%% MONTE CARLO LOOP

parfor i_MC = 1:settings.n_MC
    
%----------------------------------------------------------------
% Simulate Data
%----------------------------------------------------------------

[data] = simulate_data(MC_model,settings);

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------
    
VAR_OLS = estimateVAR(data.x,settings); 
VAR_sim = VAR_OLS;

%----------------------------------------------------------------
% OLS Point Estimate
%----------------------------------------------------------------

yzt_aux    = get2ndmoments_VAR(VAR_OLS,MC_model,settings);
bounds_OLS = get_IS(yzt_aux,MC_model,settings);

%----------------------------------------------------------------
% Bootstrap VAR
%----------------------------------------------------------------

VAR_boot = bootstrapVAR(VAR_OLS,MC_model,data,settings);

%----------------------------------------------------------------
% Pre-Assignment
%----------------------------------------------------------------

bounds_boot = struct;
for j=1:length(settings.fields)
    bounds_boot.(settings.fields{j}) = zeros([size(bounds_OLS.(settings.fields{j})) settings.n_boot]);
end

%----------------------------------------------------------------
% Get Identified Sets
%----------------------------------------------------------------

for i_boot = 1:settings.n_boot
    VAR_sim.VAR_coeff = VAR_boot.VAR_coeff(:,:,i_boot);
    VAR_sim.Sigma_u   = VAR_boot.Sigma_u(:,:,i_boot);
    
    yzt_aux = get2ndmoments_VAR(VAR_sim,MC_model,settings);
    bounds = get_IS(yzt_aux,MC_model,settings);
    
    for j=1:length(settings.fields)
        bounds_boot.(settings.fields{j})(:,:,i_boot) = bounds.(settings.fields{j}); % Store bounds
    end   
end

%----------------------------------------------------------------
% Construct CIs
%----------------------------------------------------------------

[bounds_CI_IS,bounds_CI_para] = CI_fun(bounds_boot,bounds_OLS,settings);

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

% alpha

alpha_IS_LB(:,:,i_MC)    = bounds_CI_IS.lower.alpha_LB;
alpha_IS_UB(:,:,i_MC)    = bounds_CI_IS.upper.alpha_UB;
alpha_para_LB(:,:,i_MC)  = bounds_CI_para.lower.alpha;
alpha_para_UB(:,:,i_MC)  = bounds_CI_para.upper.alpha;

% R2: invertibility

R2_inv_IS_LB(:,:,i_MC)    = bounds_CI_IS.lower.R2_inv_LB;
R2_inv_IS_UB(:,:,i_MC)    = bounds_CI_IS.upper.R2_inv_UB;
R2_inv_para_LB(:,:,i_MC)  = bounds_CI_para.lower.R2_inv;
R2_inv_para_UB(:,:,i_MC)  = bounds_CI_para.upper.R2_inv;

% R2: recoverability

R2_recov_IS_LB(:,:,i_MC)    = bounds_CI_IS.lower.R2_recov_LB;
R2_recov_IS_UB(:,:,i_MC)    = bounds_CI_IS.upper.R2_recov_UB;
R2_recov_para_LB(:,:,i_MC)  = bounds_CI_para.lower.R2_recov;
R2_recov_para_UB(:,:,i_MC)  = bounds_CI_para.upper.R2_recov;

% FVR

FVR_IS_LB(:,:,i_MC)    = bounds_CI_IS.lower.FVR_LB;
FVR_IS_UB(:,:,i_MC)    = bounds_CI_IS.upper.FVR_UB;
FVR_para_LB(:,:,i_MC)  = bounds_CI_para.lower.FVR;
FVR_para_UB(:,:,i_MC)  = bounds_CI_para.upper.FVR;

% FVD

FVD_IS_LB(:,:,i_MC) = bounds_CI_IS.lower.FVD_LB;
FVD_IS_UB(:,:,i_MC) = 1;

end

% clean up results

CI_IS   = struct;
CI_para = struct;

CI_IS.alpha_LB   = alpha_IS_LB;
CI_IS.alpha_UB   = alpha_IS_UB;
CI_para.alpha_LB = alpha_para_LB;
CI_para.alpha_UB = alpha_para_UB;

CI_IS.R2_inv_LB   = R2_inv_IS_LB;
CI_IS.R2_inv_UB   = R2_inv_IS_UB;
CI_para.R2_inv_LB = R2_inv_para_LB;
CI_para.R2_inv_UB = R2_inv_para_UB;

CI_IS.R2_recov_LB   = R2_recov_IS_LB;
CI_IS.R2_recov_UB   = R2_recov_IS_UB;
CI_para.R2_recov_LB = R2_recov_para_LB;
CI_para.R2_recov_UB = R2_recov_para_UB;

CI_IS.FVR_LB   = FVR_IS_LB;
CI_IS.FVR_UB   = FVR_IS_UB;
CI_para.FVR_LB = FVR_para_LB;
CI_para.FVR_UB = FVR_para_UB;

CI_IS.FVD_LB   = FVD_IS_LB;
CI_IS.FVD_UB   = FVD_IS_UB;

clear bounds bounds_boot bounds_CI_IS bounds_CI_para bounds_OLS i_boot i_MC j VAR_boot VAR_OLS VAR_sim yzt_aux ...
    alpha_IS_LB alpha_IS_UB alpha_para_LB alpha_para_UB R2_inv_IS_LB R2_inv_IS_UB R2_inv_para_LB R2_inv_para_UB ...
    R2_recov_IS_LB R2_recov_IS_UB R2_recov_para_LB R2_recov_para_UB FVR_IS_LB FVR_IS_UB FVR_para_LB FVR_para_UB ...
    FVD_IS_LB FVD_IS_UB    

%% RESULTS

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

coverage_para.alpha = (CI_para.alpha_LB <= MC_model.alpha & CI_para.alpha_UB >= MC_model.alpha);

disp(['The fraction of confidence intervals for alpha covering the truth is ' num2str(sum(coverage_para.alpha)/settings.n_MC)])

coverage_para_all.alpha{i_exp} = coverage_para.alpha;

coverage_IS.alpha = (CI_IS.alpha_LB <= bounds_pop.alpha_LB & CI_IS.alpha_UB >= bounds_pop.alpha_UB);

disp(['The fraction of confidence intervals for the identified set of alpha covering the truth is ' num2str(sum(coverage_IS.alpha)/settings.n_MC)])

coverage_IS_all.alpha{i_exp} = coverage_IS.alpha;

%----------------------------------------------------------------
% R2
%----------------------------------------------------------------

% invertibility

if settings.CI_for_R2_inv == 1

coverage_para.R2_inv = (CI_para.R2_inv_LB <= bounds_pop.R2_inv & CI_para.R2_inv_UB >= bounds_pop.R2_inv);

disp(['The fraction of confidence intervals for R2_inv covering the truth is ' num2str(sum(coverage_para.R2_inv)/settings.n_MC)])

coverage_para_all.R2_inv{i_exp} = coverage_para.R2_inv;

coverage_IS.R2_inv = (CI_IS.R2_inv_LB <= bounds_pop.R2_inv_LB & CI_IS.R2_inv_UB >= bounds_pop.R2_inv_UB);

disp(['The fraction of confidence intervals for the identified set of R2_inv covering the truth is ' num2str(sum(coverage_IS.R2_inv)/settings.n_MC)])

coverage_IS_all.R2_inv{i_exp} = coverage_IS.R2_inv;

end

% recoverability

if settings.CI_for_R2_recov == 1

if settings.bnd_recov == 0

coverage_para.R2_recov = (CI_para.R2_recov_LB <= bounds_pop.R2_recov & CI_para.R2_recov_UB >= bounds_pop.R2_recov);
coverage_IS.R2_recov = (CI_IS.R2_recov_LB <= bounds_pop.R2_recov_LB & CI_IS.R2_recov_UB >= bounds_pop.R2_recov_UB);

else
    
coverage_para.R2_recov = (CI_para.R2_recov_LB <= bounds_pop.R2_recov);
coverage_IS.R2_recov = (CI_IS.R2_recov_LB <= bounds_pop.R2_recov_LB);

end

coverage_para_all.R2_recov{i_exp} = coverage_para.R2_recov;
coverage_IS_all.R2_recov{i_exp} = coverage_IS.R2_recov;

disp(['The fraction of confidence intervals for R2_recov covering the truth is ' num2str(sum(coverage_para.R2_recov)/settings.n_MC)])
disp(['The fraction of confidence intervals for the identified set of R2_recov covering the truth is ' num2str(sum(coverage_IS.R2_recov)/settings.n_MC)])

end

%----------------------------------------------------------------
% FVR
%----------------------------------------------------------------

if settings.CI_for_FVR
    
% get rid of first variable

CI_para.FVR_LB = squeeze(permute(CI_para.FVR_LB(:,settings.FVR_var,:),[3 1 2]));
CI_para.FVR_UB = squeeze(permute(CI_para.FVR_UB(:,settings.FVR_var,:),[3 1 2]));
CI_IS.FVR_LB   = squeeze(permute(CI_IS.FVR_LB(:,settings.FVR_var,:),[3 1 2]));
CI_IS.FVR_UB   = squeeze(permute(CI_IS.FVR_UB(:,settings.FVR_var,:),[3 1 2]));

% report coverage

for hor_ind = 1:length(settings.FVR_hor)
    
    hor = settings.FVR_hor(hor_ind);
    coverage_para.FVR{hor_ind} = (CI_para.FVR_LB(:,hor_ind) <= MC_model.FVR(hor_ind,settings.FVR_var) & CI_para.FVR_UB(:,hor_ind) >= MC_model.FVR(hor_ind,settings.FVR_var));
    disp(['The fraction of confidence intervals for FVR at horizon ' num2str(hor) ' covering the truth is ' num2str(sum(coverage_para.FVR{hor_ind})/settings.n_MC)])
    coverage_IS.FVR{hor_ind} = (CI_IS.FVR_LB(:,hor_ind) <= bounds_pop.FVR_LB(hor_ind,settings.FVR_var) & CI_IS.FVR_UB(:,hor_ind) >= bounds_pop.FVR_UB(hor_ind,settings.FVR_var));
    disp(['The fraction of confidence intervals for the identified set of FVR at horizon ' num2str(hor) ' covering the truth is ' num2str(sum(coverage_IS.FVR{hor_ind})/settings.n_MC)])
    
end

coverage_para_all.FVR{i_exp} = coverage_para.FVR;
coverage_IS_all.FVR{i_exp} = coverage_para.FVR;

end

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

if settings.CI_for_FVD
    
% get rid of first variable
    
CI_IS.FVD_LB   = squeeze(permute(CI_IS.FVD_LB(:,settings.FVD_var,:),[3 1 2]));
CI_IS.FVD_UB   = squeeze(permute(CI_IS.FVD_UB(:,settings.FVD_var,:),[3 1 2]));

% report coverage

for hor_ind = 1:length(settings.FVR_hor)
    
    hor = settings.FVD_hor(hor_ind);
    coverage_IS.FVD{hor_ind} = (CI_IS.FVD_LB(:,hor_ind) <= MC_model.FVD(hor_ind,settings.FVD_var));
    disp(['The fraction of confidence intervals for the identified set of FVD at horizon ' num2str(hor) ' covering the truth is ' num2str(sum(coverage_IS.FVD{hor_ind})/settings.n_MC)])
    
end

coverage_IS_all.FVD{i_exp} = coverage_para.FVD;

end

end