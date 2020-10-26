%% EMPIRICAL APPLICATION: GK
% Mikkel Plagborg-Moller and Christian Wolf
% This version: 05/21/2018

%% HOUSEKEEPING

clc
clear all
close all

addpath('../Auxiliary Functions')

rng(2018)

%% SETTINGS

%----------------------------------------------------------------
% Simulation and VAR Estimation
%----------------------------------------------------------------

settings.select_VAR_simlaglength = 1;  % estimated VAR: pre-select or choose lag length?
settings.VAR_simlaglength        = 12; % estimated VAR lag length (if pre-set)
settings.max_simlaglength        = 24; % maximal lag length for estimated VAR (when chosen)
settings.use_KF                  = 1; % use Kalman filter for FVR computations?
settings.penalty                 = @(T) 2/T; % penalty term for lag length selection

%----------------------------------------------------------------
% Bootstrap
%----------------------------------------------------------------

settings.n_boot          = 10000; % bootstrap draws
settings.signif_level    = 0.1; % significance level
settings.optimopts       = optimoptions('fmincon', 'Display', 'notify'); % options for Stoye CI construction
settings.CI_for_R2_inv   = 1; % construct CI for R2_inv?
settings.CI_for_R2_recov = 1; % construct CI for R2_recov?
settings.CI_for_FVR      = 1; % construct CI for FVR?
settings.CI_for_FVD      = 0; % construct CI for FVD?

settings.fields          = {'alpha_LB', 'alpha_UB', 'R2_inv_LB', 'R2_inv_UB', 'R2_recov_LB', 'R2_recov_UB', 'FVR_LB', 'FVR_UB', 'FVD_LB'};
settings.fields_param    = {'alpha', 'R2_inv', 'R2_recov', 'FVR'};

%----------------------------------------------------------------
% Identified Set Characterization
%----------------------------------------------------------------

settings.VMA_hor        = 100; % maximal horizon in Wold/structural VMA representation; horizon M for bounds is set as function of that
settings.alpha_ngrid    = 1000; % grid points for lower bound on alpha
settings.bnd_recov      = 1; % naive recoverability-based lower bound on alpha only?
settings.FVR_hor        = 1:24; % horizons for FVR analysis
settings.FVD_hor        = 1:24; % horizons for FVD analysis

%% VAR REPRESENTATION

%----------------------------------------------------------------
% Get Data
%----------------------------------------------------------------

data.data_file = 'gk2015.csv'; % Data file
data.sample = [datetime('1990-1-1'), datetime('2012-06-01')]; % Sample end points
data.endo_vars = {'ff', 'dlogip', 'dlogcpi', 'ebp'}; % Endogenous (y) variables
data.iv_var = 'ff4_tc'; % External IV (z)

data.dat_table = readtable(data.data_file); % Read from file
data.sample_bool = datetime(data.dat_table.date) >= data.sample(1) & datetime(data.dat_table.date) <= data.sample(2); % Sample marker

data.time = data.dat_table.date(data.sample_bool); % Timestamps
data.Y = data.dat_table{data.sample_bool, data.endo_vars}; % Endogenous variables
data.Z = data.dat_table{data.sample_bool, data.iv_var}; % External IV

GK.data.x = [data.Y data.Z]; % Data matrix
GK.data.y = data.Y;
GK.data.z = data.Z;

clear data

%----------------------------------------------------------------
% Model Size
%----------------------------------------------------------------

GK.n_x   = size(GK.data.x,2);
GK.n_z   = 1;
GK.n_y   = GK.n_x - 1;

settings.T = size(GK.data.x,1);

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

disp('Getting the VAR representation...')
    
VAR_OLS = estimateVAR(GK.data.x,settings); 
VAR_sim = VAR_OLS;

disp('...done!')

%----------------------------------------------------------------
% Bootstrap VAR
%----------------------------------------------------------------

disp('Bootstrapping the VAR...')

VAR_boot = bootstrapVAR(VAR_OLS,GK,GK.data,settings);

disp('...done!')

%% CONSTRUCTION OF CI

%----------------------------------------------------------------
% OLS Point Estimate
%----------------------------------------------------------------

disp('Getting the OLS point estimates of the identified sets...')

yzt_aux    = get2ndmoments_VAR(VAR_OLS,GK,settings);
bounds_OLS = get_IS(yzt_aux,GK,settings);

disp('...done!')

%----------------------------------------------------------------
% Pre-Assignment
%----------------------------------------------------------------

for j=1:length(settings.fields)
    bounds_boot.(settings.fields{j}) = zeros([size(bounds_OLS.(settings.fields{j})) settings.n_boot]);
end

%----------------------------------------------------------------
% Get Identified Sets
%----------------------------------------------------------------

disp('Mapping each bootstrap draw into objects of interest...')

for i_boot = 1:settings.n_boot
    if mod(i_boot,50) == 0
        disp(['I am at bootstrap iteration ' int2str(i_boot)]);
    end
    VAR_sim.VAR_coeff = VAR_boot.VAR_coeff(:,:,i_boot);
    VAR_sim.Sigma_u   = VAR_boot.Sigma_u(:,:,i_boot);
    
    yzt_aux = get2ndmoments_VAR(VAR_sim,GK,settings);
    bounds = get_IS(yzt_aux,GK,settings);
    
    for j=1:length(settings.fields)
        bounds_boot.(settings.fields{j})(:,:,i_boot) = bounds.(settings.fields{j}); % Store bounds
    end
    
end

clear i_boot j VAR_sim yzt_aux bounds

disp('...done!')

%----------------------------------------------------------------
% Construct CIs
%----------------------------------------------------------------

disp('Constructing the confidence intervals...')

[CI.bounds_CI_IS,CI.bounds_CI_para] = CI_fun(bounds_boot,bounds_OLS,settings);

disp('...done!')

%% RESULTS

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

CI.alpha_CI_IS_LB = CI.bounds_CI_IS.lower.alpha_LB;
CI.alpha_CI_IS_UB = CI.bounds_CI_IS.upper.alpha_UB;

disp(['The confidence interval for the identified set of alpha is [' num2str(CI.alpha_CI_IS_LB) ',' num2str(CI.alpha_CI_IS_UB) ']'])

CI.alpha_CI_para_LB = CI.bounds_CI_para.lower.alpha;
CI.alpha_CI_para_UB = CI.bounds_CI_para.upper.alpha;

disp(['The confidence interval for alpha is [' num2str(CI.alpha_CI_para_LB) ',' num2str(CI.alpha_CI_para_UB) ']'])

%----------------------------------------------------------------
% R2
%----------------------------------------------------------------

% invertibility

if settings.CI_for_R2_inv == 1

CI.R2_inv_CI_IS_LB = CI.bounds_CI_IS.lower.R2_inv_LB;
CI.R2_inv_CI_IS_UB = CI.bounds_CI_IS.upper.R2_inv_UB;

disp(['The confidence interval for the identified set of R_0^2 is [' num2str(CI.R2_inv_CI_IS_LB) ',' num2str(CI.R2_inv_CI_IS_UB) ']'])
disp(['The point estimate for the identified set of R_0^2 is [' num2str(CI.bounds_CI_IS.OLS_biascorr.R2_inv_LB) ',' num2str(CI.bounds_CI_IS.OLS_biascorr.R2_inv_UB) ']'])

CI.R2_inv_CI_para_LB = CI.bounds_CI_para.lower.R2_inv;
CI.R2_inv_CI_para_UB = CI.bounds_CI_para.upper.R2_inv;

disp(['The confidence interval for R_0^2 is [' num2str(CI.R2_inv_CI_para_LB) ',' num2str(CI.R2_inv_CI_para_UB) ']'])

end

% recoverability

if settings.CI_for_R2_recov == 1

CI.R2_recov_CI_IS_LB = CI.bounds_CI_IS.lower.R2_recov_LB;
CI.R2_recov_CI_IS_UB = CI.bounds_CI_IS.upper.R2_recov_UB;

disp(['The confidence interval for the identified set of R_infty^2 is [' num2str(CI.R2_recov_CI_IS_LB) ',' num2str(CI.R2_recov_CI_IS_UB) ']'])
disp(['The point estimate for the identified set of R_infty^2 is [' num2str(CI.bounds_CI_IS.OLS_biascorr.R2_recov_LB) ',' num2str(CI.bounds_CI_IS.OLS_biascorr.R2_recov_UB) ']'])

CI.R2_recov_CI_para_LB = CI.bounds_CI_para.lower.R2_recov;
CI.R2_recov_CI_para_UB = CI.bounds_CI_para.upper.R2_recov;

disp(['The confidence interval for R_infty^2 is [' num2str(CI.R2_recov_CI_para_LB) ',' num2str(CI.R2_recov_CI_para_UB) ']'])

end

%----------------------------------------------------------------
% FVR
%----------------------------------------------------------------

if settings.CI_for_FVR
    
CI.FVR_CI_IS_LB = CI.bounds_CI_IS.lower.FVR_LB;
CI.FVR_CI_IS_UB = CI.bounds_CI_IS.upper.FVR_UB;

CI.FVR_CI_para_LB = CI.bounds_CI_para.lower.FVR;
CI.FVR_CI_para_UB = CI.bounds_CI_para.upper.FVR;

end

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

if settings.CI_for_FVD
    
CI.FVD_CI_para_LB   = CI.bounds_CI_IS.lower.FVD_LB;
CI.FVD_CI_para_UB   = CI.bounds_CI_IS.upper.FVD_LB;

end

%% PLOTS

% Plotting
plots.xticks = 3:3:24; % X axis ticks for FVR plot
plots.titles = {'FVR of Federal Funds Rate', 'FVR of Industrial Production Growth', 'FVR of CPI Growth', 'FVR of Excess Bond Premium'};
plots.xlabel = 'Horizon (Months)'; % X axis label for FVR plot
plots.ylabel = ''; % Y axis label for FVR plot

for i=1:GK.n_y
    
    % CI for identified set
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    plot_band(CI.FVR_CI_IS_LB(:,i), CI.FVR_CI_IS_UB(:,i), CI.bounds_CI_IS.OLS_biascorr.FVR_LB(:,i), CI.bounds_CI_IS.OLS_biascorr.FVR_UB(:,i), ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Estimate of identif. set', sprintf('%d%s', 100*(1-settings.signif_level), '\% conf. interval for identif. set')}, ...
              'YLim', [0 1], 'XLim', [1 24], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    grid on;

end