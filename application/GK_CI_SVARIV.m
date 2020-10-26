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

settings.select_VAR_simlaglength = 0;  % estimated VAR: pre-select or choose lag length?
settings.VAR_simlaglength        = 6; % estimated VAR lag length (if pre-set)
settings.max_simlaglength        = 24; % maximal lag length for estimated VAR (when chosen)
settings.penalty                 = @(T) 2/T; % penalty term for lag length selection

%----------------------------------------------------------------
% Bootstrap
%----------------------------------------------------------------

settings.n_boot          = 10000; % bootstrap draws
settings.signif_level    = 0.1; % significance level

%----------------------------------------------------------------
% Identified Set Characterization
%----------------------------------------------------------------

settings.VMA_hor        = 50; % maximal horizon in Wold/structural VMA representation
settings.FVR_hor        = 1:24; % horizons for FVR analysis
settings.FVD_hor        = settings.FVR_hor; % horizons for FVD analysis
settings.IRF_hor        = 1:24; % horizons for IRF analysis

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
    
VAR_OLS = estimateVAR_IV(GK.data.y,GK.data.z,settings);  
VAR_sim = VAR_OLS;

disp('...done!')

%----------------------------------------------------------------
% Bootstrap VAR
%----------------------------------------------------------------

disp('Bootstrapping the VAR...')

VAR_boot = bootstrapVAR_IV(VAR_OLS,GK,GK.data,settings);

disp('...done!')

%% CONSTRUCTION OF CI

%----------------------------------------------------------------
% OLS Point Estimate
%----------------------------------------------------------------

disp('Getting the OLS point estimates of the identified sets...')

[SVARIV_OLS.IRF,SVARIV_OLS.FVD] = SVARIV_analysis(VAR_OLS,GK,settings);

disp('...done!')

%----------------------------------------------------------------
% Pre-Assignment
%----------------------------------------------------------------

SVARIV_boot.IRF = NaN(length(settings.IRF_hor),GK.n_y,settings.n_boot);
SVARIV_boot.FVD = NaN(length(settings.FVD_hor),GK.n_y,settings.n_boot);

%----------------------------------------------------------------
% Get Identified Sets
%----------------------------------------------------------------

disp('Mapping each bootstrap draw into objects of interest...')

for i_boot = 1:settings.n_boot
    if mod(i_boot,50) == 0
        disp(['I am at bootstrap iteration ' int2str(i_boot)]);
    end
    VAR_sim.VAR_coeff_y = VAR_boot.VAR_coeff_y(:,:,i_boot);
    VAR_sim.Sigma_u_y   = VAR_boot.Sigma_u_y(:,:,i_boot);
    VAR_sim.gamma       = VAR_boot.gamma(:,i_boot);
    
    [SVARIV_boot.IRF(:,:,i_boot),SVARIV_boot.FVD(:,:,i_boot)] = SVARIV_analysis(VAR_sim,GK,settings);
    
end

clear i_boot VAR_sim

disp('...done!')

%----------------------------------------------------------------
% Construct CIs
%----------------------------------------------------------------

disp('Constructing the confidence intervals...')

[IRF_CI,FVD_CI] = CI_SVARIV_fun(SVARIV_OLS,SVARIV_boot,settings);

disp('...done!')

%% PLOTS

plots.xticks = 3:3:24; % X axis ticks for FVR plot
plots.titles = {'FVD of Federal Funds Rate', 'FVD of Industrial Production Growth', 'FVD of CPI Growth', 'FVD of Excess Bond Premium'};
plots.xlabel = 'Horizon (Months)'; % X axis label for FVR plot
plots.ylabel = ''; % Y axis label for FVR plot

for i=1:GK.n_y
    
    % CI
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    plot_band(FVD_CI.lower(:,i), FVD_CI.upper(:,i), FVD_CI.biascorr(:,i), [], ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Point estimate', sprintf('%d%s', 100*(1-settings.signif_level), '\% conf. interval')}, ...
              'YLim', [0 1], 'XLim', [1 24], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    grid on;
          
end