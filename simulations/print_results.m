clear all;

% Display results of Monte Carlo study


%% Settings

% Output
FVR_var         = 2;            % Variable for which to report FVR

% File names
results_folder  = 'results';    % Folder with results
results_name    = 'dgp';        % Initial part of file name for .mat files


%% Load data

% Find .mat files in results folder
files = dir(fullfile(results_folder, strcat(results_name, '*.mat')));
n_dgp = length(files);
dgps = cell(1,n_dgp);
for i_dgp=1:n_dgp
    dgps{i_dgp} = load(fullfile(results_folder, files(i_dgp).name));
end


%% Create results table

FVR_hor = dgps{1}.settings.FVR_hor;
n_hor = length(FVR_hor);

dgp_no = nan(n_dgp,1);
true_R2_inv = nan(n_dgp,1);
cov_par_R2_inv = nan(n_dgp,1);
cov_set_R2_inv = nan(n_dgp,1);
true_FVR = nan(n_dgp,n_hor);
cov_par_FVR = nan(n_dgp,n_hor);
cov_set_FVR = nan(n_dgp,n_hor);
cov_par_FVR_svar = nan(n_dgp,n_hor);

cover1 = @(cis,param_lb,param_ub) mean(cis(1,:)<=param_lb & param_ub<=cis(2,:), 2);
cover2 = @(cis,param_lb,param_ub) mean(reshape(cis(:,FVR_var,1,:)<=param_lb(:,FVR_var) & param_ub(:,FVR_var)<=cis(:,FVR_var,2,:), n_hor, []), 2);

for i_dgp=1:n_dgp % For each DGP...
    
    % DGP number
    dgp_no(i_dgp) = dgps{i_dgp}.model.dgp;
    
    % Degree of invertibility
    true_R2_inv(i_dgp) = dgps{i_dgp}.model.R2_inv;
    cov_par_R2_inv(i_dgp) = cover1(dgps{i_dgp}.svma_R2_inv_cis, dgps{i_dgp}.model.R2_inv, dgps{i_dgp}.model.R2_inv);
    cov_set_R2_inv(i_dgp) = cover1(dgps{i_dgp}.svma_R2_inv_cis, dgps{i_dgp}.bounds_pop.R2_inv_LB, dgps{i_dgp}.bounds_pop.R2_inv_UB);
    
    % FVR
    true_FVR(i_dgp,:) = dgps{i_dgp}.model.FVR(:,FVR_var);
    cov_par_FVR(i_dgp,:) = cover2(dgps{i_dgp}.svma_FVR_cis, dgps{i_dgp}.model.FVR, dgps{i_dgp}.model.FVR);
    cov_set_FVR(i_dgp,:) = cover2(dgps{i_dgp}.svma_FVR_cis, dgps{i_dgp}.bounds_pop.FVR_LB, dgps{i_dgp}.bounds_pop.FVR_UB);
    cov_par_FVR_svar(i_dgp,:) = cover2(dgps{i_dgp}.svar_FVR_cis, dgps{i_dgp}.model.FVR, dgps{i_dgp}.model.FVR);
    
end 

% Create table
results = table;
results.dgp = dgp_no;
results.true_R2_inv = true_R2_inv;
results.cov_par_R2_inv = cov_par_R2_inv;
results.cov_set_R2_inv = cov_set_R2_inv;
for i_h=1:n_hor
    results.(sprintf('%s%d', 'true_FVR', i_h)) = true_FVR(:,i_h);
    results.(sprintf('%s%d', 'cov_par_FVR', i_h)) = cov_par_FVR(:,i_h);
    results.(sprintf('%s%d', 'cov_set_FVR', i_h)) = cov_set_FVR(:,i_h);
    results.(sprintf('%s%d%s', 'cov_par_FVR', i_h, '_svar')) = cov_par_FVR_svar(:,i_h);
end

% Print to screen
disp(results);


clearvars i_dgp i_h;