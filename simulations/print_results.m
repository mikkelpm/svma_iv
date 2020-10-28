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
R2_inv_true = nan(n_dgp,1);
R2_inv_cov = nan(n_dgp,1);
FVR_true = nan(n_dgp,n_hor);
FVR_cov = nan(n_dgp,n_hor);
FVR_cov_svar = nan(n_dgp,n_hor);

cover1 = @(cis,param) nanmean(cis(1,:)<=param & param<=cis(2,:), 2);
cover2 = @(cis,param) nanmean(reshape(cis(:,FVR_var,1,:)<=param & param<=cis(:,FVR_var,2,:), n_hor, []), 2);

for i_dgp=1:n_dgp % For each DGP...
    
    % DGP number
    dgp_no(i_dgp) = dgps{i_dgp}.model.dgp;
    
    % Degree of invertibility
    R2_inv_true(i_dgp) = dgps{i_dgp}.model.R2_inv;
    R2_inv_cov(i_dgp) = cover1(dgps{i_dgp}.svma_R2_inv_cis, R2_inv_true(i_dgp));
    
    % FVR
    FVR_true(i_dgp,:) = dgps{i_dgp}.model.FVR(:,FVR_var);
    FVR_cov(i_dgp,:) = cover2(dgps{i_dgp}.svma_FVR_cis, FVR_true(i_dgp,:)');
    FVR_cov_svar(i_dgp,:) = cover2(dgps{i_dgp}.svar_FVR_cis, FVR_true(i_dgp,:)');
    
end 

% Create table
results = table;
results.dgp = dgp_no;
results.R2_inv_true = R2_inv_true;
results.R2_inv_cov = R2_inv_cov;
for i_h=1:n_hor
    results.(sprintf('%s%d%s', 'FVR', i_h, '_true')) = FVR_true(:,i_h);
    results.(sprintf('%s%d%s', 'FVR', i_h, '_cov')) = FVR_cov(:,i_h);
    results.(sprintf('%s%d%s', 'FVR', i_h, '_cov_svar')) = FVR_cov_svar(:,i_h);
end

% Print to screen
disp(results);


clearvars i_dgp i_h;