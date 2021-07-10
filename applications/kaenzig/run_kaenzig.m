clc
clear all;
close all;

% Empirical application: Kaenzig (AER, 2021)


%% Load data

% Settings
data.endo_vars  = {'oil_price','oil_production','oil_inventories','world_ip','neer','ip','cpi','ffr','vxo','tot'}; % Abbreviations of endogenous variables
data.diff       = [0,1,1,1,1,1,1,0,0,0]; % =1: take first differences of endogenous variable
data.iv         = 'oil_surprise';
data.smpl_start = '1983M03'; % Sample start point
data.smpl_end   = '2017M12'; % Sample end point

% Endogenous variables
dat_endo = load('dataQuantM'); % Endogenous variables
data.sampleDates = dat_endo.sampleDates;
data.smplStartInd = find(strcmp(data.sampleDates,data.smpl_start));
data.smplEndInd   = find(strcmp(data.sampleDates,data.smpl_end));
data.Y = dat_endo.data(data.smplStartInd:data.smplEndInd,:);

% First differences
for i = 1:size(data.Y,2)
    if data.diff(i) == 1
        data.Y(:,i) = [NaN;12 * (data.Y(2:end,i) - data.Y(1:end-1,i))];
    end
end

% IV
dat_iv=load('OilSurprisesMLog'); % IV
proxyRaw = dat_iv.oilProxiesWTIM(:,14);
data.Z = proxyRaw(data.smplStartInd:data.smplEndInd,:);

% Adjust sample
if max(data.diff) == 1
    data.Y = data.Y(2:end,:);
    data.Z = data.Z(2:end,:);
end

clear sampleDates dataEndo oilProxiesWTIM proxyRaw;


%% SVMA-IV inference

disp('*** SVMA-IV analysis ***');

% Preliminaries
addpath('../../functions');     % Add folder with SVMA-IV analysis functions
rng(2018);                      % Seed random number generator (for bootstrap)

% Estimation settings (see other optional settings in "functions/SVMAIV_estim.m")
settings = {'p', 12;            % Lag length
            'n_boot', 1e3;      % Number of bootstrap samples
            'signif', 0.1;      % Significance level
            'horiz', 1:48}';    % Horizons of FVR to report

% Run inference routines
[bounds, id_recov, inv_test, settings_struct] = SVMAIV_estim(data.Y, data.Z, settings{:});


%% Display pre-test for invertibility

disp('Invertibility pre-test p-value: all equations jointly');
disp(inv_test.pval.all);

disp('Invertibility pre-test p-value: each equation separately');
disp(inv_test.pval.eqns);


%% Display bounds on alpha and degree of invertibility/recoverability

% Scale parameter
disp('Bound estimates: alpha');
disp([bounds.estim.lower.alpha bounds.estim.upper.alpha]);
disp('Confidence interval: alpha');
disp([bounds.ci.lower.alpha bounds.ci.upper.alpha]);

% Degree of invertibility
disp('Bound estimates: degree of invertibility');
disp([bounds.estim.lower.R2_inv bounds.estim.upper.R2_inv]);
disp('Confidence interval: degree of invertibility');
disp([bounds.ci.lower.R2_inv bounds.ci.upper.R2_inv]);

% Degree of recoverability
disp('Bound estimates: degree of recoverability');
disp([bounds.estim.lower.R2_recov bounds.estim.upper.R2_recov]);
disp('Confidence interval: degree of recoverability');
disp([bounds.ci.lower.R2_recov bounds.ci.upper.R2_recov]);


%% Report FVR bounds

% figure

plots.xticks = 0:12:48; % X axis ticks for FVR plot
plots.titles = {'FVR of Oil Price', 'FVR of Oil Production Growth', 'FVR of Oil Inventories Growth', 'FVR of World IP Growth', 'FVR of U.S. NEER Growth', ...
    'FVR of U.S. IP Growth', 'FVR of U.S. CPI Growth', 'FVR of FFR', 'FVR of VXO', 'FVR of U.S. TOT'};
plots.xlabel = 'Horizon (Months)'; % X axis label for FVR plot
plots.ylabel = ''; % Y axis label for FVR plot
mkdir('figures'); % Figure output folder

for i=1:size(data.Y,2) % For each macro variable...
    
    % Plot bound estimates and CI for identified set
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    if i == 1 || i == 5
    plot_band(bounds.ci.lower.FVR(:,i), bounds.ci.upper.FVR(:,i), bounds.estim.lower.FVR(:,i), bounds.estim.upper.FVR(:,i), ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Estimate of identif. set', sprintf('%d%s', 100*(1-settings_struct.signif_level), '\% conf. interval for identif. set')}, ...
              'YLim', [0 1], 'XLim', [1 max(settings_struct.FVR_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    else
    plot_band(bounds.ci.lower.FVR(:,i), bounds.ci.upper.FVR(:,i), bounds.estim.lower.FVR(:,i), bounds.estim.upper.FVR(:,i), ...
          plots.titles{i}, plots.xlabel, plots.ylabel, [], ...
          'YLim', [0 1], 'XLim', [1 max(settings_struct.FVR_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    end
    grid on;
    drawnow;
    save_fig('figures', strcat('svma_kaenzig_', data.endo_vars{i}));

end

clear i;


%% SVAR-IV analysis for comparison (assumes invertibility)

disp('*** SVAR-IV analysis ***');

% Run analysis
[~, SVARIV_FVD, SVARIV_settings_struct] = SVARIV_estim(data.Y, data.Z, settings{:});

% FVD figure

for i=1:size(data.Y,2) % For each macro variable...
    
    % Plot point estimates and CIs
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    if i == 1 || i == 5
    plot_band(SVARIV_FVD.ci.lower(:,i), SVARIV_FVD.ci.upper(:,i), SVARIV_FVD.estim(:,i), [], ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Point estimate', sprintf('%d%s', 100*(1-SVARIV_settings_struct.signif_level), '\% conf. interval')}, ...
              'YLim', [0 1], 'XLim', [1 max(SVARIV_settings_struct.FVD_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    else
    plot_band(SVARIV_FVD.ci.lower(:,i), SVARIV_FVD.ci.upper(:,i), SVARIV_FVD.estim(:,i), [], ...
          plots.titles{i}, plots.xlabel, plots.ylabel, [], ...
          'YLim', [0 1], 'XLim', [1 max(SVARIV_settings_struct.FVD_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    end
    grid on;
    drawnow;
    save_fig('figures', strcat('svar_kaenzig_', data.endo_vars{i}));
    
end

clear i;
