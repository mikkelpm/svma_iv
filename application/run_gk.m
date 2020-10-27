clear all;

% Empirical application: Gertler & Karadi (AEJ Macro, 2015)


%% Load data

% Settings
data.file = 'gk2015.csv';                               % Data file
data.sample = datetime({'1990-1-1','2012-06-01'});      % Sample end points
data.endo_vars = {'ff', 'dlogip', 'dlogcpi', 'ebp'};    % Endogenous (y) variables
data.iv_var = 'ff4_tc';                                 % External IV (z)

% Load
data.table = readtable(data.file); % Read from file
data.sample_bool = datetime(data.table.date) >= data.sample(1) & datetime(data.table.date) <= data.sample(2); % Sample marker

% Select variables
data.Y = data.table{data.sample_bool, data.endo_vars};  % Endogenous variable data
data.Z = data.table{data.sample_bool, data.iv_var};     % External IV data


%% SVMA-IV inference

disp('*** SVMA-IV analysis ***');

% Preliminaries
addpath('../functions');        % Add folder with SVMA-IV analysis functions
rng(2018);                      % Seed random number generator (for bootstrap)

% Estimation settings (see other optional settings in "functions/SVMAIV_estim.m")
settings = {'ic', 'aic';        % Information criterion
            'n_boot', 500;      % Number of bootstrap samples
            'signif', 0.1;      % Significance level
            'horiz', 1:24}';    % Horizons of FVR to report

% Run inference routines
[bounds, id_recov, settings_struct] = SVMAIV_estim(data.Y, data.Z, settings{:});


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


%% Figure of FVR bounds

% Settings
plots.xticks = 3:3:24; % X axis ticks for FVR plot
plots.titles = {'FVR of Federal Funds Rate', 'FVR of Industrial Production Growth', 'FVR of CPI Growth', 'FVR of Excess Bond Premium'};
plots.xlabel = 'Horizon (Months)'; % X axis label for FVR plot
plots.ylabel = ''; % Y axis label for FVR plot

for i=1:size(data.Y,2) % For each macro variable...
    
    % Plot bound estimates and CI for identified set
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    plot_band(bounds.ci.lower.FVR(:,i), bounds.ci.upper.FVR(:,i), bounds.estim.lower.FVR(:,i), bounds.estim.upper.FVR(:,i), ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Estimate of identif. set', sprintf('%d%s', 100*(1-settings_struct.signif_level), '\% conf. interval for identif. set')}, ...
              'YLim', [0 1], 'XLim', [1 max(settings_struct.FVR_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    grid on;

end

drawnow;


%% SVAR-IV analysis for comparison (assumes invertibility)

disp('*** SVAR-IV analysis ***');

% Run analysis
[~, SVARIV_FVD, SVARIV_settings_struct] = SVARIV_estim(data.Y, data.Z, settings{:});

% Plot FVD
for i=1:size(data.Y,2) % For each macro variable...
    
    % Plot point estimates and CIs
    figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);
    plot_band(SVARIV_FVD.ci.lower(:,i), SVARIV_FVD.ci.upper(:,i), SVARIV_FVD.estim(:,i), [], ...
              plots.titles{i}, plots.xlabel, plots.ylabel, {'Point estimate', sprintf('%d%s', 100*(1-SVARIV_settings_struct.signif_level), '\% conf. interval')}, ...
              'YLim', [0 1], 'XLim', [1 max(SVARIV_settings_struct.FVD_hor)], 'XTick', plots.xticks, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.2);
    grid on;
    
end

