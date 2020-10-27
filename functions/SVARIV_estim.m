function [IRF, FVD, settings] = SVARIV_estim(Y, Z, varargin)

    % SVAR-IV inference on variance decompositions
    
    
    % Inputs: see below
    
    % Outputs:
    % IRF       struct  Estimation results for absolute impulse responses
    %                   - field "estim": point estimates
    %                   - field "ci": bootstrap confidence intervals
    % FVD       struct  
    % settings  struct  Settings structure (see below)
    
    
    %% Inputs
    
    ip = inputParser;
    
    % Required inputs
    addRequired(ip, 'Y', @isnumeric);           % T x n_y   endogenous variable data matrix
    addRequired(ip, 'Z', @isnumeric);           % T x 1     instrument data vector
    
    % Optional inputs: VAR specification
    addParameter(ip, 'p', [], @isnumeric);      % 1 x 1     VAR lag length, [] means use information criterion (default: [])
    addParameter(ip, 'ic', 'aic', @ischar);     % 1 x 1     information criterion, 'aic' or 'bic' (default: 'aic')
    addParameter(ip, 'ic_max', 24, @isnumeric); % 1 x 1     max lag length for information criterion (default: 24)
    
    % Optional inputs: output requested
    addParameter(ip, 'horiz', 1:24, @isnumeric);            % 1 x k     Horizons of IRF/FVD to report (default: 1:24)
    
    % Optional inputs: inference/bootstrap
    addParameter(ip, 'signif', 0.1, @isnumeric);            % 1 x 1     Significance level (default: 10%)
    addParameter(ip, 'n_boot', 1000, @isnumeric);           % 1 x 1     No. of bootstrap repetitions (default: 1000)
    
    parse(ip, Y, Z, varargin{:}); % Parse inputs
    
    
    %% Create settings structure
    
    settings.select_VAR_simlaglength = isempty(ip.Results.p);   % Use information criterion?
    settings.VAR_simlaglength = ip.Results.p;                   % Lag length (if pre-set)
    settings.max_simlaglength = ip.Results.ic_max;              % Max lag length for information criterion
    switch ip.Results.ic
        case 'bic'
            settings.penalty = @(T) log(T); % BIC
        otherwise
            settings.penalty = @(T) 2/T;    % AIC
    end
    
    settings.IRF_hor        = ip.Results.horiz;         % Horizons for IRF
    settings.FVD_hor        = ip.Results.horiz;         % Horizons for FVD
    settings.VMA_hor        = max(settings.IRF_hor);    % Maximal VMA horizon to compute
    
    settings.signif_level   = ip.Results.signif;        % Significance level
    settings.n_boot         = ip.Results.n_boot;        % No. of bootstrap iterations
    
    
    %% Estimate and bootstrap reduced-form VAR
    
    % Data
    Y = ip.Results.Y;
    Z = ip.Results.Z;
    dataobj.data.x = [Y Z]; % Data matrix
    dataobj.data.y = Y;
    dataobj.data.z = Z;

    % Model dimensions
    dataobj.n_x   = size(dataobj.data.x,2);
    dataobj.n_z   = 1;
    dataobj.n_y   = dataobj.n_x - 1;
    settings.T = size(dataobj.data.x,1);

    % Estimate VAR
    disp('Estimating the VAR...');
    VAR_OLS = estimateVAR_IV(dataobj.data.y,dataobj.data.z,settings); 
    disp('...done!');

    % Bootstrap VAR
    disp('Bootstrapping the VAR...');
    VAR_boot = bootstrapVAR_IV(VAR_OLS,dataobj,dataobj.data,settings);
    disp('...done!');
    
    
    %% FVD point estimates
    
    disp('Getting the OLS point estimates of the FVD...');
    [SVARIV_OLS.IRF,SVARIV_OLS.FVD] = SVARIV_analysis(VAR_OLS,dataobj,settings);
    disp('...done!');
    
    
    %% Compute FVD for each bootstrap iteration
    
    VAR_sim = VAR_OLS;
    
    SVARIV_boot.IRF = NaN(length(settings.IRF_hor),dataobj.n_y,settings.n_boot);
    SVARIV_boot.FVD = NaN(length(settings.FVD_hor),dataobj.n_y,settings.n_boot);

    disp('Mapping each bootstrap draw into objects of interest...')
    fprintf(strcat(repmat('%4d',1,10), '%s\n'), 10:10:100, '%');
    progress_markers = 1/40:1/40:1;

    for i_boot = 1:settings.n_boot
        
        if sum(i_boot/settings.n_boot>=progress_markers)>sum((i_boot-1)/settings.n_boot>=progress_markers)
            fprintf('x');
        end
        
        VAR_sim.VAR_coeff_y = VAR_boot.VAR_coeff_y(:,:,i_boot);
        VAR_sim.Sigma_u_y   = VAR_boot.Sigma_u_y(:,:,i_boot);
        VAR_sim.gamma       = VAR_boot.gamma(:,i_boot);

        [SVARIV_boot.IRF(:,:,i_boot),SVARIV_boot.FVD(:,:,i_boot)] = SVARIV_analysis(VAR_sim,dataobj,settings);

    end

    disp(' ');
    disp('...done!')
    
    
    %% Construct CIs
    
    disp('Constructing the confidence intervals...');
    [IRF_CI,FVD_CI] = CI_SVARIV_fun(SVARIV_OLS,SVARIV_boot,settings);
    disp('...done!');
    
    
    %% Collect results
    
    IRF = struct;
    FVD = struct;
    
    IRF.estim = IRF_CI.biascorr;
    FVD.estim = FVD_CI.biascorr;
    
    IRF.ci = rmfield(IRF_CI, 'biascorr');
    FVD.ci = rmfield(FVD_CI, 'biascorr');
    
end