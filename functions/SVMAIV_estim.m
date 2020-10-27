function [bounds, id_recov, settings] = SVMAIV_estim(Y, Z, varargin)

    % Inference routines for SVMA-IV analysis
    % Point estimates and bootstrap confidence intervals for identification bounds
    
    % Reference: Mikkel Plagborg-Moller & Christian K. Wolf (2020)
    % "Instrumental Variable Identification of Dynamic Variance Decompositions"
    % https://scholar.princeton.edu/mikkelpm/decomp_iv
    
    
    % Inputs: see below
    
    % Outputs:
    % bounds    struct  Partial identification results:
    %                   - field "estim" contains estimates of bounds (bootstrap bias-corrected)
    %                   - field "ci" contains confidence intervals for identified intervals
    %                   - field "ci_param" contains Stoye (2009) confidence intervals for parameters (if option 'ci_param'=true)
    % id_recov  struct  Point identification results under assumption of recoverability:
    %                   - field "estim" contains parameter estimates (bootstrap bias-corrected)
    %                   - field "ci" contains confidence intervals
    % settings  struct  Settings structure (see below)
    
    % Parameter names in output:
    % alpha     scale parameter
    % R2_inv    degree of invertibility
    % R2_recov  degree of recoverability
    % FVR       Forecast Variance Ratio
    % FVD       Forecast Variance Decomposition
    
    
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
    addParameter(ip, 'compute_R2_inv', true, @islogical);   % bool      Compute degree of invertibility? (default: yes)
    addParameter(ip, 'compute_R2_recov', true, @islogical); % bool      Compute degree of recoverability? (default: yes)
    addParameter(ip, 'compute_FVR', true, @islogical);      % bool      Compute Forecast Variance Ratio? (default: yes)
    addParameter(ip, 'compute_FVD', true, @islogical);      % bool      Compute Forecast Variance Decomposition? (default: yes)
    addParameter(ip, 'horiz', 1:24, @isnumeric);            % 1 x k     Horizons of FVR/FVD to report (default: 1:24)
    addParameter(ip, 'ci_param', false, @islogical);        % bool      Compute confidence intervals for parameters themselves (not identified sets)? (default: no)
    
    % Optional inputs: inference/bootstrap
    addParameter(ip, 'signif', 0.1, @isnumeric);            % 1 x 1     Significance level (default: 10%)
    addParameter(ip, 'n_boot', 1000, @isnumeric);           % 1 x 1     No. of bootstrap repetitions (default: 1000)
    addParameter(ip, 'optim_opts', optimoptions('fmincon', 'Display', 'notify'), @isobject); % obj  Numerical options for Stoye CI construction
    
    % Optional inputs: numerical settings
    addParameter(ip, 'use_kalman', true, @islogical);       % bool      Use Kalman filter for conditional variance calculations? (default: true)
    addParameter(ip, 'VMA_hor', 100, @isnumeric);           % 1 x 1     Truncation horizon for VMA representation of VAR
    
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
    
    settings.CI_for_R2_inv   = ip.Results.compute_R2_inv;   % CI for R2_inv?
    settings.CI_for_R2_recov = ip.Results.compute_R2_recov; % CI for R2_rec?
    settings.CI_for_FVR      = ip.Results.compute_FVR;      % CI for FVR?
    settings.CI_for_FVD      = ip.Results.compute_FVD;      % CI for FVD?
    settings.CI_para         = ip.Results.ci_param;         % Stoye (2009) CI?
    
    settings.FVR_hor        = ip.Results.horiz;         % Horizons for FVR
    settings.FVD_hor        = ip.Results.horiz;         % Horizons for FVD
    
    settings.signif_level   = ip.Results.signif;        % Significance level
    settings.n_boot         = ip.Results.n_boot;        % No. of bootstrap iterations
    
    settings.VMA_hor        = ip.Results.VMA_hor;       % Maximal horizon in Wold/structural VMA representation
    settings.alpha_ngrid    = [];                       % No. of grid points for sharp lower bound on alpha (not used in estimation)
    settings.bnd_recov      = 1;                        % Use weaker/practical lower bound on alpha
    settings.use_KF         = ip.Results.use_kalman;    % Use Kalman filter for computations?
    settings.optimopts      = ip.Results.optim_opts;    % fmincon options for Stoye CI
    
    
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
    VAR_OLS = estimateVAR(dataobj.data.x,settings); 
    disp('...done!');

    % Bootstrap VAR
    disp('Bootstrapping the VAR...');
    VAR_boot = bootstrapVAR(VAR_OLS,dataobj,dataobj.data,settings);
    disp('...done!');
    
    
    %% Bound point estimates
    
    disp('Getting the OLS point estimates of the identified sets...');
    yzt_aux    = get2ndmoments_VAR(VAR_OLS,dataobj,settings);
    bounds_OLS = get_IS(yzt_aux,dataobj,settings);
    disp('...done!');
    
    
    %% Compute bounds for each bootstrap iteration
    
    n_boot = ip.Results.n_boot;
    VAR_sim = VAR_OLS;
    fields = fieldnames(bounds_OLS);
    
    for j=1:length(fields)
        bounds_boot.(fields{j}) = zeros([size(bounds_OLS.(fields{j})) settings.n_boot]);
    end
    
    disp('Mapping each bootstrap draw into objects of interest...');

    for i_boot = 1:n_boot
        
        if mod(i_boot,ceil(n_boot/50)) == 0
            fprintf('%3d%s\n', round(100*i_boot/n_boot), '%');
        end
        
        VAR_sim.VAR_coeff = VAR_boot.VAR_coeff(:,:,i_boot);
        VAR_sim.Sigma_u   = VAR_boot.Sigma_u(:,:,i_boot);

        the_yzt_aux = get2ndmoments_VAR(VAR_sim,dataobj,settings);
        the_bounds = get_IS(the_yzt_aux,dataobj,settings);
        
        for j=1:length(fields)
            bounds_boot.(fields{j})(:,:,i_boot) = the_bounds.(fields{j}); % Store bounds
        end

    end

    disp('...done!');
    
    
    %% Construct CIs
    
    disp('Constructing the confidence intervals...');
    [CI.bounds_CI_IS,CI.bounds_CI_para] = CI_fun(bounds_boot,bounds_OLS,settings);
    disp('...done!');
    
    
    %% Collect results
    
    bounds = struct;
    id_recov = struct;
    params = {'alpha', 'R2_inv', 'R2_recov', 'FVR', 'FVD'};
    
    for ip=1:length(params)
        
        % Parameter name
        the_param = params{ip};
        the_param_LB = sprintf('%s%s',the_param,'_LB');
        the_param_UB = sprintf('%s%s',the_param,'_UB');
        
        % Point estimates of bounds
        bounds.estim.lower.(the_param) = CI.bounds_CI_IS.OLS_biascorr.(the_param_LB);
        bounds.estim.upper.(the_param) = CI.bounds_CI_IS.OLS_biascorr.(the_param_UB);
        
        % Confidence intervals for identified set
        bounds.ci.lower.(the_param) = CI.bounds_CI_IS.lower.(the_param_LB);
        bounds.ci.upper.(the_param) = CI.bounds_CI_IS.upper.(the_param_UB);
        
        % Estimates/CIs under recoverability
        switch the_param
            case 'alpha'
                id_recov.estim.(the_param) = bounds.estim.lower.(the_param);
                id_recov.ci.lower.(the_param) = CI.bounds_CI_IS.lower.(the_param_LB);
                id_recov.ci.upper.(the_param) = CI.bounds_CI_IS.upper.(the_param_LB);
            case 'FVD'
                % Do nothing; FVD is not point-identified even under recoverability
            otherwise
                id_recov.estim.(the_param) = bounds.estim.upper.(the_param);
                id_recov.ci.lower.(the_param) = CI.bounds_CI_IS.lower.(the_param_UB);
                id_recov.ci.upper.(the_param) = CI.bounds_CI_IS.upper.(the_param_UB);
        end
        
    end
    
    bounds.ci_param = CI.bounds_CI_para; % Stoye (2009) confidence intervals for parameters
    
end