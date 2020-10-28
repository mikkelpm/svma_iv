function [bounds_CI_IS,bounds_CI_para] = CI_fun(bounds_boot,bounds_OLS,settings)

%----------------------------------------------------------------
% Get Inputs
%----------------------------------------------------------------

fields       = fieldnames(bounds_OLS);
signif_level = settings.signif_level;
optimopts    = settings.optimopts;

%----------------------------------------------------------------
% Quantiles of Bootstrap Draws
%----------------------------------------------------------------

for j=1:length(fields)
    bounds_boot_mean.(fields{j}) = squeeze(mean(bounds_boot.(fields{j}),3)); % Average
    bounds_boot_plow.(fields{j}) = squeeze(quantile(bounds_boot.(fields{j})-bounds_OLS.(fields{j}),signif_level/2,3)); % Lower quantile
    bounds_boot_phigh.(fields{j}) = squeeze(quantile(bounds_boot.(fields{j})-bounds_OLS.(fields{j}),1-signif_level/2,3)); % Upper quantile
end

%----------------------------------------------------------------
% CI for IS
%----------------------------------------------------------------

for j=1:length(fields)
    bounds_CI_IS.OLS_biascorr.(fields{j}) = 2*bounds_OLS.(fields{j}) - bounds_boot_mean.(fields{j}); % Bias correction
    bounds_CI_IS.lower.(fields{j}) = bounds_OLS.(fields{j}) - bounds_boot_phigh.(fields{j}); % Hall's bootstrap percentile interval, lower bound
    bounds_CI_IS.upper.(fields{j}) = bounds_OLS.(fields{j}) - bounds_boot_plow.(fields{j}); % Hall's bootstrap percentile interval, upper bound
end

%----------------------------------------------------------------
% CI for Parameter
%----------------------------------------------------------------

bounds_CI_para = struct;
if ~settings.CI_para
    return;
end

for j=1:length(fields)
    
    lb_pos = strfind(fields{j},'_LB');
    if isempty(lb_pos)
        continue;
    else
        the_param = extractBefore(fields{j},lb_pos); % Name of parameter
    end
    
    field_LB = fields{j}; % Lower bound field
    field_UB = sprintf('%s_UB', the_param); % Upper bound field
    bounds_CI_para.lower.(the_param) = nan(size(bounds_OLS.(field_LB)));
    bounds_CI_para.upper.(the_param) = nan(size(bounds_OLS.(field_LB)));
    
    for l=1:size(bounds_OLS.(field_LB),1)
        for m=1:size(bounds_OLS.(field_LB),2)
            
            % Bootstrap var-cov matrix of estimated lower and upper bounds
            varcov = cov([squeeze(bounds_boot.(field_LB)(l,m,:)) squeeze(bounds_boot.(field_UB)(l,m,:))]);
            
            % Enforce parameter in [0,1] (except alpha)
            if ~strcmp(the_param,'alpha')
                  bounds_CI_IS.OLS_biascorr.(field_LB)(l,m) = max(0,bounds_CI_IS.OLS_biascorr.(field_LB)(l,m));
                  bounds_CI_IS.OLS_biascorr.(field_UB)(l,m) = max(0,bounds_CI_IS.OLS_biascorr.(field_UB)(l,m));
            end
            
            if ~any(strcmp(field_LB, {'FVD_LB', 'R2_recov_LB'}))
                % Compute Stoye (2009) confidence interval
                CI = stoye_CI(bounds_CI_IS.OLS_biascorr.(field_LB)(l,m), ...
                              bounds_CI_IS.OLS_biascorr.(field_UB)(l,m), ...
                              varcov, ...
                              signif_level, ...
                              optimopts);
                bounds_CI_para.lower.(the_param)(l,m) = CI(1);
                bounds_CI_para.upper.(the_param)(l,m) = CI(2);
            else
                % FVD and R2_recov: one-sided lower confidence interval,
                % since upper bound is always 1
                bounds_CI_para.lower.(the_param)(l,m) = bounds_CI_IS.OLS_biascorr.(field_LB)(l,m)+norminv(signif_level)*sqrt(varcov(1,1));
                bounds_CI_para.upper.(the_param)(l,m) = 1;
            end
            
        end    
    end
    
end

end