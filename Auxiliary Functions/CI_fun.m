function [bounds_CI_IS,bounds_CI_para] = CI_fun(bounds_boot,bounds_OLS,settings)

%----------------------------------------------------------------
% Get Inputs
%----------------------------------------------------------------

fields       = settings.fields;
signif_level = settings.signif_level;
optimopts    = settings.optimopts;
fields_param = settings.fields_param;

%----------------------------------------------------------------
% Quantiles of Bootstrap Draws
%----------------------------------------------------------------

for j=1:length(fields)
    bounds_boot_mean.(fields{j}) = squeeze(mean(bounds_boot.(fields{j}),3)); % Average
    bounds_boot_plow.(fields{j}) = squeeze(quantile(bsxfun(@minus, bounds_boot.(fields{j}), bounds_OLS.(fields{j})),signif_level/2,3)); % Lower quantile
    bounds_boot_phigh.(fields{j}) = squeeze(quantile(bsxfun(@minus, bounds_boot.(fields{j}), bounds_OLS.(fields{j})),1-signif_level/2,3)); % Upper quantile
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

for j=1:length(fields_param)
    
    field_LB = sprintf('%s_LB', fields_param{j});
    field_UB = sprintf('%s_UB', fields_param{j});
    bounds_CI_para.lower.(fields_param{j}) = zeros(size(bounds_OLS.(field_LB)));
    bounds_CI_para.upper.(fields_param{j}) = zeros(size(bounds_OLS.(field_LB)));
    
    for l=1:size(bounds_OLS.(field_LB),1)
        for m=1:size(bounds_OLS.(field_LB),2)
            
            % Bootstrap var-cov matrix of estimated lower and upper bounds
            varcov = cov([squeeze(bounds_boot.(field_LB)(l,m,:)) squeeze(bounds_boot.(field_UB)(l,m,:))]);
            
            % adjust for boundary cases in FVR/FVD
            
%             if strcmp(sprintf(fields_param{j}),'FVR') || strcmp(sprintf(fields_param{j}),'FVD')
%                 if bounds_CI_IS.OLS_biascorr.(field_LB)(l,m) <= 0
%                     bounds_CI_IS.OLS_biascorr.(field_LB)(l,m) = bounds_OLS.(field_LB)(l,m);
%                 end
%                 if bounds_CI_IS.OLS_biascorr.(field_UB)(l,m) <= 0
%                     bounds_CI_IS.OLS_biascorr.(field_UB)(l,m) = bounds_OLS.(field_UB)(l,m);
%                 end
%             end    
            if strcmp(sprintf(fields_param{j}),'FVR') || strcmp(sprintf(fields_param{j}),'FVD')
                  bounds_CI_IS.OLS_biascorr.(field_LB)(l,m) = max(0,bounds_CI_IS.OLS_biascorr.(field_LB)(l,m));
                  bounds_CI_IS.OLS_biascorr.(field_UB)(l,m) = max(0,bounds_CI_IS.OLS_biascorr.(field_UB)(l,m));
            end                

            
            % Compute Stoye (2009) confidence interval
            CI = stoye_CI(bounds_CI_IS.OLS_biascorr.(field_LB)(l,m), ...
                              bounds_CI_IS.OLS_biascorr.(field_UB)(l,m), ...
                              varcov, ...
                              signif_level, ...
                              optimopts);
            bounds_CI_para.lower.(fields_param{j})(l,m) = CI(1);
            bounds_CI_para.upper.(fields_param{j})(l,m) = CI(2);
            
        end    
    end
    
end

end