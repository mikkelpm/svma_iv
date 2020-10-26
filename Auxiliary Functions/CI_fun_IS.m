function bounds_CI_IS = CI_fun_IS(bounds_boot,bounds_OLS,settings)

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

end