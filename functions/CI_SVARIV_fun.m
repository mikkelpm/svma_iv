function [IRF_CI,FVD_CI] = CI_SVARIV_fun(SVARIV_OLS,SVARIV_boot,settings)

% Bootstrap confidence intervals, SVAR-IV

%----------------------------------------------------------------
% IRF
%----------------------------------------------------------------

IRF_boot_mean = mean(SVARIV_boot.IRF,3); % Average
IRF_boot_plow = quantile(bsxfun(@minus, SVARIV_boot.IRF, SVARIV_OLS.IRF),settings.signif_level/2,3); % Lower quantile
IRF_boot_phigh = quantile(bsxfun(@minus, SVARIV_boot.IRF, SVARIV_OLS.IRF),1-settings.signif_level/2,3); % Upper quantile

IRF_CI.biascorr = 2*SVARIV_OLS.IRF - IRF_boot_mean; % Bias correction
IRF_CI.lower = SVARIV_OLS.IRF - IRF_boot_phigh; % Hall's bootstrap percentile interval, lower bound
IRF_CI.upper = SVARIV_OLS.IRF - IRF_boot_plow; % Hall's bootstrap percentile interval, upper bound

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

FVD_boot_mean = mean(SVARIV_boot.FVD,3); % Average
FVD_boot_plow = quantile(bsxfun(@minus, SVARIV_boot.FVD, SVARIV_OLS.FVD),settings.signif_level/2,3); % Lower quantile
FVD_boot_phigh = quantile(bsxfun(@minus, SVARIV_boot.FVD, SVARIV_OLS.FVD),1-settings.signif_level/2,3); % Upper quantile

FVD_CI.biascorr = 2*SVARIV_OLS.FVD - FVD_boot_mean; % Bias correction
FVD_CI.lower = SVARIV_OLS.FVD - FVD_boot_phigh; % Hall's bootstrap percentile interval, lower bound
FVD_CI.upper = SVARIV_OLS.FVD - FVD_boot_plow; % Hall's bootstrap percentile interval, upper bound

end
