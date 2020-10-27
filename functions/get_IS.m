function bounds = get_IS(yzt_aux,model,settings)

bounds = struct;

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

[bounds.alpha_LB,bounds.alpha_UB,bounds.alpha_true,bounds.alpha_plot] = alpha_IS(yzt_aux,model,settings);
alpha.alpha_LB = bounds.alpha_LB;
alpha.alpha_UB = bounds.alpha_UB;
if isfield(model,'alpha')
    alpha.alpha_true = model.alpha;
else
    alpha.alpha_true = [];
end

%----------------------------------------------------------------
% R2
%----------------------------------------------------------------

% invertibility

if settings.CI_for_R2_inv == 1

[bounds.R2_inv_LB,bounds.R2_inv_UB,bounds.R2_inv_true] = R2_IS(yzt_aux,model,settings,1,alpha);

else
    
bounds.R2_inv_LB = [];
bounds.R2_inv_UB = [];

end

% recoverability

if settings.CI_for_R2_recov == 1

[bounds.R2_recov_LB,bounds.R2_recov_UB,bounds.R2_recov_true] = R2_IS(yzt_aux,model,settings,round(settings.VMA_hor/2)-1,alpha); % use exactly same bound as for two-sided alpha recoverability computation

else
    
bounds.R2_recov_LB = [];
bounds.R2_recov_UB = [];

end


%----------------------------------------------------------------
% FVR
%----------------------------------------------------------------

if settings.CI_for_FVR == 1
    
[bounds.FVR_LB,bounds.FVR_UB,bounds.FVR_true] = FVR_IS(yzt_aux,model,settings,alpha);

else
    
bounds.FVR_LB = [];
bounds.FVR_UB = [];

end

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

if settings.CI_for_FVD == 1
    
bounds.FVD_LB = FVD_IS(yzt_aux,model,settings,alpha);
bounds.FVD_UB = ones(size(bounds.FVD_LB));

else
    
bounds.FVD_LB = [];
bounds.FVD_UB = [];
    
end

end