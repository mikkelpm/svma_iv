function bounds = get_IS(yzt_aux,model,settings)

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

[alpha_LB,alpha_UB,~,~] = alpha_IS(yzt_aux,model,settings);
alpha.alpha_LB = alpha_LB;
alpha.alpha_UB = alpha_UB;
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

[R2_inv_LB,R2_inv_UB,~] = R2_IS(yzt_aux,model,settings,1,alpha);

else
    
R2_inv_LB = [];
R2_inv_UB = [];

end

% recoverability

if settings.CI_for_R2_recov == 1

[R2_recov_LB,R2_recov_UB,~] = R2_IS(yzt_aux,model,settings,round(settings.VMA_hor/2)-1,alpha); % use exactly same bound as for two-sided alpha recoverability computation

else
    
R2_recov_LB = [];
R2_recov_UB = [];

end


%----------------------------------------------------------------
% FVR
%----------------------------------------------------------------

if settings.CI_for_FVR == 1
    
[FVR_LB,FVR_UB,~] = FVR_IS(yzt_aux,model,settings,alpha);

else
    
FVR_LB = [];
FVR_UB = [];

end

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

if settings.CI_for_FVD == 1
    
FVD_LB = FVD_IS(yzt_aux,model,settings,alpha);

else
    
FVD_LB = [];
    
end

%----------------------------------------------------------------
% Package Results
%----------------------------------------------------------------

bounds = struct;
for i=1:length(settings.fields)
    bounds.(settings.fields{i}) = eval(settings.fields{i});
end

end