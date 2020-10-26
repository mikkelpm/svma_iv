function [SVARIV_IRF,SVARIV_FVD,SVARIV_weights] = SVARIV_analysis(VAR,model,settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_y         = model.n_y;
VMA_hor     = settings.VMA_hor;
laglength   = VAR.laglength;
VAR_coeff_y = VAR.VAR_coeff_y;
gamma       = VAR.gamma;
Sigma_u_y   = VAR.Sigma_u_y;
FVD_hor     = settings.FVD_hor;

%----------------------------------------------------------------
% IRFs
%----------------------------------------------------------------

% compute IRFs

IRF_Wold = zeros(n_y,n_y,VMA_hor);
IRF_Wold(:,:,1) = eye(n_y); % Reduced-form impact impulse response matrix
SVARIV_IRF = zeros(n_y,1,VMA_hor);

for l = 1:VMA_hor
    
    if l<VMA_hor
        for j=1:min(l,laglength)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + VAR_coeff_y(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
    SVARIV_IRF(:,:,l) = IRF_Wold(:,:,l)*gamma;
    
end

SVARIV_IRF = squeeze(SVARIV_IRF)';

SVARIV_IRF = SVARIV_IRF(FVD_hor,:);

%----------------------------------------------------------------
% FVDs
%----------------------------------------------------------------

SVARIV_FVD = zeros(VMA_hor,n_y);
FVD_numer = zeros(1,n_y);
FVD_denom = zeros(1,n_y);

for l=1:VMA_hor % For each forecast horizon...
    
    % Add to numerator of FVD
    FVD_numer = FVD_numer + ((IRF_Wold(:,:,l)*gamma).^2)';
    
    % Add to denominator of FVD
    forec_var_contrib = IRF_Wold(:,:,l)*Sigma_u_y*IRF_Wold(:,:,l)';
    FVD_denom = FVD_denom + diag(forec_var_contrib)';
    
    % FVD at horizon l
    SVARIV_FVD(l,:) = FVD_numer./FVD_denom;
    
end

SVARIV_FVD = SVARIV_FVD(FVD_hor,:);

%----------------------------------------------------------------
% Weights
%----------------------------------------------------------------

if isfield(model,'M')

n_xi       = model.n_xi;
SVARIV_weights = NaN(VMA_hor,n_xi);

for l = 1:VMA_hor
    SVARIV_weights(l,:) = gamma' * Sigma_u_y^(-1) * model.M(:,:,l);
end

else
    
SVARIV_weights = [];

end

end