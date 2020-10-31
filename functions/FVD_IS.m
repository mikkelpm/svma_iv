function FVD_LB = FVD_IS(yzt_aux,model,settings,alpha)

% Identified set for forecast variance decomposition

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_y       = model.n_y;
FVD_hor   = settings.FVD_hor;
FVD_LB    = NaN(size(FVD_hor,2),n_y);
VMA_hor   = settings.VMA_hor;
alpha_UB  = alpha.alpha_UB;
Sigma_yzt = yzt_aux.Sigma_yzt;

Cov_y_shock1_LB   = cell(VMA_hor,1);
Cov_yt_LB         = cell(VMA_hor,1);
for i = 1:VMA_hor
    Cov_y_shock1_LB{i} = 0;
    for j = 1:(VMA_hor-i+1)
        Cov_y_shock1_LB{i} = Cov_y_shock1_LB{i} + (1/alpha_UB)*Sigma_yzt(j,:,1)' * (1/alpha_UB)*Sigma_yzt(j+i-1,:,1);
    end
    Cov_yt_LB{i} = yzt_aux.Cov_y{i} - Cov_y_shock1_LB{i};
end

Sigma_yt_big_LB = NaN(n_y*VMA_hor,n_y*VMA_hor);
for i = 1:VMA_hor
    for j = 1:VMA_hor
        if i > j
            Sigma_yt_big_LB(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_yt_LB{1+abs(i-j)};
        else
            Sigma_yt_big_LB(1+(i-1)*n_y:i*n_y,1+(j-1)*n_y:j*n_y) = Cov_yt_LB{1+abs(i-j)}';
        end
    end
end

yzt_aux.maxVar_yt_ythor = NaN(n_y,n_y,FVD_hor(end));

for hor = 1:FVD_hor(end)
    Sigma_yt_ythor  = Sigma_yt_big_LB;
    Cov_yt          = Cov_yt_LB;
    if hor >= 2
        for j = 2:VMA_hor-(hor-1)
            Sigma_yt_ythor(1:n_y,1+(j-1)*n_y:j*n_y) = Cov_yt{1+abs(1-(j+(hor-1)))}';
        end
        Sigma_yt_ythor(1:n_y,1+(VMA_hor-hor+2-1)*n_y:end) = 0;
        for i = 2:VMA_hor-(hor-1)
            Sigma_yt_ythor(1+(i-1)*n_y:i*n_y,1:n_y) = Cov_yt{1+abs(i+(hor-1)-1)};
        end
        Sigma_yt_ythor(1+(VMA_hor-hor+2-1)*n_y:end,1:n_y) = 0;
    end
    yzt_aux.maxVar_yt_ythor(:,:,hor) = Sigma_yt_ythor(1:n_y,1:n_y) ...
    - Sigma_yt_ythor(1:n_y,n_y+1:end) * Sigma_yt_ythor(n_y+1:end,n_y+1:end)^(-1) * Sigma_yt_ythor(n_y+1:end,1:n_y);
end

%----------------------------------------------------------------
% Lower Bound
%----------------------------------------------------------------

for hor_indx = 1:size(FVD_hor,2)
    hor = FVD_hor(hor_indx);
    for j = 1:n_y
        FVD_LB(hor_indx,j) = FVD_fun(j,hor,yzt_aux,alpha.alpha_UB);
    end
end

end