function FVD = FVD_fun(var,hor,yzt_aux,alpha);

Sigma_yzt       = yzt_aux.Sigma_yzt;

maxVar_yt_ythor = yzt_aux.maxVar_yt_ythor(:,:,hor);

FVD = sum(((1/alpha) * Sigma_yzt(1:hor,var)).^2)/(sum(((1/alpha) * Sigma_yzt(1:hor,var,1)).^2) + maxVar_yt_ythor(var,var));

end