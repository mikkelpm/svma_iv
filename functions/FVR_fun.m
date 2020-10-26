function FVR = FVR_fun(var,hor,yzt_aux,alpha);

Sigma_yzt     = yzt_aux.Sigma_yzt;

Var_y_yhor = yzt_aux.Var_y_yhor(:,:,hor);

FVR = sum(((1/alpha) * Sigma_yzt(1:hor,var)).^2)/Var_y_yhor(var,var);

end