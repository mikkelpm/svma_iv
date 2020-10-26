function cross_sd = cross_sd_fun(omega,Sigma_yztilde)

% preparations
hor = size(Sigma_yztilde,1);

% derive lambda_yz_star(omega)
cross_sd = 0;
for l = 1:hor
    cross_sd = cross_sd + Sigma_yztilde(l,:) * exp(-1i * omega * l);
end

cross_sd = reshape(cross_sd,size(cross_sd,2),size(cross_sd,1));

end