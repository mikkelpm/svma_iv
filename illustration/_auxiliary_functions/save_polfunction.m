function save_polfunction()

% Save decision rule in readable format
% Similar to Dynare's function disp_dr.m when order=1

global M_ oo_ options_;
dr = oo_.dr;

% Preliminaries
if options_.block
    k1 = 1:M_.endo_nbr;
else
    k1 = dr.order_var;
end
[~,ivar] = sort(k1);

% Construct decision rule matrix in correct order
decision = [dr.ys(k1(ivar))'; % constant
            dr.ghx(ivar,:)';  % ghx
            dr.ghu(ivar,:)']; % ghu

% Save
save polfunction decision;

end