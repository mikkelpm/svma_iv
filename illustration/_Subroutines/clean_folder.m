delete sw_*_shock.log
delete sw_*_shock.m
delete sw_*_shock_dynamic.m
delete sw_*_shock_results.mat
delete sw_*_shock_set_auxiliary_variables.m
delete sw_*_shock_static.m
delete(strcat('+sw_', plots.shock, '_shock/*')) % Dynare 4.6+
try
	rmdir(strcat('sw_', plots.shock, '_shock/Output'))
	rmdir(strcat('sw_', plots.shock, '_shock'),'s')
	rmdir(strcat('+sw_', plots.shock, '_shock')) % Dynare 4.6+
catch
end