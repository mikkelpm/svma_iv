delete SW_Model.log
delete SW_Model.m
delete SW_Model_dynamic.m
delete SW_Model_results.mat
delete SW_Model_set_auxiliary_variables.m
delete SW_Model_static.m
rmdir SW_Model/Output
rmdir('SW_Model','s')
load polfunction
data_raw = [ewma epinfma yf y r a b g qs ms spinf sw kpf kp cf invef c inve pinf w lab zcapf rkf kf pkf labf wf rrf mc zcap rk k pk z_ext];
clc