%% APPLICATION: SW (2007)
% Mikkel Plagborg-Moller and Christian Wolf
% This version: 05/21/2017

%% HOUSEKEEPING

clc
clear all
close all

addpath('../../Auxiliary Functions')
addpath('../_Subroutines')
shock = 'FG_shock';

%% SETTINGS

%----------------------------------------------------------------
% Model Specification
%----------------------------------------------------------------

settings.set_obsvars      = 1; % set of observables (see below for details)

%----------------------------------------------------------------
% Simulation and VAR Estimation
%----------------------------------------------------------------

settings.VAR_poplaglength        = 250; % population VAR lag length

%----------------------------------------------------------------
% SVAR-IV Identification
%----------------------------------------------------------------

settings.VMA_hor        = 350; % maximal horizon in Wold/structural VMA representation
settings.FVD_hor        = 1:1:11; % horizon for FVD analysis
settings.IRF_hor        = settings.FVD_hor; % horizon for IRF analysis
settings.weight_hor     = 6; % horizon for weight analysis
settings.co             = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]; % color for plots

%----------------------------------------------------------------
% Extra Stuff for Population FVR
%----------------------------------------------------------------

settings.use_KF          = 1; % use Kalman filter for FVR computations?
settings.alpha_ngrid     = 1000; % grid points for lower bound on alpha
settings.bnd_recov       = 0; % naive recoverability-based lower bound on alpha only?
settings.FVR_hor         = settings.FVD_hor; % horizons for FVD analysis
settings.CI_for_R2_inv   = 0; % construct CI for R2_inv?
settings.CI_for_R2_recov = 0; % construct CI for R2_recov?
settings.CI_for_FVR      = 1; % construct CI for FVR?
settings.CI_for_FVD      = 0; % construct CI for FVD?

settings.fields          = {'alpha_LB', 'alpha_UB', 'R2_inv_LB', 'R2_inv_UB', 'R2_recov_LB', 'R2_recov_UB', 'FVR_LB', 'FVR_UB', 'FVD_LB'};

%% VAR REPRESENTATION

%----------------------------------------------------------------
% Solve SW (2007) Model
%----------------------------------------------------------------

% model run

dynare SW_Model noclearall
clean_folder

disp('I have solved and simulated the model.')

disp('Collecting model properties...')

% specify observables

if settings.set_obsvars == 1
    SW_model.obs_y = [5 4 19]; % (r,y,pi)
    SW_model.series = ['Interest Rate';'Real Output  ';'Inflation    '];
elseif settings.set_obsvars == 2
    SW_model.obs_y = [5 4 19 17 18]; % (r,y,pi,c,i)
    SW_model.series = ['Interest Rate';'Real Output  ';'Inflation    ';'Consumption  ';'Investment   '];
elseif settings.set_obsvars == 3
    SW_model.obs_y = [5 4 19 21]; % (r,y,pi,lab)
    SW_model.series = ['Interest Rate';'Real Output  ';'Inflation    ';'Hours        '];
elseif settings.set_obsvars == 4
    SW_model.obs_y = [5 4 19 17 18 20 21]; % (r,y,pi,c,i,w,lab)
    SW_model.series = ['Interest Rate';'Real Output  ';'Inflation    ';'Consumption  ';'Investment   ';'Wage         ';'Hours        '];
end

SW_model.obs_x = [SW_model.obs_y, 34];

% specify shock

SW_model.shock = shock;
clear shock

% get IV parameters

SW_model.alpha    = alpha_ext;
SW_model.sigma_nu = sigma_nu_ext;

% size indicators

SW_model.n_y   = size(SW_model.obs_y,2);
SW_model.n_z   = 1;
SW_model.n_xi  = M_.exo_nbr; % 7 structural shocks + instrument noise
SW_model.n_eps = SW_model.n_xi - 1;
SW_model.n_nu  = 1;
SW_model.n_x   = SW_model.n_y + SW_model.n_z;
SW_model.n_s   = M_.nspred;

% get law of motion for all model variables

SW_model.decision = decision(2:end,:);

% ABCD representations

SW_model.ABCD = ABCD_fun(SW_model);

% delete superfluous variables

clean_workspace

%----------------------------------------------------------------
% Get Population IRFs + FVDs + Shock Sequences + FVRs
%----------------------------------------------------------------

[SW_model.IRF,SW_model.FVD,SW_model.M,SW_model.tot_weights] = pop_analysis(SW_model,settings);

VAR_pop                         = popVAR(SW_model,settings);
[SVARIV_pop.IRF,SVARIV_pop.FVD] = SVARIV_analysis(VAR_pop,SW_model,settings);
yzt_aux                         = get2ndmoments_VAR(VAR_pop,SW_model,settings);
bounds_pop                      = get_IS(yzt_aux,SW_model,settings);
SW_model.FVR                    = bounds_pop.FVR_UB * bounds_pop.alpha_LB^2/SW_model.alpha^2;

disp('...done!')

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

disp('Getting the VAR representation...')

VAR     = popVAR(SW_model,settings);

disp('...done!')

%% ANALYSIS

disp('Doing the SVAR-IV analysis...')

[SVARIV.IRF,SVARIV.FVD,SVARIV.weights] = SVARIV_analysis(VAR,SW_model,settings);

disp('...done!')

%% RESULTS

disp('Plotting results...')

%----------------------------------------------------------------
% IRF
%----------------------------------------------------------------

plotwidth = 0.275;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2*gapsize + 2*plotwidth];
for j = 1:3
    figure(1)
    subplot(1,3,j)
    pos = get(gca, 'Position');
    pos(1) = left_pos(j);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    set(gca,'FontSize',18);
    hold on
    plot(settings.IRF_hor-1,SW_model.IRF(1:settings.IRF_hor(end),j),'linewidth',2,'linestyle','-','color',[0 0 0])
    plot(settings.IRF_hor-1,SVARIV.IRF(1:settings.IRF_hor(end),j),'linewidth',2,'linestyle',':','color',[0 0 0])
    set(gcf,'color','w')
    xlim([0 settings.IRF_hor(end)-1])
    xlabel('Horizon (Quarters)','FontSize',22,'interpreter','latex')
    title(['IRF of ',SW_model.series(j,:)],'fontsize',25,'interpreter','latex')
    if j == 2
        legend({'Truth','SVAR-IV'},'Location','South','fontsize',16,'interpreter','latex')
    end
    grid on
    hold off
end
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.1*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

clear gapsize gapsize_edges j left_pos plotwidth pos

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

plotwidth = 0.275;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2*gapsize + 2*plotwidth];
for j = 1:3
    figure(2)
    subplot(1,3,j)
    pos = get(gca, 'Position');
    pos(1) = left_pos(j);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    set(gca,'FontSize',18);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    plot(settings.FVD_hor-1,SW_model.FVR(1:settings.FVD_hor(end),j),'linewidth',2,'linestyle','-','color',[0 0 0])
    plot(settings.FVD_hor-1,SVARIV.FVD(1:settings.FVD_hor(end),j),'linewidth',2,'linestyle',':','color',[0 0 0])
    set(gcf,'color','w')
    xlim([0 settings.FVD_hor(end)-1])
    ylim([0 0.8])
    xlabel('Horizon (Quarters)','FontSize',20,'interpreter','latex')
    title(['FVR of ',SW_model.series(j,:)],'fontsize',22,'interpreter','latex')
    if j == 2
        legend({'Truth','SVAR-IV'},'Location','North','fontsize',20,'interpreter','latex')
    end
    grid on
    hold off
end
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.1*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('SW_FG_SVARIV_FVR','-deps');
% print('SW_FG_SVARIV_FVR','-dpng');

clear gapsize gapsize_edges j left_pos plotwidth pos

%----------------------------------------------------------------
% Shock Weights
%----------------------------------------------------------------

figure(3)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',18);
hold on
for i_eps = 1:SW_model.n_eps
    plot(0:1:settings.weight_hor-1,SVARIV.weights(1:settings.weight_hor,i_eps),'linewidth',2,'linestyle','-','color',settings.co(i_eps,:))
end
set(gcf,'color','w')
xlim([0 settings.weight_hor-1])
xlabel('Horizon (Quarters)','FontSize',22,'interpreter','latex')
title('SVAR-IV Shock Weights','fontsize',25,'interpreter','latex')
legend('Monetary Policy','Technology','Risk Premium','Government Spending','Investment Efficiency','Price Cost-Push','Wage Cost-Push')
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 1.5*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

clear i_eps pos

disp('...done!')