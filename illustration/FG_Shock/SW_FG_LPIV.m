%% APPLICATION: SW (2007)
% Mikkel Plagborg-Moller and Christian Wolf
% This version: 05/21/2018

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

settings.VAR_poplaglength        = 350; % population VAR lag length
settings.use_KF                  = 1; % use Kalman filter for FVR computations?

%----------------------------------------------------------------
% Identified Set Characterization
%----------------------------------------------------------------

settings.VMA_hor        = 350; % maximal horizon in Wold/structural VMA representation; horizon M for bounds is set as function of that
settings.alpha_ngrid    = 1000; % grid points for lower bound on alpha
settings.bnd_recov      = 0; % naive recoverability-based lower bound on alpha only?
settings.FVR_hor        = [1:1:11]; % horizons for FVR analysis
settings.FVD_hor        = [1:1:11]; % horizons for FVD analysis

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

SW_model.obs_x = [SW_model.obs_y,34];

% specify shock

SW_model.shock = shock;

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
% Get Population IRFs + FVDs + Shock Sequences
%----------------------------------------------------------------

[SW_model.IRF,SW_model.FVD,SW_model.M,SW_model.tot_weights] = pop_analysis(SW_model,settings);

disp('...done!')

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

disp('Getting the VAR representation...')

VAR_pop     = popVAR(SW_model,settings);

disp('...done!')

%% IDENTIFIED SETS

%----------------------------------------------------------------
% Collect Relevant Second-Moment Properties
%----------------------------------------------------------------

disp('Collecting implied second-moment properties...')

yzt_aux = get2ndmoments_VAR(VAR_pop,SW_model,settings);

disp('...done!')

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

disp('Getting the identified set for alpha...')

[alpha.alpha_LB,alpha.alpha_UB,alpha.alpha_true,alpha.alpha_plot] = alpha_IS(yzt_aux,SW_model,settings);

disp('...done!')

%----------------------------------------------------------------
% R^2
%----------------------------------------------------------------

disp('Getting the identified set for R2...')

% invertibility

[R2.R2_inv_LB,R2.R2_inv_UB,R2.R2_inv_true] = R2_IS(yzt_aux,SW_model,settings,1,alpha);

% recoverability

[R2.R2_recov_LB,R2.R2_recov_UB,R2.R2_recov_true] = R2_IS(yzt_aux,SW_model,settings,round(settings.VMA_hor/2)-1,alpha);

disp('...done!')

%----------------------------------------------------------------
% FVR
%----------------------------------------------------------------

disp('Getting the identified set for the FVR...')

[FVR.FVR_LB,FVR.FVR_UB,FVR.FVR_true] = FVR_IS(yzt_aux,SW_model,settings,alpha);

disp('...done!')

%----------------------------------------------------------------
% FVD
%----------------------------------------------------------------

disp('Getting the identified set for the FVD...')

FVD.FVD_LB          = FVD_IS(yzt_aux,SW_model,settings,alpha);

disp('...done!')

%% RESULTS

disp('Plotting results...')

%----------------------------------------------------------------
% Alpha
%----------------------------------------------------------------

figure(1)
hold on
plot(alpha.alpha_plot.omega_grid,alpha.alpha_plot.alpha_LB_vals.^2,'linewidth',2,'linestyle','-','color',[0 0 0])
set(gcf,'color','w')
set(gca,'FontSize',18);
title('Monetary Shock: Spectral Density of 2-Sided Predictor','fontsize',30,'interpreter','latex')
xlabel('Frequency','FontSize',26,'interpreter','latex')
xlim([0 pi])
ylim([0.5 1])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.8*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

%----------------------------------------------------------------
% FVR
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
    hold on
    plot(settings.FVR_hor-1,FVR.FVR_true(:,j),'linewidth',2,'linestyle','-','color',[0 0 0])
    plot(settings.FVR_hor-1,FVR.FVR_UB(:,j),'linewidth',1,'linestyle',':','color',[0 0 0])
    plot(settings.FVR_hor-1,FVR.FVR_LB(:,j),'linewidth',1,'linestyle','--','color',[0 0 0])
    set(gcf,'color','w')
    xlim([0 size(settings.FVR_hor,2)-1])
    limsy=get(gca,'YLim');
    ylim([0 limsy(2)])
    xlabel('Horizon (Quarters)','FontSize',22,'interpreter','latex')
    title(['FVR of ',SW_model.series(j,:)],'fontsize',25,'interpreter','latex')
    if j == 2
        legend({'Truth','Upper Bound','Lower Bound'},'Location','South','fontsize',16,'interpreter','latex')
    end
    grid on
    hold off
end
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 1.5*pos(4)]);
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
    figure(3)
    subplot(1,3,j)
    pos = get(gca, 'Position');
    pos(1) = left_pos(j);
    pos(3) = plotwidth;
    set(gca,'Position', pos)
    set(gca,'FontSize',18);
    hold on
    plot(settings.FVD_hor-1,SW_model.FVD(:,j),'linewidth',2,'linestyle','-','color',[0 0 0])
    plot(settings.FVD_hor-1,FVD.FVD_LB(:,j),'linewidth',1,'linestyle','--','color',[0 0 0])
    set(gcf,'color','w')
    xlim([0 size(settings.FVD_hor,2)-1])
    limsy=get(gca,'YLim');
    ylim([0 limsy(2)])
    xlabel('Horizon (Quarters)','FontSize',22,'interpreter','latex')
    title(['FVD of ',SW_model.series(j,:)],'fontsize',25,'interpreter','latex')
    if j == 2
        legend({'Truth','Lower Bound'},'Location','South','fontsize',16,'interpreter','latex')
    end
    grid on
    hold off
end
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 1.5*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

clear gapsize gapsize_edges j left_pos plotwidth pos

disp('...done!')