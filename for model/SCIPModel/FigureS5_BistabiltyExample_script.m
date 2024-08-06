% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate bistability in NXYP model
% INPUTS:   none
% METHOD:   Simulate model with different initial conditions
% OUTPUTS:  Time series of all state variables

%%% Script for epidemics in DDDT model
clear all 
format long

    newcolors = [0   0.447000000000000   0.741000000000000;
   0.850000000000000   0.325000000000000   0.098000000000000;
   0.929000000000000   0.694000000000000   0.125000000000000;
   0.494000000000000   0.184000000000000   0.556000000000000;
   0.466000000000000   0.674000000000000   0.188000000000000;
   0.301000000000000   0.745000000000000   0.933000000000000;
   0.635000000000000   0.078000000000000   0.184000000000000];


figure('Position',[10 10 1000 400]);

% ORIGINAL Model parameter values
    fac = 1;
    r1 = 0.2*fac;           % Susceptible reproduction rate
    r2 = 0.2*fac;           % Compromised reproduction rate 
    ri = 0.2*fac;          % Infected reproduction rate
    K = 100;            % Carrying Capacity
    m1 = 0.01*fac;          % non-disease mortality
    mi = 0.06*fac;            % disease induced mortality
    mc = 0.06*fac;          % Compromised mortality
    u = 0.05*fac;         % filtering rate
    chi = 8100*fac;          % Shedding rate
    beta1 = 4*10^(-7)*fac;    % Susceptible transmission rate 
    beta2 = 45*10^(-7)*fac;      % Compromised transmission rate
    delta = 0.5*fac;        % Spore degradation rate
    tf = 2000;  
    tf2 = 100;
    tf3 = 20000;  
    times = 0:0.01:tf;

 
% Initial Condition 1: Convergence to disease free equilibrium 
    parms = [r1 r2 ri K m1 mc mi beta1 beta2 chi u delta];
    x0 = [95,0,0,1000]; % initial density values
    [t1,x1] = ode23s(@model_NXYP,times,x0,[],parms);
    
subplot(131) 
    hold on
    plot(t1,x1(:,1)/100,'LineWidth',2,...
        'DisplayName','Total Density') 
    plot(t1,x1(:,2),'LineWidth',2,...
        'DisplayName','Infection Prevalence') 
    plot(t1,x1(:,3),'LineWidth',2,...
        'DisplayName','Proportion Compromised') 
    plot(t1,x1(:,4)/x0(4),'LineWidth',2,...
        'DisplayName','Infectious Propagules') 
    
    xlim([0 20])
    xlabel('Time')
    ylabel('Scaled Variables')
    title('A');
    legend('show','Location','west',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
% Initial Condition 2: Convergence to endemic equilibrium 
    parms = [r1 r2 ri K m1 mc mi beta1 beta2 chi u delta];
    x0 = [95 0 0 1200000];
    times = 0:0.01:tf3;
    [t2,x2] = ode23s(@model_NXYP,times,x0,[],parms);
    
subplot(132) 
    hold on
    plot(t2,x2(:,1)/100,'LineWidth',2,...
        'DisplayName','Total Density') 
    plot(t2,x2(:,2),'LineWidth',2,...
        'DisplayName','Infection Prevalence') 
    plot(t2,x2(:,3),'LineWidth',2,...
        'DisplayName','Proportion Compromised') 
    plot(t2,x2(:,4)/x0(4)*10,'LineWidth',2,...
        'DisplayName','Infectious Propagules') 
    
    xlim([0 25])
    ylim([0 1])
    xlabel('Time')
    ylabel('Scaled Variables')
    title('B');
    legend('show','Location','west',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')

subplot(133) 
    hold on
    plot(t2,x2(:,1)/100,'LineWidth',2,...
        'DisplayName','Total Density') 
    plot(t2,x2(:,2),'LineWidth',2,...
        'DisplayName','Infection Prevalence') 
    plot(t2,x2(:,3),'LineWidth',2,...
        'DisplayName','Proportion Compromised') 
    plot(t2,x2(:,4)/x0(4)*10,'LineWidth',2,...
        'DisplayName','Infectious Propagules') 
    
    xlim([0 tf3])
    xlabel('Time')
    ylim([0 1])
    ylabel('Scaled Variables')
    title('C');
    legend('show','Location','west',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
