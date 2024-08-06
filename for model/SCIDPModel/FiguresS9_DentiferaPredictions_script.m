% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate predictions for D. Dentifera system
% INPUTS:   none
% METHOD:   1: Simulate model with and without TGV
%           2: Plot results
% OUTPUTS:  Time series of all state variables 

clear all 
format long

    newcolors = [0   0.447000000000000   0.741000000000000;
   0.850000000000000   0.325000000000000   0.098000000000000;
   0.929000000000000   0.694000000000000   0.125000000000000;
   0.494000000000000   0.184000000000000   0.556000000000000;
   0.466000000000000   0.674000000000000   0.188000000000000;
   0.301000000000000   0.745000000000000   0.933000000000000;
   0.635000000000000   0.078000000000000   0.184000000000000];
colororder(newcolors)

% Estimated Parameter values
    rs = 0.3;           % Susceptible reproduction rate
    rC = 0.3;           % Compromised reproduction rate 
    rI = 0.3;           % Infected reproduction rate
    rD = 0.3;           % Decimated reproduction rate
    K = 100;            % Carrying Capacity
    mS = 1/30;          % Non-disease mortality rate
    mI = mS;            % Infected induced mortality rate
    mC = 5.66*mS;       % Compromised mortality rate
    mD = mC;            % Decimated mortality rate
    u = 0.0348;         % Uptake rate
    chiI = 666;         % Infected shedding rate
    chiD = 666;         % Decimated shedding rate
    betaS = 0.9*0.1/chiI/(1-exp(-u*2))*u;    % Susceptible transmission rate 
    betaC = betaS;      % Compromised transmission rate
    delta = 10;        % Spore degradation rate
    tf = 200;  
    tf2 = 100;
    times = 0:0.01:tf;
    R0 = betaS*chiI*K/mI/(delta+u*K);
    
% Simulation with TGV
    parms = [rs rC rI rD K mS mC mI mD betaS betaC chiI chiD u delta];
    x0 = [K,0,0,0,100]; % initial density values
    [t,x] = ode45(@model_NXYZP,times,x0,[],parms);
    
    
% Simulation without TGV
    parms_alt = [rs rs rI rI K mS mS mI mI betaS betaS chiI chiI u delta];
    [t_alt,x_alt] = ode45(@model_NXYZP,times,x0,[],parms_alt);
    
% Plot Results   
figure('Position',[10 10 666 840]);

% With Compromised class
    subplot(221) 
    hold on
    plot(t,x(:,1)/K,'-','LineWidth',2,'DisplayName','Total Density')
    plot(t,x(:,3),'-','LineWidth',2,'DisplayName','Infection Prevalence')
    plot(t,x(:,2),'-','LineWidth',2,'DisplayName','Proportion Compromised')
    plot(t,x(:,4),'-','LineWidth',2,'DisplayName','Proportion Decimated')
    plot(t,x(:,5)/max(x(:,5)),'-','LineWidth',2,'DisplayName','Infectious Propagules')
    xlim([0 tf])
    xlabel('Time')
    ylabel('Scaled Variables')
    title('A: With TGV');
    legend('show','Location','west',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'XTick',[0 50 100 150 200])
    
    
% Without Compromised class
    subplot(222) 
    hold on
    plot(t_alt,x_alt(:,1)/K,'--','LineWidth',2,'DisplayName','Total Density')
    plot(t_alt,x_alt(:,3),'--','LineWidth',2,'DisplayName','Infection Prevalence')
    plot(t_alt,x_alt(:,2),'--','LineWidth',2,'DisplayName','Proportion Compromised')
    plot(t_alt,x_alt(:,4),'--','LineWidth',2,'DisplayName','Proportion Decimated')
    plot(t_alt,x_alt(:,5)/max(x_alt(:,5)),'--','LineWidth',2,'DisplayName','Infectious Propagules')
    xlim([0 tf])
    xlabel('Time')
    ylabel('Scaled Variables')
    title('B: Without TGV');
    legend('show','Location','west',...
    'FontSize', 10)
    set(legend,'NumColumns',1) 
    set(gca,'Box','on')
    set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'XTick',[0 50 100 150 200])
    
     
% Comparing Total Density and Infected Density
    subplot(223) 
    hold on
    plot(t,x(:,1)/K,'-','LineWidth',2,'DisplayName','Total Density (TGV)',...
        'color',[0 0.447 0.741])
    plot(t,x(:,3).*x(:,1)/K,'-','LineWidth',2,'DisplayName','Total Infected Density (TGV)',...
        'color', [0.85 0.325 0.098])
    plot(t,x_alt(:,1)/K,'--','LineWidth',2,'DisplayName','Total Density',...
        'color',[0 0.447 0.741])
    plot(t,x_alt(:,3).*x_alt(:,1)/K,'--','LineWidth',2,'DisplayName','Total Infected Density',...
        'color', [0.85 0.325 0.098])
    xlim([0 tf])
    xlabel('Time')
    ylabel('Scaled Variables')
    title('C: Total and Infected Density');
    legend('show','Location','southwest',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'XTick',[0 50 100 150 200])
    
    
% Comparing Infection Prevalence
    subplot(224) 
   hold on
    plot(t,x(:,3),'-','LineWidth',2,'DisplayName','Infection Prevalence (TGV)',...
        'color', [0.85 0.325 0.098])
    plot(t,x_alt(:,3),'--','LineWidth',2,'DisplayName','Infection Prevalence',...
        'color', [0.85 0.325 0.098])
    xlim([0 tf])
    xlabel('Time')
    ylabel('Infection Prevalence')
    title('D: Infection Prevalence');
    legend('show','Location','southwest',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    ylim([0 1])
    set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'XTick',[0 50 100 150 200])
    
     
