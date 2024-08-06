% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate dependence of infected density on parameter values
% INPUTS:   none
% METHOD:   Simulate model for varying parameter values and plot
% OUTPUTS:  Curves showing relationships between state variable and parameters


clear all 

figure('Position',[10 10 1000 400]);

    newcolors = [0 0 255 
             60 0 190
             130 0 130
             190 0 60
             255 0 0]/255;
    colororder(newcolors)

  

%% Sensitivity of Infected Density to transmission parameter (m_c)
% Model parameter values
    rS = 6;
    rC = 6;
    rI = 6;
    K = 100;
    mS = 2;
    mI = 10;
    betaS = 2; 
    betaC = [1 2 3 4 5]; 
    chiI = 6;
    u = 1;
    delta = 2;
    mC = [2 2.5 3 3.5 4 4.5 5 5.5 6];
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(mC),length(betaC));
    Yvals = zeros(length(mC),length(betaC));
    Ivals = zeros(length(mC),length(betaC));
    
    for ii=1:length(mC)
        for jj = 1:length(betaC)
            parms = [rS rC rI K mS mC(ii) mI betaS betaC(jj) chiI u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(131) 
    hold on

    for jj = length(betaC):-1:1
        plot(mC,Ivals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(betaC(jj)))) 
    end
    xlim([mC(1) mC(end)])
    xlabel('Compromised Mortality Rate (m_C)')
    ylabel('Infected Density (I^*)')
    title('A');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
    
%% Sensitivity of Infected Density to reproduction parameter (4_c)
% Model parameter values
    rS = 6;
    rC = [0 1 2 3 4 5 6];
    rI = 3;
    K = 100;
    mS = 2;
    mI = 10;
    mC = 4;
    betaS = 2; 
    betaC = [1 10 20 30 40]; 
    chiI = 6;
    u = 1;
    delta = 1;
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(rC),length(betaC));
    Yvals = zeros(length(rC),length(betaC));
    Ivals = zeros(length(rC),length(betaC));
    
    for ii=1:length(rC)
        for jj = 1:length(betaC)
            parms = [rS rC(ii) rI K mS mC mI betaS betaC(jj) chiI u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(132) 
    hold on

%     (255-140)/4/255*(jj-1)+140
%     (194-1)/4/255*(jj-1)+1
%     (158-3)/4/255*(jj-1)+158
    
    for jj = 1:length(betaC)
        plot(rC,Ivals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(betaC(jj))),...
        'color',[(158-3)/4/255*(jj-1)+3/255 (194-1)/4/255*(jj-1)+1/255 (255-140)/4/255*(jj-1)+140/255]) 
    end
    xlim([rC(1) rC(end)])
    xlabel('Compromised Reproduction Rate (r_C)')
    ylabel('Infected Density (I^*)')
    title('B');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')    
    
    
%% Sensitivity of Infected Density to transmission parameter (beta_c)
% Model parameter values
    rS = 6;
    rC = 6;
    rI = 6;
    K = 100;
    mS = 2;
    mI = 4;
    mC = [2 10 20 30 40];
    betaS = 4; 
    betaC = [1 2 3 4 5]; 
    chiI = 3;
    u = 1.5;
    delta = 1.5;
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(mC),length(betaC));
    Yvals = zeros(length(mC),length(betaC));
    Ivals = zeros(length(mC),length(betaC));
    
    for ii=1:length(mC)
        for jj = 1:length(betaC)
            parms = [rS rC rI K mS mC(ii) mI betaS betaC(jj) chiI u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(133) 
    hold on

    for ii = 1:length(mC)
        plot(betaC,Ivals(ii,:),'LineWidth',2,...
        'DisplayName',strcat('m_c=',num2str(mC(ii)))) 
    end
    xlim([betaC(1) betaC(end)])
    xlabel('Compromised Transmission Rate (\beta_C)')
    ylabel('Infected Density (I^*)')
    title('C');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
