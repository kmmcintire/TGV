% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate dependence of infected density on parameter values
% INPUTS:   txt files of eq values and parameters generated by Maple
% METHOD:   A: Simulate model for varying parameter values and plot
%           B: Import data from txt files and plot
% OUTPUTS:  Curves showing relationships between state variable and parameters

clear all 

figure('Position',[10 10 1000 600]);

    newcolors = [0 0 255 
             60 0 190
             130 0 130
             190 0 60
             255 0 0]/255;
    colororder(newcolors)


%% Panel A: Sensitivity to mC
% Model parameter values
    rS = 6;
    rC = 6;
    rI = 6;
    rD = 6;
    K = 100;
    mS = 2;
    mI = 10;
    mD = 10;
    betaS = 2; 
    betaC = 2; 
    chiI = 6;
    chiD = [1 3 5 7 9];
    u = 1;
    delta = 2;
    mC = [2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];
    tf = 500;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(mC),length(chiD));
    Thetavals = zeros(length(mC),length(chiD));
    Itotalvals = zeros(length(mC),length(chiD));
    
    for ii=1:length(mC)
        for jj = 1:length(chiD)
            parms = [rS rC rI rD K mS mC(ii) mI mD betaS betaC chiI chiD(jj) u delta];
            x0 = [K,0,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYZP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Thetavals(ii,jj) = x(end,3);
            Itotalvals(ii,jj) = x(end,1)*x(end,3);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(231) 
    hold on

    for jj = length(chiD):-1:1
        plot(mC,Itotalvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\chi_D=',num2str(chiD(jj)))) 
    end
    xlim([mC(1) mC(end)])
    xlabel('Compromised Mortality Rate (m_C)')
    ylabel('Total Infected Density (I^*+D^*)')
    title('A');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    

%% Panel B: Sensitivity to rC
% Model parameter values
    rS = 7;
    rC = [2 3 4 5 6 7];
    rI = 4;
    rD = 4;
    K = 1000;
    mS = 2;
    mC =3;
    mI = 10;
    mD = 10;
    betaS = 2; 
    betaC = 2; 
    chiI = 6;
    chiD = [0 3 6 9 12];
    u = 1;
    delta = 2;
    tf = 750;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(rC),length(chiD));
    Thetavals = zeros(length(rC),length(chiD));
    Itotalvals = zeros(length(rC),length(chiD));
    
    for ii=1:length(rC)
        for jj = 1:length(chiD)
            parms = [rS rC(ii) rI rD K mS mC mI mD betaS betaC chiI chiD(jj) u delta];
            x0 = [K,0,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYZP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Thetavals(ii,jj) = x(end,3);
            Itotalvals(ii,jj) = x(end,1)*x(end,3);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(232) 
    hold on

    for jj = 1:length(chiD)
        plot(rC,Itotalvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\chi_D=',num2str(chiD(jj))),...
        'color',[(158-3)/4/255*(jj-1)+3/255 (194-1)/4/255*(jj-1)+1/255 (255-140)/4/255*(jj-1)+140/255]) 
    end
    xlim([rC(1) rC(end)])
    xlabel('Compromised Reproduction Rate (r_C)')
    ylabel('Total Infected Density (I^*+D^*)','color','white')
    title('B');
    legend('show','Location','southeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')    
    

%% Panel C: Sensitivity to betaC
% Model parameter values
    rS = 6;
    rC = 6;
    rI = 6;
    rD = 6;
    K = 100;
    mS = 2;
    mC =[2 5 10 15 20];
    mI = 4;
    mD = 4;
    betaS = 4; 
    betaC = 0:1:9; 
    chiI = 3;
    chiD = 3;
    u = 1.5;
    delta = 1.5;
    tf = 500;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(betaC),length(mC));
    Thetavals = zeros(length(betaC),length(mC));
    Itotalvals = zeros(length(betaC),length(mC));
    
    for ii=1:length(betaC)
        for jj = 1:length(mC)
            parms = [rS rC rI rD K mS mC(jj) mI mD betaS betaC(ii) chiI chiD u delta];
            x0 = [K,0,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYZP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Thetavals(ii,jj) = x(end,3);
            Itotalvals(ii,jj) = x(end,1)*x(end,3);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(233) 
    hold on

    for jj = 1:length(mC)
        plot(betaC,Itotalvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('m_C=',num2str(mC(jj)))) 
    end
    xlim([betaC(1) betaC(end)])
    xlabel('Compromised Infection Rate (\beta_C)')
    ylabel('Total Infected Density (I^*+D^*)','color','white')
    title('C');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')    
    
    
%% Panel D: Sensitivity to mD
% Model parameter values
    rS = 4;
    rC = 4;
    rI = 2;
    rD = 2;
    K = 100;
    mS = 1;
    mC =4;
    mI = 2;
    mD = 7:20;
    betaS = 4; 
    betaC = 4; 
    chiI = 3;
    chiD = [1 3 5 7 9];
    u = 1.5;
    delta = 1.5;
    tf = 500;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(mD),length(chiD));
    Thetavals = zeros(length(mD),length(chiD));
    Itotalvals = zeros(length(mD),length(chiD));
    
    for ii=1:length(mD)
        for jj = 1:length(chiD)
            parms = [rS rC rI rD K mS mC mI mD(ii) betaS betaC chiI chiD(jj) u delta];
            x0 = [K,0,0,0,100]; % initial density values
            [t,x] = ode45(@model_NXYZP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Thetavals(ii,jj) = x(end,3);
            Itotalvals(ii,jj) = x(end,1)*x(end,3);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(234) 
    hold on

    for jj = 1:length(chiD)
        plot(mD,Itotalvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\chi_D=',num2str(chiD(jj)))) 
    end
    xlim([mD(1) mD(end)])
    xlabel('Decimated Mortality Rate (m_D)')
    ylabel('Total Infected Density (I^*+D^*)')
    title('D');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')       

%% Panel E: Sensitivity to rD - Imported from Maple
% Model parameter values
    rS = 4;
    rC = 4;
    rI = 2;
    rD = 0:0.25:2;
    K = 100;
    mS = 1;
    mC =4;
    mI = 4;
    mD = 4;
    betaS = 4; 
    betaC = 6; 
    chiI = 3;
    chiD = [1 3 5 7 9];
    u = 1.5;
    delta = 2;
    tf = 500;
    times = 0:.01:tf;
    
    Curve1 = importdata('NumericalEq/dIdrD_chiD1_Isave1.txt');
    Parmsave1 = importdata('NumericalEq/dIdrD_chiD1_Isave1_parm.txt');
    Curve2 = importdata('NumericalEq/dIdrD_chiD3_Isave1.txt');
    Parmsave2 = importdata('NumericalEq/dIdrD_chiD3_Isave1_parm.txt');
    Curve3 = importdata('NumericalEq/dIdrD_chiD5_Isave1.txt');
    Parmsave3 = importdata('NumericalEq/dIdrD_chiD5_Isave1_parm.txt');
    Curve4 = importdata('NumericalEq/dIdrD_chiD7_Isave1.txt');
    Parmsave4 = importdata('NumericalEq/dIdrD_chiD7_Isave1_parm.txt');
    Curve5 = importdata('NumericalEq/dIdrD_chiD9_Isave1.txt');
    Parmsave5 = importdata('NumericalEq/dIdrD_chiD9_Isave1_parm.txt');
    
    Curves = [Curve1 Curve2 Curve3 Curve4 Curve5];
    
% Plot Equilibrium values    
     
    subplot(235) 
    hold on

    for jj = 1:length(chiD)
        plot(Parmsave1,Curves(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\chi_D=',num2str(chiD(jj)))) 
    end
    xlim([rD(1) rD(end)])
    xlabel('Decimated Reproduction Rate (r_D)')
    ylabel('Total Infected Density (I^*+D^*)','color','white')
    title('E');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')       

%% Panel F: Sensitivity to chiD - Imported from Maple
% Model parameter values
    rS = 7;
    rC = [3,4,5,6,7];
    rI = 7;
    rD = 2;
    K = 1000;
    mS = 2;
    mC =3;
    mI = 8;
    mD = 8;
    betaS = 2; 
    betaC = 2; 
    chiI = 6;
    chiD = [2 3 4 5 6 7];
    u = 1;
    delta = 2;
    tf = 500;
    times = 0:.01:tf;

    Curve1 = importdata('NumericalEq/dIdchiD_rS3_Isave1.txt');
    Parmsave1 = importdata('NumericalEq/dIdchiD_rS3_Isave1_parm.txt');
    Curve2 = importdata('NumericalEq/dIdchiD_rS4_Isave1.txt');
    Parmsave2 = importdata('NumericalEq/dIdchiD_rS4_Isave1_parm.txt');
    Curve3 = importdata('NumericalEq/dIdchiD_rS5_Isave1.txt');
    Parmsave3 = importdata('NumericalEq/dIdchiD_rS5_Isave1_parm.txt');
    Curve4 = importdata('NumericalEq/dIdchiD_rS6_Isave1.txt');
    Parmsave4 = importdata('NumericalEq/dIdchiD_rS6_Isave1_parm.txt');
    Curve5 = importdata('NumericalEq/dIdchiD_rS7_Isave1.txt');
    Parmsave5 = importdata('NumericalEq/dIdchiD_rS7_Isave1_parm.txt');
    
    Curves = [Curve1 Curve2 Curve3 Curve4 Curve5];
       
% Plot Equilibrium values    
     
    subplot(236) 
    hold on

    for jj = length(rC):-1:1
        plot(Parmsave1(1:15),Curves(1:15,jj),'LineWidth',2,...
        'DisplayName',strcat('r_C=',num2str(rC(jj)))) 
    end
    for jj = length(rC):-1:1
        plot(Parmsave1(20:end),Curves(20:end,jj),'LineWidth',2,...
        'HandleVisibility','off') 
    end
    xlim([Parmsave1(1) Parmsave1(end)])
    xlabel('Decimated Shedding Rate (\chi_D)')
    ylabel('Total Infected Density (I^*+D^*)','color','white')
    title('F');
    legend('show','Location','southwest',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')      