% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate dependence of total density on parameter values
% INPUTS:   none
% METHOD:   Simulate model for varying parameter values and plot
% OUTPUTS:  Curves showing relationships between state variable and parameters

figure('Position',[10 10 1000 400]);

    newcolors = [255 0 0 
                190 0 60
                130 0 130
                60 0 190    
                0 0 255]/255;
    colororder(newcolors)

%% Panel A
    r = 6;
    r_c = 6;
    ri = 3;
    K = 1000;
    m = 2;
    mi = 10;
    beta = 2; 
    beta_c = [1 2 3 4 5]; 
    chi = 6;
    u = 1;
    delta = 2;
    m_c = [2 3 4 5 6 7 8 9 10];
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(m_c),length(beta_c));
    Yvals = zeros(length(m_c),length(beta_c));
    Ivals = zeros(length(m_c),length(beta_c));
    
    for ii=1:length(m_c)
        for jj = 1:length(beta_c)
            parms = [r r_c ri K m m_c(ii) mi beta beta_c(jj) chi u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode23s(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(131) 
    hold on

    newcolors = [255 0 0 
                190 0 60
                130 0 130
                60 0 190    
                0 0 255]/255;
    colororder(newcolors)

    for jj = 1:length(beta_c)
        plot(m_c,Nvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(beta_c(jj)))) 
    end
    xlim([m_c(1) m_c(end)])
    xlabel('Compromised Mortality Rate (m_C)')
    ylabel('Total Density (N^*)')
    title('A');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
%% Panel B
    r = 10;
    r_c = [1 2 3 4 5];
    ri = 4;
    K = 9000/9;
    m = .01;
    mi = 2.25;
    beta = 0.3; 
    beta_c = [0.1 0.3 0.6 0.9 1.2]; 
    chi = 4/9;
    u = .25/9;
    delta = .25;
    m_c = .75;
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(r_c),length(beta_c));
    Yvals = zeros(length(r_c),length(beta_c));
    Ivals = zeros(length(r_c),length(beta_c));
    
    for ii=1:length(r_c)
        for jj = 1:length(beta_c)
            parms = [r r_c(ii) ri K m m_c mi beta beta_c(jj) chi u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode23s(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(132) 
    hold on

    newcolors = [0 0 255
             60 0 190
             130 0 130
             190 0 60
             255 0 0 ]/255;
    colororder(newcolors)

    for jj = length(beta_c):-1:1
        max_b_c = length(beta_c);
        plot(r_c,Nvals(:,max_b_c-jj+1),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(beta_c(max_b_c-jj+1))),...
        'color',[(158-3)/4/255*(jj-1)+3/255 (194-1)/4/255*(jj-1)+1/255 (255-140)/4/255*(jj-1)+140/255]) 
    end
    xlim([r_c(1) r_c(end)])
    xlabel('Compromised Growth Rate (r_C)')
    ylabel('Total Density (N^*)')
    title('B');
    legend('show','Location','southeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
    %% Panel C
    r = 6;
    r_c = [0.1 1 2 3 4];
    ri = 4;
    K = 1000;
    m = 1.5;
    mi = 3;
    beta = 6; 
    beta_c = [6 7 8 9 10 11]; 
    chi = 9.5;
    u = 1;
    delta = 1;
    m_c = 1;
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
  Nvals = zeros(length(beta_c),length(r_c));
    Yvals = zeros(length(beta_c),length(r_c));
    Ivals = zeros(length(beta_c),length(r_c));
    
    for ii=1:length(beta_c)
        for jj = 1:length(r_c)
            parms = [r r_c(jj) ri K m m_c mi beta beta_c(ii) chi u delta];
            x0 = [K,0,0,100]; % initial density values
            [t,x] = ode23s(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
% Plot Equilibrium values    
     
    subplot(133) 
    hold on

    for kk = length(r_c):-1:1
        plot(beta_c,Nvals(:,kk),'LineWidth',2,...
        'DisplayName',strcat('r_c=',num2str(r_c(kk)))) 
    end
    xlim([beta_c(1) beta_c(end)])
    xlabel('Compromised Transmission Rate (\beta_C)')
    ylabel('Total Density (N^*)')
    title('C');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
    
