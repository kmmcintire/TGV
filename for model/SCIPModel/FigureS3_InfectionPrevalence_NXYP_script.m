% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Illustrate dependence of infection prevalence on parameter values
% INPUTS:   none
% METHOD:   Simulate model for varying parameter values and plot
% OUTPUTS:  Curves showing relationships between state variable and parameters


figure('Position',[10 10 1000 400]);

    newcolors = [0 0 255
             60 0 190
             130 0 130
             190 0 60
             255 0 0 ]/255;
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

    for jj = length(beta_c):-1:1
        plot(m_c,Yvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(beta_c(jj)))) 
    end
    xlim([m_c(1) m_c(end)])
    xlabel('Compromised Mortality Rate (m_C)')
    ylabel('Infection Prevalence (I^*/N^*)')
    title('A');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
%% Panel B
    r = 7;
    r_c = [3 4 5 6 7];
    ri = 4;
    m = 2;
    m_c = 3;
    mi = 9.5;
    K = 1000;
    beta = 2; 
    beta_c = [1 2 3 4 5]; 
    chi = 6;
    u = 1;
    delta = 2;
    
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

    for jj = 1:length(beta_c)
        max_b_c = length(beta_c);
        plot(r_c,Yvals(:,max_b_c-jj+1),'LineWidth',2,...
        'DisplayName',strcat('\beta_c=',num2str(beta_c(max_b_c-jj+1))),...
        'color', newcolors(length(beta_c)-jj+1,:)) 
    end
    xlim([r_c(1) r_c(end)])
    xlabel('Compromised Growth Rate (r_C)')
    ylabel('Infection Prevalence (I^*/N^*)')
    title('B');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    
    %% Panel C
    r = 6;
    r_c = 5;  
    ri = 3;
    m = 4;
    m_c = [4 8 12 16 20];
    mi = 8;
    K = 10;
    beta = 9;
    beta_c = [1 2 3 4 5]; 
    chi = 2;
    u = 0.5;
    delta = 1;
    
    tf = 500;
    tf2 = 1.5;
    times = 0:.01:tf;
    
% Simulate model for each parameter combination
    Nvals = zeros(length(beta_c),length(m_c));
    Yvals = zeros(length(beta_c),length(m_c));
    Ivals = zeros(length(beta_c),length(m_c));
    
    for ii=1:length(beta_c)
        for jj = 1:length(m_c)
            parms = [r r_c ri K m m_c(jj) mi beta beta_c(ii) chi u delta];
            x0 = [100,0,0,100]; % initial density values
            [t,x] = ode23s(@model_NXYP,times,x0,[],parms);
            Nvals(ii,jj) = x(end,1);
            Yvals(ii,jj) = x(end,2);
            Ivals(ii,jj) = x(end,1)*x(end,2);
        end
    end
    
% Plot Equilibrium values    
     
    subplot(133) 
    hold on

    for jj = 1:length(m_c)
        kk = jj;
        plot(beta_c,Yvals(:,jj),'LineWidth',2,...
        'DisplayName',strcat('m_c=',num2str(m_c(jj))),...
        'color',[(158-3)/4/255*(kk-1)+3/255 (194-1)/4/255*(kk-1)+1/255 (255-140)/4/255*(kk-1)+140/255]) 
%         plot(beta_c,Yvals(:,jj),'LineWidth',2,...
%         'DisplayName',strcat('m_c,K=',num2str(m_c(jj)),',',num2str(K(jj))))
    end
    xlim([beta_c(1) beta_c(end)])
    xlabel('Compromised Transmission Rate (\beta_C)')
    ylabel('Infection Prevalence (I^*/N^*)')
    title('C');
    legend('show','Location','northeast',...
    'FontSize', 10)
    set(legend,'NumColumns',1)
    set(gca,'Box','on')
    ylim([.15 0.25])
    
    
