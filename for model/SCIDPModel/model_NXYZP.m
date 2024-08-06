% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Simulate frequency form of SCIDP model
% INPUTS:   time (t), state variables (x), and parameters (parms)
% OUTPUTS:  vectors of derivatives (xdot)

function [xdot]=model_NXYZP(t,x,parms)
% State variables
    N = x(1); % Total density
    X = x(2); % Proportion compromised
    Theta = x(3); % Propotion infected
    Z = x(4); % Proportion decimated
    P = x(5); % Spore Density
    W = 1-Theta-X; % Proportion susceptible
    
% parameters = parms 
    rS = parms(1);
    rC = parms(2);
    rI = parms(3);
    rD = parms(4);
    K = parms(5);
    mS = parms(6); 
    mC = parms(7); 
    mI = parms(8); 
    mD = parms(9); 
    betaS = parms(10); 
    betaC = parms(11); 
    chiI = parms(12); 
    chiD = parms(13); 
    u = parms(14); 
    delta = parms(15); 

% Equations
    xdot = zeros(5,1);
    xdot(1) = N*((rS*W+rC*X+rI*(Theta-Z)+rD*Z)*(1-N/K)-mS*W-mC*X-mI*(Theta-Z)-mD*Z);%dN/dt
    xdot(2) = (Theta-Z)*rI*(1-N/K)+Z*rD*(1-N/K)-betaC*X*P-mC*X-X/N*xdot(1);         %dX/dt
    xdot(3) = betaS*W*P+betaC*X*P-mI*(Theta-Z)-mD*Z-Theta/N*xdot(1);                %dTheta/dt
    xdot(4) = betaC*X*P-mD*Z-Z/N*xdot(1);                                           %dZ/dt
    xdot(5) = chiI*(Theta-Z)*N+chiD*Z*N-u*P*N-delta*P;                              %dP/dt
end
