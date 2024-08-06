% AUTHOR:   Michael Cortez, mcortez@fsu.esu
% DATE:     February 2023
% PURPOSE:  Model to simulate frequency dependent form of SCIP model
% INPUTS:   time (t), state variables (x), and parameters (parms)
% OUTPUTS:  vector of derivatives

function [xdot]=model_NXYP(t,x,parms)
% A function to simulate the NYWZ Model with compromised individuals
% State variables
    N = x(1); % Total density
    Y = x(2); % Propotion infected
    X = x(3); % Proportion compromised
    P = x(4); % Spore Density
    W = 1-Y-X; % Proportion susceptible
    
% parameters = parms 
    rS = parms(1);
    rC = parms(2);
    rI = parms(3);
    K = parms(4);
    mS = parms(5); 
    mC = parms(6); 
    mI = parms(7); 
    betaS = parms(8); 
    betaC = parms(9); 
    chiI = parms(10); 
    u = parms(11); 
    delta = parms(12); 

% Equations
    xdot = zeros(4,1);
    xdot(1) = N*((rS*W+rC*X+rI*Y)*(1-N/K)-mS*W-mC*X-mI*Y);  %dN/dt
    xdot(2) = betaS*W*P+betaC*X*P-mI*Y-Y/N*xdot(1);         %dY/dt
    xdot(3) = Y*rI*(1-N/K)-betaC*X*P-mC*X-X/N*xdot(1);      %dC/dt
    xdot(4) = chiI*Y*N-u*P*N-delta*P;                       %dP/dt
end
