# AUTHOR:   Michael Cortez, mcortez@fsu.esu
# DATE:     February 2022
# PURPOSE:  Simulate frequency-form of model with transgenerational virulence
# INPUTS:   None
# METHOD:   Step 1: Define parameter values and model
#           Step 2: Simulate model with and without transgenerational virulence
#           Step 3: Plot simulated time series
# OUTPUTS:  Time series of total density, infected density, infection prevalence

rm(list=ls())
library(deSolve)

#########################################################################
### STEP 1: DEFINE MODEL AND PARAMETER VALUES

## Frequency-form of model
NXYZP_Model=function(t,V,parms){
  
  #parameters
  rS=parms[1]
  rC=parms[2]
  rI=parms[3]
  rD=parms[4]
  K=parms[5]
  mS=parms[6]
  mC=parms[7]
  mI=parms[8]
  mD=parms[9]
  u=parms[10]
  chiI=parms[11]
  chiD=parms[12]
  betaS=parms[13]
  betaC=parms[14]
  delta=parms[15]
  
  #variables
  N=V[1]
  X=V[2]
  Theta=V[3]
  Z=V[4]
  P=V[5]
  W=1-X-Theta
  
  #equations
  dN = N*((rS*W+rC*X+rI*(Theta-Z)+rD*Z)*(1-N/K)-mS*W-mC*X-mI*(Theta-Z)-mD*Z)
  dX = (Theta-Z)*rI*(1-N/K) + Z*rD*(1-N/K)-betaC*X*P-mC*X-X/N*dN
  dTheta = betaS*W*P+betaC*X*P-(Theta-Z)*mI-Z*mD-Theta/N*dN 
  dZ = betaC*X*P-mD*Z - Z/N*dN
  dP = chiI*(Theta-Z)*N+chiD*Z*N-u*P*N-delta*P
  dV = c(dN,dX,dTheta,dZ,dP)
  
  #output
  return(list(dV))
}

## Parameters 
  rS = 0.3;             # Susceptible reproduction rate
  rC = 0.3;             # Compromised reproduction rate 
  rI = 0.3;             # Infected reproduction rate
  rD = 0.3;             # Decimated reproduction rate
  K = 100;              # Carrying Capacity
  mS = 1/30;            # non-disease mortality
  mC = 5.66*mS;         # Compromised mortality
  mI = mS;              # disease induced mortality
  mD = 5.66*mS;         # Decimated mortality
  u = 0.0348;           # filtering rate
  chiI = 666;           # Shedding rate
  chiD = 666;           # Shedding rate
  betaS = 0.9*0.1/chiI/(1-exp(-u*2))*u;    # Susceptible transmission rate 
  betaC = betaS;        # Compromised transmission rate
  delta = 10;   

  parms_yesTGV = c(rS,rC,rI,rD,K,mS,mC,mI,mD,u,chiI,chiD,betaS,betaC,delta) # with TGV
  parms_noTGV = c(rS,rS,rI,rI,K,mS,mS,mI,mI,u,chiI,chiI,betaS,betaS,delta)# w/out TGV

#########################################################################
### STEP 2: SIMULATE MODEL
  init_conds=c(K,0,0,0,100)
  times=seq(0,200,by=0.01)
  
  # Simulate model with transgenerational virulence
  out_yesTGV=ode(init_conds,times,NXYZP_Model,parms_yesTGV,method="lsoda")
  
  # Simulate model without transgenerational virulence
  out_noTGV=ode(init_conds,times,NXYZP_Model,parms_noTGV,method="lsoda")
  

#########################################################################
### STEP 3: PLOT TIME SERIES
  tiff(file="plot.tiff",
       width=8, height=3, units="in", res=400)
  
  par(mfrow=c(1,3))
  
  # Total Host Density
  matplot(out_yesTGV[,1],out_yesTGV[,2],type='l',lty="solid",col="#0583D2",lwd=3,
          xlab="Time (days)",ylab="Total Density (indv./L)",main="Total Host Density",
          font=2,font.lab=2,cex.lab = 1.2,
          xlim=c(0,100),ylim=c(0,100))
  lines(out_noTGV[,1],out_noTGV[,2],type='l',lty="solid",col="black",lwd=2)
  legend("center", legend=c("With TGV","Without TGV"),col=c("#0583D2","black"), lty=1, lwd=c(3,2), cex=1)
  
  # Infected Density
  matplot(out_yesTGV[,1],out_yesTGV[,4]*out_yesTGV[,2],type='l',lty="solid",col="#0583D2",lwd=3,
          xlab="Time (days)",ylab="Infected Density (indv./L)",main="Infected Density",
          font=2,font.lab=2,cex.lab = 1.2,
          xlim=c(0,100),ylim=c(0,100))
  lines(out_noTGV[,1],out_noTGV[,4]*out_noTGV[,2],type='l',lty="solid",col="black",lwd=2)
  
  # Infection Prevalence
  matplot(out_yesTGV[,1],out_yesTGV[,4],type='l',lty="solid",col="#0583D2",lwd=3,
          xlab="Time (days)",ylab="Infection Prevalence",main="Infection Prevalence",
          font=2,font.lab=2,cex.lab = 1.2,
          xlim=c(0,100),ylim=c(0,1))
  lines(out_noTGV[,1],out_noTGV[,4],type='l',lty="solid",col="black",lwd=2)
  


  dev.off()  
  
  
  
  