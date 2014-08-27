
# File structure
File = "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2014 -- Paulik diagram/"

# Load libraries
library(R2jags)
library(abind)

# Settings
Nyears = 30
Nstages = 6
SigmaP = rep(0.2, Nstages-1)
SigmaM = rep(0.5, Nstages-1)
Alpha = c(0.5,0.5,0.5,0.5,0.5)
Beta = c(0.0,0.0,0.1,0.1,0.3)

# Generate data
Results = array(NA, dim=c(2,Nstages,Nyears), dimnames=list(c("True","Obs"),paste("Stage",1:Nstages),paste("Year",1:Nyears)))
Results['True',1,] = runif(Nyears,min=0,max=1000)
  Results['Obs',1,] = Results[1,1,] * exp(rnorm(Nyears, mean=0, sd=SigmaM[1]))
for(StageI in 1:(Nstages-1)){
  Results['True',StageI+1,] = Results['True',StageI,] * exp(-Alpha[StageI]) * exp(-Beta[StageI]*Results['True',StageI,]) * exp(rnorm(Nyears, mean=0, sd=SigmaP[StageI]))
  Results['Obs',StageI+1,] = Results['True',StageI+1,] * exp(rnorm(Nyears, mean=0, sd=SigmaM[StageI]))
}

# or load real data
if(FALSE){
  CSV = read.csv( file=paste(File,"Indices 2014_v2.csv",sep=""), header=TRUE)
  Results = t(CSV[,-c(1:2)])
  Results = Results[,which(apply(Results,MARGIN=2,FUN=function(Vec){any(!is.na(Vec))}))]
  Results = aperm(abind( array(NA, dim=dim(Results)), Results, along=3),c(3,1,2))
  dimnames(Results) = list(c("True","Obs"), paste("Stage",1:dim(Results)[2]), paste("Year",1:dim(Results)[3]))
  Nstages = dim(Results)[2]
  Nyears = dim(Results)[3]
}

#############################
# Single-stage model
#  Only uses spawning biomass and recruits
#############################

# Plot data -- as single-state model
Which = which( !is.na(Results['Obs',1,]) & !is.na(Results['Obs',6,]) )
Col = colorRampPalette(colors=c("blue","purple","red"))
par(mfrow=c(1,1), mar=c(3,3,2,0))
  plot( x=Results['Obs',1,], y=Results['Obs',6,], type="b", pch=20, col=Col(length(Which)))
  #points( x=Results['True',StageI,], y=Results['True',StageI+1,], col="red")

# Define Ricker model for first to last quadrant
Ricker_Fn = function(){
  # Priors
  Alpha ~ dunif(-20,20)
  Beta ~ dunif(-20,20)
  SigmaM ~ dunif(0,10)
  TauM <- pow( SigmaM, -2)
  # Projection and likelihood
  for(t in 1:Nyears){
    R_exp[t] <- B_obs[t] * exp(-Alpha) * exp(-Beta*B_obs[t])
    R_obs[t] ~ dlnorm( log(R_exp[t]), TauM )
  }
  # Predictive distribution
  for(p in 1:Npred){
    Rpred_exp[p] <- Bpred[p] * exp(-Alpha) * exp(-Beta*Bpred[p])
    Rpred_hat[p] ~ dlnorm( log(Rpred_exp[p]), TauM )
  }
}
# Define Ricker model for first to last quadrant -- Using multiple regimes
RickerRegime_Fn = function(){
  # Priors
  for(reg in 1:Nregime){
    Alpha[reg] ~ dunif(-20,20)
  }
  Beta ~ dunif(-20,20)
  SigmaM ~ dunif(0,10)
  TauM <- pow( SigmaM, -2)
  # Transition matrix
  Regime[1] <- 1
  for(reg in 1:Nregime){
    Prior[reg] <- 1/Nregime
  }
  for(reg in 1:Nregime){
    Trans_Mat[reg,1:Nregime] ~ ddirich( Prior )
  }
  for(t in 2:Nyears){
    Regime[t] ~ dcat( Trans_Mat[Regime[t-1],1:Nregime] )
  }
  # Projection and likelihood
  for(t in 1:Nyears){
    R_exp[t] <- B_obs[t] * exp(-Alpha[Regime[t]]) * exp(-Beta*B_obs[t])
    R_obs[t] ~ dlnorm( log(R_exp[t]), TauM )
  }
  # Predictive distribution
  for(reg in 1:Nregime){
    for(p in 1:Npred){
      Rpred_exp[reg,p] <- Bpred[p] * exp(-Alpha[reg]) * exp(-Beta*Bpred[p])
      Rpred_hat[reg,p] ~ dlnorm( log(Rpred_exp[reg,p]), TauM )
    }
  }
}

# MCMC settings
Nsim = Nburnin = 5e3
  Nthin = 5e0
Which = which( !is.na(Results['Obs',1,]) & !is.na(Results['Obs',6,]) )

# Run jags for 1-regime 
Params = c("Alpha","Beta","SigmaM","Npred_hat","Npred_exp","Rpred_hat")
Data = list(Nyears=length(Which), B_obs=Results['Obs',1,Which], R_obs=Results['Obs',6,Which], Bpred=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=100), Npred=100 )
Ricker <- jags(model.file=Ricker_Fn, working.directory=NULL, data=Data, parameters.to.save=Params, n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
Ricker$BUGSoutput$summary

# Run jags for 2-regime 
Params = c("Alpha","Beta","SigmaM","Npred_hat","Npred_exp","Rpred_hat","Trans_Mat","Regime")
Data = list(Nregime=2, Nyears=length(Which), B_obs=Results['Obs',1,Which], R_obs=Results['Obs',6,Which], Bpred=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=100), Npred=100 )
RickerRegime <- jags(model.file=RickerRegime_Fn, working.directory=NULL, data=Data, parameters.to.save=Params, n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
Ricker$BUGSoutput$summary

# Predictive distribution
png( file=paste(File,"SRcurve_Ricker.png",sep=""), width=6, height=6, res=200, units="in")
  Pred = apply( Ricker$BUGSoutput$sims.list$Rpred_hat, MARGIN=2, FUN=mean)
  plot( x=Data$Bpred, y=Pred, type="n", ylim=c(0,max(Results['Obs',6,],quantile(Ricker$BUGSoutput$sims.list$Rpred_hat,prob=0.99),na.rm=TRUE)) )
  for(p in 1:Data$Npred){
    lines( x=rep(Data$Bpred[p],2), y=quantile(Ricker$BUGSoutput$sims.list$Rpred_hat[,p], prob=c(0.1,0.9)) )
    points( x=Data$Bpred[p], y=mean(Ricker$BUGSoutput$sims.list$Rpred_hat[,p]) )
  }
  points( x=Results['Obs',1,], y=Results['Obs',6,], pch=20, type="b", col=Col(length(Which)))
dev.off()

#############################
# State-space Paulik approach
#############################

# Plot data -- as multi-state model
Which = which( !is.na(Results['Obs',1,]) & !is.na(Results['Obs',6,]) )
Col = colorRampPalette(colors=c("blue","purple","red"))
  Ncol=ceiling(sqrt(Nstages)); Nrow=ceiling(Nstages/Ncol)
par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0))
for(StageI in 1:(Nstages-1)){
  plot( x=Results['Obs',StageI,], y=log(Results['Obs',StageI+1,]/Results['Obs',StageI,]), col=Col(length(Which)), type="p" )
  text(x=Results['Obs',StageI,], y=log(Results['Obs',StageI+1,]/Results['Obs',StageI,]), labels=1:dim(Results)[3], col="red")
}

# Define simple Paulik diagram
Paulik_Fn = function(){
  # Priors
  #Alpha_p ~ dunif(0,20)
  #Beta_p ~ dunif(0,20)
  #SigmaP_p ~ dunif(0,10)
  SigmaM_p ~ dunif(0,10)
  for(s in 1:(Nstages-1)){
    #Alpha[s] <- Alpha_p
    #Beta[s] <- Beta_p
    Alpha[s] ~ dunif(0,20)
    Beta[s] ~ dunif(-20,20)
    SigmaP[s] ~ dunif(0,20)
    TauP[s] <- pow( SigmaP[s], -2)
  }
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  }
  # Projection and likelihood
  for(t in 1:Nyears){
    N_hat[1,t] ~ dlnorm(0,0.000001)
    N_obs[1,t] ~ dlnorm( log(N_hat[1,t]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,t] <- N_hat[s-1,t] * Alpha[s-1] * exp(-Beta[s-1]*N_hat[s-1,t])
      N_hat[s,t] ~ dlnorm( log(N_exp[s,t]), TauP[s-1] )
      N_obs[s,t] ~ dlnorm( log(N_hat[s,t]), TauM[s] )
    }
  }
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[s-1,p] * Alpha[s-1] * exp(-Beta[s-1]*Npred_hat[s-1,p])
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[s-1] )
    }
  }
}
# Define simple Paulik diagram
PaulikRegime_Fn = function(){
  # Priors
  #Alpha_p ~ dunif(0,20)
  #Beta_p ~ dunif(0,20)
  #SigmaP_p ~ dunif(0,10)
  SigmaM_p ~ dunif(0,10)
  for(s in 1:(Nstages-1)){
    #Alpha[s] <- Alpha_p
    #Beta[s] <- Beta_p
    Alpha[1,s] ~ dunif(-10,10)
    Beta[s] ~ dunif(-20,20)
    SigmaP[s] ~ dunif(0,20)
    TauP[s] <- pow( SigmaP[s], -2)
  }
  Alpha[2,1] <- Alpha[1,1]
  Alpha[2,2] ~ dunif(-10,10)
  Alpha[2,3] <- Alpha[1,3]
  Alpha[2,4] <- Alpha[1,4]
  Alpha[2,5] <- Alpha[1,5]
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  }
  # Transition matrix
  Regime[1] <- 1
  for(reg in 1:Nregime){
    Prior[reg] <- 1/Nregime
  }
  for(reg in 1:Nregime){
    Trans_Mat[reg,1:Nregime] ~ ddirich( Prior )
  }
  for(t in 2:Nyears){
    Regime[t] ~ dcat( Trans_Mat[Regime[t-1],1:Nregime] )
  }
  # Projection and likelihood
  for(t in 1:Nyears){
    N_hat[1,t] ~ dlnorm(0,0.000001)
    N_obs[1,t] ~ dlnorm( log(N_hat[1,t]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,t] <- N_hat[s-1,t] * exp(-Alpha[Regime[t],s-1]) * exp(-Beta[s-1]*N_hat[s-1,t])
      N_hat[s,t] ~ dlnorm( log(N_exp[s,t]), TauP[s-1] )
      N_obs[s,t] ~ dlnorm( log(N_hat[s,t]), TauM[s] )
    }
  }
  # Predictive distribution
  for(reg in 1:Nregime){
    for(p in 1:Npred){
      Npred_hat[reg,1,p] <- Npred_init[p]
      for(s in 2:Nstages){
        Npred_exp[reg,s,p] <- Npred_hat[reg,s-1,p] * exp(-Alpha[reg,s-1]) * exp(-Beta[s-1]*Npred_hat[reg,s-1,p])
        Npred_hat[reg,s,p] ~ dlnorm( log(Npred_exp[reg,s,p]), TauP[s-1] )
      }
    }
  }  
}

# MCMC settings
Nsim = Nburnin = 1e4
  Nthin = 1e1

# Run jags  -- Single regime
Data = list(N_obs=Results['Obs',,], Nstages=Nstages, Nyears=Nyears, Npred_init=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=100), Npred=100 )
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data, parameters.to.save=c("Alpha","Beta","SigmaM_p","SigmaP","Npred_hat"), n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
Paulik$BUGSoutput$summary

# Run jags -- 2-regime
Data = list(Nregime=2, N_obs=Results['Obs',,], Nstages=Nstages, Nyears=Nyears, Npred_init=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=100), Npred=100 )
PaulikRegime <- jags(model.file=PaulikRegime_Fn, working.directory=NULL, data=Data, parameters.to.save=c("Alpha","Beta","SigmaM_p","SigmaP","Trans_Mat","Regime"), n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates                                                                                                                         # "Npred_hat",
PaulikRegime$BUGSoutput$summary

### Plot estimates vs. truth
if(FALSE){
  dev.new()
  Ncol=ceiling(sqrt(Nyears)); Nrow=ceiling(Nyears/Ncol)
  par(mfrow=c(Nrow,Ncol), mar=c(2,2,1,0))
  for(t in 1:Nyears){
    plot( x=col(Results[,,t]), y=Results[,,t], col=c("black","blue")[row(Results[,,t])], log="y", pch=20, cex=2, main=paste("Year",t))
    CI = apply(Jags$BUGSoutput$sims.list$N_hat[,,t], MARGIN=2, FUN=quantile, prob=c(0.1,0.5,0.9))
    polygon( x=c(1:Nstages,Nstages:1), y=c(CI[1,],rev(CI[3,])), col=rgb(1,0,0,alpha=0.2))
  }
}

# Posterior vs. true values
if(FALSE){
  dev.new()
  Names = c("SigmaP_p","SigmaM_p","Alpha","Beta")
  True_List = list(SigmaP[1], SigmaM[1], Alpha, Beta)
  par( mfrow=c(3,4), mar=c(3,3,2,0))
  for(ParI in 1:4){
    Mat = Jags$BUGSoutput$sims.list[[Names[ParI]]]
    for(ColI in 1:ncol(Mat)){
      hist( Mat[,ColI], breaks=50, main=paste(Names[ParI],ColI) )
      abline(v=True_List[[ParI]][ColI], lwd=3)
    }
  }
}

# Estimated states for Paulik diagram
dev.new()
Ncol=ceiling(sqrt(Nstages)); Nrow=ceiling(Nstages/Ncol)
par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0))
N_hat = apply( Jags$BUGSoutput$sims.list$N_hat, MARGIN=2:3, FUN=mean)
for(StageI in 1:(Nstages-1)){
  plot( x=N_hat[StageI,], y=N_hat[StageI+1,])
  points( x=N_hat[StageI,], y=N_hat[StageI+1,], col="red")
}

# Predictive distribution
png( file=paste(File,"SRcurve_Paulik.png",sep=""), width=6, height=6, res=200, units="in")
  Pred = apply( Paulik$BUGSoutput$sims.list$Npred_hat, MARGIN=2:3, FUN=mean)
  plot( x=Data$Npred_init, y=Pred[6,], type="n", ylim=c(0,quantile(as.vector(Paulik$BUGSoutput$sims.list$Npred_hat[,6,]),prob=0.99)) )
  for(p in 1:Data$Npred){
    lines( x=rep(Data$Npred_init[p],2), y=quantile(Paulik$BUGSoutput$sims.list$Npred_hat[,6,p], prob=c(0.1,0.9)) )
    points( x=Data$Npred_init[p], y=mean(Paulik$BUGSoutput$sims.list$Npred_hat[,6,p]) )
  }
dev.off()

# Predictive distribution -- Variance
pdf( file=paste(File,"SRcurve_Paulik_variance.pdf",sep=""), width=10, height=10, onefile=TRUE)
  par(mfrow=c(4,5), mar=c(2,2,1,0), mgp=c(2,0.5,0.0), tck=-0.02)
  for(p in 1:100){
    Pred = apply( Paulik$BUGSoutput$sims.list$Npred_hat, MARGIN=2:3, FUN=mean)
    plot( x=1:Nstages, y=Pred[,s], type="n", ylim=quantile(as.vector(Paulik$BUGSoutput$sims.list$Npred_hat[,,]),prob=c(0.001,0.999)), log="y" )
    for(s in 1:Nstages){
      lines( x=rep(s,2), y=quantile(Paulik$BUGSoutput$sims.list$Npred_hat[,s,p], prob=c(0.1,0.9)) )
      points( x=s, y=mean(Paulik$BUGSoutput$sims.list$Npred_hat[,s,p]) )
    }
  }
dev.off()

# Predictive distribution -- Variance
pdf( file=paste(File,"SRcurve_Paulik_variance_2.pdf",sep=""), width=10, height=10, onefile=TRUE)
  par(mfrow=c(4,5), mar=c(2,2,1,0), mgp=c(2,0.5,0.0), tck=-0.02)
  for(p in 1:100){
    plot( 1, type="n", xlim=c(1,6), ylim=c(0,2) )
    CV = apply(Paulik$BUGSoutput$sims.list$Npred_hat[,,p], MARGIN=2, FUN=function(Vec){(quantile(Vec,0.75)-quantile(Vec,0.25))/mean(Vec)})
    lines( x=1:6, y=CV )
  }
dev.off()

####################################
# TRY RUNNING IN TMB
# Most stable (converges!): Version=Paulik_v1; Config=NA; EstMethod=1; UseREML=TRUE
# Paul_v2c doesn't fully converge in any configuration
####################################

# File structure
File = "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2014 -- Paulik diagram/"

# Load libraries
library(abind)
library(TMB)
library(numDeriv)   # grad()

options(digits=8)
Version = c("Paulik_v1", "Paulik_v2", "Paulik_v2b", "Paulik_v2c", "Paulik_v2d")[1]
  # Paulik_v1 -- state-space
  # Paulik_v2 -- state-space with multiple regimes
    # Paulik_v2b -- Switched notation for log-likelihood and loop structure
    # Paulik_v2c -- Separated process and measurement errors
    # DEPRECATED: Paulik_v2d -- Added 1e-10 to various parts (DIDN'T FIX PROBLEM! -- so don't use it)
Config = c("1-stage in all transitions", "2-regime in 2to3 transition", "2-regime in all transitions")[2] 
  # ONLY AFFECTS "Paulik_v2c"
  # "1-stage in all transitions" -- Single regime Paulik state-space model
  # "2-stage in all transitions" -- 2-regime Paulik state-space model
  # "2-stage in 2to3 transitions" -- 2-regime Paulik state-space model for transition from 2nd to 3rd stage, and 1-regime for other transitions
    # NOTE Paulik_v1 using importance sampler is more stable than Paulik_v2c for 1-stage model
EstMethod = 1
  # 1 : Importance sampler (doesn't work for Paulik_v2c)
  # 2 : Laplace approximation
UseREML = TRUE
  # REML : Integrates across all parameters except sigmas
  # ML : Integrates only across states 
      
# Compile spatial
setwd( File )
  dyn.unload( paste(File,dynlib(Version),sep="") )
  #file.remove( paste(File,Version,c(".o",".dll"),sep="") )
  compile( paste(Version,".cpp",sep="") )

# Re-load data
if(TRUE){
  CSV = read.csv( file=paste(File,"Indices 2014.csv",sep=""), header=TRUE)
  Results = t(CSV[,-c(1:2)])
  Results = Results[,which(apply(Results,MARGIN=2,FUN=function(Vec){any(!is.na(Vec))}))]
  Results = Results / outer(rowMeans(Results,na.rm=TRUE), rep(1,ncol(Results)))
  Results = aperm(abind( array(NA, dim=dim(Results)), Results, along=3),c(3,1,2))
  dimnames(Results) = list(c("True","Obs"), paste("Stage",1:dim(Results)[2]), paste("Year",1:dim(Results)[3]))
  Nstages = dim(Results)[2]
  Nyears = dim(Results)[3]
}
      
  # Run spatial model
  dyn.load( paste(File,dynlib(Version),sep="") )
  if(Version=="Paulik_v1"){
    Data = list( "Nstages"=Nstages, "Nyears"=Nyears, "N_obs"=ifelse(is.na(Results['Obs',,]),0,Results['Obs',,]) )
    Parameters = list("log_Alpha_s"=rep(log(1),Nstages-1), "Beta_s"=rep(0,Nstages-1), "log_SigmaP_s"=rep(log(2),Nstages-1), "log_SigmaM_s"=rep(log(2),Nstages), "log_N_hat_1"=rep(log(1),Nyears), "log_N_hat_2plus"=matrix(log(1),ncol=Nyears,nrow=Nstages-1))
    Map = list("log_SigmaM_s" = factor(rep(1,Nstages)) )
    Random = c("log_N_hat_2plus", "log_N_hat_1") #, "Alpha_s", "Beta_s")
    if(UseREML==TRUE) Random = c(Random, "log_Alpha_s", "Beta_s") 
  }
  if(Version%in%c("Paulik_v2")){
    Data = list( "Nstages"=Nstages, "Nyears"=Nyears, "N_obs"=ifelse(is.na(Results['Obs',,]),0,Results['Obs',,]) )
    Parameters = list("log_Alpha_s"=rep(log(1),Nstages-1), "Beta_s"=rep(0,Nstages-1), "log_SigmaP_s"=rep(log(2),Nstages-1), "log_SigmaM_s"=rep(log(2),Nstages), "log_N_hat_1"=rep(log(1),Nyears), "TransPars"=matrix(0,nrow=Nstages-1,ncol=2), "log_N_hat_2plus"=matrix(log(1),ncol=Nyears,nrow=Nstages-1))
    Map = list("log_SigmaM_s"=factor(rep(1,Nstages)) )
    Random = c("log_N_hat_2plus", "log_N_hat_1") #, "Alpha_s", "Beta_s")
  }
  if(Version%in%c("Paulik_v2b","Paulik_v2c","Paulik_v2d")){
    Data = list( "Nstages"=Nstages, "Nyears"=Nyears, "N_obs"=ifelse(is.na(Results['Obs',,]),0,Results['Obs',,]) )
    Parameters = list("log_Alpha_sr"=matrix(log(c(1,3)),nrow=Nstages-1,ncol=2,byrow=TRUE), "Beta_s"=rep(0,Nstages-1), "log_SigmaP_s"=rep(log(200),Nstages-1), "log_SigmaM_s"=rep(log(200),Nstages), "log_N_hat_1"=rep(log(1),Nyears), "TransPars"=matrix(c(0.1,0.9),nrow=Nstages-1,ncol=2,byrow=TRUE), "log_N_hat_2plus"=matrix(log(1),ncol=Nyears,nrow=Nstages-1))
    Map = list("log_SigmaM_s"=factor(rep(1,Nstages)) )
    Random = c("log_N_hat_2plus", "log_N_hat_1") #
    if(Config=="1-stage in all transitions"){
      Map[["log_Alpha_sr"]] = factor( matrix(c(1:(Nstages-1),c(NA,NA,NA,NA,NA)),nrow=Nstages-1,ncol=2) )
      Parameters[["log_Alpha_sr"]] = matrix(log(c(rep(1,Nstages-1),c(1e-5,1e-5,1e-5,1e-5,1e-5))),nrow=Nstages-1,ncol=2)
      Map[["TransPars"]] = factor( matrix(c(NA,NA,NA,NA,NA, NA,NA,NA,NA,NA),nrow=Nstages-1,ncol=2) )
    }
    if(Config=="2-regime in all transitions"){
      # No changes needed to above code
    }
    if(Config=="2-regime in 2to3 transition"){
      Map[["log_Alpha_sr"]] = factor( matrix(c(1:(Nstages-1),c(NA,Nstages,NA,NA,NA)),nrow=Nstages-1,ncol=2) )
      Parameters[["log_Alpha_sr"]] = matrix(log(c(rep(1,Nstages-1),c(1e-5,1,1e-5,1e-5,1e-5))),nrow=Nstages-1,ncol=2)
      Map[["TransPars"]] = factor( matrix(c(NA,1,NA,NA,NA, NA,2,NA,NA,NA),nrow=Nstages-1,ncol=2) )
    }
    if(UseREML==TRUE) Random = c(Random, "log_Alpha_sr", "Beta_s") 
  }
  MCcontrol = list( list(doMC=TRUE,seed=124,n=100), list(doMC=FALSE,seed=123,n=100) )[[EstMethod]]
  obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, map=Map, MCcontrol=MCcontrol)
  #obj$gr( obj$par ) - obj$gr( obj$par )             #  
  # Robustify
  if(TRUE){
    obj$gr_orig <- obj$gr
    obj$gr <- function(...){
      Return = obj$gr_orig(...)
      if( any(is.na(Return)) ) Return = grad(obj$fn,...)
      return(Return)
    }
  }
  # Settings                                       #                                           # , inner.method="BFGS", inner.control=list(trace=1)
  obj$control <- list(trace=1,parscale=rep(1,13),REPORT=1,reltol=1e-12,maxit=100)
  obj$hessian <- FALSE
  obj$fn( obj$par )
  # Explore guts of Importance Sampler MC() in TMB.R
  if(FALSE){
    obj$env$ff(obj$par,order=0)     # obj$env$last.par[-obj$env$random]
    par <- obj$env$last.par
    obj$env$MC(obj$par, n=obj$env$MCcontrol$n, seed=obj$env$MCcontrol$seed, order=0)
  }
  # Run optimizer                                #  
  opt = nlminb(start=obj$par, objective=obj$fn, gradient=obj[[c("gr_orig","gr")[2]]], lower=-20, upper=20, control=list(eval.max=1e4, iter.max=1e4, rel.tol=1e-10))
  ParHat = obj$env$parList( obj$env$last.par.best )
  Report = obj$report()
  Sdreport = sdreport(obj)

#################
# Diagnostics
#################
  # Inspect results
  Report$TransParams
  Report$TransArray[2,,]
  Report$N_exp_str[3,,]  
  Results['Obs',3,]
  # Plot transition
  par(mfrow=c(1,2))
  matplot( log(cbind(Results['Obs',3,],Report$N_exp_str[3,,])), type="l", col=c("black","red","blue") )
  plot( Report$p_st[2,], type="l", ylim=c(0,1) )
  
  # re-run 
  dyn.load( paste(File,dynlib(Version),sep="") )
  MCcontrol = list( list(doMC=TRUE,seed=123,n=1e2), list(doMC=FALSE,seed=123,n=100) )[[2]]
  obj <- MakeADFun(data=Data, parameters=ParHat, random=Random, hessian=FALSE, map=Map, MCcontrol=MCcontrol)
  obj$fn( obj$par )
  opt = nlminb(start=obj$par, objective=obj$fn, gradient=obj[[c("gr","gr_mod")[2]]], lower=-20, upper=20, control=list(eval.max=1e4, iter.max=1e4))
  
  