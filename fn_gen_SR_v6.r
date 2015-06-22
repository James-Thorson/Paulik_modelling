## Function for multi-stage Paulik diagram
# Functions:  Generate_Data(); Fit_Jags(); plot.fit()
# 30 April 2015
# initial Jags code by Jim Thorson, modifications by Liz Brooks
#
########################################################################################
########################################################################################
Generate_Data = function(rseed=rseed, gen.fn=gen.fn, out.folder=out.folder, n.stages=n.stages, n.years=n.years,
   params=params) {

set.seed(rseed)
Results = array(NA, dim=c(2,n.stages,n.years), dimnames=list(c("True","Obs"),paste("Stage",1:n.stages),paste("Year",1:n.years)))
Results[1,1,] <- runif(n.years, 0, spawn.max)
Results[2,1,] <- Results[1,1,]*exp(rnorm(n.years, mean=0, sd=Sigma_M[1]))

if (gen.fn==1) {    #Generalized form
for(s in 1:(n.stages-1)){
  Results[1,s+1,] = Results[1,s,]*exp((-M[s]-R[s]*Results[1,s,])*t[s])/(  1 + (B[s]/(M[s]+R[s]*Results[1,s,]))*Results[1,s,]*(1-exp((-M[s]-R[s]*Results[1,s,])*t[s]) ) )* exp(rnorm(n.years, mean=0, sd=Sigma_P[s]))
  Results[2,s+1,] = Results[1,(s+1),] * exp(rnorm(n.years, mean=0, sd=Sigma_M[(s+1)]))
     } #end s loop
  } #end if 1

if (gen.fn==2) {    #Shepherd form
for(s in 1:(n.stages-1)){
  Results[1,(s+1),] = Results[1,s,] *exp(-M[s]*t[s]) /(1+(K[s]/M[s])*(1-exp(-M[s]*t[s]))*Results[1,s,]^G[s]) * exp(rnorm(n.years, mean=0, sd=Sigma_P[s]))
  Results[2,(s+1),] = Results[1,(s+1),] * exp(rnorm(n.years, mean=0, sd=Sigma_M[(s+1)]))
     } #end s loop
  } #end if 2


if (gen.fn==3) {    #BH form
 for(s in 1:(n.stages-1)){
    Results[1,(s+1),] = Results[1,s,] * exp(-M[s]*t[s])/(1+(B[s]/M[s])*(1-exp(-M[s]*t[s]))*Results[1,s,]) * exp(rnorm(n.years, mean=0, sd=Sigma_P[s]))
    Results[2,(s+1),] = Results[1,(s+1),] * exp(rnorm(n.years, mean=0, sd=Sigma_M[(s+1)]))
     } #end s loop
 } #end if  3



if (gen.fn==4) {    #Ricker form
 for(s in 1:(n.stages-1)){
  Results[1,(s+1),] = exp(-M[s]*t[s])*Results[1,s,]* exp(-R[s]*Results[1,s,]*t[s])* exp(rnorm(n.years, mean=0, sd=Sigma_P[s]))
    Results[2,(s+1),] = Results[1,(s+1),] * exp(rnorm(n.years, mean=0, sd=Sigma_M[(s+1)]))
     } #end s loop
 } #end if  4


 return(Results)
} #end fn Generate_Data




########################################################################################
########################################################################################
Fit_Jags = function(est.fn=est.fn, od=od, n.stages=n.stages, n.years=n.years, Data=Data,
    params=params, nchains=nchains, Nsim=Nsim, Nburnin=Nburnin, Nthin=Nthin, ...) {
     # inits=jags.inits,

d1 <- date()

if (est.fn==1) {    #Generalized form

Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,10)
    Ricker[s] ~ dunif(0,10)
    Beta[s] ~ dunif(0,10)
    SigmaP[s] ~ dunif(0,10)     #stage specific process error
    TauP[s] <- pow( SigmaP[s], -2)
  } #end s loop

    SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    #SigmaM[s] ~ dunif(0,1)
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } #end s loop


  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[s-1,y] * exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)])/(  1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*N_hat[(s-1),y]))*N_hat[(s-1),y]*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)]) ) )
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } #end s loop
  } #end y loop
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[(s-1),p] * exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)])/(  1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*Npred_hat[(s-1),p]))*Npred_hat[(s-1),p]*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)]) ) )
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } # end s
  } # end p loop
} # end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data, inits=jags.inits,
    parameters.to.save=c("Mu", "Ricker","Beta", "SigmaM_p", "SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 1


if (est.fn==2) {    #Shepherd form

# Define simple Paulik diagram
Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,20)
    Kappa[s] ~ dunif(0,20)
    Gamma[s] ~ dunif(0,5)
    SigmaP[s] ~ dunif(0,20)
    TauP[s] <- pow( SigmaP[s], -2)
  } #end s loop

  SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } #end s loop

  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[(s-1),y] * exp(-Mu[(s-1)]*t[(s-1)]) / (1+(Kappa[(s-1)]/Mu[(s-1)])*(1-exp(-Mu[(s-1)]*t[(s-1)]))*N_hat[(s-1),y]^Gamma[(s-1)])
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } #end s loop
  } #end y loop
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[(s-1),p] * exp(-Mu[(s-1)]*t[(s-1)]) /(1+(Kappa[(s-1)]/Mu[(s-1)])*(1-exp(-Mu[(s-1)]*t[(s-1)]))*Npred_hat[(s-1),p]^Gamma[(s-1)])
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } #end s loop
  } #end p loop
} #end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data, inits=jags.inits,
    parameters.to.save=c( "Mu","Kappa","Gamma", "SigmaM_p","SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 2


if (est.fn==3) {    #BH form
# Define simple Paulik diagram
Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,20)
    Beta[s] ~ dunif(0,20)
    SigmaP[s] ~ dunif(0,20)
    TauP[s] <- pow( SigmaP[s], -2)
  } # end s loop

  SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } # end s loop

  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[(s-1),y] * exp(-Mu[(s-1)]*t[(s-1)]) / (1+(Beta[(s-1)]/Mu[(s-1)])*(1-exp(-Mu[(s-1)]*t[(s-1)]))*N_hat[(s-1),y])
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } # end s loop
  } #end y loop

  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[s-1,p] * exp(-Mu[(s-1)]*t[(s-1)]) / (1+(Beta[(s-1)]/Mu[(s-1)])*(1-exp(-Mu[(s-1)]*t[(s-1)]))*Npred_hat[(s-1),p])
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } #end s loop
  } #end p loop
} #end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data,
    parameters.to.save=c( "Mu","Beta", "SigmaM_p","SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 3


if (est.fn==4) {    #Ricker form
# Define simple Paulik diagram
Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,20)
    Ricker[s] ~ dunif(-20,20)
    SigmaP[s] ~ dunif(0,20)
    TauP[s] <- pow( SigmaP[s], -2)
  } #end s loop

  SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } #end s loop

  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] ) #initializing?     exp(-Beta[s-1]*N_hat[s-1,y]/100)
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[s-1,y] * exp(-Mu[(s-1)]*t[(s-1)]) * exp(-Ricker[(s-1)]*N_hat[(s-1),y]*t[(s-1)])
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } #end s loop
  } #end y loop
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[(s-1),p] * exp(-Mu[(s-1)]*t[(s-1)]) * exp(-Ricker[(s-1)]*Npred_hat[(s-1),p]*t[(s-1)])
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } #end s loop
  } #end p loop
} #end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data,
    parameters.to.save=c("Mu","Ricker", "SigmaM_p","SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 4



if (est.fn==13) {    #Generalized form

Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,2)
    Ricker[s]<- 0
    Beta[s] ~ dunif(0,3)
    SigmaP[s] ~ dunif(0,10)     #stage specific process error
    TauP[s] <- pow( SigmaP[s], -2)
  } #end s loop

    SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } #end s loop


  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[s-1,y] * exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)])/( ( 1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*N_hat[(s-1),y]))*N_hat[(s-1),y])*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)]) ) )
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } #end s loop
  } #end y loop
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[(s-1),p] * exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)])/( ( 1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*Npred_hat[(s-1),p]))*Npred_hat[(s-1),p])*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)]) ) )
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } # end s
  } # end p loop
} # end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data,
    parameters.to.save=c("Mu", "Ricker","Beta", "SigmaM_p", "SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 13



if (est.fn==14) {    #Generalized form

Paulik_Fn = function(){

  for(s in 1:(Nstages-1)){
    Mu[s] ~ dunif(0,2)
    Ricker[s] ~ dunif(0,1)
    Beta[s] <-  0
    SigmaP[s] ~ dunif(0,10)     #stage specific process error
    TauP[s] <- pow( SigmaP[s], -2)
  } #end s loop

    SigmaM_p ~ dunif(0,1) #stage specific measurement error
  for(s in 1:Nstages){
    SigmaM[s] <- SigmaM_p
    TauM[s] <- pow( SigmaM[s], -2)
  } #end s loop


  # Projection and likelihood
  for(y in 1:Nyears){
    N_hat[1,y] ~ dlnorm(0,0.000001) #initializing?
    N_obs[1,y] ~ dlnorm( log(N_hat[1,y]), TauM[1] )
    for(s in 2:Nstages){
      N_exp[s,y] <- N_hat[s-1,y] * exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)])/( ( 1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*N_hat[(s-1),y]))*N_hat[(s-1),y])*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*N_hat[(s-1),y])*t[(s-1)]) ) )
      N_hat[s,y] ~ dlnorm( log(N_exp[s,y]), TauP[(s-1)] )
      N_obs[s,y] ~ dlnorm( log(N_hat[s,y]), TauM[s] )
    } #end s loop
  } #end y loop
  # Predictive distribution
  for(p in 1:Npred){
    Npred_hat[1,p] <- Npred_init[p]
    for(s in 2:Nstages){
      Npred_exp[s,p] <- Npred_hat[(s-1),p] * exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)])/( ( 1 + (Beta[(s-1)]/(Mu[(s-1)]+Ricker[(s-1)]*Npred_hat[(s-1),p]))*Npred_hat[(s-1),p])*(1-exp((-Mu[(s-1)]-Ricker[(s-1)]*Npred_hat[(s-1),p])*t[(s-1)]) ) )
      Npred_hat[s,p] ~ dlnorm( log(Npred_exp[s,p]), TauP[(s-1)] )
    } # end s
  } # end p loop
} # end function

# Run jags
Paulik <- jags(model.file=Paulik_Fn, working.directory=NULL, data=Data,
    parameters.to.save=c("Mu", "Ricker","Beta", "SigmaM_p", "SigmaP","Npred_hat"),
    n.chains=nchains, n.thin=Nthin, n.iter=(Nthin*Nsim+Nburnin), n.burnin=Nburnin)

   } #end if 14




d2<-date()


 return(Paulik)
} #end fn Fit_Jags




########################################################################################
########################################################################################

plot.fit=function(jags.out=jags.out, wd=wd, od=od, n.stages=n.stages, Data=Data, file.end.txt=file.end.txt) {

  dev.new()
#check.dir <- shell("dir plots2", intern=T)
check.dir <-shell(paste("dir ", od,  sep=""), intern=T )
if ( is.null(attr(check.dir, "status")) )  attr(check.dir, "status") <- 0
#check.dir <-system(paste('"dir "', od,   sep=""), intern=T )
flag.dir <- which(check.dir=="File Not Found")
if (length(flag.dir)>0 | attr(check.dir, "status", exact=T)==1)  shell(paste("mkdir ",  od, sep=""), intern=T )



beta.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,4)=="Beta" )
gamma.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,5)=="Gamma" )
kappa.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,5)=="Kappa" )
mu.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,2)=="Mu" )
ricker.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,6)=="Ricker" )
sigmaM.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,8)=="SigmaM_p" )
sigmaP.index<- which(substr(rownames(jags.out$BUGSoutput$summary),1,6)=="SigmaP" )
deviance.index<- which(rownames(jags.out$BUGSoutput$summary)=="deviance" )
estim.par.index <- c(beta.index,  mu.index,  ricker.index, gamma.index,  kappa.index,
             sigmaM.index, sigmaP.index,deviance.index)

n.par <- length(estim.par.index)
est.par <- matrix(NA, nrow=n.par, ncol=dim(jags.out$BUGSoutput$summary)[2])
est.par <- jags.out$BUGSoutput$summary[c(estim.par.index),]
est.par.head <- colnames(jags.out$BUGSoutput$summary)
est.par.rows <- rownames(jags.out$BUGSoutput$summary) [estim.par.index]


# write parameter estimates to csv file
write.csv(est.par, file=paste(od,"Param.est.csv",sep="") )


   ### NOTE: coda may have this stuff
# plot Traces of each chain , look at autocorrelation, calculate PSRF
chain.len = jags.out$BUGSoutput$n.sims/jags.out$BUGSoutput$n.chains
 # autocorrelation
jags.mcmc <-  as.mcmc(jags.out$BUGSoutput$sims.matrix[,estim.par.index])
ac.jags <- autocorr.diag(jags.mcmc, lags=seq(0,10), relative=F )
write.csv(ac.jags, file=paste(od,"Autocorr.pars.csv",sep="") )
par.corr <-crosscorr(jags.mcmc)
write.csv(par.corr, file=paste(od,"Cross_corr.pars.csv",sep="") )

 #Traces and density plots
for (p in 1:ceiling(n.par/3) ) {
png( file=paste(od,"Trace_density_plots",p,".png",sep=""), width=8, height=10, res=400, units="in")

plot(jags.mcmc[,( 1+(p-1)*3):min(3*p, n.par)], ask=F)

 dev.off()
  } # end p loop




#potential scale reduction factor
png( file=paste(od,"PSRF_est_pars.png",sep=""), width=8, height=10, res=400, units="in")
psrf <- jags.out$BUGSoutput$summary[,8]
plot(seq(1,n.par), seq(0,2, length.out=n.par), type='n', axes=F, xlab="", 
     ylab="Potential Scale Reduction Factor")
abline(h=1, lty=2, col='red')
abline(h=1.1, lty=1, col='red')
segments(x0=seq(1,n.par), y0=rep(0,n.par), x1=seq(1,n.par), y1=psrf[estim.par.index], col="#1111ccaa")
axis(side=1, at=seq(1,n.par), labels=est.par.rows, las=2, cex.axis=0.75 )    
axis(side=2, at=pretty(seq(0,2, length.out=n.par)), labels=pretty(seq(0,2, length.out=n.par)))
box()
 dev.off()


psrf.Nhat <- psrf[-estim.par.index]
for (s in 1:jags.out$n.stages ) {
png( file=paste(od,"PSRF_Npred_hat_stage_", s, "_.png",sep=""), width=8, height=10, res=400, units="in")
tmp.psrf <- seq(s,jags.out$n.years*jags.out$n.stages, by= jags.out$n.stages)
plot(seq(1,jags.out$n.years), seq(1,jags.out$n.years), type='n', axes=F, xlab="", 
     ylab="Potential Scale Reduction Factor", ylim=c(0,2))
abline(h=1, lty=2, col='red')
abline(h=1.1, lty=1, col='red')
title(main=paste("Nhat_pred in stage ", s, sep="")  )
segments(x0=seq(1,jags.out$n.years), y0=rep(0,jags.out$n.years), x1=seq(1,jags.out$n.years), y1=psrf.Nhat[tmp.psrf], col="#1111ccaa")
axis(side=1, at=seq(1,jags.out$n.years), labels=seq(1,jags.out$n.years), las=2, cex.axis=0.75 )    
axis(side=2, at=pretty(seq(0,2)), labels=pretty(seq(0,2)) )
box()
 dev.off()
  } # end s loop


#Calculate Relative Errors
re.mat <- matrix(NA, nrow=chain.len*nchains, ncol= (n.par-1))
re.count <- 1

if (length(beta.index)>0) {
re.mat[,re.count:length(beta.index)] <- (jags.out$BUGSoutput$sims.matrix[,beta.index] - jags.out$true.params$B)/jags.out$true.params$B 
true.zero <- which(jags.out$true.params$B==0)
re.mat[,c(re.count+true.zero -1)] <- NA

re.count <- re.count+length(beta.index)
  } #end if

  
if (length(gamma.index)>0) {
re.mat[,re.count:(re.count+length(gamma.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,gamma.index] - jags.out$true.params$G)/jags.out$true.params$G 
re.count <- re.count+length(gamma.index)
  } #end if

if (length(kappa.index)>0) {
re.mat[,re.count:(re.count+length(kappa.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,kappa.index] - jags.out$true.params$G)/jags.out$true.params$K 
re.count <- re.count+length(kappa.index)
  } #end if


if (length(mu.index)>0) {
re.mat[,re.count:(re.count+length(mu.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,mu.index] - jags.out$true.params$M)/jags.out$true.params$M 
re.count <- re.count+length(mu.index)
  } #end if


if (length(ricker.index)>0) {
re.mat[,re.count:(re.count+length(ricker.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,ricker.index] - jags.out$true.params$R)/jags.out$true.params$R 
true.zero <- which(jags.out$true.params$R==0)
re.mat[,c(re.count+true.zero -1)] <- NA
re.count <- re.count+length(ricker.index)
  } #end if


if (length(sigmaM.index)>0) {
re.mat[,re.count:(re.count+length(sigmaM.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,sigmaM.index] - jags.out$true.params$Sigma_M)/jags.out$true.params$Sigma_M 
re.count <- re.count+length(sigmaM.index)
  } #end if


if (length(sigmaP.index)>0) {
re.mat[,re.count:(re.count+length(sigmaP.index)-1)] <- (jags.out$BUGSoutput$sims.matrix[,sigmaP.index] - jags.out$true.params$Sigma_P)/jags.out$true.params$Sigma_P 
re.count <- re.count+length(sigmaP.index)
  } #end if






# Create Violin plot
png( file=paste(od,"Relative_Error_ViolinPlots.png",sep=""), width=8, height=10, res=400, units="in")
x.ticks <- seq(min(re.mat, na.rm=T), min(3, max(re.mat, na.rm=T)), length.out=(n.par-1) )
plot(x=x.ticks, y=seq(1,(n.par-1)), type='n',xlab="Relative Error",ylab="", axes=F) 
for (i in 1:(n.par-1)) {
if (!is.na(re.mat[,i]) ) vioplot(re.mat[,i], horizontal=T, add=T, na.rm=T, at=(n.par-i),
   col="#bbbbffaa" )

}#end i loop
abline(v=0, col='red', lwd=2)
axis(side=1, at=pretty(x.ticks))
axis(side=2, at=( (n.par-1):1), label= est.par.rows[1:(n.par-1)], cex.axis=0.75, las=1 )
dev.off()




## plot data (True , Obs AND estimated fit... 
for (stage in 2: n.stages) {
png( file=paste(od,"Pred_SR_", (stage-1), ".png", sep=""), width=8, height=10, res=400, units="in")
  Pred.x = apply( jags.out$BUGSoutput$sims.list$Npred_hat[,(stage-1),], MARGIN=2, FUN=mean)
  Pred.y = apply( jags.out$BUGSoutput$sims.list$Npred_hat[,stage,], MARGIN=2, FUN=mean)  
  plot( x=Pred.x, y=Pred.y, type="n", xlab=paste("Predicted Stage ", (stage-1),sep=""), ylab=paste("Predicted Stage ", stage,sep=""),
    ylim=c(0,1.02*max(Results[1,stage,], Results[2,stage,],quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,],prob=0.995), na.rm=T) )   )
  for(p in 1:length(Pred.y)){
    lines( x=rep(Pred.x[p],2), y=quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p], prob=c(0.01,0.99)), col='grey65' )
    lines( x=rep(Pred.x[p],2), y=quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p], prob=c(0.1,0.9)), lwd=2, col='grey25' )
    points( x=Pred.x[p], y=mean(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p]), pch=1 )
  }
 #overlay true points (blue)
  points(Results[1,(stage-1),], Results[1,stage,], col='blue', pch=3)
  # one-off
  #points(Results[1,(stage-1),], Results[1,(stage+1),], col='blue', pch=3)  
  #overlay obs points (red)
  points(Results[2,(stage-1),], Results[2,stage,], col='red', pch=4)
  #one-off
  #points(Results[2,(stage-1),], Results[2,(stage+1),], col='red', pch=4)  
  legend('topleft', legend=c("True", "Obs", "Mean", "90% PI", "99% PI"), col=c("Blue", "Red", "Black", "grey25", 'grey65'),
     pch=c(3,4,1,NA, NA), lwd=c(NA,NA,NA,2,1), cex=0.75, ncol=2)
dev.off()


} #end loop over stage



## plot  just the observed values
for (p in 2:n.stages)  {
png( file=paste(od,"SR_",p,"_vs_", (p-1),".png",sep=""), width=8, height=10, res=400, units="in")
#png( file=paste(od,"SR_3_vs_1.png",sep=""), width=8, height=10, res=400, units="in")
plot(Results[2,(p-1),], Results[2,p,], ylim=c(0, max(Results[2,p,]) ), 
  xlab= paste("Stage_", (p-1), sep=""), ylab=paste("Stage_",p,sep="")   )
#plot(Results[2,1,], Results[2,3,], ylim=c(0, max(Results[2,3,]) ), 
#  xlab= paste("Stage_1", "", sep=""), ylab=paste("Stage_3","",sep="")   )

dev.off()
   }  # end p-loop over stages
   

 # include correlation plots from Paulik_test_SR_form_v3.r, same rules as above
    ## not included yet
    
    
 # not sure this plot is worthwhile
png( file=paste(od,"Cross_corr_plot",p,".png",sep=""), width=8, height=10, res=400, units="in")
 crosscorr.plot(jags.mcmc[,1:(n.par-1)], cex.axis=0.75)
dev.off()



#---------  Create PDF  --------------------------------------------------


graphics.off()
pdf(file=paste(od, "Paulik_Plots_gen.fn.", jags.out$gen.fn, "_est.fn.", jags.out$est.fn,
      "_", file.end.txt,".pdf",sep=""), onefile=TRUE)


### Summary cover page for pdf
plot(seq(1,20), seq(1,20), type='n', xlab="", ylab="", axes=F, ylim=c(-1,12) )
title(main="Simulation Summary")
if (jags.out$gen.fn==1) gname <- "Generalized"
if (jags.out$gen.fn==2) gname <- "Shepherd"
if (jags.out$gen.fn==3) gname <- "Beverton-Holt"
if (jags.out$gen.fn==4) gname <- "Ricker"
if (jags.out$est.fn==1) ename <- "Generalized"
if (jags.out$est.fn==2) ename <- "Shepherd"
if (jags.out$est.fn==3) ename <- "Beverton-Holt"
if (jags.out$est.fn==4) ename <- "Ricker"
if (jags.out$est.fn==13) ename <- "Generalized (B-H form)"
if (jags.out$est.fn==13) ename <- "Generalized (Ricker form)"
text(1,12, label=paste("Generation model = ", gname, sep=""), pos=4)
text(1,11, label=paste("Estimation model = ", ename, sep=""), pos=4)
text(1,10, label=paste("Nyears = ", jags.out$n.years, sep=""), pos=4)
text(1,9, label=paste("Nstages = ", jags.out$n.stages, sep=""), pos=4)
text(1,8, label="True Process error = ", pos=4)
n.sigmaP <- length(jags.out$true.params$Sigma_P)
text(x=seq(1,8, length.out=n.sigmaP), y=rep(7.5, n.sigmaP), label= jags.out$true.params$Sigma_P, pos=4, cex=0.8)
text(1,7, label="True Measurement error = ",  pos=4)
n.sigmaM <- length(jags.out$true.params$Sigma_M)
text(x=seq(1,8, length.out=n.sigmaM), y=rep(6.5, n.sigmaM), label= jags.out$true.params$Sigma_M, pos=4, cex=0.8)

text(1,5, label=paste("DIC = ", jags.out$BUGSoutput$DIC, sep=""), pos=4)
text(1,4, label=paste("pD = ", jags.out$BUGSoutput$pD, sep=""), pos=4)

text(1,2, label=paste("Run Start = ", jags.out$d1, sep=""), pos=4)
text(1,1, label=paste("Run End = ", jags.out$d2, sep=""), pos=4)
text(1,-0.5, label=paste("Niter= ", (jags.out$n.thin*jags.out$n.sims+jags.out$n.burnin), 
   "  NBurnin=",jags.out$n.burnin, "  Nthin=", jags.out$n.thin, sep=""), pos=4)


 #Traces and density plots
for (p in 1:ceiling(n.par/3) ) {
plot(jags.mcmc[,( 1+(p-1)*3):min(3*p, n.par)], ask=F)
  } # end p loop

#Correlation plot (crappy)
 crosscorr.plot(jags.mcmc[,1:(n.par-1)], cex.axis=0.75)



#  Create Violin Plot of Relative Errors
x.ticks <- seq(min(re.mat, na.rm=T), min(3, max(re.mat, na.rm=T)), length.out=(n.par-1) )
plot(x=x.ticks, y=seq(1,(n.par-1)), type='n',xlab="Relative Error",ylab="", axes=F) 
for (i in 1:(n.par-1)) {
if (!is.na(re.mat[,i]) ) vioplot(re.mat[,i], horizontal=T, add=T, na.rm=T, at=(n.par-i),
   col="#bbbbffaa" )

}#end i loop
abline(v=0, col='red', lwd=2)
axis(side=1, at=pretty(x.ticks))
axis(side=2, at=( (n.par-1):1), label= est.par.rows[1:(n.par-1)], cex.axis=0.75, las=1 )




## plot data (True , Obs AND estimated fit... 
for (stage in 2: n.stages) {

  Pred.x = apply( jags.out$BUGSoutput$sims.list$Npred_hat[,(stage-1),], MARGIN=2, FUN=mean)
  Pred.y = apply( jags.out$BUGSoutput$sims.list$Npred_hat[,stage,], MARGIN=2, FUN=mean)  
  plot( x=Pred.x, y=Pred.y, type="n", xlab=paste("Predicted Stage ", (stage-1),sep=""), ylab=paste("Predicted Stage ", stage,sep=""),
    ylim=c(0,1.02*max(Results[1,stage,], Results[2,stage,],quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,],prob=0.995), na.rm=T) )   )
  for(p in 1:length(Pred.y)){
    lines( x=rep(Pred.x[p],2), y=quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p], prob=c(0.01,0.99)), col='grey65' )
    lines( x=rep(Pred.x[p],2), y=quantile(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p], prob=c(0.1,0.9)), lwd=2, col='grey25' )
    points( x=Pred.x[p], y=mean(jags.out$BUGSoutput$sims.list$Npred_hat[,stage,p]), pch=1 )
  }
 #overlay true points (blue)
  points(Results[1,(stage-1),], Results[1,stage,], col='blue', pch=3)
  #overlay obs points (red)
  points(Results[2,(stage-1),], Results[2,stage,], col='red', pch=4)
  legend('topleft', legend=c("True", "Obs", "Mean", "90% PI", "99% PI"), col=c("Blue", "Red", "Black", "grey25", 'grey65'),
     pch=c(3,4,1,NA, NA), lwd=c(NA,NA,NA,2,1), cex=0.75, ncol=2)


} #end loop over stage





#potential scale reduction factor

psrf <- jags.out$BUGSoutput$summary[,8]
plot(seq(1,n.par), seq(0,2, length.out=n.par), type='n', axes=F, xlab="", 
     ylab="Potential Scale Reduction Factor")
abline(h=1, lty=2, col='red')
abline(h=1.1, lty=1, col='red')
segments(x0=seq(1,n.par), y0=rep(0,n.par), x1=seq(1,n.par), y1=psrf[estim.par.index], col="#1111ccaa")
axis(side=1, at=seq(1,n.par), labels=est.par.rows, las=2, cex.axis=0.75 )    
axis(side=2, at=pretty(seq(0,2, length.out=n.par)), labels=pretty(seq(0,2, length.out=n.par)))
box()




psrf.Nhat <- psrf[-estim.par.index]
for (s in 1:jags.out$n.stages ) {

tmp.psrf <- seq(s,jags.out$n.years*jags.out$n.stages, by= jags.out$n.stages)
plot(seq(1,jags.out$n.years), seq(1,jags.out$n.years), type='n', axes=F, xlab="", 
     ylab="Potential Scale Reduction Factor", ylim=c(0,2))
abline(h=1, lty=2, col='red')
abline(h=1.1, lty=1, col='red')
title(main=paste("Nhat_pred in stage ", s, sep="")  )
segments(x0=seq(1,jags.out$n.years), y0=rep(0,jags.out$n.years), x1=seq(1,jags.out$n.years), y1=psrf.Nhat[tmp.psrf], col="#1111ccaa")
axis(side=1, at=seq(1,jags.out$n.years), labels=seq(1,jags.out$n.years), las=2, cex.axis=0.75 )    
axis(side=2, at=pretty(seq(0,2)), labels=pretty(seq(0,2)) )
box()

  } # end s loop



dev.off()      
graphics.off()
#---------  Close PDF  --------------------------------------------------


} #end function plot.fit




