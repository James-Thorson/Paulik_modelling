## Function for multi-stage Paulik diagram
#  23 april 2015
# code initiated by jim thorson 2014-06-19
# modifications by liz brooks 23 april 2015
# last modified: 30 april 2015


rm(list=ls())      # Remove all variables, etc. from the session memory
graphics.off()     # close any open graphics windows

# Load libraries
library(R2jags)
library(abind)
library(coda)
library(vioplot)


########################################################################################
########################################################################################
###   Generate data


# parameters used for each function
# M  Dens-indep instantaneous natural mortality between stages  (vector of length (n.stages-1))
# t  fraction of year for stage duration (vector of length (n.stages-1))
# B  dens-dep term  (vector of length (n.stages-1))
# R  dens-dep mortality of ricker-type (vector of length (n.stages-1))
# G - exponent on dens-dep term, only used for Shepherd (G=1 ==> BH; G>1 ==> Ricker-like; G<1 ==> Cushing) (vector of length (n.stages-1)
# Sigma_P   ==> process error between stages, vector of length (n.stages-1)
# Sigma_M   ==> measurement error for each stage, vector of length (n.stages)
t <- c(0.5, 0.25)
R <- c(0.0, 0.0)
B <- c(0.002, 0.009)
M <- c(0.4, 0.2)
K <- c(0.0,0.00) #only used if gen.fn=2 (Shepherd)
G <- c(1.0,1.0) #only used if gen.fn=2
Sigma_M <- c(0.15, 0.15, 0.15)
Sigma_P <- c(0.20, 0.20)

gen.fn <- 1 # (1=generalized; 2= shepherd, 3=BH, 4= Ricker)
est.fn <- 1  #same numbering as gen.fn
n.stages <- 3  #number of life stages (2=typical S-R function, >2 includes multiple S-R fns
file.end.txt <- "runGen3_3stage_case"  # text pasted-on to end of pdf file and list and output directory


n.years <- 100  #number of years of data to be generated
spawn.max <- 6e3  #maximum for generating 1st stage observations (spawners)

rseed<-43771  #seed for the random number generator


#MCMC specs
Nburnin=1e5
Nthin=1e2
Nsim=1e3    #number of saved iterations per chain (after thinning and burnin)
nchains=3

params=list(M=M, B=B, G=G, R=R, K=K, t=t, Sigma_P=Sigma_P, Sigma_M=Sigma_M)

#wd <- "C:\\liz\\R\\Paulik\\paulik_sim\\"    #working directory (desktop)
wd <- "C:\\liz_temp\\R\\Paulik\\paulik_sim\\for_ppt\\"    #working directory (laptop)
fd <- "C:\\liz_temp\\R\\Paulik\\paulik_sim\\"   #directory where Paulik functions are
od <- paste(wd,"gen.fn.", gen.fn, "_est.fn.", est.fn, "_Nstages.",n.stages,"\\", file.end.txt, "\\", sep="")  #output directory

source(paste(fd, "fn_gen_SR_v6.r", sep="") )   #the file with functions




Results <- Generate_Data(rseed=rseed, gen.fn=gen.fn, out.folder=out.folder, n.stages=n.stages, n.years=n.years,
   params=params)

for (p in 2:n.stages)  {
plot(Results[2,(p-1),], Results[2,p,], ylim=c(0, max(Results[2,p,]) ), 
  xlab= paste("Stage_", (p-1), sep=""), ylab=paste("Stage_",p,sep="")   )
   }  #end p loop 


Data = list( N_obs=Results['Obs',,], Nstages=n.stages, Nyears=n.years, t=t,
      Npred_init=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=n.years),
      Npred=n.years )
#Data1 <- Data
#Data <- list(N_obs = Results['Obs', c(1,3),], Nstages=2, Nyears=n.years, t=t,
#       Npred_init=seq(min(Results['Obs',1,],na.rm=TRUE),max(Results['Obs',1,],na.rm=TRUE),length=n.years),
#      Npred=n.years )
#file.end.txt <- "runBH1_sim3stage_fit2stage_case"  # text pasted-on to end of pdf file and list and output directory      
#od <- paste(wd,"gen.fn.", gen.fn, "_est.fn.", est.fn, "_Nstages.",n.stages,"\\", file.end.txt, "\\", sep="")  #output directory
#rseed<-43771  #seed for the random number generator
#n.stages1 <- n.stages
#n.stages <- 2
d1=date()
##----- Fit Jags, specifying initial conditions ----
#generalized inits
inits.1 <- list( "Mu"=rep(0.3, (n.stages-1)), "Beta"= rep(3.5e-3, (n.stages-1)), "Ricker"=rep(1e-5, (n.stages-1)),
          "SigmaP"=rep(1e-2, (n.stages-1)), "SigmaM_p"=1e-2)
#shepherd inits
#inits.1 <- list( "Mu"=rep(0.3, (n.stages-1)), "Kappa"= rep(3.5e-3, (n.stages-1)), "Gamma"=rep(1.1, (n.stages-1)),
#          "SigmaP"=rep(1e-2, (n.stages-1)), "SigmaM_p"=1e-2)
inits.2 <- inits.1
inits.3 <- inits.1
jags.inits <- list(inits.1, inits.2, inits.3)
#jags.out <- Fit_Jags( est.fn=est.fn, od=od, n.stages=n.stages, n.years=n.years, Results=Results,
#   inits=jags.inits, params=params, nchains=nchains, Nsim=Nsim, Nburnin=Nburnin, Nthin=Nthin )

#---- OR .... Fit Jags, without specifying initial conditions -----
jags.out <- Fit_Jags( est.fn=est.fn, od=od, n.stages=n.stages, n.years=n.years, Data=Data,
    params=params, nchains=nchains, Nsim=Nsim, Nburnin=Nburnin, Nthin=Nthin , inits=jags.inits)
d2=date()
jags.out$gen.fn=gen.fn
jags.out$est.fn=est.fn
jags.out$seed = rseed
jags.out$d1=d1
jags.out$d2=d2
jags.out$true.params<- params
jags.out$n.years <- n.years
jags.out$n.stages <- n.stages
jags.out$n.chains <- nchains
jags.out$n.stages <- n.stages
jags.out$n.sims <- Nsim
jags.out$n.burnin <- Nburnin
jags.out$n.thin <- Nthin


list.out <- list(est=jags.out$BUGSoutput$summary, DIC=jags.out$BUGSoutput$DIC, pd=jags.out$BUGSoutput$pD,
   d1=jags.out$d1, d2=jags.out$d2, rseed=rseed, gen.fn=jags.out$gen.fn, est.fn=jags.out$est.fn,
   true.params=params, Niter=jags.out$BUGSoutput$n.iter, 
   Nsim=Nsim, Nburnin=Nburnin, Nthin=Nthin, nchains=nchains)




## function to evaluate and plot results
graphics.off()
plot.fit(jags.out=jags.out, wd=wd, od=od, n.stages=n.stages, Data=Data, file.end.txt=file.end.txt )

dput(list.out, file=paste(od,"list_", file.end.txt, ".out", sep="")  )

########################################################################################
########################################################################################
