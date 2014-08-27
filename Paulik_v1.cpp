
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(Nstages);
  DATA_INTEGER(Nyears);
  DATA_MATRIX(N_obs);

  // Fixed effects
  PARAMETER_VECTOR(log_Alpha_s);
  PARAMETER_VECTOR(Beta_s);
  PARAMETER_VECTOR(log_SigmaP_s);
  PARAMETER_VECTOR(log_SigmaM_s);
  PARAMETER_VECTOR(log_N_hat_1)

  // Random effects
  PARAMETER_MATRIX(log_N_hat_2plus)

  // Variables
  matrix<Type> N_hat(Nstages,Nyears);
  matrix<Type> N_exp(Nstages,Nyears);
  vector<Type> SigmaP_s(Nstages-1);
  vector<Type> SigmaM_s(Nstages);
  vector<Type> Alpha_s(Nstages-1);

  // Indices and objective function
  int s,t;
  Type g = 0;

  // Derived variables
  for(s=0;s<(Nstages-1);s++) SigmaP_s(s) = exp(log_SigmaP_s(s));
  for(s=0;s<Nstages;s++) SigmaM_s(s) = exp(log_SigmaM_s(s));
  for(s=0;s<(Nstages-1);s++) Alpha_s(s) = exp(log_Alpha_s(s));
  
  // Projection and likelihood
  for(t=0;t<Nyears;t++){
    N_exp(0,t) = exp(log_N_hat_1(t));
    N_hat(0,t) = exp(log_N_hat_1(t));
    if(N_obs(0,t)>0){ g -= dnorm( log(N_obs(0,t)), log(N_exp(0,t)), SigmaM_s(0), 1 ); }
    for(s=1;s<Nstages;s++){
      N_hat(s,t) = exp(log_N_hat_2plus(s-1,t));
      //if(s==1){ log_N_exp(s,t) = log_N_hat_1(t) + Alpha_s(s-1) - Beta_s(s-1)*exp(log_N_hat_1(t)); }
      //if(s>=2){ log_N_exp(s,t) = log_N_hat_2plus(s-2,t) + Alpha_s(s-1) - Beta_s(s-1)*exp(log_N_hat_2plus(s-2,t)); }
      N_exp(s,t) = N_hat(s-1,t) * Alpha_s(s-1) * exp(-Beta_s(s-1)*N_hat(s-1,t));
      // Process error
      g -= dnorm( log(N_hat(s,t)), log(N_exp(s,t)), SigmaP_s(s-1), 1 );
      // Measurement error
      if(N_obs(s,t)>0){ g -= dnorm( log(N_obs(s,t)), log(N_hat(s,t)), SigmaM_s(s), 1 ); }
    }
  }
  //Reporting
  ADREPORT( SigmaP_s );
  ADREPORT( SigmaM_s );
  REPORT( SigmaP_s );
  REPORT( SigmaM_s );
  REPORT( N_exp )
  return g;
}
