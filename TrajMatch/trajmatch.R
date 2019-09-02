
library(pomp)
setwd("~/Git/TrajMatch")
start_time <- Sys.time()
######################################################  Model Snippet
seasonal.sir.ode <- Csnippet("
                             double beta0;
                             double va, tt;
                             double seas;
                             
                             
                             // term-time seasonality
                             tt = (t-floor(t))*365.25;
                             if ((tt>=7&&tt<=100) || (tt>=115&&tt<=199) || (tt>=252&&tt<=300) || (tt>=308&&tt<=356))
                             seas = 1.0+amplitude*0.2411/0.7589;
                             else
                             seas = 1.0-amplitude;
                             
                             // transmission rate
                             beta0 = R0*(gamma+mu)*(sigma+mu)/sigma;
                             // expected force of infection
                             double Beta = beta0*seas/pop;
                             
                             //Vacination uptake
                             //Vacination uptake
                             if (t < 1968)
                             va = 0;
                             else if (t>=1968 && t<=1969)
                             va = 0.33;
                             else if (t>=1969 && t<=1970)
                             va = 0.46;
                             else if (t>=1970 && t<=1971)
                             va = 0.51;
                             else if (t>=1971 && t<=1972)
                             va = 0.53;
                             else if (t>=1972 && t<=1973)
                             va = 0.52;
                             else if (t>=1973 && t<=1974)
                             va = 0.46;
                             else if (t>=1974 && t<=1975)
                             va = 0.46;
                             else if (t>=1975 && t<=1976)
                             va = 0.48;
                             else if (t>=1976 && t<=1977)
                             va = 0.48;
                             else if (t>=1977 && t<=1978)
                             va = 0.51;
                             else if (t>=1978 && t<=1979)
                             va = 0.53;
                             else if (t>=1979 && t<=1980)
                             va = 0.55;
                             else if (t>=1980 && t<=1981)
                             va = 0.58;
                             else if (t>=1981 && t<=1982)
                             va = 0.60;
                             else if (t>=1982 && t<=1983)
                             va = 0.63;
                             else if (t>=1983 && t<=1984)
                             va = 0.68;
                             else if (t>=1984 && t<=1985)
                             va = 0.71;
                             else if (t>=1985 && t<=1988)
                             va = 0.76;
                             else if (t>=1988 && t<=1989)
                             va = 0.814;
                             else if (t>=1989 && t<=1990)
                             va = 0.9488;
                             else if (t>=1990 && t<=1991)
                             va = 0.9818;
                             else if (t>=1991 && t<=1992)
                             va = 0.90;
                             else if (t>=1992 && t<=1993)
                             va = 0.92;
                             else if (t>=1993 && t<=1994)
                             va = 0.91;
                             else if (t>=1994 && t<=1995)
                             va = 0.91;
                             else if (t>=1995 && t<=1996)
                             va = 0.92;
                             else if (t>=1996 && t<=1997)
                             va = 0.92;
                             else if (t>=1997 && t<=1998)
                             va = 0.91;
                             else if (t>=1998 && t<=1999)
                             va = 0.88;
                             else if (t>=1999 && t<=2000)
                             va = 0.88;
                             else if (t>=2000 && t<=2001)
                             va = 0.87;
                             else if (t>=2001 && t<=2002)
                             va = 0.84;
                             else if (t>=2002 && t<=2003)
                             va = 0.82;
                             else if (t>=2003 && t<=2004)
                             va = 0.80;
                             else if (t>=2004 && t<=2005)
                             va = 0.81;
                             else if (t>=2005 && t<=2006)
                             va = 0.84;
                             else if (t>=2006 && t<=2007)
                             va = 0.85;
                             else if (t>=2007 && t<=2008)
                             va = 0.85;
                             else if (t>=2008 && t<=2009)
                             va = 0.85;
                             else if (t>=2009 && t<=2010)
                             va = 0.88;
                             else
                             va = 0.89;
                             
                             
                             DS =  birthrate*(1-va) -  Beta*S*I  - mu*S;
                             DE =  Beta*S*I - (sigma+mu)*E;
                             DI =  sigma*E - (gamma+mu)*I;
                             DR =  gamma*I  - mu*R + birthrate*va;
                             DH =  gamma*I;
                             ") #Seasonally forced snippet



# The above uses true incidence for reporting

initz <- Csnippet("
                  double m = pop/(S_0+E_0+I_0+R_0);
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  H = 0;
                  //printf(\"%f %f %f %f %f \\n\",S,E,R,I,H);
                  
                  ")
# Sampling from the normal approximation of the binomial distribution
#' Measurement process, calling on both rmeasure and dmeasure on H
#' The Data is  binomial distributed
dmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  
                  if (R_FINITE(cases)) {
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  
                  if (give_log) lik = log(lik);
                  } else {
                  lik = (give_log) ? 0 : 1;
                  }
                  
                  ")


rmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")




toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tpsi = log(psi);
                  Tgamma = log(gamma);
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tpsi = exp(psi);
                    Tgamma = exp(gamma);
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    ")

params_mod <- c("R0","amplitude","gamma","mu","sigma","rho","psi")

params_ic <- c("S_0","E_0","R_0","I_0")

statenames <- c("S","E","I","R","H")

zeronames <- c("H")
#########################################            data
setwd("~/Git/TrajMatch")
cureentfolder <-"~/Git/TrajMatch"
London_BiData <- read.csv(file.path(cureentfolder, "London_BiData.csv"))
London_covar <- read.csv(file.path(cureentfolder, "London_covar.csv"))
# covart <- covariate_table(London_covar,times="time")
#########################################       make pomp
pomp(
  data = London_BiData,
  times="time",
  t0=1940,#with(get(paste0(name,"_BiData")),2*time[1]-time[2]),
  skeleton=vectorfield(seasonal.sir.ode),
  rmeasure=rmeas,
  covar=London_covar, #covart,
  tcovar="time",
  dmeasure=dmeas,
  zeronames=zeronames,
  initializer=initz,
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  statenames=statenames,
  paramnames=c(params_mod,params_ic)
) -> m1


params_fixed <- "R0"
params_fixed <- c("gamma","sigma",params_fixed)
params_mod_fit <- c("amplitude","mu", "rho" ,"psi")
params_ic_fit <- character(0)

ParamSetFile <- paste0("ParamSet_DeterministicSEIR_run26.csv")
full_set <- read.csv(file=ParamSetFile,header=T)
current_params <- unlist(full_set[1,])

coef(m1) <- c(current_params)

traj.match(m1,
           
           start=current_params,
           transform=T,
           est=c()
) -> sets_traj
 coef(sets_traj)
 logLik(sets_traj)
 end_time <- Sys.time()
 end_time - start_time
 