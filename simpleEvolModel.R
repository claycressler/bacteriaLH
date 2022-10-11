 pick_individuals_multivariate <-function(N, traitmeans, traitsds, corr, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  if (any(traitmeans < 1e-12)) {
    traitmeans[which(traitmeans < 1e-12)] <- 1e-12
  }
  print(traitmeans) 
   
  p <- length(traitmeans);
  mu <- matrix(0,p,1);
  sigma <- matrix(0,p,p);
  ## Generate Lognormal random variables
  ## To produce a distribution with mean equal to traitmean
  ## and variance equal to traitsd^2, the parameters of
  ## the corresponding lognormal distribution are:
  for (i in 1:p) {
    mu[i] = log(traitmeans[i]^2 / sqrt(traitsds[i]^2+traitmeans[i]^2))
    sigma[i,i] = log(1 + traitsds[i]^2/traitmeans[i]^2);
  }
  ## generate the correlation matrix
  corr <- matrix(c(1,corr,corr,1),nrow=2,byrow=TRUE)
  ## generate the covariance matrix
  Cov <- sqrt(sigma)%*%corr%*%sqrt(sigma);
  Cov.Eigen = eigen(Cov);
  Cov.EigenVal <- Cov.Eigen$values*(abs(Cov.Eigen$values) >= 1e-12);
  RootCov = (Cov.Eigen$vectors)%*%diag(sqrt(Cov.EigenVal))%*%t(Cov.Eigen$vectors);
  X <- matrix(rnorm(p*N,0,1), nrow = N, ncol = p);
  for (i in 1:N){
    X[i,] <- exp(mu + RootCov%*%X[i,])
  }
  colnames(X) <- names(traitmeans) ## give the columns names
  return(X);
}

simulate_model <- function(params) {
  P0 = params[1] ## initial parasite abundance
  R0 = params[2] ## initial resource aubndance
  Rconc = params[3] ## influx concentration of resource
  Rin = params[4] ## influx rate
  Rout = params[5] ## outflux rate for use in the body for other purposes
  g = params[6] ## maximum per-capita growth rate of the parasite
  K = params[7] ## half-saturation constant of growth rate equation
  m = params[8] ## maintenance rate
  Y = params[9] ## yield
  gSD = params[10] ## standard deviation in per-capita growth rate
  mSD = params[11] ## standard deviation in maintenance
  corr = params[12] ## positive correlation (trade-off) between growth rate and maintenance
  tmax = params[13] ## maximum time to run
  Rinterval = params[14] ## how often to cycle food on and off
  test = params[15] ## whether to allow evolution or not
  
  ## Food will be influxing only in every other Rinterval, with it being on at t=0
  ## E.g. if Rinterval = 5, then food influxes from 0-5, 10-15, 20-25, etc.
  if (Rinterval!=0) { # Rinterval = 0 implies resources are always inflowing
    if ((length(seq(0,tmax,by=Rinterval))%%2)==0)
      influxIntervals = seq(0,tmax,by=Rinterval) %>% matrix(.,ncol=2,byrow=TRUE)
    else 
      influxIntervals = seq(0,tmax-Rinterval,by=Rinterval) %>% matrix(.,ncol=2,byrow=TRUE)
  }
  
  ## Set up storage of parasite traits and success
  storage <- data.frame(ID=1:1e5, g=rep(NA,1e5), m=rep(NA,1e5), numRep=rep(0,1e5), tBorn=rep(NA,1e5), tDie=rep(Inf,1e5))
  ## Set the traits of the first P0 individuals in the population
  if (test) {
    storage[1:P0,"g"] <- g
    storage[1:P0,"m"] <- m
  }
  else {
    new <- pick_individuals_multivariate(P0, c(g,m), c(gSD,mSD), corr)
    storage[1:P0,"g"] <- new[,1]
    storage[1:P0,"m"] <- new[,2]
  }
  storage[1:P0,'tBorn'] <- 0
  ## popn-level stat storage
  stats <- data.frame(t=rep(NA,1e5), R=rep(NA,1e5), P=rep(NA,1e5), gmean=rep(NA,1e5), mmean=rep(NA,1e5))
  
  t = 0 
  i = 1
  R = R0
  
  while (t < tmax & P > 0) {
    #print(t)
    popn = filter(storage, tBorn <= t & is.infinite(tDie))
    P = nrow(popn)
    ## store the current population-level stats
    stats[i,] = c(t, R, P, mean(popn$g), mean(popn$m)); i = i + 1
    
    ## compute rates
    ## are resources inflowing?
    if (Rinterval!=0) {
      if (any(apply(influxIntervals, 1, function(int) (t >= int[1] & t < int[2]))))
        Rinflow = Rin*Rconc
      else Rinflow = 0
    }
    else Rinflow = Rin*Rconc
    Routflow = Rout*R
    Rconsump = sum(popn$g/Y * R/(K + R))
    Pgrowth = popn$g * R/(K+R) 
    Pdeath = popn$m 
    ## wheel of fortune
    rates = c(Pgrowth, Pdeath, Rinflow, Routflow+Rconsump)
    wheel = cumsum(rates)/sum(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    ## update t 
    t <- t + dt
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:P) {## first P events are reproductions
      ## get the ID of this new individual
      id = min(which(is.na(storage$tBorn)))
      if (test) storage[id,c('g','m')] <- c(g,m)
      else 
        ## set its traits based on those of its "parent"
        storage[id,c('g','m')] = pick_individuals_multivariate(1, c(popn$g[event],popn$m[event]), c(gSD,mSD), corr)
      storage[id,'tBorn'] = t
      ## update the fitness of the individual that reproduced
      storage[popn$ID[event],'numRep'] = storage[popn$ID[event],'numRep'] + 1
    } 
    else if (event%in%((P+1):(2*P))) { ## next P events are deaths
      ## get the ID of the individual who dies and set its time of death
      storage$tDie[popn$ID[event-P]] = t
    }
    else if (event==(2*P+1))  ## increment resources 
      R = R + 1
    else R = R - 1 ## decrement resources
  }
  ## drop all the missing information from storage and popn
  stats = stats[1:(i-1),]
  storage = filter(storage, !is.na(tBorn))
  return(list(stats,storage))
}

params0 = c(P0 = 20, R0 = 50, Rconc = 50, Rin=1, Rout=1, g=1, K=10, m=0.5, Y=2, gSD=0.1, mSD=0.05, corr=0.75, tmax=50, Rinterval=0, test=TRUE)
sim0 = simulate_model(params0)
mean(tail(sim0[[1]],4000)[,2]) ## should be 10
mean(tail(sim0[[1]],4000)[,3]) ## should be 160

params0 = c(P0 = 20, R0 = 50, Rconc = 50, Rin=1, Rout=1, g=1, K=10, m=0.5, Y=2, gSD=0.1, mSD=0.05, corr=0.75, tmax=50, Rinterval=0, test=FALSE)
sim0 = simulate_model(params0)
plot(sim0[[1]]$gmean, sim0[[1]]$mmean)

plot(sim0[[1]][,c('t','gmean')], type='l')
plot(sim0[[1]][,c('t','mmean')], type='l')

plot(sim0[[1]][,c('t','R')], type='l')
plot(sim0[[1]][,c('t','P')], type='l')


params1 = c(P0 = 20, R0 = 50, Rconc = 50, Rin=1, Rout=1, g=1, K=10, m=0.5, Y=2, gSD=0.1, mSD=0.05, corr=0.75, tmax=100, Rinterval=5, test=FALSE)
sim1 = simulate_model(params1)
plot(sim1[[1]][,c('t','gmean')], type='l')
plot(sim1[[1]][,c('t','mmean')], type='l')
plot(sim1[[1]][,c('t','R')], type='l')
plot(sim1[[1]][,c('t','P')], type='l')

plot(sim1[[1]]$gmean, sim1[[1]]$mmean, type='l')

   
params0 = c(P0 = 20, R0 = 50, Rconc = 50, Rin=1, Rout=1, 
            g=1, K=10, m=0.5, Y=2, gSD=0.1, mSD=0.05, corr=0.95, 
            tmax=200, Rinterval=0, test=FALSE)
sim0 = simulate_model(params0)

params1 = c(P0 = 20, R0 = 50, Rconc = 50, Rin=1, Rout=1, 
            g=1, K=10, m=0.5, Y=2, gSD=0.1, mSD=0.05, corr=0.95, 
            tmax=200, Rinterval=10, test=FALSE)
sim1 = simulate_model(params1)


    

