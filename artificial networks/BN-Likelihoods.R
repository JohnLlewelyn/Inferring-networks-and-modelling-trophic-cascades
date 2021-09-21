## Get likelihood given observations of the extant species
GetLogLikelihoodExtantGivenMarginals <- function(Extant,MarginalProbs){
  NReps <- dim(Extant)[1]
  SumExtant <- apply(Extant,2,sum)
  LogLikeli <- SumExtant*log(MarginalProbs)+(NReps-SumExtant)*log(1.0-MarginalProbs)
  LogLikeli[is.nan(LogLikeli)] <- 0
  return(sum(LogLikeli))
}

## Maximum Log Likelihood
GetMaximumLogLikelihood <-function(Extant){
  cs <- apply(Extant,2,sum)
  NReps <- dim(Extant)[1]
  phat<-cs/NReps
  return(GetLogLikelihoodExtantGivenMarginals(Extant,phat))
}

GetBayesNetMarginals <-function(M,Extant,PiVector,Type="linear",alpha=10,beta=10,Label){
  ## Read the files
  S <- dim(M)[1]
  ## Remove redundant
  Functional <- M
  ## Functional <- Functional2(RootFW(M))
  ## Remove the root node as it is no longer needed
  ## Functional <- Functional[-1,-1]
  BayesNet <- list(M=M,Functional=Functional,PiVector=PiVector,Extant=Extant)
  ## Build gRain table
  MyTable <- BuildTable(Functional,PiVector,Type,alpha,beta)
  ##DBGprint(MyTable,Label,2)
  ## Get Marginals
  ResultsgRain <- getMarginals_gRain(MyTable, file=paste("gRainTable-",Label,".R",sep=""))
  ##DBGprint(ResultsgRain,Label,2)
  ## MarginalExtant
  MarginalExtant <- rep(0,S)
  for (i in 1:S){
    MarginalExtant[i] <- ResultsgRain[[i]][1]
  }
  BayesNet$MarginalsExtant<-MarginalExtant
  BayesNet$Likelihood<-GetLogLikelihoodExtantGivenMarginals(Extant,MarginalExtant)
  BayesNet$ML<-GetMaximumLogLikelihood(Extant)
  BayesNet$BetaDistralpha <- alpha
  BayesNet$BetaDistrbeta <- beta
  return(BayesNet)
}

InterfaceForNelderMead <- function(Pars,M,Extant,PiVector,Type="nonlinear",Label){
  alpha <- exp(Pars[1])
  beta <- exp(Pars[2])
  return(-GetBayesNetMarginals(M,Extant,PiVector,Type,alpha,beta,Label)$Likelihood)
}

FindMLAlphaAndBeta <- function(M,Extant,PiVector,Type="nonlinear",Label, HowManyPoints=5){
  ## First starting point
  Pars <- log(c(0.01,0.01))
  Current <- optim(Pars,InterfaceForNelderMead,M=M,Extant=Extant,PiVector=PiVector,Type=Type,Label=Label)
  CurrentLik <- Current$value
  CurrentPars <- exp(Current$par)
  BestPars <- CurrentPars
  BestLik <- CurrentLik
  print(c(BestLik,0,BestPars))
  ## Second starting point
  Pars <- log(c(1.0,1.0))
  Current <- optim(Pars,InterfaceForNelderMead,M=M,Extant=Extant,PiVector=PiVector,Type=Type,Label=Label)
  CurrentLik <- Current$value
  CurrentPars <- exp(Current$par)
  if (CurrentLik<BestLik){
    BestPars <- CurrentPars
    BestLik <- CurrentLik
    print(c(BestLik,1,BestPars))
  }
  for (i in 1:HowManyPoints){
    Pars <- rnorm(2,0,5)
    Current <-  optim(Pars,InterfaceForNelderMead,M=M,Extant=Extant,PiVector=PiVector,Type=Type,Label=Label)
    CurrentLik <- Current$value
    CurrentPars <- exp(Current$par)
    if (CurrentLik<BestLik){
      BestPars <- CurrentPars
      BestLik <- CurrentLik
      print(c(BestLik,1+i,BestPars))
    }
  }
  return(BestPars)
}
