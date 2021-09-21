
################ RESPONSE TO THE LOSS OF RESOURCES
GetLinearResponse <- function(frac,p){
  return(p + (1-p)*frac)
}

GetNonLinearResponse <- function(frac,p,a,b){
  return(p + (1-p)*(pbeta(frac, a, b))) 
}


############### BUILD TABLE FOR gRain
getNames <- function(prey){
  n <- length(prey)
  ch <- matrix("", 2^n, n);
  status <- rep("", 2^n);
  for(i in 0:((2^n)-1 )){
    x <- i;
    for(j in 1:n){
      if(x%%2==0){
        ch[i+1, j] <- "-"
        x <- floor(x/2);
      } else{
        ch[i+1, j] <- as.character(prey[j])
        x <- floor(x/2);
      }
      status[i+1] <- paste(status[i+1], " ", ch[i+1, j], sep='')
    }
  }
  return(status)
}

BuildTable <- function(M,Pi,model,alpha=1,beta=1){
  PiVector <- Pi
  S <- length(PiVector)
  Table <- list()
  Table$PiVector <- PiVector
  Table$M <- M;
  for (i in 1:S){
    variable_name <- paste("V",i,sep='')
    prey <- which(M[,i]!=0)
    n <- length(prey)
    if (n==0){
      ## Species is basal
      X <- matrix(0,1,1)
      X[1,1] <- PiVector[i]
      Table[[variable_name]] <- X
    }
    else{
      ## Species is consumer
      X<- matrix(0,2^sum(M[, i]>0),1)
      rownames(X) <- getNames(prey)
      if(model!="topo"){       
        for(j in 1:(2^n)){
          tmp <- j-1
          v <- numeric(n);
          for(k in 1:n){
            if((tmp %% 2) == 0){
              v[k]=0
              tmp <- floor(tmp/2)
            } else{
              v[k]=1
              tmp <- floor(tmp/2)
            }
          }
          ## Fraction resource lost. Works for binary and flow
          frac_loss <- 1-(M[M[, i]!=0, i] %*% v / (M[M[, i]!=0, i] %*% rep(1, n)))
          if(model=="linear"){
            X[j,1] <- GetLinearResponse (frac_loss, PiVector[i])
          } 
          if(model=="nonlinear"){
            X[j,1] <- GetNonLinearResponse(frac_loss, PiVector[i], alpha, beta)
          }	
        }
      } else if(model=="topo"){
        for(j in 1:(2^n)){
          tmp <- j-1
          v <- numeric(n)
          X[j, 1] <- 1
          for(k in 1:n){
            if(tmp %%2 == 0){
              v[k]=0;
              tmp <- floor(tmp/2);
            } else{
              v[k]=1
              X[j, 1] <- PiVector[i]
              break
            }
          } 
        }
      } else {
        stop("You must specify a response form: nonlinear, linear or topo(logical)")
      }	
      Table[[variable_name]] <- X
    }
  } 
  Table$model <- model
  return(Table)
}

############## Run gRain AND GET MARGINAL PROBS
getPreys <- function(v){
  n <- length(v)
  p <- sum(v>0)
  prey_str <- ''
  for(i in 1:n){
    if(p > 1){
      if(v[i]>0){
        prey_str <- paste(prey_str, 'V', i, "+", sep='')
        p <- p-1
      }
    } else{
      break
    }	
  }
  for(j in i:n){
    if(v[j]>0){
      prey_str <- paste(prey_str, 'V', j, sep='');
    }
  }
  return(prey_str)
}

getMarginals_gRain <- function(Table,file="test.R"){
  ## write down the script file
  write("require(gRain)", file=file);
  write("temp_fn <- function(){", file, append=TRUE)
  sp_levels <- "levels=c('extant', 'extinct'))"
  S <- length(Table$PiVector)
  for(i in 1:S){
    sp_name <- paste("V", i, sep='')
    if( dim(Table[[sp_name]])[1] == 1 ){ # Basal Species (Producers)
      X <- Table[[sp_name]]
      command <- paste("V", i, "<-cptable(~V", i, ", values=c(", 1-X[1,1],",", X[1,1], "), ",  sp_levels, sep="")
      write(command, file=file, append=TRUE)
    } else{ # Consumers 
      X <- Table[[sp_name]];
      nX <- dim(X)[1]
      command <- paste("V", i, "<-cptable(~V", i, "|", sep='')
      prey_str <- getPreys(Table$M[, i])
      command <- paste(command, prey_str, ", values=c(", sep='')
      for(j in 1:(nX-1)){
        command <- paste(command, 1-X[nX-j+1, 1], ",", X[nX-j+1, 1], ",", sep='')
      }
      command <- paste(command, 1-X[1, 1], ",", X[1, 1], "), ", sp_levels, sep='')
      write(command, file, append=TRUE)
    }
  }
  command2 <- paste("plist<-compileCPT(list(", sep='')
  for(i in 1:(S-1)){
    command2 <- paste(command2, "V", i, ",", sep='')
  }
  command2 <- paste(command2, "V", S, "))", sep='')
  write(command2, file, append=TRUE)
  write("BN <- grain(plist)", file, append=TRUE);
  write("Results_tmp <- querygrain(BN)", file, append=TRUE);
  write("return(Results_tmp)", file, append=TRUE);
  write("}", file, append=TRUE)
  
  ## Run the script
  source(file)
  r <- temp_fn()
  results <- list()
  ## Re-orgainze the sequence in the output
  for(i in 1:S){
    sp_name <- paste("V", i, sep='');
    results[[sp_name]] <- r[[sp_name]]
  }
  return(results)	
}
