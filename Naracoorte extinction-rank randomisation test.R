## Network extinction-rank randomisation test
## CJA Bradshaw
## March 2020

## remove everything
rm(list = ls())

# set working directory
#setwd("~/Documents/Papers/Palaeo/Sahul/Networks/Naracoorte co-extinction/data/")

# load network from file
dat <- read.table("extinction_rank.csv", header=T, sep=",")
head(dat)
summary(dat)

# choose random sample of extinct species
#ext.spp1 <- subset(dat, extinction_status == 1 & replicate == 1) # number of extinct species in set - but differs between replicates
#lext.spp <- dim(ext.spp1)[1]

iter <- 10000
itdiv <- iter/100
rep.rand <- sample(min(dat$replicate):max(dat$replicate), size=iter, replace=T) # random replicate selection vector
obs.ran.diff <- lte <- rep(NA,iter)

for (i in 1:iter) {
  ext.spp.rep <- subset(dat, extinction_status == 1 & replicate == rep.rand[i])           #get extinct species in replicate [i]
  ext.spp.ran <- ext.spp.rep[sample(1:dim(ext.spp.rep)[1], size=dim(ext.spp.rep)[1], replace=T),]    #From the extinct subset, sample (with replacement) rows as many times are there are rows #
  extn.spp.rep <- subset(dat, extinction_status == 0 & replicate == rep.rand[i])          #get extant species in replicate [i]
  extn.spp.ran <- extn.spp.rep[sample(1:dim(extn.spp.rep)[1], size=dim(ext.spp.rep)[1], replace=T),] #From the extant subset, sample (with replacement) rows as many times are there are rows #
  ext.extn.diff <- mean(ext.spp.ran$extinction_rank) - mean(extn.spp.ran$extinction_rank) #difference in mean extinction ranks between the two sampled groups
  
  ran.spp.rep <- subset(dat, replicate == rep.rand[i])                                    #take one replicate
  ran.spp.ran <- ran.spp.rep[sample(1:dim(ran.spp.rep)[1], size=2*dim(ext.spp.rep)[1], replace=T),]  #for that replicate, sample twice as many rows as there are extinct species
  ran.spp.ran1 <- ran.spp.ran[1:dim(ext.spp.rep)[1],]                                                #subset, take the first half...
  ran.spp.ran2 <- ran.spp.ran[(dim(ext.spp.rep)[1]+1):(2*dim(ext.spp.rep)[1]),]                                 #and the second half of the data set
  ran.diff <- mean(ran.spp.ran1$extinction_rank) - mean(ran.spp.ran2$extinction_rank)     #calculate the difference between teh two halves

  obs.ran.diff[i] <- ext.extn.diff - ran.diff                                             #calculate difference between ext/extant versus difference between random samples
  lte[i] <- ifelse(obs.ran.diff[i] <= 0, 1, 0)                                            #if extinct species don't have higher rank than extant in the comparison, then 1, otherwise 0.
  if (i %% itdiv==0) print(i)                                                             #track where comparisons are up to
}

pr.lte <- sum(lte, na.rm=T)/(iter - length(which(is.na(lte)==T)))                         #number of times extinct species didn't have higher extinction ranks divided by number of iterations
pr.lte
hist(obs.ran.diff, col="black", border="white")
abline(v=mean(obs.ran.diff, na.rm=T), col="red", lty=2, lwd=3)
abline(v=quantile(obs.ran.diff, probs=0.025, na.rm=T), col="grey", lty=2, lwd=2)
abline(v=quantile(obs.ran.diff, probs=0.975, na.rm=T), col="grey", lty=2, lwd=2)
