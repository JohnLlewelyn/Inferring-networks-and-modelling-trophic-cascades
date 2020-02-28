#Infer links for predators of vertebrates
#Add links for invertivores
#Add links for herbivores 
#Remove unrealistic links for large birds (above 5 kg)

#Libraries and wd
library(igraph)
library(gtools)
library(adegenet)
library(plyr)
library(quantreg)
library(dplyr)

#save file path
pp <- getwd()

#get data
keep80<-readRDS(paste(pp,"/data/keep80.rds",sep=""))

#remove megafauna that were unlikely to have been there, and fix some names
source(paste(pp,"/scripts/megafauna vetting.R", sep=""))

#Get quantile regressions
qrup<-readRDS(paste(pp,"/data/qrup.rds", sep="")) #generated in Trophic niche space validation.R script
qrlow<-readRDS(paste(pp,"/data/qrlow.rds", sep=""))

##Function for making edge list from data frames structured as keep80.rds
make.edge<-function(keep){
  #list of species that prey on vertebrates
  vpred<-keep[keep$eat.vert==TRUE,]
  
  #get all the species permutations
  options(expressions=100000)
  perms<-data.frame(permutations(n=length(keep$SciName),r=2,v=keep$SciName,repeats=TRUE))
  colnames(perms)<-c("from","to")
  
  #merge diet and body mass data
  perms<-merge(perms,keep[,c("SciName","mean.mass.kg","eat.plant","eat.invert","eat.vert")],by.x="from",by.y="SciName",all.x=T)
  colnames(perms)<-c("from","to","fmasskg","fplant","finvert","fvert")
  perms<-merge(perms,keep[,c("SciName","mean.mass.kg","eat.plant","eat.invert","eat.vert","Class")],by.x="to",by.y="SciName",all.x=T)
  colnames(perms)<-c("to","from","fmasskg","fplant","finvert","fvert","tmasskg","tplant","tinvert","tvert","Cpred")
  
  vpred<-perms[perms$tvert==TRUE,]  #just predators of vertebrates in the 'to' column
  vpred$Bprey<-log10(vpred$fmasskg*1000)  #log10 of mass, in g rather than kg
  vpred$Bpred<-log10(vpred$tmasskg*1000)
  vpred$Cpred<-as.character(vpred$Cpred)  #Fix class name and categories
  vpred$Cpred[vpred$Cpred=="Aves"]<-"aves"
  vpred$Cpred[vpred$Cpred=="Mammalia"]<-"mammals"
  vpred$Cpred[vpred$Cpred=="Amphibia"|vpred$Cpred=="Reptilia"]<-"amphibia and reptilia"
  
  #predict upper and lower limits
  vpred$up<-predict.rq(qrup,vpred)
  vpred$low<-predict.rq(qrlow,vpred)
  
  #then just take predicted links
  vpred<-subset(vpred, vpred$Bprey>vpred$low & vpred$Bprey<vpred$up)
  plot(vpred$Bpred,vpred$Bprey)
  #make dataframe of number of prey predicted per predator
  np<-data.frame(predator=as.character(data.frame(table(vpred$to))[,1]),number.prey.predicted=as.numeric(data.frame(table(vpred$to))[,2]),stringsAsFactors = FALSE)
  #get rid of unrealistic links for some of the bigger birds
  vpred<-vpred[vpred$to!="Grus rubicunda"|vpred$from%in%keep$SciName[keep$mean.mass.kg<(keep$mean.mass.kg[keep$SciName=="Grus rubicunda"])/3],]
  vpred<-vpred[vpred$to!="Ardeotis australis"|vpred$from%in%keep$SciName[keep$mean.mass.kg<(keep$mean.mass.kg[keep$SciName=="Ardeotis australis"])/3],]
  vpred<-vpred[vpred$to!="Threskiornis spinicollis"|vpred$from%in%keep$SciName[keep$mean.mass.kg<(keep$mean.mass.kg[keep$SciName=="Threskiornis spinicollis"])/3],]
  vpred<-vpred[vpred$to!="Pelecanus conspicillatus"|vpred$from%in%keep$SciName[keep$mean.mass.kg<(keep$mean.mass.kg[keep$SciName=="Pelecanus conspicillatus"])/3],]
  vpred$to<-as.character(vpred$to)
  vpred<-vpred[!duplicated(vpred[c("from","to")]),]
  
  #add probability, based on normal distribution
  vpred$prob<-dnorm(vpred$Bprey,mean=(vpred$up+vpred$low)/2,sd = (vpred$up - ((vpred$up+vpred$low)/2))/2  )
  preds<-vpred[,c("from","to","prob")]
  
  #now get invertivore links
  inv<-data.frame(cbind("invertebrates", keep$SciName[keep$eat.invert==TRUE]))
  inv$prob<-1
  colnames(inv)<-c("from","to","prob")
  
  #now get herbivore links
  herb<-data.frame(cbind("plants", keep$SciName[keep$eat.plant==TRUE]))
  colnames(herb)<-c("from","to")
  herb$prob<-1
  
  #now get links to fish and aquatic invertebrates
  pisc<-data.frame(cbind("FishInv", keep$SciName[keep$eat.FishInv==TRUE]))
  colnames(pisc)<-c("from","to")
  pisc$prob<-1
  
  #Stick them all together to make an edge list
  links<-rbind(preds,inv,herb,pisc)
  
  #node list
  nl<-unique(c(links$from,links$to))
  return(list(links,np))
}

#apply function to species list
ik80<-make.edge(keep80)

#split edge lists from number of prey predicted
num80<-data.frame(ik80[[2]],stringsAsFactors = F)#data.frame(predator=as.character(ik80[[2]][,2]), number.prey.predicted=as.numeric(ik80[[2]][,1]),stringsAsFactors = FALSE)
ed80<-data.frame(ik80[[1]],stringsAsFactors = FALSE)

#########################################################
#Build networks##########################################
#For each predator
#1) randomly select how many links to keep/remove
#2) round off
#3) randomly keep/remove that many links, taking into account size-based probability
#4) add links to plants, invertebrates and aquatic animas
#5) save network and make the next

#get Serengeti overestimation distribution (from de Visser's data) to randomly select how many to delete for each predator
#The Serengeti food web: empirical quantification and analysis of topological changes under increasing human impact.
#de Visser SN1, Freymann BP, Olff H. 2011
#Also de Visser unpublished + relevant GloBI records
d<-readRDS(paste(pp,"/data/Serengeti.density.rds",sep=""))
pred.ob<-readRDS(paste(pp,"/data/pred.ob.rds", sep=""))

#get herbivore/plant density data
#Baskerville, Edward B., Andy P. Dobson, Trevor Bedford, Stefano Allesina, T. Michael Anderson, and Mercedes Pascual. “Spatial Guilds in the Serengeti Food Web Revealed by a Bayesian Group Model.” Edited by Lauren Ancel Meyers. PLoS Computational Biology 7, no. 12 (December 29, 2011): e1002321. https://doi.org/10.1371/journal.pcbi.1002321.
#rds file from Herbivores and plant diversity
bask<-readRDS(paste(pp,"/data/mammal.herb.density.rds", sep=""))
bask.d<-bask[[1]]
bask.r<-bask[[2]]

#Get distributions for invertivore and piscivore diet breadth
#take the mean and sd and turn into normal distribution (because only got a small sampel size)
diet<-read.csv(paste(pp,"/tables/Table S5 diet breadth.csv",sep=""))
samp<-seq(0, 300, length=1000)
invs<-dnorm(samp,mean(diet$breadth[diet$diet.group=="inv"]),sd(diet$breadth[diet$diet.group=="inv"]))
piscs<-dnorm(samp,mean(diet$breadth[diet$diet.group=="pisc"]),sd(diet$breadth[diet$diet.group=="pisc"])) 

#Function to calculate how many links to keep for predator
lkp<-function(y) {                 #y = the number of predicted links before unrealistic ones were removed
  lkp<-sample(pred.ob$prop, dim(y)[1], replace=TRUE) + rnorm(dim(y)[1], 0, d$bw) #sample overestimations from Serengeti data
  lkp<-ifelse(lkp<0.023,0.023,lkp) #lower bound - must have > 2.3% of predicted links, because all theSerengeti predators did
  lkp<-ifelse(lkp>0.63,0.63,lkp)   #upper bound - must have less than 63% of links, because all theSerengeti predators did
  lkp<-lkp*y[,2]
  lkp<-ifelse(lkp<1,1,round(lkp))
  lkp<-data.frame(pred=y[,1],links=lkp,stringsAsFactors = FALSE)
  lkp$pred<-as.character(lkp$pred)
  return(lkp)
}

#Function to select # of rows based on links column and probability column
rsel<-function(x){
  x[sample(nrow(x),mean(x$links),prob = x$prob,replace = F),]
}

#function to randomly select number of rows specified by lkp with probability taken into account
lkp.ap<-function(y, net) #y = number predicted; net = edge list  
{     
  oth<-net[net$from=="plants"|net$from=="invertebrates"|net$from=="FishInv",]  #separate consumption of plants, inverts, fish and inverts from predation on vertebrates
  pred<-net[net$from!="plants"&net$from!="invertebrates"&net$from!="FishInv",] #predation on vertebrates
  pred$to<-as.character(pred$to)
  pred$from<-as.character(pred$from)
  kp<-lkp(y)                                                                   #Get how many links to keep
  pred<-merge(pred,kp,by.x="to",by.y="pred",all.x=T)                           #add column for number of links to keep for each predator
  out <- split( pred , f = pred$to)                                            #Split data frame
  rand<-lapply(out,rsel)                                                       #randomly select rows for each df, taking into account probability (based on body size)
  rand<-ldply(rand, rbind)                                                     #put dfs in rand together as one df
  pls<-oth[oth$from=="plants",]                                                #Add links for herbivores, invertivores and piscivores
  inv<-oth[oth$from=="invertebrates",]
  fishi<-oth[oth$from=="FishInv",]
  plinks<-round(sample(bask.r$diversity,dim(pls)[1], replace=TRUE) + rnorm(dim(pls)[1], 0, bask.d$bw))
  plinks<-ifelse(plinks<3,3,ifelse(plinks>73,73,plinks))
  pls$links<-plinks
  pls<-data.frame(pls[rep(seq_len(dim(pls)[1]), pls$links), ], row.names=NULL)
  pls<-pls[,c("to","from")]
  pls<-pls[pls$to!="invertebrates"&pls$to!="FishInv",]                          #get rid of links from plant nodes to the lumped invertebrate and FishInv nodes
  
  ilinks<-round(sample(samp,dim(inv)[1],replace=T,prob = invs))                 #Add invertebrate links; mean and sd of links based on studies
  ilinks<-ifelse(ilinks<3,3,ifelse(ilinks>300,300,ilinks))                      #make upper limit double limit from records
  inv$links<-ilinks
  inv<-data.frame(inv[rep(seq_len(dim(inv)[1]), inv$links), ], row.names=NULL)
  inv<-inv[,c("to","from")]
  
  flinks<-round(sample(samp,dim(fishi)[1],replace=T,prob = piscs))                 #Add FishInv links; diversity based on Bloon lagoon report; mean and sd of links based on studies
  flinks<-ifelse(flinks<3,3,ifelse(flinks>300,300,flinks))
  fishi$links<-flinks
  fishi<-data.frame(fishi[rep(seq_len(dim(fishi)[1]), fishi$links), ], row.names=NULL)
  fishi<-fishi[,c("to","from")]
  fishi<-fishi[fishi$to!="invertebrates"&fishi$to!="FishInv",]                     #get rid of links from plant nodes to the lumped invertebrate and FishInv nodes
  
  
  four<-keep80$SciName[keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="TRUE"]           #omnivores, adjust number of links because feed from more than one category
  three<-unique(keep80$SciName[keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="FALSE"|
                              keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.vert=="TRUE"|
                              keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="TRUE"|
                              keep80$eat.FishInv=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="TRUE"])
  two<-unique(keep80$SciName[  keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.vert=="FALSE"|
                                keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="FALSE"&keep80$eat.vert=="TRUE"|
                                keep80$eat.FishInv=="FALSE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="TRUE"|
                                keep80$eat.FishInv=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="FALSE"|
                                keep80$eat.FishInv=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.vert=="TRUE"|
                                keep80$eat.FishInv=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.vert=="FALSE"])
  omnis<-c(four,three,two)                                                                               #species that feed from more than one group
  rand<-rand[,c("to","from")]
  four<-list(rand[rand$to%in%four,],pls[pls$to%in%four,],inv[inv$to%in%four,],fishi[fishi$to%in%four,])  #split animals depending on how many groups they feed from
  three<-list(rand[rand$to%in%three,],pls[pls$to%in%three,],inv[inv$to%in%three,],fishi[fishi$to%in%three,])
  two<-list(rand[rand$to%in%two,],pls[pls$to%in%two,],inv[inv$to%in%two,],fishi[fishi$to%in%two,])
  one<-list(rand[!(rand$to%in%omnis),],pls[!(pls$to%in%omnis),],inv[!(inv$to%in%omnis),],fishi[!(fishi$to%in%omnis),])
  
  fourfun<-function(y){ddply(y,.(to),function(x) x[sample(nrow(x),round(0.25*nrow(x)),replace=FALSE),])}# Select rows by category, do it for the 4 lists
  threefun<-function(y){ddply(y,.(to),function(x) x[sample(nrow(x),round((1/3)*nrow(x)),replace=FALSE),])}
  twofun<-function(y){ddply(y,.(to),function(x) x[sample(nrow(x),round(0.5*nrow(x)),replace=FALSE),])}
  four<-do.call("rbind",(lapply(four,fourfun))) #now, apply to list and stack them 
  three<-do.call("rbind",(lapply(three,threefun)))
  two<-do.call("rbind",(lapply(two,twofun)))
  one<-do.call("rbind",one)
  nodes<-data.frame(to=c("invertebrates","invertebrates","invertebrates","FishInv","FishInv","FishInv"),from=c("plants","invertebrates","FishInv","plants","invertebrates","FishInv")) #add links between lumped nodes at the end
  rand<-do.call("rbind",list(one,two,three,four,nodes))
  
  return(data.frame(rand))
}

#Now apply it to assemblage 1000 network dataframes
set.seed(10)
x80 <- t(replicate(1000, lkp.ap(num80,ed80),simplify = FALSE))
saveRDS(x80,paste(pp,"/networks/x80nets.rds",sep=""))


