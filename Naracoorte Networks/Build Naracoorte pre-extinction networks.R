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
#setwd("~/###")
pp <- getwd() 

#get data
keep80<-readRDS(paste(pp,"/data/keep80.rds",sep=""))

#remove megafauna that were unlikely to have been there, fix some names, update thylacine traits
source(paste(pp,"/script/megafauna vetting and fix thylacine.R", sep=""))

#fix some other stuff in keep80 (diet categories, some diet info and body mass)
source(paste(pp,"/script/fix diet categories.R", sep=""))


#Get quantile regressions
qrup<-readRDS(paste(pp,"/data/qrup.rds", sep="")) #generated in Trophic niche space validation.R script
qrlow<-readRDS(paste(pp,"/data/qrlow.rds", sep=""))

##Function for making edge list from data frames structured as keep80.rds
make.edge<-function(keep){
  #list of species that prey on vertebrates
  vpred<-keep[keep$eat.terr_vert==TRUE,]
  
  #get all the species permutations
  options(expressions=100000)
  perms<-data.frame(permutations(n=length(keep$SciName),r=2,v=keep$SciName,repeats=TRUE))
  colnames(perms)<-c("from","to")
  
  #merge diet and body mass data
  perms<-merge(perms,keep[,c("SciName","mean.mass.kg","eat.plant","eat.invert","eat.terr_vert","eat.fish")],by.x="from",by.y="SciName",all.x=T)
  colnames(perms)<-c("from","to","fmasskg","fplant","finvert","fvert","ffish")
  perms<-merge(perms,keep[,c("SciName","mean.mass.kg","eat.plant","eat.invert","eat.terr_vert","eat.fish","Class")],by.x="to",by.y="SciName",all.x=T)
  colnames(perms)<-c("to","from","fmasskg","fplant","finvert","fvert","ffish","tmasskg","tplant","tinvert","tvert","tfish","Cpred")
  
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
  #plot it to check
  library(RColorBrewer)
  cols = brewer.pal(4, "Blues")
  pal = colorRampPalette(c("blue", "red"))
  pal = colorRampPalette(cols)
  ord <- findInterval(vpred$prob,sort(vpred$prob))
  plot(vpred$Bpred,vpred$Bprey, pch=19, col=pal(nrow(vpred))[ord])
  preds<-vpred[,c("from","to","prob")]
  
  #now get invertivore links
  inv<-data.frame(cbind("invertebrates", keep$SciName[keep$eat.invert==TRUE]))
  inv$prob<-1
  colnames(inv)<-c("from","to","prob")
  
  #now get herbivore links
  herb<-data.frame(cbind("plants", keep$SciName[keep$eat.plant==TRUE]))
  colnames(herb)<-c("from","to")
  herb$prob<-1
  
  #now get links to fish 
  pisc<-data.frame(cbind("fish", keep$SciName[keep$eat.fish==TRUE]))
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
num80<-data.frame(ik80[[2]],stringsAsFactors = F) #data.frame(predator=as.character(ik80[[2]][,2]), number.prey.predicted=as.numeric(ik80[[2]][,1]),stringsAsFactors = FALSE)
ed80<-data.frame(ik80[[1]],stringsAsFactors = FALSE)

#plot number of suitably sized prey versus predator mass
pm <- merge(num80,keep80[,c("SciName","mean.mass.kg")],by.x="predator",by.y="SciName",all.x=T)

pdf("###.pdf",height = 6, width = 8)
plot(log10(pm$mean.mass.kg),pm$number.prey.predicted,pch=19,xlab="log10 predator mass",ylab="# potential prey")
dev.off()
#########################################################
#Build networks##########################################
#For each predator
#1) randomly sample from carnivore diet breadth (literature) the number of links for each predator #the proportion of links to keep/remove
#2) randomly keep/remove that many links, 
#3) add links to plants, invertebrates and aquatic animals in the same way
#4) save network and make the next

#get Serengeti overestimation distribution (from de Visser's data) to randomly select how many to delete for each predator
#The Serengeti food web: empirical quantification and analysis of topological changes under increasing human impact.
#de Visser SN1, Freymann BP, Olff H. 2011
#Also de Visser unpublished + relevant GloBI records
#d<-readRDS(paste(pp,"/data/Serengeti.density.rds",sep=""))
#pred.ob<-readRDS(paste(pp,"/data/pred.ob.rds", sep=""))

#get herbivore/plant density data - see Table S5 in Supplementary Material
diet <- read.csv(paste(pp,"/data/Table S5 diet breadth.csv", sep=""))
#don't include the DNA metabarcoding studies
diet <- diet[diet$DNA_study!="1",]

#fix diet group based on 'eat' columns
diet$diet.group <- paste(diet$eat.plants,diet$eat.invert,diet$eat.Fish,diet$eat.vert, sep ="")
diet$diet.group <- ifelse(diet$diet.group=="FALSEFALSEFALSETRUE", "carn",
                        ifelse(diet$diet.group=="FALSETRUEFALSEFALSE","invert",
                               ifelse(diet$diet.group=="TRUEFALSEFALSEFALSE","herb",
                                      ifelse(diet$diet.group=="TRUETRUEFALSETRUE","herb-invert-vert",
                                             ifelse(diet$diet.group=="FALSEFALSETRUEFALSE","pisc",
                                                    ifelse(diet$diet.group=="FALSEFALSETRUETRUE","pisc-vert",
                                                           ifelse(diet$diet.group=="FALSETRUETRUEFALSE","inv-pisc",
                                                                  ifelse(diet$diet.group=="FALSETRUETRUETRUE","inv-pisc-vert",
                                                                         ifelse(diet$diet.group=="TRUEFALSEFALSETRUE","herb-vert",
                                                                                ifelse(diet$diet.group=="TRUETRUEFALSEFALSE","herb-inv",
                                                                                       ifelse(diet$diet.group=="TRUETRUEFALSETRUE","herb-inv-vert",
                                                                                              ifelse(diet$diet.group=="FALSETRUEFALSETRUE","inv-vert","all"))))))))))))

herb_b <- diet[diet$diet.group=="herb",]
hd <- density(herb_b$breadth, adjust = 2)
inv_b <- diet[diet$diet.group=="invert",]
id <- density(inv_b$breadth, adjust = 2)
carn_b <- diet[diet$diet.group=="carn",]
cd <- density(carn_b$breadth, adjust = 2)
#get piscivore data
fi<-read.csv(paste(pp,"/data/Table S5b fish diversity in Australian inland piscivores.csv",sep=""))
fi$numb_cat <- apply(fi[,c("eat.plants","eat.invert","eat.vert","eat.Fish")],1,table)[2,]
fish_breadth <- fi$fish_div*fi$numb_cat
fd <- density(fish_breadth, adjust = 2)

#Function to calculate how many links to keep for predator
lkp<-function(y) {                 #y = the number of predicted links before adjustment to get realistic number of links
  lkp<-rpois(length(y),lambda=y)
  return(lkp)
}

#Function to select # of rows based on links column and probability column
rsel<-function(x){
  x[sample(nrow(x),mean(x$links),prob = x$prob,replace = F),]
}

#function to randomly select number of rows specified by lkp with probability taken into account
lkp.ap<-function(y, net) #y = number prey predicted; net = edge list  ; could add omi_adj = "yes" or "no" indicating if adjustment should be made to reduce number of links an omnivore has
{     
  oth<-net[net$from=="plants"|net$from=="invertebrates"|net$from=="fish",]  #separate consumption of plants, inverts, fish and inverts from predation on vertebrates
  pred<-net[net$from!="plants"&net$from!="invertebrates"&net$from!="fish",] #predation on vertebrates
  pred$to<-as.character(pred$to)
  pred$from<-as.character(pred$from)
  kp<-cbind.data.frame(y$predator, lkp(y$number.prey.predicted))                                                                  #Get how many links to keep
  names(kp) <- c("pred","pois_prey_num")                                    #small shuffle of number of potential prey by sampling from Poisson distribution
  kp <- kp[order(kp$pois_prey_num),]
  clinks <- round(sample(carn_b$breadth, 2*length(unique(pred$to)), replace=TRUE) + rnorm(2*length(unique(pred$to)), 0, hd$bw)) #sample from distribution, take extra samples so can remove ones outside the realistic range
  clinks <- clinks[!(clinks<1)&!(clinks>2*max(carn_b$breadth))]                  #get rid of unrealistic numbers
  clinks<-clinks[1:length(unique(pred$to))]                                                #take the number required
  kp <- cbind.data.frame(kp,clinks[order(clinks)])
  names(kp) <- c("pred","pois_prey_num","links")                               #links is the column indicating how many links to keep
  pred<-merge(pred,kp,by.x="to",by.y="pred",all.x=T)                           #add column for number of links to keep for each predator
  out <- split( pred , f = pred$to)                                            #Split data frame
  rand<-lapply(out,rsel)                                                       #randomly select rows for each df, taking into account probability (based on body size)
  rand<-ldply(rand, rbind)                                                     #put dfs in rand together as one df
  
  pls<-oth[oth$from=="plants",]                                                #Add links for herbivores, invertivores and piscivores
  inv<-oth[oth$from=="invertebrates",]
  fishi<-oth[oth$from=="fish",]
  
  plinks<-round(sample(herb_b$breadth, 2*dim(pls)[1], replace=TRUE) + rnorm(2*dim(pls)[1], 0, hd$bw)) #sample from distribution, take extra samples so can remove ones outside the realistic range 
  plinks<-plinks[!(plinks<1)&!(plinks>2*max(herb_b$breadth))]                  #get rid of unrealistic numbers
  plinks<-plinks[1:dim(pls)[1]]                                                #take the number required
  pls$links<-plinks
  pls<-data.frame(pls[rep(seq_len(dim(pls)[1]), pls$links), ], row.names=NULL)
  pls<-pls[,c("to","from")]
  #pls<-pls[pls$to!="invertebrates"&pls$to!="FishInv",]                          #get rid of links from plant nodes to the lumped invertebrate and FishInv nodes
  
  ilinks<-round(sample(inv_b$breadth, 2*dim(inv)[1], replace=TRUE) + rnorm(2*dim(inv)[1], 0, id$bw))
  ilinks<-ilinks[!(ilinks<1)&!(ilinks>2*max(inv_b$breadth))] 
  ilinks<-ilinks[1:dim(inv)[1]] 
  inv$links<-ilinks
  inv<-data.frame(inv[rep(seq_len(dim(inv)[1]), inv$links), ], row.names=NULL)
  inv<-inv[,c("to","from")]
  
  flinks<-round(sample(fish_breadth, 2*dim(fishi)[1], replace=TRUE) + rnorm(2*dim(fishi)[1], 0, fd$bw))
  flinks<-flinks[!(flinks<1)&!(flinks>2*max(fish_breadth))] 
  flinks<-flinks[1:dim(fishi)[1]] 
  fishi$links<-flinks
  fishi<-data.frame(fishi[rep(seq_len(dim(fishi)[1]), fishi$links), ], row.names=NULL)
  fishi<-fishi[,c("to","from")]
  #fishi<-fishi[fishi$to!="invertebrates"&fishi$to!="FishInv",]                     #get rid of links from plant nodes to the lumped invertebrate and FishInv nodes
  
  
  four<-keep80$SciName[keep80$eat.fish=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="TRUE"]           #omnivores, adjust number of links because feed from more than one category
  three<-unique(keep80$SciName[keep80$eat.fish=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="FALSE"|
                              keep80$eat.fish=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.terr_vert=="TRUE"|
                              keep80$eat.fish=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="TRUE"|
                              keep80$eat.fish=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="TRUE"])
  two<-unique(keep80$SciName[  keep80$eat.fish=="TRUE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.terr_vert=="FALSE"|
                                keep80$eat.fish=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="FALSE"&keep80$eat.terr_vert=="TRUE"|
                                keep80$eat.fish=="FALSE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="TRUE"|
                                keep80$eat.fish=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="FALSE"|
                                keep80$eat.fish=="FALSE"&keep80$eat.invert=="TRUE"&keep80$eat.plant=="FALSE"&keep80$eat.terr_vert=="TRUE"|
                                keep80$eat.fish=="TRUE"&keep80$eat.invert=="FALSE"&keep80$eat.plant=="TRUE"&keep80$eat.terr_vert=="FALSE"])
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
  nodes<-data.frame(to=c("invertebrates","invertebrates","invertebrates","fish","fish","fish"),from=c("plants","invertebrates","fish","plants","invertebrates","fish")) #add links between lumped nodes at the end
  rand<-do.call("rbind",list(one,two,three,four,nodes))
  
  return(data.frame(rand))
}

#Now apply it to assemblage 1000 network dataframes, after testing how long it takes to build 1 network
set.seed(10)
start_time <- Sys.time()
x_time <- t(replicate(1, lkp.ap(num80,ed80),simplify = FALSE))
end_time <- Sys.time()
end_time - start_time
set.seed(11)
x80 <- t(replicate(1000, lkp.ap(num80,ed80),simplify = FALSE))
saveRDS(x80,(paste(pp,"/data/networks/x80.rds",sep="")))
#saveRDS(x80,paste("/###/x80nets.rds",sep=""))


