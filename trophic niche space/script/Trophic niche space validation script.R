#Get predator/prey interaction data
#Determine body size relationship, and whether this relationship is affected by predator class
library(quantreg)

#project path
pp <- getwd()

#get and tidy data 
#Extracted from GloBI -- April 2019
data<-read.csv(paste((pp),"/data/interactions_globi_cat.csv",sep="")) #Mass in grams
#merge reptiles and amphibians in class category
data$class_resource<-ifelse(data$class_resource=="amphibians","amphibia and reptilia",as.character(data$class_resource))
data$class_resource<-ifelse(data$class_resource=="reptilia","amphibia and reptilia",as.character(data$class_resource))
data$class_consumer<-ifelse(data$class_consumer=="amphibians","amphibia and reptilia",as.character(data$class_consumer))
data$class_consumer<-ifelse(data$class_consumer=="reptilia","amphibia and reptilia",as.character(data$class_consumer))

#remove marine mammals
data<-subset(data,data$category_consumer!="marmam")

#Function to fit different quantile regressions# Bprey = prey mass, Bpred = predator mass, Cpred = predator class, ta = quantile
QR <- function(Bprey,Bpred,Cpred,ta) {           
  ft1 = rq(Bprey~1,tau=ta)    #Null model
  ft2 = rq(Bprey~Cpred,tau=ta)  	
  ft3 = rq(Bprey~Bpred,tau=ta)   	
  ft4 = rq(Bprey~Bpred+Cpred,tau=ta)  
  ft5 = rq(Bprey~Bpred*Cpred,tau=ta)  
  
  return(list(AIC(ft1,k=log(length(Bprey))),AIC(ft2,k=log(length(Bprey))),   #BIC = AIC with penalty k=log(n) #n = number of observations
              AIC(ft3, k=log(length(Bprey))),AIC(ft4, k=log(length(Bprey))),
              AIC(ft5, k=log(length(Bprey)))#
              ))
}

#assign variables
Bprey<-jitter(log10(data$body_mass_resource)) #jitter for quntile regression, so not lots of repeats of exact same x or y value# mass in grams
Bpred<-jitter(log10(data$body_mass_consumer))
Cprey<-data$class_resource
Cpred<-data$class_consumer

#apply function
qr.fts<-QR(Bprey,Bpred,Cpred,0.95) 
qr.fts2<-QR(Bprey,Bpred,Cpred,0.05)

qr.fts #Model ft5 
qr.fts2#Model ft4 

#View relationship between lm and quantiles for each class
brds<-which(Cpred%in%"aves")
mams<-which(Cpred%in%"mammals")
herps<-which(Cpred%in%"amphibia and reptilia")

Cpred[Cpred=="aves"]<-"birds"
Cpred[Cpred=="amphibia and reptilia"]<-"amphibians and reptiles"

class.plot<-function(x){ #x=class
plot(Bpred[x],Bprey[x],pch=20, xlab = "log10 predator mass (g)", ylab = "log10 prey mass (g)",main = unique(Cpred[x]), font.lab=2)
abline(lm(Bprey[x]~Bpred[x]),col="green",lwd=3)
abline(rq(Bprey[x]~Bpred[x],tau=0.05),col="red",lwd=3 )
abline(rq(Bprey[x]~Bpred[x],tau=0.95),col="red",lwd=3 )
}

#plot it
par(mfrow=c(2,2))
par(mgp = c(2, 0.5, 0))
class.plot(brds)
class.plot(mams)
class.plot(herps)

#Model validation
#Some functions for validation########################################
#Function to prepare training data from GloBI
fix.dat<- function(x){
  x$class_consumer<-as.character(x$class_consumer)
  x$class_consumer[x$class_consumer=="amphibians"]<-"amphibia and reptilia"
  x$class_consumer[x$class_consumer=="reptilia"]<-"amphibia and reptilia"
  x$class_resource<-as.character(x$class_resource)
  x$class_resource[x$class_resource=="amphibians"]<-"amphibia and reptilia"
  x$class_resource[x$class_resource=="reptilia"]<-"amphibia and reptilia"
  Bprey<-as.numeric(jitter(log10(x$body_mass_resource))) 
  Bpred<-as.numeric(jitter(log10(x$body_mass_consumer)))
  Cprey<-as.character(x$class_resource)
  Cpred<-as.character(x$class_consumer)
  x<-cbind.data.frame(Bprey,Bpred,Cprey,Cpred)
}

#Function to prepare validation data
fix.ndat<- function(x){
  x$class_consumer<-as.character(x$class_consumer)
  x$class_consumer[x$class_consumer=="amphibians"]<-"amphibia and reptilia"
  x$class_consumer[x$class_consumer=="reptilia"]<-"amphibia and reptilia"
  x$class_resource<-as.character(x$class_resource)
  x$class_resource[x$class_resource=="amphibians"]<-"amphibia and reptilia"
  x$class_resource[x$class_resource=="reptilia"]<-"amphibia and reptilia"
  Bprey<-as.numeric(jitter(log10(x$body_mass_resource))) 
  Bpred<-as.numeric(jitter(log10(x$body_mass_consumer)))
  Cprey<-as.character(x$class_resource)
  Cpred<-as.character(x$class_consumer)
  obs<-as.character(x$occ)
  x<-cbind.data.frame(Bprey,Bpred,Cprey,Cpred,obs)
}

#Functions for getting predicted and observed presence/absence of interactions
#simple body mass model
Mpred<-function(data,newd)
{
  #fix training data
  data<-fix.dat(data)
  #quantile regressions
  qrup<-rq(Bprey~Bpred,data=data,tau=0.95)
  qrlow<-rq(Bprey~Bpred,data=data,tau=0.05)
  #fix validation data
  newd<-fix.ndat(newd)
  #make predictions using regressions from training data
  newd$low<-predict.rq(qrlow, newdata= newd)
  newd$up<-predict.rq(qrup, newdata= newd)
  newd$predicted<-newd$Bprey<=newd$up&newd$Bprey>=newd$low
  newd$predicted[newd$predicted=="TRUE"]<-1
  print(table(paste(newd$obs,newd$predicted)))
}


#More complex model, diff models for upper vs lower quantiles (mass*class for upper, mass + class for lower)
Mpred2<-function(data,newd)
{
  #fix training data
  data<-fix.dat(data)
  #quantile regressions
  qrup<-rq(Bprey~Bpred*Cpred,data=data,tau=0.95, method='fn')
  qrlow<-rq(Bprey~Bpred+Cpred,data=data,tau=0.05, method='fn')
  #fix validation data
  newd<-fix.ndat(newd)
  #make predictions using regressions from training data
  newd$low<-predict.rq(qrlow, newdata= newd)
  newd$up<-predict.rq(qrup, newdata= newd)
  newd$predicted<-newd$Bprey<=newd$up&newd$Bprey>=newd$low
  newd$predicted[newd$predicted=="TRUE"]<-1
  print(table(paste(newd$obs,newd$predicted)))
}


#true skill statistic functions
#for multiple sets of results
tssF<-function(x) {
  ((x[,4]*x[,1])-((x[,2]*x[,3])))/((x[,4]+x[,3])*(x[,2]+x[,1])) 
}

#for one set of results
tssF1<-function(x) {
  ((x[4]*x[1])-((x[2]*x[3])))/((x[4]+x[3])*(x[2]+x[1])) 
}

#######################################################################
#get training and validation data (GloBi data split, with unobserved interactions generated)
setwd("data/training/")
temp = list.files(pattern="*.csv")
myTfiles = lapply(temp, read.csv)
setwd(paste(pp,"/data/validation/",sep=""))
temp = list.files(pattern="*.csv")
myVfiles = lapply(temp, read.csv)

#Simple body size model 
#Function for calculating predicted/observed etc
#First, make empty data frame
cats1 <- data.frame('00' = numeric(length(temp)),
                    '01' = numeric(length(temp)),
                    '10' = numeric(length(temp)),
                    '11' = numeric(length(temp)))
#Make predictions
run.all<-for (i in 1:length(temp)) {
  tdat<-myTfiles[[i]]
  vdat<-myVfiles[[i]]
  cats1[i,]<-Mpred(tdat,vdat)
}


#calculate the True Skill Statistic
colnames(cats1)<-c("d","b","c","a")
cats1$tss<-(cats1$a*cats1$d-cats1$b*cats1$c)/((cats1$a + cats1$c)*(cats1$b+cats1$d))

#Compare to model with diff regression for upper vs lower (best regressions according to BIC)
cats2 <- data.frame('00' = numeric(length(temp)),
                    '01' = numeric(length(temp)),
                    '10' = numeric(length(temp)),
                    '11' = numeric(length(temp)))

run.all<-for (i in 1:length(temp)) {
  tdat<-myTfiles[[i]]
  vdat<-myVfiles[[i]]
  cats2[i,]<-Mpred2(tdat,vdat)
}

#calculate the True Skill Statistic
colnames(cats2)<-c("d","b","c","a")
cats2$tss<-(cats2$a*cats2$d-cats2$b*cats2$c)/((cats2$a + cats2$c)*(cats2$b+cats2$d))

#Obs:Pred compare
colMeans(cats1) #simple body size model
colMeans(cats2) #predator class*body size interaction for upper regression, addititive for lower
#The second model performs best

#########################################################################
#Have a look at the model 
#assign variables
Bprey<-jitter(log10(data$body_mass_resource)) 
Bpred<-jitter(log10(data$body_mass_consumer))
Cprey<-data$class_resource
Cpred<-data$class_consumer

#Best model, different regs for upper vs lower quantiles 
qrup2<-rq(Bprey~Bpred*Cpred,tau=0.95)  
qrlow2<-rq(Bprey~Bpred+Cpred,tau=0.05) 
#simple body size model
qrup<-rq(Bprey~Bpred,tau=0.95)  
qrlow<-rq(Bprey~Bpred,tau=0.05) 

###Where are the quantiles for these records?
#plot predictions
#pdf("/Class panels and best model slopes.pdf")
par(mfrow=c(2,2))
par(mgp = c(2, 0.5, 0))
class.plot(brds)
class.plot(mams)
class.plot(herps)
plot(x=NULL,y=NULL,xlim=c(0.3,6),ylim=c(0,8),xlab="log10 predator mass (g)",ylab="log10 prey mass (g)", font.lab=2, main = "best model")
abline(a = qrlow2$coefficients[1], b = qrlow2$coefficients[2],lty = 1, lwd = 1, col="green") #reptiles and amphibians
abline(a = qrlow2$coefficients[1]+qrlow2$coefficients[3],b = qrlow2$coefficients[2],lty = 1, lwd = 1,col="red") #birds
abline(a = qrlow2$coefficients[1]+qrlow2$coefficients[4],b = qrlow2$coefficients[2],lty = 1, lwd = 1,col="brown") #mammals
abline(qrup2$coefficients[1:2], lty = 1, lwd = 1, col="green")
abline(qrup2$coefficients[1]+qrup2$coefficients[3],qrup2$coefficients[2]+qrup2$coefficients[5],lty = 1, lwd = 1,col="red")
abline(qrup2$coefficients[1]+qrup2$coefficients[4],qrup2$coefficients[2]+qrup2$coefficients[6],lty = 1, lwd = 1,col="brown")
legend(0.25, 8, legend=c("amphibians and reptiles", "birds","mammals"),
       col=c("green", "red", "brown"), lty=1, cex=0.8, box.lty=0, bg = "transparent")
#dev.off()

###############################################################################################
#Test on Serengeti data - detailed terrestrial food web including large vertebrates (>20kg), as well as small vertebrates.
#Need all the potential interactions and all the observed interactions
#Partly from: The Serengeti food web: empirical quantification and analysis of topological changes under increasing human impact.
#de Visser SN1, Freymann BP, Olff H. 2011 Journal of Animal Ecology 80(2)
#Also de Visser unpublished + relevant GloBI records
ser<-read.csv(paste(pp,"/data/Predator-Prey Data.csv",sep=""))
ser<-subset(ser,ser$source=="tog")

#Fix black mumba mass
ser$PredMassKg[ser$Predator=="Dendroaspis polylepis"]<-1.6
ser$PreyMassKg[ser$Prey=="Dendroaspis polylepis"]<-1.6

#change mass to log10(g)
ser$LogPred<-log10(ser$PredMassKg*1000)
ser$LogPrey<-log10(ser$PreyMassKg*1000)

#add column for prey diet
ser$prey.diet<-ifelse(ser$Prey%in%ser$Predator,"Carni","?")
preyCh<-unique(ser$Prey[ser$prey.diet=="?"])

#remove strictly aquatic species 
ser<-subset(ser,ser$Prey!="Perciformes")
ser<-subset(ser,ser$Prey!="Xenopus muelleri")

#add general diet for prey from literature #But could also do this based on Elton traits database 
ser$prey.diet[ser$Prey=="Phrynobatrachus mababiensis"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Crocidura sp."]<-"Omnivore"
ser$prey.diet[ser$Prey=="Pipistrellus nanus"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Empidonax wrightii"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Apalis flavida"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Eremomela icteropygialis"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Sylvietta whytii"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Charadrius tricollaris"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Micropteropus pusillus"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Steatomys pratensis"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Anas capensis"]<-"Herbivore/Invertivore"
ser$prey.diet[ser$Prey=="Lygodactylus capensis"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Pelomys fallax"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Pterocles gutturalis"]<-"Herbivore/Invertivore"
ser$prey.diet[ser$Prey=="Francolinus coqui"]<-"Herbivore/Invertivore"
ser$prey.diet[ser$Prey=="Heterohyrax brucei"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Pronolagus rupestris"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Lutra maculicollis"]<-"Carnivore"
ser$prey.diet[ser$Prey=="Colobus guereza"]<-"Herbivore/Invertivore"
ser$prey.diet[ser$Prey=="Manis temminckii"]<-"Invertivore"
ser$prey.diet[ser$Prey=="Madoqua kirkii"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Hystrix cristata"]<-"Omnivore"
ser$prey.diet[ser$Prey=="Ourebia ourebi"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Redunca redunca"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Damaliscus lunatus"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Syncerus caffer"]<-"Herbivore"
ser$prey.diet[ser$Prey=="Tragelaphus oryx"]<-"Herbivore"

ser$prey.diet[ser$prey.diet=="Carni"|ser$prey.diet=="Carnivore"|ser$prey.diet=="Omni"|ser$prey.diet=="Omnivore"]<-"vert.pred"
ser$prey.diet[ser$prey.diet!="vert.pred"]<-"not.pred"
ser$pred.diet<-"vert.pred" #Anything recorded to have eaten a vertebrate should be a vertebrate predator (though scavangers may be included)

#fix class data
ser$PredClass<-as.character(ser$PredClass)
ser$PredClass[ser$PredClass=="Reptilia"]<-"amphibia and reptilia"
ser$PredClass[ser$PredClass=="Aves"]<-"aves"
ser$PredClass[ser$PredClass=="Mammalia"]<-"mammals"

ser$PreyClass<-as.character(ser$PreyClass)
ser$PreyClass[ser$PreyClass=="Reptilia"|ser$PreyClass=="Amphibia"]<-"amphibia and reptilia"
ser$PreyClass[ser$PreyClass=="Aves"]<-"aves"
ser$PreyClass[ser$PreyClass=="Mammalia"]<-"mammals"

#make long list of all species, diet, and body mass
prey<-unique(ser[,c("Prey","LogPrey","prey.diet","PreyClass")])
pred<-unique(ser[,c("Predator","LogPred","pred.diet","PredClass")])
cn<-c("species","logmass","diet","Cpred")
colnames(prey)<-cn
colnames(pred)<-cn

ll<-rbind.data.frame(prey,pred)
ll<- ll[!duplicated(ll$species),] 
ll$numb<-1:dim(ll)[1]
species<-ll[,c("species","numb")]
species$species<-as.character(species$species)

library(gtools)
options(expressions=100000)
pll<-data.frame(permutations(n=length(species$species),r=2,v=species$species,repeats=TRUE))
colnames(pll)<-c("cons","res")

#merge body size and class data
cons<-ll[,c("species","logmass","Cpred")]
colnames(cons)<-c("cons","cons.logmass","Cpred")
pll<-merge(pll,cons,by="cons",all.x=T)
res<-ll[,c("species","logmass")]
colnames(res)<-c("res","res.logmass")
pll<-merge(pll,res,by="res",all.x=T)

#add col for observed interactions
pll$obs<-(paste(pll$res,pll$cons)%in%paste(ser$Prey,ser$Predator))
pll$obs<-ifelse(pll$obs=="TRUE",1,0)

#Check in rGlobi - is the data set missing trophic links?
#library(rglobi)
#library(devtools)
#library(igraph)
#dd<-list()
#sp<-unique(ll$species)
#for (i in 1:length(sp)){
#dd[[i]]<-as.data.frame(get_interactions(sp[i], interaction.type = "preysOn"))
#}
#now stick the objects together
#dd<-do.call("rbind", dd)
#write.csv(dd,paste(pp, "/data/ll.prey.csv", sep="))
ll.prey<-read.csv(paste(pp, "/data/ll.prey.csv", sep=""))

#restrict to only those with predator and prey in the Serengeti data set, and add 1s for links on GloBi that weren't in the data set yet
sll<-subset(ll.prey,ll.prey$target_taxon_name%in%ll$species)
sll<-unique(sll[,c("source_taxon_name","target_taxon_name")])
table(paste(sll$source_taxon_name,sll$target_taxon_name)%in%paste(pll$cons[pll$obs=="1"],pll$res[pll$obs=="1"])) #4 missing links
sll<-sll[paste(sll$source_taxon_name,sll$target_taxon_name)%in%paste(pll$cons[pll$obs=="1"],pll$res[pll$obs=="1"])=="FALSE",] 
pll$obs<-ifelse(paste(pll$cons,pll$res)%in%paste(sll$source_taxon_name,sll$target_taxon_name),"1",pll$obs)


#Make Serengeti predictions based on GloBI data, only using body size (no diet filter)
#Simple body size model 
newd<-pll
colnames(newd)<-c("prey","pred","Bpred","Cpred","Bprey","obs")
newd$low<-predict.rq(qrlow, newdata= newd)
newd$up<-predict.rq(qrup, newdata= newd)
newd$predicted<-newd$Bprey<=newd$up&newd$Bprey>=newd$low
newd$predicted<-ifelse(newd$predicted=="TRUE",1,0)
newd$outcome<-paste(newd$obs,newd$predicted)
#True skill statistic
tssF1(table(newd$outcome)) #0.36

#Compare to best model (upper qr = mass*class, lower qr = mass + class)
newd<-pll
colnames(newd)<-c("prey","pred","Bpred","Cpred","Bprey","obs")
newd$low<-predict.rq(qrlow2, newdata= newd)
newd$up<-predict.rq(qrup2, newdata= newd)
newd$predicted<-newd$Bprey<=newd$up&newd$Bprey>=newd$low
newd$predicted<-ifelse(newd$predicted=="TRUE",1,0)
newd$outcome<-paste(newd$obs,newd$predicted)
#True skill statistic
tssF1(table(newd$outcome)) #This model is still better

#add information on diet (whether a species preys on vertebrates), and recalculate TSS
colnames(newd)[colnames(newd)=="prey"]<-"target"
colnames(newd)[colnames(newd)=="pred"]<-"consumer"
colnames(ll)[colnames(ll)=="species"]<-"consumer"
newd<-merge(newd,ll[,c("consumer","diet")],all.x=TRUE)
#fix predicted and recalculate outcome
newd$predicted<-ifelse(newd$diet=="not.pred",0,newd$predicted)
newd$outcome<-paste(newd$obs,newd$predicted)
#True skill statistic, taking diet into consideration
tssF1(table(newd$outcome))  #0.6 with diff reg model, 0.58 for simple mass model

#plot model
pred1<-subset(newd,newd$diet=="vert.pred")
color_transparent <- adjustcolor(c("red","blue","black","green"), alpha.f = 0.3) 
pred1$outcome[pred1$outcome=="1 1"] <- color_transparent[1]
pred1$outcome[pred1$outcome=="0 1"] <- color_transparent[2]
pred1$outcome[pred1$outcome=="0 0"] <- color_transparent[3]
pred1$outcome[pred1$outcome=="1 0"] <- color_transparent[4]

par(mfrow = c(1, 1))
par(mgp = c(2, 0.5, 0))
par(xpd=TRUE)
par(mar=c(4,4,4,6))
plot(pred1$Bpred,pred1$Bprey,col=pred1$outcome,pch=19,xlab="log10 predator mass (g)",
     ylab="log10 prey mass (g)",main="Serengeti food web", font.lab=2, bty='L')
legend(6.05, 6, legend=c("absent  absent", "present present", "absent  present","present absent"),
       col=c(color_transparent[3], color_transparent[1], color_transparent[4], color_transparent[2]),
       pch=19, cex=0.8, box.lty=1, bg = "transparent", title = expression(bold("predicted:observed")))


###############################################################################################################
#Spatial Guilds in the Serengeti Food Web Revealed by a Bayesian Group Model
#Baskerville, Edward B., Andy P. Dobson, Trevor Bedford, Stefano Allesina, T. Michael Anderson, and Mercedes Pascual. “Spatial Guilds in the Serengeti Food Web Revealed by a Bayesian Group Model.” Edited by Lauren Ancel Meyers. PLoS Computational Biology 7, no. 12 (December 29, 2011): e1002321. https://doi.org/10.1371/journal.pcbi.1002321.
#Get the interaction, species code, and animal list data (to differentiate animals and plants)
net<-read.csv(paste(pp,"/data/Baskerville et al Serengeti.csv",sep=""),skip=1)
sp<-read.csv(paste(pp,"/data/Baskerville et al species.csv",sep=""),skip=1)
animals<-read.csv(paste(pp,"/data/Baskerville animals.csv",sep=""),header = FALSE,colClasses = "character")
animals$V1<-trimws(animals$V1,which="both")
net<-merge(net,sp,by.x = 1, by.y = 1,all.x=TRUE)
net<-merge(net,sp,by.x = 2, by.y = 1,all.x=TRUE)
net<-net[,3:4]
colnames(net)<-c("pred","prey")
net<-apply(net,2,as.character)
net<-as.data.frame(net)
#only keep predator/prey interactions (get rid of plants)
net<-subset(net,net$prey%in%animals$V1)

#Now make a data frame with species, mass and diet data.
colnames(animals)<-"species"
#Get EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals data
#Wilman, Hamish, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, Marcelo M. Rivadeneira, and Walter Jetz. 
#“EltonTraits 1.0: Species-Level Foraging Attributes of the World’s Birds and Mammals: Ecological Archives E095-178.” Ecology 95, no. 7 (July 2014): 2027–2027. https://doi.org/10.1890/13-1917.1.
et<-read.csv(paste(pp,"/data/MamFuncDat.csv",sep=""))
et$vert.pred<-et$Diet.Vect>0|et$Diet.Vend>0|et$Diet.Vunk
et<-et[,c("Scientific","vert.pred","BodyMass.Value")]
animals<-merge(animals, et, by.x=1, by.y=1, all.x=TRUE)

#get all the species permutations
perms<-as.data.frame(permutations(n=32,r=2,v=animals$species,repeats.allowed=T))
colnames(perms)<-c("source","target")
#merge body size and diet data
perms<-merge(perms, animals, by.x=1, by.y=1, all.x=TRUE)
colnames(perms)[colnames(perms)=="vert.pred"]<-"source.vert.pred"
colnames(perms)[colnames(perms)=="BodyMass.Value"]<-"source.mass"
perms<-merge(perms, animals, by.x=2, by.y=1, all.x=TRUE)
perms$vert.pred<-NULL
colnames(perms)[colnames(perms)=="BodyMass.Value"]<-"target.mass"

#add column for observed interactions
perms$obs<-(paste(perms$target,perms$source)%in%paste(net$prey,net$pred))
perms$obs<-ifelse(perms$obs=="TRUE",1,0)

#Check if GLOBI adds any links
#bd<-list()
#sp<-unique(perms$source)
#for (i in 1:length(sp)){
#bd[[i]]<-as.data.frame(get_interactions(sp[i], interaction.type = "preysOn"))
#}
#now stick the objects together
#bd<-do.call("rbind", bd)
#write.csv(bd,"GLOBI.BS.prey.csv")
BS.prey<-read.csv(paste(pp,"/data/GLOBI.BS.prey.csv",sep=""))

#restrict to only those with predator and prey in dataset, and add 1s for links in GloBi that weren't already in the network
bll<-subset(BS.prey,BS.prey$target_taxon_name%in%perms$target)
bll<-unique(bll[,c("source_taxon_name","target_taxon_name")])
table(paste(bll$source_taxon_name,bll$target_taxon_name)%in%paste(perms$source[perms$obs=="1"],perms$target[perms$obs=="1"])) #3 missing links
bll<-bll[paste(bll$source_taxon_name,bll$target_taxon_name)%in%paste(perms$source[perms$obs=="1"],perms$target[perms$obs=="1"])=="FALSE",] 
perms$obs<-ifelse(paste(perms$source,perms$target)%in%paste(bll$source_taxon_name,bll$target_taxon_name),"1",perms$obs)

#log transform mass
perms$source.mass<-log10(perms$source.mass)
perms$target.mass<-log10(perms$target.mass)

#Predator class col
perms$Cpred<-"mammals"

#Make predictions#
colnames(perms)<-c("prey","pred","vert.pred","Bpred","Bprey","obs","Cpred")
perms$low<-predict.rq(qrlow2, newdata= perms)
perms$up<-predict.rq(qrup2, newdata= perms)
perms$predicted<-perms$Bprey<=perms$up&perms$Bprey>=perms$low
perms$predicted<-ifelse(perms$predicted=="TRUE",1,0)
perms$outcome<-paste(perms$obs,perms$predicted)

#True skill statistic, not taking diet into consideration
tssF1(table(perms$outcome))  #Class*mass model TSS = 0.18 - does poorly because so many herbivores 

#add information on diet (whether a species preys on vertebrates), and recalculate TSS
#fix predicted and recalculate outcome
perms$predicted<-ifelse(perms$vert.pred=="FALSE",0,perms$predicted)
perms$outcome<-paste(perms$obs,perms$predicted)
#True skill statistic, taking diet into consideration
tssF1(table(perms$outcome)) #Class*mass model TSS = 0.73 - 
#remove rows of herbivores as predators
perms2<-perms[perms$vert.pred=="TRUE",]
tssF1(table(perms2$outcome))

#plot it
#restrict to preds
pred<-subset(perms,perms$vert.pred=="TRUE")
#Make outcomes colours that are partly transparent
pred$outcome[pred$outcome=="1 1"] <- color_transparent[1]
pred$outcome[pred$outcome=="0 1"] <- color_transparent[2]
pred$outcome[pred$outcome=="0 0"] <- color_transparent[3]
pred$outcome[pred$outcome=="1 0"] <- color_transparent[4]

par(mgp = c(2, 0.5, 0))
par(xpd=TRUE)
par(mar=c(4,4,4,8))
plot(pred$Bpred,pred$Bprey,col=pred$outcome,pch=19,xlab="log10 predator mass (g)",ylab="log10 prey mass (g)",
     main="Serengeti spatial guilds")
legend(5.3, 6.8, legend=c("absent  absent", "present present", "absent  present","present absent"),
       col=c(color_transparent[3], color_transparent[1], color_transparent[4], color_transparent[2]),
       pch=19, cex=0.8, box.lty=1, bg = "transparent", title = expression(bold("predicted:observed")))

##################################################
#Test for a relationship between predator body size/number predicted prey and proportion of appropriately sized prey species that a predator takes - data from de Visser et al
#The Serengeti food web: empirical quantification and analysis of topological changes under increasing human impact.
#de Visser SN1, Freymann BP, Olff H. 2011
#Also de Visser unpublished + relevant GloBI records
pred.ob<-data.frame(cbind(table(pred1$predicted,pred1$consumer)[2,],table(pred1$obs,pred1$consumer)[2,]))
colnames(pred.ob)<-c("predicted","observed")
pred.ob$species<-rownames(pred.ob)
pred.ob<-merge(pred.ob,unique(pred1[,c("consumer","Bpred")]), by.x=3, by.y=1, all.x=TRUE)
pred.ob$prop<-pred.ob$observed/(pred.ob$predicted)
mass.prop<-lm(pred.ob$prop~pred.ob$Bpred)
summary(mass.prop) #no sign of a correlation
#get rid of the 0 observed - doesn't add anything
pred.ob<-pred.ob[pred.ob$observed!=0,]
plot(pred.ob$Bpred,pred.ob$prop,xlab="log10 predator body mass (g)",ylab="prey diversity taken ÷ predicted", main = "proportion by predator mass")
abline(mass.prop)
summary(lm(pred.ob$prop~pred.ob$Bpred)) #No indication of a relationship

s.lm <- lm(pred.ob$prop~pred.ob$predicted)
s.lm.n <- lm(pred.ob$prop~1)
BIC(s.lm)-BIC(s.lm.n) #Null model has lower BIC
delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
mods.dAIC <- delta.AIC(c(BIC(s.lm),BIC(s.lm.n)))
weight.AIC(mods.dAIC)

d <- density(pred.ob$prop) # returns the density data
saveRDS(d,paste(pp,"/data/Serengeti.density.rds", sep=""))
saveRDS(pred.ob,paste(pp,"/data/pred.ob.rds", sep=""))

plot(pred.ob$Bpred,pred.ob$prop,xlab="log10 predator body mass (g)",ylab="prey diversity taken ÷ predicted", main = "proportion by predator mass")
plot(pred.ob$predicted,pred.ob$prop,xlab="diversity of prey predicted",ylab="prey diversity taken ÷ predicted", main = "proportion by diversity predicted")

num.ln<-summary(lm(pred.ob$observed~pred.ob$predicted)) #not significant
plot(pred.ob$predicted,pred.ob$observed,xlab="Number of prey predicted",ylab="Number observed")
abline(num.ln$coefficients[1],num.ln$coefficients[2])
mean(pred.ob$prop)

