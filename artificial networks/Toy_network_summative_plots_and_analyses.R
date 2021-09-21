#Toy network summative plots and analyses

library(tidyverse)

setwd("###")
#Get the results
df <- list.files(path = "###", pattern = "orig_vuln.rds") %>% #adjust patterns as required
  map_dfr(readRDS, .id = 'ID') %>% 
  bind_rows()


df$ID <- as.numeric(df$ID)

#remova basal nodes from data set
df <- df[df$basal_cons!=0,]

#get some node metrics variation
min(df$tl)
max(df$tl)
min(df$diet_breadth)
max(df$diet_breadth)
min(df$basal_cons)
max(df$basal_cons)

#remove basal_cons = 7 because only one record
df <- df[df$basal_cons!=7,]

#plot both measures for TL, breadth and basal connections, on one pdf
#pdf(paste("###",tre,".pdf",sep=""))
#par(mfrow=c(3,2))
boxplot(df$bay_vuln~df$tl, ylab="bayes vuln", xlab = "trophic level")
boxplot(df$sim_vuln~df$tl, ylab = "sim vuln", xlab = "trophic level")
boxplot(df$bay_vuln~df$diet_breadth, ylab = "bayes vuln", xlab = "diet breadth")
boxplot(df$sim_vuln~df$diet_breadth, ylab = "sim vuln", xlab = "diet breadth")
boxplot(df$bay_vuln~df$basal_cons, ylab = "bayes vuln", xlab = "basal connections")
boxplot(df$sim_vuln~df$basal_cons, ylab = "sim vuln", xlab = "basal connections")
#dev.off()


#fit the models to Eklöf's analytical and simulation results
library(lme4)
library(MuMIn)
library(MASS)
library(car)


#scale, but not centre, data
df$tl <- scale(df$tl, center = FALSE, scale = TRUE)
df$basal_cons <- scale(df$basal_cons, center = FALSE, scale = TRUE)
df$diet_breadth <- scale(df$diet_breadth, center = FALSE, scale = TRUE)



#model sets
ek <- lmer(bay_vuln ~ tl*diet_breadth*basal_cons +(1|ID), data = df,na.action = "na.fail")

sim <- lmer(sim_vuln ~ tl*diet_breadth*basal_cons + (1|ID), data = df,na.action = "na.fail")

#QQplot of residuals - they look okay
qqnorm(resid(sim))
qqnorm(resid(ek))
#other assumptions
plot(resid(ek)~df$tl)
plot(resid(ek)~df$diet_breadth)
plot(resid(ek)~df$basal_cons) #basal_cons = 7 is problematic, only one record and it appears to be an outlier - remove it.
plot(resid(sim)~df$tl)
plot(resid(sim)~df$diet_breadth)
plot(resid(sim)~df$basal_cons)  #Basal_cons = 7 a problem here too

#dredge
ekb <- lmer(bay_vuln ~ tl*diet_breadth*basal_cons +(1|ID), data = df,na.action = "na.fail",REML=F)
ekcomb <- data.frame(dredge(ekb,extra = "R^2"))
simb <- lmer(sim_vuln ~ tl*diet_breadth*basal_cons + (1|ID), data = df,na.action = "na.fail",REML=F)
simcomb <- data.frame(dredge(simb))

#add marginal (variance explained by fixed effects) and conditional (variance explained by random effects) R2 to tables. 
library(performance)
r2_nakagawa(simb)
cond1 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + tl*basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg1 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + tl*basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond2 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg2 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond3 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth +  tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg3 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth +  tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond4 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg4 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond5 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg5 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond6 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg6 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond7 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg7 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[2]))
cond8 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg8 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[2]))
cond9 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg9 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond10 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg10 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond11 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg11 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2]))
cond12 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg12 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond13 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg13 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond14 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg14 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond15 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg15 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond16 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg16 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond17 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg17 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond18 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[1]))
marg18 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond19 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ 1 + (1|ID), data = df,na.action = "na.fail"))[1]))
marg19 <- unname(unlist(r2_nakagawa(lmer(sim_vuln ~ 1 + (1|ID), data = df,na.action = "na.fail"))[2])) 
#add it to results
simcomb$condR2 <- c(cond1,cond2,cond3,cond4,cond5,cond6,cond7,cond8,cond9,cond10,cond11,cond12,cond13,cond14,cond15,cond16,cond17,cond18,cond19)
simcomb$margR2 <- c(marg1,marg2,marg3,marg4,marg5,marg6,marg7,marg8,marg9,marg10,marg11,marg12,marg13,marg14,marg15,marg16,marg17,marg18,marg19)

#now do the same for ekcomb
cond1 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + tl*basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg1 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + tl*basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond2 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg2 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond3 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg3 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond4 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg4 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond5 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg5 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond6 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg6 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + basal_cons*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond7 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg7 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond8 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg8 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond9 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg9 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond10 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg10 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + tl*basal_cons + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond11 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg11 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond12 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg12 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond13 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg13 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond14 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + (1|ID), data = df,na.action = "na.fail"))[1]))
marg14 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond15 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg15 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + diet_breadth + tl*diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond16 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg16 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond17 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[1]))
marg17 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond18 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[1]))
marg18 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[2])) 
cond19 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ 1 + (1|ID), data = df,na.action = "na.fail"))[1]))
marg19 <- unname(unlist(r2_nakagawa(lmer(bay_vuln ~ 1 + (1|ID), data = df,na.action = "na.fail"))[2])) 
#add to results
ekcomb$condR2 <- c(cond1,cond2,cond3,cond4,cond5,cond6,cond7,cond8,cond9,cond10,cond11,cond12,cond13,cond14,cond15,cond16,cond17,cond18,cond19)
ekcomb$margR2 <- c(marg1,marg2,marg3,marg4,marg5,marg6,marg7,marg8,marg9,marg10,marg11,marg12,marg13,marg14,marg15,marg16,marg17,marg18,marg19)


#save tables
#write.csv(ekcomb,"###.csv")
#write.csv(simcomb,"###.csv")

#calculate model-averaged coefficients
#mAVco <- function(x,var){
#  wc <- sum(x[,var]*x$weight,na.rm=T)/sum(x$weight[!(is.na(x[,var]))])
#  return(wc)
#} #Or just use model.avg
simav <- model.avg(dredge(simb))$coefficients[1,]
ekav <- model.avg(dredge(ekb))$coefficients[1,]
#and 95%CIs
CIs <- as.data.frame(confint(model.avg(dredge(simb)), full=T))
CIe <- as.data.frame(confint(model.avg(dredge(ekb)), full=T))
#s95 <- model.avg(dredge(simb), cumsum(weight) <= .95)$coefficients[1,] #very small, not worth plotting
#e95 <- model.avg(dredge(ekb), cumsum(weight) <= .95)$coefficients[1,] #no interval because only one models has weight (much batter than the rest)

#data frame of model averages for plotting
sc <- simav[c("basal_cons","diet_breadth","tl")]
ec <- ekav[c("basal_cons","diet_breadth","tl")]
CIe <- CIe[c("basal_cons","diet_breadth","tl"),]
CIs <- CIs[c("basal_cons","diet_breadth","tl"),]
CIe$method <- "Bayesian"
CIs$method <- "simulation"
sc <- cbind(sc,CIs)
ec <- cbind(ec,CIe)
names(sc)[names(sc)%in%"sc"]<- "coe"
names(ec)[names(ec)%in%"ec"]<- "coe"
dt <- rbind(sc,ec)
dt$var <- gsub('[[:digit:]]+', '', rownames(dt))
names(dt)[2:3] <- c("low","up")
#plot it
  #pdf("###.pdf")
  coefp <- ggplot(data = dt, aes(x = var, y = coe, fill = factor(method) )) + theme_classic() +
  geom_bar(stat = "identity", position=position_dodge()) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
  axis.title.x = element_text(vjust=-0.5),legend.text=element_text(size=12), 
  legend.justification=c(1,0), legend.position=c(0.9,0.3),
  legend.box.background = element_rect(colour = "black"),panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),panel.border = element_blank(),
  panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept=0, color = "black", lwd=1) + labs (x= "variable",y= "coefficient") +
  scale_x_discrete(labels = c("basal connections","diet breadth", "trophic level")) +
  guides(fill=guide_legend(title="method",title.theme = element_text(size = 14, face = "bold"))) + 
  geom_errorbar(aes(ymin=low, ymax=up), colour="black", width=0.2, lwd=1,position=position_dodge(width = 0.9)) +
  #annotate(geom="text", x=2.6,y=-0.22, label="weighted means\nand 95CIs") + 
  scale_fill_manual(values=c("#999999", "#E69F00"))
  #scale_y_continuous(limits = c(-0.3, 0.3))
  #dev.off()

#R2 for main effects, and plot
emarg <- c(unname(unlist(r2_nakagawa(lmer(bay_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])),unname(unlist(r2_nakagawa(lmer(bay_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])),unname(unlist(r2_nakagawa(lmer(bay_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[2]))) 
smarg <- c(unname(unlist(r2_nakagawa(lmer(sim_vuln ~ basal_cons + (1|ID), data = df,na.action = "na.fail"))[2])),unname(unlist(r2_nakagawa(lmer(sim_vuln ~ diet_breadth + (1|ID), data = df,na.action = "na.fail"))[2])),unname(unlist(r2_nakagawa(lmer(sim_vuln ~ tl + (1|ID), data = df,na.action = "na.fail"))[2]))) 
emarg <- data.frame(var = c("basal_cons","diet_breadth","tl"),r2 = emarg, method = "Bayesian")
smarg <- data.frame(var = c("basal_cons","diet_breadth","tl"),r2 = smarg, method = "simulation")
marg <- rbind(emarg,smarg)
#now plot it
  #pdf("###.pdf")
  r2 <- 
  ggplot(data = marg, aes(x = var, y = r2, fill = factor(method) )) + theme_classic() + 
  geom_bar(stat = "identity", position=position_dodge()) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.title.x = element_text(vjust=-0.5),
  legend.text=element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank()) + 
  labs (x= "variable",y= bquote(bold("marginal r"^2))) +
  scale_x_discrete(labels = c("basal connections","diet breadth", "trophic level")) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  guides(fill=guide_legend(title="method",title.theme = element_text(size = 14, face = "bold"))) +
  scale_fill_manual(values=c("#999999", "#E69F00")) #+
    #geom_errorbar(aes(ymin=low, ymax=up), colour="black", width=0.2, lwd=1,position=position_dodge(width = 0.9))
  #dev.off()

  #stick coef and r2 plots together
  library(ggpubr)
  #pdf("###.pdf",width = 10, height = 5)
  ggarrange(coefp,r2, 
            labels = c("A", "B"),
            ncol = 2, nrow = 1)
  #dev.off()
  
#network metrics
nm <- list.files(path = "~/Dropbox/trophic links/Docs/Full ms/Ecography/revision 1/script/Giovanni/for_john/toy_network_mets", pattern = ".rds") %>%
  map_dfr(readRDS, .id = 'ID') %>% 
  bind_rows()

#N = number of nodes, Ltot = number of links, LD = link density, C = connectance
min(nm$N) # 3
max(nm$N) # 20
min(nm$Ltot) #2
max(nm$Ltot) #92
min(nm$C) #0.07
max(nm$C) #0.33


