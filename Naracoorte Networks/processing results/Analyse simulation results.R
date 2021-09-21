#Get simulation results and TL (from Trophic level and results.R) and fit models

#set working dir and get data
setwd("~/###")
sim <- readRDS("/data/simTL.rds")

#express vulnerability as the proportion of plants needed to be removed to trigger coextinction
sim$vulnerability <- sim$vulnerability/sim$plant_n

#names extinction status categories
sim$extinct_status <- ifelse(sim$extinct_status==0, "extant","extinct")

#get some node metrics variation
min(sim$TL)
max(sim$TL)
min(sim$deg_in)
max(sim$deg_in)
min(sim$paths_to_bas)
max(sim$paths_to_bas)

#restrict to runs with low tre (need to lose almost all food)
sim <- sim[sim$tre>0.98,]

##get residual vulnerability
resV <- resid(lm(sim$vulnerability~sim$tre))

#log transform
sim$deg_in <- log(sim$deg_in)
sim$basal_cons <- log(sim$paths_to_bas)
sim$TL <- log(sim$TL)

#scale but don't centre
sim$TL <- scale(sim$TL, center = FALSE, scale = TRUE)
#sim$basal_cons <- scale(sim$basal_cons, center = FALSE, scale = TRUE)
sim$deg_in <- scale(sim$deg_in, center = FALSE, scale = TRUE)
sim$paths_to_bas <- scale(sim$paths_to_bas, center = FALSE, scale = TRUE)

#aggregate for each species (means)
sim$resV <- resV
aggs <- aggregate(.~species, sim[,c("species","paths_to_bas","TL","deg_in","vulnerability","resV")], mean)
names(aggs) <- c("species","paths_to_bas","TL","deg_in","vulnerability","resV")
aggs$extinct_status <- ifelse(aggs$species%in%sim$species[sim$extinct_status=="extant"],"extant","extinct")


#pdf("###.pdf")
#par(mfrow=c(2,2))
plot(aggs$resV~aggs$paths_to_bas, xlab="basal connections",ylab="resid coext vuln")
abline(lm(aggs$resV~aggs$paths_to_bas),col="red")
plot(aggs$resV~aggs$deg_in, xlab="diet breadth",ylab="resid coext vuln")
abline(lm(aggs$resV~aggs$deg_in),col="red")
plot(aggs$resV~aggs$TL, xlab="trophic level",ylab="resid coext vuln")
abline(lm(aggs$resV~aggs$TL),col="red")
#dev.off()

library(lme4)
fts <- lmer(resV ~ paths_to_bas*deg_in*TL + (1|species), data = sim,na.action = "na.fail")
fts2 <- lmer(resV ~ paths_to_bas*deg_in*TL + (1|species) + (1|ID), data = sim,na.action = "na.fail",REML = F)
fts3 <- lmer(resV ~ paths_to_bas*deg_in*TL + (1|ID), data = sim,na.action = "na.fail",REML = F) #species explains none of the variance, and AIC/BIC indicate support models without species, so remove it
#AIC and BIC support fts3 so use that random effect and see if extinct versus extant differ
fts_ex <- lmer(resV ~ extinct_status + (1|ID), data = sim,na.action = "na.fail",REML = F) # t = -2.235 but lots of overlap
#check if residuals different for extinct versus extant
boxplot(resid(fts)~sim$extinct_status)

#weight of evidence for model with versus without extinction status
with <- lmer(resV~extinct_status + (1|ID), data = sim, REML = T, na.action = "na.fail")
without <- lmer(resV~1 + (1|ID), data = sim, REML = T, na.action = "na.fail")
BIC(with,without) 

#violin plot of vulnerability
library(ggplot2)
library(dplyr)
v <- ggplot(sim, aes(x=extinct_status, y=resV, fill=extinct_status)) + geom_violin(show.legend = FALSE) + 
  theme_bw() + labs(x = "extinction status", y = "vulnerability (resid)") #+
  #geom_point(shape = 21,size=0.5, position = position_jitterdodge(0.5), color="black",alpha=0.2)
  
b <- ggplot(sim, aes(x=extinct_status, y=paths_to_bas, fill=extinct_status)) + geom_violin(show.legend = FALSE) + 
  theme_bw() + labs(x = "extinction status", y = "basal connections (scaled)") #+
  #geom_point(shape = 21,size=0.5, position = position_jitterdodge(0.5), color="black",alpha=0.2)

d <- ggplot(sim, aes(x=extinct_status, y=deg_in, fill=extinct_status)) + geom_violin(show.legend = FALSE) + 
  theme_bw() + labs(x = "extinction status", y = "diet breadth (scaled)") #+
  #geom_point(shape = 21,size=0.5, position = position_jitterdodge(0.5), color="black",alpha=0.2)

t <- ggplot(sim, aes(x=extinct_status, y=TL, fill=extinct_status)) + geom_violin(show.legend = FALSE) + 
  theme_bw() + labs(x = "extinction status", y = "trophic level (scaled)") #+
#geom_point(shape = 21,size=0.5, position = position_jitterdodge(0.5), color="black",alpha=0.2)
library(ggpubr)
#pdf("###.pdf",width = 10, height = 10)
ggarrange(v, b, d, t,
labels = c("A", "B","C","D"),
ncol = 2, nrow = 2)
#dev.off()

#density plot of vulnerability
vul <- ggplot(data=aggs, aes(x=resV, group=as.factor(extinct_status), fill=as.factor(extinct_status))) +
  geom_density(adjust=1.5, alpha=.4) + guides(fill=guide_legend(title="extinction status")) +xlim(-0.14,0.15) +
  annotate(geom="text", x=0.075, y=50, label="tre > 0.9",color="red")

#bin trophic levels
breaks <- c(0,0.8,1.2,1.6)
tags <- c("[0-0.8]","[0.8-1.2]","[1.2-1.6]")
tl <- cut(aggs$TL, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)
tlp <- ggplot(data=aggs, aes(x=resV, group=tl, fill=tl)) +
  geom_density(adjust=1.5, alpha=.4) + guides(fill=guide_legend(title="trophic level\n(scaled)"))

#bin diet breadth
breaks <- c(0,0.5,1,1.5)
tags <- c("[0-0.5]","[0.5-1]","[1-1.5]")
db <- cut(aggs$deg_in, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)
dbp <- ggplot(data=aggs, aes(x=resV, group=db, fill=db)) +
  geom_density(adjust=1.5, alpha=.4) + guides(fill=guide_legend(title="diet breadth\n(scaled)"))

#bin basal cons
breaks <- c(0,0.7,1.4,2)
tags <- c("[0-0.7]","[0.7-1.4]","[1.4-2]")
bc <- cut(aggs$paths_to_bas, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)
bcp <- ggplot(data=aggs, aes(x=resV, group=bc, fill=bc)) +
  geom_density(adjust=1.5, alpha=.4) + guides(fill=guide_legend(title="basal cons\n(scaled)"))

#pdf("###.pdf",width = 10, height = 10)
ggarrange(vul, tlp, dbp, bcp,
          labels = c("A", "B","C","D"),
          ncol = 2, nrow = 2)
#dev.off()

#vulnerability by threshold
ci <- function(n) qnorm(0.975)*sd(n)/sqrt(length(n))
out <- list()
for(i in 1:10){
  out[[i]] <- data.frame(int = paste(i*0.1, "to",i*0.1-0.1),mean = aggregate(sim[sim$tre<i*0.1&sim$tre>(i*0.1-0.1),]$resV,by=list(sim[sim$tre<i*0.1&sim$tre>(i*0.1-0.1),]$extinct_status),FUN=mean), ci = ci(sim[sim$tre<i*0.1&sim$tre>(i*0.1-0.1),]$resV))
}
  

dt <- do.call("rbind",out)
names(dt) <- c("interval","extinct_status","mean","ci")
dt$up <- dt$mean + dt$ci
dt$down <- dt$mean - dt$ci
dt$ord <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)

pdf("###.pdf")
plot(dt$mean~dt$ord,col=as.factor(dt$extinct_status),pch=19,xlab="threshold",ylab="resid vuln",xaxt='n', bty='L')
axis(side = 1, at = c(1,2,3,4,5,6,7,8,9,10), labels = c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1"),las=1)
arrows(x0 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, x1 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, y0 = dt$mean, y1 = dt$up,angle = 0,col=as.factor(dt$extinct_status))
arrows(x0 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, x1 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, y0 = dt$mean, y1 = dt$down,angle = 0,col=as.factor(dt$extinct_status))
par(xpd=TRUE)
legend(1,0.1,c("extant", "extinct"), pch = c(19,19), col = as.factor(unique(dt$extinct_status)))
dev.off()


sim2 <- readRDS("###data/simTL.rds")
  
  out <- list()
for(i in 1:10){
  out[[i]] <- data.frame(int = paste(i*0.1, "to",i*0.1-0.1),mean = aggregate(sim2[sim2$tre<i*0.1&sim2$tre>(i*0.1-0.1),]$vulnerability,by=list(sim2[sim2$tre<i*0.1&sim2$tre>(i*0.1-0.1),]$extinct_status),FUN=mean), ci = ci(sim2[sim2$tre<i*0.1&sim2$tre>(i*0.1-0.1),]$vulnerability))
}


dt <- do.call("rbind",out)
names(dt) <- c("interval","extinct_status","mean","ci")
dt$up <- dt$mean + dt$ci
dt$down <- dt$mean - dt$ci
dt$ord <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)

plot(dt$mean~dt$ord,col=as.factor(dt$extinct_status),pch=19,xlab="threshold",ylab="n plants remain",xaxt='n', bty='L')
axis(side = 1, at = c(1,2,3,4,5,6,7,8,9,10), labels = c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1"),las=1)
arrows(x0 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, x1 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, y0 = dt$mean, y1 = dt$up,angle = 0,col=as.factor(dt$extinct_status))
arrows(x0 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, x1 = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)+1, y0 = dt$mean, y1 = dt$down,angle = 0,col=as.factor(dt$extinct_status))
