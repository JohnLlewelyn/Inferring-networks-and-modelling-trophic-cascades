## PCA for Naracoorte network metrics
## CJA Bradshaw and J Llewelyn
## February 2020

## remove everything
rm(list = ls())

## libraries
library(Rcmdr)
library(vegan)
library(viridis)

## import data
#setwd("~/Documents/Papers/Palaeo/Sahul/Networks/data")
pcadat <- read.table("network_data_pca.csv", header=T, sep=",")

# create average pca first i.e., each species' mean for each variable
species.vec <- attr(table(pcadat$species), "names")
lspp <- length(species.vec)

mean.mat <- matrix(data=NA, ncol=dim(pcadat)[2]-1, nrow=lspp)
ext.vec <- rep(NA, lspp)
for (i in 1:lspp) {
  sp.dat <- subset(pcadat, species==species.vec[i])
  ext.vec[i] <- as.character(sp.dat$extinct[1])
  mean.mat[i,3:13] <- as.numeric(apply(sp.dat[4:14], MARGIN=2, mean, na.rm=T))
  print(i)
}
pca.mean <- as.data.frame(mean.mat)
pca.mean[,1] <- ext.vec
pca.mean[,2] <- species.vec
colnames(pca.mean) <- colnames(pcadat[,-1])

## correlation matrix
cor.mat <- round(cor(na.omit(pca.mean[,-c(1,2)]), method="spearman"), 3)
cor.mat
cor.mat <- cor(na.omit(pca.mean[,-c(1,2)]), method="pearson")
round(cor.mat*cor.mat, 3)
pairs(pca.mean[,-c(1,2)], pch=19, cex=0.7)

# remove eigenvector_centrality, closeness_in, eccentricity_in, coreness_in, closeness_out (as redundant)
pca.mean.clip <- pca.mean[,-c(5,6,7,9,10)]
cor.mat <- round(cor(na.omit(pca.mean.clip[,-c(1,2)]), method="spearman"), 3)
cor.mat

# PCA
z <- rda(pca.mean.clip[,3:8], data = pca.mean, scale=T, center=T, na.action="na.omit")
plot(z, type="text")
summ.out <- summary(z)
pca.out <- summ.out$sites

# merge
pca.out2 <- data.frame(pca.mean, pca.out)
pca.ext <- subset(pca.out2, extinct=="ext")
pca.srv <- subset(pca.out2, extinct=="surv")

# extinct vs. survived
plot(z, type="none", ylim=c(min(pca.out2$PC2),max(pca.out2$PC2)), xlim=c(min(pca.out2$PC1),max(pca.out2$PC1)))
points(pca.ext$PC1, pca.ext$PC2, pch="E", col="red")
points(pca.srv$PC1, pca.srv$PC2, pch="S", col="blue")

## ggbiplot
library(devtools)
library(ggbiplot)
library(factoextra)

#fix some labels/categories in the dfs
colnames(pca.mean)[colnames(pca.mean)=="extinct"] <- "status"
pca.mean$status <- ifelse(pca.mean$status=="ext","extinct", "extant")

colnames(pca.mean.clip)[colnames(pca.mean.clip)=="degree_out"] <- "degree out"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="eccentricity_out"] <- "eccentricity out"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="coreness_out"] <- "coreness out"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="degree_in"] <- "degree in"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="pagerank"] <- "pageRank"


#run pca
z2 <- prcomp(pca.mean.clip[,3:8], center=T, scale=T)

#colours for plots
cols <- ifelse(pca.mean$status=="extinct","green","yellow")

#biplot
pca_plot <- fviz_pca_biplot(z2,axes = c(1, 2), label="var", labelsize = 8, col.ind = as.factor(pca.mean$status), addEllipses=TRUE, 
                            ellipse.level=0.95, legend.title = "status",xlim=c(-7.5,7.5),col.var = "black",repel=T,palette = c("#CC4678FF","#0D0887FF")) + theme_classic() + ggtitle("") + 
                            theme(axis.title = element_text(size=28), axis.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=18),
                            legend.key=element_rect(fill='white'),legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
                            legend.justification = "top") + guides(colour = guide_legend(override.aes = list(size=1, stroke=1))) + 
                            annotate("text", -7.8, Inf, label = "a)", hjust = 0, vjust = 1,cex=10,fontface=2)#

#fix categories
pcadat$extinct <- ifelse(pcadat$extinct=="ext", "extinct", "extant")

preds <- aggregate(pcadat$degree_out,by=list(pcadat$species,pcadat$extinct),FUN=mean)
colnames(preds) <- c("species", "status", "score")

box_fun <- function(data, y_lab, x_lab, panel) {
  ggplot(data, aes(x=status, y=score, fill=status)) + 
    geom_boxplot(notch = F, alpha=0.5) + #scale_fill_manual(values=c("deepskyblue3","orange")) +
    #stat_summary(fun.y=mean, geom="point", shape=18, size=4) +
    theme_classic() + theme(axis.title = element_text(siz=28), axis.text = element_text(size=20)) + ylab(y_lab) +xlab(x_lab) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    annotate("text", -Inf, Inf, label = panel, hjust = 0, vjust = 1,cex=10,fontface=2)+scale_fill_manual(values=c("#CC4678FF","#0D0887FF"))
}
pred_box <- box_fun(preds,"number of predators","species' status\n(end Pleistocene)", "  b)")

#Combine plots
library(gridExtra)
library(ggpubr)

pdf("~/Dropbox/trophic links/Docs/Figures/PCA and number of predators.pdf", height=8,width=18)
ggarrange(pca_plot, pred_box, ncol = 2, nrow = 1, align = "h", heights = c(1,1))
dev.off()


#plots for supp mat
pdf("~/Dropbox/trophic links/Docs/supp material/keepers/figures/biplot dimensions 3 and 4.pdf")
fviz_pca_biplot(z2,axes = c(3, 4), label="var", habillage=as.factor(pca.mean$status), addEllipses=TRUE, 
                ellipse.level=0.95, legend.title = "status",xlim=c(-7.5,7.5),col.var = "black",repel=T) + ggtitle("") + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=18), legend.title = element_text(size=20), legend.text = element_text(size=18),
        legend.key=element_rect(fill='white'),legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.justification = "top") + guides(colour = guide_legend(override.aes = list(size=1, stroke=1)))  
dev.off()

pdf("~/Dropbox/trophic links/Docs/supp material/keepers/figures/contributions to PC1.pdf")
fviz_contrib(z2, choice = "var", axes = 1)
dev.off()

#eigenvalues/variances of dimensions
pdf("~/Dropbox/trophic links/Docs/supp material/keepers/figures/variance of dimension (PCs).pdf")
fviz_eig(z2, addlabels = TRUE)
dev.off()