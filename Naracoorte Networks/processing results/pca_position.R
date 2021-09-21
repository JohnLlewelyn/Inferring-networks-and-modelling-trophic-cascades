#PCA on position metrics, uses 'tog' from networks and positions metrics.R
library(Rcmdr)
library(vegan)
library(viridis)

setwd("###")
tog <- readRDS("network_position_metrics.rds")

#Check correlation between variables, remove those that are too highly correlated (>0.8 spearman)
res <- cor(tog[,2:(ncol(tog)-1)],method = "spearman")
View(round(res, 3))
pairs(tog[,-c(1,14)], pch=19, cex=0.7)

tog_cut <- tog[,-which(names(tog) %in% c("closeIn","eigenvector","eccentricityIn","trophicLevel", "corenessOut","degreeOut"))]
res <- cor(tog_cut[,2:(ncol(tog_cut)-1)],method = "spearman")
res

pca.mean <- tog
pca.mean.clip <- tog_cut

#PCA
data.pca <- princomp(pca.mean.clip[,-c(1,8)], cor=TRUE)
summary(data.pca, loadings=TRUE)

library(klaR)
z <- rda(pca.mean.clip[,2:(ncol(pca.mean.clip)-1)], data = pca.mean, scale=T, center=T, na.action="na.omit")
plot(z, type="text")
summ.out <- summary(z)
pca.out <- summ.out$sites

#merge
pca.out2 <- data.frame(pca.mean, pca.out)
pca.ext <- subset(pca.out2, extinct=="extinct")
pca.srv <- subset(pca.out2, extinct=="extant")

# extinct vs. survived
plot(z, type="none", ylim=c(min(pca.out2$PC2),max(pca.out2$PC2)), xlim=c(min(pca.out2$PC1),max(pca.out2$PC1)))
points(pca.ext$PC1, pca.ext$PC2, pch="E", col="red")
points(pca.srv$PC1, pca.srv$PC2, pch="S", col="blue")

## ggbiplot
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(factoextra)

#fix some labels/categories in the dfs
colnames(pca.mean)[colnames(pca.mean)=="extinct"] <- "status"
#pca.mean$status <- ifelse(pca.mean$status=="ext","extinct", "extant")

#colnames(pca.mean.clip)[colnames(pca.mean.clip)=="degreeOut"] <- "degree out"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="eccentricityOut"] <- "eccentricity out"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="corenessIn"] <- "coreness in"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="degreeIn"] <- "degree in"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="pageRank"] <- "pageRank"
colnames(pca.mean.clip)[colnames(pca.mean.clip)=="closeOut"] <- "closeness out"

#run pca
z2 <- prcomp(pca.mean.clip[,2:(ncol(pca.mean.clip)-1)], center=T, scale=T)

#colours for plots
cols <- ifelse(pca.mean$status=="extinct","red","green3")

#biplot
pca_plot <- fviz_pca_biplot(z2,alpha.var=0.5, axes = c(1, 2), label="var",ylab="", xlab="", labelsize = 8, col.ind = as.factor(pca.mean$status), 
                            addEllipses=TRUE, ellipse.level=0.95, legend.title = "status",xlim=c(-9.5,9.5),col.var = "black",repel=T,
                            palette = c("green3", "red")) + 
  theme_classic() + ggtitle("") + theme(axis.title = element_text(size=28), axis.text = element_text(size=20), legend.title = element_text(size=20),
                            legend.text = element_text(size=18), legend.key=element_rect(fill='white'), legend.background = element_rect(linetype = 1, 
                            size = 0.5, colour = 1), legend.position = c(0.9,0.9)) + 
  guides(colour = guide_legend(override.aes = list(size=1, stroke=1))) + 
  annotate("text", -10, Inf, label = "A", hjust = 0, vjust = 1,cex=10,fontface=2) + xlab("dimension 1 (44.1%)") + ylab("dimension 2 (23%)")#

#Plot number of preds
preds <- tog[,c("species","extinct","degreeOut")]
names(preds) <- c("species","status","score")
#average number of predators
aggregate(preds$score,by=list(preds$status),FUN=mean)
#diff in number of predators across all models


box_fun <- function(data, y_lab, x_lab, panel) {
  ggplot(data, aes(x=status, y=score, fill=status, col =status)) + 
    geom_violin(alpha = 0.5) + 
    geom_point(shape = 21,size=1., position = position_jitterdodge(0.5),alpha=0.5) +
    scale_x_discrete(name = "status\n(end Pliestocene)") +
    theme_classic() +
    ylab("number of predators") +
    theme(plot.title = element_text(size = 14),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 28),   #axis.title = element_text(size=28), axis.text = element_text(size=20),
          legend.position = "none")  +
    scale_fill_manual(values=c("green3","red")) +
    scale_colour_manual(values=c("green3","red")) +
    annotate("text", -Inf, Inf, label = panel, hjust = 0, vjust = 1,cex=10,fontface=2)
  }

pred_box <- box_fun(preds,"number of predators","species' status\n(end Pleistocene)", "  B")



#Plot number of resources
food <- tog[,c("species","extinct","degreeIn")]
names(food) <- c("species","status","score")

box_fun <- function(data, y_lab, x_lab, panel) {
  ggplot(data, aes(x=status, y=score, fill=status)) + 
    geom_boxplot(notch = F, alpha=0.5) + #scale_fill_manual(values=c("deepskyblue3","orange")) +
    #stat_summary(fun.y=mean, geom="point", shape=18, size=4) +
    theme_classic() + theme(axis.title = element_text(siz=28), axis.text = element_text(size=20)) + ylab(y_lab) +xlab(x_lab) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    annotate("text", -Inf, Inf, label = panel, hjust = 0, vjust = 1,cex=10,fontface=2)+scale_fill_manual(values=c("green3","red"))
}
food_box <- box_fun(food,"number of resource nodes","species' status\n(end Pleistocene)", "  C")

#Combine plots
library(gridExtra)
library(ggpubr)

pdf("###.pdf", height=8,width=18)
ggarrange(pca_plot,NULL, pred_box, nrow = 1, widths = c(1,0.05,1), align = "h", heights = c(1,1))
dev.off()

pdf("###.pdf", height=8,width=22)
ggarrange(pca_plot,NULL, pred_box, NULL, food_box, nrow = 1, widths = c(1,0.05,1,0.05,1), align = "h", heights = c(1,1,1))
dev.off()

#plots for supp mat
pdf("###.pdf")
fviz_pca_biplot(z2,axes = c(3, 4), label="var", habillage=as.factor(pca.mean$status), addEllipses=TRUE, 
                ellipse.level=0.95, legend.title = "status",xlim=c(-7.5,7.5),col.var = "black",repel=T,palette = c("green3", "red")) + ggtitle("") + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=18), legend.title = element_text(size=20), legend.text = element_text(size=18),
        legend.key=element_rect(fill='white'),legend.background = element_rect(linetype = 1, size = 0.5, colour = 1), 
        legend.justification = "top") + guides(colour = guide_legend(override.aes = list(size=1, stroke=1)))  
dev.off()

pdf("###.pdf")
fviz_contrib(z2, choice = "var", axes = 1)
dev.off()

pdf("###.pdf")
fviz_contrib(z2, choice = "var", axes = 2)
dev.off()

#eigenvalues/variances of dimensions
pdf("###.pdf")
fviz_eig(z2, addlabels = TRUE)
dev.off()