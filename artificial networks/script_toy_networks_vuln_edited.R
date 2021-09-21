### Code for analyzing extinctions in food webs using Bayesian networks. Code is written by Anna Eklöf, Si Tang and Stefano Allesina. 
  # Checked 2014-04-30 for R version 3.0.2.
  # Checked 2016-05-05 for R version 3.2.4.

### REFERENCE: Eklöf, Anna, Si Tang, and Stefano Allesina. "Secondary extinctions in food webs: a Bayesian network approach." Methods in Ecology and Evolution 4.8 (2013): 760-770.

### Please cite and acknowledge properly!
# source("http://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
# 
# install.packages('gRain')


source("BN-BayesNet.R")
source("BN-Likelihoods.R")
source("BN-DFS.R")

LaunchAnalysis <- function(Adjacency,
                           Extant,
                           ProbExtinction,
                           Functional = "topo",
                           Alpha = 0.0,
                           Beta = 0.0,
                           Label){
    ## Wrapper for the analysis.
    ##############################
    ## INPUT
    ##############################
    # Adjacency: square, binary matrix
    #            Adjacency[i,j] = 1 means that j consumes i
    # Extant: observed presence of the species
    #         binary, Reps x S matrix, where
    #         Reps is the number of observations
    #         and each row is a different replicate
    #         marking whether the species was present (1) or absent (0)
    # ProbExtinction: vector of length S, measuring
    #                 the probability that each species goes extinct due
    #                 to external causes
    # Functional: functional form, as in the paper.
    #             Options are "topo", "linear", and "nonlinear"
    # Alpha, Beta: parameters for nonlinear functional response
    # Label: label for the system
    #
    # Results returned as "results.Rdata" stored in "FW"
    ##############################

    # First, we need to make the adjacency matrix acyclic
    # For this, we use BN-DFS.R
    A <- DFS(A=Adjacency)
    if (sum(A) < sum(Adjacency)){
        print("The network is not acyclic, and thus some links were removed!")
    }
    ## Call the function
    FW <- GetBayesNetMarginals(A,
                              Extant,
                              ProbExtinction,
                              Functional,
                              Alpha,
                              Beta,
                              Label)
    FW$Functional <- NULL
    FW$OriginalAdjacency <- Adjacency
    FW$RemovedLinks <- Adjacency - A
    
    save(FW, file='results.RData')
    return(FW)
}

# 
# Adjacency<-as.matrix(read.csv('mat.csv',header=F))
# S<-nrow(Adjacency)
# diag(Adjacency) <- 0
# Reps <- 50
# Extant <- (matrix(runif(Reps * S), Reps, S) < 0.2) * 1 #matrix(rep(1,S),1)#matrix(rep(1,Reps * S), Reps, S) 
# res<-LaunchAnalysis(Adjacency, Extant, rep(0.1, S), "topo", Label="FW1")
# #print(LaunchAnalysis(Adjacency, Extant, rep(0.1, S), "nonlinear", 0.5, 3, "FW1"))
# m<-res$M
# vuln<-res$MarginalsExtant
# vuln<-1-vuln
# 
# v_topo<-c(0.748591,
#            0.688881,
#            0.7664486666666668,
#            0.6423376666666666,
#            0.6190663333333334,
#            0.5665310000000001,
#            0.566428)
# 
# v_an<-c(0.18461010582681336,
#          0.13531924899999992,
#          0.13531924899999992,
#          0.10899999999999999,
#          0.10899999999999999,
#          0.09999999999999998,
#          0.09999999999999998)
# 
# pdf('comparison_plots.pdf',width=15,height=5)
# par(mfrow=c(1,3))
# plot(vuln,v_topo,xlab='bayes',ylab='simulated',pch=16,cex=2)
# #lines(vuln,v_topo)
# plot(v_an,v_topo,xlab='analytical',ylab='simulated',pch=16,cex=2)
# #lines(v_an,v_topo)
# plot(vuln,v_an,xlab='bayes',ylab='analytical',pch=16,cex=2)
# #lines(vuln,v_an)
# dev.off()
# 

library(igraph)
library(NetIndices)

tre = 0 #set threshold
pdf(paste('bottom_up_bayesian_comp_thresh_new',tre,'.pdf',sep=''),height=4,width=12)
#pdf('equal_vuln_bayesian.pdf',height=5,width=10)

par(mfrow=c(1,3))
tl_vs_vuln<-c()
for (mat_n in 0:195){
  #Adjacency <- (matrix(runif(S*S), S, S) < 0.5) * 1
  Adjacency<-as.matrix(read.csv(paste0('./mats/',mat_n,'_mat.csv'),header=F))
  #get the graph/network details
  deets <- unlist(GenInd(Adjacency))
  saveRDS(deets,paste0('./toy_network_mets/',mat_n,'_metrics.rds'))
  basal<-which(colSums(Adjacency)==0)
  non_basal<-which(colSums(Adjacency)>0)
  g<-graph_from_adjacency_matrix(Adjacency)
  d<-distances(g,mode='IN')
  if (length(basal)>1){tl<-apply(d[,basal],FUN='min',1)} else {tl<-as.numeric(d[,basal])} #get trophic level of each node
  tl<-tl+1
  diet_breadth <- colSums(Adjacency) #number of 'in' links
  basal_cons <- table(unlist(ego(graph=graph_from_adjacency_matrix(Adjacency), order=20, nodes=basal,mode="out"))) #number of basal nodes each node is directly or indirectly connected to
  basal_cons[basal] <- 0
  S<-nrow(Adjacency)
  diag(Adjacency) <- 0
  Reps <- 50
  Extant <- (matrix(runif(Reps * S), Reps, S) < 0.2) * 1 #matrix(rep(1,S),1)#matrix(rep(1,Reps * S), Reps, S) 
  #base_vuln<-rep(0.1,S)
  base_vuln<-rep(0.0,S)  #base vulnerability for all nodes, set to 0
  base_vuln[basal]<-0.1  #give basal nodes a base vulnerability of 0.1
  res<-LaunchAnalysis(Adjacency, Extant, base_vuln, "topo", Label="FW1")  #To run Eklöf's vulnerability measure
  #print(LaunchAnalysis(Adjacency, Extant, rep(0.1, S), "nonlinear", 0.5, 3, "FW1"))
  m<-res$M
  vuln<-res$MarginalsExtant #the vulnerabilities for each node using Eklöf's method
  coe_<-c()
  all_res <- data.frame()
  for (sim_rep in 1:1000){
    mat_ext<-Adjacency
    sc<-1 #start with removing 1 basal node
    ext<-basal
    rbasal<-sample(basal) #shuffle order of basal nodes (for extinction order)
    for (b in rbasal){
      mat_ext[b,]<-0  #make row representing extinct basal resource all '0'
      coe<-setdiff(which(colSums(mat_ext)<=diet_breadth*tre),ext) #which nodes don't have above the threshold 'in' links now and are not basal resources themselves
      if (length(coe)>0){coe_<-rbind(coe_,cbind(coe,rep(sc,length(coe))))} #add information on how many primary extinctions before this secondary extinction(s)
      ext<-c(ext,coe)
      while (length(coe)>0){  #further trophic cascades
        mat_ext[coe,]<-0
        mat_ext[,coe]<-0
        coe<-setdiff(which(colSums(mat_ext)<=diet_breadth*tre),ext) #following the previous secondary extinctions, which, if any, more nodes lose all their resources?
        ext<-c(ext,coe)
        if (length(coe)>0){coe_<-rbind(coe_,cbind(coe,rep(sc,length(coe))))}
      }
      sc<-sc+1
    }}
  coe_m<-aggregate(coe_[,2]~coe_[,1],FUN='mean')   #coe_ = matrix with order each (non-basal) node went extinct in 1000 simulations, coe_m aggregated across sims
  coe_m[,2]<-1-coe_m[,2]/length(basal) #scale 
  #coe_m[,2]<-1-sapply(coe_m[,2],pbinom,size=length(basal),prob=0.1)
  sim_vuln<-rep(0,nrow(mat_ext)) #assign vulnerability to basal and higher trophic levels
  sim_vuln[coe_m[,1]]<-coe_m[,2]
  sim_vuln[basal]<-0.1
  #write.table(m,paste0('./results/',rep,'_mat.csv'),quote=F,sep=',',row.names=F,col.names=F)
  # write.table(vuln,paste0('./vulns/',rep,'_vuln.csv'),quote=F,sep=',',row.names=F,col.names=F)
  tl_vs_vuln<-rbind(tl_vs_vuln,cbind(tl,1-vuln))
  tl_db_Bay_sim <- data.frame(cbind(tl,diet_breadth,basal_cons,1-vuln,sim_vuln))
  names(tl_db_Bay_sim)[names(tl_db_Bay_sim)=="V4"] <- "bay_vuln"
  #saveRDS(tl_db_Bay_sim,paste0('./toy_vulns/',mat_n,tre,'orig_vuln.rds'))
  saveRDS(tl_db_Bay_sim,paste0('./toy_vulns_new/',mat_n,'orig_vuln.rds'))
  print (mat_n)
  plot(tl,1-vuln,xlab='trophic level',ylab='vulnerability bayes',las=1,cex.axis=1.2,cex.lab=1.2,pch=16,main=mat_n,cex=3)
 r2<-round(cor(sim_vuln[non_basal],1-vuln[non_basal],method='spearman'),2)
 #rank_sim<-rank(round(sim_vuln[non_basal],3),ties.method='max')
 #rank_bay<-rank(round(1-vuln[non_basal],3),ties.method='max')
 rank_sim<-sim_vuln[non_basal]
 rank_bay<-1-vuln[non_basal]
 plot(rank_sim,rank_bay,xlab='vulnerability simulations',
     ylab='vulnerability bayes',las=1,cex.axis=1.2,cex.lab=1.2,pch=16,main=paste("Spearman's rho = ",r2),cex=3,
    )
  V(g)$label<-round(1-vuln,3)
  x_pos<-c()
  sc<-rep(0,max(tl))
  for (i in seq(1:S)){
    ntl<-sum(tl==tl[i])
    to_add<-((1/ntl)*sc[tl[i]])
    x_pos<-c(x_pos,to_add)
    sc[tl[i]]<-sc[tl[i]]+1
     }
  lay<-cbind(x_pos,tl)
  V(g)$color<-'white'
  V(g)$size<-30
  V(g)$shape<-'square'
  plot(g,layout=lay)
}  
dev.off()

#plot mean vulnerability using Eklöf's method
library(viridis)

b<-aggregate(tl_vs_vuln[,2]~tl_vs_vuln[,1],FUN='mean')
b_sd<-aggregate(tl_vs_vuln[,2]~tl_vs_vuln[,1],FUN='sd')
pdf('average_tl_vs_bayes_bottom_up_new.pdf',height=5,width=5)

plot(b[,1],b[,2],type='l',xlab='trophic level',ylab='vulnerability bayes',ylim=c(0,max(b[,2]+b_sd[,2])),
     las=1,cex.axis=1.2,cex.lab=1.2,lwd=2)
points(b[,1],b[,2],pch=16,cex=2)
ll<-b[,2]-b_sd[,2]
ll[which(ll<0)]<-0
arrows(b[,1], ll, b[,1], b[,2]+b_sd[,2], length=0.05, angle=90, code=3)
dev.off()


