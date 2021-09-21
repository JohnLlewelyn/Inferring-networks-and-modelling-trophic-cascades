#get full networks
#save as igraph networks
#extract network position metrics

filenames <- list.files("###", pattern="*.csv", full.names=TRUE)
nts <- lapply(filenames, read.csv)

#igraph objects
library(igraph)
library(NetIndices)
nts_basic <- lapply(nts, graph_from_data_frame,directed=TRUE)

#list of vertebrates
keep80<-readRDS("~/Dropbox/trophic links/Docs/supp material/keepers/data/keep80.rds")
extinct <- keep80$SciName[keep80$Status=="Extinct"]

#get the following: , pageRank, betweenness centrality, eigenvector centrality, closeness centrality (in), 
#coreness (in), degree (in), eccentricity (in), closeness centrality (out), coreness (out), degree (out), 
#and eccentricity (out), trophic level 

#make function that adapts pageRank function for our purpose
page_rank_vec <- function(x) {
  dt <- page_rank(x)$vector
  dt <- data.frame(species = names(dt),pageRank = dt)
 # colnames(dt) <- c("species","pageRank")
  return(dt)}
#adapt betweenness function
betweenness_vec <- function(x) {
  dt <- betweenness(x,directed = TRUE)
  dt <- data.frame(species = names(dt), betweenness = dt)
  return(dt)
}
#eigenvector centrality adapted
eigen_vec <- function (x) {
  dt <- eigen_centrality(x, directed = TRUE)$vector
  dt <- data.frame(species = names(dt), eigenvector = dt)
  return(dt)
}
#closeness centrality (in) adapted
closeIn_vec <- function (x) {
  dt <- closeness(x, mode = c("in"))#$vector
  dt <- data.frame(species = names(dt), closeIn = dt)
  return(dt)
}
#closeness centrality (out) adapted
closeOut_vec <- function (x) {
  dt <- closeness(x, mode = c("out"))#$vector
  dt <- data.frame(species = names(dt), closeOut = dt)
  return(dt)
}
#coreness (In) adapted
corenessIn_vec <- function (x) {
  dt <- coreness(x, mode = c("in"))#$vector
  dt <- data.frame(species = names(dt), corenessIn = dt)
  return(dt)
}
#coreness (out) adapted
corenessOut_vec <- function (x) {
  dt <- coreness(x, mode = c("out"))#$vector
  dt <- data.frame(species = names(dt), corenessOut = dt)
  return(dt)
}

#degree (in) adapted
degreeIn_vec <- function (x) {
  dt <- degree(x, mode = c("in"))
  dt <- data.frame(species = names(dt), degreeIn = dt)
  return(dt)
}

#degree (out) adapted
degreeOut_vec <- function (x) {
  dt <- degree(x, mode = c("out"))
  dt <- data.frame(species = names(dt), degreeOut = dt)
  return(dt)
}

#eccentricity (in) adapted
eccentricityIn_vec <- function (x) {
  dt <- eccentricity(x, mode = c("in"))
  dt <- data.frame(species = names(dt), eccentricityIn = dt)
  return(dt)
}

#eccentricity (out) adapted
eccentricityOut_vec <- function (x) {
  dt <- eccentricity(x, mode = c("out"))
  dt <- data.frame(species = names(dt), eccentricityOut = dt)
  return(dt)
}

pos_fun <- function(fun,n_list) {
  aa <- lapply(n_list, fun)
  ab <- do.call("rbind",aa)
  ac <- ab[ab$species%in%keep80$SciName,]
  ad <- aggregate(ac[,2],by=list(ac$species),FUN=mean)
  names(ad) <- c("species",names(ac)[2])
  return(ad)
}

#get the position metrics
pageRank_nara <- pos_fun(page_rank_vec, nts_basic)
betweenness_nara <- pos_fun(betweenness_vec, nts_basic)
eigenvector_nara <- pos_fun(eigen_vec, nts_basic)
closeIn_nara <- pos_fun(closeIn_vec, nts_basic)
closeOut_nara <- pos_fun(closeOut_vec, nts_basic)
corenessIn_nara <- pos_fun(corenessIn_vec, nts_basic)
corenessOut_nara <- pos_fun(corenessOut_vec, nts_basic)
degreeIn_nara <- pos_fun(degreeIn_vec, nts_basic)
degreeOut_nara <- pos_fun(degreeOut_vec, nts_basic)
eccentricityIn_nara <- pos_fun(eccentricityIn_vec, nts_basic)
eccentricityOut_nara <- pos_fun(eccentricityOut_vec, nts_basic)

#Get the trophic level data
#get the simulations results
sim <- read.csv("results_30_06_21.csv",header=T)
trophicLevel_nara <- aggregate (sim$mean_tl,by=list(sim$species),FUN=mean)
names(trophicLevel_nara) <- c("species","trophicLevel")

#all degOut data
all_d <- lapply(nts_basic,degreeOut_vec)
#stick them together with model id
all_d <- bind_rows(all_d, .id = "net_id")
all_d$net_id <- as.numeric(all_d$net_id)
#keep only verts
all_d <- all_d[all_d$species%in%keep80$SciName,]

#add if extinct
all_d$extinct <- ifelse(all_d$species%in%extinct, "extinct","extant")
#get all the diffs
df_fun <- function(x) data.frame(cbind(aggregate(x$degreeOut, by=list(x$extinct),mean)[1,2],aggregate(x$degreeOut, by=list(x$extinct),mean)[2,2]))
dfs <- ddply(all_d,~net_id,df_fun)
min(dfs$X1-dfs$X2) #in all model, extant species had more predators than extinct species

#merge and save position metrics
library(tidyverse)
tog <- list(pageRank_nara, betweenness_nara, eigenvector_nara, closeIn_nara, closeOut_nara, corenessIn_nara, corenessOut_nara, degreeIn_nara, degreeOut_nara, eccentricityIn_nara, eccentricityOut_nara, trophicLevel_nara) %>% reduce(left_join, by = "species")
#add extinction stats# $extinct
tog$extinct <- ifelse(tog$species%in%extinct, "extinct","extant")

#save network position metrics
saveRDS(tog,"network_position_metrics.rds")

#Next, run correlations and pca on position metrics #pca_position.R