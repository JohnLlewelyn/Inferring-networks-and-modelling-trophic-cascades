#make keep80 have similar diet categories to the Elton traits databases (i.e., invertebrates = aquatic and terrestrial; Fish instead of FishInv)
#Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M.M. and Jetz, W., 2014. EltonTraits 1.0: Species‐level foraging attributes of the world's birds and mammals: Ecological Archives E095‐178. Ecology, 95(7), pp.2027-2027.
#setwd("~/###")
#pp <- getwd() 
#get data
#keep80<-readRDS(paste(pp,"/data/keep80.rds",sep=""))

#Elton traits databases
mam<-read.csv(paste(pp,"/data/MamFuncDatEltonTraits.csv",sep=""))
ave<-read.csv(paste(pp,"/data/BirdFuncDatEltonTraits.csv",sep=""))

#cut to columns needed
cl <- c("Scientific","Diet.Inv","Diet.Vend","Diet.Vect","Diet.Vfish","Diet.Vunk","Diet.Scav",
        "Diet.Fruit","Diet.Nect","Diet.Seed","Diet.PlantO")
mam <- mam[,cl]
ave <- ave[,cl]
elt <- rbind(mam,ave)

#collapse the columns to 4 columns: plants, invertebrates, terrestrial vertebrates, fish
elt$eat.plant <- ifelse(elt$Diet.Fruit>0|elt$Diet.Nect>0|elt$Diet.Seed>0|elt$Diet.PlantO>0,"TRUE","FALSE")
elt$eat.invert <- ifelse(elt$Diet.Inv>0,"TRUE","FALSE")
elt$eat.fish <- ifelse(elt$Diet.Vfish>0,"TRUE","FALSE")
elt$eat.terr_vert <- ifelse(elt$Diet.Vend>0|elt$Diet.Vect>0|elt$Diet.Vunk>0|elt$Diet.Scav>0,"TRUE","FALSE")
elt <- elt[,c("Scientific","eat.plant","eat.invert","eat.fish","eat.terr_vert")]

#fix some scientific names
elt$Scientific[elt$Scientific=="Casmerodius albus"] <- "Ardea alba"
elt$Scientific[elt$Scientific=="Mesophoyx intermedia"] <- "Ardea intermedia"
elt$Scientific[elt$Scientific=="Alcedo azurea"] <- "Ceyx azureus"
elt$Scientific[elt$Scientific=="Larus novaehollandiae"] <- "Chroicocephalus novaehollandiae"
elt$Scientific[elt$Scientific=="Sterna nilotica"] <- "Gelochelidon nilotica"
elt$Scientific[elt$Scientific=="Sterna caspia"] <- "Hydroprogne caspia"
elt$Scientific[elt$Scientific=="Ixobrychus minutus"] <- "Ixobrychus dubius"
elt$eat.invert[elt$Scientific=="Ixobrychus dubius"] <- "TRUE"
elt$eat.terr_vert[elt$Scientific=="Ixobrychus dubius"] <- "TRUE"
elt$eat.fish[elt$Scientific=="Ixobrychus dubius"] <- "TRUE"
elt$Scientific[elt$Scientific=="Sterna bergii"] <- "Thalasseus bergii"
elt$Scientific[elt$Scientific=="Gallinula ventralis"] <- "Tribonyx ventralis"

#split keep80 into endotherms and ectotherms
kend <- keep80[keep80$Class=="Mammalia"|keep80$Class=="Aves",]
kect <- keep80[keep80$Class=="Amphibia"|keep80$Class=="Reptilia",]

#add updated diet columns to kend
kend_inc <- kend[kend$SciName%in%elt$Scientific,]
kend_inc[,c("eat.plant","eat.invert","eat.vert","eat.FishInv")] <- NULL
kend_inc <- merge(kend_inc,elt,by.x="SciName",by.y="Scientific")

#fix the species that don't have scientific names that line up with the Elton Traits database
kend_exc <- kend[!(kend$SciName%in%elt$Scientific),]
kend_exc1 <- kend_exc[kend_exc$eat.FishInv!="TRUE",]   #the ones that have FALSE for FishInv don't need new diet info, just fix some column names
names(kend_exc1)[names(kend_exc1)=="eat.vert"] <- "eat.terr_vert"
names(kend_exc1)[names(kend_exc1)=="eat.FishInv"] <- "eat.fish"
kend_exc1 <- kend_exc1[,names(kend_inc)]
#add back in the species that didn't have FishInv
kend_inc <- rbind(kend_inc,kend_exc1)
kend <- kend_inc

#now check kect
names(kect)[names(kect)=="eat.vert"] <- "eat.terr_vert"
names(kect)[names(kect)=="eat.FishInv"] <- "eat.fish"
kect$eat.terr_vert[kect$SciName=="Amphibolurus muricatus"] <- "FALSE"  #Fitzsimons, J.A. and Thomas, J.L., 2018. Diet and foraging strategies of the Jacky Lizard'Amphibolurus muricatus'. Australian Zoologist, 39(3), p.440.
kect$eat.terr_vert[kect$SciName=="Amphibolurus norrisi"] <- "FALSE"  #Geiser, F. and Learmonth, R.P., 1994. Dietary fats, selected body temperature and tissue fatty acid composition of agamid lizards (Amphibolurus nuchalis). Journal of Comparative Physiology B, 164(1), pp.55-61.
kect$eat.plant[kect$SciName=="Amphibolurus norrisi"] <- "FALSE"  #
kect$eat.terr_vert[kect$SciName=="Aprasia striolata"] <- "FALSE"  #Webb, J.K. and Shine, R., 1994. Feeding habits and reproductive biology of Australian pygopodid lizards of the genus Aprasia. Copeia, pp.390-398.
kect$eat.terr_vert[kect$SciName=="Ctenophorus pictus"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Ctenotus orientalis"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Ctenotus spaldingi"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Ctenotus uber"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Delma impar"] <- "FALSE" #Kutt, A.S., Coulson, G. and Wainer, J., 1998. Diet of the Striped Legless Lizard Delma impar (Squamata: pygopodidae) in a western (basalt) plains grassland, Victoria. Australian Zoologist, 30(4), pp.412-418.
kect$eat.terr_vert[kect$SciName=="Delma inornata"] <- "FALSE" #Patchell, F.C. and Shine, R., 1986. Food habits and reproductive biology of the Australian legless lizards (Pygopodidae). Copeia, pp.30-39.
kect$eat.fish[kect$SciName=="Delma inornata"] <- "FALSE"
kect$eat.fish[kect$SciName=="Emydura macquarii"] <- "TRUE"
kect$eat.terr_vert[kect$SciName=="Limnodynastes dumerilii"] <- "FALSE"
kect$mean.mass.kg[kect$SciName=="Limnodynastes peronii"] <- 0.03
kect$eat.terr_vert[kect$SciName=="Limnodynastes peronii"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Litoria peronii"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Neobatrachus pictus"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Neobatrachus sudellae"] <- "FALSE"
kect$eat.terr_vert[kect$SciName=="Pygopus lepidopodus"] <- "FALSE"  #Patchell, F.C. and Shine, R., 1986. Food habits and reproductive biology of the Australian legless lizards (Pygopodidae). Copeia, pp.30-39.
kect$eat.terr_vert[kect$SciName=="Tympanocryptis lineata"] <- "FALSE" #MacMillen, R.E., Augee, M.L. and Ellis, B.A., 1989. Thermal ecology and diet of some xerophilous lizards from western New South Wales. Journal of Arid Environments, 16(2), pp.193-201.
kect$eat.plant[kect$SciName=="Tympanocryptis lineata"] <- "TRUE"
kect$eat.terr_vert[kect$SciName=="Underwoodisaurus milii"] <- "FALSE"

kect <-kect[,names(kend)]
keep80 <- rbind(kend,kect)

#fix platypus - they rarely eat fish
keep80$eat.fish[keep80$SciName=="Ornithorhynchus anatinus"] <- "FALSE"
keep80$eat.fish[keep80$SciName=="Rhipidura leucophrys"] <- "FALSE"
#saveRDS(keep80,paste(pp,"/data/keep80_updated.rds",sep=""))


