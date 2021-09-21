#Vetting megafauna

#Fix grey kangaroo
keep80$SciName[keep80$SciName=="Macropus fuliginosus/giganteus/titan"] <- "Macropus fuliginosus"	
keep80$mean.mass.kg[keep80$SciName=="Macropus fuliginosus"] <- 25.6  #From Elton traits
keep80$min.admass[keep80$SciName=="Macropus fuliginosus"] <- NA
keep80$max.admass[keep80$SciName=="Macropus fuliginosus"] <- NA
keep80$Status[keep80$SciName=="Macropus fuliginosus"] <- "Extant"
keep80$ComName[keep80$SciName=="Macropus fuliginosus"] <- "Western grey kangaroo"

#fix Latigallina naracoortensis name
keep80$SciName[keep80$SciName=="Leipoa gallinacea"] <- "Latigallina naracoortensis"

#megafauna to remove
rem <- c("Palorchestes azael",
         "Simosthenurus pales")

keep80 <- keep80[!(keep80$SciName%in%rem),]

#fix thylacine mass and diet
keep80$mean.mass.kg[keep80$SciName=="Thylacinus cynocephalus"] <- 23.35
#keep80$eat.invert[keep80$SciName=="Thylacinus cynocephalus"] <- "TRUE"  #according to Elton traits, coyotes and devils don't eat invertebrates, so keep invertebrates as false
