a<-readRDS('coextinction sim/data/x18lump.rds')


for (i in 1:1000){
	net<-a[i][[1]]  
	#for (cat in c('invertebrates','plants','FishInv')){
	#	net<-net[net[,1]!=cat,]
	#	net<-net[net[,2]!=cat,]
	#	}
	write.table(net,paste0('coextinction sim/data/nets_18/',i,'.csv'),quote=F,sep=',')}


a<-readRDS('coextinction sim/data/x80nets.rds')
for (i in 1:1000){
	net<-a[i][[1]]
	#for (cat in c('invertebrates','plants','FishInv')){
	#	net<-net[net[,1]!=cat,]
	#	net<-net[net[,2]!=cat,]
	#	}
	write.table(net,paste0('coextinction sim/data/nets_80_ok/',i,'.csv'),quote=F,sep=',')}



a<-readRDS('coextinction sim/data/xIntlump.rds')
for (i in 1:1000){
	net<-a[i][[1]]
	for (cat in c('invertebrates','plants','FishInv')){
		net<-net[net[,1]!=cat,]
		net<-net[net[,2]!=cat,]
		}
	write.table(net,paste0('coextinction sim/data/nets_int/',i,'.csv'),quote=F,sep=',')}



