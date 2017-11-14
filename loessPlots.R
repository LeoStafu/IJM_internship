

#function for plot smooth curves for all lines and compare
#pos1 = 564424
#pos2 = 2253652


multiLineViewer = function(RTdf, pos1, pos2, cellLines = c("293T", "HeLa", "K562", "MRC5C", "RKO", "U2OS"), loessSpan = 1200000 , threshold = 0.8){

	colCell = rainbow(length(cellLines))
	nbCell = 1

	#If pos not found in data, the smaller range containing the two given pos will be chosen
	x1 = max(which(RTdf$Position <= pos1))
	if(x1 == -Inf){ x1 = 1 }
	x2 = min(which(RTdf$Position >= pos2))
	if(x2 == Inf){ x2 = dim(RTdf)[1] }
	
	#subDataFrame of the window
	RTdf = RTdf[x1:x2,]

	plot(NULL,xlim=c(pos1,pos2), ylim=c(-3,3), 
		ylab="Timing (Z-score)", xlab="Position on chromosome (bp)",
		main=paste("Replication timing on",RTdf[1,"Chromosome"]), 
		sub=paste("Position (arb. unit)",x1,"to",x2))

	for(cell in cellLines){
		#Here Loess computed on the given region beetween pos1 et pos2 arguments
		loesschr = loess(RTdf[,cell] ~ RTdf$Position, span = loessSpan/(max(RTdf$Position)-min(RTdf$Position)))
		
		lines(loesschr$fitted~RTdf$Position, col=colCell[nbCell], lwd = 3)
		nbCell = nbCell +1
	}

	legend("bottomleft", legend = as.character(cellLines), col=colCell, lwd  =3)

	#Add legend with result from other function
	conservCell = conservAtPos(rawRT_to_loessRT(RTdf), pos1, pos2, cellLines, threshold)
	legend("topright",inset = 0.02, legend =conservCell, bty="n", title=paste("Threshold=", threshold), lty=1, lwd = 3, col=colCell[cellLines %in% conservCell])

	abline(h=0, col="black", lty = "dotted")
}


#Transform raw RT data in loess RT
rawRT_to_loessRT = function(RTdf, loessSpan = 1200000, cellLines = c("293T", "HeLa", "K562", "MRC5C", "RKO", "U2OS")){

	for(cell in cellLines){
		#Here Loess computed on the whole chromosome
		loesschr = loess(RTdf[,cell] ~ RTdf$Position, span = loessSpan/(max(RTdf$Position)-min(RTdf$Position)))
		RTdf[,cell] = loesschr$fitted
	}
	return(RTdf)
}






RTchr = RTchr1
size = dim(RTchr[which(RTchr$Centro == 1), ])[1]
RTchr = RTchr[1:size,]
range1 = 1
range2 = 100

pdf("test.pdf")
while(range2 < dim(RTchr)[1]){

	multiLineViewer(RTchr, RTchr[range1,"Position"], RTchr[range2,"Position"], threshold = 0.80)

	range1 = range1 + 50
	range2 = range2 + 50
	if(range2 > dim(RTchr)){
		range2  = dim(RTchr)
	}
}
dev.off()










range = 100
plot(NULL, xlim=c(1,1.1), ylim=c(-0.1,1.2))
points(y = loessChr1[1:range,"MRC5C"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "lightblue", pch = 16, cex = 0.7)
points(y = loessChr1[1:range,"(range-1)3T"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "red", pch = 16, cex = 0.7)
points(y = loessChr1[1:range,"RKO"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "blue", pch = 16, cex = 0.7)
points(y = loessChr1[1:range,"HeLa"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "yellow", pch = 16, cex = 0.7)
points(y = loessChr1[1:range,"K562"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "green", pch = 16, cex = 0.7)
points(y = loessChr1[1:range,"U2OS"],x = seq(from=1, to = 1.1, by = 0.1/(range-1)), col = "purple", pch = 16, cex = 0.7)











RTchr = RTchr1
size = dim(RTchr[which(RTchr$Centro == 1), ])[1]
RTchr = RTchr[1:size,]
range1 = 1
range2 = 100

res= NULL

while(range2 < dim(RTchr)[1]){

	x  =hclust(dist(t(loessChr1[range1:range2,3:8]), method="euclidian"), method = "complete")$height
	
	res = c(res, sort(x))

	range1 = range1 + 50
	range2 = range2 + 50
	if(range2 > dim(RTchr)){
		range2  = dim(RTchr)
	}
}


res = sort(res)
plot(res)


for




