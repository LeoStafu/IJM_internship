#######################################################################################################################
#Matrix construction for replication timing conservation analysis

#Aim : contruct a matrix for replication timing conservation analysis at different genomic regions
#######################################################################################################################



#get data, already normalized
RTall = read.table("rawData.txt", header=T)

#define list of cell Lines and list of chromosomes names
cellLines = c("HEK293T", "HeLa", "K562", "MRC5C", "RKO", "U2OS")
allChr  = c(paste("chr",c(1:22), sep=""),"chrX")


#Transform a data frame of raw Replication Timing data into the same data frame with loess Data
rawRT_to_loessRT = function(RTdf, loessSpan = 1200000, cellLines = c("HEK293T", "HeLa", "K562", "MRC5C", "RKO", "U2OS")){

	for(cell in cellLines){
		#Here Loess computed on the whole chromosome
		loesschr = loess(RTdf[,cell] ~ RTdf$Position, span = loessSpan/(max(RTdf$Position)-min(RTdf$Position)))
		RTdf[,cell] = loesschr$fitted
	}
	return(RTdf)
}



#function generating subMatrix
#window size in number of points <=> number of probes, avg. 20kb per probe
getSubMatrix = function(RT, chr, arm, wSize = 30 ,cellLines = c("HEK293T", "HeLa", "K562", "MRC5C", "RKO", "U2OS") ){

	subMatrix = data.frame(NULL,matrix(ncol = 14, nrow = 0))
	colnames(subMatrix) = c("Chromosome", "posToCentro", "Begin", "End", cellLines, "Conservation", "whichConserv", "mean", "sd")
	
	#Get RT for the given chromosome
	RTchr = RT[which(RT$Chromosome == chr), ]
	#get the proper arm
	RTchr = RTchr[which(RTchr$posToCentro == arm), ]

	#Transform raw data of replication timing with loess
	RTchr = rawRT_to_loessRT(RTchr, loessSpan = 1200000 )
	
	#For moving window
	range1 = 1
	range2 = wSize

	colCell = rainbow(length(cellLines))

	#analysis, sweep all position in the chromosome arm
	#generate plots in a pdf file
	pdf(paste("RT_",chr,"arm",arm,".pdf", sep=))

	while(range2 < dim(RTchr)[1]){

		#Plot the loess curve of all cell Lines for the active window
		nbCell = 1
		plot(NULL,xlim=c(RTchr[range1,"Position"],RTchr[range2,"Position"]), ylim=c(-3,3), 
			ylab="Timing (Z-score)", xlab="Position on chromosome (bp)",
			main=paste("Replication timing on",RTchr[1,"Chromosome"]), 
			sub=paste("Position (arb. unit)",range1,"to",range2))
		for(cell in cellLines){
			lines(RTchr[,cell]~RTchr$Position, col=colCell[nbCell], lwd = 3	)
			nbCell = nbCell +1
		}
		abline(h=0, col = "black", lty = "dotted")

		#get RT values of each cell lines for the given position
		dataSeg = t(RTchr[range1:range2, cellLines])
		#Build a hierarchical clustering tree based on euclidian distance
		#Method complete <=> dist beetween cluster = max distance
		tree = hclust(dist(dataSeg, method="euclidian"), method="complete")
		
		#Cutoff found by studying the elbow af the curve of all height of all tree (all position).
		cutoff = 4.043117

		#Is there a curve diffenrent from another based on the distance matrix ?
		if( max(tree$height) > cutoff ){
			
			cellGroups = try(rect.hclust(tree, 	h = cutoff))
			if(!("try-error" %in% class(cellGroups))){

				sizeGrp = lapply(cellGroups, length)

				pbGrp = cellGroups[which(sizeGrp == max(unlist(sizeGrp)))]

				#Which groups of similar curves have the same size ?
				if(length(pbGrp) > 1){
					#If more than one group have the max number of cell lines in it, find the group whose curves are the most similar
					minDist = max(dist(dataSeg[names(pbGrp[[1]]), ]))
					grp = 1
					for(i in 2:length(pbGrp)){
						dist  = max(dist(dataSeg[names(pbGrp[[i]]), ]))
					 	if(minDist > dist ){
					 		minDist = dist
					 		grp = i
				 		}
					}
					segInfo = data.frame(chr, arm, RTchr[range1,"Position"], RTchr[range2,"Position"], t(as.matrix(cutree(tree, h = cutoff))), length(pbGrp[[grp]]), 
						paste(names(pbGrp[[grp]]), collapse="-") , mean(dataSeg[names(pbGrp[[grp]]), ]), sd(dataSeg[names(pbGrp[[grp]]), ]) )
				}else{
					#One and only group with the maximum number of cell Lines
					segInfo = data.frame(chr, arm, RTchr[range1,"Position"], RTchr[range2,"Position"], t(as.matrix(cutree(tree, h = cutoff))), length(pbGrp[[1]]), 
						paste(names(pbGrp[[1]]), collapse="-"), mean(dataSeg[names(pbGrp[[1]]), ]), sd(dataSeg[names(pbGrp[[1]]), ])  )
				}	
				#Legend on plot, wich cell lines are conserved
				legend("topright",inset = 0.02, legend =names(pbGrp[[1]]), bty="n", lty=1, lwd = 3, col = colCell[cellLines %in% names(pbGrp[[1]])])
			}else{
				#error, 0 groups, conserved into 0 cell lines
				segInfo = data.frame(chr, arm, RTchr[range1,"Position"], RTchr[range2,"Position"], t(as.matrix(seq(1,6))), 0, "NULL", NA, NA)
			}

		}else{
			#Conserved into all cell Lines
			segInfo = data.frame(chr, arm, RTchr[range1,"Position"], RTchr[range2,"Position"], t(as.matrix(rep(1,6))), 6, 
				paste(cellLines, collapse="-"), mean(dataSeg), sd(dataSeg)  )
			
			legend("topright",inset = 0.02, legend =cellLines, bty="n", lty=1, lwd = 3, col=colCell)

		}

		colnames(segInfo) = c("Chromosome", "posToCentro", "Begin", "End", cellLines, "Conservation", "whichConserv", "mean", "sd")
		subMatrix = rbind(subMatrix, segInfo)
		range1 = range1 + wSize/2
		range2 = range2 + wSize/2

	}
	dev.off()
	
	return(subMatrix)
}
	


#Concatenate when the conserved portion of the genome span across multiple windows in the raw matrix
simplifyMatrix = function(RTmat){

	newMat = RTmat[,1:12]
	suppr = c()

	i = 1
	while(i < dim(RTmat)[1] ){

		if(all( newMat[i,5:12] == newMat[i+1, 5:12])){

			newMat[i+1,"Begin"] = newMat[i,"Begin"]
			suppr = c(suppr,i)	
		}
		i = i+1
	}

	if(!is.null(suppr))
		newMat = newMat[-suppr,]
	return(newMat)
}

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################









#PAS MAL DE VRAC A PARTIR DE LA

#################################################################     RUN     #############################################################

swagMatrix2 = data.frame(NULL,matrix(ncol = 14, nrow = 0))
colnames(swagMatrix2) = c("Chromosome", "posToCentro", "Begin", "End", cellLines, "Conservation", "whichConserv","mean", "sd")
	

for(chr in allChr){

	if(chr %in% c("chr13","chr14","chr15","chr22")){
		x =  getSubMatrix(RTall, chr, arm = 2, wSize = 30)
		#x = simplifyMatrix(x)
		swagMatrix2 = rbind(swagMatrix2, x)
	}else{
		x =  getSubMatrix(RTall, chr, arm = 1, wSize = 30)
		#x = simplifyMatrix(x)
		swagMatrix2 = rbind(swagMatrix2, x)
		x =  getSubMatrix(RTall, chr, arm = 2, wSize = 30)
		#x = simplifyMatrix(x)
		swagMatrix2 = rbind(swagMatrix2, x)
	}
}





##########################################################     BARPLOT     ##########################################################



#Ce que djihad appelle Godzilla


chrBarlpot = function(chr, ybottom, ytop){

	x = swagMatrix[which(swagMatrix$Chromosome == chr & swagMatrix$posToCentro == 1), ]
	if(dim(x)[1] !=0){
		for(i in 1:dim(x)[1]){

			if(i == 1){
				xl = x[i,"Begin"]
				xr = (x[i+1,"Begin"] +x[i,"Begin"]) / 2
			}else if(i == dim(x)[1]){
				xl = (x[i-1,"Begin"] +x[i,"Begin"]) / 2
				xr = x[i,"End"] 
			}else{
				xl = (x[i-1,"Begin"] +x[i,"Begin"]) / 2
				xr = (x[i+1,"Begin"] +x[i,"Begin"]) / 2
			}
			rect(xleft = xl, xright = xr, ybottom = ybottom, ytop = ytop, col = colConserv[x[i,"Conservation"]], lty = 0)
		}
	}

	rect(xleft = centromere(chr,build="hg18")[1], xright = centromere(chr,build="hg38")[2], ybottom = ybottom, ytop = ytop, col = "black", lty = 0)

	x = swagMatrix[which(swagMatrix$Chromosome ==  chr & swagMatrix$posToCentro == 2), ]
	for(i in 1:dim(x)[1]){

		if(i == 1){
			xl = x[i,"Begin"]
			xr = (x[i+1,"Begin"] +x[i,"Begin"]) / 2
		}else if(i == dim(x)[1]){
			xl = (x[i-1,"Begin"] +x[i,"Begin"]) / 2
			xr = x[i,"End"] 
		}else{
			xl = (x[i-1,"Begin"] +x[i,"Begin"]) / 2
			xr = (x[i+1,"Begin"] +x[i,"Begin"]) / 2
		}
		rect(xleft = xl, xright = xr, ybottom = ybottom, ytop = ytop, col = colConserv[x[i,"Conservation"]], lty = 0)
	}	
}

#RUN DU GODZILLA
library(SNPchip)
colConserv = c("grey","purple", "red","yellow", "green","blue")

svg("genomeRT.svg", width=80, height=15)
plot(NULL, ylim = c(0,33), xlim = c(0, max(swagMatrix2$End)), yaxt = "n", ylab = "", xlab="" )

yb = 0
yt = 1

for(chr in allChr){
	chrBarlpot(chr, yb, yt)
	yb = yb +1.5
	yt = yt +1.5
}
axis(2, at=seq(from=0.5, to = 33.5, by = 1.5), labels = allChr, las= 2)


legend("topright",inset = 0.02, title="Conserved into :", legend = seq(2,6), bty="n", fill = c("purple", "red","yellow", "green","blue"), cex = 4 )


dev.off()







############################################################################################################################################################
############################################################################################################################################################
#RTloess

RTloess = data.frame(NULL,matrix(ncol = 14, nrow = 0))
colnames(RTloess) = c("Chromosome", "posToCentro", "Begin", "End", cellLines)
	

for(chr in allChr){

	if(chr %in% c("chr13","chr14","chr15","chr22")){
		x2 = RTall[which(RTall$Chromosome == chr & RTall$posToCentro == 2), ]
		x2 =  rawRT_to_loessRT(x2)
		RTloess = rbind(RTloess, x2)
	}else{
		x1 = RTall[which(RTall$Chromosome == chr & RTall$posToCentro == 1), ]
		x1 =  rawRT_to_loessRT(x1)
		RTloess = rbind(RTloess, x1)

		x2 = RTall[which(RTall$Chromosome == chr & RTall$posToCentro == 2), ]
		x2 =  rawRT_to_loessRT(x2)
		RTloess = rbind(RTloess, x2)
	}
}


pdf("RTtrees.pdf")
	plot(hclust(dist(t(RTloess[,3:8]))), main = "Distance of loess RT all genome", sub ="Distance euclidian, hclust complete", xlab ="")
	for(chr in allChr)
		plot(hclust(dist(t(RTloess[which(RTloess$Chromosome == chr),3:8])), method ="ward.D1"), main = paste("Distance of loess RT", chr), sub ="Distance euclidian, hclust complete", xlab ="")
dev.off()





####################################################################################################################
#									VRAC DATA TO FIND THREHOLD OF DIST FOR CUTTREE FUNCTION
####################################################################################################################


allHeight =c()
for(chr in allChr){

	if(chr %in% c("chr13","chr14","chr15","chr22")){
		x = getSubMatrix(RTall, chr, arm = 2, wSize = 30)
		allHeight  =c(allHeight, x)
	}

	x = getSubMatrix(RTall, chr, arm = 1, wSize = 30)
	allHeight  =c(allHeight, x)
	x = getSubMatrix(RTall, chr, arm = 2, wSize = 30)
	allHeight  =c(allHeight, x)
}



# fake data
x <- c(1,41750)
y <- c(min(allHeight), max(allHeight))
myData <- data.frame(x, y)

# fitted curve (I used a simpler example)
result <- lm(y ~ x, data=myData)
p <- predict(result)

# line connecting endpoints of fitted curve
i1 <- which.min(x)
i2 <- which.max(x)
slope <- (p[i2] - p[i1]) / (x[i2] - x[i1])
int <- p[i1] - slope*x[i1]

# for every point on the predicted curve (xi, pi), the perpendicular line that goes through that point has
perpslope <- -1/slope
perpint <- p - perpslope*x

# the intersection of the perp line(s) with the connecting line is
xcross <- (int - perpint) / (perpslope - slope)
ycross <- slope*xcross + int

# the distance between the intersection and the point(s) is
dists <- sqrt((x - xcross)^2 + (y - ycross)^2)
# the index of the farthest point
elbowi <- which.max(dists)

