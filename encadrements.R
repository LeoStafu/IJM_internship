
#Script which generated all figures in the encadrement folder 

#IN : the name of a chr as "chr1", "chrX" etc...
#		matrix name 

chrBarlpot2 = function(chr, mat = mcortMat){
	
	library(SNPchip)

	plot(NULL, ylim = c(0,20), xlim = c(0, max(mcortMat[which(mcortMat$Chromosome == chr), "End"])*10 ), yaxt = "n", ylab = "", xlab="" )

	colConserv = c("purple","grey", "red","yellow", "green","blue")
	colGene = c("grey20","grey50","grey80")

	#Loop for the first chromosomic arm
	x = mcortMat[which(mcortMat$Chromosome == chr & mcortMat$posToCentro == 1), ]
	if(dim(x)[1] !=0){
		for(i in 1:dim(x)[1]){

			#Period and conservation PLOT
			xl = x[i,"Begin"] *10
			xr = x[i,"End"] *10
			if(x[i,"periodClust"] == "early"){
				ytop = 9
			}else if(x[i,"periodClust"] == "mid"){
				ytop = 6
			}else{
				ytop = 3
			}
			rect(xleft = xl, xright = xr, ybottom = 0, ytop = ytop, col = colConserv[x[i,"Conservation"]], lty = 1, lwd = 2)

			#Gene coverage and GC% PLOT
			if(x[i,"GCclust"] == "Low"){
				ytop = 12
			}else if(x[i,"GCclust"] == "midLow"){
				ytop = 14
			}else if(x[i,"GCclust"] == "midHigh"){
				ytop = 16
			}else{
				ytop = 18
			}
			rect(xleft = xl, xright = xr, ybottom = 10, ytop = ytop, col = colGene[x[i,"protGeneClust"]], lty = 1, lwd = 2)
			#points(x= (xl+xr)/2, y = ytop, pch = 16, col = "black", cex = 1 )
			#abline(h=10)
		}
	}

	#Loop for the second chromosomic arm
	x = mcortMat[which(mcortMat$Chromosome ==  chr & mcortMat$posToCentro == 2), ]
	for(i in 1:dim(x)[1]){

			xl = x[i,"Begin"] *10
			xr = x[i,"End"] *10
			if(x[i,"periodClust"] == "early"){
				ytop = 9
			}else if(x[i,"periodClust"] == "mid"){
				ytop = 6
			}else{
				ytop = 3
			}
			rect(xleft = xl, xright = xr, ybottom = 0, ytop = ytop, col = colConserv[x[i,"Conservation"]], lty = 1, lwd = 2)

			if(x[i,"GCclust"] == "Low"){
				ytop = 12
			}else if(x[i,"GCclust"] == "midLow"){
				ytop = 14
			}else if(x[i,"GCclust"] == "midHigh"){
				ytop = 16
			}else{
				ytop = 18
			}
			#points(x= (xl+xr)/2, y = ytop, pch = 16, col = "black", cex = 1 )
			#abline(h=10)
			rect(xleft = xl, xright = xr, ybottom = 10, ytop = ytop, col = colGene[x[i,"protGeneClust"]], lty = 1, lwd = 2)
	}	
	legend("bottomright",inset = 0.01, title="Conserved\ninto :", legend = seq(2,6), bty="n", fill = c("grey", "red","yellow", "green","blue"), cex = 4 )
	legend("topright",inset = 0.01, title="Gene coverage:", legend = c("High", "Mid", "Low"), bty="n", fill = c("grey20","grey50","grey80"), cex = 4 )
	mtext("Timing E/M/L\t\tGC percent", side = 2, line = -10, cex = 5)
}



#Run the aove function for all chromosome from 1 to 22 + X
allChr  = c(paste("chr",c(1:22), sep=""),"chrX")

for(chr in allChr){
	svg(paste(chr, ".svg", sep=""), width=150, height=15)

	chrBarlpot2(chr)
	dev.off()
}
