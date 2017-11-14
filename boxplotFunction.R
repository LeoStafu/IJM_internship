#function plotting boxplot of a given descriptor of the MCORT matrix
#boxplot grouped by Conservation of RT and by RT class (early, mid, late)
#generate random samples of early, mid and late segment
#IN : name of the column of the descriptor, optionnal graphical parameter could be added


boxplotSwag = function(colname, mat = mcortMat, label = as.character(colname), pVal = F, ...){

	#define margins	
	par(mar = c(6,4,3,2) +0.1)
	#Define x positions where boxplots will stand, make place for random sample + margin beetween groups
	seqX = c(1,2,3,4,5, 8,9,10,11,12, 15,16,17,18,19)

	#Boxplot plotted and its data saved in the object b (countaining size, means...)
	b = boxplot(mat[,colname] ~ mat$Conservation + mat$periodClust,
		at = seqX, xlim = c(1,20),
		col = rainbow(5), main=paste(colname,"grouped by Conservation and Replication Timing"), ylab = label, cex.axis = 0.8, lwd = 1.4 ,  ...)
	axis(label = paste("n=",b$n), at = seqX, side = 1, line = 1, cex.axis = 0.8, lty = 0, ...)

	#Random sample for each group (early, mid, late). Size = mean size of the group
	boxplot(mat[sample(which(mat$periodClust == "early"), mean(b$n[1:5])), colname], at=6, add = T, yaxt = "n", col ="grey70" )
	boxplot(mat[sample(which(mat$periodClust == "mid"), mean(b$n[6:10])), colname], at=13, add = T, yaxt = "n", col ="grey70" )
	boxplot(mat[sample(which(mat$periodClust == "late"), mean(b$n[11:15])), colname], at=20, add = T, yaxt = "n", col ="grey70" )

	axis(label = c(paste("Random\nearly\nn=", mean(b$n[1:5])), paste("Random\nmid\nn=", mean(b$n[6:10])), paste("Random\nlate\nn=", mean(b$n[11:15]))),
		at = c(6,13,20), side = 1, line = 1, cex.axis = 0.8, lty = 0)

	#Generate empirical pvalues (parameter of the function, by default = FALSE)
	if(pVal){
		pValues = c()
		for(w in 1:15){

			means = c()
			#generate 10 000 samples of the size of a given sample
			for(i in 1:10000){
				means = c(means, mean(mat[sample(nrow(mat), b$n[w]), colname])	)
			}

			#pcal 5% = * or 1% = **
			if(is.na(table(b$stats[3,w] < means)["TRUE"]))
				pValues = c(pValues, "**")
			else if ((table(b$stats[3,w] < means)["TRUE"] /10000 ) < 0.025 | ((table(b$stats[3,w] < means)["TRUE"] /10000 ) > 0.975))
				pValues = c(pValues, "*")
			else if ((table(b$stats[3,w] < means)["TRUE"] /10000 ) < 0.005 | ((table(b$stats[3,w] < means)["TRUE"] /10000 ) > 0.995))
				pValues = c(pValues, "**")
			else
				pValues = c(pValues, " ")		
		}
		axis(label =pValues, at = seqX, side = 1, line = 2, cex.axis = 0.8, lty = 0, ...)
	}
}