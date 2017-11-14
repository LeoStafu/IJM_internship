#function plotting variability explained by clustering depending on how much clusters are build
#point in red =  adding one more cluster will add less than 5% of explained variability
#IN : any data

#Arguments :
#nbPts = how many numbers of groups should be tested
#line = add a line beetween points
#elbowStep = color in red optimum number of clusters <=> adding a groups will add less than elbowStep %  variability

plotElbowKmeans = function(data, nbPts = 10, line = T, elbowStep = 5){
	
	res = c()

	for(i in 1:nbPts){
		#kmeans with 100 start so as to avoid multiple converging states 
		x = kmeans(data, i, nstart = 100)
		#calcul percent of variability explained by clusters
		score = 100 * (x$betweenss / x$totss)
		res = c(res, score)
	}

	plot(res, pch = 16, cex = 1.2, yaxt="n", xaxt="n", ylab = "between_SS / total_SS", xlab = "Number of cluster")
	
	elbow = which(diff(res) < elbowStep)[1]
	
	points(y = res[elbow], x = elbow, col = "red", pch = 16, cex = 1.2)
	axis(side = 1, at=seq(1,nbPts, by=1), labels = seq(1,nbPts,by = 1))
	axis(side = 2, at=seq(0,100, by=10), labels = paste(seq(0,100,by = 10), "%"), las = 2)
	if(line)
		lines(res)

}
