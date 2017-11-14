



setwd(dir = "/home/laboratoire/Téléchargements/dbSNP/")


for(i in 1:dim(swagMatrix4)[1]){

	name = paste(swagMatrix4[i,"Chromosome"],"SNP", sep="")

	if(!exists(name)){
		filename = paste(swagMatrix4[i,"Chromosome"],"ok.txt", sep="")
		assign(name, read.table(filename, sep=";"))
	}
}





allChr  = c(paste("chr",c(1:22), sep=""),"chrX")
resF  = matrix(0, ncol = 7, nrow = 0)
for(chr in allChr){

	cat(chr, "\n")

	filename = paste(chr,"ok.txt", sep="")
	db = read.table(filename, sep=";", stringsAsFactors = F)

	for(arm in 1:2){
		cat("arm", arm, "\n")

		miniSwag = swagMatrix4[which(swagMatrix4$Chromosome == chr & swagMatrix4$posToCentro == arm),c("Chromosome", "Begin", "End")]
		res = matrix(0, ncol = 7, nrow = dim(miniSwag)[1])
		colnames(res) = c("deletion", "indel", "insertion", "SNV", "substitution", "transition","transversion")
		if(dim(miniSwag)[1] !=0 ){
			for(i in 1:dim(miniSwag)[1]){
				cat(i, "\n")

				lol = which((db$V3 >= miniSwag[i,"Begin"]) & (db$V4 <= miniSwag[i,"End"]))
				x = db[ lol, ]
				db = db[-lol, ]


				if(!is.na(table(x$V2)["deletion"])) res[i, "deletion"] = table(x$V2)["deletion"]
				if(!is.na(table(x$V2)["indel"])) res[i, "indel"] = table(x$V2)["indel"]
				if(!is.na(table(x$V2)["insertion"])) res[i, "insertion"] = table(x$V2)["insertion"]
				if(!is.na(table(x$V2)["SNV"])) res[i, "SNV"] = table(x$V2)["SNV"]
				if(!is.na(table(x$V2)["substitution"])) res[i, "substitution"] = table(x$V2)["substitution"]

				if(dim(x)[1] != 0){
				#study details of SNV
					for(j in 1:dim(x)[1]){

						if(x[j,"V2"] != "SNV") next
						w =	strsplit(x[j,"V5"], split="=")
						w = unlist(w)
						w = w[c(FALSE, TRUE)]
						w = unlist(strsplit(w, split=","))
						if(length(w) == 1) next

						if( (w[1] == "A" && w[2] == "G") ||(w[1] == "C" && w[2] == "T") || (w[1] == "G" && w[2] == "A") ||(w[1] == "T" && w[2] == "C")){
							res[i,"transition"] = res[i,"transition"] + 1
						}else{
							res[i,"transversion"] = res[i,"transversion"] + 1
						}
					}
				}

			}
		}

		resF = rbind(resF, res)

	}#end study arm

}
