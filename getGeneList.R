
#modify condition of this list for different gene lists
#exemple 6 early : mcortMat[which(mcortMat$Conservation == 6 & mcortMat$periodClust == "early"),"protCodingGene"]
listGene = mcortMat[which(mcortMat$Conservation == 6),"protCodingGene"]

listGene = paste(listGene, collapse = "/")
listGene = strsplit(listGene, split = "/")
listGene = unlist(listGene)
listGene = unique(listGene)
listGene = paste(listGene, collapse="\n")
write.table(listGene, file="geneInteret.txt")
