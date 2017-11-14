

#Get ontologies on ensembl BiomaRt
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
martGO = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

GOswag = c()

for(i in 1:dim(swagMatrix4)[1]){

	cat(i, "\n")

	repeat{
		x <- try(getBM(attributes = c('hgnc_symbol','chromosome_name','start_position','end_position','gene_biotype','status', "go_id", "name_1006","definition_1006", "namespace_1003"), 
	                    filters = c('chromosome_name', 'start', 'end','status', 'biotype'), 
	                    values = list(substr(swagMatrix4[i,"Chromosome"], start = 4, stop = 5), swagMatrix4[i,"Begin"], swagMatrix4[i,"End"], "KNOWN", "protein_coding"), mart = ensembl)
				)
		if(!("try-error" %in% class(x))) break
	}

	x = x[which(x$namespace_1003 == "biological_process"),]

	if(dim(x)[1] != 0){
		GOswag = rbind(GOswag, cbind(x$name_1006, paste("seg", i), swagMatrix4[i,"Conservation"], swagMatrix4[i,"posToCentro"], swagMatrix4[i,"periodClust"]))
	}

}
colnames(GOswag) = c("OntologyBP", "Segment", "Conservation", "PosToCentro", "periodClust")

#ontologies counted only one time by segment even if present multiple times
GOswagU = unique(GOswag)

##################################################################################################################################################

#Study of ontologies

#Define groups
goCons6 = sort(table(GOswagU[which(GOswagU$Conservation == 6), "OntologyBP" ]))
goCons6 = 100 * (goCons6 / sum(goCons6))

goCons2 = sort(table(GOswagU[which(GOswagU$Conservation == 2), "OntologyBP" ]))
goCons2 = 100 * (goCons2 / sum(goCons2))

goCons6E = sort(table(GOswagU[which(GOswagU$periodClust == "early" & GOswagU$Conservation == 6), "OntologyBP" ]))
goCons6E = 100 * (goCons6E / sum(goCons6E))

goCons6L = sort(table(GOswagU[which(GOswagU$periodClust == "late" & GOswagU$Conservation == 6), "OntologyBP" ]))
goCons6L = 100 * (goCons6L / sum(goCons6L))

goCons2L = sort(table(GOswagU[which(GOswagU$periodClust == "late" & GOswagU$Conservation == 2), "OntologyBP" ]))
goCons2L = 100 * (goCons2L / sum(goCons2L))

goCons2E = sort(table(GOswagU[which(GOswagU$periodClust == "early" & GOswagU$Conservation == 2), "OntologyBP" ]))
goCons2E = 100 * (goCons2E / sum(goCons2E))

goConsE = sort(table(GOswagU[which(GOswagU$periodClust == "early"), "OntologyBP" ]))
goConsE = 100 * (goConsE / sum(goConsE))

goConsL = sort(table(GOswagU[which(GOswagU$periodClust == "late"), "OntologyBP" ]))
goConsL = 100 * (goConsL / sum(goConsL))

goAll = sort(table(GOswagU$OntologyBP))
goAll = 100 * (goAll / sum(goAll))


#NOT UNIQUE IN SEGMENTS
if(TRUE){

	goCons6 = sort(table(GOswag[which(GOswag$Conservation == 6), "OntologyBP" ]))
	goCons6 = 100 * (goCons6 / sum(goCons6))

	goCons2 = sort(table(GOswag[which(GOswag$Conservation == 2), "OntologyBP" ]))
	goCons2 = 100 * (goCons2 / sum(goCons2))

	goCons6E = sort(table(GOswag[which(GOswag$periodClust == "early" & GOswag$Conservation == 6), "OntologyBP" ]))
	goCons6E = 100 * (goCons6E / sum(goCons6E))

	goCons6L = sort(table(GOswag[which(GOswag$periodClust == "late" & GOswag$Conservation == 6), "OntologyBP" ]))
	goCons6L = 100 * (goCons6L / sum(goCons6L))

	goCons2L = sort(table(GOswag[which(GOswag$periodClust == "late" & GOswag$Conservation == 2), "OntologyBP" ]))
	goCons2L = 100 * (goCons2L / sum(goCons2L))

	goCons2E = sort(table(GOswag[which(GOswag$periodClust == "early" & GOswag$Conservation == 2), "OntologyBP" ]))
	goCons2E = 100 * (goCons2E / sum(goCons2E))

	goConsE = sort(table(GOswag[which(GOswag$periodClust == "early"), "OntologyBP" ]))
	goConsE = 100 * (goConsE / sum(goConsE))

	goConsL = sort(table(GOswag[which(GOswag$periodClust == "late"), "OntologyBP" ]))
	goConsL = 100 * (goConsL / sum(goConsL))

	goAll = sort(table(GOswag$OntologyBP))
	goAll = 100 * (goAll / sum(goAll))

}


#Compute the difference in proportions of ontologies in a given group "goTable" (groups defined above) and the group goAll (all the segments)
diffPropOnto = function(goTable){
	res = data.frame(matrix(nrow = length(goAll), ncol = 2))
	colnames(res) = c("OntologyBP", "diff")
	for(i in 1:length(goAll)){

		w = goTable[names(goAll[i])]

		if(is.na(w)){w=0}
			res[i, "diff"] = w - goAll[i]
			res[i, "OntologyBP"] = names(goAll[i])
	}	
	#res  =res[complete.cases(res),]
	#keep the 0.5% with the biggest difference
	res = res[which(res$diff >= quantile(res$diff, 0.995) | res$diff <= quantile(res$diff, 0.005)),]
	rownames(res) = res$OntologyBP
	return(res)
}

#Get all ontologies and their difference
res6E = diffPropOnto(goCons6E)
res6L = diffPropOnto(goCons6L)
res2E = diffPropOnto(goCons2E)
res2L = diffPropOnto(goCons2L)
res6 = diffPropOnto(goCons6)
res2 = diffPropOnto(goCons2)
resE = diffPropOnto(goConsE)
resL = diffPropOnto(goConsL)

#concatenate all ontologies, one that is in the 1% most different for one group might be absent for another group but still interresting to study
#old allOnto avec les 1% max => allOnto1pc
allOnto1pc = c(res6E$OntologyBP, res6L$OntologyBP, res2E$OntologyBP, res2L$OntologyBP, res6$OntologyBP, res2$OntologyBP, resE$OntologyBP, resL$OntologyBP)
allOnto1pc = unique(allOnto1pc)

allOnto1pcEL = c(resE$OntologyBP, resL$OntologyBP)
allOnto1pcEL = unique(allOnto1pcEL)

allOnto1pc62 = c(res6$OntologyBP, res2$OntologyBP)
allOnto1pc62 = unique(allOnto1pc62)

allOnto1pc6EL = c(res6E$OntologyBP, res6L$OntologyBP)
allOnto1pc6EL = unique(allOnto1pc6EL)

allOnto1pc2EL = c(res2E$OntologyBP, res2L$OntologyBP)
allOnto1pc2EL = unique(allOnto1pc2EL)

allOnto1pc2L6L = c(res2L$OntologyBP, res6L$OntologyBP)
allOnto1pc2L6L = unique(allOnto1pc2L6L)

allOnto1pc2E6E = c(res2E$OntologyBP, res6E$OntologyBP)
allOnto1pc2E6E = unique(allOnto1pc2E6E)



#re get all ontologies proportions differences based on the concatenated list
diffPropOnto = function(goTable){
	res = data.frame(matrix(nrow = length(goAll), ncol = 2))
	colnames(res) = c("OntologyBP", "diff")
	for(i in 1:length(goAll)){
		w = goTable[names(goAll[i])]
		if(is.na(w)){w=0}
			res[i, "diff"] = w - goAll[i]
			res[i, "OntologyBP"] = names(goAll[i])
	}	
	rownames(res) = res$OntologyBP
	return(res)
}

res6E = diffPropOnto(goCons6E)
res6L = diffPropOnto(goCons6L)
res2E = diffPropOnto(goCons2E)
res2L = diffPropOnto(goCons2L)
res6 = diffPropOnto(goCons6)
res2 = diffPropOnto(goCons2)
resE = diffPropOnto(goConsE)
resL = diffPropOnto(goConsL)


#########################################################
#Plots of the heatmap after scaling

lol = data.frame(matrix(nrow = length(allOnto1pc), ncol = 0))
rownames(lol) = allOnto1pc
lol[,"6E"] = res6E[allOnto1pc,"diff"]
lol[,"6L"] = res6L[allOnto1pc,"diff"]
lol[,"2E"] = res2E[allOnto1pc,"diff"]
lol[,"2L"] = res2L[allOnto1pc,"diff"]
lol[,"6"] = res6[allOnto1pc,"diff"]
lol[,"2"] = res2[allOnto1pc,"diff"]
lol[,"E"] = resE[allOnto1pc,"diff"]
lol[,"L"] = resL[allOnto1pc,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

library(gplots)
#heatmap a la main
pdf("allHM.pdf", height = 20, width = 20)
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), density.info="none", trace="none", cexRow = 0.4, col = colorsHM, margins = c(2,10), keysize = 1, cexCol = 0.8)
dev.off()
 
################################### ADD JC DERNIERE MINUTE

lol = data.frame(matrix(nrow = length(allOnto1pc2L6L), ncol = 0))
rownames(lol) = allOnto1pc2L6L
lol[,"2L"] = res2L[allOnto1pc2L6L,"diff"]
lol[,"6L"] = res6L[allOnto1pc2L6L,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("late62_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1)
dev.off()

###################################

lol = data.frame(matrix(nrow = length(allOnto1pc2E6E), ncol = 0))
rownames(lol) = allOnto1pc2E6E
lol[,"2E"] = res2E[allOnto1pc2E6E,"diff"]
lol[,"6E"] = res6E[allOnto1pc2E6E,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("early62_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1)
dev.off()


###################################

lol = data.frame(matrix(nrow = length(allOnto1pcEL), ncol = 0))
rownames(lol) = allOnto1pcEL
lol[,"E"] = resE[allOnto1pcEL,"diff"]
lol[,"L"] = resL[allOnto1pcEL,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("EL_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1)
dev.off()




lol = data.frame(matrix(nrow = length(allOnto1pc62), ncol = 0))
rownames(lol) = allOnto1pc62
lol[,"6"] = res6[allOnto1pc62,"diff"]
lol[,"2"] = res2[allOnto1pc62,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("62_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1,)
dev.off()



lol = data.frame(matrix(nrow = length(allOnto1pc6EL), ncol = 0))
rownames(lol) = allOnto1pc6EL
lol[,"6E"] = res6E[allOnto1pc6EL,"diff"]
lol[,"6L"] = res6L[allOnto1pc6EL,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("6EL_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), dendrogram="row", density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1)
dev.off()


lol = data.frame(matrix(nrow = length(allOnto1pc2EL), ncol = 0))
rownames(lol) = allOnto1pc2EL
lol[,"2E"] = res2E[allOnto1pc2EL,"diff"]
lol[,"2L"] = res2L[allOnto1pc2EL,"diff"]

lolScale = apply(lol,2, scale)
rownames(lolScale) = rownames(lol)

pdf("2EL_HM.pdf")
colorsHM = colorRampPalette(c("blue4","blue1","white","red", "red4"))(20)
heatmap.2(as.matrix(lolScale), dendrogram="row", density.info="none", trace="none", cexRow = 0.2, col = colorsHM, margins = c(7,10), keysize = 1)
dev.off()








