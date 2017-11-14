



save(swagMatrix, swagMatrix2, swagMatrix3, swagMatrix4, file="swagMatrix.RData")





library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

x <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','strand','band','percentage_gene_gc_content','gene_biotype','status'), 
                    filters = c('chromosome_name', 'start', 'end','status'), 
                    values = list("1","564424", "1109858","KNOWN"), mart = ensembl)




##########################################################################################################################################

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

posSnp = getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
               filters = c('chr_name','start','end'), 
               values = list(1,564424,1109858), 
               mart = snpmart)

res = matrix(0, ncol = 5, nrow = dim(swagMatrix3)[1])
colnames(res) = c("nbMut", "nbInsert", "nbDel", "nbTransition", "nbTransversion")
for(i in 1:dim(swagMatrix3)[1]){

	repeat {
		x = try(getBM(attributes = c('refsnp_id','allele'), 
	               filters = c('chr_name','start','end'), 
	               values = list(substr(swagMatrix3[i,"Chromosome"], start = 4, stop = 5), swagMatrix3[i,"Begin"], swagMatrix3[i,"End"]), 
	               mart = snpmart)
			)
		if(!("try-error" %in% class(x))) break
	}
	res[i,"nbMut"] = dim(x)[1]

	for(j in 1:dim(x)[1]){

		w = unlist(strsplit(x[j,"allele"], split="/"))
		
		if(is.na(w[1]) || is.na(w[2]))
			next
		if(w[1] == "-"){
			res[i,"nbInsert"] = res[i,"nbInsert"] + 1
		}else if(w[2] == "-"){
			res[i,"nbDel"] = res[i,"nbDel"] + 1
		}else if( (w[1] == "A" && w[2] == "G") ||(w[1] == "C" && w[2] == "T")){
			res[i,"nbTransition"] = res[i,"nbTransition"] + 1
		}else{
			res[i,"nbTransversion"] = res[i,"nbTransversion"] + 1
		}

	}
	cat(i,"\n")
}	
save(res, file="mutCountAll.RData")




##########################################################################################################################################
res = c()
for(i in 1:dim(swagMatrix3)[1]){

	repeat {
		x = try(getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','strand','gene_biotype','status'), 
                    filters = c('chromosome_name', 'start', 'end','status'), 
                    values = list(substr(swagMatrix3[i,"Chromosome"], start = 4, stop = 5), swagMatrix3[i,"Begin"], swagMatrix3[i,"End"], "KNOWN"), mart = ensembl)
			)
		if(!("try-error" %in% class(x))) break
	}
	
	res  =c (res, dim(x)[1])
	
	cat(i, "/11505\n")

}

#count non coding RNA

res = c()
for(i in 1:dim(swagMatrix3)[1]){

	repeat {
		x = try(getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','gene_biotype','status'), 
                    filters = c('chromosome_name', 'start', 'end','status'), 
                    values = list(substr(swagMatrix3[i,"Chromosome"], start = 4, stop = 5), swagMatrix3[i,"Begin"], swagMatrix3[i,"End"], "KNOWN"), mart = ensembl)
			)
		if(!("try-error" %in% class(x))) break
	}
	
	res  =c (res, length(grep("RNA", x$gene_biotype)))
	
	cat(i, "/11505\n")

}





##########################################################################################################################################

res = c()
for(i in 1:dim(swagMatrix4)[1]){

	repeat {
		x = try(getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','strand','gene_biotype','status'), 
                    filters = c('chromosome_name', 'start', 'end','status', 'biotype'), 
                    values = list(substr(swagMatrix4[i,"Chromosome"], start = 4, stop = 5), swagMatrix4[i,"Begin"], swagMatrix4[i,"End"], "KNOWN", "protein_coding"), mart = ensembl)
			)
		if(!("try-error" %in% class(x))) break
	}

	x = x[which(x$hgnc_symbol != ""),"hgnc_symbol"]

	res  =c (res, paste(x, collapse = "/"))
	
	cat(i, "\n")

}	






##########################################################################################################################################

#GC percent

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

res = c()

for(i in 1:dim(swagMatrix4)[1]){

	x = getSeq(BSgenome.Hsapiens.UCSC.hg3, swagMatrix4[i,"Chromosome"], swagMatrix4[i,"Begin"], swagMatrix4[i,"End"])
	cat((i/11505)*100,"\n")
	res = c(res, letterFrequency(x, "GC", as.prob = TRUE))

}


##########################################################################################################################################

#Def period

res = c()
for(i in 1:dim(swagMatrix3)[1]){

	if(swagMatrix3[i,"mean"] < quantile(swagMatrix3$mean, 0.3333)){
		res = c(res, "late")
	}else if(swagMatrix3[i,"mean"] < quantile(swagMatrix3$mean, 0.6666)){
		res = c(res, "mid")
	}else{
		res = c(res,"early") 
	}
}	


res = c()
for(i in 1:dim(swagMatrix4)[1]){

	if(swagMatrix4[i,"freqLINE"] < quantile(swagMatrix4$freqLINE, 0.20)){
		res = c(res, "Low")
	}else if(swagMatrix4[i,"freqLINE"] > quantile(swagMatrix4$freqLINE, 0.20) & swagMatrix4[i,"freqLINE"] < quantile(swagMatrix4$freqLINE, 0.40 )) {
		res = c(res, "midLow")
	}else if(swagMatrix4[i,"freqLINE"] > quantile(swagMatrix4$freqLINE, 0.40) & swagMatrix4[i,"freqLINE"] < quantile(swagMatrix4$freqLINE, 0.60 )) {
		res = c(res, "mid")
	}else if(swagMatrix4[i,"freqLINE"] > quantile(swagMatrix4$freqLINE, 0.60) & swagMatrix4[i,"freqLINE"] < quantile(swagMatrix4$freqLINE, 0.80 )) {
		res = c(res, "midHigh")
	}else{
		res = c(res,"high") 
	}
}	



##########################################################################################################################################
##############################################################################################################################
#BArplot period by bp

res = matrix(ncol = 5, nrow = 3)
rownames(res) = c("early", "late", "mid")
colnames(res) = c("2","3","4","5","6")

for(z in 2:6){
	for(i in c("early", "late", "mid") )
		res[i,z-1] = sum(swagMatrix4[which(swagMatrix4$Conservation == z & swagMatrix4$periodClust == i), "lengthSeg"])
}


barplot(prop.table(res, margin = 2), 
	col=c("red","darkgrey", "orange"), ylim = c(0,1.2), 
	main="Proportion of periods of replication timing", 
	sub="Proportion based on number of bp\ncell Lines : HEK293T HeLa K562 MRC5C RKO U2OS")
legend("top",title="RT period", c("early","mid","late"), bty="n",fill=c("red","darkgrey", "orange"), horiz = T)
##############################################################################################################################
##############################################################################################################################
##########################################################################################################################################

#swagMatrix3 = swagMatrix2 but no overlapping seg, cut at the middle


swagMatrix3 = swagMatrix2[,1:12]

for(i in 1:(dim(swagMatrix3)[1]-1)){

	cat(i,"\n")
	if((swagMatrix3[i,"Chromosome"] == swagMatrix3[(i+1), "Chromosome"]) && (swagMatrix3[i,"posToCentro"] == swagMatrix3[(i+1), "posToCentro"])){

		x = (swagMatrix3[i,"End"] + swagMatrix3[(i+1),"Begin"]) / 2
		swagMatrix3[i,"End"] = x
		swagMatrix3[(i+1),"Begin"] = x
	}

}

##########################################################################################################################################
#get mean and sd for swagMatrix3


res = matrix(ncol = 2, nrow = 11505)
colnames(res) =c("mean", "sd")

for(i in 1:dim(swagMatrix3)[1]){

	RTloc = RTloess[which(RTloess$Chromosome == swagMatrix3[i,"Chromosome"]   & RTloess$posToCentro == swagMatrix3[i,"posToCentro"]),]

	j = min(which(RTloc$Position >= swagMatrix3[i,"Begin"]))
	k = min(which(RTloc$Position >= swagMatrix3[i,"End"])) - 1

	x = unlist(strsplit(as.character(swagMatrix3[i,"whichConserv"]), split="-"))

	res[i,"mean"] = mean(as.matrix(RTloc[j:k,x]))
	res[i,"sd"] = sd(as.matrix(RTloc[j:k,x]))

	cat("Seg",i, "From", j, "to", k, "\n")

}



#########################################################################################################################


hsaps_bindMotif = useMart(biomart = "ENSEMBL_MART_FUNCGEN", dataset = "hsapiens_motif_feature")


res = matrix(0, ncol = 2, nrow = dim(swagMatrix4)[1])
colnames(res) = c("nbBinding", "binding")
for(i in 1:dim(swagMatrix4)[1]){

	repeat {
		x = try(getBM(attributes = listAttributes(hsaps_bindMotif)$name, 
          				filters = c('chromosome_name', 'start', 'end'), 
          				values = list(substr(swagMatrix4[i,"Chromosome"], start = 4, stop = 5), swagMatrix4[i,"Begin"], swagMatrix4[i,"End"]),
          				mart = hsaps_bindMotif)
			)
		if(!("try-error" %in% class(x))) break
	}


	res[i,"nbBinding"] = length(unique(x$feature_type_name))
	res[i,"binding"]  =  paste(unique(x[which(x$feature_type_name != ""),"feature_type_name"]), collapse = "/")
	
	cat(i,"\n")
}	
save(res, file="bindDataSeg.RData")



#########################################################################################################################


hsaps_strVar = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_structvar")

res = matrix(0, ncol = 4, nrow = dim(swagMatrix4)[1])
colnames(res) = c("CNV", "insert","mobileEltInsert","shortTandemRepeatVar")
for(i in 1:dim(swagMatrix4)[1]){

	repeat {
			x = try(getBM(attributes = c("chr_name", "chrom_start","chrom_end","sv_variant_type"), 
	          				filters = c('chr_name', 'start', 'end'), 
	          				values = list(substr(swagMatrix4[i,"Chromosome"], start = 4, stop = 5), swagMatrix4[i,"Begin"], swagMatrix4[i,"End"]),
	          				mart = hsaps_strVar)
				)
			if(!("try-error" %in% class(x))) break
		}

	if(is.na(table(x$sv_variant_type)["copy_number_variation"])){
		res[i,"CNV"] = table(x$sv_variant_type)["copy_number_variation"]
	}

	if(is.na(table(x$sv_variant_type)["insertion"])){
		res[i,"insert"] = table(x$sv_variant_type)["insertion"]
	}

	if(is.na(table(x$sv_variant_type)["mobile_element_insertion"])){
		res[i,"mobileEltInsert"] = table(x$sv_variant_type)["mobile_element_insertion"]
	}

	if(is.na(table(x$sv_variant_type)["short_tandem_repeat_variation"])){
		res[i,"shortTandemRepeatVar"] = table(x$sv_variant_type)["short_tandem_repeat_variation"]
	}




	cat(i,"\n")
}	

#########################################################################################################################

#COUNT SNP HOME

chr = "chr1"

countMut = function(chr, strand = 1){

	miniSwag = swagMatrix4[which(swagMatrix4$Chromosome == chr & swagMatrix4$posToCentro == strand),c("Chromosome", "Begin", "End")]

	repeat {
		db = try( getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_end'), 
                filters = c('chr_name', "start", "end"), 
                values = list(substr(chr, start = 4, stop = 5), min(miniSwag$Begin), max(miniSwag$End)), 
                mart = snpmart)
				)

    	if(!("try-error" %in% class(db))) break
	}


	res = matrix(0, ncol = 4, nrow = dim(miniSwag)[1])
	colnames(res) = c("nbInsert", "nbDel", "nbTransition","nbTransversion")


	for(i in 1:dim(miniSwag)[1]){

		cat(i,"\n")

		lol = which((db$chrom_start >= miniSwag[i,"Begin"]) & (db$chrom_end <= miniSwag[i,"End"]))

		x = db[ lol, ]

		db = db[-lol, ]

		for(j in 1:dim(x)[1]){

			w = unlist(strsplit(x[j,"allele"], split="/"))
			if(is.null(w))
				next
			if(is.na(w[1]) || is.na(w[2]))
				next
			if(w[1] == "-"){
				res[i,"nbInsert"] = res[i,"nbInsert"] + 1
			}else if(w[2] == "-"){
				res[i,"nbDel"] = res[i,"nbDel"] + 1
			}else if( (w[1] == "A" && w[2] == "G") ||(w[1] == "C" && w[2] == "T")){
				res[i,"nbTransition"] = res[i,"nbTransition"] + 1
			}else{
				res[i,"nbTransversion"] = res[i,"nbTransversion"] + 1
			}

		}

	}	

	return(res)
}


#########################################################################################################################



res = matrix(0, ncol = 8, nrow = dim(swagMatrix4)[1])
colnames(res) = c("nbLINE", "nbSINE","nbLTR","nbSimpleRepeat", "nbLowCmplx", "nbDNArep","nbSatell", "nbRetrop")
for(i in 1:dim(swagMatrix4)[1]){

	x = repp[which(as.character(repp$genoName) == as.character(swagMatrix4[i,"Chromosome"]) & repp$genoStart >= swagMatrix4[i,"Begin"] & repp$genoEnd <= swagMatrix4[i,"End"]), ]


	if(!is.na(table(x$repClass)["LINE"])){
		res[i,"nbLINE"] = table(x$repClass)["LINE"]
	}

	if(!is.na(table(x$repClass)["SINE"])){
		res[i,"nbSINE"] = table(x$repClass)["SINE"]
	}

	if(!is.na(table(x$repClass)["LTR"])){
		res[i,"nbLTR"] = table(x$repClass)["LTR"]
	}

	if(!is.na(table(x$repClass)["Simple_repeat"])){
		res[i,"nbSimpleRepeat"] = table(x$repClass)["Simple_repeat"]
	}

	if(!is.na(table(x$repClass)["Low_complexity"])){
		res[i,"nbLowCmplx"] = table(x$repClass)["Low_complexity"]
	}

	if(!is.na(table(x$repClass)["DNA"])){
		res[i,"nbDNArep"] = table(x$repClass)["DNA"]
	}

	if(!is.na(table(x$repClass)["Satellite"])){
		res[i,"nbSatell"] = table(x$repClass)["Satellite"]
	}

	if(!is.na(table(x$repClass)["Retroposon"])){
		res[i,"nbRetrop"] = table(x$repClass)["Retroposon"]
	}



	cat(i,"\n")
}	