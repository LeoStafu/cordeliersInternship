#Main function
#Argument : two cells population, a MPC counter output Matrix and annotation data

twoCellsToKM = function (pop1, pop2, MCPcMatrix, annotData){
	
	#Test pop1 pop2 and annot Data
	#Error messages
	#cat(pop1, pop2, "corners\n", sep=" : ")
	#try(
	#	plotCentroidAndKM(pop1, pop2, MCPcMatrix, annotData,"corners"),
	#	silent = FALSE
	#);

	cat(pop1, pop2, "mean\n", sep=" : ")
	try(
		plotCentroidAndKM(pop1, pop2, MCPcMatrix, annotData,"mean"),
		silent = FALSE
	);

	#cat(pop1, pop2, "majority\n", sep=" : ")
	#try(
	#	plotCentroidAndKM(pop1, pop2, MCPcMatrix, annotData,"majority"),
	#	silent = FALSE
	#);
}

#===========================================================================================================================================

#Initials centroid are based on extrem values of pop1 and pop2
methodCorners = function(pop1, pop2, MCPcMatrix){
	
	initCenters = matrix(ncol=2, nrow=4)
	colnames(initCenters)=c(pop1, pop2);
	initCenters[1,] = c( min(MCPcMatrix[pop1,]), min(MCPcMatrix[pop2,]) )
	initCenters[2,] = c( min(MCPcMatrix[pop1,]), max(MCPcMatrix[pop2,]) )
	initCenters[3,] = c( max(MCPcMatrix[pop1,]), min(MCPcMatrix[pop2,]) )
	initCenters[4,] = c( max(MCPcMatrix[pop1,]), max(MCPcMatrix[pop2,]) )
	x = as.matrix(MCPcMatrix[c(pop1,pop2),])
	x = t(x)
	res = kmeans(x, 4, centers=initCenters, iter = 1000)
	#print(res$centers)
	return(res);
}

#===========================================================================================================================================

#Mean of 200 kmeans method with random initial centroids
methodMean = function(pop1, pop2, MCPcMatrix){
	
	i=0
	j=0
	sumCentroid=matrix(data=0, nrow=4,ncol=2)
	rownames(sumCentroid) = c("HH","HL","LH","LL")
	colnames(sumCentroid) = c(pop1, pop2)
	for(j in 1:200){
		x = as.matrix(MCPcMatrix[c(pop1,pop2),])
		x = t(x)
		res = kmeans(x, 4, iter=1000)
		
		#Find which point is which (HH, HL, LH, LL) ?
		#=======================================================
		x = as.data.frame(res$centers)
		x$rank= NA
		x$sort1 = NA
		x$sort2 = NA
		for(i in 1:4){
			x[which(x[,pop1]==sort(x[,pop1])[i]),"sort1"] = i;
			x[which(x[,pop2]==sort(x[,pop2])[i]),"sort2"] = i;
		}
		#Complete rank column with HH, LL, HL, LH
		for(i in 1:4){
			if(x[i,"sort1"] <=2 & x[i,"sort2"] <= 2)
				x[i,"rank"] = "LL"
			else if(x[i,"sort1"] >=3 & x[i,"sort2"] >= 3)
				x[i,"rank"] = "HH"
			else if(x[i,"sort1"] <= 2)
				x[i,"rank"] = "LH"
			else x[i,"rank"] = "HL"
		}
		#=======================================================		
		#How much HH ? LL ? HL? LH ?
		#If 2 HH then 2 LL
		#If 2 LH then 2 HL
		#If centroids aligned, more than one center is HH
		if( dim(x[which(x$rank == "HH"),])[1] != 0){

			if( dim(x[which(x$rank == "HH"),])[1] >= 2 ){
				#At least 2 HH
				#case centroids aligned (increase)
				sumCentroid["LL", pop1] = sumCentroid["LL", pop1] + sum(x[which(x$rank == "LL"),pop1]) / 2 ;
				sumCentroid["LL", pop2] = sumCentroid["LL", pop2] + sum(x[which(x$rank == "LL"),pop2]) / 2 ;
				sumCentroid["HH", pop1] = sumCentroid["HH", pop1] + sum(x[which(x$rank == "HH"),pop1]) / 2 ;
				sumCentroid["HH", pop2] = sumCentroid["HH", pop2] + sum(x[which(x$rank == "HH"),pop2]) / 2 ;

				sumCentroid["LH", pop1] = sumCentroid["LH", pop1] + sum(x[which(x$rank == "LL"),pop1]) / 2 ;
				sumCentroid["LH", pop2] = sumCentroid["LH", pop2] + sum(x[which(x$rank == "HH"),pop2]) / 2 ;
				sumCentroid["HL", pop1] = sumCentroid["HL", pop1] + sum(x[which(x$rank == "LL"),pop2]) / 2 ;
				sumCentroid["HL", pop2] = sumCentroid["HL", pop2] + sum(x[which(x$rank == "HH"),pop1]) / 2 ;
			}else{
				# else exactly 1 HH, case centroids not aligned - good case :)
				sumCentroid["HH", pop1] = sumCentroid["HH", pop1] + x[which(x$rank == "HH"),pop1];
				sumCentroid["HH", pop2] = sumCentroid["HH", pop2] + x[which(x$rank == "HH"),pop2];
				sumCentroid["LL", pop1] = sumCentroid["LL", pop1] + x[which(x$rank == "LL"),pop1];
				sumCentroid["LL", pop2] = sumCentroid["LL", pop2] + x[which(x$rank == "LL"),pop2];
				sumCentroid["HL", pop1] = sumCentroid["HL", pop1] + x[which(x$rank == "HL"),pop1];
				sumCentroid["HL", pop2] = sumCentroid["HL", pop2] + x[which(x$rank == "HL"),pop2];
				sumCentroid["LH", pop1] = sumCentroid["LH", pop1] + x[which(x$rank == "LH"),pop1];
				sumCentroid["LH", pop2] = sumCentroid["LH", pop2] + x[which(x$rank == "LH"),pop2];
			}
		}else{
			#no HH, only LH, LH = case centroids aligned (decreasing)
			sumCentroid["LH", pop1] = sumCentroid["LH", pop1] + sum(x[which(x$rank == "LH"),pop1]) / 2 ;
			sumCentroid["LH", pop2] = sumCentroid["LH", pop2] + sum(x[which(x$rank == "LH"),pop2]) / 2 ;
			sumCentroid["HL", pop1] = sumCentroid["HL", pop1] + sum(x[which(x$rank == "HL"),pop1]) / 2 ;
			sumCentroid["HL", pop2] = sumCentroid["HL", pop2] + sum(x[which(x$rank == "HL"),pop2]) / 2 ;

			sumCentroid["HH", pop1] = sumCentroid["HH", pop1] + sum(x[which(x$rank == "HL"),pop1]) / 2 ;
			sumCentroid["HH", pop2] = sumCentroid["HH", pop2] + sum(x[which(x$rank == "LH"),pop2]) / 2 ;
			sumCentroid["LL", pop1] = sumCentroid["LL", pop1] + sum(x[which(x$rank == "LH"),pop1]) / 2 ;
			sumCentroid["LL", pop2] = sumCentroid["LL", pop2] + sum(x[which(x$rank == "HL"),pop2]) / 2 ;
		}
	}
	#Mean of 200 centroids
	sumCentroid = sumCentroid /200
	#Final kmeans call
	x = as.matrix(MCPcMatrix[c(pop1,pop2),])
	x = t(x)
	res = kmeans(x, 4,centers=sumCentroid, iter=100)
	return(res)
}

#===========================================================================================================================================

#PROBLEM, RESULT NOT UNIQUE, CANNOT USE THIS FUNCTION, RANDOM RESULTS
methodMajority = function(pop1, pop2, MCPcMatrix){

	i=0
	j=1
	allCenters=list(NULL)
	majorityCenters = c(1)
	x = as.matrix(MCPcMatrix[c(pop1,pop2),])
	x = t(x)
	res = kmeans(x, 4, iter = 2000)
	allCenters[[1]] = res$centers

	for(i in 1:1000){
		
		res = kmeans(x, 4, iter = 1000)
		predicted = 0
		
		#check if centers already found in a previous iteration
		for(j in 1:length(allCenters)){
			if(identical(sort(res$centers), sort(allCenters[[j]])) )
				predicted = j
		}
		#if yes, then increase counter
		if(predicted != 0){
			majorityCenters[j] = majorityCenters[j] + 1
		}
		#else, put new centers in the possibilities list
		else{
			allCenters[[length(allCenters)+1]] = res$centers
			majorityCenters =c(majorityCenters,1)
		}
	}

	#print(majorityCenters)
	#print(sum(majorityCenters))
	if(max(majorityCenters) > 2*sort(majorityCenters, T)[2] ){
		res = kmeans(x, 4,centers =allCenters[[which.max(majorityCenters)]],   iter = 100)
	}else{
		res = methodMajority(pop1, pop2, MCPcMatrix)		
	}
	return(res)
}

#===========================================================================================================================================

plotCentroidAndKM = function(pop1, pop2, MCPcMatrix, annotData, method){

	if(method =="corners"){
		clustData = methodCorners(pop1, pop2, MCPcMatrix)
	}else if(method =="mean"){
		clustData = methodMean(pop1, pop2, MCPcMatrix)
	}else if(method =="majority"){
		clustData = methodMajority(pop1, pop2, MCPcMatrix)		
	}

	annotDataMethod = cbind(clustData$cluster, annotData)
	colnames(annotDataMethod)[1] ="clust"

	#annotDataMethod = annotDataMethod[!duplicated(annotDataMethod$bcr_patient_barcode) & !duplicated(annotDataMethod$bcr_patient_barcode, fromLast=TRUE),]
	plotData = cbind(as.data.frame(clustData$cluster), t(MCPcMatrix)[,c(pop1, pop2)])
	colnames(plotData)[1] = "clust"

	#================== PREPARE KM ==========================
	chi = survdiff(Surv(annotDataMethod$OS.delay, annotDataMethod$OS.event)~annotDataMethod$clust)$chisq
	p = 1- pchisq(chi, 3) #4 groups, 4curves so 3 degrees of freedom
	#H0 : All survival curve are the same
	#P value at 3 and 5 years
	annotSurv3ans = annotDataMethod;
	annotSurv3ans[which(annotSurv3ans$OS.delay >=1095),"OS.event"] = 0 ;
	chi3ans = survdiff(Surv(annotSurv3ans$OS.delay, annotSurv3ans$OS.event)~annotSurv3ans$clust)$chisq
	p3ans = 1- pchisq(chi3ans, 3) #4 groups, 4curves so 3 degrees of freedom

	annotSurv5ans = annotDataMethod;
	annotSurv5ans[which(annotSurv3ans$OS.delay >=1825),"OS.event"] = 0 ;
	chi5ans = survdiff(Surv(annotSurv5ans$OS.delay, annotSurv5ans$OS.event)~annotSurv5ans$clust)$chisq
	p5ans = 1- pchisq(chi5ans, 3) #4 groups, 4curves so 3 degrees of freedom

	#================== PLOT KAPLAN MEIER ==========================
	#p <= 0.05 || p3ans <= 0.05 || p5ans <= 0.05
	if(1){

		#================== PLOT CENTROIDS ==========================
		plot(x=MCPcMatrix[pop1,], y=MCPcMatrix[pop2,], 
				xlab=pop1, ylab=pop2, col="white", pch=1,
				main=paste("GSE3538 -",pop1,"-",pop2,"\nMethod", method),
				sub=paste("Correlation", signif(cor(MCPcMatrix[pop1,], MCPcMatrix[pop2,]), digits = 5),
				" Pvalue", signif(cor.test(MCPcMatrix[pop1,], MCPcMatrix[pop2,])$p.value, digits = 5))
			);
		#Plot patients
		points(plotData[which(plotData$clust == 1),2:3], col="red", pch=1);
		points(plotData[which(plotData$clust == 2),2:3], col="black", pch=1);
		points(plotData[which(plotData$clust == 3),2:3], col="lightblue", pch=1);
		points(plotData[which(plotData$clust == 4),2:3], col="blue", pch=1);
		#Plot centroids
		points(t(clustData$centers[1,1:2]), col="red", pch=16, cex=1.5);
		points(t(clustData$centers[2,1:2]), col="black", pch=16, cex=1.5);
		points(t(clustData$centers[3,1:2]), col="lightblue", pch=16, cex=1.5);
		points(t(clustData$centers[4,1:2]), col="blue", pch=16, cex=1.5);

		plot(survfit(Surv(annotDataMethod[(annotDataMethod$clust==1),"OS.delay"], annotDataMethod[(annotDataMethod$clust==1),"OS.event"])~1 ),
				mark.time = T, conf.int=F,
				col="red", lwd=1.5,
				main=paste("GSE3538 -",pop1,"with", pop2, "\nMethod", method),
				ylab="Overall Survival", xlab = "Time in months",
				sub= paste("pValue", signif(p, digits=5)),
				cex.sub=0.9,
				xlim = c(0,3000)
				);
		lines(survfit(Surv(annotDataMethod[(annotDataMethod$clust==2),"OS.delay"], annotDataMethod[(annotDataMethod$clust==2),"OS.event"])~1 ),
				mark.time = T, conf.int=F,
				col="black", lwd=1.5) ;
		lines(survfit(Surv(annotDataMethod[(annotDataMethod$clust==3),"OS.delay"], annotDataMethod[(annotDataMethod$clust==3),"OS.event"])~1 ),
			mark.time = T, conf.int=F,
			col="lightblue", lwd=1.5); 
		lines(survfit(Surv(annotDataMethod[(annotDataMethod$clust==4),"OS.delay"], annotDataMethod[(annotDataMethod$clust==4),"OS.event"])~1 ),
			mark.time = T, conf.int=F,
			col="blue", lwd=1.5) ;
		legend(
			x = 1095,
			y=1,
			legend=paste("p at 3y\n",signif(p3ans, digits=4)),
			cex=0.7,
			box.lty = 0	)
		legend(
			x = 1825,
			y=1,
			legend=paste("p at 5y\n",signif(p5ans, digits=4)),
			cex=0.7,
			box.lty = 0)
		abline(v=1095, lty="dashed") #36 months, OS at 3years
		abline(v=1825, lty="dashed") #60 months, OS at 5years	

		legend(x="bottomleft",
				legend= c(paste("n1 =", clustData$size[1]),
						paste("n2 =", clustData$size[2]),
						paste("n3 =", clustData$size[3]),
						paste("n4 =", clustData$size[4])),	
		inset=0.01, col=c("red", "black", "lightblue", "blue"),lty=1, lwd=1.5,cex=0.8,box.lty = 0)
		
	}
}

#===========================================================================================================================================


#Recup arg : pop1 pop2 estimate et annot

# Faire le Clust avec chaque methode et ploter les points (etiquettes : 1 2 3 4, 4 couleurs)
# kmean a 1000 iter
#	Methode corners : initial centers = coins du nuage de point
#	Methode mean : faire la moyenne sur 200 iterations
#	Methode majority : etat d'equilibre majoritaire

#Les trois methodes renvoient un objet kmean


#	!!!		Penser à faire des modifs dans le titre des plot selon si la cohorte est TCGA KIRC ou GSE3538 	 !!!!

#===========================================================================================================================================
#===========================================================================================================================================
#	TCGA RUN

annotData = annotTCGA_KIRC[,c("OS.delay","OS.event", "MFS.delay", "MFS.event","bcr_patient_barcode")];

pdf("CouplesImmunoCells_TCGA-KIRC_CornersOnly.pdf")
for( i in 1:(dim(estimateTCGA_KIRC)[1]-1) ){
	for( j in (i+1):dim(estimateTCGA_KIRC)[1] ){

		pop1 = rownames(estimateTCGA_KIRC)[i];
		pop2 = rownames(estimateTCGA_KIRC)[j];
		twoCellsToKM(pop1, pop2, estimateTCGA_KIRC, annotData)
	}
}
dev.off()

#===========================================================================================================================================
#===========================================================================================================================================
#GSE RUN

pdf("coupleImmuno_GSE3538_CornersOnly_OSSevent.pdf")

estimateGSE3538save = estimateGSE3538
estimateGSE3538[estimateGSE3538=="NaN"] = NA
estimateGSE3538 = estimateGSE3538[,complete.cases(t(estimateGSE3538))]
annotDataSave = annotSurv

for(i in 1:(dim(estimateGSE3538save)[1]-1) ){
	for(j in (i+1):dim(estimateGSE3538save)[1]){

		pop1 = rownames(estimateGSE3538save)[i];
		pop2 = rownames(estimateGSE3538save)[j];

		estimateGSE3538 = estimateGSE3538save[, !is.na(estimateGSE3538save[pop1,]) & !is.na(estimateGSE3538save[pop2,])]
		annotData = annotDataSave[colnames(estimateGSE3538),]
	
		twoCellsToKM(pop1, pop2, estimateGSE3538, annotData)
	}
}
dev.off()

#===========================================================================================================================================
#===========================================================================================================================================
#	TCGA RUN

annotData = annotTCGA_LIHC[,c("OS.delay","OS.event","bcr_patient_barcode")];
estimateTCGA_LIHC = MCPcounter.estimate(expTCGA_LIHC, featuresType="HUGO_symbols")

pdf("CouplesImmunoCells_TCGA-LIHC.pdf")
for( i in 1:(dim(estimateTCGA_LIHC)[1]-1) ){
	for( j in (i+1):dim(estimateTCGA_LIHC)[1] ){

		pop1 = rownames(estimateTCGA_LIHC)[i];
		pop2 = rownames(estimateTCGA_LIHC)[j];
		twoCellsToKM(pop1, pop2, estimateTCGA_LIHC, annotData)
	}
}
dev.off()
