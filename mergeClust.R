estimateTCGA_KIRC = MCPcounter.estimate(expTCGA_KIRC, featuresType = "HUGO_symbols");
estimateGSE3538 = MCPcounter.estimate(expGSE3538, featuresType = "HUGO_symbols");
estimateGSE3538[estimateGSE3538=="NaN"] = NA
estimateGSE3538save = estimateGSE3538

#============================== START ========================

pop1 = "Endothelial cells"
pop2 = "Fibroblasts"


resTCGA = methodMean(pop1, pop2, estimateTCGA_KIRC )
estimateGSE3538 = estimateGSE3538save[, !is.na(estimateGSE3538save[pop1,]) & !is.na(estimateGSE3538save[pop2,])]
resGSE = methodMean(pop1, pop2, estimateGSE3538 )


resTCGA = cbind(resTCGA$cluster, annotTCGA_KIRC[,c("OS.delay","OS.event")])
colnames(resTCGA)[1] = "clust"
resGSE = cbind(resGSE$cluster, annotGSE3538[colnames(estimateGSE3538),c("OS.delay","OS.event")])
colnames(resGSE)[1] = "clust"

#T cells + Fibro
#resTot = rbind(resTCGA, resGSE[which(resGSE$clust == 3 | resGSE$clust == 4 ),])

#B + Endo
#resGSE[which(resGSE$clust == 4),"clust"] = 3
#resTCGA[which(resTCGA$clust == 4),"clust"] = 2
#resTot = rbind(resTCGA[which(resTCGA$clust != 4),], resGSE[which(resGSE$clust != 4),])


resTot = rbind(resTCGA, resGSE)
annotDataMethod = resTot

#pdf("T_Myelo_KidneyMerged.pdf")
#====================================================================================================
	annotDataMethod[which(annotDataMethod$OS.delay >=80),"OS.event"] = 0 ;
	chi = survdiff(Surv(annotDataMethod$OS.delay, annotDataMethod$OS.event)~annotDataMethod$clust)$chisq
	p = 1- pchisq(chi, 3) #4 groups, 4curves so 3 degrees of freedom
	#H0 : All survival curve are the same

	#================== PLOT KAPLAN MEIER ==========================

		plot(survfit(Surv(annotDataMethod[(annotDataMethod$clust==1),"OS.delay"], annotDataMethod[(annotDataMethod$clust==1),"OS.event"])~1 ),
				mark.time = T, conf.int=F,
				col="red", lwd=1.5,
				main=paste("Merged Cohorts Kidney\n",pop1,"with", pop2, "\nMethod Coners"),
				ylab="Overall Survival", xlab = "Time in months",
				sub= paste("pValue", signif(p, digits=5)),
				cex.sub=0.9,
				xlim = c(0,100)
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

		abline(v=80, lty="dashed") #end of study at 80 months
		legend(
			x = 81,
			y=1,
			legend="End of study\nAll censored",
			cex=0.7,
			box.lty = 0)

		legend(x="bottomleft",
				legend= c(paste("n1 =", sum(annotDataMethod$clust == 1), ""),
						paste("n2 =", sum(annotDataMethod$clust == 2), ""),
						paste("n3 =", sum(annotDataMethod$clust == 3), ""),
						paste("n4 =", sum(annotDataMethod$clust == 4), "")),
		inset=0.02, col=c("red", "black", "lightblue", "blue"),lty=1, lwd=1.5,cex=0.8,box.lty = 0)


#dev.off()








