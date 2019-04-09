# Original code wrote by Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html
# Based on version: Tuesday, April 19, 2011

# modified by Guanglong Jiang
# Last updated: 06/19/2012


# R code for making manhattan plots and QQ plots from plink output files. 
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now. 

###############################################################################################################################
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
is.loaded <- function(mypkg) is.element(mypkg, loadedNamespaces()) 

## simulate data for testing.
set.seed(42)
nchr=26
nsnps=1000
testdata=data.frame(
    SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
    CHR=rep(1:nchr,each=nsnps), 
    BP=rep(1:nsnps,nchr), 
	P=runif(nchr*nsnps)
)
annotatesnps <- testdata$SNP[7550:7750]

###############################################################################################################################
# manhattan plot using base graphics
manhattan = function(dataframe, myylab="", colors=c("gray10", "gray50"), ymin=0, ymax="max", cex.x.axis=0.7, limitchromosomes=1:26, suggestiveline=-log10(1e-6), genomewideline=-log10(5e-8), annotate=NULL, ...) {
    cat(myylab)
	d=dataframe
    if (any(limitchromosomes)) 
		d=d[d$CHR %in% limitchromosomes, ] #remove SNPs whose CHR outside of 1:26
	d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
	cat(nrow(d), "SNPs are plotted\n")
	
	labchr = sort(unique(d$CHR))
	labchr = ifelse(labchr<=22, labchr, ifelse(labchr==23, "X", ifelse(labchr==24, "Y", ifelse(labchr==25,"XY", ifelse(labchr==26, "MT",labchr)))))

    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) 
		stop("Make sure your data frame contains columns CHR, BP, and P")

	d$logp = -log10(d$P)
	
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors, length=length(unique(d$CHR)))
    if (ymax=="max") 
		ymax<-ceiling(max(d$logp))
    if (ymin=="min")
		ymin=floor(min(d$logp))
		
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        #ticks=floor(length(d$pos))/2+1
        ticks = max(d$pos)/2+1
    } else {
        for (i in 1:numchroms) {
			mychr = sort(unique(d$CHR))[i]
        	if (i==1) {
    			d[d$CHR==mychr, ]$pos = d[d$CHR==mychr, ]$BP
				ticks = max(d[d$CHR==mychr,]$pos)/2
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==(unique(d$CHR)[i-1]))$BP, 1)
    			d[d$CHR==mychr, ]$pos = d[d$CHR==mychr, ]$BP+lastbase
				ticks = c(ticks, lastbase+tail(d[d$CHR==mychr,]$BP,1)/2)
    		}
    	}
    }
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), cex=0.4,cex.lab=1,cex.axis=1, ...))
    }	else {
        with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", cex=0.4,cex.lab=1,cex.axis=1, ...))
        axis(1, at=ticks, lab=labchr, las=2, lwd=0, lwd.ticks=1, line=-0.3, cex.axis=cex.x.axis)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
        }
    }
   
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...)) 
    }
    
    if (suggestiveline)
		abline(h=suggestiveline, col="gray40", lty=5)
    if (genomewideline)
		abline(h=genomewideline, col="gray40", lty=5)
}

###############################################################################################################################
# Base graphics qq plot
qq = function(pvector, ...) {
	if (!is.numeric(pvector)) stop("Error! P value vector is not numeric.")
	pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
	
	if(!is.installed("gap")){
		o = -log10(sort(pvector,decreasing=F))
		e = -log10( ppoints(length(pvector) )) # e = -log10( 1:length(o)/length(o) )
		plot(e,o,pch=19,cex=1, xlab=expression(-log[10](italic(Expected))), ylab=expression(-log[10](italic(Observed))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
		abline(0,1,col="red")
	}else{
		if(!is.loaded("gap")){
			library(gap)
		}
		r = gcontrol2(pvector)
		title(paste("lamda=",round(r$lambda, digits=3)))
		return(r$lambda)
	}
}

###############################################################################################################################
#par(ask=T)

#manhattan plot
# manhattan(testdata) 

#plot options
# manhattan(testdata, colors=1:8, pch=20, ymax=8) 

#highlighting
# snps_to_highlight <- testdata[testdata$P<1e-4,1] 
# manhattan(testdata, annotate=snps_to_highlight, pch=20, main="Manhattan Plot") 

#zoom in
# manhattan(subset(testdata, CHR==2), pch=20, main="Chromosome 2") 

#QQ plot of P 
# qq(testdata$P) 

#par(ask=F)

###############################################################################################################################
locuszoom = function(filename, refsnp, markercol="SNP", pvalcol="P", flank="400kb", delim=",", ...){
    if(missing("refsnp")){
        stop("Reference SNP(refsnp) can not be blank for locuszoom:\n \t\t\t\tExample: locuszom(input-file, refsnp=\"rs123\")")
	}
    cat("load data from", filename,"...\n")
    fileheader = read.csv(filename, header=T, nrow=1)
    if(!("SNP" %in% names(fileheader))){
        stop('Column header for SNP id must be "SNP" in file ', filename,"!")
    }else if(!("P" %in% names(fileheader))){
        stop('Column header for P value must be "P" in file ', filename,"!")
    }
    
	cat("Locuszoom will plot 400kb flanking variants by default. \nUse the argument flank=\"600kb\" if need to show further flanking variants.\n")
    
    system(paste("locuszoom --metal", filename, "--markercol ", markercol, "--pvalcol", pvalcol, 
    "--refsnp", refsnp, "--flank", flank,
    "--delim", delim,
    ...))
}


