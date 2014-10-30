## input:
## data is a data.frame with columns
## p=single snp p value
## MAF=minor allele frequency
## Immuno=ImmunoChip SNP id
## dbSNP=dbSNP SNP id
## Start=SNP position, build 37

## requires tab2.csv: Table 2 extracted from
## Eyre et al Nat Genet. 2012 Dec;44(12):1336-40. doi: 10.1038/ng.2462. Epub 2012 Nov 11.
## High-density genetic mapping identifies new susceptibility loci for rheumatoid arthritis.


## calculate Bayes Factors
library(randomFunctions) ## get from https://github.com/chr1swallace/random-functions
data$BF <- -abf(p=data$p, maf=data$MAF, n0=13838, n1=33742)

## sanity check
library(ggplot2)
ggplot(data,aes(x=-log10(data$p),y=BF)) + geom_point()

## top hits from RA paper
tab2 = read.table("table2-mod.csv",header=TRUE,sep="\t",as.is=TRUE)
## add chr, position
tab2$index.snp <- gsub("[abcdef].*$|/.*","",tab2$SNP)
tab2$chr <- gsub("[pq].*","",tab2$Chr)
table(immuno=tab2$index.snp %in% data$Immuno,dbsnp=tab2$index.snp %in% data$dbSNP)
tab2$position <- NA
m <- match(tab2$index.snp,data$Immuno)
tab2$position[!is.na(m)] <- data$Start[m[!is.na(m)]]
m <- match(tab2$index.snp,data$dbSNP)
tab2$position[!is.na(m)] <- data$Start[m[!is.na(m)]]
tab2$region.start <- as.numeric(gsub(",|\342.*","",tab2$LD))
tab2$region.end <- as.numeric(gsub(",|.*\223","",tab2$LD))

## exclude chr X
tab2 <- subset(tab2,chr!="X")

## step through, line by line, to create credible sets
library(coloc) # logsum function, get library from CRAN

## some functions

##' Summarize credible set for a region
##'
##' @title dosumm
##' @param dsub a data.frame with one row per SNP in a region, and a
##' cpp column giving the cumulative posterior probabilities that SNPs
##' on or below a given line contain the causal variant for this
##' region
##' @return a vector giving the number of SNPs within cpp<0.95 or
##' 0.99, the number of SNPs in the region, and the length of the
##' region.
dosumm <- function(dsub) { # dsub is sorted on cpp
  n95 <- which(dsub$cpp>=0.95)[1]
  n99 <- which(dsub$cpp>=0.99)[1]
  return(c(n95=n95,n99=n99,
           region.snps=nrow(dsub),region.length=diff(range(dsub$Start))))
}
##' process a single row in tab2 to generate PP for SNPs in the region indexed by that SNP
##'
##' @title pproc
##' @param data data.frame described at top of this file
##' @param i row of tab2
##' @param window window (in bp) by which region in tab2 should be extended
##' @return a data.frame with one row per SNP giving posterior probs
##' @author chris
    pproc <- function(data,i,window=0) {
        dsub=subset(data,Chromosome==tab2[i,"chr"] & Start>=tab2[i,"region.start"] &
                   Start<=tab2[i,"region.end"] & !is.na(MAF))
        dsub$index <- dsub$Start==tab2[i,"position"]
        if(sum(dsub$index)!=1)
            stop()
        dsub$pp <- exp(dsub$BF - coloc:::logsum(dsub$BF))
        dsub <- dsub[order(dsub$pp,decreasing=TRUE),]
        dsub$cpp <- cumsum(dsub$pp)
        return(dsub)
    }

## objects to hold results
cred99 <- cred95 <- credw99 <- credw95 <- summ <- creddata <- vector("list",nrow(tab2))

## bp to add to regions defined in paper for "wide" results
window <- 50000

for(i in seq_along(cred99)) {

 ## SNPs from regions defined in paper
 dsub <- pproc(data,i)
 ## wide version, add <window> to regions defined in paper
 dsubw <- pproc(data,i,window=window)
 
 ## collate results
 sm <- dosumm(dsub)
 smw <- dosumm(dsubw)
 cred95[[i]] <- dsub[1:sm["n95"],
                     c("Immunochip.ID","dbSNP.ID","Chromosome","Start","MAF","p","index","pp","cpp")]
 credw95[[i]] <- dsubw[1:smw["n95"],
                       c("Immunochip.ID","dbSNP.ID","Chromosome","Start","MAF","p","index","pp","cpp")]
 cred99[[i]] <- dsub[1:sm["n99"],
                     c("Immunochip.ID","dbSNP.ID","Chromosome","Start","MAF","p","index","pp","cpp")]
 credw99[[i]] <- dsubw[1:smw["n99"],
                       c("Immunochip.ID","dbSNP.ID","Chromosome","Start","MAF","p","index","pp","cpp")]
 summ[[i]] <- c(narrow=sm,wide=smw)
 print(summ[[i]])
 creddata[[i]] <- dsubw
}

## create result table
result <- cbind(tab2,do.call("rbind",summ))
result$cred95.start <- sapply(cred95,function(x) min(x$Start))
result$cred95.end <- sapply(cred95,function(x) max(x$Start))
result$cred99.start <- sapply(cred99,function(x) min(x$Start))
result$cred99.end <- sapply(cred99,function(x) max(x$Start))
result$credw95.start <- sapply(credw95,function(x) min(x$Start))
result$credw95.end <- sapply(credw95,function(x) max(x$Start))
result$credw99.start <- sapply(credw99,function(x) min(x$Start))
result$credw99.end <- sapply(credw99,function(x) max(x$Start))
result$P <- sub("\303\227","x",result$P)
result$P <- sub("\342\210\222","-",result$P)
result <- result[,c("index.snp","Gene","Chr.","MAF","Risk.allele","P","OR","position",
                    "region.start","region.end",
                    "narrow.region.snps","narrow.n95","narrow.n99",
                    "wide.region.snps","wide.n95","wide.n99",
                    "cred95.start","cred95.end","cred99.start","cred99.end",
                    "credw95.start","credw95.end","credw99.start","credw99.end")]

colnames(result) <- sub("credw","wide.credible.",colnames(result))
colnames(result) <- sub("cred9","narrow.credible.9",colnames(result))

## create details
for(i in seq_along(credw99)) {
  credw99[[i]]$index.snp <- result[i,"index.snp"]
}
credw99.long <- do.call("rbind",credw99)
for(i in seq_along(cred99)) {
  cred99[[i]]$index.snp <- result[i,"index.snp"]
}
cred99.long <- do.call("rbind",cred99)

## write to disk
write.table(cred99.long,"credible-sets.tab",col.names=TRUE,row.names=FALSE,sep="\t")
write.table(credw99.long,"credible-sets-wide.tab",col.names=TRUE,row.names=FALSE,sep="\t")
write.table(result,"target-regions.tab",col.names=TRUE,row.names=FALSE,sep="\t")

## Optional: granges objects
library(GenomicRanges)
gr95 <- GRanges(seqnames=paste0("chr",sub("pq.*","",result$Chr.)),
                IRanges(start=result$narrow.credible.95.start,end=result$narrow.credible.95.end))
wgr95 <- GRanges(seqnames=paste0("chr",sub("pq.*","",result$Chr.)),
                IRanges(start=result$wide.credible.95.start,end=result$wide.credible.95.end))
gr99 <- GRanges(seqnames=paste0("chr",sub("pq.*","",result$Chr.)),
                IRanges(start=result$narrow.credible.99.start,end=result$narrow.credible.99.end))
wgr99 <- GRanges(seqnames=paste0("chr",sub("pq.*","",result$Chr.)),
                IRanges(start=result$wide.credible.99.start,end=result$wide.credible.99.end))

sum(width(gr95))==sum(width(reduce(gr95))) ## all TRUE
sum(width(wgr99))==sum(width(reduce(wgr99))) ## all TRUE
