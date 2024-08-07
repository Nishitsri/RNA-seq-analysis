##load relevant objects from differential expression analysis

suppressPackageStartupMessages(library(edgeR))
load("Robjects/DE.Rdata")

results <- as.data.frame(topTags(lrt.BvsL,n = Inf))
results

##visualise DE analysis using plotsmear
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)

##Add annotation to edgeR results
source("http://www.bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
# For Human
biocLite("org.Hs.eg.db")
library(org.Mm.eg.db)

#extract columns from annotation database
columns(org.Mm.eg.db)

##filter the database
keytypes(org.Mm.eg.db)

##check the ENTREZid and make sure that keys are valid
keys(org.Mm.eg.db, keytype="ENTREZID")[1:10]
## Build up the query step-by-step
my.keys <- c("50916", "110308","12293")
my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID")
all(my.keys %in% keys(org.Mm.eg.db, keytype="ENTREZID"))

#Get genename as well
ann <- select(org.Mm.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
ann

#bind annotation to results dataframe

results.annotated <- cbind(results, ann)
results.annotated

#save the results so that it can be opened in excel
write.csv(results.annotated,file="B.PregVsLacResults.csv",row.names=FALSE)


##Add genomic location to the annotation table
source("http://www.bioconductor.org/biocLite.R")
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")

## For Humans
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

#Load the library and create a new object
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
tx <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(tx)

#Genomic intervals
library(GenomicRanges)
simple.range <-GRanges("1", IRanges(start=1000,end=2000))
simple.range

keys <- c("50916","110308","12293")
genePos <- select(tx, keys=keys,
                  keytype = "GENEID",
                  columns=c("EXONCHROM","EXONSTART","EXONEND")
)

geneRanges <- GRanges(genePos$EXONCHROM, IRanges(genePos$EXONSTART,genePos$EXONEND), GENEID=genePos$GENEID)
geneRanges


#find overlaps

findOverlaps(my.ranges,geneRanges)
seqlevelsStyle(geneRanges)
seqlevelsStyle(simple.range)


#retreiving gene coordinates as genomic ranges
exo <- exonsBy(tx,"gene")
exo


#exporting tracks
sigGenes <- results.annotated[detags,]
sigGenes

exoRanges <- unlist(range(exo))
sigRegions <- exoRanges[na.omit(match(sigGenes$ENTREZID, names(exoRanges)))]
sigRegions

#Restricting data to genes located between chr 1 and 19 and sex chromosomes
sigRegions <- keepSeqlevels(sigRegions, paste0("chr", c(1:19,"X","Y")))
