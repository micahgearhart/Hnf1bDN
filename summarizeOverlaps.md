``` r
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(DESeq2)
```

``` r
load("ens84_mouse.rdata")
ts<-format(Sys.time(), "%a_%b_%d_%Y_%H%M")

(fls <- list.files("GRCm38", pattern=".bam$",full=TRUE))
bamlst <- BamFileList(fls,yieldSize=1e7)
register(MulticoreParam(workers=24))
genehits <- summarizeOverlaps(ens84,bamlst,mode="Union",singleEnd=FALSE,ignore.strand=TRUE)
save(genehits,file=paste0("igarash0_Project_004_genehits_ens84",ts,".rdata"))
(n<-apply(assays(genehits)$counts,2,sum))
stopifnot(sum(n) > 0)
```

``` r
load("igarash0_Project_004_genehits_ens84Thu_Jun_30_2016_2351.rdata")
cds <- DESeqDataSet(genehits,design=~1)
cds <- estimateSizeFactors(cds)
cds$filename<-rownames(colData(cds))
```

``` r
#Define Function for system calls to convert Bam files to bigwig files
x<-colData(cds)
bam2bw <- function (i) {
 # system2("touch",paste0("myfile","_",x,".txt"))
  system2("bedtools", paste0("genomecov -scale ",1/x[i,"sizeFactor"],
                 " -split -bg -ibam GRCm38/",
                 x[i,"filename"]),stdout=paste0(rownames(x)[i],".bedGraph"),
                 ,stderr=paste0(i,"_testerr.out"))

  system2("/home/bardwell/shared/bedGraphToBigWig",
          paste0(rownames(x)[i],".bedGraph ",
                 "/home/bardwell/shared/STAR_GENOME/GRCm38/chrNameLength.txt ",
                 rownames(x)[i],".bigWig"),
                 stdout=paste0(i,"_test2.out"),stderr=paste0(i,"_testerr2.out"))
}

bplapply(1:nrow(colData(cds)),bam2bw)
```

    ## [[1]]
    ## [1] 0
    ## 
    ## [[2]]
    ## [1] 0
    ## 
    ## [[3]]
    ## [1] 0
    ## 
    ## [[4]]
    ## [1] 0
    ## 
    ## [[5]]
    ## [1] 0
    ## 
    ## [[6]]
    ## [1] 0
    ## 
    ## [[7]]
    ## [1] 0
    ## 
    ## [[8]]
    ## [1] 0
    ## 
    ## [[9]]
    ## [1] 0
    ## 
    ## [[10]]
    ## [1] 0
    ## 
    ## [[11]]
    ## [1] 0
    ## 
    ## [[12]]
    ## [1] 0
