Simulate GWAS data
================

Download packages required to create simulation of GWAS data.

``` r
library(VariantAnnotation)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomeInfoDb

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: Rsamtools

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## 
    ## Attaching package: 'VariantAnnotation'

    ## The following object is masked from 'package:base':
    ## 
    ##     tabulate

``` r
library(gprofiler2)
```

Generate test data - this is a simple generation of random data to
replicate GWAS data. This is not any form of genetic data and the
simulation does not adjust for linkeage disequilibrium or heritability.

``` r
ID <- paste0("rs", floor(runif(100, min=1, max=1e8))) #generate random rsid
base <- c("A","T","C","G") #list of bases
Allele1 <- sample(base, 100, replace=TRUE) #generate random base
Allele2 <- c() #empty list
for (i in Allele1){
if (i == "A"){
Allele2 <- append(Allele2, "T")}
else if (i == "T"){
Allele2 <- append(Allele2, "A")}
else if (i == "C"){
Allele2 <- append(Allele2, "G")}
else Allele2 <- append(Allele2, "C")} #generate base dependening on Allele1 base
chr <- floor(runif(100, min=1, max=22)) #generate chromsome location
bp <- floor(runif(100, min=1, max=14e7)) #generate base position
MarkerName= paste(chr, bp, Allele1, Allele2, sep=":") #create marker name
AF <- abs(runif(100, 0.05, 0.5)) #generate minor allele freqency
ES <- runif(100, 0, 1)#generate effect estimate
SS <- 100 #sample size of 100
p <- runif(100, min=0, max=1) #generate p-value
SE <- abs(ES)/qnorm(1-p/2) #this equation was confirmed using Z <- ES/SE and check <- 2*pnorm( Z, lower=F )
gene <- random_query(organism = "hsapiens")[1] #generate random gene name
gene_name <- gene
df <- data.frame(gene_name, ID, Allele1, Allele2, chr, bp, MarkerName, AF, ES, SE, SS, p) #create dataframe of all values for GWAS
```

See output of simulated GWAS data

``` r
head(df)
```

    ##         gene_name         ID Allele1 Allele2 chr        bp       MarkerName
    ## 1 ENSG00000123364 rs68538093       A       T   8  78536096   8:78536096:A:T
    ## 2 ENSG00000123364 rs69461568       C       G  20 104510935 20:104510935:C:G
    ## 3 ENSG00000123364  rs7745278       T       A   1 121902851  1:121902851:T:A
    ## 4 ENSG00000123364 rs59582645       A       T   9  16970128   9:16970128:A:T
    ## 5 ENSG00000123364 rs73077609       T       A  12   3627487   12:3627487:T:A
    ## 6 ENSG00000123364 rs70192785       T       A  18 119587853 18:119587853:T:A
    ##          AF        ES         SE  SS          p
    ## 1 0.4138154 0.3123727  0.1629037 100 0.05517061
    ## 2 0.3626379 0.6030055 49.7473572 100 0.99032879
    ## 3 0.4353666 0.6086147  0.3860676 100 0.11492310
    ## 4 0.4875335 0.2382613  0.2189072 100 0.27641314
    ## 5 0.3842818 0.6843204  1.0713205 100 0.52297682
    ## 6 0.4503532 0.9103380  2.5134714 100 0.71721487
