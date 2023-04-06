Simulate GWAS data
================

Download packages required to create simulation of GWAS data.

``` r
library(VariantAnnotation)
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

    ##         gene_name         ID Allele1 Allele2 chr        bp      MarkerName
    ## 1 ENSG00000102081 rs76943034       A       T   9  39627072  9:39627072:A:T
    ## 2 ENSG00000102081 rs16312248       T       A   3  52399547  3:52399547:T:A
    ## 3 ENSG00000102081 rs29618759       G       C   5 108153621 5:108153621:G:C
    ## 4 ENSG00000102081 rs91797184       A       T   8  49485561  8:49485561:A:T
    ## 5 ENSG00000102081 rs26432326       G       C  18  84874009 18:84874009:G:C
    ## 6 ENSG00000102081 rs93848964       C       G   2  40100752  2:40100752:C:G
    ##          AF        ES         SE  SS         p
    ## 1 0.3710578 0.2091585 0.23896056 100 0.3814191
    ## 2 0.1430606 0.4828864 1.77042508 100 0.7850441
    ## 3 0.1263383 0.0529422 0.04633079 100 0.2531631
    ## 4 0.4529944 0.6188569 8.55818183 100 0.9423538
    ## 5 0.3251773 0.5006739 0.45518604 100 0.2713616
    ## 6 0.3389784 0.3565727 0.42132797 100 0.3973817
