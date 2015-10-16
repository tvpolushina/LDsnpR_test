#' Example  of gene-based scoring for SCZ GWAS on chromosome 18.
#' 
#'@description This function wraps the example of a SNP analysis on chromosome 18. 
#'@details   example.genome.url = system.file("data/ENSEMBL66_gene_bld37_example.txt", package="LDsnpR")   example.snpdata.url = system.file("data/pgc2_scz_chr18_example.txt", package="LDsnpR"); example.ld.file  = system.file("data/chr18_LD1kgEUR.h5", package="LDsnpR") 
#'@export

example<-function(){
  example.genome.url <-   system.file("data/ENSEMBL66_gene_bld37_example.txt", package="LDsnpR")
  example.snpdata.url <-  system.file("data/pgc2_scz_chr18_example.txt", package="LDsnpR") 
  example.ld.file  <- system.file("data/chr18_LD1kgEUR.h5", package="LDsnpR") 
  result <- snp.ld.analysis(snpdata.url = example.snpdata.url,
                            genome.url = example.genome.url, 
                            ld.data.hdf.url = example.ld.file,
                            use.position = TRUE,
                            ld.rho.cutoff=0.8,
                            comparator=">=",
                            full.match=FALSE,
                            flank.genes.left=10000,
                            flank.genes.right=10000,
                            population = "EUR",
                            scoring.function = c("slowscore"),
                            ld.structure = TRUE,
                            correction.type = c("none"), 
                            genome.name="hg19",
                            outfile.path=NULL) 
return(result)                            
}

  