\name{LDsnpR-package}
\alias{LDsnpR-package}
\docType{package}
\title{
\packageTitle{LDsnpR}
}
\description{
\packageDescription{LDsnpR}
}
\details{

\tabular{ll}{
Package: \tab LDsnpR \cr
Version: \tab 1.0.1 \cr
Date: \tab 2011-09-02 \cr
License: \tab FreeBSD \cr
Depends:  \tab IRanges (>= 1.14.0), rhdf5 (>=2.12.0) \cr
URL:  \tab http://services.cbu.uib.no/software/ldsnpr/  \cr
}

LDsnpR is  a flexible  open-source R-package, which efficiently assigns SNPs to genes, or user-defined bins, based both on their physical position and on their pairwise LD with other SNPs.

LDsnpR allows to compute gene-wise scores based on the combined p-values of all the markers assigned to the gene bins with or without LD binning, and numerous summary statistics for computation of joint p-values, which can be extended by user-defined functions. The package supports formats of the widely-used \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{Plink} genotyping software as input data and can also generate PLINK-set format of SNP identifiers for subsequent analysis. LD calculations are based on the four initial HapMap populations. The population code with which to carry out the LD-binning is a parameter to the analysis.

The design and execution of genome-wide association studies (GWASs) are critically dependent on detailed knowledge of the pattern of linkage disequilibrium (LD) in the human genome, as characterised by the HapMap project and The 1000 Genomes project. A GWAS generates a list of variants, usually single-nucleotide polymorphisms (SNPs), ranked according to the significance of their association to the trait of interest. Downstream analyses generally focus on the gene(s) that are physically closest to the aforementioned SNPs, failing to account for their LD profile with other SNPs.  Modern array-based platforms allow for the parallel genotyping of millions of SNPs. 

For selecting the markers present on array platforms, LD information is often taken into account, to allow for example interrogating genotypes of regions in high LD using an optimally reduced set of markers. Therefore, it can be beneficial to account for LD in the subsequent analysis process (Christoforou et. al., 2012).  

The method is implemented in the \href{http://www.r-project.org/}{R-language}  and uses the \href{http://www.bioconductor.org/packages/release/bioc/html/IRanges.html}{IRanges} package from Bioconductor for fast interval overlapping. To facilitate efficient storage and retrieval of LD-data, we have converted the compressed text files found in  \href{http://hapmap.ncbi.nlm.nih.gov/}{HapMap} into a compressed \href{http://www.hdfgroup.org/HDF5/}{HDF5} binary data format. That way, we are able to provide LD data for all chromosomes and four populations in a single file of only 4.9 GB (about 20\% of the original size).  Computation of LD-based binning of 1 Million SNPs into all human Ensembl genes can be carried out with less than 2 GB of memory.
Thanks to the efficiency and maturity of HDF5 by the \href{http://www.hdfgroup.org/HDF5/}{the HDF5 group}, it lends itself to efficient storage of LD-data on a much larger scale. 

As an example, we have presented gene-based analysis of public avaliable GWAS: \href{http://www.med.unc.edu/pgc/}{schizophrenia GWAS}. We have converted data on linkage-disequilibrium  from Europien cohort  \href{http://www.1000genomes.org/}{the 1000 Genomes project}  to HDF5 format, contributed the converted data, and integrated them into our analysis work-flow.

}
\author{
\packageAuthor{LDsnpR}

Maintainer: \packageMaintainer{LDsnpR}
}
\references{
Andrea Christoforou, Michael Dondrup, Morten Mattingsdal, Manuel Mattheisen, Sudheer Giddaluru, Markus M. N?then, Marcella Rietschel, Sven Cichon, Srdjan Djurovic, Ole A. Andreassen, Inge Jonassen, Vidar M. Steen, Pal Puntervoll and Stephanie Le Hellard.
Linkage-Disequilibrium-Based Binning Affects the Interpretation of GWASs. The American Journal of Human Genetics  (2012) 10.1016/j.ajhg.2012.02.025
}

\keyword{ GWAS
SNP
LD }
\seealso{
\code{\link[IRanges:IRanges]{IRanges}},
\code{\link[rhdf5:rhdf5]{rhdf5}},
\code{\link[LDsnpR:snp.ld.analysis]{snp.ld.analysis}} 

}
