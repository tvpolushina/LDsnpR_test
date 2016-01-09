#' Implementation of the SNP and LD data analysis using IRanges 
#' 
#' @description This function wraps the work-flow of a SNP analysis. SNPs are assigned to genes using findOverlaps. Then, gene-specific scores are computed as combined scores from the SNP-pvalues.
#' 
#' 
#' @param snpdata.url  character: A URL or file to read the SNPs in tab separated format.
#' @param genome.url  character: A URL or file to read the Genome in tab separated format.
#' @param ld.data.hdf.url	character: A URL or file to read the optional LD-data as a HDF5 file.
#' @param include.ld.data	logical: Should LD data be included in scoring? Takes more time and memory (default: FALSE)
#' @param population	Choose a population code, character vector of length one. The population code corresponds to the population codes in the LD data file. (see Details) (default: "EUR")
#' @param use.position logical: Should markers be matched by genomic postions or rs ids? 
#' @param full.match	Do a full crossmatching of LD data, only if include.ld.data == TRUE. See Details. (default: TRUE)
#' @param ld.rho.cutoff	numeric: Cutoff value for the LD correlation read from HDF file (Default: 0.8)
#' @param comparator	Binary comparator (as string) for ld cutoff, used with the LD score in the ld file. one of ">","<",">=","<=","==", "!=" (default: ">=")
#' @param flank.genes.left	integer: Expand gene regions by this on the left (default: 1000)
#' @param flank.genes.right	integer: Expand gene regions by this on the right (default: 1000)
#' @param multiple.hits	logical: Allow multiple overlaps (Default: TRUE) 
#' @param scoring.function	character: the name of the scoring function to use. One of: "max", "min", "mean", "median", "product", "sum", "count",  "p.ratio", "snp.ratio.score", "saccone", "sidak", "bonferroni", "get.snps", "w.fisher", "brown", "forge", "vegas"
#' @param correction.type character: the name of the correction function to use. One of: "Nyholt", "Moskvina", "Gao"
#' @param genome.name	character: optional genome.name
#' @param outfile.path	character: A file path to save the output to. Default "./LDsnpRout.txt" If NULL, no file is generated
#' @param add.call.header	logical: if TRUE a comment header is added to the output file containing the function call parameters
#' @param generate.plink.set	logical: Generate output file in Plik set format (Default: FALSE) If TRUE, scoring.function is ignored
#' @param ld.structure logical:  If TRUE, LD stucture from hdf5 file is loaded. If FALSE, LD structure is estimated based on p-values (Default: FALSE)
#' @param  ...	Additional parameters passed to the scoring function. (e.g. pMax for the p.ratio)
#' @export 
#' @return gene.score RangedData: A hidden object of type RangedData containing the genomic ranges and a score column with the gene scores resulting from the scoring function. The result is writen to a tab separated txt file unless the output file parameter is NULL.
#' @keywords GWAS,  SNP, LD
#' @references Fisher R. (1925), Statistical methods for research workers. London: Oliver and Lloyd.
#' @references Peirce JL, Broman KW, Lu L, Williams RW.  (2007)  A simple method for combining genetic mapping data from multiple crosses and experimental designs. PLoS One. 2(10):e1036.  PMID: 17940600
#' @references Saccone SF, Hinrichs AL, Saccone NL, Chase GA, Konvicka K, et al. (2007) Cholinergic nicotinic receptor genes implicated in a nicotine dependence association study targeting 348 candidate genes with 3713 SNPs. Hum Mol Genet 16(1).  PMID: 17135278
#' @exportClass H5IdComponent
#' @exportMethod show
#' @details This function wraps the work-flow of a SNP analysis. SNPs are assigned to genes using findOverlaps. Then, gene-specific scores are computed as combined scores from the SNP-pvalues. If LD data is provided (not NULL), the scores of SNPs in high LD are added to  the gene-specific list of scoring genes. Due to the large data-volume, the LD-data must be stored in a hdf5 file in a defined format, grouped by chromosome with 3 or 5 datasets (1-dimensional) per group: snp.id.one(the rs-identifier, first column), snp.id.two (second column), value.set (numeric score value). For full crossmatching the additional data-sets marker.pos.one, marker.pos.two are required. Missing chromosomes are ignored.   Group names in the hdf5 must match chromosome names exactly, otherwise the data will not be loaded. See the example data file for the structure. The rs-identifiers must be stored as integers without the rs prefix.   
#' @import IRanges
#' @import rhdf5
#' @import fastmatch



snp.ld.analysis <- function (snpdata.url = stop("missing snpdata file"),
                             genome.url = stop("missing genome file"),
                             ld.data.hdf.url = NULL,
                             include.ld.data = !is.null( ld.data.hdf.url ),
                             population="CEU",
                             full.match=TRUE,
                             ld.rho.cutoff=0.8,
                             comparator=c(">=","<=",">","<","==", "!="),
                             flank.genes.left=10000,
                             flank.genes.right=10000,
                             multiple.hits=TRUE,
                             
                             scoring.function=c("max",
                               "min",
                               "mean",
                               "median",
                               "product",
                               "sum",
                               "count",
                               "p.ratio",
                               "snp.ratio.score",
                               "fisher",
                               "saccone",
                               "sidak",
                               "bonferroni",
                               "slowscore",#   "simes", # "vegas",    "forge",
                               "get.snps",
                               "get.scores"
                               ),
                             correction.type=c( "Moskvina",
                                                "none"),
                             ld.structure = FALSE,
                             use.position=(scoring.function != "get.snps"),
                             genome.name="unknown genome",
                             outfile.path="./LDsnpRout.txt",
                             add.call.header = TRUE,
                             generate.plink.set=FALSE,
                             ...
                             ){
 
  population <- match.arg(population, c("ASW",
                                        "CEU",
                                        "CHB",
                                        "CHD",
                                        "EUR",
                                        "GIH",
                                        "JPT",
                                        "LWK",
                                        "MEX",
                                        "MKK",
                                        "TSI",
                                        "YRI"))
 
  scoring.function <- match.arg(scoring.function)
  correction.type<-match.arg(correction.type)
  comparator <-  match.arg(comparator)
  comparator <- getMethod(comparator)
  flank.genes.right <- as.integer(flank.genes.right)
  flank.genes.left <- as.integer(flank.genes.left)

  if (! include.ld.data )
    ld.data.hdf.url <- NULL

  cat("loading snps\n")
  snps <-  tab2RangedData(snpdata.url, genome.name=genome.name, name.col=1, start.col=3, chrom.col=2, score.col=4, constwidth=1, header=TRUE, remove.duplicates=FALSE, sep="\t", use.position=use.position)
  cat ("loading genome\n")
  genes <- parse.cds(genome.url, filetype="delim", flank.at.start=flank.genes.left, flank.at.end=flank.genes.right )
  cat("finding Overlaps\n")
  ol <- findOverlaps(query=snps, subject=genes, select="all") # removed obsolete argument multiple =TRUE
  fun <- if (generate.plink.set)
    get.snps2
    else
      switch(scoring.function,
             max = max,
             min = min,
             mean = mean,
             median = median,
             sum = sum,
             product = prod,
             count =  function(x, ...) sum(!is.na(x)),
             p.ratio = ratio.score,
             snp.ratio.score = snp.ratio.score,
             fisher = fisher.combined.p,
             saccone = saccone,
             sidak = sidak,
             slowscore = slowscore,
             bonferroni = bonferroni,            #vegas = vegas,             #forge = forge,
             get.snps = get.snps,
             get.scores = get.scores
             )
  
   
  cat ("computing weights\n")
  
  comb.weights = compute.weights(matched.ranges=ol, snps=snps, genes=genes, ld.hdf.file=ld.data.hdf.url , ld.cutoff=ld.rho.cutoff, population=population, fun=fun, full.match=full.match, comparator=comparator, use.position=use.position, ld.structure = ld.structure, correction.type = correction.type, ...)
 #  reformat the data a bit:
 

  if (! is.null(outfile.path)) {

    
    if (add.call.header)
      writeLines(paste("#", deparse(match.call())), con=outfile.path)

    if (generate.plink.set) {
      write.plink.set(comb.weights, outfile.path, append=add.call.header)
    } else {
        
      comb.weights.frame = as.data.frame(comb.weights)
      if ('gpos' %in% colnames(comb.weights.frame))
        comb.weights.frame = comb.weights.frame[,-grep('gpos',colnames(comb.weights.frame))]
      print(paste("writing output to ", outfile.path))
      
      suppressWarnings(write.table(file=outfile.path, comb.weights.frame, quote=FALSE, sep="\t", append=add.call.header))
    }
  }
  print ("done")
  attr(comb.weights, "call") <- match.call()
  invisible(comb.weights)
}



# dependencies: IRanges, rhdf5
#uire(IRanges)
#require(rhdf5)

## get the names

get.snps <- function(x, sep=";", ...) {
  12398274783298775;
  return(paste(x, sep=sep, collapse=sep))
}

get.scores <- function(x, sep=";", ...) {
  
  return(paste(as.numeric(x), sep=sep, collapse=sep))
}

get.snps2 <- function(x, sep="\n", ...) {
  paste(x, sep=sep, collapse=sep)
}

annotate.snps <- function(x, ...) {
  return(list(as.factor(x)))
}


# define some gene scoring functions according to Andrea Christoforou's notes
## (1) Minimum P-value 
## > minP = min(gene1_pvals)  	#already implemented in the package.

## proportion of significant scores
## (2) SNP ratio – the percentage of the SNPs in the bin that meet a user-defined P-value threshold.

ratio.score <- function(pvals, maxP=0.05, ...) {
  sum (pvals <= maxP, na.rm=TRUE) / sum(!is.na(pvals))
}

## (2a) SNP ratio score ^2 - ratio weighted by the number of total snps assigned 

snp.ratio.score <-  function(pvals, maxP=0.05, ...) {
  rs <- ratio.score(pvals, maxP=0.05, ...)
  len <- sum(!is.na(pvals))
  rs^2 * len  
}
## (3) Fisher’s combination test P-value

## > fishers_stat=(-2)*(sum(log(gene1_pvals),na.rm=TRUE))
## > fishers_pval=pchisq(fishers_stat,2*(length(gene1_pvals)),lower.tail=FALSE)

fisher.combined.p <- function (pvals, ...) {
  fishers.stat= -2 * (sum(log(pvals), na.rm=TRUE))
  
  return (pchisq(fishers.stat, 2*(sum(!is.na(pvals))), lower.tail=FALSE))
}

# Sidak_Saccone2007 = 1-((1-minP)^((N+1)/2))
saccone <- function(pvals, ...) {
  N <- sum(!is.na(pvals))
  minP <- min(pvals, na.rm=T)
  1-((1-minP)^((N+1)/2))
}

sidak <- function(pvals, cor.array=NULL, correction.type = "none", ...){
  
  if (correction.type!="none"){
    fun2 <- switch(correction.type,
                   Nyholt = Nyholt,
                   Moskvina = Moskvina,
                   Gao = Gao
    )
    
    N<-fun2(cor.array)
    
  }else{N <- sum(!is.na(pvals))}
  
  
  minP <- min(pvals, na.rm=T)
  1-(1-minP)^(N)  
}



bonferroni <- function(pvals, ...){
  N <- sum(!is.na(pvals))
  minP <- min(pvals, na.rm=T)
  minP*N
}


simes <- function(pvals,   ...){
  N <- sum(!is.na(pvals))
  r<-rep((1/N), N)
  index<-sort(pvals[!is.na(pvals)], index.return=TRUE)
  N_e<- c(1:N)
  
  min((N*index$x)/N_e)
}


###################################################################################################
###correction types


Moskvina<-function(cor.array)
{
  
  if (is.matrix(cor.array)){
    r_M<-rep(0,nrow(cor.array)-1) 
    k_M<-rep(0,nrow(cor.array)-1)
    
    for (k1 in 1:(nrow(cor.array)-1))
    {
      r_M[k1]<-max(sqrt(cor.array[1:k1,k1+1]))
      k_M[k1]<-sqrt(1-(r_M[k1])^(-1.31*log10(0.01)))
    }
    K_Mosk<-1+sum(k_M)
    
  }else{K_Mosk<-1}
  
  return(K_Mosk)
  
}


#######combination approarch
slowscore<-function(pvals, cor.array, ...)
{
  T_0 <- -2 * (sum(log(pvals), na.rm=TRUE))
  sigma2<-4*length(pvals) + 2*sum(cor.array*(3.25 + 0.75*cor.array))
  c <- sigma2/(4*length(pvals))
  
  pchisq((T_0/c), (2*length(pvals)/c), lower.tail=FALSE)
  
}

write.plink.set <- function(res, file=stdout(), append) {
  df= as.data.frame(res[!is.na(res$'score.data'),])[,c('names','score.data')]
  cat(paste(df$names, df$score.data, "END\n\n", sep="\n"), sep="", file=file, append=append)
  
}


parse.cds <- function (cdsfile, filetype=NULL, flank.at.start=NULL, flank.at.end=NULL) {

if (missing(filetype)) {
    gr <- import(cdsfile)
    return(gr)
  } else if (filetype=="gff") {
    gr <- import.gff(cdsfile, asRangedData=TRUE)
    return (gr)
  } else if(filetype=="delim") {
  
    cdstable = read.delim(cdsfile)
    if (ncol(cdstable) < 4 )
      stop("need at least four columns")
  
    if (! setequal (tolower(colnames(cdstable)[1:4]), c("names","chr", "start", "end")))
       colnames(cdstable)[1:4] <-  c("names","chr", "start", "end")
    genenames <- make.names(cdstable$names, unique=TRUE)
    startstop = t ( apply(cbind(as.integer(cdstable$start), as.integer(cdstable$end)), 1, sort) )
    if (! is.null(flank.at.start))
      startstop[,1] = startstop[,1] - flank.at.start
    if (! is.null(flank.at.end))
      startstop[,2] = startstop[,2] + flank.at.end
    if (any(duplicated(genenames)))
      stop(genenames[(duplicated(genenames))])
    #if (! is.null(flank.at.start))
    ranges = IRanges(start=startstop[,1], end=startstop[,2], names=genenames)
    
    if (ncol(cdstable) > 4 ) {
      meta <- cdstable[,-c(1:4)]
      #row.names(meta) <- genenames
      rd = RangedData(ranges, space=cdstable$"chr", meta)
    } else {
      rd = RangedData(ranges, space=cdstable$"chr")
    }
    
    return(rd)
  } else {
    gr <- import(cdsfile, type=type)
    return(gr)
  }
}



# utility function to load input data into RangedData object, not to be used externally

tab2RangedData <- function(filename, genome.name=NULL, name.col=1, start.col=2, end.col=3, chrom.col=4, score.col=NULL, flank.at.start=NULL, flank.at.end=NULL, constwidth=NULL, remove.duplicates=FALSE, header=TRUE, na.strings=c("", "NA","na"), sep="\t", use.position, ...) {
  
  genes = read.table(filename, sep=sep, na.strings=na.strings, nrows=4, ...)
  what = as.list(rep(NULL, ncol(genes)))
  what[[name.col]] = character(0)
  what[[chrom.col]] = character(0)
  what[[start.col]] = integer(0)
  if (! is.null(end.col) )
    what[[end.col]] = integer(0)
  if (! is.null(score.col) )
    what[[score.col]] = double(0)
  
  genes = scan (filename, what=what, sep=sep, na.strings=na.strings, skip=(as.numeric(header) ))
  ## remove possibly incomplete records
 ## genes = genes[complete.cases(genes),]
  if (remove.duplicates) {
    # TODO: not neccessary
    
  #  genes = genes[! duplicated(genes[[name.col]])]
  }
  
  ##  if (any(duplicated(genes[[name.col]]))) {
  ##    stop (genes[duplicated(genes[[name.col]])])
  ##  }

  
  ## make new marker names based on possition: like [chr.pos] e.g. 22.55100
  ## or keep the marker names as they are
  range.names = if (use.position)
    as.character(paste(genes[[chrom.col]], genes[[start.col]], sep=".")) 
  else
    make.names(as.character(genes[[name.col]]))
  print (head(range.names))  
  if (is.null(constwidth)) {
    stpos = genes[[start.col]]
    endpos = genes[[end.col]]
    if (! is.null(flank.at.start))
      stpos = stpos - flank.at.start
    if (! is.null(flank.at.end))
      endpos = endpos + flank.at.end
    ir =  IRanges(start=stpos, end=endpos, names=range.names)
    
    genome = RangedData(ranges=ir, space=genes[[chrom.col]], universe=genome.name)
  }
  else {
    genome <- if (is.null (score.col)) {
      RangedData(ranges =
                 IRanges(start=genes[[start.col]], width=constwidth, names=range.names),
                 space=genes[[chrom.col]], universe=genome.name) }
    else {
      RangedData(ranges =
                 IRanges(start=genes[[start.col]], width=constwidth, names=range.names),
                 space=genes[[chrom.col]], score=genes[[score.col]], universe=genome.name)
    }
  }
  
  #rownames(genome) <-  genome[,"names"]
  return (genome)
}

# simulate some SNP data if you don't have them

simulate.snpdata <-  function (n, genome.size, nchr=1){

  ir = IRanges(start=as.integer(runif(n, min=1, max=genome.size)), width=1)
  rd = RangedData(ir , p.value = pt(rt(n,df=100), df=100), space=rep(paste( 1:nchr, sep=""), length=n ), universe="Mygenome")
  row

  rownames(rd) = paste ("rs",sprintf(paste("%0",ceiling(log10(n))+1,".0f", sep=""),1:n), sep="")
  return(rd)

}



# the LD data for a whole genome are too big to load into memory
# the LD data is loaded from a hdf5 file with a certain layout:

# top level groups - the population code
# the chromosomes must be saved on the second  as groups
# the group name must match the chromosome name
# on the third level, there is a lists with elements: "rs1" "rs2" "score"
## example /CEU/1/rs1

# cutoff: cutoff value for the LD list of score. The list is reduced immediately
# to shrink it as soon as possible to save memory
# the content is returned as a list with elements $rs1 $rs2 $score
# if full.match == TRUE then 'pos1' and 'pos2' are added containing
# the snp positions in the LD files
# if sim=TRUE, simulated snp data of max. length max.pairs is returned




load.ld <- function(chrom.name=NULL, population=population, snps=NULL, hdf5.file=NULL, max.pairs=1E5, cutoff=0.8, comparator=get(">="), rsprefix="rs", sim=FALSE, full.match=FALSE, use.position) {
  if (sim) {
    # generate a random snp list of size max.pairs
    nams = rownames(snps)
    max.pairs = min(max.pairs, length(nams))
    if (is.null(nams)) stop ("need all rownames to simulate ld data")
    print (paste("sampling random ld space data of size", max.pairs))
    return (list(rs1= sample(nams, size=max.pairs),  rs2= sample(nams, size=max.pairs) ))    
  } 
                                        # this is supported now

  if (is.na(hdf5.file)) {
    print("NA hdf5 file")
    
    return (NULL)
  }
  
  dsnams = as.vector(h5ls(hdf5.file)[[1]])       ##hdf5dir(file=hdf5.file, paste("/", population, sep=""))
  if (!(paste("/", population,"/", chrom.name, sep="") %in% dsnams) ) {
    cat(file=stderr(),paste("Chromosome", chrom.name, "in population", population, "not matched\n"))
    return(NULL)
  }
  
  chr <- chrom.name # save the original chrom name
  chrom.name <-   paste("/", population, "/", chrom.name , sep="")
  # get the data with low memory profile
  print ( paste("loading ld data for ", chrom.name))
  tmp0 = h5read(file=hdf5.file,  chrom.name)
  global.tmp <<- tmp0
  global.ls.slice.name <<- chrom.name
  
  ret = tmp0$value.set  ##(hdf5load(file=hdf5.file, dataset="value.set", group=chrom.name,load=FALSE, verb=2))
  #print ("gc()")
  #gc() # free the mem as soon as possible
                                        # make the filter:
  print ("make filter")
  fint = ( comparator(ret[[1]], cutoff) )
  stopifnot ( is.logical(fint) )
  print (paste("got", sum(fint), "entries"))
  print("loading the data")                                     # filter the list:

  if (!use.position) {
  
    ret = list(
      rs1= tmp0[[3]], ###(hdf5load(file=hdf5.file, dataset="snp.id.one", group=chrom.name,load=FALSE, verb=2))[[1]][fint], 
      rs2= tmp0[[4]] ###(hdf5load(file=hdf5.file, dataset="snp.id.two", group=chrom.name,load=FALSE, verb=2))[[1]][fint]
      )
  
  #gc()
    ## add the prefix again, not so nice
    ret[[1]] = paste(rsprefix, ret[[1]], sep="")
    ret[[2]] = paste(rsprefix, ret[[2]], sep="")
  } else {
    ## make marker names based on possition: like [chr.pos] e.g. 22.55100
    ret = list(
      rs1= tmp0[[1]], ###       (hdf5load(file=hdf5.file, dataset="marker.pos.one", group=chrom.name,load=FALSE, verb=2)[[1]][fint]),
      rs2= tmp0[[2]] ### (hdf5load(file=hdf5.file, dataset="marker.pos.two", group=chrom.name,load=FALSE, verb=2))[[1]][fint]
      )
    ret[[1]] = paste(chr, ret[[1]], sep=".")
    ret[[2]] = paste(chr, ret[[2]], sep=".")
    
    stopifnot(all(!is.na(ret)))
    
  }


  
  if (full.match) {
                                        # load the other contents of the data    
    ret = append (ret, list(
      loc1=tmp0[[1]][fint], ### (hdf5load(file=hdf5.file, dataset="marker.pos.one", group=chrom.name,load=FALSE, verb=2))[[1]][fint], 
      loc2=tmp0[[2]][fint] ###(hdf5load(file=hdf5.file, dataset="marker.pos.two", group=chrom.name,load=FALSE, verb=2))[[1]][fint]
      )
      )
  }
  rm(tmp0)
  
  return(ret)
}


match.inbound <- function(snp.pos, names, snpnames, filter, genes, space) {
          
          # match the =right= column of snp pos, using the names from the left
          # why? because:
          # snppos.1 snppos.2 snpid.1 snpid.2
          # 10       1000     rs1(on chip)     rs2(not on chip)
          # rs1 ~ rs2, but rs2 is not on the chip, take the name of rs2 to get its p-value
  # snps must be unique but the names are not, because they replace snps with
  # on in linkage disequilibrium

  # if there is nothing match then we can get off here
  if (length(filter) < 1 || sum(filter) < 1)
    return ( NULL )
    
  # we to keep the snps in line in the names vector
  names <- names[filter]
  # for reasons of debugging make unique names
  snpsadd = RangedData(
    IRanges(start=snp.pos[filter], width=1, names=make.names(names, unique=TRUE) ), 
    space=space,
    nu.names=names
    )
  #print (snpsadd)
          # the scoring function expects matches per chromosome, even if there is only
          # one space in the matched list 
          # the indices of the genes should be ok then
  ol1 <- as.matrix (findOverlaps(query=snpsadd, subject=genes, select="all")[[space]])
 
          # add these matches to the initial match matrix,
          # therefore, we have to match the snpnames to the positions in the
          # chip snps, the genes are the same as before so we can leave them

          # 1. get the names of the snps:
  rn <- names[ol1[,"queryHits"]] 
          # 2. match the names with the global snps
          # and replace them, there cannot be any missing because we tested for them before 
  ol1[,"queryHits"] <- match(rn, snpnames)
  # ol1 now must be made unique again
  ol1 <- unique(ol1)
 cat(paste("number of inbound overlaps unique:",nrow(ol1),"\n"))
  return (ol1)
}




# compute the combined weights,
# internal function

compute.weights <- function(matched.ranges, snps, genes, ld.hdf.file=NULL, population=NULL, ld.cutoff=0.8, fun, full.match=FALSE, comparator = get(">="), use.position, ld.structure = ld.structure, correction.type = correction.type, ...) {
  # the datasets in the ld data list must have the same names as the spaces
  if (! is.null(ld.hdf.file)) {
    # add a global position vector
    genes$"gpos" <-  1:nrow(genes)
    gscore = rep(NA, nrow(genes))
    
    for (m in 1:length(matched.ranges)) {
      # cannot load all ld data at once, so we have to go space by space
      
      space = names(matched.ranges[m])
      print (paste("processing space",space))
      if (! space %in% names(genes)) {
        cat (paste("space", space, "not present in gene space, skipping\n"))
        next
      }
      snp.space = snps[space]
      gene.space = genes[space]
      matchmat = as.matrix(matched.ranges[[m]])
      print(paste("got", dim(matchmat), "matches"))
      snpnames = as.character(rownames(snp.space)) # a bit faster not to load the names all the time
      if (is.null(snpnames))
        stop ("All SNPs need valid rs-snp ids to run LD analysis")
      
      ldmatch = load.ld(chrom.name=space, population=population, snp.space, hdf5.file=ld.hdf.file, max.pairs=1E6, cutoff=ld.cutoff, sim=FALSE, full.match=full.match, comparator=comparator, use.position=use.position )        
     # ldmatch=NULL
      if (!is.null(ldmatch)) {
        # something is found
        # exchange the identifiers of the ld data with the numeric match indices
        print ("matching ids on chip:")
        print (length(ldmatch[[1]]))
        # print(head(ldmatch[[1]]))
        # print(head(ldmatch[[2]]))
        print(head(snpnames))
        rs1 = match(as.character(ldmatch[[1]]), snpnames)
        rs2 = match(as.character(ldmatch[[2]]), snpnames)
        # print(cbind(rs1,rs2))
        # remove NA from non-matching snps,
        # they have to be matched back by crossmatching
        tmp =  ! (is.na(rs1) | is.na(rs2))
         if ( sum((tmp)) <= 0 ) warning("Nothing matched! Incompatible SNP coordinates?")
        rs1 = rs1[tmp]
        rs2 = rs2[tmp]
       
        print("reduced to")
        print(length(rs1))
                                        # boil down:
                                        # remove those, which are not in "query"
        print("boiling down to:")
        query = as.numeric(matchmat[,"queryHits"])
        ind =  ( ( rs1 %in% query ) | (rs2 %in% query) )
        print (length(ind))
        rs1 = as.numeric(rs1[ind])
        rs2 = as.numeric(rs2[ind])
     
        print("making unique pairs...")
        crossmat = unique( rbind (cbind(rs1,rs2), cbind(rs2,rs1) ) )
        print (dim(crossmat))
        print("added value...")
        ## make two crossmatching lists for fast query:
        addedv = split(f=crossmat[,1], crossmat[,2])
        print("value, not NA:")
        ## extend the match matrix: this needs some dirty tricks:
        addedv2 = addedv[match(as.character(query), names(addedv))]      
        # need to remove data later
        fint = which( ! is.na(names(addedv2)))
        print(length(fint))
       # list has length of 'query' and the same order as the matchmatrix
       # with some missing values, but the
       # names are still wrong, so replace them with the gene names 
        genenams =  as.character(matchmat[,"subjectHits"])    
        names(addedv2) <- genenams
        

       # now we can remove the empty entries,
       # the list now associates genes and added snps 
       addedv3 = addedv2[fint]
        print (paste("length of added value: ", length(addedv3)))
      # the list names are geneids and not unique
       print ("double split") 
      spl2 = split(addedv3, names(addedv3)) # split on listnames to get these into the list, BAD TRICK!
      # make a match gene matrix by repeating the gene
       print ("extend match matrix") 
        ret = lapply(spl2, function(x) {
          as.matrix( cbind(unlist(x), rep(names(x)[1], length(unlist(x))) ) )
        }
          )
        # add crossmatching, therefore we need to overlap backwards again
        # This step is equivalent to simple crossmatch, if and only if: 
        # All LD snps in the LD data are also on the chip
        if (full.match) {
          cat("doing inbound matching\n")
          # we need both the matches to filter
          match.back.1 <- ldmatch[[1]] %in% snpnames 
          match.back.2 <- ldmatch[[2]] %in% snpnames
          # we seek the pairs where exactly one partner is on the chip,
          # if both are on the chip, we have already matched outbound:
          # Done: if only inbound matching is applied, then
          # crossmatch has to be OR because otherwise we would loose these cases
          # Done: implement inbound/outbound strategy separatly        
          crossmatch = xor(match.back.1, match.back.2)
          # this gives us a reduced overlap index
          # reduce the data further:
          match.back.1 <- match.back.1 & crossmatch
          match.back.2 <- match.back.2 & crossmatch
          ol1 <-  match.inbound(ldmatch[[4]], ldmatch[[1]], snpnames, match.back.1, genes, space)
          ol2 <-  match.inbound(ldmatch[[3]], ldmatch[[2]], snpnames, match.back.2, genes, space) 
          # this can be added to the match.matrix:
          matchmat <- rbind(matchmat,ol1, ol2)
          
          
        }

        print ("extend step2")
     #lapply(ret, function(x) { #matchmat <<- rbind(matchmat, x);})
        # add all results to the match matrix
        
        matchmat <-  rbind(matchmat, do.call(rbind, ret))
        print(paste("after extension got", nrow(matchmat), "matches"))
        #
        matchmat <- unique(matchmat)
        print (paste("after uniquifying got", nrow(matchmat), "matches"))
      }
      print("add to scoring")

      gscore =  do.score2(matchmat, genes, gene.space, snp.space, fun, gscore, correction.type, hdf5.file=ld.hdf.file, ld.structure, ...)

    }
 
    # assign score vector to genome:
    genes$"score.data" <- unlist(gscore)
    return(genes)
    
  } else {
   # print(matched.ranges)
    tmp = as.matrix(matched.ranges)   
   # print(dim(tmp))
    return (do.score(tmp, genes, snps, fun, correction.type, hdf5.file=ld.hdf.file, ld.structure, ...))

  }
} 


# rbinds a list of matrices together, all matrices must have
# the same number of columns
# this is faster than using an intermediate dataframe

# a simple test:
# bla = list(a=(matrix(1:55,nrow=1E4,ncol=2,by=T)), b=(matrix(1:7,nrow=1E4,ncol=2,by=T)))  
# all (list.rbind(bla) == rbind(bla[[1]],bla[[2]]))
#should be true

list.rbind <- function(matlist) {
  ncol = ncol(matlist[[1]])
  
  matlist = lapply(matlist,t)
 return ( t(matrix(unlist(matlist), nrow=ncol, byrow=F)) )
}



# compute combined scores given a whole match matrix object
# for the non-LD case
         
do.score <- function (match.mat, genes, snps, fun, correction.type, hdf5.file, ld.structure, ...) {
  ## claculate the combined score
  
  ## get the p-values
  print ("get scores....")
  scorevect = rep(NA, nrow(genes))
  
  match.mat= as.matrix(match.mat)
  
  nu.nams = rownames(snps)
  ##add chr and bp 
  chr.nams = as.numeric(as.character(snps$space))
  bp.nams = data.frame(snps$ranges)$start
  
  #tmp = data.frame(match.mat, score=score( snps[as.numeric(match.mat[,"queryHits"]),]), snp.ids=( rownames( snps[as.numeric(match.mat[,"queryHits"]),]) ) )
  #####12.11.2014. include chr and bp
  tmp = data.frame(match.mat, score=score( snps[as.numeric(match.mat[,"queryHits"]),]),
                   snp.ids = nu.nams[as.numeric(match.mat[,"queryHits"])], snp.chr = chr.nams[as.numeric(match.mat[,"queryHits"])], snp.bp = bp.nams[as.numeric(match.mat[,"queryHits"])], row.names=NULL )
  
  
  
  ## split p-values by gene (subject)
  print("split...")
  if ( ( ! is.null (body(fun)))
      && ( length (body(fun)) == length(body(get.snps)) )
      && all(body(fun) == body(get.snps)))
    {
      print ("using snp names instead of score")
      splitti = split(tmp$snp.ids, f=tmp$subjectHits)
    } else {
      
      splitti = split(tmp$score, f=tmp$subjectHits)
      splitti.chr = split(tmp$snp.chr, f=tmp$subjectHits)
      splitti.bp = split(tmp$snp.bp, f=tmp$subjectHits)
    }
  
  rm (tmp)
  ## apply the combination function
   print("before sapply")
  
  if (is.na(pmatch("cor.array",formalArgs(fun))))
  { 
    tmp2 = lapply(splitti, fun, ...)
    
  }
  else{
    ##splitti.chr and splitti.bp to function func.r2
    ##create r2 matrix and merge list of pvals with the list of matrixes
    ##
    ## print("merge coordinates...")
    coor<-cbind(splitti.chr, splitti.bp)
    
   cor.array<-list()
   
   if (ld.structure){
     cor.array<-mapply(func_r2, splitti.chr, splitti.bp, hdf5.file)
   }else{
     cor.array<-mapply(func_r2_est, splitti, splitti.bp)##mapply(func_r2, splitti.chr, splitti.bp)
     ##i use this option because the exact ld loading is very long. 
     ##
   }
     
   tmp2<-mapply(fun, splitti, cor.array, correction.type)#list()
      
    
  }
  
  
  
  rm(splitti)
  print("getting index")
  ## assign it to the right positions in the scorevector
  ind = as.numeric(names(tmp2))
  print("add to scorevect")
  scorevect[ind] <- tmp2
  
  stopifnot(length(ind) == length(tmp2))
  
  genes$"score.data" <- unlist(scorevect)
  ## the scorevector contains all gene scores  and can be attached to the genes
  return (genes)

}



## do the combined score and return a score vector
## applied for each chromosome/space and a match matrix made
## from ranges lists
## the gene space needs a gpos field denoting the absolute position
## of the Range
## takes and returns a scorevector
         
do.score2 <- function (match.mat, global.genes, local.genes, local.snps, fun, score.vector, correction.type, hdf5.file=ld.hdf.file, ld.structure, ld.hdf.file,...) {
                                        # claculate the combined score
  
                                        # get the p-values
  print ("get scores 2....")
  stopifnot( nrow(global.genes) == length(score.vector))
  if (any (is.na(local.genes$gpos)) || any (is.null(local.genes$gpos)))
    stop("all ranges need a global position attribute(gpos)")
  
             

  match.mat = as.matrix(match.mat)
  nu.nams = rownames(local.snps)
  ##add chr and bp 
  chr.nams = as.numeric(as.character(local.snps$space))
  bp.nams = data.frame(local.snps$ranges)$start
  
  
  #tmp = data.frame(match.mat, score=score( local.snps[as.numeric(match.mat[,"queryHits"]),]),
  #  snp.ids = nu.nams[as.numeric(match.mat[,"queryHits"])], row.names=NULL )
  #                                      # split p-values by gene (subject)

  #####12.11.2014. include chr and bp
  tmp = data.frame(match.mat, score=score( local.snps[as.numeric(match.mat[,"queryHits"]),]),
                   snp.ids = nu.nams[as.numeric(match.mat[,"queryHits"])], snp.chr = chr.nams[as.numeric(match.mat[,"queryHits"])], snp.bp = bp.nams[as.numeric(match.mat[,"queryHits"])], row.names=NULL )
  # split p-values by gene (subject)
  
  print("split...")
 if ( ( ! is.null (body(fun)))
     && ( length (body(fun)) == length(body(get.snps)) )
     && all(body(fun) == body(get.snps))) {
    print ("using snp names instead of score")
    
    splitti = split(tmp$snp.ids, f=tmp$subjectHits)
  } else {
    
    splitti = split(tmp$score, f=tmp$subjectHits)
    splitti.chr = split(tmp$snp.chr, f=tmp$subjectHits)
    splitti.bp = split(tmp$snp.bp, f=tmp$subjectHits)
    
  }

    
 
  #rm (tmp)
                                        # apply the combination function
 # print("before sapply")
   
 if (is.na(pmatch("cor.array",formalArgs(fun))))
   {
     tmp2 = lapply(splitti, fun, ...)
   } else{
   ##splitti.chr and splitti.bp to function func.r2
   ##create r2 matrix and merge list of pvals with the list of matrixes
   ##
   coor<-cbind(splitti.chr, splitti.bp)
   
   if (ld.structure){
     print ('retreiving ld structure for all genes')
     ### SLOW, this is the point that needs to be improved
     print(system.time(cor.array<-mapply(func_r2, splitti.chr, splitti.bp, hdf5.file)))
     }else{
     cor.array<-mapply(func_r2_est, splitti, splitti.bp)##mapply(func_r2, splitti.chr, splitti.bp)
     }
                                                                            
   tmp2<-mapply(fun, splitti, cor.array, correction.type)
  
   
    
 }
  
  #rm(splitti)
  print("getting index")
    # assign it to the right positions in the scorevector
  ind = as.integer(names(tmp2))
  # the index is the names of the local.genes in chromosome
  # as split might have reordered completely, we have to re-order it
  global.ind = local.genes[ind,]$gpos
  if (length(global.ind) != length(tmp2))
    stop (paste("different length in scores: gind: ",length(global.ind),"tmp ", length(tmp2) ))
  
  print("add to scorevect")
  score.vector[global.ind] <- tmp2
 # genes$"score.data" <- scorevect
     
  # the scorevector contains all gene scores  and can be attached to the genes
  return (score.vector)

}


#######################################################################################################################################
func_r2<-function(chr, bp, hdf5.file){
  chr<-as.numeric(chr)
  bp<-as.numeric(bp)
  
  if (length(bp)>1){
    print(bp[1])
      mat<-load.ld.matrix.2(hdf5.file=hdf5.file,chrom.name=chr[1], positions=bp)
     # index<-match(bp, row.names(matrix)) This is already done by the load ld matrix function
     # cor.array<-matrix[index,index]
  }else{
    mat<-matrix(1)
  }
  
  return(mat)
}




#####################################################################################################


load.ld.matrix.bak <- function(chrom.name=NULL, population="EUR", hdf5.file=NA, positions=NA, verbose=F, rsprefix="rs") {
  
  positions <- sort(unique(positions)) # just in case
  print("tmp1...")
  if (is.na(hdf5.file)) {  print("NA hdf5 file"); return (NULL) } 
  #if (anyNA(positions)) { print("NA positions"); return(NULL) }
  
  dsnams = as.vector(h5ls(hdf5.file)[[1]])       ##hdf5dir(file=hdf5.file, paste("/", population, sep=""))
  if (!(paste("/", population,"/", chrom.name, sep="") %in% dsnams) ) {
    cat(file=stderr(),paste("Chromosome", chrom.name, "in population", population, "not matched\n"))
    return(NULL)
  }
  
  ld.slice.name <- paste("/", population, "/", chrom.name , sep="")

  ## try to avoid loading same ld data over and over again
  if (ld.slice.name != global.ls.slice.name) {
    rm (global.tmp)
    rm (global.ls.slice.name)
    global.ld.slice.name <<- ld.slice.name
    ## get the data with low memory profile
    print ( paste("Loading ld data positions for", ld.slice.name))
    tmp <- h5read(file=hdf5.file, ld.slice.name)
    global.tmp <<- tmp
  } else {
    print ("Not loading, using global object!")
    tmp <- global.tmp
  }

    
    snp.pairs = list(
      pos1 = tmp[[1]],
      pos2 = tmp[[2]]
      )
  
  print ( paste("Selecting associations for",length(positions),'snps on chromosome',chrom.name) )
  ### SLOW:
  ### making the selection is slow
  selection <- (( snp.pairs[[1]] %in% positions & snp.pairs[[2]] %in% positions ))
  ### if selection is empty, we can return immediately
  if (sum (selection) < 1) {
    return(diag(1,nrow=length(positions), ncol=length(positions) ))
  }
  ### SLOW: change this into some vectorized operation
  print (" make matrix of indexes, this is slow" )
  pos1.sel = snp.pairs[[1]][selection]
  pos2.sel = snp.pairs[[2]][selection]
  #### This is better written as: 
  #m<-matrix(rep(NA,length(positions)**2),ncol=length(positions))
  m<-matrix(NA, nrow = length(positions), ncol=length(positions))
  
  for (i in 1:length(positions)) {

    ### SLOW: if the matrix is symmetric, we don't need to start at
    ### j for each i, instead we could start at i+1 ???? 
    for (j in 1:length(positions)) {
      if (verbose) {
        #cat(length(new.row),'.',fill=F,sep='')
        cat('.',fill=F,sep='')
      }
      ### SLOW: this is the diagonale,isnt it set to 1 anyway?
      ### This will no longer happen and the matrix was already initialized as such
      if (i == j) { 
        m[i,j] <- NA
        next
      } 
      
      if ( sum(positions[i] == pos1.sel & positions[j] == pos2.sel) >= 1 ) {
        if (sum(positions[i] == pos1.sel & positions[j] == pos2.sel) > 1) { print("Multiple associations!!!") 
                                                                            tmp0 <- which(selection)[which(pos1.sel == positions[i] & pos2.sel == positions[j])]
                                                                            m[i,j]<-tmp0[1]}else{
                                                                              m[i,j] <- which(selection)[which(pos1.sel == positions[i] & pos2.sel == positions[j])]}
        m[j,i]<-m[i,j]
      }
    }
  }
  print ("...done")
  #cross matrix.selection with assoc
  assoc = tmp$value.set[as.vector(m)]
  #  rm(tmp)
  if (length(pos1.sel)>0){
    ret = matrix(assoc, nrow=length(positions), ncol=length(positions))  
  }else{ret<-mat.or.vec(length(positions),length(positions))}
  # make result matrix
  
  ret[is.na(ret)]<-0
  
  rownames(ret)<-positions
  colnames(ret)<-positions
  # set the diagonal
  diag(ret) <- 1.0 
  
  return(ret)
}

load.ld.matrix.2 <- function(chrom.name=NULL, population="EUR", hdf5.file=NA, positions=NA, verbose=F, rsprefix="rs") {
  
  positions <- sort(unique(as.integer(positions)))
   mti <- rep(positions, 1000) # don't ask :D
  # mti <- positions
                                        # just in case, we don't need anything twice
  #print("tmp1...")
   if (is.na(hdf5.file)) {  print("NA hdf5 file"); return (NULL) } 
  #if (anyNA(positions)) { print("NA positions"); return(NULL) }
  
  dsnams = as.vector(h5ls(hdf5.file)[[1]])       ##hdf5dir(file=hdf5.file, paste("/", population, sep=""))
  if (!(paste("/", population,"/", chrom.name, sep="") %in% dsnams) ) {
    cat(file=stderr(),paste("Chromosome", chrom.name, "in population", population, "not matched\n"))
    return(NULL)
  }
  
  ld.slice.name <- paste("/", population, "/", chrom.name , sep="")

  ## try to avoid loading same ld data over and over again
  if (ld.slice.name != global.ls.slice.name) {
    rm (global.tmp)
    rm (global.ls.slice.name)
    global.ld.slice.name <<- ld.slice.name
    ## get the data with low memory profile
    print ( paste("Loading ld data positions for", ld.slice.name))
    tmp <- h5read(file=hdf5.file, ld.slice.name)
    global.tmp <<- tmp
  } else {
    print ("Not loading, using global object!")
    tmp <- global.tmp
  }
  
  ## redundant assignment  
  #snp.pairs <- NULL
  ## check what the fastest processing would be
  
  print ( paste("Selecting associations for",length(positions),'snps on chromosome',chrom.name) )
  ### SLOW:
  ### making the selection is slow, and incorrect, because we wanted all snps contained in genes, right?
  ### select all snps where at least one is inside the gene:
  # selection <- ( tmp[[1]] %in% positions | tmp[[2]] %in% positions )
### select all snps where both snps are inside the gene:
  ## leaving it like this for timing comparison
  selection <- intersect(
                         which(fmatch(tmp[[1]], mti, nomatch = 0 ) != 0),
                         which(fmatch(tmp[[2]], mti, nomatch = 0 ) != 0)
  )
  
  #browser()
  

  #### selection <- ( tmp[[1]] %in% positions) & (tmp[[2]] %in% positions )
  ### if selection is empty, we can return immediately,
  ### Question, what should be return then????
  
  ### Make the return matrix, it's a diagonale matrix with all other entries 0
  ### dimension is the same as length positions 
  mat <- diag(1,nrow=length(positions), ncol=length(positions) )
  ### assing proper dimnames
  dimnames(mat) <- list(positions, positions)
  
  if (length (selection) < 1) {
    return(mat)
  }
  ### SLOW: change this into some vectorized operation
  print ("reduce to selected entries only" )
  pos1.sel <- tmp$marker.pos.one[selection]
  pos2.sel <- tmp$marker.pos.two[selection]
  val <- tmp$value.set[selection] # the r2 values


  
  #### This is better written as: 
  #m<-matrix(rep(NA,length(positions)**2),ncol=length(positions))
  #m<-matrix(NA, nrow = length(positions), ncol=length(positions))

  ### assign sparse values to matrix by indexing using a matrix:
  ### create the match index:
  print("matching positions")
  #ind1 <- match(pos1.sel, positions)
  #ind2 <- match(pos2.sel, positions)

  mat[cbind(as.character(pos1.sel), as.character(pos2.sel))] <- val
  return (mat)
  # print ("computing eigen values")
  ### that's it, compute the eigen values and return them only:
  ## does not work 
  #  return (eigen(mat, only.values=T))
  #### END!


}






############################
###dummy temporal function 
func_r2_est<-function(p,bp){
  p<-as.numeric(p)
  bp<-as.numeric(bp)
  
   
  if (length(p)>1){
  c<-mat.or.vec(length(p),length(p))
  diag(c)<-1
  for (i in 1:(length(p)-1)){
    for (j in (i+1):(length(p)))
    {
      if (abs(bp[i]-bp[j])<500000){
      c[i,j]<-min(p[i],p[j])/max(p[i],p[j])
      c[j,i]<-c[i,j]
      }
    }
  }
  }else{c<-1}
  
  return(c)
}


