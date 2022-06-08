##################################################
## Project: sc-STITCHIT
## Script purpose: Gene Expression Discretization
## Author: Laura Rumpf
##################################################


suppressPackageStartupMessages(require(mclust))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(EDASeq))
suppressPackageStartupMessages(require(biomaRt))
suppressPackageStartupMessages(require(DESeq2))

# input
read_input <- function(file){

  file_name <- unlist(file)

  if(strsplit(file_name, "\\.")[[1]][2] == "csv"){
    read.csv(file, header = TRUE,  check.names=FALSE)
  } else {
    read.table(file, header = TRUE,  check.names=FALSE)
  }


}

# Remove duplicate genes and set gene ids as row names
remove_duplicates <- function(df){

  if (is.character(df[1,1])){

    #print("gene names in first column")
    #print(df[1,1])
    # filter gene id duplicates in count matrix: upper-/lowercase
    df[,1] <- toupper(df[,1])
    df <- df[df[,1] %in% names(which(table(df[,1]) < 2)), ]

    rownames(df) <- df[,1]
    df <- df[,-1]

  }else{

    #print("gene names as rownames")
    df[,ncol(df)+1] <- toupper(rownames(df))
    df <- df[df[,ncol(df)] %in% names(which(table(df[,ncol(df)]) < 2)), ]
    rownames(df) <- toupper(rownames(df))
    df <- df[,-(ncol(df))]



  }


}


#' Gene lengths query function
#'
#' @description
#' 'get_genelengths' takes an ASCII input file containing a count matrix of expression
#' values with ensemble gene ids as rownames
#'
#' @param file (Path) input ASCII file containing count values (Columns=Samples, Rows=Genes)
#' @return output file containing gene length vector that can be used in tpm normalization
#'
#' @export
#'
get_genelengths <- function(file, output, genomeflag){

  print(genomeflag)

  data <- read_input(file)

  # log_file <- file("deleted_genes.log")

  # get path of input file
  target.directory <- dirname(normalizePath(file))

  # create log file for removed genes in translation gene symbols -> ensemble ids
  log_file <- file.path(target.directory, "deleted_genes.log")
  cat("Removed Genes When Mapping Gene Symbols To Ensemble Id\n\n", file = log_file)

  # set path for output file
  if (is.na(output)){
    output.folder <- target.directory
  }

  # get gene names

  data <- remove_duplicates(data)

  gene.names <- rownames(data)

  #print("gene ids: ")
  #print(gene.names)
  #print("data without duplicates: ")
  #print(head(data))

 # human genome
 if (genomeflag == 'h'){

   print("get human gene lengths")

   if(!mapply(grepl,'ensg', gene.names[1], ignore.case=TRUE)){
   # if (!grepl('ens', gene.names[1])){

     ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

     mart_res <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                       filters= "hgnc_symbol", values= gene.names, mart= ensembl)


     # multiple ensemble ids for one gene symbol -> log entry
     mult_ens <- unique(mart_res$hgnc_symbol[duplicated(mart_res$hgnc_symbol)])

     print("multiple ensemble ids")
     print(mult_ens)

     for (entry in mult_ens){
       cat(entry, ": more than one ensemble id for gene symbol\n", file = log_file, append = TRUE)

     }

     # multiple gene symbols get one ensemble id -> log entry -> obsolete?
     mult_hgnc <- unique(mart_res$hgnc_symbol[duplicated(mart_res$ensembl_gene_id)])


     for (entry in mult_hgnc){
       cat(entry, ": same ensemble id for multiple gene symbols\n", file = log_file, append = TRUE)
     }


     #no ensemble id for gene symbol -> log entry
     no_ens <- gene.names[which(!(toupper(gene.names) %in% toupper(mart_res$hgnc_symbol)))]

     # print("gene symols with no ensemble ids")
     # print(no_ens)

     #if (length(no_ens) !=0){
     for (entry in no_ens){
       cat(entry, ": no ensemble id for gene symbol\n", file = log_file, append = TRUE)
       # }
     }

     # keep only 1 to 1 mapping
     reduced_mart <- mart_res[mart_res$ensembl_gene_id %in% names(which(table(mart_res$ensembl_gene_id) < 2)), ]
     # obsolete?
     reduced_mart2 <- reduced_mart[reduced_mart$hgnc_symbol %in% names(which(table(reduced_mart$hgnc_symbol) < 2)),]

     #change rownames in counts to ensembl ids
     data <- data[toupper(rownames(data)) %in% toupper(reduced_mart2$hgnc_symbol), ] #keep only genes with 1 to 1 mapping in count matrix

     #create hashmap ensemble id -> gene symbol

      print("creating hashmap gene symbol -> ensembl id")
      keys <- toupper(reduced_mart2$hgnc_symbol) #reduced_mart2[,1]
      H <- new.env(hash = TRUE, size = nrow(reduced_mart2))
      for (i in 1:length(keys)) {
       x <- toString(keys[i])
       #print(i)
       #print(x)
       H[[x]] <- reduced_mart2[i,1]
     }

       #print("Hash: ")
       #print(ls.str(H))

     for (key in keys){
       id <- toString(key)
       rownames(data)[rownames(data) == id] <- toString(H[[id]])
     }

     info <- getGeneLengthAndGCContent(id=reduced_mart2$ensembl_gene_id , org="hsa")

     # reduced_mart2[reduced_mart2$ensembl_gene_id %in% rownames(info),3] <- info[reduced_mart2$ensembl_gene_id %in% rownames(info),1]

     # ensembl ids necessary for stitchIt
     # lengths <- reduced_mart2[, ! names(reduced_mart2) %in% "ensembl_gene_id", drop = f]
     lengths <- info[,-2,drop=F]

   } else {

     info <- getGeneLengthAndGCContent(id=gene.names , org="hsa")
     lengths <- info[,-2,drop=F]
     #lengths <- data.frame(gene.names)
     #lengths[lengths$gene.names %in% rownames(info),2] <- info[lengths$gene.names %in% rownames(info),1]

   }

 } else if (genomeflag == 'mm') {

    print("get mmusculus gene lengths")

   if(!mapply(grepl,'ensmus', gene.names[1], ignore.case=TRUE)){

      ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

      #get ensemble ids
      mart_res <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"),
                        filters= "mgi_symbol", values= gene.names, mart= ensembl)


      # multiple ensemble ids for one gene symbol -> log entry
      mult_ens <- unique(mart_res$mgi_symbol[duplicated(mart_res$mgi_symbol)])

      print("multiple ensemble ids")
      print(mult_ens)

      #if (length(mult_ens) != 0){
        for (entry in mult_ens){
          cat(entry, ": more than one ensemble id for gene symbol\n", file = log_file, append = TRUE)

        }

      # multiple gene symbols get one ensemble id -> log entry -> obsolete?
      mult_mgi <- mart_res$mgi_symbol[duplicated(mart_res$ensembl_gene_id)]


        for (entry in mult_mgi){
          cat(entry, ": same ensemble id for multiple gene symbols\n", file = log_file, append = TRUE)
        }

      #no ensemble id for gene symbol -> log entry
      no_ens <- gene.names[which(!(toupper(gene.names) %in% toupper(mart_res$mgi_symbol)))]

      # print("gene symols with no ensemble ids")
      # print(no_ens)

      #if (length(no_ens) !=0){
        for (entry in no_ens){
          cat(entry, ": no ensemble id for gene symbol\n", file = log_file, append = TRUE)
       # }
      }

      # keep only 1 to 1 mapping
      reduced_mart <- mart_res[mart_res$ensembl_gene_id %in% names(which(table(mart_res$ensembl_gene_id) < 2)), ]
      reduced_mart2 <- reduced_mart[reduced_mart$mgi_symbol %in% names(which(table(reduced_mart$mgi_symbol) < 2)),]

      # print(head(reduced_mart2))

      #change rownames in counts to ensembl ids
      data <- data[toupper(rownames(data)) %in% toupper(reduced_mart2$mgi_symbol), ] #keep only genes with 1 to 1 mapping in count matrix
      #print(counts)

        #create hashmap ensemble id -> gene symbol
        keys <- toupper(reduced_mart2$mgi_symbol) #reduced_mart2[,1]
        H <- new.env(hash = TRUE, size = nrow(reduced_mart2))
        for (i in 1:length(keys)) {
          x <- toString(keys[i])
          #print(i)
          #print(x)
          H[[x]] <- reduced_mart2[i,1]
        }

        #print("Hash: ")
        #print(ls.str(H))

        for (key in keys){
          id <- toString(key)
          rownames(data)[rownames(data) == id] <- toString(H[[id]])
        }

      info <- getGeneLengthAndGCContent(id=reduced_mart2$ensembl_gene_id , org="mmusculus_gene_ensembl")

      #reduced_mart2[reduced_mart2$ensembl_gene_id %in% rownames(info),3] <- info[reduced_mart2$ensembl_gene_id %in% rownames(info),1]

      # only keep ensembl ids
      #lengths <- reduced_mart2[, ! names(reduced_mart2) %in% "ensembl_gene_id", drop = F]
      lengths <- info[,-2,drop=F]

    } else {

      info <- getGeneLengthAndGCContent(id=gene.names , org="mmusculus_gene_ensembl")
      lengths <- info[,-2,drop=F]
      #lengths <- data.frame(gene.names)
      #lengths[lengths$gene.names %in% rownames(info),2] <- info[lengths$gene.names %in% rownames(info),1]

    }

  } else if (genomeflag != 'h' & genomeflag != 'mm') {
       stop("supported genome options -g for the gene lengths are human 'h' and mouse 'mm'")
  }

    #print("gene lengths: ")
    #print(head(lengths))

    #print("counts: ")
    #print(head(data))

  print("save counts with ensemble gene ids and gene lengths")

  if(is.na(output)){
    write.table(data, file=file.path(output.folder, "ensembl_counts.txt"), sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
    write.table(lengths, file=file.path(output.folder, "gene_lengths.txt"), sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
    gl <- file.path(output.folder, "gene_lengths.txt")
    c <- file.path(output.folder, "ensembl_counts.txt")
    glList <- list("infile" = c, "outfile" = gl)
    return(glList)
  }else{
    out <- dirname(normalizePath(output))
    write.table(data, file=file.path(out, "ensembl_counts.txt"), sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
    write.table(lengths, file=output, row.names=TRUE, sep = "\t", quote = FALSE, col.names=TRUE)
    gl <- output
    c <- file.path(out, "ensembl_counts.txt")
    glList <- list("infile" = c, "outfile" = gl)
    return(glList)
  }



}

#' Normalize function
#'
#' @description
#' 'normalize' takes an ASCII input file containing a count matrix of expression values and a vector of the same length
#' containing the corresponding gene lengths and returns an output file with the DeSeq2 and TPM normalized values
#'
#' @param file (Path) input ASCII file containing count values (Columns=Samples, Rows=Genes)
#' @param lenvector file containing the gene names and corresponding lengths
#' @param output.norm
#' @return output file containing a DESeq and TPM normalized value for every count value of the input file
#'
#' @export
#'
normalize <- function(file, lenvector, output){

  print("normalizing ...")

  # get path of input file
  target.directory <- dirname(normalizePath(file))

  # set path for output file
  if (is.na(output)){
    output.folder <- target.directory
  }

  # read input files

  count.table <- read_input(file)

  # print(head(count.table))

  count.table <- remove_duplicates(count.table)

  # print(head(count.table))

  length.table <- read.table(lenvector, header = TRUE)

  #length.table <- remove_duplicates(length.table)

  #if (is.character(length.table[1,1])){
  #  rownames(length.table) <- toupper(length.table[,1])
  #  length.table <- length.table[,-1] #length vector
  #}else{
  #  rownames(length.table) <- toupper(rownames(length.table))
  #}

  if(!(mapply(grepl,'ens', rownames(count.table), ignore.case=TRUE)
       && mapply(grepl,'ens', rownames(length.table), ignore.case=TRUE))){
    warning("count matrix and length vector have to contain ensembl ids as gene names (rownames) when used for stitchit input!")
  }


  # not transformed ensembl ids: keep only genes with existing length value
  count.table <- count.table[rownames(count.table) %in% rownames(length.table),]

  print("first step: normalize across samples using DEseq")

  ##############################################################################
  # first step: DEseq2 for normalization across samples
  ##############################################################################

  # create dummy meta df as only normalized count matrix is needed
  # rownames in same order as colnames count matrix
  # dummy value: celltypes (-> not used for creating normalized count matrix)
  sample.num <- ncol(count.table)
  meta.df <- data.frame(matrix(NA, nrow = sample.num, ncol = 1))
  rownames(meta.df) <- colnames(count.table)
  colnames(meta.df) <- "celltypes"
  meta.df[,1] <- rownames(meta.df)

  dds <- DESeqDataSetFromMatrix(countData = round(count.table), colData = meta.df, design = ~ celltypes)
  dds <- estimateSizeFactors(dds)
  #sizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)

  counts <- as.matrix(normalized_counts)

  print("second step: tpm normalization")

  ##############################################################################
  #second step: tpm normalization
  ##############################################################################

  # create hash table for gene lengths

  keys <- rownames(length.table)

  H <- new.env(hash = TRUE, size = nrow(length.table))

  vapply(keys, function(x) {
    H[[x]] <- length.table[x,1]/1000 #kilobases
    # print(x)
    # print(length.table[x,2])
    logical(0)
  }, FUN.VALUE = logical(0))

  # print(all.equal(sort(names(H)), keys))
  # print(ls.str(H))


  # check if all genes from input file are in hash table
  g <- row.names(normalized_counts)
  if(!all(g %in% names(H))){
    stop("Every gene in the count matrix must have a specified length")

  }

  rpk <- t(sapply(1:nrow(counts), function(i) {
     counts[i,]/H[[rownames(normalized_counts)[i]]]
   }))

  # normalize for sequencing depth
  tpm <- apply(rpk, MARGIN = 2, function(x) {x/sum(as.numeric(x)) * 10^6})

  #print("tpm values: ")
  #print(head(tpm))

  rownames(tpm) <- rownames(normalized_counts)

  if(is.na(output)){
    write.table(tpm, file=file.path(output.folder, "tpm.txt"), sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
    return(file.path(output.folder, "tpm.txt"))
  }else{
    write.table(tpm, file=output, sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
    return(output)
  }

}

#' Discretize function
#'
#' @description
#' 'discretize' takes an ASCII input file containing a matrix of normalized
#' gene expression values and returns an output file with the discretized values (either classes 0,1 or -1,0,1)
#'
#' @param file (Path) input ASCII file containing normalized gene expression values (Columns=Samples, Rows=Genes)
#' @param figflag Binary parameter for plotting (TRUE = plot output)
#' @param binaryflag Binary parameter for binary discretization (TRUE = 0,1 classes)
#' @param output (Path) name output file
#' @return output file containing a discrete value for every value of the input file (0,1 values for binaryflag = TRUE, else -1,0,1)
#'
#' @export
#'
discretize <- function(file, figflag, binaryflag, output){

  print("discretizing ...")

 # read input file
 data <- read_input(file) #read.table(file, header = TRUE)

 # path settings for output files

   # get path of input file
   target.directory <- dirname(dirname(normalizePath(file)))

   # set path for plot output files
   if (figflag){

    # set path
    plot.folder <- file.path(target.directory, "Plots") # folder name for plot outputs

    # if folder "Sample Plots" does not exist, create folder
    if (!file.exists(plot.folder)){
      dir.create(plot.folder, showWarnings = TRUE)
     }
   }

   # set path for output file
   if (is.na(output)){

    output.folder <- file.path(target.directory, "DGE") # folder name for discretized matrix

    if (!file.exists(output.folder)){
      dir.create(output.folder, showWarnings = TRUE)
     }
   }

   tpm.counts <- as.matrix(data)

   signal <- log2(tpm.counts)  # log-transform the tpm values
   signal[is.infinite(signal)] = -10000  # transform all -inf to -10000

   discretized.keep <- matrix(nrow = nrow(tpm.counts), ncol = ncol(tpm.counts))
   rownames(discretized.keep) <- rownames(data)
   colnames(discretized.keep) <- colnames(data)

 # discretize each sample
 for (j in 1 : dim(signal)[2]){  # for each sample

  # console output: status
  write(paste("Calculation of Sample:", colnames(signal)[j]), "")

  signal.sample <- signal[,j]  # save current sample for density estimation
  data.keep <- signal[,j]  # all values of current sample for discretization
  signal.sample <- signal.sample[signal.sample > -10000]  # remove sample values with low expression

  # signal density estimation
  est <- density(signal.sample, bw = "nrd0", adjust = 1, kernel="gaussian", n = 100)


  if (figflag){

   # plot signal density
   est.df <- data.frame(X = est$x, Y = est$y)

   p <- ggplot2::ggplot(est.df, mapping = aes(x=est$x, y=est$y))
   p <- p + geom_line(aes(y=est$y, color = "signal", linetype = "signal"), size = 1.2)
   p <- p + labs(title=paste("Density Plot for Sample:", colnames(signal)[j]),
                x ="log2(tpm)", y = "density")

  }


   # gaussian mixture model with 2 mixture components

    signalfit <- mclust::densityMclust(signal.sample, G=2)

    # first gaussian component of mixture model
    norm1 <- signalfit$parameters$pro[1]*dnorm(est$x, signalfit$parameters$mean[1], sqrt(signalfit$parameters$variance$sigmasq[1]))

    # Second gaussian component of mixture model
    if (length(signalfit$parameters$variance$sigmasq) == 1){ # both components have same variance -> only one entry
      norm2 <- signalfit$parameters$pro[2]*dnorm(est$x, signalfit$parameters$mean[2], sqrt(signalfit$parameters$variance$sigmasq[1]))
    }else{
      norm2 <- signalfit$parameters$pro[2]*dnorm(est$x, signalfit$parameters$mean[2], sqrt(signalfit$parameters$variance$sigmasq[2]))
    }


   # plot gaussian mixture components
   if(figflag){

     p <- p + geom_line(aes(y=norm1, color = "gauss1", linetype = "gauss1"), size = 1.2)
     p <- p + geom_line(aes(y=norm2,  color = "gauss2", linetype = "gauss2"), size = 1.2)

   }

  # treshold calculation

   # starting point for index search: mu1 - 2*sigma1
    low.idx.gauss1 <- signalfit$parameters$mean[1]-2*sqrt(signalfit$parameters$variance$sigmasq[1])

   # binary index: intersection of two gaussian components
    intersection.idx.bin <- (norm1 < norm2) & (est$x >= low.idx.gauss1)  # y-value first gaussian is lower than second gaussian
   # intersection.idx2 <- (norm1 < est$y) & (est$x >= min(est$x[est$y > 0.001])) #(est$x >= low.idx.gauss1)  # y-value first gaussian is lower than signal-curve

   # maximum signal density estimation curve
    peak.max.y <- which.max(est$y)  # y-value
    peak.max.x <- est$x[peak.max.y]  # x-value
    # peak.max.idx <- c(peak.max.x,peak.max.y)

    if(binaryflag){ # calculation expression treshold

      # intersection two gaussian curves
      if (any(intersection.idx.bin)){
          intersection.point.bin <- min(est$x[intersection.idx.bin]) #leftmost intersection

         # lower bound for intersection point
         if(intersection.point.bin < min(est$x[est$y > 0.001])){
            intersection.point.bin <- peak.max.x
         }
      } else {
         intersection.point.bin <- peak.max.x
      }
    } else { # no binary flag: calculation inexpression and expression treshold

       # means of gaussian mixture components as tresholds
       intersection.point.inexp <- signalfit$parameters$mean[1]
       intersection.point.exp <- signalfit$parameters$mean[2]

       # lower bound inexpression treshold
       if(intersection.point.inexp < min(est$x[est$y > 0.001])){
          intersection.point.inexp <- min(est$x[est$y > 0.001])
       }
    }


    # plot tresholds
    if (figflag){
      if (binaryflag){
        p <- p + geom_line(aes(x =intersection.point.bin, linetype = "treshold", color = "treshold"), size = 1)
      }else{
        p <- p + geom_line(aes(x=intersection.point.inexp, color = "inexpression treshold", linetype = "inexpression treshold"), size = 1)
        p <- p + geom_line(aes(x=intersection.point.exp, color = "expression treshold",  linetype = "expression treshold"), size = 1)
      }
    }


   # set expression tresholds
    if(binaryflag){
      expression.treshold <- intersection.point.bin
    }else{
      inexpression.treshold <- intersection.point.inexp
      expression.treshold <- intersection.point.exp
    }

    # indices expression classes according to tresholds
    if (binaryflag){
      activated.exp <- which(data.keep >= expression.treshold)
      not.exp <- which(data.keep < expression.treshold)
    }else{
     activated.exp <- which(data.keep >= expression.treshold)
     not.exp <- which(data.keep <= inexpression.treshold)
     normal.exp <- which(data.keep > inexpression.treshold & data.keep < expression.treshold)
    }

    # log2(tpm) values to discrete values
    if(binaryflag){
      data.keep[activated.exp] <- 1
      data.keep[not.exp] <- 0
    }else{
      data.keep[activated.exp] <- 1
      data.keep[not.exp] <- -1
      data.keep[normal.exp] <- 0
    }

    # save discretized gene values in output matrix
    discretized.keep[,j] <- data.keep

    # save plot
    if (figflag){

      # add legend
      if (binaryflag){
        legend.title <- "Legend"
        p <- p +
          scale_color_manual(name = legend.title,
                             values = c("treshold"="black", "gauss1"="blue", "gauss2"="red", "signal"="black")) +
          # labels = c("treshold", "gauss1", "gauss2", "signal")) +
          scale_linetype_manual(name = legend.title, values=c("treshold"="dotted",  "gauss1"="solid", "gauss2"="solid", "signal"="solid")) +
          theme(legend.direction = "vertical")
        #guides(linetype=guide_legend(keywidth = 3, keyheight = 1), colour=guide_legend(keywidth = 3, keyheight = 1))

      }else{
        legend.title <- "Legend"
        p <- p +
          scale_color_manual(name = legend.title,
                             # breaks = c("inexpression treshold", "expression treshold", "gauss1", "gauss2", "signal"),
                             values = c( "inexpression treshold" = "blue", "expression treshold" = "red", "gauss1" = "blue", "gauss2" = "red", "signal" = "black")) +
          scale_linetype_manual(name = legend.title, values = c("inexpression treshold" = "dotted", "expression treshold" = "dotted", "gauss1" = "solid", "gauss2"="solid", "signal"="solid" )) +
          theme(legend.direction = "vertical") # legend.key = element_blank()) +
        # guides(linetype=guide_legend(keywidth = 3, keyheight = 1),colour=guide_legend(keywidth = 3, keyheight = 1))

      }

      # name of individual plot file
      plot.file <- paste("Sample_",colnames(signal)[j], ".png", sep = "") # plot file names

      # set file path
      fin.path <- file.path(plot.folder, plot.file)

      # save plot
      ggplot2::ggsave(filename = fin.path)

    }

 }

  print("Save output file...")
  if(is.na(output)){
    write.table(discretized.keep, file=file.path(output.folder, "dge.txt"), sep = "\t", quote = FALSE, row.names=TRUE, col.names=TRUE)
  }else{
    write.table(discretized.keep, file=output,  quote = FALSE, sep = "\t", row.names=TRUE, col.names=TRUE)
  }

}

main <- function() {

  # Command line arguments

  args <- commandArgs(trailingOnly = TRUE)

  option_list <- list(
    make_option(c("-f", "--file"),
                type ="character",
                dest ="filename",
                default = NULL,
                help = "dataset file name",
                metavar = "character"),
   make_option(c("-n", "--normalize"),
                dest="normflag",
                default = FALSE,
                action = "store_true",
                help = "Should the expression data be TPM normalized first? [default %default]",
                metavar="character"),
   make_option(c("-g", "--genome"),
               dest="genomeflag",
               type = "character",
               default = 'h',
               help = "Options 'h' for Hsapiens, 'mm' for Mmusculus genome when no gene length file is provided for TPM normalisation. [default %default]",
               metavar="character"),
   make_option(c("-v", "--vector"),
               type ="character",
               dest ="lenvector",
               default = NULL,
               help = "gene length file name",
               metavar = "character"),
    make_option(c("-p", "--plot"),
                dest = "figflag" ,
                default = FALSE,
                action = "store_true",
                help="Should the program plot the resulting densities of each sample? [default %default]",
                metavar="character"),
    make_option(c("-b", "--binary"),
                dest = "binaryflag" ,
                default = FALSE,
                action = "store_true",
                help="Should the program's output be a binary classification? [default %default]",
                metavar="character"),
    make_option(c("-o", "--out"),
                type="character",
                dest = "output",
                default=NA,
                help="output file name discretization [default= %default]",
                metavar="character"),
   make_option(c("-t", "--tpmout"),
               type="character",
               dest = "output.norm",
               default=NA,
               help="output file name normalization [default= %default]",
               metavar="character"),
   make_option(c("-l", "--lout"),
               type="character",
               dest = "output.len",
               default=NA,
               help="output file name gene lengths [default= %default]",
               metavar="character")
  )


  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);

  # input file name mandatory
  if (is.null(opt$file)){
    print_help(opt_parser)
    stop("The argument 'input file' needs to be supplied", call.=FALSE)
  }

  if (opt$normflag){

    if (is.null(opt$lenvector)){

        #if(is.null(opt$genomeflag)){
        #
        # print_help(opt_parser)
        #  stop("The argument 'genomeflag' -g needs to be supplied for the gene length query", call.=FALSE)
        #
        #}

      gl <- get_genelengths(opt$filename, opt$output.len, opt$genomeflag)
      norm <- normalize(gl[[1]], gl[[2]], opt$output.norm)

      discretize(norm, opt$figflag, opt$binaryflag, opt$output)

      #discretize(normalize(opt$filename, get_genelengths(opt$filename, opt$output.len, opt$genomeflag),
      #                     opt$output.norm), opt$figflag, opt$binaryflag, opt$output)

      #discretize(normalize(opt$filename, get_genelengths(opt$filename, opt$output.len),
      #                     opt$output.norm), opt$figflag, opt$binaryflag, opt$output)

    # print_help(opt_parser)
    # stop("The argument 'gene length file' needs to be supplied for the normalization step", call.=FALSE)
     } else {

    # normalize(opt$filename, opt$lenvector, opt$output.norm)
     discretize(normalize(opt$filename, opt$lenvector, opt$output.norm), opt$figflag, opt$binaryflag, opt$output)
     }

  } else {

     discretize(opt$filename, opt$figflag, opt$binaryflag, opt$output)
  }


}

main()









