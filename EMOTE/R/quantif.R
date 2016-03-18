

#' Parse reads generated with EMOTE protocol
#' @param sr a ShortRead object to parse
#' @param valid.barcodes a character vector listing valid barcodes. If NULL, all
#'        barcodes are considered valid.
#' @export
#' @import S4Vectors
#' @import Biostrings
EMOTE_parse_reads <- function(sr,max.mismatch=1,valid.barcodes=DNAStringSet(c("TACA","GTAT","CGTC","AAGT","ACAC","GGTA","GCCT","TCGG","CAAG","TTGA","GCTG","CCGA","CTCG","AGGA","ATTG","GACG","TGTT"))) {
  X <- DataFrame(first_nt = narrow(sread(sr),1,1),
       barcode = narrow(sread(sr),2,5),
       recognition = narrow(sread(sr),6,20),
       umi = narrow(sread(sr),21,27),
       control = narrow(sread(sr),28,30),
       read = narrow(sr,31)
  )
  X <- within(X,{
    is_valid_recognition <- vcountPattern("CGGCACCAACCGAGG",recognition,max.mismatch=.(max.mismatch))>0
    is_valid_control <- control==DNAString("CGC")
    is_valid_umi <- as.vector(letterFrequency(umi,"ACG")==width(umi))
    is_valid_barcode <- (is.null(.(valid.barcodes)) | (barcode %in% .(valid.barcodes)))
    is_valid_EMOTE_read <- is_valid_recognition & is_valid_control & is_valid_umi & is_valid_barcode
  })
  return(X)
}



#' Parse and demultiplex a EMOTE FASTQ file according to barcode sequence of the reads.
#'
#' Parse the input reads and store valid mapping sequences in separate FASTQ files according to the barcode sequence.
#' Additionnaly, the UMI sequence is prepended to the read names to keep a trace during mapping and use it during quantification
#'
#' @param fq.file a character with the path of the FASTQ file (eventually gzipped) containing single-end EMOTE reads
#' @param out.dir a character with the path of the output directory where the demultiplexed FASTQ files will be stored
#' @param force when TRUE recursively remove the content of out.dir instead of stopping with an error message
#' @param invalid.reads logical saying whether invalid reads should be outputed
#' @param yieldSize number of read to read at once
#' @param ... additional parameters are passed to EMOTE_parse_reads
#' @export
#' @import ShortRead
EMOTE_demultiplex_fastq <- function(fq.file,out.dir=paste0(fq.file,".demux"),force=FALSE,yieldSize=1e6,invalid.reads=TRUE,...) {
  if (dir.exists(out.dir)) {
    if (!force) {
      cat("The output directory already exists => load demultiplexing report from cache. Use force=TRUE to force building the output.\n")
      parsing_report <- read.table(file.path(out.dir,"parsing_report.txt"),sep="\t",header=TRUE)
      demultiplex_report <- read.table(file.path(out.dir,"demultiplex_report.txt"),sep="\t",header=TRUE)
      return(merge(parsing_report,demultiplex_report,by="source",all=TRUE))
    }
    unlink(out.dir,recursive=TRUE)
  }
  dir.create(out.dir)

  f <- FastqStreamer(fq.file,yieldSize)
  on.exit(close(f))


  parsing_report <- data.frame(source = fq.file,read.count = 0,valid.read = 0,invalid.contol = 0,invalid.barcode = 0,invalid.umi = 0,invalid.recognition = 0)
  demultiplex_report <- data.frame(source=character(),barcode=character(),mapseq.fastq=character(),num_mapseq=integer())

  while (length(fq <- yield(f))) {
    # parse the reads
    X <- EMOTE_parse_reads(fq,...)
    X$barcode_group <- as.character(X$barcode)
    X$barcode_group[!X$is_valid_EMOTE_read] <- "invalid_EMOTE_read"
    X$barcode_group <- factor(X$barcode_group)

    # output read sequences into separate files
    local({
      FQ <- split(fq,X$barcode_group)
      if (!invalid.reads) FQ$invalid_EMOTE_read <- NULL
      fn <- file.path(out.dir,paste0(names(FQ),".fastq.gz"))
      mapply(writeFastq, FQ, fn, MoreArgs = list(compress=TRUE,mode="a"))
    })

    # output valid mapping sequences
    Z <- subset(X,is_valid_EMOTE_read)
    Z$read@id <- xscat(BStringSet(Z$umi),":",id(Z$read))
    Z <- split(Z$read,Z$barcode_group[drop=TRUE])
    mapseq.fastq <- file.path(out.dir,paste0(names(Z),"_mapseq.fastq.gz"))
    mapply(writeFastq, Z, mapseq.fastq, MoreArgs = list(compress=TRUE,mode="a"))

    # update parsing_report
    parsing_report$read.count <- parsing_report$read.count + length(fq)
    parsing_report$valid.read <- parsing_report$valid.read + sum(X$is_valid_EMOTE_read)
    parsing_report$invalid.contol <- parsing_report$invalid.contol + sum(!X$is_valid_control)
    parsing_report$invalid.barcode <- parsing_report$invalid.barcode + sum(!X$is_valid_barcode)
    parsing_report$invalid.umi <- parsing_report$invalid.umi + sum(!X$is_valid_umi)
    parsing_report$invalid.recognition <- parsing_report$invalid.recognition + sum(!X$is_valid_recognition)

    # update demultiplex_report
    demultiplex_report <- merge(all=TRUE,demultiplex_report,data.frame(source=fq.file,barcode=names(Z),mapseq.fastq=mapseq.fastq,num_mapseq=elementLengths(Z)),by=c("source","barcode","mapseq.fastq"))
    demultiplex_report$num_mapseq <- rowSums(cbind(demultiplex_report$num_mapseq.x,demultiplex_report$num_mapseq.y),na.rm=TRUE)
    demultiplex_report$num_mapseq.x <- demultiplex_report$num_mapseq.y <- NULL
  }
  write.table(parsing_report,file=file.path(out.dir,"parsing_report.txt"),sep="\t",row.names=FALSE)
  write.table(demultiplex_report,file=file.path(out.dir,"demultiplex_report.txt"),sep="\t",row.names=FALSE)
  return(merge(parsing_report,demultiplex_report,by="source",all=TRUE))
}

#' Helper function that unzip a gzfile to the given destination
gunzip <- function(gz.file,dest.file) {
  gz.file <- gzfile(gz.file,"rb");on.exit(close(gz.file))
  dest.file <- file(dest.file,"wb");on.exit(close(dest.file))
  while (length(data <- readBin(gz.file,raw(),n=20e6)))
    writeBin(data,dest.file)
}

#' Map reads to a genome with Rbowtie
#' @export
#' @import Rbowtie
#' @import Rsamtools
EMOTE_map <- function(bowtie_index,fq.file,bam.file=sub("(.fastq|.fq)(.gz)?$",".bam",fq.file),force=FALSE,threads=3) {
  if (!all(grepl("\\.bam$",bam.file))) stop("bam.file must end with .bam suffix")
  if (!force && file.exists(bam.file)) return(bam.file)
  fq.file <- as.character(fq.file)
  bowtie_index <- as.character(bowtie_index)

  if (all(grepl(".gz$",fq.file))) {
    gz.file <- fq.file
    gunzip(gz.file,fq.file <- tempfile(fileext=".fq"))
  }
  sam.file <- tempfile(fileext=".sam")

  # map with Bowtie
  #     -v 1: no more than V mismatches in the alignment
  #   --best: guarantee that reported alignments are best in number of mismatches
  #     -M 1: report only one random alignment when more than M alignment for the read
  bowtie(sam=TRUE,best=TRUE,M=1,sequences=fq.file,index=bowtie_index,outfile=sam.file,threads=threads,v=1)
  asBam(sam.file,sub(".bam$","",bam.file),indexDestination=TRUE,overwrite=TRUE)

  return(bam.file)
}


#' Quantify number of read starting at each position of the genome
#'
#' Quantify number of read starting at each position of the genome with optional
#' UMI correction
#'
#' @param bam.files the input bam files to process
#' @param quantif.files names of the output RData file that store the coverage
#' @param remove.umi.duplicate a logical telling whether the umi correction must be done
#' @param mode the quantification mode to use: either all reads, either only ambiguous reads or only unambiguous reads.
#'   Reads with a mapping quality of 0 are considered ambiguously mapped.
#' @param yieldSize yieldSize value passed to BamFile to process input file by chunks when limited amount of memory is available
#' @export
#' @import rtracklayer
#' @import GenomicAlignments
EMOTE_quantify <- function(
    bam.file,
    remove.umi.duplicate=TRUE,
    mode=c("unambiguous","ambiguous","all"),
    yieldSize=NA_integer_,
    pos.bw.file=paste0(bam.file,".pos.bw"),
    neg.bw.file=paste0(bam.file,".neg.bw"),
    quantif.report.file=paste0(bam.file,".quantif_report.txt"),
    force=FALSE
) {
    # check parameters
    mode <- match.arg(mode)

    if (!force & file.exists(quantif.report.file)) {
      quantif_report <- read.table(quantif.report.file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
      return(quantif_report)
    }

    # open the BAM file and set the number of read to retreive at each iteration
    bf <- BamFile(bam.file,yieldSize=yieldSize)
    open(bf);on.exit(close(bf))

    # initialize the report and the coverage variables that will be updated at each iteration
    covPos <- covNeg <- 0
    quantif_report <- data.frame(bam.file=bam.file,pos.bw.file=pos.bw.file,neg.bw.file=neg.bw.file,totalMap=0,totalAmbiguous=0,totalUnambiguous=0,stringsAsFactors=FALSE)

    # iterate over all the reads of the BAM file by chunk
    while(length(gr <- readGAlignments(bf,use.names = TRUE,param = ScanBamParam(what="mapq")))>0) {

        # update quantif_report statistics with the new chunk
        quantif_report$totalMap <- quantif_report$totalMap + length(gr)
        quantif_report$totalAmbiguous <- quantif_report$totalAmbiguous + sum(mcols(gr)$mapq==0)
        quantif_report$totalUnambiguous <- quantif_report$totalUnambiguous + sum(mcols(gr)$mapq>0)

        # subset the reads depending on the mode used
        names(gr) <- sub(":.*","",names(gr)) # parse UMI from read name prefix
        gr <- switch(mode,
               unambiguous = subset(gr,mapq>0),
               ambiguous = subset(gr,mapq==0),
               all = gr
        )
        gr <- granges(gr)

        # remove PCR duplicates (with identical position for a given UMI)
        if (remove.umi.duplicate) {
          gr <- split(resize(gr,1,use.names=FALSE),names(gr))
          gr <- unlist(unique(gr))
        }

        # update coverage values
        covPos <- coverage(subset(gr,strand!="-")) + covPos
        covNeg <- coverage(subset(gr,strand=="-")) + covNeg
    }

    export.bw(covPos,pos.bw.file)
    export.bw(covNeg,neg.bw.file)
    write.table(quantif_report,file=quantif.report.file,sep="\t",row.names=FALSE)

    return(quantif_report)
}



#' Merge several bigwig files into a SummarizedExperiment
#'
#' Merge several bigwig files into a SummarizedExperiment
#' @param pos.bw the input bam files to process
#' @param neg.bw names of the output RData file that store the coverage
#' @export
EMOTE_bw_merge <- function(pos.bw.files,neg.bw.files) {
  Q <- mapply(function(p,n) {
    pos <- GRanges(import.bw(p),strand="+")
    neg <- GRanges(import.bw(p),strand="-")
    subset(c(pos,neg),score>0)
  },pos.bw.files,neg.bw.files)
  Q <- stack(GRangesList(unname(Q)))

  i <- disjoin(Q)
  h <- findOverlaps(i,Q)

  N <- matrix(0,length(i),nlevels(Q$sample))
  N[cbind(queryHits(h),as.integer(Q$sample)[subjectHits(h)])] <- score(Q)[subjectHits(h)]

  E <- SummarizedExperiment(N,rowRanges=i,colData=DataFrame(pos.bw.file=pos.bw.files,neg.bw.file=neg.bw.files))
  return(E)
}



