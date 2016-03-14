

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
#' @param invalid.reads path of the FASTQ file that will contains all invalid reads
#' @param yieldSize number of read to read at once
#' @param ... additional parameters are passed to EMOTE_parse_reads
#' @export
#' @import ShortRead
EMOTE_demultiplex_fastq <- function(fq.file,out.dir=paste0(fq.file,".demux"),force=FALSE,invalid.reads=NULL,yieldSize=1e6,...) {
  if (dir.exists(out.dir)) {
    if (!force) stop("out.dir already exists")
    unlink(out.dir,recursive=TRUE)
  }
  dir.create(out.dir)

  if (!is.null(invalid.reads) && file.exists(invalid.reads)) {
    stop("the file invalid.reads already exists")
  }

  f <- FastqStreamer(fq.file,yieldSize)
  on.exit(close(f))

  demux.report <- list(
    summary = data.frame(source=fq.file,total.read.count = 0,total.valid = 0),
    read.counts = data.frame(source=character(),index=character(),demultiplexed.fastq=character(),n=integer())
  )
  while (length(fq <- yield(f))) {
    X <- EMOTE_parse_reads(fq,...)
    X$read@id <- xscat(BStringSet(X$umi),":",id(X$read))
    if (!is.null(invalid.reads)) writeFastq(fq[!X$is_valid_EMOTE_read],invalid.reads,mode="a")
    X <- subset(X,is_valid_EMOTE_read)
    X <- split(X$read,as.character(X$barcode))
    fn <- file.path(out.dir,paste0(names(X),".fq.gz"))
    mapply(writeFastq, X, fn, MoreArgs = list(compress=TRUE,mode="a"))

    demux.report <- within(demux.report,{
      summary$total.read.count <- summary$total.read.count + length(fq)
      summary$total.valid <- summary$total.valid + sum(elementLengths(X))
      read.counts <- merge(all=TRUE,read.counts,data.frame(source=fq.file,index=names(X),demultiplexed.fastq=fn,n=elementLengths(X)),by=c("source","index","demultiplexed.fastq"))
      read.counts$n <- rowSums(cbind(read.counts$n.x,read.counts$n.y),na.rm=TRUE)
      read.counts$n.x <- read.counts$n.y <- NULL
    })
  }
  demux.report
}


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
EMOTE_map <- function(bowtie_index,fq.files,bam.files=sub("(.fastq|.fq)(.gz)?$",".bam",fq.files),threads=3) {
  fq.files <- as.character(fq.files)
  if (!all(grepl("\\.bam$",bam.files))) stop("bam.files must end with .bam suffix")

  mapply(function(fq.file,bam.file) {
    if (grepl(".gz$",fq.file)) {
      gz.file <- fq.file
      gunzip(gz.file,fq.file <- tempfile(fileext=".fq"))
    }
    sam.file <- tempfile(fileext=".sam")
    bowtie(sam=TRUE,best=TRUE,M=1,sequences=fq.file,index=bowtie_index,outfile=sam.file,threads=threads,v=1)
    asBam(sam.file,sub(".bam$","",bam.file),indexDestination=TRUE,overwrite=TRUE)
    bam.file
  },fq.files,bam.files)

  data.frame(reference=bowtie_index,demultiplexed.fastq=fq.files,bam.file=bam.files,stringsAsFactors = FALSE)
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
EMOTE_quantify <- function(
    bam.files,
    quantif.files=sub(".bam$",".quantif.RData",bam.files),
    remove.umi.duplicate=TRUE,
    mode=c("unambiguous","ambiguous","all"),
    yieldSize=NA_integer_
  ) {

  # helper function that do the job for a single BAM file
  quantify_single_bam <- function(bam.file,remove.umi.duplicate,mode,yieldSize) {
    bf <- BamFile(bam.file,yieldSize=yieldSize)
    open(bf);on.exit(close(bf))

    totalAmbiguous <- totalUnambiguous <- totalMap <- covPos <- covNeg <- 0
    while(length(gr <- readGAlignments(bf,use.names = TRUE,param = ScanBamParam(what="mapq")))>0) {
      # update summary statistics
      totalMap <- totalMap + length(gr)
      totalAmbiguous <- totalAmbiguous + sum(mcols(gr)$mapq==0)
      totalUnambiguous <- totalUnambiguous + sum(mcols(gr)$mapq>0)

      # subset the reads depending on the mode used
      names(gr) <- sub(":.*","",names(gr)) # parse UMI from read name prefix
      switch(mode,
             unambiguous = {gr <- subset(gr,mapq>0)},
             ambiguous = {gr <- subset(gr,mapq==0)},
             all = {}
      )
      gr <- granges(gr)

      # remove PCR duplicates (with identical position for a given UMI)
      if (remove.umi.duplicate) {
        gr <- split(resize(gr,1,use.names=FALSE),names(gr))
        gr <- unlist(unique(gr))
      }

      # update coverage values
      covPos <- coverage(subset(gr,strand!="+")) + covPos
      covNeg <- coverage(subset(gr,strand=="-")) + covNeg
    }
    covPos <- subset(GRanges(covPos,strand="+"),score>0)
    covNeg <- subset(GRanges(covNeg,strand="-"),score>0)

    cov <- c(covPos,covNeg)
    metadata(cov) <- data.frame(stringsAsFactors=FALSE,
      TotalMap=totalMap,
      TotalAmbiguous=totalAmbiguous,
      TotalUnambiguous=totalUnambiguous
    )
    cov
  }

  # check parameters
  mode <- match.arg(mode)
  if (length(quantif.files) != length(bam.files)) stop("plus.bw.files must be of the same length than bam.files")

  out <- mapply(function(bam.file,quantif.file) {
    x <- quantify_single_bam(bam.file,remove.umi.duplicate,mode,yieldSize)
    save(x,file=quantif.file)
    metadata(x)$quantif.file <- quantif.file
    metadata(x)
  },bam.files,quantif.files,SIMPLIFY = FALSE)
  stack(SplitDataFrameList(out),"bam.file")
}




EMOTE_bw_merge <- function() {

}



