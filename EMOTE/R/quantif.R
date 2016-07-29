

#' Parse reads generated with EMOTE protocol
#'
#' The function takes the given ShortRead object, split the reads into elements
#' expected to by found in EMOTE reads, and test for their validity.
#'
#'
#' @param sr a ShortRead object to parse
#' @param max.recognition.mismatch number of mismatch allowed in the recognition
#'        sequence
#' @param valid.barcodes a character vector listing valid barcodes. If NULL, all
#'        barcodes are considered valid.
#' @param first.nuc.len length of the first nucleotide sequence before EMOTE
#'        recognition sequence start
#' @param recognition.sequence a character with the recognition sequence to
#'        detect in the reads
#' @param umi.length length of the UMI sequence
#' @param control.sequence a character specifying the control sequence to look
#'        for in the read
#' @return A DataFrame where the columns contains the elements parsed in the
#'         input reads or a logical checking for the validity of the elements.
#'         The fields are:
#'  \itemize{
#'        \item first_nt: the first nucleotide of the read;
#'        \item barcode: the 4 nucleotides at position 2-5 in the read, which is
#'              expected to be the the barcode sequence;
#'        \item recognition: the 15 nucleotides at position 6-20 in the read,
#'              which is expected to be the EMOTE recognition sequence equal to
#'              CGGCACCAACCGAGG;
#'        \item umi: the 7 nucleotides at position 21-27 in the read, which is
#'              expected to be an universal molecular identifier (UMI);
#'        \item control: the 3 nucleotides at position 28-30 in the read, which
#'              is expected to be equal to the EMOTE control sequence ACG;
#'        \item read: the end of the read from position 31 which is the read
#'              sequence expected to map on the genome
#'        \item is_valid_recognition: a logical which is TRUE when the
#'              reconigiton sequence match with the sequence CGGCACCAACCGAGG
#'              with max.mismatch error tolerance;
#'        \item is_valid_control: a logical which is TRUE when the control
#'              sequence perfectly match with the sequence CGC;
#'        \item is_valid_umi: a logical which is TRUE when the UMI sequence only
#'              contains A,C or G nucleotides;
#'        \item is_valid_barcode: a logical which is TRUE when the barcode
#'              perfectly match with one of the given valid.barcodes, or when
#'              valid.barcodes is NULL;
#'        \item is_valid_EMOTE_read: a logical which is TRUE when
#'              is_valid_recognition, is_valid_control, is_valid_umi, and
#'              is_valid_barcode are all TRUE.
#' }
#' @export
#' @import S4Vectors
#' @import Biostrings
EMOTE_parse_reads <- function(
  sr,
  max.recognition.mismatch=1L,
  valid.barcodes=DNAStringSet(c("TACA","GTAT","CGTC","AAGT","ACAC","GGTA","GCCT","TCGG","CAAG","TTGA","GCTG","CCGA","CTCG","AGGA","ATTG","GACG","TGTT")),
  first.nuc.len=1L,
  recognition.sequence="CGGCACCAACCGAGG",
  umi.length=7L,
  control.sequence="CGC"
) {
  barcode.length <- unique(nchar(valid.barcodes))
  if (length(barcode.length)>1) stop("All barcode sequences must be of the same length")

  bar.start <- 1L + first.nuc.len
  rec.start <- bar.start + barcode.length
  umi.start <- rec.start + nchar(recognition.sequence)
  control.start <- umi.start + umi.length
  read.start <- control.start + nchar(control.sequence)

  X <- DataFrame(first_nt = narrow(sread(sr),1L,first.nuc.len),
       barcode = narrow(sread(sr),bar.start,rec.start-1L),
       recognition = narrow(sread(sr),rec.start,umi.start-1L),
       umi = narrow(sread(sr),umi.start,control.start-1L),
       control = narrow(sread(sr),control.start,read.start-1L),
       read = narrow(sr,read.start)
  )
  X$is_valid_recognition <- vcountPattern(recognition.sequence,X$recognition,max.mismatch=max.recognition.mismatch)>0
  X$is_valid_control <- X$control==DNAString(control.sequence)
  X$is_valid_umi <- as.vector(letterFrequency(X$umi,"ACG")==width(X$umi))
  X$is_valid_barcode <- (is.null(valid.barcodes) | (X$barcode %in% valid.barcodes))
  X$is_valid_EMOTE_read <- X$is_valid_recognition & X$is_valid_control & X$is_valid_umi & X$is_valid_barcode
  return(X)
}



#' Demultiplex a EMOTE FASTQ file according to barcode sequence of the reads.
#'
#' Parse the reads from the input FASTQ file and store valid mapping sequences
#' in separate FASTQ files in the given output directory. The name of the file
#' in the ouput directory depends on the barcode sequence. At the end of the
#' process demultiplexing statistics are reported in the file
#' demultiplex_report.txt.
#'
#' When the file demultiplex_report.txt already exists in the output directory
#' and force=FALSE, no demultiplexing is performed, but the report is load and
#' the function return with a warning message. When force=TRUE, the output
#' directory is recursively deleted and the demultiplex process if forced to
#' run.
#'
#' By default the function processes the input fastq by chunks of 1e6 reads.
#' The number of read processed by the function can be tuned with the parameter
#' yieldSize to control memory usage.
#'
#' @param fq.file a character with the path of the FASTQ file (eventually
#'        gzipped) containing single-end EMOTE reads
#' @param out.dir a character with the path of the output directory where the
#'        demultiplexed FASTQ files will be stored
#' @param force when TRUE recursively remove the content of out.dir before
#'        running instead of stopping the function with an error message
#' @param invalid.reads logical saying whether invalid reads should be outputed
#'        in the file INVALID.fastq.gz
#' @param yieldSize number of read to process per iteration
#' @param ... additional parameters are passed to EMOTE_parse_reads
#' @return a data.frame with demultiplexing statistics
#' @export
#' @import ShortRead
EMOTE_demultiplex_fastq <- function(fq.file,out.dir=paste0(fq.file,".demux"),force=FALSE,yieldSize=1e6,invalid.reads=TRUE,...) {
  demux.report.file <- file.path(out.dir,"demultiplex_report.txt")
  if (file.exists(demux.report.file) && !force) {
      warning("The output directory already exists => load demultiplexing report from cache. Use force=TRUE to force building the output.\n")
      X <- read.table(demux.report.file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
      if (!invalid.reads) X <- subset(X,barcode!="INVALID")
      return(X)
  }

  # create the ouptut directory
  dir.create(out.dir,showWarnings=FALSE)

  # initilialize a streamer on the input file
  f <- FastqStreamer(fq.file,yieldSize)
  on.exit(close(f))

  # initialize reporting data
  R <- data.frame(source=fq.file,source.read.count=0,source.valid.read=0,source.invalid.contol=0,source.invalid.barcode=0,source.invalid.umi=0,source.invalid.recognition=0,stringsAsFactors=FALSE)
  R.demux <- data.frame(source=character(),barcode=character(),demux.fastq=character(),demux.read.count=integer(),stringsAsFactors=FALSE)

  # loop over all reads of the input FASTQ file
  while (length(fq <- yield(f))) {
    # parse the reads
    X <- EMOTE_parse_reads(fq,...)
    X$barcode_group <- as.character(X$barcode)
    X$barcode_group[!X$is_valid_EMOTE_read] <- "INVALID"
    X$barcode_group <- factor(X$barcode_group)

    # output read sequences into separate files
    FQ <- split(fq,X$barcode_group)
    if (!invalid.reads) FQ$INVALID <- NULL
    fn <- file.path(out.dir,paste0(names(FQ),".fastq.gz"))
    mapply(writeFastq, FQ, fn, MoreArgs = list(compress=TRUE,mode="a"))

    # update report
    R$source.read.count <- R$source.read.count + length(fq)
    R$source.valid.read <- R$source.valid.read + sum(X$is_valid_EMOTE_read)
    R$source.invalid.contol <- R$source.invalid.contol + sum(!X$is_valid_control)
    R$source.invalid.barcode <- R$source.invalid.barcode + sum(!X$is_valid_barcode)
    R$source.invalid.umi <- R$source.invalid.umi + sum(!X$is_valid_umi)
    R$source.invalid.recognition <- R$source.invalid.recognition + sum(!X$is_valid_recognition)

    # update demultiplex_report
    R.demux <- merge(all=TRUE,R.demux,data.frame(source=fq.file,barcode=names(FQ),demux.fastq=fn,demux.read.count=elementNROWS(FQ)),by=c("source","barcode","demux.fastq"))
    R.demux$demux.read.count <- rowSums(cbind(R.demux$demux.read.count.x,R.demux$demux.read.count.y),na.rm=TRUE)
    R.demux$demux.read.count.x <- R.demux$demux.read.count.y <- NULL
  }
  X <- merge(R,R.demux,by="source",all=TRUE)
  write.table(X,file=demux.report.file,sep="\t",row.names=FALSE)
  return(X)
}


#' Map reads to a genome with Rbowtie
#'
#' Parse the EMOTE reads in the given FASTQ file to extract the mapping sequence
#' and the UMI. Create an intermediate FASTQ file with this mapping sequence and
#' the UMI sequence added as a prefix of the reads headers. Map the intermediate
#' FASTQ file with Rbowtie and generate a BAM.
#'
#' When the BAM file to generate already exists, the function simply return its name
#'
#' @param bowtie_index a length-one character whith the name of the bowtie index
#'        onto which the reads have to be mapped
#' @param fq.file name of the EMOTE FASTQ file to map (usualy a demultiplexed FASTQ)
#' @param bam.file name of the generated BAM file
#' @param intermediate.fastq.file name of the intermediate FASTQ file that is generated
#' @param force when TRUE the BAM file is generated wven if it already exists
#' @param threads number of thread to be used by Rbowtie during the mapping process
#' @param yieldSize the number of read processed per iteration when parsing the input FASTQ file to generate the intermediate file
#' @param ... additional paramaters passed to EMOTE_parse_reads()
#' @return the name of the generated BAM file
#' @seealso bowtie_build
#' @export
#' @import Rbowtie
#' @import Rsamtools
EMOTE_map <- function(
  bowtie_index,
  fq.file,bam.file=sub("(.fastq|.fq)(.gz)?$",".bam",fq.file),
  intermediate.fastq.file=tempfile(fileext=".fastq"),
  force=FALSE,threads=3,yieldSize=1e6,...
) {
  if (!all(grepl("\\.bam$",bam.file))) stop("bam.file must end with .bam suffix")
  if (!force && file.exists(bam.file)) {
    warning("The BAM file already exists. Use force=TRUE to force recomputing")
    return(bam.file)
  }
  fq.file <- as.character(fq.file)
  bowtie_index <- as.character(bowtie_index)

  # Process the fastq file to add the UMI as header prefix of each read
  local({
    f <- FastqStreamer(fq.file,yieldSize)
    on.exit(close(f))
    while (length(fq <- yield(f))) {
      X <- EMOTE_parse_reads(fq,...)
      X$read@id <- xscat(BStringSet(X$umi),":",id(X$read))
      writeFastq(X$read,intermediate.fastq.file,compress=FALSE,mode="a")
    }
  })

  # map with Bowtie
  #     -v 1: no more than V mismatches in the alignment
  #   --best: guarantee that reported alignments are best in number of mismatches
  #     -M 1: report only one random alignment when more than M alignment for the read
  sam.file <- tempfile(fileext=".sam")
  bowtie(sam=TRUE,best=TRUE,M=1,sequences=intermediate.fastq.file,index=bowtie_index,outfile=sam.file,threads=threads,v=1)
  asBam(sam.file,sub(".bam$","",bam.file),indexDestination=TRUE,overwrite=TRUE)

  return(bam.file)
}


#' Helper function to write an EMOTE quantification file
#' @import rtracklayer
EMOTE_write_quantif <- function(covPos,covNeg,R,file) {
  dir.create(file,showWarnings=FALSE)
  report.file <- file.path(file,"report.txt")
  pos.bw.file <- file.path(file,"positive_strand.bw")
  neg.bw.file <- file.path(file,"negative_strand.bw")
  write.table(R,file=report.file,sep="\t",row.names=FALSE)
  export.bw(covPos,pos.bw.file)
  export.bw(covNeg,neg.bw.file)
}


#' Load one or multiple EMOTE quantification files
#'
#' @param quantif.files a character vector listing the quantification directories to load
#' @return a SummarizedExperiment object containing the quantification values for all samples
#'         at each non-empty genomic location
#' @export
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import rtracklayer
EMOTE_read_quantif <- function(quantif.files) {
  report.files <- file.path(quantif.files,"report.txt")
  pos.bw.files <- file.path(quantif.files,"positive_strand.bw")
  neg.bw.files <- file.path(quantif.files,"negative_strand.bw")

  # load thre report statistics
  R <- lapply(report.files,read.table,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  R <- unlist(SplitDataFrameList(setNames(R,quantif.files)))

  # load the counts from the BigWig files
  cov <- local({
    read.bw <- function(f,strand="*") {
      x <- RleList(import.bw(f,as="RleList"))
      pos <- which(x>0)
      seq <- Rle(factor(names(pos),names(pos)),elementNROWS(pos))
      n <- unlist(x[x>0],use.names = FALSE)
      gr <- GRanges(seq,IRanges(unlist(pos),width=1),strand=strand,score=n)
      seqlengths(gr) <- elementNROWS(x)
      gr
    }
    covPos <- GRangesList(lapply(pos.bw.files,read.bw,strand="+"))
    covNeg <- GRangesList(lapply(neg.bw.files,read.bw,strand="-"))
    setNames(pc(covPos,covNeg),quantif.files)
  })

  # build the SummarizedExperiment object
  gr <- unique(unlist(cov,use.names = FALSE))
  mcols(gr) <- NULL
  gr <- sort(gr)
  col <- togroup(PartitioningByEnd(cov))
  row <- unlist(match(cov,gr))
  n <- matrix(0,length(gr),length(cov))
  n[cbind(row,col)] <- score(unlist(cov))

  SummarizedExperiment(n,rowRanges=gr,colData=R)
}


#' Quantify number of read starting at each position of the genome
#'
#' Quantify number of read starting at each position of the genome with optional
#' UMI correction, and store the result in an output directory
#'
#' @param bam.file the input bam files to process
#' @param quantif.file name of the output directory that will contain the quantification results
#' @param remove.umi.duplicate a logical telling whether the umi correction must be done
#' @param mode the quantification mode to use: either all reads, either only ambiguous reads or only unambiguous reads.
#'   Reads with a mapping quality of 0 are considered ambiguously mapped.
#' @param yieldSize yieldSize value passed to BamFile to process input file by chunks when limited amount of memory is available
#' @return name of the generated quantification directory
#' @export
#' @import GenomicAlignments
#' @import Rsamtools
#' @import SummarizedExperiment
EMOTE_quantify <- function(
    bam.file,
    quantif.file=paste0(bam.file,".quantif"),
    remove.umi.duplicate=TRUE,
    mode=c("unambiguous","ambiguous","all"),
    yieldSize=NA_integer_,
    force=FALSE
) {
    # check parameters
    mode <- match.arg(mode)

    # If the BED file already exists, load it and return
    if (!force && file.exists(quantif.file)) {
      warning("quantification file already exists. Use force=TRUE to force recomputing")
      return(quantif.file)
    }
    # open the BAM file and set the number of read to retreive at each iteration
    bf <- BamFile(bam.file,yieldSize=yieldSize)
    open(bf);on.exit(close(bf))

    # initialize the report and the coverage variables that will be updated at each iteration
    covPos <- covNeg <- 0
    R <- DataFrame(totalMap=0,totalAmbiguous=0,totalUnambiguous=0,totalDedup=0)

    # iterate over all the reads of the BAM file by chunk
    while(length(gr <- readGAlignments(bf,use.names=TRUE,param=ScanBamParam(what="mapq")))>0) {

      # update quantif_report statistics with the new chunk
      R$totalMap <- R$totalMap + length(gr)
      R$totalAmbiguous <- R$totalAmbiguous + sum(mcols(gr)$mapq==0)
      R$totalUnambiguous <- R$totalUnambiguous + sum(mcols(gr)$mapq>0)

      # subset the reads depending on the mode used
      names(gr) <- sub(":.*","",names(gr)) # parse UMI from read name prefix
      gr <- switch(mode,
                   unambiguous = gr[mcols(gr)$mapq>0],
                   ambiguous = gr[mcols(gr)$mapq==0],
                   all = gr
      )
      gr <- granges(gr)

      # remove PCR duplicates (with identical position for a given UMI)
      if (remove.umi.duplicate) {
        gr <- split(resize(gr,1,use.names=FALSE),names(gr))
        gr <- unlist(unique(gr))
      }
      R$totalDedup <- R$totalDedup + length(gr)

      # update coverage values
      covPos <- coverage(gr[strand(gr)!="-"]) + covPos
      covNeg <- coverage(gr[strand(gr)=="-"]) + covNeg
    }
    EMOTE_write_quantif(covPos,covNeg,R,file=quantif.file)
    return(quantif.file)
}



