

library(EMOTE)
library(Rbowtie)
library(Rsamtools)
library(ShortRead)
library(rtracklayer)

# index the genome
genome.file <- "test/SA564.fa"
genome.id <- sub("(.fa|.fasta)$","",basename(genome.file))
genome.index <- paste0(genome.file,".bowtie/index")
bowtie_build(genome.file,dirname(genome.index),prefix=basename(genome.index))

# demultiplex EMOTE FASTQ file
demux.report <- EMOTE_demultiplex_fastq("test/EMOTE_example.fq.gz",force=TRUE)

# map demultiplexed reads
bam.report <- EMOTE_map(genome.index,demux.report$read.counts$demultiplexed.fastq)

# quantify alignments from BAM file
quantif.report <- EMOTE_quantify(bam.report$bam.file)


Q <- GRangesList(lapply(quantif.report$quantif.file,function(f) {load(f);x}))
names(Q) <- quantif.report$quantif.file

# synchronize counts
q <- local({
  i <- disjoin(unlist(Q))
  N <- lapply(Q,function(q) {q <- score(q)[findOverlaps(i,q,select="first")];replace(q,is.na(q),0)})
  values(i) <- as(N,"DataFrame")
  i
})
names(mcols(q)) <- sub(".quantif.RData$","",basename(names(mcols(q))))




# apply beta-binomial model to compute probability of the positions to be a TSS
library(VGAM)
quantif.report$TotalMap[togroup(Q)]


Q$pval <- pmax(1 - pbetabinom.ab(Q$pRppH,metadata(pRppH)$TotalMap,1+Q$mRppH,1+metadata(mRppH)$TotalMap-Q$mRppH),0)


# compute the probability that the position is a TSS according to a beta-binomial distribution









#
# Process Curran data
#
ss <- read.table("test/EMOTE16.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
ss$RppH <- ifelse(grepl("[+]RppH",ss$SampleId),"+","-")
ss$sampleId <- sub(" [+-]RppH$","",ss$SampleId)
ss <- subset(ss,species=="human")

genome.index <- "/Users/prados/Downloads/Homo_sapiens_NCBI_GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/BowtieIndex/genome"
demux.report <- EMOTE_demultiplex_fastq("/Users/prados/Downloads/EMOTE16_L4_R1.fq.gz",force=TRUE)
bam.report <- EMOTE_map(genome.index,demux.report$read.counts$demultiplexed.fastq)
quantif.report <- EMOTE_quantify(bam.report$bam.file)

Q <- GRangesList(lapply(quantif.report$quantif.file,function(f) {load(f);x}))
M <- SplitDataFrameList(lapply(quantif.report$quantif.file,function(f) {load(f);metadata(x)}))
names(M) <- names(Q) <- quantif.report$quantif.file
M <- stack(M,"quantif.file")

Z <- merge(demux.report$read.counts,ss)
Z <- merge(Z,bam.report,by.x="demultiplexed.fastq",by.y="fq.file")
Z <- merge(Z,quantif.report)
Q <- Q[Z$quantif.file]


q <- local({
  i <- disjoin(unlist(Q))
  N <- lapply(Q,function(q) {q <- score(q)[findOverlaps(i,q,select="first")];replace(q,is.na(q),0)})
  values(i) <- as(N,"DataFrame")
  i
})
names(mcols(q)) <- sub(".quantif.RData$","",basename(names(mcols(q))))
q$wRppH <- q$ATTG + q$TGTT
q$woRppH <- q$AGGA + q$GACG
q <- subset(q,wRppH>=5)

wRppH_totalmap <- sum(subset(Z,index %in% c("ATTG","TGTT"))$TotalMap)
woRppH_totalmap <- sum(subset(Z,index %in% c("AGGA","GACG"))$TotalMap)


library(VGAM)
q$pval <- pmax(1 - pbetabinom.ab(q$wRppH,wRppH_totalmap,1+q$woRppH,1+woRppH_totalmap-q$woRppH),0)
q$fdr <- p.adjust(q$pval,"fdr")
tss <- subset(q,fdr<0.01)
save(tss,file="test/curran_tss_GRCh38.RData")
export.bed(tss,"test/curran_tss_GRCh38.bed")
score(tss) <- pmin(-log10(tss$fdr),30)
export.bw(tss,"test/curran_tss_GRCh38_FDR.bw")





