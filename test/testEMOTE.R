

library(EMOTE)
library(Rbowtie)
library(Rsamtools)
library(ShortRead)


# index the genome
genome.file <- "test/SA564.fa"
genome.id <- sub("(.fa|.fasta)$","",basename(genome.file))
genome.index <- paste0(genome.file,".bowtie/index")
bowtie_build(genome.file,dirname(genome.index),prefix=basename(genome.index))


# demultiplex EMOTE FASTQ file
demux.report <- EMOTE_demultiplex_fastq("test/EMOTE_example.fq.gz",force=TRUE)

# map demultiplexed reads
bams <- EMOTE_map(genome.index,demux.report$read.counts$demultiplexed.fastq)

# quantify alignments from BAM file
Q <- EMOTE_quantify(bams)

# synchronize counts
q <- local({
  i <- disjoin(unlist(GRangesList(Q)))
  N <- lapply(Q,function(q) {q <- score(q)[findOverlaps(i,q,select="first")];replace(q,is.na(q),0)})
  values(i) <- as(N,"DataFrame")
  i
})

# apply beta-binomial model to compute probability of the positions to be a TSS
library(VGAM)
Q$pval <- pmax(1 - pbetabinom.ab(Q$pRppH,metadata(pRppH)$TotalMap,1+Q$mRppH,1+metadata(mRppH)$TotalMap-Q$mRppH),0)


# compute the probability that the position is a TSS according to a beta-binomial distribution








