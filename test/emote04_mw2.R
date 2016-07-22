

library(EMOTE)
library(Rbowtie)
library(rtracklayer)
library(SummarizedExperiment)

# Initialize the sample sheet for Sepi organism
sample.sheet <- read.table(header=TRUE,sep=",",stringsAsFactors=FALSE,text="
source,barcode,library,assay,RppH,sample
test/EMOTE04_1M.fastq.gz,TTGA,A1-,A1,-,A
test/EMOTE04_1M.fastq.gz,GCTG,A1+,A1,+,A
test/EMOTE04_1M.fastq.gz,ACAC,A2-,A2,-,A
test/EMOTE04_1M.fastq.gz,GGTA,A2+,A2,+,A
")

# Build the Bowtie index for the Sepi Genome
bowtie_build("test/Sepi12228.fa","test/Sepi12228.fa.bowtie")

# Run EMOTE demultiplexing on the pooled FASTQ file
demux.report <- EMOTE_demultiplex_fastq("test/EMOTE04_1M.fastq.gz")

# Merge resulting demultiplex report with the informations from the sample sheet
demux.report <- merge(sample.sheet,demux.report,by=c("source","barcode"))
demux.report <- demux.report[order(demux.report$sample,demux.report$assay,demux.report$RppH),]


# Map demultiplexed reads and generate BAM files
demux.report$bam.file <- sapply(demux.report$demux.fastq,EMOTE_map,bowtie_index="test/Sepi12228.fa.bowtie/index")

# Quantify reads at each genomic position from the BAM files
demux.report$quantif.file <- sapply(demux.report$bam.file,EMOTE_quantify)

# Load the quantification results
q <- EMOTE_read_quantif(demux.report$quantif.file)
colData(q) <- cbind(colData(q),demux.report)


# filter the candidate position by requiring at least 5 read in one of the +RppH condition
q <- q[rowSums(assay(q[,q$RppH=="+"])>=5)>0]


# compute p-value for each assay according to a beta binomial model
mcols(q)$betabino.p <- local({
  cmp <- tapply(colnames(q),list(q$assay,q$RppH),identity)
  pos <- q[,cmp[,"+"]]
  neg <- q[,cmp[,"-"]]
  p <- tss.model(assay(pos),assay(neg),colData(pos)$totalMap[col(pos)],colData(neg)$totalMap[col(neg)])
  colnames(p) <- rownames(cmp)
  p
})

# combine p values of several replicates using fisher methods
mcols(q)$fisher <- local({
  cmp <- tapply(q$assay,q$sample,unique)
  sapply(cmp,function(v) p.combine(mcols(q)$betabino.p[,v]))
})

# compute FDR of combined p-values
mcols(q)$fdr <- apply(mcols(q)$fisher,2,p.adjust,method="fdr")

# Select as TSS the positions where FDR<0.01 for one of the sample
q <- q[rowSums(mcols(q)$fdr<0.01)>0]










#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# TSS clustering
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
k <- reduce(rowRanges(q),min.gapwidth=5,with.revmap=TRUE)
k$rank <- rank(extractList(matrixStats::rowMins(mcols(q)$betabino.p),k$revmap),ties.method="first")
mcols(q)$cluster.id[unlist(k$revmap)] <- togroup(PartitioningByEnd(k$revmap))
mcols(q)$cluster.size <- elementNROWS(k[mcols(q)$cluster.id]$revmap)
mcols(q)$cluster.rank[unlist(k$revmap)] <- unlist(k$rank)




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# TSS annotation with surrounding sequence
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

mcols(q)$seq <- getPromoterSeq(Rsamtools::FaFile("test/Sepi12228.fa"),rowRanges(q))


#     # Test for presence of the TATAAT box in the sequence between pos -5 to -15
#     rowRanges(tss)$TATAAT <- local({
#       pos <- c(-11,-10,-12,-9,-13)
#       ed <- t(neditAt("TATAAT",rowRanges(tss)$seq,at=45 + pos))
#       i <- cbind(1:nrow(ed),max.col(-ed))
#       cbind(pos = pos[i[,2]], err = ed[i])
#     })
#
#     # Test for presence of the TATAAT box in the sequence between pos -25 to -44
#     rowRanges(tss)$TTGACA <- local({
#       pos <- c(-37,-36,-38,-35,-39,-34,-40,-33,-41,-32,-42)
#       ed <- t(neditAt("TTGACA",rowRanges(tss)$seq,at=45 + pos))
#       i <- cbind(1:nrow(ed),max.col(-ed))
#       cbind(pos = pos[i[,2]],err = ed[i])
#     })
#
#     # Generate the pattern that show matching sequences with uppercase letters in a lowercase letter DNA sequence
#     rowRanges(tss)$pattern <- tolower(rowRanges(tss)$seq)
#     subseq(rowRanges(tss)$pattern,45+rowRanges(tss)$TATAAT[,"pos"],width=6) <- toupper(subseq(rowRanges(tss)$pattern,45+rowRanges(tss)$TATAAT[,"pos"],width=6))
#     subseq(rowRanges(tss)$pattern,45+rowRanges(tss)$TTGACA[,"pos"],width=6) <- toupper(subseq(rowRanges(tss)$pattern,45+rowRanges(tss)$TTGACA[,"pos"],width=6))
#     subseq(rowRanges(tss)$pattern,46,width=5) <- toupper(subseq(rowRanges(tss)$pattern,46,width=5))
#     rowRanges(tss)$TATAAT.seq <- subseq(rowRanges(tss)$seq,45+rowRanges(tss)$TATAAT[,"pos"],width=6)
#     rowRanges(tss)$TTGACA.seq <- subseq(rowRanges(tss)$seq,45+rowRanges(tss)$TTGACA[,"pos"],width=6)
#     rowRanges(tss)$TGT_TATAAT.seq <- subseq(rowRanges(tss)$seq,45+rowRanges(tss)$TATAAT[,"pos"]-3,width=3)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# TSS annotation with surrounding genes
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
gff <- import.gff("test/GCF_000007645.1_ASM764v1_genomic.gff.gz",feature.type="gene")
seqlevels(gff) <- sub("\\.[0-9]*$","",seqlevels(gff))

h <- as(findOverlaps(rowRanges(q),gff),"IntegerList")
mcols(q)$overlapping.genes <- extractList(gff$Name,h)
mcols(q)$overlapping.genes <- unstrsplit(mcols(q)$overlapping.genes,";")

nn <- precede(q,gff)
mcols(q)$followingGeneID[!is.na(nn)] <- gff[nn[!is.na(nn)]]$ID
mcols(q)$followingGeneName[!is.na(nn)] <- gff[nn[!is.na(nn)]]$Name
mcols(q)$followingGeneDistance[!is.na(nn)] <- distance(q[!is.na(nn)],gff[nn[!is.na(nn)]])



mcols(q)





