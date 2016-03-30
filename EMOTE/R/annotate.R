



#' Retreive the promoter sequence around given genomic ranges
#' @param fa a FaFile to search the sequence in
#' @param gr the genomic range used as reference to compute the promoter range
#' @param upstream number of required bases upstream
#' @param downstream number of required bases downstream
#' @return a DNAStringSet with the promoters sequences
#' @export
getPromoterSeq <- function(fa,gr,upstream=45,downstream=5) {
  seq.gr <- promoters(granges(gr),upstream,downstream)
  valid.gr <- trim(seq.gr)
  seq <- getSeq(fa,valid.gr)
  seqTemplate <- paste0(rep("N",upstream+downstream),collapse="")
  seqTemplate <- DNAStringSet(rep(seqTemplate,length(seq)))
  subseq(seqTemplate,
         1 + abs(start(resize(valid.gr,1)) - start(resize(seq.gr,1))),
         nchar(seqTemplate) - abs(end(resize(valid.gr,1,fix="end")) - end(resize(seq.gr,1,fix="end")))
  ) <- seq
  seqTemplate
}
