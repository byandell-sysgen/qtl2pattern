# Read genotype probability object from file
read_probs_calc <- function(chr, datapath, allele = TRUE, probdir = "genoprob") {

  ## Redesign
  ##   allele option may be logical or one of c("allele","pair","snp")
  ##   allele probs may be "apr", "aprobs", "alleleprobs"
  ##   allele pair probs may be "pr", "probs", "genoprobs", "pairprobs"
  ##   SNP probs may be "snp_pr", "snpr", "snpprobs"  
  if(is.logical(allele)) {
    allele <- ifelse(allele, "allele", "pair")
  }
  allele <- match.arg(allele, c("allele", "pair", "snp"))
  
  if(allele == "snp")
    stop("not designed for SNP probabilities yet")
  if((allele == "pair") & length(chr) > 1)
    stop("must supply at most one chr")

  ## Read in probs for selected chromosomes and cbind.
  prefix <- switch(allele,
    allele = "aprobs_",
    pair = "probs_")
  ## Note broader choices in `read_probs_fast()`.
  probs <- convert_probs(
    readRDS(file.path(datapath, probdir,
                      paste0(prefix, chr[1], ".rds"))))

  if(allele == "pair") {

    ## Fix rownames of probs. Begin with "DO-".
    pr <- probs
    tmp <- substring(rownames(pr[[chr]]), 4)
    rownames(pr[[chr]]) <- tmp
    ## Sort them in increasing number order.
    pr[[chr]] <- pr[[chr]][order(as.numeric(tmp)),,, drop = FALSE]
    probs <- modify_object(probs, pr)
  }

  if(length(chr) > 1) {
    for(chri in chr[-1]) {
      probs <- cbind(probs,
                     convert_probs(
                       readRDS(file.path(datapath, probdir,
                                         paste0(prefix, chri, ".rds")))))
    }
  }
  probs
}
#' @export
#' @method dimnames calc_genoprob
dimnames.calc_genoprob <- function (x)
{
  dnames <- lapply(x, dimnames)
  list(ind = dnames[[1]][[1]],
       gen = lapply(dnames, "[[", 2),
       mar = lapply(dnames, "[[", 3))
}
