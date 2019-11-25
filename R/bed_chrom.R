#' bed_chrom
#'
#' For each chromosome create file that let's the user plot the chromosome
#' profile for the chromatin access and tumor mutation burden
#'
#' @param bed BED file with fully processed with 7 columns
#' ("Chromosome","start","end","type","score","strand","signal").
#' @param seg.size Size of segments to be used.
#' @return dat chromosome, position and their
#' estimated chromatin access.
#'
#'@export
#'

bed_chrom <- function(bed,seg.size){

  final <- data.table(do.call('rbind',lapply(c(1:22,"X"),function(x){ #c(1:22,"X")
    print(x)
    ############
    # subset to that chromosome #
    bed.sub <- bed %>%
      filter(Chromosome == paste0("chr",x)) %>%
      select(start,end,signal) %>%
      mutate(Position = trunc(start/seg.size)) %>%
      mutate_all(as.numeric) %>%
      group_by(Position) %>%
      summarise(DNaseI = -sum(signal)) %>% #-mean(signal)
      ungroup() %>%
      select(Position, DNaseI)

    full.dat <- full.dat[complete.cases(full.dat),]
    full.dat$chrom <- rep(x,nrow(full.dat))
    return(full.dat)
  })))
  return(final)
}
