#' chrom_profile
#'
#' For each chromosome create file that let's the user plot the chromosome
#' profile for the chromatin access and tumor mutation burden
#'
#' @param bed BED file with fully processed with 7 columns
#' ("Chromosome","start","end","type","score","strand","signal").
#' @param gen complete maf file with chromosome, chromosome start/end, sample ID
#' @param seg.size Size of segments to be used.
#' @return dat chromosome, position and their
#' estimated chromatin access and TMB (scaled).
#'
#'@export
#'

chrom_profile <- function(bed,gen,seg.size){

  n <- length(unique(gen$submitted_sample_id))
  final <- data.table(do.call('rbind',lapply(c(1:22,"X"),function(x){ #c(1:22,"X")
    print(x)
    ############

    mut.summary <- gen %>%
      filter(chromosome == x) %>%
      mutate(Position = trunc(chromosome_start/(seg.size))) %>%
      group_by(Position) %>%
      summarise(N = n()) %>%
      rename(MutCount = N) %>%
      select(Position,MutCount)

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

    full.dat <- full_join(mut.summary,bed.sub,"Position")
    # full.dat$DNaseI[abs(full.dat$DNaseI) > 2000] <- NA
    full.dat <- full.dat[complete.cases(full.dat),]
    full.dat$chrom <- rep(x,nrow(full.dat))
    # scaled.full <- as.data.frame(scale(full.dat %>%
    #                                      select(MutCount,DNaseI)))
    # scaled.full$Position <- full.dat$Position
    return(full.dat)
  })))
  return(final)
}
