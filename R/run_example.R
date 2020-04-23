#' run_example
#'
#' Recreate chromosome 2 plot from 'Cell-of-origin chromatin organization shapes the mutational landscape of cancer'
#' (https://www.nature.com/articles/nature14221)
#'
#' @param chrom Numeric integer between 1 and 22 specifying which chromosome to plot
#' @param span Numeric value, for the span to be used for smoothing
#' @return Specified chromosome TMB and chromatin access signla overlay
#'
#' @export
#'
#' @import
#' ggplot2
#' cowplot
#' dplyr
#' dtplyr
#' data.table


run_example <- function(chrom = 2, span = 0.025){
  dat <- rbind(bed1,bed2) %>%
    group_by(Position,chrom) %>%
    summarise(TMB = mean(MutCount),
              DNaseI = mean(DNaseI)) %>%
    ungroup() %>%
    filter(DNaseI > -10^6)
  i = 2
  scaled.full <- as.data.frame(scale(dat %>%
                                       filter(chrom == i) %>%
                                       select(TMB,DNaseI)))
  scaled.full$Position <- 0:(nrow(dat %>%
                                    filter(chrom == i))-1)

  span = 0.025
  p <- scaled.full %>%
    mutate(
      DNaseI_Smooth = predict(loess(DNaseI~Position,span = span)),
      MutCount_Smooth = predict(loess(TMB~Position,span = span))
    ) %>%
    ggplot(aes(x=Position)) +
    geom_line(aes(y = MutCount_Smooth, colour = "Scaled mutation count"),size=1.5) +
    geom_line(aes(y = DNaseI_Smooth, colour = "DNaseI"),size=1.5) +
    scale_y_continuous(sec.axis = sec_axis(~.*(-1), name = "Scaled DNaseI")) +
    xlab("Chromosomal position")
  return(p)
}
