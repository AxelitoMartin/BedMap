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

  final <- data.table(do.call('rbind',lapply(c(1:22,"X"),function(x){
    print(x)
    ############

    mut.summary <- gen %>%
      filter(chromosome == x) %>%
      mutate(Position = trunc(chromosome_start/(seg.size))) %>%
      group_by(Position) %>%
      summarise(N = n()/length(unique(submitted_sample_id))) %>%
      rename(MutCount = N)





    # subset to that chromosome #
    bed.sub <- bed %>%
      filter(Chromosome == paste0("chr",x)) %>%
      select(start,end,signal) %>%
      mutate(start_mut = trunc(start/seg.size)) %>%
      mutate_all(as.numeric)


    start.pos <- 0
    all.pos <- data.frame()

    while(start.pos < max(bed.sub$end)){
      ## chrom access ##
      temp <- bed.sub %>%
        filter(start >= start.pos & start <= (start.pos + seg.size))
      all.pos <- rbind(all.pos,c(start.pos,mean(temp$signal)))

      ## TMB ##


      #################################
      start.pos <- start.pos + seg.size
    }



  })))
}
