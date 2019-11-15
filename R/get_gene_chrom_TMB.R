#' get_gene_chrom_TMB
#'
#' Combines bed files with WGS/WXS/IMPACT genes
#' for chromatin accessibility and TMB
#'
#' @param bed BED file with fully processed with 7 columns
#' ("Chromosome","start","end","type","score","strand","signal").
#' @param gen complete maf file with chromosome, chromosome start/end, sample ID
#' @param map Start and end positions of all genes in the genome (or interest)
#' @param seg.size Size of segments to be used.
#' @return dat a data frame with the name of genes, their positions and their
#' estimated chromatin access and TMB.
#'
#'@export
#'

get_gene_chrom_TMB <- function(bed,gen,map,seg.size,cores=1){

  final <- data.frame()

  for(x in c(1:22,"X")){ # 1:22,"X"
    print(x)
    bed.sub <- bed[Chromosome == paste0("chr",x),c("start","end","signal")]
    gen.sub <- gen[chromosome == x,]
    map.sub <- map[chrom == paste0("chr",x),]
    print("Data processed")
    info <- data.table(do.call('rbind',lapply(1:nrow(map.sub),function(y,map.sub,bed.sub,gen.sub,seg.size){
      y <- map.sub[y,]
      # print(as.character(y[1]))
      range <- c(as.numeric(y[,3])-seg.size,as.numeric(y[,4])+seg.size)
      mean.bed <- mean(unlist(bed.sub[start >= range[1] & end <= range[2],"signal"]))
      tmb <- mean(unlist(gen.sub %>%
                           filter(chromosome_start >= range[1], chromosome_end <= range[2]) %>%
                           group_by(as.character(submitted_sample_id)) %>%
                           summarise(N = n()) %>%
                           select(N)))

      out <- c(unlist(y),mean.bed,tmb)
      out[3:6] <- as.numeric(out[3:6])
      # gc()
      return(out)
    },map = map.sub, bed = bed.sub, gen = gen.sub, seg.size = seg.size)))
    colnames(info)[5:6] <- c("ChromAccess","TMB")
    final <- rbind(final,info)
    gc()
  }


  # for(x in c(1:2)){ # 1:22,"X"
  #   print(x)
  #   bed.sub <- bed[Chromosome == paste0("chr",x),c("start","end","signal")]
  #   gen.sub <- gen[chromosome == x,]
  #   map.sub <- map[chrom == paste0("chr",x),]
  #   print("Data processed")
  #   cl <- makeCluster(2)
  #   info <- data.table(do.call('rbind',parLapply(cl,1:5,function(y,map.sub,bed.sub,gen.sub,seg.size){
  #     y <- map.sub[y,]
  #
  #     range <- c(as.numeric(y[,3])-seg.size,as.numeric(y[,4])+seg.size)
  #     mean.bed <- mean(unlist(bed.sub[start >= range[1] & end <= range[2],"signal"]))
  #     tmb <- mean(unlist(gen.sub %>%
  #                          filter(chromosome_start >= range[1], chromosome_end <= range[2]) %>%
  #                          group_by(as.character(submitted_sample_id)) %>%
  #                          summarise(N = n()) %>%
  #                          select(N)))
  #
  #     out <- c(unlist(y),mean.bed,tmb)
  #     out[3:6] <- as.numeric(out[3:6])
  #     # gc()
  #     return(out)
  #   },map = map.sub, bed = bed.sub, gen = gen.sub, seg.size = seg.size)))
  #   stopCluster(cl)
  #   colnames(info)[5:6] <- c("ChromAccess","TMB")
  #   final <- rbind(final,info)
  # }
  return(final)
}
