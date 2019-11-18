#' gene_plot
#'
#' quick plots
#'
#' @param dat
#' @return dat a data frame with the name of genes, their positions and their
#' estimated chromatin access and TMB.
#'
#'@export
#'

gene_plot <- function(dat,lim.C,lim.T){
  final.all <- dat %>%
    tbl_df() %>%
    mutate(ChromAccess = -as.numeric(DNaseI),
           TMB = as.numeric(MutCount)) %>%
    group_by(gene) %>%
    summarise(ChromAccess = mean(ChromAccess),
              TMB = mean(TMB)) %>%
    ungroup()

  p1 <- final.all %>%
    ggplot(aes(x = ChromAccess,y = TMB,fill = chrom,color = chrom)) + geom_point() +  #fill = chrom,color = chrom
    theme(legend.position = "none") + ggtitle("Genome wide association per chromosome")

  fit <- lm(ChromAccess~TMB, data=final.all)
  corre <- cor(final.all %>% select(ChromAccess,TMB))
  p2 <- final.all %>%
    ggplot(aes(x = ChromAccess,y = TMB)) + geom_point() +
    theme(legend.position = "none") + geom_smooth(method = "lm")  +
    ggtitle(paste0("Gene association (R2 = ",round(as.numeric(summary(fit)[8]),digits = 3),
                   " cor = ",round(corre[1,2],digits = 3),")"))

  fit <- lm(ChromAccess~TMB, data=final.all%>%
              filter(ChromAccess < lim.C, TMB < lim.T))
  corre <- cor(final.all %>%
                 filter(ChromAccess < lim.C, TMB < lim.T) %>%
                 select(ChromAccess,TMB))
  p3 <- final.all %>%
    filter(ChromAccess < lim.C, TMB < lim.T) %>%
    ggplot(aes(x = ChromAccess,y = TMB)) + geom_point() +
    theme(legend.position = "none") + geom_smooth(method = "lm") +
    ggtitle(paste0("Gene association (R2 = ",round(as.numeric(summary(fit)[8]),digits = 3),
                   " cor = ",round(corre[1,2],digits = 3),")"))
  # plot_grid(p1, p2, p3,ncol = 1)

  return(list(p1=p1,p2=p2,p3=p3))
}


