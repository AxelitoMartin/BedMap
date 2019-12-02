#' process_cancer
#'
#' For each chromosome create file that let's the user plot the chromosome
#' profile for the chromatin access and tumor mutation burden
#'
#' @param cancer name of the cancer to process
#' @param path path to the folder containing the data
#' @param output.path path to save output
#' @return dat chromosome, position and their
#' estimated mean chromatin access.
#'
#'@export
#'


process_cancer <- function(cancer,path,output.path){

  # chrom #
  files <- list.files(paste0(path,"/data/",cancer,"/Chrom/"))
  print(files)
  if(length(files) > 0){
    count <- 1
    for(i in files){
      bed <- readRDS(paste0(path,"data/",cancer,"/Chrom/",i))
      gc()
      if(i == files[1]){
        first <- bed_chrom(bed,seg.size=10^6)
        colnames(first)[match("DNaseI",colnames(first))] <- paste0("file",count)
      }
      else if(i == files[2]){
        chrom <- bed_chrom(bed,seg.size=10^6)
        colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
        temp <- full_join(first,chrom,
                          by = c("chrom","Position"))
      }
      else{
        chrom <- bed_chrom(bed,seg.size=10^6)
        colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
        temp <- full_join(temp,chrom ,by = c("chrom","Position"))
      }
      count <- count + 1
      gc()
    }
    if(length(files) > 1){
      temp <- temp %>%
        select(chrom,Position,paste0("file",1:length(files)))
      temp$Chromatin <- apply(temp,1,function(x){
        median(as.numeric(as.character(x[3:ncol(temp)])),rm.na=T)
      })
    }
    if(length(files) == 1){ colnames(temp)[ncol(temp)] <- "Chromatin"}
    Chrom <- temp %>% select(chrom,Position,Chromatin)
    print("Chrom done")
  }
  else{Chrom <- NULL}
  ######################################################
  # H3K4me1 #
  files <- list.files(paste0(path,"/data/",cancer,"/H3K4me1/"))
  print(files)
  count <- 1
  for(i in files){
    bed <- readRDS(paste0(path,"data/",cancer,"/H3K4me1/",i))
    gc()
    if(i == files[1]){
      first <- bed_chrom(bed,seg.size=10^6)
      colnames(first)[match("DNaseI",colnames(first))] <- paste0("file",count)
    }
    else if(i == files[2]){
      chrom <- bed_chrom(bed,seg.size=10^6)
      colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
      temp <- full_join(first,chrom,
                        by = c("chrom","Position"))
    }
    else{
      chrom <- bed_chrom(bed,seg.size=10^6)
      colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
      temp <- full_join(temp,chrom ,by = c("chrom","Position"))
    }
    count <- count + 1
    gc()
  }

  if(length(files) == 1) temp <- first
  if(length(files) > 1){
    temp <- temp %>%
      select(chrom,Position,paste0("file",1:length(files)))
    temp$H3K4me1 <- apply(temp,1,function(x){
      median(as.numeric(as.character(x[3:ncol(temp)])),rm.na=T)
    })
  }
  if(length(files) == 1){ colnames(temp)[ncol(temp)] <- "H4K3me1"}
  H3K4me1 <- temp %>% select(chrom,Position,H3K4me1)

  if(is.null(Chrom)) fulldat <- H3K4me1
  print("H3K4me1 done")

  #######################################################
  # H3K36me3 #
  files <- list.files(paste0(path,"/data/",cancer,"/H3K36me3/"))
  print(files)
  count <- 1
  for(i in files){
    bed <- readRDS(paste0(path,"data/",cancer,"/H3K36me3/",i))
    gc()
    if(i == files[1]){
      first <- bed_chrom(bed,seg.size=10^6)
      colnames(first)[match("DNaseI",colnames(first))] <- paste0("file",count)
    }
    else if(i == files[2]){
      chrom <- bed_chrom(bed,seg.size=10^6)
      colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
      temp <- full_join(first,chrom,
                        by = c("chrom","Position"))
    }
    else{
      chrom <- bed_chrom(bed,seg.size=10^6)
      colnames(chrom)[match("DNaseI",colnames(chrom))] <- paste0("file",count)
      temp <- full_join(temp,chrom ,by = c("chrom","Position"))
    }
    count <- count + 1
    gc()
  }

  if(length(files) == 1) temp <- first
  if(length(files) > 1){
    temp <- temp %>%
      select(chrom,Position,paste0("file",1:length(files)))
    temp$H3K36me3 <- apply(temp,1,function(x){
      median(as.numeric(as.character(x[3:ncol(temp)])),rm.na=T)
    })
  }
  if(length(files) == 1){ colnames(temp)[ncol(temp)] <- "H3K36me3"}
  H3K36me3 <- temp %>% select(chrom,Position,H3K36me3)
  fulldat <- full_join(fulldat,H3K36me3 ,by = c("chrom","Position"))
  print("H3K36me3 done")

  saveRDS(fulldat,file = paste0(output.path,cancer,"_mean.rds"))

}


