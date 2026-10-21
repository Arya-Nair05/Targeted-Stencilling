#--------------------------------------------------------------------------------
## Function to find the overlap between 2 or more bed files against a reference bed file

## Arguments:  reference: A reference bed file
#              tables: A list of all bed files to be used
#              names: a vector of column names
#--------------------------------------------------------------------------------


## laoding libraries
  library(tidyverse)
  library(GenomicRanges)


## defining the function

  find_overlaps_reference <- function(reference, tables, names= NULL){
    
    gr_list <- lapply(tables, function(x){
      
      colnames(x)[1:3] <- c("chr", "start", "end")
      makeGRangesFromDataFrame(x, 
                               seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end",
                               keep.extra.columns = FALSE)
      
    })
    
    colnames(reference)[1:3] <- c("chr", "start", "end")
    reference_GR <- makeGRangesFromDataFrame(reference, 
                                             seqnames.field = "chr", 
                                             start.field = "start", 
                                             end.field = "end")
    
    overlap_list <- lapply(gr_list, function(x){
      as.integer(countOverlaps(reference_GR, x)) > 0
      
    })
    
    binary_mat <- data.frame(
      chr   = as.character(seqnames(reference_GR)),
      start = start(reference_GR),
      end   = end(reference_GR),
      do.call(cbind,overlap_list)
    )
    

    colnames(binary_mat)[4:ncol(binary_mat)] <- names
    binary_mat[4:ncol(binary_mat)] <-  binary_mat[4:ncol(binary_mat)] * 1
    
    
    binary_mat <- binary_mat %>%
      mutate(peak = paste(chr, start, end, sep = "_")) %>%
      select(-chr,-start, -end) %>%
      select(peak, everything())
    
    return(binary_mat)
    
  }
  

 