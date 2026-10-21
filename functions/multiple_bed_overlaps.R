#--------------------------------------------------------------------------------
## Function to find the overlap between 2 or more bed files

## Arguments:  tables: A list of all bed files to be used
#              names: a vector of column names
# Output:     binary_matrix: A binary matrix showing the overlaps
#             upset_plot: An upset Plot
#--------------------------------------------------------------------------------


## laoding libraries
  library(tidyverse)
  library(GenomicRanges)
  library(UpSetR)
  

## defining the function
  
  find_overlaps <- function(tables, names= NULL){
    
    gr_list <- lapply(tables, function(x){
      
      colnames(x)[1:3] <- c("chr", "start", "end")
      makeGRangesFromDataFrame(x, 
                               seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end",
                               keep.extra.columns = FALSE)
      
    })
    
    all_peaks <- do.call(c, gr_list)
    all_peaks <- GenomicRanges::reduce(all_peaks)
   
  
    overlap_list <- lapply(gr_list, function(x){
      as.integer(countOverlaps(all_peaks, x)) > 0
      
    })
    
    binary_mat <- data.frame(
      chr   = as.character(seqnames(all_peaks)),
      start = start(all_peaks),
      end   = end(all_peaks),
      do.call(cbind,overlap_list)
    )
    
    
    colnames(binary_mat)[4:ncol(binary_mat)] <- names
    binary_mat[4:ncol(binary_mat)] <-  binary_mat[4:ncol(binary_mat)] * 1
    
    
    binary_mat <- binary_mat %>%
      mutate(peak = paste(chr, start, end, sep = "_")) %>%
      select(-chr,-start, -end) %>%
      select(peak, everything())
    
    max_scale <- nrow(binary_mat)/ 0.9
    
    p <-  upset(binary_mat, order.by = "freq", set_size.show = TRUE,
                         matrix.color = "black", main.bar.color = "#065c8c",
                         sets.bar.color = "#065c8c", shade.color = "white", shade.alpha = 0.2,
                         point.size = 4, mb.ratio = c(0.75,0.25),
                         text.scale=c(1.5,3,1,1.25,1.5,2), set_size.scale_max = max_scale)
    
    
    
    return(list(binary_matrix= binary_mat,
                upset_plot= p))
    
  }
  
  
  