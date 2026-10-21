#--------------------------------------------------------------------------------------------
## Function to find the overlap between 2  bed files and return the overlapping coords

## Arguments:  bed_1, bed_2: 2 BED files 
#              names:        a vector of 2 column names corresponding to each of the BED
#              overlap_type: type argument of findOverlapPairs function
#--------------------------------------------------------------------------------------------


## laoding libraries
  library(tidyverse)
  library(GenomicRanges)


## defining the function

bed_overlaps <- function(bed.1= bed_1,
                           bed.2=bed_2,
                           names= NULL,
                           overlap_type= "any"){
  
    tables <- list(bed.1, bed.2)
    gr_list <- lapply(tables, function(x){
      
      colnames(x)[1:3] <- c("chr", "start", "end")
      makeGRangesFromDataFrame(x, 
                               seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end",
                               keep.extra.columns = FALSE)
      
    })
    
    overlaps <- findOverlapPairs(gr_list[[1]], gr_list[[2]], type= overlap_type)
      
    first_gr <- S4Vectors::first(overlaps)
    second_gr <- S4Vectors::second(overlaps)
    
    df <- data.frame(first=
                       paste(seqnames(first_gr),
                             start(first_gr),
                             end(first_gr),
                             sep = "-"),
                     second=
                       paste(seqnames(second_gr),
                              start(second_gr),
                              end(second_gr),
                              sep="-"))
    
    colnames(df) <- names
    
    return(df)
    
}
  
  
