##-------------------------------------------------------------------------------
##. Function to find and add modification rates to a given data frame
## Rate = read counts of a mod in group A/ total counts of the mod in group A
##-------------------------------------------------------------------------------

## Library:
  library(tidyverse)
  library(tibble)

## Arguments:
#     1. Data frame output from smgstats diff function


## function:

  add_mod_rate <- function(file= data_frame){
    
    
  #mod6ma rate:
    # group A
    file$mod6mA_rateA <- file$mod6mA_countA/ (file$mod6mA_totalA +1e-300)
    
    # group B
    file$mod6mA_rateB <- file$mod6mA_countB/ (file$mod6mA_totalB +1e-300)
    
  # fire rate
    # group A
    file$fire_rateA <- file$fire_countA/ (file$fire_totalA +1e-300)
    
    # group B
    file$fire_rateB <- file$fire_countB/ (file$fire_totalB +1e-300)
    
  # 5mC rate
    # group A
    file$mod5mC_rateA <- file$mod5mC_countA/ (file$mod5mC_totalA +1e-300)
    
    # group B
    file$mod5mC_rateB <- file$mod5mC_countB/ (file$mod5mC_totalB +1e-300)
    
    
  # 5hmC rate
    # group A
    file$mod5hmC_rateA <- file$mod5hmC_countA/ (file$mod5hmC_totalA +1e-300)
    
    # group B
    file$mod5hmC_rateB <- file$mod5hmC_countB/ (file$mod5hmC_totalB +1e-300)
    
    
    return(file)
    
  }