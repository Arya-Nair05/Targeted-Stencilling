##-------------------------------------------------------------------------------------
##     Function to perform FDR correction 
##-------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(tibble)
  library(DiscreteFDR)

## Arguments:
#       1. file: Diff peaks smgstats diff --mods. Can be a file path or a data frame
#       2. fisher_p_val: string specifyting which specific fisher p value to correct for FDR
#.      3. discrete_fdr : A logical argument. If TRUE will do D-FDR on FIREs if not will do BH correction

## Function:

  diff_peaks_fdr_corr <- function(file= "file_path",
                                  fisher_p_val = "fisher_bpv",
                                  discrete_fdr = TRUE){
 
    
    ## checking the format of input file:
      if (class(file) == "character"){
        
        diff_df <- read.delim(file.path(file), header = TRUE, sep = '\t')
        diff_df <- diff_df %>%
          filter(rowSums(across(ends_with(c("countA","countB","totalA","totalB")))) != 0)
        
      } else if (class(file) == "data.frame"){
        
        diff_df <- file %>%
          filter(rowSums(across(ends_with(c("countA","countB","totalA","totalB")))) != 0)
        
      } else {
      
        print("Wrong format of file input")
        
      }
  
   
    ## finding the columns with the given input string of fisher:
      col_names <- colnames(diff_df)[grep(fisher_p_val, colnames(diff_df))]
      
      if (length(col_names)==0){
        stop("No matching fisher_p_val string. Check the input argument")
      }
    
    ## FDR correction on all Modifications
    # 6mA:
      p_val_col <- col_names[grep("mod6mA", col_names)]
      p_val_6ma <- diff_df[[p_val_col]]
      bh_6ma <- p.adjust(p_val_6ma, method = "BH")
      
      diff_df$bh_6mA <- bh_6ma
      
      
    # 5hmC:
      p_val_col <- col_names[grep("mod5mC", col_names)]
      p_val_5mC <- diff_df[[p_val_col]]
      bh_5mc <- p.adjust(p_val_5mC, method = "BH")
      
      diff_df$bh_5mC <- bh_5mc
      
    # 5hmC:
      p_val_col <- col_names[grep("mod5hmC", col_names)]
      p_val_5hmC <- diff_df[[p_val_col]]
      bh_5hmc <- p.adjust(p_val_5hmC, method = "BH")
      
      diff_df$bh_5hmC <- bh_5hmc
      
      
    ## D-FDR correction on FIRE:
      
      if (isTRUE(discrete_fdr)){
      fire_contingency <- diff_df[,grep("^fire_contingency", colnames(diff_df))] %>%
        setNames(c("a", "b", "c", "d")) %>%
        as.matrix()
      
      #D-FDR correction:
      fire_dfdr <- direct.discrete.BH(fire_contingency,
                                      "fisher",
                                      direction = "sd")
      
      fire_dfdr_summary <- summary(fire_dfdr)
      mod_dfdr_df <- fire_dfdr_summary$Table
      
      diff_df$Index <- 1:nrow(diff_df)
      diff_df <- left_join(diff_df, mod_dfdr_df[, c("Index", "Adjusted")], by = "Index")  
      
      } else {
        # Fire BH correction:
        p_val_col <- col_names[grep("fire", col_names)]
        p_val_fire <- diff_df[[p_val_col]]
        bh_fire <- p.adjust(p_val_fire, method = "BH")
        
        diff_df$bh_fire <- bh_fire
        
      }
      
      return(diff_df)
      
      
  }   
      
      
  
  
  