##-------------------------------------------------------------------------------------
## Function to make a bar plot and upset plot of significant regions
## The p value thresholds can be specified
##-------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(ggplot2)
  library(tibble)
  library(GenomicRanges)
  library(patchwork)
  library(UpSetR)
  library(ggplotify)
  library(ggthemes)


##Arguments:
#     1.file: Data frame with FDR corrected values. Output from diff_peaks_fdr_corr() function
#     2. mods_thresh : Threshold to select regions with significant modification changes
#     3. fire_thresh: Threshold to select regions with signficant FIRE differences
#     4. discrete_dfdr: Logical value denoting if D-FDR was run on FIRE regions
#     5. mod6mA_dfdr: Logical value denoting if D-FDR was run on mod6mA regions

#     The file should have columns named chrom, start and stop

## Function:
  filter_and_plot_mods <- function(file_df, 
                                   mods_thresh = 1e-20,
                                   fire_thresh= 0.01,
                                   discrete_fdr= TRUE,
                                   mod6mA_dfdr= TRUE) {
  
    
    ## columns denoting significant modifications
    # 5mC
      file_df$mod5mC <- file_df$bh_5mC < mods_thresh
    
    # 5hmC
      file_df$mod5hmC <- file_df$bh_5hmC < mods_thresh
    
    # mod6mA
      if(isTRUE(mod6mA_dfdr)){
        
        file_df$mod6mA <- file_df$mod6mA_Adjusted < mods_thresh
        
      } else {
        
        file_df$mod6mA <- file_df$bh_6mA < mods_thresh
        
      }
    
    # fire
      if(isTRUE(discrete_fdr)){
        
        file_df$Fire <- file_df$Adjusted < fire_thresh
        
      } else {
        
        file_df$Fire <- file_df$bh_fire < fire_thresh
        
      }
    
   
    ## making a frequency table
      sig_df <- file_df %>%
        select(peak, mod6mA,mod5mC,mod5hmC,Fire) %>%
        pivot_longer(cols = c(mod6mA,mod5mC,mod5hmC,Fire), names_to = "mods", values_to = "is_sig") %>%
        mutate(mods= sub("mod", "", mods))
      
      num_sig <- data.frame(table(mods= sig_df$mods, sig= sig_df$is_sig)) %>%
        filter(sig == TRUE) 
    
    ## Generating a bar plot
      p <-  ggplot(num_sig, aes(x = mods, y = Freq)) +
        geom_col(fill = "#065c8c",color = "black", width = 0.35,linewidth = 0.3) +
        geom_text(aes(label = Freq), vjust = -0.4, size = 4.5) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs( title = "Number of significant peaks per modification", 
              x = "Modification",
              y= "Frequency") +
        theme_few(base_size = 14) +
        theme(axis.text = element_text(size = 12),
              plot.title = element_text(face = "bold",size = 14),
              panel.border = element_rect(colour = "black",linewidth = 0.5))
      
    
    ## Upset plot
      upset_df <- file_df[,c("peak","mod6mA","Fire","mod5mC","mod5hmC")] %>%
        column_to_rownames("peak")
      upset_df <- upset_df *1
     
      # setting the max scale of the sets:
      set_max_scale <- max(num_sig$Freq) * 1.3
      
      ## upset plot cannot be saved to an obj if we use a matrix for plotting
      ## instead use a data frame as the initial input
      
      upset_plot <-  upset(upset_df, order.by = "freq", set_size.show = TRUE,
                           matrix.color = "black", main.bar.color = "#065c8c",
                           sets.bar.color = "#065c8c", shade.color = "steelblue", shade.alpha = 0.2,
                           point.size = 2.5, mb.ratio = c(0.75,0.25),
                           text.scale=c(1.5,1.5,1,1.25,1.5,2), set_size.scale_max = set_max_scale)
      
    ## defining return objectes 
      return(list(all_sig_mods = file_df,
                  upset_df = upset_df,
                  bar_plot = p,
                  upset_plot = upset_plot))
    
    
  }
  
  
  
 
  
 
  