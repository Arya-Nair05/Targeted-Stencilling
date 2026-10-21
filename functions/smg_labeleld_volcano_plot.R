##-------------------------------------------------------------------------------------
## Function to make volcano plots from SMGstats diff peaks per modification:
## Giving a reference file-- marking points that overlap with the regions in the bed file
##-------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(ggplot2)
  library(tibble)
  library(GenomicRanges)
  library(patchwork)
  library(ggthemes)

  source("../../R_Projects/Functions/Function_Find_Overlaps.R")

## Arguments
#     1. file: Diff peaks smgstats diff --mods. Can be a file path or a data frame
#     2. bed _file: Bed file of regions that needs to be labelled. Can be a file path or a data frame
#     2. y axis data : For each modification we can specify which column to plot on the y-axis
#     3. y axis limits: for each modification we can set the y axis limits


## Example:
#      results <- diff_peaks_mods_labelled(file= diff_peaks_df,
#                               bed_file= ref_bed,
#                               y_axis_6ma = "BH_corrected_6ma",
#                               y_limits_6ma = c(0,50))

# access the results: 
#        results$volcano_plots


## function:

  diff_peaks_mods_labelled <- function(file= "file_path",
                                    bed_file= "file_path",
                                    label_names= c("diff_peak", "bed_peak"),
                                    y_axis_6ma = "mod6mA_fisher_bpv",
                                    y_axis_fire = "fire_fisher_bpv",
                                    y_axis_5mc = "mod5mC_fisher_bpv",
                                    y_axis_5hmc = "mod5hmC_fisher_bpv",
                                    y_limits_6ma = c(0,0),
                                    y_limits_fire = c(0,0),
                                    y_limits_5mc = c(0,0),
                                    y_limits_5hmc = c(0,0)){
    
    
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
    
    ## checking the format of the input reference bed file:
      if (class(bed_file) == "character"){
        
        ref_bed <- read.delim(file.path(bed_file), header = FALSE, sep = "\t")
        ref_bed <- ref_bed[,1:3]
        
      } else if (class(bed_file) == "data.frame"){
        
        ref_bed <- bed_file[,1:3]
      } else {
        
        print("Wrong format of bed_file input")
      }
      
    
    
    ## finding overlaps of diff peaks with reference peaks
      peaks_overlap <-  bed_overlaps(bed.1 = diff_df[,1:3],
                                     bed.2 = ref_bed,
                                     names = label_names)
      
      col_name <- label_names[2]
      
      diff_df <-  diff_df %>%
        mutate(peak= paste(chrom, start, stop, sep = "-"))
      
      diff_df[col_name] <- diff_df$peak %in% peaks_overlap[[1]]
      
   
    ## Volcano plots per modifications:
    # 6mA
      p1 <- ggplot(diff_df, 
                   aes(x= -log2(mod6mA_or), y= -log10(.data[[y_axis_6ma]]), colour = .data[[col_name]]))+
        geom_point(size=1, alpha=0.8)+
        scale_color_manual(name = col_name, 
                           values = c("FALSE" = "#69a2c2", "TRUE" = "#D55E00"))+
        theme_few(base_size = 12)+
        geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
        geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
        labs(title= "6mA",
             x = bquote(Log[2] ~ "Odds Ratio (6mA)"))+
        theme(axis.title.y = element_text(size = 10))
    
    
    # fire:
      p2 <- ggplot(diff_df, 
                   aes(x= -log2(fire_or), y= -log10(.data[[y_axis_fire]]),colour = .data[[col_name]]))+
        geom_point(size=1, alpha=0.8)+
        scale_color_manual(name = col_name, 
                           values = c("FALSE" = "#69a2c2", "TRUE" = "#D55E00"))+
        theme_few(base_size = 12)+
        geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
        geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
        labs(title= "FIRE",
             x = bquote(Log[2] ~ "Odds Ratio (FIRE)"))+
        theme(axis.title.y = element_text(size = 10))
    
    # 5mC
      p3 <- ggplot(diff_df,
                   aes(x= -log2(mod5mC_or), y= -log10(.data[[y_axis_5mc]]),colour = .data[[col_name]]))+
        geom_point(size=1, alpha=0.8)+
        scale_color_manual(name = col_name, 
                           values = c("FALSE" = "#69a2c2", "TRUE" = "#D55E00"))+
        theme_few(base_size = 12)+
        geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
        geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
        labs(title= "5mC",
             x = bquote(Log[2] ~ "Odds Ratio (5mC)"))+
        theme(axis.title.y = element_text(size = 10))
    
    # 5hmC
      p4 <- ggplot(diff_df,
                   aes(x= -log2(mod5hmC_or), y= -log10(.data[[y_axis_5hmc]]),colour = .data[[col_name]]))+
        geom_point(size=1, alpha=0.8)+
        scale_color_manual(name = col_name, 
                           values = c("FALSE" = "#69a2c2", "TRUE" = "#D55E00"))+
        theme_few(base_size = 12)+
        geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
        geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
        labs(title= "5hmC",
             x = bquote(Log[2] ~ "Odds Ratio (5hmC)"))+
        theme(axis.title.y = element_text(size = 10))
    
    ## checking if the y limits argument has been provided:
    
      if (max(y_limits_6ma) != 0){
        p1 <- p1 + ylim(y_limits_6ma) + ylab( bquote(-Log[10] ~ y_axis_6ma))
      } else {
        p1 <-  p1 + ylab(bquote(-log[10] ~ .(y_axis_6ma)))
      }
      
      if (max(y_limits_fire) != 0){
        p2 <- p2 + ylim(y_limits_fire)+ ylab( bquote(-Log[10] ~ y_axis_fire))
      } else {
        p2 <-  p2 + ylab( bquote(-Log[10] ~ .(y_axis_fire)))
      }
      
      if (max(y_limits_5mc) != 0){
        p3 <- p3 + ylim(y_limits_5mc) + ylab( bquote(-Log[10] ~ y_axis_5mc))
      } else {
        p3 <-  p3 + ylab( bquote(-Log[10] ~ .(y_axis_5mc)))
      }
      
      if (max(y_limits_5hmc) != 0){
        p4 <- p4 + ylim(y_limits_5hmc) + ylab( bquote(-Log[10] ~ y_axis_5hmc))
      } else {
        p4 <-  p4 + ylab( bquote(-Log[10] ~ .(y_axis_5hmc)))
      }
      
    
     plots <- list(p1,p2,p3,p4)
     final_plots <-  wrap_plots(plots, ncol = 2, nrow = 2)+
        plot_layout(axes = "collect", axis_titles = "collect_y", guides = "collect")
      
    
      return(list(volcano_plots = final_plots,
        diff_peaks_anno = diff_df,
        mod6ma_plot= p1,
        fire_plot= p2,
        mod5mc_plot= p3,
        mod5hmc_plot= p4))
   
  }
