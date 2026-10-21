##-------------------------------------------------------------------------------------
## Function to make volcano plots from SMGstats diff peaks per modification:
##-------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(ggplot2)
  library(tibble)
  library(GenomicRanges)
  library(patchwork)
  library(ggthemes)

## Arguments
#     1. file: Can be a Path to a file that has the output from smgstats diff --mods or can be a data frame
#     2. y axis data : For each modification we can specify which column to plot on the y-axis
#     3. y axis limits: for each modification we can set the y axis limits

## Example:
#      diff_peaks_mods_plots(file= diff_peaks_df,
#                            y_axis_6ma = "BH_corrected_6ma",
#                            y_limits_6ma = c(0,50))

## function:

  diff_peaks_mods_plots <- function(file= "file_path",
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
  
    
    
    ## Volcano plots per modifications:
    # 6mA
    p1 <- ggplot(diff_df, aes(x= -log2(mod6mA_or), y= -log10(.data[[y_axis_6ma]])))+
      geom_point(size=0.8, alpha=0.6, color="#0072b2")+
      theme_few(base_size = 12)+
      geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
      geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
      labs(title= "6mA",
           x = bquote(Log[2] ~ "Odds Ratio (6mA)"))+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    
    # fire:
    p2 <- ggplot(diff_df, aes(x= -log2(fire_or), y= -log10(.data[[y_axis_fire]])))+
      geom_point(size=0.8, alpha=0.6,color="#0072b2")+
      theme_few(base_size = 12)+
      geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
      geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
      labs(title= "FIRE",
           x = bquote(log[2] ~ "Odds Ratio (FIRE)"))+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    # 5mC
    p3 <- ggplot(diff_df, aes(x= -log2(mod5mC_or), y= -log10(.data[[y_axis_5mc]])))+
      geom_point(size=0.8, alpha=0.6,color="#0072b2")+
      theme_few(base_size = 12)+
      geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
      geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
      labs(title= "5mC",
           x = bquote(log[2] ~ "Odds Ratio (5mC)"))+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    # 5hmC
    p4 <- ggplot(diff_df, aes(x= -log2(mod5hmC_or), y= -log10(.data[[y_axis_5hmc]])))+
      geom_point(size=0.8, alpha=0.6,color="#0072b2")+
      theme_few(base_size = 12)+
      geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
      geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
      labs(title= "5hmC",
           x = bquote(log[2] ~ "Odds Ratio (5hmC)"))+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    
    
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
    wrap_plots(plots, ncol = 2, nrow = 2)+
      plot_layout( guides = "collect")
    

  }
  