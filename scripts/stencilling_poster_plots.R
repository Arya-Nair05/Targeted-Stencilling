##-----------------------------------------------------------------------------------------------------------------
## Understanding changes in 5mC DNA methylation at MER11s upon the loss of ZNF808
## Script to generate plots used for the poster:
## Dissecting ZNF808 and MER11 elements during pancreatic development using targeted chromatin stencilling
##-----------------------------------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(ggplot2)
  library(tibble)
  library(GenomicRanges)
  library(DiscreteFDR)
  library(UpSetR)
  library(grid)
  library(DiscreteFDR)
  library(pheatmap)
  library(parallel)
  library(sparseMatrixStats)

## function
  source("functions/bed_overlaps.R")
  source("functions/multiple_bed_overlaps.R")
  source("functions/multi_overlaps_with_reference.R")
  source("functions/smg_volcano_plot.R")
  source("functions/smg_labeleld_volcano_plot.R")
  source("functions/smg_FDR_correction.R")
  source("functions/add_modification_rate.R")
  source("functions/filter_sig_and_plot.R")

  set.seed(1234)

## loading MER11 BED file:
  # all MER11
  all_mer11 <- read.delim("data/all_mer11_hg38.bed", header = FALSE)
  chr_names <- c(paste0("chr", 1:22), "chrX", "chrY")
  all_mer11 <- all_mer11[all_mer11$V1 %in%chr_names,]
  
  # A subset of 220 MER11 that lost H3K9me3 and gained H3K27ac in KO
  mer11_220 <- read.delim("data/mer11_220_clusters.hg19lohg38.bed", header = FALSE)
  mer11_220 <- mer11_220[,1:3]

# mer11_220 overlap
  mer11_overlap <- bed_overlaps(bed.1 = all_mer11,
                                bed.2 = mer11_220,
                                names = c("mer11", "mer11_220"))
  
## loading the DA peaks from SMGstats:
  mods_df <- read.delim("data/znf808_diff_access_with_mods_CORRECTED.tsv", sep = "\t") 
  
## annotating the DA peaks as mer11/non-mer11
  peaks_mer11 <- bed_overlaps(bed.1 = mods_df[,1:3],
                              bed.2 = all_mer11,
                              names = c("peaks", "mer11"))
  
  peaks_mer11 <- left_join(peaks_mer11, mer11_overlap, by = "mer11")
  peaks_mer11_220 <- peaks_mer11[!is.na(peaks_mer11$mer11_220),]
  
## annotating the main DA peaks data frame if it's mer11 or not
  mods_df <- mods_df %>%
    mutate(peak= paste(chrom, start, stop, sep = "-"),
           is_mer11= peak %in% peaks_mer11$peaks,
           is_mer11_220= peak %in% peaks_mer11_220$peaks,
           category= case_when(is_mer11_220  ~ "mer11_220",
                               is_mer11 ~ "mer11", 
                               .default = "non_mer11"))

## volcano plots:
  # All modifacations and FIRE
  all_volcano <- diff_peaks_mods_plots(file = mods_df)
  all_volcano +
    plot_annotation(title = "Volcano plots for all mods and FIRE",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
  
  
## FDR correction: Benjamini Hochberg for 5mC/6mA/5hmC modifications and Discrete-FDR for FIRE
  mods_df <- diff_peaks_fdr_corr(mods_df, 
                                 discrete_fdr = TRUE)
  
## volcano plots with corrected p values:
  fdr_correted <- diff_peaks_mods_plots(file = mods_df,
                        y_axis_6ma = "bh_6mA",
                        y_axis_fire = "Adjusted",
                        y_axis_5mc = "bh_5mC",
                        y_axis_5hmc = "bh_5hmC")
  
  fdr_correted+
    plot_annotation(title = "Volcano plots for all mods and FIRE-after FDR corrrection",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))

## labelled plot: All MER11
  fdr_labelled_plot <- diff_peaks_mods_labelled(file = mods_df,
                                                bed_file = all_mer11,
                                                label_names = c("peaks", "Mer11"),
                                                y_axis_6ma = "bh_6mA",
                                                y_axis_fire = "Adjusted",
                                                y_axis_5mc = "bh_5mC",
                                                y_axis_5hmc = "bh_5hmC")
  
  fdr_labelled_plot$volcano_plots + plot_annotation(title = "After FDR: All MER11 marked",
                                                    theme = theme(plot.title = 
                                                                    element_text(hjust = 0.5, face = "bold", size = 16)))
  
  # Mer11 220
  fdr_labelled_plot_2 <- diff_peaks_mods_labelled(file = mods_df,
                                                bed_file = mer11_220,
                                                label_names = c("peaks", "Mer11"),
                                                y_axis_6ma = "bh_6mA",
                                                y_axis_fire = "Adjusted",
                                                y_axis_5mc = "bh_5mC",
                                                y_axis_5hmc = "bh_5hmC")
  
  fdr_labelled_plot_2$volcano_plots + plot_annotation(title = "After FDR: 220 MER11 marked",
                                                    theme = theme(plot.title = 
                                                                    element_text(hjust = 0.5, face = "bold", size = 16)))
  
## Finding the peaks/regions with significant mods /fire:
  mods_thresh <- 1e-4
  fire_thresh <- 0.01
  
  sig_mods <- filter_and_plot_mods(file = mods_df,
                                   mods_thresh = mods_thresh,
                                   fire_thresh = fire_thresh,
                                   mod6mA_dfdr = FALSE)
  
  ## upset plot:
  sig_mods$upset_plot
  grid.text("Upset plot of significant regions per modification", just = "top", 
            x = 0.65, y = 0.98, gp = gpar(fontsize=15))
  
  ## bar plot:
  sig_mods$bar_plot
  
## finding peaks with higher accessibility in KO (in our case -ve odds ratio since first group was WT)
  all_sig_df <- sig_mods$all_sig_mods
  mod_names <- c("mod6mA", "mod5mC", "mod5hmC")
  
  mods_thresh <- 1e-3
  log2_or_thresh <- 0.1
  
## adding a new column for each of the modification
  for (i in mod_names) {
    
    bh_col <- paste0("bh_",sub("mod", "",i))
    or_col <- paste0(i, "_or")
    
    log2_or_col <- paste0("log2_", or_col)
    new_col <- paste0(i, "_sig_gain")
    
    all_sig_df[[log2_or_col]] <- -log2(all_sig_df[[or_col]])
    all_sig_df[[new_col]] <- all_sig_df[[bh_col]] < mods_thresh & all_sig_df[[log2_or_col]] > log2_or_thresh
    
  }
  
## adding sig gain column for the fire:
  all_sig_df$log2_fire_or <- -log2(all_sig_df$fire_or)
  all_sig_df$fire_sig_gain <- all_sig_df$Adjusted < fire_thresh & all_sig_df$log2_fire_or > log2_or_thresh
  
## plotting bar plot and upset plot for these sig gain peaks:
  sig_gain_peaks <- all_sig_df %>%
    filter(if_any(contains("sig_gain"), ~. ==TRUE))
  
  mods_thresh <- 1e-4
  fire_thresh <- 0.01
  
  sig_mods_gain <- filter_and_plot_mods(file = sig_gain_peaks,
                                   mods_thresh = mods_thresh,
                                   fire_thresh = fire_thresh)
  
  
  sig_mods_gain$bar_plot
  sig_mods_gain$upset_plot
  grid.text("Upset plot for all peaks gaining accessibility in the KO", just = "top", 
            x= 0.65, y = 0.98, gp = gpar(fontsize=15))
  
  upset_df <- sig_mods_gain$upset_df
  
  
##-------------------------------------------------------------------------
## D-FDR correction for the 6mA modified regions 
##-------------------------------------------------------------------------
## getting the contingency table
  mod6ma_mat <- mods_df[, grep("^mod6mA_contingency", colnames(mods_df))] %>%
    setNames(c("a", "b", "c", "d")) %>%
    as.matrix()
  
## D-FDR correction:ao
  mod_dfdr <- direct.discrete.BH(mod6ma_mat,
                                  "fisher",
                                  direction = "sd")
    
  mod_dfdr_summary <- summary(mod_dfdr)
  mod_dfdr_df <- mod_dfdr_summary$Table
  
  colnames(mod_dfdr_df)[3] <- "mod6mA_Adjusted"
  
## merging this table to the data frame with all the mod regions:
  mods_df <- left_join(mods_df,mod_dfdr_df[, c("Index", "mod6mA_Adjusted")], by = "Index" )
  
## volcano plots using D-FDR corrected 6mA regions:
  fdr_labelled_plot <- diff_peaks_mods_labelled(file = mods_df,
                                                bed_file = all_mer11,
                                                label_names = c("peaks", "Mer11"),
                                                y_axis_6ma = "mod6mA_Adjusted",
                                                y_axis_fire = "Adjusted",
                                                y_axis_5mc = "bh_5mC",
                                                y_axis_5hmc = "bh_5hmC")
  
  fdr_labelled_plot$volcano_plots+
    plot_annotation(title = "After D-FDR on 6mA & FIRE",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))
  
  
  
## Extracting individual plots for poster:
 # fdr_labelled_plot$mod6ma_plot
  
  labelled_df <- fdr_labelled_plot$diff_peaks_anno
  
  # 6mA
  ggplot(labelled_df, 
         aes(x= -log2(mod6mA_or), y= -log10(mod6mA_Adjusted), colour = Mer11))+
    geom_point(size=2, alpha=0.8)+
    scale_color_manual(name = "Mer11", 
                       values = c("FALSE" = "grey40", "TRUE" = "#D55E00"))+
    theme_few(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
    geom_vline(xintercept = c(-0.5,0.5), linetype="dashed", color="grey40",linewidth=0.4)+
    labs(title= "6mA",
         x = bquote(Log[2] ~ "Odds Ratio (6mA)"),
         y= bquote(-Log[10] ~"Adj p value"))+
    theme(axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(colour = "black", linewidth = 0.4))
  
  # 5mC
  labelled_df <-  labelled_df %>%
    arrange(Mer11)
    
  ggplot(labelled_df, 
         aes(x= -log2(mod5mC_or), y= -log10(bh_5mC), colour = Mer11))+
    geom_point(size=2.8, alpha=0.8)+
    scale_color_manual(name = "Mer11", 
                       values = c("FALSE" = "grey40", "TRUE" = "#D55E00"))+
    theme_few(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
    geom_vline(xintercept = c(-0.5,0.5), linetype="dashed", color="grey40",linewidth=0.4)+
    labs(title= "5mC",
         x = bquote(Log[2] ~ "Odds Ratio (5mC)"),
         y= bquote(-Log[10] ~"Adj p value"))+
    theme(axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(colour = "black", linewidth = 0.4))
  
  
  # fire
  ggplot(labelled_df, 
         aes(x= -log2(fire_or), y= -log10(Adjusted), colour = Mer11))+
    geom_point(size=2.8, alpha=0.8)+
    scale_color_manual(name = "Mer11", 
                       values = c("FALSE" = "steelblue", "TRUE" = "#D55E00"))+
    theme_few(base_size = 14)+
    geom_hline(yintercept = 0, linetype="dashed", color="grey40", linewidth=0.4)+
    geom_vline(xintercept = c(-0.1,0.1), linetype="dashed", color="grey40",linewidth=0.4)+
    labs(title= "FIRE",
         x = bquote(Log[2] ~ "Odds Ratio (Fire)"),
         y= bquote(-Log[10] ~"Adj p value"))+
    theme(axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(colour = "black", linewidth = 0.4))
  
## finding the number of MER11 regions that show sig gain of accessibility/methylation
  # labels on volcano plot:
  
  #5mC:
  length(unique(all_sig_df$mer11[all_sig_df$log2_mod5mC_or > 0.5 & all_sig_df$bh_5mC < 1e-4 & all_sig_df$is_mer11 == TRUE]))
  length(unique(all_sig_df$mer11[all_sig_df$log2_mod5mC_or< -0.5 & all_sig_df$bh_5mC < 1e-4 & all_sig_df$is_mer11 == TRUE]))
  
  #6mA
  length(unique(all_sig_df$mer11[all_sig_df$log2_mod6mA_or< -0.5 & all_sig_df$mod6mA_Adjusted < 1e-4 & all_sig_df$is_mer11 == TRUE]))
  length(unique(all_sig_df$mer11[all_sig_df$log2_mod6mA_or > 0.5 & all_sig_df$mod6mA_Adjusted < 1e-4 & all_sig_df$is_mer11 == TRUE]))
  
  
  
## finding regions with a significant change of modifcations:
  mods_thresh <- 1e-4
  fire_thresh <- 0.01
  
  sig_mods <- filter_and_plot_mods(file = mods_df,
                                   mods_thresh = mods_thresh,
                                   fire_thresh = fire_thresh,
                                   discrete_fdr = TRUE,
                                   mod6mA_dfdr= TRUE)
  
  ## upset plot:
  sig_mods$upset_plot
  grid.text("Upset plot of significant regions per modification", just = "top", 
            x = 0.65, y = 0.98, gp = gpar(fontsize=15))
  
  ## bar plot:
  sig_mods$bar_plot
  
## finding peaks with higher accessibility in KO (in our case -ve log2 (odds ratio) since first group was WT)
  all_sig_df <- sig_mods$all_sig_mods
  mod_names <- c("mod5mC", "mod5hmC")
  
  mods_thresh <- 1e-4
  log2_or_thresh <- 0.5
  
  ## adding a new column for each of the modification
    for (i in mod_names) {
      
      bh_col <- paste0("bh_",sub("mod", "",i))
      or_col <- paste0(i, "_or")
      
      log2_or_col <- paste0("log2_", or_col)
      new_col <- paste0(i, "_sig_gain")
      
      all_sig_df[[log2_or_col]] <- -log2(all_sig_df[[or_col]])
      all_sig_df[[new_col]] <- all_sig_df[[bh_col]] < mods_thresh & all_sig_df[[log2_or_col]] > log2_or_thresh
      
    }
  
## adding sig gain column for the fire and mod6mA:
  # fIRE:
    all_sig_df$log2_fire_or <- -log2(all_sig_df$fire_or)
    all_sig_df$fire_sig_gain <- all_sig_df$Adjusted < fire_thresh & all_sig_df$log2_fire_or > log2_or_thresh
    
  # mod6mA:
    all_sig_df$log2_mod6mA_or <- -log2(all_sig_df$mod6mA_or)
    all_sig_df$mod6mA_sig_gain <- all_sig_df$mod6mA_Adjusted < mods_thresh & all_sig_df$log2_mod6mA_or > log2_or_thresh
    
## plotting bar plot and upset plot for these sig gain peaks:
  sig_gain_peaks <- all_sig_df %>%
    filter(if_any(contains("sig_gain"), ~. ==TRUE)) %>%
    select(peak, colnames(.)[grep("sig_gain", colnames(.))]) %>%
    column_to_rownames("peak")
      
  sig_gain_peaks <- sig_gain_peaks *1
  set_max_scale <- nrow(sig_gain_peaks)*1.1
  
  upset(sig_gain_peaks, order.by = "freq", set_size.show = TRUE,
                       matrix.color = "black", main.bar.color = "#065c8c",
                       sets.bar.color = "#065c8c", shade.color = "white", shade.alpha = 0.2,
                       point.size = 4, mb.ratio = c(0.75,0.25),
                       text.scale=c(1.5,1.5,1,1.25,1.5,2), set_size.scale_max = set_max_scale)

  grid.text("Upset plot for all peaks gaining accessibility in the KO", just = "top", 
              x= 0.65, y = 0.98, gp = gpar(fontsize=15))
  
   
## adding the MER11 coordinates to the table
  all_sig_df <- left_join(all_sig_df, peaks_mer11, by = c("peak"="peaks"))
  all_sig_df <-  left_join(all_sig_df, mer11_overlap, by = "mer11")
##----------------------------------------------------------------------------------------------------------
# Few stats about the data
##----------------------------------------------------------------------------------------------------------
# how many sig mod6mA are MER11?
  # 311 (175 from mer11 220) (dfdr pval- 1e-3 and log2: 0.1 thresh)
  # 288 (167 from mer11 220) (dfdr pval 1e-4 and log2: 0.5 thresh)
  length(unique(all_sig_df[all_sig_df$mod6mA_sig_gain==TRUE &all_sig_df$is_mer11==TRUE,]$mer11))
  length(unique(all_sig_df[all_sig_df$mod6mA_sig_gain==TRUE &all_sig_df$is_mer11==TRUE,]$mer11_220))
  
# how many regions show significant loss of 5mC
  # 1341 (< -0.1)
  # 1173 (< -0.5)
  length(which(all_sig_df$bh_5mC < mods_thresh & all_sig_df$log2_mod5mC_or < -0.5))
  
# how many of these regions losing 5mC are MER11
  # 218 (149 from mer11 2220)- < -0.1
  # 193 (143 from mer11 220)- < -0.5
  length(unique(all_sig_df[all_sig_df$bh_5mC < mods_thresh & all_sig_df$log2_mod5mC_or < -0.5,]$mer11))
  length(unique(all_sig_df[all_sig_df$bh_5mC < mods_thresh & all_sig_df$log2_mod5mC_or < -0.5,]$mer11_220))
  
  mer_220_meth_loss_no_access <- all_sig_df %>%
    filter(mod6mA_sig_gain==TRUE & is_mer11==TRUE) %>%
    distinct(mer11) %>%
    drop_na()
  
  all_sig_df %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    distinct(mer11) %>%
    drop_na() %>%
    nrow()
  
  all_sig_df %>%
    filter(mod6mA_Adjusted < mods_thresh & log2_mod6mA_or > 0.5) %>%
    distinct(mer11) %>%
    drop_na() %>%
    nrow()
  
trial <- all_sig_df %>%
    filter(mod6mA_Adjusted < mods_thresh & log2_mod6mA_or > 0.5) %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    filter(!is.na(mer11))
  
# how many of the MER11 gaining accessibility lose methylation
 # 208 (146 from mer11 220) : < -0.1
  # 183 (136 from mer11 220): < -0.5

  all_sig_df %>%
    filter(mod6mA_sig_gain==TRUE & is_mer11==TRUE) %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    distinct(mer11) %>%
    drop_na() %>%
    nrow()
  
  no_access_meth_loss <- all_sig_df %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    filter(mod6mA_sig_gain==FALSE & is_mer11==TRUE) %>%
    distinct(mer11) %>%
    drop_na()
  
  mer_220peakmer_220_meth_loss <- all_sig_df %>%
    filter(mod6mA_sig_gain==TRUE & is_mer11==TRUE) %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    distinct(mer11) %>%
    drop_na()
  
## how many regions lose methylation but do not gain accessibility:
  # 12 MER11---however these are the mer11 that has multiple (2) sten peaks overlapping with the MER11
  # So one of it shows sig 6mA and sig 5mC while the other show only sig 5mC loss
  
  no_access_meth_loss <- all_sig_df %>%
    filter(bh_5mC < mods_thresh & log2_mod5mC_or < -0.5) %>%
    filter(mod6mA_sig_gain==FALSE) %>%
    distinct(mer11) %>%
    drop_na()

  mer_220_only_meth_loss <- anti_join(mer_220_meth_loss_no_access,mer_220_meth_loss)
  
#--------------------------------------------------------------------------------------------------
# Comparison with scATAC markers:
#--------------------------------------------------------------------------------------------------
## loading scATAC markers:
  sc_atac <- read.csv("SC_KO_vs_WT_markers.csv")
  
## filtering significant scATAC markers:
  padj <- 1e-3
  logfc <- 0.2
  
  sc_atac_peaks <- sc_atac %>%
    separate(X, into = c("chrom", "start", "end"), sep = "-") %>%
    filter(p_val_adj < padj & avg_log2FC > logfc) %>%
    select(c("chrom", "start", "end"))
  
  write.table(sc_atac_peaks, "SC_atac_KO_gained_markers.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
## finding regions from stencilling with significant mod6mA
  da_peaks_sten <- sig_gain_peaks %>%
    filter(mod6mA_sig_gain==1) %>%
    rownames_to_column("peak") %>%
    select(peak) %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-")
  
  
## finding overlaps:
  sten_atac_overlap <- bed_overlaps(bed.1 = da_peaks_sten,
                                    bed.2 = sc_atac_peaks,
                                    names = c("peak", "atac_peak"))

## upset plot of scATAC and stencilling markers:
  peaks_list <- list(da_peaks_sten,sc_atac_peaks)
 
  atac_sten_overlap <- find_overlaps(tables = peaks_list, 
                                     names = c("Stencilling Markers", "scATAC Markers"))
  
  atac_sten_overlap$upset_plot
  grid.text("Stencilling markers Vs pseudo-bulked scATAC markers", x = 0.68, y=0.97, gp = gpar(fontsize=15))
  
  
## finding the width of the peaks per catgeory
  # getting the matrix:
  #overlap_mat <- atac_sten_overlap$binary_matrix
  #overlap_mat <- overlap_mat %>%
  #  mutate(category= case_when(`Stencilling Markers`==1 & `scATAC Markers`==1 ~ "Both",
   #                            `Stencilling Markers`==1 & `scATAC Markers`==0 ~ "Stencilling_Unique",
    #                           `Stencilling Markers`==0 & `scATAC Markers`==1 ~ "scATAC_Unique",
     #                          `Stencilling Markers`==0 & `scATAC Markers`==0 ~ "none")) %>%
  #  separate(peak, into = c("chr", "start", "end"), sep = "_", remove = FALSE) %>%
  #  mutate(width= as.numeric(end)- as.numeric(start))
  
  #ggplot(overlap_mat, aes(x=category, y = width, fill = category))+
    #geom_boxplot(stat = "boxplot")+
    #theme_few()+
    #scale_fill_tableau(direction = -1)+
    #scale_y_log10()
  
  
## since the matrix was formed by merging the peaks--this approach of finding the width of the peaks might be wrong:
  overlap_mat <- atac_sten_overlap$binary_matrix
  merged_peaks_bed <- overlap_mat %>%
    select(peak) %>%
    separate(peak, into = c("chr", "start", "end"))
  
  # stencilling:
  merged_sten <- bed_overlaps(bed.1 = merged_peaks_bed,
                              bed.2 = da_peaks_sten,
                              names = c("peak", "stencilling"))
  
  merged_sten$peak <- gsub("-", "_", merged_sten$peak)
  
  overlap_mat <- left_join(overlap_mat, merged_sten, by = "peak") %>%
    group_by(peak) %>%
    mutate(stencilling= paste(stencilling, collapse = ",")) %>%
    ungroup() %>%
    distinct()
  
  #scATAC
  merged_atac <- bed_overlaps(bed.1 = merged_peaks_bed,
                              bed.2 = sc_atac_peaks,
                              names = c("peak", "scATAC"))
  
  merged_atac$peak <- gsub("-", "_", merged_atac$peak)
  
  overlap_mat <- left_join(overlap_mat, merged_atac, by = "peak") %>%
    group_by(peak) %>%
    mutate(scATAC= paste(scATAC, collapse = ",")) %>%
    ungroup() %>%
    distinct()
  
## finding categories of peaks:
  overlap_mat <- overlap_mat %>%
    mutate(category= case_when(`Stencilling Markers`==1 & `scATAC Markers`==1 ~ "Both",
                               `Stencilling Markers`==1 & `scATAC Markers`==0 ~ "Stencilling_Unique",
                               `Stencilling Markers`==0 & `scATAC Markers`==1 ~ "scATAC_Unique",
                               TRUE ~ "none")) %>%
    separate_rows(c("stencilling", "scATAC"), sep = ",") %>%
    separate(stencilling, into = c("sten_chr", "sten_start", "sten_end"), sep = "-", remove = FALSE) %>%
    mutate(sten_width = as.numeric(sten_end)- as.numeric(sten_start)) %>%
    select(-c("sten_chr", "sten_start", "sten_end")) %>%
    separate(scATAC, into = c("atac_chr", "atac_start", "atac_end"), sep = "-", remove = FALSE) %>%
    mutate(atac_width = as.numeric(atac_end)- as.numeric(atac_start)) %>%
    select(-c("atac_chr", "atac_start", "atac_end"))
    
  
## making long format df for grouped boxplot:
  boxplot_df <- overlap_mat %>%
    pivot_longer(cols = c("atac_width", "sten_width"), names_to = "width_category", values_to = "width")
  
  boxplot_df$category <- factor(boxplot_df$category, levels = c("Stencilling_Unique","scATAC_Unique","Both" ))
  
  ggplot(boxplot_df, aes(x=category, y = width, fill = width_category))+
    geom_boxplot(stat = "boxplot", width=0.5, outlier.size = 1, outlier.alpha = 0.2)+
    theme_few(base_size = 14)+
    scale_fill_tableau(name= "Category", labels= c("atac_width"= "scATAC", "sten_width"="Stencilling"))+
    scale_y_log10()+
    labs(title = "Peak width of Stencilling and scATAC Markers",
         x= "", 
         y= "Width (bp)")+
    scale_x_discrete(labels= c(1,2,3))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", size=18, margin = margin(r = 12)),
          legend.position = "top",
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank())
  
  
## finding overlaps with ALL mer11:
  atac_sten_mer11 <-  find_overlaps_reference(reference = all_mer11,
                                              tables = peaks_list,
                                              names = c("Stencilling", "scATAC"))
  
  upset(atac_sten_mer11, order.by = "freq")
  
  
  
##----------------------------------------------------------------------------------------------------------
# Qudrant scatter plot of accessibility and methylation:
##----------------------------------------------------------------------------------------------------------
  
## plotting:
  
  thresh <- 0.5
  padj_thresh <- 1e-4
  
  
  #quad_plot <- all_sig_df %>%
    #select(bh_5mC, mod6mA_Adjusted, mod6mA_or, mod5mC_or, log2_mod5mC_or, log2_mod6mA_or) %>%
    #filter(mod6mA_Adjusted < padj_thresh) %>%
    #filter(bh_5mC < padj_thresh) %>%
    #filter(is.finite(log2_mod5mC_or), 
      #     is.finite(log2_mod6mA_or)) %>%
    #filter(log2_mod6mA_or > thresh | log2_mod6mA_or < -thresh) %>%
    #filter(log2_mod5mC_or > thresh | log2_mod5mC_or < -thresh)

## Re-doing Fishers test for 5mC data to remove all Infinite odds ratio values:
    col_names_5mc <- colnames(all_sig_df)[grep("mod5mC_contingency", colnames(all_sig_df))]
    
    contin_5mc <- all_sig_df %>%
      select(col_names_5mc, Index)
    
    # replacing all 0 with 1
    contin_5mc[contin_5mc==0] <- 1
    
    # splitting by rows
    contin_5mc_split <- split(contin_5mc, rownames(contin_5mc))
  
    # fishers test
    fishers_results <- lapply(contin_5mc_split, function(x){
    
      index_val <- x$Index
      df <- x %>%
        select(-Index) %>%
        setNames(c("a", "b", "c", "d"))
      
      df_mat <- matrix(as.numeric(df), nrow = 2, ncol = 2, byrow = TRUE)
      df_fish <- fisher.test(df_mat)
      
       data.frame(a = df$a,
                  b= df$b,
                  c = df$c,
                  d = df$d,
                  p_val = df_fish$p.value,
                  or_5mC_corrected = df_fish$estimate,
                  Index= index_val)
      
    })
    
    fishers_5mc <- do.call(rbind,fishers_results)
    
## merging with the main data frame:
    all_sig_df <- left_join(all_sig_df,fishers_5mc[, c("Index", "or_5mC_corrected")], by = "Index" )
    all_sig_df$log2_mod5mC_or_CORRECTED <- -log2(all_sig_df$or_5mC_corrected)
    
## plotting all categories of peaks: MER11 and non-MER11
  quad_plot <- all_sig_df %>% 
    select(bh_5mC, mod6mA_Adjusted, mod6mA_or, mod5mC_or, log2_mod5mC_or_CORRECTED, log2_mod6mA_or, is_mer11) %>% 
    mutate(sig_6ma= mod6mA_Adjusted < padj_thresh & abs(log2_mod6mA_or) > thresh,
           sig_5mc= bh_5mC < padj_thresh & abs(log2_mod5mC_or_CORRECTED) > thresh,
           is_mer11 = gsub("FALSE", "Non_MER11", is_mer11),
           is_mer11 = gsub("TRUE", "MER11", is_mer11)) %>%
    mutate(category = case_when(sig_6ma==TRUE & sig_5mc==TRUE ~ "Both",
                                sig_6ma==TRUE & sig_5mc==FALSE ~ "6mA",
                                sig_6ma==FALSE & sig_5mc==TRUE ~ "5mC",
                                TRUE ~ "Non_sig")) %>%
    filter(category != "Non_sig") %>%
  mutate(quad_category = case_when(
    log2_mod6mA_or < -0.5 & log2_mod5mC_or_CORRECTED > 0.5 ~ "Q1",
    log2_mod6mA_or > 0.5 & log2_mod5mC_or_CORRECTED > 0.5 ~ "Q2",
    log2_mod6mA_or < -0.5 & log2_mod5mC_or_CORRECTED < -0.5 ~ "Q3",
    log2_mod6mA_or > 0.5 & log2_mod5mC_or_CORRECTED < -0.5 ~ "Q4",
    TRUE ~ "None"
  )) %>%
    arrange((category))
  

  
  ##  accessibility vs methylation
  # Mer11
  ggplot(quad_plot[quad_plot$is_mer11=="MER11",], aes(x= log2_mod6mA_or, y= log2_mod5mC_or_CORRECTED, colour = category))+
    scale_color_manual(name= "", values = c("5mC"= "lightblue", "6mA"= "#0072a1", "Both" = "#d55e00"))+
    geom_point(size=1.8, alpha=0.8)+
    theme_few(base_size = 18)+
    xlim(-6,6)+
    ylim(-9,9)+
    #facet_wrap(~is_mer11) +
    geom_hline(yintercept = 0, color="grey")+
    geom_vline(xintercept = 0, color= "grey")+
    labs(title = "Correlation scatter plot between Accessibility and Methylation",
         x= bquote("Accessibility (" ~Log[2] ~ "Odds Ratio)"),
         y= bquote("5mC Methylation (" ~ Log[2] ~ "Odds Ratio)"))+
    theme(axis.text = element_text(face = "bold", size = 18),
          plot.title = element_text(face = "bold", size = 16),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 16)) 

  #Non mer11
  ggplot(quad_plot[quad_plot$is_mer11=="Non_MER11",], aes(x= log2_mod6mA_or, y= log2_mod5mC_or_CORRECTED, colour = category))+
    scale_color_manual(name= "", values = c("5mC"= "lightblue", "6mA"= "#0072a1", "Both" = "#d55e00"))+
    geom_point(size=1.3, alpha=0.8)+
    theme_few(base_size = 18)+
    xlim(-6,6)+
    ylim(-9,9)+
    #facet_wrap(~is_mer11) +
    geom_hline(yintercept = 0, color="grey")+
    geom_vline(xintercept = 0, color= "grey")+
    labs(title = "Correlation scatter plot between Accessibility and Methylation",
         x= bquote("Accessibility (" ~Log[2] ~ "Odds Ratio)"),
         y= bquote("5mC Methylation (" ~ Log[2] ~ "Odds Ratio)"))+
    theme(axis.text = element_text(face = "bold", size = 18),
          plot.title = element_text(face = "bold", size = 16),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 16))

##----------------------------------------------------------------------------------------------------------
# Heatmap:
##----------------------------------------------------------------------------------------------------------
## Finding the proportion per modification:
  
  mods_df_2 <-  all_sig_df
  mods_df_2 <-  add_mod_rate(mods_df_2)
  
  sc_fire_mat <- mods_df_2 %>%
    filter(mod6mA_sig_gain == TRUE & is_mer11==TRUE) %>%
    select(contains("rate"), peak) %>%
    column_to_rownames("peak") %>%
    distinct()
  
  colnames(sc_fire_mat) <- gsub("_rateA", "_SC_WT", colnames(sc_fire_mat))
  colnames(sc_fire_mat) <- gsub("_rateB", "_SC_KO", colnames(sc_fire_mat))
  
  num_peaks <- nrow(sc_fire_mat)
  sc_fire_heatmap <- pheatmap(sc_fire_mat, 
                              show_rownames = FALSE,
                              clustering_distance_rows = "euclidean",
                              clustering_method = "ward.D2",
                              cluster_rows = T, 
                              cluster_cols = F,
                              color = viridis(100),
                              border_color = NA,
                              cutree_rows= 5,
                              labels_row = TRUE,
                              gaps_col = seq(0,8,2),
                              fontsize_col = 8,
                              angle_col = "315",
                              main = paste("SC stage: FIRE", num_peaks ,"peaks"),
                              breaks = seq(0, 1, length.out = 101))
  
  sc_clus <- data.frame("Clusters"= factor(cutree(sc_fire_heatmap$tree_row, k = 5)))
  
  sc_fire_heatmap <- pheatmap(sc_fire_mat, 
                              show_rownames = FALSE,
                              clustering_distance_rows = "euclidean",
                              clustering_method = "ward.D2",
                              cluster_rows = T, 
                              cluster_cols = F,
                              color = viridis(100),
                              border_color = NA,
                              cutree_rows= 5,
                              labels_row = TRUE,
                              annotation_row = sc_clus,
                              gaps_col = seq(0,8,2),
                              fontsize_col = 8,
                              angle_col = "315",
                              main = paste("SC stage: 6mA MER11-", num_peaks ,"peaks"),
                              breaks = seq(0, 1, length.out = 101))
  
## Saving the peaks as a bed file to extract reads from Chip seq data:
  mer11_6ma_bed <- sc_fire_mat %>%
    rownames_to_column("peak") %>%
    select(peak) %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-")
  
  write.table(mer11_6ma_bed,"stencilling_320_MER11.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)

##---------------------------------------------------------------------------------------------------------------------
## Adding Chip seq signals and ATAC signals to the heatmap:
##---------------------------------------------------------------------------------------------------------------------
## Loading the scATAC read counts
  atac_reads <- read.delim("sten_mer11_all_ATAC_reads.bed", header = FALSE, sep = "\t")
  
  atac_reads <- atac_reads %>%
    setNames(c("chr", "start","end", "counts", "sample"))
  
## Loading the ChIP reads
  k27_reads <- read.delim("sorted_sten_mer11_H3K27ac.bed" , 
                          header = FALSE,
                          sep = "\t") %>%
    setNames(c("chr", "start","end", "length","counts", "sample")) %>%
    select(-length)
  
  
  k9_reads <- read.delim("sorted_sten_mer11_H3K9me3.bed" ,
                         header = FALSE, 
                         sep = "\t") %>%
    setNames(c("chr", "start","end", "length","counts", "sample"))%>%
    select(-length)
  
  
## laoding the library sizes of each assay:
  # K9me3
  k9_lib_size <- read.delim("../../R_Projects/Mer11_cluster_methylation/library_sizes_H3K9me3", sep = "\t", header = FALSE)
  k9_lib_size <- k9_lib_size %>%
    mutate(V1= gsub("\\s*\\+.*\\)", "", V1)) %>%
    separate(V1, into = c("size", "sample"), sep = "\\s+")
  
  # K27ac
  k27_lib_size <- read.delim("../../R_Projects/Mer11_cluster_methylation/library_sizes_H3K27ac", sep = "\t", header = FALSE)
  k27_lib_size <- k27_lib_size %>%
    mutate(V1= gsub("\\s*\\+.*\\)", "", V1)) %>%
    separate(V1, into = c("size", "sample"), sep = "\\s+")
  
  # Atac
  atac_lib_size <- read.delim("../../R_Projects/Mer11_cluster_methylation/atac_fragments_size.bed", sep = "\t", header = FALSE) %>%
    setNames(c("size", "sample"))
 
## making a raw matrix of all the data:
  # Filtering to onlu contain SC stage values:
  
  df_list <- list(atac_reads, k9_reads, k27_reads)
  df_list_raw <- lapply(df_list, function(x){
    
    x <- x %>%
      mutate(Mer11= paste(chr, start, end, sep = "-")) %>%
      select(4:6) %>%
      mutate(sample= gsub("stage_0_", "SC_", sample)) %>%
      pivot_wider(names_from = sample, values_from = counts) %>%
      column_to_rownames("Mer11")%>%
      select(starts_with("SC")) %>%
      as.matrix()
    
    return(x)
    
  })
  
  names(df_list_raw) <- c("atac_raw_df","k9_raw_df", "k27_raw_df")
  list2env(df_list_raw, envir =.GlobalEnv)
  
 
## normalising each assay separately:
  ## ATAC and ChIP reads:
  # CPM method: 
  # Raw count/ library size (num of properly aligend pairs in the BAM file), * by scaling factor = 1e6
  
  # K9 reads
    lib_size <- as.numeric(k9_lib_size$size)
    names(lib_size) <- k9_lib_size$sample
    lib_size <- lib_size[colnames(k9_raw_df)] # ensuring that the colnames are in the same order as the lib size vector
    
    k9_norm_df <- sweep(k9_raw_df,2,lib_size, FUN = "/") * 1e+6
    
    #max_score <- rowMaxs(k9_norm_df)
    #k9_norm_df <- sweep(k9_norm_df, 1, max_score, FUN = "/")
    
    colnames(k9_norm_df) <- paste0("K9_", colnames(k9_norm_df))
    
  # K27 reads
    lib_size <- as.numeric(k27_lib_size$size)
    names(lib_size) <- k27_lib_size$sample
    lib_size <- lib_size[colnames(k27_raw_df)]
    
    k27_norm_df <- sweep(k27_raw_df,2,lib_size, FUN = "/") *1e+6
    
    #max_score <- rowMaxs(k27_norm_df)
    #k27_norm_df <- sweep(k27_norm_df, 1, max_score, FUN = "/")
    
    colnames(k27_norm_df) <- paste0("K27_", colnames(k27_norm_df))
  
  # atac reads:
    lib_size <- as.numeric(atac_lib_size$size)
    names(lib_size) <- atac_lib_size$sample
    lib_size <- lib_size[colnames(atac_raw_df)]
    
    atac_norm <- sweep(atac_raw_df,2,lib_size, FUN = "/") *1e+6
  
  
## Collapsing the replicates for the ATAC data: finding mean normalised counts
    samples <- sub("_[1:2]$", "", colnames(atac_norm))
    
    atac_norm_df <- sapply(unique(samples), function(x){
      rowMeans(atac_norm[,samples == x, drop= TRUE])
      
    })
    
    #max_score <- rowMaxs(atac_norm_df)
    #atac_norm_df <- sweep(atac_norm_df, 1, max_score, FUN = "/")
    
    colnames(atac_norm_df) <- paste0("ATAC_", colnames(atac_norm_df))

  
## combining these normalised matrices togetehr:
  chip_atac_mat <- as.matrix(cbind(atac_norm_df, k9_norm_df, k27_norm_df, sc_fire_mat))
  
## Scaling the matrix:
  ## Scaling row MAX:
  max_score <- rowMaxs(chip_atac_mat)
  max_score[max_score==0] <- 1 # to remove NaN values after divsion
  
  chip_atac_mat <- sweep(chip_atac_mat, 1, max_score, FUN = "/")
  chip_atac_mat <- chip_atac_mat %>%
    as.data.frame() %>%
    relocate(ends_with("KO"), .after = ends_with("WT")) %>%
    as.matrix()
  
  # Z score scaling
  chip_atac_mat_scaled <- t(scale(t(chip_atac_mat)))
  chip_atac_mat_scaled <- chip_atac_mat_scaled %>%
    as.data.frame() %>%
    relocate(ends_with("KO"), .after = ends_with("WT")) %>%
    as.matrix()
  
  
## Heatmap:
  num_peaks <- nrow(chip_atac_mat_scaled)
  sc_chip_heatmap <- pheatmap(chip_atac_mat_scaled, 
                              show_rownames = FALSE,
                              clustering_distance_rows = "euclidean",
                              clustering_method = "ward.D2",
                              cluster_rows = T, 
                              cluster_cols = F,
                              color = viridis(100),
                              border_color = NA,
                              cutree_rows= 8,
                              labels_row = TRUE,
                              gaps_col = seq(0,14,7),
                              fontsize_col = 8,
                              angle_col = "315",
                              main = paste("SC stage: ChIP", num_peaks ,"peaks"),
                              breaks = seq(-1, 1, length.out = 101))
  
  sc_clus <- data.frame("Clusters"= factor(cutree(sc_chip_heatmap$tree_row, k = 5)))
  
  sc_chip_heatmap <- pheatmap(chip_atac_mat_scaled, 
                              show_rownames = FALSE,
                              clustering_distance_rows = "euclidean",
                              clustering_method = "ward.D2",
                              cluster_rows = T, 
                              cluster_cols = F,
                              color = viridis(100),
                              border_color = NA,
                              cutree_rows= 5,
                              labels_row = TRUE,
                              annotation_row = sc_clus,
                              gaps_col = seq(0,14,7),
                              fontsize_col = 8,
                              angle_col = "315",
                              main = paste("SC stage: 6mA MER11-", num_peaks ,"peaks"),
                              breaks = seq(-1, 1, length.out = 101))

  
##--------------------------------------------------------------------------------------------- 
# Haplotype specific diff peaks from KO WGS   
##--------------------------------------------------------------------------------------------- 
## loading the file:
  hap_wgs <- read.delim("znf808_KO_AS_hap_CORRECTED.tsv", sep = "\t")
  
## BH adjustemnet of p values--- 5mC and 6mA
  bh_6ma <- p.adjust(hap_wgs$mod6mA_fisher_bpv, method = "BH")
  hap_wgs$bh_6ma <- bh_6ma
  
  bh_5mc <- p.adjust(hap_wgs$mod5mC_fisher_bpv, method = "BH")
  hap_wgs$bh_5mc <- bh_5mc
  
  hap_wgs$log2_6ma <- log2(hap_wgs$mod6mA_or)
  hap_wgs$log2_5mc <- log2(hap_wgs$mod5mC_or)
  
## overlap with mer11:
  peaks_mer11 <- bed_overlaps(bed.1 = hap_wgs[,1:3],
                              bed.2 = all_mer11,
                              names = c("peaks", "mer11"))
  
  hap_wgs <- hap_wgs %>%
    mutate(peaks= paste(chrom, start, stop, sep = "-"),
           is_mer11= peaks %in% peaks_mer11$peaks) 
  hap_wgs <- left_join(hap_wgs, peaks_mer11, by = "peaks")
  hap_wgs <-  left_join(hap_wgs, mer11_overlap, by = c("mer11"="mer11_220"))
  
## filtering only the sig peaks:
  hap_thresh <- 1e-4
  hap_thresh_or <- 0.5
  
  sig_hap_6mA <- hap_wgs %>%
    filter(bh_6ma < hap_thresh & abs(log2_6ma) > hap_thresh_or) %>%
    filter(is_mer11==TRUE)
  
  sig_hap_5mc <- hap_wgs %>%
    filter(bh_5mc < hap_thresh & abs(log2_5mc) > hap_thresh_or) %>%
    filter(is_mer11==TRUE)
  
  sig_hap_6ma_5mc <- hap_wgs %>%
    filter(bh_5mc < hap_thresh & abs(log2_5mc) > hap_thresh_or & bh_6ma < hap_thresh & abs(log2_6ma) > hap_thresh_or) %>%
    filter(is_mer11==TRUE)
  
  