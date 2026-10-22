# Targeted-Stencilling

This repo contains the scripts used for analysing the results from targeted chromatin stencilling done on the ZNF808 KO cell line.

## Background

ZNF808 is a primate specific transcription factor that functions by binding to and repressing a specific class of transposable elements called the MER11 elements. Loss of ZNF808 causes pancreatic agenesis, a rare developmental disorder in which the pancreas fails to develop. When ZNF808 is lost, the MER11 elements which have the potential to act as an enhancer, are no longer repressed, and therefore gains accessibility, causing aberrant gene activation. In vitro stem cell studies have shown higher expression of liver like gene expression in the ZNF808 KO cell lines during pancreatic differentiation.

## Stencilling/Fiber-seq

[Fiber-seq](https://www.science.org/doi/10.1126/science.aaz1646) also known as stencilling, is a long-read sequencing based assay that simultaneously profiles chromatin accessibility, DNA methylation, nucleosome positioning and transcription factor footprints at a single molecule resolution.

## Targeted Stencilling: Stencilling + ONT's Adaptive Sampling

To understand how the chromatin accessibility and DNA methylation at the MER11 elements change in the absence of ZNF808 we performed targeted stencilling. This was done by using [Oxford Nanopore's (ONT) Adaptive Sampling (AS)](https://nanoporetech.com/document/adaptive-sampling) technology that selectively enriches for reads mapping to the MER11 elements.

Through this targeted approach we get a \>100x read depth (combining WGS and AS) at the MER11 elements that allows us to understand the epigenetic states at a single molecule resolution. The increased coverage also allows the reads to be phased enabling us to resolve haplotype-specific patterns of regulation.

## Analysis:

All the pre-processing of the files were done using [FIBERTOOLS](https://fiberseq.github.io/fibertools/fibertools.html) and the [EPI2ME pipeline](https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-human-variation/) by ONT. For further processing, including finding the differentially accessible and methylated regions [SMGStats](https://github.com/exetergrdd/SMGStats.jl/tree/main), an in-house developed package was used. For plotting the reads along with the linker, nucleosome and FIRE labelling, [SMGDataPlotting](https://github.com/exetergrdd/SMGDataPlotting.jl) package was used.

-   `data`: This folder has all the data required to perform the downstream analysis.

-   `functions`: This folder has all the functions used for running the script

-   `scripts`: This folder contains a markdown script along with the outputs and the plots
