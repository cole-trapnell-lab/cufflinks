---
layout: page
permalink: /papers/
description: "Papers describing the Cufflinks suite."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

Cufflinks is an ongoing research project as well as a suite of tools.  Here are the papers that describe the science behind the programs.  If you use Cufflinks, **please cite these papers** in your work!

# [Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation](http://dx.doi.org/10.1038/nbt.1621) 

Cole Trapnell, Brian Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Jeltje van Baren, Steven Salzberg, Barbara Wold, Lior Pachter. 

***Nature Biotechnology***, 2010 

High-throughput mRNA sequencing (RNA-Seq) promises simultaneous transcript discovery and abundance estimation. However, this would require algorithms that are not restricted by prior gene annotations and that account for alternative transcription and splicing. Here we introduce such algorithms in an open-source software program called Cufflinks. To test Cufflinks, we sequenced and analyzed >430 million paired 75-bp RNA-Seq reads from a mouse myoblast cell line over a differentiation time series. We detected 13,692 known transcripts and 3,724 previously unannotated ones, 62% of which are supported by independent expression data or by homologous genes in other species. Over the time series, 330 genes showed complete switches in the dominant transcription start site (TSS) or splice isoform, and we observed more subtle shifts in 1,304 other genes. These results suggest that Cufflinks can illuminate the substantial regulatory flexibility and complexity in even this well-studied model of muscle development and that it can improve transcriptome-based genome annotation.

doi:10.1038/nbt.1621

***Note:*** *This is the original Cufflinks paper.  Please cite this paper if you use Cufflinks in your work.*


# [Improving RNA-Seq expression estimates by correcting for fragment bias](http://genomebiology.com/2011/12/3/R22/abstract)

Adam Roberts, Cole Trapnell, Julie Donaghey, John L. Rinn, Lior Pachter.  

***Genome Biology***, 2011 

The biochemistry of RNA-Seq library preparation results in cDNA fragments that are not uniformly distributed within the transcripts they represent. This non-uniformity must be accounted for when estimating expression levels, and we show how to perform the needed corrections using a likelihood based approach. We find improvements in expression estimates as measured by correlation with independently performed qRT-PCR and show that correction of bias leads to improved replicability of results across libraries and sequencing technologies.

doi:10.1186/gb-2011-12-3-r22

***Note:*** *This paper describes improvements made to Cufflinks to handle bias in RNA-Seq read coverage.  Please cite this paper if you use Cufflinks with the -b option in your work.*


# [Identification of novel transcripts in annotated genomes using RNA-Seq](http://bioinformatics.oxfordjournals.org/content/27/17/2325)

Adam Roberts, Harold Pimentel, Cole Trapnell, Lior Pachter.  

***Bioinformatics***, 2011 

*Summary:* We describe a new ‘reference annotation based transcript assembly’ problem for RNA-Seq data that involves assembling novel transcripts in the context of an existing annotation. This problem arises in the analysis of expression in model organisms, where it is desirable to leverage existing annotations for discovering novel transcripts. We present an algorithm for reference annotation-based transcript assembly and show how it can be used to rapidly investigate novel transcripts revealed by RNA-Seq in comparison with a reference annotation.

*Availability:* The methods described in this article are implemented in the Cufflinks suite of software for RNA-Seq, freely available from http://bio.math.berkeley.edu/cufflinks. The software is released under the BOOST license.

doi:10.1093/bioinformatics/btr355

***Note:*** *This paper describes the RABT assembly algorithm.  Please cite this paper if you use Cufflinks in RABT mode in your work.*


# [Differential analysis of gene regulation at transcript resolution with RNA-seq](http://dx.doi.org/10.1038/nbt.2450)

Cole Trapnell, David Hendrickson, Martin Sauvageau, Loyal Goff, John L. Rinn, Lior Pachter  

***Nature Biotechnology***, 2012 

Differential analysis of gene and transcript expression using high-throughput RNA sequencing (RNA-seq) is complicated by several sources of measurement variability and poses numerous statistical challenges. We present Cuffdiff 2, an algorithm that estimates expression at transcript-level resolution and controls for variability evident across replicate libraries. Cuffdiff 2 robustly identifies differentially expressed transcripts and genes and reveals differential splicing and promoter-preference changes. We demonstrate the accuracy of our approach through differential analysis of lung fibroblasts in response to loss of the developmental transcription factor HOXA1, which we show is required for lung fibroblast and HeLa cell cycle progression. Loss of HOXA1 results in significant expression level changes in thousands of individual transcripts, along with isoform switching events in key regulators of the cell cycle. Cuffdiff 2 performs robust differential analysis in RNA-seq experiments at transcript resolution, revealing a layer of regulation not readily observable with other high-throughput technologies.

doi:10.1038/nbt.2450

***Note:*** *This paper describes Cuffdiff 2.  Please cite this paper if you use Cuffdiff in your work.*

