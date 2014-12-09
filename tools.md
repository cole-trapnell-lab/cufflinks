---
layout: page
permalink: /tools/
description: "Other tools for analysis high-throughput experiments."
modified: 2013-09-11
tags: [Bowtie, Tophat, Monocle, CummeRbund]
---

# Bowtie: ultrafast short read alignment

[Bowtie](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

Bowtie is provided under the OSI-approved Artistic License 2.0.

# TopHat: alignment of short RNA-Seq reads

[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons. 

TopHat is provided under the OSI-approved Artistic License 2.0.

# CummeRbund: visualization of RNA-Seq differential analysis

[CummeRbund](http://compbio.mit.edu/cummeRbund/) is an R package that is designed to aid and simplify the task of analyzing Cufflinks RNA-Seq output.

CummeRbund is provided under the OSI-approved Artistic License 2.0.

# Monocle: Differential expression for single-cell RNA-Seq and qPCR.

[Monocle](http://monocle-bio.sourceforge.net/) is a toolkit for analyzing single-cell gene expression experiments. Monocle was designed for RNA-Seq, but can also work with single cell qPCR. It performs differential expression analysis, and can find genes that differ between cell types or between cell states. When used to study an ongoing biological process such as cell differentiation, Monocle learns that process and places cells in order according to their progress through it. Monocle finds genes that are dynamically regulated during that process. 

Monocle is provided under the OSI-approved Artistic License (version 2.0)
