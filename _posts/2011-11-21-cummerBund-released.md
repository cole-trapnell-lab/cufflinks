---
layout: post
title: CummeRbund released
description: "CummeRbund released"
modified: 2011-11-21
category: announcements
tags: [CummeRbund]
comments: false
share: true
---
Extracting biological insight from transcript-level RNA-Seq analysis can be very challenging. Due to the volume and complexity of output from RNA-Seq analysis pipelines, many users may choose to focus only on gene-level results, and thus miss crucial biological insights that a transcript-level analysis can provide. We are happy to present CummeRbund, an R/Bioconductor package that simplifies the organization, access, exploration, and visualization of the various output files of a Cuffdiff differential expression analysis. CummeRbund begins by re-organizing the Cuffdiff output files, and storing these data in a local SQLite database. During this process, CummeRbund indexes the data to speed up access to specific feature data (genes, isoforms, TSS, CDS, etc.), and preserves the various relationships between these features. Access to data elements is done through R via the RSQLite package and data are presented in appropriately structured R classes with various convenience functions designed to streamline your workflow.

CummeRbund simplifies the way in which you access and analyze your RNA-Seq data. Features include:

- Numerous plotting methods to allow for rapid visualization of RNA-Seq data quality and global statistics such as FPKM distribution, as well as simple routines for plotting expression levels for one or more genes, their isoforms, TSS groups, or CDS groups.
	- Plots are generally publication-ready. However, because they are built with [ggplot2](http://had.co.nz/ggplot2/), plot objects returned by the plotting routines are easy to manipulate in a live R session so you can tweak them to your specification.
- Persistent storage and indexing of cuffdiff data in a relational database.
	- Quicker searching for feature-specific information
	- Seamless aggregation of all related data points.
	- Direct SQL/RSQLite querying of cuffdiff data for access to complex datasets using specialized queries.
- Formalized R classes for data access and manipulation.
	- 'Pointer' classes call data directly from SQL tables.
	- 'Data' classes contain queried results and can be manipulated directly in R.
	- Methods for direct access to FPKM values, differential expression data, and additional annotation of all gene-, isoform-, TSS-, and CDS-level features.
	- Output formats allow for standard browsing/analysis in R (data.frame, list, etc).
- Geneset-level data access:
	- Relevant data for meaningful subsets of genes (e.g. 'significantly regulated') can be quickly aggregated and stored in a geneset object.
	- plot wrappers for heatmaps, expression profiles, barplots, volcano plots, scatterplots, etc.
- Gene-centric data access:
	- All interconnected data for a single gene are aggregated into a single data object.
	- Access to all relevant FPKM, differential expression, and annotation data for all features of a single gene.
	- Gene-specific plots (expression profiles)
- Individual Feature-level data access as well:
	- Individual Isoforms, TSS, CDS, etc. are also available.
	- Direct access to distribution test results for alternative splicing, TSS groups, and CDS usage.
- CummeRbund is being made freely available under the OSI approved Artistic License 2.0 and can be downloaded from [http://compbio.mit.edu/cummeRbund/](http://compbio.mit.edu/cummeRbund/). CummeRbund has also been included as part of the new R/[Bioconductor](http://www.bioconductor.org/) v2.9 release and can be installed in a similar manner to standard Bioconductor packages.


