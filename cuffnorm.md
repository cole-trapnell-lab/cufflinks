---
layout: page
permalink: /cuffnorm/
description: "Documentation for Cuffnorm."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

# Running Cuffnorm

Cufflinks includes a program, "Cuffnorm", that you can use to generate tables of expression values that are properly normalized for library size. From the command line, run cuffnorm as follows:

**<tt>cuffnorm [options] \<transcripts.gtf\> \<sample1_replicate1.sam[,...,sample1_replicateM.sam]\> \\</tt>**

**<tt>\<sample2_replicate1.sam[,...,sample2_replicateM.sam]\>... [sampleN.sam_replicate1.sam[,...,sample2_replicateM.sam]]</tt>**

## Cuffnorm Input

Running Cuffnorm is very similar to running Cuffdiff. Cuffnorm takes a GTF2/GFF3 file of transcripts as input, along with two or more SAM, BAM, or CXB files for two or more samples. It produces a number of output files that contain expression levels and normalized fragment counts at the level of transcripts, primary transcripts, and genes. It also tracks changes in the relative abundance of transcripts sharing a common transcription start site, and in the relative abundances of the primary transcripts of each gene. Tracking the former allows one to see changes in splicing, and the latter lets one see changes in relative promoter use within a gene.

If you have more than one replicate for a sample, supply the SAM files for the sample as a single comma-separated list. It is not necessary to have the same number of replicates for each sample.

Note that Cuffnorm can also accepted BAM files (which are binary, compressed SAM files). It can also accept CXB files produced by Cuffquant. Note that mixing SAM and BAM files is supported, but you cannot currently mix CXB and SAM/BAM. If one of the samples is supplied as a CXB file, all of the samples must be supplied as CXB files.

Cuffnorm also requires a GFF/GTF file, conforming to the same specifications as needed for Cuffdiff.

## Cuffnorm arguments	

**<tt><transcripts.(gtf/gff)></tt>**	

A transcript annotation file produced by cufflinks, cuffcompare, or other source.

**<tt><sample1.(sam/bam/cxb)></tt>**	

A SAM file of aligned RNA-Seq reads. If more than two are provided, Cuffdiff tests for differential expression and regulation between all pairs of samples.

## Cuffnorm options	

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-o/--output-dir \<string\></tt>**

Sets the name of the directory in which Cuffdiff will write all of its output. The default is "./".

**<tt>-L/--labels \<label1,label2,...,labelN\></tt>**

Specify a label for each sample, which will be included in various output files produced by Cuffdiff.

**<tt>-p/--num-threads \<int\></tt>**

Use this many threads to align reads. The default is 1.

**<tt>--total-hits-norm</tt>**

With this option, Cuffquant counts all fragments, including those not compatible with any reference transcript, towards the number of mapped fragments used in the FPKM denominator. It is inactive by default.

**<tt>--compatible-hits-norm</tt>**

With this option, Cuffnorm counts only those fragments compatible with some reference transcript towards the number of mapped fragments used in the FPKM denominator. It is active by default.


**<tt>--library-type</tt>**

See [Library Types](#library-types)

**<tt>--library-norm-method</tt>**	

See [Library Normalization Methods](#library-normalization-methods)

**<tt>--output-format</tt>**

See [Output format options](#output-format-options)

## Cuffnorm advanced options	

**<tt>-v/--verbose</tt>**	

Print lots of status updates and other diagnostic information.

**<tt>-q/--quiet</tt>**

Suppress messages other than serious warnings and errors.

**<tt>--no-update-check</tt>**

Turns off the automatic routine that contacts the Cufflinks server to check for a more recent version.

## Cuffnorm output

Cuffnorm outputs a set of files containing normalized expression levels for each gene, transcript, TSS group, and CDS group in the experiment. It does not perform differential expression analysis. To assess the significance of changes in expression for genes and transcripts between conditions, use Cuffdiff. Cuffnorm's output files are useful when you have many samples and you simply want to cluster them or plot expression levels of genes important in your study.

By default, Cuffnorm reports expression levels in the "simple-table" tab-delimted text files. The program also reports information about your samples and about the genes, transcripts, TSS groups, and CDS groups as tab delimited text files. Note that these files have a different format than the files used by Cuffdiff. However, you can direct Cuffnorm to report its output in the same format as used by Cuffdiff if you wish. Simply supply the option <tt>--output-format cuffdiff</tt> at the command line.

Cuffnorm will report both FPKM values and **normalized**, estimates for the number of fragments that originate from each gene, transcript, TSS group, and CDS group. Note that because these counts are already normalized to account for differences in library size, they should not be used with downstream differential expression tools that require **raw** counts as input.

To see the details of the simple table format used by Cuffnorm, refer to the simple table expression format, simple table sample attribute format, and simple table feature (e.g. gene) attribute format sections below. 

# Library Types

In cases where Cufflinks cannot determine the platform and protocol used to generate input reads, you can supply this information manually, which will allow Cufflinks to infer source strand information with certain protocols. The available options are listed below. For paired-end data, we currently only support protocols where reads are point towards each other.

Library Type|	Examples|	Description
---|----|----
fr-unstranded (default)|	Standard Illumina|	Reads from the left-most end of the fragment (in transcript coordinates) map to the transcript strand, and the right-most end maps to the opposite strand.
fr-firststrand|	dUTP, NSR, NNSR|	Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.
fr-secondstrand|	Directional Illumina (Ligation), Standard SOLiD|	Same as above except we enforce the rule that the left-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during second strand synthesis is sequenced.
{: class="table"}

Please contact tophat.cufflinks@gmail.com to request support for a new protocol.

# Library Normalization Methods

You can control how library sizes (i.e. sequencing depths) are normalized in Cufflinks and Cuffdiff. Cuffdiff has several methods that require multiple libraries in order to work. Library normalization methods supported by Cufflinks work on one library at a time.

Normalization Method|	Supported by Cufflinks|	Supported by Cuffdiff|	Description
---|----|----
classic-fpkm|	Yes|	Yes|	Library size factor is set to 1 - no scaling applied to FPKM values or fragment counts. (default for Cufflinks)
geometric|	No|	Yes|	FPKMs and fragment counts are scaled via the median of the geometric means of fragment counts across all libraries, as described in Anders and Huber (Genome Biology, 2010). This policy identical to the one used by DESeq. (default for Cuffdiff)
quartile|	No|	Yes|	FPKMs and fragment counts are scaled via the ratio of the 75 quartile fragment counts to the average 75 quartile value across all libraries.
{: class="table"}


