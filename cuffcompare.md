---
layout: page
permalink: /cuffcompare/
description: "Documentation for Cuffcompare."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

<h1>Running Cuffcompare</h1>


Cufflinks includes a program that you can use to help analyze the transfrags you assemble. The program cuffcompare helps you:

 - Compare your assembled transcripts to a reference annotation
 - Track Cufflinks transcripts across multiple experiments (e.g. across a time course)

From the command line, run cuffcompare as follows:

**<tt>cuffcompare [options]* \<cuff1.gtf\> [cuff2.gtf] ... [cuffN.gtf]</tt>**


##Cuffcompare input files

Cuffcompare takes Cufflinks' GTF output as input, and optionally can take a "reference" annotation (such as from [Ensembl](ftp://ftp.ensembl.org/pub/current_gtf/) 

###Cuffcompare arguments

**<tt>\<cuff1.gtf\> [cuff2.gtf] ... [cuffN.gtf]</tt>**	

GTF files produced by cufflinks.

###Cuffcompare options	

**<tt>-h</tt>**

Prints the help message and exits

**<tt>-o \<outprefix\></tt>**

All output files created by Cuffcompare will have this prefix (e.g. <outprefix>.loci, <outprefix>.tracking, etc.). If this option is not provided the default output prefix being used is: "cuffcmp"

**<tt>-r</tt>**
An optional "reference" annotation GFF file. Each sample is matched against this file, and sample isoforms are tagged as overlapping, matching, or novel where appropriate. See the refmap and tmap output file descriptions below.

**<tt>-R</tt>**

If -r was specified, this option causes cuffcompare to ignore reference transcripts that are not overlapped by any transcript in one of cuff1.gtf,...,cuffN.gtf. Useful for ignoring annotated transcripts that are not present in your RNA-Seq samples and thus adjusting the "sensitivity" calculation in the accuracy report written in the <outprefix> file

**<tt>-s \<seq_dir\></tt>**

Causes cuffcompare to look into for fasta files with the underlying genomic sequences (one file per contig) against which your reads were aligned for some optional classification functions. For example, Cufflinks transcripts consisting mostly of lower-case bases are classified as repeats. Note that <seq_dir> must contain one fasta file per reference chromosome, and each file must be named after the chromosome, and have a .fa or .fasta extension.

**<tt>-C</tt>**

Enables the "contained" transcripts to be also written in the <outprefix>.combined.gtffile, with the attribute "contained_in" showing the first container transfrag found. By default, without this option, cuffcompare does not write in that file isoforms that were found to be fully contained/covered (with the same compatible intron structure) by other transfrags in the same locus.

**<tt>-V</tt>**

Cuffcompare is a little more verbose about what it's doing, printing messages to stderr, and it will also show warning messages about any inconsistencies or potential issues found while reading the given GFF file(s).

##Cuffcompare output files

Cuffcompare produces the following output files:

### Overall summary statistics: \<outprefix\>.stats

Cuffcompare reports various statistics related to the "accuracy" of the transcripts in each sample when compared to the reference annotation data. The typical gene finding measures of "sensitivity" and "specificity" (as defined in Burset, M., Guigó, R. : **Evaluation of gene structure prediction programs** (1996) Genomics, 34 (3), pp. 353-367. [doi: 10.1006/geno.1996.0298](http://dx.doi.org/10.1006/geno.1996.0298)) are calculated at various levels (nucleotide, exon, intron, transcript, gene) for each input file and reported in this file. The **Sn** and **Sp** columns show specificity and sensitivity values at each level, while the **fSn** and **fSp** columns are "fuzzy" variants of these same accuracy calculations, allowing for a very small variation in exon boundaries to still be counted as a "match". (If the -o option was not given the default prefix "cuffcmp" is used and these stats will be printed into a file named cuffcmp.stats in the current directory) 

### The "union" of all transfrags in all assemblies: \<outprefix\>.combined.gtf

Cuffcompare reports a GTF file containing the "union" of all transfrags in each sample. If a transfrag is present in both samples, it is thus reported once in the combined gtf. 

### Transfrags matching to each reference transcript: \<cuff_in\>.refmap
This tab delimited file lists, for each reference transcript, which cufflinks transcripts either fully or partially match it. There is one row per reference transcript, and the columns are as follows:

Column name|Example|Description
------|-----|-------
1|Reference gene name|	Myog|	The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.
2|Reference transcript id|	uc007crl.1|	The transcript_id attribute of the reference GTF record for this transcript
3|Class code|c|	The type of match between the Cufflinks transcripts in column 4 and the reference transcript. One of either 'c' for partial match, or '=' for full match.
4|Cufflinks matches|	CUFF.23567.0,CUFF.24689.0|	A comma separated list of Cufflinks transcript ids matching the reference transcript
{: class="table"}

### Best reference transcript for each transfrag: \<cuff_in\>.tmap

This tab delimited file lists the most closely matching reference transcript for each Cufflinks transcript. There is one row per Cufflinks transcript, and the columns are as follows:

Column number|Column name|Example|Description|
------|-----|-------|-----------
1|Reference gene name|	Myog|	The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.
2|	Reference transcript id|	uc007crl.1|	The transcript_id attribute of the reference GTF record for this transcript
3|	Class code|	c|	The type of relationship between the Cufflinks transcripts in column 4 and the reference transcript (as described in the Class Codes section below)
4|	Cufflinks gene id|	CUFF.23567|	The Cufflinks internal gene id
5|	Cufflinks transcript id|	CUFF.23567.0|	The Cufflinks internal transcript id
6|	Fraction of major isoform (FMI)|	100|	The expression of this transcript expressed as a fraction of the major isoform for the gene. Ranges from 1 to 100.
7|	FPKM|	1.4567|	The expression of this transcript expressed in FPKM
8|	FPKM_conf_lo|	0.7778|	The lower limit of the 95% FPKM confidence interval
9|	FPKM_conf_hi|	1.9776|	The upper limit of the 95% FPKM confidence interval
10|	Coverage|	3.2687|	The estimated average depth of read coverage across the transcript.
11|	Length|	1426|	The length of the transcript
12|	Major isoform ID|	CUFF.23567.0|	The Cufflinks ID of the gene's major isoform
{: class="table"}

### Tracking transfrags through multiple samples: \<outprefix\>.tracking

This file matches transcripts up between samples. Each row contains a transcript structure that is present in one or more input GTF files. Because the transcripts will generally have different IDs (unless you assembled your RNA-Seq reads against a reference transcriptome), cuffcompare examines the structure of each the transcripts, matching transcripts that agree on the coordinates and order of all of their introns, as well as strand. Matching transcripts are allowed to differ on the length of the first and last exons, since these lengths will naturally vary from sample to sample due to the random nature of sequencing.

If you ran <tt>cuffcompare</tt> with the <tt>-r</tt> option, the first and second columns contain the closest matching reference transcript to the one described by each row.

Here's an example of a line from the tracking file:

<pre>
TCONS_00000045 XLOC_000023 Tcea|uc007afj.1	j	\
     q1:exp.115|exp.115.0|100|3.061355|0.350242|0.350207 \
     q2:60hr.292|60hr.292.0|100|4.094084|0.000000|0.000000
</pre>

In this example, a transcript present in the two input files, called <tt>exp.115.0</tt> in the first and <tt>60hr.292.0</tt> in the second, doesn't match any reference transcript exactly, but shares exons with <tt>uc007afj.1</tt>, an isoform of the gene Tcea, as indicated by the [class code](#transfrag-class-codes) <tt>j</tt>. The first three columns are as follows:

Column number|Column name|Example|Description
------|-----|-------|-----------
1|Cufflinks transfrag id|TCONS_00000045|A unique internal id for the transfrag
2|Cufflinks locus id|XLOC_000023|A unique internal id for the locus
3|Reference gene id||Tcea|The gene_name attribute of the reference GTF record for this transcript, or '-' if no reference transcript overlaps this Cufflinks transcript
4|Reference transcript id|uc007afj.1|The transcript_id attribute of the reference GTF record for this transcript, or '-' if no reference transcript overlaps this Cufflinks transcript
5|Class code|c|The type of match between the Cufflinks transcripts in column 6 and the reference transcript. See [class codes](#transfrag-class-codes)
{: class="table"}

Each of the columns after the fifth have the following format:

<tt>qJ:\<gene_id\>|\<transcript_id\>|\<FMI\>|\<FPKM\>|\<conf_lo\>|\<conf_hi>\|\<cov\>|\<len\><tt>

A transcript need not be present in all samples to be reported in the tracking file. A sample not containing a transcript will have a "-" in its entry in the row for that transcript. 

(The following output files are created for each of the \<cuff_in\> file given, in the same directories where the \<cuff_in\> files reside) 


##Transfrag class codes

If you ran cuffcompare with the -r option, tracking rows will contain the following values. If you did not use -r, the rows will all contain "-" in their class code column.

Priority|Code|Description
-|-|-
1|	=|	Complete match of intron chain
2|	c|	Contained	
3|	j|	Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript	
4|	e|	Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.	
5|	i|	A transfrag falling entirely within a reference intron	
6|	o|	Generic exonic overlap with a reference transcript	
7|	p|	Possible polymerase run-on fragment (within 2Kbases of a reference transcript)	
8|	r|	Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case	
9|	u|	Unknown, intergenic transcript	
10|	x|	Exonic overlap with reference on the opposite strand	
11|	s|	An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)	
12|	.|	(.tracking file only, indicates multiple classifications)
{: class="table"}