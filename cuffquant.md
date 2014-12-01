---
layout: page
permalink: /cuffquant/
description: "Documentation for Cuffquant."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

<h1>Running Cuffquant</h1>

Run cuffquant from the command line as follows:

**<tt>cuffquant [options]* \<annotation.(gtf/gff)\> \<aligned_reads.(sam/bam)\></tt>**

### Cuffquant arguments

**<tt>\<annotation.(gtf/gff)\></tt>**	

Tells Cuffquant to use the supplied reference annotation ([a GFF file](#gffgtf-files)) to estimate isoform expression. The program will ignore alignments not structurally compatible with any reference transcript.

**<tt>\<aligned_reads.(sam/bam)\></tt>**

A file of RNA-Seq read alignments in the SAM format. SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cuffquant requires that the alignments you supply have some of these tags. Please see Input formats for more details.

### Cuffquant general options

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-o/--output-dir \<string\></tt>**

Sets the name of the directory in which Cuffquant will write all of its output. The default is "./".

**<tt>-p/--num-threads \<int\></tt>**

Use this many threads to align reads. The default is 1.

**<tt>-M/--mask-file \<mask.(gtf/gff)\></tt>**

Tells Cuffquant to ignore all reads that could have come from transcripts in this GTF file. We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment methods and rRNA depletion kits, masking these transcripts often improves the overall robustness of transcript abundance estimates.

**<tt>-b/--frag-bias-correct \<genome.fa\></tt>**

Providing Cuffquant with a multifasta file via this option instructs it to run bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. See [How Cufflinks Works]({{ site.url}}/how_it_works/index.html) for more details.

**<tt>-u/--multi-read-correct</tt>**

Tells Cuffquant to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. See [How Cufflinks Works]({{ site.url}}/how_it_works/index.html) for more details.

**<tt>--library-type</tt>**

See [Library Types](#library-types)

### Cuffquant advanced abundance estimation options

**<tt>-m/--frag-len-mean \<int\></tt>**

This is the expected (mean) fragment length. The default is 200bp. 

*Note: Cuffquant learns the fragment length mean for each SAM file, so using this option is no longer recommended with paired-end reads.*

**<tt>-s/--frag-len-std-dev \<int\></tt>**

The standard deviation for the distribution on fragment lengths. The default is 80bp. 

*Note: Cuffquant learns the fragment length standard deviation for each SAM file, so using this option is no longer recommended with paired-end reads.*

**<tt>--max-mle-iterations \<int\></tt>**

Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000

**<tt>--max-bundle-frags \<int\></tt>**

Sets the maximum number of fragments a locus may have before being skipped. Default: 1000000

**<tt>--no-effective-length-correction</tt>**

Cuffquant will not employ its "effective" length normalization to transcript FPKM.

**<tt>--no-length-correction</tt>**

Cuffquant will not normalize fragment counts by transcript length at all. Use this option when fragment count is independent of the size of the features being quantified (e.g. for small RNA libraries, where no fragmentation takes place, or 3 prime end sequencing, where sampled RNA fragments are all essentially the same length). Experimental option, use with caution.


###Advanced Program Behavior Options:	

**<tt>-v/--verbose</tt>**

Print lots of status updates and other diagnostic information.

**<tt>-q/--quiet</tt>**

Suppress messages other than serious warnings and errors.

**<tt>--no-update-check</tt>**

Turns off the automatic routine that contacts the Cufflinks server to check for a more recent version.

##Cuffquant input

Cuffquant takes as input a single SAM/BAM file of aligned reads and a single GTF/GFF file of gene annotations.

##Cuffquant output

Cuffquant produces writes a single output file, abundances.cxb, into the output directory. CXB files are binary files, and can be passed to Cuffnorm or Cuffdiff for further processing.
