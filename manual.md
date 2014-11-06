---
layout: page
permalink: /manual/
description: "Documentation for Cufflinks."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}


#Prerequisites

Cufflinks runs on intel-based computers running Linux or Mac OS X and that have GCC 4.0 or greater installed. You can install pre-compiled binaries or build Cufflinks from the source code. If you wish to build Cufflinks yourself, you will need to install the [Boost C++ libraries](http://www.boost.org/). See [Installing Boost]({{ site.url }}/getting_started/index.html#boost), on the [getting started]({{ site.url }}/getting_started/index.html) page. You will also need to build and install the [SAM tools](http://samtools.sourceforge.net/), but you should take a look at the getting started page for detailed instructions, because the headers and libbam must be accessible to the Cufflinks build scripts.

#Running Cufflinks
Run cufflinks from the command line as follows:

**<tt>cufflinks [options] <aligned_reads.(sam/bam)></tt>**

The following is a detailed description of the options used to control Cufflinks:

###Arguments
	
**<tt><aligned_reads.(sam/bam)></tt>**

A file of RNA-Seq read alignments in the [SAM format](http://samtools.sourceforge.net/). SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cufflinks requires that the alignments you supply have some of these tags. Please see [Input formats](#cufflinks-input-files) for more details.|

###Cufflinks General Options

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-g/--GTF-guide <reference_annotation.(gtf/gff)></tt>**

Tells Cufflinks to use the supplied reference annotation [a GFF file](#gffgtf-files) to guide [RABT]({{ site.url }}/howitworks/index.html#how-does-reference-annotation-based-transcript-rabt-assembly-work) assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

**<tt>-M/--mask-file <mask.(gtf/gff)></tt>**

Tells Cufflinks to ignore all reads that could have come from transcripts in this GTF file. We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment methods and rRNA depletion kits, masking these transcripts often improves the overall robustness of transcript abundance estimates.

**<tt>-b/--frag-bias-correct \<genome.fa\></tt>**

Providing Cufflinks with a multifasta file via this option instructs it to run our new bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. See [How Cufflinks Works]({{ site.url }}/how_it_works/index.html) for more details.

**<tt>-u/--multi-read-correct</tt>**	

Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. See [How Cufflinks Works]({{ site.url }}/how_it_works/index.html) for more details.

**<tt>--library-type</tt>**

See [Library Types](#library-types)

**<tt>--library-norm-method</tt>**	

See [Library Normalization Methods](#library-normalization-methods)

###Advanced Abundance Estimation Options

**<tt>-m/--frag-len-mean \<int\></tt>**	

This is the expected (mean) fragment length. The default is 200bp. 
*Note: Cufflinks now learns the fragment length mean for each SAM file, so using this option is no longer recommended with paired-end reads.*

**<tt>-s/--frag-len-std-dev \<int\></tt>**

The standard deviation for the distribution on fragment lengths. The default is 80bp. 
*Note: Cufflinks now learns the fragment length standard deviation for each SAM file, so using this option is no longer recommended with paired-end reads.*

**<tt>-N/--upper-quartile-norm</tt>**

With this option, Cufflinks normalizes by the upper quartile of the number of fragments mapping to individual loci instead of the total number of sequenced fragments. This can improve robustness of differential expression calls for less abundant genes and transcripts.

**<tt>--total-hits-norm</tt>**

With this option, Cufflinks counts all fragments, including those not compatible with any reference transcript, towards the number of mapped hits used in the FPKM denominator. This option can be combined with -N/--upper-quartile-norm. It is active by default.

**<tt>--compatible-hits-norm</tt>**

With this option, Cufflinks counts only those fragments compatible with some reference transcript towards the number of mapped hits used in the FPKM denominator. This option can be combined with -N/--upper-quartile-norm. It is inactive by default, and can only be used in combination with 

**<tt>--GTF</tt>**

Use with either RABT or ab initio assembly is not supported

**<tt>--max-mle-iterations \<int\></tt>**

Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000

**<tt>--max-bundle-frags \<int></tt>**

Sets the maximum number of fragments a locus may have before being skipped. Skipped loci are listed in skipped.gtf. Default: 1000000

**<tt>--no-effective-length-correction</tt>**

Cufflinks will not employ its "effective" length normalization to transcript FPKM.

**<tt>--no-length-correction</tt>**

Cufflinks will not normalize fragment counts by transcript length at all. Use this option when fragment count is independent of the size of the features being quantified (e.g. for small RNA libraries, where no fragmentation takes place, or 3 prime end sequencing, where sampled RNA fragments are all essentially the same length). Experimental option, use with caution.

###Advanced Assembly Options	

**<tt>-L/--label</tt>**

Cufflinks will report transfrags in GTF format, with a prefix given by this option. The default prefix is "CUFF".

**<tt>-F/--min-isoform-fraction \<0.0-1.0\></tt>**

After calculating isoform abundance for a gene, Cufflinks filters out transcripts that it believes are very low abundance, because isoforms expressed at extremely low levels often cannot reliably be assembled, and may even be artifacts of incompletely spliced precursors of processed transcripts. This parameter is also used to filter out introns that have far fewer spliced alignments supporting them. The default is 0.1, or 10% of the most abundant isoform (the major isoform) of the gene.

**<tt>-j/--pre-mrna-fraction \<0.0-1.0\></tt>**

Some RNA-Seq protocols produce a significant amount of reads that originate from incompletely spliced transcripts, and these reads can confound the assembly of fully spliced mRNAs. Cufflinks uses this parameter to filter out alignments that lie within the intronic intervals implied by the spliced alignments. The minimum depth of coverage in the intronic region covered by the alignment is divided by the number of spliced reads, and if the result is lower than this parameter value, the intronic alignments are ignored. The default is 15%.

**<tt>-I/--max-intron-length \<int\></tt>**

The maximum intron length. Cufflinks will not report transcripts with introns longer than this, and will ignore SAM alignments with REF_SKIP CIGAR operations longer than this. The default is 300,000.

**<tt>-a/--junc-alpha \<0.0-1.0\></tt>**

The alpha value for the binomial test used during false positive spliced alignment filtration. Default: 0.001

**<tt>-A/--small-anchor-fraction \<0.0-1.0\></tt>**

Spliced reads with less than this percent of their length on each side of the junction are considered suspicious and are candidates for filtering prior to assembly. Default: 0.09.

**<tt>--min-frags-per-transfrag \<int\></tt>**

Assembled transfrags supported by fewer than this many aligned RNA-Seq fragments are not reported. Default: 10.

**<tt>--overhang-tolerance \<int\></tt>**

The number of bp allowed to enter the intron of a transcript when determining if a read or another transcript is mappable to/compatible with it. The default is 8 bp based on the default bowtie/TopHat parameters.

**<tt>--max-bundle-length \<int\></tt>**

Maximum genomic length allowed for a given bundle. The default is 3,500,000 bp.

**<tt>--min-intron-length \<int\></tt>**

Minimum intron size allowed in genome. The default is 50 bp.

**<tt>--trim-3-avgcov-thresh \<int\></tt>**

Minimum average coverage required to attempt 3' trimming. The default is 10.

**<tt>--trim-3-dropoff-frac \<int\></tt>**

The fraction of average coverage below which to trim the 3' end of an assembled transcript. The default is 0.1.

**<tt>--max-multiread-fraction \<0.0-1.0\></tt>**

The fraction a transfrag's supporting reads that may be multiply mapped to the genome. A transcript composed of more than this fraction will not be reported by the assembler. Default: 0.75 (75% multireads or more is suppressed).

**<tt>--overlap-radius \<int\></tt>**

Transfrags that are separated by less than this distance get merged together, and the gap is filled. Default: 50 (in bp).

###Advanced Reference Annotation Based Transcript (RABT) Assembly Options

*Note: These options have an affect only when used in conjuction with <tt>-g/--GTF-guide</tt>.*

**<tt>--3-overhang-tolerance \<int\></tt>**	

The number of bp allowed to overhang the 3' end of a reference transcript when determining if an assembled transcript should be merged with it (ie, the assembled transcript is not novel). The default is 600 bp.

**<tt>--intron-overhang-tolerance \<int\></tt>**

The number of bp allowed to enter the intron of a reference transcript when determining if an assembled transcript should be merged with it (ie, the assembled transcript is not novel). The default is 50 bp.

**<tt>--no-faux-reads</tt>**

This option disables tiling of the reference transcripts with faux reads. Use this if you only want to use sequencing reads in assembly but do not want to output assembled transcripts that lay within reference transcripts. All reference transcripts in the input annotation will also be included in the output.

###Advanced Program Behavior Options
	
**<tt>-v/--verbose</tt>**

Print lots of status updates and other diagnostic information.

**<tt>-q/--quiet</tt>**

Suppress messages other than serious warnings and errors.

**<tt>--no-update-check</tt>**

Turns off the automatic routine that contacts the Cufflinks server to check for a more recent version.

##Cufflinks Input Files

Cufflinks takes a text file of SAM alignments, or a binary SAM (BAM) file as input. For more details on the SAM format, see the specification. The RNA-Seq read mapper TopHat produces output in this format, and is recommended for use with Cufflinks. However Cufflinks will accept SAM alignments generated by any read mapper. Here's an example of an alignment Cufflinks will accept:

<pre>
s6.25mer.txt-913508	16	chr1 4482736 255 14M431N11M * 0 0 \
   CAAGATGCTAGGCAAGTCTTGGAAG IIIIIIIIIIIIIIIIIIIIIIIII NM:i:0 XS:A:-
</pre>

Note the use of the custom tag XS. This attribute, which must have a value of "+" or "-", indicates which strand the RNA that produced this read came from. While this tag can be applied to any alignment, including unspliced ones, it must be present for all spliced alignment records (those with a 'N' operation in the CIGAR string).
The SAM file supplied to Cufflinks must be sorted by reference position. If you aligned your reads with TopHat, your alignments will be properly sorted already. If you used another tool, you may want to make sure they are properly sorted as follows:

{% highlight bash %}
sort -k 3,3 -k 4,4n hits.sam > hits.sam.sorted
{% endhighlight %}


##Cufflinks Output Files

Cufflinks produces three output files:

###Transcriptome assembly: transcripts.gtf

This GTF file contains Cufflinks' assembled isoforms. The first 7 columns are [standard GTF](#gffgtf-files), and the last column contains attributes, some of which are also standardized ("gene_id", and "transcript_id"). There one GTF record per row, and each record represents either a transcript or an exon within a transcript. The columns are defined as follows:

|Column number|Column name|Example|Description
|:-:|-|-|-|
|1|seqname|chrX|Chromosome or contig name
|2|source|Cufflinks|The name of the program that generated this file (always 'Cufflinks')
|3|feature|exon|The type of record (always either "transcript" or "exon".
|4|start|77696957|The leftmost coordinate of this record (where 1 is the leftmost possible coordinate)
|5|end|77712009|The rightmost coordinate of this record, inclusive.
|6|score|77712009|The most abundant isoform for each gene is assigned a score of 1000. Minor isoforms are scored by the ratio (minor FPKM/major FPKM)
|7|strand|+|Cufflinks' guess for which strand the isoform came from. Always one of "+", "-", "."
|7|frame|.|Cufflinks does not predict where the start and stop codons (if any) are located within each transcript, so this field is not used.|
|8|attributes|...|See below.
{: class="table"}


Each GTF record is decorated with the following attributes:

|Attribute|Example|Description|
|:-------|--------|------|
|gene_id|CUFF.1|Cufflinks gene id|
|transcript_id|CUFF.1.1|Cufflinks transcript id
|FPKM|101.267|Isoform-level relative abundance in Fragments Per Kilobase of exon model per Million mapped fragments
|frac|0.7647|Reserved. Please ignore, as this attribute may be deprecated in the future
|conf_lo|0.07|Lower bound of the 95% confidence interval of the abundance of this isoform, as a fraction of the isoform abundance. That is, lower bound = FPKM * (1.0 - conf_lo)
|conf_hi|0.1102|Upper bound of the 95% confidence interval of the abundance of this isoform, as a fraction of the isoform abundance. That is, upper bound = FPKM * (1.0 + conf_lo)
|cov|100.765|Estimate for the absolute depth of read coverage across the whole transcript
|full_read_support|yes|When RABT assembly is used, this attribute reports whether or not all introns and internal exons were fully covered by reads from the data.
{: class="table"}

### Transcript-level expression: isoforms.fpkm_tracking

This file contains the estimated isoform-level expression values in the generic [FPKM Tracking Format](#fpkm-tracking-format). Note, however that as there is only one sample, the "q" format is not used.

###Gene-level expression: genes.fpkm_tracking

This file contains the estimated gene-level expression values in the generic [FPKM Tracking Format](#fpkm-tracking-format). Note, however that as there is only one sample, the "q" format is not used.

#Running Cuffcompare

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

#Merging assemblies with cuffmerge

Cufflinks includes a script called cuffmerge that you can use to merge together several Cufflinks assemblies. It handles also handles running Cuffcompare for you, and automatically filters a number of transfrags that are probably artfifacts. If you have a reference GTF file available, you can provide it to the script in order to gracefully merge novel isoforms and known isoforms and maximize overall assembly quality. The main purpose of this script is to make it easier to make an assembly GTF file suitable for use with Cuffdiff. From the command line, run <tt>cuffmerge</tt> as follows:

**<tt>cuffmerge [options]* \<assembly_GTF_list.txt\></tt>**

##Cuffmerge input files

cuffmerge takes several assembly GTF files from Cufflinks' as input. Input GTF files must be specified in a "manifest" file listing full paths to the files. 

###Cuffmerge arguments	

**<tt>\<assembly_list.txt\></tt>**

Text file "manifest" with a list (one per line) of GTF files that you'd like to merge together into a single GTF file.

###Cuffmerge options	

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-o \<outprefix\></tt>**

Write the summary stats into the text output file \<outprefix\>(instead of stdout)

**<tt>-g/--ref-gtf</tt>**

An optional "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.

**<tt>-p/--num-threads \<int\></tt>**

Use this many threads to align reads. The default is 1.

**<tt>-s/--ref-sequence \<seq_dir\>/\<seq_fasta\></tt>**

This argument should point to the genomic DNA sequences for the reference. If a directory, it should contain one fasta file per contig. If a multifasta file, all contigs should be present. The merge script will pass this option to <tt>cuffcompare</tt>, which will use the sequences to assist in classifying transfrags and excluding artifacts (e.g. repeats). For example, Cufflinks transcripts consisting mostly of lower-case bases are classified as repeats. 

*Note that \<seq_dir\> **must** contain one fasta file per reference chromosome, and each file must be named after the chromosome, and have a <tt>.fa</tt> or <tt>.fasta</tt> extension*.

##Cuffmerge output files

Cuffmerge produces a GTF file, <tt>merged.gtf</tt> that merges together the input assemblies.

#Running Cuffquant

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

Tells Cuffquant to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. See How Cufflinks Works for more details.

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

#Running Cuffdiff

Cufflinks includes a program, "Cuffdiff", that you can use to find significant changes in transcript expression, splicing, and promoter use. From the command line, run cuffdiff as follows:

**<tt>cuffdiff [options]* \<transcripts.gtf\> \\</tt>** 

**<tt>\<sample1_replicate1.sam[,...,sample1_replicateM.sam]\> \\</tt>**

**<tt>\<sample2_replicate1.sam[,...,sample2_replicateM.sam]\> ... \\</tt>**

**<tt>[sampleN.sam_replicate1.sam[,...,sample2_replicateM.sam]]</tt>**

###Cuffdiff arguments	

**<tt>\<transcripts.(gtf/gff)\></tt>**

A transcript annotation file produced by cufflinks, cuffcompare, or other source.

**<tt>\<sample1.(sam/bam/cxb)\></tt>**

A SAM file of aligned RNA-Seq reads. If more than two are provided, Cuffdiff tests for differential expression and regulation between all pairs of samples.

Options:	

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-o/--output-dir \<string\></tt>**	

Sets the name of the directory in which Cuffdiff will write all of its output. The default is "./".

**<tt>-L/--labels \<label1,label2,...,labelN\></tt>**	

Specify a label for each sample, which will be included in various output files produced by Cuffdiff.

**<tt>-p/--num-threads \<int\></tt>**	

Use this many threads to align reads. The default is 1.

**<tt>-T/--time-series</tt>**	

Instructs Cuffdiff to analyze the provided samples as a time series, rather than testing for differences between all pairs of samples. Samples should be provided in increasing time order at the command line (e.g first time point SAM, second timepoint SAM, etc.)

**<tt>--total-hits-norm</tt>**	

With this option, Cufflinks counts all fragments, including those not compatible with any reference transcript, towards the number of mapped fragments used in the FPKM denominator. This option can be combined with -N/--upper-quartile-norm. It is inactive by default.

**<tt>--compatible-hits-norm</tt>**	

With this option, Cufflinks counts only those fragments compatible with some reference transcript towards the number of mapped fragments used in the FPKM denominator. This option can be combined with -N/--upper-quartile-norm. Using this mode is generally recommended in Cuffdiff to reduce certain types of bias caused by differential amounts of ribosomal reads which can create the impression of falsely differentially expressed genes. It is active by default.

**<tt>-b/--frag-bias-correct \<genome.fa\></tt>**	

Providing Cufflinks with the multifasta file your reads were mapped to via this option instructs it to run our bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. See How Cufflinks Works for more details.

**<tt>-u/--multi-read-correct</tt>**	

Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. See How Cufflinks Works for more details.

**<tt>-c/--min-alignment-count \<int\></tt>**	

The minimum number of alignments in a locus for needed to conduct significance testing on changes in that locus observed between samples. If no testing is performed, changes in the locus are deemed not signficant, and the locus' observed changes don't contribute to correction for multiple testing. The default is 10 fragment alignments.

**<tt>-M/--mask-file \<mask.(gtf/gff)\></tt>**	

Tells Cuffdiff to ignore all reads that could have come from transcripts in this GTF file. We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment methods and rRNA depletion kits, masking these transcripts often improves the overall robustness of transcript abundance estimates.

**<tt>--FDR \<float\></tt>**	

The allowed false discovery rate. The default is 0.05.

**<tt>--library-type</tt>**	

See Library Types

**<tt>--library-norm-method</tt>**	

See Library Normalization Methods

**<tt>--dispersion-method</tt>**	

See Cross-replicate dispersion estimation methods

Advanced Options:	

**<tt>-m/--frag-len-mean \<int\></tt>**	

This is the expected (mean) fragment length. The default is 200bp. 
Note: Cufflinks now learns the fragment length mean for each SAM file, so using this option is no longer recommended with paired-end reads.

**<tt>-s/--frag-len-std-dev \<int\></tt>**	

The standard deviation for the distribution on fragment lengths. The default is 80bp. 
Note: Cufflinks now learns the fragment length standard deviation for each SAM file, so using this option is no longer recommended with paired-end reads.

**<tt>--max-mle-iterations \<int\></tt>**	

Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000

**<tt>-v/--verbose</tt>**	

Print lots of status updates and other diagnostic information.

**<tt>-q/--quiet</tt>**	

Suppress messages other than serious warnings and errors.

**<tt>--no-update-check</tt>**	

Turns off the automatic routine that contacts the Cufflinks server to check for a more recent version.

**<tt>--poisson-dispersion</tt>**	

Use the Poisson fragment dispersion model instead of learning one in each condition.

**<tt>--emit-count-tables</tt>**	

Cuffdiff will output a file for each condition (called \<sample\>_counts.txt) containing the fragment counts, fragment count variances, and fitted variance model. For internal debugging only. This option will be removed in a future version of Cuffdiff.

**<tt>-F/--min-isoform-fraction \<0.0-1.0\></tt>**	

Cuffdiff will round down to zero the abundance of alternative isoforms quantified at below the specified fraction of the major isoforms. This is done after MLE estimation but before MAP estimation to improve robustness of confidence interval generation and differential expression analysis. The default is 1e-5, and we recommend you not alter this parameter.

**<tt>--max-bundle-frags \<int\></tt>**	

Sets the maximum number of fragments a locus may have before being skipped. Skipped loci are marked with status HIDATA. Default: 1000000

**<tt>--max-frag-count-draws \<int\></tt>**	

Cuffdiff will make this many draws from each transcript's predicted negative binomial random numbder generator. Each draw is a number of fragments that will be probabilistically assigned to the transcripts in the transcriptome. Used to estimate the variance-covariance matrix on assigned fragment counts. Default: 100.

**<tt>--max-frag-assign-draws \<int\></tt>**	

For each fragment drawn from a transcript, Cuffdiff will assign it this many times (probabilistically), thus estimating the assignment uncertainty for each transcript. Used to estimate the variance-covariance matrix on assigned fragment counts. Default: 50.

**<tt>--min-reps-for-js-test \<int\></tt>**	

Cuffdiff won't test genes for differential regulation unless the conditions in question have at least this many replicates. Default: 3.

**<tt>--no-effective-length-correction	

Cuffdiff will not employ its "effective" length normalization to transcript FPKM.

**<tt>--no-length-correction</tt>**

Cuffdiff will not normalize fragment counts by transcript length at all. Use this option when fragment count is independent of the size of the features being quantified (e.g. for small RNA libraries, where no fragmentation takes place, or 3 prime end sequencing, where sampled RNA fragments are all essentially the same length). Experimental option, use with caution.


##Cuffdiff input Files

Cuffdiff takes a GTF2/GFF3 file of transcripts as input, along with two or more SAM files containing the fragment alignments for two or more samples. It produces a number of output files that contain test results for changes in expression at the level of transcripts, primary transcripts, and genes. It also tracks changes in the relative abundance of transcripts sharing a common transcription start site, and in the relative abundances of the primary transcripts of each gene. Tracking the former allows one to see changes in splicing, and the latter lets one see changes in relative promoter use within a gene.

If you have more than one **replicate** for a sample, supply the SAM files for the sample as a single **comma-separated** list. It is not necessary to have the same number of replicates for each sample.

Note that Cuffdiff can also accepted BAM files (which are binary, compressed SAM files). It can also accept CXB files produced by Cuffquant. Note that mixing SAM and BAM files is supported, but you cannot currently mix CXB and SAM/BAM. If one of the samples is supplied as a CXB file, all of the samples must be supplied as CXB files.

Cuffdiff requires that transcripts in the input GTF be annotated with certain attributes in order to look for changes in primary transcript expression, splicing, coding output, and promoter use. These attributes are:

**<tt>tss_id</tt>**

The ID of this transcript's inferred start site. Determines which primary transcript this processed transcript is believed to come from. Cuffcompare appends this attribute to every transcript reported in the .combined.gtf file.

**<tt>p_id</tt>**	

The ID of the coding sequence this transcript contains. This attribute is attached by Cuffcompare to the .combined.gtf records only when it is run with a reference annotation that include CDS records. Further, differential CDS analysis is only performed when all isoforms of a gene have p_id attributes, because neither Cufflinks nor Cuffcompare attempt to assign an open reading frame to transcripts.

*Note: If an arbitrary GTF/GFF3 file is used as input (instead of the .combined.gtf file produced by Cuffcompare), these attributes will not be present, but Cuffcompare can still be used to obtain these attributes with a command like this:* 

{% highlight bash %}
cuffcompare -s /path/to/genome_seqs.fa -CG -r annotation.gtf annotation.gtf
{% endhighlight %}

The resulting cuffcmp.combined.gtf file created by this command will have the <tt>tss_id</tt> and <tt>p_id</tt> attributes added to each record and this file can be used as input for cuffdiff.

##Cuffdiff output Files

###FPKM tracking

###Count tracking

###Read group tracking

###Differential expression

###Differential splicing

###Differential coding sequence output

###Differential promoter use

###Read group info

###Run info

#Running Cuffnorm

##Cuffnorm input Files

##Cuffnorm output Files

##Sample sheets for Cuffdiff and Cuffnorm

#Contrast files for Cuffdiff

#Output formats

##FPKM Tracking format

##Count Tracking format

##Read Group Tracking format

##Simple-table expression format

##Simple-table gene attributes format

##Simple-table sample attributes format

#Library Types

#Library Normalization Methods

#Cross-replicate dispersion estimation methods

#GFF/GTF Files

