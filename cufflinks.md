---
layout: page
permalink: /cufflinks/
description: "Documentation for Cufflinks."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

<h1>Running Cufflinks</h1>
Run cufflinks from the command line as follows:

**<tt>cufflinks [options] \<aligned_reads.(sam/bam)\></tt>**

The following is a detailed description of the options used to control Cufflinks:

###Arguments
	
**<tt><aligned_reads.(sam/bam)></tt>**

A file of RNA-Seq read alignments in the [SAM format](http://samtools.sourceforge.net/). SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cufflinks requires that the alignments you supply have some of these tags. Please see [Input formats](#cufflinks-input-files) for more details.|

###Cufflinks General Options

**<tt>-h/--help</tt>**

Prints the help message and exits

**<tt>-g/--GTF-guide \<reference_annotation.(gtf/gff)\></tt>**

Tells Cufflinks to use the supplied reference annotation [a GFF file](#gffgtf-files) to guide [RABT]({{ site.url }}/howitworks/index.html#how-does-reference-annotation-based-transcript-rabt-assembly-work) assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

**<tt>-M/--mask-file \<mask.(gtf/gff)\></tt>**

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

###Library Types

In cases where Cufflinks cannot determine the platform and protocol used to generate input reads, you can supply this information manually, which will allow Cufflinks to infer source strand information with certain protocols. The available options are listed below. For paired-end data, we currently only support protocols where reads are point towards each other.

Library Type|	Examples|	Description
---|----|----
fr-unstranded (default)|	Standard Illumina|	Reads from the left-most end of the fragment (in transcript coordinates) map to the transcript strand, and the right-most end maps to the opposite strand.
fr-firststrand|	dUTP, NSR, NNSR|	Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.
fr-secondstrand|	Directional Illumina (Ligation), Standard SOLiD|	Same as above except we enforce the rule that the left-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during second strand synthesis is sequenced.
{: class="table"}

Please contact tophat.cufflinks@gmail.com to request support for a new protocol.
