---
layout: page
permalink: /cuffdiff/
description: "Documentation for Cuffdiff."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

<h1>Running Cuffdiff</h1>

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

## Cuffdiff options:	

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

With this option, Cufflinks counts all fragments, including those not compatible with any reference transcript, towards the number of mapped fragments used in the FPKM denominator. It is inactive by default.

**<tt>--compatible-hits-norm</tt>**	

With this option, Cufflinks counts only those fragments compatible with some reference transcript towards the number of mapped fragments used in the FPKM denominator. Using this mode is generally recommended in Cuffdiff to reduce certain types of bias caused by differential amounts of ribosomal reads which can create the impression of falsely differentially expressed genes. It is active by default.

**<tt>-b/--frag-bias-correct \<genome.fa\></tt>**	

Providing Cufflinks with the multifasta file your reads were mapped to via this option instructs it to run our bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. See [How Cufflinks Works]({{ site.url}}/how_it_works/index.html)how_it_works/index.html#) for more details.

**<tt>-u/--multi-read-correct</tt>**	

Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. See [How Cufflinks Works]({{ site.url}}/how_it_works/index.html) for more details.

**<tt>-c/--min-alignment-count \<int\></tt>**	

The minimum number of alignments in a locus for needed to conduct significance testing on changes in that locus observed between samples. If no testing is performed, changes in the locus are deemed not signficant, and the locus' observed changes don't contribute to correction for multiple testing. The default is 10 fragment alignments.

**<tt>-M/--mask-file \<mask.(gtf/gff)\></tt>**	

Tells Cuffdiff to ignore all reads that could have come from transcripts in this GTF file. We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment methods and rRNA depletion kits, masking these transcripts often improves the overall robustness of transcript abundance estimates.

**<tt>--FDR \<float\></tt>**	

The allowed false discovery rate. The default is 0.05.

**<tt>--library-type</tt>**	

See [Library Types](#library-types)

**<tt>--library-norm-method</tt>**	

See [Library Normalization Methods](#library-normalization-methods)

**<tt>--dispersion-method</tt>**	

See [Cross-replicate dispersion estimation methods](#cross-replicate-dispersion-estimation-methods)

### Cuffdiff advanced options:	

**<tt>-m/--frag-len-mean \<int\></tt>**	

This is the expected (mean) fragment length. The default is 200bp. 

*Note: Cuffdiff now learns the fragment length mean for each SAM file, so using this option is no longer recommended with paired-end reads.*

**<tt>-s/--frag-len-std-dev \<int\></tt>**	

The standard deviation for the distribution on fragment lengths. The default is 80bp. 

*Note: Cuffdiff now learns the fragment length standard deviation for each SAM file, so using this option is no longer recommended with paired-end reads.*

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


### FPKM tracking files

Cuffdiff calculates the FPKM of each transcript, primary transcript, and gene in each sample. Primary transcript and gene FPKMs are computed by summing the FPKMs of transcripts in each primary transcript group or gene group. The results are output in FPKM tracking files in the format described [here](#fpkm-tracking-format). There are four FPKM tracking files:

|isoforms.fpkm_tracking|Transcript FPKMs|
genes.fpkm_tracking|Gene FPKMs. Tracks the summed FPKM of transcripts sharing each gene_id
cds.fpkm_tracking|Coding sequence FPKMs. Tracks the summed FPKM of transcripts sharing each p_id, independent of tss_id
tss_groups.fpkm_tracking|Primary transcript FPKMs. Tracks the summed FPKM of transcripts sharing each tss_id
{: class="table"}

### Count tracking files

Cuffdiff estimates the number of fragments that originated from each transcript, primary transcript, and gene in each sample. Primary transcript and gene counts are computed by summing the counts of transcripts in each primary transcript group or gene group. The results are output in count tracking files in the format described [here](#count-tracking-format). There are four Count tracking files:

isoforms.count_tracking|	Transcript counts
genes.count_tracking| Gene counts. Tracks the summed counts of transcripts sharing each gene_id
cds.count_tracking|	Coding sequence counts. Tracks the summed counts of transcripts sharing each p_id, independent of tss_id
tss_groups.count_tracking|	Primary transcript counts. Tracks the summed counts of transcripts sharing each tss_id
{: class="table"}

### Read group tracking files

Cuffdiff calculates the expression and fragment count for each transcript, primary transcript, and gene in each replicate. The results are output in per-replicate tracking files in the format described [here](#read-group-tracking-format). There are four read group tracking files:

isoforms.read_group_tracking|Transcript read group tracking
genes.read_group_tracking|Gene read group tracking. Tracks the summed expression and counts of transcripts sharing each gene_id in each replicate
cds.read_group_tracking|Coding sequence FPKMs. Tracks the summed expression and counts of transcripts sharing each p_id, independent of tss_id in each replicate
tss_groups.read_group_tracking|Primary transcript FPKMs. Tracks the summed expression and counts of transcripts sharing each tss_id in each replicate
{: class="table"}

### Differential expression tests

This tab delimited file lists the results of differential expression testing between samples for spliced transcripts, primary transcripts, genes, and coding sequences. Four files are created:

isoform_exp.diff|	Transcript-level differential expression.
gene_exp.diff|	Gene-level differential expression. Tests differences in the summed FPKM of transcripts sharing each gene_id
tss_group_exp.diff|	Primary transcript differential expression. Tests differences in the summed FPKM of transcripts sharing each tss_id
cds_exp.diff|	Coding sequence differential expression. Tests differences in the summed FPKM of transcripts sharing each p_id independent of tss_id
{: class="table"}

Each of the above files has the following format:

Column number|	Column name|	Example|	Description
---|---|---|---
1|	Tested id|	XLOC_000001|	A unique identifier describing the transcipt, gene, primary transcript, or CDS being tested
2|	gene|	Lypla1|	The gene_name(s) or gene_id(s) being tested
3|	locus|	chr1:4797771-4835363|	Genomic coordinates for easy browsing to the genes or transcripts being tested.
4|	sample 1|	Liver|	Label (or number if no labels provided) of the first sample being tested
5|	sample 2|	Brain|	Label (or number if no labels provided) of the second sample being tested
6|	Test status|	NOTEST|	Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.
7|	FPKMx|	8.01089|	FPKM of the gene in sample x
8|	FPKMy|	8.551545|	FPKM of the gene in sample y
9|	log2(FPKMy/FPKMx)|	0.06531|	The (base 2) log of the fold change y/x	
10|	test stat|	0.860902|	The value of the test statistic used to compute significance of the observed change in FPKM
11|	p| value	0.389292|	The uncorrected p-value of the test statistic
12|	q| value	0.985216|	The FDR-adjusted p-value of the test statistic
13|	significant|	no|	Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing
{: class="table"}

### Differential splicing tests - splicing.diff
This tab delimited file lists, for each primary transcript, the amount of isoform switching detected among its isoforms, i.e. how much differential splicing exists between isoforms processed from a single primary transcript. Only primary transcripts from which two or more isoforms are spliced are listed in this file.

Column number|	Column name|	Example|	Description|
---|---|---|---
1|	Tested id|	TSS10015|	A unique identifier describing the primary transcript being tested.
2|	gene name|	Rtkn|	The gene_name or gene_id that the primary transcript being tested belongs to
3|	locus|	chr6:83087311-83102572|	Genomic coordinates for easy browsing to the genes or transcripts being tested.
4|	sample 1|	Liver|	Label (or number if no labels provided) of the first sample being tested
5|	sample 2|	Brain|	Label (or number if no labels provided) of the second sample being tested
6|	Test status|	OK|	Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.
7|	Reserved|	0|	
8|	Reserved|	0|	
9|	√JS(x,y)|	0.22115|	The amount of isoform switching between the isoforms originating from this TSS, as measured by the square root of the Jensen-Shannon divergence computed on the relative abundances of the splice variants	
10|	test stat|	0.22115|	The value of the test statistic used to compute significance of the observed overloading, equal to √JS(x,y)
11|	p value|	0.000174982|	The uncorrected p-value of the test statistic.
12|	q value|	0.985216|	The FDR-adjusted p-value of the test statistic
13|	significant|	no|	Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing
{: class="table"}

### Differential coding output - cds.diff

This tab delimited file lists, for each gene, the amount of overloading detected among its coding sequences, i.e. how much differential CDS output exists between samples. Only genes producing two or more distinct CDS (i.e. multi-protein genes) are listed here.

Column number|	Column name|	Example|	Description
---|---|---|---
1|	Tested id|	XLOC_000002 |	A unique identifier describing the gene being tested.
2|	gene name|	Atp6v1h|	The gene_name or gene_id
3|	locus|	chr1:5073200-5152501|	Genomic coordinates for easy browsing to the genes or transcripts being tested.
4|	sample 1|	Liver|	Label (or number if no labels provided) of the first sample being tested
5|	sample 2|	Brain|	Label (or number if no labels provided) of the second sample being tested
6|	Test status|	OK|	Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.
7|	Reserved|	0|	
8|	Reserved|	0|	
9|	√JS(x,y)|	0.0686517|	The CDS overloading of the gene, as measured by the square root of the Jensen-Shannon divergence computed on the relative abundances of the coding sequences	
10|	test stat|	0.0686517|	The value of the test statistic used to compute significance of the observed overloading, equal to √JS(x,y)
11|	p value|	0.00546783|	The uncorrected p-value of the test statistic
12|	q value|	0.985216|	The FDR-adjusted p-value of the test statistic
13|	significant|	no|	Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing
{: class="table"}

### Differential promoter use - promoters.diff

This tab delimited file lists, for each gene, the amount of overloading detected among its primary transcripts, i.e. how much differential promoter use exists between samples. Only genes producing two or more distinct primary transcripts (i.e. multi-promoter genes) are listed here.

Column number|	Column name|	Example|	Description
---|---|---|---
1|	Tested id|	XLOC_000019|	A unique identifier describing the gene being tested.
2|	gene name|	Tmem70|	The gene_name or gene_id
3|	locus|	chr1:16651657-16668357|	Genomic coordinates for easy browsing to the genes or transcripts being tested.
4|	sample 1|	Liver|	Label (or number if no labels provided) of the first sample being tested
5|	sample 2|	Brain|	Label (or number if no labels provided) of the second sample being tested
6|	Test status|	OK|	Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.
7|	Reserved|	0|	
8|	Reserved|	0|	
9|	√JS(x,y)|	0.0124768|	The promoter overloading of the gene, as measured by the square root of the Jensen-Shannon divergence computed on the relative abundances of the primary transcripts	
10|	test stat|	0.0124768|	The value of the test statistic used to compute significance of the observed overloading, equal to √JS(x,y)
11|	p value|	0.394327|	The uncorrected p-value of the test statistic
12|	q value|	0.985216|	The FDR-adjusted p-value of the test statistic
13|	significant|	no|	Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing
{: class="table"}

### Read group info - read_groups.info

This tab delimited file lists, for each replicate, key properties used by Cuffdiff during quantification, such as library normalization factors. The read_groups.info file has the following format:

Column number|	Column name|	Example|	Description
---|---|---|---
1|	file|	mCherry_rep_A/accepted_hits.bam|	BAM or SAM file containing the data for the read group
2|	condition|	mCherry|	Condition to which the read group belongs
3|	replicate_num|	0|	Replicate number of the read group
4|	total_mass|	4.72517e+06|	Total number of fragments for the read group
5|	norm_mass|	4.72517e+06|	Fragment normalization constant used during calculation of FPKMs.
6|	internal_scale|	1.23916|	Scaling factor used to normalize for library size
7|	external_scale|	1.0|	Currently unused, and always equal to 1.0. 
{: class="table"}

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

#Cross-replicate dispersion estimation methods

Cuffdiff works by modeling the variance in fragment counts across replicates as a function of the mean fragment count across replicates. Strictly speaking, models a quantitity called dispersion - the variance present in a group of samples beyond what is expected from a simple Poisson model of RNA_Seq. You can control how Cuffdiff constructs its model of dispersion in locus fragment counts. Each condition that has replicates can receive its own model, or Cuffdiff can use a global model for all conditions. All of these policies are identical to those used by DESeq (Anders and Huber, Genome Biology, 2010)

Dispersion Method|	Description
----|----
pooled|	Each replicated condition is used to build a model, then these models are averaged to provide a single global model for all conditions in the experiment. (Default)
per-condition|	Each replicated condition receives its own model. Only available when all conditions have replicates.
blind|	All samples are treated as replicates of a single global "condition" and used to build one model.
poisson|	The Poisson model is used, where the variance in fragment count is predicted to equal the mean across replicates. Not recommended.
{: class="table"}

Which method you choose largely depends on whether you expect variability in each group of samples to be similar. For example, if you are comparing two groups, A and B, where A has low cross-replicate variability and B has high variability, it may be best to choose per-condition. However, if the conditions have similar levels of variability, you might stick with the default, which sometimes provides a more robust model, especially in cases where each group has few replicates. Finally, if you only have a single replicate in each condition, you must use blind, which treats all samples in the experiment as replicates of a single condition. This method works well when you expect the samples to have very few differentially expressed genes. If there are many differentially expressed genes, Cuffdiff will construct an overly conservative model and you may not get any significant calls. In this case, you will need more replicates in your experiment.
