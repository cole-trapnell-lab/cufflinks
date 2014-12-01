---
layout: page
permalink: /file_formats/
description: "Documentation for the file formats used by the Cufflinks suite."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}


#GFF/GTF Files

Some of the Cufflinks modules take as input a file (or more) containing known gene annotations or other transcript data in GFF format (General Feature Format). GFF has many versions, but the two most popular that are supported by Cufflinks (and other programs in the Tuxedo suite, like [Tophat](http://tophat.cbcb.umd.edu)) are GTF2 (Gene Transfer Format, described [here](http://mblab.wustl.edu/GTF2.html)) and GFF3 (defined [here](http://www.sequenceontology.org/gff3.shtml)). Here are a few notes about the way these formats are interpreted by the Cufflinks programs.

## GTF2
As seen in the GTF2 specification, the **transcript_id** attribute is also required by our GFF parser, and a **gene_id** attribute, though not strictly required in our programs, is very useful for grouping alternative transcripts under a gene/locus identifier. An optional gene_name attribute, if found, will be taken and shown as  a symbolic gene name or short-form abbreviation (e.g. gene symbols from HGNC or Entrez Gene). Some annotation sources (e.g. Ensembl) place a "human readable" gene name/symbol in the **gene_name** attribute, like a HUGO symbol (while gene_id might be just an automatically generated numeric identifier for the gene). 

TopHat and Cufflinks generally expect exon features to define a transcript structure, with optional CDS features to specify the coding segments. Our GFF reader will ignore redundant features like start_codon, stop_codon when whole CDS features were provided, or *UTR features when whole exon features were also given. However, it is still possible to provide only CDS and *UTR features and our GFF parser will reassemble the exonic structure accordingly.

Example of a GTF2 transcript record with minimal attributes:

<pre>
20	protein_coding	exon	9873504	9874841	.	+	.	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	CDS	9873504	9874841	.	+	0	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	exon	9877488	9877679	.	+	.	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	CDS	9877488	9877679	.	+	0	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	exon	9888412	9888586	.	+	.	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	CDS	9888412	9888586	.	+	0	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	exon	9891475	9891998	.	+	.	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
20	protein_coding	CDS	9891475	9891995	.	+	2	gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";
</pre>

## GFF3

As defined by the [GFF3](http://www.sequenceontology.org/gff3.shtml) specification, the parent features (usually transcripts, i.e. "mRNA" features) are required to have an **ID** attribute, but here again an optional gene_name attribute can be used to specify a common gene name abbreviation. If **gene_name** is not given, it can be also inferred from the Name or ID attributes of the parent gene feature of the current parent mRNA feature (if given in the input file). Exon or CDS features arerequired to have a Parent attribute whose value must match the value of the ID attribute of a parent transcript feature (usually a "mRNA" feature).

**Feature restrictions**

For various reasons we currently assume the following limits (maximum values) for the genomic length (span) of gene and transcript features:

- Genes and transcripts cannot span more than 7 Megabases on the genomic sequence
- Exons cannot be longer than 30 Kilobases
- Introns cannot be larger than 6 Megabases

Also, transcript IDs are expected to be unique per GFF input file (though we relaxed this restriction by limiting it to a chromosome/contig scope).

Due to these requirements, Cufflinks programs may fail to load the user provided GFF file, and an error message should specify the offending GFF record. The user is expected to remove or correct such GFF records in order to continue the analysis.

Example of a GFF3 transcript record with minimal attributes:

<pre>
ctg123	example	mRNA	1300	9950	.	+	.	ID=t_012143;gene_name=EDEN
ctg123	example	exon	1300	1500	.	+	.	Parent=t_012143
ctg123	example	exon	3000	3902	.	+	.	Parent=t_012143
ctg123	example	exon	5000	5500	.	+	.	Parent=t_012143
ctg123	example	exon	7000	9000	.	+	.	Parent=t_012143
ctg123	example	exon	9400	9950	.	+	.	Parent=t_012143
ctg123	example	CDS	3301	3902	.	+	0	Parent=t_012143
ctg123	example	CDS	5000	5500	.	+	1	Parent=t_012143
ctg123	example	CDS	7000	7600	.	+	1	Parent=t_012143
</pre>


## The gffread utility
A program called gffread is included with the Cufflinks package and it can be used to verify or perform various operations on GFF files (use gffread -h to see the various usage options). Because the program shares the same GFF parser code with Cufflinks and other programs in the Tuxedo suite, it could be used to verify that a GFF file from a certain annotation source is correctly "understood" by these programs. Thus the gffread utility can be used to simply read the transcripts from the file and print these transcripts back, in either GFF3 (default) or GTF2 format (-T option), discarding any extra attributes and keeping only the essential ones, so the user can quickly verify if the transcripts in that file are properly parsed by the GFF reader code. The command line for such a quick cleanup and visual inspection of a given GFF file could be:

{% highlight bash %}
gffread -E annotation.gff -o- | more
{% endhighlight %}


This will show the minimalist GFF3 re-formatting of the transcript records given in the input file (annotation.gff in this example). The -E option directs gffread to "expose" (display warnings about) any potential issues encountered while parsing the given GFF file.
In order to see the GTF2 version of the same transcripts the -T option should be added:

{% highlight bash %}
gffread -E annotation.gff -T -o- | more
{% endhighlight %}

From these examples it can be seen that gffread can also be used to convert a file between GTF2 and GFF3 formats. 


## Extracting transcript sequences

The gffread utility can be used to generate a FASTA file with the DNA sequences for all transcripts in a GFF file. For this operation a fasta file with the genomic sequences have to be provided as well. For example, one might want to extract the sequence of all transfrags assembled from a Cufflinks assembly session. This can be accomplished with a command line like this:

{% highlight bash %}
gffread -w transcripts.fa -g /path/to/genome.fa transcripts.gtf
{% endhighlight %}

The file genome.fa in this example would be a multi fasta file with the genomic sequences of the target genome. This also requires that every contig or chromosome name found in the 1st column of the input GFF file (transcript.gtf in this example) must have a corresponding sequence entry in chromosomes.fa. This should be the case in our example if genome.fa is the file corresponding to the same genome (index) that was used for mapping the reads with Tophat. Note that the retrieval of the transcript sequences this way is going to be quicker if a fasta index file (genome.fa.fai in this example) is found in the same directory with the genomic fasta file. Such an index file can be created with samtools prior to running gffread, like this: 

{% highlight bash %}
samtools faidx genome.fa
{% endhighlight %}

Then in subsequent runs using the -g option gffread will find the fasta index and use it to speed up the extraction of the transcript sequences. 

#Sample sheets for Cuffdiff and Cuffnorm

Both Cuffdiff and Cuffnorm can be run by specifying a list of SAM, BAM, or CXB files at the command line. For analysis with many such files, specifying them in this way can be cumbersome and error-prone. Both programs also allow you to specify these inputs in a simple, tab-delimited table. Create a file called sample_sheet.txt or another name of your choice, and specify samples as follows, one per line:

Column number|	Column name|	Example|	Description
-|-|-|-
1|	sample_id|	C1_R1.sam|	the path to the SAM/BAM/CXB file for this sample
2|	group_label| C1|	The condition label for this sample. Replicates of a condition should share the same label
{: class="table"}

To run Cuffdiff or Cuffnorm with a sample sheet, create one and then at the command line, provide the --use-sample-sheet option and replace the list of SAM/BAM/CXB files with the name of your sample sheet file, as follows:

**<tt>cuffdiff --use-sample-sheet \<transcripts.gtf\> \<sample_sheet.txt\> </tt>**

An example sample sheet might look like this:

<pre>
sample_id	group_label
C1_R1.sam	C1
C1_R2.sam	C1
C2_R1.sam	C2
C2_R2.sam	C2
</pre>

# Output formats used in the Cufflinks suite

##FPKM Tracking format

FPKM tracking files use a generic format to output estimated expression values. Each FPKM tracking file has the following format:

Column number|	Column name|	Example|	Description
-|-|-|-
1|	tracking_id|	TCONS_00000001|	A unique identifier describing the object (gene, transcript, CDS, primary transcript)
2|	class_code|	=|	The class_code attribute for the object, or "-" if not a transcript, or if class_code isn't present
3|	nearest_ref_id|	NM_008866.1|	The reference transcript to which the class code refers, if any
4|	gene_id|	NM_008866|	The gene_id(s) associated with the object
5|	gene_short_name|	Lypla1|	The gene_short_name(s) associated with the object
6|	tss_id|	TSS1|	The tss_id associated with the object, or "-" if not a transcript/primary transcript, or if tss_id isn't present
7|	locus|	chr1:4797771-4835363|	Genomic coordinates for easy browsing to the object
8|	length|	2447|	The number of base pairs in the transcript, or '-' if not a transcript/primary transcript
9|	coverage|	43.4279|	Estimate for the absolute depth of read coverage across the object
10|	q0_FPKM|	8.01089|	FPKM of the object in sample 0
11|	q0_FPKM_lo|	7.03583|	the lower bound of the 95% confidence interval on the FPKM of the object in sample 0
12|	q0_FPKM_hi|	8.98595|	the upper bound of the 95% confidence interval on the FPKM of the object in sample 0
13|	q0_status|	OK|	Quantification status for the object in sample 0. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
14|	q1_FPKM|	8.55155|	FPKM of the object in sample 1
15|	q1_FPKM_lo|	7.77692|	the lower bound of the 95% confidence interval on the FPKM of the object in sample 0
16|	q1_FPKM_hi|	9.32617|	the upper bound of the 95% confidence interval on the FPKM of the object in sample 1
17|	q1_status|	9.32617|	the upper bound of the 95% confidence interval on the FPKM of the object in sample 1. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
3N + 12|	qN_FPKM|	7.34115|	FPKM of the object in sample N
3N + 13|	qN_FPKM_lo|	6.33394|	the lower bound of the 95% confidence interval on the FPKM of the object in sample N
3N + 14|	qN_FPKM_hi|	8.34836|	the upper bound of the 95% confidence interval on the FPKM of the object in sample N
3N + 15|	qN_status|	OK|	Quantification status for the object in sample N. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
{: class="table"}

##Count Tracking format

Count tracking files use a generic format to output estimated fragment count values. Each Count tracking file has the following format:

Column number|	Column name|	Example|	Description
-|-|-|-
1|	tracking_id|	TCONS_00000001|	A unique identifier describing the object (gene, transcript, CDS, primary transcript)
2|	q0_count|	201.334|	Estimated (externally scaled) number of fragments generated by the object in sample 0
3|	q0_count_variance|	5988.24|	Estimated variance in the number of fragments generated by the object in sample 0
4|	q0_count_uncertainty_var|	170.21|	Estimated variance in the number of fragments generated by the object in sample 0 due to fragment assignment uncertainty.
5|	q0_count_dispersion_var|	4905.63|	Estimated variance in the number of fragments generated by the object in sample 0 due to cross-replicate variability.
6|	q0_status|	OK|	Quantification status for the object in sample 0. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
7|	q1_count|	201.334|	Estimated (externally scaled) number of fragments generated by the object in sample 1
8|	q1_count_variance|	5988.24|	Estimated variance in the number of fragments generated by the object in sample 1
9|	q1_count_uncertainty_var|	170.21|	Estimated variance in the number of fragments generated by the object in sample 1 due to fragment assignment uncertainty.
10|	q1_count_dispersion_var|	4905.63|	Estimated variance in the number of fragments generated by the object in sample 1 due to cross-replicate variability.
11|	q1_status|	OK|	Quantification status for the object in sample 1. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
7|	qN_count|	201.334|	Estimated (externally scaled) number of fragments generated by the object in sample N
4N + 6|	qN_count_variance|	7.34115|	Estimated variance in the number of fragments generated by the object in sample N
4N + 7|	qN_count_uncertainty_var|	6.33394|	Estimated variance in the number of fragments generated by the object in sample N due to fragment assignment uncertainty.
4N + 8|	qN_count_dispersion_var|	8.34836|	Estimated variance in the number of fragments generated by the object in sample N due to cross-replicate variability.
4N + 9|	qN_status|	OK|	Quantification status for the object in sample N. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
{: class="table"}

##Read Group Tracking format

Read group tracking files record per-replicate expression and count data. Each Count tracking file has the following format:

Column number|	Column name|	Example|	Description
-|-|-|-
1|	tracking_id|	TCONS_00000001|	A unique identifier describing the object (gene, transcript, CDS, primary transcript)
2|	condition|	Fibroblasts|	Name of the condition
3|	replicate|	1|	Name of the replicate of the condition
4|	raw_frags|	170.21|	The estimate number of (unscaled) fragments originating from the object in this replicate
5|	internal_scaled_frags|	4905.63|	Estimated number of fragments originating from the object, after transforming to the internal common count scale (for comparison between replicates of this condition.)
6|	external_scaled_frags|	99.21|	Estimated number of fragments originating from the object, after transforming to the external common count scale (for comparison between conditions)
7|	FPKM|	201.334|	FPKM of this object in this replicate
8|	effective_length|	5988.24|	The effective length of the object in this replicate. Currently a reserved, unreported field.
9|	status|	OK|	Quantification status for the object. Can be one of OK (deconvolution successful), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents deconvolution.
{: class="table"}

## Read group info - read_groups.info

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

## Simple-table expression format

Cuffnorm reports two different types of files with this format: *.fpkm_table files and *.count_table files for each group of features: genes, transcripts, TSS groups, and CDS groups. The files start with a column indicating the feature ID for each row. There is one subsequent column for each sample in the analysis:

Column number| Column name | Example | Description
---|---|---|---
1|	tracking_id|	TCONS_00000001|	A unique identifier describing the object (gene, transcript, CDS, primary transcript)
2|	conditionX_N|	=|	The FPKM value (for *.fpkm_table files) or normalized fragment count (for *.count_table files) for this feature in replicate N of conditionX
{: class="table"}

## Simple-table gene attributes format

Cuffdiff reports metadata for each gene, transcript, TSS group, and CDS group in the following tab delimited format:

Column number|	Column name|	Example|	Description
---|---|---|---
1|	tracking_id|	TCONS_00000001|	A unique identifier describing the object (gene, transcript, CDS, primary transcript)
2|	class_code|	=|	The class_code attribute for the object, or "-" if not a transcript, or if class_code isn't present
3|	nearest_ref_id|	NM_008866.1|	The reference transcript to which the class code refers, if any
4|	gene_id|	NM_008866|	The gene_id(s) associated with the object
5|	gene_short_name|	Lypla1|	The gene_short_name(s) associated with the object
6|	tss_id|	TSS1|	The tss_id associated with the object, or "-" if not a transcript/primary transcript, or if tss_id isn't present
7|	locus|	chr1:4797771-4835363|	Genomic coordinates for easy browsing to the object
8|	length|	2447|	The number of base pairs in the transcript, or '-' if not a transcript/primary transcript
{: class="table"}

## Simple-table sample attributes format

Cuffnorm reports some information about each sample (i.e. each SAM, BAM, or CXB file) in the analysis in the following format:

Column number|	Column name|	Example|	Description
---|---|---|---
1|	sample_id|	q1_0|	A unique identifier describing the sample. Has the format conditionX_N, meaning replicate N of conditionX.
2|	file|	C1_R1.sam|	The path to the file (SAM/BAM/CXB) attribute for the sample
3|	total_mass|	94610|	The total (un-normalized) number of fragment alignments for this sample
4|	internal_scale|	1.0571|	The scaling factor used to normalize this sample library size.
5|	external_scale|	1|	Reserved
{: class="table"}



