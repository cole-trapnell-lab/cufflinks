---
layout: page
permalink: /cuffmerge/
description: "Documentation for Cuffmerge."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

* Table of Contents
{:toc}

<h1>Merging assemblies with cuffmerge</h1>

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
