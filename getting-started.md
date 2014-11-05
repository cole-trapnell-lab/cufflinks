---
layout: page
permalink: /getting-started/
description: "Instructions on how to install Monocle."
modified: 2013-09-11
tags: [monocle, install, setup]
---

* Table of Contents
{:toc}


##Install quick-start

###Required software
Monocle runs in the [R statistical computing environment](http://www.r-project.org/). You will need R version 3.0 or higher. You will need to install the following packages through CRAN:

- VGAM
- irlba
- matrixStats
- igraph (version >= 0.7.0)
- combinat
- fastICA
- grid
- ggplot2
- reshape2
- plyr
- parallel
- methods

You can install these packages by starting an R session and typing: 

{% highlight R %}
> install.packages(c("VGAM", "irlba", "matrixStats", "igraph", 
"combinat", "fastICA", "grid", "ggplot2", 
"reshape2", "plyr", "parallel", "methods"))
{% endhighlight %}

You will also need to install [Bioconductor](http://bioconductor.org/install/): 

{% highlight R %}
> source("http://bioconductor.org/biocLite.R") 
> biocLite()
{% endhighlight %}

###Installing the myoblast example data
Monocle includes a detailed documentation vignette and snippets of example code that depend on the skeletal muscle myoblast data described in Trapnell and Cacchiarelli et al. You'll need to install the data package containing it before installing Monocle. To do so, type:

{% highlight bash %}
$ R CMD INSTALL HSMMSingleCell_0.99.0.tar.gz 
{% endhighlight %}

###Installing Monocle from source
To install the Monocle package, download the source tarball, change to the directory in which you saved it, and type

{% highlight bash %}
$ R CMD INSTALL monocle_0.99.0.tar.gz 
{% endhighlight %}

###Testing the installation
To ensure that Monocle was installed correctly, start a new R session and type:

{% highlight R %}
> library(monocle)
{% endhighlight %}

##Computing expression values for single cells

To use Monocle, you must first compute the expression of each gene in each cell for your experiment. There are a number of ways to do this for RNA-Seq. We recommend using Cufflinks, but you could also use [RSEM](http://deweylab.biostat.wisc.edu/rsem/), [eXpress](http://bio.math.berkeley.edu/eXpress/), Sailfish, or another tool for estimating gene and transcript expression levels from aligned reads. Here, we'll show a simplified workflow for using TopHat and Cufflinks to estimate expression. You can read more about how to use TopHat and Cufflinks to calculate expression [here](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html).

To estimate gene and transcript expression levels for single-cell RNA-Seq using [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) and [Cufflinks](http://cufflinks.cbcb.umd.edu/), you must have a file of RNA-Seq reads for each cell you captured. If you performed paired-end RNA-Seq, you should have two files for each cell. Depending on how the base calling was performed, the naming conventions for these files may differ. In the examples below, we assume that each file follows the format:

CELL_TXX_YYY.RZ.fastq.gz

Where XX is the time point at which the cell was collected in our experiment, YY is the well of the 96-well plate used during library prep, and Z is either 1 or 2 depending on whether we are looking at the left mate or the right mate in a paired end sequencing run. So CELL_T24_A01.R1.fastq.gz means we are looking at the left mate file for a cell collected 24 hours into our experiment and which was prepped in well A01 of the 24-hour capture plate. 

###Aligning reads to the genome with TopHat
We begin by aligning each cell's reads separately, so we will have one BAM file for each cell. The commands below show how to run each cell's reads through TopHat. These alignment commands can take a while, but they can be run in parallel if you have access to a compute cluster. If so, contact your cluster administrator for more information on how to run TopHat in a cluster environment. 

{% highlight bash %}
tophat -o CELL_T24_A01_thout -G GENCODE.gtf bowtie-hg19-idx CELL_T24_A01.R1.fastq.gz CELL_T24_A01.R2.fastq.gz 
tophat -o CELL_T24_A02_thout -G GENCODE.gtf bowtie-hg19-idx CELL_T24_A02.R1.fastq.gz CELL_T24_A02.R2.fastq.gz 
tophat -o CELL_T24_A03_thout -G GENCODE.gtf bowtie-hg19-idx CELL_T24_A03.R1.fastq.gz CELL_T24_A03.R2.fastq.gz 
{% endhighlight %}

The commands above show how to align the reads for each of three cells in the experiment. You will need to run a similar command for each cell you wish to include in your analysis. These TopHat alignment commands are simplified for brevity - there are options to control the number of CPUs used by TopHat and otherwise control how TopHat aligns reads that you may want to explore on the TopHat manual. The key components of the above commands are:

- The -o option, which sets the directory in which each cell's output will be written.
- The gene annotation file, specified with -G, which tells TopHat where to look for splice junctions.
- The Bowtie index for genome of your organism, in this case build hg19 of the human genome.
- The read files for each cell as mentioned above.

When the commands finish, there will be a BAM file in each cell's TopHat output directory. For example, CELL_T24_A01_thout/accepted_hits.bam will contain the alignments for cell T24_A01.

###Computing gene expression using Cufflinks
Now, we will use Cufflinks to estimate gene expression levels for each cell in your study. 

{% highlight bash %}
cuffquant -o CELL_T24_A01_cuffquant_out GENCODE.gtf CELL_T24_A01_thout/accepted_hits.bam 
cuffquant -o CELL_T24_A02_cuffquant_out GENCODE.gtf CELL_T24_A02_thout/accepted_hits.bam 
cuffquant -o CELL_T24_A03_cuffquant_out GENCODE.gtf CELL_T24_A03_thout/accepted_hits.bam 
{% endhighlight %}

The commands above show how to convert aligned reads for each cell into gene expression values for that cell. You will need to run a similar command for each cell you wish to include in your analysis. These commands are simplified for brevity - there are options to control the number of CPUs used by the cuffquant utility and otherwise control how cuffquant estimates expression that you may want to explore on the [Cufflinks](http://cufflinks.cbcb.umd.edu/) manual. The key components of the above commands are:

- The -o option, which sets the directory in which each cell's output will be written.
- The gene annotation file, which tells cuffquant what the gene structures are in the genome.
- The BAM file containing the aligned reads.

Next, you will need to merge the expression estimates into a single table for use with Monocle. You can do this with the following command: 

{% highlight bash %}
cuffnorm --use-sample-sheet -o sc_expr_out GENCODE.gtf sample_sheet.txt
{% endhighlight %}

The option --use-sample-sheet tells cuffnorm that it should look in the file sample_sheet.txt for the expression files, to make the above command simpler. If you choose not to use a sample sheet, you will need to specify the expression files on the command line directly. The sample sheet is a tab-delimited file that looks like this: 

| sample_name | group |
|:--------|:-------:|--------:|
| CELL_T24_A01_cuffquant_out/abundances.cxb   | T24_A01 |
| CELL_T24_A02_cuffquant_out/abundances.cxb   | T24_A02 |
| CELL_T24_A03_cuffquant_out/abundances.cxb   | T24_A03 |

Now, you are ready to load the expression data into Monocle and start analyzing your experiment. 

##Analyzing data with Monocle

Monocle provides a number of tools you can use to analyze your single cell expression experiments. To get started, we must create a CellDataSet object. You can do this with the commands below:

{% highlight R %}
> library(monocle)
> sample_sheet <- read.delim("sc_expr_out/samples.table", row.names=1)
> gene_annotations <- read.delim("sc_expr_out/genes.attr_table", row.names=1)
> fpkm_matrix <- read.delim("sc_expr_out/genes.fpkm_table", row.names=1)
> pd <- new("AnnotatedDataFrame", data = sample_sheet)
> my_data <- new("CellDataSet", exprs = as.matrix(fpkm_matrix), phenoData = pd, featureData = fd)
{% endhighlight %}

Now, you have created an object named "my_data" that stores your single-cell expression data. This object is the central object in Monocle. You will use it to identify differentially expressed genes and perform other analyses. To see what Monocle can do for you and how to proceed, please have a look at the [vignette (PDF)](http://www.bioconductor.org/packages/devel/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf). Good luck! 