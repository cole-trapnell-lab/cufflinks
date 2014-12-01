---
layout: page
permalink: /benchmarks/
description: "Benchmarks of Cuffdiff accuracy."
modified: 2014-11-19
tags: [cufflinks, cuffdiff, performance, benchmarks]
---


# Cuffdiff performance assessments

Cuffdiff is a robust tool for differential analysis of RNA-Seq experiments. It performs as well as or better than other commonly used tools over a large range of RNA-Seq experimental designs. Here, we provide up-to-date performance data to show that the current release performs at or exceeds the standards described in the Cuffdiff 2 manuscript. The experiments here include both simulated and real data.

**PLEASE NOTE:** The data presented here are **not intended to provide a comprehensive evaluation of Cuffdiff's performance**. These experiments are part of a testing suite we use to make sure that changes made between versions work as expected, and do not hurt performance or cause other problems. While the simulation experiments reflect a wide range of experimental scenarios that we commonly encounter when performing RNA-Seq, there may be use cases not captured here.

# Real-world datasets

In Trapnell and Hendrickson et al, we analyzed several real RNA-Seq datasets and compared them with matched data from other expression measurement technologies. These comparisons showed that Cuffdiff 2 is robust and accurate. The studies below include a re-analysis of those datasets.

- [HOXA1 knockdown in fibroblasts]({{ site.url }}/benchmarks/hoxkd/index.html)
- [Griffith et al]({{ site.url }}/benchmarks/hoxkd/index.html)
- [MAQC]({{ site.url }}/benchmarks/maqc/index.html)

# Simulation studies

The simulations studies presented below assess Cuffdiff's performance under various hypothetical experimental scenarios. They are all intended to be “realistic”, in the sense that users may encounter these experimental settings during routine RNA-Seq analysis. We have performed extensive sequencing using our in-house RNA-Seq read simulator, TuxSim.

Most of the simulations are designed around the fibroblast control data presented in Trapnell and Hendrickson et al. That is, the number of replicates, sequencing depth, cross-replicate variability, and other properties have be set to match the real experiment in the paper. Except where otherwise noted, the experiments are designed as follows:

Experimental parameter|	Value
---|---
Genome|	Homo sapiens (hg19)
Transcriptome|	UCSC protein coding genes
Experimental base data set|	Trapnell and Hendrickson et al Lung fibroblast
Sequencing depth|	~15,000,000 fragments per replicate
Replicates|	3 per condition
Read length|	100bp
Paired-end reads|	Yes
Isoform-switching|	Minimal
Isoform perturbed|	Most abundant of each gene
Number of genes perturbed|	1000
Minimum depth for gene perturbation|	10 fragments
{: class="table"}

You can find the results for each simulation in the sections listed. The “basic” simulation covers what we feel will be a representative experiment for many users performing in-vitro studies (e.g. knockdowns, overexpression, drug treatment in a dish, etc). It includes a full report of Cuffdiff's accuracy. Some simulations do not show results to the same level of detail as the “basic” simulation, because the results are largely the same as that setting. We thus exclude some panels for brevity.

- [Basic performance]({{ site.url }}/benchmarks/basic/index.html)
	- Estimation of expression values in a single condition
	- Estimation of expression values between conditions
	- True and false positives
- [Varying replication]({{ site.url }}/benchmarks/replication_series/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Varying depth]({{ site.url }}/benchmarks/depth_curve/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Uneven depth]({{ site.url }}/benchmarks/uneven_libs/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Uneven replication]({{ site.url }}/benchmarks/uneven_replication/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Biased DE genes]({{ site.url }}/benchmarks/biased_de_genes/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Many DE genes]({{ site.url }}/benchmarks/many_de_genes/index.html)
	- Estimation of expression values between conditions
	- True and false positives
- [Increased isoform switching]({{ site.url }}/benchmarks/multi_iso/index.html)
	- Estimation of expression values in a single condition
	- Estimation of expression values between conditions
	- True and false positives
- [High cross-replicate variability]({{ site.url }}/benchmarks/high_variability/index.html)
	- Estimation of expression values in a single condition
	- Estimation of expression values between conditions
	- True and false positives