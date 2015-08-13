# CUFFLINKS

Cufflinks is a reference-guided assembler for RNA-Seq experiments. It
simultaneously assembles transcripts from reads and estimates their relative
abundances, without using a reference annotation.  The software expects as 
input RNA-Seq read alignments in SAM format (http://samtools.sourceforge.net).

Here's an example spliced read alignment record:
<pre>
s6.25mer.txt-913508	16	chr1	4482736	255	14M431N11M	*	0	0	CAAGATGCTAGGCAAGTCTTGGAAG	IIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	XS:A:-
</pre>

This record includes a custom tag used by Cufflinks to determine the strand 
of the mRNA from which this read originated.  Often, RNA-Seq experiments lose
strand information, and so reads will map to either strand of the genome.  
However, strandedness of spliced alignments can often be inferred from the 
orientation of consensus splice sites, and Cufflinks requires that spliced 
aligments have the custom strand tag XS, which has SAM attribute type "A", and
can take values of "+" or "-".  If your RNA-Seq protocol is strand specific,
including this tag for all alignments, including unspliced alignments, will 
improve the assembly quality.

The SAM records MUST BE SORTED by reference coordinate, like so:

```bash
sort -k 3,3 -k 4,4n hits.sam 
```

The program is fully threaded, and when running with multiple threads, should
be run on a machine with plenty of RAM. 4 GB per thread is probably reasonable
for most experiments.  Since many experiments feature a handful of genes that
are very abundantly transcribed, Cufflinks will spend much of its time 
assembling" a few genes.  When using more than one thread, Cufflinks may 
appear to "hang" while these genes are being assembled.

Cufflinks assumes that library fragment lengths are size selected and normally
distributed. When using paired end RNA-Seq reads, you must take care to supply
Cufflinks with the mean and variance on the inner distances between mate 
pairs. For the moment, Cufflinks doesn't support assembling mixtures of paired
end reads from different fragment size distributions.  Mixtures of single 
ended reads (of varying lengths) with paired ends are supported.

Cufflinks also assumes that the donor does not contain major structural 
variations with respect to the reference.  Tumor genomes are often highly 
rearranged, and while Cufflinks may eventually identify gene fusions and 
gracefully handle genome breakpoints, users are encouraged to be careful when
using Cufflinks on tumor RNA-Seq experiments. 

The full manual may be found at http://cole-trapnell-lab.github.io/cufflinks/

# CUFFCOMPARE

Please see http://cole-trapnell-lab.github.io/cufflinks/manual


# REQUIREMENTS


Cufflinks is a standalone tool that requires gcc 4.0 or greater, and runs on
Linux and OS X.  It depends on Boost (http://www.boost.org) version 1.38 or 
higher.

# REFERENCES

Cufflinks is an ongoing research project as well as a suite of tools. Here are the papers that describe the science behind the programs. If you use Cufflinks, please cite these papers in your work!

**Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation**
Cole Trapnell, Brian Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Jeltje van Baren, Steven Salzberg, Barbara Wold, Lior Pachter.
*Nature Biotechnology*, 2010, doi:10.1038/nbt.1621
*Note: This is the original Cufflinks paper. Please cite this paper if you use Cufflinks in your work.*

**Improving RNA-Seq expression estimates by correcting for fragment bias**
Adam Roberts, Cole Trapnell, Julie Donaghey, John L. Rinn, Lior Pachter.
*Genome Biology*, 2011, doi:10.1186/gb-2011-12-3-r22
*Note: This paper describes improvements made to Cufflinks to handle bias in RNA-Seq read coverage. Please cite this paper if you use Cufflinks with the <tt>-b</tt> option in your work.*

**Identification of novel transcripts in annotated genomes using RNA-Seq**
Adam Roberts, Harold Pimentel, Cole Trapnell, Lior Pachter.
*Bioinformatics*, 2011, doi:10.1093/bioinformatics/btr355
*Note: This paper describes the RABT assembly algorithm. Please cite this paper if you use Cufflinks in RABT mode in your work.*

**Differential analysis of gene regulation at transcript resolution with RNA-seq**
Cole Trapnell, David Hendrickson, Martin Sauvageau, Loyal Goff, John L. Rinn, Lior Pachter
*Nature Biotechnology*, 2012, doi:10.1038/nbt.2450
*Note: This paper describes Cuffdiff 2. Please cite this paper if you use Cuffdiff in your work.*

Cufflinks builds on many ideas, including some
proposed in the following papers:

1. Ali Mortazavi, Brian A Williams, Kenneth McCue, Lorian Schaeffer and Barbara 
Wold, "Mapping and quantifying mammalian transcriptomes by RNA-Seq",Nature 
Methods, volume 5, 621 - 628 (2008)
2. Hui Jiang and Wing Hung Wong, "Statistical Inferences for isoform expression", 
Bioinformatics, 2009 25(8):1026-1032=
3.Nicholas Eriksson, Lior Pachter, Yumi Mitsuya, Soo-Yon Rhee, Chunlin Wang, 
Baback Gharizadeh, Mostafa Ronaghi, Robert W. Shafer, Niko Beerenwinkel, "Viral 
population estimation using pyrosequencing", PLoS Computational Biology, 
4(5):e1000074

