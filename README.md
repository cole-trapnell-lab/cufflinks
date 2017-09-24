
# Cufflinks

The main *website* for cufflinks is [here](http://cole-trapnell-lab.github.io/cufflinks/)

*NOTE*: If you're looking for old releases of Cufflinks, including source, you can find them [here](http://cole-trapnell-lab.github.io/cufflinks/install/).

Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one, taking into account biases in library preparation protocols. 

Cufflinks was originally developed as part of a collaborative effort between the [Laboratory for Mathematical and Computational Biology](http://bio.math.berkeley.edu/), led by Lior Pachter at UC Berkeley, Steven Salzberg's [computational genomics group](http://ccb.jhu.edu/people/salzberg/) at the Institute of Genetic Medicine at Johns Hopkins University, and [Barbara Wold's lab](http://woldlab.caltech.edu/) at Caltech. The project is now maintained by [Cole Trapnell's](http://cole-trapnell-lab.github.io/) lab at the University of Washington.

Cufflinks is provided under the OSI-approved [Boost License](http://en.wikipedia.org/wiki/Boost_Software_License)

# News

*To get the latest updates on the Cufflinks project and the rest of the "Tuxedo tools", please subscribe to our [**mailing list**](https://lists.sourceforge.net/lists/listinfo/bowtie-bio-announce)* 

# Install quick-start

## Installing a pre-compiled binary release

In order to make it easy to install Cufflinks, we provide a few binary packages [here](http://cole-trapnell-lab.github.io/cufflinks/install/) to save users from the occasionally frustrating process of building Cufflinks, which requires that you install the Boost libraries. To use the binary packages, simply download the appropriate one for your machine, untar it, and make sure the cufflinks,cuffdiff and cuffcompare binaries are in a directory in your PATH environment variable.

# Building Cufflinks from source

In order to build Cufflinks, you must have the [Boost C++ libraries](http://www.boost.org/) (version 1.47 or higher) installed on your system. See below for instructions on installing Boost.

## Installing Boost

1. Download Boost and the bjam build engine. **WARNING:** Due to a serious issue with Boost Serlialization library introduced in version 1.56, Cufflinks currently can only be built with Boost version 1.55 or lower.  The issue is expected to be fixed in the upcoming Boost v1.59.
2. Unpack bjam and add it to your PATH.
3. Unpack the Boost tarball and cd to the Boost source directory. This directory is called the BOOST_ROOT in some Boost installation instructions.
4. Build Boost. Note that you can specify where to put Boost with the --prefix option. The default Boost installation directory is /usr/local. Take note of the boost installation directory, because you will need to tell the Cufflinks installer where to find Boost later on.

- If you are on Mac OS X, type (all on one line): 
```bash
bjam --prefix=<YOUR_BOOST_INSTALL_DIRECTORY> --toolset=darwin architecture=x86 address_model=32_64 link=static runtime-link=static --layout=versioned stage install
```

- If you are on a 32-bit Linux system, type (all on one line): 
```bash
bjam --prefix=<YOUR_BOOST_INSTALL_DIRECTORY> --toolset=gcc architecture=x86 address_model=32 link=static runtime-link=static stage install
```

- If you are on a 64-bit Linux system, type (all on one line): 
```bash
bjam --prefix=<YOUR_BOOST_INSTALL_DIRECTORY> --toolset=gcc architecture=x86 address_model=64 link=static runtime-link=static stage install
```

## Installing the SAM tools

1. [Download the SAM tools](http://samtools.sourceforge.net/)
2. Unpack the SAM tools tarball and cd to the SAM tools source directory.
3. Build the SAM tools by typing make at the command line.
4. Choose a directory into which you wish to copy the SAM tools binary, the included library <tt>libbam.a</tt>, and the library headers. A common choice is <tt>/usr/local/</tt>.
5. Copy libbam.a to the lib/ directory in the folder you've chosen above (e.g. <tt>/usr/local/lib/</tt>)
6. Create a directory called "bam" in the <tt>include/</tt> directory (e.g. <tt>/usr/local/include/bam</tt>)
7. Copy the headers (files ending in <tt>.h</tt>) to the include/bam directory you've created above (e.g. <tt>/usr/local/include/</tt>bam)
8. Copy the samtools binary to some directory in your <tt>PATH</tt>.

## Installing the Eigen libraries

1. [Download Eigen](http://eigen.tuxfamily.org/)
2. Unpack the Eigen tarball and cd to the Eigen source directory.
3. Copy the Eigen/ subdirectory someplace on your system where you keep header files (e.g. /usr/local/include)

## Building Cufflinks

### If you are starting from a source tarball downloaded from [here](http://cole-trapnell-lab.github.io/cufflinks/install/):

Unpack the Cufflinks source tarball (in this example for version 2.2.1):
```bash
tar zxvf cufflinks-2.2.1.tar.gz
```
Change to the Cufflinks directory:
```bash
cd cufflinks-2.2.1
```

### If you want to clone the Cufflinks github repo:
```bash
git clone https://github.com/cole-trapnell-lab/cufflinks.git
cd cufflinks
autoreconf --install
```
The above will generate the configure script.

### To configure Cufflinks prior to the build

If Boost is installed somewhere other than /usr/local, you will need to tell the installer where to find it using the --with-boost option. Specify where to install Cufflinks using the --prefix option.
```bash
./configure --prefix=/path/to/cufflinks/install --with-boost=/path/to/boost --with-eigen=/path/to/eigen
```

If you see any errors during configuration, verify that you are using Boost version 1.47 or higher, and that the directory you specified via --with-boost contains the boost header files and libraries. See the Boost Getting started page for more details. If you copied the SAM tools binaries to someplace other than /usr/local/, you may need to supply the --with-bam configuration option.
Finally, make and install Cufflinks.
```bash
make
make install
```

## Testing the installation

1. [Download](http://cufflinks.cbcb.umd.edu/downloads/test_data.sam) the test data
2. In the directory where you placed the test file, type:

```bash
cufflinks ./test_data.sam
```

You should see the following output:

<pre>
[bam_header_read] EOF marker is absent. The input is probably truncated.
[bam_header_read] invalid BAM binary header (this is not a BAM file).
File ./test_data.sam doesn't appear to be a valid BAM file, trying SAM...
[13:23:15] Inspecting reads and determining fragment length distribution.
> Processed 1 loci.                            [*************************] 100%
Warning: Using default Gaussian distribution due to insufficient paired-end reads in open ranges.  
It is recommended that correct paramaters (--frag-len-mean and --frag-len-std-dev) be provided.
> Map Properties:
>       Total Map Mass: 102.50
>       Read Type: 75bp x 75bp
>       Fragment Length Distribution: Truncated Gaussian (default)
>                     Estimated Mean: 200
>                  Estimated Std Dev: 80
[13:23:15] Assembling transcripts and estimating abundances.
> Processed 1 loci.                            [*************************] 100%
</pre>

Verify that the file transcripts.gtf is in the current directory and looks like this (your file will have GTF attributes, omitted here for clarity)

<pre>
test_chromosome Cufflinks       exon    53      250     1000    +       . 
test_chromosome Cufflinks       exon    351     400     1000    +       . 
test_chromosome Cufflinks       exon    501     550     1000    +       .
</pre>	

# Common uses of the Cufflinks package

Cufflinks includes a number of tools for analyzing RNA-Seq experiments. Some of these tools can be run on their own, while others are pieces of a larger workflow. The complexity of your workflow depends on what you want to achieve with your analysis. For a complete discussion of how Cufflinks can help you with your analysis, please [see our protocol paper](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html). The paper includes a diagram (Figure 2) describing how the various parts of the Cufflinks package (and its companion tool TopHat) fit together. As of version 2.2.0, you can also run Cuffquant and Cuffnorm to make large scale analyses easier to handle. The figure below is an updated version of Figure 2 showing how the two utilities released after the protocol paper appeared fit into the workflow: 

<div style="text-align:center">
![Workflow]({{ site.url }}/images/tuxedo_workflow.png)
</div>

You can use Cuffquant to pre-compute gene expression levels for each of your samples, which can save time if you have to re-run part of your analysis. Using Cuffquant also makes it easier to spread the load of computation for lots of samples across multiple computers. If you don't want to perform differential expression analysis, you can run Cuffnorm instead of Cuffdiff. Cuffnorm produces simple tables of expression values that you can look at in R (for example) to cluster samples and perform other follow up analysis.	

# Using pre-built annotation packages

A number of steps in the Tuxedo package work better if you have pre-existing gene annotations. How you can use these annotations is detailed in our [protocol paper](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html). Illumina has kindly provided a large number of annotation packages for commonly used model organisms and humans. You can find these packages [here]({{ site.url }}/igenome_table/index.html).	

# References

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
