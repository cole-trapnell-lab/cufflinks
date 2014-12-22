---
layout: page
permalink: /
tags: [cufflinks, RNA-Seq]
modified: 2013-09-13
---

Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one, taking into account biases in library preparation protocols. 

Cufflinks was originally developed as part of a collaborative effort between the [Laboratory for Mathematical and Computational Biology](http://bio.math.berkeley.edu/), led by Lior Pachter at UC Berkeley, Steven Salzberg's [computational genomics group](http://ccb.jhu.edu/) at the Institute of Genetic Medicine at Johns Hopkins University, and [Barbara Wold's lab](http://woldlab.caltech.edu/) at Caltech. The project is now maintained by [Cole Trapnell's](http://cole-trapnell-lab.github.io/) lab at the University of Washington.

Cufflinks is provided under the OSI-approved [Boost License](http://en.wikipedia.org/wiki/Boost_Software_License)

# News

*To get the latest updates on the Cufflinks project and the rest of the "Tuxedo tools", please subscribe to our [**mailing list**](https://lists.sourceforge.net/lists/listinfo/bowtie-bio-announce)* 

<ul class="post-list">
{% for post in site.posts %} 
  <li><article><a href="{{ site.url }}{{ post.url }}">{{ post.title }} <span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></a></article></li>
{% endfor %}
</ul>
