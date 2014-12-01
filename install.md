---
layout: page
permalink: /install/
description: "Documentation for the file formats used by the Cufflinks suite."
modified: 2013-09-11
tags: [cufflinks, cuffdiff, cuffquant, cuffnorm, formats, manual]
---

Cufflinks is available for Linux and Mac OS X. You can find the full list of releases below. 

The Cufflinks source code for each point release is available below as well. If you want to grab the current code, check out the <a href="https://github.com/cole-trapnell-lab/cufflinks">Cufflinks GitHub repository</a>.

<center>
<iframe src="http://ghbtns.com/github-btn.html?user=cole-trapnell-lab&repo=cufflinks&type=watch&count=true" allowtransparency="true" frameborder="0" scrolling="0" width="110" height="20"></iframe>

<iframe src="http://ghbtns.com/github-btn.html?user=cole-trapnell-lab&repo=cufflinks&type=fork&count=true" allowtransparency="true" frameborder="0" scrolling="0" width="110" height="20"></iframe>
</center>

<h1>Cufflinks Releases</h1>
<table class="table">
  <thead>
    <tr>
      <th style="text-align: left">Version</th>
      <th>Date</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>

{% for post in site.categories.releases %}
    <tr>
    	<td><a href="{{ site.url }}{{ post.url }}">{{ post.version }}</a></td>
    	<td><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></td>
    	<td><a href="{{ site.url }}/assets/downloads/cufflinks-{{ post.version }}.Linux_x86_64.tar.gz">Linux</a></td>
    	<td><a href="{{ site.url }}/assets/downloads/cufflinks-{{ post.version }}.OSX_x86_64.tar.gz">Mac OS X</a></td>
    	<td><a href="{{ site.url }}/assets/downloads/cufflinks-{{ post.version }}.tar.gz">Source</a></td>
    </tr>
{% endfor %} 
</table>