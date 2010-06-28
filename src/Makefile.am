AM_CPPFLAGS = -I$(top_srcdir)/src 

EXTRA_DIST = $(top_srcdir)/AUTHORS $(top_srcdir)/make_bin.sh

AM_CXXFLAGS =

bin_PROGRAMS = \
	cufflinks \
	cuffcompare \
	cuffdiff 
	#	gtf_reads

noinst_HEADERS = \
	codons.h \
	clustering.h \
	differential.h \
	GArgs.h \
	GBase.h \
	gdna.h \
	GFaSeqGet.h \
	GFastaFile.h \
	gff.h \
	GHash.hh \
	GList.hh \
	GStr.h \
	assemble.h \
	bundles.h \
	filters.h \
	hits.h \
	tokenize.h \
	scaffolds.h \
	abundances.h \
	genes.h \
	transitive_closure.h \
    transitive_reduction.h \
	gtf_tracking.h \
	bundles.h \
    scaffold_graph.h \
    matching_merge.h \
    graph_optimize.h \
	lemon/maps.h \
	lemon/bin_heap.h \
	lemon/fib_heap.h \
	lemon/bipartite_matching.h \
	lemon/topology.h \
	lemon/bucket_heap.h \
	lemon/smart_graph.h \
	lemon/list_graph.h \
	lemon/error.h \
	lemon/graph_utils.h \
	lemon/bfs.h \
	lemon/dfs.h \
	lemon/math.h \
	lemon/tolerance.h \
	lemon/graph_adaptor.h \
	lemon/concept_check.h \
	lemon/bits/default_map.h \
	lemon/bits/array_map.h \
	lemon/bits/vector_map.h \
	lemon/bits/debug_map.h \
	lemon/bits/map_extender.h \
	lemon/bits/base_extender.h \
	lemon/bits/graph_adaptor_extender.h \
	lemon/bits/graph_extender.h \
	lemon/bits/invalid.h \
	lemon/bits/path_dump.h \
	lemon/bits/traits.h \
	lemon/bits/utility.h \
	lemon/bits/alteration_notifier.h \
	lemon/bits/variant.h \
	lemon/concepts/graph.h \
	lemon/concepts/bpugraph.h \
	lemon/concepts/graph_components.h \
	lemon/concepts/ugraph.h \
	lemon/concepts/heap.h \
	lemon/concepts/maps.h \
	lemon/concepts/path.h \
	lemon/concepts/matrix_maps.h \
	common.h 

noinst_LIBRARIES = libcufflinks.a libgc.a

libcufflinks_a_SOURCES = \
	clustering.cpp \
	differential.cpp \
	common.cpp \
	assemble.cpp \
	tokenize.cpp \
	abundances.cpp \
	scaffolds.cpp \
	hits.cpp \
	genes.cpp \
	bundles.cpp \
	filters.cpp \
    scaffold_graph.cpp \
    matching_merge.cpp \
    graph_optimize.cpp 
	
libgc_a_SOURCES = \
	codons.cpp \
	GArgs.cpp \
	GBase.cpp \
	gdna.cpp \
	GStr.cpp \
	GFaSeqGet.cpp \
	gff.cpp \
	gtf_tracking.cpp 

cufflinks_SOURCES = cufflinks.cpp  
cufflinks_LDADD = libcufflinks.a libgc.a $(BOOST_THREAD_LIB) $(BAM_LIB) 
cufflinks_LDFLAGS =  $(BOOST_LDFLAGS) $(BAM_LDFLAGS) #$(ZLIB_LDFLAGS)

cuffcompare_SOURCES = cuffcompare.cpp
cuffcompare_LDADD = libgc.a

cuffdiff_SOURCES = cuffdiff.cpp
cuffdiff_LDADD = libcufflinks.a libgc.a $(BOOST_THREAD_LIB) $(BAM_LIB)
cuffdiff_LDFLAGS =  $(BOOST_LDFLAGS) $(BAM_LDFLAGS) $(ZLIB_LDFLAGS)

#gtf_reads_SOURCES = gtf_reads.cpp
#gtf_reads_LDADD = libcufflinks.a libgc.a $(BOOST_THREAD_LIB) $(BAM_LIB) 
#gtf_reads_LDFLAGS =  $(BOOST_LDFLAGS) $(BAM_LDFLAGS) #$(ZLIB_LDFLAGS)
