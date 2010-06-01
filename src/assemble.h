#ifndef ASSEMBLE_H
#define ASSEMBLE_H

/*
 *  assemble.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <map>

//#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
//#include <lemon/concepts/bpugraph.h>
#include <boost/graph/adjacency_list.hpp>

#include "hits.h"
#include "bundles.h"
#include "scaffolds.h"


bool assemble_hits(BundleFactory& bundle_factory);

//bool intron_compatible(const MateHit& lhs, const MateHit& rhs);
bool read_hits_overlap(const ReadHit* lhs, const ReadHit* rhs);
bool read_hits_intron_agree(const ReadHit* h1, const ReadHit* h2);

int  match_length(const MateHit& m, int left, int right);

bool distance_compatible(const MateHit& lhs, 
						 const MateHit& rhs, 
						 int max_inner_dist);

typedef boost::adjacency_list<boost::vecS, 
							  boost::vecS, 
							  boost::bidirectionalS, 
							  boost::property<boost::vertex_name_t, Scaffold*> > DAG;

typedef boost::graph_traits<DAG>::vertex_descriptor DAGNode;

bool make_scaffolds(int bundle_left,
					int bundle_length,
					vector<Scaffold>& hits,
					vector<Scaffold>& scaffolds);

#endif
