#ifndef SCAFFOLD_GRAPH_H
#define SCAFFOLD_GRAPH_H
/*
 *  scaffold_graph.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/2/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/version.hpp>

#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include <boost/property_map/vector_property_map.hpp>
#endif

class Scaffold;

typedef boost::adjacency_list<boost::vecS, 
                              boost::vecS, 
                              boost::bidirectionalS, 
                              boost::property<boost::vertex_name_t, Scaffold*> > DAG;

typedef boost::graph_traits<DAG>::vertex_descriptor DAGNode;

typedef boost::property_map<DAG, boost::vertex_name_t>::type HitsForNodeMap;

bool create_overlap_dag(std::vector<Scaffold>& hits, DAG& bundle_dag);
std::pair<DAGNode, DAGNode> add_terminal_nodes(DAG& bundle_dag);

#endif
