#ifndef GRAPH_OPTIMIZE_H
#define GRAPH_OPTIMIZE_H
/*
 *  graph_optimize.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>

#include "bundles.h"
#include "scaffold_graph.h"
#include "scaffolds.h"

using namespace std;

using namespace boost;

template < typename PredecessorMap, 
           typename PathIDMap > 
class path_compress_visitor : public default_dfs_visitor 
{    
public:
    path_compress_visitor(PathIDMap pm) : curr_path_id(0), path_map(pm) {}
    
    template < typename Vertex, typename Graph >
    void initialize_vertex(Vertex u, const Graph & g) const
    {
        put(predecessor, u, u);
        put(path_map, u, u);
    }
    
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g)
    {
        //fprintf(stderr, "node %d has indegree %d, outdegree %d\n",u,in_degree(u, g),out_degree(u, g)); 
        if (in_degree(u, g) == 1)
        {
            
            Vertex v = get(predecessor, u);
            
            assert(v != u);
            
            if (out_degree(v, g) == 1)
            {
                // compress into predecessor's path 
                typename PathIDMap::value_type path = get(path_map, v);
                put(path_map, u, path);
                //fprintf(stderr, "\told path for node %d = %d\n", u, path);
                
                return;
            }
        }
        // start a new path
        curr_path_id++;
        put(path_map, u, curr_path_id);
        //fprintf(stderr, "\tnew path for node %d = %d\n", u, curr_path_id);
        
    }
    
    template < typename Edge, typename Graph >
    void tree_edge(Edge e, const Graph & g) const
    {
        put(predecessor, target(e, g), source(e, g));
    }
    
    size_t last_path_id() const { return curr_path_id; }
    
    PredecessorMap predecessor;
    
    size_t curr_path_id;
	PathIDMap path_map;
};

void fill_gaps(vector<Scaffold>& scaffolds, int fill_size);

void compress_fragments(vector<Scaffold>& hits);

bool collapse_equivalent_transfrags(vector<Scaffold>& scaffolds, 
                                    uint32_t max_rounds = 0xFFFFFFFF);

bool collapse_contained_transfrags(vector<Scaffold>& scaffolds, 
                                   uint32_t max_rounds = 0xFFFFFFFF);

void compress_overlap_dag_paths(DAG& bundle_dag,
                                vector<Scaffold>& hits);

#endif
