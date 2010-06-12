/*
 *  graph_optimize.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>

#include "graph_optimize.h"
// for graph optimization only

#include <lemon/topology.h>
#include <lemon/smart_graph.h>
#include <lemon/bipartite_matching.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "scaffold_graph.h"
#include "scaffolds.h"
#include "matching_merge.h"

using namespace std;
using namespace boost;

namespace ublas = boost::numeric::ublas;

void fill_gaps(vector<Scaffold>& scaffolds, int fill_size)
{
	for (size_t i = 0; i < scaffolds.size(); ++i)
		scaffolds[i].fill_gaps(fill_size);
}

enum ConflictState { UNKNOWN_CONFLICTS = 0, SAME_CONFLICTS, DIFF_CONFLICTS };

void same_conflicts(const vector<Scaffold>& scaffolds,
                    ublas::mapped_matrix<ConflictState>& conflict_states,
                    vector<vector<size_t> >* conf_out)
{       
    vector<vector<size_t> > conflicts(scaffolds.size());
    
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = i; j < scaffolds.size(); ++j)
		{
			if (i == j)
				continue;
			if (Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], 0))
			{
				if (!Scaffold::compatible(scaffolds[i], scaffolds[j]))
                {
                    conflicts[i].push_back(j);
                    conflicts[j].push_back(i);
                }
			}
		}
        
	}
    
    for (size_t i = 0; i < scaffolds.size(); ++i)
	{   
        sort(conflicts[i].begin(), conflicts[i].end());
    }
    
    for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = i; j < scaffolds.size(); ++j)
		{
			if (i == j)
				continue;
			if (Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], 0))
			{
				vector<size_t> diff;
                set_symmetric_difference(conflicts[i].begin(), 
                                         conflicts[i].end(),
                                         conflicts[j].begin(), 
                                         conflicts[j].end(),
                                         back_inserter(diff));
                if (!diff.empty())
                {
                    conflict_states(i,j) = DIFF_CONFLICTS;
                    conflict_states(j,i) = DIFF_CONFLICTS;
                }
                else
                {
                    conflict_states(i,j) = SAME_CONFLICTS;
                    conflict_states(j,i) = SAME_CONFLICTS;
                }
			}
		}
	}
    
    if (conf_out)
    {
        *conf_out = conflicts;
    }
}

void add_non_constitutive_to_scaffold_mask(const vector<Scaffold>& scaffolds,
										   vector<bool>& scaffold_mask)
{	
	//scaffold_mask = vector<bool>(scaffolds.size(), 0);
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
        for (size_t j = 0; j < scaffolds.size(); ++j)
        {
            if (Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], 0) &&
                !Scaffold::compatible(scaffolds[i], scaffolds[j]))
            {
                scaffold_mask[i] = true;
            }
        }
	}
}

bool collapse_contained_transfrags(vector<Scaffold>& scaffolds, 
                                   uint32_t max_rounds)
{
	// The containment graph is a bipartite graph with an edge (u,v) when
	// u is (not necessarily properly) contained in v and is the two are
	// compatible.
	typedef lemon::SmartBpUGraph ContainmentGraph;
	normal norm(0, 0.1);
	bool performed_collapse = false;
    
	while (max_rounds--)
	{
		
#if ASM_VERBOSE
		fprintf(stderr, "\tStarting new collapse round\n");
#endif
        
        ublas::mapped_matrix<ConflictState> conflict_states;
        conflict_states = ublas::mapped_matrix<ConflictState>(scaffolds.size(), 
                                                              scaffolds.size(), 
                                                              scaffolds.size() * scaffolds.size());
        
        vector<vector<size_t> > conf_out;
        same_conflicts(scaffolds, conflict_states, &conf_out);
		
		ContainmentGraph containment;
		
		
		typedef pair<ContainmentGraph::ANode, ContainmentGraph::BNode> NodePair;
		vector<NodePair> node_ids;
		vector<size_t> A_to_scaff(scaffolds.size());
		vector<size_t> B_to_scaff(scaffolds.size());
		
		for (size_t n = 0; n < scaffolds.size(); ++n)
		{
			NodePair p = make_pair(containment.addANode(),
								   containment.addBNode());
			node_ids.push_back(p);
			A_to_scaff[containment.aNodeId(p.first)] = n;
			B_to_scaff[containment.bNodeId(p.second)] = n;
		}
		
		bool will_perform_collapse = false;
		for (size_t i = 0; i < scaffolds.size(); ++i)
		{
			for (size_t j = 0; j < scaffolds.size(); ++j)
			{
				if (i == j)
					continue;
                
				if (scaffolds[i].contains(scaffolds[j]) &&
					Scaffold::compatible(scaffolds[i], scaffolds[j]))
				{
                    bool same = (conflict_states(i,j) == SAME_CONFLICTS);
                    if (!same)
                    {
                        //bool X = Scaffold::compatible(scaffolds[0], scaffolds[4]);
                        continue;
                    }   
					// To gaurd against the identity collapse, which won't 
					// necessary reduce the total number of scaffolds.
					if (scaffolds[j].contains(scaffolds[i]) && i < j)
						continue;
					const NodePair& nj = node_ids[j];
					const NodePair& ni = node_ids[i];
					assert (nj.first != ni.second);
                    
					will_perform_collapse = true;
					ContainmentGraph::UEdge e = containment.addEdge(nj.first,
																	ni.second);						
                    
				}
			}
		}
		
		if (will_perform_collapse == false)
			return performed_collapse;
		
		lemon::MaxBipartiteMatching<ContainmentGraph>  matcher(containment);
        
#if ASM_VERBOSE
		fprintf(stderr, "\tContainment graph has %d nodes, %d edges\n", containment.aNodeNum(), containment.uEdgeNum());
		fprintf(stderr, "\tFinding a maximum matching to collapse scaffolds\n");
#endif
        
		matcher.run();
        
#if ASM_VERBOSE
		fprintf(stderr, "\t\tWill collapse %d scaffolds\n", matcher.matchingSize());
#endif
		
		ContainmentGraph::UEdgeMap<bool> matched_edges(containment);
		
		matcher.matching(matched_edges);
		
        merge_from_matching(containment, matcher, scaffolds);
        
		performed_collapse = true;
	}
	return performed_collapse;
}

void compress_consitutive(vector<Scaffold>& hits)
{
    vector<bool> scaffold_mask;
    
#ifdef ASM_VERBOSE
    fprintf(stderr, "Building constitutivity mask..."); 
#endif
    
    scaffold_mask = vector<bool>(hits.size(), false);
    add_non_constitutive_to_scaffold_mask(hits, scaffold_mask);
    
#ifdef ASM_VERBOSE
    fprintf(stderr, "done\n"); 
#endif
    
    vector<Scaffold> constitutive;
    vector<Scaffold> non_constitutive;
    
    for (size_t i = 0; i < scaffold_mask.size(); ++i)
    {
        if (!scaffold_mask[i])
            constitutive.push_back(hits[i]);
        else
            non_constitutive.push_back(hits[i]);
    }
    
#ifdef ASM_VERBOSE
    size_t pre_compress = hits.size();
#endif
    hits.clear();
    if (!constitutive.empty())
    {
        Scaffold compressed = Scaffold(constitutive); 
        vector<Scaffold> completes; 
        compressed.get_complete_subscaffolds(completes);
        hits.insert(hits.end(), completes.begin(), completes.end()); 
    }
    
    hits.insert(hits.end(), non_constitutive.begin(), non_constitutive.end());
    sort(hits.begin(), hits.end(), scaff_lt);
#ifdef ASM_VERBOSE
    size_t post_compress = hits.size();
    size_t delta = pre_compress - post_compress;
    double collapse_ratio = delta / (double) pre_compress; 
    fprintf(stderr, 
            "Compressed %lu of %lu constitutive fragments (%lf percent)\n", 
            delta, 
            pre_compress, 
            collapse_ratio);
#endif
}

void compress_fragments(vector<Scaffold>& fragments)
{
    compress_consitutive(fragments);
    
#if ASM_VERBOSE
    fprintf(stderr,"\tPerforming preliminary containment collapse on %lu fragments\n", fragments.size());
    size_t pre_hit_collapse_size = fragments.size();
#endif
    if (perform_full_collapse)
    {
        collapse_contained_transfrags(fragments);
    }
    
#if ASM_VERBOSE
    size_t post_hit_collapse_size = fragments.size();
    fprintf(stderr,"\tIgnoring %lu strictly contained fragments\n", pre_hit_collapse_size - post_hit_collapse_size);
#endif   
}

void compress_overlap_dag_paths(DAG& bundle_dag,
                                vector<Scaffold>& hits)
{
    HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);
    
    vector_property_map<size_t> path_for_scaff;
    path_compress_visitor<vector_property_map<DAGNode>,
    vector_property_map<size_t> > vis(path_for_scaff);
    depth_first_search(bundle_dag,
                       visitor(vis));
    
    vector<vector<Scaffold> > compressed_paths(hits.size()+1);

    vector<Scaffold> new_scaffs;
    
    for (size_t i = 0; i < num_vertices(bundle_dag); ++i)
    {
        size_t path_id = path_for_scaff[i];
        assert (path_id < compressed_paths.size());
        const Scaffold* h = hits_for_node[i];
        if (h)
        {
            compressed_paths[path_id].push_back(*h);
        }
    }
    for (size_t i = 0; i < compressed_paths.size(); ++i)
    {
        if (!compressed_paths[i].empty())
        {
            //fprintf(stderr, "Path %d has %d fragments in it\n", i, compressed_paths[i].size());
            new_scaffs.push_back(Scaffold(compressed_paths[i]));
        }
    }
    //hits = new_scaffs;
    fprintf(stderr, "Compressed overlap graph from %d to %d fragments (%f percent)\n", 
            hits.size(), 
            new_scaffs.size(), 
            (hits.size() - new_scaffs.size())/(double)hits.size());
    hits = new_scaffs;
    sort(hits.begin(), hits.end(), scaff_lt);
    create_overlap_dag(hits, bundle_dag);
}

