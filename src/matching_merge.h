#ifndef MATCHING_MERGE_H
#define MATCHING_MERGE_H

/*
 *  matching_merge.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>

//#include <lemon/topology.h>
#include <lemon/smart_graph.h>
#include <lemon/bipartite_matching.h>

#include "scaffold_graph.h"
#include "scaffolds.h"

using namespace std;

typedef lemon::SmartBpUGraph ReachGraph;

template<class Matcher>
void make_chains_from_matching(const ReachGraph& bp, 
							   const Matcher& matcher,
							   vector<vector<DAGNode> >& chains)
{
	// Make chains out of the matching
	ReachGraph::ANodeMap<ReachGraph::UEdge> matched_a_nodes(bp);
	matcher.aMatching(matched_a_nodes);
	
	ReachGraph::BNodeMap<ReachGraph::UEdge> matched_b_nodes(bp);
	matcher.bMatching(matched_b_nodes);
	
	set<ReachGraph::ANode> chain_heads;
	
	for (ReachGraph::ANodeIt i(bp); i!=lemon::INVALID; ++i)
	{
		int a_id = bp.aNodeId(i);
		ReachGraph::ANode a = bp.nodeFromANodeId(a_id);
		if (matched_a_nodes[a] == lemon::INVALID)
			chain_heads.insert(bp.nodeFromANodeId(bp.aNodeId(i)));
	}
	
	for (set<ReachGraph::ANode>::iterator i = chain_heads.begin();
		 i != chain_heads.end();
		 ++i)
	{
		vector<DAGNode> chain;
		ReachGraph::ANode n = *i;
		chain.push_back(bp.aNodeId(*i));
		while(true)
		{
			//int a_id = bp.aNodeId(n);
			int b_id = bp.bNodeId(n);
			
			ReachGraph::BNode b = bp.nodeFromBNodeId(b_id);
			
			if (matched_b_nodes[b] == lemon::INVALID)
			{
				break;
			}
			else
			{
				ReachGraph::ANode a_match_to_b = bp.source(matched_b_nodes[b]);
				chain.push_back(bp.aNodeId(a_match_to_b));
				n = a_match_to_b;
			}
		}
		chains.push_back(chain);
	}
	assert (chains.size() == chain_heads.size());
}

template<class Matcher>
void merge_from_matching(const ReachGraph& bp_graph,
                         const Matcher& matcher,
                         vector<Scaffold>& scaffolds)
{
    vector<Scaffold> merged_scaffolds;
    
    vector<vector<DAGNode> > chains;
    
    make_chains_from_matching(bp_graph, matcher, chains);
    for (size_t i = 0; i < chains.size(); ++i)
    {
        vector<Scaffold> chain;
        for (size_t j = 0; j < chains[i].size(); ++j)
        {
            chain.push_back(scaffolds[chains[i][j]]);
        }
        merged_scaffolds.push_back(Scaffold(chain));
    }
    
    sort (merged_scaffolds.begin(), merged_scaffolds.end(), scaff_lt);
    
    scaffolds = merged_scaffolds;   
}

void extend_chains_to_paths(const DAG& bundle_dag,
							vector<vector<DAGNode> >& chains,
							boost::adjacency_list<>& TC,
							DAGNode source,
							DAGNode sink,
							vector<vector<DAGNode> >& paths);

void make_scaffolds_from_paths(DAG& bundle_dag,
							   const vector<vector<DAGNode> >& paths,
							   vector<Scaffold>& scaffolds);
#endif
