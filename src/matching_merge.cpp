/*
 *  matching_merge.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "matching_merge.h"

using namespace std;
using namespace boost;

void find_path(const DAG& bundle_dag,
			   const adjacency_list<>& TC,
			   const DAGNode& source,
			   const DAGNode& target,
			   vector<DAGNode>& path)
{
	if (source == target)
		return;
	bool done = false;
	DAGNode curr = source;
	while(!done)
	{
		graph_traits<DAG>::adjacency_iterator i, iend;
		for (tie(i,iend) = adjacent_vertices(curr, bundle_dag); i != iend; ++i)
		{
			DAGNode I = *i;
			pair<adjacency_list<>::edge_descriptor, bool> p;
			p = edge(I, target, TC);
			if (p.second)
			{
				path.push_back(*i);
				curr = *i;
				break;
			}
			if (*i == target)
			{
				path.push_back(*i);
				done = true;
				break;
			}
		}
	}
}

void extend_chains_to_paths(const DAG& bundle_dag,
							vector<vector<DAGNode> >& chains,
							adjacency_list<>& TC,
							DAGNode source,
							DAGNode sink,
							vector<vector<DAGNode> >& paths)
{	
	//Extend each chain to a path
	for(size_t c = 0; c < chains.size(); ++c)
	{
		vector<DAGNode>& chain = chains[c];
		assert (!chain.empty());
		reverse(chain.begin(), chain.end());
		vector<DAGNode> path;
		find_path(bundle_dag, TC, source, chain[0], path);
		for (size_t n = 1; n < chain.size(); ++n)
		{
			assert (path.back() == chain[n - 1]);
			DAGNode last = chain[n-1];
			DAGNode next = chain[n];
			find_path(bundle_dag, TC, last, next, path);
		}
		find_path(bundle_dag, TC, chain.back(), sink, path);
		assert (path.back() == sink);
		path.pop_back();
		paths.push_back(path);
	}
}

void make_scaffolds_from_paths(DAG& bundle_dag,
							   const vector<vector<DAGNode> >& paths,
							   vector<Scaffold>& scaffolds)
{
	HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);
	for (size_t p = 0; p < paths.size(); ++p)
	{
		const vector<DAGNode>& path = paths[p];
		
		vector<Scaffold> path_alignments;
		for (size_t m = 0; m < path.size(); ++m)
		{
			//fprintf(stderr, "%d ", scaff_id);
			path_alignments.push_back(*(hits_for_node[path[m]]));
		}
		//fprintf(stderr,"\n");
		//fprintf(stderr, "\tMerging path %d into scaffold\n", p);
		Scaffold s(path_alignments);
		
		//fprintf(stderr, "PATH %d\n-------------------------------------\n",s.get_id());
		scaffolds.push_back(s);
	}
	sort(scaffolds.begin(), scaffolds.end(), scaff_lt);
}


