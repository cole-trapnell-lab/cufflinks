/*
 *  assemble.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>
#include <list>
#include <vector>
//#include <deque>
#include <iostream>
//#include <numeric>

#include <lemon/topology.h>

//#include <lemon/graph_utils.h>
#include <lemon/bipartite_matching.h>

//#include <boost/thread.hpp>
#include <boost/ref.hpp>
// DON'T move this, or mystery compiler errors will result. Affects gcc >= 4.1
#include <boost/graph/vector_as_graph.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <boost/vector_property_map.hpp>

#include <boost/math/distributions/normal.hpp> // for normal_distribution
using boost::math::normal; // typedef provides default type is double.

#include "transitive_closure.h"
#include "transitive_reduction.h"

#include "common.h"
#include "assemble.h"
//#include "abundances.h"
#include "bundles.h"
//#include "filters.h"
//#include "genes.h"
#include "clustering.h"

// for graph optimization only
#include <boost/numeric/ublas/matrix_sparse.hpp>


using namespace boost;
using namespace std;

bool mate_graphs(const HitBundle& bundle, BundleStats* stats);

typedef lemon::SmartBpUGraph ReachGraph;
void reachability_graph(const DAG& dag, ReachGraph& reach_graph);

#define ASM_VERBOSE 0

struct HitBufBasket
{
	HitBufBasket(int coord, Scaffold* h, DAGNode d)
		: expiration_coord(coord), hit(h), node(d) {}
	int expiration_coord;
	Scaffold* hit;
	DAGNode node;
};

bool right_lt (const HitBufBasket& lhs, 
			   const HitBufBasket& rhs)
{
	return lhs.expiration_coord < rhs.expiration_coord;
}

struct Expired
{
	Expired(int ec) : expiration_coord(ec) {}
	bool operator()(const HitBufBasket& lhs)
	{
		return lhs.expiration_coord <= expiration_coord;
	}
	
	int expiration_coord;
};

enum ConnectState { UNKNOWN, CONNECT, DONT_CONNECT };

template <class CompatibilityMap, class ConnectMap,  class Tag>
struct connect_visitor
	: public base_visitor<connect_visitor<CompatibilityMap, ConnectMap,  Tag> >
{
    typedef Tag event_filter;
    connect_visitor(CompatibilityMap compatibility, ConnectMap connect, DAGNode t) 
		: _connect(connect),_compatibility(compatibility), _target(t) { }
	
    template <class Vertex, class Graph>
    void operator()(Vertex u, const Graph& g)
	{
		typedef graph_traits<Graph> GraphTraits;
		
		typename GraphTraits::adjacency_iterator v, vend;
		
		if (_compatibility[u] == true)
		{
			for (tie(v,vend) = adjacent_vertices(u, g); v != vend; ++v)
			{
				if (_compatibility[*v])
				{
					//fprintf(stderr, "Avoiding a redundant edge from %d to %d\n", u, *v);
					_connect[u] = DONT_CONNECT;
					return;
				}
			}
			
			// If we get here, u is compatible with the target, but has no
			// compatible successors, so it's safe to add the edge after the DFS
			_connect[u] = CONNECT;
		}
		else
		{
			_connect[u] = DONT_CONNECT;
		}
		//put(_compat, v, compatible);
    }
	
    ConnectMap _connect;
	CompatibilityMap _compatibility;
	DAGNode _target;
};

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
	//Extend each chain to paths
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

template <class CompatibilityMap, class ConnectMap, class Tag>
connect_visitor<CompatibilityMap, ConnectMap, Tag>
record_connections(CompatibilityMap compatibility, 
				   ConnectMap connect, 
				   DAGNode target, 
				   Tag) 
{
    return connect_visitor<CompatibilityMap, ConnectMap, Tag> (compatibility, connect, target);
}


//template <class PredecessorMap, class Tag>
//struct predecessor_recorder
//: public base_visitor<predecessor_recorder<PredecessorMap, Tag> >
//{
//    typedef Tag event_filter;
//    predecessor_recorder(PredecessorMap pa) : m_predecessor(pa) { }
//    template <class Edge, class Graph>
//    void operator()(Edge e, const Graph& g) {
//        put(m_predecessor, target(e, g), source(e, g));
//    }
//    PredecessorMap m_predecessor;
//};

using namespace boost;
template < typename PredecessorMap, 
           typename PathIDMap > 
class path_compress_visitor : public default_dfs_visitor {
    //typedef typename property_traits < PredecessorMap >::value_type T;
public:
    path_compress_visitor(PathIDMap pm) : curr_path_id(0), path_map(pm) {}
    
    template < typename Vertex, typename Graph >
    void initialize_vertex(Vertex u, const Graph & g) const
    {
        //put(m_dtimemap, u, m_time++);
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
        //put(m_dtimemap, u, m_time++);
        put(predecessor, target(e, g), source(e, g));
    }
    
    size_t last_path_id() const { return curr_path_id; }
    
    //TimeMap m_dtimemap;
    //TimeMap m_ftimemap;
    
    PredecessorMap predecessor;
    PathIDMap path_map;
    
    size_t curr_path_id;
    //T & m_time;
};

//template <class PredecessorMap, class Tag>
//predecessor_recorder<PredecessorMap, Tag>
//record_predecessors(PredecessorMap pa, Tag) {
//    return predecessor_recorder<PredecessorMap, Tag> (pa);
//}


long long weight_of_merge(Scaffold& lhs, 
						  Scaffold& rhs, 
						  double source_psi,
						  double target_psi)
{
	//double expected_cov_diff = max(1.0, 0.1 * (source_doc + target_doc));
	//normal cov_norm(0, expected_cov_diff);
	
	normal cov_test_norm(0, 1.0);

	double score = 0.0;
	
	// HACK: This early breakout prevents spliced reads that cross exactly one
	// intron from being matched up if they both cross the same intron.
	// Otherwise, we get phasing problems, as introns aren't matched up long 
	// distance.  Ugh..
	if (Scaffold::overlap_in_genome(lhs, rhs, 0))
	{
		bool lh_intron = lhs.has_intron();
		bool rh_intron = rhs.has_intron();
		
		if (lh_intron && rh_intron)
		{
			vector<pair<int, int> > lh_gaps = lhs.gaps();
			vector<pair<int, int> > rh_gaps = rhs.gaps();
			if (lh_gaps.size() == 1 && lh_gaps == rh_gaps)
				return 999999999;
		}

		//return 999999999;
	}
//	if (source_psi < 1 || target_psi < 1)
//	{
//		int a =4;
//	}
	if (source_psi > 0 && target_psi > 0 )
	{
		double test_stat = log(1.0 - abs(source_psi - target_psi));
		score = test_stat;
		assert (score <= 0.0);
	}
	else
	{
		return 999999999;
	}
	
	assert (score <= 0.0);
	if (score >= -1e-6)
		score = -1e-6;
	int weight = (int)(score * -1e6);
	
	//assert (a_scaff != b_scaff);
//#if ASM_VERBOSE
//	if (lh_intron && rh_intron)
//		fprintf(stderr, "Score between reads (%lf, %lf: %lf : %lf) = %d, distance = %d\n", source_doc, target_doc, cov_stat, cov_score, weight, rhs.left() - lhs.right());
//#endif
	return weight;
	
}

enum ConflictState { UNKNOWN_CONFLICTS = 0, SAME_CONFLICTS, DIFF_CONFLICTS };
void same_conflicts(const vector<Scaffold>& scaffolds,
                    ublas::mapped_matrix<ConflictState>& conflict_states)
{       
    vector<vector<size_t> > conflicts(scaffolds.size());
    
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = i; j < scaffolds.size(); ++j)
		{
			if (i == j)
				continue;
			if (Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], olap_radius))
			{
				if (!Scaffold::compatible(scaffolds[i], scaffolds[j], true))
                {
                    conflicts[i].push_back(j);
                    conflicts[j].push_back(i);
                }
			}
		}
	}
    
    for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = i; j < scaffolds.size(); ++j)
		{
			if (i == j)
				continue;
			if (Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], olap_radius))
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
    
}

void strict_containment_collapse(vector<Scaffold>& scaffolds)
{
    
    ublas::mapped_matrix<ConflictState> conflict_states;
    conflict_states = ublas::mapped_matrix<ConflictState>(scaffolds.size(), 
                                                          scaffolds.size(), 
                                                          scaffolds.size() * scaffolds.size());
    same_conflicts(scaffolds, conflict_states);
    
	vector<char> contained_by_somebody(scaffolds.size(),0);
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = 0; j < scaffolds.size(); ++j)
		{
			if (i == j)
				continue;
			if (scaffolds[i].strictly_contains(scaffolds[j]))
			{
				if (scaffolds[j].strictly_contains(scaffolds[i]) && i < j)
					continue; // shouldn't happen
                ConflictState c = conflict_states(i,j);
                //fprintf(stderr, "%d,%d = %d\n", i, j, c); 
                if (conflict_states(i,j) == SAME_CONFLICTS)
                    contained_by_somebody[j] = true;
			}
		}
	}
	
	vector<Scaffold> merged_scaffolds;
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		if (!contained_by_somebody[i])
		{
			merged_scaffolds.push_back(scaffolds[i]);
		}
	}
	
	scaffolds = merged_scaffolds;
}

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

bool collapse_contained_scaffolds(vector<Scaffold>& scaffolds, uint32_t max_rounds = 0xFFFFFFFF)
{
	// The containment graph is a bipartite graph with an edge (u,v) when
	// u is (not necessarily properly) contained in v and is the two are
	// compatible.
	typedef lemon::SmartBpUGraph ContainmentGraph;
	normal norm(0, 0.1);
	bool performed_collapse = false;
	//sort (scaffolds.begin(), scaffolds.end(), scaff_lt);
	
    size_t original_size = scaffolds.size();
    
	while (max_rounds--)
	{
		
#if ASM_VERBOSE
		fprintf(stderr, "\tStarting new collapse round\n");
#endif
				
		// Remove identical or strictly contained elements
		
		//strict_containment_collapse(scaffolds);
        
        ublas::mapped_matrix<ConflictState> conflict_states;
        conflict_states = ublas::mapped_matrix<ConflictState>(scaffolds.size(), 
                                                              scaffolds.size(), 
                                                              scaffolds.size() * scaffolds.size());
        same_conflicts(scaffolds, conflict_states);
        
//		size_t prelim_collapse_size = scaffolds.size();
//#if ASM_VERBOSE
//		fprintf (stderr, "\tPreliminary collapse had %lu scaffolds, left %lu scaffolds\n", original_size, prelim_collapse_size);
//#endif
//		if (original_size != prelim_collapse_size)
//		{
//			performed_collapse = true;
//		}
		
		ContainmentGraph containment;
		//ReachGraph::UEdgeMap<long long> weights(containment);
		
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
					Scaffold::compatible(scaffolds[i], scaffolds[j], true) &&
                    conflict_states(i,j) == SAME_CONFLICTS)
				{
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
		
		lemon::MaxBipartiteMatching<ReachGraph>  matcher(containment);

#if ASM_VERBOSE
		
		fprintf(stderr, "\tContainment graph has %d nodes, %d edges\n", containment.aNodeNum(), containment.uEdgeNum());
		fprintf(stderr, "\tFinding a maximum matching to collapse scaffolds\n");
#endif
		matcher.run();
#if ASM_VERBOSE
		fprintf(stderr, "\t\tWill collapse %d scaffolds\n", matcher.matchingSize());
#endif
		
		ContainmentGraph::UEdgeMap<bool> matched_edges(containment);
		vector<char> keep_scaffold(scaffolds.size(), 1);
		matcher.matching(matched_edges);
		
		vector<Scaffold> merged_scaffolds;
		
		vector<vector<DAGNode> > chains;
		
		make_chains_from_matching(containment,matcher, chains);
		for (size_t i = 0; i < chains.size(); ++i)
		{
			vector<Scaffold> chain;
			for (size_t j = 0; j < chains[i].size(); ++j)
			{
				chain.push_back(scaffolds[chains[i][j]]);
			}
			merged_scaffolds.push_back(Scaffold(chain));
		}

#if ASM_VERBOSE
		fprintf (stderr, "\t\tpre = %d, post = %d\n", (int)scaffolds.size(), (int)merged_scaffolds.size());
#endif
		sort (merged_scaffolds.begin(), merged_scaffolds.end(), scaff_lt);
		scaffolds = merged_scaffolds;
		performed_collapse = true;
	}
	return performed_collapse;
}

typedef property_map<DAG, vertex_name_t>::type HitsForNodeMap;

bool create_dag_for_bundle(vector<Scaffold>& hits,
						   DAG& bundle_dag)
{
    bundle_dag = DAG();
	vector<Scaffold>::iterator hi =  hits.begin();
	bool found_compatible_scaffolds = false;
	
	typedef list<HitBufBasket> HitBuf;
	HitBuf hit_buf;
	
	HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);
	
	while (hi != hits.end())
	{
		int new_left = hi->left();
		int new_right = hi->right();

        //fprintf(stderr, "Adding to hit buffer: [%d, %d)\n", new_left, new_right);
        
		HitBufBasket new_basket(new_right, &(*hi), add_vertex(bundle_dag));
		hits_for_node[new_basket.node] = new_basket.hit;
		
		HitBuf::iterator new_end = remove_if(hit_buf.begin(), 
											 hit_buf.end(), 
											 Expired(new_left));
		
		hit_buf.erase(new_end, hit_buf.end());

		// Now check the each hit in the buffer for compatibility with this
		// new one
		
		vector<const Scaffold*> containing_hits;
		
		boost::vector_property_map<bool> c(num_vertices(bundle_dag));
		boost::vector_property_map<ConnectState> connected(num_vertices(bundle_dag));

		for (HitBuf::iterator bi = hit_buf.begin();
			 bi != hit_buf.end();
			 ++bi)
			
		{
			const Scaffold& lhs = *(bi->hit);
			const Scaffold& rhs = *(new_basket.hit);

			assert (lhs.left() <= rhs.left());
			if (!lhs.contains(rhs))
			{
                //fprintf(stderr, "Checking [%d, %d) and [%d, %d)\n", lhs.left(), lhs.right(), rhs.left(), rhs.right());
				if (Scaffold::compatible(lhs, rhs, true))
				{
					c[bi->node] = true;
				}
			}
		}
		
		for (HitBuf::iterator bi = hit_buf.begin();
			 bi != hit_buf.end();
			 ++bi)
			
		{
			if (connected[bi->node] == UNKNOWN)
			{
				depth_first_search(bundle_dag,
								   root_vertex(bi->node).
								   visitor(make_dfs_visitor(make_pair(record_connections(c, connected, new_basket.node, on_finish_vertex()), null_visitor()))));
			}
		}
		
		for (HitBuf::iterator bi = hit_buf.begin();
			 bi != hit_buf.end();
			 ++bi)
		{
			if (connected[bi->node] == CONNECT)
			{ 
				add_edge(bi->node, new_basket.node, bundle_dag);
				found_compatible_scaffolds = true;
			}
		}
		
		hit_buf.push_back(new_basket);
		
		++hi;
	}
	
	vector<bool> has_parent(num_vertices(bundle_dag), false);
	vector<bool> has_child (num_vertices(bundle_dag), false);
	
	graph_traits < DAG >::vertex_iterator u, uend;
	for (tie(u, uend) = vertices(bundle_dag); u != uend; ++u)
	{
		graph_traits < DAG >::adjacency_iterator v, vend;
		for (tie(v,vend) = adjacent_vertices(*u, bundle_dag); v != vend; ++v)
		{
			DAGNode U = *u;
			DAGNode V = *v;
			has_parent[V] = true;
			has_child[U] = true;
		}
	}
	
#ifdef DEBUG
	set<const Scaffold*> introns;
#endif
	for (size_t i = 0; i < num_vertices(bundle_dag); ++i)
	{
		if (has_child[i])
			continue;
		const Scaffold* hit_i = hits_for_node[i];

		for (size_t j = 0; j < num_vertices(bundle_dag); ++j)
		{
			if (has_parent[j])
				continue;
			const Scaffold* hit_j = hits_for_node[j];
			if (hit_i->right() < hit_j->left() &&
				hit_j->left() - hit_i->right() < olap_radius)
			{
				add_edge(i, j, bundle_dag);
			}
		}
	}
//    
#if DEBUG
    DAG tr;
    boost::vector_property_map<DAGNode> G_to_TR;
    property_map<DAG, vertex_index_t>::type w = get(vertex_index, bundle_dag);
    transitive_reduction(bundle_dag, 
                         tr, 
                         G_to_TR,
                         w);
    fprintf(stderr, "dag has %d edges, tr has %d edges\n", num_edges(bundle_dag), num_edges(tr));
    assert (num_edges(bundle_dag) == num_edges(tr));
#endif
	
//    bundle_dag = tr;
	return found_compatible_scaffolds;
}

inline bool partial_scaffold(const Scaffold& scaff)
{
	return scaff.has_unknown();
}


typedef map<DAGNode, pair<int, int> > DagToBp;
void reachability_graph(DAG& dag, 
						DAGNode src,
						DAGNode sink,
						ReachGraph& reach_graph,
						vector<ReachGraph::BNode> b_to_a,
						DagToBp& dag_to_bp,
						const adjacency_list<>& TC,
						const vector<bool>& scaffold_mask)
{	
	//typedef graph_traits<DAG>::vertex_descriptor Vertex;
	HitsForNodeMap hits_for_node = get(vertex_name, dag);
	
	//fprintf (stdout, "\tclosure edges:\t\t\t\%d\n", num_edges(TC));
	graph_traits < adjacency_list<> >::vertex_iterator v, vend;
	
	b_to_a.resize(num_vertices(TC));
	
	for (tie(v, vend) = vertices(TC); v != vend; ++v)
	{
		DagToBp::iterator itr = dag_to_bp.find(*v);
		if (itr == dag_to_bp.end())
		{
			ReachGraph::ANode A = reach_graph.addANode();
			int a = reach_graph.aNodeId(A);
			
			ReachGraph::BNode B = reach_graph.addBNode();
			int b = reach_graph.bNodeId(B);
			b_to_a[b] = A;
			dag_to_bp[*v] = make_pair(a, b);
			
		}
	}
	
	//size_t node_size = num_vertices(TC) * 
	//	(sizeof(ReachGraph::ANode) + sizeof(ReachGraph::BNode));
	
	//size_t edge_size = num_edges(TC) * sizeof(ReachGraph::UEdge);
			
	
	reach_graph.reserveEdge(num_edges(TC));
	reach_graph.reserveANode(num_vertices(TC));
	reach_graph.reserveBNode(num_vertices(TC));
	
	graph_traits < adjacency_list<> >::edge_iterator i, end;
	for (tie(i, end) = edges(TC); i != end; ++i)
	{
		int a_id = -1;
		int b_id = -1;
		DAGNode s = source(*i, TC);
		DAGNode t = target(*i, TC);

		DagToBp::iterator itr = dag_to_bp.find(s);
		if (itr == dag_to_bp.end())
		{
			assert (false);
		}
		else
		{
			a_id = itr->second.first;
		}
		
		itr = dag_to_bp.find(t);
		
		if (itr == dag_to_bp.end())
		{
			assert(false);
		}
		else
		{
			b_id = itr->second.second;
		}
		
		if (s == src || t == sink)
			continue;
		
		assert (a_id != -1);
		assert (b_id != -1);
		
		ReachGraph::ANode a_node = reach_graph.nodeFromANodeId(a_id);
		ReachGraph::BNode b_node = reach_graph.nodeFromBNodeId(b_id);
		ReachGraph::ANode a_for_b = b_to_a[b_id];
		
		assert (a_for_b != a_node);
		
		if (scaffold_mask[a_id] && 
			scaffold_mask[b_id])
		{
			reach_graph.addEdge(a_node, 
								b_node);
		}
	}
}

void add_non_constitutive_to_scaffold_mask(const vector<Scaffold>& scaffolds,
										   vector<bool>& scaffold_mask)
{	
	//scaffold_mask = vector<bool>(scaffolds.size(), 0);
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		//if (!intron_scaffold_mask[i])
		{
			for (size_t j = 0; j < scaffolds.size(); ++j)
			{
				if (/*intron_scaffold_mask[j] && */
					Scaffold::overlap_in_genome(scaffolds[i], scaffolds[j], 0) &&
					!Scaffold::compatible(scaffolds[i], scaffolds[j], true))
				{
					scaffold_mask[i] = true;
				}
			}
		}
	}
}

void add_weights_to_reachability_graph(ReachGraph& bp,
									   const HitsForNodeMap& hits_for_node,
									   const vector<Scaffold>& hits,
									   const vector<Scaffold>& scaffolds,
									   ReachGraph::UEdgeMap<long long>& weights)
{
	// number of reads of reads spliced in to each scaffold
	vector<double> spliced_in(scaffolds.size(), 0.0);
	
	// number spliced in / length of scaffold
	vector<double> density(scaffolds.size(), 0.0);
	
	// number of reads spliced in to scaffold + those that just overlap it
	//vector<double> overlapping(scaffolds.size(), 0.0);
	
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		for (size_t j = 0; j < hits.size(); ++j)
		{
			if (scaffolds[i].contains(hits[j]))
			{
				if (Scaffold::compatible(scaffolds[i],hits[j],true))
				{
					spliced_in[i]++;
				}
				//overlapping[i]++;
			}
		}
		density[i] = spliced_in[i] / scaffolds[i].length();
	}
	
	// percent spliced in = density / (total density of overlapping scaffolds) 
	vector<double> psi(scaffolds.size(), 0.0);
	vector<double> local_density(scaffolds.size(), 0.0);
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		double total_density = 0.0;
		int num_overlaps = 0;
		double compatible_density = 0.0;
		for (size_t j = 0; j < scaffolds.size(); ++j)
		{
			if (Scaffold::overlap_in_genome(scaffolds[i],scaffolds[j], 0))
			{
				total_density += density[j];
				num_overlaps++;
				if (Scaffold::compatible(scaffolds[i], scaffolds[j], true))
				{
					compatible_density += density[j];
				}
			}
		}
		if (total_density)
			psi[i] = compatible_density / total_density;
		local_density[i] = compatible_density;
	}
	
//	for (size_t i = 0; i < scaffolds.size(); ++i)
//	{
//		fprintf(stderr, "%d - %d: count = %lf, density = %lf, local = %lf psi = %lf\n", 
//				scaffolds[i].left(),
//				scaffolds[i].right(),
//				spliced_in[i],
//				density[i],
//				local_density[i],
//				psi[i]);
//	}
	
	for (ReachGraph::UEdgeIt i(bp); i!=lemon::INVALID; ++i)
	{
		ReachGraph::ANode a = bp.source(i);
		ReachGraph::BNode b = bp.target(i);
		DAGNode a_dag = bp.aNodeId(a);
		int a_id_for_b = bp.aNodeId(b);
		ReachGraph::ANode a_for_b = bp.nodeFromANodeId(a_id_for_b);
		assert (a_for_b != lemon::INVALID);
		DAGNode b_dag = a_id_for_b;
		
		Scaffold* a_scaff = hits_for_node[a_dag];
		Scaffold* b_scaff = hits_for_node[b_dag];

		size_t aidx = a_scaff - &scaffolds[0];
		size_t bidx = b_scaff - &scaffolds[0];
		
		double a_psi = psi[aidx];
		
		double b_psi = psi[bidx];
		
		
		//assert (a_psi != 1.0);
		//assert (b_psi != 1.0);
		
		long long weight = weight_of_merge(*a_scaff, 
										   *b_scaff, 
										   a_psi,
										   b_psi);
		
		if (weight < 0)
			weight = 10000000;
        
		//fprintf(stderr, "Similarity between %d, %d = %.20lf\n", bp.aNodeId(a), a_id_for_b, weight);
		weights[i] = weight;
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

struct ScaffBelowThresh
{
	ScaffBelowThresh(double thresh) : _thresh(thresh) {}
	bool operator()(const Scaffold& s)
	{
		return s.score()/(double)s.mate_hits().size() < _thresh;
	}
	
	double _thresh;
};

struct WorstMateScoreBelowThresh
{
	WorstMateScoreBelowThresh(double thresh) : _thresh(thresh) {}
	bool operator()(const Scaffold& s)
	{
		return s.worst_mate_score() < _thresh;
	}
	
	double _thresh;
};

struct ScaffHasIntron
{
	ScaffHasIntron(HitsForNodeMap m)
	:_m(m) {}
	bool operator()(const DAGNode& n)
	{
		const Scaffold* s = _m[n];
		if (!s) 
			return false;
		return s->has_intron();
	}
	
	HitsForNodeMap _m;
};

struct ScaffBelowIntronDoc
{
	ScaffBelowIntronDoc(const map<pair<int, int>, int>& docs,
						double doc_thresh)
	: intron_doc(docs), bundle_doc_thresh(doc_thresh) {}
	
	const map<pair<int,int>, int >& intron_doc;
	double bundle_doc_thresh;
	
	bool operator()(const Scaffold& scaff)
	{
		const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
		for (size_t i = 0; i < ops.size(); ++i)
		{
			if (ops[i].opcode == CUFF_INTRON)
			{
				map<pair<int, int>, int>::const_iterator itr;
				itr = intron_doc.find(make_pair(ops[i].g_left(), ops[i].g_right()));
				assert (itr != intron_doc.end());
				double doc = itr->second;
				//fprintf(stderr, "doc = %f\n", doc);
				if (doc < bundle_doc_thresh)
					return true;
			}
		}
		return false;
	}
};	

void fill_gaps(vector<Scaffold>& scaffolds, int fill_size)
{
	for (size_t i = 0; i < scaffolds.size(); ++i)
		scaffolds[i].fill_gaps(fill_size);
}

void holdout_transitivity_hazards(vector<Scaffold>& hits, 
								  vector<Scaffold>& hazards)
{
	vector<pair<int,int> > introns;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const vector<AugmentedCuffOp>& ops = hits[i].augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			const AugmentedCuffOp& op = ops[j];
			if (op.opcode == CUFF_INTRON)
				introns.push_back(make_pair(op.g_left(), op.g_right()));
		}
	}

	sort(introns.begin(), introns.end());
	vector<pair<int, int> >::iterator new_end = unique(introns.begin(), 
													   introns.end());
	introns.erase(new_end, introns.end());
	
	vector<pair<int, int> > evil_introns;
	for (size_t i = 0; i < introns.size(); ++i)
	{
		for (size_t j = i + 1; j < introns.size(); ++j)
		{
			if (overlap_in_genome(introns[i].first, introns[i].second,
								  introns[j].first, introns[j].second))
			{
				{
					evil_introns.push_back(introns[i]);
					evil_introns.push_back(introns[j]);
				}
			}
		}
	}
	
	sort(evil_introns.begin(), evil_introns.end());
	new_end = unique(evil_introns.begin(), evil_introns.end());
	evil_introns.erase(new_end, evil_introns.end());
	
	vector<Scaffold> filtered_hits;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		bool overlaps_evil_intron = false;
		const vector<AugmentedCuffOp>& ops = hits[i].augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			const AugmentedCuffOp& op = ops[j];
			if (op.opcode == CUFF_UNKNOWN)
			{
				for (size_t k = 0; k < evil_introns.size(); ++k)
				{
					
					if (overlap_in_genome(op.g_left(), op.g_right(),
										  evil_introns[k].first, evil_introns[k].second))
					{
						overlaps_evil_intron = true;
//						int left_diff = abs(op.g_left() - evil_introns[j].first);
//						int right_diff = abs(op.g_right() - evil_introns[j].second);
//						if (left_diff + right_diff > max_inner_dist)
//						{
//							overlaps_evil_intron = true;
//							break;
//						}
					}
				}				  
			}
		}
		if (overlaps_evil_intron)
		{
			hazards.push_back(hits[i]);
		}	
		else
		{
			filtered_hits.push_back(hits[i]);
		}
	}
	
#if ASM_VERBOSE
	fprintf(stderr, "Held out %lu scaffolds as transitivity hazards\n", hazards.size());
#endif
	hits = filtered_hits;
}

bool collapse_redundant_scaffolds(vector<Scaffold>& scaffolds, 
                                  uint32_t max_rounds = 0xFFFFFFFF, 
                                  uint32_t max_chain_size = 2)
{
	// The containment graph is a bipartite graph with an edge (u,v) when
	// u is (not necessarily properly) contained in v and is the two are
	// compatible.
	typedef lemon::SmartBpUGraph ContainmentGraph;
	
	bool performed_collapse = false;
	//sort (scaffolds.begin(), scaffolds.end(), scaff_lt);
	
	while (max_rounds--)
	{
		
#if ASM_VERBOSE
		fprintf(stderr, "\tStarting new collapse round\n");
#endif
        
		// Remove identical or strictly contained elements
		size_t original_size = scaffolds.size();
		//strict_containment_collapse(scaffolds);
		size_t prelim_collapse_size = scaffolds.size();
#if ASM_VERBOSE
		fprintf (stderr, "\tPreliminary collapse had %lu scaffolds, left %lu scaffolds\n", original_size, prelim_collapse_size);
#endif
		if (original_size != prelim_collapse_size)
		{
			performed_collapse = true;
		}
		
		ContainmentGraph containment;
		//ReachGraph::UEdgeMap<long long> weights(containment);
		
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
					Scaffold::compatible(scaffolds[i], scaffolds[j], true))
				{
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
		
		lemon::MaxBipartiteMatching<ReachGraph>  matcher(containment);
        
#if ASM_VERBOSE
		
		fprintf(stderr, "\tContainment graph has %d nodes, %d edges\n", containment.aNodeNum(), containment.uEdgeNum());
		fprintf(stderr, "\tFinding a maximum matching to collapse scaffolds\n");
#endif
		matcher.run();
#if ASM_VERBOSE
		fprintf(stderr, "\t\tWill collapse %d scaffolds\n", matcher.matchingSize());
#endif
		
		ContainmentGraph::UEdgeMap<bool> matched_edges(containment);
		vector<char> keep_scaffold(scaffolds.size(), 1);
		matcher.matching(matched_edges);
		
		vector<Scaffold> merged_scaffolds;
		
		vector<vector<DAGNode> > chains;
		
		// Make chains out of the matching
        ReachGraph::ANodeMap<ReachGraph::UEdge> matched_a_nodes(containment);
        matcher.aMatching(matched_a_nodes);
        
        ReachGraph::BNodeMap<ReachGraph::UEdge> matched_b_nodes(containment);
        matcher.bMatching(matched_b_nodes);
        
        set<ReachGraph::ANode> chain_heads;
        
        for (ReachGraph::ANodeIt i(containment); i!=lemon::INVALID; ++i)
        {
            int a_id = containment.aNodeId(i);
            ReachGraph::ANode a = containment.nodeFromANodeId(a_id);
            if (matched_a_nodes[a] == lemon::INVALID)
                chain_heads.insert(containment.nodeFromANodeId(containment.aNodeId(i)));
        }
        vector<DAGNode> held_out_of_collapse;
        for (set<ReachGraph::ANode>::iterator i = chain_heads.begin();
             i != chain_heads.end();
             ++i)
        {
            vector<DAGNode> chain;
            ReachGraph::ANode n = *i;
            chain.push_back(containment.aNodeId(*i));
            while(true)
            {
                //int a_id = bp.aNodeId(n);
                int b_id = containment.bNodeId(n);
                
                ReachGraph::BNode b = containment.nodeFromBNodeId(b_id);
                
                if (matched_b_nodes[b] == lemon::INVALID)
                {
                    break;
                }
                else
                {
                    ReachGraph::ANode a_match_to_b = containment.source(matched_b_nodes[b]);
                    if (chain.size() < max_chain_size)
                    {
                        chain.push_back(containment.aNodeId(a_match_to_b));
                    }
                    else
                    {
                        held_out_of_collapse.push_back(containment.aNodeId(a_match_to_b));  
                    }
                    n = a_match_to_b;
                }
            }
            chains.push_back(chain);
        }
        assert (chains.size() == chain_heads.size());
        
		for (size_t i = 0; i < chains.size(); ++i)
		{
			vector<Scaffold> chain;
			for (size_t j = 0; j < chains[i].size(); ++j)
			{
				chain.push_back(scaffolds[chains[i][j]]);
			}
			merged_scaffolds.push_back(Scaffold(chain));
		}
        
        for (size_t j = 0; j < held_out_of_collapse.size(); ++j)
        {
            merged_scaffolds.push_back(scaffolds[held_out_of_collapse[j]]);
        }
        
#if ASM_VERBOSE
		fprintf (stderr, "\t\tpre = %d, post = %d\n", (int)scaffolds.size(), (int)merged_scaffolds.size());
#endif
		sort (merged_scaffolds.begin(), merged_scaffolds.end(), scaff_lt);
		performed_collapse |= (scaffolds.size() != merged_scaffolds.size());
        scaffolds = merged_scaffolds;
		
        
	}
	return performed_collapse;
}

void compress_consitutive(vector<Scaffold>& hits)
{
    vector<bool> scaffold_mask;
    
    scaffold_mask = vector<bool>(hits.size(), false);
    add_non_constitutive_to_scaffold_mask(hits, scaffold_mask);
    
    vector<Scaffold> constitutive;
    vector<Scaffold> non_constitutive;
    
    for (size_t i = 0; i < scaffold_mask.size(); ++i)
    {
        if (!scaffold_mask[i])
            constitutive.push_back(hits[i]);
        else
            non_constitutive.push_back(hits[i]);
    }
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
}

bool make_scaffolds(int bundle_left,
					int bundle_length,
					vector<Scaffold>& hits,
					vector<Scaffold>& scaffolds)
{
	if (hits.empty())
		return true;

	int intron_hits = 0;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].has_intron())
		{
			intron_hits++;
		}
	}
	
	if (!intron_hits)
	{
#if ASM_VERBOSE
		fprintf(stderr, "No introns in bundle, collapsing all hits to single transcript\n");
#endif
		scaffolds.push_back(Scaffold(hits));
		fill_gaps(scaffolds, 2 * olap_radius);
	}
	else
	{
#if ASM_VERBOSE
		fprintf(stderr, "Bundle has %d spliced reads\n", intron_hits);
#endif
//#if ASM_VERBOSE
		FILE* f_verbose = stderr;
//#endif

		vector<Scaffold> hazards;
		holdout_transitivity_hazards(hits, hazards);
		
        vector<Scaffold> split_hazards;
        // Cleave the partials at their unknowns to minimize FPKM dilation on  
        // the low end of the expression profile. 
        for (size_t i = 0; i < hazards.size(); ++i) 
        { 
            vector<Scaffold> c; 
            hazards[i].get_complete_subscaffolds(c); 
            split_hazards.insert(split_hazards.end(), c.begin(), c.end()); 
        } 
        
		vector<Scaffold> orig_hits = hits;
		
        hits.insert(hits.end(), split_hazards.begin(), split_hazards.end());

        compress_consitutive(hits);
        
#if ASM_VERBOSE
        fprintf(f_verbose,"\tPerforming preliminary containment collapse on %lu fragments\n", hits.size());
        size_t pre_hit_collapse_size = hits.size();
#endif
		if (perform_full_collapse)
		{
            collapse_contained_scaffolds(hits);
		}
			
        
#if ASM_VERBOSE
		size_t post_hit_collapse_size = hits.size();
		fprintf(f_verbose,"\tIgnoring %lu strictly contained scaffolds\n", pre_hit_collapse_size - post_hit_collapse_size);
#endif
        
        fprintf(stderr, "Assembling bundle with %lu hits\n", hits.size());
		
		vector<int> depth_of_coverage(bundle_length,0);
		map<pair<int,int>, int> intron_depth_of_coverage;
		compute_doc(bundle_left, 
					hits, 
					depth_of_coverage, 
					intron_depth_of_coverage,
					false);

		normal norm(0, 0.1);
		
		bool first_round = false;
		//bool first_round = false;
		vector<const MateHit*> prev_chaff;
		while (true)
		{		
			static size_t MAX_BUNDLE_ALIGNMENTS = 0xFFFF;
			
			sort(hits.begin(), hits.end(), scaff_lt);
			
			if (hits.size() >= MAX_BUNDLE_ALIGNMENTS)
			{
#if ASM_VERBOSE
				fprintf(f_verbose, "Warning: bundle too large, skipping assembly\n");
#endif
				return false;
			}
           
            fprintf(stderr, "\tCurrent bundle has %d non-constitutive fragments\n", hits.size());            
			
			DAG bundle_dag;
			
			if (hits.empty())
				return true;
#if ASM_VERBOSE
			fprintf(f_verbose, "\tCalculating scaffold densities\n");
#endif
			vector<double> scaff_doc;
			record_doc_for_scaffolds(bundle_left, 
									 hits, 
									 depth_of_coverage, 
									 intron_depth_of_coverage,
									 scaff_doc);
#if ASM_VERBOSE
			fprintf(f_verbose, "\tCreating compatibility graph\n");
#endif
			if (!create_dag_for_bundle(hits, bundle_dag))
			{
				//scaffolds = hits;
				//create_dag_for_bundle(hits, bundle_dag, first_round ? -2 : -999999);
				break;
			}
            
            HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);

            
            vector_property_map<size_t> path_for_scaff;
            path_compress_visitor<vector_property_map<DAGNode>,
                                  vector_property_map<size_t> > vis(path_for_scaff);
            depth_first_search(bundle_dag,
                               visitor(vis));
            
           
            vector<vector<Scaffold> > compressed_paths(hits.size()+1);
            {
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
                create_dag_for_bundle(hits, bundle_dag);
            }
			
			
			vector<char> has_parent(num_vertices(bundle_dag) + 2, false);
			vector<char> has_child (num_vertices(bundle_dag) + 2, false);
			
			graph_traits < DAG >::vertex_iterator u, uend;
			for (tie(u, uend) = vertices(bundle_dag); u != uend; ++u)
			{
				graph_traits < DAG >::adjacency_iterator v, vend;
				for (tie(v,vend) = adjacent_vertices(*u, bundle_dag); v != vend; ++v)
				{
					DAGNode U = *u;
					DAGNode V = *v;
					has_parent[V] = true;
					has_child[U] = true;
				}
			}
			
			DAGNode source = add_vertex(bundle_dag);
			DAGNode sink = add_vertex(bundle_dag);
			
#ifdef DEBUG
			set<const Scaffold*> introns;
#endif
			for (size_t i = 0; i < num_vertices(bundle_dag); ++i)
			{
//				
//#if ASM_VERBOSE
				const Scaffold* hit = hits_for_node[i];
				if (hit && hit->has_intron())
				{
//					//add_edge(source, i, bundle_dag);
//					//add_edge(i, sink, bundle_dag);
//
//					introns.insert(hit);
//					fprintf(stderr, "%d, %d\n", hit->left(), hit->right());
//
				}
//#endif
				if (!has_parent[i] && i != sink && i != source)
					add_edge(source, i, bundle_dag);
				
				if (!has_child[i] && i != source && i != sink)
					add_edge(i, sink, bundle_dag);
			}
			
			ReachGraph bp;
			
//#if ASM_VERBOSE
			fprintf(f_verbose, "\tConstructing reachability graph\n");
//#endif
			vector<ReachGraph::BNode> b_to_a;
			adjacency_list<> TC;
			
			transitive_closure(bundle_dag, TC);
			DagToBp dag_to_bp;
			vector<bool> scaffold_mask(num_vertices(bundle_dag), true);
			            // Add ends to scaffold mask
//            for (size_t i = 0; i < num_vertices(bundle_dag); ++i)
//			{
//				if (!has_parent[i] && i != sink && i != source)
//					scaffold_mask[i] = true;
//				
//				if (!has_child[i] && i != source && i != sink)
//					scaffold_mask[i] = true;
//			}
            
			reachability_graph(bundle_dag, source, sink,  bp, b_to_a, dag_to_bp, TC, scaffold_mask);
			
			ReachGraph::UEdgeMap<long long> cov_weights(bp);
			add_weights_to_reachability_graph(bp, hits_for_node, orig_hits, hits, cov_weights);				
			
			vector<vector<DAGNode> > chains;
			if (first_round)
			{
#if ASM_VERBOSE
				fprintf(f_verbose, "\tPerforming initial intronic matching\n");
#endif
				typedef lemon::MinCostMaxBipartiteMatching<ReachGraph,ReachGraph::UEdgeMap<long long> > InitialMatcher;
				InitialMatcher init_matcher(bp, cov_weights);
				init_matcher.run();

#if ASM_VERBOSE
				fprintf(f_verbose, "\tConverting matching into chain decomposition\n");
#endif
				make_chains_from_matching<InitialMatcher>(bp, init_matcher, chains);

#if ASM_VERBOSE
				fprintf(f_verbose, "\tConverting matching into chain decomposition\n");
#endif
			}
			else
			{
//#if ASM_VERBOSE
				fprintf(f_verbose, "\tPerforming weighted matching\n");
//#endif
				typedef lemon::MinCostMaxBipartiteMatching<ReachGraph,ReachGraph::UEdgeMap<long long> > Matcher;
				Matcher matcher(bp, cov_weights);
				matcher.run();
				make_chains_from_matching<Matcher>(bp, matcher, chains);
			}
//#if ASM_VERBOSE
			fprintf(f_verbose, "\tFound %d distinct chains\n", (int)chains.size());
//#endif		
			vector<vector<DAGNode> > paths;
			extend_chains_to_paths(bundle_dag, chains, TC, source, sink, paths);
			
	#ifdef DEBUG
			set<const Scaffold*> path_introns;
			for (size_t p = 0; p < paths.size(); ++p)
			{
				const vector<DAGNode>& path = paths[p];
				for (size_t m = 0; m < path.size(); ++m)
				{
					const Scaffold* hit = hits_for_node[path[m]];
					if (hit && first_round && hit->has_intron())
					{
						path_introns.insert(hit);
					}
				}
			}
			
			//assert (introns == path_introns);
	#endif

//#if ASM_VERBOSE
			fprintf(f_verbose, "\tCreating scaffolds for %d paths\n", (int)paths.size());
//#endif
			vector<Scaffold> new_scaffs;
			make_scaffolds_from_paths(bundle_dag, paths, new_scaffs);
//#if ASM_VERBOSE
			fprintf(f_verbose, "\tCollapsing scaffolds\n");
//#endif
			collapse_contained_scaffolds(new_scaffs);
			hits = new_scaffs;
			first_round = false;
		}
	
		scaffolds = hits;
		
		// One last collapse attempt...
		vector<Scaffold> new_scaffs = scaffolds;
#if ASM_VERBOSE
		fprintf(stderr, "Performing final collapse round\n");
#endif
		
		//remove_improbable_scaffolds(new_scaffs, prev_chaff);
		fill_gaps(new_scaffs, 2 * olap_radius);
        
        scaffolds = new_scaffs;
        
        // Cleave the partials at their unknowns to minimize FPKM dilation on  
        // the low end of the expression profile. 
        vector<Scaffold> completes; 
        for (size_t i = 0; i < scaffolds.size(); ++i) 
        { 
            vector<Scaffold> c; 
            scaffolds[i].get_complete_subscaffolds(c); 
            completes.insert(completes.end(), c.begin(), c.end()); 
        } 
        new_scaffs = completes;
        sort(new_scaffs.begin(), new_scaffs.end(), scaff_lt);
        
		collapse_contained_scaffolds(new_scaffs);
		sort(new_scaffs.begin(), new_scaffs.end(), scaff_lt);
		scaffolds = new_scaffs;
	}
    
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		//assert(!scaffolds[i].has_unknown());
		
		const vector<const MateHit*>& supporting = scaffolds[i].mate_hits();
		CuffStrand s = CUFF_STRAND_UNKNOWN;
		for (size_t j = 0; j < supporting.size(); ++j)
		{

//			assert (supporting[j]->strand() == CUFF_STRAND_UNKNOWN || 
//					s == supporting[j]->strand());
			if (supporting[j]->strand() != CUFF_STRAND_UNKNOWN)
				s = supporting[j]->strand();
		}
        if (scaffolds[i].strand() == CUFF_STRAND_UNKNOWN)
            scaffolds[i].strand(s);
		
//		if (!scaffolds[i].has_intron())
//			scaffolds[i].strand(CUFF_STRAND_UNKNOWN);
	}
	
	return true;
}



