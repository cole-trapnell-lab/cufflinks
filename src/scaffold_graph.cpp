/*
 *  scaffold_graph.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/2/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>
#include "scaffold_graph.h"
#include "scaffolds.h"

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>

#ifndef NDEBUG
#include "transitive_reduction.h"
#endif

using namespace std;
using namespace boost;

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

template <class CompatibilityMap, class ConnectMap, class Tag>
connect_visitor<CompatibilityMap, ConnectMap, Tag>
record_connections(CompatibilityMap compatibility, 
				   ConnectMap connect, 
				   DAGNode target, 
				   Tag) 
{
    return connect_visitor<CompatibilityMap, ConnectMap, Tag> (compatibility, connect, target);
}


bool create_overlap_dag(vector<Scaffold>& hits,
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
				if (Scaffold::compatible(lhs, rhs))
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
    
#ifndef NDEBUG
    DAG tr;
    boost::vector_property_map<DAGNode> G_to_TR;
    property_map<DAG, vertex_index_t>::type w = get(vertex_index, bundle_dag);
    transitive_reduction(bundle_dag, 
                         tr, 
                         G_to_TR,
                         w);
    verbose_msg("dag has %lu edges, tr has %lu edges\n", num_edges(bundle_dag), num_edges(tr));
    
	//assert (num_edges(bundle_dag) == num_edges(tr));
#endif
    
	return found_compatible_scaffolds;
}

pair<DAGNode, DAGNode> add_terminal_nodes(DAG& bundle_dag)
{
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
    
    int num_attached_to_source = 0;
    int num_attached_to_sink = 0;
    
    for (size_t i = 0; i < num_vertices(bundle_dag); ++i)
    {
        if (!has_parent[i] && i != sink && i != source)
        {
            num_attached_to_source++;
            add_edge(source, i, bundle_dag);
        }
        if (!has_child[i] && i != source && i != sink)
        {
            num_attached_to_sink++;
            add_edge(i, sink, bundle_dag);
        }
    }
    
#if verbose_msg
    HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);
    DAG::vertex_iterator ki, ke;
	for (tie(ki, ke) = vertices(bundle_dag); ki != ke; ++ki)
	{
        if (edge(source, *ki, bundle_dag).second)
        {
            const Scaffold* pS = hits_for_node[*ki];
            fprintf(stderr, "%d-%d has edge from source\n", pS->left(), pS->right());
        }
        
        if (edge(*ki, sink, bundle_dag).second)
        {
            const Scaffold* pS = hits_for_node[*ki];
            fprintf(stderr, "%d-%d has edge to sink\n", pS->left(), pS->right());
        }
    }
    verbose_msg("%d source nodes, %d sink nodes\n", num_attached_to_source, num_attached_to_sink);
#endif
    return make_pair(source, sink);
}
