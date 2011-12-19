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
#include <iostream>

#include <lemon/topology.h>

//#include <lemon/graph_utils.h>
#include <lemon/bipartite_matching.h>

#include <boost/ref.hpp>
// DON'T move this, or mystery compiler errors will result. Affects gcc >= 4.1
#include <boost/graph/vector_as_graph.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/graph_traits.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <boost/version.hpp>

#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include <boost/property_map/vector_property_map.hpp>
#endif

#include <boost/math/distributions/normal.hpp> // for normal_distribution
using boost::math::normal; // typedef provides default type is double.

#include "transitive_closure.h"
#include "transitive_reduction.h"

#include "common.h"
#include "assemble.h"

#include "bundles.h"
//#include "filters.h"
//#include "genes.h"
#include "scaffolds.h"
#include "clustering.h"
#include "matching_merge.h"
#include "graph_optimize.h"

using namespace boost;
using namespace std;

bool mate_graphs(const HitBundle& bundle, BundleStats* stats);

typedef lemon::SmartBpUGraph ReachGraph;

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
	}

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
	return weight;
}


typedef map<DAGNode, pair<int, int> > DagToBp;
void create_reachability_bp_graph(DAG& dag,
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
		
		if (in_degree(s, dag) == 0  /* virtual "source"? */|| 
            out_degree(t, dag) == 0 /* virtual "sink"?*/)
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

void add_weights_to_reachability_bp_graph(ReachGraph& bp,
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
			if (!hits[j].is_ref() && scaffolds[i].contains(hits[j]))
			{
				if (Scaffold::compatible(scaffolds[i],hits[j]))
				{
					spliced_in[i]++;
				}
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
				if (Scaffold::compatible(scaffolds[i], scaffolds[j]))
				{
					compatible_density += density[j];
				}
			}
		}
		if (total_density)
			psi[i] = compatible_density / total_density;
		local_density[i] = compatible_density;
	}
	
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
				evil_introns.push_back(introns[i]);
				evil_introns.push_back(introns[j]);
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
					}
				}				  
			}
		}
		if (overlaps_evil_intron)
		{
//            if (hits[i].has_intron())
//            {
//                fprintf(stderr, "&&& Holding out intron-containing hazard at %d-%d\n", hits[i].left(), hits[i].right());
//            }
			hazards.push_back(hits[i]);
		}	
		else
		{
			filtered_hits.push_back(hits[i]);
		}
	}
	
	verbose_msg( "%s\tHeld out %lu scaffolds as transitivity hazards\n", bundle_label->c_str(), hazards.size());
	
	hits = filtered_hits;
}

bool make_scaffolds(int bundle_left,
					int bundle_length,
					vector<Scaffold>& hits,
					vector<Scaffold>& scaffolds)
{
	if (hits.empty())
		return true;

	bool intron_hits = false;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].has_intron())
		{
			intron_hits = true;
			break;
		}
	}
	
	if (!intron_hits)
	{
		verbose_msg( "%s\tNo introns in bundle, collapsing all hits to single transcript\n", bundle_label->c_str());
		scaffolds.push_back(Scaffold(hits));
		fill_gaps(scaffolds, 2 * olap_radius);
	}
	else
	{
		verbose_msg( "%s\tBundle has spliced reads\n", bundle_label->c_str());

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
		

        compress_fragments(hits);

        verbose_msg( "%s\tAssembling bundle with %lu hits\n", bundle_label->c_str(), hits.size());
       
		vector<float> depth_of_coverage(bundle_length,0);
		map<pair<int,int>, float> intron_depth_of_coverage;
		compute_doc(bundle_left, 
					hits, 
					depth_of_coverage, 
					intron_depth_of_coverage,
					false);

		normal norm(0, 0.1);
		
		vector<const MateHit*> prev_chaff;
		while (true)
		{		
			static size_t MAX_BUNDLE_ALIGNMENTS = 0xFFFF;
			
			sort(hits.begin(), hits.end(), scaff_lt);
           
            //fprintf(stderr, "\tCurrent bundle has %d non-constitutive fragments\n", hits.size());            
			
			DAG bundle_dag;
			
			if (hits.empty())
				return true;
			verbose_msg( "%s\tCalculating scaffold densities\n", bundle_label->c_str());
			vector<double> scaff_doc;
			record_doc_for_scaffolds(bundle_left, 
									 hits, 
									 depth_of_coverage, 
									 intron_depth_of_coverage,
									 scaff_doc);
			verbose_msg( "%s\tCreating compatibility graph\n", bundle_label->c_str());
			
			if (!create_overlap_dag(hits, bundle_dag))
			{			
				break;
			}
            
            HitsForNodeMap hits_for_node = get(vertex_name, bundle_dag);

            compress_overlap_dag_paths(bundle_dag, hits);
            
            if (hits.size() >= MAX_BUNDLE_ALIGNMENTS)
			{
				verbose_msg( "%s\tWarning: bundle too large, skipping assembly\n", bundle_label->c_str());
				return false;
			}
            
            pair<DAGNode, DAGNode> terminal = add_terminal_nodes(bundle_dag);
            DAGNode source = terminal.first;
            DAGNode sink = terminal.second;
            
			ReachGraph bp;
			
			verbose_msg( "%s\tConstructing reachability graph\n", bundle_label->c_str());
			
			vector<ReachGraph::BNode> b_to_a;
			adjacency_list<> TC;
			
			transitive_closure(bundle_dag, TC);
			DagToBp dag_to_bp;
            
            // TODO: deprecate dependence of create_reachability_bp_graph() on scaffold_mask
			vector<bool> scaffold_mask(num_vertices(bundle_dag), true);
            
			create_reachability_bp_graph(bundle_dag, bp, b_to_a, dag_to_bp, TC, scaffold_mask);
			
			ReachGraph::UEdgeMap<long long> cov_weights(bp);
			add_weights_to_reachability_bp_graph(bp, hits_for_node, orig_hits, hits, cov_weights);				

            verbose_msg( "%s\tPerforming weighted matching\n", bundle_label->c_str());

            typedef lemon::MinCostMaxBipartiteMatching<ReachGraph,ReachGraph::UEdgeMap<long long> > Matcher;
            Matcher matcher(bp, cov_weights);
            matcher.run();
            
            vector<vector<DAGNode> > chains;
            make_chains_from_matching<Matcher>(bp, matcher, chains);
			
			verbose_msg( "%s\tFound %d distinct chains\n", bundle_label->c_str(), (int)chains.size());
			
			vector<vector<DAGNode> > paths;
			extend_chains_to_paths(bundle_dag, chains, TC, source, sink, paths);

			verbose_msg( "%s\tCreating scaffolds for %d paths\n", bundle_label->c_str(), (int)paths.size());
			
			vector<Scaffold> new_scaffs;
			make_scaffolds_from_paths(bundle_dag, paths, new_scaffs);
			
			verbose_msg( "%s\tCollapsing scaffolds\n", bundle_label->c_str());
           
			collapse_contained_transfrags(new_scaffs);
			hits = new_scaffs;
		}
	
		scaffolds = hits;
		
		// One last collapse attempt...
		vector<Scaffold> new_scaffs = scaffolds;
		
		verbose_msg( "%s\tPerforming final collapse round\n", bundle_label->c_str());
		
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
        
        verbose_msg( "Extracted %lu contiguous transfrags from %lu scaffolds\n", completes.size(), scaffolds.size());
        
        new_scaffs = completes;
        sort(new_scaffs.begin(), new_scaffs.end(), scaff_lt);
        
		collapse_contained_transfrags(new_scaffs);
		sort(new_scaffs.begin(), new_scaffs.end(), scaff_lt);
		scaffolds = new_scaffs;
	}
    
    // TODO: refactor into subroutine.  This routine shouldn't actually be 
    // necessary, and should be eliminated after thorough testing that it really
    // isn't needed
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
	}
    // end refactor
	
	return true;
}
