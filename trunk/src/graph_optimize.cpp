/*
 *  graph_optimize.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 6/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>
#include <algorithm>

#include "graph_optimize.h"
// for graph optimization only

#include <lemon/topology.h>
#include <lemon/smart_graph.h>
#include <lemon/bipartite_matching.h>

#include "scaffold_graph.h"
#include "scaffolds.h"
#include "filters.h"
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

bool scaff_left_lt_right_gt(const Scaffold& lhs, const Scaffold& rhs)
{
    if (lhs.left() != rhs.left())
        return lhs.left() < rhs.left();
    return lhs.right() > rhs.right();
}

bool op_left_lt_right_lt(const AugmentedCuffOp& lhs, const AugmentedCuffOp& rhs)
{
	if (lhs.genomic_offset != rhs.genomic_offset)
	{
		return lhs.genomic_offset < rhs.genomic_offset;
	}
	if (lhs.genomic_length != rhs.genomic_length)
	{
		return lhs.genomic_length < rhs.genomic_length;
	}
	return false;
}

void extract_conflicting_ops(const vector<AugmentedCuffOp>& ops,
                                       vector<AugmentedCuffOp>& conflict_ops)
{	
	for (size_t i = 0; i < ops.size(); ++i)
	{
		for (size_t j = i+1; j < ops.size(); ++j)
		{
			if (AugmentedCuffOp::overlap_in_genome(ops[i], ops[j]))
			{
				if (!AugmentedCuffOp::compatible(ops[i], ops[j]))
				{
					if (!binary_search(conflict_ops.begin(), conflict_ops.end(), ops[i]))
					{
						conflict_ops.push_back(ops[i]);
						sort(conflict_ops.begin(), conflict_ops.end());
					}
					
					if (!binary_search(conflict_ops.begin(), conflict_ops.end(), ops[j]))
					{
						conflict_ops.push_back(ops[j]);
						sort(conflict_ops.begin(), conflict_ops.end());
					}
				}
			}
			else
			{
				break;
			}
            
		}
	}    
}

void collect_non_redundant_ops(const vector<Scaffold>& scaffolds,
                               vector<AugmentedCuffOp>& ops)
{

	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		ops.insert(ops.end(), 
				   scaffolds[i].augmented_ops().begin(), 
				   scaffolds[i].augmented_ops().end());
	}
    sort(ops.begin(), ops.end());
    vector<AugmentedCuffOp>::iterator new_end = unique(ops.begin(), ops.end());
	ops.erase(new_end, ops.end());
	
	sort (ops.begin(), ops.end(), op_left_lt_right_lt);   
}

void fill_unambiguous_unknowns(vector<Scaffold>& to_fill, 
                               const vector<Scaffold>& constitutive)
{
//    vector<AugmentedCuffOp> conflict_ops;
//    vector<AugmentedCuffOp> ops;

    for (size_t i = 0; i < to_fill.size(); ++i)
    {
        if (to_fill[i].has_unknown())
        {
            for( size_t j = 0; j < constitutive.size(); ++j)
            {
                const Scaffold& cons = constitutive[j];
                if (Scaffold::overlap_in_genome(to_fill[i], cons, 0) &&
                    Scaffold::compatible(to_fill[i], cons))
                {
                    if (cons.strand() != CUFF_STRAND_UNKNOWN)
                        to_fill[i].strand(cons.strand());

                    to_fill[i].fill_gaps(cons.augmented_ops());
                    if (!to_fill[i].has_unknown())
                    {
                        break;
                    }
                }
            }
        }
    }
}

// WARNING: scaffolds MUST be sorted by scaff_lt_rt() in order for this routine
// to work correctly.
void add_non_constitutive_to_scaffold_mask(const vector<Scaffold>& scaffolds,
										   vector<bool>& scaffold_mask)
{	
	
	// First, we filter out all fragments that are entirely contained in a 
	// constitutive exon of the "gene".  If such fragments were non-constitutive,
	// neither would that exon.  Also, we can examine all fragments at most
	// once here, because if we don't look at them here, we'll look hard in the
	// next loop
	Scaffold smashed_gene;
	
	// setting introns_overwrite_matches in a gene smash takes only it's
	// constitutive regions
	Scaffold::merge(scaffolds, smashed_gene, true); 
	
	vector<bool> smash_filter(scaffolds.size(), false);
	
	const vector<AugmentedCuffOp>& cig = smashed_gene.augmented_ops();
	size_t next_frag = 0; 
	
	size_t num_filtered = 0;
	for (size_t j = 0; j < cig.size(); ++j)
	{
		if (cig[j].opcode == CUFF_MATCH)
		{
			for (;next_frag < scaffolds.size(); ++next_frag)
			{
				const Scaffold& frag = scaffolds[next_frag];
				
				if (frag.left() >= cig[j].g_left() && 
					frag.right() <= cig[j].g_right())
				{
					smash_filter[next_frag] = true;
					//scaffold_mask[next_frag] = true;
					num_filtered++;
				}
				if (frag.left() >= cig[j].g_right())
				{
					break;
				}
			}
		}
	}
    
	verbose_msg("%lu constitutive reads of %lu smash-filtered from further consideration\n", num_filtered, smash_filter.size());
    
    vector<AugmentedCuffOp> ops;
    collect_non_redundant_ops(scaffolds, ops);
	
    vector<AugmentedCuffOp> conflict_ops;
    extract_conflicting_ops(ops, conflict_ops);
	
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		if (smash_filter[i])
			continue;
		const vector<AugmentedCuffOp>& s_ops = scaffolds[i].augmented_ops();
		for (size_t j = 0; j < s_ops.size(); ++j)
		{
			if (binary_search(conflict_ops.begin(), conflict_ops.end(), s_ops[j]))
			{
				scaffold_mask[i] = true;
				break;
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
		
		verbose_msg("%s\tStarting new collapse round\n", bundle_label->c_str());
        
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
        
		verbose_msg("%s\tContainment graph has %d nodes, %d edges\n", bundle_label->c_str(), containment.aNodeNum(), containment.uEdgeNum());
		verbose_msg("%s\tFinding a maximum matching to collapse scaffolds\n", bundle_label->c_str());
        
		matcher.run();
        
		verbose_msg( "%s\tWill collapse %d scaffolds\n", bundle_label->c_str(), matcher.matchingSize());
		
		ContainmentGraph::UEdgeMap<bool> matched_edges(containment);
		
		matcher.matching(matched_edges);
		
        merge_from_matching(containment, matcher, scaffolds);
        
		performed_collapse = true;
	}
	return performed_collapse;
}

bool scaff_smaller_lt_rt(const Scaffold& lhs, const Scaffold& rhs)
{
	size_t lhs_len = lhs.right() - lhs.left();
	size_t rhs_len = rhs.right() - rhs.left();
	
	if (lhs_len != rhs_len)
	{
		return lhs_len < rhs_len;
	}
	else 
	{
		return scaff_lt_rt(lhs, rhs);
	}
	return false;
}

struct FragIndexSortSmallerLR
{
	FragIndexSortSmallerLR(const vector<Scaffold>& frags) : fragments(frags) {}
	
	const vector<Scaffold>& fragments;
	
	bool operator()(size_t lhs_frag_idx, size_t rhs_frag_idx)
	{
		const Scaffold& lhs = fragments[lhs_frag_idx];
		const Scaffold& rhs = fragments[rhs_frag_idx];
		
		size_t lhs_len = lhs.right() - lhs.left();
		size_t rhs_len = rhs.right() - rhs.left();
		
		if (lhs_len != rhs_len)
		{
			return lhs_len > rhs_len;
		}
		else 
		{
			return scaff_lt_rt(lhs, rhs);
		}
		return false;
	}
};

bool collapse_equivalent_transfrags(vector<Scaffold>& fragments, 
                                    uint32_t max_rounds)
{
	// The containment graph is a bipartite graph with an edge (u,v) when
	// u is (not necessarily properly) contained in v and is the two are
	// compatible.
	typedef lemon::SmartBpUGraph ContainmentGraph;
	normal norm(0, 0.1);
	bool performed_collapse = false;
	
	//double last_size = -1;
    //long leftmost = 9999999999;
    long leftmost = LONG_MAX-1;
    long rightmost = -1;
    
    for (size_t i = 0; i < fragments.size(); ++i)
    {
        leftmost = std::min((long)fragments[i].left(), leftmost);
        rightmost = std::max((long)fragments[i].right(), rightmost);
    }
    
    //long bundle_length = rightmost - leftmost;
    
	while (max_rounds--)
	{
		
		sort (fragments.begin(), fragments.end(), scaff_lt_rt);
		
		
		vector<size_t> smaller_idx_array;
		for (size_t i = 0; i < fragments.size(); ++i)
        {
            smaller_idx_array.push_back(i);
        }
		
		sort(smaller_idx_array.begin(), 
			 smaller_idx_array.end(), 
			 FragIndexSortSmallerLR(fragments));
		
		verbose_msg("%s\tStarting new collapse round\n", bundle_label->c_str());
        verbose_msg("%s\tFinding fragment-level conflicts\n", bundle_label->c_str());

        bool will_perform_collapse = false;
		
        verbose_msg( "%s\tAssessing overlaps between %lu fragments for identical conflict sets\n", 
                bundle_label->c_str(), 
                fragments.size());
        vector<size_t> replacements;
        for (size_t i = 0; i < fragments.size(); ++i)
        {
            replacements.push_back(i);
        }
		
		size_t curr_frag = 0;
		vector<int> curr_conflicts;
		
//		for (int i = 0; i < fragments.size(); ++i)
//		{
//			if (Scaffold::overlap_in_genome(fragments[0], fragments[i], 0))
//			{
//				if (!Scaffold::compatible(fragments[0], fragments[i]))
//				{
//					curr_conflicts.push_back(i);
//				}
//			}
//		}
		
        double mean_length = 0;
        for (size_t i = 0; i < fragments.size(); ++i)
        {
            mean_length += fragments[i].length();
        }
        
        mean_length /= fragments.size();
        
        double variance = 0.0;
        for (size_t i = 0; i < fragments.size(); ++i)
        {
            double v = fragments[i].length() - mean_length;
            v *= v;
            variance += v;
        }
        
        variance /= fragments.size();
        double frag_len_std_dev = sqrt(variance);
        
		int num_merges = 0;
		
		while (curr_frag < smaller_idx_array.size())
		{	
			size_t curr_frag_native_idx = smaller_idx_array[curr_frag];
			if (replacements[curr_frag_native_idx] == curr_frag_native_idx)
			{
				size_t lhs = curr_frag;
				
				size_t lhs_native_idx = smaller_idx_array[lhs];
				
				const Scaffold& lhs_scaff = fragments[lhs_native_idx];
				curr_conflicts.clear();
				
				double lhs_len = lhs_scaff.right() - lhs_scaff.left();
				
				for (size_t i = 0; i < smaller_idx_array.size(); ++i)
				{
					size_t j_scaff_idx = smaller_idx_array[i];
                    if (replacements[j_scaff_idx] == j_scaff_idx)
                    {
                        if (Scaffold::overlap_in_genome(lhs_scaff, fragments[j_scaff_idx], 0))
                        {
                            if (!Scaffold::compatible(lhs_scaff, fragments[j_scaff_idx]))
                            {
                                curr_conflicts.push_back(j_scaff_idx);
                            }
                        }
                    }
				}
				sort(curr_conflicts.begin(), curr_conflicts.end());
				
				//bool advanced_curr = false;
				for (size_t c = lhs + 1; c < smaller_idx_array.size(); ++c)
				{
					size_t c_native_idx = smaller_idx_array[c];
					const Scaffold& c_scaff = fragments[c_native_idx];
					if (replacements[c_native_idx] == c_native_idx &&
						lhs_scaff.contains(c_scaff))
					{
						double c_len = c_scaff.right() - c_scaff.left();

						if (lhs_len - c_len > frag_len_std_dev)
							break;

						if (c_scaff.augmented_ops() == lhs_scaff.augmented_ops())
						{
							if (num_merges % 100 == 0)
							{
								verbose_msg("%s\tCollapsing frag # %d\n", 
										bundle_label->c_str(), 
										num_merges);
							}
							vector<Scaffold> s;
							s.push_back(c_scaff);
							s.push_back(lhs_scaff);
							fragments[lhs_native_idx] = Scaffold(s);
							replacements[c_native_idx] = lhs_native_idx;
							//fragments[c_native_idx] = Scaffold();
							//curr_conflicts = c_conflicts;
							fragments[c_native_idx].clear_hits();
							///lhs = c;
							//advanced_curr = true;
							will_perform_collapse = true;
							num_merges++;
							continue;
						}

						if (!Scaffold::compatible(lhs_scaff, c_scaff))
							continue;
						vector<int> c_conflicts;
						// Find c's conflicts
						
						// If c fails to overlap lhs's conflicts, or if it's
						// compatible with any of them, they aren't equivalent
						bool not_equivalent = false;
						for (size_t j = 0; j < curr_conflicts.size(); ++j)
						{
							if (!Scaffold::overlap_in_genome(fragments[curr_conflicts[j]], c_scaff, 0) ||
								 Scaffold::compatible(fragments[curr_conflicts[j]], c_scaff))
							{
								not_equivalent = true;
								break;
							}
						}
						
						if (not_equivalent)
							continue;
						
						// If we get here, then c disagrees with at least all 
						// the guys lhs does.
					
						// Now check that c doesn't have any additional conflicts
						// of it's own
						for (size_t i = lhs_native_idx + 1; i < fragments.size(); ++i)
						{
                            if (replacements[i] == i)
                            {
                                if (Scaffold::overlap_in_genome(fragments[i], lhs_scaff, 0))
                                {
                                    if (Scaffold::overlap_in_genome(fragments[i], c_scaff, 0))
                                    {
                                        if (!Scaffold::compatible(fragments[i], c_scaff))
                                        {
                                            //c_conflicts.push_back(i);
                                            if (!binary_search(curr_conflicts.begin(), curr_conflicts.end(), i))
                                            {
                                                not_equivalent = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                            }
						}
						
						
						if (not_equivalent)
							continue;
					
						// merge
						if (num_merges % 100 == 0)
						{
							verbose_msg("%s\tCollapsing frag # %d\n", 
									bundle_label->c_str(), 
									num_merges);
						}
						
						vector<Scaffold> s;
						s.push_back(c_scaff);
						s.push_back(lhs_scaff);
						fragments[lhs_native_idx] = Scaffold(s);
						replacements[c_native_idx] = lhs_native_idx;
						//fragments[c_native_idx] = Scaffold();
						fragments[c_native_idx].clear_hits();
						//curr_conflicts = c_conflicts;
						//advanced_curr = true;
						will_perform_collapse = true;
						num_merges++;
						//break;
					}
					else 
					{
						continue;
					}
					
				}
			}
			
			//if (!advanced_curr)
			{
				++curr_frag;
			}
		}
		
		if (will_perform_collapse == false)
			return performed_collapse;
		
		
        vector<Scaffold> replaced;
        for (size_t i = 0; i < fragments.size(); ++i)
        {
            if (replacements[i] == i)
            {
                replaced.push_back(fragments[i]);
            }
        }
        
        fragments = replaced;
        sort(fragments.begin(), fragments.end(), scaff_lt_rt);
		performed_collapse = true;
	}
	return performed_collapse;
}

void compress_consitutive(vector<Scaffold>& hits)
{
    vector<bool> scaffold_mask;
    
    verbose_msg("%s\tBuilding constitutivity mask\n", bundle_label->c_str()); 
    
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
    
    size_t pre_compress = hits.size();
    hits.clear();
    if (!constitutive.empty())
    {
        Scaffold compressed = Scaffold(constitutive); 
        vector<Scaffold> completes; 
        compressed.fill_gaps(2 * olap_radius);
        compressed.get_complete_subscaffolds(completes);
        
        hits.insert(hits.end(), completes.begin(), completes.end()); 
    }
    
    fill_unambiguous_unknowns(non_constitutive, hits);
    
    hits.insert(hits.end(), non_constitutive.begin(), non_constitutive.end());
    sort(hits.begin(), hits.end(), scaff_lt);
    
    
    size_t post_compress = hits.size();
    size_t delta = pre_compress - post_compress;
    double collapse_ratio = delta / (double) pre_compress; 
    verbose_msg("%s\tCompressed %lu of %lu constitutive fragments (%lf percent)\n", 
            bundle_label->c_str(),
            delta, 
            pre_compress, 
            collapse_ratio);
}


void compress_redundant(vector<Scaffold>& fragments)
{
    double last_size = -1;
    //long leftmost = 9999999999;
    long leftmost = LONG_MAX-1;
    long rightmost = -1;
    
    for (size_t i = 0; i < fragments.size(); ++i)
    {
        leftmost = std::min((long)fragments[i].left(), leftmost);
        rightmost = std::max((long)fragments[i].right(), rightmost);
    }
    
    size_t pre_compress = fragments.size();
    
    while (true)
    {
        if (last_size == -1 || 0.9 * last_size > fragments.size())
        {
            last_size = fragments.size();
            if (!collapse_equivalent_transfrags(fragments, 1))
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
    
    size_t post_compress = fragments.size();
    size_t delta = pre_compress - post_compress;
    double collapse_ratio = delta / (double) pre_compress; 
    verbose_msg("%s\tCompressed %lu of %lu redundant fragments (%lf percent)\n", 
                bundle_label->c_str(),
                delta, 
                pre_compress, 
                collapse_ratio);
}

void compress_fragments(vector<Scaffold>& fragments)
{
    verbose_msg("%s\tPerforming preliminary containment collapse on %lu fragments\n", bundle_label->c_str(), fragments.size());
    size_t pre_hit_collapse_size = fragments.size();
    sort(fragments.begin(), fragments.end(), scaff_lt_rt);
    
	compress_consitutive(fragments);
	
	compress_redundant(fragments);
    
    size_t post_hit_collapse_size = fragments.size();
    verbose_msg("%s\tIgnoring %lu strictly contained fragments\n", bundle_label->c_str(), pre_hit_collapse_size - post_hit_collapse_size);
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
            Scaffold s(compressed_paths[i]);
            verbose_msg("Path over %d-%d has %lu fragments in it\n", s.left(), s.right(), compressed_paths[i].size());
            new_scaffs.push_back(s);
        }
    }
    //hits = new_scaffs;
	
    verbose_msg("%s\tCompressed overlap graph from %lu to %lu fragments (%f percent)\n",
            bundle_label->c_str(), 
            hits.size(), 
            new_scaffs.size(), 
            (hits.size() - new_scaffs.size())/(double)hits.size());

    hits = new_scaffs;
    sort(hits.begin(), hits.end(), scaff_lt);
    create_overlap_dag(hits, bundle_dag);
}

