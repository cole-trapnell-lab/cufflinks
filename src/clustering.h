#ifndef CLUSTERING_H
#define CLUSTERING_H
/*
 *  abundances.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 4/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// DON'T move this, or mystery compiler errors will result. Affects gcc >= 4.1
#include <boost/graph/vector_as_graph.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>

#ifdef DEBUG
#include <boost/numeric/ublas/io.hpp>
#endif

#include <boost/version.hpp>

#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include <boost/property_map/vector_property_map.hpp>
#endif


#include "abundances.h"

using namespace boost;

typedef adjacency_list <vecS, vecS, undirectedS> AbundanceGraph;

struct ConnectByExonOverlap
{
	void operator()(const AbundanceGroup& cluster,
					AbundanceGraph& G);
};

struct ConnectByAnnotatedGeneId
{
	void operator()(const AbundanceGroup& cluster,
					AbundanceGraph& G);
};

struct ConnectByAnnotatedTssId
{
	void operator()(const AbundanceGroup& cluster,
					AbundanceGraph& G);
};

struct ConnectByAnnotatedProteinId
{
	void operator()(const AbundanceGroup& cluster,
					AbundanceGraph& G);
};

struct ConnectByStrand
{
	void operator()(const AbundanceGroup& cluster,
					AbundanceGraph& G);
};

// A "transcript cluster is a set of transcripts whose projections into the
// genome overlap on the same strand.  They may thus share fragment alignments,
// and so they need to be quantitated together.  After quantitation, they
// can be picked apart.
template<class cluster_policy>
void cluster_transcripts(const AbundanceGroup& transfrags,
						 vector<AbundanceGroup>& transfrags_by_cluster,
						 ublas::matrix<double>* new_gamma = NULL,
                         ublas::matrix<double>* new_iterated_count = NULL,
                         ublas::matrix<double>* new_count = NULL,
                         ublas::matrix<double>* new_fpkm = NULL,
                         ublas::matrix<double>* new_gamma_bootstrap = NULL)
{
	adjacency_list <vecS, vecS, undirectedS> G;
	
	transfrags_by_cluster.clear();
	
	cluster_policy cp;
	
	cp(transfrags, G);
	
	std::vector<int> component(num_vertices(G));
	connected_components(G, &component[0]);
	
	vector<vector<bool> > clusters(transfrags.abundances().size(), 
								   vector<bool>(transfrags.abundances().size(), false));
	
	vector<vector<size_t> > cluster_indices(transfrags.abundances().size());
	for (size_t i = 0; i < transfrags.abundances().size(); ++i)
	{
		clusters[component[i]][i] = true;
		cluster_indices[component[i]].push_back(i);
	}
	for (size_t i = 0; i < cluster_indices.size(); ++i)
	{
		if (cluster_indices[i].empty())
		{
			cluster_indices.resize(i);
			break;
		}
	}
	for (size_t i = 0; i < clusters.size(); ++i)
	{
		AbundanceGroup cluster;
		transfrags.filter_group(clusters[i], cluster);
		if (!cluster.abundances().empty())
			transfrags_by_cluster.push_back(cluster);
	}
	
	if (new_gamma != NULL)
	{
		const ublas::matrix<double>& trans_gamma_cov = transfrags.gamma_cov();
        const ublas::matrix<double>& trans_gamma_bootstrap_cov = transfrags.gamma_bootstrap_cov();
        const ublas::matrix<double>& trans_iterated_count_cov = transfrags.iterated_count_cov();
        const ublas::matrix<double>& trans_count_cov = transfrags.count_cov();
        const ublas::matrix<double>& trans_fpkm_cov = transfrags.fpkm_cov();
        
		ublas::matrix<double>& cov = *new_gamma;
        ublas::matrix<double>& boot_cov = *new_gamma_bootstrap;
        ublas::matrix<double>& iterated_count_cov = *new_iterated_count;
        ublas::matrix<double>& count_cov = *new_count;
        ublas::matrix<double>& fpkm_cov = *new_fpkm;
        
		// number of primary transcripts for this gene
		size_t num_pt = cluster_indices.size();
		cov = ublas::zero_matrix<double>(num_pt, num_pt);
        boot_cov = ublas::zero_matrix<double>(num_pt, num_pt);
        count_cov = ublas::zero_matrix<double>(num_pt, num_pt);
        iterated_count_cov = ublas::zero_matrix<double>(num_pt, num_pt);
        fpkm_cov = ublas::zero_matrix<double>(num_pt, num_pt);
		//cerr << "combined " << combined << endl;
		
		//cerr << "locus isoform gamma cov" << gamma_cov << endl;
		for (size_t L = 0; L < cluster_indices.size(); ++L)
		{
			const vector<size_t>& L_isos = cluster_indices[L];
			for (size_t K = 0; K < cluster_indices.size(); ++K)
			{
				const vector<size_t>& K_isos = cluster_indices[K];
				for (size_t l = 0; l < L_isos.size(); ++l)
				{
					for (size_t k = 0; k < K_isos.size(); ++k)
					{
						cov(L,K) += trans_gamma_cov(L_isos[l],K_isos[k]);
                        boot_cov(L,K) += trans_gamma_bootstrap_cov(L_isos[l],K_isos[k]);
                        count_cov(L,K) += trans_count_cov(L_isos[l],K_isos[k]);
                        iterated_count_cov(L,K) += trans_iterated_count_cov(L_isos[l],K_isos[k]);
                        fpkm_cov(L,K) += trans_fpkm_cov(L_isos[l],K_isos[k]);
					}
				}
			}
		}
	}
}

#endif

