/*
 *  clustering.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/15/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "clustering.h"

void ConnectByExonOverlap::operator()(const AbundanceGroup& cluster,
									  AbundanceGraph& G)
{
	const vector<shared_ptr<Abundance> >& abundances = cluster.abundances();
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		shared_ptr<Scaffold> scaff_i = abundances[i]->transfrag();
		assert (scaff_i);
		
		for (size_t j = i + 1; j < abundances.size(); ++j)
		{
			shared_ptr<Scaffold> scaff_j = abundances[j]->transfrag();
			assert (scaff_j);
			
			if (Scaffold::exons_overlap(*scaff_i, *scaff_j))
				add_edge(i, j, G);
		}
	}
}

void ConnectByAnnotatedGeneId::operator()(const AbundanceGroup& cluster,
										  AbundanceGraph& G)
{
	const vector<shared_ptr<Abundance> >& abundances = cluster.abundances();
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		set<string> i_gene_id = abundances[i]->gene_id();
		for (size_t j = i + 1; j < abundances.size(); ++j)
		{
			set<string> j_gene_id = abundances[j]->gene_id();
			if (i_gene_id == j_gene_id)
			{
				add_edge(i, j, G);
			}
		}
	}
}

void ConnectByAnnotatedTssId::operator()(const AbundanceGroup& cluster,
										  AbundanceGraph& G)
{
	const vector<shared_ptr<Abundance> >& abundances = cluster.abundances();
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		set<string> i_tss_id = abundances[i]->tss_id();
		for (size_t j = i + 1; j < abundances.size(); ++j)
		{
			set<string> j_tss_id = abundances[j]->tss_id();
			if (i_tss_id == j_tss_id)
			{
				add_edge(i, j, G);
			}
		}
	}
}

void ConnectByAnnotatedProteinId::operator()(const AbundanceGroup& cluster,
											 AbundanceGraph& G)
{
	const vector<shared_ptr<Abundance> >& abundances = cluster.abundances();
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		set<string> i_p_id = abundances[i]->protein_id();
		for (size_t j = i + 1; j < abundances.size(); ++j)
		{
			set<string> j_p_id = abundances[j]->protein_id();
			if (i_p_id == j_p_id)
			{
				add_edge(i, j, G);
			}
		}
	}
}

void ConnectByStrand::operator()(const AbundanceGroup& cluster,
								 AbundanceGraph& G)
{
	const vector<shared_ptr<Abundance> >& abundances = cluster.abundances();
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		shared_ptr<Scaffold> scaff_i = abundances[i]->transfrag();
		assert (scaff_i);
		
		for (size_t j = i + 1; j < abundances.size(); ++j)
		{
			shared_ptr<Scaffold> scaff_j = abundances[j]->transfrag();
			assert (scaff_j);
			if (scaff_i->strand() == scaff_j->strand())
			{
				add_edge(i, j, G);
			}
		}
	}
}
