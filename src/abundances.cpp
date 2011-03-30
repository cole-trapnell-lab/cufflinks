/*
 *  abundances.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 4/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 *  NOTE: some of the code in this file was derived from (Eriksson et al, 2008)
 */

#include "abundances.h"
#include <numeric>
#include <limits>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/constants/constants.hpp>

#include "filters.h"
#include "replicates.h"

//#define USE_LOG_CACHE

void compute_compatibilities(vector<shared_ptr<Abundance> >& transcripts,
						 const vector<MateHit>& alignments,
							 vector<vector<char> >& compatibilities)
{
	int M = alignments.size();
	int N = transcripts.size();
	
	vector<Scaffold> alignment_scaffs;
	
	for (size_t i = 0; i < alignments.size(); ++i)
	{
		const MateHit& hit = alignments[i];
		alignment_scaffs.push_back(Scaffold(hit));
	} 
	
    for (int j = 0; j < N; ++j) 
    {
		shared_ptr<Scaffold> transfrag_j = transcripts[j]->transfrag();
		for (int i = 0; i < M; ++i) 
        {
			if (transfrag_j->contains(alignment_scaffs[i])
				&& Scaffold::compatible(*transfrag_j, alignment_scaffs[i]))
            {
				compatibilities[j][i] = 1;
            }
        }
	}
}

AbundanceStatus AbundanceGroup::status() const
{
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->status() != NUMERIC_OK)
		{
			return NUMERIC_FAIL;
		}
	}
	return NUMERIC_OK;
}

double AbundanceGroup::num_fragments() const
{
	double num_f = 0;
	
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		num_f += ab->num_fragments();
	}
	return num_f;
}

double AbundanceGroup::mass_fraction() const
{
	double mass = 0;
	
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		mass += ab->mass_fraction();
	}
	return mass;
}

double AbundanceGroup::mass_variance_fraction() const
{
	double mass_var = 0;
	
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		mass_var += ab->mass_variance_fraction();
	}
	return mass_var;
}

double AbundanceGroup::FPKM() const
{
	double fpkm = 0;
	
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		fpkm += ab->FPKM();
	}
	
	return fpkm;
}

double AbundanceGroup::gamma() const
{
	double gamma = 0;
	
	foreach(shared_ptr<Abundance> ab, _abundances)
	{
		gamma += ab->gamma();
	}
	
	return gamma;
}

void AbundanceGroup::filter_group(const vector<bool>& to_keep, 
								  AbundanceGroup& filtered_group) const
{
	//filtered_group = AbundanceGroup();
	
	assert (to_keep.size() == _abundances.size());
	
	size_t num_kept = 0;
	foreach(bool keeper, to_keep)
	{
		num_kept += keeper;
	}
	
	ublas::matrix<double> new_cov = ublas::zero_matrix<double>(num_kept,num_kept);
	vector<shared_ptr<Abundance> > new_ab;
	
	// rebuild covariance matrix and abundance vector after filtration
	
	size_t next_cov_row = 0;
	for (size_t i = 0; i < _abundances.size(); ++i)
	{
		if (to_keep[i])
		{
			new_ab.push_back(_abundances[i]);
			size_t next_cov_col = 0;
			for (size_t j = 0; j < _abundances.size(); ++j)
			{
				if (to_keep[j])
				{
					new_cov(next_cov_row,next_cov_col) = _gamma_covariance(i, j);
					next_cov_col++;
				}
			}
			next_cov_row++;
		}
	}
	
	filtered_group = AbundanceGroup(new_ab, new_cov);
}

void AbundanceGroup::get_transfrags(vector<shared_ptr<Abundance> >& transfrags) const
{
	transfrags.clear();
	foreach(shared_ptr<Abundance> pA, _abundances)
	{
		shared_ptr<Scaffold> pS = pA->transfrag();
		if (pS)
		{
			transfrags.push_back(pA);
		}
	}
}

set<string> AbundanceGroup::gene_id() const	
{
	set<string> s;
	
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_id();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}

set<string> AbundanceGroup::gene_name() const	
{
	set<string> s;
	
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_name();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}


set<string> AbundanceGroup::tss_id() const	
{
	set<string> s;

	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->tss_id();
		s.insert(sub.begin(), sub.end());
	}

	return s;
}

set<string> AbundanceGroup::protein_id() const	
{
	set<string> s;
	
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->protein_id();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}

const string& AbundanceGroup::locus_tag() const	
{
	static string default_locus_tag = "-";
	const string* pLast = NULL;
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		if (pLast)
		{
			if (pA->locus_tag() != *pLast)
			{
				assert (false);
				return default_locus_tag;
			}
		}
		pLast = &(pA->locus_tag());
	}
	if (pLast)
	{
		return *pLast;
	}
	assert (false);
	return default_locus_tag;
}

const string& AbundanceGroup::reference_tag() const	
{
	static string default_reference_tag = "-";
	const string* pLast = NULL;
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		if (pLast)
		{
			if (pA->reference_tag() != *pLast)
			{
				assert (false);
				return default_reference_tag;
			}
		}
		pLast = &(pA->reference_tag());
	}
	if (pLast)
	{
		return *pLast;
	}
	assert (false);
	return default_reference_tag;
}

double AbundanceGroup::effective_length() const
{
	double eff_len = 0.0;
	double group_fpkm = FPKM();
	if (group_fpkm == 0)
		return 0;
	foreach (shared_ptr<Abundance> ab, _abundances)
	{
		eff_len += (ab->effective_length() * (ab->FPKM() / group_fpkm));
	}
	return eff_len;
}

void AbundanceGroup::calculate_counts(const vector<MateHit>& alignments,
                                      const vector<shared_ptr<Abundance> >& transcripts)
{
	size_t M = alignments.size();
	size_t N = transcripts.size();
	
	if (transcripts.empty())
		return;
    
    map<shared_ptr<ReadGroupProperties const>, double> count_per_replicate;
    
	for (size_t i = 0; i < M; ++i)
	{	
		if (!alignments[i].left_alignment())
			continue;
		
		bool mapped = false;
		for (size_t j = 0; j < N; ++j)
        {
			if (_abundances[j]->cond_probs()->at(i) > 0)
            {
				mapped = true;
				break;
			}
		}
		if (mapped)
        {
            shared_ptr<ReadGroupProperties const> rg_props = alignments[i].read_group_props();
            //assert (parent != NULL);
            pair<map<shared_ptr<ReadGroupProperties const>, double>::iterator, bool> inserted;
            inserted = count_per_replicate.insert(make_pair(rg_props, 0.0));

            double more_mass = alignments[i].collapse_mass();
            inserted.first->second += more_mass;
        }
    }
    
    double avg_X_g = 0.0;
    double avg_mass_fraction = 0.0;
    double avg_variance_fraction = 0.0;
    
    
    for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = count_per_replicate.begin();
         itr != count_per_replicate.end();
         ++itr)
    {
        shared_ptr<ReadGroupProperties const> rg_props = itr->first;
        double scaled_mass = rg_props->scale_mass(itr->second);
        double scaled_total_mass = rg_props->scale_mass(rg_props->normalized_map_mass());
        avg_X_g += scaled_mass;
        shared_ptr<MassDispersionModel const> disperser = rg_props->mass_dispersion_model();
        double scaled_variance = disperser->scale_mass_variance(scaled_mass);
        avg_variance_fraction += scaled_variance / scaled_total_mass;
        assert (scaled_total_mass != 0.0);
        avg_mass_fraction += (scaled_mass / scaled_total_mass);
    }
    
    double num_replicates = count_per_replicate.size();
    
    if (num_replicates)
    {
        avg_X_g /= num_replicates;
        avg_mass_fraction /= num_replicates;
        avg_variance_fraction /= num_replicates;
    }
    
	for (size_t j = 0; j < N; ++j)
	{
		_abundances[j]->num_fragments(_abundances[j]->gamma() * avg_X_g);
        
        double j_avg_mass_fraction = _abundances[j]->gamma() * avg_mass_fraction;
        _abundances[j]->mass_fraction(j_avg_mass_fraction);
        
        double j_avg_variance_fraction = _abundances[j]->gamma() * avg_variance_fraction;
        _abundances[j]->mass_variance_fraction(j_avg_variance_fraction);
        
        if (j_avg_mass_fraction > 0)
        {
            double FPKM = j_avg_mass_fraction * 1000000000/ _abundances[j]->effective_length();
            _abundances[j]->FPKM(FPKM);
        }
        else 
        {
            _abundances[j]->FPKM(0);
            _abundances[j]->mass_variance_fraction(0);
            _abundances[j]->mass_fraction(0);
        }
	}
}

int total_cond_prob_calls = 0;
void collapse_equivalent_hits(const vector<MateHit>& alignments,
                              vector<shared_ptr<Abundance> >& transcripts,
                              vector<shared_ptr<Abundance> >& mapped_transcripts,
                              vector<MateHit>& nr_alignments,
                              vector<double>& log_conv_factors)
{
    int N = transcripts.size();
	int M = alignments.size();
    
    nr_alignments.clear();
    
	vector<vector<char> > compatibilities(N, vector<char>(M,0));
	compute_compatibilities(transcripts, alignments, compatibilities);
    
    vector<vector<double> > cached_cond_probs (M, vector<double>());
    
    vector<bool> replaced(M, false);
    int num_replaced = 0;
    
    vector<BiasCorrectionHelper> bchs;
    for (size_t j = 0; j < N; ++j)
    {
        bchs.push_back(BiasCorrectionHelper(transcripts[j]->transfrag()));   
    }
    
    for(int i = 0 ; i < M; ++i)
    {
        vector<double> cond_probs_i(N,0);
        if (replaced[i] == true)
            continue;
        
        if (cached_cond_probs[i].empty())
        {
            for (int j = 0; j < N; ++j)
            {
                shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
                
                if (compatibilities[j][i]==1)
                {
                    total_cond_prob_calls++;
                    cond_probs_i[j] = bchs[j].get_cond_prob(alignments[i]);
                }
                
            }
            cached_cond_probs[i] = cond_probs_i;
        }
        else
        {
            cond_probs_i = cached_cond_probs[i];
        }
        
        MateHit* curr_align = NULL;
        
        nr_alignments.push_back(alignments[i]);
        curr_align = &nr_alignments.back();
        log_conv_factors.push_back(0);
        
        if (alignments[i].is_multi()) // don't reduce other hits into multihits
            continue;
        
        bool seen_olap = false;
        
        for(int k = i + 1 ; k < M; ++k)
        {
            if (replaced[k] || alignments[k].is_multi())
                continue;
            if (!::overlap_in_genome(curr_align->left(), curr_align->right(),
                                     alignments[k].left(), alignments[k].right()))
            {
                if (seen_olap) 
                    break;
                else
                    continue;
            }
            else
            {
                seen_olap = true;   
            }
            
            vector<double>* cond_probs_k;
            double last_cond_prob = -1;
            
            bool equiv = true;
            
//            for (int j = 0; j < N; ++j)
//            {
//                if (compatibilities[j][k] != compatibilities[j][i])
//                {
//                    equiv = false;
//                    break;
//                }
//            }
//            
//            if (!equiv)
//                continue;
            
            if (cached_cond_probs[k].empty())
            {
                cached_cond_probs[k] = vector<double>(N, 0.0);
                cond_probs_k = &cached_cond_probs[k];
                for (int j = 0; j < N; ++j)
                {
                    shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
                    
                    if (compatibilities[j][k]==1)
                    {
                        total_cond_prob_calls++;
                        (*cond_probs_k)[j] = bchs[j].get_cond_prob(alignments[k]);
                    }
                }
                //cached_cond_probs[k] = cond_probs_k;
            }
            else
            {
                cond_probs_k = &cached_cond_probs[k];
            }
               
            
            for (int j = 0; j < N; ++j)
            {
                if ((*cond_probs_k)[j] != 0 && cond_probs_i[j] != 0)
                {
                    double ratio =  (*cond_probs_k)[j] / cond_probs_i[j];
                    if (last_cond_prob == -1)
                    {
                        last_cond_prob = ratio;
                    }
                    else
                    {
                        if (last_cond_prob != ratio)
                        {
                            equiv = false;
                            break;
                        }
                    }
                }
                else if ((*cond_probs_k)[j] == 0 && cond_probs_i[j] == 0)
                {
                    // just do nothing in this iter.
                    // last_cond_prob = 0.0;
                }
                else
                {
                    equiv = false;
                    break;
                }
            }
            
            // cond_prob_i vector is a scalar multiple of cond_prob_k, so we
            // can collapse k into i via the mass.
            if (equiv && last_cond_prob > 0.0)
            {
                assert (last_cond_prob > 0);
                //double mass_muliplier = sqrt(last_cond_prob);
                double mass_multiplier = log(last_cond_prob);
                assert (!isinf(mass_multiplier) && !isnan(mass_multiplier));
                log_conv_factors[log_conv_factors.size() - 1] += mass_multiplier; 
                replaced[k] = true;
                cached_cond_probs[k].clear();
                vector<double>(cached_cond_probs[k]).swap(cached_cond_probs[k]);
                num_replaced++;
                double more_mass = alignments[k].common_scale_mass() * alignments[k].collapse_mass() ;
                
                curr_align->incr_collapse_mass(more_mass);
            }
        }
    }
    
    N = transcripts.size();
	//M = nr_alignments.size();
        
	for (int j = 0; j < N; ++j) 
    {
		shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
		vector<double>& cond_probs = *(new vector<double>(nr_alignments.size(),0));
		
		BiasCorrectionHelper& bch = bchs[j];
		
        size_t last_cond_prob_idx = 0;
		for(int i = 0 ; i < M; ++i)
		{
            if (!cached_cond_probs[i].empty())
            {
                if (compatibilities[j][i]==1)
                {
                    assert (cached_cond_probs[i].size() > j);
                    cond_probs[last_cond_prob_idx] = cached_cond_probs[i][j];
                }
                last_cond_prob_idx++;
            }
        }
		
        assert (last_cond_prob_idx == nr_alignments.size());
        
		transcripts[j]->effective_length(bch.get_effective_length());
		transcripts[j]->cond_probs(&cond_probs);
		
		if (bch.is_mapped()) 
			mapped_transcripts.push_back(transcripts[j]);
	}
    if (nr_alignments.size())
    {
        asm_verbose("\nReduced %lu frags to %lu (%lf percent)\n", alignments.size(), nr_alignments.size(), 100.0 * nr_alignments.size()/(double)alignments.size());
    }
}

#define PERFORM_EQUIV_COLLAPSE 1

void AbundanceGroup::calculate_abundance(const vector<MateHit>& alignments)
{
	vector<shared_ptr<Abundance> > transcripts;
	get_transfrags(transcripts);
	vector<shared_ptr<Abundance> > mapped_transcripts; // This collects the transcripts that have alignments mapping to them
	
	vector<MateHit> nr_alignments;
	collapse_hits(alignments, nr_alignments);
    
#if PERFORM_EQUIV_COLLAPSE
    vector<MateHit> non_equiv_alignments;
    vector<double> log_conv_factors;
    collapse_equivalent_hits(nr_alignments, transcripts, mapped_transcripts, non_equiv_alignments, log_conv_factors);
    assert (non_equiv_alignments.size() == log_conv_factors.size());
    nr_alignments.clear();
#else
    vector<MateHit> non_equiv_alignments = nr_alignments;
    vector<double> log_conv_factors(nr_alignments.size(), 0);
	compute_cond_probs_and_effective_lengths(non_equiv_alignments, transcripts, mapped_transcripts);
#endif
    
	calculate_gammas(non_equiv_alignments, log_conv_factors, transcripts, mapped_transcripts);		
	
    //non_equiv_alignments.clear();
	//collapse_hits(alignments, nr_alignments);
    // This will also compute the transcript level FPKMs
    calculate_counts(non_equiv_alignments, transcripts);  

	if(corr_multi && !final_est_run)
	{
		update_multi_reads(non_equiv_alignments, mapped_transcripts);
	}
	
	if (final_est_run) // Only on last estimation run
	{
        calculate_conf_intervals();
        calculate_kappas();
    }
    
    //fprintf(stderr, "Total calls to get_cond_prob = %d\n", total_cond_prob_calls);
}

void AbundanceGroup::update_multi_reads(const vector<MateHit>& alignments, vector<shared_ptr<Abundance> > transcripts)
{
	size_t M = alignments.size();
	size_t N = transcripts.size();
	
	if (transcripts.empty())
		return;
    
    for (size_t i = 0; i < M; ++i)
	{
		if (alignments[i].is_multi())
		{
			double expr = 0.0;
			for (size_t j = 0; j < N; ++j)
			{
				expr += _abundances[j]->cond_probs()->at(i) * _abundances[j]->FPKM() * _abundances[j]->effective_length();
			}
			alignments[i].read_group_props()->multi_read_table()->add_expr(alignments[i], expr);
		}
	}
}

void AbundanceGroup::calculate_conf_intervals()
{        
	if (status() == NUMERIC_OK)
	{
		// This will compute the transcript level FPKM confidence intervals
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
			if (_abundances[j]->effective_length() > 0.0 && mass_fraction() > 0)
			{
                double iso_fpkm_var = 0.0;
                
                double norm_frag_density = 1000000000;
                norm_frag_density /= _abundances[j]->effective_length();
                
                norm_frag_density *= mass_variance_fraction();
                
                iso_fpkm_var = norm_frag_density * _abundances[j]->gamma();
                assert (!isnan(_gamma_covariance(j,j)));
                
                iso_fpkm_var += norm_frag_density * norm_frag_density * _gamma_covariance(j,j);
                
				double FPKM_hi = _abundances[j]->FPKM() + 2 * sqrt(iso_fpkm_var);
				double FPKM_lo = max(0.0, _abundances[j]->FPKM() - 2 * sqrt(iso_fpkm_var));
				assert (FPKM_lo <= _abundances[j]->FPKM() && _abundances[j]->FPKM() <= FPKM_hi);
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				_abundances[j]->FPKM_conf(conf);
				_abundances[j]->FPKM_variance(iso_fpkm_var);
			}
			else
			{
				_abundances[j]->FPKM_conf(ConfidenceInterval(0.0, 0.0));
				_abundances[j]->FPKM_variance(0.0);
			}
		}
		
		double group_fpkm = FPKM();
		if (group_fpkm > 0.0)
		{
			calculate_FPKM_variance();
			double FPKM_hi = FPKM() + 2 * sqrt(FPKM_variance());
			double FPKM_lo = max(0.0, FPKM() - 2 * sqrt(FPKM_variance()));
			ConfidenceInterval conf(FPKM_lo, FPKM_hi);
			FPKM_conf(conf);
		}
		else
		{
			_FPKM_variance = 0.0;
			ConfidenceInterval conf(0.0, 0.0);
			FPKM_conf(conf);
		}
	}
	else
	{
		double sum_transfrag_FPKM_hi = 0;
        
		foreach(shared_ptr<Abundance> pA, _abundances)
		{
			double FPKM_hi;
			double FPKM_lo;
			if (pA->effective_length() > 0)
			{
                double norm_frag_density = 1000000000;
                norm_frag_density /= pA->effective_length();
                
                norm_frag_density *= mass_variance_fraction();
                double fpkm_high = norm_frag_density;
                
                double var_fpkm = fpkm_high; 
				
				FPKM_hi = fpkm_high + 2 * sqrt(var_fpkm);
				FPKM_lo = 0.0;
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				assert (FPKM_lo <= pA->FPKM() && pA->FPKM() <= FPKM_hi);
				pA->FPKM_conf(conf);
                pA->FPKM_variance(var_fpkm);
				sum_transfrag_FPKM_hi = max(sum_transfrag_FPKM_hi, FPKM_hi);
			}
			else
			{
				FPKM_hi = 0.0;
				FPKM_lo = 0.0;
				ConfidenceInterval conf(0.0, 0.0);
				pA->FPKM_conf(conf);
                pA->FPKM_variance(0.0);
			}
				
		}
		
		// In the case of a numeric failure, the groups error bars need to be 
		// set such that 
		FPKM_conf(ConfidenceInterval(0.0, sum_transfrag_FPKM_hi));
	}
}


void AbundanceGroup::calculate_FPKM_variance()
{
	if (mass_fraction() == 0)
	{
		_FPKM_variance = 0.0;
		return;
	}
    
    double var_cumul_gamma = 0;
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        for (size_t j = 0; j < _abundances.size(); ++j)
        {
            if (_abundances[i]->effective_length() && _abundances[j]->effective_length())
            {
                assert(!isnan(_gamma_covariance(i,j)));
                double L = _abundances[i]->effective_length() * _abundances[j]->effective_length();
                assert(!isnan(L)); 
                if (L != 0.0)
                {
                    double g = _gamma_covariance(i,j) / L;
                    var_cumul_gamma += g;
                }
            }
        }    
    }
    
    double mass_coeff = (1000000000 * mass_variance_fraction());
	//double left = (mass_coeff * gamma());
    double eff_len = effective_length(); 
    //double right = (mass_coeff * mass_coeff * var_cumul_gamma);
    _FPKM_variance = ((mass_coeff * gamma()) / eff_len) + (mass_coeff * mass_coeff * var_cumul_gamma);
    assert (!isinf(_FPKM_variance) && !isnan(_FPKM_variance));
}


void AbundanceGroup::compute_cond_probs_and_effective_lengths(const vector<MateHit>& alignments,
															  vector<shared_ptr<Abundance> >& transcripts,
															  vector<shared_ptr<Abundance> >& mapped_transcripts)
{																	
	int N = transcripts.size();
	int M = alignments.size();

	vector<vector<char> > compatibilities(N, vector<char>(M,0));
	compute_compatibilities(transcripts, alignments, compatibilities);

	for (int j = 0; j < N; ++j) {
		shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
		vector<double>& cond_probs = *(new vector<double>(M,0));
		
		BiasCorrectionHelper bch(transfrag);
		
		for(int i = 0 ; i < M; ++i)
		{
			if (compatibilities[j][i]==1)
            {
                total_cond_prob_calls++;
				cond_probs[i] = bch.get_cond_prob(alignments[i]);
            }
		}
		
		transcripts[j]->effective_length(bch.get_effective_length());
		transcripts[j]->cond_probs(&cond_probs);
		
		if (bch.is_mapped()) 
			mapped_transcripts.push_back(transcripts[j]);
	}
}

// FIXME: This function doesn't really need to copy the transcripts out of 
// the cluster.  Needs refactoring
bool AbundanceGroup::calculate_gammas(const vector<MateHit>& nr_alignments, 
                                      const vector<double>& log_conv_factors,
									  const vector<shared_ptr<Abundance> >& transcripts, 
									  const vector<shared_ptr<Abundance> >& mapped_transcripts)
{
	if (mapped_transcripts.empty())
    {
		//gammas = vector<double>(transfrags.size(), 0.0);
		foreach (shared_ptr<Abundance> ab, _abundances)
		{
			ab->gamma(0);
		}
		_gamma_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                       transcripts.size());
		return true;
    }
	
	vector<double> gammas;
	
	asm_verbose( "Calculating intial MLE\n");
	
	gamma_mle(mapped_transcripts,
			  nr_alignments,
              log_conv_factors,
			  gammas);
	
	asm_verbose( "Tossing likely garbage isoforms\n");
	
	for (size_t i = 0; i < gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			verbose_msg("Warning: isoform abundance is NaN!\n");
		}
	}
	
	vector<shared_ptr<Abundance> > filtered_transcripts = mapped_transcripts;
	vector<double> filtered_gammas = gammas;
	filter_junk_isoforms(filtered_transcripts, filtered_gammas, mapped_transcripts);
	
	if (filtered_transcripts.empty())
	{
		//gammas = vector<double>(transfrags.size(), 0.0);
		foreach (shared_ptr<Abundance> ab, _abundances)
		{
			ab->gamma(0);
		}
		_gamma_covariance = ublas::zero_matrix<double>(transcripts.size(), 
													  transcripts.size());
		return true;
	}
	
	filtered_gammas.clear();
	
	asm_verbose( "Revising MLE\n");
	
	gamma_mle(filtered_transcripts,
			  nr_alignments,
              log_conv_factors, 
			  filtered_gammas);
	

	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(filtered_gammas[i]))
		{
			verbose_msg("Warning: isoform abundance is NaN!\n");
		}
	}
	
	asm_verbose( "Importance sampling posterior distribution\n");
	
	bool success = true;
	size_t N = transcripts.size();
	
	if (final_est_run) // Only on last estimation run.
	{
		success = gamma_map(filtered_transcripts,
							nr_alignments,
                            log_conv_factors,
							filtered_gammas,
							_gamma_covariance);
	}
	else 
	{
		_gamma_covariance = ublas::zero_matrix<double>(N, N);
	}

	
	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			verbose_msg( "Warning: isoform abundance is NaN!\n");
			success = false;
		}
	}
	
	// Now we need to fill in zeros for the isoforms we filtered out of the 
	// MLE/MAP calculation
	vector<double> updated_gammas = vector<double>(N, 0.0);
	ublas::matrix<double> updated_gamma_cov;
	updated_gamma_cov = ublas::zero_matrix<double>(N, N);
	
	size_t cfs = 0;
	shared_ptr<Scaffold> curr_filtered_scaff = filtered_transcripts[cfs]->transfrag();
	StructurallyEqualScaffolds se;
	vector<size_t> scaff_present(N, N);
	
	for (size_t i = 0; i < N; ++i)
	{
		shared_ptr<Scaffold> scaff_i = transcripts[i]->transfrag();
		if (cfs < filtered_transcripts.size())
		{
			curr_filtered_scaff = filtered_transcripts[cfs]->transfrag();
			if (se(scaff_i, curr_filtered_scaff))
			{
				scaff_present[i] = cfs;
				cfs++;
			}
		}
	}
	
	for (size_t i = 0; i < N; ++i)
	{
		if (scaff_present[i] != N)
		{
			// then scaffolds[i] has a non-zero abundance, we need to fill
			// that in along with relevant cells from the covariance matrix
			updated_gammas[i] = filtered_gammas[scaff_present[i]];
			for (size_t j = 0; j < N; ++j)
			{
				if (scaff_present[j] != N)
				{
					updated_gamma_cov(i,j) = _gamma_covariance(scaff_present[i],
															   scaff_present[j]);
					assert (!isinf(updated_gamma_cov(i,j)));
					assert (!isnan(updated_gamma_cov(i,j)));
				}
			}
		}
	}
	
	// All scaffolds that go in get abundances, but those that get "filtered"
	// from the calculation get zeros.
	//gammas = updated_gammas;
	for (size_t i = 0; i < _abundances.size(); ++i)
	{
		_abundances[i]->gamma(updated_gammas[i]);
		_abundances[i]->status(success ? NUMERIC_OK : NUMERIC_FAIL);
	}
	_gamma_covariance = updated_gamma_cov;
	
	return (status() == NUMERIC_OK);
}

void AbundanceGroup::calculate_kappas()
{
	size_t num_members = _abundances.size();
	_kappa_covariance = ublas::matrix<double>(num_members, 
											  num_members);
	//cerr << gamma_cov <<endl;
	
	assert (_gamma_covariance.size1() == num_members);
	assert (_gamma_covariance.size2() == num_members);
	
	//tss_group.sub_quants = vector<QuantGroup>(isos_in_tss);
	
	double Z_kappa = 0.0;
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		if (pA->effective_length() > 0)
		{
			Z_kappa += pA->gamma() / pA->effective_length();
		}
	}
	
	foreach (shared_ptr<Abundance> pA, _abundances)
	{
		if (pA->effective_length() > 0)
		{
			pA->kappa((pA->gamma() / pA->effective_length()) / Z_kappa);
		}
		else
		{
			pA->kappa(0); 
		}
	}
	
	for (size_t k = 0; k < num_members; ++k)
	{
		for (size_t m = 0; m < num_members; ++m)
		{
			_kappa_covariance(k,m) = _gamma_covariance(k, m);
			double L = _abundances[k]->effective_length() * 
					   _abundances[m]->effective_length();
			if (L > 0.0)
			{
				_kappa_covariance(k,m) /= (L * Z_kappa * Z_kappa);
			}
			else
			{
				_kappa_covariance(k,m) = 0.0;
			}
		}
	}
}

void get_alignments_from_scaffolds(const vector<shared_ptr<Abundance> >& abundances,
								   vector<MateHit>& alignments)
{
	set<const MateHit*> hits_in_gene_set;

	foreach(shared_ptr<Abundance> pA, abundances)
	{			
		shared_ptr<Scaffold> pS = pA->transfrag();
		assert (pS);
		hits_in_gene_set.insert(pS->mate_hits().begin(),
								pS->mate_hits().end());
	}
	
	for(set<const MateHit*>::iterator itr = hits_in_gene_set.begin();
		itr != hits_in_gene_set.end();
		++itr)
	{
		alignments.push_back(**itr);
	}
	
	sort(alignments.begin(), alignments.end(), mate_hit_lt);
}

void round(vector<double> & p) {
	
	double KILLP = 0; // kill all probabilities below this
	
	for (vector<double>::iterator i = p.begin(); i != p.end(); ++i) {
		if ((*i) < KILLP) 
			*i = 0;
	}
}

void Estep (int N, 
			int M, 
			vector<double> const & p,
			vector<vector<double> >& U,
			const vector<vector<double> >& cond_probs,
			const vector<double>& u) {
	// given p, fills U with expected frequencies
	int i,j;	
    
    vector<double> frag_prob_sums(M, 0.0);
    
    for (j = 0; j < N; ++j) 
    {
        for (i = 0; i < M; ++i) 
        {
            frag_prob_sums [i] += cond_probs[j][i] * p[j];
        }
    }
    
    for (i = 0; i < M; ++i) 
    {
        frag_prob_sums[i] = frag_prob_sums[i] ? (1.0 / frag_prob_sums[i]) : 0.0;
    }
    
    for (j = 0; j < N; ++j) 
    {
        for (i = 0; i < M; ++i) 
        {
            double ProbY = frag_prob_sums[i];
            double exp_i_j = u[i] * cond_probs[j][i] * p[j] * ProbY;
            U[j][i] = exp_i_j;
        }
    }
}


void Mstep (int N, int M, vector<double> & p, vector<vector<double> > const & U) {
	vector<double> v(N,0);
	double m = 0;
	int i,j;
	
	//#pragma omp parallel for
	for (j = 0; j < N; ++j) {
		//cout << "." <<  v[j] << ".\n";
		for (i = 0; i < M; ++i) {
			//	cout << U[i][j] << " \n";
			v[j] += U[j][i];
		}
		m += v[j];
	}
	
	if (m)
	{
		for (j = 0; j < N; ++j) {
			p[j] = v[j] / m;
		}
	}
	else
	{
		assert(false);
	}
}
 

double logLike (int N, 
				int M, 
				vector<double> & p,
				const vector<vector<double> >& cond_prob, 
				const vector<double>& u,
                const vector<double>& log_conv_factors) {
	int i,j;
	
	double ell = accumulate(log_conv_factors.begin(), log_conv_factors.end(), 0.0);
	double Prob_Y;
	for (i= 0; i < M; i++) {
		Prob_Y = 0;
		for (j= 0; j < N; j++) {
			Prob_Y += cond_prob[j][i] * p[j];
		}
		if (Prob_Y > 0) {
			ell += (u[i] * log(Prob_Y));
		}
	}
	return ell;
}

void grad_ascent_step (int N, 
                       int M, 
                       vector<double> const & p,
                       vector<vector<double> >& U,
                       const vector<vector<double> >& cond_probs,
                       const vector<double>& u,
                       vector<double>& newP,
                       double& epsilon) 
{
	// given p, fills U with expected frequencies
	//int i,j;	
    
    vector<double> dLL_dj(N, 0.0);
    
    for (size_t i = 0; i < M; ++i)
    {
        double denom = 0.0;
        for (size_t j = 0; j < N; ++j)
        {
            denom += p[j] * cond_probs[j][i];
        }
        
        for (size_t j = 0; j < N; ++j)
        {
            if (denom > 0)
            {
                dLL_dj[j] += u[i] * p[j] * cond_probs[j][i] / denom;
            }
        }
    }
    
    for (size_t j = 0; j < N; ++j)
    {
        newP[j] = p[j] + epsilon * dLL_dj[j];
    }
    
    double m = accumulate(newP.begin(), newP.end(), 0.0);
    if (m > 0)
    {    
        for (int j = 0; j < N; ++j) {
            newP[j] = newP[j] / m;
        }
    }
    else
    {
        return;
    }
}

double grad_ascent (int N, int M, vector<double> & newP, 
                    const vector<vector<double> >& cond_prob, 
                    vector<double> const & u,
                    vector<double> const & log_conv_factors) 
{
    double sum = 0;
	double newEll = 0;
	vector<double> p(N,0);
	vector<vector<double> > U(N, vector<double>(M,0));
	double ell = 0; 
	int iter = 0;
	int j;
	
	for (j = 0; j < N; ++j) {
		p[j] = rand();
		sum += p[j];
	}
	for (j = 0; j < N; ++j) {
		p[j] = p[j] / sum;
	}
	
    ell = logLike(N, M, p, cond_prob, u, log_conv_factors);
    
    double epsilon = 1e-2;
    
    static const double ACCURACY = .001; // convergence criteria
	
	while (iter <= 2 || iter < max_mle_iterations) 
    {
        grad_ascent_step(N, M, p, U, cond_prob, u, newP, epsilon);
		
		newEll = logLike(N, M, newP, cond_prob,u, log_conv_factors);
		
        double delta = newEll - ell;
        //fprintf (stderr, "%g\n", delta);
        if (delta > 0)
        {
            //round(newP);
			p = newP;
			ell = newEll;
            if (abs(delta) < ACCURACY)
            {
                break;
            }
        }
        else
        {
            verbose_msg("Reducing EPSILON \n");
            epsilon /= 10;
        }
		iter++;
	}
	if (iter == max_mle_iterations)
		verbose_msg("Warning: ITERMAX reached in abundance estimation, estimation hasn't fully converged\n");
	verbose_msg("Convergence reached in %d iterations \n", iter);
	return newEll;

}

double EM (int N, int M, vector<double> & newP, 
		   const vector<vector<double> >& cond_prob, 
		   vector<double> const & u,
           vector<double> const & log_conv_factors) 
{
	double sum = 0;
	double newEll = 0;
	vector<double> p(N,0);
	vector<vector<double> > U(N, vector<double>(M,0));
	double ell = 0; 
	int iter = 0;
	int j;
	
	for (j = 0; j < N; ++j) {
		p[j] = rand();
		sum += p[j];
	}
	for (j = 0; j < N; ++j) {
		p[j] = p[j] / sum;
	}
	
	//#ifdef DEBUG
	//	for (j = 0; j < N; ++j) {
	//		cout << p[j] << " ";
	//	}
	//	cout << endl;
	//#endif

	static const double ACCURACY = .000001; // convergence for EM
	
	while (((iter <= 2) || (abs(ell - newEll) > ACCURACY)) && (iter < max_mle_iterations)) {
		if (iter > 0) {
			round(newP);
			p = newP;
			ell = newEll;
		}
		
		Estep(N, M, p, U, cond_prob, u); //  fills U
		Mstep(N, M, newP,U); // fills p
		
		newEll = logLike(N, M, newP, cond_prob,u, log_conv_factors);
		
		//fprintf(stderr, "%d\t%lf\n", iter, newEll);
		
		//printf("%.3f %.3f %.3f ", newP[0], newP[1], newP[2]);
		//printf("%.3f %.3f %.3f ", newP[3], newP[4], newP[5]);
		//printf("%.3f %.3f %.3f\n", newP[6], newP[7], newP[8]);
		iter++;
	}
	if (iter == max_mle_iterations)
		verbose_msg("Warning: ITERMAX reached in abundance estimation, estimation hasn't fully converged\n");
	verbose_msg("Convergence reached in %d iterations \n", iter);
	return newEll;
}

void compute_fisher(const vector<shared_ptr<Abundance> >& transcripts,
					const vector<double>& abundances,
					const vector<MateHit>& alignments,
					const vector<double>& u,
					boost::numeric::ublas::matrix<double>& fisher)
{
	int M = alignments.size();
	int N = transcripts.size();  
	
	vector<long double> denoms(M, 0.0);
	vector<vector<long double> > P(M,vector<long double>(N,0));

	for (int j = 0; j < N; ++j)
	{
		const vector<double>& cond_probs_j = *(transcripts[j]->cond_probs());
		for (int x = 0; x < M; ++x)
		{
			if (cond_probs_j[x]==0)
				continue;
			long double alpha = 0.0;
			alpha = cond_probs_j[x];
			alpha *= abundances[j];
			denoms[x] += alpha;
		}
	}
	
	for (int x = 0; x < M; ++x)
		denoms[x] *= denoms[x];
	
	
	for (int j = 0; j < N; ++j)
	{
		const vector<double>& cond_probs_j = *(transcripts[j]->cond_probs());
		for (int k = 0; k < N; ++k)
		{
			
			const vector<double>& cond_probs_k = *(transcripts[k]->cond_probs());

			for (int x = 0; x < M; ++x)
			{
				if (cond_probs_j[x]==0 && cond_probs_k[x]==0)
					continue;
				
				assert(denoms[x] != 0.0);
				
				double fisher_x_j_k = cond_probs_j[x] * cond_probs_k[x] / denoms[x];
				
				fisher(j,k) += u[x] * fisher_x_j_k;
			}
		}
	}
}

//static long double EPSILON = 1e-10;

// Boost Cholesky factorizations in the spirit of lu.hpp
// Written by Robbie Vogt, found at: 
// http://lists.boost.org/MailArchives/ublas/2005/07/0568.php

namespace boost { namespace numeric { namespace ublas {
	
	// Cholesky factorization
	template<class M>
	double cholesky_factorize (M &m) 
	{
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::value_type value_type;
		
		BOOST_UBLAS_CHECK (m.size1() == m.size2(), external_logic("Cholesky decomposition is only valid for a square, positive definite matrix."));
		
		size_type size = m.size1();
		vector<value_type> d(size);
		//bool positive_definite = true;
		for (size_type i = 0; i < size; ++ i) {
			matrix_row<M> mri (row (m, i));
			for (size_type j = i; j < size; ++ j) {
				matrix_row<M> mrj (row (m, j));
				
				value_type elem = m(i,j) - inner_prod(project(mri,range(0,i)), project(mrj,range(0,i)));
				
				if (i == j) {
					if (elem <= 0.0) {
						// matrix after rounding errors is not positive definite
						return elem;
					}
					else {
						d(i) = sqrtl(elem);
					}
				}
				else {
					m(j,i) = elem / d(i);
				}
			}
		}
		
		// put the diagonal back in
		for (size_type i = 0; i < size; ++ i) {
			m(i,i) = d(i);
		}
		
		//cerr << m << endl;
		for (size_type i = 0; i < size; ++i) {
			for (size_type j = 0; j < i; ++j)
			{
				m(j,i) = 0;
			}
		}
		//cerr << m << endl;
		// decomposition succeeded
		return 0.0;
	}
	
	
	// Cholesky substitution 
	template<class M, class E> 
	void cholesky_substitute (const M &m, vector_expression<E> &e) { 
		typedef const M const_matrix_type; 
		typedef vector<typename E::value_type> vector_type; 
		inplace_solve (m, e, lower_tag ()); 
		inplace_solve (trans(m), e, upper_tag ()); 
	} 
	template<class M, class E> 
	void cholesky_substitute (const M &m, matrix_expression<E> &e) { 
		typedef const M const_matrix_type; 
		typedef matrix<typename E::value_type> matrix_type; 
		inplace_solve (m, e, lower_tag ()); 
		inplace_solve (trans(m), e, upper_tag ()); 
	} 
	template<class E, class M> 
	void cholesky_substitute_left (vector_expression<E> &e, const M &m) { 
		typedef const M const_matrix_type; 
		typedef vector<typename E::value_type> vector_type; 
		inplace_solve (trans(m), e, upper_tag ()); 
		inplace_solve (m, e, lower_tag ()); 
	} 
	template<class E, class M> 
	void cholesky_substitute_left (matrix_expression<E> &e, const M &m) { 
		typedef const M const_matrix_type; 
		typedef matrix<typename E::value_type> matrix_type; 
		inplace_solve (trans(m), e, upper_tag ()); 
		inplace_solve (m, e, lower_tag ()); 
	} 
	// Cholesky matrix inversion 
	template<class M> 
	void cholesky_invert (M &m) 
	{ 
		typedef typename M::size_type size_type; 
		typedef typename M::value_type value_type; 
		size_type size = m.size1(); 
		// determine the inverse of the lower traingular matrix 
		for (size_type i = 0; i < size; ++ i) { 
			m(i,i) = 1 / m(i,i); 
			for (size_type j = i+1; j < size; ++ j) { 
				value_type elem(0); 
				for (size_type k = i; k < j; ++ k) { 
					elem -= m(j,k)*m(k,i); 
				} 
				m(j,i) = elem / m(j,j); 
			} 
		} 
		// multiply the upper and lower inverses together 
		m = prod(trans(triangular_adaptor<M,lower>(m)), triangular_adaptor<M,lower>(m)); 
	} 
	
	
}}}

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool lu_invert_matrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());
	
 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
	if( res != 0 ) return false;
	
 	// create identity matrix of "inverse"
 	inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));
	
 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
	
 	return true;
}

///* Matrix inversion routine.
// Expects input to be PRE-FACTORIZED */
template<class T>
bool chol_invert_matrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
	
 	using namespace boost::numeric::ublas;
	inverse = input;
	
	cholesky_invert(inverse);
	
 	return true;
}


// Adapted for Boost from Numerical Recipes
template<typename ValueType>
class multinormal_generator
{
	typedef boost::mt19937 base_generator_type;
	typedef boost::normal_distribution<> distribution_type;

public:
	// expects the mean vector and the *CHOLESKY* factorization of the covariance
	multinormal_generator(const ublas::vector<ValueType>& mean,
						  const ublas::matrix<ValueType>& chol_cov)
		:	
			_engine(time(0)),
			_distribution(),
			_generator(boost::variate_generator<base_generator_type&, 
			     distribution_type >(_engine, 
									 _distribution))
	{
		_rand = ublas::zero_vector<ValueType>(mean.size());
		_mean = mean; 
		_cholesky = chol_cov;
	}
	
	const ublas::vector<ValueType>& next_rand()
	{
		ublas::vector<ValueType> temp(_mean.size());
		for (size_t i = 0; i < _mean.size(); ++i)
		{
			double r = _generator();
			temp(i) = r;
			_rand(i) = 0.0;
		}
		
		//cerr << "rand ="<<temp << endl;
		//_rand = prod(ublas::triangular_adaptor<ublas::matrix<ValueType>,ublas::lower>(_cholesky), temp);
		for (size_t i = 0; i < _cholesky.size1(); ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
				_rand(i) += _cholesky(i,j) * temp(j);
			}
		}
		//cerr <<_rand << " + " << _mean << "=";
		_rand = _rand + _mean;
		//cerr <<_rand <<endl;
		
		return _rand;
	}
	
private:
	ublas::vector<ValueType>						_rand;
	ublas::vector<ValueType>						_mean;
	ublas::matrix<ValueType>						_cholesky;
	
	base_generator_type								_engine;
	distribution_type								_distribution;
	boost::variate_generator<base_generator_type&, 
							 distribution_type>		_generator;
};

// expects a cholesky factorized covariance matrix
template<class matrix_T>
double determinant(ublas::matrix_expression<matrix_T> const& mat_r)
{
	double det = 1.0;
	
	matrix_T chol(mat_r());
	
	for (size_t i = 0; i < chol.size1(); ++i)
	{
		det *= chol(i,i);
	}
	
	return det * det;
}


// Given log(p) and log(q) returns log(p+q)double
template<class float_type>
float_type log_space_add(float_type log_p, float_type log_q)
{
	if (log_p < log_q)
	{
		float_type tmp = log_p;
		log_p = log_q;
		log_q = tmp;
	}

#if 0
	if (!(log_p >= log_q))
	{
		fprintf(stderr, "ERROR: %lf >= %lf\n", log_p, log_q);
	}
#endif
	assert (log_p >= log_q);
	return log (1.0 + exp(log_q - log_p)) + log_p;
}

void compute_sample_weights(const ublas::matrix<double>& proposed_cov,
							const vector<vector<double> >& cond_probs,
							const vector<ublas::vector<double> >& samples,
							const vector<double>& u,
                            const vector<double>& log_conv_factors,
							double scale,
							const ublas::vector<double>& MLE,
							vector<ublas::vector<double> >& weighted_samples,
							vector<pair<size_t, double> >& sample_weights)
{
	if (cond_probs.empty())
		return;
	
	int M = cond_probs.front().size();
	int N = cond_probs.size();
	
	//cerr << "Cov^-1"<<inv_cov << endl;
	for (size_t i = 0; i < samples.size(); ++i)
	{
		vector<double> sample(samples[i].begin(), samples[i].end()); 
		
		//cerr << "s: "<<samples[i] << endl;
		

		
		double ell = logLike(N,
                             M,
                             sample, 
                             cond_probs, 
							 u,
                             log_conv_factors);
		
		ublas::vector<double> diff = (samples[i] - MLE);
		//cerr << "diff: "<<diff << endl;
		
		ublas::vector<double> diff_transpose = ublas::trans(diff);
		//cerr << "diff^T" << diff_transpose << endl;
		ublas::vector<double> P = prod(proposed_cov, diff);
		//cerr << "Prod: "<< P << endl;
		double X = inner_prod(diff_transpose,P);
		
		//cerr << diff_transpose << " "<< P << " " << X << endl;
		
		double sample_prob = exp(-0.5 * X) / scale;
		
		if (sample_prob == 0.0)
		{
			//			fprintf(stderr, "Error: sample_prob == 0, %lf after rounding. \n", X);
			//			cerr << "diff: "<<diff << endl;//cerr << covariance << endl;
			//			cerr << "Prod: "<< P << endl;
			//			cerr << "s: "<<samples[i] << endl;
			//			return false;
			continue; // prob is zero after rounding, skip this sample
		}
		
		assert (sample_prob);
		assert (!isinf(sample_prob));
		assert (!isnan(sample_prob));
		//cerr << "Prob(sample) = " << sample_prob << endl;
		double log_weight;
		
		
		if (sample_prob == 0)
		{
			continue;
		}
		else
		{
			//assert (sample_prob > 0.0 && sample_prob <= 1.0);
			//sample_prob *= scale;
			double e_p = ell - log(sample_prob);
			log_weight = e_p;
		}
		
		ublas::vector<double> scaled_sample(N);
		for (size_t v = 0; v < scaled_sample.size(); ++v)
		{
			assert (samples[i][v]);
			scaled_sample(v) = log_weight + log(samples[i][v]);
			assert (scaled_sample(v));
			assert (!isinf(scaled_sample(v)));
			assert (!isnan(scaled_sample(v)));
		}
		
		//cerr << scaled_sample << endl;
		weighted_samples.push_back(scaled_sample);
		
		sample_weights.push_back(make_pair(i, log_weight));
	}
}

void compute_posterior_expectation(const vector<ublas::vector<double> >& weighted_samples,
								   const vector<pair<size_t, double> >& sample_weights,
								   ublas::vector<double>& log_expectation,
								   long double& log_total_weight)
{
	for (size_t i = 0; i < weighted_samples.size(); ++i)
	{
		const ublas::vector<double>& scaled_sample = weighted_samples[i];
		double log_weight = sample_weights[i].second;
		if (log_total_weight == 0.0)
		{
			log_expectation = weighted_samples[i];
			log_total_weight = log_weight;
		}
		else
		{			
			for (size_t e = 0; e < log_expectation.size(); ++e)
			{
				log_expectation(e) = log_space_add<long double>(log_expectation[e], scaled_sample[e]);
			}
			log_total_weight = log_space_add<long double>(log_total_weight, log_weight);
		}	
	}
}

bool gamma_map(const vector<shared_ptr<Abundance> >& transcripts,
			   const vector<MateHit>& nr_alignments,
               const vector<double>& log_conv_factors,
			   vector<double>& gamma_map_estimate,
			   ublas::matrix<double>& gamma_covariance)
{	
	size_t N = transcripts.size();	
	size_t M = nr_alignments.size();
	
	gamma_covariance = ublas::zero_matrix<double>(N);
	
	if (N == 1 || M == 0)
	{
		return true;
	}

	typedef ublas::matrix<double> matrix_type;
	matrix_type fisher = ublas::zero_matrix<double>(N,N);
	
	vector<double> u(M);
	for (size_t i = 0; i < M; ++i)
	{
		u[i] = nr_alignments[i].collapse_mass();
	}
	compute_fisher(transcripts,
				   gamma_map_estimate,
				   nr_alignments,
				   u,
				   fisher);
	
	ublas::matrix<double> epsilon = ublas::zero_matrix<double>(N,N);
	for (size_t i = 0; i < N; ++i)
	{
		epsilon(i,i) = 1e-6;
	}
	
	fisher += epsilon; // modify matrix to avoid problems during inverse
	
	ublas::matrix<double> fisher_chol = fisher;
	
	//matrix_type test_fisher(fisher);
	
	//cerr << "FISHER" << fisher << endl << endl;
	
	double ch = cholesky_factorize(fisher_chol);
	if (ch != 0.0)
	{
		verbose_msg("Warning: Fisher matrix is not positive definite (bad element: %lg)\n", ch);
		return false;
	}
	
	//cerr << "FISHER" << fisher << endl << endl;
	
	//cerr << "CHOLESKY" << endl;
	//cerr << test_fisher << endl << endl;
	
	ublas::matrix<double> inv_fisher = ublas::zero_matrix<double>(N,N);
	bool invertible = chol_invert_matrix(fisher_chol, inv_fisher);
	//bool invertible = lu_invert_matrix(fisher, inv_fisher);
	
	//cerr << "FISHER^-1"<< inv_fisher << endl << endl;
	
	//cerr << "FISHER * FISHER^-1" << prod(fisher,inv_fisher) << endl;
	
	//inv_fisher += epsilon;
	gamma_covariance = inv_fisher;
	
	ublas::matrix<double> test_fisher = inv_fisher;
	ch = cholesky_factorize(test_fisher);
	if (ch != 0.0 || !invertible)
	{
		verbose_msg("Warning: Inverse fisher matrix is not positive definite (bad element: %lg)\n", ch);
		return false;
	}
	
	ublas::vector<double> MLE(N);
	for (size_t i = 0; i < N; ++i)
	{
		MLE(i) = gamma_map_estimate[i];
	}
	
	//cerr << "MEAN: "<< MLE << endl;
	
	ublas::matrix<double> covariance = inv_fisher;
	
//	ublas::matrix_vector_range<matrix_type> diag(covariance, ublas::range (0,N), ublas::range (0,N));
//	covariance = 4 * (covariance + (sum(diag) * ublas::identity_matrix<double>(N)));
	
	
	ublas::matrix<double> inv_cov = ublas::zero_matrix<double>(N,N);
	
	ublas::matrix<double> covariance_chol = covariance;
	ch = cholesky_factorize(covariance_chol);
	
	if (ch != 0.0)
	{
		verbose_msg("Warning: Covariance matrix is not positive definite (bad element: %lg)\n", ch);
		return false;
	}
	
	chol_invert_matrix(covariance_chol, inv_cov);
	
	//cerr << "COV" << endl << covariance << endl;
	//cerr << "COV^-1" << inv_cov << endl;
	//cerr << "COV * COV^-1" << prod(covariance,inv_cov) << endl;
	
	multinormal_generator<double> generator(MLE, covariance_chol);
	
	vector<ublas::vector<double> > samples;
	//int num_samples = 1000;
	for (int i = 0; i < num_importance_samples; ++i)
	{
		ublas::vector<double> r = generator.next_rand();
		ublas::vector<double> scaled_sample = r;
		
		for (size_t j = 0; j < N; ++j) {
			if (scaled_sample(j) < 0)
				scaled_sample(j) = 1e-10;
		}
		
		double m = sum(scaled_sample);
		if (m && !isnan(m))
		{
			for (size_t j = 0; j < N; ++j) {
				scaled_sample(j) = scaled_sample(j) / m;
			}
			
			bool has_zero = false;
			for (size_t j = 0; j < scaled_sample.size(); ++j)
			{
				if (scaled_sample[j] == 0)
				{
					has_zero = true;
					break;
				}
			}
			
			if (has_zero)
				continue;
			samples.push_back(scaled_sample);
		}
		else
		{
			//cerr << r << endl;
			//cerr << scaled_sample << endl;
		}
	}
	
	if (samples.size() < 100)
	{
		verbose_msg("Warning: not-enough samples for MAP re-estimation\n");
		return false;
	}
	
	double det = determinant(covariance_chol);
	
	double denom = pow(2.0*boost::math::constants::pi<double>(), N/2.0);
	double s = sqrt(det);
	denom *= s;
	
	//assert (det);
	if (s == 0.0)
	{
		verbose_msg("Error: sqrt(det(cov)) == 0, %lf after rounding. \n", det);
		//cerr << covariance << endl;
		return false;
	}
	assert (s);
	assert (denom);
	
	//double scale = 1.0 / (denom); 
	
	if (!invertible)
	{
		verbose_msg("Warning: covariance matrix is not invertible, probability interval is not available\n");
		return false;
	}
	
	long double log_total_weight = 0.0;
	vector<pair<size_t, double> > sample_weights;
	
	ublas::vector<double> log_expectation(N);
	vector<ublas::vector<double> > weighted_samples;
	
	vector<vector<double> > cond_probs(N, vector<double>());
	for(size_t j = 0; j < N; ++j)
	{
		cond_probs[j]= *(transcripts[j]->cond_probs());
	}
	
	compute_sample_weights(fisher,
						   cond_probs,
						   samples,
						   u,
                           log_conv_factors,
						   denom,
						   MLE,
						   weighted_samples,
						   sample_weights);
	
	compute_posterior_expectation(weighted_samples,
								  sample_weights,
								  log_expectation,
								  log_total_weight);
	
	if (log_total_weight == 0 || sample_weights.size() < 100)
	{
		verbose_msg("Warning: restimation failed, importance samples have zero weight.\n\tResorting to MLE and observed Fisher\n");
		return false;
	}
	
	ublas::vector<long double> expectation(N);
	for (size_t e = 0; e < expectation.size(); ++e)
	{
		expectation(e) = (long double)log_expectation(e) - log_total_weight;
		expectation(e) = exp(expectation(e));
	}
	
	for (size_t e = 0; e < expectation.size(); ++e)
	{
		if (isinf(expectation(e)) || isnan(expectation(e)))
		{
			verbose_msg("Warning: isoform abundance is NaN, restimation failed.\n\tResorting to MLE and observed Fisher.");
			return false;
		}
	}
	
	//fprintf(stderr, "%d valid samples\n", num_valid_samples);
	//cerr << "e = " << log_expectation << " e / w = ";
	//expectation /= total_weight;
	//cerr << expectation << endl;
	
	for (size_t j = 0; j < N; ++j) 
	{
		if (expectation(j) < 0)
			expectation(j) = 0;
	}
	
	long double m = sum(expectation);
	
	if (m == 0 || isinf(m) || isnan(m))
	{
		verbose_msg("Warning: restimation failed, could not renormalize MAP estimate\n\tResorting to MLE and observed Fisher.");
		return false;
	}
	
	for (size_t j = 0; j < N; ++j) {
		expectation(j) = expectation(j) / m;
	}
	
	// revise gamma by setting it to the posterior expectation computed via the
	// importance sampling
	gamma_map_estimate = vector<double>(expectation.begin(), expectation.end());
	
	// calculate the sample - mean vectors, store them in log space
	vector<ublas::vector<double> > sample_expectation_diffs;
	
	for (size_t j = 0; j < weighted_samples.size(); ++j)
	{
		ublas::vector<double> sample = weighted_samples[j];
		double log_sample_weight = sample_weights[j].second;
		for (size_t e = 0; e < expectation.size(); ++e)
		{
			// sample is already log transformed after it was weighted, so we
			// need to divide by the sample weight to recover the original sample
			// value, then undo the log transform, then subtract the mean from it
			sample(e) = (exp(sample(e) - log_sample_weight) - expectation(e));
			//sample(e) *= exp((log_sample_weight - log_total_weight));
		}
		//cerr << sample << endl;
		sample_expectation_diffs.push_back(sample);
	}
	
	// We want to revise the covariance matrix from the samples, since we'll 
	// need it later for the CIs.
	ublas::matrix<double> revised_cov = ublas::zero_matrix<double>(N,N);
	
	// initialize the revised covariance with the values from the first sample
	for (size_t k = 0; k < N; ++k)
	{
		for (size_t i = 0; i < N; ++i)
		{
			double log_sample_weight = sample_weights[0].second;
			double x = sample_expectation_diffs[0](k);
			double z = sample_expectation_diffs[0](i);
			revised_cov(k,i) = x * z;
			//cerr << x * z<<",";
			revised_cov(k,i) *= exp((log_sample_weight - log_total_weight));
			//cerr << exp((log_sample_weight - log_total_weight))<<endl;
		}
	}
	
	// accumulate the contributions from the other samples (doing one cell of 
	// covariance matrix per outer (i x j) loop iteration.
	for (size_t k = 0; k < N; ++k)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 1; j < sample_expectation_diffs.size(); ++j)
			{
				double x = sample_expectation_diffs[j](k);
				double z = sample_expectation_diffs[j](i);
				double log_sample_weight = sample_weights[j].second;
				revised_cov(k,i) += exp((log_sample_weight - log_total_weight)) * x * z;
			}			
		}
	}

	//cerr << "Revised COV" << endl;
	//cerr << revised_cov << endl;
	gamma_covariance = revised_cov;
	
	return true;
}

void gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
			   const vector<MateHit>& nr_alignments,
               const vector<double>& log_conv_factors,
			   vector<double>& gammas)
{
	gammas.clear();
	if (transcripts.empty())
		return;
	
	//long double bundle_mass_fraction = bundle_mass / (long double) map_mass;
	if (transcripts.size() == 1)
	{
		gammas.push_back(1.0);
		return;
	}
	
	size_t M = nr_alignments.size();
	size_t N = transcripts.size();

	
	if (M > 0)
	{
		
		//vector<vector<double> > saliencies (M,vector<double>(N,0));
		
		
		//compute_saliencies(cond_probs, saliencies, saliency_weight);
		
		vector<double> prob(N,0);

		double logL;
		
		
		vector<vector<double> > cond_probs(N, vector<double>());
		for (size_t j = 0; j < N; ++j)
		{
			cond_probs[j] = *(transcripts[j]->cond_probs());
		}
		vector<double> u(M);
		for (size_t i = 0; i < M; ++i)
		{
			u[i] = nr_alignments[i].collapse_mass();
		}
		
		//logL = EM(N, M, prob, cond_probs, u);
        logL = grad_ascent(N, M, prob, cond_probs, u, log_conv_factors);
        
		gammas = prob;
		
		for (size_t i = 0; i < gammas.size(); ++i)
		{
			assert (!isnan(gammas[i]));
		}
	}
	else
	{
		gammas = vector<double>(N, 0.0);
	}
}

void calc_isoform_fpkm_conf_intervals(double FPKM,
									  double variance,
									  ConfidenceInterval& FPKM_conf)
{
	double FPKM_lo = 0.0;
	double FPKM_hi = 0.0;
	FPKM_hi = FPKM + 2 * sqrt(variance);
	FPKM_lo = max(0.0, FPKM - 2 * sqrt(variance));
	FPKM_conf = ConfidenceInterval(FPKM_lo, FPKM_hi);
}

double compute_doc(int bundle_origin, 
				   const vector<Scaffold>& scaffolds,
				   vector<int>& depth_of_coverage,
				   map<pair<int, int>, int>& intron_depth_of_coverage,
				   bool exclude_intra_intron)
{
	vector<bool> intronic(depth_of_coverage.size(), false);
	depth_of_coverage = vector<int>(depth_of_coverage.size(), 0);
		
	set<const MateHit*> hits_in_gene_set;
	
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		hits_in_gene_set.insert(scaffolds[i].mate_hits().begin(),
								scaffolds[i].mate_hits().end());
	}
	
	vector<Scaffold> hits;
	
	for(set<const MateHit*>::iterator itr = hits_in_gene_set.begin();
		itr != hits_in_gene_set.end();
		++itr)
	{
		hits.push_back(Scaffold(**itr));
        hits.back().fpkm((**itr).mass());
	}
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const vector<AugmentedCuffOp>& aug_ops = hits[i].augmented_ops();
		for (size_t j = 0; j < aug_ops.size(); ++j)
		{
			const AugmentedCuffOp& op = aug_ops[j];
			if (op.opcode == CUFF_INTRON)
			{
				for (int K = op.g_left(); K < op.g_right(); ++K)
				{
					intronic[K - bundle_origin] = true; 
				}
			}
		}
	}
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const vector<AugmentedCuffOp>& aug_ops = hits[i].augmented_ops();
		for (size_t j = 0; j < aug_ops.size(); ++j)
		{
			const AugmentedCuffOp& op = aug_ops[j];
			if (op.opcode == CUFF_MATCH)
			{
				for (int K = op.g_left(); K < op.g_right(); ++K)
				{
					depth_of_coverage[K - bundle_origin] += hits[i].fpkm();
				}
			}
			else if (op.opcode == CUFF_INTRON)
			{
				pair<map<pair<int,int>,int>::iterator, bool> is = intron_depth_of_coverage.insert(make_pair(make_pair(op.g_left(), op.g_right()), 0));
				is.first->second += hits[i].fpkm();
			}
		}
	}
	
	vector<int> knockout(depth_of_coverage);
	
	int total_doc = 0;
	int total_len = 0;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const vector<AugmentedCuffOp>& aug_ops = hits[i].augmented_ops();
		for (size_t j = 0; j < aug_ops.size(); ++j)
		{
			const AugmentedCuffOp& op = aug_ops[j];
			if (op.opcode == CUFF_MATCH)
			{
				for (int K = op.g_left(); K < op.g_right(); ++K)
				{
					if (!exclude_intra_intron || !intronic[K - bundle_origin])
					{
						total_doc += knockout[K - bundle_origin];
						total_len += (knockout[K - bundle_origin] != 0);
						knockout[K - bundle_origin] = 0;
					}
				}
			}
		}
	}
	
	return (double)total_doc / (double)total_len;
}

double major_isoform_intron_doc(map<pair<int, int>, int>& intron_doc)
{
	double major_isoform_intron_doc = 0;
	int num_major_introns = 0;
	for(map<pair<int, int>, int>::const_iterator itr = intron_doc.begin();
		itr != intron_doc.end(); 
		++itr)
	{
		bool heaviest = true;
		
		for (map<pair<int,int>, int>::const_iterator itr2 = intron_doc.begin();
			 itr2 != intron_doc.end();
			 ++itr2)
		{	
			if (itr != itr2 &&
				itr->second < itr2->second &&
				overlap_in_genome(itr->first.first,
								  itr->first.second,
								  itr2->first.first,
								  itr2->first.second))
			{
				heaviest = false;
				break;
			}
		}
		
		if (heaviest)
		{
			major_isoform_intron_doc += itr->second;
			num_major_introns++;
		}
	}
	if (num_major_introns)
	{
		return major_isoform_intron_doc / num_major_introns;
	}
	else
	{
		return 0.0;
	}	
}

void record_min_doc_for_scaffolds(int bundle_origin, 
								  const vector<Scaffold>& hits,
								  const vector<int>& depth_of_coverage,
								  const map<pair<int, int>, int>& intron_depth_of_coverage,
								  vector<double>& scaff_doc)
{
	for (size_t h = 0; h < hits.size(); ++h)
	{
		double doc = 99999999.0;
		if (hits[h].has_intron())
			doc = get_intron_doc(hits[h], intron_depth_of_coverage);
		
		doc = min(doc, get_scaffold_min_doc(bundle_origin, 
											hits[h], 
											depth_of_coverage));
		scaff_doc.push_back(doc);
	}
}

void record_doc_for_scaffolds(int bundle_origin, 
							  const vector<Scaffold>& hits,
							  const vector<int>& depth_of_coverage,
							  vector<double>& scaff_doc)
{
	for (size_t h = 0; h < hits.size(); ++h)
	{
		double doc;
		doc = get_scaffold_doc(bundle_origin, 
							   hits[h], 
							   depth_of_coverage);
		scaff_doc.push_back(doc);
	}
}

void record_doc_for_scaffolds(int bundle_origin, 
							  const vector<Scaffold>& hits,
							  const vector<int>& depth_of_coverage,
							  const map<pair<int, int>, int>& intron_depth_of_coverage,
							  vector<double>& scaff_doc)
{
	for (size_t h = 0; h < hits.size(); ++h)
	{
		double doc;
		if (hits[h].has_intron())
			doc = get_intron_doc(hits[h], intron_depth_of_coverage);
		else
			doc = get_scaffold_doc(bundle_origin, 
								   hits[h], 
								   depth_of_coverage);
		scaff_doc.push_back(doc);
	}
}

double get_intron_doc(const Scaffold& s,
					  const map<pair<int, int>, int >& intron_depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	int num_introns = 0;
	int doc = 0;
	for (size_t j = 0; j < aug_ops.size(); ++j)
	{
		const AugmentedCuffOp& op = aug_ops[j];
		if (op.opcode == CUFF_INTRON)
		{
			num_introns++;
			pair<int,int> op_intron(op.g_left(), op.g_right());
			map<pair<int, int>, int >::const_iterator itr = intron_depth_of_coverage.find(op_intron);
			//			assert (itr != intron_depth_of_coverage.end());
			if (itr == intron_depth_of_coverage.end())
			{
				map<pair<int, int>, int >::const_iterator zi;
				for (zi = intron_depth_of_coverage.begin();
					 zi != intron_depth_of_coverage.end();
					 ++zi)
				{
					asm_verbose( "intron: [%d-%d], %d\n", zi->first.first, zi->first.second, zi->second);
				}
			}

			doc += itr->second;
		}
	}	
	return doc / (double)num_introns;
}

double get_scaffold_doc(int bundle_origin, 
						const Scaffold& s,
						const vector<int>& depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	int m_len = 0;
	int doc = 0;
	for (size_t j = 0; j < aug_ops.size(); ++j)
	{
		const AugmentedCuffOp& op = aug_ops[j];
		if (op.opcode == CUFF_MATCH)
		{
			for (int K = op.g_left(); K < op.g_right(); ++K)
			{
				m_len++;
				doc += depth_of_coverage[K - bundle_origin];
			}
		}
	}
	
	return (double) doc / (double) m_len; 
}

double get_scaffold_min_doc(int bundle_origin, 
							const Scaffold& s,
							const vector<int>& depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	int min_doc = 99999999;
	
	for (size_t j = 0; j < aug_ops.size(); ++j)
	{
		const AugmentedCuffOp& op = aug_ops[j];
		if (op.opcode == CUFF_MATCH)
		{
			for (int K = op.g_left(); K < op.g_right(); ++K)
			{
				if (min_doc > depth_of_coverage[K - bundle_origin])
					min_doc = depth_of_coverage[K - bundle_origin];
			}
		}
	}
	
	return min_doc;
}
