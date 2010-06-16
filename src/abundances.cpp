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

//#define USE_LOG_CACHE

void get_compatibilities(const vector<Scaffold>& transcripts,
						 const vector<MateHit>& alignments,
						 vector<vector<char> >& compatibilities,
						 vector<int>& transcripts_per_alignment)
{
	int M = alignments.size();
	int N = transcripts.size();
	
	vector<Scaffold> alignment_scaffs;
	
	for (size_t i = 0; i < alignments.size(); ++i)
	{
		const MateHit& hit = alignments[i];
		//CuffStrand hs = hit.strand();
		alignment_scaffs.push_back(Scaffold(hit));
	} 
	
	for (int i = 0; i < M; ++i) 
	{
		for (int j = 0; j < N; ++j) 
		{
			if (Scaffold::overlap_in_genome(transcripts[j],alignment_scaffs[i], 0) &&
				Scaffold::compatible(transcripts[j], alignment_scaffs[i]))
			{
				compatibilities[i][j] = 1;
				transcripts_per_alignment[i]++;
			}
		}
	}
}


// This matrix computes the conditional probability that a read came from a 
// given transcripts, indenpendent of abundnace.
void compute_alignment_cond_probs(const vector<Scaffold>& transcripts,
								  const vector<MateHit>& alignments,
								  vector<vector<double> >& cond_probs)
{
	int M = alignments.size();
	int N = transcripts.size();
	// count how many transcripts an alignment is compatible with
	vector<int> transcripts_per_alignment(M,0);
	
	// the M x N matrix of 0/1 incidences
	vector<vector<char> > A(M,vector<char>(N,0));
	
	get_compatibilities(transcripts, alignments, A, transcripts_per_alignment);	
	
	normal norm(inner_dist_mean, inner_dist_std_dev);
	
	//Z[i][j] == Prob(Y = y_i | X = x_j)
	
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			if (A[i][j] == 1) {
				const MateHit& mate = alignments[i];
				pair<int,int> span = mate.genomic_inner_span();
				double inner_dist_prob = 0.0;
				if (span.first != -1 && span.second != -1)
				{
					int inner_dist = transcripts[j].match_length(span.first, 
																 span.second);
					inner_dist_prob = pdf(norm, inner_dist);
				}
				else
				{
					//inner_dist_prob = 1e-3;
					inner_dist_prob = pdf(norm, inner_dist_mean);
				}
				
				int transcript_len = transcripts[j].length();
				double cond_prob = (1.0 / transcript_len);
				cond_prob *= inner_dist_prob; 
				cond_probs[i][j] = cond_prob;
			}
		}
	}
}


// FIXME: use library-adjusted effective_length
double TranscriptAbundance::effective_length() const
{
	if (_transfrag)
	{
		return _transfrag->length();
	}
	else
	{
		return 0;
	}
}

// NOTE: Do not call this function directly - it needs to be called by 
// a containing abundance group, which can probabilistically assign fragment
// counts to this transcript before calling this function.
void TranscriptAbundance::calculate_abundance(const vector<MateHit>& alignments,
											   long double map_mass)
{
	_FPKM = num_fragments();
	if (_FPKM == 0.0 || map_mass == 0.0)
	{
		return;
	}
	_FPKM /= (effective_length() / 1000.0);
	_FPKM /= (map_mass / 1000000.0);
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
	//filtered_group._gamma_covariance = ublas::zero_matrix<double>(num_kept,num_kept);
	//filtered_group._sample_mass = _sample_mass;
	//vector<shared_ptr<Abundance> >& new_ab = filtered_group._abundances;
	//ublas::matrix<double>& new_cov = filtered_group._gamma_covariance;
	
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
	
	filtered_group = AbundanceGroup(new_ab, new_cov, _sample_mass);
	
	//filtered_group.calculate_FPKM_variance();
	//filtered_group.calculate_conf_intervals();
	//filtered_group.calculate_kappas();
}

void AbundanceGroup::get_transfrags(vector<Scaffold>& scaffolds) const
{
	scaffolds.clear();
	foreach(shared_ptr<Abundance> pA, _abundances)
	{
		const Scaffold* pS = pA->transfrag();
		if (pS)
		{
			scaffolds.push_back(*pS);
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

void AbundanceGroup::calculate_counts(const vector<MateHit>& alignments)
{
	int M = alignments.size();
	
	vector<Scaffold> transcripts;
	get_transfrags(transcripts);
	
	if (transcripts.empty())
		return;
	
	int N = transcripts.size();
	
	vector<vector<double> > cond_probs (M,vector<double>(N,0));
	
	// FIXME: should operate on an AbundanceGroup, not on a vector of Scaffolds
	compute_alignment_cond_probs(transcripts, alignments, cond_probs);
	
	vector<double> counts(N, 0);
	
	int total_gene_mass = 0;	
	for (size_t i = 0; i < alignments.size(); ++i)
	{	
		double prob_i = 0.0;
		
		int mate_mass = 0;
		if (alignments[i].left_alignment() || alignments[i].right_alignment())
			mate_mass = 1;
		double assigned_mass = 0;
		
		for (size_t j = 0; j < transcripts.size(); ++j)
		{
			prob_i += _abundances[j]->gamma() * cond_probs[i][j];
		}
		
		if (prob_i)
		{
			total_gene_mass += mate_mass;
			
			for (size_t j = 0; j < transcripts.size(); ++j)
			{
				// each transcript takes a fraction of this reads density proportional 
				// to P(R_i|T_j)/P(R_i)
				
				double mass_frac = mate_mass * (1.0 - alignments[i].error_prob());
				assigned_mass += mate_mass;
				counts[j] += ((_abundances[j]->gamma() *  cond_probs[i][j])/prob_i) * mass_frac;
			}
		}
	}
	
	for (size_t j = 0; j < transcripts.size(); ++j)
	{
		_abundances[j]->num_fragments(counts[j]);
	}
}

void AbundanceGroup::calculate_abundance(const vector<MateHit>& alignments,
										 long double map_mass)
{
	_sample_mass = map_mass;
	calculate_gammas(alignments);	
	calculate_counts(alignments);
	
	// This will compute the transcript level FPKMs
	foreach(shared_ptr<Abundance> pA, _abundances)
	{
		pA->calculate_abundance(alignments, map_mass);
	}
	
	calculate_conf_intervals();
	calculate_kappas();
}

void AbundanceGroup::calculate_conf_intervals()
{
	if (status() == NUMERIC_OK)
	{
		// This will compute the transcript level FPKM confidence intervals
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
			double iso_fpkm_variance = 0.0;
			if (_abundances[j]->effective_length() > 0.0 && _sample_mass > 0)
			{
				iso_fpkm_variance = (1000000000 / (_abundances[j]->effective_length() * _sample_mass));
				iso_fpkm_variance *= iso_fpkm_variance;
				double frac = _gamma_covariance(j,j)*(1 + num_fragments()) + 
				(_abundances[j]->gamma() * _abundances[j]->gamma());
				iso_fpkm_variance *= frac;
				iso_fpkm_variance *= num_fragments();
				double FPKM_hi = _abundances[j]->FPKM() + 2 * sqrt(iso_fpkm_variance);
				double FPKM_lo = max(0.0, _abundances[j]->FPKM() - 2 * sqrt(iso_fpkm_variance));
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				_abundances[j]->FPKM_conf(conf);
				_abundances[j]->FPKM_variance(iso_fpkm_variance);
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
			//TODO: calculate group FPKM variances.
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
		double max_transfrag_FPKM_hi = 0;
		double min_transfrag_FPKM_hi = numeric_limits<double>::max();
		foreach(shared_ptr<Abundance> pA, _abundances)
		{
			double FPKM_hi;
			double FPKM_lo;
			if (pA->effective_length() > 0)
			{
				double fpkm_high = num_fragments();
				fpkm_high /= (pA->effective_length() / 1000.0);
				fpkm_high /= (_sample_mass / 1000000.0);
				
				double var_fpkm = (1000000000 / (_sample_mass * pA->effective_length()));
				var_fpkm *= var_fpkm;
				var_fpkm *= fpkm_high;
				
				FPKM_hi = fpkm_high + 2 * sqrt(var_fpkm);
				FPKM_lo = 0.0;
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				pA->FPKM_conf(conf);
			}
			else
			{
				FPKM_hi = 0.0;
				FPKM_lo = 0.0;
				ConfidenceInterval conf(0.0, 0.0);
				pA->FPKM_conf(conf);
			}
			max_transfrag_FPKM_hi = max(max_transfrag_FPKM_hi, FPKM_hi);
			min_transfrag_FPKM_hi = min(min_transfrag_FPKM_hi, FPKM_hi);
				
		}
		
		// In the case of a numeric failure, the groups error bars need to be 
		// set such that 
		FPKM_conf(ConfidenceInterval(0, max_transfrag_FPKM_hi));
	}
}

void AbundanceGroup::calculate_FPKM_variance()
{
	if (_sample_mass == 0)
	{
		_FPKM_variance = 0.0;
		return;
	}
	
	double variance = 0;
			
	double var_left = 0;
	double var_right = 0;
	double var_gamma = 0;
	for (size_t i = 0; i < _abundances.size(); ++i)
	{
		if (_abundances[i]->effective_length() <= 0.0)
			continue;
		double t = _abundances[i]->effective_length() * _abundances[i]->effective_length();
		t = 1.0/t;
		var_left += (t * _gamma_covariance(i,i));
		
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
			if (_abundances[j]->effective_length() <= 0.0)
				continue;
			double t = _abundances[i]->effective_length() * _abundances[j]->effective_length();
			t = 1.0/t;
			assert (!isinf(t) && !isnan(t));
			var_right += (t * _gamma_covariance(i,j));
		}
		
		var_gamma += (_abundances[i]->gamma() / _abundances[i]->effective_length());
	}
	
	var_gamma *= var_gamma;
	
	variance = (1000000000.0 / _sample_mass);
	variance *= variance;
	double frac = ((1.0 + num_fragments())*(var_left + var_right) + (var_gamma));
	variance *= frac;
	variance *= num_fragments();
	
	assert (!isinf(variance) && !isnan(variance));
	_FPKM_variance = variance;
}

// FIXME: This function doesn't really need to copy the transcripts out of 
// the cluster.  Needs refactoring
bool AbundanceGroup::calculate_gammas(const vector<MateHit>& hits_in_cluster)
{
	vector<Scaffold> transfrags;
	get_transfrags(transfrags);
	
	if (transfrags.empty())
		return true;
	
	vector<int> collapse_counts;
	vector<MateHit> non_redundant_hits_in_gene;
	
	vector<double> gammas;
	
	collapse_hits(hits_in_cluster, non_redundant_hits_in_gene, collapse_counts);
	
#if ASM_VERBOSE
	fprintf(stderr, "%s\tCalculating initial MLE\n", bundle_label->c_str());
#endif	
	
	gamma_mle(transfrags,
			  non_redundant_hits_in_gene,
			  collapse_counts,
			  gammas);
	
#if ASM_VERBOSE
	fprintf(stderr, "%s\tTossing likely garbage isoforms\n", bundle_label->c_str());
#endif
	
	for (size_t i = 0; i < gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			fprintf(stderr, "%s\tWarning: isoform abundance is NaN!\n", bundle_label->c_str());
		}
	}
	
	vector<Scaffold> filtered_scaffolds = transfrags;
	vector<double> filtered_gammas = gammas;
	
	filter_junk_isoforms(filtered_scaffolds, filtered_gammas);
	
	if (filtered_scaffolds.empty())
	{
		//gammas = vector<double>(transfrags.size(), 0.0);
		foreach (shared_ptr<Abundance> ab, _abundances)
		{
			ab->gamma(0);
		}
		_gamma_covariance = ublas::zero_matrix<double>(transfrags.size(), 
													  transfrags.size());
		return true;
	}
	
	filtered_gammas.clear();
	
#if ASM_VERBOSE
	fprintf(stderr, "%s\tRevising MLE\n");
#endif
	
	gamma_mle(filtered_scaffolds,
			  non_redundant_hits_in_gene,
			  collapse_counts,
			  filtered_gammas);
	
	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(filtered_gammas[i]))
		{
			fprintf(stderr, "%s\tWarning: isoform abundance is NaN!\n", bundle_label->c_str());
		}
	}
	
#if ASM_VERBOSE
	fprintf(stderr, "%s\tImportance sampling posterior distribution\n", bundle_label->c_str());
#endif
	
	bool success = gamma_map(filtered_scaffolds,
							 non_redundant_hits_in_gene,
							 collapse_counts,
							 filtered_gammas,
							 _gamma_covariance);
	
	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			fprintf(stderr, "%s\tWarning: isoform abundance is NaN!\n", bundle_label->c_str());
			success = false;
		}
	}
	
	// Now we need to fill in zeros for the isoforms we filtered out of the 
	// MLE/MAP calculation
	vector<double> updated_gammas = vector<double>(transfrags.size(), 0.0);
	ublas::matrix<double> updated_gamma_cov;
	updated_gamma_cov = ublas::zero_matrix<double>(transfrags.size(), 
												   transfrags.size());
	
	size_t curr_filtered_scaff = 0;
	StructurallyEqualScaffolds se;
	vector<size_t> scaff_present(transfrags.size(), transfrags.size());
	
	for (size_t i = 0; i < transfrags.size(); ++i)
	{
		if (curr_filtered_scaff < filtered_scaffolds.size())
		{
			if (se(transfrags[i], filtered_scaffolds[curr_filtered_scaff]))
			{
				scaff_present[i] = curr_filtered_scaff;
				curr_filtered_scaff++;
			}
		}
	}
	
	for (size_t i = 0; i < transfrags.size(); ++i)
	{
		if (scaff_present[i] != transfrags.size())
		{
			// then scaffolds[i] has a non-zero abundance, we need to fill
			// that in along with relevant cells from the covariance matrix
			updated_gammas[i] = filtered_gammas[scaff_present[i]];
			for (size_t j = 0; j < transfrags.size(); ++j)
			{
				if (scaff_present[j] != transfrags.size())
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
		const Scaffold* pS = pA->transfrag();
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
			vector<vector<double> > const & cond_probs,
			const vector<double>& u) {
	// given p, fills U with expected frequencies
	int i,j;
	double ProbY;

	for (i = 0; i < M; ++i) {
		ProbY = 0;
		for (j = 0; j < N; ++j) {
			ProbY += cond_probs[i][j] * p[j];
					//cout << "ProbY = " << ProbY << endl;
		}
		ProbY = ProbY ? (1.0 / ProbY) : 0.0;
		
		for (j = 0; j < N; ++j) 
		{
			U[i][j] = u[i]* cond_probs[i][j] * p[j] * ProbY;
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
			v[j] += U[i][j];
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
//#ifdef DEBUG
//	for (j = 0; j < N; ++j) {
//		cout << p[j] << " ";
//	}
//	cout << endl;
//#endif
	
}
 
//double logLike (int N, 
//				int M, 
//				vector<double> & p,
//				vector<vector<double> > const & cond_prob, 
//				vector<double> const & u) {
//	//int i,j;
//	
//	double ell = 0;
//	double Prob_Y;
//	
//	double* tmp_prob = (float*)calloc(M, sizeof(float));
//	
//	for (int i= 0; i < M; i++) 
//	{
//		//tmp_prob[i] = 0;
//		for (int j= 0; j < N; j++) {
//			tmp_prob[i] += (float)(cond_prob[i][j] * p[j]);
//		}
//	} 
//	
//	for (int i= 0; i < M; i++) 
//	{
//		float l = tmp_prob[i]; 
//		tmp_prob[i] = logf(l); 
//	}
//	
//	for (int i= 0; i < M; i++) 
//	{
//		if (!isinf(tmp_prob[i]) && !isnan(tmp_prob[i]))
//		{
//			ell += tmp_prob[i];
//		}
//	}
//	
//	free(tmp_prob);
//	return ell;
//}

double logLike (int N, 
				int M, 
				vector<double> & p,
				vector<vector<double> > const & cond_prob, 
				vector<double> const & u) {
	int i,j;
	
	double ell = 0;
	double Prob_Y;
	for (i= 0; i < M; i++) {
		Prob_Y = 0;
		for (j= 0; j < N; j++) {
			Prob_Y += cond_prob[i][j] * p[j];
		}
		if (Prob_Y > 0) {
			ell += (u[i] * log(Prob_Y));
		}
	}
	return ell;
}

//#define SEEPROB

double EM (int N, int M, vector<double> & newP, 
		   vector<vector<double> > const & cond_prob, 
		   vector<double> const & u) 
{
	double sum = 0;
	double newEll = 0;
	vector<double> p(N,0);
	vector<vector<double> > U(M,vector<double>(N,0));
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
		
		newEll = logLike(N, M, newP, cond_prob,u);
		
		//fprintf(stderr, "%d\t%lf\n", iter, newEll);
		
		//printf("%.3f %.3f %.3f ", newP[0], newP[1], newP[2]);
		//printf("%.3f %.3f %.3f ", newP[3], newP[4], newP[5]);
		//printf("%.3f %.3f %.3f\n", newP[6], newP[7], newP[8]);
		iter++;
	}
	
	if (iter == max_mle_iterations)
		fprintf(stderr, "%s\tWARNING: ITERMAX reached in abundance estimation, estimation hasn't fully converged\n", bundle_label->c_str());
	
	//fprintf(stderr, "Convergence reached in %d iterations \n", iter);
	return newEll;
}

void compute_fisher(const vector<Scaffold>& transcripts,
					const vector<double>& abundances,
					const vector<MateHit>& alignments,
					const vector<double>& u,
					const vector<vector<char> >& compatibilities,
					const vector<vector<double> >& cond_probs,
					boost::numeric::ublas::matrix<double>& fisher)
{
	int M = alignments.size();
	int N = transcripts.size();  
	
	normal norm(inner_dist_mean, inner_dist_std_dev);
	
	vector<long double> denoms(M, 0.0);
	vector<vector<long double> > P(M,vector<long double>(N,0));

	for (int x = 0; x < M; ++x)
	{
		for (int j = 0; j < N; ++j)
		{
			if (!compatibilities[x][j])
				continue;
			long double alpha = 0.0;
			alpha = cond_probs[x][j];
			alpha *= abundances[j];
			denoms[x] += alpha;
		}
		denoms[x] *= denoms[x];
	}
	
	for (int j = 0; j < N; ++j)
	{
		for (int k = 0; k < N; ++k)
		{
			for (int x = 0; x < M; ++x)
			{
				if (!(compatibilities[x][j] && compatibilities[x][k]))
					continue;
				
				assert(denoms[x] != 0.0);
				
				double fisher_x_j_k = cond_probs[x][j] * cond_probs[x][k] / denoms[x];
				
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
							const vector<double>& collapse_counts,
							double scale,
							const ublas::vector<double>& MLE,
							vector<ublas::vector<double> >& weighted_samples,
							vector<pair<size_t, double> >& sample_weights)
{
	if (cond_probs.empty())
		return;
	
	int N = cond_probs.front().size();
	int M = cond_probs.size();
	
	//cerr << "Cov^-1"<<inv_cov << endl;
	for (size_t i = 0; i < samples.size(); ++i)
	{
		vector<double> sample(samples[i].begin(), samples[i].end()); 
		
		//cerr << "s: "<<samples[i] << endl;
		
		double ell = logLike(N,
								 M,
								 sample, 
								 cond_probs, 
								 collapse_counts);
		
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

bool gamma_map(const vector<Scaffold>& transcripts,
			   const vector<MateHit>& alignments,
			   const vector<int>& collapse_counts,
			   vector<double>& gamma_map_estimate,
			   ublas::matrix<double>& gamma_covariance)
{	
	int N = transcripts.size();	
	int M = alignments.size();
	
	gamma_covariance = ublas::zero_matrix<double>(N);
	
	if (N == 1 || M == 0.0)
	{
		return true;
	}

	typedef ublas::matrix<double> matrix_type;
	matrix_type fisher = ublas::zero_matrix<double>(N,N);
	
	// the M x N matrix of 0/1 incidences
	vector<vector<char> > A(M,vector<char>(N,0));
	vector<int> transcripts_per_alignment(M, 0);
	
	get_compatibilities(transcripts, alignments, A, transcripts_per_alignment);	
	
	vector<double> u;
	for (size_t i = 0; i < collapse_counts.size(); ++i)
	{
		u.push_back(collapse_counts[i]);
	}
	
	vector<vector<double> > cond_probs (M,vector<double>(N,0));
	
	compute_alignment_cond_probs(transcripts, alignments, cond_probs);
	
	compute_fisher(transcripts,
				   gamma_map_estimate,
				   alignments,
				   u,
				   A,
				   cond_probs,
				   fisher);
	
	ublas::matrix<double> epsilon = ublas::zero_matrix<double>(N,N);
	for (int i = 0; i < N; ++i)
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
		fprintf(stderr, "%s\tWarning: Fisher matrix is not positive definite (bad element: %lg)\n", bundle_label->c_str(), ch);
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
		fprintf(stderr, "%s\tWarning: Inverse fisher matrix is not positive definite (bad element: %lg)\n", bundle_label->c_str(), ch);
		return false;
	}
	
	ublas::vector<double> MLE(N);
	for (int i = 0; i < N; ++i)
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
		fprintf(stderr, "%s\tWarning: Covariance matrix is not positive definite (bad element: %lg)\n", bundle_label->c_str(), ch);
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
		
		for (int j = 0; j < N; ++j) {
			if (scaled_sample(j) < 0)
				scaled_sample(j) = 1e-10;
		}
		
		double m = sum(scaled_sample);
		if (m && !isnan(m))
		{
			for (int j = 0; j < N; ++j) {
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
		fprintf(stderr, "%s\tWarning: not-enough samples for MAP re-estimation\n", bundle_label->c_str());
		return false;
	}
	
	double det = determinant(covariance_chol);
	
	double denom = pow(2.0*boost::math::constants::pi<double>(), N/2.0);
	double s = sqrt(det);
	denom *= s;
	
	//assert (det);
	if (s == 0.0)
	{
		fprintf(stderr, "%s\tError: sqrt(det(cov)) == 0, %lf after rounding. \n", bundle_label->c_str(), det);
		//cerr << covariance << endl;
		return false;
	}
	assert (s);
	assert (denom);
	
	//double scale = 1.0 / (denom); 
	
	if (!invertible)
	{
		fprintf(stderr, "%s\tWarning: covariance matrix is not invertible, probability interval is not available\n", bundle_label->c_str());
		return false;
	}
	
	long double log_total_weight = 0.0;
	vector<pair<size_t, double> > sample_weights;
	
	ublas::vector<double> log_expectation(N);
	vector<ublas::vector<double> > weighted_samples;
	
	compute_sample_weights(fisher,
						   cond_probs,
						   samples,
						   u,
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
		fprintf(stderr, "%s\tWarning: restimation failed, importance samples have zero weight.\n\tResorting to MLE and observed Fisher\n", bundle_label->c_str());
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
			fprintf(stderr, "%s\tWarning: isoform abundance is NaN, restimation failed.\n\tResorting to MLE and observed Fisher.", bundle_label->c_str());
			return false;
		}
	}
	
	//fprintf(stderr, "%d valid samples\n", num_valid_samples);
	//cerr << "e = " << log_expectation << " e / w = ";
	//expectation /= total_weight;
	//cerr << expectation << endl;
	
	for (int j = 0; j < N; ++j) 
	{
		if (expectation(j) < 0)
			expectation(j) = 0;
	}
	
	long double m = sum(expectation);
	
	if (m == 0 || isinf(m) || isnan(m))
	{
		fprintf(stderr, "%s\tWarning: restimation failed, could not renormalize MAP estimate\n", bundle_label->c_str());
		return false;
	}
	
	for (int j = 0; j < N; ++j) {
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
	for (int k = 0; k < N; ++k)
	{
		for (int i = 0; i < N; ++i)
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
	for (int k = 0; k < N; ++k)
	{
		for (int i = 0; i < N; ++i)
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

void gamma_mle(const vector<Scaffold>& transcripts,
				   const vector<MateHit>& alignments,
				   const vector<int>& collapse_counts,
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
	
	int M = alignments.size();
	int N = transcripts.size();

	
	if (M > 0)
	{
		
		vector<vector<double> > cond_probs (M,vector<double>(N,0));
		//vector<vector<double> > saliencies (M,vector<double>(N,0));
		
		compute_alignment_cond_probs(transcripts, alignments, cond_probs);
		
		//compute_saliencies(cond_probs, saliencies, saliency_weight);
		
		vector<double> prob(N,0);

		double logL;
		
		vector<int> transcripts_per_alignment(M,0);
		
		// the M x N matrix of 0/1 incidences
		vector<vector<char> > A(M,vector<char>(N,0));
		
		get_compatibilities(transcripts, alignments, A, transcripts_per_alignment);
		
		vector<double> u;
		for (size_t i = 0; i < collapse_counts.size(); ++i)
		{
			double saliency = collapse_counts[i];
			u.push_back(saliency);
		}
		
		logL = EM(N, M, prob, cond_probs, u);

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
				   bool exclude_intra_intron,
                   bool use_non_redundant)
{
	vector<bool> intronic(depth_of_coverage.size(), false);
	depth_of_coverage = vector<int>(depth_of_coverage.size(), 0);
	
	vector<Scaffold> hits;
	
    if (!use_non_redundant)
    {
        for (size_t i = 0; i < scaffolds.size(); ++i)
        {
            const vector<const MateHit*>& scaff_hits = scaffolds[i].mate_hits(); 
            for (vector<const MateHit*>::const_iterator itr = scaff_hits.begin();
                 itr != scaff_hits.end();
                 ++itr)
            {
                hits.push_back(Scaffold(**itr));
            }
        }
    }
    else
    {
        hits = scaffolds;
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
					depth_of_coverage[K - bundle_origin]++;
				}
			}
			else if (op.opcode == CUFF_INTRON)
			{
				pair<map<pair<int,int>,int>::iterator, bool> is = intron_depth_of_coverage.insert(make_pair(make_pair(op.g_left(), op.g_right()), 0));
				is.first->second++;
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
#if ASM_VERBOSE
			//fprintf (stderr, "[%d-%d]: %d might be a major isoform intron\n", itr->first.first, itr->first.second, itr->second);
#endif
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
//#if ASM_VERBOSE
//			if (itr == intron_depth_of_coverage.end())
//			{
//				map<pair<int, int>, int >::const_iterator zi;
//				for (zi = intron_depth_of_coverage.begin();
//					 zi != intron_depth_of_coverage.end();
//					 ++zi)
//				{
//					fprintf(stderr, "intron: [%d-%d], %d\n", zi->first.first, zi->first.second, zi->second);
//				}
//			}
//#endif
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

long double get_map_mass(BundleFactory& bundle_factory)
{
	HitBundle bundle;
	long double map_mass = 0;
	
	while(bundle_factory.next_bundle(bundle))
	{
		const vector<MateHit>& hits = bundle.hits();
		for (size_t i = 0; i < bundle.hits().size(); ++i)
		{
			double mate_len = 0;
			if (hits[i].left_alignment() || hits[i].right_alignment())
				mate_len = 1.0;
			map_mass += mate_len * (1.0 - hits[i].error_prob()); 
		}
	}
	
	bundle_factory.reset();
	return map_mass;
}
