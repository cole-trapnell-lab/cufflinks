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

//#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <complex>

#include <boost/math/distributions/chi_squared.hpp>

#include "filters.h"
#include "replicates.h"
#include "sampling.h"
#include "jensen_shannon.h"
#include "rounding.h"


#include "negative_binomial_distribution.h"


#include <Eigen/Dense>
//using Eigen::MatrixXd;

//#ifdef  __USE_ISOC99
///* INFINITY and NAN are defined by the ISO C99 standard */
//#else
//double my_infinity(void) {
//    double zero = 0.0;
//    return 1.0/zero;
//}
//double my_nan(void) {
//    double zero = 0.0;
//    return zero/zero;
//}
//#define INFINITY my_infinity()
//#define NAN my_nan()
//#endif

#define EulerGamma 0.57721566490153286060651209008240243104215933593992

/* The digamma function is the derivative of gammaln.
 
 Reference:
 J Bernardo,
 Psi ( Digamma ) Function,
 Algorithm AS 103,
 Applied Statistics,
 Volume 25, Number 3, pages 315-317, 1976.
 
 From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
 (with modifications for negative numbers and extra precision)
 */
double digamma(double x)
{
    double neginf = -INFINITY;
    static const double c = 12,
    d1 = -EulerGamma,
    d2 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    /*static const double s8 = 691/32760.0, s9 = 1/12.0, s10 = 3617/8160.0;*/
    double result;
#if 0
    static double cache_x = 0;
    static int hits = 0, total = 0;
    total++;
    if(x == cache_x) {
        hits++;
    }
    if(total % 1000 == 1) {
        printf("hits = %d, total = %d, hits/total = %g\n", hits, total,
               ((double)hits)/total);
    }
    cache_x = x;
#endif
    if( x==1.0 )
        return d1;
    
    /* Illegal arguments */
    if((x == neginf) || isnan(x)) {
        return NAN;
    }
    /* Singularities */
    if((x <= 0) && (floor(x) == x)) {
        return neginf;
    }
    /* Negative values */
    /* Use the reflection formula (Jeffrey 11.1.6):
     * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
     *
     * This is related to the identity
     * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
     * where z is the fractional part of x
     * For example:
     * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
     *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
     * Then we use
     * digamma(1-z) - digamma(z) = pi*cot(pi*z)
     */
    if(x < 0) {
        return digamma(1-x) + M_PI/tan(-M_PI*x);
    }
    /* Use Taylor series if argument <= S */
    if(x <= s) return d1 - 1/x + d2*x;
    /* Reduce to digamma(X + N) where (X + N) >= C */
    result = 0;
    while(x < c) {
        result -= 1/x;
        x++;
    }
    /* Use de Moivre's expansion if argument >= C */
    /* This expansion can be computed in Maple via asympt(Psi(x),x) */
    if(x >= c) {
        double r = 1/x, t;
        result += log(x) - 0.5*r;
        r *= r;
#if 0
        result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
#else
        /* this version for lame compilers */
        t = (s5 - r * (s6 - r * s7));
        result -= r * (s3 - r * (s4 - r * t));
#endif
    }
    return result;
}


/* The trigamma function is the derivative of the digamma function.
 
 Reference:
 
 B Schneider,
 Trigamma Function,
 Algorithm AS 121,
 Applied Statistics,
 Volume 27, Number 1, page 97-99, 1978.
 
 From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
 (with modification for negative arguments and extra precision)
 */
double trigamma(double x)
{
    double neginf = -INFINITY,
    small = 1e-4,
    large = 8,
    c = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    c1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,
    b4 = -1./30,
    b6 =  1./42,
    b8 = -1./30,
    b10 = 5./66;
    double result;
    /* Illegal arguments */
    if((x == neginf) || isnan(x)) {
        return NAN;
    }
    /* Singularities */
    if((x <= 0) && (floor(x) == x)) {
        return neginf;
    }
    /* Negative values */
    /* Use the derivative of the digamma reflection formula:
     * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
     */
    if(x < 0) {
        result = M_PI/sin(-M_PI*x);
        return -trigamma(1-x) + result*result;
    }
    /* Use Taylor series if argument <= small */
    if(x <= small) {
        return 1/(x*x) + c + c1*x;
    }
    result = 0;
    /* Reduce to trigamma(x+n) where ( X + N ) >= B */
    while(x < large) {
        result += 1/(x*x);
        x++;
    }
    /* Apply asymptotic formula when X >= B */
    /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
    if(x >= large) {
        double r = 1/(x*x), t;
#if 0
        result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
#else
        t = (b4 + r*(b6 + r*(b8 + r*b10)));
        result += 0.5*r + (1 + r*(b2 + r*t))/x;
#endif
    }
    return result;
}


bool fit_gamma_dist(const vector<double> samples, double& k, double& theta_hat)
{
    
    if (samples.size() == 0)
    {
        k = 0;
        theta_hat = 0;
        return true;
    }
    
    double s_1 = 0.0;
    double s_2 = 0.0;
    double s = 0.0;
    
    double N = samples.size();
    
    BOOST_FOREACH(double sample, samples)
    {
        assert (isnan(sample) == false);
        if (sample > 0)
        {
            s_1 += sample;
            s_2 += log(sample);
        }
    }
    
    if (s_1 > 0)
    {
        s = log(s_1/N) - (s_2/N);
    }
    else
    {
        k = 0;
        theta_hat = 0;
        return true;
    }
    
   
    
    k = (3 - s + sqrt(((s - 3) * (s - 3)) + 24*s)) / (12 * s);
    double k_next = 0;
    
    static const double k_newton_raphson_conv_thresh = 0.001;
    static const int max_iterations = 100;
    int num_iter = 0;
    
    while (abs(k_next - k) < k_newton_raphson_conv_thresh && num_iter < max_iterations)
    {
        if (k > 0)
        {
            k_next = k - ( (log(k) - digamma(k) - s) / ( (1.0/k)  - trigamma(k) ) );
        }
        else
        {
            k = 0;
            theta_hat = 0;
            return true;
        }
        num_iter++;
    }
    
    if (k != 0)
    {
        theta_hat = (1.0/(k*N)) * s_1;
    }
    
    return true;
}

//#define USE_LOG_CACHE

void compute_compatibilities(const vector<shared_ptr<Abundance> >& transcripts,
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

AbundanceGroup::AbundanceGroup(const vector<shared_ptr<Abundance> >& abundances,
                               const ublas::matrix<double>& gamma_covariance,
                               const ublas::matrix<double>& iterated_exp_count_covariance,
                               const ublas::matrix<double>& count_covariance,
                               const ublas::matrix<double>& fpkm_covariance,
                               const long double max_mass_variance,
                               const set<shared_ptr<ReadGroupProperties const> >& rg_props) :
    _abundances(abundances), 
    _iterated_exp_count_covariance(iterated_exp_count_covariance),
    _count_covariance(count_covariance),
    _fpkm_covariance(fpkm_covariance),
    _gamma_covariance(gamma_covariance),
    _max_mass_variance(max_mass_variance),
    _salient_frags(0.0),
    _total_frags(0.0),
    _read_group_props(rg_props)
{
    // Calling calculate_FPKM_covariance() also estimates cross-replicate
    // count variances
    // calculate_FPKM_covariance();
    double fpkm_var = 0.0;
    for (size_t i = 0; i < _fpkm_covariance.size1(); ++i)
    {
        for (size_t j = 0; j < _fpkm_covariance.size2(); ++j)
        {
            assert (!isnan(_fpkm_covariance(i,j)) && !isinf(_fpkm_covariance(i,j)));
            fpkm_var += _fpkm_covariance(i,j);
        }
    }
    
    _FPKM_variance = fpkm_var;
    
    if (FPKM() > 0 && final_est_run && library_type != "transfrags")
    {
        
        ublas::matrix<double> test = _fpkm_covariance;
        double ret = cholesky_factorize(test);
        if (ret != 0 || (_FPKM_variance < 0 && status() == NUMERIC_OK))
        {
            //fprintf(stderr, "Warning: total FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            for (size_t j = 0; j < _abundances.size(); ++j)
            {
                _abundances[j]->status(NUMERIC_FAIL);
            }
        }
        
       if(!(FPKM() == 0 || fpkm_var > 0 || status() != NUMERIC_OK))
       {
           //cerr << _count_covariance << endl;
           //cerr << _fpkm_covariance << endl;
       }
        
        assert (FPKM() == 0 || fpkm_var > 0 || status() != NUMERIC_OK);
    }
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        if (_fpkm_samples.empty())
            _fpkm_samples = vector<double>(_abundances[i]->fpkm_samples().size(), 0);
        for (size_t j = 0; j < _abundances[i]->fpkm_samples().size(); ++j)
        {
            _fpkm_samples[j] += _abundances[i]->fpkm_samples()[j];
        }
    }
    
    fit_gamma_distributions();
    
    calculate_conf_intervals();
    
    if (no_js_tests == false && _read_group_props.size() >= min_reps_for_js_test)
    {
        calculate_kappas();
    }
}

AbundanceStatus AbundanceGroup::status() const
{
    bool has_lowdata_member = false;
    bool has_ok_member = false;
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->status() == NUMERIC_FAIL)
		{
			return NUMERIC_FAIL;
		}
        else if (ab->status() == NUMERIC_LOW_DATA)
		{
			has_lowdata_member = true;
            //return NUMERIC_LOW_DATA;
		}
        else if (ab->status() == NUMERIC_HI_DATA)
		{
			return NUMERIC_HI_DATA;
		}
        else if (ab->status() == NUMERIC_OK)
		{
			has_ok_member = true;
		}
	}
    
    if (has_ok_member == false)
        return NUMERIC_LOW_DATA;
    
	return NUMERIC_OK;
}

void TranscriptAbundance::FPKM_variance(double v)
{ 
    assert (v >= 0); 
    assert(!isnan(v));
    _FPKM_variance = v; 
}

bool AbundanceGroup::has_member_with_status(AbundanceStatus member_status) const
{
    BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->status() == member_status)
		{
			return true;
		}
	}
    return false;
}

double AbundanceGroup::num_fragments() const
{
	double num_f = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		num_f += ab->num_fragments();
	}
    assert (!isnan(num_f));
	return num_f;
}

CountPerReplicateTable AbundanceGroup::num_fragments_by_replicate() const
{
	CountPerReplicateTable cpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (cpr.empty())
        {
            cpr = ab->num_fragments_by_replicate();
        }
        else
        {
            CountPerReplicateTable ab_cpr = ab->num_fragments_by_replicate();
            for (CountPerReplicateTable::const_iterator itr = ab_cpr.begin(); 
                 itr != ab_cpr.end();
                 ++itr)
            {
                CountPerReplicateTable::iterator cpr_itr = cpr.find(itr->first);
                assert (cpr_itr != cpr.end());
                cpr_itr->second += itr->second;
            }
        }
	}
    
    //assert (cpr.empty() != false);
	return cpr;
}

FPKMPerReplicateTable AbundanceGroup::FPKM_by_replicate() const
{
	FPKMPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
        FPKMPerReplicateTable ab_fpr = ab->FPKM_by_replicate();
        
        for (FPKMPerReplicateTable::const_iterator itr = ab_fpr.begin();
             itr != ab_fpr.end();
             ++itr)
        {
            FPKMPerReplicateTable::iterator fpr_itr = fpr.find(itr->first);
            if (fpr_itr != fpr.end())
                fpr_itr->second += itr->second;
            else
                fpr[itr->first] = itr->second;
        }
	}
    
    //assert (cpr.empty() != false);
	return fpr;
}

StatusPerReplicateTable AbundanceGroup::status_by_replicate() const
{
	StatusPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (fpr.empty())
        {
            fpr = ab->status_by_replicate();
        }
        else
        {
            StatusPerReplicateTable ab_fpr = ab->status_by_replicate();
            for (StatusPerReplicateTable::const_iterator itr = ab_fpr.begin(); 
                 itr != ab_fpr.end();
                 ++itr)
            {
                StatusPerReplicateTable::iterator fpr_itr = fpr.find(itr->first);
                assert (fpr_itr != fpr.end());
                
                AbundanceStatus s = itr->second;
                
                if (s == NUMERIC_FAIL)
                {
                    fpr_itr->second = NUMERIC_FAIL;
                }
                else if (s == NUMERIC_LOW_DATA && (fpr_itr->second != NUMERIC_HI_DATA && fpr_itr->second != NUMERIC_FAIL && fpr_itr->second != NUMERIC_OK))
                {
                    fpr_itr->second = NUMERIC_LOW_DATA;
                }
                else if (s == NUMERIC_HI_DATA)
                {
                    fpr_itr->second = NUMERIC_HI_DATA;
                }
                else if (s == NUMERIC_OK && (fpr_itr->second != NUMERIC_HI_DATA && fpr_itr->second != NUMERIC_FAIL))
                {
                    fpr_itr->second = NUMERIC_OK;
                }
            }
        }
	}
    
    //assert (cpr.empty() != false);
	return fpr;
}

double AbundanceGroup::mass_fraction() const
{
	double mass = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		mass += ab->mass_fraction();
	}
	return mass;
}

double AbundanceGroup::mass_variance() const
{
    double mass_var = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		mass_var += ab->mass_variance();
	}
	return mass_var;
}

// This tracks the final modeled variance in the assigned counts.
double AbundanceGroup::num_fragment_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        for (size_t j = 0; j < _abundances.size(); ++j)
        {
            frag_var += _count_covariance(i,j);
        }
    }
    return frag_var;
}

// This tracks the final modeled variance in the assigned counts.
double AbundanceGroup::num_fragment_uncertainty_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        for (size_t j = 0; j < _abundances.size(); ++j)
        {
            frag_var += _iterated_exp_count_covariance(i,j);
        }
    }
    return frag_var;
}

double AbundanceGroup::FPKM() const
{
	double fpkm = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		fpkm += ab->FPKM();
	}
	
	return fpkm;
}

double AbundanceGroup::gamma() const
{
	double gamma = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
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
	BOOST_FOREACH(bool keeper, to_keep)
	{
		num_kept += keeper;
	}
	
	ublas::matrix<double> new_cov = ublas::zero_matrix<double>(num_kept,num_kept);
    ublas::matrix<double> new_iterated_em_count_cov = ublas::zero_matrix<double>(num_kept,num_kept);
    ublas::matrix<double> new_count_cov = ublas::zero_matrix<double>(num_kept,num_kept);
    ublas::matrix<double> new_fpkm_cov = ublas::zero_matrix<double>(num_kept,num_kept);
    
	vector<shared_ptr<Abundance> > new_ab;
    
	vector<vector<double> > new_fpkm_samples(_fpkm_samples.size(), vector<double>(num_kept, 0));
	
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
                    new_iterated_em_count_cov(next_cov_row,next_cov_col) = _iterated_exp_count_covariance(i, j);
                    new_count_cov(next_cov_row,next_cov_col) = _count_covariance(i, j);
                    new_fpkm_cov(next_cov_row,next_cov_col) = _fpkm_covariance(i, j);
					next_cov_col++;
				}
			}
			next_cov_row++;
		}
	}
    
    
//    size_t curr_abundance_idx = 0;
//    for (size_t i = 0; i < _abundances.size(); ++i)
//    {
//        if (to_keep[i])
//		{
//            for (size_t j = 0; j < _fpkm_samples.size(); ++j)
//            {
//                new_fpkm_samples[j][curr_abundance_idx] = _fpkm_samples[j][i];
//            }
//            curr_abundance_idx++;
//        }
//        
//    }
    
   
    
	filtered_group = AbundanceGroup(new_ab, 
                                    new_cov, 
                                    new_iterated_em_count_cov, 
                                    new_count_cov, 
                                    new_fpkm_cov,
                                    _max_mass_variance,
                                    _read_group_props);
    //assert (filtered_group.FPKM() == 0 || new_fpkm_samples.size() > 0);
    
    filtered_group.description(_description);
}

void AbundanceGroup::get_transfrags(vector<shared_ptr<Abundance> >& transfrags) const
{
	transfrags.clear();
	BOOST_FOREACH(shared_ptr<Abundance> pA, _abundances)
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
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_id();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}

set<string> AbundanceGroup::gene_name() const	
{
	set<string> s;
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_name();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}


set<string> AbundanceGroup::tss_id() const	
{
	set<string> s;

	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->tss_id();
		s.insert(sub.begin(), sub.end());
	}

	return s;
}

set<string> AbundanceGroup::protein_id() const	
{
	set<string> s;
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
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
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
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
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
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
	BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
	{
		eff_len += (ab->effective_length() * (ab->FPKM() / group_fpkm));
	}
	return eff_len;
}

//void AbundanceGroup::collect_read_group_props()
//{
//    size_t M = alignments.size();
//    
//	for (size_t i = 0; i < M; ++i)
//	{	
//		if (!alignments[i].left_alignment())
//			continue;
//        shared_ptr<ReadGroupProperties const> rg_props = alignments[i].read_group_props();
//           
//        _read_group_props.insert(rg_props;
//    }
//}

void AbundanceGroup::collect_per_replicate_mass(const vector<MateHit>& alignments,
                                                vector<shared_ptr<Abundance> >& transcripts)
{
    size_t M = alignments.size();
	size_t N = transcripts.size();
	
    //_count_per_replicate.clear();
    
    for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin(); 
         itr != _count_per_replicate.end();
         ++itr)
    {
        itr->second = 0.0;
    }
    
	if (transcripts.empty())
		return;
    
    //map<shared_ptr<ReadGroupProperties const>, double> count_per_replicate;

    vector<shared_ptr<Abundance> > mapped_transcripts; // This collects the transcripts that have alignments mapping to them
	compute_cond_probs_and_effective_lengths(alignments, transcripts, mapped_transcripts);
    
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
            inserted = _count_per_replicate.insert(make_pair(rg_props, 0.0));
            _read_group_props.insert(rg_props);
            
            // these are the *internally* scaled masses.
            double more_mass = alignments[i].collapse_mass();
            inserted.first->second += more_mass;
        }
    }
}

void AbundanceGroup::calculate_locus_scaled_mass_and_variance(const vector<MateHit>& alignments,
                                                              const vector<shared_ptr<Abundance> >& transcripts)
{
	size_t N = transcripts.size();
	
	if (transcripts.empty())
		return;
    
    double avg_X_g = 0.0;
    double avg_mass_fraction = 0.0;
    
    // as long as all the read groups share the same dispersion model (currently true)
    // then all the variances from each read group will be the same, so this
    // averaging step isn't strictly necessary.  Computing it this way is simply
    // convenient.
    vector<double> avg_mass_variances(N, 0.0);
    
    double max_mass_var = 0.0;
    
    double external_scale_factor = -1.0;
    for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin();
         itr != _count_per_replicate.end();
         ++itr)
    {
        shared_ptr<ReadGroupProperties const> rg_props = itr->first;
        
        if (external_scale_factor < 0)
        {
            external_scale_factor = rg_props->external_scale_factor();
        }
        else
        {
            assert (external_scale_factor == rg_props->external_scale_factor());
        }
        
        // Since the _count_per_replicate table stores internally scaled
        // fragment counts, we need to scale the fragment counts up so we 
        // can compare between conditions, rather than just between replicates
        // of this condition.
        double scaled_mass = itr->second;
        double scaled_total_mass = rg_props->normalized_map_mass();
        avg_X_g += scaled_mass;
        shared_ptr<MassDispersionModel const> disperser = rg_props->mass_dispersion_model();
        for (size_t j = 0; j < N; ++j)
        {
            double scaled_variance;
            //scaled_variance = disperser->scale_mass_variance(scaled_mass * _abundances[j]->gamma());
            scaled_variance = _abundances[j]->gamma() * disperser->scale_mass_variance(scaled_mass);
            avg_mass_variances[j] += scaled_variance;
        }
        assert (disperser->scale_mass_variance(scaled_mass) != 0 || scaled_mass == 0); 
        max_mass_var += disperser->scale_mass_variance(scaled_mass);
        assert (scaled_total_mass != 0.0);
        avg_mass_fraction += (scaled_mass / scaled_total_mass);
    }
    
    // Set the maximum mass variance in case we get an identifiability failure
    // and need to bound the group expression.
    if (!_count_per_replicate.empty())
        max_mass_var /= _count_per_replicate.size();
    
    
    double num_replicates = _count_per_replicate.size();
    
    if (num_replicates)
    {
        avg_X_g /= num_replicates;
        avg_mass_fraction /= num_replicates;
        for (size_t j = 0; j < N; ++j)
        {
            avg_mass_variances[j] /= num_replicates;
        }
    }
    
    assert (max_mass_var != 0 || avg_X_g == 0);
    max_mass_variance(max_mass_var);
    
    for (size_t j = 0; j < _abundances.size(); ++j)
	{
		_abundances[j]->num_fragments(_abundances[j]->gamma() * avg_X_g);
        
        double j_avg_mass_fraction = _abundances[j]->gamma() * avg_mass_fraction;
        _abundances[j]->mass_fraction(j_avg_mass_fraction);
        _abundances[j]->mass_variance(avg_mass_variances[j]);
        
        if (j_avg_mass_fraction > 0)
        {
            double FPKM = j_avg_mass_fraction * 1000000000/ _abundances[j]->effective_length();
            FPKM *= 1.0 / external_scale_factor;
            
            _abundances[j]->FPKM(FPKM);
        }
        else 
        {
            _abundances[j]->FPKM(0);
            _abundances[j]->mass_variance(0);
            _abundances[j]->mass_fraction(0);
        }
	}

}

int total_cond_prob_calls = 0;
void collapse_equivalent_hits(const vector<MateHit>& alignments,
                              vector<shared_ptr<Abundance> >& transcripts,
                              vector<shared_ptr<Abundance> >& mapped_transcripts,
                              vector<MateHit>& nr_alignments,
                              vector<double>& log_conv_factors, 
                              bool require_overlap = true)
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
        
        if (corr_multi && alignments[i].is_multi()) // don't reduce other hits into multihits
            continue;
        
        bool seen_olap = false;
        
        for(int k = i + 1 ; k < M; ++k)
        {
            if (replaced[k] || (corr_multi && alignments[k].is_multi()) || alignments[i].read_group_props() != alignments[k].read_group_props())
                continue;
            if (require_overlap && !::overlap_in_genome(curr_align->left(), curr_align->right(),
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
                    double cp_j = (*cond_probs_k)[j];
                    double cp_i = cond_probs_i[j];
                    double ratio =  cp_j / cp_i;
                    if (last_cond_prob == -1)
                    {
                        //assert(ratio < 5);
                        last_cond_prob = ratio;
                    }
                    else
                    {
                        if (last_cond_prob != ratio)
                        //if (abs(last_cond_prob - ratio) > 0.001)
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
            if (equiv)
            {
                if (last_cond_prob > 0.0)
                {
                    //assert(curr_align->read_group_props() == alignments[k].read_group_props());
                    assert (last_cond_prob > 0);
                    //double mass_muliplier = sqrt(last_cond_prob);
                    double mass_multiplier = log(last_cond_prob);
                    //assert(last_cond_prob < 5);
                    assert (!isinf(mass_multiplier) && !isnan(mass_multiplier));
                    log_conv_factors[log_conv_factors.size() - 1] += mass_multiplier; 
                    replaced[k] = true;
                    cached_cond_probs[k].clear();
                    vector<double>(cached_cond_probs[k]).swap(cached_cond_probs[k]);
                    num_replaced++;
                    double more_mass = alignments[k].collapse_mass();
                    curr_align->incr_collapse_mass(more_mass);
                }
                else
                {
                    replaced[k] = true;
                    num_replaced++;
                    cached_cond_probs[k].clear();
                    vector<double>(cached_cond_probs[k]).swap(cached_cond_probs[k]);
                }
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
        verbose_msg("\nReduced %lu frags to %lu (%lf percent)\n", alignments.size(), nr_alignments.size(), 100.0 * (1 - nr_alignments.size()/(double)alignments.size()));
    }
}

void collapse_equivalent_hits_helper(const vector<MateHit>& alignments,
                                     vector<shared_ptr<Abundance> >& transcripts,
                                     vector<shared_ptr<Abundance> >& mapped_transcripts,
                                     vector<MateHit>& nr_alignments,
                                     vector<double>& log_conv_factors)
{
    int N = transcripts.size();
	int M = alignments.size();
    
    if (N == 1)
    {
        nr_alignments = alignments;
        log_conv_factors = vector<double>(M, 0.0);
        return;
    }
    // TODO: Remove this short cut after verifying that it doesn't really make sense
    // for large bundles.  The collapse is almost certainly more efficient.
    // If there's a lot of transcripts, just use the old, overlap constrained 
    // version of the equivalence collapse.
    if (N > 24)
    {
        collapse_equivalent_hits(alignments,
                                 transcripts,
                                 mapped_transcripts,
                                 nr_alignments,
                                 log_conv_factors, 
                                 true);
        return;
    }
    
    vector<vector<const MateHit*> > compat_table(1 << N);
    vector<vector<char> > compatibilities(N, vector<char>(M,0));
    compute_compatibilities(transcripts, alignments, compatibilities);
    
    for(int i = 0; i < M; ++i)
    {
        size_t compat_mask = 0;
        for (int j = 0; j < N; ++j)
        {
            compat_mask |= ((compatibilities[j][i] !=0) << j);
        }
        assert (compat_mask < compat_table.size());
        compat_table[compat_mask].push_back(&(alignments[i]));
    }
    
    for (size_t i = 0; i < compat_table.size(); ++i)
    {
        vector<MateHit> tmp_hits;
        vector<MateHit> tmp_nr_hits;
        vector<double> tmp_log_conv_factors;
        vector<shared_ptr<Abundance> > tmp_mapped_transcripts;
        for (size_t j = 0; j < compat_table[i].size(); ++j)
        {
            tmp_hits.push_back(*(compat_table[i][j]));
        }
        if (tmp_hits.empty())
            continue;
        collapse_equivalent_hits(tmp_hits,
                                 transcripts,
                                 tmp_mapped_transcripts,
                                 tmp_nr_hits,
                                 tmp_log_conv_factors, 
                                 false);
        copy(tmp_nr_hits.begin(), tmp_nr_hits.end(), back_inserter(nr_alignments));
        copy(tmp_log_conv_factors.begin(), tmp_log_conv_factors.end(), back_inserter(log_conv_factors));
    }
}

#define PERFORM_EQUIV_COLLAPSE 1

void AbundanceGroup::calculate_abundance(const vector<MateHit>& alignments, bool perform_collapse, bool calculate_per_replicate, const ublas::matrix<double>* max_count_assign_covariance)
{
	vector<shared_ptr<Abundance> > transcripts;
    
	get_transfrags(transcripts);
	vector<shared_ptr<Abundance> > mapped_transcripts; // This collects the transcripts that have alignments mapping to them
    
	vector<MateHit> nr_alignments;
    
    if (perform_collapse)
    {
        collapse_hits(alignments, nr_alignments);
    }
    else
    {
        nr_alignments = alignments;
    }
    
    collect_per_replicate_mass(nr_alignments, transcripts);
    
    vector<MateHit> non_equiv_alignments;
    vector<double> log_conv_factors;
    if (perform_collapse)
    {
        collapse_equivalent_hits_helper(nr_alignments, transcripts, mapped_transcripts, non_equiv_alignments, log_conv_factors);
        assert (non_equiv_alignments.size() == log_conv_factors.size());
        log_conv_factors = vector<double>(nr_alignments.size(), 0);
        nr_alignments.clear();
        mapped_transcripts.clear();
        compute_cond_probs_and_effective_lengths(non_equiv_alignments, transcripts, mapped_transcripts);
    }
    else
    {
        non_equiv_alignments = nr_alignments;
        log_conv_factors = vector<double>(nr_alignments.size(), 0);
        compute_cond_probs_and_effective_lengths(non_equiv_alignments, transcripts, mapped_transcripts);
    }
    
    vector<double> joint_mle_gammas;
    if (final_est_run || corr_multi || corr_bias) 
    {
        calculate_gammas(non_equiv_alignments, log_conv_factors, transcripts, mapped_transcripts);
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            joint_mle_gammas.push_back(_abundances[i]->gamma());
        }
    }

//    std::map<shared_ptr<ReadGroupProperties const >, ublas::vector<double> > mles_for_read_groups;
//    std::map<shared_ptr<ReadGroupProperties const >, AbundanceStatus > status_per_replicate;
//    std::map<shared_ptr<ReadGroupProperties const >, std::vector<double> > fpkm_samples_for_read_groups;
//    std::map<shared_ptr<ReadGroupProperties const >, double > fpkm_per_replicate;
    
    std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> > ab_group_per_replicate;
    
//	for (std::map<shared_ptr<ReadGroupProperties const >, double >::const_iterator itr =_count_per_replicate.begin(); itr != _count_per_replicate.end(); ++itr)
//	{
//        mles_for_read_groups.insert(make_pair(itr->first, ublas::vector<double>(_abundances.size(), 0)));
//	}
    mapped_transcripts.clear();
    compute_cond_probs_and_effective_lengths(non_equiv_alignments, transcripts, mapped_transcripts);
    
    //non_equiv_alignments.clear();
	//collapse_hits(alignments, nr_alignments);
    //This will also compute the transcript level FPKMs
    
    calculate_iterated_exp_count_covariance(joint_mle_gammas, non_equiv_alignments, transcripts, _iterated_exp_count_covariance);
    
    if (final_est_run || (!corr_multi && !corr_bias))
    {
        if (calculate_per_replicate)
        {
            calculate_per_replicate_abundances(transcripts,
                                               non_equiv_alignments,
                                               log_conv_factors,
                                               ab_group_per_replicate,
                                               &_iterated_exp_count_covariance);
            
            for (size_t i = 0; i < _abundances.size(); ++i)
            {
                CountPerReplicateTable cpr;
                FPKMPerReplicateTable fpr;
                StatusPerReplicateTable spr;
                for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
                     itr != ab_group_per_replicate.end();
                     ++itr)
                {
                    assert(itr->second->abundances().size() == _abundances.size());
                    cpr[itr->first] = itr->second->abundances()[i]->num_fragments();
                    //fprintf(stderr, "FPKM = %lg\n", itr->second->abundances()[i]->FPKM());
                    fpr[itr->first] = itr->second->abundances()[i]->FPKM();
                    spr[itr->first] = itr->second->abundances()[i]->status();
                }
                
                _abundances[i]->num_fragments_by_replicate(cpr);
                _abundances[i]->FPKM_by_replicate(fpr);
                _abundances[i]->status_by_replicate(spr);
            }
        }
    }

    
    // Calculate the initial estimates for the number of fragments originating
    // from each transcript, and set the NB variances
    calculate_locus_scaled_mass_and_variance(non_equiv_alignments, transcripts);  
    
    // Refresh the variances to match the new gammas computed during iterated
    // expectation
    // calculate_locus_scaled_mass_and_variance(non_equiv_alignments, transcripts);  
    
	if(corr_multi && !final_est_run)
	{
		update_multi_reads(non_equiv_alignments, mapped_transcripts);
	}
	
	if (final_est_run) // Only on last estimation run
	{
        // Simulate NB draws and fragment assignment under uncertainty to sample
        // from the BNBs.
        vector<double> frags_per_transcript;
        vector<double> frag_variances;
        
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            frags_per_transcript.push_back(_abundances[i]->num_fragments());
            frag_variances.push_back(_abundances[i]->mass_variance());
        }
        
        ublas::matrix<double> count_assign_covariance = _iterated_exp_count_covariance;
        if (calculate_per_replicate == false)
        {
            //cerr << "+++++++++++++++" << endl;
            
            //cerr << "pooled:" << endl;
            //cerr << _iterated_exp_count_covariance << endl;
            
            //cerr << "before:" << endl;
            //cerr << count_assign_covariance << endl;
            
            if (max_count_assign_covariance)
            {
//                for (size_t i = 0; i < _abundances.size(); ++i)
//                {
//                    count_assign_covariance(i,i) = std::max<double>(count_assign_covariance(i,i), (*max_count_assign_covariance)(i,i));
//                }
                count_assign_covariance = *max_count_assign_covariance;
            }
            //cerr << "after:" << endl;
            //cerr << count_assign_covariance << endl;
        }
        simulate_count_covariance(frags_per_transcript, frag_variances, count_assign_covariance, non_equiv_alignments, transcripts, _count_covariance, _assigned_count_samples);
    
        // Calling calculate_FPKM_covariance() also estimates cross-replicate
        // count variances
        
        
        if (calculate_per_replicate == false)
        {
            generate_fpkm_samples();
            fit_gamma_distributions();
        }
        else
        {

            ublas::vector<double> joint_mle = ublas::zero_vector<double>(_abundances.size());
            for (size_t i = 0; i < joint_mle.size(); ++i)
            {
                joint_mle(i) = _abundances[i]->gamma() * num_fragments();
            }
            
            //cerr << _iterated_exp_count_covariance << endl;
            
            ublas::matrix<double> chol_covariance = _iterated_exp_count_covariance;
            
            ublas::matrix<double> epsilon = ublas::zero_matrix<double>(_abundances.size(),_abundances.size());
            for (size_t i = 0; i < _abundances.size(); ++i)
            {
                epsilon(i,i) = 1e-6;
            }
            
            chol_covariance  += epsilon; // modify matrix to avoid problems during inverse
            
            
            double ret = cholesky_factorize(chol_covariance);
            if (ret != 0)
            {
                // set status fail?
            }
            
            ublas::matrix<double> inverse_covariance = _iterated_exp_count_covariance;
            bool invertible = chol_invert_matrix(chol_covariance, inverse_covariance);
            if (invertible == false)
            {
                // set status fail?
            }
            
            //cerr << inverse_covariance << endl;
            
            // calculate the mahalanobis distance between each replicate and the joint mle to flag outliers
            // TODO: should re-do this to work on the kappas (or the gammas)
            
            std::map<shared_ptr<ReadGroupProperties const >, ublas::vector<double> > mle_per_replicate;
            for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
                 itr != ab_group_per_replicate.end();
                 ++itr)
            {
                ublas::vector<double> mle = ublas::zero_vector<double>(_abundances.size());
                for (size_t i = 0; i < mle.size(); ++i)
                {
                    mle(i) = itr->second->abundances()[i]->gamma() * num_fragments();
                }
                
                ublas::vector<double> diff = (mle - joint_mle);
                ublas::vector<double> p_diff = prod(diff, inverse_covariance);
                double mahalonobis_dist = inner_prod(p_diff, diff);
                
                boost::math::chi_squared_distribution<double> csd(_abundances.size());
                double tail_1 = cdf(csd, mahalonobis_dist);
                double p_val = 1.0 - (tail_1);
                
                if (p_val < min_outlier_p)
                {
                    for (size_t i = 0; i < _abundances.size(); ++i)
                    {
                        StatusPerReplicateTable st = _abundances[i]->status_by_replicate();
                        st[itr->first] = NUMERIC_LOW_DATA;
                        _abundances[i]->status_by_replicate(st);
                    }
                }
            }
            
            for (size_t i = 0; i < _abundances.size(); ++i)
            {
                vector<double> ab_i_fpkm_samples;
                
                for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
                     itr != ab_group_per_replicate.end();
                     ++itr)
                {
                    StatusPerReplicateTable st = _abundances[i]->status_by_replicate();
                    StatusPerReplicateTable::const_iterator si = st.find(itr->first);
                    if (si == st.end() || si->second == NUMERIC_LOW_DATA)
                        continue;
                    ab_i_fpkm_samples.insert(ab_i_fpkm_samples.end(), itr->second->abundances()[i]->fpkm_samples().begin(), itr->second->abundances()[i]->fpkm_samples().end());
                }
                _abundances[i]->fpkm_samples(ab_i_fpkm_samples);
            }
            
            vector<double> ab_group_fpkm_samples;
            for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
                 itr != ab_group_per_replicate.end();
                 ++itr)
            {
                ab_group_fpkm_samples.insert(ab_group_fpkm_samples.end(), itr->second->fpkm_samples().begin(), itr->second->fpkm_samples().end());
            }
            
            fpkm_samples(ab_group_fpkm_samples);

            calculate_FPKM_covariance();
            
            fit_gamma_distributions();
            
            // Derive confidence intervals from the FPKM variance/covariance matrix
            calculate_conf_intervals();
            
            // Calculate the inter-group relative abundances and variances
            if (no_js_tests == false && _read_group_props.size() >= min_reps_for_js_test)
            {
                calculate_kappas();
            }
        }
    }
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        for (size_t j = 0; j < _abundances.size(); ++j)
        {
            if (i != j)
            {
                assert(!isinf(_fpkm_covariance(i,j)) && !isnan(_fpkm_covariance(i,j)));
                if (_abundances[i]->transfrag()->contains(*_abundances[j]->transfrag()) &&
                    Scaffold::compatible(*_abundances[i]->transfrag(),*_abundances[j]->transfrag()))
                {
                    _abundances[j]->status(NUMERIC_LOW_DATA);
                }
            }
        }
    }
    
    assert (FPKM() == 0 || _assigned_count_samples.size() > 0);
    
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


long double solve_beta(long double A, long double B, long double C)
{
    long double a = -C/B;
    long double b = (A + 4*A*C/(B*B) - (4*C/B));
    long double c = -A + B - 5*A*A*C/(B*B*B) + 10*A*C/(B*B) - 5*C/B;
    long double d = 2*A*A*A*C/(B*B*B*B) - 6*A*A*C/(B*B*B) + 6*A*C/(B*B) - 2*C/B;
    complex<long double> q((3*a*c - b*b)/(a*a*9.0));
    complex<long double> r((9.0*a*c*b - 27.0*a*a*d - 2.0*b*b*b)/(a*a*a*54.0));    
    complex<long double> s1 = std::pow((r + std::sqrt(q*q*q + r*r)),complex<long double>(1/3.0));
    complex<long double> s2 = std::pow((r - std::sqrt(q*q*q + r*r)),complex<long double>(1/3.0));
    complex<long double> R1 = s1 + s2 - complex<long double>(b/(a*3.0));   
    complex<long double> R2 = -(s1+s2)/complex<long double>(2.0) - complex<long double>(b/(a*3.0)) + (s1-s2) * complex<long double>(0, sqrtl(3.0)/2.0);  
    complex<long double> R3 = -(s1+s2)/complex<long double>(2.0) - complex<long double>(b/(a*3.0)) - (s1-s2) * complex<long double>(0, sqrtl(3.0)/2.0);
  
    vector<long double> roots;
    if (R1.imag() == 0)
        roots.push_back(R1.real());
    if (R2.imag() == 0)
        roots.push_back(R2.real());
    if (R3.imag() == 0)
        roots.push_back(R3.real());
    sort(roots.begin(), roots.end());

    if (roots.empty())
        return 0;
    
    long double root = roots.back();
    return root;
}


// This function takes the point estimate of the number of fragments from
// a transcript, the iterated expection count matrix, and the locus level
// cross replicate variance, and calculates the transcript-level cross-replicate
// count variance
bool estimate_count_variance(long double& variance,
                             double gamma_t, 
                             double psi_t_count_var, 
                             double X_g, 
                             double V_X_g_t,
                             double l_t)
{
    if (l_t == 0)
    {
        return true;
    }

    long double A = X_g * gamma_t;
    
    long double B = V_X_g_t;
    
    long double C = psi_t_count_var;
    
    variance = 0.0;
    bool numeric_ok = true;
    
    long double dispersion = V_X_g_t - (X_g * gamma_t);
    
    if (psi_t_count_var < 0)
    {
        //fprintf (stderr, "Warning: psi_t is negative! (psi_t = %lf)\n", psi_t);
        psi_t_count_var = 0;
    }
    assert (psi_t_count_var >= 0);
    
    // we multiply A with the constants here to make things work out 
    // at the end of the routine when we multiply by the square of those
    // constants
    long double poisson_variance = A + psi_t_count_var;
    long double alpha = 0.0;
    long double beta = 0.0;
    long double bnb_mean = 0.0;
    long double r = 0.0;
    
    if (dispersion < -1 || abs(dispersion) < 1)
    {
        // default to poisson dispersion
        variance = poisson_variance;
    }
    else // there's some detectable overdispersion here, use mixture of negative binomials
    {
        if (psi_t_count_var < 1) 
        {
            // default to regular negative binomial case.
            variance = V_X_g_t;
        }
        else
        {
            r = ceil((A * A) / (B - A));
            
            if (r < 0)
            {
                numeric_ok = false;
            }

            // exact cubic
            beta = solve_beta(A,B,C);
            alpha = 1.0 - (A/(A-B)) * beta;
            //fprintf(stdout, "%Lg\t%Lg\t%Lg\n",A,alpha,beta);
            //if (beta <= 2 || alpha <= 1)
            if (alpha <= 2)
            {
                //printf ("Warning: beta for is %Lg\n", beta);
                numeric_ok = false;
                variance = V_X_g_t;
            }
            else
            {                
                bnb_mean = r * beta / (alpha - 1.0);
                variance = r * (alpha + r - 1.0) * beta * (alpha + beta - 1);
                variance /= (alpha - 2.0) * (alpha - 1.0) * (alpha - 1.0);
            }
            if (variance < 0)
            {
                numeric_ok = false;
                variance = V_X_g_t;
            }
            
            if (variance == 0 && A != 0)
            {
                variance = poisson_variance;
            }
            
            assert (!numeric_ok || variance >= poisson_variance);
            assert (!numeric_ok || variance >= V_X_g_t);
            
            if (variance < poisson_variance)
                variance = poisson_variance;
            
            if (variance < V_X_g_t)
                variance = V_X_g_t;
            
            //assert (abs(FPKM - mean) < 1e-3);
        }
    }
    
    if (variance < 0)
        variance = 0;
    
    variance = ceil(variance);
    
    assert (!numeric_ok || (!isinf(variance) && !isnan(variance)));
    assert (!numeric_ok || variance != 0 || A == 0);
    return numeric_ok;
}

// This version of the function draws directly from the iterated expectation count covariance matrix
// and treats the betas as normally distributed (which is a good approximation), allowing good capture
// of their covariance without resorting to lots of expensive fragment level sampling.
bool simulate_count_covariance(const vector<double>& num_fragments,
                               const vector<double>& frag_variances,
                               const ublas::matrix<double>& iterated_exp_count_covariance,
                               const vector<MateHit>& nr_alignments,
                               const vector<shared_ptr<Abundance> >& transcripts,
                               ublas::matrix<double>& count_covariance,
                               vector<Eigen::VectorXd>& assigned_count_samples)
{
    count_covariance = ublas::zero_matrix<double>(transcripts.size(), transcripts.size());
    
    if (frag_variances.size() <= 1)
    {
        count_covariance(0,0) = frag_variances[0];
        //return;
    }
    
    double total_frag_counts = accumulate(num_fragments.begin(), num_fragments.end(), 0.0);
    double total_var = accumulate(frag_variances.begin(), frag_variances.end(), 0.0);
    
    if (total_frag_counts == 0)
    {
        count_covariance = ublas::zero_matrix<double>(transcripts.size(), transcripts.size());
        assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * num_frag_assignments, Eigen::VectorXd::Zero(transcripts.size()));
        return true;
    }
    
    //size_t num_frag_count_draws = 1000;
    //const int num_multinomial_samples = 1;
    
    boost::mt19937 rng;
    
    
    vector<Eigen::VectorXd > generated_counts (num_frag_count_draws, Eigen::VectorXd::Zero(transcripts.size()));
    ublas::vector<double> mle_frag_counts = ublas::zero_vector<double>(transcripts.size());
    
    for (size_t j = 0; j < transcripts.size(); ++j)
    {
        mle_frag_counts(j) = num_fragments[j];
    }
    
//    cerr << "********************" << endl;
//    cerr << "initial MLE counts: " << mle_frag_counts << endl;

    ublas::matrix<double> mle_count_covar = iterated_exp_count_covariance;
    
    ublas::matrix<double> epsilon = ublas::zero_matrix<double>(transcripts.size(),transcripts.size());
	for (size_t i = 0; i < transcripts.size(); ++i)
	{
		epsilon(i,i) = 1e-6;
	}
	
	mle_count_covar += epsilon; // modify matrix to avoid problems during inverse
    
    double ret = cholesky_factorize(mle_count_covar);
    if (ret != 0)
    {
        fprintf(stderr, "Warning: Iterated expectation count covariance matrix cannot be cholesky factorized!\n");
        //fprintf(stderr, "Warning: FPKM covariance is not positive definite (ret = %lg)!\n", ret);
//        for (size_t j = 0; j < _abundances.size(); ++j)
//        {
//            _abundances[j]->status(NUMERIC_FAIL);
//        }
        return false;
    }
    
//    cerr << endl << "Cholesky factored covariance matrix: " << endl;
//    for (unsigned i = 0; i < _count_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (mle_count_covar, i);
//        cerr << i << " : " << _abundances[i]->num_fragments() << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }
    //cerr << "======" << endl;
    
    multinormal_generator<double> generator(mle_frag_counts, mle_count_covar);
    //vector<Eigen::VectorXd> multinormal_samples;
    
    
    assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * num_frag_assignments, Eigen::VectorXd::Zero(transcripts.size()));
    
    Eigen::VectorXd expected_generated_counts = Eigen::VectorXd::Zero(transcripts.size());
    
    
    boost::uniform_01<> uniform_dist;
    boost::mt19937 null_rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uniform_gen(null_rng, uniform_dist);

   
    for (size_t assign_idx = 0; assign_idx < num_frag_assignments; ++assign_idx)
    {
        
        boost::numeric::ublas::vector<double> random_count_assign = generator.next_rand();
        //cerr << random_count_assign << endl;

        for (size_t r_idx = 0; r_idx < random_count_assign.size(); ++r_idx)
        {
            if (random_count_assign(r_idx) < 0)
                random_count_assign(r_idx) = 0;
        }
        
        double total_sample_counts = accumulate(random_count_assign.begin(), random_count_assign.end(), 0);
        if (total_sample_counts > 0)
            random_count_assign = total_frag_counts * (random_count_assign / total_sample_counts);
        else
            random_count_assign = boost::numeric::ublas::zero_vector<double>(transcripts.size());
        

        
        //cerr << "*** sample around MLE: " << random_count_assign << endl;
        
        //double r = _abundances[j]->num_fragments();
        
        for (size_t gen_idx = 0; gen_idx < num_frag_count_draws; ++gen_idx)
        {
            Eigen::VectorXd generated_and_assigned_counts = Eigen::VectorXd::Zero(transcripts.size());
            
            for (size_t j = 0; j < transcripts.size(); ++j)
            {
                double r =  random_count_assign(j);
                
                //r = r < 1 ? 1 : r;
                
                if (r > 0)
                {
                    //double fit_var = _abundances[j]->mass_variance();
                    
                    double fit_var = total_var * (random_count_assign(j) / total_frag_counts);
                    double frags = random_count_assign(j);
                    
                    if (fit_var - frags > 1e-1)
                    {
                        r *= r;
                        double over_disp_scale = fit_var - frags;
                        r /= over_disp_scale;
                        
                        if (r == 0)
                        {
                            generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        //double p = _abundances[j]->num_fragments() / fit_var;
                        double p = r / (r + frags);
                        
                        negative_binomial_distribution<int, double> nb_j(r, p);

                        generated_and_assigned_counts(j) = nb_j(rng);

                    }
                    else
                    {
                        double after_decimal = r - (long)r;
                        //fprintf( stderr, "after decimal = %lg\n", after_decimal);
                        if (uniform_gen() < after_decimal)
                            r = floor(r);
                        else
                            r = ceil(r);

                        if (r == 0)
                        {
                            generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        boost::random::poisson_distribution<int, double> nb_j(r);

                        generated_and_assigned_counts(j) = nb_j(rng);
                        
                    }
                }
                else
                {
                    generated_and_assigned_counts(j) = 0;
                }
            }
            //cerr << "     assigned count sample: " << generated_and_assigned_counts.transpose() << endl;
            assigned_count_samples[assign_idx*num_frag_count_draws + gen_idx] = generated_and_assigned_counts;
        }
    }
        
    Eigen::VectorXd expected_counts = Eigen::VectorXd::Zero(transcripts.size());
    Eigen::VectorXd expected_relative_abundances = Eigen::VectorXd::Zero(transcripts.size());
    
    for (size_t i = 0; i < assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(assigned_count_samples[i](j)) && !isinf(assigned_count_samples[i](j)));
        }
        
        
        expected_counts += assigned_count_samples[i];
        //
        //expected_relative_abundances += relative_abundances[i];
    }
    if (assigned_count_samples.size() > 0)
    {
        expected_counts /= assigned_count_samples.size();
        //expected_generated_counts /= assigned_counts.size();
        //expected_relative_abundances /= assigned_counts.size();
    }
    
//    if (num_frag_assignments > 0)
//    {
//        expected_generated_counts /= num_frag_assignments;
//        //expected_relative_abundances /= assigned_counts.size();
//    }
    
    
    //    cerr << "======" << endl;
    //    cerr << "updated expected counts #1: " << endl;
    //    std::cerr << expected_counts << std::endl;
    //    cerr << "updated expected generated counts #1: " << endl;
    //    std::cerr << expected_generated_counts << std::endl;
    //    cerr << "======" << endl;
    
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            for (size_t k = 0 ; k < assigned_count_samples.size(); ++k)
            {
                double c = (assigned_count_samples[k](i) - expected_counts(i)) * (assigned_count_samples[k](j) - expected_counts(j));
                count_covariance(i,j) += c;
                
                assert (!isinf(count_covariance(i,j)) && !isnan(count_covariance(i,j)));
                //double r = (relative_abundances[k](i) - expected_relative_abundances(i)) * (relative_abundances[k](j) - expected_relative_abundances(j));
                //_kappa_covariance(i,j) +=
            }
        }
    }
    
    count_covariance /= assigned_count_samples.size();
        
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        // Make sure we aren't below the fit for the single isoform case
        if (count_covariance(i,i) < ceil(frag_variances[i]))
        {
            //fprintf(stderr, "Counts for %d (var = %lg) are underdispersed, reverting to fitted variance model (%lg)\n", i, _count_covariance(i,i), ceil(_abundances[i]->mass_variance()));
            count_covariance(i,i) = ceil(frag_variances[i]);
        }
        
        // Check that we aren't below what the Poisson model says we ought to be at
        if (count_covariance(i,i) < ceil(num_fragments[i] + iterated_exp_count_covariance(i,i)))
        {
            //fprintf(stderr, "Counts for %d (var = %lg) are underdispersed, reverting to additive variance model (%lg)\n", i, _count_covariance(i,i),  ceil(_abundances[i]->num_fragments() + _iterated_exp_count_covariance(i,i)));
            count_covariance(i,i) = ceil(num_fragments[i] + iterated_exp_count_covariance(i,i));
        }
        
//        long double count_var = 0.0;
        
//        // Check that we aren't below what the BNB model says we ought to be at
//        bool numerics_ok = estimate_count_variance(count_var,
//                                                   _abundances[i]->gamma(),
//                                                   _iterated_exp_count_covariance(i,i),
//                                                   num_fragments(),
//                                                   _abundances[i]->mass_variance(),
//                                                   _abundances[i]->effective_length());
//        //        if (numerics_ok == false)
//        //        {
//        //            fprintf(stderr, "Warning: BNB has no analytic solution\n");
//        //        }
//        
//        if (numerics_ok && _count_covariance(i,i) < ceil(count_var))
//        {
//            //fprintf(stderr, "Counts for %d (var = %lg) are underdispersed, reverting to additive variance model (%lg)\n", i, _count_covariance(i,i),  ceil(_abundances[i]->num_fragments() + _iterated_exp_count_covariance(i,i)));
//            _count_covariance(i,i) = ceil(count_var);
//        }
    }
    
    //    for (size_t i = 0; i < _abundances.size(); ++i)
    //    {
    //        _count_covariance(i,i) = ceil(_count_covariance(i,i));
    //    }
    
//    cerr << endl << "simulated count covariance: " << endl;
//    for (unsigned i = 0; i < _count_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (_count_covariance, i);
//        cerr << i << " : " << _abundances[i]->num_fragments() << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }
//    cerr << "======" << endl;
//    cerr << "updated expected counts: " << endl;
//    std::cerr << expected_counts << std::endl;
//    cerr << "======" << endl;

    //    for (size_t i = 0; i < num_count_draws; ++i)
    //    {
    //        cerr << generated_counts[i] << endl;
    //        
    //    }
    
    return true;
}

void AbundanceGroup::generate_fpkm_samples()
{
    double external_scale_factor = -1.0;
    
    double M = 0;
    
    for (set<shared_ptr<ReadGroupProperties const> >::iterator itr = _read_group_props.begin();
         itr != _read_group_props.end();
         ++itr)
    {
        shared_ptr<ReadGroupProperties const> rg_props = *itr;
        M += rg_props->normalized_map_mass();
        
        if (external_scale_factor < 0)
        {
            external_scale_factor = rg_props->external_scale_factor();
        }
        else
        {
            assert (external_scale_factor == rg_props->external_scale_factor());
        }
    }
    
    M /= _read_group_props.size();
    
    // set up individual vectors of FPKM samples for each abundance object in this group.
    vector<vector<double> > fpkm_sample_vectors(_abundances.size());
    vector<double> group_sum_fpkm_samples;
    
    vector<double> fpkm_means(_abundances.size(), 0);
    
    for(size_t i = 0; i < _assigned_count_samples.size(); ++i)
    {
        const Eigen::VectorXd sample = _assigned_count_samples[i];
        double total_fpkm = 0;
        
        for (size_t j = 0; j < sample.size(); ++j)
        {
            double fpkm_sample = sample[j] / M;
            
            if (_abundances[j]->effective_length() > 0)
            {
                fpkm_sample *= 1000000000;
                fpkm_sample /= _abundances[j]->effective_length();
                fpkm_sample /= external_scale_factor;
                //double standard_fpkm = _abundances[j]->FPKM();
                //fprintf(stderr, "count = %lg, fpkm = %lg, standard fpkm = %lg\n", sample[j], fpkm_sample, standard_fpkm);
            }
            else
            {
                fpkm_sample = 0;
            }
            
            assert (isnan(fpkm_sample) == false);
            
            fpkm_sample_vectors[j].push_back(fpkm_sample);
            fpkm_means[j] += fpkm_sample;
            total_fpkm += fpkm_sample;
        }
        
        if (effective_length() > 0)
        {
            group_sum_fpkm_samples.push_back(total_fpkm);
        }
        else
        {
            group_sum_fpkm_samples.push_back(0);
        }
    }
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        fpkm_means[i] /= _assigned_count_samples.size();
        _abundances[i]->fpkm_samples(fpkm_sample_vectors[i]);
        //fprintf(stderr, "standard fpkm = %lg, sample mean = %lg\n", _abundances[i]->FPKM(), fpkm_means[i]);
    }
    
    fpkm_samples(group_sum_fpkm_samples);
}

void AbundanceGroup::fit_gamma_distributions()
{
    // Now fit a gamma distribution to the FPKM samples for each abundance object in this group.
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        double i_k = 0;
        double i_theta = 0;
        
        const vector<double> fpkm_samples = _abundances[i]->fpkm_samples();
        
        bool good_fit = fit_gamma_dist(fpkm_samples, i_k, i_theta);
        if (fpkm_samples.empty() || good_fit == false)
        {
            _abundances[i]->status(NUMERIC_LOW_DATA);
        }
        else
        {
            _abundances[i]->fpkm_gamma_dist_k(i_k);
            _abundances[i]->fpkm_gamma_dist_theta(i_theta);
            //double fpkm_gamma_mean = i_k * i_theta;
            //double fpkm_gamma_var = i_k * i_theta * i_theta;
            //fprintf (stderr, "standard mean = %lg, standard var = %lg, sample mean = %lg, gamma_mean = %lg, gamma_var = %lg\n", _abundances[i]->FPKM(),_abundances[i]->FPKM_variance(), fpkm_means[i], fpkm_gamma_mean, fpkm_gamma_var);
        }
    }
    
    double group_k = 0;
    double group_theta = 0;
    bool good_fit = fit_gamma_dist(fpkm_samples(), group_k, group_theta);
    if (good_fit == false)
    {
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            _abundances[i]->status(NUMERIC_FAIL);
        }
    }
    else
    {
        fpkm_gamma_dist_k(group_k);
        fpkm_gamma_dist_theta(group_theta);
    }
}

void AbundanceGroup::calculate_FPKM_covariance()
{
	if (mass_fraction() == 0 || effective_length() == 0)
	{
		_fpkm_covariance = ublas::zero_matrix<double>(_abundances.size(), _abundances.size());
		return;
	}
    
    //long double M = num_fragments()/mass_fraction();
    
    //estimate_count_covariance();
    
    long double total_var = 0.0;
    
    double abundance_weighted_length = 0.0;
    double total_abundance = 0.0;
    
    double external_scale_factor = -1.0;
    
    double M = 0;
    
    for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin();
         itr != _count_per_replicate.end();
         ++itr)
    {
        shared_ptr<ReadGroupProperties const> rg_props = itr->first;
        M += rg_props->normalized_map_mass();
        
        if (external_scale_factor < 0)
        {
            external_scale_factor = rg_props->external_scale_factor();
        }
        else
        {
            assert (external_scale_factor == rg_props->external_scale_factor());
        }
    }
    
    M /= _count_per_replicate.size();
    
    for (size_t j = 0; j < _abundances.size(); ++j)
    {
        abundance_weighted_length += _abundances[j]->effective_length() * _abundances[j]->FPKM();
        total_abundance += _abundances[j]->FPKM();
        
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            _fpkm_covariance(i,j) = _count_covariance(i,j);
            
            // FPKMs need to be on the external scale, so we can compare them
            // between conditions.  Counts are internally scaled up until here.
            _fpkm_covariance(i,j) *= 1.0 / (external_scale_factor * external_scale_factor);
            assert (!isinf(_count_covariance(i,j)) && !isnan(_fpkm_covariance(i,j)));
            
            long double length_i = _abundances[i]->effective_length();
            long double length_j = _abundances[j]->effective_length();
            assert (!isinf(length_i) && !isnan(length_i));
            assert (!isinf(length_j) && !isnan(length_j));
            if (length_i > 0 && length_j > 0)
            {
                _fpkm_covariance(i,j) *=
                    ((1000000000.0 / (length_j *M)))*((1000000000.0 / (length_i *M)));
                assert (!isinf(_fpkm_covariance(i,j)) && !isnan(_fpkm_covariance(i,j)));
            }
            else
            {
                _fpkm_covariance(i,j) = 0.0;
            }
            
            if (i == j)
            {
                double fpkm = _abundances[i]->FPKM();
                double fpkm_var = _fpkm_covariance(i,j);
                if (fpkm > 0 && fpkm_var == 0 )
                {
                    cerr << _count_covariance << endl;
                }
                assert (fpkm == 0 || fpkm_var > 0 || _abundances[i]->status() != NUMERIC_OK);
                assert (!isinf(fpkm_var) && !isnan(fpkm_var));
                _abundances[i]->FPKM_variance(fpkm_var);
                _abundances[i]->num_fragment_var(_count_covariance(i,j));
                
            }
            
            total_var += _fpkm_covariance(i,j);

            assert (!isinf(_fpkm_covariance(i,j)) && !isnan(_fpkm_covariance(i,j)));
        }
    }
    
    _FPKM_variance = total_var;
    
    if (final_est_run && library_type != "transfrags")
    {
        ublas::matrix<double> test = _fpkm_covariance;
        double ret = cholesky_factorize(test);
        if (ret != 0 || (_FPKM_variance < 0 && status() == NUMERIC_OK))
        {
            //fprintf(stderr, "Warning: FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            for (size_t j = 0; j < _abundances.size(); ++j)
            {
                _abundances[j]->status(NUMERIC_FAIL);
            }
        }

        assert (FPKM() == 0 || _FPKM_variance > 0 || status() != NUMERIC_OK);
    }
    
    assert (!isinf(_FPKM_variance) && !isnan(_FPKM_variance));
}

void AbundanceGroup::calculate_conf_intervals()
{    
    // We only really ever call this function for primary abundance groups
    // (i.e. the transcript groups and read bundles with which we calculate
    // transcript MLE expression levels.  Genes, TSS groups, etc get broken
    // off of primary bundles, so we should not call this function on those
    // secondary groups.  The group splitting code needs to manage the task
    // of splitting up all the variout covariance matrices we're calculating
    // here.
	if (status() == NUMERIC_OK)
	{
		// This will compute the transcript level FPKM confidence intervals
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
            long double fpkm_var = _abundances[j]->FPKM_variance();
            double FPKM_hi = 0.0;      
            double FPKM_lo = 0.0;
            if (_abundances[j]->status() != NUMERIC_FAIL)
            {
                FPKM_hi = _abundances[j]->FPKM() + 2 * sqrt(fpkm_var);
                FPKM_lo = max(0.0, (double)(_abundances[j]->FPKM() - 2 * sqrt(fpkm_var)));
                if (!(FPKM_lo <= _abundances[j]->FPKM() && _abundances[j]->FPKM() <= FPKM_hi))
                {
                    //fprintf(stderr, "Error: confidence intervals are illegal! var = %Lg, fpkm = %lg, lo = %lg, hi %lg, status = %d\n", fpkm_var, _abundances[j]->FPKM(), FPKM_lo, FPKM_hi, _abundances[j]->status());
                }
                assert (FPKM_lo <= _abundances[j]->FPKM() && _abundances[j]->FPKM() <= FPKM_hi);
                ConfidenceInterval conf(FPKM_lo, FPKM_hi);
                _abundances[j]->FPKM_conf(conf);
                //_abundances[j]->FPKM_variance(fpkm_var);
            }
            else
            {
                // we shouldn't be able to get here
                assert(false);
                // TODO: nothing to do here?
            }
		}
		
        // Now build a confidence interval for the whole abundance group
		double group_fpkm = FPKM();
		if (group_fpkm > 0.0)
		{
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
        double max_fpkm = 0.0;
        //double min_fpkm = 1e100;
		BOOST_FOREACH(shared_ptr<Abundance> pA, _abundances)
		{
			double FPKM_hi;
			double FPKM_lo;
			if (pA->effective_length() > 0)
			{
                double norm_frag_density = 1000000000;
                norm_frag_density /= pA->effective_length();
                
                norm_frag_density *= mass_fraction();
                double fpkm_high = norm_frag_density;
                
                double var_fpkm = fpkm_high; 
                
				FPKM_hi = fpkm_high + 2 * sqrt(var_fpkm);
				FPKM_lo = 0.0;
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				assert (FPKM_lo <= pA->FPKM() && pA->FPKM() <= FPKM_hi);
				pA->FPKM_conf(conf);
                //pA->FPKM_variance(var_fpkm);
				max_fpkm = max(sum_transfrag_FPKM_hi, FPKM_hi);
			}
			else
			{
				FPKM_hi = 0.0;
				FPKM_lo = 0.0;
				ConfidenceInterval conf(0.0, 0.0);
				pA->FPKM_conf(conf);
                //pA->FPKM_variance(0.0);
			}
            
		}
		// In the case of a numeric failure, the groups error bars need to be 
		// set such that 
		FPKM_conf(ConfidenceInterval(0.0, max_fpkm + 2 * sqrt(FPKM_variance())));
	}
}

void compute_cond_probs_and_effective_lengths(const vector<MateHit>& alignments,
															  vector<shared_ptr<Abundance> >& transcripts,
															  vector<shared_ptr<Abundance> >& mapped_transcripts)
{		
	int N = transcripts.size();
	int M = alignments.size();

	vector<vector<char> > compatibilities(N, vector<char>(M,0));
	compute_compatibilities(transcripts, alignments, compatibilities);

	for (int j = 0; j < N; ++j) 
    {
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


double trace(const ublas::matrix<double>& m)
{
    
    double t = 0.0;
    for (size_t i = 0.0; i < m.size1(); ++i)
    {
        t += m(i,i);
    }
    
    return t;
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
		BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
		{
			ab->gamma(0);
		}
		_gamma_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                       transcripts.size());
        _count_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                       transcripts.size());
        _iterated_exp_count_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                       transcripts.size());
        _fpkm_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                      transcripts.size());
        
		return true;
    }
	
	vector<double> gammas;
	
	verbose_msg( "Calculating intial MLE\n");
	
	AbundanceStatus mle_success = gamma_mle(mapped_transcripts,
                                                  nr_alignments,
                                                  log_conv_factors,
                                                  gammas);
	
	verbose_msg( "Tossing likely garbage isoforms\n");
	
	for (size_t i = 0; i < gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			verbose_msg("Warning: isoform abundance is NaN!\n");
		}
	}
	
    double locus_mass = 0.0;
   
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        const MateHit& alignment = nr_alignments[i];
        locus_mass += alignment.collapse_mass();
    }
    
	vector<shared_ptr<Abundance> > filtered_transcripts = mapped_transcripts;
	vector<double> filtered_gammas = gammas;
	filter_junk_isoforms(filtered_transcripts, filtered_gammas, mapped_transcripts, locus_mass);
	
	if (filtered_transcripts.empty())
	{
		//gammas = vector<double>(transfrags.size(), 0.0);
		BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
		{
			ab->gamma(0);
		}
		_gamma_covariance = ublas::zero_matrix<double>(transcripts.size(), 
													  transcripts.size());
        _count_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                       transcripts.size());
        _iterated_exp_count_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                                    transcripts.size());
        _fpkm_covariance = ublas::zero_matrix<double>(transcripts.size(), 
                                                      transcripts.size());
        
		return true;
	}
	
    if (filtered_transcripts.size() != mapped_transcripts.size())
    {    
        filtered_gammas.clear();
        
        verbose_msg( "Revising MLE\n");
        
        mle_success = gamma_mle(filtered_transcripts,
                                        nr_alignments,
                                        log_conv_factors, 
                                        filtered_gammas);
    }

	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(filtered_gammas[i]))
		{
			verbose_msg("Warning: isoform abundance is NaN!\n");
		}
	}
	
	size_t N = transcripts.size();
	
    set<shared_ptr<ReadGroupProperties const> > rg_props;
	for (size_t i = 0; i < nr_alignments.size(); ++i)
	{
        rg_props.insert(nr_alignments[i].read_group_props());
	}
    
    AbundanceStatus map_success = NUMERIC_OK;
	if (final_est_run) // Only on last estimation run.
	{
        ublas::vector<double> gamma_mle(filtered_gammas.size());
        std::copy(filtered_gammas.begin(), filtered_gammas.end(), gamma_mle.begin());
    }
	
	for (size_t i = 0; i < filtered_gammas.size(); ++i)
	{
		if (isnan(gammas[i]))
		{
			verbose_msg( "Warning: isoform abundance is NaN!\n");
			map_success = NUMERIC_FAIL;
		}
	}
	
	// Now we need to fill in zeros for the isoforms we filtered out of the 
	// MLE/MAP calculation
	vector<double> updated_gammas = vector<double>(N, 0.0);
    
    
	ublas::matrix<double> updated_gamma_cov;
	updated_gamma_cov = ublas::zero_matrix<double>(N, N);
    
    ublas::matrix<double> updated_count_cov;
    updated_count_cov = ublas::zero_matrix<double>(N, N);
    ublas::matrix<double> updated_iterated_exp_count_cov;
    updated_iterated_exp_count_cov = ublas::zero_matrix<double>(N, N);
    ublas::matrix<double> updated_fpkm_cov;
    updated_fpkm_cov = ublas::zero_matrix<double>(N, N);
	
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
            //cerr << updated_gammas[i] << ",";
            
			for (size_t j = 0; j < N; ++j)
			{
				if (scaff_present[j] != N)
				{
					updated_gamma_cov(i,j) = _gamma_covariance(scaff_present[i],
															   scaff_present[j]);
                    
                    updated_iterated_exp_count_cov(i,j) = _iterated_exp_count_covariance(scaff_present[i],
                                                                                         scaff_present[j]);
                    // Should still be empty but let's do these for consistency:
                    updated_count_cov(i,j) = _count_covariance(scaff_present[i],
                                                               scaff_present[j]);
                    updated_fpkm_cov(i,j) = _fpkm_covariance(scaff_present[i],
                                                             scaff_present[j]);
					assert (!isinf(updated_gamma_cov(i,j)));
					assert (!isnan(updated_gamma_cov(i,j)));
				}
			}
		}
	}
    
    //cerr << endl;
	
    AbundanceStatus numeric_status = NUMERIC_OK;
    if (mle_success == NUMERIC_LOW_DATA)
    {
        numeric_status = NUMERIC_LOW_DATA;
    }
    else if (mle_success == NUMERIC_FAIL)
    {
        numeric_status = NUMERIC_FAIL;
    }
    else
    {
        assert (mle_success == NUMERIC_OK);
        if (map_success == NUMERIC_FAIL)
        {
            numeric_status = NUMERIC_FAIL;
        }
        else if (map_success == NUMERIC_LOW_DATA)
        {
            numeric_status = NUMERIC_LOW_DATA;
        }
        // otherwise, we're cool.
    }
    
    
    
	// All scaffolds that go in get abundances, but those that get "filtered"
	// from the calculation get zeros.
	//gammas = updated_gammas;
	for (size_t i = 0; i < _abundances.size(); ++i)
	{
		_abundances[i]->gamma(updated_gammas[i]);
		_abundances[i]->status(numeric_status);
	}
	_gamma_covariance = updated_gamma_cov;
    _count_covariance = updated_count_cov;
    _iterated_exp_count_covariance = updated_iterated_exp_count_cov;
    _fpkm_covariance = updated_fpkm_cov;
	
	return (status() == NUMERIC_OK);
}

void calculate_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                const Eigen::MatrixXd& transcript_cond_probs,
                                const Eigen::VectorXd& proposed_gammas,
                                Eigen::MatrixXd& assignment_probs)
{
//    vector<double> u(nr_alignments.size());
//    for (size_t i = 0; i < alignment_multiplicity.size(); ++i)
//    {
//        u[i] = nr_alignments[i].collapse_mass();
//    }
    
    //ublas::vector<double> total_cond_prob = ublas::prod(proposed_gammas,transcript_cond_probs);
    Eigen::VectorXd total_cond_prob = proposed_gammas.transpose() * transcript_cond_probs ;
    
    // Compute the marginal conditional probability for each fragment against each isoform
    //ublas::matrix<double>  marg_cond_prob = ublas::zero_matrix<double>(transcript_cond_probs.size1(), transcript_cond_probs.size2());
    Eigen::MatrixXd marg_cond_prob(transcript_cond_probs.rows(), transcript_cond_probs.cols());
    
    for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
        marg_cond_prob.array().col(i) = proposed_gammas.array() * transcript_cond_probs.array().col(i);
        
        if (total_cond_prob(i) > 0)
        {
            marg_cond_prob.array().col(i) /= total_cond_prob(i);
            //column(marg_cond_prob,i) /= total_cond_prob(i);
        }
    }
    
    assignment_probs = marg_cond_prob;
}

void calculate_average_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                        const Eigen::MatrixXd& transcript_cond_probs,
                                        const Eigen::VectorXd& proposed_gammas,
                                        Eigen::MatrixXd& assignment_probs)
{
    //ublas::vector<double> total_cond_prob = ublas::prod(proposed_gammas,transcript_cond_probs);
    Eigen::VectorXd total_cond_prob = proposed_gammas.transpose() * transcript_cond_probs ;
    
    // Compute the marginal conditional probability for each fragment against each isoform
    //ublas::matrix<double>  marg_cond_prob = ublas::zero_matrix<double>(transcript_cond_probs.size1(), transcript_cond_probs.size2());
    Eigen::MatrixXd marg_cond_prob(transcript_cond_probs.rows(), transcript_cond_probs.cols());
    
    for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
        marg_cond_prob.array().col(i) = proposed_gammas.array() * transcript_cond_probs.array().col(i);
        
        if (total_cond_prob(i) > 0)
        {
            marg_cond_prob.array().col(i) /= total_cond_prob(i);
            //column(marg_cond_prob,i) /= total_cond_prob(i);
        }
    }
    
    assignment_probs = Eigen::MatrixXd::Zero(proposed_gammas.size(), proposed_gammas.size());
    
    vector<double> num_compatible(proposed_gammas.size(), 0);
    
    for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
        for (size_t j = 0; j < proposed_gammas.size(); ++j)
        {
            if (marg_cond_prob(j,i) > 0)
            {
                assignment_probs.col(j) += marg_cond_prob.col(i) * alignment_multiplicities[i];
                num_compatible[j] += alignment_multiplicities[i];
            }
        }
    }
    
//    cerr << "marg cond prob:" << endl;
//    cerr << marg_cond_prob << endl;
    
    for (size_t j = 0; j < proposed_gammas.size(); ++j)
    {
        if (num_compatible[j] > 0)
            assignment_probs.col(j) /= num_compatible[j];
    }
    
//    cerr << "multiplicities:" << endl;
//    cerr << alignment_multiplicities << endl;
//    cerr << "avg matrix:" << endl;
//    cerr << assignment_probs << endl;
    //assignment_probs = marg_cond_prob;
}



void calculate_iterated_exp_count_covariance(const vector<double>& gammas,
                                             const vector<MateHit>& nr_alignments,
                                             const vector<shared_ptr<Abundance> >& transcripts,
                                             ublas::matrix<double>& count_covariance)
{
    // Now calculate the _iterated_exp_count_covariance matrix via iterated expectation
    vector<vector<double> > cond_probs(transcripts.size(), vector<double>());
    for(size_t j = 0; j < transcripts.size(); ++j)
    {
        cond_probs[j]= *(transcripts[j]->cond_probs());
    }
    
    vector<double> u(nr_alignments.size());
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        u[i] = nr_alignments[i].collapse_mass();
    }
    
    count_covariance = ublas::zero_matrix<double>(transcripts.size(), transcripts.size());
    
    ublas::vector<double> total_cond_prob = ublas::zero_vector<double>(nr_alignments.size());
    
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        // the replicate gamma mles might not be available, if one of the
        // replicates returned an error, we'll consider all to be unreliable
        for (size_t j = 0; j < cond_probs.size(); ++j)
        {
            if (cond_probs[j][i] > 0)
            {
                total_cond_prob(i) += gammas[j] * cond_probs[j][i];
                assert (!isnan(total_cond_prob(i) && ! isinf(total_cond_prob(i))));
            }
        }
    }
    
    // Compute the marginal conditional probability for each fragment against each isoform
    ublas::matrix<double>  marg_cond_prob = ublas::zero_matrix<double>(transcripts.size(), nr_alignments.size());
    
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        // the replicate gamma mles might not be available, if one of the
        // replicates returned an error, we'll consider all to be unreliable
        for (size_t j = 0; j < cond_probs.size(); ++j)
        {
            if (total_cond_prob(i))
            {
                if (cond_probs[j][i] > 0)
                {
                    marg_cond_prob(j,i) = (gammas[j] * cond_probs[j][i])/total_cond_prob(i);
                }
            }
        }
    }
    
    double total_var = 0.0;
    
    ublas::vector<double> expected_counts = ublas::zero_vector<double>(cond_probs.size());

    //iterate over fragments
    for (size_t i = 0; i < marg_cond_prob.size2(); ++i)
    {
        
        // iterate over transcripts
        for (size_t j = 0; j < marg_cond_prob.size1(); ++j)
        {
            double c_j_i = marg_cond_prob(j,i);
            
            expected_counts(j) += u[i] * marg_cond_prob(j,i);
            
            //if (c_j_i == 0 || c_j_i == 1.0)
            //    continue;
            
            for (size_t k = 0; k < marg_cond_prob.size1(); ++k)
            {
                double c_k_i = marg_cond_prob(k,i);
                //if (c_k_i == 0 || c_k_i == 1.0)
                //    continue;
                
                if (j == k)
                {
                    if (c_k_i != 0 && c_k_i != 1.0)
                    {
                        double var = u[i] * c_k_i * (1.0 - c_k_i);
                        count_covariance(k,k) += var;
                        assert (var >= 0);
                        assert (!isnan(var) && !isinf(var));
                        total_var += var;
                    }
                }
                else
                {
                    if (c_k_i != 0 && c_k_i != 1.0 &&
                        c_j_i != 0 && c_j_i != 1.0)
                    {
                        double covar = -u[i] * c_k_i * c_j_i;
                        assert (covar <= 0);
                        assert (!isnan(covar) && !isinf(covar));
                        count_covariance(k,j) += covar;
                    }
                } 
            }
        }

    }
    
//    // take care of little rounding errors
//    for (size_t i = 0; i < _iterated_exp_count_covariance.size1(); ++i)
//    {
//        for (size_t j = 0; j < _iterated_exp_count_covariance.size2(); ++j)
//        {
//            if (i == j)
//            {
//                double c = _iterated_exp_count_covariance(i,j);
//                if (c < 0)
//                    _iterated_exp_count_covariance(i,j) = 0;
//                //assert(c >= 0);
//            }
//            else
//            {
//                double c = _iterated_exp_count_covariance(i,j);
//                if (c > 0)
//                    _iterated_exp_count_covariance(i,j) = 0;
//                //assert(c <= 0);
//            }
//        }
//        _abundances[i]->num_fragment_uncertainty_var(_iterated_exp_count_covariance(i,i));
//    }
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
	
	double S_FPKM = 0.0;
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		if (pA->effective_length() > 0)
		{
			S_FPKM += pA->FPKM();
		}
	}
	
    //fprintf (stderr, "*********\n");
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		if (S_FPKM > 0)
		{
			pA->kappa(pA->FPKM() / S_FPKM);
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
            double L = _abundances[k]->effective_length() * 
            _abundances[m]->effective_length();
            if (L == 0.0)
            {
                _kappa_covariance(k,m) = 0.0;
            }
            else if (m == k)
            {
                // Use the modeled count variance here instead
                double kappa_var;
                if (S_FPKM)
                {
                    kappa_var = _fpkm_covariance(k,k) / (S_FPKM * S_FPKM);
                }
                else
                {
                    kappa_var = 0.0;
                }
                
                if (isnan(kappa_var) || isinf(kappa_var)) // to protect against underflow
                    kappa_var = 0;
                _kappa_covariance(k,m) = kappa_var;
            }
            else
            {
                double kappa_covar;
                if (S_FPKM)
                {
                    kappa_covar = _fpkm_covariance(k,m) / (S_FPKM * S_FPKM);
                }
                else
                {
                    kappa_covar = 0.0;
                }
                _kappa_covariance(k,m) = kappa_covar;
            }
		}
	}

    
//	size_t num_members = _abundances.size();
//	_kappa_covariance = ublas::zero_matrix<double>(num_members, 
//											  num_members);
//    if (FPKM() == 0)
//    {
//        for (size_t k = 0; k < num_members; ++k)
//        {
//            _abundances[k]->kappa(0);
//        }
//        return;
//    }
//
    
    // FIXME:
    size_t num_count_draws = _assigned_count_samples.size();
    vector<Eigen::VectorXd > relative_abundances (num_count_draws, Eigen::VectorXd::Zero(num_members));
    
    // We'll use the effective lengths to transform counts into relative abundances,
    // and then use that to calculate the kappa variances and covariances.
    Eigen::VectorXd effective_length_recip = Eigen::VectorXd::Zero(_abundances.size());
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        if (_abundances[i]->effective_length() > 0)
            effective_length_recip(i) = 1.0 / _abundances[i]->effective_length();
    }
    
    for (size_t i = 0; i < num_count_draws; ++i)
    {
        
        Eigen::VectorXd relative_abundance = effective_length_recip.array() * _assigned_count_samples[i].array();
        double total = relative_abundance.sum();
        if (total > 0)
            relative_abundance /= total;
        //cerr << relative_abundance.transpose() << endl;
        relative_abundances[i] = relative_abundance;
    }
    
//    
//    Eigen::VectorXd expected_relative_abundances = Eigen::VectorXd::Zero(_abundances.size());
//    
//    for (size_t i = 0; i < relative_abundances.size(); ++i)
//    {
//        expected_relative_abundances += relative_abundances[i];
//    }
//    
//    if (relative_abundances.size() > 0)
//    {
//        expected_relative_abundances /= relative_abundances.size();
//    }
//    
//    for (size_t k = 0; k < num_members; ++k)
//    {
//        _abundances[k]->kappa(expected_relative_abundances(k));
//    }
    
//    cerr << "======" << endl;
//    cerr << "updated expected relative abundances: " << endl;
//    std::cerr << expected_relative_abundances << std::endl;
//    cerr << "======" << endl;
//    
//    cerr << "simulated kappa deviations: " << endl;
//    for (unsigned i = 0; i < _count_covariance.size1 (); ++ i) 
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (_kappa_covariance, i);
//        cerr << i << " : " << _abundances[i]->kappa() << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }
    
//    for (size_t i = 0; i < _abundances.size(); ++i)
//    {
//        for (size_t j = 0; j < _abundances.size(); ++j)
//        {
//            for (size_t k = 0 ; k < relative_abundances.size(); ++k)
//            {
//                double r = (relative_abundances[k](i) - expected_relative_abundances(i)) * (relative_abundances[k](j) - expected_relative_abundances(j));
//                assert (r <= 1.0 && r >= -1);
//                _kappa_covariance(i,j) += r;
//                //assert (_kappa_covariance(i,j) >= -1 * relative_abundances.size() && _kappa_covariance(i,j) <= relative_abundances.size());
//            }
//        }
//    }
    
//    cerr << "simulated kappa deviations: " << endl;
//    for (unsigned i = 0; i < _count_covariance.size1 (); ++ i) 
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (_kappa_covariance, i);
//        cerr << i << " : " << _abundances[i]->kappa() << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }

    
//    _kappa_covariance /= relative_abundances.size();
    
    vector<double> js_samples;
    
    ublas::vector<double> kappa_mean(_abundances.size());
    for (size_t j = 0; j < _abundances.size(); ++j)
    {
        kappa_mean(j) = _abundances[j]->kappa();
    }
    
    ublas::matrix<double> kappa_cov_chol = _kappa_covariance;
    double ret = cholesky_factorize(kappa_cov_chol);
    if (ret == 0)
    {
        multinormal_generator<double> generator(kappa_mean, kappa_cov_chol);
        vector<Eigen::VectorXd> multinormal_samples;
        
        generate_importance_samples(generator, multinormal_samples, 10000, false);

        // We used to sample the JS using the real assigned count samples, but
        // that's not quite as accurate as simulating from a multinomial built from
        // the bounded covariance matrices.
        
        //generate_null_js_samples(relative_abundances, 100000, js_samples);
        generate_null_js_samples(multinormal_samples, 100000, js_samples);
        
        _null_js_samples = js_samples;
        //if (_null_js_samples.size() > 0)
        //    fprintf(stderr, "Max JS from null: %lg\n",_null_js_samples.back()); 
    }
    else
    {
        _null_js_samples.clear();
    }
}

void get_alignments_from_scaffolds(const vector<shared_ptr<Abundance> >& abundances,
								   vector<MateHit>& alignments)
{
	set<const MateHit*> hits_in_gene_set;

	BOOST_FOREACH(shared_ptr<Abundance> pA, abundances)
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

//void round(Eigen::VectorXd&  p) {
//	
//	double KILLP = 0; // kill all probabilities below this
//	
//	for (size_t i = 0; i < p.size(); ++i) {
//		if (p(i) < KILLP) 
//			p(i) = 0;
//	}
//}

void Estep (int N, 
			int M, 
			const Eigen::VectorXd& p,
			Eigen::MatrixXd& U,
			const Eigen::MatrixXd& cond_probs,
			const Eigen::VectorXd& alignment_multiplicities) {
	// given p, fills U with expected frequencies
	int i,j;	
    
//    Eigen::VectorXd frag_prob_sums = Eigen::VectorXd::Zero(M);
//    
//    for (j = 0; j < N; ++j) 
//    {
//        for (i = 0; i < M; ++i) 
//        {
//            frag_prob_sums(i) += cond_probs(j,i) * p(j);
//        }
//    }
//    
//    for (i = 0; i < M; ++i) 
//    {
//        frag_prob_sums(i) = frag_prob_sums(i) ? (1.0 / frag_prob_sums(i)) : 0.0;
//    }
//    
//    for (j = 0; j < N; ++j) 
//    {
//        for (i = 0; i < M; ++i) 
//        {
//            double ProbY = frag_prob_sums(i);
//            double exp_i_j = alignment_multiplicities(i) * cond_probs(j,i) * p(j) * ProbY;
//            U(j,i) = exp_i_j;
//        }
//    }
    
    Eigen::VectorXd frag_prob_sums;// = Eigen::VectorXd::Zero(M);
    
//    for (j = 0; j < N; ++j) 
//    {
//        for (i = 0; i < M; ++i) 
//        {
//            frag_prob_sums(i) += cond_probs(j,i) * p(j);
//            assert (!isnan(cond_probs(j,i)) && !isinf(cond_probs(j,i)));
//        }
//    }
    
    frag_prob_sums = p.transpose() * cond_probs;
    
    for (i = 0; i < M; ++i) 
    {
        assert (!isnan(frag_prob_sums(i)) && !isinf(frag_prob_sums(i)));
        //double x = frag_prob_sums(i);
        frag_prob_sums(i) = frag_prob_sums(i) ? (1.0 / frag_prob_sums(i)) : 0.0;
        if (isnan(frag_prob_sums(i)) || isinf(frag_prob_sums(i)))
        {
            frag_prob_sums(i) = 0; // protect against overflow/underflow
        }
    }
    
    Eigen::VectorXd x = frag_prob_sums.array() * alignment_multiplicities.array();
//    Eigen::MatrixXd y = p.transpose() * cond_probs;
//    Eigen::MatrixXd UU = y.array() * x.array();
    
//    cerr << UU <<endl;
//    for (j = 0; j < N; ++j) 
//    {
//        for (i = 0; i < M; ++i) 
//        {
//            //double ProbY = frag_prob_sums(i);
//            double exp_i_j = cond_probs(j,i) * p(j) * x(i);
//            U(j,i) = exp_i_j;
//        }
//    }
    
    U = Eigen::MatrixXd(N,M);

    for (i = 0; i < M; ++i) 
    {
        //double ProbY = frag_prob_sums(i);
        //double exp_i_j = cond_probs(j,i) * p(j) * x(i);
        //U.r = exp_i_j;
        U.col(i) = cond_probs.col(i).array() * p.array();
    }
    for (j = 0; j < N; ++j) 
    {
        U.array().row(j) *= x.transpose().array();
    }
    
    //cerr << UU << endl;
    //cerr << "==========" << endl;
    //cerr << U << endl;
    //cerr << "**********" << endl;
    //Eigen::ArrayXXd a_cond_probs = cond_probs.array();
    //Eigen::ArrayXXd a_p = p.array();
    //Eigen::ArrayXXd a_U = a_p * a_cond_probs.colwise();
    //U = (.colwise() * p.transpose().array()).matrix();
    //U = cond_probs.colwise() * p.array();
    //U = alignment_multiplicities.array() * p.array() * cond_probs.array() * frag_prob_sums.array();
    //U = (frag_prob_sums.array() * alignment_multiplicities.array()).rowwise() * p.array().colwise() * cond_probs.array().colwise();
}


void Mstep (int N, 
            int M, 
            Eigen::VectorXd& p, 
            const Eigen::MatrixXd& U) 
{
	Eigen::VectorXd v;// = Eigen::VectorXd::Zero(N);
    
	double m = 0;
    
	v = U.rowwise().sum();
    m = v.colwise().sum()(0);
    
	if (m)
	{
		p = v / m;
	}
	else
	{
        p = Eigen::VectorXd::Zero(N);
	}
}
 

double logLike (int N, 
				int M, 
				Eigen::VectorXd& p,
				const Eigen::MatrixXd& cond_prob, 
				const Eigen::VectorXd& alignment_multiplicities,
                const vector<double>& log_conv_factors) {
	//int i,j;
	
	double ell = accumulate(log_conv_factors.begin(), log_conv_factors.end(), 0.0);
	//double Prob_Y;
//	for (int i= 0; i < M; i++) {
//		Prob_Y = 0;
//		for (int j= 0; j < N; j++) {
//			Prob_Y += cond_prob(j,i) * p(j);
//		}
//		if (Prob_Y > 0) {
//			ell += (alignment_multiplicities(i) * log(Prob_Y));
//		}
//	}
    Eigen::VectorXd Prob_Ys = p.transpose() * cond_prob;
    Eigen::VectorXd log_Prob_Ys = Prob_Ys.array().log();
    //log_Prob_Ys *= alignment_multiplicities;
    for (int i = 0; i < M; i++) {
        if (Prob_Ys(i) > 0)
        {
            ell += alignment_multiplicities(i) * log_Prob_Ys(i);
        }
    }
	return ell;
}

double EM(int N, int M, 
          Eigen::VectorXd&  newP, 
          const Eigen::MatrixXd& cond_prob, 
		  const Eigen::VectorXd& alignment_multiplicities,
          vector<double> const & log_conv_factors,
          bool& converged) 
{
    converged = true;
	//double sum = 0;
	double newEll = 0;
    Eigen::VectorXd p = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd  U = Eigen::MatrixXd::Zero(N,M);
	double ell = 0; 
	int iter = 0;
	int j;
    
    if (N == 0 || M == 0)
        return NUMERIC_OK;
    
    for (j = 0; j < N; ++j) {
        //p[j] = drand48();
        //sum += p[j];
        p(j) = 1.0/(double)N;
    }
  
    

    static const double ACCURACY = mle_accuracy; // convergence for EM
    
	while (((iter <= 2) || (abs(ell - newEll) > ACCURACY)) && (iter < max_mle_iterations)) {
		if (iter > 0) {
			//round(newP);
			p = newP;
			ell = newEll;
		}
		
		Estep(N, M, p, U, cond_prob, alignment_multiplicities); //  fills U
		Mstep(N, M, newP,U); // fills p
		
		newEll = logLike(N, M, newP, cond_prob,alignment_multiplicities, log_conv_factors);
		
		//fprintf(stderr, "%d\t%lf\n", iter, newEll);
		
		//printf("%.3f %.3f %.3f ", newP[0], newP[1], newP[2]);
		//printf("%.3f %.3f %.3f ", newP[3], newP[4], newP[5]);
		//printf("%.3f %.3f %.3f\n", newP[6], newP[7], newP[8]);
		iter++;
	}
	if (iter >= max_mle_iterations)
    {
		verbose_msg("Warning: ITERMAX reached in abundance estimation, estimation hasn't fully converged\n");
        converged = false;
    }
    verbose_msg("Convergence reached in %d iterations \n", iter);
	return newEll;
}

void compute_fisher(const vector<shared_ptr<Abundance> >& transcripts,
					const ublas::vector<double>& abundances,
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
			alpha *= abundances(j);
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

AbundanceStatus compute_posterior_expectation(const vector<ublas::vector<double> >& weighted_samples,
                                              const vector<pair<size_t, double> >& sample_weights,
                                              ublas::vector<double>& expectation,
                                              long double& log_total_weight)
{
    ublas::vector<double> log_expectation(expectation.size());
    log_total_weight = 0.0;
    
    // compute the weighted sum of the samples from the proposal distribution,
    // but keep the result in log space to avoid underflow.
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
    
    if (log_total_weight == 0 || sample_weights.size() < 100)
	{
		verbose_msg("Warning: restimation failed, importance samples have zero weight.\n\tResorting to MLE and observed Fisher\n");
		return NUMERIC_FAIL;
	}
	
    // compute the weighted mean, and transform back out of log space
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
			return NUMERIC_FAIL;
		}
	}
	
	for (size_t j = 0; j < expectation.size(); ++j) 
	{
		if (expectation(j) < 0)
			expectation(j) = 0;
	}
	
	long double m = sum(expectation);
	
	if (m == 0 || isinf(m) || isnan(m))
	{
		verbose_msg("Warning: restimation failed, could not renormalize MAP estimate\n\tResorting to MLE and observed Fisher.");
		return NUMERIC_FAIL;
	}
	
	for (size_t j = 0; j < expectation.size(); ++j) {
		expectation(j) = expectation(j) / m;
	}
    return NUMERIC_OK;
}

AbundanceStatus AbundanceGroup::calculate_per_replicate_abundances(vector<shared_ptr<Abundance> >& transcripts,
                                                                   const vector<MateHit>& nr_alignments,
                                                                   const vector<double>& log_conv_factors,
                                                                   std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >& ab_group_per_replicate,
                                                                   const ublas::matrix<double>* max_count_assign_covariance)
{
    for(std::set<shared_ptr<ReadGroupProperties const > >::iterator itr = _read_group_props.begin();
        itr != _read_group_props.end(); 
        ++itr)
    {
        vector<shared_ptr<Abundance> > new_transcripts;
        BOOST_FOREACH(shared_ptr<Abundance> ab, transcripts)
        {
            boost::shared_ptr<TranscriptAbundance> d = boost::static_pointer_cast<TranscriptAbundance>(ab);
            //new_transcripts.push_back(shared_ptr<Abundance>(new TranscriptAbundance(*boost::static_pointer_cast<TranscriptAbundance>(ab))));
            TranscriptAbundance* pT = new TranscriptAbundance;
            pT->transfrag(d->transfrag());
            shared_ptr<Abundance> ab_new(pT);
            ab_new->description(ab_new->description());
            ab_new->locus_tag("");
            new_transcripts.push_back(ab_new);
        }
        shared_ptr<AbundanceGroup> ab_group(new AbundanceGroup(new_transcripts));
        std::set<shared_ptr<ReadGroupProperties const > > rg_props;
        rg_props.insert(*itr);
        ab_group->init_rg_props(rg_props);
        
        vector<MateHit> rep_hits;
        vector<double> rep_log_conv_factors;

        for (size_t i = 0; i < nr_alignments.size(); ++i)
        {
            if (nr_alignments[i].read_group_props() == *itr)
            {
                rep_hits.push_back(nr_alignments[i]);
                rep_log_conv_factors.push_back(log_conv_factors[i]);

            }
	
        }	
        //rep_hit_counts.push_back(count_per_replicate.find(*itr)->second);
        
        ab_group->calculate_abundance(rep_hits, false, false, max_count_assign_covariance);
        //fprintf (stderr, "FPKM = %lg\n", ab_group->FPKM());
        ab_group_per_replicate[*itr] = ab_group;
    }

    return NUMERIC_OK;
}

AbundanceStatus calculate_inverse_fisher(const vector<shared_ptr<Abundance> >& transcripts,
                                         const vector<MateHit>& alignments,
                                         const ublas::vector<double>& gamma_mean,
                                         ublas::matrix<double>& inverse_fisher)
{
//    size_t N = gamma_covariance.size1();	
	
//	gamma_map_covariance = ublas::zero_matrix<double>(N);
    
	typedef ublas::matrix<double> matrix_type;
	matrix_type fisher = ublas::zero_matrix<double>(gamma_mean.size(),gamma_mean.size());
	
	vector<double> u(alignments.size());
	for (size_t i = 0; i < alignments.size(); ++i)
	{
		u[i] = alignments[i].collapse_mass();
	}
    
	compute_fisher(transcripts,
				   gamma_mean,
				   alignments,
				   u,
				   fisher);
	
	ublas::matrix<double> epsilon = ublas::zero_matrix<double>(gamma_mean.size(),gamma_mean.size());
	for (size_t i = 0; i < gamma_mean.size(); ++i)
	{
		epsilon(i,i) = 1e-6;
	}
	
	fisher += epsilon; // modify matrix to avoid problems during inverse
    
	ublas::matrix<double> fisher_chol = fisher;
	
	double ch = cholesky_factorize(fisher_chol);
	if (ch != 0.0)
	{
		verbose_msg("Warning: Fisher matrix is not positive definite (bad element: %lg)\n", ch);
		return NUMERIC_FAIL;
	}
	
	inverse_fisher = ublas::zero_matrix<double>(gamma_mean.size(),gamma_mean.size());
	bool invertible = chol_invert_matrix(fisher_chol, inverse_fisher);
    
    ublas::matrix<double> test_fisher = inverse_fisher;
	ch = cholesky_factorize(test_fisher);
	if (ch != 0.0 || !invertible)
	{
		verbose_msg("Warning: Fisher matrix is not inverible\n", ch);
		return NUMERIC_FAIL;
	}
    
    return NUMERIC_OK;
}

AbundanceStatus revise_map_mean_and_cov_estimate(double log_total_weight,
                                                 const ublas::vector<double>& expectation,
                                                 const vector<pair<size_t, double> >& sample_weights,
                                                 const vector<ublas::vector<double> >& weighted_samples,
                                                 ublas::vector<double>& gamma_map_estimate,
                                                 ublas::matrix<double>& gamma_map_covariance)
{
    int N = expectation.size();
    
	// revise gamma by setting it to the posterior expectation computed via the
	// importance sampling
	gamma_map_estimate = expectation;
    
	// calculate the sample - mean vectors, store them in log space
	vector<ublas::vector<double> > sample_expectation_diffs;
	
    ublas::vector<double> check_expectation = ublas::zero_vector<double>(expectation.size());
    
	for (size_t j = 0; j < weighted_samples.size(); ++j)
	{
		ublas::vector<double> sample = weighted_samples[j];
		double log_sample_weight = sample_weights[j].second;
        
		for (size_t e = 0; e < expectation.size(); ++e)
		{
			// sample is already log transformed after it was weighted, so we
			// need to divide by the sample weight to recover the original sample
			// value, then undo the log transform, then subtract the mean from it
			sample(e) = (exp(((long double)sample(e) - log_sample_weight)) - expectation(e));
			//sample(e) *= exp((log_sample_weight - log_total_weight));
		}
		//cerr << sample << endl;
		sample_expectation_diffs.push_back(sample);
	}
     
    // We want to revise the covariance matrix from the samples, since we'll 
    // need it later for the CIs.
    ublas::matrix<double> revised_cov = ublas::zero_matrix<double>(N,N);
    
    // accumulate the contributions from the other samples (doing one cell of 
	// covariance matrix per outer (i x j) loop iteration.
    
    for (size_t j = 0; j < sample_expectation_diffs.size(); ++j)
    {
        double log_sample_weight = sample_weights[j].second;
        double w = exp((log_sample_weight - log_total_weight));
        ublas::vector<double> sample = weighted_samples[j];
        
        for (size_t e = 0; e < expectation.size(); ++e)
		{
			// sample is already log transformed after it was weighted, so we
			// need to divide by the sample weight to recover the original sample
			// value, then undo the log transform, then subtract the mean from it
			sample(e) = exp(sample(e) - log_sample_weight);
			//sample(e) *= exp((log_sample_weight - log_total_weight));
		}
        
        revised_cov += w * (outer_prod(sample,sample));
	}
    
    revised_cov -= outer_prod(expectation,expectation);
    
	//cerr << "Revised COV" << endl;
	//cerr << revised_cov << endl;
	gamma_map_covariance = revised_cov;
	
    //cerr << "Revised MAP estimate: " << expectation << endl;
    //cerr << "Revised Covariance matrix:" << endl;
    //cerr << gamma_map_covariance << endl;
    //cerr << "*************" << endl;
    
    return NUMERIC_OK;
}

AbundanceStatus calc_is_scale_factor(const ublas::matrix<double>& covariance_chol,
                                     double& is_scale_factor)
{
    double det = determinant(covariance_chol);
    is_scale_factor = pow(2.0*boost::math::constants::pi<double>(), covariance_chol.size1()/2.0);
	double s = sqrt(det);
	is_scale_factor *= s;
	
	//assert (det);
	if (s == 0.0)
	{
		verbose_msg("Error: sqrt(det(cov)) == 0, %lf after rounding. \n", det);
		//cerr << covariance << endl;
		return NUMERIC_FAIL;
	}
	assert (s);
	assert (is_scale_factor);
    return NUMERIC_OK;
}

template<class M, class PM>
bool is_identifiable(M &m, PM &pm)
{
    using namespace ublas;
    typedef M matrix_type;
    typedef typename M::size_type size_type;
    typedef typename M::value_type value_type;

    int singular = 0;
    size_type size1 = m.size1 ();
    size_type size2 = m.size2 ();
    size_type size = (std::min) (size1, size2);
    for (size_type i = 0; i < size; ++ i) {
        matrix_column<M> mci (column (m, i));
        matrix_row<M> mri (row (m, i));
        size_type i_norm_inf = i + index_norm_inf (project (mci, range (i, size1)));
        if (m (i_norm_inf, i) != value_type/*zero*/()) {
            if (i_norm_inf != i) {
                pm (i) = i_norm_inf;
                row (m, i_norm_inf).swap (mri);
            } else {
                //BOOST_UBLAS_CHECK (pm (i) == i_norm_inf, external_logic ());
            }
            project (mci, range (i + 1, size1)) *= value_type (1) / m (i, i);
        } else if (singular == 0) {
            singular = i + 1;
        }
        project (m, range (i + 1, size1), range (i + 1, size2)).minus_assign (outer_prod (project (mci, range (i + 1, size1)),
                                                                              project (mri, range (i + 1, size2))));
    }
    return singular == 0;
}

AbundanceStatus gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                          const vector<MateHit>& nr_alignments,
                          const vector<double>& log_conv_factors,
                          vector<double>& gammas,
                          bool check_identifiability)
{
	//gammas.clear();
    gammas = vector<double>(transcripts.size(), 0);
	if (transcripts.empty())
		return NUMERIC_OK;
	
	//long double bundle_mass_fraction = bundle_mass / (long double) map_mass;
	if (transcripts.size() == 1)
	{
		gammas = vector<double>(1,1.0);
		return NUMERIC_OK;
	}
	
	size_t M = nr_alignments.size();
	size_t N = transcripts.size();

    bool converged = true;
    bool identifiable = true;
    
	if (M > 0)
	{
		
		//vector<vector<double> > saliencies (M,vector<double>(N,0));
		
		
		//compute_saliencies(cond_probs, saliencies, saliency_weight);
		
        Eigen::VectorXd prob = Eigen::VectorXd::Zero(N);
        Eigen::VectorXd _phint = Eigen::VectorXd::Zero(N);
        
		double logL;
		
        Eigen::MatrixXd cond_probs(N, M);
		
        for (size_t j = 0; j < cond_probs.rows(); ++j)
        {
            for (size_t i = 0; i < cond_probs.cols(); ++i)
            {
                cond_probs(j,i) = (*(transcripts[j]->cond_probs()))[i];
            }
        }
        
        if (check_identifiability)
        {
            ublas::matrix<double> compat = ublas::zero_matrix<double>(M,N);
            
            for (size_t j = 0; j < N; ++j)
            {
                for (size_t i = 0; i < M; ++i)
                {
                    if (cond_probs(j,i))
                    {
                        //compat(i,j) = cond_probs[j][i];
                        compat(i,j) = 1;
                    }
                }
            }
            
            vector<size_t> transcripts_with_frags;
            for (size_t j = 0; j < N; ++j)
            {
                bool has_fragment = false;
                for (size_t i = 0; i < M; ++i)
                {
                    if (compat(i,j))
                    {
                        has_fragment = true;
                        break;
                    }
                }
                if (has_fragment)
                    transcripts_with_frags.push_back(j);
            }
            ublas::matrix<double> reduced_compat = ublas::zero_matrix<double>(M,transcripts_with_frags.size());
            for (size_t j = 0; j < transcripts_with_frags.size(); ++j)
            {
                column(reduced_compat, j) = column(compat, transcripts_with_frags[j]);
            }
            
            
            typedef ublas::permutation_matrix<std::size_t> pmatrix;
            
            // create a permutation matrix for the LU-factorization
            pmatrix pm(reduced_compat.size1());
            
            // cerr << compat.size2() <<endl;
            // perform LU-factorization
            identifiable = is_identifiable<ublas::matrix<double>,pmatrix>(reduced_compat,pm);
        }
        
        Eigen::VectorXd alignment_multiplicities(M);
		for (size_t i = 0; i < M; ++i)
		{
			alignment_multiplicities[i] = nr_alignments[i].collapse_mass();
		}
        		
        logL = EM(N, M, prob, cond_probs, alignment_multiplicities, log_conv_factors, converged);
        
		gammas = vector<double>(N, 0.0);
		
		for (size_t i = 0; i < gammas.size(); ++i)
		{
            gammas[i] = prob(i);
			if (isnan(gammas[i]) || isinf(gammas[i]))
            {
                return NUMERIC_FAIL;
            }
		}
	}
	else
	{
		gammas = vector<double>(N, 0.0);
	}
    
    double round_err = 0.0;
    double num_good = 0;
    BOOST_FOREACH (double& g, gammas)
    {
        if (g < min_isoform_fraction)
        {
            round_err += g;
            g = 0.0;
        }
        else
        {
            num_good += 1;
        }
    }
    BOOST_FOREACH (double& g, gammas)
    {
        if (g != 0)
        {
            g += (round_err/num_good);
        }
    }
    
    if (converged && identifiable)
        return NUMERIC_OK;
    else
    {
        if (!identifiable)
            //return NUMERIC_LOW_DATA;
            return NUMERIC_OK;
        else
            return NUMERIC_FAIL;
    }
    
    return NUMERIC_OK;
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

bool not_intronic(int p, vector<float>& depth_of_coverage, vector<float>& intronic_cov, float min_intra_intron_fraction,
      int& intronic_status) {
  bool not_an_intron = (intronic_cov[p]==0 ||
        depth_of_coverage[p]/intronic_cov[p] >= min_intra_intron_fraction);
 if (not_an_intron) intronic_status--;
               else intronic_status++;
 return not_an_intron;
}


double compute_doc(int bundle_origin, 
				   const vector<Scaffold>& scaffolds,
				   vector<float>& depth_of_coverage,
				   map<pair<int, int>, float>& intron_depth_of_coverage,
				   bool exclude_intra_intron,
				   vector<float>* intronic_cov,
				   vector<int>* scaff_intronic_status)
{
  vector<int> i_status;
  if (scaff_intronic_status==NULL)
    scaff_intronic_status=&i_status;
  *scaff_intronic_status = vector<int>(scaffolds.size(), 0);
  vector<float> intronic;
  if (intronic_cov==NULL)
    intronic_cov=&intronic;
  *intronic_cov = vector<float>(depth_of_coverage.size(), 0);
  //vector<bool> intronic(depth_of_coverage.size(), false);
	depth_of_coverage = vector<float>(depth_of_coverage.size(), 0);
		
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

	/*
	//no need for this here, we do it below with depth_of_coverage
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
	*/
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
        for (int K = op.g_left(); K < op.g_right(); ++K)
        {
          (*intronic_cov)[K - bundle_origin] += hits[i].fpkm();
          //intronic[K - bundle_origin] = true;
        }

				pair<map<pair<int,int>,float>::iterator, bool> is = intron_depth_of_coverage.insert(make_pair(make_pair(op.g_left(), op.g_right()), 0));
				is.first->second += hits[i].fpkm();
			}
		}
	}
	
	vector<float> knockout(depth_of_coverage);
	
	double total_doc = 0;
	int total_len = 0;
	float min_intra_intron_fraction = min(pre_mrna_fraction, min_isoform_fraction);
	//for (size_t i = 0; i < hits.size(); ++i)
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		//const vector<AugmentedCuffOp>& aug_ops = hits[i].augmented_ops();
	  const vector<AugmentedCuffOp>& aug_ops = scaffolds[i].augmented_ops();
		for (size_t j = 0; j < aug_ops.size(); ++j)
		{
			const AugmentedCuffOp& op = aug_ops[j];
			if (op.opcode == CUFF_MATCH)
			{
				for (int K = op.g_left(); K < op.g_right(); ++K)
				{
					//if (!exclude_intra_intron || !intronic[K - bundle_origin])
				  if (!exclude_intra_intron ||
				       not_intronic(K-bundle_origin, depth_of_coverage, *intronic_cov, min_intra_intron_fraction,
				           (*scaff_intronic_status)[i]) )
					{
						total_doc += knockout[K - bundle_origin];
						total_len += (knockout[K - bundle_origin] != 0);
						knockout[K - bundle_origin] = 0;
					}
				}
			}
		}
	}
	
	return total_doc/(double)total_len;
}

double major_isoform_intron_doc(map<pair<int, int>, float>& intron_doc)
{
	double major_isoform_intron_doc = 0;
	int num_major_introns = 0;
	for(map<pair<int, int>, float>::const_iterator itr = intron_doc.begin();
		itr != intron_doc.end(); 
		++itr)
	{
		bool heaviest = true;
		
		for (map<pair<int,int>, float>::const_iterator itr2 = intron_doc.begin();
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
								  const vector<float>& depth_of_coverage,
								  const map<pair<int, int>, float>& intron_depth_of_coverage,
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
							  const vector<float>& depth_of_coverage,
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
							  const vector<float>& depth_of_coverage,
							  const map<pair<int, int>, float>& intron_depth_of_coverage,
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
					  const map<pair<int, int>, float >& intron_depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	int num_introns = 0;
	double doc = 0;
	for (size_t j = 0; j < aug_ops.size(); ++j)
	{
		const AugmentedCuffOp& op = aug_ops[j];
		if (op.opcode == CUFF_INTRON)
		{
			num_introns++;
			pair<int,int> op_intron(op.g_left(), op.g_right());
			map<pair<int, int>, float >::const_iterator itr = intron_depth_of_coverage.find(op_intron);
			//			assert (itr != intron_depth_of_coverage.end());
			if (itr == intron_depth_of_coverage.end())
			{
				map<pair<int, int>, float >::const_iterator zi;
				for (zi = intron_depth_of_coverage.begin();
					 zi != intron_depth_of_coverage.end();
					 ++zi)
				{
					verbose_msg( "Warning: intron not within scaffold ([%d-%d], %d)\n", zi->first.first, zi->first.second, zi->second);
				}
			}

			doc += itr->second;
		}
	}
	return doc / (double)num_introns;
}

double get_scaffold_doc(int bundle_origin, 
						const Scaffold& s,
						const vector<float>& depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	int m_len = 0;
	double doc = 0;
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
	
	return doc/(double)m_len;
}

double get_scaffold_min_doc(int bundle_origin, 
							const Scaffold& s,
							const vector<float>& depth_of_coverage)
{
	const vector<AugmentedCuffOp>& aug_ops = s.augmented_ops();
	float min_doc = 99999999;
	
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

bool generate_null_js_samples(const vector<Eigen::VectorXd>& rel_abundances,
                              size_t num_js_samples,
                              vector<double>& js_samples)
{
    if (rel_abundances.empty())
        return true;
    
    size_t num_abundances = rel_abundances.front().size();
    
    if (num_abundances <= 1)
        return true;
    //    ublas::vector<double> null_kappa_mean(num_abundances);
    //    for (size_t i = 0; i < num_abundances; ++i)
    //    {
    //        null_kappa_mean(i) = rel_abundances[i];
    //    }
    
    //    cerr << endl << null_kappa_mean << endl;
    //    for (unsigned i = 0; i < null_kappa_cov.size1 (); ++ i) 
    //    {
    //        ublas::matrix_row<ublas::matrix<double> > mr (null_kappa_cov, i);
    //        std::cerr << i << " : " << mr << std::endl;
    //    }
    //    cerr << "======" << endl;
    
    //vector<ublas::vector<double> > null_samples;
    
    js_samples.clear();
    
    //size_t num_samples = std::min(prev_samples.size(), curr_samples.size());
    size_t num_samples = num_js_samples;
    vector<Eigen::VectorXd> sample_kappas(2);
    
    boost::uniform_int<> null_uniform_dist(0,rel_abundances.size()-1);
    boost::mt19937 null_rng; 
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > null_uniform_gen(null_rng, null_uniform_dist); 
    
    for (size_t i = 0; i < num_samples; ++i)
    {
		sample_kappas[0] = rel_abundances[null_uniform_gen()];
        sample_kappas[1] = rel_abundances[null_uniform_gen()];
        
		double js = jensen_shannon_distance(sample_kappas);  
        assert(!isnan(js));
        //cerr << sample_kappas[0].transpose() << " vs. " <<  sample_kappas[1].transpose() << " = " << js << endl;
        js_samples.push_back(js);
    }
    
    sort(js_samples.begin(), js_samples.end());
    
    return true;
}

