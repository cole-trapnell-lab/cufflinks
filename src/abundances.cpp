/*
 *  abundances.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 4/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 *  NOTE: some of the code in this file was derived from (Eriksson et al, 2008)
 */


//#define BOOST_MATH_INSTRUMENT 1
//#include <iostream>
//#include <iomanip>

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
#include <boost/math/tools/tuple.hpp>

#include <boost/math/distributions/chi_squared.hpp>

#include "filters.h"
#include "replicates.h"
#include "sampling.h"
#include "jensen_shannon.h"
#include "rounding.h"
#include "clustering.h"


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

// Returns the log likelihood of the negative binomial with a given r and the mle value of p
long double negbin_log_likelihood(const vector<double>& samples, long double r, long double p)
{
    if (samples.empty() || r == 0 || p == 0)
    {
        return 1.0;
    }
    
    long double T1 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T1 += lgamma(samples[i] + r);
    }
    
    long double T2 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T2 += lgamma(samples[i] + 1);
    }
    
    long double T3 = samples.size() * lgamma(r);
    long double T4 = samples.size() * r * logl(1 - p);
    
    long double T5 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T5 += samples[i] * logl(p);
    }
    
    assert ((isnan(T1) || isnan(T2) || isnan(T3) || isnan(T4) || isnan(T5)) == false);
    long double log_likelihood = T1 - T2 - T3 + T4 + T5;
    return log_likelihood;
}


// Returns the log likelihood of the negative binomial with a given r and the mle value of p
long double poisson_log_likelihood(const vector<double>& samples, long double lambda)
{
    if (samples.empty() || lambda == 0.0)
    {
        return 1.0;
    }
    
    long double T1 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T1 += samples[i] * logl(lambda);
    }
    
    long double T2 = samples.size() * lambda;
        
    assert ((isnan(T1) || isnan(T2)) == false);
    long double log_likelihood = T1 - T2;
    return log_likelihood;
}



// Returns the log likelihood of the negative binomial with a given r and the mle value of p
long double negbin_log_likelihood_helper(const vector<double>& samples, long double r)
{
    long double p = 0;
    long double mean_count = accumulate(samples.begin(), samples.end(), 0.0);
    
    if (r == 0 || mean_count == 0)
    {
        fprintf(stderr, "Error: r must be > 0!\n");
        return 0;
    }
    
    if (samples.empty() == false)
        mean_count /= samples.size();
    
    p = mean_count / (r + mean_count);
    
    long double T1 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T1 += lgamma(samples[i] + r);
    }

    long double T2 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T2 += lgamma(samples[i] + 1);
    }

    long double T3 = samples.size() * lgamma(r);
    long double T4 = samples.size() * r * log(1 - p);
    
    long double T5 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T5 += samples[i] * log(p);
    }
    
    assert ((isnan(T1) || isnan(T2) || isnan(T3) || isnan(T4) || isnan(T5)) == false);
    long double log_likelihood = T1 - T2 - T3 + T4 + T5;
    return log_likelihood;
}


// Returns the log likelihood of the negative binomial with a given r and the mle value of p
long double negbin_log_likelihood_prime_helper(const vector<double>& samples, long double r)
{
    long double T1 = 0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        T1 += digamma(samples[i] + r);
    }
    
    long double T2 = samples.size() * digamma(r);
    
    long double mean_count = accumulate(samples.begin(), samples.end(), 0.0);
    if (samples.empty() == false)
        mean_count /= samples.size();

    long double T3 = samples.size() * log(r / (r + mean_count));
    
    assert ((isnan(T1) || isnan(T2) || isnan(T3)) == false);
    
    long double log_likelihood_prime = T1 - T2 + T3;
    return log_likelihood_prime;
}


struct negbin_ll_functor
{
    negbin_ll_functor(const vector<double>& count_samples) : samples(count_samples){}
    boost::math::tuple<long double, long double> operator()(long double r)
    {
        long double llk = negbin_log_likelihood_helper(samples, r);
        long double llk_d = negbin_log_likelihood_prime_helper(samples, r);
        return boost::math::make_tuple(llk, llk_d);
    }
private:
    const vector<double>& samples;
};


bool fit_negbin_dist(const vector<double> samples, double& r, double& p)
{
    
    if (samples.size() == 0)
    {
        r = 0;
        p = 0;
        return true;
    }
    
    long double guess = accumulate(samples.begin(), samples.end(), 0.0); 
    if (samples.empty() == false)
        guess /= samples.size();
    
    if (guess == 0)
    {
        r = 0;
        p = 0;
        return true;
    }
    
    long double min = 0;
    long double max = *std::max_element(samples.begin(), samples.end());
    max *= max;
    int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
    
    boost::uintmax_t max_iters = 50;
    
    r = boost::math::tools::newton_raphson_iterate(negbin_ll_functor(samples), guess, min, max, digits, max_iters);
    
    long double mean_count = accumulate(samples.begin(), samples.end(), 0.0);
    if (samples.empty() == false)
        mean_count /= samples.size();
    
    p = mean_count / (r + mean_count);
    
    //fprintf(stderr, "r = %lg, p = %lg, max_r = %Lg\n", r, p, max);
    
    if (isnan(r) || isnan(p))
    {
        fprintf(stderr, "warning: negative binomial parameters are Nan!\n");
        return false;
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
                               const set<shared_ptr<ReadGroupProperties const> >& rg_props) :
    _abundances(abundances), 
    _iterated_exp_count_covariance(iterated_exp_count_covariance),
    _count_covariance(count_covariance),
    _fpkm_covariance(fpkm_covariance),
    _gamma_covariance(gamma_covariance),
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
        {
            _fpkm_samples = vector<double>(_abundances[i]->fpkm_samples().size(), 0);
            _member_fpkm_samples = vector<Eigen::VectorXd>(_abundances[i]->fpkm_samples().size(), Eigen::VectorXd::Zero(_abundances.size()));
        }
        for (size_t j = 0; j < _abundances[i]->fpkm_samples().size(); ++j)
        {
            _fpkm_samples[j] += _abundances[i]->fpkm_samples()[j];
            _member_fpkm_samples[j](i) = _abundances[i]->fpkm_samples()[j];
        }
    }
    
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

void TranscriptAbundance::clear_non_serialized_data()
{
    _fpkm_samples.clear();
    std::vector<double>().swap(_fpkm_samples);
    
    if (_cond_probs)
    {
        _cond_probs->clear();
        std::vector<double>().swap(*_cond_probs);
    }
    
    if (_transfrag)
    {
        _transfrag->clear_hits();
        _transfrag = shared_ptr<Scaffold>();
    }
}

void AbundanceGroup::clear_non_serialized_data()
{
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        _abundances[i]->clear_non_serialized_data();
    }
    
    _fpkm_samples.clear();
    std::vector<double>().swap(_fpkm_samples);
    _member_fpkm_samples.clear();
    std::vector<Eigen::VectorXd>().swap(_member_fpkm_samples);
    _assigned_count_samples.clear();
    std::vector<Eigen::VectorXd>().swap(_assigned_count_samples);
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
    
	//vector<vector<double> > new_fpkm_samples(_fpkm_samples.size(), vector<double>(num_kept, 0));
	
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
	//assert (false);
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
	//assert (false);
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
            double curr_mass = inserted.first->second;
            assert (isnan(more_mass) == false);
            inserted.first->second += more_mass;
        }
    }
}

void AbundanceGroup::calculate_locus_scaled_mass_and_variance(const vector<shared_ptr<Abundance> >& transcripts)
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
        
        //assert (scaled_total_mass != 0.0);
        avg_mass_fraction += (scaled_mass / scaled_total_mass);
    }
    
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
    
    for (size_t j = 0; j < _abundances.size(); ++j)
	{
		_abundances[j]->num_fragments(_abundances[j]->gamma() * avg_X_g);
        
        double j_avg_mass_fraction = _abundances[j]->gamma() * avg_mass_fraction;
        
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
        }
	}

}

int total_cond_prob_calls = 0;
void collapse_equivalent_hits(const vector<MateHit>& alignments,
                              vector<shared_ptr<Abundance> >& transcripts,
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
    
    double total_mass_pre_collapse = 0.0;
    
    for(int i = 0 ; i < M; ++i)
    {
        total_mass_pre_collapse += alignments[i].collapse_mass();
        
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
	}
    
    double total_mass_post_collapse = 0.0;
    for(int i = 0 ; i < nr_alignments.size(); ++i)
    {
        total_mass_post_collapse += nr_alignments[i].collapse_mass();
    }
    
    //assert(abs(total_mass_pre_collapse - total_mass_post_collapse) < 1);
    
    if (nr_alignments.size())
    {
        verbose_msg("\nReduced %lu frags to %lu (%lf percent)\n", alignments.size(), nr_alignments.size(), 100.0 * (1 - nr_alignments.size()/(double)alignments.size()));
    }
}

void collapse_equivalent_hits_helper(const vector<MateHit>& alignments,
                                     vector<shared_ptr<Abundance> >& transcripts,
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
        
        for (size_t j = 0; j < compat_table[i].size(); ++j)
        {
            tmp_hits.push_back(*(compat_table[i][j]));
        }
        if (tmp_hits.empty())
            continue;
        collapse_equivalent_hits(tmp_hits,
                                 transcripts,
                                 tmp_nr_hits,
                                 tmp_log_conv_factors, 
                                 false);
        copy(tmp_nr_hits.begin(), tmp_nr_hits.end(), back_inserter(nr_alignments));
        copy(tmp_log_conv_factors.begin(), tmp_log_conv_factors.end(), back_inserter(log_conv_factors));
    }
}


AbundanceStatus bootstrap_gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                                    const vector<MateHit>& nr_alignments,
                                    const vector<double>& log_conv_factors,
                                    ublas::vector<double>& gamma_map_estimate,
                                    ublas::matrix<double>& gamma_covariance,
                                    double& cross_replicate_js)
{
    size_t N = transcripts.size();
	size_t M = nr_alignments.size();
    
    if (M == 0)
    {
        gamma_map_estimate = ublas::vector<double>(1);
        gamma_map_estimate(0) = 0.0;
        gamma_covariance = ublas::matrix<double>(1,1);
        gamma_covariance(0,0) = 0.0;
        return NUMERIC_OK;
    }
    
    if (N == 1)
    {
        gamma_map_estimate = ublas::vector<double>(1);
        gamma_map_estimate(0) = 1.0;
        gamma_covariance = ublas::matrix<double>(1,1);
        gamma_covariance(0,0) = 0.0;
        return NUMERIC_OK;
    }
    
    vector<MateHit> alignments = nr_alignments;
    vector<double>  scaled_masses;
    vector<double>  unscaled_masses;
    double num_uncollapsed_frags = 0.0;
    for (size_t i = 0; i < M; ++i)
    {
        double uncollapsed_mass = alignments[i].collapse_mass();
        num_uncollapsed_frags += (uncollapsed_mass);
        scaled_masses.push_back(alignments[i].collapse_mass());
        unscaled_masses.push_back(uncollapsed_mass);
        alignments[i].collapse_mass(uncollapsed_mass);
    }
    
    // FIXME: this has already been computed above, so just pass it in.
    vector<double> orig_gammas(0.0, transcripts.size());
    gamma_mle(transcripts,
              nr_alignments,
              log_conv_factors,
              orig_gammas,
              false);
    
    std::vector<ublas::vector<double> > mle_gammas;
    
    boost::uniform_int<> uniform_dist(0,num_uncollapsed_frags-1);
    boost::mt19937 rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniform_gen(rng, uniform_dist);
    
    int num_sample_frags = floor(num_uncollapsed_frags * bootstrap_fraction);
    
    if (num_sample_frags <= 0)
    {
        return NUMERIC_FAIL;
    }
    
    for (size_t i = 0; i < num_bootstrap_samples; ++i)
    {
        vector<int> sample_idxs;
        for (size_t j = 0; j < num_sample_frags; ++j)
        {
            sample_idxs.push_back(uniform_gen());
        }
        sort (sample_idxs.begin(), sample_idxs.end());
        assert (sample_idxs.empty() == false);
        
        size_t curr_sample = 0;
        size_t processed_hits = 0;
        vector<double> adjusted_masses(alignments.size(), 0);
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            int adjusted_mass = 0.0;
            while (curr_sample < sample_idxs.size() &&
                   sample_idxs[curr_sample] >= processed_hits &&
                   sample_idxs[curr_sample] < processed_hits + alignments[j].collapse_mass())
            {
                adjusted_mass++;
                curr_sample++;
            }
            processed_hits += alignments[j].collapse_mass();
            alignments[j].collapse_mass(adjusted_mass);
            adjusted_masses[j] = adjusted_mass;
        }
        
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            alignments[j].collapse_mass(alignments[j].collapse_mass());
        }
        
        vector<double> bs_gammas(0.0, transcripts.size());
        
//        AbundanceStatus mle_success = gamma_mle(transcripts,
//                                                alignments,
//                                                log_conv_factors,
//                                                bs_gammas,
//                                                false,
//                                                &orig_gammas);
        
        AbundanceStatus mle_success = gamma_mle(transcripts,
                                                alignments,
                                                log_conv_factors,
                                                bs_gammas,
                                                false);
        if (mle_success == NUMERIC_OK)
        {
            ublas::vector<double> mle = ublas::zero_vector<double>(N);
            for(size_t j = 0; j < N; ++j)
            {
                mle(j) = bs_gammas[j];
            }
            mle_gammas.push_back(mle);
        }
        
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            alignments[j].collapse_mass(unscaled_masses[j]);
        }
    }
    
    //fprintf(stderr, "Ran %lu bootstrap samples succesfully\n", mle_gammas.size());
    
    if (mle_gammas.empty())
        return NUMERIC_FAIL;
    
    gamma_covariance = ublas::zero_matrix<double>(N,N);
    ublas::vector<double> expected_mle_gamma = ublas::zero_vector<double>(N);
    
    BOOST_FOREACH(ublas::vector<double>& mle, mle_gammas)
    {
        //cerr << "MLE # "<< MLENUM++ << endl;
        //cerr << mle << endl;
        expected_mle_gamma += mle;
    }
    expected_mle_gamma /= mle_gammas.size();
    
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            for (size_t k = 0 ; k < mle_gammas.size(); ++k)
            {
                double c = (mle_gammas[k](i) - expected_mle_gamma(i)) * (mle_gammas[k](j) - expected_mle_gamma(j));
                gamma_covariance(i,j) += c;
            }
        }
    }
    
    gamma_covariance /= mle_gammas.size();
    gamma_map_estimate = expected_mle_gamma;
    
    //cerr << "MLE: " << expected_mle_gamma << endl;
    //cerr << "COV:" << endl;
    //cerr << gamma_covariance << endl;
    //cerr << "*************" << endl;
    return NUMERIC_OK;
}



AbundanceStatus bootstrap_gammas(const vector<shared_ptr<Abundance> >& transcripts,
                                 const vector<MateHit>& alignments,
                                 const vector<double>& log_conv_factors,
                                 ublas::vector<double>& gamma_estimate,
                                 ublas::matrix<double>& gamma_covariance,
                                 double& cross_replicate_js)
{
    ublas::vector<double> empirical_gamma_mle = gamma_estimate;
    ublas::matrix<double> empirical_gamma_covariance = gamma_covariance;
    
    // Calculate the mean gamma MLE and covariance matrix across replicates, so
    // we can use it as the proposal distribution for importance sampling.  This will
    // make the Bayesian prior more conservative than using the inverse of the
    // Fisher Information matrix on the mixed likelihood function.
    AbundanceStatus empirical_mle_status = bootstrap_gamma_mle(transcripts,
                                                               alignments,
                                                               log_conv_factors,
                                                               empirical_gamma_mle,
                                                               empirical_gamma_covariance,
                                                               cross_replicate_js);
    
    if (empirical_mle_status != NUMERIC_OK)
        return empirical_mle_status;
    
    gamma_estimate = empirical_gamma_mle;
    gamma_covariance = empirical_gamma_covariance;
    
    
    
    
    
    return NUMERIC_OK;
}


 // Sample from a multivariate normal to estimate the gamma distribution
bool generate_count_assignment_samples(int num_draws,
                                       const vector<double>& count_mean,
                                       const ublas::matrix<double>& count_covariance,
                                       vector<ublas::vector<double> >& assigned_count_samples)
{
    double total_frag_counts = accumulate(count_mean.begin(), count_mean.end(), 0.0);
    
    if (total_frag_counts == 0)
    {
        assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(count_mean.size()));
        return true;
    }
    
    boost::mt19937 rng;
    
    ublas::vector<double> mle_frag_counts = ublas::zero_vector<double>(count_mean.size());
    
    for (size_t j = 0; j < count_mean.size(); ++j)
    {
        mle_frag_counts(j) = count_mean[j];
    }
    
    //    cerr << "********************" << endl;
    //    cerr << "initial MLE counts: " << mle_frag_counts << endl;
    
    ublas::matrix<double> mle_count_covar = count_covariance;
    
    ublas::matrix<double> epsilon = ublas::zero_matrix<double>(count_mean.size(),count_mean.size());
    for (size_t i = 0; i < count_mean.size(); ++i)
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
    
    assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(count_mean.size()));
    
    boost::uniform_01<> uniform_dist;
    boost::mt19937 null_rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uniform_gen(null_rng, uniform_dist);
    
    
    for (size_t assign_idx = 0; assign_idx < num_draws; ++assign_idx)
    {
        
        boost::numeric::ublas::vector<double> random_count_assign = generator.next_rand();
        //cerr << random_count_assign << endl;
        
        for (size_t r_idx = 0; r_idx < random_count_assign.size(); ++r_idx)
        {
            if (random_count_assign(r_idx) < 0)
                random_count_assign(r_idx) = 0;
        }
        
        double total_sample_counts = accumulate(random_count_assign.begin(), random_count_assign.end(), 0.0);
        if (total_sample_counts > 0)
            random_count_assign = total_frag_counts * (random_count_assign / total_sample_counts);
        else
            random_count_assign = ublas::zero_vector<double>(count_mean.size());
        assigned_count_samples[assign_idx] = random_count_assign;
    }
    
    return true;
}
 
void calculate_gamma_mle_covariance(const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >& ab_group_per_replicate,
                                    ublas::vector<double>& estimated_gamma_mean,
                                    ublas::matrix<double>& estimated_gamma_covariance)
{
    vector<ublas::vector<double> > all_assigned_count_samples;
    
    if (ab_group_per_replicate.empty())
        return;
    int num_transcripts = ab_group_per_replicate.begin()->second->abundances().size();
    if (num_transcripts == 1)
    {
        estimated_gamma_mean = ublas::vector<double>(1);
        estimated_gamma_mean(0) = 1.0;
        
        estimated_gamma_covariance = ublas::matrix<double>(1,1);
        estimated_gamma_covariance(0,0) = 0.0;
        return;
    }
    
    for(std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
        itr != ab_group_per_replicate.end();
        ++itr)
    {
        ublas::vector<double> count_mean = ublas::zero_vector<double>(itr->second->abundances().size());
        for (size_t i = 0; i < itr->second->abundances().size(); ++i)
        {
            count_mean(i) = itr->second->abundances()[i]->gamma();
        }
        
        all_assigned_count_samples.push_back(count_mean);
    }
    
    //    for (size_t i = 0; i < all_assigned_count_samples.size(); ++i)
    //    {
    //        double total = accumulate(all_assigned_count_samples[i].begin(), all_assigned_count_samples[i].end(), 0.0);
    //        if (total > 0)
    //            all_assigned_count_samples[i] /= total;
    //    }
    
    estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
    estimated_gamma_covariance = ublas::zero_matrix<double>(num_transcripts, num_transcripts);
    
    for (size_t i = 0; i < all_assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < all_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(all_assigned_count_samples[i](j)) && !isinf(all_assigned_count_samples[i](j)));
        }
        
        estimated_gamma_mean += all_assigned_count_samples[i];
        //
        //expected_relative_abundances += relative_abundances[i];
    }
    
    if (all_assigned_count_samples.size() > 1)
    {
        estimated_gamma_mean /= (all_assigned_count_samples.size() - 1);
    }
    
    for (size_t i = 0; i < num_transcripts; ++i)
    {
        for (size_t j = 0; j < num_transcripts; ++j)
        {
            for (size_t k = 0 ; k < all_assigned_count_samples.size(); ++k)
            {
                double c = (all_assigned_count_samples[k](i) - estimated_gamma_mean(i)) * (all_assigned_count_samples[k](j) - estimated_gamma_mean(j));
                estimated_gamma_covariance(i,j) += c;
                
                assert (!isinf(estimated_gamma_covariance(i,j)) && !isnan(estimated_gamma_covariance(i,j)));
            }
        }
    }
    
    if (all_assigned_count_samples.size() > 1)
    {
        estimated_gamma_covariance /= (all_assigned_count_samples.size() - 1);
    }
}

void calculate_fragment_assignment_distribution(const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >& ab_group_per_replicate,
                                                ublas::vector<double>& estimated_gamma_mean,
                                                ublas::matrix<double>& estimated_gamma_covariance,
                                                vector<ublas::vector<double> >& all_assigned_count_samples)
{
    
    all_assigned_count_samples.clear();
    
    if (ab_group_per_replicate.empty())
        return;
    int num_transcripts = ab_group_per_replicate.begin()->second->abundances().size();
    if (num_transcripts == 1)
    {
        estimated_gamma_mean = ublas::vector<double>(1);
        estimated_gamma_mean(0) = 1.0;
        
        estimated_gamma_covariance = ublas::matrix<double>(1,1);
        estimated_gamma_covariance(0,0) = 0.0;
        return;
    }
    
    for(std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
        itr != ab_group_per_replicate.end();
        ++itr)
    {
        vector<double> count_mean;
        for (size_t i = 0; i < itr->second->abundances().size(); ++i)
        {
            count_mean.push_back(itr->second->abundances()[i]->num_fragments());
        }
        
        ublas::matrix<double> count_covariance = itr->second->iterated_count_cov();
        ublas::matrix<double> mle_error = ublas::zero_matrix<double>(count_mean.size(), count_mean.size());
        shared_ptr<const MleErrorModel> mle_model = itr->first->mle_error_model();
        if (mle_model != NULL)
        {
            for (size_t i = 0; i < count_mean.size(); ++i)
            {
                double mle_var = mle_model->scale_mle_variance(count_mean[i]);
                mle_error(i,i) = max(0.0, mle_var);
            }
            count_covariance += mle_error;
//            cerr << endl << "MLE error correction: " << endl;
//            for (unsigned i = 0; i < mle_error.size1 (); ++ i)
//            {
//                ublas::matrix_row<ublas::matrix<double> > mr (mle_error, i);
//                cerr << i << " : " << count_mean[i] << " : ";
//                std::cerr << i << " : " << mr << std::endl;
//            }
//            cerr << "======" << endl;
        }
        
        
        vector<ublas::vector<double> > assigned_count_samples;
        generate_count_assignment_samples(num_frag_assignments,
                                          count_mean,
                                          count_covariance,
                                          assigned_count_samples);

        all_assigned_count_samples.insert(all_assigned_count_samples.end(), assigned_count_samples.begin(), assigned_count_samples.end());
    }
    
    for (size_t i = 0; i < all_assigned_count_samples.size(); ++i)
    {
        double total = accumulate(all_assigned_count_samples[i].begin(), all_assigned_count_samples[i].end(), 0.0);
        if (total > 0)
            all_assigned_count_samples[i] /= total;
        //cerr << all_assigned_count_samples[i] << endl;
    }
    
    estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
    estimated_gamma_covariance = ublas::zero_matrix<double>(num_transcripts, num_transcripts);
    
    for (size_t i = 0; i < all_assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < all_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(all_assigned_count_samples[i](j)) && !isinf(all_assigned_count_samples[i](j)));
        }
        
        estimated_gamma_mean += all_assigned_count_samples[i];
        //
        //expected_relative_abundances += relative_abundances[i];
    }
    
    if (all_assigned_count_samples.size() > 1)
    {
        estimated_gamma_mean /= (all_assigned_count_samples.size() - 1);
    }
       
    for (size_t i = 0; i < num_transcripts; ++i)
    {
        for (size_t j = 0; j < num_transcripts; ++j)
        {
            for (size_t k = 0 ; k < all_assigned_count_samples.size(); ++k)
            {
                double c = (all_assigned_count_samples[k](i) - estimated_gamma_mean(i)) * (all_assigned_count_samples[k](j) - estimated_gamma_mean(j));
                estimated_gamma_covariance(i,j) += c;
                
                assert (!isinf(estimated_gamma_covariance(i,j)) && !isnan(estimated_gamma_covariance(i,j)));
            }
        }
    }
    
    if (all_assigned_count_samples.size() > 1)
    {
        estimated_gamma_covariance /= (all_assigned_count_samples.size() - 1);
    }
    
//    cerr << "Gamma covariance: " << endl;
//    for (unsigned i = 0; i < estimated_gamma_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (estimated_gamma_covariance, i);
//        cerr << i << " : " << estimated_gamma_mean[i] << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }

}


#define PERFORM_EQUIV_COLLAPSE 1

//void AbundanceGroup::calculate_abundance_for_replicate(const vector<MateHit>& alignments)
//{
//    
//}


void AbundanceGroup::calculate_abundance_group_variance(const vector<shared_ptr<Abundance> >& transcripts,
                                                        const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >& ab_group_per_replicate)
{
	if (final_est_run) // Only on last estimation run
	{
        // Simulate NB draws and fragment assignment under uncertainty to sample
        // from the BNBs.
        
        ublas::matrix<double> count_assign_covariance;
        
        ublas::vector<double> estimated_gamma_mean;
        ublas::matrix<double> estimated_gamma_covariance;
        
        vector<ublas::vector<double> > gamma_samples;
        calculate_fragment_assignment_distribution(ab_group_per_replicate,
                                                   estimated_gamma_mean,
                                                   estimated_gamma_covariance,
                                                   gamma_samples);
        ublas::vector<double> estimated_count_mean = estimated_gamma_mean * num_fragments();
        ublas::matrix<double>  estimated_count_covariance = estimated_gamma_covariance * num_fragments() * num_fragments();
        
        //cerr << estimated_count_mean << endl;
        //cerr << estimated_count_covariance << endl;
        vector<double> frags_per_transcript;
        vector<double> frag_variances;
        
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            assert (estimated_count_mean.size() > i);
            frags_per_transcript.push_back(estimated_count_mean[i]);
            //frags_per_transcript.push_back(_abundances[i]->num_fragments());
            
            frag_variances.push_back(_abundances[i]->mass_variance());
        }
        
        simulate_count_covariance(frags_per_transcript, frag_variances, estimated_count_covariance, transcripts, _count_covariance, _assigned_count_samples, &gamma_samples);
        
        generate_fpkm_samples();
        
        calculate_FPKM_covariance();
        
        // Derive confidence intervals from the FPKM variance/covariance matrix
        calculate_conf_intervals();
        
        // Calculate the inter-group relative abundances and variances
        if (no_js_tests == false && _read_group_props.size() >= min_reps_for_js_test)
        {
            calculate_kappas();
        }
        
    }
    
//    //cerr << _count_covariance << endl;
//    for (size_t i = 0; i < _abundances.size(); ++i)
//    {
//        for (size_t j = 0; j < _abundances.size(); ++j)
//        {
//            if (i != j)
//            {
//                assert(!isinf(_fpkm_covariance(i,j)) && !isnan(_fpkm_covariance(i,j)));
//                if (_abundances[i]->transfrag()->contains(*_abundances[j]->transfrag()) &&
//                    Scaffold::compatible(*_abundances[i]->transfrag(),*_abundances[j]->transfrag()))
//                {
//                    _abundances[j]->status(NUMERIC_LOW_DATA);
//                }
//            }
//        }
//    }
    
    //assert (FPKM() == 0 || _assigned_count_samples.size() > 0);
    
    //fprintf(stderr, "Total calls to get_cond_prob = %d\n", total_cond_prob_calls);

}

void AbundanceGroup::calculate_abundance_for_replicate(const vector<MateHit>& alignments, bool perform_collapse)
{
    vector<shared_ptr<Abundance> > transcripts;
    
	get_transfrags(transcripts);
	vector<shared_ptr<Abundance> > mapped_transcripts; // This collects the transcripts that have alignments mapping to them
    
	vector<MateHit> nr_alignments;
    vector<double> joint_mle_gammas;
    
    vector<MateHit> non_equiv_alignments;
    vector<double> log_conv_factors;
    
    if (final_est_run || corr_multi || corr_bias)
    {
        compute_cond_probs_and_effective_lengths(alignments, transcripts, mapped_transcripts);
        collect_per_replicate_mass(alignments, transcripts);
        calculate_gammas(alignments, log_conv_factors, transcripts, mapped_transcripts);
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            joint_mle_gammas.push_back(_abundances[i]->gamma());
        }
    }
    
    calculate_iterated_exp_count_covariance(joint_mle_gammas, alignments, transcripts, _iterated_exp_count_covariance);
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        _abundances[i]->num_fragment_uncertainty_var(_iterated_exp_count_covariance(i,i));
    }
    
    if (corr_multi && !final_est_run)
    {
        update_multi_reads(alignments, mapped_transcripts);
    }
    
    // Calculate the initial estimates for the number of fragments originating
    // from each transcript, the FPKMs, and set the NB variances
    calculate_locus_scaled_mass_and_variance(transcripts);
}

void AbundanceGroup::aggregate_replicate_abundances(const map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >& ab_group_per_replicate)
{
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        CountPerReplicateTable cpr;
        FPKMPerReplicateTable fpr;
        StatusPerReplicateTable spr;
        
        double avg_fpkm = 0.0;
        double avg_num_frags = 0.0;
        double avg_gamma = 0.0;
        double avg_mass_variance = 0.0;
        double avg_effective_length = 0.0;
        
        map<AbundanceStatus, int> status_table;
        status_table[NUMERIC_OK] = 0;
        status_table[NUMERIC_LOW_DATA] = 0;
        status_table[NUMERIC_FAIL] = 0;
        status_table[NUMERIC_HI_DATA] = 0;
        
        for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
             itr != ab_group_per_replicate.end();
             ++itr)
        {
            status_table[itr->second->abundances()[i]->status()] += 1;
        }
        
        for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
             itr != ab_group_per_replicate.end();
             ++itr)
        {
            const vector<shared_ptr<Abundance> >& sc_ab = itr->second->abundances();
            assert(itr->second->abundances().size() == _abundances.size());
            cpr[itr->first] = itr->second->abundances()[i]->num_fragments();
            //fprintf(stderr, "FPKM = %lg\n", itr->second->abundances()[i]->FPKM());
            fpr[itr->first] = itr->second->abundances()[i]->FPKM();
            spr[itr->first] = itr->second->abundances()[i]->status();
            
            /*
            if (itr->second->abundances()[i]->status() == NUMERIC_OK)
            {
                avg_fpkm += itr->second->abundances()[i]->FPKM() / (double)status_table[NUMERIC_OK];
                avg_num_frags += itr->second->abundances()[i]->num_fragments() / (double)status_table[NUMERIC_OK];
                avg_gamma += itr->second->abundances()[i]->gamma() / (double)status_table[NUMERIC_OK];
                avg_mass_variance += itr->second->abundances()[i]->mass_variance() / (double)status_table[NUMERIC_OK];
                avg_effective_length += itr->second->abundances()[i]->effective_length() / (double)status_table[NUMERIC_OK];
            }
            */
            avg_fpkm += itr->second->abundances()[i]->FPKM() / (double)ab_group_per_replicate.size();
            avg_num_frags += itr->second->abundances()[i]->num_fragments() / (double)ab_group_per_replicate.size();
            avg_gamma += itr->second->abundances()[i]->gamma() / (double)ab_group_per_replicate.size();
            avg_mass_variance += itr->second->abundances()[i]->mass_variance() / (double)ab_group_per_replicate.size();
            avg_effective_length += itr->second->abundances()[i]->effective_length() / (double)ab_group_per_replicate.size();
            
        }
        
        _abundances[i]->FPKM(avg_fpkm);
        _abundances[i]->gamma(avg_gamma);
        _abundances[i]->num_fragments(avg_num_frags);
        _abundances[i]->mass_variance(avg_mass_variance);
        _abundances[i]->effective_length(avg_effective_length);
        
        
        // if there was at least one good replicate, set the status to OK.  The reduced power will be reflected
        // during testing
        if (status_table[NUMERIC_OK] >= 1)
        {
            _abundances[i]->status(NUMERIC_OK);
        }
        else
        {
            if (status_table[NUMERIC_LOW_DATA] >= status_table[NUMERIC_FAIL])
            {
                _abundances[i]->status(NUMERIC_LOW_DATA);
            }
            else if (status_table[NUMERIC_HI_DATA] >= status_table[NUMERIC_FAIL])
            {
                _abundances[i]->status(NUMERIC_HI_DATA);
            }
//            else if (status_table[NUMERIC_HI_DATA] >= status_table[NUMERIC_LOW_DATA]) // not sure this ever happens in practice
//            {
//                _abundances[i]->status(NUMERIC_FAIL);
//            }
            else
            {
                _abundances[i]->status(NUMERIC_FAIL);
            }
        }
        
        _abundances[i]->num_fragments_by_replicate(cpr);
        _abundances[i]->FPKM_by_replicate(fpr);
        _abundances[i]->status_by_replicate(spr);
    }
}

void AbundanceGroup::calculate_abundance(const vector<MateHit>& alignments, bool perform_collapse)
{
	vector<shared_ptr<Abundance> > transcripts;
    
	get_transfrags(transcripts);
    
    map<shared_ptr<ReadGroupProperties const >, vector<MateHit> > alignments_per_read_group;
    
    for(std::set<shared_ptr<ReadGroupProperties const > >::iterator itr = _read_group_props.begin();
        itr != _read_group_props.end();
        ++itr)
    {
        alignments_per_read_group[*itr] = vector<MateHit>();
        
        vector<MateHit>& rep_hits = alignments_per_read_group[*itr];
        
        for (size_t i = 0; i < alignments.size(); ++i)
        {
            if (alignments[i].read_group_props() == *itr)
            {
                rep_hits.push_back(alignments[i]);
            }
        }
    }
    
    collect_per_replicate_mass(alignments, transcripts);
    
    std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> > ab_group_per_replicate;
    
    calculate_per_replicate_abundances(transcripts,
                                       alignments_per_read_group,
                                       ab_group_per_replicate);
    
    std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AbundanceGroup> > const_ab_group_per_replicate;
    
    for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >::iterator itr = ab_group_per_replicate.begin();
         itr != ab_group_per_replicate.end(); ++itr)
    {
        const_ab_group_per_replicate[itr->first] = itr->second;
    }
    
    aggregate_replicate_abundances(const_ab_group_per_replicate);

    calculate_abundance_group_variance(transcripts, const_ab_group_per_replicate);
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
                               const vector<shared_ptr<Abundance> >& transcripts,
                               ublas::matrix<double>& count_covariance,
                               vector<Eigen::VectorXd>& assigned_count_samples,
                               vector<ublas::vector<double> >* gamma_samples = NULL)
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
        
        size_t num_frags = num_frag_assignments;
        if (gamma_samples && gamma_samples->empty() == false)
            num_frags = gamma_samples->size();
        
        assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * num_frags, Eigen::VectorXd::Zero(transcripts.size()));
        return true;
    }
    
    vector<ublas::vector<double> > assigned_gamma_samples;
    
    boost::mt19937 rng;
    
    if (gamma_samples == NULL || gamma_samples->empty())
    {
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
        
        multinormal_generator<double> generator(mle_frag_counts, mle_count_covar);
        
        for (size_t assign_idx = 0; assign_idx < num_frag_assignments; ++assign_idx)
        {
            boost::numeric::ublas::vector<double> random_count_assign;
            double total_sample_counts = 0;
            do {
                random_count_assign = generator.next_rand();
                //cerr << random_count_assign << endl;
                
                for (size_t r_idx = 0; r_idx < random_count_assign.size(); ++r_idx)
                {
                    if (random_count_assign(r_idx) < 0)
                        random_count_assign(r_idx) = 0;
                }
                
                total_sample_counts = accumulate(random_count_assign.begin(), random_count_assign.end(), 0.0);
                if (total_sample_counts > 0)
                    random_count_assign = total_frag_counts * (random_count_assign / total_sample_counts);
                else
                    random_count_assign = boost::numeric::ublas::zero_vector<double>(transcripts.size());
            } while(total_sample_counts <= 0);
            
            assigned_gamma_samples.push_back(random_count_assign);
        }
    }
    else
    {
        for (size_t assign_idx = 0; assign_idx < gamma_samples->size(); ++assign_idx)
        {
            boost::numeric::ublas::vector<double> random_count_assign = total_frag_counts * (*gamma_samples)[assign_idx];
            assigned_gamma_samples.push_back(random_count_assign);
        }
    }
    
    assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * assigned_gamma_samples.size(), Eigen::VectorXd::Zero(transcripts.size()));
    
    Eigen::VectorXd expected_generated_counts = Eigen::VectorXd::Zero(transcripts.size());
    
    
    boost::uniform_01<> uniform_dist;
    boost::mt19937 null_rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uniform_gen(null_rng, uniform_dist);
    
    
    for (size_t assign_idx = 0; assign_idx < assigned_gamma_samples.size(); ++assign_idx)
    {
        boost::numeric::ublas::vector<double>& random_count_assign = assigned_gamma_samples[assign_idx];
        
        //cerr << random_count_assign << endl;
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
    
    if (assigned_count_samples.empty() == false)
        count_covariance /= assigned_count_samples.size();
    
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        // Make sure we aren't below the fit for the single isoform case
        if (count_covariance(i,i) < ceil(frag_variances[i]))
        {
            //fprintf(stderr, "Counts for %d (var = %lg) are underdispersed, reverting to fitted variance model (%lg)\n", i, _count_covariance(i,i), ceil(_abundances[i]->mass_variance()));
            count_covariance(i,i) = ceil(frag_variances[i]);
            assert (!isinf(count_covariance(i,i)) && !isnan(count_covariance(i,i)));
        }
        
        // Check that we aren't below what the Poisson model says we ought to be at
        if (count_covariance(i,i) < ceil(num_fragments[i] + iterated_exp_count_covariance(i,i)))
        {
            //fprintf(stderr, "Counts for %d (var = %lg) are underdispersed, reverting to additive variance model (%lg)\n", i, _count_covariance(i,i),  ceil(_abundances[i]->num_fragments() + _iterated_exp_count_covariance(i,i)));
            count_covariance(i,i) = ceil(num_fragments[i] + iterated_exp_count_covariance(i,i));
            assert (!isinf(count_covariance(i,i)) && !isnan(count_covariance(i,i)));
        }
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            assert (!isinf(count_covariance(i,j)) && !isnan(count_covariance(i,j)));
        }
    }
    
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
            double fpkm_sample =  sample[j] / M;
            
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
    
    vector<Eigen::VectorXd> mem_fpkm_samples;
    for (size_t i = 0; i < fpkm_sample_vectors.size(); ++i)
    {
        Eigen::VectorXd sample(fpkm_sample_vectors[i].size());
        mem_fpkm_samples.push_back(sample);
        for (size_t j = 0; j < fpkm_sample_vectors[i].size(); ++j)
        {
            mem_fpkm_samples[i](j) = fpkm_sample_vectors[i][j];
        }
    }
    
    member_fpkm_samples(mem_fpkm_samples);
    
    fpkm_samples(group_sum_fpkm_samples);
    
    assert (group_sum_fpkm_samples.empty() == false);
}

void AbundanceGroup::calculate_FPKM_covariance()
{
	if (effective_length() == 0)
	{
		_fpkm_covariance = ublas::zero_matrix<double>(_abundances.size(), _abundances.size());
		return;
	}
    
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
            if (length_i > 0 && length_j > 0 & M > 0)
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
//                if (fpkm > 0 && fpkm_var == 0 )
//                {
//                    cerr << _count_covariance << endl;
//                }
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
    
//    cerr << "FPKM covariance: " << endl;
//    for (unsigned i = 0; i < _fpkm_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (_fpkm_covariance, i);
//        cerr << i << " : " << _abundances[i]->FPKM() << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }
    
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
        
        double avg_X_g = 0.0;
        double avg_mass_fraction = 0.0;
        
        int N = _abundances.size();
        
        vector<double> avg_mass_variances(_abundances.size(), 0.0);
        
        for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin();
             itr != _count_per_replicate.end();
             ++itr)
        {
            shared_ptr<ReadGroupProperties const> rg_props = itr->first;
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
            
            //assert (scaled_total_mass != 0.0);
            avg_mass_fraction += (scaled_mass / scaled_total_mass);
        }
        
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
                
		BOOST_FOREACH(shared_ptr<Abundance> pA, _abundances)
		{
			double FPKM_hi;
			double FPKM_lo;
			if (pA->effective_length() > 0 && avg_mass_fraction > 0)
			{
                double norm_frag_density = 1000000000;
                norm_frag_density /= pA->effective_length();
                
                norm_frag_density *= avg_mass_fraction;
                double fpkm_high = norm_frag_density;
                
                double total_mass = (num_fragments() / avg_mass_fraction);
                double fpkm_constant = 1000000000 / pA->effective_length() / total_mass;
                double var_fpkm = mass_variance() * (fpkm_constant * fpkm_constant);
                
				FPKM_hi = fpkm_high + 2 * sqrt(var_fpkm);
				FPKM_lo = 0.0;
				ConfidenceInterval conf(FPKM_lo, FPKM_hi);
				//assert (FPKM_lo <= pA->FPKM() && pA->FPKM() <= FPKM_hi);
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


/*
void AbundanceGroup::calculate_conf_intervals()
{
    for (size_t j = 0; j < _abundances.size(); ++j)
    {
        double FPKM_hi = 0.0;
        double FPKM_lo = numeric_limits<double>::max();
        const vector<double> ab_j_samples = _abundances[j]->fpkm_samples();
        vector<pair<double, double> > fpkm_samples;
        double target_fpkm = _abundances[j]->FPKM();
        for (size_t i = 0; i < ab_j_samples.size(); ++i)
            fpkm_samples.push_back(make_pair(abs(ab_j_samples[i] - _abundances[j]->FPKM()), ab_j_samples[i]));
        
        sort(fpkm_samples.begin(), fpkm_samples.end());
        
        for (size_t i = 0; i < 0.95*fpkm_samples.size(); ++i)
        {
            if (FPKM_lo > fpkm_samples[i].second)
                FPKM_lo = fpkm_samples[i].second;
            if (FPKM_hi < fpkm_samples[i].second)
                FPKM_hi = fpkm_samples[i].second;
        }
        
        if (fpkm_samples.size() > 0 && FPKM_lo > target_fpkm)
        {
            fprintf(stderr, "Warning: transcript confidence interval lower bound is > FPKM\n");
        }
        
        if (fpkm_samples.size() > 0 && FPKM_hi < target_fpkm)
        {
            fprintf(stderr, "Warning: transcript confidence interval upper bound is < FPKM\n");
        }
        
        ConfidenceInterval conf(FPKM_lo, FPKM_hi);
        _abundances[j]->FPKM_conf(conf);
    }
    
    double FPKM_hi = 0.0;
    double FPKM_lo = numeric_limits<double>::max();
    const vector<double>& ab_j_samples = _fpkm_samples;
    vector<pair<double, double> > fpkm_samples;
    for (size_t i = 0; i < ab_j_samples.size(); ++i)
        fpkm_samples.push_back(make_pair(abs(ab_j_samples[i] - FPKM()), ab_j_samples[i]));
    
    sort(fpkm_samples.begin(), fpkm_samples.end());
    
    for (size_t i = 0; i < 0.95*fpkm_samples.size(); ++i)
    {
        if (FPKM_lo > fpkm_samples[i].second)
            FPKM_lo = fpkm_samples[i].second;
        if (FPKM_hi < fpkm_samples[i].second)
            FPKM_hi = fpkm_samples[i].second;
    }
    
    double target_fpkm = FPKM();
    
    if (fpkm_samples.size() > 0 && FPKM_lo > target_fpkm)
    {
        fprintf(stderr, "Warning: group confidence interval lower bound is > FPKM\n");
    }
    
    if (fpkm_samples.size() > 0 && FPKM_hi < target_fpkm)
    {
        fprintf(stderr, "Warning: group confidence interval upper bound is < FPKM\n");
    }
    
    ConfidenceInterval conf(FPKM_lo, FPKM_hi);
    FPKM_conf(conf);
    
    return;
}

*/
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
    //cerr << "expected counts" << endl;
    //cerr << expected_counts << endl;
    
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
                                                                   const map<shared_ptr<ReadGroupProperties const >, vector<MateHit> >& alignments_per_read_group,
                                                                   std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AbundanceGroup> >& ab_group_per_replicate,
                                                                   bool perform_collapse)
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
        
        map<shared_ptr<ReadGroupProperties const >, vector<MateHit> >::const_iterator al_itr =
            alignments_per_read_group.find(*itr);
        
        assert(al_itr != alignments_per_read_group.end());
        const vector<MateHit>& rep_hits = alignments_per_read_group.find(*itr)->second;
       
        vector<MateHit> nr_alignments;
        vector<MateHit> non_equiv_alignments;
        vector<double> log_conv_factors;
        vector<shared_ptr<Abundance> > mapped_transcripts;
        
        if (perform_collapse)
        {
            collapse_hits(rep_hits, nr_alignments);
            collapse_equivalent_hits_helper(nr_alignments, transcripts, non_equiv_alignments, log_conv_factors);
            assert (non_equiv_alignments.size() == log_conv_factors.size());
            log_conv_factors = vector<double>(nr_alignments.size(), 0);
            nr_alignments.clear();
        }
        else
        {
            nr_alignments = rep_hits;
            non_equiv_alignments = nr_alignments;
            log_conv_factors = vector<double>(nr_alignments.size(), 0);
        }
        
        //rep_hit_counts.push_back(count_per_replicate.find(*itr)->second);
        
        ab_group->calculate_abundance_for_replicate(non_equiv_alignments, false);
        
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
        size_type i_norm_inf = i + index_norm_inf (project (mci, boost::numeric::ublas::range (i, size1)));
        if (m (i_norm_inf, i) != value_type/*zero*/()) {
            if (i_norm_inf != i) {
                pm (i) = i_norm_inf;
                row (m, i_norm_inf).swap (mri);
            } else {
                //BOOST_UBLAS_CHECK (pm (i) == i_norm_inf, external_logic ());
            }
            project (mci, boost::numeric::ublas::range (i + 1, size1)) *= value_type (1) / m (i, i);
        } else if (singular == 0) {
            singular = i + 1;
        }
        project (m, boost::numeric::ublas::range (i + 1, size1), boost::numeric::ublas::range (i + 1, size2)).minus_assign (outer_prod (project (mci, boost::numeric::ublas::range (i + 1, size1)),
                                                                              project (mri, boost::numeric::ublas::range (i + 1, size2))));
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

void tss_analysis(const string& locus_tag, SampleAbundances& sample)
{
    // Cluster transcripts by start site (TSS)
    vector<AbundanceGroup> transcripts_by_tss;
    
    ublas::matrix<double> tss_gamma_cov;
    ublas::matrix<double> tss_count_cov;
    ublas::matrix<double> tss_iterated_exp_count_cov;
    ublas::matrix<double> tss_fpkm_cov;
    vector<Eigen::VectorXd> tss_assigned_counts;
    
    vector<bool> mask(sample.transcripts.abundances().size(), true);
    for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
    {
        if (*(sample.transcripts.abundances()[i]->tss_id().begin()) == "")
        {
            mask[i] = false;
        }
    }
    
    AbundanceGroup trans_with_tss;
    sample.transcripts.filter_group(mask, trans_with_tss);
    
    cluster_transcripts<ConnectByAnnotatedTssId>(trans_with_tss,
                                                 transcripts_by_tss,
                                                 &tss_gamma_cov,
                                                 &tss_iterated_exp_count_cov,
                                                 &tss_count_cov,
                                                 &tss_fpkm_cov);
    
    
    BOOST_FOREACH(AbundanceGroup& ab_group, transcripts_by_tss)
    {
        ab_group.locus_tag(locus_tag);
        set<string> tss_ids = ab_group.tss_id();
        assert (tss_ids.size() == 1);
        string desc = *(tss_ids.begin());
        assert (desc != "");
        ab_group.description(*(tss_ids.begin()));
    }
    
    sample.primary_transcripts = transcripts_by_tss;
    
    // Group TSS clusters by gene
    vector<shared_ptr<Abundance> > primary_transcript_abundances;
    set<shared_ptr<ReadGroupProperties const> > rg_props;
    BOOST_FOREACH (AbundanceGroup& ab_group, sample.primary_transcripts)
    {
        primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
        rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
    }
    
    AbundanceGroup primary_transcripts(primary_transcript_abundances,
                                       tss_gamma_cov,
                                       tss_iterated_exp_count_cov,
                                       tss_count_cov,
                                       tss_fpkm_cov,
                                       rg_props);
    
    vector<AbundanceGroup> primary_transcripts_by_gene;
    
    cluster_transcripts<ConnectByAnnotatedGeneId>(primary_transcripts,
                                                  primary_transcripts_by_gene);
    
    BOOST_FOREACH(AbundanceGroup& ab_group, primary_transcripts_by_gene)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        if (gene_ids.size() > 1)
        {
            BOOST_FOREACH (string st, gene_ids)
            {
                fprintf(stderr, "%s\n", st.c_str());
            }
            ab_group.gene_id();
        }
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
    
    sample.gene_primary_transcripts = primary_transcripts_by_gene;
}

void cds_analyis(const string& locus_tag, SampleAbundances& sample)
{
    // Cluster transcripts by CDS
    vector<AbundanceGroup> transcripts_by_cds;
    ublas::matrix<double> cds_gamma_cov;
    ublas::matrix<double> cds_count_cov;
    ublas::matrix<double> cds_iterated_exp_count_cov;
    ublas::matrix<double> cds_fpkm_cov;
    
    vector<bool> mask(sample.transcripts.abundances().size(), true);
    for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
    {
        if (*(sample.transcripts.abundances()[i]->protein_id().begin()) == "")
        {
            mask[i] = false;
        }
    }
    
    AbundanceGroup trans_with_p_id;
    sample.transcripts.filter_group(mask, trans_with_p_id);
    
    cluster_transcripts<ConnectByAnnotatedProteinId>(trans_with_p_id,
                                                     transcripts_by_cds,
                                                     &cds_gamma_cov,
                                                     &cds_iterated_exp_count_cov,
                                                     &cds_count_cov,
                                                     &cds_fpkm_cov);
    
    BOOST_FOREACH(AbundanceGroup& ab_group, transcripts_by_cds)
    {
        ab_group.locus_tag(locus_tag);
        set<string> protein_ids = ab_group.protein_id();
        assert (protein_ids.size() == 1);
        string desc = *(protein_ids.begin());
        //if (desc != "")
        //{
        assert (desc != "");
        ab_group.description(*(protein_ids.begin()));
        //}
    }
    
    sample.cds = transcripts_by_cds;
    
    // Group the CDS clusters by gene
    vector<shared_ptr<Abundance> > cds_abundances;
    
    set<shared_ptr<ReadGroupProperties const> > rg_props;
    BOOST_FOREACH (AbundanceGroup& ab_group, sample.cds)
    {
        //if (ab_group.description() != "")
        {
            cds_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
            rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
        }
    }
    AbundanceGroup cds(cds_abundances,
                       cds_gamma_cov,
                       cds_iterated_exp_count_cov,
                       cds_count_cov,
                       cds_fpkm_cov,
                       rg_props);
    
    vector<AbundanceGroup> cds_by_gene;
    
    cluster_transcripts<ConnectByAnnotatedGeneId>(cds,
                                                  cds_by_gene);
    
    BOOST_FOREACH(AbundanceGroup& ab_group, cds_by_gene)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
    
    sample.gene_cds = cds_by_gene;

}

void sample_abundance_worker(const string& locus_tag,
                             const set<shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAbundances& sample,
                             shared_ptr<HitBundle> sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis)
{
    vector<shared_ptr<Abundance> > abundances;
    
    BOOST_FOREACH(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        TranscriptAbundance* pT = new TranscriptAbundance;
        pT->transfrag(s);
        shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }
    
    sample.transcripts = AbundanceGroup(abundances);
    
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
        get_alignments_from_scaffolds(sample.transcripts.abundances(),
                                      hits_in_cluster);
        
        // Compute the individual transcript FPKMs via each sample's
        // AbundanceGroup for this locus.
        
        sample.transcripts.calculate_abundance(hits_in_cluster);
    }
    else
    {
        BOOST_FOREACH(shared_ptr<Abundance>  ab, abundances)
        {
            ab->status(NUMERIC_HI_DATA);
            
            CountPerReplicateTable cpr;
            FPKMPerReplicateTable fpr;
            StatusPerReplicateTable spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                cpr[*itr] = 0;
                fpr[*itr] = 0;
                spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_fragments_by_replicate(cpr);
            ab->FPKM_by_replicate(fpr);
            ab->status_by_replicate(spr);
        }
    }
    
    // Cluster transcripts by gene_id
    vector<AbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	BOOST_FOREACH(AbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
    
    if (perform_cds_analysis)
    {
        cds_analyis(locus_tag, sample);
    }
    
    if (perform_tss_analysis)
    {
        tss_analysis(locus_tag, sample);
    }
}

// This function applies library size factors to pre-computed expression entries
void AbundanceGroup::apply_normalization_to_abundances(const map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AbundanceGroup> >& unnormalized_ab_group_per_replicate,
                                                       map<shared_ptr<const ReadGroupProperties>, shared_ptr<AbundanceGroup> >& normalized_ab_group_per_replicate)
{
    for (map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AbundanceGroup> >::const_iterator itr = unnormalized_ab_group_per_replicate.begin();
         itr != unnormalized_ab_group_per_replicate.end(); ++itr)
    {
        shared_ptr<AbundanceGroup> norm_ab = shared_ptr<AbundanceGroup>(new AbundanceGroup(*itr->second));
        shared_ptr<const ReadGroupProperties> rg_props = itr->first;
        shared_ptr<const MassDispersionModel> disp_model = rg_props->mass_dispersion_model();
        
        shared_ptr<const ReadGroupProperties> old_rg_props = *(itr->second->rg_props().begin());
        
        double fpkm_correction_factor = old_rg_props->normalized_map_mass() / rg_props->normalized_map_mass();
        double internal_scale_factor = rg_props->internal_scale_factor();
        
        double total_mass = 0.0;
        
        for (size_t i = 0; i < norm_ab->_abundances.size(); ++i)
        {
            norm_ab->_abundances[i]->num_fragments(itr->second->_abundances[i]->num_fragments() / internal_scale_factor);
            
            total_mass += norm_ab->_abundances[i]->num_fragments();
            
            norm_ab->_abundances[i]->FPKM(fpkm_correction_factor * itr->second->_abundances[i]->FPKM() / internal_scale_factor);
            norm_ab->_iterated_exp_count_covariance = norm_ab->iterated_count_cov() / (internal_scale_factor*internal_scale_factor);
            norm_ab->_fpkm_covariance = norm_ab->_fpkm_covariance * (fpkm_correction_factor * fpkm_correction_factor)/ (internal_scale_factor*internal_scale_factor);
            norm_ab->_count_covariance = norm_ab->_count_covariance/ (internal_scale_factor*internal_scale_factor);
        }
        
        double locus_mass_variance = disp_model->scale_mass_variance(total_mass);
        
        for (size_t i = 0; i < norm_ab->_abundances.size(); ++i)
        {
            norm_ab->_abundances[i]->mass_variance(locus_mass_variance * norm_ab->_abundances[i]->gamma());
        }
        
        normalized_ab_group_per_replicate[itr->first] = norm_ab;
    }
}

void merge_precomputed_expression_worker(const string& locus_tag,
                                         const vector<shared_ptr<PrecomputedExpressionBundleFactory> >& expression_factories,
                                         SampleAbundances& sample,
                                         shared_ptr<HitBundle> sample_bundle,
                                         bool perform_cds_analysis,
                                         bool perform_tss_analysis)
{
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AbundanceGroup> > unnormalized_ab_group_per_replicate;
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<AbundanceGroup> > normalized_ab_group_per_replicate;
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AbundanceGroup> > const_ab_group_per_replicate;
    
    if (locus_tag == "chr6:31126318-31148508")
    {
        int a = 4;
    }
    
    set<shared_ptr<const ReadGroupProperties> > rg_props;
    for (size_t i = 0; i < expression_factories.size(); ++i)
    {
        shared_ptr<PrecomputedExpressionBundleFactory> pBundleFac = expression_factories[i];
		shared_ptr<const PrecomputedExpressionHitFactory> pHitFac = dynamic_pointer_cast<const PrecomputedExpressionHitFactory> (pBundleFac->hit_factory());
        assert (pHitFac);
        
        shared_ptr<const ReadGroupProperties> rg_prop = pBundleFac->read_group_properties();
        rg_props.insert(rg_prop);
		shared_ptr<const AbundanceGroup> ab = pBundleFac->get_abundance_for_locus(sample_bundle->id());
		pBundleFac->clear_abundance_for_locus(sample_bundle->id());
        if (!ab)
        {
            fprintf(stderr, "Error: no bundle with id %d in precomputed expression file\n", sample_bundle->id());
        }
        else if(ab->abundances().size() != sample_bundle->ref_scaffolds().size())
        {
            fprintf(stderr, "Error: bad bundle merge %s != %s\n", ab->description().c_str(), locus_tag.c_str());
        }
        unnormalized_ab_group_per_replicate[rg_prop] = ab;
    }
    
    AbundanceGroup::apply_normalization_to_abundances(unnormalized_ab_group_per_replicate, normalized_ab_group_per_replicate);
    
    for (map<shared_ptr<const ReadGroupProperties>, shared_ptr<AbundanceGroup> >::const_iterator itr = normalized_ab_group_per_replicate.begin();
         itr != normalized_ab_group_per_replicate.end(); ++itr)
    {
        const_ab_group_per_replicate[itr->first] = itr->second;
    }
    
    vector<shared_ptr<Abundance> > abundances;
    
    BOOST_FOREACH(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        TranscriptAbundance* pT = new TranscriptAbundance;
        pT->transfrag(s);
        shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }
    
    sample.transcripts = AbundanceGroup(abundances);
    
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
        sample.transcripts.collect_per_replicate_mass(const_ab_group_per_replicate);
        sample.transcripts.aggregate_replicate_abundances(const_ab_group_per_replicate);
        sample.transcripts.calculate_abundance_group_variance(abundances, const_ab_group_per_replicate);
    }
    else
    {
        BOOST_FOREACH(shared_ptr<Abundance>  ab, abundances)
        {
            ab->status(NUMERIC_HI_DATA);
            
            CountPerReplicateTable cpr;
            FPKMPerReplicateTable fpr;
            StatusPerReplicateTable spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                cpr[*itr] = 0;
                fpr[*itr] = 0;
                spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_fragments_by_replicate(cpr);
            ab->FPKM_by_replicate(fpr);
            ab->status_by_replicate(spr);
        }
    }
    
    // Cluster transcripts by gene_id
    vector<AbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	BOOST_FOREACH(AbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
    
    if (perform_cds_analysis)
    {
        cds_analyis(locus_tag, sample);
    }
    
    if (perform_tss_analysis)
    {
        tss_analysis(locus_tag, sample);
    }
}

//allele implementation
void compute_compatibilities(const vector<shared_ptr<Abundance> >& transcripts,
							 const vector<MateHit>& alignments,
							 vector<vector<char> >& paternalCompatibilities,
							 vector<vector<char> >& maternalCompatibilities)
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
				AlleleInfo allele = alignments[i].allele();
				if(allele == ALLELE_PATERNAL)
				{
					paternalCompatibilities[j][i] = 1;
				}
				else if(allele == ALLELE_MATERNAL)
				{
					maternalCompatibilities[j][i] = 1;
				}	
				else
				{
					paternalCompatibilities[j][i] = 1;
					maternalCompatibilities[j][i] = 1;
				}
			}
		}	
	}
}

AlleleAbundanceGroup::AlleleAbundanceGroup(const vector<shared_ptr<Abundance> >& abundances,
                               const ublas::matrix<double>& parental_gamma_covariance,
                               const ublas::matrix<double>& iterated_exp_parental_count_covariance,
                               const ublas::matrix<double>& parental_count_covariance,
                               const ublas::matrix<double>& parental_fpkm_covariance,
                               const set<shared_ptr<ReadGroupProperties const> >& rg_props) :
    _abundances(abundances), 
    _iterated_exp_parental_count_covariance(iterated_exp_parental_count_covariance),
    _parental_count_covariance(parental_count_covariance),
    _parental_fpkm_covariance(parental_fpkm_covariance),
    _parental_gamma_covariance(parental_gamma_covariance),
    _salient_frags(0.0),
    _total_frags(0.0),
    _read_group_props(rg_props)
{
	
    // Calling calculate_FPKM_covariance() also estimates cross-replicate
    // count variances
    // calculate_FPKM_covariance();
    double paternal_fpkm_var = 0.0;
    for (size_t i = 0; i < _parental_fpkm_covariance.size1()/2; ++i)
    {
        for (size_t j = 0; j < _parental_fpkm_covariance.size2(); ++j)
        {
            assert (!isnan(_parental_fpkm_covariance(i,j)) && !isinf(_parental_fpkm_covariance(i,j)));
            paternal_fpkm_var += _parental_fpkm_covariance(i,j);
        }
    }
    
    _paternal_FPKM_variance = paternal_fpkm_var;
	
	double maternal_fpkm_var = 0.0;
    for (size_t i = _parental_fpkm_covariance.size1()/2; i < _parental_fpkm_covariance.size1(); ++i)
    {
        for (size_t j = 0; j < _parental_fpkm_covariance.size2(); ++j)
        {
            assert (!isnan(_parental_fpkm_covariance(i,j)) && !isinf(_parental_fpkm_covariance(i,j)));
            maternal_fpkm_var += _parental_fpkm_covariance(i,j);
        }
    }
    
    _maternal_FPKM_variance = maternal_fpkm_var;
    
    if (paternal_FPKM()+maternal_FPKM() > 0 && final_est_run && library_type != "transfrags")
    {
        
        ublas::matrix<double> test = _parental_fpkm_covariance;
				
        double ret = cholesky_factorize(test);
        if (ret != 0 || ((_paternal_FPKM_variance < 0 || _maternal_FPKM_variance < 0) && (paternal_status() == NUMERIC_OK && maternal_status() == NUMERIC_OK)))
        {
            fprintf(stderr, "Warning: total FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            for (size_t j = 0; j < _abundances.size(); ++j)
            {
                _abundances[j]->paternal_status(NUMERIC_FAIL);
				_abundances[j]->maternal_status(NUMERIC_FAIL);
            }
        }
        
		if(!(paternal_FPKM()+maternal_FPKM() == 0 || (paternal_fpkm_var > 0 || maternal_fpkm_var > 0) || (paternal_status() != NUMERIC_OK || maternal_status() != NUMERIC_OK)))
       {
           //cerr << _parental_count_covariance << endl;
           //cerr << _parental_fpkm_covariance << endl;
       }
        
        assert (paternal_FPKM()+maternal_FPKM() == 0 || (paternal_fpkm_var > 0 || maternal_fpkm_var > 0) ||(paternal_status() != NUMERIC_OK || maternal_status() != NUMERIC_OK));
    }
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        if (_paternal_fpkm_samples.empty())
        {
			_paternal_fpkm_samples = vector<double>(_abundances[i]->paternal_fpkm_samples().size(), 0);
            _paternal_member_fpkm_samples = vector<Eigen::VectorXd>(_abundances[i]->paternal_fpkm_samples().size(), Eigen::VectorXd::Zero(_abundances.size()));
        }
		if (_maternal_fpkm_samples.empty())
        {
			_maternal_fpkm_samples = vector<double>(_abundances[i]->maternal_fpkm_samples().size(), 0);
            _maternal_member_fpkm_samples = vector<Eigen::VectorXd>(_abundances[i]->maternal_fpkm_samples().size(), Eigen::VectorXd::Zero(_abundances.size()));
        }
		assert(_abundances[i]->paternal_fpkm_samples().size() == _abundances[i]->maternal_fpkm_samples().size());
		for (size_t j = 0; j < _abundances[i]->paternal_fpkm_samples().size(); ++j)
        {
            _paternal_fpkm_samples[j] += _abundances[i]->paternal_fpkm_samples()[j];
            _paternal_member_fpkm_samples[j](i) = _abundances[i]->paternal_fpkm_samples()[j];
			_maternal_fpkm_samples[j] += _abundances[i]->maternal_fpkm_samples()[j];
            _maternal_member_fpkm_samples[j](i) = _abundances[i]->maternal_fpkm_samples()[j];
        }
    }
		
    calculate_conf_intervals();
    
    if (no_js_tests == false && _read_group_props.size() >= min_reps_for_js_test)
    {
        calculate_kappas();
    }
}

AbundanceStatus AlleleAbundanceGroup::paternal_status() const
{
    bool has_lowdata_member = false;
    bool has_ok_member = false;
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->paternal_status() == NUMERIC_FAIL)
		{
			return NUMERIC_FAIL;
		}
        else if (ab->paternal_status() == NUMERIC_LOW_DATA)
		{
			has_lowdata_member = true;
            //return NUMERIC_LOW_DATA;
		}
        else if (ab->paternal_status() == NUMERIC_HI_DATA)
		{
			return NUMERIC_HI_DATA;
		}
        else if (ab->paternal_status() == NUMERIC_OK)
		{
			has_ok_member = true;
		}
	}
    
    if (has_ok_member == false)
        return NUMERIC_LOW_DATA;
    
	return NUMERIC_OK;
}

AbundanceStatus AlleleAbundanceGroup::maternal_status() const
{
    bool has_lowdata_member = false;
    bool has_ok_member = false;
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->maternal_status() == NUMERIC_FAIL)
		{
			return NUMERIC_FAIL;
		}
        else if (ab->maternal_status() == NUMERIC_LOW_DATA)
		{
			has_lowdata_member = true;
            //return NUMERIC_LOW_DATA;
		}
        else if (ab->maternal_status() == NUMERIC_HI_DATA)
		{
			return NUMERIC_HI_DATA;
		}
        else if (ab->maternal_status() == NUMERIC_OK)
		{
			has_ok_member = true;
		}
	}
    
    if (has_ok_member == false)
        return NUMERIC_LOW_DATA;
    
	return NUMERIC_OK;
}

void AlleleTranscriptAbundance::paternal_FPKM_variance(double v)
{ 
    assert (v >= 0); 
    assert(!isnan(v));
    _paternal_FPKM_variance = v; 
}

void AlleleTranscriptAbundance::maternal_FPKM_variance(double v)
{ 
    assert (v >= 0); 
    assert(!isnan(v));
    _maternal_FPKM_variance = v; 
}

bool AlleleAbundanceGroup::has_paternal_member_with_status(AbundanceStatus member_status) const
{
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->paternal_status() == member_status)
		{
			return true;
		}
	}
    return false;
}

bool AlleleAbundanceGroup::has_maternal_member_with_status(AbundanceStatus member_status) const
{
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (ab->maternal_status() == member_status)
		{
			return true;
		}
	}
    return false;
}

double AlleleAbundanceGroup::num_paternal_fragments() const
{
	double num_f = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		num_f += ab->num_paternal_fragments();
	}
    assert (!isnan(num_f));
	return num_f;
}

double AlleleAbundanceGroup::num_maternal_fragments() const
{
	double num_f = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		num_f += ab->num_maternal_fragments();
	}
    assert (!isnan(num_f));
	return num_f;
}

CountPerReplicateTable AlleleAbundanceGroup::num_paternal_fragments_by_replicate() const
{
	CountPerReplicateTable cpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (cpr.empty())
        {
            cpr = ab->num_paternal_fragments_by_replicate();
        }
        else
        {
            CountPerReplicateTable ab_cpr = ab->num_paternal_fragments_by_replicate();
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

CountPerReplicateTable AlleleAbundanceGroup::num_maternal_fragments_by_replicate() const
{
	CountPerReplicateTable cpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (cpr.empty())
        {
            cpr = ab->num_maternal_fragments_by_replicate();
        }
        else
        {
            CountPerReplicateTable ab_cpr = ab->num_maternal_fragments_by_replicate();
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

FPKMPerReplicateTable AlleleAbundanceGroup::paternal_FPKM_by_replicate() const
{
	FPKMPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
        FPKMPerReplicateTable ab_fpr = ab->paternal_FPKM_by_replicate();
        
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

FPKMPerReplicateTable AlleleAbundanceGroup::maternal_FPKM_by_replicate() const
{
	FPKMPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
        FPKMPerReplicateTable ab_fpr = ab->maternal_FPKM_by_replicate();
        
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

StatusPerReplicateTable AlleleAbundanceGroup::paternal_status_by_replicate() const
{
	StatusPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (fpr.empty())
        {
            fpr = ab->paternal_status_by_replicate();
        }
        else
        {
            StatusPerReplicateTable ab_fpr = ab->paternal_status_by_replicate();
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

StatusPerReplicateTable AlleleAbundanceGroup::maternal_status_by_replicate() const
{
	StatusPerReplicateTable fpr;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		if (fpr.empty())
        {
            fpr = ab->maternal_status_by_replicate();
        }
        else
        {
            StatusPerReplicateTable ab_fpr = ab->maternal_status_by_replicate();
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

double AlleleAbundanceGroup::paternal_mass_variance() const
{
    double mass_var = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		mass_var += ab->paternal_mass_variance();
	}
	return mass_var;
}

double AlleleAbundanceGroup::maternal_mass_variance() const
{
    double mass_var = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		mass_var += ab->maternal_mass_variance();
	}
	return mass_var;
}


// This tracks the final modeled variance in the assigned counts.
double AlleleAbundanceGroup::num_paternal_fragment_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = 0; i < _parental_count_covariance.size1()/2; ++i)
    {
        for (size_t j = 0; j < _parental_count_covariance.size2(); ++j)
        {
            frag_var += _parental_count_covariance(i,j);
        }
    }
    return frag_var;
}

double AlleleAbundanceGroup::num_maternal_fragment_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = _parental_count_covariance.size1()/2; i < _parental_count_covariance.size1(); ++i)
    {
        for (size_t j = 0; j < _parental_count_covariance.size2(); ++j)
        {
            frag_var += _parental_count_covariance(i,j);
        }
    }
    return frag_var;
}

// This tracks the final modeled variance in the assigned counts.
double AlleleAbundanceGroup::num_paternal_fragment_uncertainty_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = 0; i < _iterated_exp_parental_count_covariance.size1()/2; ++i)
    {
        for (size_t j = 0; j < _iterated_exp_parental_count_covariance.size2(); ++j)
        {
            frag_var += _iterated_exp_parental_count_covariance(i,j);
        }
    }
    return frag_var;
}

double AlleleAbundanceGroup::num_maternal_fragment_uncertainty_var() const			
{ 
    double frag_var = 0.0;
    for (size_t i = _iterated_exp_parental_count_covariance.size1()/2; i < _iterated_exp_parental_count_covariance.size1(); ++i)
    {
        for (size_t j = 0; j < _iterated_exp_parental_count_covariance.size2(); ++j)
        {
            frag_var += _iterated_exp_parental_count_covariance(i,j);
        }
    }
	return frag_var;
}

double AlleleAbundanceGroup::paternal_FPKM() const
{
	double fpkm = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		fpkm += ab->paternal_FPKM();
	}
	
	return fpkm;
}

double AlleleAbundanceGroup::maternal_FPKM() const
{
	double fpkm = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		fpkm += ab->maternal_FPKM();
	}
	
	return fpkm;
}

double AlleleAbundanceGroup::paternal_gamma() const
{
	double gamma = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		gamma += ab->paternal_gamma();
	}
	
	return gamma;
}

double AlleleAbundanceGroup::maternal_gamma() const
{
	double gamma = 0;
	
	BOOST_FOREACH(shared_ptr<Abundance> ab, _abundances)
	{
		gamma += ab->maternal_gamma();
	}
	
	return gamma;
}

void AlleleTranscriptAbundance::clear_non_serialized_data()
{
    _paternal_fpkm_samples.clear();
	_maternal_fpkm_samples.clear();
    std::vector<double>().swap(_paternal_fpkm_samples);
	std::vector<double>().swap(_maternal_fpkm_samples);
    
    if (_paternal_cond_probs)
    {
        _paternal_cond_probs->clear();
        std::vector<double>().swap(*_paternal_cond_probs);
    }
    if (_maternal_cond_probs)
    {
        _maternal_cond_probs->clear();
        std::vector<double>().swap(*_maternal_cond_probs);
    }
    
    if (_transfrag)
    {
        _transfrag->clear_hits();
        _transfrag = shared_ptr<Scaffold>();
    }
}

void AlleleAbundanceGroup::clear_non_serialized_data()
{
    
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        _abundances[i]->clear_non_serialized_data();
    }
    
    _paternal_fpkm_samples.clear();
	_maternal_fpkm_samples.clear();
    std::vector<double>().swap(_paternal_fpkm_samples);
	std::vector<double>().swap(_maternal_fpkm_samples);
    _paternal_member_fpkm_samples.clear();
	_maternal_member_fpkm_samples.clear();
    std::vector<Eigen::VectorXd>().swap(_paternal_member_fpkm_samples);
	std::vector<Eigen::VectorXd>().swap(_maternal_member_fpkm_samples);
    _paternal_assigned_count_samples.clear();
	_maternal_assigned_count_samples.clear();
    std::vector<Eigen::VectorXd>().swap(_paternal_assigned_count_samples);
	std::vector<Eigen::VectorXd>().swap(_maternal_assigned_count_samples);
}

void AlleleAbundanceGroup::filter_group(const vector<bool>& to_keep, 
								  AlleleAbundanceGroup& filtered_group) const
{
	//filtered_group = AlleleAbundanceGroup();
	
	assert (to_keep.size() == _abundances.size());
	
	size_t num_kept = 0;
	BOOST_FOREACH(bool keeper, to_keep)
	{
		num_kept += keeper;
	}
	
	ublas::matrix<double> new_parental_cov = ublas::zero_matrix<double>(2*num_kept,2*num_kept);
    ublas::matrix<double> new_iterated_em_parental_count_cov = ublas::zero_matrix<double>(2*num_kept,2*num_kept);
    ublas::matrix<double> new_parental_count_cov = ublas::zero_matrix<double>(2*num_kept,2*num_kept);
    ublas::matrix<double> new_parental_fpkm_cov = ublas::zero_matrix<double>(2*num_kept,2*num_kept);
    
	vector<shared_ptr<Abundance> > new_ab;
    
	//vector<vector<double> > new_paternal_fpkm_samples(_paternal_fpkm_samples.size(), vector<double>(num_kept, 0));
	//vector<vector<double> > new_maternal_fpkm_samples(_maternal_fpkm_samples.size(), vector<double>(num_kept, 0));
	
    // rebuild covariance matrix and abundance vector after filtration
	
	size_t next_cov_row = 0;
	for (size_t i = 0; i < to_keep.size(); ++i)
	{
		if (to_keep[i])
		{
			new_ab.push_back(_abundances[i]);
			size_t next_cov_col = 0;
			for (size_t j = 0; j < to_keep.size(); ++j)
			{
				if (to_keep[j])
				{
					new_parental_cov(next_cov_row,next_cov_col) = _parental_gamma_covariance(i, j);
					new_parental_cov(next_cov_row,next_cov_col+num_kept) = _parental_gamma_covariance(i, j+to_keep.size());
					new_parental_cov(next_cov_row+num_kept,next_cov_col+num_kept) = _parental_gamma_covariance(i+to_keep.size(), j+to_keep.size());
					new_parental_cov(next_cov_row+num_kept,next_cov_col) = _parental_gamma_covariance(i+to_keep.size(), j);
					
					new_iterated_em_parental_count_cov(next_cov_row,next_cov_col) = _iterated_exp_parental_count_covariance(i, j);
					new_iterated_em_parental_count_cov(next_cov_row,next_cov_col+num_kept) = _iterated_exp_parental_count_covariance(i, j+to_keep.size());
					new_iterated_em_parental_count_cov(next_cov_row+num_kept,next_cov_col+num_kept) = _iterated_exp_parental_count_covariance(i+to_keep.size(), j+to_keep.size());
					new_iterated_em_parental_count_cov(next_cov_row+num_kept,next_cov_col) = _iterated_exp_parental_count_covariance(i+to_keep.size(), j);
					
					new_parental_count_cov(next_cov_row,next_cov_col) = _parental_count_covariance(i, j);
					new_parental_count_cov(next_cov_row,next_cov_col+num_kept) = _parental_count_covariance(i, j+to_keep.size());
					new_parental_count_cov(next_cov_row+num_kept,next_cov_col+num_kept) = _parental_count_covariance(i+to_keep.size(), j+to_keep.size());
					new_parental_count_cov(next_cov_row+num_kept,next_cov_col) = _parental_count_covariance(i+to_keep.size(), j);
					
					new_parental_fpkm_cov(next_cov_row,next_cov_col) = _parental_fpkm_covariance(i, j);
					new_parental_fpkm_cov(next_cov_row,next_cov_col+num_kept) = _parental_fpkm_covariance(i, j+to_keep.size());
					new_parental_fpkm_cov(next_cov_row+num_kept,next_cov_col+num_kept) = _parental_fpkm_covariance(i+to_keep.size(), j+to_keep.size());
					new_parental_fpkm_cov(next_cov_row+num_kept,next_cov_col) = _parental_fpkm_covariance(i+to_keep.size(), j);
					
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
//            for (size_t j = 0; j < _paternal_fpkm_samples.size(); ++j)
//            {
//                new_paternal_fpkm_samples[j][curr_abundance_idx] = _paternal_fpkm_samples[j][i];
//            }
//            curr_abundance_idx++;
//        }
//        
//    }
//    curr_abundance_idx = 0;
//    for (size_t i = 0; i < _abundances.size(); ++i)
//    {
//        if (to_keep[i])
//		{
//            for (size_t j = 0; j < _maternal_fpkm_samples.size(); ++j)
//            {
//                new_maternal_fpkm_samples[j][curr_abundance_idx] = _maternal_fpkm_samples[j][i];
//            }
//            curr_abundance_idx++;
//        }
//        
//    }   
   
    
	filtered_group = AlleleAbundanceGroup(new_ab, 
                                    new_parental_cov, 
                                    new_iterated_em_parental_count_cov, 
                                    new_parental_count_cov, 
                                    new_parental_fpkm_cov,
                                    _read_group_props);
    //assert (filtered_group.paternal_FPKM() == 0 || new_paternal_fpkm_samples.size() > 0);
	//assert (filtered_group.maternal_FPKM() == 0 || new_maternal_fpkm_samples.size() > 0);
	
    filtered_group.description(_description);	
}

void AlleleAbundanceGroup::get_transfrags(vector<shared_ptr<Abundance> >& transfrags) const
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

set<string> AlleleAbundanceGroup::gene_id() const	
{
	set<string> s;
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_id();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}

set<string> AlleleAbundanceGroup::gene_name() const	
{
	set<string> s;
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->gene_name();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}


set<string> AlleleAbundanceGroup::tss_id() const	
{
	set<string> s;

	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->tss_id();
		s.insert(sub.begin(), sub.end());
	}

	return s;
}

set<string> AlleleAbundanceGroup::protein_id() const	
{
	set<string> s;
	
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		set<string> sub = pA->protein_id();
		s.insert(sub.begin(), sub.end());
	}
	
	return s;
}

const string& AlleleAbundanceGroup::locus_tag() const	
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
	//assert (false);
	return default_locus_tag;
}

const string& AlleleAbundanceGroup::reference_tag() const	
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
	//assert (false);
	return default_reference_tag;
}

double AlleleAbundanceGroup::paternal_effective_length() const
{
	double eff_len = 0.0;
	double group_fpkm = paternal_FPKM();
	if (group_fpkm == 0)
		return 0;
	BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
	{
		eff_len += (ab->paternal_effective_length() * (ab->paternal_FPKM() / group_fpkm));
	}
	return eff_len;
}

double AlleleAbundanceGroup::maternal_effective_length() const
{
	double eff_len = 0.0;
	double group_fpkm = maternal_FPKM();
	if (group_fpkm == 0)
		return 0;
	BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
	{
		eff_len += (ab->maternal_effective_length() * (ab->maternal_FPKM() / group_fpkm));
	}
	return eff_len;
}

void AlleleAbundanceGroup::collect_per_replicate_mass(const vector<MateHit>& alignments,
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
	compute_cond_probs_and_effective_lengths_allele(alignments, transcripts, mapped_transcripts);
    
	for (size_t i = 0; i < M; ++i)
	{	
		if (!alignments[i].left_alignment())
			continue;
		
		bool mapped = false;
		for (size_t j = 0; j < N; ++j)
        {
			if (_abundances[j]->paternal_cond_probs()->at(i) > 0 || _abundances[j]->maternal_cond_probs()->at(i) > 0)
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
			double curr_mass = inserted.first->second;
			assert (isnan(more_mass) == false);
            inserted.first->second += more_mass;
        }
    }
}

void AlleleAbundanceGroup::calculate_locus_scaled_mass_and_variance(const vector<shared_ptr<Abundance> >& transcripts)
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
    vector<double> avg_mass_paternal_variances(N, 0.0);
	vector<double> avg_mass_maternal_variances(N, 0.0);
    
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
            double scaled_paternal_variance,scaled_maternal_variance;
            //scaled_variance = disperser->scale_mass_variance(scaled_mass * _abundances[j]->gamma());
            scaled_paternal_variance = _abundances[j]->paternal_gamma() * disperser->scale_mass_variance(scaled_mass);
			scaled_maternal_variance = _abundances[j]->maternal_gamma() * disperser->scale_mass_variance(scaled_mass);
            avg_mass_paternal_variances[j] += scaled_paternal_variance;
			avg_mass_maternal_variances[j] += scaled_maternal_variance;
        }
        assert (disperser->scale_mass_variance(scaled_mass) != 0 || scaled_mass == 0); 
        
        //assert (scaled_total_mass != 0.0);
        avg_mass_fraction += (scaled_mass / scaled_total_mass);
    }
    
    double num_replicates = _count_per_replicate.size();
    
    if (num_replicates)
    {
        avg_X_g /= num_replicates;
        avg_mass_fraction /= num_replicates;
        for (size_t j = 0; j < N; ++j)
        {
            avg_mass_paternal_variances[j] /= num_replicates;
			avg_mass_maternal_variances[j] /= num_replicates;
        }
    }
    
    for (size_t j = 0; j < _abundances.size(); ++j)
	{
		_abundances[j]->num_paternal_fragments(_abundances[j]->paternal_gamma() * avg_X_g);
		_abundances[j]->num_maternal_fragments(_abundances[j]->maternal_gamma() * avg_X_g);
        
        double j_avg_mass_paternal_fraction = _abundances[j]->paternal_gamma() * avg_mass_fraction;
		double j_avg_mass_maternal_fraction = _abundances[j]->maternal_gamma() * avg_mass_fraction;
        
        _abundances[j]->paternal_mass_variance(avg_mass_paternal_variances[j]);
		_abundances[j]->maternal_mass_variance(avg_mass_maternal_variances[j]);
        
        if (j_avg_mass_paternal_fraction > 0)
        {
            double paternal_FPKM = j_avg_mass_paternal_fraction * 1000000000/ _abundances[j]->paternal_effective_length();
            paternal_FPKM *= 1.0 / external_scale_factor;
			_abundances[j]->paternal_FPKM(paternal_FPKM);
        }
        else 
        {
            _abundances[j]->paternal_FPKM(0);
            _abundances[j]->paternal_mass_variance(0);
        }
		if (j_avg_mass_maternal_fraction > 0)
        {
            double maternal_FPKM = j_avg_mass_maternal_fraction * 1000000000/ _abundances[j]->maternal_effective_length();
            maternal_FPKM *= 1.0 / external_scale_factor;
			_abundances[j]->maternal_FPKM(maternal_FPKM);
        }
        else 
        {
            _abundances[j]->maternal_FPKM(0);
            _abundances[j]->maternal_mass_variance(0);
        }
	}

}

int total_paternal_cond_prob_calls = 0;
int total_maternal_cond_prob_calls = 0;

//new implementation
void collapse_equivalent_hits_allele(const vector<MateHit>& alignments,
                              vector<shared_ptr<Abundance> >& transcripts,
                              vector<MateHit>& nr_alignments,
                              vector<double>& log_conv_factors, 
                              bool require_overlap = true)
{
	
    int N = transcripts.size();
	int M = alignments.size();
    
    nr_alignments.clear();
    
	vector<vector<char> > paternal_compatibilities(N, vector<char>(M,0));
	vector<vector<char> > maternal_compatibilities(N, vector<char>(M,0));
	compute_compatibilities(transcripts, alignments, paternal_compatibilities, maternal_compatibilities);
    
	//vector<vector<double> > cached_cond_probs (M, vector<double>());
    vector<vector<double> > paternal_cached_cond_probs (M, vector<double>());
	vector<vector<double> > maternal_cached_cond_probs (M, vector<double>());
    
    vector<bool> replaced(M, false);
    int num_replaced = 0;
    
    vector<BiasCorrectionHelper> bchs;
    for (size_t j = 0; j < N; ++j)
    {
        bchs.push_back(BiasCorrectionHelper(transcripts[j]->transfrag()));   
    }
    
    double total_mass_pre_collapse = 0.0;
    
    for(int i = 0 ; i < M; ++i)
    {
        total_mass_pre_collapse += alignments[i].collapse_mass();
        
		//vector<double> cond_probs_i(N,0);
        vector<double> paternal_cond_probs_i(N,0);
		vector<double> maternal_cond_probs_i(N,0);
        if (replaced[i] == true)
            continue;
		if (paternal_cached_cond_probs[i].empty() || maternal_cached_cond_probs[i].empty()) //equivalent to if (cached_cond_probs[i].empty())
        {
            for (int j = 0; j < N; ++j)
            {
                shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
                
                if (paternal_compatibilities[j][i]==1)
                {
                    total_paternal_cond_prob_calls++;
                    paternal_cond_probs_i[j] = bchs[j].get_cond_prob(alignments[i]);
					//cond_probs_i[j] = bchs[j].get_cond_prob(alignments[i]);
                }
				if (maternal_compatibilities[j][i]==1)
                {
                    total_maternal_cond_prob_calls++;
                    maternal_cond_probs_i[j] = bchs[j].get_cond_prob(alignments[i]);
					//cond_probs_i[j] = bchs[j].get_cond_prob(alignments[i]);
                }
                
            }
			//cached_cond_probs[i] = cond_probs_i;
            paternal_cached_cond_probs[i] = paternal_cond_probs_i;
			maternal_cached_cond_probs[i] = maternal_cond_probs_i;
        }
        else
        {
			//cond_probs_i = cached_cond_probs[i];
            paternal_cond_probs_i = paternal_cached_cond_probs[i];
			maternal_cond_probs_i = maternal_cached_cond_probs[i];
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
			if (require_overlap && (!::overlap_in_genome(curr_align->left(), curr_align->right(),
														 alignments[k].left(), alignments[k].right())))
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
            
			//vector<double>* cond_probs_k;
            //double last_cond_prob = -1;
            vector<double>* paternal_cond_probs_k;
            double last_paternal_cond_prob = -1;
			vector<double>* maternal_cond_probs_k;
            double last_maternal_cond_prob = -1;
            
            bool equiv = true;
            
			if (paternal_cached_cond_probs[k].empty() || maternal_cached_cond_probs[k].empty())// equivalent to if (cached_cond_probs[k].empty())
            {
				//cached_cond_probs[k] = vector<double>(N, 0.0);
                //cond_probs_k = &cached_cond_probs[k];
				paternal_cached_cond_probs[k] = vector<double>(N, 0.0);
                paternal_cond_probs_k = &paternal_cached_cond_probs[k];
                maternal_cached_cond_probs[k] = vector<double>(N, 0.0);
                maternal_cond_probs_k = &maternal_cached_cond_probs[k];
                for (int j = 0; j < N; ++j)
                {
                    shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
                    
                    if (paternal_compatibilities[j][k]==1)
                    {
                        total_paternal_cond_prob_calls++;
                        (*paternal_cond_probs_k)[j] = bchs[j].get_cond_prob(alignments[k]);
						//(*cond_probs_k)[j] = bchs[j].get_cond_prob(alignments[k]);
                    }
					if (maternal_compatibilities[j][k]==1)
                    {
                        total_maternal_cond_prob_calls++;
                        (*maternal_cond_probs_k)[j] = bchs[j].get_cond_prob(alignments[k]);
						//(*cond_probs_k)[j] = bchs[j].get_cond_prob(alignments[k]);
                    }
                }
				//paternal_cached_cond_probs[k] = paternal_cond_probs_k;
				//maternal_cached_cond_probs[k] = maternal_cond_probs_k;
            }
            else
            {
				//cond_probs_k = &cached_cond_probs[k];
				paternal_cond_probs_k = &paternal_cached_cond_probs[k];
                maternal_cond_probs_k = &maternal_cached_cond_probs[k];
            }
			
            for (int j = 0; j < N; ++j)
            {
				if (((*paternal_cond_probs_k)[j] != 0 && paternal_cond_probs_i[j] != 0) && ((*maternal_cond_probs_k)[j] != 0 && maternal_cond_probs_i[j] != 0))
                {
					//double cp_j = (*cond_probs_k)[j];
                    //double cp_i = cond_probs_i[j];
                    //double ratio =  cp_j / cp_i;
                    double paternal_cp_j = (*paternal_cond_probs_k)[j];
                    double paternal_cp_i = paternal_cond_probs_i[j];
                    double paternal_ratio =  paternal_cp_j / paternal_cp_i;
					double maternal_cp_j = (*maternal_cond_probs_k)[j];
                    double maternal_cp_i = maternal_cond_probs_i[j];
                    double maternal_ratio =  maternal_cp_j / maternal_cp_i;
					
					if (last_paternal_cond_prob == -1 && last_maternal_cond_prob == -1) //equivalent if (last_cond_prob == -1)
					{
                        //assert(paternal_ratio < 5 && maternal_ratio < 5);
						//last_cond_prob = ratio;
                        last_paternal_cond_prob = paternal_ratio;
						last_maternal_cond_prob = maternal_ratio;
                    }
                    else
                    {
                        if (last_paternal_cond_prob != paternal_ratio || last_maternal_cond_prob != maternal_ratio)
                        //if (abs(last_paternal_cond_prob - paternal_ratio) > 0.001 || abs(last_maternal_cond_prob - maternal_ratio) > 0.001)
                        {
                            equiv = false;
                            break;
                        }
                    }
				}
				else if (((*paternal_cond_probs_k)[j] != 0 && paternal_cond_probs_i[j] != 0) && ((*maternal_cond_probs_k)[j] == 0 && maternal_cond_probs_i[j] == 0))
                {
					//double cp_j = (*cond_probs_k)[j];
                    //double cp_i = cond_probs_i[j];
                    //double ratio =  cp_j / cp_i;
                    double paternal_cp_j = (*paternal_cond_probs_k)[j];
                    double paternal_cp_i = paternal_cond_probs_i[j];
                    double paternal_ratio =  paternal_cp_j / paternal_cp_i;
					if (last_paternal_cond_prob == -1 && last_maternal_cond_prob == -1) //equivalent if (last_cond_prob == -1)
                    {
                        //assert(paternal_ratio < 5);
						//last_cond_prob = ratio;
                        last_paternal_cond_prob = paternal_ratio;
					}
                    else
                    {
                        if (last_paternal_cond_prob != paternal_ratio)
                        //if (abs(last_paternal_cond_prob - paternal_ratio) > 0.001)
                        {
                            equiv = false;
                            break;
                        }
                    }
                }
				else if (((*paternal_cond_probs_k)[j] != 0 && paternal_cond_probs_i[j] != 0) && ((*maternal_cond_probs_k)[j] == 0 && maternal_cond_probs_i[j] == 0))
                {
					//double cp_j = (*cond_probs_k)[j];
                    //double cp_i = cond_probs_i[j];
                    //double ratio =  cp_j / cp_i;
                    double maternal_cp_j = (*maternal_cond_probs_k)[j];
                    double maternal_cp_i = maternal_cond_probs_i[j];
                    double maternal_ratio =  maternal_cp_j / maternal_cp_i;
					if (last_paternal_cond_prob == -1 && last_maternal_cond_prob == -1) //equivalent if (last_cond_prob == -1)
                    {
                        //assert(maternal_ratio < 5);
						//last_cond_prob = ratio;
                        last_maternal_cond_prob = maternal_ratio;
					}
                    else
                    {
                        if (last_maternal_cond_prob != maternal_ratio)
                        //if (abs(last_maternal_cond_prob - maternal_ratio) > 0.001)
                        {
                            equiv = false;
                            break;
                        }
                    }
                }
                else if (((*paternal_cond_probs_k)[j] == 0 && paternal_cond_probs_i[j] == 0) && ((*maternal_cond_probs_k)[j] == 0 && maternal_cond_probs_i[j] == 0))
                {
                    // just do nothing in this iter.
					// last_cond_prob = 0.0;
                    // last_paternal_cond_prob = 0.0;
					// last_maternal_cond_prob = 0.0;
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
                if (last_paternal_cond_prob > 0.0 || last_maternal_cond_prob > 0.0)
                {
                    //assert(curr_align->read_group_props() == alignments[k].read_group_props());
                    assert (last_paternal_cond_prob > 0 || last_maternal_cond_prob > 0);
                    //double paternal_mass_muliplier = sqrt(last_paternal_cond_prob);
                    double paternal_mass_multiplier = log(last_paternal_cond_prob);
                    //assert(last_paternal_cond_prob < 5);
					//double maternal_mass_muliplier = sqrt(last_maternal_cond_prob);
                    double maternal_mass_multiplier = log(last_maternal_cond_prob);
                    //assert(last_maternal_cond_prob < 5);
                    assert (!isinf(paternal_mass_multiplier) && !isnan(paternal_mass_multiplier));
					assert (!isinf(maternal_mass_multiplier) && !isnan(maternal_mass_multiplier));
                    log_conv_factors[log_conv_factors.size() - 1] += paternal_mass_multiplier+maternal_mass_multiplier; 
                    replaced[k] = true;
					//cached_cond_probs[k].clear();
                    paternal_cached_cond_probs[k].clear();
					maternal_cached_cond_probs[k].clear();
					//vector<double>(cached_cond_probs[k]).swap(cached_cond_probs[k]);
                    vector<double>(paternal_cached_cond_probs[k]).swap(paternal_cached_cond_probs[k]);
					vector<double>(maternal_cached_cond_probs[k]).swap(maternal_cached_cond_probs[k]);
                    num_replaced++;
                    double more_mass = alignments[k].collapse_mass();
                    curr_align->incr_collapse_mass(more_mass);
                }
                else
                {
                    replaced[k] = true;
                    num_replaced++;
					//cached_cond_probs[k].clear();
                    paternal_cached_cond_probs[k].clear();
					maternal_cached_cond_probs[k].clear();
					//vector<double>(cached_cond_probs[k]).swap(cached_cond_probs[k]);
                    vector<double>(paternal_cached_cond_probs[k]).swap(paternal_cached_cond_probs[k]);
					vector<double>(maternal_cached_cond_probs[k]).swap(maternal_cached_cond_probs[k]);
                }
            }            
        }
    }
    
    N = transcripts.size();
	//M = nr_alignments.size();
        
	for (int j = 0; j < N; ++j) 
    {
		shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
		//vector<double>& cond_probs = *(new vector<double>(nr_alignments.size(),0));
		vector<double>& paternal_cond_probs = *(new vector<double>(nr_alignments.size(),0));
		vector<double>& maternal_cond_probs = *(new vector<double>(nr_alignments.size(),0));
		
		BiasCorrectionHelper& bch = bchs[j];
		
        size_t last_cond_prob_idx = 0;
		for(int i = 0 ; i < M; ++i)
		{
            if (!paternal_cached_cond_probs[i].empty() || !maternal_cached_cond_probs[i].empty())
            {
                if (paternal_compatibilities[j][i]==1 || maternal_compatibilities[j][i]==1)
                {
                    assert (paternal_cached_cond_probs[i].size() > j || maternal_cached_cond_probs[i].size() > j);
					//cond_probs[last_cond_prob_idx] = cached_cond_probs[i][j];
                    paternal_cond_probs[last_cond_prob_idx] = paternal_cached_cond_probs[i][j];
					maternal_cond_probs[last_cond_prob_idx] = maternal_cached_cond_probs[i][j];
                }
                last_cond_prob_idx++;
            }
        }
		
        assert (last_cond_prob_idx == nr_alignments.size());
        double bch_eff_len = bch.get_effective_length();
		if(transcripts[j]->paternal_effective_length() == transcripts[j]->maternal_effective_length()){
			transcripts[j]->paternal_effective_length(bch_eff_len);
			transcripts[j]->maternal_effective_length(bch_eff_len);
		}
		else{//update according to their proportions
			if(transcripts[j]->paternal_effective_length() > transcripts[j]->maternal_effective_length()){
				transcripts[j]->paternal_effective_length(bch_eff_len);
				transcripts[j]->maternal_effective_length((transcripts[j]->maternal_effective_length()/transcripts[j]->paternal_effective_length())*bch_eff_len);
			}
			else{
				transcripts[j]->paternal_effective_length((transcripts[j]->paternal_effective_length()/transcripts[j]->maternal_effective_length())*bch_eff_len);
				transcripts[j]->maternal_effective_length(bch_eff_len);
			}
		}
		transcripts[j]->paternal_cond_probs(&paternal_cond_probs);
		transcripts[j]->maternal_cond_probs(&maternal_cond_probs);
	}
    
    double total_mass_post_collapse = 0.0;
    for(int i = 0 ; i < nr_alignments.size(); ++i)
    {
        total_mass_post_collapse += nr_alignments[i].collapse_mass();
    }
    
    //assert(abs(total_mass_pre_collapse - total_mass_post_collapse) < 1);
    
    if (nr_alignments.size())
    {
        verbose_msg("\nReduced %lu frags to %lu (%lf percent)\n", alignments.size(), nr_alignments.size(), 100.0 * (1 - nr_alignments.size()/(double)alignments.size()));
    }
}

void collapse_equivalent_hits_helper_allele(const vector<MateHit>& alignments,
                                     vector<shared_ptr<Abundance> >& transcripts,
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
        collapse_equivalent_hits_allele(alignments,
                                 transcripts,
                                 nr_alignments,
                                 log_conv_factors, 
                                 true);
        return;
    }
    
    vector<vector<const MateHit*> > compat_table(1 << N);
    vector<vector<char> > paternal_compatibilities(N, vector<char>(M,0));
	vector<vector<char> > maternal_compatibilities(N, vector<char>(M,0));
    compute_compatibilities(transcripts, alignments, paternal_compatibilities, maternal_compatibilities);
    
    for(int i = 0; i < M; ++i)
    {
        size_t compat_mask = 0;
        for (int j = 0; j < N; ++j)
        {
            compat_mask |= ((paternal_compatibilities[j][i] !=0 || maternal_compatibilities[j][i] !=0) << j);
        }
        assert (compat_mask < compat_table.size());
        compat_table[compat_mask].push_back(&(alignments[i]));
    }
    
    for (size_t i = 0; i < compat_table.size(); ++i)
    {
        vector<MateHit> tmp_hits;
        vector<MateHit> tmp_nr_hits;
        vector<double> tmp_log_conv_factors;
        
        for (size_t j = 0; j < compat_table[i].size(); ++j)
        {
            tmp_hits.push_back(*(compat_table[i][j]));
        }
        if (tmp_hits.empty())
            continue;
        collapse_equivalent_hits_allele(tmp_hits,
                                 transcripts,
                                 tmp_nr_hits,
                                 tmp_log_conv_factors, 
                                 false);
        copy(tmp_nr_hits.begin(), tmp_nr_hits.end(), back_inserter(nr_alignments));
        copy(tmp_log_conv_factors.begin(), tmp_log_conv_factors.end(), back_inserter(log_conv_factors));
    }
}

AbundanceStatus bootstrap_gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
						 const vector<MateHit>& nr_alignments,
						 const vector<double>& log_conv_factors,
						 ublas::vector<double>& paternal_gamma_map_estimate,
						 ublas::vector<double>& maternal_gamma_map_estimate,
						 ublas::matrix<double>& parental_gamma_covariance,
						 double& cross_replicate_js)
{
    size_t N = transcripts.size();
	size_t M = nr_alignments.size();
    
    if (M == 0)
    {
        paternal_gamma_map_estimate = ublas::vector<double>(1);
        paternal_gamma_map_estimate(0) = 0.0;
		maternal_gamma_map_estimate = ublas::vector<double>(1);
        maternal_gamma_map_estimate(0) = 0.0;
        parental_gamma_covariance = ublas::matrix<double>(2,2);
        parental_gamma_covariance(0,0) = 0.0;
		parental_gamma_covariance(0,1) = 0.0;
		parental_gamma_covariance(1,0) = 0.0;
		parental_gamma_covariance(1,1) = 0.0;
		return NUMERIC_OK;
	}
        
    vector<MateHit> alignments = nr_alignments;
    vector<double>  scaled_masses;
    vector<double>  unscaled_masses;
    double num_uncollapsed_frags = 0.0;
    for (size_t i = 0; i < M; ++i)
    {
        double uncollapsed_mass = alignments[i].collapse_mass();
        num_uncollapsed_frags += (uncollapsed_mass);
        scaled_masses.push_back(alignments[i].collapse_mass());
        unscaled_masses.push_back(uncollapsed_mass);
        alignments[i].collapse_mass(uncollapsed_mass);
    }
    
    // FIXME: this has already been computed above, so just pass it in.
    vector<double> paternal_orig_gammas(0.0, transcripts.size());
	vector<double> maternal_orig_gammas(0.0, transcripts.size());
    gamma_mle(transcripts,
              nr_alignments,
              log_conv_factors,
              paternal_orig_gammas,
			  maternal_orig_gammas,
              false);
    
    std::vector<ublas::vector<double> > paternal_mle_gammas,maternal_mle_gammas;
    
    boost::uniform_int<> uniform_dist(0,num_uncollapsed_frags-1);
    boost::mt19937 rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniform_gen(rng, uniform_dist);
    
    int num_sample_frags = floor(num_uncollapsed_frags * bootstrap_fraction);
    
    if (num_sample_frags <= 0)
    {
		return NUMERIC_FAIL;
	}
    
    for (size_t i = 0; i < num_bootstrap_samples; ++i)
    {
        vector<int> sample_idxs;
        for (size_t j = 0; j < num_sample_frags; ++j)
        {
            sample_idxs.push_back(uniform_gen());
        }
        sort (sample_idxs.begin(), sample_idxs.end());
        assert (sample_idxs.empty() == false);
        
        size_t curr_sample = 0;
        size_t processed_hits = 0;
        vector<double> adjusted_masses(alignments.size(), 0);
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            int adjusted_mass = 0.0;
            while (curr_sample < sample_idxs.size() &&
                   sample_idxs[curr_sample] >= processed_hits &&
                   sample_idxs[curr_sample] < processed_hits + alignments[j].collapse_mass())
            {
                adjusted_mass++;
                curr_sample++;
            }
            processed_hits += alignments[j].collapse_mass();
            alignments[j].collapse_mass(adjusted_mass);
            adjusted_masses[j] = adjusted_mass;
        }
        
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            alignments[j].collapse_mass(alignments[j].collapse_mass());
        }
        
        vector<double> paternal_bs_gammas(0.0, transcripts.size());
		vector<double> maternal_bs_gammas(0.0, transcripts.size());
        
//        AbundanceStatus mle_success = gamma_mle(transcripts,
//                                                alignments,
//                                                log_conv_factors,
//                                                paternal_bs_gammas,
//		                                          maternal_bs_gammas,
//                                                false,
//                                                &paternal_orig_gammas,
//		                                          &maternal_orig_gammas);
        
        AbundanceStatus mle_success = gamma_mle(transcripts,
                                                alignments,
                                                log_conv_factors,
                                                paternal_bs_gammas,
												maternal_bs_gammas,
                                                false);
        if (mle_success == NUMERIC_OK)
        {
            ublas::vector<double> paternal_mle = ublas::zero_vector<double>(N);
			ublas::vector<double> maternal_mle = ublas::zero_vector<double>(N);
            for(size_t j = 0; j < N; ++j)
            {
                paternal_mle(j) = paternal_bs_gammas[j];
				maternal_mle(j) = maternal_bs_gammas[j];
            }
            paternal_mle_gammas.push_back(paternal_mle);
			maternal_mle_gammas.push_back(maternal_mle);
        }
        
        for (size_t j = 0; j < alignments.size(); ++j)
        {
            alignments[j].collapse_mass(unscaled_masses[j]);
        }
    }
    
    //fprintf(stderr, "Ran %lu bootstrap samples succesfully\n", paternal_mle_gammas.size());
    
    if (paternal_mle_gammas.empty() || maternal_mle_gammas.empty()){
		return NUMERIC_FAIL;
	}
    
    parental_gamma_covariance = ublas::zero_matrix<double>(2*N,2*N);
    ublas::vector<double> paternal_expected_mle_gamma = ublas::zero_vector<double>(N);
	ublas::vector<double> maternal_expected_mle_gamma = ublas::zero_vector<double>(N);
    
    BOOST_FOREACH(ublas::vector<double>& paternal_mle, paternal_mle_gammas)
    {
        //cerr << "PATERNAL MLE # "<< MLENUM++ << endl;
        //cerr << paternal_mle << endl;
        paternal_expected_mle_gamma += paternal_mle;
    }
    paternal_expected_mle_gamma /= paternal_mle_gammas.size();
	BOOST_FOREACH(ublas::vector<double>& maternal_mle, maternal_mle_gammas)
    {
        //cerr << "MATERNAL MLE # "<< MLENUM++ << endl;
        //cerr << maternal_mle << endl;
        maternal_expected_mle_gamma += maternal_mle;
    }
    maternal_expected_mle_gamma /= maternal_mle_gammas.size();
    
	assert(paternal_mle_gammas.size() == maternal_mle_gammas.size());
	
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            for (size_t k = 0 ; k < paternal_mle_gammas.size(); ++k)
            {
                double cpp = (paternal_mle_gammas[k](i) - paternal_expected_mle_gamma(i)) * (paternal_mle_gammas[k](j) - paternal_expected_mle_gamma(j));
                parental_gamma_covariance(i,j) += cpp;
				double cpm = (paternal_mle_gammas[k](i) - paternal_expected_mle_gamma(i)) * (maternal_mle_gammas[k](j) - maternal_expected_mle_gamma(j));
				parental_gamma_covariance(i,j+N) += cpm;
				double cmm = (maternal_mle_gammas[k](i) - maternal_expected_mle_gamma(i)) * (maternal_mle_gammas[k](j) - maternal_expected_mle_gamma(j));
                parental_gamma_covariance(i+N,j+N) += cmm;
				double cmp = (maternal_mle_gammas[k](i) - maternal_expected_mle_gamma(i)) * (paternal_mle_gammas[k](j) - paternal_expected_mle_gamma(j));
                parental_gamma_covariance(i+N,j) += cmp;
            }
        }
    }
    
    parental_gamma_covariance /= (paternal_mle_gammas.size()+maternal_mle_gammas.size());
    paternal_gamma_map_estimate = paternal_expected_mle_gamma;
	maternal_gamma_map_estimate = maternal_expected_mle_gamma;
    
    //cerr << "PATERNAL MLE: " << paternal_expected_mle_gamma << endl;
	//cerr << "MATERNAL MLE: " << maternal_expected_mle_gamma << endl;
    //cerr << "COV:" << endl;
    //cerr << parental_gamma_covariance << endl;
    //cerr << "*************" << endl;
	return NUMERIC_OK;
}

AbundanceStatus bootstrap_gammas(const vector<shared_ptr<Abundance> >& transcripts,
					   const vector<MateHit>& alignments,
					   const vector<double>& log_conv_factors,
					   ublas::vector<double>& paternal_gamma_estimate,
					   ublas::vector<double>& maternal_gamma_estimate,
					   ublas::matrix<double>& parental_gamma_covariance,
					   double& cross_replicate_js)
{
    ublas::vector<double> paternal_empirical_gamma_mle = paternal_gamma_estimate;
	ublas::vector<double> maternal_empirical_gamma_mle = maternal_gamma_estimate;
    ublas::matrix<double> parental_empirical_gamma_covariance = parental_gamma_covariance;
    
    // Calculate the mean gamma MLE and covariance matrix across replicates, so
    // we can use it as the proposal distribution for importance sampling.  This will
    // make the Bayesian prior more conservative than using the inverse of the
    // Fisher Information matrix on the mixed likelihood function.
    AbundanceStatus empirical_mle_status = bootstrap_gamma_mle(transcripts,
                                                               alignments,
                                                               log_conv_factors,
                                                               paternal_empirical_gamma_mle,
															   maternal_empirical_gamma_mle,
                                                               parental_empirical_gamma_covariance,
                                                               cross_replicate_js);
    
    if (empirical_mle_status != NUMERIC_OK)
        return empirical_mle_status;
    
    paternal_gamma_estimate = paternal_empirical_gamma_mle;
	maternal_gamma_estimate = maternal_empirical_gamma_mle;
    parental_gamma_covariance = parental_empirical_gamma_covariance;
    
    return NUMERIC_OK;
}

bool generate_count_assignment_samples(int num_draws,
                                       const vector<double>& paternal_count_mean,
									   const vector<double>& maternal_count_mean,
                                       const ublas::matrix<double>& count_covariance,
                                       vector<ublas::vector<double> >& paternal_assigned_count_samples,
									   vector<ublas::vector<double> >& maternal_assigned_count_samples)
{
	
    double total_paternal_frag_counts = accumulate(paternal_count_mean.begin(), paternal_count_mean.end(), 0.0);
	double total_maternal_frag_counts = accumulate(maternal_count_mean.begin(), maternal_count_mean.end(), 0.0);
    
	bool paternal_pass = true;
	bool maternal_pass = true;

    if (total_paternal_frag_counts == 0)
    {
        paternal_assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(paternal_count_mean.size()));
        paternal_pass = false;
	}
	if (total_maternal_frag_counts == 0)
    {
        maternal_assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(maternal_count_mean.size()));
        maternal_pass = false;
	}
    if(!paternal_pass && !maternal_pass) return true;
	 
    boost::mt19937 rng;
    
	ublas::vector<double> paternal_mle_frag_counts = ublas::zero_vector<double>(paternal_count_mean.size());
	ublas::vector<double> maternal_mle_frag_counts = ublas::zero_vector<double>(maternal_count_mean.size());
    
    for (size_t j = 0; j < paternal_count_mean.size(); ++j)
    {
        paternal_mle_frag_counts(j) = paternal_count_mean[j];
    }
	for (size_t j = 0; j < maternal_count_mean.size(); ++j)
    {
        maternal_mle_frag_counts(j) = maternal_count_mean[j];
    }
	
    //    cerr << "********************" << endl;
    //    cerr << "initial MLE counts: " << mle_frag_counts << endl;
    
    ublas::matrix<double> parental_mle_count_covar = count_covariance;
	
    ublas::matrix<double> parental_epsilon = ublas::zero_matrix<double>(paternal_count_mean.size()+maternal_count_mean.size(),paternal_count_mean.size()+maternal_count_mean.size());
    for (size_t i = 0; i < paternal_count_mean.size()+maternal_count_mean.size(); ++i)
    {
        parental_epsilon(i,i) = 1e-6;
    }
	parental_mle_count_covar += parental_epsilon; // modify matrix to avoid problems during inverse
	double ret = cholesky_factorize(parental_mle_count_covar);
    if (ret != 0)
    {
        fprintf(stderr, "Warning: Iterated expectation parental count covariance matrix cannot be cholesky factorized!\n");
        //fprintf(stderr, "Warning: parental FPKM covariance is not positive definite (ret = %lg)!\n", ret);
        //        for (size_t j = 0; j < _abundances.size(); ++j)
        //        {
        //            _abundances[j]->paternal_status(NUMERIC_FAIL);
		//            _abundances[j]->maternal_status(NUMERIC_FAIL);
        //        }
        return false;
    }
    
    //    cerr << endl << "Cholesky factored parental covariance matrix: " << endl;
    //    for (unsigned i = 0; i < _parental_count_covariance.size1 (); ++ i)
    //    {
    //        ublas::matrix_row<ublas::matrix<double> > mr (parental_mle_count_covar, i);
    //        cerr << i << " : " << _abundances[i]->num_paternal_fragments() << " : "<< _abundances[i]->num_maternal_fragments() << " : ";
    //        std::cerr << i << " : " << mr << std::endl;
    //    }
    //cerr << "======" << endl;
    ublas::vector<double> parental_mle_frag_counts(paternal_mle_frag_counts.size()+maternal_mle_frag_counts.size());
	for(size_t c = 0;c < paternal_mle_frag_counts.size();++c)
		parental_mle_frag_counts[c] = paternal_mle_frag_counts[c];
	for(size_t c = 0;c < maternal_mle_frag_counts.size();++c)
		parental_mle_frag_counts[c+paternal_mle_frag_counts.size()] = maternal_mle_frag_counts[c];
    multinormal_generator<double> generator(parental_mle_frag_counts, parental_mle_count_covar);
	//vector<Eigen::VectorXd> multinormal_samples;
    
    paternal_assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(paternal_count_mean.size()));
	maternal_assigned_count_samples = vector<ublas::vector<double> > (num_draws, ublas::zero_vector<double>(maternal_count_mean.size()));
    	
    boost::uniform_01<> uniform_dist;
    boost::mt19937 null_rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uniform_gen(null_rng, uniform_dist);
			
	for (size_t assign_idx = 0; assign_idx < num_draws; ++assign_idx)
    {
		boost::numeric::ublas::vector<double> random_count_assign = generator.next_rand();
		boost::numeric::ublas::vector<double> paternal_random_count_assign(random_count_assign.size()/2,0);
		boost::numeric::ublas::vector<double> maternal_random_count_assign(random_count_assign.size()/2,0);
		//cerr << random_count_assign << endl;
        
        for (size_t r_idx = 0; r_idx < random_count_assign.size(); ++r_idx)
        {
			if(r_idx < random_count_assign.size()/2){
				paternal_random_count_assign(r_idx) = random_count_assign(r_idx);
				if (paternal_random_count_assign(r_idx) < 0)
					paternal_random_count_assign(r_idx) = 0;
			}
			else{
				maternal_random_count_assign(r_idx-random_count_assign.size()/2) = random_count_assign(r_idx);
				if (maternal_random_count_assign(r_idx-random_count_assign.size()/2) < 0)
					maternal_random_count_assign(r_idx-random_count_assign.size()/2) = 0;
			}			
        }
        		
        double total_paternal_sample_counts = accumulate(paternal_random_count_assign.begin(), paternal_random_count_assign.end(), 0.0);
		double total_maternal_sample_counts = accumulate(maternal_random_count_assign.begin(), maternal_random_count_assign.end(), 0.0);
		if (total_paternal_sample_counts > 0)
			paternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (paternal_random_count_assign / (total_paternal_sample_counts+total_maternal_sample_counts));
        else
            paternal_random_count_assign = ublas::zero_vector<double>(paternal_count_mean.size());
		if (total_maternal_sample_counts > 0)
			maternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (maternal_random_count_assign / (total_paternal_sample_counts+total_maternal_sample_counts));
        else
            maternal_random_count_assign = ublas::zero_vector<double>(maternal_count_mean.size());
		paternal_assigned_count_samples[assign_idx] = paternal_random_count_assign;
		maternal_assigned_count_samples[assign_idx] = maternal_random_count_assign;
	}
    	
    return true;
}

void calculate_gamma_mle_covariance(const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AlleleAbundanceGroup> >& ab_group_per_replicate,
                                    ublas::vector<double>& paternal_estimated_gamma_mean,
									ublas::vector<double>& maternal_estimated_gamma_mean,
                                    ublas::matrix<double>& parental_estimated_gamma_covariance)
{
    vector<ublas::vector<double> > paternal_all_assigned_count_samples,maternal_all_assigned_count_samples;
    
    if (ab_group_per_replicate.empty())
        return;
    int num_transcripts = ab_group_per_replicate.begin()->second->abundances().size();

    for(std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AlleleAbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
        itr != ab_group_per_replicate.end();
        ++itr)
    {
        ublas::vector<double> paternal_count_mean = ublas::zero_vector<double>(itr->second->abundances().size());
		ublas::vector<double> maternal_count_mean = ublas::zero_vector<double>(itr->second->abundances().size());
        for (size_t i = 0; i < itr->second->abundances().size(); ++i)
        {
            paternal_count_mean(i) = itr->second->abundances()[i]->paternal_gamma();
			maternal_count_mean(i) = itr->second->abundances()[i]->maternal_gamma();
        }
        
        paternal_all_assigned_count_samples.push_back(paternal_count_mean);
		maternal_all_assigned_count_samples.push_back(maternal_count_mean);
    }
    
    //    for (size_t i = 0; i < paternal_all_assigned_count_samples.size(); ++i)
    //    {
    //        double paternal_total = accumulate(paternal_all_assigned_count_samples[i].begin(), paternal_all_assigned_count_samples[i].end(), 0.0);
    //        if (paternal_total > 0)
    //            paternal_all_assigned_count_samples[i] /= paternal_total;
    //    }
	//    for (size_t i = 0; i < maternal_all_assigned_count_samples.size(); ++i)
    //    {
    //        double maternal_total = accumulate(maternal_all_assigned_count_samples[i].begin(), maternal_all_assigned_count_samples[i].end(), 0.0);
    //        if (maternal_total > 0)
    //            maternal_all_assigned_count_samples[i] /= maternal_total;
    //    }
    
    paternal_estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
	maternal_estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
    parental_estimated_gamma_covariance = ublas::zero_matrix<double>(2*num_transcripts, 2*num_transcripts);
    assert(paternal_all_assigned_count_samples.size() == maternal_all_assigned_count_samples.size());
	
    for (size_t i = 0; i < paternal_all_assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < paternal_all_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(paternal_all_assigned_count_samples[i](j)) && !isinf(paternal_all_assigned_count_samples[i](j)) && !isnan(maternal_all_assigned_count_samples[i](j)) && !isinf(maternal_all_assigned_count_samples[i](j)));
        }
        
        paternal_estimated_gamma_mean += paternal_all_assigned_count_samples[i];
		maternal_estimated_gamma_mean += maternal_all_assigned_count_samples[i];
        //
        //paternal_expected_relative_abundances += paternal_relative_abundances[i];
		//maternal_expected_relative_abundances += maternal_relative_abundances[i];
    }
    
    if (paternal_all_assigned_count_samples.size() > 1)
    {
        paternal_estimated_gamma_mean /= (paternal_all_assigned_count_samples.size() - 1);
    }
	if (maternal_all_assigned_count_samples.size() > 1)
    {
        maternal_estimated_gamma_mean /= (maternal_all_assigned_count_samples.size() - 1);
    }
    
    for (size_t i = 0; i < num_transcripts; ++i)
    {
        for (size_t j = 0; j < num_transcripts; ++j)
        {
            for (size_t k = 0 ; k < paternal_all_assigned_count_samples.size(); ++k)
            {
                double cpp = (paternal_all_assigned_count_samples[k](i) - paternal_estimated_gamma_mean(i)) * (paternal_all_assigned_count_samples[k](j) - paternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i,j) += cpp;
                
                assert (!isinf(parental_estimated_gamma_covariance(i,j)) && !isnan(parental_estimated_gamma_covariance(i,j)));
				
				double cpm = (paternal_all_assigned_count_samples[k](i) - paternal_estimated_gamma_mean(i)) * (maternal_all_assigned_count_samples[k](j) - maternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i,j+num_transcripts) += cpm;
                
                assert (!isinf(parental_estimated_gamma_covariance(i,j+num_transcripts)) && !isnan(parental_estimated_gamma_covariance(i,j+num_transcripts)));
				
				double cmm = (maternal_all_assigned_count_samples[k](i) - maternal_estimated_gamma_mean(i)) * (maternal_all_assigned_count_samples[k](j) - maternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i+num_transcripts,j+num_transcripts) += cmm;
                
                assert (!isinf(parental_estimated_gamma_covariance(i+num_transcripts,j+num_transcripts)) && !isnan(parental_estimated_gamma_covariance(i+num_transcripts,j+num_transcripts)));
				double cmp = (maternal_all_assigned_count_samples[k](i) - maternal_estimated_gamma_mean(i)) * (paternal_all_assigned_count_samples[k](j) - paternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i+num_transcripts,j) += cmp;
                
                assert (!isinf(parental_estimated_gamma_covariance(i+num_transcripts,j)) && !isnan(parental_estimated_gamma_covariance(i+num_transcripts,j)));
				
            }
        }
    }
    
    if (paternal_all_assigned_count_samples.size() > 1 && maternal_all_assigned_count_samples.size() > 1)
    {
        parental_estimated_gamma_covariance /= (paternal_all_assigned_count_samples.size() - 1 + maternal_all_assigned_count_samples.size() - 1);
    }
}

void calculate_fragment_assignment_distribution(const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >& ab_group_per_replicate,
                                                ublas::vector<double>& paternal_estimated_gamma_mean,
												ublas::vector<double>& maternal_estimated_gamma_mean,
                                                ublas::matrix<double>& parental_estimated_gamma_covariance,
                                                vector<ublas::vector<double> >& paternal_all_assigned_count_samples,
												vector<ublas::vector<double> >& maternal_all_assigned_count_samples)
{
    paternal_all_assigned_count_samples.clear();
	maternal_all_assigned_count_samples.clear();
    
    if (ab_group_per_replicate.empty())
        return;
    int num_transcripts = ab_group_per_replicate.begin()->second->abundances().size();

    for(std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
        itr != ab_group_per_replicate.end();
        ++itr)
    {
        vector<double> paternal_count_mean,maternal_count_mean;
        for (size_t i = 0; i < itr->second->abundances().size(); ++i)
        {
            paternal_count_mean.push_back(itr->second->abundances()[i]->num_paternal_fragments());
			maternal_count_mean.push_back(itr->second->abundances()[i]->num_maternal_fragments());
        }
        
        ublas::matrix<double> parental_count_covariance = itr->second->iterated_count_cov();
				
        ublas::matrix<double> parental_mle_error = ublas::zero_matrix<double>(2*itr->second->abundances().size(), 2*itr->second->abundances().size());
        shared_ptr<const MleErrorModel> parental_mle_model = itr->first->mle_error_model();
			
		//???
        if (parental_mle_model != NULL)
        {
            for (size_t i = 0; i < itr->second->abundances().size(); ++i)
            {
				double paternal_mle_var = parental_mle_model->scale_mle_variance(paternal_count_mean[i]);//assuming that the mle model will be constructed based on allele-specific mean:var this holds
				double maternal_mle_var = parental_mle_model->scale_mle_variance(maternal_count_mean[i]);//assuming that the mle model will be constructed based on allele-specific mean:var this holds
                parental_mle_error(i,i) = max(0.0, paternal_mle_var);
				parental_mle_error(i+itr->second->abundances().size(),i+itr->second->abundances().size()) = max(0.0, maternal_mle_var);
			}
			parental_count_covariance += parental_mle_error;
//            cerr << endl << "MLE error correction: " << endl;
//            for (unsigned i = 0; i < parental_mle_error.size1 (); ++ i)
//            {
//                ublas::matrix_row<ublas::matrix<double> > mr (parental_mle_error, i);
//                if(i < itr->second->abundances().size())
//                     cerr << i << " : " << paternal_count_mean[i] << " : "<< maternal_count_mean[i] << " : ";
//                else
//                     cerr << i << " : " << paternal_count_mean[i-itr->second->abundances().size()] << " : "<< maternal_count_mean[i-itr->second->abundances().size()] << " : ";
//                std::cerr << i << " : " << mr << std::endl;
//            }
//            cerr << "======" << endl;
        }
				
		vector<ublas::vector<double> > paternal_assigned_count_samples;
		vector<ublas::vector<double> > maternal_assigned_count_samples;
		generate_count_assignment_samples(num_frag_assignments,
                                          paternal_count_mean,
										  maternal_count_mean,
                                          parental_count_covariance,
                                          paternal_assigned_count_samples,
										  maternal_assigned_count_samples);
		
        paternal_all_assigned_count_samples.insert(paternal_all_assigned_count_samples.end(), paternal_assigned_count_samples.begin(), paternal_assigned_count_samples.end());
		maternal_all_assigned_count_samples.insert(maternal_all_assigned_count_samples.end(), maternal_assigned_count_samples.begin(), maternal_assigned_count_samples.end());
    }
    assert(paternal_all_assigned_count_samples.size() == maternal_all_assigned_count_samples.size());
	
    for (size_t i = 0; i < paternal_all_assigned_count_samples.size(); ++i)
    {
        double paternal_total = accumulate(paternal_all_assigned_count_samples[i].begin(), paternal_all_assigned_count_samples[i].end(), 0.0);
		double maternal_total = accumulate(maternal_all_assigned_count_samples[i].begin(), maternal_all_assigned_count_samples[i].end(), 0.0);
        if (paternal_total+maternal_total > 0){
			paternal_all_assigned_count_samples[i] /= (paternal_total+maternal_total);
			maternal_all_assigned_count_samples[i] /= (paternal_total+maternal_total);
		}
		//cerr << paternal_all_assigned_count_samples[i] << endl;
		//cerr << maternal_all_assigned_count_samples[i] << endl;
    }
    
    paternal_estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
	maternal_estimated_gamma_mean = ublas::zero_vector<double>(num_transcripts);
    parental_estimated_gamma_covariance = ublas::zero_matrix<double>(2*num_transcripts, 2*num_transcripts);
    
    for (size_t i = 0; i < paternal_all_assigned_count_samples.size(); ++i)
    {
		for (int j = 0; j < paternal_all_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(paternal_all_assigned_count_samples[i](j)) && !isinf(paternal_all_assigned_count_samples[i](j)));
        }
        
        paternal_estimated_gamma_mean += paternal_all_assigned_count_samples[i];
        //
        //paternal_expected_relative_abundances += paternal_relative_abundances[i];
    }
	for (size_t i = 0; i < maternal_all_assigned_count_samples.size(); ++i)
    {
		for (int j = 0; j < maternal_all_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(maternal_all_assigned_count_samples[i](j)) && !isinf(maternal_all_assigned_count_samples[i](j)));
        }
        
		maternal_estimated_gamma_mean += maternal_all_assigned_count_samples[i];
        //
		//maternal_expected_relative_abundances += maternal_relative_abundances[i];
    }
    
    if (paternal_all_assigned_count_samples.size() > 1)
    {
        paternal_estimated_gamma_mean /= (paternal_all_assigned_count_samples.size() - 1);
    }
	if (maternal_all_assigned_count_samples.size() > 1)
    {
        maternal_estimated_gamma_mean /= (maternal_all_assigned_count_samples.size() - 1);
    }
	
    for (size_t i = 0; i < num_transcripts; ++i)
    {
        for (size_t j = 0; j < num_transcripts; ++j)
        {
            for (size_t k = 0 ; k < paternal_all_assigned_count_samples.size(); ++k)
            {
                double cpp = (paternal_all_assigned_count_samples[k](i) - paternal_estimated_gamma_mean(i)) * (paternal_all_assigned_count_samples[k](j) - paternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i,j) += cpp;
                
                assert (!isinf(parental_estimated_gamma_covariance(i,j)) && !isnan(parental_estimated_gamma_covariance(i,j)));
				
				double cpm = (paternal_all_assigned_count_samples[k](i) - paternal_estimated_gamma_mean(i)) * (maternal_all_assigned_count_samples[k](j) - maternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i,j+num_transcripts) += cpm;
                
                assert (!isinf(parental_estimated_gamma_covariance(i,j+num_transcripts)) && !isnan(parental_estimated_gamma_covariance(i,j+num_transcripts)));
				
				double cmm = (maternal_all_assigned_count_samples[k](i) - maternal_estimated_gamma_mean(i)) * (maternal_all_assigned_count_samples[k](j) - maternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(j+num_transcripts,j+num_transcripts) += cmm;
                
                assert (!isinf(parental_estimated_gamma_covariance(i+num_transcripts,j+num_transcripts)) && !isnan(parental_estimated_gamma_covariance(i+num_transcripts,j+num_transcripts)));
				double cmp = (maternal_all_assigned_count_samples[k](i) - maternal_estimated_gamma_mean(i)) * (paternal_all_assigned_count_samples[k](j) - paternal_estimated_gamma_mean(j));
                parental_estimated_gamma_covariance(i+num_transcripts,j) += cmp;
                
                assert (!isinf(parental_estimated_gamma_covariance(i+num_transcripts,j)) && !isnan(parental_estimated_gamma_covariance(i+num_transcripts,j)));
				
            }
        }
    }
    if (paternal_all_assigned_count_samples.size() > 1 && maternal_all_assigned_count_samples.size() > 1)
    {
        parental_estimated_gamma_covariance /= (ab_group_per_replicate.size() - 1);
    }
	
//    cerr << "Gamma covariance: " << endl;
//    for (unsigned i = 0; i < parental_estimated_gamma_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (parental_estimated_gamma_covariance, i);
//        if(i < num_transcripts)
//            cerr << i << " : " << paternal_estimated_gamma_mean[i] << " : "<< maternal_estimated_gamma_mean[i] << " : ";
//        else
//            cerr << i << " : " << paternal_estimated_gamma_mean[i-num_transcripts] << " : "<< maternal_estimated_gamma_mean[i-num_transcripts] << " : ";
//        std::cerr << i << " : " << mr << std::endl;
//    }

}


#define PERFORM_EQUIV_COLLAPSE 1

void AlleleAbundanceGroup::calculate_abundance_group_variance(const vector<shared_ptr<Abundance> >& transcripts,
                                                        const std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >& ab_group_per_replicate)
{
	if (final_est_run) // Only on last estimation run
	{
        // Simulate NB draws and fragment assignment under uncertainty to sample
        // from the BNBs.
        
        ublas::matrix<double> parental_count_assign_covariance;
        
        ublas::vector<double> paternal_estimated_gamma_mean,maternal_estimated_gamma_mean;
        ublas::matrix<double> parental_estimated_gamma_covariance;
        
        vector<ublas::vector<double> > paternal_gamma_samples,maternal_gamma_samples;
        calculate_fragment_assignment_distribution(ab_group_per_replicate,
                                                   paternal_estimated_gamma_mean,
												   maternal_estimated_gamma_mean,
                                                   parental_estimated_gamma_covariance,
                                                   paternal_gamma_samples,
												   maternal_gamma_samples);
		
        ublas::vector<double> paternal_estimated_count_mean = paternal_estimated_gamma_mean * (num_paternal_fragments()+num_maternal_fragments());
		ublas::vector<double> maternal_estimated_count_mean = maternal_estimated_gamma_mean * (num_paternal_fragments()+num_maternal_fragments());
        ublas::matrix<double> parental_estimated_count_covariance = parental_estimated_gamma_covariance * (num_paternal_fragments()+num_maternal_fragments()) * (num_paternal_fragments()+num_maternal_fragments());
        
        //cerr << paternal_estimated_count_mean << endl;
		//cerr << maternal_estimated_count_mean << endl;
        //cerr << parental_estimated_count_covariance << endl;
        vector<double> paternal_frags_per_transcript,maternal_frags_per_transcript;
        vector<double> paternal_frag_variances,maternal_frag_variances;
        
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            assert (paternal_estimated_count_mean.size() > i && maternal_estimated_count_mean.size() >);
            paternal_frags_per_transcript.push_back(paternal_estimated_count_mean[i]);
			maternal_frags_per_transcript.push_back(maternal_estimated_count_mean[i]);
            //paternal_frags_per_transcript.push_back(_abundances[i]->num_paternal_fragments());
			//maternal_frags_per_transcript.push_back(_abundances[i]->num_maternal_fragments());
            
            paternal_frag_variances.push_back(_abundances[i]->paternal_mass_variance());
			maternal_frag_variances.push_back(_abundances[i]->maternal_mass_variance());
        }
        
        simulate_count_covariance(paternal_frags_per_transcript, maternal_frags_per_transcript, paternal_frag_variances, maternal_frag_variances, parental_estimated_count_covariance, transcripts, _parental_count_covariance, _paternal_assigned_count_samples, _maternal_assigned_count_samples, &paternal_gamma_samples, &maternal_gamma_samples);
		generate_fpkm_samples();
		calculate_FPKM_covariance();
		// Derive confidence intervals from the FPKM variance/covariance matrix
        calculate_conf_intervals();
		// Calculate the inter-group relative abundances and variances
        if (no_js_tests == false && _read_group_props.size() >= min_reps_for_js_test)
        {
            calculate_kappas();
        }
        
    }
    
    //cerr << _count_covariance << endl;
//    for (size_t i = 0; i < _abundances.size(); ++i)
//    {
//        for (size_t j = 0; j < _abundances.size(); ++j)
//        {
//            if (i != j)
//            {
//                assert(!isinf(_parental_fpkm_covariance(i,j)) && !isnan(_parental_fpkm_covariance(i,j)));
//				assert(!isinf(_parental_fpkm_covariance(i,j+_abundances.size())) && !isnan(_parental_fpkm_covariance(i,j+_abundances.size())));
//				assert(!isinf(_parental_fpkm_covariance(i+_abundances.size(),j)) && !isnan(_parental_fpkm_covariance(i+_abundances.size(),j)));
//				assert(!isinf(_parental_fpkm_covariance(i+_abundances.size(),j+_abundances.size())) && !isnan(_parental_fpkm_covariance(i+_abundances.size(),j+_abundances.size())));
//                if (_abundances[i]->transfrag()->contains(*_abundances[j]->transfrag()) &&
//                    Scaffold::compatible(*_abundances[i]->transfrag(),*_abundances[j]->transfrag())) 
//                {
//					_abundances[j]->paternal_status(NUMERIC_LOW_DATA);
//					_abundances[j]->maternal_status(NUMERIC_LOW_DATA);
//                }
//            }
//		}
//    }
    
    //assert (paternal_FPKM() == 0 || _paternal_assigned_count_samples.size() > 0);
	//assert (maternal_FPKM() == 0 || _maternal_assigned_count_samples.size() > 0);
    
    //fprintf(stderr, "Total calls to get_cond_prob = %d\n", total_cond_prob_calls);
}

void AlleleAbundanceGroup::calculate_abundance_for_replicate(const vector<MateHit>& alignments, bool perform_collapse)
{
    vector<shared_ptr<Abundance> > transcripts;
    
	get_transfrags(transcripts);
	vector<shared_ptr<Abundance> > mapped_transcripts; // This collects the transcripts that have alignments mapping to them
    
	vector<MateHit> nr_alignments;
    vector<double> paternal_joint_mle_gammas,maternal_joint_mle_gammas;
    
    vector<MateHit> non_equiv_alignments;
    vector<double> log_conv_factors;
    
    if (final_est_run || corr_multi || corr_bias)
    {
        compute_cond_probs_and_effective_lengths_allele(alignments, transcripts, mapped_transcripts);
        collect_per_replicate_mass(alignments, transcripts);
        calculate_gammas(alignments, log_conv_factors, transcripts, mapped_transcripts);
        for (size_t i = 0; i < _abundances.size(); ++i)
        {
            paternal_joint_mle_gammas.push_back(_abundances[i]->paternal_gamma());
			maternal_joint_mle_gammas.push_back(_abundances[i]->maternal_gamma());
		}
    }
		
	calculate_iterated_exp_count_covariance(paternal_joint_mle_gammas, maternal_joint_mle_gammas, alignments, transcripts, _iterated_exp_parental_count_covariance);
	
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        _abundances[i]->num_paternal_fragment_uncertainty_var(_iterated_exp_parental_count_covariance(i,i));
		_abundances[i]->num_maternal_fragment_uncertainty_var(_iterated_exp_parental_count_covariance(i+_abundances.size(),i+_abundances.size()));
    }
    
    if (corr_multi && !final_est_run)
    {
        update_multi_reads(alignments, mapped_transcripts);
    }
    
    // Calculate the initial estimates for the number of fragments originating
    // from each transcript, the FPKMs, and set the NB variances

    calculate_locus_scaled_mass_and_variance(transcripts);
}

void AlleleAbundanceGroup::set_allele_informative()
{
	for(size_t i = 0;i < _abundances.size();++i)
		boost::static_pointer_cast<AlleleTranscriptAbundance>(_abundances[i])->set_allele_informative();
}

void AlleleAbundanceGroup::set_allele_informative(const map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >& ab_group_per_replicate)
{
	for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
		 itr != ab_group_per_replicate.end();
		 ++itr)
	{
		for(size_t i = 0;i < itr->second->abundances().size();++i)
			boost::static_pointer_cast<AlleleTranscriptAbundance>(_abundances[i])->set_allele_informative(boost::static_pointer_cast<AlleleTranscriptAbundance>(itr->second->abundances()[i])->is_allele_informative());
	}
}

bool AlleleAbundanceGroup::is_allele_informative() const
{
	bool result = false;
	for(size_t i = 0;i < _abundances.size();++i)
	{
		if(boost::static_pointer_cast<AlleleTranscriptAbundance>(_abundances[i])->is_allele_informative())
		{
			result = true;
			break;
		}
	}
	return (result);
}


void AlleleAbundanceGroup::aggregate_replicate_abundances(const map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >& ab_group_per_replicate)
{
	for (size_t i = 0; i < _abundances.size(); ++i)
    {
        CountPerReplicateTable paternal_cpr,maternal_cpr;
        FPKMPerReplicateTable paternal_fpr,maternal_fpr;
        StatusPerReplicateTable paternal_spr,maternal_spr;
        
        double avg_paternal_fpkm = 0.0;
		double avg_maternal_fpkm = 0.0;
        double avg_paternal_num_frags = 0.0;
		double avg_maternal_num_frags = 0.0;
        double avg_paternal_gamma = 0.0;
		double avg_maternal_gamma = 0.0;
        double avg_paternal_mass_variance = 0.0;
		double avg_maternal_mass_variance = 0.0;
        double avg_paternal_effective_length = 0.0;
		double avg_maternal_effective_length = 0.0;
        
        map<AbundanceStatus, int> paternal_status_table;
		map<AbundanceStatus, int> maternal_status_table;
        paternal_status_table[NUMERIC_OK] = 0;
		maternal_status_table[NUMERIC_OK] = 0;
        paternal_status_table[NUMERIC_LOW_DATA] = 0;
		maternal_status_table[NUMERIC_LOW_DATA] = 0;
        paternal_status_table[NUMERIC_FAIL] = 0;
		maternal_status_table[NUMERIC_FAIL] = 0;
        paternal_status_table[NUMERIC_HI_DATA] = 0;
		maternal_status_table[NUMERIC_HI_DATA] = 0;
		
		
        for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
             itr != ab_group_per_replicate.end();
             ++itr)
        {
			paternal_status_table[itr->second->abundances()[i]->paternal_status()] += 1;
			maternal_status_table[itr->second->abundances()[i]->maternal_status()] += 1;
        }
        		
        for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
             itr != ab_group_per_replicate.end();
             ++itr)
        {
			const vector<shared_ptr<Abundance> >& sc_ab = itr->second->abundances();
            assert(itr->second->abundances().size() == _abundances.size());
            paternal_cpr[itr->first] = itr->second->abundances()[i]->num_paternal_fragments();
			maternal_cpr[itr->first] = itr->second->abundances()[i]->num_maternal_fragments();
            //fprintf(stderr, "Paternal FPKM = %lg\n", itr->second->abundances()[i]->paternal_FPKM());
			//fprintf(stderr, "Maternal FPKM = %lg\n", itr->second->abundances()[i]->maternal_FPKM());
            paternal_fpr[itr->first] = itr->second->abundances()[i]->paternal_FPKM();
			maternal_fpr[itr->first] = itr->second->abundances()[i]->maternal_FPKM();
            paternal_spr[itr->first] = itr->second->abundances()[i]->paternal_status();
			maternal_spr[itr->first] = itr->second->abundances()[i]->maternal_status();
			/*
			if (itr->second->abundances()[i]->paternal_status() == NUMERIC_OK)
            {
                avg_paternal_fpkm += itr->second->abundances()[i]->paternal_FPKM() / (double)paternal_status_table[NUMERIC_OK];
                avg_paternal_num_frags += itr->second->abundances()[i]->num_paternal_fragments() / (double)paternal_status_table[NUMERIC_OK];
                avg_paternal_gamma += itr->second->abundances()[i]->paternal_gamma() / (double)paternal_status_table[NUMERIC_OK];
                avg_paternal_mass_variance += itr->second->abundances()[i]->paternal_mass_variance() / (double)paternal_status_table[NUMERIC_OK];
				avg_paternal_effective_length += itr->second->abundances()[i]->paternal_effective_length() / (double)paternal_status_table[NUMERIC_OK];
            }
			*/
			avg_paternal_fpkm += itr->second->abundances()[i]->paternal_FPKM() / (double)ab_group_per_replicate.size();
            avg_paternal_num_frags += itr->second->abundances()[i]->num_paternal_fragments() / (double)ab_group_per_replicate.size();
            avg_paternal_gamma += itr->second->abundances()[i]->paternal_gamma() / (double)ab_group_per_replicate.size();
            avg_paternal_mass_variance += itr->second->abundances()[i]->paternal_mass_variance() / (double)ab_group_per_replicate.size();
            avg_paternal_effective_length += itr->second->abundances()[i]->paternal_effective_length() / (double)ab_group_per_replicate.size();
			/*
			if (itr->second->abundances()[i]->maternal_status() == NUMERIC_OK)
            {
                avg_maternal_fpkm += itr->second->abundances()[i]->maternal_FPKM() / (double)maternal_status_table[NUMERIC_OK];
                avg_maternal_num_frags += itr->second->abundances()[i]->num_maternal_fragments() / (double)maternal_status_table[NUMERIC_OK];
                avg_maternal_gamma += itr->second->abundances()[i]->maternal_gamma() / (double)maternal_status_table[NUMERIC_OK];
                avg_maternal_mass_variance += itr->second->abundances()[i]->maternal_mass_variance() / (double)maternal_status_table[NUMERIC_OK];
				avg_maternal_effective_length += itr->second->abundances()[i]->maternal_effective_length() / (double)maternal_status_table[NUMERIC_OK];
            }
			*/
			avg_maternal_fpkm += itr->second->abundances()[i]->maternal_FPKM() / (double)ab_group_per_replicate.size();
            avg_maternal_num_frags += itr->second->abundances()[i]->num_maternal_fragments() / (double)ab_group_per_replicate.size();
            avg_maternal_gamma += itr->second->abundances()[i]->maternal_gamma() / (double)ab_group_per_replicate.size();
            avg_maternal_mass_variance += itr->second->abundances()[i]->maternal_mass_variance() / (double)ab_group_per_replicate.size();
            avg_maternal_effective_length += itr->second->abundances()[i]->maternal_effective_length() / (double)ab_group_per_replicate.size();
        }
		_abundances[i]->paternal_FPKM(avg_paternal_fpkm);
        _abundances[i]->paternal_gamma(avg_paternal_gamma);
        _abundances[i]->num_paternal_fragments(avg_paternal_num_frags);
        _abundances[i]->paternal_mass_variance(avg_paternal_mass_variance);
		_abundances[i]->paternal_effective_length(avg_paternal_effective_length);
		_abundances[i]->maternal_FPKM(avg_maternal_fpkm);
        _abundances[i]->maternal_gamma(avg_maternal_gamma);
        _abundances[i]->num_maternal_fragments(avg_maternal_num_frags);
        _abundances[i]->maternal_mass_variance(avg_maternal_mass_variance);
		_abundances[i]->maternal_effective_length(avg_maternal_effective_length);
        if (paternal_status_table[NUMERIC_OK] >= 1)
        {
            _abundances[i]->paternal_status(NUMERIC_OK);
        }
        else
        {
            if (paternal_status_table[NUMERIC_LOW_DATA] >= paternal_status_table[NUMERIC_FAIL])
            {
                _abundances[i]->paternal_status(NUMERIC_LOW_DATA);
            }
            else if (paternal_status_table[NUMERIC_HI_DATA] >= paternal_status_table[NUMERIC_FAIL])
            {
                _abundances[i]->paternal_status(NUMERIC_HI_DATA);
            }
//            else if (paternal_status_table[NUMERIC_HI_DATA] >= paternal_status_table[NUMERIC_LOW_DATA]) // not sure this ever happens in practice
//            {
//                _abundances[i]->paternal_status(NUMERIC_FAIL);
//            }
            else
            {
                _abundances[i]->paternal_status(NUMERIC_FAIL);
            }
        }
		if (maternal_status_table[NUMERIC_OK] >= 1)
        {
            _abundances[i]->maternal_status(NUMERIC_OK);
        }
        else
        {
            if (maternal_status_table[NUMERIC_LOW_DATA] >= maternal_status_table[NUMERIC_FAIL])
            {
                _abundances[i]->maternal_status(NUMERIC_LOW_DATA);
            }
            if (maternal_status_table[NUMERIC_HI_DATA] >= maternal_status_table[NUMERIC_FAIL])
            {
                _abundances[i]->maternal_status(NUMERIC_HI_DATA);
            }
            if (maternal_status_table[NUMERIC_HI_DATA] >= maternal_status_table[NUMERIC_LOW_DATA]) // not sure this ever happens in practice
            {
                _abundances[i]->maternal_status(NUMERIC_FAIL);
            }
            else
            {
                _abundances[i]->maternal_status(NUMERIC_FAIL);
            }
        }
		_abundances[i]->num_paternal_fragments_by_replicate(paternal_cpr);
		_abundances[i]->num_maternal_fragments_by_replicate(maternal_cpr);
        _abundances[i]->paternal_FPKM_by_replicate(paternal_fpr);
		_abundances[i]->maternal_FPKM_by_replicate(maternal_fpr);
        _abundances[i]->paternal_status_by_replicate(paternal_spr);
		_abundances[i]->maternal_status_by_replicate(maternal_spr);
	}
}

void AlleleAbundanceGroup::calculate_abundance(const vector<MateHit>& alignments, bool perform_collapse)
{
	vector<shared_ptr<Abundance> > transcripts;
	get_transfrags(transcripts);
    map<shared_ptr<ReadGroupProperties const >, vector<MateHit> > alignments_per_read_group;
    for(std::set<shared_ptr<ReadGroupProperties const > >::iterator itr = _read_group_props.begin();
        itr != _read_group_props.end();
        ++itr)
    {
        alignments_per_read_group[*itr] = vector<MateHit>();
        
        vector<MateHit>& rep_hits = alignments_per_read_group[*itr];
        
        for (size_t i = 0; i < alignments.size(); ++i)
        {
            if (alignments[i].read_group_props() == *itr)
            {
                rep_hits.push_back(alignments[i]);
            }
        }
		
	}
	
	collect_per_replicate_mass(alignments, transcripts);
	std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AlleleAbundanceGroup> > ab_group_per_replicate;
	calculate_per_replicate_abundances(transcripts,
                                       alignments_per_read_group,
                                       ab_group_per_replicate);
	std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<const AlleleAbundanceGroup> > const_ab_group_per_replicate;
	for (std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AlleleAbundanceGroup> >::iterator itr = ab_group_per_replicate.begin();
         itr != ab_group_per_replicate.end(); ++itr)
    {
        const_ab_group_per_replicate[itr->first] = itr->second;
    }
    
    aggregate_replicate_abundances(const_ab_group_per_replicate);

    calculate_abundance_group_variance(transcripts, const_ab_group_per_replicate);
}

void AlleleAbundanceGroup::update_multi_reads(const vector<MateHit>& alignments, vector<shared_ptr<Abundance> > transcripts)
{
	size_t M = alignments.size();
	size_t N = transcripts.size();
	
	if (transcripts.empty())
		return;
    
    for (size_t i = 0; i < M; ++i)
	{
		if (alignments[i].is_multi())
		{
			double paternal_expr = 0.0;
			double maternal_expr = 0.0;
			for (size_t j = 0; j < N; ++j)
			{
				paternal_expr += _abundances[j]->paternal_cond_probs()->at(i) * _abundances[j]->paternal_FPKM() * _abundances[j]->paternal_effective_length();
				maternal_expr += _abundances[j]->maternal_cond_probs()->at(i) * _abundances[j]->maternal_FPKM() * _abundances[j]->maternal_effective_length();
			}
			alignments[i].read_group_props()->multi_read_table()->add_expr(alignments[i], paternal_expr+maternal_expr, paternal_expr, maternal_expr);
		}
	}
}

bool simulate_count_covariance(const vector<double>& num_paternal_fragments,
							   const vector<double>& num_maternal_fragments,
                               const vector<double>& paternal_frag_variances,
							   const vector<double>& maternal_frag_variances,
                               const ublas::matrix<double>& iterated_exp_parental_count_covariance,
                               const vector<shared_ptr<Abundance> >& transcripts,
                               ublas::matrix<double>& parental_count_covariance,
                               vector<Eigen::VectorXd>& paternal_assigned_count_samples,
							   vector<Eigen::VectorXd>& maternal_assigned_count_samples,
                               vector<ublas::vector<double> >* paternal_gamma_samples = NULL,
							   vector<ublas::vector<double> >* maternal_gamma_samples = NULL)
{
    parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 2*transcripts.size());
    
    
    double total_paternal_frag_counts = accumulate(num_paternal_fragments.begin(), num_paternal_fragments.end(), 0.0);
	double total_maternal_frag_counts = accumulate(num_maternal_fragments.begin(), num_maternal_fragments.end(), 0.0);
    double total_paternal_var = accumulate(paternal_frag_variances.begin(), paternal_frag_variances.end(), 0.0);
	double total_maternal_var = accumulate(maternal_frag_variances.begin(), maternal_frag_variances.end(), 0.0);
    
    if (total_paternal_frag_counts == 0 && total_maternal_frag_counts == 0)
    {
        parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 2*transcripts.size());
        
        size_t num_paternal_frags = num_frag_assignments;
		size_t num_maternal_frags = num_frag_assignments;
        if (paternal_gamma_samples && paternal_gamma_samples->empty() == false)
            num_paternal_frags = paternal_gamma_samples->size();
		if (maternal_gamma_samples && maternal_gamma_samples->empty() == false)
            num_maternal_frags = maternal_gamma_samples->size();
        
        paternal_assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * num_paternal_frags, Eigen::VectorXd::Zero(transcripts.size()));
		maternal_assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * num_maternal_frags, Eigen::VectorXd::Zero(transcripts.size()));
        return true;
    }
    
    vector<ublas::vector<double> > paternal_assigned_gamma_samples,maternal_assigned_gamma_samples;
    
    boost::mt19937 rng;
    
    if (paternal_gamma_samples == NULL || paternal_gamma_samples->empty() || maternal_gamma_samples == NULL || maternal_gamma_samples->empty()) //???not sure about the || term between alleles
    {
        vector<Eigen::VectorXd > paternal_generated_counts (num_frag_count_draws, Eigen::VectorXd::Zero(transcripts.size()));
		vector<Eigen::VectorXd > maternal_generated_counts (num_frag_count_draws, Eigen::VectorXd::Zero(transcripts.size()));
        ublas::vector<double> paternal_mle_frag_counts = ublas::zero_vector<double>(transcripts.size());
		ublas::vector<double> maternal_mle_frag_counts = ublas::zero_vector<double>(transcripts.size());
        
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            paternal_mle_frag_counts(j) = num_paternal_fragments[j];
			maternal_mle_frag_counts(j) = num_maternal_fragments[j];
        }
        
        //    cerr << "********************" << endl;
        //    cerr << "initial paternal MLE counts: " << paternal_mle_frag_counts << endl;
        
        ublas::matrix<double> parental_mle_count_covar = iterated_exp_parental_count_covariance;
        
        ublas::matrix<double> parental_epsilon = ublas::zero_matrix<double>(2*transcripts.size(),2*transcripts.size());
        for (size_t i = 0; i < 2*transcripts.size(); ++i)
        {
            parental_epsilon(i,i) = 1e-6;
		}
        
        parental_mle_count_covar += parental_epsilon; // modify matrix to avoid problems during inverse
        
        double ret = cholesky_factorize(parental_mle_count_covar);
        if (ret != 0)
        {
            fprintf(stderr, "Warning: Iterated expectation parental count covariance matrix cannot be cholesky factorized!\n");
            //fprintf(stderr, "Warning: FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            //        for (size_t j = 0; j < _abundances.size(); ++j)
            //        {
            //            _abundances[j]->status(NUMERIC_FAIL);
            //        }
            return false;
        }
        ublas::vector<double> parental_mle_frag_counts(paternal_mle_frag_counts.size()+maternal_mle_frag_counts.size());
		for(size_t c = 0;c < paternal_mle_frag_counts.size();++c)
			parental_mle_frag_counts[c] = paternal_mle_frag_counts[c];
		for(size_t c = 0;c < maternal_mle_frag_counts.size();++c)
			parental_mle_frag_counts[c+maternal_mle_frag_counts.size()] = maternal_mle_frag_counts[c];
		
        multinormal_generator<double> generator(parental_mle_frag_counts, parental_mle_count_covar);
        
        for (size_t assign_idx = 0; assign_idx < num_frag_assignments; ++assign_idx)
        {
            boost::numeric::ublas::vector<double> random_count_assign,paternal_random_count_assign,maternal_random_count_assign;
            double total_paternal_sample_counts = 0;
			double total_maternal_sample_counts = 0;
            do {
                random_count_assign = generator.next_rand();
                //cerr << random_count_assign << endl;
                
                for (size_t r_idx = 0; r_idx < random_count_assign.size(); ++r_idx)
                {
					if(r_idx < random_count_assign.size()/2){
						paternal_random_count_assign(r_idx) = random_count_assign(r_idx);
						if(paternal_random_count_assign(r_idx) < 0)
							paternal_random_count_assign(r_idx) = 0;
					}
					else{
						maternal_random_count_assign(r_idx-random_count_assign.size()/2) = random_count_assign(r_idx);
						if(maternal_random_count_assign(r_idx-random_count_assign.size()/2) < 0)
							maternal_random_count_assign(r_idx-random_count_assign.size()/2) = 0;
					}						
                }
                
                total_paternal_sample_counts = accumulate(paternal_random_count_assign.begin(), paternal_random_count_assign.end(), 0.0);
				total_maternal_sample_counts = accumulate(maternal_random_count_assign.begin(), maternal_random_count_assign.end(), 0.0);
                if (total_paternal_sample_counts > 0)
                    paternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (paternal_random_count_assign / (total_paternal_sample_counts+total_maternal_sample_counts));
                else
                    paternal_random_count_assign = boost::numeric::ublas::zero_vector<double>(transcripts.size());
				if (total_maternal_sample_counts > 0)
                    maternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (maternal_random_count_assign / (total_paternal_sample_counts+total_maternal_sample_counts));
                else
                    maternal_random_count_assign = boost::numeric::ublas::zero_vector<double>(transcripts.size());
            } while(total_paternal_sample_counts+total_maternal_sample_counts <= 0);
            
            paternal_assigned_gamma_samples.push_back(paternal_random_count_assign);
			maternal_assigned_gamma_samples.push_back(maternal_random_count_assign);
        }
    }
    else
    {
        for (size_t assign_idx = 0; assign_idx < paternal_gamma_samples->size(); ++assign_idx)
        {
            boost::numeric::ublas::vector<double> paternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (*paternal_gamma_samples)[assign_idx];
            paternal_assigned_gamma_samples.push_back(paternal_random_count_assign);
        }
		for (size_t assign_idx = 0; assign_idx < maternal_gamma_samples->size(); ++assign_idx)
        {
            boost::numeric::ublas::vector<double> maternal_random_count_assign = (total_paternal_frag_counts+total_maternal_frag_counts) * (*maternal_gamma_samples)[assign_idx];
            maternal_assigned_gamma_samples.push_back(maternal_random_count_assign);
        }
    }
    
    paternal_assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * paternal_assigned_gamma_samples.size(), Eigen::VectorXd::Zero(transcripts.size()));
	maternal_assigned_count_samples = vector<Eigen::VectorXd> (num_frag_count_draws * maternal_assigned_gamma_samples.size(), Eigen::VectorXd::Zero(transcripts.size()));
    
    Eigen::VectorXd paternal_expected_generated_counts = Eigen::VectorXd::Zero(transcripts.size());
	Eigen::VectorXd maternal_expected_generated_counts = Eigen::VectorXd::Zero(transcripts.size());
    
    
    boost::uniform_01<> uniform_dist;
    boost::mt19937 null_rng;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uniform_gen(null_rng, uniform_dist);
    
    assert(paternal_assigned_gamma_samples.size() == maternal_assigned_gamma_samples.size());
	for (size_t assign_idx = 0; assign_idx < paternal_assigned_gamma_samples.size(); ++assign_idx)
    {
        boost::numeric::ublas::vector<double>& paternal_random_count_assign = paternal_assigned_gamma_samples[assign_idx];
		boost::numeric::ublas::vector<double>& maternal_random_count_assign = maternal_assigned_gamma_samples[assign_idx];
        
        //cerr << random_count_assign << endl;
        for (size_t gen_idx = 0; gen_idx < num_frag_count_draws; ++gen_idx)
        {
            Eigen::VectorXd paternal_generated_and_assigned_counts = Eigen::VectorXd::Zero(transcripts.size());
			Eigen::VectorXd maternal_generated_and_assigned_counts = Eigen::VectorXd::Zero(transcripts.size());
            
            for (size_t j = 0; j < transcripts.size(); ++j)
            {
                double pr =  paternal_random_count_assign(j);
				double mr =  maternal_random_count_assign(j);
                
                //pr = pr < 1 ? 1 : pr;
				//mr = mr < 1 ? 1 : mr;
                
                if (pr > 0)
                {
                    //double paternal_fit_var = _abundances[j]->paternal_mass_variance();
                    
                    double paternal_fit_var = total_paternal_var * (paternal_random_count_assign(j) / (total_paternal_frag_counts+total_maternal_frag_counts));
                    double paternal_frags = paternal_random_count_assign(j);
                    
                    if (paternal_fit_var - paternal_frags > 1e-1)
                    {
                        pr *= pr;
                        double paternal_over_disp_scale = paternal_fit_var - paternal_frags;
                        pr /= paternal_over_disp_scale;

                        double paternal_after_decimal = pr - (long)pr;
                        //fprintf( stderr, "paternal after decimal = %lg\n", paternal_after_decimal);
                        if (uniform_gen() < paternal_after_decimal)
                            pr = floor(pr);
                        else
                            pr = ceil(pr);
                        
                        if (pr == 0)
                        {
                            paternal_generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        //double p = _abundances[j]->num_paternal_fragments() / paternal_fit_var;
                        double p = pr / (pr + paternal_frags);
                        
                        negative_binomial_distribution<int, double> nb_j(pr, p);
                        
                        paternal_generated_and_assigned_counts(j) = nb_j(rng);
                        
                    }
                    else
                    {
                        double paternal_after_decimal = pr - (long)pr;
                        //fprintf( stderr, "paternal after decimal = %lg\n", paternal_after_decimal);
                        if (uniform_gen() < paternal_after_decimal)
                            pr = floor(pr);
                        else
                            pr = ceil(pr);
                        
                        if (pr == 0)
                        {
                            paternal_generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        boost::random::poisson_distribution<int, double> nb_j(pr);
                        
                        paternal_generated_and_assigned_counts(j) = nb_j(rng);
                        
                    }
                }
                else
                {
                    paternal_generated_and_assigned_counts(j) = 0;
                }
				//
				if (mr > 0)
                {
                    //double maternal_fit_var = _abundances[j]->maternal_mass_variance();
                    
                    double maternal_fit_var = total_maternal_var * (maternal_random_count_assign(j) / (total_paternal_frag_counts+total_maternal_frag_counts));
                    double maternal_frags = maternal_random_count_assign(j);
                    
                    if (maternal_fit_var - maternal_frags > 1e-1)
                    {
                        mr *= mr;
                        double maternal_over_disp_scale = maternal_fit_var - maternal_frags;
                        mr /= maternal_over_disp_scale;
                        
						double maternal_after_decimal = mr - (long)mr;
                        //fprintf( stderr, "maternal after decimal = %lg\n", maternal_after_decimal);
                        if (uniform_gen() < maternal_after_decimal)
                            mr = floor(mr);
                        else
                            mr = ceil(mr);
                        
                        if (mr == 0)
                        {
                            maternal_generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        //double p = _abundances[j]->num_maternal_fragments() / maternal_fit_var;
                        double p = mr / (mr + maternal_frags);
                        
                        negative_binomial_distribution<int, double> nb_j(mr, p);
                        
                        maternal_generated_and_assigned_counts(j) = nb_j(rng);
                        
                    }
                    else
                    {
                        double maternal_after_decimal = mr - (long)mr;
                        //fprintf( stderr, "maternal after decimal = %lg\n", maternal_after_decimal);
                        if (uniform_gen() < maternal_after_decimal)
                            mr = floor(mr);
                        else
                            mr = ceil(mr);
                        
                        if (mr == 0)
                        {
                            maternal_generated_and_assigned_counts(j) = 0;
                            continue;
                        }
                        
                        boost::random::poisson_distribution<int, double> nb_j(mr);
                        
                        maternal_generated_and_assigned_counts(j) = nb_j(rng);
                        
                    }
                }
                else
                {
                    maternal_generated_and_assigned_counts(j) = 0;
                }
            }
			//cerr << "     paternal assigned count sample: " << paternal_generated_and_assigned_counts.transpose() << endl;
            paternal_assigned_count_samples[assign_idx*num_frag_count_draws + gen_idx] = paternal_generated_and_assigned_counts;
            //cerr << "     maternal assigned count sample: " << maternal_generated_and_assigned_counts.transpose() << endl;
            maternal_assigned_count_samples[assign_idx*num_frag_count_draws + gen_idx] = maternal_generated_and_assigned_counts;
        }
    }
    
    Eigen::VectorXd paternal_expected_counts = Eigen::VectorXd::Zero(transcripts.size());
	Eigen::VectorXd maternal_expected_counts = Eigen::VectorXd::Zero(transcripts.size());
    Eigen::VectorXd paternal_expected_relative_abundances = Eigen::VectorXd::Zero(transcripts.size());
	Eigen::VectorXd maternal_expected_relative_abundances = Eigen::VectorXd::Zero(transcripts.size());
    
    for (size_t i = 0; i < paternal_assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < paternal_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(paternal_assigned_count_samples[i](j)) && !isinf(paternal_assigned_count_samples[i](j)));
        }
        
        
        paternal_expected_counts += paternal_assigned_count_samples[i];
        //
        //paternal_expected_relative_abundances += paternal_relative_abundances[i];
    }
	for (size_t i = 0; i < maternal_assigned_count_samples.size(); ++i)
    {
        for (int j = 0; j < maternal_assigned_count_samples[i].size(); ++j)
        {
            assert (!isnan(maternal_assigned_count_samples[i](j)) && !isinf(maternal_assigned_count_samples[i](j)));
        }
        
        
        maternal_expected_counts += maternal_assigned_count_samples[i];
        //
        //maternal_expected_relative_abundances += maternal_relative_abundances[i];
    }
	assert(paternal_assigned_count_samples.size() == maternal_assigned_count_samples.size());
    if (paternal_assigned_count_samples.size() > 0)
    {
        paternal_expected_counts /= paternal_assigned_count_samples.size();
        //paternal_expected_generated_counts /= paternal_assigned_counts.size();
        //paternal_expected_relative_abundances /= paternal_assigned_counts.size();
    }
	if (maternal_assigned_count_samples.size() > 0)
    {
        maternal_expected_counts /= maternal_assigned_count_samples.size();
        //maternal_expected_generated_counts /= maternal_assigned_counts.size();
        //maternal_expected_relative_abundances /= maternal_assigned_counts.size();
    }
    
    //    if (num_frag_assignments > 0)
    //    {
    //        paternal_expected_generated_counts /= num_frag_assignments;
	//        maternal_expected_generated_counts /= num_frag_assignments;
    //        //paternal_expected_relative_abundances /= paternal_assigned_counts.size();
	//        //maternal_expected_relative_abundances /= maternal_assigned_counts.size();
    //    }
    
    
    //    cerr << "======" << endl;
    //    cerr << "updated paternal expected counts #1: " << endl;
    //    std::cerr << paternal_expected_counts << std::endl;
    //    cerr << "updated paternal expected generated counts #1: " << endl;
    //    std::cerr << paternal_expected_generated_counts << std::endl;
    //    cerr << "updated maternal expected counts #1: " << endl;
    //    std::cerr << maternal_expected_counts << std::endl;
    //    cerr << "updated maternal expected generated counts #1: " << endl;
    //    std::cerr << maternal_expected_generated_counts << std::endl;
    //    
    //    cerr << "======" << endl;
    
	
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            for (size_t k = 0 ; k < paternal_assigned_count_samples.size(); ++k)
            {
                double cpp = (paternal_assigned_count_samples[k](i) - paternal_expected_counts(i)) * (paternal_assigned_count_samples[k](j) - paternal_expected_counts(j));
                parental_count_covariance(i,j) += cpp;
				assert (!isinf(parental_count_covariance(i,j)) && !isnan(parental_count_covariance(i,j)));
                //double r = (paternal_relative_abundances[k](i) - paternal_expected_relative_abundances(i)) * (paternal_relative_abundances[k](j) - paternal_expected_relative_abundances(j));
				//_parental_kappa_covariance(i,j) +=
				
				double cpm = (paternal_assigned_count_samples[k](i) - paternal_expected_counts(i)) * (maternal_assigned_count_samples[k](j) - maternal_expected_counts(j));
                parental_count_covariance(i,j+transcripts.size()) += cpm;
				assert (!isinf(parental_count_covariance(i,j+transcripts.size())) && !isnan(parental_count_covariance(i,j+transcripts.size())));
                //double r = (paternal_relative_abundances[k](i) - paternal_expected_relative_abundances(i)) * (maternal_relative_abundances[k](j) - maternal_expected_relative_abundances(j));
                //_parental_kappa_covariance(i,j+transcripts.size()) +=
				
				double cmm = (maternal_assigned_count_samples[k](i) - maternal_expected_counts(i)) * (maternal_assigned_count_samples[k](j) - maternal_expected_counts(j));
                parental_count_covariance(i+transcripts.size(),j+transcripts.size()) += cmm;
				assert (!isinf(parental_count_covariance(i+transcripts.size(),j+transcripts.size())) && !isnan(parental_count_covariance(i+transcripts.size(),j+transcripts.size())));
                //double r = (maternal_relative_abundances[k](i) - maternal_expected_relative_abundances(i)) * (maternal_relative_abundances[k](j) - maternal_expected_relative_abundances(j));
				//_parental_kappa_covariance(i+transcripts.size(),j+transcripts.size()) +=
				
				double cmp = (maternal_assigned_count_samples[k](i) - maternal_expected_counts(i)) * (paternal_assigned_count_samples[k](j) - paternal_expected_counts(j));
                parental_count_covariance(i+transcripts.size(),j) += cmp;
				assert (!isinf(parental_count_covariance(i+transcripts.size(),j)) && !isnan(parental_count_covariance(i+transcripts.size(),j)));
                //double r = (maternal_relative_abundances[k](i) - maternal_expected_relative_abundances(i)) * (paternal_relative_abundances[k](j) - paternal_expected_relative_abundances(j));
				//_parental_kappa_covariance(i+transcripts.size(),j) +=
            }
        }
    }
    
    if (paternal_assigned_count_samples.empty() == false && maternal_assigned_count_samples.empty() == false)
        parental_count_covariance /= paternal_assigned_count_samples.size();
    
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        // Make sure we aren't below the fit for the single isoform case
        if (parental_count_covariance(i,i) < ceil(paternal_frag_variances[i]))
        {
            fprintf(stderr, "Counts for paternal %d (var = %lg) are underdispersed, reverting to fitted variance model (%lg)\n", i, parental_count_covariance(i,i), ceil(paternal_frag_variances[i]));
            parental_count_covariance(i,i) = ceil(paternal_frag_variances[i]);
            assert (!isinf(parental_count_covariance(i,i)) && !isnan(parental_count_covariance(i,i)));
        }
        if (parental_count_covariance(i+transcripts.size(),i+transcripts.size()) < ceil(maternal_frag_variances[i]))
        {
            fprintf(stderr, "Counts for maternal %d (var = %lg) are underdispersed, reverting to fitted variance model (%lg)\n", i, parental_count_covariance(i+transcripts.size(),i+transcripts.size()), ceil(maternal_frag_variances[i]));
            parental_count_covariance(i+transcripts.size(),i+transcripts.size()) = ceil(maternal_frag_variances[i]);
            assert (!isinf(parental_count_covariance(i+transcripts.size(),i+transcripts.size())) && !isnan(parental_count_covariance(i+transcripts.size(),i+transcripts.size())));
        }
		
        // Check that we aren't below what the Poisson model says we ought to be at
        if (parental_count_covariance(i,i) < ceil(num_paternal_fragments[i] + iterated_exp_parental_count_covariance(i,i)))
        {
            fprintf(stderr, "Counts for paternal %d (var = %lg) are underdispersed, reverting to additive variance model (%lg)\n", i, parental_count_covariance(i,i),  ceil(num_paternal_fragments[i] + iterated_exp_parental_count_covariance(i,i)));
            parental_count_covariance(i,i) = ceil(num_paternal_fragments[i] + iterated_exp_parental_count_covariance(i,i));
            assert (!isinf(parental_count_covariance(i,i)) && !isnan(parental_count_covariance(i,i)));
        }
		if (parental_count_covariance(i+transcripts.size(),i+transcripts.size()) < ceil(num_maternal_fragments[i] + iterated_exp_parental_count_covariance(i+transcripts.size(),i+transcripts.size())))
        {
            fprintf(stderr, "Counts for maternal %d (var = %lg) are underdispersed, reverting to additive variance model (%lg)\n", i, parental_count_covariance(i+transcripts.size(),i+transcripts.size()),  ceil(num_maternal_fragments[i] + iterated_exp_parental_count_covariance(i+transcripts.size(),i+transcripts.size())));
            parental_count_covariance(i+transcripts.size(),i+transcripts.size()) = ceil(num_maternal_fragments[i] + iterated_exp_parental_count_covariance(i+transcripts.size(),i+transcripts.size()));
            assert (!isinf(parental_count_covariance(i+transcripts.size(),i+transcripts.size())) && !isnan(parental_count_covariance(i+transcripts.size(),i+transcripts.size())));
        }
		
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            assert (!isinf(parental_count_covariance(i,j)) && !isnan(parental_count_covariance(i,j)));
			assert (!isinf(parental_count_covariance(i,j++transcripts.size())) && !isnan(parental_count_covariance(i,j++transcripts.size())));
			assert (!isinf(parental_count_covariance(i++transcripts.size(),j++transcripts.size())) && !isnan(parental_count_covariance(i++transcripts.size(),j++transcripts.size())));
			assert (!isinf(parental_count_covariance(i++transcripts.size(),j)) && !isnan(parental_count_covariance(i++transcripts.size(),j)));
        }
    }
    
    return true;
}

void AlleleAbundanceGroup::generate_fpkm_samples()
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
    vector<vector<double> > paternal_fpkm_sample_vectors(_abundances.size());
	vector<vector<double> > maternal_fpkm_sample_vectors(_abundances.size());
    vector<double> paternal_group_sum_fpkm_samples,maternal_group_sum_fpkm_samples;
    
    vector<double> paternal_fpkm_means(_abundances.size(), 0);
	vector<double> maternal_fpkm_means(_abundances.size(), 0);
    assert(_paternal_assigned_count_samples.size() == _maternal_assigned_count_samples.size());
	
    for(size_t i = 0; i < _paternal_assigned_count_samples.size(); ++i)
    {
        const Eigen::VectorXd paternal_sample = _paternal_assigned_count_samples[i];
        double total_paternal_fpkm = 0;
		const Eigen::VectorXd maternal_sample = _maternal_assigned_count_samples[i];
        double total_maternal_fpkm = 0;
        
        for (size_t j = 0; j < paternal_sample.size(); ++j)
        {
            double paternal_fpkm_sample = paternal_sample[j] / M;
            
            if (_abundances[j]->paternal_effective_length() > 0)
            {
                paternal_fpkm_sample *= 1000000000;
                paternal_fpkm_sample /= _abundances[j]->paternal_effective_length();
                paternal_fpkm_sample /= external_scale_factor;
                //double paternal_standard_fpkm = _abundances[j]->paternal_FPKM();
                //fprintf(stderr, "count = %lg, paternal fpkm = %lg, paternal standard fpkm = %lg\n", paternal_sample[j], paternal_fpkm_sample, paternal_standard_fpkm);
            }
            else
            {
                paternal_fpkm_sample = 0;
            }
            
            assert (isnan(paternal_fpkm_sample) == false);
            
            paternal_fpkm_sample_vectors[j].push_back(paternal_fpkm_sample);
            paternal_fpkm_means[j] += paternal_fpkm_sample;
            total_paternal_fpkm += paternal_fpkm_sample;
        }
		if (paternal_effective_length() > 0)
        {
			paternal_group_sum_fpkm_samples.push_back(total_paternal_fpkm);
        }
        else
        {
            paternal_group_sum_fpkm_samples.push_back(0);
        }
		for (size_t j = 0; j < maternal_sample.size(); ++j)
        {
            double maternal_fpkm_sample = maternal_sample[j] / M;
            
            if (_abundances[j]->maternal_effective_length() > 0)
            {
                maternal_fpkm_sample *= 1000000000;
                maternal_fpkm_sample /= _abundances[j]->maternal_effective_length();
                maternal_fpkm_sample /= external_scale_factor;
                //double maternal_standard_fpkm = _abundances[j]->maternal_FPKM();
                //fprintf(stderr, "count = %lg, maternal fpkm = %lg, maternal standard fpkm = %lg\n", maternal_sample[j], maternal_fpkm_sample, maternal_standard_fpkm);
            }
            else
            {
                maternal_fpkm_sample = 0;
            }
            
            assert (isnan(maternal_fpkm_sample) == false);
            
            maternal_fpkm_sample_vectors[j].push_back(maternal_fpkm_sample);
            maternal_fpkm_means[j] += maternal_fpkm_sample;
            total_maternal_fpkm += maternal_fpkm_sample;
        }
		if (maternal_effective_length() > 0)
        {
			maternal_group_sum_fpkm_samples.push_back(total_maternal_fpkm);
        }
        else
        {
            maternal_group_sum_fpkm_samples.push_back(0);
        }
    }
    	
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        paternal_fpkm_means[i] /= _paternal_assigned_count_samples.size();
        _abundances[i]->paternal_fpkm_samples(paternal_fpkm_sample_vectors[i]);
        //fprintf(stderr, "standard paternal fpkm = %lg, paternal sample mean = %lg\n", _abundances[i]->paternal_FPKM(), paternal_fpkm_means[i]);
		maternal_fpkm_means[i] /= _maternal_assigned_count_samples.size();
        _abundances[i]->maternal_fpkm_samples(maternal_fpkm_sample_vectors[i]);
        //fprintf(stderr, "standard maternal fpkm = %lg, maternal sample mean = %lg\n", _abundances[i]->maternal_FPKM(), maternal_fpkm_means[i]);
    }
    
    vector<Eigen::VectorXd> paternal_mem_fpkm_samples,maternal_mem_fpkm_samples;
	assert(paternal_fpkm_sample_vectors.size() == maternal_fpkm_sample_vectors.size());
	
    for (size_t i = 0; i < paternal_fpkm_sample_vectors.size(); ++i)
    {
        Eigen::VectorXd paternal_sample(paternal_fpkm_sample_vectors[i].size());
		Eigen::VectorXd maternal_sample(maternal_fpkm_sample_vectors[i].size());
        paternal_mem_fpkm_samples.push_back(paternal_sample);
		maternal_mem_fpkm_samples.push_back(maternal_sample);
        for (size_t j = 0; j < paternal_fpkm_sample_vectors[i].size(); ++j)
        {
            paternal_mem_fpkm_samples[i](j) = paternal_fpkm_sample_vectors[i][j];
        }
		for (size_t j = 0; j < maternal_fpkm_sample_vectors[i].size(); ++j)
        {
            maternal_mem_fpkm_samples[i](j) = maternal_fpkm_sample_vectors[i][j];
        }
    }
	
    paternal_member_fpkm_samples(paternal_mem_fpkm_samples);
	maternal_member_fpkm_samples(maternal_mem_fpkm_samples);

    paternal_fpkm_samples(paternal_group_sum_fpkm_samples);
    maternal_fpkm_samples(maternal_group_sum_fpkm_samples);
    assert (paternal_group_sum_fpkm_samples.empty() == false);
	assert (maternal_group_sum_fpkm_samples.empty() == false);
}

void AlleleAbundanceGroup::calculate_FPKM_covariance()
{
	
	if (paternal_effective_length() == 0 || maternal_effective_length() == 0)
	{
		_parental_fpkm_covariance = ublas::zero_matrix<double>(2*_abundances.size(), 2*_abundances.size());
		return;
	}
    
    //estimate_count_covariance();
    
    long double total_paternal_var = 0.0;
	long double total_maternal_var = 0.0;
    
    double paternal_abundance_weighted_length = 0.0;
	double maternal_abundance_weighted_length = 0.0;
    double total_paternal_abundance = 0.0;
	double total_maternal_abundance = 0.0;
    
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
    
    for (size_t j = 0; j < 2*_abundances.size(); ++j)
    {
		if(j < _abundances.size())
		{
			paternal_abundance_weighted_length += _abundances[j]->paternal_effective_length() * _abundances[j]->paternal_FPKM();
			total_paternal_abundance += _abundances[j]->paternal_FPKM();
		}
		else{
			maternal_abundance_weighted_length += _abundances[j-_abundances.size()]->maternal_effective_length() * _abundances[j-_abundances.size()]->maternal_FPKM();
			total_maternal_abundance += _abundances[j-_abundances.size()]->maternal_FPKM();
		}
		
		for (size_t i = 0; i < 2*_abundances.size(); ++i)
        {
            _parental_fpkm_covariance(i,j) = _parental_count_covariance(i,j);
            
            // FPKMs need to be on the external scale, so we can compare them
            // between conditions.  Counts are internally scaled up until here.
            _parental_fpkm_covariance(i,j) *= 1.0 / (external_scale_factor * external_scale_factor);
            assert (!isinf(_parental_count_covariance(i,j)) && !isnan(_parental_fpkm_covariance(i,j)));
            long double length_i,length_j;
			if(j < _abundances.size())
			{
				length_j = _abundances[j]->paternal_effective_length();
			}
			else{
				length_j = _abundances[j-_abundances.size()]->maternal_effective_length();
			}
			if(i < _abundances.size())
			{
				length_i = _abundances[i]->paternal_effective_length();
			}
			else
			{
				length_i = _abundances[i-_abundances.size()]->maternal_effective_length();
			}
            
            assert (!isinf(length_i) && !isnan(length_i));
            assert (!isinf(length_j) && !isnan(length_j));
            if (length_i > 0 && length_j > 0 & M > 0)
            {
                _parental_fpkm_covariance(i,j) *=
                    ((1000000000.0 / (length_j *M)))*((1000000000.0 / (length_i *M)));
                assert (!isinf(_parental_fpkm_covariance(i,j)) && !isnan(_parental_fpkm_covariance(i,j)));
            }
            else
            {
                _parental_fpkm_covariance(i,j) = 0.0;
            }
            
            if (i == j)
            {
				if(j < _abundances.size() && i < _abundances.size())
				{
					double fpkm = _abundances[i]->paternal_FPKM();
					double fpkm_var = _parental_fpkm_covariance(i,j);
//                if (fpkm > 0 && fpkm_var == 0 )
//                {
//                    cerr << _parental_count_covariance << endl;
//                }
					assert (fpkm == 0 || fpkm_var > 0 || _abundances[i]->paternal_status() != NUMERIC_OK);
					assert (!isinf(fpkm_var) && !isnan(fpkm_var));
					_abundances[i]->paternal_FPKM_variance(fpkm_var);
					_abundances[i]->num_paternal_fragment_var(_parental_count_covariance(i,j));
					total_paternal_var += _parental_fpkm_covariance(i,j);
				}
				else{
					double fpkm = _abundances[i-_abundances.size()]->maternal_FPKM();
					double fpkm_var = _parental_fpkm_covariance(i,j);
//                if (fpkm > 0 && fpkm_var == 0 )
//                {
//                    cerr << _parental_count_covariance << endl;
//                }
					assert (fpkm == 0 || fpkm_var > 0 || _abundances[i-_abundances.size()]->maternal_status() != NUMERIC_OK);
					assert (!isinf(fpkm_var) && !isnan(fpkm_var));
					_abundances[i-_abundances.size()]->maternal_FPKM_variance(fpkm_var);
					_abundances[i-_abundances.size()]->num_maternal_fragment_var(_parental_count_covariance(i,j));
					total_maternal_var += _parental_fpkm_covariance(i,j);
				}
            }
			assert (!isinf(_parental_fpkm_covariance(i,j)) && !isnan(_parental_fpkm_covariance(i,j)));
        }
    }
    
    _paternal_FPKM_variance = total_paternal_var;
	_maternal_FPKM_variance = total_maternal_var;
    
    if (final_est_run && library_type != "transfrags")
    {
        ublas::matrix<double> test = _parental_fpkm_covariance;
        double ret = cholesky_factorize(test);
        if (ret != 0 || (_paternal_FPKM_variance < 0 && paternal_status() == NUMERIC_OK))
        {
            fprintf(stderr, "Warning: FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            for (size_t j = 0; j < _abundances.size(); ++j)
            {
                _abundances[j]->paternal_status(NUMERIC_FAIL);
            }
        }
		if (ret != 0 || (_maternal_FPKM_variance < 0 && maternal_status() == NUMERIC_OK))
        {
            fprintf(stderr, "Warning: FPKM covariance is not positive definite (ret = %lg)!\n", ret);
            for (size_t j = 0; j < _abundances.size(); ++j)
            {
                _abundances[j]->maternal_status(NUMERIC_FAIL);
            }
        }
		assert(paternal_status() == maternal_status());
		assert (paternal_FPKM()+maternal_FPKM() == 0 || (_paternal_FPKM_variance > 0 || paternal_status() != NUMERIC_OK) || (_maternal_FPKM_variance > 0 || maternal_status() != NUMERIC_OK));
    }
    
//    cerr << "FPKM covariance: " << endl;
//    for (unsigned i = 0; i < _parental_fpkm_covariance.size1 (); ++ i)
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (_parental_fpkm_covariance, i);
//        cerr << i << " : " << _abundances[i]->paternal_FPKM() << " : "<< _abundances[i]->maternal_FPKM() ; ";
//        std::cerr << i << " : " << mr << std::endl;
//    }
    assert (!isinf(_paternal_FPKM_variance) && !isnan(_paternal_FPKM_variance));
	assert (!isinf(_maternal_FPKM_variance) && !isnan(_maternal_FPKM_variance));    
}

void AlleleAbundanceGroup::calculate_conf_intervals()
{
	// We only really ever call this function for primary abundance groups
    // (i.e. the transcript groups and read bundles with which we calculate
    // transcript MLE expression levels.  Genes, TSS groups, etc get broken
    // off of primary bundles, so we should not call this function on those
    // secondary groups.  The group splitting code needs to manage the task
    // of splitting up all the variout covariance matrices we're calculating
    // here.
	if (paternal_status() == NUMERIC_OK)
	{
		// This will compute the paternal transcript level FPKM confidence intervals
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
            long double paternal_fpkm_var = _abundances[j]->paternal_FPKM_variance();
            double paternal_FPKM_hi = 0.0;
            double paternal_FPKM_lo = 0.0;
            if (_abundances[j]->paternal_status() != NUMERIC_FAIL)
            {
                paternal_FPKM_hi = _abundances[j]->paternal_FPKM() + 2 * sqrt(paternal_fpkm_var);
                paternal_FPKM_lo = max(0.0, (double)(_abundances[j]->paternal_FPKM() - 2 * sqrt(paternal_fpkm_var)));
                if (!(paternal_FPKM_lo <= _abundances[j]->paternal_FPKM() && _abundances[j]->paternal_FPKM() <= paternal_FPKM_hi))
                {
                    //fprintf(stderr, "Error: paternal confidence intervals are illegal! var = %Lg, fpkm = %lg, lo = %lg, hi %lg, status = %d\n", paternal_fpkm_var, _abundances[j]->paternal_FPKM(), paternal_FPKM_lo, paternal_FPKM_hi, _abundances[j]->paternal_status());
                }
                assert (paternal_FPKM_lo <= _abundances[j]->paternal_FPKM() && _abundances[j]->paternal_FPKM() <= paternal_FPKM_hi);
                ConfidenceInterval paternal_conf(paternal_FPKM_lo, paternal_FPKM_hi);
                _abundances[j]->paternal_FPKM_conf(paternal_conf);
                //_abundances[j]->paternal_FPKM_variance(paternal_fpkm_var);
            }
            else
            {
                // we shouldn't be able to get here
                assert(false);
                // TODO: nothing to do here?
            }
		}
		
		// Now build a confidence interval for the whole abundance group
		double paternal_group_fpkm = paternal_FPKM();
		if (paternal_group_fpkm > 0.0)
		{
			double paternal_FPKM_hi = paternal_FPKM() + 2 * sqrt(paternal_FPKM_variance());
			double paternal_FPKM_lo = max(0.0, paternal_FPKM() - 2 * sqrt(paternal_FPKM_variance()));
			ConfidenceInterval paternal_conf(paternal_FPKM_lo, paternal_FPKM_hi);
			paternal_FPKM_conf(paternal_conf);
		}
		else
		{
			_paternal_FPKM_variance = 0.0;
			ConfidenceInterval paternal_conf(0.0, 0.0);
			paternal_FPKM_conf(paternal_conf);
		}
	}
	else
	{
		double sum_paternal_transfrag_FPKM_hi = 0;
        double max_paternal_fpkm = 0.0;
        //double min_paternal_fpkm = 1e100;
        
        double avg_X_g = 0.0;
        double avg_mass_fraction = 0.0;
        
        int N = _abundances.size();
        
        vector<double> avg_paternal_mass_variances(_abundances.size(), 0.0);
        
        for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin();
             itr != _count_per_replicate.end();
             ++itr)
        {
            shared_ptr<ReadGroupProperties const> rg_props = itr->first;
            double scaled_mass = itr->second;
            double scaled_total_mass = rg_props->normalized_map_mass();
            avg_X_g += scaled_mass;
            shared_ptr<MassDispersionModel const> disperser = rg_props->mass_dispersion_model();
            for (size_t j = 0; j < N; ++j)
            {
                double scaled_paternal_variance;
				scaled_paternal_variance = _abundances[j]->paternal_gamma() * disperser->scale_mass_variance(scaled_mass);
                avg_paternal_mass_variances[j] += scaled_paternal_variance;
            }
            assert (disperser->scale_mass_variance(scaled_mass) != 0 || scaled_mass == 0);
            
            //assert (scaled_total_mass != 0.0);
            avg_mass_fraction += (scaled_mass / scaled_total_mass);
        }
        
        double num_replicates = _count_per_replicate.size();
        
        if (num_replicates)
        {
            avg_X_g /= num_replicates;
            avg_mass_fraction /= num_replicates;
            for (size_t j = 0; j < N; ++j)
            {
                avg_paternal_mass_variances[j] /= num_replicates;
            }
        }
                
		BOOST_FOREACH(shared_ptr<Abundance> pA, _abundances)
		{
			double paternal_FPKM_hi;
			double paternal_FPKM_lo;
			if (pA->paternal_effective_length() > 0 && avg_mass_fraction > 0)
			{
                double norm_frag_density = 1000000000;
                norm_frag_density /= pA->paternal_effective_length();
                
                norm_frag_density *= avg_mass_fraction;
                double paternal_fpkm_high = norm_frag_density;
                
                double total_mass = ((num_paternal_fragments()+num_maternal_fragments()) / avg_mass_fraction); //hope it is num_paternal_fragments()+num_maternal_fragments() rather than num_paternal_fragments()
                double paternal_fpkm_constant = 1000000000 / pA->paternal_effective_length() / total_mass;
                double paternal_var_fpkm = paternal_mass_variance() * (paternal_fpkm_constant * paternal_fpkm_constant);
                
				paternal_FPKM_hi = paternal_fpkm_high + 2 * sqrt(paternal_var_fpkm);
				paternal_FPKM_lo = 0.0;
				ConfidenceInterval paternal_conf(paternal_FPKM_lo, paternal_FPKM_hi);
				//assert (paternal_FPKM_lo <= pA->paternal_FPKM() && pA->paternal_FPKM() <= paternal_FPKM_hi);
				pA->paternal_FPKM_conf(paternal_conf);
                //pA->paternal_FPKM_variance(paternal_var_fpkm);
				max_paternal_fpkm = max(sum_paternal_transfrag_FPKM_hi, paternal_FPKM_hi);
			}
			else
			{
				paternal_FPKM_hi = 0.0;
				paternal_FPKM_lo = 0.0;
				ConfidenceInterval paternal_conf(0.0, 0.0);
				pA->paternal_FPKM_conf(paternal_conf);
                //pA->paternal_FPKM_variance(0.0);
			}
            
		}
		// In the case of a numeric failure, the groups error bars need to be
		// set such that
		paternal_FPKM_conf(ConfidenceInterval(0.0, max_paternal_fpkm + 2 * sqrt(paternal_FPKM_variance())));
	}
	//maternal
	if (maternal_status() == NUMERIC_OK)
	{
		// This will compute the maternal transcript level FPKM confidence intervals
		for (size_t j = 0; j < _abundances.size(); ++j)
		{
            long double maternal_fpkm_var = _abundances[j]->maternal_FPKM_variance();
            double maternal_FPKM_hi = 0.0;
            double maternal_FPKM_lo = 0.0;
            if (_abundances[j]->maternal_status() != NUMERIC_FAIL)
            {
                maternal_FPKM_hi = _abundances[j]->maternal_FPKM() + 2 * sqrt(maternal_fpkm_var);
                maternal_FPKM_lo = max(0.0, (double)(_abundances[j]->maternal_FPKM() - 2 * sqrt(maternal_fpkm_var)));
                if (!(maternal_FPKM_lo <= _abundances[j]->maternal_FPKM() && _abundances[j]->maternal_FPKM() <= maternal_FPKM_hi))
                {
                    //fprintf(stderr, "Error: maternal confidence intervals are illegal! var = %Lg, fpkm = %lg, lo = %lg, hi %lg, status = %d\n", maternal_fpkm_var, _abundances[j]->maternal_FPKM(), maternal_FPKM_lo, maternal_FPKM_hi, _abundances[j]->maternal_status());
                }
                assert (maternal_FPKM_lo <= _abundances[j]->maternal_FPKM() && _abundances[j]->maternal_FPKM() <= maternal_FPKM_hi);
                ConfidenceInterval maternal_conf(maternal_FPKM_lo, maternal_FPKM_hi);
                _abundances[j]->maternal_FPKM_conf(maternal_conf);
                //_abundances[j]->maternal_FPKM_variance(maternal_fpkm_var);
            }
            else
            {
                // we shouldn't be able to get here
                assert(false);
                // TODO: nothing to do here?
            }
		}
		
		// Now build a confidence interval for the whole abundance group
		double maternal_group_fpkm = maternal_FPKM();
		if (maternal_group_fpkm > 0.0)
		{
			double maternal_FPKM_hi = maternal_FPKM() + 2 * sqrt(maternal_FPKM_variance());
			double maternal_FPKM_lo = max(0.0, maternal_FPKM() - 2 * sqrt(maternal_FPKM_variance()));
			ConfidenceInterval maternal_conf(maternal_FPKM_lo, maternal_FPKM_hi);
			paternal_FPKM_conf(maternal_conf);
		}
		else
		{
			_maternal_FPKM_variance = 0.0;
			ConfidenceInterval maternal_conf(0.0, 0.0);
			maternal_FPKM_conf(maternal_conf);
		}
	}
	else
	{
		double sum_maternal_transfrag_FPKM_hi = 0;
        double max_maternal_fpkm = 0.0;
        //double min_maternal_fpkm = 1e100;
        
        double avg_X_g = 0.0;
        double avg_mass_fraction = 0.0;
        
        int N = _abundances.size();
        
        vector<double> avg_maternal_mass_variances(_abundances.size(), 0.0);
        
        for (map<shared_ptr<ReadGroupProperties const>, double>::iterator itr = _count_per_replicate.begin();
             itr != _count_per_replicate.end();
             ++itr)
        {
            shared_ptr<ReadGroupProperties const> rg_props = itr->first;
            double scaled_mass = itr->second;
            double scaled_total_mass = rg_props->normalized_map_mass();
            avg_X_g += scaled_mass;
            shared_ptr<MassDispersionModel const> disperser = rg_props->mass_dispersion_model();
            for (size_t j = 0; j < N; ++j)
            {
                double scaled_maternal_variance;
				scaled_maternal_variance = _abundances[j]->maternal_gamma() * disperser->scale_mass_variance(scaled_mass);
                avg_maternal_mass_variances[j] += scaled_maternal_variance;
            }
            assert (disperser->scale_mass_variance(scaled_mass) != 0 || scaled_mass == 0);
            
            //assert (scaled_total_mass != 0.0);
            avg_mass_fraction += (scaled_mass / scaled_total_mass);
        }
        
        double num_replicates = _count_per_replicate.size();
        
        if (num_replicates)
        {
            avg_X_g /= num_replicates;
            avg_mass_fraction /= num_replicates;
            for (size_t j = 0; j < N; ++j)
            {
                avg_maternal_mass_variances[j] /= num_replicates;
            }
        }
                
		BOOST_FOREACH(shared_ptr<Abundance> pA, _abundances)
		{
			double maternal_FPKM_hi;
			double maternal_FPKM_lo;
			if (pA->maternal_effective_length() > 0 && avg_mass_fraction > 0)
			{
                double norm_frag_density = 1000000000;
                norm_frag_density /= pA->maternal_effective_length();
                
                norm_frag_density *= avg_mass_fraction;
                double maternal_fpkm_high = norm_frag_density;
                
                double total_mass = ((num_paternal_fragments()+num_maternal_fragments()) / avg_mass_fraction); //hope it is num_paternal_fragments()+num_maternal_fragments() rather than num_maternal_fragments()
                double maternal_fpkm_constant = 1000000000 / pA->maternal_effective_length() / total_mass;
                double maternal_var_fpkm = maternal_mass_variance() * (maternal_fpkm_constant * maternal_fpkm_constant);
                
				maternal_FPKM_hi = maternal_fpkm_high + 2 * sqrt(maternal_var_fpkm);
				maternal_FPKM_lo = 0.0;
				ConfidenceInterval maternal_conf(maternal_FPKM_lo, maternal_FPKM_hi);
				//assert (maternal_FPKM_lo <= pA->maternal_FPKM() && pA->maternal_FPKM() <= maternal_FPKM_hi);
				pA->maternal_FPKM_conf(maternal_conf);
                //pA->maternal_FPKM_variance(maternal_var_fpkm);
				max_maternal_fpkm = max(sum_maternal_transfrag_FPKM_hi, maternal_FPKM_hi);
			}
			else
			{
				maternal_FPKM_hi = 0.0;
				maternal_FPKM_lo = 0.0;
				ConfidenceInterval maternal_conf(0.0, 0.0);
				pA->maternal_FPKM_conf(maternal_conf);
                //pA->maternal_FPKM_variance(0.0);
			}
            
		}
		// In the case of a numeric failure, the groups error bars need to be
		// set such that
		maternal_FPKM_conf(ConfidenceInterval(0.0, max_maternal_fpkm + 2 * sqrt(maternal_FPKM_variance())));
	}
}

/*
void AlleleAbundanceGroup::calculate_conf_intervals()
{
    for (size_t j = 0; j < _abundances.size(); ++j)
    {
        double paternal_FPKM_hi = 0.0;
		double maternal_FPKM_hi = 0.0;
        double paternal_FPKM_lo = numeric_limits<double>::max();
		double maternal_FPKM_lo = numeric_limits<double>::max();
        const vector<double> paternal_ab_j_samples = _abundances[j]->paternal_fpkm_samples();
		const vector<double> maternal_ab_j_samples = _abundances[j]->maternal_fpkm_samples();
        vector<pair<double, double> > paternal_fpkm_samples,maternal_fpkm_samples;
        for (size_t i = 0; i < paternal_ab_j_samples.size(); ++i)
            paternal_fpkm_samples.push_back(make_pair(abs(paternal_ab_j_samples[i] - _abundances[j]->paternal_FPKM()), paternal_ab_j_samples[i]));
		for (size_t i = 0; i < maternal_ab_j_samples.size(); ++i)
            maternal_fpkm_samples.push_back(make_pair(abs(maternal_ab_j_samples[i] - _abundances[j]->maternal_FPKM()), maternal_ab_j_samples[i]));
        
        sort(paternal_fpkm_samples.begin(), paternal_fpkm_samples.end());
		sort(maternal_fpkm_samples.begin(), maternal_fpkm_samples.end());
        
        for (size_t i = 0; i < 0.95*paternal_fpkm_samples.size(); ++i)
        {
            if (paternal_FPKM_lo > paternal_fpkm_samples[i].second)
                paternal_FPKM_lo = paternal_fpkm_samples[i].second;
            if (paternal_FPKM_hi < paternal_fpkm_samples[i].second)
                paternal_FPKM_hi = paternal_fpkm_samples[i].second;
        }
		for (size_t i = 0; i < 0.95*maternal_fpkm_samples.size(); ++i)
        {
            if (maternal_FPKM_lo > maternal_fpkm_samples[i].second)
                maternal_FPKM_lo = maternal_fpkm_samples[i].second;
            if (maternal_FPKM_hi < maternal_fpkm_samples[i].second)
                maternal_FPKM_hi = maternal_fpkm_samples[i].second;
        }       
        
        ConfidenceInterval paternal_conf(paternal_FPKM_lo, paternal_FPKM_hi);
        _abundances[j]->paternal_FPKM_conf(paternal_conf);
		ConfidenceInterval maternal_conf(maternal_FPKM_lo, maternal_FPKM_hi);
        _abundances[j]->maternal_FPKM_conf(maternal_conf);
    }
    
    double paternal_FPKM_hi = 0.0;
	double maternal_FPKM_hi = 0.0;
    double paternal_FPKM_lo = numeric_limits<double>::max();
	double maternal_FPKM_lo = numeric_limits<double>::max();
    const vector<double> paternal_ab_j_samples = _paternal_fpkm_samples;
	const vector<double> maternal_ab_j_samples = _maternal_fpkm_samples;
    vector<pair<double, double> > paternal_fpkm_samples,maternal_fpkm_samples;
    for (size_t i = 0; i < paternal_ab_j_samples.size(); ++i)
        paternal_fpkm_samples.push_back(make_pair(abs(paternal_ab_j_samples[i] - paternal_FPKM()), paternal_ab_j_samples[i]));
    for (size_t i = 0; i < maternal_ab_j_samples.size(); ++i)
        maternal_fpkm_samples.push_back(make_pair(abs(maternal_ab_j_samples[i] - maternal_FPKM()), maternal_ab_j_samples[i]));
    
    sort(paternal_fpkm_samples.begin(), paternal_fpkm_samples.end());
	sort(maternal_fpkm_samples.begin(), maternal_fpkm_samples.end());
    
    for (size_t i = 0; i < 0.95*paternal_fpkm_samples.size(); ++i)
    {
        if (paternal_FPKM_lo > paternal_fpkm_samples[i].second)
            paternal_FPKM_lo = paternal_fpkm_samples[i].second;
        if (paternal_FPKM_hi < paternal_fpkm_samples[i].second)
            paternal_FPKM_hi = paternal_fpkm_samples[i].second;
    }
    for (size_t i = 0; i < 0.95*maternal_fpkm_samples.size(); ++i)
    {
        if (maternal_FPKM_lo > maternal_fpkm_samples[i].second)
            maternal_FPKM_lo = maternal_fpkm_samples[i].second;
        if (maternal_FPKM_hi < maternal_fpkm_samples[i].second)
            maternal_FPKM_hi = maternal_fpkm_samples[i].second;
    }
    ConfidenceInterval paternal_conf(paternal_FPKM_lo, paternal_FPKM_hi);
    paternal_FPKM_conf(paternal_conf);
	ConfidenceInterval maternal_conf(maternal_FPKM_lo, maternal_FPKM_hi);
    maternal_FPKM_conf(maternal_conf);
    
    return;
}
*/

void compute_cond_probs_and_effective_lengths_allele(const vector<MateHit>& alignments,
													 vector<shared_ptr<Abundance> >& transcripts,
													 vector<shared_ptr<Abundance> >& mapped_transcripts)
{		
	int N = transcripts.size();
	int M = alignments.size();

	vector<vector<char> > paternal_compatibilities(N, vector<char>(M,0));
	vector<vector<char> > maternal_compatibilities(N, vector<char>(M,0));
	compute_compatibilities(transcripts, alignments, paternal_compatibilities, maternal_compatibilities);

	for (int j = 0; j < N; ++j) 
    {
		shared_ptr<Scaffold> transfrag = transcripts[j]->transfrag();
		vector<double>& paternal_cond_probs = *(new vector<double>(M,0));
		vector<double>& maternal_cond_probs = *(new vector<double>(M,0));
		
		BiasCorrectionHelper bch(transfrag);
		
		for(int i = 0 ; i < M; ++i)
		{
			if (paternal_compatibilities[j][i]==1)
            {
                total_paternal_cond_prob_calls++;
				paternal_cond_probs[i] = bch.get_cond_prob(alignments[i]);
            }
			if (maternal_compatibilities[j][i]==1)
            {
                total_maternal_cond_prob_calls++;
				maternal_cond_probs[i] = bch.get_cond_prob(alignments[i]);
            }
		}
		double bch_eff_len = bch.get_effective_length();
		if(transcripts[j]->paternal_effective_length() == transcripts[j]->maternal_effective_length()){
			transcripts[j]->paternal_effective_length(bch_eff_len);
			transcripts[j]->maternal_effective_length(bch_eff_len);
		}
		else{//update according to their proportions
			if(transcripts[j]->paternal_effective_length() > transcripts[j]->maternal_effective_length()){
				transcripts[j]->paternal_effective_length(bch_eff_len);
				transcripts[j]->maternal_effective_length((transcripts[j]->maternal_effective_length()/transcripts[j]->paternal_effective_length())*bch_eff_len);
			}
			else{
				transcripts[j]->paternal_effective_length((transcripts[j]->paternal_effective_length()/transcripts[j]->maternal_effective_length())*bch_eff_len);
				transcripts[j]->maternal_effective_length(bch_eff_len);
			}
		}

		transcripts[j]->paternal_cond_probs(&paternal_cond_probs);
		transcripts[j]->maternal_cond_probs(&maternal_cond_probs);
		
		if (bch.is_mapped()) 
			mapped_transcripts.push_back(transcripts[j]);
	}
}

bool AlleleAbundanceGroup::calculate_gammas(const vector<MateHit>& nr_alignments, 
                                      const vector<double>& log_conv_factors,
									  const vector<shared_ptr<Abundance> >& transcripts, 
									  const vector<shared_ptr<Abundance> >& mapped_transcripts)
{
	if (mapped_transcripts.empty())
    {
		//gammas = vector<double>(transfrags.size(), 0.0);
		BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
		{
			ab->paternal_gamma(0);
			ab->maternal_gamma(0);
		}
		_parental_gamma_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                       2*transcripts.size());
        _parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                       2*transcripts.size());
        _iterated_exp_parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                       2*transcripts.size());
        _parental_fpkm_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                      2*transcripts.size());
        
		return true;
    }
	
	vector<double> paternal_gammas,maternal_gammas;
	
	verbose_msg( "Calculating intial MLE\n");
	
	AbundanceStatus mle_success = gamma_mle(mapped_transcripts,
											nr_alignments,
											log_conv_factors,
											paternal_gammas,
											maternal_gammas);
	
	verbose_msg( "Tossing likely garbage isoforms\n");
	
	for (size_t i = 0; i < paternal_gammas.size(); ++i)
	{
		if (isnan(paternal_gammas[i]))
		{
			verbose_msg("Warning: paternal isoform abundance is NaN!\n");
		}
	}
	for (size_t i = 0; i < maternal_gammas.size(); ++i)
	{
		if (isnan(maternal_gammas[i]))
		{
			verbose_msg("Warning: maternal isoform abundance is NaN!\n");
		}
	}
	
    double locus_mass = 0.0;
   
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        const MateHit& alignment = nr_alignments[i];
        locus_mass += alignment.collapse_mass();
    }
    
	vector<shared_ptr<Abundance> > filtered_transcripts = mapped_transcripts;
	vector<double> paternal_filtered_gammas = paternal_gammas;
	vector<double> maternal_filtered_gammas = maternal_gammas;
	filter_junk_isoforms(filtered_transcripts, paternal_filtered_gammas, maternal_filtered_gammas, mapped_transcripts, locus_mass);
	
	if (filtered_transcripts.empty())
	{
		//gammas = vector<double>(transfrags.size(), 0.0);
		BOOST_FOREACH (shared_ptr<Abundance> ab, _abundances)
		{
			ab->paternal_gamma(0);
			ab->maternal_gamma(0);
		}
		_parental_gamma_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
													  2*transcripts.size());
        _parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                       2*transcripts.size());
        _iterated_exp_parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                                    2*transcripts.size());
        _parental_fpkm_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 
                                                      2*transcripts.size());
        
		return true;
	}
	
    if (filtered_transcripts.size() != mapped_transcripts.size())
    {    
        paternal_filtered_gammas.clear();
		maternal_filtered_gammas.clear();
        
        verbose_msg( "Revising MLE\n");
        
        mle_success = gamma_mle(filtered_transcripts,
								nr_alignments,
								log_conv_factors, 
								paternal_filtered_gammas,
								maternal_filtered_gammas);
    }

	for (size_t i = 0; i < paternal_filtered_gammas.size(); ++i)
	{
		if (isnan(paternal_filtered_gammas[i]))
		{
			verbose_msg("Warning: paternal isoform abundance is NaN!\n");
		}
	}
	for (size_t i = 0; i < maternal_filtered_gammas.size(); ++i)
	{
		if (isnan(maternal_filtered_gammas[i]))
		{
			verbose_msg("Warning: maternal isoform abundance is NaN!\n");
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
        ublas::vector<double> paternal_gamma_mle(paternal_filtered_gammas.size());
		std::copy(paternal_filtered_gammas.begin(), paternal_filtered_gammas.end(), paternal_gamma_mle.begin());
		ublas::vector<double> maternal_gamma_mle(maternal_filtered_gammas.size());
        std::copy(maternal_filtered_gammas.begin(), maternal_filtered_gammas.end(), maternal_gamma_mle.begin());
    }
	
	for (size_t i = 0; i < paternal_filtered_gammas.size(); ++i)
	{
		if (isnan(paternal_gammas[i]))
		{
			verbose_msg( "Warning: paternal isoform abundance is NaN!\n");
			map_success = NUMERIC_FAIL;
		}
	}
	for (size_t i = 0; i < maternal_filtered_gammas.size(); ++i)
	{
		if (isnan(maternal_gammas[i]))
		{
			verbose_msg( "Warning: maternal isoform abundance is NaN!\n");
			map_success = NUMERIC_FAIL;
		}
	}
	
	// Now we need to fill in zeros for the isoforms we filtered out of the 
	// MLE/MAP calculation
	vector<double> updated_paternal_gammas = vector<double>(N, 0.0);
	vector<double> updated_maternal_gammas = vector<double>(N, 0.0);
    
    
	ublas::matrix<double> updated_parental_gamma_cov;
	updated_parental_gamma_cov = ublas::zero_matrix<double>(2*N, 2*N);
    
    ublas::matrix<double> updated_parental_count_cov;
    updated_parental_count_cov = ublas::zero_matrix<double>(2*N, 2*N);
    ublas::matrix<double> updated_iterated_exp_parental_count_cov;
    updated_iterated_exp_parental_count_cov = ublas::zero_matrix<double>(2*N, 2*N);
    ublas::matrix<double> updated_parental_fpkm_cov;
    updated_parental_fpkm_cov = ublas::zero_matrix<double>(2*N, 2*N);
	
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
			updated_paternal_gammas[i] = paternal_filtered_gammas[scaff_present[i]];
			updated_maternal_gammas[i] = maternal_filtered_gammas[scaff_present[i]];
            //cerr << updated_gammas[i] << ",";
            
			for (size_t j = 0; j < N; ++j)
			{
				if (scaff_present[j] != N)
				{
					updated_parental_gamma_cov(i,j) = _parental_gamma_covariance(scaff_present[i],
															   scaff_present[j]);
					updated_iterated_exp_parental_count_cov(i,j) = _iterated_exp_parental_count_covariance(scaff_present[i],
                                                                                         scaff_present[j]);
                    // Should still be empty but let's do these for consistency:
                    updated_parental_count_cov(i,j) = _parental_count_covariance(scaff_present[i],
                                                               scaff_present[j]);
                    updated_parental_fpkm_cov(i,j) = _parental_fpkm_covariance(scaff_present[i],
                                                             scaff_present[j]);
					assert (!isinf(updated_parental_gamma_cov(i,j)));
					assert (!isnan(updated_parental_gamma_cov(i,j)));
					
					updated_parental_gamma_cov(i,j+N) = _parental_gamma_covariance(scaff_present[i],
															   scaff_present[j]+N);
					updated_iterated_exp_parental_count_cov(i,j+N) = _iterated_exp_parental_count_covariance(scaff_present[i],
                                                                                         scaff_present[j]+N);
                    // Should still be empty but let's do these for consistency:
                    updated_parental_count_cov(i,j+N) = _parental_count_covariance(scaff_present[i],
                                                               scaff_present[j]+N);
                    updated_parental_fpkm_cov(i,j+N) = _parental_fpkm_covariance(scaff_present[i],
                                                             scaff_present[j]+N);
					assert (!isinf(updated_parental_gamma_cov(i,j+N)));
					assert (!isnan(updated_parental_gamma_cov(i,j+N)));
					
					updated_parental_gamma_cov(i+N,j+N) = _parental_gamma_covariance(scaff_present[i]+N,
															   scaff_present[j]+N);
					updated_iterated_exp_parental_count_cov(i+N,j+N) = _iterated_exp_parental_count_covariance(scaff_present[i]+N,
                                                                                         scaff_present[j]+N);
                    // Should still be empty but let's do these for consistency:
                    updated_parental_count_cov(i+N,j+N) = _parental_count_covariance(scaff_present[i]+N,
                                                               scaff_present[j]+N);
                    updated_parental_fpkm_cov(i+N,j+N) = _parental_fpkm_covariance(scaff_present[i]+N,
                                                             scaff_present[j]+N);
					assert (!isinf(updated_parental_gamma_cov(i+N,j+N)));
					assert (!isnan(updated_parental_gamma_cov(i+N,j+N)));
					
					updated_parental_gamma_cov(i+N,j) = _parental_gamma_covariance(scaff_present[i]+N,
															   scaff_present[j]);
					updated_iterated_exp_parental_count_cov(i+N,j) = _iterated_exp_parental_count_covariance(scaff_present[i]+N,
                                                                                         scaff_present[j]);
                    // Should still be empty but let's do these for consistency:
                    updated_parental_count_cov(i+N,j) = _parental_count_covariance(scaff_present[i]+N,
                                                               scaff_present[j]);
                    updated_parental_fpkm_cov(i+N,j) = _parental_fpkm_covariance(scaff_present[i]+N,
                                                             scaff_present[j]);
					assert (!isinf(updated_parental_gamma_cov(i+N,j)));
					assert (!isnan(updated_parental_gamma_cov(i+N,j)));
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
		_abundances[i]->paternal_gamma(updated_paternal_gammas[i]);
		_abundances[i]->maternal_gamma(updated_maternal_gammas[i]);
		_abundances[i]->paternal_status(numeric_status);
		_abundances[i]->maternal_status(numeric_status);
	}
	_parental_gamma_covariance = updated_parental_gamma_cov;
    _parental_count_covariance = updated_parental_count_cov;
    _iterated_exp_parental_count_covariance = updated_iterated_exp_parental_count_cov;
    _parental_fpkm_covariance = updated_parental_fpkm_cov;
	
	return ((paternal_status() == NUMERIC_OK) && (maternal_status() == NUMERIC_OK));
}

void calculate_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                const Eigen::MatrixXd& paternal_transcript_cond_probs,
								const Eigen::MatrixXd& maternal_transcript_cond_probs,
                                const Eigen::VectorXd& paternal_proposed_gammas,
								const Eigen::VectorXd& maternal_proposed_gammas,
                                Eigen::MatrixXd& paternal_assignment_probs,
								Eigen::MatrixXd& maternal_assignment_probs)
{
//    vector<double> u(nr_alignments.size());
//    for (size_t i = 0; i < alignment_multiplicity.size(); ++i)
//    {
//        u[i] = nr_alignments[i].collapse_mass();
//    }
    
    //ublas::vector<double> total_cond_prob = ublas::prod(proposed_gammas,transcript_cond_probs);
	Eigen::VectorXd total_cond_prob = (paternal_proposed_gammas.transpose() * paternal_transcript_cond_probs) + (maternal_proposed_gammas.transpose() * maternal_transcript_cond_probs);
    // Compute the marginal conditional probability for each fragment against each isoform
    //ublas::matrix<double>  marg_cond_prob = ublas::zero_matrix<double>(transcript_cond_probs.size1(), transcript_cond_probs.size2());
	Eigen::MatrixXd paternal_marg_cond_prob(paternal_transcript_cond_probs.rows(), paternal_transcript_cond_probs.cols());
	Eigen::MatrixXd maternal_marg_cond_prob(maternal_transcript_cond_probs.rows(), maternal_transcript_cond_probs.cols());
    
    for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
        paternal_marg_cond_prob.array().col(i) = paternal_proposed_gammas.array() * paternal_transcript_cond_probs.array().col(i);
		maternal_marg_cond_prob.array().col(i) = maternal_proposed_gammas.array() * maternal_transcript_cond_probs.array().col(i);
        
        if (total_cond_prob(i) > 0)
        {
			paternal_marg_cond_prob.array().col(i) /= total_cond_prob(i);
			//column(paternal_marg_cond_prob,i) /= total_cond_prob(i);
			maternal_marg_cond_prob.array().col(i) /= total_cond_prob(i);
			//column(maternal_marg_cond_prob,i) /= total_cond_prob(i);
        }
    }
    paternal_assignment_probs = paternal_marg_cond_prob;
	maternal_assignment_probs = maternal_marg_cond_prob;
}

void calculate_average_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                        const Eigen::MatrixXd& paternal_transcript_cond_probs,
										const Eigen::MatrixXd& maternal_transcript_cond_probs,
                                        const Eigen::VectorXd& proposed_paternal_gammas,
										const Eigen::VectorXd& proposed_maternal_gammas,
                                        Eigen::MatrixXd& paternal_assignment_probs,
										Eigen::MatrixXd& maternal_assignment_probs)
{
	//ublas::vector<double> total_cond_prob = ublas::prod(proposed_gammas,transcript_cond_probs);
	Eigen::VectorXd total_cond_prob = (proposed_paternal_gammas.transpose() * paternal_transcript_cond_probs) +  (proposed_maternal_gammas.transpose() * maternal_transcript_cond_probs);

	// Compute the marginal conditional probability for each fragment against each isoform
    //ublas::matrix<double>  paternal_marg_cond_prob = ublas::zero_matrix<double>(paternal_transcript_cond_probs.size1(), paternal_transcript_cond_probs.size2());
	//ublas::matrix<double>  maternal_marg_cond_prob = ublas::zero_matrix<double>(maternal_transcript_cond_probs.size1(), maternal_transcript_cond_probs.size2());
    Eigen::MatrixXd marg_paternal_cond_prob(paternal_transcript_cond_probs.rows(), paternal_transcript_cond_probs.cols());
	Eigen::MatrixXd marg_maternal_cond_prob(maternal_transcript_cond_probs.rows(), maternal_transcript_cond_probs.cols());
	
	for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
		marg_paternal_cond_prob.array().col(i) = proposed_paternal_gammas.array() * paternal_transcript_cond_probs.array().col(i);
		marg_maternal_cond_prob.array().col(i) = proposed_maternal_gammas.array() * maternal_transcript_cond_probs.array().col(i);
        
        if (total_cond_prob(i) > 0)
        {
            marg_paternal_cond_prob.array().col(i) /= total_cond_prob(i);
			marg_maternal_cond_prob.array().col(i) /= total_cond_prob(i);
            //column(marg_paternal_cond_prob,i) /= total_cond_prob(i);
			//column(marg_maternal_cond_prob,i) /= total_cond_prob(i);
        }
    }
    paternal_assignment_probs = Eigen::MatrixXd::Zero(proposed_paternal_gammas.size(), proposed_paternal_gammas.size());
	maternal_assignment_probs = Eigen::MatrixXd::Zero(proposed_maternal_gammas.size(), proposed_maternal_gammas.size());
    
	vector<double> num_paternal_compatible(proposed_paternal_gammas.size(), 0);
	vector<double> num_maternal_compatible(proposed_maternal_gammas.size(), 0);
          
    for (size_t i = 0; i < alignment_multiplicities.size(); ++i)
    {
        for (size_t j = 0; j < proposed_paternal_gammas.size(); ++j)
        {
            if (marg_paternal_cond_prob(j,i) > 0)
            {
                paternal_assignment_probs.col(j) += marg_paternal_cond_prob.col(i) * alignment_multiplicities[i];
                num_paternal_compatible[j] += alignment_multiplicities[i];
            }
        }
		for (size_t j = 0; j < proposed_maternal_gammas.size(); ++j)
        {
            if (marg_maternal_cond_prob(j,i) > 0)
            {
                maternal_assignment_probs.col(j) += marg_maternal_cond_prob.col(i) * alignment_multiplicities[i];
                num_maternal_compatible[j] += alignment_multiplicities[i];
            }
        }
    }
    
//    cerr << "marg cond prob:" << endl;
//    cerr << marg_cond_prob << endl;
    
    for (size_t j = 0; j < proposed_paternal_gammas.size(); ++j)
    {
        if (num_paternal_compatible[j] > 0)
            paternal_assignment_probs.col(j) /= num_paternal_compatible[j];
    }
	for (size_t j = 0; j < proposed_maternal_gammas.size(); ++j)
    {
        if (num_maternal_compatible[j] > 0)
            maternal_assignment_probs.col(j) /= num_maternal_compatible[j];
    }
    
//    cerr << "multiplicities:" << endl;
//    cerr << alignment_multiplicities << endl;
//    cerr << "avg matrix:" << endl;
//    cerr << assignment_probs << endl;
    //assignment_probs = marg_cond_prob;
}

void calculate_iterated_exp_count_covariance(const vector<double>& paternal_gammas,
											 const vector<double>& maternal_gammas,
                                             const vector<MateHit>& nr_alignments,
                                             const vector<shared_ptr<Abundance> >& transcripts,
                                             ublas::matrix<double>& parental_count_covariance)
{
    // Now calculate the _iterated_exp_count_covariance matrix via iterated expectation
    vector<vector<double> > paternal_cond_probs(transcripts.size(), vector<double>());
	vector<vector<double> > maternal_cond_probs(transcripts.size(), vector<double>());
    for(size_t j = 0; j < transcripts.size(); ++j)
    {
        paternal_cond_probs[j]= *(transcripts[j]->paternal_cond_probs());
		maternal_cond_probs[j]= *(transcripts[j]->maternal_cond_probs());
	}
    	
    vector<double> u(nr_alignments.size());
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        u[i] = nr_alignments[i].collapse_mass();
    }
    
    parental_count_covariance = ublas::zero_matrix<double>(2*transcripts.size(), 2*transcripts.size());
    
    ublas::vector<double> total_cond_prob = ublas::zero_vector<double>(nr_alignments.size());
    
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        // the replicate gamma mles might not be available, if one of the
        // replicates returned an error, we'll consider all to be unreliable
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
			if (paternal_cond_probs[j][i] > 0 || maternal_cond_probs[j][i] > 0)
            {
                total_cond_prob(i) += (transcripts[j]->paternal_gamma() * paternal_cond_probs[j][i]) + (transcripts[j]->maternal_gamma() * maternal_cond_probs[j][i]);
                assert (!isnan(total_cond_prob(i) && ! isinf(total_cond_prob(i))));
            }
        }
    }
	
    // Compute the marginal conditional probability for each fragment against each isoform
    ublas::matrix<double>  marg_paternal_cond_prob = ublas::zero_matrix<double>(transcripts.size(), nr_alignments.size());
	ublas::matrix<double>  marg_maternal_cond_prob = ublas::zero_matrix<double>(transcripts.size(), nr_alignments.size());
    
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
        // the replicate gamma mles might not be available, if one of the
        // replicates returned an error, we'll consider all to be unreliable
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
            if (total_cond_prob(i))
            {
				marg_paternal_cond_prob(j,i) = (transcripts[j]->paternal_gamma() * paternal_cond_probs[j][i])/total_cond_prob(i);
				marg_maternal_cond_prob(j,i) = (transcripts[j]->maternal_gamma() * maternal_cond_probs[j][i])/total_cond_prob(i);
			}
        }
    }
    
    double total_paternal_paternal_var = 0.0;
	double total_maternal_maternal_var = 0.0;
    
    ublas::vector<double> paternal_expected_counts = ublas::zero_vector<double>(paternal_cond_probs.size());
	ublas::vector<double> maternal_expected_counts = ublas::zero_vector<double>(maternal_cond_probs.size());

    //iterate over fragments
    for (size_t i = 0; i < nr_alignments.size(); ++i)
    {
		// iterate over transcripts
        for (size_t j = 0; j < transcripts.size(); ++j)
        {
			double paternal_c_j_i = marg_paternal_cond_prob(j,i);
			double maternal_c_j_i = marg_maternal_cond_prob(j,i);
            
			paternal_expected_counts(j) += u[i] * marg_paternal_cond_prob(j,i);
			maternal_expected_counts(j) += u[i] * marg_maternal_cond_prob(j,i);
            
            //if (c_j_i == 0 || c_j_i == 1.0)
            //    continue;
            
            for (size_t k = 0; k < transcripts.size(); ++k)
            {
				double paternal_c_k_i = marg_paternal_cond_prob(k,i);
				double maternal_c_k_i = marg_maternal_cond_prob(k,i);
								
				//if (c_k_i == 0 || c_k_i == 1.0)
                //    continue;
				if (j == k)
                {
					//fill in the diagonal elements (paternal[k]:paternal[k] and maternal[k]:maternal[k])
					if (paternal_c_k_i != 0 && paternal_c_k_i != 1.0)
                    {
                        double paternal_paternal_var = u[i] * paternal_c_k_i * (1.0 - paternal_c_k_i);
                        parental_count_covariance(k,k) += paternal_paternal_var;
                        assert (paternal_paternal_var >= 0);
                        assert (!isnan(paternal_paternal_var) && !isinf(paternal_paternal_var));
                        total_paternal_paternal_var += paternal_paternal_var;
                    }
					if (maternal_c_k_i != 0 && maternal_c_k_i != 1.0)
                    {
                        double maternal_maternal_var = u[i] * maternal_c_k_i * (1.0 - maternal_c_k_i);
                        parental_count_covariance(k+transcripts.size(),k+transcripts.size()) += maternal_maternal_var;
                        assert (maternal_maternal_var >= 0);
                        assert (!isnan(maternal_maternal_var) && !isinf(maternal_maternal_var));
                        total_maternal_maternal_var += maternal_maternal_var;
                    }
					//fill in the off diagonal elements (paternal[k]:maternal[k] and maternal[k]:paternal[k])
					if (paternal_c_k_i != 0 && paternal_c_k_i != 1.0 && maternal_c_k_i != 0 && maternal_c_k_i != 1.0)
                    {
                        double covar = -u[i] * paternal_c_k_i * maternal_c_k_i;
						assert (covar >= 0);
                        assert (!isnan(covar) && !isinf(covar));
                        parental_count_covariance(k,k+transcripts.size()) += covar;
						parental_count_covariance(k+transcripts.size(),k) += covar;
					}					
                }
                else
                {
					if (paternal_c_k_i != 0 && paternal_c_k_i != 1.0 && paternal_c_j_i != 0 && paternal_c_j_i != 1.0)
					{
						double paternal_paternal_covar = -u[i] * paternal_c_k_i * paternal_c_j_i;
						assert (paternal_paternal_covar <= 0);
						assert (!isnan(paternal_paternal_covar));
						parental_count_covariance(k,j) += paternal_paternal_covar;
					}
					if (paternal_c_k_i != 0 && paternal_c_k_i != 1.0 && maternal_c_j_i != 0 && maternal_c_j_i != 1.0)
					{
						double paternal_maternal_covar = -u[i] * paternal_c_k_i * maternal_c_j_i;
						assert (paternal_maternal_covar <= 0);
						assert (!isnan(paternal_maternal_covar));
						parental_count_covariance(k,j+transcripts.size()) += paternal_maternal_covar;
					}
					if (maternal_c_k_i != 0 && maternal_c_k_i != 1.0 && maternal_c_j_i != 0 && maternal_c_j_i != 1.0)
					{
						double maternal_maternal_covar = -u[i] * maternal_c_k_i * maternal_c_j_i;
						assert (maternal_maternal_covar <= 0);
						assert (!isnan(maternal_maternal_covar));
						parental_count_covariance(k+transcripts.size(),j+transcripts.size()) += maternal_maternal_covar;
					}
					if (maternal_c_k_i != 0 && maternal_c_k_i != 1.0 && paternal_c_j_i != 0 && paternal_c_j_i != 1.0)
					{
						double maternal_paternal_covar = -u[i] * maternal_c_k_i * paternal_c_j_i;
						assert (maternal_paternal_covar <= 0);
						assert (!isnan(maternal_paternal_covar));
						parental_count_covariance(k+transcripts.size(),j) += maternal_paternal_covar;
					}
                } 
            }
        }

    }
	
	
    //cerr << "paternal expected counts" << endl;
    //cerr << paternal_expected_counts << endl;
	//cerr << "maternal expected counts" << endl;
    //cerr << maternal_expected_counts << endl;
    
//    // take care of little rounding errors
//    for (size_t i = 0; i < _iterated_exp_parental_count_covariance.size1(); ++i)
//    {
//        for (size_t j = 0; j < _iterated_exp_parental_count_covariance.size2(); ++j)
//        {
//            if (i == j)
//            {
//                double c = _iterated_exp_parental_count_covariance(i,j);
//                if (c < 0)
//                    _iterated_exp_parental_count_covariance(i,j) = 0;
//                //assert(c >= 0);
//            }
//            else
//            {
//                double c = _iterated_exp_parental_count_covariance(i,j);
//                if (c > 0)
//                    _iterated_exp_parental_count_covariance(i,j) = 0;
//                //assert(c <= 0);
//            }
//        }
//	      if(i < _iterated_exp_parental_count_covariance.size1()/2)
//		  {
//			  _abundances[i]->num_paternal_fragment_uncertainty_var(_iterated_exp_parental_count_covariance(i,i));
//		  }
//		  else
//		  {
//			  _abundances[i-_iterated_exp_parental_count_covariance.size1()/2]->num_maternal_fragment_uncertainty_var(_iterated_exp_parental_count_covariance(i,i));
//		  }
//    }
}

void AlleleAbundanceGroup::calculate_kappas()
{
    size_t num_members = 2*_abundances.size();
    _parental_kappa_covariance = ublas::matrix<double>(num_members, 
											  num_members);
	//cerr << gamma_cov <<endl;
	
	assert (_parental_gamma_covariance.size1() == num_members);
	assert (_parental_gamma_covariance.size2() == num_members);
		
	//tss_group.sub_quants = vector<QuantGroup>(isos_in_tss);
	
	double paternal_S_FPKM = 0.0;
	double maternal_S_FPKM = 0.0;
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		if (pA->paternal_effective_length() > 0 || pA->maternal_effective_length() > 0)
		{
			paternal_S_FPKM += pA->paternal_FPKM();
			maternal_S_FPKM += pA->maternal_FPKM();
		}
	}
	
    //fprintf (stderr, "*********\n");
	BOOST_FOREACH (shared_ptr<Abundance> pA, _abundances)
	{
		if (paternal_S_FPKM+maternal_S_FPKM > 0)
		{
			pA->paternal_kappa(pA->paternal_FPKM() / (paternal_S_FPKM+maternal_S_FPKM));
			pA->maternal_kappa(pA->maternal_FPKM() / (paternal_S_FPKM+maternal_S_FPKM));
        }
		else
		{
			pA->paternal_kappa(0); 
			pA->maternal_kappa(0);
		}
	}
    
	for (size_t k = 0; k < num_members; ++k)
	{
		for (size_t m = 0; m < num_members; ++m)
		{
			double L = 1.0;
			if(k < num_members/2)
			{
				if(m < num_members/2)
				{
					L = _abundances[k]->paternal_effective_length() * _abundances[m]->paternal_effective_length();
				}
				else{
					L = _abundances[k]->paternal_effective_length() * _abundances[m-(num_members/2)]->maternal_effective_length();
				}
			}
			else{
				if(m < num_members/2)
				{
					L = _abundances[k-(num_members/2)]->maternal_effective_length() * _abundances[m]->paternal_effective_length();
				}
				else{
					L = _abundances[k-(num_members/2)]->maternal_effective_length() * _abundances[m-(num_members/2)]->maternal_effective_length();
				}
			}
            if (L == 0.0)
            {
                _parental_kappa_covariance(k,m) = 0.0;
            }
            else if (m == k)
            {
                // Use the modeled count variance here instead
                double kappa_var;
                if (paternal_S_FPKM+maternal_S_FPKM)
                {
                    kappa_var = _parental_fpkm_covariance(k,k) / ((paternal_S_FPKM+maternal_S_FPKM) * (paternal_S_FPKM+maternal_S_FPKM));
                }
                else
                {
                    kappa_var = 0.0;
                }
                
                if (isnan(kappa_var) || isinf(kappa_var)) // to protect against underflow
                    kappa_var = 0;
                _parental_kappa_covariance(k,m) = kappa_var;
            }
            else
            {
                double kappa_covar;
                if (paternal_S_FPKM+maternal_S_FPKM)
                {
                    kappa_covar = _parental_fpkm_covariance(k,m) / ((paternal_S_FPKM+maternal_S_FPKM) * (paternal_S_FPKM+maternal_S_FPKM));
                }
                else
                {
                    kappa_covar = 0.0;
                }
                _parental_kappa_covariance(k,m) = kappa_covar;
            }
		}
	}

    
//	size_t num_members = 2*_abundances.size();
//	_parental_kappa_covariance = ublas::zero_matrix<double>(num_members, 
//											  num_members);
//    bools flag = true;
//    if (paternal_FPKM() == 0)
//    {
//        for (size_t k = 0; k < _abundances.size(); ++k)
//        {
//            _abundances[k]->paternal_kappa(0);
//        }
//        flag = 0;
//    }
//    if (maternal_FPKM() == 0)
//    {
//        for (size_t k = 0; k < _abundances.size(); ++k)
//        {
//            _abundances[k]->maternal_kappa(0);
//        }
//        flag = 0;
//    }
//	if(!flag) return;   


    // FIXME:
    size_t num_paternal_count_draws = _paternal_assigned_count_samples.size();
	size_t num_maternal_count_draws = _maternal_assigned_count_samples.size();
    vector<Eigen::VectorXd > relative_paternal_abundances (num_paternal_count_draws, Eigen::VectorXd::Zero(_abundances.size()));
	vector<Eigen::VectorXd > relative_maternal_abundances (num_maternal_count_draws, Eigen::VectorXd::Zero(_abundances.size()));
    
    // We'll use the effective lengths to transform counts into relative abundances,
    // and then use that to calculate the kappa variances and covariances.
    Eigen::VectorXd paternal_effective_length_recip = Eigen::VectorXd::Zero(_abundances.size());
	Eigen::VectorXd maternal_effective_length_recip = Eigen::VectorXd::Zero(_abundances.size());
    for (size_t i = 0; i < _abundances.size(); ++i)
    {
        if (_abundances[i]->paternal_effective_length() > 0)
            paternal_effective_length_recip(i) = 1.0 / _abundances[i]->paternal_effective_length();
		if (_abundances[i]->maternal_effective_length() > 0)
            maternal_effective_length_recip(i) = 1.0 / _abundances[i]->maternal_effective_length();
    }
    assert(num_paternal_count_draws == num_maternal_count_draws);
    for (size_t i = 0; i < num_paternal_count_draws; ++i)
    {
        
        Eigen::VectorXd relative_paternal_abundance = paternal_effective_length_recip.array() * _paternal_assigned_count_samples[i].array();
        double paternal_total = relative_paternal_abundance.sum();
		Eigen::VectorXd relative_maternal_abundance = maternal_effective_length_recip.array() * _maternal_assigned_count_samples[i].array();
        double maternal_total = relative_maternal_abundance.sum();
        if (paternal_total+maternal_total > 0){
			relative_paternal_abundance /= (paternal_total+maternal_total);
			relative_maternal_abundance /= (paternal_total+maternal_total);
		}
		
        //cerr << relative_paternal_abundance.transpose() << endl;
		//cerr << relative_maternal_abundance.transpose() << endl;
        relative_paternal_abundances[i] = relative_paternal_abundance;
		relative_maternal_abundances[i] = relative_maternal_abundance;
    }
}

void Estep (int N, 
			int M, 
			const Eigen::VectorXd& paternal_p,
			const Eigen::VectorXd& maternal_p,
			Eigen::MatrixXd& paternal_U,
			Eigen::MatrixXd& maternal_U,
			const Eigen::MatrixXd& paternal_cond_probs,
			const Eigen::MatrixXd& maternal_cond_probs,
			const Eigen::VectorXd& alignment_multiplicities,
			const bool first) {
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
    Eigen::VectorXd paternal_frag_prob_sums;// = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd maternal_frag_prob_sums;// = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd frag_prob_sums;// = Eigen::VectorXd::Zero(M);
    
//    for (j = 0; j < N; ++j) 
//    {
//        for (i = 0; i < M; ++i) 
//        {
//            frag_prob_sums(i) += cond_probs(j,i) * p(j);
//            assert (!isnan(cond_probs(j,i)) && !isinf(cond_probs(j,i)));
//        }
//    }
	
	paternal_frag_prob_sums = paternal_p.transpose() * paternal_cond_probs;//this is a vector of size M, where each element is the col sum over all paternal N isoforms
	maternal_frag_prob_sums = maternal_p.transpose() * maternal_cond_probs;//this is a vector of size M, where each element is the col sum over all maternal N isoforms
	if(first)
	{
		for(int i = 0; i < M; ++i){
			if(paternal_cond_probs.col(i).sum() > 0 && maternal_cond_probs.col(i).sum() > 0)
			{
				paternal_frag_prob_sums[i] /= 2.0;
				maternal_frag_prob_sums[i] /= 2.0;
			}
		}
	}	
	frag_prob_sums = paternal_frag_prob_sums+maternal_frag_prob_sums;//this is a vector of size M, where each element is the combined cols sum over all paternal+maternal N isoforms
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
    
    paternal_U = Eigen::MatrixXd(N,M);
	maternal_U = Eigen::MatrixXd(N,M);

	if(first){
		for (i = 0; i < M; ++i) 
		{
			for (j = 0; j < N; ++j) 
			{
				if(paternal_cond_probs.col(i).sum() > 0 && maternal_cond_probs.col(i).sum() > 0)
				{
					paternal_U(j,i) = paternal_cond_probs(j,i) * (paternal_p[j]/2.0);
					maternal_U(j,i) = maternal_cond_probs(j,i) * (maternal_p[j]/2.0);
				}
				else{
					paternal_U(j,i) = paternal_cond_probs(j,i) * paternal_p[j];
					maternal_U(j,i) = maternal_cond_probs(j,i) * maternal_p[j];
				}
			}
		}
	}
	else{
		for (i = 0; i < M; ++i) 
		{
			//double ProbY = frag_prob_sums(i);
			//double exp_i_j = cond_probs(j,i) * p(j) * x(i);
			//U.r = exp_i_j;
			paternal_U.col(i) = paternal_cond_probs.col(i).array() * paternal_p.array();
			maternal_U.col(i) = maternal_cond_probs.col(i).array() * maternal_p.array();
		}
	}
	
	for (j = 0; j < N; ++j) 
    {
        paternal_U.array().row(j) *= x.transpose().array();
		maternal_U.array().row(j) *= x.transpose().array();
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
            Eigen::VectorXd& paternal_p, 
			Eigen::VectorXd& maternal_p, 
            const Eigen::MatrixXd& paternal_U,
			const Eigen::MatrixXd& maternal_U) 
{
	Eigen::VectorXd paternal_v;// = Eigen::VectorXd::Zero(N);
	Eigen::VectorXd maternal_v;// = Eigen::VectorXd::Zero(N);
    
	double m = 0; //m is the number of reads (paternal+maternal)
    
	paternal_v = paternal_U.rowwise().sum();
	maternal_v = maternal_U.rowwise().sum();
    m = paternal_v.colwise().sum()(0)+maternal_v.colwise().sum()(0);
    
	if (m)
	{
		paternal_p = paternal_v / m;
		maternal_p = maternal_v / m;
	}
	else
	{
        paternal_p = Eigen::VectorXd::Zero(N);
		maternal_p = Eigen::VectorXd::Zero(N);
	}
}
 

double logLike (int N, 
				int M, 
				Eigen::VectorXd& paternal_p,
				Eigen::VectorXd& maternal_p,
				const Eigen::MatrixXd& paternal_cond_prob, 
				const Eigen::MatrixXd& maternal_cond_prob, 
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
	Eigen::VectorXd Prob_Ys = paternal_p.transpose() * paternal_cond_prob + maternal_p.transpose() * maternal_cond_prob;
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
          Eigen::VectorXd&  newPaternalP, 
		  Eigen::VectorXd&  newMaternalP, 
          const Eigen::MatrixXd& paternal_cond_prob, 
		  const Eigen::MatrixXd& maternal_cond_prob, 
		  const Eigen::VectorXd& alignment_multiplicities,
          vector<double> const & log_conv_factors,
          bool& converged) 
{
	converged = true;
	//double sum = 0;
	double newEll = 0;
    //Eigen::VectorXd paternal_p = Eigen::VectorXd::Zero(N);
	//Eigen::VectorXd maternal_p = Eigen::VectorXd::Zero(N);
	Eigen::VectorXd paternal_p = Eigen::VectorXd::Zero(N,M);
	Eigen::VectorXd maternal_p = Eigen::VectorXd::Zero(N,M);
    Eigen::MatrixXd paternal_U = Eigen::MatrixXd::Zero(N,M);
	Eigen::MatrixXd maternal_U = Eigen::MatrixXd::Zero(N,M);
	double ell = 0; 
	int iter = 0;
	int j,i;
    bool first = true;
	
    if (N == 0 || M == 0)
        return NUMERIC_OK;
    
    for (j = 0; j < N; ++j) {
        //p[j] = drand48();
        //sum += p[j];
		paternal_p(j) = 1.0/(double)N;
		maternal_p(j) = 1.0/(double)N;
    }
		
    static const double ACCURACY = mle_accuracy; // convergence for EM
	
	while (((iter <= 2) || (abs(ell - newEll) > ACCURACY)) && (iter < max_mle_iterations)) {
		if (iter > 0) {
			//round(newP);
			paternal_p = newPaternalP;
			maternal_p = newMaternalP;
			ell = newEll;
		}
		
		Estep(N, M, paternal_p, maternal_p, paternal_U, maternal_U, paternal_cond_prob, maternal_cond_prob, alignment_multiplicities, first); //  fills parental U's
		first = false;
		Mstep(N, M, newPaternalP, newMaternalP, paternal_U, maternal_U); // fills parental p's
		newEll = logLike(N, M, newPaternalP, newMaternalP, paternal_cond_prob, maternal_cond_prob, alignment_multiplicities, log_conv_factors);
						
		//nimrod
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

AbundanceStatus AlleleAbundanceGroup::calculate_per_replicate_abundances(vector<shared_ptr<Abundance> >& transcripts,
																		 const map<shared_ptr<ReadGroupProperties const >, vector<MateHit> >& alignments_per_read_group,
																		 std::map<shared_ptr<ReadGroupProperties const >, shared_ptr<AlleleAbundanceGroup> >& ab_group_per_replicate,
                                                                   bool perform_collapse)
{
	for(std::set<shared_ptr<ReadGroupProperties const > >::iterator itr = _read_group_props.begin();
        itr != _read_group_props.end(); 
        ++itr)
    {
        vector<shared_ptr<Abundance> > new_transcripts;
        BOOST_FOREACH(shared_ptr<Abundance> ab, transcripts)
        {
            boost::shared_ptr<AlleleTranscriptAbundance> d = boost::static_pointer_cast<AlleleTranscriptAbundance>(ab);
            //new_transcripts.push_back(shared_ptr<Abundance>(new AlleleTranscriptAbundance(*boost::static_pointer_cast<AlleleTranscriptAbundance>(ab))));
            AlleleTranscriptAbundance* pT = new AlleleTranscriptAbundance;
            pT->transfrag(d->transfrag());
			pT->set_allele_informative();
			shared_ptr<Abundance> ab_new(pT);
            ab_new->description(ab_new->description());
            ab_new->locus_tag("");
            new_transcripts.push_back(ab_new);
        }
        shared_ptr<AlleleAbundanceGroup> ab_group(new AlleleAbundanceGroup(new_transcripts));
        std::set<shared_ptr<ReadGroupProperties const > > rg_props;
        rg_props.insert(*itr);
        ab_group->init_rg_props(rg_props);
        
        map<shared_ptr<ReadGroupProperties const >, vector<MateHit> >::const_iterator al_itr =
            alignments_per_read_group.find(*itr);
        
        assert(al_itr != alignments_per_read_group.end());
        const vector<MateHit>& rep_hits = alignments_per_read_group.find(*itr)->second;
				
        vector<MateHit> nr_alignments;
        vector<MateHit> non_equiv_alignments;
        vector<double> log_conv_factors;
        vector<shared_ptr<Abundance> > mapped_transcripts;
				
        if (perform_collapse)
        {
			collapse_hits(rep_hits, nr_alignments,true);
            collapse_equivalent_hits_helper_allele(nr_alignments, transcripts, non_equiv_alignments, log_conv_factors);
            assert (non_equiv_alignments.size() == log_conv_factors.size());
            log_conv_factors = vector<double>(nr_alignments.size(), 0);
            nr_alignments.clear();
        }
        else
        {
			nr_alignments = rep_hits;
            non_equiv_alignments = nr_alignments;
            log_conv_factors = vector<double>(nr_alignments.size(), 0);
        }
		
		ab_group->calculate_abundance_for_replicate(non_equiv_alignments, false);
        //fprintf (stderr, "FPKM = %lg\n", ab_group->FPKM());
        ab_group_per_replicate[*itr] = ab_group;			
    }
	
    return NUMERIC_OK;
}

AbundanceStatus gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                          const vector<MateHit>& nr_alignments,
                          const vector<double>& log_conv_factors,
                          vector<double>& paternal_gammas,
						  vector<double>& maternal_gammas,
                          bool check_identifiability)
{
	//paternal_gammas.clear();
	//maternal_gammas.clear();
    paternal_gammas = vector<double>(transcripts.size(), 0);
	maternal_gammas = vector<double>(transcripts.size(), 0);
	if (transcripts.empty())
		return NUMERIC_OK;
	
	/*This condition doesn't apply in the allele case
	if (transcripts.size() == 1)
	{
		gammas = vector<double>(1,1.0);
		return NUMERIC_OK;
	}
	*/
	size_t M = nr_alignments.size();
	size_t N = transcripts.size();
	
    bool converged = true;
    bool identifiable = true;
    
	if (M > 0)
	{
		
		//vector<vector<double> > saliencies (M,vector<double>(N,0));
		
		
		//compute_saliencies(cond_probs, saliencies, saliency_weight);
		
        Eigen::VectorXd paternal_prob = Eigen::VectorXd::Zero(N);
		Eigen::VectorXd maternal_prob = Eigen::VectorXd::Zero(N);
        Eigen::VectorXd _phint = Eigen::VectorXd::Zero(N);
        
		double logL;
		
        Eigen::MatrixXd paternal_cond_probs(N, M);
		Eigen::MatrixXd maternal_cond_probs(N, M);
		
        for (size_t j = 0; j < paternal_cond_probs.rows(); ++j)
        {
            for (size_t i = 0; i < paternal_cond_probs.cols(); ++i)
            {
                paternal_cond_probs(j,i) = (*(transcripts[j]->paternal_cond_probs()))[i];
            }
        }
		for (size_t j = 0; j < maternal_cond_probs.rows(); ++j)
        {
            for (size_t i = 0; i < maternal_cond_probs.cols(); ++i)
            {
                maternal_cond_probs(j,i) = (*(transcripts[j]->maternal_cond_probs()))[i];
            }
        }
        
        if (check_identifiability)
        {
            ublas::matrix<double> compat = ublas::zero_matrix<double>(M,N);
            
            for (size_t j = 0; j < N; ++j)//transcripts
            {
                for (size_t i = 0; i < M; ++i)//alignments
                {
                    if (paternal_cond_probs(j,i) || maternal_cond_probs(j,i))
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
		
        logL = EM(N, M, paternal_prob, maternal_prob, paternal_cond_probs, maternal_cond_probs, alignment_multiplicities, log_conv_factors, converged);
        
		paternal_gammas = vector<double>(N, 0.0);
		maternal_gammas = vector<double>(N, 0.0);
		
		for (size_t i = 0; i < paternal_gammas.size(); ++i)
		{
            paternal_gammas[i] = paternal_prob(i);
			if (isnan(paternal_gammas[i]) || isinf(paternal_gammas[i]))
            {
                return NUMERIC_FAIL;
            }
		}
		for (size_t i = 0; i < maternal_gammas.size(); ++i)
		{
            maternal_gammas[i] = maternal_prob(i);
			if (isnan(maternal_gammas[i]) || isinf(maternal_gammas[i]))
            {
                return NUMERIC_FAIL;
            }
		}
	}
	else
	{
		paternal_gammas = vector<double>(N, 0.0);
		maternal_gammas = vector<double>(N, 0.0);
	}
	
    double paternal_round_err = 0.0;
	double maternal_round_err = 0.0;
    double num_good = 0;
    for(size_t g = 0;g < N;++g)
    {
        if (paternal_gammas[g]+maternal_gammas[g] < min_isoform_fraction)
        {
            paternal_round_err += paternal_gammas[g];
			maternal_round_err += maternal_gammas[g];
			paternal_gammas[g] = 0.0;
			maternal_gammas[g] = 0.0;
        }
        else
        {
            num_good += 1;
        }
    }
	for(size_t g = 0;g < N;++g)
    {
        if(paternal_gammas[g]+maternal_gammas[g] != 0)
        {
            paternal_gammas[g] += (paternal_round_err/num_good);
			maternal_gammas[g] += (maternal_round_err/num_good);
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

double compute_doc_allele(int bundle_origin, 
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
		double paternal_mass,maternal_mass;
		(**itr).parental_masses(paternal_mass,maternal_mass);
		hits.back().paternal_fpkm(paternal_mass);
		hits.back().maternal_fpkm(maternal_mass);
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
					depth_of_coverage[K - bundle_origin] += (hits[i].paternal_fpkm()+hits[i].maternal_fpkm());
				}
			}
			else if (op.opcode == CUFF_INTRON)
			{
				for (int K = op.g_left(); K < op.g_right(); ++K)
				{
					(*intronic_cov)[K - bundle_origin] += (hits[i].paternal_fpkm()+hits[i].maternal_fpkm());
					//intronic[K - bundle_origin] = true;
				}
				
				pair<map<pair<int,int>,float>::iterator, bool> is = intron_depth_of_coverage.insert(make_pair(make_pair(op.g_left(), op.g_right()), 0));
				is.first->second += (hits[i].paternal_fpkm()+hits[i].maternal_fpkm());
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

void tss_analysis(const string& locus_tag, SampleAlleleAbundances& sample)
{
    // Cluster transcripts by start site (TSS)
    vector<AlleleAbundanceGroup> transcripts_by_tss;
    
    ublas::matrix<double> tss_parental_gamma_cov;
    ublas::matrix<double> tss_parental_count_cov;
    ublas::matrix<double> tss_iterated_exp_parental_count_cov;
    ublas::matrix<double> tss_parental_fpkm_cov;
    vector<Eigen::VectorXd> tss_paternal_assigned_counts;
	vector<Eigen::VectorXd> tss_maternal_assigned_counts;
    
    vector<bool> mask(sample.transcripts.abundances().size(), true);
    for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
    {
        if (*(sample.transcripts.abundances()[i]->tss_id().begin()) == "")
        {
            mask[i] = false;
        }
    }
    
    AlleleAbundanceGroup trans_with_tss;
    sample.transcripts.filter_group(mask, trans_with_tss);
    
    cluster_transcripts<ConnectByAnnotatedTssId>(trans_with_tss,
                                                 transcripts_by_tss,
                                                 &tss_parental_gamma_cov,
                                                 &tss_iterated_exp_parental_count_cov,
                                                 &tss_parental_count_cov,
                                                 &tss_parental_fpkm_cov);
    
    
    BOOST_FOREACH(AlleleAbundanceGroup& ab_group, transcripts_by_tss)
    {
        ab_group.locus_tag(locus_tag);
        set<string> tss_ids = ab_group.tss_id();
        assert (tss_ids.size() == 1);
        string desc = *(tss_ids.begin());
        assert (desc != "");
        ab_group.description(*(tss_ids.begin()));
    }
    
    sample.primary_transcripts = transcripts_by_tss;
    
    // Group TSS clusters by gene
    vector<shared_ptr<Abundance> > primary_transcript_abundances;
    set<shared_ptr<ReadGroupProperties const> > rg_props;
    BOOST_FOREACH (AlleleAbundanceGroup& ab_group, sample.primary_transcripts)
    {
        primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AlleleAbundanceGroup(ab_group)));
        rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
    }
    
    AlleleAbundanceGroup primary_transcripts(primary_transcript_abundances,
                                       tss_parental_gamma_cov,
                                       tss_iterated_exp_parental_count_cov,
                                       tss_parental_count_cov,
                                       tss_parental_fpkm_cov,
                                       rg_props);
    
    vector<AlleleAbundanceGroup> primary_transcripts_by_gene;
    
    cluster_transcripts<ConnectByAnnotatedGeneId>(primary_transcripts,
                                                  primary_transcripts_by_gene);
    
    BOOST_FOREACH(AlleleAbundanceGroup& ab_group, primary_transcripts_by_gene)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        if (gene_ids.size() > 1)
        {
            BOOST_FOREACH (string st, gene_ids)
            {
                fprintf(stderr, "%s\n", st.c_str());
            }
            ab_group.gene_id();
        }
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
    
    sample.gene_primary_transcripts = primary_transcripts_by_gene;
}

void cds_analyis(const string& locus_tag, SampleAlleleAbundances& sample)
{
    // Cluster transcripts by CDS
    vector<AlleleAbundanceGroup> transcripts_by_cds;
    ublas::matrix<double> cds_parental_gamma_cov;
    ublas::matrix<double> cds_parental_count_cov;
    ublas::matrix<double> cds_iterated_exp_parental_count_cov;
    ublas::matrix<double> cds_parental_fpkm_cov;
    
    vector<bool> mask(sample.transcripts.abundances().size(), true);
    for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
    {
        if (*(sample.transcripts.abundances()[i]->protein_id().begin()) == "")
        {
            mask[i] = false;
        }
    }
    
    AlleleAbundanceGroup trans_with_p_id;
    sample.transcripts.filter_group(mask, trans_with_p_id);
    
    cluster_transcripts<ConnectByAnnotatedProteinId>(trans_with_p_id,
                                                     transcripts_by_cds,
                                                     &cds_parental_gamma_cov,
                                                     &cds_iterated_exp_parental_count_cov,
                                                     &cds_parental_count_cov,
                                                     &cds_parental_fpkm_cov);
    
    BOOST_FOREACH(AlleleAbundanceGroup& ab_group, transcripts_by_cds)
    {
        ab_group.locus_tag(locus_tag);
        set<string> protein_ids = ab_group.protein_id();
        assert (protein_ids.size() == 1);
        string desc = *(protein_ids.begin());
        //if (desc != "")
        //{
        assert (desc != "");
        ab_group.description(*(protein_ids.begin()));
        //}
    }
    
    sample.cds = transcripts_by_cds;
    
    // Group the CDS clusters by gene
    vector<shared_ptr<Abundance> > cds_abundances;
    
    set<shared_ptr<ReadGroupProperties const> > rg_props;
    BOOST_FOREACH (AlleleAbundanceGroup& ab_group, sample.cds)
    {
        //if (ab_group.description() != "")
        {
            cds_abundances.push_back(shared_ptr<Abundance>(new AlleleAbundanceGroup(ab_group)));
            rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
        }
    }
    AlleleAbundanceGroup cds(cds_abundances,
                       cds_parental_gamma_cov,
                       cds_iterated_exp_parental_count_cov,
                       cds_parental_count_cov,
                       cds_parental_fpkm_cov,
                       rg_props);
    
    vector<AlleleAbundanceGroup> cds_by_gene;
    
    cluster_transcripts<ConnectByAnnotatedGeneId>(cds,
                                                  cds_by_gene);
    
    BOOST_FOREACH(AlleleAbundanceGroup& ab_group, cds_by_gene)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
    
    sample.gene_cds = cds_by_gene;

}

void sample_abundance_worker(const string& locus_tag,
                             const set<shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAlleleAbundances& sample,
                             shared_ptr<HitBundle> sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis)
{
    vector<shared_ptr<Abundance> > abundances;
    
    BOOST_FOREACH(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        AlleleTranscriptAbundance* pT = new AlleleTranscriptAbundance;
        pT->transfrag(s);
		shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }
    
    sample.transcripts = AlleleAbundanceGroup(abundances);
    
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
        get_alignments_from_scaffolds(sample.transcripts.abundances(),
                                      hits_in_cluster);
        
        // Compute the individual transcript FPKMs via each sample's
        // AlleleAbundanceGroup for this locus.
        
        sample.transcripts.calculate_abundance(hits_in_cluster);
    }
    else
    {
        BOOST_FOREACH(shared_ptr<Abundance>  ab, abundances)
        {
            ab->paternal_status(NUMERIC_HI_DATA);
			ab->maternal_status(NUMERIC_HI_DATA);
            
            CountPerReplicateTable paternal_cpr,maternal_cpr;
            FPKMPerReplicateTable paternal_fpr,maternal_fpr;
            StatusPerReplicateTable paternal_spr,maternal_spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                paternal_cpr[*itr] = 0;
				maternal_cpr[*itr] = 0;
                paternal_fpr[*itr] = 0;
				maternal_fpr[*itr] = 0;
                paternal_spr[*itr] = NUMERIC_HI_DATA;
				maternal_spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_paternal_fragments_by_replicate(paternal_cpr);
			ab->num_maternal_fragments_by_replicate(maternal_cpr);
            ab->paternal_FPKM_by_replicate(paternal_fpr);
			ab->maternal_FPKM_by_replicate(maternal_fpr);
            ab->paternal_status_by_replicate(paternal_spr);
			ab->maternal_status_by_replicate(maternal_spr);
        }
    }
	sample.transcripts.set_allele_informative();
	// Cluster transcripts by gene_id
    vector<AlleleAbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	BOOST_FOREACH(AlleleAbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
    
    if (perform_cds_analysis)
    {
        cds_analyis(locus_tag, sample);
    }
    
    if (perform_tss_analysis)
    {
        tss_analysis(locus_tag, sample);
    }
}

// This function applies library size factors to pre-computed expression entries
void AlleleAbundanceGroup::apply_normalization_to_abundances(const map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AlleleAbundanceGroup> >& unnormalized_ab_group_per_replicate,
                                                       map<shared_ptr<const ReadGroupProperties>, shared_ptr<AlleleAbundanceGroup> >& normalized_ab_group_per_replicate)
{
    for (map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AlleleAbundanceGroup> >::const_iterator itr = unnormalized_ab_group_per_replicate.begin();
         itr != unnormalized_ab_group_per_replicate.end(); ++itr)
    {
        shared_ptr<AlleleAbundanceGroup> norm_ab = shared_ptr<AlleleAbundanceGroup>(new AlleleAbundanceGroup(*itr->second));
        shared_ptr<const ReadGroupProperties> rg_props = itr->first;
        shared_ptr<const MassDispersionModel> disp_model = rg_props->mass_dispersion_model();
        
        shared_ptr<const ReadGroupProperties> old_rg_props = *(itr->second->rg_props().begin());
        
        double fpkm_correction_factor = old_rg_props->normalized_map_mass() / rg_props->normalized_map_mass();
        double internal_scale_factor = rg_props->internal_scale_factor();
        
        double total_mass = 0.0;
        
        for (size_t i = 0; i < norm_ab->_abundances.size(); ++i)
        {
            norm_ab->_abundances[i]->num_paternal_fragments(itr->second->_abundances[i]->num_paternal_fragments() / internal_scale_factor);
			norm_ab->_abundances[i]->num_maternal_fragments(itr->second->_abundances[i]->num_maternal_fragments() / internal_scale_factor);
            
            total_mass += norm_ab->_abundances[i]->num_paternal_fragments()+norm_ab->_abundances[i]->num_maternal_fragments();
            
            norm_ab->_abundances[i]->paternal_FPKM(fpkm_correction_factor * itr->second->_abundances[i]->paternal_FPKM() / internal_scale_factor);
			norm_ab->_abundances[i]->maternal_FPKM(fpkm_correction_factor * itr->second->_abundances[i]->maternal_FPKM() / internal_scale_factor);
            norm_ab->_iterated_exp_parental_count_covariance = norm_ab->iterated_count_cov() / (internal_scale_factor*internal_scale_factor);
            norm_ab->_parental_fpkm_covariance = norm_ab->_parental_fpkm_covariance * (fpkm_correction_factor * fpkm_correction_factor)/ (internal_scale_factor*internal_scale_factor);
            norm_ab->_parental_count_covariance = norm_ab->_parental_count_covariance/ (internal_scale_factor*internal_scale_factor);
        }
        
        double locus_mass_variance = disp_model->scale_mass_variance(total_mass);
        
        for (size_t i = 0; i < norm_ab->_abundances.size(); ++i)
        {
            norm_ab->_abundances[i]->paternal_mass_variance(locus_mass_variance * norm_ab->_abundances[i]->paternal_gamma());
			norm_ab->_abundances[i]->maternal_mass_variance(locus_mass_variance * norm_ab->_abundances[i]->maternal_gamma());
        }
        
        normalized_ab_group_per_replicate[itr->first] = norm_ab;
    }
}

void merge_precomputed_allele_expression_worker(const string& locus_tag,
                                         const vector<shared_ptr<PrecomputedAlleleExpressionBundleFactory> >& expression_factories,
                                         SampleAlleleAbundances& sample,
                                         shared_ptr<HitBundle> sample_bundle,
                                         bool perform_cds_analysis,
                                         bool perform_tss_analysis)
{
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AlleleAbundanceGroup> > unnormalized_ab_group_per_replicate;
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<AlleleAbundanceGroup> > normalized_ab_group_per_replicate;
    map<shared_ptr<const ReadGroupProperties>, shared_ptr<const AlleleAbundanceGroup> > const_ab_group_per_replicate;
    //?
    if (locus_tag == "chr6:31126318-31148508")
    {
        int a = 4;
    }
    
    set<shared_ptr<const ReadGroupProperties> > rg_props;
    for (size_t i = 0; i < expression_factories.size(); ++i)
    {
        shared_ptr<PrecomputedAlleleExpressionBundleFactory> pBundleFac = expression_factories[i];
        shared_ptr<const PrecomputedAlleleExpressionHitFactory> pHitFac = dynamic_pointer_cast<const PrecomputedAlleleExpressionHitFactory> (pBundleFac->hit_factory());
        assert (pHitFac);
        
        shared_ptr<const ReadGroupProperties> rg_prop = pBundleFac->read_group_properties();
        rg_props.insert(rg_prop);
        
        shared_ptr<const AlleleAbundanceGroup> ab = pBundleFac->get_abundance_for_locus(sample_bundle->id());
        pBundleFac->clear_abundance_for_locus(sample_bundle->id());
        if (!ab)
        {
            fprintf(stderr, "Error: no bundle with id %d in precomputed expression file\n", sample_bundle->id());
        }
        else if(ab->abundances().size() != sample_bundle->ref_scaffolds().size())
        {
            fprintf(stderr, "Error: bad bundle merge %s != %s\n", ab->description().c_str(), locus_tag.c_str());
        }
        unnormalized_ab_group_per_replicate[rg_prop] = ab;
    }
    
    AlleleAbundanceGroup::apply_normalization_to_abundances(unnormalized_ab_group_per_replicate, normalized_ab_group_per_replicate);
    
    for (map<shared_ptr<const ReadGroupProperties>, shared_ptr<AlleleAbundanceGroup> >::const_iterator itr = normalized_ab_group_per_replicate.begin();
         itr != normalized_ab_group_per_replicate.end(); ++itr)
    {
        const_ab_group_per_replicate[itr->first] = itr->second;
    }
    
    vector<shared_ptr<Abundance> > abundances;
    
    BOOST_FOREACH(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        AlleleTranscriptAbundance* pT = new AlleleTranscriptAbundance;
        pT->transfrag(s);
		shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
	}

    sample.transcripts = AlleleAbundanceGroup(abundances);
    
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
        sample.transcripts.collect_per_replicate_mass(const_ab_group_per_replicate);
        sample.transcripts.aggregate_replicate_abundances(const_ab_group_per_replicate);
        sample.transcripts.calculate_abundance_group_variance(abundances, const_ab_group_per_replicate);
		sample.transcripts.set_allele_informative(const_ab_group_per_replicate);
	}
    else
    {
        BOOST_FOREACH(shared_ptr<Abundance>  ab, abundances)
        {
            ab->paternal_status(NUMERIC_HI_DATA);
			ab->maternal_status(NUMERIC_HI_DATA);
            
            CountPerReplicateTable paternal_cpr,maternal_cpr;
            FPKMPerReplicateTable paternal_fpr,maternal_fpr;
            StatusPerReplicateTable paternal_spr,maternal_spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                paternal_cpr[*itr] = 0;
				maternal_cpr[*itr] = 0;
                paternal_fpr[*itr] = 0;
				maternal_fpr[*itr] = 0;
                paternal_spr[*itr] = NUMERIC_HI_DATA;
				maternal_spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_paternal_fragments_by_replicate(paternal_cpr);
			ab->num_maternal_fragments_by_replicate(maternal_cpr);
            ab->paternal_FPKM_by_replicate(paternal_fpr);
			ab->maternal_FPKM_by_replicate(maternal_fpr);
            ab->paternal_status_by_replicate(paternal_spr);
			ab->maternal_status_by_replicate(maternal_spr);
		}
    }
	
    // Cluster transcripts by gene_id
    vector<AlleleAbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	BOOST_FOREACH(AlleleAbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
	
    if (perform_cds_analysis)
    {
        cds_analyis(locus_tag, sample);
    }
    
    if (perform_tss_analysis)
    {
        tss_analysis(locus_tag, sample);
    }
}

//BOOST_CLASS_EXPORT(Abundance)
BOOST_CLASS_EXPORT(TranscriptAbundance)
BOOST_SERIALIZATION_SHARED_PTR(TranscriptAbundance);
BOOST_CLASS_EXPORT(AlleleTranscriptAbundance)
BOOST_SERIALIZATION_SHARED_PTR(AlleleTranscriptAbundance);
BOOST_CLASS_EXPORT(AbundanceGroup);
BOOST_SERIALIZATION_SHARED_PTR(AbundanceGroup);
BOOST_CLASS_EXPORT(AlleleAbundanceGroup);
BOOST_SERIALIZATION_SHARED_PTR(AlleleAbundanceGroup);
BOOST_SERIALIZATION_ASSUME_ABSTRACT(Abundance);

