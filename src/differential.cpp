/*
 *  differential.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/15/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <algorithm>
#include <functional>
#include <numeric>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include "abundances.h"
#include "differential.h"
#include "clustering.h"
#include "differential.h"
#include "sampling.h"

using namespace std;

double min_read_count = 10;

#if ENABLE_THREADS
mutex _launcher_lock;
mutex locus_thread_pool_lock;
int locus_curr_threads = 0;
int locus_num_threads = 0;

void decr_pool_count()
{
	locus_thread_pool_lock.lock();
	locus_curr_threads--;
	locus_thread_pool_lock.unlock();	
}
#endif

SampleDifference::SampleDifference():
    sample_1(-1),
    sample_2(-1),
    value_1(0.0),
    value_2(0.0),
    test_stat(0.0),
    p_value(1.0),
    corrected_p(1.0),
    test_status(NOTEST),
    significant(false) {}


void add_to_tracking_table(size_t sample_index,
                           Abundance& ab,
						   FPKMTrackingTable& track)

{
	pair<FPKMTrackingTable::iterator,bool> inserted;
	pair<string, FPKMTracking > p;
	p = make_pair(ab.description(), FPKMTracking());
	inserted = track.insert(p);
	
	FPKMTracking& fpkm_track = inserted.first->second;
	
	set<string> tss = ab.tss_id();
    set<string> gene_ids = ab.gene_id();
	set<string> genes = ab.gene_name();
	set<string> proteins = ab.protein_id();
	
	fpkm_track.tss_ids.insert(tss.begin(), tss.end());
    fpkm_track.gene_ids.insert(gene_ids.begin(), gene_ids.end());
	fpkm_track.gene_names.insert(genes.begin(), genes.end());
	fpkm_track.protein_ids.insert(proteins.begin(), proteins.end());
	
	if (inserted.second)
	{
		fpkm_track.locus_tag = ab.locus_tag();
		fpkm_track.description = ab.description();
		shared_ptr<Scaffold> transfrag = ab.transfrag();
		if (transfrag && transfrag->nearest_ref_id() != "")
		{
			fpkm_track.classcode = transfrag->nearest_ref_classcode();
			fpkm_track.ref_match = transfrag->nearest_ref_id();
		}
		else
		{
			fpkm_track.classcode = 0;
			fpkm_track.ref_match = "-";
		}
        if (transfrag)
        {
            fpkm_track.length = transfrag->length(); 
        }
        else
        {
            fpkm_track.length = 0;
        }
	}
	
	FPKMContext r1 = FPKMContext(ab.num_fragments(), 
                                 ab.num_fragment_var(),
                                 ab.num_fragment_uncertainty_var(),
                                 ab.mass_variance(),
                                 ab.num_fragments_by_replicate(),
								 ab.FPKM(), 
								 ab.FPKM_by_replicate(),
                                 ab.FPKM_variance(),
                                 ab.FPKM_conf().low,
                                 ab.FPKM_conf().high,
                                 ab.status(),
                                 ab.status_by_replicate(),
                                 ab.fpkm_samples(),
                                 ab.gamma());
    
    
	
    vector<FPKMContext>& fpkms = inserted.first->second.fpkm_series;
    if (sample_index < fpkms.size())
    {
        // if the fpkm series already has an entry matching this description
        // for this sample index, then we are dealing with a group of transcripts
        // that occupies multiple (genomically disjoint) bundles.  We need
        // to add this bundle's contribution to the FPKM, fragments, and variance 
        // to whatever's already there.  
        
        // NOTE: we can simply sum the FKPM_variances, because we are currently
        // assuming that transcripts in disjoint bundles share no alignments and 
        // thus have FPKM covariance == 0;  This assumption will no longer be
        // true if we decide to do multireads the right way.
        
        FPKMContext& existing = fpkms[sample_index];
        existing.FPKM += r1.FPKM;
        existing.count_mean += r1.count_mean;
        existing.FPKM_variance += r1.FPKM_variance;
        if (existing.status == NUMERIC_FAIL || r1.status == NUMERIC_FAIL)
        {
            existing.status = NUMERIC_FAIL;
        }
        else 
        {
            existing.status = NUMERIC_OK;
        }
        
    }
    else 
    {
        fpkms.push_back(r1);
    }
}


TestLauncher::launcher_sample_table::iterator TestLauncher::find_locus(const string& locus_id)
{
    launcher_sample_table::iterator itr = _samples.begin();
    for(; itr != _samples.end(); ++itr)
    {
        if (itr->first == locus_id)
            return itr;
    }
    return _samples.end();
}

void TestLauncher::register_locus(const string& locus_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        pair<launcher_sample_table::iterator, bool> p;
        vector<shared_ptr<SampleAbundances> >abs(_orig_workers);
        _samples.push_back(make_pair(locus_id, abs));
    }
}

void TestLauncher::abundance_avail(const string& locus_id, 
                                   shared_ptr<SampleAbundances> ab, 
                                   size_t factory_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        assert(false);
    }
    itr->second[factory_id] = ab;
    //itr->second(factory_id] = ab;
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
bool TestLauncher::all_samples_reported_in(vector<shared_ptr<SampleAbundances> >& abundances)
{    
    BOOST_FOREACH (shared_ptr<SampleAbundances> ab, abundances)
    {
        if (!ab)
        {
            return false;
        }
    }
    return true;
}

#if ENABLE_THREADS
mutex test_storage_lock; // don't modify the above struct without locking here
#endif

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void TestLauncher::perform_testing(vector<shared_ptr<SampleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
    
    // Just verify that all the loci from each factory match up.
    for (size_t i = 1; i < abundances.size(); ++i)
    {
        const SampleAbundances& curr = *(abundances[i]);
        const SampleAbundances& prev = *(abundances[i-1]);
        
        assert (curr.locus_tag == prev.locus_tag);
        
        const AbundanceGroup& s1 = curr.transcripts;
        const AbundanceGroup& s2 =  prev.transcripts;
        
        assert (s1.abundances().size() == s2.abundances().size());
        
        for (size_t j = 0; j < s1.abundances().size(); ++j)
        {
            assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
        }
    }
    
    test_differential(abundances.front()->locus_tag, abundances, _contrasts, *_tests, *_tracking);
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void TestLauncher::record_tracking_data(vector<shared_ptr<SampleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
    
    // Just verify that all the loci from each factory match up.
    for (size_t i = 1; i < abundances.size(); ++i)
    {
        const SampleAbundances& curr = *(abundances[i]);
        const SampleAbundances& prev = *(abundances[i-1]);
        
        assert (curr.locus_tag == prev.locus_tag);
        
        const AbundanceGroup& s1 = curr.transcripts;
        const AbundanceGroup& s2 =  prev.transcripts;
        
        assert (s1.abundances().size() == s2.abundances().size());
        
        for (size_t j = 0; j < s1.abundances().size(); ++j)
        {
            assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
        }
    }
    
#if ENABLE_THREADS
	test_storage_lock.lock();
#endif
    
    // Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		const AbundanceGroup& ab_group = abundances[i]->transcripts;
        //fprintf(stderr, "[%d] count = %lg\n",i,  ab_group.num_fragments());
		BOOST_FOREACH (shared_ptr<Abundance> ab, ab_group.abundances())
		{
			add_to_tracking_table(i, *ab, _tracking->isoform_fpkm_tracking);
            //assert (_tracking->isoform_fpkm_tracking.num_fragments_by_replicate().empty() == false);
		}
		
		BOOST_FOREACH (AbundanceGroup& ab, abundances[i]->cds)
		{
			add_to_tracking_table(i, ab, _tracking->cds_fpkm_tracking);
		}
		
		BOOST_FOREACH (AbundanceGroup& ab, abundances[i]->primary_transcripts)
		{
			add_to_tracking_table(i, ab, _tracking->tss_group_fpkm_tracking);
		}
		
		BOOST_FOREACH (AbundanceGroup& ab, abundances[i]->genes)
		{
			add_to_tracking_table(i, ab, _tracking->gene_fpkm_tracking);
		}
	}
    
#if ENABLE_THREADS
    test_storage_lock.unlock();
#endif
    
}

void TestLauncher::test_finished_loci()
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif  

    launcher_sample_table::iterator itr = _samples.begin(); 
    while(itr != _samples.end())
    {
        if (all_samples_reported_in(itr->second))
        {
            // In some abundance runs, we don't actually want to perform testing 
            // (eg initial quantification before bias correction).
            // _tests and _tracking will be NULL in these cases.
            if (_tests != NULL && _tracking != NULL)
            {
                if (_p_bar)
                {
                    verbose_msg("Testing for differential expression and regulation in locus [%s]\n", itr->second.front()->locus_tag.c_str());
                    _p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
                }
                record_tracking_data(itr->second);
                perform_testing(itr->second);
                
            }
            else
            {
                if (_p_bar)
                {
                    //verbose_msg("Testing for differential expression and regulation in locus [%s]\n", abundances.front()->locus_tag.c_str());
                    _p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
                }
                record_tracking_data(itr->second);
            }
            
            // Removes the samples that have already been tested and transferred to the tracking tables,
            itr = _samples.erase(itr);
        }
        else
        {
            
            ++itr;
        }
    }
}

long double wrapped_lgamma(double gamma_k)
{
    long double lgamma_k = 1;
    try
    {
        lgamma_k = boost::math::lgamma<long double>(gamma_k);
    }
    catch (boost::exception_detail::error_info_injector<std::overflow_error>& e)
    {
        lgamma_k = numeric_limits<long double>::max();
    }
    catch (boost::exception_detail::error_info_injector<std::underflow_error>& e)
    {
        lgamma_k = numeric_limits<long double>::max();
    }
    catch (boost::exception_detail::error_info_injector<std::domain_error>& e)
    {
        lgamma_k = numeric_limits<long double>::max();
    }
    return lgamma_k;
}

struct TruncNormalMoments
{
    double m1;
    double m2;
    double m3;
};

void calculate_sample_moments(const vector<double>& samples, TruncNormalMoments& moments)
{
    double m1 = accumulate(samples.begin(), samples.end(), 0.0);
    if (samples.size())
        m1 /= samples.size();

    double m2 = 0.0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        m2 += samples[i]*samples[i];
    }
    if (samples.size())
        m2 /= samples.size();
    
    double m3 = 0.0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        m3 += samples[i]*samples[i]*samples[i];
    }
    if (samples.size())
        m3 /= samples.size();
    
    moments.m1 = m1;
    moments.m2 = m2;
    moments.m3 = m3;
}

// Estimate the location and scale parameters for the normal distribution
// that's truncated at zero.
void estimate_trunc_normal_params(const TruncNormalMoments& m, double& location, double& scale)
{
    double denom = 2.0*m.m1*m.m1 - m.m2;
    if (denom > 0)
    {
        location = (2.0*m.m1*m.m2 - m.m3)/denom;
        scale = (m.m1*m.m3 - m.m2*m.m2)/denom;
    }
    else
    {
        location = 0;
        scale = 0;
    }
    
}

long double trunc_normal_log_likelihood(const vector<double>& samples, double mean, double variance)
{
    boost::math::normal norm; //standard normal
    
    if (mean < 0 || variance <= 0.0)
        return 0;
    
    double denom_sigma = 1.0/sqrt(variance);
    
    long double log_likelihood = 0.0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        long double sL = 0.0;
        sL = denom_sigma * pdf(norm, (samples[i] - mean) * denom_sigma);
        double sL_denom = (1.0 - cdf(norm, (0.0 - mean) * denom_sigma));
        if (sL_denom != 0)
            sL /= sL_denom;
        else
            sL = 0.0;
        
        sL = logl(sL);
        log_likelihood += sL;
    }
    return log_likelihood;
}

//// Sampling-based test:
SampleDifference test_diffexp(const FPKMContext& curr,
                              const FPKMContext& prev)
{
	bool performed_test = false;
    
    SampleDifference test = SampleDifference();
    
    test.p_value = 1.0;
    test.differential = 0.0;
    test.test_stat = 0.0;
    test.test_status = NOTEST;
    test.value_1 = 0;
    test.value_2 = 0;
    test.significant = 0;
    test.corrected_p = 1.0;
    
    double p_value = 1.0;
        
    double differential = 0.0;
    
    // static const long double min_gamma_params = 1e-20;
    
    vector<double> null_log_ratio_samples;
    
    static const size_t num_null_ratio_samples = 10000;
    
    boost::random::mt19937 rng;
    
    if ((curr.FPKM != 0 || prev.FPKM != 0) && (prev.fpkm_samples.size() > 0 && curr.fpkm_samples.size() > 0))
    {
        boost::random::uniform_int_distribution<> prev_sampler(0, prev.fpkm_samples.size()-1);
        boost::random::uniform_int_distribution<> curr_sampler(0, curr.fpkm_samples.size()-1);
    
        vector<double> prev_rep_samples;
        vector<double> curr_rep_samples;
        
        
        for (FPKMPerReplicateTable::const_iterator itr = curr.fpkm_per_rep.begin();
             itr != curr.fpkm_per_rep.end(); ++itr)
        {
            StatusPerReplicateTable::const_iterator si = curr.status_per_rep.find(itr->first);
            if (si == curr.status_per_rep.end() || si->second == NUMERIC_LOW_DATA)
                continue;
            curr_rep_samples.push_back(itr->second);
        }
        
        for (FPKMPerReplicateTable::const_iterator itr = prev.fpkm_per_rep.begin();
             itr != prev.fpkm_per_rep.end(); ++itr)
        {
            StatusPerReplicateTable::const_iterator si = prev.status_per_rep.find(itr->first);
            if (si == prev.status_per_rep.end() || si->second == NUMERIC_LOW_DATA)
                continue;
            prev_rep_samples.push_back(itr->second);
        }
        
        double curr_fpkm = accumulate(curr_rep_samples.begin(), curr_rep_samples.end(), 0.0);
        if (curr_rep_samples.size() > 0)
            curr_fpkm /= curr_rep_samples.size();
        
        double prev_fpkm = accumulate(prev_rep_samples.begin(), prev_rep_samples.end(), 0.0);
        if (prev_rep_samples.size() > 0)
            prev_fpkm /= prev_rep_samples.size();
        
        
        if (curr_fpkm > 0.0 && prev_fpkm > 0.0)
            differential = log2(curr_fpkm) - log2(prev_fpkm);
        else if (curr_fpkm)
            differential = numeric_limits<double>::infinity();
        else if (prev_fpkm)
            differential = -numeric_limits<double>::infinity();

        // set the asymptotic delta method test stat for backward compatibility
        double curr_log_fpkm_var = (curr.FPKM_variance) / (curr.FPKM * curr.FPKM);
        double prev_log_fpkm_var = (prev.FPKM_variance) / (prev.FPKM * prev.FPKM);
        double numerator = log(curr.FPKM / prev.FPKM);
        double denominator = sqrt(prev_log_fpkm_var + curr_log_fpkm_var);
        if (denominator > 0.0)
            test.test_stat = numerator / denominator;
        else if (numerator > 0)
            test.test_stat = numeric_limits<double>::infinity();
        else if (numerator < 0)
            test.test_stat = -numeric_limits<double>::infinity();
        else
            test.test_stat = 0;

        
        test.test_stat = numerator / denominator;
        
        
        // Draw from prev fpkm_samples to make the first half of the null
        for (size_t i = 0; i < num_null_ratio_samples; ++i)
        {
            double curr_set_sample = 0.0;
            for (size_t k = 0; k < curr_rep_samples.size(); ++k)
            {
                int next_sample_idx = prev_sampler(rng);
                if (next_sample_idx >= 0 && next_sample_idx < prev.fpkm_samples.size())
                    curr_set_sample += prev.fpkm_samples[next_sample_idx] / (double)curr_rep_samples.size();
            }
            
            double prev_set_sample = 0.0;
            for (size_t k = 0; k < prev_rep_samples.size(); ++k)
            {
                int next_sample_idx = prev_sampler(rng);
                if (next_sample_idx >= 0 && next_sample_idx < prev.fpkm_samples.size())
                    prev_set_sample += prev.fpkm_samples[next_sample_idx] / (double)prev_rep_samples.size();
            }
            
            double null_ratio_sample = 0.0;
            if (curr_set_sample > 0.0 && prev_set_sample > 0.0)
                null_ratio_sample = log2(curr_set_sample) - log2(prev_set_sample);
            else if (curr_set_sample > 0.0)
                null_ratio_sample = numeric_limits<double>::infinity();
            else if (prev_set_sample)
                null_ratio_sample = -numeric_limits<double>::infinity();
            
            null_log_ratio_samples.push_back(null_ratio_sample);
        }
        
        // Draw from curr's fpkm_samples to make the other half of the null
        for (size_t i = 0; i < num_null_ratio_samples; ++i)
        {
            double curr_set_sample = 0.0;
            for (size_t k = 0; k < curr_rep_samples.size(); ++k)
            {
                int next_sample_idx = curr_sampler(rng);
                if (next_sample_idx >= 0 && next_sample_idx < curr.fpkm_samples.size())
                    curr_set_sample += curr.fpkm_samples[next_sample_idx] / (double)curr_rep_samples.size();
            }
            
            double prev_set_sample = 0.0;
            for (size_t k = 0; k < prev_rep_samples.size(); ++k)
            {
                int next_sample_idx = curr_sampler(rng);
                if (next_sample_idx >= 0 && next_sample_idx < curr.fpkm_samples.size())
                    prev_set_sample += curr.fpkm_samples[next_sample_idx] / (double)prev_rep_samples.size();
            }
            
            double null_ratio_sample = 0.0;
            if (curr_set_sample > 0.0 && prev_set_sample > 0.0)
                null_ratio_sample = log2(curr_set_sample) - log2(prev_set_sample);
            else if (curr_set_sample > 0.0)
                null_ratio_sample = numeric_limits<double>::infinity();
            else if (prev_set_sample)
                null_ratio_sample = -numeric_limits<double>::infinity();
            
            null_log_ratio_samples.push_back(null_ratio_sample);
        }

        
        sort(null_log_ratio_samples.begin(), null_log_ratio_samples.end());
        
        double lower_tail_val;
        double upper_tail_val;
        if (differential > 0)
        {
            upper_tail_val = differential;
            lower_tail_val = -differential;
        }
        else
        {
            upper_tail_val = -differential;
            lower_tail_val = differential;
        }
        
        vector<double>::iterator lower_tail_null_range_iter = upper_bound(null_log_ratio_samples.begin(), null_log_ratio_samples.end(), lower_tail_val);
        vector<double>::iterator upper_tail_null_range_iter = upper_bound(null_log_ratio_samples.begin(), null_log_ratio_samples.end(), upper_tail_val);
        size_t num_smaller_lower_tail = lower_tail_null_range_iter - null_log_ratio_samples.begin();
        size_t num_smaller_upper_tail = upper_tail_null_range_iter - null_log_ratio_samples.begin();
        
        size_t num_samples_more_extreme = (null_log_ratio_samples.size() - num_smaller_upper_tail) + num_smaller_lower_tail;
        
        p_value = num_samples_more_extreme / (double)null_log_ratio_samples.size();
        if (p_value == 0)
            p_value = 1.0 / (double)null_log_ratio_samples.size();
        if (p_value > 1)
            p_value = 1.0;
        
        
        
        performed_test = true;

        //test = SampleDifference(sample1, sample2, prev.FPKM, curr.FPKM, stat, p_value, transcript_group_id);
        test.p_value = p_value;
        test.differential = differential;
        //test.test_stat = stat;
        test.value_1 = prev.FPKM;
        test.value_2 = curr.FPKM;
    }
    else
    {
        test.p_value = 1.0;
        test.test_stat = 0.0;
        test.value_1 = prev.FPKM;
        test.value_2 = curr.FPKM;
        test.differential = 0;
        performed_test = false;
    }
    
	test.test_status = performed_test ? OK : NOTEST;
	return test;
}

SampleDiffMetaDataTable meta_data_table;
#if ENABLE_THREADS
boost::mutex meta_data_lock;
#endif

shared_ptr<SampleDifferenceMetaData> get_metadata(const string description)
{
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(meta_data_lock);
#endif
    pair<SampleDiffMetaDataTable::iterator, bool> p;
    p = meta_data_table.insert(make_pair(description, new SampleDifferenceMetaData()));
    return p.first->second;
}

// This performs between-group tests on isoforms or TSS groupings in a single
// locus, on two different samples.
SampleDifference get_de_tests(const string& description,
                 const FPKMContext& prev_abundance,
				 const FPKMContext& curr_abundance,
				 //SampleDiffs& de_tests,
				 bool enough_reads)
{
	int total_iso_de_tests = 0;
			
	SampleDifference test;
    
    const FPKMContext& r1 = curr_abundance;
    const FPKMContext& r2 = prev_abundance;
    
	if (curr_abundance.status == NUMERIC_FAIL || 
        prev_abundance.status == NUMERIC_FAIL ||
        prev_abundance.status == NUMERIC_HI_DATA ||
        curr_abundance.status == NUMERIC_HI_DATA)
    {
        test = test_diffexp(r1, r2);
        test.test_stat = 0;
		test.p_value = 1.0;
		test.differential = 0.0;
        if (curr_abundance.status == NUMERIC_FAIL || 
            prev_abundance.status == NUMERIC_FAIL)
        {
            test.test_status = FAIL;
        }
        else if (prev_abundance.status == NUMERIC_HI_DATA ||
                 curr_abundance.status == NUMERIC_HI_DATA)
        {
            test.test_status = HIDATA;
        }
    }
    else if (curr_abundance.status == NUMERIC_LOW_DATA || 
             prev_abundance.status == NUMERIC_LOW_DATA)
    {
        // perform the test, but mark it as not significant and don't add it to the 
        // pile. This way we don't penalize for multiple testing, but users can still
        // see the fold change.
		test = test_diffexp(r1, r2);
        test.test_stat = 0;
        test.p_value = 1.0;
        //test.differential = 0.0;
		
		test.test_status = LOWDATA;
    }
    else // at least one is OK, the other might be LOW_DATA
	{
        test = test_diffexp(r1, r2);
        if (test.test_status == OK && enough_reads)
        {
            total_iso_de_tests++;
        }
        else
        {
            test.test_status = NOTEST;
			test.test_stat = 0;
			test.p_value = 1.0;
			//test.differential = 0.0;
        }
//		if (test_diffexp(r1, r2, test))
//		{
//			total_iso_de_tests++;
//		}
//		else
//		{
//			test.test_stat = 0;
//			test.p_value = 1.0;
//			test.differential = 0.0;
//		}
//		if (enough_reads)
//			test.test_status = OK;
//		else
//			test.test_status = NOTEST;
		
	}
	
    
	//inserted.first->second = test;
	
	//return make_pair(total_iso_de_tests, inserted.first);
    return test;
}

bool test_js(const AbundanceGroup& prev_abundance,
             const AbundanceGroup& curr_abundance,
             double& js,
             double& p_val)
{
    vector<Eigen::VectorXd> sample_kappas;
//    Eigen::VectorXd curr_kappas(Eigen::VectorXd::Zero(curr_abundance.abundances().size()));
//    for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
//    {
//        curr_kappas(i) = curr_abundance.abundances()[i]->kappa();
//    }
//    
//    Eigen::VectorXd prev_kappas(Eigen::VectorXd::Zero(prev_abundance.abundances().size()));
//    for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
//    {
//        prev_kappas(i) = prev_abundance.abundances()[i]->kappa();
//    }
    
    
    ////////////
    
    Eigen::VectorXd prev_kappas = Eigen::VectorXd::Zero(prev_abundance.abundances().size());
    for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
    {
        FPKMPerReplicateTable ab_i_by_rep = prev_abundance.abundances()[i]->FPKM_by_replicate();
        for (FPKMPerReplicateTable::const_iterator itr = ab_i_by_rep.begin(); itr != ab_i_by_rep.end(); ++itr)
        {
            prev_kappas(i) += itr->second;
        }
        prev_kappas(i) /= ab_i_by_rep.size();
    }
    
    if (prev_kappas.sum() > 0)
        prev_kappas /= prev_kappas.sum();

    Eigen::VectorXd curr_kappas = Eigen::VectorXd::Zero(curr_abundance.abundances().size());
    for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
    {
        FPKMPerReplicateTable ab_i_by_rep = curr_abundance.abundances()[i]->FPKM_by_replicate();
        for (FPKMPerReplicateTable::const_iterator itr = ab_i_by_rep.begin(); itr != ab_i_by_rep.end(); ++itr)
        {
            curr_kappas(i) += itr->second;
        }
        curr_kappas(i) /= ab_i_by_rep.size();
    }
    
    if (curr_kappas.sum() > 0)
        curr_kappas /= curr_kappas.sum();
    
    ////////////
    
    sample_kappas.push_back(prev_kappas);
    sample_kappas.push_back(curr_kappas);
    
    js = jensen_shannon_distance(sample_kappas);
    
    if (isinf(js) || isnan(js))
        return false;
    
    vector<double> null_js_samples;
    
    static const size_t num_null_js_samples = 10000;
    
    boost::random::mt19937 rng;
    boost::random::uniform_int_distribution<> prev_sampler(0, prev_abundance.member_fpkm_samples().size()-1);
    boost::random::uniform_int_distribution<> curr_sampler(0, curr_abundance.member_fpkm_samples().size()-1);
    
    for (size_t i = 0; i < num_null_js_samples; ++i)
    {
        // Draw k values, from prev's fpkm_samples to make the first half of the null, where k is the number of replicates in *curr*, and compute their average
        Eigen::VectorXd curr_set_sample = Eigen::VectorXd::Zero(curr_abundance.abundances().size());
        for (size_t k = 0; k < curr_abundance.rg_props().size(); ++k)
        {
            int next_sample_idx = prev_sampler(rng);
            if (next_sample_idx >= 0 && next_sample_idx < prev_abundance.member_fpkm_samples().size())
                curr_set_sample += prev_abundance.member_fpkm_samples()[next_sample_idx] / (double)curr_abundance.rg_props().size();
        }
        double total_curr_set_sample_fpkm = curr_set_sample.sum();
        if (total_curr_set_sample_fpkm > 0.0)
            curr_set_sample /= total_curr_set_sample_fpkm;
        
        // Draw k values, from prev's fpkm_samples to make the first half of the null, where k is the number of replicates in *prev*, and compute their average
        Eigen::VectorXd prev_set_sample = Eigen::VectorXd::Zero(prev_abundance.abundances().size());
        for (size_t k = 0; k < prev_abundance.rg_props().size(); ++k)
        {
            int next_sample_idx = prev_sampler(rng);
            if (next_sample_idx >= 0 && next_sample_idx < prev_abundance.member_fpkm_samples().size())
                prev_set_sample += prev_abundance.member_fpkm_samples()[next_sample_idx] / (double)prev_abundance.rg_props().size();
        }
        double total_prev_set_sample_fpkm = prev_set_sample.sum();
        if (total_prev_set_sample_fpkm > 0.0)
            prev_set_sample /= total_prev_set_sample_fpkm;
        
        vector<Eigen::VectorXd> sample_kappas;
        
        sample_kappas.push_back(curr_set_sample);
        sample_kappas.push_back(prev_set_sample);
        
        double js_sample = jensen_shannon_distance(sample_kappas);
        
        //cerr << curr_set_sample.transpose() << " vs. " << prev_set_sample.transpose() << " : " << js_sample << endl;
        
        null_js_samples.push_back(js_sample);
    }
    
    for (size_t i = 0; i < num_null_js_samples; ++i)
    {
        // Draw k values, from curr's fpkm_samples to make the first half of the null, where k is the number of replicates in *curr*, and compute their average
        Eigen::VectorXd curr_set_sample = Eigen::VectorXd::Zero(curr_abundance.abundances().size());
        for (size_t k = 0; k < curr_abundance.rg_props().size(); ++k)
        {
            int next_sample_idx = curr_sampler(rng);
            if (next_sample_idx >= 0 && next_sample_idx < curr_abundance.member_fpkm_samples().size())
                curr_set_sample += curr_abundance.member_fpkm_samples()[next_sample_idx] / (double)curr_abundance.rg_props().size();
        }
        double total_curr_set_sample_fpkm = curr_set_sample.sum();
        if (total_curr_set_sample_fpkm > 0.0)
            curr_set_sample /= total_curr_set_sample_fpkm;
        
        // Draw k values, from curr's fpkm_samples to make the first half of the null, where k is the number of replicates in *prev*, and compute their average
        Eigen::VectorXd prev_set_sample = Eigen::VectorXd::Zero(curr_abundance.abundances().size());
        for (size_t k = 0; k < prev_abundance.rg_props().size(); ++k)
        {
            int next_sample_idx = curr_sampler(rng);
            if (next_sample_idx >= 0 && next_sample_idx < curr_abundance.member_fpkm_samples().size())
                prev_set_sample += curr_abundance.member_fpkm_samples()[next_sample_idx] / (double)prev_abundance.rg_props().size();
        }
        double total_prev_set_sample_fpkm = prev_set_sample.sum();
        if (total_prev_set_sample_fpkm > 0.0)
            prev_set_sample /= total_prev_set_sample_fpkm;
        
        vector<Eigen::VectorXd> sample_kappas;
        
        sample_kappas.push_back(curr_set_sample);
        sample_kappas.push_back(prev_set_sample);
        
        double js_sample = jensen_shannon_distance(sample_kappas);
        
        null_js_samples.push_back(js_sample);
    }    
    
    sort(null_js_samples.begin(), null_js_samples.end());
    
    double upper_tail_val = js;
    
    vector<double>::iterator upper_tail_null_range_iter = upper_bound(null_js_samples.begin(), null_js_samples.end(), upper_tail_val);
    
    size_t num_smaller_upper_tail = upper_tail_null_range_iter - null_js_samples.begin();
    
    size_t num_samples_more_extreme = (null_js_samples.size() - num_smaller_upper_tail);
    
    p_val = num_samples_more_extreme / (double)null_js_samples.size();
    if (p_val == 0)
        p_val = 1.0 / (double)null_js_samples.size();
    if (p_val > 1)
        p_val = 1.0;
    
    return true;
}


// This performs within-group tests on a set of isoforms or a set of TSS groups.
// This is a way of looking for meaningful differential splicing or differential
// promoter use.
SampleDifference get_ds_tests(const AbundanceGroup& prev_abundance,
                              const AbundanceGroup& curr_abundance,
//                              SampleDiffs& diff_tests,
                              bool enough_reads)
{
	SampleDifference test;
    
	test.p_value = 1.0;
    test.differential = 0.0;
    test.test_stat = 0.0;
    test.test_status = NOTEST;
    test.value_1 = 0;
    test.value_2 = 0;
    test.significant = 0;
    test.corrected_p = 1.0;

	
	AbundanceStatus prev_status = curr_abundance.status();
	AbundanceStatus curr_status = prev_abundance.status();
    
    if (prev_abundance.abundances().size() == 1 ||
        prev_abundance.num_fragments() == 0 ||
        curr_abundance.num_fragments() == 0 ||
        prev_abundance.member_fpkm_samples().size() == 0 ||
        curr_abundance.member_fpkm_samples().size() == 0)
    {
        test.p_value = 1;
        test.value_1 = 0;
        test.value_2 = 0;
        test.differential = 0;
        test.test_status = NOTEST;
    }
	else if (prev_abundance.abundances().size() > 1 &&
        /*prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
        filtered_curr.has_member_with_status(NUMERIC_LOW_DATA) == false &&*/ 
        prev_status == NUMERIC_OK && prev_abundance.num_fragments() > 0 &&
		curr_status == NUMERIC_OK && curr_abundance.num_fragments() > 0)
	{
		vector<ublas::vector<double> > sample_kappas;
		ublas::vector<double> curr_kappas(curr_abundance.abundances().size());
        
        for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->kappa();
		}
		
		ublas::vector<double> prev_kappas(prev_abundance.abundances().size());
        for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->kappa();
		}
		
		sample_kappas.push_back(prev_kappas);
		sample_kappas.push_back(curr_kappas);
		
        double js = 0.0;
        double p_val = 1.0;
        
        bool success;
        success = test_js(prev_abundance, curr_abundance, js, p_val);
        
		if (js == 0.0 || success == false)
		{
			test.test_stat = 0;
			test.p_value = 1.0;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = 0;
			test.test_status = NOTEST;
		}
		else
		{
            //test.test_stat = 0;
			test.p_value = p_val;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = js;
			test.test_status = enough_reads ? OK : NOTEST;
        }
	}
	else // we won't even bother with the JS-based testing in LOWDATA cases.
	{
        if (prev_status == NUMERIC_OK && curr_status == NUMERIC_OK && 
            prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
            curr_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false)
            test.test_status = NOTEST;
        else if (prev_status == NUMERIC_FAIL || curr_status == NUMERIC_FAIL)
            test.test_status = FAIL;
        else
            test.test_status = LOWDATA;
            
		test.test_stat = 0;
		test.p_value = 0.0;
		test.differential = 0.0;
	}
    
    return test;
}

string make_ref_tag(const string& ref, char classcode)
{
	char tag_buf[1024];

	sprintf(tag_buf, 
			"%s(%c)",
			ref.c_str(),
			classcode);
	
	return string(tag_buf);
}

string bundle_locus_tag(const RefSequenceTable& rt, 
						const HitBundle& bundle)
{
	char locus_buf[1024];
	RefID bundle_chr_id = bundle.ref_id();
	assert (bundle_chr_id != 0);
	const char* chr_name = rt.get_name(bundle_chr_id);
	
	sprintf(locus_buf, 
			"%s:%d-%d",
			chr_name,
			bundle.left(),
			bundle.right());
	
	return string(locus_buf);
}

void sample_abundance_worker(const string& locus_tag,
                             const set<shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAbundances& sample,
                             HitBundle* sample_bundle,
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
    
    if (perform_tss_analysis)
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
//            if (gene_ids.size() > 1)
//            {
//                BOOST_FOREACH (string st, gene_ids)
//                {
//                    fprintf(stderr, "%s\n", st.c_str());
//                }
//                ab_group.gene_id();
//            }
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_primary_transcripts = primary_transcripts_by_gene;
    }
}

struct LocusVarianceInfo
{
    int factory_id;
    double count_mean;
    double count_empir_var;
    double locus_count_fitted_var;
    double isoform_fitted_var_sum;
    double cross_replicate_js;
    int num_transcripts;
    double bayes_gamma_trace;
    double empir_gamma_trace;
    vector<double> gamma;
    vector<double> gamma_var;
    vector<double> gamma_bootstrap_var;
    vector<string> transcript_ids;
    vector<double> count_sharing;
    double locus_salient_frags;
    double locus_total_frags;

};

#if ENABLE_THREADS
mutex variance_info_lock; // don't modify the above struct without locking here
#endif

vector<LocusVarianceInfo> locus_variance_info_table;

// We'll use this tracking table to collect per replicate counts for each 
// transcript, so we can re-fit the variance model.
FPKMTrackingTable transcript_count_tracking; 

void sample_worker(const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   shared_ptr<SampleAbundances> abundance,
                   size_t factory_id,
                   shared_ptr<TestLauncher> launcher)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    HitBundle bundle;
    bool non_empty = sample_factory.next_bundle(bundle);
    
    char bundle_label_buf[2048];
    sprintf(bundle_label_buf,
            "%s:%d-%d",
            rt.get_name(bundle.ref_id()),
            bundle.left(),
            bundle.right());
    string locus_tag = bundle_label_buf;
    
    if (!non_empty || (bias_run && bundle.ref_scaffolds().size() != 1)) // Only learn on single isoforms
    {
#if !ENABLE_THREADS
        // If Cuffdiff was built without threads, we need to manually invoke 
        // the testing functor, which will check to see if all the workers
        // are done, and if so, perform the cross sample testing.
        launcher->abundance_avail(locus_tag, abundance, factory_id);
        launcher->test_finished_loci();
        //launcher();
#endif
    	return;
    }
    
    abundance->cluster_mass = bundle.mass();
    
    launcher->register_locus(locus_tag);
    
    abundance->locus_tag = locus_tag;
    
    bool perform_cds_analysis = false;
    bool perform_tss_analysis = false;

    BOOST_FOREACH(shared_ptr<Scaffold> s, bundle.ref_scaffolds())
    {
        if (s->annotated_tss_id() != "")
        {
            perform_tss_analysis = final_est_run;
        }
        if (s->annotated_protein_id() != "")
        {
            perform_cds_analysis = final_est_run;
        }
    }

    set<shared_ptr<ReadGroupProperties const> > rg_props;
    for (size_t i = 0; i < sample_factory.factories().size(); ++i)
    {
        shared_ptr<BundleFactory> bf = sample_factory.factories()[i];
        rg_props.insert(bf->read_group_properties());
    }
    
    sample_abundance_worker(boost::cref(locus_tag),
                            boost::cref(rg_props),
                            boost::ref(*abundance),
                            &bundle,
                            perform_cds_analysis,
                            perform_tss_analysis);
    
#if ENABLE_THREADS
    variance_info_lock.lock();
#endif
    
#if ENABLE_THREADS
    variance_info_lock.unlock();
#endif
    
    ///////////////////////////////////////////////
    
    
    BOOST_FOREACH(shared_ptr<Scaffold> ref_scaff,  bundle.ref_scaffolds())
    {
        ref_scaff->clear_hits();
    }
    
    launcher->abundance_avail(locus_tag, abundance, factory_id);
    launcher->test_finished_loci();
    
#if !ENABLE_THREADS
    // If Cuffdiff was built without threads, we need to manually invoke 
    // the testing functor, which will check to see if all the workers
    // are done, and if so, perform the cross sample testing.
    //launcher->test_finished_loci();
#endif
}

void dump_locus_variance_info(const string& filename)
{
#if ENABLE_THREADS
    variance_info_lock.lock();
#endif
    
    FILE* fdump = fopen(filename.c_str(), "w");
    
    fprintf(fdump, 
            "condition\tdescription\tlocus_counts\tempir_var\tlocus_fit_var\tsum_iso_fit_var\tcross_replicate_js\tnum_transcripts\tbayes_gamma_trace\tempir_gamma_trace\tcount_mean\tgamma_var\tlocus_salient_frags\tlocus_total_frags\tcount_sharing\n");
    BOOST_FOREACH (LocusVarianceInfo& L, locus_variance_info_table)
    {
        for (size_t i = 0; i < L.gamma.size(); ++i)
        {
            fprintf(fdump, "%d\t%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", L.factory_id, L.transcript_ids[i].c_str(), L.count_mean, L.count_empir_var, L.locus_count_fitted_var, L.isoform_fitted_var_sum, L.cross_replicate_js, L.num_transcripts, L.bayes_gamma_trace, L.empir_gamma_trace,L.gamma[i],L.gamma_var[i], L.locus_salient_frags, L.locus_total_frags, L.count_sharing[i]);
        }
        
    }
    
#if ENABLE_THREADS
    variance_info_lock.unlock();
#endif
}

void filter_group_for_js_testing(vector<vector<AbundanceGroup> >& source_groups)
{
    if (source_groups.empty())
        return;
    
    // iterate over transcript groups
    for (size_t ab_group_idx = 0; ab_group_idx < source_groups[0].size(); ++ab_group_idx)
    {
        vector<bool> to_keep(source_groups[0][ab_group_idx].abundances().size(), false);
        // iterate over samples
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            // iterate over each member of the current abundance group
            AbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
            {
                const Abundance& ab = *(ab_group.abundances()[ab_idx]);
                if (ab.num_fragments() && ab.effective_length())
                {
                    double frags_per_kb = ab.num_fragments() / (ab.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        to_keep[ab_idx] = true;
                }
                
            }
        }
        
        // Now that we know which ones we want to keep, get rid of the rest
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            AbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            AbundanceGroup f = ab_group;
            ab_group.filter_group(to_keep, f);
            ab_group = f;
        }
    }
}

// The two functions below just clear the FPKM samples and other data used to estimate
// the FPKM distributions during testing.  Once testing for a locus is complete,
// we can throw that stuff out.  We do so to save memory over the run.
void clear_samples_from_fpkm_tracking_table(const string& locus_desc, FPKMTrackingTable& fpkm_tracking)
{
    FPKMTrackingTable::iterator itr = fpkm_tracking.find(locus_desc);
    if (itr == fpkm_tracking.end())
        return;
    
    for (size_t i = 0; i < itr->second.fpkm_series.size(); ++i)
    {
        itr->second.fpkm_series[i].fpkm_samples.clear();
        std::vector<double>().swap(itr->second.fpkm_series[i].fpkm_samples);
        //itr->second.fpkm_series[i].fpkm_samples.swap(vector<double>(itr->second.fpkm_series[i].fpkm_samples));
    }
}

void clear_samples_from_tracking_table(shared_ptr<SampleAbundances> sample, Tracking& tracking)
{
    for (size_t k = 0; k < sample->transcripts.abundances().size(); ++k)
    {
        const Abundance& abundance = *(sample->transcripts.abundances()[k]);
        const string& desc = abundance.description();
        clear_samples_from_fpkm_tracking_table(desc, tracking.isoform_fpkm_tracking);
    }

    for (size_t k = 0; k < sample->cds.size(); ++k)
    {
        const Abundance& abundance = sample->cds[k];
        const string& desc = abundance.description();
        clear_samples_from_fpkm_tracking_table(desc, tracking.cds_fpkm_tracking);
    }
    
    for (size_t k = 0; k < sample->primary_transcripts.size(); ++k)
    {
        const Abundance& abundance = sample->primary_transcripts[k];
        const string& desc = abundance.description();
        clear_samples_from_fpkm_tracking_table(desc, tracking.tss_group_fpkm_tracking);
    }

    for (size_t k = 0; k < sample->genes.size(); ++k)
    {
        const Abundance& abundance = sample->genes[k];
        const string& desc = abundance.description();
        clear_samples_from_fpkm_tracking_table(desc, tracking.gene_fpkm_tracking);
    }
}

int total_tests = 0;
void test_differential(const string& locus_tag,
					   const vector<shared_ptr<SampleAbundances> >& samples,
                       const vector<pair<size_t, size_t> >& contrasts,
					   Tests& tests,
					   Tracking& tracking)
{
	if (samples.empty())
		return;
    
    if (no_differential == true)
    {
//#if ENABLE_THREADS
//        test_storage_lock.unlock();
//#endif
        for (size_t i = 0; i < samples.size(); ++i)
        {
            clear_samples_from_tracking_table(samples[i], tracking);
        }
        return;
    }
	
    vector<vector<AbundanceGroup> > filtered_primary_trans_groups;
    vector<vector<AbundanceGroup> > filtered_promoter_groups;
    vector<vector<AbundanceGroup> > filtered_cds_groups;

    for (size_t i = 0; i < samples.size(); ++i)
    {
        filtered_primary_trans_groups.push_back(samples[i]->primary_transcripts);
        filtered_promoter_groups.push_back(samples[i]->gene_primary_transcripts);
        filtered_cds_groups.push_back(samples[i]->gene_cds);
    }
    
    filter_group_for_js_testing(filtered_primary_trans_groups);
    filter_group_for_js_testing(filtered_promoter_groups);
    filter_group_for_js_testing(filtered_cds_groups);
    
    
    // Perform pairwise significance testing between samples. If this is a
    // time series, only test between successive pairs of samples, as supplied 
    // by the user.
    
    for (size_t contrast_idx = 0; contrast_idx < contrasts.size(); ++contrast_idx)
    {
        size_t i = contrasts[contrast_idx].first;
        size_t j = contrasts[contrast_idx].second;
        //            bool enough_reads = (samples[i]->cluster_mass >= min_read_count ||
        //                                 samples[j]->cluster_mass >= min_read_count);
        assert (samples[i]->transcripts.abundances().size() ==
                samples[j]->transcripts.abundances().size());
        for (size_t k = 0; k < samples[i]->transcripts.abundances().size(); ++k)
        {
            const Abundance& curr_abundance = *(samples[j]->transcripts.abundances()[k]);
            const Abundance& prev_abundance = *(samples[i]->transcripts.abundances()[k]);
            const string& desc = curr_abundance.description();
            FPKMTrackingTable::iterator itr = tracking.isoform_fpkm_tracking.find(desc);
            assert (itr != tracking.isoform_fpkm_tracking.end());
            
            bool enough_reads = false;
            if (curr_abundance.num_fragments() && curr_abundance.effective_length())
            {
                double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            if (prev_abundance.num_fragments() && prev_abundance.effective_length())
            {
                double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            
//            if (enough_reads)
//            {
//                if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
//                    enough_reads = false;
//            }
            
            SampleDifference test;
            test = get_de_tests(desc,
                                itr->second.fpkm_series[j],
                                itr->second.fpkm_series[i],
                                //tests.isoform_de_tests[i][j],
                                enough_reads);
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.isoform_de_tests[i][j].insert(make_pair(desc,
                                                                     test));
            
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = curr_abundance.gene_id();
            meta_data->gene_names = curr_abundance.gene_name();
            meta_data->protein_ids = curr_abundance.protein_id();
            meta_data->locus_desc = curr_abundance.locus_tag();
            meta_data->description = curr_abundance.description();
            inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        for (size_t k = 0; k < samples[i]->cds.size(); ++k)
        {
            const Abundance& curr_abundance = samples[i]->cds[k];
            const Abundance& prev_abundance = samples[j]->cds[k];
            
            const string& desc = curr_abundance.description();
            FPKMTrackingTable::iterator itr = tracking.cds_fpkm_tracking.find(desc);
            assert (itr != tracking.cds_fpkm_tracking.end());
            
            bool enough_reads = false;
            if (curr_abundance.num_fragments() && curr_abundance.effective_length())
            {
                double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            if (prev_abundance.num_fragments() && prev_abundance.effective_length())
            {
                double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            
//            if (enough_reads)
//            {
//                if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
//                    enough_reads = false;
//            }
            
            
            SampleDifference test;
            test = get_de_tests(desc,
                                itr->second.fpkm_series[j],
                                itr->second.fpkm_series[i],
                                //tests.cds_de_tests[i][j],
                                enough_reads);
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.cds_de_tests[i][j].insert(make_pair(desc,
                                                                 test));
            
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = curr_abundance.gene_id();
            meta_data->gene_names = curr_abundance.gene_name();
            meta_data->protein_ids = curr_abundance.protein_id();
            meta_data->locus_desc = curr_abundance.locus_tag();
            meta_data->description = curr_abundance.description();
            inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
        {
            const Abundance& curr_abundance = samples[i]->primary_transcripts[k];
            const Abundance& prev_abundance = samples[j]->primary_transcripts[k];
            
            const string& desc = curr_abundance.description();
            FPKMTrackingTable::iterator itr = tracking.tss_group_fpkm_tracking.find(desc);
            assert (itr != tracking.tss_group_fpkm_tracking.end());
            
            bool enough_reads = false;
            if (curr_abundance.num_fragments() && curr_abundance.effective_length())
            {
                double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            if (prev_abundance.num_fragments() && prev_abundance.effective_length())
            {
                double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            
//            if (enough_reads)
//            {
//                if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
//                    enough_reads = false;
//            }
            
            
            SampleDifference test;
            test = get_de_tests(desc,
                                itr->second.fpkm_series[j],
                                itr->second.fpkm_series[i],
                                //tests.tss_group_de_tests[i][j],
                                enough_reads);
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.tss_group_de_tests[i][j].insert(make_pair(desc,
                                                                       test));
            
            
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = curr_abundance.gene_id();
            meta_data->gene_names = curr_abundance.gene_name();
            meta_data->protein_ids = curr_abundance.protein_id();
            meta_data->locus_desc = curr_abundance.locus_tag();
            meta_data->description = curr_abundance.description();
            inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        for (size_t k = 0; k < samples[i]->genes.size(); ++k)
        {
            const AbundanceGroup& curr_abundance = samples[i]->genes[k];
            const AbundanceGroup& prev_abundance = samples[j]->genes[k];
            const string& desc = curr_abundance.description();
            FPKMTrackingTable::iterator itr = tracking.gene_fpkm_tracking.find(desc);
            assert (itr != tracking.gene_fpkm_tracking.end());
            
            bool enough_reads = false;
            if (curr_abundance.num_fragments() && curr_abundance.effective_length())
            {
                double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            if (prev_abundance.num_fragments() && prev_abundance.effective_length())
            {
                double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                if (frags_per_kb >= min_read_count)
                    enough_reads = true;
            }
            
            SampleDifference test;
            test = get_de_tests(desc,
                                itr->second.fpkm_series[j],
                                itr->second.fpkm_series[i],
                                //tests.gene_de_tests[i][j],
                                enough_reads);
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.gene_de_tests[i][j].insert(make_pair(desc,
                                                                  test));
            
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = curr_abundance.gene_id();
            meta_data->gene_names = curr_abundance.gene_name();
            meta_data->protein_ids = curr_abundance.protein_id();
            meta_data->locus_desc = curr_abundance.locus_tag();
            meta_data->description = curr_abundance.description();
            inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        // Skip all the JS based testing for genes with an isoform switch?
        if (no_js_tests)
            continue;
        
        // FIXME: the code below will not properly test for differential
        // splicing/promoter use when a gene (e.g.) occupies two
        // disjoint bundles.  We need to store the covariance matrices (etc)
        // in the FPKMContexts to handle that case properly.
        
        // Differential promoter use
        for (size_t k = 0; k < samples[i]->gene_primary_transcripts.size(); ++k)
        {
            const AbundanceGroup& curr_abundance = filtered_promoter_groups[i][k];
            const AbundanceGroup& prev_abundance = filtered_promoter_groups[j][k];
            const string& desc = curr_abundance.description();
            
            bool enough_reads = (curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test &&
                                 prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);
            SampleDifference test;
            test = get_ds_tests(prev_abundance,
                                curr_abundance,
                                enough_reads);
            
            // The filtered group might be empty, so let's grab metadata from
            // the unfiltered group
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = samples[i]->gene_primary_transcripts[k].gene_id();
            meta_data->gene_names = samples[i]->gene_primary_transcripts[k].gene_name();
            meta_data->protein_ids = samples[i]->gene_primary_transcripts[k].protein_id();
            meta_data->locus_desc = samples[i]->gene_primary_transcripts[k].locus_tag();
            meta_data->description = samples[i]->gene_primary_transcripts[k].description();
            
            test.meta_data = meta_data;
            
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.diff_promoter_tests[i][j].insert(make_pair(desc,test));
            inserted.first->second = test;
            
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        // Differential coding sequence output
        for (size_t k = 0; k < samples[i]->gene_cds.size(); ++k)
        {
            const AbundanceGroup& curr_abundance = filtered_cds_groups[i][k];
            const AbundanceGroup& prev_abundance = filtered_cds_groups[j][k];
            const string& desc = curr_abundance.description();
            
            bool enough_reads =  (curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test &&
                                  prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);
            SampleDifference test;
            test = get_ds_tests(prev_abundance,
                                curr_abundance,
                                enough_reads);
            
            // The filtered group might be empty, so let's grab metadata from
            // the unfiltered group
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = samples[i]->gene_cds[k].gene_id();
            meta_data->gene_names = samples[i]->gene_cds[k].gene_name();
            meta_data->protein_ids = samples[i]->gene_cds[k].protein_id();
            meta_data->locus_desc = samples[i]->gene_cds[k].locus_tag();
            meta_data->description = samples[i]->gene_cds[k].description();
            
            test.meta_data = meta_data;
            
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.diff_cds_tests[i][j].insert(make_pair(desc,test));
            inserted.first->second = test;
            
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
        // Differential splicing of primary transcripts
        for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
        {
            const AbundanceGroup& curr_abundance = filtered_primary_trans_groups[i][k];
            const AbundanceGroup& prev_abundance = filtered_primary_trans_groups[j][k];
            const string& desc = curr_abundance.description();
            
            bool enough_reads = (curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test &&
                                 prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);
            
            SampleDifference test;
            test = get_ds_tests(prev_abundance, 
                                curr_abundance,
                                enough_reads);
            
            // The filtered group might be empty, so let's grab metadata from
            // the unfiltered group
            shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
            
            meta_data->gene_ids = samples[i]->primary_transcripts[k].gene_id();
            meta_data->gene_names = samples[i]->primary_transcripts[k].gene_name();
            meta_data->protein_ids = samples[i]->primary_transcripts[k].protein_id();
            meta_data->locus_desc = samples[i]->primary_transcripts[k].locus_tag();
            meta_data->description = samples[i]->primary_transcripts[k].description();
            
            test.meta_data = meta_data;
            
#if ENABLE_THREADS
            test_storage_lock.lock();
#endif
            pair<SampleDiffs::iterator, bool> inserted;
            inserted = tests.diff_splicing_tests[i][j].insert(make_pair(desc,test)); 
            inserted.first->second = test;
            
#if ENABLE_THREADS
            test_storage_lock.unlock();
#endif
        }
        
    }
    
    for (size_t i = 0; i < samples.size(); ++i)
    {
        clear_samples_from_tracking_table(samples[i], tracking);
    }
}
