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

#include "abundances.h"
#include "differential.h"
#include "clustering.h"
#include "differential.h"

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
    foreach (shared_ptr<SampleAbundances> ab, abundances)
    {
        if (!ab)
        {
            return false;
        }
    }
    return true;
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void TestLauncher::perform_testing(vector<shared_ptr<SampleAbundances> >& abundances)
{
    if (_p_bar)
    {
        verbose_msg("Testing for differential expression and regulation in locus [%s]\n", abundances.front()->locus_tag.c_str());
        _p_bar->update(abundances.front()->locus_tag.c_str(), 1);
    }
    
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
    
    test_differential(abundances.front()->locus_tag, abundances, *_tests, *_tracking, _samples_are_time_series);
}

void TestLauncher::test_finished_loci()
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif  
    // In some abundance runs, we don't actually want to perform testing 
    // (eg initial quantification before bias correction).
    // _tests and _tracking will be NULL in these cases.
    if (_tests != NULL && _tracking != NULL && !_samples.empty())
    {
        launcher_sample_table::iterator itr = _samples.begin(); 
        while(itr != _samples.end())
        {
            if (all_samples_reported_in(itr->second))
            {
                perform_testing(itr->second);
                itr = _samples.erase(itr);
            }
            else
            {
                ++itr;
            }
        }
    }
}

// This performs a between-group test on an isoform or TSS grouping, on two 
// different samples.
bool test_diffexp(const FPKMContext& curr,
				  const FPKMContext& prev,
				  SampleDifference& test)
{
	bool performed_test = false;
	if (curr.FPKM > 0.0 && prev.FPKM > 0.0)
	{
		//assert (curr.FPKM_variance > 0.0 && prev.FPKM_variance > 0.0);
//		double log_curr = log(curr.counts);
//		double log_prev = log(prev.counts);
        
        double stat = 0.0;
        double p_value = 1.0;
        
        if (curr.FPKM_variance > 0.0 || prev.FPKM_variance > 0.0)
        {
            double curr_log_fpkm_var = (curr.FPKM_variance) / (curr.FPKM * curr.FPKM);
            double prev_log_fpkm_var = (prev.FPKM_variance) / (prev.FPKM * prev.FPKM);
            
            double numerator = log(prev.FPKM / curr.FPKM);
            
            double denominator = sqrt(prev_log_fpkm_var + curr_log_fpkm_var);
            stat = numerator / denominator;
        
		
            normal norm;
            double t1, t2;
            if (stat > 0.0)
            {
                t1 = stat;
                t2 = -stat;
            }
            else
            {
                t1 = -stat;
                t2 = stat;
            }
            
            if (isnan(t1) || isinf(t1) || isnan(t2) || isnan(t2))
            {
                
                //fprintf(stderr, "Warning: test statistic is NaN! %s (samples %lu and %lu)\n", test.locus_desc.c_str(), test.sample_1, test.sample_2);
                p_value = 1.0;
            }
            else
            {
                double tail_1 = cdf(norm, t1);
                double tail_2 = cdf(norm, t2);
                p_value = 1.0 - (tail_1 - tail_2);                
            }
        }
		
		double differential = log(curr.FPKM) - log(prev.FPKM);
		
		//test = SampleDifference(sample1, sample2, prev.FPKM, curr.FPKM, stat, p_value, transcript_group_id);
		test.p_value = p_value;
		test.differential = differential;
		test.test_stat = stat;
		test.value_1 = prev.FPKM;
		test.value_2 = curr.FPKM;
		
		performed_test = true;
	}
	else
	{
		if (curr.FPKM > 0.0)
		{
			//test = SampleDifference(sample1, sample2, 0, curr.FPKM, DBL_MAX, 0, transcript_group_id); 
			test.p_value = 0;
			test.differential = numeric_limits<double>::max();;
			test.test_stat = numeric_limits<double>::max();
			test.value_1 = 0;
			test.value_2 = curr.FPKM;
			performed_test = true;
		}
		else if (prev.FPKM > 0.0)
		{
			//test = SampleDifference(sample1, sample2, prev.FPKM, 0, -DBL_MAX, 0, transcript_group_id); 
			test.p_value = 0;
			test.differential = -numeric_limits<double>::max();
			test.test_stat = -numeric_limits<double>::max();
			test.value_1 = prev.FPKM;
			test.value_2 = 0;
			performed_test = true;
		}
	}	
	
	test.test_status = performed_test ? OK : NOTEST;
	return performed_test;
}

// This performs between-group tests on isoforms or TSS groupings in a single
// locus, on two different samples.
pair<int, SampleDiffs::iterator>  get_de_tests(const string& description,
                 const FPKMContext& prev_abundance,
				 const FPKMContext& curr_abundance,
				 SampleDiffs& de_tests,
				 bool enough_reads)
{
	int total_iso_de_tests = 0;
			
	SampleDifference test;
    
	pair<SampleDiffs::iterator, bool> inserted;
//	inserted = de_tests.insert(make_pair(curr_abundance.description(),
//										 SampleDifference())); 
    inserted = de_tests.insert(make_pair(description,
    									 SampleDifference())); 
    
    const FPKMContext& r1 = curr_abundance;
    const FPKMContext& r2 = prev_abundance;
    
	if (curr_abundance.status == NUMERIC_FAIL || 
        prev_abundance.status == NUMERIC_FAIL)
    {
        test.test_stat = 0;
		test.p_value = 1.0;
		test.differential = 0.0;
		test.test_status = FAIL;
    }
    else if (curr_abundance.status == NUMERIC_LOW_DATA && 
             prev_abundance.status == NUMERIC_LOW_DATA)
    {
        // perform the test, but mark it as not significant and don't add it to the 
        // pile. This way we don't penalize for multiple testing, but users can still
        // see the fold change.
		test_diffexp(r1, r2, test);
        test.test_stat = 0;
        test.p_value = 1.0;
        //test.differential = 0.0;
		
		test.test_status = LOWDATA;
    }
    else // at least one is OK, the other might be LOW_DATA
	{
		test.test_status = FAIL;

		if (test_diffexp(r1, r2, test))
		{
			total_iso_de_tests++;
		}
		else
		{
			test.test_stat = 0;
			test.p_value = 1.0;
			test.differential = 0.0;
		}
		if (enough_reads) 
			test.test_status = OK;
		else
			test.test_status = NOTEST;
		
	}
	
    
	inserted.first->second = test;
	
	return make_pair(total_iso_de_tests, inserted.first);
}


// This performs within-group tests on a set of isoforms or a set of TSS groups.
// This is a way of looking for meaningful differential splicing or differential
// promoter use.
void get_ds_tests(const AbundanceGroup& prev_abundance,
				  const AbundanceGroup& curr_abundance,
				  SampleDiffs& diff_tests,
				  bool enough_reads)
{	
	const string& name = curr_abundance.description();
	
	pair<SampleDiffs::iterator, bool> inserted;
	inserted = diff_tests.insert(make_pair(name,SampleDifference())); 
	SampleDifference test;
	
    test.gene_ids = curr_abundance.gene_id();
	test.gene_names = curr_abundance.gene_name();
	test.protein_ids = curr_abundance.protein_id();
	test.locus_desc = curr_abundance.locus_tag();
	test.description = curr_abundance.description();
	
	test.test_status = NOTEST;
	
	bool prev_status = curr_abundance.status();
	bool curr_status = prev_abundance.status();
	
	if (prev_status == NUMERIC_OK && prev_abundance.num_fragments() > 0 &&
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
		
		double js = jensen_shannon_div(sample_kappas);
		if (js == 0.0 || isnan(js) || isinf(js))
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
			ublas::vector<double> js_gradient;
			jensen_shannon_gradient(sample_kappas, js, js_gradient);
			
			vector<ublas::matrix<double> > covariances;
			
			covariances.push_back(prev_abundance.kappa_cov());
			covariances.push_back(curr_abundance.kappa_cov());
			
			ublas::matrix<double> js_covariance;
			assert (covariances.size() > 0);
			for (size_t i = 0; i < covariances.size(); ++i)
			{
				assert (covariances[i].size1() > 0 && covariances[i].size2() > 0);
			}
			make_js_covariance_matrix(covariances,js_covariance);
			assert (js_covariance.size1() > 0 && js_covariance.size2() > 0);
			
			//cout << "grad: " << js_gradient << endl;
			//cout << "js_cov: " << js_covariance << endl;
			//cout << prod(js_covariance, js_gradient) << endl;
			
			double js_var = inner_prod(js_gradient, 
									   prod(js_covariance, js_gradient));
#ifdef DEBUG
			if (isinf(js_var) || isnan(js_var))
			{
				cerr << "grad: " << js_gradient << endl;
				cerr << "js_cov: " << js_covariance << endl;
				cerr << prod(js_covariance, js_gradient) << endl;	
			}
#endif
			if (js_var <= 0.0)
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
                // We're dealing with a standard normal that's been truncated below zero, so
                // we need to multiply the baseline pdf and cdf by 2.0.
                
				normal test_dist(0,1.0);
				//double denom = sqrt(js_var);
				double p = js/sqrt(js_var);
                test.test_stat = 2 * pdf(test_dist, p);
				test.p_value = cdf(test_dist, p);
				test.value_1 = 0;
				test.value_2 = 0;
				test.differential = js;
				test.test_status = enough_reads ? OK : NOTEST;
			}
			if (isinf(test.test_stat) || isnan(test.test_stat))
			{
				fprintf(stderr, "Warning: test stat is invalid!\n");
				exit(1);
			}
		}

		inserted.first->second = test;
	}
	else // we won't even bother with the JS-based testing in LOWDATA cases.
	{
        if (prev_status == NUMERIC_OK && curr_status == NUMERIC_OK)
            test.test_status = NOTEST;
        else if (prev_status == NUMERIC_FAIL || curr_status == NUMERIC_FAIL)
            test.test_status = FAIL;
        else
            test.test_status = LOWDATA;
            
		test.test_stat = 0;
		test.p_value = 0.0;
		test.differential = 0.0;
		inserted.first->second = test;
	}
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
								 ab.FPKM(), 
								 ab.FPKM_variance(),
                                 ab.status());
    
    
	
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
        existing.counts += r1.counts;
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

     // TODO: remove this assert
    //assert (inserted.first->second.fpkm_series.size() <= 2);
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

#if ENABLE_THREADS
mutex test_storage_lock; // don't modify the above struct without locking here
#endif

void sample_abundance_worker(const string& locus_tag,
                             SampleAbundances& sample,
                             HitBundle* sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis)
{
    vector<shared_ptr<Abundance> > abundances;
    
    foreach(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        TranscriptAbundance* pT = new TranscriptAbundance;
        pT->transfrag(s);
        shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }
    
    sample.transcripts = AbundanceGroup(abundances);
    
    vector<MateHit> hits_in_cluster;
    
    get_alignments_from_scaffolds(sample.transcripts.abundances(),
                                  hits_in_cluster);
    
    // Compute the individual transcript FPKMs via each sample's 
    // AbundanceGroup for this locus.
    
    sample.transcripts.calculate_abundance(hits_in_cluster);
    
    // Cluster transcripts by gene_id
    vector<AbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	foreach(AbundanceGroup& ab_group, transcripts_by_gene_id)
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
        cluster_transcripts<ConnectByAnnotatedProteinId>(sample.transcripts,
                                                         transcripts_by_cds,
                                                         &cds_gamma_cov);
        foreach(AbundanceGroup& ab_group, transcripts_by_cds)
        {
            ab_group.locus_tag(locus_tag);
            set<string> protein_ids = ab_group.protein_id();
            assert (protein_ids.size() == 1);
            string desc = *(protein_ids.begin()); 
            assert (desc != "");
            ab_group.description(*(protein_ids.begin()));
        }
        
        sample.cds = transcripts_by_cds;
        
        // Group the CDS clusters by gene
        vector<shared_ptr<Abundance> > cds_abundances;
        foreach (AbundanceGroup& ab_group, sample.cds)
        {
            cds_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
        }
        AbundanceGroup cds(cds_abundances,
                           cds_gamma_cov);
        
        vector<AbundanceGroup> cds_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(cds,
                                                      cds_by_gene);
        
        foreach(AbundanceGroup& ab_group, cds_by_gene)
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
        cluster_transcripts<ConnectByAnnotatedTssId>(sample.transcripts,
                                                     transcripts_by_tss,
                                                     &tss_gamma_cov);
        
        foreach(AbundanceGroup& ab_group, transcripts_by_tss)
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
        foreach (AbundanceGroup& ab_group, sample.primary_transcripts)
        {
            primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
        }
        
        AbundanceGroup primary_transcripts(primary_transcript_abundances,
                                           tss_gamma_cov);
        
        vector<AbundanceGroup> primary_transcripts_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(primary_transcripts,
                                                      primary_transcripts_by_gene);
        
        foreach(AbundanceGroup& ab_group, primary_transcripts_by_gene)
        {
            ab_group.locus_tag(locus_tag);
            set<string> gene_ids = ab_group.gene_id();
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_primary_transcripts = primary_transcripts_by_gene;
    }
}

void sample_worker(const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   shared_ptr<SampleAbundances> abundance,
                   size_t factory_id,
                   shared_ptr<bool> non_empty,
                   shared_ptr<TestLauncher> launcher)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    HitBundle bundle;
    *non_empty = sample_factory.next_bundle(bundle);
    
    if (!*non_empty || (!corr_multi && !final_est_run && bundle.ref_scaffolds().size() != 1)) // Only learn on single isoforms
    {
#if !ENABLE_THREADS
        // If Cuffdiff was built without threads, we need to manually invoke 
        // the testing functor, which will check to see if all the workers
        // are done, and if so, perform the cross sample testing.
        launcher();
#endif
    	return;
    }
    
    abundance->cluster_mass = bundle.mass();
    
    char bundle_label_buf[2048];
    sprintf(bundle_label_buf, 
            "%s:%d-%d", 
            rt.get_name(bundle.ref_id()),
            bundle.left(),
            bundle.right());
    string locus_tag = bundle_label_buf;
    
    launcher->register_locus(locus_tag);
    
    abundance->locus_tag = locus_tag;
    
    bool perform_cds_analysis = final_est_run;
    bool perform_tss_analysis = final_est_run;

    foreach(shared_ptr<Scaffold> s, bundle.ref_scaffolds())
    {
        if (s->annotated_tss_id() == "")
        {
            perform_tss_analysis = false;
        }
        if (s->annotated_protein_id() == "")
        {
            perform_cds_analysis = false;
        }
    }

    sample_abundance_worker(boost::cref(locus_tag),
                            boost::ref(*abundance),
                            &bundle,
                            perform_cds_analysis,
                            perform_tss_analysis);
    
    foreach(shared_ptr<Scaffold> ref_scaff,  bundle.ref_scaffolds())
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

int total_tests = 0;
void test_differential(const string& locus_tag,
					   const vector<shared_ptr<SampleAbundances> >& samples,
					   Tests& tests,
					   Tracking& tracking,
                       bool samples_are_time_series)
{
	if (samples.empty())
		return;
    
#if ENABLE_THREADS
	test_storage_lock.lock();
    total_tests++;
#endif
    
    //fprintf(stderr, "\nTesting in %s (%d total tests)\n", locus_tag.c_str(), total_tests);
    
	// Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
	for (size_t i = 0; i < samples.size(); ++i)
	{
		const AbundanceGroup& ab_group = samples[i]->transcripts;
		foreach (shared_ptr<Abundance> ab, ab_group.abundances())
		{
			add_to_tracking_table(i, *ab, tracking.isoform_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i]->cds)
		{
			add_to_tracking_table(i, ab, tracking.cds_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i]->primary_transcripts)
		{
			add_to_tracking_table(i, ab, tracking.tss_group_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i]->genes)
		{
			add_to_tracking_table(i, ab, tracking.gene_fpkm_tracking);
		}
	}
	
    // Perform pairwise significance testing between samples. If this is a
    // time series, only test between successive pairs of samples, as supplied 
    // by the user.
	for (size_t i = 1; i < samples.size(); ++i)
	{
		//bool multi_transcript_locus = samples[i]->transcripts.abundances().size() > 1;
		
        int sample_to_start_test_against = 0;
        if (samples_are_time_series)
            sample_to_start_test_against = i - 1;
        
        for (size_t j = sample_to_start_test_against; j < i; ++j)
        {
            bool enough_reads = (samples[i]->cluster_mass >= min_read_count ||
                                 samples[j]->cluster_mass >= min_read_count);
            assert (samples[i]->transcripts.abundances().size() == 
                    samples[j]->transcripts.abundances().size());
            for (size_t k = 0; k < samples[i]->transcripts.abundances().size(); ++k)
            {
                const Abundance& curr_abundance = *(samples[j]->transcripts.abundances()[k]);
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.isoform_fpkm_tracking.find(desc);
                assert (itr != tracking.isoform_fpkm_tracking.end());
                
                pair<int, SampleDiffs::iterator> result;
                result = get_de_tests(desc,
                                      itr->second.fpkm_series[j], 
                                      itr->second.fpkm_series[i],
                                      tests.isoform_de_tests[i][j],
                                      enough_reads);
                
                result.second->second.gene_ids = curr_abundance.gene_id();
                result.second->second.gene_names = curr_abundance.gene_name();
                result.second->second.protein_ids = curr_abundance.protein_id();
                result.second->second.locus_desc = curr_abundance.locus_tag();
                result.second->second.description = curr_abundance.description();
            }
            
            for (size_t k = 0; k < samples[i]->cds.size(); ++k)
            {
                const Abundance& curr_abundance = samples[j]->cds[k];
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.cds_fpkm_tracking.find(desc);
                assert (itr != tracking.cds_fpkm_tracking.end());
                
                pair<int, SampleDiffs::iterator> result;
                result = get_de_tests(desc,
                             itr->second.fpkm_series[j], 
                             itr->second.fpkm_series[i],
                             tests.cds_de_tests[i][j],
                             enough_reads);
                
                result.second->second.gene_ids = curr_abundance.gene_id();
                result.second->second.gene_names = curr_abundance.gene_name();
                result.second->second.protein_ids = curr_abundance.protein_id();
                result.second->second.locus_desc = curr_abundance.locus_tag();
                result.second->second.description = curr_abundance.description();
            }
            
            for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
            {
                const Abundance& curr_abundance = samples[j]->primary_transcripts[k];
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.tss_group_fpkm_tracking.find(desc);
                assert (itr != tracking.tss_group_fpkm_tracking.end());
                
                pair<int, SampleDiffs::iterator> result;
                result = get_de_tests(desc,
                             itr->second.fpkm_series[j], 
                             itr->second.fpkm_series[i],
                             tests.tss_group_de_tests[i][j],
                             enough_reads);
                
                result.second->second.gene_ids = curr_abundance.gene_id();
                result.second->second.gene_names = curr_abundance.gene_name();
                result.second->second.protein_ids = curr_abundance.protein_id();
                result.second->second.locus_desc = curr_abundance.locus_tag();
                result.second->second.description = curr_abundance.description();
            }
            
            for (size_t k = 0; k < samples[i]->genes.size(); ++k)
            {
                const Abundance& curr_abundance = samples[j]->genes[k];
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.gene_fpkm_tracking.find(desc);
                assert (itr != tracking.gene_fpkm_tracking.end());
                
                pair<int, SampleDiffs::iterator> result;
                result = get_de_tests(desc,
                             itr->second.fpkm_series[j], 
                             itr->second.fpkm_series[i],
                             tests.gene_de_tests[i][j],
                             enough_reads);
                
                result.second->second.gene_ids = curr_abundance.gene_id();
                result.second->second.gene_names = curr_abundance.gene_name();
                result.second->second.protein_ids = curr_abundance.protein_id();
                result.second->second.locus_desc = curr_abundance.locus_tag();
                result.second->second.description = curr_abundance.description();
            }
            
            // FIXME: the code below will not properly test for differential
            // splicing/promoter use when a gene (e.g.) occupies two
            // disjoint bundles.  We need to store the covariance matrices (etc)
            // in the FPKMContexts to handle that case properly.
            
            // Differential promoter use
            for (size_t k = 0; k < samples[i]->gene_primary_transcripts.size(); ++k)
            {
                get_ds_tests(samples[j]->gene_primary_transcripts[k], 
                             samples[i]->gene_primary_transcripts[k],
                             tests.diff_promoter_tests[i][j],
                             enough_reads);
            }
            
            // Differential coding sequence output
            for (size_t k = 0; k < samples[i]->gene_cds.size(); ++k)
            {
                get_ds_tests(samples[j]->gene_cds[k], 
                             samples[i]->gene_cds[k],
                             tests.diff_cds_tests[i][j],
                             enough_reads);
            }
            
            // Differential splicing of primary transcripts
            for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
            {
                get_ds_tests(samples[j]->primary_transcripts[k], 
                             samples[i]->primary_transcripts[k],
                             tests.diff_splicing_tests[i][j],
                             enough_reads);
            }
        }
	}
	
#if ENABLE_THREADS
	test_storage_lock.unlock();
#endif
}
