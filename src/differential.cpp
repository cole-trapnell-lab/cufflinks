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

using namespace std;

double min_read_count = 1000;

// This performs a between-group test on an isoform or TSS grouping, on two 
// different samples.
bool test_diffexp(const FPKMContext& curr,
				  const FPKMContext& prev,
				  SampleDifference& test)
{
	bool performed_test = false;
	if (curr.FPKM > 0.0 && prev.FPKM > 0.0)
	{
		assert (curr.FPKM_variance > 0.0 && prev.FPKM_variance > 0.0);
		double log_curr = log(curr.counts);
		double log_prev = log(prev.counts);
		
		double curr_log_fpkm_var = (curr.FPKM_variance) / (curr.FPKM * curr.FPKM);
		double prev_log_fpkm_var = (prev.FPKM_variance) / (prev.FPKM * prev.FPKM);
		
		// Note: this is written a little differently in the supplement, 
		// where the numerator of the tests includes an explicit sum of
		// logs to get the isoform counts in log space.  We have already
		// performed X_g * gamma_j here, so we just take the log of 
		// that product.
		
		double numerator = (log_curr + prev.log_sample_mass - log_prev - curr.log_sample_mass);
		double denominator = sqrt(prev_log_fpkm_var + curr_log_fpkm_var);
		double stat = numerator / denominator;
		
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
		double tail_1 = cdf(norm, t1);
		double tail_2 = cdf(norm, t2);
		
		double differential;
		if (curr.FPKM == 0.0 && prev.FPKM != 0.0)
		{
			differential = numeric_limits<double>::max();;
		}
		else
		{
			differential = log(curr.FPKM) - log(prev.FPKM);
		}
		
		double p_value = 1.0 - (tail_1 - tail_2);
		
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
			test.test_stat = numeric_limits<double>::max();;
			test.value_1 = 0;
			test.value_2 = curr.FPKM;
			performed_test = true;
		}
		else if (prev.FPKM > 0.0)
		{
			//test = SampleDifference(sample1, sample2, prev.FPKM, 0, -DBL_MAX, 0, transcript_group_id); 
			test.p_value = 0;
			test.test_stat = numeric_limits<double>::min();
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
int get_de_tests(const Abundance& prev_abundance,
				 long double prev_mass,
				 const Abundance& curr_abundance,
				 long double curr_mass,
				 SampleDiffs& de_tests,
				 bool enough_reads)
{
	int total_iso_de_tests = 0;
			
	SampleDifference test;
	pair<SampleDiffs::iterator, bool> inserted;
	inserted = de_tests.insert(make_pair(curr_abundance.description(),
										 SampleDifference())); 
	if (curr_abundance.status() == NUMERIC_OK && 
		prev_abundance.status() == NUMERIC_OK && 
		curr_mass > 0 &&
		prev_mass > 0)
	{
		double log_mass_curr = log(curr_mass);
		double log_mass_prev = log(prev_mass);
		
		FPKMContext r1(curr_abundance.num_fragments(), 
					   curr_abundance.FPKM(), 
					   curr_abundance.FPKM_variance(), 
					   log_mass_curr);
		
		FPKMContext r2(prev_abundance.num_fragments(), 
					   prev_abundance.FPKM(), 
					   prev_abundance.FPKM_variance(), 
					   log_mass_prev);
		
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
	else
	{
		test.test_stat = 0;
		test.test_stat = 1.0;
		test.differential = 0.0;
		test.test_status = NOTEST;
	}
	
	test.gene_names = curr_abundance.gene_name();
	test.protein_ids = curr_abundance.protein_id();
	test.locus_desc = curr_abundance.locus_tag();
	test.description = curr_abundance.description();
	inserted.first->second = test;
	
	return total_iso_de_tests;
}

double entropy(const ublas::vector<double>& p)
{
	double e = 0;  
	for (size_t i = 0; i < p.size(); ++i)
	{
		double P = p[i];
		e -= (P * log(P));
	}
	return e;
}

double jensen_shannon_div(vector<ublas::vector<double> >& sample_kappas)
{
	assert (sample_kappas.size() > 1);
	
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		//cerr << sample_kappas[i] << endl;
		double kappa_sum = accumulate(sample_kappas[i].begin(), 
									  sample_kappas[i].end(), 0.0);
		if (abs(kappa_sum - 1.0) > 1e-10)
		{
			cerr << kappa_sum << " " << sample_kappas[i] << endl;
		}
		assert (abs(kappa_sum - 1.0) < 1e-10);
	}
	
	size_t kappa_length = 0;
	for (size_t i = 1; i < sample_kappas.size(); ++i)
	{
		assert (sample_kappas[i].size() == sample_kappas[i-1].size());
		kappa_length = sample_kappas[i].size();
	}
	
	ublas::vector<double> avg_kappas = ublas::zero_vector<double>(kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		//cout << "kappa " << i<< " "<< sample_kappas[i] << endl;
		avg_kappas += sample_kappas[i];
	}
	avg_kappas /= sample_kappas.size();
	//cout << avg_kappas << endl;
	
	double avg_entropy = 0.0;
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		avg_entropy += entropy(sample_kappas[i]);
	}
	avg_entropy /= sample_kappas.size();
	//cout << avg_entropy << endl;
	
	double entropy_avg = entropy(avg_kappas);
	
	double js = entropy_avg - avg_entropy;
	return sqrt(js);
}

void jensen_shannon_gradient(vector<ublas::vector<double> >& sample_kappas,
							 double js,
							 ublas::vector<double>& gradient)
{
	assert (sample_kappas.size() > 1);
	size_t kappa_length = sample_kappas.front().size();
	for (size_t i = 1; i < sample_kappas.size(); ++i)
	{
		assert (sample_kappas[i].size() == sample_kappas[i-1].size());
		kappa_length = sample_kappas[i].size();
	}
	
	if (kappa_length == 0)
		return;
	
	gradient = ublas::zero_vector<double>(sample_kappas.size() * kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		for (size_t j = 0; j < kappa_length; ++j)
		{
			gradient(i*kappa_length + j) = sample_kappas[i](j);
		}
	}
	
	//cout << "t1: " << gradient<< endl;
	
	ublas::vector<double> denoms = ublas::zero_vector<double>(kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		denoms += sample_kappas[i];
	}
	denoms /= sample_kappas.size();
	
	
	//cout << "t2 " << denoms << endl;
	
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		for (size_t j = 0; j < kappa_length; ++j)
		{
			if (denoms(j) == 0.0 || gradient(i*kappa_length + j) == 0.0)
			{
				gradient(i*kappa_length + j) = 0.0;
			}
			else
			{
#ifdef DEBUG
				ublas::vector<double>& grad_tmp = gradient;
#endif
				gradient(i*kappa_length + j) /= denoms(j);
				gradient(i*kappa_length + j) = log(gradient(i*kappa_length + j));
				gradient(i*kappa_length + j) /= sample_kappas.size();
				gradient(i*kappa_length + j) *= (1.0/(2.0*js));
				
#ifdef DEBUG
				if(isinf(gradient(i*kappa_length + j)))
				{
					cerr << grad_tmp << endl;
					cerr << sample_kappas[i] << endl;
					assert (false);
				}
#endif
				
			}
		}
	}
}

void make_js_covariance_matrix(vector<ublas::matrix<double> >& kappa_covariances,
							   ublas::matrix<double>& js_covariance)
{
	size_t kappa_length = 0;
	for (size_t i = 1; i < kappa_covariances.size(); ++i)
	{
		assert (kappa_covariances[i].size1() == kappa_covariances[i-1].size1());
		assert (kappa_covariances[i].size2() == kappa_covariances[i-1].size2());
		
		kappa_length = kappa_covariances[i].size1();
	}
	
	if (kappa_length == 0)
		return;
	
	js_covariance = ublas::zero_matrix<double>(kappa_covariances.size() * kappa_length,
											   kappa_covariances.size() * kappa_length);
	for (size_t i = 0; i < kappa_covariances.size(); ++i)
	{
		for (size_t j = 0; j < kappa_length; ++j)
		{
			for (size_t k = 0; k < kappa_length; ++k)
			{
				js_covariance(i*kappa_length + j, i*kappa_length + k) = 
				kappa_covariances[i](j,k);
			}
		}
	}
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
				normal test_dist(0,1.0);
				//double denom = sqrt(js_var);
				test.test_stat = js;
				double p = js/sqrt(js_var);
				test.p_value = 1.0 - cdf(test_dist, p);
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
	else
	{
		test.test_status = FAIL;
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
void add_to_tracking_table(Abundance& ab,
						   double log_sample_mass,
						   FPKMTrackingTable& track)

{
	pair<FPKMTrackingTable::iterator,bool> inserted;
	pair<string, FPKMTracking > p;
	p = make_pair(ab.description(), FPKMTracking());
	inserted = track.insert(p);
	
	FPKMTracking& fpkm_track = inserted.first->second;
	
	set<string> tss = ab.tss_id();
	set<string> genes = ab.gene_name();
	set<string> proteins = ab.protein_id();
	
	fpkm_track.tss_ids.insert(tss.begin(), tss.end());
	fpkm_track.gene_names.insert(genes.begin(), genes.end());
	fpkm_track.protein_ids.insert(proteins.begin(), proteins.end());
	
	if (inserted.second)
	{
		fpkm_track.locus_tag = ab.locus_tag();
		fpkm_track.description = ab.description();
		const Scaffold* transfrag = ab.transfrag();
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
	}
	
	FPKMContext r1 = FPKMContext(ab.num_fragments(), 
								 ab.FPKM(), 
								 ab.FPKM_variance(), 
								 log_sample_mass);
	
	inserted.first->second.fpkm_series.push_back(r1);
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

struct SampleAbundances
{
	AbundanceGroup transcripts;
	vector<AbundanceGroup> primary_transcripts;
	vector<AbundanceGroup> gene_primary_transcripts;
	vector<AbundanceGroup> cds;
	vector<AbundanceGroup> gene_cds;
	vector<AbundanceGroup> genes;
	double sample_mass;
	double cluster_mass;
};

#if ENABLE_THREADS
mutex test_storage_lock; // don't modify the above struct without locking here
#endif

void test_differential(const RefSequenceTable& rt, 
					   const vector<HitBundle*>& sample_bundles,
					   const vector<long double>& sample_masses,
					   Tests& tests,
					   Tracking& tracking)
{
	if (sample_bundles.empty())
		return;
	
	string locus_tag = bundle_locus_tag(rt, *(sample_bundles.front()));
	
	bool perform_cds_analysis = true;
	bool perform_tss_analysis = true;
	for (size_t i = 0; i < sample_bundles.size(); ++i)
	{
		foreach(const Scaffold& s, sample_bundles[i]->ref_scaffolds())
		{
			if (s.annotated_tss_id() == "")
			{
				perform_tss_analysis = false;
			}
			if (s.annotated_protein_id() == "")
			{
				perform_cds_analysis = false;
			}
		}
	}
	
	vector<SampleAbundances> samples(sample_bundles.size());

	// set up the transfrag abundance group for each samples
	for (size_t i = 0; i < sample_bundles.size(); ++i)
	{
		samples[i].cluster_mass = sample_bundles[i]->hits().size();
		samples[i].sample_mass = sample_masses[i];
		vector<shared_ptr<Abundance> > abundances;
		
		foreach(const Scaffold& s, sample_bundles[i]->ref_scaffolds())
		{
			TranscriptAbundance* pT = new TranscriptAbundance;
			pT->transfrag(&s);
			shared_ptr<Abundance> ab(pT);
			ab->description(s.annotated_trans_id());
			ab->locus_tag(locus_tag);
			abundances.push_back(ab);
		}
		
		samples[i].transcripts = AbundanceGroup(abundances);
		
		vector<MateHit> hits_in_cluster;
		
		get_alignments_from_scaffolds(samples[i].transcripts.abundances(),
									  hits_in_cluster);
		
		// Compute the individual transcript FPKMs
		samples[i].transcripts.calculate_abundance(hits_in_cluster,
												   samples[i].sample_mass);
		
		// Cluster transcripts by gene_id
		vector<AbundanceGroup> transcripts_by_gene_id;
		cluster_transcripts<ConnectByAnnotatedGeneId>(samples[i].transcripts,
													  transcripts_by_gene_id);
		foreach(AbundanceGroup& ab_group, transcripts_by_gene_id)
		{
			ab_group.locus_tag(locus_tag);
			set<string> gene_ids = ab_group.gene_id();
			assert (gene_ids.size() == 1);
			ab_group.description(*(gene_ids.begin()));
		}
		
		samples[i].genes = transcripts_by_gene_id;
		
		if (perform_cds_analysis)
		{
			// Cluster transcripts by CDS
			vector<AbundanceGroup> transcripts_by_cds;
			ublas::matrix<double> cds_gamma_cov;
			cluster_transcripts<ConnectByAnnotatedProteinId>(samples[i].transcripts,
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
			
			samples[i].cds = transcripts_by_cds;
			
			// Group the CDS clusters by gene
			vector<shared_ptr<Abundance> > cds_abundances;
			foreach (AbundanceGroup& ab_group, samples[i].cds)
			{
				cds_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
			}
			AbundanceGroup cds(cds_abundances,
							   cds_gamma_cov, 
							   samples[i].sample_mass);
			
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
			
			samples[i].gene_cds = cds_by_gene;
		}
		
		if (perform_tss_analysis)
		{
			// Cluster transcripts by start site (TSS)
			vector<AbundanceGroup> transcripts_by_tss;
			
			ublas::matrix<double> tss_gamma_cov;
			cluster_transcripts<ConnectByAnnotatedTssId>(samples[i].transcripts,
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
			
			samples[i].primary_transcripts = transcripts_by_tss;
			
			// Group TSS clusters by gene
			vector<shared_ptr<Abundance> > primary_transcript_abundances;
			foreach (AbundanceGroup& ab_group, samples[i].primary_transcripts)
			{
				primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
			}
			AbundanceGroup primary_transcripts(primary_transcript_abundances,
											   tss_gamma_cov,
											   samples[i].sample_mass);
			
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
			
			samples[i].gene_primary_transcripts = primary_transcripts_by_gene;
		}
	}

#if ENABLE_THREADS
	test_storage_lock.lock();
#endif
	
	for (size_t i = 0; i < samples.size(); ++i)
	{
		const AbundanceGroup& ab_group = samples[i].transcripts;
		samples[i].sample_mass = sample_masses[i];
		double log_mass_curr = log(sample_masses[i]);
		foreach (shared_ptr<Abundance> ab, ab_group.abundances())
		{
			add_to_tracking_table(*ab, log_mass_curr, tracking.isoform_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i].cds)
		{
			add_to_tracking_table(ab, log_mass_curr, tracking.cds_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i].primary_transcripts)
		{
			add_to_tracking_table(ab, log_mass_curr, tracking.tss_group_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, samples[i].genes)
		{
			add_to_tracking_table(ab, log_mass_curr, tracking.gene_fpkm_tracking);
		}
	}
	
	for (size_t i = 1; i < samples.size(); ++i)
	{
		bool multi_transcript_locus = samples[i].transcripts.abundances().size() > 1;
		bool enough_reads = !(multi_transcript_locus &&
							 (samples[i].cluster_mass < min_read_count ||
							  samples[i-1].cluster_mass < min_read_count));
		
		assert (samples[i].transcripts.abundances().size() == 
				samples[i-1].transcripts.abundances().size());
		for (size_t j = 0; j < samples[i].transcripts.abundances().size(); ++j)
		{
			get_de_tests(*(samples[i-1].transcripts.abundances()[j]), 
						 sample_masses[i-1],
						 *(samples[i].transcripts.abundances()[j]),
						 sample_masses[i],
						 tests.isoform_de_tests[i-1],
						 enough_reads);
		}
		
		for (size_t j = 0; j < samples[i].cds.size(); ++j)
		{
			get_de_tests(samples[i-1].cds[j], 
						 sample_masses[i-1],
						 samples[i].cds[j],
						 sample_masses[i],
						 tests.cds_de_tests[i-1],
						 enough_reads);
		}
		
		for (size_t j = 0; j < samples[i].primary_transcripts.size(); ++j)
		{
			get_de_tests(samples[i-1].primary_transcripts[j], 
						 sample_masses[i-1],
						 samples[i].primary_transcripts[j],
						 sample_masses[i],
						 tests.tss_group_de_tests[i-1],
						 enough_reads);
		}
		
		for (size_t j = 0; j < samples[i].genes.size(); ++j)
		{
			get_de_tests(samples[i-1].genes[j], 
						 sample_masses[i-1],
						 samples[i].genes[j],
						 sample_masses[i],
						 tests.gene_de_tests[i-1],
						 enough_reads);
		}
		
		// Differential promoter use
		for (size_t j = 0; j < samples[i].gene_primary_transcripts.size(); ++j)
		{
			get_ds_tests(samples[i-1].gene_primary_transcripts[j], 
						 samples[i].gene_primary_transcripts[j],
						 tests.diff_promoter_tests[i-1],
						 enough_reads);
		}
		
		// Differential coding sequence output
		for (size_t j = 0; j < samples[i].gene_cds.size(); ++j)
		{
			get_ds_tests(samples[i-1].gene_cds[j], 
						 samples[i].gene_cds[j],
						 tests.diff_cds_tests[i-1],
						 enough_reads);
		}
		
		// Differential splicing of primary transcripts
		for (size_t j = 0; j < samples[i].primary_transcripts.size(); ++j)
		{
			get_ds_tests(samples[i-1].primary_transcripts[j], 
						 samples[i].primary_transcripts[j],
						 tests.diff_splicing_tests[i-1],
						 enough_reads);
		}
	}
	
#if ENABLE_THREADS
	test_storage_lock.unlock();
#endif
}
