/*
 *  cuffdiff.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 10/21/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif


#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <numeric>
#include <cfloat>
#include <iostream>

#include "common.h"
#include "hits.h"
#include "bundles.h"
#include "abundances.h"
#include "tokenize.h"
#include "biascorrection.h"
#include "update_check.h"

#include <boost/thread.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "gtf_tracking.h"
#include "differential.h"

// Need at least this many reads in a locus to do any testing on it

vector<string> sample_labels;

double FDR = 0.05; 
bool samples_are_time_series = false;
using namespace std;
using namespace boost;

// We leave out the short codes for options that don't take an argument
#if ENABLE_THREADS
const char *short_options = "m:p:s:F:c:I:j:Q:L:M:o:r:TNqv";
#else
const char *short_options = "m:s:F:c:I:j:Q:L:M:o:r:TNqv";
#endif



static struct option long_options[] = {
{"frag-len-mean",			required_argument,       0,          'm'},
{"frag-len-std-dev",			required_argument,       0,          's'},
{"transcript-score-thresh", required_argument,       0,          't'},
{"min-isoform-fraction",    required_argument,       0,          'F'},
{"pre-mrna-fraction",		required_argument,		 0,			 'j'},
{"max-intron-length",		required_argument,		 0,			 'I'},
{"min-map-qual",			required_argument,		 0,			 'Q'},
{"labels",					required_argument,		 0,			 'L'},
{"min-alignment-count",     required_argument,		 0,			 'c'},
{"FDR",					    required_argument,		 0,			 OPT_FDR},
{"mask-gtf",                required_argument,		 0,			 'M'},
{"output-dir",			    required_argument,		 0,			 'o'},
{"verbose",			    	no_argument,			 0,			 'v'},
{"quiet",			    	no_argument,			 0,			 'q'},
{"reference-seq",			required_argument,		 0,			 'r'},
{"time-series",             no_argument,             0,			 'T'},
{"quartile-normalization",  no_argument,	 		 0,	         'N'},
#if ENABLE_THREADS
{"num-threads",				required_argument,       0,          'p'},
#endif
{"library-type",		    required_argument,		 0,			 OPT_LIBRARY_TYPE},
{"num-importance-samples",  required_argument,		 0,			 OPT_NUM_IMP_SAMPLES},
{"max-mle-iterations",		 required_argument,		 0,			 OPT_MLE_MAX_ITER},
#if ADAM_MODE
{"bias-mode",		    required_argument,		 0,			 OPT_BIAS_MODE},
#endif
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "cuffdiff v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]\n");
	fprintf(stderr, "   Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam\n");
    fprintf(stderr, "Options:\n\n");
    fprintf(stderr, "  -T/--time-series             treat samples as a time-series                        [ default:  FALSE ]\n");
   	fprintf(stderr, "  -N/--quartile-normalization  use upper-quartile normalization                      [ default:  FALSE ]\n");
	fprintf(stderr, "  -Q/--min-map-qual            ignore alignments with lower than this mapping qual   [ default:      0 ]\n");
	fprintf(stderr, "  -c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   1000 ]\n");
	fprintf(stderr, "  --FDR                        False discovery rate used in testing                  [ default:   0.05 ]\n");
	fprintf(stderr, "  -M/--mask-file               ignore all alignment within transcripts in this file  [ default:   NULL ]\n");
    fprintf(stderr, "  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
	fprintf(stderr, "  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
	fprintf(stderr, "  -o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");
	fprintf(stderr, "  -r/--reference-seq           reference fasta file for sequence bias correction     [ default:   NULL ]\n");
    fprintf(stderr, "  -L/--labels                  comma-separated list of condition labels\n");
#if ENABLE_THREADS
	fprintf(stderr, "  -p/--num-threads             number of threads used during quantification          [ default:      1 ]\n");
#endif
	fprintf(stderr, "\nAdvanced Options:\n\n");
    fprintf(stderr, "  --library-type               Library prep used for input reads                     [ default:  below ]\n");
    fprintf(stderr, "  -m/--frag-len-mean           the average fragment length (use with unpaired reads only)               \n");
	fprintf(stderr, "  -s/--frag-len-std-dev        the fragment length standard deviation  (use with unpaired reads only)   \n");
	fprintf(stderr, "  --num-importance-samples     number of importance samples for MAP restimation      [ default:   1000 ]\n");
	fprintf(stderr, "  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]\n");
    print_library_table();
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    string sample_label_list;
    
    do {
        next_option = getopt_long_only(argc, argv, short_options, long_options, &option_index);
        if (next_option == -1)     /* Done with options. */
            break;
        switch (next_option) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                break;
                
			case 'm':
				user_provided_fld = true;
				def_frag_len_mean = (uint32_t)parseInt(0, "-m/--frag-len-mean arg must be at least 0", print_usage);
				break;
			case 'c':
				min_read_count = (uint32_t)parseInt(0, "-c/--min-alignment-count arg must be at least 0", print_usage);
				break;
			case 's':
				user_provided_fld = true;
				def_frag_len_std_dev = (uint32_t)parseInt(0, "-s/--frag-len-std-dev arg must be at least 0", print_usage);
				break;
			case 'p':
				num_threads = (uint32_t)parseInt(1, "-p/--num-threads arg must be at least 1", print_usage);
				break;
            case 'L':
				sample_label_list = optarg;
				break;
			case OPT_FDR:
				FDR = (double)parseFloat(0.00, 1.00, "--FDR arg must be between 0 and 1", print_usage);
				break;
			case OPT_NUM_IMP_SAMPLES:
				num_importance_samples = parseInt(1, "--num-importance-samples must be at least 1", print_usage);
				break;
			case OPT_MLE_MAX_ITER:
				max_mle_iterations = parseInt(1, "--max-mle-iterations must be at least 1", print_usage);
				break;
			case OPT_BIAS_MODE:
				if (!strcmp(optarg, "site"))
					bias_mode = SITE;
				else if (!strcmp(optarg, "pos"))
					bias_mode = POS;
				else if (!strcmp(optarg, "pos_vlmm"))
					bias_mode = POS_VLMM;
				else if (!strcmp(optarg, "vlmm"))
					bias_mode = VLMM;
				else
				{
					fprintf(stderr, "Unknown bias mode.\n");
					exit(1);
				}
				break;
			case 'Q':
			{
				int min_map_qual = parseInt(0, "-Q/--min-map-qual must be at least 0", print_usage);
				if (min_map_qual > 0)
				{
					long double p = (-1.0 * min_map_qual) / 10.0;
					max_phred_err_prob = pow(10.0L, p);
				}
				else
				{
					max_phred_err_prob = 1.0;
				}
				break;
			}
			case 'M':
			{
				mask_gtf_filename = optarg;
				break;
			}
			case 'v':
			{
				if (cuff_quiet)
				{
					fprintf(stderr, "Warning: Can't be both verbose and quiet!  Setting verbose only.\n");
				}
				cuff_quiet = false;
				cuff_verbose = true;
				break;
			}
			case 'q':
			{
				if (cuff_verbose)
				{
					fprintf(stderr, "Warning: Can't be both verbose and quiet!  Setting quiet only.\n");
				}
				cuff_verbose = false;
				cuff_quiet = true;
				break;
			}
            case 'o':
			{
				output_dir = optarg;
				break;
			}
			case 'r':
			{
				fasta_dir = optarg;
				corr_bias = true;
				break;
            }    
                
            case 'T':
			{
                samples_are_time_series = true;
				break;
            }
            case 'N':
            {
            	use_quartile_norm = true;
            	break;
            }
            case OPT_LIBRARY_TYPE:
			{
				library_type = optarg;
				break;
			}

			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    tokenize(sample_label_list, ",", sample_labels);
    
	allow_junk_filtering = false;
	
	return 0;
}


#if ENABLE_THREADS
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

template<typename T>
string cat_strings(const T& container, const char* delimiter=",")
{
	string cat;
	if (container.empty())
	{
		cat = "";
	}
	else
	{
		typename T::const_iterator itr = container.begin();
		//cat = *(itr);
		for (; itr != container.end(); itr++)
		{
			if (!(*itr).empty()) {
				if (!cat.empty()) cat += delimiter;
				cat += *itr; 
				}
		}
	}

	return cat;
}

void print_tests(FILE* fout,
                 const char* sample_1_label,
                 const char* sample_2_label,
				 const SampleDiffs& de_tests)
{
	for (SampleDiffs::const_iterator itr = de_tests.begin(); 
		 itr != de_tests.end(); 
		 ++itr)
	{
		const SampleDifference& test = itr->second;
		
		string all_gene_names = cat_strings(test.gene_names);
		if (all_gene_names == "")
			all_gene_names = "-";
		
		string all_protein_ids = cat_strings(test.protein_ids);	
		if (all_protein_ids == "")
			all_protein_ids = "-";
		
		fprintf(fout, "%s\t%s\t%s\t%s\t%s", 
                itr->first.c_str(), 
                all_gene_names.c_str(), 
                test.locus_desc.c_str(),
                sample_1_label,
                sample_2_label);
		
		if (test.test_status != FAIL)
		{
			double t = test.test_stat;
			double r1 = test.value_1;
			double r2 = test.value_2;
			double d = test.differential;
			double p = test.p_value;
			const char* sig;
			if (test.significant && test.test_status == OK)
				sig = "yes";
			else
				sig = "no";
			
			const char* status;
			if (test.test_status == OK)
				status = "OK";
			else
				status = "NOTEST";
			
			fprintf(fout, "\t%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%s", status, r1, r2, d, t, p, sig);
			fprintf(fout, "\n");
		}
		else
		{
			fprintf(fout, "\tFAIL\t0.0\t0.0\t0.0\t0.0\t1.0\tno\n");
		}
	}
}

void print_FPKM_tracking(FILE* fout, 
						 const FPKMTrackingTable& tracking)
{
	fprintf(fout,"tracking_id\tclass_code\tnearest_ref_id\tgene_short_name\ttss_id\tlocus");
	FPKMTrackingTable::const_iterator first_itr = tracking.begin();
	if (first_itr != tracking.end())
	{
		const FPKMTracking& track = first_itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
			fprintf(fout, "\t%s_FPKM\t%s_conf_lo\t%s_conf_hi", sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str());
		}
	}
	fprintf(fout, "\n");
	for (FPKMTrackingTable::const_iterator itr = tracking.begin(); itr != tracking.end(); ++itr)
	{
		const string& description = itr->first;
		const FPKMTracking& track = itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		
		string all_gene_names = cat_strings(track.gene_names);
		if (all_gene_names == "")
			all_gene_names = "-";
		
		string all_tss_ids = cat_strings(track.tss_ids);
		if (all_tss_ids == "")
			all_tss_ids = "-";
		
		fprintf(fout, "%s\t%c\t%s\t%s\t%s\t%s", 
				description.c_str(),
				track.classcode ? track.classcode : '-',
				track.ref_match.c_str(),
				all_gene_names.c_str(), 
				all_tss_ids.c_str(),
				track.locus_tag.c_str());
		
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
			double fpkm = fpkms[i].FPKM;
			double std_dev = sqrt(fpkms[i].FPKM_variance);
			double fpkm_conf_hi = fpkm + 2.0 * std_dev;
			double fpkm_conf_lo = max(0.0, fpkm - 2.0 * std_dev);
			fprintf(fout, "\t%lg\t%lg\t%lg", fpkm, fpkm_conf_lo, fpkm_conf_hi);
		}
		
		fprintf(fout, "\n");
	}
}

bool p_value_lt(const SampleDifference* lhs, const SampleDifference* rhs)
{
	return lhs->p_value < rhs->p_value;
}

// Benjamani-Hochberg procedure
int fdr_significance(double fdr, 
					  vector<SampleDifference*>& tests)
{
	sort(tests.begin(), tests.end(), p_value_lt);
	vector<SampleDifference*> passing;

	for (int k = 0; k < (int)tests.size(); ++k)
	{
		if (tests[k]->test_status == OK)
		{
			passing.push_back(tests[k]);
		}
		else
		{
			tests[k]->significant = false;
		}
	}
	
	for (int k = 0; k < (int)passing.size(); ++k)
	{
		double r = (double)passing.size() / (k + 1);
		double corrected_p = passing[k]->p_value * r;
		passing[k]->significant = (corrected_p <= fdr);
	}
	return passing.size();
}

void extract_sample_diffs(SampleDiffs& diff_map,
						  vector<SampleDifference*>& diffs)
{
	for (SampleDiffs::iterator itr = diff_map.begin();
		 itr != diff_map.end();
		 ++itr)
	{
		diffs.push_back(&(itr->second));
	}
}

#if ENABLE_THREADS
mutex inspect_lock;
#endif

void inspect_map_worker(ReplicatedBundleFactory& fac,
                        int& tmp_min_frag_len, 
                        int& tmp_max_frag_len)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    int min_f = std::numeric_limits<int>::max();
    int max_f = 0;
    
    fac.inspect_replicate_maps(min_f, max_f);
    
#if ENABLE_THREADS
    inspect_lock.lock();
#endif
    tmp_min_frag_len = min(min_f, tmp_min_frag_len);
    tmp_max_frag_len = max(max_f, tmp_max_frag_len);
#if ENABLE_THREADS
    inspect_lock.unlock();
#endif
}

void learn_bias_worker(shared_ptr<BundleFactory> fac)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
	shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
	BiasLearner* bl = new BiasLearner(rg_props->frag_len_dist());
	learn_bias(*fac, *bl, false);
	rg_props->bias_learner(shared_ptr<BiasLearner const>(bl));
	fac->reset();
}

bool quantitate_next_locus(const RefSequenceTable& rt,
                           vector<ReplicatedBundleFactory>& bundle_factories,
                           vector<shared_ptr<SampleAbundances> >& abundances)
{
    vector<shared_ptr<bool> > non_empty_bundle_flags;
    
    for (size_t i = 0; i < bundle_factories.size(); ++i)
    {
        shared_ptr<bool> sample_non_empty = shared_ptr<bool>(new bool);
        shared_ptr<SampleAbundances> s_ab = shared_ptr<SampleAbundances>(new SampleAbundances);
#if ENABLE_THREADS					
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                locus_thread_pool_lock.unlock();
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
            
        }
        
        locus_thread_pool_lock.lock();
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        thread quantitate(sample_worker,
                          boost::ref(rt),
                          boost::ref(bundle_factories[i]),
                          s_ab,
                          sample_non_empty);  
#else
        sample_worker(boost::ref(rt),
                      boost::ref(bundle_factories[i]),
                      s_ab,
                      sample_non_empty);
#endif
		
        abundances.push_back(s_ab);
        non_empty_bundle_flags.push_back(sample_non_empty);
    }
    
    // wait for the workers to finish up before doing the cross-sample testing.
#if ENABLE_THREADS	
    while(1)
    {
        locus_thread_pool_lock.lock();
        if (locus_curr_threads == 0)
        {
            locus_thread_pool_lock.unlock();
            break;
        }
        
        locus_thread_pool_lock.unlock();
        
        boost::this_thread::sleep(boost::posix_time::milliseconds(5));
        
    }
#endif
    
    bool more_loci_remain = false;
    foreach (shared_ptr<bool> sample_flag, non_empty_bundle_flags)
    {
        if (*sample_flag)
        {
            more_loci_remain = true;
            break;
        }
    }
    
    // There's no more loci to be processed.
    return more_loci_remain;

}

void driver(FILE* ref_gtf, FILE* mask_gtf, vector<string>& sam_hit_filename_lists, Outfiles& outfiles)
{
	check_version(PACKAGE_VERSION);

	ReadTable it;
	RefSequenceTable rt(true, false);
    
	vector<shared_ptr<Scaffold> > ref_mRNAs;
	
	vector<ReplicatedBundleFactory> bundle_factories;    
    vector<shared_ptr<HitFactory> > all_hit_factories;
    
	for (size_t i = 0; i < sam_hit_filename_lists.size(); ++i)
	{
        vector<string> sam_hit_filenames;
        tokenize(sam_hit_filename_lists[i], ",", sam_hit_filenames);
        
        vector<shared_ptr<BundleFactory> > replicate_factories;
        for (size_t j = 0; j < sam_hit_filenames.size(); ++j)
        {
            shared_ptr<HitFactory> hs;
            try
            {
                hs = shared_ptr<HitFactory>(new BAMHitFactory(sam_hit_filenames[j], it, rt));
            }
            catch (std::runtime_error& e) 
            {
                try
                {
                    fprintf(stderr, "File %s doesn't appear to be a valid BAM file, trying SAM...\n",
                            sam_hit_filename_lists[i].c_str());
                    hs = shared_ptr<HitFactory>(new SAMHitFactory(sam_hit_filenames[j], it, rt));
                }
                catch (std::runtime_error& e)
                {
                    fprintf(stderr, "Error: cannot open alignment file %s for reading\n",
                            sam_hit_filenames[j].c_str());
                    exit(1);
                }
            }
            
            all_hit_factories.push_back(hs);
            
            shared_ptr<BundleFactory> hf(new BundleFactory(hs, REF_DRIVEN));
            shared_ptr<ReadGroupProperties> rg_props(new ReadGroupProperties);
            
            if (global_read_properties)
            {
                *rg_props = *global_read_properties;
            }
            else 
            {
                *rg_props = hs->read_group_properties();
            }
            
            hf->read_group_properties(rg_props);
            
            replicate_factories.push_back(hf);
            //replicate_factories.back()->set_ref_rnas(ref_mRNAs);
        }
        
        ReplicatedBundleFactory rep_factory(replicate_factories);
        bundle_factories.push_back(rep_factory);
	}
    
    ::load_ref_rnas(ref_gtf, rt, ref_mRNAs, corr_bias, false);
    if (ref_mRNAs.empty())
        return;
    
    vector<shared_ptr<Scaffold> > mask_rnas;
    if (mask_gtf)
    {
        ::load_ref_rnas(mask_gtf, rt, mask_rnas, false, false);
    }
    
    foreach (ReplicatedBundleFactory& fac, bundle_factories)
    {
        fac.set_ref_rnas(ref_mRNAs);
        if (mask_gtf) fac.set_mask_rnas(mask_rnas);
    }
    
#if ENABLE_THREADS
    locus_num_threads = num_threads;
#endif
    
	int tmp_min_frag_len = numeric_limits<int>::max();
	int tmp_max_frag_len = 0;
	
	ProgressBar p_bar("Inspecting maps and determining fragment length distributions.",0);
	foreach (ReplicatedBundleFactory& fac, bundle_factories)
    {
#if ENABLE_THREADS	
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                locus_thread_pool_lock.unlock();
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
        }
        locus_thread_pool_lock.lock();
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        thread inspect(inspect_map_worker,
                       boost::ref(fac),
                       boost::ref(tmp_min_frag_len),
                       boost::ref(tmp_max_frag_len));  
#else
        inspect_map_worker(boost::ref(fac),
                           boost::ref(tmp_min_frag_len),
                           boost::ref(tmp_max_frag_len));
#endif
    }
    
    // wait for the workers to finish up before reporting everthing.
#if ENABLE_THREADS	
    while(1)
    {
        locus_thread_pool_lock.lock();
        if (locus_curr_threads == 0)
        {
            locus_thread_pool_lock.unlock();
            break;
        }
        
        locus_thread_pool_lock.unlock();
        
        boost::this_thread::sleep(boost::posix_time::milliseconds(5));
    }
#endif

	min_frag_len = tmp_min_frag_len;
    max_frag_len = tmp_max_frag_len;
	
	final_est_run = false;
	
	int num_bundles = bundle_factories[0].num_bundles();
	
	if (corr_bias) // Only run initial estimation if correcting bias
	{
		p_bar = ProgressBar("Calculating initial abundance estimates.", num_bundles);
		
		while (1) 
		{
			p_bar.update("",1);
			vector<shared_ptr<SampleAbundances> > abundances;
			bool more_loci_remain = quantitate_next_locus(rt, bundle_factories, abundances);
			
			if (!more_loci_remain)
				break;
		}
		p_bar.complete();
	
		p_bar = ProgressBar("Learning bias parameters.", 0);
		foreach (ReplicatedBundleFactory& rep_fac, bundle_factories)
		{
			rep_fac.reset();
			foreach (shared_ptr<BundleFactory> fac, rep_fac.factories())
			{
#if ENABLE_THREADS	
				while(1)
				{
					locus_thread_pool_lock.lock();
					if (locus_curr_threads < locus_num_threads)
					{
						locus_thread_pool_lock.unlock();
						break;
					}
					
					locus_thread_pool_lock.unlock();
					
					boost::this_thread::sleep(boost::posix_time::milliseconds(5));
				}
				locus_thread_pool_lock.lock();
				locus_curr_threads++;
				locus_thread_pool_lock.unlock();
				
				thread bias(learn_bias_worker, fac);
#else
				learn_bias_worker(fac);
#endif
			}
    	}
    
    // wait for the workers to finish up before reporting everthing.
#if ENABLE_THREADS	
		while(1)
		{
			locus_thread_pool_lock.lock();
			if (locus_curr_threads == 0)
			{
				locus_thread_pool_lock.unlock();
				break;
			}
			
			locus_thread_pool_lock.unlock();
			
			boost::this_thread::sleep(boost::posix_time::milliseconds(5));
		}
#endif
	}
	
	Tests tests;
    
    int N = (int)sam_hit_filename_lists.size();
    
    tests.isoform_de_tests = vector<vector<SampleDiffs> >(N);
    tests.tss_group_de_tests = vector<vector<SampleDiffs> >(N);
    tests.gene_de_tests = vector<vector<SampleDiffs> >(N);
    tests.cds_de_tests = vector<vector<SampleDiffs> >(N);
    tests.diff_splicing_tests = vector<vector<SampleDiffs> >(N);
    tests.diff_promoter_tests = vector<vector<SampleDiffs> >(N);
    tests.diff_cds_tests = vector<vector<SampleDiffs> >(N);
    
	for (int i = 1; i < N; ++i)
    {
        tests.isoform_de_tests[i] = vector<SampleDiffs>(i);
        tests.tss_group_de_tests[i] = vector<SampleDiffs>(i);
        tests.gene_de_tests[i] = vector<SampleDiffs>(i);
        tests.cds_de_tests[i] = vector<SampleDiffs>(i);
        tests.diff_splicing_tests[i] = vector<SampleDiffs>(i);
        tests.diff_promoter_tests[i] = vector<SampleDiffs>(i);
        tests.diff_cds_tests[i] = vector<SampleDiffs>(i);
    }

	Tracking tracking;
	
	final_est_run = true;
	p_bar = ProgressBar("Testing for differential expression and regulation in locus.", num_bundles);
	while (true)
	{
        vector<shared_ptr<SampleAbundances> > abundances;
        bool more_loci_remain = quantitate_next_locus(rt, bundle_factories, abundances);
        
        if (!more_loci_remain)
            break;
        
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
                assert (s1.abundances()[j]->description() == s1.abundances()[j]->description());
            }
        }

        verbose_msg("Testing for differential expression and regulation in locus [%s]\n", abundances.front()->locus_tag.c_str());
		p_bar.update(abundances.front()->locus_tag.c_str(), 1);
		test_differential(abundances.front()->locus_tag, abundances, tests, tracking, samples_are_time_series);
	}
	
	p_bar.complete();

	//double FDR = 0.05;
	int total_iso_de_tests = 0;
	
	vector<SampleDifference*> isoform_exp_diffs;
	for (size_t i = 1; i < tests.isoform_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            total_iso_de_tests += tests.isoform_de_tests[i][j].size();
            extract_sample_diffs(tests.isoform_de_tests[i][j], isoform_exp_diffs);
        }
	}
	int iso_exp_tests = fdr_significance(FDR, isoform_exp_diffs);
	fprintf(stderr, "Performed %d isoform-level transcription difference tests\n", iso_exp_tests);
    fprintf(outfiles.isoform_de_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tln(fold_change)\ttest_stat\tp_value\tsignificant\n");
	for (size_t i = 1; i < tests.isoform_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.isoform_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.isoform_de_tests[i][j]);
        }
	}
	
	int total_group_de_tests = 0;
	vector<SampleDifference*> tss_group_exp_diffs;
	for (size_t i = 1; i < tests.tss_group_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            extract_sample_diffs(tests.tss_group_de_tests[i][j], tss_group_exp_diffs);
            total_group_de_tests += tests.tss_group_de_tests[i][j].size();
        }
	}
	
	int tss_group_exp_tests = fdr_significance(FDR, tss_group_exp_diffs);
	fprintf(stderr, "Performed %d tss-level transcription difference tests\n", tss_group_exp_tests);
    fprintf(outfiles.group_de_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tln(fold_change)\ttest_stat\tp_value\tsignificant\n");
	for (size_t i = 1; i < tests.tss_group_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.group_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.tss_group_de_tests[i][j]);
        }
	}
	
	int total_gene_de_tests = 0;
	vector<SampleDifference*> gene_exp_diffs;
	for (size_t i = 1; i < tests.gene_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            total_gene_de_tests += tests.gene_de_tests[i][j].size();
            extract_sample_diffs(tests.gene_de_tests[i][j], gene_exp_diffs);
        }
	}
	
	int gene_exp_tests = fdr_significance(FDR, gene_exp_diffs);
	fprintf(stderr, "Performed %d gene-level transcription difference tests\n", gene_exp_tests);
	fprintf(outfiles.gene_de_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tln(fold_change)\ttest_stat\tp_value\tsignificant\n");
    for (size_t i = 1; i < tests.gene_de_tests.size(); ++i)
	{        
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.gene_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.gene_de_tests[i][j]);
        }
	}	

	int total_cds_de_tests = 0;
	vector<SampleDifference*> cds_exp_diffs;
	for (size_t i = 1; i < tests.cds_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            total_cds_de_tests += tests.cds_de_tests[i][j].size();
            extract_sample_diffs(tests.cds_de_tests[i][j], cds_exp_diffs);
        }
	}
	int cds_exp_tests = fdr_significance(FDR, cds_exp_diffs);
	fprintf(stderr, "Performed %d CDS-level transcription difference tests\n", cds_exp_tests);
	fprintf(outfiles.cds_de_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tln(fold_change)\ttest_stat\tp_value\tsignificant\n");
    for (size_t i = 1; i < tests.cds_de_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.cds_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.cds_de_tests[i][j]);
        }
	}
	
	int total_diff_splice_tests = 0;
	vector<SampleDifference*> splicing_diffs;
	for (size_t i = 1; i < tests.diff_splicing_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            total_diff_splice_tests += tests.diff_splicing_tests[i][j].size();
            extract_sample_diffs(tests.diff_splicing_tests[i][j], splicing_diffs);
        }
	}
	
	int splicing_tests = fdr_significance(FDR, splicing_diffs);
	fprintf(stderr, "Performed %d splicing tests\n", splicing_tests);
	fprintf(outfiles.diff_splicing_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tsignificant\n");
    for (size_t i = 1; i < tests.diff_splicing_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            const SampleDiffs& diffs = tests.diff_splicing_tests[i][j];
            print_tests(outfiles.diff_splicing_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), diffs);
        }
	}
	
	int total_diff_promoter_tests = 0;
	vector<SampleDifference*> promoter_diffs;
	for (size_t i = 1; i < tests.diff_splicing_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            total_diff_promoter_tests += tests.diff_promoter_tests[i][j].size();
            extract_sample_diffs(tests.diff_promoter_tests[i][j], promoter_diffs);
        }
	}
	int promoter_tests = fdr_significance(FDR, promoter_diffs);
	fprintf(stderr, "Performed %d promoter preference tests\n", promoter_tests);
    fprintf(outfiles.diff_promoter_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tsignificant\n");
    for (size_t i = 1; i < tests.diff_promoter_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.diff_promoter_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.diff_promoter_tests[i][j]);
        }
	}

	int total_diff_cds_tests = 0;
	vector<SampleDifference*> cds_use_diffs;
	for (size_t i = 1; i < tests.diff_cds_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            extract_sample_diffs(tests.diff_cds_tests[i][j], cds_use_diffs);
            total_diff_cds_tests += tests.diff_cds_tests[i][j].size();
        }
	}
	int cds_use_tests = fdr_significance(FDR, cds_use_diffs);
	fprintf(stderr, "Performing %d relative CDS output tests\n", cds_use_tests);
	fprintf(outfiles.diff_cds_outfile, "test_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tsignificant\n");
    for (size_t i = 1; i < tests.diff_cds_tests.size(); ++i)
	{
        for (size_t j = 0; j < i; ++j)
        {
            print_tests(outfiles.diff_cds_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.diff_cds_tests[i][j]);
        }
	}
	
	FILE* fiso_fpkm_tracking =  outfiles.isoform_fpkm_tracking_out;
	fprintf(stderr, "Writing isoform-level FPKM tracking\n");
	print_FPKM_tracking(fiso_fpkm_tracking,tracking.isoform_fpkm_tracking); 
	
	FILE* ftss_fpkm_tracking =  outfiles.tss_group_fpkm_tracking_out;
	fprintf(stderr, "Writing TSS group-level FPKM tracking\n");
	print_FPKM_tracking(ftss_fpkm_tracking,tracking.tss_group_fpkm_tracking);
	
	FILE* fgene_fpkm_tracking =  outfiles.gene_fpkm_tracking_out;
	fprintf(stderr, "Writing gene-level FPKM tracking\n");
	print_FPKM_tracking(fgene_fpkm_tracking,tracking.gene_fpkm_tracking);
	
	FILE* fcds_fpkm_tracking =  outfiles.cds_fpkm_tracking_out;
	fprintf(stderr, "Writing CDS-level FPKM tracking\n");
	print_FPKM_tracking(fcds_fpkm_tracking,tracking.cds_fpkm_tracking);
}

int main(int argc, char** argv)
{
    init_library_table();
    
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string ref_gtf_filename = argv[optind++];
	
	vector<string> sam_hit_filenames;
    while(optind < argc)
    {
        string sam_hits_file_name = argv[optind++];
		sam_hit_filenames.push_back(sam_hits_file_name);
    }
    	
	while (sam_hit_filenames.size() < 2)
    {
        fprintf(stderr, "Error: cuffdiff requires at least 2 SAM files\n");
        exit(1);
    }
	
    if (sample_labels.size() == 0)
    {
        for (size_t i = 1; i < sam_hit_filenames.size() + 1; ++i)
        {
            char buf[256];
            sprintf(buf, "q%lu", i);
            sample_labels.push_back(buf);
        }   
    }
    
    if (sam_hit_filenames.size() != sample_labels.size())
    {
        fprintf(stderr, "Error: number of labels must match number of conditions\n");
        exit(1);
    }
    
	// seed the random number generator - we'll need it for the importance
	// sampling during MAP estimation of the gammas
	srand48(time(NULL));
	
	FILE* ref_gtf = NULL;
	if (ref_gtf_filename != "")
	{
		ref_gtf = fopen(ref_gtf_filename.c_str(), "r");
		if (!ref_gtf)
		{
			fprintf(stderr, "Error: cannot open reference GTF file %s for reading\n",
					ref_gtf_filename.c_str());
			exit(1);
		}
	}
	
	FILE* mask_gtf = NULL;
	if (mask_gtf_filename != "")
	{
		mask_gtf = fopen(mask_gtf_filename.c_str(), "r");
		if (!mask_gtf)
		{
			fprintf(stderr, "Error: cannot open mask GTF file %s for reading\n",
					mask_gtf_filename.c_str());
			exit(1);
		}
	}
	
    
	
	// Note: we don't want the assembly filters interfering with calculations 
	// here
	min_isoform_fraction = 0.0;
	pre_mrna_fraction = 0.0;
    olap_radius = 0;
	
	Outfiles outfiles;
	
    if (output_dir != "")
    {
        int retcode = mkpath(output_dir.c_str(), 0777);
        if (retcode == -1)
        {
            if (errno != EEXIST)
            {
                fprintf (stderr, 
                         "Error: cannot create directory %s\n", 
                         output_dir.c_str());
                exit(1);
            }
        }
    }
    
    char out_file_prefix[64];
    sprintf(out_file_prefix, "%s/", output_dir.c_str());
    char iso_out_file_name[256];
    sprintf(iso_out_file_name, "%sisoform_exp.diff", out_file_prefix);
    FILE* iso_out = fopen(iso_out_file_name, "w");
    if (!iso_out)
    {
        fprintf(stderr, "Error: cannot open differential isoform transcription file %s for writing\n",
                iso_out_file_name);
        exit(1);
    }
    
    char group_out_file_name[256];
    sprintf(group_out_file_name, "%stss_group_exp.diff", out_file_prefix);
    FILE* group_out = fopen(group_out_file_name, "w");
    if (!group_out)
    {
        fprintf(stderr, "Error: cannot open differential TSS group transcription file %s for writing\n",
                group_out_file_name);
        exit(1);
    }
    
    char gene_out_file_name[256];
    sprintf(gene_out_file_name, "%sgene_exp.diff", out_file_prefix);
    FILE* gene_out = fopen(gene_out_file_name, "w");
    if (!group_out)
    {
        fprintf(stderr, "Error: cannot open gene expression file %s for writing\n",
                gene_out_file_name);
        exit(1);
    }
    
    char cds_out_file_name[256];
    sprintf(cds_out_file_name, "%scds_exp.diff", out_file_prefix);
    FILE* cds_out = fopen(cds_out_file_name, "w");
    if (!cds_out)
    {
        fprintf(stderr, "Error: cannot open cds expression file %s for writing\n",
                cds_out_file_name);
        exit(1);
    }
    
    char diff_splicing_out_file_name[256];
    sprintf(diff_splicing_out_file_name, "%ssplicing.diff", out_file_prefix);
    FILE* diff_splicing_out = fopen(diff_splicing_out_file_name, "w");
    if (!diff_splicing_out)
    {
        fprintf(stderr, "Error: cannot open differential splicing file %s for writing\n",
                diff_splicing_out_file_name);
        exit(1);
    }
    
    char diff_promoter_out_file_name[256];
    sprintf(diff_promoter_out_file_name, "%spromoters.diff", out_file_prefix);
    FILE* diff_promoter_out = fopen(diff_promoter_out_file_name, "w");
    if (!diff_promoter_out)
    {
        fprintf(stderr, "Error: cannot open differential transcription start file %s for writing\n",
                diff_promoter_out_file_name);
        exit(1);
    }
    
    char diff_cds_out_file_name[256];
    sprintf(diff_cds_out_file_name, "%scds.diff", out_file_prefix);
    FILE* diff_cds_out = fopen(diff_cds_out_file_name, "w");
    if (!diff_cds_out)
    {
        fprintf(stderr, "Error: cannot open differential relative CDS file %s for writing\n",
                diff_cds_out_file_name);
        exit(1);
    }
    
    outfiles.isoform_de_outfile = iso_out;
    outfiles.group_de_outfile = group_out;
    outfiles.gene_de_outfile = gene_out;
    outfiles.cds_de_outfile = cds_out;
    outfiles.diff_splicing_outfile = diff_splicing_out;
    outfiles.diff_promoter_outfile = diff_promoter_out;
    outfiles.diff_cds_outfile = diff_cds_out;
	
	char isoform_fpkm_tracking_name[256];
	sprintf(isoform_fpkm_tracking_name, "%s/isoforms.fpkm_tracking", output_dir.c_str());
	FILE* isoform_fpkm_out = fopen(isoform_fpkm_tracking_name, "w");
	if (!isoform_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level FPKM tracking file %s for writing\n",
				isoform_fpkm_tracking_name);
		exit(1);
	}
	outfiles.isoform_fpkm_tracking_out = isoform_fpkm_out;

	char tss_group_fpkm_tracking_name[256];
	sprintf(tss_group_fpkm_tracking_name, "%s/tss_groups.fpkm_tracking", output_dir.c_str());
	FILE* tss_group_fpkm_out = fopen(tss_group_fpkm_tracking_name, "w");
	if (!tss_group_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level FPKM tracking file %s for writing\n",
				tss_group_fpkm_tracking_name);
		exit(1);
	}
	outfiles.tss_group_fpkm_tracking_out = tss_group_fpkm_out;

	char cds_fpkm_tracking_name[256];
	sprintf(cds_fpkm_tracking_name, "%s/cds.fpkm_tracking", output_dir.c_str());
	FILE* cds_fpkm_out = fopen(cds_fpkm_tracking_name, "w");
	if (!cds_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open CDS level FPKM tracking file %s for writing\n",
				cds_fpkm_tracking_name);
		exit(1);
	}
	outfiles.cds_fpkm_tracking_out = cds_fpkm_out;
	
	char gene_fpkm_tracking_name[256];
	sprintf(gene_fpkm_tracking_name, "%s/genes.fpkm_tracking", output_dir.c_str());
	FILE* gene_fpkm_out = fopen(gene_fpkm_tracking_name, "w");
	if (!gene_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open gene-level FPKM tracking file %s for writing\n",
				gene_fpkm_tracking_name);
		exit(1);
	}
	outfiles.gene_fpkm_tracking_out = gene_fpkm_out;
	
    driver(ref_gtf, mask_gtf, sam_hit_filenames, outfiles);
	
	return 0;
}

