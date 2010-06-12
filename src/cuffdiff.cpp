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

double FDR = 0.05; 
using namespace std;
using namespace boost;

#if ENABLE_THREADS
const char *short_options = "m:p:s:F:c:I:j:Q:L:G:";
#else
const char *short_options = "m:s:F:c:I:j:Q:L:G:";
#endif

static struct option long_options[] = {
{"inner-dist-mean",			required_argument,       0,          'm'},
{"inner-dist-stddev",		required_argument,       0,          's'},
{"transcript-score-thresh", required_argument,       0,          't'},
{"min-isoform-fraction",    required_argument,       0,          'F'},
{"pre-mrna-fraction",		required_argument,		 0,			 'j'},
{"max-intron-length",		required_argument,		 0,			 'I'},
{"min-map-qual",			required_argument,		 0,			 'Q'},
{"label",					required_argument,		 0,			 'L'},
{"min-alignment-count",     required_argument,		 0,			 'c'},
{"GTF",					    required_argument,		 0,			 'G'},
{"FDR",					    required_argument,		 0,			 OPT_FDR},
{"output-dir",			    required_argument,		 0,			 'o'},
#if ENABLE_THREADS
{"num-threads",				required_argument,       0,          'p'},
#endif
{"num-importance-samples",  required_argument,		 0,			 OPT_NUM_IMP_SAMPLES},
{"max-mle-iterations",		 required_argument,		 0,			 OPT_MLE_MAX_ITER},
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "cuffdiff v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-m/--inner-dist-mean         the average inner distance between mates              [ default:     45 ]\n");
	fprintf(stderr, "-s/--inner-dist-std-dev      the inner distance standard deviation                 [ default:     20 ]\n");
	fprintf(stderr, "-Q/--min-map-qual            ignore alignments with lower than this mapping qual   [ default:      0 ]\n");
	fprintf(stderr, "-c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   1000 ]\n");
	fprintf(stderr, "--FDR						  False discovery rate used in testing   [ default:   0.05 ]\n");
	fprintf(stderr, "-o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");

#if ENABLE_THREADS
	fprintf(stderr, "-p/--num-threads             number of threads used during assembly                [ default:      1 ]\n");
#endif
	fprintf(stderr, "\nAdvanced Options:\n\n");
	fprintf(stderr, "--num-importance-samples     number of importance samples for MAP restimation      [ default:   1000 ]\n");
	fprintf(stderr, "--max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]\n");
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'm':
				inner_dist_mean = (uint32_t)parseInt(-200, "-m/--inner-dist-mean arg must be at least -200", print_usage);
				break;
			case 'c':
				min_read_count = (uint32_t)parseInt(0, "-c/--min-alignment-count arg must be at least 0", print_usage);
				break;
			case 's':
				inner_dist_std_dev = (uint32_t)parseInt(0, "-s/--inner-dist-std-dev arg must be at least 0", print_usage);
				break;
			case 'p':
				num_threads = (uint32_t)parseInt(1, "-p/--num-threads arg must be at least 1", print_usage);
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
			case 'G':
			{
				ref_gtf_filename = optarg;
				break;
			}
            case 'o':
			{
				output_dir = optarg;
				break;
			}
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	max_inner_dist = inner_dist_mean + 3 * inner_dist_std_dev;
	inner_dist_norm = normal(inner_dist_mean, inner_dist_std_dev);
	allow_junk_filtering = false;
	
	return 0;
}

// like a normal bundle factory, except the user must explicitly ask for 
// a range of alignments for the bunde.  NOTE: this factory does NOT seek 
// backwards - it is up to you to ask for monotonically increasing loci
class LocusBundleFactory
{
public:
	LocusBundleFactory(SAMHitFactory& fac, FILE* hfile)
		: sam_hit_fac(fac), hit_file(hfile), _next_line_num(0) {}
	
	bool next_bundle(HitBundle& bundle_out, 
					 RefID ref_id, 
					 int left_boundary, 
					 int right_boundary);
	
	SAMHitFactory& hit_factory() { return sam_hit_fac; } 
	
	void reset() { rewind(hit_file); _next_line_num = 0; }
	
private:
	SAMHitFactory sam_hit_fac;
	FILE* hit_file;
	int _next_line_num;
};

bool LocusBundleFactory::next_bundle(HitBundle& bundle_out, 
									 RefID ref_id, 
									 int left_boundary, 
									 int right_boundary)
{
	HitBundle bundle = bundle_out;
	
	if (feof(hit_file))
	{
		return false;
	}
	char bwt_buf[2048];
	
	RefID last_ref_id_seen = 0;
	int last_pos_seen = 0;
	
	off_t curr_pos = ftello(hit_file);
	
	
	while (fgets(bwt_buf, 2048, hit_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		_next_line_num++;
		
		shared_ptr<ReadHit> bh(new ReadHit());
		
		if (!sam_hit_fac.get_hit_from_buf(_next_line_num, bwt_buf, *bh, false))
		{
			continue;
		}
		
		if (bh->ref_id() == 84696373) // corresponds to SAM "*" under FNV hash. unaligned read record  
			continue;
		
		bool hit_within_boundary = false;
		
		if (bh->ref_id() != ref_id)
		{
			RefSequenceTable& rt = sam_hit_fac.ref_table();
			const char* gtf_name = rt.get_name(ref_id);
			const char* sam_name = rt.get_name(bh->ref_id());
			if (!gtf_name)
			{
				// The GTF name isn't even inthe SAM file, we'll never find
				// records by advancing the SAM.  Give up.  hit_within_boundary
				// will remain false, and we'll break out of the loop and 
				// reset the SAM file pointer
			}
			else
			{
				// the LocusBundleFactory's ref table should be populated with all 
				// values from the SAM
				assert (sam_name); 
				int c = strcmp(gtf_name, sam_name);
				
				if (c > 0)
				{
					// we need to keep advancing the SAM file, to catch up
					// with the GTF file
					continue;
				}
				else
				{
					// else the GTF file is before the SAM file, we'll never find
					// records by advancing the SAM.  Give up.  hit_within_boundary
					// will remain false, and we'll break out of the loop and 
					// reset the SAM file pointer
				}
			}
		}
		else
		{
			if (bh->right() < left_boundary)
				continue;
			
			if (bh->error_prob() > max_phred_err_prob)
				continue;
			
			if (bh->left() <= right_boundary)
				hit_within_boundary = true;
		}
		
		if (hit_within_boundary)
		{
			if (bh->left() < last_pos_seen)
			{
				fprintf(stderr, "Error: this SAM file doesn't appear to be correctly sorted!\n");
				fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
						sam_hit_fac.ref_table().get_name(bh->ref_id()),
						bh->left(),
						sam_hit_fac.ref_table().get_name(last_ref_id_seen),
						last_pos_seen);
				
				exit(1);
			}
			
			bundle.add_open_hit(bh);
		}
		else
		{
			fseeko(hit_file, curr_pos, SEEK_SET);
			break;
		}
		
		last_ref_id_seen = bh->ref_id();
		last_pos_seen = bh->left();
		
		curr_pos = ftello(hit_file);
	}
	
	bundle.finalize_open_mates();
	bundle_out = bundle;
	bundle_out.finalize();
	assert (!bundle_out.ref_scaffolds().empty());
	return true;
}

void make_sample_bundles(vector<Scaffold>& ref_mrnas,
						 vector<LocusBundleFactory>& bundle_factories,
						 vector<HitBundle*>& locus_bundles)
{
	RefID ref_id = 0;
	int left_boundary = 999999999;
	int right_boundary = -1;
	for (size_t i = 0; i < ref_mrnas.size(); ++i)
	{
		assert (ref_id == 0 || ref_mrnas[i].ref_id());
		ref_id = ref_mrnas[i].ref_id();
		
		if (ref_mrnas[i].left() < left_boundary)
			left_boundary = ref_mrnas[i].left();
		if (ref_mrnas[i].right() > right_boundary)
			right_boundary = ref_mrnas[i].right();
	}
	
	for (size_t i = 0; i < bundle_factories.size(); ++i)
	{
		LocusBundleFactory& fac = bundle_factories[i];
		HitBundle* bundle = new HitBundle;
		for (size_t j = 0; j < ref_mrnas.size(); ++j)
		{
			bundle->add_ref_scaffold(ref_mrnas[j]);
		}
		
		fac.next_bundle(*bundle, ref_id, left_boundary, right_boundary);
		assert (!bundle->ref_scaffolds().empty());
		locus_bundles.push_back(bundle);
	}
}	

#if ENABLE_THREADS

mutex thread_pool_lock;
int curr_threads = 0;

void decr_pool_count()
{
	thread_pool_lock.lock();
	curr_threads--;
	thread_pool_lock.unlock();	
}
#endif

void quantitation_worker(const RefSequenceTable& rt,
						 vector<HitBundle*>* sample_bundles,
						 const vector<long double>& sample_masses,
						 Tests& tests,
						 Tracking& tracking)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif

	test_differential(rt, *sample_bundles, sample_masses, tests, tracking);
	
	for (size_t i = 0; i < sample_bundles->size(); ++i)
	{
		delete (*sample_bundles)[i];
	}
	
	delete sample_bundles;
}

template<typename T>
string cat_strings(const T& container)
{
	string cat;
	if (container.empty())
	{
		cat = "";
	}
	else
	{
		typename T::const_iterator itr = container.begin();
		cat = *(itr);
		for (++itr; itr != container.end(); ++itr)
		{
			cat += "," + *itr;
		}
	}

	return cat;
}

void print_tests(FILE* fout,
				 const char* label,
				 const SampleDiffs& de_tests)
{
	fprintf(fout, "test_id\tgene\tlocus\tstatus\tvalue_1\tvalue_2\t%s\ttest_stat\tp_value\tsignificant\n", label);
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
		
		fprintf(fout, "%s\t%s\t%s", itr->first.c_str(), all_gene_names.c_str(), test.locus_desc.c_str());
		
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
			fprintf(fout, "\tq%lu_FPKM\tq%lu_conf_lo\tq%lu_conf_hi", i, i, i);
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

void driver(FILE* ref_gtf, vector<FILE*>& sam_hit_files, Outfiles& outfiles)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<Scaffold> ref_mRNAs;
	load_ref_rnas(ref_gtf, rt, ref_mRNAs);
	if (ref_mRNAs.empty())
		return;
	
	vector<LocusBundleFactory> bundle_factories;
	vector<long double> map_masses;
	for (size_t i = 0; i < sam_hit_files.size(); ++i)
	{
		SAMHitFactory hs(it, rt);
		LocusBundleFactory lf(hs, sam_hit_files[i]);
		bundle_factories.push_back(lf);
		BundleFactory standard_factory(hs, sam_hit_files[i], NULL);
		
		fprintf(stderr, "Counting hits in sample %lu\n", i);
		long double map_mass = get_map_mass(standard_factory);
		map_masses.push_back(map_mass);
	}
	
	vector<Scaffold>::iterator curr_ref_scaff = ref_mRNAs.begin();
	
	RefID last_ref_seen = curr_ref_scaff->ref_id();
	int left_boundary = curr_ref_scaff->left();
	int right_boundary = curr_ref_scaff->right();
	
	vector<Scaffold>::iterator ref_scaff_bundle_start = curr_ref_scaff;
	Tests tests;
	
	tests.isoform_de_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.tss_group_de_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.gene_de_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.cds_de_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.diff_splicing_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.diff_promoter_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	tests.diff_cds_tests = vector<SampleDiffs>((int)sam_hit_files.size() - 1);
	
	Tracking tracking;
	
	while (curr_ref_scaff != ref_mRNAs.end())
	{
		if (curr_ref_scaff->left() > right_boundary ||
			curr_ref_scaff->ref_id() != last_ref_seen)
		{
#if ENABLE_THREADS			
			while(1)
			{
				thread_pool_lock.lock();
				if (curr_threads < num_threads)
				{
					thread_pool_lock.unlock();
					break;
				}
				
				thread_pool_lock.unlock();
				
				boost::this_thread::sleep(boost::posix_time::milliseconds(5));
				
			}
#endif
			
			vector<Scaffold> ref_mrnas;
			
			ref_mrnas.insert(ref_mrnas.end(), 
							 ref_scaff_bundle_start, 
							 curr_ref_scaff);
			
			vector<HitBundle*>* sample_bundles = new vector<HitBundle*>();
			
			// grab the alignments for this locus in each of the samples
			make_sample_bundles(ref_mrnas, bundle_factories, *sample_bundles);
			bool non_empty_bundle = false;
			for (size_t i = 0; i < sample_bundles->size(); ++i)
			{
				if (!(*sample_bundles)[i]->hits().empty())
				{
					non_empty_bundle = true;
					break;
				}
			}
			if (non_empty_bundle)
			{
				RefID bundle_chr_id = sample_bundles->front()->ref_id();
				assert (bundle_chr_id != 0);
				const char* chr_name = rt.get_name(bundle_chr_id);
				assert (chr_name);
				fprintf(stderr, "Quantitating samples in locus [ %s:%d-%d ] \n", 
						chr_name,
						sample_bundles->front()->left(),
						sample_bundles->front()->right());
			
			
#if ENABLE_THREADS			
				thread_pool_lock.lock();
				curr_threads++;
				thread_pool_lock.unlock();
				
				thread quantitate(quantitation_worker,
								  boost::cref(rt), 
								  sample_bundles, 
								  boost::cref(map_masses), 
								  boost::ref(tests),
								  boost::ref(tracking));

#else
				quantitation_worker(boost::cref(rt), 
									sample_bundles, 
									boost::cref(map_masses), 
									boost::ref(tests),
									boost::ref(tracking));
#endif
			}
			left_boundary = curr_ref_scaff->left();
			if (curr_ref_scaff->ref_id() != last_ref_seen)
			{
				right_boundary = curr_ref_scaff->right();
			}
			last_ref_seen = curr_ref_scaff->ref_id();
			ref_scaff_bundle_start = curr_ref_scaff;
		}
		right_boundary = max(right_boundary, curr_ref_scaff->right());
		++curr_ref_scaff;
	}
	
	if (ref_scaff_bundle_start != ref_mRNAs.end())
	{
		vector<Scaffold> ref_mrnas;
		
		ref_mrnas.insert(ref_mrnas.end(), 
						 ref_scaff_bundle_start, 
						 curr_ref_scaff);
		
		vector<HitBundle*> sample_bundles;
		
		// grab the alignments for this locus in each of the samples
		make_sample_bundles(ref_mrnas, bundle_factories, sample_bundles);
		bool non_empty_bundle = false;
		for (size_t i = 0; i < sample_bundles.size(); ++i)
		{
			if (!sample_bundles[i]->hits().empty())
			{
				non_empty_bundle = true;
				break;
			}
		}
		if (non_empty_bundle)
		{
			RefID bundle_chr_id = sample_bundles.front()->ref_id();
			assert (bundle_chr_id != 0);
			const char* chr_name = rt.get_name(bundle_chr_id);
			assert (chr_name);
			fprintf(stderr, "Quantitating samples in locus [ %s:%d-%d ] \n", 
					chr_name,
					sample_bundles.front()->left(),
					sample_bundles.front()->right());
			test_differential(rt, sample_bundles, map_masses, tests, tracking);
			
			for (size_t i = 0; i < sample_bundles.size(); ++i)
			{
				delete sample_bundles[i];
			}
		}
	}
	
	// wait for the workers to finish up before reporting everthing.
#if ENABLE_THREADS	
	while(1)
	{
		thread_pool_lock.lock();
		if (curr_threads == 0)
		{
			thread_pool_lock.unlock();
			break;
		}
		
		thread_pool_lock.unlock();
		//fprintf(stderr, "waiting to for all workers to finish\n");
		boost::this_thread::sleep(boost::posix_time::milliseconds(5));
	}
#endif
	
	//double FDR = 0.05;
	int total_iso_de_tests = 0;
	
	vector<SampleDifference*> isoform_exp_diffs;
	for (size_t i = 0; i < tests.isoform_de_tests.size(); ++i)
	{
		total_iso_de_tests += tests.isoform_de_tests[i].size();
		extract_sample_diffs(tests.isoform_de_tests[i], isoform_exp_diffs);
	}
	int iso_exp_tests = fdr_significance(FDR, isoform_exp_diffs);
	fprintf(stderr, "Performed %d isoform-level transcription difference tests\n", iso_exp_tests);
	for (size_t i = 0; i < tests.isoform_de_tests.size(); ++i)
	{
		FILE* fout = outfiles.isoform_de_outfiles[i];
		print_tests(fout, "log(fold_change)", tests.isoform_de_tests[i]);
	}
	
	int total_group_de_tests = 0;
	vector<SampleDifference*> tss_group_exp_diffs;
	for (size_t i = 0; i < tests.tss_group_de_tests.size(); ++i)
	{
		extract_sample_diffs(tests.tss_group_de_tests[i], tss_group_exp_diffs);
		total_group_de_tests += tests.tss_group_de_tests[i].size();
	}
	
	int tss_group_exp_tests = fdr_significance(FDR, tss_group_exp_diffs);
	fprintf(stderr, "Performed %d tss-level transcription difference tests\n", tss_group_exp_tests);
	for (size_t i = 0; i < tests.tss_group_de_tests.size(); ++i)
	{
		FILE* fout = outfiles.group_de_outfiles[i];
		print_tests(fout, "log(fold_change)", tests.tss_group_de_tests[i]);
	}
	
	int total_gene_de_tests = 0;
	vector<SampleDifference*> gene_exp_diffs;
	for (size_t i = 0; i < tests.gene_de_tests.size(); ++i)
	{
		total_gene_de_tests += tests.gene_de_tests[i].size();
		extract_sample_diffs(tests.gene_de_tests[i], gene_exp_diffs);
	}
	
	int gene_exp_tests = fdr_significance(FDR, gene_exp_diffs);
	fprintf(stderr, "Performed %d gene-level transcription difference tests\n", gene_exp_tests);
	for (size_t i = 0; i < tests.gene_de_tests.size(); ++i)
	{
		FILE* fout = outfiles.gene_de_outfiles[i];
		print_tests(fout, "log(fold_change)", tests.gene_de_tests[i]);
	}	

	int total_cds_de_tests = 0;
	vector<SampleDifference*> cds_exp_diffs;
	for (size_t i = 0; i < tests.cds_de_tests.size(); ++i)
	{
		total_cds_de_tests += tests.cds_de_tests[i].size();
		extract_sample_diffs(tests.cds_de_tests[i], cds_exp_diffs);
	}
	int cds_exp_tests = fdr_significance(FDR, cds_exp_diffs);
	fprintf(stderr, "Performed %d CDS-level transcription difference tests\n", cds_exp_tests);
	for (size_t i = 0; i < tests.cds_de_tests.size(); ++i)
	{
		FILE* fout = outfiles.cds_de_outfiles[i];
		print_tests(fout, "log(fold_change)", tests.cds_de_tests[i]);
	}
	
	int total_diff_splice_tests = 0;
	vector<SampleDifference*> splicing_diffs;
	for (size_t i = 0; i < tests.diff_splicing_tests.size(); ++i)
	{
		total_diff_splice_tests += tests.diff_splicing_tests[i].size();
		extract_sample_diffs(tests.diff_splicing_tests[i], splicing_diffs);
	}
	
	int splicing_tests = fdr_significance(FDR, splicing_diffs);
	fprintf(stderr, "Performed %d splicing tests\n", splicing_tests);
	for (size_t i = 0; i < tests.diff_splicing_tests.size(); ++i)
	{
		FILE* fout = outfiles.diff_splicing_outfiles[i];
		const SampleDiffs& diffs = tests.diff_splicing_tests[i];
		print_tests(fout, "sqrt(JS)", diffs);
	}
	
	int total_diff_promoter_tests = 0;
	vector<SampleDifference*> promoter_diffs;
	for (size_t i = 0; i < tests.diff_splicing_tests.size(); ++i)
	{
		total_diff_promoter_tests += tests.diff_promoter_tests[i].size();
		extract_sample_diffs(tests.diff_promoter_tests[i], promoter_diffs);

	}
	int promoter_tests = fdr_significance(FDR, promoter_diffs);
	fprintf(stderr, "Performed %d promoter preference tests\n", promoter_tests);
	for (size_t i = 0; i < tests.diff_promoter_tests.size(); ++i)
	{
		FILE* fout = outfiles.diff_promoter_outfiles[i];
		print_tests(fout, "sqrt(JS)", tests.diff_promoter_tests[i]);
	}

	int total_diff_cds_tests = 0;
	vector<SampleDifference*> cds_use_diffs;
	for (size_t i = 0; i < tests.diff_cds_tests.size(); ++i)
	{
		extract_sample_diffs(tests.diff_cds_tests[i], cds_use_diffs);
		total_diff_cds_tests += tests.diff_cds_tests[i].size();
	}
	int cds_use_tests = fdr_significance(FDR, cds_use_diffs);
	fprintf(stderr, "Performing %d relative CDS output tests\n", cds_use_tests);
	for (size_t i = 0; i < tests.diff_cds_tests.size(); ++i)
	{
		FILE* fout = outfiles.diff_cds_outfiles[i];
		print_tests(fout, "sqrt(JS)", tests.diff_cds_tests[i]);
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
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string ref_gtf_filename = argv[optind++];
	
	vector<FILE*> sam_hit_files;
	vector<string> sam_hit_file_names;
    while(optind < argc)
    {
        string sam_hits_file_name = argv[optind++];
		// Open the approppriate files
		FILE* sam_hits_file = fopen(sam_hits_file_name.c_str(), "r");
		if (sam_hits_file == NULL)
		{
			fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
					sam_hits_file_name.c_str());
			exit(1);
		}
		
		sam_hit_files.push_back(sam_hits_file);
		sam_hit_file_names.push_back(sam_hits_file_name);
    }
    	
	while (sam_hit_files.size() < 2)
    {
        fprintf(stderr, "Error: cuffdiff requires at least 2 SAM files\n");
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
			fprintf(stderr, "Error: cannot open GTF file %s for reading\n",
					ref_gtf_filename.c_str());
			exit(1);
		}
	}
	
	
	// Note: we don't want the assembly filters interfering with calculations 
	// here
	min_isoform_fraction = 0.0;
	pre_mrna_fraction = 0.0;
	
	Outfiles outfiles;
	
	for (size_t i = 1; i < sam_hit_files.size(); ++i)
	{
		char out_file_prefix[64];
		sprintf(out_file_prefix, "%s/%lu_%lu", output_dir.c_str(), i - 1, i);
		char iso_out_file_name[256];
		sprintf(iso_out_file_name, "%s_isoform_exp.diff", out_file_prefix);
		FILE* iso_out = fopen(iso_out_file_name, "w");
		if (!iso_out)
		{
			fprintf(stderr, "Error: cannot open differential isoform transcription file %s for writing\n",
					iso_out_file_name);
			exit(1);
		}
		
		char group_out_file_name[256];
		sprintf(group_out_file_name, "%s_tss_group_exp.diff", out_file_prefix);
		FILE* group_out = fopen(group_out_file_name, "w");
		if (!group_out)
		{
			fprintf(stderr, "Error: cannot open differential TSS group transcription file %s for writing\n",
					group_out_file_name);
			exit(1);
		}
		
		char gene_out_file_name[256];
		sprintf(gene_out_file_name, "%s_gene_exp.diff", out_file_prefix);
		FILE* gene_out = fopen(gene_out_file_name, "w");
		if (!group_out)
		{
			fprintf(stderr, "Error: cannot open gene expression file %s for writing\n",
					gene_out_file_name);
			exit(1);
		}
		
		char cds_out_file_name[256];
		sprintf(cds_out_file_name, "%s_cds_exp.diff", out_file_prefix);
		FILE* cds_out = fopen(cds_out_file_name, "w");
		if (!cds_out)
		{
			fprintf(stderr, "Error: cannot open cds expression file %s for writing\n",
					cds_out_file_name);
			exit(1);
		}
		
		char diff_splicing_out_file_name[256];
		sprintf(diff_splicing_out_file_name, "%s_splicing.diff", out_file_prefix);
		FILE* diff_splicing_out = fopen(diff_splicing_out_file_name, "w");
		if (!diff_splicing_out)
		{
			fprintf(stderr, "Error: cannot open differential splicing file %s for writing\n",
					diff_splicing_out_file_name);
			exit(1);
		}
		
		char diff_promoter_out_file_name[256];
		sprintf(diff_promoter_out_file_name, "%s_promoters.diff", out_file_prefix);
		FILE* diff_promoter_out = fopen(diff_promoter_out_file_name, "w");
		if (!diff_promoter_out)
		{
			fprintf(stderr, "Error: cannot open differential transcription start file %s for writing\n",
					diff_promoter_out_file_name);
			exit(1);
		}
		
		char diff_cds_out_file_name[256];
		sprintf(diff_cds_out_file_name, "%s_cds.diff", out_file_prefix);
		FILE* diff_cds_out = fopen(diff_cds_out_file_name, "w");
		if (!diff_cds_out)
		{
			fprintf(stderr, "Error: cannot open differential relative CDS file %s for writing\n",
					diff_cds_out_file_name);
			exit(1);
		}
		
		outfiles.isoform_de_outfiles.push_back(iso_out);
		outfiles.group_de_outfiles.push_back(group_out);
		outfiles.gene_de_outfiles.push_back(gene_out);
		outfiles.cds_de_outfiles.push_back(cds_out);
		outfiles.diff_splicing_outfiles.push_back(diff_splicing_out);
		outfiles.diff_promoter_outfiles.push_back(diff_promoter_out);
		outfiles.diff_cds_outfiles.push_back(diff_cds_out);
	}
	
	char isoform_fpkm_tracking_name[256];
	sprintf(isoform_fpkm_tracking_name, "isoforms.fpkm_tracking");
	FILE* isoform_fpkm_out = fopen(isoform_fpkm_tracking_name, "w");
	if (!isoform_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level FPKM tracking file %s for writing\n",
				isoform_fpkm_tracking_name);
		exit(1);
	}
	outfiles.isoform_fpkm_tracking_out = isoform_fpkm_out;

	char tss_group_fpkm_tracking_name[256];
	sprintf(tss_group_fpkm_tracking_name, "tss_groups.fpkm_tracking");
	FILE* tss_group_fpkm_out = fopen(tss_group_fpkm_tracking_name, "w");
	if (!tss_group_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level FPKM tracking file %s for writing\n",
				tss_group_fpkm_tracking_name);
		exit(1);
	}
	outfiles.tss_group_fpkm_tracking_out = tss_group_fpkm_out;

	char cds_fpkm_tracking_name[256];
	sprintf(cds_fpkm_tracking_name, "cds.fpkm_tracking");
	FILE* cds_fpkm_out = fopen(cds_fpkm_tracking_name, "w");
	if (!cds_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open CDS level FPKM tracking file %s for writing\n",
				cds_fpkm_tracking_name);
		exit(1);
	}
	outfiles.cds_fpkm_tracking_out = cds_fpkm_out;
	
	char gene_fpkm_tracking_name[256];
	sprintf(gene_fpkm_tracking_name, "genes.fpkm_tracking");
	FILE* gene_fpkm_out = fopen(gene_fpkm_tracking_name, "w");
	if (!gene_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open gene-level FPKM tracking file %s for writing\n",
				gene_fpkm_tracking_name);
		exit(1);
	}
	outfiles.gene_fpkm_tracking_out = gene_fpkm_out;
	
    driver(ref_gtf, sam_hit_files, outfiles);
	
	return 0;
}

