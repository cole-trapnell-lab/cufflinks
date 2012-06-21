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
#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "differential.h"

// Need at least this many reads in a locus to do any testing on it

vector<string> sample_labels;

double FDR = 0.05; 
bool samples_are_time_series = false;
using namespace std;
using namespace boost;

bool use_geometric_norm = false;
bool use_raw_mapped_norm = false;
bool use_isoform_count_dispersion = true;

// We leave out the short codes for options that don't take an argument
#if ENABLE_THREADS
const char *short_options = "m:p:s:c:I:j:L:M:o:b:TNqvuF:";
#else
const char *short_options = "m:s:c:I:j:L:M:o:b:TNqvuF:";
#endif



static struct option long_options[] = {
{"frag-len-mean",			required_argument,       0,          'm'},
{"frag-len-std-dev",        required_argument,       0,          's'},
{"transcript-score-thresh", required_argument,       0,          't'},
{"pre-mrna-fraction",		required_argument,		 0,			 'j'},
{"max-intron-length",		required_argument,		 0,			 'I'},
{"labels",					required_argument,		 0,			 'L'},
{"min-alignment-count",     required_argument,		 0,			 'c'},
{"FDR",					    required_argument,		 0,			 OPT_FDR},
{"seed",                    required_argument,		 0,			 OPT_RANDOM_SEED},
{"mask-file",                required_argument,		 0,			 'M'},
{"output-dir",			    required_argument,		 0,			 'o'},
{"verbose",			    	no_argument,			 0,			 'v'},
{"quiet",			    	no_argument,			 0,			 'q'},
{"frag-bias-correct",       required_argument,		 0,			 'b'},
{"multi-read-correct",      no_argument,			 0,			 'u'},
{"time-series",             no_argument,             0,			 'T'},
{"upper-quartile-norm",     no_argument,	 		 0,	         'N'},
{"geometric-norm",          no_argument,	 		 0,	         OPT_GEOMETRIC_NORM},
{"raw-mapped-norm",         no_argument,	 		 0,	         OPT_RAW_MAPPED_NORM},
{"min-isoform-fraction",    required_argument,       0,          'F'},
#if ENABLE_THREADS
{"num-threads",				required_argument,       0,          'p'},
#endif
{"library-type",		    required_argument,		 0,			 OPT_LIBRARY_TYPE},
{"seed",                    required_argument,		 0,			 OPT_RANDOM_SEED},
{"no-collapse-cond-prob",   no_argument,             0,			 OPT_COLLAPSE_COND_PROB},
{"num-importance-samples",  required_argument,		 0,			 OPT_NUM_IMP_SAMPLES},
{"max-mle-iterations",		required_argument,		 0,			 OPT_MLE_MAX_ITER},
{"min-mle-accuracy",		required_argument,		 0,			 OPT_MLE_MIN_ACC},
{"poisson-dispersion",		no_argument,             0,		     OPT_POISSON_DISPERSION},
{"bias-mode",               required_argument,		 0,			 OPT_BIAS_MODE},
{"no-update-check",         no_argument,             0,          OPT_NO_UPDATE_CHECK},
{"emit-count-tables",       no_argument,             0,          OPT_EMIT_COUNT_TABLES},
{"compatible-hits-norm",    no_argument,	 		 0,	         OPT_USE_COMPAT_MASS},
{"total-hits-norm",         no_argument,	 		 0,	         OPT_USE_TOTAL_MASS},
//{"analytic-diff",           no_argument,	 		 0,	         OPT_ANALYTIC_DIFF},
{"no-diff",                 no_argument,	 		 0,	         OPT_NO_DIFF},
{"num-frag-count-draws",	required_argument,		 0,			 OPT_NUM_FRAG_COUNT_DRAWS},
{"num-frag-assign-draws",	required_argument,		 0,			 OPT_NUM_FRAG_ASSIGN_DRAWS},
    
// Some options for testing different stats policies
{"fisher-covariance",       no_argument,	 		 0,	         OPT_USE_FISHER_COVARIANCE},
{"empirical-covariance",    no_argument,	 		 0,	         OPT_USE_EMPIRICAL_COVARIANCE},
{"split-mass",              no_argument,	 		 0,	         OPT_SPLIT_MASS},
{"split-variance",          no_argument,	 		 0,	         OPT_SPLIT_VARIANCE},
{"num-bootstrap-samples",   required_argument,	 	 0,          OPT_NUM_BOOTSTRAP_SAMPLES},
{"bootstrap-fraction",      required_argument,	 	 0,          OPT_BOOTSTRAP_FRACTION},
{"max-bundle-frags",        required_argument,       0,          OPT_MAX_FRAGS_PER_BUNDLE}, 
{"read-skip-fraction",      required_argument,	     0,          OPT_READ_SKIP_FRACTION},
{"no-read-pairs",           no_argument,	 		 0,          OPT_NO_READ_PAIRS},
{"trim-read-length",        required_argument,	     0,          OPT_TRIM_READ_LENGTH},
{"cov-delta",               required_argument,	     0,          OPT_MAX_DELTA_GAP},
{"locus-count-dispersion",  no_argument,             0,          OPT_LOCUS_COUNT_DISPERSION},
{"max-frag-multihits",      required_argument,       0,          OPT_FRAG_MAX_MULTIHITS},
{"min-outlier-p",           required_argument,       0,          OPT_MIN_OUTLIER_P},
{"min-reps-for-js-test",      required_argument,     0,          OPT_MIN_REPS_FOR_JS_TEST},
{"no-effective-length-correction",  no_argument,     0,          OPT_NO_EFFECTIVE_LENGTH_CORRECTION},
{"no-length-correction",  no_argument,     0,          OPT_NO_LENGTH_CORRECTION},
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "cuffdiff v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]\n");
	fprintf(stderr, "   Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "  -o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");
    fprintf(stderr, "  --seed                       value of random number generator seed                 [ default:      0 ]\n");
    fprintf(stderr, "  -T/--time-series             treat samples as a time-series                        [ default:  FALSE ]\n");
	fprintf(stderr, "  -c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   10 ]\n");
	fprintf(stderr, "  --FDR                        False discovery rate used in testing                  [ default:   0.05 ]\n");
	fprintf(stderr, "  -M/--mask-file               ignore all alignment within transcripts in this file  [ default:   NULL ]\n");
    fprintf(stderr, "  -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]\n");
    fprintf(stderr, "  -u/--multi-read-correct      use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]\n");
    fprintf(stderr, "  -N/--upper-quartile-norm     use upper-quartile normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, "  --geometric-norm             use geometric mean normalization                      [ default:  TRUE ]\n");
    fprintf(stderr, "  --raw-mapped-norm            use raw mapped count normalized (classic FPKM)        [ default:  FALSE ]\n");
    fprintf(stderr, "  -L/--labels                  comma-separated list of condition labels\n");
#if ENABLE_THREADS
	fprintf(stderr, "  -p/--num-threads             number of threads used during quantification          [ default:      1 ]\n");
#endif
    fprintf(stderr, "  --no-diff                    Don't generate differential analysis files            [ default:  FALSE ]\n");
	fprintf(stderr, "\nAdvanced Options:\n");
    fprintf(stderr, "  --library-type               Library prep used for input reads                     [ default:  below ]\n");
    fprintf(stderr, "  -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]\n");
    fprintf(stderr, "  -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]\n");
    fprintf(stderr, "  --num-importance-samples     number of importance samples for MAP restimation      [ DEPRECATED      ]\n");
	fprintf(stderr, "  --num-bootstrap-samples      Number of bootstrap replications                      [ DEPRECATED      ]\n");
    fprintf(stderr, "  --bootstrap-fraction         Fraction of fragments in each bootstrap sample        [ DEPRECATED      ]\n");
    fprintf(stderr, "  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]\n");
    fprintf(stderr, "  --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:   TRUE ]\n");
    fprintf(stderr, "  --total-hits-norm            count all hits for normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, "  --poisson-dispersion         Don't fit fragment counts for overdispersion          [ default:  FALSE ]\n");
    fprintf(stderr, "  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
	fprintf(stderr, "  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
    fprintf(stderr, "  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]\n");    
    fprintf(stderr, "  --emit-count-tables          print count tables used to fit overdispersion         [ default:  FALSE ]\n");
    fprintf(stderr, "  --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]\n");
    fprintf(stderr, "  --num-frag-count-draws       Number of fragment generation samples                 [ default:   1000 ]\n");
    fprintf(stderr, "  --num-frag-assign-draws      Number of fragment assignment samples per generation  [ default:      1 ]\n");
    fprintf(stderr, "  --max-frag-multihits         Maximum number of alignments allowed per fragment     [ default: unlim  ]\n");
    fprintf(stderr, "  --min-outlier-p              Min replicate p value to admit for testing            [ default:   0.01 ]\n");
    fprintf(stderr, "  --min-reps-for-js-test       Replicates needed for relative isoform shift testing  [ default:      3 ]\n");
    fprintf(stderr, "  --no-effective-length-correction   No effective length correction                  [ default:  FALSE ]\n");
    fprintf(stderr, "  --no-length-correction       No effective length correction                        [ default:  FALSE ]\n");
    fprintf(stderr, "\nDebugging use only:\n");
    fprintf(stderr, "  --read-skip-fraction         Skip a random subset of reads this size               [ default:    0.0 ]\n");
    fprintf(stderr, "  --no-read-pairs              Break all read pairs                                  [ default:  FALSE ]\n");
    fprintf(stderr, "  --trim-read-length           Trim reads to be this long (keep 5' end)              [ default:   none ]\n");
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
            case 'F':
				min_isoform_fraction = parseFloat(0, 1.0, "-F/--min-isoform-fraction must be between 0 and 1.0", print_usage);
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
                else if (!strcmp(optarg, "pos_site"))
					bias_mode = POS_SITE;
				else
				{
					fprintf(stderr, "Unknown bias mode.\n");
					exit(1);
				}
				break;
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
			case 'b':
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
            case 'u':
            {
                corr_multi = true;
                break;
            }
            case OPT_LIBRARY_TYPE:
			{
				library_type = optarg;
				break;
			}
            case OPT_POISSON_DISPERSION:
			{
				poisson_dispersion = true;
				break;
			}
            case OPT_NO_UPDATE_CHECK:
            {
                no_update_check = true;
                break;
            }
            case OPT_RANDOM_SEED:
            {
                random_seed = parseInt(0, "--seed must be at least 0", print_usage);
                break;
            }
            case OPT_EMIT_COUNT_TABLES:
            {
                emit_count_tables = true;
                break;
            }    
            case OPT_COLLAPSE_COND_PROB:
            {
                cond_prob_collapse = false;
                break;
            }
            case OPT_USE_COMPAT_MASS:
            {
                use_compat_mass = true;
                break;
            }
            case OPT_USE_TOTAL_MASS:
            {
                use_total_mass = true;
                break;
            }
            case OPT_USE_FISHER_COVARIANCE:
            {
                use_fisher_covariance = true;
                break;
            }
            case OPT_USE_EMPIRICAL_COVARIANCE:
            {
                use_fisher_covariance = false;
                break;
            }
            case OPT_SPLIT_MASS:
            {
                split_variance = false;
                break;
            }
            case OPT_SPLIT_VARIANCE:
            {
                split_variance = true;
                break;
            }
            case OPT_NUM_BOOTSTRAP_SAMPLES:
            {
                //num_bootstrap_samples = parseInt(1, "--num-bootstrap-samples must be at least 1", print_usage);
                break;
            }
            case OPT_BOOTSTRAP_FRACTION:
            {
                bootstrap_fraction = parseFloat(0, 1.0, "--bootstrap-fraction must be between 0 and 1.0", print_usage);
                break;
            }
            case OPT_MAX_FRAGS_PER_BUNDLE:
            {
                max_frags_per_bundle = parseInt(0, "--max-bundle-frags must be at least 0", print_usage);
                break;
            }
            case OPT_READ_SKIP_FRACTION:
            {
                read_skip_fraction = parseFloat(0, 1.0, "--read-skip-fraction must be between 0 and 1.0", print_usage);
                break;
            }
            case OPT_NO_READ_PAIRS:
            {
                no_read_pairs = true;
                break;
            }
            case OPT_TRIM_READ_LENGTH:
            {
                trim_read_length = parseInt(0, "--trim-read-length must be at least 1", print_usage);
                break;
            }
            case OPT_MAX_DELTA_GAP:
            {
                bootstrap_delta_gap = parseFloat(0, 10000000.0, "--read-skip-fraction must be between 0 and 10000000.0", print_usage);
                break;
            }
            case OPT_MLE_MIN_ACC:
            {
                bootstrap_delta_gap = parseFloat(0, 10000000.0, "--read-skip-fraction must be between 0 and 10000000.0", print_usage);
                break;
            }
//            case OPT_ANALYTIC_DIFF:
//            {
//                analytic_diff = true;
//                break;
//            }
            case OPT_NO_DIFF:
            {
                no_differential = true;
                break;
            }
            case OPT_GEOMETRIC_NORM:
            {
                use_geometric_norm = true;
                break;
            } 
            case OPT_RAW_MAPPED_NORM:
            {
                use_raw_mapped_norm = true;
                break;
            } 
            case OPT_NUM_FRAG_COUNT_DRAWS:
            {
                num_frag_count_draws = parseInt(1, "--num-frag-count-draws must be at least 1", print_usage);
                break;
            }
            case OPT_NUM_FRAG_ASSIGN_DRAWS:
            {
                num_frag_assignments = parseInt(1, "--num-frag-assign-draws must be at least 1", print_usage);
                break;
            }
            case OPT_LOCUS_COUNT_DISPERSION:
            {
                use_isoform_count_dispersion = false;
                break;
            }
            case OPT_FRAG_MAX_MULTIHITS:
            {
                max_frag_multihits = parseInt(1, "--max-frag-multihits must be at least 1", print_usage);
                break;
            }
            case OPT_MIN_OUTLIER_P:
            {
                min_outlier_p = parseFloat(0, 1.0, "--min-outlier-p must be between 0 and 1.0", print_usage);
                break;
            }
            case OPT_MIN_REPS_FOR_JS_TEST:
            {
                min_reps_for_js_test = parseInt(1, "--min-reps-for-js-test must be at least 1", print_usage);
                break;
            }
            case OPT_NO_EFFECTIVE_LENGTH_CORRECTION:
            {
                no_effective_length_correction = true;
                break;
            }
            case OPT_NO_LENGTH_CORRECTION:
            {
                no_length_correction = true;
                break;
            }
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	if (library_type != "")
    {
        map<string, ReadGroupProperties>::iterator lib_itr = 
		library_type_table.find(library_type);
        if (lib_itr == library_type_table.end())
        {
            fprintf(stderr, "Error: Library type %s not supported\n", library_type.c_str());
            exit(1);
        }
        else 
        {
            if (library_type == "transfrags")
            {
                allow_junk_filtering = false;
            }
            global_read_properties = &lib_itr->second;
        }
    }
    
    if (use_total_mass && use_compat_mass)
    {
        fprintf (stderr, "Error: please supply only one of --compatibile-hits-norm and --total-hits-norm\n");
        exit(1);
    }
    
    if (use_raw_mapped_norm + use_geometric_norm + use_quartile_norm > 1)
    {
        fprintf(stderr, "Error: Choose one of upper-quartile-norm, geometric-norm, or raw-mapped-norm\n");
        exit(1);
    }
    if (use_raw_mapped_norm + use_geometric_norm + use_quartile_norm == 0)
    {
        use_geometric_norm = true;
    }
    
    tokenize(sample_label_list, ",", sample_labels);
    
	allow_junk_filtering = false;
	
	return 0;
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
        
        string all_gene_ids = cat_strings(test.meta_data->gene_ids);
		if (all_gene_ids == "")
			all_gene_ids = "-";
        
		string all_gene_names = cat_strings(test.meta_data->gene_names);
		if (all_gene_names == "")
			all_gene_names = "-";
		
		string all_protein_ids = cat_strings(test.meta_data->protein_ids);	
		if (all_protein_ids == "")
			all_protein_ids = "-";
		
		fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s", 
                itr->first.c_str(), 
                all_gene_ids.c_str(),
                all_gene_names.c_str(), 
                test.meta_data->locus_desc.c_str(),
                sample_1_label,
                sample_2_label);
		
        double t = test.test_stat;
        double r1 = test.value_1;
        double r2 = test.value_2;
        double d = test.differential;
        double p = test.p_value;
        double q = test.corrected_p;
        const char* sig;
        if (test.significant && test.test_status == OK)
            sig = "yes";
        else
            sig = "no";
        
        const char* status;
        if (test.test_status == OK)
            status = "OK";
        else if (test.test_status == LOWDATA)
            status = "LOWDATA";
        else if (test.test_status == HIDATA)
            status = "HIDATA";
        else if (test.test_status == NOTEST)
            status = "NOTEST";
        else
            status = "FAIL";
        
        fprintf(fout, "\t%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%s", status, r1, r2, d, t, p, q, sig);
        fprintf(fout, "\n");
	}
}

void print_FPKM_tracking(FILE* fout, 
						 const FPKMTrackingTable& tracking)
{
	fprintf(fout,"tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage");
	FPKMTrackingTable::const_iterator first_itr = tracking.begin();
	if (first_itr != tracking.end())
	{
		const FPKMTracking& track = first_itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
			fprintf(fout, "\t%s_FPKM\t%s_conf_lo\t%s_conf_hi\t%s_status", sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str());
		}
	}
	fprintf(fout, "\n");
	for (FPKMTrackingTable::const_iterator itr = tracking.begin(); itr != tracking.end(); ++itr)
	{
		const string& description = itr->first;
		const FPKMTracking& track = itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		
        AbundanceStatus status = NUMERIC_OK;
        foreach (const FPKMContext& c, fpkms)
        {
            if (c.status == NUMERIC_FAIL)
                status = NUMERIC_FAIL;
        }
        
        string all_gene_ids = cat_strings(track.gene_ids);
		if (all_gene_ids == "")
			all_gene_ids = "-";
        
		string all_gene_names = cat_strings(track.gene_names);
		if (all_gene_names == "")
			all_gene_names = "-";
		
		string all_tss_ids = cat_strings(track.tss_ids);
		if (all_tss_ids == "")
			all_tss_ids = "-";
		
        char length_buff[33] = "-";
        if (track.length)
            sprintf(length_buff, "%d", track.length);
        
        fprintf(fout, "%s\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
                description.c_str(),
                track.classcode ? track.classcode : '-',
                track.ref_match.c_str(),
                all_gene_ids.c_str(),
                all_gene_names.c_str(), 
                all_tss_ids.c_str(),
                track.locus_tag.c_str(),
                length_buff,
                "-");
       		
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
			double fpkm = fpkms[i].FPKM;
			double std_dev = sqrt(fpkms[i].FPKM_variance);
			double fpkm_conf_hi = fpkm + 2.0 * std_dev;
			double fpkm_conf_lo = max(0.0, fpkm - 2.0 * std_dev);
            const char* status_str = "OK";
            
            if (fpkms[i].status == NUMERIC_OK)
            {
                status_str = "OK";
            }
            else if (fpkms[i].status == NUMERIC_FAIL)
            {
                status_str = "FAIL";
            }
            else if (fpkms[i].status == NUMERIC_LOW_DATA)
            {
                status_str = "LOWDATA";
            }
            else if (fpkms[i].status == NUMERIC_HI_DATA)
            {
                status_str = "HIDATA";
            }
            else
            {
                assert(false);
            }
            
			fprintf(fout, "\t%lg\t%lg\t%lg\t%s", fpkm, fpkm_conf_lo, fpkm_conf_hi, status_str);
		}
		
		fprintf(fout, "\n");
	}
}

void print_count_tracking(FILE* fout, 
						  const FPKMTrackingTable& tracking)
{
	fprintf(fout,"tracking_id");
	FPKMTrackingTable::const_iterator first_itr = tracking.begin();
	if (first_itr != tracking.end())
	{
		const FPKMTracking& track = first_itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
			fprintf(fout, "\t%s_count\t%s_count_variance\t%s_count_uncertainty_var\t%s_count_dispersion_var\t%s_status", sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str());
		}
	}
	fprintf(fout, "\n");
	for (FPKMTrackingTable::const_iterator itr = tracking.begin(); itr != tracking.end(); ++itr)
	{
		const string& description = itr->first;
		const FPKMTracking& track = itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
		
        AbundanceStatus status = NUMERIC_OK;
        foreach (const FPKMContext& c, fpkms)
        {
            if (c.status == NUMERIC_FAIL)
                status = NUMERIC_FAIL;
        }
        
        fprintf(fout, "%s", 
                description.c_str());
        
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
            const char* status_str = "OK";
            
            if (fpkms[i].status == NUMERIC_OK)
            {
                status_str = "OK";
            }
            else if (fpkms[i].status == NUMERIC_FAIL)
            {
                status_str = "FAIL";
            }
            else if (fpkms[i].status == NUMERIC_LOW_DATA)
            {
                status_str = "LOWDATA";
            }
            else if (fpkms[i].status == NUMERIC_HI_DATA)
            {
                status_str = "HIDATA";
            }
            else
            {
                assert(false);
            }
            
            double external_counts = fpkms[i].count_mean;
            double external_count_var = fpkms[i].count_var;
            double uncertainty_var = fpkms[i].count_uncertainty_var;
            double dispersion_var = fpkms[i].count_dispersion_var;
			fprintf(fout, "\t%lg\t%lg\t%lg\t%lg\t%s", external_counts, external_count_var, uncertainty_var, dispersion_var, status_str);
		}
		
		fprintf(fout, "\n");
	}
}

void print_read_group_tracking(FILE* fout, 
                               const FPKMTrackingTable& tracking)
{
	fprintf(fout,"tracking_id\tcondition\treplicate\traw_frags\tinternal_scaled_frags\texternal_scaled_frags\tFPKM\teffective_length\tstatus");
	
	fprintf(fout, "\n");
	for (FPKMTrackingTable::const_iterator itr = tracking.begin(); itr != tracking.end(); ++itr)
	{
		const string& description = itr->first;
		const FPKMTracking& track = itr->second;
		const vector<FPKMContext>& fpkms = track.fpkm_series;
                
		for (size_t i = 0; i < fpkms.size(); ++i)
		{
            for (CountPerReplicateTable::const_iterator itr = fpkms[i].count_per_rep.begin(); 
                 itr != fpkms[i].count_per_rep.end(); 
                 ++itr)
            { 
                FPKMPerReplicateTable::const_iterator f_itr = fpkms[i].fpkm_per_rep.find(itr->first);
                StatusPerReplicateTable::const_iterator s_itr = fpkms[i].status_per_rep.find(itr->first);
                
                
                if (f_itr == fpkms[i].fpkm_per_rep.end())
                {
                    fprintf(stderr, "Error: missing per-replicate FPKM data\n");
                }
                
                double FPKM = f_itr->second;
                double internal_count = itr->second;
                double external_count = internal_count / itr->first->external_scale_factor();
                double raw_count = internal_count * itr->first->internal_scale_factor();
                const  string& condition_name = itr->first->condition_name();
                AbundanceStatus status = s_itr->second;
                
                int rep_num = itr->first->replicate_num();
                
                const char* status_str = "OK";
                
                if (status == NUMERIC_OK)
                {
                    status_str = "OK";
                }
                else if (status == NUMERIC_FAIL)
                {
                    status_str = "FAIL";
                }
                else if (status == NUMERIC_LOW_DATA)
                {
                    status_str = "LOWDATA";
                }
                else if (status == NUMERIC_HI_DATA)
                {
                    status_str = "HIDATA";
                }
                else
                {
                    assert(false);
                }
                
                fprintf(fout, "%s\t%s\t%d\t%lg\t%lg\t%lg\t%lg\t%s\t%s\n",
                        description.c_str(),
                        condition_name.c_str(),
                        rep_num,
                        raw_count,
                        internal_count,
                        external_count,
                        FPKM,
                        "-",
                        status_str);
            }
		}
	}
}

void print_read_group_info(FILE* fout, 
                           const vector<shared_ptr<ReadGroupProperties> >& all_read_groups)
{
    fprintf(fout, "file\tcondition\treplicate_num\ttotal_mass\tnorm_mass\tinternal_scale\texternal_scale\n");
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        shared_ptr<ReadGroupProperties const> rg_props = all_read_groups[i];
        fprintf(fout, "%s\t%s\t%d\t%Lg\t%Lg\t%lg\t%lg\n",
                rg_props->file_path().c_str(),
                rg_props->condition_name().c_str(),
                rg_props->replicate_num(),
                rg_props->total_map_mass(),
                rg_props->normalized_map_mass(),
                rg_props->internal_scale_factor(),
                rg_props->external_scale_factor());
                
    }
}

void print_run_info(FILE* fout)
{
    fprintf(fout, "param\tvalue\n");
    fprintf(fout, "cmd_line\t%s\n", cmd_str.c_str());
    fprintf(fout, "version\t%s\n", PACKAGE_VERSION);
    fprintf(fout, "SVN_revision\t%s\n",SVN_REVISION); 
    fprintf(fout, "boost_version\t%d\n", BOOST_VERSION);
}

bool p_value_lt(const SampleDifference* lhs, const SampleDifference* rhs)
{
	return lhs->p_value < rhs->p_value;
}

//// Benjamani-Hochberg procedure
//int fdr_significance(double fdr, 
//					  vector<SampleDifference*>& tests)
//{
//	sort(tests.begin(), tests.end(), p_value_lt);
//	vector<SampleDifference*> passing;
//
//	for (int k = 0; k < (int)tests.size(); ++k)
//	{
//		if (tests[k]->test_status == OK)
//		{
//			passing.push_back(tests[k]);
//		}
//		else
//		{
//			tests[k]->significant = false;
//		}
//	}
//	int significant = 0;
//	for (int k = 0; k < (int)passing.size(); ++k)
//	{
//		double r = (double)passing.size() / ((double) k + 1);
//		double corrected_p = passing[k]->p_value * r;
//		passing[k]->corrected_p = corrected_p;
//		passing[k]->significant = (corrected_p <= fdr);
//        significant += passing[k]->significant;
//	}
//    
//	return passing.size();
//}

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
	int significant = 0;
	float pmin=1;
	int n = (int) passing.size();
    //use the same procedure as p.adjust(...,"BH") in R
	for (int k = n-1; k >= 0; k--)
	{
		double corrected_p = (double) passing[k]->p_value * ((double) n/(double) (k+1));
        //make sure that no entry with lower p-value will get higher q-value than any entry with higher p-value
		if (corrected_p < pmin) 
		{
			pmin = corrected_p;
		}
		else
		{
			corrected_p = pmin;
		}
        // make sure that the q-value is always <= 1 
		passing[k]->corrected_p = (corrected_p < 1 ? corrected_p : 1); 
		passing[k]->significant = (corrected_p <= fdr);
        significant += passing[k]->significant;
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
boost::mutex inspect_lock;
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
}


shared_ptr<TestLauncher> test_launcher;

bool quantitate_next_locus(const RefSequenceTable& rt,
                           vector<shared_ptr<ReplicatedBundleFactory> >& bundle_factories,
                           shared_ptr<TestLauncher> launcher)
{
    for (size_t i = 0; i < bundle_factories.size(); ++i)
    {
        shared_ptr<SampleAbundances> s_ab = shared_ptr<SampleAbundances>(new SampleAbundances);
        
#if ENABLE_THREADS					
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
            
        }
        
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        thread quantitate(sample_worker,
                          boost::ref(rt),
                          boost::ref(*(bundle_factories[i])),
                          s_ab,
                          i,
                          launcher);  
#else
        sample_worker(boost::ref(rt),
                      boost::ref(*(bundle_factories[i])),
                      s_ab,
                      i,
                      launcher);
#endif
    }
    return true;
}

void normalize_as_pool(vector<shared_ptr<ReadGroupProperties> >& all_read_groups)
{
    vector<LocusCountList> sample_count_table;
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        const vector<LocusCount>& raw_count_table = rg_props->raw_counts();
        
        for (size_t j = 0; j < raw_count_table.size(); ++j)
        {
            if (sample_count_table.size() == j)
            {
                const string& locus_id = raw_count_table[j].locus_desc;
                int num_transcripts = raw_count_table[j].num_transcripts;
                sample_count_table.push_back(LocusCountList(locus_id,all_read_groups.size(), num_transcripts));
            }
            double scaled = raw_count_table[j].count;
            //sample_count_table[j].counts[i] = scaled * unscaling_factor;
            sample_count_table[j].counts[i] = scaled;
            assert(sample_count_table[j].counts[i] >= 0 && !isinf(sample_count_table[j].counts[i]));
        }
    }
    
    vector<double> scale_factors(all_read_groups.size(), 0.0);
    
    // TODO: needs to be refactored - similar code exists in replicates.cpp
    calc_scaling_factors(sample_count_table, scale_factors);
    
    for (size_t j = 0; j < scale_factors.size(); ++j)
    {
        double total = 0.0;
        double sf = scale_factors[j];
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            total += sample_count_table[i].counts[j];
        }
        //fprintf(stderr, "SF: %lg, Total: %lg\n", sf, total);
    }
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        rg_props->internal_scale_factor(scale_factors[i]);
        //rg_props->external_scale_factor(scale_factors[i]);
    }
    
    // Transform raw counts to the common scale
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        LocusCountList& p = sample_count_table[i];
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            assert (scale_factors.size() > j);
            p.counts[j] *= (1.0 / scale_factors[j]);
        }
    }
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        vector<LocusCount> scaled_counts;
        for (size_t j = 0; j < sample_count_table.size(); ++j)
        {
            string& locus_id = sample_count_table[j].locus_desc;
            double count = sample_count_table[j].counts[i];
            int num_transcripts = sample_count_table[j].num_transcripts;
            LocusCount locus_count(locus_id, count, num_transcripts);
            scaled_counts.push_back(locus_count);
        }
        rg_props->common_scale_counts(scaled_counts);
        // revert each read group back to native scaling to avoid a systematic fold change toward the mean.
        
        // rg_props->internal_scale_factor(1.0);         
    }
    
    if (poisson_dispersion == false)
    {
        shared_ptr<MassDispersionModel const> disperser;
        disperser = fit_dispersion_model("pooled", scale_factors, sample_count_table, false);
        
        foreach (shared_ptr<ReadGroupProperties> rg_props, all_read_groups)
        {
            rg_props->mass_dispersion_model(disperser);
        }
    }
    
    double avg_total_common_scaled_count = 0.0;
    
    for (size_t fac_idx = 0; fac_idx < all_read_groups.size(); ++fac_idx)
    {
        //shared_ptr<ReadGroupProperties> rg = bundle_factories[fac_idx];
        //double scaled_mass = scale_factors[fac_idx] * rg->total_map_mass();
        double total_common = 0.0;
        for (size_t j = 0; j < sample_count_table.size(); ++j)
        {
            total_common += sample_count_table[j].counts[fac_idx];
        }
        
        avg_total_common_scaled_count += (1.0/all_read_groups.size()) * total_common;
    }
    
    foreach(shared_ptr<ReadGroupProperties> rg, all_read_groups)
    {
        rg->normalized_map_mass(avg_total_common_scaled_count);
        //bf->read_group_properties()->normalized_map_mass(scale_factors[fac_idx]);
    }

}

void fit_isoform_level_count_dispersion(const FPKMTrackingTable& isoform_fpkm_tracking,
                                        vector<shared_ptr<ReplicatedBundleFactory> >& bundle_factories)
{
    map<shared_ptr<MassDispersionModel const>, vector<LocusCountList> > sample_count_table_for_disp_model;
    
    for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
    {
        shared_ptr<ReplicatedBundleFactory> rep_fac = bundle_factories[fac_idx];
        vector<shared_ptr<BundleFactory> > replicates = rep_fac->factories();
        //vector<double> count_table;
        
        vector<LocusCountList> sample_count_table;
        
        set<shared_ptr<ReadGroupProperties> > reps_for_condition;
        
        for (size_t i = 0; i < replicates.size(); ++i)
        {
            shared_ptr<ReadGroupProperties> rg_props = replicates[i]->read_group_properties();
            reps_for_condition.insert(rg_props);
        }
    }

    size_t iso_num = 0;
    
    for (FPKMTrackingTable::const_iterator itr = isoform_fpkm_tracking.begin(); itr != isoform_fpkm_tracking.end(); ++itr)
    {
        const string& description = itr->first;
		const FPKMTracking& track = itr->second;
        const vector<FPKMContext>& fpkms = track.fpkm_series;
		
        bool skip_iso = false;
        for (size_t i = 0; i < fpkms.size(); ++i)
		{
            if (fpkms[i].count_per_rep.empty())
            {
                skip_iso = true;
                break;
            }
        }
        if (skip_iso)
            continue;
        
        for (size_t i = 0; i < fpkms.size(); ++i)
		{
            size_t fac_idx = 0;
            for (CountPerReplicateTable::const_iterator itr = fpkms[i].count_per_rep.begin(); 
                 itr != fpkms[i].count_per_rep.end(); 
                 ++itr)
            {
                shared_ptr<ReadGroupProperties const> rg_props = itr->first;
                pair<map<shared_ptr<MassDispersionModel const>, vector<LocusCountList> >::iterator, bool >  p;
                p = sample_count_table_for_disp_model.insert(make_pair(rg_props->mass_dispersion_model(),
                                                                   vector<LocusCountList>()));
                map<shared_ptr<MassDispersionModel const>, vector<LocusCountList> >::iterator lc_itr = p.first;
                vector<LocusCountList>& lc_list = lc_itr->second;
                double count = itr->second;
                
                if (iso_num >= lc_list.size())
                {
                    LocusCountList locus_count(description, fpkms[i].count_per_rep.size(), 1); 
                    lc_list.push_back(locus_count);
                    lc_list[lc_list.size() - 1].counts[0] = count;
                }
                else
                {
                    const string& ld = lc_list[iso_num].locus_desc;
                    if (ld != description)
                    {
                        fprintf (stderr, "Error: bundle boundaries don't match across replicates!\n");
                        exit(1);
                    }
                    lc_list[iso_num].counts[fac_idx] = count;
                }
                fac_idx++;
            }
        }
        iso_num++;
    }
    
    for (map<shared_ptr<MassDispersionModel const>, vector<LocusCountList> >::const_iterator itr = sample_count_table_for_disp_model.begin(); 
         itr != sample_count_table_for_disp_model.end(); 
         ++itr)
    {
        if (itr->second.empty())
            continue;
        
        vector<double> scale_factors(itr->second.front().counts.size(), 1);
        shared_ptr<MassDispersionModel const> model = fit_dispersion_model(itr->first->name()+"iso", scale_factors, itr->second, true);
        for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
        {
            shared_ptr<ReplicatedBundleFactory> rep_fac = bundle_factories[fac_idx];
            vector<shared_ptr<BundleFactory> > replicates = rep_fac->factories();
            
            for (size_t i = 0; i < replicates.size(); ++i)
            {
                shared_ptr<ReadGroupProperties> rg_props = replicates[i]->read_group_properties();
                if (rg_props->mass_dispersion_model() == itr->first)
                {
                    rg_props->mass_dispersion_model(model);
                }
            }
        }
    }
}

void driver(FILE* ref_gtf, FILE* mask_gtf, vector<string>& sam_hit_filename_lists, Outfiles& outfiles)
{

	ReadTable it;
	RefSequenceTable rt(true, false);
    
	vector<shared_ptr<Scaffold> > ref_mRNAs;
	
	vector<shared_ptr<ReplicatedBundleFactory> > bundle_factories;
    vector<shared_ptr<ReadGroupProperties> > all_read_groups;
    vector<shared_ptr<HitFactory> > all_hit_factories;
    
	for (size_t i = 0; i < sam_hit_filename_lists.size(); ++i)
	{
        vector<string> sam_hit_filenames;
        tokenize(sam_hit_filename_lists[i], ",", sam_hit_filenames);
        
        vector<shared_ptr<BundleFactory> > replicate_factories;
        
        string condition_name = sample_labels[i];
        
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
                            sam_hit_filenames[j].c_str());
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
            
            rg_props->condition_name(condition_name);
            rg_props->replicate_num(j);
            rg_props->file_path(sam_hit_filenames[j]);
            
            all_read_groups.push_back(rg_props);
            
            hf->read_group_properties(rg_props);
            
            replicate_factories.push_back(hf);
            //replicate_factories.back()->set_ref_rnas(ref_mRNAs);
        }
        
        
        bundle_factories.push_back(shared_ptr<ReplicatedBundleFactory>(new ReplicatedBundleFactory(replicate_factories, condition_name)));
	}
    
    ::load_ref_rnas(ref_gtf, rt, ref_mRNAs, corr_bias, false);
    if (ref_mRNAs.empty())
        return;
    
    vector<shared_ptr<Scaffold> > mask_rnas;
    if (mask_gtf)
    {
        ::load_ref_rnas(mask_gtf, rt, mask_rnas, false, false);
    }
    
    foreach (shared_ptr<ReplicatedBundleFactory> fac, bundle_factories)
    {
        fac->set_ref_rnas(ref_mRNAs);
        if (mask_gtf) 
            fac->set_mask_rnas(mask_rnas);
    }
    
#if ENABLE_THREADS
    locus_num_threads = num_threads;
#endif
    
	int tmp_min_frag_len = numeric_limits<int>::max();
	int tmp_max_frag_len = 0;
	
	ProgressBar p_bar("Inspecting maps and determining fragment length distributions.",0);
	foreach (shared_ptr<ReplicatedBundleFactory> fac, bundle_factories)
    {
#if ENABLE_THREADS	
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
        }
        
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        thread inspect(inspect_map_worker,
                       boost::ref(*fac),
                       boost::ref(tmp_min_frag_len),
                       boost::ref(tmp_max_frag_len));  
#else
        inspect_map_worker(boost::ref(*fac),
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
    
    int most_reps = -1;
    int most_reps_idx = 0;
    
    bool single_replicate_fac = false;
    
    for (size_t i = 0; i < bundle_factories.size(); ++i)
    {
        ReplicatedBundleFactory& fac = *(bundle_factories[i]);
        if (fac.num_replicates() > most_reps)
        {
            most_reps = fac.num_replicates();
            most_reps_idx = i;
        }
        if (most_reps == 1)
        {
            single_replicate_fac = true;
        }
    }
    
    if (most_reps == 1)
    {
        normalize_as_pool(all_read_groups);
    }
    
    if (most_reps != 1 && (use_quartile_norm || use_geometric_norm))
    {
        vector<LocusCountList> sample_count_table;
        
        //vector<shared_ptr<ReplicatedBundleFactory> > bundle_factories;
        
        for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
        {
            shared_ptr<ReplicatedBundleFactory> rep_fac = bundle_factories[fac_idx];
            vector<shared_ptr<BundleFactory> > replicates = rep_fac->factories();
            vector<double> count_table;
            for (size_t j = 0; j < replicates.size(); ++j)
            {
                shared_ptr<ReadGroupProperties> rg = replicates[j]->read_group_properties();
                const vector<LocusCount>& rep_count_table = rg->common_scale_counts();
                if (count_table.empty())
                    count_table = vector<double>(rep_count_table.size(), 0);
                
                for (size_t i = 0; i < rep_count_table.size(); ++i)
                {
                    const LocusCount& c = rep_count_table[i];
                    double count = c.count;
                    count_table[i] += (count / replicates.size());;
                }
                
            }
            
            for (size_t i = 0; i < count_table.size(); ++i)
            {
                
                const LocusCount& c = replicates.front()->read_group_properties()->common_scale_counts()[i];
                double count = count_table[i];
                
                if (i >= sample_count_table.size())
                {
                    LocusCountList locus_count(c.locus_desc, bundle_factories.size(), c.num_transcripts); 
                    sample_count_table.push_back(locus_count);
                    sample_count_table.back().counts[0] = count;
                }
                else
                {
                    if (sample_count_table[i].locus_desc != c.locus_desc)
                    {
                        fprintf (stderr, "Error: bundle boundaries don't match across replicates!\n");
                        exit(1);
                    }
                    sample_count_table[i].counts[fac_idx] = count;
                }
                
            }
        }
        
        vector<double> scale_factors(bundle_factories.size(), 0.0);
        
        calc_scaling_factors(sample_count_table, scale_factors);
        
        for (size_t j = 0; j < scale_factors.size(); ++j)
        {
            shared_ptr<ReplicatedBundleFactory> rep_fac = bundle_factories[j];
            vector<shared_ptr<BundleFactory> > replicates = rep_fac->factories();
            
            for (size_t i = 0; i < replicates.size(); ++i)
            {
                shared_ptr<ReadGroupProperties> rg = replicates[i]->read_group_properties();
                rg->external_scale_factor(scale_factors[j]);
            }
            
            double total = 0.0;
            double sf = scale_factors[j];
            for (size_t i = 0; i < sample_count_table.size(); ++i)
            {
                total += sample_count_table[i].counts[j];
            }
            //fprintf(stderr, "SF: %lg, Total: %lg\n", sf, total);
        }
        
        if (use_quartile_norm)
        {
            vector<double> upper_quartiles(bundle_factories.size(), 0);
            vector<double> total_common_masses(bundle_factories.size(), 0);
            
            for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
            {
                //shared_ptr<ReadGroupProperties> rg = bundle_factories[fac_idx];
                //double scaled_mass = scale_factors[fac_idx] * rg->total_map_mass();
                vector<double> common_scaled_counts;
                double total_common = 0.0;
                
                for (size_t j = 0; j < sample_count_table.size(); ++j)
                {
                    total_common += sample_count_table[j].counts[fac_idx];
                    common_scaled_counts.push_back(sample_count_table[j].counts[fac_idx]);
                }
                
                sort(common_scaled_counts.begin(), common_scaled_counts.end());
                if (common_scaled_counts.empty())
                    continue;
                
                int upper_quart_index = common_scaled_counts.size() * 0.75;
                double upper_quart_count = common_scaled_counts[upper_quart_index];
                upper_quartiles[fac_idx] = upper_quart_count;
                total_common_masses[fac_idx] = total_common;
            }
            
            long double total_mass = accumulate(total_common_masses.begin(), total_common_masses.end(), 0.0);
            long double total_norm_mass = accumulate(upper_quartiles.begin(), upper_quartiles.end(), 0.0);
            
            for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
            {
                if (total_mass > 0)
                {
                    //double scaling_factor = total_mass / total_norm_mass;
                    double external_scaling_factor = upper_quartiles[fac_idx] / (total_norm_mass / upper_quartiles.size());
                    foreach(shared_ptr<BundleFactory> bf, bundle_factories[fac_idx]->factories())
                    {
                        //double scaled_mass = scaling_factor * upper_quartiles[fac_idx];
                        bf->read_group_properties()->normalized_map_mass(total_mass / total_common_masses.size());
                        //bf->read_group_properties()->external_scale_factor(1.0);
                        bf->read_group_properties()->external_scale_factor(external_scaling_factor);
                    }
                }
            }
        }
        else
        {
            transform_counts_to_common_scale(scale_factors, sample_count_table);
            
            double avg_total_common_scaled_count = 0.0;
            
            for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
            {
                //shared_ptr<ReadGroupProperties> rg = bundle_factories[fac_idx];
                //double scaled_mass = scale_factors[fac_idx] * rg->total_map_mass();
                double total_common = 0.0;
                for (size_t j = 0; j < sample_count_table.size(); ++j)
                {
                    total_common += sample_count_table[j].counts[fac_idx];
                }
                
                avg_total_common_scaled_count += (1.0/bundle_factories.size()) * total_common;
                
                //rg->normalized_map_mass(scale_factors[fac_idx])
            }
            
            for (size_t fac_idx = 0; fac_idx < bundle_factories.size(); ++fac_idx)
            {
                foreach(shared_ptr<BundleFactory> bf, bundle_factories[fac_idx]->factories())
                {
                    bf->read_group_properties()->normalized_map_mass(avg_total_common_scaled_count);
                }
            }
        }  
        
        foreach (shared_ptr<ReplicatedBundleFactory> fac, bundle_factories)
        {
            // for now, "borrow" the dispersion model for the condition with the most replicates
            size_t borrowed_disp_model_idx = most_reps_idx;
            if (fac->num_replicates() == 1)
            {
                fac->mass_dispersion_model(bundle_factories[borrowed_disp_model_idx]->mass_dispersion_model());
                double borrowed_internal_size_factor = scale_factors[borrowed_disp_model_idx];
                double borrowed_external_size_factor = scale_factors[borrowed_disp_model_idx];
                double borrowed_norm_map_mass = bundle_factories[borrowed_disp_model_idx]->factories().front()->read_group_properties()->normalized_map_mass();
                foreach(shared_ptr<BundleFactory> bf, fac->factories())
                {
                    // we need to adjust the scaling factors so that the FPKMs aren't skewed
                    // and the variance function from the dispersion model is correct.
                    //bf->read_group_properties()->normalized_map_mass(avg_total_common_scaled_count);
                    bf->read_group_properties()->internal_scale_factor(bf->read_group_properties()->external_scale_factor()/borrowed_internal_size_factor);
                    bf->read_group_properties()->normalized_map_mass(borrowed_norm_map_mass);
                    bf->read_group_properties()->external_scale_factor(borrowed_external_size_factor);
                }
            }
        }
    }
    else if (use_raw_mapped_norm)
    {
        // no need to do anything beyond what's already being done during 
        // per-condition map inspection.  Counts are common-scale-transformed 
        // on a per condition basis.  External scale factors are set to 1.0 
        // by default
    }


//    if (use_quartile_norm || use_raw_mapped_norm)
//    {
//        // scale the normalized masses so that both quantile total count normalization
//        // are roughly on the same numerical scale
//        foreach (shared_ptr<ReadGroupProperties> rg_props, all_read_groups)
//        {
//            long double new_norm = rg_props->normalized_map_mass() * (total_mass / total_norm_mass);
//            rg_props->normalized_map_mass(new_norm);
//        }
//    }
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        shared_ptr<ReadGroupProperties> rg = all_read_groups[i];
        fprintf(stderr, "> Map Properties:\n");
        
        fprintf(stderr, ">\tNormalized Map Mass: %.2Lf\n", rg->normalized_map_mass());
        fprintf(stderr, ">\tRaw Map Mass: %.2Lf\n", rg->total_map_mass());
        if (corr_multi)
            fprintf(stderr,">\tNumber of Multi-Reads: %zu (with %zu total hits)\n", rg->multi_read_table()->num_multireads(), rg->multi_read_table()->num_multihits()); 
        
        if (rg->frag_len_dist()->source() == LEARNED)
        {
            fprintf(stderr, ">\tFragment Length Distribution: Empirical (learned)\n");
            fprintf(stderr, ">\t              Estimated Mean: %.2f\n", rg->frag_len_dist()->mean());
            fprintf(stderr, ">\t           Estimated Std Dev: %.2f\n", rg->frag_len_dist()->std_dev());
        }
        else
        {
            if (rg->frag_len_dist()->source() == USER)
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (user-specified)\n");
            else //rg->frag_len_dist()->source == FLD::DEFAULT
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (default)\n");
            fprintf(stderr, ">\t              Default Mean: %d\n", def_frag_len_mean);
            fprintf(stderr, ">\t           Default Std Dev: %d\n", def_frag_len_std_dev);
        }
    }
    
    long double total_norm_mass = 0.0;
    long double total_mass = 0.0;
    foreach (shared_ptr<ReadGroupProperties> rg_props, all_read_groups)
    {
        total_norm_mass += rg_props->normalized_map_mass();
        total_mass += rg_props->total_map_mass();
    }

	min_frag_len = tmp_min_frag_len;
    max_frag_len = tmp_max_frag_len;
	
	final_est_run = false;
	
	double num_bundles = (double)bundle_factories[0]->num_bundles();
	
//    if (corr_bias && corr_multi)
//        p_bar = ProgressBar("Calculating initial abundance estimates for bias and multi-read correction.", num_bundles);
//    else if (corr_bias)
//        p_bar = ProgressBar("Calculating initial abundance estimates for bias correction.", num_bundles);
//    else if (corr_multi)
//        p_bar = ProgressBar("Calculating initial abundance estimates for multi-read correction.", num_bundles);
//    else
    
    p_bar = ProgressBar("Calculating preliminary abundance estimates", num_bundles);
    
    Tracking tracking;
    
    test_launcher = shared_ptr<TestLauncher>(new TestLauncher(bundle_factories.size(), NULL, &tracking, samples_are_time_series, &p_bar));
    
	if (corr_bias || corr_multi) // Only run initial estimation if correcting bias or multi-reads
	{
		while (1) 
		{
			//p_bar.update("",1);
            //test_launcher = shared_ptr<TestLauncher>(new TestLauncher((int)bundle_factories.size(), NULL, &tracking, samples_are_time_series, &p_bar));
                                                     
			shared_ptr<vector<shared_ptr<SampleAbundances> > > abundances(new vector<shared_ptr<SampleAbundances> >());
			quantitate_next_locus(rt, bundle_factories, test_launcher);
			bool more_loci_remain = false;
            foreach (shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories) 
            {
                if (rep_fac->bundles_remain())
                {
                    more_loci_remain = true;
                    break;
                }
            }
            
			if (!more_loci_remain)
            {
                // wait for the workers to finish up before breaking out.
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
				break;
            }
		}
        
        foreach (shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories)
		{
			rep_fac->reset();
        }
        
		p_bar.complete();
	}
    if (corr_bias)
    {
        bias_run = true;
        p_bar = ProgressBar("Learning bias parameters.", 0);
		foreach (shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories)
		{
			foreach (shared_ptr<BundleFactory> fac, rep_fac->factories())
			{
#if ENABLE_THREADS	
				while(1)
				{
					locus_thread_pool_lock.lock();
					if (locus_curr_threads < locus_num_threads)
					{
						break;
					}
					
					locus_thread_pool_lock.unlock();
					
					boost::this_thread::sleep(boost::posix_time::milliseconds(5));
				}
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
        foreach (shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories)
		{
			rep_fac->reset();
        }
        bias_run = false;
	}
    
//    if (use_isoform_count_dispersion)
//    {
//        fit_isoform_level_count_dispersion(tracking.isoform_fpkm_tracking, bundle_factories);
//    }
    
    // Allow the multiread tables to do their thing...
    foreach (shared_ptr<ReadGroupProperties> rg_props, all_read_groups)
    {
        rg_props->multi_read_table()->valid_mass(true);
    }
    
    test_launcher->clear_tracking_data();
	
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
	
	final_est_run = true;
	p_bar = ProgressBar("Testing for differential expression and regulation in locus.", num_bundles);
                                                     
    test_launcher = shared_ptr<TestLauncher>(new TestLauncher(bundle_factories.size(), &tests, &tracking, samples_are_time_series, &p_bar));
                                                                                              
	while (true)
	{
        //shared_ptr<vector<shared_ptr<SampleAbundances> > > abundances(new vector<shared_ptr<SampleAbundances> >());
        quantitate_next_locus(rt, bundle_factories, test_launcher);
        bool more_loci_remain = false;
        foreach (shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories) 
        {
            if (rep_fac->bundles_remain())
            {
                more_loci_remain = true;
                break;
            }
        }
        if (!more_loci_remain)
        {
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
            break;
        }
    }
	
	p_bar.complete();

	//double FDR = 0.05;
	int total_iso_de_tests = 0;
	
	vector<SampleDifference*> isoform_exp_diffs;
    fprintf(outfiles.isoform_de_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    //if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.isoform_de_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.isoform_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.isoform_de_tests[i][j]);
            }
        }
	}
    
	int total_group_de_tests = 0;
	vector<SampleDifference*> tss_group_exp_diffs;
    fprintf(outfiles.group_de_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
	//if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.tss_group_de_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.group_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.tss_group_de_tests[i][j]);
            }
        }
    }
	
	int total_gene_de_tests = 0;
	vector<SampleDifference*> gene_exp_diffs;
    fprintf(outfiles.gene_de_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    if (no_differential == false)
    {
        for (size_t i = 1; i < tests.gene_de_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                total_gene_de_tests += tests.gene_de_tests[i][j].size();
                extract_sample_diffs(tests.gene_de_tests[i][j], gene_exp_diffs);
            }
        }
    
        //fprintf(stderr, "***There are %lu difference records in gene_exp_diffs\n", gene_exp_diffs.size());
        int gene_exp_tests = fdr_significance(FDR, gene_exp_diffs);
        fprintf(stderr, "Performed %d gene-level transcription difference tests\n", gene_exp_tests);
        
        for (size_t i = 1; i < tests.gene_de_tests.size(); ++i)
        {        
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.gene_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.gene_de_tests[i][j]);
            }
        }	
    }
    
	int total_cds_de_tests = 0;
	vector<SampleDifference*> cds_exp_diffs;
    fprintf(outfiles.cds_de_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    //if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.cds_de_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.cds_de_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.cds_de_tests[i][j]);
            }
        }
	}
    
	int total_diff_splice_tests = 0;
	vector<SampleDifference*> splicing_diffs;
    fprintf(outfiles.diff_splicing_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    //if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.diff_splicing_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                const SampleDiffs& diffs = tests.diff_splicing_tests[i][j];
                print_tests(outfiles.diff_splicing_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), diffs);
            }
        }
	}
    
	int total_diff_promoter_tests = 0;
	vector<SampleDifference*> promoter_diffs;
    fprintf(outfiles.diff_promoter_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    //if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.diff_promoter_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.diff_promoter_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.diff_promoter_tests[i][j]);
            }
        }
    }
    
	int total_diff_cds_tests = 0;
	vector<SampleDifference*> cds_use_diffs;
    fprintf(outfiles.diff_cds_outfile, "test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tsqrt(JS)\ttest_stat\tp_value\tq_value\tsignificant\n");
    
    //if (no_differential == false)
    {
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
        
        for (size_t i = 1; i < tests.diff_cds_tests.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                print_tests(outfiles.diff_cds_outfile, sample_labels[j].c_str(), sample_labels[i].c_str(), tests.diff_cds_tests[i][j]);
            }
        }
    }
	
    // FPKM tracking
    
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

    // Count tracking
    
    FILE* fiso_count_tracking =  outfiles.isoform_count_tracking_out;
	fprintf(stderr, "Writing isoform-level count tracking\n");
	print_count_tracking(fiso_count_tracking,tracking.isoform_fpkm_tracking); 
	
	FILE* ftss_count_tracking =  outfiles.tss_group_count_tracking_out;
	fprintf(stderr, "Writing TSS group-level count tracking\n");
	print_count_tracking(ftss_count_tracking,tracking.tss_group_fpkm_tracking);
	
	FILE* fgene_count_tracking =  outfiles.gene_count_tracking_out;
	fprintf(stderr, "Writing gene-level count tracking\n");
	print_count_tracking(fgene_count_tracking,tracking.gene_fpkm_tracking);
	
	FILE* fcds_count_tracking =  outfiles.cds_count_tracking_out;
	fprintf(stderr, "Writing CDS-level count tracking\n");
	print_count_tracking(fcds_count_tracking,tracking.cds_fpkm_tracking);
    
    // Read group tracking
    
    FILE* fiso_rep_tracking =  outfiles.isoform_rep_tracking_out;
	fprintf(stderr, "Writing isoform-level read group tracking\n");
	print_read_group_tracking(fiso_rep_tracking,tracking.isoform_fpkm_tracking); 
	
	FILE* ftss_rep_tracking =  outfiles.tss_group_rep_tracking_out;
	fprintf(stderr, "Writing TSS group-level read group tracking\n");
	print_read_group_tracking(ftss_rep_tracking,tracking.tss_group_fpkm_tracking);
	
	FILE* fgene_rep_tracking =  outfiles.gene_rep_tracking_out;
	fprintf(stderr, "Writing gene-level read group tracking\n");
	print_read_group_tracking(fgene_rep_tracking,tracking.gene_fpkm_tracking);
	
	FILE* fcds_rep_tracking =  outfiles.cds_rep_tracking_out;
	fprintf(stderr, "Writing CDS-level read group tracking\n");
	print_read_group_tracking(fcds_rep_tracking,tracking.cds_fpkm_tracking);
    
    FILE* fread_group_info =  outfiles.read_group_info_out;
	fprintf(stderr, "Writing read group info\n");
	print_read_group_info(fread_group_info,all_read_groups);

    FILE* frun_info =  outfiles.run_info_out;
	fprintf(stderr, "Writing run info\n");
	print_run_info(frun_info);
}

int main(int argc, char** argv)
{
    for (int i = 0; i < argc; ++i)
    {
        cmd_str += string(argv[i]) + " ";
    }
    
    init_library_table();
    
    min_isoform_fraction = 1e-5;
    
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
    if (!use_total_mass && !use_compat_mass)
    {
        use_total_mass = false;
        use_compat_mass = true;   
    }
    
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    if (!no_update_check)
        check_version(PACKAGE_VERSION);
        
    
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
    
    if (random_seed == -1)
        random_seed = time(NULL);
    
	// seed the random number generator - we'll need it for the importance
	// sampling during MAP estimation of the gammas
	srand48(random_seed);
	
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
    
    static const int filename_buf_size = 2048;
    
    char out_file_prefix[filename_buf_size];
    sprintf(out_file_prefix, "%s/", output_dir.c_str());
    char iso_out_file_name[filename_buf_size];
    sprintf(iso_out_file_name, "%sisoform_exp.diff", out_file_prefix);
    FILE* iso_out = fopen(iso_out_file_name, "w");
    if (!iso_out)
    {
        fprintf(stderr, "Error: cannot open differential isoform transcription file %s for writing\n",
                iso_out_file_name);
        exit(1);
    }
    
    char group_out_file_name[filename_buf_size];
    sprintf(group_out_file_name, "%stss_group_exp.diff", out_file_prefix);
    FILE* group_out = fopen(group_out_file_name, "w");
    if (!group_out)
    {
        fprintf(stderr, "Error: cannot open differential TSS group transcription file %s for writing\n",
                group_out_file_name);
        exit(1);
    }
    
    char gene_out_file_name[filename_buf_size];
    sprintf(gene_out_file_name, "%sgene_exp.diff", out_file_prefix);
    FILE* gene_out = fopen(gene_out_file_name, "w");
    if (!group_out)
    {
        fprintf(stderr, "Error: cannot open gene expression file %s for writing\n",
                gene_out_file_name);
        exit(1);
    }
    
    char cds_out_file_name[filename_buf_size];
    sprintf(cds_out_file_name, "%scds_exp.diff", out_file_prefix);
    FILE* cds_out = fopen(cds_out_file_name, "w");
    if (!cds_out)
    {
        fprintf(stderr, "Error: cannot open cds expression file %s for writing\n",
                cds_out_file_name);
        exit(1);
    }
    
    char diff_splicing_out_file_name[filename_buf_size];
    sprintf(diff_splicing_out_file_name, "%ssplicing.diff", out_file_prefix);
    FILE* diff_splicing_out = fopen(diff_splicing_out_file_name, "w");
    if (!diff_splicing_out)
    {
        fprintf(stderr, "Error: cannot open differential splicing file %s for writing\n",
                diff_splicing_out_file_name);
        exit(1);
    }
    
    char diff_promoter_out_file_name[filename_buf_size];
    sprintf(diff_promoter_out_file_name, "%spromoters.diff", out_file_prefix);
    FILE* diff_promoter_out = fopen(diff_promoter_out_file_name, "w");
    if (!diff_promoter_out)
    {
        fprintf(stderr, "Error: cannot open differential transcription start file %s for writing\n",
                diff_promoter_out_file_name);
        exit(1);
    }
    
    char diff_cds_out_file_name[filename_buf_size];
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
	
	char isoform_fpkm_tracking_name[filename_buf_size];
	sprintf(isoform_fpkm_tracking_name, "%s/isoforms.fpkm_tracking", output_dir.c_str());
	FILE* isoform_fpkm_out = fopen(isoform_fpkm_tracking_name, "w");
	if (!isoform_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level FPKM tracking file %s for writing\n",
				isoform_fpkm_tracking_name);
		exit(1);
	}
	outfiles.isoform_fpkm_tracking_out = isoform_fpkm_out;

	char tss_group_fpkm_tracking_name[filename_buf_size];
	sprintf(tss_group_fpkm_tracking_name, "%s/tss_groups.fpkm_tracking", output_dir.c_str());
	FILE* tss_group_fpkm_out = fopen(tss_group_fpkm_tracking_name, "w");
	if (!tss_group_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level FPKM tracking file %s for writing\n",
				tss_group_fpkm_tracking_name);
		exit(1);
	}
	outfiles.tss_group_fpkm_tracking_out = tss_group_fpkm_out;

	char cds_fpkm_tracking_name[filename_buf_size];
	sprintf(cds_fpkm_tracking_name, "%s/cds.fpkm_tracking", output_dir.c_str());
	FILE* cds_fpkm_out = fopen(cds_fpkm_tracking_name, "w");
	if (!cds_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open CDS level FPKM tracking file %s for writing\n",
				cds_fpkm_tracking_name);
		exit(1);
	}
	outfiles.cds_fpkm_tracking_out = cds_fpkm_out;
	
	char gene_fpkm_tracking_name[filename_buf_size];
	sprintf(gene_fpkm_tracking_name, "%s/genes.fpkm_tracking", output_dir.c_str());
	FILE* gene_fpkm_out = fopen(gene_fpkm_tracking_name, "w");
	if (!gene_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open gene-level FPKM tracking file %s for writing\n",
				gene_fpkm_tracking_name);
		exit(1);
	}
	outfiles.gene_fpkm_tracking_out = gene_fpkm_out;

    char isoform_count_tracking_name[filename_buf_size];
	sprintf(isoform_count_tracking_name, "%s/isoforms.count_tracking", output_dir.c_str());
	FILE* isoform_count_out = fopen(isoform_count_tracking_name, "w");
	if (!isoform_count_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level count tracking file %s for writing\n",
				isoform_count_tracking_name);
		exit(1);
	}
	outfiles.isoform_count_tracking_out = isoform_count_out;
    
	char tss_group_count_tracking_name[filename_buf_size];
	sprintf(tss_group_count_tracking_name, "%s/tss_groups.count_tracking", output_dir.c_str());
	FILE* tss_group_count_out = fopen(tss_group_count_tracking_name, "w");
	if (!tss_group_count_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level count tracking file %s for writing\n",
				tss_group_count_tracking_name);
		exit(1);
	}
	outfiles.tss_group_count_tracking_out = tss_group_count_out;
    
	char cds_count_tracking_name[filename_buf_size];
	sprintf(cds_count_tracking_name, "%s/cds.count_tracking", output_dir.c_str());
	FILE* cds_count_out = fopen(cds_count_tracking_name, "w");
	if (!cds_count_out)
	{
		fprintf(stderr, "Error: cannot open CDS level count tracking file %s for writing\n",
				cds_count_tracking_name);
		exit(1);
	}
	outfiles.cds_count_tracking_out = cds_count_out;
	
	char gene_count_tracking_name[filename_buf_size];
	sprintf(gene_count_tracking_name, "%s/genes.count_tracking", output_dir.c_str());
	FILE* gene_count_out = fopen(gene_count_tracking_name, "w");
	if (!gene_count_out)
	{
		fprintf(stderr, "Error: cannot open gene-level count tracking file %s for writing\n",
				gene_count_tracking_name);
		exit(1);
	}
	outfiles.gene_count_tracking_out = gene_count_out;
    
    char isoform_rep_tracking_name[filename_buf_size];
	sprintf(isoform_rep_tracking_name, "%s/isoforms.read_group_tracking", output_dir.c_str());
	FILE* isoform_rep_out = fopen(isoform_rep_tracking_name, "w");
	if (!isoform_rep_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level read group tracking file %s for writing\n",
				isoform_rep_tracking_name);
		exit(1);
	}
	outfiles.isoform_rep_tracking_out = isoform_rep_out;
    
	char tss_group_rep_tracking_name[filename_buf_size];
	sprintf(tss_group_rep_tracking_name, "%s/tss_groups.read_group_tracking", output_dir.c_str());
	FILE* tss_group_rep_out = fopen(tss_group_rep_tracking_name, "w");
	if (!tss_group_rep_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level read group tracking file %s for writing\n",
				tss_group_rep_tracking_name);
		exit(1);
	}
	outfiles.tss_group_rep_tracking_out = tss_group_rep_out;
    
	char cds_rep_tracking_name[filename_buf_size];
	sprintf(cds_rep_tracking_name, "%s/cds.read_group_tracking", output_dir.c_str());
	FILE* cds_rep_out = fopen(cds_rep_tracking_name, "w");
	if (!cds_rep_out)
	{
		fprintf(stderr, "Error: cannot open CDS level read group tracking file %s for writing\n",
				cds_rep_tracking_name);
		exit(1);
	}
	outfiles.cds_rep_tracking_out = cds_rep_out;
	
	char gene_rep_tracking_name[filename_buf_size];
	sprintf(gene_rep_tracking_name, "%s/genes.read_group_tracking", output_dir.c_str());
	FILE* gene_rep_out = fopen(gene_rep_tracking_name, "w");
	if (!gene_rep_out)
	{
		fprintf(stderr, "Error: cannot open gene-level read group tracking file %s for writing\n",
				gene_rep_tracking_name);
		exit(1);
	}
	outfiles.gene_rep_tracking_out = gene_rep_out;
    
    char read_group_info_name[filename_buf_size];
	sprintf(read_group_info_name, "%s/read_groups.info", output_dir.c_str());
	FILE* read_group_out = fopen(read_group_info_name, "w");
	if (!read_group_out)
	{
		fprintf(stderr, "Error: cannot open read group info file %s for writing\n",
				read_group_info_name);
		exit(1);
	}
	outfiles.read_group_info_out = read_group_out;
    
    char run_info_name[filename_buf_size];
	sprintf(run_info_name, "%s/run.info", output_dir.c_str());
	FILE* run_info_out = fopen(run_info_name, "w");
	if (!run_info_out)
	{
		fprintf(stderr, "Error: cannot open run info file %s for writing\n",
				run_info_name);
		exit(1);
	}
	outfiles.run_info_out = run_info_out;
    
    driver(ref_gtf, mask_gtf, sam_hit_filenames, outfiles);
	
#if 0
    if (emit_count_tables)
    {
        dump_locus_variance_info(output_dir + string("/locus_var.txt"));
    }
#endif
    
	return 0;
}

