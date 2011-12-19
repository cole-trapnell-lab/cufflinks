/*
 *  cufflinks.cpp
 *  Cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
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

#include "common.h"
#include "hits.h"

#include <boost/thread.hpp>

#include "update_check.h"
#include "clustering.h"
#include "abundances.h"
#include "bundles.h"
#include "filters.h"
#include "genes.h"
#include "assemble.h"
#include "biascorrection.h"
#include "multireads.h"

using namespace std;

#if ENABLE_THREADS
const char *short_options = "m:p:s:F:I:j:Q:L:G:g:o:M:b:a:A:Nqvu";
#else
const char *short_options = "m:s:F:I:j:Q:L:G:g:o:M:b:a:A:Nqvu";
#endif

static struct option long_options[] = {
// general options
{"GTF",					    required_argument,		 0,			 'G'},
{"GTF-guide",			    required_argument,		 0,			 'g'},
{"mask-gtf",                required_argument,		 0,			 'M'},
{"library-type",		    required_argument,		 0,			 OPT_LIBRARY_TYPE},
{"seed",                    required_argument,		 0,			 OPT_RANDOM_SEED},

// program behavior
{"output-dir",			    required_argument,		 0,			 'o'},
{"verbose",			    	no_argument,		 	 0,			 'v'},
{"quiet",			    	no_argument,			 0,			 'q'},
{"no-update-check",         no_argument,             0,          OPT_NO_UPDATE_CHECK},
#if ENABLE_THREADS
    {"num-threads",			required_argument,       0,          'p'},
#endif    
{"output-fld",              no_argument,             0,          OPT_OUTPUT_FLD},
{"output-bias-params",      no_argument,             0,          OPT_OUTPUT_BIAS_PARAMS},
    
// abundance estimation
{"frag-len-mean",			required_argument,       0,          'm'},
{"frag-len-std-dev",		required_argument,       0,          's'},
{"min-isoform-fraction",    required_argument,       0,          'F'},
{"upper-quartile-normalization",  no_argument,	 		 0,	         'N'},
{"frag-bias-correct",       required_argument,		 0,			 'b'},
{"multi-read-correct",      no_argument,			 0,			 'u'},

{"num-importance-samples",  required_argument,		 0,			 OPT_NUM_IMP_SAMPLES},
{"max-mle-iterations",		required_argument,		 0,			 OPT_MLE_MAX_ITER},
{"bias-mode",               required_argument,		 0,			 OPT_BIAS_MODE},
{"use-grad-ascent",         no_argument,             0,			 OPT_USE_EM},
{"no-collapse-cond-prob",   no_argument,             0,			 OPT_COLLAPSE_COND_PROB},
{"compatible-hits-norm",    no_argument,	 		 0,	         OPT_USE_COMPAT_MASS},
{"total-hits-norm",         no_argument,	 		 0,	         OPT_USE_TOTAL_MASS},
    
// assembly
{"pre-mrna-fraction",		required_argument,		 0,			 'j'},
{"junc-alpha",				required_argument,		 0,			 'a'},	
{"small-anchor-fraction",	required_argument,		 0,			 'A'},
{"max-intron-length",		required_argument,		 0,			 'I'},
{"label",					required_argument,		 0,			 'L'},
{"overhang-tolerance",      required_argument,		 0,			 OPT_OVERHANG_TOLERANCE},
{"min-frags-per-transfrag",required_argument,		 0,			 OPT_MIN_FRAGS_PER_TRANSFRAG},
{"min-intron-length",       required_argument,	     0,			 OPT_MIN_INTRON_LENGTH},
{"max-bundle-length",       required_argument,		 0,			 OPT_MAX_BUNDLE_LENGTH},
{"trim-3-dropoff-frac",     required_argument,		 0,			 OPT_3_PRIME_AVGCOV_THRESH},
{"trim-3-avgcov-thresh",	required_argument,		 0,			 OPT_3_PRIME_AVGCOV_THRESH},
    
{"3-overhang-tolerance",	required_argument,		 0,			 OPT_3_OVERHANG_TOLERANCE},
{"intron-overhang-tolerance",	required_argument,		 0,		 OPT_INTRON_OVERHANG_TOLERANCE},
{"no-faux-reads",           no_argument,             0,          OPT_NO_FAUX_READS},
{"no-5-extend",             no_argument,             0,          OPT_NO_5_EXTEND},
{"tile-read-len",           required_argument,       0,          OPT_TILE_LEN}, 
{"tile-read-sep",           required_argument,       0,          OPT_TILE_SEP}, 
    
{"max-bundle-frags",        required_argument,        0,          OPT_MAX_FRAGS_PER_BUNDLE}, 
{0, 0, 0, 0} // terminator
};

void print_usage()
{
    //NOTE: SPACES ONLY, bozo
    fprintf(stderr, "cufflinks v%s\n", PACKAGE_VERSION);
    fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
    fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "Usage:   cufflinks [options] <hits.sam>\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "  -o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");
#if ENABLE_THREADS
    fprintf(stderr, "  -p/--num-threads             number of threads used during analysis                [ default:      1 ]\n");
#endif    
    fprintf(stderr, "  --seed                       value of random number generator seed                 [ default:      0 ]\n");
    fprintf(stderr, "  -G/--GTF                     quantitate against reference transcript annotations                      \n");
    fprintf(stderr, "  -g/--GTF-guide               use reference transcript annotation to guide assembly                   \n");
    fprintf(stderr, "  -M/--mask-file               ignore all alignment within transcripts in this file                     \n");
    fprintf(stderr, "  -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]\n");
    fprintf(stderr, "  -u/--multi-read-correct      use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]\n");
    fprintf(stderr, "  --library-type               library prep used for input reads                     [ default:  below ]\n");
    
    fprintf(stderr, "\nAdvanced Abundance Estimation Options:\n");
    fprintf(stderr, "  -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]\n");
    fprintf(stderr, "  -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]\n");
    fprintf(stderr, "  --upper-quartile-norm        use upper-quartile normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, "  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]\n");
    fprintf(stderr, "  --num-importance-samples     number of importance samples for MAP restimation      [ default:   1000 ]\n");
    fprintf(stderr, "  --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:  FALSE ]\n");
    fprintf(stderr, "  --total-hits-norm            count all hits for normalization                      [ default:  TRUE  ]\n");
    
    fprintf(stderr, "\nAdvanced Assembly Options:\n");
    fprintf(stderr, "  -L/--label                   assembled transcripts have this ID prefix             [ default:   CUFF ]\n");
    fprintf(stderr, "  -F/--min-isoform-fraction    suppress transcripts below this abundance level       [ default:   0.10 ]\n");
    fprintf(stderr, "  -j/--pre-mrna-fraction       suppress intra-intronic transcripts below this level  [ default:   0.15 ]\n");
    fprintf(stderr, "  -I/--max-intron-length       ignore alignments with gaps longer than this          [ default: 300000 ]\n");
    fprintf(stderr, "  -a/--junc-alpha              alpha for junction binomial test filter               [ default:  0.001 ]\n");
    fprintf(stderr, "  -A/--small-anchor-fraction   percent read overhang taken as 'suspiciously small'   [ default:   0.09 ]\n");
    fprintf(stderr, "  --min-frags-per-transfrag    minimum number of fragments needed for new transfrags [ default:     10 ]\n");
    fprintf(stderr, "  --overhang-tolerance         number of terminal exon bp to tolerate in introns     [ default:      8 ]\n");
    fprintf(stderr, "  --max-bundle-length          maximum genomic length allowed for a given bundle     [ default:3500000 ]\n");
    fprintf(stderr, "  --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]\n");
    fprintf(stderr, "  --min-intron-length          minimum intron size allowed in genome                 [ default:     50 ]\n");
    fprintf(stderr, "  --trim-3-avgcov-thresh       minimum avg coverage required to attempt 3' trimming  [ default:     10 ]\n");
    fprintf(stderr, "  --trim-3-dropoff-frac        fraction of avg coverage below which to trim 3' end   [ default:    0.1 ]\n");
    
    fprintf(stderr, "\nAdvanced Reference Annotation Guided Assembly Options:\n");
//    fprintf(stderr, "  --tile-read-len              length of faux-reads                                  [ default:    405 ]\n");
//    fprintf(stderr, "  --tile-read-sep              distance between faux-reads                           [ default:     15 ]\n");
    fprintf(stderr, "  --no-faux-reads              disable tiling by faux reads                          [ default:  FALSE ]\n");
    fprintf(stderr, "  --3-overhang-tolerance       overhang allowed on 3' end when merging with reference[ default:    600 ]\n");
    fprintf(stderr, "  --intron-overhang-tolerance  overhang allowed inside reference intron when merging [ default:     30 ]\n");
    
    fprintf(stderr, "\nAdvanced Program Behavior Options:\n");
    fprintf(stderr, "  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
    fprintf(stderr, "  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
    fprintf(stderr, "  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]\n");
    print_library_table();
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
	bool F_set = false;
	
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'm':
				user_provided_fld = true;
				def_frag_len_mean = (uint32_t)parseInt(0, "-m/--frag-len-mean arg must be at least 0", print_usage);
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
				F_set = true;
				break;
			case 'I':
				max_intron_length = parseInt(1, "-I/--max-intron-length must be at least 1", print_usage);
				break;
			case 'j':
				pre_mrna_fraction = parseFloat(0, 1.0, "-I/--pre-mrna-fraction must be at least 0", print_usage);
				break;
				
			case 'a':
				binomial_junc_filter_alpha = parseFloat(0, 1.0, "-a/--junc-alpha must be  between 0 and 1.0", print_usage);
				break;
			case 'A':
				small_anchor_fraction = parseFloat(0, 1.0, "-A/--small-anchor-fraction must be  between 0 and 1.0", print_usage);
				break;
            case OPT_OVERHANG_TOLERANCE:
				bowtie_overhang_tolerance = parseInt(0, "--overhang-tolerance must be at least 0", print_usage);
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
			case 'L':
			{
				user_label = optarg;
				break;
			}
			case 'G':
			{
				ref_gtf_filename = optarg;
				bundle_mode = REF_DRIVEN;
                init_bundle_mode = REF_DRIVEN;
				break;
			}
			case 'g':
			{
				ref_gtf_filename = optarg;
				bundle_mode = REF_GUIDED;
                init_bundle_mode = REF_GUIDED;
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
			case 'N':
            {
            	use_quartile_norm = true;
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
            case OPT_MAX_BUNDLE_LENGTH:
			{
				max_gene_length = parseInt(1, "--max-bundle-length must be at least 1", print_usage);;
				break;
			}
            case OPT_MIN_FRAGS_PER_TRANSFRAG:
			{
				min_frags_per_transfrag = parseInt(0, "--min-frags-per-transfrag must be at least 0", print_usage);;
				break;
			}
            case OPT_MIN_INTRON_LENGTH:
			{
				min_intron_length = parseInt(0, "--min-intron-length must be at least 0", print_usage);
				break;
			}
			case OPT_3_PRIME_AVGCOV_THRESH:
			{
				trim_3_avgcov_thresh = parseFloat(0, 9999999, "--trim-3-avgcov-thresh must be at least 0", print_usage);
				break;
			}
            case OPT_3_PRIME_DROPOFF_FRAC:
			{
				trim_3_dropoff_frac = parseFloat(0, 1.0, "--trim-3-dropoff-frac must be between 0 and 1.0", print_usage);
				break;
			}
            case OPT_NO_UPDATE_CHECK:
            {
                no_update_check = true;
                break;
            }
            case OPT_OUTPUT_FLD:
            {
                output_fld = true;
                break;
            }
            case OPT_OUTPUT_BIAS_PARAMS:
            {
                output_bias_params = true;
                break;
            }
            case OPT_USE_EM:
            {
                use_em = false;
                break;
            }
            case OPT_COLLAPSE_COND_PROB:
            {
                cond_prob_collapse = false;
                break;
            }
            case OPT_NO_FAUX_READS:
            {
                enable_faux_reads = false;
                break;
            }
            case OPT_NO_5_EXTEND:
            {
                enable_5_extend = false;
                break;
            }
            case OPT_3_OVERHANG_TOLERANCE:
            {
                overhang_3 = parseInt(0, "--3-overhang-tolernace must be at least 0", print_usage);
                break;
            }
            case OPT_TILE_LEN:
            {
                tile_len = parseInt(0, "--tile-read-len must be at least 0", print_usage);
                break;
            }
            case OPT_TILE_SEP:
            {
                tile_off = parseInt(0, "--tile-read-sep must be at least 0", print_usage);
                break;
            }
            case OPT_INTRON_OVERHANG_TOLERANCE:
            {
                ref_merge_overhang_tolerance = parseInt(0, "--intron-overhang-tolerance must be at least 0", print_usage);
                break;
            }
            case OPT_RANDOM_SEED:
            {
                random_seed = parseInt(0, "--seed must be at least 0", print_usage);
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
            case OPT_MAX_FRAGS_PER_BUNDLE:
            {
                max_frags_per_bundle = parseInt(0, "--max-bundle-frags must be at least 0", print_usage);
                break;
            }
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    
	if (bundle_mode == REF_DRIVEN)
	{
        if (!F_set)
        {
            min_isoform_fraction = 0.0;
        }
	}
	
	if (bundle_mode == REF_DRIVEN)
	{
		allow_junk_filtering = false;
	}
	
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
//            if (library_type == "transfrags")
//            {
//                allow_junk_filtering = false;
//            }
            global_read_properties = &lib_itr->second;
        }
    }
    
    if (use_total_mass && use_compat_mass)
    {
        fprintf (stderr, "Error: please supply only one of --compatibile-hits-norm and --total-hits-norm\n");
        exit(1);
    }
    if (use_compat_mass && bundle_mode != REF_DRIVEN)
    {
        fprintf (stderr, "Error: cannot use --compatible-hits-norm without --GTF\n");
        exit(1);
    }
	
    return 0;
}

void combine_strand_assemblies(vector<Scaffold>& lhs, 
						   vector<Scaffold>& rhs,
						   vector<Scaffold>& scaffolds,
						   vector<shared_ptr<Scaffold> >* ref_scaffs)
{
	// first check for strand support
    for (size_t l = 0; l < lhs.size(); ++l)
    {
		if (!lhs[l].has_strand_support(ref_scaffs))
			lhs[l].strand(CUFF_STRAND_UNKNOWN);
	}
	for (size_t r = 0; r < rhs.size(); ++r)
    {
		if (!rhs[r].has_strand_support(ref_scaffs))
			rhs[r].strand(CUFF_STRAND_UNKNOWN);
	}
	
    vector<bool> kept_lhs(lhs.size(), true);
    vector<bool> kept_rhs(rhs.size(), true);
    
	// next filter both lists based on reference transcripts (if available)
	if (ref_scaffs != NULL)
	{
		for(size_t l = 0; l < lhs.size(); ++l)
		{			
			foreach(shared_ptr<Scaffold> ref_scaff, *ref_scaffs)
			{
                // if we're past all the overlaps, just stop
				if (ref_scaff->left() >= lhs[l].right() + overhang_3)
				{
					//break;
				}
                // don't emit assembled transfrags that are contained within reference ones
				else if (ref_scaff->contains(lhs[l], 0, overhang_3) && Scaffold::compatible(*ref_scaff, lhs[l], ref_merge_overhang_tolerance))
				{
					kept_lhs[l] = false;
				}
                // if they're compatible but not equal, let's check a few more criteria before 
                // we decide to emit the assembled guy
				else if (ref_scaff->overlapped_3(lhs[l], 0, overhang_3) && Scaffold::compatible(*ref_scaff, lhs[l], ref_merge_overhang_tolerance))
				{
					if (ref_scaff->gaps() == lhs[l].gaps())
                    {
                        kept_lhs[l] = false;
                    }
                    else
                    {
//                        if (enable_5_extend)
//                        {
//                            ref_scaff->extend_5(lhs[l]);
//                            kept_lhs[l] = false;
//                        }
                    }
				}
			}
		}
		for(size_t r = 0; r < rhs.size(); ++r)
		{			
			foreach(shared_ptr<Scaffold> ref_scaff, *ref_scaffs)
			{
				if (ref_scaff->left() >= rhs[r].right() + overhang_3)
				{
					//break;
				}
				else if (ref_scaff->contains(rhs[r], 0, overhang_3) && Scaffold::compatible(*ref_scaff, rhs[r], ref_merge_overhang_tolerance))
				{
					kept_rhs[r] = false;
				}
				else if (ref_scaff->overlapped_3(rhs[r], 0, overhang_3) && Scaffold::compatible(*ref_scaff, rhs[r], ref_merge_overhang_tolerance))
				{
                    if (ref_scaff->gaps() == rhs[r].gaps())
                    {
                        kept_rhs[r] = false;
                    }
                    else
                    {
//                        if (enable_5_extend)
//                        {
//                            ref_scaff->extend_5(rhs[r]);
//                            kept_rhs[r] = false;
//                        }
                    }
				}
			}
		}
	}
	
    // We want to keep all fwd, all reverse, and only the non-redundant unknowns
    // if two unknown strand frags olap, merge them.
    for (size_t l = 0; l < lhs.size(); ++l)
    {		 
		if (!kept_lhs[l])
			continue;
		bool lhs_support = (lhs[l].strand() != CUFF_STRAND_UNKNOWN);
        
		for (size_t r = 0; r < rhs.size(); ++r)
		{
			if (!kept_rhs[r])
				continue;
			if (Scaffold::overlap_in_genome(lhs[l], rhs[r], 0))
			{
				if (Scaffold::compatible(lhs[l], rhs[r]))
				{
					bool rhs_support = (rhs[r].strand() != CUFF_STRAND_UNKNOWN);
					if (!lhs_support && !rhs_support)
					{
						Scaffold merged;
						Scaffold::merge(lhs[l],rhs[r],merged, true);
						scaffolds.push_back(merged);
						kept_lhs[l] = false;
						kept_rhs[r] = false;
						break;
					}
					else if (lhs_support && !rhs_support)
					{
						kept_rhs[r] = false;
					}
					else if (!lhs_support && rhs_support)
					{
						kept_lhs[l] = false;
						break;
					}					
				}
			}
		}		
    }
    
    // first trim off any polymerase run-ons, and make 3' ends consistent
    clip_by_3_prime_dropoff(lhs);
    clip_by_3_prime_dropoff(rhs);

    for (size_t i = 0; i < lhs.size(); ++i)
    {
        if (kept_lhs[i])
            scaffolds.push_back(lhs[i]);
    }
    
    for (size_t i = 0; i < rhs.size(); ++i)
    {
        if (kept_rhs[i])
            scaffolds.push_back(rhs[i]);
    }
}

void guess_strand(int bundle_origin, 
				  const vector<Scaffold>& hits,
				  vector<uint8_t>& strand_guess)
{
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].strand() == CUFF_STRAND_UNKNOWN)
			continue;
        
		for (int K = hits[i].left(); K < hits[i].right(); ++K)
			strand_guess[K - bundle_origin] |= hits[i].strand();
        
	}	
}

CuffStrand guess_strand_for_interval(const vector<uint8_t>& strand_guess, 
									 int left, 
									 int right)
{
	uint8_t guess = CUFF_STRAND_UNKNOWN;
	
	for (int i = left; i < right; ++i)
	{
		if (guess == CUFF_BOTH)
			return (CuffStrand)guess;
		guess |= strand_guess[i];
	}
	return (CuffStrand)guess;
}


bool scaffolds_for_bundle(const HitBundle& bundle, 
						  vector<shared_ptr<Scaffold> >& scaffolds,
						  vector<shared_ptr<Scaffold> >* ref_scaffs = NULL,
						  BundleStats* stats = NULL)
{
    if (bundle.hits().size() >= max_frags_per_bundle)
        return false;
    
	bool ref_guided = (ref_scaffs != NULL);
	
	vector<Scaffold> hits;
	vector<Scaffold> tmp_scaffs;
	
	for (size_t i = 0; i < bundle.hits().size(); ++i)
	{
		const MateHit& hit = bundle.hits()[i];
		hits.push_back(Scaffold(hit));
	}
    
  vector<float> depth_of_coverage(bundle.length(),0);
	vector<double> scaff_doc;
	map<pair<int,int>, float> intron_doc;
	
	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
	double bundle_avg_doc = compute_doc(bundle.left(), 
										hits, 
										depth_of_coverage, 
										intron_doc,
										true);
    
    if (bundle_avg_doc > 3000)
    {
        filter_introns(bundle.length(), 
                       bundle.left(), 
                       hits, 
                       min_isoform_fraction, 
                       false,
                       true);
    }
    
	if (ref_guided && enable_faux_reads && !hits.empty())
	{
		vector<Scaffold> pseudohits;
		foreach(shared_ptr<Scaffold const> ref_scaff, *ref_scaffs)
		{
			ref_scaff->tile_with_scaffs(pseudohits, tile_len, tile_off);
		}
		hits.insert(hits.end(),
					pseudohits.begin(),
					pseudohits.end());
		inplace_merge(hits.begin(),hits.end()-pseudohits.size(), hits.end(), scaff_lt);
	}
	
	vector<uint8_t> strand_guess(bundle.length(), CUFF_STRAND_UNKNOWN);
	guess_strand(bundle.left(),
				 hits,
				 strand_guess);
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].strand() == CUFF_STRAND_UNKNOWN)
		{
            assert (!hits[i].has_intron());
			uint8_t guess = CUFF_STRAND_UNKNOWN;
			Scaffold& hit = hits[i];
			const vector<AugmentedCuffOp>& ops = hit.augmented_ops();
            
			for (size_t j = 0; j < ops.size(); ++j)
			{
				const AugmentedCuffOp& op = ops[j];
				if (op.opcode == CUFF_UNKNOWN && op.genomic_length > (int)min_intron_length)
				{
					guess |= guess_strand_for_interval(strand_guess, 
													   hit.left() - bundle.left(),
													   hit.right() - bundle.left());
					
					break;
				}
			}

            
			if (guess != CUFF_BOTH && guess != CUFF_STRAND_UNKNOWN)
				hits[i].strand((CuffStrand)guess);
			//else
			//	fprintf(stderr, "Unknown strand for pair [%d-%d]\n", hit.left(), hit.right());
		}
	}
	
	bool saw_fwd = false;
	bool saw_rev = false;
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const Scaffold& hit = hits[i];
		CuffStrand hs = hit.strand();
						
		if (hs == CUFF_FWD)
			saw_fwd = true;
		if (hs == CUFF_REV)
			saw_rev = true;
		
//		if (hs != CUFF_REV) 
//			fwd_hits.push_back(hit);
//		if (hs != CUFF_FWD)
//			rev_hits.push_back(hit);
	}

    
	vector<Scaffold> fwd_scaffolds;
	vector<Scaffold> rev_scaffolds;
	
	bool assembled_successfully = false;
	
	if (saw_fwd && saw_rev)
	{
        // Forward strand hits
        {
            vector<Scaffold> fwd_hits;
            for (size_t i = 0; i < hits.size(); ++i)
            {
                const Scaffold& hit = hits[i];
                CuffStrand hs = hit.strand();
                if (hs != CUFF_REV) 
                    fwd_hits.push_back(hit);
            }
            
            verbose_msg ("%s\tFiltering forward strand\n", bundle_label->c_str());
            filter_hits(bundle.length(), bundle.left(), fwd_hits);
            assembled_successfully |= make_scaffolds(bundle.left(), 
                                                     bundle.length(), 
                                                     fwd_hits, 
                                                     fwd_scaffolds);
        }
        
        // Reverse strand hits
        {
            vector<Scaffold> rev_hits;
            for (size_t i = 0; i < hits.size(); ++i)
            {
                const Scaffold& hit = hits[i];
                CuffStrand hs = hit.strand();
                if (hs != CUFF_FWD)
                    rev_hits.push_back(hit);
            }
            
            verbose_msg ("%s\tFiltering reverse strand\n", bundle_label->c_str());
            filter_hits(bundle.length(), bundle.left(), rev_hits);
            assembled_successfully |= make_scaffolds(bundle.left(), 
                                                     bundle.length(), 
                                                     rev_hits, 
                                                     rev_scaffolds);
        }
	}
	else
	{
		if (saw_fwd || (!saw_fwd && !saw_rev))
		{
            // Forward strand hits
            {
                vector<Scaffold> fwd_hits;
                for (size_t i = 0; i < hits.size(); ++i)
                {
                    const Scaffold& hit = hits[i];
                    CuffStrand hs = hit.strand();
                    if (hs != CUFF_REV) 
                        fwd_hits.push_back(hit);
                }
                
                verbose_msg ("%s\tFiltering forward strand\n", bundle_label->c_str());
                filter_hits(bundle.length(), bundle.left(), fwd_hits);
                assembled_successfully |= make_scaffolds(bundle.left(), 
                                                         bundle.length(), 
                                                         fwd_hits, 
                                                         fwd_scaffolds);
            
            }
		}
		else
		{
            // Reverse strand hits
            {
                vector<Scaffold> rev_hits;
                for (size_t i = 0; i < hits.size(); ++i)
                {
                    const Scaffold& hit = hits[i];
                    CuffStrand hs = hit.strand();
                    if (hs != CUFF_FWD)
                        rev_hits.push_back(hit);
                }
                
                verbose_msg ("%s\tFiltering reverse strand\n", bundle_label->c_str());
                filter_hits(bundle.length(), bundle.left(), rev_hits);
                assembled_successfully |= make_scaffolds(bundle.left(), 
                                                         bundle.length(), 
                                                         rev_hits, 
                                                         rev_scaffolds);
            }
		}
	}
	
	combine_strand_assemblies(fwd_scaffolds, rev_scaffolds, tmp_scaffs, ref_scaffs);

	
	// Make sure all the reads are accounted for, including the redundant ones...
	for (size_t i = 0; i < tmp_scaffs.size(); ++i)
	{
		tmp_scaffs[i].clear_hits();
		for (size_t j = 0; j < bundle.hits().size(); ++j)
		{
			const MateHit& h = bundle.hits()[j];
			tmp_scaffs[i].add_hit(&h);
		}
	}
	
	if (ref_guided)
	{
		scaffolds = *ref_scaffs;
	}
	if (assembled_successfully)
	{
		foreach(Scaffold& scaff, tmp_scaffs)
		{
			scaffolds.push_back(shared_ptr<Scaffold>(new Scaffold(scaff)));
		}
	}
	sort(scaffolds.begin(), scaffolds.end(), scaff_lt_sp);
	
	return assembled_successfully;
}

//static long double min_abundance = 0.000001;

#if ENABLE_THREADS
boost::mutex out_file_lock;
boost::mutex thread_pool_lock;
int curr_threads = 0;

void decr_pool_count()
{
	thread_pool_lock.lock();
	curr_threads--;
	thread_pool_lock.unlock();	
}
#endif

void quantitate_transcript_cluster(AbundanceGroup& transfrag_cluster,
								   //const RefSequenceTable& rt,
                                   double total_map_mass,
                                   vector<Gene>& genes,
                                   bool bundle_too_large)
{
	if (transfrag_cluster.abundances().empty())
		return;
	
	vector<double> gammas;
    
	vector<MateHit> hits_in_cluster;
	
	get_alignments_from_scaffolds(transfrag_cluster.abundances(),
								  hits_in_cluster);
	
	
	// need the avg read length for depth of coverage calculation 
	double avg_read_length = 0;
	foreach (MateHit& hit, hits_in_cluster)
	{
		if (hit.left_alignment())
			avg_read_length += hit.left_alignment()->read_len(); 
		if (hit.right_alignment())
			avg_read_length += hit.right_alignment()->read_len(); 
	}
	
    if (hits_in_cluster.size())
        avg_read_length /= hits_in_cluster.size();
	
    if (library_type != "transfrags")
    {
        if (bundle_too_large == false)
        {
            transfrag_cluster.calculate_abundance(hits_in_cluster);
        }
        else
        {
            foreach(shared_ptr<Abundance>  ab, transfrag_cluster.abundances())
            {
                ab->status(NUMERIC_HI_DATA);
            }
        }
	}
    else
    {
        vector<shared_ptr<Abundance> >& abundances = transfrag_cluster.abundances();
        
        int N = abundances.size();
        double total_fpkm = 0.0;
        vector<double> gammas;
        for (size_t j = 0; j < N; ++j)
        {
            double FPKM = abundances[j]->transfrag()->fpkm();
            abundances[j]->FPKM(FPKM);
            total_fpkm += FPKM;
            gammas.push_back(FPKM);
        }
        
        for (size_t j = 0; j < N; ++j)
        {
            if (total_fpkm)
                gammas[j] /= total_fpkm;
        }
        
        vector<shared_ptr<Abundance> > filtered_transcripts = abundances;
        filter_junk_isoforms(filtered_transcripts, gammas, abundances, 0);
        vector<bool> to_keep (abundances.size(), false);
        for(size_t i = 0; i < abundances.size(); ++i)
        {
            shared_ptr<Abundance> ab_i = abundances[i];
            bool found = false;
            foreach (shared_ptr<Abundance> ab_j, filtered_transcripts)
            {
                if (ab_i == ab_j)
                {
                    found = true;
                    break;
                }
            }
            if (found)
                to_keep[i] = true;
        }
        
        AbundanceGroup kept;
        transfrag_cluster.filter_group(to_keep, kept);
        transfrag_cluster = kept;
    }
    
	vector<AbundanceGroup> transfrags_by_strand;
	cluster_transcripts<ConnectByStrand>(transfrag_cluster,
										 transfrags_by_strand);
	
	
	foreach (const AbundanceGroup& strand_group, transfrags_by_strand)
	{	
		vector<AbundanceGroup> transfrags_by_gene;
		
		if (bundle_mode == REF_DRIVEN)
		{
			cluster_transcripts<ConnectByAnnotatedGeneId>(strand_group, transfrags_by_gene);
		}
		else
		{
			cluster_transcripts<ConnectByExonOverlap>(strand_group, transfrags_by_gene);
		}

		foreach(const AbundanceGroup& gene, transfrags_by_gene)
		{
			const vector<shared_ptr<Abundance> >& iso_abundances = gene.abundances();
			vector<Isoform> isoforms;
			
			int gene_id = -1;
			int num_ref_gene_ids = 0;
            bool has_novel_isoform = false;
			string ref_gene_id = "";
			
			double major_isoform_FPKM = 0;
			foreach (shared_ptr<Abundance> iso_ab, iso_abundances)
			{
				if (iso_ab->transfrag()->is_ref())
				{
					if (iso_ab->transfrag()->annotated_gene_id() != ref_gene_id)
					{
						ref_gene_id = iso_ab->transfrag()->annotated_gene_id();
						num_ref_gene_ids++;
					}
				}
                else
                {
                    has_novel_isoform = true;
                }
				major_isoform_FPKM = max(iso_ab->FPKM(), major_isoform_FPKM);
			}
			
			foreach (shared_ptr<Abundance> iso_ab, iso_abundances)
			{
				// Calculate transcript depth of coverage and FMI from FPKM
				double FPKM = iso_ab->FPKM();
				double density_score = major_isoform_FPKM ? (FPKM / major_isoform_FPKM) : 0;
				double density_per_bp = FPKM;
				
				shared_ptr<Scaffold> transfrag = iso_ab->transfrag();
				assert(transfrag);
				
				double s_len = transfrag->length();
				
				density_per_bp *= (total_map_mass / 1000000.0); // yields (mass/(length/1000))
				density_per_bp *= (s_len/ 1000.0);
                double estimated_count = density_per_bp;
				density_per_bp /= s_len;
				density_per_bp *= avg_read_length;
				//double density_per_bp = (FPKM * (map_mass / 1000000.0) * 1000.0);
				
				if (!allow_junk_filtering || transfrag->is_ref() || density_score > min_isoform_fraction)
				{
					if (gene_id == -1 && (has_novel_isoform || num_ref_gene_ids > 1))
						gene_id = get_next_gene_id();
					
					isoforms.push_back(Isoform(*transfrag,
											   gene_id,
											   (int)isoforms.size() + 1,
											   FPKM,
											   iso_ab->effective_length(),
											   iso_ab->gamma(),
											   iso_ab->FPKM_conf(),
											   density_per_bp, 
                                               estimated_count,
											   density_score,
											   iso_ab->status(),
											   ref_gene_id));
				}
			}
			
			if (!isoforms.empty())
			{
				Gene g(isoforms, gene.FPKM(), gene.FPKM_conf(), gene.status());
				genes.push_back(g);	
			}
		}
		
	}
    
}

void quantitate_transcript_clusters(vector<shared_ptr<Scaffold> >& scaffolds,
									long double total_map_mass,
									vector<Gene>& genes,
                                    bool bundle_too_large)
{	
	//vector<shared_ptr<Scaffold> > partials;
	//vector<shared_ptr<Scaffold> > completes;
    
    vector<shared_ptr<Scaffold> > split_partials;
    // Cleave the partials at their unknowns to minimize FPKM dilation on  
    // the low end of the expression profile. 
    for (size_t i = 0; i < scaffolds.size(); ++i) 
    { 
        vector<Scaffold> c; 
        scaffolds[i]->get_complete_subscaffolds(c); 
        foreach (Scaffold& s, c)
        {
            split_partials.push_back(shared_ptr<Scaffold>(new Scaffold(s))); 
        }
    } 
    
    scaffolds = split_partials;
	
	vector<shared_ptr<Abundance> > abundances;
	foreach(shared_ptr<Scaffold> s, scaffolds)
	{
		TranscriptAbundance* pT = new TranscriptAbundance;
		pT->transfrag(s);
		shared_ptr<Abundance> ab(pT);
		abundances.push_back(ab);
	}
	
	AbundanceGroup transfrags = AbundanceGroup(abundances);
	
	vector<AbundanceGroup> transfrags_by_cluster;
	
	cluster_transcripts<ConnectByExonOverlap>(transfrags,
                                              transfrags_by_cluster);
	
	foreach(AbundanceGroup& cluster, transfrags_by_cluster)
	{
		quantitate_transcript_cluster(cluster, total_map_mass, genes, bundle_too_large);
	}
    verbose_msg( "%s\tBundle quantitation complete\n", bundle_label->c_str());
}

void assemble_bundle(const RefSequenceTable& rt,
					 HitBundle* bundle_ptr, 
					 BiasLearner* bl_ptr,
					 long double map_mass,
					 FILE* ftranscripts,
					 FILE* fgene_abundances,
					 FILE* ftrans_abundances,
					 FILE* fskipped)
{
	HitBundle& bundle = *bundle_ptr;
    
    char bundle_label_buf[2048];
    sprintf(bundle_label_buf, 
            "%s:%d-%d", 
            rt.get_name(bundle.ref_id()),
            bundle.left(),
            bundle.right());

#if ENABLE_THREADS
    bundle_label.reset(new string(bundle_label_buf));
#else
    bundle_label = shared_ptr<string>(new string(bundle_label_buf));
#endif

    verbose_msg( "%s\tProcessing new bundle with %d alignments\n", 
            bundle_label->c_str(),
            (int)bundle.hits().size());

#if ENABLE_THREADS	
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
	
	vector<shared_ptr<Scaffold> > scaffolds;
	
    bool successfully_assembled = true;
    
	switch(bundle_mode)
	{
		case REF_DRIVEN:
			scaffolds = bundle.ref_scaffolds();
			if (!final_est_run && scaffolds.size() != 1) // Only learn bias on single isoforms
			{
				delete bundle_ptr;
				return;
			}
			break;
		case REF_GUIDED:
			successfully_assembled = scaffolds_for_bundle(bundle, scaffolds, &bundle.ref_scaffolds());
			break;
		case HIT_DRIVEN:
			successfully_assembled = scaffolds_for_bundle(bundle, scaffolds);
			break;
		default:
			assert(false);
	}
	
    if (successfully_assembled == false)
    {

#if ENABLE_THREADS	
        out_file_lock.lock();
#endif
        
        int mask_region_id = get_next_skipped_region_id();
        fprintf(fskipped, 
                "%s\tCufflinks\texon\t%d\t%d\t%d\t%s\t.\tgene_id \"mask_%d\"; transcript_id \"mask_id%d/+\";\n",
                rt.get_name(bundle.ref_id()),
                bundle.left() + 1,
                bundle.right(), // GTF intervals are inclusive on both ends, but ours are half-open
                0,
                "+",
                mask_region_id,
                mask_region_id);
        
        fprintf(fskipped, 
                "%s\tCufflinks\texon\t%d\t%d\t%d\t%s\t.\tgene_id \"mask_%d\"; transcript_id \"mask_id%d/-\";\n",
                rt.get_name(bundle.ref_id()),
                bundle.left() + 1,
                bundle.right(), // GTF intervals are inclusive on both ends, but ours are half-open
                0,
                "-",
                mask_region_id,
                mask_region_id);
        
        
#if ENABLE_THREADS	
        out_file_lock.unlock();
#endif
        delete bundle_ptr;
		return;
    }
    
	if (scaffolds.empty())
	{
		delete bundle_ptr;
		return;
	}
		
	vector<Gene> genes;
    
    bool bundle_too_large = bundle_ptr->hits().size() >= max_frags_per_bundle;
    
    // FIXME: this routine does more than just quantitation, and should be 
    // renamed or refactored.
    quantitate_transcript_clusters(scaffolds, 
                                   map_mass,
                                   genes,
                                   bundle_too_large);
    
    verbose_msg( "%s\tFiltering bundle assembly\n", bundle_label->c_str());
    
    if (allow_junk_filtering)
        filter_junk_genes(genes);

	
	if (!final_est_run && bundle_mode==REF_DRIVEN) // Bias needs to be learned
	{
		for (size_t i = 0; i < genes.size(); ++i)
		{
			if (genes[i].isoforms().size() == 1)
			{
				bl_ptr -> preProcessTranscript(genes[i].isoforms()[0].scaffold()); 
			}
		}
	}
	
#if ENABLE_THREADS	
	out_file_lock.lock();
#endif
    
    // Get hit_introns for full_read_support test if ref-guided
    set<AugmentedCuffOp>* hit_introns = NULL;
    if (init_bundle_mode == REF_GUIDED)
    {
        hit_introns = new set<AugmentedCuffOp>();
        foreach(const MateHit& h, bundle.non_redundant_hits())
        {
            Scaffold s(h);
            foreach (AugmentedCuffOp a, s.augmented_ops())
            {
                if (a.opcode == CUFF_INTRON)
                {
                    hit_introns->insert(a);
                }
            }
        }
    }
    
	
	size_t num_scaffs_reported = 0;
	for (size_t i = 0; i < genes.size(); ++i)
	{
		const Gene& gene = genes[i];
		const vector<Isoform>& isoforms = gene.isoforms();
        set<string> annotated_gene_names;
        set<string> annotated_tss_ids;
		for (size_t j = 0; j < isoforms.size(); ++j)
		{
			const Isoform& iso = isoforms[j];
			
			vector<const MateHit*> H(iso.scaffold().mate_hits().size(), 0);
			copy(iso.scaffold().mate_hits().begin(), 
				 iso.scaffold().mate_hits().end(),
				 H.begin());
			
			vector<string> isoform_exon_recs;
            
			iso.get_gtf(isoform_exon_recs, rt, hit_introns);
			
			for (size_t g = 0; g < isoform_exon_recs.size(); ++g)
			{
				fprintf(ftranscripts, "%s", isoform_exon_recs[g].c_str());
			}
			
			fflush(ftranscripts);
			
			const char* status;
			if (iso.status()==NUMERIC_OK) 
				status = "OK";
			else if (iso.status() == NUMERIC_LOW_DATA)
                status = "LOWDATA";
            else if (iso.status() == NUMERIC_HI_DATA)
                status = "HIDATA";
            else if (iso.status() == NUMERIC_FAIL)
				status = "FAIL";
            else
                assert (false);
			
			fprintf(ftrans_abundances,"%s\t%c\t%s\t%s\t%s\t%s\t%s:%d-%d\t%d\t%lg\t%lg\t%lg\t%lg\t%s\n", 
					iso.trans_id().c_str(),
                    (iso.scaffold().nearest_ref_classcode() == 0 ? '-' : iso.scaffold().nearest_ref_classcode()),
                    (iso.scaffold().nearest_ref_id() == "" ? "-" : iso.scaffold().nearest_ref_id().c_str()),
                    gene.gene_id().c_str(),
                    (iso.scaffold().annotated_gene_name() == "" ? "-" : iso.scaffold().annotated_gene_name().c_str()), 
                    (iso.scaffold().annotated_tss_id() == "" ? "-" : iso.scaffold().annotated_tss_id().c_str()),
					rt.get_name(bundle.ref_id()),
					iso.scaffold().left(),
					iso.scaffold().right(),
                    iso.scaffold().length(),
                    iso.coverage(),
                    iso.FPKM(),
					iso.confidence().low,
					iso.confidence().high,
                    status);
			fflush(ftrans_abundances);
			
            annotated_gene_names.insert(iso.scaffold().annotated_gene_name());
            annotated_tss_ids.insert(iso.scaffold().annotated_tss_id());
            
			num_scaffs_reported++;
		}
		
		const char* status = "OK";
        if (gene.status()==NUMERIC_OK) 
            status = "OK";
        else if (gene.status() == NUMERIC_LOW_DATA)
            status = "LOWDATA";
        else if (gene.status() == NUMERIC_HI_DATA)
            status = "HIDATA";
        else if (gene.status() == NUMERIC_FAIL)
            status = "FAIL";
        else
            assert (false);

        string gene_names = cat_strings(annotated_gene_names);
        if (gene_names == "") gene_names = "-";
        string tss_ids = cat_strings(annotated_tss_ids);
        if (tss_ids == "") tss_ids = "-";
        
        fprintf(fgene_abundances,"%s\t%c\t%s\t%s\t%s\t%s\t%s:%d-%d\t%s\t%s\t%lg\t%lg\t%lg\t%s\n",
                gene.gene_id().c_str(),
                '-',
                "-",
                gene.gene_id().c_str(),
                gene_names.c_str(), 
                tss_ids.c_str(),
                rt.get_name(bundle.ref_id()),
                gene.left(),
                gene.right(),
                "-",
                "-",
                gene.FPKM(),
                gene.confidence().low,
                gene.confidence().high,
                status);
		fflush(fgene_abundances);
	}
    delete hit_introns;
	//fprintf(fbundle_tracking, "CLOSE %d\n", bundle.id());
	
	if (bundle_mode==REF_DRIVEN && num_scaffs_reported > bundle.ref_scaffolds().size())
    {
		fprintf(stderr, "Error: reported more isoforms than in reference!\n");
		exit(1);
	}
	
    verbose_msg( "%s\tBundle complete\n", bundle_label->c_str());
    
#if ENABLE_THREADS
	out_file_lock.unlock();
#endif

    genes.clear();
    scaffolds.clear();
	delete bundle_ptr;
}

bool assemble_hits(BundleFactory& bundle_factory, BiasLearner* bl_ptr)
{
	//srand(time(0));
		
	RefSequenceTable& rt = bundle_factory.ref_table();
    
	//FILE* fbundle_tracking = fopen("open_bundles", "w");
    
	//FILE* fstats = fopen("bundles.stats", "w");
	FILE* ftrans_abundances = fopen(string(output_dir + "/" + "isoforms.fpkm_tracking").c_str(), "w");
	//fprintf(ftrans_abundances,"trans_id\tbundle_id\tchr\tleft\tright\tFPKM\tFMI\tfrac\tFPKM_conf_lo\tFPKM_conf_hi\tcoverage\tlength\teffective_length\tstatus\n");
	fprintf(ftrans_abundances,"tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n");
	FILE* fgene_abundances = fopen(string(output_dir + "/" + "genes.fpkm_tracking").c_str(), "w");
	//fprintf(fgene_abundances,"gene_id\tbundle_id\tchr\tleft\tright\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tstatus\n");
    fprintf(fgene_abundances,"tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n");
    
	FILE* ftranscripts = fopen(string(output_dir + "/" + "transcripts.gtf").c_str(), "w");
    FILE* fskipped = fopen(string(output_dir + "/" + "skipped.gtf").c_str(), "w");
    
	string process;
	if (corr_bias && corr_multi && final_est_run)
		process = "Re-estimating abundances with bias and multi-read correction.";
	else if (corr_multi && final_est_run)
		process = "Re-estimating abundances with multi-read correction.";
	else if (corr_bias && final_est_run)
		process = "Re-estimating abundances with bias correction.";
	else if (bundle_mode==REF_DRIVEN && final_est_run)
		process = "Estimating transcript abundances."; 
	else if (bundle_mode==REF_DRIVEN && corr_bias)
		process = "Learning bias parameters.";
	else if (bundle_mode==REF_DRIVEN && corr_multi)
		process = "Initializing transcript abundances for multi-read correction.";
	else if (corr_multi)
		process = "Assembling transcripts and initializing abundances for multi-read correction.";
	else
		process = "Assembling transcripts and estimating abundances.";
		
	ProgressBar p_bar(process, bundle_factory.read_group_properties()->total_map_mass());

	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		if (!bundle_factory.next_bundle(*bundle_ptr))
		{
			delete bundle_ptr;
			break;
		}
		
		HitBundle& bundle = *bundle_ptr;
        
		char bundle_label_buf[2048];
		sprintf(bundle_label_buf, 
				"%s:%d-%d", 
				rt.get_name(bundle.ref_id()),
				bundle.left(),
				bundle.right());

		if (bundle.right() - bundle.left() > max_gene_length)
		{
			fprintf(stderr, "\n%s\tWarning: Skipping large bundle.\n", bundle_label_buf);
			delete bundle_ptr;
			continue;
		}

		BundleStats stats;
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
		p_bar.update(bundle_label_buf, bundle.raw_mass());	

#if ENABLE_THREADS			
		thread_pool_lock.lock();
		curr_threads++;
		thread_pool_lock.unlock();
		
		thread asmbl(assemble_bundle,
					 boost::cref(rt), 
					 bundle_ptr, 
					 bl_ptr,
					 bundle_factory.read_group_properties()->normalized_map_mass(),
					 ftranscripts, 
					 fgene_abundances,
					 ftrans_abundances,
                     fskipped);
#else
		assemble_bundle(boost::cref(rt), 
						bundle_ptr, 
						bl_ptr,
						bundle_factory.read_group_properties()->normalized_map_mass(),
						ftranscripts,
						fgene_abundances,
						ftrans_abundances,
                        fskipped);
#endif			
		
	}

#if ENABLE_THREADS	
	while(1)
	{
		thread_pool_lock.lock();
		if (curr_threads == 0)
		{
			thread_pool_lock.unlock();
			break;
		}
		p_bar.remaining(curr_threads);
		
		thread_pool_lock.unlock();
		//fprintf(stderr, "waiting to exit\n");
		boost::this_thread::sleep(boost::posix_time::milliseconds(5));
	}
#endif
	
	p_bar.complete();
	
	if(!final_est_run && bundle_mode==REF_DRIVEN) // We are learning bias
	{
		bl_ptr->normalizeParameters();
        if (output_bias_params)
            bl_ptr->output();
	}
	
	fclose(ftranscripts);
	fclose(ftrans_abundances);
	fclose(fgene_abundances);
    fclose(fskipped);
	return true;
}
	
void driver(const string& hit_file_name, FILE* ref_gtf, FILE* mask_gtf)
{	
    ReadTable it;
	RefSequenceTable rt(true, false);
	    
	shared_ptr<HitFactory> hit_factory;

    try
	{
		hit_factory = shared_ptr<BAMHitFactory>(new BAMHitFactory(hit_file_name, it, rt));
	}
	catch (std::runtime_error& e)
	{
		fprintf(stderr, "File %s doesn't appear to be a valid BAM file, trying SAM...\n",
				hit_file_name.c_str());
	
        try
        {
            hit_factory = shared_ptr<SAMHitFactory>(new SAMHitFactory(hit_file_name, it, rt));
        }
        catch (std::runtime_error& e)
        {
            fprintf(stderr, "Error: cannot open alignment file %s for reading\n",
                    hit_file_name.c_str());
            exit(1);
        }
	}
	
	BundleFactory& bundle_factory = *(new BundleFactory(hit_factory, bundle_mode));
	shared_ptr<ReadGroupProperties> rg_props =bundle_factory.read_group_properties();
	BadIntronTable bad_introns;
    
    rt.print_rec_ordering();
    
    vector<shared_ptr<Scaffold> > ref_mRNAs;
    if (ref_gtf)
    {
        ::load_ref_rnas(ref_gtf, bundle_factory.ref_table(), ref_mRNAs, corr_bias && bundle_mode == REF_DRIVEN, false);
        bundle_factory.set_ref_rnas(ref_mRNAs);
    }
    rt.print_rec_ordering();
    vector<shared_ptr<Scaffold> > mask_rnas;
    if (mask_gtf)
    {
        ::load_ref_rnas(mask_gtf, bundle_factory.ref_table(), mask_rnas, false, false);
        bundle_factory.set_mask_rnas(mask_rnas);
    }
    
    vector<LocusCount> count_table;
    if (bundle_mode != HIT_DRIVEN)
        inspect_map(bundle_factory, NULL, count_table);
    else 
        inspect_map(bundle_factory, &bad_introns, count_table);
    
    
    verbose_msg("%d ReadHits still live\n", num_deleted);
    verbose_msg("Found %lu reference contigs\n", rt.size());
    
    foreach(shared_ptr<Scaffold> ref_scaff, ref_mRNAs)
    {
        ref_scaff->clear_hits();
    }
    
    //fprintf(stderr, "ReadHit delete count is %d\n", num_deleted);
    
	BiasLearner* bl_ptr = new BiasLearner(rg_props->frag_len_dist());
    bundle_factory.read_group_properties(rg_props);

	//if (ref_gtf) -- why? bad introns are bad
		bundle_factory.bad_intron_table(bad_introns);
	
	max_frag_len = rg_props->frag_len_dist()->max();
	min_frag_len = rg_props->frag_len_dist()->min();
	verbose_msg("\tTotal map density: %Lf\n", rg_props->total_map_mass());

	if (corr_bias || corr_multi) final_est_run = false;

	assemble_hits(bundle_factory, bl_ptr);
    
	if (final_est_run) 
    {
        ref_mRNAs.clear();
        return;
    }
    
	hit_factory->reset();
	delete &bundle_factory;
	BundleFactory bundle_factory2(hit_factory, REF_DRIVEN);
	rg_props->bias_learner(shared_ptr<BiasLearner const>(bl_ptr));
	rg_props->multi_read_table()->valid_mass(true);
	bundle_factory2.read_group_properties(rg_props);

    if (bundle_mode==HIT_DRIVEN || bundle_mode==REF_GUIDED)
    {
		ref_gtf = fopen(string(output_dir + "/transcripts.gtf").c_str(), "r");
        ref_mRNAs.clear();
        ::load_ref_rnas(ref_gtf, bundle_factory2.ref_table(), ref_mRNAs, corr_bias, true);
    }    
	bundle_factory2.set_ref_rnas(ref_mRNAs);
    if (mask_gtf)
    {
        mask_rnas.clear();
        ::load_ref_rnas(mask_gtf, bundle_factory2.ref_table(), mask_rnas, false, false);
        bundle_factory2.set_mask_rnas(mask_rnas);
    }    
	bundle_factory2.reset();
	
	if(corr_bias && (bundle_mode==HIT_DRIVEN || bundle_mode==REF_GUIDED)) 
	{
        // We still need to learn the bias since we didn't have the sequences before assembly
		learn_bias(bundle_factory2, *bl_ptr);
		bundle_factory2.reset();
	}

    bundle_mode = REF_DRIVEN;
	final_est_run = true;
	assemble_hits(bundle_factory2, bl_ptr);
	ref_mRNAs.clear();
}

int main(int argc, char** argv)
{	
    init_library_table();
  string cmdline;
  for (int i=0;i<argc;i++) {
    cmdline+=argv[i];
    cmdline+=" ";
    }
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
    if (!use_total_mass && !use_compat_mass)
    {
        use_total_mass = true;
        use_compat_mass = false;   
    }
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    if (!no_update_check)
        check_version(PACKAGE_VERSION);
    
    if (cuff_quiet || cuff_verbose)
      fprintf(stderr, "Command line:\n%s\n", cmdline.c_str());
    string sam_hits_file_name = argv[optind++];
	

	if (random_seed == -1)
        random_seed = time(NULL);
    
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
    
    driver(sam_hits_file_name, ref_gtf, mask_gtf);
	
	return 0;
}
