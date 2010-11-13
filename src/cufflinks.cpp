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

#include "clustering.h"
#include "abundances.h"
#include "bundles.h"
#include "filters.h"
#include "genes.h"
#include "assemble.h"
#include "biascorrection.h"
#include "gtf_tracking.h"

using namespace std;

#if ENABLE_THREADS
const char *short_options = "m:p:s:F:c:I:j:Q:L:G:f:o:M:r:a:A:Nqv";
#else
const char *short_options = "m:s:F:c:I:j:Q:L:G:f:o:M:r:a:A:Nqv";
#endif

static struct option long_options[] = {
{"inner-dist-mean",			required_argument,       0,          'm'},
{"inner-dist-stddev",		required_argument,       0,          's'},
{"transcript-score-thresh", required_argument,       0,          't'},
{"min-isoform-fraction",    required_argument,       0,          'F'},
{"min-intron-fraction",     required_argument,       0,          'f'},
{"pre-mrna-fraction",		required_argument,		 0,			 'j'},
{"junc-alpha",				required_argument,		 0,			 'a'},	
{"small-anchor-fraction",	required_argument,		 0,			 'A'},
{"max-intron-length",		required_argument,		 0,			 'I'},
{"min-map-qual",			required_argument,		 0,			 'Q'},
{"label",					required_argument,		 0,			 'L'},
{"collapse-thresh",         required_argument,		 0,			 'c'},
{"GTF",					    required_argument,		 0,			 'G'},
{"mask-gtf",                required_argument,		 0,			 'M'},
{"output-dir",			    required_argument,		 0,			 'o'},
{"quartile-normalization",  no_argument,	 		 0,	         'N'},
{"verbose",			    	no_argument,		 	0,			 'v'},
{"quiet",			    	no_argument,		 	0,			 'q'},
{"reference-seq",			required_argument,		 0,			 'r'},	
#if ENABLE_THREADS
{"num-threads",				required_argument,       0,          'p'},
#endif
{"overhang-tolerance",      required_argument,		 0,			 OPT_OVERHANG_TOLERANCE},

{"num-importance-samples",  required_argument,		 0,			 OPT_NUM_IMP_SAMPLES},
{"max-mle-iterations",		required_argument,		 0,			 OPT_MLE_MAX_ITER},
{"library-type",		    required_argument,		 0,			 OPT_LIBRARY_TYPE},
{"max-bundle-length",       required_argument,		 0,			 OPT_MAX_BUNDLE_LENGTH},
{"min-frags-per-transfrags",required_argument,		 0,			 OPT_MIN_FRAGS_PER_TRANSFRAG},
{"min-intron-length",required_argument,		         0,			 OPT_MIN_INTRON_LENGTH},

#if ADAM_MODE
{"bias-mode",		    required_argument,		 0,			 OPT_BIAS_MODE},
#endif

{0, 0, 0, 0} // terminator
};

void print_usage()
{
    //NOTE: SPACES ONLY, bozo (I agree -- everywhere else in the code too, please don't use tabs)
    fprintf(stderr, "cufflinks v%s\n", PACKAGE_VERSION);
    fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
    fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "Usage:   cufflinks [options] <hits.sam>\n");
    fprintf(stderr, "Options:\n\n");
#if ENABLE_THREADS
    fprintf(stderr, "  -p/--num-threads             number of threads used during analysis                [ default:      1 ]\n");
#endif
    
    fprintf(stderr, "  -L/--label                   assembled transcripts have this ID prefix             [ default:   CUFF ]\n");
    fprintf(stderr, "  -G/--GTF                     quantitate against reference transcript annotations                      \n");
    fprintf(stderr, "  -F/--min-isoform-fraction    suppress transcripts below this abundance level       [ default:   0.15 ]\n");
    fprintf(stderr, "  -f/--min-intron-fraction     filter spliced alignments below this level            [ default:   0.05 ]\n");
    fprintf(stderr, "  -j/--pre-mrna-fraction       suppress intra-intronic transcripts below this level  [ default:   0.15 ]\n");
    fprintf(stderr, "  -I/--max-intron-length       ignore alignments with gaps longer than this          [ default: 300000 ]\n");
    fprintf(stderr, "  -Q/--min-map-qual            ignore alignments with lower than this mapping qual   [ default:      0 ]\n");
    fprintf(stderr, "  -M/--mask-file               ignore all alignment within transcripts in this file                     \n");
    fprintf(stderr, "  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
	fprintf(stderr, "  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
    fprintf(stderr, "  -o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");
    fprintf(stderr, "  -r/--reference-seq           reference fasta file for sequence bias correction     [ default:   NULL ]\n");
    fprintf(stderr, "\nAdvanced Options:\n\n");
    fprintf(stderr, "  -N/--quartile-normalization  use quartile normalization instead of total counts    [ default:  FALSE ]\n");
    fprintf(stderr, "  -a/--junc-alpha              alpha for junction binomial test filter               [ default:   0.01 ]\n");
    fprintf(stderr, "  -A/--small-anchor-fraction   percent read overhang taken as 'suspiciously small'   [ default:   0.12 ]\n");
    fprintf(stderr, "  -m/--frag-len-mean           the average fragment length                           [ default:    200 ]\n");
    fprintf(stderr, "  -s/--frag-len-std-dev        the fragment length standard deviation                [ default:     80 ]\n");
    fprintf(stderr, "  --min-frags-per-transfrag    minimum number of fragments needed for new transfrags [ default:     10 ]\n");
    fprintf(stderr, "  --overhang-tolerance         number of terminal exon bp to tolerate in introns     [ default:      8 ]\n");
    fprintf(stderr, "  --num-importance-samples     number of importance samples for MAP restimation      [ default:   1000 ]\n");
    fprintf(stderr, "  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]\n");
    fprintf(stderr, "  --library-type               Library prep used for input reads                     [ default:  below ]\n");
    fprintf(stderr, "  --max-bundle-length          maximum genomic length allowed for a given bundle     [ default:3500000 ]\n");
    fprintf(stderr, "  --min-intron-length          minimum intron size allowed in genome                 [ default:     50 ]\n");
    
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
				def_frag_len_mean = (uint32_t)parseInt(0, "-m/--frag-len-mean arg must be at least 0", print_usage);
				break;
			case 's':
				def_frag_len_std_dev = (uint32_t)parseInt(0, "-s/--frag-len-std-dev arg must be at least 0", print_usage);
				break;
			case 't':
				transcript_score_thresh = parseFloat(-99999999, 0, "-t/--transcript-score-thresh must less than or equal to 0", print_usage);
				break;
			case 'p':
				num_threads = (uint32_t)parseInt(1, "-p/--num-threads arg must be at least 1", print_usage);
				break;
			case 'F':
				min_isoform_fraction = parseFloat(0, 1.0, "-F/--min-isoform-fraction must be between 0 and 1.0", print_usage);
				F_set = true;
				break;
			case 'f':
				min_isoform_fraction = parseFloat(0, 1.0, "-f/--min-intron-fraction must be between 0 and 1.0", print_usage);
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
				bias_mode = optarg;
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
			case 'L':
			{
				user_label = optarg;
				break;
			}
			case 'G':
			{
				ref_gtf_filename = optarg;
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
            case 'r':
			{
				fasta_dir = optarg;
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
                
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    
	if (ref_gtf_filename != "")
	{
        if (!F_set)
        {
            min_isoform_fraction = 0.0;
        }
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
            if (library_type == "transfrags")
            {
                allow_junk_filtering = false;
            }
            global_read_properties = &lib_itr->second;
        }
    }
	
#if ADAM_MODE
	if (fasta_dir != "")
		output_dir = output_dir + "/" + bias_mode;
#endif
	
    
    //inner_dist_norm = normal(0, inner_dist_std_dev);
	return 0;
}

void add_non_shadow_scaffs(const vector<Scaffold>& lhs, 
						   const vector<Scaffold>& rhs,
						   vector<Scaffold>& scaffolds,
						   bool include_unknown_strand)
{
	for (size_t i = 0; i < lhs.size(); ++i)
	{
		bool add_to_asm = true;
		if (lhs[i].strand() == CUFF_STRAND_UNKNOWN)
		{
			for (size_t j = 0; j < rhs.size(); ++j)
			{
				if (include_unknown_strand || 
					rhs[j].strand() != CUFF_STRAND_UNKNOWN)
				{
					if (Scaffold::overlap_in_genome(lhs[i], rhs[j], 0) &&
						Scaffold::compatible(lhs[i], rhs[j]))
					{
						add_to_asm = false;
						break;
					}
				}
			}
		}
		if (add_to_asm)
		{
			scaffolds.push_back(lhs[i]);	
		}
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
						  BundleStats* stats = NULL)
{
	vector<Scaffold> hits;
	vector<Scaffold> tmp_scaffs;
	
    //	for (size_t i = 0; i < bundle.non_redundant_hits().size(); ++i)
    //	{
    //		const MateHit& hit = bundle.non_redundant_hits()[i];
    //		hits.push_back(Scaffold(hit));
    //	}
    
	for (size_t i = 0; i < bundle.hits().size(); ++i)
	{
		const MateHit& hit = bundle.hits()[i];
		hits.push_back(Scaffold(hit));
	}
    
    vector<int> depth_of_coverage(bundle.length(),0);
	vector<double> scaff_doc;
	map<pair<int,int>, int> intron_doc;
	
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
    
	vector<uint8_t> strand_guess(bundle.length(), CUFF_STRAND_UNKNOWN);
	guess_strand(bundle.left(),
				 hits,
				 strand_guess);
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].strand() == CUFF_STRAND_UNKNOWN)
		{
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
	
	vector<Scaffold> fwd_hits, rev_hits;
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
		
		if (hs != CUFF_REV) 
			fwd_hits.push_back(hit);
		if (hs != CUFF_FWD)
			rev_hits.push_back(hit);
	}
	
	{
		asm_verbose ("%s\tFiltering forward strand\n", bundle_label->c_str());
		filter_hits(bundle.length(), bundle.left(), fwd_hits);
		asm_verbose ("%s\tFiltering reverse strand\n", bundle_label->c_str());
		filter_hits(bundle.length(), bundle.left(), rev_hits);
	}
    
	vector<Scaffold> fwd_scaffolds;
	vector<Scaffold> rev_scaffolds;
	
	bool assembled_successfully = false;
	
	if (saw_fwd && saw_rev)
	{
		assembled_successfully |= make_scaffolds(bundle.left(), 
                                                 bundle.length(), 
                                                 fwd_hits, 
                                                 fwd_scaffolds);
        
		assembled_successfully |= make_scaffolds(bundle.left(), 
                                                 bundle.length(), 
                                                 rev_hits, 
                                                 rev_scaffolds);
		
		add_non_shadow_scaffs(fwd_scaffolds, rev_scaffolds, tmp_scaffs, true);
		add_non_shadow_scaffs(rev_scaffolds, fwd_scaffolds, tmp_scaffs, false);
	}
	else
	{
		if (saw_fwd || (!saw_fwd && !saw_rev))
		{
			assembled_successfully |= make_scaffolds(bundle.left(),
													 bundle.length(),
													 fwd_hits,
													 fwd_scaffolds);
			tmp_scaffs.insert(tmp_scaffs.end(),fwd_scaffolds.begin(), fwd_scaffolds.end());
		}
		else
		{
			assembled_successfully |= make_scaffolds(bundle.left(), 
													 bundle.length(), 
													 rev_hits, 
													 rev_scaffolds);
			tmp_scaffs.insert(tmp_scaffs.end(),rev_scaffolds.begin(), rev_scaffolds.end());
		}
	}
	
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
	
	sort(tmp_scaffs.begin(), tmp_scaffs.end(), scaff_lt);
	
	foreach(Scaffold& scaff, tmp_scaffs)
	{
		scaffolds.push_back(shared_ptr<Scaffold>(new Scaffold(scaff)));
	}
	
	return assembled_successfully;
}


//static long double min_abundance = 0.000001;

#if ENABLE_THREADS
mutex out_file_lock;
mutex thread_pool_lock;
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
                                   vector<Gene>& genes)
{
	if (transfrag_cluster.abundances().empty())
		return;
	
	vector<double> gammas;
	ublas::matrix<double> gamma_covariance;
    
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
        transfrag_cluster.calculate_abundance(hits_in_cluster);
	}
    
	vector<AbundanceGroup> transfrags_by_strand;
	cluster_transcripts<ConnectByStrand>(transfrag_cluster,
										 transfrags_by_strand);
	
	foreach (const AbundanceGroup& strand_group, transfrags_by_strand)
	{	
		vector<AbundanceGroup> transfrags_by_gene;
		
		cluster_transcripts<ConnectByAnnotatedGeneId>(strand_group,
													  transfrags_by_gene);
		
		foreach (const AbundanceGroup& gene, transfrags_by_gene)
		{
			const vector<shared_ptr<Abundance> >& iso_abundances = gene.abundances();
			vector<Isoform> isoforms;
			
            int gene_id = get_next_gene_id();
            
			double major_isoform_FPKM = 0;
			foreach (shared_ptr<Abundance> iso_ab, iso_abundances)
			{
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
				density_per_bp /= s_len;
				density_per_bp *= avg_read_length;
				//double density_per_bp = (FPKM * (map_mass / 1000000.0) * 1000.0);
				
				if (!allow_junk_filtering || density_score > min_isoform_fraction || major_isoform_FPKM == 0.0)
				{
					isoforms.push_back(Isoform(*transfrag,
											   gene_id,
											   (int)isoforms.size() + 1,
											   FPKM,
											   iso_ab->effective_length(),
											   iso_ab->gamma(),
											   iso_ab->FPKM_conf(),
											   density_per_bp, 
											   density_score,
											   iso_ab->status()));
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
									vector<Gene>& genes)
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
		quantitate_transcript_cluster(cluster, total_map_mass, genes);
	}
    asm_verbose( "%s\tBundle quantitation complete\n", bundle_label->c_str());
}

void assemble_bundle(const RefSequenceTable& rt,
					 HitBundle* bundle_ptr, 
					 long double map_mass,
					 FILE* ftranscripts,
					 FILE* fgene_abundances,
					 FILE* ftrans_abundances)
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
	
	if (ref_gtf_filename != "")
	{
		scaffolds = bundle.ref_scaffolds();
	}
	else 
	{
		bool success = scaffolds_for_bundle(bundle, scaffolds, NULL);
		if (!success)
		{
			delete bundle_ptr;
			return;
		}
	}
	
	vector<Gene> genes;
    
    // FIXME: this routine does more than just quantitation, and should be 
    // renamed or refactored.
    quantitate_transcript_clusters(scaffolds, 
                                   map_mass,
                                   genes);
    
    asm_verbose( "%s\tFiltering bundle assembly\n", bundle_label->c_str());
    
    if (allow_junk_filtering)
        filter_junk_genes(genes);

    
#if ENABLE_THREADS	
	out_file_lock.lock();
#endif
	
	size_t num_scaffs_reported = 0;
	for (size_t i = 0; i < genes.size(); ++i)
	{
		const Gene& gene = genes[i];
		const vector<Isoform>& isoforms = genes[i].isoforms();
		for (size_t j = 0; j < isoforms.size(); ++j)
		{
			const Isoform& iso = isoforms[j];
			
			vector<const MateHit*> H(iso.scaffold().mate_hits().size(), 0);
			copy(iso.scaffold().mate_hits().begin(), 
				 iso.scaffold().mate_hits().end(),
				 H.begin());
			
			vector<string> isoform_exon_recs;
            
			iso.get_gtf(isoform_exon_recs, rt);
			
			for (size_t g = 0; g < isoform_exon_recs.size(); ++g)
			{
				fprintf(ftranscripts, "%s", isoform_exon_recs[g].c_str());
			}
			
			fflush(ftranscripts);
			
			const char* status;
			if (iso.status()==NUMERIC_OK) 
				status = "OK";
			else 
				status = "FAIL";
			
			fprintf(ftrans_abundances,"%s\t%d\t%s\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\t%lg\t%s\n", 
					iso.trans_id().c_str(),
					bundle.id(),
					rt.get_name(bundle.ref_id()),
					iso.scaffold().left(),
					iso.scaffold().right(),
					iso.FPKM(),
					iso.FMI(),
					iso.fraction(),
					iso.confidence().low,
					iso.confidence().high,
					iso.coverage(),
					iso.scaffold().length(),
					iso.effective_length(),
					status);
			fflush(ftrans_abundances);
			
			num_scaffs_reported++;
		}
		
		const char* status;
		if (gene.status()==NUMERIC_OK)
			status = "OK";
		else
			status = "FAIL";

		fprintf(fgene_abundances,"%s\t%d\t%s\t%d\t%d\t%lg\t%lg\t%lg\t%s\n",
				gene.gene_id().c_str(),
				bundle.id(),
				rt.get_name(bundle.ref_id()),
				gene.left(),
				gene.right(),
				gene.FPKM(),
				gene.confidence().low,
				gene.confidence().high,
				status);
		fflush(fgene_abundances);
	}
    
	//fprintf(fbundle_tracking, "CLOSE %d\n", bundle.id());
	
	if (ref_gtf_filename != "" && num_scaffs_reported > bundle.ref_scaffolds().size())
    {
		fprintf(stderr, "Error: reported more isoforms than in reference!\n");
		exit(1);
	}
	
    asm_verbose( "%s\tBundle complete\n", bundle_label->c_str());
    
#if ENABLE_THREADS
	out_file_lock.unlock();
#endif

    genes.clear();
    scaffolds.clear();
	delete bundle_ptr;
}

bool assemble_hits(BundleFactory& bundle_factory)
{
	srand(time(0));
	
	int num_bundles = 0;
	
	RefSequenceTable& rt = bundle_factory.ref_table();
    
	//FILE* fbundle_tracking = fopen("open_bundles", "w");
    
	//FILE* fstats = fopen("bundles.stats", "w");
	FILE* ftrans_abundances = fopen(string(output_dir + "/" + "transcripts.expr").c_str(), "w");
	fprintf(ftrans_abundances,"trans_id\tbundle_id\tchr\tleft\tright\tFPKM\tFMI\tfrac\tFPKM_conf_lo\tFPKM_conf_hi\tcoverage\tlength\teffective_length\tstatus\n");
	
	FILE* fgene_abundances = fopen(string(output_dir + "/" + "genes.expr").c_str(), "w");
	fprintf(fgene_abundances,"gene_id\tbundle_id\tchr\tleft\tright\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tstatus\n");
	FILE* ftranscripts = fopen(string(output_dir + "/" + "transcripts.gtf").c_str(), "w");
    
	string process;
	if (ref_gtf_filename != "" && fasta_dir != "" && final_est_run)
		process = "Re-estimating abundances with bias correction.";
	else if (ref_gtf_filename != "")
		process = "Calculating estimated abundances.";
	else
		process = "Assembling transcripts and estimating abundances.";
	ProgressBar p_bar(process, bundle_factory.num_bundles());

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

		if (bundle.right() - bundle.left() > 3000000)
		{
			verbose_msg( "%s\tWarning: large bundle encountered\n", bundle_label_buf);
		}

		BundleStats stats;
		num_bundles++;
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
		p_bar.update(bundle_label_buf, 1);	

#if ENABLE_THREADS			
		thread_pool_lock.lock();
		curr_threads++;
		thread_pool_lock.unlock();
		
		thread asmbl(assemble_bundle,
					 boost::cref(rt), 
					 bundle_ptr, 
					 bundle_factory.read_group_properties()->total_map_mass(),
					 ftranscripts, 
					 fgene_abundances,
					 ftrans_abundances);
#else
		assemble_bundle(boost::cref(rt), 
						bundle_ptr, 
						bundle_factory.read_group_properties()->total_map_mass(),
						ftranscripts,
						fgene_abundances,
						ftrans_abundances);
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
		
		thread_pool_lock.unlock();
		//fprintf(stderr, "waiting to exit\n");
		boost::this_thread::sleep(boost::posix_time::milliseconds(5));
	}
#endif
	
	p_bar.complete();
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
	BundleFactory& bundle_factory = *(new BundleFactory(hit_factory));
	
	shared_ptr<EmpDist> frag_len_dist(new EmpDist);
	long double map_mass = 0.0;
	BadIntronTable bad_introns;
    
    rt.print_rec_ordering();
    
    vector<shared_ptr<Scaffold> > ref_mRNAs;
    if (ref_gtf)
    {
        ::load_ref_rnas(ref_gtf, bundle_factory.ref_table(), ref_mRNAs, false, false);
        bundle_factory.set_ref_rnas(ref_mRNAs);
    }
    rt.print_rec_ordering();
    vector<shared_ptr<Scaffold> > mask_rnas;
    if (mask_gtf)
    {
        ::load_ref_rnas(mask_gtf, bundle_factory.ref_table(), mask_rnas, false, false);
        bundle_factory.set_mask_rnas(mask_rnas);
    }
    
    if (ref_gtf)
    {
        inspect_map(bundle_factory, map_mass, NULL, *frag_len_dist);
    }
    else 
    {
        inspect_map(bundle_factory, map_mass, &bad_introns, *frag_len_dist);
    }
    
    asm_verbose("%d ReadHits still live\n", num_deleted);
    asm_verbose("Found %lu reference contigs\n", rt.size());
    
    foreach(shared_ptr<Scaffold> ref_scaff, ref_mRNAs)
    {
        ref_scaff->clear_hits();
    }
    
    //fprintf(stderr, "ReadHit delete count is %d\n", num_deleted);
    shared_ptr<ReadGroupProperties> rg_props(new ReadGroupProperties);
    
    if (global_read_properties)
    {
        *rg_props = *global_read_properties;
    }
    else 
    {
        *rg_props = hit_factory->read_group_properties();
    }
    
    rg_props->frag_len_dist(frag_len_dist);
    rg_props->total_map_mass(map_mass);
    
    bundle_factory.read_group_properties(rg_props);

	if (ref_gtf_filename == "")
	{
		bundle_factory.bad_intron_table(bad_introns);
	}
	
	max_frag_len = frag_len_dist->max();
	min_frag_len = frag_len_dist->min();
	verbose_msg("\tTotal map density: %Lf\n", map_mass);

	if (fasta_dir != "") final_est_run = false;
#if ADAM_MODE
	if (fasta_dir == "")  assemble_hits(bundle_factory);
#else
	assemble_hits(bundle_factory);
#endif
    if (fasta_dir == "") 
    {
        ref_mRNAs.clear();
        return;
    }
    
	hit_factory->reset();
	int num_bundles = bundle_factory.num_bundles();
	delete &bundle_factory;
	BundleFactory bundle_factory2(hit_factory);
	bundle_factory2.num_bundles(num_bundles);

#if ADAM_MODE
	ref_gtf = fopen(string(output_dir + "/../init/transcripts.gtf").c_str(), "r");
#else	
	ref_gtf = fopen(string(output_dir + "/transcripts.gtf").c_str(), "r");
#endif

    if (ref_gtf)
    {
        ref_mRNAs.clear();
        ::load_ref_rnas(ref_gtf, bundle_factory2.ref_table(), ref_mRNAs, true, true);
        bundle_factory2.set_ref_rnas(ref_mRNAs);
    }
    
    if (mask_gtf)
    {
        mask_rnas.clear();
        ::load_ref_rnas(mask_gtf, bundle_factory2.ref_table(), mask_rnas, false, false);
        bundle_factory2.set_mask_rnas(mask_rnas);
    }
    
    bundle_factory2.read_group_properties(rg_props);
    
	BiasLearner* bl = new BiasLearner(rg_props->frag_len_dist());
	learn_bias(bundle_factory2, *bl);
	rg_props->bias_learner(shared_ptr<BiasLearner const>(bl));
	
	bundle_factory2.reset();

	final_est_run = true;
	assemble_hits(bundle_factory2);
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
	
    string sam_hits_file_name = argv[optind++];
	

	
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
