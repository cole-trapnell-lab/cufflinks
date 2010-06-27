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
#include "gtf_tracking.h"

using namespace std;

#if ENABLE_THREADS
const char *short_options = "m:p:s:F:c:I:j:Q:L:G:f:o:";
#else
const char *short_options = "m:s:F:c:I:j:Q:L:G:f:o:";
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
	//NOTE: SPACES ONLY, bozo
	fprintf(stderr, "cufflinks v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "Usage:   cufflinks [options] <hits.sam>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-m/--inner-dist-mean         the average inner distance between mates              [ default:     45 ]\n");
	fprintf(stderr, "-s/--inner-dist-std-dev      the inner distance standard deviation                 [ default:     20 ]\n");
	fprintf(stderr, "-c/--collapse-thresh         Median depth of cov needed for preassembly collapse   [ default:  10000 ]\n");
	fprintf(stderr, "-F/--min-isoform-fraction    suppress transcripts below this abundance level       [ default:   0.15 ]\n");
	fprintf(stderr, "-f/--min-intron-fraction     filter spliced alignments below this level            [ default:   0.05 ]\n");
	fprintf(stderr, "-a/--junc-alpha              alpha for junction binomial test filter               [ default:   0.01 ]\n");
	fprintf(stderr, "-A/--small-anchor-fraction   percent read overhang taken as 'suspiciously small'   [ default:   0.12 ]\n");
	fprintf(stderr, "-j/--pre-mrna-fraction       suppress intra-intronic transcripts below this level  [ default:   0.15 ]\n");
	fprintf(stderr, "-I/--max-intron-length       ignore alignments with gaps longer than this          [ default: 300000 ]\n");
	fprintf(stderr, "-Q/--min-map-qual            ignore alignments with lower than this mapping qual   [ default:      0 ]\n");
	fprintf(stderr, "-L/--label                   all transcripts have this prefix in their IDs         [ default:   CUFF ]\n");
	fprintf(stderr, "-G/--GTF                     quantitate against reference transcript annotations                      \n");
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
			case 's':
				inner_dist_std_dev = (uint32_t)parseInt(0, "-s/--inner-dist-std-dev arg must be at least 0", print_usage);
				break;
			case 't':
				transcript_score_thresh = parseFloat(-99999999, 0, "-t/--transcript-score-thresh must less than or equal to 0", print_usage);
				break;
			case 'p':
				num_threads = (uint32_t)parseInt(1, "-p/--num-threads arg must be at least 1", print_usage);
				break;
			case 'c':
				collapse_thresh = (uint32_t)parseInt(0, "-c/--collapse-thresh arg must be at least 0", print_usage);
				break;
			case 'F':
				min_isoform_fraction = parseFloat(0, 1.0, "-F/--min-isoform-fraction must be between 0 and 1.0", print_usage);
				break;
			case 'f':
				min_intron_fraction = parseFloat(0, 1.0, "-f/--min-intron-fraction must be between 0 and 1.0", print_usage);
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
	
	max_inner_dist = inner_dist_mean + 7 * inner_dist_std_dev;
	inner_dist_norm = normal(inner_dist_mean, inner_dist_std_dev);
	min_intron_fraction = min_isoform_fraction;
    
	if (ref_gtf_filename != "")
	{
		allow_junk_filtering = false;	
	}
	
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
						  vector<Scaffold>& scaffolds,
						  BundleStats* stats = NULL)
{
	vector<Scaffold> hits;
	
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
    
	filter_introns(bundle.length(), 
				   bundle.left(), 
				   hits, 
				   min_intron_fraction, 
				   false,
				   true);
	
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
#if ASM_VERBOSE
		fprintf (stderr, "%s\tFiltering forward strand\n", bundle_label->c_str());
#endif
		filter_hits(bundle.length(), bundle.left(), fwd_hits);
#if ASM_VERBOSE
		fprintf (stderr, "%s\tFiltering reverse strand\n", bundle_label->c_str());
#endif
		filter_hits(bundle.length(), bundle.left(), rev_hits);
	}
    
	vector<Scaffold> fwd_scaffolds;
	vector<Scaffold> rev_scaffolds;
	
	bool assembled_successfully = false;
	
	if (saw_fwd && saw_rev)
	{
		assembled_successfully |= make_scaffolds(bundle.left(), bundle.length(), fwd_hits, fwd_scaffolds);
        
		assembled_successfully |= make_scaffolds(bundle.left(), bundle.length(), rev_hits, rev_scaffolds);
		
		add_non_shadow_scaffs(fwd_scaffolds, rev_scaffolds, scaffolds, true);
		add_non_shadow_scaffs(rev_scaffolds, fwd_scaffolds, scaffolds, false);
	}
	else
	{
		if (saw_fwd || (!saw_fwd && !saw_rev))
		{
			assembled_successfully |= make_scaffolds(bundle.left(),
													 bundle.length(),
													 fwd_hits,
													 fwd_scaffolds);
			scaffolds.insert(scaffolds.end(),fwd_scaffolds.begin(), fwd_scaffolds.end());
		}
		else
		{
			assembled_successfully |= make_scaffolds(bundle.left(), 
													 bundle.length(), 
													 rev_hits, 
													 rev_scaffolds);
			scaffolds.insert(scaffolds.end(),rev_scaffolds.begin(), rev_scaffolds.end());
		}
	}
	
	// Make sure all the reads are accounted for, including the redundant ones...
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		scaffolds[i].clear_hits();
		for (size_t j = 0; j < bundle.hits().size(); ++j)
		{
			const MateHit& h = bundle.hits()[j];
			scaffolds[i].add_hit(&h);
		}
	}
	
	sort(scaffolds.begin(), scaffolds.end(), scaff_lt);
	
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
								   long double map_mass,
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
	
	avg_read_length /= hits_in_cluster.size();
	
	transfrag_cluster.calculate_abundance(hits_in_cluster, map_mass);
	
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
				
				const Scaffold* transfrag = iso_ab->transfrag();
				assert(transfrag);
				
				double s_len = transfrag->length();
				
				density_per_bp *= (map_mass / 1000000.0); // yields (mass/(length/1000))
				density_per_bp *= (s_len/ 1000.0);
				density_per_bp /= s_len;
				density_per_bp *= avg_read_length;
				//double density_per_bp = (FPKM * (map_mass / 1000000.0) * 1000.0);
                
				
				if (density_score > min_isoform_fraction)
				{
					isoforms.push_back(Isoform(*transfrag,
											   gene_id,
											   (int)isoforms.size() + 1,
											   FPKM, 
											   iso_ab->gamma(),
											   iso_ab->FPKM_conf(),
											   density_per_bp, 
											   density_score));
				}
			}
			
			if (!isoforms.empty())
			{
				Gene g(isoforms, gene.FPKM(), gene.FPKM_conf());
				genes.push_back(g);	
			}
		}
	}
    
}

void quantitate_transcript_clusters(vector<Scaffold>& scaffolds,
									long double map_mass,
									vector<Gene>& genes)
{	
	vector<Scaffold> partials;
	vector<Scaffold> completes;
    
	for (size_t i = 0; i < scaffolds.size(); ++i)
	{
		if (scaffolds[i].has_unknown())
		{
			partials.push_back(scaffolds[i]);
		}
		else
		{
			completes.push_back(scaffolds[i]);
		}
	}
	
	scaffolds = completes;
	
	vector<shared_ptr<Abundance> > abundances;
	foreach(Scaffold& s, scaffolds)
	{
		TranscriptAbundance* pT = new TranscriptAbundance;
		pT->transfrag(&s);
		shared_ptr<Abundance> ab(pT);
		abundances.push_back(ab);
	}
	
	AbundanceGroup transfrags = AbundanceGroup(abundances);
	
	vector<AbundanceGroup> transfrags_by_cluster;
	
	cluster_transcripts<ConnectByExonOverlap>(transfrags,
                                              transfrags_by_cluster);
	
	foreach(AbundanceGroup& cluster, transfrags_by_cluster)
	{
		quantitate_transcript_cluster(cluster, map_mass, genes);
	}
#if ASM_VERBOSE
    fprintf(stderr, "%s\tBundle quantitation complete\n", bundle_label->c_str());
#endif
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
    bundle_label.reset(new string(bundle_label_buf));
    
    fprintf(stderr, "%s\tProcessing new bundle with %d alignments\n", 
            bundle_label->c_str(),
            (int)bundle.hits().size());
    
#if ENABLE_THREADS	
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
	
	vector<Scaffold> scaffolds;
	
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
	quantitate_transcript_clusters(scaffolds, 
								   map_mass, 
								   genes);
    
#if ASM_VERBOSE
    fprintf(stderr, "%s\tFiltering bundle assembly\n", bundle_label->c_str());
#endif
    
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
			
			fprintf(ftrans_abundances,"%s\t%d\t%s\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\n", 
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
					iso.scaffold().length());
			fflush(ftrans_abundances);
			
			num_scaffs_reported++;
		}
		
		fprintf(fgene_abundances,"%s\t%d\t%s\t%d\t%d\t%lg\t%lg\t%lg\n",
				gene.gene_id().c_str(),
				bundle.id(),
				rt.get_name(bundle.ref_id()),
				gene.left(),
				gene.right(),
				gene.FPKM(),
				gene.confidence().low,
				gene.confidence().high);
		fflush(fgene_abundances);
	}
    
	//fprintf(fbundle_tracking, "CLOSE %d\n", bundle.id());
	
	if (ref_gtf_filename != "" && num_scaffs_reported > bundle.ref_scaffolds().size())
	{
		fprintf(stderr, "Error: reported more isoforms than in reference!\n");
		exit(1);
	}
	
#if ASM_VERBOSE
    fprintf(stderr, "%s\tBundle complete\n", bundle_label->c_str());
#endif
    
#if ENABLE_THREADS
	out_file_lock.unlock();
#endif
	
	delete bundle_ptr;
}

bool assemble_hits(BundleFactory& bundle_factory)
{
	fprintf(stderr, "Counting hits in map\n");
	long double map_mass = 0.0;
	
	BadIntronTable bad_introns;
	
	inspect_map(bundle_factory, map_mass, bad_introns);
	
	bundle_factory.load_ref_rnas();
	
	if (ref_gtf_filename == "")
	{
		bundle_factory.bad_intron_table(bad_introns);
	}
	
	fprintf(stderr, "\tTotal map density: %Lf\n", map_mass);
	
	srand(time(0));
	
	HitBundle bundle;
	
	int num_bundles = 0;
	vector<size_t> bundle_sizes;
	vector<size_t> nr_bundle_sizes;
	vector<int> bundle_compatible;
	vector<int> bundle_uncollapsible;
	vector<int> bundle_spans;
	
	RefSequenceTable& rt = bundle_factory.hit_factory().ref_table();
    
	//FILE* fbundle_tracking = fopen("open_bundles", "w");
    
	//FILE* fstats = fopen("bundles.stats", "w");
	FILE* ftrans_abundances = fopen(string(output_dir + "/" + "transcripts.expr").c_str(), "w");
	fprintf(ftrans_abundances,"trans_id\tbundle_id\tchr\tleft\tright\tFPKM\tFMI\tfrac\tFPKM_conf_lo\tFPKM_conf_hi\tcoverage\tlength\n");
	
	FILE* fgene_abundances = fopen(string(output_dir + "/" + "genes.expr").c_str(), "w");
	fprintf(fgene_abundances,"gene_id\tbundle_id\tchr\tleft\tright\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\n");
	FILE* ftranscripts = fopen(string(output_dir + "/" + "transcripts.gtf").c_str(), "w");
    
	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		if (!bundle_factory.next_bundle(*bundle_ptr))
		{
			delete bundle_ptr;
			break;
		}
		
		HitBundle& bundle = *bundle_ptr;
		
		if (bundle.right() - bundle.left() > 3000000)
		{
			fprintf(stderr, "%s\tWarning: large bundle encountered\n", bundle_label->c_str());
		}
		if (bundle.hits().size())
		{
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
            
#if ENABLE_THREADS			
			thread_pool_lock.lock();
			curr_threads++;
			thread_pool_lock.unlock();
            
			thread asmbl(assemble_bundle,
						 boost::cref(rt), 
						 bundle_ptr, 
						 map_mass, 
						 ftranscripts, 
						 fgene_abundances,
						 ftrans_abundances
						 );
#else
			assemble_bundle(boost::cref(rt), 
							bundle_ptr, 
							map_mass, 
							ftranscripts,
							fgene_abundances,
							ftrans_abundances);
#endif			
            
		}
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
	
	return true;
}


void driver(FILE* sam_hit_file, FILE* ref_gtf)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	SAMHitFactory hit_factory(it, rt);
	
	BundleFactory bundle_factory(hit_factory, sam_hit_file, ref_gtf);
	
#if ENABLE_THREDS
	boost::thread asm_thread(assemble_hits,
							 bundle_factory);
	asm_thread.join();
#else	
	assemble_hits(bundle_factory);
#endif
	
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
	
    string sam_hits_file_name = argv[optind++];
	
    // Open the approppriate files
    FILE* sam_hits_file = fopen(sam_hits_file_name.c_str(), "r");
    if (sam_hits_file == NULL)
    {
        fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
                sam_hits_file_name.c_str());
        exit(1);
    }
	
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
    
    driver(sam_hits_file, ref_gtf);
	
	return 0;
}
