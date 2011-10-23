/*
 *  common.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <stdarg.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <libgen.h>
#include <string.h>

#include "getopt.h"
#include "common.h"
#include "replicates.h"

using namespace std;


// Non-option globals
bool final_est_run = true;
bool allow_junk_filtering = true;
bool user_provided_fld = false;

// Behavior options
int num_threads = 1;
bool no_update_check = false;
bool cuff_quiet = false;
#if ASM_VERBOSE
bool cuff_verbose = true;
#else
bool cuff_verbose = false;
#endif
bool output_fld = false;
bool output_bias_params = false;

// General options
BundleMode bundle_mode = HIT_DRIVEN;
BundleMode init_bundle_mode = HIT_DRIVEN;
int max_partner_dist = 50000;
uint32_t max_gene_length = 3500000;
std::string ref_gtf_filename = "";
std::string mask_gtf_filename = "";
std::string output_dir = "./";
std::string fasta_dir;
string default_library_type = "fr-unstranded";
string library_type = default_library_type;


// Abundance estimation options
bool corr_bias = false;
bool corr_multi = false;
bool use_quartile_norm = false;
bool poisson_dispersion = false;
BiasMode bias_mode = POS_VLMM;
int def_frag_len_mean = 200;
int def_frag_len_std_dev = 80;
int def_max_frag_len = 800;
int max_frag_len = 800;
int min_frag_len = 1;
float min_isoform_fraction = 0.1;
int max_mle_iterations = 5000;
int num_importance_samples = 10000;
bool use_compat_mass = false;
bool use_total_mass = false;


// Ref-guided assembly options
int overhang_3 = 600;
int ref_merge_overhang_tolerance = 30;
int tile_len = 405;
int tile_off = 15;
bool enable_faux_reads = true;
bool enable_5_extend = true;

// Assembly options
uint32_t min_intron_length = 50;
uint32_t max_intron_length = 300000;
int olap_radius = 50;
int bowtie_overhang_tolerance = 8; // Typically don't need to change this, except in special cases, such as meta-assembly.
int min_frags_per_transfrag = 10;
int microexon_length = 25;
float pre_mrna_fraction = 0.15;
float high_phred_err_prob = 0.50; // about MAPQ = 3
double small_anchor_fraction = 7 / 75.0;
double binomial_junc_filter_alpha = 0.001;
double trim_3_dropoff_frac = .1;
double trim_3_avgcov_thresh = 10.0;
std::string user_label = "CUFF";

bool use_em = true;
bool cond_prob_collapse = true;

bool emit_count_tables = false;
bool use_fisher_covariance = true;
bool split_variance = false;
bool bootstrap = true;
int num_bootstrap_samples = 20;
double bootstrap_fraction = 1.0;
double bootstrap_delta_gap = 0.001;
int max_frags_per_bundle = 1000000;

// SECRET OPTIONS: 
// These options are just for instrumentation and benchmarking code

float read_skip_fraction = 0.0;
bool no_read_pairs = false;
int trim_read_length = -1;
double mle_accuracy = 1e-6;

// END SECRET OPTIONS


map<string, ReadGroupProperties> library_type_table;
const ReadGroupProperties* global_read_properties = NULL;

#if ENABLE_THREADS
boost::thread_specific_ptr<std::string> bundle_label;
#else
boost::shared_ptr<std::string> bundle_label;
#endif

long random_seed = 0;

extern void print_usage();

bool gaurd_assembly()
{
	return ref_gtf_filename == "";
}

void asm_verbose(const char* fmt,...)
{
#if !ASM_VERBOSE
	return;
#endif
     va_list argp;
     va_start(argp, fmt);
     vfprintf(stderr, fmt, argp);
     va_end(argp);
}

void verbose_msg(const char* fmt,...) {

	if (!cuff_verbose)
		return;
	
   va_list argp;
   va_start(argp, fmt);
   vfprintf(stderr, fmt, argp);
   va_end(argp);
}


/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
 

int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
float parseFloat(float lower, float upper, const char *errmsg, void (*print_usage)()) {
    float l;
    l = (float)atof(optarg);
	
    if (l < lower) {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    if (l > upper)
    {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    return l;
	
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/* Function with behaviour like `mkdir -p'  */
/* found at: http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/ */

int mkpath(const char *s, mode_t mode)
{
    char *q, *r = NULL, *path = NULL, *up = NULL;
    int rv;
    
    rv = -1;
    if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
        return (0);
    
    if ((path = strdup(s)) == NULL)
        exit(1);
    
    if ((q = strdup(s)) == NULL)
        exit(1);
    
    if ((r = dirname(q)) == NULL)
        goto out;
    
    if ((up = strdup(r)) == NULL)
        exit(1);
    
    if ((mkpath(up, mode) == -1) && (errno != EEXIST))
        goto out;
    
    if ((mkdir(path, mode) == -1) && (errno != EEXIST))
        rv = -1;
    else
        rv = 0;
    
out:
    if (up != NULL)
        free(up);
    free(q);
    free(path);
    return (rv);
}

void init_library_table()
{
    ReadGroupProperties fr_unstranded;
    fr_unstranded.platform(UNKNOWN_PLATFORM);
	fr_unstranded.mate_strand_mapping(FR);
    fr_unstranded.std_mate_orientation(MATES_POINT_TOWARD);
    fr_unstranded.strandedness(UNSTRANDED_PROTOCOL);
    
    library_type_table["fr-unstranded"] = fr_unstranded;
        	
	ReadGroupProperties fr_firststrand;
    fr_firststrand.platform(UNKNOWN_PLATFORM);
	fr_firststrand.mate_strand_mapping(RF);
    fr_firststrand.std_mate_orientation(MATES_POINT_TOWARD);
    fr_firststrand.strandedness(STRANDED_PROTOCOL);
	
    library_type_table["fr-firststrand"] = fr_firststrand;

	ReadGroupProperties fr_secondstrand;
    fr_secondstrand.platform(UNKNOWN_PLATFORM);
	fr_secondstrand.mate_strand_mapping(FR);
    fr_secondstrand.std_mate_orientation(MATES_POINT_TOWARD);
    fr_secondstrand.strandedness(STRANDED_PROTOCOL);
	
    library_type_table["fr-secondstrand"] = fr_secondstrand;
	
	ReadGroupProperties ff_unstranded;
    ff_unstranded.platform(UNKNOWN_PLATFORM);
	ff_unstranded.mate_strand_mapping(FF);
    ff_unstranded.std_mate_orientation(MATES_POINT_TOWARD);
    ff_unstranded.strandedness(UNSTRANDED_PROTOCOL);
    
    library_type_table["ff-unstranded"] = ff_unstranded;
	
	ReadGroupProperties ff_firststrand;
    ff_firststrand.platform(UNKNOWN_PLATFORM);
	ff_firststrand.mate_strand_mapping(FF);
    ff_firststrand.std_mate_orientation(MATES_POINT_TOWARD);
    ff_firststrand.strandedness(STRANDED_PROTOCOL);
	
    library_type_table["ff-firststrand"] = ff_firststrand;
	
	ReadGroupProperties ff_secondstrand;
    ff_secondstrand.platform(UNKNOWN_PLATFORM);
	ff_secondstrand.mate_strand_mapping(RR);
    ff_secondstrand.std_mate_orientation(MATES_POINT_TOWARD);
    ff_secondstrand.strandedness(STRANDED_PROTOCOL);
	
    library_type_table["ff-secondstrand"] = ff_secondstrand;
    
    ReadGroupProperties transfrags;
    transfrags.platform(UNKNOWN_PLATFORM);
	transfrags.mate_strand_mapping(FR);
    transfrags.std_mate_orientation(MATES_POINT_TOWARD);
    transfrags.strandedness(UNSTRANDED_PROTOCOL);
	transfrags.complete_fragments(true);
    
    library_type_table["transfrags"] = transfrags;
	
    //global_read_properties = &(library_type_table.find(default_library_type)->second);
}

void print_library_table()
{
    fprintf (stderr, "\nSupported library types:\n");
    for (map<string, ReadGroupProperties>::const_iterator itr = library_type_table.begin();
         itr != library_type_table.end();
         ++itr)
    {
        if (itr->first == default_library_type)
        {
            fprintf(stderr, "\t%s (default)\n", itr->first.c_str());
        }
        else
        {
            fprintf(stderr, "\t%s\n", itr->first.c_str());
        }
    }
}


// c_seq is complement, *NOT* REVERSE complement
void encode_seq(const string seqStr, char* seq, char* c_seq)
{
    
	for (size_t i = 0; i < seqStr.length(); ++i)
	{
		switch(seqStr[i])
		{
			case 'A' : 
			case 'a' : seq[i] = 0; c_seq[i] = 3; break;
			case 'c' : 
			case 'C' : seq[i] = 1; c_seq[i] = 2; break;
			case 'G' :
			case 'g' : seq[i] = 2; c_seq[i] = 1; break;
			case 'T' :
			case 't' : seq[i] = 3; c_seq[i] = 0; break;
			default  : seq[i] = 4; c_seq[i] = 4; break; // N
		}
	}
}


ReadGroupProperties::ReadGroupProperties() : 
    _strandedness(UNKNOWN_STRANDEDNESS), 
    _std_mate_orient(UNKNOWN_MATE_ORIENTATION),
    _platform(UNKNOWN_PLATFORM),
    _total_map_mass(0.0),
    _norm_map_mass(0.0),
    _mass_scaling_factor(1.0),
    _complete_fragments(false)
{
    _mass_dispersion_model = boost::shared_ptr<MassDispersionModel const>(new PoissonDispersionModel);
} 
