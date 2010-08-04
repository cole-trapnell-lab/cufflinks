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

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <libgen.h>
#include <string.h>

#include "getopt.h"
#include "common.h"

using namespace std;

//int insert_len = 250;
//int insert_len_std_dev = 20;
//
//uint32_t min_anchor_len = 5;
uint32_t min_intron_length = 50;
uint32_t max_intron_length = 300000;
//uint32_t min_exon_length = 100; 

uint32_t max_gene_length = 3000000;
int max_partner_dist = 50000;
int def_frag_len_mean = 195;
int def_frag_len_std_dev = 80;
int def_max_frag_len = 1000;
int max_frag_len = 1000;
int olap_radius = 50;

float min_isoform_fraction = 0.1;
//float min_isoform_fraction = 0.1;
float pre_mrna_fraction = 0.25;
float high_phred_err_prob = 0.50; // about MAPQ = 3

double transcript_score_thresh = -0.693;

int num_threads = 1;

float max_phred_err_prob = 1.0;

std::string user_label = "CUFF";
std::string ref_gtf_filename = "";
std::string mask_gtf_filename = "";
std::string output_dir = "./";
std::string fasta_dir;

int collapse_thresh = 10;

int microexon_length = 25;

bool perform_full_collapse = true;

bool allow_junk_filtering = true;

int max_mle_iterations = 5000;
int num_importance_samples = 1000;

double small_anchor_fraction = 7 / 75.0;
double binomial_junc_filter_alpha = 0.001;

string default_library_type = "illumina-unstranded-paired-end";
string library_type = "";

map<string, ReadGroupProperties> library_type_table;
const ReadGroupProperties* global_read_properties = NULL;

#if ENABLE_THREADS
boost::thread_specific_ptr<std::string> bundle_label;
#else
boost::shared_ptr<std::string> bundle_label;
#endif

extern void print_usage();

bool gaurd_assembly()
{
	return ref_gtf_filename == "";
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
    ReadGroupProperties std_illumina_paired;
    std_illumina_paired.platform(ILLUMINA);
    std_illumina_paired.std_mate_orientation(MATES_POINT_TOWARD);
    std_illumina_paired.strandedness(UNSTRANDED_PROTOCOL);
    
    library_type_table["illumina-unstranded-paired-end"] = std_illumina_paired;
   
    ReadGroupProperties solid_fragment;
    solid_fragment.platform(SOLID);
    solid_fragment.std_mate_orientation(UNPAIRED);
    solid_fragment.strandedness(STRANDED_PROTOCOL);
    
    library_type_table["solid-fragment"] = solid_fragment;
    
    ReadGroupProperties illumina_fragment;
    illumina_fragment.platform(ILLUMINA);
    illumina_fragment.std_mate_orientation(UNPAIRED);
    illumina_fragment.strandedness(UNSTRANDED_PROTOCOL);
    
    library_type_table["illumina-fragment"] = illumina_fragment;
    
    global_read_properties = &(library_type_table.find(default_library_type)->second);
}

void print_library_table()
{
    fprintf (stderr, "Supported library types:\n");
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
    
	for (int i = 0; i < seqStr.length(); ++i)
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

