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
//int max_mate_inner_dist = -1; 
//
//uint32_t min_anchor_len = 5;
uint32_t min_intron_length = 50;
uint32_t max_intron_length = 300000;
//uint32_t min_exon_length = 100; 

uint32_t max_gene_length = 3000000;
int max_partner_dist = 50000;
int inner_dist_mean = 45;
int inner_dist_std_dev = 40;
int max_inner_dist = 105;
int olap_radius = 50;

float min_isoform_fraction = 0.05;
float min_intron_fraction = 0.05;
float pre_mrna_fraction = 0.15;
float high_phred_err_prob = 0.50; // about MAPQ = 3

double transcript_score_thresh = -0.693;

normal inner_dist_norm;

int num_threads = 1;

float max_phred_err_prob = 1.0;

std::string user_label = "CUFF";
std::string ref_gtf_filename = "";
std::string output_dir = "./";

int collapse_thresh = 10;

int microexon_length = 25;

bool perform_full_collapse = true;

bool allow_junk_filtering = true;

int max_mle_iterations = 5000;
int num_importance_samples = 1000;

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
