#ifndef COMMON_H
#define COMMON_H
/*
 *  common.h
 *  Cufflinks
 *
 *  Created by Cole Trapnell on 11/26/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <cassert>
#include <string>

#include <boost/math/distributions/normal.hpp> 
using boost::math::normal;

#include <boost/foreach.hpp>
#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH


extern uint32_t max_intron_length;
extern uint32_t min_intron_length;

extern uint32_t max_gene_length;
extern int max_partner_dist;
extern int inner_dist_mean;
extern int inner_dist_std_dev;
extern int max_inner_dist;
extern normal inner_dist_norm;
extern double transcript_score_thresh;
extern int olap_radius;
extern float pre_mrna_fraction;

extern int num_threads;

static const int bowtie_overhang_tolerance = 8;
extern float min_isoform_fraction;
extern float min_intron_fraction;
extern float max_phred_err_prob;
extern float high_phred_err_prob;

extern std::string user_label;
extern std::string ref_gtf_filename;
extern std::string output_dir;

extern int collapse_thresh;
extern int microexon_length;

extern bool perform_full_collapse;

extern bool allow_junk_filtering;

extern int max_mle_iterations;
extern int num_importance_samples;

#define ENABLE_THREADS 1
#define ASM_VERBOSE 0

bool gaurd_assembly();

int parseInt(int lower, 
			 const char *errmsg, 
			 void (*print_usage)());

float parseFloat(float lower, 
				 float upper, 
				 const char *errmsg, 
				 void (*print_usage)());

template<typename InputIterator,
		 typename OutputIterator,
		 typename Predicate>
OutputIterator copy_if(InputIterator begin,
					   InputIterator end,
					   OutputIterator destBegin,
					   Predicate p)
{
	while (begin != end)
	{
		if (p(*begin)) *destBegin++ = *begin;
		++begin;
	}
	return destBegin;
}


#define OPT_NUM_IMP_SAMPLES		260
#define OPT_MLE_MAX_ITER		261
#define OPT_FDR					262

#endif
