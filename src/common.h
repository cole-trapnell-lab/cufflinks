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

#include <boost/thread.hpp>

extern uint32_t max_intron_length;
extern uint32_t min_intron_length;

extern uint32_t max_gene_length;
extern int max_partner_dist;
extern int def_frag_len_mean;
extern int def_frag_len_std_dev;
extern int def_max_frag_len;
extern int max_frag_len;
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
extern std::string mask_gtf_filename;
extern std::string output_dir;
extern std::string fasta_dir;

extern int collapse_thresh;
extern int microexon_length;

extern bool perform_full_collapse;

extern bool allow_junk_filtering;

extern int max_mle_iterations;
extern int num_importance_samples;

extern double small_anchor_fraction;

extern double binomial_junc_filter_alpha;

extern std::string library_type;

#define ENABLE_THREADS 1

#if ENABLE_THREADS
extern boost::thread_specific_ptr<std::string> bundle_label; // for consistent, traceable logging
#else
extern boost::shared_ptr<std::string> bundle_label;
#endif

#define ASM_VERBOSE 0

bool gaurd_assembly();

int parseInt(int lower, 
			 const char *errmsg, 
			 void (*print_usage)());

float parseFloat(float lower, 
				 float upper, 
				 const char *errmsg, 
				 void (*print_usage)());

void encode_seq(const std::string seqStr, char* seq, char* c_seq);
int mkpath(const char *s, mode_t mode);


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

enum Strandedness 
{
    UNKNOWN_STRANDEDNESS,
    STRANDED_PROTOCOL,
    UNSTRANDED_PROTOCOL
};

enum StandardMateOrientation
{
    UNKNOWN_MATE_ORIENTATION,
    MATES_POINT_TOWARD,
    MATES_POINT_SAME,
    MATES_POINT_AWAY,
    UNPAIRED,
};

enum Platform
{
    UNKNOWN_PLATFORM,
    ILLUMINA,
    SOLID
};

class ReadGroupProperties
{
public:
    
    ReadGroupProperties() : 
    _strandedness(UNKNOWN_STRANDEDNESS), 
    _std_mate_orient(UNKNOWN_MATE_ORIENTATION),
    _platform(UNKNOWN_PLATFORM){} 
    
    Strandedness strandedness() const { return _strandedness; }
    void strandedness(Strandedness s) { _strandedness = s; }
    
    StandardMateOrientation std_mate_orientation() const { return _std_mate_orient; }
    void std_mate_orientation(StandardMateOrientation so)  { _std_mate_orient = so; }
    
    Platform platform() const { return _platform; }
    void platform(Platform p)  { _platform = p; }    
private:
    
    Strandedness _strandedness;
    StandardMateOrientation _std_mate_orient;
    Platform _platform;
};

extern std::map<std::string, ReadGroupProperties> library_type_table;

extern const ReadGroupProperties* global_read_properties;

void print_library_table();
void init_library_table();

#define OPT_NUM_IMP_SAMPLES		260
#define OPT_MLE_MAX_ITER		261
#define OPT_FDR					262
#define OPT_LIBRARY_TYPE        263

#endif
