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

using namespace std;

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

#include <boost/shared_ptr.hpp>

extern bool final_est_run;
extern bool corr_bias;
extern bool ref_driven;

extern uint32_t max_intron_length;
extern uint32_t min_intron_length;

extern uint32_t max_gene_length;

extern int max_partner_dist;
extern int def_frag_len_mean;
extern int def_frag_len_std_dev;
extern int def_max_frag_len;
extern int max_frag_len;
extern int min_frag_len;

extern double transcript_score_thresh;
extern int olap_radius;
extern float pre_mrna_fraction;

extern int num_threads;

extern int bowtie_overhang_tolerance;
extern float min_isoform_fraction;
//extern float min_isoform_fraction;
extern float max_phred_err_prob;
extern float high_phred_err_prob;

extern std::string user_label;
extern std::string ref_gtf_filename;
extern std::string mask_gtf_filename;
extern std::string output_dir;
extern std::string fasta_dir;
extern std::string bias_mode;

extern int microexon_length;
extern bool cuff_verbose;
extern bool cuff_quiet;
extern bool perform_full_collapse;

extern bool allow_junk_filtering;

extern bool use_quartile_norm;

extern int max_mle_iterations;
extern int num_importance_samples;

extern double small_anchor_fraction;

extern double binomial_junc_filter_alpha;

extern std::string library_type;

extern int min_frags_per_transfrag;

#define ADAM_MODE 0
#define ASM_VERBOSE 0
#define ENABLE_THREADS 0

#if ENABLE_THREADS
extern boost::thread_specific_ptr<std::string> bundle_label; // for consistent, traceable logging
#else
extern boost::shared_ptr<std::string> bundle_label;
#endif


bool gaurd_assembly();

void asm_verbose(const char* fmt,...);
void verbose_msg(const char* fmt,...); 

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

enum MateStrandMapping
{
	FF,
	FR,
	RF, // This is really FR with first-strandedness
	RR // This is really FF with first-strandedness
};

enum Platform
{
    UNKNOWN_PLATFORM,
    ILLUMINA,
    SOLID
};

class EmpDist
{
	//Vectors only valid between min and max!
	vector<double> _pdf;
	vector<double> _cdf;
	int _mode;
	int _max;
	int _min;
    double _mean;
    double _std_dev;
	
public:
	
	void pdf(vector<double>& pdf)	{ _pdf = pdf; }
	double pdf(int l) const
	{
		if (!valid_len(l))
			return 0.0;
		return _pdf[l];
	}
	
	// pdf renomalized over the lengths <= r
	double npdf(int l, int r) const
 	{
		if (!valid_len(l))
			return 0.0;
		
		if (r > _max || r == 0)
			return pdf(l);
		
		return pdf(l)/cdf(r);
	}
	
	void cdf(vector<double>& cdf)	{ _cdf = cdf; }
	double cdf(int l) const
	{
		if (l > _max)
			return 1.0;
        if (l < 0)
            return 0.0;
		return _cdf[l];
	}
	
	bool valid_len(int l) const { return (l >= _min && l <= _max); }
	bool too_short(int l) const { return (l < _min); }
	
	void mode(int mode)				{ _mode = mode; }
	int mode() const				{ return _mode; }
	
	void max(int max)				{ _max = max;  }
	int max() const					{ return _max; }
	
	void min(int min)				{ _min = min;  }
	int min() const					{ return _min; }
    
    void mean(double mean)				{ _mean = mean;  }
	double mean() const					{ return _mean; }
    
    void std_dev(double std_dev)				{ _std_dev = std_dev;  }
	double std_dev() const					{ return _std_dev; }
};

class BiasLearner;

class ReadGroupProperties
{
public:
    
    ReadGroupProperties() : 
    _strandedness(UNKNOWN_STRANDEDNESS), 
    _std_mate_orient(UNKNOWN_MATE_ORIENTATION),
    _platform(UNKNOWN_PLATFORM),
    _total_map_mass(0.0) {} 
    
    Strandedness strandedness() const { return _strandedness; }
    void strandedness(Strandedness s) { _strandedness = s; }
    
    StandardMateOrientation std_mate_orientation() const { return _std_mate_orient; }
    void std_mate_orientation(StandardMateOrientation so)  { _std_mate_orient = so; }
    
	MateStrandMapping mate_strand_mapping() const { return _mate_strand_mapping; }
	void mate_strand_mapping(MateStrandMapping msm) { _mate_strand_mapping = msm; }
	
    Platform platform() const { return _platform; }
    void platform(Platform p)  { _platform = p; }   
    
    long double total_map_mass() const { return _total_map_mass; }
    void total_map_mass(long double p)  { _total_map_mass = p; }  
    
    boost::shared_ptr<EmpDist const> frag_len_dist() const { return _frag_len_dist; }
    void frag_len_dist(boost::shared_ptr<EmpDist const> p)  { _frag_len_dist = p; }  
    
	boost::shared_ptr<BiasLearner const> bias_learner() const { return _bias_learner; }
    void bias_learner(boost::shared_ptr<BiasLearner const> bl)  { _bias_learner = bl; } 
	
private:
    
    Strandedness _strandedness;
    StandardMateOrientation _std_mate_orient;
	MateStrandMapping _mate_strand_mapping;
    Platform _platform;
    long double _total_map_mass;
    boost::shared_ptr<EmpDist const> _frag_len_dist;
	boost::shared_ptr<BiasLearner const> _bias_learner;
};

extern std::map<std::string, ReadGroupProperties> library_type_table;

extern const ReadGroupProperties* global_read_properties;

void print_library_table();
void init_library_table();

#define OPT_NUM_IMP_SAMPLES         260
#define OPT_MLE_MAX_ITER            261
#define OPT_FDR                     262
#define OPT_LIBRARY_TYPE            263
#define OPT_OVERHANG_TOLERANCE      264
#define OPT_MAX_BUNDLE_LENGTH       265
#define OPT_MIN_FRAGS_PER_TRANSFRAG 266
#define OPT_BIAS_MODE               267
#define OPT_MIN_INTRON_LENGTH       268

#endif
