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
#include <utility>

#include <boost/math/distributions/normal.hpp> 
using boost::math::normal;

#include <boost/foreach.hpp>
#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

// Non-option globals
extern bool final_est_run;
extern bool allow_junk_filtering;
extern bool user_provided_fld;
extern int def_max_frag_len;
extern int max_frag_len;
extern int min_frag_len;

// Behavior options
extern int num_threads;
extern bool no_update_check;
extern bool cuff_quiet;
extern bool cuff_verbose;
extern bool output_fld;
extern bool output_bias_params;

// General options
extern int max_partner_dist;
extern uint32_t max_gene_length;
extern std::string ref_gtf_filename;
extern std::string mask_gtf_filename;
extern std::string output_dir;
extern std::string fasta_dir;
extern std::string library_type;

// Abundance estimation options
extern bool corr_bias;
extern bool corr_multi;
extern bool use_quartile_norm;
extern bool poisson_dispersion;
extern int def_frag_len_mean;
extern int def_frag_len_std_dev;
extern int max_mle_iterations;
extern int num_importance_samples;
extern float min_isoform_fraction;
extern bool use_em;
extern bool cond_prob_collapse;
extern bool use_compat_mass;
extern bool use_total_mass;

// Ref-guided assembly options
extern int overhang_3;
extern int ref_merge_overhang_tolerance;
extern int tile_len;
extern int tile_off;
extern bool enable_faux_reads;
extern bool enable_5_extend;

// Assembly options
extern uint32_t min_intron_length;
extern uint32_t max_intron_length;
extern int olap_radius;
extern int bowtie_overhang_tolerance;
extern int min_frags_per_transfrag;
extern int microexon_length;
extern float pre_mrna_fraction;
extern float high_phred_err_prob;
extern double trim_3_dropoff_frac;
extern double trim_3_avgcov_thresh;
extern double small_anchor_fraction;
extern double binomial_junc_filter_alpha;
extern std::string user_label;
extern long random_seed;
extern bool emit_count_tables;
extern bool use_fisher_covariance;
extern bool split_variance;
extern bool bootstrap;
extern int num_bootstrap_samples;
extern double bootstrap_fraction;
extern double bootstrap_delta_gap;
extern int max_frags_per_bundle;

// SECRET OPTIONS: 
// These options are just for instrumentation and benchmarking code

extern bool no_read_pairs;
extern float read_skip_fraction;
extern int trim_read_length;
extern double mle_accuracy;

// END SECRET OPTIONS

#define ASM_VERBOSE 0
#define ENABLE_THREADS 1

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

enum BundleMode
{
	HIT_DRIVEN,
	REF_DRIVEN,
	REF_GUIDED
};
extern BundleMode bundle_mode;
extern BundleMode init_bundle_mode;

enum BiasMode
{
	SITE,
	VLMM,
	POS,
	POS_VLMM,
    POS_SITE
};
extern BiasMode bias_mode;

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
	std::vector<double> _pdf;
	std::vector<double> _cdf;
	int _mode;
	double _mean;
    double _std_dev;
	int _min;
	int _max;
	
public:
	EmpDist(std::vector<double>& pdf, std::vector<double>& cdf, int mode, double mean, double std_dev, int min, int max)
	: _pdf(pdf), _cdf(cdf), _mode(mode), _mean(mean), _std_dev(std_dev), _min(min), _max(max) {}
	
	void pdf(std::vector<double>& pdf)	{ _pdf = pdf; }
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
	
	void cdf(std::vector<double>& cdf)	{ _cdf = cdf; }
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
class MultiReadTable;

class MassDispersionModel;

struct LocusCount
{
    LocusCount(std::string ld, double c, int nt) : 
        locus_desc(ld), count(c), num_transcripts(nt) {}
    std::string locus_desc;
    double count;
    int num_transcripts;
};

class ReadGroupProperties
{
public:
    
    ReadGroupProperties(); 
    
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
    
    long double normalized_map_mass() const { return _norm_map_mass; }
    void normalized_map_mass(long double p)  { _norm_map_mass = p; }  
    
    boost::shared_ptr<EmpDist const> frag_len_dist() const { return _frag_len_dist; }
    void frag_len_dist(boost::shared_ptr<EmpDist const> p)  { _frag_len_dist = p; }  
    
	boost::shared_ptr<BiasLearner const> bias_learner() const { return _bias_learner; }
    void bias_learner(boost::shared_ptr<BiasLearner const> bl)  { _bias_learner = bl; } 
	
    void mass_scale_factor(double sf) { _mass_scaling_factor = sf; }
    double mass_scale_factor() const  { return _mass_scaling_factor; }
    
    void complete_fragments(bool c)  { _complete_fragments = c; }
    bool complete_fragments() const { return _complete_fragments; }
    
    double scale_mass(double unscaled_mass) const 
    { 
        if (_mass_scaling_factor == 0)
            return unscaled_mass;
        
        return unscaled_mass * (1.0 / _mass_scaling_factor);
    }
    
    boost::shared_ptr<const MassDispersionModel> mass_dispersion_model() const 
    { 
        return _mass_dispersion_model; 
    };
    
    void mass_dispersion_model(boost::shared_ptr<const MassDispersionModel> nm) 
    { 
        _mass_dispersion_model = nm; 
    }
    
    const std::vector<LocusCount>& common_scale_counts() { return _common_scale_counts; }
    void common_scale_counts(const std::vector<LocusCount>& counts) { _common_scale_counts = counts; }
    
	boost::shared_ptr<MultiReadTable> multi_read_table() const {return _multi_read_table; }	
	void multi_read_table(boost::shared_ptr<MultiReadTable> mrt) { _multi_read_table = mrt;	}
	
private:
    
    Strandedness _strandedness;
    StandardMateOrientation _std_mate_orient;
	MateStrandMapping _mate_strand_mapping;
    Platform _platform;
    long double _total_map_mass;
    long double _norm_map_mass;
    boost::shared_ptr<EmpDist const> _frag_len_dist;
	boost::shared_ptr<BiasLearner const> _bias_learner;
	boost::shared_ptr<MultiReadTable> _multi_read_table;
    
    double _mass_scaling_factor;
    boost::shared_ptr<const MassDispersionModel> _mass_dispersion_model;
    std::vector<LocusCount> _common_scale_counts;
    
    bool _complete_fragments;
};

extern std::map<std::string, ReadGroupProperties> library_type_table;

extern const ReadGroupProperties* global_read_properties;

void print_library_table();
void init_library_table();


template<typename T>
std::string cat_strings(const T& container, const char* delimiter=",")
{
    std::string cat;
	if (container.empty())
	{
		cat = "";
	}
	else
	{
		typename T::const_iterator itr = container.begin();
		//cat = *(itr);
		for (; itr != container.end(); itr++)
		{
			if (!(*itr).empty()) {
				if (!cat.empty()) cat += delimiter;
				cat += *itr; 
            }
		}
	}
    
	return cat;
}

#define OPT_NUM_IMP_SAMPLES         260
#define OPT_MLE_MAX_ITER            261
#define OPT_FDR                     262
#define OPT_LIBRARY_TYPE            263
#define OPT_OVERHANG_TOLERANCE      264
#define OPT_MAX_BUNDLE_LENGTH       265
#define OPT_MIN_FRAGS_PER_TRANSFRAG 266
#define OPT_BIAS_MODE               267
#define OPT_MIN_INTRON_LENGTH       268
#define OPT_3_PRIME_AVGCOV_THRESH	269
#define OPT_3_PRIME_DROPOFF_FRAC    270
#define OPT_POISSON_DISPERSION      271
#define OPT_NO_UPDATE_CHECK         272
#define OPT_OUTPUT_FLD              273
#define OPT_OUTPUT_BIAS_PARAMS      274
#define OPT_USE_EM                  275
#define OPT_COLLAPSE_COND_PROB      276
#define OPT_RANDOM_SEED             277
#define OPT_NO_FAUX_READS           278
#define OPT_3_OVERHANG_TOLERANCE    279
#define OPT_INTRON_OVERHANG_TOLERANCE 280
#define OPT_EMIT_COUNT_TABLES       281
#define OPT_USE_COMPAT_MASS         282
#define OPT_USE_TOTAL_MASS          283
#define OPT_USE_FISHER_COVARIANCE   284
#define OPT_USE_EMPIRICAL_COVARIANCE   285
#define OPT_SPLIT_MASS              286
#define OPT_SPLIT_VARIANCE          287
#define OPT_BOOTSTRAP               288
#define OPT_NUM_BOOTSTRAP_SAMPLES   289
#define OPT_BOOTSTRAP_FRACTION      290
#define OPT_TILE_LEN                291
#define OPT_TILE_SEP                292
#define OPT_NO_5_EXTEND             293
#define OPT_MAX_FRAGS_PER_BUNDLE    294
#define OPT_READ_SKIP_FRACTION      295
#define OPT_NO_READ_PAIRS           296
#define OPT_TRIM_READ_LENGTH        297
#define OPT_MAX_DELTA_GAP           298
#define OPT_MLE_MIN_ACC             299
#endif
