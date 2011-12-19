#ifndef ABUNDANCES_H
#define ABUNDANCES_H
/*
 *  abundances.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 4/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "hits.h"
#include "scaffolds.h"
#include "bundles.h"
#include "biascorrection.h"

namespace ublas = boost::numeric::ublas;

struct ConfidenceInterval
{
	ConfidenceInterval(double Low = 0.0, double High = 0.0) 
	: low(Low), high(High) {}
	double low;
	double high;
};

enum AbundanceStatus { NUMERIC_OK, NUMERIC_FAIL, NUMERIC_LOW_DATA, NUMERIC_HI_DATA };

class Abundance
{
public:
	virtual ~Abundance() {}
	
	// Status of the numerical calculation performed on this object.  Safe to
	// do testing only if status == NUMERIC_OK
	virtual AbundanceStatus status() const = 0;
	virtual void status(AbundanceStatus s) = 0;
	
	// Fragments Per Kilbase of transcript per Million fragments mapped
	virtual double				FPKM() const = 0;
    virtual void				FPKM(double fpkm) = 0;
	virtual double				FPKM_variance() const = 0;
	virtual void				FPKM_variance(double v) = 0;
	
	virtual ConfidenceInterval	FPKM_conf() const = 0;
	virtual void FPKM_conf(const ConfidenceInterval& cf) = 0;
	
	// gamma is a fixed property of each transcript or transcript group.  It's
	// the probability that one would draw a fragment from this object, scaled
	// to an arbitrary locus' worth of fragments.
	virtual double			gamma() const = 0;
	virtual void			gamma(double g) = 0;
	
	// Kappa is only really meaningful when this Abundance record is part of a
	// group - it's the relative abundance of this object within the larger
	// group.
	virtual double			kappa() const	= 0;
	virtual void			kappa(double k)	= 0;
	
	virtual double			num_fragments() const = 0;
	virtual void			num_fragments(double nf) = 0;
    
    virtual double          mass_fraction() const = 0;
	virtual void            mass_fraction(double mf) = 0;
    
    virtual double          mass_variance() const = 0;
	virtual void            mass_variance(double mv) = 0;
    
	virtual double			effective_length() const= 0;
	virtual void			effective_length(double el) = 0;
	
	virtual const vector<double>*	cond_probs() const		{ return NULL; }
	virtual void					cond_probs(vector<double>* cp) = 0;
	
	// The structural information for the object, if defined.
	virtual shared_ptr<Scaffold> transfrag() const		{ return shared_ptr<Scaffold>(); }
	
	
	virtual set<string>		gene_id() const = 0;
	virtual set<string>		gene_name() const = 0;
	virtual set<string>		tss_id() const = 0;
	virtual set<string>		protein_id() const = 0;
	
	virtual const string&	description() const = 0;
	virtual void			description(const string& d) = 0;
	
	virtual const string&	locus_tag() const = 0;
	virtual void			locus_tag(const string& L) = 0;
	
	virtual const string&	reference_tag() const = 0;
	virtual void			reference_tag(const string& r) = 0;
};

class TranscriptAbundance : public Abundance
{
public:
	
	TranscriptAbundance() : 
		_status(NUMERIC_OK), 
		_transfrag(shared_ptr<Scaffold>()), 
		_FPKM(0), 
		_FPKM_variance(0),
		_gamma(0), 
		_kappa(1.0), 
		_num_fragments(0),
		_eff_len(0),
		_cond_probs(NULL),
        _sample_mass_fraction(0.0),
        _sample_mass_variance(0.0){}
	
	TranscriptAbundance(const TranscriptAbundance& other)
	{
		_status = other._status;
		_transfrag = other._transfrag;
		_FPKM = other._FPKM;
		_FPKM_conf = other._FPKM_conf;
		_gamma = other._gamma;
		_num_fragments = other._num_fragments;
		_eff_len = other._eff_len;
		_cond_probs = other._cond_probs;
        _sample_mass_fraction = other._sample_mass_fraction;
        _sample_mass_variance = other._sample_mass_variance;
	}
	
	~TranscriptAbundance()
	{
		if (_cond_probs != NULL)
		{
			delete _cond_probs;
			_cond_probs = NULL;
		}
	}
	
	AbundanceStatus status() const			{ return _status; }
	void status(AbundanceStatus s)			{ _status = s; }
	
	double FPKM() const						{ return _FPKM; }
    void FPKM(double fpkm)                  
	{ 
		_FPKM = fpkm;
		_transfrag->fpkm(fpkm);
	}
	double FPKM_variance() const			{ return _FPKM_variance; }
	void   FPKM_variance(double v);			
	
	ConfidenceInterval	FPKM_conf() const   { return _FPKM_conf; }
	void FPKM_conf(const ConfidenceInterval& cf) { _FPKM_conf = cf; }
	
	double gamma() const					{ return _gamma; }
	void gamma(double g)					{ assert(!isnan(g)); _gamma = g; };
	
	double kappa() const					{ return _kappa; }
	void kappa(double k)					{ _kappa = k; }
	
	double num_fragments() const			{ return _num_fragments; }
	void num_fragments(double nf)			{ assert (!isnan(nf)); _num_fragments = nf; }
    
	double mass_fraction() const			{ return _sample_mass_fraction; }
	void mass_fraction(double mf)			{ _sample_mass_fraction = mf; }
	
    double mass_variance() const			{ return _sample_mass_variance; }
	void mass_variance(double mv)			{ _sample_mass_variance = mv; }
	
	void transfrag(shared_ptr<Scaffold> tf)		{ _transfrag = tf; }
	shared_ptr<Scaffold> transfrag() const		{ return _transfrag; }
	
	double effective_length() const			{ return _eff_len; }
	void effective_length(double el)		{ _eff_len = el; }
	
	const vector<double>* cond_probs() const	{ return _cond_probs; }
	void cond_probs(vector<double>* cp) 	
	{ 
		if(_cond_probs != NULL) { delete _cond_probs; };
		_cond_probs = cp;
	}
	
	
	set<string> gene_id() const	
	{
		if (_transfrag)
		{
			set<string> s;
			s.insert(_transfrag->annotated_gene_id());
			return s;
		}
		else 
		{
			assert (false);
			return set<string>();
		}
	}
	
	set<string> gene_name() const	
	{
		if (_transfrag)
		{
			set<string> s;
			s.insert(_transfrag->annotated_gene_name());
			return s;
		}
		else 
		{
			assert (false);
			return set<string>();
		}
	}
	
	set<string> tss_id() const	
	{
		if (_transfrag)
		{
			set<string> s;
			s.insert(_transfrag->annotated_tss_id());
			return s;
		}
		else 
		{
			assert (false);
			return set<string>();
		}
	}
	
	set<string> protein_id() const	
	{
		if (_transfrag)
		{
			set<string> s;
			s.insert(_transfrag->annotated_protein_id());
			return s;
		}
		else 
		{
			assert (false);
			return set<string>();
		}
	}
	
	virtual const string&	description() const	{ return _description; }
	virtual void			description(const string& d) { _description = d; }
	
	virtual const string&	locus_tag() const { return _locus_tag; }
	virtual void			locus_tag(const string& L) { _locus_tag = L; }
	
	virtual const string&	reference_tag() const { return _ref_tag; }
	virtual void			reference_tag(const string& r) { _ref_tag = r; } 
	
private:
	
	void calculate_FPKM_err_bar(double variance);
	
	AbundanceStatus _status;
	shared_ptr<Scaffold> _transfrag;
	double _FPKM;
	double _FPKM_variance;
	ConfidenceInterval _FPKM_conf;
	double _gamma;
	double _kappa;
	double _num_fragments;
	double _eff_len;
	vector<double>* _cond_probs;
	
	string _description;
	string _locus_tag;
	string _ref_tag;
	
    long double _sample_mass_fraction;
    long double _sample_mass_variance;
};

class AbundanceGroup : public Abundance
{
public:
	AbundanceGroup() : _kappa(1.0), _FPKM_variance(0.0), _max_mass_variance(0.0), _salient_frags(0.0), _total_frags(0.0) {}
	
	AbundanceGroup(const AbundanceGroup& other) 
	{
		_abundances = other._abundances;
		_iterated_exp_count_covariance = other._iterated_exp_count_covariance;
        _count_covariance = other._count_covariance;
        _fpkm_covariance = other._fpkm_covariance;
        _gamma_covariance = other._gamma_covariance;
        _gamma_bootstrap_covariance = other._gamma_bootstrap_covariance;
		_FPKM_conf = other._FPKM_conf;
		_kappa = other._kappa;
		_kappa_covariance = other._kappa_covariance;
		_FPKM_variance = other._FPKM_variance;
		_description = other._description;
        _max_mass_variance = other._max_mass_variance;
        _salient_frags = other._salient_frags;
        _total_frags = other._total_frags;
        _read_group_props = other._read_group_props;
    }
	
	AbundanceGroup(const vector<shared_ptr<Abundance> >& abundances) : 
		_abundances(abundances), 
        _iterated_exp_count_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
        _count_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
        _fpkm_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
		_gamma_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
        _gamma_bootstrap_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
		_kappa_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())),
		_kappa(1.0),
		_FPKM_variance(0.0), 
        _max_mass_variance(0.0),
        _salient_frags(0.0),
        _total_frags(0.0) {}
    
	AbundanceGroup(const vector<shared_ptr<Abundance> >& abundances,
				   const ublas::matrix<double>& gamma_covariance,
                   const ublas::matrix<double>& gamma_bootstrap_covariance,
                   const ublas::matrix<double>& iterated_exp_count_covariance,
                   const ublas::matrix<double>& count_covariance,
                   const ublas::matrix<double>& fpkm_covariance,
                   const long double max_mass_variance,
                   const std::set<shared_ptr<ReadGroupProperties const > >& rg_props); 
	
	AbundanceStatus status() const;
	void status(AbundanceStatus s)			{ }
	bool has_member_with_status(AbundanceStatus member_status);
    
	double FPKM() const;
    void   FPKM(double fpkm)                { }
    
	double FPKM_variance() const			{ return _FPKM_variance; }
	void   FPKM_variance(double v)			{ }
	
	ConfidenceInterval	FPKM_conf() const   { return _FPKM_conf; }

	
	double gamma() const;
	void gamma(double g)					{  };
	
	double kappa() const					{ return _kappa; }
	void kappa(double k)					{ _kappa = k; }
	
	double num_fragments() const;
	void num_fragments(double nf)			{  }
    
    double mass_fraction() const;
	void mass_fraction(double mf)			{  }
    
    double mass_variance() const;
	void mass_variance(double mf)	{  }
	
	set<string> gene_id() const;	
	set<string> gene_name() const;
	set<string> tss_id() const;
	set<string> protein_id() const;
	
	virtual const string&	description() const	{ return _description; }
	virtual void			description(const string& d) { _description = d; }
	
	virtual const string&	locus_tag() const;
	virtual void			locus_tag(const string& L) { }
	
	virtual const string&	reference_tag() const;
	virtual void			reference_tag(const string& r) { } 
	
	double effective_length() const;
	
	//DUMMY FUNCTIONS
	void effective_length(double ef) {}
	void cond_probs(vector<double>* cp) {}


	void filter_group(const vector<bool>& to_keep, 
					  AbundanceGroup& filtered_group) const;
	
	void get_transfrags(vector<shared_ptr<Abundance> >& transfrags) const;
												 
	vector<shared_ptr<Abundance> >& abundances() { return _abundances; }
	const vector<shared_ptr<Abundance> >& abundances() const { return _abundances; }
	
	const ublas::matrix<double>& gamma_cov() const { return _gamma_covariance; }
    
    const ublas::matrix<double>& gamma_bootstrap_cov() const { return _gamma_bootstrap_covariance; }
    
    const ublas::matrix<double>& iterated_count_cov() const { return _iterated_exp_count_covariance; }
    
    const ublas::matrix<double>& count_cov() const { return _count_covariance; }
    
	const ublas::matrix<double>& kappa_cov() const { return _kappa_covariance; }
    
    const ublas::matrix<double>& fpkm_cov() const { return _kappa_covariance; }
	
	
	void calculate_abundance(const vector<MateHit>& alignments);
	
    void max_mass_variance(double mmv) { _max_mass_variance = mmv; }
    double max_mass_variance() const { return _max_mass_variance; }
    
    double salient_frags() const { return _salient_frags; }
    void salient_frags(double nf) { _salient_frags = nf; }
    
    double total_frags() const { return _total_frags; }
    void total_frags(double nf) { _total_frags = nf; }
    
    const std::set<shared_ptr<ReadGroupProperties const> >& rg_props() const { return _read_group_props; }
    
private:
	
	void FPKM_conf(const ConfidenceInterval& cf)  { _FPKM_conf = cf; }
	
	bool calculate_gammas(const vector<MateHit>& nr_alignments, 
                          const vector<double>& log_conv_factors,
						  const vector<shared_ptr<Abundance> >& transcripts,
						  const vector<shared_ptr<Abundance> >& mapped_transcripts);
	void calculate_FPKM_covariance();
    void estimate_count_covariance();
    void calculate_conf_intervals();
	void calculate_locus_scaled_mass_and_variance(const vector<MateHit>& nr_alignments,
                          const vector<shared_ptr<Abundance> >& transcripts);
    void calculate_iterated_exp_count_covariance(const vector<MateHit>& nr_alignments, 
                                                 const vector<shared_ptr<Abundance> >& transcripts);
	void calculate_kappas();
    
    
	void update_multi_reads(const vector<MateHit>& alignments, vector<shared_ptr<Abundance> > transcripts);

    
	void compute_cond_probs_and_effective_lengths(const vector<MateHit>& alignments, 
												  vector<shared_ptr<Abundance> >& transcripts,
												  vector<shared_ptr<Abundance> >& mapped_transcripts);
    
    void update_transcript_expression(double locus_mass, double locus_mass_fraction);
    
    
    
    //void collect_read_group_props();
	
	vector<shared_ptr<Abundance> > _abundances;
	
    // _count_covariance is the final count covariance matrix.  It's includes our estimates
    // of transcript-level biological variability on counts
    ublas::matrix<double> _count_covariance;
    
    // _iterated_exp_count_covariance is the ITERATED EXPECTATION count covariance matrix.  It's not the 
    // estimated count covariance matrix (i.e. it doesn't include biological variability from
    // the fitted model.
    ublas::matrix<double> _iterated_exp_count_covariance;
    ublas::matrix<double> _fpkm_covariance;
	ublas::matrix<double> _gamma_covariance;
    ublas::matrix<double> _gamma_bootstrap_covariance;
    
	ConfidenceInterval _FPKM_conf;
	
	ublas::matrix<double> _kappa_covariance;
	double _kappa;
	double _FPKM_variance;
	string _description;
    double _max_mass_variance;  // upper bound on the count variance that could come from this group.
    double _salient_frags;
    double _total_frags;
    
    std::set<shared_ptr<ReadGroupProperties const > > _read_group_props;
    //std::map<shared_ptr<ReadGroupProperties const >, ublas::vector<double> > _mles_for_read_groups;
};

void compute_compatibilities(vector<shared_ptr<Abundance> >& transcripts,
							 const vector<MateHit>& alignments,
							 vector<vector<char> >& compatibilities);

void get_alignments_from_scaffolds(const vector<shared_ptr<Abundance> >& abundances,
								   vector<MateHit>& alignments);

AbundanceStatus empirical_mean_replicate_gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                                                   const vector<MateHit>& nr_alignments,
                                                   const vector<double>& log_conv_factors,
                                                   ublas::vector<double>& gamma_map_estimate,
                                                   ublas::matrix<double>& gamma_covariance,
                                                   std::map<shared_ptr<ReadGroupProperties const >, ublas::vector<double> >& mles_for_read_groups);

AbundanceStatus empirical_replicate_gammas(const vector<shared_ptr<Abundance> >& transcripts,
                                           const vector<MateHit>& nr_alignments,
                                           const vector<double>& log_conv_factors,
                                           ublas::vector<double>& gamma_map_estimate,
                                           ublas::matrix<double>& gamma_map_covariance,
                                           std::map<shared_ptr<ReadGroupProperties const >, ublas::vector<double> >& mles_for_read_groups);

AbundanceStatus bootstrap_gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                                    const vector<MateHit>& nr_alignments,
                                    const vector<double>& log_conv_factors,
                                    ublas::vector<double>& gamma_map_estimate,
                                    ublas::matrix<double>& gamma_covariance,
                                    double& cross_replicate_js);

AbundanceStatus bootstrap_gammas(const vector<shared_ptr<Abundance> >& transcripts,
                                 const vector<MateHit>& nr_alignments,
                                 const vector<double>& log_conv_factors,
                                 ublas::vector<double>& gamma_map_estimate,
                                 ublas::matrix<double>& gamma_map_covariance,
                                 double& cross_replicate_js);

AbundanceStatus bayesian_gammas(const vector<shared_ptr<Abundance> >& transcripts,
                                 const vector<MateHit>& nr_alignments,
                                 const vector<double>& log_conv_factors,
                                 const ublas::vector<double>& gamma_mle,
                                 ublas::vector<double>& gamma_map_estimate,
                                 ublas::matrix<double>& gamma_map_covariance);

AbundanceStatus map_estimation(const vector<shared_ptr<Abundance> >& transcripts,
                               const vector<MateHit>& alignments,
                               const vector<double>& log_conv_factors,
                               const ublas::vector<double>&  proposal_gamma_mean,
                               const ublas::matrix<double>& proposal_gamma_covariance,
                               ublas::vector<double>& gamma_map_estimate,
                               ublas::matrix<double>& gamma_map_covariance);

AbundanceStatus gamma_mle(const vector<shared_ptr<Abundance> >& transcripts,
                          const vector<MateHit>& nr_alignments,
                          const vector<double>& log_conv_factors,
                          vector<double>& gammas,
                          bool check_identifiability = true,
                          vector<double>* p_hint = NULL);

double compute_doc(int bundle_origin, 
				   const vector<Scaffold>& scaffolds,
				   vector<float>& depth_of_coverage,
				   map<pair<int, int>, float>& intron_depth_of_coverage,
				   bool exclude_intra_intron=false,
				   vector<float>* intronic_cov=NULL,
				   vector<int>* scaff_intronic_status=NULL);

double major_isoform_intron_doc(map<pair<int, int>, float>& intron_doc);

void record_doc_for_scaffolds(int bundle_origin, 
							  const std::vector<Scaffold>& hits,
							  const std::vector<float>& depth_of_coverage,
							  const std::map<std::pair<int, int>, float>& intron_depth_of_coverage,
							  std::vector<double>& scaff_doc);

void record_doc_for_scaffolds(int bundle_origin, 
							  const std::vector<Scaffold>& hits,
							  const std::vector<float>& depth_of_coverage,
							  std::vector<double>& scaff_doc);

void record_min_doc_for_scaffolds(int bundle_origin, 
								  const std::vector<Scaffold>& hits,
								  const std::vector<float>& depth_of_coverage,
								  const std::map<std::pair<int, int>, float>& intron_depth_of_coverage,
								  std::vector<double>& scaff_doc);


double get_intron_doc(const Scaffold& s,
					  const map<pair<int, int>, float>& intron_depth_of_coverage);

double get_scaffold_doc(int bundle_origin, 
						const Scaffold& s,
						const vector<float>& depth_of_coverage);

double get_scaffold_min_doc(int bundle_origin, 
							const Scaffold& s,
							const vector<float>& depth_of_coverage);

AbundanceStatus calculate_inverse_fisher(const vector<shared_ptr<Abundance> >& transcripts,
                                         const vector<MateHit>& alignments,
                                         const ublas::vector<double>& gamma_mean,
                                         ublas::matrix<double>& inverse_fisher);
#endif
