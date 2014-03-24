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

#include <Eigen/Dense>

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
    
private:
    friend std::ostream & operator<<(std::ostream &os, const ConfidenceInterval &gp);
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){
        ar & low & high;
    }
};

enum AbundanceStatus { NUMERIC_OK, NUMERIC_FAIL, NUMERIC_LOW_DATA, NUMERIC_HI_DATA };

typedef map<boost::shared_ptr<ReadGroupProperties const>, double> CountPerReplicateTable;
typedef map<boost::shared_ptr<ReadGroupProperties const>, double> FPKMPerReplicateTable;
typedef map<boost::shared_ptr<ReadGroupProperties const>, AbundanceStatus> StatusPerReplicateTable;

bool fit_negbin_dist(const vector<double> samples, double& r, double& p);
long double negbin_log_likelihood(const vector<double>& samples, long double r, long double p);
long double poisson_log_likelihood(const vector<double>& samples, long double lambda);

namespace bser = boost::serialization;

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
	
    // This tracks the final modeled variance in the assigned counts.
	virtual double			num_fragments() const = 0;
	virtual void			num_fragments(double nf) = 0;
    
    // This tracks the final modeled variance in the assigned counts.
    virtual double num_fragment_var() const	= 0;
	virtual void num_fragment_var(double nfv) = 0;
    
    virtual CountPerReplicateTable num_fragments_by_replicate() const = 0;
    virtual void            num_fragments_by_replicate(CountPerReplicateTable& cpr) = 0;
    
    virtual FPKMPerReplicateTable FPKM_by_replicate() const = 0;
    virtual void            FPKM_by_replicate(CountPerReplicateTable& cpr) = 0;
    
    virtual StatusPerReplicateTable status_by_replicate() const = 0;
    virtual void status_by_replicate(StatusPerReplicateTable& fpr) = 0;
    
    // This tracks the fitted variance from the overdispersion model,
    // and does not include the fragment assignment uncertainty.
    virtual double          mass_variance() const = 0;
	virtual void            mass_variance(double mv) = 0;
    
    // This tracks the fragment assignment uncertainty,
    // and does not include the overdispersion.
    virtual double          num_fragment_uncertainty_var() const = 0;
	virtual void            num_fragment_uncertainty_var(double mv) = 0;
    
	virtual double			effective_length() const= 0;
	virtual void			effective_length(double el) = 0;
    
	virtual const vector<double>*	cond_probs() const		{ return NULL; }
	virtual void					cond_probs(vector<double>* cp) = 0;
	
	// The structural information for the object, if defined.
	virtual boost::shared_ptr<Scaffold> transfrag() const		{ return boost::shared_ptr<Scaffold>(); }
	
    virtual const vector<double>& fpkm_samples() const = 0;
    virtual void  fpkm_samples(const vector<double>& s) = 0;
    
	
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
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        
    }
    
    virtual void clear_non_serialized_data() = 0;
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Abundance);

class TranscriptAbundance : public Abundance
{
public:
	
	TranscriptAbundance() : 
		_status(NUMERIC_OK), 
		_transfrag(boost::shared_ptr<Scaffold>()), 
		_FPKM(0), 
		_FPKM_variance(0),
		_gamma(0), 
		_kappa(1.0), 
		_num_fragments(0),
        _num_fragment_var(0),
        _num_fragment_uncertainty_var(0),
		_eff_len(0),
		_cond_probs(NULL),
        _sample_mass_variance(0.0){}
	
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
        if (_transfrag)
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
	
    // This returns the estimated number of fragments
	double num_fragments() const			{ return _num_fragments; }
	void num_fragments(double nf)
    {
        assert (!isnan(nf));
        _num_fragments = nf;
        if (_transfrag)
            _transfrag->num_fragments(nf);
    }
    
    // This tracks the final modeled variance in the assigned counts.
    double num_fragment_var() const			{ return _num_fragment_var; }
	void num_fragment_var(double nfv)		{ assert (!isnan(nfv)); _num_fragment_var = nfv; }
    
    // This tracks the fragment assignment uncertainty,
    // and does not include the overdispersion.
    double          num_fragment_uncertainty_var() const { return _num_fragment_uncertainty_var; }
	void            num_fragment_uncertainty_var(double nfv) { assert (!isnan(nfv)); _num_fragment_uncertainty_var = nfv; }
    
    CountPerReplicateTable num_fragments_by_replicate() const { return _num_fragments_per_replicate; }
    void num_fragments_by_replicate(CountPerReplicateTable& cpr) { _num_fragments_per_replicate = cpr; }
    
    FPKMPerReplicateTable FPKM_by_replicate() const { return _fpkm_per_replicate; }
    void FPKM_by_replicate(FPKMPerReplicateTable& fpr) { _fpkm_per_replicate = fpr; }
    
    StatusPerReplicateTable status_by_replicate() const { return _status_per_replicate; }
    void status_by_replicate(StatusPerReplicateTable& fpr) { _status_per_replicate = fpr; }
    
    double mass_variance() const			{ return _sample_mass_variance; }
	void mass_variance(double mv)			{ _sample_mass_variance = mv; }
	
	void transfrag(boost::shared_ptr<Scaffold> tf)		{ _transfrag = tf; }
	boost::shared_ptr<Scaffold> transfrag() const		{ return _transfrag; }
	
	double effective_length() const			{ return _eff_len; }
	void effective_length(double el)		{ _eff_len = el; }
	
	const vector<double>* cond_probs() const	{ return _cond_probs; }
	void cond_probs(vector<double>* cp) 	
	{ 
		if(_cond_probs != NULL) { delete _cond_probs; };
		_cond_probs = cp;
	}
	
    const vector<double>& fpkm_samples() const { return _fpkm_samples; }
    void  fpkm_samples(const vector<double>& s) { _fpkm_samples = s; }
	
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
    
	void clear_non_serialized_data();
private:
	
    friend std::ostream & operator<<(std::ostream &os, const TranscriptAbundance &gp);
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){
        
        ar & _status;
        //ar & _transfrag;
        ar & _FPKM;
        ar & _FPKM_conf;
        ar & _gamma;
        ar & _kappa;
        ar & _num_fragments;
        ar & _num_fragment_var;
        ar & _num_fragment_uncertainty_var;
        ar & _eff_len;
        //ar & _cond_probs;
        ar & _description;
        ar & _locus_tag;
        ar & _ref_tag;
        ar & _sample_mass_variance;
        ar & _num_fragments_per_replicate;
        ar & _fpkm_per_replicate;
        ar & _status_per_replicate;
    }
    
	void calculate_FPKM_err_bar(double variance);
	
	AbundanceStatus _status;
	boost::shared_ptr<Scaffold> _transfrag;
	double _FPKM;
	double _FPKM_variance;
	ConfidenceInterval _FPKM_conf;
	double _gamma;
	double _kappa;
	double _num_fragments;
    double _num_fragment_var;
    double _num_fragment_uncertainty_var;
	double _eff_len;
	vector<double>* _cond_probs;
    
    vector<double> _fpkm_samples;
	
	string _description;
	string _locus_tag;
	string _ref_tag;
	
    long double _sample_mass_variance;
    
    CountPerReplicateTable _num_fragments_per_replicate;
    FPKMPerReplicateTable _fpkm_per_replicate;
    StatusPerReplicateTable _status_per_replicate;
};

//BOOST_CLASS_EXPORT_GUID(TranscriptAbundance, "TranscriptAbundance");
//BOOST_SERIALIZATION_boost::shared_ptr(TranscriptAbundance)

class AbundanceGroup : public Abundance
{
public:
	AbundanceGroup() :
        _kappa(1.0),
        _FPKM_variance(0.0),
        _salient_frags(0.0),
        _total_frags(0.0) {}
	
	AbundanceGroup(const vector<boost::shared_ptr<Abundance> >& abundances) :
		_abundances(abundances), 
        _iterated_exp_count_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
        _count_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
        _fpkm_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
		_gamma_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
		_kappa_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())),
		_kappa(1.0),
		_FPKM_variance(0.0),
        _salient_frags(0.0),
        _total_frags(0.0) {}
    
	AbundanceGroup(const vector<boost::shared_ptr<Abundance> >& abundances,
				   const ublas::matrix<double>& gamma_covariance,
                   const ublas::matrix<double>& iterated_exp_count_covariance,
                   const ublas::matrix<double>& count_covariance,
                   const ublas::matrix<double>& fpkm_covariance,
                   const std::set<boost::shared_ptr<ReadGroupProperties const > >& rg_props);
	
	AbundanceStatus status() const;
	void status(AbundanceStatus s)			{ }
	bool has_member_with_status(AbundanceStatus member_status) const;
    
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
    
    // This tracks the final modeled variance in the assigned counts.
    double num_fragment_var() const;
	void num_fragment_var(double nfv)		{  }
    
    // This tracks the uncertainty in the assigned counts
    double num_fragment_uncertainty_var() const;
	void num_fragment_uncertainty_var(double nfv)		{  }
    
    CountPerReplicateTable num_fragments_by_replicate() const;
    void num_fragments_by_replicate(CountPerReplicateTable& cpr) { }
    
    FPKMPerReplicateTable FPKM_by_replicate() const;
    void FPKM_by_replicate(FPKMPerReplicateTable& fpr) { }
    
    StatusPerReplicateTable status_by_replicate() const;
    void status_by_replicate(StatusPerReplicateTable& fpr) { }
    
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
	
	void get_transfrags(vector<boost::shared_ptr<Abundance> >& transfrags) const;
												 
	vector<boost::shared_ptr<Abundance> >& abundances() { return _abundances; }
	const vector<boost::shared_ptr<Abundance> >& abundances() const { return _abundances; }
	
	const ublas::matrix<double>& gamma_cov() const { return _gamma_covariance; }
    
    const ublas::matrix<double>& iterated_count_cov() const { return _iterated_exp_count_covariance; }
    
    const ublas::matrix<double>& count_cov() const { return _count_covariance; }
    
	const ublas::matrix<double>& kappa_cov() const { return _kappa_covariance; }
    
    const ublas::matrix<double>& fpkm_cov() const { return _fpkm_covariance; }
	
	const vector<Eigen::VectorXd>& assigned_counts() const { return _assigned_count_samples; }
    
    const vector<double>& fpkm_samples() const { return _fpkm_samples; }
    void  fpkm_samples(const vector<double>& s) { _fpkm_samples = s; }
    void clear_fpkm_samples() {
        _fpkm_samples.clear();
        std::vector<double>().swap(_fpkm_samples);
        for (size_t i = 0; i < _member_fpkm_samples.size(); ++i)
        {
            _member_fpkm_samples.clear();
            std::vector<Eigen::VectorXd>().swap(_member_fpkm_samples);
        }
    }
    
    const vector<Eigen::VectorXd>& member_fpkm_samples() const { return _member_fpkm_samples; }
    void  member_fpkm_samples(const vector<Eigen::VectorXd>& s) { _member_fpkm_samples = s; }
    
	void calculate_abundance(const vector<MateHit>& alignments,
                             bool perform_collapse = true,
                             bool calculate_variance = true);
	
    double salient_frags() const { return _salient_frags; }
    void salient_frags(double nf) { _salient_frags = nf; }
    
    double total_frags() const { return _total_frags; }
    void total_frags(double nf) { _total_frags = nf; }
    
    const std::set<boost::shared_ptr<ReadGroupProperties const> >& rg_props() const { return _read_group_props; }
    
    void init_rg_props(const std::set<boost::shared_ptr<ReadGroupProperties const> >& rg) 
    {
        _read_group_props = rg;
        _count_per_replicate.clear();
        for ( std::set<boost::shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg.begin();
             itr != rg.end();
             ++itr)
        {
            _count_per_replicate[*itr] = 0;
        }
    }
    
    void fit_gamma_distributions();
    
    void clear_non_serialized_data();
    
    void aggregate_replicate_abundances(const std::map<boost::shared_ptr<ReadGroupProperties const >, boost::shared_ptr<const AbundanceGroup> >& ab_group_per_replicate);
    
    void calculate_abundance_group_variance(const vector<boost::shared_ptr<Abundance> >& transcripts,
                                            const std::map<boost::shared_ptr<ReadGroupProperties const >, boost::shared_ptr<const AbundanceGroup> >& ab_group_per_replicate);

    void collect_per_replicate_mass(std::map<boost::shared_ptr<ReadGroupProperties const >, boost::shared_ptr<const AbundanceGroup> >& ab_group_per_replicate)
    {
        _count_per_replicate.clear();
        for (std::map<boost::shared_ptr<ReadGroupProperties const >, boost::shared_ptr<const AbundanceGroup> >::const_iterator itr = ab_group_per_replicate.begin();
             itr != ab_group_per_replicate.end(); ++itr)
        {
            _count_per_replicate[itr->first] = itr->second->num_fragments();
        }
    }
    
    static void apply_normalization_to_abundances(const map<boost::shared_ptr<const ReadGroupProperties>, boost::shared_ptr<const AbundanceGroup> >& unnormalized_ab_group_per_replicate,
                                                  map<boost::shared_ptr<const ReadGroupProperties>, boost::shared_ptr<AbundanceGroup> >& normalized_ab_group_per_replicate);
    
private:
    
    void calculate_abundance_for_replicate(const vector<MateHit>& alignments, bool perform_collapse);
    
    
	
	void FPKM_conf(const ConfidenceInterval& cf)  { _FPKM_conf = cf; }
	
	bool calculate_gammas(const vector<MateHit>& nr_alignments, 
                          const vector<double>& log_conv_factors,
						  const vector<boost::shared_ptr<Abundance> >& transcripts,
						  const vector<boost::shared_ptr<Abundance> >& mapped_transcripts);
	void calculate_FPKM_covariance();

    void calculate_conf_intervals();
	void calculate_locus_scaled_mass_and_variance(const vector<boost::shared_ptr<Abundance> >& transcripts);
   
    
    AbundanceStatus calculate_per_replicate_abundances(vector<boost::shared_ptr<Abundance> >& transcripts,
                                                       const std::map<boost::shared_ptr<ReadGroupProperties const >, vector<MateHit> >& alignments_per_read_group,
                                                       std::map<boost::shared_ptr<ReadGroupProperties const >, boost::shared_ptr<AbundanceGroup> >& ab_group_per_replicate,
                                                       bool perform_collapse = true);
    
        
	void calculate_kappas();
    
    
	void update_multi_reads(const vector<MateHit>& alignments, vector<boost::shared_ptr<Abundance> > transcripts);
    
    void collect_per_replicate_mass(const vector<MateHit>& alignments,
                                    vector<boost::shared_ptr<Abundance> >& transcripts);
    
    void generate_fpkm_samples();
    
    friend std::ostream & operator<<(std::ostream &os, const AbundanceGroup &gp);
    friend class boost::serialization::access;
    
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        //ar & _abundances;
        vector<TranscriptAbundance*> tmp;
        BOOST_FOREACH(boost::shared_ptr<Abundance> ab, _abundances)
        {
            tmp.push_back((TranscriptAbundance*)&(*ab));
        }
        ar & tmp;
        ar & _iterated_exp_count_covariance;
        ar & _count_covariance;
        ar & _fpkm_covariance;
        ar & _gamma_covariance;
        
        //ar & _FPKM_conf;
        //ar & _kappa_covariance;
        //ar & _assign_probs;
        
        ar & _kappa;
        ar & _FPKM_variance;
        ar & _description;
        //ar & _salient_frags;
        ar & _total_frags;
        //ar & _fpkm_samples; // don't save the samples
        ar & _read_group_props;
        //ar & _member_fpkm_samples // don't save the member samples either
        ar & _count_per_replicate;

    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {

        vector<TranscriptAbundance*> tmp;
        ar & tmp;
        BOOST_FOREACH(TranscriptAbundance* ab, tmp)
        {
            _abundances.push_back(boost::shared_ptr<Abundance>(ab));
        }
        
        //ar & _abundances;
        
        ar & _iterated_exp_count_covariance;
        ar & _count_covariance;
        ar & _fpkm_covariance;
        ar & _gamma_covariance;
        
        //ar & _FPKM_conf;
        //ar & _kappa_covariance;
        //ar & _assign_probs;
        
        ar & _kappa;
        ar & _FPKM_variance;
        ar & _description;
        //ar & _salient_frags;
        ar & _total_frags;
        //ar & _fpkm_samples; // don't save the samples
        ar & _read_group_props;
        //ar & _member_fpkm_samples // don't save the member samples either
        ar & _count_per_replicate;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    //void collect_read_group_props();
	
	vector<boost::shared_ptr<Abundance> > _abundances;
	
    // _iterated_exp_count_covariance is the ITERATED EXPECTATION count covariance matrix.  It's not the 
    // estimated count covariance matrix (i.e. it doesn't include biological variability from
    // the fitted model.
    ublas::matrix<double> _iterated_exp_count_covariance;
    
    // _count_covariance is the final count covariance matrix.  It's includes our estimates
    // of transcript-level biological variability on counts
    ublas::matrix<double> _count_covariance;
    
    ublas::matrix<double> _fpkm_covariance;
	ublas::matrix<double> _gamma_covariance;
    
	ConfidenceInterval _FPKM_conf;
	
	ublas::matrix<double> _kappa_covariance;
    ublas::matrix<double> _assign_probs;
    
	double _kappa;
	double _FPKM_variance;
	string _description;
    double _salient_frags;
    double _total_frags;
    
    vector<double> _fpkm_samples;
    vector<Eigen::VectorXd> _member_fpkm_samples;
    
    std::set<boost::shared_ptr<ReadGroupProperties const > > _read_group_props;
    vector<Eigen::VectorXd> _assigned_count_samples;
    
    map<boost::shared_ptr<ReadGroupProperties const>, double> _count_per_replicate;
    //std::map<boost::shared_ptr<ReadGroupProperties const >, ublas::vector<double> > _mles_for_read_groups;
};

struct SampleAbundances
{
    string locus_tag;
	AbundanceGroup transcripts;
	vector<AbundanceGroup> primary_transcripts;
	vector<AbundanceGroup> gene_primary_transcripts;
	vector<AbundanceGroup> cds;
	vector<AbundanceGroup> gene_cds;
	vector<AbundanceGroup> genes;
	double cluster_mass;
};

void compute_cond_probs_and_effective_lengths(const vector<MateHit>& alignments, 
                                              vector<boost::shared_ptr<Abundance> >& transcripts,
                                              vector<boost::shared_ptr<Abundance> >& mapped_transcripts);

void compute_compatibilities(const vector<boost::shared_ptr<Abundance> >& transcripts,
							 const vector<MateHit>& alignments,
							 vector<vector<char> >& compatibilities);

void get_alignments_from_scaffolds(const vector<boost::shared_ptr<Abundance> >& abundances,
								   vector<MateHit>& alignments);

AbundanceStatus empirical_replicate_gammas(vector<boost::shared_ptr<Abundance> >& transcripts,
                                           const vector<MateHit>& nr_alignments,
                                           const vector<double>& log_conv_factors,
                                           ublas::vector<double>& gamma_map_estimate,
                                           ublas::matrix<double>& gamma_map_covariance,
                                           std::map<boost::shared_ptr<ReadGroupProperties const >, ublas::vector<double> >& mles_for_read_groups);

AbundanceStatus gamma_mle(const vector<boost::shared_ptr<Abundance> >& transcripts,
                          const vector<MateHit>& nr_alignments,
                          const vector<double>& log_conv_factors,
                          vector<double>& gammas,
                          bool check_identifiability = true);

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

AbundanceStatus calculate_inverse_fisher(const vector<boost::shared_ptr<Abundance> >& transcripts,
                                         const vector<MateHit>& alignments,
                                         const ublas::vector<double>& gamma_mean,
                                         ublas::matrix<double>& inverse_fisher);

void calculate_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                const Eigen::MatrixXd& transcript_cond_probs,
                                const Eigen::VectorXd& proposed_gammas,
                                Eigen::MatrixXd& assignment_probs);

void calculate_average_assignment_probs(const Eigen::VectorXd& alignment_multiplicities, 
                                        const Eigen::MatrixXd& transcript_cond_probs,
                                        const Eigen::VectorXd& proposed_gammas,
                                        Eigen::MatrixXd& assignment_probs);

void calculate_iterated_exp_count_covariance(const vector<double>& gammas,
                                             const vector<MateHit>& nr_alignments,
                                             const vector<boost::shared_ptr<Abundance> >& transcripts,
                                             ublas::matrix<double>& count_covariance);

bool simulate_count_covariance(const vector<double>& num_fragments,
                               const vector<double>& frag_variances,
                               const ublas::matrix<double>& iterated_exp_count_covariances,
                               const vector<boost::shared_ptr<Abundance> >& transcripts,
                               ublas::matrix<double>& count_covariances,
                               vector<Eigen::VectorXd>& assigned_count_samples,
                               vector<ublas::vector<double> >* gamma_samples);

void sample_abundance_worker(const string& locus_tag,
                             const set<boost::shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAbundances& sample,
                             boost::shared_ptr<HitBundle> sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis,
                             bool calculate_variance);

void merge_precomputed_expression_worker(const string& locus_tag,
                                         const vector<boost::shared_ptr<PrecomputedExpressionBundleFactory> >& expression_factories,
                                         SampleAbundances& sample,
                                         boost::shared_ptr<HitBundle> sample_bundle,
                                         bool perform_cds_analysis,
                                         bool perform_tss_analysis,
                                         bool calculate_variance);
#endif
