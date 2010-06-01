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

namespace ublas = boost::numeric::ublas;

struct ConfidenceInterval
{
	ConfidenceInterval(double Low = 0.0, double High = 0.0) 
	: low(Low), high(High) {}
	double low;
	double high;
};

enum AbundanceStatus { NUMERIC_OK, NUMERIC_FAIL };

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
	
	virtual double			effective_length() const = 0;
	
	// The structural information for the object, if defined.
	virtual const Scaffold* transfrag() const		{ return NULL; }
	
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
	
	virtual	void calculate_abundance(const vector<MateHit>& alignments,
									 long double total_map_mass) = 0;
};

class TranscriptAbundance : public Abundance
{
public:
	
	TranscriptAbundance() : 
		_status(NUMERIC_OK), 
		_transfrag(NULL), 
		_FPKM(0), 
		_FPKM_variance(0),
		_gamma(0), 
		_kappa(1.0), 
		_num_fragments(0) {}
	
	TranscriptAbundance(const TranscriptAbundance& other)
	{
		_status = other._status;
		_transfrag = other._transfrag;
		_FPKM = other._FPKM;
		_FPKM_conf = other._FPKM_conf;
		_gamma = other._gamma;
		_num_fragments = other._num_fragments;
	}
	
	AbundanceStatus status() const			{ return _status; }
	void status(AbundanceStatus s)			{ _status = s; }
	
	double FPKM() const						{ return _FPKM; }
	double FPKM_variance() const			{ return _FPKM_variance; }
	void   FPKM_variance(double v)			{ _FPKM_variance = v; }
	
	ConfidenceInterval	FPKM_conf() const   { return _FPKM_conf; }
	void FPKM_conf(const ConfidenceInterval& cf) { _FPKM_conf = cf; }
	
	double gamma() const					{ return _gamma; }
	void gamma(double g)					{ _gamma = g; };
	
	double kappa() const					{ return _kappa; }
	void kappa(double k)					{ _kappa = k; }
	
	double num_fragments() const			{ return _num_fragments; }
	void num_fragments(double nf)			{ _num_fragments = nf; }
	
	void transfrag(const Scaffold* tf)		{ _transfrag = tf; }
	const Scaffold* transfrag() const		{ return _transfrag; }
	
	double effective_length() const;
	
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
	
	void calculate_abundance(const vector<MateHit>& alignments,
							 long double total_map_mass);
	
private:
	
	void calculate_FPKM_err_bar(double variance);
	
	AbundanceStatus _status;
	const Scaffold* _transfrag;
	double _FPKM;
	double _FPKM_variance;
	ConfidenceInterval _FPKM_conf;
	double _gamma;
	double _kappa;
	double _num_fragments;
	
	string _description;
	string _locus_tag;
	string _ref_tag;
};

class AbundanceGroup : public Abundance
{
public:
	AbundanceGroup() : _kappa(1.0), _FPKM_variance(0.0), _sample_mass(0.0) {}
	
	AbundanceGroup(const AbundanceGroup& other) 
	{
		_abundances = other._abundances;
		_gamma_covariance = other._gamma_covariance;
		_FPKM_conf = other._FPKM_conf;
		_kappa = other._kappa;
		_kappa_covariance = other._kappa_covariance;
		_FPKM_variance = other._FPKM_variance;
		_description = other._description;
		_sample_mass = other._sample_mass;
	}
	
	AbundanceGroup(const vector<shared_ptr<Abundance> >& abundances) : 
		_abundances(abundances), 
		_gamma_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())), 
		_kappa_covariance(ublas::zero_matrix<double>(abundances.size(), abundances.size())),
		_kappa(1.0),
		_FPKM_variance(0.0),
		_sample_mass(0.0) {}
	
	AbundanceGroup(const vector<shared_ptr<Abundance> >& abundances,
				   const ublas::matrix<double>& gamma_covariance,
				   long double map_mass) :
		_abundances(abundances), 
		_gamma_covariance(gamma_covariance),
		_sample_mass(map_mass)
	{
		calculate_conf_intervals();
		calculate_kappas();
	}
	
	AbundanceStatus status() const;
	void status(AbundanceStatus s)			{ }
	
	double FPKM() const;
	double FPKM_variance() const			{ return _FPKM_variance; }
	void   FPKM_variance(double v)			{ }
	
	ConfidenceInterval	FPKM_conf() const   { return _FPKM_conf; }

	
	double gamma() const;
	void gamma(double g)					{  };
	
	double kappa() const					{ return _kappa; }
	void kappa(double k)					{ _kappa = k; }
	
	double num_fragments() const;
	void num_fragments(double nf)			{  }
	
	//const Scaffold* transfrag() const		{ return NULL; }
	
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
	
	void filter_group(const vector<bool>& to_keep, 
					  AbundanceGroup& filtered_group) const;
	
	void get_transfrags(vector<Scaffold>& scaffolds) const;
												 
	vector<shared_ptr<Abundance> >& abundances() { return _abundances; }
	const vector<shared_ptr<Abundance> >& abundances() const { return _abundances; }
	
	const ublas::matrix<double>& gamma_cov() const { return _gamma_covariance; }
	const ublas::matrix<double>& kappa_cov() const { return _kappa_covariance; }
	
	
	void calculate_abundance(const vector<MateHit>& alignments,
							 long double total_map_mass);
	
private:
	
	void FPKM_conf(const ConfidenceInterval& cf)  { _FPKM_conf = cf; }
	
	bool calculate_gammas(const vector<MateHit>& alignments);
	void calculate_FPKM_variance();
	void calculate_conf_intervals();
	void calculate_counts(const vector<MateHit>& alignments);
	void calculate_kappas();
	
	vector<shared_ptr<Abundance> > _abundances;
	
	ublas::matrix<double> _gamma_covariance;
	ConfidenceInterval _FPKM_conf;
	
	ublas::matrix<double> _kappa_covariance;
	double _kappa;
	double _FPKM_variance;
	string _description;
	double _sample_mass;
};

void get_alignments_from_scaffolds(const vector<shared_ptr<Abundance> >& abundances,
								   vector<MateHit>& alignments);

bool gamma_map(const vector<Scaffold>& transcripts,
			   const vector<MateHit>& alignments,
			   const vector<int>& collapse_counts,
			   vector<double>& gamma_map_estimate,
			   ublas::matrix<double>& gamma_covariance);

void gamma_mle(const vector<Scaffold>& transcripts,
			   const vector<MateHit>& alignments,
			   const vector<int>& collapse_counts,
			   vector<double>& gammas);

double compute_doc(int bundle_origin, 
				   const vector<Scaffold>& scaffolds,
				   vector<int>& depth_of_coverage,
				   map<pair<int, int>, int>& intron_depth_of_coverage,
				   bool exclude_intra_intron,
                   bool use_non_redundant = false);

double major_isoform_intron_doc(map<pair<int, int>, int>& intron_doc);

void record_doc_for_scaffolds(int bundle_origin, 
							  const std::vector<Scaffold>& hits,
							  const std::vector<int>& depth_of_coverage,
							  const std::map<std::pair<int, int>, int>& intron_depth_of_coverage,
							  std::vector<double>& scaff_doc);

void record_doc_for_scaffolds(int bundle_origin, 
							  const std::vector<Scaffold>& hits,
							  const std::vector<int>& depth_of_coverage,
							  std::vector<double>& scaff_doc);

void record_min_doc_for_scaffolds(int bundle_origin, 
								  const std::vector<Scaffold>& hits,
								  const std::vector<int>& depth_of_coverage,
								  const std::map<std::pair<int, int>, int>& intron_depth_of_coverage,
								  std::vector<double>& scaff_doc);


double get_intron_doc(const Scaffold& s,
					  const map<pair<int, int>, int >& intron_depth_of_coverage);

double get_scaffold_doc(int bundle_origin, 
						const Scaffold& s,
						const vector<int>& depth_of_coverage);

double get_scaffold_min_doc(int bundle_origin, 
							const Scaffold& s,
							const vector<int>& depth_of_coverage);

long double get_map_mass(BundleFactory& bundle_factory);

#endif
