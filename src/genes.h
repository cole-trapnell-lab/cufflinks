#ifndef ISOFORM_H
#define ISOFORM_H

/*
 *  genes.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scaffolds.h"
#include "abundances.h"
#include "common.h"

extern int next_isoform_id; 

int get_next_isoform_id();

extern int next_gene_id; 

int get_next_gene_id();

extern int next_skipped_region_id; 

int get_next_skipped_region_id();

class Isoform
{
public:
	Isoform(const Scaffold& s,
			int gid,
			int tid,
			double FPKM = 0.0, 
			double eff_len = 0.0,
			double fraction = 0.0,
			ConfidenceInterval ci = ConfidenceInterval(),
			double cov = 0.0,
            double est_frag_count = 0.0,
			double fmi = 0.0,
			AbundanceStatus status = NUMERIC_FAIL,
			string ref_gene_id = "") :
		_scaffold(s),
		_FPKM(FPKM),
		_eff_len(eff_len),
		_fraction(fraction),
		_confidence(ci),
		_coverage(cov),
        _estimated_count(est_frag_count),
		_FMI(fmi),
		_status(status),
		//allele
		_paternal_FPKM(0.0),
		_maternal_FPKM(0.0),
		_paternal_eff_len(eff_len),
		_maternal_eff_len(eff_len),
		_paternal_fraction(0.0),
		_maternal_fraction(0.0),
		_paternal_confidence(ci),
		_maternal_confidence(ci),
		_paternal_estimated_count(0.0),
		_maternal_estimated_count(0.0),
		_paternal_FMI(0.0),
		_maternal_FMI(0.0),
		_paternal_status(status),
		_maternal_status(status),
		_allele_informative(false)
	{
		_id = get_next_isoform_id();
		
		char trans_id_str[256];
		if (_scaffold.annotated_trans_id() != "")
			strncpy(trans_id_str, _scaffold.annotated_trans_id().c_str(), 255);
		else if (gid == -1)
			sprintf(trans_id_str, "%s.%s.%d", user_label.c_str(), ref_gene_id.c_str(), tid);
		else
			sprintf(trans_id_str, "%s.%d.%d", user_label.c_str(), gid, tid);
		
		_trans_id = trans_id_str;
		
		char gene_id_str[256];
		if(gid == -1)
			strncpy(gene_id_str, ref_gene_id.c_str(), 255);
		else
			sprintf(gene_id_str, "%s.%d", user_label.c_str(), gid);
		_gene_id = gene_id_str;
	}
	//allele
	Isoform(const Scaffold& s,
			int gid,
			int tid,
			double paternal_FPKM = 0.0, 
			double maternal_FPKM = 0.0, 	
			double paternal_eff_len = 0.0,
			double maternal_eff_len = 0.0,
			double paternal_fraction = 0.0,
			double maternal_fraction = 0.0,	
			ConfidenceInterval paternal_ci = ConfidenceInterval(),
			ConfidenceInterval maternal_ci = ConfidenceInterval(),	
			double paternal_cov = 0.0,
			double maternal_cov = 0.0,	
            double paternal_est_frag_count = 0.0,
			double maternal_est_frag_count = 0.0,
			double paternal_fmi = 0.0,
			double maternal_fmi = 0.0,	
			AbundanceStatus paternal_status = NUMERIC_FAIL,
			AbundanceStatus maternal_status = NUMERIC_FAIL,	
			string ref_gene_id = "",
			bool allele_informative = false) :
		_scaffold(s),
		_paternal_FPKM(paternal_FPKM),
		_maternal_FPKM(maternal_FPKM),
		_paternal_eff_len(paternal_eff_len),
		_maternal_eff_len(maternal_eff_len),
		_paternal_fraction(paternal_fraction),
		_maternal_fraction(maternal_fraction),
		_paternal_confidence(paternal_ci),
		_maternal_confidence(maternal_ci),
		_paternal_coverage(paternal_cov),
		_maternal_coverage(maternal_cov),
        _paternal_estimated_count(paternal_est_frag_count),
		_maternal_estimated_count(maternal_est_frag_count),
		_paternal_FMI(paternal_fmi),
		_maternal_FMI(maternal_fmi),
		_paternal_status(paternal_status),
		_maternal_status(maternal_status),
		_FPKM(paternal_FPKM+maternal_FPKM),
		_eff_len(paternal_eff_len),
		_fraction(paternal_fraction+maternal_fraction),
		_confidence(paternal_ci),
		_estimated_count(0.0),
		_FMI(0.0),
		_status(paternal_status),
		_allele_informative(allele_informative)
		{
			_id = get_next_isoform_id();
			
			char trans_id_str[256];
			if (_scaffold.annotated_trans_id() != "")
				strncpy(trans_id_str, _scaffold.annotated_trans_id().c_str(), 255);
			else if (gid == -1)
				sprintf(trans_id_str, "%s.%s.%d", user_label.c_str(), ref_gene_id.c_str(), tid);
			else
				sprintf(trans_id_str, "%s.%d.%d", user_label.c_str(), gid, tid);
			
			_trans_id = trans_id_str;
			
			char gene_id_str[256];
			if(gid == -1)
				strncpy(gene_id_str, ref_gene_id.c_str(), 255);
			else
				sprintf(gene_id_str, "%s.%d", user_label.c_str(), gid);
			_gene_id = gene_id_str;
			//This is performed cause _scaffold member is not allele specific
			_scaffold.fpkm(paternal_FPKM+maternal_FPKM);
			_scaffold.paternal_fpkm(paternal_FPKM);
			_scaffold.maternal_fpkm(maternal_FPKM);
		}
	const Scaffold& scaffold() const { return _scaffold; }
	double FPKM() const { return _FPKM; } 
	void   FPKM(double fpkm) { _FPKM = fpkm; }

	//allele
	double paternal_FPKM() const { return _paternal_FPKM; } 
	double maternal_FPKM() const { return _maternal_FPKM; } 
	void   paternal_FPKM(double paternal_fpkm) { _paternal_FPKM = paternal_fpkm; }
	void   maternal_FPKM(double maternal_fpkm) { _maternal_FPKM = maternal_fpkm; }	
	double effective_length() const { return _eff_len; } 
	//allele
	double paternal_effective_length() const { return _paternal_eff_len; } 
	double maternal_effective_length() const { return _maternal_eff_len; }
	void   effective_length(double eff_len) { _eff_len = eff_len; }
	//allele
	void   paternal_effective_length(double paternal_eff_len) { _paternal_eff_len = paternal_eff_len; }
	void   maternal_effective_length(double maternal_eff_len) { _maternal_eff_len = maternal_eff_len; }		
	AbundanceStatus status() const { return _status; } 
	//allele
	AbundanceStatus paternal_status() const { return _paternal_status; } 
	AbundanceStatus maternal_status() const { return _maternal_status; } 	
	void   status(AbundanceStatus status) { _status = status; }
	//allele
	void   paternal_status(AbundanceStatus paternal_status) { _paternal_status = paternal_status; }
	void   maternal_status(AbundanceStatus maternal_status) { _maternal_status = maternal_status; }	
	double fraction() const {return _fraction; }
	//allele
	double paternal_fraction() const {return _paternal_fraction; }
	double maternal_fraction() const {return _maternal_fraction; }
	void fraction(double f) { _fraction = f; }
	//allele
	void paternal_fraction(double paternal_f) { _paternal_fraction = paternal_f; }
	void maternal_fraction(double maternal_f) { _maternal_fraction = maternal_f; }	
	ConfidenceInterval confidence() const { return _confidence; }
	//allele
	ConfidenceInterval paternal_confidence() const { return _paternal_confidence; }
	ConfidenceInterval maternal_confidence() const { return _maternal_confidence; }
	void   confidence(ConfidenceInterval c) { _confidence = c; }
	//allele
	void   paternal_confidence(ConfidenceInterval paternal_c) { _paternal_confidence = paternal_c; }
	void   maternal_confidence(ConfidenceInterval maternal_c) { _maternal_confidence = maternal_c; }	
	double coverage() const { return _coverage; }
	//allele
	double paternal_coverage() const { return _paternal_coverage; }
	double maternal_coverage() const { return _maternal_coverage; }
	void   coverage(double cov) { _coverage = cov; }
	//allele
	void   paternal_coverage(double paternal_cov) { _paternal_coverage = paternal_cov; }
	void   maternal_coverage(double maternal_cov) { _maternal_coverage = maternal_cov; }	
	// fraction of major isoform expression
	double FMI() const { return _FMI; }
	//allele
	double paternal_FMI() const { return _paternal_FMI; }
	double maternal_FMI() const { return _maternal_FMI; }	
	void   FMI(double fmi) { _FMI = fmi; }
	//allele
	void   paternal_FMI(double paternal_fmi) { _paternal_FMI = paternal_fmi; }
	void   maternal_FMI(double maternal_fmi) { _maternal_FMI = maternal_fmi; }		
	int ID() const { return _id; }

	void get_gtf(vector<string>& gtf_recs, 
				 const RefSequenceTable& rt,
                 set<AugmentedCuffOp>* hit_introns=NULL) const;
	//allele
	void get_allele_gtf(vector<string>& gtf_recs, 
						const RefSequenceTable& rt,
						set<AugmentedCuffOp>* hit_introns=NULL) const;	
	void gene_id(string& gid) { _gene_id = gid; }
	const string& gene_id() const { return _gene_id; }
	const string& trans_id() const {return _trans_id; }
	
	bool is_ref_trans() const { return _scaffold.is_ref(); }
	
    double estimated_count() const { return _estimated_count; }
	//allele
	double paternal_estimated_count() const { return _paternal_estimated_count; }
	double maternal_estimated_count() const { return _maternal_estimated_count; }
    void estimated_count(double est) { _estimated_count = est; }
	//allele
	void paternal_estimated_count(double paternal_est) { _paternal_estimated_count = paternal_est; }
	void maternal_estimated_count(double maternal_est) { _maternal_estimated_count = maternal_est; }
	bool get_allele_informativeness () const { return _allele_informative; }
	void set_allele_informativeness (bool allele_informative)  {_allele_informative = allele_informative; }
private:
	
	Scaffold _scaffold;
	double _FPKM;
	//allele
	double _paternal_FPKM;
	double _maternal_FPKM;
	double _eff_len;
	//allele
	double _paternal_eff_len;
	double _maternal_eff_len;
	double _fraction;
	//allele
	double _paternal_fraction;
	double _maternal_fraction;
	ConfidenceInterval _confidence;
	//allele
	ConfidenceInterval _paternal_confidence;
	ConfidenceInterval _maternal_confidence;
	double _coverage;
	//allele
	double _paternal_coverage;
	double _maternal_coverage;
    double _estimated_count;
	//allele
	double _paternal_estimated_count;
	double _maternal_estimated_count;
	double _FMI;
	//allele
	double _paternal_FMI;
	double _maternal_FMI;
	int _id;
	string _gene_id;
	string _trans_id;
	AbundanceStatus _status;
	//allele
	AbundanceStatus _paternal_status;
	AbundanceStatus _maternal_status;
	bool _allele_informative;
};

class Gene
{
public:
	Gene(const vector<Isoform>& isoforms, 
		 double FPKM = 0.0,
		 const ConfidenceInterval& ci = ConfidenceInterval(),
		 AbundanceStatus status=NUMERIC_FAIL) : 
		_isoforms(isoforms), 
		_FPKM(FPKM),
		_confidence(ci),
		_status(status),
		//allele
		_paternal_FPKM(0.0),
		_maternal_FPKM(0.0),
		_paternal_confidence(ci),
		_maternal_confidence(ci),
		_paternal_status(status),
		_maternal_status(status)
	{		
		vector<Scaffold> scaffolds;
		for (size_t i = 0; i < isoforms.size(); ++i)
			scaffolds.push_back(isoforms[i].scaffold());
		
		// Now compute FPKM for the whole gene
		Scaffold smashed_gene;
		Scaffold::merge(scaffolds, smashed_gene, false);
		_left = smashed_gene.left();
		_right = smashed_gene.right();
		
		_gene_id = _isoforms.front().gene_id();
	}
	//allele
	Gene(const vector<Isoform>& isoforms, 
		 double paternal_FPKM = 0.0,
		 double maternal_FPKM = 0.0,
		 const ConfidenceInterval& paternal_ci = ConfidenceInterval(),
		 const ConfidenceInterval& maternal_ci = ConfidenceInterval(),
		 AbundanceStatus paternal_status=NUMERIC_FAIL,
		 AbundanceStatus maternal_status=NUMERIC_FAIL) : 
		_isoforms(isoforms), 
		_paternal_FPKM(paternal_FPKM),
		_maternal_FPKM(maternal_FPKM),
		_paternal_confidence(paternal_ci),
		_maternal_confidence(maternal_ci),
		_paternal_status(paternal_status),
		_maternal_status(maternal_status),
		_FPKM(paternal_FPKM+maternal_FPKM),
		_confidence(paternal_ci),
		_status(paternal_status)
		{		
			vector<Scaffold> scaffolds;
			for (size_t i = 0; i < isoforms.size(); ++i)
				scaffolds.push_back(isoforms[i].scaffold());
			
			// Now compute FPKM for the whole gene
			//Nimrod - since the whole gene FPKM isn't really used i didn't separate to an allele cases
			Scaffold smashed_gene;
			Scaffold::merge(scaffolds, smashed_gene, false);
			_left = smashed_gene.left();
			_right = smashed_gene.right();
			
			_gene_id = _isoforms.front().gene_id();
		}
	
	const vector<Isoform>& isoforms() const { return _isoforms; }
	double FPKM() const { return _FPKM; } 
	//allele
   	double paternal_FPKM() const { return _paternal_FPKM; } 
	double maternal_FPKM() const { return _maternal_FPKM; }	
	ConfidenceInterval confidence() const { return _confidence; }
	//allele
	ConfidenceInterval paternal_confidence() const { return _paternal_confidence; }
	ConfidenceInterval maternal_confidence() const { return _maternal_confidence; }
	void   confidence(ConfidenceInterval c) { _confidence = c; }
	//allele
	void   paternal_confidence(ConfidenceInterval paternal_c) { _paternal_confidence = paternal_c; }
	void   maternal_confidence(ConfidenceInterval maternal_c) { _maternal_confidence = maternal_c; }	
	AbundanceStatus status() const { return _status; }
	//allele
	AbundanceStatus paternal_status() const { return _paternal_status; }
	AbundanceStatus maternal_status() const { return _maternal_status; }
	void   status(AbundanceStatus status) { _status = status; }
	//allele
	void   paternal_status(AbundanceStatus paternal_status) { _paternal_status = paternal_status; }
	void   maternal_status(AbundanceStatus maternal_status) { _maternal_status = maternal_status; }	
	int left() const { return _left; }
	int right() const { return _right; }
	
	const string& gene_id() const { return _gene_id; }
	
	bool has_ref_trans() const
	{
		BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			if (iso.is_ref_trans())
				return true;
		}
		return false;
	}
    
    double estimated_count() const 
    {
        double est = 0.0;
        BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			est += iso.estimated_count(); 
		}
		return est;
    }
	//allele
	double paternal_estimated_count() const 
    {
        double est = 0.0;
		BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			est += iso.paternal_estimated_count(); 
		}
		return est;
    }
	
	double maternal_estimated_count() const 
    {
        double est = 0.0;
		BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			est += iso.maternal_estimated_count(); 
		}
		return est;
    }    
    double effective_length() const 
    {
        double eff = 0.0;
        double total_fpkm = 0;
        BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			eff += iso.FPKM() * iso.effective_length();
            total_fpkm += iso.FPKM();
		}
        if (total_fpkm)
            return eff / total_fpkm;
		else
            return 0;
    }
	//allele
	double paternal_effective_length() const 
    {
        double eff = 0.0;
        double total_fpkm = 0;
		BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			eff += iso.paternal_FPKM() * iso.paternal_effective_length();
            total_fpkm += iso.paternal_FPKM();
		}
        if (total_fpkm)
            return eff / total_fpkm;
		else
            return 0;
    }
	
	double maternal_effective_length() const 
    {
        double eff = 0.0;
        double total_fpkm = 0;
		BOOST_FOREACH (const Isoform& iso, _isoforms)
		{
			eff += iso.maternal_FPKM() * iso.maternal_effective_length();
            total_fpkm += iso.maternal_FPKM();
		}
        if (total_fpkm)
            return eff / total_fpkm;
		else
            return 0;
    }
	
	//allele
	//fills all isoform IDs that have at least n reads with allele information (does not map to any other isoform of this gene) 
	void set_allele_informative_isoforms(const int n = min_allele_reads);
	bool has_allele_informative_isoforms() const;	
private:
	
	vector<Isoform> _isoforms;
	double _FPKM;
	//allele
	double _paternal_FPKM;
	double _maternal_FPKM;
	ConfidenceInterval _confidence;
	//allele
	ConfidenceInterval _paternal_confidence;
	ConfidenceInterval _maternal_confidence;
	int _id;
	int _left;
	int _right;
	string _gene_id;
	AbundanceStatus _status;
	//allele
	AbundanceStatus _paternal_status;
	AbundanceStatus _maternal_status;
};

#endif
