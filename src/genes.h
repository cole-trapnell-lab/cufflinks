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
		_status(status)
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
	
	const Scaffold& scaffold() const { return _scaffold; }
	double FPKM() const { return _FPKM; } 
	void   FPKM(double fpkm) { _FPKM = fpkm; }
	
	double effective_length() const { return _eff_len; } 
	void   effective_length(double eff_len) { _eff_len = eff_len; }
	
	AbundanceStatus status() const { return _status; } 
	void   status(AbundanceStatus status) { _status = status; }
	
	double fraction() const {return _fraction; }
	void fraction(double f) { _fraction = f; }
	
	ConfidenceInterval confidence() const { return _confidence; }
	void   confidence(ConfidenceInterval c) { _confidence = c; }
	
	double coverage() const { return _coverage; }
	void   coverage(double cov) { _coverage = cov; }
	
	// fraction of major isoform expression
	double FMI() const { return _FMI; }
	void   FMI(double fmi) { _FMI = fmi; }
	
	int ID() const { return _id; }

	void get_gtf(vector<string>& gtf_recs, 
				 const RefSequenceTable& rt,
                 set<AugmentedCuffOp>* hit_introns=NULL) const;
	
	void gene_id(string& gid) { _gene_id = gid; }
	const string& gene_id() const { return _gene_id; }
	const string& trans_id() const {return _trans_id; }
	
	bool is_ref_trans() const { return _scaffold.is_ref(); }
	
    double estimated_count() const { return _estimated_count; }
    void estimated_count(double est) { _estimated_count = est; }
private:
	
	Scaffold _scaffold;
	double _FPKM;
	double _eff_len;
	double _fraction;
	ConfidenceInterval _confidence;
	double _coverage;
    double _estimated_count;
	double _FMI;
	int _id;
	string _gene_id;
	string _trans_id;
	AbundanceStatus _status;
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
		_status(status)
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
	
	const vector<Isoform>& isoforms() const { return _isoforms; }
	double FPKM() const { return _FPKM; } 
	
	ConfidenceInterval confidence() const { return _confidence; }
	void   confidence(ConfidenceInterval c) { _confidence = c; }
	
	AbundanceStatus status() const { return _status; }
	void   status(AbundanceStatus status) { _status = status; }
	
	int left() const { return _left; }
	int right() const { return _right; }
	
	const string& gene_id() const { return _gene_id; }
	
	bool has_ref_trans() const
	{
		foreach (const Isoform& iso, _isoforms)
		{
			if (iso.is_ref_trans())
				return true;
		}
		return false;
	}
    
    double estimated_count() const 
    {
        double est = 0.0;
        foreach (const Isoform& iso, _isoforms)
		{
			est += iso.estimated_count(); 
		}
		return est;
    }
    
    double effective_length() const 
    {
        double eff = 0.0;
        double total_fpkm = 0;
        foreach (const Isoform& iso, _isoforms)
		{
			eff += iso.FPKM() * iso.effective_length();
            total_fpkm += iso.FPKM();
		}
        if (total_fpkm)
            return eff / total_fpkm;
		else
            return 0;
    }
	
private:
	
	vector<Isoform> _isoforms;
	double _FPKM;
	ConfidenceInterval _confidence;
	int _id;
	int _left;
	int _right;
	string _gene_id;
	AbundanceStatus _status;
};

#endif
