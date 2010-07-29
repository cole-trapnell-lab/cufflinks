#ifndef BUNDLES_H
#define BUNDLES_H
/*
 *  bundles.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <boost/bind.hpp>
#include <vector>
#include "hits.h"
#include "scaffolds.h"
#include "GList.hh"
#include "gtf_tracking.h"

class EmpDist
{
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
		if (l >= _pdf.size() || l < 0)
			return 0;
		return _pdf[l];
	}
	
	void cdf(vector<double>& cdf)	{ _cdf = cdf; }
	double cdf(int l) const
	{
		if (l >= _cdf.size())
			return 1;
        if (l < 0)
            return 0;
		return _cdf[l];
	}
	
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

struct BundleStats
{		
	BundleStats() : 
	compatible(0), 
	uncollapsible(0),
	closure_edges(0),
	matched_edges(0)
	{}
	
	int compatible;
	int uncollapsible;
	int closure_edges;
	int matched_edges;
};

typedef map<RefID, vector<AugmentedCuffOp> > BadIntronTable;

/*******************************************************************************
 HitBundle is a set of MateHit objects that, were you to look at the interval
 graph of their spanning intervals in genomic coordinates, you'd see a single
 connected component. Note that bundles do not correspond to single transcripts,
 or even single genes
 *******************************************************************************/
class HitBundle
{
public:
	HitBundle() : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id) {}
	~HitBundle() 
	{
	}
	
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
	
	const std::vector<MateHit>& hits() const { return _hits; } 
	const std::vector<MateHit>& non_redundant_hits() const { return _non_redundant; } 
	const std::vector<double>& collapse_counts() const { return _collapse_counts; } 
	
	RefID ref_id()  const
	{
		if (!_hits.empty())
			return _hits.front().ref_id();
		else if (!_ref_scaffs.empty())
			return _ref_scaffs.front().ref_id();
		else
			return 0;
	}
	
	int id() const { return _id; }
	
	void add_ref_scaffold(const Scaffold& scaff)
	{
		_ref_scaffs.push_back(scaff);
	}
	
	const vector<Scaffold>& ref_scaffolds() const { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	void add_open_hit(shared_ptr<ReadHit> bh);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	// Sorts the hits, and performs other functions needed to prepare this 
	// bundle for assembly and quantitation
	void finalize();
	
	void remove_hitless_scaffolds();
	
private:
	int _leftmost;
	int _rightmost;
	std::vector<MateHit> _hits;
	std::vector<MateHit> _non_redundant;
	std::vector<double> _collapse_counts;
	std::vector<Scaffold> _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	OpenMates _open_mates;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<Scaffold>& ref_mRNAs,
				   bool loadSeqs=false,
				   bool loadFPKM=false);

class BundleFactory
{
public:
	BundleFactory(HitFactory& fac, FILE* ref_rna_file, FILE* maskf)
	: _hit_fac(fac), ref_mRNA_file(ref_rna_file), mask_file(maskf){}

	bool next_bundle(HitBundle& bundle_out);
	
	HitFactory& hit_factory() { return _hit_fac; } 
	
	void reset() 
	{ 
		//rewind(hit_file); 
		_hit_fac.reset();
		next_ref_scaff = ref_mRNAs.begin(); 
	}
	
	void load_ref_rnas(bool loadSeqs=false, bool loadFPKM=false) 
    {
        if (ref_mRNA_file)
        {
            ::load_ref_rnas(ref_mRNA_file, _hit_fac.ref_table(), ref_mRNAs, loadSeqs, loadFPKM);
            RefID last_id = 0;
            for (vector<Scaffold>::iterator i = ref_mRNAs.begin(); i < ref_mRNAs.end(); ++i)
            {
                if (i->ref_id() != last_id)
                {
                    _ref_scaff_offsets.push_back(make_pair(i->ref_id(), i));
                }
                last_id = i->ref_id();
            }
            
            next_ref_scaff = ref_mRNAs.begin();
        }
        
        if (mask_file)
        {
            ::load_ref_rnas(mask_file, _hit_fac.ref_table(), mask_gtf_recs, false, false);
            RefID last_id = 0;
            for (vector<Scaffold>::iterator i = mask_gtf_recs.begin(); i < mask_gtf_recs.end(); ++i)
            {
                if (i->ref_id() != last_id)
                {
                    _mask_scaff_offsets.push_back(make_pair(i->ref_id(), i));
                }
                last_id = i->ref_id();
            }
            
            next_mask_scaff = mask_gtf_recs.begin();
        }
    }
        
	void bad_intron_table(const BadIntronTable& bad_introns) 
	{ 
		_bad_introns = bad_introns;
	}
    
    void frag_len_dist(shared_ptr<const EmpDist> fld) 
	{ 
		_frag_len_dist = fld;
	}
    
    shared_ptr<const EmpDist> frag_len_dist() const 
	{ 
		return _frag_len_dist;
	}
	
	void map_mass(long double mm) 
	{
		_map_mass = mm;
	}
	
	long double map_mass() const
	{
		return _map_mass;
	}
	
	bool spans_bad_intron(const ReadHit& read);
	
private:
	HitFactory& _hit_fac;

	vector<Scaffold> ref_mRNAs;
	FILE* ref_mRNA_file;
	vector<pair<RefID, vector<Scaffold>::iterator> > _ref_scaff_offsets;
	vector<Scaffold>::iterator next_ref_scaff;
    
    vector<Scaffold> mask_gtf_recs;
	FILE* mask_file;
	vector<pair<RefID, vector<Scaffold>::iterator> > _mask_scaff_offsets;
	vector<Scaffold>::iterator next_mask_scaff;
	
	BadIntronTable _bad_introns;
    shared_ptr<const EmpDist> _frag_len_dist;
	long double _map_mass;
};

void inspect_map(BundleFactory& bundle_factory, 
				 long double& map_mass, 
				 BadIntronTable& bad_introns,
                 EmpDist& frag_len_dist);

#endif
