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
	const std::vector<int>& collapse_counts() const { return _collapse_counts; } 
	
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
	std::vector<int> _collapse_counts;
	std::vector<Scaffold> _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	
	OpenMates _open_mates;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<Scaffold>& ref_mRNAs);

class BundleFactory
{
public:
	BundleFactory(HitFactory& fac, FILE* ref_rna_file)
	: _hit_fac(fac), ref_mRNA_file(ref_rna_file), _next_hit_num(0) {}

	bool next_bundle(HitBundle& bundle_out);
	
	HitFactory& hit_factory() { return _hit_fac; } 
	
	void reset() 
	{ 
		//rewind(hit_file); 
		_hit_fac.reset();
		
		_next_hit_num = 0;
		next_ref_scaff = ref_mRNAs.begin(); 
	}
	
	void load_ref_rnas() 
	{
		::load_ref_rnas(ref_mRNA_file, _hit_fac.ref_table(), ref_mRNAs);
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
	
	void bad_intron_table(const BadIntronTable& bad_introns) 
	{ 
		_bad_introns = bad_introns;
	}
	
	
	bool spans_bad_intron(const ReadHit& read);
	
private:
	HitFactory& _hit_fac;

	vector<Scaffold> ref_mRNAs;
	FILE* ref_mRNA_file;
	vector<pair<RefID, vector<Scaffold>::iterator> > _ref_scaff_offsets;
	vector<Scaffold>::iterator next_ref_scaff;
	int _next_hit_num;
	
	BadIntronTable _bad_introns;
};

void inspect_map(BundleFactory& bundle_factory, 
				 long double& map_mass, 
				 BadIntronTable& bad_introns);

#endif
