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
#include <boost/random.hpp>
#include <vector>
#include <numeric>
#include "common.h"
#include "hits.h"
#include "scaffolds.h"
#include "gtf_tracking.h"
#include "progressbar.h"

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
private:
    HitBundle(const HitBundle& rhs) {} 
public:
	HitBundle() 
    : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id), _ref_id(0), _raw_mass(0.0), _num_replicates(1), _compatible_mass(0.0) {}
	
    ~HitBundle()
    {
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        foreach(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            // This bundle and the factory that actually owns the ref_mRNAs
            // are the only objects that should have access to these scaffolds
            // so if the use count is 2, we can clear these guys.
			// Updated to 3 since the bias learner now uses them.
            if (ref_scaff.use_count() <= 3)
            {
                ref_scaff->clear_hits();
            }
            else if (ref_scaff->mate_hits().size() > 0)
            {
                fprintf(stderr, "Warning: bundle %d-%d shared reference scaffolds with others.  Possible soft memory leak.\n", left(), right());
            }
        }   
        
        foreach (MateHit& hit, _hits)
		{
			delete hit.left_alignment();
			delete hit.right_alignment();
		}
        
        for(OpenMates::iterator itr = _open_mates.begin(); itr != _open_mates.end(); ++itr)
		{
			foreach (MateHit& hit,  itr->second)
            {
                delete hit.left_alignment();
                delete hit.right_alignment();
            }
		}
		
    }
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
    
	// This is to keep track of mass of all hits, including
	// thosethat are not added to any bundle
	// but are skipped during the creation of this bundle
	void add_raw_mass(double rm) { _raw_mass += rm; }
	void rem_raw_mass(double rm) { _raw_mass -= rm; }
	double raw_mass() { return _raw_mass; }
    
    double compatible_mass() const 
    {
        return _compatible_mass;
    }
    
	void clear_hits() 
    {
        _hits.clear(); 
        _non_redundant.clear();
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        foreach(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            if (ref_scaff.use_count() <= 3)
            {
                ref_scaff->clear_hits();
            }
            else 
            {
                fprintf(stderr, "Warning: bundle %d-%d shared reference scaffolds with others.  Possible soft memory leak.\n", left(), right());
            }
        } 
    }
	
    const std::vector<MateHit>& hits() const { return _hits; } 
	const std::vector<MateHit>& non_redundant_hits() const { return _non_redundant; } 
	
	RefID ref_id()  const {return _ref_id; }
	
	int id() const { return _id; }
	
	void add_ref_scaffold(shared_ptr<Scaffold> scaff)
	{
		if (scaff->left() < _leftmost)
			_leftmost = scaff->left();
		if (scaff->right() > _rightmost)
			_rightmost = scaff->right();
		_ref_scaffs.push_back(scaff);
        _ref_id = scaff->ref_id();
	}
	
	vector<shared_ptr<Scaffold> >& ref_scaffolds() { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	bool add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                      const ReadHit* bh,
					  bool expand_by = true);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	// Sorts the hits, and performs other functions needed to prepare this 
	// bundle for assembly and quantitation
	void finalize(bool is_combined=false);
	
	void remove_hitless_scaffolds();

	void collapse_hits();
	
    int num_replicates() const { return _num_replicates; }
    
	double mass() const
	{
		double mass = 0;
		for(size_t i = 0; i < _non_redundant.size(); i++)
		{
			mass += _non_redundant[i].collapse_mass();
		}
		return mass;
	}
	
    static void combine(const vector<HitBundle*>& in_bundles,
                        HitBundle& out_bundle);
    
private:
    int _leftmost;
	int _rightmost;
	std::vector<MateHit> _hits;
	std::vector<MateHit> _non_redundant;
	std::vector<shared_ptr<Scaffold> > _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
    RefID _ref_id;
	double _raw_mass;

	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	OpenMates _open_mates;
    int _num_replicates;
    double _compatible_mass;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs=false,
				   bool loadFPKM=false);

class BundleFactory
{
public:
    
	BundleFactory(shared_ptr<HitFactory> fac, BundleMode bm)
	: _hit_fac(fac), _bundle_mode(bm), _prev_pos(0), _prev_ref_id(0), _curr_bundle(0),  _zeroone(rng)
	{
		_rg_props = shared_ptr<ReadGroupProperties>(new ReadGroupProperties(fac->read_group_properties()));
        
        
       
	}

    bool bundles_remain()  
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        return _curr_bundle < num_bundles();
    }
    
	bool next_bundle(HitBundle& bundle_out);
	bool next_bundle_hit_driven(HitBundle& bundle_out);
	bool next_bundle_ref_driven(HitBundle& bundle_out);
	bool next_bundle_ref_guided(HitBundle& bundle_out);

    
    RefSequenceTable& ref_table() { return _hit_fac->ref_table(); }
    
    // Not available until after inspect_bundle
    int num_bundles() const { return _num_bundles; }
    void num_bundles(int n) { _num_bundles = n; }
    
	void reset() 
	{ 
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _curr_bundle = 0;
		//rewind(hit_file); 
		_hit_fac->reset();
		next_ref_scaff = ref_mRNAs.begin(); 
        next_mask_scaff = mask_gtf_recs.begin();
        
        foreach(shared_ptr<Scaffold> ref_scaff, ref_mRNAs)
        {
            ref_scaff->clear_hits();
        }
        
        _prev_pos = 0;
        _prev_ref_id = 0;
	}
	
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other
    void set_ref_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        ref_mRNAs.clear();
        for (vector<shared_ptr<Scaffold> >::const_iterator i = mRNAs.begin(); i < mRNAs.end(); ++i)
        {
            ref_mRNAs.push_back(shared_ptr<Scaffold>(new Scaffold(**i)));
        }
        
        RefID last_id = 0;
        for (vector<shared_ptr<Scaffold> >::iterator i = ref_mRNAs.begin(); i < ref_mRNAs.end(); ++i)
        {
            if ((*i)->ref_id() != last_id)
            {
                _ref_scaff_offsets.push_back(make_pair((*i)->ref_id(), i));
            }
            last_id = (*i)->ref_id();
        }
        
        next_ref_scaff = ref_mRNAs.begin();
    }
    
    void set_mask_rnas(const vector<shared_ptr<Scaffold> >& masks)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        mask_gtf_recs = masks;
        RefID last_id = 0;
        for (vector<shared_ptr<Scaffold> >::iterator i = mask_gtf_recs.begin(); i < mask_gtf_recs.end(); ++i)
        {
            if ((*i)->ref_id() != last_id)
            {
                _mask_scaff_offsets.push_back(make_pair((*i)->ref_id(), i));
            }
            last_id = (*i)->ref_id();
        }
        
        next_mask_scaff = mask_gtf_recs.begin();
    }
        
	void bad_intron_table(const BadIntronTable& bad_introns) 
	{ 
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
		_bad_introns = bad_introns;
	}
    
    void read_group_properties(shared_ptr<ReadGroupProperties> rg)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _rg_props = rg;
    }
    
    shared_ptr<ReadGroupProperties> read_group_properties()
    {
        return _rg_props;
    }
	
	bool spans_bad_intron(const ReadHit& read);
	
private:
	
	bool _expand_by_hits(HitBundle& bundle);
	bool _expand_by_refs(HitBundle& bundle);
	
	shared_ptr<HitFactory> _hit_fac;
    
	vector<shared_ptr<Scaffold> > ref_mRNAs;
	//FILE* ref_mRNA_file;
	vector<pair<RefID, vector<shared_ptr<Scaffold> >::iterator> > _ref_scaff_offsets;
	vector<shared_ptr<Scaffold> >::iterator next_ref_scaff;
    
    vector<shared_ptr<Scaffold> > mask_gtf_recs;
	//FILE* mask_file;
	vector<pair<RefID, vector<shared_ptr<Scaffold> >::iterator> > _mask_scaff_offsets;
	vector<shared_ptr<Scaffold> >::iterator next_mask_scaff;
	
	BadIntronTable _bad_introns;
    
    shared_ptr<ReadGroupProperties> _rg_props;
    
	// Sets nva to point to the next valid alignment
	// Returns the mass of any alignments that are seen, valid or not
    double next_valid_alignment(const ReadHit*& nva);
	
	// Backs up the factory to before the last valid alignment
	// and returns the mass of that alignment (rh)
	double rewind_hit(const ReadHit* rh);

    BundleMode _bundle_mode;
    int _prev_pos;
    RefID _prev_ref_id;
    int _num_bundles;
    int _curr_bundle;
    
    boost::mt19937 rng;
    boost::uniform_01<boost::mt19937> _zeroone;
    
#if ENABLE_THREADS    
    boost::mutex _factory_lock;
#endif
};

void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

void inspect_map(BundleFactory& bundle_factory,
					BadIntronTable* bad_introns,
					vector<LocusCount>& count_table,
					bool progress_bar = true,
					bool show_stats = true);

#endif
