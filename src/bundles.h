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
#include <boost/unordered_map.hpp>
#include <vector>
#include <numeric>
#include <cmath>
#include "common.h"
#include "hits.h"
#include "scaffolds.h"
#include "gtf_tracking.h"
#include "progressbar.h"
#include "rounding.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

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
		: _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id), _ref_id(0), _raw_mass(0.0), _num_replicates(1), _compatible_mass(0.0) {
		// Shrink the open mates table when it is less than 10% full and larger than 10,000 items.
		_rehash_threshold = _open_mates.max_load_factor() / 10;
		_rehash_size_threshold = 10000;
	}
	
    ~HitBundle()
    {
        vector<boost::shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        BOOST_FOREACH(boost::shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
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
        
        BOOST_FOREACH (MateHit& hit, _hits)
		{
			delete hit.left_alignment();
			delete hit.right_alignment();
		}
        
        for(OpenMates::iterator itr = _open_mates.begin(); itr != _open_mates.end(); ++itr)
		{
			delete itr->second.left_alignment();
			delete itr->second.right_alignment();
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
    
    void compatible_mass(double c) { _compatible_mass = c; }
    
	void clear_hits() 
    {
        _hits.clear(); 
        _non_redundant.clear();
        vector<boost::shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        BOOST_FOREACH(boost::shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
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
    void id(int i) { _id = i; }
	
	void add_ref_scaffold(boost::shared_ptr<Scaffold> scaff)
	{
		if (scaff->left() < _leftmost)
			_leftmost = scaff->left();
		if (scaff->right() > _rightmost)
			_rightmost = scaff->right();
		_ref_scaffs.push_back(scaff);
        _ref_id = scaff->ref_id();
	}
	
	vector<boost::shared_ptr<Scaffold> >& ref_scaffolds() { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	bool add_open_hit(boost::shared_ptr<ReadGroupProperties const> rg_props,
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
	std::vector<boost::shared_ptr<Scaffold> > _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
    RefID _ref_id;
	double _raw_mass;

	
	static int _next_id;
	
	// InsertIDs are already hashes of SAM record names; no need
	// to hash them again. Note on 32-bit this will take the truncate
	// the larger InsertID hash to size_t.
	struct identity_hash {
		size_t operator()(const InsertID& in) const { return (size_t)in; }
	};

	typedef boost::unordered_multimap<uint64_t, MateHit, identity_hash> OpenMates;
	OpenMates _open_mates;
    int _num_replicates;
    double _compatible_mass;
	float _rehash_threshold;
	OpenMates::size_type _rehash_size_threshold;

};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<boost::shared_ptr<Scaffold> >& ref_mRNAs,
                   boost::crc_32_type& gtf_crc_result,
				   bool loadSeqs=false,
				   bool loadFPKM=false);

class BundleFactory
{
public:
    
	BundleFactory(boost::shared_ptr<HitFactory> fac, BundleMode bm)
	: _hit_fac(fac), _bundle_mode(bm), _prev_pos(0), _prev_ref_id(0), _curr_bundle(0),  _zeroone(rng)
	{
		_rg_props = boost::shared_ptr<ReadGroupProperties>(new ReadGroupProperties(fac->read_group_properties()));
        next_ref_scaff = ref_mRNAs.begin(); 
        next_mask_scaff = mask_gtf_recs.begin();
	}

    virtual ~BundleFactory() {} 
    boost::shared_ptr<const HitFactory> hit_factory() const { return _hit_fac; }
    
    bool bundles_remain()  
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        return _curr_bundle < num_bundles();
    }
    
	virtual bool next_bundle(HitBundle& bundle_out, bool cache_bundle);
	bool next_bundle_hit_driven(HitBundle& bundle_out);
	bool next_bundle_ref_driven(HitBundle& bundle_out);
	bool next_bundle_ref_guided(HitBundle& bundle_out);

    
    RefSequenceTable& ref_table() { return _hit_fac->ref_table(); }
    
    // Not available until after inspect_bundle
    int num_bundles() const { return _num_bundles; }
    void num_bundles(int n) { _num_bundles = n; }
    
	virtual void reset()
	{ 
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _curr_bundle = 0;
		//rewind(hit_file); 
		_hit_fac->reset();
		next_ref_scaff = ref_mRNAs.begin(); 
        next_mask_scaff = mask_gtf_recs.begin();
        
        BOOST_FOREACH(boost::shared_ptr<Scaffold> ref_scaff, ref_mRNAs)
        {
            ref_scaff->clear_hits();
        }
        
        _prev_pos = 0;
        _prev_ref_id = 0;
	}
	
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other. Deep copying is necessary whenever
    // you need to set fpkm or num fragments in Scaffold objects (e.g. during bias correction or multiread correction)
    void set_ref_rnas(const vector<boost::shared_ptr<Scaffold> >& mRNAs, bool deep_copy = true)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        ref_mRNAs.clear();
        if (deep_copy)
        {
            for (vector<boost::shared_ptr<Scaffold> >::const_iterator i = mRNAs.begin(); i < mRNAs.end(); ++i)
            {
                ref_mRNAs.push_back(boost::shared_ptr<Scaffold>(new Scaffold(**i)));
            }
        }
        else
        {
            ref_mRNAs = mRNAs;
        }
        
        RefID last_id = 0;
        for (vector<boost::shared_ptr<Scaffold> >::iterator i = ref_mRNAs.begin(); i < ref_mRNAs.end(); ++i)
        {
            if ((*i)->ref_id() != last_id)
            {
                _ref_scaff_offsets.push_back(make_pair((*i)->ref_id(), i));
            }
            last_id = (*i)->ref_id();
        }
        
        next_ref_scaff = ref_mRNAs.begin();
    }
    
    void set_mask_rnas(const vector<boost::shared_ptr<Scaffold> >& masks)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        mask_gtf_recs = masks;
        RefID last_id = 0;
        for (vector<boost::shared_ptr<Scaffold> >::iterator i = mask_gtf_recs.begin(); i < mask_gtf_recs.end(); ++i)
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
    
    void read_group_properties(boost::shared_ptr<ReadGroupProperties> rg)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_factory_lock);
#endif
        _rg_props = rg;
    }
    
    boost::shared_ptr<ReadGroupProperties> read_group_properties()
    {
        return _rg_props;
    }
	
	bool spans_bad_intron(const ReadHit& read);
	
    virtual double mode_transcript_coverage(int num_bins = 1000) { return -1; }
    
private:
	
	bool _expand_by_hits(HitBundle& bundle);
	bool _expand_by_refs(HitBundle& bundle);
	
	boost::shared_ptr<HitFactory> _hit_fac;
    
	vector<boost::shared_ptr<Scaffold> > ref_mRNAs;
	//FILE* ref_mRNA_file;
	vector<pair<RefID, vector<boost::shared_ptr<Scaffold> >::iterator> > _ref_scaff_offsets;
	vector<boost::shared_ptr<Scaffold> >::iterator next_ref_scaff;
    
    vector<boost::shared_ptr<Scaffold> > mask_gtf_recs;
	//FILE* mask_file;
	vector<pair<RefID, vector<boost::shared_ptr<Scaffold> >::iterator> > _mask_scaff_offsets;
	vector<boost::shared_ptr<Scaffold> >::iterator next_mask_scaff;
	
	BadIntronTable _bad_introns;
    
    boost::shared_ptr<ReadGroupProperties> _rg_props;
    
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

class PrecomputedExpressionBundleFactory : public BundleFactory
{
public:
    PrecomputedExpressionBundleFactory(boost::shared_ptr<PrecomputedExpressionHitFactory> fac)
	: BundleFactory(fac, REF_DRIVEN), _hit_fac(fac)
	{
		
	}
    
    bool next_bundle(HitBundle& bundle_out, bool cache_bundle);
    
    boost::shared_ptr<const AbundanceGroup> get_abundance_for_locus(int locus_id);
    void clear_abundance_for_locus(int locus_id);
    
    void reset() { BundleFactory::reset(); transcript_coverages.clear(); }
    
//    double median_transcript_coverage() {
//        
//        if (transcript_coverages.size() > 0)
//        {
//            size_t median_idx = transcript_coverages.size() * 0.5;
//            //double mean_transcript_coverage = accumulate(transcript_coverages.begin(), transcript_coverages.end(), 0.0) /transcript_coverages.size();
//            return transcript_coverages[median_idx];
//            //return mean_transcript_coverage;
//        }
//        else
//        {
//            return 0.0;
//        }
//    }

//    double mode_transcript_coverage(int num_bins = 1000) {
//        
//        if (transcript_coverages.size() > 0)
//        {
//            
//            using namespace boost;
//            using namespace boost::accumulators;
//            
//            typedef accumulator_set<double, features<tag::density> > acc;
//            typedef iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type;
//            
//            double min_cov = 99999999;
//            double max_cov = -1;
//            
//            for (size_t i = 0; i < transcript_coverages.size(); ++i)
//            {
//                if (min_cov > transcript_coverages[i])
//                    min_cov = transcript_coverages[i];
//                if (max_cov < transcript_coverages[i])
//                    max_cov = transcript_coverages[i];
//            }
//            
//            if (min_cov == max_cov)
//            {
//                return 0.0;
//            }
//            
//            double cache_size = transcript_coverages.size();
//            
//            // we want to have bin resolution at around 1e-5 between 0-1 frags per KB.
//            num_bins = max_cov / 1e-3;
//            
////            if (transcript_coverages.size() > num_bins)
////            {
////                num_bins = sqrt(transcript_coverages.size());
////            }
////            
//            //create an accumulator
//            acc myAccumulator( tag::density::num_bins = num_bins, tag::density::cache_size = cache_size);
//            
//            //fill accumulator
//            for (size_t i = 0; i < transcript_coverages.size(); ++i)
//            {
//                if (transcript_coverages[i] > 0)
//                {
//                    float nearest = floorf(transcript_coverages[i] * 1e4 + 0.5) / 1e3;
//                    myAccumulator(transcript_coverages[i]);
//                }
//            }
//            
//            histogram_type hist = density(myAccumulator);
//            
//            //double total = 0.0;
//            
//            int mode_idx = -1;
//            double mode_count = -1;
//            for(int i = 0; i < hist.size(); i++ )
//            {
//                //std::cout << "Bin lower bound: " << hist[i].first << ", Value: " << hist[i].second << std::endl;
//                //total += hist[i].second;
//                
//                double lower_b = hist[i].first;
//                double count = hist[i].second;
//                
//                if (hist[i].second > mode_count)
//                {
//                    mode_idx = i;
//                    mode_count = hist[i].second;
//                }
//            }
//            
//            if (mode_idx == -1 || mode_idx >= hist.size() - 1)
//                return 0;
//            
//            double mode_val = (hist[mode_idx+1].first - hist[mode_idx].first) / 2.0;
//            return mode_val;
//            //std::cout << "Total: " << total << std::endl; //should be 1 (and it is)
//
//        }
//        else
//        {
//            return 0.0;
//        }
//    }

    
    double mode_transcript_coverage(int num_bins = 1000) {
        
        if (transcript_coverages.size() > 0)
        {
            size_t max_score_i = 0;
            int max_score = 0;
            for (size_t i = 0; i < transcript_coverages.size(); ++i)
            {
                double score_i = 0;
                for (size_t j = 0; j < transcript_coverages.size(); ++j)
                {
                    if (rounding::roundhalfeven(transcript_coverages[j]/transcript_coverages[i]) == 1)
                    {
                        score_i++;
                    }
                }
                if (score_i > max_score)
                {
                    max_score = score_i;
                    max_score_i = i;
                }
            }
            
            return transcript_coverages[max_score_i];
            
        }
        else
        {
            return 0.0;
        }
    }

    
private:
    
    boost::shared_ptr<PrecomputedExpressionHitFactory> _hit_fac;
#if ENABLE_THREADS
    boost::mutex _factory_lock;
#endif
    
    vector<double> transcript_coverages;
};

void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

class IdToLocusMap
{
    
#if ENABLE_THREADS
	boost::mutex _bl_lock;
#endif
    
    IdToLocusMap(const IdToLocusMap& rhs) {}
public:
    
    IdToLocusMap(boost::shared_ptr<map<string, set<string> > > id_to_locus_map) :
        _id_to_locus_map(id_to_locus_map) {}
    
    void register_locus_to_id(const string& ID, const string& locus_desc)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_bl_lock);
#endif

        pair< map<string, set<string> >::iterator, bool> p = _id_to_locus_map->insert(make_pair(ID, set<string>()));
        p.first->second.insert(locus_desc);
    }
    
    int unregister_locus_from_id(const string& ID, const string& locus_desc)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_bl_lock);
#endif

        map<string, set<string> >::iterator itr = _id_to_locus_map->find(ID);
        if (itr != _id_to_locus_map->end()){
            itr->second.erase(locus_desc);
            return itr->second.size();
        }
        return 0;
    }
    
    int num_locus_registered(const string& ID) {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_bl_lock);
#endif

        map<string, set<string> >::const_iterator itr = _id_to_locus_map->find(ID);
        if (itr == _id_to_locus_map->end())
            return 0;
        return itr->second.size();
    }
private:
    boost::shared_ptr<map<string, set<string> > > _id_to_locus_map;
};

void inspect_map(boost::shared_ptr<BundleFactory> bundle_factory,
                 BadIntronTable* bad_introns,
                 vector<LocusCount>& compatible_count_table,
                 vector<LocusCount>& total_count_table,
                 IdToLocusMap& id_to_locus_map,
                 bool progress_bar = true,
                 bool show_stats = true);

#endif
