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
	void remove_unmapped_hits();
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
    boost::mt19937 rng;
    boost::uniform_01<boost::mt19937> _zeroone;
    
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
#if ENABLE_THREADS    
    boost::mutex _factory_lock;
#endif
};

void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

template<class BundleFactoryType>
void inspect_map(BundleFactoryType& bundle_factory,
                 BadIntronTable* bad_introns,
                 vector<LocusCount>& count_table,
                 bool progress_bar = true)
{

	ProgressBar p_bar;
	if (progress_bar)
		p_bar = ProgressBar("Inspecting reads and determining fragment length distribution.",bundle_factory.ref_table().size());
	RefID last_chrom = 0;

	long double map_mass = 0.0;
    long double norm_map_mass = 0.0;
	
	int min_len = numeric_limits<int>::max();
	int max_len = def_max_frag_len;
	vector<double> frag_len_hist(def_max_frag_len+1,0);
	bool has_pairs = false;

	int num_bundles = 0;
	size_t total_hits = 0;
	size_t total_non_redundant_hits = 0;
	
	//To be used for quartile normalization
	vector<long double> mass_dist; 	
	
	// Store the maximum read length for "first" and "second" reads to report to user.
	int max_1 = 0;
	int max_2 = 0;
	
	shared_ptr<MultiReadTable> mrt(new MultiReadTable());
	
	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		bool valid_bundle = bundle_factory.next_bundle(*bundle_ptr);
		HitBundle& bundle = *bundle_ptr;

        if (use_compat_mass) //only count hits that are compatible with ref transcripts
        {
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.compatible_mass();
            if (use_quartile_norm && bundle.compatible_mass() > 0) 
            {
                mass_dist.push_back(bundle.compatible_mass());
            }
        }
        else if (use_total_mass) //use all raw mass
        { 
            
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.raw_mass();
            if (use_quartile_norm && bundle.raw_mass() > 0) 
            {
                mass_dist.push_back(bundle.raw_mass());
            }
        }
        else
        {
            fprintf(stderr, "Error: hit counting scheme for normalization is not set!\n");
            assert(false);
            exit(1);
        }
		
		const RefSequenceTable& rt = bundle_factory.ref_table();
		const char* chrom = rt.get_name(bundle.ref_id());
		char bundle_label_buf[2048];
        if (chrom)
        {
            sprintf(bundle_label_buf, "%s:%d-%d", chrom, bundle.left(), bundle.right());
            verbose_msg("Inspecting bundle %s with %lu reads\n", bundle_label_buf, bundle.hits().size());
            count_table.push_back(LocusCount(bundle_label_buf, bundle.raw_mass(), bundle.ref_scaffolds().size()));
		}
        
        if (!valid_bundle)
		{
			delete bundle_ptr;
			break;
		}
		num_bundles++;
        
        if (progress_bar) 
        {
			double inc_amt = last_chrom == bundle.ref_id() ? 0.0 : 1.0;
			p_bar.update(bundle_label_buf, inc_amt);
			last_chrom = bundle.ref_id();
        }
        
        if (bad_introns != NULL)
		{
			identify_bad_splices(bundle, *bad_introns);
		}
		
		const vector<MateHit>& hits = bundle.non_redundant_hits();
		if (hits.empty())
		{
			delete bundle_ptr;
			continue;
		}
		
		list<pair<int, int> > open_ranges;
		int curr_range_start = hits[0].left();
		int curr_range_end = numeric_limits<int>::max();
		int next_range_start = -1;
		
		total_non_redundant_hits += bundle.non_redundant_hits().size();
		total_hits += bundle.hits().size();
		
		// This first loop calclates the map mass and finds ranges with no introns
		// Note that we are actually looking at non-redundant hits, which is why we use collapse_mass
		// This loop will also add multi-reads to the MultiReads table 
		for (size_t i = 0; i < hits.size(); ++i) 
		{
			assert(hits[i].left_alignment());
            
            // Add to table if multi-read
			if (hits[i].is_multi())
			{
				mrt->add_hit(hits[i]);
			}
			
			// Find left length
			int left_len = hits[i].left_alignment()->right()-hits[i].left_alignment()->left();
			min_len = min(min_len, left_len);
			if (!hits[i].left_alignment()->contains_splice())
            {
				if (hits[i].left_alignment()->is_first())
                    max_1 = max(max_1, left_len);
                else
                    max_2 = max(max_2, left_len);
            }
			
			// Find right length
			if (hits[i].right_alignment())
			{
				int right_len = hits[i].right_alignment()->right()-hits[i].right_alignment()->left();
				min_len = min(min_len, right_len);
				if (!hits[i].right_alignment()->contains_splice())
                {
                    if (hits[i].right_alignment()->is_first())
                        max_1 = max(max_1, right_len);
                    else
                        max_2 = max(max_2, right_len);
                }
                has_pairs = true;
			}
			
			// Find fragment length
			if (bundle.ref_scaffolds().size()==1 && hits[i].is_pair())
			// Annotation provided and single isoform gene
			{
				int start, end, mate_length;
				shared_ptr<Scaffold> scaff = bundle.ref_scaffolds()[0];
				if (scaff->map_frag(hits[i], start, end, mate_length))
				{
					if (mate_length >= min_len && mate_length <= max_len)
						frag_len_hist[mate_length] += hits[i].collapse_mass();
				}
			}
			else if (bundle.ref_scaffolds().empty())
			// No annotation provided.  Look for ranges.
			{
				if (hits[i].left() > curr_range_end)
				{
					if (curr_range_end - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, curr_range_end));
					curr_range_start = next_range_start;
					curr_range_end = numeric_limits<int>::max();
				}
				if (hits[i].left_alignment()->contains_splice())
				{
					if (hits[i].left() - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, hits[i].left()-1));
					curr_range_start = max(next_range_start, hits[i].left_alignment()->right());
				}
				if (hits[i].right_alignment() && hits[i].right_alignment()->contains_splice())
				{
					assert(hits[i].right_alignment()->left() >= hits[i].left());
					curr_range_end = min(curr_range_end, hits[i].right_alignment()->left()-1);
					next_range_start = max(next_range_start, hits[i].right());
				}
			}
		}
        
        if (bundle.ref_scaffolds().empty() && has_pairs) // No annotation provided
		{
			pair<int, int> curr_range(-1,-1);
			
			// This second loop uses the ranges found above to find the estimated frag length distribution
			// It also finds the minimum read length to use in the linear interpolation
			for (size_t i = 0; i < hits.size(); ++i)
			{
				if (hits[i].left() > curr_range.second && open_ranges.empty())
					break;
				
				if (hits[i].left() > curr_range.second)
				{
					curr_range = open_ranges.front();
					open_ranges.pop_front();
				}
				
				if (hits[i].left() >= curr_range.first && hits[i].right() <= curr_range.second && hits[i].is_pair())
				{
					int mate_len = hits[i].right()-hits[i].left();
					if (mate_len <= max_len)
						frag_len_hist[mate_len] += hits[i].collapse_mass();
				}
			}
		}
		
        open_ranges.clear();
		delete bundle_ptr;
	}
	
    norm_map_mass = map_mass;
    
	if (use_quartile_norm && mass_dist.size() > 0)
	{
		sort(mass_dist.begin(),mass_dist.end());
		int upper_quart_index = mass_dist.size() * 0.75;
		norm_map_mass = mass_dist[upper_quart_index];
	}

    if (bad_introns != NULL)
    {
        size_t alloced = 0;
        size_t used = 0;
        size_t num_introns = 0;
        for (BadIntronTable::const_iterator itr = bad_introns->begin();
             itr != bad_introns->end();
             ++itr)
        {
            alloced += itr->second.capacity() * sizeof(AugmentedCuffOp);
            used += itr->second.size() * sizeof(AugmentedCuffOp);
            num_introns += itr->second.size();
        }
        
        verbose_msg( "Bad intron table has %lu introns: (%lu alloc'd, %lu used)\n", num_introns, alloced, used);
    	verbose_msg( "Map has %lu hits, %lu are non-redundant\n", total_hits, total_non_redundant_hits);
    } 
    
	if (progress_bar)
		p_bar.complete();
	
	vector<double> frag_len_pdf(max_len+1, 0.0);
	vector<double> frag_len_cdf(max_len+1, 0.0);
    long double tot_count = accumulate(frag_len_hist.begin(), frag_len_hist.end(), 0.0 );
    bool empirical = false;
	
	if (user_provided_fld && has_pairs && tot_count >= 10000)
	{
		fprintf(stderr, "Warning: Overriding empirical fragment length distribution with user-specified parameters is not recommended.\n");
	}
	
	if (!has_pairs || tot_count < 10000)
	{
		if (has_pairs && !user_provided_fld)
		{
			fprintf(stderr, "Warning: Using default Gaussian distribution due to insufficient paired-end reads in open ranges.  It is recommended that correct parameters (--frag-len-mean and --frag-len-std-dev) be provided.\n");
		}
		tot_count = 0;
		normal frag_len_norm(def_frag_len_mean, def_frag_len_std_dev);
		max_len = def_frag_len_mean + 3*def_frag_len_std_dev;
		for(int i = min_len; i <= max_len; i++)
		{
			frag_len_hist[i] = cdf(frag_len_norm, i+0.5)-cdf(frag_len_norm, i-0.5);
			tot_count += frag_len_hist[i];
		}
	}
	else
	// Calculate the max frag length and interpolate all zeros between min read len and max frag len
	{	
		empirical = true;
		double curr_total = 0;
		size_t last_nonzero = min_len-1;
		for(size_t i = last_nonzero+1; i < frag_len_hist.size(); i++)
		{
			if (frag_len_hist[i] > 0)
			{
				if (last_nonzero != i-1)
				{
					double b = frag_len_hist[last_nonzero];
					double m = (frag_len_hist[i] - b)/(i-last_nonzero);
					for (size_t x = 1; x < i - last_nonzero; x++)
					{
						frag_len_hist[last_nonzero+x] = m * x + b;
						tot_count += frag_len_hist[last_nonzero+x];
						curr_total += frag_len_hist[last_nonzero+x];
					}	
				}
				last_nonzero = i;
			}
			
			curr_total += frag_len_hist[i];
			
			if (curr_total/tot_count > 0.9999)
			{
				max_len = i; 
				tot_count = curr_total;
				break;
			}
		}
	}
	
    double mean = 0.0;

    if (output_fld)
    {
        FILE* fhist = fopen(string(output_dir + "/frag_len_hist.csv").c_str(),"w");
        fprintf(fhist, "Length,Count\n");
        for(size_t i = 1; i < frag_len_hist.size(); i++)
        {
            fprintf(fhist, "%zu,%f\n", i, frag_len_hist[i]);
        }
        fclose(fhist);
    }

	// Convert histogram to pdf and cdf, calculate mean
	int frag_len_mode = 0;
	for(size_t i = min_len; i <= (size_t)max_len; i++)
	{
		frag_len_pdf[i] = frag_len_hist[i]/tot_count;
		frag_len_cdf[i] = frag_len_cdf[i-1] + frag_len_pdf[i];
        
		if (frag_len_pdf[i] > frag_len_pdf[frag_len_mode])
			frag_len_mode = i;
        mean += frag_len_pdf[i] * i;
	}
    
    double std_dev =  0.0;
    for(size_t i = 1; i < frag_len_hist.size(); i++)
    {
        std_dev += frag_len_pdf[i] * ((i - mean) * (i - mean));
    }
    
    std_dev = sqrt(std_dev);
	
	shared_ptr<ReadGroupProperties> rg_props = bundle_factory.read_group_properties();
	shared_ptr<EmpDist const> fld(new EmpDist(frag_len_pdf, frag_len_cdf, frag_len_mode, mean, std_dev, min_len, max_len));
	rg_props->multi_read_table(mrt);
	rg_props->frag_len_dist(fld);
	rg_props->normalized_map_mass(norm_map_mass);
    rg_props->total_map_mass(map_mass);

	fprintf(stderr, "> Map Properties:\n");
	if (use_quartile_norm)
		fprintf(stderr, ">\tUpper Quartile: %.2Lf\n", norm_map_mass);
	else
		fprintf(stderr, ">\tTotal Map Mass: %.2Lf\n", norm_map_mass);
	if (corr_multi)
		fprintf(stderr,">\tNumber of Multi-Reads: %zu (with %zu total hits)\n", mrt->num_multireads(), mrt->num_multihits()); 
//	if (has_pairs)
//		fprintf(stderr, ">\tRead Type: %dbp x %dbp\n", max_1, max_2);
//	else
//		fprintf(stderr, ">\tRead Type: %dbp single-end\n", max(max_1,max_2));

	if (empirical)
	{
		fprintf(stderr, ">\tFragment Length Distribution: Empirical (learned)\n");
		fprintf(stderr, ">\t              Estimated Mean: %.2f\n", mean);
		fprintf(stderr, ">\t           Estimated Std Dev: %.2f\n", std_dev);
	}
	else
	{
		if (user_provided_fld)
			fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (user-specified)\n");
		else
			fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (default)\n");
		fprintf(stderr, ">\t              Default Mean: %d\n", def_frag_len_mean);
		fprintf(stderr, ">\t           Default Std Dev: %d\n", def_frag_len_std_dev);
	}

	bundle_factory.num_bundles(num_bundles);
	bundle_factory.reset(); 
	return;
}

#endif
