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
#include <numeric>
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
	HitBundle() 
    : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id), _num_replicates(1) {}
    
	~HitBundle() {}
	
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
	
	const std::vector<MateHit>& hits() const { return _hits; } 
	const std::vector<MateHit>& non_redundant_hits() const { return _non_redundant; } 
	
	RefID ref_id()  const
	{
		if (!_hits.empty())
			return _hits.front().ref_id();
		else if (!_ref_scaffs.empty())
			return _ref_scaffs.front()->ref_id();
		else
			return 0;
	}
	
	int id() const { return _id; }
	
	void add_ref_scaffold(shared_ptr<Scaffold> scaff)
	{
		_ref_scaffs.push_back(scaff);
	}
	
	vector<shared_ptr<Scaffold> >& ref_scaffolds() { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	void add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                      shared_ptr<ReadHit> bh);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	// Sorts the hits, and performs other functions needed to prepare this 
	// bundle for assembly and quantitation
	void finalize(bool is_combined=false);
	
	void remove_hitless_scaffolds();
	
    int num_replicates() const { return _num_replicates; }
    
    static void combine(const vector<HitBundle>& in_bundles,
                        HitBundle& out_bundle)
    {
        out_bundle._hits.clear();
		out_bundle._non_redundant.clear();
        out_bundle._ref_scaffs.clear();
        
        for (size_t i = 1; i < in_bundles.size(); ++i)
        {
            assert(in_bundles[i].ref_id() == in_bundles[i-1].ref_id());
        }

		// Merge  hits
		vector<int> indices(in_bundles.size(),0);
		while(true)
		{
			int next_bundle = -1;
			const MateHit* next_hit; 
			for(int i = 0; i < in_bundles.size(); ++i)
			{
				const vector<MateHit>& curr_hits = in_bundles[i].hits();
				
				if (indices[i] == curr_hits.size())
					continue;
				
				const MateHit* curr_hit = &curr_hits[indices[i]];
				
				if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
				{
					next_bundle = i;
					next_hit = curr_hit;
				}
			}
			
			if(next_bundle==-1)
				break;
			
			out_bundle._hits.push_back(*next_hit);
			indices[next_bundle]++;
		}
		
		// Merge collapsed hits
		indices = vector<int>(in_bundles.size(), 0);
		while(true)
		{
			int next_bundle = -1;
			const MateHit* next_hit; 
			for(int i = 0; i < in_bundles.size(); ++i)
			{
				const vector<MateHit>& curr_non_redundant_hits = in_bundles[i].non_redundant_hits();
				
				if (indices[i] == curr_non_redundant_hits.size())
					continue;
				
				const MateHit* curr_hit = &curr_non_redundant_hits[indices[i]];
				
				if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
				{
					next_bundle = i;
					next_hit = curr_hit;
				}
			}
			
			if(next_bundle==-1)
				break;
			
			out_bundle._non_redundant.push_back(*next_hit);
			indices[next_bundle]++;
		}
        
		// Merge ref scaffolds
		indices = vector<int>(in_bundles.size(), 0);
		while(true)
		{
			int next_bundle = -1;
			shared_ptr<Scaffold> next_scaff; 
			for(int i = 0; i < in_bundles.size(); ++i)
			{
				const vector<shared_ptr<Scaffold> >& curr_scaffs = in_bundles[i]._ref_scaffs;
				
				if (indices[i] == curr_scaffs.size())
					continue;
				
				shared_ptr<Scaffold> curr_scaff = curr_scaffs[indices[i]];
				
				if (next_bundle == -1 || scaff_lt_rt_oplt(*curr_scaff, *next_scaff))
				{
					next_bundle = i;
					next_scaff = curr_scaff;
				}
			}
			
			if(next_bundle==-1)
				break;
			
			if (out_bundle._ref_scaffs.size()==0 || out_bundle._ref_scaffs.back()->annotated_trans_id() != next_scaff->annotated_trans_id()) 
				out_bundle._ref_scaffs.push_back(next_scaff);
			indices[next_bundle]++;
		}
		        
        out_bundle.finalize(true); // true means everything is already sorted, etc.
        out_bundle._num_replicates = (int)in_bundles.size();
    }
    
private:
    int _leftmost;
	int _rightmost;
	std::vector<MateHit> _hits;
	std::vector<MateHit> _non_redundant;
	std::vector<shared_ptr<Scaffold> > _ref_scaffs; // user-supplied reference annotations overlapping the bundle
	bool _final;
	int _id;
	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	OpenMates _open_mates;
    int _num_replicates;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs=false,
				   bool loadFPKM=false);

class BundleFactory
{
public:
    
	BundleFactory(HitFactory& fac)
	: _hit_fac(fac) {}

	bool next_bundle(HitBundle& bundle_out);
    
    RefSequenceTable& ref_table() { return _hit_fac.ref_table(); }
    
	void reset() 
	{ 
		//rewind(hit_file); 
		_hit_fac.reset();
		next_ref_scaff = ref_mRNAs.begin(); 
        next_mask_scaff = mask_gtf_recs.begin();
	}
	
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other
    void set_ref_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
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
		_bad_introns = bad_introns;
	}
    
    void read_group_properties(shared_ptr<ReadGroupProperties> rg)
    {
        _rg_props = rg;
    }
    
    shared_ptr<ReadGroupProperties> read_group_properties()
    {
        return _rg_props;
    }
	
	bool spans_bad_intron(const ReadHit& read);
	
private:
	HitFactory& _hit_fac;
    
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
    
    shared_ptr<ReadHit> next_valid_alignment();
};

void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

template<class BundleFactoryType>
void inspect_map(BundleFactoryType& bundle_factory,
                 long double& map_mass, 
                 BadIntronTable& bad_introns,
                 EmpDist& frag_len_dist)
{
	fprintf(stderr, "Inspecting reads and determining empirical fragment length distribution.\n");
    
	HitBundle bundle;
    map_mass = 0.0;
    int min_len = numeric_limits<int>::max();
	int max_len = def_max_frag_len;
    vector<double> frag_len_hist(def_max_frag_len+1,0);
	list<pair<int, int> > open_ranges;
    bool has_pairs = false;
    
    size_t total_hits = 0;
    size_t total_non_redundant_hits = 0;
    
	while(bundle_factory.next_bundle(bundle))
	{
//#if ASM_VERBOSE
        const RefSequenceTable& rt = bundle_factory.ref_table();
        const char* chrom = rt.get_name(bundle.ref_id());
        char bundle_label_buf[2048];
        sprintf(bundle_label_buf, 
                "%s:%d-%d",
                chrom,
                bundle.left(),
                bundle.right());
        fprintf(stderr, "Inspecting bundle %s with %lu reads\n", bundle_label_buf, bundle.hits().size());
//#endif
		
		identify_bad_splices(bundle, bad_introns);
        
        const vector<MateHit>& hits = bundle.non_redundant_hits();
		if (hits.empty())
            continue;
        
		int curr_range_start = hits[0].left();
		int curr_range_end = numeric_limits<int>::max();
		int next_range_start = -1;
		
        total_non_redundant_hits += bundle.non_redundant_hits().size();
        total_hits += bundle.hits().size();

		// This first loop calclates the map mass and finds ranges with no introns
		for (size_t i = 0; i < hits.size(); ++i) 
		{
  			map_mass += hits[i].collapse_mass(); 
            
			min_len = min(min_len, hits[i].left_alignment()->right()-hits[i].left_alignment()->left());
			if (hits[i].right_alignment())
				min_len = min(min_len, hits[i].right_alignment()->right()-hits[i].right_alignment()->left());
            
			
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
				has_pairs = true;
				assert(hits[i].right_alignment()->left() > hits[i].left());
				curr_range_end = min(curr_range_end, hits[i].right_alignment()->left()-1);
				next_range_start = max(next_range_start, hits[i].right());
			}
		}
        
        if (has_pairs)
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
					if (mate_len < max_len)
						frag_len_hist[mate_len] += hits[i].collapse_mass();
					min_len = min(min_len, hits[i].left_alignment()->right()-hits[i].left_alignment()->left());
					min_len = min(min_len, hits[i].right_alignment()->right()-hits[i].right_alignment()->left());
				}
			}
		}
        
        
        foreach(shared_ptr<Scaffold>& ref_scaff, bundle.ref_scaffolds())
        {
            ref_scaff->clear_hits();
        }
	}
	
    size_t alloced = 0;
    size_t used = 0;
    size_t num_introns = 0;
    for (BadIntronTable::const_iterator itr = bad_introns.begin();
         itr != bad_introns.end();
         ++itr)
    {
        alloced += itr->second.capacity() * sizeof(AugmentedCuffOp);
        used += itr->second.size() * sizeof(AugmentedCuffOp);
        num_introns += itr->second.size();
    }
    
#if ASM_VERBOSE
    fprintf(stderr, "Bad intron table has %lu introns: (%lu alloc'd, %lu used)\n", num_introns, alloced, used);
#endif
    fprintf(stderr, "Map has %lu hits, %lu are non-redundant\n", total_hits, total_non_redundant_hits);
            
    long double tot_count = 0;
	vector<double> frag_len_pdf(max_len+1, 0.0);
	vector<double> frag_len_cdf(max_len+1, 0.0);
	normal frag_len_norm(def_frag_len_mean, def_frag_len_std_dev);
        
    tot_count = accumulate(frag_len_hist.begin(), frag_len_hist.end(), 0.0 );
    
	// Calculate the max frag length and interpolate all zeros between min read len and max frag len
	if (!has_pairs || tot_count == 0)
	{
		max_len = def_frag_len_mean + 3*def_frag_len_std_dev;
		for(int i = min_len; i < max_len; i++)
		{
			frag_len_hist[i] = cdf(frag_len_norm, i+0.5)-cdf(frag_len_norm, i-0.5);
			tot_count += frag_len_hist[i];
		}
	}
	else 
	{	
		double curr_total = 0;
		int last_nonzero = min_len-1;
		for(int i = 1; i < frag_len_hist.size(); i++)
		{
			if (frag_len_hist[i] > 0)
			{
				if (last_nonzero > 0 && last_nonzero != i-1)
				{
					double b = frag_len_hist[last_nonzero];
					double m = (frag_len_hist[i] - b)/(i-last_nonzero);
					for (int x = 1; x < i - last_nonzero; x++)
					{
						frag_len_hist[last_nonzero+x] = m * x + b;
						tot_count += frag_len_hist[last_nonzero+x];
						curr_total += frag_len_hist[last_nonzero+x];
					}	
				}
				last_nonzero = i;
			}
			
			curr_total += frag_len_hist[i];
			
			if (curr_total/tot_count > .9999)
			{
				max_len = i; 
				break;
			}
		}
	}
	
    double mean = 0.0;
    
	// Convert histogram to pdf and cdf, calculate mean
	int frag_len_mode = 0;
	for(int i = 1; i < frag_len_hist.size(); i++)
	{
		frag_len_pdf[i] = frag_len_hist[i]/tot_count;
		frag_len_cdf[i] = frag_len_cdf[i-1] + frag_len_pdf[i];
		//fprintf(stderr, "%f\n", frag_len_hist[i]);
        
		if (frag_len_pdf[i] > frag_len_pdf[frag_len_mode])
			frag_len_mode = i;
        mean += frag_len_pdf[i] * i;
	}
    
    double std_dev =  0.0;
    for(int i = 1; i < frag_len_hist.size(); i++)
    {
        std_dev += frag_len_pdf[i] * ((i - mean) * (i - mean));
    }
    
    std_dev = sqrt(std_dev);
	
	frag_len_dist.pdf(frag_len_pdf);
	frag_len_dist.cdf(frag_len_cdf);
	frag_len_dist.mode(frag_len_mode);
	frag_len_dist.max(max_len);
	frag_len_dist.min(min_len);
    frag_len_dist.mean(mean);
    frag_len_dist.std_dev(std_dev);
    
    bundle_factory.reset();
    
   	return;
}


#endif
