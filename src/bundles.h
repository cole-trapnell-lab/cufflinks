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
    : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id), _num_replicates(1) {}
	
    ~HitBundle()
    {
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        foreach(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            // This bundle and the factory that actually owns the ref_mRNAs
            // are the only objects that should have access to these scaffolds
            // so if the use count is 2, we can clear these guys.
            if (ref_scaff.use_count() <= 2)
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
    }
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
    
    
	void clear_hits() 
    {
        _hits.clear(); 
        _non_redundant.clear();
        vector<shared_ptr<Scaffold> >& bundle_ref_scaffs = ref_scaffolds();
        foreach(shared_ptr<Scaffold>& ref_scaff, bundle_ref_scaffs)
        {
            if (ref_scaff.use_count() <= 2)
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
        //scaff->clear_hits();
		_ref_scaffs.push_back(scaff);
	}
	
	vector<shared_ptr<Scaffold> >& ref_scaffolds() { return _ref_scaffs; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	void add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                      const ReadHit* bh);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	// Sorts the hits, and performs other functions needed to prepare this 
	// bundle for assembly and quantitation
	void finalize(bool is_combined=false);
	
	void remove_hitless_scaffolds();
	void remove_unmapped_hits();
	
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
    
	BundleFactory(shared_ptr<HitFactory> fac)
	: _hit_fac(fac), _ref_driven(false), _prev_pos(0), _prev_ref_id(0) {}

	bool next_bundle(HitBundle& bundle_out);
    
    RefSequenceTable& ref_table() { return _hit_fac->ref_table(); }
    
    // Not available until after inspect_bundle
    int num_bundles() const { return _num_bundles; }
    void num_bundles(int n) { _num_bundles = n; }
    
	void reset() 
	{ 
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
        _ref_driven = true;
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
    
    const ReadHit* next_valid_alignment();
    bool _ref_driven;
    int _prev_pos;
    RefID _prev_ref_id;
    int _num_bundles;
};

void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops);

template<class BundleFactoryType>
void inspect_map(BundleFactoryType& bundle_factory,
                 long double& map_mass, 
                 BadIntronTable* bad_introns,
                 EmpDist& frag_len_dist,
                 bool progress_bar = true)
{

	ProgressBar p_bar;
	if (progress_bar)
	{
		p_bar = ProgressBar("Inspecting reads and determining fragment length distribution.",bundle_factory.ref_table().size());
	}

	char last_chrom[100];
    map_mass = 0.0;
    int min_len = numeric_limits<int>::max();
	int max_len = def_max_frag_len;
    vector<double> frag_len_hist(def_max_frag_len+1,0);
    bool has_pairs = false;
    
    int num_bundles = 0;
    size_t total_hits = 0;
    size_t total_non_redundant_hits = 0;
	
	vector<long double> mass_dist; //To be used for quartile normalization
	
	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		if (!bundle_factory.next_bundle(*bundle_ptr))
		{
			delete bundle_ptr;
			break;
		}
		num_bundles++;

		HitBundle& bundle = *bundle_ptr;

		const RefSequenceTable& rt = bundle_factory.ref_table();
		const char* chrom = rt.get_name(bundle.ref_id());	
		char bundle_label_buf[2048];
		sprintf(bundle_label_buf, "%s:%d-%d", chrom, bundle.left(), bundle.right());

		asm_printf("Inspecting bundle %s with %lu reads\n", bundle_label_buf, bundle.hits().size());

		if (progress_bar)
		{
			int inc_amt = (strncmp(last_chrom, chrom, 100)==0) ? 0 : 1;
			p_bar.update(bundle_label_buf, inc_amt);
			strncpy(last_chrom, chrom, 100);
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
		double bundle_mass = 0;
		for (size_t i = 0; i < hits.size(); ++i) 
		{
			bundle_mass += hits[i].collapse_mass();
            
            assert(hits[i].left_alignment());
			min_len = min(min_len, hits[i].left_alignment()->right()-hits[i].left_alignment()->left());
			if (hits[i].right_alignment())
				min_len = min(min_len, hits[i].right_alignment()->right()-hits[i].right_alignment()->left());
            
			if (bundle.ref_scaffolds().size()==1 && hits[i].is_pair()) // Annotation provided and single isoform gene.
			{
				int start, end, mate_length;
				shared_ptr<Scaffold> scaff = bundle.ref_scaffolds()[0];
				if (scaff->map_frag(hits[i], start, end, mate_length))
				{
					if (mate_length >= min_len && mate_length <= max_len)
						frag_len_hist[mate_length] += hits[i].collapse_mass();
				}
			}
			else
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
					has_pairs = true;
					assert(hits[i].right_alignment()->left() > hits[i].left());
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
		
		map_mass += bundle_mass;
		if (use_quartile_norm && bundle_mass > 0) mass_dist.push_back(bundle_mass);
		
		foreach(shared_ptr<Scaffold>& ref_scaff, bundle.ref_scaffolds())
		{
			ref_scaff->clear_hits();
		}
        open_ranges.clear();
        
		delete bundle_ptr;
	}
	
	if (use_quartile_norm)
	{
		sort(mass_dist.begin(),mass_dist.end());
		int num_included = mass_dist.size() * 0.75;
		fprintf(stderr,"%Lf\n",map_mass);
		map_mass = accumulate(mass_dist.begin(), mass_dist.begin()+num_included, 0.0);
		fprintf(stderr,"%Lf\n",map_mass);

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
        
        asm_printf( "Bad intron table has %lu introns: (%lu alloc'd, %lu used)\n", num_introns, alloced, used);
    	asm_printf( "Map has %lu hits, %lu are non-redundant\n", total_hits, total_non_redundant_hits);
    } 
    
    long double tot_count = 0;
	vector<double> frag_len_pdf(max_len+1, 0.0);
	vector<double> frag_len_cdf(max_len+1, 0.0);
        
    tot_count = accumulate(frag_len_hist.begin(), frag_len_hist.end(), 0.0 );
    
    
    string distr_type;
	// Calculate the max frag length and interpolate all zeros between min read len and max frag len
	if (!has_pairs || tot_count == 0)
	{
		distr_type  = "Gaussian (default)";
		normal frag_len_norm(def_frag_len_mean, def_frag_len_std_dev);
		max_len = def_frag_len_mean + 3*def_frag_len_std_dev;
		for(int i = min_len; i < max_len; i++)
		{
			frag_len_hist[i] = cdf(frag_len_norm, i+0.5)-cdf(frag_len_norm, i-0.5);
			tot_count += frag_len_hist[i];
		}
	}
	else 
	{	
		distr_type = "Empirical (learned)";
		double curr_total = 0;
		size_t last_nonzero = min_len-1;
		for(size_t i = last_nonzero+1; i < frag_len_hist.size(); i++)
		{
			if (frag_len_hist[i] > 0)
			{
				if (last_nonzero > 0 && last_nonzero != i-1)
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
			
			if (curr_total/tot_count > .9999)
			{
				max_len = i; 
				break;
			}
		}
	}
	
    double mean = 0.0;
    
#if ADAM_MODE    
    FILE* fhist = fopen(string(output_dir + "/frag_len_hist.csv").c_str(),"w");
    fprintf(fhist, "Length,Count\n");
	for(int i = 1; i < frag_len_hist.size(); i++)
	{
		fprintf(fhist, "%d,%f\n", i, frag_len_hist[i]);
	}
	fclose(fhist);
#endif	

	// Convert histogram to pdf and cdf, calculate mean
	int frag_len_mode = 0;
	for(size_t i = 1; i < frag_len_hist.size(); i++)
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
	
	frag_len_dist.pdf(frag_len_pdf);
	frag_len_dist.cdf(frag_len_cdf);
	frag_len_dist.mode(frag_len_mode);
	frag_len_dist.max(max_len);
	frag_len_dist.min(min_len);
    frag_len_dist.mean(mean);
    frag_len_dist.std_dev(std_dev);
    
   	if (progress_bar)
   	{
   		p_bar.complete();
		fprintf(stderr, "> Map Properties:\n");
		fprintf(stderr, ">\tTotal Map Mass: %.2Lf\n", map_mass);
		string type = (has_pairs) ? "paired-end" : "single-end";
		fprintf(stderr, ">\tRead Type: %dbp %s\n", min_len, type.c_str());
		fprintf(stderr, ">\tFragment Length Distribution: %s\n", distr_type.c_str());
		fprintf(stderr, ">\t                        Mean: %.2f\n", mean);
		fprintf(stderr, ">\t                     Std Dev: %.2f\n", std_dev);
   	}
   	
   	bundle_factory.num_bundles(num_bundles);
    bundle_factory.reset(); 
   	return;
}

#endif
