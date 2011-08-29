//
//  replicates.h
//  cufflinks
//
//  Created by Cole Trapnell on 3/11/11.
//  Copyright 2011 Cole Trapnell. All rights reserved.
//

#include "common.h"
#include "bundles.h"
#include <vector>
#include <algorithm>
#include <map>

class MassDispersionModel
{
public:
    MassDispersionModel() {}
    MassDispersionModel(const std::vector<double>& scaled_mass_means, 
                        const std::vector<double>& scaled_raw_variances,
                        const std::vector<double>& scaled_mass_variances);

    
    virtual double scale_mass_variance(double scaled_mass) const;
    
    const vector<double>& scaled_mass_means() const { return _scaled_mass_means; }
    const vector<double>& scaled_raw_variances() const { return _scaled_raw_variances; }
    const vector<double>& scaled_mass_variances() const { return _scaled_mass_variances; }
    
    std::pair<double, double> get_raw_mean_and_var(const std::string& locus_desc) const
    {
        std::map<std::string, std::pair<double, double> >::const_iterator itr;
        itr = _raw_mv_by_locus.find(locus_desc);
        std::pair<double, double> p = make_pair<double, double>(0.0,0.0);
        if (itr != _raw_mv_by_locus.end())
        {
            p = itr->second;
        }
        return p;
    }
    
    void set_raw_mean_and_var(const std::string& locus_desc, const std::pair<double, double>& p) 
    {
        _raw_mv_by_locus[locus_desc] = p;
    }
    
private:
    std::vector<double> _scaled_mass_means;
    std::vector<double> _scaled_raw_variances;
    std::vector<double> _scaled_mass_variances;
    
    std::map<std::string, std::pair<double, double> > _raw_mv_by_locus;
};

class PoissonDispersionModel : public MassDispersionModel
{
public:
    
    virtual double scale_mass_variance(double scaled_mass) const 
    { 
        return scaled_mass; 
    }
};

struct LocusCountList
{
    LocusCountList(std::string ld, int num_reps, int nt) : 
    locus_desc(ld), counts(std::vector<double>(num_reps, 0)), num_transcripts(nt) {}
    std::string locus_desc;
    std::vector<double> counts;
    int num_transcripts;
};

void calc_scaling_factors(const std::vector<LocusCountList>& sample_count_table,
                          std::vector<double>& scale_factors);


boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model(const string& condition_name,
                     const std::vector<double>& scale_factors,
                     const std::vector<LocusCountList>& sample_count_table);

// This factory merges bundles in a requested locus from several replicates
class ReplicatedBundleFactory
{
public:
	ReplicatedBundleFactory(const std::vector<shared_ptr<BundleFactory> >& factories, 
                            const string& condition_name)
    : _factories(factories), _condition_name(condition_name) {}
	
	int num_bundles() { return _factories[0]->num_bundles(); }
    std::vector<boost::shared_ptr<BundleFactory> > factories() { return _factories; }
	
    const string& condition_name() const { return _condition_name; }
    void condition_name(const string& cn) { _condition_name = cn; }
    
    bool bundles_remain() 
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        foreach (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            if (fac->bundles_remain())
                return true;
        }
        return false;
    }
    
	bool next_bundle(HitBundle& bundle_out)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        std::vector<HitBundle*> bundles;
        
        bool non_empty_bundle = false;
        foreach (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            bundles.push_back(new HitBundle());
            if (fac->next_bundle(*(bundles.back())))
            {
                non_empty_bundle = true;
            }
        }
        
        if (non_empty_bundle == false)
        {
            foreach (HitBundle* in_bundle, bundles)
            {
                in_bundle->ref_scaffolds().clear();
                in_bundle->clear_hits();
                delete in_bundle;
            }
            return false;
        }
        
        for (size_t i = 1; i < bundles.size(); ++i)
        {
            const vector<shared_ptr<Scaffold> >& s1 = bundles[i]->ref_scaffolds();
            const vector<shared_ptr<Scaffold> >& s2 =  bundles[i-1]->ref_scaffolds();
            assert (s1.size() == s2.size());
            for (size_t j = 0; j < s1.size(); ++j)
            {
                assert (s1[j]->annotated_trans_id() == s2[j]->annotated_trans_id());
            }
        }
        
        // Merge the replicates into a combined bundle of hits.
        HitBundle::combine(bundles, bundle_out);
        
        foreach (HitBundle* in_bundle, bundles)
        {
            in_bundle->ref_scaffolds().clear();
            in_bundle->clear_hits();
            delete in_bundle;
        }
        return true;
    }
	
	void reset() 
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        foreach (shared_ptr<BundleFactory> fac, _factories)
        {
            fac->reset();
        }
    }
    
    void inspect_replicate_maps(int& min_len, int& max_len)
    {
        vector<LocusCountList> sample_count_table;
        vector<double> sample_masses;
        
        for (size_t fac_idx = 0; fac_idx < _factories.size(); ++fac_idx)
        {
            shared_ptr<BundleFactory> fac = _factories[fac_idx];        
            BadIntronTable bad_introns;
            
            vector<LocusCount> count_table;
            inspect_map(*fac, NULL, count_table, false);
            
            shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            
            for (size_t i = 0; i < count_table.size(); ++i)
            {
                LocusCount& c = count_table[i];
                double raw_count = c.count;
                
                if (i >= sample_count_table.size())
                {
                    LocusCountList locus_count(c.locus_desc, _factories.size(), c.num_transcripts); 
                    sample_count_table.push_back(locus_count);
                    sample_count_table.back().counts[0] = raw_count;
                }
                else
                {
                    if (sample_count_table[i].locus_desc != c.locus_desc)
                    {
                        fprintf (stderr, "Error: bundle boundaries don't match across replicates!\n");
                        exit(1);
                    }
                    sample_count_table[i].counts[fac_idx] = raw_count;
                }
            }
            sample_masses.push_back(rg_props->total_map_mass());
			min_len = min(min_len, rg_props->frag_len_dist()->min());
			max_len = max(max_len, rg_props->frag_len_dist()->max());
        }
        
        vector<double> scale_factors(_factories.size(), 0.0);
        
        calc_scaling_factors(sample_count_table, scale_factors);
        
        for (size_t i = 0; i < scale_factors.size(); ++i)
        {
            shared_ptr<ReadGroupProperties> rg_props = _factories[i]->read_group_properties();
            assert (scale_factors[i] != 0);
            rg_props->mass_scale_factor(scale_factors[i]);
        }
        
        // Transform raw counts to the common scale
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            LocusCountList& p = sample_count_table[i];
            for (size_t j = 0; j < p.counts.size(); ++j)
            {
                assert (scale_factors.size() > j);
                p.counts[j] *= (1.0 / scale_factors[j]);
            }
        }
        
        for (size_t i = 0; i < _factories.size(); ++i)
        {
            shared_ptr<ReadGroupProperties> rg_props = _factories[i]->read_group_properties();
            vector<LocusCount> scaled_counts;
            for (size_t j = 0; j < sample_count_table.size(); ++j)
            {
                scaled_counts.push_back(LocusCount(sample_count_table[j].locus_desc, sample_count_table[j].counts[i], sample_count_table[j].num_transcripts));
            }
            rg_props->common_scale_counts(scaled_counts);
        }
        
        shared_ptr<MassDispersionModel const> disperser;
        disperser = fit_dispersion_model(_condition_name,scale_factors, sample_count_table);
        
        foreach (shared_ptr<BundleFactory> fac, _factories)
        {
            shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            rg_props->mass_dispersion_model(disperser);
        }
    }
    
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other
    void set_ref_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        foreach(shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_ref_rnas(mRNAs);
        }
    }
    
    void set_mask_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        foreach(shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_mask_rnas(mRNAs);
        }
    }
    
    int num_replicates() const { return _factories.size(); }
    
    void mass_dispersion_model(shared_ptr<MassDispersionModel const> disperser)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        foreach(shared_ptr<BundleFactory>& fac, _factories)
        {
            fac->read_group_properties()->mass_dispersion_model(disperser);
        }
    }
    
    shared_ptr<MassDispersionModel const> mass_dispersion_model() const
    {
        return _factories.front()->read_group_properties()->mass_dispersion_model();
    }
    
private:
	vector<shared_ptr<BundleFactory> > _factories;
#if ENABLE_THREADS
    boost::mutex _rep_factory_lock;
#endif
    string _condition_name; 
};
