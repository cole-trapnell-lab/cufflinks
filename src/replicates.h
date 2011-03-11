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

class MassDispersionModel
{
public:
    MassDispersionModel() {}
    MassDispersionModel(const std::vector<double>& scaled_mass_means, 
                        const std::vector<double>& scaled_mass_variances) :
    _scaled_mass_means(scaled_mass_means),
    _scaled_mass_variances(scaled_mass_variances) 
    {
        std::vector<double>::iterator new_end = unique(_scaled_mass_means.begin(), _scaled_mass_means.end());
        _scaled_mass_means.erase(new_end, _scaled_mass_means.end());
        
        new_end = unique(_scaled_mass_variances.begin(), _scaled_mass_variances.end());
        _scaled_mass_variances.erase(new_end, _scaled_mass_variances.end());
    }
    
    virtual double scale_mass_variance(double scaled_mass) const;
    
private:
    std::vector<double> _scaled_mass_means;
    std::vector<double> _scaled_mass_variances;
};

class PoissonDispersionModel : public MassDispersionModel
{
public:
    
    virtual double scale_mass_variance(double scaled_mass) const 
    { 
        return scaled_mass; 
    }
};

void calc_scaling_factors(const std::vector<pair<std::string, std::vector<double> > >& sample_count_table,
                          std::vector<double>& scale_factors);

boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model(const std::vector<double>& scale_factors,
                     const std::vector<std::pair<std::string, std::vector<double> > >& sample_count_table);

// This factory merges bundles in a requested locus from several replicates
class ReplicatedBundleFactory
{
public:
	ReplicatedBundleFactory(const std::vector<shared_ptr<BundleFactory> >& factories)
    : _factories(factories) {}
	
	int num_bundles() { return _factories[0]->num_bundles(); }
    std::vector<boost::shared_ptr<BundleFactory> > factories() { return _factories; }
	
	bool next_bundle(HitBundle& bundle_out)
    {
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
        foreach (shared_ptr<BundleFactory> fac, _factories)
        {
            fac->reset();
        }
    }
    
    void inspect_replicate_maps(int& min_len, int& max_len)
    {
        vector<pair<string, vector<double> > > sample_count_table;
        vector<double> sample_masses;
        
        foreach (shared_ptr<BundleFactory> fac, _factories)
        {
            BadIntronTable bad_introns;
            
            vector<pair<string, double> > count_table;
            inspect_map(*fac, NULL, count_table, false);
            
            shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            
            for (size_t i = 0; i < count_table.size(); ++i)
            {
                pair<string, double>& c = count_table[i];
                double raw_count = c.second;
                
                
                if (i >= sample_count_table.size())
                {
                    sample_count_table.push_back(make_pair(c.first, vector<double>()));
                    sample_count_table.back().second.push_back(raw_count);
                }
                else
                {
                    const string& label = sample_count_table[i].first;
                    assert (label == c.first);
                    sample_count_table[i].second.push_back(raw_count);
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
            rg_props->mass_scale_factor(scale_factors[i]);
        }
        
        // Transform raw counts to the common scale
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            pair<string, vector<double> >& p = sample_count_table[i];
            for (size_t j = 0; j < p.second.size(); ++j)
            {
                assert (scale_factors.size() > j);
                p.second[j] *= (1.0 / scale_factors[j]);
            }
        }
        
        for (size_t i = 0; i < _factories.size(); ++i)
        {
            shared_ptr<ReadGroupProperties> rg_props = _factories[i]->read_group_properties();
            vector<pair<string, double> > scaled_counts;
            for (size_t j = 0; j < sample_count_table.size(); ++j)
            {
                scaled_counts.push_back(make_pair(sample_count_table[j].first, sample_count_table[j].second[i]));
            }
            rg_props->common_scale_counts(scaled_counts);
        }
        
        shared_ptr<MassDispersionModel const> disperser;
        disperser = fit_dispersion_model(scale_factors, sample_count_table);
        
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
        foreach(shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_ref_rnas(mRNAs);
        }
    }
    
    void set_mask_rnas(const vector<shared_ptr<Scaffold> >& mRNAs)
    {
        foreach(shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_mask_rnas(mRNAs);
        }
    }
    
    int num_replicates() const { return _factories.size(); }
    
private:
	vector<shared_ptr<BundleFactory> > _factories;
};
