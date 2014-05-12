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
    MassDispersionModel(const std::string& name,
                        const std::vector<double>& scaled_compatible_mass_means, 
                        const std::vector<double>& scaled_compatible_variances,
                        const std::vector<double>& scaled_mass_variances);

    virtual const std::string& name() const { return _name; }
    
    virtual double scale_mass_variance(double scaled_mass) const;
    
    const vector<double>& scaled_compatible_mass_means() const { return _scaled_compatible_mass_means; }
    const vector<double>& scaled_compatible_variances() const { return _scaled_compatible_variances; }
    const vector<double>& scaled_mass_variances() const { return _scaled_mass_variances; }
    
    std::pair<double, double> get_compatible_mean_and_var(const std::string& locus_desc) const
    {
        std::map<std::string, std::pair<double, double> >::const_iterator itr;
        itr = _compatible_mv_by_locus.find(locus_desc);
        std::pair<double, double> p = make_pair<double, double>(0.0,0.0);
        if (itr != _compatible_mv_by_locus.end())
        {
            p = itr->second;
        }
        return p;
    }
    
    std::pair<double, double> get_total_mean_and_var(const std::string& locus_desc) const
    {
        std::map<std::string, std::pair<double, double> >::const_iterator itr;
        itr = _total_mv_by_locus.find(locus_desc);
        std::pair<double, double> p = make_pair<double, double>(0.0,0.0);
        if (itr != _total_mv_by_locus.end())
        {
            p = itr->second;
        }
        return p;
    }
    
    void set_compatible_mean_and_var(const std::string& locus_desc, const std::pair<double, double>& p)
    {
        _compatible_mv_by_locus[locus_desc] = p;
    }
    
    void set_total_mean_and_var(const std::string& locus_desc, const std::pair<double, double>& p)
    {
        _total_mv_by_locus[locus_desc] = p;
    }
    
    const std::map<std::string, std::pair<double, double> >& total_mv_by_locus() const { return _total_mv_by_locus; }
    const std::map<std::string, std::pair<double, double> >& compatible_mv_by_locus() const { return _compatible_mv_by_locus; }
    
private:
    std::string         _name;
    std::vector<double> _scaled_compatible_mass_means;
    std::vector<double> _scaled_compatible_variances;
    std::vector<double> _scaled_mass_variances;
    
    std::map<std::string, std::pair<double, double> > _compatible_mv_by_locus;
    std::map<std::string, std::pair<double, double> > _total_mv_by_locus;
};

class PoissonDispersionModel : public MassDispersionModel
{
    std::string         _name;
    
public:
    PoissonDispersionModel(const std::string& name) : _name(name) {}
        
        virtual const std::string& name() const { return _name; }
    virtual double scale_mass_variance(double scaled_mass) const 
    { 
        return scaled_mass; 
    }
};



class MleErrorModel
{
public:
    MleErrorModel() {}
    MleErrorModel(const std::string& name,
                        const std::vector<double>& scaled_compatible_mass_means,
                        const std::vector<double>& scaled_mle_variances);
    
    virtual const std::string& name() const { return _name; }
    
    virtual double scale_mle_variance(double scaled_mass) const;
    
//    const vector<double>& scaled_compatible_mass_means() const { return _scaled_compatible_mass_means; }
//    const vector<double>& scaled_compatible_variances() const { return _scaled_compatible_variances; }
//    const vector<double>& scaled_mass_variances() const { return _scaled_mass_variances; }
    
private:
    std::string         _name;
    std::vector<double> _scaled_compatible_mass_means;
    std::vector<double> _scaled_mle_variances;
};

struct LocusCountList
{
    LocusCountList(std::string ld, int num_reps, int nt, const std::vector<std::string>& gids, const std::vector<std::string>& gnms) :
    locus_desc(ld), counts(std::vector<double>(num_reps, 0)), num_transcripts(nt), gene_ids(gids), gene_short_names(gnms) {}
    std::string locus_desc;
    std::vector<double> counts;
    int num_transcripts;
    vector<std::string> gene_ids;
    vector<std::string> gene_short_names;
};

void transform_counts_to_common_scale(const vector<double>& scale_factors,
                                      vector<LocusCountList>& sample_count_table);

void calc_geometric_scaling_factors(const std::vector<LocusCountList>& sample_count_table,
                                    std::vector<double>& scale_factors);

void calc_classic_fpkm_scaling_factors(const std::vector<LocusCountList>& sample_count_table,
                                       std::vector<double>& scale_factors);

void calc_quartile_scaling_factors(const std::vector<LocusCountList>& sample_count_table,
                                        std::vector<double>& scale_factors);

void calc_tmm_scaling_factors(const std::vector<LocusCountList>& sample_count_table,
                                        std::vector<double>& scale_factors);


boost::shared_ptr<MassDispersionModel>
fit_dispersion_model(const string& condition_name,
                     const std::vector<double>& scale_factors,
                     const std::vector<LocusCountList>& sample_count_table);

void calculate_count_means_and_vars(const vector<LocusCountList>& sample_compatible_count_table,
                                    vector<pair<double, double> >& means_and_vars);

// This factory merges bundles in a requested locus from several replicates
class ReplicatedBundleFactory
{
public:
	ReplicatedBundleFactory(const std::vector<boost::shared_ptr<BundleFactory> >& factories, 
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
        BOOST_FOREACH (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            if (fac->bundles_remain())
                return true;
        }
        return false;
    }
    
	bool next_bundle(HitBundle& bundle_out, bool cache_bundle)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        std::vector<HitBundle*> bundles;
        
        bool non_empty_bundle = false;
        BOOST_FOREACH (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            bundles.push_back(new HitBundle());
            if (fac->next_bundle(*(bundles.back()), cache_bundle))
            {
                non_empty_bundle = true;
            }
        }
        
        int locus_id = -1;
        for (size_t i = 0; i < bundles.size(); ++i)
        {
            if (locus_id == -1)
                locus_id = bundles[i]->id();
            if (locus_id != bundles[i]->id())
            {
                fprintf(stderr, "Error: locus id mismatch!\n");
                exit(1);
            }
        }
        
        
//        if (non_empty_bundle == false)
//        {
//            bundle_out.id(locus_id);
//            BOOST_FOREACH (HitBundle* in_bundle, bundles)
//            {
//                in_bundle->ref_scaffolds().clear();
//                in_bundle->clear_hits();
//                delete in_bundle;
//            }
//            return false;
//        }
        
        for (size_t i = 1; i < bundles.size(); ++i)
        {
            const vector<boost::shared_ptr<Scaffold> >& s1 = bundles[i]->ref_scaffolds();
            const vector<boost::shared_ptr<Scaffold> >& s2 =  bundles[i-1]->ref_scaffolds();
            assert (s1.size() == s2.size());
            for (size_t j = 0; j < s1.size(); ++j)
            {
                assert (s1[j]->annotated_trans_id() == s2[j]->annotated_trans_id());
            }
        }
        
        double total_compatible_mass = 0.0;
        double total_raw_mass = 0.0;
        
        for (size_t i = 0; i < bundles.size(); ++i)
        {
            total_compatible_mass += bundles[i]->compatible_mass();
            total_raw_mass += bundles[i]->raw_mass();
        }
        
        // Merge the replicates into a combined bundle of hits.
        HitBundle::combine(bundles, bundle_out);
        
        bundle_out.compatible_mass(total_compatible_mass);
        bundle_out.add_raw_mass(total_raw_mass);
        
        bundle_out.id(locus_id);
        
        BOOST_FOREACH (HitBundle* in_bundle, bundles)
        {
            in_bundle->ref_scaffolds().clear();
            in_bundle->clear_hits();
            delete in_bundle;
        }
        return non_empty_bundle;
    }
	
	void reset() 
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        BOOST_FOREACH (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            fac->reset();
        }
    }
    
    void inspect_replicate_maps(int& min_len, int& max_len, IdToLocusMap& id_to_locus_map)
    {
        vector<LocusCountList> sample_compatible_count_table;
        vector<LocusCountList> sample_total_count_table;
        
        vector<double> sample_masses;
        
        for (size_t fac_idx = 0; fac_idx < _factories.size(); ++fac_idx)
        {
            boost::shared_ptr<BundleFactory> fac = _factories[fac_idx];        
            BadIntronTable bad_introns;
            
            vector<LocusCount> compatible_count_table;
            vector<LocusCount> total_count_table;
            
            inspect_map(fac, NULL, compatible_count_table, total_count_table, id_to_locus_map, false, false);
            
            boost::shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            
            assert (compatible_count_table.size() == total_count_table.size());
            
            for (size_t i = 0; i < compatible_count_table.size(); ++i)
            {
                LocusCount& c = compatible_count_table[i];
                double raw_count = c.count;
                
                if (i >= sample_compatible_count_table.size())
                {
                    LocusCountList locus_count(c.locus_desc, _factories.size(), c.num_transcripts, c.gene_ids, c.gene_short_names);
                    sample_compatible_count_table.push_back(locus_count);
                    sample_compatible_count_table.back().counts[0] = raw_count;
                    sample_total_count_table.push_back(locus_count);
                    sample_total_count_table.back().counts[0] = total_count_table[i].count;
                }
                else
                {
                    if (sample_compatible_count_table[i].locus_desc != c.locus_desc)
                    {
                        fprintf (stderr, "Error: bundle boundaries don't match across replicates!\n");
                        exit(1);
                    }
                    sample_compatible_count_table[i].counts[fac_idx] = raw_count;
                    sample_total_count_table[i].counts[fac_idx] = total_count_table[i].count;
                }
            }
            
            rg_props->raw_compatible_counts(compatible_count_table);
            rg_props->raw_total_counts(total_count_table);
            
            sample_masses.push_back(rg_props->total_map_mass());
			min_len = min(min_len, rg_props->frag_len_dist()->min());
			max_len = max(max_len, rg_props->frag_len_dist()->max());
        }
        
        /*
        vector<double> scale_factors(_factories.size(), 0.0);
        
        calc_scaling_factors(sample_compatible_count_table, scale_factors);
        
        for (size_t i = 0; i < scale_factors.size(); ++i)
        {
            boost::shared_ptr<ReadGroupProperties> rg_props = _factories[i]->read_group_properties();
            assert (scale_factors[i] != 0);
            rg_props->internal_scale_factor(scale_factors[i]);
        }
        
        transform_counts_to_common_scale(scale_factors, sample_compatible_count_table);
        transform_counts_to_common_scale(scale_factors, sample_total_count_table);
        
        for (size_t fac_idx = 0; fac_idx < _factories.size(); ++fac_idx)
        {
            boost::shared_ptr<ReadGroupProperties> rg_props = _factories[fac_idx]->read_group_properties();
            assert (scale_factors[fac_idx] != 0);
            vector<LocusCount> common_scaled_compatible_counts;
            vector<LocusCount> common_scaled_total_counts;
            for (size_t j = 0; j < sample_compatible_count_table.size(); ++j)
            {
                common_scaled_compatible_counts.push_back(LocusCount(sample_compatible_count_table[j].locus_desc, sample_compatible_count_table[j].counts[fac_idx], sample_compatible_count_table[j].num_transcripts));
            }
            for (size_t j = 0; j < sample_total_count_table.size(); ++j)
            {
                common_scaled_total_counts.push_back(LocusCount(sample_total_count_table[j].locus_desc, sample_total_count_table[j].counts[fac_idx], sample_total_count_table[j].num_transcripts));
            }
            rg_props->common_scale_compatible_counts(common_scaled_compatible_counts);
            rg_props->common_scale_total_counts(common_scaled_total_counts);
        }
        fit_dispersion_model();
        */
        
    }

    void fit_dispersion_model()
    {
        vector<LocusCountList> sample_compatible_count_table;
        vector<LocusCountList> sample_total_count_table;
        
        vector<double> scale_factors;
        
        for (size_t fac_idx = 0; fac_idx < _factories.size(); ++fac_idx)
        {
            boost::shared_ptr<BundleFactory> fac = _factories[fac_idx];        
            
            boost::shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            const vector<LocusCount>& compatible_count_table = rg_props->common_scale_compatible_counts();
            const vector<LocusCount>& total_count_table = rg_props->common_scale_total_counts();
            
            assert(compatible_count_table.size() == compatible_count_table.size());
            
            for (size_t i = 0; i < compatible_count_table.size(); ++i)
            {
                const LocusCount& c = compatible_count_table[i];
                double common_scale_compatible_count = c.count;
                double common_scale_total_count = total_count_table[i].count;
                
                if (i >= sample_compatible_count_table.size())
                {
                    LocusCountList locus_count(c.locus_desc, _factories.size(), c.num_transcripts, c.gene_ids, c.gene_short_names);
                    sample_compatible_count_table.push_back(locus_count);
                    sample_compatible_count_table.back().counts[0] = common_scale_compatible_count;
                    sample_total_count_table.push_back(locus_count);
                    sample_total_count_table.back().counts[0] = common_scale_total_count;
                }
                else
                {
                    if (sample_compatible_count_table[i].locus_desc != c.locus_desc)
                    {
                        fprintf (stderr, "Error: bundle boundaries don't match across replicates!\n");
                        exit(1);
                    }
                    sample_compatible_count_table[i].counts[fac_idx] = common_scale_compatible_count;
                    sample_total_count_table[i].counts[fac_idx] = common_scale_total_count;
                }
            }
            scale_factors.push_back(rg_props->internal_scale_factor());
        }
        
        boost::shared_ptr<MassDispersionModel> disperser;
        disperser = ::fit_dispersion_model(_condition_name, scale_factors, sample_compatible_count_table);
        
        vector<pair<double, double> > compatible_means_and_vars;
        calculate_count_means_and_vars(sample_compatible_count_table,
                                       compatible_means_and_vars);
        
        for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
        {
            const LocusCountList& p = sample_compatible_count_table[i];
            double mean = compatible_means_and_vars[i].first;
            double var = compatible_means_and_vars[i].second;
            disperser->set_compatible_mean_and_var(p.locus_desc, make_pair(mean, var));
            //labeled_mv_table[p.locus_desc] = make_pair(mean, var);
        }
        
        vector<pair<double, double> > total_means_and_vars;
        calculate_count_means_and_vars(sample_total_count_table,
                                       total_means_and_vars);
        
        for (size_t i = 0; i < sample_total_count_table.size(); ++i)
        {
            const LocusCountList& p = sample_total_count_table[i];
            double mean = total_means_and_vars[i].first;
            double var = total_means_and_vars[i].second;
            disperser->set_total_mean_and_var(p.locus_desc, make_pair(mean, var));
            //labeled_mv_table[p.locus_desc] = make_pair(mean, var);
        }

        
        BOOST_FOREACH (boost::shared_ptr<BundleFactory> fac, _factories)
        {
            boost::shared_ptr<ReadGroupProperties> rg_props = fac->read_group_properties();
            rg_props->mass_dispersion_model(disperser);
        }
    }
    
    
    // This function NEEDS to deep copy the ref_mRNAs, otherwise cuffdiff'd
    // samples will clobber each other
    void set_ref_rnas(const vector<boost::shared_ptr<Scaffold> >& mRNAs, bool deep_copy = true)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        BOOST_FOREACH(boost::shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_ref_rnas(mRNAs, deep_copy);
        }
    }
    
    void set_mask_rnas(const vector<boost::shared_ptr<Scaffold> >& mRNAs)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        BOOST_FOREACH(boost::shared_ptr<BundleFactory> fac, _factories)
        {
            fac->set_mask_rnas(mRNAs);
        }
    }
    
    int num_replicates() const { return _factories.size(); }
    
    void mass_dispersion_model(boost::shared_ptr<MassDispersionModel const> disperser)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        BOOST_FOREACH(boost::shared_ptr<BundleFactory>& fac, _factories)
        {
            fac->read_group_properties()->mass_dispersion_model(disperser);
        }
    }
    
    boost::shared_ptr<MassDispersionModel const> mass_dispersion_model() const
    {
        return _factories.front()->read_group_properties()->mass_dispersion_model();
    }
    
    void mle_error_model(boost::shared_ptr<MleErrorModel const> mle_model)
    {
#if ENABLE_THREADS
        boost::mutex::scoped_lock lock(_rep_factory_lock);
#endif
        BOOST_FOREACH(boost::shared_ptr<BundleFactory>& fac, _factories)
        {
            fac->read_group_properties()->mle_error_model(mle_model);
        }
    }
    
    boost::shared_ptr<MleErrorModel const> mle_error_model() const
    {
        return _factories.front()->read_group_properties()->mle_error_model();
    }
    
private:
	vector<boost::shared_ptr<BundleFactory> > _factories;
#if ENABLE_THREADS
    boost::mutex _rep_factory_lock;
#endif
    string _condition_name; 
};

void normalize_counts(std::vector<boost::shared_ptr<ReadGroupProperties> > & all_read_groups);
