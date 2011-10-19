//
//  replicates.cpp
//  cufflinks
//
//  Created by Cole Trapnell on 3/11/11.
//  Copyright 2011 Cole Trapnell. All rights reserved.
//

#include <boost/thread.hpp>

extern "C" {
#include "locfit/local.h"
}

#include "replicates.h"

#if ENABLE_THREADS	
boost::mutex _locfit_lock;
#endif

MassDispersionModel::MassDispersionModel(const std::vector<double>& scaled_mass_means, 
                                         const std::vector<double>& scaled_raw_variances,
                                         const std::vector<double>& scaled_mass_variances) 
{
    if (scaled_mass_means.size() != scaled_mass_variances.size())
    {
        fprintf (stderr, "Error: dispersion model table is malformed\n");
    }
    
    double last_val = 0;
    for (size_t i = 0; i < scaled_mass_means.size(); i++)
    {
        
        if (last_val > scaled_mass_means[i])
        {
            fprintf (stderr, "Error: DispersionModel input is malformed\n");
        }
        
        if ( i == 0 || last_val < scaled_mass_means[i])
        {
            _scaled_mass_means.push_back(scaled_mass_means[i]);
            _scaled_raw_variances.push_back(scaled_raw_variances[i]);
            _scaled_mass_variances.push_back(scaled_mass_variances[i]);
        }
        else
        {
            // skip this element if it's equal to what we've already seen
        }
        
        last_val = scaled_mass_means[i];
    }
}

double MassDispersionModel::scale_mass_variance(double scaled_mass) const 
{
    if (scaled_mass <= 0)
        return 0.0;
    
    if (_scaled_mass_means.size() < 2 || _scaled_mass_variances.size() < 2)
    {
        return scaled_mass; // revert to poisson.
    }
    if (scaled_mass > _scaled_mass_means.back())
    {
        // extrapolate to the right
        // off the right end
        double x1_mean = _scaled_mass_means[_scaled_mass_means.size()-2];
        double x2_mean = _scaled_mass_means[_scaled_mass_means.size()-1];
        
        double y1_var = _scaled_mass_variances[_scaled_mass_means.size()-2];
        double y2_var = _scaled_mass_variances[_scaled_mass_means.size()-1];
        double slope = 0.0;                
        if (x2_mean != x1_mean)
        {
            slope = (y2_var - y1_var) / (x2_mean-x1_mean);
        }
        else if (y1_var == y2_var)
        {
            assert (false); // should have a unique'd table
        }
        double mean_interp = _scaled_mass_variances[_scaled_mass_means.size()-1] -
        slope*(scaled_mass - _scaled_mass_means.size()-1);
        if (mean_interp < scaled_mass)
            mean_interp = scaled_mass;
        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    else if (scaled_mass < _scaled_mass_means.front())
    {
        // extrapolate to the left
        // off the left end?
        double x1_mean = _scaled_mass_means[0];
        double x2_mean = _scaled_mass_means[1];
        
        double y1_var = _scaled_mass_variances[0];
        double y2_var = _scaled_mass_variances[1];
        double slope = 0.0;                
        if (x2_mean != x1_mean)
        {
            slope = (y2_var - y1_var) / (x2_mean-x1_mean);
        }
        else if (y1_var == y2_var)
        {
            assert (false); // should have a unique'd table
        }
        double mean_interp = _scaled_mass_variances[0] - slope*(_scaled_mass_means[0] - scaled_mass);
        if (mean_interp < scaled_mass)
            mean_interp = scaled_mass;

        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    
    vector<double>::const_iterator lb;
    lb = lower_bound(_scaled_mass_means.begin(), 
                     _scaled_mass_means.end(), 
                     scaled_mass);
    if (lb < _scaled_mass_means.end())
    {
        int d = lb - _scaled_mass_means.begin();
        if (*lb == scaled_mass || lb == _scaled_mass_means.begin())
        {
            double var = _scaled_mass_variances[d];
            if (var < scaled_mass) // revert to poisson if underdispersed
                var = scaled_mass;
            assert (!isnan(var) && !isinf(var));
            return var;
        }
        
        
        //in between two points on the scale.
        d--;
        
        if (d < 0)
        {
            fprintf(stderr, "ARG d < 0, d = %d \n", d);
        }
        
        if (d >= _scaled_mass_means.size())
        {
            fprintf(stderr, "ARG d >= _scaled_mass_means.size(), d = %d\n", d);
        }
        if (d >= _scaled_mass_variances.size())
        {
            fprintf(stderr, "ARG d >= _scaled_mass_variances.size(), d = %d\n", d);
        }
        
        double x1_mean = _scaled_mass_means[d];
        double x2_mean = _scaled_mass_means[d + 1];
        
        double y1_var = _scaled_mass_variances[d];
        double y2_var = _scaled_mass_variances[d + 1];
        double slope = 0.0;                
        if (x2_mean != x1_mean)
        {
            slope = (y2_var - y1_var) / (x2_mean-x1_mean);
        }
        else if (y1_var == y2_var)
        {
            assert (false); // should have a unique'd table
        }
        double mean_interp = _scaled_mass_variances[d] + slope*(scaled_mass - _scaled_mass_means[d]);
        if (mean_interp < scaled_mass) // revert to poisson if underdispersed
            mean_interp = scaled_mass;
 
        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    else
    {
        assert (!isnan(scaled_mass) && !isinf(scaled_mass));
        return scaled_mass; // revert to poisson assumption
    }
}

void calc_scaling_factors(const vector<LocusCountList>& sample_count_table,
                          vector<double>& scale_factors)
{
    vector<double> geom_means(sample_count_table.size(), 0.0);
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_count_table[i];
        
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            //assert (geom_means.size() > j);
            if (geom_means[i] > 0  && p.counts[j] > 0)
            {
                geom_means[i] *= p.counts[j];
            }
            else if (p.counts[j] > 0)
            {
                geom_means[i] = p.counts[j];
            }
            
        }
        geom_means[i] = pow(geom_means[i], 1.0/(double)p.counts.size());
    }
    
    for (size_t j = 0; j < scale_factors.size(); ++j)
    {
        vector<double> tmp_counts;
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            if (geom_means[i] && !isinf(geom_means[i]) && !isnan(geom_means[i]) && sample_count_table[i].counts[j])
            {
                double gm = (double)sample_count_table[i].counts[j] / geom_means[i];
                assert (!isinf(gm));
                tmp_counts.push_back(gm);
            }
        }
        sort(tmp_counts.begin(), tmp_counts.end());
        if (!tmp_counts.empty())
            scale_factors[j] = tmp_counts[tmp_counts.size()/2];
        else
            scale_factors[j] = 1.0;
    }
}

static const int min_loci_for_fitting = 30;

boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model_helper(const string& condition_name,
                            const vector<double>& scale_factors,
                            const vector<LocusCountList>& sample_count_table)
{
    vector<pair<double, double> > raw_means_and_vars;
    map<string, pair<double, double> > labeled_mv_table;
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_count_table[i];
        double mean = accumulate(p.counts.begin(), p.counts.end(), 0.0);
        if (mean > 0.0 && p.counts.size() > 0)
            mean /= p.counts.size();
        
        double var = 0.0;
        foreach (double d, p.counts)
        {
            var += (d - mean) * (d - mean);
        }
        if (var > 0.0 && p.counts.size())
            var /= p.counts.size();
        labeled_mv_table[p.locus_desc] = make_pair(mean, var);
        if (mean > 0 && var > 0.0)
        {
            //fprintf(stderr, "%s\t%lg\t%lg\n", p.locus_desc.c_str(), mean, var);
            raw_means_and_vars.push_back(make_pair(mean, var));
        }
    }
    
    if (raw_means_and_vars.size() < min_loci_for_fitting)
    {
        shared_ptr<MassDispersionModel> disperser;
        disperser = shared_ptr<MassDispersionModel>(new PoissonDispersionModel);
        
        for (map<string, pair<double, double> >::iterator itr = labeled_mv_table.begin();
             itr != labeled_mv_table.end();
             ++itr)
        {
            string label = itr->first;
            pair<double, double> p = itr->second;
            disperser->set_raw_mean_and_var(itr->first, itr->second);
        }
        //fprintf(stderr, "Warning: fragment count variances between replicates are all zero, reverting to Poisson model\n");
        return disperser;
    }
    
    sort(raw_means_and_vars.begin(), raw_means_and_vars.end());
    
    vector<double> raw_means(raw_means_and_vars.size(), 0.0);
    vector<double> raw_variances(raw_means_and_vars.size(), 0.0);
    
    for(size_t i = 0; i < raw_means_and_vars.size(); ++i)
    {
        raw_means[i] = raw_means_and_vars[i].first;
        raw_variances[i] = raw_means_and_vars[i].second;
    }
    
    vector<double> fitted_values(raw_means_and_vars.size(), 0.0);
    
    setuplf();  
    
    // WARNING: locfit doesn't like undescores - need camel case for 
    // variable names
    
    char namebuf[256];
    sprintf(namebuf, "countMeans");
    vari* cm = createvar(namebuf,STREGULAR,raw_means.size(),VDOUBLE);
    for (size_t i = 0; i < raw_means.size(); ++i)
    {
        cm->dpr[i] = log(raw_means[i]);
    }
    
    sprintf(namebuf, "countVariances");
    vari* cv = createvar(namebuf,STREGULAR,raw_variances.size(),VDOUBLE);
    for (size_t i = 0; i < raw_variances.size(); ++i)
    {
        cv->dpr[i] = raw_variances[i]; 
    }
    
    char locfit_cmd[2048];
    sprintf(locfit_cmd, "locfit countVariances~countMeans family=gamma");
    
    locfit_dispatch(locfit_cmd);
    
    sprintf(locfit_cmd, "fittedVars=predict countMeans");
    locfit_dispatch(locfit_cmd);
    
    //sprintf(locfit_cmd, "prfit x fhat h nlx");
    //locfit_dispatch(locfit_cmd);
    
    int n = 0;
    sprintf(namebuf, "fittedVars");
    vari* cp = findvar(namebuf, 1, &n);
    assert(cp != NULL);
    for (size_t i = 0; i < cp->n; ++i)
    {
        fitted_values[i] = cp->dpr[i];
    }
    
    shared_ptr<MassDispersionModel> disperser;
    disperser = shared_ptr<MassDispersionModel>(new MassDispersionModel(raw_means, raw_variances, fitted_values));
    if (poisson_dispersion)
        disperser = shared_ptr<MassDispersionModel>(new PoissonDispersionModel);
    
    for (map<string, pair<double, double> >::iterator itr = labeled_mv_table.begin();
         itr != labeled_mv_table.end();
         ++itr)
    {
        string label = itr->first;
        pair<double, double> p = itr->second;
        disperser->set_raw_mean_and_var(itr->first, itr->second);
    }
    
    //pair<double, double> p = disperser->get_raw_mean_and_var("chr1:11873-29961");
    
    return disperser;
}

boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model(const string& condition_name,
                     const vector<double>& scale_factors,
                     const vector<LocusCountList>& sample_count_table)
{
//    
//#if ENABLE_THREADS
//	boost::mutex::scoped_lock lock(_locfit_lock);
//#endif
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        if (sample_count_table[i].counts.size() <= 1)
        {
            // only one replicate - no point in fitting variance
            return shared_ptr<MassDispersionModel const>(new PoissonDispersionModel);
        }
    }
#if ENABLE_THREADS
    _locfit_lock.lock();
#endif
    
    ProgressBar p_bar("Modeling fragment count overdispersion.",0);
    
    int max_transcripts = 0;
    foreach(const LocusCountList& L, sample_count_table)
    {
        if (L.num_transcripts > max_transcripts)
        {
            max_transcripts = L.num_transcripts;
        }
    }
    
    // This vector holds a dispersion model for each transcript multiplicity. 
    // The model for multiplicity is fitted to all the data, and lives at index 0 in the 
    // vector below.
    vector<boost::shared_ptr<MassDispersionModel const> > disp_models(max_transcripts+1);

    for (size_t i = 0; i < max_transcripts; i++)
    {
        boost::shared_ptr<MassDispersionModel const> model;
        if (i != 0)
        {
//            vector<LocusCountList> sample_count_subtable;
//            foreach(const LocusCountList& L, sample_count_table)
//            {
//                if (L.num_transcripts == i)
//                {
//                    sample_count_subtable.push_back(L);
//                }
//            }
//            model = fit_dispersion_model_helper(condition_name, scale_factors, sample_count_subtable);
        }
        else
        {
            model = fit_dispersion_model_helper(condition_name, scale_factors, sample_count_table);
        }
        disp_models[i] = model;
    }
    
    if (emit_count_tables)
    {
//        string cond_count_filename = output_dir + "/" + condition_name + "_counts.txt";
//        
//        FILE* sample_count_file = fopen(cond_count_filename.c_str(), "w");
//        
//        if (sample_count_file)
//        {
//            fprintf(sample_count_file, "count_mean\tcount_var\tfitted_var\tnum_transcripts\n");
//            for (size_t j = 0; j < max_transcripts; j++)
//            {
//                boost::shared_ptr<MassDispersionModel const> model = disp_models[j];
//                const vector<double>& means = model->scaled_mass_means();
//                const vector<double>& raw_vars  = model->scaled_raw_variances();
//                
//                for (size_t i = 0; i < means.size(); ++i)
//                {
//                    fprintf(sample_count_file, "%lg\t%lg\t%lg\t%lu\n", 
//                            means[i], 
//                            raw_vars[i],
//                            model->scale_mass_variance(means[i]),
//                            j);
//                }
//            }
//            fclose(sample_count_file);
//        } 

        string cond_count_filename = output_dir + "/" + condition_name + "_counts.txt";
        
        FILE* sample_count_file = fopen(cond_count_filename.c_str(), "w");
        
        if (sample_count_file)
        {
            fprintf(sample_count_file, "count_mean\tcount_var\tfitted_var\n");
            
            boost::shared_ptr<MassDispersionModel const> model = disp_models[0];
            const vector<double>& means = model->scaled_mass_means();
            const vector<double>& raw_vars  = model->scaled_raw_variances();
            
            for (size_t i = 0; i < means.size(); ++i)
            {
                fprintf(sample_count_file, "%lg\t%lg\t%lg\n", 
                        means[i], 
                        raw_vars[i],
                        model->scale_mass_variance(means[i]));
            }
            fclose(sample_count_file);
        } 
    }

#if ENABLE_THREADS
    _locfit_lock.unlock();
#endif
    return disp_models[0];
}
