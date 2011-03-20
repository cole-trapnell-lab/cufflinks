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
    if (lb != _scaled_mass_means.end())
    {
        int d = lb - _scaled_mass_means.begin();
        if (*lb == scaled_mass)
        {
            assert (!isnan(_scaled_mass_variances[d]) && !isinf(_scaled_mass_variances[d]));
            return _scaled_mass_variances[d];
        }
        
        //in between two points on the scale.
        d--;
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
        if (mean_interp < scaled_mass)
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

void calc_scaling_factors(const vector<pair<string, vector<double> > >& sample_count_table,
                          vector<double>& scale_factors)
{
    vector<double> geom_means(sample_count_table.size(), 0.0);
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const pair<string, vector<double> >& p = sample_count_table[i];
        
        for (size_t j = 0; j < p.second.size(); ++j)
        {
            assert (geom_means.size() > j);
            if (geom_means[i] > 0  && p.second[j] > 0)
            {
                geom_means[i] *= p.second[j];
            }
            else if (p.second[j] > 0)
            {
                geom_means[i] = p.second[j];
            }
        }
        geom_means[i] = pow(geom_means[i], 1.0/(double)p.second.size());
    }
    
    for (size_t j = 0; j < scale_factors.size(); ++j)
    {
        vector<double> tmp_counts;
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            if (geom_means[i] && !isinf(geom_means[i]) && !isnan(geom_means[i]) && sample_count_table[i].second[j])
            {
                double gm = (double)sample_count_table[i].second[j] / geom_means[i];
                assert (!isinf(gm));
                tmp_counts.push_back(gm);
            }
        }
        sort(tmp_counts.begin(), tmp_counts.end());
        if (!tmp_counts.empty())
            scale_factors[j] = tmp_counts[tmp_counts.size()/2];
        else
            scale_factors[j] = 0.0;
    }
}

boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model(const vector<double>& scale_factors,
                     const vector<pair<string, vector<double> > >& sample_count_table)
{
//    
//#if ENABLE_THREADS
//	boost::mutex::scoped_lock lock(_locfit_lock);
//#endif
    
    vector<pair<double, double> > raw_means_and_vars;
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        if (sample_count_table[i].second.size() <= 1)
        {
            // only one replicate - no point in fitting variance
            return shared_ptr<MassDispersionModel const>(new PoissonDispersionModel);
        }
    }
    
    _locfit_lock.lock();
    
    ProgressBar p_bar("Modeling fragment count overdispersion.",0);
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const pair<string, vector<double> >& p = sample_count_table[i];
        double mean = accumulate(p.second.begin(), p.second.end(), 0.0);
        if (mean > 0.0 && p.second.size() > 0)
            mean /= p.second.size();
        
        double var = 0.0;
        foreach (double d, p.second)
        {
            var += (d - mean) * (d - mean);
        }
        if (var > 0.0 && p.second.size())
            var /= p.second.size();
        if (mean > 0 && var > 0.0)
            raw_means_and_vars.push_back(make_pair(mean, var));
    }
    
    if (raw_means_and_vars.empty())
    {
        fprintf(stderr, "Warning: fragment count variances between replicates are all zero, reverting to Poisson model\n");
        _locfit_lock.unlock();
        return shared_ptr<MassDispersionModel const>(new PoissonDispersionModel);
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
    disperser = shared_ptr<MassDispersionModel>(new MassDispersionModel(raw_means, fitted_values));
    if (poisson_dispersion)
        disperser = shared_ptr<MassDispersionModel>(new PoissonDispersionModel);

    char sample_name_buf[256];
    int sample_id = rand();
    sprintf(sample_name_buf, "%d_counts.txt", sample_id);
    FILE* sample_count_file = fopen(sample_name_buf, "w");
    
    if (sample_count_file)
    {
        fprintf(sample_count_file, "count_mean\tcount_var\tfitted_var\n");
        for (size_t i = 0; i < raw_means_and_vars.size(); ++i)
        {
            fprintf(sample_count_file, "%lg\t%lg\t%lg\n", 
                    raw_means_and_vars[i].first, 
                    raw_means_and_vars[i].second,
                    fitted_values[i]);
        }
        fclose(sample_count_file);
    } 
 
    _locfit_lock.unlock();
    
    return disperser;
}
