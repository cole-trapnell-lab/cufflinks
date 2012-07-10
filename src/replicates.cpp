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

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>

#include "replicates.h"
//#include "rounding.h"

#if ENABLE_THREADS	
boost::mutex _locfit_lock;
#endif

MassDispersionModel::MassDispersionModel(const std::string& name,
                                         const std::vector<double>& scaled_mass_means, 
                                         const std::vector<double>& scaled_raw_variances,
                                         const std::vector<double>& scaled_mass_variances) 
{
    _name = name;
    
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

void transform_counts_to_common_scale(const vector<double>& scale_factors,
                                      vector<LocusCountList>& sample_count_table)
{
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
}

void calc_scaling_factors(const vector<LocusCountList>& sample_count_table,
                          vector<double>& scale_factors)
{
    
    vector<double> log_geom_means(sample_count_table.size(), 0.0);
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_count_table[i];
        
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            //assert (log_geom_means.size() > j);
            if (floor(p.counts[j]) > 0)
            {
                log_geom_means[i] += (1.0/p.counts.size()) * log(floor(p.counts[j]));
            }
            
        }
        //log_geom_means[i] = pow(log_geom_means[i], 1.0/(double)p.counts.size());
    }
    
    for (size_t j = 0; j < scale_factors.size(); ++j)
    {
        vector<double> tmp_counts;
        for (size_t i = 0; i < sample_count_table.size(); ++i)
        {
            if (log_geom_means[i] && !isinf(log_geom_means[i]) && !isnan(log_geom_means[i]) && floor(sample_count_table[i].counts[j]))
            {
                double gm = (double)log(floor(sample_count_table[i].counts[j])) - log_geom_means[i];
                assert (!isinf(gm));
                tmp_counts.push_back(gm);
            }
        }
        sort(tmp_counts.begin(), tmp_counts.end());
        if (!tmp_counts.empty())
            scale_factors[j] = exp(tmp_counts[tmp_counts.size()/2]);
        else
            scale_factors[j] = 1.0;
    }
}

static const int min_loci_for_fitting = 30;

struct SCVInterpolator
{
    void add_scv_pair(double est_scv, double true_scv)
    {
        true_scvs.push_back(true_scv);
        est_scvs.push_back(est_scv);
    }
    
    void finalize() 
    {
        vector<pair<double, double> > pv;
        for (size_t i =0; i < true_scvs.size(); ++i)
        {
            pv.push_back(make_pair(est_scvs[i], true_scvs[i]));
        }
        sort(pv.begin(), pv.end());
    }
    
    // This was built from the dispersion model interpolator - we should refactor these 
    // into a single routine.
    double interpolate_scv(double est_scv)
    {
        if (est_scv <= 0)
            return 0.0;
        
        if (est_scvs.size() < 2 || true_scvs.size() < 2)
        {
            return est_scv; // revert to poisson.
        }
        if (est_scv > est_scvs.back())
        {
            // extrapolate to the right
            // off the right end
            double x1_mean = est_scvs[est_scvs.size()-2];
            double x2_mean = est_scvs[est_scvs.size()-1];
            
            double y1_var = true_scvs[est_scvs.size()-2];
            double y2_var = true_scvs[est_scvs.size()-1];
            double slope = 0.0;                
            if (x2_mean != x1_mean)
            {
                slope = (y2_var - y1_var) / (x2_mean-x1_mean);
            }
            else if (y1_var == y2_var)
            {
                assert (false); // should have a unique'd table
            }
            double mean_interp = true_scvs[est_scvs.size()-1] -
            slope*(est_scv - est_scvs.size()-1);
            if (mean_interp < est_scv)
                mean_interp = est_scv;
            assert (!isnan(mean_interp) && !isinf(mean_interp));
            return mean_interp;
        }
        else if (est_scv < est_scvs.front())
        {
            // extrapolate to the left
            // off the left end?
            double x1_mean = est_scvs[0];
            double x2_mean = est_scvs[1];
            
            double y1_var = true_scvs[0];
            double y2_var = true_scvs[1];
            double slope = 0.0;                
            if (x2_mean != x1_mean)
            {
                slope = (y2_var - y1_var) / (x2_mean-x1_mean);
            }
            else if (y1_var == y2_var)
            {
                assert (false); // should have a unique'd table
            }
            double mean_interp = true_scvs[0] - slope*(est_scvs[0] - est_scv);
            if (mean_interp < est_scv)
                mean_interp = est_scv;
            
            assert (!isnan(mean_interp) && !isinf(mean_interp));
            return mean_interp;
        }
        
        vector<double>::const_iterator lb;
        lb = lower_bound(est_scvs.begin(), 
                         est_scvs.end(), 
                         est_scv);
        if (lb < est_scvs.end())
        {
            int d = lb - est_scvs.begin();
            if (*lb == est_scv || lb == est_scvs.begin())
            {
                double var = true_scvs[d];
                if (var < est_scv) // revert to poisson if underdispersed
                    var = est_scv;
                assert (!isnan(var) && !isinf(var));
                return var;
            }
            
            
            //in between two points on the scale.
            d--;
            
            if (d < 0)
            {
                fprintf(stderr, "ARG d < 0, d = %d \n", d);
            }
            
            if (d >= est_scvs.size())
            {
                fprintf(stderr, "ARG d >= est_scvs.size(), d = %d\n", d);
            }
            if (d >= true_scvs.size())
            {
                fprintf(stderr, "ARG d >= true_scvs.size(), d = %d\n", d);
            }
            
            double x1_mean = est_scvs[d];
            double x2_mean = est_scvs[d + 1];
            
            double y1_var = true_scvs[d];
            double y2_var = true_scvs[d + 1];
            double slope = 0.0;                
            if (x2_mean != x1_mean)
            {
                slope = (y2_var - y1_var) / (x2_mean-x1_mean);
            }
            else if (y1_var == y2_var)
            {
                assert (false); // should have a unique'd table
            }
            double mean_interp = true_scvs[d] + slope*(est_scv - est_scvs[d]);
            if (mean_interp < est_scv) // revert to poisson if underdispersed
                mean_interp = est_scv;
            
            assert (!isnan(mean_interp) && !isinf(mean_interp));
            return mean_interp;
        }
        else
        {
            assert (!isnan(est_scv) && !isinf(est_scv));
            return est_scv; // revert to poisson assumption
        }
        
        return 0.0;
    }
    
private:
    vector<double> true_scvs;
    vector<double> est_scvs;
};

void build_scv_correction_fit(int nreps, int ngenes, int mean_count, SCVInterpolator& true_to_est_scv_table)
{    
    setuplf();  
    
    vector<pair<double, double> > alpha_vs_scv;
   
    boost::mt19937 rng;
    vector<boost::random::negative_binomial_distribution<int, double> > nb_gens;
    
    vector<double> alpha_range;
    for(double a = 0.02; a < 2.0; a += (2.0 / 100.0))
    {
        alpha_range.push_back(a);
    }
    
    for (double a = 2; a < 10.0; a += (8.0 / 20.0))
    {
        alpha_range.push_back(a);
    }
    
    foreach (double alpha, alpha_range)
    {
        double k = 1.0/alpha;
        double p = k / (k + mean_count);
        double r = (mean_count * p) / (1-p);

        //boost::random::negative_binomial_distribution<int, double> nb(r, p);
        
        boost::random::gamma_distribution<double> gamma(r, (1-p)/p);
        
        
        vector<double> scvs_for_alpha;
        vector<double> draws;
        for (size_t i = 0; i < ngenes; ++i)
        {
            LocusCountList locus_count("", nreps, 1); 
            for (size_t rep_idx = 0; rep_idx < nreps; ++rep_idx)
            {
                boost::random::poisson_distribution<long, double> poisson(gamma(rng));
                locus_count.counts[rep_idx] = poisson(rng);
                draws.push_back(locus_count.counts[rep_idx]);
                //fprintf(stderr, "%lg\t", locus_count.counts[rep_idx]);
            }
            
            double mean = accumulate(locus_count.counts.begin(), locus_count.counts.end(), 0);
            if (mean == 0)
                continue;
            mean /= locus_count.counts.size();
            double var = 0.0;
            foreach(double c,  locus_count.counts)
            {
                var += (c-mean)*(c-mean);
            }
            var /= locus_count.counts.size();
            var *= locus_count.counts.size() / (locus_count.counts.size() - 1);
            
            double scv = var / (mean*mean);
            scvs_for_alpha.push_back(scv);
            //fprintf(stderr, " : mean = %lg, var = %lg, scv = %lg\n", mean, var, scv);
            
            //fprintf(stderr, "\n");
        }
        
        double mean = accumulate(draws.begin(), draws.end(), 0);
        mean /= draws.size();
        double var = 0.0;
        foreach(int c,  draws)
        {
            var += (c - mean)*(c-mean);
        }
        var /= draws.size();
        var *= draws.size() / (draws.size() - 1);
        
       
        
        //fprintf(stderr, "##########\n");
        //fprintf(stderr, "mean = %lf, var = %lg\n", mean, var);
        if (scvs_for_alpha.size() > 0)
        {
            double mean_scv = accumulate(scvs_for_alpha.begin(),scvs_for_alpha.end(), 0.0);
            mean_scv /= scvs_for_alpha.size();
            //fprintf(stderr, "alpha = %lg scv = %lg\n", alpha, mean_scv);
            alpha_vs_scv.push_back(make_pair(alpha, mean_scv));
        }
    }
    
    //fprintf(stderr, "$$$$$$$$$\n");
    
    //sort (alpha_range.begin(), alpha_range.end());
    
    char namebuf[256];
    sprintf(namebuf, "trueSCV");
    vari* cm = createvar(namebuf,STREGULAR,alpha_vs_scv.size(),VDOUBLE);
    for (size_t i = 0; i < alpha_vs_scv.size(); ++i)
    {
        cm->dpr[i] = alpha_vs_scv[i].first;
    }
    
    sprintf(namebuf, "estSCV");
    vari* cv = createvar(namebuf,STREGULAR,alpha_vs_scv.size(),VDOUBLE);
    for (size_t i = 0; i < alpha_vs_scv.size(); ++i)
    {
        cv->dpr[i] = alpha_vs_scv[i].second;
    }
    
    char locfit_cmd[2048];
    sprintf(locfit_cmd, "locfit trueSCV~estSCV");
    
    locfit_dispatch(locfit_cmd);
    
    sprintf(namebuf, "domainSCV");
    vari* cd = createvar(namebuf,STREGULAR,alpha_vs_scv.size(),VDOUBLE);
    for (size_t i = 0; i < alpha_vs_scv.size(); ++i)
    {
        cd->dpr[i] = alpha_vs_scv[i].second;
    }
    
    sprintf(locfit_cmd, "fittedSCV=predict domainSCV");
    locfit_dispatch(locfit_cmd);
    
    int n = 0;
    sprintf(namebuf, "fittedSCV");
    vari* cp = findvar(namebuf, 1, &n);
    assert(cp != NULL);
    
    for (size_t i = 0; i < cp->n; ++i)
    {
        //fprintf(stderr, "%lg\t%lg\n",alpha_range[i], cp->dpr[i]);
        true_to_est_scv_table.add_scv_pair(alpha_range[i], cp->dpr[i]);
    }
}

boost::shared_ptr<MassDispersionModel const> 
fit_dispersion_model_helper(const string& condition_name,
                            const vector<double>& scale_factors,
                            const vector<LocusCountList>& sample_count_table,
                            bool exclude_zero_samples)
{
    vector<pair<double, double> > raw_means_and_vars;
    map<string, pair<double, double> > labeled_mv_table;
    
    SCVInterpolator true_to_est_scv_table;
    
    build_scv_correction_fit(scale_factors.size(), 10000, 100000, true_to_est_scv_table);
    
    setuplf();  
    
    double xim = 0;
    foreach(double s, scale_factors)
    {
        if (s)
            xim += 1.0 / s;
    }
    xim /= scale_factors.size();
    
    for (size_t i = 0; i < sample_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_count_table[i];
        double mean = accumulate(p.counts.begin(), p.counts.end(), 0.0);
        if (mean > 0.0 && p.counts.size() > 0)
            mean /= p.counts.size();
        
        double var = 0.0;
        double num_non_zero = 0;
        foreach (double d, p.counts)
        {
            if (d > 0)
                num_non_zero++;
            var += (d - mean) * (d - mean);
        }
        if (var > 0.0 && p.counts.size())
        {
            var /= p.counts.size();
            var *= p.counts.size() / (p.counts.size() - 1);
        }
        labeled_mv_table[p.locus_desc] = make_pair(mean, var);
        //double scv = (var - xim * mean) / (mean * var);
        if (mean > 0 && var > 0.0 && /*scv > 0 && */ (!exclude_zero_samples || num_non_zero == p.counts.size()))
        {
            //fprintf(stderr, "%s\t%lg\t%lg\n", p.locus_desc.c_str(), mean, var);
            raw_means_and_vars.push_back(make_pair(mean, var));
        }
    }
    
    if (raw_means_and_vars.size() < min_loci_for_fitting)
    {
        shared_ptr<MassDispersionModel> disperser;
        disperser = shared_ptr<MassDispersionModel>(new PoissonDispersionModel(condition_name));
        
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
    vector<double> raw_scvs(raw_means_and_vars.size(), 0.0);
    
    for(size_t i = 0; i < raw_means_and_vars.size(); ++i)
    {
        raw_means[i] = raw_means_and_vars[i].first;
        raw_variances[i] = raw_means_and_vars[i].second;
        raw_scvs[i] = (raw_variances[i] - xim * raw_means[i]) / (raw_means[i] * raw_means[i]);
    }
    
    vector<double> fitted_values(raw_means_and_vars.size(), 0.0);
    
    // WARNING: locfit doesn't like undescores - need camel case for 
    // variable names
    
    char namebuf[256];
    sprintf(namebuf, "countMeans");
    vari* cm = createvar(namebuf,STREGULAR,raw_means.size(),VDOUBLE);
    for (size_t i = 0; i < raw_means.size(); ++i)
    {
        cm->dpr[i] = log(raw_means[i]);
    }
    
    //sprintf(namebuf, "countSCV");
    sprintf(namebuf, "countVariances");
    vari* cv = createvar(namebuf,STREGULAR,raw_variances.size(),VDOUBLE);
    for (size_t i = 0; i < raw_variances.size(); ++i)
    {
        cv->dpr[i] = raw_variances[i]; 
        //cv->dpr[i] = raw_scvs[i];
    }
    
    char locfit_cmd[2048];
    //sprintf(locfit_cmd, "locfit countVariances~countMeans family=gamma");
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
        // var[NB] = mu + mu^2 / (fitted_disp^-1)
        // var[NB] = mu + mu^2 * fitted_disp
        //double fitted_dispersion = 0;
        //if (cp->dpr[i] >0)
        //    fitted_dispersion = 1.0 / fitted_dispersion;
        if (cp->dpr[i] >= 0)
        {
            //fprintf (stderr, "%lg\t%lg\n", raw_means[i], cp->dpr[i]);
            double mean = exp(cm->dpr[i]);
            double fitted_scv = (cp->dpr[i]) / (mean * mean);
            double corrected_scv = true_to_est_scv_table.interpolate_scv(fitted_scv);
            //double corrected_scv = cp->dpr[i];
            
            // uncorrected fitted variance:
            fitted_values[i] = mean + (cp->dpr[i] - xim * mean);
            
//            double k = 1.0/corrected_scv;
//            double p = k / (k + 100000);
//            double r = (100000 * p) / (1-p);
//            
//            double hypothetical_mean = p*r/ (1-p);
//            double hypothetical_var = p*r/((1-p)*(1-p));
//            double raw_var = hypothetical_var - hypothetical_mean;
            // bias corrected fitted_variance:
            //fitted_values[i] = mean + raw_var;
            //fprintf(stderr, "mean = %lg, variance = %lg, uncorrected scv = %lg, corrected scv = %lg, var = %lg\n", mean, cp->dpr[i], fitted_scv, corrected_scv, fitted_values[i]);
        }
        else
        {
            fitted_values[i] = raw_means[i];
        }
    }
    
    shared_ptr<MassDispersionModel> disperser;
    disperser = shared_ptr<MassDispersionModel>(new MassDispersionModel(condition_name, raw_means, raw_variances, fitted_values));
    if (poisson_dispersion)
        disperser = shared_ptr<MassDispersionModel>(new PoissonDispersionModel(condition_name));
    
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
                     const vector<LocusCountList>& sample_count_table,
                     bool exclude_zero_samples)
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
            return shared_ptr<MassDispersionModel const>(new PoissonDispersionModel(condition_name));
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
            model = fit_dispersion_model_helper(condition_name, scale_factors, sample_count_table, exclude_zero_samples);
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
