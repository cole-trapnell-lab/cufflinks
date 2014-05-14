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
                                         const std::vector<double>& scaled_compatible_mass_means, 
                                         const std::vector<double>& scaled_compatible_variances,
                                         const std::vector<double>& scaled_mass_variances) 
{
    _name = name;
    
    if (scaled_compatible_mass_means.size() != scaled_mass_variances.size())
    {
        fprintf (stderr, "Error: dispersion model table is malformed\n");
    }
    
    double last_val = 0;
    for (size_t i = 0; i < scaled_compatible_mass_means.size(); i++)
    {
        
        if (last_val > scaled_compatible_mass_means[i])
        {
            fprintf (stderr, "Error: DispersionModel input is malformed\n");
        }
        
        if ( i == 0 || last_val < scaled_compatible_mass_means[i])
        {
            _scaled_compatible_mass_means.push_back(scaled_compatible_mass_means[i]);
            _scaled_compatible_variances.push_back(scaled_compatible_variances[i]);
            _scaled_mass_variances.push_back(scaled_mass_variances[i]);
        }
        else
        {
            // skip this element if it's equal to what we've already seen
        }
        
        last_val = scaled_compatible_mass_means[i];
    }
}

double MassDispersionModel::scale_mass_variance(double scaled_mass) const 
{
    if (scaled_mass <= 0)
        return 0.0;
    
    if (_scaled_compatible_mass_means.size() < 2 || _scaled_mass_variances.size() < 2)
    {
        return scaled_mass; // revert to poisson.
    }
    if (scaled_mass > _scaled_compatible_mass_means.back())
    {
        // extrapolate to the right
        // off the right end
        double x1_mean = _scaled_compatible_mass_means[_scaled_compatible_mass_means.size()-2];
        double x2_mean = _scaled_compatible_mass_means[_scaled_compatible_mass_means.size()-1];
        
        double y1_var = _scaled_mass_variances[_scaled_compatible_mass_means.size()-2];
        double y2_var = _scaled_mass_variances[_scaled_compatible_mass_means.size()-1];
        double slope = 0.0;                
        if (x2_mean != x1_mean)
        {
            slope = (y2_var - y1_var) / (x2_mean-x1_mean);
        }
        else if (y1_var == y2_var)
        {
            assert (false); // should have a unique'd table
        }
        double mean_interp = _scaled_mass_variances[_scaled_compatible_mass_means.size()-1] -
        slope*(scaled_mass - _scaled_compatible_mass_means.size()-1);
        if (mean_interp < scaled_mass)
            mean_interp = scaled_mass;
        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    else if (scaled_mass < _scaled_compatible_mass_means.front())
    {
        // extrapolate to the left
        // off the left end?
        double x1_mean = _scaled_compatible_mass_means[0];
        double x2_mean = _scaled_compatible_mass_means[1];
        
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
        double mean_interp = _scaled_mass_variances[0] - slope*(_scaled_compatible_mass_means[0] - scaled_mass);
        if (mean_interp < scaled_mass)
            mean_interp = scaled_mass;

        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    
    vector<double>::const_iterator lb;
    lb = lower_bound(_scaled_compatible_mass_means.begin(), 
                     _scaled_compatible_mass_means.end(), 
                     scaled_mass);
    if (lb < _scaled_compatible_mass_means.end())
    {
        int d = lb - _scaled_compatible_mass_means.begin();
        if (*lb == scaled_mass || lb == _scaled_compatible_mass_means.begin())
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
        
        if (d >= _scaled_compatible_mass_means.size())
        {
            fprintf(stderr, "ARG d >= _scaled_compatible_mass_means.size(), d = %d\n", d);
        }
        if (d >= _scaled_mass_variances.size())
        {
            fprintf(stderr, "ARG d >= _scaled_mass_variances.size(), d = %d\n", d);
        }
        
        double x1_mean = _scaled_compatible_mass_means[d];
        double x2_mean = _scaled_compatible_mass_means[d + 1];
        
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
        double mean_interp = _scaled_mass_variances[d] + slope*(scaled_mass - _scaled_compatible_mass_means[d]);
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


MleErrorModel::MleErrorModel(const std::string& name,
                             const std::vector<double>& scaled_compatible_mass_means,
                             const std::vector<double>& scaled_mle_variances)
{
    _name = name;
    
    if (scaled_compatible_mass_means.size() != scaled_mle_variances.size())
    {
        fprintf (stderr, "Error: dispersion model table is malformed\n");
    }
    
    double last_val = 0;
    for (size_t i = 0; i < scaled_compatible_mass_means.size(); i++)
    {
        
        if (last_val > scaled_compatible_mass_means[i])
        {
            fprintf (stderr, "Error: MLEError input is malformed\n");
        }
        
        if ( i == 0 || last_val < scaled_compatible_mass_means[i])
        {
            _scaled_compatible_mass_means.push_back(scaled_compatible_mass_means[i]);
            _scaled_mle_variances.push_back(scaled_mle_variances[i]);
        }
        else
        {
            // skip this element if it's equal to what we've already seen
        }
        
        last_val = scaled_compatible_mass_means[i];
    }
}

double MleErrorModel::scale_mle_variance(double scaled_mass) const
{
    if (scaled_mass <= 0)
        return 0.0;
    
    if (_scaled_compatible_mass_means.size() < 2 || _scaled_mle_variances.size() < 2)
    {
        return 0; // revert to poisson.
    }
    if (scaled_mass > _scaled_compatible_mass_means.back())
    {
        return 0;
    }
    else if (scaled_mass < _scaled_compatible_mass_means.front())
    {
        return 0; // we won't add anything if we're out of the range of the table
    }
    
    vector<double>::const_iterator lb;
    lb = lower_bound(_scaled_compatible_mass_means.begin(),
                     _scaled_compatible_mass_means.end(),
                     scaled_mass);
    if (lb < _scaled_compatible_mass_means.end())
    {
        int d = lb - _scaled_compatible_mass_means.begin();
        if (*lb == scaled_mass || lb == _scaled_compatible_mass_means.begin())
        {
            double var = _scaled_mle_variances[d];
            assert (!isnan(var) && !isinf(var));
            if (var < 0)
                var = 0;
            return var;
        }
        
        //in between two points on the scale.
        d--;
        
        if (d < 0)
        {
            fprintf(stderr, "ARG d < 0, d = %d \n", d);
        }
        
        if (d >= _scaled_compatible_mass_means.size())
        {
            fprintf(stderr, "ARG d >= _scaled_compatible_mass_means.size(), d = %d\n", d);
        }
        if (d >= _scaled_mle_variances.size())
        {
            fprintf(stderr, "ARG d >= _scaled_mass_variances.size(), d = %d\n", d);
        }
        
        double x1_mean = _scaled_compatible_mass_means[d];
        double x2_mean = _scaled_compatible_mass_means[d + 1];
        
        double y1_var = _scaled_mle_variances[d];
        double y2_var = _scaled_mle_variances[d + 1];
        double slope = 0.0;
        if (x2_mean != x1_mean)
        {
            slope = (y2_var - y1_var) / (x2_mean-x1_mean);
        }
        else if (y1_var == y2_var)
        {
            assert (false); // should have a unique'd table
        }
        double mean_interp = _scaled_mle_variances[d] + slope*(scaled_mass - _scaled_compatible_mass_means[d]);
        if (mean_interp < 0)
            mean_interp = 0;
        
        assert (!isnan(mean_interp) && !isinf(mean_interp));
        return mean_interp;
    }
    else
    {
        assert (!isnan(scaled_mass) && !isinf(scaled_mass));
        return 0; // revert to poisson assumption
    }
}

void transform_counts_to_common_scale(const vector<double>& scale_factors,
                                      vector<LocusCountList>& sample_compatible_count_table)
{
    // Transform raw counts to the common scale
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        LocusCountList& p = sample_compatible_count_table[i];
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            assert (scale_factors.size() > j);
            p.counts[j] *= (1.0 / scale_factors[j]);
        }
    }
}

void calc_geometric_scaling_factors(const vector<LocusCountList>& sample_compatible_count_table,
                                    vector<double>& scale_factors)
{
    
    vector<double> log_geom_means(sample_compatible_count_table.size(), 0.0);
    
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_compatible_count_table[i];
        
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            //assert (log_geom_means.size() > j);
            //if (floor(p.counts[j]) > 0)
            //{
                log_geom_means[i] += (1.0/p.counts.size()) * log(floor(p.counts[j]));
            //}
            
        }
        //log_geom_means[i] = pow(log_geom_means[i], 1.0/(double)p.counts.size());
    }
    
    for (size_t j = 0; j < scale_factors.size(); ++j)
    {
        vector<double> tmp_counts;
        for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
        {
            if (log_geom_means[i] && !isinf(log_geom_means[i]) && !isnan(log_geom_means[i]) && floor(sample_compatible_count_table[i].counts[j]))
            {
                double gm = (double)log(floor(sample_compatible_count_table[i].counts[j])) - log_geom_means[i];
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

void calc_classic_fpkm_scaling_factors(const vector<LocusCountList>& sample_compatible_count_table,
                                       vector<double>& scale_factors)
{
    vector<double> total_counts(sample_compatible_count_table.size(), 0.0);
    
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_compatible_count_table[i];
        
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            total_counts[j] += floor(p.counts[j]);
        }
    }
    
    double all_library_counts = accumulate(total_counts.begin(), total_counts.end(), 0.0);
    if (all_library_counts == 0.0)
    {
        for (size_t i = 0; i < scale_factors.size(); ++i)
        {
            scale_factors[i] = 1.0;
        }
        return;
    }
    for (size_t i = 0; i < scale_factors.size(); ++i)
    {
        scale_factors[i] = total_counts[i] / all_library_counts;
    }
    
    double avg_scaling_factor = accumulate(scale_factors.begin(), scale_factors.end(), 0.0);
    avg_scaling_factor /= scale_factors.size();
    
    for (size_t i = 0; i < scale_factors.size(); ++i)
    {
        scale_factors[i] = scale_factors[i] / avg_scaling_factor;
    }
}

void calc_estimated_absolute_scaling_factors(vector<boost::shared_ptr<ReadGroupProperties> > & all_read_groups,
                                             vector<double>& scale_factors)
{
    assert (scale_factors.size() == all_read_groups.size());
    for (size_t i = 0; i < scale_factors.size(); ++i)
    {
        double med_cov = all_read_groups[i]->mode_transcript_coverage();
        if (med_cov == 0)
        {
            scale_factors[i] = 0.0;
        }
        else if (med_cov == -1)
        {
            fprintf(stderr, "Error: estimated-absolute requires pre-calculated CXB files\n");
            exit(1);
        }
        
        scale_factors[i] =  med_cov * 1000000000 / (all_read_groups[i]->normalized_map_mass()) ;
    }
    
//
//    void adjust_fpkms_by_median_coverage(const vector<vector<double > >& median_coverage,
//                                         FPKMTrackingTable& tracking)
//    {
//        for (FPKMTrackingTable::iterator itr = tracking.begin(); itr != tracking.end(); ++itr)
//        {
//            if (itr != tracking.end())
//            {
//                FPKMTracking& track = itr->second;
//                vector<FPKMContext>& fpkms = track.fpkm_series;
//                
//                for (size_t i = 0; i < fpkms.size(); ++i)
//                {
//                    for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
//                         ++j)
//                    {
//                        double& FPKM = fpkms[i].tracking_info_per_rep[j].fpkm;
//                        FPKM /= median_coverage[i][j];
//                        FPKM *= fpkms[i].tracking_info_per_rep[j].rg_props->normalized_map_mass();
//                        FPKM /= 1000000000;
//                    }
//                }
//            }
//        }
//    }
//
//    
}



void calc_quartile_scaling_factors(const vector<LocusCountList>& sample_compatible_count_table,
                                   vector<double>& scale_factors)
{
    
    if (sample_compatible_count_table.empty())
        return;
    
    vector<double> upper_quartiles(sample_compatible_count_table.front().counts.size(), 0.0);
    vector<double> total_common_masses(sample_compatible_count_table.front().counts.size(), 0.0);
    
    for (size_t i = 0; i < sample_compatible_count_table.front().counts.size(); ++i)
    {
        vector<double> common_scaled_counts;
        double total_common = 0.0;

        for (size_t j = 0; j < sample_compatible_count_table.size(); ++j)
        {
        
            //boost::shared_ptr<ReadGroupProperties> rg = bundle_factories[fac_idx];
            //double scaled_mass = scale_factors[fac_idx] * rg->total_map_mass();
            
            total_common += sample_compatible_count_table[j].counts[i];
            common_scaled_counts.push_back(sample_compatible_count_table[j].counts[i]);
        }
    
        sort(common_scaled_counts.begin(), common_scaled_counts.end());
        if (common_scaled_counts.empty())
            continue;
        
        int upper_quart_index = common_scaled_counts.size() * 0.75;
        double upper_quart_count = common_scaled_counts[upper_quart_index];
        upper_quartiles[i] = upper_quart_count;
        total_common_masses[i] = total_common;
    }
    
    long double total_mass = accumulate(total_common_masses.begin(), total_common_masses.end(), 0.0);
    long double total_norm_mass = accumulate(upper_quartiles.begin(), upper_quartiles.end(), 0.0);
    
    for (size_t i = 0; i < sample_compatible_count_table.front().counts.size(); ++i)
    {
        if (total_mass > 0)
        {
            double scaling_factor = upper_quartiles[i];
            scaling_factor /= (total_norm_mass / upper_quartiles.size());
            scale_factors[i] = scaling_factor;
        }
        else
        {
            scale_factors[i] = 1.0;
        }
    }
}

void calc_tmm_scaling_factors(const vector<LocusCountList>& sample_compatible_count_table,
                              vector<double>& scale_factors)
{
    scale_factors = vector<double>(sample_compatible_count_table.size(), 1.0);
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
//        if (est_scv <= 0)
//            return 0.0;
        
        if (est_scvs.size() < 2 || true_scvs.size() < 2)
        {
            return est_scv; // revert to poisson.
        }
        if (est_scv > est_scvs.back())
        {
            //fprintf(stderr, "Warning: extrapolating to the right\n");
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
//            if (mean_interp < est_scv)
//                mean_interp = est_scv;
            assert (!isnan(mean_interp) && !isinf(mean_interp));
            return mean_interp;
        }
        else if (est_scv < est_scvs.front())
        {
            //fprintf(stderr, "Warning: extrapolating to the left\n");
            
            // If we're extrapolating to the left, our fit is too coarse, but
            // that probably means we don't need SCV bias correction at all.
            return est_scv;
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
//                if (var < est_scv) // revert to poisson if underdispersed
//                    var = est_scv;
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
                fprintf(stderr, "Warning: SCV table does not have unique keys\n!");
                assert (false); // should have a unique'd table
            }
            double mean_interp = true_scvs[d] + slope*(est_scv - est_scvs[d]);
//            if (mean_interp < est_scv) // revert to poisson if underdispersed
//                mean_interp = est_scv;
            
            assert (!isnan(mean_interp) && !isinf(mean_interp));
            return mean_interp;
        }
        else
        {
            assert (!isnan(est_scv) && !isinf(est_scv));
            return est_scv; // revert to poisson assumption
        }
        
        return est_scv;
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
    for(double a = 0.0002; a < 2.0; a += 0.0002)
    {
        alpha_range.push_back(a);
    }
    
    for (double a = 2; a < 100.0; a += 1)
    {
        alpha_range.push_back(a);
    }
    
    BOOST_FOREACH (double alpha, alpha_range)
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
            LocusCountList locus_count("", nreps, 1, vector<string>(), vector<string>());
            for (size_t rep_idx = 0; rep_idx < nreps; ++rep_idx)
            {
                double gamma_draw = gamma(rng);
                if (gamma_draw == 0)
                {
                    locus_count.counts[rep_idx] = 0;
                    draws.push_back(0);
                }
                else
                {
                    boost::random::poisson_distribution<long, double> poisson(gamma_draw);
                    locus_count.counts[rep_idx] = poisson(rng);
                    draws.push_back(locus_count.counts[rep_idx]);
                    //fprintf(stderr, "%lg\t", locus_count.counts[rep_idx]);
                }
            }
            
            double mean = accumulate(locus_count.counts.begin(), locus_count.counts.end(), 0.0);
            if (mean == 0)
                continue;
            mean /= locus_count.counts.size();
            double var = 0.0;
            BOOST_FOREACH(double c,  locus_count.counts)
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
        
        double mean = accumulate(draws.begin(), draws.end(), 0.0);
        mean /= draws.size();
        double var = 0.0;
        BOOST_FOREACH(int c,  draws)
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
    true_to_est_scv_table.finalize();
}

void calculate_count_means_and_vars(const vector<LocusCountList>& sample_compatible_count_table,
                                    vector<pair<double, double> >& means_and_vars)
{
    
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        const LocusCountList& p = sample_compatible_count_table[i];
        double mean = accumulate(p.counts.begin(), p.counts.end(), 0.0);
        if (mean > 0.0 && p.counts.size() > 0)
            mean /= p.counts.size();
        
        double var = 0.0;
        double num_non_zero = 0;
        BOOST_FOREACH (double d, p.counts)
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
        means_and_vars.push_back(make_pair(mean, var));
    }
}
                              
boost::shared_ptr<MassDispersionModel>
fit_dispersion_model_helper(const string& condition_name,
                            const vector<double>& scale_factors,
                            const vector<LocusCountList>& sample_compatible_count_table)
{
    vector<pair<double, double> > compatible_means_and_vars;
    
    SCVInterpolator true_to_est_scv_table;
    
    int num_samples = sample_compatible_count_table.front().counts.size();
    if (no_scv_correction == false)
    {
        build_scv_correction_fit(num_samples, 10000, 100000, true_to_est_scv_table);
    }
    
    setuplf();  
    
    calculate_count_means_and_vars(sample_compatible_count_table, compatible_means_and_vars);

    sort(compatible_means_and_vars.begin(), compatible_means_and_vars.end());
    
    vector<double> compatible_count_means;
    vector<double> raw_variances;
    
    for(size_t i = 0; i < compatible_means_and_vars.size(); ++i)
    {
        if (compatible_means_and_vars[i].first > 0 && compatible_means_and_vars[i].second > 0.0)
        {
            compatible_count_means.push_back(compatible_means_and_vars[i].first);
            raw_variances.push_back(compatible_means_and_vars[i].second);
        }
    }
    
    if (compatible_count_means.size() < min_loci_for_fitting)
    {
        boost::shared_ptr<MassDispersionModel> disperser;
        disperser = boost::shared_ptr<MassDispersionModel>(new PoissonDispersionModel(condition_name));
        
        return disperser;
    }
    
    
    vector<double> fitted_values;
    
    // WARNING: locfit doesn't like undescores - need camel case for 
    // variable names
    
    char namebuf[256];
    sprintf(namebuf, "countMeans");
    vari* cm = createvar(namebuf,STREGULAR,compatible_count_means.size(),VDOUBLE);
    for (size_t i = 0; i < compatible_count_means.size(); ++i)
    {
        cm->dpr[i] = log(compatible_count_means[i]);
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
    
    double xim = 0;
    BOOST_FOREACH(double s, scale_factors)
    {
        if (s)
            xim += 1.0 / s;
    }
    xim /= scale_factors.size();
    
    int n = 0;
    sprintf(namebuf, "fittedVars");
    vari* cp = findvar(namebuf, 1, &n);
    assert(cp != NULL);
    for (size_t i = 0; i < cp->n; ++i)
    {
//        if (cp->dpr[i] >= 0)
//        {
            double mean = exp(cm->dpr[i]);
            double fitted_scv = (cp->dpr[i] - mean) / (mean * mean);
            double corrected_scv = true_to_est_scv_table.interpolate_scv(fitted_scv);
            double corrected_variance = mean + (corrected_scv * (mean * mean));
            double uncorrected_variance = mean + (fitted_scv * (mean * mean));
            //fitted_values.push_back(mean + (cp->dpr[i] - xim * mean));
            if (no_scv_correction == false && corrected_variance > uncorrected_variance)
                fitted_values.push_back(corrected_variance);
            else if (uncorrected_variance > 0)
                fitted_values.push_back(uncorrected_variance);
            else
                fitted_values.push_back(compatible_count_means[i]);
        
            
//        }
//        else
//        {
//            fitted_values.push_back(compatible_count_means[i]);
//        }
    }
    
    boost::shared_ptr<MassDispersionModel> disperser;
    disperser = boost::shared_ptr<MassDispersionModel>(new MassDispersionModel(condition_name, compatible_count_means, raw_variances, fitted_values));
    if (dispersion_method == POISSON)
        disperser = boost::shared_ptr<MassDispersionModel>(new PoissonDispersionModel(condition_name));
    
//    for (map<string, pair<double, double> >::iterator itr = labeled_mv_table.begin();
//         itr != labeled_mv_table.end();
//         ++itr)
//    {
//        string label = itr->first;
//        disperser->set_compatible_mean_and_var(itr->first, itr->second);
//    }
    
    return disperser;
}

boost::shared_ptr<MassDispersionModel>
fit_dispersion_model(const string& condition_name,
                     const vector<double>& scale_factors,
                     const vector<LocusCountList>& sample_compatible_count_table)
{
//    
//#if ENABLE_THREADS
//	boost::mutex::scoped_lock lock(_locfit_lock);
//#endif
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        if (sample_compatible_count_table[i].counts.size() <= 1)
        {
            // only one replicate - no point in fitting variance
            return boost::shared_ptr<MassDispersionModel>(new PoissonDispersionModel(condition_name));
        }
    }
#if ENABLE_THREADS
    _locfit_lock.lock();
#endif
    
    ProgressBar p_bar("Modeling fragment count overdispersion.",0);
    
    int max_transcripts = 0;
    BOOST_FOREACH(const LocusCountList& L, sample_compatible_count_table)
    {
        if (L.num_transcripts > max_transcripts)
        {
            max_transcripts = L.num_transcripts;
        }
    }
    
    boost::shared_ptr<MassDispersionModel>  model = fit_dispersion_model_helper(condition_name, scale_factors, sample_compatible_count_table);

#if ENABLE_THREADS
    _locfit_lock.unlock();
#endif
    return model;
}

void build_norm_table(const vector<LocusCountList>& full_count_table,
                      boost::shared_ptr<const map<string, LibNormStandards> > normalizing_standards,
                      vector<LocusCountList>& norm_table)
{
    // If we're using housekeeping genes or spike-in controls, select the rows we'll be using from the full count table.
    if (normalizing_standards)
    {
        for (size_t i = 0; i < full_count_table.size(); ++i)
        {
            const vector<string>& gene_ids = full_count_table[i].gene_ids;
            const vector<string>& gene_short_names = full_count_table[i].gene_short_names;
            
            // If the row has an ID that's in the table, take it.
            map<string, LibNormStandards>::const_iterator g_id_itr = normalizing_standards->end();
            map<string, LibNormStandards>::const_iterator g_name_itr = normalizing_standards->end();
            
            for (size_t j = 0; j < gene_ids.size(); ++j)
            {
                g_id_itr = normalizing_standards->find(gene_ids[j]);
                if (g_id_itr != normalizing_standards->end())
                {
                    break;
                }
            }
            
            if (g_id_itr != normalizing_standards->end())
            {
                norm_table.push_back(full_count_table[i]);
                continue;
            }

            for (size_t j = 0; j < gene_short_names.size(); ++j)
            {
                g_name_itr = normalizing_standards->find(gene_short_names[j]);
                if (g_name_itr != normalizing_standards->end())
                {
                    break;
                }
            }
            
            if (g_name_itr != normalizing_standards->end())
            {
                norm_table.push_back(full_count_table[i]);
                continue;
            }

        }
    }
    else // otherwise, just take all rows.
    {
        norm_table = full_count_table;
    }
}

void normalize_counts(vector<boost::shared_ptr<ReadGroupProperties> > & all_read_groups)
{
    vector<LocusCountList> sample_compatible_count_table;
    vector<LocusCountList> sample_total_count_table;
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        boost::shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        const vector<LocusCount>& raw_compatible_counts = rg_props->raw_compatible_counts();
        const vector<LocusCount>& raw_total_counts = rg_props->raw_total_counts();
        
        for (size_t j = 0; j < raw_compatible_counts.size(); ++j)
        {
            if (sample_compatible_count_table.size() == j)
            {
                const string& locus_id = raw_compatible_counts[j].locus_desc;
                int num_transcripts = raw_compatible_counts[j].num_transcripts;
                
                const vector<string>& gene_ids = raw_compatible_counts[j].gene_ids;
                const vector<string>& gene_short_names = raw_compatible_counts[j].gene_short_names;
                
                sample_compatible_count_table.push_back(LocusCountList(locus_id,all_read_groups.size(), num_transcripts, gene_ids, gene_short_names));
                sample_total_count_table.push_back(LocusCountList(locus_id,all_read_groups.size(), num_transcripts, gene_ids, gene_short_names));
            }
            double scaled = raw_compatible_counts[j].count;
            //sample_compatible_count_table[j].counts[i] = scaled * unscaling_factor;
            sample_compatible_count_table[j].counts[i] = floor(scaled);
            sample_total_count_table[j].counts[i] = floor(raw_total_counts[j].count);
            
            assert(sample_compatible_count_table[j].counts[i] >= 0 && !isinf(sample_compatible_count_table[j].counts[i]));
        }
    }
    
    vector<double> scale_factors(all_read_groups.size(), 0.0);
    
    vector<LocusCountList> norm_table;
    
    if (use_compat_mass)
    {
        build_norm_table(sample_compatible_count_table, lib_norm_standards, norm_table);
    }
    else // use_total_mass
    {
        assert(use_total_mass);
        build_norm_table(sample_total_count_table, lib_norm_standards, norm_table);
    }
    
    if (lib_norm_method == GEOMETRIC)
    {
        calc_geometric_scaling_factors(norm_table, scale_factors);
    }
    else if (lib_norm_method == CLASSIC_FPKM)
    {
        calc_classic_fpkm_scaling_factors(norm_table, scale_factors);
    }
    else if (lib_norm_method == QUARTILE)
    {
        calc_quartile_scaling_factors(norm_table, scale_factors);
    }
    else if (lib_norm_method == TMM)
    {
        calc_tmm_scaling_factors(norm_table, scale_factors);
    }
    else if (lib_norm_method == ESTIMATED_ABSOLUTE)
    {
        calc_estimated_absolute_scaling_factors(all_read_groups, scale_factors);
    }
    else
    {
        assert (false);
    }
    
    
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        boost::shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        rg_props->internal_scale_factor(scale_factors[i]);
    }
    
    assert(sample_compatible_count_table.size() == sample_total_count_table.size());
    
    // Transform raw counts to the common scale
    for (size_t i = 0; i < sample_compatible_count_table.size(); ++i)
    {
        LocusCountList& p = sample_compatible_count_table[i];
        for (size_t j = 0; j < p.counts.size(); ++j)
        {
            assert (scale_factors.size() > j);
            p.counts[j] *= (1.0 / scale_factors[j]);
        }
        
        LocusCountList& t = sample_total_count_table[i];
        for (size_t j = 0; j < t.counts.size(); ++j)
        {
            assert (scale_factors.size() > j);
            t.counts[j] *= (1.0 / scale_factors[j]);
        }
    }
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        boost::shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        vector<LocusCount> scaled_compatible_counts;
        for (size_t j = 0; j < sample_compatible_count_table.size(); ++j)
        {
            string& locus_id = sample_compatible_count_table[j].locus_desc;
            double count = sample_compatible_count_table[j].counts[i];
            int num_transcripts = sample_compatible_count_table[j].num_transcripts;
            
            const vector<string>& gids = sample_compatible_count_table[j].gene_ids;
            const vector<string>& gnms = sample_compatible_count_table[j].gene_short_names;
            
            LocusCount locus_count(locus_id, count, num_transcripts, gids, gnms);
            scaled_compatible_counts.push_back(locus_count);
        }
        rg_props->common_scale_compatible_counts(scaled_compatible_counts);
    }
    
    for (size_t i = 0; i < all_read_groups.size(); ++i)
    {
        boost::shared_ptr<ReadGroupProperties> rg_props = all_read_groups[i];
        vector<LocusCount> scaled_total_counts;
        for (size_t j = 0; j < sample_total_count_table.size(); ++j)
        {
            string& locus_id = sample_total_count_table[j].locus_desc;
            double count = sample_total_count_table[j].counts[i];
            int num_transcripts = sample_total_count_table[j].num_transcripts;
            
            const vector<string>& gids = sample_total_count_table[j].gene_ids;
            const vector<string>& gnms = sample_total_count_table[j].gene_short_names;
            
            LocusCount locus_count(locus_id, count, num_transcripts, gids, gnms);
            scaled_total_counts.push_back(locus_count);
        }
        rg_props->common_scale_total_counts(scaled_total_counts);
    }
    
    double avg_total_common_scaled_count = 0.0;
    
    for (size_t fac_idx = 0; fac_idx < all_read_groups.size(); ++fac_idx)
    {
        double total_common = 0.0;
        if (use_compat_mass)
        {
            for (size_t j = 0; j < sample_compatible_count_table.size(); ++j)
            {
                total_common += sample_compatible_count_table[j].counts[fac_idx];
            }
        }
        else
        {
            for (size_t j = 0; j < sample_compatible_count_table.size(); ++j)
            {
                total_common += sample_total_count_table[j].counts[fac_idx];
            }

        }
        
        avg_total_common_scaled_count += (1.0/all_read_groups.size()) * total_common;
    }
    
    BOOST_FOREACH(boost::shared_ptr<ReadGroupProperties> rg, all_read_groups)
    {
        rg->normalized_map_mass(avg_total_common_scaled_count);
    }
}

