/*
 *  sorting_hat.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/30/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <stdlib.h>
#include <getopt.h>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


#include "common.h"
#include "tokenize.h"

#include "differential.h"

using namespace std;
using namespace boost;

bool compute_row_matrix = false;
bool log_transform_fpkm = false;
bool output_row_density = false;

int k_clusters = 0;
int max_iterations = 1000;

string row_matrix_out_filename = "";
string row_density_out_filename = "";

#if ENABLE_THREADS
const char *short_options = "o:p:d:P:k:I:l";
#else
const char *short_options = "o:d:P:k:I:l";
#endif

static struct option long_options[] = {
    {"output-dir",			    required_argument,		 0,			 'o'},
    {"row-matrix",			    required_argument,       0,			 'd'},
    {"row-densities",			required_argument,       0,			 'P'},
    {"log-fpkm",			    no_argument,             0,			 'l'},
    {"k-means",                 required_argument,       0,			 'k'},
    {"max-iterations",          required_argument,       0,			 'I'},
#if ENABLE_THREADS
    {"num-threads",				required_argument,       0,          'p'},
#endif
    {0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "sorting_hat v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffdiff [options] <input.fpkm_tracking> <output.shout>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");  
    fprintf(stderr, "-d/--row-matrix              compute row distance matrix                           [ default:     off ]\n");
    fprintf(stderr, "-l/--log-fpkm                JS on log(fpkm+1) instead of                          [ default:     off ]\n");
#if ENABLE_THREADS
	fprintf(stderr, "-p/--num-threads             number of threads used during assembly                [ default:      1 ]\n");
#endif
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;

            case 'd':
			{
				compute_row_matrix = true;
                row_matrix_out_filename = optarg;
				break;
			}
                
            case 'P':
			{
				output_row_density = true;
                row_density_out_filename = optarg;
				break;
			}
                
            case 'l':
			{
				log_transform_fpkm = true;
				break;
			}
                
            case 'k':
			{
				k_clusters = (int)parseInt(1, "-k/--k-means arg must be at least 1", print_usage);
				break;
			}
                
            case 'I':
			{
				max_iterations = (int)parseInt(1, "-I/--max-iterations arg must be at least 1", print_usage);
				break;
			}
                
            case 'o':
			{
				output_dir = optarg;
				break;
			}
  
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	allow_junk_filtering = false;
	
	return 0;
}


struct ExprRecord
{
    ExprRecord() : 
        total_FPKM(0.0),
        total_FPKM_conf_hi(0.0),
        total_FPKM_conf_lo(0.0),
        max_FPKM(0.0),
        total_log_FPKM(std::numeric_limits<double>::infinity()),
        total_log_FPKM_conf_hi(std::numeric_limits<double>::infinity()),
        total_log_FPKM_conf_lo(std::numeric_limits<double>::infinity()),
        max_log_FPKM(std::numeric_limits<double>::infinity()),
        cluster_id(-1) {} 
    
    string tracking_id;
    string class_code;
    string nearest_ref_id;
    string gene_id;
    string gene_short_name;
    string tss_id;
    string locus; 
    
    string length;
    string coverage;
    string status;
    
    vector<double> FPKMs;
    vector<double> FPKM_conf_los;
    vector<double> FPKM_conf_his;
    
    vector<double> log_FPKMs;
    vector<double> log_FPKM_conf_los;
    vector<double> log_FPKM_conf_his;
    
    double total_FPKM;
    double total_FPKM_conf_hi;
    double total_FPKM_conf_lo;
    double max_FPKM;
    
    double total_log_FPKM;
    double total_log_FPKM_conf_hi;
    double total_log_FPKM_conf_lo;
    double max_log_FPKM;
    
    ublas::vector<double> cond_density;
    vector<double> cond_specificities;
    
    int cluster_id;
};

struct ClusterStats
{
    ClusterStats(int dim) 
    {
        mean = ublas::zero_vector<double>(dim);
        variance = 0.0;
    }
    
    ublas::vector<double> mean;
    double variance;
    vector<size_t> members;
};

struct SortByVariance
{
    bool operator()(const ClusterStats& lhs, const ClusterStats& rhs)
    {
        if (lhs.variance != rhs.variance)
            return lhs.variance < rhs.variance;
        else if (lhs.members.size() != rhs.members.size())
            return lhs.members.size() < rhs.members.size();
        return false;
    }
};


void assign_to_nearest_cluster(vector<ExprRecord>& expr_records, vector<ClusterStats>& clusters)
{
    
    int num_clusters = clusters.size();
    int num_conditions = expr_records.front().cond_density.size();
    
    BOOST_FOREACH(ClusterStats& c, clusters)
    {
        c.members.clear();
    }
    
    // E step - assign each record to nearest cluster
    for (size_t r = 0; r < expr_records.size(); ++r)
    {
        ExprRecord& rec = expr_records[r];
        
        // Skip unexpressed records;
        if (log_transform_fpkm)
        {
            if (rec.total_log_FPKM == 0 ||
                rec.total_log_FPKM == std::numeric_limits<double>::infinity())
            {
                continue;
            }
        }
        else
        {
            if (rec.total_FPKM == 0)
            {
                continue;
            }
        }    

        double min_dist = std::numeric_limits<double>::max();
        size_t closest_mean = 0;
        vector<ublas::vector<double> >  kappas;
        
        kappas.push_back(rec.cond_density);
        //cerr << rec.cond_density << endl;
        kappas.push_back(ublas::zero_vector<double>());
        
        for (size_t m = 0; m < clusters.size(); ++m)
        {
            kappas[1] = clusters[m].mean;
            double js = jensen_shannon_distance(kappas);
            if (js < min_dist)
            {
                closest_mean = m;
                min_dist = js;
            }
        }
        
        clusters[closest_mean].members.push_back(r);
    }
    
    // M step - calculate new means
    
    vector<ublas::vector<double> > tmp_means(num_clusters, 
                                             ublas::zero_vector<double>(num_conditions));
    // Calculate mean
    
    for (size_t c = 0; c < clusters.size(); ++c)
    {
        ClusterStats& cluster = clusters[c];
        BOOST_FOREACH(int m, cluster.members)
        {
            tmp_means[c] += expr_records[m].cond_density;
        }
    }
    
    for (size_t c = 0; c < clusters.size(); ++c)
    {
        if (clusters[c].members.size() > 0)
        {
            clusters[c].mean = tmp_means[c] / clusters[c].members.size();
        }
    }
    
    // And now the variance
    for (size_t c = 0; c < clusters.size(); ++c)
    {
        ClusterStats& cluster = clusters[c];
        if (cluster.members.empty())
            continue;
        
        vector<ublas::vector<double> >  kappas;
        kappas.push_back(cluster.mean);
        kappas.push_back(ublas::zero_vector<double>());
        BOOST_FOREACH (int m, cluster.members)
        {
            kappas[1] = expr_records[m].cond_density;
            double js = jensen_shannon_distance(kappas);
            cluster.variance += (js * js);
        }
        
        cluster.variance /= cluster.members.size();
    }
}

void get_assignments(const vector<ClusterStats>& clusters, 
                     vector<int>& assignments)
{
    BOOST_FOREACH (int& i, assignments)
    {
        i = -1;
    }
    
    for (size_t c = 0; c < clusters.size(); ++c)
    {
        for (size_t r = 0; r < clusters[c].members.size(); ++r)
        {
            assignments[clusters[c].members[r]] = c;
        }
    }
}


void split_cluster(const vector<ExprRecord>& expr_records, 
                   ClusterStats& to_split, 
                   vector<ClusterStats>& new_clusters)
{
    
    fprintf (stderr, "Splitting cluster...\n");
    
    int num_conditions = expr_records.front().cond_density.size();
    
    vector<ublas::vector<double> >  kappas;
    kappas.push_back(to_split.mean);
    
    vector<double> dist_from_mean;
    
    for (size_t m = 0; m < to_split.members.size(); ++m)
    {
        kappas.push_back(expr_records[to_split.members[m]].cond_density);
        double js = jensen_shannon_distance(kappas);
        dist_from_mean.push_back(js);
    }
    
    // pick the point farthest from the mean and the point closest.
    
    double largest_dist = -1.0;
    size_t farthest_point = 0;
    
    double smallest_dist = std::numeric_limits<double>::max();
    size_t closest_point = 0;
    
    for (size_t i = 0; i < dist_from_mean.size(); ++i)
    {
        if (largest_dist > dist_from_mean[i])
        {
            largest_dist = dist_from_mean[i];
            farthest_point = i;
        }
        
        if (smallest_dist < dist_from_mean[i])
        {
            smallest_dist = dist_from_mean[i];
            closest_point = i;
        }
    }
    
    vector<ublas::vector<double> > a_kappas;
    a_kappas.push_back(expr_records[to_split.members[closest_point]].cond_density);
    a_kappas.push_back(ublas::zero_vector<double>());
    
    vector<ublas::vector<double> > b_kappas;
    b_kappas.push_back(expr_records[to_split.members[farthest_point]].cond_density);
    b_kappas.push_back(ublas::zero_vector<double>());
    
    ClusterStats a_cluster(num_conditions);
    ClusterStats b_cluster(num_conditions);
    
    a_cluster.mean = ublas::zero_vector<double>(num_conditions);
    b_cluster.mean = ublas::zero_vector<double>(num_conditions);
    
    for (size_t m = 0; m < to_split.members.size(); ++m)
    {
        a_kappas[1] = expr_records[to_split.members[m]].cond_density;
        double a_js = jensen_shannon_distance(a_kappas);
        
        b_kappas[1] = expr_records[to_split.members[m]].cond_density;
        double b_js = jensen_shannon_distance(b_kappas);
        
        if (a_js < b_js)
        {
            a_cluster.members.push_back(to_split.members[m]);
            a_cluster.mean += expr_records[to_split.members[m]].cond_density;
        }
        else 
        {
            b_cluster.members.push_back(to_split.members[m]);    
            b_cluster.mean += expr_records[to_split.members[m]].cond_density;
        }
    }  
    
    a_cluster.mean /= a_cluster.members.size();
    b_cluster.mean /= b_cluster.members.size();
    
    // Calculate the variances of the split clusters.
    a_kappas[0] = a_cluster.mean;
    b_kappas[0] = b_cluster.mean;
    
    for (size_t m = 0; m < a_cluster.members.size(); ++m)
    {
        a_kappas[1] = expr_records[a_cluster.members[m]].cond_density;
        double js = jensen_shannon_distance(a_kappas);
        a_cluster.variance += (js * js);
    }
    
    a_cluster.variance /= a_cluster.members.size();
    
    for (size_t m = 0; m < b_cluster.members.size(); ++m)
    {
        b_kappas[1] = expr_records[b_cluster.members[m]].cond_density;
        double js = jensen_shannon_distance(b_kappas);
        b_cluster.variance += (js * js);
    }
    
    b_cluster.variance /= b_cluster.members.size();
    
    new_clusters.push_back(a_cluster);
    new_clusters.push_back(b_cluster);
}

void kmeans(vector<ExprRecord>& expr_records, int num_clusters, int num_iterations)
{
    // generate k random means
    size_t num_conditions = 0;
    for (size_t i = 1; i < expr_records.size(); ++i)
    {
        if (expr_records[i - 1].FPKMs.size() != expr_records[i].FPKMs.size())
        {
            fprintf(stderr, "Error: not all records have the same number of conditions!\n");
            exit(1);
        }
        num_conditions = expr_records[i].FPKMs.size();
    }
    
    // This is a typedef for a random number generator.
    typedef boost::minstd_rand base_generator_type;
    base_generator_type generator(time(NULL));

    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    
    // each mean gets an accumulation vector and a counter.  For now,
    // let's just keep ClusterStats thin and do all this outside the class. 
    vector<ClusterStats> clusters(num_clusters, ClusterStats(num_conditions));
    
    for (size_t i = 0; i < clusters.size(); ++i)
    {
        ClusterStats& c = clusters[i];
        ublas::vector<double>& mean_i = c.mean;
        
        for (size_t j = 0; j < mean_i.size(); ++j)
        {
            mean_i(j) = uni();
        }
        double total = accumulate(mean_i.begin(), mean_i.end(), 0.0);
        mean_i /= total;
        //cerr << mean_i << endl;
    }
    
    vector<int> assignments(expr_records.size());
    vector<int> prev_assignments;
    
    size_t iter = 0;
    for (;iter < num_iterations; ++iter)
    {        
        if (iter % 5 == 0)
        {
            fprintf(stderr, "Iteration # %d\n", iter);
        }
        
        assign_to_nearest_cluster(expr_records, clusters);
        get_assignments(clusters, assignments);
        
        sort(clusters.begin(), clusters.end(), SortByVariance());
        
        vector<ClusterStats> empty_clusters;
        vector<ClusterStats> full_clusters;
        
        // Look for empty clusters, and if found split others (in order of 
        // largest variance)
        for(size_t p = 0; p < clusters.size(); ++p)
        {
            const ClusterStats& P = clusters[p];
            if (P.members.size() == 0) // found empty cluster
            {
                empty_clusters.push_back(P);
            }
            else
            {
                full_clusters.push_back(P);   
            }
        }
        
        if (!empty_clusters.empty())
        {
            vector<ClusterStats> new_clusters;
            for (size_t c = 0; c < empty_clusters.size(); ++c)
            {
                //size_t full = 0;
                // split this cluster.
                ClusterStats to_split = full_clusters.back();
                full_clusters.pop_back();
                split_cluster(expr_records, to_split, new_clusters);
                //get_assignments(clusters, assignments);
                //prev_assignments = assignments;
            }
            
            full_clusters.insert(full_clusters.end(), new_clusters.begin(), new_clusters.end());
            clusters = full_clusters;
            get_assignments(clusters, assignments);
            prev_assignments = assignments;
        }
        
        else 
        {
            if (assignments == prev_assignments)
            {
                break;
            }
            else 
            {
                prev_assignments = assignments;
            }
            
        }
    }
    
    for (size_t a = 0; a < assignments.size(); ++a)
    {
        expr_records[a].cluster_id = assignments[a];
    }
    
    if (iter == num_iterations)
    {
        fprintf(stderr, "Warning: clustering did not converge after %d iterations\n", iter);
    }
}

void driver(FILE* fpkm_file, FILE* spec_out, FILE* row_matrix_out, FILE* row_density_out)
{
    char buf[10 * 1024];
    
    vector<string> sample_names;
    
    int line_num = 1;
    
    while(fgets(buf, sizeof(buf), fpkm_file))
    {
        if (buf[0])
        {
            // Chomp the newline
            char* nl = strrchr(buf, '\n');
            if (nl) *nl = 0;
            
            vector<string> tokens;
            tokenize(buf, "\t", tokens);
            
            if (tokens.size() < 16)
            {
                fprintf(stderr, "Error:  FPKM tracking files must have at least 12 columns\n");
                exit(1);
            }
            
            if ((tokens.size() - 10) % 3 != 0)
            {
                fprintf(stderr, "Error:  FPKM tracking files must have FPKM, FPKM_lo, and FPKM_hi columns for each sample\n");
                exit(1);
            }
            
            static const size_t first_sample_idx = 10;
            
            for (size_t i = first_sample_idx; i < tokens.size(); i += 3)
            {
                string FPKM_label = tokens[i];
                string::size_type name_end = FPKM_label.rfind("_FPKM");
                if (name_end == string::npos)
                {
                    fprintf(stderr, "Error:  Malformed FPKM column header %s.  Should end in \"_FPKM\".\n", FPKM_label.c_str());
                    exit(1);
                }
                string name = FPKM_label.substr(0, name_end);
                sample_names.push_back(name);
            }
            
            break;
        }
        
        line_num++;
    }
    
    fprintf(spec_out, "tracking_id\tclass_code\tnearest_ref\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tstatus\ttotal_FPKM\ttotal_FPKM_lo\ttotal_FPKM_hi\tmax_FPKM\tcluster_id");

    for (size_t i = 0; i < sample_names.size(); ++i)
    {
        fprintf(spec_out, "\t%s", sample_names[i].c_str());
    }
    fprintf(spec_out, "\n");
    
    vector<ExprRecord> expr_records;
    
    while(fgets(buf, sizeof(buf), fpkm_file))
    {
        // Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
        
        vector<string> tokens;
        tokenize(buf, "\t", tokens);
        
        if (((tokens.size() - 10) / 3) != sample_names.size())
        {
            fprintf(stderr, "Error:  Line %d has %lu columns, should have %lu\n", line_num, tokens.size(), (sample_names.size() * 3) + 6);
            exit(1);
        }
        
        ExprRecord rec;
        
        rec.tracking_id = tokens[0];
        rec.class_code = tokens[1];
        rec.nearest_ref_id = tokens[2];
        rec.gene_id = tokens[3];
        rec.gene_short_name = tokens[4];
        rec.tss_id = tokens[5];
        rec.locus = tokens[6];
        
        rec.length = tokens[7];
        rec.coverage = tokens[8];
        rec.status = tokens[9];
        
        static const size_t first_sample_idx = 10;
        
        size_t num_samples = sample_names.size();
        ublas::vector<double> u1 = ublas::unit_vector<double>(sample_names.size(), 0);
        ublas::vector<double> u2 = ublas::unit_vector<double>(sample_names.size(), 1);
        
        vector<ublas::vector<double> > norm_kappas;
        norm_kappas.push_back(u1);
        norm_kappas.push_back(u2);
        
        double norm_js = jensen_shannon_distance(norm_kappas);
        double max_FPKM = -1;
        double max_log_FPKM = -1;
        
        for (size_t i = first_sample_idx; i < tokens.size() - 2; i += 3)
        {
            string FPKM_string = tokens[i];
            double fpkm = atof(FPKM_string.c_str());
            
            string FPKM_conf_lo_string = tokens[i+1];
            double fpkm_conf_lo = atof(FPKM_conf_lo_string.c_str());
            
            string FPKM_conf_hi_string = tokens[i+2];
            double fpkm_conf_hi = atof(FPKM_conf_hi_string.c_str());
            
            double log_fpkm = log10(fpkm + 1.0);
            double log_fpkm_conf_lo = log10(fpkm_conf_lo + 1.0);
            double log_fpkm_conf_hi = log10(fpkm_conf_hi + 1.0);
            
            if (isnan(fpkm) || isnan(fpkm_conf_lo) || isnan(fpkm_conf_hi))
            {
                fprintf (stderr, "Warning: gene %s (%s) on line %d has FPKM = NaN\n", 
                         rec.tracking_id.c_str(), rec.gene_short_name.c_str(), line_num); 
                fpkm = 0.0;
                fpkm_conf_lo = 0;
                fpkm_conf_hi = 0;
                
                log_fpkm = std::numeric_limits<double>::infinity();
                log_fpkm_conf_lo = std::numeric_limits<double>::infinity();
                log_fpkm_conf_hi = std::numeric_limits<double>::infinity();
            }
            else
            {
                if (log_fpkm > max_log_FPKM)
                {
                    max_log_FPKM = log_fpkm;
                }
                
                if (fpkm > max_FPKM)
                {
                    max_FPKM = fpkm;
                }
            }
            
            rec.FPKMs.push_back(fpkm);
            rec.FPKM_conf_los.push_back(fpkm_conf_lo);
            rec.FPKM_conf_his.push_back(fpkm_conf_hi);

            rec.log_FPKMs.push_back(log_fpkm);
            rec.log_FPKM_conf_los.push_back(log_fpkm_conf_lo);
            rec.log_FPKM_conf_his.push_back(log_fpkm_conf_hi);
        }
        
        rec.total_FPKM = accumulate(rec.FPKMs.begin(), rec.FPKMs.end(), 0.0);
        rec.total_FPKM_conf_lo = accumulate(rec.FPKM_conf_los.begin(), rec.FPKM_conf_los.end(), 0.0);
        rec.total_FPKM_conf_hi = accumulate(rec.FPKM_conf_his.begin(), rec.FPKM_conf_his.end(), 0.0);
        rec.max_FPKM = max_FPKM;
        
        rec.total_log_FPKM = accumulate(rec.log_FPKMs.begin(), rec.log_FPKMs.end(), 0.0);
        rec.total_log_FPKM_conf_lo = accumulate(rec.log_FPKM_conf_los.begin(), rec.log_FPKM_conf_los.end(), 0.0);
        rec.total_log_FPKM_conf_hi = accumulate(rec.log_FPKM_conf_his.begin(), rec.log_FPKM_conf_his.end(), 0.0);
        rec.max_log_FPKM = max_log_FPKM;
        
        assert (!isnan(rec.total_FPKM) && !isinf(rec.total_FPKM));
        
        rec.cond_density = ublas::vector<double>(sample_names.size());
        rec.cond_specificities = vector<double>(sample_names.size(), std::numeric_limits<double>::max());
        
                
        if (rec.total_FPKM == 0.0 || (log_transform_fpkm && (rec.total_log_FPKM == 0.0 || rec.total_log_FPKM == std::numeric_limits<double>::infinity())))
        {
             rec.cond_density = ublas::zero_vector<double>(sample_names.size());
        }
        else 
        {
            for (size_t i = 0; i < rec.cond_density.size(); ++i)
            {
                if (log_transform_fpkm)
                {

                    rec.cond_density(i) = rec.log_FPKMs[i] / rec.total_log_FPKM;
                }
                else 
                {
                    rec.cond_density(i) = rec.FPKMs[i] / rec.total_FPKM;
                }

                
            }
        
            //fpkm_dists.push_back(FPKM_dist);
            
            //cerr << tracking_id << FPKM_dist<< endl;
            
            const size_t N = sample_names.size();
            
            assert (N >= 2);
            
            vector<ublas::vector<double> > kappas;
            kappas.push_back(rec.cond_density);
            kappas.push_back(ublas::zero_vector<double>(sample_names.size()));
            
            for (size_t i = 0; i < sample_names.size(); ++i)
            {
                ublas::vector<double> specific_vec = ublas::unit_vector<double>(N, i);
                kappas[1] = specific_vec;
                
                double js = jensen_shannon_distance(kappas);
                js /= norm_js;
                
                rec.cond_specificities[i] = 1.0 - js;
            }
        }
        
        expr_records.push_back(rec);
        
        line_num++;
    }
    
    if (k_clusters > 0)
    {
        kmeans(expr_records, k_clusters, max_iterations);
    }
    
    
    for (size_t i = 0; i < expr_records.size(); ++i)
    {
        const ExprRecord& rec = expr_records[i];
        
        char cluster_str[256];
        
        if (k_clusters == 0)
        {
            sprintf(cluster_str, "-");
        }
        else 
        {
            sprintf(cluster_str, "%d", rec.cluster_id);
        }
        
        fprintf(spec_out, 
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%lg\t%s",
                rec.tracking_id.c_str(),
                rec.class_code.c_str(),
                rec.nearest_ref_id.c_str(),
                rec.gene_id.c_str(),
                rec.gene_short_name.c_str(),
                rec.tss_id.c_str(),
                rec.locus.c_str(),
                rec.length.c_str(),
                rec.coverage.c_str(),
                rec.status.c_str(),
                rec.total_FPKM,
                rec.total_FPKM_conf_lo,
                rec.total_FPKM_conf_hi,
                rec.max_FPKM,
                cluster_str);
        
        for (size_t i = 0; i < rec.cond_density.size(); ++i)
        {
            fprintf(spec_out, "\t%g", rec.cond_specificities[i]);
        }
        
        fprintf(spec_out, "\n");
    }
    
//    if (row_matrix_out)
//    {
//         for (size_t i = 0; i < fpkm_dists.size(); ++i)
//         {
//             const ublas::vector<double>& i_vec = fpkm_dists[i];
//             vector<ublas::vector<double> > kappas;
//             kappas.push_back(i_vec);
//             kappas.push_back(ublas::zero_vector<double>(i_vec.size()));
//             if (i % 100 == 0)
//             {
//                 fprintf(stderr, "%lu of %lu (%g percent)\n", i, fpkm_dists.size(), (double)i/(double)fpkm_dists.size()); 
//             }
//             
//             for(size_t j = 0; j < fpkm_dists.size(); ++j)
//             {
//                 const ublas::vector<double>& j_vec = fpkm_dists[j];
//                 kappas[1] = j_vec;
//                 double i_j_js = jensen_shannon_distance(kappas);
//                 if (j == fpkm_dists.size() - 1)
//                 {
//                     fprintf(row_matrix_out, "%g\n", i_j_js);
//                 }
//                 else 
//                 {
//                     fprintf(row_matrix_out, "%g\t", i_j_js);
//                 }
//
//             }
//         }
//    }
//    
    if (row_density_out)
    {
        for (size_t i = 0; i < expr_records.size(); ++i)
        {
            const ExprRecord& rec = expr_records[i];
            const ublas::vector<double>& i_vec = rec.cond_density;
            fprintf(row_density_out, "%s\t", rec.tracking_id.c_str());                    
            for (size_t j = 0; j < i_vec.size(); ++j)
            {
                if (j == i_vec.size() - 1)
                {
                    fprintf(row_density_out, "%g\n", i_vec(j));
                }
                else 
                {
                    fprintf(row_density_out, "%g\t", i_vec(j));
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string fpkm_filename = argv[optind++];

    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string spec_out_filename = argv[optind++];
    
    
    FILE* fpkm_file = fopen(fpkm_filename.c_str(), "r");
    if (!fpkm_file)
    {
        fprintf(stderr, "Error: cannot open FPKM tracking file %s for reading\n",
                fpkm_filename.c_str());
        exit(1);
    }
    
    FILE* spec_out_file = fopen(spec_out_filename.c_str(), "w");
    if (!spec_out_file)
    {
        fprintf(stderr, "Error: cannot open output file %s for writing\n",
                spec_out_filename.c_str());
        exit(1);
    }
    
    FILE* row_matrix_out = NULL;
    if (compute_row_matrix)
    {
        row_matrix_out = fopen(row_matrix_out_filename.c_str(), "w");
        if (!row_matrix_out)
        {
            fprintf(stderr, "Error: cannot open output file %s for writing\n",
                    row_matrix_out_filename.c_str());
            exit(1);
        }
    }
    
    FILE* row_density_out = NULL;
    if (output_row_density)
    {
        row_density_out = fopen(row_density_out_filename.c_str(), "w");
        if (!row_density_out)
        {
            fprintf(stderr, "Error: cannot open output file %s for writing\n",
                    row_density_out_filename.c_str());
            exit(1);
        }
    }
    
    driver(fpkm_file, spec_out_file, row_matrix_out, row_density_out);
	
	return 0;
}


