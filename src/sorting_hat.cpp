/*
 *  sorting_hat.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/30/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif

#include <stdlib.h>
#include <getopt.h>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "common.h"
#include "tokenize.h"

#include "differential.h"

using namespace std;
using namespace boost;

bool compute_row_matrix = false;
bool log_transform_fpkm = false;
bool output_row_density = false;

string row_matrix_out_filename = "";
string row_density_out_filename = "";

#if ENABLE_THREADS
const char *short_options = "o:p:d:P:";
#else
const char *short_options = "o:d";
#endif

static struct option long_options[] = {
    {"output-dir",			    required_argument,		 0,			 'o'},
    {"row-matrix",			    required_argument,       0,			 'd'},
    {"row-densities",			required_argument,       0,			 'P'},
    {"log-fpkm",			    no_argument,             0,			 'l'},
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
				compute_row_matrix = true;
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
            
            if (tokens.size() < 12)
            {
                fprintf(stderr, "Error:  FPKM tracking files must have at least 12 columns\n");
                exit(1);
            }
            
            if ((tokens.size() - 6) % 3 != 0)
            {
                fprintf(stderr, "Error:  FPKM tracking files must have FPKM, FPKM_lo, and FPKM_hi columns for each sample\n");
                exit(1);
            }
            
            static const size_t first_sample_idx = 6;
            
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
    
    if (log_transform_fpkm)
    {
        fprintf(spec_out, "tracking_id\tclass_code\tnearest_ref\tgene_short_name\ttss_id\tlocus\t1og10_total_FPKM\t1og10_total_FPKM_lo\t1og10_total_FPKM_hi");   
    }
    else 
    {
        fprintf(spec_out, "tracking_id\tclass_code\tnearest_ref\tgene_short_name\ttss_id\tlocus\ttotal_FPKM\ttotal_FPKM_lo\ttotal_FPKM_hi");
    }

    for (size_t i = 0; i < sample_names.size(); ++i)
    {
        fprintf(spec_out, "\t%s", sample_names[i].c_str());
    }
    fprintf(spec_out, "\n");
    
    vector<ublas::vector<double> > fpkm_dists;
    
    while(fgets(buf, sizeof(buf), fpkm_file))
    {
        // Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
        
        vector<string> tokens;
        tokenize(buf, "\t", tokens);
        
        if (((tokens.size() - 6) / 3) != sample_names.size())
        {
            fprintf(stderr, "Error:  Line %d has %lu columns, should have %lu\n", line_num, tokens.size(), (sample_names.size() * 3) + 6);
            exit(1);
        }
        
        string tracking_id = tokens[0];
        string class_code = tokens[1];
        string nearest_ref_id = tokens[2];
        string gene_short_name = tokens[3];
        string tss_id = tokens[4];
        string locus = tokens[5];
        static const size_t first_sample_idx = 6;
        
        vector<double> FPKMs;
        vector<double> FPKM_conf_los;
        vector<double> FPKM_conf_his;
        for (size_t i = first_sample_idx; i < tokens.size() - 2; i += 3)
        {
            string FPKM_string = tokens[i];
            double fpkm = atof(FPKM_string.c_str());
            
            string FPKM_conf_lo_string = tokens[i+1];
            double fpkm_conf_lo = atof(FPKM_conf_lo_string.c_str());
            
            string FPKM_conf_hi_string = tokens[i+2];
            double fpkm_conf_hi = atof(FPKM_conf_hi_string.c_str());
            
            if (isnan(fpkm))
            {
                fprintf (stderr, "Warning: gene %s (%s) on line %d has FPKM = NaN\n", 
                         tracking_id.c_str(), gene_short_name.c_str(), line_num); 
                fpkm = 0.0;
                fpkm_conf_lo = 0;
                fpkm_conf_hi = 0;
            }
            
            if (!log_transform_fpkm)
            {
                FPKMs.push_back(fpkm);
                FPKM_conf_los.push_back(fpkm_conf_lo);
                FPKM_conf_his.push_back(fpkm_conf_hi);
            }
            else 
            {
                FPKMs.push_back(log10(fpkm + 1.0));
                FPKM_conf_los.push_back(log10(fpkm_conf_lo + 1.0));
                FPKM_conf_his.push_back(log10(fpkm_conf_hi + 1.0));
            }

        }
        
        double total_FPKM = accumulate(FPKMs.begin(), FPKMs.end(), 0.0);
        double total_FPKM_lo = accumulate(FPKM_conf_los.begin(), FPKM_conf_los.end(), 0.0);
        double total_FPKM_hi = accumulate(FPKM_conf_his.begin(), FPKM_conf_his.end(), 0.0);
        
        assert (!isnan(total_FPKM) && !isinf(total_FPKM));
        
        ublas::vector<double> FPKM_dist(sample_names.size());
        vector<double> specificity_js (sample_names.size(), std::numeric_limits<double>::max());
        
        if (total_FPKM == 0.0)
        {
            FPKM_dist = ublas::zero_vector<double>(sample_names.size());
        }
        else 
        {
            for (size_t i = 0; i < FPKM_dist.size(); ++i)
            {
                if (log_transform_fpkm)
                {
                    FPKM_dist(i) = log10(FPKMs[i]) / total_FPKM;
                }
                else
                {
                    FPKM_dist(i) = FPKMs[i] / total_FPKM;
                }
            }
        
            fpkm_dists.push_back(FPKM_dist);
            
            //cerr << tracking_id << FPKM_dist<< endl;
            
            const size_t N = sample_names.size();
            
            assert (N >= 2);
            
            vector<ublas::vector<double> > kappas;
            kappas.push_back(FPKM_dist);
            kappas.push_back(ublas::zero_vector<double>(sample_names.size()));
            
            for (size_t i = 0; i < sample_names.size(); ++i)
            {
                ublas::vector<double> specific_vec = ublas::unit_vector<double>(N, i);
                kappas[1] = specific_vec;
                
                double js = jensen_shannon_div(kappas);
                specificity_js[i] = js;
            }
        }

        fprintf(spec_out, 
                "%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g",
                tracking_id.c_str(),
                class_code.c_str(),
                nearest_ref_id.c_str(),
                gene_short_name.c_str(),
                tss_id.c_str(),
                locus.c_str(), 
                total_FPKM,
                total_FPKM_lo,
                total_FPKM_hi);
        
        for (size_t i = 0; i < specificity_js.size(); ++i)
        {
            fprintf(spec_out, "\t%g", specificity_js[i]);
        }
        
        fprintf(spec_out, "\n");
        
        line_num++;
    }
    
    if (row_matrix_out)
    {
         for (size_t i = 0; i < fpkm_dists.size(); ++i)
         {
             const ublas::vector<double>& i_vec = fpkm_dists[i];
             vector<ublas::vector<double> > kappas;
             kappas.push_back(i_vec);
             kappas.push_back(ublas::zero_vector<double>(i_vec.size()));
             if (i % 100 == 0)
             {
                 fprintf(stderr, "%lu of %lu (%g percent)\n", i, fpkm_dists.size(), (double)i/(double)fpkm_dists.size()); 
             }
             
             for(size_t j = 0; j < fpkm_dists.size(); ++j)
             {
                 const ublas::vector<double>& j_vec = fpkm_dists[j];
                 kappas[1] = j_vec;
                 double i_j_js = jensen_shannon_div(kappas);
                 if (j == fpkm_dists.size() - 1)
                 {
                     fprintf(row_matrix_out, "%g\n", i_j_js);
                 }
                 else 
                 {
                     fprintf(row_matrix_out, "%g\t", i_j_js);
                 }

             }
         }
    }
    
    if (row_density_out)
    {
        for (size_t i = 0; i < fpkm_dists.size(); ++i)
        {
            const ublas::vector<double>& i_vec = fpkm_dists[i];
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


