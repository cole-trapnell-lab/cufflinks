/*
 *  cuffdiff.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 10/21/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <numeric>
#include <cfloat>
#include <iostream>

#include "common.h"
#include "hits.h"
#include "bundles.h"
#include "abundances.h"
#include "tokenize.h"
#include "biascorrection.h"
#include "update_check.h"

#include <boost/thread.hpp>
#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>

#include "differential.h"

extern "C" {
#include "locfit/local.h"
}

// Need at least this many reads in a locus to do any testing on it

vector<string> sample_labels;

double FDR = 0.05; 
bool samples_are_time_series = false;
using namespace std;
using namespace boost;

// We leave out the short codes for options that don't take an argument
#if ENABLE_THREADS
const char *short_options = "m:p:s:c:I:j:L:M:o:b:TNqvuF:C:";
#else
const char *short_options = "m:s:c:I:j:L:M:o:b:TNqvuF:C:";
#endif

static struct option long_options[] = {
{"labels",					required_argument,		 0,			 'L'},
{"seed",                    required_argument,		 0,			 OPT_RANDOM_SEED},
{"norm-standards-file",     required_argument,		 0,			 OPT_NORM_STANDARDS_FILE},
{"use-sample-sheet",        no_argument,             0,			 OPT_USE_SAMPLE_SHEET},
{"output-dir",			    required_argument,		 0,			 'o'},
{"verbose",			    	no_argument,			 0,			 'v'},
{"quiet",			    	no_argument,			 0,			 'q'},
#if ENABLE_THREADS
{"num-threads",				required_argument,       0,          'p'},
#endif
{"library-type",		    required_argument,		 0,			 OPT_LIBRARY_TYPE},
{"no-update-check",         no_argument,             0,          OPT_NO_UPDATE_CHECK},
{"compatible-hits-norm",    no_argument,	 		 0,	         OPT_USE_COMPAT_MASS},
{"total-hits-norm",         no_argument,	 		 0,	         OPT_USE_TOTAL_MASS},
    
// Some options for testing different stats policies
{"max-bundle-frags",        required_argument,       0,          OPT_MAX_FRAGS_PER_BUNDLE}, 
{"library-norm-method",     required_argument,       0,          OPT_LIB_NORM_METHOD},
{"output-format",           required_argument,       0,          OPT_OUTPUT_FORMAT},
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "cuffnorm v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION);
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffnorm [options] <transcripts.gtf> <sample1_expr.cxb> <sample2_expr.cxb> [... sampleN_expr.cxb]\n");
	fprintf(stderr, "   Supply replicate CXB files as comma separated lists for each condition: sample1_rep1.cxb,sample1_rep2.cxb,...sample1_repM.cxb\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "  -o/--output-dir              write all output files to this directory              [ default:     ./ ]\n");
    fprintf(stderr, "  -L/--labels                  comma-separated list of condition labels\n");
    fprintf(stderr, "  --norm-standards-file        Housekeeping/spike genes to normalize libraries       [ default:   NULL ]\n"); // NOT YET DOCUMENTED, keep secret for now
#if ENABLE_THREADS
	fprintf(stderr, "  -p/--num-threads             number of threads used during quantification          [ default:      1 ]\n");
#endif
    fprintf(stderr, "  --library-type               Library prep used for input reads                     [ default:  below ]\n");
    fprintf(stderr, "  --library-norm-method        Method used to normalize library sizes                [ default:  below ]\n");
    fprintf(stderr, "  --output-format              Format for output tables                              [ default:  below ]\n");
    
    fprintf(stderr, "\nAdvanced Options:\n");
    fprintf(stderr, "  --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:   TRUE ]\n");
    fprintf(stderr, "  --total-hits-norm            count all hits for normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, "  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
	fprintf(stderr, "  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
    fprintf(stderr, "  --seed                       value of random number generator seed                 [ default:      0 ]\n");
    fprintf(stderr, "  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]\n");
    print_library_table();
    print_lib_norm_method_table();
    print_output_format_table();
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    string sample_label_list;
    string dispersion_method_str;
    string lib_norm_method_str;
    string output_format_str;

    do {
        next_option = getopt_long_only(argc, argv, short_options, long_options, &option_index);
        if (next_option == -1)     /* Done with options. */
            break;
        switch (next_option) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                break;
                
			case 'p':
				num_threads = (uint32_t)parseInt(1, "-p/--num-threads arg must be at least 1", print_usage);
				break;
            case 'F':
				min_isoform_fraction = parseFloat(0, 1.0, "-F/--min-isoform-fraction must be between 0 and 1.0", print_usage);
				break;
            case 'L':
				sample_label_list = optarg;
				break;


			case OPT_NORM_STANDARDS_FILE:
			{
				norm_standards_filename = optarg;
				break;
			}
            case OPT_USE_SAMPLE_SHEET:
			{
                use_sample_sheet = true;
				break;
			}
            case 'v':
			{
				if (cuff_quiet)
				{
					fprintf(stderr, "Warning: Can't be both verbose and quiet!  Setting verbose only.\n");
				}
				cuff_quiet = false;
				cuff_verbose = true;
				break;
			}
			case 'q':
			{
				if (cuff_verbose)
				{
					fprintf(stderr, "Warning: Can't be both verbose and quiet!  Setting quiet only.\n");
				}
				cuff_verbose = false;
				cuff_quiet = true;
				break;
			}
            case 'o':
			{
				output_dir = optarg;
				break;
			}


            case OPT_LIBRARY_TYPE:
			{
				library_type = optarg;
				break;
			}

            case OPT_NO_UPDATE_CHECK:
            {
                no_update_check = true;
                break;
            }
            case OPT_RANDOM_SEED:
            {
                random_seed = parseInt(0, "--seed must be at least 0", print_usage);
                break;
            }
            case OPT_USE_COMPAT_MASS:
            {
                use_compat_mass = true;
                break;
            }
            case OPT_USE_TOTAL_MASS:
            {
                use_total_mass = true;
                break;
            }
            case OPT_LIB_NORM_METHOD:
			{
				lib_norm_method_str = optarg;
				break;
			}
            case OPT_OUTPUT_FORMAT:
			{
				output_format_str = optarg;
				break;
			}
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	if (library_type != "")
    {
        map<string, ReadGroupProperties>::iterator lib_itr = 
		library_type_table.find(library_type);
        if (lib_itr == library_type_table.end())
        {
            fprintf(stderr, "Error: Library type %s not supported\n", library_type.c_str());
            exit(1);
        }
        else 
        {
            if (library_type == "transfrags")
            {
                allow_junk_filtering = false;
            }
            global_read_properties = &lib_itr->second;
        }
    }
    else
    {
        
    }
    
    // Set the count dispersion method to use
    if (dispersion_method_str == "")
    {
        dispersion_method_str = default_dispersion_method;
    }
    
    map<string, DispersionMethod>::iterator disp_itr = 
    dispersion_method_table.find(dispersion_method_str);
    if (disp_itr == dispersion_method_table.end())
    {
        fprintf(stderr, "Error: Dispersion method %s not supported\n", dispersion_method_str.c_str());
        exit(1);
    }
    else 
    {
        dispersion_method = disp_itr->second;
    }

    // Set the library size normalization method to use
    if (lib_norm_method_str == "")
    {
        lib_norm_method_str = default_lib_norm_method;
    }
    
    map<string, LibNormalizationMethod>::iterator lib_norm_itr =
    lib_norm_method_table.find(lib_norm_method_str);
    if (lib_norm_itr == lib_norm_method_table.end())
    {
        fprintf(stderr, "Error: Normalization method %s not supported\n", lib_norm_method_str.c_str());
        exit(1);
    }
    else
    {
        lib_norm_method = lib_norm_itr->second;
    }

    // Set the count dispersion method to use
    if (output_format_str == "")
    {
        output_format_str = default_output_format;
    }
    
    map<string, OutputFormat>::iterator output_itr =
    output_format_table.find(output_format_str);
    if (output_itr == output_format_table.end())
    {
        fprintf(stderr, "Error: Output format %s not supported\n", output_format_str.c_str());
        exit(1);
    }
    else
    {
        output_format = output_itr->second;
    }

    if (use_total_mass && use_compat_mass)
    {
        fprintf (stderr, "Error: please supply only one of --compatibile-hits-norm and --total-hits-norm\n");
        exit(1);
    }
    
    tokenize(sample_label_list, ",", sample_labels);
    
	allow_junk_filtering = false;
	
	return 0;
}

#if ENABLE_THREADS
boost::mutex inspect_lock;
#endif

void inspect_map_worker(ReplicatedBundleFactory& fac,
                        int& tmp_min_frag_len, 
                        int& tmp_max_frag_len,
                        IdToLocusMap& id_to_locus_map)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    int min_f = std::numeric_limits<int>::max();
    int max_f = 0;
    
    fac.inspect_replicate_maps(min_f, max_f, id_to_locus_map);
    
#if ENABLE_THREADS
    inspect_lock.lock();
#endif
    tmp_min_frag_len = min(min_f, tmp_min_frag_len);
    tmp_max_frag_len = max(max_f, tmp_max_frag_len);
#if ENABLE_THREADS
    inspect_lock.unlock();
#endif
}

boost::shared_ptr<TrackingDataWriter> tracking_data_writer;

bool quantitate_next_locus(const RefSequenceTable& rt,
                           vector<boost::shared_ptr<ReplicatedBundleFactory> >& bundle_factories,
                           boost::shared_ptr<TestLauncher> launcher)
{
    for (size_t i = 0; i < bundle_factories.size(); ++i)
    {
        boost::shared_ptr<SampleAbundances> s_ab = boost::shared_ptr<SampleAbundances>(new SampleAbundances);
        
#if ENABLE_THREADS					
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
            
        }
        
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        boost::shared_ptr<HitBundle> pBundle = boost::shared_ptr<HitBundle>(new HitBundle());
        bool non_empty = bundle_factories[i]->next_bundle(*pBundle, true);
        
        if (pBundle->compatible_mass() > 0)
        {
            thread quantitate(sample_worker,
                              non_empty,
                              pBundle,
                              boost::ref(rt),
                              boost::ref(*(bundle_factories[i])),
                              s_ab,
                              i,
                              launcher,
                              false);
        }
        else
        {
            sample_worker(non_empty,
                          pBundle,
                          boost::ref(rt),
                          boost::ref(*(bundle_factories[i])),
                          s_ab,
                          i,
                          launcher,
                          false);
            locus_thread_pool_lock.lock();
            locus_curr_threads--;
            locus_thread_pool_lock.unlock();
        }
#else
        HitBundle bundle;
        bool non_empty = sample_factory.next_bundle(bundle);
        
        sample_worker(non_emtpy,
                      pBundle,
                      boost::ref(rt),
                      boost::ref(*(bundle_factories[i])),
                      s_ab,
                      i,
                      launcher,
                      false);
#endif
    }
    return true;
}

void parse_sample_sheet_file(FILE* sample_sheet_file,
                             vector<string>& sample_labels,
                             vector<string>& sam_hit_filename_lists)
{
    
    char pBuf[10 * 1024];
    size_t non_blank_lines_read = 0;
    
    sample_labels.clear();
    
    map<string, vector<string> > sample_groups;
    
    while (fgets(pBuf, 10*1024, sample_sheet_file))
    {
        if (strlen(pBuf) > 0)
        {
            char* nl = strchr(pBuf, '\n');
            if (nl)
                *nl = 0;
            
            string pBufstr = pBuf;
            string trimmed = boost::trim_copy(pBufstr);
            
            if (trimmed.length() > 0 && trimmed[0] != '#')
            {
                non_blank_lines_read++;
                vector<string> columns;
                tokenize(trimmed, "\t", columns);
                
                if (non_blank_lines_read == 1)
                    continue;
                
                if (columns.size() < 2)
                {
                    if (columns.size() > 0)
                        fprintf(stderr, "Malformed record in sample sheet: \n   >  %s\n", pBuf);
                    else
                        continue;
                }
                
                string sam_file = columns[0];
                string sample_group = columns[1];
                
                pair<map<string, vector<string> >::iterator, bool> inserted = sample_groups.insert(make_pair(sample_group, vector<string>()));
                inserted.first->second.push_back(sam_file);
            }
        }
    }
    
    for (map<string, vector<string> >::iterator itr = sample_groups.begin();
         itr != sample_groups.end(); ++itr)
    {
        sample_labels.push_back(itr->first);
        string sam_list = boost::join(itr->second, ",");
        sam_hit_filename_lists.push_back(sam_list);
    }
}

void parse_norm_standards_file(FILE* norm_standards_file)
{
    char pBuf[10 * 1024];
    size_t non_blank_lines_read = 0;
    
    boost::shared_ptr<map<string, LibNormStandards> > norm_standards(new map<string, LibNormStandards>);
    
    while (fgets(pBuf, 10*1024, norm_standards_file))
    {
        if (strlen(pBuf) > 0)
        {
            char* nl = strchr(pBuf, '\n');
            if (nl)
                *nl = 0;
            
            string pBufstr = pBuf;
            string trimmed = boost::trim_copy(pBufstr);
            
            if (trimmed.length() > 0 && trimmed[0] != '#')
            {
                non_blank_lines_read++;
                vector<string> columns;
                tokenize(trimmed, "\t", columns);
                
                if (non_blank_lines_read == 1)
                    continue;
                
                if (columns.size() < 1) // 
                {
                    continue;
                }
                
                string gene_id = columns[0];
                LibNormStandards L;
                norm_standards->insert(make_pair(gene_id, L));
            }
        }
    }
    lib_norm_standards = norm_standards;
}


void print_variability_models(FILE* var_model_out, const vector<boost::shared_ptr<ReplicatedBundleFactory> >& factories)
{

    fprintf(var_model_out, "condition\tlocus\tcompatible_count_mean\tcompatible_count_var\ttotal_count_mean\ttotal_count_var\tfitted_var\n");
    
    for (size_t i = 0; i < factories.size(); ++i)
    {
        string factor_name = factories[i]->condition_name();
        boost::shared_ptr<ReadGroupProperties> rg = factories[i]->factories()[0]->read_group_properties();
        boost::shared_ptr<MassDispersionModel const> model = rg->mass_dispersion_model();
//        const vector<double>& means = model->scaled_compatible_mass_means();
//        const vector<double>& raw_vars  = model->scaled_compatible_variances();
        
        const vector<LocusCount>& common_scale_compatible_counts = rg->common_scale_compatible_counts();
        for (size_t j = 0; j < common_scale_compatible_counts.size(); ++j)
        {
            string locus_desc = common_scale_compatible_counts[j].locus_desc;
            pair<double, double> compat_mean_and_var = model->get_compatible_mean_and_var(locus_desc);
            pair<double, double> total_mean_and_var = model->get_total_mean_and_var(locus_desc);
//            double total_compat_count = 0;
//            if (itr != locus_to_total_count_table.end())
//                total_compat_count = itr->second.count;
            
            
            fprintf(var_model_out, "%s\t%s\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    factor_name.c_str(),
                    locus_desc.c_str(),
                    compat_mean_and_var.first,
                    compat_mean_and_var.second,
                    total_mean_and_var.first,
                    total_mean_and_var.second,
                    model->scale_mass_variance(compat_mean_and_var.first));
        }
    }
    fclose(var_model_out);

}

void driver(FILE* ref_gtf, FILE* mask_gtf, FILE* contrast_file, FILE* norm_standards_file, vector<string>& sam_hit_filename_lists, Outfiles& outfiles)
{

	ReadTable it;
	RefSequenceTable rt(true, false);
    
	vector<boost::shared_ptr<Scaffold> > ref_mRNAs;
	
	vector<boost::shared_ptr<ReplicatedBundleFactory> > bundle_factories;
    vector<boost::shared_ptr<ReadGroupProperties> > all_read_groups;
    set<boost::shared_ptr<ReadGroupProperties> > serialized_read_groups;
    
	for (size_t i = 0; i < sam_hit_filename_lists.size(); ++i)
	{
        vector<string> sam_hit_filenames;
        tokenize(sam_hit_filename_lists[i], ",", sam_hit_filenames);
        
        vector<boost::shared_ptr<BundleFactory> > replicate_factories;
        
        string condition_name = sample_labels[i];
        
        for (size_t j = 0; j < sam_hit_filenames.size(); ++j)
        {
            boost::shared_ptr<HitFactory> hs;
            boost::shared_ptr<BundleFactory> hf;
            try
            {
                hs = boost::shared_ptr<HitFactory>(new PrecomputedExpressionHitFactory(sam_hit_filenames[j], it, rt));
                hf = boost::shared_ptr<BundleFactory>(new PrecomputedExpressionBundleFactory(static_pointer_cast<PrecomputedExpressionHitFactory>(hs)));
            }
            
            catch(boost::archive::archive_exception & e)
            {
                try
                {
                    hs = boost::shared_ptr<HitFactory>(new BAMHitFactory(sam_hit_filenames[j], it, rt));
                }
                catch (std::runtime_error& e) 
                {
                    try
                    {
//                        fprintf(stderr, "File %s doesn't appear to be a valid BAM file, trying SAM...\n",
//                                sam_hit_filenames[j].c_str());
                        hs = boost::shared_ptr<HitFactory>(new SAMHitFactory(sam_hit_filenames[j], it, rt));
                    }
                    catch (std::runtime_error& e)
                    {
                        fprintf(stderr, "Error: cannot open file %s for reading. Unrecognized file type\n",
                                sam_hit_filenames[j].c_str());
                        exit(1);
                    }
                }
                hf = boost::shared_ptr<BundleFactory>(new BundleFactory(hs, REF_DRIVEN));
            }
            
            
            boost::shared_ptr<ReadGroupProperties> rg_props(new ReadGroupProperties);
            
            if (global_read_properties)
            {
                *rg_props = *global_read_properties;
            }
            else 
            {
                *rg_props = hs->read_group_properties();
            }
            
            rg_props->checked_parameters(hs->read_group_properties().checked_parameters());
            rg_props->condition_name(condition_name);
            rg_props->replicate_num(j);
            rg_props->file_path(sam_hit_filenames[j]);
            
            all_read_groups.push_back(rg_props);
            
            hf->read_group_properties(rg_props);
            
            replicate_factories.push_back(hf);
            //replicate_factories.back()->set_ref_rnas(ref_mRNAs);
        }
        
        bundle_factories.push_back(boost::shared_ptr<ReplicatedBundleFactory>(new ReplicatedBundleFactory(replicate_factories, condition_name)));
	}
    
    boost::crc_32_type ref_gtf_crc_result;
    ::load_ref_rnas(ref_gtf, rt, ref_mRNAs, ref_gtf_crc_result, corr_bias, false);
    if (ref_mRNAs.empty())
        return;
    
    vector<boost::shared_ptr<Scaffold> > mask_rnas;
    if (mask_gtf)
    {
        boost::crc_32_type mask_gtf_crc_result;
        ::load_ref_rnas(mask_gtf, rt, mask_rnas, mask_gtf_crc_result, false, false);
    }
    
    BOOST_FOREACH (boost::shared_ptr<ReplicatedBundleFactory> fac, bundle_factories)
    {
        // NOTE: we aren't deep copying here to save memory.  We can pull this off
        // because unless bias correction or multiread correction is enabled,
        // we don't need to write into Scaffold objects.
        fac->set_ref_rnas(ref_mRNAs, false);
        if (mask_gtf) 
            fac->set_mask_rnas(mask_rnas);
    }
    
    if (norm_standards_file != NULL)
    {
        parse_norm_standards_file(norm_standards_file);
    }

    validate_cross_sample_parameters(all_read_groups);
    
    vector<pair<size_t, size_t > > contrasts;
    
#if ENABLE_THREADS
    locus_num_threads = num_threads;
#endif

	int tmp_min_frag_len = numeric_limits<int>::max();
	int tmp_max_frag_len = 0;
	
    boost::shared_ptr<IdToLocusMap> id_to_locus_map(new IdToLocusMap(boost::shared_ptr<map<string, set<string> > >(new map<string, set<string> >())));
    
	ProgressBar p_bar("Inspecting maps and determining fragment length distributions.",0);
	BOOST_FOREACH (boost::shared_ptr<ReplicatedBundleFactory> fac, bundle_factories)
    {
#if ENABLE_THREADS	
        while(1)
        {
            locus_thread_pool_lock.lock();
            if (locus_curr_threads < locus_num_threads)
            {
                break;
            }
            
            locus_thread_pool_lock.unlock();
            
            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
        }
        
        locus_curr_threads++;
        locus_thread_pool_lock.unlock();
        
        thread inspect(inspect_map_worker,
                       boost::ref(*fac),
                       boost::ref(tmp_min_frag_len),
                       boost::ref(tmp_max_frag_len),
                       boost::ref(*id_to_locus_map));
#else
        inspect_map_worker(boost::ref(*fac),
                           boost::ref(tmp_min_frag_len),
                           boost::ref(tmp_max_frag_len),
                           boost::ref(*id_to_locus_map));
#endif
    }
    
    // wait for the workers to finish up before reporting everthing.
#if ENABLE_THREADS	
    while(1)
    {
        locus_thread_pool_lock.lock();
        if (locus_curr_threads == 0)
        {
            locus_thread_pool_lock.unlock();
            break;
        }
        locus_thread_pool_lock.unlock();
        
        boost::this_thread::sleep(boost::posix_time::milliseconds(5));
    }
#endif
    
    normalize_counts(all_read_groups);
    
    long double total_norm_mass = 0.0;
    long double total_mass = 0.0;
    BOOST_FOREACH (boost::shared_ptr<ReadGroupProperties> rg_props, all_read_groups)
    {
        total_norm_mass += rg_props->normalized_map_mass();
        total_mass += rg_props->total_map_mass();
        rg_props->clear_count_tables();
    }

	min_frag_len = tmp_min_frag_len;
    max_frag_len = tmp_max_frag_len;
	
	final_est_run = false;
	
	double num_bundles = (double)bundle_factories[0]->num_bundles();
	
    p_bar = ProgressBar("Calculating preliminary abundance estimates", num_bundles);
    
    Tracking tracking;
    
	final_est_run = true;
	p_bar = ProgressBar("Normalizing expression levels for locus", num_bundles);
                                                     
    tracking_data_writer = boost::shared_ptr<TrackingDataWriter>(new TrackingDataWriter(bundle_factories.size(), &outfiles, &tracking, all_read_groups, sample_labels, &p_bar, id_to_locus_map));
    
	while (true)
	{
        //boost::shared_ptr<vector<boost::shared_ptr<SampleAbundances> > > abundances(new vector<boost::shared_ptr<SampleAbundances> >());
        quantitate_next_locus(rt, bundle_factories, tracking_data_writer);
        bool more_loci_remain = false;
        BOOST_FOREACH (boost::shared_ptr<ReplicatedBundleFactory> rep_fac, bundle_factories) 
        {
            if (rep_fac->bundles_remain())
            {
                more_loci_remain = true;
                break;
            }
        }
        if (!more_loci_remain)
        {
            // wait for the workers to finish up before doing the cross-sample testing.
#if ENABLE_THREADS	
            while(1)
            {
                locus_thread_pool_lock.lock();
                if (locus_curr_threads == 0)
                {
                    locus_thread_pool_lock.unlock();
                    break;
                }
                
                locus_thread_pool_lock.unlock();
                
                boost::this_thread::sleep(boost::posix_time::milliseconds(5));
                
            }
#endif
            break;
        }
    }
	
	p_bar.complete();
    
}

void open_outfiles_for_writing_cuffdiff_format(Outfiles& outfiles)
{
    
    static const int filename_buf_size = 2048;
    
	char isoform_fpkm_tracking_name[filename_buf_size];
	sprintf(isoform_fpkm_tracking_name, "%s/isoforms.fpkm_tracking", output_dir.c_str());
	FILE* isoform_fpkm_out = fopen(isoform_fpkm_tracking_name, "w");
	if (!isoform_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level FPKM tracking file %s for writing\n",
				isoform_fpkm_tracking_name);
		exit(1);
	}
	outfiles.isoform_fpkm_tracking_out = isoform_fpkm_out;
    
	char tss_group_fpkm_tracking_name[filename_buf_size];
	sprintf(tss_group_fpkm_tracking_name, "%s/tss_groups.fpkm_tracking", output_dir.c_str());
	FILE* tss_group_fpkm_out = fopen(tss_group_fpkm_tracking_name, "w");
	if (!tss_group_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level FPKM tracking file %s for writing\n",
				tss_group_fpkm_tracking_name);
		exit(1);
	}
	outfiles.tss_group_fpkm_tracking_out = tss_group_fpkm_out;
    
	char cds_fpkm_tracking_name[filename_buf_size];
	sprintf(cds_fpkm_tracking_name, "%s/cds.fpkm_tracking", output_dir.c_str());
	FILE* cds_fpkm_out = fopen(cds_fpkm_tracking_name, "w");
	if (!cds_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open CDS level FPKM tracking file %s for writing\n",
				cds_fpkm_tracking_name);
		exit(1);
	}
	outfiles.cds_fpkm_tracking_out = cds_fpkm_out;
	
	char gene_fpkm_tracking_name[filename_buf_size];
	sprintf(gene_fpkm_tracking_name, "%s/genes.fpkm_tracking", output_dir.c_str());
	FILE* gene_fpkm_out = fopen(gene_fpkm_tracking_name, "w");
	if (!gene_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open gene-level FPKM tracking file %s for writing\n",
				gene_fpkm_tracking_name);
		exit(1);
	}
	outfiles.gene_fpkm_tracking_out = gene_fpkm_out;
    
    char isoform_count_tracking_name[filename_buf_size];
	sprintf(isoform_count_tracking_name, "%s/isoforms.count_tracking", output_dir.c_str());
	FILE* isoform_count_out = fopen(isoform_count_tracking_name, "w");
	if (!isoform_count_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level count tracking file %s for writing\n",
				isoform_count_tracking_name);
		exit(1);
	}
	outfiles.isoform_count_tracking_out = isoform_count_out;
    
	char tss_group_count_tracking_name[filename_buf_size];
	sprintf(tss_group_count_tracking_name, "%s/tss_groups.count_tracking", output_dir.c_str());
	FILE* tss_group_count_out = fopen(tss_group_count_tracking_name, "w");
	if (!tss_group_count_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level count tracking file %s for writing\n",
				tss_group_count_tracking_name);
		exit(1);
	}
	outfiles.tss_group_count_tracking_out = tss_group_count_out;
    
	char cds_count_tracking_name[filename_buf_size];
	sprintf(cds_count_tracking_name, "%s/cds.count_tracking", output_dir.c_str());
	FILE* cds_count_out = fopen(cds_count_tracking_name, "w");
	if (!cds_count_out)
	{
		fprintf(stderr, "Error: cannot open CDS level count tracking file %s for writing\n",
				cds_count_tracking_name);
		exit(1);
	}
	outfiles.cds_count_tracking_out = cds_count_out;
	
	char gene_count_tracking_name[filename_buf_size];
	sprintf(gene_count_tracking_name, "%s/genes.count_tracking", output_dir.c_str());
	FILE* gene_count_out = fopen(gene_count_tracking_name, "w");
	if (!gene_count_out)
	{
		fprintf(stderr, "Error: cannot open gene-level count tracking file %s for writing\n",
				gene_count_tracking_name);
		exit(1);
	}
	outfiles.gene_count_tracking_out = gene_count_out;
    
    char isoform_rep_tracking_name[filename_buf_size];
	sprintf(isoform_rep_tracking_name, "%s/isoforms.read_group_tracking", output_dir.c_str());
	FILE* isoform_rep_out = fopen(isoform_rep_tracking_name, "w");
	if (!isoform_rep_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level read group tracking file %s for writing\n",
				isoform_rep_tracking_name);
		exit(1);
	}
	outfiles.isoform_rep_tracking_out = isoform_rep_out;
    
	char tss_group_rep_tracking_name[filename_buf_size];
	sprintf(tss_group_rep_tracking_name, "%s/tss_groups.read_group_tracking", output_dir.c_str());
	FILE* tss_group_rep_out = fopen(tss_group_rep_tracking_name, "w");
	if (!tss_group_rep_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level read group tracking file %s for writing\n",
				tss_group_rep_tracking_name);
		exit(1);
	}
	outfiles.tss_group_rep_tracking_out = tss_group_rep_out;
    
	char cds_rep_tracking_name[filename_buf_size];
	sprintf(cds_rep_tracking_name, "%s/cds.read_group_tracking", output_dir.c_str());
	FILE* cds_rep_out = fopen(cds_rep_tracking_name, "w");
	if (!cds_rep_out)
	{
		fprintf(stderr, "Error: cannot open CDS level read group tracking file %s for writing\n",
				cds_rep_tracking_name);
		exit(1);
	}
	outfiles.cds_rep_tracking_out = cds_rep_out;
	
	char gene_rep_tracking_name[filename_buf_size];
	sprintf(gene_rep_tracking_name, "%s/genes.read_group_tracking", output_dir.c_str());
	FILE* gene_rep_out = fopen(gene_rep_tracking_name, "w");
	if (!gene_rep_out)
	{
		fprintf(stderr, "Error: cannot open gene-level read group tracking file %s for writing\n",
				gene_rep_tracking_name);
		exit(1);
	}
	outfiles.gene_rep_tracking_out = gene_rep_out;
    
    char read_group_info_name[filename_buf_size];
	sprintf(read_group_info_name, "%s/read_groups.info", output_dir.c_str());
	FILE* read_group_out = fopen(read_group_info_name, "w");
	if (!read_group_out)
	{
		fprintf(stderr, "Error: cannot open read group info file %s for writing\n",
				read_group_info_name);
		exit(1);
	}
	outfiles.read_group_info_out = read_group_out;
    
    char run_info_name[filename_buf_size];
	sprintf(run_info_name, "%s/run.info", output_dir.c_str());
	FILE* run_info_out = fopen(run_info_name, "w");
	if (!run_info_out)
	{
		fprintf(stderr, "Error: cannot open run info file %s for writing\n",
				run_info_name);
		exit(1);
	}
	outfiles.run_info_out = run_info_out;

}

void open_outfiles_for_writing_simple_table_format(Outfiles& outfiles)
{
    static const int filename_buf_size = 2048;
    
	char isoform_fpkm_tracking_name[filename_buf_size];
	sprintf(isoform_fpkm_tracking_name, "%s/isoforms.fpkm_table", output_dir.c_str());
	FILE* isoform_fpkm_out = fopen(isoform_fpkm_tracking_name, "w");
	if (!isoform_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level FPKM table %s for writing\n",
				isoform_fpkm_tracking_name);
		exit(1);
	}
	outfiles.isoform_fpkm_tracking_out = isoform_fpkm_out;
    
	char tss_group_fpkm_tracking_name[filename_buf_size];
	sprintf(tss_group_fpkm_tracking_name, "%s/tss_groups.fpkm_table", output_dir.c_str());
	FILE* tss_group_fpkm_out = fopen(tss_group_fpkm_tracking_name, "w");
	if (!tss_group_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level FPKM table %s for writing\n",
				tss_group_fpkm_tracking_name);
		exit(1);
	}
	outfiles.tss_group_fpkm_tracking_out = tss_group_fpkm_out;
    
	char cds_fpkm_tracking_name[filename_buf_size];
	sprintf(cds_fpkm_tracking_name, "%s/cds.fpkm_table", output_dir.c_str());
	FILE* cds_fpkm_out = fopen(cds_fpkm_tracking_name, "w");
	if (!cds_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open CDS level FPKM table %s for writing\n",
				cds_fpkm_tracking_name);
		exit(1);
	}
	outfiles.cds_fpkm_tracking_out = cds_fpkm_out;
	
	char gene_fpkm_tracking_name[filename_buf_size];
	sprintf(gene_fpkm_tracking_name, "%s/genes.fpkm_table", output_dir.c_str());
	FILE* gene_fpkm_out = fopen(gene_fpkm_tracking_name, "w");
	if (!gene_fpkm_out)
	{
		fprintf(stderr, "Error: cannot open gene-level FPKM table %s for writing\n",
				gene_fpkm_tracking_name);
		exit(1);
	}
	outfiles.gene_fpkm_tracking_out = gene_fpkm_out;
    
    char isoform_count_tracking_name[filename_buf_size];
	sprintf(isoform_count_tracking_name, "%s/isoforms.count_table", output_dir.c_str());
	FILE* isoform_count_out = fopen(isoform_count_tracking_name, "w");
	if (!isoform_count_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level count table %s for writing\n",
				isoform_count_tracking_name);
		exit(1);
	}
	outfiles.isoform_count_tracking_out = isoform_count_out;
    
	char tss_group_count_tracking_name[filename_buf_size];
	sprintf(tss_group_count_tracking_name, "%s/tss_groups.count_table", output_dir.c_str());
	FILE* tss_group_count_out = fopen(tss_group_count_tracking_name, "w");
	if (!tss_group_count_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level count table %s for writing\n",
				tss_group_count_tracking_name);
		exit(1);
	}
	outfiles.tss_group_count_tracking_out = tss_group_count_out;
    
	char cds_count_tracking_name[filename_buf_size];
	sprintf(cds_count_tracking_name, "%s/cds.count_table", output_dir.c_str());
	FILE* cds_count_out = fopen(cds_count_tracking_name, "w");
	if (!cds_count_out)
	{
		fprintf(stderr, "Error: cannot open CDS level count table %s for writing\n",
				cds_count_tracking_name);
		exit(1);
	}
	outfiles.cds_count_tracking_out = cds_count_out;
	
	char gene_count_tracking_name[filename_buf_size];
	sprintf(gene_count_tracking_name, "%s/genes.count_table", output_dir.c_str());
	FILE* gene_count_out = fopen(gene_count_tracking_name, "w");
	if (!gene_count_out)
	{
		fprintf(stderr, "Error: cannot open gene-level count table %s for writing\n",
				gene_count_tracking_name);
		exit(1);
	}
	outfiles.gene_count_tracking_out = gene_count_out;
 
    char isoform_attr_name[filename_buf_size];
	sprintf(isoform_attr_name, "%s/isoforms.attr_table", output_dir.c_str());
	FILE* isoform_attr_out = fopen(isoform_attr_name, "w");
	if (!isoform_attr_out)
	{
		fprintf(stderr, "Error: cannot open isoform-level attribute table %s for writing\n",
				isoform_attr_name);
		exit(1);
	}
	outfiles.isoform_attr_out = isoform_attr_out;
    
	char tss_group_attr_name[filename_buf_size];
	sprintf(tss_group_attr_name, "%s/tss_groups.attr_table", output_dir.c_str());
	FILE* tss_group_attr_out = fopen(tss_group_attr_name, "w");
	if (!tss_group_attr_out)
	{
		fprintf(stderr, "Error: cannot open TSS group-level attribute table %s for writing\n",
				tss_group_attr_name);
		exit(1);
	}
	outfiles.tss_group_attr_out = tss_group_attr_out;
    
	char cds_attr_name[filename_buf_size];
	sprintf(cds_attr_name, "%s/cds.attr_table", output_dir.c_str());
	FILE* cds_attr_out = fopen(cds_attr_name, "w");
	if (!cds_attr_out)
	{
		fprintf(stderr, "Error: cannot open CDS level attribute table %s for writing\n",
				cds_attr_name);
		exit(1);
	}
	outfiles.cds_attr_out = cds_attr_out;
	
	char gene_attr_name[filename_buf_size];
	sprintf(gene_attr_name, "%s/genes.attr_table", output_dir.c_str());
	FILE* gene_attr_out = fopen(gene_attr_name, "w");
	if (!gene_attr_out)
	{
		fprintf(stderr, "Error: cannot open gene-level attribute table %s for writing\n",
				gene_attr_name);
		exit(1);
	}
	outfiles.gene_attr_out = gene_attr_out;
    
    char read_group_info_name[filename_buf_size];
	sprintf(read_group_info_name, "%s/samples.table", output_dir.c_str());
	FILE* read_group_out = fopen(read_group_info_name, "w");
	if (!read_group_out)
	{
		fprintf(stderr, "Error: cannot open read group info file %s for writing\n",
				read_group_info_name);
		exit(1);
	}
	outfiles.read_group_info_out = read_group_out;
    
    char run_info_name[filename_buf_size];
	sprintf(run_info_name, "%s/run.info", output_dir.c_str());
	FILE* run_info_out = fopen(run_info_name, "w");
	if (!run_info_out)
	{
		fprintf(stderr, "Error: cannot open run info file %s for writing\n",
				run_info_name);
		exit(1);
	}
	outfiles.run_info_out = run_info_out;
    
}

void open_outfiles_for_writing(Outfiles& outfiles)
{
    if (output_format == CUFFDIFF_OUTPUT_FMT)
    {
        open_outfiles_for_writing_cuffdiff_format(outfiles);
    }
    else if (output_format == SIMPLE_TABLE_OUTPUT_FMT)
    {
        open_outfiles_for_writing_simple_table_format(outfiles);
    }
    else{
        fprintf(stderr, "Error: unrecognized output format!\n");
        exit(1);
    }

}

int main(int argc, char** argv)
{
    for (int i = 0; i < argc; ++i)
    {
        cmd_str += string(argv[i]) + " ";
    }
    
    init_library_table();
    init_dispersion_method_table();
    init_lib_norm_method_table();
    init_output_format_table();
    
    min_isoform_fraction = 1e-5;
    
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
    if (!use_total_mass && !use_compat_mass)
    {
        use_total_mass = false;
        use_compat_mass = true;   
    }
    
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    if (!no_update_check)
        check_version(PACKAGE_VERSION);
    
    string ref_gtf_filename = argv[optind++];
    vector<string> sam_hit_filenames;
    
    if (use_sample_sheet)
    {
        if  (optind < argc)
        {
            
            string sample_sheet_filename = argv[optind++];
            FILE* sample_sheet_file = NULL;
            if (sample_sheet_filename != "")
            {
                sample_sheet_file = fopen(sample_sheet_filename.c_str(), "r");
                if (!sample_sheet_file)
                {
                    fprintf(stderr, "Error: cannot open sample sheet file %s for reading\n",
                            sample_sheet_filename.c_str());
                    exit(1);
                }
            }
            parse_sample_sheet_file(sample_sheet_file, sample_labels, sam_hit_filenames);
        }
        else
        {
            fprintf(stderr, "Error: option --use-sample-sheet requires a single sample sheet filename instead of a list of SAM/BAM files\n");
        }
    }
    else
    {
        while(optind < argc)
        {
            string sam_hits_file_name = argv[optind++];
            sam_hit_filenames.push_back(sam_hits_file_name);
        }
        
        if (sample_labels.size() == 0)
        {
            for (size_t i = 1; i < sam_hit_filenames.size() + 1; ++i)
            {
                char buf[256];
                sprintf(buf, "q%lu", i);
                sample_labels.push_back(buf);
            }
        }
    }
    	
	while (sam_hit_filenames.size() < 2)
    {
        fprintf(stderr, "Error: cuffdiff requires at least 2 SAM files\n");
        exit(1);
    }
	
    
    if (sam_hit_filenames.size() != sample_labels.size())
    {
        fprintf(stderr, "Error: number of labels must match number of conditions\n");
        exit(1);
    }
    
    if (random_seed == -1)
        random_seed = time(NULL);
    
	// seed the random number generator - we'll need it for the importance
	// sampling during MAP estimation of the gammas
	srand48(random_seed);
	
	FILE* ref_gtf = NULL;
	if (ref_gtf_filename != "")
	{
		ref_gtf = fopen(ref_gtf_filename.c_str(), "r");
		if (!ref_gtf)
		{
			fprintf(stderr, "Error: cannot open reference GTF file %s for reading\n",
					ref_gtf_filename.c_str());
			exit(1);
		}
	}
	
	FILE* mask_gtf = NULL;
	if (mask_gtf_filename != "")
	{
		mask_gtf = fopen(mask_gtf_filename.c_str(), "r");
		if (!mask_gtf)
		{
			fprintf(stderr, "Error: cannot open mask GTF file %s for reading\n",
					mask_gtf_filename.c_str());
			exit(1);
		}
	}

	
    FILE* contrast_file = NULL;
	if (contrast_filename != "")
	{
		contrast_file = fopen(contrast_filename.c_str(), "r");
		if (!contrast_file)
		{
			fprintf(stderr, "Error: cannot open contrast file %s for reading\n",
					contrast_filename.c_str());
			exit(1);
		}
	}

    FILE* norm_standards_file = NULL;
	if (norm_standards_filename != "")
	{
		norm_standards_file = fopen(norm_standards_filename.c_str(), "r");
		if (!norm_standards_file)
		{
			fprintf(stderr, "Error: cannot open contrast file %s for reading\n",
					norm_standards_filename.c_str());
			exit(1);
		}
	}
    

	// Note: we don't want the assembly filters interfering with calculations 
	// here
	
	pre_mrna_fraction = 0.0;
    olap_radius = 0;
	
	Outfiles outfiles;
	
    if (output_dir != "")
    {
        int retcode = mkpath(output_dir.c_str(), 0777);
        if (retcode == -1)
        {
            if (errno != EEXIST)
            {
                fprintf (stderr, 
                         "Error: cannot create directory %s\n", 
                         output_dir.c_str());
                exit(1);
            }
        }
    }

    open_outfiles_for_writing(outfiles);
    
    driver(ref_gtf, mask_gtf, contrast_file, norm_standards_file, sam_hit_filenames, outfiles);
	
#if 0
    if (emit_count_tables)
    {
        dump_locus_variance_info(output_dir + string("/locus_var.txt"));
    }
#endif
    
	return 0;
}

