#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H
/*
 *  differential.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/15/10.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <string>

#include <boost/thread.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include "abundances.h"
#include "jensen_shannon.h"
#include "replicates.h"
#include "tracking.h"

using namespace std;

enum TestStatus {
	NOTEST,  // successful calculation, test not performed
    LOWDATA, // unsuccessful calculation due to low data, test not performed
    HIDATA,  // skipped calculation due to too many reads data, test not performed
	OK,      // successful numerical calc, test performed
	FAIL     // numerical exception, test not performed
}; 

struct SampleDifferenceMetaData
{
    string locus_desc;
    set<string> gene_ids;
	set<string> gene_names;
	set<string> protein_ids;
	string description; // isoforms or tss groups (e.g.) involved in this test
};

// Stores the differential expression of an isoform or set of isoforms in two
// different samples, along with a significance test statistic for the difference.
class SampleDifference
{
public:
	SampleDifference();
	
	size_t sample_1;
	size_t sample_2;
	
	double value_1;
	double value_2;
	double differential;
	double test_stat;
	double p_value;
	double corrected_p;
	
    boost::shared_ptr<SampleDifferenceMetaData> meta_data;

	TestStatus test_status;
	bool significant;
    
//    double sample_1_param_1;
//    double sample_1_param_2;
//    
//    
//    double sample_2_param_1;
//    double sample_2_param_2;
//    
//    double null_param_1;
//    double null_param_2;
};

typedef map<string, SampleDifference > SampleDiffs;
typedef map<string, boost::shared_ptr<SampleDifferenceMetaData> > SampleDiffMetaDataTable;

struct Outfiles
{
	FILE* isoform_de_outfile;
	FILE* group_de_outfile;
	FILE* gene_de_outfile;
	FILE* cds_de_outfile;
	
	FILE* diff_splicing_outfile;
	FILE* diff_promoter_outfile;
	FILE* diff_cds_outfile;
	
	FILE* isoform_fpkm_tracking_out;
	FILE* tss_group_fpkm_tracking_out;
	FILE* gene_fpkm_tracking_out;
	FILE* cds_fpkm_tracking_out;
    
    FILE* isoform_count_tracking_out;
	FILE* tss_group_count_tracking_out;
	FILE* gene_count_tracking_out;
	FILE* cds_count_tracking_out;
    
    FILE* isoform_rep_tracking_out;
	FILE* tss_group_rep_tracking_out;
	FILE* gene_rep_tracking_out;
	FILE* cds_rep_tracking_out;

    FILE* isoform_attr_out;
	FILE* tss_group_attr_out;
	FILE* gene_attr_out;
	FILE* cds_attr_out;
    
    FILE* run_info_out;
    FILE* read_group_info_out;
    FILE* bias_out;
    FILE* var_model_out;
};

struct Tests
{
	vector<vector<SampleDiffs> > isoform_de_tests;
	vector<vector<SampleDiffs> > tss_group_de_tests;
	vector<vector<SampleDiffs> > gene_de_tests;
	vector<vector<SampleDiffs> > cds_de_tests;
	
	vector<vector<SampleDiffs> > diff_splicing_tests; // to be performed on the isoforms of a single tss group
	vector<vector<SampleDiffs> > diff_promoter_tests; // to be performed on the tss groups of a single gene
	vector<vector<SampleDiffs> > diff_cds_tests; // to be performed on the cds groups of a single gene
};

#if ENABLE_THREADS
    extern boost::mutex _launcher_lock;
#endif

#if ENABLE_THREADS
extern boost::mutex test_storage_lock; 
#endif

struct TestLauncher
{
private:
    TestLauncher(TestLauncher& rhs) {}
    
public:
    TestLauncher(int num_samples,
                 const vector<pair<size_t, size_t> >& contrasts,
                 Tests* tests,
                 Tracking* tracking,
                 ProgressBar* p_bar,
                 bool track_only = false)
    :
    _orig_workers(num_samples),
    _contrasts(contrasts),
    _tests(tests),
    _tracking(tracking),
    _p_bar(p_bar),
    _track_only(track_only)
    {
    }
    
    virtual ~TestLauncher() {} 
    
    void operator()();
    
    void register_locus(const string& locus_id);
    void abundance_avail(const string& locus_id, 
                         boost::shared_ptr<SampleAbundances> ab, 
                         size_t factory_id);
    virtual vector<vector<boost::shared_ptr<SampleAbundances> > > test_finished_loci();
    virtual void perform_testing(vector<boost::shared_ptr<SampleAbundances> > abundances);
    void record_tracking_data(vector<boost::shared_ptr<SampleAbundances> >& abundances);
    bool all_samples_reported_in(vector<boost::shared_ptr<SampleAbundances> >& abundances);
    bool all_samples_reported_in(const string& locus_id);
    
    void clear_tracking_data() { _tracking->clear(); }
    
    typedef list<pair<string, vector<boost::shared_ptr<SampleAbundances> > > > launcher_sample_table;
    
private:
    
    launcher_sample_table::iterator find_locus(const string& locus_id);
    
    int _orig_workers;
    vector<pair<size_t, size_t> > _contrasts;
    launcher_sample_table _samples;
    Tests* _tests;
    Tracking* _tracking;
    ProgressBar* _p_bar;
    bool _track_only;
    
};

// This routine doesn't actually do any tests.  Use this when you just want to
// Flush tracking data once all samples have reported in.
class TrackingDataWriter : public TestLauncher
{
    
public:
    TrackingDataWriter(int num_samples,
                       Outfiles* outfiles,
                       Tracking* tracking,
                       const vector<boost::shared_ptr<ReadGroupProperties> >& all_read_groups,
                       vector<string> sls,
                       ProgressBar* p_bar,
                       boost::shared_ptr<IdToLocusMap> id_to_locus_map)
    : TestLauncher(num_samples, vector<pair<size_t, size_t> >(), NULL, tracking, p_bar),
    _tracking(tracking),
    _outfiles(outfiles),
    _all_read_groups(all_read_groups),
    sample_labels(sls),
    headers_written(false),
    _id_to_locus_map(id_to_locus_map)
    {
        
    }
    
    void perform_testing(vector<boost::shared_ptr<SampleAbundances> > abundances);
    
private:
    Tracking* _tracking;
    Outfiles* _outfiles;
    const vector<boost::shared_ptr<ReadGroupProperties> >& _all_read_groups;
    bool headers_written; // this flag records whether we've written out the file headers yet.
    boost::shared_ptr<IdToLocusMap> _id_to_locus_map;
    vector<string> sample_labels;
    
    void print_FPKM_tracking_header(FILE* fout,
                                    const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage");
        FPKMTrackingTable::const_iterator first_itr = tracking.begin();
        if (first_itr != tracking.end())
        {
            const FPKMTracking& track = first_itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                fprintf(fout, "\t%s_FPKM\t%s_conf_lo\t%s_conf_hi\t%s_status", sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str());
            }
        }
        fprintf(fout, "\n");
        fflush(fout);
    }
    
    void print_count_tracking_header(FILE* fout,
                                     const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id");
        FPKMTrackingTable::const_iterator first_itr = tracking.begin();
        if (first_itr != tracking.end())
        {
            const FPKMTracking& track = first_itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                fprintf(fout, "\t%s_count\t%s_count_variance\t%s_count_uncertainty_var\t%s_count_dispersion_var\t%s_status", sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str(), sample_labels[i].c_str());
            }
        }
        fprintf(fout, "\n");
        fflush(fout);
    }
    
    void print_FPKM_simple_table_header(FILE* fout,
                                        const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id");
        FPKMTrackingTable::const_iterator first_itr = tracking.begin();
        if (first_itr != tracking.end())
        {
            const FPKMTracking& track = first_itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
                     ++j)
                {
                    const  string& condition_name = fpkms[i].tracking_info_per_rep[j].rg_props->condition_name();
                    int rep_num = fpkms[i].tracking_info_per_rep[j].rg_props->replicate_num();
                    fprintf(fout, "\t%s_%d", condition_name.c_str(), rep_num);
                }
            }
            fflush(fout);
        }
    }
    
    void print_count_simple_table_header(FILE* fout,
                                         const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id");
        FPKMTrackingTable::const_iterator first_itr = tracking.begin();
        if (first_itr != tracking.end())
        {
            const FPKMTracking& track = first_itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
                     ++j)
                {
                    const  string& condition_name = fpkms[i].tracking_info_per_rep[j].rg_props->condition_name();
                    int rep_num = fpkms[i].tracking_info_per_rep[j].rg_props->replicate_num();
                    fprintf(fout, "\t%s_%d", condition_name.c_str(), rep_num);
                }
            }
            fflush(fout);
        }
    }
    
    void print_FPKM_tracking(const string& target_desc,
                             FILE* fout,
                             const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            const FPKMTracking& track = itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            
            AbundanceStatus status = NUMERIC_OK;
            BOOST_FOREACH (const FPKMContext& c, fpkms)
            {
                if (c.status == NUMERIC_FAIL)
                    status = NUMERIC_FAIL;
            }
            
            string all_gene_ids = cat_strings(track.gene_ids);
            if (all_gene_ids == "")
                all_gene_ids = "-";
            
            string all_gene_names = cat_strings(track.gene_names);
            if (all_gene_names == "")
                all_gene_names = "-";
            
            string all_tss_ids = cat_strings(track.tss_ids);
            if (all_tss_ids == "")
                all_tss_ids = "-";
            
            char length_buff[33] = "-";
            if (track.length)
                sprintf(length_buff, "%d", track.length);
            
            fprintf(fout, "%s\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    description.c_str(),
                    track.classcode ? track.classcode : '-',
                    track.ref_match.c_str(),
                    all_gene_ids.c_str(),
                    all_gene_names.c_str(),
                    all_tss_ids.c_str(),
                    track.locus_tag.c_str(),
                    length_buff,
                    "-");
       		
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                double fpkm = fpkms[i].FPKM;
                //double std_dev = sqrt(fpkms[i].FPKM_variance);
                double fpkm_conf_hi = fpkms[i].FPKM_conf_hi;
                double fpkm_conf_lo = fpkms[i].FPKM_conf_lo;
                const char* status_str = "OK";
                
                if (fpkms[i].status == NUMERIC_OK)
                {
                    status_str = "OK";
                }
                else if (fpkms[i].status == NUMERIC_FAIL)
                {
                    status_str = "FAIL";
                }
                else if (fpkms[i].status == NUMERIC_LOW_DATA)
                {
                    status_str = "LOWDATA";
                }
                else if (fpkms[i].status == NUMERIC_HI_DATA)
                {
                    status_str = "HIDATA";
                }
                else
                {
                    assert(false);
                }
                
                fprintf(fout, "\t%lg\t%lg\t%lg\t%s", fpkm, fpkm_conf_lo, fpkm_conf_hi, status_str);
            }
            
            fprintf(fout, "\n");
            fflush(fout);
        }
    }
    
    void print_count_tracking(const string& target_desc,
                              FILE* fout,
                              const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            const FPKMTracking& track = itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            
            AbundanceStatus status = NUMERIC_OK;
            BOOST_FOREACH (const FPKMContext& c, fpkms)
            {
                if (c.status == NUMERIC_FAIL)
                    status = NUMERIC_FAIL;
            }
            
            fprintf(fout, "%s",
                    description.c_str());
            
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                const char* status_str = "OK";
                
                if (fpkms[i].status == NUMERIC_OK)
                {
                    status_str = "OK";
                }
                else if (fpkms[i].status == NUMERIC_FAIL)
                {
                    status_str = "FAIL";
                }
                else if (fpkms[i].status == NUMERIC_LOW_DATA)
                {
                    status_str = "LOWDATA";
                }
                else if (fpkms[i].status == NUMERIC_HI_DATA)
                {
                    status_str = "HIDATA";
                }
                else
                {
                    assert(false);
                }
                
                double external_counts = fpkms[i].count_mean;
                double external_count_var = fpkms[i].count_var;
                double uncertainty_var = fpkms[i].count_uncertainty_var;
                double dispersion_var = fpkms[i].count_dispersion_var;
                fprintf(fout, "\t%lg\t%lg\t%lg\t%lg\t%s", external_counts, external_count_var, uncertainty_var, dispersion_var, status_str);
            }
            
            fprintf(fout, "\n");
            fflush(fout);
        }
    }
    
    void print_FPKM_simple_table(const string& target_desc,
                                 FILE* fout,
                                 const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            fprintf(fout, "\n%s", description.c_str());

            const FPKMTracking& track = itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
                     ++j)
                {
                    double FPKM = fpkms[i].tracking_info_per_rep[j].fpkm;
                    fprintf(fout, "\t%lg", FPKM);
                }
            }
            fflush(fout);
        }
    }
    
    void print_count_simple_table(const string& target_desc,
                                  FILE* fout,
                                  const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            fprintf(fout, "\n%s", description.c_str());

            const FPKMTracking& track = itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
                     ++j)
                {
                    double count = fpkms[i].tracking_info_per_rep[j].count;
                    fprintf(fout, "\t%lg", count);
                }
            }
            fflush(fout);
        }
    }
    
    void remove_abundance_from_tracking_table(const string& target_desc,
                                              FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            tracking.erase(itr);
        }
    }
    
    void print_read_group_tracking_header(FILE* fout,
                                          const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id\tcondition\treplicate\traw_frags\tinternal_scaled_frags\texternal_scaled_frags\tFPKM\teffective_length\tstatus");
        fprintf(fout, "\n");
    }
    
    void print_read_group_tracking(const string& target_desc,
                                   FILE* fout,
                                   const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            const FPKMTracking& track = itr->second;
            const vector<FPKMContext>& fpkms = track.fpkm_series;
            
            for (size_t i = 0; i < fpkms.size(); ++i)
            {
                for (size_t j = 0; j != fpkms[i].tracking_info_per_rep.size();
                     ++j)
                {
                    double FPKM = fpkms[i].tracking_info_per_rep[j].fpkm;
                    double internal_count = fpkms[i].tracking_info_per_rep[j].count;
                    double external_count = internal_count / fpkms[i].tracking_info_per_rep[j].rg_props->external_scale_factor();
                    double raw_count = internal_count * fpkms[i].tracking_info_per_rep[j].rg_props->internal_scale_factor();
                    const  string& condition_name = fpkms[i].tracking_info_per_rep[j].rg_props->condition_name();
                    AbundanceStatus status = fpkms[i].tracking_info_per_rep[j].status;
                    
                    int rep_num = fpkms[i].tracking_info_per_rep[j].rg_props->replicate_num();
                    
                    const char* status_str = "OK";
                    
                    if (status == NUMERIC_OK)
                    {
                        status_str = "OK";
                    }
                    else if (status == NUMERIC_FAIL)
                    {
                        status_str = "FAIL";
                    }
                    else if (status == NUMERIC_LOW_DATA)
                    {
                        status_str = "LOWDATA";
                    }
                    else if (status == NUMERIC_HI_DATA)
                    {
                        status_str = "HIDATA";
                    }
                    else
                    {
                        assert(false);
                    }
                    
                    fprintf(fout, "%s\t%s\t%d\t%lg\t%lg\t%lg\t%lg\t%s\t%s\n",
                            description.c_str(),
                            condition_name.c_str(),
                            rep_num,
                            raw_count,
                            internal_count,
                            external_count,
                            FPKM,
                            "-",
                            status_str);
                }
            }
            fflush(fout);
        }
    }
    
    void print_read_group_cuffdiff_info(FILE* fout,
                                        const vector<boost::shared_ptr<ReadGroupProperties> >& all_read_groups)
    {
        fprintf(fout, "file\tcondition\treplicate_num\ttotal_mass\tnorm_mass\tinternal_scale\texternal_scale\n");
        for (size_t i = 0; i < all_read_groups.size(); ++i)
        {
            boost::shared_ptr<ReadGroupProperties const> rg_props = all_read_groups[i];
            fprintf(fout, "%s\t%s\t%d\t%Lg\t%Lg\t%lg\t%lg\n",
                    rg_props->file_path().c_str(),
                    rg_props->condition_name().c_str(),
                    rg_props->replicate_num(),
                    rg_props->total_map_mass(),
                    rg_props->normalized_map_mass(),
                    rg_props->internal_scale_factor(),
                    rg_props->external_scale_factor());
            
        }
    }
    
    void print_read_group_simple_table_info(FILE* fout,
                                            const vector<boost::shared_ptr<ReadGroupProperties> >& all_read_groups)
    {
        //fprintf(fout, "file\tcondition\treplicate_num\ttotal_mass\tnorm_mass\tinternal_scale\texternal_scale\n");
        fprintf(fout, "sample_id\tfile\ttotal_mass\tinternal_scale\texternal_scale\tmedian_transcript_frags\n");
        for (size_t i = 0; i < all_read_groups.size(); ++i)
        {
            boost::shared_ptr<ReadGroupProperties const> rg_props = all_read_groups[i];
            fprintf(fout, "%s_%d\t%s\t%Lg\t%lg\t%lg\t%lg\n",
                    rg_props->condition_name().c_str(),
                    rg_props->replicate_num(),
                    rg_props->file_path().c_str(),
                    rg_props->total_map_mass(),
                    rg_props->internal_scale_factor(),
                    rg_props->external_scale_factor(),
                    rg_props->mode_transcript_coverage());
            
        }
    }

    void print_feature_attr_simple_table_header(FILE* fout,
                                         const FPKMTrackingTable& tracking)
    {
        fprintf(fout,"tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength");
        fprintf(fout, "\n");
    }

    
    void print_feature_attr_simple_table(const string& target_desc,
                                         FILE* fout,
                                         const FPKMTrackingTable& tracking)
    {
        FPKMTrackingTable::const_iterator itr = tracking.find(target_desc);
        if (itr != tracking.end())
        {
            const string& description = itr->first;
            const string& locus_tag = itr->second.locus_tag;
            const string& ref_match = itr->second.ref_match;
            int length = itr->second.length;
            char length_buff[33] = "-";
            if (length)
                sprintf(length_buff, "%d", length);
            const set<string>& gene_names = itr->second.gene_names;
            const set<string>& gene_ids = itr->second.gene_ids;
            const set<string>& tss_ids = itr->second.tss_ids;
            char class_code = itr->second.classcode ? itr->second.classcode : '-';
            string all_gene_names = cat_strings(gene_names);
            if (all_gene_names == "")
                all_gene_names = "-";
            string all_gene_ids = cat_strings(gene_ids);
            if (all_gene_ids == "")
                all_gene_ids = "-";
            string all_tss_ids = cat_strings(tss_ids);
            if (all_tss_ids == "")
                all_tss_ids = "-";
            fprintf(fout, "%s\t%c\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    description.c_str(),
                    class_code,
                    ref_match.c_str(),
                    all_gene_ids.c_str(),
                    all_gene_names.c_str(),
                    all_tss_ids.c_str(),
                    locus_tag.c_str(),
                    length_buff);
        }
    }
    
    void print_run_info(FILE* fout)
    {
        fprintf(fout, "param\tvalue\n");
        fprintf(fout, "cmd_line\t%s\n", cmd_str.c_str());
        fprintf(fout, "version\t%s\n", PACKAGE_VERSION);
        fprintf(fout, "SVN_revision\t%s\n",SVN_REVISION);
        fprintf(fout, "boost_version\t%d\n", BOOST_VERSION);
    }
    
    void write_header_output_cuffdiff_format()
    {
        // FPKM tracking
        
        FILE* fiso_fpkm_tracking =  _outfiles->isoform_fpkm_tracking_out;
        //fprintf(stderr, "Writing isoform-level FPKM tracking\n");
        print_FPKM_tracking_header(fiso_fpkm_tracking,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_fpkm_tracking =  _outfiles->tss_group_fpkm_tracking_out;
        //fprintf(stderr, "Writing TSS group-level FPKM tracking\n");
        print_FPKM_tracking_header(ftss_fpkm_tracking,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_fpkm_tracking =  _outfiles->gene_fpkm_tracking_out;
        //fprintf(stderr, "Writing gene-level FPKM tracking\n");
        print_FPKM_tracking_header(fgene_fpkm_tracking,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_fpkm_tracking =  _outfiles->cds_fpkm_tracking_out;
        //fprintf(stderr, "Writing CDS-level FPKM tracking\n");
        print_FPKM_tracking_header(fcds_fpkm_tracking,_tracking->cds_fpkm_tracking);
        
        // Count tracking
        
        FILE* fiso_count_tracking =  _outfiles->isoform_count_tracking_out;
        //fprintf(stderr, "Writing isoform-level count tracking\n");
        print_count_tracking_header(fiso_count_tracking,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_count_tracking =  _outfiles->tss_group_count_tracking_out;
        //fprintf(stderr, "Writing TSS group-level count tracking\n");
        print_count_tracking_header(ftss_count_tracking,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_count_tracking =  _outfiles->gene_count_tracking_out;
        //fprintf(stderr, "Writing gene-level count tracking\n");
        print_count_tracking_header(fgene_count_tracking,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_count_tracking =  _outfiles->cds_count_tracking_out;
        //fprintf(stderr, "Writing CDS-level count tracking\n");
        print_count_tracking_header(fcds_count_tracking,_tracking->cds_fpkm_tracking);
        
        // Read group tracking
        
        FILE* fiso_rep_tracking =  _outfiles->isoform_rep_tracking_out;
        //fprintf(stderr, "Writing isoform-level read group tracking\n");
        print_read_group_tracking_header(fiso_rep_tracking,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_rep_tracking =  _outfiles->tss_group_rep_tracking_out;
        //fprintf(stderr, "Writing TSS group-level read group tracking\n");
        print_read_group_tracking_header(ftss_rep_tracking,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_rep_tracking =  _outfiles->gene_rep_tracking_out;
        //fprintf(stderr, "Writing gene-level read group tracking\n");
        print_read_group_tracking_header(fgene_rep_tracking,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_rep_tracking =  _outfiles->cds_rep_tracking_out;
        //fprintf(stderr, "Writing CDS-level read group tracking\n");
        print_read_group_tracking_header(fcds_rep_tracking,_tracking->cds_fpkm_tracking);
        
        FILE* fread_group_info =  _outfiles->read_group_info_out;
        //fprintf(stderr, "Writing read group info\n");
        print_read_group_cuffdiff_info(fread_group_info,_all_read_groups);
        
        FILE* frun_info =  _outfiles->run_info_out;
        //fprintf(stderr, "Writing run info\n");
        print_run_info(frun_info);
        
    }
    
    void write_header_output_simple_table_format()
    {
        // FPKM tracking
        
        FILE* fiso_fpkm_tracking =  _outfiles->isoform_fpkm_tracking_out;
        //fprintf(stderr, "Writing isoform-level FPKM tracking\n");
        print_FPKM_simple_table_header(fiso_fpkm_tracking,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_fpkm_tracking =  _outfiles->tss_group_fpkm_tracking_out;
        //fprintf(stderr, "Writing TSS group-level FPKM tracking\n");
        print_FPKM_simple_table_header(ftss_fpkm_tracking,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_fpkm_tracking =  _outfiles->gene_fpkm_tracking_out;
        //fprintf(stderr, "Writing gene-level FPKM tracking\n");
        print_FPKM_simple_table_header(fgene_fpkm_tracking,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_fpkm_tracking =  _outfiles->cds_fpkm_tracking_out;
        //fprintf(stderr, "Writing CDS-level FPKM tracking\n");
        print_FPKM_simple_table_header(fcds_fpkm_tracking,_tracking->cds_fpkm_tracking);
        
        // Count tracking
        
        FILE* fiso_count_tracking =  _outfiles->isoform_count_tracking_out;
        //fprintf(stderr, "Writing isoform-level count tracking\n");
        print_count_simple_table_header(fiso_count_tracking,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_count_tracking =  _outfiles->tss_group_count_tracking_out;
        //fprintf(stderr, "Writing TSS group-level count tracking\n");
        print_count_simple_table_header(ftss_count_tracking,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_count_tracking =  _outfiles->gene_count_tracking_out;
        //fprintf(stderr, "Writing gene-level count tracking\n");
        print_count_simple_table_header(fgene_count_tracking,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_count_tracking =  _outfiles->cds_count_tracking_out;
        //fprintf(stderr, "Writing CDS-level count tracking\n");
        print_count_simple_table_header(fcds_count_tracking,_tracking->cds_fpkm_tracking);
        
        // We can also take care of all this metadata about samples and genes
        FILE* fiso_attr =  _outfiles->isoform_attr_out;
        //fprintf(stderr, "Writing isoform-level attributes\n");
        print_feature_attr_simple_table_header(fiso_attr,_tracking->isoform_fpkm_tracking);
        
        FILE* ftss_attr =  _outfiles->tss_group_attr_out;
        //fprintf(stderr, "Writing TSS group-level attributes\n");
        print_feature_attr_simple_table_header(ftss_attr,_tracking->tss_group_fpkm_tracking);
        
        FILE* fgene_attr =  _outfiles->gene_attr_out;
        //fprintf(stderr, "Writing gene-level attributes\n");
        print_feature_attr_simple_table_header(fgene_attr,_tracking->gene_fpkm_tracking);
        
        FILE* fcds_attr =  _outfiles->cds_attr_out;
        //fprintf(stderr, "Writing CDS-level attributes\n");
        print_feature_attr_simple_table_header(fcds_attr,_tracking->cds_fpkm_tracking);
        
        FILE* fread_group_info =  _outfiles->read_group_info_out;
        //fprintf(stderr, "Writing read group info\n");
        print_read_group_simple_table_info(fread_group_info,_all_read_groups);
        
        FILE* frun_info =  _outfiles->run_info_out;
        //fprintf(stderr, "Writing run info\n");
        print_run_info(frun_info);
        
    }
    
    void write_output(const vector<boost::shared_ptr<SampleAbundances> >& abundances)
    {
        if (output_format == CUFFDIFF_OUTPUT_FMT)
        {
            write_cuffdiff_output(abundances);
        }
        else if (output_format == SIMPLE_TABLE_OUTPUT_FMT)
        {
            write_simple_table_output(abundances);
        }
        else{
            fprintf(stderr, "Error: unrecognized output format!\n");
            exit(1);
        }
    }
    
    void write_header_output()
    {
        if (output_format == CUFFDIFF_OUTPUT_FMT)
        {
            write_header_output_cuffdiff_format();
        }
        else if (output_format == SIMPLE_TABLE_OUTPUT_FMT)
        {
            write_header_output_simple_table_format();
        }
        else{
            fprintf(stderr, "Error: unrecognized output format!\n");
            exit(1);
        }
    }
    
    void write_cuffdiff_output(const vector<boost::shared_ptr<SampleAbundances> >& abundances)
    {
        const AbundanceGroup& ab_group = abundances.front()->transcripts;
        //fprintf(stderr, "[%d] count = %lg\n",i,  ab_group.num_fragments());
        BOOST_FOREACH (boost::shared_ptr<Abundance> ab, ab_group.abundances())
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab->description(), ab->locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_tracking(ab->description(), _outfiles->isoform_fpkm_tracking_out, _tracking->isoform_fpkm_tracking);
                print_count_tracking(ab->description(), _outfiles->isoform_count_tracking_out, _tracking->isoform_fpkm_tracking);
                print_read_group_tracking(ab->description(), _outfiles->isoform_rep_tracking_out, _tracking->isoform_fpkm_tracking);
                remove_abundance_from_tracking_table(ab->description(), _tracking->isoform_fpkm_tracking);
            }
        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->cds)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_tracking(ab.description(), _outfiles->cds_fpkm_tracking_out, _tracking->cds_fpkm_tracking);
                print_count_tracking(ab.description(), _outfiles->cds_count_tracking_out, _tracking->cds_fpkm_tracking);
                print_read_group_tracking(ab.description(), _outfiles->cds_rep_tracking_out, _tracking->cds_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->cds_fpkm_tracking);
                
            }
        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->primary_transcripts)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_tracking(ab.description(), _outfiles->tss_group_fpkm_tracking_out, _tracking->tss_group_fpkm_tracking);
                print_count_tracking(ab.description(), _outfiles->tss_group_count_tracking_out, _tracking->tss_group_fpkm_tracking);
                print_read_group_tracking(ab.description(), _outfiles->tss_group_rep_tracking_out, _tracking->tss_group_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->tss_group_fpkm_tracking);
            }
        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->genes)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_tracking(ab.description(), _outfiles->gene_fpkm_tracking_out, _tracking->gene_fpkm_tracking);
                print_count_tracking(ab.description(), _outfiles->gene_count_tracking_out, _tracking->gene_fpkm_tracking);
                print_read_group_tracking(ab.description(), _outfiles->gene_rep_tracking_out, _tracking->gene_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->gene_fpkm_tracking);
            }
        }
    }
    
    
    void write_simple_table_output(const vector<boost::shared_ptr<SampleAbundances> >& abundances)
    {
        const AbundanceGroup& ab_group = abundances.front()->transcripts;
        //fprintf(stderr, "[%d] count = %lg\n",i,  ab_group.num_fragments());
        BOOST_FOREACH (boost::shared_ptr<Abundance> ab, ab_group.abundances())
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab->description(), ab->locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_simple_table(ab->description(), _outfiles->isoform_fpkm_tracking_out, _tracking->isoform_fpkm_tracking);
                print_count_simple_table(ab->description(), _outfiles->isoform_count_tracking_out, _tracking->isoform_fpkm_tracking);
                print_feature_attr_simple_table(ab->description(), _outfiles->isoform_attr_out, _tracking->isoform_fpkm_tracking);
                remove_abundance_from_tracking_table(ab->description(), _tracking->isoform_fpkm_tracking);
            }

        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->cds)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_simple_table(ab.description(), _outfiles->cds_fpkm_tracking_out, _tracking->cds_fpkm_tracking);
                print_count_simple_table(ab.description(), _outfiles->cds_count_tracking_out, _tracking->cds_fpkm_tracking);
                print_feature_attr_simple_table(ab.description(), _outfiles->cds_attr_out, _tracking->cds_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->cds_fpkm_tracking);
            }

        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->primary_transcripts)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_simple_table(ab.description(), _outfiles->tss_group_fpkm_tracking_out, _tracking->tss_group_fpkm_tracking);
                print_count_simple_table(ab.description(), _outfiles->tss_group_count_tracking_out, _tracking->tss_group_fpkm_tracking);
                print_feature_attr_simple_table(ab.description(), _outfiles->tss_group_attr_out, _tracking->tss_group_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->tss_group_fpkm_tracking);
            }
        }
        
        BOOST_FOREACH (AbundanceGroup& ab, abundances.front()->genes)
        {
            int num_remaining = _id_to_locus_map->unregister_locus_from_id(ab.description(), ab.locus_tag());
            if (num_remaining == 0)
            {
                print_FPKM_simple_table(ab.description(), _outfiles->gene_fpkm_tracking_out, _tracking->gene_fpkm_tracking);
                print_count_simple_table(ab.description(), _outfiles->gene_count_tracking_out, _tracking->gene_fpkm_tracking);
                print_feature_attr_simple_table(ab.description(), _outfiles->gene_attr_out, _tracking->gene_fpkm_tracking);
                remove_abundance_from_tracking_table(ab.description(), _tracking->gene_fpkm_tracking);
            }
        }
    }
};


extern double min_read_count;

void sample_worker(bool non_empty,
                   boost::shared_ptr<HitBundle> bundle,
                   const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   boost::shared_ptr<SampleAbundances> abundance,
                   size_t factory_id,
                   boost::shared_ptr<TestLauncher> launcher,
                   bool calculate_variance);

void test_differential(const string& locus_tag,
					   const vector<boost::shared_ptr<SampleAbundances> >& samples,
                       const vector<pair<size_t, size_t> >& constrasts,
					   Tests& tests,
					   Tracking& tracking);

void dump_locus_variance_info(const string& filename);

#if ENABLE_THREADS
void decr_pool_count();
extern boost::mutex locus_thread_pool_lock;
extern int locus_curr_threads;
extern int locus_num_threads;
#endif

#endif

void validate_cross_sample_parameters(const vector<boost::shared_ptr<ReadGroupProperties> >& all_read_groups);
