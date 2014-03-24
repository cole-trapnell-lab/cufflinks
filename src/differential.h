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

struct TestLauncher
{
private:
    TestLauncher(TestLauncher& rhs) {}
    
public:
    TestLauncher(int num_samples,
                 const vector<pair<size_t, size_t> >& contrasts,
                 Tests* tests,
                 Tracking* tracking,
                 ProgressBar* p_bar) 
    :
    _orig_workers(num_samples),
    _contrasts(contrasts),
    _tests(tests),
    _tracking(tracking),
    _p_bar(p_bar)
    {
    }
    
    void operator()();
    
    void register_locus(const string& locus_id);
    void abundance_avail(const string& locus_id, 
                         boost::shared_ptr<SampleAbundances> ab, 
                         size_t factory_id);
    void test_finished_loci();
    void perform_testing(vector<boost::shared_ptr<SampleAbundances> > abundances);
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
