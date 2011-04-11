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

using namespace std;

enum TestStatus {
	NOTEST,  // successful calculation, test not performed
    LOWDATA, // unsuccessful calculation due to low data, test not performed
	OK,      // successful numerical calc, test performed
	FAIL     // numerical exception, test not performed
}; 

// Stores the differential expression of an isoform or set of isoforms in two
// different samples, along with a significance test statistic for the difference.
struct SampleDifference
{
	SampleDifference() :
	sample_1(-1), 
	sample_2(-1), 
	value_1(0.0),
	value_2(0.0),
	test_stat(0.0),
    p_value(1.0),
	corrected_p(1.0),
	tested_group_id(-1),
	test_status(NOTEST),
	significant(false){}
	
	size_t sample_1;
	size_t sample_2;
	
	double value_1;
	double value_2;
	double differential;
	double test_stat;
	double p_value;
	double corrected_p;
	
	size_t tested_group_id; // which scaffolds' FPKMs contribute
	
	string locus_desc;
	set<string> gene_names;
	set<string> protein_ids;
	string description; // isoforms or tss groups (e.g.) involved in this test
	
	TestStatus test_status;
	bool significant;
};

typedef map<string, SampleDifference > SampleDiffs;

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

struct FPKMContext
{
	FPKMContext(double c, double r, double v, AbundanceStatus s)
		: counts(c), FPKM(r), FPKM_variance(v), status(s) {}
	double counts;
	double FPKM;
	double FPKM_variance;
    AbundanceStatus status;
};

struct FPKMTracking
{
	string locus_tag;
	char classcode;
	set<string> tss_ids; // for individual isoforms only
	set<string> gene_names;
	set<string> protein_ids;
	string description; // isoforms or tss groups (e.g.) involved in this test
	string ref_match;
	
	TestStatus test_status;
	
	vector<FPKMContext> fpkm_series;
};

typedef map<string,  FPKMTracking> FPKMTrackingTable;

struct Tracking
{
	FPKMTrackingTable isoform_fpkm_tracking;
	FPKMTrackingTable tss_group_fpkm_tracking;
	FPKMTrackingTable gene_fpkm_tracking;
	FPKMTrackingTable cds_fpkm_tracking;
};

struct SampleAbundances
{
    string locus_tag;
	AbundanceGroup transcripts;
	vector<AbundanceGroup> primary_transcripts;
	vector<AbundanceGroup> gene_primary_transcripts;
	vector<AbundanceGroup> cds;
	vector<AbundanceGroup> gene_cds;
	vector<AbundanceGroup> genes;
	double cluster_mass;
};

extern double min_read_count;

void sample_worker(const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   shared_ptr<SampleAbundances> abundance,
                   shared_ptr<bool> non_empty);

void test_differential(const string& locus_tag,
					   const vector<shared_ptr<SampleAbundances> >& samples,
					   Tests& tests,
					   Tracking& tracking,
                       bool samples_are_time_series);

#if ENABLE_THREADS
void decr_pool_count();
#endif

#endif
