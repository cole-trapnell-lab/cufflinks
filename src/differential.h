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
#include <vector>
#include <string>

#include <boost/thread.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "abundances.h"

using namespace std;

enum TestStatus {
	NOTEST,  // successful calculation, test not performed
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
	vector<FILE*> isoform_de_outfiles;
	vector<FILE*> group_de_outfiles;
	vector<FILE*> gene_de_outfiles;
	vector<FILE*> cds_de_outfiles;
	
	vector<FILE*> diff_splicing_outfiles;
	vector<FILE*> diff_promoter_outfiles;
	vector<FILE*> diff_cds_outfiles;
	
	FILE* isoform_fpkm_tracking_out;
	FILE* tss_group_fpkm_tracking_out;
	FILE* gene_fpkm_tracking_out;
	FILE* cds_fpkm_tracking_out;
};

struct Tests
{
	vector<SampleDiffs> isoform_de_tests;
	vector<SampleDiffs> tss_group_de_tests;
	vector<SampleDiffs> gene_de_tests;
	vector<SampleDiffs> cds_de_tests;
	
	vector<SampleDiffs> diff_splicing_tests; // to be performed on the isoforms of a single tss group
	vector<SampleDiffs> diff_promoter_tests; // to be performed on the tss groups of a single gene
	vector<SampleDiffs> diff_cds_tests; // to be performed on the cds groups of a single gene
};

struct FPKMContext
{
	FPKMContext(double c, double r, double v, double m)
		: counts(c), FPKM(r), FPKM_variance(v), log_sample_mass(m) {}
	double counts;
	double FPKM;
	double FPKM_variance;
	double log_sample_mass;
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

extern double min_read_count;


#if 0
// This performs a between-group test on an isoform or TSS grouping, on two 
// different samples.
bool test_diffexp(size_t sample1,
				  size_t sample2,
				  size_t transcript_group_id,
				  const FPKMContext& curr,
				  const FPKMContext& prev,
				  SampleDifference& test);

// This performs between-group tests on isoforms or TSS groupings in a single
// locus, on two different samples.
int get_de_tests(size_t curr_sample_idx,
				 const QuantGroup& curr_quant_group,
				 long double curr_mass,
				 size_t prev_sample_idx,
				 const QuantGroup& prev_quant_group,
				 long double prev_mass,
				 SampleDiffs& isoform_de_tests,
				 bool enough_reads);

// This performs within-group tests on a set of isoforms or a set of TSS groups.
// This is a way of looking for meaningful differential splicing or differential
// promoter use.
void get_ds_tests(size_t curr_sample_idx,
				  const QuantGroup& curr_quant_group,
				  size_t prev_sample_idx,
				  const QuantGroup& prev_quant_group,
				  SampleDiffs& diff_tests,
				  bool enough_reads);

#endif

void test_differential(const RefSequenceTable& rt, 
					   const vector<HitBundle*>& sample_bundles,
					   const vector<long double>& sample_masses,
                       const vector<EmpDist>& frag_len_dists,
					   Tests& tests,
					   Tracking& tracking);


#endif
