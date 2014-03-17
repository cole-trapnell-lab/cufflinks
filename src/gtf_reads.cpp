/*
 *  reads_of_interest.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 11/29/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#if ENABLE_THREADS
const char *short_options = "";
#else
const char *short_options = "";
#endif

static struct option long_options[] = {
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   cuffdiff <transcripts.gtf> <sample_hits.sam>\n");
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
			
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	return 0;
}

void driver(FILE* ref_gtf, FILE* sam_hit_file)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<Scaffold> ref_mRNAs;
	
	SAMHitFactory hs(it, rt);
    
	boost::shared_ptr<HitFactory> hit_factory;
    
    try
	{
		hit_factory = boost::shared_ptr<BAMHitFactory>(new BAMHitFactory(hit_file_name, it, rt));
	}
	catch (std::runtime_error& e)
	{
		fprintf(stderr, "File %s doesn't appear to be a valid BAM file, trying SAM...\n",
				hit_file_name.c_str());
        
        try
        {
            hit_factory = boost::shared_ptr<SAMHitFactory>(new SAMHitFactory(hit_file_name, it, rt));
        }
        catch (std::runtime_error& e)
        {
            fprintf(stderr, "Error: cannot open alignment file %s for reading\n",
                    hit_file_name.c_str());
            exit(1);
        }
	}
	BundleFactory& bundle_factory = *(new BundleFactory(hit_factory, bundle_mode));
    
    boost::crc_32_type ref_gtf_crc_result;
    vector<boost::shared_ptr<Scaffold> > ref_mRNAs;
    if (ref_gtf)
    {
        ::load_ref_rnas(ref_gtf, bundle_factory.ref_table(), ref_mRNAs, ref_gtf_crc_result, false, false);
        bundle_factory.set_ref_rnas(ref_mRNAs);
    }
    rt.print_rec_ordering();
    vector<boost::shared_ptr<Scaffold> > mask_rnas;
    boost::crc_32_type mask_gtf_crc_result;
    if (mask_gtf)
    {
        ::load_ref_rnas(mask_gtf, bundle_factory.ref_table(), mask_rnas, mask_gtf_crc_result, false, false);
        bundle_factory.set_mask_rnas(mask_rnas);
    }
    
	long num_fragments = 0;
	long num_reads = 0;
	
    const vector<MateHit>& hits = bundle.hits();
    for (size_t i = 0; i < bundle.hits().size(); ++i)
    {
        bool compatible = false;
        for (vector<Scaffold>::iterator ri = ref_mrnas.begin();
             ri != ref_mrnas.end();
             ++ri)
        {
            if (Scaffold::overlap_in_genome(*ri, hits[i], 0) &&
                Scaffold::compatible(*ri,hits[i]))
            {
                compatible = true;
                break;
            }
        }
        if (!compatible)
            continue;
        
        if (hits[i].left_alignment() || hits[i].right_alignment())
            num_fragments++;
        if (hits[i].left_alignment())
        {
            printf("%s\n", hits[i].left_alignment()->hitfile_rec().c_str());
            num_reads++;
        }
        if (hits[i].right_alignment())
        {
            printf("%s\n", hits[i].right_alignment()->hitfile_rec().c_str());
            num_reads++;
        }
    }
	
	fprintf(stderr, "Extracted %ld fragments, %ld reads\n", num_fragments, num_reads); 
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
	
    string ref_gtf_filename = argv[optind++];

	string sam_hits_file_name = argv[optind++];
	// Open the approppriate files
	FILE* sam_hits_file = fopen(sam_hits_file_name.c_str(), "r");
	if (sam_hits_file == NULL)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
				sam_hits_file_name.c_str());
		exit(1);
	}
	
	FILE* ref_gtf = NULL;
	if (ref_gtf_filename != "")
	{
		ref_gtf = fopen(ref_gtf_filename.c_str(), "r");
		if (!ref_gtf)
		{
			fprintf(stderr, "Error: cannot open GTF file %s for reading\n",
					ref_gtf_filename.c_str());
			exit(1);
		}
	}
	
    driver(ref_gtf, sam_hits_file);
	
	return 0;
}
