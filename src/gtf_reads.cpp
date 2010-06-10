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

// like a normal bundle factory, except the user must explicitly ask for 
// a range of alignments for the bunde.  NOTE: this factory does NOT seek 
// backwards - it is up to you to ask for monotonically increasing loci
// Attaches SAM alignment strings to each hit.
class FilteringBundleFactory
	{
	public:
		FilteringBundleFactory(SAMHitFactory& fac, FILE* hfile)
		: sam_hit_fac(fac), hit_file(hfile), _next_line_num(0) {}
		
		bool next_bundle(HitBundle& bundle_out, 
						 RefID ref_id, 
						 int left_boundary, 
						 int right_boundary);
		
		SAMHitFactory& hit_factory() { return sam_hit_fac; } 
		
		void reset() { rewind(hit_file); _next_line_num = 0; }
		
	private:
		SAMHitFactory sam_hit_fac;
		FILE* hit_file;
		int _next_line_num;
	};

bool FilteringBundleFactory::next_bundle(HitBundle& bundle_out, 
										 RefID ref_id, 
										 int left_boundary, 
										 int right_boundary)
{
	HitBundle bundle = bundle_out;
	
	if (feof(hit_file))
	{
		return false;
	}
	char bwt_buf[2048];
	
	RefID last_ref_id_seen = 0;
	int last_pos_seen = 0;
	
	off_t curr_pos = ftello(hit_file);
	
	
	while (fgets(bwt_buf, 2048, hit_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		_next_line_num++;
		
		shared_ptr<ReadHit> bh(new ReadHit());
		
		if (!sam_hit_fac.get_hit_from_buf(_next_line_num, bwt_buf, *bh, false))
		{
			continue;
		}
		
		if (bh->ref_id() == 84696373) // corresponds to SAM "*" under FNV hash. unaligned read record  
			continue;
		
		bh->hitfile_rec(bwt_buf);
		
		bool hit_within_boundary = false;
		
		if (bh->ref_id() != ref_id)
		{
			RefSequenceTable& rt = sam_hit_fac.ref_table();
			const char* gtf_name = rt.get_name(ref_id);
			const char* sam_name = rt.get_name(bh->ref_id());
			if (!gtf_name)
			{
				// The GTF name isn't even inthe SAM file, we'll never find
				// records by advancing the SAM.  Give up.  hit_within_boundary
				// will remain false, and we'll break out of the loop and 
				// reset the SAM file pointer
			}
			else
			{
				// the LocusBundleFactory's ref table should be populated with all 
				// values from the SAM
				assert (sam_name); 
				int c = strcmp(gtf_name, sam_name);
				
				if (c > 0)
				{
					// we need to keep advancing the SAM file, to catch up
					// with the GTF file
					continue;
				}
				else
				{
					// else the GTF file is before the SAM file, we'll never find
					// records by advancing the SAM.  Give up.  hit_within_boundary
					// will remain false, and we'll break out of the loop and 
					// reset the SAM file pointer
				}
			}
		}
		else
		{
			if (bh->right() < left_boundary)
				continue;
			
			if (bh->error_prob() > max_phred_err_prob)
				continue;
			
			if (bh->left() <= right_boundary)
				hit_within_boundary = true;
		}
		
		if (hit_within_boundary)
		{
			if (bh->left() < last_pos_seen)
			{
				fprintf(stderr, "Error: this SAM file doesn't appear to be correctly sorted!\n");
				fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
						sam_hit_fac.ref_table().get_name(bh->ref_id()),
						bh->left(),
						sam_hit_fac.ref_table().get_name(last_ref_id_seen),
						last_pos_seen);
				
				exit(1);
			}
			
			bundle.add_open_hit(bh);
		}
		else
		{
			fseeko(hit_file, curr_pos, SEEK_SET);
			break;
		}
		
		last_ref_id_seen = bh->ref_id();
		last_pos_seen = bh->left();
		
		curr_pos = ftello(hit_file);
	}
	
	bundle.finalize_open_mates();
	bundle_out = bundle;
	bundle_out.finalize();
	//assert (!bundle_out.ref_scaffolds().empty());
	return true;
}


void driver(FILE* ref_gtf, FILE* sam_hit_file)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<Scaffold> ref_mRNAs;
	load_ref_rnas(ref_gtf, rt, ref_mRNAs);
	if (ref_mRNAs.empty())
		return;
	
	SAMHitFactory hs(it, rt);
	FilteringBundleFactory lf(hs, sam_hit_file);
	HitBundle bundle;
	long num_fragments = 0;
	long num_reads = 0;
	
	vector<Scaffold>::iterator curr_ref_scaff = ref_mRNAs.begin();
	RefID last_ref_seen = curr_ref_scaff->ref_id();
	int left_boundary = curr_ref_scaff->left();
	int right_boundary = curr_ref_scaff->right();
	vector<Scaffold>::iterator ref_scaff_bundle_start = curr_ref_scaff;
	
	while (curr_ref_scaff != ref_mRNAs.end())
	{
		if (curr_ref_scaff->left() > right_boundary ||
			curr_ref_scaff->ref_id() != last_ref_seen)
		{
			
			vector<Scaffold> ref_mrnas;
			
			ref_mrnas.insert(ref_mrnas.end(), 
							 ref_scaff_bundle_start, 
							 curr_ref_scaff);
			
			HitBundle bundle;
			lf.next_bundle(bundle, last_ref_seen, left_boundary, right_boundary);
			
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
			
			left_boundary = curr_ref_scaff->left();
			if (curr_ref_scaff->ref_id() != last_ref_seen)
			{
				right_boundary = curr_ref_scaff->right();
			}
			last_ref_seen = curr_ref_scaff->ref_id();
			ref_scaff_bundle_start = curr_ref_scaff;
		}
		right_boundary = max(right_boundary, curr_ref_scaff->right());
		++curr_ref_scaff;
	}
	
	if (ref_scaff_bundle_start != ref_mRNAs.end())
	{
		vector<Scaffold> ref_mrnas;
		
		ref_mrnas.insert(ref_mrnas.end(), 
						 ref_scaff_bundle_start, 
						 curr_ref_scaff);
		
		HitBundle bundle;
		lf.next_bundle(bundle, last_ref_seen, left_boundary, right_boundary);
		
		const vector<MateHit>& hits = bundle.hits();
		for (size_t i = 0; i < bundle.hits().size(); ++i)
		{
			bool compatible = false;
			for (vector<Scaffold>::iterator ri = ref_mrnas.begin();
				 ri != ref_mrnas.end();
				 ++ri)
			{
				
				if (Scaffold::overlap_in_genome(*ri, hits[i], 0) &&
					Scaffold::compatible(*ri, hits[i]))
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
	}
	
	lf.reset();
	
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
