/*
 *  gtf_to_sam.cpp
 *  Cufflinks
 *
 *  Created by Cole Trapnell on 8/1/10.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#endif

#include <stdlib.h>
#include <getopt.h>
#include <string>

#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>

#include "common.h"
#include "hits.h"
#include "bundles.h"

#include "scaffolds.h"
#include "tokenize.h"

using namespace boost;
using namespace std;

#if ENABLE_THREADS
const char *short_options = "r:F";
#else
const char *short_options = "r:F";
#endif

bool raw_fpkm = false;

static struct option long_options[] = {
{"reference-seq",		required_argument,		 0,			 'r'},
{"raw-fpkm",            no_argument,             0,			 'F'},
{0, 0, 0, 0} // terminator
};

void print_usage()
{
	//NOTE: SPACES ONLY, bozo
	fprintf(stderr, "gtf_to_sam v%s\n", PACKAGE_VERSION);
	fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
	fprintf(stderr, "-----------------------------\n"); 
	fprintf(stderr, "Usage:   cufflinks [options] <transcripts1.gtf,...,transcriptsN.gtf> <out.sam>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-r/--reference-seq			  reference fasta file                     [ default:   NULL ]\n");
    fprintf(stderr, "-F/--raw-fpkm			      use FPKM instead of isoform fraction                        \n");
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
            case 'r':
			{
				fasta_dir = optarg;
				break;
            }    
            case 'F':
			{
				raw_fpkm = true;
				break;
            }   
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);

	return 0;
}

void print_scaff_as_sam(FILE* sam_out,
                        const RefSequenceTable& rt, 
                        const Scaffold& scaff)
{
	string seq;
	string quals;

    seq = "*";
    quals = "*";	

	uint32_t sam_flag = 0;
	if (scaff.strand() == CUFF_REV)
	{
		sam_flag |= 0x0010; // BAM_FREVERSE
//		if (sequence)
//		{
//			reverse_complement(seq);
//			reverse(quals.begin(), quals.end());
//		}
	}
	
	uint32_t sam_pos = scaff.left() + 1;
	uint32_t map_quality = 255;
	char cigar[8192];
	cigar[0] = 0;
	string mate_ref_name = "*";
	uint32_t mate_pos = 0;
	uint32_t insert_size = 0;
	
	const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
	for (size_t c = 0; c < ops.size(); ++c)
	{
		char ibuf[64];
		sprintf(ibuf, "%d", ops[c].genomic_length);
		switch(ops[c].opcode)
		{
			case CUFF_MATCH:
				strcat(cigar, ibuf);
				strcat(cigar, "M");
				break;
			case CUFF_INTRON:
				strcat(cigar, ibuf);
				strcat(cigar, "N");
				break;
			default:
                fprintf(stderr, "Warning: Transcript %s contains an unconvertible alignment operator, skipping\n", scaff.annotated_trans_id().c_str());
                return;
                break;
		}
	}
	
	//string q = string(bh.read_len(), '!');
	//string s = string(bh.read_len(), 'N');
	
    const char* ref_name = rt.get_name(scaff.ref_id());
    if (!ref_name)
    {
        fprintf(stderr, "Warning: Could not find contig name for ID %d, skipping\n", scaff.ref_id());
        return;
    }
    
    if (scaff.annotated_trans_id() == "")
    {
        fprintf(stderr, "Warning: transcript_id attribute is empty, skipping\n");
        return;
    }
    
	fprintf(sam_out,
			"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			scaff.annotated_trans_id().c_str(),
			sam_flag,
			ref_name,
			sam_pos,
			map_quality,
			cigar,
			mate_ref_name.c_str(),
			mate_pos,
			insert_size,
			seq.c_str(),
			quals.c_str());
	
	if (scaff.strand() != CUFF_STRAND_UNKNOWN)
	{
		fprintf(sam_out,
				"\tXS:A:%c",
				scaff.strand() == CUFF_REV ? '-' : '+');
	}
    
    if (scaff.fpkm() != 0)
	{
		fprintf(sam_out,
				"\tZF:f:%f",
				scaff.fpkm());
	}
	
	fprintf(sam_out, "\n");
    
}
	
void set_relative_fpkms(vector<shared_ptr<Scaffold> >& ref_mRNAs)
{
    adjacency_list <vecS, vecS, undirectedS> G;
	
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		add_vertex(G);
	}
	
    map<string, vector<int> > gene_id_idxs;
    
    for (size_t i = 0; i < ref_mRNAs.size(); ++i)
    {
        pair<map<string, vector<int> >::iterator, bool> inserted;
        inserted = gene_id_idxs.insert(make_pair(ref_mRNAs[i]->annotated_gene_id(), vector<int>()));
        inserted.first->second.push_back(i);
    }
    
    for (map<string, vector<int> >::iterator itr = gene_id_idxs.begin();
         itr != gene_id_idxs.end();
         ++itr)
    {
        vector<int>& gene = itr->second;
        for (size_t i = 0; i < gene.size(); ++i)
        {
            for (size_t j = 0; j < gene.size(); ++j)
            {
                {
                    add_edge(gene[i], gene[j], G);
                }
            }
        }
    }
    
    std::vector<int> component(num_vertices(G));
	connected_components(G, &component[0]);
	
	vector<vector<bool> > clusters(ref_mRNAs.size(), 
								   vector<bool>(ref_mRNAs.size(), false));
	
	//vector<vector<size_t> > cluster_indices(three_prime_ends.size());
    
    vector<vector<shared_ptr<Scaffold> > > grouped_scaffolds(ref_mRNAs.size());
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		clusters[component[i]][i] = true;
		grouped_scaffolds[component[i]].push_back(ref_mRNAs[i]);
	}
    
    for (size_t i = 0; i < grouped_scaffolds.size(); ++i)
    {
        vector<shared_ptr<Scaffold> >& gene = grouped_scaffolds[i];
        
        double total_fpkm = 0.0;
        foreach(shared_ptr<Scaffold> scaff, gene)
        {
            total_fpkm += scaff->fpkm();
        }
        if (total_fpkm > 0)
        {
            foreach (shared_ptr<Scaffold> scaff, gene)
            {
                scaff->fpkm(scaff->fpkm() / total_fpkm);
            }
        }
    }
}

void driver(vector<FILE*> ref_gtf_files, FILE* sam_out)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<vector<shared_ptr<Scaffold> > > ref_mRNA_table;
	vector<pair<string, vector<double> > > sample_count_table;
    
    foreach (FILE* ref_gtf, ref_gtf_files)
    {
        vector<shared_ptr<Scaffold> > ref_mRNAs;
        ::load_ref_rnas(ref_gtf, rt, ref_mRNAs, false, true);
        ref_mRNA_table.push_back(ref_mRNAs);
    }
    
    for (size_t j = 0; j < ref_mRNA_table.size(); ++j)
    {
        vector<shared_ptr<Scaffold> > ref_mRNAs = ref_mRNA_table[j];
        
        if (!raw_fpkm)
            set_relative_fpkms(ref_mRNAs);
        
        for (size_t i = 0; i < ref_mRNAs.size(); ++i)
        {
            print_scaff_as_sam(sam_out, rt, *ref_mRNA_table[j][i]);
        }
    }
}

int main(int argc, char** argv)
{
    init_library_table();
    
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string ref_gtf_in_filenames = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string sam_out_filename = argv[optind++];
    
    vector<string> ref_gtf_filenames;
    tokenize(ref_gtf_in_filenames, ",", ref_gtf_filenames);
    
    vector<FILE*> ref_gtf_files;
    
    foreach (const string& ref_gtf_in_filename, ref_gtf_filenames)
    {
        FILE* ref_gtf = NULL;
        if (ref_gtf_in_filename != "")
        {
            ref_gtf = fopen(ref_gtf_in_filename.c_str(), "r");
            if (!ref_gtf)
            {
                fprintf(stderr, "Error: cannot open GTF file %s for reading\n",
                        ref_gtf_in_filename.c_str());
                exit(1);
            }
            ref_gtf_files.push_back(ref_gtf);
        }
    }
    
    FILE* sam_out = NULL;
	if (sam_out_filename != "")
	{
		sam_out = fopen(sam_out_filename.c_str(), "w");
		if (!sam_out)
		{
			fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
					sam_out_filename.c_str());
			exit(1);
		}
	}
    
    driver(ref_gtf_files, sam_out);
	
	return 0;
}
