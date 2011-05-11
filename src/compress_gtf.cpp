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
#include <algorithm>

#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>

#include "common.h"
#include "hits.h"
#include "bundles.h"

#include "gtf_tracking.h"
#include "scaffolds.h"
#include "tokenize.h"
#include "genes.h"

using namespace boost;
using namespace std;

#if ENABLE_THREADS
const char *short_options = "r:F";
#else
const char *short_options = "r:F";
#endif

bool raw_fpkm = false;
bool proj_union = false;
bool proj_intersection = false;

static struct option long_options[] = {
{"reference-seq",		required_argument,		 0,			 'r'},
{"raw-fpkm",            no_argument,             0,			 'F'},
{"union",               no_argument,             0,			 'U'},
{"intersection",        no_argument,             0,			 'I'},


{0, 0, 0, 0} // terminator
};

void print_usage()
{
	//NOTE: SPACES ONLY, bozo
	fprintf(stderr, "compress_gtf v%s\n", PACKAGE_VERSION);
	fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
	fprintf(stderr, "-----------------------------\n"); 
	fprintf(stderr, "Usage:   compress_gtf [options] <reference.gtf> <compressed_reference.gtf>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-r/--reference-seq			  reference fasta file                     [ default:   NULL ]\n");
    fprintf(stderr, "-F/--raw-fpkm			      use FPKM instead of isoform fraction                        \n");
    fprintf(stderr, "-U/--union                   report projective union                  [ default:   OFF  ]\n");
    fprintf(stderr, "-I/--intersection            report projective intersection           [ default:   ON   ]\n");
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
            case 'U':
			{
				proj_union = true;
				break;
            } 
            case 'I':
			{
				proj_intersection = true;
				break;
            }  
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);

    if (proj_union && proj_intersection)
    {
        fprintf (stderr, "Error: please specify only one of --union and --intersection");
        exit(1);
    }
    
//    if (!proj_union && !proj_intersection)
//        proj_intersection =  true;
	return 0;
}

void compress_genes(FILE* ftranscripts,
                    RefSequenceTable& rt,
                    vector<shared_ptr<Scaffold> >& ref_mRNAs)
{
    adjacency_list <vecS, vecS, undirectedS> G;
    
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		add_vertex(G);
	}
	
    for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
        shared_ptr<Scaffold> scaff_i = ref_mRNAs[i];
        for (size_t j = 0; j < ref_mRNAs.size(); ++j)
        {
            shared_ptr<Scaffold> scaff_j = ref_mRNAs[j];
			if (scaff_i->annotated_gene_id() == scaff_j->annotated_gene_id())
				add_edge(i, j, G);
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
        vector<Scaffold> gene_scaffs;
        string gene_id;
        foreach (shared_ptr<Scaffold> s, gene)
        {
            if (gene_id == "")
                gene_id = s->annotated_gene_id();
            
            gene_scaffs.push_back(*s);
        }
        
        if (gene_scaffs.empty())
            continue;
        
        next_gene_id++;
        
        Scaffold smashed_gene;
        if (!proj_intersection && !proj_union)
        {
            foreach (shared_ptr<Scaffold> s, gene)
            {
                /*
                 *transfrag,
                 gene_id,
                 (int)isoforms.size() + 1,
                 FPKM,
                 iso_ab->effective_length(),
                 iso_ab->gamma(),
                 iso_ab->FPKM_conf(),
                 density_per_bp, 
                 estimated_count,
                 density_score,
                 iso_ab->status(),
                 ref_gene_id)*/
                
                Isoform iso(*s,
                            -1,
                            1,
                            0.0,
                            s->length(),
                            0.0,
                            ConfidenceInterval(0.0,0.0),
                            0,
                            0,
                            0,
                            NUMERIC_OK,
                            gene_id);
                vector<string> isoform_exon_recs;
                
                iso.get_gtf(isoform_exon_recs, rt);
                
                for (size_t g = 0; g < isoform_exon_recs.size(); ++g)
                {
                    fprintf(ftranscripts, "%s", isoform_exon_recs[g].c_str());
                }
            }
        }
        else
        {
            if (proj_union)
                Scaffold::merge(gene_scaffs, smashed_gene, false);
            else if (proj_intersection)
            {
                vector<AugmentedCuffOp> iso_ops;
                
                int gmax = -1;
                int gmin = numeric_limits<int>::max();
                
                foreach (shared_ptr<Scaffold> s, gene)
                {
                    //iso_ops.push_back(s->augmented_ops());
                    //sort (iso_ops.back().begin(), iso_ops.back().end());
                    if (s->left() < gmin)
                        gmin = s->left();
                    if (s->right() > gmax)
                        gmax = s->right();
                }
                
                foreach (shared_ptr<Scaffold> s, gene)
                {
                    if (s->left() > gmin)
                    {
                        iso_ops.push_back(AugmentedCuffOp(CUFF_INTRON, gmin, s->left() - gmin)); 
                    }
                    if (s->right() < gmax)
                    {
                        iso_ops.push_back(AugmentedCuffOp(CUFF_INTRON, s->right(), gmax - s->right())); 
                    }
                    iso_ops.insert(iso_ops.end(), s->augmented_ops().begin(), s->augmented_ops().end());
                }
//                vector<AugmentedCuffOp> intersect = iso_ops.front();
//                for (size_t j = 1; j < iso_ops.size(); ++j)
//                {
//                    vector<AugmentedCuffOp> tmp;
//                    const vector<AugmentedCuffOp>& iso_ops_j = iso_ops[j];
//                    //set_intersection(intersect.begin(), intersect.end(), iso_ops_j.begin(), iso_ops_j.end(), back_inserter(tmp));
//                    intersect.insert(intersect.end(), iso_ops_j.begin(), iso_ops_j.end());
//                    
//                    intersect.push_back(
//                    assert (tmp.size() <= intersect.size());
//                    //intersect = tmp;
//                    //sort(intersect.begin(), intersect.end());
//                }
//                
                sort(iso_ops.begin(), iso_ops.end(), AugmentedCuffOp::g_left_lt);
//                
//                while (!intersect.empty() && intersect.front().opcode != CUFF_MATCH)
//                {
//                    intersect.erase(intersect.begin());
//                }
//                
//                while (!intersect.empty() && intersect.back().opcode != CUFF_MATCH)
//                {
//                    intersect.pop_back();
//                }
//                
//                if (intersect.empty())
//                    continue;
                
                vector<AugmentedCuffOp> merged_ops;
                AugmentedCuffOp::merge_ops(iso_ops, merged_ops, true, true);
                vector<AugmentedCuffOp>::iterator first_match = merged_ops.begin();
                vector<AugmentedCuffOp>::iterator last_match = merged_ops.end();
                last_match--;
                while(first_match < merged_ops.end())
                {
                    if (first_match->opcode == CUFF_MATCH)
                        break;
                    first_match++;
                }
                while(last_match >= merged_ops.begin() && last_match< merged_ops.end())
                {
                    if (last_match->opcode == CUFF_MATCH)
                        break;
                    last_match--;
                }
                
                vector<AugmentedCuffOp> internal_matches;
                if (last_match >= first_match && last_match < merged_ops.end())
                {
                    last_match++;
                    
                    internal_matches.insert(internal_matches.end(), first_match, last_match);
                    smashed_gene = Scaffold(gene.front()->ref_id(), gene.front()->strand(), internal_matches);
                }
                else
                {
                    
                    fprintf(stderr, "Could not find consitutive region for %s\n", gene_id.c_str());
                    continue;
                }
                
            }
            else
                assert(false);
            assert (smashed_gene.ref_id());
            
            Isoform iso(smashed_gene,
                        -1,
                        1,
                        0.0,
                        smashed_gene.length(),
                        0.0,
                        ConfidenceInterval(0.0,0.0),
                        0, 
                        0,
                        0,
                        NUMERIC_OK,
                        gene_id);
            vector<string> isoform_exon_recs;
            
            iso.get_gtf(isoform_exon_recs, rt);
            
            for (size_t g = 0; g < isoform_exon_recs.size(); ++g)
            {
                fprintf(ftranscripts, "%s", isoform_exon_recs[g].c_str());
            }
        }
        
        fflush(ftranscripts);
    }
}

void driver(vector<FILE*> ref_gtf_files, FILE* gtf_out)
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
            compress_genes(gtf_out, rt, ref_mRNAs);
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
	
    string gtf_out_filename = argv[optind++];
    
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
    
    FILE* gtf_out = NULL;
	if (gtf_out_filename != "")
	{
		gtf_out = fopen(gtf_out_filename.c_str(), "w");
		if (!gtf_out)
		{
			fprintf(stderr, "Error: cannot open GTF file %s for writing\n",
					gtf_out_filename.c_str());
			exit(1);
		}
	}
    
    driver(ref_gtf_files, gtf_out);
	
	return 0;
}
