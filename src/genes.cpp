/*
 *  genes.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */
#include <boost/thread.hpp>
#include "genes.h"

using namespace boost;

#if ENABLE_THREADS
mutex gene_id_lock;
#endif

int next_isoform_id = 1;

int get_next_isoform_id()
{
#if ENABLE_THREADS
	gene_id_lock.lock();
#endif
	int next = next_isoform_id++;
#if ENABLE_THREADS
	gene_id_lock.unlock();
#endif
	return next;
}

void Isoform::get_gtf(vector<string>& gff_recs, 
					  const RefSequenceTable& rt,
                      set<AugmentedCuffOp>* hit_introns) const
{
	const char* ref_name = rt.get_name(_scaffold.ref_id());
	
    assert (ref_name != NULL);
    
	const char* strand_str = NULL;
	if (_scaffold.strand() == CUFF_STRAND_UNKNOWN)
		strand_str = ".";
	else if (_scaffold.strand() == CUFF_FWD)
		strand_str = "+";
	else
		strand_str = "-";
	
	int score = (int)(_FMI * 1000);
	score = min(1000, score);
	if (score == 0)
		score = 1;
	
	char buf[2048];	

	if (hit_introns != NULL)
    {
        sprintf(buf, 
			"%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%10.10lf\"; frac \"%lf\"; conf_lo \"%lf\"; conf_hi \"%lf\"; cov \"%lf\"; full_read_support \"%s\";\n",
			ref_name,
			_scaffold.left() + 1,
			_scaffold.right(), // GTF intervals are inclusive on both ends, but ours are half-open
			score,
			strand_str,
			gene_id().c_str(),
			trans_id().c_str(),
			_FPKM,
			_fraction,
			_confidence.low, 
			_confidence.high,
			_coverage,
            (_scaffold.has_struct_support(*hit_introns)) ? "yes":"no");
    }
    else
    {
        sprintf(buf, 
                "%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%10.10lf\"; frac \"%lf\"; conf_lo \"%lf\"; conf_hi \"%lf\"; cov \"%lf\";\n",
                ref_name,
                _scaffold.left() + 1,
                _scaffold.right(), // GTF intervals are inclusive on both ends, but ours are half-open
                score,
                strand_str,
                gene_id().c_str(),
                trans_id().c_str(),
                _FPKM,
                _fraction,
                _confidence.low, 
                _confidence.high,
                _coverage);
    }
    
    
	gff_recs.push_back(buf);
	
	int exon_num = 1;
	for (size_t op_id = 0; op_id < _scaffold.augmented_ops().size(); ++op_id)
	{
		const AugmentedCuffOp& op = _scaffold.augmented_ops()[op_id];
		if (op.opcode == CUFF_MATCH || op.opcode == CUFF_UNKNOWN)
		{
			const char* type = op.opcode == CUFF_MATCH ? "exon" : "missing_data";
			
			sprintf(buf, 
					"%s\tCufflinks\t\%s\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%10.10lf\"; frac \"%lf\"; conf_lo \"%lf\"; conf_hi \"%lf\"; cov \"%lf\";\n",
					ref_name,
					type,
					op.g_left() + 1,
					op.g_right(), // GTF intervals are inclusive on both ends, but ours are half-open
					score,
					strand_str,
					gene_id().c_str(),
					trans_id().c_str(),
					exon_num,
					_FPKM,
					_fraction,
					_confidence.low, 
					_confidence.high,
					_coverage);
			gff_recs.push_back(buf);
			
			exon_num++;
		}
		//gff_recs.push_back(buf);
	}
	
}


int next_gene_id = 1;

int get_next_gene_id()
{
#if ENABLE_THREADS
	gene_id_lock.lock();
#endif
	int next = next_gene_id++;
#if ENABLE_THREADS
	gene_id_lock.unlock();
#endif
	return next;
}


int next_skipped_region_id = 1;

int get_next_skipped_region_id()
{
#if ENABLE_THREADS
	gene_id_lock.lock();
#endif
	int next = next_skipped_region_id++;
#if ENABLE_THREADS
	gene_id_lock.unlock();
#endif
	return next;
}
