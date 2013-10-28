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

//allele
void Isoform::get_allele_gtf(vector<string>& gff_recs, 
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
	
	int paternal_score = (int)(_paternal_FMI * 1000);
	int maternal_score = (int)(_maternal_FMI * 1000);
	paternal_score = min(1000, paternal_score);
	maternal_score = min(1000, maternal_score);
	if (paternal_score == 0)
		paternal_score = 1;
	if (maternal_score == 0)
		maternal_score = 1;
	
	char buf[2048];	

	if (hit_introns != NULL)
    {
        sprintf(buf,
				"%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; allele_informative \"%d\"; paternal_FPKM \"%10.10lf\"; maternal_FPKM \"%10.10lf\"; paternal_frac \"%lf\"; maternal_frac \"%lf\"; paternal_conf_lo \"%lf\"; paternal_conf_hi \"%lf\"; maternal_conf_lo \"%lf\"; maternal_conf_hi \"%lf\"; paternal_cov \"%lf\"; maternal_cov \"%lf\"; full_read_support \"%s\";\n",
			ref_name,
			_scaffold.left() + 1,
			_scaffold.right(), // GTF intervals are inclusive on both ends, but ours are half-open
			paternal_score,
			maternal_score,
			strand_str,
			gene_id().c_str(),
			trans_id().c_str(),
			_allele_informative,	
			_paternal_FPKM,
            _maternal_FPKM,
			_paternal_fraction,
            _maternal_fraction,
			_paternal_confidence.low, 
			_paternal_confidence.high,
            _maternal_confidence.low, 
			_maternal_confidence.high,
			_paternal_coverage,
            _maternal_coverage,
            (_scaffold.has_struct_support(*hit_introns)) ? "yes":"no");
    }
    else
    {
        sprintf(buf, 
                "%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; allele_informative \"%d\"; paternal_FPKM \"%10.10lf\"; maternal_FPKM \"%10.10lf\"; paternal_frac \"%lf\"; maternal_frac \"%lf\"; paternal_conf_lo \"%lf\"; paternal_conf_hi \"%lf\"; maternal_conf_lo \"%lf\"; maternal_conf_hi \"%lf\"; paternal_cov \"%lf\"; maternal_cov \"%lf\";\n",
                ref_name,
                _scaffold.left() + 1,
                _scaffold.right(), // GTF intervals are inclusive on both ends, but ours are half-open
                paternal_score,
				maternal_score,
                strand_str,
                gene_id().c_str(),
                trans_id().c_str(),
				_allele_informative,
                _paternal_FPKM,
				_maternal_FPKM,
                _paternal_fraction,
				_maternal_fraction,
                _paternal_confidence.low, 
                _paternal_confidence.high,
				_maternal_confidence.low, 
                _maternal_confidence.high,
                _paternal_coverage,
				_maternal_coverage);
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
					"%s\tCufflinks\t\%s\t%d\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; allele_informative \"%d\"; exon_number \"%d\"; paternal_FPKM \"%10.10lf\"; maternal_FPKM \"%10.10lf\"; paternal_frac \"%lf\"; maternal_frac \"%lf\"; paternal_conf_lo \"%lf\"; paternal_conf_hi \"%lf\"; maternal_conf_lo \"%lf\"; maternal_conf_hi \"%lf\"; paternal_cov \"%lf\"; maternal_cov \"%lf\";\n",
					ref_name,
					type,
					op.g_left() + 1,
					op.g_right(), // GTF intervals are inclusive on both ends, but ours are half-open
					paternal_score,
					maternal_score,
					strand_str,
					gene_id().c_str(),
					trans_id().c_str(),
					_allele_informative,
					exon_num,
					_paternal_FPKM,
					_maternal_FPKM,
					_paternal_fraction,
					_maternal_fraction,
					_paternal_confidence.low, 
					_paternal_confidence.high,
					_maternal_confidence.low, 
					_maternal_confidence.high,
					_paternal_coverage,
					_maternal_coverage);
			
			gff_recs.push_back(buf);
			
			exon_num++;
			
		}
		//gff_recs.push_back(buf);
	}
}

void Gene::set_allele_informative_isoforms(const int n)
{
	for (size_t i = 0; i < _isoforms.size(); ++i)
	{
		int informative = 0;
		vector<const MateHit*> mates = _isoforms[i].scaffold().mate_hits();
		for (size_t j = 0; j < mates.size(); ++j)
		{
			if(mates[j]->allele() == ALLELE_PATERNAL || mates[j]->allele() == ALLELE_MATERNAL)
			{
				informative += 1;
				if(informative >= n) break;
			}
		}
		if(informative >= n)
		{
			_isoforms[i].set_allele_informativeness(true);
		}
		else{
			_isoforms[i].set_allele_informativeness(false);
		}
	}
}

bool Gene::has_allele_informative_isoforms() const
{
	bool result = false;
	for (size_t i = 0; i < _isoforms.size(); ++i)
	{
		if(_isoforms[i].get_allele_informativeness())
		{
			result = true;
			break;
		}
	}
	return result;
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
