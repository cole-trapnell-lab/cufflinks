/*
 *  bundles.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/6/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <list>
#include <map>
#include <numeric>
#include <boost/math/distributions/binomial.hpp>

#include "common.h"
#include "bundles.h"
#include "scaffolds.h"
#include "gff.h"
#include "GFastaFile.h"
#include "GStr.h"

using namespace std;
using boost::math::binomial;

struct ScaffoldSorter
{
	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
	bool operator()(const shared_ptr<Scaffold>& lhs, const shared_ptr<Scaffold>& rhs)
	{
		const char* lhs_name = rt.get_name(lhs->ref_id());
		const char* rhs_name = rt.get_name(rhs->ref_id());
		int c = strcmp(lhs_name, rhs_name);
		if (c != 0)
		{
			return c < 0;
		}
		if (lhs->left() != rhs->left())
		{
			return lhs->left() < rhs->left();
		}
        return false;
	}
	
	RefSequenceTable& rt;
};

char* getGSeqName(int gseq_id)
{
	return GffObj::names->gseqs.getName(gseq_id);
}
char* getFastaFile(int gseq_id) 
{
	if (fasta_dir == "") return NULL;
	GStr s(fasta_dir.c_str());
	s.trimR('/');
	s.appendfmt("/%s",getGSeqName(gseq_id));
	GStr sbase=s;
	if (!fileExists(s.chars())) 
		s.append(".fa");
	if (!fileExists(s.chars())) s.append("sta");
	if (fileExists(s.chars())) return Gstrdup(s.chars());
	else
	{
		GMessage("Warning: cannot find genomic sequence file %s{.fa,.fasta}\n",sbase.chars());
		return NULL;
	}
}

//FIXME: needs refactoring
void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs,
				   bool loadFPKM) 
{
	if (loadSeqs)
		fprintf(stderr,"Loading reference annotation and sequence\n");
	else
		fprintf(stderr,"Loading reference annotation\n");
    
	GList<GSeqData> ref_rnas;
	
	if (ref_mRNA_file)
	{
		read_transcripts(ref_mRNA_file, ref_rnas);
	}
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			char* name = GffObj::names->gseqs.getName(ref_rnas[j]->gseq_id);
			
			int f = 0;
			int r = 0;
			int u  = 0;
			GffObj* rna_p;
			RefID ref_id = rt.get_id(name, NULL);
			int f_count = ref_rnas[j]->mrnas_f.Count();
			int r_count = ref_rnas[j]->mrnas_r.Count();
			int u_count = ref_rnas[j]->umrnas.Count();
			
			
			while(!(f==f_count && r==r_count && u==u_count))
			{	
				CuffStrand strand;
				
				if (f < f_count)
                {
					rna_p = ref_rnas[j]->mrnas_f[f++];
					strand = CUFF_FWD;
                }
				else if (r < r_count) 
				{
					rna_p = ref_rnas[j]->mrnas_r[r++];
					strand = CUFF_REV;
				}
				else 
				{
					rna_p = ref_rnas[j]->umrnas[u++];
					strand = CUFF_STRAND_UNKNOWN;
				}
				
				
				GffObj& rna = *rna_p;
				
                
				if (loadSeqs && rna.gseq_id != last_gseq_id) //next chromosome
				{
					delete faseq;
					faseq = NULL;
					last_gseq_id = rna.gseq_id;
					char* sfile = getFastaFile(last_gseq_id);
                    if (sfile != NULL) 
                    {
						if (verbose)
							GMessage("Processing sequence from fasta file '%s'\n",sfile);
						faseq = new GFaSeqGet(sfile, false);
						faseq->loadall();
						GFREE(sfile);
                    }
                    else 
                    {
						assert (false);
                    }
                }
                
				vector<AugmentedCuffOp> ops;
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
					ops.push_back(AugmentedCuffOp(CUFF_MATCH, ex.start - 1, ex.end - ex.start + 1));
					
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
						ops.push_back(AugmentedCuffOp(CUFF_INTRON, ex.end, next_ex.start - ex.end - 1));
					}
				}
				
				Scaffold ref_scaff(ref_id, strand, ops);
				
				char* rna_seq;
				int seqlen=0;
				if (loadSeqs){ 
					assert (faseq);
					rna_seq = rna.getSpliced(faseq, false, &seqlen);
				}
                
				if (rna.getID())
					ref_scaff.annotated_trans_id(rna.getID());
				
				
				if (rna.getGene())
					ref_scaff.annotated_gene_id(rna.getGene());
				
				char* short_name = rna.getAttr("gene_name");
				if (short_name)
					ref_scaff.annotated_gene_name(short_name);
				
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
					ref_scaff.annotated_protein_id(protein_id);
				
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
					ref_scaff.annotated_tss_id(tss_id);
				
				
				if (loadFPKM)
				{
					const char* expr = rna.getAttr("FPKM");
					if (expr!=NULL) {
						if (expr[0]=='"') expr++;
						ref_scaff.fpkm(strtod(expr, NULL));
					}
				}
				
				if (loadSeqs)
                {
					string rs = rna_seq;
					std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
					ref_scaff.seq(rs);
					GFREE(rna_seq);
				}
				
				ref_mRNAs.push_back(shared_ptr<Scaffold>(new Scaffold(ref_scaff))); 
			}
		}
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
	}
	delete faseq;
}


int HitBundle::_next_id = 0;

bool HitBundle::add_hit(const MateHit& hit)
{
	if (_final)
		return false;
	
	// Update the bounds on the span
	if (hit.left() < _leftmost)
		_leftmost = hit.left();
	if (hit.right() > _rightmost)
		_rightmost = hit.right();
	
	_hits.push_back(hit);
	return true;
}

struct HitLessScaffold
{
	bool operator()(shared_ptr<Scaffold> x)
	{
		return x->mate_hits().empty();
	}
};

void HitBundle::add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                             shared_ptr<ReadHit const> bh)
{
	if (bh->partner_ref_id() == 0 
		|| bh->partner_ref_id() != bh->ref_id() ||
		abs((int)bh->partner_pos() - (int)bh->left()) > max_partner_dist)
	{
		// This is a singleton, so just make a closed MateHit and
		// continue
		MateHit m(rg_props, bh->ref_id(), bh, shared_ptr<ReadHit const>());
		add_hit(m);
	}
	else
	{
		OpenMates::iterator mi = _open_mates.find(bh->left());
		
		// Does this read hit close an open mate?
		if (mi == _open_mates.end())
		{
			// No, so add it to the list of open mates, unless we would
			// already have seen it's partner
			if(bh->left() < bh->partner_pos())
			{
				MateHit open_hit(rg_props,
                                 bh->ref_id(), 
                                 bh, 
                                 shared_ptr<ReadHit const>());
				
				pair<OpenMates::iterator, bool> ret;
				ret = _open_mates.insert(make_pair(bh->partner_pos(), 
												  list<MateHit>()));
				
				ret.first->second.push_back(open_hit);
			}
			else
			{
				add_hit(MateHit(rg_props,bh->ref_id(), bh, shared_ptr<ReadHit const>()));
			}
		}
		else
		{
			
			bool found_partner = false;
			// Maybe, see if we can find an ID match in the list of
			// open mates expecting a partner at this position
			for (list<MateHit>::iterator pi = mi->second.begin();
				 pi != mi->second.end();
				 ++pi)
			{
				MateHit& pm = *pi;
				
				if (pm.insert_id() == bh->insert_id())
				{
					// Found a partner?
					
					Scaffold L(MateHit(rg_props, bh->ref_id(), pm.left_alignment(), shared_ptr<ReadHit const>()));
					Scaffold R(MateHit(rg_props, bh->ref_id(), bh, shared_ptr<ReadHit const>()));
					
					bool strand_agree = L.strand() == CUFF_STRAND_UNKNOWN ||
					R.strand() == CUFF_STRAND_UNKNOWN ||
					L.strand() == R.strand();
					
					//bool orientation_agree = pm.left_alignment()->antisense_align() != bh->antisense_align();
					
					if (strand_agree && 
                        (!Scaffold::overlap_in_genome(L, R, olap_radius) ||
                         Scaffold::compatible(L,R)))
					{					
						pm.right_alignment(bh);
						add_hit(pm);
						mi->second.erase(pi);
						if (mi->second.empty())
							_open_mates.erase(mi);
						
						found_partner = true;
						break;
					}
				}
			}
			
			if (!found_partner)
			{
				// If we got here, couldn't actually close any mates with
				// this read hit, so open a new one, unless we can never
				// close this one
				if(bh->left() < bh->partner_pos())
				{
					MateHit open_hit(rg_props, bh->ref_id(), bh, shared_ptr<ReadHit const>());
					
					pair<OpenMates::iterator, bool> ret;
					ret = _open_mates.insert(make_pair(bh->partner_pos(), 
													  list<MateHit>()));
					
					ret.first->second.push_back(open_hit);
				}
				else
				{
					add_hit(MateHit(rg_props, bh->ref_id(), bh, shared_ptr<ReadHit const>()));
				}
			}
		}
	}
}

void HitBundle::finalize_open_mates()
{
	for (OpenMates::iterator itr = _open_mates.begin(); 
		 itr != _open_mates.end(); 
		 ++itr)
	{
		for (list<MateHit>::iterator mi = itr->second.begin(); mi != itr->second.end(); ++mi)
		{
			add_hit(*mi);
		}
        itr->second.clear();
	}
    _open_mates.clear();
}

void HitBundle::remove_hitless_scaffolds()
{
	vector<shared_ptr<Scaffold> >::iterator new_end = remove_if(_ref_scaffs.begin(),
												   _ref_scaffs.end(),
												   HitLessScaffold());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());	
}

void HitBundle::finalize(bool is_combined)
{
	_final = true;
    
	if (!is_combined)
	{
		sort(_hits.begin(), _hits.end(), mate_hit_lt);
		collapse_hits(_hits, _non_redundant);
		sort(_ref_scaffs.begin(), _ref_scaffs.end(), scaff_lt_rt_oplt_sp);
		vector<shared_ptr<Scaffold> >::iterator new_end = unique(_ref_scaffs.begin(), 
												_ref_scaffs.end(),
												StructurallyEqualScaffolds());
		_ref_scaffs.erase(new_end, _ref_scaffs.end());
        vector<shared_ptr<Scaffold> >(_ref_scaffs).swap(_ref_scaffs);
	}
	
    for (size_t j = 0; j < _ref_scaffs.size(); ++j)
	{
		_ref_scaffs[j]->clear_hits();
	}
    
    finalize_open_mates();
	
	for (size_t i = 0; i < _hits.size(); ++i)
	{
		const MateHit* hit = &(_hits[i]);
		
		Scaffold hs(*hit);
		
        if (i >= 1)
        {
            assert (hit->ref_id() == _hits[i-1].ref_id());
        }
        
		for (size_t j = 0; j < _ref_scaffs.size(); ++j)
		{
			// add hit only adds if the hit is structurally compatible
			if (_ref_scaffs[j]->Scaffold::contains(hs))
			{
				_ref_scaffs[j]->add_hit(hit);
			}
		}
	}
		
	for (size_t i = 0; i < _ref_scaffs.size(); ++i)
	{
		if (_ref_scaffs[i]->left() < _leftmost)
			_leftmost = _ref_scaffs[i]->left();
		if (_ref_scaffs[i]->right() > _rightmost)
			_rightmost = _ref_scaffs[i]->right();
	}
}

void print_sort_error(const char* last_chr_name, 
                      int last_chr_pos, 
                      const char* bh_name, 
                      int bh_pos)
{
    fprintf(stderr, "Error: this SAM file doesn't appear to be correctly sorted!\n");
    fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
            bh_name,
            bh_pos,
            last_chr_name,
            last_chr_pos);
    fprintf(stderr, "You may be able to fix this by running:\n\t$ LC_ALL=\"C\" sort -k 3,3 -k 4,4n input.sam > fixed.sam\n");
}

shared_ptr<ReadHit const> BundleFactory::next_valid_alignment()
{
    const char* hit_buf;
	size_t hit_buf_size = 0;
    shared_ptr<ReadHit const> bh;
    
    while (true)
    {
    
        if (!_hit_fac->next_record(hit_buf, hit_buf_size))
            break;
        
        ReadHit tmp;
        if (!_hit_fac->get_hit_from_buf(hit_buf, tmp, false))
            continue;
        
        if (tmp.error_prob() > max_phred_err_prob)
            continue;
        
        if (tmp.ref_id() == 84696373) // corresponds to SAM "*" under FNV hash. unaligned read record  
            continue;
        
        if (spans_bad_intron(tmp))
            continue;
        
        bool hit_within_mask = false;
        
        // We want to skip stuff that overlaps masked GTF records, so 
        // sync up the masking chromosome
        if (!mask_gtf_recs.empty() && 
            next_mask_scaff != mask_gtf_recs.end() &&
            (*next_mask_scaff)->ref_id() != tmp.ref_id())
        {
            bool found_scaff = false;
            vector<shared_ptr<Scaffold> >::iterator curr_mask_scaff = mask_gtf_recs.begin();
            for (size_t i = 0; i < _mask_scaff_offsets.size(); ++i)
            {
                if (_mask_scaff_offsets[i].first == tmp.ref_id())
                {
                    curr_mask_scaff = _mask_scaff_offsets[i].second;
                    found_scaff = true;
                    break;
                }
            }
            
            next_mask_scaff = curr_mask_scaff;
        }
        
        //check that we aren't sitting in the middle of a masked scaffold
        while (next_mask_scaff != mask_gtf_recs.end() && 
               (*next_mask_scaff)->ref_id() == tmp.ref_id() &&
               (*next_mask_scaff)->right() <= tmp.left())
        {
            if ((*next_mask_scaff)->left() >= tmp.left())
            {
                break;
            }
            
            next_mask_scaff++;
        }
        
        if (next_mask_scaff != mask_gtf_recs.end() &&
            (*next_mask_scaff)->ref_id() == tmp.ref_id() &&
            (*next_mask_scaff)->left() <= tmp.left() &&
            (*next_mask_scaff)->right() >= tmp.right())
        {
            hit_within_mask = true;
        }
        
        if (hit_within_mask)
            continue;
        bh = shared_ptr<ReadHit const>(new ReadHit(tmp));
        break;
    }
    return bh;
}

// This is my least favorite function in Cufflinks.  It should be banished to
// Hell and re-written.  It is utterly loathesome.
bool BundleFactory::next_bundle(HitBundle& bundle)
{
	//char bwt_buf[2048];
	
    // The most recent RefID and position we've seen in the hit stream
	RefID last_hit_ref_id_seen = 0;
    int last_hit_pos_seen = 0;
    
    // The most recent RefID and position we've seen in the reference RNA stream
	RefID last_ref_rna_ref_id_seen = 0;
    int last_ref_rna_pos_seen = 0;
    
    // The first RefID we saw in either stream
	RefID first_ref_id_seen = 0;

	
	int right_bundle_boundary = 0;
	int left_bundle_boundary = -1;
	
    
	while(true)
	{
        shared_ptr<ReadHit const> bh = next_valid_alignment();
        

		// Initialize the bundle boundaries using the hit or the next 
        // reference transcript, whichever is first in the stream
		if (left_bundle_boundary == -1)
        {
            if (_ref_driven && next_ref_scaff != ref_mRNAs.end() && bh != NULL) // We're not done
            {
                if (bh->ref_id() != (*next_ref_scaff)->ref_id())
                {
                    const char* bh_chr_name = _hit_fac->ref_table().get_name(bh->ref_id());
                    const char* ref_rna_chr_name = _hit_fac->ref_table().get_name((*next_ref_scaff)->ref_id());
                    const char* last_bh_chr_name = _hit_fac->ref_table().get_name(last_hit_ref_id_seen);
                    
                    // if the BAM/SAM file isn't lexicographically sorted by chromosome, print an error
                    // and exit.
                    if (last_bh_chr_name && strcmp(last_bh_chr_name, bh_chr_name) > 0)
                    { 
                        print_sort_error(last_bh_chr_name, 
                                         last_hit_pos_seen, 
                                         bh_chr_name, 
                                         bh->left());
                        exit(1);
                    }
                    
                    if (bh_chr_name && ref_rna_chr_name)
                    {
                        if (strcmp(bh_chr_name, ref_rna_chr_name) < 0)
                        {
                            // hit is lexicographically less than the reference,
                            // so skip the hit.
                            continue;
                        }
                        else 
                        {
                            // reference is lexicographically less than
                            // the hit, so wait on the hit and emit this bundle.
                            left_bundle_boundary = (*next_ref_scaff)->left();
                            right_bundle_boundary = (*next_ref_scaff)->right();
                            first_ref_id_seen = (*next_ref_scaff)->ref_id();
                            
                        }

                    }
                    // Then these hits are on a chromosome not in the reference annotation, but which are lexicographically
                    // less than the next reference record, so just skip them
                    
                    continue;
                }
				
                if ((bh->ref_id() == (*next_ref_scaff)->ref_id() &&
                     bh->right() <= (*next_ref_scaff)->left()))
                {
                    // if this hit doesn't overlap the next ref transcript
                    // and falls to its left, we should skip it.  We need
                    // this to maintain synchronization of multiple
                    // streams driven by the same set of reference transcripts
                    
                    // if the hit stream disagrees with the reference transcript stream
                    // on which RefID is next, just skip the hit.
                    continue;
//                    left_bundle_boundary = bh->left();
//                    right_bundle_boundary = bh->left();
//                    first_ref_id_seen = bh->ref_id();
                }
                else
                {
                    left_bundle_boundary = (*next_ref_scaff)->left();
                    right_bundle_boundary = (*next_ref_scaff)->right();
                    first_ref_id_seen = (*next_ref_scaff)->ref_id();
                }

            }
            else if (_ref_driven)
            {
                if (next_ref_scaff != ref_mRNAs.end())
                {
                    left_bundle_boundary = (*next_ref_scaff)->left();
                    right_bundle_boundary = (*next_ref_scaff)->right();
                    first_ref_id_seen = (*next_ref_scaff)->ref_id();
                }
                else 
                {
                    return false;
                }

            }
            else if (bh != NULL)
            {
                left_bundle_boundary = bh->left();
                right_bundle_boundary = bh->right();
                first_ref_id_seen = bh->ref_id();
            }
            else
            {
                if (!_hit_fac->records_remain() && next_ref_scaff == ref_mRNAs.end())
                {
                    return false;
                }
                else 
                {
                    // if we get here, there are more records, but this one is
                    // null, and we've reached the end of the reference trancripts
                    // or there weren't any
                    continue;
                }
            }
        }

        assert(left_bundle_boundary != -1);

        // if the hit here overlaps the current bundle interval,
        // we have to include it, and expand the bundle interval
        if (bh && !_ref_driven
            && overlap_in_genome(bh->left(),bh->right(),left_bundle_boundary, right_bundle_boundary + olap_radius))
        {
            right_bundle_boundary = max(right_bundle_boundary, bh->right());
        }
        
        // expand the bundle interval as far as needed to include the overlapping
        // chain of reference transcripts that also overlap the initial bundle
        // interval, whether the initial bundle interval came from a hit or the
        // current reference transcript in the GTF record stream.
		bool hit_within_boundary = false;
		bool olaps_reference = false;
		if (_ref_driven)
		{
            while(next_ref_scaff < ref_mRNAs.end())
            {
                if (first_ref_id_seen == (*next_ref_scaff)->ref_id())
                {
                    // if the ref_scaff here overlaps the current bundle interval,
                    // we have to include it, and expand the bundle interval
                    if (overlap_in_genome((*next_ref_scaff)->left(),(*next_ref_scaff)->right(),left_bundle_boundary, right_bundle_boundary))
                    {
                        right_bundle_boundary = max(right_bundle_boundary, (*next_ref_scaff)->right());
                        olaps_reference = true;
                    }
                    
                    // If the ref transcript is to the right of the bundle interval, 
                    // then we can leave this the ref transcript out for now
                    if ((*next_ref_scaff)->left() >= right_bundle_boundary)
                    {
                        // we've gone beyond the bundle and this hit, so we're 
                        // not going to expand the right boundary any further
                        // until we see more hits
                        break;
                    }
                }
                else 
                {
                    break;
                }

                
                last_ref_rna_ref_id_seen = (*next_ref_scaff)->ref_id();
                last_ref_rna_pos_seen = (*next_ref_scaff)->left();
                
                next_ref_scaff++;
            }
		}

        if (bh == NULL)
        {
            break;
        }
        
        bool past_right_end = false;
        if (last_hit_ref_id_seen == 0 || bh->ref_id() == first_ref_id_seen)
        {
            if (bh->left() >= left_bundle_boundary && bh->right() <= right_bundle_boundary)
            {
                hit_within_boundary = true;
            }
            else if (bh->left() >= right_bundle_boundary)
            {
                past_right_end = true;
            }
        }
        else
        {
            past_right_end = true;
            const char* bh_chr_name = _hit_fac->ref_table().get_name(bh->ref_id());
            const char* last_chr_name = _hit_fac->ref_table().get_name(last_hit_ref_id_seen);
            
            if (strcmp(last_chr_name, bh_chr_name) > 0)
            { 
                print_sort_error(last_chr_name, 
                                 last_hit_pos_seen, 
                                 bh_chr_name, 
                                 bh->left());
                exit(1);
            }
        }
        
		if (hit_within_boundary)
		{
			if (bh->left() < last_hit_pos_seen)
			{
                const char* bh_chr_name = _hit_fac->ref_table().get_name(bh->ref_id());
                const char* last_chr_name = _hit_fac->ref_table().get_name(last_hit_ref_id_seen);
				print_sort_error(last_chr_name, 
                                 last_hit_pos_seen,
                                 bh_chr_name,
                                 bh->left());
				exit(1);
			}
			
            // We want to drive the bundling entirely by the ref_mRNAs, if they
            // are there.
            if (!_ref_driven)
            {
                if (bh->right() + olap_radius - left_bundle_boundary < (int)max_gene_length)
                {
                    right_bundle_boundary = max(right_bundle_boundary, bh->right() + olap_radius);
                }
                
                bool singleton = (bh->partner_ref_id() == 0 
                                  || bh->partner_ref_id() != bh->ref_id() ||
                                  abs(bh->partner_pos() - bh->left()) > max_partner_dist);
                if (!singleton)
                {
                    if ((int)bh->partner_pos() + olap_radius - (int)bh->left() < (int)max_partner_dist)
                    {
                        right_bundle_boundary = max(right_bundle_boundary, bh->partner_pos() + olap_radius);
                    }
                }
            }
			
			bundle.add_open_hit(read_group_properties(), bh);
		}
		else if (past_right_end)
		{
			_hit_fac->undo_hit();
			break;
		}
		
		last_hit_ref_id_seen = bh->ref_id();
		last_hit_pos_seen = bh->left();
        
		//curr_pos = ftello(hit_file);
	}
	
    assert(left_bundle_boundary != -1);
    
	bundle.finalize_open_mates();
	
	if (_ref_driven)
	{
		vector<shared_ptr<Scaffold> >::iterator itr = next_ref_scaff;
		while(itr < ref_mRNAs.end())
		{
			if ((*itr)->ref_id() != last_ref_rna_ref_id_seen || (*itr)->left() >= right_bundle_boundary)
				break;
			itr++;
		}
		
		bool past_first_id = false;
		while(itr >= ref_mRNAs.begin())
		{
			if (itr < ref_mRNAs.end())
			{
				// if we haven't backed up to the scaffold we need to be on
				// keep scanning through the reference annotations
				if ((*itr)->ref_id() == last_ref_rna_ref_id_seen)
				{
					// Now we're on the right scaffold, gotta get to the right
					// coords

					if (overlap_in_genome((*itr)->left(), 
                                          (*itr)->right(), 
                                          left_bundle_boundary, 
                                          right_bundle_boundary))

					{	
						bundle.add_ref_scaffold(*itr);
					}
					else if ((*itr)->right() < left_bundle_boundary)
					{	
						// This reference record is now to the left of 
						// the bundle of alignments
						break;
					}
				}
				else if (past_first_id)
				{
					break;
				}
				
				// If this reference record is on the same scaffold as we were 
				// on when we first entered this routine, going back any further
				// will be unproductive, so it's safe to break out of the loop
				// as long as the bundle is on a different scaffold.
				if ((*itr)->ref_id() == first_ref_id_seen)
				{
					past_first_id = true;
				}
			}
			--itr;
		}
	}
	
	assert (left_bundle_boundary != -1);
	bundle.finalize();
    assert(bundle.right() != -1);
	//bundle_out.remove_hitless_scaffolds();
    
    return true;
}

struct IntronSpanCounter
{
	IntronSpanCounter() : left_reads(0), little_reads(0), total_reads(0), multimap_reads(0), fwd_strand_frags(0) {}
	size_t left_reads;
	size_t little_reads; // small span overhang
	size_t total_reads;
	size_t multimap_reads;
    size_t fwd_strand_frags;
	vector<size_t> hist;
};

typedef map<AugmentedCuffOp, IntronSpanCounter> IntronCountTable;

void count_introns_in_read(const ReadHit& read,
						   IntronCountTable& intron_counts)
{
	const vector<CigarOp>& cig = read.cigar();
	
	int read_len = read.read_len();
	int small_anchor = (int)floor(read_len * small_anchor_fraction);
	
	int r_left = 0;
	int g_left = read.left();
	
	for (size_t i = 0; i < cig.size(); ++i)
	{
		assert(cig[i].length >= 0);
		switch(cig[i].opcode)
		{
			case MATCH:
				//ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, cig[i].length));
				g_left += cig[i].length;
				r_left += cig[i].length;
				break;
				
			case REF_SKIP:
			{	
				AugmentedCuffOp intron(CUFF_INTRON, g_left, cig[i].length);
				pair<IntronCountTable::iterator, bool> ins_itr;
				ins_itr = intron_counts.insert(make_pair(intron, IntronSpanCounter()));
				IntronCountTable::iterator itr = ins_itr.first;
				itr->second.total_reads++;
				
				if (read.error_prob() > 0.9)
				{
					itr->second.multimap_reads++;
				}
				
				if ( r_left <= small_anchor || (read_len - r_left) < small_anchor)
				{
					itr->second.little_reads++;
				}
                
                if (read.source_strand() == CUFF_FWD)
                {
                    itr->second.fwd_strand_frags;
                }
                else 
                {
                    assert(read.source_strand() == CUFF_REV);
                }

				
				vector<size_t>& hist = itr->second.hist;
				if (hist.size() < (size_t)read_len)
				{
					size_t num_new_bins = read_len - hist.size();
					size_t new_left_bins = (size_t)floor(num_new_bins / 2.0);
					size_t new_right_bins = (size_t)ceil(num_new_bins / 2.0);
					hist.insert(hist.begin(), new_left_bins, 0);
					hist.insert(hist.end(), new_right_bins, 0);
				}
				
				assert (r_left < hist.size());
				hist[r_left]++;
				//ops.push_back(AugmentedCuffOp(CUFF_INTRON, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
			}
				
			case SOFT_CLIP:
				g_left += cig[i].length;
				break;
            case HARD_CLIP:
				break;
			default:
				assert(false);
				break;
		}
	}
}

void minor_introns(int bundle_length,
				   int bundle_left,
				   const IntronCountTable& intron_counts,
				   vector<AugmentedCuffOp>& bad_introns,
				   double fraction)

{
	for(IntronCountTable::const_iterator itr = intron_counts.begin();
		itr != intron_counts.end(); 
		++itr)
	{
		pair<AugmentedCuffOp, IntronSpanCounter> itr_cnt_pair = *itr;
		const IntronSpanCounter itr_spans = itr_cnt_pair.second;
		
		double doc = itr_spans.total_reads;
		
		for (IntronCountTable::const_iterator itr2 = intron_counts.begin();
			 itr2 != intron_counts.end(); 
			 ++itr2)
		{	
			if (itr == itr2 ||
				!AugmentedCuffOp::overlap_in_genome(itr->first, itr2->first))
			{
				continue;
			}
			
			pair<AugmentedCuffOp, IntronSpanCounter> itr2_cnt_pair = *itr2;
			const IntronSpanCounter itr2_spans = itr2_cnt_pair.second;
			
			double thresh = itr2_spans.total_reads * fraction;
			if (doc < thresh)
			{
				//#if ASM_VERBOSE
				//							fprintf(stderr, "\t Filtering intron (due to overlap) %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
				//#endif	
				bool exists = binary_search(bad_introns.begin(), 
											bad_introns.end(), 
											itr->first);
				if (!exists)
				{
#if ASM_VERBOSE
					fprintf(stderr, "Filtering intron %d-%d spanned by %lu reads based on overlap with much more abundant intron: %d-%d spanned by %lu reads\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads,
							itr2->first.g_left(), 
							itr2->first.g_right(), 
							itr2->second.total_reads);
					
#endif
					bad_introns.push_back(itr->first);
					sort(bad_introns.begin(), bad_introns.end());
				}
			}
            
//            if ((itr->second.fwd_strand_frags == 0 &&
//                 itr2->second.fwd_strand_frags != 0) ||
//                (itr2->second.fwd_strand_frags == 0 &&
//                 itr->second.fwd_strand_frags != 0))
//            {
//                int itr1_L = itr->first.g_left();
//                int itr1_R = itr->first.g_right();
//                int itr2_L = itr2->first.g_left();
//                int itr2_R = itr2->first.g_right();
//                
//                if (abs(itr1_L - itr2_L) < 25 && abs(itr1_R - itr2_R) < 25)
//                {
//                    int a = 3;
//                }
//            }
		}
	}
}

void multimapping_introns(int bundle_length,
						  int bundle_left,
						  const IntronCountTable& intron_counts,
						  vector<AugmentedCuffOp>& bad_introns,
						  double fraction)

{
	for(IntronCountTable::const_iterator itr = intron_counts.begin();
		itr != intron_counts.end(); 
		++itr)
	{
		pair<AugmentedCuffOp, IntronSpanCounter> itr_cnt_pair = *itr;
		const IntronSpanCounter itr_spans = itr_cnt_pair.second;
		
		double doc = itr_spans.total_reads;
		double multi = itr_spans.multimap_reads;
		
		double multi_fraction = multi / doc;
		
		if (multi_fraction > fraction)
		{
			bool exists = binary_search(bad_introns.begin(), 
										bad_introns.end(), 
										itr->first);
			if (!exists)
			{
	#if ASM_VERBOSE
				fprintf(stderr, "Filtering intron %d-%d spanned by %lu reads because %lg percent are multireads.\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads,
						multi_fraction * 100);
				
	#endif
				bad_introns.push_back(itr->first);
				sort(bad_introns.begin(), bad_introns.end());
			}
		}
	}
}


void identify_bad_splices(const HitBundle& bundle, 
						  BadIntronTable& bad_splice_ops)
{
	// Tracks, for each intron, how many reads
	IntronCountTable intron_counts;
	
	RefID ref_id = bundle.ref_id();
	
	pair<BadIntronTable::iterator, bool> ins_itr;
	ins_itr = bad_splice_ops.insert(make_pair(ref_id, vector<AugmentedCuffOp>()));
	vector<AugmentedCuffOp>& bad_introns = ins_itr.first->second;
	
	foreach (const MateHit& hit, bundle.hits())
	{
		if (hit.left_alignment())
		{
			count_introns_in_read(*hit.left_alignment(), intron_counts);
		}
		if (hit.right_alignment())
		{
			count_introns_in_read(*hit.right_alignment(), intron_counts);
		}
	}
	
	minor_introns(bundle.length(), bundle.left(), intron_counts, bad_introns, min_isoform_fraction);
	multimapping_introns(bundle.length(), bundle.left(), intron_counts, bad_introns, 0.5);
	for (IntronCountTable::iterator itr = intron_counts.begin();
		 itr != intron_counts.end();
		 ++itr)
	{
		if (binary_search(bad_introns.begin(), 
						  bad_introns.end(), 
						  itr->first))
		{
			continue;
		}
		pair<AugmentedCuffOp, IntronSpanCounter> cnt_pair = *itr;
		try
		{
			const IntronSpanCounter spans = cnt_pair.second;
			
			//			binomial read_half_dist(spans.total_reads, success_fraction);
			//			double left_side_p = cdf(read_half_dist, spans.total_reads - spans.left_reads);
			//			double right_side_p = cdf(complement(read_half_dist, spans.left_reads));
			
			
			double success = 2 * small_anchor_fraction;
			
			binomial read_half_dist(spans.total_reads, success);
			double right_side_p;
			
			// right_side_p describes the chance that we'd observe at least 
			// this many small overhang reads by chance with an unbiased 
			// distribution over a normal (e.g. non-artifact) junction
			if (spans.little_reads > 0)
			{
				right_side_p = 1.0 - cdf(read_half_dist, spans.little_reads - 1);
			}
			else 
			{
				right_side_p = 1.0;
			}
			
			double left_side_p = 0;
#if ASM_VERBOSE
			double expected = success * spans.total_reads;
#endif
            //double excess = spans.little_reads - expected;
			
			// left_side_p describes the chance that we'd observe this few or
			// fewer small overhang reads by chance with an unbiased 
			// distribution over a normal (e.g. non-artifact) junction
			if (spans.little_reads > 0)
			{
				left_side_p = cdf(read_half_dist, spans.little_reads);
			}
			else 
			{
				left_side_p = cdf(read_half_dist, 0);
			}
			
			//double alpha = 0.05;
			//double right_side_p = 0;
			
			// Two-tailed binomial test:
//			if (left_side_p < (binomial_junc_filter_alpha / 2.0) || 
//				right_side_p < (binomial_junc_filter_alpha / 2.0))
			// One-tailed binomial test
			
			bool filtered = false;
			
			const IntronSpanCounter& counter = itr->second;
			
			if (right_side_p < (binomial_junc_filter_alpha))
			{
				double overhang_ratio = counter.little_reads / (double) counter.total_reads;
				if (counter.total_reads < 100 || overhang_ratio >= 0.50)
				{
#if ASM_VERBOSE
					fprintf(stderr, "Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
					
#endif
					filtered = true;
					
					bool exists = binary_search(bad_introns.begin(), 
												bad_introns.end(), 
												itr->first);
					if (!exists)
					{
						bad_introns.push_back(itr->first);
						sort(bad_introns.begin(), bad_introns.end());
					}
				}
			}
			
			vector<size_t> hist = itr->second.hist;
			if (itr->second.total_reads > 1000)
			{
				sort(hist.begin(), hist.end());
				size_t median = (size_t)floor(hist.size() / 2);
				if (median <= hist.size() && hist[median] == 0)
				{
#if ASM_VERBOSE
					fprintf(stderr, "Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
					
#endif
					filtered = true;
					
					bool exists = binary_search(bad_introns.begin(), 
												bad_introns.end(), 
												itr->first);
					if (!exists)
					{
						bad_introns.push_back(itr->first);
						sort(bad_introns.begin(), bad_introns.end());
					}
				}
			}
			
			if (!filtered)
			{
#if ASM_VERBOSE
				fprintf(stderr, "Accepting intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads, 
						itr->second.little_reads, 
						expected,
						left_side_p,
						right_side_p);
#endif
				
			}
		}
		
		
		catch(const std::exception& e)
		{
			//
			/*`
			 [#coinflip_eg_catch]
			 It is always essential to include try & catch blocks because
			 default policies are to throw exceptions on arguments that
			 are out of domain or cause errors like numeric-overflow.
			 
			 Lacking try & catch blocks, the program will abort, whereas the
			 message below from the thrown exception will give some helpful
			 clues as to the cause of the problem.
			 */
			std::cout <<
			"\n""Message from thrown exception was:\n   " << e.what() << std::endl;
		}
		
	}
}

bool BundleFactory::spans_bad_intron(const ReadHit& read)
{

	const vector<CigarOp>& cig = read.cigar();
	
	size_t g_left = read.left();
	BadIntronTable::const_iterator itr = _bad_introns.find(read.ref_id());
	if (itr == _bad_introns.end())
		return false;
	
	const vector<AugmentedCuffOp>& bi = itr->second; 
	for (size_t i = 0; i < cig.size(); ++i)
	{
		assert(cig[i].length >= 0);
		switch(cig[i].opcode)
		{
			case MATCH:
				//ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
				
			case REF_SKIP:
			{	
				AugmentedCuffOp intron(CUFF_INTRON, g_left, cig[i].length);
				if (binary_search(bi.begin(), bi.end(), intron))
				{
					return true;
				}
				
				//ops.push_back(AugmentedCuffOp(CUFF_INTRON, g_left, cig[i].length));
				g_left += cig[i].length;
				break;
			}
				
			case SOFT_CLIP:
				g_left += cig[i].length;
				break;
                
            case HARD_CLIP:
				break;
			default:
				assert(false);
				break;
		}
	}
	
	return false;
}
