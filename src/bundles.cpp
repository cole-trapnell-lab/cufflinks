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
#include <boost/crc.hpp>

#include "common.h"
#include "bundles.h"
#include "scaffolds.h"

#include "abundances.h"

using namespace std;
using boost::math::binomial;

//struct ScaffoldSorter
//{
//	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
//	bool operator()(boost::shared_ptr<Scaffold const> lhs, boost::shared_ptr<Scaffold const> rhs)
//	{
//        assert (lhs);
//        assert (rhs);
//		const char* lhs_name = rt.get_name(lhs->ref_id());
//		const char* rhs_name = rt.get_name(rhs->ref_id());
//		int c = strcmp(lhs_name, rhs_name);
//		if (c != 0)
//		{
//			return c < 0;
//		}
//		if (lhs->left() != rhs->left())
//		{
//			return lhs->left() < rhs->left();
//		}
//        return false;
//	}
//	
//	RefSequenceTable& rt;
//};

struct ScaffoldSorter
{
	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
	bool operator()(boost::shared_ptr<Scaffold const> lhs, boost::shared_ptr<Scaffold const> rhs)
	{
        //assert (lhs);
        //assert (rhs);
        if (!lhs || !rhs)
            return false;
		int lhs_order = rt.observation_order(lhs->ref_id());
        assert (lhs_order != -1);
		int rhs_order = rt.observation_order(rhs->ref_id());
        assert (rhs_order != -1);

		if (lhs_order != rhs_order)
		{
			return lhs_order < rhs_order;
		}
		if (lhs->left() != rhs->left())
		{
			return lhs->left() < rhs->left();
		}
        return false;
	}
	
	RefSequenceTable& rt;
};

//FIXME: needs refactoring
void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<boost::shared_ptr<Scaffold> >& ref_mRNAs,
                   boost::crc_32_type& gtf_crc_result,
				   bool loadSeqs,
				   bool loadFPKM) 
{
	if (loadSeqs)
		ProgressBar p_bar("Loading reference annotation and sequence.",0);
    else
        ProgressBar p_bar("Loading reference annotation.",0);

	GList<GSeqData> ref_rnas;
	
	// If the RefSequenceTable already has entries, we will sort the GTF records
	// according to their observation order.  Otherwise, we will sort the 
	// RefSequenceTable's records lexicographically.
	bool reorder_GTF_recs_lexicographically = false;
	if (rt.size() == 0)
	{
	  reorder_GTF_recs_lexicographically = true;
	}

	if (ref_mRNA_file)
	{
		gtf_tracking_verbose=cuff_verbose;
		read_transcripts(ref_mRNA_file, ref_rnas, gtf_crc_result, true);
	}
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
	GFastaHandler gfasta(fasta_dir.c_str());
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			//const char* name = ref_rnas[j]->gseq_name;
			
			int f = 0;
			int r = 0;
			int u  = 0;
			GffObj* rna_p;
			RefID ref_id = rt.get_id(ref_rnas[j]->gseq_name, NULL);
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
					faseq = gfasta.fetch(last_gseq_id);
					if (faseq==NULL)
					{
						fprintf(stderr,"This contig will not be bias corrected.\n");
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
				
				Scaffold ref_scaff(ref_id, strand, ops, true);
				
				char* rna_seq = 0;
				int seqlen=0;
				if (loadSeqs && faseq){ 
					rna_seq = rna.getSpliced(faseq, false, &seqlen);
				}

				if (rna.getID())
					ref_scaff.annotated_trans_id(rna.getID());
				
				
				if (rna.getGeneID())
					ref_scaff.annotated_gene_id(rna.getGeneID());
				
				if (rna.getGeneName())
					ref_scaff.annotated_gene_name(rna.getGeneName());
				
				
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
					string rs = (rna_seq) ? rna_seq:"";
					std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
					ref_scaff.seq(rs);
					GFREE(rna_seq);
				}
                
				boost::shared_ptr<Scaffold> scaff(new Scaffold());
                *scaff = ref_scaff;
                assert (scaff);
				ref_mRNAs.push_back(scaff); 
			}
		}
        
        BOOST_FOREACH (boost::shared_ptr<Scaffold> s, ref_mRNAs)
        {
            assert (s);
        }
        
        if (reorder_GTF_recs_lexicographically)
        {
            rt.order_recs_lexicographically();
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
    {
		return false;
    }
	
	// Update the bounds on the span
	if (hit.left() < _leftmost)
		_leftmost = hit.left();
	if (hit.right() > _rightmost)
		_rightmost = hit.right();
	
	
	_hits.push_back(hit);
	return true;
}

struct HitlessScaffold
{
	bool operator()(boost::shared_ptr<Scaffold> x)
	{
		return x->mate_hits().empty();
	}
};

bool unmapped_hit(const MateHit& x)
{
	return !(x.is_mapped());
}


bool HitBundle::add_open_hit(boost::shared_ptr<ReadGroupProperties const> rg_props,
                             const ReadHit* bh,
							 bool expand_by_partner)
{
    assert (bh != NULL);
    
	_leftmost = min(_leftmost, bh->left());
	_ref_id = bh->ref_id();
    
	if (bh->is_singleton() || no_read_pairs)
	{
		_rightmost = max(_rightmost, bh->right());
		MateHit m(rg_props, bh->ref_id(), bh, NULL);
        if (m.right() - m.left() > max_gene_length)
        {
            fprintf(stderr, "Warning: hit is longer than max_gene_length, skipping\n");
            return false;
        }
		add_hit(m);
	}
	else
	{
        if (abs(bh->right() - bh->partner_pos()+1) > max_gene_length)
        {
            fprintf(stderr, "Warning: hit is longer than max_gene_length, skipping\n");
            return false;
        }
		if (expand_by_partner)
			_rightmost = max(max(_rightmost, bh->right()), bh->partner_pos()+1);
		
		uint64_t search_key = bh->insert_id() ^ bh->left();
		std::pair<OpenMates::iterator, OpenMates::iterator> its = _open_mates.equal_range(search_key);

		// Does this hit close an open mate?
		bool found_partner = false;
		for(OpenMates::iterator it = its.first; it != its.second; ++it) {
		
			MateHit& pm = it->second;

			if(pm.left_alignment()->partner_pos() != bh->left())
				continue;

			if(pm.insert_id() != bh->insert_id())
				continue;

			{
					
				Scaffold L(MateHit(rg_props, bh->ref_id(), pm.left_alignment(), NULL));
				Scaffold R(MateHit(rg_props, bh->ref_id(), bh, NULL));
				
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
						_open_mates.erase(it);
						// Boost unordered_multimap ordinarily never shrinks.
						// Compel it to do so if significantly underloaded.
						if(_open_mates.size() > _rehash_size_threshold && _open_mates.load_factor() < _rehash_threshold) {
						  _open_mates.rehash(0);
						  if(_open_mates.load_factor() < _rehash_threshold) {
							// It didn't shrink -- looks like this Boost implementation has a minimum
							// size limit, or other reason to keep the table large. Don't keep rehashing with every remove op:
							_rehash_size_threshold = _open_mates.size();
						  }
						}
						found_partner = true;
						break;
					}

			}

		}

		if(!found_partner) {

			// Add it to the list of open mates, unless we would
			// already have seen it's partner
			if(bh->left() <= bh->partner_pos())
				{
					MateHit open_hit(rg_props,
									 bh->ref_id(), 
									 bh, 
									 NULL);
					
					uint64_t insert_key = bh->insert_id() ^ bh->partner_pos();
					_open_mates.insert(make_pair(insert_key, open_hit));
				}
			else
				{
					// This should never happen during hit_driven or ref_guided bundling, and in the case of
					// ref_driven, this read clearly shouldn't map to any of the transcripts anyways.
					// Adding this hit would cause problems with multi-reads that straddle boundaries after assembly.
					// add_hit(MateHit(rg_props,bh->ref_id(), bh, NULL));
					return false;
				}

		}

	}
    return true;
}

void HitBundle::collapse_hits()
{
	::collapse_hits(_hits, _non_redundant);
    //_non_redundant = _hits;
}

void HitBundle::finalize_open_mates()
{
    // We don't want to split reads accross boundaries since this would only occur
    // in ref_driven mode and the read shouldn't map to any of the references in this case.

    for(OpenMates::iterator itr = _open_mates.begin(); itr != _open_mates.end(); ++itr)
    {
		delete itr->second.left_alignment();
		delete itr->second.right_alignment();
    }
    _open_mates.clear();
}

void HitBundle::remove_hitless_scaffolds()
{
	vector<boost::shared_ptr<Scaffold> >::iterator new_end = remove_if(_ref_scaffs.begin(),
												   _ref_scaffs.end(),
												   HitlessScaffold());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());	
}



void HitBundle::combine(const vector<HitBundle*>& in_bundles,
                        HitBundle& out_bundle)
{
    out_bundle._hits.clear();
    out_bundle._non_redundant.clear();
    out_bundle._ref_scaffs.clear();
    
    for (size_t i = 1; i < in_bundles.size(); ++i)
    {
        assert(in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id());
    }
    
    // Merge  hits
    vector<size_t> indices(in_bundles.size(),0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit=NULL; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<MateHit>& curr_hits = in_bundles[i]->hits();
            
            if (indices[i] == curr_hits.size())
                continue;
            
            const MateHit* curr_hit = &curr_hits[indices[i]];
            
            if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
            {
                next_bundle = i;
                next_hit = curr_hit;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        out_bundle._hits.push_back(*next_hit);
        indices[next_bundle]++;
    }
    
    // Merge collapsed hits
    indices = vector<size_t>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit = NULL; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<MateHit>& curr_non_redundant_hits = in_bundles[i]->non_redundant_hits();
            
            if (indices[i] == curr_non_redundant_hits.size())
                continue;
            
            const MateHit* curr_hit = &curr_non_redundant_hits[indices[i]];
            
            if (next_bundle == -1 || mate_hit_lt(*curr_hit, *next_hit))
            {
                next_bundle = i;
                next_hit = curr_hit;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        out_bundle._non_redundant.push_back(*next_hit);
        indices[next_bundle]++;
    }
    
    for(size_t i = 0; i < in_bundles.size(); ++i)
    {
        for (size_t j = 0; j < in_bundles[i]->_ref_scaffs.size(); ++j)
        {
            in_bundles[i]->_ref_scaffs[j]->clear_hits();
        }
    }
    
    /*
    // Merge ref scaffolds
    indices = vector<size_t>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        boost::shared_ptr<Scaffold> next_scaff; 
        for(size_t i = 0; i < in_bundles.size(); ++i)
        {
            const vector<boost::shared_ptr<Scaffold> >& curr_scaffs = in_bundles[i]->_ref_scaffs;
            
            if (indices[i] == curr_scaffs.size())
                continue;
            
            boost::shared_ptr<Scaffold> curr_scaff = curr_scaffs[indices[i]];
            
            if (next_bundle == -1 || scaff_lt_rt_oplt(*curr_scaff, *next_scaff))
            {
                next_bundle = i;
                next_scaff = curr_scaff;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        if (out_bundle._ref_scaffs.size()==0 || out_bundle._ref_scaffs.back()->annotated_trans_id() != next_scaff->annotated_trans_id()) 
            out_bundle.add_ref_scaffold(next_scaff);
        indices[next_bundle]++;
    }
	*/
    
    for (size_t i = 0; i < in_bundles.size(); ++i)
    {
        for (size_t j = 0; j < in_bundles[i]->ref_scaffolds().size(); ++j)
        {
            out_bundle.add_ref_scaffold(in_bundles[i]->ref_scaffolds()[j]);
        }
    }
    
    sort(out_bundle._ref_scaffs.begin(), out_bundle._ref_scaffs.end(), scaff_lt_rt_oplt_sp);
    vector<boost::shared_ptr<Scaffold> >::iterator new_end = unique(out_bundle._ref_scaffs.begin(),
                                                             out_bundle._ref_scaffs.end(),
                                                             StructurallyEqualScaffolds());
    out_bundle._ref_scaffs.erase(new_end, out_bundle._ref_scaffs.end());
    vector<boost::shared_ptr<Scaffold> >(out_bundle._ref_scaffs).swap(out_bundle._ref_scaffs);
    
    out_bundle.finalize(true); // true means everything is already sorted, etc.
    out_bundle._num_replicates = (int)in_bundles.size();
}


void HitBundle::finalize(bool is_combined)
{
	_final = true;
    if (!is_combined)
	{
        // only perform read skipping on primary bundles 
        // (i.e. don't do it on bundles we're making by combining two or more other bundles)
        size_t num_skipped = _hits.size() * read_skip_fraction;
        if (num_skipped > 0 && num_skipped < _hits.size())
        {
            random_shuffle(_hits.begin(), _hits.end());
            for (int i = (int)_hits.size() - num_skipped; i >= 0 && i < (int)_hits.size(); ++i)
            {
                delete _hits[i].left_alignment();
                _hits[i].left_alignment(NULL);
                
                delete _hits[i].right_alignment();
                _hits[i].right_alignment(NULL);
            }
            _hits.resize(_hits.size() - num_skipped);
            is_combined = false;
        }
        else if (num_skipped >= _hits.size())
        {
            for (size_t i = 0; i < _hits.size(); ++i)
            {
                delete _hits[i].left_alignment();
                delete _hits[i].right_alignment();
            }
            _hits.clear();
        }

		sort(_hits.begin(), _hits.end(), mate_hit_lt);
        if (cond_prob_collapse)
        {
            collapse_hits();
        }
        else
        {
            BOOST_FOREACH (MateHit& hit, _hits)
            {
                hit.incr_collapse_mass(hit.internal_scale_mass());
            }
            _non_redundant = _hits;
            
        }
		sort(_ref_scaffs.begin(), _ref_scaffs.end(), scaff_lt_rt_oplt_sp);
		vector<boost::shared_ptr<Scaffold> >::iterator new_end = unique(_ref_scaffs.begin(), 
												_ref_scaffs.end(),
												StructurallyEqualScaffolds());
		_ref_scaffs.erase(new_end, _ref_scaffs.end());
        vector<boost::shared_ptr<Scaffold> >(_ref_scaffs).swap(_ref_scaffs);
	}
	
    for (size_t j = 0; j < _ref_scaffs.size(); ++j)
	{
		_ref_scaffs[j]->clear_hits();
	}
    
    _compatible_mass = 0.0;
    
	for (size_t i = 0; i < _hits.size(); ++i)
	{
		MateHit& hit = _hits[i];
		
		Scaffold hs(hit);
		
        if (i >= 1)
        {
            assert (hit.ref_id() == _hits[i-1].ref_id());
        }
		hit.is_mapped(false);
		for (size_t j = 0; j < _ref_scaffs.size(); ++j)
		{
			// add hit only adds if the hit is structurally compatible
			if (_ref_scaffs[j]->contains(hs))
			{
				bool added = _ref_scaffs[j]->add_hit(&hit);
                if (added)
                    hit.is_mapped(true);
			}
		}
        if (hit.is_mapped())
        {
            _compatible_mass += hit.internal_scale_mass();
        }
	}
    
}

void print_sort_error(const char* last_chr_name, 
                      int last_chr_pos, 
                      const char* bh_name, 
                      int bh_pos)
{
    fprintf(stderr, "\nError: this SAM file doesn't appear to be correctly sorted!\n");
    fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
            bh_name,
            bh_pos,
            last_chr_name,
            last_chr_pos);
    fprintf(stderr, "Cufflinks requires that if your file has SQ records in\nthe SAM header that they appear in the same order as the chromosomes names \nin the alignments.\nIf there are no SQ records in the header, or if the header is missing,\nthe alignments must be sorted lexicographically by chromsome\nname and by position.\n \n");
}


double BundleFactory::next_valid_alignment(const ReadHit*& bh)
{
    const char* hit_buf;
	size_t hit_buf_size = 0;
    bh = NULL;
    
	// Keep track of mass of hits we skip
	double raw_mass = 0; 
	
    while (true)
    {
    
        if (!_hit_fac->next_record(hit_buf, hit_buf_size))
            break;
        
        ReadHit tmp;
        if (!_hit_fac->get_hit_from_buf(hit_buf, tmp, false))
            continue;
        
		if (tmp.ref_id() == 12638153115695167477)  // corresponds to SAM "*" under FNV hash. unaligned read record 
            continue;
        
		raw_mass += tmp.mass();
		
        if (_hit_fac->ref_table().get_name(tmp.ref_id())==NULL) // unaligned read record (!?)
            continue;
            
        if (spans_bad_intron(tmp))
            continue;
        
        int order = _hit_fac->ref_table().observation_order(tmp.ref_id());
        if (_prev_pos != 0)
        {
            int prev_order = _hit_fac->ref_table().observation_order(_prev_ref_id);
            
            if (prev_order > order || (prev_order == order && _prev_pos > tmp.left()))
            {
                const char* bh_chr_name = _hit_fac->ref_table().get_name(tmp.ref_id());
                const char* last_bh_chr_name = _hit_fac->ref_table().get_name(_prev_ref_id);
                                
                print_sort_error(last_bh_chr_name, 
                                 _prev_pos, 
                                 bh_chr_name, 
                                 tmp.left());
                exit(1);
            }
        }
        
        _prev_ref_id = tmp.ref_id();
        _prev_pos = tmp.left();
        
        bool hit_within_mask = false;
        
        // We want to skip stuff that overlaps masked GTF records, so 
        // sync up the masking chromosome
        if (!mask_gtf_recs.empty() && 
            next_mask_scaff != mask_gtf_recs.end() &&
            (*next_mask_scaff)->ref_id() != tmp.ref_id())
        {
            bool found_scaff = false;
            vector<boost::shared_ptr<Scaffold> >::iterator curr_mask_scaff = mask_gtf_recs.begin();
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
        
        // if the user's asked for read trimming, do it here.
        if (trim_read_length > 0)
        {
            tmp.trim(trim_read_length);
        }
        
        bh = new ReadHit(tmp);
        
        break;
    }
    
    return raw_mass;
}

double BundleFactory::rewind_hit(const ReadHit* rh)
{
	double mass = rh->mass();
	delete rh;
	_hit_fac->undo_hit();
	return mass;
}

bool BundleFactory::next_bundle_hit_driven(HitBundle& bundle)
{
	const ReadHit* bh = NULL;
    
    bool skip_read = false;
    
	while(bh == NULL)
	{
		if (!_hit_fac->records_remain())
		{
			return false;
		}
        
        // If we are randomly throwing out reads, check to see
        // whether this one should be kept.
        if (bundle.hits().size() >= max_frags_per_bundle)
        {
            skip_read = true;
            next_valid_alignment(bh);
        }
        else
        {
            double raw_mass = next_valid_alignment(bh);
            if (bh && bh->num_hits() > max_frag_multihits)
            {
                skip_read = true;
            }
            else
            {
                bundle.add_raw_mass(raw_mass);
            }
        }
	}
	
	if ((skip_read || !bundle.add_open_hit(read_group_properties(), bh)) && bh != NULL)
    {
        delete bh;
        bh = NULL;
    }
	_expand_by_hits(bundle);

    assert(bundle.left() != -1);    
	bundle.finalize_open_mates();
	bundle.finalize();
    assert(bundle.right() != -1);
    
    return true;
}

bool BundleFactory::next_bundle_ref_driven(HitBundle& bundle)
{
	if (next_ref_scaff == ref_mRNAs.end())
	{
		const ReadHit* bh = NULL;
		while(_hit_fac->records_remain())
		{
            double raw_mass = next_valid_alignment(bh);
            if (bundle.hits().size() < max_frags_per_bundle)
            {
                if (bh && bh->num_hits() > max_frag_multihits)
                {
                    
                }
                else
                {
                    bundle.add_raw_mass(raw_mass);
                }
                if (bh) { delete bh; }
            }
            else
            {
                delete bh;
                bh = NULL;
            }
		    
		}
		bundle.finalize();
		return false;
	}
	
	bundle.add_ref_scaffold(*next_ref_scaff);
    
	++next_ref_scaff;
    
	_expand_by_refs(bundle);
    
//    for (size_t i = 0; i < bundle.ref_scaffolds().size(); ++i)
//    {
//        boost::shared_ptr<Scaffold> s = bundle.ref_scaffolds()[i];
//        if (s->annotated_gene_id() == "ENSG00000268467.1")
//        {
//                int a = 4;
//        }
//    }

    
	// The most recent RefID and position we've seen in the hit stream
	RefID last_hit_ref_id_seen = 0;
	int last_hit_pos_seen = 0;
	
	// include hits that lay within the bundle interval
	while(true)
	{		
		const ReadHit* bh = NULL;
        
        bool skip_read = false;
		// If we are randomly throwing out reads, check to see
        // whether this one should be kept.
        if (bundle.hits().size() >= max_frags_per_bundle)
        {
            next_valid_alignment(bh);
            skip_read = true;
        }
        else
        {
            double raw_mass = next_valid_alignment(bh);
            if (bh && bh->num_hits() > max_frag_multihits)
            {
                skip_read = true;
            }
            else
            {
                bundle.add_raw_mass(raw_mass);
            }
        }
        
        if (bh == NULL)
        {
			if (_hit_fac->records_remain())
				continue;
			else
				break;
        }
        
		last_hit_ref_id_seen = bh->ref_id();
		last_hit_pos_seen = bh->left();
		
		// test if the hit stream needs to catch up or has gone too far based on ref_id
		if (bh->ref_id() != bundle.ref_id())
		{
			int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
			int bundle_chr_order = _hit_fac->ref_table().observation_order(bundle.ref_id());
			
			if (bh_chr_order < bundle_chr_order) // the hit stream has not caught up, skip
			{
				delete bh;
                bh = NULL;
				continue; 
			}
			else // the hit stream has gone too far, rewind and break
			{
                rewind_hit(bh);
                bh = NULL;
                break;
			}
		}
        
        if (bh == NULL) // the hit stream has gone too far, break
            break;
		
        if (bh->left() >= bundle.left() && bh->right() <= bundle.right())
		{
            if (skip_read)
            {
                delete bh;
                bh = NULL;
            }
			else
            {
                if (!bundle.add_open_hit(read_group_properties(), bh, false))
                {
                    delete bh;
                    bh = NULL;
                }
            }
		}
		else if (bh->left() >= bundle.right())
		{
            if (skip_read == false)
            {
                bundle.rem_raw_mass(rewind_hit(bh));
                bh = NULL;
            }
            else
            {
                delete bh;
                bh = NULL;
            }
			break;
		}
	    else
        {
            // It's not within the bundle bounds, but it's also not past the 
            // right end, so skip it.
            delete bh;
            bh = NULL;
        }
        
        if (skip_read == true && bh != NULL)
        {
            delete bh;
            bh = NULL;
        }
	}
	
    assert(bundle.left() != -1);
    bundle.finalize_open_mates();
	bundle.finalize();
    assert(bundle.right() != -1);
    
    return true;
}

// NOTE: does not support read skipping yet or max hits per bundle yet.
bool BundleFactory::next_bundle_ref_guided(HitBundle& bundle)
{
	
	if (next_ref_scaff == ref_mRNAs.end())
	{
		return next_bundle_hit_driven(bundle);
	}
	
	const ReadHit* bh = NULL;
	while(bh == NULL)
	{
		if (!_hit_fac->records_remain())
		{
			return next_bundle_ref_driven(bundle);
		}
		bundle.add_raw_mass(next_valid_alignment(bh));
	}
	
	if (bh->ref_id() != (*next_ref_scaff)->ref_id())
	{
		int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
		int scaff_chr_order = _hit_fac->ref_table().observation_order((*next_ref_scaff)->ref_id());
		
		bundle.rem_raw_mass(rewind_hit(bh));
		bh = NULL;
        
		if (bh_chr_order < scaff_chr_order)
		{
			return next_bundle_hit_driven(bundle);
		}
		else
		{
			return next_bundle_ref_driven(bundle);
		}
	}
		
	if (bh->left() < (*next_ref_scaff)->left())
	{
		if (!bundle.add_open_hit(read_group_properties(), bh))
        {
            delete bh;
            bh = NULL;
        }
	}
	else 
	{
		bundle.rem_raw_mass(rewind_hit(bh));
        bh = NULL;
        
		bundle.add_ref_scaffold(*next_ref_scaff);
		next_ref_scaff++;
		_expand_by_refs(bundle);
	}
	
	while(_expand_by_hits(bundle) || 
		  _expand_by_refs(bundle)) {}
	
	assert(bundle.left() != -1);    
	bundle.finalize_open_mates();
	bundle.finalize();
	assert(bundle.right() != -1);
	
	return true; 
}

// expand the bundle interval as far as needed to include the overlapping
// chain of reference transcripts that also overlap the initial bundle
// interval
bool BundleFactory::_expand_by_refs(HitBundle& bundle)
{
	int initial_right = bundle.right();
	while(next_ref_scaff < ref_mRNAs.end())
	{		
		assert(bundle.ref_id() != (*next_ref_scaff)->ref_id() || (*next_ref_scaff)->left() >= bundle.left());
//        if (*next_ref_scaff && (*next_ref_scaff)->annotated_gene_id() == "XLOC_009372")
//        {
//            int a = 5;
//        }
		if (bundle.ref_id() == (*next_ref_scaff)->ref_id()
			&& overlap_in_genome((*next_ref_scaff)->left(),(*next_ref_scaff)->right(),bundle.left(), bundle.right()))
		{
			bundle.add_ref_scaffold(*next_ref_scaff);
            next_ref_scaff++;
		}
		else 
		{
			break;
		}		
	}

	
	return (bundle.right() > initial_right);
}

// expand bundle by chaining overlapping hits
bool BundleFactory::_expand_by_hits(HitBundle& bundle)
{
	int initial_right = bundle.right();
	while(true)
	{
        bool skip_read = false;
        const ReadHit* bh = NULL;
        
        double raw_mass = next_valid_alignment(bh);
        if (bh && bh->num_hits() > max_frag_multihits)
        {
            skip_read = true;
        }
        else
        {
            bundle.add_raw_mass(raw_mass);
        }

		if (bh == NULL)
		{
			if (_hit_fac->records_remain())
			{
				continue;
			}
			else
			{
				break;
			}	
		}
		
		if (bh->ref_id() == bundle.ref_id() && bh->left() < bundle.right() + olap_radius)
		{			
			if (skip_read || !bundle.add_open_hit(read_group_properties(), bh))
            {
                delete bh;
                bh = NULL;
            }
		}
		else
		{
			bundle.rem_raw_mass(rewind_hit(bh));

			break;
		}
	}
	
	return (bundle.right() > initial_right);
}

bool BundleFactory::next_bundle(HitBundle& bundle, bool cache_bundle)
{    
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(_factory_lock);
#endif
    bool got_bundle = false;
	switch(_bundle_mode)
	{
		case HIT_DRIVEN:
            _curr_bundle++;
			got_bundle = next_bundle_hit_driven(bundle);
            bundle.id(_curr_bundle);
			break;
		case REF_DRIVEN:
            _curr_bundle++;
			got_bundle = next_bundle_ref_driven(bundle);
            bundle.id(_curr_bundle);
			break;
		case REF_GUIDED:
            _curr_bundle++;
			got_bundle = next_bundle_ref_guided(bundle);
            bundle.id(_curr_bundle);
			break;
	}
	return got_bundle;
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
				
				if (read.num_hits() > 10)
				{
					itr->second.multimap_reads++;
				}
				
				if ( r_left <= small_anchor || (read_len - r_left) < small_anchor)
				{
					itr->second.little_reads++;
				}
                
                if (read.source_strand() == CUFF_FWD)
                {
                    //itr->second.fwd_strand_frags;
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
            case INS:
                g_left -= cig[i].length;
                break;
            case DEL:
                g_left += cig[i].length;
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
				//#if verbose_msg
				//							fprintf(stderr, "\t Filtering intron (due to overlap) %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
				//#endif	
				bool exists = binary_search(bad_introns.begin(), 
											bad_introns.end(), 
											itr->first);
				if (!exists)
				{
					verbose_msg("Filtering intron %d-%d spanned by %lu reads based on overlap with much more abundant intron: %d-%d spanned by %lu reads\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads,
							itr2->first.g_left(), 
							itr2->first.g_right(), 
							itr2->second.total_reads);
					
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
				verbose_msg("Filtering intron %d-%d spanned by %lu reads because %lg percent are multireads.\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads,
						multi_fraction * 100);
				
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
	
	BOOST_FOREACH (const MateHit& hit, bundle.hits())
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
	// [Geo] disable filtering of multi-mapped introns:
  // multimapping_introns(bundle.length(), bundle.left(), intron_counts, bad_introns, 0.5);
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
			double expected = success * spans.total_reads;

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
					verbose_msg("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
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
					verbose_msg("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
							itr->first.g_left(), 
							itr->first.g_right(), 
							itr->second.total_reads, 
							itr->second.little_reads, 
							expected,
							left_side_p,
							right_side_p);
					
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
				verbose_msg("Accepting intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
						itr->first.g_left(), 
						itr->first.g_right(), 
						itr->second.total_reads, 
						itr->second.little_reads, 
						expected,
						left_side_p,
						right_side_p);
				
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
            case INS:
                g_left -= cig[i].length;
                break;
            case DEL:
                g_left += cig[i].length;
                break;
			default:
				assert(false);
				break;
		}
	}
	
	return false;
}

void inspect_map(boost::shared_ptr<BundleFactory> bundle_factory,
                 BadIntronTable* bad_introns,
                 vector<LocusCount>& compatible_count_table,
                 vector<LocusCount>& total_count_table,
                 IdToLocusMap& id_to_locus_map,
                 bool progress_bar,
                 bool show_stats)
{

	ProgressBar p_bar;
	if (progress_bar)
		p_bar = ProgressBar("Inspecting reads and determining fragment length distribution.",bundle_factory->ref_table().size());
	RefID last_chrom = 0;

	long double map_mass = 0.0;
    long double norm_map_mass = 0.0;
	
	int min_len = numeric_limits<int>::max();
	int max_len = def_max_frag_len;
	vector<double> frag_len_hist(def_max_frag_len+1,0);
	bool has_pairs = false;

	int num_bundles = 0;
	size_t total_hits = 0;
	size_t total_non_redundant_hits = 0;
	
	//To be used for quartile normalization
	vector<long double> mass_dist; 	
	
	// Store the maximum read length for "first" and "second" reads to report to user.
	int max_1 = 0;
	int max_2 = 0;
	
	boost::shared_ptr<MultiReadTable> mrt(new MultiReadTable());
	
	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		bool valid_bundle = bundle_factory->next_bundle(*bundle_ptr, false);
		HitBundle& bundle = *bundle_ptr;

        if (use_compat_mass) //only count hits that are compatible with ref transcripts
        {
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.compatible_mass();
//            if (lib_norm_method == QUARTILE && bundle.compatible_mass() > 0)
//            {
//                mass_dist.push_back(bundle.compatible_mass());
//            }
        }
        else if (use_total_mass) //use all raw mass
        { 
            
            // Take raw mass even if bundle is "empty", since we could be out of refs
            // with remaining hits
            map_mass += bundle.raw_mass();
//            if (lib_norm_method == QUARTILE && bundle.raw_mass() > 0)
//            {
//                mass_dist.push_back(bundle.raw_mass());
//            }
        }
        else
        {
            fprintf(stderr, "Error: hit counting scheme for normalization is not set!\n");
            assert(false);
            exit(1);
        }
		
		const RefSequenceTable& rt = bundle_factory->ref_table();
		const char* chrom = rt.get_name(bundle.ref_id());
		char bundle_label_buf[2048];
        if (chrom)
        {
            sprintf(bundle_label_buf, "%s:%d-%d", chrom, bundle.left(), bundle.right());
            verbose_msg("Inspecting bundle %s with %lu reads\n", bundle_label_buf, bundle.hits().size());
            
            vector<string> gene_ids;
            vector<string> gene_short_names;
            BOOST_FOREACH(boost::shared_ptr<Scaffold> s, bundle.ref_scaffolds())
            {
                if (s->annotated_gene_id() != "")
                    gene_ids.push_back(s->annotated_gene_id());
                if (s->annotated_gene_name() != "")
                    gene_short_names.push_back(s->annotated_gene_name());
                
                if (s->annotated_gene_id() == "ENSG00000268467.1")
                {
                    int a = 4;
                }
                

                
            }
            compatible_count_table.push_back(LocusCount(bundle_label_buf, floor(bundle.compatible_mass()), bundle.ref_scaffolds().size(), gene_ids, gene_short_names));
            total_count_table.push_back(LocusCount(bundle_label_buf, floor(bundle.raw_mass()), bundle.ref_scaffolds().size(), gene_ids, gene_short_names));
		}
        
        if (!valid_bundle)
		{
			delete bundle_ptr;
			break;
		}
		num_bundles++;
        
        BOOST_FOREACH(boost::shared_ptr<Scaffold> s, bundle.ref_scaffolds()){
            id_to_locus_map.register_locus_to_id(s->annotated_trans_id(), bundle_label_buf);
            id_to_locus_map.register_locus_to_id(s->annotated_gene_id(), bundle_label_buf);
            
            if (s->annotated_tss_id().empty() == false)
            {
                id_to_locus_map.register_locus_to_id(s->annotated_tss_id(), bundle_label_buf);
            }
            
            if (s->annotated_protein_id().empty() == false)
            {
                id_to_locus_map.register_locus_to_id(s->annotated_protein_id(), bundle_label_buf);
            }
        }
        
        if (progress_bar) 
        {
			double inc_amt = last_chrom == bundle.ref_id() ? 0.0 : 1.0;
			p_bar.update(bundle_label_buf, inc_amt);
			last_chrom = bundle.ref_id();
        }
        
        if (bad_introns != NULL)
		{
			identify_bad_splices(bundle, *bad_introns);
		}
		
		const vector<MateHit>& hits = bundle.non_redundant_hits();
		if (hits.empty())
		{
			delete bundle_ptr;
			continue;
		}
		
		list<pair<int, int> > open_ranges;
		int curr_range_start = hits[0].left();
		int curr_range_end = numeric_limits<int>::max();
		int next_range_start = -1;
		
		total_non_redundant_hits += bundle.non_redundant_hits().size();
		total_hits += bundle.hits().size();
		
		// This first loop calclates the map mass and finds ranges with no introns
		// Note that we are actually looking at non-redundant hits, which is why we use collapse_mass
		// This loop will also add multi-reads to the MultiReads table 
		for (size_t i = 0; i < hits.size(); ++i) 
		{
			assert(hits[i].left_alignment());
            
            // Add to table if multi-read
			if (hits[i].is_multi())
			{
				mrt->add_hit(hits[i]);
			}
			
			// Find left length
			int left_len = hits[i].left_alignment()->right()-hits[i].left_alignment()->left();
			min_len = min(min_len, left_len);
			if (!hits[i].left_alignment()->contains_splice())
            {
				if (hits[i].left_alignment()->is_first())
                    max_1 = max(max_1, left_len);
                else
                    max_2 = max(max_2, left_len);
            }
			
			// Find right length
			if (hits[i].right_alignment())
			{
				int right_len = hits[i].right_alignment()->right()-hits[i].right_alignment()->left();
				min_len = min(min_len, right_len);
				if (!hits[i].right_alignment()->contains_splice())
                {
                    if (hits[i].right_alignment()->is_first())
                        max_1 = max(max_1, right_len);
                    else
                        max_2 = max(max_2, right_len);
                }
                has_pairs = true;
			}
			
			// Find fragment length
			if (bundle.ref_scaffolds().size()==1 && hits[i].is_pair())
			// Annotation provided and single isoform gene
			{
				int start, end, mate_length;
				boost::shared_ptr<Scaffold> scaff = bundle.ref_scaffolds()[0];
				if (scaff->map_frag(hits[i], start, end, mate_length))
				{
					if (mate_length >= min_len && mate_length <= max_len)
						frag_len_hist[mate_length] += hits[i].collapse_mass();
				}
			}
			else if (bundle.ref_scaffolds().empty())
			// No annotation provided.  Look for ranges.
			{
				if (hits[i].left() > curr_range_end)
				{
					if (curr_range_end - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, curr_range_end));
					curr_range_start = next_range_start;
					curr_range_end = numeric_limits<int>::max();
				}
				if (hits[i].left_alignment()->contains_splice())
				{
					if (hits[i].left() - curr_range_start > max_len)
						open_ranges.push_back(make_pair(curr_range_start, hits[i].left()-1));
					curr_range_start = max(next_range_start, hits[i].left_alignment()->right());
				}
				if (hits[i].right_alignment() && hits[i].right_alignment()->contains_splice())
				{
					assert(hits[i].right_alignment()->left() >= hits[i].left());
					curr_range_end = min(curr_range_end, hits[i].right_alignment()->left()-1);
					next_range_start = max(next_range_start, hits[i].right());
				}
			}
		}
        
        if (bundle.ref_scaffolds().empty() && has_pairs) // No annotation provided
		{
			pair<int, int> curr_range(-1,-1);
			
			// This second loop uses the ranges found above to find the estimated frag length distribution
			// It also finds the minimum read length to use in the linear interpolation
			for (size_t i = 0; i < hits.size(); ++i)
			{
				if (hits[i].left() > curr_range.second && open_ranges.empty())
					break;
				
				if (hits[i].left() > curr_range.second)
				{
					curr_range = open_ranges.front();
					open_ranges.pop_front();
				}
				
				if (hits[i].left() >= curr_range.first && hits[i].right() <= curr_range.second && hits[i].is_pair())
				{
					int mate_len = hits[i].right()-hits[i].left();
					if (mate_len <= max_len)
						frag_len_hist[mate_len] += hits[i].collapse_mass();
				}
			}
		}
		
        open_ranges.clear();
		delete bundle_ptr;
	}
	
    norm_map_mass = map_mass;
    
//	if (lib_norm_method == QUARTILE && mass_dist.size() > 0)
//	{
//		sort(mass_dist.begin(),mass_dist.end());
//		int upper_quart_index = mass_dist.size() * 0.75;
//		norm_map_mass = mass_dist[upper_quart_index];
//	}

    if (bad_introns != NULL)
    {
        size_t alloced = 0;
        size_t used = 0;
        size_t num_introns = 0;
        for (BadIntronTable::const_iterator itr = bad_introns->begin();
             itr != bad_introns->end();
             ++itr)
        {
            alloced += itr->second.capacity() * sizeof(AugmentedCuffOp);
            used += itr->second.size() * sizeof(AugmentedCuffOp);
            num_introns += itr->second.size();
        }
        
        verbose_msg( "Bad intron table has %lu introns: (%lu alloc'd, %lu used)\n", num_introns, alloced, used);
    	verbose_msg( "Map has %lu hits, %lu are non-redundant\n", total_hits, total_non_redundant_hits);
    } 
    
	if (progress_bar)
		p_bar.complete();
	
	vector<double> frag_len_pdf(max_len+1, 0.0);
	vector<double> frag_len_cdf(max_len+1, 0.0);
    long double tot_count = accumulate(frag_len_hist.begin(), frag_len_hist.end(), 0.0 );
    bool empirical = false;
	
	if (user_provided_fld && has_pairs && tot_count >= 10000)
	{
		fprintf(stderr, "Warning: Overriding empirical fragment length distribution with user-specified parameters is not recommended.\n");
	}
	
	if (!has_pairs || tot_count < 10000)
	{
		if (has_pairs && !user_provided_fld)
		{
			fprintf(stderr, "Warning: Using default Gaussian distribution due to insufficient paired-end reads in open ranges.  It is recommended that correct parameters (--frag-len-mean and --frag-len-std-dev) be provided.\n");
		}
		tot_count = 0;
		normal frag_len_norm(def_frag_len_mean, def_frag_len_std_dev);
		max_len = def_frag_len_mean + 3*def_frag_len_std_dev;
		for(int i = min_len; i <= max_len; i++)
		{
			frag_len_hist[i] = cdf(frag_len_norm, i+0.5)-cdf(frag_len_norm, i-0.5);
			tot_count += frag_len_hist[i];
		}
	}
	else
	// Calculate the max frag length and interpolate all zeros between min read len and max frag len
	{	
		empirical = true;
		double curr_total = 0;
		size_t last_nonzero = min_len-1;
		for(size_t i = last_nonzero+1; i < frag_len_hist.size(); i++)
		{
			if (frag_len_hist[i] > 0)
			{
				if (last_nonzero != i-1)
				{
					double b = frag_len_hist[last_nonzero];
					double m = (frag_len_hist[i] - b)/(i-last_nonzero);
					for (size_t x = 1; x < i - last_nonzero; x++)
					{
						frag_len_hist[last_nonzero+x] = m * x + b;
						tot_count += frag_len_hist[last_nonzero+x];
						curr_total += frag_len_hist[last_nonzero+x];
					}	
				}
				last_nonzero = i;
			}
			
			curr_total += frag_len_hist[i];
			
			if (curr_total/tot_count > 0.9999)
			{
				max_len = i; 
				tot_count = curr_total;
				break;
			}
		}
	}
	
    double mean = 0.0;

    if (output_fld)
    {
        FILE* fhist = fopen(string(output_dir + "/frag_len_hist.csv").c_str(),"w");
        fprintf(fhist, "Length,Count\n");
        for(size_t i = 1; i < frag_len_hist.size(); i++)
        {
            fprintf(fhist, "%zu,%f\n", i, frag_len_hist[i]);
        }
        fclose(fhist);
    }

	// Convert histogram to pdf and cdf, calculate mean
	int frag_len_mode = 0;
	for(size_t i = min_len; i <= (size_t)max_len; i++)
	{
		frag_len_pdf[i] = frag_len_hist[i]/tot_count;
		frag_len_cdf[i] = frag_len_cdf[i-1] + frag_len_pdf[i];
        
		if (frag_len_pdf[i] > frag_len_pdf[frag_len_mode])
			frag_len_mode = i;
        mean += frag_len_pdf[i] * i;
	}
    
    double std_dev =  0.0;
    for(size_t i = 1; i < frag_len_hist.size(); i++)
    {
        std_dev += frag_len_pdf[i] * ((i - mean) * (i - mean));
    }
    
    std_dev = sqrt(std_dev);
	
	boost::shared_ptr<ReadGroupProperties> rg_props = bundle_factory->read_group_properties();

    FLDSource source = DEFAULT;
    if (empirical)
    {
        source = LEARNED;
    }
    else if (user_provided_fld)
    {
        source = USER;
    }

	boost::shared_ptr<EmpDist const> fld(new EmpDist(frag_len_pdf, frag_len_cdf, frag_len_mode, mean, std_dev, min_len, max_len, source));
	rg_props->multi_read_table(mrt);
	rg_props->frag_len_dist(fld);
	rg_props->normalized_map_mass(norm_map_mass);
    rg_props->total_map_mass(map_mass);
    rg_props->mode_transcript_coverage(bundle_factory->mode_transcript_coverage());
    
    if (show_stats)
    {
        fprintf(stderr, "> Map Properties:\n");
        //if (lib_norm_method == QUARTILE)
        //    fprintf(stderr, ">\tUpper Quartile: %.2Lf\n", norm_map_mass);
        fprintf(stderr, ">\tNormalized Map Mass: %.2Lf\n", norm_map_mass);
        fprintf(stderr, ">\tRaw Map Mass: %.2Lf\n", map_mass);
        if (corr_multi)
            fprintf(stderr,">\tNumber of Multi-Reads: %zu (with %zu total hits)\n", mrt->num_multireads(), mrt->num_multihits()); 
    //	if (has_pairs)
    //		fprintf(stderr, ">\tRead Type: %dbp x %dbp\n", max_1, max_2);
    //	else
    //		fprintf(stderr, ">\tRead Type: %dbp single-end\n", max(max_1,max_2));

        if (empirical)
        {
            fprintf(stderr, ">\tFragment Length Distribution: Empirical (learned)\n");
            fprintf(stderr, ">\t              Estimated Mean: %.2f\n", mean);
            fprintf(stderr, ">\t           Estimated Std Dev: %.2f\n", std_dev);
        }
        else
        {
            if (user_provided_fld)
            {
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (user-specified)\n");
            }
            else
            {
                fprintf(stderr, ">\tFragment Length Distribution: Truncated Gaussian (default)\n");
            }
            fprintf(stderr, ">\t              Default Mean: %d\n", def_frag_len_mean);
            fprintf(stderr, ">\t           Default Std Dev: %d\n", def_frag_len_std_dev);
        }
    }
	bundle_factory->num_bundles(num_bundles);
	bundle_factory->reset();
	return;
}

//////////////////////


bool PrecomputedExpressionBundleFactory::next_bundle(HitBundle& bundle, bool cache_bundle)
{
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(_factory_lock);
#endif
    bool got_bundle = BundleFactory::next_bundle(bundle, cache_bundle);
    if (got_bundle)
    {
        RefSequenceTable& rt = ref_table();
        
        char bundle_label_buf[2048];
        sprintf(bundle_label_buf, "%s:%d-%d", rt.get_name(bundle.ref_id()),	bundle.left(), bundle.right());
        
        boost::shared_ptr<const AbundanceGroup> ab = _hit_fac->next_locus(bundle.id(), cache_bundle);
        if (ab)
        {
            double compatible_mass = _hit_fac->get_compat_mass(bundle.id());
            double total_mass  = _hit_fac->get_total_mass(bundle.id());
            
            /*
            double compatible_mass = _hit_fac->get_compat_mass(bundle_label_buf);
            double total_mass  = _hit_fac->get_total_mass(bundle_label_buf);
            */
            
            bundle.finalize();
            bundle.add_raw_mass(total_mass);
            bundle.compatible_mass(compatible_mass);
            
            // calculate the depth of coverage for each transcript and store it in the factory
            // in case we need it later for library size normalization
//            size_t max_FPKM_idx = 0;
//            double max_FPKM = -1;
            for (size_t i = 0; i < ab->abundances().size(); ++i)
            {
                double num_frags = ab->abundances()[i]->num_fragments();
                double length = ab->abundances()[i]->effective_length();
                double FPKM =ab-> abundances()[i]->FPKM();
                // FIXME: Should probably set these "magic" params as options somewhere.
                if (num_frags >= 1 && length != 0 && FPKM >= 1e-5)
                {
                    double coverage = num_frags / length;
                    transcript_coverages.push_back(coverage);
                }

            }
            
//            if (max_FPKM_idx < ab->abundances().size())
//            {
//            }
            
            //fprintf (stderr, "Reconstituting bundle %s (%d) with mass %lf\n", bundle_label_buf, bundle.id(), compatible_mass);
            if (bundle.ref_scaffolds().size() != ab->abundances().size())
            {
                fprintf (stderr, "Error in file %s: reconstituted expression bundle %s (%lu transcripts)  does not match GTF (%d transcripts):\n", read_group_properties()->file_path().c_str(),  bundle_label_buf, ab->abundances().size(), bundle.ref_scaffolds().size());
                fprintf(stderr, "Reconstituted:\n");
                for (size_t i = 0; i < ab->abundances().size(); ++i)
                {
                    fprintf(stderr, "%s\n", ab->abundances()[i]->description().c_str());
                }
                fprintf(stderr, "GTF:\n");
                for (size_t i = 0; i < bundle.ref_scaffolds().size(); ++i)
                {
                    fprintf(stderr, "%s\n", bundle.ref_scaffolds()[i]->annotated_trans_id().c_str());
                }
                exit(1);
            }
        }
        else
        {
            fprintf (stderr, "Error: no abundance info for locus %s\n", bundle_label_buf);
        }
        
    }
    return got_bundle;
}

boost::shared_ptr<const AbundanceGroup> PrecomputedExpressionBundleFactory::get_abundance_for_locus(int locus_id)
{
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(_factory_lock);
#endif
    return _hit_fac->get_abundance_for_locus(locus_id);
}

void PrecomputedExpressionBundleFactory::clear_abundance_for_locus(int locus_id)
{
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(_factory_lock);
#endif
    _hit_fac->clear_abundance_for_locus(locus_id);
}


