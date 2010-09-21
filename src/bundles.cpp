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

using namespace std;
using boost::math::binomial;

//struct ScaffoldSorter
//{
//	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
//	bool operator()(shared_ptr<Scaffold const> lhs, shared_ptr<Scaffold const> rhs)
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
	bool operator()(shared_ptr<Scaffold const> lhs, shared_ptr<Scaffold const> rhs)
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

/*
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
*/

//FIXME: needs refactoring
void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<shared_ptr<Scaffold> >& ref_mRNAs,
				   bool loadSeqs,
				   bool loadFPKM) 
{

	
	if (loadSeqs)
		ProgressBar p_bar("Loading reference annotation and sequence",0);
    
	GList<GSeqData> ref_rnas;
	
	if (ref_mRNA_file)
	{
		read_transcripts(ref_mRNA_file, ref_rnas);
        //read_mRNAs(ref_mRNA_file, ref_rnas, &ref_rnas, true, -1, "", false);
	}
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
	GFastaHandler gfasta(fasta_dir.c_str());
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
				
				Scaffold ref_scaff(ref_id, strand, ops);
				
				char* rna_seq;
				int seqlen=0;
				if (loadSeqs && faseq){ 
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
					string rs = (rna_seq) ? rna_seq:"";
					std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
					ref_scaff.seq(rs);
					GFREE(rna_seq);
				}
                
				shared_ptr<Scaffold> scaff(new Scaffold());
                *scaff = ref_scaff;
                assert (scaff);
				ref_mRNAs.push_back(scaff); 
			}
		}
        
        foreach (shared_ptr<Scaffold> s, ref_mRNAs)
        {
            assert (s);
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
	bool operator()(shared_ptr<Scaffold> x)
	{
		return x->mate_hits().empty();
	}
};

bool unmapped_hit(const MateHit& x)
{
	return !(x.is_mapped());
}


void HitBundle::add_open_hit(shared_ptr<ReadGroupProperties const> rg_props,
                             const ReadHit* bh)
{
	if (bh->partner_ref_id() == 0 
		|| bh->partner_ref_id() != bh->ref_id() ||
		abs((int)bh->partner_pos() - (int)bh->left()) > max_partner_dist)
	{
		// This is a singleton, so just make a closed MateHit and
		// continue
		MateHit m(rg_props, bh->ref_id(), bh, NULL);
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
                                 NULL);
				
				pair<OpenMates::iterator, bool> ret;
				ret = _open_mates.insert(make_pair(bh->partner_pos(), 
												  list<MateHit>()));
				
				ret.first->second.push_back(open_hit);
			}
			else
			{
				add_hit(MateHit(rg_props,bh->ref_id(), bh, NULL));
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
					MateHit open_hit(rg_props, bh->ref_id(), bh, NULL);
					
					pair<OpenMates::iterator, bool> ret;
					ret = _open_mates.insert(make_pair(bh->partner_pos(), 
													  list<MateHit>()));
					
					ret.first->second.push_back(open_hit);
				}
				else
				{
					add_hit(MateHit(rg_props, bh->ref_id(), bh, NULL));
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
												   HitlessScaffold());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());	
}

void HitBundle::remove_unmapped_hits()
{
	
	foreach (MateHit& hit, _hits)
	{
		if (unmapped_hit(hit))
		{
			delete hit.left_alignment();
			delete hit.right_alignment();
		} 
	}
	
	vector<MateHit>::iterator new_end = remove_if(_hits.begin(),
												  _hits.end(),
												  unmapped_hit);

	_hits.erase(new_end, _hits.end());	

	new_end = remove_if(_non_redundant.begin(),
						_non_redundant.end(),
						unmapped_hit);
	_non_redundant.erase(new_end, _non_redundant.end());
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
    vector<int> indices(in_bundles.size(),0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit; 
        for(int i = 0; i < in_bundles.size(); ++i)
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
    indices = vector<int>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        const MateHit* next_hit; 
        for(int i = 0; i < in_bundles.size(); ++i)
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
    
    for(int i = 0; i < in_bundles.size(); ++i)
    {
        for (size_t j = 0; j < in_bundles[i]->_ref_scaffs.size(); ++j)
        {
            in_bundles[i]->_ref_scaffs[j]->clear_hits();
        }
    }
    
    // Merge ref scaffolds
    indices = vector<int>(in_bundles.size(), 0);
    while(true)
    {
        int next_bundle = -1;
        shared_ptr<Scaffold> next_scaff; 
        for(int i = 0; i < in_bundles.size(); ++i)
        {
            const vector<shared_ptr<Scaffold> >& curr_scaffs = in_bundles[i]->_ref_scaffs;
            
            if (indices[i] == curr_scaffs.size())
                continue;
            
            shared_ptr<Scaffold> curr_scaff = curr_scaffs[indices[i]];
            
            if (next_bundle == -1 || scaff_lt_rt_oplt(*curr_scaff, *next_scaff))
            {
                next_bundle = i;
                next_scaff = curr_scaff;
            }
        }
        
        if(next_bundle==-1)
            break;
        
        if (out_bundle._ref_scaffs.size()==0 || out_bundle._ref_scaffs.back()->annotated_trans_id() != next_scaff->annotated_trans_id()) 
            out_bundle._ref_scaffs.push_back(next_scaff);
        indices[next_bundle]++;
    }
    
    out_bundle.finalize(true); // true means everything is already sorted, etc.
    out_bundle._num_replicates = (int)in_bundles.size();
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
	}
	
	for (size_t i = 0; i < _non_redundant.size(); ++i)
	{
		MateHit& hit = _non_redundant[i];
		
		Scaffold hs(hit);
		
        if (i >= 1)
        {
            assert (hit.ref_id() == _non_redundant[i-1].ref_id());
        }
		hit.is_mapped(false);
		for (size_t j = 0; j < _ref_scaffs.size(); ++j)
		{
			// add hit only adds if the hit is structurally compatible
			if (_ref_scaffs[j]->contains(hs))
			{
                hit.is_mapped(true);
			}
		}
	}

	remove_unmapped_hits();
    
	if (_ref_scaffs.size() > 0)
    {
        _leftmost = INT_MAX;
        _rightmost = -1;
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

const ReadHit* BundleFactory::next_valid_alignment()
{
    const char* hit_buf;
	size_t hit_buf_size = 0;
    const ReadHit* bh = NULL;
    
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
        
        int order = _hit_fac->ref_table().observation_order(tmp.ref_id());
        if (_prev_pos != 0)
        {
            int prev_order = _hit_fac->ref_table().observation_order(_prev_ref_id);
            
            if (prev_order > order || (prev_order == order && _prev_pos > tmp.left()))
            {
                const char* bh_chr_name = _hit_fac->ref_table().get_name(tmp.ref_id());
                const char* last_bh_chr_name = _hit_fac->ref_table().get_name(_prev_ref_id);
                
                _hit_fac->ref_table().observation_order(_prev_ref_id);
                
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
        bh = new ReadHit(tmp);
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
        const ReadHit* bh = next_valid_alignment();
        

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
                    
                    int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
                    int ref_rna_chr_order = _hit_fac->ref_table().observation_order((*next_ref_scaff)->ref_id());
                    int last_bh_chr_order = _hit_fac->ref_table().observation_order(last_hit_ref_id_seen);
                    
                    // if the BAM/SAM file isn't lexicographically sorted by chromosome, print an error
                    // and exit.
                    if (last_bh_chr_name && last_bh_chr_order > bh_chr_order)
                    { 
                        print_sort_error(last_bh_chr_name, 
                                         last_hit_pos_seen, 
                                         bh_chr_name, 
                                         bh->left());
                        exit(1);
                    }
                    
                    if (bh_chr_name && ref_rna_chr_name)
                    {
                        if (bh_chr_order < ref_rna_chr_order)
                        {
                            // hit is lexicographically less than the reference,
                            // so skip the hit.
                            delete bh;
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
                    delete bh;
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
                    delete bh;
                    continue;
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
                    delete bh;
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
                    delete bh;
                    return false;
                }
                else 
                {
                    // if we get here, there are more records, but this one is
                    // null, and we've reached the end of the reference trancripts
                    // or there weren't any
                    delete bh;
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
            
            int bh_chr_order = _hit_fac->ref_table().observation_order(bh->ref_id());
            int last_chr_order = _hit_fac->ref_table().observation_order(last_hit_ref_id_seen);
            
            if (last_chr_order > bh_chr_order)
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
                delete bh;
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
            
            last_hit_ref_id_seen = bh->ref_id();
            last_hit_pos_seen = bh->left();
            
			bundle.add_open_hit(read_group_properties(), bh);
		}
		else if (past_right_end)
		{
            delete bh;
			_hit_fac->undo_hit();
			break;
		}
        else 
        {
            // It's not within the bundle bounds, but it's also not past the 
            // right end, so skip it.
            last_hit_ref_id_seen = bh->ref_id();
            last_hit_pos_seen = bh->left();
            delete bh;
        }

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
					asm_printf("Filtering intron %d-%d spanned by %lu reads based on overlap with much more abundant intron: %d-%d spanned by %lu reads\n", 
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
				asm_printf("Filtering intron %d-%d spanned by %lu reads because %lg percent are multireads.\n", 
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
					asm_printf("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
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
					asm_printf("Filtering intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
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
				asm_printf("Accepting intron %d-%d spanned by %lu reads (%lu low overhang, %lg expected) left P = %lg, right P = %lg\n", 
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
			default:
				assert(false);
				break;
		}
	}
	
	return false;
}
