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

#include "common.h"
#include "bundles.h"
#include "scaffolds.h"

using namespace std;

struct ScaffoldSorter
{
	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
	bool operator()(const Scaffold& lhs, const Scaffold& rhs)
	{
		const char* lhs_name = rt.get_name(lhs.ref_id());
		const char* rhs_name = rt.get_name(rhs.ref_id());
		int c = strcmp(lhs_name, rhs_name);
		if (c != 0)
		{
			return c < 0;
		}
		else
		{
			return lhs.left() < rhs.left();
		}
	}
	
	RefSequenceTable& rt;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<Scaffold>& ref_mRNAs) 
{
	GList<GSeqData> ref_rnas;
	
	if (ref_mRNA_file)
	{
		//read_mRNAs(ref_mRNA_file, false, ref_rnas, ref_rnas, NULL, -1, false);
		read_mRNAs(ref_mRNA_file, ref_rnas);
	}
	
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			char* name = GffObj::names->gseqs.getName(ref_rnas[j]->gseq_id);
			RefID ref_id = rt.get_id(name, NULL);
			for (int i = 0; i < ref_rnas[j]->mrnas_f.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_f[i]);
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
				
				Scaffold ref_scaff(ref_id, CUFF_FWD, ops);
				if (rna.getID())
				{
					ref_scaff.annotated_trans_id(rna.getID());
				}
				if (rna.getGene())
					ref_scaff.annotated_gene_id(rna.getGene());
				char* short_name = rna.getAttr("gene_name");
				if (short_name)
				{
					ref_scaff.annotated_gene_name(short_name);
				}
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
				{
					ref_scaff.annotated_protein_id(protein_id);
				}
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
				{
					ref_scaff.annotated_tss_id(tss_id);
				}
				
				ref_mRNAs.push_back(ref_scaff); 
			}
			
			for (int i = 0; i < ref_rnas[j]->mrnas_r.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_r[i]);
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
				
				Scaffold ref_scaff(ref_id, CUFF_REV, ops);
				if (rna.getID())
				{
					ref_scaff.annotated_trans_id(rna.getID());
				}
				if (rna.getGene())
					ref_scaff.annotated_gene_id(rna.getGene());
				char* short_name = rna.getAttr("gene_name");
				if (short_name)
				{
					ref_scaff.annotated_gene_name(short_name);
				}
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
				{
					ref_scaff.annotated_protein_id(protein_id);
				}
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
				{
					ref_scaff.annotated_tss_id(tss_id);
				}
				
				ref_mRNAs.push_back(ref_scaff); 
			}
		}
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
	}
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
	bool operator()(const Scaffold& x)
	{
		return x.mate_hits().empty();
	}
};

void HitBundle::add_open_hit(shared_ptr<ReadHit> bh)
{
	if (bh->partner_ref_id() == 0 
		|| bh->partner_ref_id() != bh->ref_id() ||
		abs((int)bh->partner_pos() - (int)bh->left()) > max_partner_dist)
	{
		// This is a singleton, so just make a closed MateHit and
		// continue
		MateHit m(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0);
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
				MateHit open_hit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), inner_dist_mean, max_inner_dist);
				
				pair<OpenMates::iterator, bool> ret;
				ret = _open_mates.insert(make_pair(bh->partner_pos(), 
												  list<MateHit>()));
				
				ret.first->second.push_back(open_hit);
			}
			else
			{
				add_hit(MateHit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), inner_dist_mean, max_inner_dist));
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
					
					Scaffold L(MateHit(bh->ref_id(), pm.left_alignment(), shared_ptr<ReadHit const>(), 0, 0));
					Scaffold R(MateHit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0));
					
					bool strand_agree = L.strand() == CUFF_STRAND_UNKNOWN ||
					R.strand() == CUFF_STRAND_UNKNOWN ||
					L.strand() == R.strand();
					
					bool orientation_agree = pm.left_alignment()->antisense_align() != bh->antisense_align();
					
					if (strand_agree && 
						orientation_agree && (!Scaffold::overlap_in_genome(L, R, olap_radius) ||
											  Scaffold::compatible(L,R, false)))
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
					MateHit open_hit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), inner_dist_mean, max_inner_dist);
					
					pair<OpenMates::iterator, bool> ret;
					ret = _open_mates.insert(make_pair(bh->partner_pos(), 
													  list<MateHit>()));
					
					ret.first->second.push_back(open_hit);
				}
				else
				{
					add_hit(MateHit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), inner_dist_mean, max_inner_dist));
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
	}
}

void HitBundle::remove_hitless_scaffolds()
{
	vector<Scaffold>::iterator new_end = remove_if(_ref_scaffs.begin(),
												   _ref_scaffs.end(),
												   HitLessScaffold());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());	
}

void HitBundle::finalize()
{
	_final = true;
	
	finalize_open_mates();
	
	sort(_hits.begin(), _hits.end(), mate_hit_lt);
	collapse_hits(_hits, _non_redundant, _collapse_counts);
	
	sort(_ref_scaffs.begin(), _ref_scaffs.end(), scaff_lt);
	vector<Scaffold>::iterator new_end = unique(_ref_scaffs.begin(), 
												_ref_scaffs.end(),
												StructurallyEqualScaffolds());
	_ref_scaffs.erase(new_end, _ref_scaffs.end());
	
	for (size_t i = 0; i < _hits.size(); ++i)
	{
		const MateHit* hit = &(_hits[i]);
		
		Scaffold hs(*hit);
		
		for (size_t j = 0; j < _ref_scaffs.size(); ++j)
		{
			// add hit only adds if the hit is structurally compatible
			if (_ref_scaffs[j].Scaffold::contains(hs) &&
				Scaffold::compatible(_ref_scaffs[j], hs, true))
			{
				_ref_scaffs[j].add_hit(hit);
			}
		}
	}
	
	for (size_t i = 0; i < _ref_scaffs.size(); ++i)
	{
		if (_ref_scaffs[i].left() < _leftmost)
			_leftmost = _ref_scaffs[i].left();
		if (_ref_scaffs[i].right() > _rightmost)
			_rightmost = _ref_scaffs[i].right();
	}
}

bool BundleFactory::next_bundle(HitBundle& bundle_out)
{
	HitBundle bundle;
	
	if (feof(hit_file))
	{
		return false;
	}
	char bwt_buf[2048];
	
	RefID last_ref_id_seen = 0;
	RefID first_ref_id_seen = 0;
	
	int last_pos_seen = 0;

	off_t curr_pos = ftello(hit_file);
	
	int right_bundle_boundary = 0;
	
	int left_boundary = -1;
	while (fgets(bwt_buf, 2048, hit_file))
	{
		_next_line_num++;
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		shared_ptr<ReadHit> bh(new ReadHit());
		
		if (!sam_hit_fac.get_hit_from_buf(_next_line_num, bwt_buf, *bh, false))
		{
			continue;
		}
		
		if (bh->error_prob() > max_phred_err_prob)
			continue;
		
		if (bh->ref_id() == 84696373) // corresponds to SAM "*" under FNV hash. unaligned read record  
			continue;
		
		if (left_boundary == -1)
			left_boundary = bh->left();
		
		if (!ref_mRNAs.empty() && 
			next_ref_scaff != ref_mRNAs.end() &&
			next_ref_scaff->ref_id() != bh->ref_id())
		{
			for (size_t i = 0; i < _ref_scaff_offsets.size(); ++i)
			{
				if (_ref_scaff_offsets[i].first == bh->ref_id())
				{
					next_ref_scaff = _ref_scaff_offsets[i].second;
				}
			}
		}
		
		// break the bundle if there's a coverage gap or if alignments for a 
		// new contig are encountered.
		
		bool hit_within_boundary = false;
		
		if (!ref_mRNAs.empty())
		{
			//check that we aren't sitting in the middle of an annotated scaffold
			while (next_ref_scaff != ref_mRNAs.end() && 
				   next_ref_scaff->ref_id() == bh->ref_id() &&
				   next_ref_scaff->right() <= bh->left())
			{
				if (next_ref_scaff->left() >= bh->left())
				{
					break;
				}
				next_ref_scaff++;
			}
			
			while (next_ref_scaff != ref_mRNAs.end() && 
				   (!last_ref_id_seen || bh->ref_id() == last_ref_id_seen) &&
				   next_ref_scaff->ref_id() == bh->ref_id() &&
                   right_bundle_boundary > next_ref_scaff->left() && 
				   next_ref_scaff->left() <= bh->left() &&
				   next_ref_scaff->right() >= bh->right())
			{
				hit_within_boundary = true;
				right_bundle_boundary = max(right_bundle_boundary, next_ref_scaff->right());
				next_ref_scaff++;
			}
		}
		
		if (last_ref_id_seen == 0)
			hit_within_boundary = true;
		else if (bh->ref_id() == last_ref_id_seen && bh->left() <= right_bundle_boundary)
			hit_within_boundary = true;

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
			
			if (bh->right() + olap_radius - left_boundary < (int)max_gene_length)
			{
				right_bundle_boundary = max(right_bundle_boundary, bh->right() + olap_radius);
			}
			
			bool singleton = (bh->partner_ref_id() == 0 
							  || bh->partner_ref_id() != bh->ref_id() ||
							  abs(bh->partner_pos() - bh->left()) > max_partner_dist);
			if (!singleton)
			{
				if ((int)bh->partner_pos() + olap_radius - (int)bh->left() < (int)max_partner_dist)
					right_bundle_boundary = max(right_bundle_boundary, bh->partner_pos() + olap_radius);
			}
			
			bundle.add_open_hit(bh);
		}
		else
		{
			fseeko(hit_file, curr_pos, SEEK_SET);
			break;
		}
		
		if (!last_ref_id_seen)
		{
			first_ref_id_seen = bh->ref_id();
		}
		
		last_ref_id_seen = bh->ref_id();
		last_pos_seen = bh->left();
		
		curr_pos = ftello(hit_file);
	}
	
	bundle.finalize_open_mates();
	
	if (!ref_mRNAs.empty())
	{
		vector<Scaffold>::iterator itr = next_ref_scaff;
		while(itr < ref_mRNAs.end())
		{
			if (itr->ref_id() != last_ref_id_seen || itr->left() >= right_bundle_boundary)
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
				if (itr->ref_id() == last_ref_id_seen)
				{
					// Now we're on the right scaffold, gotta get to the right
					// coords
					if (overlap_in_genome(itr->left(), 
									  itr->right(), 
									  bundle.left(), 
									  bundle.right()))
					{	
						bundle.add_ref_scaffold(*itr);
					}
					else if (itr->right() < bundle.left())
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
				if (itr->ref_id() == first_ref_id_seen)
				{
					past_first_id = true;
				}
			}
			--itr;
		}
	}
	
	
	bundle_out = bundle;
	bundle_out.finalize();
	bundle_out.remove_hitless_scaffolds();
	return true;
}
