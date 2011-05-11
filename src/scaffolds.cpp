/*
 *  scaffolds.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/30/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <list>
#include <algorithm>
#include <numeric>
#include "common.h"
#include "scaffolds.h"

using namespace std;

bool AugmentedCuffOp::compatible(const AugmentedCuffOp& lhs,
								 const AugmentedCuffOp& rhs,
								 int overhang_tolerance)
{
	if (rhs.opcode == CUFF_INTRON)
	{
		if (lhs.opcode == CUFF_INTRON)
		{ 
			if (lhs != rhs /*&& !(lhs.genomic_length == rhs.genomic_length && 
							   abs(lhs.genomic_offset - rhs.genomic_offset) < 2)*/)
				return false;
		}
		else if (lhs.opcode == CUFF_UNKNOWN)
		{
            //int left_diff = abs(lhs.g_left() - rhs.g_left());
            //int right_diff = abs(lhs.g_right() - rhs.g_right());
			//if (left_diff + right_diff > max_frag_len)
            //    return false;
        }
        else
        {
            int l_match = match_length(lhs, rhs.g_left(), rhs.g_right());
            if (l_match > overhang_tolerance)
            {
                return false;
            }
        }
	}
	else if (rhs.opcode == CUFF_UNKNOWN)
	{
        int l_match = match_length(lhs, rhs.g_left(), rhs.g_right());
		if (l_match > max_frag_len)
			return false;
	}

	if (lhs.opcode == CUFF_INTRON)
	{
		if (rhs.opcode == CUFF_INTRON)
		{
			if (lhs != rhs /*&& !(lhs.genomic_length == rhs.genomic_length && 
							  abs(lhs.genomic_offset - rhs.genomic_offset) < 2)*/)
				return false;
		}
		else if (lhs.opcode == CUFF_UNKNOWN)
		{
            //int left_diff = abs(lhs.g_left() - rhs.g_left());
            //int right_diff = abs(lhs.g_right() - rhs.g_right());
            //if (left_diff + right_diff > max_frag_len)
            //    return false;
		}
		else 
        {
            int r_match = match_length(rhs, lhs.g_left(), lhs.g_right());
            if (r_match > overhang_tolerance)
            {
                return false;
            }
        }
	}	
	else if (lhs.opcode == CUFF_UNKNOWN)
	{
        int r_match = match_length(rhs, lhs.g_left(), lhs.g_right());
		if (r_match > max_frag_len)
			return false;
	}
	
	return true;
}

bool AugmentedCuffOp::g_left_lt(const AugmentedCuffOp& lhs,
                                const AugmentedCuffOp& rhs)
{
	return lhs.g_left() < rhs.g_left();
}

void disjoint_ops(vector<AugmentedCuffOp>& to_reduce)
{
	if (to_reduce.empty())
		return;
	
	vector<AugmentedCuffOp> reduced;
	reduced.push_back(to_reduce.front());
	for (int i = 1; i < (int)to_reduce.size(); ++i)
	{
		assert (to_reduce[i].opcode == to_reduce[i - 1].opcode);
		if (reduced.back().g_right() >= to_reduce[i].g_left())
		{
			int delta = to_reduce[i].g_right() - reduced.back().g_right();
			if (delta > 0)
				reduced.back().genomic_length += delta; 
		}
		else
		{
			reduced.push_back(to_reduce[i]);
		}
	}
	
	to_reduce = reduced;
}

// Adds open intervals not covered in the genomic coordinate covered by to_fill
// to the vector gaps.  DOES NOT CLEAR gaps.
void record_gaps(const vector<AugmentedCuffOp>& to_fill, 
				 vector<pair<int, int> >& gaps)
{
	for (size_t i = 1; i < to_fill.size(); ++i)
	{
		if (to_fill[i].g_left() - to_fill[i-1].g_right() > 0)
		{
			gaps.push_back(make_pair(to_fill[i-1].g_right(), to_fill[i].g_left()));
		}
	}	
}

// This function "fills" the gaps in to_fill with 
// AugmentedCuffOps from filler. The intersection of the gaps in both vectors
// remains as gaps in the modified to_fill.

// IMPORTANT: both vectors MUST be disjoint (see disjoint_matches) before calling
// this function
void AugmentedCuffOp::fill_interstices(vector<AugmentedCuffOp>& to_fill,
                                       const vector<AugmentedCuffOp>& filler,
                                       bool allow_flank_fill,
                                       bool allow_flanking_introns)
{
	vector<AugmentedCuffOp> filled = to_fill;
	vector<pair<int, int> > gaps;
	
	
	size_t j = 0;
	
	if (to_fill.empty())
	{
		to_fill = filler;
		sort(to_fill.begin(), to_fill.end(), g_left_lt);
		return;
	}
	
	// This first loop could scan from either end and bail on hitting
	// the first gap, but this is straightforward, and probably just as
	// fast in practice, since these arrays are generally tiny
	if (allow_flank_fill)
	{
		gaps.push_back(make_pair(0, to_fill.front().g_left()));
		record_gaps(to_fill, gaps);
		gaps.push_back(make_pair(to_fill.back().g_right(), INT_MAX));
	}
	else
	{
		record_gaps(to_fill, gaps);
	}
	
	size_t i = 0; 
	
	while (i < gaps.size())
	{
		pair<int, int>& gap = gaps[i];
		
		// a break in this loop will advance the gap index
		while (j < filler.size())
		{
			const AugmentedCuffOp& op = filler[j];
			
			if (op.g_left() == gap.first && op.g_right() == gap.second)
			{
                // CASE 1
				//  gap     [     )
				//  op      [     ) 
				filled.push_back(op);
				
				// advance both indexes (i is advanced after the break);
				++j;
                //fprintf (stderr, "CASE 1: advancing both indexes\n"); 
				break;
			}
			
			else if (op.g_right() <= gap.first)
			{
                // CASE 2
				//  gap   [     )
				//  op  [ ) 
				
                //fprintf (stderr, "CASE 2: skipping op %d:%d-%d due to gap %d-%d\n", op.opcode, op.g_left(), op.g_right(), gap.first, gap.second);
				// now just move to the next op, we can't add this one
			}
			else if (op.g_left() >= gap.second)
			{
                // CASE 3
				//  gap     [     )
				//  op            [ ) 
                //fprintf (stderr, "CASE 3: advancing gap from %d-%d due to op %d:%d-%d\n", gap.first, gap.second, op.opcode, op.g_left(), op.g_right());
				break; // not safe to add yet, we've gone beyond the current gap
				//  advance the gap index
			}
			else if (op.g_left() < gap.first && op.g_right() > gap.second)
			{
                // CASE 4
				//  gap     [     )
				//  op    [         ) 
				
				// create a new op to fill the gap
				AugmentedCuffOp gap_op(op.opcode, 
									   gap.first, 
									   gap.second - gap.first);
				assert (gap_op.genomic_length > 0);
				
				filled.push_back(gap_op);
                
                //fprintf (stderr, "CASE 4: advancing gap from %d-%d due to op %d:%d-%d\n", gap.first, gap.second, op.opcode, op.g_left(), op.g_right());
				// advance the gap index
				break;
			}
			else if (op.g_left() > gap.first && op.g_right() < gap.second)
			{
                // CASE 5
				//  gap     [     )
				//  op        [ ) 
				
				//  just add this op
				filled.push_back(op);
				//  advance the op index
                //fprintf (stderr, "CASE 5: adding %d:%d-%d, advancing op index\n", op.opcode, op.g_left(), op.g_right());
			}
			else if (op.g_right() >= gap.second && op.g_left() >= gap.first)
			{
                // CASE 6
				//  gap     [     )
				//  op           [ ) 
				
				//  create a new op from the left part of this one and add it
				AugmentedCuffOp gap_op(op.opcode, 
									   op.g_left(), 
									   gap.second - op.g_left());
				assert (gap_op.genomic_length > 0);
				filled.push_back(gap_op);
                
				break;
			}
			else if (op.g_left() <= gap.first && op.g_right() >= gap.first)
			{
                // CASE 7
				//  gap     [     )
				//  op     [ ) 
				
				//  create a new op from the right part of this one and add it
				AugmentedCuffOp gap_op(op.opcode, 
									   gap.first, 
									   op.g_right() - gap.first);
				assert (gap_op.genomic_length > 0);
				
				filled.push_back(gap_op);
				
				//  advance the op index
                //fprintf (stderr, "CASE 7: advancing op\n");
			}
			else
			{
				assert(false);
			}
            
			++j;
		}
		
		++i;
	}
    
	sort(filled.begin(), filled.end(), g_left_lt);
    if (!allow_flanking_introns)
    {
        for (size_t i = 0; i < filled.size(); ++i)
        {
            if (filled[i].opcode == CUFF_INTRON)
            {
                assert (i > 0);
                assert (i < filled.size() -1);
                assert (filled[i-1].opcode == CUFF_MATCH);
                assert (filled[i+1].opcode == CUFF_MATCH);
                assert (filled[i-1].g_right() == filled[i].g_left());
                assert (filled[i+1].g_left() == filled[i].g_right());
            }
        }
    }
	to_fill = filled;
}



// ops is assumed to be sorted
void AugmentedCuffOp::merge_ops(const vector<AugmentedCuffOp>& ops, 
                                vector<AugmentedCuffOp>& merged,
                                bool introns_overwrite_matches,
                                bool allow_flank_introns)
{	
#if DEBUG
	//assert(std::adjacent_find(ops.begin(), ops.end(), g_left_lt) == ops.end());
#endif
	
	if (ops.size() < 2)
	{
		merged = ops;
		return;
	}
	
	size_t g_max = 0;
	size_t g_min = 0xFFFFFFFF;
	
	vector<AugmentedCuffOp> matches;
	vector<AugmentedCuffOp> introns;
	
	vector<AugmentedCuffOp> unknowns;
	
	for (size_t i = 0; i < ops.size(); ++i)
	{
		//if (ops[i].opcode == CUFF_INTRON)
		//	fprintf (stderr, "[%d] %d, %d (%d)\n", i, ops[i].g_left(),ops[i].g_right(), ops[i].g_right() - ops[i].g_left() );
		assert (ops[i].g_left() < ops[i].g_right());
		
		if ((size_t)ops[i].g_left() < g_min)
			g_min = ops[i].g_left();
		if ((size_t)ops[i].g_right() > g_max)
			g_max = ops[i].g_right();

		switch(ops[i].opcode)
		{
			case CUFF_MATCH:
			{
				matches.push_back(ops[i]);
			}break;
                
			case CUFF_INTRON:
			{
				introns.push_back(ops[i]);
			}break;
                
			case CUFF_UNKNOWN:
			{
                
			}break;
                
			default:
				fprintf(stderr, "Unknown opcode, exiting\n");
				exit(1);
				break;
		}
	}
	
	int merged_length = (int)g_max - (int)g_min + 1;
	if (merged_length < 0)
    {
        fprintf(stderr, "Error: nonsense gene merge - negative length product\n");
		exit(1);
    }
    if (gaurd_assembly() && merged_length > (int)max_gene_length)
	{
		fprintf(stderr, "Error: nonsense gene merge - product > max_gene_length\n");
		exit(1);
	}
	
	unknowns.push_back(AugmentedCuffOp(CUFF_UNKNOWN, g_min, g_max - g_min + 1));
	
    vector<AugmentedCuffOp> pre_merge = matches;
	
	disjoint_ops(matches);
	disjoint_ops(introns);
	
	// below, the flank fill flag is set to preserve the invariant that 
	// all scaffolds must start with CUFF_MATCH
	if (introns_overwrite_matches)
	{
		merged = introns;
		fill_interstices(merged, matches, true, allow_flank_introns);
		vector<pair<int, int> > gaps;
		record_gaps(merged, gaps);
		if (!gaps.empty())
			fill_interstices(merged, unknowns, false, allow_flank_introns); 
	}
	else
	{
		merged = matches;
		fill_interstices(merged, introns, false, allow_flank_introns);
		vector<pair<int, int> > gaps;
		record_gaps(merged, gaps);
		if (!gaps.empty())
			fill_interstices(merged, unknowns, false, allow_flank_introns); 
	}

    if (!allow_flank_introns)
    {
        assert (merged.front().opcode == CUFF_MATCH);
        assert (merged.back().opcode == CUFF_MATCH);
    }
    for (size_t i = 1; i < merged.size(); ++i)
    {
		assert(merged[i-1].g_right() == merged[i].g_left());
        if (merged[i].opcode == CUFF_INTRON)
        {
            assert (merged[i-1].opcode == CUFF_MATCH);
        }
        else if (merged[i].opcode == CUFF_UNKNOWN)
        {
            assert (merged[i-1].opcode != CUFF_INTRON);
        }
    }
}
    

bool is_known(const AugmentedCuffOp& op)
{
	return op.opcode != CUFF_UNKNOWN;
}

// verifies that no matter how the merge goes, the result wont be an insanely 
// long gene.
bool check_merge_length(const vector<AugmentedCuffOp>& ops)
{
	size_t g_max = 0;
	size_t g_min = 0xFFFFFFFF;
	
	for (size_t i = 0; i < ops.size(); ++i)
	{
		//if (ops[i].opcode == CUFF_INTRON)
		//	fprintf (stderr, "[%d] %d, %d (%d)\n", i, ops[i].g_left(),ops[i].g_right(), ops[i].g_right() - ops[i].g_left() );
		assert (ops[i].g_left() < ops[i].g_right());
		
		if ((size_t)ops[i].g_left() < g_min)
			g_min = ops[i].g_left();
		if ((size_t)ops[i].g_right() > g_max)
			g_max = ops[i].g_right();
	}	
	int merged_length = (int)g_max - (int)g_min + 1;
	if (merged_length < 0 || merged_length > (int)max_gene_length)
	{
		return false;
			}
	return true;
		}

inline bool has_intron(const Scaffold& scaff)
{
			
	const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
	for (size_t j = 0; j < ops.size(); ++j)
	{
		if (ops[j].opcode == CUFF_INTRON)
			return true;
	}
	
	return false;
}

//inline bool has_intron(const Scaffold& scaff)
//{
//	
//	const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
//	for (size_t j = 0; j < ops.size(); ++j)
//	{
//		if (ops[j].opcode == CUFF_INTRON)
//			return true;
//	}
//	
//	return false;
//}

void Scaffold::extend_5(const Scaffold& other)
{
	assert(Scaffold::compatible(*this, other, ref_merge_overhang_tolerance));
	
	if (strand() == CUFF_FWD)
	{
		AugmentedCuffOp& my_first_op = _augmented_ops.front();
		const AugmentedCuffOp& other_first_op = other.augmented_ops().front();

		my_first_op.g_left(other_first_op.g_left());
	}
	else if (strand() == CUFF_REV)
	{
		AugmentedCuffOp& my_last_op = _augmented_ops.back();
		const AugmentedCuffOp& other_last_op = other.augmented_ops().back();
		
		my_last_op.g_right(other_last_op.g_right());
	}
	else 
	{
		assert(false);
	}
	
	int id_length = annotated_trans_id().length();
	if (id_length < 4 || _annotated_trans_id.substr(id_length-4)!="_ext")
		_annotated_trans_id += "_ext";
}

// Clip final 3' exon by given amount
void Scaffold::trim_3(int to_remove)
{
	if(strand() == CUFF_FWD)
	{
		AugmentedCuffOp& exon_3 = _augmented_ops.back();
		assert (exon_3.genomic_length > to_remove);
		exon_3.genomic_length -= to_remove;
	}
	else if(strand() == CUFF_REV)
	{
		AugmentedCuffOp& exon_3 = _augmented_ops.front();
		assert (exon_3.genomic_length > to_remove);
		exon_3.genomic_offset += to_remove;
		exon_3.genomic_length -= to_remove;

	}
	else
	{
		assert(false);
	}
}

// Clip final 3' exon by given amount
void Scaffold::extend_3(int to_add)
{
	if(strand() == CUFF_FWD)
	{
		AugmentedCuffOp& exon_3 = _augmented_ops.back();
		//assert (exon_3.genomic_length > to_remove);
		exon_3.genomic_length += to_add;
	}
	else if(strand() == CUFF_REV)
	{
		AugmentedCuffOp& exon_3 = _augmented_ops.front();
		//assert (exon_3.genomic_length > to_remove);
		exon_3.genomic_offset -= to_add;
		exon_3.genomic_length += to_add;
        
	}
	else
	{
		assert(false);
	}
}

void Scaffold::tile_with_scaffs(vector<Scaffold>& tile_scaffs, int max_len, int tile_offset) const
{
    int min_len = tile_offset;
	assert(min_len % tile_offset == 0 && max_len % tile_offset == 0);
	
	if (length() < max_len)
		return;
	
	size_t l = 0;
	size_t r = 0;
	int l_off = 0;
	int r_off = min_len; // One past end
	int curr_len = min_len;
	int remaining_len = length();
	
    const vector<AugmentedCuffOp>& orig_ops = augmented_ops();
    
	while(true)
	{
		while (l < orig_ops.size() && (orig_ops[l].opcode != CUFF_MATCH || l_off >= orig_ops[l].genomic_length))
		{
			if (orig_ops[l].opcode == CUFF_MATCH)
				l_off -= orig_ops[l].genomic_length;
            ++l;
			//if(++l == augmented_ops().size())
			//	assert(false);
		}
		while (r < orig_ops.size() && (orig_ops[r].opcode != CUFF_MATCH || r_off > orig_ops[r].genomic_length)) // Strictly > because r_off is one past
		{
			if (orig_ops[r].opcode == CUFF_MATCH)
				r_off -= orig_ops[r].genomic_length;
            ++r;
			//if(++r == augmented_ops().size())
			//	assert(false);
		}
		
		vector<AugmentedCuffOp> ops;
		
        //if (l >= orig_ops.size() && r >= orig_ops.size())
        //    break;
        
		if (l==r)
		{
            assert (l < orig_ops.size());
			ops.push_back(AugmentedCuffOp(CUFF_MATCH, orig_ops[l].g_left() + l_off, r_off - l_off));
		}
		else
		{
            assert(orig_ops.size());
            //assert(orig_ops[l].genomic_offset != 0);
            
			ops.push_back(AugmentedCuffOp(CUFF_MATCH, orig_ops[l].g_left() + l_off, orig_ops[l].genomic_length - l_off));
			assert(ops.back().g_right() > ops.back().g_left());
			for(size_t i = l+1; i < r; i++)
			{
				ops.push_back(orig_ops[i]);
			}
			if (r_off > 0 && r < orig_ops.size())
            {
                assert (r < orig_ops.size());
                //assert (orig_ops[r].genomic_offset != 0);
				ops.push_back(AugmentedCuffOp(CUFF_MATCH, orig_ops[r].g_left(), r_off));
            }
		}
		assert(ops.back().g_right() > ops.back().g_left());

        // genomic_offset actually could be zero - from an exon starting at coord
        // 1 in some chromosome of the ref.
//		foreach(const AugmentedCuffOp& op, ops)
//        {
//            assert (op.genomic_offset != 0);
//        }
        
		tile_scaffs.push_back(Scaffold(this->ref_id(), this->strand(), ops, true)); 
		assert(tile_scaffs.back().length() == curr_len );
		assert(tile_scaffs.back().left() >= left() &&
			   tile_scaffs.back().right() <= right());
		
		if (l==0 && l_off == 0 && curr_len < max_len) // On left end, not yet full length
		{
			curr_len += tile_offset;
			r_off += tile_offset;
		}
		else if(remaining_len - tile_offset < max_len) // On the right end of transcript, decreasing in length
		{
			curr_len += augmented_ops().back().genomic_length - r_off - tile_offset; 
			r_off = augmented_ops().back().genomic_length;
			l_off += tile_offset;
			remaining_len -= tile_offset;
			if (curr_len < min_len)
				return;
		}
		else // In the middle of the transcript, full length
		{
			l_off += tile_offset;
			r_off += tile_offset;
			remaining_len -= tile_offset;
		}
			
	}
	
}

bool Scaffold::sub_scaffold(Scaffold& sub_scaff, int g_left, int match_length) const	
{
	size_t i;
	for(i = 0; i < augmented_ops().size(); ++i)
	{
		if (augmented_ops()[i].g_right() > g_left)
			break;
	}
	const AugmentedCuffOp& l_op = augmented_ops()[i++];
	
	assert(l_op.opcode == CUFF_MATCH);
	
	vector<AugmentedCuffOp> sub_ops;
	sub_ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, min(match_length, l_op.g_right()-g_left)));					 
	int len_so_far = sub_ops.back().genomic_length;
						 
	while(len_so_far < match_length && i < augmented_ops().size())
	{
		const AugmentedCuffOp& op = augmented_ops()[i++];
		if(op.opcode==CUFF_MATCH)
		{
			sub_ops.push_back(AugmentedCuffOp(CUFF_MATCH, op.g_left(), min(match_length-len_so_far, op.genomic_length)));
			len_so_far += sub_ops.back().genomic_length;
			if (len_so_far >= match_length)
				break;
		}
		else 
		{
			sub_ops.push_back(op);
		}
	}
	
	//assert(len_so_far == match_length);
	sub_scaff = Scaffold(this->ref_id(), this->strand(), sub_ops, true);
	return (len_so_far == match_length);
}

void Scaffold::merge(const Scaffold& lhs, 
					 const Scaffold& rhs, 
					 Scaffold& merged,
					 bool introns_overwrite_matches)
{
	OpList ops;
	
	assert (merged.ref_id() == 0);
	assert (lhs.ref_id() == rhs.ref_id());
	
	ops.insert(ops.end(), lhs._augmented_ops.begin(), lhs._augmented_ops.end());
	ops.insert(ops.end(), rhs._augmented_ops.begin(), rhs._augmented_ops.end());
	
	sort(ops.begin(), ops.end(), AugmentedCuffOp::g_left_lt);
	
	AugmentedCuffOp::merge_ops(ops, merged._augmented_ops, introns_overwrite_matches);
	
	assert (ops.empty() || !(merged._augmented_ops.empty()));
	
	merged._ref_id = lhs.ref_id();
	merged._mates_in_scaff.insert(merged._mates_in_scaff.end(),
								  lhs._mates_in_scaff.begin(), 
								  lhs._mates_in_scaff.end());
	merged._mates_in_scaff.insert(merged._mates_in_scaff.end(),
								  rhs._mates_in_scaff.begin(), 
								  rhs._mates_in_scaff.end());
	
	sort(merged._mates_in_scaff.begin(),
		 merged._mates_in_scaff.end());
	vector<const MateHit*>::iterator new_end = unique(merged._mates_in_scaff.begin(), 
													  merged._mates_in_scaff.end());
	merged._mates_in_scaff.erase(new_end, merged._mates_in_scaff.end());
	
	
	int r_check = merged.left();
	for (size_t i = 0; i < merged._augmented_ops.size(); ++i)
		r_check += merged._augmented_ops[i].genomic_length;
	
	assert(r_check == merged.right());
	
	if (lhs.strand() != CUFF_STRAND_UNKNOWN)
		merged._strand = lhs.strand();
	if (rhs.strand() != CUFF_STRAND_UNKNOWN)
		merged._strand = rhs.strand();
	
	merged._has_intron = has_intron(merged);
	assert(merged._strand != CUFF_STRAND_UNKNOWN || !merged.has_intron() );
	assert(!merged.is_ref());
    
    if (library_type == "transfrags")
    {
        double avg_fpkm = lhs.fpkm() + rhs.fpkm();
        avg_fpkm /= 2;
        merged.fpkm(avg_fpkm);
    }
}



void Scaffold::merge(const vector<Scaffold>& s, 
					 Scaffold& merged, 
					 bool introns_overwrite_matches)
{
	OpList ops;
	
	CuffStrand strand = CUFF_STRAND_UNKNOWN;
	
	for (size_t i = 0; i < s.size(); ++i)
	{
		ops.insert(ops.end(), s[i]._augmented_ops.begin(), s[i]._augmented_ops.end());
		merged._mates_in_scaff.insert(merged._mates_in_scaff.end(),
									  s[i]._mates_in_scaff.begin(), 
									  s[i]._mates_in_scaff.end());

		if (s[i].strand() != CUFF_STRAND_UNKNOWN)
		{
			//assert (strand == CUFF_STRAND_UNKNOWN || strand == s[i].strand());
			strand = s[i].strand();
		}
	}
	
	sort(merged._mates_in_scaff.begin(),merged._mates_in_scaff.end());
	vector<const MateHit*>::iterator new_end = unique(merged._mates_in_scaff.begin(), 
													  merged._mates_in_scaff.end());
	merged._mates_in_scaff.erase(new_end, merged._mates_in_scaff.end());
	
	sort(ops.begin(), ops.end(), AugmentedCuffOp::g_left_lt);
	//merged._contigs.push_back(CuffAlign(ops.front().g_left(), vector<CuffOp>()));
	
	if (ops.empty())
		return;
	
	AugmentedCuffOp::merge_ops(ops, merged._augmented_ops, introns_overwrite_matches);
	
	
	//assert (ops.empty() || !(merged._augmented_ops.empty()));
#ifdef DEBUG
	if (merged._augmented_ops.empty() || 
		merged._augmented_ops.front().opcode != CUFF_MATCH ||
		merged._augmented_ops.back().opcode != CUFF_MATCH)
	{
		AugmentedCuffOp::merge_ops(ops, merged._augmented_ops, introns_overwrite_matches);
	}
#endif
	
	merged._ref_id = s.front().ref_id();
	merged._strand = strand;
		
	int r_check = merged.left();
	for (size_t i = 0; i < merged._augmented_ops.size(); ++i)
		r_check += merged._augmented_ops[i].genomic_length;
	
#ifdef DEBUG
	if (r_check != merged.right())
	{
		AugmentedCuffOp::merge_ops(ops, merged._augmented_ops, introns_overwrite_matches);
	}
#endif
	
	assert (r_check == merged.right());
	merged._has_intron = has_intron(merged);

	assert(merged._strand != CUFF_STRAND_UNKNOWN || !merged.has_strand_support());
    
    if (library_type == "transfrags")
    {
        double avg_fpkm = 0.0;
        foreach (const Scaffold& sc, s)
        {
            avg_fpkm += sc.fpkm();
        }
        avg_fpkm /= s.size();
        merged.fpkm(avg_fpkm);
    }
}

void Scaffold::fill_gaps(int filled_gap_size)
{
	OpList ops;
	
	const vector<AugmentedCuffOp>& orig_ops = augmented_ops();
	for (size_t i = 0; i < orig_ops.size(); ++i)
	{
		if (orig_ops[i].opcode == CUFF_UNKNOWN && 
			orig_ops[i].genomic_length < filled_gap_size)
		{
			ops.push_back(AugmentedCuffOp(CUFF_MATCH, 
										  orig_ops[i].genomic_offset, 
										  orig_ops[i].genomic_length)); 
		}
		else
		{
			ops.push_back(orig_ops[i]);
		}
	}
	
	sort(ops.begin(), ops.end(), AugmentedCuffOp::g_left_lt);
	
	AugmentedCuffOp::merge_ops(ops, _augmented_ops, false);
	_has_intron = has_intron(*this);
	assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);

}

void Scaffold::fill_gaps(const vector<AugmentedCuffOp>& filler)
{
	OpList ops;
	
	const vector<AugmentedCuffOp>& orig_ops = augmented_ops();
    
	const vector<AugmentedCuffOp> og_ops = augmented_ops();
	
    vector<AugmentedCuffOp> unknowns;
	
    size_t g_max = 0;
	size_t g_min = 0xFFFFFFFF;
    
    vector<AugmentedCuffOp> tmp_filler = filler;
    
    foreach(const AugmentedCuffOp& op, orig_ops)
    {
		assert (op.g_left() < op.g_right());
		
		if ((size_t)op.g_left() < g_min)
			g_min = op.g_left();
		if ((size_t)op.g_right() > g_max)
			g_max = op.g_right();
	}

	tmp_filler.push_back(AugmentedCuffOp(CUFF_UNKNOWN, g_min, g_max - g_min + 1));
    sort(tmp_filler.begin(), tmp_filler.end(), AugmentedCuffOp::g_left_lt);
    
    vector<AugmentedCuffOp> padded_filler;
    AugmentedCuffOp::merge_ops(tmp_filler, padded_filler, false);
    	
    vector<AugmentedCuffOp> overlapping;
    foreach (const AugmentedCuffOp& op, padded_filler)
    {
		//if (left() <= op.g_left() && right() >= op.g_right()
		if(::overlap_in_genome(op.g_left(),op.g_right(), left(), right())
            && (op.opcode != CUFF_UNKNOWN || !overlapping.empty()))
        {
            overlapping.push_back(op);
        }
    }
    
    overlapping.insert(overlapping.end(), _augmented_ops.begin(),_augmented_ops.end()); 
    sort(overlapping.begin(), overlapping.end(), AugmentedCuffOp::g_left_lt);
	
	// Trim first in last in case they hang over the scaffold boundaries
	if (overlapping.front().g_left() < left())
	{
		AugmentedCuffOp& first_op = overlapping.front();
		first_op.genomic_length -= (left() - first_op.g_left());
		first_op.genomic_offset = left();
	}
	if (overlapping.back().g_right() > right())
	{
		AugmentedCuffOp& last_op = overlapping.back();
		last_op.genomic_length -= (last_op.g_right() - right());
	}
	if (overlapping.back().opcode == CUFF_INTRON && overlapping.back().genomic_length <= bowtie_overhang_tolerance)
	{
		overlapping.pop_back();
	}
	if (overlapping.front().opcode == CUFF_INTRON && overlapping.front().genomic_length <= bowtie_overhang_tolerance)
	{
		overlapping.erase(overlapping.begin(), overlapping.begin()+1);
	}
	
    // we don't want either the terminal ops in the filler to be unknowns,
    // because they'll just propogate to the scaffold we're trying to fill.
    if (!overlapping.empty() && overlapping.back().opcode == CUFF_UNKNOWN)
        overlapping.pop_back();
    
    if (overlapping.empty())
        return;
    
	AugmentedCuffOp::merge_ops(overlapping, _augmented_ops, true);
	_has_intron = has_intron(*this);
	assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);
}

bool Scaffold::overlap_in_genome(const Scaffold& lhs, 
								 const Scaffold& rhs, 
								 int overlap_radius)
{
	int ll = lhs.left() - overlap_radius;
	int rr = rhs.right() + overlap_radius;
	int lr = lhs.right() + overlap_radius;
	int rl = rhs.left() - overlap_radius;
	
	if (ll >= rl && ll < rr)
		return true;
	if (lr >rl && lr < rr)
		return true;
	if (rl >= ll && rl < lr)
		return true;
	if (rr > ll && rr < lr)
		return true;
	
	return false;
}

bool intron_op(const AugmentedCuffOp& op)
{
	return op.opcode == CUFF_INTRON;
}

bool Scaffold::strand_agree(const Scaffold& lhs, 
							const Scaffold& rhs)
{
	bool strand = (lhs.strand() == CUFF_STRAND_UNKNOWN || 
				   rhs.strand() == CUFF_STRAND_UNKNOWN ||
				   lhs.strand() == rhs.strand());
	return strand;
}

bool Scaffold::exons_overlap(const Scaffold& lhs, 
							 const Scaffold& rhs)
{
	const vector<AugmentedCuffOp>& lhs_ops = lhs.augmented_ops();
	const vector<AugmentedCuffOp>& rhs_ops = rhs.augmented_ops();

	for (size_t l = 0; l < lhs_ops.size(); l++)
	{
		if (lhs_ops[l].opcode != CUFF_MATCH)
			continue;
		
		for (size_t r = 0; r < rhs_ops.size(); r++)
		{
			if (rhs_ops[r].opcode == CUFF_MATCH &&
				AugmentedCuffOp::overlap_in_genome(lhs_ops[l], rhs_ops[r]))
			{
				return true;
			}
		}
	}
	return false;
}

bool Scaffold::compatible(const Scaffold& lhs, 
						  const Scaffold& rhs,
						  int overhang_tolerance)
{	
	if (!strand_agree(lhs, rhs))
		return false;
	
	if (lhs.left() <= rhs.left())
	{
		if (overlap_in_genome(lhs, rhs, olap_radius))
		{
			// check compatibility
			if (!compatible_contigs(lhs, rhs, overhang_tolerance))
				return false;
		}
	}
	else if (rhs.left() < lhs.left())
	{
		if (overlap_in_genome(rhs, lhs, olap_radius))
		{
			// check compatibility
			if (!compatible_contigs(rhs, lhs, overhang_tolerance))
				return false;
		}
	}

	
	return true;
}

bool Scaffold::compatible_contigs(const Scaffold& lhs, 
								  const Scaffold& rhs,
								  int overhang_tolerance)
{
	const vector<AugmentedCuffOp>& l_aug = lhs._augmented_ops;
	const vector<AugmentedCuffOp>& r_aug = rhs._augmented_ops;
	
	size_t curr_l_op = 0;
	size_t curr_r_op = 0;
	
	while (curr_l_op != l_aug.size() &&
		   curr_r_op != r_aug.size())
	{
		const AugmentedCuffOp& l_op = l_aug[curr_l_op];
		const AugmentedCuffOp& r_op = r_aug[curr_r_op];
		
		if (l_op.g_left() <= r_op.g_left())
		{
			if (AugmentedCuffOp::overlap_in_genome(l_op, r_op))
			{
				// check compatibility
				if (!AugmentedCuffOp::compatible(l_op, r_op, overhang_tolerance))
					return false;
//				if (l_op.opcode == CUFF_UNKNOWN && 
//					r_op.opcode == CUFF_MATCH)
//				{
//					//int m_len = AugmentedCuffOp::match_length(r_op, l_op.g_left(), l_op.g_right());
//					if (l_op.properly_contains(r_op))
//						return false;
//				}
			}
			if (l_op.g_right() < r_op.g_right())
				++curr_l_op;
			else if (r_op.g_right() < l_op.g_right())
				++curr_r_op;
			else // Indentical boundaries, advance both
			{
				++curr_l_op;
				++curr_r_op;
			}
		}
		else if (r_op.g_left() < l_op.g_left())
		{
			if (AugmentedCuffOp::overlap_in_genome(r_op, l_op))
			{
				// check compatibility
				if (!AugmentedCuffOp::compatible(r_op, l_op, overhang_tolerance))
						return false;
//				if (r_op.opcode == CUFF_UNKNOWN && 
//					l_op.opcode == CUFF_MATCH)
//				{
//					//int m_len = AugmentedCuffOp::match_length(l_op, r_op.g_left(), r_op.g_right());
//					if (r_op.properly_contains(l_op))
//						return false;
//				}
			}
			if (r_op.g_right() < l_op.g_right())
				++curr_r_op;
			else if (l_op.g_right() < r_op.g_right())
				++curr_l_op;
			else // Indentical boundaries, advance both
			{
				++curr_l_op;
				++curr_r_op;
			}
		}
		
	}
	
	return true;
	
}

bool overlap_in_genome(int ll, int lr, int rl, int rr)
{
	if (ll >= rl && ll < rr)
		return true;
	if (lr > rl && lr < rr)
		return true;
	if (rl >= ll && rl < lr)
		return true;
	if (rr > ll && rr < lr)
		return true;
	return false;
}

double Scaffold::score_overlap(const AugmentedCuffOp& op, 
							   int inner_left_edge, 
							   int inner_right_edge)
{
	assert(false);
	return 0; 
	/*
     assert (op.opcode == CUFF_MATCH);
     int mlen = AugmentedCuffOp::match_length(op, 
     inner_left_edge, 
     inner_right_edge);
     if (mlen > inner_dist_mean)
     {
     double c = cdf(inner_dist_norm, mlen) - 0.5;
     //	if (c == 1.0)
     //		return -1000.0;
     double weight = log(1 - c);
     assert (weight <= 0.0);
     return weight;
     }
     else
     {
     return 0;
     }
	 */
}

double Scaffold::min_score(const Scaffold& contig, 
						   const vector<pair<int, int> >& inners)
{
	double score = 0.0;
	
	size_t curr_inner = 0;
	size_t curr_op = 0;
	
	while (curr_op != contig._augmented_ops.size() &&
		   curr_inner != inners.size())
	{
		const pair<int, int>   inner = inners[curr_inner];
		const AugmentedCuffOp&  op = contig._augmented_ops[curr_op];
		
		if (inner.first <= op.g_left())
		{
			if (op.opcode == CUFF_MATCH &&
				::overlap_in_genome(inner.first, 
									inner.second, 
									op.g_left(), 
									op.g_right()))
			{
				score = min(score, score_overlap(op, inner.first, inner.second));
			}
			if (inner.second < op.g_right())
				++curr_inner;
			else if (op.g_right() < inner.second)
				++curr_op;
			else // Indentical boundaries, advance both
			{
				++curr_inner;
				++curr_op;
			}
		}
		else if (op.g_left() < inner.first)
		{
			if (op.opcode == CUFF_MATCH &&
				::overlap_in_genome(op.g_left(), 
									op.g_right(), 
									inner.first, 
									inner.second))
			{
				score = min(score, score_overlap(op, inner.first, inner.second));
			}
			if (op.g_right() < inner.second)
				++curr_op;
			else if (inner.second < op.g_right())
				++curr_inner;
			else // Indentical boundaries, advance both
			{
				++curr_inner;
				++curr_op;
			}
		}
	}
	return score;
}

double Scaffold::score_contigs(const Scaffold& contig, 
							   const vector<pair<int, int> >& inners)
{
	double score = 0.0;
	
	size_t curr_inner = 0;
	size_t curr_op = 0;
	
	while (curr_op != contig._augmented_ops.size() &&
		   curr_inner != inners.size())
	{
		const pair<int, int>   inner = inners[curr_inner];
		const AugmentedCuffOp&  op = contig._augmented_ops[curr_op];
		
		if (inner.first <= op.g_left())
		{
			if (op.opcode == CUFF_MATCH &&
				::overlap_in_genome(inner.first, 
								  inner.second, 
								  op.g_left(), 
								  op.g_right()))
			{
				score += score_overlap(op, inner.first, inner.second);
			}
			if (inner.second < op.g_right())
				++curr_inner;
			else if (op.g_right() < inner.second)
				++curr_op;
			else // Indentical boundaries, advance both
			{
				++curr_inner;
				++curr_op;
			}
		}
		else if (op.g_left() < inner.first)
		{
			if (op.opcode == CUFF_MATCH &&
				::overlap_in_genome(op.g_left(), 
								  op.g_right(), 
								  inner.first, 
								  inner.second))
			{
				score += score_overlap(op, inner.first, inner.second);
			}
			if (op.g_right() < inner.second)
				++curr_op;
			else if (inner.second < op.g_right())
				++curr_inner;
			else // Indentical boundaries, advance both
			{
				++curr_inner;
				++curr_op;
			}
		}
	}
	return score;
}

double Scaffold::score_merge(const Scaffold& lhs, const Scaffold& rhs)
{
	double score = 0.0;
		
	vector<pair<int, int> > l_inners, r_inners;
	
	vector<const MateHit*> lhs_but_not_rhs;
	vector<const MateHit*> rhs_but_not_lhs;
	
	set_difference(lhs.mate_hits().begin(),
				   lhs.mate_hits().end(),
				   rhs.mate_hits().begin(),
				   rhs.mate_hits().end(),
				   back_inserter(lhs_but_not_rhs));
	
	set_difference(rhs.mate_hits().begin(),
				   rhs.mate_hits().end(),
				   lhs.mate_hits().begin(),
				   lhs.mate_hits().end(),
				   back_inserter(rhs_but_not_lhs));
	
	for (vector<const MateHit*>::iterator i = lhs_but_not_rhs.begin();
		 i != lhs_but_not_rhs.end();
		 ++i)
	{
		pair<int,int> p = (*i)->genomic_inner_span();
		if (p.first != -1 && p.second != -1)
			l_inners.push_back(p);
	}
	
	for (vector<const MateHit*>::iterator i = rhs_but_not_lhs.begin();
		 i != rhs_but_not_lhs.end();
		 ++i)
	{
		pair<int,int> p = (*i)->genomic_inner_span();
		if (p.first != -1 && p.second != -1)
			r_inners.push_back(p);
	}
	
	score += Scaffold::score_contigs(lhs, r_inners);
	score += Scaffold::score_contigs(rhs, l_inners);
		
	return score;
}

double Scaffold::score() const
{	
	double score = 0.0;
	
	vector<pair<int, int> > inners;
	
	for (vector<const MateHit*>::const_iterator i = mate_hits().begin();
		 i != mate_hits().end();
		 ++i)
	{
		pair<int,int> p = (*i)->genomic_inner_span();
		if (p.first != -1 && p.second != -1)
			inners.push_back(p);
	}

	double contig_score = Scaffold::score_contigs(*this, inners);
	score += contig_score;

	return score;	
}

double Scaffold::worst_mate_score() const
{
	double min_score = 0.0;
	
	vector<pair<int, int> > inners;
	
	for (vector<const MateHit*>::const_iterator i = mate_hits().begin();
		 i != mate_hits().end();
		 ++i)
	{
		pair<int,int> p = (*i)->genomic_inner_span();
		if (p.first != -1 && p.second != -1)
			inners.push_back(p);
	}
	
	//return overlap_in_genome(lhs._contigs[0], rhs._contigs[0]);
	
	
	min_score = Scaffold::min_score(*this, inners);
	
	return min_score;	
}


int Scaffold::genomic_to_transcript_coord(int g_coord) const
{
	int s_coord = 0;
	size_t curr_op = 0;
	const AugmentedCuffOp* op=NULL;

	while (curr_op != _augmented_ops.size())
	{
		op = &_augmented_ops[curr_op];
		
		if (op->g_right() > g_coord)
			break;
		if (op->opcode == CUFF_MATCH)
			s_coord += op->genomic_length;
		
		++curr_op;
	}
	
	int remainder = g_coord - (*op).g_left();
	
	if (remainder <= bowtie_overhang_tolerance || op->opcode == CUFF_MATCH)
		s_coord += remainder; 
	else
		s_coord -= (op->genomic_length-remainder);
	
	if(strand()==CUFF_REV)
		s_coord = length() - 1 - s_coord;
	
	return s_coord;
}

// start and end (first and second) are the first and final coordinates of the span
// Does not handle overhangs (can return negative values)
pair <int,int> Scaffold::genomic_to_transcript_span(pair<int,int> g_span) const
{
	int s_start;
	int s_end;
	
	int s_coord = 0;
	size_t curr_op = 0;
	const AugmentedCuffOp* op = NULL;
	// First, find start
	while (curr_op != _augmented_ops.size())
	{
		op = &_augmented_ops[curr_op];
		
		if (op->g_right() > g_span.first)
			break;
		if (op->opcode == CUFF_MATCH)
			s_coord += op->genomic_length;
		
		++curr_op;
	}

	int remainder = g_span.first - (*op).g_left();
	
	if (remainder <= bowtie_overhang_tolerance || op->opcode == CUFF_MATCH)
		s_start = s_coord + remainder; 
	else
		s_start = s_coord - (op->genomic_length-remainder);

	
	// Now, find end
	while (curr_op != _augmented_ops.size())
	{
		op = &_augmented_ops[curr_op];
		
		if (op->g_right() > g_span.second)
			break;
		if (op->opcode == CUFF_MATCH)
			s_coord += op->genomic_length;
		++curr_op;
	}
	
	remainder = g_span.second - op->g_left();
	
	if (remainder < bowtie_overhang_tolerance || op->opcode == CUFF_MATCH)
		s_end = s_coord + remainder; 
	else
		s_end = s_coord - (op->genomic_length-remainder);
	
	if(strand()==CUFF_REV)
	{
		int scaff_len = length();
		s_start = scaff_len - 1 - s_start;
		s_end = scaff_len - 1 - s_end;
		swap(s_start, s_end);
	}
	
	return make_pair(s_start, s_end);
}

// Returns true only when both the start and end are found (ie, frag_len is known), which
// can only happen if the read is paired.  Returned values that are equal to the trans_len 
// should be ignored, as they are invalid.
// We can't use EmpDist unless it is unpaired since we call this function in inspect_bundle
// End is inclusive
// Returned indices are oriented with the transript!
bool Scaffold::map_frag(const MateHit& hit, int& start, int& end, int& frag_len) const
{
	
	int trans_len = length();
    
//    if (Scaffold(hit).augmented_ops() == augmented_ops())
//    {
//        int a = 4;
//    }
    
	start = trans_len;
	end = trans_len;
	
    Scaffold h(hit);
    if(!(contains(h) && Scaffold::compatible(*this, h)))
       return false;
    
    if (hit.read_group_props()->complete_fragments())
    {
        pair<int,int> g_span = make_pair(hit.left(), hit.right() - 1);
		pair<int,int> t_span = genomic_to_transcript_span(g_span);
		start = t_span.first;
		end = t_span.second;
		frag_len = abs(end-start)+1;
    }
	else if (hit.is_pair())
	{
		pair<int,int> g_span = hit.genomic_outer_span();
		pair<int,int> t_span = genomic_to_transcript_span(g_span);
		start = t_span.first;
		end = t_span.second;
		frag_len = abs(end-start)+1;		
	}
	else if (hit.read_group_props()->mate_strand_mapping()==FF)
	{
		shared_ptr<const EmpDist> frag_len_dist = hit.read_group_props()->frag_len_dist();
		frag_len = min(frag_len_dist->mode(), trans_len);
	}
	else
	{
		shared_ptr<const EmpDist> frag_len_dist = hit.read_group_props()->frag_len_dist();

		if ((hit.left_alignment()->antisense_align() && strand() != CUFF_REV) 
			|| (!hit.left_alignment()->antisense_align() && strand() == CUFF_REV))
		{
			int g_end  = (strand()!=CUFF_REV) ? hit.right()-1:hit.left();
			end = genomic_to_transcript_coord(g_end);
			frag_len = min(frag_len_dist->mode(), end);
		}
		else
		{
			int g_start = (strand()!=CUFF_REV) ? hit.left():hit.right()-1;
			start = genomic_to_transcript_coord(g_start);
			if (start == trans_len) // Overhang
				frag_len = min(frag_len_dist->mode(), trans_len);
			else
				frag_len = min(frag_len_dist->mode(), trans_len-start);
		}
	}

    if (start <= 0 || start > trans_len)
		start = trans_len;
	if (end <= 0 || end > trans_len)
		end = trans_len;
		
	return (start != trans_len && end != trans_len);
}

int Scaffold::match_length(int left, int right) const
{
	int len = 0;

	
	size_t curr_op = 0;
	
	while (curr_op != _augmented_ops.size())
	{
		const AugmentedCuffOp&  op = _augmented_ops[curr_op];
		
		if (op.opcode == CUFF_MATCH &&
			::overlap_in_genome(left, 
								right, 
								op.g_left(), 
								op.g_right()))
		{
			len += AugmentedCuffOp::match_length(op, left, right);
		}
		if (op.g_left() >= right)
			break;
		++curr_op;
	}
	
	return len;
}


void Scaffold::clear_hits()
{
    _mates_in_scaff.clear();
	vector<const MateHit*>(_mates_in_scaff).swap(_mates_in_scaff);
    //_mates_in_scaff.clear();
}

bool Scaffold::add_hit(const MateHit* hit)
{
	Scaffold hs(*hit);
	if (contains(hs) &&
		Scaffold::compatible(*this, hs))
	{
		if (!binary_search(_mates_in_scaff.begin(),
						  _mates_in_scaff.end(),
						  hit))
		{
			_mates_in_scaff.push_back(hit);
		}
		return true;
	}

	return false;
}


void Scaffold::get_complete_subscaffolds(vector<Scaffold>& complete)
{
	if (!has_unknown())
	{
		complete.push_back(*this);
	}
	else
	{
		int last_unknown = -1;
        int leftmost_known_op = -1;
        int rightmost_known_op = -1;
		for (size_t i = 0; i < _augmented_ops.size(); ++i)
		{
			const AugmentedCuffOp& op = _augmented_ops[i]; 
            if (op.opcode != CUFF_UNKNOWN)
            {
                if (leftmost_known_op == -1)
                {
                    leftmost_known_op = i;
                    assert (_augmented_ops[leftmost_known_op].opcode != CUFF_INTRON);
                }
                rightmost_known_op = i;
            }
			if (op.opcode == CUFF_UNKNOWN || i == _augmented_ops.size() - 1)
			{
				int left_known;
				int right_known;
				
				if (i == _augmented_ops.size() - 1)
					right_known = right();
				else
					right_known = op.g_left() - 1;
				
				if (last_unknown == -1)
					left_known = left();
				else
					left_known = _augmented_ops[last_unknown].g_right();
				
                //fprintf(stderr, "excluding unknown region between %d-%d\n", left_known, right_known);
                if (leftmost_known_op != -1 && rightmost_known_op != -1)
                {
                    
                    vector<AugmentedCuffOp> known_ops;
                    known_ops.insert(known_ops.end(), 
                                     &_augmented_ops[leftmost_known_op],
                                     &_augmented_ops[rightmost_known_op] + 1);
                    Scaffold known(ref_id(), strand(), known_ops);
                    
                    for (vector<const MateHit*>::iterator mitr = _mates_in_scaff.begin();
                         mitr != _mates_in_scaff.end();
                         ++mitr)
                    {

                        const MateHit& m = **mitr;
                        if (left_known <= m.left() && m.right() <= right_known)
                        {
                            known.add_hit(&m);
                        }
                    }
                    
                   // if (known.has_intron())
                        known.strand(_strand);
                    
                    //assert (!known.mate_hits().empty());
                    //assert(!known.has_intron()|| known.strand() != CUFF_STRAND_UNKNOWN);
                    
                   // const vector<const MateHit*>& hits = known.mate_hits();
                   // bool contains_spliced_hit = false;
                   // foreach (const MateHit* h, hits)
                   // {
                   //     const ReadHit* left = h->left_alignment();
                   //     const ReadHit* right = h->right_alignment();
                   //     if (left && left->contains_splice())
                   //     {
                   //         contains_spliced_hit = true;
                   //         break;
                   //     }
                   //     if (right && right->contains_splice())
                   //     {
                   //         contains_spliced_hit = true;
                   //         break;
                   //     }
				   // }
                    
                   // if (!hits.empty() && !contains_spliced_hit && !known.is_ref())
                   // {
                   //     known.strand(CUFF_STRAND_UNKNOWN);
                   // }
                    
                    complete.push_back(known);
                }
                    
				last_unknown = i;
                leftmost_known_op = -1;
			}
		}
	}
}

double Scaffold::internal_exon_coverage() const
{
    // First check if there are internal exons
    if (augmented_ops().size() < 5) 
        return 0.0;
    
    int left = augmented_ops()[2].g_left();
    int right = augmented_ops()[augmented_ops().size() - 3].g_right();
    vector<bool> covered(right-left, 0);
    foreach(const MateHit* h, mate_hits())
    {
        if (::overlap_in_genome(h->left(),h->right(), left, right))
        {
            for (int i = max(h->left()-left, 0); i < min(h->right()-left, right-left); ++i)
            {
                assert(i < covered.size());
                covered[i] = 1;
            }
        }
    }
    double percent_covered = accumulate(covered.begin(),covered.end(),0.0)/(right-left);
    return percent_covered;
}

bool Scaffold::has_strand_support(vector<shared_ptr<Scaffold> >* ref_scaffs) const
{
	if (strand() == CUFF_STRAND_UNKNOWN)
		return false;
	
	// FIXME:  This is not true for non-canonical splicing
	if (has_intron())
		return true;
	
	foreach (const MateHit* h, mate_hits())
	{
		if (h->strand() == strand())
			return true;
		assert(h->strand() == CUFF_STRAND_UNKNOWN && !h->contains_splice());
	}
	
	if (ref_scaffs == NULL)
		return false;
	
	foreach (shared_ptr<Scaffold const> ref_scaff, *ref_scaffs)
	{
		if (ref_scaff->strand() == strand() && exons_overlap(*this, *ref_scaff))
			return true;
	}
		
	return false;
}

bool Scaffold::has_struct_support(set<AugmentedCuffOp>& hit_introns) const
{
    if(augmented_ops().size() == 1)
        return mate_hits().size() > 0;

    if(augmented_ops().size() == 3)
        return hits_support_introns(hit_introns);
    
    return (hits_support_introns(hit_introns) && internal_exon_coverage() == 1.0);
}


bool Scaffold::hits_support_introns() const
{
    set<AugmentedCuffOp> hit_introns;
    set<AugmentedCuffOp> scaffold_introns;
    foreach(const MateHit* h, _mates_in_scaff)
    {
        Scaffold s(*h);
        foreach (AugmentedCuffOp a, s.augmented_ops())
        {
            if (a.opcode == CUFF_INTRON)
            {
                hit_introns.insert(a);
            }
        }
    }
    foreach (AugmentedCuffOp a, _augmented_ops)
    {
        if (a.opcode == CUFF_INTRON)
        {
            scaffold_introns.insert(a);
        }
    }
    
    if (hit_introns != scaffold_introns)
    {
        fprintf(stderr, "********************\n");
        foreach(const AugmentedCuffOp& a, hit_introns)
        {
            fprintf(stderr, "%d - %d\n", a.g_left(), a.g_right());
        }
        
        fprintf(stderr, "####################\n");    
        foreach(const AugmentedCuffOp& a, scaffold_introns)
        {
            fprintf(stderr, "%d - %d\n", a.g_left(), a.g_right());
        }
    }
    
    return hit_introns == scaffold_introns;
}

bool Scaffold::hits_support_introns(set<AugmentedCuffOp>& hit_introns) const
{
    set<AugmentedCuffOp> scaffold_introns;

    foreach (AugmentedCuffOp a, _augmented_ops)
    {
        if (a.opcode == CUFF_INTRON)
        {
            scaffold_introns.insert(a);
        }
    }
    
    return includes(hit_introns.begin(),hit_introns.end(), scaffold_introns.begin(), scaffold_introns.end());
}

bool scaff_lt(const Scaffold& lhs, const Scaffold& rhs)
{
	return lhs.left() < rhs.left();
}

bool scaff_lt_rt(const Scaffold& lhs, const Scaffold& rhs)
{
    if (lhs.left() != rhs.left())
        return lhs.left() < rhs.left();
    return lhs.right() < rhs.right();
}

bool scaff_lt_rt_oplt(const Scaffold& lhs, const Scaffold& rhs)
{
    if (scaff_lt_rt(lhs, rhs))
        return true;
    if (scaff_lt_rt(rhs, lhs))
        return false;
    
    // Now we need to actually compare the alignment ops
    const vector<AugmentedCuffOp>& lhs_ops = lhs.augmented_ops();
    const vector<AugmentedCuffOp>& rhs_ops = rhs.augmented_ops();
    
    if (lhs_ops.size() != rhs_ops.size())
        return lhs_ops.size() < rhs_ops.size();
    
    for (size_t i = 0; i < lhs_ops.size(); ++i)
    {
        if (lhs_ops[i] != rhs_ops[i])
        {
            return lhs_ops[i] < rhs_ops[i];
        }
    }
    
    return false;
}

bool scaff_lt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs)
{
	return scaff_lt(*lhs,*rhs);
}

bool scaff_lt_rt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs)
{
	return scaff_lt_rt(*lhs,*rhs);
}

bool scaff_lt_rt_oplt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs)
{
	return scaff_lt_rt_oplt(*lhs,*rhs);
}
