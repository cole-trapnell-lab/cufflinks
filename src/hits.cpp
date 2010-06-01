/*
 *  hits.cpp
 *  Cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <vector>

#include "common.h"
#include "hits.h"
#include "tokenize.h"

using namespace std;

//static const int max_read_length = 1024;

bool hit_insert_id_lt(const ReadHit& h1, const ReadHit& h2)
{
	return h1.insert_id() < h2.insert_id();
}

bool hits_eq_mod_id(const ReadHit& lhs, const ReadHit& rhs)
{
	return (lhs.ref_id() == rhs.ref_id() &&
			lhs.antisense_align() == rhs.antisense_align() &&
			lhs.left() == rhs.left() && 
			lhs.source_strand() == rhs.source_strand() &&
			lhs.cigar() == rhs.cigar());
}

bool hits_equals(const MateHit& lhs, const MateHit& rhs) 
{
	if (lhs.ref_id() != rhs.ref_id())
		return false;
	
	if ((lhs.left_alignment() == NULL) != (rhs.left_alignment() == NULL))
		return false;
	if ((lhs.right_alignment() == NULL) != (rhs.right_alignment() == NULL))
		return false;
	if (lhs.left_alignment())
	{
		if (!(hits_eq_mod_id(*lhs.left_alignment(),*(rhs.left_alignment()))))
			return false;
	}
	if (lhs.right_alignment())
	{
		if (!(hits_eq_mod_id(*lhs.right_alignment(),*(rhs.right_alignment()))))
			return false;
	}
	return true;
}

// Assumes hits are sorted by mate_hit_lt
void collapse_hits(const vector<MateHit>& hits,
				   vector<MateHit>& non_redundant,
				   vector<int>& collapse_counts)
{
	copy(hits.begin(), hits.end(), back_inserter(non_redundant));
	vector<MateHit>::iterator new_end = unique(non_redundant.begin(), 
											   non_redundant.end(), 
											   hits_equals);
	non_redundant.erase(new_end, non_redundant.end());
	
	collapse_counts = vector<int>(non_redundant.size(), 1);
	size_t curr_aln = 0;
	size_t curr_unique_aln = 0;
	while (curr_aln < collapse_counts.size())
	{
		if (hits_equals(non_redundant[curr_unique_aln], hits[curr_aln]))
			collapse_counts[curr_unique_aln]++;
		else
		{
			curr_unique_aln++;
		}
		
		++curr_aln;
	}
	
}

bool mate_hit_lt(const MateHit& lhs, const MateHit& rhs)
{
	if (lhs.left() != rhs.left())
		return lhs.left() < rhs.left();
	if (lhs.right() != rhs.right())
		return lhs.right() > rhs.right();
	
	if ((lhs.left_alignment() == NULL) != (rhs.left_alignment() == NULL))
		return (lhs.left_alignment() == NULL) < (rhs.left_alignment() == NULL);
	
	if ((lhs.right_alignment() == NULL) != (rhs.right_alignment() == NULL))
		return (lhs.right_alignment() == NULL) < (rhs.right_alignment() == NULL);
	
	assert ((lhs.right_alignment() == NULL) == (rhs.right_alignment() == NULL));
	assert ((lhs.left_alignment() == NULL) == (rhs.left_alignment() == NULL));
	
	shared_ptr<const ReadHit> lhs_l = lhs.left_alignment();
	shared_ptr<const ReadHit> lhs_r = lhs.right_alignment();
	
	shared_ptr<const ReadHit> rhs_l = rhs.left_alignment();
	shared_ptr<const ReadHit> rhs_r = rhs.right_alignment();
	
	if (lhs_l && rhs_l)
	{
		if (lhs_l->cigar().size() != rhs_l->cigar().size())
			return lhs_l->cigar().size() < rhs_l->cigar().size(); 
		for (size_t i = 0; i < lhs_l->cigar().size(); ++i)
		{
			if (lhs_l->cigar()[i].opcode != rhs_l->cigar()[i].opcode)
				return lhs_l->cigar()[i].opcode < rhs_l->cigar()[i].opcode;
			if (lhs_l->cigar()[i].length != rhs_l->cigar()[i].length)
				return lhs_l->cigar()[i].length < rhs_l->cigar()[i].length;
		}
	}
	
	if (lhs_r && rhs_r)
	{
		if (lhs_r->cigar().size() != rhs_r->cigar().size())
			return lhs_r->cigar().size() < rhs_r->cigar().size(); 
		for (size_t i = 0; i < lhs_r->cigar().size(); ++i)
		{
			if (lhs_r->cigar()[i].opcode != rhs_r->cigar()[i].opcode)
				return lhs_r->cigar()[i].opcode < rhs_r->cigar()[i].opcode;
			if (lhs_r->cigar()[i].length != rhs_r->cigar()[i].length)
				return lhs_r->cigar()[i].length < rhs_r->cigar()[i].length;
		}
	}
	
	return false;
}


ReadHit HitFactory::create_hit(const string& insert_name, 
							   const string& ref_name,
							   int left,
							   const vector<CigarOp>& cigar,
							   bool antisense_aln,
							   CuffStrand source_strand,
							   const string& partner_ref,
							   int partner_pos, 
							   double error_prob,
							   unsigned int edit_dist)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL);
	uint32_t partner_ref_id = _ref_table.get_id(partner_ref, NULL);
	
	return ReadHit(reference_id,
				   insert_id, 
				   left, 
				   cigar, 
				   antisense_aln,
				   source_strand,
				   partner_ref_id,
				   partner_pos,
				   error_prob,
				   edit_dist);	
}

ReadHit HitFactory::create_hit(const string& insert_name, 
							   const string& ref_name,
							   uint32_t left,
							   uint32_t read_len,
							   bool antisense_aln,
							   CuffStrand source_strand,
							   const string& partner_ref,
							   int partner_pos,
							   double error_prob,
							   unsigned int edit_dist)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL);
	uint32_t partner_ref_id = _ref_table.get_id(partner_ref, NULL);
	
	return ReadHit(reference_id,
				   insert_id, 
				   left,
				   read_len,
				   antisense_aln,
				   source_strand,
				   partner_ref_id,
				   partner_pos,
				   error_prob,
				   edit_dist);	
}

int anchor_mismatch = 0;

bool SAMHitFactory::get_hit_from_buf(int line_num, 
									 const char* orig_bwt_buf, 
									 ReadHit& bh,
									 bool strip_slash,
									 char* name_out,
									 char* name_tags)
{
	char bwt_buf[2048];
	strcpy(bwt_buf, orig_bwt_buf);
	// Are we still in the header region?
	if (bwt_buf[0] == '@')
		return false;
	
	const char* buf = bwt_buf;
	const char* _name = strsep((char**)&buf,"\t");
	if (!_name)
		return false;
	char name[2048];
	strncpy(name, _name, 2047); 
	
	const char* sam_flag_str = strsep((char**)&buf,"\t");
	if (!sam_flag_str)
		return false;
	
	const char* text_name = strsep((char**)&buf,"\t");
	if (!text_name)
		return false;
	
	const char* text_offset_str = strsep((char**)&buf,"\t");
	if (!text_offset_str)
		return false;
	
	const char* map_qual_str =  strsep((char**)&buf,"\t");
	if (!map_qual_str)
		return false;
	
	const char* cigar_str = strsep((char**)&buf,"\t");
	if (!cigar_str)
		return false;
	
	const char* mate_ref_name =  strsep((char**)&buf,"\t");
	if (!mate_ref_name)
		return false;
	
	const char* mate_pos_str =  strsep((char**)&buf,"\t");
	if (!mate_pos_str)
		return false;
	
	const char* inferred_insert_sz_str =  strsep((char**)&buf,"\t");
	if (!inferred_insert_sz_str)
		return false;
	
	const char* seq_str =  strsep((char**)&buf,"\t");
	if (!seq_str)
		return false;
	
	const char* qual_str =  strsep((char**)&buf,"\t");
	if (!qual_str)
		return false;
	
	
	int sam_flag = atoi(sam_flag_str);
	int text_offset = atoi(text_offset_str);
	int text_mate_pos = atoi(mate_pos_str);
	
	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
			strcpy(name_tags, pipe);
		*pipe = 0;
	}
	// Stripping the slash and number following it gives the insert name
	char* slash = strrchr(name, '/');
	if (strip_slash && slash)
		*slash = 0;
	
	const char* p_cig = cigar_str;
	//int len = strlen(sequence);
	vector<CigarOp> cigar;
	bool spliced_alignment = false;
	
	double mapQ = atoi(map_qual_str);
	long double error_prob;
	if (mapQ > 0)
	{
		long double p = (-1.0 * mapQ) / 10.0;
		error_prob = pow(10.0L, p);
	}
	else
	{
		error_prob = 1.0;
	}
	
	if (!strcmp(text_name, "*"))
	{
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(name,
						text_name,
						0, // SAM files are 1-indexed
						0,
						false,
						CUFF_STRAND_UNKNOWN,
						"*",
						0,
						1.0,
						0);
		return true;
	}
	// Mostly pilfered direct from the SAM tools:
	while (*p_cig) 
	{
		char* t;
		int length = (int)strtol(p_cig, &t, 10);
		if (length <= 0)
		{
			fprintf (stderr, "SAM error on line %d: CIGAR op has zero length\n", line_num);
			return false;
		}
		char op_char = toupper(*t);
		CigarOpCode opcode;
		if (op_char == 'M') 
		{
			/*if (length > max_read_length)
			{
				fprintf(stderr, "SAM error on line %d:  %s: MATCH op has length > %d\n", line_num, name, max_read_length);
				return false;
			}*/
			opcode = MATCH;
		}
		else if (op_char == 'I') opcode = INS;
		else if (op_char == 'D') opcode = DEL;
		else if (op_char == 'N')
		{
			opcode = REF_SKIP;
			spliced_alignment = true;
			if (length > (int)max_intron_length)
			{
				//fprintf(stderr, "Encounter REF_SKIP > max_gene_length, skipping\n");
				return false;
			}
		}
		else if (op_char == 'S') opcode = SOFT_CLIP;
		else if (op_char == 'H') opcode = HARD_CLIP;
		else if (op_char == 'P') opcode = PAD;
		else
		{
			fprintf (stderr, "SAM error on line %d: invalid CIGAR operation\n", line_num);
			return false;
		}
		p_cig = t + 1;
		//i += length;
		cigar.push_back(CigarOp(opcode, length));
	}
	if (*p_cig)
	{
		fprintf (stderr, "SAM error on line %d: unmatched CIGAR operation\n", line_num);
		return false;
	}
    
	string mrnm;
	if (strcmp(mate_ref_name, "*"))
	{
		if (!strcmp(mate_ref_name, "=") || !strcmp(mate_ref_name, text_name))
		{
			mrnm = text_name;
			if (abs((int)text_mate_pos - (int)text_offset) > max_intron_length)
			{
				//fprintf (stderr, "Mates are too distant, skipping\n");
				return false;
			}
		}
		else
		{
			//fprintf(stderr, "Trans-spliced mates are not currently supported, skipping\n");
			return false;
		}
	}
	else
	{
		text_mate_pos = 0;
	}
	
	CuffStrand source_strand = CUFF_STRAND_UNKNOWN;
	unsigned char num_mismatches = 0;
	
	const char* tag_buf = buf;
	
	while((tag_buf = strsep((char**)&buf,"\t")))
	{
		
		char* first_colon = strchr(tag_buf, ':');
		if (first_colon)
		{
			*first_colon = 0;
			++first_colon;
			char* second_colon = strchr(first_colon, ':');
			if (second_colon)
			{
				*second_colon = 0;
				++second_colon;
				const char* first_token = tag_buf;
				//const char* second_token = first_colon;
				const char* third_token = second_colon;
				if (!strcmp(first_token, "XS"))
				{				
					if (*third_token == '-')
						source_strand = CUFF_REV;
					else if (*third_token == '+')
						source_strand = CUFF_FWD;
				}
				else if (!strcmp(first_token, "NM"))
				{
					num_mismatches = atoi(third_token);
				}
				else 
				{
					
				}
			}
		}
	}
	
	if (!spliced_alignment)
	{
		
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(name,
						text_name,
						text_offset - 1, // SAM files are 1-indexed
						cigar[0].length,
						sam_flag & 0x0010,
						source_strand,
						mrnm,
						text_mate_pos - 1,
						error_prob,
						num_mismatches);
		return true;
		
	}
	else
	{	
		if (source_strand == CUFF_STRAND_UNKNOWN)
		{
			fprintf(stderr, "SAM error on line %d: found spliced alignment without XS attribute\n", line_num);
		}
		
		bh = create_hit(name,
						text_name,
						text_offset - 1,
						cigar,
						sam_flag & 0x0010,
						source_strand,
						mrnm,
						text_mate_pos - 1,
						error_prob,
						num_mismatches);
		return true;
	}
	return false;
}


//void add_hits_to_coverage(const HitList& hits, vector<short>& DoC)
//{
//	int max_hit_pos = -1;
//	for (size_t i = 0; i < hits.size(); ++i)
//	{
//		max_hit_pos = max((int)hits[i].right(),max_hit_pos);
//	}
//	
//	if ((int)DoC.size() < max_hit_pos)
//		DoC.resize(max_hit_pos);
//	
//	for (size_t i = 0; i < hits.size(); ++i)
//	{
//		const ReadHit& bh = hits[i];
//		
//		if (!bh.accepted())
//			continue;
//		// split up the coverage contibution for this reads
//		size_t j = bh.left();
//		const vector<CigarOp>& cigar = bh.cigar();
//
//		for (size_t c = 0 ; c < cigar.size(); ++c)
//		{
//			switch(cigar[c].opcode)
//			{
//				case MATCH:
//					for (size_t m = 0; m < cigar[c].length; ++m)
//						DoC[j + m]++;
//				//fall through this case to REF_SKIP is intentional
//				case REF_SKIP:
//					j += cigar[c].length;
//					break;
//				default:
//					break;
//			}
//			 
//		}
//	}
//}
