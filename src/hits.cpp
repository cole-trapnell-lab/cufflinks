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

#if ENABLE_THREADS
boost::mutex RefSequenceTable::table_lock;
#endif

int num_deleted = 0;

void ReadHit::trim(int trimmed_length)
{
    bool antisense_aln = _sam_flag & 0x10;
    
    vector<CigarOp> new_cigar;
    int new_left = 0;
    
    if (!antisense_aln)
    {
        int pos = _left;
        new_left = _left;
        int length = 0;
		for (vector<CigarOp>::iterator i = _cigar.begin(); i < _cigar.end(); ++i)
		{
			const CigarOp& op = *i;
            
            if (length < trimmed_length)
            {
                switch(op.opcode)
                {
                    case REF_SKIP:
                        //gaps_out.push_back(make_pair(pos, pos + op.length - 1));
                        pos += op.length;
                        new_cigar.push_back(op);
                        break;
                    case SOFT_CLIP:
                        assert(false); // not sure if this case is right
                        pos += op.length;
                        length += op.length;
                        new_cigar.push_back(op);
                        break;
                    case HARD_CLIP:
                        new_cigar.push_back(op);
                        break;
                    case MATCH:
                        if (length + op.length < trimmed_length)
                        {
                            pos += op.length;
                            length += op.length;
                            new_cigar.push_back(op);
                        }
                        else
                        {
                            new_cigar.push_back(CigarOp(MATCH, trimmed_length - length));
                            pos += trimmed_length - length;
                            length += trimmed_length - length;
                        }
                        break;
                    case INS:
                        assert(false); // not sure if this case is right
                        pos -= op.length;
                        length -= op.length;
                        new_cigar.push_back(op);
                        break;
                    case DEL:
                        assert(false); // not sure if this case is right
                        pos += op.length;
                        length += op.length;
                        new_cigar.push_back(op);
                        break;
                    default:
                        break;
                }
            }
		}
    }
    else
    {
        int pos = _right;
        int length = 0;
		for (vector<CigarOp>::reverse_iterator i = _cigar.rbegin(); i < _cigar.rend(); ++i)
		{
			const CigarOp& op = *i;
            
            if (length < trimmed_length)
            {
                switch(op.opcode)
                {
                    case REF_SKIP:
                        //gaps_out.push_back(make_pair(pos, pos + op.length - 1));
                        pos -= op.length;
                        new_cigar.push_back(op);
                        break;
                    case SOFT_CLIP:
                        assert(false); // not sure if this case is right
                        pos -= op.length;
                        length += op.length;
                        new_cigar.push_back(op);
                        break;
                    case HARD_CLIP:
                        new_cigar.push_back(op);
                        break;
                    case MATCH:
                        if (length + op.length < trimmed_length)
                        {
                            pos -= op.length;
                            length += op.length;
                            new_cigar.push_back(op);
                        }
                        else
                        {
                            new_cigar.push_back(CigarOp(MATCH, trimmed_length - length));
                            pos -= trimmed_length - length;
                            length += trimmed_length - length;
                        }
                        break;
                    case INS:
                        assert(false); // not sure if this case is right
                        pos += op.length;
                        length -= op.length;
                        new_cigar.push_back(op);
                        break;
                    case DEL:
                        assert(false); // not sure if this case is right
                        pos -= op.length;
                        length += op.length;
                        new_cigar.push_back(op);
                        break;
                    default:
                        break;
                }
            }
		}
        _left = pos;
    }
    _cigar = new_cigar;
    _right = get_right();
    assert (trimmed_length == read_len());
}

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

// Compares for structural equality, but won't declare multihits equal to one another
bool hits_eq_non_multi(const MateHit& lhs, const MateHit& rhs)
{
	if ((lhs.is_multi() || rhs.is_multi() ) && lhs.insert_id() != rhs.insert_id())
		return false;
	return hits_equals(lhs, rhs);
}

// Compares for structural equality, but won't declare multihits equal to one another
// and won't return true for hits from different read groups (e.g. replicate samples)
bool hits_eq_non_multi_non_replicate(const MateHit& lhs, const MateHit& rhs)
{
	if ((lhs.is_multi() || rhs.is_multi() || lhs.read_group_props() != rhs.read_group_props()) && lhs.insert_id() != rhs.insert_id())
		return false;
	return hits_equals(lhs, rhs);
}
    
// Does NOT care about the read group this hit came from.
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

bool has_no_collapse_mass(const MateHit& hit)
{
	return hit.collapse_mass() == 0;
}

// Assumes hits are sorted by mate_hit_lt
// Does not collapse hits that are multi-reads
void collapse_hits(const vector<MateHit>& hits,
				   vector<MateHit>& non_redundant)
{
	copy(hits.begin(), hits.end(), back_inserter(non_redundant));
	vector<MateHit>::iterator new_end = unique(non_redundant.begin(), 
											   non_redundant.end(), 
											   hits_eq_non_multi_non_replicate);
	non_redundant.erase(new_end, non_redundant.end());
    non_redundant.resize(non_redundant.size());
	
	foreach(MateHit& hit, non_redundant)
		hit.collapse_mass(0);
	
	size_t curr_aln = 0;
	size_t curr_unique_aln = 0;
	while (curr_aln < hits.size())
	{
		if (hits_eq_non_multi_non_replicate(non_redundant[curr_unique_aln], hits[curr_aln]) || hits_eq_non_multi_non_replicate(non_redundant[++curr_unique_aln], hits[curr_aln]))
		{
            double more_mass = hits[curr_aln].common_scale_mass();
			//assert(non_redundant[curr_unique_aln].collapse_mass() == 0 || !non_redundant[curr_unique_aln].is_multi());
			non_redundant[curr_unique_aln].incr_collapse_mass(more_mass);
		}
		else
			assert(false);
		
		++curr_aln;
	}
	
	//foreach(MateHit& hit, non_redundant)
		//assert(hit.collapse_mass() <= 1 || !hit.is_multi());
	
	//non_redundant.erase(remove_if(non_redundant.begin(),non_redundant.end(),has_no_collapse_mass), non_redundant.end()); 
	
}

// Places multi-reads to the right of reads they match
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
	
	const ReadHit* lhs_l = lhs.left_alignment();
	const ReadHit* lhs_r = lhs.right_alignment();
	
	const ReadHit* rhs_l = rhs.left_alignment();
	const ReadHit* rhs_r = rhs.right_alignment();
	
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
	
	if (lhs.is_multi() != rhs.is_multi())
	{
		return rhs.is_multi();
	}
	
	return false;
}


ReadHit HitFactory::create_hit(const string& insert_name, 
							   const string& ref_name,
							   int left,
							   const vector<CigarOp>& cigar,
							   CuffStrand source_strand,
							   const string& partner_ref,
							   int partner_pos, 
							   unsigned int edit_dist,
							   int num_hits,
                               float base_mass,
                               uint32_t sam_flag)
{
	InsertID insert_id = _insert_table.get_id(insert_name);
	RefID reference_id = _ref_table.get_id(ref_name, NULL);
	RefID partner_ref_id = _ref_table.get_id(partner_ref, NULL);
	
	return ReadHit(reference_id,
				   insert_id,
				   left, 
				   cigar, 
				   source_strand,
				   partner_ref_id,
				   partner_pos,
				   edit_dist,
				   num_hits,
                   base_mass,
                   sam_flag);	
}

ReadHit HitFactory::create_hit(const string& insert_name, 
							   const string& ref_name,
							   uint32_t left,
							   uint32_t read_len,
							   CuffStrand source_strand,
							   const string& partner_ref,
							   int partner_pos,
							   unsigned int edit_dist,
							   int num_hits,
                               float base_mass,
                               uint32_t sam_flag)
{
	InsertID insert_id = _insert_table.get_id(insert_name);
	RefID reference_id = _ref_table.get_id(ref_name, NULL);
	RefID partner_ref_id = _ref_table.get_id(partner_ref, NULL);
	
	return ReadHit(reference_id,
				   insert_id, 
				   left,
				   read_len,
				   source_strand,
				   partner_ref_id,
				   partner_pos,
				   edit_dist,
				   num_hits,
                   base_mass,
                   sam_flag);	
}

// populate a bam_t This will 
bool BAMHitFactory::next_record(const char*& buf, size_t& buf_size)
{
    if (_next_hit.data)
    {
        free(_next_hit.data);
        _next_hit.data = NULL;
    }
    
    if (records_remain() == false)
        return false;
    
	mark_curr_pos();
    
    memset(&_next_hit, 0, sizeof(_next_hit));
    
	int bytes_read = samread(_hit_file, &_next_hit);
	if (bytes_read < 0)
    {
        _eof_encountered = true;
		return false;
    }
	buf = (const char*)&_next_hit;
	buf_size = bytes_read;
	
	return true;
}

CuffStrand use_stranded_protocol(uint32_t sam_flag,  MateStrandMapping msm)
{
    bool antisense_aln = sam_flag & 0x10;
	if (((sam_flag & BAM_FPAIRED) && (sam_flag & BAM_FREAD1)) || !(sam_flag & BAM_FPAIRED)) // first-in-pair or single-end
	{
		switch(msm)
		{
			case FF:
			case FR:
				return (antisense_aln) ? CUFF_REV : CUFF_FWD;
				break;
			case RF:
			case RR:
				return (antisense_aln) ? CUFF_FWD : CUFF_REV;
				break;
		}
	}
	else // second-in-pair read
	{
		switch (msm)
		{
			case FF:
			case RF:
				return (antisense_aln) ? CUFF_REV : CUFF_FWD;
				break;
			case FR:
			case RR:
				return (antisense_aln) ? CUFF_FWD : CUFF_REV;
				break;
		}
	}
	assert(false);
    return CUFF_STRAND_UNKNOWN;
}


bool BAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
									 ReadHit& bh,
									 bool strip_slash,
									 char* name_out,
									 char* name_tags)
{
	const bam1_t* hit_buf = (const bam1_t*)orig_bwt_buf;
	
	uint32_t sam_flag = hit_buf->core.flag;
	
	int text_offset = hit_buf->core.pos;
	int text_mate_pos = hit_buf->core.mpos;
	int target_id = hit_buf->core.tid;
	int mate_target_id = hit_buf->core.mtid;
	
	vector<CigarOp> cigar;
	bool spliced_alignment = false;
	int num_hits = 1;
	
	//header->target_name[c->tid]
	
	if (sam_flag & 0x4 || target_id < 0)
	{
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(bam1_qname(hit_buf),
						"*",
						0, // SAM files are 1-indexed
						0,
						CUFF_STRAND_UNKNOWN,
						"*",
						0,
						0,
						1,
                        1.0,
                        sam_flag);
		return true;
	}
	if (target_id >= _hit_file->header->n_targets)
    {
        fprintf (stderr, "BAM error: file contains hits to sequences not in header SQ records (%s)\n", bam1_qname(hit_buf));
        return false;
    }
    
	string text_name = _hit_file->header->target_name[target_id];
	
	for (int i = 0; i < hit_buf->core.n_cigar; ++i) 
	{
		//char* t;

		int length = bam1_cigar(hit_buf)[i] >> BAM_CIGAR_SHIFT;
		if (length <= 0)
		{
		  fprintf (stderr, "BAM error: CIGAR op has zero length (%s)\n", bam1_qname(hit_buf));
			return false;
		}
		
		CigarOpCode opcode;
		switch(bam1_cigar(hit_buf)[i] & BAM_CIGAR_MASK)
		{
			case BAM_CMATCH: opcode  = MATCH; break; 
			case BAM_CINS: opcode  = INS; break;
			case BAM_CDEL: opcode  = DEL; break; 
			case BAM_CSOFT_CLIP: opcode  = SOFT_CLIP; break;
			case BAM_CHARD_CLIP: opcode  = HARD_CLIP; break;
			case BAM_CPAD: opcode  = PAD; break; 
			case BAM_CREF_SKIP:
                opcode = REF_SKIP;
				spliced_alignment = true;
				if (length > (int)max_intron_length)
				{
					//fprintf(stderr, "Encounter REF_SKIP > max_gene_length, skipping\n");
					return false;
				}
				break;
			default:
				//fprintf (stderr, "SAM error on line %d: invalid CIGAR operation\n", _line_num);
				return false;
		}
		if (opcode != HARD_CLIP)
			cigar.push_back(CigarOp(opcode, length));
	}
	
	string mrnm;
	if (mate_target_id >= 0)
	{
		if (mate_target_id == target_id)
		{
			mrnm = _hit_file->header->target_name[mate_target_id];
//			if (abs((int)text_mate_pos - (int)text_offset) > (int)max_intron_length)
//			{
//				//fprintf (stderr, "Mates are too distant, skipping\n");
//				return false;
//			}
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

	uint8_t* ptr = bam_aux_get(hit_buf, "XS");
	if (ptr)
	{
		char src_strand_char = bam_aux2A(ptr);
		if (src_strand_char == '-')
			source_strand = CUFF_REV;
		else if (src_strand_char == '+')
			source_strand = CUFF_FWD;
	}
	
	ptr = bam_aux_get(hit_buf, "NM");
	if (ptr)
	{
		num_mismatches = bam_aux2i(ptr);
	}

	ptr = bam_aux_get(hit_buf, "NH");
	if (ptr)
	{
		num_hits = bam_aux2i(ptr);
	}
    
    double mass = 1.0;
    ptr = bam_aux_get(hit_buf, "ZF");
	if (ptr)
	{
		mass = bam_aux2i(ptr);
        if (mass <= 0.0)
            mass = 1.0;
	}
	    
    if (_rg_props.strandedness() == STRANDED_PROTOCOL && source_strand == CUFF_STRAND_UNKNOWN)
		source_strand = use_stranded_protocol(sam_flag, _rg_props.mate_strand_mapping());
    
	if (!spliced_alignment)
	{
		//assert(_rg_props.strandedness() == STRANDED_PROTOCOL || source_strand == CUFF_STRAND_UNKNOWN);

		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(bam1_qname(hit_buf),
						text_name,
						text_offset,  // BAM files are 0-indexed
						cigar,
						source_strand,
						mrnm,
						text_mate_pos,
						num_mismatches,
						num_hits,
                        mass,
                        sam_flag);
		return true;
		
	}
	else
	{	
		if (source_strand == CUFF_STRAND_UNKNOWN)
		{
			fprintf(stderr, "BAM record error: found spliced alignment without XS attribute\n");
		}
		
		bh = create_hit(bam1_qname(hit_buf),
						text_name,
						text_offset,  // BAM files are 0-indexed
						cigar,
						source_strand,
						mrnm,
						text_mate_pos,
						num_mismatches,
						num_hits,
                        mass,
                        sam_flag);
		return true;
	}
	
	
	return true;
}



Platform str_to_platform(const string pl_str)
{
    if (pl_str == "SOLiD")
    {
        return SOLID;
    }
    else if (pl_str == "Illumina") 
    {
        return ILLUMINA;
    }
    else 
    {
        return UNKNOWN_PLATFORM;
    }
}

// Parses the header to determine platform and other properties
bool HitFactory::parse_header_string(const string& header_rec,
                                     ReadGroupProperties& rg_props)
{
    vector<string> columns;
    tokenize(header_rec, "\t", columns); 
    
    if (columns[0] == "@RG")
    {
        for (size_t i = 1; i < columns.size(); ++i)
        {
            vector<string> fields;
            tokenize(columns[i], ":", fields);
            if (fields[0] == "PL")
            {
                if (rg_props.platform() == UNKNOWN_PLATFORM)
                {
                    Platform p = str_to_platform(fields[1]);
                    rg_props.platform(p);
                }
                else 
                {
                    Platform p = str_to_platform(fields[1]);
                    if (p != rg_props.platform())
                    {
                        fprintf(stderr, "Error: Processing reads from different platforms is not currently supported\n");
                        return false;
                    }
                }

            }
        }
    }
    else if (columns[0] == "@SQ")
    {
        _num_seq_header_recs++;
        for (size_t i = 1; i < columns.size(); ++i)
        {
            vector<string> fields;
            tokenize(columns[i], ":", fields);
            if (fields[0] == "SN")
            {
                // Populate the RefSequenceTable with the sequence dictionary, 
                // to ensure that (for example) downstream GTF files are sorted
                // in an order consistent with the header, and to enforce that
                // BAM records appear in the order implied by the header
                RefID _id = _ref_table.get_id(fields[1], NULL);

                const RefSequenceTable::SequenceInfo* info = _ref_table.get_info(_id);
				
                if (info->observation_order != _num_seq_header_recs)
                {
                    if (info->name != fields[1])
                    {
                        fprintf(stderr, "Error: Hash collision between references '%s' and '%s'.\n", info->name, fields[1].c_str());
                    }
                    else
                    {
                        fprintf(stderr, "Error: sort order of reads in BAMs must be the same\n");
					}
                    exit(1);
                }
            }
        }
    }

    return true;
}

void HitFactory::finalize_rg_props()
{
    if (_rg_props.platform() == SOLID)
    {
        _rg_props.strandedness(STRANDED_PROTOCOL);
        _rg_props.std_mate_orientation(MATES_POINT_TOWARD);
    }
    else
    {
        // Default to Illumina's unstranded protocol params for strandedness and
        // mate orientation
        _rg_props.strandedness(UNSTRANDED_PROTOCOL);
        _rg_props.std_mate_orientation(MATES_POINT_TOWARD);
    }
}

static const unsigned MAX_HEADER_LEN = 4 * 1024 * 1024; // 4 MB

bool BAMHitFactory::inspect_header()
{
    bam_header_t* header = _hit_file->header;
    
    if (header == NULL)
    {
        fprintf(stderr, "Warning: No BAM header\n");
        return false;
    }
    
    if (header->l_text >= MAX_HEADER_LEN)
    {
        fprintf(stderr, "Warning: BAM header too large\n");
        return false;
    }

	if (header->l_text == 0)
	{
		fprintf(stderr, "Warning: BAM header has 0 length or is corrupted.  Try using 'samtools reheader'.\n");
        return false;
	}
	
	
    if (header->text != NULL)
    {
        char* h_text = strdup(header->text);
        char* pBuf = h_text;
        while(pBuf - h_text < header->l_text)
        {
            char* nl = strchr(pBuf, '\n');
            if (nl) 
            {
                *nl = 0; 
                parse_header_string(pBuf, _rg_props);
                pBuf = ++nl;
            }
            else 
            {
                pBuf = h_text + header->l_text;
            }
        }
        
        free(h_text);
    }
    
    finalize_rg_props();
    return true;
}


bool SAMHitFactory::next_record(const char*& buf, size_t& buf_size)
{
	mark_curr_pos();
	
	bool new_rec = fgets(_hit_buf,  _hit_buf_max_sz - 1, _hit_file);
	if (!new_rec)
		return false;
	++_line_num;
	char* nl = strrchr(_hit_buf, '\n');
	if (nl) *nl = 0;
	buf = _hit_buf;
	buf_size = _hit_buf_max_sz - 1;
	return true;
}

bool SAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
									 ReadHit& bh,
									 bool strip_slash,
									 char* name_out,
									 char* name_tags)
{	
	char bwt_buf[10*2048];

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
	int num_hits = 1;
	    
    if ((sam_flag & 0x4) ||!strcmp(text_name, "*"))
	{
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(name,
						"*",
						0, // SAM files are 1-indexed
						0,
						CUFF_STRAND_UNKNOWN,
						"*",
						0,
						0,
						1,
                        1.0,
                        sam_flag);
		return true;
	}
	// Mostly pilfered direct from the SAM tools:
	while (*p_cig) 
	{
		char* t;
		int length = (int)strtol(p_cig, &t, 10);
		if (length <= 0)
		{
			fprintf (stderr, "SAM error on line %d: CIGAR op has zero length\n", _line_num);
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
		else if (op_char == 'D') 
        {
            opcode = DEL;
        }
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
			fprintf (stderr, "SAM error on line %d: invalid CIGAR operation\n", _line_num);
			return false;
		}
		p_cig = t + 1;
		//i += length;
		if (opcode != HARD_CLIP)
			cigar.push_back(CigarOp(opcode, length));
	}
	if (*p_cig)
	{
		fprintf (stderr, "SAM error on line %d: unmatched CIGAR operation\n", _line_num);
		return false;
	}
    
	string mrnm;
	if (strcmp(mate_ref_name, "*"))
	{
		if (!strcmp(mate_ref_name, "=") || !strcmp(mate_ref_name, text_name))
		{
			mrnm = text_name;
//			if (abs((int)text_mate_pos - (int)text_offset) > (int)max_intron_length)
//			{
//				//fprintf (stderr, "Mates are too distant, skipping\n");
//				return false;
//			}
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
	
    double mass = 1.0;
    
	while((tag_buf = strsep((char**)&buf,"\t")))
	{
		
		char* first_colon = (char*)strchr(tag_buf, ':');
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
				else if (!strcmp(first_token, "NH"))
				{
                    num_hits = atoi(third_token);
				}
                else if (!strcmp(first_token, "ZF"))
				{
					mass = atof(third_token);
                    if (mass <= 0.0)
                        mass = 1.0;
				}
				else 
				{
					
				}
			}
		}
	}
        
    // Don't let the protocol setting override explicit XS tags
	if (_rg_props.strandedness() == STRANDED_PROTOCOL && source_strand == CUFF_STRAND_UNKNOWN)
		source_strand = use_stranded_protocol(sam_flag, _rg_props.mate_strand_mapping());
	
	if (!spliced_alignment)
	{		
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(name,
						text_name,
						text_offset - 1,
						cigar,
						source_strand,
						mrnm,
						text_mate_pos - 1,
						num_mismatches,
						num_hits,
                        mass,
                        sam_flag);
		return true;
		
	}
	else
	{	
		if (source_strand == CUFF_STRAND_UNKNOWN)
		{
			fprintf(stderr, "SAM error on line %d: found spliced alignment without XS attribute\n", _line_num);
		}
		
		bh = create_hit(name,
						text_name,
						text_offset - 1,
						cigar,
						source_strand,
						mrnm,
						text_mate_pos - 1,
						num_mismatches,
						num_hits,
                        mass,
                        sam_flag);
		return true;
	}
	return false;
}

bool SAMHitFactory::inspect_header()
{
    char pBuf[10 * 1024];
    
    off_t curr_pos = ftello(_hit_file);
    rewind(_hit_file);
    
    while (fgets(pBuf, 10*1024, _hit_file))
    {
        if (pBuf[0] != '@')
        {
            break; // done with the header.
        }
        char* nl = strchr(pBuf, '\n');
        if (nl) 
        {
            *nl = 0; 
            parse_header_string(pBuf, _rg_props);
        }
    }
    
    fseek(_hit_file, curr_pos, SEEK_SET);
    
    finalize_rg_props();
    return true;
}
