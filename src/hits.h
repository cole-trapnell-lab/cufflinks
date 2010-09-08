#ifndef BWT_MAP_H
#define BWT_MAP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <cassert>

#include <boost/shared_ptr.hpp>

#include <bam/sam.h>

#include "common.h"

using namespace std;
using boost::shared_ptr;

/*
 *  hits.h
 *  Cufflinks
 *
 *  Created by Cole Trapnell on 3/23/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

enum CuffStrand { CUFF_STRAND_UNKNOWN = 0, CUFF_FWD = 1, CUFF_REV = 2, CUFF_BOTH = 3 };


enum CigarOpCode 
{ 
	MATCH = BAM_CMATCH, 
	INS = BAM_CINS, 
	DEL = BAM_CDEL, 
	REF_SKIP = BAM_CREF_SKIP,
	SOFT_CLIP = BAM_CSOFT_CLIP, 
	HARD_CLIP = BAM_CHARD_CLIP, 
	PAD = BAM_CPAD
};

struct CigarOp
{
	CigarOp(CigarOpCode o, uint32_t l) : opcode(o), length(l) {}
	CigarOpCode opcode : 3;
	uint32_t length : 29;
	
	bool operator==(const CigarOp& rhs) const { return opcode == rhs.opcode && length == rhs.length; }
	
};

typedef uint64_t InsertID;
typedef uint32_t RefID;

extern int num_deleted;

/*  Stores the information from a single record of the bowtie map. A given read
    may have many of these.  Reads up to 255bp are supported. 
*/
struct ReadHit
{
	ReadHit() : 
		_ref_id(0),
		_insert_id(0),
		_error_prob(1.0),
		_edit_dist(0xFFFFFFFF)
    {
        //num_deleted++;
    }
	
	ReadHit(RefID ref_id,
			InsertID insert_id, 
			int left, 
			int read_len, 
			bool antisense,
			CuffStrand source_strand,
			RefID partner_ref,
			int partner_pos,
			double error_prob,
			unsigned int edit_dist) :
		_ref_id(ref_id),
		_insert_id(insert_id), 
		_left(left), 
		_partner_ref_id(partner_ref),
		_partner_pos(partner_pos),
		_cigar(vector<CigarOp>(1,CigarOp(MATCH,read_len))),
		_source_strand(source_strand),
		_antisense_aln(antisense),
		_error_prob(error_prob),
		_edit_dist(edit_dist)
	{
		assert(_cigar.capacity() == _cigar.size());
		_right = get_right();
        //num_deleted++;
	}
    
	ReadHit(RefID ref_id,
			InsertID insert_id, 
			int left,  
			const vector<CigarOp>& cigar,
			bool antisense_aln,
			CuffStrand source_strand,
			RefID partner_ref,
			int partner_pos, 
			double error_prob,
			unsigned int  edit_dist) : 
		_ref_id(ref_id),
		_insert_id(insert_id), 	
		_left(left),
		_partner_ref_id(partner_ref),
		_partner_pos(partner_pos),
		_cigar(cigar),
		_source_strand(source_strand),
		_antisense_aln(antisense_aln),
		_error_prob(error_prob),
		_edit_dist(edit_dist)
	{
		assert(_cigar.capacity() == _cigar.size());
		_right = get_right();
        //num_deleted++;
	}
	~ReadHit() 
    {
        //num_deleted--;
    }
	int read_len() const
	{
		int len = 0;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			switch(op.opcode)
			{
				case MATCH:
				case INS:
				case SOFT_CLIP:
					len += op.length;
					break;
				default:
					break;
			}
		}
		
		return len;
	}
	
	bool contains_splice() const
	{
                // FIXME: This won't be true once we start supporting INDEL operators in CIGAR strings. See Ticket #183
		return (_cigar.size() > 1);
	}
	
	bool operator==(const ReadHit& rhs) const
	{
	    return (_insert_id == rhs._insert_id &&
	            _ref_id == rhs._ref_id &&
	            _antisense_aln == rhs._antisense_aln &&
	            _left == rhs._left && 
	            _source_strand == rhs._source_strand &&
	            /* DO NOT USE ACCEPTED IN COMPARISON */
	            _cigar == rhs._cigar);
    }
	
	RefID ref_id() const				{ return _ref_id;			}
	InsertID insert_id() const			{ return _insert_id;		}
	
	RefID partner_ref_id() const		{ return _partner_ref_id;	}
	int	  partner_pos()	const			{ return _partner_pos;		}
	
	int left() const					{ return _left;				}
	int right() const					{ return _right;			}
	CuffStrand source_strand()	const	{ return _source_strand; }
	bool antisense_align() const		{ return _antisense_aln;	}
		
	double error_prob() const			{ return _error_prob;		}
	
	// For convenience, if you just want a copy of the gap intervals
	// for this hit.
	void gaps(vector<pair<int,int> >& gaps_out) const
	{
		gaps_out.clear();
		int pos = _left;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			
			switch(op.opcode)
			{
				case REF_SKIP:
					gaps_out.push_back(make_pair(pos, pos + op.length - 1));
					pos += op.length;
					break;
				case MATCH:
					pos += op.length;
					break;
				default:
					break;
			}
		}
	}
	
	const vector<CigarOp>& cigar() const { return _cigar; }
	
	bool contiguous() const 
	{ 
		return _cigar.size() == 1 && _cigar[0].opcode == MATCH;
	}
	
	unsigned int  edit_dist() const { return _edit_dist; }
	
	const string& hitfile_rec() const { return _hitfile_rec; }
	void hitfile_rec(const string& rec) { _hitfile_rec = rec; }
	
private:
	
	int get_right() const	
	{
		int r = _left;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			
			switch(op.opcode)
			{
				case MATCH:
				case REF_SKIP:
				case DEL:
					r += op.length;
					break;
				default:
					break;
			}
		}
		return r;			
	}
	
	RefID _ref_id;
	InsertID _insert_id;   // Id of the sequencing insert
	int _left;        // Position in the reference of the left side of the alignment
	int _right;
	
	RefID _partner_ref_id;  // Reference contig on which we expect the mate 
	int _partner_pos;     // Position at which we expect the mate of this hit
	
	
	vector<CigarOp> _cigar;
	
	CuffStrand _source_strand;    // Which strand the read really came from, if known
	bool _antisense_aln;       // Whether the alignment is to the reverse strand
	double _error_prob;		   // Probability that this alignment is incorrect
	unsigned int  _edit_dist;            // Number of mismatches
	string _hitfile_rec; // Points to the buffer for the record from which this hit came
};

class ReadTable
{
public:
	
	ReadTable() {}
	
	// This function should NEVER return zero
	InsertID get_id(const string& name)
	{
		uint64_t _id = hash_string(name.c_str());
		assert(_id);
		return _id;
	}
	
private:
	
	// This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
	inline uint64_t hash_string(const char* __s)
	{
		uint64_t hash = 0xcbf29ce484222325ull;
		for ( ; *__s; ++__s)
		{
			hash *= 1099511628211ull;
			hash ^= *__s;
		}
		return hash;
	}
};

class RefSequenceTable
{
public:
	
	typedef std::string Sequence;
	
	struct SequenceInfo
	{
		SequenceInfo(uint32_t _order, 
					 char* _name, 
					 Sequence* _seq) :
		observation_order(_order),
		name(_name),
		seq(_seq) {}
		uint32_t observation_order;
		char* name;
		Sequence* seq;
	};
	
	typedef map<string, RefID> IDTable;
	typedef map<RefID, SequenceInfo> InvertedIDTable;
	typedef InvertedIDTable::iterator iterator;
	typedef InvertedIDTable::const_iterator const_iterator;
	
	RefSequenceTable(bool keep_names, bool keep_seqs = false) : 
	_next_id(1), 
	_keep_names(keep_names) {}
	
	~RefSequenceTable()
	{
		for (InvertedIDTable::iterator itr = _by_id.begin();
			 itr != _by_id.end();
			 ++itr)
		{
			free(itr->second.name);
		}
	}
	
	RefID get_id(const string& name,
				 Sequence* seq)
	{
		if (name.empty())
			return 0;
#if ENABLE_THREADS
        table_lock.lock();
#endif
		uint32_t _id = hash_string(name.c_str());
		pair<InvertedIDTable::iterator, bool> ret = 
		_by_id.insert(make_pair(_id, SequenceInfo(_next_id, NULL, NULL)));
		if (ret.second == true)
		{			
			char* _name = NULL;
			if (_keep_names)
				_name = strdup(name.c_str());
			ret.first->second.name = _name;
			ret.first->second.seq	= seq;
			++_next_id;
		}
		assert (_id);
#if ENABLE_THREADS
        table_lock.unlock();
#endif
		return _id;
	}
	
	// You must call invert() before using this function
	const char* get_name(RefID ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.name;
		else
			return NULL;
	}
	
	Sequence* get_seq(RefID ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.seq;
		else
			return NULL;
	}
	
	const SequenceInfo* get_info(RefID ID) const
	{
		
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
		{
			return &(itr->second);
		}
		else
			return NULL;
	}
	
	int observation_order(RefID ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
		{
			return itr->second.observation_order;
		}
		else
			return -1;
	}
	
	iterator begin() { return _by_id.begin(); }
	iterator end() { return _by_id.end(); }
	
	const_iterator begin() const { return _by_id.begin(); }
	const_iterator end() const { return _by_id.end(); }
	
	size_t size() const { return _by_id.size(); }
	
	void clear()
	{
		//_by_name.clear();
		_by_id.clear();
	}
	
private:
	
	// This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
	inline uint32_t hash_string(const char* __s)
	{
		uint32_t hash = 0x811c9dc5;
		for ( ; *__s; ++__s)
		{
			hash *= 16777619;
			hash ^= *__s;
		}
		return hash;
	}
	
	//IDTable _by_name;
	RefID _next_id;
	bool _keep_names;
	InvertedIDTable _by_id;
#if ENABLE_THREADS
    static boost::mutex table_lock;
#endif
};


bool hit_insert_id_lt(const ReadHit& h1, const ReadHit& h2);

/******************************************************************************
 The HitFactory abstract class is responsible for returning a single ReadHit 
 from an alignment file.  The only class that actually implements this interface
 right now in Cufflinks is SAMHitFactory
*******************************************************************************/
class HitFactory
{
public:
    
	HitFactory(ReadTable& insert_table, 
			   RefSequenceTable& reference_table) : 
		_insert_table(insert_table), 
		_ref_table(reference_table) {}
	
	HitFactory& operator=(const HitFactory& rhs) 
	{
		if (this != &rhs)
		{
			//_hit_file = rhs._hit_file;
			_insert_table = rhs._insert_table;
			_ref_table = rhs._ref_table;
		}
		return *this;
	}
	virtual ~HitFactory() {}
	
	ReadHit create_hit(const string& insert_name, 
					   const string& ref_name,
					   int left,
					   const vector<CigarOp>& cigar,
					   bool antisense_aln,
					   CuffStrand source_strand,
					   const string& partner_ref,
					   int partner_pos,
					   double error_prob,
					   unsigned int  edit_dist);
	
	ReadHit create_hit(const string& insert_name, 
					   const string& ref_name,
					   uint32_t left,
					   uint32_t read_len,
					   bool antisense_aln,
					   CuffStrand source_strand,
					   const string& partner_ref,
					   int partner_pos,
					   double error_prob,
					   unsigned int  edit_dist);
	
	virtual void reset() = 0;
	
	virtual void undo_hit() = 0;
	
	// next_record() should always set _curr_pos before reading the
	// next_record so undo_hit() will work properly.
	virtual bool next_record(const char*& buf, size_t& buf_size) = 0;
	
	virtual bool records_remain() const = 0;
	
	virtual bool get_hit_from_buf(const char* bwt_buf, 
								  ReadHit& bh,
								  bool strip_slash,
								  char* name_out = NULL,
								  char* name_tags = NULL) = 0;
	
	RefSequenceTable& ref_table() { return _ref_table; }
	
	//FILE* hit_file() { return _hit_file; }
	
    virtual bool inspect_header() = 0;
    
protected:
    
    bool parse_header_string(const string& header_rec,
                             ReadGroupProperties& rg_props);
    
    void finalize_rg_props();
    
    // TODO: We want to keep a collection of these, indexed by RG ID.  See #180
    ReadGroupProperties _rg_props; 
    
private:

    
	ReadTable& _insert_table;
	RefSequenceTable& _ref_table;

};

/******************************************************************************
 SAMHitFactory turns SAM alignments into ReadHits
*******************************************************************************/
class SAMHitFactory : public HitFactory
{
public:
	SAMHitFactory(const string& hit_file_name, 
				  ReadTable& insert_table, 
				  RefSequenceTable& reference_table) : 
		HitFactory(insert_table, reference_table), 
		_line_num(0), 
		_curr_pos(0) 
	{
		_hit_file = fopen(hit_file_name.c_str(), "r");
		if (_hit_file == NULL)
		{
			throw std::runtime_error("Error: could not open file for reading");
		}
        
        if (inspect_header() == false)
        {
            throw std::runtime_error("Error: could not parse SAM header");
        }
        
        // Override header-inferred read group properities with whatever
        // the user supplied.
        if (global_read_properties != NULL)
        {
            _rg_props = *global_read_properties;
        }
	}
	
	~SAMHitFactory() 
	{
		if (_hit_file)
		{
			fclose(_hit_file);
		}
	}
	
	virtual void undo_hit() 
	{ 
		fseeko(_hit_file, _curr_pos, SEEK_SET); 
		--_line_num;
	}
	
	void reset() { rewind(_hit_file); }
	
	void mark_curr_pos() { _curr_pos = ftell(_hit_file); }
	bool records_remain() const { return !feof(_hit_file); }
	
	bool next_record(const char*& buf, size_t& buf_size);
	
	bool get_hit_from_buf(const char* bwt_buf, 
						  ReadHit& bh,
						  bool strip_slash,
						  char* name_out = NULL,
						  char* name_tags = NULL);
    
    bool inspect_header();
    
private:
	static const size_t _hit_buf_max_sz = 10 * 1024;
	char _hit_buf[_hit_buf_max_sz];
	int _line_num;
	
	FILE* _hit_file;
	off_t _curr_pos;
};

/******************************************************************************
 BAMHitFactory turns SAM alignments into ReadHits
 *******************************************************************************/
class BAMHitFactory : public HitFactory
{
public:
	BAMHitFactory(const string& hit_file_name, 
				  ReadTable& insert_table, 
				  RefSequenceTable& reference_table) : 
		HitFactory(insert_table, reference_table) 
	{
		_hit_file = samopen(hit_file_name.c_str(), "rb", 0);
        
        memset(&_next_hit, 0, sizeof(_next_hit));
        
		if (_hit_file == NULL || _hit_file->header == NULL) 
		{
			throw std::runtime_error("Fail to open BAM file");
		}
		
		_beginning = bgzf_tell(_hit_file->x.bam);
        _eof_encountered = false;
        
        if (inspect_header() == false)
        {
            throw std::runtime_error("Error: could not parse BAM header");
        }
        
        // Override header-inferred read group properities with whatever
        // the user supplied.
        if (global_read_properties != NULL)
        {
            _rg_props = *global_read_properties;
        }
	}	
	
	~BAMHitFactory() 
	{
		if (_hit_file)
		{
			samclose(_hit_file);
		}
	}
	
	void mark_curr_pos() 
	{ 
		_curr_pos = bgzf_tell(_hit_file->x.bam);
	}

	
	void undo_hit() 
	{ 
		bgzf_seek(_hit_file->x.bam, _curr_pos, SEEK_SET);
		//--_line_num;
	}
    
    bool records_remain() const 
    { 
        return !_eof_encountered;
    }
	
	void reset() 
	{ 
		if (_hit_file && _hit_file->x.bam)
		{
			bgzf_seek(_hit_file->x.bam, _beginning, SEEK_SET);
            _eof_encountered = false;
		}
	}
	
	
	bool next_record(const char*& buf, size_t& buf_size);
	
	bool get_hit_from_buf(const char* bwt_buf, 
						  ReadHit& bh,
						  bool strip_slash,
						  char* name_out = NULL,
						  char* name_tags = NULL);
    
    bool inspect_header();
	
private:
	samfile_t* _hit_file; 
	int64_t _curr_pos;
	int64_t _beginning;
	
	bam1_t _next_hit; 
    bool _eof_encountered;
};

// Forward declaration of BundleFactory, because MateHit will need a pointer
// back to the Factory that created.  Ultimately, we should replace this
// with a pointer back to the ReadGroupProperty object corresponding to each 
// MateHit.  That, however, requires that we link fragment length distributions
// and bias models, etc, with each read group, and that's more than we can 
// afford to implement right now.

/*******************************************************************************
 MateHit is a class that encapsulates a paired-end alignment as a single object.
 MateHits can be "open" when one hit has been read from a stream of individual
 read alignments, but the other hasn't.  A "closed" MateHit is one where either
 both read alignments have been installed in the MateHit, or one read hit has,
 but the other will never come (i.e. singletons)
*******************************************************************************/
class MateHit
{
public:
    MateHit() : 
    _refid(0), 
	_collapse_mass(0.0) {}
    
	MateHit(shared_ptr<ReadGroupProperties const> rg_props,
            uint32_t refid, 
			shared_ptr<ReadHit const> left_alignment, 
			shared_ptr<ReadHit const> right_alignment) : 
    _rg_props(rg_props),
	_refid(refid), 
	_left_alignment(left_alignment),
	_right_alignment(right_alignment),
	_collapse_mass(0.0)
	{
		//_expected_inner_dist = min(genomic_inner_dist(), _expected_inner_dist);
	}
	~MateHit()
	{
		//fprintf(stderr, "Killing hit %lx\n",this);
	}

	//bool closed() {return _closed;}
	
    shared_ptr<ReadGroupProperties const> read_group_props() const { return _rg_props; }
    
	shared_ptr<ReadHit const> left_alignment() const {return _left_alignment;}
	void left_alignment(shared_ptr<ReadHit const> left_alignment) 
	{
		_left_alignment = left_alignment;
	}
	
	shared_ptr<ReadHit const> right_alignment() const {return _right_alignment;}					
	void right_alignment(shared_ptr<ReadHit const> right_alignment)  
	{
		_right_alignment = right_alignment;
	}
	
	bool is_pair() const
	{
		return (_left_alignment && _right_alignment);
	}
	

	int left() const 
	{
		if (_right_alignment && _left_alignment)
		{
			return min(_right_alignment->left(),_left_alignment->left());
		}
		if (_left_alignment)
			return _left_alignment->left();
		else if (_right_alignment)
			return _right_alignment->left(); 
		return -1;
	}
	
	int right() const 
	{
		if (_right_alignment && _left_alignment)
		{
			return max(_right_alignment->right(),_left_alignment->right());
		}
		if (_right_alignment)
			return _right_alignment->right();
		else if (_left_alignment)
			return _left_alignment->right(); 
		return -1;
	}
	
	CuffStrand strand() const 
	{
		CuffStrand left_strand = CUFF_STRAND_UNKNOWN;
		CuffStrand right_strand = CUFF_STRAND_UNKNOWN;
		if (_left_alignment)
		{
			left_strand = _left_alignment->source_strand();
		}
		if (_right_alignment)
		{
			right_strand = _right_alignment->source_strand();
			//assert ( s != CUFF_STRAND_UNKNOWN ? s == r : true);
		}
		assert (left_strand == right_strand || 
				left_strand == CUFF_STRAND_UNKNOWN || 
				right_strand == CUFF_STRAND_UNKNOWN);
		
		return max(left_strand, right_strand);
	}
	
	
	InsertID insert_id() const
	{
		if (_left_alignment) return _left_alignment->insert_id();
		if (_right_alignment) return _right_alignment->insert_id();
		return 0;
	}
	
	RefID ref_id() const { return _refid; }
	
	int genomic_inner_dist() const 
	{
		if (_left_alignment && _right_alignment)
		{
			return _right_alignment->left() - _left_alignment->right();
		}
		else
		{
			return -1;
		}
		return -1;
	}
	
	pair<int,int> genomic_outer_span() const
	{
		if (_left_alignment && _right_alignment)
		{
			return make_pair(left(),
							 right() - 1);
		}

		return make_pair(-1,-1);
	}	
	
	pair<int,int> genomic_inner_span() const 
	{
		if (_left_alignment && _right_alignment)
		{
			return make_pair(_left_alignment->right(),
							 _right_alignment->left() - 1);
		}

		return make_pair(-1,-1);
	}
	
	double error_prob() const
	{
		if (_left_alignment)
			return _left_alignment->error_prob();
		else if (_right_alignment)
			return _right_alignment->error_prob();
		return 1.0;
	}
	
	unsigned int  edit_dist() const
	{
		unsigned int edits = 0;
		if (_left_alignment)
			edits += _left_alignment->edit_dist();
		if (_right_alignment)
			edits += _right_alignment->edit_dist();
		return edits;
	}
	
	double mass() const { return (1.0 - error_prob()); }
	double collapse_mass() const { return _collapse_mass; }
	void collapse_mass(double m) { _collapse_mass = m; }
	void incr_collapse_mass(double incr) { _collapse_mass += incr; }
	
private:
	
    shared_ptr<ReadGroupProperties const> _rg_props;
	RefID _refid;
	shared_ptr<ReadHit const> _left_alignment;
	shared_ptr<ReadHit const> _right_alignment;
	double _collapse_mass;
	//bool _closed;
};

bool mate_hit_lt(const MateHit& lhs, const MateHit& rhs);

bool hits_eq_mod_id(const ReadHit& lhs, const ReadHit& rhs);

bool hits_equals(const MateHit& lhs, const MateHit& rhs);

// Assumes hits are sorted by mate_hit_lt
void collapse_hits(const vector<MateHit>& hits,
				   vector<MateHit>& non_redundant);



#endif
