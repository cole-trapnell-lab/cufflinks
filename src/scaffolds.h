#ifndef SCAFFOLDS_H
#define SCAFFOLDS_H
/*
 *  scaffolds.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/30/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <set>
#include <vector>

#include "common.h"
#include "hits.h"

#include <boost/thread.hpp>

using namespace std;

enum CuffOpCode { CUFF_MATCH, CUFF_INTRON, CUFF_UNKNOWN };

struct AugmentedCuffOp 
{
	AugmentedCuffOp(const CuffOpCode& O, int g_off, int g_len) 
	: opcode(O), 
	genomic_offset(g_off),
	genomic_length(g_len) 
	{
		assert (genomic_length >= 0);
	}
	
	void g_left(int left) 
	{ 
		int right = genomic_offset + genomic_length;
		genomic_offset = left;
		genomic_length = right - left;
	}
	int g_left() const { return genomic_offset; }
	
	void g_right(int right) 
	{ 
		genomic_length = right - genomic_offset;
	}
	int g_right() const { return genomic_offset + genomic_length; }
	
	static bool overlap_in_genome(const AugmentedCuffOp& lhs,
								  const AugmentedCuffOp& rhs)
	{	
		if (lhs.g_left() >= rhs.g_left() && lhs.g_left() < rhs.g_right())
			return true;
		if (lhs.g_right() > rhs.g_left() && lhs.g_right() < rhs.g_right())
			return true;
		if (rhs.g_left() >= lhs.g_left() && rhs.g_left() < lhs.g_right())
			return true;
		if (rhs.g_right() > lhs.g_left() && rhs.g_right() < lhs.g_right())
			return true;
		return false;
	}
	
	bool contains(const AugmentedCuffOp& other) const
	{
		if (g_left() <= other.g_left() && g_right() >= other.g_right())
			return true;
		return false;
	}
	
	bool properly_contains(const AugmentedCuffOp& other) const
	{
		if ((g_left() < other.g_left() && g_right() >= other.g_right()) ||
			(g_left() <= other.g_left() && g_right() > other.g_right()))
			return true;
		return false;
	}
	
	static int match_length(const AugmentedCuffOp& op, int left, int right)
	{
		int len = 0;
		int left_off = op.g_left();
		if (op.opcode == CUFF_MATCH)
		{
			if (left_off + op.genomic_length > left && left_off < right)
			{
				if (left_off > left)
				{
					if (left_off + op.genomic_length <= right + 1)
						len += op.genomic_length;
					else
						len += right - left_off;
				}
				else
				{
					if (left_off + op.genomic_length <= right + 1)
						len += (left_off + op.genomic_length - left);
					else
						return right - left;
				}
			}
		}
		return len;
	}
	
	static bool compatible(const AugmentedCuffOp& lhs,
						   const AugmentedCuffOp& rhs,
						   int overhang_tolerance = bowtie_overhang_tolerance);
	
	bool operator==(const AugmentedCuffOp& rhs) const
	{
		return (opcode == rhs.opcode && 
				genomic_offset == rhs.genomic_offset &&
				genomic_length == rhs.genomic_length);
	}
	
	bool operator<(const AugmentedCuffOp& rhs) const
	{
		if (opcode != rhs.opcode)
		{
			return opcode < rhs.opcode;
		}
		if (genomic_offset != rhs.genomic_offset)
		{
			return genomic_offset < rhs.genomic_offset;
		}
		if (genomic_length != rhs.genomic_length)
		{
			return genomic_length < rhs.genomic_length;
		}
		return false;
	}
    
    static bool g_left_lt(const AugmentedCuffOp& lhs,
						  const AugmentedCuffOp& rhs);
	
	bool operator!=(const AugmentedCuffOp& rhs) const
	{
		return !(*this == rhs);
	}
    
    static void merge_ops(const std::vector<AugmentedCuffOp>& ops, 
						  vector<AugmentedCuffOp>& merged,
						  bool introns_overwrite_matches,
                          bool allow_flank_introns = false);
    
    static void fill_interstices(vector<AugmentedCuffOp>& to_fill,
								 const vector<AugmentedCuffOp>& filler,
								 bool allow_flank_fill,
                                 bool allow_flank_introns);
	
	CuffOpCode opcode;
	int genomic_offset;
	int genomic_length;
};

class Scaffold
{
	void cuff_ops_from_cigar(vector<AugmentedCuffOp>& ops,
							 const vector<CigarOp>& cig,
							 int& g_left)
	{
		for (size_t i = 0; i < cig.size(); ++i)
		{
			assert(cig[i].length >= 0);
			switch(cig[i].opcode)
			{
				case MATCH:
					ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, cig[i].length));
					g_left += cig[i].length;
					break;
				case REF_SKIP:
					ops.push_back(AugmentedCuffOp(CUFF_INTRON, g_left, cig[i].length));
					g_left += cig[i].length;
					break;
				case SOFT_CLIP:
					g_left += cig[i].length;
					break;
                case HARD_CLIP:
					//g_left += cig[i].length;
					break;
                case INS:
                    //ops.back().genomic_length -= cig[i].length;
                    //g_left -= cig[i].length;
					break;
                case DEL:
					if (!ops.empty())
						ops.back().genomic_length += cig[i].length;
                    g_left += cig[i].length;
					break;
				default:
					assert(false);
					break;
			}
            assert (ops.empty() || ops.back().genomic_length >= 1);
		}
	}

	RefID _ref_id;
	
public:
	
	Scaffold() :
		_ref_id(0), 
		_is_ref(false),
		_strand(CUFF_STRAND_UNKNOWN), 
		_classcode(0),
        _fpkm(0.0) {}
	
	Scaffold(const MateHit& mate) :
		_ref_id(mate.ref_id()),
		_is_ref(false),
	    _classcode(0)
	{
		const ReadHit* left_hit = mate.left_alignment();
		//CuffAlign a;
		_strand = mate.strand();
		
		vector<AugmentedCuffOp> aug_ops;
		
		if (left_hit)
		{
			const vector<CigarOp>& l_cig = left_hit->cigar();
			int g_left = left_hit->left();
			cuff_ops_from_cigar(aug_ops, l_cig, g_left);
			
			const ReadHit* right_hit = mate.right_alignment();
			if (right_hit)
			{
				const vector<CigarOp>& r_cig = right_hit->cigar();
				int gap = (right_hit->left() - g_left);
				
				if (gap < 0)
				{
					g_left += gap;
					cuff_ops_from_cigar(aug_ops, r_cig, g_left);
				}
				else 
				{
					if (gap > 0)
					{
						//if (gap < (int)min_intron_length)
						//	aug_ops.push_back(AugmentedCuffOp(CUFF_MATCH, g_left, gap));
						//else
                        aug_ops.push_back(AugmentedCuffOp(CUFF_UNKNOWN, g_left, gap));
						g_left += gap;
					
					}
					cuff_ops_from_cigar(aug_ops, r_cig, g_left);
				}
				
			}
		}
		else
		{
			assert(false);
		}
		_mates_in_scaff.push_back(&mate);
		sort(aug_ops.begin(), aug_ops.end(), AugmentedCuffOp::g_left_lt);

        for(size_t i = 0; i < aug_ops.size(); ++i)
        {
            assert (aug_ops[i].genomic_length >= 1);
        }
        
        AugmentedCuffOp::merge_ops(aug_ops, _augmented_ops, false);
		
		int r_check = left();
		for (size_t i = 0; i < _augmented_ops.size(); ++i)
			r_check += _augmented_ops[i].genomic_length;
		
#ifdef DEBUG
		if (r_check != right())
		{
            AugmentedCuffOp::merge_ops(aug_ops, _augmented_ops, false);
		}
#endif
		assert (r_check == right());
		
		_has_intron = has_intron(*this);
		
		assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);

        for(size_t i = 1; i < _augmented_ops.size(); ++i)
		{
			assert(_augmented_ops[i-1].g_right() == _augmented_ops[i].g_left());
		}
		
		assert (_augmented_ops.front().opcode == CUFF_MATCH);
        assert (_augmented_ops.back().opcode == CUFF_MATCH);
        
        if (library_type == "transfrags")
        {
            double avg_fpkm = mate.mass();
            fpkm(avg_fpkm);
        }
	}
	
	Scaffold(const vector<Scaffold>& hits, bool introns_overwrite_matches = true) 
		:  _is_ref(false), _classcode(0)
	{
		assert (!hits.empty());
		_ref_id = hits[0].ref_id();
		
		Scaffold::merge(hits, *this, introns_overwrite_matches);
		
		assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);
        
        assert (_augmented_ops.front().opcode == CUFF_MATCH);
        assert (_augmented_ops.back().opcode == CUFF_MATCH);
        
        if (library_type == "transfrags")
        {
            double avg_fpkm = 0.0;
            foreach (const Scaffold& sc, hits)
            {
                avg_fpkm += sc.fpkm();
            }
            avg_fpkm /= hits.size();
            fpkm(avg_fpkm);
        }
	}
	
	// For manually constructing scaffolds, for example when a reference is 
	// available
	Scaffold(RefID ref_id, CuffStrand strand, const vector<AugmentedCuffOp>& ops, bool is_ref = false)
	: _ref_id(ref_id), 
	  _augmented_ops(ops), 
	  _strand(strand),
	  _classcode(0)
	{
		_has_intron = has_intron(*this);
		_is_ref = is_ref;
		
		assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);

        assert (_augmented_ops.front().opcode == CUFF_MATCH);
        assert (_augmented_ops.back().opcode == CUFF_MATCH);
	}
	
	//int get_id() const { return _id; }
		
	int left() const { return _augmented_ops.front().g_left(); }
	int right() const  { return _augmented_ops.back().g_right(); }
	
	const vector<const MateHit*>& mate_hits() const { return _mates_in_scaff; }
	
	RefID ref_id() const { return _ref_id; }
	void ref_id(RefID rid) { _ref_id = rid; }
	
	const string& annotated_trans_id() const { return _annotated_trans_id; }
	void annotated_trans_id(const string& ann_name) { _annotated_trans_id = ann_name; }
	
	const string& annotated_gene_id() const { return _annotated_gene_id; }
	void annotated_gene_id(const string& ann_name) { _annotated_gene_id = ann_name; }
	
	const string& annotated_gene_name() const { return _annotated_gene_name; }
	void annotated_gene_name(const string& ann_name) { _annotated_gene_name = ann_name; }
	
	const string& annotated_protein_id() const { return _annotated_protein_id; }
	void annotated_protein_id(const string& ann_name) { _annotated_protein_id = ann_name; }	
	
	const string& annotated_tss_id() const { return _annotated_tss_id; }
	void annotated_tss_id(const string& ann_name) { _annotated_tss_id = ann_name; }	
	
	const string& nearest_ref_id() const { return _nearest_ref_id; }
	void nearest_ref_id(const string& ann_name) { _nearest_ref_id = ann_name; }	
	
	double fpkm() const {return _fpkm; }
	void fpkm(double fpkm) { _fpkm = fpkm; }
 
	const string& seq() const { return _seq; } 
	void seq(const string& s) {	_seq = s; } 
	
	double gc_content() const
	{
		if (_seq != "")
		{
			int count = 0;
			for(size_t i = 0; i < _seq.length(); ++i)
			{
				if (_seq[i] == 'G' or _seq[i] == 'g' or _seq[i] == 'C' or _seq[i] == 'c')
					count ++;
			}
			return count/double(_seq.length());
		}
		return -1.0;
	}
	
	char nearest_ref_classcode() const { return _classcode; }
	void nearest_ref_classcode(char cc) { _classcode = cc; }
	
	bool has_intron() const { return _has_intron; }
	bool has_suspicious_unknown() const { return has_suspicious_unknown(*this); }

    // returns the fraction coverage of internal exons, returns 0 if no internal exons
	double internal_exon_coverage() const;
    
    // returns true if the scaffold strand is supported with reads or exon overlap with
    // a reference scaffold of known strand (since the scaffold may have been created with faux reads)
    bool has_strand_support(vector<shared_ptr<Scaffold> >* ref_scaffs = NULL) const;
    
    // returns true if all introns are supported with splice reads, false ow
    bool hits_support_introns() const; 
    bool hits_support_introns(set<AugmentedCuffOp>& hit_introns) const; 

    // returns true if all internal exons are fully covered and hits support introns, false ow
    bool has_struct_support(set<AugmentedCuffOp>& hit_introns) const;
    
	bool is_ref() const { return _is_ref; }
	void is_ref(bool ir) { _is_ref = ir; }
	
	CuffStrand strand() const 
    { 
        return _strand; 
    }
    
	void strand(CuffStrand strand) 
    { 
		assert(!has_strand_support() || _strand != CUFF_STRAND_UNKNOWN);
        _strand = strand; 
    }
	
	// Could we merge lhs and rhs together?
	static bool compatible(const Scaffold& lhs, 
						   const Scaffold& rhs,
						   int overhang_tolerance = bowtie_overhang_tolerance);
	
	static bool strand_agree(const Scaffold& lhs, 
							 const Scaffold& rhs);
	
	static bool exons_overlap(const Scaffold& lhs,
							  const Scaffold& rhs);
	
	// Incorporate Scaffold chow into this one.
	static void merge(const Scaffold& lhs, 
					  const Scaffold& rhs, 
					  Scaffold& merged,
					  bool introns_overwrite_matches);
	
	static void merge(const vector<Scaffold>& s, 
					  Scaffold& merged,
					  bool introns_overwrite_matches);

	// Extend 5' end using beginning of other scaffold without adding new exons.
	void extend_5(const Scaffold& other);
	// Clip final 3' exon by given amount
	void trim_3(int to_remove);
    
    // Extend final 3' exon by given amount
    void extend_3(int to_add);

	void tile_with_scaffs(vector<Scaffold>& tile_scaffs, int max_len, int tile_offset) const;

	// Creates a scaffold that matches this one but only covers the section from g_left for
	// a distance of match_length.  It is assumed that this region is contained in the scaffold.
	// sub_scaff should be an empty Scaffold object.
	bool sub_scaffold(Scaffold& sub_scaff, int g_left, int match_length) const;

 	// Tests whether the other scaffold is contained allowing the given overhang
	bool contains(const Scaffold& other, int ohang_5 = 0, int ohang_3 = 0) const
	{
		if (left() <= other.left() && right()>= other.right())
			return true;
		
        if (!(ohang_5 || ohang_3))
            return false;
        
		int left_hang;
		int right_hang;
		switch(strand())
		{
			case CUFF_FWD:
				left_hang = ohang_5;
				right_hang = ohang_3;
				break;
			case CUFF_REV:
				left_hang = ohang_3;
				right_hang = ohang_5;
				break;
			default:
				left_hang = max(ohang_3, ohang_5);
				right_hang = left_hang;
		}
				
		// Test to see if it is contained within the relaxed boundaries
		if ((left()-left_hang) <= other.left() && (right() + right_hang) >= other.right())
		{
			// Ensure that there are no exons outside of the strict boundaries
			return (other.augmented_ops().front().g_right() > left() && other.augmented_ops().back().g_left() < right()); 
		}
		
		return false;

	}
	
	// Tests whether the other scaffold contains the 5' end and is contained (allowing some overhang) on the 3' end
	// There can be no additional exons on the 5' end
	bool overlapped_3(const Scaffold& other, int ohang_5 = 0, int ohang_3 = 0) const
	{
		switch(strand())
		{
			case CUFF_FWD:
				return((left() + ohang_5 >= other.left() && right() + ohang_3 >= other.right()) && other.augmented_ops().front().g_right() > left());
			case CUFF_REV:
				return ((right() - ohang_5 <= other.right() && left() - ohang_3 <= other.left()) && other.augmented_ops().back().g_left() < right());
			default:
				return false;
		}
	}
	
	int match_length(int left, int right) const;
	
	int length() const
	{
		
		if(_seq != "")
			return _seq.length();
		
		int len = 0;
		
		// FIXME: this estimate really should include estimates of the CUFF_UNKNOWN lengths
		// for better abundance estimation.

		for (size_t j = 0; j < _augmented_ops.size(); ++j)
		{
			if (_augmented_ops[j].opcode == CUFF_MATCH)
				len += _augmented_ops[j].genomic_length; 
		}

		
		return len;
	}
	
	//void augmented_ops(vector<AugmentedCuffOp>& augmented) const;
	
	// The "score" of a merge is the log odds of lhs and rhs coming from
	// different transcripts.  We want to minimize the total score of 
	// all of our merges.
	static double score_merge(const Scaffold& lhs, const Scaffold& rhs);
	
	double score() const;
	
	double worst_mate_score() const;
	
	pair<int,int> genomic_to_transcript_span(pair<int,int> g_span) const;
	int genomic_to_transcript_coord(int g_coord) const;
	bool map_frag(const MateHit& hit, int& start, int& end, int& frag_len) const;

	
	static bool g_left_lt(const AugmentedCuffOp& lhs,
						  const AugmentedCuffOp& rhs);
	
	const vector<AugmentedCuffOp>& augmented_ops() const { return _augmented_ops; }
	void augmented_ops(vector<AugmentedCuffOp> aug_ops) { _augmented_ops = aug_ops; }

	
	static bool overlap_in_genome(const Scaffold& lhs, 
								  const Scaffold& rhs, 
								  int overlap_radius);
	
	vector<pair< int, int> > gaps() const
	{
		vector<pair<int,int> > g;
		const vector<AugmentedCuffOp>& ops = augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			if (ops[j].opcode == CUFF_INTRON)
			{
				g.push_back(make_pair(ops[j].g_left(), ops[j].g_right()));
			}
		}
		return g;
	}
	
	inline bool has_unknown() const
	{
		const vector<AugmentedCuffOp>& ops = augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			if (ops[j].opcode == CUFF_UNKNOWN)
				return true;
		}
		return false;
	}
	
    // Fills in CUFF_UNKNOWNs up to filled_gap_size long
	void fill_gaps(int filled_gap_size);
    
    // Fills in CUFF_UNKNOWNs with the contents of filler.  Filler must be
    // a sortted, contiguous, non-overlapping vector of AugmentedCuffOps
    void fill_gaps(const vector<AugmentedCuffOp>& filler);
	
	void clear_hits();
	bool add_hit(const MateHit*);
	
	void get_complete_subscaffolds(vector<Scaffold>& complete);
private: 
	
	void initialize_exon_lists();
	
	static bool has_intron(const Scaffold& scaff)
	{
		
		const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			if (ops[j].opcode == CUFF_INTRON)
				return true;
		}
		
		return false;
	}
	
	static bool has_suspicious_unknown(const Scaffold& scaff)
	{
		
		const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
		for (size_t j = 0; j < ops.size(); ++j)
		{
			if (ops[j].opcode == CUFF_UNKNOWN && 
				ops[j].genomic_length > max_frag_len)
				return true;
		}
		
		return false;
	}
	
	static double score_overlap(const AugmentedCuffOp& op, 
								int inner_left_edge, 
								int inner_right_edge);
	
	static double score_contigs(const Scaffold& contig,
								const vector<pair<int, int> >& inners);
	
	static double min_score(const Scaffold& contig, 
							const vector<pair<int, int> >& inners);
	
	static bool compatible_contigs(const Scaffold& lhs, 
											const Scaffold& rhs,
											int overhang_tolerance = bowtie_overhang_tolerance);
	
	
	typedef vector<AugmentedCuffOp> OpList;
	
	bool check_merge_length(const vector<AugmentedCuffOp>& ops);
	
	static void fill_interstices(vector<AugmentedCuffOp>& to_fill,
								 const vector<AugmentedCuffOp>& filler,
								 bool allow_flank_fill);
	
	static void merge_ops(const std::vector<AugmentedCuffOp>& ops, 
						  vector<AugmentedCuffOp>& merged,
						  bool introns_overwrite_matches);
	
	vector<const MateHit*> _mates_in_scaff;
	
	bool _has_intron; 
	bool _is_ref;
	
	vector<AugmentedCuffOp> _augmented_ops;
	CuffStrand _strand;
	
	string _annotated_trans_id;
	string _annotated_gene_id;
	string _annotated_gene_name;
	string _annotated_protein_id;
	string _annotated_tss_id;
	string _nearest_ref_id;
	char _classcode;
	
	string _seq;
	double _fpkm;

};

bool scaff_lt(const Scaffold& lhs, const Scaffold& rhs);
bool scaff_lt_rt(const Scaffold& lhs, const Scaffold& rhs);
bool scaff_lt_rt_oplt(const Scaffold& lhs, const Scaffold& rhs);
bool scaff_lt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs);
bool scaff_lt_rt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs);
bool scaff_lt_rt_oplt_sp(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs);


bool overlap_in_genome(int ll, int lr, int rl, int rr);

struct StructurallyEqualScaffolds
{
	bool operator()(shared_ptr<Scaffold> lhs, shared_ptr<Scaffold> rhs)
	{
		return lhs->ref_id() == rhs->ref_id() && 
		lhs->augmented_ops() == rhs->augmented_ops();
	}
};

#endif
