#ifndef MULTIREADS_H
#define MULTIREADS_H

#include <boost/thread.hpp>

typedef uint64_t RefID;
typedef uint64_t InsertID;

struct MultiHit
{
    MultiHit(RefID id, int l, int r)
    :   r_id(id),
        left(l),
        right(r),
        expr(0),
		//allele
		paternal_expr(0),
		maternal_expr(0){}
    RefID r_id;
    int left;
    int right;
    double expr;
	//allele
	double paternal_expr;
	double maternal_expr;
};

class MultiRead
{
	size_t _curr_index;
	std::vector<MultiHit> _hits;
	double _tot_expr;
	//allele
	double _tot_paternal_expr;
	double _tot_maternal_expr;
	InsertID _id;
    
    MultiHit* get_hit(RefID r_id, int left, int right);
	
public:
//allele		
	MultiRead(InsertID id, int exp_num_hits, double exp_num_paternal_hits, double exp_num_maternal_hits)
	:	_curr_index(0),
		_tot_expr(0.0),
		//allele
		_tot_paternal_expr(0.0),
		_tot_maternal_expr(0.0),
		_id(id)
	{
		_hits.reserve(exp_num_hits);
	}
	
	size_t num_hits() { return (int)_hits.size(); }
	void add_hit(RefID r_id, int left, int right);
//allele
	void add_expr(RefID r_id, int left, int right, double expr, double paternal_expr, double maternal_expr);
	double get_mass(RefID r_id, int left, int right, bool valid_mass);
	//allele
	void get_mass(RefID r_id, int left, int right, bool valid_mass, double& paternal_mass, double& maternal_mass);
};

class MateHit;

class MultiReadTable
{
	typedef std::map<InsertID, MultiRead> MultiReadMap;
	MultiReadMap _read_map;
	bool _valid_mass;
	MultiRead* get_read(InsertID mr_id);
#if ENABLE_THREADS
	boost::mutex _lock;
#endif
public:
	MultiReadTable(): _valid_mass(false) {}
	
	void valid_mass(bool vm) { _valid_mass = vm; }
	void add_hit(const MateHit& hit);
	//allele
	void add_hit(RefID r_id, int left, int right, InsertID mr_id, int exp_num_hits, double exp_num_paternal_hits, double exp_num_maternal_hits);
	void add_expr(const MateHit& hit, double expr, double paternal_expr = 0, double maternal_expr = 0);
	void add_expr(RefID r_id, int left, int right, InsertID mr_id, double expr, double paternal_expr, double maternal_expr);
	double get_mass(const MateHit& hit);
	//allele: get_mass(const MateHit& hit) is insensitive to allele info and will return the combined parental mass
	//get_mass(const MateHit& hit, double& paternal_mass, double& maternal_mass) is allele sensitive and will return each parental mass separately
	void get_mass(const MateHit& hit, double& paternal_mass, double& maternal_mass);
	size_t num_multireads();
	size_t num_multihits();
	
};

#endif
