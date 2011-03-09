#ifndef MULTIREADS_H
#define MULTIREADS_H

#include <boost/thread.hpp>

typedef uint64_t RefID;
typedef uint64_t InsertID;

class MultiRead
{
	typedef std::pair<RefID, int> MultiHit ;
	size_t _curr_index;
	std::vector<MultiHit> _hits;
	std::vector<double> _expr;
	double _tot_expr;
	InsertID _id;
	bool* _valid_mass;

	
public:
	
	MultiRead(InsertID id, bool* valid_mass, int exp_num_hits)
	:	_curr_index(0),
		_tot_expr(0.0),
		_id(id),
		_valid_mass(valid_mass)
	{
		_hits.reserve(exp_num_hits);
		_expr.reserve(exp_num_hits);
	}
	
	size_t num_hits() { return (int)_hits.size(); }
	void add_hit(RefID r_id, int left);
	void add_expr(RefID r_id, int left, double expr);
	double get_mass(RefID r_id, int left);
	

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
	void add_hit(RefID r_id, int left, InsertID mr_id, int exp_num_hits);
	void add_expr(const MateHit& hit, double expr);
	void add_expr(RefID r_id, int left, InsertID mr_id, double expr);
	double get_mass(const MateHit& hit);
	size_t num_multireads();
	size_t num_multihits();
	
};

#endif
