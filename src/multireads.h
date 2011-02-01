#ifndef MULTIREADS_H
#define MULTIREADS_H

#include "hits.h"

class MultiRead
{
	typedef pair<RefID, int> MultiHit ;
	size_t _curr_index;
	vector<MultiHit > _hits;
	vector<double> _expr;
	double _tot_expr;
	InsertID _id;
	
public:
	
	MultiRead(InsertID id, int exp_number=0)
	:	_curr_index(0),
		_tot_expr(0.0),
		_id(id)
	{
		if (exp_number > 0)
		{
			_hits.reserve(exp_number);
			_expr.reserve(exp_number);
		}
	}
	
	size_t num_hits() { return (int)_hits.size(); }
	
	void add_hit(RefID r_id, int left)
	{
		_hits.push_back(make_pair(r_id, left));
		_expr.push_back(0.0);
	}
	
	void add_expr(RefID r_id, int left, double expr)
	{
		while (true)
		{
			MultiHit hit = _hits[_curr_index];
			if (hit.first == r_id && hit.second == left)
			{
				_expr[_curr_index] += expr;
				break;
			}
			_curr_index = (_curr_index + 1) % num_hits();
		}
		_tot_expr += expr;
	}
	
	double get_mass(RefID r_id, int left)
	{
		if (_tot_expr == 0)
		{
			return 1.0/num_hits();
		}
		
		double this_expr;
		while (true)
		{
			MultiHit& hit = _hits[_curr_index];
			if (hit.first == r_id && hit.second == left)
			{
				this_expr = _expr[_curr_index];
				break;
			}
			_curr_index = (_curr_index + 1) % num_hits();
		}
		
		return this_expr/_tot_expr;
	}
	

};


class MultiReadTable
{
	typedef map<InsertID, MultiRead> MultiReadMap;
	MultiReadMap _read_map;
	
public:
	MultiRead* add_hit(RefID r_id, int left, InsertID mr_id, int exp_number=0)
	{
		MultiRead* mr = get_read(mr_id);
		if (!mr)
		{
			mr = &((_read_map.insert(make_pair(mr_id, MultiRead(mr_id, exp_number)))).first->second); 
		}
		mr->add_hit(r_id, left);
		return mr;
	}
										
	MultiRead* get_read(InsertID mr_id)
	{
		MultiReadMap::iterator it;
		it = _read_map.find(mr_id);
		if (it == _read_map.end())
		{
			return NULL;
		}
		else 
		{
			return &(it->second);
		}
	}
};

#endif
