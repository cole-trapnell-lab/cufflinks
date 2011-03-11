/*
 *  multireads.cpp
 *  cufflinks
 *
 *  Created by Adam Roberts on 3/6/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */


#include "hits.h"
#include "multireads.h"

void MultiRead::add_hit(RefID r_id, int left)
{
	_hits.push_back(std::make_pair(r_id, left));
	_expr.push_back(0.0);
}

void MultiRead::add_expr(RefID r_id, int left, double expr)
{
	while (true)
	{
		MultiHit& hit = _hits[_curr_index];
		if (hit.first == r_id && hit.second == left)
		{
			_expr[_curr_index] += expr;
			break;
		}
		_curr_index = (_curr_index + 1) % num_hits();
	}
	_tot_expr += expr;
}

double MultiRead::get_mass(RefID r_id, int left, bool valid_mass)
{
	if (!valid_mass)
	{
		return 1.0/num_hits();
	}
	
	if (_tot_expr == 0.0)
		return 0.0;
	
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

MultiRead* MultiReadTable::get_read(InsertID mr_id)
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

void MultiReadTable::add_hit(const MateHit& hit)
{
	add_hit(hit.ref_id(), hit.left(), hit.insert_id(), hit.num_hits());
}

void MultiReadTable::add_hit(RefID r_id, int left, InsertID mr_id, int exp_num_hits)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif
	MultiRead* mr = get_read(mr_id);
	if (!mr)
	{
		mr = &((_read_map.insert(std::make_pair(mr_id, MultiRead(mr_id, exp_num_hits)))).first->second); 
	}
	mr->add_hit(r_id, left);
}

void MultiReadTable::add_expr(const MateHit& hit, double expr)
{
	add_expr(hit.ref_id(), hit.left(), hit.insert_id(), expr);
}

void MultiReadTable::add_expr(RefID r_id, int left, InsertID mr_id, double expr)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif
	MultiRead* mr = get_read(mr_id);
	mr->add_expr(r_id, left, expr);
}

double MultiReadTable::get_mass(const MateHit& hit)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif	
	MultiRead* mr = get_read(hit.insert_id());
	if(!mr)
		return 1.0;
	return mr->get_mass(hit.ref_id(), hit.left(), _valid_mass);
}

size_t MultiReadTable::num_multireads()
{
	return (int)_read_map.size();
}

size_t MultiReadTable::num_multihits()
{
	size_t count = 0;
	for (MultiReadMap::iterator it=_read_map.begin() ; it != _read_map.end(); it++ )
	{
		count += it->second.num_hits();
	}
	return count;
}
