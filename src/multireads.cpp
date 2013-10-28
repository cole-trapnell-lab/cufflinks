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

void MultiRead::add_hit(RefID r_id, int left, int right)
{
	_hits.push_back(MultiHit(r_id, left, right));
}

MultiHit* MultiRead::get_hit(RefID r_id, int left, int right)
{
    for (size_t i = 0; i < num_hits(); ++i)
	{
		MultiHit& hit = _hits[_curr_index];
		if (hit.r_id == r_id && hit.left == left && hit.right == right)
		{
			return &hit;
        }
		_curr_index = (_curr_index + 1) % num_hits();
	}  
    fprintf(stderr, "\nWARNING: Multi-Hit not found (%d,%d).\n", left, right); 
    return NULL;
}

//allele
void MultiRead::add_expr(RefID r_id, int left, int right, double expr, double paternal_expr, double maternal_expr)
{
	MultiHit* hit = get_hit(r_id, left, right);
    if (hit)
    {
        hit->expr += expr;
        _tot_expr += expr;
		//allele
		hit->paternal_expr += paternal_expr;
		hit->maternal_expr += maternal_expr;
		_tot_paternal_expr += paternal_expr;
		_tot_maternal_expr += maternal_expr;
    }
}
double MultiRead::get_mass(RefID r_id, int left, int right, bool valid_mass)
{
	if (!valid_mass)
	{
		return 1.0/num_hits();
	}
	
	if (_tot_expr == 0.0)
		return 0.0;
	
	MultiHit* hit = get_hit(r_id, left, right);
    if (hit)
        return hit->expr/_tot_expr;
    else
        return 1.0/num_hits();
}

//allele
void MultiRead::get_mass(RefID r_id, int left, int right, bool valid_mass, double& paternal_mass, double& maternal_mass)
{
	if (!valid_mass)
	{
		paternal_mass = 0.5/num_hits();
		maternal_mass = 0.5/num_hits();
	}
	else if(_tot_paternal_expr == 0.0 && _tot_maternal_expr == 0.0)
	{
		paternal_mass = 0.0;
		maternal_mass = 0.0;
	}
	else
	{
		MultiHit* hit = get_hit(r_id, left, right);
		if (hit)
		{
			paternal_mass = hit->paternal_expr/_tot_expr;
			maternal_mass = hit->maternal_expr/_tot_expr;
		}
		else
		{
			paternal_mass = 0.5/num_hits();
			maternal_mass = 0.5/num_hits();
		}
	}
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
	//allele
	double num_paternal_hits,num_maternal_hits;
	int num_hits = hit.num_hits();
	if(hit.allele() == ALLELE_PATERNAL)
	{
		num_paternal_hits += 1.0*hit.num_hits();
		
	}
	else if(hit.allele() == ALLELE_MATERNAL)
	{
		num_maternal_hits += 1.0*hit.num_hits();
		
	}
	else
	{
		num_paternal_hits = num_hits/2.0;
		num_paternal_hits = num_hits/2.0;
	}
	//allele
	add_hit(hit.ref_id(), hit.left(), hit.right(), hit.insert_id(), hit.num_hits(), num_paternal_hits, num_maternal_hits);
}
//allele
void MultiReadTable::add_hit(RefID r_id, int left, int right, InsertID mr_id, int exp_num_hits, double exp_num_paternal_hits, double exp_num_maternal_hits)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif
	MultiRead* mr = get_read(mr_id);
	if (!mr)
	{
		//allele
		mr = &((_read_map.insert(std::make_pair(mr_id, MultiRead(mr_id, exp_num_hits, exp_num_paternal_hits, exp_num_maternal_hits)))).first->second); 
	}
	mr->add_hit(r_id, left, right);
}
//allele
void MultiReadTable::add_expr(const MateHit& hit, double expr, double paternal_expr, double maternal_expr)
{
	//allele
		if(paternal_expr == 0 && maternal_expr == 0)
	{
		paternal_expr = expr/2.0;
		maternal_expr = expr/2.0;
	}
		//allele
	add_expr(hit.ref_id(), hit.left(), hit.right(), hit.insert_id(), expr, paternal_expr, maternal_expr);
}
//allele
void MultiReadTable::add_expr(RefID r_id, int left, int right, InsertID mr_id, double expr, double paternal_expr, double maternal_expr)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif
	MultiRead* mr = get_read(mr_id);
    if (mr)
		//allele
        mr->add_expr(r_id, left, right, expr, paternal_expr, maternal_expr);
}

double MultiReadTable::get_mass(const MateHit& hit)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif	
	MultiRead* mr = get_read(hit.insert_id());
	if(!mr)
		return 1.0;
	return mr->get_mass(hit.ref_id(), hit.left(), hit.right(), _valid_mass);
}

void MultiReadTable::get_mass(const MateHit& hit,double& paternal_mass, double& maternal_mass)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_lock);
#endif	
	MultiRead* mr = get_read(hit.insert_id());
	if(!mr)
	{
		paternal_mass = 0.5;
		maternal_mass = 0.5;
	}
	else
	{
		mr->get_mass(hit.ref_id(), hit.left(), hit.right(), _valid_mass, paternal_mass, maternal_mass);
	}
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
