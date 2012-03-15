/*
 *  jensen_shannon.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/30/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <numeric>

#include "jensen_shannon.h"

using namespace std;
using namespace boost;


double entropy(const Eigen::VectorXd& p)
{
	double e = 0;  
	for (size_t i = 0; i < p.size(); ++i)
	{
		double P = p[i];
        if (P != 0.0)
        {
            e -= (P * log(P));
        }
	}
    if (e < 0)
        e = 0;
	return e;
}

double jensen_shannon_distance(std::vector<Eigen::VectorXd>& sample_kappas)
{
	assert (sample_kappas.size() > 1);
	
	size_t kappa_length = 0;
	for (size_t i = 1; i < sample_kappas.size(); ++i)
	{
		assert (sample_kappas[i].size() == sample_kappas[i-1].size());
		kappa_length = sample_kappas[i].size();
	}
	
    Eigen::VectorXd avg_kappas = Eigen::VectorXd::Zero(kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		//cout << "kappa " << i<< " "<< sample_kappas[i] << endl;
		avg_kappas += sample_kappas[i];
	}
	avg_kappas /= sample_kappas.size();
	//cout << avg_kappas << endl;
	
	double avg_entropy = 0.0;
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		avg_entropy += entropy(sample_kappas[i]);
	}
	avg_entropy /= sample_kappas.size();
	//cout << avg_entropy << endl;
	
	double entropy_avg = entropy(avg_kappas);
	
	double js = entropy_avg - avg_entropy;
    if (js < 0) // can happen due to underflow or rounding errors.
        return 0;
	return sqrt(js);
}

void jensen_shannon_gradient(vector<Eigen::VectorXd>& sample_kappas,
							 double js,
							 ublas::vector<double>& gradient)
{
	assert (sample_kappas.size() > 1);
	size_t kappa_length = sample_kappas.front().size();
	for (size_t i = 1; i < sample_kappas.size(); ++i)
	{
		assert (sample_kappas[i].size() == sample_kappas[i-1].size());
		kappa_length = sample_kappas[i].size();
	}
	
	if (kappa_length == 0)
		return;
	
	gradient = ublas::zero_vector<double>(sample_kappas.size() * kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		for (size_t k = 0; k < kappa_length; ++k)
		{
            assert (!isinf(sample_kappas[i](k)) && !isnan(sample_kappas[i](k)));
			gradient(i*kappa_length + k) = sample_kappas[i](k);
		}
	}
	
	//cout << "t1: " << gradient<< endl;
	
    ublas::vector<double> p_bar = ublas::zero_vector<double>(kappa_length);
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		for (size_t j = 0; j < sample_kappas[i].size(); ++j)
        {
            p_bar[j] += sample_kappas[i](j);   
        }
	}
	p_bar /= sample_kappas.size();
	
	
	//cout << "t2 " << denoms << endl;
	
	for (size_t i = 0; i < sample_kappas.size(); ++i)
	{
		for (size_t k = 0; k < kappa_length; ++k)
		{
			if (p_bar(k) == 0.0 || gradient(i*kappa_length + k) == 0.0 || js == 0.0)
			{
				gradient(i*kappa_length + k) = 0.0;
			}
			else
			{
                double alt_grad = 0.0;
                double m = 2.0;
                alt_grad = js / (2.0 * m);
                double A = log(gradient(i*kappa_length + k)) + (1.0 / gradient(i*kappa_length + k));
                double B = log(p_bar[k]) + (1.0 / p_bar[k]);
                alt_grad *= (A - B);
                
				gradient(i*kappa_length + k) /= p_bar(k);
				gradient(i*kappa_length + k) = log(gradient(i*kappa_length + k));
				gradient(i*kappa_length + k) /= sample_kappas.size(); // m in paper notation
				gradient(i*kappa_length + k) *= (1.0/(2.0 * js)); // This is supposed to use the square root of the distance (it's not a typo)
				
                double curr_grad = gradient(i*kappa_length + k);
                
                assert (!isinf(curr_grad) && !isnan(curr_grad));
                //fprintf(stderr, "Curr gradient: %lg, alternate gradient %lg\n", curr_grad, alt_grad);
            
			}
		}
	}
}

void make_js_covariance_matrix(vector<ublas::matrix<double> >& kappa_covariances,
							   ublas::matrix<double>& js_covariance)
{
	size_t kappa_length = 0;
	for (size_t i = 1; i < kappa_covariances.size(); ++i)
	{
		assert (kappa_covariances[i].size1() == kappa_covariances[i-1].size1());
		assert (kappa_covariances[i].size2() == kappa_covariances[i-1].size2());
		
		kappa_length = kappa_covariances[i].size1();
	}
	
	if (kappa_length == 0)
		return;
	
	js_covariance = ublas::zero_matrix<double>(kappa_covariances.size() * kappa_length,
											   kappa_covariances.size() * kappa_length);
	for (size_t i = 0; i < kappa_covariances.size(); ++i)
	{
		for (size_t j = 0; j < kappa_length; ++j)
		{
			for (size_t k = 0; k < kappa_length; ++k)
			{
				js_covariance(i*kappa_length + j, i*kappa_length + k) = 
				kappa_covariances[i](j,k);
                assert (!isinf(js_covariance(i*kappa_length + j, i*kappa_length + k) && ! isnan(js_covariance(i*kappa_length + j, i*kappa_length + k))));
			}
		}
	}
}
