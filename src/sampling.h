#ifndef SAMPLING_H
#define SAMPLING_H
//
//  sampling.h
//  cufflinks
//
//  Created by Cole Trapnell on 12/19/11.
//  Copyright 2011 Cole Trapnell. All rights reserved.
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "common.h"

#include <stdint.h>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>

//#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/constants/constants.hpp>

// Boost Cholesky factorizations in the spirit of lu.hpp
// Written by Robbie Vogt, found at: 
// http://lists.boost.org/MailArchives/ublas/2005/07/0568.php

namespace boost { namespace numeric { namespace ublas {
	
	// Cholesky factorization
	template<class M>
	double cholesky_factorize (M &m) 
	{
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::value_type value_type;
		
		BOOST_UBLAS_CHECK (m.size1() == m.size2(), external_logic("Cholesky decomposition is only valid for a square, positive definite matrix."));
		
		size_type size = m.size1();
		vector<value_type> d(size);
		//bool positive_definite = true;
		for (size_type i = 0; i < size; ++ i) {
			matrix_row<M> mri (row (m, i));
			for (size_type j = i; j < size; ++ j) {
				matrix_row<M> mrj (row (m, j));
				
				value_type elem = m(i,j) - inner_prod(project(mri,range(0,i)), project(mrj,range(0,i)));
				
				if (i == j) {
					if (elem <= 0.0) {
						// matrix after rounding errors is not positive definite
						return elem;
					}
					else {
						d(i) = sqrtl(elem);
					}
				}
				else {
					m(j,i) = elem / d(i);
				}
			}
		}
		
		// put the diagonal back in
		for (size_type i = 0; i < size; ++ i) {
			m(i,i) = d(i);
		}
		
		//cerr << m << endl;
		for (size_type i = 0; i < size; ++i) {
			for (size_type j = 0; j < i; ++j)
			{
				m(j,i) = 0;
			}
		}
		//cerr << m << endl;
		// decomposition succeeded
		return 0.0;
	}
	
	
	// Cholesky substitution 
	template<class M, class E> 
	void cholesky_substitute (const M &m, vector_expression<E> &e) { 
		typedef const M const_matrix_type; 
		typedef vector<typename E::value_type> vector_type; 
		inplace_solve (m, e, lower_tag ()); 
		inplace_solve (trans(m), e, upper_tag ()); 
	} 
	template<class M, class E> 
	void cholesky_substitute (const M &m, matrix_expression<E> &e) { 
		typedef const M const_matrix_type; 
		typedef matrix<typename E::value_type> matrix_type; 
		inplace_solve (m, e, lower_tag ()); 
		inplace_solve (trans(m), e, upper_tag ()); 
	} 
	template<class E, class M> 
	void cholesky_substitute_left (vector_expression<E> &e, const M &m) { 
		typedef const M const_matrix_type; 
		typedef vector<typename E::value_type> vector_type; 
		inplace_solve (trans(m), e, upper_tag ()); 
		inplace_solve (m, e, lower_tag ()); 
	} 
	template<class E, class M> 
	void cholesky_substitute_left (matrix_expression<E> &e, const M &m) { 
		typedef const M const_matrix_type; 
		typedef matrix<typename E::value_type> matrix_type; 
		inplace_solve (trans(m), e, upper_tag ()); 
		inplace_solve (m, e, lower_tag ()); 
	} 
	// Cholesky matrix inversion 
	template<class M> 
	void cholesky_invert (M &m) 
	{ 
		typedef typename M::size_type size_type; 
		typedef typename M::value_type value_type; 
		size_type size = m.size1(); 
		// determine the inverse of the lower traingular matrix 
		for (size_type i = 0; i < size; ++ i) { 
			m(i,i) = 1 / m(i,i); 
			for (size_type j = i+1; j < size; ++ j) { 
				value_type elem(0); 
				for (size_type k = i; k < j; ++ k) { 
					elem -= m(j,k)*m(k,i); 
				} 
				m(j,i) = elem / m(j,j); 
			} 
		} 
		// multiply the upper and lower inverses together 
		m = prod(trans(triangular_adaptor<M,lower>(m)), triangular_adaptor<M,lower>(m)); 
	} 
	
	
}}}

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool lu_invert_matrix (const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());
	
 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
	if( res != 0 ) return false;
	
 	// create identity matrix of "inverse"
 	inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));
	
 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
	
 	return true;
}

///* Matrix inversion routine.
// Expects input to be PRE-FACTORIZED */
template<class T>
bool chol_invert_matrix (const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse) {
	
 	using namespace boost::numeric::ublas;
	inverse = input;
	
	cholesky_invert(inverse);
	
 	return true;
}


// Adapted for Boost from Numerical Recipes
template<typename ValueType>
class multinormal_generator
{
	typedef boost::mt19937 base_generator_type;
	typedef boost::normal_distribution<> distribution_type;
    
public:
	// expects the mean vector and the *CHOLESKY* factorization of the covariance
	multinormal_generator(const boost::numeric::ublas::vector<ValueType>& mean,
						  const boost::numeric::ublas::matrix<ValueType>& chol_cov)
    :	
    _engine(random_seed),
    _distribution(),
    _generator(boost::variate_generator<base_generator_type&, 
               distribution_type >(_engine, 
                                   _distribution))
	{
		_rand = boost::numeric::ublas::zero_vector<ValueType>(mean.size());
		_mean = mean; 
		_cholesky = chol_cov;
	}
	
	const boost::numeric::ublas::vector<ValueType>& next_rand()
	{
		boost::numeric::ublas::vector<ValueType> temp(_mean.size());
		for (size_t i = 0; i < _mean.size(); ++i)
		{
			double r = _generator();
			temp(i) = r;
			_rand(i) = 0.0;
		}
		
		//cerr << "rand ="<<temp << endl;
		//_rand = prod(ublas::triangular_adaptor<ublas::matrix<ValueType>,ublas::lower>(_cholesky), temp);
		for (size_t i = 0; i < _cholesky.size1(); ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
				_rand(i) += _cholesky(i,j) * temp(j);
			}
		}
		//cerr <<_rand << " + " << _mean << "=";
		_rand = _rand + _mean;
		//cerr <<_rand <<endl;
		
		return _rand;
	}
    
    void set_parameters(const boost::numeric::ublas::vector<ValueType>& mean,
                        const boost::numeric::ublas::matrix<ValueType>& chol_cov)
    {
        _rand = boost::numeric::ublas::zero_vector<ValueType>(mean.size());
		_mean = mean; 
		_cholesky = chol_cov;
    }
	
private:
	boost::numeric::ublas::vector<ValueType>		_rand;
	boost::numeric::ublas::vector<ValueType>		_mean;
	boost::numeric::ublas::matrix<ValueType>		_cholesky;
	
	base_generator_type								_engine;
	distribution_type								_distribution;
	boost::variate_generator<base_generator_type&, 
    distribution_type>		_generator;
};

// expects a cholesky factorized covariance matrix
template<class matrix_T>
double determinant(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_r)
{
	double det = 1.0;
	
	matrix_T chol(mat_r());
	
	for (size_t i = 0; i < chol.size1(); ++i)
	{
		det *= chol(i,i);
	}
	
	return det * det;
}


// Given log(p) and log(q) returns log(p+q)double
template<class float_type>
float_type log_space_add(float_type log_p, float_type log_q)
{
	if (log_p < log_q)
	{
		float_type tmp = log_p;
		log_p = log_q;
		log_q = tmp;
	}
    
	assert (log_p >= log_q);
	return log (1.0 + exp(log_q - log_p)) + log_p;
}

void generate_importance_samples(multinormal_generator<double>& generator,
                                 std::vector<boost::numeric::ublas::vector<double> >& samples, 
                                 int num_samples,
                                 bool no_zeros = true);

#endif
