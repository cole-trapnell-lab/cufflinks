/*
 *  jensen_shannon.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 8/30/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>

namespace ublas = boost::numeric::ublas;

double entropy(const ublas::vector<double>& p);

double jensen_shannon_distance(std::vector<ublas::vector<double> >& sample_kappas);

//void alt_jensen_shannon_gradient(std::vector<ublas::vector<double> >& sample_kappas,
//                                 double js,
//                                 ublas::vector<double>& gradient);

void jensen_shannon_gradient(std::vector<ublas::vector<double> >& sample_kappas,
							 double js,
							 ublas::vector<double>& gradient);

void make_js_covariance_matrix(std::vector<ublas::matrix<double> >& kappa_covariances,
							   ublas::matrix<double>& js_covariance);
