#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

/*
 *  biascorrection.h
 *  cufflinks
 *
 *  Created by Adam Roberts on 5/20/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "bundles.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <list>
#include <string>


namespace ublas = boost::numeric::ublas;
using namespace std;

void get_compatibility_list(const vector<Scaffold>& transcripts,
							const vector<MateHit>& alignments,
							vector<list<int> >& compatibilities);


class BiasLearner{
	
	static const int pow4[];
	static const int paramTypes[];
	static const int MAX_SLICE;
	static const int CENTER;
	static const int _M;
	static const int _N;
	
	static const int lengthBins[];
	static const double positionBins[];
	
	EmpDist _frag_len_dist;
	ublas::matrix<long double> _startParams;
	ublas::matrix<long double> _startExp;
	ublas::matrix<long double> _endParams;
	ublas::matrix<long double> _endExp;
	ublas::matrix<long double> _posParams;
	ublas::matrix<long double> _posExp;
	vector<long double> posParam_sums;

	int seqToInt(const char* seqSlice, int n);
	void getSlice(const char* seq, char* slice, int start, int end);
	void genNList(const char* seqSlice, int start, int n, list<int>& nList);
	
public:
	
	BiasLearner(const EmpDist& frag_len_dist);
	void processTranscript(const vector<long double>& startHist, const vector<long double>& endHist, const Scaffold& transcript);
	void normalizeParameters();
	void output();
	
	void getBias(const Scaffold& transcript, vector<double>& startBiases, vector<double>& endBiases);

};

bool learn_bias(BundleFactory& bundle_factory, BiasLearner& bl);

#endif
