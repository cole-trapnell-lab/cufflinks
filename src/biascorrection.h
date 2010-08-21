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
void map_frag_to_transcript(const Scaffold& transcript, 
							const MateHit& hit, 
							int& start, int& end, int& frag_len);


class BiasLearner{
	
	static const int pow4[];
	static const int paramTypes[];
	static const int MAX_SLICE;
	static const int CENTER;
	static const int _M;
	static const int _N;
	
	static const int lengthBins[];
	static const double positionBins[];
	
	shared_ptr<EmpDist const> _frag_len_dist;
	ublas::matrix<long double> _startParams;
	ublas::matrix<long double> _startExp;
	ublas::matrix<long double> _endParams;
	ublas::matrix<long double> _endExp;
	ublas::matrix<long double> _posParams;
	ublas::matrix<long double> _posExp;

	int seqToInt(const char* seqSlice, int n) const;
	void getSlice(const char* seq, char* slice, int start, int end) const;
	void genNList(const char* seqSlice, int start, int n, list<int>& nList) const;
	
public:
	
	BiasLearner(shared_ptr<EmpDist const> frag_len_dist);
	void processTranscript(const vector<long double>& startHist, const vector<long double>& endHist, const Scaffold& transcript);
	void normalizeParameters();
	void output();
	
	void getBias(const Scaffold& transcript, vector<double>& startBiases, vector<double>& endBiases, vector<double>& posBiases) const;

};

void learn_bias(BundleFactory& bundle_factory, BiasLearner& bl);
void process_bundle(HitBundle& bundle, BiasLearner& bl);

// Helps with the complexities of bias correction with replicates in cond_probs and eff_lens
class BiasCorrectionHelper{
	
	shared_ptr<Scaffold> _transcript;
	map<shared_ptr<ReadGroupProperties const>, int> _rg_index;
	int _size;
	bool _mapped;
	
	vector<vector<double> > _start_biases;
	vector<vector<double> > _end_biases;
	vector<vector<double> > _pos_biases;
	vector<vector<double> > _tot_biases_for_len;
	vector<double> _eff_lens;
	vector<double> _mean_start_biases;
	vector<double> _mean_end_biases;
	vector<double> _rg_masses;
	
	int add_read_group(shared_ptr<ReadGroupProperties const> rgp);	
	int get_index(shared_ptr<ReadGroupProperties const> rgp);
	
public:
	
	BiasCorrectionHelper(shared_ptr<Scaffold> transcript) 
	{ 
		if (transcript->annotated_trans_id()=="NM_000029")
			fprintf(stderr,"HERE\n");
		_transcript = transcript;
		_mapped = false;
		_size = 0; 
	}
	

	double get_cond_prob(const MateHit& hit);
	
	double get_effective_length();
	bool is_mapped() { return _mapped; }

};



#endif
