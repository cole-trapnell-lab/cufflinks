#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

/*
 *  biascorrection.h
 *  cufflinks
 *
 *  Created by Adam Roberts on 5/20/10.
 *  Copyright 2010 Adam Roberts. All rights reserved.
 *
 */


#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <list>
#include <string>
#include <boost/tr1/unordered_map.hpp>
#include <boost/thread.hpp>
#include "common.h"

class MateHit;
class Scaffold;
class BundleFactory;
class HitBundle;

namespace ublas = boost::numeric::ublas;

void get_compatibility_list(const std::vector<Scaffold>& transcripts,
							const std::vector<MateHit>& alignments,
							std::vector<std::list<int> >& compatibilities);

class BiasLearner{
	static const int pow4[];
	static const int MAX_SLICE;
	static const int CENTER;
	static const int _m;
	static const int _n;
	
	static const int lengthBins[];
	static const double positionBins[];
	static const int siteSpec[];
	static const int vlmmSpec[];
	
	const int* paramTypes;
	boost::shared_ptr<EmpDist const> _frag_len_dist;
	ublas::matrix<long double> _startSeqParams;
	ublas::matrix<long double> _startSeqExp;
	ublas::matrix<long double> _endSeqParams;
	ublas::matrix<long double> _endSeqExp;
	ublas::matrix<long double> _startPosParams;
	ublas::matrix<long double> _startPosExp;
	ublas::matrix<long double> _endPosParams;
	ublas::matrix<long double> _endPosExp;

	int seqToInt(const char* seqSlice, int n) const;
	void getSlice(const char* seq, char* slice, int start, int end) const;
	void genNList(const char* seqSlice, int start, int n, std::list<int>& nList) const;

#if ENABLE_THREADS	
	boost::mutex _bl_lock;
#endif
	
public:
	
	BiasLearner(boost::shared_ptr<EmpDist const> frag_len_dist);
	void preProcessTranscript(const Scaffold& transcript);
	
	void processTranscript(const std::vector<double>& startHist, const std::vector<double>& endHist, const Scaffold& transcript);
	void normalizeParameters();
	void output();
	
	void getBias(const Scaffold& transcript, std::vector<double>& startBiases, std::vector<double>& endBiases) const;

};

void learn_bias(BundleFactory& bundle_factory, BiasLearner& bl, bool progress_bar = true);
void process_bundle(HitBundle& bundle, BiasLearner& bl);

// Helps with the complexities of bias correction with replicates in cond_probs and eff_lens
class BiasCorrectionHelper{
	
	boost::shared_ptr<Scaffold> _transcript;
    boost::unordered_map<boost::shared_ptr<ReadGroupProperties const>, int> _rg_index;
	int _size;
	bool _mapped;
	
	std::vector<std::vector<double> > _start_biases;
	std::vector<std::vector<double> > _end_biases;
	std::vector<std::vector<double> > _pos_biases;
	std::vector<std::vector<double> > _tot_biases_for_len;
	std::vector<std::vector<double> > _start_biases_for_len;
	std::vector<std::vector<double> > _end_biases_for_len;
	
	std::vector<double> _eff_lens;
	std::vector<double> _rg_masses;
	
	int add_read_group(boost::shared_ptr<ReadGroupProperties const> rgp);	
	int get_index(boost::shared_ptr<ReadGroupProperties const> rgp);
	
public:
	
	BiasCorrectionHelper(boost::shared_ptr<Scaffold> transcript) 
	{ 
		_transcript = transcript;
		_mapped = false;
		_size = 0; 
	}
	

	double get_cond_prob(const MateHit& hit);
	
	double get_effective_length();
	bool is_mapped() { return _mapped; }

};

#endif
