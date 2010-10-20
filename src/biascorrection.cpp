/*
 *  biascorrection.cpp
 *  cufflinks
 *
 *  Created by Adam Roberts on 5/20/10.
 *  Copyright 2010 Adam Roberts. All rights reserved.
 *
 */

#include "biascorrection.h"
#include "abundances.h"
#include "progressbar.h"

#include <iostream>
#include <fstream>

void output_vector(vector<double>& v, char* fname)
{
	ofstream myfile1;
	string filename = output_dir + "/" + fname;
	myfile1.open (filename.c_str());
	
	for (size_t i = 0; i < v.size(); ++i)
		myfile1 << v[i] <<",";
	myfile1 << endl;
	
	myfile1.close();
}

double colSums(const ublas::matrix<long double>& A, vector<long double>& sums)
{
	long double total = 0.0;
	sums = vector<long double>(A.size2(),0.0);
	for (size_t i = 0; i < A.size1(); ++i)
		for (size_t j = 0; j < A.size2(); ++j)
		{
			sums[j] += A(i,j);
			total += A(i,j);
		}
	return total;
}

double fourSums(const ublas::matrix<long double>& A, ublas::matrix<long double>& sums)
{
	long double total = 0.0;
	sums = ublas::zero_matrix<long double>(A.size1(), A.size2()/4);
	for (size_t i = 0; i < A.size1(); ++i)
		for (size_t j = 0; j < A.size2(); ++j)
		{
			sums(i,j/4) += A(i,j);
			total += A(i,j);
		}
	
	return total;
}

void ones(ublas::matrix<long double>& A)
{
	for (size_t i = 0; i < A.size1(); ++i)
		for (size_t j = 0; j < A.size2(); ++j)
			A(i,j) = 1;
}

void get_compatibility_list(const vector<shared_ptr<Scaffold> >& transcripts,
                            const vector<MateHit>& alignments,
                            vector<list<int> >& compatibilities)
{
	int M = alignments.size();
	int N = transcripts.size();
	
	vector<Scaffold> alignment_scaffs;
	
	for (size_t i = 0; i < alignments.size(); ++i)
	{
		const MateHit& hit = alignments[i];
		alignment_scaffs.push_back(Scaffold(hit));
	} 
	
	for (int i = 0; i < M; ++i) 
	{
		for (int j = 0; j < N; ++j) 
		{
			if (transcripts[j]->strand() != CUFF_STRAND_UNKNOWN
				&& transcripts[j]->contains(alignment_scaffs[i]) 
				&& Scaffold::compatible(*transcripts[j],alignment_scaffs[i]))
			{
				compatibilities[i].push_back(j);
			}
		}
	}
}

void learn_bias(BundleFactory& bundle_factory, BiasLearner& bl)
{
	HitBundle bundle;
	RefSequenceTable& rt = bundle_factory.ref_table();

	ProgressBar p_bar("Learning bias parameters.", bundle_factory.num_bundles());

	while(true)
	{
		HitBundle* bundle_ptr = new HitBundle();
		
		if (!bundle_factory.next_bundle(*bundle_ptr))
		{
			delete bundle_ptr;
			break;
		}
		
		HitBundle& bundle = *bundle_ptr;
		
		char bundle_label_buf[2048];
		sprintf(bundle_label_buf, "%s:%d-%d", rt.get_name(bundle.ref_id()),	bundle.left(), bundle.right());
		p_bar.update(bundle_label_buf, 1);
		

		if (bundle.non_redundant_hits().size()==0)
		{
			delete bundle_ptr;
			continue;
		}

		process_bundle(bundle, bl);
			
		delete bundle_ptr;
	}
	p_bar.complete();
	bl.normalizeParameters();
#if ADAM_MODE
	bl.output();
#endif
}

void process_bundle(HitBundle& bundle, BiasLearner& bl)
{
	const vector<shared_ptr<Scaffold> >& transcripts = bundle.ref_scaffolds();
	const vector<MateHit>& nr_alignments = bundle.non_redundant_hits();
	
	int M = nr_alignments.size();
	int N = transcripts.size();
	
	vector<list<int> > compatibilities(M,list<int>());
	get_compatibility_list(transcripts, nr_alignments, compatibilities);	
	
	vector<vector<double> > startHists(N, vector<double>());
	vector<vector<double> > endHists(N, vector<double>());
	vector<double> useable_mass (N, 0.0);
	
	for (int j = 0; j < N; ++j)
	{
		int size = transcripts[j]->length();
		startHists[j].resize(size+1); // +1 catches overhangs
		endHists[j].resize(size+1);
	}
	
	for (int i = 0; i < M; ++i)
	{
		const MateHit& hit = nr_alignments[i];
		if (!hit.left_alignment() && !hit.right_alignment())
			continue;

		double num_hits = nr_alignments[i].collapse_mass();
		long double locus_fpkm = 0;
		
		for (list<int>::iterator it=compatibilities[i].begin(); it!=compatibilities[i].end(); ++it)
		{
			locus_fpkm += transcripts[*it]->fpkm(); 
		}
		
		if (locus_fpkm==0)
			continue;
		
		for (list<int>::iterator it=compatibilities[i].begin(); it!=compatibilities[i].end(); ++it)
		{	
			const Scaffold& transcript = *transcripts[*it];
			assert(transcript.strand()!=CUFF_STRAND_UNKNOWN); //Filtered in compatibility list
			
			int start;
			int end;
			int frag_len;
			
			double mass_fraction = num_hits*transcript.fpkm()/locus_fpkm;
			
			transcript.map_frag(hit, start, end, frag_len);
			
			if (start != transcript.length() || end != transcript.length())
				useable_mass[*it] += mass_fraction;
			
			startHists[*it][start] += mass_fraction;
			endHists[*it][end] += mass_fraction;
		}
	}
 	for (int j = 0; j < N; ++j)
 	{
 		if (transcripts[j]->strand()!=CUFF_STRAND_UNKNOWN && 
			transcripts[j]->seq()!="" &&  
//			useable_mass[j]/transcripts[j]->length() >= 0.1 &&
			transcripts[j]->fpkm() >= 1)
 			bl.processTranscript(startHists[j], endHists[j], *transcripts[j]);
 	}
}

const int BiasLearner::pow4[] = {1,4,16,64};
//const int BiasLearner::paramTypes[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
const int BiasLearner::paramTypes[] = {1,1,1,1,1,2,2,2,3,3,3,3,3,3,3,3,2,2,2,1,1}; //Length of connections at each position in the window
const int BiasLearner::MAX_SLICE = 3; // Maximum connection length
const int BiasLearner::CENTER = 8; //Index in paramTypes[] of first element in read
const int BiasLearner::_M = 21; //Number of positions spanned by window
const int BiasLearner::_N = 64; //Length of maximum connection in VLMM
const int BiasLearner::lengthBins[] = {791,1265,1707,2433}; //Quantiles derived from human mRNA length distribution in UCSC genome browser
const double BiasLearner::positionBins[] = {.02,.04,.06,.08,.10,.15,.2,.3,.4,.5,.6,.7,.8,.85,.9,.92,.94,.96,.98,1};
//const double BiasLearner::positionBins[] = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1};

BiasLearner::BiasLearner(shared_ptr<EmpDist const> frag_len_dist)
{
	_frag_len_dist = frag_len_dist;
	_startSeqParams = ublas::zero_matrix<long double>(_M,_N);
	_startSeqExp = ublas::zero_matrix<long double>(_M,_N);
	_endSeqParams = ublas::zero_matrix<long double>(_M,_N);
	_endSeqExp = ublas::zero_matrix<long double>(_M,_N);
	_startPosParams = ublas::zero_matrix<long double>(20,5);
	_startPosExp = ublas::zero_matrix<long double>(20,5);
	_endPosParams = ublas::zero_matrix<long double>(20,5);
	_endPosExp = ublas::zero_matrix<long double>(20,5);
}


inline int BiasLearner::seqToInt(const char* seqSlice, int n) const
{
	int c = 0;
	for(int i = 0; i < n; i++)
	{	
		if (seqSlice[i] == 4) return -1;//N
		c += (seqSlice[i])*pow4[n-i-1];
	}
	return c;
}

inline void BiasLearner::getSlice(const char* seq, char* slice, int start, int end) const// INCLUSIVE!
{
	if (end >= start)
	{
		for (int i = start; i <= end; ++i)
		{
			slice[i-start] = seq[i];
		}
	}
	else 
	{
		for(int i = start; i >= end; --i)
		{
			slice[start-i] = seq[i];
		}
	}
}


void BiasLearner::processTranscript(const std::vector<double>& startHist, const std::vector<double>& endHist, const Scaffold& transcript)
{
	double fpkm = transcript.fpkm();
	int seqLen = transcript.length();
	
	char seq[seqLen];
	char c_seq[seqLen];
	encode_seq(transcript.seq(), seq, c_seq);
	
	char seqSlice[MAX_SLICE];
	
	int lenClass=0;
	while (seqLen > lengthBins[lenClass] && lenClass < 4)
	{
		lenClass++;
	}
	
	// We want to only use the portion of the transcript where fragments can start/end
	int min_frag_len = _frag_len_dist->min();
	int currStartBin = 0;
	int startBinCutoff = positionBins[currStartBin]*(seqLen - min_frag_len);
	int currEndBin = 0;
	int endBinCutoff = positionBins[currStartBin]*(seqLen - min_frag_len);
		
	for (int i=0; i < seqLen; i++)
	{
		
		//Position Bias
		if (i > startBinCutoff && currStartBin < _startPosParams.size1()-1)
			startBinCutoff=positionBins[++currStartBin]*(seqLen - min_frag_len);
		if (i - min_frag_len > endBinCutoff)
			endBinCutoff = positionBins[++currEndBin]*(seqLen - min_frag_len);
			
		_startPosParams(currStartBin, lenClass) += startHist[i]/fpkm;
		_startPosExp(currStartBin, lenClass) += !(_frag_len_dist->too_short(seqLen-i));
		_endPosParams(currEndBin, lenClass) += endHist[i]/fpkm;
		_endPosExp(currEndBin, lenClass) += !(_frag_len_dist->too_short(i+1));

		
		bool start_in_bounds = i-CENTER >= 0 && i+(_M-1)-CENTER < seqLen;
		bool end_in_bounds = i+CENTER-(_M-1) >= 0 && i+CENTER < seqLen;
		
		if (!start_in_bounds && !end_in_bounds) // Make sure we are in bounds of the sequence
			continue;
		
		//Sequence Bias
		for(int j=0; j < _M; j++)
		{
			// Start Bias
			if (start_in_bounds) // Make sure we are in bounds of the sequence
			{
				int k = i+j-CENTER;
				getSlice(seq, seqSlice, k-(paramTypes[j]-1), k);
				int v = seqToInt(seqSlice,paramTypes[j]);
				if (v >= 0)
				{
					_startSeqParams(j,v) += startHist[i]/fpkm;
					_startSeqExp(j,v) += !(_frag_len_dist->too_short(seqLen-i));
				}
				else // There is an N.  Average over all possible values of N
				{
					list<int> nList(1,0);
					genNList(seqSlice, 0, paramTypes[j],nList);
					for (list<int>::iterator it=nList.begin(); it!=nList.end(); ++it)
					{
						_startSeqParams(j,*it) += startHist[i]/(fpkm * (double)nList.size());
						_startSeqExp(j,*it) += !(_frag_len_dist->too_short(seqLen-i))/(double)nList.size();
					}
				}
			}
			// End Bias
			if (end_in_bounds) // Make sure we are in bounds of the sequence
			{
				int k = i+CENTER-j;
				getSlice(c_seq, seqSlice, k+(paramTypes[j]-1), k); 
				int v = seqToInt(seqSlice, paramTypes[j]);
				if (v >= 0)
				{
					_endSeqParams(j,v) += endHist[i]/fpkm;
					_endSeqExp(j,v) += !(_frag_len_dist->too_short(seqLen-i));
				}
				else // There is an N.  Average over all possible values of N
				{
					list<int> nList(1,0);
					genNList(seqSlice, 0, paramTypes[j], nList);
					for (list<int>::iterator it=nList.begin(); it!=nList.end(); ++it)
					{
						_endSeqParams(j,*it) += endHist[i]/(fpkm * (double)nList.size());
						_endSeqExp(j,*it) += !(_frag_len_dist->too_short(seqLen-i))/(double)nList.size();
					}
				}
			}
		}
	}
}

void BiasLearner::getBias(const Scaffold& transcript, vector<double>& startBiases, vector<double>& endBiases) const
{
	if (transcript.seq()=="")
		return;
		
	int seqLen = transcript.length();
	
	char seq[seqLen];
	char c_seq[seqLen];
	encode_seq(transcript.seq(), seq, c_seq);
	
	char seqSlice[MAX_SLICE];
	
	int lenClass=0;
	while (seqLen > lengthBins[lenClass] && lenClass < 4)
	{
		lenClass++;
	}
	
	int min_frag_len = _frag_len_dist->min();
	int currStartBin = 0;
	int startBinCutoff = positionBins[currStartBin]*(seqLen - min_frag_len);
	int currEndBin = 0;
	int endBinCutoff = positionBins[currEndBin]*(seqLen - min_frag_len);
	
	for (int i=0; i < seqLen; i++)
	{
		//Position Bias
		if (i > startBinCutoff && currStartBin < _startPosParams.size1()-1)
			startBinCutoff=positionBins[++currStartBin]*(seqLen - min_frag_len);
		if (i - min_frag_len > endBinCutoff)
			endBinCutoff = positionBins[++currEndBin]*(seqLen - min_frag_len);

		double startBias = _startPosParams(currStartBin, lenClass);
		double endBias = _endPosParams(currEndBin,lenClass);
		
		//Sequence Bias
		
		bool start_in_bounds = i-CENTER >= 0 && i+(_M-1)-CENTER < seqLen;
		bool end_in_bounds = i+CENTER-(_M-1) >= 0 && i+CENTER < seqLen - _frag_len_dist->mean(); // don't count bias near end since we're over-counting these fragments
		
		if (start_in_bounds || end_in_bounds) // Make sure we are in bounds of the sequence
		{
			for(int j=0; j < _M; j++)
			{
				// Start Bias
				if (start_in_bounds) // Make sure we are in bounds of the sequence
				{
					int k = i+j-CENTER;
					getSlice(seq, seqSlice, k-(paramTypes[j]-1), k);
					int v = seqToInt(seqSlice, paramTypes[j]);
					if (v >= 0)
					{
						startBias *= _startSeqParams(j,v);
					}
					else // There is an N.  Average over all possible values of N
					{
						list<int> nList(1,0);
						double tot = 0;
						genNList(seqSlice, 0, paramTypes[j],nList);

						for (list<int>::iterator it=nList.begin(); it!=nList.end(); ++it)
						{
							tot += _startSeqParams(j,*it);
						}
						startBias *= tot/nList.size();
					}
				}
				
				// End Bias
				if (end_in_bounds) // Make sure we are in bounds of the sequence
				{
					int k = i+CENTER-j;
					getSlice(c_seq, seqSlice, k+(paramTypes[j]-1), k);
					int v = seqToInt(seqSlice,paramTypes[j]);
					if (v >= 0)
					{
						endBias *= _endSeqParams(j,v);
					}
					else // There is an N.  Average over all possible values of N
					{
						list<int> nList(1,0);
						double tot = 0;
						genNList(seqSlice, 0, paramTypes[j],nList);
						for (list<int>::iterator it=nList.begin(); it!=nList.end(); ++it)
						{
							tot += _endSeqParams(j,*it);
						}
						endBias *= tot/nList.size();
					}
				}
			}
		}
		startBiases[i] = startBias;
		endBiases[i] = endBias;
	}
}

void BiasLearner::genNList(const char* seqSlice, int start, int n, list<int>& nList) const
{

	if (n > 1) 
		genNList(seqSlice, start+1, n-1, nList);
	
	
	if (n==1 && seqSlice[start]==4)
	{
		for (int j=0; j<4; ++j) 
			nList.push_back(j);
	}
	else if (n==1)
	{
		nList.push_back(seqSlice[start]);
	}
	else if (seqSlice[start]==4)
	{
		for (int i = nList.size()-1; i>=0; --i)
		{
			for (int j=0; j<4; ++j) 
				nList.push_back(nList.front()+j*pow4[n-1]);
			nList.pop_front();
		}
	}
	else
	{
		for (list<int>::iterator it=nList.begin(); it!=nList.end(); ++it)
			(*it)+=seqSlice[start]*pow4[n-1];
	}
	
}

void BiasLearner::normalizeParameters()
{
	//Normalize position parameters	
	vector<long double> startPosParam_sums;
	vector<long double> startPosExp_sums;
	double start_tot = colSums(_startPosParams, startPosParam_sums); // Total starts for each length class
	colSums(_startPosExp, startPosExp_sums); // Total FPKM for each length class
	
	vector<long double> endPosParam_sums;
	vector<long double> endPosExp_sums;
	double end_tot = colSums(_endPosParams, endPosParam_sums); // Total starts for each length class
	colSums(_endPosExp, endPosExp_sums); // Total FPKM for each length class

	for(size_t i=0; i < _startPosParams.size1(); i++)
	{
		for(size_t j=0; j < _startPosParams.size2(); j++)
		{
			_startPosParams(i,j) /= startPosParam_sums[j];
			_startPosExp(i,j) /= startPosExp_sums[j];
			if (_startPosExp(i,j) == 0)
				_startPosParams(i,j) = numeric_limits<long double>::max();
			else
				_startPosParams(i,j) /= _startPosExp(i,j);
			
			_endPosParams(i,j) /= endPosParam_sums[j];
			_endPosExp(i,j) /= endPosExp_sums[j];
			if (_endPosExp(i,j) == 0)
				_endPosParams(i,j) = numeric_limits<long double>::max();
			else
				_endPosParams(i,j) /= _endPosExp(i,j);
		}
	}
	
	if (start_tot == 0.0)
		ones(_startPosParams);
	if (end_tot == 0.0)
		ones(_endPosParams);
	
	ublas::matrix<long double> startSeqExp_sums;
	ublas::matrix<long double> startParam_sums;
	start_tot = fourSums(_startSeqParams, startParam_sums);
	fourSums(_startSeqExp, startSeqExp_sums);
	
	ublas::matrix<long double> endSeqExp_sums;
	ublas::matrix<long double> endParam_sums;
	end_tot = fourSums(_endSeqParams, endParam_sums);
	fourSums(_endSeqExp, endSeqExp_sums); 
	
	//Normalize sequence parameters
	for(int i=0; i < _M; i++)
	{
		for(int j=0; j < pow4[paramTypes[i]]; j++)
		{
			if (startParam_sums(i,j/4)==0)
			{
				_startSeqParams(i,j) = 1;
			}
			else 
			{
				_startSeqParams(i,j) /= startParam_sums(i,j/4);
				_startSeqExp(i,j) /= startSeqExp_sums(i,j/4);
				_startSeqParams(i,j) /= _startSeqExp(i,j);
			}
			if (endParam_sums(i,j/4)==0)
			{
				_endSeqParams(i,j) = 1;
			}
			else 
			{
				_endSeqParams(i,j) /= endParam_sums(i,j/4);
				_endSeqExp(i,j) /= endSeqExp_sums(i,j/4);
				_endSeqParams(i,j) /= _endSeqExp(i,j);
			}

		}
	}
	
	if (start_tot==0.0)
		ones(_startSeqParams);
	if (end_tot==0.0)
		ones(_endSeqParams);
	
	if (bias_mode=="vlmm" || bias_mode=="site")
	{
		ones(_startPosParams);
		ones(_endPosParams);
	}
	else if (bias_mode =="pos")	
	{
		ones(_startSeqParams);
		ones(_endSeqParams);
	}
}

void BiasLearner::output()
{
	ofstream myfile1;
	ofstream myfile2;
	ofstream myfile3;
	ofstream myfile4;
	string startfile = output_dir + "/startSeqBias.csv";
	myfile1.open (startfile.c_str());
	string endfile = output_dir + "/endSeqBias.csv";
	myfile2.open (endfile.c_str());
	string posfile = output_dir + "/startPosBias.csv";
	myfile3.open (posfile.c_str());
	string posfile2 = output_dir + "/endPosBias.csv";
	myfile4.open (posfile2.c_str());


	for (int i = 0; i < _N; ++i)
	{
		for(int j = 0; j < _M; ++j)
		{
			myfile1 << _startSeqParams(j,i) <<",";
			myfile2 << _endSeqParams(j,i) << ",";
		}
		myfile1 << endl;
		myfile2 << endl;
	}
	
	for (size_t i = 0; i < _startPosParams.size2(); ++i)
	{
		for(size_t j = 0; j < _startPosParams.size1(); ++j)
		{
			myfile3 << _startPosParams(j,i) <<",";
			myfile4 << _endPosParams(j,i) <<",";
		}
		myfile3 <<endl;
		myfile4 <<endl;
	}
	
	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
}



int BiasCorrectionHelper::add_read_group(shared_ptr<ReadGroupProperties const> rgp)
{
	int trans_len = _transcript->length();
	_rg_index.insert(make_pair(rgp, _size));
	
	// Defaults are values for a run not using bias correction
	vector<double> start_bias(trans_len+1, 1.0);
	vector<double> end_bias(trans_len+1, 1.0);
	double eff_len = 0.0;
	
	if (rgp->bias_learner()!=NULL && _transcript->strand()!=CUFF_STRAND_UNKNOWN)
	{
		rgp->bias_learner()->getBias(*_transcript, start_bias, end_bias);
	}
	
	shared_ptr<EmpDist const> frag_len_dist = rgp->frag_len_dist();
	
	vector<double> tot_bias_for_len(trans_len+1, 0);
	vector<double> start_bias_for_len(trans_len+1, 0);
	vector<double> end_bias_for_len(trans_len+1, 0);

	tot_bias_for_len[trans_len] = trans_len;
	start_bias_for_len[trans_len] = trans_len;
	end_bias_for_len[trans_len] = trans_len;

	
	for(int l = rgp->frag_len_dist()->min(); l <= trans_len; l++)
	{
		//double tot = 0;
		//double start = 0;
		//double end = 0;
		for(int i = 0; i <= trans_len - l; i++)
		{
			double tot_bias = start_bias[i]*end_bias[i+l-1];
			tot_bias_for_len[l] += tot_bias;
			start_bias_for_len[l] += start_bias[i];
			end_bias_for_len[l] += end_bias[i+l-1];
			eff_len += tot_bias * frag_len_dist->npdf(l, trans_len - i);
		}
        //assert(tot != 0);
		//tot_bias_for_len[l] = tot;
		//eff_len += tot * p_len;
		//mean_start_bias += start * p_len;
		//mean_end_bias += end * p_len; 
	}
	
	_start_biases.push_back(start_bias);
	_end_biases.push_back(end_bias);
	_tot_biases_for_len.push_back(tot_bias_for_len);
	_eff_lens.push_back(eff_len);
	_start_biases_for_len.push_back(start_bias_for_len);
	_end_biases_for_len.push_back(end_bias_for_len);
	_rg_masses.push_back(0.0);
	
	return _size++; // Index of new element
}

int num_adds = 0;
int BiasCorrectionHelper::get_index(shared_ptr<ReadGroupProperties const> rgp)
{
	map<shared_ptr<ReadGroupProperties const>, int>::iterator iter;
	iter = _rg_index.find(rgp);
	
	if (iter==_rg_index.end()) //This rg is not yet in the index, so add it.
	{
        num_adds++;
		return add_read_group(rgp);
	}
	
	return iter->second;
}

// Hit needs to be collapsed
double BiasCorrectionHelper::get_cond_prob(const MateHit& hit)
{
	shared_ptr<ReadGroupProperties const> rgp = hit.read_group_props();
	
	int i = get_index(rgp);
	
	int start;
	int end;
	int frag_len;
	int trans_len = _transcript->length();
	
	_transcript->map_frag(hit, start, end, frag_len);

	shared_ptr<const EmpDist> fld = rgp->frag_len_dist();
	
	double cond_prob = 1.0;
	cond_prob *= _start_biases[i][start];
	cond_prob *= _end_biases[i][end];
	cond_prob *= fld->npdf(frag_len, trans_len-start); //defaults to pdf if trans_len==start
	
	if (cond_prob==0.0)
		return 0.0;
	
	if (hit.is_pair())
		cond_prob /= _tot_biases_for_len[i][frag_len];
	else if (start!=trans_len && end==trans_len) // The hit is a singleton at the start of a fragment
		cond_prob /= _start_biases_for_len[i][frag_len];
	else if (start==trans_len && end!=trans_len) // The hit is a singleton at the end of a fragment
		cond_prob /= _end_biases_for_len[i][frag_len];
	else if (frag_len==trans_len)  // We don't actually know where we start or end and can't subtract off the frag_len or we'll get inf
		cond_prob /= trans_len;
	else // Single-end read w/ library type FF or RR
		cond_prob /= trans_len-frag_len;

	if (cond_prob > 0 && hit.collapse_mass() > 0)
	{
		_rg_masses[i] += hit.collapse_mass();
		_mapped = true;
	}
	
	assert(!isinf(cond_prob));
	assert(!isnan(cond_prob));
	return cond_prob;
}

double BiasCorrectionHelper::get_effective_length()
{
	
	if (_size==0)
		return _transcript->length();

	
	double tot_mass = accumulate( _rg_masses.begin(), _rg_masses.end(), 0.0 );
	double eff_len = 0.0;
	
	if (tot_mass==0)
		return _transcript->length();
	
    for (map<shared_ptr<ReadGroupProperties const>, int>::iterator itr = _rg_index.begin();
         itr != _rg_index.end();
         ++itr)
	{
		int i = itr->second;
		double rg_eff_len = _eff_lens[i];
		eff_len += rg_eff_len * (_rg_masses[i]/tot_mass);
	}
	
	assert(eff_len>0);
	assert(!isnan(eff_len));
	return eff_len;
}
																		
									

																		
