/*
 *  filters.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 10/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include "filters.h"
#include <algorithm>
#include <numeric>

using namespace std;

void filter_introns(int bundle_length,
					int bundle_left,
					vector<Scaffold>& hits, 
					double fraction,
					bool filter_on_intron_overlap,
					bool filter_with_intron_doc)
{
	vector<int> depth_of_coverage(bundle_length,0);
	vector<double> scaff_doc;
	map<pair<int,int>, int> intron_doc;
	vector<Scaffold> filtered_hits;
	vector<bool> toss(hits.size(), false);
	
	double bundle_avg_doc = compute_doc(bundle_left, 
										hits, 
										depth_of_coverage, 
										intron_doc,
										false);
	
	double bundle_avg_thresh = bundle_avg_doc * fraction;
	
	if (filter_with_intron_doc && !intron_doc.empty())
	{
		bundle_avg_doc = major_isoform_intron_doc(intron_doc);
		bundle_avg_thresh = fraction * bundle_avg_doc;
#if ASM_VERBOSE
		fprintf(stderr, "\tFiltering bundle introns, avg (intron) doc = %lf, thresh = %f\n", bundle_avg_doc, bundle_avg_thresh);
#endif
	}
	else
	{
#if ASM_VERBOSE
		fprintf(stderr, "\tFiltering bundle introns, avg bundle doc = %lf, thresh = %f\n", bundle_avg_doc, bundle_avg_thresh);
#endif	
	}
	
	for(map<pair<int, int>, int>::const_iterator itr = intron_doc.begin();
		itr != intron_doc.end(); 
		++itr)
	{
		for (size_t j = 0; j < hits.size(); ++j)
		{
			//fprintf(stderr, "considering read [%d-%d] with min doc = %lf contained in intron with doc = %lf\n", hits[j].left(), hits[j].right(), doc, idoc);
			const vector<AugmentedCuffOp>& ops = hits[j].augmented_ops();
			
			for (size_t i = 0; i < ops.size(); ++i)
			{
				if (ops[i].opcode == CUFF_INTRON)
				{
					map<pair<int, int>, int>::const_iterator itr;
					itr = intron_doc.find(make_pair(ops[i].g_left(), ops[i].g_right()));
					
					double doc = itr->second;
					if (doc < bundle_avg_thresh)
					{
						toss[j] = true;
#if ASM_VERBOSE
						fprintf(stderr, "\t Filtering intron %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
#endif	
						continue; 
					}
					
					if (!filter_on_intron_overlap)
						continue;
					
					for (map<pair<int,int>, int>::const_iterator itr2 = intron_doc.begin();
						 itr2 != intron_doc.end();
						 ++itr2)
					{	
						if (itr == itr2 ||
							!overlap_in_genome(itr->first.first,
											   itr->first.second,
											   itr2->first.first,
											   itr2->first.second))
							continue;
						
						double thresh = itr2->second * fraction;
						if (doc < thresh)
						{
#if ASM_VERBOSE
							fprintf(stderr, "\t Filtering intron (due to overlap) %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
#endif	
							toss[j] = true;
						}
					}
				}
			}
		}
	}
	
	for (size_t j = 0; j < hits.size(); ++j)
	{	
		if (!toss[j])
		{
			filtered_hits.push_back(hits[j]);
//#if ASM_VERBOSE
//			if (hits[j].has_intron())
//			{
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif
		}
		else
		{
#if ASM_VERBOSE
			if (hits[j].has_intron())
			{
				
				fprintf(stderr, "\tFiltering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
#endif	
		}
	}
	
#if ASM_VERBOSE
	
	fprintf(stderr, "\tIntron filtering pass finished: excluded %d fragments\n", (int)hits.size() - (int)filtered_hits.size());
#endif	
	hits = filtered_hits;
}

double background_rate(const vector<int> depth_of_coverage,
                       int left, 
                       int right)
{
    vector<int> tmp;
    
    size_t r_bound = (size_t)min(right, (int) depth_of_coverage.size());
    size_t l_bound = (size_t)max(left, 0);
    
    tmp.insert(tmp.end(), 
               depth_of_coverage.begin() + l_bound, 
               depth_of_coverage.begin() + r_bound);
    
    if (tmp.empty())
        return 0;
    
    vector<int>::iterator new_end =  remove(tmp.begin(), tmp.end(), 0);
    tmp.erase(new_end, tmp.end());
    sort(tmp.begin(), tmp.end());
    
    size_t median = (size_t)floor(tmp.size() / 2);
    double median_doc = tmp[median];
    return median_doc;
}

void pre_mrna_filter(int bundle_length,
					 int bundle_left,
					 vector<Scaffold>& hits)
{
	vector<int> depth_of_coverage(bundle_length,0);
	vector<double> scaff_doc;
	map<pair<int,int>, int> intron_doc;
	vector<Scaffold> filtered_hits;
	vector<bool> toss(hits.size(), false);
	
	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
    double bundle_avg_doc = compute_doc(bundle_left, 
										hits, 
										depth_of_coverage, 
										intron_doc,
										true);
    
	compute_doc(bundle_left, 
                hits, 
                depth_of_coverage, 
                intron_doc,
                false);
    
	record_doc_for_scaffolds(bundle_left, 
							 hits, 
							 depth_of_coverage, 
							 intron_doc,
							 scaff_doc);
    
    
	for(map<pair<int, int>, int >::const_iterator itr = intron_doc.begin();
		itr != intron_doc.end(); 
		++itr)
	{
		int i_left = itr->first.first;
		int i_right = itr->first.second;
		
        double intron_background = background_rate(depth_of_coverage, 
                                                   i_left - bundle_left, 
                                                   i_right - bundle_left);
        
        double cumul_cov = 0;
        for (size_t i = 0; i < i_right - i_left; ++i)
        {
            size_t pos = (i_left - bundle_left) + i;
            cumul_cov += depth_of_coverage[pos];
        }
        cumul_cov /= i_right - i_left;
#if ASM_VERBOSE
        fprintf(stderr, "retained intron %d-%d background: %lf\n", i_left, i_right, intron_background);
#endif
        if (cumul_cov / bundle_avg_doc >= pre_mrna_fraction)
        {
            //fprintf(stderr, "\tskipping\n");
            
            continue;
        }
        
        for (size_t j = 0; j < hits.size(); ++j)
        {
            //if (hits[j].has_intron())
            //    continue;
            double thresh = (1.0/pre_mrna_fraction) * intron_background;
            
            int len = 0;
            double doc = 0.0;
            size_t curr_op = 0;
            const vector<AugmentedCuffOp>& ops = hits[j].augmented_ops();
            while (curr_op != ops.size())
            {
                const AugmentedCuffOp&  op = ops[curr_op];
                
                if (op.opcode == CUFF_MATCH)
                {
                    int op_len = 0;
                    double op_doc = 0.0;
                    int left_off = op.g_left();
                    if (left_off + op.genomic_length > i_left && left_off < i_right)
                    {
                        if (left_off > i_left)
                        {
                            if (left_off + op.genomic_length <= i_right + 1)
                            {
                                op_len += op.genomic_length;
                                int L = left_off - bundle_left;
                                int R = L + op.genomic_length;
                                op_doc += accumulate(depth_of_coverage.begin() + L, depth_of_coverage.begin() + R, 0); 
                            }
                            else
                            {
                                op_len += i_right - left_off;
                                int L = left_off - bundle_left;
                                int R = L + (i_right - left_off);
                                op_doc += accumulate(depth_of_coverage.begin() + L, depth_of_coverage.begin() + R, 0);
                            }
                        }
                        else
                        {
                            if (left_off + op.genomic_length <= i_right + 1)
                            {
                                op_len += (left_off + op.genomic_length - i_left);
                                int L = left_off - bundle_left;
                                int R = L + (left_off + op.genomic_length - i_left);
                                op_doc += accumulate(depth_of_coverage.begin() + L, depth_of_coverage.begin() + R, 0);
                            }
                            else
                            {
                                op_len = i_right - i_left;
                                int L = left_off - bundle_left;
                                int R = L + (i_right - i_left);
                                op_doc = accumulate(depth_of_coverage.begin() + L, depth_of_coverage.begin() + R, 0);
                            }
                        }
                    }
                    
                    len += op_len;
                    doc += op_doc;
                }
                
                if (op.g_left() >= i_right)
                    break;
                ++curr_op;
            }
            
            if (len)
            {
                double hit_doc_in_region = doc / len;
                if (hit_doc_in_region < thresh)
                {
                    toss[j] = true;
//                    if (hits[j].has_intron())
//                    {
//                        fprintf(stderr, "\t$$$ Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//                    }
                }
            }
		}
	}
    
	for (size_t j = 0; j < hits.size(); ++j)
	{	
		if (!toss[j])
		{
			filtered_hits.push_back(hits[j]);
			//#if ASM_VERBOSE
			//			if (hits[j].has_intron())
			//			{
			//				
			//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			//			}
			//#endif	
		}
		else
		{
#if ASM_VERBOSE
			if (hits[j].has_intron())
			{
				
                fprintf(stderr, "\t@@@ Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
#endif	
		}
	}
    
#if ASM_VERBOSE
	fprintf(stderr, "\tPre-mRNA filter pass complete, excluded %lu fragments\n", hits.size() - filtered_hits.size());
#endif
	
	hits = filtered_hits;
}

//void pre_mrna_filter(int bundle_length,
//					 int bundle_left,
//					 vector<Scaffold>& hits)
//{
//	vector<int> depth_of_coverage(bundle_length,0);
//	vector<double> scaff_doc;
//	map<pair<int,int>, int> intron_doc;
//	vector<Scaffold> filtered_hits;
//	vector<bool> toss(hits.size(), false);
//	
//	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
//	double bundle_avg_doc = compute_doc(bundle_left, 
//										hits, 
//										depth_of_coverage, 
//										intron_doc,
//										true);
//	
//	// recompute the real DoCs
//	compute_doc(bundle_left, 
//				hits, 
//				depth_of_coverage, 
//				intron_doc,
//				false);
//	
//	record_doc_for_scaffolds(bundle_left, 
//							 hits, 
//							 depth_of_coverage, 
//							 intron_doc,
//							 scaff_doc);
//	
////	vector<int>::iterator new_end = remove(depth_of_coverage.begin(), depth_of_coverage.end(), 0);
////	depth_of_coverage.erase(new_end, depth_of_coverage.end());
////	sort(depth_of_coverage.begin(), depth_of_coverage.end());
////	
////	size_t median = floor(depth_of_coverage.size() / 2);
//	
//	
//	Scaffold smashed_gene;
//	vector<Scaffold> spliced_hits;
//    for (size_t i = 0; i < hits.size(); ++i)
//    {
//        const vector<const MateHit*>& m_hits = hits[i].mate_hits();
//        for (size_t j = 0; j < m_hits.size(); ++j)
//        {
//            if (m_hits[j]->left_alignment() && !m_hits[j]->left_alignment()->contiguous())
//            {
//                spliced_hits.push_back(Scaffold(MateHit(m_hits[j]->ref_id(),
//                                                        m_hits[j]->left_alignment(),
//                                                        shared_ptr<const ReadHit>(),
//                                                        0,
//                                                        0)));
//
//            }
//            if (m_hits[j]->right_alignment() && !m_hits[j]->right_alignment()->contiguous())
//            {
//                spliced_hits.push_back(Scaffold(MateHit(m_hits[j]->ref_id(),
//                                                        m_hits[j]->right_alignment(),
//                                                        shared_ptr<const ReadHit>(),
//                                                        0,
//                                                        0)));
//                
//            }
//        }
//    }
//    
//	Scaffold::merge(spliced_hits, smashed_gene, false);
//	vector<bool> constitutive_introns(intron_doc.size(), true);
//	
//    vector<pair<int, int> > gaps = smashed_gene.gaps();
//    
//	//size_t intron_idx = 0;
//    
//	for(map<pair<int, int>, int >::const_iterator itr = intron_doc.begin();
//		itr != intron_doc.end(); 
//		++itr)
//	{
//		int i_left = itr->first.first;
//		int i_right = itr->first.second;
//		
//        double cumul_cov = 0;
//        for (size_t i = 0; i < i_right - i_left; ++i)
//        {
//            size_t pos = (i_left - bundle_left) + i;
//            cumul_cov += depth_of_coverage[pos];
//        }
//        cumul_cov /= i_right - i_left;
//#if ASM_VERBOSE
//        fprintf(stderr, "retained intron %d-%d depth of cov: %lf\n", i_left, i_right, cumul_cov);
//#endif
//        if (cumul_cov / bundle_avg_doc >= pre_mrna_fraction)
//        {
//            continue;
//        }
//        
//        for (size_t j = 0; j < hits.size(); ++j)
//        {
//            //if (hits[j].has_intron())
//            //    continue;
//            if (hits[j].match_length(i_left, i_right))
//            {
//                toss[j] = true;
//            }
//		}
//	}
//    
//	for (size_t j = 0; j < hits.size(); ++j)
//	{	
//		if (!toss[j])
//		{
//			filtered_hits.push_back(hits[j]);
//			//#if ASM_VERBOSE
//			//			if (hits[j].has_intron())
//			//			{
//			//				
//			//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			//			}
//			//#endif	
//		}
//		else
//		{
//#if ASM_VERBOSE
//			//			if (hits[j].has_intron())
//			//			{
//			//				
//			//				fprintf(stderr, "\tFiltering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			//			}
//#endif	
//		}
//	}
//	
//	//#if ASM_VERBOSE
//	//	fprintf(stderr, "\tInitial filter pass complete\n");
//	//#endif
//	
//	hits = filtered_hits;
//}

void filter_hits(int bundle_length,
				 int bundle_left,
				 vector<Scaffold>& hits)
{
	
	pre_mrna_filter(bundle_length, bundle_left, hits);
	
	vector<int> depth_of_coverage(bundle_length,0);
	vector<double> scaff_doc;
	map<pair<int,int>, int> intron_doc;
	vector<Scaffold> filtered_hits;
	vector<bool> toss(hits.size(), false);
	
	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
	double bundle_avg_doc = compute_doc(bundle_left, 
										hits, 
										depth_of_coverage, 
										intron_doc,
										true);
	
	// recompute the real DoCs
	compute_doc(bundle_left, 
				hits, 
				depth_of_coverage, 
				intron_doc,
				false);
	
	
	record_min_doc_for_scaffolds(bundle_left, 
								 hits, 
								 depth_of_coverage, 
								 intron_doc,
								 scaff_doc);
	
	//double bundle_avg_thresh = min_isoform_fraction * bundle_avg_doc;
	
	if (!intron_doc.empty())
	{
		double intron_avg_doc = major_isoform_intron_doc(intron_doc);
		double intron_multiplier = intron_avg_doc / bundle_avg_doc;
		
		// we don't want this to be more than 1.0 ...
		intron_multiplier = min(intron_avg_doc, 1.0);
		//bundle_avg_thresh = min_isoform_fraction * bundle_avg_doc;
		
		set<pair<int, int> > tossed_introns;
		for(map<pair<int, int>, int>::const_iterator itr = intron_doc.begin();
			itr != intron_doc.end(); 
			++itr)
		{
			for (size_t j = 0; j < hits.size(); ++j)
			{
				int i_left = itr->first.first;
				int i_right = itr->first.second;
				int j_match_len = hits[j].match_length(i_left, i_right); 
				if (j_match_len > 0)
				{
					double idoc = itr->second;
					double doc = scaff_doc[j];
					
					if (!hits[j].has_intron() && 
						doc < pre_mrna_fraction * (idoc * intron_multiplier))
                    {
						toss[j] = true;
                    }
					
					const vector<AugmentedCuffOp>& ops = hits[j].augmented_ops();
					
					unsigned int num_mismatches = 0;
					assert (hits[j].mate_hits().size() == 1);
					const MateHit& hit = **(hits[j].mate_hits().begin());
					num_mismatches = hit.edit_dist();
					
					double percent_mismatches = num_mismatches / (double)hits[j].length();
					
					bool intron_pokin_read = false;
					
					const AugmentedCuffOp& first = ops.front();
					// intron                    =========================
					// hit                     ******************
					if (first.g_left() < i_left && first.g_right() > i_left  && first.g_right() < i_right)
					{
						intron_pokin_read = true;
					}
					
					// intron          =========================
					// hit                      ******************
					if (first.g_left() < i_right && first.g_right() > i_right  && first.g_left() > i_left)
					{
						intron_pokin_read = true;
					}
					
					const AugmentedCuffOp& last = ops.back();
					// intron          =========================
					// hit                     ******************
					if (last.g_left() < i_left && last.g_right() > i_left  && last.g_right() < i_right)
					{
						intron_pokin_read = true;
					}
					
					// intron          =========================
					// hit                            ******************
					if (last.g_left() < i_right && last.g_right() > i_right  && last.g_left() > i_left)
					{
						intron_pokin_read = true;
					}
					
					if (intron_pokin_read)
					{
						double fraction;
//						if (!hits[j].has_intron())
//						{ 
//							fraction = (3 * pre_mrna_fraction) + percent_mismatches;
//						}
//						else
						{
							fraction = pre_mrna_fraction + percent_mismatches;
						}
						double thresh = fraction * (intron_avg_doc * intron_multiplier);
						if (doc < thresh)
						{
							toss[j] = true; 
//							if (hits[j].has_intron())
//							{
//								fprintf(stderr, "\t^^^Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//							}
						}
					}
				}
			}
		}
	}
	
	for (size_t j = 0; j < hits.size(); ++j)
	{	
		if (!toss[j])
		{
			filtered_hits.push_back(hits[j]);
//#if ASM_VERBOSE
//			if (hits[j].has_intron())
//			{
//				
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif	
		}
		else
		{
#if ASM_VERBOSE
			if (hits[j].has_intron())
			{
				
				fprintf(stderr, "\t!!!Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
#endif	
		}
	}
	
//#if ASM_VERBOSE
//	fprintf(stderr, "\tInitial filter pass complete\n");
//#endif
	
	hits = filtered_hits;
	
	scaff_doc.clear();
	filtered_hits.clear();
	
	toss = vector<bool>(hits.size(), false);
	
	map<pair<int, int>, int> dummy;
	bundle_avg_doc = compute_doc(bundle_left, 
								 hits, 
								 depth_of_coverage, 
								 dummy,
								 false);
	
//#if ASM_VERBOSE
//	fprintf(stderr, "\tUpdated avg bundle doc = %lf\n", bundle_avg_doc);
//#endif
	
	record_doc_for_scaffolds(bundle_left, 
							 hits, 
							 depth_of_coverage, 
							 intron_doc,
							 scaff_doc);
	
	
	
//#if ASM_VERBOSE
//    double bundle_thresh = pre_mrna_fraction * bundle_avg_doc;
//	fprintf(stderr, "\tthreshold is = %lf\n", bundle_thresh);
//#endif
	
	if (!intron_doc.empty())
	{
//		filter_introns(bundle_length, 
//					   bundle_left, 
//					   hits, 
//					   min_isoform_fraction, 
//					   true,
//					   true);
		if (bundle_avg_doc > 3000)
		{
			filter_introns(bundle_length, 
						   bundle_left, 
						   hits, 
						   min_isoform_fraction, 
						   true,
						   false);
		}
	}
	
	for (size_t j = 0; j < hits.size(); ++j)
	{	
		if (!toss[j])
		{
			filtered_hits.push_back(hits[j]);
//#if ASM_VERBOSE
//			if (hits[j].has_intron())
//			{
//				
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif
		}
		else
		{
#if ASM_VERBOSE
			if (hits[j].has_intron())
			{
				
				fprintf(stderr, "\t***Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
#endif	
		}
	}
	
	//fprintf(stderr, "\tTossed %d hits as noise\n", (int)hits.size() - filtered_hits.size());
	
	hits = filtered_hits;
}


void filter_junk_isoforms(vector<shared_ptr<Abundance> >& transcripts,
						  vector<double>& abundances)
{
	//	vector<double>::iterator max_ab = std::max_element(abundances.begin(),
	//													   abundances.end());
	double max_fwd_ab = -1.0;
	double max_rev_ab = -1.0;
	
	for (size_t t = 0; t < transcripts.size(); ++t)
	{
		shared_ptr<Scaffold> scaff = transcripts[t]->transfrag();
		if (scaff->strand() == CUFF_FWD || scaff->strand() == CUFF_STRAND_UNKNOWN)
		{
			if (abundances[t] > max_fwd_ab)
				max_fwd_ab = abundances[t];
		}
		if (scaff->strand() == CUFF_REV || scaff->strand() == CUFF_STRAND_UNKNOWN)
		{			
			if (abundances[t] > max_rev_ab)
				max_rev_ab = abundances[t];
		}
	}
	
	// Try to categorize the crap transcripts for suppression
	vector<bool> pre_mrna_junk(transcripts.size(), false); //intra-intron, much lower abundance than container
	vector<bool> chaff(transcripts.size(), false); // only a single MateHit, impossible to reliably quantitate
	vector<bool> repeats(transcripts.size(), false); // too many low-quality hits
	vector<bool> too_rare(transcripts.size(), false); // too rare to be reliably quantitated, could be error
	
	vector<bool> illegal_microexon(transcripts.size(), false); // initial or terminal exons are too short
	
	//cerr << "Chucked : ";
	for (size_t t = 0; t < transcripts.size(); ++t)
	{
		shared_ptr<Scaffold> scaff = transcripts[t]->transfrag();

		if (allow_junk_filtering)
		{
			const vector<const MateHit*> hits = scaff->mate_hits();
			
			const vector<AugmentedCuffOp>& ops = scaff->augmented_ops();
			
			if (ops.size() == 1 && ops[0].opcode == CUFF_MATCH)
			{
				for (size_t j = 0; j < transcripts.size(); ++j)
				{
					const vector<AugmentedCuffOp>& j_ops = scaff->augmented_ops();
					for (size_t L = 0; L < j_ops.size(); L++)
					{
						if (AugmentedCuffOp::overlap_in_genome(ops[0], j_ops[L]) &&
							j_ops[L].opcode == CUFF_INTRON)
						{
							pre_mrna_junk[t] = true;
						}
					}
				}
			}
			
			double low_qual_hits = 0.0;
			static const double low_qual_err_prob = high_phred_err_prob; // hits with error_prob() above this are low quality;
			static const double low_qual_thresh = 0.75; // hits with more than this fraction of low qual hits are repeats
			for (vector<const MateHit*>::const_iterator itr = hits.begin();
				 itr != hits.end();
				 ++itr)
			{
				double e = (*itr)->error_prob();
				if (e >= low_qual_err_prob)
					low_qual_hits += 1.0;
			}
			double low_qual_frac = low_qual_hits / (double)hits.size();
			if (low_qual_frac > low_qual_thresh)
				repeats[t] = true;
		}
		
		if (scaff->strand() == CUFF_FWD &&
			(abundances[t] / max_fwd_ab) < min_isoform_fraction)
			too_rare[t] = true;
		if ((scaff->strand() == CUFF_REV ||  scaff->strand() == CUFF_STRAND_UNKNOWN) &&
			(abundances[t] / max_rev_ab) < min_isoform_fraction)
			too_rare[t] = true;
		
	}
	
	vector<shared_ptr<Abundance> > non_junk_transcripts;
	vector<double> non_junk_abundances;
	for (size_t t = 0; t < transcripts.size(); ++t)
	{
		if (!repeats[t] && !pre_mrna_junk[t] && !too_rare[t])
		{
			non_junk_transcripts.push_back(transcripts[t]);
			non_junk_abundances.push_back(abundances[t]);
		}
        else
        {
#if ASM_VERBOSE
            fprintf(stderr, "Filtering isoform %d-%d\n", transcripts[t]->transfrag()->left(), transcripts[t]->transfrag()->right());
#endif
        }
	}
	
	transcripts = non_junk_transcripts;
	abundances = non_junk_abundances;
}

// Designed to strip out remaining pre-mrna genes, assembled repeats, and 
// fragments from isoforms too short to be reliably quantitated.
void filter_junk_genes(vector<Gene>& genes)
{
	vector<Gene> good_genes;
	vector<Isoform> all_isoforms;
	for (size_t i = 0; i < genes.size(); ++i)
	{
		all_isoforms.insert(all_isoforms.end(), 
							genes[i].isoforms().begin(), 
							genes[i].isoforms().end());
	}
	
	for (size_t i = 0; i < genes.size(); ++i)
	{
		const Gene& g = genes[i];
		bool good_gene = true;
		for (size_t j = 0; j < all_isoforms.size(); ++j)
		{
			vector<pair<int, int> > introns = all_isoforms[j].scaffold().gaps();
			for (size_t k = 0; k < introns.size(); ++k)
			{
				if (g.left() > introns[k].first && g.right() < introns[k].second &&
					g.FPKM() / all_isoforms[j].FPKM() < pre_mrna_fraction)
				{
					good_gene = false;
				}
			}
		}
		if (good_gene)
        {
			good_genes.push_back(g);
        }
        else
        {
#if ASM_VERBOSE
            fprintf(stderr, "Filtering transfrags from gene %d-%d\n", g.left(), g.right());
#endif
        }
	}
	
	genes = good_genes;
	
}
