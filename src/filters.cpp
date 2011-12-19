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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

void filter_introns(int bundle_length,
					int bundle_left,
					vector<Scaffold>& hits, 
					double fraction,
					bool filter_on_intron_overlap,
					bool filter_with_intron_doc)
{
	vector<float> depth_of_coverage(bundle_length,0);
	vector<double> scaff_doc;
	map<pair<int,int>, float> intron_doc;
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
		verbose_msg("\tFiltering bundle introns, avg (intron) doc = %lf, thresh = %f\n", bundle_avg_doc, bundle_avg_thresh);
	}
	else
	{
		verbose_msg("\tFiltering bundle introns, avg bundle doc = %lf, thresh = %f\n", bundle_avg_doc, bundle_avg_thresh);
	}
	
	for(map<pair<int, int>, float>::const_iterator itr = intron_doc.begin();
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
					map<pair<int, int>, float>::const_iterator itr;
					itr = intron_doc.find(make_pair(ops[i].g_left(), ops[i].g_right()));
					
					double doc = itr->second;
					if (doc < bundle_avg_thresh)
					{
						toss[j] = true;
						verbose_msg("\t Filtering intron %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
						continue; 
					}
					
					if (!filter_on_intron_overlap)
						continue;
					
					for (map<pair<int,int>, float>::const_iterator itr2 = intron_doc.begin();
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
							verbose_msg("\t Filtering intron (due to overlap) %d - %d: %f thresh %f\n", itr->first.first, itr->first.second, doc, bundle_avg_thresh);
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
//#if verbose_msg
//			if (hits[j].has_intron())
//			{
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif
		}
		else
		{
			if (hits[j].has_intron())
			{
				
				verbose_msg("\tFiltering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
		}
	}
	
	
	verbose_msg("\tIntron filtering pass finished: excluded %d fragments\n", (int)hits.size() - (int)filtered_hits.size());
	hits = filtered_hits;
}

double background_rate(const vector<float> depth_of_coverage,
                       int left, 
                       int right)
{
    vector<float> tmp;
    
    size_t r_bound = (size_t)min(right, (int) depth_of_coverage.size());
    size_t l_bound = (size_t)max(left, 0);
    
    tmp.insert(tmp.end(), 
               depth_of_coverage.begin() + l_bound, 
               depth_of_coverage.begin() + r_bound);
    
    if (tmp.empty())
        return 0;
    
    vector<float>::iterator new_end =  remove(tmp.begin(), tmp.end(), 0);
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
	vector<float> depth_of_coverage(bundle_length,0);
	vector<double> scaff_doc;
	map<pair<int,int>, float> intron_doc;
	vector<Scaffold> filtered_hits;
	vector<bool> toss(hits.size(), false);
	vector<float> through_introns; //for each location, how many introns pass through

	vector<int> scaff_intron_status;
	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
  double bundle_avg_doc = compute_doc(bundle_left,
										hits, 
										depth_of_coverage, 
										intron_doc,
										true,
										&through_introns,
										&scaff_intron_status);
   verbose_msg("Pre-mRNA flt: bundle average doc = %lf\n", bundle_avg_doc);
  /*
   //2nd call not needed, the vectors won't change, only the return value
	compute_doc(bundle_left, 
                hits, 
                depth_of_coverage, 
                intron_doc,
                false);
  */
	record_doc_for_scaffolds(bundle_left, 
							 hits, 
							 depth_of_coverage, 
							 intron_doc,
							 scaff_doc);

	for(map<pair<int, int>, float >::const_iterator itr = intron_doc.begin();
		itr != intron_doc.end(); 
		++itr)
	{
		int i_left = itr->first.first;
		int i_right = itr->first.second;
		int i_doc   = itr->second;
    double intron_background = background_rate(depth_of_coverage,
                                               i_left - bundle_left,
                                               i_right - bundle_left);
    double cumul_cov = 0;
    for (int i = 0; i < i_right - i_left; ++i)
    {
        size_t pos = (i_left - bundle_left) + i;
        cumul_cov += depth_of_coverage[pos];
    }
    cumul_cov /= i_right - i_left;
    verbose_msg("Pre-mRNA flt: intron %d-%d : background: %lf, inner coverage: %lf, junction coverage: %f\n",
                 i_left, i_right, intron_background, cumul_cov, i_doc);
    if (cumul_cov / bundle_avg_doc >= pre_mrna_fraction)
    {
        //fprintf(stderr, "\tskipping\n");
        continue;
    }
        
		////double thresh = (1.0/pre_mrna_fraction) * intron_background;
		double thresh = pre_mrna_fraction * intron_background;
    float min_flt_fraction = min(pre_mrna_fraction, min_isoform_fraction);
    //double thresh = min_flt_fraction * i_doc;

    for (size_t j = 0; j < hits.size(); ++j)
        {
           if (hits[j].left()>i_right) break;
            if (hits[j].is_ref())
                continue;
			      if (toss[j])
			          continue;
			      //find maximum intron support in the hit region

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
                    if (hits[j].has_intron())
                    {
                    //    fprintf(stderr, "\t$$$ Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());

                       verbose_msg("\t@@@ Filtering intron scaff [%d-%d] (scaff_doc=%lf, doc_in_region=%lf)\n",
                            hits[j].left(), hits[j].right(), scaff_doc[j], hit_doc_in_region);
                    }
                }
            }
		} //for each scaffold
	} //for each intron
	for (size_t j = 0; j < hits.size(); ++j)
	{	
		if (!toss[j])
		{
			filtered_hits.push_back(hits[j]);
		}
		/*else
		{
			if (hits[j].has_intron())
			{
				
                verbose_msg( "\t@@@ Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
		}
		*/
	}
    
	if (cuff_verbose && hits.size()>filtered_hits.size())
	  verbose_msg("\tPre-mRNA flt tossed %lu fragments\n", hits.size() - filtered_hits.size());
	
	hits = filtered_hits;
}

void filter_hits(int bundle_length,
				 int bundle_left,
				 vector<Scaffold>& hits)
{
	
	pre_mrna_filter(bundle_length, bundle_left, hits);
	
	vector<float> depth_of_coverage(bundle_length+1,0);
	vector<double> scaff_doc;
	map<pair<int,int>, float> intron_doc;
	vector<Scaffold> filtered_hits;
	vector<bool> toss(hits.size(), false);
	
	// Make sure the avg only uses stuff we're sure isn't pre-mrna fragments
	double bundle_avg_doc = compute_doc(bundle_left, 
										hits, 
										depth_of_coverage, 
										intron_doc,
										true);
	
	// recompute the real DoCs
	/* not needed, vectors are not changed
	compute_doc(bundle_left, 
				hits, 
				depth_of_coverage, 
				intron_doc,
				false);
	*/
	
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
		for(map<pair<int, int>, float>::const_iterator itr = intron_doc.begin();
			itr != intron_doc.end(); 
			++itr)
		{
			for (size_t j = 0; j < hits.size(); ++j)
			{
				if (hits[j].is_ref())
                {
					continue;
                }
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
//#if verbose_msg
//			if (hits[j].has_intron())
//			{
//				
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif	
		}
		else
		{
			if (hits[j].has_intron())
			{
				
				verbose_msg("\t!!!Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
		}
	}
	
//#if verbose_msg
//	fprintf(stderr, "\tInitial filter pass complete\n");
//#endif
	
	hits = filtered_hits;
	
	scaff_doc.clear();
	filtered_hits.clear();
	
	toss = vector<bool>(hits.size(), false);
	
	map<pair<int, int>, float> dummy;
	bundle_avg_doc = compute_doc(bundle_left, 
								 hits, 
								 depth_of_coverage, 
								 dummy,
								 false);
	
//#if verbose_msg
//	fprintf(stderr, "\tUpdated avg bundle doc = %lf\n", bundle_avg_doc);
//#endif
	
	record_doc_for_scaffolds(bundle_left, 
							 hits, 
							 depth_of_coverage, 
							 intron_doc,
							 scaff_doc);
	
	
	
//#if verbose_msg
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
//#if verbose_msg
//			if (hits[j].has_intron())
//			{
//				
//				fprintf(stderr, "KEEPING intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
//			}
//#endif
		}
		else
		{
			if (hits[j].has_intron())
			{
				verbose_msg("\t***Filtering intron scaff [%d-%d]\n", hits[j].left(), hits[j].right());
			}
		}
	}
	
	//fprintf(stderr, "\tTossed %d hits as noise\n", (int)hits.size() - filtered_hits.size());
	
	hits = filtered_hits;
}


void filter_junk_isoforms(vector<shared_ptr<Abundance> >& transcripts,
						  vector<double>& abundances,
                          const vector<shared_ptr<Abundance> >& mapped_transcripts,
                          double locus_mass)
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
	
	//cerr << "Chucked : ";
	for (size_t t = 0; t < transcripts.size(); ++t)
	{
		shared_ptr<Scaffold> scaff = transcripts[t]->transfrag();
        
		if (!(scaff->is_ref()) && allow_junk_filtering)
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
			
            if (library_type != "transfrags")
            {
                double low_qual_hits = 0.0;
                static const double low_qual_err_prob = high_phred_err_prob; // hits with error_prob() above this are low quality;
                static const double low_qual_thresh = 0.75; // hits with more than this fraction of low qual hits are repeats
                for (vector<const MateHit*>::const_iterator itr = hits.begin();
                     itr != hits.end();
                     ++itr)
                {
                    double e = 1-(*itr)->mass();
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

            const vector<double>* cond_probs = (mapped_transcripts[t]->cond_probs());
            if (cond_probs)
            {
                assert (library_type != "transfrags");
                double supporting_hits = abundances[t] * locus_mass;
                if (supporting_hits < min_frags_per_transfrag)
                    chaff[t] = true;
            }
		}
//        else // we should still filter things that are zero to improve robustness of MAP estimation
//        {
//            if (abundances[t] == 0.0)
//                too_rare[t] = true;
//        }
	}
	
	vector<shared_ptr<Abundance> > non_junk_transcripts;
	vector<double> non_junk_abundances;
	for (size_t t = 0; t < transcripts.size(); ++t)
	{
		if (!repeats[t] && !pre_mrna_junk[t] && !too_rare[t] && !chaff[t])
		{
			non_junk_transcripts.push_back(transcripts[t]);
			non_junk_abundances.push_back(abundances[t]);
		}
        else
        {
            verbose_msg( "Filtering isoform %d-%d\n", transcripts[t]->transfrag()->left(), transcripts[t]->transfrag()->right());
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
		
		if(g.has_ref_trans())
		{
			good_genes.push_back(g);
			continue;
		}
		
		bool good_gene = true;
		for (size_t j = 0; j < all_isoforms.size(); ++j)
		{
			vector<pair<int, int> > introns = all_isoforms[j].scaffold().gaps();
            
            //assert (!allow_junk_filtering || all_isoforms[j].scaffold().mate_hits().size() >= min_frags_per_transfrag);
			for (size_t k = 0; k < introns.size(); ++k)
			{
				if (g.left() > introns[k].first && g.right() < introns[k].second &&
					g.FPKM() / all_isoforms[j].FPKM() < pre_mrna_fraction)
				{
					good_gene = false;
				}
			}
		}
        if (allow_junk_filtering)
        {
            if (g.FPKM() == 0)
            {
                good_gene = false;
            }
        }
		if (good_gene)
        {
			good_genes.push_back(g);
        }
        else
        {
            verbose_msg("Filtering transfrags from gene %d-%d\n", g.left(), g.right());
        }
	}
	
	genes = good_genes;
	
}

void clip_by_3_prime_dropoff(vector<Scaffold>& scaffolds)
{
    vector<pair<double, Scaffold*> > three_prime_ends;
    
    if (library_type != "transfrags")
    {
        foreach (Scaffold& scaff, scaffolds)
        {
            if (!(scaff.strand() == CUFF_FWD || scaff.strand() == CUFF_REV))
                continue;

            int scaff_len = scaff.length();
            vector<double> coverage(scaff_len, 0.0);
            
            double total = 0;
            foreach(const MateHit* hit, scaff.mate_hits())
            {
                int start, end, frag_len;
                if (!scaff.map_frag(*hit, start, end, frag_len)) continue;
                
                if (scaff.strand() == CUFF_REV)
                {
                    start = scaff_len - 1 - start;
                    end = scaff_len - 1 - end;
                    swap(start, end);
                }
                
                for(int i = start; i <= end; ++i)
                {
                    coverage[i] += hit->mass();
                    total += hit->mass();
                }
            }
            double avg_cov = total/scaff_len;
    //        if (avg_cov < trim_3_avgcov_thresh)
    //            continue;
            
            const AugmentedCuffOp* exon_3 = NULL;
            int mult;
            int offset;
            
            if (scaff.strand() == CUFF_REV)
            {
                mult = 1;
                offset = 0;
                exon_3 = &scaff.augmented_ops().front();
            }
            else if (scaff.strand() == CUFF_FWD)
            {
                mult = -1;
                offset = scaff_len - 1;
                exon_3 = &scaff.augmented_ops().back();
            }
            else
            {
                continue;
            }

            int to_remove;
            double min_cost = numeric_limits<double>::max();
            double mean_to_keep = 0.0;
            double mean_to_trim = 0.0;
            double tmp_mean_to_trim = 0.0;
            double tmp_mean_to_keep = 0.0;
            double tmp_mean_3prime = 0.0;
            for (int i = 0; i < exon_3->genomic_length; i++)
            {
                tmp_mean_3prime += coverage[offset + mult*i];
            }
            tmp_mean_3prime /= exon_3->genomic_length;
            
            double base_cost = 0.0;
            for (int i = 0; i < exon_3->genomic_length; i++)
            {
                double d = (coverage[offset + mult*i] - tmp_mean_3prime);
                d *= d;
                base_cost += d;
            }
            base_cost /= exon_3->genomic_length;
            
            size_t min_cost_x = -1;
            for (to_remove = 1; to_remove < exon_3->genomic_length - 1; to_remove++)
            {
                tmp_mean_to_trim = 0.0;
                tmp_mean_to_keep = 0.0;
                for (size_t i = 0; i < exon_3->genomic_length; i++)
                {
                    if (i <= to_remove)
                    {
                        tmp_mean_to_trim += coverage[offset + mult*i];
                    }
                    else 
                    {
                        tmp_mean_to_keep += coverage[offset + mult*i];
                    }
                }
                
                tmp_mean_to_trim /= to_remove;
                tmp_mean_to_keep /= (exon_3->genomic_length - to_remove);
                
                double tmp_mean_trim_cost = 0.0;
                double tmp_mean_keep_cost = 0.0;
                for (int i = 0; i < exon_3->genomic_length; i++)
                {
                    if (i <= to_remove)
                    {
                        double d = (coverage[offset + mult*i] - tmp_mean_to_trim);
                        d *= d;
                        tmp_mean_trim_cost += d;
                    }
                    else 
                    {
                        double d = (coverage[offset + mult*i] - tmp_mean_to_keep);
                        d *= d;
                        tmp_mean_keep_cost += d;
                    }
                }
                
                tmp_mean_trim_cost /= to_remove;
                tmp_mean_keep_cost /= (exon_3->genomic_length - to_remove);
                
                double new_cost = tmp_mean_trim_cost + tmp_mean_keep_cost;
                
                if (new_cost < min_cost && trim_3_dropoff_frac * tmp_mean_to_keep > tmp_mean_to_trim && new_cost < base_cost && to_remove > scaff_len * 0.05)
                {
                    min_cost = tmp_mean_trim_cost + tmp_mean_keep_cost;
                    min_cost_x = to_remove;
                    mean_to_keep = tmp_mean_to_keep;
                    mean_to_trim = tmp_mean_to_trim;
                }
            }
            
            // If trimming reduces the overall mean squared error of the coverage
            // do it
            if (avg_cov >= trim_3_avgcov_thresh && min_cost_x < exon_3->genomic_length)
            {
                scaff.trim_3(min_cost_x);
            }
            
            // store the mean squared error for this exon
            tmp_mean_3prime = 0.0;
            for (int i = 0; i < exon_3->genomic_length; i++)
            {
                tmp_mean_3prime += coverage[offset + mult*i];
            }
            tmp_mean_3prime /= exon_3->genomic_length;
            
            base_cost = 0.0;
            for (int i = 0; i < exon_3->genomic_length; i++)
            {
                double d = (coverage[offset + mult*i] - tmp_mean_3prime);
                d *= d;
                base_cost += d;
            }
            base_cost /= exon_3->genomic_length;
            three_prime_ends.push_back(make_pair(base_cost, &scaff));
        }
    }
    else
    {
        foreach (Scaffold& scaff, scaffolds)
        {
            if (!(scaff.strand() == CUFF_FWD || scaff.strand() == CUFF_REV))
                continue;
            
            int scaff_len = scaff.length();
            vector<double> coverage(scaff_len, 0.0);
            
            double total = 0;
            foreach(const MateHit* hit, scaff.mate_hits())
            {
                int start, end, frag_len;
                if (!scaff.map_frag(*hit, start, end, frag_len)) continue;
                
                if (scaff.strand() == CUFF_REV)
                {
                    start = scaff_len - 1 - start;
                    end = scaff_len - 1 - end;
                    swap(start, end);
                }
                
                for(int i = start; i <= end; ++i)
                {
                    coverage[i] += hit->mass();
                    total += hit->mass();
                }
            }
            double avg_cov = total/scaff_len;
            //        if (avg_cov < trim_3_avgcov_thresh)
            //            continue;
            
            const AugmentedCuffOp* exon_3 = NULL;
            int mult;
            int offset;
            
            if (scaff.strand() == CUFF_REV)
            {
                mult = 1;
                offset = 0;
                exon_3 = &scaff.augmented_ops().front();
            }
            else if (scaff.strand() == CUFF_FWD)
            {
                mult = -1;
                offset = scaff_len - 1;
                exon_3 = &scaff.augmented_ops().back();
            }
            else
            {
                continue;
            }
            
            three_prime_ends.push_back(make_pair(scaff.fpkm(), &scaff));
        }
        
    }
    
    adjacency_list <vecS, vecS, undirectedS> G;
	
	for (size_t i = 0; i < three_prime_ends.size(); ++i)
	{
		add_vertex(G);
	}
	
	for (size_t i = 0; i < three_prime_ends.size(); ++i)
	{
		Scaffold* scaff_i = three_prime_ends[i].second;
		//assert (scaff_i);
		
        const AugmentedCuffOp* scaff_i_exon_3 = NULL;
        
        if (scaff_i->strand() == CUFF_REV)
        {
            scaff_i_exon_3 = &(scaff_i->augmented_ops().front());
        }
        else if (scaff_i->strand() == CUFF_FWD)
        {
            scaff_i_exon_3 = &(scaff_i->augmented_ops().back());
        }
        
		for (size_t j = i + 1; j < three_prime_ends.size(); ++j)
		{
			Scaffold* scaff_j = three_prime_ends[j].second;
            
            if (scaff_i->strand() != scaff_j->strand())
                continue;
            
            const AugmentedCuffOp* scaff_j_exon_3 = NULL;
            
            if (scaff_j->strand() == CUFF_REV)
            {
                scaff_j_exon_3 = &(scaff_j->augmented_ops().front());
            }
            else if (scaff_j->strand() == CUFF_FWD)
            {
                scaff_j_exon_3 = &(scaff_j->augmented_ops().back());
            }
			
			if (AugmentedCuffOp::overlap_in_genome(*scaff_j_exon_3, *scaff_i_exon_3) && 
                AugmentedCuffOp::compatible(*scaff_j_exon_3, *scaff_i_exon_3, 0))
				add_edge(i, j, G);
		}
	}
	
	std::vector<int> component(num_vertices(G));
	connected_components(G, &component[0]);
	
	vector<vector<bool> > clusters(three_prime_ends.size(), 
								   vector<bool>(three_prime_ends.size(), false));
	
	//vector<vector<size_t> > cluster_indices(three_prime_ends.size());
    
    vector<vector<pair<double, Scaffold*> > > grouped_scaffolds(three_prime_ends.size());
	for (size_t i = 0; i < three_prime_ends.size(); ++i)
	{
		clusters[component[i]][i] = true;
		grouped_scaffolds[component[i]].push_back(three_prime_ends[i]);
	}
    
    for (size_t i = 0; i < grouped_scaffolds.size(); ++i)
    {
        vector<pair<double, Scaffold*> >& group = grouped_scaffolds[i];
        sort(group.begin(), group.end());
        if (group.empty())
            continue;
        
        Scaffold* group_leader = NULL;
        int trim_point = -1;
        
        const AugmentedCuffOp* group_exon_3 = NULL;
        vector<pair<double, Scaffold*> >::iterator l_itr = group.begin();
        while (l_itr != group.end())
        {
            Scaffold* possible_leader = l_itr->second;
            bool ok_clip_leader = true;
            vector<pair<double, Scaffold*> >::iterator g_itr = group.begin();
            const AugmentedCuffOp* l_exon_3 = NULL;
            CuffStrand s = possible_leader->strand();

            if (s != CUFF_STRAND_UNKNOWN)
            {  
                if (s == CUFF_REV)
                    l_exon_3 = &(possible_leader->augmented_ops().front());
                else 
                    l_exon_3 = &(possible_leader->augmented_ops().back());
                for (; g_itr != group.end(); ++g_itr)
                {
                    const AugmentedCuffOp* g_exon_3 = NULL;
                    if (s == CUFF_REV)
                    {
                        //  bad:
                        //              leader  
                        //    follower
                        g_exon_3 = &(g_itr->second->augmented_ops().front());
                        if (g_exon_3->g_right() <= l_exon_3->g_left())
                            ok_clip_leader = false;
                        
                        // for meta-assembly libraries, don't ever allow clipping, just extension
                        // bad:
                        //             leader
                        //         follower
                        if (library_type == "transfrags" && 
                            g_exon_3->g_left() < l_exon_3->g_left())
                            ok_clip_leader = false;
                    }
                    else 
                    {
                        //  bad:
                        //          follower  
                        //  leader
                        g_exon_3 = &(g_itr->second->augmented_ops().back());
                        if (g_exon_3->g_left() >= l_exon_3->g_right())
                            ok_clip_leader = false;
                        
                        // for meta-assembly libraries, don't ever allow clipping, just extension
                        // bad:
                        //             leader
                        //                follower
                        if (library_type == "transfrags" && 
                            g_exon_3->g_right() > l_exon_3->g_right())
                            ok_clip_leader = false;
                    }
                }
            }
            else
            {
                ok_clip_leader = false;
            }
            
            if (ok_clip_leader)
            {
                if (s == CUFF_REV)
                {
                    if (trim_point == -1)
                        trim_point = l_exon_3->g_left();
                    else if (l_exon_3->g_left() < trim_point)
                        ok_clip_leader = false;
                }
                else 
                {
                    if (trim_point == -1)
                        trim_point = l_exon_3->g_right();
                    else if (l_exon_3->g_right() > trim_point)
                        ok_clip_leader = false;
                }
            }
            
            if (ok_clip_leader)
            {
                group_leader = possible_leader;
                group_exon_3 = l_exon_3;
                break;
            }
            ++l_itr;
        }
        
        if (!group_leader || !group_exon_3)
            continue;
        
        for (size_t j = 0; j < group.size(); ++j)
        {
            const AugmentedCuffOp* exon_3 = NULL;
            int end_diff = 0;
            if (group_leader->strand() == CUFF_REV)
            {
                exon_3 = &(group[j].second->augmented_ops().front());
                end_diff = group_exon_3->g_left() - exon_3->g_left();
            }
            else 
            {
                exon_3 = &(group[j].second->augmented_ops().back());                
                end_diff = exon_3->g_right() - group_exon_3->g_right();
            }
            
            if (end_diff > 0)
            {
                // leader
                //     follower
                group[j].second->trim_3(end_diff);
            }
            else if (end_diff < 0)
            {
                //        leader
                //   follower
                group[j].second->extend_3(-end_diff);
            }
        }
    }
    
    
	return;
	
}
