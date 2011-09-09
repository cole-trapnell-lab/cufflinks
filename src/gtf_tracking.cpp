/*
 *  gtf_tracking.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/5/09.
 *  Copyright 2009 Geo Pertea. All rights reserved.
 *
 */

#include "gtf_tracking.h"

bool gtf_tracking_verbose = false;
bool gtf_tracking_largeScale=false; //many input Cufflinks files processed at once by cuffcompare, discard exon attributes

int GXConsensus::count=0;

char* getGSeqName(int gseq_id) {
 return GffObj::names->gseqs.getName(gseq_id);
}

int cmpByPtr(const pointer p1, const pointer p2) {
  return (p1>p2) ? 1: ((p1==p2)? 0 : -1);
  }

bool betterRef(GffObj* a, GffObj* b) {
 if (a==NULL || b==NULL) return (a!=NULL);
 if (a->exons.Count()!=b->exons.Count()) return (a->exons.Count()>b->exons.Count());
 if (a->hasCDS() && !b->hasCDS())
        return true;
   else {
     if (b->hasCDS() && !a->hasCDS()) return false;
     return (a->covlen>b->covlen);
     }
 }

GffObj* is_RefDup(GffObj* m, GList<GffObj>& mrnas, int& dupidx) {
 //mrnas MUST be sorted by start coordinate
  int ovlen=0;
  dupidx=-1;
  if (mrnas.Count()==0) return NULL;
  int nidx=qsearch_mrnas(m->end, mrnas);
  if (nidx==0) return NULL;
  if (nidx==-1) nidx=mrnas.Count();//all can overlap
  for (int i=nidx-1;i>=0;i--) {
      GffObj& omrna=*mrnas[i];
      if (m->start>omrna.end) {
           if (m->start-omrna.start>GFF_MAX_EXON) break; //give up already
           continue;
           }
      if (omrna.start>m->end) continue; //this should never be the case if nidx was found correctly
      //locus overlap here:
      if (tMatch(*m, omrna, ovlen, false, true)) {
             dupidx=i;
             return mrnas[i];
             }
      }
  return NULL;
}


bool intronRedundant(GffObj& ti, GffObj&  tj) {
 //two transcripts are "intron redundant" iff one transcript's intron chain
  // is a sub-chain of the other's
 int imax=ti.exons.Count()-1;
 int jmax=tj.exons.Count()-1;
 if (imax==0 || jmax==0) return false; //don't deal with single-exon transcripts here
 if (ti.exons[imax]->start<tj.exons[0]->end ||
     tj.exons[jmax]->start<ti.exons[0]->end )
         return false; //intron chains do not overlap at all
 
 uint eistart=0, eiend=0, ejstart=0, ejend=0; //exon boundaries
 int i=1; //exon idx to the right of the current intron of ti
 int j=1; //exon idx to the right of the current intron of tj
 //find the first intron overlap:
 while (i<=imax && j<=jmax) {
    eistart=ti.exons[i-1]->end;
    eiend=ti.exons[i]->start;
    ejstart=tj.exons[j-1]->end;
    ejend=tj.exons[j]->start;
    if (ejend<eistart) { j++; continue; }
    if (eiend<ejstart) { i++; continue; }
    //we found an intron overlap
    break;
    }
 if ((i>1 && j>1) || i>imax || j>jmax) {
     return false; //either no intron overlaps found at all
                  //or it's not the first intron for at least one of the transcripts
     }
 if (eistart!=ejstart || eiend!=ejend) return false; //not an exact intron match
 //we have the first matching intron on the left
 if (j>i) {
   //i==1, ti's start must not conflict with the previous intron of tj
   if (ti.start<tj.exons[j-1]->start) return false;
   //so i's first intron starts AFTER j's first intron
   // then j must contain i, so i's last intron must end with or before j's last intron
   if (ti.exons[imax]->start>tj.exons[jmax]->start) return false;
      //comment out the line above if you just want "intron compatibility" (i.e. extension of intron chains )
   }
  else if (i>j) {
   //j==1, tj's start must not conflict with the previous intron of ti
   if (tj.start<ti.exons[i-1]->start) return false;
   //so j's intron chain starts AFTER i's
   // then i must contain j, so j's last intron must end with or before j's last intron
   if (tj.exons[jmax]->start>ti.exons[imax]->start) return false;
      //comment out the line above for just "intronCompatible()" check
   }
 //now check if the rest of the introns overlap, in the same sequence
 i++;
 j++;
 while (i<=imax && j<=jmax) {
  if (ti.exons[i-1]->end!=tj.exons[j-1]->end ||
      ti.exons[i]->start!=tj.exons[j]->start) return false;
  i++;
  j++;
  }
 i--;
 j--;
 if (i==imax && j<jmax) {
   // tj has more introns to the right, check if ti's end doesn't conflict with the current tj exon boundary
   if (ti.end>tj.exons[j]->end) return false;
   }
 else if (j==jmax && i<imax) {
   if (tj.end>ti.exons[i]->end) return false;
   }
 return true;
}

bool t_contains(GffObj& a, GffObj& b) {
 //returns true if b's intron chain (or single exon) is included in a
 if (b.exons.Count()>=a.exons.Count()) return false;
 if (b.exons.Count()==1) {
    //check if b is contained in any of a's exons:
    for (int i=0;i<a.exons.Count();i++) {
       if (b.start>=a.exons[i]->start && b.end<=a.exons[i]->end) return true;
       }
    return false;
    }
 if (intronRedundant(a,b)) {
    //intronRedudant allows b's initial/terminal exons to extend beyond a's boundaries
    //but we don't allow this kind of behavior here
    return (b.start>=a.start && b.end<=a.end);
    }
  else return false;
 }

int is_Redundant(GffObj*m, GList<GffObj>* mrnas) {
 //first locate the list index of the mrna starting just ABOVE
 //the end of this mrna
 if (mrnas->Count()==0) return -1;
 int nidx=qsearch_mrnas(m->end, *mrnas);
 if (nidx==0) return -1;
 if (nidx==-1) nidx=mrnas->Count();//all can overlap
 for (int i=nidx-1;i>=0;i--) {
     GffObj& omrna=*mrnas->Get(i);
     if (m->start>omrna.end) {
          if (m->start-omrna.start>GFF_MAX_EXON) break; //give up already
          continue;
          }
     if (omrna.start>m->end) continue; //this should never be the case if nidx was found correctly
     
     if (intronRedundant(*m, omrna)) return i;
     }
 return -1;
}

bool t_dominates(GffObj* a, GffObj* b) {
 // for redundant / intron compatible transfrags:
 // returns true if a "dominates" b, i.e. a has more exons or is longer
 if (a->exons.Count()==b->exons.Count())
         return (a->covlen>b->covlen);
    else return (a->exons.Count()>b->exons.Count());
}

bool betterDupRef(GffObj* a, GffObj* b) {
  if (a->exons.Count()!=b->exons.Count())
    return (a->exons.Count()>b->exons.Count());
  if (a->hasCDS()!=b->hasCDS())
     return (a->hasCDS()>b->hasCDS());
   //for annotation purposes, it's more important to keep the 
   //longer transcript, instead of the one that was loaded first
  if (a->covlen != b->covlen)
         return (a->covlen > b->covlen);
    else return (a->track_id < b->track_id);
}

int parse_mRNAs(GfList& mrnas,
				 GList<GSeqData>& glstdata,
				 bool is_ref_set,
				 bool check_for_dups,
				 int qfidx, bool only_multiexon) {
	int refdiscarded=0; //ref duplicates discarded
	int tredundant=0; //cufflinks redundant transcripts discarded
	for (int k=0;k<mrnas.Count();k++) {
		GffObj* m=mrnas[k];
		int i=-1;
		GSeqData f(m->gseq_id);
		GSeqData* gdata=NULL;
		uint tlen=m->len();
		if (m->hasErrors() || (tlen+500>GFF_MAX_LOCUS)) { //should probably report these in a file too..
			if (gtf_tracking_verbose) 
			      GMessage("Warning: transcript %s discarded (structural errors found, length=%d).\n", m->getID(), tlen);
			continue;
			}
		if (only_multiexon && m->exons.Count()<2) {
			continue;
			}
		//GStr feature(m->getFeatureName());
		//feature.lower();
		//bool gene_or_locus=(feature.endsWith("gene") ||feature.index("loc")>=0);
		//if (m->exons.Count()==0 && gene_or_locus) {
		if (m->isDiscarded()) {
			//discard generic "gene" or "locus" features with no other detailed subfeatures
			//if (gtf_tracking_verbose)
			//   GMessage("Warning: discarding GFF generic gene/locus container %s\n",m->getID());
			continue;
			}
		if (m->exons.Count()==0) {
				//if (gtf_tracking_verbose)
				// GMessage("Warning: %s %s found without exon segments (adding default exon).\n",m->getFeatureName(), m->getID());
				m->addExon(m->start,m->end);
				}
		if (glstdata.Found(&f,i)) gdata=glstdata[i];
		else {
			gdata=new GSeqData(m->gseq_id);
			glstdata.Add(gdata);
			}
		
		double fpkm=0;
		double cov=0;
		double conf_hi=0;
		double conf_lo=0;

		GList<GffObj>* target_mrnas=NULL;
		if (is_ref_set) { //-- ref transcripts
		   if (m->strand=='.') {
		     //unknown strand - discard from reference set (!)
		     continue;
		     }
		   target_mrnas=(m->strand=='+') ? &(gdata->mrnas_f) : &(gdata->mrnas_r);
		   if (check_for_dups) {
		     //check all gdata->mrnas_r (ref_data) for duplicate ref transcripts
		     int rpidx=-1;
		     GffObj* rp= is_RefDup(m, *target_mrnas, rpidx);
		     if (rp!=NULL) { //duplicate found
		      //discard one of them
		      //but let's keep the gene_name if present
		      //DEBUG:
		      //GMessage("Ref duplicates: %s = %s\n", rp->getID(), m->getID());
		      refdiscarded++;
		      if (betterDupRef(rp, m)) {
		           if (rp->getGeneName()==NULL && m->getGeneName()!=NULL) {
		                  rp->setGeneName(m->getGeneName());
		                  }
		           continue;
		           }
		         else {
		           if (m->getGeneName()==NULL && rp->getGeneName()!=NULL) {
		                  m->setGeneName(rp->getGeneName());
		                  }
		           ((CTData*)(rp->uptr))->mrna=NULL;
		           rp->isUsed(false);
		           target_mrnas->Forget(rpidx);
		           target_mrnas->Delete(rpidx);
		           }
		       }
		     } //check for duplicate ref transcripts
		   } //ref transcripts
		else { //-- transfrags
		   if (m->strand=='+') { target_mrnas = &(gdata->mrnas_f); }
		     else if (m->strand=='-') { target_mrnas=&(gdata->mrnas_r); }
		        else { m->strand='.'; target_mrnas=&(gdata->umrnas); }
		   if (check_for_dups) { //check for redundancy
		     // check if there is a redundancy between this and another already loaded Cufflinks transcript
		     int cidx =  is_Redundant(m, target_mrnas);
		     if (cidx>=0) {
		        //always discard the redundant transcript with the fewer exons OR shorter
		        if (t_dominates(target_mrnas->Get(cidx),m)) {
		            //new transcript is shorter, discard it
		            continue;
		            } 
		        else {
		            //discard the older transfrag
		            ((CTData*)(target_mrnas->Get(cidx)->uptr))->mrna=NULL;
		            target_mrnas->Get(cidx)->isUsed(false);
		            target_mrnas->Forget(cidx);
		            target_mrnas->Delete(cidx);
		            //the uptr (CTData) pointer will still be kept in gdata->ctdata and deallocated eventually
		            }
		        tredundant++;
		        }
		     }// redundant transfrag check
		   if (m->gscore==0.0)   
		     m->gscore=m->exons[0]->score; //Cufflinks exon score = isoform abundance
		   //const char* expr = (gtf_tracking_largeScale) ? m->getAttr("FPKM") : m->exons[0]->getAttr(m->names,"FPKM");
		   const char* expr = m->getAttr("FPKM");
		   if (expr!=NULL) {
		       if (expr[0]=='"') expr++;
		       fpkm=strtod(expr, NULL);
		       } else { //backward compatibility: read RPKM if FPKM not found
		       //expr=(gtf_tracking_largeScale) ? m->getAttr("RPKM") : m->exons[0]->getAttr(m->names,"RPKM");
		       expr=m->getAttr("RPKM");
		       if (expr!=NULL) {
		           if (expr[0]=='"') expr++;
		           fpkm=strtod(expr, NULL);
		           }
		       }
		   //const char* scov=(gtf_tracking_largeScale) ? m->getAttr("cov") : m->exons[0]->getAttr(m->names,"cov");
		   const char* scov=m->getAttr("cov");
		   if (scov!=NULL) {
		       if (scov[0]=='"') scov++; 
		       cov=strtod(scov, NULL);
		       }
		   //const char* sconf_hi=(gtf_tracking_largeScale) ? m->getAttr("conf_hi") : m->exons[0]->getAttr(m->names,"conf_hi");
		   const char* sconf_hi=m->getAttr("conf_hi");
		   if (sconf_hi!=NULL){
		       if (sconf_hi[0]=='"') sconf_hi++;
		       conf_hi=strtod(sconf_hi, NULL);
		       }
		   //const char* sconf_lo=(gtf_tracking_largeScale) ? m->getAttr("conf_lo") : m->exons[0]->getAttr(m->names,"conf_lo");
		   const char* sconf_lo=m->getAttr("conf_lo");
		   if (sconf_lo!=NULL) {
		       if (sconf_lo[0]=='"') sconf_lo++;
		       conf_lo=strtod(sconf_lo, NULL);
		       }
		   } //Cufflinks transfrags
		target_mrnas->Add(m);
		m->isUsed(true);
		CTData* mdata=new CTData(m);
		mdata->qset=qfidx;
		gdata->tdata.Add(mdata);
		if (!is_ref_set) {
		// Cufflinks - attributes parsing
		   mdata->FPKM=fpkm;
		   mdata->cov=cov;
		   mdata->conf_hi=conf_hi;
		   mdata->conf_lo=conf_lo;
		   }
	}//for each mrna read
 //if (mrna_deleted>0)
 //  mrnas.Pack();
 
 return (is_ref_set ? refdiscarded : tredundant);
}

bool tMatch(GffObj& a, GffObj& b, int& ovlen, bool fuzzunspl, bool contain_only) {
	//strict intron chain match, or single-exon perfect match
	int imax=a.exons.Count()-1;
	int jmax=b.exons.Count()-1;
	ovlen=0;
	if (imax!=jmax) return false; //different number of introns
	if (imax==0) { //single-exon mRNAs
		if (contain_only) {
		   return ((a.start>=b.start && a.end<=b.end) || 
		           (b.start>=a.start && b.end<=a.end));
		}
		if (fuzzunspl) {
			//fuzz match for single-exon transfrags: 
			// it's a match if they overlap at least 80% of shortest one
			ovlen=a.exons[0]->overlapLen(b.exons[0]);
			int maxlen=GMAX(a.covlen,b.covlen);
			return (ovlen>=maxlen*0.8);
		}
	  else {
			//only exact match, or strictly contained
			ovlen=a.covlen;
			return (a.exons[0]->start==b.exons[0]->start &&
					a.exons[0]->end==b.exons[0]->end);
		}
	}
	if ( a.exons[imax]->start<b.exons[0]->end ||
		b.exons[jmax]->start<a.exons[0]->end )
		return false; //intron chains do not overlap at all
	//check intron overlaps
	ovlen=a.exons[0]->end-(GMAX(a.start,b.start))+1;
	ovlen+=(GMIN(a.end,b.end))-a.exons.Last()->start;
	for (int i=1;i<=imax;i++) {
		if (i<imax) ovlen+=a.exons[i]->len();
		if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
			(a.exons[i]->start!=b.exons[i]->start)) {
			return false; //intron mismatch
		}
	}
	if (contain_only) 
		     return ((a.start>=b.start && a.end<=b.end) || 
		           (b.start>=a.start && b.end<=a.end));
		else return true;
}


void cluster_mRNAs(GList<GffObj> & mrnas, GList<GLocus> & loci, int qfidx) {
	//mrnas sorted by start coordinate
	//and so are the loci
	//int rdisc=0;
		for (int t=0;t<mrnas.Count();t++) {
		GArray<int> mrgloci(false);
		GffObj* mrna=mrnas[t];
		int lfound=0; //count of parent loci
		/*for (int l=0;l<loci.Count();l++) {
			if (loci[l]->end<mrna->exons.First()->start) continue;
			if (loci[l]->start>mrna->exons.Last()->end) break; */
		 for (int l=loci.Count()-1;l>=0;l--) {
		   if (loci[l]->end<mrna->exons.First()->start) {
		       if (mrna->exons.First()->start-loci[l]->start > GFF_MAX_LOCUS) break;
		       continue;
		       }
		   if (loci[l]->start>mrna->exons.Last()->end) continue;
			//here we have mrna overlapping loci[l]
			if (loci[l]->add_mRNA(mrna)) {
				//a parent locus was found
				lfound++;
				mrgloci.Add(l); //locus indices added here, in decreasing order
			}
		}//loci loop
		//if (lfound<0) continue; //mrna was a ref duplicate, skip it
		if (lfound==0) {
			//create a locus with only this mRNA
 			 loci.Add(new GLocus(mrna, qfidx));
		    }
		 else if (lfound>1) {
			//more than one locus found parenting this mRNA, merge loci
		     lfound--;
			 for (int l=0;l<lfound;l++) {
				  int mlidx=mrgloci[l]; //largest indices first, so it's safe to remove
				  loci[mrgloci[lfound]]->addMerge(*loci[mlidx], mrna);
				  loci.Delete(mlidx);
			    }
		    }
	}//mrnas loop
	//if (rdisc>0) mrnas.Pack();
	//return rdisc;
}

int fix_umrnas(GSeqData& seqdata, GSeqData* rdata, FILE* fdis=NULL) {
	//attempt to find the strand for seqdata.umrnas
	//based on a) overlaps with oriented reference mRNAs if present
	//         b) overlaps with oriented mRNAs from the same input set
	if (rdata!=NULL) { //we have reference mrnas
		for (int i=0;i<rdata->mrnas_f.Count();i++) {
			for (int j=0;j<seqdata.umrnas.Count();j++) {
				if (rdata->mrnas_f[i]->gseq_id!=seqdata.umrnas[j]->gseq_id) continue;
				if (seqdata.umrnas[j]->strand!='.') continue;
				uint ustart=seqdata.umrnas[j]->exons.First()->start;
				uint uend=seqdata.umrnas[j]->exons.Last()->end;
				uint rstart=rdata->mrnas_f[i]->exons.First()->start;
				uint rend=rdata->mrnas_f[i]->exons.Last()->end;
				if (ustart>rend) break;
				if (rstart>uend) continue;
				if (rdata->mrnas_f[i]->exonOverlap(ustart,uend)) {
					seqdata.umrnas[j]->strand='+';
				}
				else { //within intron
					//if (seqdata.umrnas[j]->ulink==NULL ||
					//     seqdata.umrnas[j]->ulink->covlen<rdata->mrnas_f[i]->covlen) {
					CTData* mdata=(CTData*)seqdata.umrnas[j]->uptr;
					mdata->addOvl('i',rdata->mrnas_f[i]);
				}
			}
		}
		for (int i=0;i<rdata->mrnas_r.Count();i++) {
			for (int j=0;j<seqdata.umrnas.Count();j++) {
				if (seqdata.umrnas[j]->strand!='.') continue;
				uint ustart=seqdata.umrnas[j]->exons.First()->start;
				uint uend=seqdata.umrnas[j]->exons.Last()->end;
				uint rstart=rdata->mrnas_r[i]->exons.First()->start;
				uint rend=rdata->mrnas_r[i]->exons.Last()->end;
				if (ustart>rend) break;
				if (rstart>uend) continue;
				if (rdata->mrnas_r[i]->exonOverlap(ustart,uend)) {
					seqdata.umrnas[j]->strand='-';
				}
				else { //within intron
					CTData* mdata=(CTData*)seqdata.umrnas[j]->uptr;
					mdata->addOvl('i',rdata->mrnas_r[i]);
				}
				
			}
		}
	}//we have reference transcripts
	//---- now compare to other transcripts
	for (int i=0;i<seqdata.mrnas_f.Count();i++) {
		for (int j=0;j<seqdata.umrnas.Count();j++) {
			if (seqdata.umrnas[j]->strand!='.') continue;
			uint ustart=seqdata.umrnas[j]->exons.First()->start;
			uint uend=seqdata.umrnas[j]->exons.Last()->end;
			uint rstart=seqdata.mrnas_f[i]->exons.First()->start;
			uint rend=seqdata.mrnas_f[i]->exons.Last()->end;
			if (ustart>rend) break;
			if (rstart>uend) continue;
			if (seqdata.mrnas_f[i]->exonOverlap(ustart,uend)) {
				seqdata.umrnas[j]->strand='+';
			}
		}
	}
	for (int i=0;i<seqdata.mrnas_r.Count();i++) {
		for (int j=0;j<seqdata.umrnas.Count();j++) {
			if (seqdata.umrnas[j]->strand!='.') continue;
			uint ustart=seqdata.umrnas[j]->exons.First()->start;
			uint uend=seqdata.umrnas[j]->exons.Last()->end;
			uint rstart=seqdata.mrnas_r[i]->exons.First()->start;
			uint rend=seqdata.mrnas_r[i]->exons.Last()->end;
			if (ustart>rend) break;
			if (rstart>uend) continue;
			//overlap
			if (seqdata.mrnas_r[i]->exonOverlap(ustart,uend)) {
				seqdata.umrnas[j]->strand='-';
			}
		}
    }
	int fcount=0;
	for (int i=0;i<seqdata.umrnas.Count();i++) {
		if (seqdata.umrnas[i]->strand=='+') {
			seqdata.mrnas_f.Add(seqdata.umrnas[i]);
			seqdata.umrnas.Forget(i);
		}
		else if (seqdata.umrnas[i]->strand=='-') {
		    seqdata.mrnas_r.Add(seqdata.umrnas[i]);
		    seqdata.umrnas.Forget(i);
		}
		else {  //discard mRNAs not settled
			seqdata.umrnas[i]->strand='.';
			if (fdis!=NULL) {
				seqdata.umrnas[i]->printGtf(fdis);
				}
			fcount++;
		}
	}
	seqdata.umrnas.Pack();
	return fcount;
}

//retrieve ref_data for a specific genomic sequence
GSeqData* getRefData(int gid, GList<GSeqData>& ref_data) {
	int ri=-1;
	GSeqData f(gid);
	GSeqData* r=NULL;
	if (ref_data.Found(&f,ri))
		r=ref_data[ri];
	return r;
}

void read_transcripts(FILE* f, GList<GSeqData>& seqdata, bool keepAttrs) {
	rewind(f);
	GffReader gffr(f, true); //loading only recognizable transcript features
	gffr.showWarnings(gtf_tracking_verbose);

	//          keepAttrs    mergeCloseExons   noExonAttrs
	gffr.readAll(keepAttrs,          true,        true);

	//                               is_ref?    check_for_dups,
	parse_mRNAs(gffr.gflst, seqdata, false,       false);
}

int cmpGSeqByName(const pointer p1, const pointer p2) {
 return strcmp(((GSeqData*)p1)->gseq_name, ((GSeqData*)p2)->gseq_name);
}

void sort_GSeqs_byName(GList<GSeqData>& seqdata) {
  seqdata.setSorted(&cmpGSeqByName);
}

void read_mRNAs(FILE* f, GList<GSeqData>& seqdata, GList<GSeqData>* ref_data,
	         bool check_for_dups, int qfidx, const char* fname, bool only_multiexon) {
	//>>>>> read all transcripts/features from a GTF/GFF3 file
	//int imrna_counter=0;
	int loci_counter=0;
	if (ref_data==NULL) ref_data=&seqdata;
	bool isRefData=(&seqdata==ref_data);
	                          //(f, transcripts_only)
	GffReader* gffr=new GffReader(f, true); //load only transcript annotations
	gffr->showWarnings(gtf_tracking_verbose);
	//            keepAttrs   mergeCloseExons   noExonAttrs
	gffr->readAll(!isRefData,          true,        isRefData || gtf_tracking_largeScale);
	//so it will read exon attributes only for low number of Cufflinks files
	
	int d=parse_mRNAs(gffr->gflst, seqdata, isRefData, check_for_dups, qfidx,only_multiexon);
	if (gtf_tracking_verbose && d>0) {
	  if (isRefData) GMessage(" %d duplicate reference transcripts discarded.\n",d);
	            else GMessage(" %d redundant cufflinks transfrags discarded.\n",d);
	  }
	//imrna_counter=gffr->mrnas.Count();
	delete gffr; //free the extra memory and unused GffObjs
	
	//for each genomic sequence, cluster transcripts
	int discarded=0;
	GStr bname(fname);
	GStr s;
	if (!bname.is_empty()) {
		int di=bname.rindex('.');
		if (di>0) bname.cut(di);
		int p=bname.rindex('/');
		if (p<0) p=bname.rindex('\\');
		if (p>=0) bname.remove(0,p);
	}
	FILE* fdis=NULL;
	FILE* frloci=NULL;

	for (int g=0;g<seqdata.Count();g++) {
		//find the corresponding refseqdata with the same gseq_id
		int gseq_id=seqdata[g]->get_gseqid();
		if (!isRefData) { //cufflinks data, find corresponding ref data
			GSeqData* rdata=getRefData(gseq_id, *ref_data);
			if (rdata!=NULL && seqdata[g]->umrnas.Count()>0) {
			    discarded+=fix_umrnas(*seqdata[g], rdata, fdis);
			    }
			}
		//>>>>> group mRNAs into locus-clusters (based on exon overlap)
		cluster_mRNAs(seqdata[g]->mrnas_f, seqdata[g]->loci_f, qfidx);
		cluster_mRNAs(seqdata[g]->mrnas_r, seqdata[g]->loci_r, qfidx);
		if (!isRefData) {
			cluster_mRNAs(seqdata[g]->umrnas, seqdata[g]->nloci_u, qfidx);
			}
		loci_counter+=seqdata[g]->loci_f.Count();
		loci_counter+=seqdata[g]->loci_r.Count();
//		if (refData) {
//			if (frloci==NULL) {
//				s=bname;
//				s.append(".loci.lst");
//				frloci=fopen(s.chars(), "w");
//			}
//			writeLoci(frloci, seqdata[g]->loci_f);
//			writeLoci(frloci, seqdata[g]->loci_r);
//		}//write ref loci
	}//for each genomic sequence
	if (fdis!=NULL) fclose(fdis);
	if (frloci!=NULL) fclose(frloci);
	if (discarded>0) {
		if (gtf_tracking_verbose) GMessage("Found %d transcripts with undetermined strand.\n", discarded);
	}
	else { if (fdis!=NULL) remove(s.chars()); }
}

int qsearch_mrnas(uint x, GList<GffObj>& mrnas) {
  //binary search
  //do the simplest tests first:
  if (mrnas[0]->start>x) return 0;
  if (mrnas.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=mrnas.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=mrnas[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && mrnas[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && mrnas[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int qsearch_loci(uint x, GList<GLocus>& loci) {
 // same as above, but for GSeg lists
  //binary search
  //do the simplest tests first:
  if (loci[0]->start>x) return 0;
  if (loci.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=loci.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l + h) >> 1;
     istart=loci[i]->start;
     if (istart < x) l=i+1;
                else {
                   if (istart == x) { //found matching coordinate here
                        idx=i;
                        while (idx<=maxh && loci[idx]->start==x) {
                           idx++;
                           }
                        return (idx>maxh) ? -1 : idx;
                        }
                   h=i-1;
                   }
     } //while
 idx = l;
 while (idx<=maxh && loci[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

