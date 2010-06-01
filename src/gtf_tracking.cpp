/*
 *  gtf_tracking.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/5/09.
 *  Copyright 2009 Geo Pertea. All rights reserved.
 *
 */

#include <cstdlib>
#include "gtf_tracking.h"
#include "gff.h"
#include "GStr.h"

const char* ATTR_GENE_NAME=  "gene_name";

bool verbose = false;
bool largeScale=false; //many input Cufflinks files processed at once by cuffcompare, discard exon attributes

int GXConsensus::count=0;

int cmpByPtr(const pointer p1, const pointer p2) {
  return (p1>p2) ? 1: ((p1==p2)? 0 : -1);
  }


GffObj* is_mRNADup(GffObj* m, GList<GffObj>& mrnas) {
 //mrnas MUST be sorted by start coordinate
  int ovlen=0;
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
      if (tMatch(*m, omrna, ovlen)) return mrnas[i];
      }
  return NULL;
}


bool intronRedundant(GffObj& ti, GffObj&  tj) {
 //two transcripts are "intron redundant" iff one transcript's intron chain
  // is a sub-chain of the other's
 int imax=ti.exons.Count()-1;
 int jmax=tj.exons.Count()-1;
 if (imax==0 || jmax==0) return false; //don't care about single-exon transcripts
 if (ti.exons[imax]->start<tj.exons[0]->end ||
     tj.exons[jmax]->start<ti.exons[0]->end )
         return false; //intron chains do not overlap at all
 //find the first intron overlap:
 uint istart=0, iend=0, jstart=0, jend=0;
 int i=1; //exon idx to the right of the current intron of ti
 int j=1; //exon idx to the right of the current intron of tj
 while (i<=imax && j<=jmax) {
    istart=ti.exons[i-1]->end;
    iend=ti.exons[i]->start;
    jstart=tj.exons[j-1]->end;
    jend=tj.exons[j]->start;
    if (jend<istart) { j++; continue; }
    if (iend<jstart) { i++; continue; }
    //we found an intron overlap
    break;
    }
 if ((i>1 && j>1) || i>imax || j>jmax) {
     return false; //no intron overlaps found at all
                  //or not the first intron for at least one of the transcripts
     }
 if (istart!=jstart || iend!=jend) return false; //not an exact intron match
 if (j>i) {
   //i==1, ti's start must not conflict with the previous intron of tj
   if (ti.start<tj.exons[j-1]->start) return false;
   } else if (i>j) {
   //j==1, tj's start must not conflict with the previous intron of ti
   if (tj.start<ti.exons[i-1]->start) return false;
   }
 //now check if the rest of the introns overlap, in a linear succession
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

int parse_mRNAs(GList<GffObj>& mrnas,
				 GList<GSeqData>& glstdata,
				 bool is_ref_set,
				 bool check_for_dups,
				 int qfidx) {
    int refdiscarded=0; //ref duplicates discarded
    int tredundant=0; //cufflinks redundant transcripts discarded
    int mrna_deleted=0;
	for (int k=0;k<mrnas.Count();k++) {
		GffObj* m=mrnas[k];
		int i=-1;
		GSeqData f(m->gseq_id);
		GSeqData* gdata=NULL;
		uint tlen=m->end-m->start+1;

		if (m->hasErrors || (tlen+500>GFF_MAX_LOCUS)) { //should probably report these in a file too..
			GMessage("Warning: transcript %s discarded (structural errors found, length=%d).\n", m->getID(), tlen);
			delete m;
			mrnas.Forget(k);
			mrna_deleted++;
			continue;
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
		if (is_ref_set) {
		  if (check_for_dups) {
          //check all gdata->mrnas_r (ref_data) for duplicate ref transcripts
		  GffObj* rp= (m->strand=='+') ? is_mRNADup(m, gdata->mrnas_f) :
		                                 is_mRNADup(m, gdata->mrnas_r);
		  if (rp!=NULL) {
		    //just discard this, it's a duplicate ref
		     //but let's just keep the gene_name if present
		     char* gname=m->getAttr(ATTR_GENE_NAME);
		     char* pgname=rp->getAttr(ATTR_GENE_NAME);
		     if (pgname==NULL && gname!=NULL)
		         rp->addAttr(ATTR_GENE_NAME, gname);
		     delete m;
		     mrnas.Forget(k);
		     mrna_deleted++;
		     refdiscarded++;
		     continue;
		     }
		   } //suppress ref dups
		} //is ref 
		else { // Cufflinks gtf
		if (m->strand!='.' && check_for_dups) { //oriented transfrag
			// check if there is a redundancy between this and another already loaded Cufflinks transcript
			GList<GffObj>* ckmrnas=(m->strand=='+') ? &(gdata->mrnas_f) : &(gdata->mrnas_r);
			int cidx =  is_Redundant(m, ckmrnas);
			if (cidx>=0) {
				//always discard the "shorter" transcript of the redundant pair
				if (ckmrnas->Get(cidx)->covlen>m->covlen) {
				//new transcript is shorter, discard it
					delete m;
					mrnas.Forget(k);
					mrna_deleted++;
					continue;
				} else {
					//new transcript is longer, discard the older one
					((CTData*)(ckmrnas->Get(cidx)->uptr))->mrna=NULL;
					ckmrnas->Delete(cidx);
				//the uptr (CTData) pointer will still be kept in gdata->ctdata and freed accordindly at the end
					}
				tredundant++;
				}
			// ^^^ redundant transcript check
			} //oriented transfrag
		if (m->gscore==0.0)   
		   m->gscore=m->exons[0]->score; //Cufflinks exon score = isoform abundance
		//for Cufflinks file, parse expr attribute from the 1st exon (the lowest coordinate exon)
		const char* expr = (largeScale) ? m->getAttr("FPKM") : m->exons[0]->getAttr(m->names,"FPKM");
		if (expr!=NULL) {
			if (expr[0]=='"') expr++;
			fpkm=strtod(expr, NULL);
			} else { //backward compatibility: read RPKM if FPKM not found
			expr=(largeScale) ? m->getAttr("RPKM") : m->exons[0]->getAttr(m->names,"RPKM");
			if (expr!=NULL) {
				if (expr[0]=='"') expr++;
				fpkm=strtod(expr, NULL);
				}
			}
		const char* scov=(largeScale) ? m->getAttr("cov") : m->exons[0]->getAttr(m->names,"cov");
		if (scov!=NULL) {
			if (scov[0]=='"') scov++; 
			cov=strtod(scov, NULL);
			}
		const char* sconf_hi=(largeScale) ? m->getAttr("conf_hi") : m->exons[0]->getAttr(m->names,"conf_hi");
		if (sconf_hi!=NULL){
			if (sconf_hi[0]=='"') sconf_hi++;
			conf_hi=strtod(sconf_hi, NULL);
			}
		const char* sconf_lo=(largeScale) ? m->getAttr("conf_lo") : m->exons[0]->getAttr(m->names,"conf_lo");
		if (sconf_lo!=NULL) {
			if (sconf_lo[0]=='"') sconf_lo++;
			conf_lo=strtod(sconf_lo, NULL);
			}
		} //Cufflinks transfrags

		if (m->strand=='+') gdata->mrnas_f.Add(m);
		   else {
			 if (m->strand=='-') gdata->mrnas_r.Add(m);
			 else { //unknown strand, unoriented mRNA
				if (is_ref_set) {// discard these from reference
					delete m; //just free
					mrnas.Forget(k);
					mrna_deleted++;
					continue;
					}
				else {//store them in the unknown strand pile, to be analyzed later
					m->strand=0;
					gdata->umrnas.Add(m);
					}
				} //unknown strand
			} //- or unknown strand
		CTData* mdata=new CTData(m);
		mdata->qset=qfidx;
		gdata->tdata.Add(mdata);
		//if (!is_ref_set) { //Cufflinks - attributes parsing
		mdata->FPKM=fpkm;
		mdata->cov=cov;
		mdata->conf_hi=conf_hi;
		mdata->conf_lo=conf_lo;
		//} //cufflinks transcripts
		/* for (int ex=0;ex<m->exons.Count();ex++) {
		 delete m->exons[ex]->attrs; //just free this extra memory here
		 m->exons[ex]->attrs=NULL;
		 } */
	}//for each mrna read
 if (mrna_deleted>0) {
   mrnas.Pack();
   }
 return (is_ref_set ? refdiscarded : tredundant);
}

bool tMatch(GffObj& a, GffObj& b, int& ovlen, bool equnspl) {
	//strict intron chain match, or single-exon perfect match
	int imax=a.exons.Count()-1;
	int jmax=b.exons.Count()-1;
	ovlen=0;
	if (imax!=jmax) return false; //different number of introns
	if (imax==0) { //single-exon mRNAs
		//consider match if they overlap over 50% of max len
		if (equnspl) {
			ovlen=a.exons[0]->overlapLen(b.exons[0]);
			int maxlen=GMAX(a.covlen,b.covlen);
			return (ovlen>maxlen/2);
		}
        else {
            //only exact match
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
	//GMessage("tMatch found: %s == %s\n", a.getID(),b.getID());
	return true;
}

void cluster_mRNAs(GList<GffObj> & mrnas, GList<GLocus> & loci, int qfidx, bool refData) {
	//mrnas sorted by start coordinate
	//and so are the loci
	//int rdisc=0;
   //TODO: optimize this clustering pattern, using the fact that both mrnas and loci are sorted by start coordinate
   // (so perhaps we should scan loci in reverse and give up after certain distance ? but check for problems with mrgloci )
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
			//more than one loci found parenting this mRNA, merge loci
			//if (lfound>2) GMessage(" merging %d loci \n",lfound);
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
				if (seqdata.umrnas[j]->strand>0) continue;
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
				if (seqdata.umrnas[j]->strand) continue;
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
			if (seqdata.umrnas[j]->strand) continue;
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
			if (seqdata.umrnas[j]->strand) continue;
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

void read_mRNAs(FILE* f, GList<GSeqData>& seqdata, GList<GSeqData>* ref_data,
	         bool check_for_dups, int qfidx, const char* fname, bool checkseq) {
	//>>>>> read all the mRNAs from a file
	GffReader* gffr=new GffReader(f, true);
	//int imrna_counter=0;
	int loci_counter=0;
	if (ref_data==NULL) ref_data=&seqdata;
	bool isRefData=(&seqdata==ref_data);
	//           keepAttrs   mergeCloseExons   noExonAttrs
	gffr->readAll(true,          true,        isRefData || largeScale);
	//so it will read exon attributes if low number of Cufflinks files
	
	int d=parse_mRNAs(gffr->gflst, seqdata, isRefData, check_for_dups, qfidx);
	if (verbose && d>0) {
	  if (isRefData) GMessage(" %d duplicate reference transcripts discarded.\n",d);
	             else GMessage(" %d redundant cufflinks transfrags discarded.\n",d);
	  }
	//imrna_counter=gffr->mrnas.Count();
	delete gffr; //free the extra memory
	
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
		int gseq_id=seqdata[g]->gseq_id;
		//if (verbose) GMessage("Clustering mRNAs on %s...\n", getGSeqName(gseq_id));
		if (!isRefData) { //cufflinks data, find corresponding ref data
			GSeqData* rdata=getRefData(gseq_id, *ref_data);
			if (rdata!=NULL && seqdata[g]->umrnas.Count()>0) {
				discarded+=fix_umrnas(*seqdata[g], rdata, fdis);
			    }
		    }
		//>>>>> group mRNAs into locus-clusters (based on exon overlap)
		cluster_mRNAs(seqdata[g]->mrnas_f, seqdata[g]->loci_f, qfidx, isRefData);
		cluster_mRNAs(seqdata[g]->mrnas_r, seqdata[g]->loci_r, qfidx, isRefData);
		if (!isRefData) {
			cluster_mRNAs(seqdata[g]->umrnas, seqdata[g]->nloci_u, qfidx, false);
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
		GMessage("Warning: found %d transcripts with undetermined strand.\n", discarded);
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

int qsearch_segs(uint x, GList<GSeg>& segs) {
 // same as above, but for GSeg lists
  //binary search
  //do the simplest tests first:
  if (segs[0]->start>x) return 0;
  if (segs.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=segs.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l + h) >> 1;
     istart=segs[i]->start;
     if (istart < x) l=i+1;
                else {
                   if (istart == x) { //found matching coordinate here
                        idx=i;
                        while (idx<=maxh && segs[idx]->start==x) {
                           idx++;
                           }
                        return (idx>maxh) ? -1 : idx;
                        }
                   h=i-1;
                   }
     } //while
 idx = l;
 while (idx<=maxh && segs[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

