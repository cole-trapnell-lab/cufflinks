#include "gff_utils.h"

extern bool verbose;
extern bool debugMode;

void printFasta(FILE* f, GStr& defline, char* seq, int seqlen) {
 if (seq==NULL) return;
 int len=(seqlen>0)?seqlen:strlen(seq);
 if (len<=0) return;
 if (!defline.is_empty())
     fprintf(f, ">%s\n",defline.chars());
 int ilen=0;
 for (int i=0; i < len; i++, ilen++) {
   if (ilen == 70) {
     fputc('\n', f);
     ilen = 0;
     }
   putc(seq[i], f);
   } //for
 fputc('\n', f);
}

int qsearch_gloci(uint x, GList<GffLocus>& loci) {
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
     i = (l+h)>>1;
     istart=loci[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && loci[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && loci[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int qsearch_rnas(uint x, GList<GffObj>& rnas) {
  //binary search
  //do the simplest tests first:
  if (rnas[0]->start>x) return 0;
  if (rnas.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=rnas.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=rnas[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && rnas[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && rnas[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int cmpRedundant(GffObj& a, GffObj& b) {
  if (a.exons.Count()==b.exons.Count()) {
     if (a.covlen==b.covlen) {
       return strcmp(a.getID(), b.getID());
       }
     else return (a.covlen>b.covlen)? 1 : -1;
     }
   else return (a.exons.Count()>b.exons.Count())? 1: -1;
}


bool tMatch(GffObj& a, GffObj& b) {
  //strict intron chain match, or single-exon perfect match
  int imax=a.exons.Count()-1;
  int jmax=b.exons.Count()-1;
  int ovlen=0;
  if (imax!=jmax) return false; //different number of introns
  if ( a.exons[imax]->start<b.exons[0]->end ||
   b.exons[jmax]->start<a.exons[0]->end )
   return false; //intron chains do not overlap at all

  if (imax==0) { //single-exon mRNAs
    //if (equnspl) {
      //fuzz match for single-exon transfrags: 
      // it's a match if they overlap at least 80% of max len
      ovlen=a.exons[0]->overlapLen(b.exons[0]);
      int maxlen=GMAX(a.covlen,b.covlen);
      return (ovlen>=maxlen*0.8);
    /*}
    else {
      //only exact match
      ovlen=a.covlen;
      return (a.exons[0]->start==b.exons[0]->start &&
          a.exons[0]->end==b.exons[0]->end);
      
       }*/
     }
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
  return true;
}


bool unsplContained(GffObj& ti, int imax, GffObj&  tj, int jmax) {
 //returns true only if ti (which MUST be single-exon) is "almost" contained in any of tj's exons
 //but it does not cross any intron-exon boundary of tj!
  if (imax>0) GError("Error: bad unsplContained() call, 1st param must be single-exon\n");
  int tilen=ti.exons[0]->len();
  for (int j=0;j<=jmax;j++) {
       //must NOT overlap the introns
       if ((j>0 && ti.start<tj.exons[j]->start) 
          || (j<jmax && ti.end>tj.exons[j]->end))
         return false;
       if (ti.exons[0]->overlapLen(tj.exons[j])>=tilen*0.8) {
              return true;
              }
       }
 return false;
}

bool redundantTranscripts(GffObj& ti, GffObj&  tj, bool fullmatchonly) {
  //two transcripts are "intron redundant" iff one transcript's intron chain
  // is a sub-chain of the other's
  
 int imax=ti.exons.Count()-1;
 int jmax=tj.exons.Count()-1;
  if (ti.exons[imax]->start<tj.exons[0]->end ||
     tj.exons[jmax]->start<ti.exons[0]->end )
         return false; //intron chains do not overlap at all
 if (fullmatchonly) { 
   if (imax!=jmax) return false;
   if (imax==0) {
     //single-exon transcripts: at least 80% of the longest one is overlapped by the other
     int ovlen=ti.exons[0]->overlapLen(tj.exons[0]);
     int maxlen=GMAX(ti.covlen,tj.covlen);
     return (ovlen>=maxlen*0.8);
     }
   for (int i=1;i<=imax;i++) {
        if ((ti.exons[i-1]->end!=tj.exons[i-1]->end) ||
            (ti.exons[i]->start!=tj.exons[i]->start)) {
            return false; //intron mismatch
            }
        }
   return true;
   }
 //containment is also considered redundancy
 if (imax==0) {
   //check if this single exon is contained in any of tj exons
   return unsplContained(ti, imax, tj, jmax);
   }
 if (jmax==0) {
   return unsplContained(tj,jmax,ti,imax);
   }
 //full intron chain containment is also considered redundancy
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
      //comment out the line above for just "intronCompatible()" check (allowing extension of intron chain)
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


int gseqCmpName(const pointer p1, const pointer p2) {
 return strcmp(((GenomicSeqData*)p1)->gseq_name, ((GenomicSeqData*)p2)->gseq_name);
}


void printLocus(GffLocus* loc, const char* pre=NULL) {
  if (pre!=NULL) fprintf(stderr, "%s", pre);
  GMessage(" [%d-%d] : ", loc->start, loc->end);
  GMessage("%s",loc->rnas[0]->getID());
  for (int i=1;i<loc->rnas.Count();i++) {
    GMessage(",%s",loc->rnas[i]->getID());
    }
  GMessage("\n");
}

void placeGf(GffObj* t, GenomicSeqData* gdata, bool doCluster, bool collapseRedundant, bool collapseContained) {
  //GMessage(">>Placing transcript %s\n", t->getID());
  GTData* tdata=new GTData(t);
  gdata->tdata.Add(tdata);
  int tidx=-1;
  if (t->exons.Count()>0)
              tidx=gdata->rnas.Add(t); //added it in sorted order
            else {
              gdata->gfs.Add(t);
              return; //nothing to do with these non-transcript objects
              }
  if (!doCluster) return;
  if (gdata->loci.Count()==0) {
       gdata->loci.Add(new GffLocus(t));
       //GMessage("  <<make it first locus %d-%d \n",t->start, t->end);
       return;
       }
  //DEBUG: show available loci:
  // GMessage("  [%d loci already:\n", gdata->loci.Count());
  //for (int l=0;l<gdata->loci.Count();l++) {
  //    printLocus(gdata->loci[l]);
  //    }
  int nidx=qsearch_gloci(t->end, gdata->loci); //get index of nearest locus starting just ABOVE t->end
  //GMessage("\tlooking up end coord %d in gdata->loci.. (qsearch got nidx=%d)\n", t->end, nidx);
  if (nidx==0) {
     //cannot have any overlapping loci
     //GMessage("  <<no ovls possible, create locus %d-%d \n",t->start, t->end);
     gdata->loci.Add(new GffLocus(t));
     return;
     }
  if (nidx==-1) nidx=gdata->loci.Count();//all loci start below t->end
  int lfound=0; //count of parent loci
  GArray<int> mrgloci(false);
  GList<GffLocus> tloci(true); //candidate parent loci to adopt this
  //GMessage("\tchecking all loci from %d to 0\n",nidx-1);
  for (int l=nidx-1;l>=0;l--) {
      GffLocus& loc=*(gdata->loci[l]);
      if (loc.strand!='.' && t->strand!='.'&& loc.strand!=t->strand) continue;
      if (t->start>loc.end) {
           if (t->start-loc.start>GFF_MAX_LOCUS) break; //give up already
           continue;
           }
      if (loc.start>t->end) continue;
          //this should never be the case if nidx was found correctly
      //GMessage(" !range overlap found with locus ");
      //printLocus(&loc);
      if (loc.add_RNA(t)) {
         //will add this transcript to loc
         lfound++;
         mrgloci.Add(l);
         if (collapseRedundant) {
           //compare to every single transcript in this locus
           for (int ti=0;ti<loc.rnas.Count();ti++) {
                 if (loc.rnas[ti]==t) continue;
                 GTData* odata=(GTData*)(loc.rnas[ti]->uptr);
                 //GMessage("  ..redundant check vs overlapping transcript %s\n",loc.rnas[ti]->getID());
                 if (odata->replaced_by==NULL && redundantTranscripts(*t, *(loc.rnas[ti]), !collapseContained)) {
                     if (cmpRedundant(*t, *(loc.rnas[ti]))>0) {
                        odata->replaced_by=t;
                        }
                     else {
                        tdata->replaced_by=loc.rnas[ti];
                        }
                     }
              }//for each transcript in the exon-overlapping locus
          } //if doCollapseRedundant
         } //overlapping locus
      }
  if (lfound==0) {
      //overlapping loci not found, create a locus with only this mRNA
      //GMessage("  overlapping locus not found, create locus %d-%d \n",t->start, t->end);
      int addidx=gdata->loci.Add(new GffLocus(t));
      if (addidx<0) {
         GMessage("  WARNING: new GffLocus(%s:%d-%d) not added!\n",t->getID(), t->start, t->end);
         }
      }
   else if (lfound>1) {
      //more than one loci found parenting this mRNA, merge loci
      //if (lfound>2) GMessage(" merging %d loci \n",lfound);
       lfound--;
       for (int l=0;l<lfound;l++) {
          int mlidx=mrgloci[l]; //largest indices first, so it's safe to remove
          gdata->loci[mrgloci[lfound]]->addMerge(*(gdata->loci[mlidx]), t);
          gdata->loci.Delete(mlidx);
          }
      }
}

void collectLocusData(GList<GenomicSeqData>& ref_data) {
  int locus_num=0;
  for (int g=0;g<ref_data.Count();g++) {
    GenomicSeqData* gdata=ref_data[g];
    for (int l=0;l<gdata->loci.Count();l++) {
      GffLocus& loc=*(gdata->loci[l]);
      GHash<int> gnames(true); //gene names in this locus
      GHash<int> geneids(true); //Entrez GeneID: numbers
      for (int i=0;i<loc.rnas.Count();i++) {
        GffObj& t=*(loc.rnas[i]);
        GStr gname(t.getGeneName());
        if (!gname.is_empty()) {
           gname.upper();
           int* prevg=gnames.Find(gname.chars());
           if (prevg!=NULL) (*prevg)++;
                  else gnames.Add(gname, new int(1));
           }
        //parse GeneID xrefs, if any:
        GStr xrefs(t.getAttr("xrefs"));
        if (!xrefs.is_empty()) {
          xrefs.startTokenize(",");
          GStr token;
          while (xrefs.nextToken(token)) {
            token.upper();
            if (token.startsWith("GENEID:")) {
              token.cut(0,token.index(':')+1);
              int* prevg=geneids.Find(token.chars());
              if (prevg!=NULL) (*prevg)++;
                     else geneids.Add(token, new int(1));
              }
            } //for each xref
          } //xrefs parsing
        }//for each transcript
      locus_num++;
      loc.locus_num=locus_num;
      if (gnames.Count()>0) { //collect all gene names associated to this locus
         gnames.startIterate();
         int* gfreq=NULL;
         char* key=NULL;
         while ((gfreq=gnames.NextData(key))!=NULL) {
            loc.gene_names.Add(new CGeneSym(key,*gfreq));
            }
         } //added collected gene_names
      if (loc.gene_ids.Count()>0) { //collect all GeneIDs names associated to this locus
         geneids.startIterate();
         int* gfreq=NULL;
         char* key=NULL;
         while ((gfreq=geneids.NextData(key))!=NULL) {
           loc.gene_ids.Add(new CGeneSym(key,*gfreq));
            }
          }
      } //for each locus
  }//for each genomic sequence
}


void GffLoader::load(GList<GenomicSeqData>& seqdata, GFValidateFunc* gf_validate, 
                              bool doCluster, bool doCollapseRedundant, bool collapseContained) {
   GffReader* gffr=new GffReader(f, this->transcriptsOnly, false); //not only mRNA features, not sorted
   gffr->showWarnings(this->showWarnings);
   //           keepAttrs   mergeCloseExons  noExonAttr
   gffr->readAll(this->fullAttributes,    this->mergeCloseExons,  this->noExonAttrs);
  //int redundant=0; //redundant annotation discarded
  if (verbose) GMessage("   .. loaded %d genomic features from %s\n", gffr->gflst.Count(), fname.chars());
  //int rna_deleted=0;
  //add to GenomicSeqData, adding to existing loci and identifying intron-chain duplicates
  for (int k=0;k<gffr->gflst.Count();k++) {
     GffObj* m=gffr->gflst[k];
     if (strcmp(m->getFeatureName(), "locus")==0 && 
          m->getAttr("transcripts")!=NULL) {
        continue;
        }
     
     char* rloc=m->getAttr("locus");
     if (rloc!=NULL && startsWith(rloc, "RLOC_")) {
        m->removeAttr("locus", rloc);
        }
     if (m->exons.Count()==0 && m->children.Count()==0) {
       //a non-mRNA feature with no subfeatures
       //add a dummy exon just to have the generic exon checking work
       m->addExon(m->start,m->end);
       }
     GList<GffObj> gfadd(false,false);
     if (gf_validate!=NULL && !(*gf_validate)(m, &gfadd)) {
       continue;
       }
     m->isUsed(true); //so the gffreader won't destroy it
     int i=-1;
     GenomicSeqData f(m->gseq_id);
     GenomicSeqData* gdata=NULL;
     
     if (seqdata.Found(&f,i)) gdata=seqdata[i];
         else { //entry not created yet for this genomic seq
           gdata=new GenomicSeqData(m->gseq_id);
           seqdata.Add(gdata);
           }
    for (int k=0;k<gfadd.Count();k++) {
      placeGf(gfadd[k], gdata, doCluster, doCollapseRedundant, collapseContained);
      }
    placeGf(m, gdata, doCluster, doCollapseRedundant, collapseContained);
    } //for each read gffObj
   //if (verbose) GMessage("  .. %d records from %s clustered into loci.\n", gffr->gflst.Count(), fname.chars());
   if (f!=stdin) { fclose(f); f=NULL; }
   delete gffr;
}
