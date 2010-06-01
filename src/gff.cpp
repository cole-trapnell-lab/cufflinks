#include "gff.h"

//GffNames* GffReader::names=NULL;
GffNames* GffObj::names=NULL;
//global set of feature names, attribute names etc.
// -- common for all GffObjs in current application!

const uint GFF_MAX_LOCUS = 4000000; //longest known gene in human is ~2.2M, UCSC claims a gene for mouse of ~ 3.1 M
const uint GFF_MAX_EXON  =   20000; //longest known exon in human is ~11K
const uint GFF_MAX_INTRON= 1600000;

const int gff_fid_mRNA=0;
const int gff_fid_exon=1;
const int gff_fid_CDS=2;
bool gff_warns=true;

void gffnames_ref(GffNames* &n) {
  if (n==NULL) n=new GffNames();
  n->numrefs++;
}

void gffnames_unref(GffNames* &n) {
  if (n==NULL) GError("Error: attempt to remove reference to null GffNames object!\n");
  n->numrefs--;
  if (n->numrefs==0) { delete n; n=NULL; }
}

static char fnamelc[32];

GffLine::GffLine(GffReader* reader, const char* l) {
 line=Gstrdup(l);
 skip=true;
 gseqname=NULL;
 track=NULL;
 ftype=NULL;
 info=NULL;
 Parent=NULL;
 is_cds=false;
 is_mrna=false;
 is_exon=false;
 exontype=0;
 gname=NULL;
 qstart=0;
 qend=0;
 qlen=0;
 exontype=0;
 ID=NULL;
 char* t[9];
 int i=0;
 int tidx=1;
 t[0]=line;
 
 while (line[i]!=0) {
  if (line[i]=='\t') {
   line[i]=0;
   t[tidx]=line+i+1;
   tidx++;
   if (tidx>8) break;
   }
  i++;
  }

 if (tidx<8) { // ignore non-GFF lines
  // GMessage("Warning: error parsing GFF/GTF line:\n%s\n", l);
  return;
  }
 gseqname=t[0];
 track=t[1];
 ftype=t[2];
 info=t[8];
 char* p=t[3];
 if (!parseUInt(p,fstart))
   GError("Error parsing start coordinate from GFF line:\n%s\n",l);
 p=t[4];
 if (!parseUInt(p,fend))
   GError("Error parsing end coordinate from GFF line:\n%s\n",l);
 if (fend<fstart) swap(fend,fstart);
 p=t[5];
 if (p[0]=='.' && p[1]==0) {
  score=0;
  }
 else {
  if (!parseDouble(p,score))
       GError("Error parsing feature score from GFF line:\n%s\n",l);
  }
 strand=*t[6];
 if (strand!='+' && strand!='-' && strand!='.')
     GError("Error parsing strand (%c) from GFF line:\n%s\n",strand,l);
 phase=*t[7];
 ID=NULL;
 Parent=NULL;
 // exon/CDS/mrna filter
 strncpy(fnamelc, ftype, 31);
 fnamelc[31]=0;
 strlower(fnamelc); //convert to lower case
 if (strstr(fnamelc, "utr")!=NULL) {
   exontype=exgffUTR;
   is_exon=true;
   }
  else if (strstr(fnamelc, "exon")!=NULL) {
   exontype=exgffExon;
   is_exon=true;
   }
  else if (startsWith(fnamelc, "stop") && strstr(fnamelc, "codon")!=NULL) {
   exontype=exgffStop;
   is_cds=true;
   }
 else if (strcmp(fnamelc, "cds")==0) {
   exontype=exgffCDS;
   is_cds=true;
   }
 else { //is_mrna is set only if the *current* line is a mRNA or transcript
   is_mrna=(strcmp(fnamelc,"mrna")==0 ||
          strcmp(fnamelc,"transcript")==0);
   }

 if (reader->mrnaOnly) {
   if (!is_mrna && !is_cds && !is_exon)
                     return; //skip this line, unwanted feature name
   }
 p=strstr(info,"ID=");
 if (p!=NULL) { //has ID attr
   ID=p+3;
   p=ID;
   while (*p!=';' && *p!=0) p++;
   ID=Gstrdup(ID, p-1);
   //look for a name attr too:
   p=strstr(info,"Name=");
   if (p!=NULL) {
     gname=p+5;
     p=gname;
     while (*p!=';' && *p!=0) p++;
     gname=Gstrdup(gname, p-1);
     }
   }
 p=NULL;
 if (!is_mrna)
     p=strifind(info,"Parent="); //don't care about the parent for mRNA features..
 if (p!=NULL) { //has Parent attr
   Parent=p+7;
   p=Parent;
   while (*p!=';' && *p!=0) p++;
   Parent=Gstrdup(Parent, p-1);
   }
  else if (ID==NULL) { //no "Parent=" and no "ID=", attempt GTF parsing instead
   p=strstr(info,"transcript_id");
   if (p!=NULL) { //GTF detected
     p+=13;
     //requires quotes
     while (*p!='"' && *p!=0) p++;
     if (*p==0) GError("Error parsing transcript_id (double quotes expected) at GTF line:\n%s\n",l);
     p++;
     Parent=p;
     while (*p!='"' && *p!=0) p++;
     if (*p==0) GError("Error parsing transcript_id (ending double quotes expected) at GTF line:\n%s\n",l);
     if (is_mrna) { // RGASP GTF exception: a parent "transcript" feature preceding exon/CDS subfeatures
        ID=Gstrdup(Parent, p-1);
        Parent=NULL;
        }
       else {
        Parent=Gstrdup(Parent, p-1);
        }
     //check for gene_name or gene_id
     //p=strstr(info, "gene_name");// this is preferred over gene_id
     //if (p==NULL)
     p=strstr(info,"gene_id");
     if (p!=NULL) {
       p+=7;//skip 'gene_id'
       while (*p!='"' && *p!=0) p++;
       if (*p==0) GError("Error parsing gene_id (double quotes expected) at GTF line:\n%s\n",l);
       p++;
       gname=p;
       while (*p!='"' && *p!=0) p++;
       if (*p==0) GError("Error parsing gene_id (ending double quotes expected) at GTF line:\n%s\n",l);
       gname=Gstrdup(gname, p-1);
       }
     //prepare for parseAttr by adding '=' character instead of spaces for all attributes
     //after the attribute name
     p=info;
     bool noed=true; //not edited after the last delim
     bool nsp=false; //non-space found after last delim
     while (*p!=0) {
      if (*p==' ') {
         if (nsp && noed) {
           *p='=';
            noed=false;
            p++;
            continue;
            }
         }
      else nsp=true;
      if (*p==';') { noed=true; nsp=false; }
      p++;
      }
     } //gtf detected
    else {//check for jigsaw or cufflinks format
     char* fexon=strstr(fnamelc, "exon");
     if (fexon!=NULL) {
       if (startsWith(track,"jigsaw")) {
        is_cds=true;
        strcpy(track,"jigsaw");
        p=strchr(info,';');
        if (p==NULL) Parent=Gstrdup(info);
           else { Parent=Gstrdup(info,p-1); info=p+1;  }
        }
       else if ((i=strcspn(info,"; \t\n\r"))<=(int)(strlen(info)+1)) {//one word ID
          Parent=Gstrdup(info,info+i-1);
        }

      }
      else GError("Error parsing Parent/ID at input line:\n%s\n",l);
     }
   }
 p=strstr(info,"Target=");
 if (p!=NULL) { //has Target attr
   p+=7;
   while (*p!=';' && *p!=0 && *p!=' ') p++;
   if (*p!=' ') {
      GError("Error parsing target coordinates from GFF line:\n%s\n",l);
      }
   if (!parseUInt(p,qstart))
     GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
   if (*p!=' ') {
      GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
      }
   p++;
   if (!parseUInt(p,qend))
     GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
   }
 else {
   p=strifind(info,"Qreg=");
   if (p!=NULL) { //has Qreg attr
     p+=5;
     if (!parseUInt(p,qstart))
       GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
     if (*p!='-') {
        GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
        }
     p++;
     if (!parseUInt(p,qend))
       GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
     if (*p=='|') {
       p++;
       if (!parseUInt(p,qlen))
         GError("Error parsing target length from GFF Qreg|: \n%s\n",l);
       }
     }//has Qreg attr
   }
 if (qlen==0 && (p=strifind(info,"Qlen="))!=NULL) {
   p+=5;
   if (!parseUInt(p,qlen))
       GError("Error parsing target length from GFF Qlen:\n%s\n",l);
   }
 skip=false;
}

int GffObj::addExon(GffLine* gl, bool keepAttr) {
  //this will make sure we have the right subftype_id!
  int subf_id=-1;
  if (ftype_id==gff_fid_mRNA) {
     if (subftype_id<0) subftype_id=gff_fid_exon; 
     //any recognized mRNA segment gets the generic "exon" type (also applies to CDS)
     if (!gl->is_cds && !gl->is_exon)
         //extraneous mRNA feature, will discard for now
          return -1;
     }
  else { //other kind of parent feature, check this subf type
    subf_id=names->feats.addName(gl->ftype);
    if (subftype_id<0)
       subftype_id=subf_id;
    else {
       if (subftype_id!=subf_id)
         GMessage("GFF Warning: multiple subfeatures (%s and %s) found for %s, only %s is kept\n",
             names->feats.getName(subftype_id), names->feats.getName(subf_id),
             gffID,names->feats.getName(subftype_id));
         return -1; //skip this 2nd subfeature type for this parent!
       }
    }
  if (gl->exontype==exgffUTR || gl->exontype==exgffStop) 
      udata=1;//merge 0-distance segments
  int eidx=addExon(gl->fstart, gl->fend, gl->score, gl->phase,
         gl->qstart,gl->qend, gl->is_cds, gl->exontype);
  if (eidx>=0 && keepAttr)
      parseAttrs(exons[eidx]->attrs, gl->info);
  return eidx;
}


void GffObj::expandExon(int oi, uint segstart, uint segend, char exontype, double sc, char fr, int qs, int qe) {
  //oi is the index of the *first* overlapping segment found that must be enlarged
  covlen-=exons[oi]->len();
  if (segstart<exons[oi]->start)
    exons[oi]->start=segstart;
  if (qs<exons[oi]->qstart) exons[oi]->qstart=qs;
  if (segend>exons[oi]->end)
    exons[oi]->end=segend;
  if (qe>exons[oi]->qend) exons[oi]->qend=qe;
  //warning: score cannot be properly adjusted! e.g. if it's a p-value it's just going to get worse
  if (sc!=0) exons[oi]->score=sc;
  covlen+=exons[oi]->len();
  //if (exons[oi]->exontype< exontype) -- always true
  exons[oi]->exontype = exontype;
  if (exontype==exgffCDS) exons[oi]->phase=fr; 
  //we must check if any more exons are also overlapping this 
  int ni=oi+1; //next exon index after oi
  while (ni<exons.Count() && segend>=exons[ni]->start) { // next segment overlaps new enlarged segment
     //only allow this if next segment is fully included, and a subordinate
     if (exons[ni]->exontype<exontype && exons[ni]->end<=segend) {
         if (exons[ni]->qstart<exons[oi]->qstart) exons[oi]->qstart=exons[ni]->qstart;
         if (exons[ni]->qend>exons[oi]->qend) exons[oi]->qend=exons[ni]->qend;
         exons.Delete(ni);
         }
      else {
         GMessage("GFF Warning: overlapping feature segment (%d-%d) for GFF ID %s\n", segstart, segend, gffID);
         hasErrors=true;
         break;
         }
     }
  // -- make sure any other related boundaries are updated:
  start=exons.First()->start;
  end=exons.Last()->end;
  if (uptr!=NULL) { //collect stats about the underlying genomic sequence
    GSeqStat* gsd=(GSeqStat*)uptr;
    if (start<gsd->mincoord) gsd->mincoord=start;
    if (end>gsd->maxcoord) gsd->maxcoord=end;
    }
}

int GffObj::addExon(uint segstart, uint segend, double sc, char fr, int qs, int qe, bool iscds, char exontype) {
  if (exons.Count()==0) {
      if (iscds) isCDS=true; //for now, assume CDS only if first "exon" given is a CDS
      if (subftype_id<0) {
         subftype_id = (ftype_id==gff_fid_mRNA) ? gff_fid_exon : ftype_id;
         }
      }
  if (iscds) { //update CDS anchors:
     if (CDstart==0 || segstart<CDstart)  {
           CDstart=segstart;
           if (exontype==exgffCDS && strand=='+') CDphase=fr;
           }
     if (segend>CDend) {
           if (exontype==exgffCDS && strand=='-') CDphase=fr;
           CDend=segend;
           }
     }
   else { // !iscds
     isCDS=false;
     }
  if (qs || qe) {
    if (qs>qe) swap(qs,qe);
    if (qs==0) qs=1;
    }
  int ovlen=0;
  int oi=exonOverlapIdx(segstart, segend, &ovlen);
  if (oi>=0) { //overlap existing segment
     //only allow this for CDS within exon, stop_codon within exon, stop_codon within UTR,
     //                        or stop_codon within CDS
    if (exons[oi]->exontype>exontype && 
         exons[oi]->start<=segstart && exons[oi]->end>=segend &&
         !(exons[oi]->exontype==exgffUTR && exontype==exgffCDS)) {
          //larger segment given first, now the smaller included one
          return oi; //only used to store attributes from current GffLine
          }
    if (exontype>exons[oi]->exontype && 
         segstart<=exons[oi]->start && segend>=exons[oi]->end &&
         !(exontype==exgffUTR && exons[oi]->exontype==exgffCDS)) {
           //smaller segment given first, so we have to enlarge it
          expandExon(oi, segstart, segend, exontype, sc, fr, qs, qe); 
            //this must also check for overlapping next exon (oi+1) 
          return oi;
          }
    //there is also the special case of "ribosomal slippage exception" (programmed frameshift)
    //where two CDS segments may actually overlap for 1 or 2 bases, but there should be only one encompassing exon
    if (ovlen>2 || exons[oi]->exontype!=exgffCDS || exontype!=exgffCDS) {
       GMessage("GFF Warning: discarded overlapping feature segment (%d-%d) for GFF ID %s\n", 
           segstart, segend, gffID);
       hasErrors=true;
       return false; //segment NOT added
       }
     /* -- 
     else {
       // else add the segment if the overlap is small and between two CDS segments
       //we might want to add an attribute here with the slippage coordinate and size?
       } 
     */
    }//overlap of existing segment

   // --- no overlap, or accepted micro-overlap (ribosomal slippage)
   // create & add the new segment
   GffExon* enew=new GffExon(segstart, segend, sc, fr, qs, qe, exontype);
   int eidx=exons.Add(enew);
   if (eidx<0) {
     GMessage("GFF Warning: duplicate segment %d-%d not added for GFF ID %s!\n",
            segstart, segend, gffID);
     return -1;            
     }
   covlen+=(int)(exons[eidx]->end-exons[eidx]->start)+1;
   start=exons.First()->start;
   end=exons.Last()->end;
   if (uptr!=NULL) { //collect stats about the underlying genomic sequence
       GSeqStat* gsd=(GSeqStat*)uptr;
       if (start<gsd->mincoord) gsd->mincoord=start;
       if (end>gsd->maxcoord) gsd->maxcoord=end;
       }
   return eidx;
}

void GffObj::removeExon(int idx) {
  /*
   if (idx==0 && segs[0].start==gstart)
                  gstart=segs[1].start;
   if (idx==segcount && segs[segcount].end==gend)
                  gend=segs[segcount-1].end;
  */
  if (idx<0 || idx>=exons.Count()) return;
  int segstart=exons[idx]->start;
  int segend=exons[idx]->end;
  exons.Delete(idx);
  covlen -= (int)(segend-segstart)+1;
  start=exons.First()->start;
  end=exons.Last()->end;
  if (isCDS) { CDstart=start; CDend=end; }
}

GffObj::GffObj(GffReader *gfrd, GffLine* gffline, bool keepAttr, bool noExonAttr):
     GSeg(0,0), exons(true,true,true) {
 xstart=0;
 xend=0;
 xstatus=0;
 partial=false;
 isCDS=false;
 uptr=NULL;
 ulink=NULL;
 udata=0;
 CDstart=0;
 CDend=0;
 gname=NULL;
 attrs=NULL;
 gffID=NULL;
 track_id=-1;
 gseq_id=-1;
 ftype_id=-1;
 subftype_id=-1;
 hasErrors=false;
 if (gfrd==NULL)
    GError("Cannot use this GffObj constructor with a NULL GffReader!\n");
 gffnames_ref(names);
 if (gfrd->names==NULL) gfrd->names=names;
 qlen=0;qstart=0;qend=0;
 gscore=0;
 uscore=0;
 covlen=0;
 qcov=0;
 if (gffline->Parent!=NULL) {
    //GTF style -- subfeature given directly
    if (gffline->is_cds || gffline->is_exon)
         ftype_id=gff_fid_mRNA;
      else {
        //group of other subfeatures of type ftype:
        ftype_id=names->feats.addName(gffline->ftype);
        }
    gffID=gffline->Parent;
    gffline->Parent=NULL; //just take over
    if (gffline->gname!=NULL) {
        gname=gffline->gname;
        gffline->gname=NULL;
        }
    gseq_id=names->gseqs.addName(gffline->gseqname);
    track_id=names->tracks.addName(gffline->track);
    strand=gffline->strand;
    qlen=gffline->qlen;
    start=gffline->fstart;
    end=gffline->fend;
    isCDS=gffline->is_cds; //for now
    addExon(gffline, keepAttr);
    if (keepAttr && noExonAttr) {
      //simply move the attrs from this first exon
      //to the transcript
      attrs=exons.First()->attrs;
      exons.First()->attrs=NULL;
      }
    }
 else { //GffReader made sure this is a parent line (no parent)
    //even for a mRNA with a Parent= line
    gscore=gffline->score;
    if (gffline->ID==NULL || gffline->ID[0]==0)
       GError("Error: no ID found for GFF record start\n");
    gffID=gffline->ID; //there must be an ID here
    if (gffline->is_mrna) ftype_id=gff_fid_mRNA;
        else ftype_id=names->feats.addName(gffline->ftype);
    gffline->ID=NULL; //steal it
    if (gffline->gname!=NULL) {
        gname=gffline->gname;
        gffline->gname=NULL;
        }
    start=gffline->fstart;
    end=gffline->fend;
    gseq_id=names->gseqs.addName(gffline->gseqname);
    track_id=names->tracks.addName(gffline->track);
    qlen=gffline->qlen;
    qstart=gffline->qstart;
    qend=gffline->qend;
    strand=gffline->strand;
    if (keepAttr) this->parseAttrs(attrs, gffline->info);
    }
 GSeqStat* gsd=gfrd->gseqstats.AddIfNew(new GSeqStat(gseq_id,names->gseqs.lastNameUsed()),true);
 uptr=gsd;
 gsd->gflst.Add(this);
 if (start<gsd->mincoord) gsd->mincoord=start;
 if (end>gsd->maxcoord) gsd->maxcoord=end;
 gfrd->gfoAdd(gffID, gffline->gseqname, this);
}


GffLine* GffReader::nextGffLine() {
 if (gffline!=NULL) return gffline; //caller should free gffline after processing
 while (gffline==NULL) {
    //const char* l=linebuf->getLine();
    int llen=0;
    buflen=GFF_LINELEN-1;
    char* l=fgetline(linebuf, buflen, fh, &fpos, &llen);
    if (l==NULL) {
         return NULL; //end of file
         }
    int ns=0; //first nonspace position
    while (l[ns]!=0 && isspace(l[ns])) ns++;
    if (l[ns]=='#' || llen<10) continue;
    gffline=new GffLine(this, l);
    if (gffline->skip) {
       delete gffline;
       gffline=NULL;
       }
    }
return gffline;
}


GffObj* GffReader::parse(bool keepAttr, bool noExonAttr) { //<- do not use this!
   //parses one parent feature at a time, assuming absolute grouping of subfeatures
   //!ASSUMES that subfeatures (exons) are properly grouped together
   //any interleaving exons from other transcripts will terminate the parsing of the current feature!
  GffObj* gfo=NULL;
  while (nextGffLine()!=NULL) {
    if (gfo==NULL) {//record starts fresh here
       gfo=new GffObj(this, gffline, keepAttr, noExonAttr);
       delete gffline; gffline=NULL;
       continue;
       }
 // -- gfo is not NULL from here --
    if (gffline->Parent==NULL) {// new parent feature starting here
       //new record start, return what we have so far,
       //gffline was NOT deleted, so it will be used for the next parse() call
       return ((gfo==NULL) ? NULL : gfo->finalize());
       }
    //-- has a Parent so it's a subfeature segment (exon/CDS/other subfeature)
      // is it a subfeature of the current gf?
    if (strcmp(gffline->Parent, gfo->gffID)==0) {
          //yes, add it
          gfo->addExon(gffline, !noExonAttr);
          delete gffline; gffline=NULL;
          continue;
          }
      // is it a subfeature of a previously loaded gfo?
    GffObj* prevgfo=gfoFind(gffline->Parent, gffline->gseqname);
    if (prevgfo==NULL)
           return ((gfo==NULL) ? NULL : gfo->finalize()); // new subfeature, gffline will be used for the next parse()
    //this is for an earlier parent
    prevgfo->addExon(gffline, !noExonAttr);
    delete gffline;
    gffline=NULL;
    } //while reading gfflines
  return ((gfo==NULL) ? NULL : gfo->finalize());
}
char* GffReader::gfoBuildId(const char* id, const char* ctg) {
//caller must free the returned pointer
 char* buf=NULL;
 int idlen=strlen(id);
 GMALLOC(buf, idlen+strlen(ctg)+2);
 strcpy(buf, id);
 buf[idlen]='~';
 strcpy(buf+idlen+1, ctg);
 return buf;
}

void GffReader::gfoRemove(const char* id, const char* ctg) {
 char* buf=gfoBuildId(id,ctg);
 phash.Remove(buf);
 GFREE(buf);
}

void GffReader::gfoAdd(const char* id, const char* ctg, GffObj* gfo) {
 char* buf=gfoBuildId(id,ctg);
 phash.Add(buf, gfo);
 GFREE(buf);
}
GffObj* GffReader::gfoFind(const char* id, const char* ctg) {
 char* buf=gfoBuildId(id,ctg);
 GffObj* r=phash.Find(buf);
 GFREE(buf);
 return r;
}


void GffReader::parseAll(GffRecFunc* gproc, bool keepAttr, bool noExonAttr, void* userptr1, void* userptr2) {
  //WARNING: this is all messed up if the GFF lines are NOT grouped by parent feature
  //!!! do not use this !!!
  //iterates through all mappings in the input file
  //calling gproc with each parsed mapping
    GffObj* gfo;
    while ((gfo=this->parse(keepAttr,noExonAttr))!=NULL) { //a valid gff record was parsed
        if (gfo->empty()) { //shouldn't happen!
             delete gfo;
             gfo=NULL;
             continue;
             }
     if (gproc(gfo, userptr1, userptr2)) {
             //true returned from GfProcFunc means no longer needed
              gfoRemove(gfo->gffID, gfo->getGSeqName());
              delete gfo;
              }
         else {
              gflst.Add(gfo);
              }
     gfo=NULL;
     } //while records are parsed
    phash.Clear();
}


void GffReader::readAll(bool keepAttr, bool mergeCloseExons, bool noExonAttr) {
  while (nextGffLine()!=NULL) {
    if (gffline->Parent==NULL) {//no parent, new GFF3-like record starting
       //check for uniqueness of gffline->ID in phash !
       GffObj* f=gfoFind(gffline->ID, gffline->gseqname);
       if (f!=NULL) {
            GError("Error: duplicate GFF ID '%s' encountered!\n",gffline->ID);
            }
       gflst.Add(new GffObj(this, gffline, keepAttr, noExonAttr));
       }
    else { //--- it's a subfeature (exon/CDS/other):
       GffObj* prevgfo=gfoFind(gffline->Parent, gffline->gseqname);
       if (prevgfo!=NULL) { //exon of a previously seen GffObj
                 if (gffline->strand!=prevgfo->strand) {
                    GError("Error: duplicate GFF ID '%s' (exons found on different strands of %s)\n",
                       prevgfo->gffID, prevgfo->getGSeqName());
                    }
                 int gdist=(gffline->fstart>prevgfo->end) ? gffline->fstart-prevgfo->end :
                                     ((gffline->fend<prevgfo->start)? prevgfo->start-gffline->fend : 
                                        0 );
                 if (gdist>(int)GFF_MAX_LOCUS) { //too far apart, most likely this is a duplicate ID
                   GError("Error: duplicate GFF ID '%s' (or exons too far apart)!\n",prevgfo->gffID);
                   }
                 prevgfo->addExon(gffline, !noExonAttr);
                 }
            else {//new GTF-like record starting here with a subfeature
                 gflst.Add(new GffObj(this, gffline, keepAttr, noExonAttr));
                 //even those with errors will be added here!
                 }
       } //subfeature
      //--
    delete gffline;
    gffline=NULL;
    }//while
  for (int i=0;i<gflst.Count();i++) {
    gflst[i]->finalize(mergeCloseExons); //finalize the parsing - also merge close exons/features if so requested
    }
 // all gff records are now loaded in GList gflst
 // so we can free the hash
  phash.Clear();
}

//this may be called prematurely if exon records are not grouped by parent
GffObj* GffObj::finalize(bool mergeCloseExons) {
 udata=0;
 uptr=NULL;
 //TODO: merge adjacent or close segments - for mRNAs
 //must merge adjacent features, and optionally small gaps
 if (ftype_id==gff_fid_mRNA) {
   int mindist=mergeCloseExons ? 5:1;
   for (int i=0;i<exons.Count()-1;i++) {
     int ni=i+1;
     uint mend=exons[i]->end;
     while (ni<exons.Count()) {
       int dist=(int)(exons[ni]->start-mend);
       if (dist<0 || dist>mindist) break; //no merging with next segment
       mend=exons[ni]->end;
       exons[i]->end=mend;
       if (exons[ni]->attrs!=NULL && (exons[i]->attrs==NULL || 
            exons[i]->attrs->Count()<exons[ni]->attrs->Count())) {
              //use the other exon attributes, if more
              delete(exons[i]->attrs);
              exons[i]->attrs=exons[ni]->attrs;
              exons[ni]->attrs=NULL;
              }
       exons.Delete(ni);
       } //check for merge with next exon
     } //for each exon
   }
 return this;
}

void GffObj::parseAttrs(GffAttrs*& atrlist, char* info) {
  if (names==NULL)
     GError(ERR_NULL_GFNAMES, "parseAttrs()");
  if (atrlist==NULL)
      atrlist=new GffAttrs();
  char* endinfo=info+strlen(info);
  char* start=info;
  char* pch=start;
  while (start<endinfo) {
    //skip spaces
    while (*start==' ' && start<endinfo) start++;
    pch=strchr(start, ';');
    if (pch==NULL) pch=endinfo;
       else {
            *pch='\0';
            pch++;
            }
    char* ech=strchr(start,'=');
    if (ech!=NULL) { // attr=value format found
       *ech='\0';
       /*if (strcmp(start, "ID")==0 || strcmp(start,"Target")==0 || Gstricmp(start, "Qreg")==0 ||
           Gstricmp(start, "Qlen")==0 || strcmp(start,"Parent")==0 || strcmp(start,"Name")==0 ||
           strcmp(start,"transcript_id")==0 || strcmp(start,"gene_id")==0) */
       if (strcmp(start, "ID")==0 || strcmp(start,"Parent")==0 || strcmp(start,"Name")==0 ||
            strcmp(start,"transcript_id")==0 || strcmp(start,"gene_id")==0)
          { start=pch; continue; } //skip this already recognized and stored attribute
       ech++;
       while (*ech==' ' && ech<endinfo) ech++;//skip extra spaces after the '='
       atrlist->Add(new GffAttr(names->attrs.addName(start),ech));
       }
      /*
      else { //not an attr=value format
        atrlist->Add(new GffAttr(names->attrs.addName(start),"1"));
        }
      */
    start=pch;
    }
  if (atrlist->Count()==0) { delete atrlist; atrlist=NULL; }
}

void GffObj::addAttr(const char* attrname, char* attrvalue) {
  if (this->attrs==NULL)
      this->attrs=new GffAttrs();
  this->attrs->Add(new GffAttr(names->attrs.addName(attrname),attrvalue));
}

void GffObj::getCDS_ends(uint& cds_start, uint& cds_end) {
  cds_start=0;
  cds_end=0;
  if (CDstart==0 || CDend==0) return; //no CDS info
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  cds_start=CDstart;
  cds_end=CDend;
  if (strand=='-') cds_end-=cdsadj;
              else cds_start+=cdsadj;
  }

void GffObj::mRNA_CDS_coords(uint& cds_mstart, uint& cds_mend) {
  //sets cds_start and cds_end to the CDS start,end coordinates on the spliced mRNA transcript
  cds_mstart=0;
  cds_mend=0;
  if (CDstart==0 || CDend==0) return; //no CDS info
  //restore normal coordinates, just in case
  unxcoord();
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  /*
   uint seqstart=CDstart;
   uint seqend=CDend;
  */
  uint seqstart=exons.First()->start;
  uint seqend=exons.Last()->end;
  int s=0; //resulting nucleotide counter
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend)
             sgend=seqend; //seqend within this segment
       s+=(int)(sgend-sgstart)+1;
       if (CDstart>=sgstart && CDstart<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             cds_mend=s-(int)(CDstart-sgstart);
             //GMessage("Setting cds_mend to %d\n",cds_mend);
             }
       if (CDend>=sgstart && CDend<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             cds_mstart=s-(int)(CDend-cdsadj-sgstart);
             //GMessage("Setting cds_mstart to %d\n",cds_mstart);
             }
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      s+=(int)(sgend-sgstart)+1;
      /* for (uint i=sgstart;i<=sgend;i++) {
          spliced[s]=gsubseq[i-gstart];
          s++;
          }//for each nt
          */
      if (CDstart>=sgstart && CDstart<=sgend) {
            //CDstart in this segment
            cds_mstart=s-(int)(sgend-CDstart-cdsadj);
            }
      if (CDend>=sgstart && CDend<=sgend) {
            //CDend in this segment
            cds_mend=s-(int)(sgend-CDend);
            }
      } //for each exon
    } // + strand
  //spliced[s]=0;
  //if (rlen!=NULL) *rlen=s;
  //return spliced;
}

char* GffObj::getSpliced(GFaSeqGet* faseq, bool CDSonly, int* rlen, uint* cds_start, uint* cds_end,
          GList<GSeg>* seglst) {
  if (CDSonly && CDstart==0) return NULL;
  if (faseq==NULL) { GMessage("Warning: getSpliced(NULL,.. ) called!\n");
              return NULL;
              }
  //restore normal coordinates:
  unxcoord();
  if (exons.Count()==0) return NULL;
  int fspan=end-start+1;
  const char* gsubseq=faseq->subseq(start, fspan);
  if (gsubseq==NULL) {
        GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
        }
  char* spliced=NULL;
  GMALLOC(spliced, covlen+1); //allocate more here
  uint seqstart, seqend;
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  if (CDSonly) {
     seqstart=CDstart;
     seqend=CDend;
     if (strand=='-') seqend-=cdsadj;
           else seqstart+=cdsadj;
     }
   else {
     seqstart=exons.First()->start;
     seqend=exons.Last()->end;
     }
  int s=0; //resulting nucleotide counter
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend)
             sgend=seqend; //seqend within this segment
       if (seglst!=NULL)
           seglst->Add(new GSeg(s+1,s+1+sgend-sgstart));
       for (uint i=sgend;i>=sgstart;i--) {
            spliced[s] = ntComplement(gsubseq[i-start]);
            s++;
            }//for each nt

       if (!CDSonly && cds_start!=NULL && CDstart>0) {
          if (CDstart>=sgstart && CDstart<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             *cds_end=s-(CDstart-sgstart);
             }
          if (CDend>=sgstart && CDend<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             *cds_start=s-(CDend-cdsadj-sgstart);
             }
         }//update local CDS coordinates
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      if (seglst!=NULL)
          seglst->Add(new GSeg(s+1,s+1+sgend-sgstart));
      for (uint i=sgstart;i<=sgend;i++) {
          spliced[s]=gsubseq[i-start];
          s++;
          }//for each nt
      if (!CDSonly && cds_start!=NULL && CDstart>0) {
         if (CDstart>=sgstart && CDstart<=sgend) {
            //CDstart in this segment
            //and we are getting the whole transcript
            *cds_start=s-(sgend-CDstart-cdsadj);
            }
         if (CDend>=sgstart && CDend<=sgend) {
            //CDstart in this segment
            //and we are getting the whole transcript
            *cds_end=s-(sgend-CDend);
            }
        }//update local CDS coordinates
      } //for each exon
    } // + strand
  spliced[s]=0;
  if (rlen!=NULL) *rlen=s;
  return spliced;
}

char* GffObj::getSplicedTr(GFaSeqGet* faseq, bool CDSonly, int* rlen) {
  if (CDSonly && CDstart==0) return NULL;
  //restore normal coordinates:
  unxcoord();
  if (exons.Count()==0) return NULL;
  int fspan=end-start+1;
  const char* gsubseq=faseq->subseq(start, fspan);
  if (gsubseq==NULL) {
    GError("Error getting subseq for %s (%d..%d)!\n", gffID, start, end);
    }

  char* translation=NULL;
  GMALLOC(translation, (int)(covlen/3)+1);
  uint seqstart, seqend;
  int cdsadj=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsadj=CDphase-'0';
      }
  if (CDSonly) {
     seqstart=CDstart;
     seqend=CDend;
     if (strand=='-') seqend-=cdsadj;
           else seqstart+=cdsadj;
     }
   else {
     seqstart=exons.First()->start;
     seqend=exons.Last()->end;
     }
  Codon codon;
  int nt=0; //codon nucleotide counter (0..2)
  int aa=0; //aminoacid count
  if (strand=='-') {
    for (int x=exons.Count()-1;x>=0;x--) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (seqend<sgstart || seqstart>sgend) continue;
       if (seqstart>=sgstart && seqstart<=sgend)
             sgstart=seqstart; //seqstart within this segment
       if (seqend>=sgstart && seqend<=sgend) {
             sgend=seqend; //seqend within this segment
             }
       for (uint i=sgend;i>=sgstart;i--) {
            codon.nuc[nt]=ntComplement(gsubseq[i-start]);
            nt++;
            if (nt==3) {
               nt=0;
               translation[aa]=codon.translate();
               aa++;
               }
            }//for each nt
      } //for each exon
    } // - strand
   else { // + strand
    for (int x=0;x<exons.Count();x++) {
      uint sgstart=exons[x]->start;
      uint sgend=exons[x]->end;
      if (seqend<sgstart || seqstart>sgend) continue;
      if (seqstart>=sgstart && seqstart<=sgend)
            sgstart=seqstart; //seqstart within this segment
      if (seqend>=sgstart && seqend<=sgend)
            sgend=seqend; //seqend within this segment
      for (uint i=sgstart;i<=sgend;i++) {
          codon.nuc[nt]=gsubseq[i-start];
          nt++;
          if (nt==3) {
             nt=0;
             translation[aa]=codon.translate();
             aa++;
             }
          }//for each nt
        } //for each exon
    } // + strand
 translation[aa]=0;
 if (rlen!=NULL) *rlen=aa;
 return translation;
}

void GffObj::printSummary(FILE* fout) {
 if (fout==NULL) fout=stdout;
 fprintf(fout, "%s\t%c\t%d\t%d\t%4.2f\t%4.1f\n", gffID,
          strand, start, end, gscore, (float)qcov/10.0);
}

void GffObj::printGxfLine(FILE* fout, char* tlabel, char* gseqname, bool iscds,
                             uint segstart, uint segend, int exidx, char phase, bool gff3) {
  static char scorestr[14];
  strcpy(scorestr,".");
  GffAttrs* xattrs=NULL;
  if (exidx>=0) {
     if (exons[exidx]->score) sprintf(scorestr,"%.2f", exons[exidx]->score);
     xattrs=exons[exidx]->attrs;
  }
  char* geneid=(gname!=NULL)? gname : gffID;
  if (phase==0 || !iscds) phase='.';
  const char* ftype=iscds ? "CDS" : getSubfName();
  if (gff3) {
    fprintf(fout,
      "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\tParent=%s",
      gseqname, tlabel, ftype, segstart, segend, scorestr, strand,
      phase, gffID);
    if (xattrs!=NULL) {
      for (int i=0;i<xattrs->Count();i++)
         fprintf(fout, ";%s=%s",names->attrs.getName(xattrs->Get(i)->attr_id),
                           xattrs->Get(i)->attr_val);
         }
    fprintf(fout, "\n");
    } //GFF
  else {//for GTF -- we can only print mRNAs here
    fprintf(fout, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\t",
        gseqname, tlabel, ftype, segstart, segend, scorestr, strand, phase);
    if (ismRNA())
       fprintf(fout,"gene_id \"%s\"; transcript_id \"%s\";", geneid, gffID);
    if (xattrs!=NULL) {
       for (int i=0;i<xattrs->Count();i++) {
         if (xattrs->Get(i)->attr_val==NULL) continue;
         fprintf(fout, " %s ",names->attrs.getName(xattrs->Get(i)->attr_id));
          if (xattrs->Get(i)->attr_val[0]=='"')
                  fprintf(fout, "%s;",xattrs->Get(i)->attr_val);
             else fprintf(fout, "\"%s\";",xattrs->Get(i)->attr_val);
          }
       }
    fprintf(fout, "\n");
    }//GTF
}

void GffObj::printGxf(FILE* fout, GffPrintMode gffp, char* tlabel) {
 static char tmpstr[255];
 if (tlabel==NULL) {
    tlabel=track_id>=0 ? names->tracks.Get(track_id)->name :
         (char*)"gffobj" ;
    }

 unxcoord();
 if (exons.Count()==0) return;
 char* gseqname=names->gseqs.Get(gseq_id)->name;
 bool gff3 = (gffp>=pgffAny);
 bool showCDS = (gffp==pgtfAny || gffp==pgtfCDS || gffp==pgffCDS || gffp==pgffAny || gffp==pgffBoth);
 bool showExon = (gffp<=pgtfExon || gffp==pgffAny || gffp==pgffExon || gffp==pgffBoth);
 if (gff3) {
   //print GFF3 mRNA line:
   if (gscore>0.0) sprintf(tmpstr,"%.2f", gscore);
          else strcpy(tmpstr,".");
   uint pstart, pend;
   if (gffp==pgffCDS) {
      pstart=CDstart;
      pend=CDend;
      }
   else { pstart=start;pend=end; }
   const char* ftype=ismRNA() ? "mRNA" : getFeatureName();
   fprintf(fout,
     "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\tID=%s",
     gseqname, tlabel, ftype, pstart, pend, tmpstr, strand, gffID);
   if (gname!=NULL)
       fprintf(fout, ";Name=%s",gname);
   if (CDstart>0 && !showCDS && !isCDS) fprintf(fout,";CDS=%d:%d",CDstart,CDend);
   if (attrs!=NULL) {
      for (int i=0;i<attrs->Count();i++) {
        fprintf(fout,";%s=%s", names->attrs.getName(attrs->Get(i)->attr_id),
               attrs->Get(i)->attr_val);
        }
      }
    fprintf(fout,"\n");
   }// gff3 mRNA line
 if (showExon) {
   //print exons
   for (int i=0;i<exons.Count();i++) {
     printGxfLine(fout, tlabel, gseqname, isCDS, exons[i]->start, exons[i]->end, i, exons[i]->phase, gff3);
     }
 }//printing exons
 if (showCDS && !isCDS && CDstart>0) {
    GArray<GffCDSeg> cds(true,true);
    getCDSegs(cds);
    for (int i=0;i<cds.Count();i++) {
      printGxfLine(fout, tlabel, gseqname, true, cds[i].start, cds[i].end, -1, cds[i].phase, gff3);
      }
  } //showCDS
}


void GffObj::getCDSegs(GArray<GffCDSeg>& cds) {
  GffCDSeg cdseg;
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
      }
  if (strand=='-') {
     for (int x=exons.Count()-1;x>=0;x--) {
        uint sgstart=exons[x]->start;
        uint sgend=exons[x]->end;
        if (CDend<sgstart || CDstart>sgend) continue;
        if (CDstart>=sgstart && CDstart<=sgend)
              sgstart=CDstart; //cdstart within this segment
        if (CDend>=sgstart && CDend<=sgend)
              sgend=CDend; //cdend within this segment
        cdseg.start=sgstart;
        cdseg.end=sgend;
        cdseg.exonidx=x;
        //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
        cdseg.phase='0'+ (3-cdsacc%3)%3;
        cdsacc+=sgend-sgstart+1;
        cds.Add(cdseg);
       } //for each exon
     } // - strand
    else { // + strand
     for (int x=0;x<exons.Count();x++) {
       uint sgstart=exons[x]->start;
       uint sgend=exons[x]->end;
       if (CDend<sgstart || CDstart>sgend) continue;
       if (CDstart>=sgstart && CDstart<=sgend)
             sgstart=CDstart; //seqstart within this segment
       if (CDend>=sgstart && CDend<=sgend)
             sgend=CDend; //seqend within this segment
       cdseg.start=sgstart;
       cdseg.end=sgend;
       cdseg.exonidx=x;
       //cdseg.phase='0'+(cdsacc>0 ? (3-cdsacc%3)%3 : 0);
       cdseg.phase='0' + (3-cdsacc%3)%3 ;
       cdsacc+=sgend-sgstart+1;
       cds.Add(cdseg);
       } //for each exon
   } // + strand
}
/*
#ifdef DEBUG
void GffObj::dbgPrint(const char* msg) {
 if (msg!=NULL) fprintf(stdout, ">> %s\n",msg);
 char* tlabel=track_id>=0 ? names->tracks.Get(track_id)->name :
       (char*)"gmapobj" ;
 char scorestr[14];
 char strand=revstrand?'-':'+';
 unxcoord();
 char* gseqname=names->gseqs.Get(gseq_id)->name;
 char* fname=f_id>=0 ? names->feats.Get(f_id)->name : (char*)"nofeatname";

 fprintf(stdout, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tID=%s;Name=%s\n",
       gseqname, tlabel, fname, start, end, strand, gffID, gffID);

 for (int fi=0;fi<features->Count();fi++) {
   GFeature* feature=features->Get(fi);
   fname=names->feats.Get(feature->name_id)->name;
   GffExon* segs=feature->segs;
   int segcount=feature->segcount;
   if (segcount==0 || segs==NULL) continue;
   for (int i=0;i<segcount;i++) {
      if (segs[i].start==0) continue;
      if (segs[i].score) sprintf(scorestr,"%.2f", segs[i].score/100.00);
                  else strcpy(scorestr,".");
      fprintf(stdout,
         "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\tParent=%s\n",
         gseqname, tlabel, fname, segs[i].start, segs[i].end, scorestr, strand, gffID);
      }
   }
 fflush(stdout);
}
#endif
*/

