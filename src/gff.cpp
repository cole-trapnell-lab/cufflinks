#include "gff.h"

//GffNames* GffReader::names=NULL;
GffNames* GffObj::names=NULL;
//global set of feature names, attribute names etc.
// -- common for all GffObjs in current application!

const uint GFF_MAX_LOCUS = 7000000; //longest known gene in human is ~2.2M, UCSC claims a gene for mouse of ~ 3.1 M
const uint GFF_MAX_EXON  =   30000; //longest known exon in human is ~11K
const uint GFF_MAX_INTRON= 6000000; //Ensembl shows a >5MB human intron 
bool gff_show_warnings = false; //global setting, set by GffReader->showWarnings()
const int gff_fid_mRNA=0;
const int gff_fid_transcript=1;
const int gff_fid_exon=2;
const int gff_fid_CDS=3; //never really used in GffObj ftype_id or subftype_id
const uint gfo_flag_HAS_ERRORS       = 0x00000001;
const uint gfo_flag_CHILDREN_PROMOTED= 0x00000002;
const uint gfo_flag_IS_GENE          = 0x00000004;
const uint gfo_flag_IS_TRANSCRIPT    = 0x00000008;
const uint gfo_flag_FROM_GFF3        = 0x00000010;
const uint gfo_flag_BY_EXON          = 0x00000020; //created by subfeature (exon) directly
const uint gfo_flag_DISCARDED        = 0x00000100;
const uint gfo_flag_LST_KEEP         = 0x00000200;
const uint gfo_flag_LEVEL_MSK        = 0x00FF0000;
const byte gfo_flagShift_LEVEL           = 16;

void gffnames_ref(GffNames* &n) {
  if (n==NULL) n=new GffNames();
  n->numrefs++;
}

void gffnames_unref(GffNames* &n) {
  if (n==NULL) GError("Error: attempt to remove reference to null GffNames object!\n");
  n->numrefs--;
  if (n->numrefs==0) { delete n; n=NULL; }
}

int gfo_cmpByLoc(const pointer p1, const pointer p2) {
 
 GffObj& g1=*((GffObj*)p1);
 GffObj& g2=*((GffObj*)p2);
 if (g1.gseq_id==g2.gseq_id) {
             if (g1.start!=g2.start)
                    return (int)(g1.start-g2.start);
               else if (g1.getLevel()!=g2.getLevel())
                        return (int)(g1.getLevel()-g2.getLevel());
                    else
                        if (g1.end!=g2.end)
                              return (int)(g1.end-g2.end);
                        else return strcmp(g1.getID(), g2.getID());
             }
             else return (int)(g1.gseq_id-g2.gseq_id);
}

char* GffLine::extractAttr(const char* pre, bool caseStrict, bool enforce_GTF2) {
 //parse a key attribute and remove it from the info string
 //(only works for attributes that have values following them after ' ' or '=')
 static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required) at GTF line:\n%s\n";
 int lpre=strlen(pre);
 char cend=pre[lpre-1];
 char* pos = (caseStrict) ? strstr(info, pre) : strifind(info, pre);
 if (pos==NULL) return NULL;
 char* findstart=info;
 //require word boundary on the left:
 while (pos!=NULL && pos!=info && *(pos-1)!=';' && *(pos-1)!=' ') {
    findstart=pos+lpre;
    pos = (caseStrict) ? strstr(findstart, pre) : strifind(findstart, pre);
    }
 if (pos==NULL) return NULL;
 if (cend!=' ' && cend!='=') {
    //require word boundary on the right:
    while (pos!=NULL && *(pos+lpre)!=' ' && *(pos+lpre)!='=') {
       findstart=pos+lpre;
       pos = (caseStrict) ? strstr(findstart, pre) : strifind(findstart, pre);
       }
    }
 if (pos==NULL) return NULL;
 char* vp=pos+lpre;
 while (*vp==' ') vp++;
 if (*vp==';' || *vp==0)
      GError("Error parsing value of GFF attribute \"%s\", line:\n%s\n", pre, dupline);
 bool dq_enclosed=false; //value string enclosed by double quotes
 if (*vp=='"') {
     dq_enclosed=true;
     vp++;
     }
 if (enforce_GTF2 && !dq_enclosed)
      GError(GTF2_ERR,pre, dupline);
 char* vend=vp;
 if (dq_enclosed) {
    while (*vend!='"' && *vend!=';' && *vend!=0) vend++;
    }
 else {
    while (*vend!=';' && *vend!=0) vend++;
    }
 if (enforce_GTF2 && *vend!='"')
     GError(GTF2_ERR, pre, dupline);
 char *r=Gstrdup(vp, vend-1);
 //-- now remove this attribute from the info string
 while (*vend!=0 && (*vend=='"' || *vend==';' || *vend==' ')) vend++;
 if (*vend==0) vend--;
 for (char *src=vend, *dest=pos;;src++,dest++) {
   *dest=*src;
   if (*src==0) break;
   }
 return r;
}

static char fnamelc[128];

GffLine::GffLine(GffReader* reader, const char* l) {
 llen=strlen(l);
 GMALLOC(line,llen+1);
 memcpy(line, l, llen+1);
 GMALLOC(dupline, llen+1);
 memcpy(dupline, l, llen+1);
 skip=true;
 gseqname=NULL;
 track=NULL;
 ftype=NULL;
 info=NULL;
 _parents=NULL;
 _parents_len=0;
 num_parents=0;
 parents=NULL;
 is_gff3=false;
 is_cds=false;
 is_transcript=false;
 is_exon=false;
 is_gene=false;
 exontype=0;
 gene_id=NULL;
 gene_name=NULL;
 qstart=0;
 qend=0;
 qlen=0;
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
 if (!parseUInt(p,fstart)) {
   //FIXME: chromosome_band entries in Flybase
   GMessage("Warning: invalid start coordinate at line:\n%s\n",l);
   return;
   }
 p=t[4];
 if (!parseUInt(p,fend)) {
   GMessage("Warning: invalid end coordinate at line:\n%s\n",l);
   return;
   }
 if (fend<fstart) swap(fend,fstart); //make sure fstart>=fend, always
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
 phase=*t[7]; // must be '.', '0', '1' or '2'
 ID=NULL;
 // exon/CDS/mrna filter
 strncpy(fnamelc, ftype, 127);
 fnamelc[127]=0;
 strlower(fnamelc); //convert to lower case
 bool is_t_data=false;
 if (strstr(fnamelc, "utr")!=NULL) {
   exontype=exgffUTR;
   is_exon=true;
   is_t_data=true;
   }
  else if (endsWith(fnamelc, "exon")) {
   exontype=exgffExon;
   is_exon=true;
   is_t_data=true;
   }
  else if (strstr(fnamelc, "stop") && 
      (strstr(fnamelc, "codon") || strstr(fnamelc, "cds"))){
   exontype=exgffStop;
   is_cds=true; //though some place it outside the last CDS segment
   is_t_data=true;
   }
  else if (strstr(fnamelc, "start") && 
      ((strstr(fnamelc, "codon")!=NULL) || strstr(fnamelc, "cds")!=NULL)){
   exontype=exgffStart;
   is_cds=true;
   is_t_data=true;
   }
 else if (strcmp(fnamelc, "cds")==0) {
   exontype=exgffCDS;
   is_cds=true;
   is_t_data=true;
   }
 else if (endsWith(fnamelc, "gene") || startsWith(fnamelc, "gene")) {
   is_gene=true;
   is_t_data=true; //because its name will be attached to parented transcripts
   }
 else if (endsWith(fnamelc,"rna") || endsWith(fnamelc,"transcript")) {
   is_transcript=true;
   is_t_data=true;
   }

if (reader->transcriptsOnly && !is_t_data) {
        char* id=extractAttr("ID=");
        if (id==NULL) id=extractAttr("transcript_id");
        //GMessage("Discarding non-transcript line:\n%s\n",l);
        if (id!=NULL) {
          reader->discarded_ids.Add(id, new int(1));
          GFREE(id);
          }
        return; //skip this line, unwanted feature name
        }
 ID=extractAttr("ID=");
 char* Parent=extractAttr("Parent=");
 is_gff3=(ID!=NULL || Parent!=NULL);
 if (is_gff3) {
   //parse as GFF3
    if (ID!=NULL) {
       //has ID attr so it's likely to be a parent feature
       //look for explicit gene name
       gene_name=extractAttr("gene_name=",false);
       if (gene_name==NULL) {
           gene_name=extractAttr("geneName=",false);
           if (gene_name==NULL) {
               gene_name=extractAttr("gene_sym=",false);
               if (gene_name==NULL) {
                   gene_name=extractAttr("gene=",false);
                   }
               }
           }
       gene_id=extractAttr("geneID=",false);
       if (gene_id==NULL) {
          gene_id=extractAttr("gene_id=",false);
          }
       if (is_gene) {
         //special case: keep the Name and ID attributes of the gene feature
         if (gene_name==NULL)
              gene_name=extractAttr("Name=");
         if (gene_id==NULL) //the ID is also gene_id in this case
              gene_id=Gstrdup(ID);
         //skip=false;
         //return;
         GFREE(Parent); //TMI, we really don't care about gene Parents?
         } //gene feature
       }// has GFF3 ID
   if (Parent!=NULL) {
        //keep Parent attr
         //parse multiple parents
         num_parents=1;
         p=Parent;
         int last_delim_pos=-1;
         while (*p!=';' && *p!=0) {
             if (*p==',' && *(p+1)!=0 && *(p+1)!=';') {
                 num_parents++;
                 last_delim_pos=(p-Parent);
                 }
             p++;
             }
         _parents_len=p-Parent+1;
         _parents=Parent;
         GMALLOC(parents, num_parents*sizeof(char*));
         parents[0]=_parents;
         int i=1;
         if (last_delim_pos>0) {
           for (p=_parents+1;p<=_parents+last_delim_pos;p++) {
              if (*p==',') {
                 char* ep=p-1;
                 while (*ep==' ' && ep>_parents) ep--;
                 *(ep+1)=0; //end the string there
                 parents[i]=p+1;
                 i++;
                 }
              }
           }
         } //has Parent field
   } //GFF3
  else { // GTF-like expected
   Parent=extractAttr("transcript_id");
   if (Parent!=NULL) { //GTF2 format detected
     if (is_transcript) {
         // atypical GTF with a parent transcript line declared
         ID=Parent;
         Parent=NULL;
         }
     gene_id=extractAttr("gene_id"); // for GTF this is the only attribute accepted as geneID
     gene_name=extractAttr("gene_name");
     if (gene_name==NULL) {
           gene_name=extractAttr("gene_sym");
           if (gene_name==NULL)
               gene_name=extractAttr("gene");
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
         else nsp=true; //non-space
       if (*p==';') { noed=true; nsp=false; }
       p++;
       }
     } //GTF2 detected (no parent line)
    else {// Parent is NULL, check for jigsaw format or other pre-GTF2 format
     //char* fexon=strstr(fnamelc, "exon");
     //if (fexon!=NULL) {
     if (exontype==exgffExon) {
       if (startsWith(track,"jigsaw")) {
          is_cds=true;
          strcpy(track,"jigsaw");
          p=strchr(info,';');
          if (p==NULL) { Parent=Gstrdup(info); info=NULL; }
           else { Parent=Gstrdup(info,p-1);
                  info=p+1;
                }
          }
        } //exon feature?
        if (Parent==NULL && exontype>=exgffCDS &&
               (i=strcspn(info,"; \t\n\r"))<=(int)(strlen(info)+1)) {
          //one word ID ? really desperate attempt to parse it here
          Parent=Gstrdup(info,info+i-1);
          info=NULL; //discard anything else on the line
          }
     }
   if (Parent!=NULL) { //GTF transcript_id for exon/CDS feature
      _parents=Parent;
      GMALLOC(parents,sizeof(char*));
      num_parents=1;
      parents[0]=_parents;
      }
   } //GTF-like

 //parse other potentially useful features
 if (is_gff3) {
   if ((p=strstr(info,"Target="))!=NULL) { //has Target attr
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
   if ((p=strifind(info,"Qreg="))!=NULL) { //has Qreg attr
       p+=5;
       if (!parseUInt(p,qstart))
         GError("Error parsing target start coordinate from GFF line:\n%s\n",l);
       if (*p!='-') {
          GError("Error parsing next target coordinate from GFF line:\n%s\n",l);
          }
       p++;
       if (!parseUInt(p,qend))
         GError("Error parsing target end coordinate from GFF line:\n%s\n",l);
       if (*p=='|' || *p==':') {
         p++;
         if (!parseUInt(p,qlen))
           GError("Error parsing target length from GFF Qreg|: \n%s\n",l);
         }
       }//has Qreg attr
   if (qlen==0 && (p=strifind(info,"Qlen="))!=NULL) {
     p+=5;
     if (!parseUInt(p,qlen))
         GError("Error parsing target length from GFF Qlen:\n%s\n",l);
     }
   }//parsing some useful attributes in GFF3 records
 if (ID==NULL && parents==NULL) {
      if (reader->gff_warns)
          GMessage("Warning: could not parse ID or Parent from GFF line:\n%s\n",dupline);
      return; //skip
      }
 skip=false;
}


void GffObj::addCDS(uint cd_start, uint cd_end, char phase) {
  if (cd_start>=this->start) {
        this->CDstart=cd_start;
        if (strand=='+') this->CDphase=phase;
        }
      else this->CDstart=this->start;
  if (cd_end<=this->end) {
      this->CDend=cd_end;
      if (strand=='-') this->CDphase=phase;
      }
     else this->CDend=this->end;
  isTranscript(true);
  exon_ftype_id=gff_fid_exon;
  if (monoFeature()) {
     if (exons.Count()==0) addExon(this->start, this->end,0,'.',0,0,false,exgffExon);
            else exons[0]->exontype=exgffExon;
     }
}

int GffObj::addExon(GffReader* reader, GffLine* gl, bool keepAttr, bool noExonAttr) {
  //this will make sure we have the right subftype_id!
  //int subf_id=-1;
  if (!isTranscript() && gl->is_cds) {
          isTranscript(true);
          exon_ftype_id=gff_fid_exon;
          if (exons.Count()==1) exons[0]->exontype=exgffExon;
          }
  if (isTranscript()) {
     if (exon_ftype_id<0) {//exon_ftype_id=gff_fid_exon;
          if (gl->exontype>0) exon_ftype_id=gff_fid_exon;
                         else exon_ftype_id=names->feats.addName(gl->ftype);
          }
     //any recognized mRNA segment gets the generic "exon" type (also applies to CDS)
     if (gl->exontype==0 && !gl->is_transcript) {
          //extraneous mRNA feature, discard
          if (reader->gff_warns)
            GMessage("Warning: discarding unrecognized transcript subfeature %s of %s\n", 
                gl->ftype, gffID);
          return -1;
          }
     }
  else { //non-mRNA parent feature, check this subf type
    int subf_id=names->feats.addName(gl->ftype);
    if (exon_ftype_id<0 || exons.Count()==0) //never assigned a subfeature type before (e.g. first exon being added)
       exon_ftype_id=subf_id;
     else {
       if (exon_ftype_id!=subf_id) {
         //
         if (exon_ftype_id==ftype_id && exons.Count()==1 && exons[0]->start==start && exons[0]->end==end) {
            //the existing exon was just a dummy one created by default, discard it
            exons.Clear();
            covlen=0;
            exon_ftype_id=subf_id; //allow the new subfeature to completely takeover
            }
         else { //multiple subfeatures, prefer those with 
             if (reader->gff_warns)
               GMessage("GFF Warning: multiple subfeatures (%s and %s) found for %s, discarding ", 
                  names->feats.getName(subf_id), names->feats.getName(exon_ftype_id),gffID);
            if (gl->exontype!=0) { //new feature is an exon, discard previously parsed subfeatures
               if (reader->gff_warns) GMessage("%s.\n", names->feats.getName(exon_ftype_id));
               exon_ftype_id=subf_id;
               exons.Clear();
               covlen=0;
               }
              else { //discard new feature
               if (reader->gff_warns) GMessage("%s.\n", names->feats.getName(subf_id));
               return -1; //skip this 2nd subfeature type for this parent!
               }
            }
         } //incoming subfeature is of different type
       } //new subfeature type
    } //non-mRNA parent
  int eidx=addExon(gl->fstart, gl->fend, gl->score, gl->phase,
         gl->qstart,gl->qend, gl->is_cds, gl->exontype);
  if (eidx<0) return eidx; //this should never happen
  if (keepAttr) {
     if (noExonAttr) { 
         if (attrs==NULL) //place the parsed attributes directly at transcript level
           parseAttrs(attrs, gl->info);
         }
       else { //need all exon-level attributes
         parseAttrs(exons[eidx]->attrs, gl->info, true);
         }
      }
  return eidx;
}


int GffObj::addExon(uint segstart, uint segend, double sc, char fr, int qs, int qe, bool iscds, char exontype) {
  if (exons.Count()==0) {
      if (iscds) isCDS=true; //for now, assume CDS only if first "exon" given is a CDS
      if (exon_ftype_id<0) {
         exon_ftype_id = isTranscript() ? gff_fid_exon : ftype_id;
         }
      }
  //special treatment of start/stop codon features, they might be broken/split between exons 
  //and in that case some providers will still give the wrong end coordinate as start+2 (e.g. UCSC)
  //so we should not trust the end coordinate for such features
  if (exontype==exgffStart || exontype==exgffStop) {
     if (strand=='-') segstart=segend;
                else  segend=segstart;
     if (exontype==exgffStart) {
           if (CDstart==0 || segstart<CDstart) CDstart=segstart;
           }
         else {
           if (segstart>CDend) CDend=segstart;
           }
     }
    else if (iscds) { //update CDS anchors:
     if (CDstart==0 || segstart<CDstart)  {
           CDstart=segstart;
           if (exontype==exgffCDS && strand=='+') CDphase=fr;
           }
     if (segend>CDend) {
           if (exontype==exgffCDS && strand=='-') CDphase=fr;
           CDend=segend;
           }
     }
   else { // not a CDS/start/stop 
     isCDS=false;
     }
  if (qs || qe) {
    if (qs>qe) swap(qs,qe);
    if (qs==0) qs=1;
	}
  int ovlen=0;
  if (exontype>0) { //check for overlaps between exon-type segments
      int oi=exonOverlapIdx(segstart, segend, &ovlen);
      if (oi>=0) { //overlap existing segment
         if (ovlen==0) {
			  //adjacent segments will be merged
			  //e.g. CDS to (UTR|exon)
			  if ((exons[oi]->exontype>=exgffUTR && exontype==exgffCDS) ||
				  (exons[oi]->exontype==exgffCDS && exontype>=exgffUTR)) {
					expandExon(oi, segstart, segend, exgffCDSUTR, sc, fr, qs, qe);
					return oi;
					}
			  //CDS adjacent to stop_codon: UCSC does (did?) this
			  if ((exons[oi]->exontype==exgffStop && exontype==exgffCDS) ||
				  (exons[oi]->exontype==exgffCDS && exontype==exgffStop)) {
					expandExon(oi, segstart, segend, exgffCDS, sc, fr, qs, qe);
					return oi;
					}
             }
		 //only allow this for CDS within exon, stop_codon within (CDS|UTR|exon),
         //                   start_codon within (CDS|exon)
        if (exons[oi]->exontype>exontype && 
             exons[oi]->start<=segstart && exons[oi]->end>=segend &&
             !(exons[oi]->exontype==exgffUTR && exontype==exgffCDS)) {
              //larger segment given first, now the smaller included one is redundant
              return oi; //only used to store attributes from current GffLine
              }
        if (exontype>exons[oi]->exontype && 
             segstart<=exons[oi]->start && segend>=exons[oi]->end &&
             !(exontype==exgffUTR && exons[oi]->exontype==exgffCDS)) {
               //smaller segment given first, so we have to enlarge it
			  expandExon(oi, segstart, segend, exontype, sc, fr, qs, qe);
				//this should also check for overlapping next exon (oi+1) ?
              return oi;
              }
        //there is also the special case of "ribosomal slippage exception" (programmed frameshift)
        //where two CDS segments may actually overlap for 1 or 2 bases, but there should be only one encompassing exon
		//if (ovlen>2 || exons[oi]->exontype!=exgffCDS || exontype!=exgffCDS) {
		// had to relax this because of some weird UCSC annotations with exons partially overlapping the CDS segments
		/*
		if (ovlen>2 && exons[oi]->exontype!=exgffUTR && exontype!=exgffUTR) {
		   if (gff_show_warnings)
			   GMessage("GFF Warning: discarding overlapping feature segment (%d-%d) (vs %d-%d (%s)) for GFF ID %s on %s\n",
			   segstart, segend, exons[oi]->start, exons[oi]->end, getSubfName(), gffID, getGSeqName());
		   hasErrors(true);
		   return -1; //segment NOT added
		   }
		*/

		 if ((ovlen>2 || ovlen==0) || exons[oi]->exontype!=exgffCDS || exontype!=exgffCDS) {
		  if (gff_show_warnings)
			 GMessage("GFF Warning: merging overlapping/adjacent feature segment (%d-%d) into (%d-%d) (%s) for GFF ID %s on %s\n",
				 segstart, segend, exons[oi]->start, exons[oi]->end, getSubfName(), gffID, getGSeqName());
			expandExon(oi, segstart, segend, exontype, sc, fr, qs, qe);
			return oi;
		 }
		// else add the segment if the overlap is small and between two CDS segments
		//TODO: we might want to add an attribute here with the slippage coordinate and size?
        covlen-=ovlen;
		}//overlap or adjacent to existing segment
	   } //check for overlap
   // --- no overlap, or accepted micro-overlap (ribosomal slippage)
   // create & add the new segment
   /*
   if (start>0 && exontype==exgffCDS && exons.Count()==0) {
      //adding a CDS directly as the first subfeature of a declared parent
      segstart=start;
      segend=end;
      } 
   */
   GffExon* enew=new GffExon(segstart, segend, sc, fr, qs, qe, exontype);
   int eidx=exons.Add(enew);
   if (eidx<0) {
    //this would actually be acceptable if the object is a "Gene" and "exons" are in fact isoforms
     if (gff_show_warnings) 
       GMessage("GFF Warning: failed adding segment %d-%d for %s (discarded)!\n",
            segstart, segend, gffID);
     delete enew;
     hasErrors(true);
     return -1;            
     }
   covlen+=(int)(exons[eidx]->end-exons[eidx]->start)+1;
   //adjust parent feature coordinates to contain this exon
   if (start==0 || start>exons.First()->start) {
     start=exons.First()->start;
     }
   if (end<exons.Last()->end) end=exons.Last()->end;
     
   if (uptr!=NULL) { //collect stats about the underlying genomic sequence
       GSeqStat* gsd=(GSeqStat*)uptr;
       if (start<gsd->mincoord) gsd->mincoord=start;
       if (end>gsd->maxcoord) gsd->maxcoord=end;
       if (this->len()>gsd->maxfeat_len) {
          gsd->maxfeat_len=this->len();
          gsd->maxfeat=this;
          }
       }
   return eidx;
}

void GffObj::expandExon(int oi, uint segstart, uint segend, char exontype, double sc, char fr, int qs, int qe) {
  //oi is the index of the *first* overlapping segment found that must be enlarged
  covlen-=exons[oi]->len();
  if (segstart<exons[oi]->start)
    exons[oi]->start=segstart;
  if (qs && qs<exons[oi]->qstart) exons[oi]->qstart=qs;
  if (segend>exons[oi]->end)
    exons[oi]->end=segend;
  if (qe && qe>exons[oi]->qend) exons[oi]->qend=qe;
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
/* I guess we have to relax this due to stupid UCSC hg18 files having a start_codon sticking out
chr1	hg18_knownGene	start_codon	69806911	69806913	0.000000	+	.
chr1	hg18_knownGene	CDS	69806911	69806912	0.000000	+	0
chr1	hg18_knownGene	exon	69805456	69806912	0.000000	+	.     
*/
         if (exons[ni]->qstart<exons[oi]->qstart) exons[oi]->qstart=exons[ni]->qstart;
         if (exons[ni]->qend>exons[oi]->qend) exons[oi]->qend=exons[ni]->qend;
         exons.Delete(ni);
         }
      else {
         if (gff_show_warnings) GMessage("GFF Warning: overlapping existing exon(%d-%d) while expanding to %d-%d for GFF ID %s\n",
                exons[ni]->start, exons[ni]->end, segstart, segend, gffID);
         //hasErrors(true);
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
    if (this->len()>gsd->maxfeat_len) {
        gsd->maxfeat_len=this->len();
        gsd->maxfeat=this;
        }
    }
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

void GffObj::removeExon(GffExon* p) {
  for (int idx=0;idx<exons.Count();idx++) {
     if (exons[idx]==p) {
        int segstart=exons[idx]->start;
        int segend=exons[idx]->end;
        exons.Delete(idx);
        covlen -= (int)(segend-segstart)+1;
        start=exons.First()->start;
        end=exons.Last()->end;
        if (isCDS) { CDstart=start; CDend=end; }
        return;
        }
     }
}



GffObj::GffObj(GffReader *gfrd, GffLine* gffline, bool keepAttr, bool noExonAttr):
     GSeg(0,0), exons(true,true,false), children(1,false) {
  xstart=0;
  xend=0;
  xstatus=0;
  partial=false;
  isCDS=false;
  uptr=NULL;
  ulink=NULL;
  parent=NULL;
  udata=0;
  flags=0;
  CDstart=0;
  CDend=0;
  CDphase=0;
  geneID=NULL;
  gene_name=NULL;
  attrs=NULL;
  gffID=NULL;
  track_id=-1;
  gseq_id=-1;
  ftype_id=-1;
  exon_ftype_id=-1;
  strand='.';
  if (gfrd==NULL)
    GError("Cannot use this GffObj constructor with a NULL GffReader!\n");
  gffnames_ref(names);
  if (gfrd->names==NULL) gfrd->names=names;
  //qlen=0;qstart=0;qend=0;
  gscore=0;
  uscore=0;
  covlen=0;
  qcov=0;
  start=gffline->fstart;
  end=gffline->fend;
  gseq_id=names->gseqs.addName(gffline->gseqname);
  track_id=names->tracks.addName(gffline->track);
  strand=gffline->strand;
  qlen=gffline->qlen;
  qstart=gffline->qstart;
  qend=gffline->qend;
  //setup flags from gffline
  isCDS=gffline->is_cds; //for now
  isGene(gffline->is_gene);
  isTranscript(gffline->is_transcript || gffline->exontype!=0);
  fromGff3(gffline->is_gff3);

  if (gffline->parents!=NULL) {
    //GTF style -- create a GffObj directly by subfeature
    //(also possible orphan GFF3 exon line, or an exon given before its parent (chado))
    if (gffline->exontype!=0) { //recognized exon-like feature
       ftype_id=gff_fid_transcript; //so this is some sort of transcript
       exon_ftype_id=gff_fid_exon; //subfeatures MUST be exons
       }
     else {//unrecognized subfeatures
       //make this GffObj of the same feature type
       ftype_id=names->feats.addName(gffline->ftype);
       }
    if (gffline->ID==NULL) { //typical GTF
        gffID=Gstrdup(gffline->parents[0]);
        this->createdByExon(true);
        //this is likely the first exon/segment of the feature
        addExon(gfrd, gffline, keepAttr, noExonAttr);
        }
      else { //a parented feature with an ID -- probably an orphan GFF3 line
        if (gffline->is_gff3 && gffline->exontype!=0) {
             //premature exon given before its parent transcript
             //create the transcript entry here
             gffID=Gstrdup(gffline->parents[0]);
             this->createdByExon(true);
             //this is the first exon/segment of the transcript
             addExon(gfrd, gffline, keepAttr, noExonAttr);
             }
            else { //unrecognized non-exon feature ? use the ID instead
             gffID=Gstrdup(gffline->ID);
             if (keepAttr) this->parseAttrs(attrs, gffline->info);
             }
        }
    } //subfeature given directly
  else { //gffline->parents==NULL
    //create a parent feature in its own right
    gscore=gffline->score;
    if (gffline->ID==NULL || gffline->ID[0]==0)
      GError("Error: no ID found for GFF record start\n");
    gffID=Gstrdup(gffline->ID); //there must be an ID here
    //if (gffline->is_transcript) ftype_id=gff_fid_mRNA;
      //else
    ftype_id=names->feats.addName(gffline->ftype);
    if (gffline->is_transcript)
      exon_ftype_id=gff_fid_exon;

    if (keepAttr) this->parseAttrs(attrs, gffline->info);
    }//no parent

  if (gffline->gene_name!=NULL) {
     gene_name=Gstrdup(gffline->gene_name);
     }
  if (gffline->gene_id!=NULL) {
     geneID=Gstrdup(gffline->gene_id);
     }

  GSeqStat* gsd=gfrd->gseqstats.AddIfNew(new GSeqStat(gseq_id,names->gseqs.lastNameUsed()),true);
  uptr=gsd;
  if (start<gsd->mincoord) gsd->mincoord=start;
  if (end>gsd->maxcoord) gsd->maxcoord=end;
    if (this->len()>gsd->maxfeat_len) {
        gsd->maxfeat_len=this->len();
        gsd->maxfeat=this;
        }
}

GffLine* GffReader::nextGffLine() {
 if (gffline!=NULL) return gffline; //caller should free gffline after processing
 while (gffline==NULL) {
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
       continue;
       }
    if (gffline->ID==NULL && gffline->parents==NULL)  { //it must have an ID
        //this might not be needed, already checked in the GffLine constructor
        if (gff_warns)
            GMessage("Warning: malformed GFF line, no parent or record Id (kipping\n");
        delete gffline;
        gffline=NULL;
        //continue;
        }
    }
return gffline;
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

//Warning: if gflst gets altered, idx becomes obsolete
GfoHolder* GffReader::gfoAdd(const char* id, const char* ctg, GffObj* gfo, int idx) {
 char* buf=gfoBuildId(id,ctg);
 GfoHolder* r=new GfoHolder(gfo,idx);
 phash.Add(buf, r);
 GFREE(buf);
 return r;
}

GfoHolder* GffReader::gfoFind(const char* id, const char* ctg) {
 char* buf=gfoBuildId(id,ctg);
 GfoHolder* r=phash.Find(buf);
 GFREE(buf);
 return r;
}

GfoHolder* GffReader::replaceGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr, int replaceidx) {
  GffObj* newgfo=new GffObj(this, gffline, keepAttr, noExonAttr);
  GfoHolder* r=NULL;
  if (replaceidx>=0) {
     gflst.Put(replaceidx,newgfo);
     r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, replaceidx);
     }
   else {
     int gfoidx=gflst.Add(newgfo);
     r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, gfoidx);
     }
  if (gff_warns) {
    int* pcount=tids.Find(newgfo->gffID);
    if (pcount!=NULL) {
       if (gff_warns) GMessage("Warning: duplicate GFF ID: %s\n", newgfo->gffID);
       (*pcount)++;
       }
     else {
       tids.Add(newgfo->gffID,new int(1));
       }
    }
  return r;
}

GfoHolder* GffReader::updateParent(GfoHolder* newgfh, GffObj* parent) {
  //assert(parent);
  //assert(newgfo);
  parent->children.Add(newgfh->gffobj);
  if (newgfh->gffobj->parent==NULL) newgfh->gffobj->parent=parent;
  newgfh->gffobj->setLevel(parent->getLevel()+1);
  if (parent->isGene()) {
      if (parent->gene_name!=NULL && newgfh->gffobj->gene_name==NULL)
        newgfh->gffobj->gene_name=Gstrdup(parent->gene_name);
      if (parent->geneID!=NULL && newgfh->gffobj->geneID==NULL)
        newgfh->gffobj->geneID=Gstrdup(parent->geneID);
      }

  return newgfh;
}

GfoHolder* GffReader::newGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr,
                          GffObj* parent, GffExon* pexon) {
  GffObj* newgfo=new GffObj(this, gffline, keepAttr, noExonAttr);
  GfoHolder* r=NULL;
  int gfoidx=gflst.Add(newgfo);
  r=gfoAdd(newgfo->gffID, gffline->gseqname, newgfo, gfoidx);
  if (parent!=NULL) {
    updateParent(r, parent);
    if (pexon!=NULL) parent->removeExon(pexon);
    }
  if (gff_warns) {
    int* pcount=tids.Find(newgfo->gffID);
    if (pcount!=NULL) {
       if (gff_warns) GMessage("Warning: duplicate GFF ID: %s\n", newgfo->gffID);
       (*pcount)++;
       }
     else {
       tids.Add(newgfo->gffID,new int(1));
       }
    }
  return r;
}

GfoHolder* GffReader::updateGffRec(GfoHolder* prevgfo, GffLine* gffline, 
                                         bool keepAttr) {
 if (prevgfo==NULL) return NULL;
 prevgfo->gffobj->createdByExon(false);
 prevgfo->gffobj->ftype_id=prevgfo->gffobj->names->feats.addName(gffline->ftype);
 prevgfo->gffobj->start=gffline->fstart;
 prevgfo->gffobj->end=gffline->fend;
 prevgfo->gffobj->isGene(gffline->is_gene);
 prevgfo->gffobj->isTranscript(gffline->is_transcript || gffline->exontype!=0);
 prevgfo->gffobj->fromGff3(gffline->is_gff3);
 if (keepAttr) {
   if (prevgfo->gffobj->attrs!=NULL) prevgfo->gffobj->attrs->Clear();
   prevgfo->gffobj->parseAttrs(prevgfo->gffobj->attrs, gffline->info);
   }
 return prevgfo;
}


bool GffReader::addExonFeature(GfoHolder* prevgfo, GffLine* gffline, GHash<CNonExon>& pex, bool noExonAttr) {
  bool r=true;
  if (gffline->strand!=prevgfo->gffobj->strand) {
  //TODO: add support for trans-splicing and even inter-chromosomal fusions
     if (prevgfo->gffobj->strand=='.') {
            prevgfo->gffobj->strand=gffline->strand;
        }
     else {
       GMessage("GFF Error at %s (%c): exon %d-%d (%c) found on different strand; discarded.\n",
       prevgfo->gffobj->gffID, prevgfo->gffobj->strand,
       gffline->fstart, gffline->fend, gffline->strand, prevgfo->gffobj->getGSeqName());
       //r=false;
       return true; //FIXME: split trans-spliced mRNAs by strand
       }
   }
  int gdist=(gffline->fstart>prevgfo->gffobj->end) ? gffline->fstart-prevgfo->gffobj->end :
                      ((gffline->fend<prevgfo->gffobj->start)? prevgfo->gffobj->start-gffline->fend :
                         0 );
  if (gdist>(int)GFF_MAX_LOCUS) { //too far apart, most likely this is a duplicate ID
    GMessage("Error: duplicate GFF ID '%s' (or exons too far apart)!\n",prevgfo->gffobj->gffID);
    //validation_errors = true;
    r=false;
    if (!gff_warns) exit(1);
    }
  int eidx=prevgfo->gffobj->addExon(this, gffline, !noExonAttr, noExonAttr);
  if (eidx>=0 && gffline->ID!=NULL && gffline->exontype==0)
      subfPoolAdd(pex, prevgfo);
  return r;
}

CNonExon* GffReader::subfPoolCheck(GffLine* gffline, GHash<CNonExon>& pex, char*& subp_name) {
  CNonExon* subp=NULL;
  subp_name=NULL;
  for (int i=0;i<gffline->num_parents;i++) {
    if (transcriptsOnly && discarded_ids.Find(gffline->parents[i])!=NULL)
        continue;
    subp_name=gfoBuildId(gffline->parents[i], gffline->gseqname); //e.g. mRNA name
    subp=pex.Find(subp_name);
    if (subp!=NULL)
       return subp;
    GFREE(subp_name);
    }
  return NULL;
}

void GffReader::subfPoolAdd(GHash<CNonExon>& pex, GfoHolder* newgfo) {
//this might become a parent feature later
if (newgfo->gffobj->exons.Count()>0) {
   char* xbuf=gfoBuildId(gffline->ID, gffline->gseqname);
   pex.Add(xbuf, new CNonExon(newgfo->idx, newgfo->gffobj,
       newgfo->gffobj->exons[0], gffline));
   GFREE(xbuf);
   }
}

GfoHolder* GffReader::promoteFeature(CNonExon* subp, char*& subp_name, GHash<CNonExon>& pex,
    bool keepAttr, bool noExonAttr) {
  GffObj* prevp=subp->parent; //grandparent of gffline (e.g. gene)
  if (prevp!=gflst[subp->idx])
    GError("Error promoting subfeature %s, gflst index mismatch?!\n", subp->gffline->ID);
  subp->gffline->discardParent();
  GfoHolder* gfoh=newGffRec(subp->gffline, keepAttr, noExonAttr, prevp, subp->exon);
  pex.Remove(subp_name); //no longer a potential parent, moved it to phash already
  prevp->promotedChildren(true);
  return gfoh; //returns the holder of newly promoted feature
}

//have to parse the whole file because exons can be scattered all over
void GffReader::readAll(bool keepAttr, bool mergeCloseExons, bool noExonAttr) {
  bool validation_errors = false;
  //loc_debug=false;
  GHash<CNonExon> pex; //keep track of any "exon"-like features that have an ID
                     //and thus could become promoted to parent features
  while (nextGffLine()!=NULL) {
       //seen this gff ID before?
     GfoHolder* prevseen=NULL;
     if (gffline->ID) //GFF3
         prevseen=gfoFind(gffline->ID, gffline->gseqname);
     if (prevseen!=NULL) {
            if (prevseen->gffobj->createdByExon()) {
                updateGffRec(prevseen, gffline, keepAttr);
                }
             else {
                GMessage("Error: duplicate GFF ID '%s' encountered!\n",gffline->ID);
                validation_errors = true;
                if (gff_warns) { 
                       delete gffline; gffline=NULL; continue;
                       }
                   else exit(1);
                }
            }
    if (gffline->parents==NULL) {//start GFF3-like record with no parent (mRNA, gene)
       if (!prevseen) newGffRec(gffline, keepAttr, noExonAttr);
       }
    else { //--- it's a parented feature (could still be a mRNA)
       bool found_parent=false;
       GfoHolder* newgfo=prevseen;
       for (int i=0;i<gffline->num_parents;i++) {
            if (transcriptsOnly && discarded_ids.Find(gffline->parents[i])!=NULL)
                continue; //skipping discarded parent feature
            GfoHolder* parentgfo=gfoFind(gffline->parents[i], gffline->gseqname);
            if (parentgfo!=NULL) { //parent GffObj parsed earlier
                   found_parent=true;
                   if (parentgfo->gffobj->isGene() && gffline->is_transcript
                                   && gffline->exontype==0) {
                       //not an exon, but a transcript parented by a gene
                       if (newgfo) {
                           updateParent(newgfo, parentgfo->gffobj);
                           }
                         else {
                           newgfo=newGffRec(gffline, keepAttr, noExonAttr, parentgfo->gffobj);
                           }
                       }
                   else { //potential exon subfeature
                       if (!addExonFeature(parentgfo, gffline, pex, noExonAttr))
                         validation_errors=true;
                       }
                   }
            } //for each parsed parent Id
       if (!found_parent) { //new GTF-like record starting here with a subfeature directly
             //or it could be some chado GFF3 barf with exons declared BEFORE their parent :(
            //check if this feature isn't parented by a previously stored "exon" subfeature
            char* subp_name=NULL;
            CNonExon* subp=subfPoolCheck(gffline, pex, subp_name);
            if (subp!=NULL) { //found a subfeature that is the parent of this gffline
               //promote that subfeature to a full GffObj
               GfoHolder* gfoh=promoteFeature(subp, subp_name, pex, keepAttr, noExonAttr);
               //add current gffline as an exon of the newly promoted subfeature
               if (!addExonFeature(gfoh, gffline, pex, noExonAttr))
                      validation_errors=true;
               }
              else { //no parent seen before, create one directly with this exon
               //loc_debug=true;
               GfoHolder* newgfo=prevseen ? prevseen : newGffRec(gffline, keepAttr, noExonAttr);
               if (gffline->ID!=NULL && gffline->exontype==0)
                     subfPoolAdd(pex, newgfo);
               //even those with errors will be added here!
               }
            GFREE(subp_name);
            } //no previous parent found
       } //parented feature
        //--
      delete gffline;
      gffline=NULL;
      }//while gff lines
  gflst.finalize(this, mergeCloseExons, keepAttr, noExonAttr); //force sorting by locus if so constructed
 // all gff records are now loaded in GList gflst
 // so we can free the hash
  phash.Clear();
  tids.Clear();
  if (validation_errors) {
    exit(1);
    }
}

GffObj* GffObj::finalize(GffReader* gfr, bool mergeCloseExons, bool keepAttrs, bool noExonAttr) {
 //merge
 //always merge adjacent or overlapping segments
 //but if mergeCloseExons then merge even when distance is up to 5 bases
 udata=0;
 uptr=NULL;
 if (gfr->transcriptsOnly && !(isTranscript() || (isGene() && children.Count()==0))) {
       isDiscarded(true);
       }
 if (ftype_id==gff_fid_transcript && CDstart>0) {
    ftype_id=gff_fid_mRNA;
    //exon_ftype_id=gff_fid_exon;
    }
 //if (ftype_id==gff_fid_mRNA || exon_ftype_id==gff_fid_exon || mergeCloseExons) {
 if (isTranscript() || exon_ftype_id==gff_fid_exon || mergeCloseExons) {
   int mindist=mergeCloseExons ? 5:1;
   for (int i=0;i<exons.Count()-1;i++) {
     int ni=i+1;
     uint mend=exons[i]->end;
     while (ni<exons.Count()) {
       int dist=(int)(exons[ni]->start-mend);
       if (dist>mindist) break; //no merging with next segment
       if (gfr!=NULL && gfr->gff_warns && dist!=0 && (exons[ni]->exontype!=exgffUTR && exons[i]->exontype!=exgffUTR)) {
          GMessage("GFF warning: merging adjacent/overlapping segments of %s on %s (%d-%d, %d-%d)\n",
               gffID, getGSeqName(), exons[i]->start, exons[i]->end,exons[ni]->start, exons[ni]->end);
          }
       mend=exons[ni]->end;
       covlen-=exons[i]->len();
       exons[i]->end=mend;
       covlen+=exons[i]->len();
       covlen-=exons[ni]->len();
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
 //attribute reduction for GTF records
 if (keepAttrs && !noExonAttr && !fromGff3() 
          && exons.Count()>0 && exons[0]->attrs!=NULL) {
   bool attrs_discarded=false;
   for (int a=0;a<exons[0]->attrs->Count();a++) {
      int attr_name_id=exons[0]->attrs->Get(a)->attr_id;
      char* attr_name=names->attrs.getName(attr_name_id);
      char* attr_val =exons[0]->attrs->Get(a)->attr_val;
      bool sameExonAttr=true;
      for (int i=1;i<exons.Count();i++) {
         char* ov=exons[i]->getAttr(attr_name_id);
         if (ov==NULL || (strcmp(ov,attr_val)!=0)) { 
             sameExonAttr=false; 
             break;
             }
         }
      if (sameExonAttr) {
             //delete this attribute from exons level
             attrs_discarded=true;
             this->addAttr(attr_name, attr_val);
             for (int i=1;i<exons.Count();i++) {
                 removeExonAttr(*(exons[i]), attr_name_id);
                 }
             exons[0]->attrs->freeItem(a);
             }
      }
   if (attrs_discarded) exons[0]->attrs->Pack();
   }
 return this;
}

void GffObj::parseAttrs(GffAttrs*& atrlist, char* info, bool isExon) {
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
       //if (noExonAttr && (strcmp(start, "exon_number")==0 || strcmp(start, "exon")==0)) { start=pch; continue; }
       if (strcmp(start, "exon_number")==0 || strcmp(start, "exon")==0 ||
              strcmp(start, "exon_id")==0)
           { start=pch; continue; }
       ech++;
       while (*ech==' ' && ech<endinfo) ech++;//skip extra spaces after the '='
       //atrlist->Add(new GffAttr(names->attrs.addName(start),ech));
       //make sure we don't add the same attribute more than once
       if (isExon && (strcmp(start, "protein_id")==0)) {
             //Ensembl special case
             this->addAttr(start, ech);
             start=pch;
             continue;
             }
       atrlist->add_or_update(names, start, ech);
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

void GffObj::addAttr(const char* attrname, const char* attrvalue) {
  if (this->attrs==NULL)
      this->attrs=new GffAttrs();
  //this->attrs->Add(new GffAttr(names->attrs.addName(attrname),attrvalue));
  this->attrs->add_or_update(names, attrname, attrvalue);
}


void GffObj::setFeatureName(const char* feature) {
 //change the feature name/type for a transcript
 int fid=names->feats.addName(feature);
 if (monoFeature() && exons.Count()>0) 
    this->exon_ftype_id=fid;
 this->ftype_id=fid;
}

void GffObj::setRefName(const char* newname) {
 //change the feature name/type for a transcript
 int rid=names->gseqs.addName(newname);
 this->gseq_id=rid;
}



int GffObj::removeAttr(const char* attrname, const char* attrval) {
  if (this->attrs==NULL || attrname==NULL || attrname[0]==0) return 0;
  int aid=this->names->attrs.getId(attrname);
  if (aid<0) return 0;
  int delcount=0;  //could be more than one ?
  for (int i=0;i<this->attrs->Count();i++) {
     if (aid==this->attrs->Get(i)->attr_id) {
       if (attrval==NULL || 
          strcmp(attrval, this->attrs->Get(i)->attr_val)==0) {
             delcount++;
             this->attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) this->attrs->Pack(); 
  return delcount;
}

int GffObj::removeAttr(int aid, const char* attrval) {
  if (this->attrs==NULL || aid<0) return 0;
  int delcount=0;  //could be more than one ?
  for (int i=0;i<this->attrs->Count();i++) {
     if (aid==this->attrs->Get(i)->attr_id) {
       if (attrval==NULL || 
          strcmp(attrval, this->attrs->Get(i)->attr_val)==0) {
             delcount++;
             this->attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) this->attrs->Pack(); 
  return delcount;
}


int GffObj::removeExonAttr(GffExon& exon, const char* attrname, const char* attrval) {
  if (exon.attrs==NULL || attrname==NULL || attrname[0]==0) return 0;
  int aid=this->names->attrs.getId(attrname);
  if (aid<0) return 0;
  int delcount=0;  //could be more than one
  for (int i=0;i<exon.attrs->Count();i++) {
     if (aid==exon.attrs->Get(i)->attr_id) {
       if (attrval==NULL || 
          strcmp(attrval, exon.attrs->Get(i)->attr_val)==0) {
             delcount++;
             exon.attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) exon.attrs->Pack(); 
  return delcount;
}

int GffObj::removeExonAttr(GffExon& exon, int aid, const char* attrval) {
  if (exon.attrs==NULL || aid<0) return 0;
  int delcount=0;  //could be more than one
  for (int i=0;i<exon.attrs->Count();i++) {
     if (aid==exon.attrs->Get(i)->attr_id) {
       if (attrval==NULL || 
          strcmp(attrval, exon.attrs->Get(i)->attr_val)==0) {
             delcount++;
             exon.attrs->freeItem(i);
             }
       }
     }
  if (delcount>0) exon.attrs->Pack(); 
  return delcount;
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
             }
       if (CDend>=sgstart && CDend<=sgend) {
             //CDstart in this segment
             //and we are getting the whole transcript
             cds_mstart=s-(int)(CDend-cdsadj-sgstart);
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

char* GffObj::getUnspliced(GFaSeqGet* faseq, int* rlen, GList<GSeg>* seglst) 
{
    if (faseq==NULL) { GMessage("Warning: getUnspliced(NULL,.. ) called!\n");
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
    char* unspliced=NULL;

    int seqstart=exons.First()->start;
    int seqend=exons.Last()->end;
    
    int unsplicedlen = 0;

    unsplicedlen += seqend - seqstart + 1;

    GMALLOC(unspliced, unsplicedlen+1); //allocate more here
    //uint seqstart, seqend;

    int s = 0; //resulting nucleotide counter
    if (strand=='-') 
    {
        if (seglst!=NULL)
            seglst->Add(new GSeg(s+1,s+1+seqend-seqstart));
        for (int i=seqend;i>=seqstart;i--) 
        {
            unspliced[s] = ntComplement(gsubseq[i-start]);
            s++;
        }//for each nt
    } // - strand
    else 
    { // + strand
        if (seglst!=NULL)
            seglst->Add(new GSeg(s+1,s+1+seqend-seqstart));
        for (int i=seqstart;i<=seqend;i++) 
        {
            unspliced[s]=gsubseq[i-start];
            s++;
        }//for each nt
    } // + strand
    //assert(s <= unsplicedlen);
    unspliced[s]=0;
    if (rlen!=NULL) *rlen=s;
    return unspliced;
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
  if (fspan<(int)(end-start+1)) { //special case: stop coordinate was extended past the gseq length, must adjust
     int endadj=end-start+1-fspan;
     uint prevend=end;
     end-=endadj;
     if (CDend>end) CDend=end;
     if (exons.Last()->end>end) {
         exons.Last()->end=end; //this could get us into trouble if exon start is also > end
         if (exons.Last()->start>exons.Last()->end) {
            GError("GffObj::getSpliced() error: improper genomic coordinate %d on %s for %s\n",
                  prevend,getGSeqName(), getID());
            }
         covlen-=endadj;
         }
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

void GffObj::printGxfLine(FILE* fout, const char* tlabel, const char* gseqname, bool iscds,
                             uint segstart, uint segend, int exidx, char phase, bool gff3) {
  static char scorestr[14];
  strcpy(scorestr,".");
  GffAttrs* xattrs=NULL;
  if (exidx>=0) {
     if (exons[exidx]->score) sprintf(scorestr,"%.2f", exons[exidx]->score);
     xattrs=exons[exidx]->attrs;
  }
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
  else {//for GTF -- we print only transcripts
    //if (isValidTranscript())
    fprintf(fout, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t%c\ttranscript_id \"%s\";",
           gseqname, tlabel, ftype, segstart, segend, scorestr, strand, phase, gffID);
    //char* geneid=(geneID!=NULL)? geneID : gffID;
    if (geneID)
      fprintf(fout," gene_id \"%s\";",geneID);
    if (gene_name!=NULL) {
       //fprintf(fout, " gene_name ");
       //if (gene_name[0]=='"') fprintf (fout, "%s;",gene_name);
       //              else fprintf(fout, "\"%s\";",gene_name);
       fprintf(fout," gene_name \"%s\";",gene_name);
       }
    if (xattrs!=NULL) {
          for (int i=0;i<xattrs->Count();i++) {
            if (xattrs->Get(i)->attr_val==NULL) continue;
            const char* attrname=names->attrs.getName(xattrs->Get(i)->attr_id);
            fprintf(fout, " %s ",attrname);
            if (xattrs->Get(i)->attr_val[0]=='"')
                     fprintf(fout, "%s;",xattrs->Get(i)->attr_val);
                else fprintf(fout, "\"%s\";",xattrs->Get(i)->attr_val);
             }
          }
    //for GTF, also append the GffObj attributes to each exon line
    if ((xattrs=this->attrs)!=NULL) {
          for (int i=0;i<xattrs->Count();i++) {
            if (xattrs->Get(i)->attr_val==NULL) continue;
            const char* attrname=names->attrs.getName(xattrs->Get(i)->attr_id);
            fprintf(fout, " %s ",attrname);
            if (xattrs->Get(i)->attr_val[0]=='"')
                     fprintf(fout, "%s;",xattrs->Get(i)->attr_val);
                else fprintf(fout, "\"%s\";",xattrs->Get(i)->attr_val);
             }
           }
    fprintf(fout, "\n");
    }//GTF
}

void GffObj::printGxf(FILE* fout, GffPrintMode gffp, 
                   const char* tlabel, const char* gfparent) {
 static char tmpstr[255];
 if (tlabel==NULL) {
    tlabel=track_id>=0 ? names->tracks.Get(track_id)->name :
         (char*)"gffobj" ;
    }
 unxcoord();
 //if (exons.Count()==0) return;
 const char* gseqname=names->gseqs.Get(gseq_id)->name;
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
   //const char* ftype=isTranscript() ? "mRNA" : getFeatureName();
   const char* ftype=getFeatureName();
   fprintf(fout,
     "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\tID=%s",
     gseqname, tlabel, ftype, pstart, pend, tmpstr, strand, gffID);
   if (CDstart>0 && !showCDS && !isCDS) fprintf(fout,";CDS=%d-%d",CDstart,CDend);
   if (gfparent!=NULL) {
      //parent override
      fprintf(fout, ";Parent=%s",gfparent);
      }
     else {
       if (parent!=NULL && !parent->isDiscarded())
           fprintf(fout, ";Parent=%s",parent->getID());
       }
   if (geneID!=NULL)
      fprintf(fout, ";geneID=%s",geneID);
   if (gene_name!=NULL)
      fprintf(fout, ";gene_name=%s",gene_name);
   if (attrs!=NULL) {
      for (int i=0;i<attrs->Count();i++) {
        const char* attrname=names->attrs.getName(attrs->Get(i)->attr_id);
        fprintf(fout,";%s=%s", attrname,
               attrs->Get(i)->attr_val);
        }
      }
    fprintf(fout,"\n");
   }// gff3 mRNA line
 if (showExon) {
   //print exons
    if (isCDS && exons.Count()>0 && 
        ((strand=='-' && exons.Last()->phase<'0') || (strand=='+' && exons.Last()->phase<'0')))
         updateExonPhase();

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

void GffObj::updateExonPhase() {
  if (!isCDS) return;
  int cdsacc=0;
  if (CDphase=='1' || CDphase=='2') {
      cdsacc+= 3-(CDphase-'0');
      }
  if (strand=='-') { //reverse strand
     for (int i=exons.Count()-1;i>=0;i--) {
         exons[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=exons[i]->end-exons[i]->start+1;
         }
     }
    else { //forward strand
     for (int i=0;i<exons.Count();i++) {
         exons[i]->phase='0'+ (3-cdsacc%3)%3;
         cdsacc+=exons[i]->end-exons[i]->start+1;
         }
     }
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
