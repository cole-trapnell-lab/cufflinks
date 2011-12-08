#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "SVN"
#endif

#include "GArgs.h"
#include <ctype.h>
#include <errno.h>
#include "gtf_tracking.h"

#ifdef HAVE_CONFIG_H
#include "update_check.h"
#endif

#define USAGE "Usage:\n\
cuffcompare [-r <reference_mrna.gtf>] [-R] [-T] [-V] [-s <seq_path>] \n\
    [-o <outprefix>] [-p <cprefix>] \n\
    {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}\n\
\n\
 Cuffcompare provides classification, reference annotation mapping and various\n\
 statistics for Cufflinks transfrags.\n\
 Cuffcompare clusters and tracks transfrags across multiple samples, writing\n\
 matching transcripts (intron chains) into <outprefix>.tracking, and a GTF\n\
 file <outprefix>.combined.gtf containing a nonredundant set of transcripts \n\
 across all input files (with a single representative transfrag chosen\n\
 for each clique of matching transfrags across samples).\n\
\n\
Options:\n\
-i provide a text file with a list of Cufflinks GTF files to process instead\n\
   of expecting them as command line arguments (useful when a large number\n\
   of GTF files should be processed)\n\
\n\
-r  a set of known mRNAs to use as a reference for assessing \n\
    the accuracy of mRNAs or gene models given in <input.gtf>\n\
\n\
-R  for -r option, reduce the set of reference transcripts to \n\
    only those found to overlap any of the input loci\n\
-M  discard (ignore) single-exon transfrags and reference transcripts\n\
-N  discard (ignore) single-exon reference transcripts\n\
\n\
-s  <seq_path> can be a multi-fasta file with all the genomic sequences or \n\
    a directory containing multiple single-fasta files (one file per contig);\n\
    lower case bases will be used to classify input transcripts as repeats\n\
\n\
-d  max distance (range) for grouping transcript start sites (100)\n\
-p  the name prefix to use for consensus transcripts in the \n\
    <outprefix>.combined.gtf file (default: 'TCONS')\n\
-C  include the \"contained\" transcripts in the .combined.gtf file\n\
-G  generic GFF input file(s) (do not assume Cufflinks GTF)\n\
-T  do not generate .tmap and .refmap files for each input file\n\
-V  verbose processing mode (showing all GFF parsing warnings)\n\
"
bool debug=false;
bool perContigStats=false; // -S to enable stats for every single contig
bool generic_GFF=false; //-G, don't assume Cufflinks GTF as input
bool showContained=false; // -C
bool reduceRefs=false;
bool checkFasta=false;
bool tmapFiles=true;
bool only_spliced_refs=false;
int debugCounter=0;

int polyrun_range=2000; //polymerase run range 2KB
double scoreThreshold=0;
double exprThreshold=0;
char* cprefix=NULL;
FILE* ffasta=NULL; //genomic seq file
FILE *f_ref=NULL; //reference mRNA GFF, if provided
FILE* f_in=NULL; //sequentially, each input GFF file
FILE* f_out=NULL; //stdout if not provided
GFastaHandler gfasta;
int xlocnum=0;
int tsscl_num=0; //for tss cluster IDs
int protcl_num=0; //for "unique" protein IDs within TSS clusters
int tssDist=100;
//int total_tcons=0;
int total_xloci_alt=0;

void openfwrite(FILE* &f, GArgs& args, char opt) {
  GStr s=args.getOpt(opt);
  if (!s.is_empty()) {
      if (s=='-')
       f=stdout;
      else {
       f=fopen(s,"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

//-- structure to keep track of data from multiple qry input files for a single genomic seq
class GSeqTrack {
  int gseq_id;
 public:
  const char* gseq_name;
  GList<GLocus>* rloci_f; //reference loci for this genomic sequence
  GList<GLocus>* rloci_r;
  GList<GXLocus> xloci_f; // extended super-loci across all qry datasets
  GList<GXLocus> xloci_r; // extended super-loci across all qry datasets
  GList<GXLocus> xloci_u; // extended super-loci across all qry datasets
  GSeqData* qdata[MAX_QFILES]; //fixed order array with GSeqData for each qry input
                 //element in array is NULL if a qry file has no transcripts on this genomic sequence
  int get_gseqid() { return gseq_id; }
  GSeqTrack(int gid=-1):xloci_f(true,true,false),
        xloci_r(true,true,false), xloci_u(true,true,false) {
    gseq_id=gid;
    if (gseq_id>=0) {
      gseq_name=GffObj::names->gseqs.getName(gseq_id);
      }
    rloci_f=NULL;
    rloci_r=NULL;
    for (int i=0;i<MAX_QFILES;i++) qdata[i]=NULL;
    }
  bool operator==(GSeqTrack& d){
      return (gseq_id==d.gseq_id);
      }
  bool operator>(GSeqTrack& d){
     return (gseq_id>d.gseq_id);
     }
  bool operator<(GSeqTrack& d){
     return (gseq_id<d.gseq_id);
     }
};

char* getFastaFile(int gseq_id);

// ref globals
bool haveRefs=false;  //true if a reference annotation (-r) is provide

GList<GSeqData> ref_data(true,true,true); //list of reference mRNAs and loci data for each genomic seq
              //each locus will keep track of any superloci which includes it, formed during the analysis

void processLoci(GSeqData& seqdata, GSeqData* refdata=NULL, int qfidx=0);

void reportStats(FILE* fout, const char* setname, GSuperLocus& stotal,
       GSeqData* seqdata=NULL, GSeqData* refdata=NULL);

GSeqData* getQryData(int gid, GList<GSeqData>& qdata);
void trackGData(int qcount, GList<GSeqTrack>& gtracks, GStr& fbasename, FILE** ftr, FILE** frs);

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)
#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)

FILE* f_mintr=NULL; //missed ref introns

bool multiexon_only=false;
bool multiexonrefs_only=false;

GHash<GStr> refdescr;
void loadRefDescr(const char* fname);

GList<GStr> qryfiles(false,true,false);

//list of GSeqTrack data, sorted by gseq_id
GList<GSeqTrack> gseqtracks(true,true,true);
GSeqTrack* findGSeqTrack(int gsid);


int cmpGTrackByName(const pointer p1, const pointer p2) {
 return strcmp(((GSeqTrack*)p1)->gseq_name, ((GSeqTrack*)p2)->gseq_name);
}


void show_usage() {
  GMessage("cuffcompare v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION);
  GMessage( "-----------------------------\n");
  GMessage("%s\n", USAGE);
  }

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "XDTMNVGSCKRLhp:c:d:s:i:n:r:o:");
  int e;
  if ((e=args.isError())>0) {
    show_usage();
    GMessage("Invalid argument: %s\n", argv[e]);
    exit(1);
    }
  if (args.getOpt('h')!=NULL){
    show_usage();
    exit(1);
    }
  showContained=(args.getOpt('C')!=NULL);
  debug=(args.getOpt('D')!=NULL);
  tmapFiles=(args.getOpt('T')==NULL);
  multiexon_only=(args.getOpt('M')!=NULL);
  multiexonrefs_only=(args.getOpt('N')!=NULL); 
  perContigStats=(args.getOpt('S')!=NULL);
  checkFasta=(args.getOpt('K')!=NULL);
  gtf_tracking_verbose=((args.getOpt('V')!=NULL) || debug);
  FILE* finlst=NULL;
  GStr s=args.getOpt('i');
  if (!s.is_empty()) {
      if (s=='-')
       finlst=stdin;
      else {
       finlst=fopen(s,"r");
       if (finlst==NULL) GError("Error opening file: %s\n", s.chars());
             }
           }
        int numqryfiles=0;
        if (finlst) {
          GLineReader* lr=new GLineReader(finlst);
          char* l=NULL;
          while ((l=lr->getLine())!=NULL) {
            if (strlen(l)<2 || startsWith(l,"# ") || isspace(*l)) continue;
            if (!fileExists(l)) GError("Error: cannot locate input file: %s\n", l);
            qryfiles.Add(new GStr(l));
            }
          delete lr;
          //if (qryfiles.Count()>10) 
          gtf_tracking_largeScale=true;
          }
         else {
          numqryfiles=args.startNonOpt();
          char *infile=NULL;
          if (numqryfiles>0) {
            while ((infile=args.nextNonOpt())!=NULL) {
              if (!fileExists(infile)) GError("Error: cannot locate input file: %s\n", infile);
              qryfiles.Add(new GStr(infile));
              } //for each argument
            }
          }
        numqryfiles=qryfiles.Count();
        if (numqryfiles==0) {
          show_usage();
          exit(1);
          }
      if (numqryfiles>MAX_QFILES) {
           GMessage("Error: too many input files (limit set to %d at compile time)\n",MAX_QFILES);
           GMessage("(if you need to raise this limit set a new value for\nMAX_QFILES in gtf_tracking.h and recompile)\n");
           exit(0x5000);
           }
  #ifdef HAVE_CONFIG_H
  check_version(PACKAGE_VERSION);
  #endif
  gfasta.init(args.getOpt('s'));
   // determine if -s points to a multi-fasta file or a directory
  s=args.getOpt('c');
  if (!s.is_empty()) scoreThreshold=s.asReal();
  s=args.getOpt('p');
  if (!s.is_empty()) cprefix=Gstrdup(s.chars());
                  else  cprefix=Gstrdup("TCONS");
  s=args.getOpt('e');
  if (!s.is_empty()) exprThreshold=s.asReal();
  s=args.getOpt('d');
  if (!s.is_empty()) {
     tssDist=s.asInt();
     }

  s=args.getOpt('n');
  if (!s.is_empty()) loadRefDescr(s.chars());
  s=args.getOpt('r');
  if (!s.is_empty()) {
    f_ref=fopen(s,"r");
    if (f_ref==NULL) GError("Error opening reference gff: %s\n",s.chars());
    haveRefs=true;
    if (gtf_tracking_verbose) GMessage("Loading reference transcripts..\n");
    read_mRNAs(f_ref, ref_data, &ref_data, true, -1, s.chars(), (multiexonrefs_only || multiexon_only));
    haveRefs=(ref_data.Count()>0);
    reduceRefs=(args.getOpt('R')!=NULL);
    if (gtf_tracking_verbose) GMessage("..reference annotation loaded\n");
    }
  bool discard_redundant=true; //discard redundant input transfrags
  generic_GFF=args.getOpt('G');
  if (generic_GFF) discard_redundant=false; //generic GTF, don't try to discard "redundant" transcripts
  //if a full pathname is given
  //the other common output files will still be created in the current directory:
  // .loci, .tracking, .stats
  GStr outbasename; //include path, if provided
  GStr outprefix; //without path and/or extension
  GStr outstats=args.getOpt('o');
  if (outstats.is_empty() || outstats=="-") {
       outstats="cuffcmp";
       }
  outbasename=outstats;
  GStr outext(getFileExt(outstats.chars()));
  if (outext.is_empty()) {
    outext="stats";
    outstats.append(".stats");
    outbasename=outstats;
    }
    else outext.lower();
  if (outext=="txt" || outext=="out" || outext=="stats" || outext=="summary") {
      outbasename.cut(outbasename.length()-outext.length()-1);
      }

  outprefix=outbasename;
  int di=outprefix.rindex(CHPATHSEP);
  if (di>=0) outprefix.cut(0,di+1);
  
  if (debug) { //create a few more files potentially useful for debugging
        s=outbasename;
        s.append(".missed_introns.gtf");
        f_mintr=fopen(s.chars(),"w");
        if (f_mintr==NULL) GError("Error creating file %s!\n",s.chars());
        /*
        s=outbasename;
        s.append(".noTP_introns.gtf");
        f_nintr=fopen(s.chars(),"w");
        s=outbasename;
        s.append(".wrong_Qintrons.gtf");
        f_qintr=fopen(s.chars(),"w");
        */
        }

  f_out=fopen(outstats, "w");
  if (f_out==NULL) GError("Error creating output file %s!\n", outstats.chars());
  if (gtf_tracking_verbose) GMessage("Prefix for output files: %s\n", outprefix.chars());
  fprintf(f_out, "# Cuffcompare v%s | Command line was:\n#", PACKAGE_VERSION);
  for (int i=0;i<argc-1;i++) 
    fprintf(f_out, "%s ", argv[i]);
  fprintf(f_out, "%s\n#\n", argv[argc-1]);
  //int qfileno=0;
  GList<GSeqData>** qrysdata=NULL;
  FILE** tfiles=NULL;
  FILE** rtfiles=NULL;
  GMALLOC(qrysdata, numqryfiles*sizeof(GList<GSeqData>*));
  if (tmapFiles) {
      GMALLOC(tfiles, numqryfiles*sizeof(FILE*));
      if (haveRefs) {
          GMALLOC(rtfiles, numqryfiles*sizeof(FILE*));
          }
      }
  for (int fi=0;fi<qryfiles.Count();fi++) {
    GStr in_file(qryfiles[fi]->chars());
    GStr infname(getFileName(qryfiles[fi]->chars())); //file name only
    GStr indir(qryfiles[fi]->chars());
    di=indir.rindex(CHPATHSEP);
    if (di>=0) indir.cut(di+1); //directory path for this input file
          else indir=""; //current directory
    
    if (debug || (gtf_tracking_verbose && !gtf_tracking_largeScale))
        GMessage("Processing qfile #%d: %s\n",fi+1, in_file.chars());
    if (in_file=="-") { f_in=stdin; in_file="stdin"; }
      else {
        f_in=fopen(in_file.chars(),"r");
        if (f_in==NULL) 
            GError("Cannot open input file %s!\n",in_file.chars());
        }
    //f_in is the query gff file to process
    
    GStr sbase(indir);
    sbase.append(outprefix);
    sbase.append(".");
    sbase.append(infname);
    if (tmapFiles) {
        //-- we should keep the infname path, otherwise the remaining file names 
        //   may be the same and clobber each other
        s=sbase;
        s.append(".tmap");
        tfiles[fi]=fopen(s.chars(),"w");
        if (tfiles[fi]==NULL)
          GError("Error creating file '%s'!\n",s.chars());
        fprintf(tfiles[fi],"ref_gene_id\tref_id\tclass_code\tcuff_gene_id\tcuff_id\tFMI\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tcov\tlen\tmajor_iso_id\tref_match_len\n");
        if (haveRefs) {
          s=sbase;
          s.append(".refmap");
          rtfiles[fi]=fopen(s.chars(),"w");
          if (rtfiles[fi]==NULL)
             GError("Error creating file '%s'!\n",s.chars());
          fprintf(rtfiles[fi],"ref_gene_id\tref_id\tclass_code\tcuff_id_list\n");
          }
        }

      GList<GSeqData>* pdata=new GList<GSeqData>(true,true,true);
      qrysdata[fi]=pdata;
      if (gtf_tracking_verbose) GMessage("Loading transcripts from %s..\n",in_file.chars());
      read_mRNAs(f_in, *pdata, &ref_data, discard_redundant, fi, in_file.chars(), multiexon_only);
      GSuperLocus gstats;
      GFaSeqGet *faseq=NULL;
      for (int g=0;g<pdata->Count();g++) { //for each seqdata related to a genomic sequence
        int gsid=pdata->Get(g)->get_gseqid();
        GSeqData* refdata=getRefData(gsid, ref_data);//ref data for this contig
        if (!gtf_tracking_largeScale)
          processLoci(*(pdata->Get(g)), refdata, fi);
        GSeqTrack* seqtrack=findGSeqTrack(gsid); //this will add a gseqtrack if it doesn't exist
        // for gsid
        if (refdata!=NULL) {
          seqtrack->rloci_f=&(refdata->loci_f);
          seqtrack->rloci_r=&(refdata->loci_r);
        }
        seqtrack->qdata[fi]=pdata->Get(g);
        //will only gather data into stats if perContig==false
        if (!gtf_tracking_largeScale) reportStats(f_out, getGSeqName(gsid), gstats,
              pdata->Get(g), refdata);
        if (faseq!=NULL) delete faseq;
      } //for each genomic sequence data
      //there could be genomic sequences with no qry transcripts
      //but with reference transcripts
      if (haveRefs && !reduceRefs && !gtf_tracking_largeScale) {
        for (int r=0;r<ref_data.Count();r++) {
          GSeqData* refdata=ref_data[r];
          int gsid=refdata->get_gseqid();
          if (getQryData(gsid, *pdata)==NULL) {
            reportStats(f_out, getGSeqName(gsid), gstats, NULL, refdata);
            }//completely missed all refdata on this contig
        }
      }
      //now report the summary:
      if (!gtf_tracking_largeScale) reportStats(f_out, in_file.chars(), gstats);
      if (f_in!=stdin) fclose(f_in);
      //qfileno++;
  }//for each input file
  if (f_mintr!=NULL) fclose(f_mintr);
  gseqtracks.setSorted(&cmpGTrackByName);
  if (gtf_tracking_verbose) GMessage("Tracking transcripts across %d query files..\n", numqryfiles);
  trackGData(numqryfiles, gseqtracks, outbasename, tfiles, rtfiles);
  fprintf(f_out, "\n Total union super-loci across all input datasets: %d \n", xlocnum);
  if (numqryfiles>1) {
      fprintf(f_out, "  (%d multi-transcript, ~%.1f transcripts per locus)\n",
           total_xloci_alt, ((double)(GXConsensus::count))/xlocnum);
      }
  if (gtf_tracking_verbose) GMessage("Cleaning up..\n");
  GFREE(cprefix);
  // clean up
  for (int i=0;i<numqryfiles;i++) {
    delete qrysdata[i];
    }
  GFREE(qrysdata);
  GFREE(tfiles);
  GFREE(rtfiles);
  gseqtracks.Clear();
  FRCLOSE(f_ref);
  FWCLOSE(f_out);
  if (gtf_tracking_verbose) GMessage("Done.\n");
  ref_data.Clear();
  //getchar();
} //main ends here

void show_exons(FILE* f, GffObj& m) {
  fprintf(f,"(");
  int imax=m.exons.Count()-1;
  for (int i=0;i<=imax;i++) {
    if (i==imax) fprintf(f,"%d-%d)",m.exons[i]->start, m.exons[i]->end);
            else fprintf(f,"%d-%d,",m.exons[i]->start, m.exons[i]->end);
    }
}

bool ichainMatch(GffObj* t, GffObj* r, bool& exonMatch, int fuzz=0) {
  //t's intron chain is considered matching to reference r
  //if r chain is the same or a subset of  t's chain
  exonMatch=false;
  int imax=r->exons.Count()-1;
  int jmax=t->exons.Count()-1;
  if (imax==0 || jmax==0) {   //single-exon mRNAs
     if (imax!=jmax) return false;
     exonMatch=r->exons[0]->coordMatch(t->exons[0],fuzz);
     /*if (exonMatch) return true;
       else return (r->exons[0]->start>=t->exons[0]->start &&
                    r->exons[0]->end<=t->exons[0]->end);*/
     return exonMatch;
     }

  if (r->exons[imax]->start<t->exons[0]->end ||
      t->exons[jmax]->start<r->exons[0]->end ) //intron chains do not overlap at all
          {
           return false;
          }
  //check intron overlaps
  int i=1;
  int j=1;
  bool exmism=false; //any mismatch
  while (i<=imax && j<=jmax) {
     uint rstart=r->exons[i-1]->end;
     uint rend=r->exons[i]->start;
     uint tstart=t->exons[j-1]->end;
     uint tend=t->exons[j]->start;
     if (tend<rstart) { j++; continue; }
     if (rend<tstart) { i++; continue; }
     break; //here we have an intron overlap
     }
  if (i>1 || i>imax || j>jmax) {
      return false; //no intron overlaps found at all
                   //or first intron of ref not overlapping
      }
  //from now on we expect intron matches up to imax
  if (i!=j || imax!=jmax) { exmism=true; if (fuzz==0) return false; }
  for (;i<=imax && j<=jmax;i++,j++) {
    if (abs((int)(r->exons[i-1]->end-t->exons[j-1]->end))>fuzz ||
        abs((int)(r->exons[i]->start-t->exons[j]->start))>fuzz) {
        return false; //just run away
        }
    }
  //if we made it here, we have matching intron chains up to MIN(imax,jmax)
  if (imax!=jmax) {
     exmism=true;
     if (jmax<imax) return false; // qry ichain included in ref ichain
                     else //ref ichain included in qry ichain
                      if (fuzz==0) return false;
     }
  if (exmism) {
          //exonMatch=false; -- it's default
          return true;
          }
  exonMatch = ( abs((int)(r->exons[0]->start-t->exons[0]->start))<=fuzz &&
               abs((int)(r->exons[imax]->end-t->exons[jmax]->end))<=fuzz );
  return true;
}


void compareLoci2R(GList<GLocus>& loci, GList<GSuperLocus>& cmpdata,
                             GList<GLocus>& refloci, int qfidx) {
 cmpdata.Clear();//a new list of superloci will be built
 if (refloci.Count()==0 || loci.Count()==0) return;
 //reset cmpovl and stats
 for (int i=0;i<refloci.Count();i++) refloci[i]->creset();
 //find loci with overlapping refloci
 //and store cmpovl links both ways for ALL loci and refloci on this strand
 for (int l=0;l<loci.Count();l++) {
   GLocus* locus=loci[l];
   locus->creset();
   for (int j=0;j<refloci.Count();j++) {
     //if (refloci[j]->start>locus->end) break;
     if (refloci[j]->start>locus->end) {
         if (refloci[j]->start-locus->end > GFF_MAX_LOCUS) break;
         continue;
         }
     if (locus->start>refloci[j]->end) continue;
     // then we must have overlap here:
     //if (locus->overlap(refloci[j]->start, refloci[j]->end)) {
        locus->cmpovl.Add(refloci[j]);
        refloci[j]->cmpovl.Add(locus);
        //}
     }//for each reflocus
   } //for each locus

 //create corresponding "superloci" from transitive overlapping between loci and ref
 for (int l=0;l<loci.Count();l++) {
  if (loci[l]->v!=0) continue; //skip, already processed
  GSuperLocus* super=new GSuperLocus();
  super->qfidx=qfidx;
  //try to find all other loci connected to this locus loci[l]
  GList<GLocus> lstack(false,false,false);  //traversal stack
  lstack.Push(loci[l]);
  while (lstack.Count()>0) {
      GLocus* locus=lstack.Pop();
      if (locus->v!=0 || locus->cmpovl.Count()==0) continue;
      super->addQlocus(*locus);
      locus->v=1;
      for (int r=0;r<locus->cmpovl.Count();r++) {
        GLocus* rloc=locus->cmpovl[r];
        if (rloc->v==0) {
          super->addRlocus(*rloc);
          rloc->v=1;
          for (int ll=0;ll<rloc->cmpovl.Count();ll++) {
              if (rloc->cmpovl[ll]->v==0) lstack.Push(rloc->cmpovl[ll]);
              }
           }
        } //for each overlapping reflocus
      } //while linking

  if (super->qloci.Count()==0) {
    delete super;
    continue; //try next query loci
    }
  //--here we have a "superlocus" region data on both qry and ref
  // -- analyze mexons matching (base level metrics)
  cmpdata.Add(super);
  //make each ref locus keep track of all superloci containing it
    for (int rl=0;rl<super->rloci.Count();rl++) {
      super->rloci[rl]->superlst->Add(super);
      }
  for (int x=0;x<super->rmexons.Count();x++) {
    super->rbases_all += super->rmexons[x].end-super->rmexons[x].start+1;
    }
  for (int x=0;x<super->qmexons.Count();x++) {
    super->qbases_all += super->qmexons[x].end-super->qmexons[x].start+1;
    }
  int i=0; //locus mexons
  int j=0; //refmexons
  while (i<super->qmexons.Count() && j<super->rmexons.Count()) {
     uint istart=super->qmexons[i].start;
     uint iend=super->qmexons[i].end;
     uint jstart=super->rmexons[j].start;
     uint jend=super->rmexons[j].end;
     if (iend<jstart) { i++; continue; }
     if (jend<istart) { j++; continue; }
     //v--overlap here:
     uint ovlstart = jstart>istart? jstart : istart;
     uint ovlend = iend<jend ? iend : jend;
     uint ovlen=ovlend-ovlstart+1;
     super->baseTP+=ovlen; //qbases_cov
     if (iend<jend) i++;
               else j++;
     } //while mexons ovl search
  /* if (reduceRefs) {
    super->baseFP=super->qbases_all-super->baseTP;
    super->baseFN=super->rbases_all-super->baseTP;
    }
  */
  // -- exon level comparison:
  int* qexovl; //flags for qry exons with ref overlap
  GCALLOC(qexovl,super->quexons.Count()*sizeof(int));
  int* rexovl; //flags for ref exons with qry overlap
  GCALLOC(rexovl,super->ruexons.Count()*sizeof(int));
  for (int i=0;i<super->quexons.Count();i++) {
    uint istart=super->quexons[i].start;
    uint iend=super->quexons[i].end;
    for (int j=0;j<super->ruexons.Count();j++) {
      uint jstart=super->ruexons[j].start;
      uint jend=super->ruexons[j].end;
      if (iend<jstart) break;
      if (jend<istart) continue;
      //--- overlap here between quexons[i] and ruexons[j]
      qexovl[i]++;
      rexovl[j]++;
      if (super->quexons[i].coordMatch(&super->ruexons[j],5)) {
         super->exonATP++;
         if (super->quexons[i].coordMatch(&super->ruexons[j])) {
             super->exonTP++;
             } //exact match
         } //fuzzy match
      } //ref uexon loop
   } //qry uexon loop
  super->m_exons=0; //ref exons with no query overlap
  super->w_exons=0; //qry exons with no ref overlap
  for (int x=0;x<super->quexons.Count();x++)
       if (qexovl[x]==0) super->w_exons++;
  for (int x=0;x<super->ruexons.Count();x++)
       if (rexovl[x]==0) super->m_exons++;
  GFREE(rexovl);
  GFREE(qexovl);

  //-- intron level stats:
  //query:
  int* qinovl=NULL; //flags for qry introns with at least some ref overlap
  int* qtpinovl=NULL; //flags for qry introns with perfect ref overlap
  if (super->qintrons.Count()>0) {
   GCALLOC(qinovl,super->qintrons.Count()*sizeof(int));
   GCALLOC(qtpinovl,super->qintrons.Count()*sizeof(int));
   }
  //-- reference:
  int* rinovl=NULL; //flags for ref introns with qry overlap
  int* rtpinovl=NULL; //ref introns with perfect qry intron overlap
  if (super->rintrons.Count()>0) {
   GCALLOC(rinovl,super->rintrons.Count()*sizeof(int));
   GCALLOC(rtpinovl,super->rintrons.Count()*sizeof(int));
   }
  for (int i=0;i<super->qintrons.Count();i++) {
    uint istart=super->qintrons[i].start;
    uint iend=super->qintrons[i].end;
    for (int j=0;j<super->rintrons.Count();j++) {
      uint jstart=super->rintrons[j].start;
      uint jend=super->rintrons[j].end;
      if (iend<jstart) break;
      if (jend<istart) continue;
      //--- overlap here between qintrons[i] and rintrons[j]
      qinovl[i]++;
      rinovl[j]++;
      if (super->qintrons[i].coordMatch(&super->rintrons[j],5)) {
         super->intronATP++;
         if (super->qintrons[i].coordMatch(&super->rintrons[j])) {
             super->intronTP++;
             qtpinovl[i]++;
             rtpinovl[j]++;
             } //exact match
         } //fuzzy match
      } //ref intron loop
   } //qry intron loop
  super->m_introns=0; //ref introns with no query overlap
  super->w_introns=0; //qry introns with no ref overlap
  for (int x=0;x<super->qintrons.Count();x++) {
       if (qinovl[x]==0) { super->w_introns++;
                 //qry introns with no ref intron overlap AT ALL
                 super->i_qwrong.Add(super->qintrons[x]);
                 }
          else
             if (qtpinovl[x]==0) {
               super->i_qnotp.Add(super->qintrons[x]);
             }
       }
  for (int x=0;x<super->rintrons.Count();x++) {
       if (rinovl[x]==0) { //no intron overlap at all
             super->m_introns++;
             super->i_missed.Add(super->rintrons[x]);
             }
       else if (rtpinovl[x]==0) { //no perfect intron match
            super->i_notp.Add(super->rintrons[x]);
            }
       }
  GFREE(rinovl);
  GFREE(rtpinovl);
  GFREE(qinovl);
  GFREE(qtpinovl);

  // ---- now intron-chain and transcript comparison
  for (int i=0;i<super->qmrnas.Count();i++) {
    uint istart=super->qmrnas[i]->exons.First()->start;
    uint iend=super->qmrnas[i]->exons.Last()->end;
    for (int j=0;j<super->rmrnas.Count();j++) {
      uint jstart=super->rmrnas[j]->exons.First()->start;
      uint jend=super->rmrnas[j]->exons.Last()->end;
      if (iend<jstart) break;
      if (jend<istart) continue;
      //--- overlap here --
      bool exonMatch=false;
      if ((super->qmrnas[i]->udata & 3) > 1) continue; //already counted a ichainTP for this qry
      if (ichainMatch(super->qmrnas[i],super->rmrnas[j],exonMatch, 5)) { //fuzzy match
         GLocus* qlocus=((CTData*)super->qmrnas[i]->uptr)->locus;
         GLocus* rlocus=((CTData*)super->rmrnas[j]->uptr)->locus;
         if (super->qmrnas[i]->exons.Count()>1) {
              super->ichainATP++;
              qlocus->ichainATP++;
              rlocus->ichainATP++;
              }
          if (exonMatch) {
                super->mrnaATP++;
                qlocus->mrnaATP++;
                rlocus->mrnaATP++;
                }
         if (ichainMatch(super->qmrnas[i],super->rmrnas[j],exonMatch)) { //exact match
             if (super->qmrnas[i]->exons.Count()>1) {
                super->qmrnas[i]->udata|=1;
                super->ichainTP++;
                qlocus->ichainTP++;
                rlocus->ichainTP++;
                }
             if (exonMatch) {
                super->qmrnas[i]->udata|=2;
                super->mrnaTP++;
                qlocus->mrnaTP++;
                rlocus->mrnaTP++;
                }
             } //exact match
         } //fuzzy match
      } //ref mrna loop
    } //qry mrna loop
  for (int ql=0;ql<super->qloci.Count();ql++) {
      if (super->qloci[ql]->ichainTP+super->qloci[ql]->mrnaTP >0 )
                 super->locusQTP++;
      if (super->qloci[ql]->ichainATP+super->qloci[ql]->mrnaATP>0)
                 super->locusAQTP++;
      }
  for (int rl=0;rl<super->rloci.Count();rl++) {
      if (super->rloci[rl]->ichainTP+super->rloci[rl]->mrnaTP >0 )
                 super->locusTP++;
      if (super->rloci[rl]->ichainATP+super->rloci[rl]->mrnaATP>0)
                 super->locusATP++;
      }

  }//for each unlinked locus

}

//look for qry data for a specific genomic sequence
GSeqData* getQryData(int gid, GList<GSeqData>& qdata) {
  int qi=-1;
  GSeqData f(gid);
  GSeqData* q=NULL;
  if (qdata.Found(&f,qi))
        q=qdata[qi];
  return q;
}

const char* findDescr(GffObj* gfobj) {
  if (refdescr.Count()==0) return NULL;
  GStr* s=refdescr.Find(gfobj->getID());
  if (s==NULL) {
       s=refdescr.Find(gfobj->getGeneName());
       if (s==NULL) s=refdescr.Find(gfobj->getGeneID());
       }
  if (s!=NULL)
     return s->chars();
  return NULL;
}

const char* getGeneID(GffObj* gfobj) {
 //returns anything that might resemble a gene identifier for this transcript
 //or, if everything fails, returns the transcript ID
 const char* s=gfobj->getGeneName();
 if (s) return s;
 if ((s=gfobj->getGeneID())!=NULL) return s;
 if ((s=gfobj->getAttr("Name"))!=NULL) return s;
 return gfobj->getID();
}

const char* getGeneID(GffObj& gfobj) {
 return getGeneID(&gfobj);
}

void writeLoci(FILE* f, GList<GLocus> & loci) {
 for (int l=0;l<loci.Count();l++) {
   GLocus& loc=*(loci[l]);
   fprintf(f,"%s\t%s[%c]%d-%d\t", loc.mrna_maxcov->getID(),
       loc.mrna_maxcov->getGSeqName(),
           loc.mrna_maxcov->strand, loc.start,loc.end);
   //now print all transcripts in this locus, comma delimited
   int printfd=0;
   for (int i=0;i<loc.mrnas.Count();i++) {
      if (loc.mrnas[i]==loc.mrna_maxcov) continue;
      if (printfd==0) fprintf(f,"%s",loc.mrnas[i]->getID());
          else fprintf(f,",%s",loc.mrnas[i]->getID());
      printfd++;
      }
   const char* rdescr=findDescr(loc.mrna_maxcov);
   if (rdescr==NULL)  fprintf(f,"\t\n");
                 else fprintf(f,"\t%s\n",rdescr);
   }
}

void printXQ1(FILE* f, int qidx, GList<GLocus>& qloci) {
  int printfd=0;
  //print
  for (int i=0;i<qloci.Count();i++) {
     if (qloci[i]->qfidx!=qidx) continue;
      for (int j=0;j<qloci[i]->mrnas.Count();j++) {
        if (printfd==0) fprintf(f,"%s",qloci[i]->mrnas[j]->getID());
            else fprintf(f,",%s",qloci[i]->mrnas[j]->getID());
        printfd++;
        }
      }
  if (printfd==0) fprintf(f,"-");
 }

void numXLoci(GList<GXLocus>& xloci, int& last_id) {
  for (int l=0;l<xloci.Count();l++) {
    if (xloci[l]->qloci.Count()==0) continue; //we never print ref-only xloci
    last_id++;
    xloci[l]->id=last_id;
    }
}


class GProtCl {
 public:
   GList<GXConsensus> protcl;
   GProtCl(GXConsensus* c=NULL):protcl(true,false,false) {
    if (c!=NULL)
       protcl.Add(c);
    }
   bool add_Pcons(GXConsensus* c) {
    if (c==NULL || c->aalen==0) return false;
    if (protcl.Count()==0) {
        protcl.Add(c);
        return true;
        }
    if (protcl[0]->aalen!=c->aalen) return false;
    if (strcmp(protcl[0]->aa,c->aa)!=0) return false;
    protcl.Add(c);
    return true;
    }

   void addMerge(GProtCl& pcl, GXConsensus* pclnk) {
    for (int i=0;i<pcl.protcl.Count();i++) {
      if (pcl.protcl[i]!=pclnk) {
          protcl.Add(pcl.protcl[i]);
          }
      }
    }

   int aalen() {
    if (protcl.Count()==0) return 0;
    return protcl[0]->aalen;
   }
   bool operator==(GProtCl& cl) {
    return this==&cl;
    }
   bool operator>(GProtCl& cl) {
     return (this>&cl);
    }
   bool operator<(GProtCl& cl) {
    return (this<&cl);
    }
};

class GTssCl:public GSeg { //experiment cluster of ref loci (isoforms)
 public:
   uint fstart; //lowest coordinate of the first exon
   uint fend; //highest coordinate of the first exon
   GList<GXConsensus> tsscl;
   GTssCl(GXConsensus* c=NULL):tsscl(true,false,false) {
     start=0;
     end=0;
     fstart=0;
     fend=0;
     if (c!=NULL) addFirst(c);
     }

   void addFirst(GXConsensus* c) {
     tsscl.Add(c);
     start=c->start;
     end=c->end;
     GffExon* fexon=(c->tcons->strand=='-') ? c->tcons->exons.Last() :
                                             c->tcons->exons.First();
     fstart=fexon->start;
     fend=fexon->end;
     }
   bool add_Xcons(GXConsensus* c) {
     if (tsscl.Count()==0) {
            addFirst(c);
            return true;
            }
     //check if it can be added to existing xconsensi
     uint nfend=0;
     uint nfstart=0;
     /*
     if (tsscl.Get(0)->tcons->getGeneID()!=NULL && 
             c->tcons->getGeneID()!=NULL && 
            strcmp(tsscl.Get(0)->tcons->getGeneID(), c->tcons->getGeneID()))
        //don't tss cluster if they don't have the same GeneID (?)
        //FIXME: we might not want this if input files are not from Cufflinks
        //       and they could simply lack proper GeneID
          return false;
      */
     if (c->tcons->strand=='-') {
        //no, the first exons don't have to overlap
        //if (!c->tcons->exons.Last()->overlap(fstart,fend)) return false;
        nfstart=c->tcons->exons.Last()->start;
        nfend=c->tcons->exons.Last()->end;
        //proximity check for the transcript start:
        if (nfend>fend+tssDist || fend>nfend+tssDist)
                return false;
        }
      else {
        //if (!c->tcons->exons.First()->overlap(fstart,fend)) return false;
        nfstart=c->tcons->exons.First()->start;
        nfend=c->tcons->exons.First()->end;
        if (nfstart>fstart+tssDist || fstart>nfstart+tssDist)
            return false;
        }
     // -- if we are here, we can add to tss cluster

     tsscl.Add(c);
     if (fstart>nfstart) fstart=nfstart;
     if (fend<nfend) fend=nfend;
     if (start>c->start) start=c->start;
     if (end<c->end) end=c->end;
     return true;
     }

   void addMerge(GTssCl& cl, GXConsensus* clnk) {
     for (int i=0;i<cl.tsscl.Count();i++) {
         if (cl.tsscl[i]==clnk) continue;
         tsscl.Add(cl.tsscl[i]);
         }
     if (fstart>cl.fstart) fstart=cl.fstart;
     if (fend<cl.fend) fend=cl.fend;
     if (start>cl.start) start=cl.start;
     if (end<cl.end) end=cl.end;
     }
};
/*
class IntArray { //two dimensional int array
    int* mem;
    int xsize;
    int ysize;
  public:
   IntArray(int xlen, int ylen) {
     xsize=xlen;
     ysize=ylen;
     GMALLOC(mem, xsize*ysize*sizeof(int));
     }
   ~IntArray() {
     GFREE(mem);
     }
   int& data(int x, int y) {
    return mem[y*xsize+x];
    }
};

int aa_diff(GXConsensus* c1, GXConsensus* c2) {
 int diflen=abs(c1->aalen-c2->aalen);
 if (diflen>=6) return diflen;
 //obvious case: same CDS
 if (diflen==0 && strcmp(c1->aa, c2->aa)==0) return 0;
 //simple edit distance calculation
 IntArray dist(c1->aalen+1, c2->aalen+1);
 for (int i=0;i<=c1->aalen;i++) {
     dist.data(i,0) = i;
     }
 for (int j = 0; j <= c2->aalen; j++) {
     dist.data(0,j) = j;
     }
 for (int i = 1; i <= c1->aalen; i++)
     for (int j = 1; j <= c2->aalen; j++) {
         dist.data(i,j) = GMIN3( dist.data(i-1,j)+1,
             dist.data(i,j-1)+1,
                 dist.data(i-1,j-1)+((c1->aa[i-1] == c2->aa[j-1]) ? 0 : 1) );
         }
 int r=dist.data(c1->aalen,c2->aalen);
 return r;
}
*/
void printConsGTF(FILE* fc, GXConsensus* xc, int xlocnum) {
 for (int i=0;i<xc->tcons->exons.Count();i++) {
   fprintf(fc,
     "%s\t%s\texon\t%d\t%d\t.\t%c\t.\tgene_id \"XLOC_%06d\"; transcript_id \"%s_%08d\"; exon_number \"%d\";",
     xc->tcons->getGSeqName(),xc->tcons->getTrackName(),xc->tcons->exons[i]->start, xc->tcons->exons[i]->end, xc->tcons->strand,
       xlocnum, cprefix, xc->id, i+1);
    //if (i==0) {
  const char* gene_name=NULL;
  if (xc->ref) {
     gene_name=xc->ref->getGeneName();
     if (gene_name==NULL) gene_name=xc->ref->getGeneID();
     if (gene_name) {
       fprintf (fc, " gene_name \"%s\";", gene_name);
       }
   }
   if (!haveRefs) {
      if (gene_name==NULL && xc->tcons->getGeneName())
        fprintf (fc, " gene_name \"%s\";", xc->tcons->getGeneName());
      char* s=xc->tcons->getAttr("nearest_ref", true);
      if (s) fprintf(fc, " nearest_ref \"%s\";",s);
      s=xc->tcons->getAttr("class_code", true);
      if (s) fprintf(fc, " class_code \"%s\";", s);
      }
   fprintf(fc, " oId \"%s\";",xc->tcons->getID());
   if (xc->contained) {
     fprintf(fc, " contained_in \"%s_%08d\";", cprefix, xc->contained->id);
     }
   if (haveRefs) {
     if (xc->ref!=NULL)
       fprintf(fc, " nearest_ref \"%s\";",xc->ref->getID());
     fprintf(fc, " class_code \"%c\";",xc->refcode ? xc->refcode : '.');
     }
   if (xc->tss_id>0) fprintf(fc, " tss_id \"TSS%d\";",xc->tss_id);
   if (xc->p_id>0) fprintf(fc, " p_id \"P%d\";",xc->p_id);
   //      }
   fprintf(fc,"\n");
   }
}

void tssCluster(GXLocus& xloc) 
{
    GList<GTssCl> xpcls(true,true,false);
    for (int i=0;i<xloc.tcons.Count();i++) 
    {
        GXConsensus* c=xloc.tcons[i];
        //if (c->tcons->exons.Count()<2) continue;  //should we skip single-exon transcripts ??
        GArray<int> mrgloci(true);
        int lfound=0;
        for (int l=0;l<xpcls.Count();l++) 
        {
            if (xpcls[l]->end<c->tcons->exons.First()->start) continue;
            if (xpcls[l]->start>c->tcons->exons.Last()->end) break;
            if (xpcls[l]->add_Xcons(c)) 
            {
                lfound++;
                mrgloci.Add(l);
                
            }
            
        } // for each xpcluster
        if (lfound==0) 
        {
            //create a xpcl with only this xconsensus
            xpcls.Add(new GTssCl(c));
            
        }
        else if (lfound>1) 
        {
            for (int l=1;l<lfound;l++) 
            {
                int mlidx=mrgloci[l]-l+1;
                xpcls[mrgloci[0]]->addMerge(*xpcls[mlidx], c);
                xpcls.Delete(mlidx);
            }
        }
        
    }//for each xconsensus in this xlocus
    for (int l=0;l<xpcls.Count();l++) 
    {
        //if (xpcls[l]->tsscl.Count()<2) continue;
        tsscl_num++;
        for (int i=0;i<xpcls[l]->tsscl.Count();i++)
            xpcls[l]->tsscl[i]->tss_id=tsscl_num;
        //processTssCl(xcds_num, xpcls[l], faseq);
    }
}

void protCluster(GXLocus& xloc, GFaSeqGet *faseq) {
  if (!faseq)
    return;
  GList<GProtCl> xpcls(true,true,false);
  for (int i=0;i<xloc.tcons.Count();i++) {
    GXConsensus* c=xloc.tcons[i];
    if (c->ref==NULL || c->ref->CDstart==0) continue;  //no ref or CDS available
    if (c->refcode!='=') continue;
    //get the CDS translation here
    if (c->aa==NULL) {
       c->aa=c->ref->getSplicedTr(faseq, true, &c->aalen);
       if (c->aalen>0 && c->aa[c->aalen-1]=='.') {
             //discard the final stop codon
             c->aalen--;
             c->aa[c->aalen]=0;
             }
       }
    GArray<int> mrgloci(true);
    int lfound=0;
    for (int l=0;l<xpcls.Count();l++) {
        if (xpcls[l]->aalen()!=c->aalen) continue;
        if (xpcls[l]->add_Pcons(c)) {
            lfound++;
            mrgloci.Add(l);
            }
        } // for each xpcluster
    if (lfound==0) {
        //create a xpcl with only this xconsensus
        xpcls.Add(new GProtCl(c));
        }
      else if (lfound>1) {
        for (int l=1;l<lfound;l++) {
              int mlidx=mrgloci[l]-l+1;
              xpcls[mrgloci[0]]->addMerge(*xpcls[mlidx], c);
              xpcls.Delete(mlidx);
              }
        }
    }//for each xconsensus in this xlocus
 for (int l=0;l<xpcls.Count();l++) {
      protcl_num++;
      for (int i=0;i<xpcls[l]->protcl.Count();i++)
           xpcls[l]->protcl[i]->p_id=protcl_num;
      }
 for (int i=0;i<xloc.tcons.Count();i++) {
   GXConsensus* c=xloc.tcons[i];
   if (c->aa!=NULL) { GFREE(c->aa); }
   }
}

void printXLoci(FILE* f, FILE* fc, int qcount, GList<GXLocus>& xloci, GFaSeqGet *faseq) {
  for (int l=0;l<xloci.Count();l++) {
    if (xloci[l]->qloci.Count()==0) continue;
    GXLocus& xloc=*(xloci[l]);
    xloc.checkContainment(); 
    tssCluster(xloc);//cluster and assign tss_id and cds_id to each xconsensus in xloc
    protCluster(xloc,faseq);
    for (int c=0;c<xloc.tcons.Count();c++) {
       if (showContained || xloc.tcons[c]->contained==NULL)
         printConsGTF(fc,xloc.tcons[c],xloc.id);
       }
    fprintf(f,"XLOC_%06d\t%s[%c]%d-%d\t", xloc.id,
        xloc.qloci[0]->mrna_maxcov->getGSeqName(),
            xloc.strand, xloc.start,xloc.end);
    //now print all transcripts in this locus, comma delimited
    //first, ref loci, if any
    int printfd=0;
    if (xloc.rloci.Count()>0) {
       for (int i=0;i<xloc.rloci.Count();i++) {
          for (int j=0;j<xloc.rloci[i]->mrnas.Count();j++) {
            if (printfd==0) fprintf(f,"%s|%s",getGeneID(xloc.rloci[i]->mrnas[j]),
                                                    xloc.rloci[i]->mrnas[j]->getID());
                else fprintf(f,",%s|%s",getGeneID(xloc.rloci[i]->mrnas[j]),
                                                    xloc.rloci[i]->mrnas[j]->getID());
            printfd++;
            }
          }
       }
    else {
      fprintf(f,"-");
      }
   //second, all the cufflinks transcripts
    for (int qi=0;qi<qcount;qi++) {
       fprintf(f,"\t");
       printXQ1(f,qi,xloc.qloci);
       }
    fprintf(f,"\n");
    }
}

void writeIntron(FILE* f, char strand, GFaSeqGet* faseq, GSeg& iseg,
                GList<GffObj>& mrnas, bool wrong=false) {
//find a ref mrna having this intron
  GffObj* rm=NULL;
  for (int i=0;i<mrnas.Count();i++) {
   GffObj* m=mrnas[i];
   if (m->start>iseg.end) break;
   if (m->end<iseg.start) continue;
   //intron coords overlaps mrna region
   for (int j=1;j<m->exons.Count();j++) {
      if (iseg.start==m->exons[j-1]->end+1 &&
            iseg.end==m->exons[j]->start-1) { rm=m; break; } //match found
      }//for each intron
   if (rm!=NULL) break;
   } //for each ref mrna in this locus
 if (rm==NULL) GError("Error: couldn't find ref mrna for intron %d-%d! (BUG)\n",
                         iseg.start,iseg.end);
 int ilen=iseg.end-iseg.start+1;
 fprintf(f,"%s\t%s\tintron\t%d\t%d\t.\t%c\t.\t",
            rm->getGSeqName(),rm->getTrackName(),iseg.start,iseg.end,strand);
 if (faseq!=NULL) {
   const char* gseq=faseq->subseq(iseg.start, ilen);
   char* cseq=Gstrdup(gseq, gseq+ilen-1);
   if (strand=='-') reverseComplement(cseq, ilen);
   fprintf(f,"spl=\"%c%c..%c%c\"; ", toupper(cseq[0]),toupper(cseq[1]),
           toupper(cseq[ilen-2]),toupper(cseq[ilen-1]));
   GFREE(cseq);
   }
  fprintf(f,"transcript_id \"%s\";", rm->getID());
  if (wrong) fprintf(f," noOvl=1;");
  fprintf(f,"\n");

}

void reportMIntrons(FILE* fm, FILE* fn, FILE* fq, char strand,
            GList<GSuperLocus>& cmpdata) {
  if (fm==NULL) return;
  for (int l=0;l<cmpdata.Count();l++) {
    GSuperLocus *sl=cmpdata[l];
    //cache the whole locus sequence if possible
    //write these introns and their splice sites into the file
    for (int i=0;i<sl->i_missed.Count();i++)
      writeIntron(fm, strand, NULL, sl->i_missed[i], sl->rmrnas);
    if (fn!=NULL) {
        for (int i=0;i<sl->i_notp.Count();i++)
          writeIntron(fn, strand, NULL, sl->i_notp[i], sl->rmrnas);
        }
    if (fq!=NULL) {
     for (int i=0;i<sl->i_qwrong.Count();i++) {
       writeIntron(fq, strand, NULL, sl->i_qwrong[i], sl->qmrnas, true);
       }
     for (int i=0;i<sl->i_qnotp.Count();i++) {
       writeIntron(fq, strand, NULL, sl->i_qnotp[i], sl->qmrnas);
       }
     }
    }
}


void processLoci(GSeqData& seqdata, GSeqData* refdata, int qfidx) {
    //GList<GSeqLoci>& glstloci, GList<GSeqCmpRegs>& cmpdata)

  if (refdata!=NULL) {
     //if (gtf_tracking_verbose) GMessage(" ..comparing to reference loci..\n") ;
     compareLoci2R(seqdata.loci_f, seqdata.gstats_f, refdata->loci_f, qfidx);
     compareLoci2R(seqdata.loci_r, seqdata.gstats_r, refdata->loci_r, qfidx);
     // -- report
     
     if (f_mintr!=NULL) {
       GMessage(" ..reporting missed ref introns..\n");
       //reportIntrons(f_mintr, f_nintr, f_qintr, faseq, '+', seqdata.gstats_f);
       //reportIntrons(f_mintr, f_nintr, f_qintr, faseq, '-', seqdata.gstats_r);
       reportMIntrons(f_mintr, NULL, NULL, '+', seqdata.gstats_f);
       reportMIntrons(f_mintr, NULL, NULL, '-', seqdata.gstats_r);
       }
     }
}

//adjust stats for a list of unoverlapped (completely missed) ref loci
void collectRLocData(GSuperLocus& stats, GLocus& loc) {
stats.total_rmrnas+=loc.mrnas.Count();
stats.total_rexons+=loc.uexons.Count();
stats.total_rintrons+=loc.introns.Count();
stats.total_rmexons+=loc.mexons.Count();
stats.total_richains+=loc.ichains;
stats.m_exons+=loc.uexons.Count();
stats.m_introns+=loc.introns.Count();
stats.total_rloci++;
for (int e=0;e<loc.mexons.Count();e++) {
   stats.rbases_all+=loc.mexons[e].end-loc.mexons[e].start+1;
   }
}

void collectRData(GSuperLocus& stats, GList<GLocus>& loci) {
 for (int l=0;l<loci.Count();l++)
      collectRLocData(stats,*loci[l]);
}

//adjust stats for a list of unoverlapped (completely "wrong" or novel) qry loci
void collectQLocData(GSuperLocus& stats, GLocus& loc) {
 stats.total_qmrnas+=loc.mrnas.Count();
 stats.total_qexons+=loc.uexons.Count();
 stats.total_qmexons+=loc.mexons.Count();
 stats.total_qintrons+=loc.introns.Count();
 stats.total_qichains+=loc.ichains;
 stats.total_qloci++;
 if (loc.ichains>0 && loc.mrnas.Count()>1)
    stats.total_qloci_alt++;
 stats.w_exons+=loc.uexons.Count();
 stats.w_introns+=loc.introns.Count();
 for (int e=0;e<loc.mexons.Count();e++) {
   stats.qbases_all+=loc.mexons[e].end-loc.mexons[e].start+1;
   }
}

void collectQData(GSuperLocus& stats, GList<GLocus>& loci, GList<GLocus>& nloci) {
 for (int l=0;l<loci.Count();l++) {
     //this is called when no refdata is given, so all these loci are nloci
     nloci.Add(loci[l]);
     collectQLocData(stats,*loci[l]);
     }
}

void collectQNOvl(GSuperLocus& stats, GList<GLocus>& loci, GList<GLocus>& nloci) {
  for (int l=0;l<loci.Count();l++) {
    if (loci[l]->cmpovl.Count()==0) {//locus with no ref loci overlaps
      stats.w_loci++; //novel/wrong loci
      nloci.Add(loci[l]);
      collectQLocData(stats,*loci[l]);
      }
  }
}

void collectQU(GSuperLocus& stats, GList<GLocus>& nloci) {
  for (int l=0;l<nloci.Count();l++) {
    stats.w_loci++; //novel/wrong loci
    collectQLocData(stats, *nloci[l]);
    }
}

void printLocus(FILE* f, GLocus& loc, const char* gseqname) {
  fprintf(f, "## Locus %s:%d-%d\n",gseqname, loc.start, loc.end);
  for (int m=0;m<loc.mrnas.Count();m++) {
    loc.mrnas[m]->printGtf(f);
    }
}

void collectRNOvl(GSuperLocus& stats, GList<GLocus>& loci) { //, const char* gseqname) {
  for (int l=0;l<loci.Count();l++) {
    if (loci[l]->cmpovl.Count()==0) {
      stats.m_loci++; //missed ref loci
      //if (f_mloci!=NULL)
      //      printLocus(f_mloci,*loci[l], gseqname);
      collectRLocData(stats,*loci[l]);
      }
  }
}


void collectCmpData(GSuperLocus& stats, GList<GSuperLocus>& cmpdata) { //, const char* gseqname) {
 for (int c=0;c<cmpdata.Count();c++) {
   stats.addStats(*cmpdata[c]);
   /*
   if (f_nloci!=NULL && cmpdata[c]->locusTP==0 && cmpdata[c]->rloci.Count()>0) {       
      fprintf(f_nloci, "# Superlocus %s:%d-%d\n",gseqname, cmpdata[c]->start, cmpdata[c]->end);
      for (int l=0;l<cmpdata[c]->rloci.Count();l++) {
         printLocus(f_nloci,*cmpdata[c]->rloci[l], gseqname);
         }
      }
   */   
   }
}

void collectStats(GSuperLocus& stats, GSeqData* seqdata, GSeqData* refdata) {
 //collect all stats for a single genomic sequence into stats
 if (seqdata==NULL) {
   if (reduceRefs || refdata==NULL) return;
   //special case with completely missed all refs on a contig/chromosome
   collectRData(stats, refdata->loci_f);
   collectRData(stats, refdata->loci_r);
   return;
   }
 if (refdata==NULL) {//reference data missing on this contig
   collectQData(stats, seqdata->loci_f, seqdata->nloci_f);
   collectQData(stats, seqdata->loci_r, seqdata->nloci_r);
   collectQU(stats, seqdata->nloci_u);
   return;
   }

 /*stats.total_qloci+=seqdata->loci_f.Count();
 stats.total_qloci+=seqdata->loci_r.Count();
 if (reduceRefs) { //only collect ref loci from superloci

    }
   else {
    stats.total_rloci+=refdata->loci_f.Count();
    stats.total_rloci+=refdata->loci_r.Count();
    }
 */
 //collect data for overlapping superloci (already in seqdata->gstats_f/_r)
 //char* gseqname=getGSeqName(seqdata->gseq_id);
 collectCmpData(stats, seqdata->gstats_f);
 collectCmpData(stats, seqdata->gstats_r);
 //for non-overlapping qry loci, always add them as false positives FP
 collectQNOvl(stats, seqdata->loci_f, seqdata->nloci_f);
 collectQNOvl(stats, seqdata->loci_r, seqdata->nloci_r);
 collectQU(stats, seqdata->nloci_u);
 if (!reduceRefs) { //find ref loci with empty cmpovl and add them
  collectRNOvl(stats, refdata->loci_f);
  collectRNOvl(stats, refdata->loci_r);
  }
}

void reportStats(FILE* fout, const char* setname, GSuperLocus& stotal,
                          GSeqData* seqdata, GSeqData* refdata) {
  GSuperLocus stats;
  bool finalSummary=(seqdata==NULL && refdata==NULL);
  GSuperLocus *ps=(finalSummary ? &stotal : &stats );
  if (!finalSummary) { //collecting contig stats
    //gather statistics for all loci/superloci here
    collectStats(stats, seqdata, refdata);
    stotal.addStats(stats);
    if (!perContigStats) return;
    }
  ps->calcF();
  if (seqdata!=NULL) fprintf(fout, "#> Genomic sequence: %s \n", setname);
                else fprintf(fout, "\n#= Summary for dataset: %s :\n", setname);

  fprintf(fout,   "#     Query mRNAs : %7d in %7d loci  (%d multi-exon transcripts)\n",
          ps->total_qmrnas, ps->total_qloci, ps->total_qichains);
  fprintf(fout, "#            (%d multi-transcript loci, ~%.1f transcripts per locus)\n",
          ps->total_qloci_alt, ((double)ps->total_qmrnas/ps->total_qloci));

  if (haveRefs) {
    fprintf(fout, "# Reference mRNAs : %7d in %7d loci  (%d multi-exon)\n",
            ps->total_rmrnas, ps->total_rloci, ps->total_richains);
    if (ps->baseTP+ps->baseFP==0 || ps->baseTP+ps->baseFN==0) return;
    fprintf(fout, "# Corresponding super-loci:        %7d\n",ps->total_superloci);

    /*if (seqdata!=NULL) {
      fprintf(fout, "          ( %d/%d on forward/reverse strand)\n",
             seqdata->gstats_f.Count(),seqdata->gstats_r.Count());
       }*/
    fprintf(fout, "#--------------------|   Sn   |  Sp   |  fSn |  fSp  \n");
    double sp=(100.0*(double)ps->baseTP)/(ps->baseTP+ps->baseFP);
    double sn=(100.0*(double)ps->baseTP)/(ps->baseTP+ps->baseFN);
    fprintf(fout, "        Base level: \t%5.1f\t%5.1f\t  - \t  - \n",sn, sp);
    sp=(100.0*(double)ps->exonTP)/(ps->exonTP+ps->exonFP);
    sn=(100.0*(double)ps->exonTP)/(ps->exonTP+ps->exonFN);
    double fsp=(100.0*(double)ps->exonATP)/(ps->exonATP+ps->exonAFP);
    double fsn=(100.0*(double)ps->exonATP)/(ps->exonATP+ps->exonAFN);
    if (fsp>100.0) fsp=100.0;
    if (fsn>100.0) fsn=100.0;
    fprintf(fout, "        Exon level: \t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sn, sp, fsn, fsp);
  if (ps->total_rintrons>0) {
    //intron level
    sp=(100.0*(double)ps->intronTP)/(ps->intronTP+ps->intronFP);
    sn=(100.0*(double)ps->intronTP)/(ps->intronTP+ps->intronFN);
    fsp=(100.0*(double)ps->intronATP)/(ps->intronATP+ps->intronAFP);
    fsn=(100.0*(double)ps->intronATP)/(ps->intronATP+ps->intronAFN);
    if (fsp>100.0) fsp=100.0;
    if (fsn>100.0) fsn=100.0;
    fprintf(fout, "      Intron level: \t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sn, sp, fsn, fsp);
    //intron chains:
    sp=(100.0*(double)ps->ichainTP)/(ps->ichainTP+ps->ichainFP);
    sn=(100.0*(double)ps->ichainTP)/(ps->ichainTP+ps->ichainFN);
    if (sp>100.0) sp=100.0;
    if (sn>100.0) sn=100.0;
    fsp=(100.0*(double)ps->ichainATP)/(ps->ichainATP+ps->ichainAFP);
    fsn=(100.0*(double)ps->ichainATP)/(ps->ichainATP+ps->ichainAFN);
    if (fsp>100.0) fsp=100.0;
    if (fsn>100.0) fsn=100.0;
    fprintf(fout, "Intron chain level: \t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sn, sp, fsn, fsp);
    }
  else {
    fprintf(fout, "      Intron level: \t  -  \t  -  \t  -  \t  -  \n");
    fprintf(fout, "Intron chain level: \t  -  \t  -  \t  -  \t  -  \n");
    }
    sp=(100.0*(double)ps->mrnaTP)/(ps->mrnaTP+ps->mrnaFP);
    sn=(100.0*(double)ps->mrnaTP)/(ps->mrnaTP+ps->mrnaFN);
    fsp=(100.0*(double)ps->mrnaATP)/(ps->mrnaATP+ps->mrnaAFP);
    fsn=(100.0*(double)ps->mrnaATP)/(ps->mrnaATP+ps->mrnaAFN);
    if (fsp>100.0) fsp=100.0;
    if (fsn>100.0) fsn=100.0;
    fprintf(fout, "  Transcript level: \t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sn, sp, fsn, fsp);
    //sp=(100.0*(double)ps->locusTP)/(ps->locusTP+ps->locusFP);
    sp=(100.0*(double)ps->locusQTP)/ps->total_qloci;
    sn=(100.0*(double)ps->locusTP)/ps->total_rloci;  //(ps->locusTP+ps->locusFN);
    fsp=(100.0*(double)ps->locusAQTP)/ps->total_qloci; //(ps->locusATP+ps->locusAFP);
    fsn=(100.0*(double)ps->locusATP)/ps->total_rloci; //(ps->locusATP+ps->locusAFN);
    fprintf(fout, "       Locus level: \t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sn, sp, fsn, fsp);
    //fprintf(fout, "                   (locus TP=%d, total ref loci=%d)\n",ps->locusTP, ps->total_rloci);
    fprintf(fout, "\nMatching intron chains: %7d\n",ps->ichainTP);
    fprintf(fout, "         Matching loci: %7d\n",ps->locusTP);
    fprintf(fout, "\n");
    sn=(100.0*(double)ps->m_exons)/(ps->total_rexons);
    fprintf(fout, "          Missed exons: %7d/%d\t(%5.1f%%)\n",ps->m_exons, ps->total_rexons, sn);
    sn=(100.0*(double)ps->w_exons)/(ps->total_qexons);
    fprintf(fout, "           Novel exons: %7d/%d\t(%5.1f%%)\n",ps->w_exons, ps->total_qexons,sn);
    if (ps->total_rintrons>0) {
    sn=(100.0*(double)ps->m_introns)/(ps->total_rintrons);
    fprintf(fout, "        Missed introns: %7d/%d\t(%5.1f%%)\n",ps->m_introns, ps->total_rintrons, sn);
    }
    if (ps->total_qintrons>0) {
    sn=(100.0*(double)ps->w_introns)/(ps->total_qintrons);
    fprintf(fout, "         Novel introns: %7d/%d\t(%5.1f%%)\n",ps->w_introns, ps->total_qintrons,sn);
    }
    if (ps->total_rloci>0) {
    sn=(100.0*(double)ps->m_loci)/(ps->total_rloci);
    fprintf(fout, "           Missed loci: %7d/%d\t(%5.1f%%)\n",ps->m_loci, ps->total_rloci, sn);
    }
    if (ps->total_qloci>0) {
    sn=(100.0*(double)ps->w_loci)/(ps->total_qloci);
    fprintf(fout, "            Novel loci: %7d/%d\t(%5.1f%%)\n",ps->w_loci, ps->total_qloci,sn);
    }

  }
}

int inbuf_len=1024; //starting inbuf capacity
char* inbuf=NULL; // incoming buffer for sequence lines.

void loadRefDescr(const char* fname) {
  if (inbuf==NULL)  { GMALLOC(inbuf, inbuf_len); }
  FILE *f=fopen(fname, "rb");
  if (f==NULL) GError("Error opening exon file: %s\n",fname);
  char* line;
  int llen=0;
  off_t fpos;
  while ((line=fgetline(inbuf, inbuf_len, f, &fpos, &llen))!=NULL) {
   if (strlen(line)<=2) continue;
   int idlen=strcspn(line,"\t ");
   char* p=line+idlen;
   if (idlen<llen && idlen>0) {
     *p=0;
      p++;
      refdescr.Add(line, new GStr(p));
      }
  }
}

GSeqTrack* findGSeqTrack(int gsid) {
  GSeqTrack f(gsid);
  int fidx=-1;
  if (gseqtracks.Found(&f,fidx))
     return gseqtracks[fidx];
  fidx=gseqtracks.Add(new GSeqTrack(gsid));
  return gseqtracks[fidx];
}



GffObj* findRefMatch(GffObj& m, GLocus& rloc, int& ovlen) {
 ovlen=0;
 CTData* mdata=((CTData*)m.uptr);
 if (mdata->eqref!=NULL && ((CTData*)(mdata->eqref->uptr))->locus==&rloc) {
	  mdata->eqref=mdata->ovls.First()->mrna; //this should be unnecessary
      //check it?
      return mdata->ovls.First()->mrna;
      }
 //if (rloc==NULL|| m==NULL) return NULL;
 GffObj* ret=NULL;
 for (int r=0;r<rloc.mrnas.Count();r++) {
    int olen=0;
	if (tMatch(m, *(rloc.mrnas[r]),olen, true)) { //return rloc->mrnas[r];
	  /*
	  if (ovlen<olen) {
		  ovlen=olen;
		  ret=rloc.mrnas[r]; //keep the longest matching ref
			 //but this is unnecessary, there can be only one matching ref
			 // because duplicate refs were discarded (unless -G was used)
		  }
	  */
	  mdata->addOvl('=',rloc.mrnas[r], olen);
	  ret=mdata->ovls.First()->mrna;
      //this must be called only for the head of an equivalency chain
      CTData* rdata=(CTData*)rloc.mrnas[r]->uptr;
      rdata->addOvl('=',&m,olen);
      //if (rdata->eqnext==NULL) rdata->eqnext=&m;
      }
    }
 if (ret!=NULL)
   mdata->eqref=ret;
 return ret;
 }


void addXCons(GXLocus* xloc, GffObj* ref, char ovlcode, GffObj* tcons, CEqList* ts) {
 GXConsensus* c=new GXConsensus(tcons, ts, ref, ovlcode);
 //xloc->tcons.Add(c);
 //this will also check c against the other tcons for containment:
 xloc->addXCons(c); 
}


const uint pre_mrna_threshold = 100;

char getOvlCode(GffObj& m, GffObj& r, int& ovlen) {
  ovlen=0;
  if (!m.overlap(r.start,r.end)) return 'u';
  int jmax=r.exons.Count()-1;
  
  if (m.exons.Count()==1) { //single-exon transfrag
     GSeg mseg(m.start, m.end);
     if (jmax==0) { //also single-exon ref
         ovlen=mseg.overlapLen(r.start,r.end);
         int lmax=GMAX(r.covlen, m.covlen);
         if (ovlen >= lmax*0.8) return '='; //fuzz matching for single-exon transcripts: 80% of the longer one
         //if (m.covlen<=ovlen+12 && m.covlen<r.covlen) return 'c'; //contained
         if (m.covlen<r.covlen && ovlen >= m.covlen*0.8) return 'c';
         return 'o'; //just plain overlapping
         }
     //single-exon qry overlaping multi-exon ref
     for (int j=0;j<=jmax;j++) {
        //check if it's contained by an exon
        if (m.start>r.exons[j]->start-8 && m.end<r.exons[j]->end+8)
            return 'c';
        if (j==jmax) break;
        //check if it's contained by an intron
        if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
           return 'i';
        // check if it's a potential pre-mRNA transcript
        // (if overlaps an intron at least 10 bases)
        uint iovlen=mseg.overlapLen(r.exons[j]->end+1, r.exons[j+1]->start-1);
        if (iovlen>=10 && mseg.len()>iovlen+10) return 'e';
        }
     return 'o'; //plain overlap, uncategorized
     } //single-exon transfrag
  //-- from here on we have a multi-exon transfrag --
  // * check if contained by a ref intron
  for (int j=0;j<jmax;j++) {
     if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end) 
           return 'i';
     }
  //> check if m's intron chain is a subset of  r's intron chain
  int imax=m.exons.Count()-1;// imax>0 here
  if (m.exons[imax]->start<r.exons[0]->end ||
      r.exons[jmax]->start<m.exons[0]->end ) //intron chains do not overlap at all
           return 'o'; //but terminal exons do, otherwise we wouldn't be here
  int i=1; //index of exon to the right of current qry intron
  int j=1; //index of exon to the right of current ref intron
  //find first intron overlap
  while (i<=imax && j<=jmax) {
     if (r.exons[j]->start<m.exons[i-1]->end) { j++; continue; }
     if (m.exons[i]->start<r.exons[j-1]->end) { i++; continue; }
     break; //here we have an intron overlap
     }
  if (i>imax || j>jmax)
      return 'o'; //no initial intron overlap found
  //from here on we check all qry introns against ref introns
  bool jmatch=false; //true if at least a junction match is found
  bool icmatch=(i==1); //intron chain match - it will be updated as introns are checked
  //bool exovli=false; // if any terminal exon of qry extends into a ref intron
  int jmstart=j; //index of first intron overlap of reference
  int jmend=0;  //index of last intron overlap of reference
  int imend=0;  //index of last intron overlap of query
  //check for intron matches
  while (i<=imax && j<=jmax) {
    uint mstart=m.exons[i-1]->end;
    uint mend=m.exons[i]->start;
    uint rstart=r.exons[j-1]->end;
    uint rend=r.exons[j]->start;
    if (rend<mstart) { j++; icmatch=false; continue; } //skipping ref intron, no ichain match
    if (mend<rstart) { i++; icmatch=false; continue; } //skipping qry intron, no ichain match
    //overlapping introns here, test junction matching
    jmend=j; //keep track of last overlapping intron
    imend=i;
    bool smatch=(mstart==rstart);
    bool ematch=(mend==rend);
    if (smatch || ematch) jmatch=true;
    if (smatch && ematch) { i++; j++; } //perfect match for this intron
                     else { //at least one junction doesn't match
                          icmatch=false;
                          if (mend>rend) j++; else i++;
                          }
    } //while checking intron overlaps
  
  if (icmatch && imend==imax) { // qry intron chain match
     if (jmstart==1 && jmend==jmax) return '='; //identical intron chains
     // -- qry intron chain is shorter than ref intron chain --
     int l_iovh=0;   // overhang of leftmost q exon left boundary beyond the end of ref intron to the left
     int r_iovh=0;   // same type of overhang through the ref intron on the right
     if (jmstart>1 && r.exons[jmstart-1]->start>m.start) 
        l_iovh = r.exons[jmstart-1]->start - m.start;
     if (jmend<jmax && m.end > r.exons[jmend]->end)
        r_iovh = m.end - r.exons[jmend]->end;
     if (l_iovh<4 && r_iovh<4) return 'c';
     //TODO? check if any x_iovl>10 and return 'e' to signal an "unspliced intron" ?
     // or we can check if any of them are >= the length of the corresponding ref intron on that side
     return 'j';
     }
   /*
  if (icmatch && (jmax>=imax)) { //all qry introns match
       //but they may overlap
     // if ((lbound && lbound > m.exons[0]->start+10) ||
     //    (j<=jmax && m.exons[i-1]->end > r.exons[j-1]->end+10)) return 'j';
     //       return 'c';
     // }
    int code = 'c';
    if (lbound)
    {
      uint ref_boundary = lbound;
      uint cuff_boundary = m.exons[0]->start;
      if (ref_boundary > (cuff_boundary + pre_mrna_threshold)) // cuff extends a lot
      {
        code = 'j';
      }
      if (ref_boundary > cuff_boundary) // cuff extends just a bit into a ref intron
      {
        code = 'e';
      }
    }
    if (j <= jmax)
    {
      uint ref_boundary = r.exons[j-1]->end;
      uint cuff_boundary = m.exons[i-1]->end;
      if (cuff_boundary > (ref_boundary + pre_mrna_threshold)) // cuff extends a lot
      {
        code = 'j';
      }
      if (cuff_boundary > ref_boundary) // cuff extends just a bit into a ref intron
      {
        code = 'e';
      }
    }
    //if ((lbound && lbound > m.exons[0]->start+10) ||
    //  (j<=jmax && m.exons[i-1]->end > r.exons[j-1]->end+10)) return 'j';
    return code;
     }
  
//  if (!ichain) // first and last exons more or less match, but there's a different intron somewhere
//  {
//    
//  }
  
  */
  return jmatch ? 'j':'o';
}

char getRefOvl(GffObj& m, GLocus& rloc, GffObj*& rovl, int& ovlen) {
  rovl=NULL;
  ovlen=0;
  if (m.start>rloc.end || m.end<rloc.start) {
     //see if it's a polymerase run ?
     /*
       if ((m.strand=='+' && m.end<=rloc.end+polyrun_range) ||
         (m.strand=='-' && m.start>=rloc.start-polyrun_range)) {
            rovl=rloc.mrna_maxcov;
            ((CTData*)m.uptr)->addOvl('p',rloc.mrna_maxcov);
            return 'p';
            }
     */
     return 0; //unknown -> intergenic space
     }
  for (int i=0;i<rloc.mrnas.Count();i++) {
     GffObj* r=rloc.mrnas[i];
     int olen=0;
     char ovlcode=getOvlCode(m,*r,olen);
     if (ovlcode!=0) { //has some sort of "overlap" with r
       ((CTData*)m.uptr)->addOvl(ovlcode,r,olen);
       if (olen>ovlen) ovlen=olen;
       if (ovlcode=='c' || ovlcode=='=') //keep match/containment for each reference transcript
          ((CTData*)r->uptr)->addOvl(ovlcode,&m,olen);
       }
     }//for each ref in rloc
  // i,j,o
  return ((CTData*)m.uptr)->getBestCode();
}

/*
void findTMatches(GTrackLocus& loctrack, int qcount) {
 //perform an all vs. all ichain-match for all transcripts across all loctrack[i]->qloci
for (int q=0;q<qcount-1;q++) { //for each qry dataset
  if (loctrack[q]==NULL) continue;
  for (int qi=0;qi<loctrack[q]->Count();qi++) { // for each transcript in q dataset
    GffObj* qi_t=loctrack[q]->Get(qi);
    CTData* qi_d=(CTData*)qi_t->uptr;
    if ((qi_d->eqdata & EQHEAD_TAG) !=0) { //this is set as an EQ chain head already
         //if (qi_t->exons.Count()>1) 
          continue; 
         }
    for (int n=q+1;n<qcount;n++) { // for every successor dataset
       if (loctrack[n]==NULL) continue;
       for (int ni=0;ni<loctrack[n]->Count();ni++) {
          GffObj* ni_t=loctrack[n]->Get(ni);
          CTData* ni_d=(CTData*)ni_t->uptr;
          //if (ni_d->eqdata!=0) continue; //already part of an EQ chain
          //single exon transfrags have special treatment:
          bool s_match=(ni_t->exons.Count()==1 && qi_t->exons.Count()==1);
          if (ni_d->eqdata!=0 && !s_match) continue; //already part of an EQ chain
          if (ni_d->eqnext!=NULL) continue;
          int ovlen=0;
          
          if (qi_d->eqnext!=NULL) {
              if (!s_match) continue;
              //test all in the EQ list for a match
              bool matchFound=false;
              CTData* next_eq_d=qi_d;
              if (tMatch(*qi_t, *ni_t, ovlen, true)) {
                 matchFound=true;
                 }
               else {
                 while (next_eq_d->eqnext!=NULL) {
                     if (tMatch(*(next_eq_d->eqnext), *ni_t, ovlen, true)) {
                         matchFound=true;
                         break;
                         }
                     next_eq_d=(CTData*)(next_eq_d->eqnext->uptr);
                     } //last in the chain
                 }
              if (matchFound) {
                 //add this to the end of the EQ chain instead
                 next_eq_d=(CTData*)(qi_d->eqnext->uptr);
                 while (next_eq_d->eqnext!=NULL) {
                        next_eq_d=(CTData*)(next_eq_d->eqnext->uptr);
                        } //last in the chain
                 next_eq_d->eqnext=ni_t;
                 ni_d->eqdata=n+1;
                 ni_d->eqnext=NULL; ///TEST
                 }
              }
           else {
              if (tMatch(*qi_t,*ni_t, ovlen, true)) {
                 qi_d->eqnext=ni_t;
                 ni_d->eqnext=NULL; ///TEST
                 if (qi_d->eqdata == 0) {//only start of chain is tagged
                   qi_d->eqdata = ((q+1) | EQHEAD_TAG);
                   }
                 ni_d->eqdata=n+1; //EQ chain member only marked with qry# (1-based)
                 }
              } //multi-exon case
          } //for each transfrag in the next qry dataset
       
       if (qi_d->eqnext!=NULL && qi_t->exons.Count()>1) break;
         //part of a chain already, skip other datasets
       } // for each successor dataset
    } //for each transcript in qry dataset
  } //for each qry dataset
}
*/

void findTMatches(GTrackLocus& loctrack, int qcount) {
 //perform an all vs. all ichain-match for all transcripts across all loctrack[i]->qloci
for (int q=0;q<qcount-1;q++) { //for each qry dataset
  if (loctrack[q]==NULL) continue;
  for (int qi=0;qi<loctrack[q]->Count();qi++) { // for each transcript in q dataset
    GffObj* qi_t=loctrack[q]->Get(qi);
    CTData* qi_d=(CTData*)qi_t->uptr;
    if (qi_d->eqlist!=NULL && qi_t->exons.Count()>1) {
         continue; //this is part of an EQ chain already
         }
    for (int n=q+1;n<qcount;n++) { // for every successor dataset
       if (loctrack[n]==NULL) continue;
       for (int ni=0;ni<loctrack[n]->Count();ni++) {
          GffObj* ni_t=loctrack[n]->Get(ni);
          CTData* ni_d=(CTData*)ni_t->uptr;
          bool singleExon=(ni_t->exons.Count()==1 && qi_t->exons.Count()==1);
          if (ni_d->eqlist!=NULL && 
                 (ni_d->eqlist==qi_d->eqlist || !singleExon)) continue;
          int ovlen=0;
          if ((ni_d->eqlist==qi_d->eqlist && qi_d->eqlist!=NULL) ||
                  tMatch(*qi_t,*ni_t, ovlen, singleExon)) {
             qi_d->joinEqList(ni_t);
             }
          }
       } // for each successor dataset
    } //for each transcript in qry dataset
  } //for each qry dataset
}


int cmpTData_qset(const pointer* p1, const pointer* p2) {
 CTData* d1=(CTData*)(((GffObj*)p1)->uptr);
 CTData* d2=(CTData*)(((GffObj*)p2)->uptr);
 return (d1->qset - d2->qset);
 }

void printITrack(FILE* ft, GList<GffObj>& mrnas, int qcount, int& cnum) {
  for (int i=0;i<mrnas.Count();i++) {
   GffObj& qt=*(mrnas[i]);
   CTData* qtdata=(CTData*)qt.uptr;
   int qfidx=qtdata->qset;
   char ovlcode=qtdata->classcode;
   //GList<GffObj> eqchain(false,false,false);
   CEqList* eqchain=qtdata->eqlist;
   GffObj* ref=NULL; //related ref -- it doesn't have to be fully matching
   GffObj* eqref=NULL; //fully ichain-matching ref
   GffObj* tcons=NULL; //"consensus" (largest) transcript for a clique
   int tmaxcov=0;
   //eqchain.Add(&qt);
   eqref=qtdata->eqref;
   if (qtdata->ovls.Count()>0 && qtdata->ovls[0]->mrna!=NULL) {
       //if it has ovlcode with a ref
       ref=qtdata->ovls[0]->mrna;
	   //consistency check: qtdata->ovls[0]->code==ovlcode
	   // -- let tcons be a transfrag, not a ref transcript
	   //tcons=eqref;
       //if (tcons!=NULL) tmaxcov=tcons->covlen;
       }
   //chain pre-check
   if (tcons==NULL || mrnas[i]->covlen>tmaxcov) {
       tcons=mrnas[i];
       tmaxcov=tcons->covlen;
       }
   if (qtdata->isEqHead()) {//head of a equivalency chain
      //check if all transcripts in this chain have the same ovlcode
      for (int k=0;k<qtdata->eqlist->Count();k++) {
         GffObj* m=qtdata->eqlist->Get(k);
        if (m->covlen>tmaxcov) {
            tmaxcov=m->covlen;
            tcons=m;
            }
        if (ovlcode!='=' && ovlcode!='.' && ((CTData*)m->uptr)->getBestCode()!=ovlcode) {
              ovlcode='.'; //non-uniform ovlcode
              }
         }
      /*
      GffObj* m=mrnas[i];
      while (((CTData*)m->uptr)->eqnext!=NULL) {
        m=((CTData*)m->uptr)->eqnext;
        eqchain.Add(m);
        if (m->covlen>tmaxcov) {
            tmaxcov=m->covlen;
            tcons=m;
            }
        if (ovlcode!='=' && ovlcode!='.' && ((CTData*)m->uptr)->getBestCode()!=ovlcode) {
              ovlcode='.'; //non-uniform ovlcode
              //break;
              }
        } //while elements in chain
       */
       
      }//chain check
   //if (ovlcode=='p') ref=NULL; //ignore polymerase runs?
   if (ovlcode==0 || ovlcode=='-') ovlcode = (ref==NULL) ? 'u' : '.';
   //-- print columns 1 and 2 as LOCUS_ID and TCONS_ID
   //bool chainHead=(qtdata->eqnext!=NULL && ((qtdata->eqdata & EQHEAD_TAG)!=0));
   bool chainHead=qtdata->isEqHead();
   //bool noChain=((qtdata->eqdata & EQCHAIN_TAGMASK)==0);
   bool noChain=(eqchain==NULL);
   if (chainHead || noChain) {
     cnum++;
     if (ft!=NULL) fprintf(ft,"%s_%08d\t",cprefix,cnum);
     GXLocus* xloc=qtdata->locus->xlocus;
     if (xloc!=NULL) {
         if (ft!=NULL) fprintf(ft, "XLOC_%06d\t",xloc->id);
         if (tcons->exons.Count()>1) {
            //! only multi-exon mRNAs are counted for multi-transcript xloci !
              xloc->num_mtcons++;
              if (xloc->num_mtcons==2)
                 total_xloci_alt++;
              }
         }
      else {
        //should NEVER happen!
        int fidx=qtdata->qset;
        GError("Error: no XLocus created for transcript %s (file %s) [%d, %d], on %s%c:%d-%d\n", qt.getID(),
                           qryfiles[qtdata->locus->qfidx]->chars(), qtdata->locus->qfidx, fidx, qt.getGSeqName(), qt.strand, qt.start, qt.end);
        }
     addXCons(xloc, ref, ovlcode, tcons, eqchain);
     } // if chain head or uniq entry (not part of a chain)
   if (ft==NULL) continue;
   if (chainHead) {
      //this is the start of a equivalence class as a printing chain
      if (ref!=NULL) fprintf(ft,"%s|%s\t%c", getGeneID(ref),ref->getID(), ovlcode);
                else fprintf(ft,"-\t%c", ovlcode);
      GffObj* m=mrnas[i];
      CTData* mdata=(CTData*)m->uptr;
      
      int lastpq=-1;
      /*
      for (int ptab=mdata->qset-lastpq; ptab>0;ptab--)
             if (ptab>1) fprintf(ft,"\t-");
                    else fprintf(ft,"\t");
      lastpq=mdata->qset;
      fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%8.6f|%d", lastpq+1, getGeneID(m), m->getID(),
          iround(m->gscore/10), mdata->FPKM, mdata->conf_lo, mdata->conf_hi, mdata->cov, m->covlen);
      //traverse linked list of matching transcripts
      while (mdata->eqnext!=NULL) {
         m=mdata->eqnext;
         mdata=(CTData*)m->uptr;
         for (int ptab=mdata->qset-lastpq;ptab>0;ptab--)
             if (ptab>1) fprintf(ft,"\t-");
                    else fprintf(ft,"\t");
         lastpq = mdata->qset;
         fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%8.6f|%d", lastpq+1, getGeneID(m), m->getID(),
            iround(m->gscore/10), mdata->FPKM,mdata->conf_lo,mdata->conf_hi,mdata->cov, m->covlen);
         } //traverse and print row
      */
      eqchain->setUnique(false);
      eqchain->setSorted((GCompareProc*) cmpTData_qset);
      
      for (int k=0;k<eqchain->Count();k++) {
         m=eqchain->Get(k);
         mdata=(CTData*)m->uptr;
         if (mdata->qset==lastpq) {
            //shouldn't happen, unless this input set is messed up (has duplicates/redundant transfrags)
            fprintf(ft,",%s|%s|%d|%8.6f|%8.6f|%8.6f|%8.6f|%d", getGeneID(m), m->getID(),
               iround(m->gscore/10), mdata->FPKM,mdata->conf_lo,mdata->conf_hi,mdata->cov, m->covlen);
            continue;
            }
         for (int ptab=mdata->qset-lastpq;ptab>0;ptab--)
             if (ptab>1) fprintf(ft,"\t-");
                    else fprintf(ft,"\t");
         lastpq = mdata->qset;
         fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%8.6f|%d", lastpq+1, getGeneID(m), m->getID(),
            iround(m->gscore/10), mdata->FPKM,mdata->conf_lo,mdata->conf_hi,mdata->cov, m->covlen);
         }
      for (int ptab=qcount-lastpq-1;ptab>0;ptab--)
            fprintf(ft,"\t-");
      fprintf(ft,"\n");
      continue;
      } //start of eq class (printing chain)

   if (eqchain!=NULL) continue; //part of a matching chain, dealt with previously

   //--------- not in an ichain-matching class, print as singleton

   if (ref!=NULL) fprintf(ft,"%s|%s\t%c",getGeneID(ref), ref->getID(), ovlcode);
             else fprintf(ft,"-\t%c",ovlcode);
   for (int ptab=qfidx;ptab>=0;ptab--)
      if (ptab>0) fprintf(ft,"\t-");
             else fprintf(ft,"\t");
   fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%8.6f|-",qfidx+1, getGeneID(qt), qt.getID(),iround(qt.gscore/10),
       qtdata->FPKM, qtdata->conf_lo,qtdata->conf_hi,qtdata->cov);
   for (int ptab=qcount-qfidx-1;ptab>0;ptab--)
         fprintf(ft,"\t-");
   fprintf(ft,"\n");
   } //for each transcript
}


void findTRMatch(GTrackLocus& loctrack, int qcount, GLocus& rloc) {
 //requires loctrack to be already populated with overlapping qloci by findTMatches()
  // which also found (and tagged) all matching qry transcripts
 for (int q=0;q<qcount;q++) { //for each qry dataset
  if (loctrack[q]==NULL) continue;
  for (int qi=0;qi<loctrack[q]->Count();qi++) { // for each transcript in q dataset
    //if (loctrack[q]->cl[qi]->exons.Count()<2) continue; //skip single-exon transcripts
    GffObj& qt=*(loctrack[q]->Get(qi));
    CTData* qtdata=(CTData*)qt.uptr;
    GffObj* rmatch=NULL; //== ref match for this row
    int rovlen=0;
    //if (qtdata->eqnext!=NULL && ((qtdata->eqdata & EQHEAD_TAG)!=0)) { 
    if (qtdata->isEqHead()) {
        //EQ chain head -- transfrag equivalency list starts here
        if (qtdata->eqref==NULL) { //find rloc overlap
           if (qt.overlap(rloc.start, rloc.end)) {
                 rmatch=findRefMatch(qt, rloc, rovlen);
                 }
           } else rmatch=qtdata->eqref;
       if (rmatch!=NULL) {
         /*
          GffObj* m=loctrack[q]->Get(qi);
          //traverse linked list of matching transcripts
          while (((CTData*)m->uptr)->eqnext!=NULL) {
            m=((CTData*)m->uptr)->eqnext;
            if (rmatch!=NULL) {
              ((CTData*)m->uptr)->addOvl('=',rmatch,rovlen);
              }
            } //traverse qry data sets
          continue;
          }
         */
          for (int k=0;k<qtdata->eqlist->Count();k++) {
            GffObj* m=qtdata->eqlist->Get(k);
            ((CTData*)m->uptr)->addOvl('=',rmatch,rovlen);
            continue;
            }
         }
        //if (rmatch!=NULL) continue;
        } //equivalence class (chain of intron-matching)
    //if ((qtdata->eqdata & EQCHAIN_TAGMASK)!=0) continue; //part of a matching chain, dealt with previously

    //--------- qry mrna not in a '=' matching clique
    if (qtdata->eqref==NULL) { //find any rloc overlap
       if (qt.overlap(rloc.start, rloc.end)) {
          rmatch=findRefMatch(qt, rloc, rovlen);
          if (rmatch==NULL) {
            //not an ichain match, look for other codes
            GffObj* rovl=NULL;
            int rovlen=0;
            //char ovlcode=
            getRefOvl(qt, rloc,rovl,rovlen);
            }
          }
       }
     else rmatch=qtdata->eqref;
    } //for each qry transcript
  }//for each qry dataset
}


bool inPolyRun(char strand, GffObj& m, GList<GLocus>* rloci, int& rlocidx) {
 //we are only here if there is no actual overlap between m and any locus in rloci
 if (rloci==NULL || rloci->Count()==0) return false; // || m.exons.Count()>1
  if (strand=='-') {
        rlocidx=qsearch_loci(m.end, *rloci);
        //returns index of locus starting just ABOVE m.end
        // or -1 if last locus start <= m.end
        GLocus* rloc=NULL;
        if (rlocidx<0) return false;
        while (rlocidx<rloci->Count()) {
           rloc=rloci->Get(rlocidx);
           if (rloc->start>m.end+polyrun_range) break;
           if (rloc->start+6>m.end) return true;
           rlocidx++;
           }
        }
      else { // strand == '+' (or '.' ?)
        rlocidx=qsearch_loci(m.end, *rloci);
        GLocus* rloc=NULL;
        //returns index of closest locus starting ABOVE m.end
        // or -1 if last locus start <= m.end
        if (rlocidx<0) rlocidx=rloci->Count(); //this may actually start below m.end
        while ((--rlocidx)>=0) {
          rloc=rloci->Get(rlocidx);
          if (m.start>rloc->start+GFF_MAX_LOCUS) break;
          if (m.start+6>rloc->end && m.start<rloc->end+polyrun_range) return true;
          }
        }
  return false;
}

CTData* getBestOvl(GffObj& m) {
 //CTData* mdata=(CTData*)m.uptr;
 //return mdata->getBestCode();
  if ( ((CTData*)m.uptr)->ovls.Count()>0)
     return (CTData*)m.uptr;
  return NULL;
}

void reclass_XStrand(GList<GffObj>& mrnas, GList<GLocus>* rloci) {
  //checking for relationship with ref transcripts on opposite strand
  if (rloci==NULL || rloci->Count()<1) return;
  int j=0;//current rloci index
  for (int i=0;i<mrnas.Count();i++) {
     GffObj& m=*mrnas[i];
     char ovlcode=((CTData*)m.uptr)->getBestCode();
     if (ovlcode>47 && strchr("=cjeo",ovlcode)!=NULL) continue;
     GLocus* rloc=rloci->Get(j);
     if (rloc->start>m.end) continue; //check next transfrag
     while (m.start>rloc->end && j+1<rloci->Count()) {
           j++;
           rloc=rloci->Get(j);
           }
     if (rloc->start>m.end) continue; //check next transfrag
     //m overlaps rloc:
     //check if m has a fuzzy intron overlap -> 's' (shadow, mapped on the wrong strand)
     //  then if m is contained within an intron -> 'i'
     //  otherwise it's just a plain cross-strand overlap: 'x'
     int jm=0;
     do { //while rloci overlap this transfrag (m)
       rloc=rloci->Get(j+jm);
       bool is_shadow=false;
       GffObj* sovl=NULL;
       bool is_intraintron=false;
       GffObj* iovl=NULL;
       if (rloc->introns.Count()>0) {
           for (int n=0;n<rloc->introns.Count();n++) {
              GISeg& rintron=rloc->introns[n];
              if (rintron.start>m.end) break;
              if (m.start>rintron.end) continue;
              //overlap between m and intron
              if (m.end<=rintron.end && m.start>=rintron.start) {
                  is_intraintron=true;
                  if (iovl==NULL || iovl->covlen<rintron.t->covlen) iovl=rintron.t;
                  continue;
                  }
              //check if any intron of m has a fuzz-match with rintron
              for (int e=1;e<m.exons.Count();e++) {
                 GSeg mintron(m.exons[e-1]->end+1,m.exons[e]->start-1);
                 if (rintron.coordMatch(&mintron,10)) {
                    is_shadow=true;
                    if (sovl==NULL || sovl->covlen<rintron.t->covlen) sovl=rintron.t;
                    break;
                    }
                 } //for each m intron
              } //for each intron of rloc
           }//rloc has introns
       bool xcode=true;
       if (is_shadow) { ((CTData*)m.uptr)->addOvl('s', sovl); xcode=false; }
             // else
       if (ovlcode!='i' && is_intraintron) { ((CTData*)m.uptr)->addOvl('i', iovl); xcode=false; }
       if (xcode) {
               // just plain overlap, find the overlapping mrna in rloc
               GffObj* maxovl=NULL;
               int ovlen=0;
               GffObj* max_lovl=NULL; //max len ref transcript
                       // having no exon overlap but simply range overlap (interleaved exons)
               for (int ri=0;ri<rloc->mrnas.Count();ri++) {
                  if (!m.overlap(*(rloc->mrnas[ri]))) continue;
                  int o=m.exonOverlapLen(*(rloc->mrnas[ri]));
                  if (o>0) {
                     if (o>ovlen) {
                        ovlen=o;
                        maxovl=rloc->mrnas[ri];
                        }
                     }
                    else { //no exon overlap, but still overlapping (interleaved exons)
                     if (max_lovl==NULL || max_lovl->covlen<rloc->mrnas[ri]->covlen)
                         max_lovl=rloc->mrnas[ri];
                     }
                  }
               if (maxovl) ((CTData*)m.uptr)->addOvl('x',maxovl);
                 else if (max_lovl) ((CTData*)m.uptr)->addOvl('x',max_lovl);
               } //'x'
       jm++;
       } while (j+jm<rloci->Count() && rloci->Get(j+jm)->overlap(m));
     } //for each transfrag
}

void reclass_mRNAs(char strand, GList<GffObj>& mrnas, GList<GLocus>* rloci, GFaSeqGet *faseq) {
  int rlocidx=-1;
  for (int i=0;i<mrnas.Count();i++) {
    GffObj& m=*mrnas[i];
    char ovlcode=((CTData*)m.uptr)->getBestCode();
    //if (ovlcode=='u' || ovlcode=='i' || ovlcode==0) {
    if (ovlcode=='u' || ovlcode<47) {
      //check for overlaps with ref transcripts on the other strand
      if (m.exons.Count()==1 && inPolyRun(strand, m, rloci, rlocidx)) {
         ((CTData*)m.uptr)->addOvl('p',rloci->Get(rlocidx)->mrna_maxcov);
         }
      else { //check for repeat content
         if (faseq!=NULL) {
            int seqlen;
            char* seq=m.getSpliced(faseq, false, &seqlen);
            //get percentage of lowercase
            int numlc=0;
            for (int c=0;c<seqlen;c++) if (seq[c]>='a') numlc++;
            if (numlc > seqlen/2)
               ((CTData*)m.uptr)->addOvl('r');
            GFREE(seq);
            }
         }
      } //for unassigned class
  }//for each mrna

}

void reclassLoci(char strand, GList<GLocus>& qloci, GList<GLocus>* rloci, GFaSeqGet *faseq) {
  for (int ql=0;ql<qloci.Count();ql++) {
    reclass_mRNAs(strand, qloci[ql]->mrnas, rloci, faseq);
    //find closest upstream ref locus for this q locus
  } //for each locus
}

//for a single genomic sequence, all qry data and ref data is stored in gtrack
//check for all 'u' transfrags if they are repeat ('r') or polymerase run 'p' or anything else
void umrnaReclass(int qcount,  GSeqTrack& gtrack, FILE** ftr, GFaSeqGet* faseq=NULL) {
    for (int q=0;q<qcount;q++) {
        if (gtrack.qdata[q]==NULL) continue; //no transcripts in this q dataset for this genomic seq
        reclassLoci('+', gtrack.qdata[q]->loci_f, gtrack.rloci_f, faseq);
        reclassLoci('-', gtrack.qdata[q]->loci_r, gtrack.rloci_r, faseq);
        reclass_mRNAs('+', gtrack.qdata[q]->umrnas, gtrack.rloci_f, faseq);
        reclass_mRNAs('-', gtrack.qdata[q]->umrnas, gtrack.rloci_r, faseq);
        //and also check for special cases with cross-strand overlaps:
        reclass_XStrand(gtrack.qdata[q]->mrnas_f, gtrack.rloci_r);
        reclass_XStrand(gtrack.qdata[q]->mrnas_r, gtrack.rloci_f);
        // print all tmap data here here:
        for (int i=0;i<gtrack.qdata[q]->tdata.Count();i++) {
            CTData* mdata=gtrack.qdata[q]->tdata[i];
            if (mdata->mrna==NULL) continue; //invalidated -- removed earlier
            //GLocus* rlocus=NULL;
            mdata->classcode='u';
            GffObj* ref=NULL;
            if (mdata->ovls.Count()>0) {
                mdata->classcode=mdata->ovls[0]->code;
                ref=mdata->ovls[0]->mrna;
            }
            //if (mdata->classcode<33) mdata->classcode='u';      
            if (mdata->classcode<47) mdata->classcode='u'; // if 0, '-' or '.'
            if (tmapFiles) {
                char ref_match_len[2048];
                if (ref!=NULL) {
                    sprintf(ref_match_len, "%d",ref->covlen);
                    fprintf(ftr[q],"%s\t%s\t",getGeneID(ref),ref->getID());
                    //rlocus=((CTData*)(ref->uptr))->locus;
                }
                else {
                    fprintf(ftr[q],"-\t-\t");
                    strcpy(ref_match_len, "-");
                }
                //fprintf(ftr[q],"%c\t%s\t%d\t%8.6f\t%8.6f\t%d\n", ovlcode, mdata->mrna->getID(),
                //    iround(mdata->mrna->gscore/10), mdata->FPKM, mdata->cov, mdata->mrna->covlen);
                const char* mlocname = (mdata->locus!=NULL) ? mdata->locus->mrna_maxcov->getID() : mdata->mrna->getID();
                fprintf(ftr[q],"%c\t%s\t%s\t%d\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%d\t%s\t%s\n", mdata->classcode, getGeneID(mdata->mrna), mdata->mrna->getID(),
                        iround(mdata->mrna->gscore/10), mdata->FPKM, mdata->conf_lo,mdata->conf_hi, mdata->cov, mdata->mrna->covlen, mlocname, ref_match_len);
            }
        } //for each tdata
    } //for each qdata
}

void buildXLoci(GTrackLocus& loctrack, int qcount, GSeqTrack& gtrack, char strand,
    GList<GXLocus>* retxloci=NULL) {
  GList<GXLocus>* dest_xloci=NULL;
  GList<GXLocus> tmpxloci(true,false,true); //local set of newly created xloci
  GList<GXLocus>* xloci=&tmpxloci;
  if (strand=='+') {
       dest_xloci=& gtrack.xloci_f;
       }
    else if (strand=='-') {
      dest_xloci = & gtrack.xloci_r;
      }
   else dest_xloci= & gtrack.xloci_u;

  if (retxloci==NULL) {
     //if no return set of build xloci was given
     //take it as a directive to work directly on the global xloci
     xloci=dest_xloci;
     dest_xloci=NULL;
   }
 for (int q=-1;q<qcount;q++) {
   GList<GLocus>* wrkloci=NULL;
   if (q<0) {
      if (loctrack.rloci.Count()==0) continue;
      //loci=new GList<GLocus>(true,false,false);
      //loci->Add(loctrack.rloc);
      wrkloci = &(loctrack.rloci);
      }
     else {
      if (loctrack[q]==NULL) continue;
      wrkloci = &(loctrack[q]->qloci);
      }
   
   for (int t=0;t<wrkloci->Count();t++) {
      GLocus* loc=wrkloci->Get(t);
      int xfound=0; //count of parent xloci
      if (loc->xlocus!=NULL) continue; //already assigned a superlocus
      GArray<int> mrgxloci(true);
      for (int xl=0;xl<xloci->Count();xl++) {
         GXLocus& xloc=*(xloci->Get(xl));
         if (xloc.start>loc->end) {
            if (xloc.start-loc->end > GFF_MAX_LOCUS) break;
            continue;
            }
         if (loc->start>xloc.end) continue;
         if (xloc.add_Locus(loc)) {
            xfound++;
            mrgxloci.Add(xl);
            }
         } //for each existing Xlocus
      if (xfound==0) {
         xloci->Add(new GXLocus(loc));
         }
      else {
         int il=mrgxloci[0];
         GXLocus& xloc=*(xloci->Get(il));
         if (xfound>1) {
            for (int l=1;l<xfound;l++) {
              int mlidx=mrgxloci[l]-l+1;
              xloc.addMerge(*(xloci->Get(mlidx)));
              GXLocus* ldel=xloci->Get(mlidx);
              xloci->Delete(mlidx);
              if (retxloci!=NULL)
                    delete ldel;
              }
            }
         //in case xloc.start was decreased, bubble-down until it's in the proper order
         while (il>0 && xloc<*(xloci->Get(il-1))) {
            il--;
            xloci->Swap(il,il+1);
            }
         } //at least one locus is being merged 
      }//for each locus
   }//for each set of loci in the region (refs and each qry set)
  //-- add xloci to the global set of xloci unless retxloci was given,
  if (retxloci!=NULL) retxloci->Add(*xloci);
                 else dest_xloci->Add(*xloci);
}

void singleQData(GList<GLocus>& qloci, GList<GTrackLocus>& loctracks) {
 for (int i=0;i<qloci.Count();i++) {
  if (qloci[i]->t_ptr==NULL) {
    GTrackLocus* tloc=new GTrackLocus();
    tloc->addQLocus(qloci[i],0);
    loctracks.Add(tloc);
    }
  }
}

void recheckUmrnas(GSeqData* gseqdata, GList<GffObj>& mrnas,
     GList<GLocus>& loci, GList<GLocus>& nloci,  GList<GLocus>& oloci) {
 GList<GLocus> reassignedLocs(false,false);
 for (int u=0;u<gseqdata->umrnas.Count();u++) {
   for (int l=0;l<oloci.Count();l++) {
     if (gseqdata->umrnas[u]==NULL) break;
     if (gseqdata->umrnas[u]->end<oloci[l]->start) break; //try next umrna
     if (oloci[l]->end<gseqdata->umrnas[u]->start) continue; //try next locus
     if (gseqdata->umrnas[u]->strand=='+' || gseqdata->umrnas[u]->strand=='-') {
       gseqdata->umrnas.Forget(u);
       continue; //already reassigned earlier
       }
     //umrna overlaps locus region
     GffObj* umrna=gseqdata->umrnas[u];
     for (int m=0;m<oloci[l]->mrnas.Count();m++) {
        if (oloci[l]->mrnas[m]->exonOverlap(umrna)) {
            gseqdata->umrnas.Forget(u);
            CTData* umdata=((CTData*)umrna->uptr);
            //must be in a Loci anyway
            if (umdata==NULL || umdata->locus==NULL)
                GError("Error: no locus pointer for umrna %s!\n",umrna->getID());
            for (int i=0;i<umdata->locus->mrnas.Count();i++) {
               GffObj* um=umdata->locus->mrnas[i];
               um->strand=oloci[l]->mrnas[m]->strand;
               }
            reassignedLocs.Add(umdata->locus);
            break;
            }
        } //for each mrna in locus
      } //for each locus
   } //for each umrna
 if (reassignedLocs.Count()>0) {
   gseqdata->umrnas.Pack();
   gseqdata->nloci_u.setFreeItem(false);
   for (int i=0;i<reassignedLocs.Count();i++) {
     GLocus* loc=reassignedLocs[i];
     for (int m=0;m<loc->mrnas.Count();m++) {
        mrnas.Add(loc->mrnas[m]);
        }
     loci.Add(loc);
     nloci.Add(loc);
     gseqdata->nloci_u.Remove(loc);
     }
   gseqdata->nloci_u.setFreeItem(true);
   }
}

void umrnasXStrand(GList<GXLocus>& xloci, GSeqTrack& gtrack) {
  //try to determine the strand of unoriented transfrags based on possible overlaps
  //with other, oriented transfrags
 for (int x=0;x<xloci.Count();x++) {
   if (xloci[x]->strand=='.') continue;
   if (xloci[x]->qloci.Count()==0) continue;
   //go through all qloci in this xlocus
   for (int l = 0; l < xloci[x]->qloci.Count(); l++) {
     char locstrand=xloci[x]->qloci[l]->mrna_maxcov->strand;
     if (locstrand=='.') {
        //this is a umrna cluster
        GLocus* qloc=xloci[x]->qloci[l];
        //we don't really need to update loci lists (loci_f, nloci_f etc.)
        /*
        if (xloci[x]->strand=='+') {
           }
         else { // - strand
           }
        */
        for (int i=0;i<qloc->mrnas.Count();i++) {
           qloc->mrnas[i]->strand=xloci[x]->strand;
           int uidx=gtrack.qdata[qloc->qfidx]->umrnas.IndexOf(qloc->mrnas[i]);
           if (uidx>=0) {
                gtrack.qdata[qloc->qfidx]->umrnas.Forget(uidx);
                gtrack.qdata[qloc->qfidx]->umrnas.Delete(uidx);
                if (xloci[x]->strand=='+')
                     gtrack.qdata[qloc->qfidx]->mrnas_f.Add(qloc->mrnas[i]);
                   else
                     gtrack.qdata[qloc->qfidx]->mrnas_r.Add(qloc->mrnas[i]);
                }
           }
        } //unknown strand
     } //for each xloci[x].qloci (l)

   } //for each xloci (x)
}

//cluster loci across all datasets
void xclusterLoci(int qcount, char strand, GSeqTrack& gtrack) {
  //gtrack holds data for all input qry datasets for a chromosome/contig
  //cluster QLoci
  GList<GTrackLocus> loctracks(true,true,false);
  //all vs all clustering across all qry data sets + ref
  //one-strand set of loci from all datasets + ref loci
  GList<GLocus>* wrkloci=NULL;
  //build xloci without references first
  //then add references only if they overlap an existing xloci
  
  int nq=0;
  for (int q=0;q<=qcount+1;q++) {
    bool refcheck=false;
    if (q==qcount) { // check the unoriented loci for each query file
       while (nq<qcount &&
              (gtrack.qdata[nq]==NULL || gtrack.qdata[nq]->nloci_u.Count()==0))
                 nq++; //skip query files with no unoriented loci
       if (nq<qcount) {
             wrkloci=&(gtrack.qdata[nq]->nloci_u);
             nq++;
             if (nq<qcount) q--; //so we can fetch the next nq in the next q cycle
             }
          else continue; //no more q files with unoriented loci
       }
    else if (q==qcount+1) { // check the reference loci
           if (strand=='+') wrkloci=gtrack.rloci_f;
                       else wrkloci=gtrack.rloci_r;
            
           if (wrkloci==NULL) break; //no ref loci here
           refcheck=true;
           }
     else  {
          if (gtrack.qdata[q]==NULL) continue;
          if (strand=='+') wrkloci=&(gtrack.qdata[q]->loci_f);
                      else wrkloci=&(gtrack.qdata[q]->loci_r);
         }
   // now do the all-vs-all clustering thing:
   for (int t=0;t<wrkloci->Count();t++) {
      GLocus* loc=wrkloci->Get(t);
      int xfound=0; //count of parent loctracks
      if (loc->t_ptr!=NULL) continue; //already assigned a loctrack
      GArray<int> mrgloctracks(true);
      for (int xl=0;xl<loctracks.Count();xl++) {
         GTrackLocus& trackloc=*loctracks[xl];
         if (trackloc.start>loc->end) break;
         if (loc->start>trackloc.end) continue;
         if (trackloc.add_Locus(loc)) {
            xfound++;
            mrgloctracks.Add(xl);
            }
         } //for each existing Xlocus
      if (xfound==0) {
         if (!refcheck) //we really don't care about ref-only clusters
           loctracks.Add(new GTrackLocus(loc));
         }
      else {
         int il=mrgloctracks[0];
         GTrackLocus& tloc=*(loctracks.Get(il));
         if (xfound>1) {
           for (int l=1;l<xfound;l++) {
             int mlidx=mrgloctracks[l]-l+1;
             tloc.addMerge(loctracks[mlidx], qcount, loc);
             loctracks.Delete(mlidx);
             }
           }
         //in case tloc.start was decreased, bubble-down 'til it's in the proper place
         while (il>0 && tloc<*(loctracks[il-1])) {
            il--;
            loctracks.Swap(il,il+1);
            }
        } //at least one locus found
      }//for each wrklocus
     } //for each set of loci (q)
   //loctracks is now set with all x-clusters on this strand
 for (int i=0;i<loctracks.Count();i++) {
   if (!loctracks[i]->hasQloci) continue; //we really don't care here about reference-only clusters
   GTrackLocus& loctrack=*loctracks[i];
   findTMatches(loctrack, qcount); //find matching transfrags in this xcluster
   for (int rl=0; rl < loctrack.rloci.Count(); rl++) {
      findTRMatch(loctrack, qcount, *(loctrack.rloci[rl]));
      //find matching reference annotation for this xcluster and assign class codes to transfrags
      }
    GList<GXLocus> xloci(false,false,false);
    buildXLoci(loctrack, qcount, gtrack, strand, &xloci);
    //the newly created xloci are in xloci
    umrnasXStrand(xloci, gtrack);
    //also merge these xloci into the global list of xloci
    for (int l=0; l < xloci.Count(); l++) {
       if (xloci[l]->strand=='+') {
           gtrack.xloci_f.Add(xloci[l]);
           }
          else if (xloci[l]->strand=='-') {
              gtrack.xloci_r.Add(xloci[l]);
              }
            else gtrack.xloci_u.Add(xloci[l]);
       }
    }//for each xcluster
}


void printRefMap(FILE** frs, int qcount, GList<GLocus>* rloci) {
  if (rloci==NULL) return;

  for (int l=0;l<rloci->Count(); l++) {
    for (int r=0;r<rloci->Get(l)->mrnas.Count(); r++) {
      GffObj& ref = *(rloci->Get(l)->mrnas[r]);
      CTData* refdata = ((CTData*)ref.uptr);
      GStr* clist = new GStr[qcount];
      GStr* eqlist = new GStr[qcount];
      for (int i = 0; i<refdata->ovls.Count(); i++) {
        GffObj* m=refdata->ovls[i]->mrna;
        char ovlcode=refdata->ovls[i]->code;
        if (m==NULL) {
          GMessage("Warning: NULL mRNA found for ref %s with ovlcode '%c'\n",
               ref.getID(), refdata->ovls[i]->code);
          continue;
        }
        int qfidx = ((CTData*)m->uptr)->qset;
        if (ovlcode == '=') {
          eqlist[qfidx].append(getGeneID(m));
          eqlist[qfidx].append('|');
          eqlist[qfidx].append(m->getID());
          eqlist[qfidx].append(',');
        }
        else if (ovlcode == 'c') {
          clist[qfidx].append(getGeneID(m));
          clist[qfidx].append('|');
          clist[qfidx].append(m->getID());
          clist[qfidx].append(',');
        }
      }//for each reference overlap
      for (int q=0;q<qcount;q++) {
        if (!eqlist[q].is_empty()) {
          eqlist[q].trimR(',');
          fprintf(frs[q],"%s\t%s\t=\t%s\n", getGeneID(ref), ref.getID(),eqlist[q].chars());
        }
        if (!clist[q].is_empty()) {
          clist[q].trimR(',');
          fprintf(frs[q],"%s\t%s\tc\t%s\n",getGeneID(ref), ref.getID(),clist[q].chars());
        }
      }
      delete[] clist;
      delete[] eqlist; 
    }// ref loop
  }//ref locus loop
}



void trackGData(int qcount, GList<GSeqTrack>& gtracks, GStr& fbasename, FILE** ftr, FILE** frs) {
  FILE* f_ltrack=NULL;
  FILE* f_itrack=NULL;
  FILE* f_ctrack=NULL;
  FILE* f_xloci=NULL;
  int cnum=0; //consensus numbering for printITrack()
  GStr s=fbasename;
  //if (qcount>1 || generic_GFF) { //doesn't make much sense for only 1 query file
    s.append(".tracking");
    f_itrack=fopen(s.chars(),"w");
    if (f_itrack==NULL) GError("Error creating file %s !\n",s.chars());
  //  }
  s=fbasename;
  s.append(".combined.gtf");
  f_ctrack=fopen(s.chars(),"w");
  if (f_ctrack==NULL) GError("Error creating file %s !\n",s.chars());
  
  s=fbasename;
  s.append(".loci");
  f_xloci=fopen(s.chars(),"w");
  if (f_xloci==NULL) GError("Error creating file %s !\n",s.chars());
  for (int g=0;g<gtracks.Count();g++) { //for each genomic sequence
    GSeqTrack& gseqtrack=*gtracks[g];

    xclusterLoci(qcount,  '+', gseqtrack);
    xclusterLoci(qcount,  '-', gseqtrack);

    //count XLoci, setting their id
    numXLoci(gseqtrack.xloci_f, xlocnum);
    numXLoci(gseqtrack.xloci_r, xlocnum);
    numXLoci(gseqtrack.xloci_u, xlocnum);
    //transcript accounting: for all those transcripts with 'u' or 0 class code
    // we have to check for polymerase runs 'p' or repeats 'r'

    GFaSeqGet *faseq=gfasta.fetch(gseqtrack.get_gseqid(), checkFasta);

    umrnaReclass(qcount, gseqtrack, ftr, faseq);

    // print transcript tracking (ichain_tracking)
    //if (qcount>1)
    for (int q=0;q<qcount;q++) {
         if (gseqtrack.qdata[q]==NULL) continue;
         printITrack(f_itrack, gseqtrack.qdata[q]->mrnas_f, qcount, cnum);
         printITrack(f_itrack, gseqtrack.qdata[q]->mrnas_r, qcount, cnum);
         //just for the sake of completion:
         printITrack(f_itrack, gseqtrack.qdata[q]->umrnas, qcount, cnum);
         }
    //print XLoci and XConsensi within each xlocus
    //also TSS clustering and protein ID assignment for XConsensi
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_f, faseq);
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_r, faseq);
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_u, faseq);
    if (tmapFiles && haveRefs) {
      printRefMap(frs, qcount, gseqtrack.rloci_f);
      printRefMap(frs, qcount, gseqtrack.rloci_r);
      }
    delete faseq;
    }
  if (tmapFiles) {
   for (int q=0;q<qcount;q++) {
        fclose(ftr[q]); 
        if (haveRefs) fclose(frs[q]); 
        }
   }
  if (f_ltrack!=NULL) fclose(f_ltrack);
  if (f_itrack!=NULL) fclose(f_itrack);
  if (f_ctrack!=NULL) fclose(f_ctrack);
  if (f_xloci!=NULL) fclose(f_xloci);

}
