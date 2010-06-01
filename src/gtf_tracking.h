#ifndef GTF_TRACKING_H
#define GTF_TRACKING_H
/*
 *  gtf_tracking.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 9/5/09.
 *  Copyright 2009 Geo Pertea. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "GList.hh"
#include "gff.h"

#define MAX_QFILES 500
extern const char* ATTR_GENE_NAME;

extern bool verbose;

extern bool largeScale;
//many input files, no accuracy stats are generated, no *.tmap
// exon attributes are discarded

int cmpByPtr(const pointer p1, const pointer p2);

class GLocus;

class COvLink {
public:
    static int coderank(char c) {
		switch (c) {
			case '=': return 0; //ichain match
			case 'c': return 2; //containment (chain fragment)
			case 'e': return 4; // overlap exon+intron of unspliced pre-mRNA
			case 'j': return 6; // overlap with at least a junction match
			case 'o': return 8; // overlap (exon)
			case 'r': return 14; //repeats
			case 'i': return 15; // intra-intron
			case 'p': return 16; //polymerase run
			case 'u': return 90; //intergenic
			case  0 : return 100;
			default: return 95;
        }
	}
    char code;
    int rank;
    GffObj* mrna;
    int ovlen;
    COvLink(char c=0,GffObj* m=NULL, int ovl=0) {
		code=c;
		mrna=m;
		ovlen=ovl;
		rank=coderank(c);
	}
    bool operator<(COvLink& b) {
		if (rank==b.rank)
			return (ovlen==b.ovlen)? (mrna>b.mrna) : (ovlen>b.ovlen);
		else return rank<b.rank;
	}
    bool operator>(COvLink& b) {
		if (rank==b.rank)
			return (ovlen==b.ovlen)? (mrna<b.mrna) : (ovlen<b.ovlen);
		else return rank>b.rank;
	}
    bool operator==(COvLink& b) {
		return (rank==b.rank && mrna==b.mrna);
	}
};

class CTData { //transcript associated data
public:
	GffObj* mrna;
	GLocus* locus;
	GList<COvLink> ovls; //overlaps with other transcripts (ref vs query)
	//-- just for ichain match tracking:
	GffObj* eqref; //ref having an ichain match
	int qset; //qry set index (qfidx)
	GffObj* eqnext;
	int eqdata;
	// Cufflinks specific data:
	double FPKM;
	double conf_hi;
	double conf_lo;
	double cov;
	char classcode; //the best/final classcode
	CTData(GffObj* m=NULL, GLocus* l=NULL):ovls(true,true,true) {
		mrna=m;
		if (mrna!=NULL) mrna->uptr=this;
		locus=l;
        classcode=0;
		eqref=NULL;
		eqnext=NULL;
		eqdata=0;
		qset=-2;
		FPKM=0;
		conf_lo=0;
		conf_hi=0;
		cov=0;
	}
	~CTData() {
		ovls.Clear();
	}
	void addOvl(char code,GffObj* target=NULL, int ovlen=0) {
		ovls.AddIfNew(new COvLink(code, target, ovlen));
	}
	char getBestCode() {
		return (ovls.Count()>0) ? ovls[0]->code : 0 ;
    }
	bool operator>(CTData& b) { return (mrna > b.mrna); }
	bool operator<(CTData& b) { return (mrna < b.mrna); }
	bool operator==(CTData& b) { return (mrna==b.mrna); }
};

class GSuperLocus;
class GTrackLocus;
class GXLocus;

//Data structure holding a query locus data (overlapping mRNAs on the same strand)
// and also the accuracy data of all mRNAs of a query locus
// (against all reference loci overlapping the same region)
class GLocus:public GSeg {
public:
    int gseq_id; //id of underlying genomic sequence
    int qfidx; // for locus tracking
    GTrackLocus* t_ptr; //for locus tracking cluster
    GffObj* mrna_maxcov;  //transcript with maximum coverage (for main "ref" transcript)
    GffObj* mrna_maxscore; //transcript with maximum gscore (for major isoform)
    GList<GffObj> mrnas; //list of transcripts (isoforms) for this locus
	GArray<GSeg> uexons; //list of unique exons (covered segments) in this region
	GArray<GSeg> mexons; //list of merged exons in this region
	GArray<GSeg> introns;
	GList<GLocus> cmpovl; //temp list of overlapping qry/ref loci to compare to (while forming superloci)
	
	//only for reference loci --> keep track of all superloci found for each qry dataset
	//                           which contain this reference locus
	GList<GSuperLocus>* superlst;
	GXLocus* xlocus; //superlocus formed by exon overlaps across all qry datasets
	// -- if genomic sequence was given:
	int spl_major; // number of GT-AG splice site consensi
	int spl_rare; // number of GC-AG, AT-AC and other rare splice site consensi
	int spl_wrong; //number of "wrong" (unrecognized) splice site consensi
	int ichains; //number of multi-exon mrnas
	int ichainTP;
	int ichainATP;
	int mrnaTP;
	int mrnaATP;
	int v; //user flag/data
	GLocus(GffObj* mrna=NULL, int qidx=-1):mrnas(true,false,false),uexons(true,true),mexons(true,true),
	introns(true,true), cmpovl(true,false,true) {
		//this will NOT free mrnas!
		ichains=0;
		gseq_id=-1;
		qfidx=qidx;
		t_ptr=NULL;
		creset();
		xlocus=NULL;
		mrna_maxcov=NULL;
		mrna_maxscore=NULL;
		superlst=new GList<GSuperLocus>(true,false,false);
		if (mrna!=NULL) {
			start=mrna->exons.First()->start;
			end=mrna->exons.Last()->end;;
			gseq_id=mrna->gseq_id;
			GSeg seg;
			for (int i=0;i<mrna->exons.Count();i++) {
				seg.start=mrna->exons[i]->start;
				seg.end=mrna->exons[i]->end;
				uexons.Add(seg);
				mexons.Add(seg);
				if (i>0) {
					seg.start=mrna->exons[i-1]->end+1;
					seg.end=mrna->exons[i]->start-1;
					introns.Add(seg);
				}
			}
			mrnas.Add(mrna);
			if (mrna->exons.Count()>1) ichains++;
			((CTData*)(mrna->uptr))->locus=this;
			mrna_maxscore=mrna;
			mrna_maxcov=mrna;
		}
	}
	~GLocus() {
		delete superlst;
	}
	void creset() {
		spl_major=0;spl_rare=0;spl_wrong=0;
		v=0; //visited/other data
		ichainTP=0;
		ichainATP=0;
		mrnaTP=0;
		mrnaATP=0;
		cmpovl.Clear();
	}
	
	void addMerge(GLocus& locus, GffObj* lnkmrna) {
		//add all the elements of the other locus (merging)
		//-- merge mexons
		GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
		int i=0; //index of first mexons with a merge
		int j=0; //index current mrna exon
		while (i<mexons.Count() && j<locus.mexons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=locus.mexons[j].start;
			uint jend=locus.mexons[j].end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//if (mexons[i].overlap(jstart, jend)) {
			//exon overlap was found :
			ovlexons.Add(j);
			//extend mexons[i] as needed
			if (jstart<istart) mexons[i].start=jstart;
			if (jend>iend) { //mexons[i] end extend
				mexons[i].end=jend;
				//now this could overlap the next mexon(s), so we have to merge them all
				while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
					uint nextend=mexons[i+1].end;
					mexons.Delete(i+1);
					if (nextend>mexons[i].end) {
						mexons[i].end=nextend;
						break; //no need to check next mexons
					}
                } //while next mexons merge
            } // mexons[i] end extend
			//  } //exon overlap
			j++; //check the next locus.mexon
		}
		//-- add the rest of the non-overlapping mexons:
		GSeg seg;
		for (int i=0;i<locus.mexons.Count();i++) {
			seg.start=locus.mexons[i].start;
			seg.end=locus.mexons[i].end;
			if (!ovlexons.Exists(i)) mexons.Add(seg);
		}
        // -- merge uexons
        //add to uexons:
		for (int i=0;i<locus.uexons.Count();i++) {
			uexons.Add(locus.uexons[i]);
		}
		for (int i=0;i<locus.introns.Count();i++) {
			introns.Add(locus.introns[i]);
		}
		
		// -- add locus.mrnas
		for (int i=0;i<locus.mrnas.Count();i++) {
			((CTData*)(locus.mrnas[i]->uptr))->locus=this;
			if (locus.mrnas[i]!=lnkmrna) {
				mrnas.Add(locus.mrnas[i]);
				if (locus.mrnas[i]->exons.Count()>1) ichains++;
            }
		  }
		// -- adjust start/end as needed
		if (start>locus.start) start=locus.start;
		if (end<locus.end) end=locus.end;
		if (mrna_maxcov->covlen<locus.mrna_maxcov->covlen)
			mrna_maxcov=locus.mrna_maxcov;
		if (mrna_maxscore->gscore<locus.mrna_maxscore->gscore)
			mrna_maxscore=locus.mrna_maxscore;
     }
	

	bool exonOverlap(GLocus& loc) {
		//check if any mexons overlap!
		int i=0;
		int j=0;
		while (i<mexons.Count() && j<loc.mexons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=loc.mexons[j].start;
			uint jend=loc.mexons[j].end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//exon overlap found if we're here:
			return true;
		}
		return false;
    }
	
	bool add_mRNA(GffObj* mrna) {
		if (mrnas.Count()>0 && mrna->gseq_id!=gseq_id) return false; //mrna must be on the same genomic seq
		//check for exon overlap with existing mexons
		//also update uexons and mexons accordingly, if mrna is added
		uint mrna_start=mrna->exons.First()->start;
		uint mrna_end=mrna->exons.Last()->end;
		if (mrna_start>end || start>mrna_end) return false;
		bool hasovl=false;
		int i=0; //index of first mexons with a merge
		int j=0; //index current mrna exon
		GArray<int> ovlexons(true,true); //list of mrna exon indexes overlapping mexons
		while (i<mexons.Count() && j<mrna->exons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=mrna->exons[j]->start;
			uint jend=mrna->exons[j]->end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//exon overlap found if we're here:
			ovlexons.Add(j);
			hasovl=true;
			//extend mexons[i] as needed
			if (jstart<istart) mexons[i].start=jstart;
			if (jend>iend) { //mexon stretch up
				mexons[i].end=jend;
				//now this could overlap the next mexon(s), so we have to merge them all
				while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
					uint nextend=mexons[i+1].end;
					mexons.Delete(i+1);
					if (nextend>mexons[i].end) {
						mexons[i].end=nextend;
						break; //no need to check next mexons
					}
				} //while next mexons merge
			} //possible mexons merge
			
			j++; //check the next mrna exon
		}//all vs all exon check loop
		if (hasovl) {
			GSeg seg;
	         //add the rest of the non-overlapping exons,
			 // and also to uexons etc.
			for (int i=0;i<mrna->exons.Count();i++) {
				seg.start=mrna->exons[i]->start;
				seg.end=mrna->exons[i]->end;
				if (!ovlexons.Exists(i)) mexons.Add(seg);
				uexons.Add(seg);
				if (i>0) {
					seg.start=mrna->exons[i-1]->end+1;
					seg.end=mrna->exons[i]->start-1;
					introns.Add(seg);
				}
			}
			
			mrnas_add(mrna);
			// add to mrnas
			((CTData*)mrna->uptr)->locus=this;
			gseq_id=mrna->gseq_id;
			if (mrna->exons.Count()>1) ichains++;
		}
		return hasovl;
	}
	
	//simpler,basic adding of a mrna
	void mrnas_add(GffObj* mrna) {
		mrnas.Add(mrna);
		// adjust start/end
		if (start>mrna->start) start=mrna->start;
		if (end<mrna->end) end=mrna->end;
		if (mrna_maxcov->covlen<mrna->covlen) mrna_maxcov=mrna;
		if (mrna_maxscore->gscore<mrna->gscore) mrna_maxscore=mrna;
    }
};

class GSuperLocus;
class GTrackLocus;

class GSuperLocus : public GSeg {
public:
    int qfidx; //index of query dataset/file for which this superlocus was created
    GList<GLocus> qloci;
    GList<GLocus> rloci;
    GList<GffObj> qmrnas; //list of transcripts (isoforms) for this locus
    GArray<GSeg> qmexons; //list of merged exons in this region
    GArray<GSeg> quexons; //list of unique exons (covered segments) in this region
    GArray<GSeg> qintrons; //list of unique exons (covered segments) in this region
    //same lists for reference:
    GList<GffObj> rmrnas; //list of transcripts (isoforms) for this locus
    GArray<GSeg> rmexons; //list of merged exons in this region
    GArray<GSeg> ruexons; //list of unique exons (covered segments) in this region
    GArray<GSeg> rintrons; //list of unique exons (covered segments) in this region
    // store problematic introns for printing:
    GArray<GSeg> i_missed; //missed reference introns (not overlapped by any qry intron)
    GArray<GSeg> i_notp;  //wrong ref introns (one or both ends not matching any qry intron)
    //
    GArray<GSeg> i_qwrong; //totally wrong qry introns (not overlapped by any ref intron)
    GArray<GSeg> i_qnotp;  //imperfect qry introns (may overlap but has no "perfect" match)
	
	
    int qbases_all;
    int rbases_all; //in fact, it's all ref bases overlapping any query loci
    int in_rmrnas; //count of ALL ref mrnas and loci given for this region
    int in_rloci; //not just those overlapping qry data
    // this will keep track of total qry loci, mrnas and exons in an area
    int total_superloci;
    int total_qloci;
    int total_qloci_alt; //total qloci with multiple transcripts

    int total_qmrnas;
    int total_qichains; //multi exon mrnas
    int total_qexons; //unique exons
    int total_qmexons;
    int total_qintrons; //unique introns
    // these ref totals are in fact only limited to data from
    // loci overlapping any of qry loci
    int total_rmexons;
    int total_richains; //multi exon mrnas
    int total_rloci;
    int total_rmrnas;
    int total_rexons;
    int total_rintrons; //unique introns
	
    //--- accuracy data after compared to ref loci:
	int locusTP; // 1 if mrnaTP>0
	int locusATP; // 1 if ichainATP + mrnaATP > 0
	int locusFP;
	int locusAFP;
	int locusAFN;
	int locusFN;
	//---transcript level accuracy -- all exon coordinates should match (most stringent)
	int mrnaTP; // number of qry mRNAs with perfect match with ref transcripts
	int mrnaFP; // number of qry mRNAs with no perfect match with a ref transcript
	int mrnaFN; // number of ref mRNAs in this region having no perfect match with a qry transcript
	int mrnaATP;
	int mrnaAFN;
	int mrnaAFP;
	//---intron level accuracy (comparing the ordered set of splice sites):
	int ichainTP; // number of qry intron chains covering a reference intron chain
	// (covering meaning that the ordered set of reference splice sites
	//  is the same with a ordered subset of the query splice sites)
	int ichainFP; // number of qry intron chains not covering a reference intron chain
	int ichainFN; // number of ref intron chains in this region not being covered by a reference intron chain
	// same as above, but approximate -- allowing a 10bp distance error for splice sites
	int ichainATP;
	int ichainAFP;
	int ichainAFN;
	//---projected features ---
	//---exon level accuracy:
	int exonTP;  //number of perfectly overlapping exons (true positives)
	int exonFP; //number of exons of query with no perfect match with a reference exon
	int exonFN; //number of exons of reference with no perfect match with a query exon
	// same as the above but with acceptable approximation (10bp error window):
	int exonATP;
	int exonAFP;
	int exonAFN;
	
	int intronTP;  //number of perfectly overlapping introns (true positives)
	int intronFP; //number of introns of query with no perfect match with a reference intron
	int intronFN; //number of introns of reference with no perfect match with a query intron
	// same as the above but with acceptable approximation (10bp error window):
	int intronATP;
	int intronAFP;
	int intronAFN;
	
	//-- EGASP added these too:
	int m_exons; //number of exons totally missed (not overlapped *at all* by any query exon)
	int w_exons; //numer of totally wrong exons (query exons not overlapping *at all* any reference exon)
	int m_introns; //number of introns totally missed (not overlapped *at all* by any query intron)
	int w_introns; //numer of totally wrong introns (query introns not overlapping *at all* any reference intron)
	int m_loci; //missed loci
	int w_loci; //novel/wrong loci
	//---base level accuracy
	int baseTP; //number of overlapping bases
	int baseFP; //number of qry bases not overlapping reference
	int baseFN; //number of ref bases not overlapping qry
	//            sorted,free,unique       sorted,unique
    GSuperLocus(uint lstart=0,uint lend=0):qloci(true,false,false),rloci(true,false,false),
	qmrnas(true,false,false), qmexons(true,false), quexons(true,false), qintrons(true,false),
	rmrnas(true,false,false), rmexons(true,false), ruexons(true,false), rintrons(true,false),
	i_missed(true,false),i_notp(true,false), i_qwrong(true,false), i_qnotp(true,false){
		qfidx=-1;
		start=lstart;
		end=lend;
		qbases_all=0;
		rbases_all=0;
		baseTP=0;baseFP=0;baseFN=0;
		locusTP=0;locusATP=0;
		locusFP=0;locusAFP=0;locusAFN=0;
		locusFN=0;
		in_rmrnas=0;
		in_rloci=0;
		w_loci=0;
		m_loci=0;
		total_superloci=0;
		mrnaTP=0;mrnaFP=0;mrnaFN=0;
		mrnaATP=0;mrnaAFP=0;mrnaAFN=0;
		ichainTP=0;ichainFP=0;ichainFN=0;
		ichainATP=0;ichainAFP=0;ichainAFN=0;
		exonTP=0;exonFP=0;exonFN=0;
		exonATP=0;exonAFP=0;exonAFN=0;
		intronTP=0;intronFP=0;intronFN=0;
		intronATP=0;intronAFP=0;intronAFN=0;
		total_rmexons=0;
		total_qmexons=0;
		total_qexons=0;total_qloci=0;total_qmrnas=0;
		total_qloci_alt=0;
		total_qintrons=0;total_qichains=0;
		total_rexons=0;total_rloci=0;total_rmrnas=0;
		total_rintrons=0;total_richains=0;
		w_exons=0;
		m_exons=0;
		w_introns=0;
		m_introns=0;
	}
    void addQlocus(GLocus& loc) {
		if (start==0 || start>loc.start) start=loc.start;
		if (end<loc.end) end=loc.end;
		qloci.Add(&loc);
		total_qloci++;
		if (loc.ichains>0 && loc.mrnas.Count()>1)
		    total_qloci_alt++;
		qmrnas.Add(loc.mrnas);
		total_qmrnas+=loc.mrnas.Count();
		total_qichains+=loc.ichains;
		qmexons.Add(loc.mexons);
		total_qmexons+=loc.mexons.Count();
		quexons.Add(loc.uexons);
		total_qexons+=loc.uexons.Count();
		qintrons.Add(loc.introns);
		total_qintrons+=loc.introns.Count();
	}
    void addRlocus(GLocus& loc) {
		if (start==0 || start>loc.start) start=loc.start;
		if (end<loc.end) end=loc.end;
		rloci.Add(&loc);
		total_rloci++;
		rmrnas.Add(loc.mrnas);
		total_rmrnas+=loc.mrnas.Count();
		total_richains+=loc.ichains;
		rmexons.Add(loc.mexons);
		total_rmexons+=loc.mexons.Count();
		ruexons.Add(loc.uexons);
		total_rexons+=loc.uexons.Count();
		rintrons.Add(loc.introns);
		total_rintrons+=loc.introns.Count();
	}
	
    void calcF() {
		// base level
		baseFP=qbases_all-baseTP;
		baseFN=rbases_all-baseTP;
		//exon level:
		exonAFP=total_qexons-exonATP;
		exonFP=total_qexons-exonTP;
		exonAFN=total_rexons-exonATP;
		exonFN=total_rexons-exonTP;
		//intron stats
		intronAFP=total_qintrons-intronATP;
		intronFP=total_qintrons-intronTP;
		intronAFN=total_rintrons-intronATP;
		intronFN=total_rintrons-intronTP;
		
		// ichain and transcript levels:
		ichainAFP=total_qichains-ichainATP;
		ichainFP=total_qichains-ichainTP;
		ichainAFN=total_richains-ichainATP;
		ichainFN=total_richains-ichainTP;
		mrnaFP=total_qmrnas-mrnaTP;
		mrnaFN=total_rmrnas-mrnaTP;
		mrnaAFP=total_qmrnas-mrnaATP;
		mrnaAFN=total_rmrnas-mrnaATP;
		// locus/gene level:
		locusAFP=total_qloci-locusATP;
		locusFP=total_qloci-locusTP;
		locusAFN=total_rloci-locusATP;
		locusFN=total_rloci-locusTP;
	}
	
    void addStats(GSuperLocus& s) {
		in_rmrnas+=s.in_rmrnas;
		in_rloci+=s.in_rloci;
		baseTP+=s.baseTP;
		exonTP+=s.exonTP;
		exonATP+=s.exonATP;
		intronTP+=s.intronTP;
		intronATP+=s.intronATP;
		ichainTP+=s.ichainTP;
		ichainATP+=s.ichainATP;
		mrnaTP+=s.mrnaTP;
		mrnaATP+=s.mrnaATP;
		locusTP+=s.locusTP;
		locusATP+=s.locusATP;
		m_exons+=s.m_exons;
		w_exons+=s.w_exons;
		m_introns+=s.m_introns;
		w_introns+=s.w_introns;
		if (s.total_superloci==0 && s.qloci.Count()>0) s.total_superloci=1;
		total_superloci+=s.total_superloci;
		qbases_all+=s.qbases_all;
		rbases_all+=s.rbases_all;
		m_loci+=s.m_loci;
		w_loci+=s.w_loci;
		total_qexons+=s.total_qexons;
		total_qintrons+=s.total_qintrons;
		total_qmexons+=s.total_qmexons;
		total_rexons+=s.total_rexons;
		total_rintrons+=s.total_rintrons;
		total_rmexons+=s.total_rmexons;
		total_qmrnas+=s.total_qmrnas;
		total_qichains+=s.total_qichains;
		total_rmrnas+=s.total_rmrnas;
		total_richains+=s.total_richains;
		total_qloci+=s.total_qloci;
		total_qloci_alt+=s.total_qloci_alt;
		total_rloci+=s.total_rloci;
    }
};

class GSeqData {
public:
	int gseq_id;
    GList<GffObj> refs_f; //forward strand mRNAs
    GList<GffObj> refs_r; //reverse strand mRNAs
	GList<GffObj> mrnas_f; //forward strand mRNAs
	GList<GffObj> mrnas_r; //reverse strand mRNAs
	GList<GLocus> loci_f; //forward strand loci
	GList<GLocus> loci_r; //reverse strand loci
	//--> the fields below are not used by reference data --
	GList<GSuperLocus> gstats_f; //stats for forward strand superloci
	GList<GSuperLocus> gstats_r; //stats for reverse strand superloci
	GList<GLocus> nloci_f; //"novel" loci on forward strand
	GList<GLocus> nloci_r; //"novel" loci on reverse strand
	GList<GffObj> umrnas; //unknown orientation mrnas
	GList<GLocus> nloci_u; //"novel" loci with no orientation found
	
	GList<CTData> tdata; //transcript data (uptr holder for all mrnas here)
	//--<
	GSeqData(int gid=-1):mrnas_f(true,true,false),mrnas_r(true,true,false),
	loci_f(true,true,true),loci_r(true,true,true),
	gstats_f(true,true,false),gstats_r(true,true,false),
	nloci_f(true,false,true), nloci_r(true,false,true),
	umrnas(true,true,false), nloci_u(true,true,true), tdata(false,true,false) {
		gseq_id=gid;
	}
	bool operator==(GSeqData& d){
		return (gseq_id==d.gseq_id);
	}
	bool operator>(GSeqData& d){
		return (gseq_id>d.gseq_id);
	}
	bool operator<(GSeqData& d){
		return (gseq_id<d.gseq_id);
	}
};


// a group of qry loci and a transcript cluster for a single qry dataset
class GQCluster : public GList<GffObj> {
 public:
   GffObj* mrna_maxcov;  //transcript with maximum coverage (for largest transcript)
   GffObj* mrna_maxscore; //transcript with maximum gscore ( = major isoform for Cufflinks)
   uint start;
   uint end;
   GList<GLocus> qloci;
   //GCluster cl; //just a more compact way of keeping all transcripts in these loci
   GQCluster(GList<GLocus>* loci=NULL):GList<GffObj>(true,false,false),
                                    qloci(true,false,false) {
     mrna_maxcov=NULL;
     mrna_maxscore=NULL;
     start=0;
     end=0;
     if (loci!=NULL)  {
          qloci.Add(*loci);
          for (int i=0;i<loci->Count();i++) {
             addLocus(loci->Get(i),false);
             }
          }
      }
   void addLocus(GLocus* loc, bool toLoci=true) {
     //check so we don't add locus duplicates
     if (toLoci) {
        for (int i=0;i<qloci.Count();i++) {
           if (loc==qloci[i]) return;
           }
        qloci.Add(loc);
        }
     for (int m=0;m<loc->mrnas.Count();m++) {
        GffObj* mrna=loc->mrnas[m];
        Add(mrna);
        if (start==0 || start>mrna->start) start=mrna->start;
        if (end<mrna->end) end=mrna->end;
        if (mrna_maxcov==NULL || mrna_maxcov->covlen<mrna->covlen) mrna_maxcov=mrna;
        if (mrna_maxscore==NULL || mrna_maxscore->gscore<mrna->gscore) mrna_maxscore=mrna;
        }
     }
};

//track a set of clustered qloci across multiple qry datasets
// the qloci in qcls[] overlap but not necessarily at exon level
// (so there can be multiple genes here in fact)
class GTrackLocus:public GSeg {
  public:
    char strand;
    bool hasQloci;
    //GLocus* rloc; //corresponding reference locus, if available
    GList<GLocus> rloci; //ref loci found overlapping this region
    GQCluster* qcls[MAX_QFILES]; //all qloci for this superlocus, grouped by dataset
    GTrackLocus(GLocus* qloc=NULL, int q=-1):GSeg(0,0),rloci(true,false,true) {
      strand='.';
      for (int i=0;i<MAX_QFILES;i++) qcls[i]=NULL;
      if (qloc!=NULL) addQLocus(qloc,q);
      }

    void addRLocus(GLocus* rl) {
      if (rl==NULL) return;
      if (rl->qfidx>=0)
          GError("Error: GTrackLocus::addRLocus called with a query locus (set# %d)\n",
                        rl->qfidx+1);
      if (strand=='.') strand=rl->mrna_maxcov->strand;
      if (start==0 || start>rl->start) start=rl->start;
      if (end==0 || end<rl->end) end=rl->end;
      rl->t_ptr=this;
      rloci.Add(rl);
      }

    void addQLocus(GLocus* loc, int q=-1) { //adding qry locus
      if (loc==NULL) return;
      if (strand=='.' && loc->mrna_maxcov->strand!='.')
           strand=loc->mrna_maxcov->strand;
      if (loc->qfidx<0 && q<0)
         GError("Error at GTrackLocus::addQLocus(): locus.qfidx not set and index not given!\n");
      if (q>=0) loc->qfidx=q;
           else q=loc->qfidx;
      if (start==0 || start>loc->start) start=loc->start;
      if (end==0 || end<loc->end) end=loc->end;
      if (qcls[q]==NULL) qcls[q]=new GQCluster();
      hasQloci=true;
      loc->t_ptr = this;
      qcls[q]->addLocus(loc);
      }

    bool add_Locus(GLocus* loc) {
      if (start==0 || overlap(*loc)) { //simple range overlap, not exon overlap
         if (loc->qfidx<0) addRLocus(loc);
                      else addQLocus(loc);
         return true;
         }
      return false;
      }


   void addQCl(int q, GQCluster* qcl, GLocus* lnkloc) {
      for (int i=0;i<qcl->qloci.Count();i++) {
         GLocus* loc=qcl->qloci[i];
         if (loc==lnkloc || loc->t_ptr==this) continue;
         hasQloci=true;
         loc->t_ptr=this;
         qcls[q]->addLocus(loc);
         }
     }

   void addMerge(GTrackLocus* loctrack, int qcount, GLocus* lnkloc) {
      if (loctrack==NULL) return;
      //merge qloci
      for (int q=0; q < qcount; q++) {
          if (qcls[q]==NULL) {
             if (loctrack->qcls[q]!=NULL) {
                 qcls[q]=loctrack->qcls[q];
                 loctrack->qcls[q]=NULL; //just move pointer here
                 //set all t_ptr pointers for moved loci
                 for (int ql = 0; ql < qcls[q]->qloci.Count(); ql++) {
                    qcls[q]->qloci[ql]->t_ptr=this;
                    }
                 hasQloci=true;
                 }
             }
           else //existing qloci at q
             if (loctrack->qcls[q]!=NULL) { //merge elements
                addQCl(q, loctrack->qcls[q], lnkloc);
              }
          }//for each qset
      //merge rloci, if any
      if (loctrack->rloci.Count()>0) {
         for (int l=0;l<loctrack->rloci.Count();l++) {
           if (loctrack->rloci[l]!=lnkloc && loctrack->rloci[l]->t_ptr!=this) {
              rloci.Add(loctrack->rloci[l]);
              loctrack->rloci[l]->t_ptr=this;
              }
           }
         }
      if (loctrack->start<start) start=loctrack->start;
      if (loctrack->end>end) end=loctrack->end;
      if (strand=='.' && loctrack->strand!='.') strand=loctrack->strand;
      }

    /*
    void add_QLoci(GList<GLocus>* loci, int q, GLocus& r) {
      // only add loci overlapping given refloc
      //rloc=&r;
      if (loci==NULL) return;
      for (int i=0;i<loci->Count();i++) {
         GLocus* loc=loci->Get(i);
         // if (!loc->exonOverlap(r)) continue;  //do we really needed exon overlap?
         if (!loc->overlap(r)) continue;
         if (start==0 || start>loc->start) start=loc->start;
         if (end==0 || end<loc->end) end=loc->end;
         loc->t_ptr=this;
         loc->qfidx=q;
         if (qcls[q]==NULL) qcls[q]=new GQCluster();
         qcls[q]->addLocus(loc);
         }
      strand=r.mrnas[0]->strand;
      }
     */
    ~GTrackLocus() {
      for (int q=0;q<MAX_QFILES;q++)
           if (qcls[q]!=NULL) { delete qcls[q]; qcls[q]=NULL; }
      }

    GQCluster* operator[](int q) {
      if (q<0 || q>=MAX_QFILES)
          GError("Error: qfidx index out of bounds (%d) for GTrackLocus!\n",q);
      return qcls[q];
      }
};

class GXConsensus:public GSeg {
 public:
   static int count;
   int id; //XConsensus ID
   int tss_id; //group id for those xconsensi with shared first exon
   int p_id; //group id for those xconsensi with "similar" protein
   GffObj* tcons; //chosen (longest) transcript to define the "consensus"
   GffObj* ref; //overlapping reference transcript 
   char refcode; // the code for ref relationship (like in the tracking file)
   char* aa;
   int aalen;
   //list of ichain-matching query (cufflinks) transcripts that contributed to this consensus
   GList<GffObj> qchain;
   GXConsensus(GffObj* c, GList<GffObj> qlst, GffObj* r=NULL, char rcode=0)
                   :qchain(false,false,false) {
      ref=r;
      refcode=rcode;
      tcons=c;
      qchain.Add(qlst);
      count++;
      tss_id=0;
      p_id=0;
      aalen=0;
      id=count;
      aa=NULL;
      start=tcons->start;
      end=tcons->end;
      }
   ~GXConsensus() {
     if (aa!=NULL) GFREE(aa);
     }
};

class GXLocus:public GSeg {
 public:
    int id;
    int num_mtcons; //number of multi-exon "consensus" transcripts in this locus
    char strand;
    GList<GLocus> rloci; //list of ref loci overlapping any of the mexons
    GList<GLocus> qloci; //loci from all qry datasets that have overlapping exons with this region
    GArray<GSeg> mexons; //list of merged exonic regions for this locus
    GList<GXConsensus> tcons;
    GXLocus(GLocus* lfirst=NULL):GSeg(0,0),
        rloci((GCompareProc*)cmpByPtr, (GFreeProc*)NULL, true),
        qloci((GCompareProc*)cmpByPtr, (GFreeProc*)NULL, true),
        mexons(true,true), tcons(true,true,false) {
      strand='.';
      num_mtcons=0;
      if (lfirst!=NULL) {
         add_Locus(lfirst);
         }
      id=0;
      }

    bool add_Locus(GLocus* loc) {
      if (mexons.Count()>0 && (end<loc->start || start > loc->end))
              return false; //no chance for overlapping exons
      if (mexons.Count()==0) {
          mexons.Add(loc->mexons);
          start=loc->start;
          end=loc->end;
          if (loc->qfidx<0) rloci.Add(loc);
                      else  qloci.Add(loc);
          strand=loc->mrna_maxcov->strand;
          loc->xlocus=this;
          return true;
          }
      int f=0;
      if (loc->qfidx<0) {
        if (rloci.Found(loc,f)) return false;
        }
      else if (qloci.Found(loc,f)) return false;

      // -- merge mexons
      GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
      int i=0; //index of first mexons with a merge
      int j=0; //index current mrna exon
      while (i<mexons.Count() && j<loc->mexons.Count()) {
          uint istart=mexons[i].start;
          uint iend=mexons[i].end;
          uint jstart=loc->mexons[j].start;
          uint jend=loc->mexons[j].end;
          if (iend<jstart) { i++; continue; }
          if (jend<istart) { j++; continue; }
          //if (mexons[i].overlap(jstart, jend)) {
          //exon overlap was found :
          ovlexons.Add(j);
          //extend mexons[i] as needed
          if (jstart<istart) mexons[i].start=jstart;
          if (jend>iend) { //mexons[i] end extend
              mexons[i].end=jend;
              //now this could overlap the next mexon(s), so we have to merge them all
              while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                  uint nextend=mexons[i+1].end;
                  mexons.Delete(i+1);
                  if (nextend>mexons[i].end) {
                      mexons[i].end=nextend;
                      break; //no need to check next mexons
                  }
              } //while next mexons merge
          } // mexons[i] end extend
          //  } //exon overlap
          j++; //check the next locus.mexon
        }//while mexons
      if (ovlexons.Count()==0) return false;
      if (strand=='.' && loc->mrna_maxcov->strand>0)
             strand=loc->mrna_maxcov->strand;
      //have exon overlap:
      //-- add the rest of the non-overlapping mexons:
       GSeg seg;
       for (int i=0;i<loc->mexons.Count();i++) {
            seg.start=loc->mexons[i].start;
            seg.end=loc->mexons[i].end;
            if (!ovlexons.Exists(i)) mexons.Add(seg);
            }
      // -- adjust start/end as needed
      if (start>loc->start) start=loc->start;
      if (end<loc->end) end=loc->end;
      loc->xlocus=this;
      if (loc->qfidx<0) rloci.Add(loc);
                  else  qloci.Add(loc);
      return true;
      }

  void addMerge(GXLocus& oxloc) {
    GArray<int> ovlexons(true,true); //list of oxloc.mexons indexes overlapping existing mexons
    int i=0; //index of first mexons with a merge
    int j=0; //index current mrna exon
    while (i<mexons.Count() && j<oxloc.mexons.Count()) {
        uint istart=mexons[i].start;
        uint iend=mexons[i].end;
        uint jstart=oxloc.mexons[j].start;
        uint jend=oxloc.mexons[j].end;
        if (iend<jstart) { i++; continue; }
        if (jend<istart) { j++; continue; }
        //if (mexons[i].overlap(jstart, jend)) {
        //exon overlap was found :
        ovlexons.Add(j);
        //extend mexons[i] as needed
        if (jstart<istart) mexons[i].start=jstart;
        if (jend>iend) { //mexons[i] end extend
            mexons[i].end=jend;
            //now this could overlap the next mexon(s), so we have to merge them all
            while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                uint nextend=mexons[i+1].end;
                mexons.Delete(i+1);
                if (nextend>mexons[i].end) {
                    mexons[i].end=nextend;
                    break; //no need to check next mexons
                }
            } //while next mexons merge
        } // mexons[i] end extend
        //  } //exon overlap
        j++; //check the next oxloc.mexon
    }
    if (ovlexons.Count()==0) {
      GError("Error: attempt to merge GXLoci with non-overlapping exons!\n");
      }
    //-- add the rest of the non-overlapping mexons:
    GSeg seg;
    for (int i=0;i<oxloc.mexons.Count();i++) {
        seg.start=oxloc.mexons[i].start;
        seg.end=oxloc.mexons[i].end;
        if (!ovlexons.Exists(i)) mexons.Add(seg);
        }
   if (start>oxloc.start) start=oxloc.start;
   if (end<oxloc.end) end=oxloc.end;
   if (strand=='.') strand=oxloc.strand;
    //-- steal all qloci and rloci
    for (int i=0;i<oxloc.qloci.Count();i++) {
         if (oxloc.qloci[i]->xlocus==this) continue;
         qloci.Add(oxloc.qloci[i]);
         oxloc.qloci[i]->xlocus=this;
         }
    for (int i=0;i<oxloc.rloci.Count();i++) {
         if (oxloc.rloci[i]->xlocus==this) continue;
         rloci.Add(oxloc.rloci[i]);
         oxloc.rloci[i]->xlocus=this;
         }
  }
};



int parse_mRNAs(GList<GffObj>& mrnas,
				 GList<GSeqData>& glstdata,
				 bool is_ref_set=true,
				 bool check_for_dups=false,
				 int qfidx=-1);

//reading a mRNAs from a gff file and grouping them into loci
void read_mRNAs(FILE* f, GList<GSeqData>& seqdata, GList<GSeqData>* ref_data=NULL, 
              bool check_for_dups=false, int qfidx=-1, const char* fname=NULL, bool checkseq=false);


bool tMatch(GffObj& a, GffObj& b, int& ovlen, bool equnspl=false);

//"position" a given coordinate x within a list of transcripts sorted by their start (lowest)
//coordinate, using quick-search; the returned int is the list index of the closest *higher*
//GffObj - i.e. starting right *ABOVE* the given coordinate
//Convention: returns -1 if there is no such GffObj (i.e. last GffObj starts below x)
int qsearch_mrnas(uint x, GList<GffObj>& mrnas);
int qsearch_segs(uint x, GList<GSeg>& segs); // same as above, but for GSeg lists

GSeqData* getRefData(int gid, GList<GSeqData>& ref_data); //returns reference GSeqData for a specific genomic sequence

#endif
