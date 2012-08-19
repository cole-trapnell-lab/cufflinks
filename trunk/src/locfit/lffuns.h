


/* FILES IN THE src DIRECTORY */

/* adap.c */
extern double afit(), aband2(), aband3();
extern INT ainitband();

/* band.c */
extern void band(), kdeselect();

/* density.c */
extern INT densinit();
extern INT fact[];
extern int likeden();
extern void prodint_resp(), prresp();

/* dens_haz.c */
extern void haz_init();
extern INT hazint();

/* dens_int.c */
extern double dens_integrate();
extern void dens_renorm(), dens_lscv(), lforder();

/* dist.c */
extern double igamma(), ibeta();
extern double pf(), pchisq(), pnorm();
extern double df(), dchisq();

/* ev_atree.c */
extern void atree_start(), atree_grow(), atree_guessnv();
extern double atree_int();

/* ev_interp.c */
extern double dointpoint(), cubintd();
extern double linear_interp(), cubic_interp(), rectcell_interp();
extern INT exvval();
extern void exvvalpv(), hermite2();

/* ev_kdtre.c */
extern void kdtre_start();
extern double kdtre_int();

/* ev_main.c */
extern void trchck(), guessnv();
extern void dataf(), gridf(), crossf(), xbarf(), preset();
extern INT newsplit();
extern int lfit_reqd(), lfit_reqi();
#ifndef CVERSION
extern vari *createvar();
#endif

/* ev_trian.c */
extern void triang_start(), triang_grow();
extern double triang_int();

/* family.c */
extern INT links(), stdlinks(), defaultlink(), validlinks();
extern double b2(), b3(), b4(), lf_link(), invlink();

/* frend.c */
extern void fitfun(), degfree(), ressumm(), makecfn();
extern INT procv(), procvraw(), procvvord();
extern double base(), cens(), prwt(), resp(), getxi(), rss();
extern INT calcp();

/* kappa0.c */
extern double critval(), critvalc(), tailp(), taild();
extern INT constants();

/* lf_dercor.c */
extern void dercor();

/* lf_fitfun.c */
extern void fitfun(), designmatrix();
extern INT calcp(), coefnumber();

/* lf_robust.c */
extern double median();
extern void lf_robust();

/* lfstr.c */
extern void setstrval();
extern INT ppwhat(), restyp();

/* lf_vari.c */
extern void comp_vari(), local_df();
extern double comp_infl();

/* linalg.c */
extern void svd(), hsvdsolve();
extern void addouter(), multmatscal();
extern void QRupd(), QR1(), bacK(), bacT(), solve(), grsc();
extern void setzero(), unitvec();
extern void transpose();
extern double innerprod(), m_trace();
extern INT svdsolve();

/* locfit.c or parfit.c (most) */
extern int ident, locfit(), lf_iter();

/* math.c */
extern double lflgamma(), lferf(), lferfc(), lfdaws();
extern double ptail(), logit(), expit();
//extern double lgamma(), erf(), erfc();
extern int factorial();

/* minmax.c */
extern double ipower(), minmax();

/* nbhd.c */
extern double kordstat(), nbhd(), rho();

/* odint.c */
extern INT onedint();
extern void recurint();

/* pcomp.c */
extern double addparcomp();
extern void compparcomp(), subparcomp(), subparcomp2(), pcchk();
extern int pc_reqd();
extern INT noparcomp(), hasparcomp();

/* preplot.c */
extern void preplot(), cpreplot();
extern INT setpppoints();

/* resid.c */
extern double resid();
extern void cfitted();
extern vari *vfitted(), *vresid();

/* scb.c */
extern void scb(), cscbsim();

/* simul.c */
extern void liksim(), scbsim(), scbmax(), regband(), rband();

/* startlf.c */
extern void bbox(), deschk(), startlf(), preproc(), fitdefault();
extern void fitoptions(), clocfit(), endfit();
extern INT nofit();

/* strings.c */
extern int stm(), pmatch(), matchlf(), matchrt(), checkltor(), checkrtol();
extern void strip();

/* wdiag.c */
extern INT wdiag(), procvhatm();
extern void cwdiag();

/* weight.c */
extern double W(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern INT iscompact(), wtaylor();

/* arith.c */
extern INT arvect(), intitem();
extern double areval(), arith(), darith(), dareval();
extern vari *varith(), *saveresult(), *arbuild();

/* c_args.c */
#define argused(v,i) (((carg *)viptr(v,i))->used)
#define setused(v,i) { ((carg *)viptr(v,i))->used = 1; }
#define setunused(v,i) { ((carg *)viptr(v,i))->used = 0; }
#define argarg(v,i) (((carg *)viptr(v,i))->arg)
#define argvalis(v,i,z) (strcmp(argval(v,i),z)==0)
extern char *argval(), *getargval();
extern int getarg(), readilist(), getlogic();

/* cmd.c */
extern int locfit_dispatch(char*);
extern void setuplf(), recondat(), cmdint();
extern double backtr(), docrit();

/* c_plot.c */
extern void plotdata(), plotfit(), plottrack(), plotopt(), setplot();

/* help.c */
extern void example();

/* lfd.c */
extern void doreaddata(), dosavedata(), dosavefit();
extern INT  setfilename();

/* main.c */
extern void SetWinDev();

/* makecmd.c */
extern vari *getcmd();
extern void makecmd(), del_clines(), inc_forvar(), dec_forvar();

/* post.c */
extern void SetPSDev();

/* pout.c */
extern INT pretty();
extern void displayplot();
extern void plotmaple(), plotmathe(), plotmatlb(), plotgnup(), plotxwin();

/* random.c */
extern double rnorm(), rexp(), runif(), rpois();
extern void rseed();

/* readfile.c */
extern void readfile();

/* scbmax.c */
extern void cscbmax();

/* vari.c */
#include "vari.hpp"

