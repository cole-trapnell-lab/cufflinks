#include "abundances.h"

struct FPKMContext
{
	FPKMContext(double cm,
                double cv,
                double cuv,
                double cdv,
                const CountPerReplicateTable& cpr,
                double r,
                const FPKMPerReplicateTable& fpr,
                double v,
                double fcl,
                double fch,
                AbundanceStatus s,
                const StatusPerReplicateTable& spr,
                const vector<double>& fs,
                double g)
    : count_mean(cm),
    count_var(cv),
    count_uncertainty_var(cuv),
    count_dispersion_var(cdv),
    count_per_rep(cpr),
    fpkm_per_rep(fpr),
    FPKM(r),
    FPKM_variance(v),
    FPKM_conf_lo(fcl),
    FPKM_conf_hi(fch),
    status(s),
    status_per_rep(spr),
    fpkm_samples(fs),
    gamma(g) {}
    
	double count_mean;
    double count_var;
    double count_uncertainty_var;
    double count_dispersion_var;
    CountPerReplicateTable count_per_rep;
    FPKMPerReplicateTable fpkm_per_rep;
    StatusPerReplicateTable status_per_rep;
	double FPKM;
	double FPKM_variance;
    double FPKM_conf_lo;
    double FPKM_conf_hi;
    AbundanceStatus status;
    vector<double> fpkm_samples;
    double gamma;
};

struct FPKMTracking
{
	string locus_tag;
	char classcode;
	set<string> tss_ids; // for individual isoforms only
    set<string> gene_ids;
	set<string> gene_names;
	set<string> protein_ids;
	string description; // isoforms or tss groups (e.g.) involved in this test
	string ref_match;
    int length;
	
	vector<FPKMContext> fpkm_series;
};

typedef map<string,  FPKMTracking> FPKMTrackingTable;

struct Tracking
{
	FPKMTrackingTable isoform_fpkm_tracking;
	FPKMTrackingTable tss_group_fpkm_tracking;
	FPKMTrackingTable gene_fpkm_tracking;
	FPKMTrackingTable cds_fpkm_tracking;
    
    void clear()
    {
        isoform_fpkm_tracking.clear();
        tss_group_fpkm_tracking.clear();
        gene_fpkm_tracking.clear();
        cds_fpkm_tracking.clear();
    }
};

void add_to_tracking_table(size_t sample_index,
                           Abundance& ab,
						   FPKMTrackingTable& track);
