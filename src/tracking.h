#include "abundances.h"

struct TrackingInfoPerRep
{
    boost::shared_ptr<const ReadGroupProperties> rg_props;
    double fpkm;
    double count;
    AbundanceStatus status;
};

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
    
    FPKM(r),
    FPKM_variance(v),
    FPKM_conf_lo(fcl),
    FPKM_conf_hi(fch),
    status(s),
    
    fpkm_samples(fs),
    gamma(g)
    {
        assert (fpr.size() == cpr.size());
        assert (fpr.size() == spr.size());
        assert (cpr.size() == spr.size());
        
        // TODO: should check for proper alignment of these tables...
        for (CountPerReplicateTable::const_iterator itr = cpr.begin(); itr != cpr.end(); ++itr)
        {
            TrackingInfoPerRep info;
            
            info.rg_props = itr->first;
            info.count = itr->second;
            
            FPKMPerReplicateTable::const_iterator f_itr = fpr.find(itr->first);
            if (f_itr != fpr.end())
                info.fpkm = f_itr->second;
            
            StatusPerReplicateTable::const_iterator s_itr = spr.find(itr->first);
            if (s_itr != spr.end())
                info.status = s_itr->second;
            
            tracking_info_per_rep.push_back(info);
        }
        
        vector<TrackingInfoPerRep>(tracking_info_per_rep).swap(tracking_info_per_rep);
    }
    
	double count_mean;
    double count_var;
    double count_uncertainty_var;
    double count_dispersion_var;
    vector<TrackingInfoPerRep> tracking_info_per_rep;
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
	
    vector<vector<boost::shared_ptr<const ReadGroupProperties> > > rg_props;
    
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
