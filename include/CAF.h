#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"
// #include "Ntuple/NtpMCEventRecord.h"

class CAF {

public:

    CAF();

    CAF( std::string filename );

    CAF(const CAF &) = default;

    ~CAF();

    bool BookTFile();

    void FillTTree();

    void WriteTTree();

    void CloseTFile();

    // event accounting
    int run, event;
    // Truth information

    // Reco information CV
    double Ev_rec[10000], Ev[10000];
    int reco_numu, reco_nue, reco_nc;

    // Gas TPC variables
    int gastpc_pi_min_mult, gastpc_pi_pl_mult, true_numu_fhc_caf, true_numu_rhc_caf;
    int nFSP;
    int truepdg[10000], recopid[10000], tpirpi_prob[10000],tpirmu_prob[10000],tpirp_prob[10000],tpirk_prob[10000],tpird_prob[10000],tpire_prob[10000];
    double pid_prob_arr[10000], trkLen[10000], trkLenPerp[10000], partEvReco[10000], mureco[10000];
    double partereco[10000], anglereco[10000], angle[10000], preco[10000], erecon[10000],truep[10000], truepx[10000], truepy[10000], truepz[10000];

private:
    std::string _outputFile;
    TFile * cafFile;
    TTree * cafMVA;
};

#endif
