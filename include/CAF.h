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

    void ClearVectors();

    //Need to be careful as they can be modified easily / better to make vectors!
    //Event-wise values
    int Run, Event, SubRun;
    //Generator values
    std::vector<int> mode, ccnc, ntype, gint, weight, tgtpdg, gt_t, intert;
    std::vector<double> q2, w, y, x, theta, t, mctime, mcnupx, mcnupy, mcnupz, vertx, verty, vertz;
    //MC Particle Values
    std::vector<int> nFSP;
    std::vector<int> mother, pdgmother, truepdg, MCPStartX, MCPStartY, MCPStartZ, MCPEndX, MCPEndY, MCPEndZ, MCPStartPX, MCPStartPY, MCPStartPZ;
    std::vector<std::string> MCProc, MCEndProc;
    std::vector<float> trkLen, trkLenPerp, truep, truepx, truepy, truepz, angle;
    //Reco values
    std::vector<int> recopid;
    std::vector<float> prob_arr, partereco, anglereco, preco, erecon;

private:
    std::string _outputFile;
    TFile * cafFile;
    TTree * cafMVA;
};

#endif
