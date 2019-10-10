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



    int nFSP , Run, verts, Event, SubRun, VertN, vert;
    int mode[10000], ccnc[10000], ntype[10000], gint[10000], weight[10000], tgtpdg[10000], gt_t[10000], mother[10000], pdgmother[10000], intert[10000];
    double q2[10000], w[10000], y[10000], x[10000], theta[10000], t[10000], mctime[10000], vertt, vertq, mcnupx[10000], mcnupy[10000], mcnupz[10000], vertx, verty, vertz;
    int truepdg[10000], recopid[10000], MCPStartX[10000], MCPStartY[10000], MCPStartZ[10000], MCPEndX[10000], MCPEndY[10000], MCPEndZ[10000], MCPStartPX[10000], MCPStartPY[10000], MCPStartPZ[10000];
    std::string MCProc[10000], MCEndProc[10000];
    double prob_arr[10000], trkLen[10000], trkLenPerp[10000];
    double partereco[10000], anglereco[10000], angle[10000], preco[10000], erecon[10000],truep[10000], truepx[10000], truepy[10000], truepz[10000];

private:
    std::string _outputFile;
    TFile * cafFile;
    TTree * cafMVA;
};

#endif
