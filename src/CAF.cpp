#include "CAF.h"

#include <iostream>

CAF::CAF()
: cafFile(nullptr), _outputFile("")
{

}

CAF::CAF( std::string filename )
: cafFile(nullptr), _outputFile(filename)
{

}

bool CAF::BookTFile()
{
    cafFile = new TFile(_outputFile.c_str(), "RECREATE");
    if(cafFile == nullptr){
        return false;
    }
    else
    {
        cafMVA = new TTree( "caf", "caf tree" );
        cafMVA->SetDirectory(0);

	cafMVA = new TTree( "caf", "caf tree" );
        cafMVA->SetDirectory(0);

        cafMVA->Branch( "Run", &Run, "Run/I" );
        cafMVA->Branch( "Event", &Event, "Event/I" );
        cafMVA->Branch( "SubRun", &SubRun, "SubRun/I" );
        cafMVA->Branch( "mode", &mode, "mode[Event]/I" );
        cafMVA->Branch( "q2", &q2, "q2[Event]/D" );
        cafMVA->Branch( "w", &w, "w[Event]/D" );
        cafMVA->Branch( "y", &y, "y[Event]/D" );
        cafMVA->Branch( "x", &x, "x[Event]/D" );
        cafMVA->Branch( "theta", &theta, "theta[Event]/D" );
        cafMVA->Branch( "t", &t, "t[Event]/D" );
        cafMVA->Branch( "mctime", &mctime, "mctime[Event]/D" );
        cafMVA->Branch( "vertt", &vertt, "vertt[Event]/D" );
        cafMVA->Branch( "vertq", &vertq, "vertq[Event]/D" );
        cafMVA->Branch( "ntype", &ntype, "ntype[Event]/I" );
        cafMVA->Branch( "ccnc", &ccnc, "ccnc[Event]/I" );
        cafMVA->Branch( "gint", &gint, "gint[Event]/I" );
        cafMVA->Branch( "tgtpdg", &tgtpdg, "tgtpdg[Event]/I" );
        cafMVA->Branch( "weight", &weight, "weight[Event]/I" );
        cafMVA->Branch( "gt_t", &gt_t, "gt_t[Event]/I" );
        cafMVA->Branch( "intert", &intert, "intert[Event]/I" );
        cafMVA->Branch( "VertN", &VertN, "VertN/I" );
        cafMVA->Branch( "verts", &verts, "verts/I" );
        cafMVA->Branch( "mcnupx", &mcnupx, "mcnupx[Event]/D" );
        cafMVA->Branch( "mcnupy", &mcnupy, "mcnupy[Event]/D" );
        cafMVA->Branch( "mcnupz", &mcnupz, "mcnupz[Event]/D" );

	cafMVA->Branch( "vertx", &vertx, "vertx[verts]/D" );
        cafMVA->Branch( "verty", &verty, "verty[verts]/D" );
        cafMVA->Branch( "vertz", &vertz, "vertz[verts]/D" );

	cafMVA->Branch( "nFSP", &nFSP, "nFSP/I" );
        cafMVA->Branch( "mother", &mother, "mother[nFSP]/I" );
        cafMVA->Branch( "pdgmother", &pdgmother, "pdgmother[nFSP]/I" );
        cafMVA->Branch( "MCPStartX", &MCPStartX, "MCPStartX[nFSP]/D" );
        cafMVA->Branch( "MCPStartY", &MCPStartY, "MCPStartY[nFSP]/D" );
        cafMVA->Branch( "MCPStartZ", &MCPStartZ, "MCPStartZ[nFSP]/D" );

        cafMVA->Branch( "MCPStartPX", &MCPStartPX, "MCPStartPX[nFSP]/D" );
        cafMVA->Branch( "MCPStartPY", &MCPStartPY, "MCPStartPY[nFSP]/D" );
        cafMVA->Branch( "MCPStartPZ", &MCPStartPZ, "MCPStartPZ[nFSP]/D" );

        cafMVA->Branch( "MCPEndX", &MCPEndX, "MCPEndX[nFSP]/D" );
        cafMVA->Branch( "MCPEndtY", &MCPEndY, "MCPEndY[nFSP]/D" );
        cafMVA->Branch( "MCPEndZ", &MCPEndZ, "MCPEndZ[nFSP]/D" );

        cafMVA->Branch( "MCProc", &MCProc, "MCProc[nFSP]/S" );
        cafMVA->Branch( "MCEndProc", &MCEndProc, "MCEndProc[nFSP]/S" );

        cafMVA->Branch( "truepdg", &truepdg, "truepdg[nFSP]/I" );
        cafMVA->Branch( "recopid", &recopid, "recopid[nFSP]/I" );
        cafMVA->Branch( "trkLen", &trkLen, "trkLen[nFSP]/D" );
        cafMVA->Branch( "trkLenPerp", &trkLenPerp, "trkLenPerp[nFSP]/D" );
        cafMVA->Branch( "preco", &preco, "preco[nFSP]/D" );
        cafMVA->Branch( "anglereco", &anglereco, "anglereco[nFSP]/D" );
        cafMVA->Branch( "angle", &angle, "angle[nFSP]/D" );
        cafMVA->Branch( "erecon", &erecon, "erecon[nFSP]/D" );
        cafMVA->Branch( "truep", &truep, "truep[nFSP]/D" );	
	cafMVA->Branch( "prob_arr", &prob_arr, "prob_arr[nFSP]/D" );
        return true;
    }
}

void CAF::CloseTFile()
{
    if(cafFile != nullptr) {
        cafFile->Close();
        return;
    }
    else{ return; }
}

void CAF::WriteTTree()
{
    if(cafMVA != nullptr) {
        cafFile->cd();
        cafMVA->Write();
        cafFile->Close();
    }

    return;
}

void CAF::FillTTree()
{
    if(cafMVA != nullptr) {
        cafFile->cd();
        cafMVA->Fill();
    }

    return;
}
