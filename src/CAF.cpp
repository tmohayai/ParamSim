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

        //Event number
        cafMVA->Branch("Event", &Event);
        //Run number
        cafMVA->Branch("Run", &Run);
        //Sub Run number
        cafMVA->Branch("SubRun", &SubRun);

        //MC Truth info (Generator)
        cafMVA->Branch("mode", &mode);
        cafMVA->Branch("q2", &q2);
        cafMVA->Branch("w", &w);
        cafMVA->Branch("y", &y);
        cafMVA->Branch("x", &x);
        cafMVA->Branch("theta", &theta);
        cafMVA->Branch("t", &t);
        cafMVA->Branch("mctime", &mctime);
        cafMVA->Branch("ntype", &ntype);
        cafMVA->Branch("ccnc", &ccnc);
        cafMVA->Branch("gint", &gint);
        cafMVA->Branch("tgtpdg", &tgtpdg);
        cafMVA->Branch("weight", &weight);
        cafMVA->Branch("gt_t", &gt_t);
        cafMVA->Branch("intert", &intert);
        cafMVA->Branch("mcnupx", &mcnupx);
        cafMVA->Branch("mcnupy", &mcnupy);
        cafMVA->Branch("mcnupz", &mcnupz);
        cafMVA->Branch("vertx", &vertx);
        cafMVA->Branch("verty", &verty);
        cafMVA->Branch("vertz", &vertz);

        //Number of final state particle (primaries)
        cafMVA->Branch("nFSP", &nFSP);
        //MC Particle info
        cafMVA->Branch("mother", &mother);
        cafMVA->Branch("pdgmother", &pdgmother);
        cafMVA->Branch("MCPStartX", &MCPStartX);
        cafMVA->Branch("MCPStartY", &MCPStartY);
        cafMVA->Branch("MCPStartZ", &MCPStartZ);
        cafMVA->Branch("MCPStartPX", &MCPStartPX);
        cafMVA->Branch("MCPStartPY", &MCPStartPY);
        cafMVA->Branch("MCPStartPZ", &MCPStartPZ);
        cafMVA->Branch("MCPEndX", &MCPEndX);
        cafMVA->Branch("MCPEndY", &MCPEndY);
        cafMVA->Branch("MCPEndZ", &MCPEndZ);
        cafMVA->Branch("MCProc", &MCProc);
        cafMVA->Branch("MCEndProc", &MCEndProc);
        cafMVA->Branch("angle", &angle);
        cafMVA->Branch("truep", &truep);
        cafMVA->Branch("truepdg", &truepdg);
        //Reco info
        cafMVA->Branch("recopid", &recopid);
        cafMVA->Branch("trkLen", &trkLen);
        cafMVA->Branch("trkLenPerp", &trkLenPerp);
        cafMVA->Branch("preco", &preco);
        cafMVA->Branch("anglereco", &anglereco);
        cafMVA->Branch("erecon", &erecon);
        cafMVA->Branch("prob_arr", &prob_arr);

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

void CAF::ClearVectors()
{
    Run = 0;
    Event = 0;
    SubRun = 0;
    //Generator values
    mode.clear(); ccnc.clear(); ntype.clear(); gint.clear(); weight.clear(); tgtpdg.clear(); gt_t.clear(); intert.clear();
    q2.clear(); w.clear(); y.clear(); x.clear(); theta.clear(); t.clear(); mctime.clear(); mcnupx.clear(); mcnupy.clear(); mcnupz.clear(); vertx.clear(); verty.clear(); vertz.clear();
    //MC Particle Values
    nFSP.clear();
    mother.clear(); pdgmother.clear(); truepdg.clear(); MCPStartX.clear(); MCPStartY.clear(); MCPStartZ.clear(); MCPEndX.clear(); MCPEndY.clear(); MCPEndZ.clear(); MCPStartPX.clear(); MCPStartPY.clear(); MCPStartPZ.clear();
    MCProc.clear(); MCEndProc.clear();
    trkLen.clear(); trkLenPerp.clear(); truep.clear(); truepx.clear(); truepy.clear(); truepz.clear(); angle.clear();
    //Reco values
    recopid.clear();
    prob_arr.clear(); partereco.clear(); anglereco.clear(); preco.clear(); erecon.clear();
}
