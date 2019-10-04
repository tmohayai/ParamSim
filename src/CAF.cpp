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

        cafMVA->Branch( "run", &run, "run/I" );
        cafMVA->Branch( "event", &event, "event/I" );
        cafMVA->Branch( "Ev[event]", &Ev, "Ev[event]/D" );
        cafMVA->Branch( "Ev_rec[event]", &Ev_rec, "Ev_rec[event]/D" );
        cafMVA->Branch( "gastpc_pi_pl_mult", &gastpc_pi_pl_mult, "gastpc_pi_pl_mult/I" );
        cafMVA->Branch( "gastpc_pi_min_mult", &gastpc_pi_min_mult, "gastpc_pi_min_mult/I" );

        cafMVA->Branch( "nFSP", &nFSP, "nFSP/I" );
        cafMVA->Branch( "truepdg", &truepdg, "truepdg[nFSP]/I" );
        cafMVA->Branch( "recopid", &recopid, "recopid[nFSP]/I" );
        cafMVA->Branch( "trkLen", &trkLen, "trkLen[nFSP]/D" );
        cafMVA->Branch( "trkLenPerp", &trkLenPerp, "trkLenPerp[nFSP]/D" );
        cafMVA->Branch( "pid_prob_arr", &pid_prob_arr, "pid_prob_arr[nFSP]/D" );
        cafMVA->Branch( "preco", &preco, "preco[nFSP]/D" );
        cafMVA->Branch( "anglereco", &anglereco, "anglereco[nFSP]/D" );

        cafMVA->Branch( "tpirpi_prob", &tpirpi_prob, "tpirpi_prob[nFSP]/D" );
        cafMVA->Branch( "tpirmu_prob", &tpirmu_prob, "tpirmu_prob[nFSP]/D" );

        cafMVA->Branch( "tpirp_prob", &tpirp_prob, "tpirp_prob[nFSP]/D" );
        cafMVA->Branch( "tpirk_prob", &tpirk_prob, "tpirk_prob[nFSP]/D" );

        cafMVA->Branch( "tpird_prob", &tpird_prob, "tpird_prob[nFSP]/D" );
        cafMVA->Branch( "tpire_prob", &tpire_prob, "tpire_prob[nFSP]/D" );

        cafMVA->Branch( "angle", &angle, "angle[nFSP]/D" );
        cafMVA->Branch( "erecon", &erecon, "erecon[nFSP]/D" );
        cafMVA->Branch( "truep", &truep, "truep[nFSP]/D" );
        cafMVA->Branch( "truepx", &truepx, "truepx[nFSP]/D" );
        cafMVA->Branch( "truepy", &truepy, "truepy[nFSP]/D" );
        cafMVA->Branch( "truepz", &truepz, "truepz[nFSP]/D" );
        cafMVA->Branch( "partEvReco", &partEvReco, "partEvReco[nFSP]/D" );
        cafMVA->Branch( "partereco", &partereco, "partereco[nFSP]/D" );
        cafMVA->Branch( "reco_nue", &reco_nue, "recon_nue/I" );
        cafMVA->Branch( "reco_numu", &reco_numu, "reco_numu/I" );
        cafMVA->Branch( "reco_nc", &reco_nc, "reco_nc/I" );

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
