#include "CAF.h"

#include "TRandom3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <functional>

CAF::CAF()
: cafFile(nullptr), _outputFile(""), _inputfile(""), _intfile(nullptr), _inttree(nullptr)
{

}

CAF::CAF( std::string infile, std::string filename )
: cafFile(nullptr), _outputFile(filename), _inputfile(infile), _intfile(nullptr), _inttree(nullptr)
{

}

bool CAF::BookTFile()
{
    _intfile = new TFile( _inputfile.c_str() );
    if(nullptr == _intfile){
        std::cout << "Cannot open file " << _inputfile.c_str() << std::endl;
        return false;
    }

    _inttree = (TTree*) _intfile->Get( "anatree/GArAnaTree" );
    if(nullptr == _inttree){
        std::cout << "Cannot find tree anatree/GArAnaTree" << std::endl;
        return false;
    }

    cafFile = new TFile(_outputFile.c_str(), "RECREATE");
    if(nullptr == cafFile){
        return false;
    }
    else
    {
        cafMVA = new TTree( "caf", "caf tree" );
        cafMVA->SetDirectory(0);

        //Event number
        cafMVA->Branch("Event", &_Event);
        //Run number
        cafMVA->Branch("Run", &_Run);
        //Sub Run number
        cafMVA->Branch("SubRun", &_SubRun);

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
        cafMVA->Branch("nFSP", &_nFSP);
        //MC Particle info
        cafMVA->Branch("mother", &mother);
        cafMVA->Branch("pdgmother", &pdgmother);
        cafMVA->Branch("MCPStartX", &_MCPStartX);
        cafMVA->Branch("MCPStartY", &_MCPStartY);
        cafMVA->Branch("MCPStartZ", &_MCPStartZ);
        cafMVA->Branch("MCPStartPX", &_MCPStartPX);
        cafMVA->Branch("MCPStartPY", &_MCPStartPY);
        cafMVA->Branch("MCPStartPZ", &_MCPStartPZ);
        cafMVA->Branch("MCPEndX", &_MCPEndX);
        cafMVA->Branch("MCPEndY", &_MCPEndY);
        cafMVA->Branch("MCPEndZ", &_MCPEndZ);
        cafMVA->Branch("MCProc", &_MCProc);
        cafMVA->Branch("MCEndProc", &_MCEndProc);
        cafMVA->Branch("angle", &_angle);
        cafMVA->Branch("truep", &truep);
        cafMVA->Branch("truepdg", &truepdg);
        //Reco info
        cafMVA->Branch("recopid", &recopid);
        cafMVA->Branch("trkLen", &trkLen);
        cafMVA->Branch("trkLenPerp", &trkLenPerp);
        cafMVA->Branch("preco", &_preco);
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
    //Generator values
    mode.clear();
    ccnc.clear();
    ntype.clear();
    gint.clear();
    weight.clear();
    tgtpdg.clear();
    gt_t.clear();
    intert.clear();
    q2.clear();
    w.clear();
    y.clear();
    x.clear();
    theta.clear();
    t.clear();
    mctime.clear();
    mcnupx.clear();
    mcnupy.clear();
    mcnupz.clear();
    vertx.clear();
    verty.clear();
    vertz.clear();

    //MC Particle Values
    _nFSP.clear();
    mother.clear();
    pdgmother.clear();
    truepdg.clear();
    _MCPStartX.clear();
    _MCPStartY.clear();
    _MCPStartZ.clear();
    _MCPEndX.clear();
    _MCPEndY.clear();
    _MCPEndZ.clear();
    _MCPStartPX.clear();
    _MCPStartPY.clear();
    _MCPStartPZ.clear();
    _MCProc.clear();
    _MCEndProc.clear();
    trkLen.clear();
    trkLenPerp.clear();
    truep.clear();
    truepx.clear();
    truepy.clear();
    truepz.clear();
    _angle.clear();

    //Reco values
    recopid.clear();
    prob_arr.clear();
    partereco.clear();
    anglereco.clear();
    _preco.clear();
    erecon.clear();
}

// main loop function
void CAF::loop()
{
    //double gastpc_len = 5.; // track length cut in cm
    float gastpc_len = 2.; // new track length cut in cm based on Thomas' study of low energy protons
    // dont care about electrons -- check momentum and see if hit ECAL
    float gastpc_B = 0.5; // B field strength in Tesla
    float gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
    float gastpc_X0 = 1300.; // cm = 13m radiation length

    //Resolution for short tracks //TODO check this numbers!
    float sigmaP_short = 0.1; //in GeV
    float sigma_angle_short = 0.01; //in radian

    // point resolution
    float sigma_x = 0.1;

    int seed = 7;
    TRandom3 *rando = new TRandom3( seed );

    std::vector<float> v = std::vector<float>();
    for (float pit = 0.040; pit < 20.0; pit += 0.001)
    {
        v.push_back(pit);
    }

    TRandom *R = new TRandom(time(0));

    //TODO as function of momentum
    float NeutronECAL_detEff = 0.4;
    float sigmaNeutronECAL_first = 0.11;
    float sigmaNeutronECAL_rescatter = 0.26;

    //ECAL energy resolution sigmaE/E
    float ECAL_stock = 0.06; //in %
    float ECAL_const = 0.02;
    TF1 *fRes = new TF1("fRes", "TMath::Sqrt ( [0]*[0]/x + [1]*[1] )", 3);
    fRes->FixParameter(0, ECAL_stock);
    fRes->FixParameter(1, ECAL_const);

    float ECAL_pi0_resolution = 0.13; //sigmaE/E in between at rest (17%) and high energy (~few %)

    TParticlePDG *neutron = TDatabasePDG::Instance()->GetParticle(2112);
    float neutron_mass = neutron->Mass(); //in GeV
    TParticlePDG *pi0 = TDatabasePDG::Instance()->GetParticle(111);
    float pi0_mass = pi0->Mass(); //in GeV

    //------------------------------------------------------------------------

    int           	 Event = 0;
    int           	 SubRun = 0;
    int           	 Run = 0;

    std::vector<float> *MC_Q2 = 0;
    std::vector<float> *MC_W = 0;
    std::vector<float> *MC_Y = 0;
    std::vector<float> *MC_X = 0;
    std::vector<float> *MC_Theta = 0;
    std::vector<float> *MC_T = 0;
    std::vector<float> *MCVertX = 0;
    std::vector<float> *MCVertY = 0;
    std::vector<float> *MCVertZ = 0;
    std::vector<float> *MCNuPx = 0;
    std::vector<float> *MCNuPy = 0;
    std::vector<float> *MCNuPz = 0;
    std::vector<int>   *NType = 0;
    std::vector<int>   *CCNC = 0;
    std::vector<int>   *Mode = 0;
    std::vector<int>   *Gint=0;
    std::vector<int>   *TgtPDG=0;
    std::vector<int>   *GT_T=0;
    std::vector<int>   *InterT=0;
    std::vector<float> *Weight=0;

    std::vector<int>     *PDG = 0;
    std::vector<int>     *Mother=0;
    std::vector<int>     *PDGMother=0;
    std::vector<int>     *MCPTrkID = 0;
    std::vector<float>   *MCPStartX = 0;
    std::vector<float>   *MCPStartY = 0;
    std::vector<float>   *MCPStartZ = 0;
    std::vector<float>   *MCPStartPX = 0;
    std::vector<float>   *MCPStartPY = 0;
    std::vector<float>   *MCPStartPZ = 0;
    std::vector<std::string>   *MCPProc = 0;
    std::vector<std::string>   *MCPEndProc = 0;
    std::vector<float>   *MCPEndX = 0;
    std::vector<float>   *MCPEndY = 0;
    std::vector<float>   *MCPEndZ = 0;
    std::vector<float>   *TrajMCPX = 0;
    std::vector<float>   *TrajMCPY = 0;
    std::vector<float>   *TrajMCPZ = 0;
    std::vector<int>     *TrajMCPTrajIndex = 0;

    _inttree->SetBranchAddress("Event", &Event);
    _inttree->SetBranchAddress("SubRun", &SubRun);
    _inttree->SetBranchAddress("Run", &Run);

    //Generator info
    _inttree->SetBranchAddress("NType", &NType);
    _inttree->SetBranchAddress("CCNC", &CCNC);
    _inttree->SetBranchAddress("MC_Q2", &MC_Q2);
    _inttree->SetBranchAddress("MC_W", &MC_W);
    _inttree->SetBranchAddress("MC_Y", &MC_Y);
    _inttree->SetBranchAddress("MC_X", &MC_X);
    _inttree->SetBranchAddress("MC_Theta", &MC_Theta);
    _inttree->SetBranchAddress("MC_T", &MC_T);
    _inttree->SetBranchAddress("Mode", &Mode);
    _inttree->SetBranchAddress("Gint", &Gint);
    _inttree->SetBranchAddress("TgtPDG", &TgtPDG);
    _inttree->SetBranchAddress("GT_T", &GT_T);
    _inttree->SetBranchAddress("MCVertX", &MCVertX);
    _inttree->SetBranchAddress("MCVertY", &MCVertY);
    _inttree->SetBranchAddress("MCVertZ", &MCVertZ);
    _inttree->SetBranchAddress("MCNuPx", &MCNuPx);
    _inttree->SetBranchAddress("MCNuPy", &MCNuPy);
    _inttree->SetBranchAddress("MCNuPz", &MCNuPz);
    _inttree->SetBranchAddress("InterT", &InterT);
    _inttree->SetBranchAddress("Weight", &Weight);

    //MC info
    _inttree->SetBranchAddress("PDG", &PDG);
    _inttree->SetBranchAddress("MCPTrkID", &MCPTrkID);
    _inttree->SetBranchAddress("MCPStartX", &MCPStartX);
    _inttree->SetBranchAddress("MCPStartY", &MCPStartY);
    _inttree->SetBranchAddress("MCPStartZ", &MCPStartZ);
    _inttree->SetBranchAddress("MCPEndX", &MCPEndX);
    _inttree->SetBranchAddress("MCPEndY", &MCPEndY);
    _inttree->SetBranchAddress("MCPEndZ", &MCPEndZ);
    _inttree->SetBranchAddress("Mother", &Mother);
    _inttree->SetBranchAddress("PDGMother", &PDGMother);
    _inttree->SetBranchAddress("MCPStartPX", &MCPStartPX);
    _inttree->SetBranchAddress("MCPStartPY", &MCPStartPY);
    _inttree->SetBranchAddress("MCPStartPZ", &MCPStartPZ);
    _inttree->SetBranchAddress("MCPProc", &MCPProc);
    _inttree->SetBranchAddress("MCPEndProc", &MCPEndProc);
    _inttree->SetBranchAddress("TrajMCPX", &TrajMCPX);
    _inttree->SetBranchAddress("TrajMCPY", &TrajMCPY);
    _inttree->SetBranchAddress("TrajMCPZ", &TrajMCPZ);
    _inttree->SetBranchAddress("TrajMCPTrajIndex", &TrajMCPTrajIndex);

    //gamma, neutron, pi0, k0L, k0S, k0, delta0
    std::vector<int> pdg_neutral = {22, 2112, 111, 130, 310, 311, 2114};

    //-------------------------------------------------------------------

    // Main event loop
    for( int entry = 0; entry < _inttree->GetEntries(); entry++ )
    {
        _inttree->GetEntry(entry);

        this->ClearVectors();

        if(Event%100 == 0)
        std::cout << "Treating Evt " << Event << std::endl;

        //Filling MCTruth values
        // std::cout << "Event " << Event << " Run " << Run << std::endl;
        _Event = Event;
        _Run = Run;
        _SubRun = SubRun;

        for(size_t i = 0; i < NType->size(); i++)
        {
            ccnc.push_back(CCNC->at(i));
            ntype.push_back(NType->at(i));
            q2.push_back(MC_Q2->at(i));
            w.push_back(MC_W->at(i));
            y.push_back(MC_Y->at(i));
            x.push_back(MC_X->at(i));
            theta.push_back(MC_Theta->at(i));
            mode.push_back(Mode->at(i));
            t.push_back(MC_T->at(i));
            intert.push_back(InterT->at(i));
            vertx.push_back(MCVertX->at(i));
            verty.push_back(MCVertY->at(i));
            vertz.push_back(MCVertZ->at(i));
            mcnupx.push_back(MCNuPx->at(i));
            mcnupy.push_back(MCNuPy->at(i));
            mcnupz.push_back(MCNuPz->at(i));
        }

        for(size_t i = 0; i < Gint->size(); i++)
        {
            gint.push_back(Gint->at(i));
            tgtpdg.push_back(TgtPDG->at(i));
            gt_t.push_back(GT_T->at(i));
            weight.push_back(Weight->at(i));
        }

        //--------------------------------------------------------------------------
        // Start of Parameterized Reconstruction
        //--------------------------------------------------------------------------
        int nFSP = 0;

        //TODO we should skip particles that we don't see or cannot reconstruct! change filling caf with i to index counting the particle
        //have to be careful with indexes and continue
        //---------------------------------------------------------------
        // all Gluckstern calculations happen in the following loop
        for(size_t i = 0; i < MCPStartPX->size(); ++i )
        {
            //Get the creating process
            std::string mcp_process = MCPProc->at(i);
            //Get ending process
            std::string mcp_endprocess = MCPEndProc->at(i);

            nFSP++;

            int pdg = PDG->at(i);
            //need to ignore neutrals for this - put the value to 0
            auto result = std::find(pdg_neutral.begin(), pdg_neutral.end(), pdg);
            bool isNeutral = (result == pdg_neutral.end()) ? false : true;

            if( isNeutral )
            {
                trkLen.push_back(0.);
                trkLenPerp.push_back(0.);
            }

            // calculate the total and the transverse track lengths and restrict the
            // tracklength to be above the gas TPC track length threshold
            double tracklen = 0.;
            double tracklen_perp = 0.;

            //TODO check if the mcp point is within the TPC volume! Skip for mcp in the ECAL (showers)
            //TODO Link showers to original mcp?
            for(int itraj = 1; itraj < TrajMCPX->size(); itraj++){
                //check that it is the correct mcp
                if(TrajMCPTrajIndex->at(itraj) == i){
                    //Traj point+1
                    TVector3 point(TrajMCPX->at(itraj), TrajMCPY->at(itraj), TrajMCPZ->at(itraj));
                    //point is not in the TPC anymore - stop traj loop
                    if(not isInTPC(point)){
                        // std::cout << "Point not within the TPC: " << point.X() << " r " << std::sqrt(point.Y()*point.Y() + point.Z()*point.Z()) << std::endl;
                        continue;
                    }
                    // find the length of the track by getting the distance between each hit
                    TVector3 diff(TrajMCPX->at(itraj) - TrajMCPX->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1), TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1));
                    // perp length
                    TVector2 tracklen_perp_vec(TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1));
                    // Summing up
                    tracklen += diff.Mag();
                    tracklen_perp += tracklen_perp_vec.Mod();
                }
            }

            trkLen.push_back(tracklen);
            trkLenPerp.push_back(tracklen_perp);

            //----------------------------------

            TVector3 mcp(MCPStartPX->at(i), MCPStartPY->at(i), MCPStartPZ->at(i));
            TVector3 xhat(1, 0, 0);

            float ptrue = (mcp).Mag();
            float pz = mcp.Z();
            float pt = (mcp.Cross(xhat)).Mag();
            float px = mcp.X();
            float py = mcp.Y();
            float mctrackid = MCPTrkID->at(i);
            // angle with respect to the incoming neutrino
            float angle  = atan(mcp.X() / mcp.Z());

            //for neutrons
            if(pdg == 2112)
            {
                //check if it can be detected by the ECAL
                //Assumes 40% efficiency to detect
                float random_number = R->Rndm();
                if(random_number >= NeutronECAL_detEff)
                {
                    //TODO random is first interaction of rescatter and smear accordingly to Chris's study
                    //Detected in the ECAL
                    recopid.push_back(2112);
                    float ereco = rando->Gaus( std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass), sigmaNeutronECAL_first );
                    erecon.push_back(ereco > 0 ? ereco : 0.);
                    // std::cout << "true part n true energy " << std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass) << " ereco " << erecon[i] << std::endl;
                    truepdg.push_back(pdg);
                    truepx.push_back(MCPStartPX->at(i));
                    truepy.push_back(MCPStartPY->at(i));
                    truepz.push_back(MCPStartPZ->at(i));
                    _MCPStartX.push_back(MCPStartX->at(i));
                    _MCPStartY.push_back(MCPStartY->at(i));
                    _MCPStartZ.push_back(MCPStartZ->at(i));
                    _MCPEndX.push_back(MCPEndX->at(i));
                    _MCPEndY.push_back(MCPEndY->at(i));
                    _MCPEndZ.push_back(MCPEndZ->at(i));
                    mother.push_back(Mother->at(i));
                    pdgmother.push_back(PDGMother->at(i));
                    // save the true momentum
                    truep.push_back(ptrue);
                    // save the true angle
                    _angle.push_back(angle);
                    //Save MC process
                    _MCProc.push_back(mcp_process);
                    _MCEndProc.push_back(mcp_endprocess);
                }
            }

            //for pi0s
            if(pdg == 111)
            {
                //TODO smear the pi0 energy (and decay vertex?) according to previous pi0 reco studies
                float ereco = rando->Gaus( std::sqrt(ptrue*ptrue + pi0_mass*pi0_mass), ECAL_pi0_resolution*std::sqrt(ptrue*ptrue + pi0_mass*pi0_mass));
                erecon.push_back(ereco);
                recopid.push_back(111);

                truepdg.push_back(pdg);
                truepx.push_back(MCPStartPX->at(i));
                truepy.push_back(MCPStartPY->at(i));
                truepz.push_back(MCPStartPZ->at(i));
                _MCPStartX.push_back(MCPStartX->at(i));
                _MCPStartY.push_back(MCPStartY->at(i));
                _MCPStartZ.push_back(MCPStartZ->at(i));
                _MCPEndX.push_back(MCPEndX->at(i));
                _MCPEndY.push_back(MCPEndY->at(i));
                _MCPEndZ.push_back(MCPEndZ->at(i));
                mother.push_back(Mother->at(i));
                pdgmother.push_back(PDGMother->at(i));
                // save the true momentum
                truep.push_back(ptrue);
                // save the true angle
                _angle.push_back(angle);
                //Save MC process
                _MCProc.push_back(mcp_process);
                _MCEndProc.push_back(mcp_endprocess);
            }

            //for gammas
            if(pdg == 22)
            {
                //TODO check if they are not from a pi0 or decayed in the TPC and hit the ECAL!
                if( PDGMother->at(i) != 111 )
                {
                    TVector3 point(MCPEndX->at(i), MCPEndY->at(i), MCPEndZ->at(i));
                    //Endpoint is not in the TPC
                    if(not isInTPC(point)){
                        //if they hit the ECAL and smear their energy
                        float ECAL_resolution = fRes->Eval(ptrue)*ptrue;
                        float ereco = rando->Gaus(ptrue, ECAL_resolution);
                        erecon.push_back(ereco);
                        recopid.push_back(22);

                        truepdg.push_back(pdg);
                        truepx.push_back(MCPStartPX->at(i));
                        truepy.push_back(MCPStartPY->at(i));
                        truepz.push_back(MCPStartPZ->at(i));
                        _MCPStartX.push_back(MCPStartX->at(i));
                        _MCPStartY.push_back(MCPStartY->at(i));
                        _MCPStartZ.push_back(MCPStartZ->at(i));
                        _MCPEndX.push_back(MCPEndX->at(i));
                        _MCPEndY.push_back(MCPEndY->at(i));
                        _MCPEndZ.push_back(MCPEndZ->at(i));
                        mother.push_back(Mother->at(i));
                        pdgmother.push_back(PDGMother->at(i));
                        // save the true momentum
                        truep.push_back(ptrue);
                        // save the true angle
                        _angle.push_back(angle);
                        //Save MC process
                        _MCProc.push_back(mcp_process);
                        _MCEndProc.push_back(mcp_endprocess);
                    }
                }
            }

            //Visible in the TPC
            if( trkLen.at(i) > gastpc_len )
            {
                //Use range instead of Gluckstern for stopping tracks
                //TODO is that correct? What if it is a scatter in the TPC? Need to check if daughter is same particle
                float preco = 0;
                TVector3 point(MCPEndX->at(i), MCPEndY->at(i), MCPEndZ->at(i));

                // save the true PDG, parametrized PID comes later
                truepdg.push_back(pdg);
                truepx.push_back(MCPStartPX->at(i));
                truepy.push_back(MCPStartPY->at(i));
                truepz.push_back(MCPStartPZ->at(i));
                _MCPStartX.push_back(MCPStartX->at(i));
                _MCPStartY.push_back(MCPStartY->at(i));
                _MCPStartZ.push_back(MCPStartZ->at(i));
                _MCPEndX.push_back(MCPEndX->at(i));
                _MCPEndY.push_back(MCPEndY->at(i));
                _MCPEndZ.push_back(MCPEndZ->at(i));
                mother.push_back(Mother->at(i));
                pdgmother.push_back(PDGMother->at(i));
                // save the true momentum
                truep.push_back(ptrue);
                // save the true angle
                _angle.push_back(angle);
                //Save MC process
                _MCProc.push_back(mcp_process);
                _MCEndProc.push_back(mcp_endprocess);

                if(isInTPC(point))
                {
                    preco = rando->Gaus( ptrue, sigmaP_short );
                    float angle_reco = rando->Gaus(angle, sigma_angle_short);
                    _preco.push_back(preco);
                    anglereco.push_back(angle_reco);
                }
                else{
                    // calculate number of trackpoints
                    float nHits = round (trkLen.at(i) / gastpc_padPitch);
                    // measurement term in Gluckstern formula
                    float fracSig_meas = sqrt(720./(nHits+4)) * ((0.01*gastpc_padPitch*ptrue) / (0.3 * gastpc_B * 0.0001 *trkLenPerp.at(i)*trkLenPerp.at(i)));
                    // multiple Coulomb scattering term in Gluckstern formula
                    float fracSig_MCS = (0.052*sqrt(1.43)) / (gastpc_B * sqrt(gastpc_X0*trkLenPerp.at(i)*0.0001));
                    // momentum resoltion from the two terms above
                    float sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
                    // now Gaussian smear the true momentum using the momentum resolution
                    preco = rando->Gaus( ptrue, sigmaP );

                    // measurement term in the Gluckstern formula for calculating the
                    // angular resolution
                    float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / trkLen.at(i)*trkLen.at(i)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
                    // scattering term in Gluckstern formula
                    float sigma_angle_2 = (0.015*0.015 / (3. * ptrue * ptrue)) * (trkLen.at(i)/gastpc_X0);
                    // angular resolution from the two terms above
                    float sigma_angle = sqrt(sigma_angle_1 + sigma_angle_2);
                    // now Gaussian smear the true angle using the angular resolution
                    float angle_reco = rando->Gaus(angle, sigma_angle);

                    // save reconstructed momentum and angle to cafanatree
                    _preco.push_back(preco);
                    anglereco.push_back(angle_reco);
                }

                //--------------------------------------------------------------------------
                // Start of PID Parametrization
                //--------------------------------------------------------------------------

                float p = preco;
                // read the PID parametrization ntuple from T. Junk
                TString filename="pid.root";
                TFile infile(filename,"READ");

                char str[10];
                std::vector<double> vec;
                std::vector<int> pdglist = {2112, 211, 13, 2212, 321, 100001020, 11};
                std::vector<std::string> pnamelist     = {"n", "#pi", "#mu", "p", "K", "d", "e"};
                std::vector<std::string> recopnamelist = {"n", "#pi", "#mu", "p", "K", "d", "e"};

                int qclosest = 0;
                float dist = 99990.;
                for (int q = 0; q < 501; ++q)
                {
                    if (v[q] < p < v[q+1])
                    {
                        sprintf(str, "%d", q);
                        std::string s = "pidmatrix";
                        s.append(str);
                        // read the 500 histograms one by one; each histogram is a
                        // 6 by 6 matrix of probabilities for a given momentum value
                        TH2F *pidinterp = (TH2F*) infile.Get(s.c_str())->Clone("pidinterp");
                        //Check the title and the reco momentum take only the one that fits
                        std::string fulltitle = pidinterp->GetTitle();
                        unsigned first = fulltitle.find("=");
                        unsigned last = fulltitle.find("GeV");
                        std::string substr = fulltitle.substr(first+1, last - first-1);
                        float pidinterp_mom = std::atof(substr.c_str());
                        //std::cout << pidinterp_mom << std::endl;
                        //calculate the distance between the bin and mom, store the q the closest
                        float disttemp = std::abs(pidinterp_mom - p);
                        //std::cout << disttemp << " " << dist << std::endl;
                        if( disttemp < dist ){
                            dist = disttemp;
                            qclosest = q;
                            // std::cout << "pid mom " << pidinterp_mom << " reco mom " << p << " dist " << dist << " qclosest " << qclosest << std::endl;
                        }
                    } // closes the pid vector loop
                } // closes the "pidmatrix" loop

                // std::cout << qclosest << std::endl;
                sprintf(str, "%d", qclosest);
                std::string mtx = "pidmatrix";
                mtx.append(str);
                // std::cout << mtx << std::endl;
                TH2F *pidinterp = (TH2F*) infile.Get(mtx.c_str())->Clone("pidinterp");

                //loop over the columns (true pid)
                for (int pidm = 1; pidm <= 6; ++pidm)
                {
                    //if it is in the pid table at index pidm
                    std::vector< P > v_prob;
                    if ( abs(pdg) == pdglist[pidm] )
                    {
                        //get true particle name
                        std::string trueparticlename = pidinterp->GetXaxis()->GetBinLabel(pidm);
                        if ( trueparticlename == pnamelist[pidm] )
                        {
                            //loop over the rows (reco pid)
                            for (int pidr=1; pidr <= 6; ++pidr)
                            {
                                std::string recoparticlename = pidinterp->GetYaxis()->GetBinLabel(pidr);
                                if (recoparticlename == recopnamelist[pidr])
                                {
                                    float prob = pidinterp->GetBinContent(pidm,pidr);
                                    prob_arr.push_back(prob);
                                    // std::cout << "true part " << trueparticlename << " true pid " << pdglist[pidm] << " reco name " << recoparticlename << " reco part list "
                                    // << recopnamelist[pidr] <<  " true mom " << ptrue << " reco mom " <<  p << " prob " << pidinterp->GetBinContent(pidm,pidr) << '\n';
                                    //Need to check random number value and prob value then associate the recopdg to the reco prob
                                    v_prob.push_back( std::make_pair(prob, recoparticlename) );
                                }
                            }

                            if(v_prob.size() > 1){
                                //Order the vector of prob
                                std::sort(v_prob.begin(), v_prob.end());
                                //Throw a random number between 0 and 1
                                float random_number = R->Rndm();
                                //Make cumulative sum to get the range
                                std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});
                                // std::cout << "rand " << random_number << std::endl;
                                // for(int ivec = 0; ivec < v_prob.size(); ivec++)
                                // {
                                //     std::cout << "Cumulative prob " << v_prob.at(ivec).first << " particle " << v_prob.at(ivec).second << std::endl;
                                // }
                                for(int ivec = 0; ivec < v_prob.size()-1; ivec++)
                                {
                                    if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first )
                                    // std::cout << "Reco pid " << v_prob.at(ivec+1).second <<std::endl;
                                    recopid.push_back( pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ) );
                                }
                            }
                            else{
                                // std::cout << v_prob.at(0).first << " " << v_prob.at(0).second << std::endl;
                                recopid.push_back( pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ) );
                            }
                        } // closes the if statement
                    } // closes the conditional statement of trueparticlename == MC true pdg
                } // closes the vertical bining loop of the pid matrix
            }//close if track_length > tpc_min_length
        } // closes the MC truth loop

        _nFSP.push_back(nFSP);
        // std::cout << "Fill CAF TTree" << std::endl;
        this->FillTTree();

    } // closes the event loop
} // closes the main loop function

bool CAF::isInTPC(TVector3 &point)
{
    //TPC Volume radius 2600 mm
    //TPC full length 5 m
    bool isInside = true;
    float r_point = std::sqrt(point.Y()*point.Y() + point.Z()*point.Z());

    //in the Barrel
    if( r_point > 260 ) isInside = false;
    //in the Endcap
    if(r_point < 260 && std::abs(point.X()) > 250) isInside = false;

    return isInside;
}
