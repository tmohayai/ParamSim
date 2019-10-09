#include "CAF.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <numeric>
#include <functional>

typedef std::pair<float, std::string> P;

////////////////////////////////////////////////////////////////////////
//// Class:       cafanatree
//// File Name:   cafanatree_module.cc
////
//// Authors: Tanaz Mohayai and Eldwan Brianne
//// To run this module:
//// 1) cd Build
//// 2) rm -rf *
//// 3) cmake ${module home direcory}
//// 4) make install
//// 5) cd ${module home direcory}
//// 6) bin/cafanatree_module --edepfile ${name of the anatree file that is output from anatree module, not to be confused with edepsim file} --outfile ${a name of your choosing for the output file}
//////////////////////////////////////////////////////////////////////////

// main loop function
void loop(CAF *caf, TTree *tree)
{
    //double gastpc_len = 5.; // track length cut in cm
    double gastpc_len = 2.; // new track length cut in cm based on Thomas' study of low energy protons
    // dont care about electrons -- check momentum and see if hit ECAL
    double gastpc_B = 0.5; // B field strength in Tesla
    double gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
    double gastpc_X0 = 1300.; // cm = 13m radiation length

    int seed = 7;
    TRandom3 *rando = new TRandom3( seed );

    std::vector<double> v = std::vector<double>();
    for (double pit=0.040; pit<20.0; pit+=0.001)
    {
        v.push_back(pit);
    }

    TRandom *R = new TRandom(time(0));

    //TODO as function of momentum
    float sigmaNeutronECAL_first = 0.11;
    float sigmaNeutronECAL_rescatter = 0.26;

    //------------------------------------------------------------------------

    Int_t           	 Event;
    Int_t           	 SubRun;
    Int_t           	 Run;
    std::vector<int>     *NType = 0;
    std::vector<int>     *CCNC = 0;
    std::vector<int>     *PDG = 0;
    std::vector<float>   *TrackLenF = 0;
    std::vector<float>   *TrackLenB = 0;
    std::vector<int>     *NTPCClustersOnTrack = 0;
    std::vector<float>   *MCNuPx = 0;
    std::vector<float>   *MCNuPy = 0;
    std::vector<float>   *MCNuPz = 0;
    std::vector<int>     *MCPTrkID = 0;
    std::vector<float>   *MCPStartX = 0;
    std::vector<float>   *MCPStartY = 0;
    std::vector<float>   *MCPStartZ = 0;
    std::vector<float>   *MCPStartPX = 0;
    std::vector<float>   *MCPStartPY = 0;
    std::vector<float>   *MCPStartPZ = 0;
    std::vector<std::string>   *MCPProc = 0;
    std::vector<float>   *TrackStartX = 0;
    std::vector<float>   *TrackStartY = 0;
    std::vector<float>   *TrackStartZ = 0;
    std::vector<float>   *TrackStartPX = 0;
    std::vector<float>   *TrackStartPY = 0;
    std::vector<float>   *TrackStartPZ = 0;
    std::vector<float>   *TrackEndX = 0;
    std::vector<float>   *TrackEndY = 0;
    std::vector<float>   *TrackEndZ = 0;
    std::vector<float>   *TrackEndPX = 0;
    std::vector<float>   *TrackEndPY = 0;
    std::vector<float>   *TrackEndPZ = 0;
    std::vector<float>   *VertX = 0;
    std::vector<float>   *VertY = 0;
    std::vector<float>   *VertZ = 0;
    std::vector<int>     *VertN = 0;
    std::vector<int> 	 *TrackIDNumber = 0;

    std::vector<float>   *MCPEndX = 0;
    std::vector<float>   *MCPEndY = 0;
    std::vector<float>   *MCPEndZ = 0;

    std::vector<float>   *TrajMCPX = 0;
    std::vector<float>   *TrajMCPY = 0;
    std::vector<float>   *TrajMCPZ = 0;
    std::vector<int>     *TrajMCPTrajIndex = 0;

    TBranch        *b_Event;
    TBranch        *b_SubRun;
    TBranch        *b_Run;
    TBranch        *b_NType;
    TBranch        *b_CCNC;
    TBranch        *b_PDG;
    TBranch	       *b_MCPTrkID;
    TBranch        *b_TrackLenB;
    TBranch 	   *b_NTPCClustersOnTrack;
    TBranch 	   *b_TrackLenF;
    TBranch        *b_MCNuPx;
    TBranch        *b_MCNuPy;
    TBranch        *b_MCNuPz;
    TBranch        *b_MCPStartX;
    TBranch        *b_MCPStartY;
    TBranch        *b_MCPStartZ;
    TBranch        *b_MCPStartPX;
    TBranch        *b_MCPStartPY;
    TBranch        *b_MCPStartPZ;
    TBranch        *b_MCPProc;
    TBranch        *b_TrackStartX;
    TBranch        *b_TrackStartY;
    TBranch        *b_TrackStartZ;
    TBranch        *b_TrackStartPX;
    TBranch        *b_TrackStartPY;
    TBranch        *b_TrackStartPZ;
    TBranch        *b_TrackEndX;
    TBranch        *b_TrackEndY;
    TBranch        *b_TrackEndZ;
    TBranch        *b_TrackEndPX;
    TBranch        *b_TrackEndPY;
    TBranch        *b_TrackEndPZ;
    TBranch        *b_VertX;
    TBranch        *b_VertY;
    TBranch        *b_VertZ;
    TBranch        *b_VertN;
    TBranch 	   *b_TrackIDNumber;

    TBranch        *b_MCPEndX;
    TBranch        *b_MCPEndY;
    TBranch        *b_MCPEndZ;

    TBranch        *b_TrajMCPX;
    TBranch        *b_TrajMCPY;
    TBranch        *b_TrajMCPZ;
    TBranch        *b_TrajMCPTrajIndex;

    NType = 0;
    CCNC = 0;
    PDG = 0;
    TrackIDNumber = 0;
    MCPTrkID = 0;
    TrackLenF = 0;
    NTPCClustersOnTrack = 0;
    TrackLenB = 0;
    MCPStartX = 0;
    MCPStartY = 0;

    TrajMCPX= 0;
    TrajMCPY = 0;
    TrajMCPZ = 0;

    MCPStartZ = 0;
    MCNuPx = 0;
    MCNuPy = 0;
    MCNuPz = 0;
    MCPStartPX = 0;
    MCPStartPY = 0;
    MCPStartPZ = 0;

    MCPEndX = 0;
    MCPEndY = 0;
    MCPEndZ = 0;

    TrackStartX = 0;
    TrackStartY = 0;
    TrackStartZ = 0;
    TrackStartPX = 0;
    TrackStartPY = 0;
    TrackStartPZ = 0;
    TrackEndX = 0;
    TrackEndY = 0;
    TrackEndZ = 0;
    TrackEndPX = 0;
    TrackEndPY = 0;
    TrackEndPZ = 0;
    VertX = 0;
    VertY = 0;
    VertZ = 0;
    VertN = 0;

    tree->SetBranchAddress("Event", &Event, &b_Event);
    tree->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
    tree->SetBranchAddress("Run", &Run, &b_Run);
    tree->SetBranchAddress("NType", &NType, &b_NType);
    tree->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
    tree->SetBranchAddress("PDG", &PDG, &b_PDG);
    tree->SetBranchAddress("MCPTrkID", &MCPTrkID, &b_MCPTrkID);
    tree->SetBranchAddress("TrackLenF", &TrackLenF, &b_TrackLenF);
    tree->SetBranchAddress("TrackLenB", &TrackLenB, &b_TrackLenB);
    tree->SetBranchAddress("NTPCClustersOnTrack", &NTPCClustersOnTrack, &b_NTPCClustersOnTrack);
    tree->SetBranchAddress("TrackIDNumber", &TrackIDNumber, &b_TrackIDNumber);

    tree->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    tree->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    tree->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);

    tree->SetBranchAddress("MCPEndX", &MCPEndX, &b_MCPEndX);
    tree->SetBranchAddress("MCPEndY", &MCPEndY, &b_MCPEndY);
    tree->SetBranchAddress("MCPEndZ", &MCPEndZ, &b_MCPEndZ);


    tree->SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    tree->SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    tree->SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    tree->SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);

    tree->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    tree->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    tree->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    tree->SetBranchAddress("TrajMCPTrajIndex", &TrajMCPTrajIndex, &b_TrajMCPTrajIndex);

    tree->SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
    tree->SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
    tree->SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);

    tree->SetBranchAddress("TrackStartX", &TrackStartX, &b_TrackStartX);
    tree->SetBranchAddress("TrackStartY", &TrackStartY, &b_TrackStartY);
    tree->SetBranchAddress("TrackStartZ", &TrackStartZ, &b_TrackStartZ);
    tree->SetBranchAddress("TrackStartPX", &TrackStartPX, &b_TrackStartPX);
    tree->SetBranchAddress("TrackStartPY", &TrackStartPY, &b_TrackStartPY);
    tree->SetBranchAddress("TrackStartPZ", &TrackStartPZ, &b_TrackStartPZ);
    tree->SetBranchAddress("TrackEndX", &TrackEndX, &b_TrackEndX);
    tree->SetBranchAddress("TrackEndY", &TrackEndY, &b_TrackEndY);
    tree->SetBranchAddress("TrackEndZ", &TrackEndZ, &b_TrackEndZ);
    tree->SetBranchAddress("TrackEndPX", &TrackEndPX, &b_TrackEndPX);
    tree->SetBranchAddress("TrackEndPY", &TrackEndPY, &b_TrackEndPY);
    tree->SetBranchAddress("TrackEndPZ", &TrackEndPZ, &b_TrackEndPZ);
    tree->SetBranchAddress("VertX", &VertX, &b_VertX);
    tree->SetBranchAddress("VertY", &VertY, &b_VertY);
    tree->SetBranchAddress("VertZ", &VertZ, &b_VertZ);
    tree->SetBranchAddress("VertN", &VertN, &b_VertN);

    //gamma, neutron, pi0, k0L, k0S, k0, delta0
    std::vector<int> pdg_neutral = {22, 2112, 111, 130, 310, 311, 2114};

    // Main event loop
    int N = tree->GetEntries();
    for( int entry = 0; entry < N; entry++ )
    // for( int entry = 0; entry < 1; entry++ )
    {
        tree->GetEntry(entry);
        caf->event = Event;

        int nprimary = 0;
        //--------------------------------------------------------------------------
        // Start of Parameterized Reconstruction
        //--------------------------------------------------------------------------

        caf->Ev_rec[entry] = 0.;
        // save the number of final state particles that are primary
        for(size_t i=0; i < MCPStartPX->size(); ++i )
        {
            std::string mcp_process = MCPProc->at(i);
            // std::cout << mcp_process << std::endl;
            if(mcp_process != "primary") continue;
            nprimary++;

            int pdg = PDG->at(i);
            //need to ignore neutrals for this - put the value to 0
            auto result = std::find(pdg_neutral.begin(), pdg_neutral.end(), pdg);
            bool isNeutral = (result == pdg_neutral.end()) ? false : true;

            if( isNeutral )
            {
                caf->trkLen[i] = 0.;
                caf->trkLenPerp[i] = 0.;
            }

            // calculate the total and the transverse track lengths and restrict the
            // tracklength to be above the gas TPC track length threshold
            double tracklen = 0.;
            double tracklen_perp = 0.;

            for(int itraj = 1; itraj < TrajMCPX->size(); itraj++){
                //check that it is the correct mcp
                if(TrajMCPTrajIndex->at(itraj) == i){
                    // find the length of the track by getting the distance between each hit
                    TVector3 diff(TrajMCPX->at(itraj) - TrajMCPX->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1), TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1));
                    // perp length
                    TVector2 tracklen_perp_vec(TrajMCPZ->at(itraj) - TrajMCPZ->at(itraj-1), TrajMCPY->at(itraj) - TrajMCPY->at(itraj-1));
                    // Summing up
                    tracklen += diff.Mag();
                    tracklen_perp += tracklen_perp_vec.Mod();
                }
            }

            std::cout << "tracklength for mcp " << i << " is: " << tracklen << " cm" << std::endl;

            caf->trkLen[i] = tracklen;
            caf->trkLenPerp[i] = tracklen_perp;
        }

        caf->nFSP = nprimary;

        //---------------------------------------------------------------

        // all Gluckstern calculations happen in the following loop
        for(size_t i=0; i< MCPStartPX->size(); ++i )
        {
            //check if mcp is primary
            std::string mcp_process = MCPProc->at(i);
            // std::cout << mcp_process << std::endl;
            if(mcp_process != "primary") continue;
            //std::cout << "Treating mcp " << i << std::endl;

            TVector3 mcp(MCPStartPX->at(i),MCPStartPY->at(i),MCPStartPZ->at(i));
            TVector3 xhat(1,0,0);
            float ptrue = (mcp).Mag();
            float pz = mcp.Z();
            float pt = (mcp.Cross(xhat)).Mag();
            float px = mcp.X();
            float py = mcp.Y();
            int pdg = PDG->at(i);
            float mctrackid = MCPTrkID->at(i);
            // point resolution
            double sigma_x = 0.1;
            // angle with respect to the incoming neutrino
            float angle  = atan(mcp.X() / mcp.Z());
            // save the true PDG, parametrized PID comes later
            caf->truepdg[i] = pdg;
            // save the true momentum
            caf->truep[i] = ptrue;
            caf->truepx[i] = MCPStartPX->at(i);
            caf->truepy[i] = MCPStartPY->at(i);
            caf->truepz[i] = MCPStartPZ->at(i);
            caf->angle[i] = angle;

            //for neutrons
            if(pdg == 2112)
            {
                //check if it can be detected by the ECAL
                //Assumes 40% efficiency to detect
                double random_number = R->Rndm();
                if(random_number > 0.4)
                {
                    //TODO random is first interaction of rescatter and smear accordingly to Chris's study
                    //Detected in the ECAL
                    caf->recopid[i] = 2112;
                    float ereco = rando->Gaus( std::sqrt(ptrue*ptrue + 0.93957*0.93957), sigmaNeutronECAL_first );
                    caf->erecon[i] = ereco > 0 ? ereco : 0.;
                    std::cout << "true part n true energy " << std::sqrt(ptrue*ptrue + 0.93957*0.93957) << " ereco " << caf->erecon[i] << std::endl;
                }
                else{
                    //not detected drop it
                    caf->erecon[i] = 0.;
                }
            }

            //for pi0s
            if(pdg == 111)
            {
                //TODO smear the pi0 energy (and decay vertex?) according to previous pi0 reco studies
            }

            if( caf->trkLen[i] > gastpc_len )
            {
                //std::cout << "tracklength is:" << trklen_arr[i] << '\n';
                // calculate number of trackpoints
                float nHits = round (caf->trkLen[i] / gastpc_padPitch);
                // measurement term in Gluckstern formula
                float fracSig_meas = sqrt(720./(nHits+4)) * ((0.01*gastpc_padPitch*ptrue) / (0.3 * gastpc_B * 0.0001 *caf->trkLenPerp[i]*caf->trkLenPerp[i]));
                // multiple Coulomb scattering term in Gluckstern formula
                float fracSig_MCS = (0.052*sqrt(1.43)) / (gastpc_B * sqrt(gastpc_X0*caf->trkLenPerp[i]*0.0001));
                // momentum resoltion from the two terms above
                float sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
                // now Gaussian smear the true momentum using the momentum resolution
                float preco = rando->Gaus( ptrue, sigmaP );

                // measurement term in the Gluckstern formula for calculating the
                // angular resolution
                float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / caf->trkLen[i]*caf->trkLen[i]*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
                // scattering term in Gluckstern formula
                float sigma_angle_2 = (0.015*0.015 / (3. * ptrue * ptrue)) * (caf->trkLen[i]/gastpc_X0);
                // angular resolution from the two terms above
                float sigma_angle = sqrt(sigma_angle_1 + sigma_angle_2);
                // now Gaussian smear the true angle using the angular resolution
                float angle_reco = rando->Gaus(angle, sigma_angle);
                // save reconstructed momentum and angle to cafanatree
                caf->preco[i] = preco;
                caf->anglereco[i] = angle_reco;

                //--------------------------------------------------------------------------
                // Start of PID Parametrization
                //--------------------------------------------------------------------------

                float p = caf->preco[i];
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
                for (int q=0; q < 501; ++q)
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
                for (int pidm=1; pidm <= 6; ++pidm)
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

                                    std::cout << "true part " << trueparticlename << " true pid " << pdglist[pidm] << " reco name " << recoparticlename << " reco part list "
                                    << recopnamelist[pidr] <<  " true mom " << ptrue << " reco mom " <<  p << " prob " << pidinterp->GetBinContent(pidm,pidr) << '\n';

                                    //Need to check random number value and prob value then associate the recopdg to the reco prob
                                    v_prob.push_back( std::make_pair(prob, recoparticlename) );
                                }
                            }

                            if(v_prob.size() > 1){
                                //Order the vector of prob
                                std::sort(v_prob.begin(), v_prob.end());
                                //Throw a random number between 0 and 1
                                double random_number = R->Rndm();
                                //Make cumulative sum to get the range
                                std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& x, const P& y){return P(x.first + y.first, y.second);});

                                // std::cout << "rand " << random_number << std::endl;
                                // for(int ivec = 0; ivec < v_prob.size(); ivec++)
                                // {
                                //     std::cout << "Cumulative prob " << v_prob.at(ivec).first << " particle " << v_prob.at(ivec).second << std::endl;
                                // }

                                for(int ivec = 0; ivec < v_prob.size()-1; ivec++)
                                {
                                    if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first )
                                    // std::cout << "Reco pid " << v_prob.at(ivec+1).second <<std::endl;
                                    caf->recopid[i] = pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) );
                                }
                            }
                            else{
                                std::cout << v_prob.at(0).first << " " << v_prob.at(0).second << std::endl;
                                caf->recopid[i] = pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) );
                            }
                        } // closes the if statement
                    } // closes the conditional statement of trueparticlename == MC true pdg
                } // closes the vertical bining loop of the pid matrix
            }//close if track_length > tpc_min_length
        } // closes the MC truth loop
    } // closes the event loop

    std::cout << "Fill CAF TTree" << std::endl;
    caf->FillTTree();
} // closes the main loop function

void ShowHelp()
{
    std::cout << "./cafanatree_module --edepfile <inputfile> --outfile <outputfile>" << std::endl;
}


int main( int argc, char const *argv[] )
{
    if( argc == 1 || ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) || argc < 5 ){
        ShowHelp();
        return 2;
    }

    if( argv[1] != std::string("--edepfile") || argv[3] != std::string("--outfile") ) {
        ShowHelp();
        return -2;
    }

    // get command line options
    std::string outfile = "";
    std::string edepfile = "";
    int p = 0;
    while( p < argc )
    {
        if( argv[p] == std::string("--edepfile") ){
            edepfile = argv[p+1];
            p++;
        }
        else if( argv[p] == std::string("--outfile") ){
            outfile = argv[p+1];
            p++;
        }
        else{
            p++;
        }
    }

    printf( "Making CAF from edep-sim tree dump: %s\n", edepfile.c_str() );
    printf( "Output CAF file: %s\n", outfile.c_str() );

    TFile * tf = new TFile( edepfile.c_str() );
    if(nullptr == tf)
    return -1;
    TTree * tree = (TTree*) tf->Get( "anatree/GArAnaTree" );
    if(nullptr == tree)
    return -1;

    CAF *caf = new CAF(outfile);
    if(not caf->BookTFile()) return -1;

    loop(caf, tree);

    caf->WriteTTree();
    printf( "-30-\n" );

    return 0;
}
