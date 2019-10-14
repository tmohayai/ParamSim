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
//// 6) bin/cafanatree_module --infile ${name of the anatree file that is output from anatree module, not to be confused with edepsim file} --outfile ${a name of your choosing for the output file}
//////////////////////////////////////////////////////////////////////////

bool isInTPC(TVector3 &point);

// main loop function
void loop(CAF *caf, TTree *tree)
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
    std::vector<float>   *MCVertX = 0;
    std::vector<float>   *MCVertY = 0;
    std::vector<float>   *MCVertZ = 0;
    std::vector<float>   *MCNuPx = 0;
    std::vector<float>   *MCNuPy = 0;
    std::vector<float>   *MCNuPz = 0;
    std::vector<int>     *NType = 0;
    std::vector<int>     *CCNC = 0;
    std::vector<int>     *Mode = 0;
    std::vector<int>     *Gint=0;
    std::vector<int>     *TgtPDG=0;
    std::vector<int>     *GT_T=0;
    std::vector<int>     *InterT=0;
    std::vector<float>   *Weight=0;

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

    tree->SetBranchAddress("Event", &Event);
    tree->SetBranchAddress("SubRun", &SubRun);
    tree->SetBranchAddress("Run", &Run);

    //Generator info
    tree->SetBranchAddress("NType", &NType);
    tree->SetBranchAddress("CCNC", &CCNC);
    tree->SetBranchAddress("MC_Q2", &MC_Q2);
    tree->SetBranchAddress("MC_W", &MC_W);
    tree->SetBranchAddress("MC_Y", &MC_Y);
    tree->SetBranchAddress("MC_X", &MC_X);
    tree->SetBranchAddress("MC_Theta", &MC_Theta);
    tree->SetBranchAddress("MC_T", &MC_T);
    tree->SetBranchAddress("Mode", &Mode);
    tree->SetBranchAddress("Gint", &Gint);
    tree->SetBranchAddress("TgtPDG", &TgtPDG);
    tree->SetBranchAddress("GT_T", &GT_T);
    tree->SetBranchAddress("MCVertX", &MCVertX);
    tree->SetBranchAddress("MCVertY", &MCVertY);
    tree->SetBranchAddress("MCVertZ", &MCVertZ);
    tree->SetBranchAddress("MCNuPx", &MCNuPx);
    tree->SetBranchAddress("MCNuPy", &MCNuPy);
    tree->SetBranchAddress("MCNuPz", &MCNuPz);
    tree->SetBranchAddress("InterT", &InterT);
    tree->SetBranchAddress("Weight", &Weight);

    //MC info
    tree->SetBranchAddress("PDG", &PDG);
    tree->SetBranchAddress("MCPTrkID", &MCPTrkID);
    tree->SetBranchAddress("MCPStartX", &MCPStartX);
    tree->SetBranchAddress("MCPStartY", &MCPStartY);
    tree->SetBranchAddress("MCPStartZ", &MCPStartZ);
    tree->SetBranchAddress("MCPEndX", &MCPEndX);
    tree->SetBranchAddress("MCPEndY", &MCPEndY);
    tree->SetBranchAddress("MCPEndZ", &MCPEndZ);
    tree->SetBranchAddress("Mother", &Mother);
    tree->SetBranchAddress("PDGMother", &PDGMother);
    tree->SetBranchAddress("MCPStartPX", &MCPStartPX);
    tree->SetBranchAddress("MCPStartPY", &MCPStartPY);
    tree->SetBranchAddress("MCPStartPZ", &MCPStartPZ);
    tree->SetBranchAddress("MCPProc", &MCPProc);
    tree->SetBranchAddress("MCPEndProc", &MCPEndProc);
    tree->SetBranchAddress("TrajMCPX", &TrajMCPX);
    tree->SetBranchAddress("TrajMCPY", &TrajMCPY);
    tree->SetBranchAddress("TrajMCPZ", &TrajMCPZ);
    tree->SetBranchAddress("TrajMCPTrajIndex", &TrajMCPTrajIndex);

    //gamma, neutron, pi0, k0L, k0S, k0, delta0
    std::vector<int> pdg_neutral = {22, 2112, 111, 130, 310, 311, 2114};

    //-------------------------------------------------------------------

    // Main event loop
    for( int entry = 0; entry < tree->GetEntries(); entry++ )
    {
        tree->GetEntry(entry);

        caf->ClearVectors();

        //Filling MCTruth values
        // std::cout << "Event " << Event << " Run " << Run << std::endl;
        caf->Event = Event;
        caf->Run = Run;
        caf->SubRun = SubRun;

        for(size_t i = 0; i < NType->size(); i++)
        {
            caf->ccnc.push_back(CCNC->at(i));
            caf->ntype.push_back(NType->at(i));
            caf->q2.push_back(MC_Q2->at(i));
            caf->w.push_back(MC_W->at(i));
            caf->y.push_back(MC_Y->at(i));
            caf->x.push_back(MC_X->at(i));
            caf->theta.push_back(MC_Theta->at(i));
            caf->mode.push_back(Mode->at(i));
            caf->t.push_back(MC_T->at(i));
            caf->intert.push_back(InterT->at(i));
            caf->vertx.push_back(MCVertX->at(i));
            caf->verty.push_back(MCVertY->at(i));
            caf->vertz.push_back(MCVertZ->at(i));
            caf->mcnupx.push_back(MCNuPx->at(i));
            caf->mcnupy.push_back(MCNuPy->at(i));
            caf->mcnupz.push_back(MCNuPz->at(i));
        }

        for(size_t i = 0; i < Gint->size(); i++)
        {
            caf->gint.push_back(Gint->at(i));
            caf->tgtpdg.push_back(TgtPDG->at(i));
            caf->gt_t.push_back(GT_T->at(i));
            caf->weight.push_back(GT_T->at(i));
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

            //Need to check what kind of particle?
            // if(mcp_process != "primary") continue;

            nFSP++;

            int pdg = PDG->at(i);
            //need to ignore neutrals for this - put the value to 0
            auto result = std::find(pdg_neutral.begin(), pdg_neutral.end(), pdg);
            bool isNeutral = (result == pdg_neutral.end()) ? false : true;

            if( isNeutral )
            {
                caf->trkLen.push_back(0.);
                caf->trkLenPerp.push_back(0.);
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
                        break;
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

            caf->trkLen.push_back(tracklen);
            caf->trkLenPerp.push_back(tracklen_perp);

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

            //Get ending process
            std::string mcp_endprocess = MCPEndProc->at(i);
            //Save MC process
            /*caf->MCProc.push_back(mcp_process);
            caf->MCEndProc.push_back(mcp_endprocess);
	        */
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
                    caf->recopid.push_back(2112);
                    float ereco = rando->Gaus( std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass), sigmaNeutronECAL_first );
                    caf->erecon.push_back(ereco > 0 ? ereco : 0.);
                    // std::cout << "true part n true energy " << std::sqrt(ptrue*ptrue + neutron_mass*neutron_mass) << " ereco " << caf->erecon[i] << std::endl;
                }
            }

            //for pi0s
            if(pdg == 111)
            {
                //TODO smear the pi0 energy (and decay vertex?) according to previous pi0 reco studies
                float ereco = rando->Gaus( std::sqrt(ptrue*ptrue + pi0_mass*pi0_mass), ECAL_pi0_resolution*std::sqrt(ptrue*ptrue + pi0_mass*pi0_mass));
                caf->erecon.push_back(ereco);
                caf->recopid.push_back(111);
            }

            //for gammas
            if(pdg == 22)
            {
                //TODO check if they are not from a pi0 or decayed in the TPC
                if( PDGMother->at(i) != 111 )
                {
                    TVector3 point(MCPEndX->at(i), MCPEndY->at(i), MCPEndZ->at(i));
                    //Endpoint is not in the TPC
                    if(not isInTPC(point)){
                        //if they hit the ECAL and smear their energy
                        float ECAL_resolution = fRes->Eval(ptrue)*ptrue;
                        float ereco = rando->Gaus(ptrue, ECAL_resolution);
                        caf->erecon.push_back(ereco);
                        caf->recopid.push_back(22);
                    }
                }
            }

            //Visible in the TPC
            if( caf->trkLen.at(i) > gastpc_len )
            {
                //Use range instead of Gluckstern for stopping tracks
                //TODO is that correct? What if it is a scatter in the TPC? Need to check if daughter is same particle
                float preco = 0;
                TVector3 point(MCPEndX->at(i), MCPEndY->at(i), MCPEndZ->at(i));
                // moved all the truth-level MC vectors to here to have same-sized
                // vectors between truth and reco
		        caf->truepdg.push_back(pdg);
            	caf->truepx.push_back(MCPStartPX->at(i));
            	caf->truepy.push_back(MCPStartPY->at(i));
            	caf->truepz.push_back(MCPStartPZ->at(i));
            	caf->MCPStartX.push_back(MCPStartX->at(i));
            	caf->MCPStartY.push_back(MCPStartY->at(i));
            	caf->MCPStartZ.push_back(MCPStartZ->at(i));
            	caf->MCPEndX.push_back(MCPEndX->at(i));
            	caf->MCPEndY.push_back(MCPEndY->at(i));
            	caf->MCPEndZ.push_back(MCPEndZ->at(i));
            	caf->mother.push_back(Mother->at(i));
            	caf->pdgmother.push_back(PDGMother->at(i));
            	// save the true momentum
                caf->truep.push_back(ptrue);
                // save the true angle
                caf->angle.push_back(angle);
		        caf->MCProc.push_back(mcp_process);
            	caf->MCEndProc.push_back(mcp_endprocess);
		        if(isInTPC(point))
                {
                    preco = rando->Gaus( ptrue, sigmaP_short );
                    float angle_reco = rando->Gaus(angle, sigma_angle_short);
                    caf->preco.push_back(preco);
                    caf->anglereco.push_back(angle_reco);
                }
                else{
                    // calculate number of trackpoints
                    float nHits = round (caf->trkLen.at(i) / gastpc_padPitch);
                    // measurement term in Gluckstern formula
                    float fracSig_meas = sqrt(720./(nHits+4)) * ((0.01*gastpc_padPitch*ptrue) / (0.3 * gastpc_B * 0.0001 *caf->trkLenPerp.at(i)*caf->trkLenPerp.at(i)));
                    // multiple Coulomb scattering term in Gluckstern formula
                    float fracSig_MCS = (0.052*sqrt(1.43)) / (gastpc_B * sqrt(gastpc_X0*caf->trkLenPerp.at(i)*0.0001));
                    // momentum resoltion from the two terms above
                    float sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
                    // now Gaussian smear the true momentum using the momentum resolution
                    preco = rando->Gaus( ptrue, sigmaP );

                    // measurement term in the Gluckstern formula for calculating the
                    // angular resolution
                    float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / caf->trkLen.at(i)*caf->trkLen.at(i)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
                    // scattering term in Gluckstern formula
                    float sigma_angle_2 = (0.015*0.015 / (3. * ptrue * ptrue)) * (caf->trkLen.at(i)/gastpc_X0);
                    // angular resolution from the two terms above
                    float sigma_angle = sqrt(sigma_angle_1 + sigma_angle_2);
                    // now Gaussian smear the true angle using the angular resolution
                    float angle_reco = rando->Gaus(angle, sigma_angle);

                    // save reconstructed momentum and angle to cafanatree
                    caf->preco.push_back(preco);
                    caf->anglereco.push_back(angle_reco);
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
                                    caf->prob_arr.push_back(prob);
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
                                    caf->recopid.push_back( pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ) );
                                }
                            }
                            else{
                                // std::cout << v_prob.at(0).first << " " << v_prob.at(0).second << std::endl;
                                caf->recopid.push_back( pdglist.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ) );
                            }
                        } // closes the if statement
                    } // closes the conditional statement of trueparticlename == MC true pdg
                } // closes the vertical bining loop of the pid matrix
            }//close if track_length > tpc_min_length
        } // closes the MC truth loop

        caf->nFSP.push_back(nFSP);
        // std::cout << "Fill CAF TTree" << std::endl;
        caf->FillTTree();

    } // closes the event loop
} // closes the main loop function

bool isInTPC(TVector3 &point)
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

void ShowHelp()
{
    std::cout << "./cafanatree_module --infile <inputfile> --outfile <outputfile>" << std::endl;
}


int main( int argc, char const *argv[] )
{
    if( argc == 1 || ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) || argc < 5 ){
        ShowHelp();
        return 2;
    }

    if( argv[1] != std::string("--infile") || argv[3] != std::string("--outfile") ) {
        ShowHelp();
        return -2;
    }

    // get command line options
    std::string outfile = "";
    std::string infile = "";
    int p = 0;
    while( p < argc )
    {
        if( argv[p] == std::string("--infile") ){
            infile = argv[p+1];
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

    printf( "Making CAF from edep-sim tree dump: %s\n", infile.c_str() );
    printf( "Output CAF file: %s\n", outfile.c_str() );

    TFile * tf = new TFile( infile.c_str() );
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
