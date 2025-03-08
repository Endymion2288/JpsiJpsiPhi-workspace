#define secCut_cxx
#include "/home/storage0/users/xingcheng/storage2/CMS-Analysis/JpsJpsPhi-workspace/secCut/secCut.h"
#include "/home/storage0/users/xingcheng/storage2/CMS-Analysis/JpsJpsPhi-workspace/includes/ParticleCand.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <memory>
#include <cstdio>

// Include the header file for the roofit.
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooChebychev.h"

#define CUT_DR
#define CUT_PRI_VTXPROB

//#define SHOW_DEBUG

void secCut::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();\
    
    printf("Entries: %lld\n", nentries);

    const unsigned int nBins = 20;
    const unsigned int nBin_cut = 20;
    const unsigned int nCandsAllowed = 20;

    // Use Roofit to draw the plot with proper error bars.
    // Define mass histograms for Jpsi, Phi and Pri passing the cut. Using Roofit.
    RooRealVar Jpsi_1_mass_var("Jpsi_1_mass_cut", "Jpsi_1_mass_cut", 3.0, 3.2);
    RooRealVar Jpsi_2_mass_var("Jpsi_2_mass_cut", "Jpsi_2_mass_cut", 3.0, 3.2);
    RooRealVar Phi_mass_var("Phi_mass_cut","Phi_mass_cut", 0.99, 1.07);
    RooRealVar Pri_mass_var("Pri_mass_cut","Pri_mass_cut", 0.0, 100.0);

    RooRealVar Jpsi_1_ctau_var("Jpsi_1_ctau_cut", "Jpsi_1_ctau_cut", 0.0, 0.1);
    RooRealVar Jpsi_2_ctau_var("Jpsi_2_ctau_cut", "Jpsi_2_ctau_cut", 0.0, 0.1);
    RooRealVar Phi_ctau_var("Phi_ctau_cut","Phi_ctau_cut", 0.0, 0.1);
    RooRealVar Pri_ctau_var("Pri_ctau_cut","Pri_ctau_cut", 0.0, 0.1);

    // Define dataset for Jpsi, Phi and Pri passing the cut. Using Roofit.
    RooDataSet Jpsi_1_mass_set_multi("Jpsi_1_mass_set_multi", "Jpsi_1_mass_set_multi", RooArgList(Jpsi_1_mass_var));
    RooDataSet Jpsi_2_mass_set_multi("Jpsi_2_mass_set_multi", "Jpsi_2_mass_set_multi", RooArgList(Jpsi_2_mass_var));
    RooDataSet Phi_mass_set_multi("Phi_mass_set_multi", "Phi_mass_set_multi", RooArgList(Phi_mass_var));
    RooDataSet Pri_mass_set_multi("Pri_mass_set_multi", "Pri_mass_set_multi", RooArgList(Pri_mass_var));

    // Define dataset for Jpsi, Phi and Pri passing the cut. Using Roofit. Removed the multiple candidates.
    RooDataSet Jpsi_1_mass_set("Jpsi_1_mass_set", "Jpsi_1_mass_set", RooArgList(Jpsi_1_mass_var));
    RooDataSet Jpsi_2_mass_set("Jpsi_2_mass_set", "Jpsi_2_mass_set", RooArgList(Jpsi_2_mass_var));
    RooDataSet Phi_mass_set("Phi_mass_set", "Phi_mass_set", RooArgList(Phi_mass_var));
    RooDataSet Pri_mass_set("Pri_mass_set", "Pri_mass_set", RooArgList(Pri_mass_var));

    RooDataSet Jpsi_1_ctau_set("Jpsi_1_ctau_set", "Jpsi_1_ctau_set", RooArgList(Jpsi_1_ctau_var));
    RooDataSet Jpsi_2_ctau_set("Jpsi_2_ctau_set", "Jpsi_2_ctau_set", RooArgList(Jpsi_2_ctau_var));
    RooDataSet Phi_ctau_set("Phi_ctau_set", "Phi_ctau_set", RooArgList(Phi_ctau_var));
    RooDataSet Pri_ctau_set("Pri_ctau_set", "Pri_ctau_set", RooArgList(Pri_ctau_var));

    // --- Register your cut parameters here  ---

    double Jpsi1_Jpsi2_DR = 0.0;
    double Jpsi1_Phi_DR   = 0.0;
    double Jpsi2_Phi_DR   = 0.0;

    double Jpsi1_Jpsi2_DR_max = 10.0;
    double Jpsi1_Phi_DR_max    = 10.0;
    double Jpsi2_Phi_DR_max    = 10.0;

    double Jpsi1_Jpsi2_DR_min = 0.0;
    double Jpsi1_Phi_DR_min = 0.0;
    double Jpsi2_Phi_DR_min  = 0.0;

    double Ups_pT_min = 6.0;
    double Ups_mu_pT_min = 4.0;

    bool Ups_mu_require_medium = true;
    bool Ups_mu_require_loose  = false;
    bool Ups_mu_require_tight  = false;


    // --- End of cut parameters registration ---

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // Marker set to show the progress.
      if(jentry % 500 == 0){
         printf(">>> Processing entry %ld <<<\n", jentry);
      }

      #ifdef TRY_E4

      if(jentry == 10000){
            break;
      }

      #endif

      // Loop over all candidates.

      std::vector<std::shared_ptr<ParticleCand> > CandList;
      ParticleCand tempCand;
      ParticleCand::PartIdxList_t tempList;
      double temp_massChi2;
      double temp_pt_abs;
      double temp_VtxProb;
      double Pri_VtxProb_cut = 0.75;

      for (unsigned int iCand = 0; iCand < Jpsi_1_mass->size(); iCand++){
            bool passCut = true;

            // Calculate the DRs
            double Jpsi1_Jpsi2_DR = sqrt((Jpsi_1_eta->at(iCand) - Jpsi_2_eta->at(iCand)) * (Jpsi_1_eta->at(iCand) - Jpsi_2_eta->at(iCand))
                                    + (Jpsi_1_phi->at(iCand) - Jpsi_2_phi->at(iCand)) * (Jpsi_1_phi->at(iCand) - Jpsi_2_phi->at(iCand)));
            double Jpsi1_Phi_DR = sqrt((Jpsi_1_eta->at(iCand) - Phi_eta->at(iCand)) * (Jpsi_1_eta->at(iCand) - Phi_eta->at(iCand)) 
                                    + (Jpsi_1_phi->at(iCand) - Phi_phi->at(iCand)) * (Jpsi_1_phi->at(iCand) - Phi_phi->at(iCand)));
            double Jpsi_2_Phi_DR  = sqrt(( Jpsi_2_eta->at(iCand) - Phi_eta->at(iCand)) * ( Jpsi_2_eta->at(iCand) - Phi_eta->at(iCand))
                                    +( Jpsi_2_phi->at(iCand) - Phi_phi->at(iCand)) * ( Jpsi_2_phi->at(iCand) - Phi_phi->at(iCand)));

            #ifdef SHOW_DEBUG
            std::cout << ">>> Calculated DRs <<<" << std::endl;
            #endif

            // Apply DR cuts

            #ifdef CUT_DR

            if (Jpsi1_Jpsi2_DR < Jpsi1_Jpsi2_DR_min || Jpsi1_Jpsi2_DR > Jpsi1_Jpsi2_DR_max){
               passCut = false;
            }

            if (Jpsi1_Phi_DR < Jpsi1_Phi_DR_min || Jpsi1_Phi_DR > Jpsi1_Phi_DR_max){
               passCut = false;
            }

            if (Jpsi2_Phi_DR < Jpsi2_Phi_DR_min || Jpsi2_Phi_DR > Jpsi2_Phi_DR_max){
               passCut = false;
            }

            #endif

            // Apply primary vertex probability cut
            #ifdef CUT_PRI_VTXPROB
            if (Pri_VtxProb->at(iCand) < Pri_VtxProb_cut){
               passCut = false;
            }
            #endif

            // // Apply Upsilon cuts
            // if (Ups_pt->at(iCand) < Ups_pT_min){
            //    passCut = false;
            // }

            // // Apply muon cuts
            // if (Ups_mu_1_pt->at(iCand) < Ups_mu_pT_min || Ups_mu_2_pt->at(iCand) < Ups_mu_pT_min){
            //    passCut = false;
            // }

            // if (Ups_mu_require_medium && (!Ups_mu_1_isPatMediumMuon->at(iCand) || !Ups_mu_2_isPatMediumMuon->at(iCand))){
            //    passCut = false;
            // }

            // if (Ups_mu_require_loose && (!Ups_mu_1_isPatLooseMuon->at(iCand) || !Ups_mu_2_isPatLooseMuon->at(iCand))){
            //    passCut = false;
            // }

            // if (Ups_mu_require_tight && (!Ups_mu_1_isPatTightMuon->at(iCand) || !Ups_mu_2_isPatTightMuon->at(iCand))){
            //    passCut = false;
            // }

            if (!passCut){
               continue;
            }

            // Fill the temporary candidate.
            tempList.clear();
            tempCand.Clear();
            tempCand.SetId(iCand);

            tempList.push_back(Jpsi_1_mu_1_Idx->at(iCand));
            tempList.push_back(Jpsi_1_mu_2_Idx->at(iCand));
            tempList.push_back(Jpsi_2_mu_1_Idx->at(iCand));
            tempList.push_back(Jpsi_2_mu_2_Idx->at(iCand));
            tempCand.AddParticle(ParticleCand::PartType::Muon, tempList);
            tempList.clear();
            tempList.push_back(Phi_K_1_Idx->at(iCand));
            tempList.push_back(Phi_K_2_Idx->at(iCand));
            tempCand.AddParticle(ParticleCand::PartType::Track, tempList);
            tempList.clear();

            temp_massChi2 = (Jpsi_1_massDiff->at(iCand) / Jpsi_1_massErr->at(iCand)) * (Jpsi_1_massDiff->at(iCand) / Jpsi_1_massErr->at(iCand)) +
                           (Jpsi_2_massDiff->at(iCand) /  Jpsi_2_massErr->at(iCand)) *  (Jpsi_2_massDiff->at(iCand) /  Jpsi_2_massErr->at(iCand)) +
                           (Phi_massDiff->at(iCand) /  Phi_massErr->at(iCand)) *  (Phi_massDiff->at(iCand) /  Phi_massErr->at(iCand));
            
            temp_pt_abs = sqrt(Jpsi_1_pt->at(iCand) * Jpsi_1_pt->at(iCand) + Jpsi_2_pt->at(iCand) * Jpsi_2_pt->at(iCand) + Phi_pt->at(iCand) * Phi_pt->at(iCand));

            temp_VtxProb = Jpsi_1_VtxProb->at(iCand) * Jpsi_2_VtxProb->at(iCand) * Phi_VtxProb->at(iCand);

            tempCand.SetScore(temp_VtxProb);
            CandList.push_back(std::make_shared<ParticleCand>(tempCand));
      }

      #ifdef SHOW_DEBUG
      std::cout << ">>> Finished filling the candidate list <<<" << std::endl;
      #endif

      if(CandList.size() == 0){
            continue;
      }
      // Save all filtered candidates to the tree.
      for (auto cand: CandList){
            filtered_Jpsi_1_mass->push_back(Jpsi_1_mass->at(cand->GetId()));
            filtered_Jpsi_1_massErr->push_back(Jpsi_1_massErr->at(cand->GetId()));
            filtered_Jpsi_1_massDiff->push_back(Jpsi_1_massDiff->at(cand->GetId()));
            filtered_Jpsi_1_ctau->push_back(Jpsi_1_ctau->at(cand->GetId()));
            filtered_Jpsi_1_ctauErr->push_back(Jpsi_1_ctauErr->at(cand->GetId()));
            filtered_Jpsi_1_Chi2->push_back(Jpsi_1_Chi2->at(cand->GetId()));
            filtered_Jpsi_1_ndof->push_back(Jpsi_1_ndof->at(cand->GetId()));
            filtered_Jpsi_1_VtxProb->push_back(Jpsi_1_VtxProb->at(cand->GetId()));
            filtered_Jpsi_1_px->push_back(Jpsi_1_px->at(cand->GetId()));
            filtered_Jpsi_1_py->push_back(Jpsi_1_py->at(cand->GetId()));
            filtered_Jpsi_1_pz->push_back(Jpsi_1_pz->at(cand->GetId()));
            filtered_Jpsi_1_phi->push_back(Jpsi_1_phi->at(cand->GetId()));
            filtered_Jpsi_1_eta->push_back(Jpsi_1_eta->at(cand->GetId()));
            filtered_Jpsi_1_pt->push_back(Jpsi_1_pt->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_Idx->push_back(Jpsi_1_mu_1_Idx->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_Idx->push_back(Jpsi_1_mu_2_Idx->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Jpsi 1 to the dataset <<<" << std::endl;
            #endif

            filtered_Phi_mass->push_back(Phi_mass->at(cand->GetId()));
            filtered_Phi_massErr->push_back(Phi_massErr->at(cand->GetId()));
            filtered_Phi_massDiff->push_back(Phi_massDiff->at(cand->GetId()));
            filtered_Phi_ctau->push_back(Phi_ctau->at(cand->GetId()));
            filtered_Phi_ctauErr->push_back(Phi_ctauErr->at(cand->GetId()));
            filtered_Phi_Chi2->push_back(Phi_Chi2->at(cand->GetId()));
            filtered_Phi_ndof->push_back(Phi_ndof->at(cand->GetId()));
            filtered_Phi_VtxProb->push_back(Phi_VtxProb->at(cand->GetId()));
            filtered_Phi_px->push_back(Phi_px->at(cand->GetId()));
            filtered_Phi_py->push_back(Phi_py->at(cand->GetId()));
            filtered_Phi_pz->push_back(Phi_pz->at(cand->GetId()));
            filtered_Phi_phi->push_back(Phi_phi->at(cand->GetId()));
            filtered_Phi_eta->push_back(Phi_eta->at(cand->GetId()));
            filtered_Phi_pt->push_back(Phi_pt->at(cand->GetId()));
            filtered_Phi_K_1_Idx->push_back(Phi_K_1_Idx->at(cand->GetId()));
            filtered_Phi_K_2_Idx->push_back(Phi_K_2_Idx->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Phi to the dataset <<<" << std::endl;
            #endif

            filtered_Pri_mass->push_back(Pri_mass->at(cand->GetId()));
            filtered_Pri_massErr->push_back(Pri_massErr->at(cand->GetId()));
            filtered_Pri_ctau->push_back(Pri_ctau->at(cand->GetId()));
            filtered_Pri_ctauErr->push_back(Pri_ctauErr->at(cand->GetId()));
            filtered_Pri_Chi2->push_back(Pri_Chi2->at(cand->GetId()));
            filtered_Pri_ndof->push_back(Pri_ndof->at(cand->GetId()));
            filtered_Pri_VtxProb->push_back(Pri_VtxProb->at(cand->GetId()));
            filtered_Pri_px->push_back(Pri_px->at(cand->GetId()));
            filtered_Pri_py->push_back(Pri_py->at(cand->GetId()));
            filtered_Pri_pz->push_back(Pri_pz->at(cand->GetId()));
            filtered_Pri_phi->push_back(Pri_phi->at(cand->GetId()));
            filtered_Pri_eta->push_back(Pri_eta->at(cand->GetId()));
            filtered_Pri_pt->push_back(Pri_pt->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Pri to the dataset <<<" << std::endl;
            #endif

            filtered_Jpsi_2_mass->push_back(Jpsi_2_mass->at(cand->GetId()));
            filtered_Jpsi_2_massErr->push_back(Jpsi_2_massErr->at(cand->GetId()));
            filtered_Jpsi_2_massDiff->push_back(Jpsi_2_massDiff->at(cand->GetId()));
            filtered_Jpsi_2_ctau->push_back(Jpsi_2_ctau->at(cand->GetId()));
            filtered_Jpsi_2_ctauErr->push_back(Jpsi_2_ctauErr->at(cand->GetId()));
            filtered_Jpsi_2_Chi2->push_back(Jpsi_2_Chi2->at(cand->GetId()));
            filtered_Jpsi_2_ndof->push_back(Jpsi_2_ndof->at(cand->GetId()));
            filtered_Jpsi_2_VtxProb->push_back(Jpsi_2_VtxProb->at(cand->GetId()));
            filtered_Jpsi_2_px->push_back(Jpsi_2_px->at(cand->GetId()));
            filtered_Jpsi_2_py->push_back(Jpsi_2_py->at(cand->GetId()));
            filtered_Jpsi_2_pz->push_back(Jpsi_2_pz->at(cand->GetId()));
            filtered_Jpsi_2_phi->push_back(Jpsi_2_phi->at(cand->GetId()));
            filtered_Jpsi_2_eta->push_back(Jpsi_2_eta->at(cand->GetId()));
            filtered_Jpsi_2_pt->push_back(Jpsi_2_pt->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_Idx->push_back(Jpsi_2_mu_1_Idx->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_Idx->push_back(Jpsi_2_mu_2_Idx->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Jpsi 2 to the dataset <<<" << std::endl;
            #endif

            // For the muons: copy the corresponding muon information.
            filtered_Jpsi_1_mu_1_px->push_back(Jpsi_1_mu_1_px->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_py->push_back(Jpsi_1_mu_1_py->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_pz->push_back(Jpsi_1_mu_1_pz->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_eta->push_back(Jpsi_1_mu_1_eta->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_phi->push_back(Jpsi_1_mu_1_phi->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_pt->push_back(Jpsi_1_mu_1_pt->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_isPatLooseMuon->push_back(Jpsi_1_mu_1_isPatLooseMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_isPatSoftMuon->push_back(Jpsi_1_mu_1_isPatSoftMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_isPatMediumMuon->push_back(Jpsi_1_mu_1_isPatMediumMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_isPatTightMuon->push_back(Jpsi_1_mu_1_isPatTightMuon->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Jpsi muons to the dataset <<<" << std::endl;
            #endif

            filtered_Jpsi_1_mu_2_px->push_back(Jpsi_1_mu_2_px->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_py->push_back(Jpsi_1_mu_2_py->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_pz->push_back(Jpsi_1_mu_2_pz->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_eta->push_back(Jpsi_1_mu_2_eta->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_phi->push_back(Jpsi_1_mu_2_phi->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_pt->push_back(Jpsi_1_mu_2_pt->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_isPatLooseMuon->push_back(Jpsi_1_mu_2_isPatLooseMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_isPatSoftMuon->push_back(Jpsi_1_mu_2_isPatSoftMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_isPatMediumMuon->push_back(Jpsi_1_mu_2_isPatMediumMuon->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_isPatTightMuon->push_back(Jpsi_1_mu_2_isPatTightMuon->at(cand->GetId()));

            filtered_Jpsi_2_mu_1_px->push_back(Jpsi_2_mu_1_px->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_py->push_back(Jpsi_2_mu_1_py->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_pz->push_back(Jpsi_2_mu_1_pz->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_eta->push_back(Jpsi_2_mu_1_eta->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_phi->push_back(Jpsi_2_mu_1_phi->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_pt->push_back(Jpsi_2_mu_1_pt->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_isPatLooseMuon->push_back(Jpsi_2_mu_1_isPatLooseMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_isPatSoftMuon->push_back(Jpsi_2_mu_1_isPatSoftMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_isPatMediumMuon->push_back(Jpsi_2_mu_1_isPatMediumMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_isPatTightMuon->push_back(Jpsi_2_mu_1_isPatTightMuon->at(cand->GetId()));

            filtered_Jpsi_2_mu_2_px->push_back(Jpsi_2_mu_2_px->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_py->push_back(Jpsi_2_mu_2_py->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_pz->push_back(Jpsi_2_mu_2_pz->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_eta->push_back(Jpsi_2_mu_2_eta->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_phi->push_back(Jpsi_2_mu_2_phi->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_pt->push_back(Jpsi_2_mu_2_pt->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_isPatLooseMuon->push_back(Jpsi_2_mu_2_isPatLooseMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_isPatSoftMuon->push_back(Jpsi_2_mu_2_isPatSoftMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_isPatMediumMuon->push_back(Jpsi_2_mu_2_isPatMediumMuon->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_isPatTightMuon->push_back(Jpsi_2_mu_2_isPatTightMuon->at(cand->GetId()));

            // Kaons from phi
            #ifdef SHOW_DEBUG
            std::cout << ">>> Phi_K_1_px size: " << Phi_K_1_px->size() << " <<<" << std::endl;
            std::cout << ">>> Phi_K_1_py size: " << Phi_K_1_py->size() << " <<<" << std::endl;
            std::cout << ">>> ID: " << cand->GetId() << " <<<" << std::endl;
            #endif
            filtered_Phi_K_1_px->push_back(Phi_K_1_px->at(cand->GetId()));
            filtered_Phi_K_1_py->push_back(Phi_K_1_py->at(cand->GetId()));
            filtered_Phi_K_1_pt->push_back(Phi_K_1_pt->at(cand->GetId()));
            filtered_Phi_K_1_eta->push_back(Phi_K_1_eta->at(cand->GetId()));
            filtered_Phi_K_1_phi->push_back(Phi_K_1_phi->at(cand->GetId()));

            filtered_Phi_K_2_px->push_back(Phi_K_2_px->at(cand->GetId()));
            filtered_Phi_K_2_py->push_back(Phi_K_2_py->at(cand->GetId()));
            filtered_Phi_K_2_pz->push_back(Phi_K_2_pz->at(cand->GetId()));
            filtered_Phi_K_2_pt->push_back(Phi_K_2_pt->at(cand->GetId()));
            filtered_Phi_K_2_eta->push_back(Phi_K_2_eta->at(cand->GetId()));
            filtered_Phi_K_2_phi->push_back(Phi_K_2_phi->at(cand->GetId()));
            #ifdef SHOW_DEBUG
            std::cout << ">>> Finished adding Phi kaons to the dataset <<<" << std::endl;
            #endif    
      }
      filteredTree->Fill();
      // To plot in parallel: multi candidates allowed and non-overlap candidate only.
      // For all candidates, store them all.
      for (auto cand: CandList){
            Jpsi_1_mass_var.setVal(Jpsi_1_mass->at(cand->GetId()));
            Jpsi_2_mass_var.setVal(Jpsi_2_mass->at(cand->GetId()));
            Phi_mass_var.setVal(Phi_mass->at(cand->GetId()));
            Pri_mass_var.setVal(Pri_mass->at(cand->GetId()));

            Jpsi_1_mass_var.setError(Jpsi_1_massErr->at(cand->GetId()));
            Jpsi_2_mass_var.setError(Jpsi_2_massErr->at(cand->GetId()));
            Phi_mass_var.setError(Phi_massErr->at(cand->GetId()));
            Pri_mass_var.setError(Pri_massErr->at(cand->GetId()));

            Jpsi_1_mass_set_multi.add(RooArgSet(Jpsi_1_mass_var));
            Jpsi_2_mass_set_multi.add(RooArgSet(Jpsi_2_mass_var));
            Phi_mass_set_multi.add(RooArgSet(Phi_mass_var));
            Pri_mass_set_multi.add(RooArgSet(Pri_mass_var));
      }

      // For non-overlap candidates, store them.
      std::vector<std::shared_ptr<ParticleCand> > CandList_nonOverlap;
      std::sort(CandList.begin(), CandList.end(),
            [](std::shared_ptr<ParticleCand> cand1, std::shared_ptr<ParticleCand> cand2){
               return cand1->GetScore() > cand2->GetScore();
            }
      ); 
      for (auto& cand: CandList){
            if(CandList_nonOverlap.size() == 0){
               CandList_nonOverlap.push_back(cand);
            }
            else{
               bool isOverlap = false;
               for (auto& cand_nonOverlap: CandList_nonOverlap){
                  if(cand->Overlap(*cand_nonOverlap)){
                        isOverlap = true;
                        break;
                  }
               }
               if(!isOverlap){
                  CandList_nonOverlap.push_back(cand);
               }
            }
      }
      for( auto cand: CandList_nonOverlap){
         Jpsi_1_mass_var.setVal(Jpsi_1_mass->at(cand->GetId()));
         Jpsi_2_mass_var.setVal(Jpsi_2_mass->at(cand->GetId()));
         Phi_mass_var.setVal(Phi_mass->at(cand->GetId()));
         Pri_mass_var.setVal(Pri_mass->at(cand->GetId()));

         Jpsi_1_mass_var.setError(Jpsi_1_massErr->at(cand->GetId()));
         Jpsi_2_mass_var.setError(Jpsi_2_massErr->at(cand->GetId()));
         Phi_mass_var.setError(Phi_massErr->at(cand->GetId()));
         Pri_mass_var.setError(Pri_massErr->at(cand->GetId()));

            Jpsi_1_mass_set.add(RooArgSet(Jpsi_1_mass_var));
            Jpsi_2_mass_set.add(RooArgSet(Jpsi_2_mass_var));
            Phi_mass_set.add(RooArgSet(Phi_mass_var));
            Pri_mass_set.add(RooArgSet(Pri_mass_var));
         
         Jpsi_1_ctau_var.setVal(Jpsi_1_ctau->at(cand->GetId()));
         Jpsi_2_ctau_var.setVal(Jpsi_2_ctau->at(cand->GetId()));
         Phi_ctau_var.setVal(Phi_ctau->at(cand->GetId()));
         Pri_ctau_var.setVal(Pri_ctau->at(cand->GetId()));

         Jpsi_1_ctau_var.setError(Jpsi_1_ctauErr->at(cand->GetId()));
         Jpsi_2_ctau_var.setError(Jpsi_2_ctauErr->at(cand->GetId()));
         Phi_ctau_var.setError(Phi_ctauErr->at(cand->GetId()));
         Pri_ctau_var.setError(Pri_ctauErr->at(cand->GetId()));

         Jpsi_1_ctau_set.add(RooArgSet(Jpsi_1_ctau_var));
         Jpsi_2_ctau_set.add(RooArgSet(Jpsi_2_ctau_var));
         Phi_ctau_set.add(RooArgSet(Phi_ctau_var));
         Pri_ctau_set.add(RooArgSet(Pri_ctau_var));
      }
      ClearBranches();
      CandList.clear();
   }

   // Draw the histograms.
   TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
   c1->Divide(2,2);
   RooPlot* Jpsi_1_mass_frame = Jpsi_1_mass_var.frame(nBins);
   RooPlot* Jpsi_2_mass_frame = Jpsi_2_mass_var.frame(nBins);
   RooPlot* Phi_mass_frame = Phi_mass_var.frame(nBins);
   RooPlot* Pri_mass_frame = Pri_mass_var.frame(nBins);
   Jpsi_1_mass_set_multi.plotOn(Jpsi_1_mass_frame);
   Jpsi_2_mass_set_multi.plotOn(Jpsi_2_mass_frame);
   Phi_mass_set_multi.plotOn(Phi_mass_frame);
   Pri_mass_set_multi.plotOn(Pri_mass_frame);
   c1->cd(1);
   Jpsi_1_mass_frame->Draw();
   c1->cd(2);
   Jpsi_2_mass_frame->Draw();
   c1->cd(3);
   Phi_mass_frame->Draw();
   c1->cd(4);
   Pri_mass_frame->Draw();
   c1->SaveAs("secCut_multi_by_VtxProb.pdf");
   c1->SaveAs("secCut_multi_by_VtxProb.png");

   // Draw the histograms for the non-overlap candidates.
   TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
   c2->Divide(2,2);
   RooPlot* Jpsi_1_mass_frame_nonOverlap = Jpsi_1_mass_var.frame(nBins);
   RooPlot* Jpsi_2_mass_frame_nonOverlap = Jpsi_2_mass_var.frame(nBins);
   RooPlot* Phi_mass_frame_nonOverlap = Phi_mass_var.frame(nBins);
   RooPlot* Pri_mass_frame_nonOverlap = Pri_mass_var.frame(nBins);
   Jpsi_1_mass_set.plotOn(Jpsi_1_mass_frame_nonOverlap);
   Jpsi_2_mass_set.plotOn(Jpsi_2_mass_frame_nonOverlap);
   Phi_mass_set.plotOn(Phi_mass_frame_nonOverlap);
   Pri_mass_set.plotOn(Pri_mass_frame_nonOverlap);
   c2->cd(1);
   Jpsi_1_mass_frame_nonOverlap->Draw();
   c2->cd(2);
   Jpsi_2_mass_frame_nonOverlap->Draw();
   c2->cd(3);
   Phi_mass_frame_nonOverlap->Draw();
   c2->cd(4);
   Pri_mass_frame_nonOverlap->Draw();
   c2->SaveAs("secCut_nonOverlap_by_VtxProb.pdf");
   c2->SaveAs("secCut_nonOverlap_by_VtxProb.png");

   // Draw the histograms for the ctau
   TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
   c3->Divide(2,2);
   RooPlot* Jpsi_1_ctau_frame = Jpsi_1_ctau_var.frame(nBins);
   RooPlot* Jpsi_2_ctau_frame = Jpsi_2_ctau_var.frame(nBins);
   RooPlot* Phi_ctau_frame = Phi_ctau_var.frame(nBins);
   RooPlot* Pri_ctau_frame = Pri_ctau_var.frame(nBins);
   Jpsi_1_ctau_set.plotOn(Jpsi_1_ctau_frame);
   Jpsi_2_ctau_set.plotOn(Jpsi_2_ctau_frame);
   Phi_ctau_set.plotOn(Phi_ctau_frame);
   Pri_ctau_set.plotOn(Pri_ctau_frame);
   c3->cd(1);
   gPad->SetLogy();
   Jpsi_1_ctau_frame->Draw();
   c3->cd(2);
   gPad->SetLogy();
   Jpsi_2_ctau_frame->Draw();
   c3->cd(3);
   gPad->SetLogy();
   Phi_ctau_frame->Draw();
   c3->cd(4);
   gPad->SetLogy();
   Pri_ctau_frame->Draw();
   c3->SaveAs("secCut_ctau.pdf");
   c3->SaveAs("secCut_ctau.png");

   // Save the output tree.
   TFile *outputFile = new TFile("filtered_data_secCut.root", "RECREATE");
   filteredTree->Write();
   outputFile->Close();
}
