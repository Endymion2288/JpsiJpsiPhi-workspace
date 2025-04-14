#define secCut_cxx
#include "/home/storage2/users/xingcheng/CMSSW_14_0_18/src/JpsiJpsiPhi-workspace/secCut/secCut.h"
#include "/home/storage2/users/xingcheng/CMSSW_14_0_18/src/JpsiJpsiPhi-workspace/includes/ParticleCand.C"
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
#include "RooAddPdf.h"
#include "RooProdPdf.h"

//#define CUT_DR
//#define CUT_PRI_VTXPROB
//#define CUT_PHI_VTXPROB
//#define CUT_PHI_CTAU
#define Stefanos_CUT

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

    RooRealVar Jpsi_1_ctau_var("Jpsi_1_ctau_cut", "Jpsi_1_ctau_cut", -0.05, 0.1);
    RooRealVar Jpsi_2_ctau_var("Jpsi_2_ctau_cut", "Jpsi_2_ctau_cut", -0.05, 0.1);
    RooRealVar Phi_ctau_var("Phi_ctau_cut","Phi_ctau_cut", -0.1, 0.1);
    RooRealVar Pri_ctau_var("Pri_ctau_cut","Pri_ctau_cut", -0.1, 0.1);

   //  RooRealVar Jpsi_1_ctau_error_var("Jpsi_1_ctau_error_cut", "Jpsi_1_ctau_error_cut", 0.0, 0.01);
   //    RooRealVar Jpsi_2_ctau_error_var("Jpsi_2_ctau_error_cut", "Jpsi_2_ctau_error_cut", 0.0, 0.01);
   //    RooRealVar Phi_ctau_error_var("Phi_ctau_error_cut","Phi_ctau_error_cut", 0.0, 0.1);
   //    RooRealVar Pri_ctau_error_var("Pri_ctau_error_cut","Pri_ctau_error_cut", 0.0, 0.1);

   RooRealVar Jpsi_1_Lxy_var("Jpsi_1_Lxy_cut", "Jpsi_1_Lxy_cut", -0.5, 0.5);
   RooRealVar Jpsi_2_Lxy_var("Jpsi_2_Lxy_cut", "Jpsi_2_Lxy_cut", -0.5, 0.5);
   RooRealVar Phi_Lxy_var("Phi_Lxy_cut","Phi_Lxy_cut", -0.5, 0.5);

    RooRealVar Jpsi_1_VtxProb_var("Jpsi_1_VtxProb_cut", "Jpsi_1_VtxProb_cut", 0.0, 1.0);
    RooRealVar Jpsi_2_VtxProb_var("Jpsi_2_VtxProb_cut", "Jpsi_2_VtxProb_cut", 0.0, 1.0);
    RooRealVar Phi_VtxProb_var("Phi_VtxProb_cut","Phi_VtxProb_cut", 0.0, 1.0);
    RooRealVar Pri_VtxProb_var("Pri_VtxProb_cut","Pri_VtxProb_cut", 0.0, 1.0);

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

   //  RooDataSet Jpsi_1_ctau_error_set("Jpsi_1_ctau_error_set", "Jpsi_1_ctau_error_set", RooArgList(Jpsi_1_ctau_error_var));
   //  RooDataSet Jpsi_2_ctau_error_set("Jpsi_2_ctau_error_set", "Jpsi_2_ctau_error_set", RooArgList(Jpsi_2_ctau_error_var));
   //  RooDataSet Phi_ctau_error_set("Phi_ctau_error_set", "Phi_ctau_error_set", RooArgList(Phi_ctau_error_var));
   //  RooDataSet Pri_ctau_error_set("Pri_ctau_error_set", "Pri_ctau_error_set", RooArgList(Pri_ctau_error_var));

   RooDataSet Jpsi_1_Lxy_set("Jpsi_1_Lxy_set", "Jpsi_1_Lxy_set", RooArgList(Jpsi_1_Lxy_var));
   RooDataSet Jpsi_2_Lxy_set("Jpsi_2_Lxy_set", "Jpsi_2_Lxy_set", RooArgList(Jpsi_2_Lxy_var));
   RooDataSet Phi_Lxy_set("Phi_Lxy_set", "Phi_Lxy_set", RooArgList(Phi_Lxy_var));

    RooDataSet Jpsi_1_VtxProb_set("Jpsi_1_VtxProb_set", "Jpsi_1_VtxProb_set", RooArgList(Jpsi_1_VtxProb_var));
    RooDataSet Jpsi_2_VtxProb_set("Jpsi_2_VtxProb_set", "Jpsi_2_VtxProb_set", RooArgList(Jpsi_2_VtxProb_var));
    RooDataSet Phi_VtxProb_set("Phi_VtxProb_set", "Phi_VtxProb_set", RooArgList(Phi_VtxProb_var));
    RooDataSet Pri_VtxProb_set("Pri_VtxProb_set", "Pri_VtxProb_set", RooArgList(Pri_VtxProb_var));

   gStyle->SetOptStat(0);
   gStyle->SetCanvasColor(33);
   gStyle->SetFrameFillColor(18);
   Int_t cancolor = 17;
   auto h2 = new TH2F("h2", "Jpsi_1_mass and Jpsi_2_mass", 20, 2.9, 3.3, 20, 2.9, 3.3);
   h2->SetFillColor(46);
   auto h3 = new TH2F("h3", "Jpsi_1_mass and Phi_mass", 20, 2.9, 3.3, 20, 0.99, 1.07);
   h3->SetFillColor(46);
   auto h4 = new TH2F("h4", "Jpsi_2_mass and Phi_mass", 20, 2.9, 3.3, 20, 0.99, 1.07);
   h4->SetFillColor(46);

    RooRealVar m_Jpsi_1("m_Jpsi_1", "J/#psi_{1} invariant mass", 3.097, 2.9, 3.3);
    RooRealVar m_Jpsi_2("m_Jpsi_2", "J/#psi_{2} invariant mass", 3.097, 2.9, 3.3);
    RooRealVar m_Phi("m_Phi", "#Phi invariant mass", 1.019, 0.99, 1.07);

    RooDataSet data("data", "Dataset with m_Jpsi_1, m_Jpsi_2, m_Phi", RooArgSet(m_Jpsi_1, m_Jpsi_2, m_Phi));

   // J/psi_1 的信号模型（Crystal Ball + Gaussian）
   RooRealVar mean_Jpsi_1("mean_Jpsi_1", "Mean of J/#psi_{1}", 3.097, 3.05, 3.15);
   RooRealVar sigma_CB_1("sigma_CB_1", "CB Sigma of J/#psi_{1}", 0.02, 0.005, 0.05);
   RooRealVar alpha_1("alpha_1", "alpha of J/#psi_{1}", 1.5, 0.5, 5.0);
   RooRealVar n_1("n_1", "n of J/#psi_{1}", 2.0, 1.0, 10.0);
   RooCBShape cb_Jpsi_1("cb_Jpsi_1", "CB function for J/#psi_{1}", m_Jpsi_1, mean_Jpsi_1, sigma_CB_1, alpha_1, n_1);

   RooRealVar sigma_Gauss_1("sigma_Gauss_1", "Gaussian Sigma of J/#psi_{1}", 0.03, 0.01, 0.1);
   RooGaussian gauss_Jpsi_1("gauss_Jpsi_1", "Gaussian for J/#psi_{1}", m_Jpsi_1, mean_Jpsi_1, sigma_Gauss_1);

   RooRealVar frac_1("frac_1", "Fraction of CB in J/#psi_{1}", 0.7, 0.0, 1.0);
   RooAddPdf signal_Jpsi_1("signal_Jpsi_1", "Signal PDF for J/#psi_{1}", RooArgList(cb_Jpsi_1, gauss_Jpsi_1), frac_1);

   // J/psi_1 的背景模型（指数函数）
   RooRealVar lambda_1("lambda_1", "Background slope for J/#psi_{1}", -0.5, -10.0, 0.0);
   RooExponential bkg_Jpsi_1("bkg_Jpsi_1", "Background PDF for J/#psi_{1}", m_Jpsi_1, lambda_1);

   // J/psi_2 的信号模型（Crystal Ball + Gaussian）
   RooRealVar mean_Jpsi_2("mean_Jpsi_2", "Mean of J/#psi_{2}", 3.097, 3.05, 3.15);
   RooRealVar sigma_CB_2("sigma_CB_2", "CB Sigma of J/#psi_{2}", 0.02, 0.005, 0.05);
   RooRealVar alpha_2("alpha_2", "alpha of J/#psi_{2}", 1.5, 0.5, 5.0);
   RooRealVar n_2("n_2", "n of J/#psi_{2}", 2.0, 1.0, 10.0);
   RooCBShape cb_Jpsi_2("cb_Jpsi_2", "CB function for J/#psi_{2}", m_Jpsi_2, mean_Jpsi_2, sigma_CB_2, alpha_2, n_2);

   RooRealVar sigma_Gauss_2("sigma_Gauss_2", "Gaussian Sigma of J/#psi_{2}", 0.03, 0.01, 0.1);
   RooGaussian gauss_Jpsi_2("gauss_Jpsi_2", "Gaussian for J/#psi_{2}", m_Jpsi_2, mean_Jpsi_2, sigma_Gauss_2);

   RooRealVar frac_2("frac_2", "Fraction of CB in J/#psi_{2}", 0.7, 0.0, 1.0);
   RooAddPdf signal_Jpsi_2("signal_Jpsi_2", "Signal PDF for J/#psi_{2}", RooArgList(cb_Jpsi_2, gauss_Jpsi_2), frac_2);

   // J/psi_2 的背景模型（指数函数）
   RooRealVar lambda_2("lambda_2", "Background slope for J/#psi_{2}", -0.5, -10.0, 0.0);
   RooExponential bkg_Jpsi_2("bkg_Jpsi_2", "Background PDF for J/#psi_{2}", m_Jpsi_2, lambda_2);

   // Phi 的信号模型（高斯分布）
   RooRealVar mean_Phi("mean_Phi", "Mean of #Phi", 1.020, 1.01, 1.03); // MeV单位
   RooRealVar sigma_Phi("sigma_Phi", "Sigma of #Phi", 0.003, 0.001, 0.01); // 3.1 MeV宽度
   RooGaussian signal_Phi("signal_Phi", "Signal PDF for #Phi", m_Phi, mean_Phi, sigma_Phi);

   // // Phi 的背景模型（4阶多项式）
   // RooRealVar c0("c0", "c0", 1.0, -1000.0, 1000.0);
   // RooRealVar c1("c1", "c1", 0.1, -1000.0, 1000.0);
   // RooRealVar c2("c2", "c2", 0.01, -1000.0, 1000.0);
   // // RooRealVar c3("c3", "c3", 0.001, -1000.0, 1000.0);
   // // RooRealVar c4("c4", "c4", 0.0001, -10000.0, 10000.0);
   // // RooPolynomial bkg_Phi("bkg_Phi", "Background PDF for #Phi", m_Phi, RooArgList(c0, c1, c2, c3, c4));
   // RooPolynomial bkg_Phi("bkg_Phi", "Background PDF for #Phi", m_Phi, RooArgList(c0, c1, c2));

   // Phi 的背景模型（误差函数）
   RooRealVar erf_c0("erf_c0", "erf_c0", 1.0, 0.0, 10.0);         // 整体归一化参数
   RooRealVar erf_c1("erf_c1", "erf_c1", 1.02, 0.98, 1.05);       // 位置参数（阈值）
   RooRealVar erf_c2("erf_c2", "erf_c2", 0.005, 0.001, 0.02);     // 宽度参数
   RooRealVar erf_c3("erf_c3", "erf_c3", 0.0, -0.5, 0.5);         // 线性分量的斜率
   RooRealVar erf_c4("erf_c4", "erf_c4", 1.0, 0.0, 10.0);         // 误差函数的比例因子
   RooRealVar erf_c5("erf_c5", "二次项系数", 0.0, -5.0, 5.0);
   // 背景模型：p0 + p1*(m-m0) + p2*erf((m-m0)/p3)
   // 使用线性项和误差函数的组合
   RooGenericPdf bkg_Phi("bkg_Phi", "Phi介子背景PDF", 
      "erf_c0 + erf_c3*(m_Phi - erf_c1) + erf_c5*(m_Phi - erf_c1)*(m_Phi - erf_c1) + erf_c4*TMath::Erf((m_Phi - erf_c1)/erf_c2)",
      RooArgSet(m_Phi, erf_c0, erf_c1, erf_c2, erf_c3, erf_c4, erf_c5));


   // // 定义三个子组件的比例
   // RooRealVar nsig_Jpsi_1("nsig_Jpsi_1", "Number of J/#psi_{1} signal events", 1000, 0, 10000);
   // RooRealVar nbkg_Jpsi_1("nbkg_Jpsi_1", "Number of J/#psi_{1} background events", 500, 0, 5000);
   // RooRealVar nsig_Jpsi_2("nsig_Jpsi_2", "Number of J/#psi_{2} signal events", 1000, 0, 10000);
   // RooRealVar nbkg_Jpsi_2("nbkg_Jpsi_2", "Number of J/#psi_{2} background events", 500, 0, 5000);
   // RooRealVar nsig_Phi("nsig_Phi", "Number of #Phi signal events", 1000, 0, 10000);
   // RooRealVar nbkg_Phi("nbkg_Phi", "Number of #Phi background events", 500, 0, 5000);

   // // 组合 J/psi_1 的信号和背景
   // RooAddPdf model_Jpsi_1("model_Jpsi_1", "Signal + Background for J/#psi_{1}", 
   //                      RooArgList(signal_Jpsi_1, bkg_Jpsi_1), 
   //                      RooArgList(nsig_Jpsi_1, nbkg_Jpsi_1));

   // // 组合 J/psi_2 的信号和背景
   // RooAddPdf model_Jpsi_2("model_Jpsi_2", "Signal + Background for J/#psi_{2}", 
   //                      RooArgList(signal_Jpsi_2, bkg_Jpsi_2), 
   //                      RooArgList(nsig_Jpsi_2, nbkg_Jpsi_2));

   // // 组合 Phi 的信号和背景
   // RooAddPdf model_Phi("model_Phi", "Signal + Background for #Phi", 
   //                   RooArgList(signal_Phi, bkg_Phi), 
   //                   RooArgList(nsig_Phi, nbkg_Phi));

   // SSS: signal+signal+signal
   RooRealVar yield_SSS("yield_SSS", "Yield of SSS", 100, 0, 10000);
   RooProdPdf pdf_SSS("pdf_SSS", "Signal+Signal+Signal PDF", 
                  RooArgList(signal_Jpsi_1, signal_Jpsi_2, signal_Phi));

   // SSB: signal+signal+background
   RooRealVar yield_SSB("yield_SSB", "Yield of SSB", 50, 0, 10000);
   RooProdPdf pdf_SSB("pdf_SSB", "Signal+Signal+Background PDF", 
                  RooArgList(signal_Jpsi_1, signal_Jpsi_2, bkg_Phi));

   // SBS: signal+background+signal
   RooRealVar yield_SBS("yield_SBS", "Yield of SBS", 50, 0, 10000);
   RooProdPdf pdf_SBS("pdf_SBS", "Signal+Background+Signal PDF", 
                  RooArgList(signal_Jpsi_1, bkg_Jpsi_2, signal_Phi));

   // BSS: background+signal+signal
   RooRealVar yield_BSS("yield_BSS", "Yield of BSS", 50, 0, 10000);
   RooProdPdf pdf_BSS("pdf_BSS", "Background+Signal+Signal PDF", 
                  RooArgList(bkg_Jpsi_1, signal_Jpsi_2, signal_Phi));

   // SBB: signal+background+background
   RooRealVar yield_SBB("yield_SBB", "Yield of SBB", 25, 0, 10000);
   RooProdPdf pdf_SBB("pdf_SBB", "Signal+Background+Background PDF", 
                  RooArgList(signal_Jpsi_1, bkg_Jpsi_2, bkg_Phi));

   // BSB: background+signal+background
   RooRealVar yield_BSB("yield_BSB", "Yield of BSB", 25, 0, 10000);
   RooProdPdf pdf_BSB("pdf_BSB", "Background+Signal+Background PDF", 
                  RooArgList(bkg_Jpsi_1, signal_Jpsi_2, bkg_Phi));

   // BBS: background+background+signal
   RooRealVar yield_BBS("yield_BBS", "Yield of BBS", 25, 0, 10000);
   RooProdPdf pdf_BBS("pdf_BBS", "Background+Background+Signal PDF", 
                  RooArgList(bkg_Jpsi_1, bkg_Jpsi_2, signal_Phi));

   // BBB: background+background+background
   RooRealVar yield_BBB("yield_BBB", "Yield of BBB", 10, 0, 10000);
   RooProdPdf pdf_BBB("pdf_BBB", "Background+Background+Background PDF", 
                  RooArgList(bkg_Jpsi_1, bkg_Jpsi_2, bkg_Phi));

   // 组合8种可能性为一个完整模型
   RooAddPdf totalModel("totalModel", "Total PDF", 
                     RooArgList(pdf_SSS, pdf_SSB, pdf_SBS, pdf_BSS, pdf_SBB, pdf_BSB, pdf_BBS, pdf_BBB),
                     RooArgList(yield_SSS, yield_SSB, yield_SBS, yield_BSS, yield_SBB, yield_BSB, yield_BBS, yield_BBB));


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
         printf(">>> Processing entry %lld <<<\n", jentry);
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
      double Pri_VtxProb_cut = 0.01;
      double Phi_VtxProb_cut = 0.05;

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

            #ifdef CUT_PHI_VTXPROB
            if (Phi_VtxProb->at(iCand) < Phi_VtxProb_cut){
               passCut = false;
            }
            #endif

            #ifdef CUT_PHI_CTAU
            if (Phi_ctau->at(iCand) < -0.02 || Phi_ctau->at(iCand) > 0.02){
               passCut = false;
            }
            #endif

            #ifdef Stefanos_CUT
            if(Jpsi_1_pt->at(iCand) < 6.0 || Jpsi_2_pt->at(iCand) < 6.0){
               passCut = false;
            }

            if(Jpsi_1_VtxProb->at(iCand) < 0.1 || Jpsi_2_VtxProb->at(iCand) < 0.1){
               passCut = false;
            }

            double Jpsi_1_E = sqrt(Jpsi_1_mass->at(iCand) * Jpsi_1_mass->at(iCand) + Jpsi_1_px->at(iCand) * Jpsi_1_px->at(iCand) + Jpsi_1_py->at(iCand) * Jpsi_1_py->at(iCand) + Jpsi_1_pz->at(iCand) * Jpsi_1_pz->at(iCand));
            double Jpsi_2_E = sqrt(Jpsi_2_mass->at(iCand) * Jpsi_2_mass->at(iCand) + Jpsi_2_px->at(iCand) * Jpsi_2_px->at(iCand) + Jpsi_2_py->at(iCand) * Jpsi_2_py->at(iCand) + Jpsi_2_pz->at(iCand) * Jpsi_2_pz->at(iCand));
            double Jpsi_1_gamma = 0.5 * log((Jpsi_1_E + Jpsi_1_pz->at(iCand)) / (Jpsi_1_E - Jpsi_1_pz->at(iCand)));
            double Jpsi_2_gamma = 0.5 * log((Jpsi_2_E + Jpsi_2_pz->at(iCand)) / (Jpsi_2_E - Jpsi_2_pz->at(iCand)));
            if(fabs(Jpsi_1_gamma) > 2.4 || fabs(Jpsi_2_gamma) > 2.4){
               passCut = false;
            }
            
            //JpsiJpsi四muon顶点拟合>0.01在ntuple里做过了

            
            Jpsi_1_Lxy->push_back(Jpsi_1_ctau->at(iCand) * Jpsi_1_pt->at(iCand) / Jpsi_1_mass->at(iCand));
            Jpsi_2_Lxy->push_back(Jpsi_2_ctau->at(iCand) * Jpsi_2_pt->at(iCand) / Jpsi_2_mass->at(iCand));
            Phi_Lxy->push_back(Phi_ctau->at(iCand) * Phi_pt->at(iCand) / Phi_mass->at(iCand));
            // if(Jpsi_1_Lxy->at(iCand) > 0.1 || Jpsi_2_Lxy->at(iCand) > 0.1 || Phi_Lxy->at(iCand) > 0.1){
            //    passCut = false;
            // }

            // double Jpsi_1_Lxy = Jpsi_1_ctau->at(iCand) * Jpsi_1_pt->at(iCand) / Jpsi_1_mass->at(iCand);
            // double Jpsi_2_Lxy = Jpsi_2_ctau->at(iCand) * Jpsi_2_pt->at(iCand) / Jpsi_2_mass->at(iCand);
            // double Phi_Lxy = Phi_ctau->at(iCand) * Phi_pt->at(iCand) / Phi_mass->at(iCand);
            // if(Jpsi_1_Lxy > 0.01 || Jpsi_2_Lxy > 0.01 || Phi_Lxy > 0.01){
            //    passCut = false;
            // }
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

            tempCand.SetScore(temp_pt_abs);
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

         // Jpsi_1_ctau_error_var.setVal(Jpsi_1_ctauErr->at(cand->GetId()));
         // Jpsi_2_ctau_error_var.setVal(Jpsi_2_ctauErr->at(cand->GetId()));
         // Phi_ctau_error_var.setVal(Phi_ctauErr->at(cand->GetId()));
         // Pri_ctau_error_var.setVal(Pri_ctauErr->at(cand->GetId()));

         // Jpsi_1_ctau_error_set.add(RooArgSet(Jpsi_1_ctau_error_var));
         // Jpsi_2_ctau_error_set.add(RooArgSet(Jpsi_2_ctau_error_var));
         // Phi_ctau_error_set.add(RooArgSet(Phi_ctau_error_var));
         // Pri_ctau_error_set.add(RooArgSet(Pri_ctau_error_var));

         Jpsi_1_Lxy_var.setVal(Jpsi_1_Lxy->at(cand->GetId()));
         Jpsi_2_Lxy_var.setVal(Jpsi_2_Lxy->at(cand->GetId()));
         Phi_Lxy_var.setVal(Phi_Lxy->at(cand->GetId()));

         Jpsi_1_Lxy_set.add(RooArgSet(Jpsi_1_Lxy_var));
         Jpsi_2_Lxy_set.add(RooArgSet(Jpsi_2_Lxy_var));
         Phi_Lxy_set.add(RooArgSet(Phi_Lxy_var));

         Jpsi_1_VtxProb_var.setVal(Jpsi_1_VtxProb->at(cand->GetId()));
         Jpsi_2_VtxProb_var.setVal(Jpsi_2_VtxProb->at(cand->GetId()));
         Phi_VtxProb_var.setVal(Phi_VtxProb->at(cand->GetId()));
         Pri_VtxProb_var.setVal(Pri_VtxProb->at(cand->GetId()));

         Jpsi_1_VtxProb_set.add(RooArgSet(Jpsi_1_VtxProb_var));
         Jpsi_2_VtxProb_set.add(RooArgSet(Jpsi_2_VtxProb_var));
         Phi_VtxProb_set.add(RooArgSet(Phi_VtxProb_var));
         Pri_VtxProb_set.add(RooArgSet(Pri_VtxProb_var));

         double jpsi1_mass = Jpsi_1_mass->at(cand->GetId()); // 获取 Jpsi_1_mass 的值
         double jpsi2_mass = Jpsi_2_mass->at(cand->GetId()); // 获取 Jpsi_2_mass 的值
         double phi_mass = Phi_mass->at(cand->GetId());      // 获取 Phi_mass 的值
         h2->Fill(jpsi1_mass, jpsi2_mass);      // 填充到二维直方图
         h3->Fill(jpsi1_mass, phi_mass);         // 填充到二维直方图
         h4->Fill(jpsi2_mass, phi_mass);         // 填充到二维直方图

         m_Jpsi_1.setVal(Jpsi_1_mass->at(cand->GetId()));
         m_Jpsi_2.setVal(Jpsi_2_mass->at(cand->GetId()));
         m_Phi.setVal(Phi_mass->at(cand->GetId()));
         data.add(RooArgSet(m_Jpsi_1, m_Jpsi_2, m_Phi));
      }
      ClearBranches();
      CandList.clear();
   }

   // Draw the histograms.
   TCanvas *Jpsi1MassCanvas = new TCanvas("Jpsi1MassCanvas", "Jpsi1MassCanvas", 1600, 1200);
   Jpsi1MassCanvas->Divide(2,2);
   RooPlot* Jpsi_1_mass_frame = Jpsi_1_mass_var.frame(nBins);
   RooPlot* Jpsi_2_mass_frame = Jpsi_2_mass_var.frame(nBins);
   RooPlot* Phi_mass_frame = Phi_mass_var.frame(nBins);
   RooPlot* Pri_mass_frame = Pri_mass_var.frame(nBins);
   Jpsi_1_mass_set_multi.plotOn(Jpsi_1_mass_frame);
   Jpsi_2_mass_set_multi.plotOn(Jpsi_2_mass_frame);
   Phi_mass_set_multi.plotOn(Phi_mass_frame);
   Pri_mass_set_multi.plotOn(Pri_mass_frame);
   Jpsi1MassCanvas->cd(1);
   Jpsi_1_mass_frame->Draw();
   Jpsi1MassCanvas->cd(2);
   Jpsi_2_mass_frame->Draw();
   Jpsi1MassCanvas->cd(3);
   Phi_mass_frame->Draw();
   Jpsi1MassCanvas->cd(4);
   Pri_mass_frame->Draw();
   Jpsi1MassCanvas->SaveAs("secCut_multi_by_pT.pdf");
   Jpsi1MassCanvas->SaveAs("secCut_multi_by_pT.png");

   // Draw the histograms for the non-overlap candidates.
   TCanvas *Jpsi2MassCanvas = new TCanvas("Jpsi2MassCanvas", "Jpsi2MassCanvas", 1600, 1200);
   Jpsi2MassCanvas->Divide(2,2);
   RooPlot* Jpsi_1_mass_frame_nonOverlap = Jpsi_1_mass_var.frame(nBins);
   RooPlot* Jpsi_2_mass_frame_nonOverlap = Jpsi_2_mass_var.frame(nBins);
   RooPlot* Phi_mass_frame_nonOverlap = Phi_mass_var.frame(nBins);
   RooPlot* Pri_mass_frame_nonOverlap = Pri_mass_var.frame(nBins);
   Jpsi_1_mass_set.plotOn(Jpsi_1_mass_frame_nonOverlap);
   Jpsi_2_mass_set.plotOn(Jpsi_2_mass_frame_nonOverlap);
   Phi_mass_set.plotOn(Phi_mass_frame_nonOverlap);
   Pri_mass_set.plotOn(Pri_mass_frame_nonOverlap);
   Jpsi2MassCanvas->cd(1);
   Jpsi_1_mass_frame_nonOverlap->Draw();
   Jpsi2MassCanvas->cd(2);
   Jpsi_2_mass_frame_nonOverlap->Draw();
   Jpsi2MassCanvas->cd(3);
   Phi_mass_frame_nonOverlap->Draw();
   Jpsi2MassCanvas->cd(4);
   Pri_mass_frame_nonOverlap->Draw();
   Jpsi2MassCanvas->SaveAs("secCut_nonOverlap_by_pT.pdf");
   Jpsi2MassCanvas->SaveAs("secCut_nonOverlap_by_pT.png");

   // Draw the histograms for the ctau
   TCanvas *PhiMassCanvas = new TCanvas("PhiMassCanvas", "PhiMassCanvas", 1600, 1200);
   PhiMassCanvas->Divide(2,2);
   RooPlot* Jpsi_1_ctau_frame = Jpsi_1_ctau_var.frame(nBins);
   RooPlot* Jpsi_2_ctau_frame = Jpsi_2_ctau_var.frame(nBins);
   RooPlot* Phi_ctau_frame = Phi_ctau_var.frame(nBins);
   RooPlot* Pri_ctau_frame = Pri_ctau_var.frame(nBins);
   Jpsi_1_ctau_set.plotOn(Jpsi_1_ctau_frame);
   Jpsi_2_ctau_set.plotOn(Jpsi_2_ctau_frame);
   Phi_ctau_set.plotOn(Phi_ctau_frame);
   Pri_ctau_set.plotOn(Pri_ctau_frame);
   PhiMassCanvas->cd(1);
   gPad->SetLogy();
   Jpsi_1_ctau_frame->Draw();
   PhiMassCanvas->cd(2);
   gPad->SetLogy();
   Jpsi_2_ctau_frame->Draw();
   PhiMassCanvas->cd(3);
   gPad->SetLogy();
   Phi_ctau_frame->Draw();
   PhiMassCanvas->cd(4);
   gPad->SetLogy();
   Pri_ctau_frame->Draw();
   PhiMassCanvas->SaveAs("secCut_ctau_by_pT.pdf");
   PhiMassCanvas->SaveAs("secCut_ctau_by_pT.png");

   // Draw the histograms for the VtxProb
   TCanvas *PriMassCanvas = new TCanvas("PriMassCanvas", "PriMassCanvas", 1600, 1200);
   PriMassCanvas->Divide(2,2);
   RooPlot* Jpsi_1_VtxProb_frame = Jpsi_1_VtxProb_var.frame(nBins);
   RooPlot* Jpsi_2_VtxProb_frame = Jpsi_2_VtxProb_var.frame(nBins);
   RooPlot* Phi_VtxProb_frame = Phi_VtxProb_var.frame(nBins);
   RooPlot* Pri_VtxProb_frame = Pri_VtxProb_var.frame(nBins);
   Jpsi_1_VtxProb_set.plotOn(Jpsi_1_VtxProb_frame);
   Jpsi_2_VtxProb_set.plotOn(Jpsi_2_VtxProb_frame);
   Phi_VtxProb_set.plotOn(Phi_VtxProb_frame);
   Pri_VtxProb_set.plotOn(Pri_VtxProb_frame);
   PriMassCanvas->cd(1);
   // gPad->SetLogy();
   Jpsi_1_VtxProb_frame->Draw();
   PriMassCanvas->cd(2);
   // gPad->SetLogy();
   Jpsi_2_VtxProb_frame->Draw();
   PriMassCanvas->cd(3);
   // gPad->SetLogy();
   Phi_VtxProb_frame->Draw();
   PriMassCanvas->cd(4);
   // gPad->SetLogy();
   Pri_VtxProb_frame->Draw();
   PriMassCanvas->SaveAs("secCut_VtxProb_by_pT.pdf");
   PriMassCanvas->SaveAs("secCut_VtxProb_by_pT.png");

   // Draw the histograms for the ctau error
   // TCanvas *c5 = new TCanvas("c5", "c5", 1600, 1200);
   // c5->Divide(2,2);
   // RooPlot* Jpsi_1_ctau_error_frame = Jpsi_1_ctau_error_var.frame(nBins);
   // RooPlot* Jpsi_2_ctau_error_frame = Jpsi_2_ctau_error_var.frame(nBins);
   // RooPlot* Phi_ctau_error_frame = Phi_ctau_error_var.frame(nBins);
   // RooPlot* Pri_ctau_error_frame = Pri_ctau_error_var.frame(nBins);
   // Jpsi_1_ctau_error_set.plotOn(Jpsi_1_ctau_error_frame);
   // Jpsi_2_ctau_error_set.plotOn(Jpsi_2_ctau_error_frame);
   // Phi_ctau_error_set.plotOn(Phi_ctau_error_frame);
   // Pri_ctau_error_set.plotOn(Pri_ctau_error_frame);
   // c5->cd(1);
   // gPad->SetLogy();
   // Jpsi_1_ctau_error_frame->Draw();
   // c5->cd(2);
   // gPad->SetLogy();
   // Jpsi_2_ctau_error_frame->Draw();
   // c5->cd(3);
   // gPad->SetLogy();
   // Phi_ctau_error_frame->Draw();
   // c5->cd(4);
   // gPad->SetLogy();
   // Pri_ctau_error_frame->Draw();
   // c5->SaveAs("secCut_ctau_error.pdf");
   // c5->SaveAs("secCut_ctau_error.png");

   // Draw the histograms for the Lxy
   TCanvas *c6 = new TCanvas("c6", "c6", 1600, 1200);
   c6->Divide(2,2);
   RooPlot* Jpsi_1_Lxy_frame = Jpsi_1_Lxy_var.frame(nBins);
   RooPlot* Jpsi_2_Lxy_frame = Jpsi_2_Lxy_var.frame(nBins);
   RooPlot* Phi_Lxy_frame = Phi_Lxy_var.frame(nBins);
   Jpsi_1_Lxy_set.plotOn(Jpsi_1_Lxy_frame);
   Jpsi_2_Lxy_set.plotOn(Jpsi_2_Lxy_frame);
   Phi_Lxy_set.plotOn(Phi_Lxy_frame);
   c6->cd(1);
   gPad->SetLogy();
   Jpsi_1_Lxy_frame->Draw();
   c6->cd(2);
   gPad->SetLogy();
   Jpsi_2_Lxy_frame->Draw();
   c6->cd(3);
   gPad->SetLogy();
   Phi_Lxy_frame->Draw();
   c6->SaveAs("secCut_Lxy_by_pT_raw.pdf");
   c6->SaveAs("secCut_Lxy_by_pT_raw.png");



   

   // auto lego = new TCanvas("lego", "lego options", 150, 150, 800, 600);
   // lego->Divide(2, 2);
   // lego->SetFillColor(cancolor);
   // lego->cd(1);
   // h2->Draw("lego");
   // lego->SaveAs("secCut_Mass_lego.pdf");
   // lego->SaveAs("secCut_Mass_lego.png");

   auto surf = new TCanvas("surfopt", "surface options", 200, 200, 800, 600);
   surf->Divide(2, 2);
   surf->SetFillColor(cancolor);
   surf->cd(1);
   h2->Draw("surf3");
   surf->cd(2);
   h3->Draw("surf3");
   surf->cd(3);
   h4->Draw("surf3");
   surf->SaveAs("secCut_Mass_surf.pdf");
   surf->SaveAs("secCut_Mass_surf.png");

   // 执行三维拟合
   //RooFitResult* fitResult = model.fitTo(data, RooFit::Extended(true), RooFit::Save(), RooFit::PrintLevel(1));
   //RooFitResult* fitResult = extendedModel.fitTo(data, RooFit::Save(), RooFit::PrintLevel(1));
   RooFitResult* fitResult = totalModel.fitTo(data, RooFit::Extended(true), RooFit::Save(), RooFit::PrintLevel(1));

   // 创建新的画布展示拟合结果
   TCanvas *fitCanvas = new TCanvas("fitCanvas", "Fit Results", 1800, 600);
   fitCanvas->Divide(3, 1);

   // 创建分布的帧
   RooPlot* frameMJpsi1 = m_Jpsi_1.frame(RooFit::Title("J/#psi_{1} Mass Distribution"), RooFit::Bins(40));
   RooPlot* frameMJpsi2 = m_Jpsi_2.frame(RooFit::Title("J/#psi_{2} Mass Distribution"), RooFit::Bins(40));
   RooPlot* frameMPhi = m_Phi.frame(RooFit::Title("#Phi Mass Distribution"), RooFit::Bins(40));

   // 绘制数据
   data.plotOn(frameMJpsi1);
   data.plotOn(frameMJpsi2);
   data.plotOn(frameMPhi);

   // 在每个变量上投影PDF
   // model.plotOn(frameMJpsi1, RooFit::Components("model_Jpsi_1"), RooFit::LineColor(kRed));
   // model.plotOn(frameMJpsi1, RooFit::Components("signal_Jpsi_1"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   // model.plotOn(frameMJpsi1, RooFit::Components("bkg_Jpsi_1"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   // model.plotOn(frameMJpsi2, RooFit::Components("model_Jpsi_2"), RooFit::LineColor(kRed));
   // model.plotOn(frameMJpsi2, RooFit::Components("signal_Jpsi_2"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   // model.plotOn(frameMJpsi2, RooFit::Components("bkg_Jpsi_2"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   // model.plotOn(frameMPhi, RooFit::Components("model_Phi"), RooFit::LineColor(kRed));
   // model.plotOn(frameMPhi, RooFit::Components("signal_Phi"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   // model.plotOn(frameMPhi, RooFit::Components("bkg_Phi"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   totalModel.plotOn(frameMJpsi1, RooFit::LineColor(kRed));
   totalModel.plotOn(frameMJpsi1, RooFit::Components("pdf_SSS"), 
                  RooFit::LineColor(kRed), RooFit::FillColor(kRed-4), 
                  RooFit::FillStyle(3004), RooFit::DrawOption("F"));
   totalModel.plotOn(frameMJpsi1, RooFit::Components("pdf_SSS,pdf_SSB,pdf_SBS,pdf_SBB"), 
                  RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   totalModel.plotOn(frameMJpsi1, RooFit::Components("pdf_BSS,pdf_BSB,pdf_BBS,pdf_BBB"), 
                  RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   totalModel.plotOn(frameMJpsi2, RooFit::LineColor(kRed));
   totalModel.plotOn(frameMJpsi2, RooFit::Components("pdf_SSS"), 
                  RooFit::LineColor(kRed), RooFit::FillColor(kRed-4), 
                  RooFit::FillStyle(3004), RooFit::DrawOption("F"));
   totalModel.plotOn(frameMJpsi2, RooFit::Components("pdf_SSS,pdf_SSB,pdf_BSS,pdf_BSB"), 
                  RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   totalModel.plotOn(frameMJpsi2, RooFit::Components("pdf_SBS,pdf_SBB,pdf_BBS,pdf_BBB"), 
                  RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   totalModel.plotOn(frameMPhi, RooFit::LineColor(kRed));
   totalModel.plotOn(frameMPhi, RooFit::Components("pdf_SSS"), 
                 RooFit::LineColor(kRed), RooFit::FillColor(kRed-4), 
                 RooFit::FillStyle(3004), RooFit::DrawOption("F"));
   totalModel.plotOn(frameMPhi, RooFit::Components("pdf_SSS,pdf_SBS,pdf_BSS,pdf_BBS"), 
                  RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
   totalModel.plotOn(frameMPhi, RooFit::Components("pdf_SSB,pdf_SBB,pdf_BSB,pdf_BBB"), 
                  RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

   // 添加拟合结果文本框
   fitCanvas->cd(1);
   frameMJpsi1->Draw();
   TPaveText *textJpsi1 = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
   textJpsi1->AddText(Form("Mean = %.3f #pm %.3f GeV", mean_Jpsi_1.getVal(), mean_Jpsi_1.getError()));
   textJpsi1->AddText(Form("CB #sigma = %.3f #pm %.3f GeV", sigma_CB_1.getVal(), sigma_CB_1.getError()));
   textJpsi1->AddText(Form("Gauss #sigma = %.3f #pm %.3f GeV", sigma_Gauss_1.getVal(), sigma_Gauss_1.getError()));
   //textJpsi1->AddText(Form("N_{sss} = %.0f #pm %.0f", nsig_Jpsi_1_SSS.getVal(), nsig_Jpsi_1_SSS.getError()));
   textJpsi1->AddText(Form("SSS events = %.0f #pm %.0f", 
      yield_SSS.getVal(), sqrt(pow(yield_SSS.getError(),2))));
   textJpsi1->SetFillColor(0);
   textJpsi1->SetBorderSize(1);
   textJpsi1->Draw();

   fitCanvas->cd(2);
   frameMJpsi2->Draw();
   TPaveText *textJpsi2 = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
   textJpsi2->AddText(Form("Mean = %.3f #pm %.3f GeV", mean_Jpsi_2.getVal(), mean_Jpsi_2.getError()));
   textJpsi2->AddText(Form("CB #sigma = %.3f #pm %.3f GeV", sigma_CB_2.getVal(), sigma_CB_2.getError()));
   textJpsi2->AddText(Form("Gauss #sigma = %.3f #pm %.3f GeV", sigma_Gauss_2.getVal(), sigma_Gauss_2.getError()));
   //textJpsi2->AddText(Form("N_{sss} = %.0f #pm %.0f", nsig_Jpsi_2_SSS.getVal(), nsig_Jpsi_2_SSS.getError()));
   textJpsi2->AddText(Form("SSS events = %.0f #pm %.0f", 
      yield_SSS.getVal(), sqrt(pow(yield_SSS.getError(),2))));
   textJpsi2->SetFillColor(0);
   textJpsi2->SetBorderSize(1);
   textJpsi2->Draw();

   fitCanvas->cd(3);
   frameMPhi->Draw();
   TPaveText *textPhi = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
   textPhi->AddText(Form("Mean = %.3f #pm %.3f GeV", mean_Phi.getVal(), mean_Phi.getError()));
   textPhi->AddText(Form("Sigma = %.3f #pm %.3f GeV", sigma_Phi.getVal(), sigma_Phi.getError()));
   //textPhi->AddText(Form("N_{sss} = %.0f #pm %.0f", nsig_Phi_SSS.getVal(), nsig_Phi_SSS.getError()));
   textPhi->AddText(Form("SSS events = %.0f #pm %.0f", 
      yield_SSS.getVal(), sqrt(pow(yield_SSS.getError(),2))));
   textPhi->SetFillColor(0);
   textPhi->SetBorderSize(1);
   textPhi->Draw();

   // 保存拟合结果
   fitCanvas->SaveAs("secCut_Fit_Results.pdf");
   fitCanvas->SaveAs("secCut_Fit_Results.png");

   // 打印拟合参数
   // std::cout << "\n===== Fit Results =====" << std::endl;
   // std::cout << "J/psi_1 Mean: " << mean_Jpsi_1.getVal() << " ± " << mean_Jpsi_1.getError() << " GeV" << std::endl;
   // std::cout << "J/psi_2 Mean: " << mean_Jpsi_2.getVal() << " ± " << mean_Jpsi_2.getError() << " GeV" << std::endl;
   // std::cout << "Phi Mean: " << mean_Phi.getVal() << " ± " << mean_Phi.getError() << " GeV" << std::endl;
   // std::cout << "Phi Width: " << sigma_Phi.getVal()*1000 << " ± " << sigma_Phi.getError()*1000 << " MeV" << std::endl;

   // 打印每个model的参数
   // RooArgSet* components = model.getComponents();
   // std::cout << "\n===== Model Components =====" << std::endl;
   // components->Print("v");

   std::cout << "\n===== Component Yields =====" << std::endl;
   std::cout << "SSS: " << yield_SSS.getVal() << " ± " << yield_SSS.getError() << std::endl;
   std::cout << "SSB: " << yield_SSB.getVal() << " ± " << yield_SSB.getError() << std::endl;
   std::cout << "SBS: " << yield_SBS.getVal() << " ± " << yield_SBS.getError() << std::endl;
   std::cout << "BSS: " << yield_BSS.getVal() << " ± " << yield_BSS.getError() << std::endl;
   std::cout << "SBB: " << yield_SBB.getVal() << " ± " << yield_SBB.getError() << std::endl;
   std::cout << "BSB: " << yield_BSB.getVal() << " ± " << yield_BSB.getError() << std::endl;
   std::cout << "BBS: " << yield_BBS.getVal() << " ± " << yield_BBS.getError() << std::endl;
   std::cout << "BBB: " << yield_BBB.getVal() << " ± " << yield_BBB.getError() << std::endl;

   // 计算总信号数和背景数
   double total_events = yield_SSS.getVal() + yield_SSB.getVal() + yield_SBS.getVal() + yield_BSS.getVal() +
                     yield_SBB.getVal() + yield_BSB.getVal() + yield_BBS.getVal() + yield_BBB.getVal();
                     
   double pure_signal_events = yield_SSS.getVal();
   double signal_err = yield_SSS.getError();

   std::cout << "\n===== Summary =====" << std::endl;
   std::cout << "Total events: " << total_events << std::endl;
   std::cout << "Pure signal events (SSS): " << pure_signal_events << " ± " << signal_err << std::endl;
   std::cout << "Pure signal fraction: " << pure_signal_events/total_events*100 << "%" << std::endl;

   // 1. 首先保存包含信号+背景的拟合结果
   RooFitResult* fitResult_SB = totalModel.fitTo(data, RooFit::Extended(true), 
   RooFit::Save(), RooFit::PrintLevel(1));
   double nll_SB = fitResult_SB->minNll(); // 带信号+背景的负对数似然值

   // 2. 强制信号产率为零，只拟合背景
   RooConstVar zero_SSS("zero_SSS", "Constrained zero for SSS", 0.0);
   RooAbsReal* yield_SSS_orig = (RooAbsReal*)yield_SSS.clone("yield_SSS_orig");
   yield_SSS.setAttribute("Constant");
   yield_SSS.setVal(0);

   // 3. 用零信号假设再次拟合
   RooFitResult* fitResult_B = totalModel.fitTo(data, RooFit::Extended(true), 
   RooFit::Save(), RooFit::PrintLevel(1));
   double nll_B = fitResult_B->minNll(); // 仅背景的负对数似然值

   // 4. 恢复原始状态
   yield_SSS.setAttribute("Constant", false);
   yield_SSS.setVal(yield_SSS_orig->getVal());
   delete yield_SSS_orig;

   // 5. 计算似然比检验统计量
   double deltaLL = nll_B - nll_SB;
   double significance = sqrt(2 * deltaLL);

   std::cout << "\n===== Significance Calculation =====" << std::endl;
   std::cout << "Log-likelihood (S+B): " << -nll_SB << std::endl;
   std::cout << "Log-likelihood (B only): " << -nll_B << std::endl;
   std::cout << "Delta Log-likelihood: " << deltaLL << std::endl;
   std::cout << "Significance: " << significance << " sigma" << std::endl;



   // Save the output tree.
   TFile *outputFile = new TFile("filtered_data_secCut.root", "RECREATE");
   filteredTree->Write();
   outputFile->Close();
}
