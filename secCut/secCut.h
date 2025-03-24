//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 20 07:24:38 2025 by ROOT version 6.24/06
// from TTree filteredTree/filteredTree
// found on file: ../preCut/preCut_2023.root
//////////////////////////////////////////////////////////

#ifndef secCut_h
#define secCut_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class secCut {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *Jpsi_1_mass;
   vector<float>   *Jpsi_1_massErr;
   vector<float>   *Jpsi_1_massDiff;
   vector<float>   *Jpsi_1_ctau;
   vector<float>   *Jpsi_1_ctauErr;
   vector<float>   *Jpsi_1_Chi2;
   vector<float>   *Jpsi_1_ndof;
   vector<float>   *Jpsi_1_VtxProb;
   vector<float>   *Jpsi_1_px;
   vector<float>   *Jpsi_1_py;
   vector<float>   *Jpsi_1_pz;
   vector<float>   *Jpsi_1_phi;
   vector<float>   *Jpsi_1_eta;
   vector<float>   *Jpsi_1_pt;
   vector<float>   *Jpsi_1_mu_1_Idx;
   vector<float>   *Jpsi_1_mu_2_Idx;
   vector<float>   *Jpsi_2_mass;
   vector<float>   *Jpsi_2_massErr;
   vector<float>   *Jpsi_2_massDiff;
   vector<float>   *Jpsi_2_ctau;
   vector<float>   *Jpsi_2_ctauErr;
   vector<float>   *Jpsi_2_Chi2;
   vector<float>   *Jpsi_2_ndof;
   vector<float>   *Jpsi_2_VtxProb;
   vector<float>   *Jpsi_2_px;
   vector<float>   *Jpsi_2_py;
   vector<float>   *Jpsi_2_pz;
   vector<float>   *Jpsi_2_phi;
   vector<float>   *Jpsi_2_eta;
   vector<float>   *Jpsi_2_pt;
   vector<float>   *Jpsi_2_mu_1_Idx;
   vector<float>   *Jpsi_2_mu_2_Idx;
   vector<float>   *Pri_mass;
   vector<float>   *Pri_massErr;
   vector<float>   *Pri_ctau;
   vector<float>   *Pri_ctauErr;
   vector<float>   *Pri_Chi2;
   vector<float>   *Pri_ndof;
   vector<float>   *Pri_VtxProb;
   vector<float>   *Pri_px;
   vector<float>   *Pri_py;
   vector<float>   *Pri_pz;
   vector<float>   *Pri_phi;
   vector<float>   *Pri_eta;
   vector<float>   *Pri_pt;
   vector<float>   *Phi_mass;
   vector<float>   *Phi_massErr;
   vector<float>   *Phi_massDiff;
   vector<float>   *Phi_ctau;
   vector<float>   *Phi_ctauErr;
   vector<float>   *Phi_Chi2;
   vector<float>   *Phi_ndof;
   vector<float>   *Phi_VtxProb;
   vector<float>   *Phi_px;
   vector<float>   *Phi_py;
   vector<float>   *Phi_pz;
   vector<float>   *Phi_phi;
   vector<float>   *Phi_eta;
   vector<float>   *Phi_pt;
   vector<float>   *Phi_K_1_Idx;
   vector<float>   *Phi_K_2_Idx;
   vector<float>   *Jpsi_1_mu_1_px;
   vector<float>   *Jpsi_1_mu_1_py;
   vector<float>   *Jpsi_1_mu_1_pz;
   vector<float>   *Jpsi_1_mu_1_phi;
   vector<float>   *Jpsi_1_mu_1_eta;
   vector<float>   *Jpsi_1_mu_1_pt;
   vector<int>     *Jpsi_1_mu_1_isPatLooseMuon;
   vector<int>     *Jpsi_1_mu_1_isPatSoftMuon;
   vector<int>     *Jpsi_1_mu_1_isPatMediumMuon;
   vector<int>     *Jpsi_1_mu_1_isPatTightMuon;
   vector<float>   *Jpsi_1_mu_2_px;
   vector<float>   *Jpsi_1_mu_2_py;
   vector<float>   *Jpsi_1_mu_2_pz;
   vector<float>   *Jpsi_1_mu_2_phi;
   vector<float>   *Jpsi_1_mu_2_eta;
   vector<float>   *Jpsi_1_mu_2_pt;
   vector<int>     *Jpsi_1_mu_2_isPatLooseMuon;
   vector<int>     *Jpsi_1_mu_2_isPatSoftMuon;
   vector<int>     *Jpsi_1_mu_2_isPatMediumMuon;
   vector<int>     *Jpsi_1_mu_2_isPatTightMuon;
   vector<float>   *Jpsi_2_mu_1_px;
   vector<float>   *Jpsi_2_mu_1_py;
   vector<float>   *Jpsi_2_mu_1_pz;
   vector<float>   *Jpsi_2_mu_1_phi;
   vector<float>   *Jpsi_2_mu_1_eta;
   vector<float>   *Jpsi_2_mu_1_pt;
   vector<int>     *Jpsi_2_mu_1_isPatLooseMuon;
   vector<int>     *Jpsi_2_mu_1_isPatSoftMuon;
   vector<int>     *Jpsi_2_mu_1_isPatMediumMuon;
   vector<int>     *Jpsi_2_mu_1_isPatTightMuon;
   vector<float>   *Jpsi_2_mu_2_px;
   vector<float>   *Jpsi_2_mu_2_py;
   vector<float>   *Jpsi_2_mu_2_pz;
   vector<float>   *Jpsi_2_mu_2_phi;
   vector<float>   *Jpsi_2_mu_2_eta;
   vector<float>   *Jpsi_2_mu_2_pt;
   vector<int>     *Jpsi_2_mu_2_isPatLooseMuon;
   vector<int>     *Jpsi_2_mu_2_isPatSoftMuon;
   vector<int>     *Jpsi_2_mu_2_isPatMediumMuon;
   vector<int>     *Jpsi_2_mu_2_isPatTightMuon;
   vector<float>   *Phi_K_1_px;
   vector<float>   *Phi_K_1_py;
   vector<float>   *Phi_K_1_pz;
   vector<float>   *Phi_K_1_phi;
   vector<float>   *Phi_K_1_eta;
   vector<float>   *Phi_K_1_pt;
   vector<float>   *Phi_K_2_px;
   vector<float>   *Phi_K_2_py;
   vector<float>   *Phi_K_2_pz;
   vector<float>   *Phi_K_2_phi;
   vector<float>   *Phi_K_2_eta;
   vector<float>   *Phi_K_2_pt;

   vector<float>   *Jpsi_1_Lxy;
   vector<float>   *Jpsi_2_Lxy;
   vector<float>   *Phi_Lxy;

   // List of branches
   TBranch        *b_Jpsi_1_mass;   //!
   TBranch        *b_Jpsi_1_massErr;   //!
   TBranch        *b_Jpsi_1_massDiff;   //!
   TBranch        *b_Jpsi_1_ctau;   //!
   TBranch        *b_Jpsi_1_ctauErr;   //!
   TBranch        *b_Jpsi_1_Chi2;   //!
   TBranch        *b_Jpsi_1_ndof;   //!
   TBranch        *b_Jpsi_1_VtxProb;   //!
   TBranch        *b_Jpsi_1_px;   //!
   TBranch        *b_Jpsi_1_py;   //!
   TBranch        *b_Jpsi_1_pz;   //!
   TBranch        *b_Jpsi_1_phi;   //!
   TBranch        *b_Jpsi_1_eta;   //!
   TBranch        *b_Jpsi_1_pt;   //!
   TBranch        *b_Jpsi_1_mu_1_Idx;   //!
   TBranch        *b_Jpsi_1_mu_2_Idx;   //!
   TBranch        *b_Jpsi_2_mass;   //!
   TBranch        *b_Jpsi_2_massErr;   //!
   TBranch        *b_Jpsi_2_massDiff;   //!
   TBranch        *b_Jpsi_2_ctau;   //!
   TBranch        *b_Jpsi_2_ctauErr;   //!
   TBranch        *b_Jpsi_2_Chi2;   //!
   TBranch        *b_Jpsi_2_ndof;   //!
   TBranch        *b_Jpsi_2_VtxProb;   //!
   TBranch        *b_Jpsi_2_px;   //!
   TBranch        *b_Jpsi_2_py;   //!
   TBranch        *b_Jpsi_2_pz;   //!
   TBranch        *b_Jpsi_2_phi;   //!
   TBranch        *b_Jpsi_2_eta;   //!
   TBranch        *b_Jpsi_2_pt;   //!
   TBranch        *b_Jpsi_2_mu_1_Idx;   //!
   TBranch        *b_Jpsi_2_mu_2_Idx;   //!
   TBranch        *b_Pri_mass;   //!
   TBranch        *b_Pri_massErr;   //!
   TBranch        *b_Pri_ctau;   //!
   TBranch        *b_Pri_ctauErr;   //!
   TBranch        *b_Pri_Chi2;   //!
   TBranch        *b_Pri_ndof;   //!
   TBranch        *b_Pri_VtxProb;   //!
   TBranch        *b_Pri_px;   //!
   TBranch        *b_Pri_py;   //!
   TBranch        *b_Pri_pz;   //!
   TBranch        *b_Pri_phi;   //!
   TBranch        *b_Pri_eta;   //!
   TBranch        *b_Pri_pt;   //!
   TBranch        *b_Phi_mass;   //!
   TBranch        *b_Phi_massErr;   //!
   TBranch        *b_Phi_massDiff;   //!
   TBranch        *b_Phi_ctau;   //!
   TBranch        *b_Phi_ctauErr;   //!
   TBranch        *b_Phi_Chi2;   //!
   TBranch        *b_Phi_ndof;   //!
   TBranch        *b_Phi_VtxProb;   //!
   TBranch        *b_Phi_px;   //!
   TBranch        *b_Phi_py;   //!
   TBranch        *b_Phi_pz;   //!
   TBranch        *b_Phi_phi;   //!
   TBranch        *b_Phi_eta;   //!
   TBranch        *b_Phi_pt;   //!
   TBranch        *b_Phi_K_1_Idx;   //!
   TBranch        *b_Phi_K_2_Idx;   //!
   TBranch        *b_Jpsi_1_mu_1_px;   //!
   TBranch        *b_Jpsi_1_mu_1_py;   //!
   TBranch        *b_Jpsi_1_mu_1_pz;   //!
   TBranch        *b_Jpsi_1_mu_1_phi;   //!
   TBranch        *b_Jpsi_1_mu_1_eta;   //!
   TBranch        *b_Jpsi_1_mu_1_pt;   //!
   TBranch        *b_Jpsi_1_mu_1_isPatLooseMuon;   //!
   TBranch        *b_Jpsi_1_mu_1_isPatSoftMuon;   //!
   TBranch        *b_Jpsi_1_mu_1_isPatMediumMuon;   //!
   TBranch        *b_Jpsi_1_mu_1_isPatTightMuon;   //!
   TBranch        *b_Jpsi_1_mu_2_px;   //!
   TBranch        *b_Jpsi_1_mu_2_py;   //!
   TBranch        *b_Jpsi_1_mu_2_pz;   //!
   TBranch        *b_Jpsi_1_mu_2_phi;   //!
   TBranch        *b_Jpsi_1_mu_2_eta;   //!
   TBranch        *b_Jpsi_1_mu_2_pt;   //!
   TBranch        *b_Jpsi_1_mu_2_isPatLooseMuon;   //!
   TBranch        *b_Jpsi_1_mu_2_isPatSoftMuon;   //!
   TBranch        *b_Jpsi_1_mu_2_isPatMediumMuon;   //!
   TBranch        *b_Jpsi_1_mu_2_isPatTightMuon;   //!
   TBranch        *b_Jpsi_2_mu_1_px;   //!
   TBranch        *b_Jpsi_2_mu_1_py;   //!
   TBranch        *b_Jpsi_2_mu_1_pz;   //!
   TBranch        *b_Jpsi_2_mu_1_phi;   //!
   TBranch        *b_Jpsi_2_mu_1_eta;   //!
   TBranch        *b_Jpsi_2_mu_1_pt;   //!
   TBranch        *b_Jpsi_2_mu_1_isPatLooseMuon;   //!
   TBranch        *b_Jpsi_2_mu_1_isPatSoftMuon;   //!
   TBranch        *b_Jpsi_2_mu_1_isPatMediumMuon;   //!
   TBranch        *b_Jpsi_2_mu_1_isPatTightMuon;   //!
   TBranch        *b_Jpsi_2_mu_2_px;   //!
   TBranch        *b_Jpsi_2_mu_2_py;   //!
   TBranch        *b_Jpsi_2_mu_2_pz;   //!
   TBranch        *b_Jpsi_2_mu_2_phi;   //!
   TBranch        *b_Jpsi_2_mu_2_eta;   //!
   TBranch        *b_Jpsi_2_mu_2_pt;   //!
   TBranch        *b_Jpsi_2_mu_2_isPatLooseMuon;   //!
   TBranch        *b_Jpsi_2_mu_2_isPatSoftMuon;   //!
   TBranch        *b_Jpsi_2_mu_2_isPatMediumMuon;   //!
   TBranch        *b_Jpsi_2_mu_2_isPatTightMuon;   //!
   TBranch        *b_Phi_K_1_px;   //!
   TBranch        *b_Phi_K_1_py;   //!
   TBranch        *b_Phi_K_1_pz;   //!
   TBranch        *b_Phi_K_1_phi;   //!
   TBranch        *b_Phi_K_1_eta;   //!
   TBranch        *b_Phi_K_1_pt;   //!
   TBranch        *b_Phi_K_2_px;   //!
   TBranch        *b_Phi_K_2_py;   //!
   TBranch        *b_Phi_K_2_pz;   //!
   TBranch        *b_Phi_K_2_phi;   //!
   TBranch        *b_Phi_K_2_eta;   //!
   TBranch        *b_Phi_K_2_pt;   //!

      // Output tree
      TTree *filteredTree;

      // Define the branches for the filtered data.
      std::vector<float> *filtered_Jpsi_1_mass;
      std::vector<float> *filtered_Jpsi_1_massErr;
      std::vector<float> *filtered_Jpsi_1_massDiff;
      std::vector<float> *filtered_Jpsi_1_ctau;
      std::vector<float> *filtered_Jpsi_1_ctauErr;
      std::vector<float> *filtered_Jpsi_1_Chi2;
      std::vector<float> *filtered_Jpsi_1_ndof;
      std::vector<float> *filtered_Jpsi_1_VtxProb;
      std::vector<float> *filtered_Jpsi_1_px;
      std::vector<float> *filtered_Jpsi_1_py;
      std::vector<float> *filtered_Jpsi_1_pz;
      std::vector<float> *filtered_Jpsi_1_phi;
      std::vector<float> *filtered_Jpsi_1_eta;
      std::vector<float> *filtered_Jpsi_1_pt;
      std::vector<float> *filtered_Jpsi_1_mu_1_Idx;
      std::vector<float> *filtered_Jpsi_1_mu_2_Idx;
   
      std::vector<float> *filtered_Jpsi_2_mass;
      std::vector<float> *filtered_Jpsi_2_massErr;
      std::vector<float> *filtered_Jpsi_2_massDiff;
      std::vector<float> *filtered_Jpsi_2_ctau;
      std::vector<float> *filtered_Jpsi_2_ctauErr;
      std::vector<float> *filtered_Jpsi_2_Chi2;
      std::vector<float> *filtered_Jpsi_2_ndof;
      std::vector<float> *filtered_Jpsi_2_VtxProb;
      std::vector<float> *filtered_Jpsi_2_px;
      std::vector<float> *filtered_Jpsi_2_py;
      std::vector<float> *filtered_Jpsi_2_pz;
      std::vector<float> *filtered_Jpsi_2_phi;
      std::vector<float> *filtered_Jpsi_2_eta;
      std::vector<float> *filtered_Jpsi_2_pt;
      std::vector<float> *filtered_Jpsi_2_mu_1_Idx;
      std::vector<float> *filtered_Jpsi_2_mu_2_Idx;
   
      std::vector<float> *filtered_Pri_mass;
      std::vector<float> *filtered_Pri_massErr;
      std::vector<float> *filtered_Pri_ctau;
      std::vector<float> *filtered_Pri_ctauErr;
      std::vector<float> *filtered_Pri_Chi2;
      std::vector<float> *filtered_Pri_ndof;
      std::vector<float> *filtered_Pri_VtxProb;
      std::vector<float> *filtered_Pri_px;
      std::vector<float> *filtered_Pri_py;
      std::vector<float> *filtered_Pri_pz;
      std::vector<float> *filtered_Pri_phi;
      std::vector<float> *filtered_Pri_eta;
      std::vector<float> *filtered_Pri_pt;
   
      std::vector<float> *filtered_Phi_mass;
      std::vector<float> *filtered_Phi_massErr;
      std::vector<float> *filtered_Phi_massDiff;
      std::vector<float> *filtered_Phi_Chi2;
      std::vector<float> *filtered_Phi_ndof;
      std::vector<float> *filtered_Phi_VtxProb;
      std::vector<float> *filtered_Phi_px;
      std::vector<float> *filtered_Phi_py;
      std::vector<float> *filtered_Phi_pz;
      std::vector<float> *filtered_Phi_phi;
      std::vector<float> *filtered_Phi_eta;
      std::vector<float> *filtered_Phi_pt;
      std::vector<float> *filtered_Phi_ctau;
      std::vector<float> *filtered_Phi_ctauErr;
      std::vector<float> *filtered_Phi_K_1_Idx;
      std::vector<float> *filtered_Phi_K_2_Idx;
   
      // Transferring the muon data to the output tree
      std::vector<float> *filtered_Jpsi_1_mu_1_px;
      std::vector<float> *filtered_Jpsi_1_mu_1_py;
      std::vector<float> *filtered_Jpsi_1_mu_1_pz;
      std::vector<float> *filtered_Jpsi_1_mu_1_phi;
      std::vector<float> *filtered_Jpsi_1_mu_1_eta;
      std::vector<float> *filtered_Jpsi_1_mu_1_pt;
      std::vector<int>   *filtered_Jpsi_1_mu_1_isPatLooseMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_1_isPatSoftMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_1_isPatMediumMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_1_isPatTightMuon;
   
      std::vector<float> *filtered_Jpsi_1_mu_2_px;
      std::vector<float> *filtered_Jpsi_1_mu_2_py;
      std::vector<float> *filtered_Jpsi_1_mu_2_pz;
      std::vector<float> *filtered_Jpsi_1_mu_2_phi;
      std::vector<float> *filtered_Jpsi_1_mu_2_eta;
      std::vector<float> *filtered_Jpsi_1_mu_2_pt;
      std::vector<int>   *filtered_Jpsi_1_mu_2_isPatLooseMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_2_isPatSoftMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_2_isPatMediumMuon;
      std::vector<int>   *filtered_Jpsi_1_mu_2_isPatTightMuon;
   
      std::vector<float> *filtered_Jpsi_2_mu_1_px;
      std::vector<float> *filtered_Jpsi_2_mu_1_py;
      std::vector<float> *filtered_Jpsi_2_mu_1_pz;
      std::vector<float> *filtered_Jpsi_2_mu_1_phi;
      std::vector<float> *filtered_Jpsi_2_mu_1_eta;
      std::vector<float> *filtered_Jpsi_2_mu_1_pt;
      std::vector<int>   *filtered_Jpsi_2_mu_1_isPatLooseMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_1_isPatSoftMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_1_isPatMediumMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_1_isPatTightMuon;
   
      std::vector<float> *filtered_Jpsi_2_mu_2_px;
      std::vector<float> *filtered_Jpsi_2_mu_2_py;
      std::vector<float> *filtered_Jpsi_2_mu_2_pz;
      std::vector<float> *filtered_Jpsi_2_mu_2_phi;
      std::vector<float> *filtered_Jpsi_2_mu_2_eta;
      std::vector<float> *filtered_Jpsi_2_mu_2_pt;
      std::vector<int>   *filtered_Jpsi_2_mu_2_isPatLooseMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_2_isPatSoftMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_2_isPatMediumMuon;
      std::vector<int>   *filtered_Jpsi_2_mu_2_isPatTightMuon;
   
      std::vector<float> *filtered_Phi_K_1_px;
      std::vector<float> *filtered_Phi_K_1_py;
      std::vector<float> *filtered_Phi_K_1_pz;
      std::vector<float> *filtered_Phi_K_1_phi;
      std::vector<float> *filtered_Phi_K_1_eta;
      std::vector<float> *filtered_Phi_K_1_pt;
   
      std::vector<float> *filtered_Phi_K_2_px;
      std::vector<float> *filtered_Phi_K_2_py;
      std::vector<float> *filtered_Phi_K_2_pz;
      std::vector<float> *filtered_Phi_K_2_phi;
      std::vector<float> *filtered_Phi_K_2_eta;
      std::vector<float> *filtered_Phi_K_2_pt;

   secCut(TTree *tree=0);
   virtual ~secCut();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

      // Initialize the filtered data tree
      virtual void     InitFilteredTree();
      virtual void     ClearBranches();
};

#endif

#ifdef secCut_cxx
secCut::secCut(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../preCut/preCut_2023.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../preCut/preCut_2023.root");
      }
      f->GetObject("filteredTree",tree);

   }
   Init(tree);
   InitFilteredTree();
}

secCut::~secCut()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t secCut::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t secCut::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void secCut::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jpsi_1_mass = 0;
   Jpsi_1_massErr = 0;
   Jpsi_1_massDiff = 0;
   Jpsi_1_ctau = 0;
   Jpsi_1_ctauErr = 0;
   Jpsi_1_Chi2 = 0;
   Jpsi_1_ndof = 0;
   Jpsi_1_VtxProb = 0;
   Jpsi_1_px = 0;
   Jpsi_1_py = 0;
   Jpsi_1_pz = 0;
   Jpsi_1_phi = 0;
   Jpsi_1_eta = 0;
   Jpsi_1_pt = 0;
   Jpsi_1_mu_1_Idx = 0;
   Jpsi_1_mu_2_Idx = 0;
   Jpsi_2_mass = 0;
   Jpsi_2_massErr = 0;
   Jpsi_2_massDiff = 0;
   Jpsi_2_ctau = 0;
   Jpsi_2_ctauErr = 0;
   Jpsi_2_Chi2 = 0;
   Jpsi_2_ndof = 0;
   Jpsi_2_VtxProb = 0;
   Jpsi_2_px = 0;
   Jpsi_2_py = 0;
   Jpsi_2_pz = 0;
   Jpsi_2_phi = 0;
   Jpsi_2_eta = 0;
   Jpsi_2_pt = 0;
   Jpsi_2_mu_1_Idx = 0;
   Jpsi_2_mu_2_Idx = 0;
   Pri_mass = 0;
   Pri_massErr = 0;
   Pri_ctau = 0;
   Pri_ctauErr = 0;
   Pri_Chi2 = 0;
   Pri_ndof = 0;
   Pri_VtxProb = 0;
   Pri_px = 0;
   Pri_py = 0;
   Pri_pz = 0;
   Pri_phi = 0;
   Pri_eta = 0;
   Pri_pt = 0;
   Phi_mass = 0;
   Phi_massErr = 0;
   Phi_massDiff = 0;
   Phi_ctau = 0;
   Phi_ctauErr = 0;
   Phi_Chi2 = 0;
   Phi_ndof = 0;
   Phi_VtxProb = 0;
   Phi_px = 0;
   Phi_py = 0;
   Phi_pz = 0;
   Phi_phi = 0;
   Phi_eta = 0;
   Phi_pt = 0;
   Phi_K_1_Idx = 0;
   Phi_K_2_Idx = 0;
   Jpsi_1_mu_1_px = 0;
   Jpsi_1_mu_1_py = 0;
   Jpsi_1_mu_1_pz = 0;
   Jpsi_1_mu_1_phi = 0;
   Jpsi_1_mu_1_eta = 0;
   Jpsi_1_mu_1_pt = 0;
   Jpsi_1_mu_1_isPatLooseMuon = 0;
   Jpsi_1_mu_1_isPatSoftMuon = 0;
   Jpsi_1_mu_1_isPatMediumMuon = 0;
   Jpsi_1_mu_1_isPatTightMuon = 0;
   Jpsi_1_mu_2_px = 0;
   Jpsi_1_mu_2_py = 0;
   Jpsi_1_mu_2_pz = 0;
   Jpsi_1_mu_2_phi = 0;
   Jpsi_1_mu_2_eta = 0;
   Jpsi_1_mu_2_pt = 0;
   Jpsi_1_mu_2_isPatLooseMuon = 0;
   Jpsi_1_mu_2_isPatSoftMuon = 0;
   Jpsi_1_mu_2_isPatMediumMuon = 0;
   Jpsi_1_mu_2_isPatTightMuon = 0;
   Jpsi_2_mu_1_px = 0;
   Jpsi_2_mu_1_py = 0;
   Jpsi_2_mu_1_pz = 0;
   Jpsi_2_mu_1_phi = 0;
   Jpsi_2_mu_1_eta = 0;
   Jpsi_2_mu_1_pt = 0;
   Jpsi_2_mu_1_isPatLooseMuon = 0;
   Jpsi_2_mu_1_isPatSoftMuon = 0;
   Jpsi_2_mu_1_isPatMediumMuon = 0;
   Jpsi_2_mu_1_isPatTightMuon = 0;
   Jpsi_2_mu_2_px = 0;
   Jpsi_2_mu_2_py = 0;
   Jpsi_2_mu_2_pz = 0;
   Jpsi_2_mu_2_phi = 0;
   Jpsi_2_mu_2_eta = 0;
   Jpsi_2_mu_2_pt = 0;
   Jpsi_2_mu_2_isPatLooseMuon = 0;
   Jpsi_2_mu_2_isPatSoftMuon = 0;
   Jpsi_2_mu_2_isPatMediumMuon = 0;
   Jpsi_2_mu_2_isPatTightMuon = 0;
   Phi_K_1_px = 0;
   Phi_K_1_py = 0;
   Phi_K_1_pz = 0;
   Phi_K_1_phi = 0;
   Phi_K_1_eta = 0;
   Phi_K_1_pt = 0;
   Phi_K_2_px = 0;
   Phi_K_2_py = 0;
   Phi_K_2_pz = 0;
   Phi_K_2_phi = 0;
   Phi_K_2_eta = 0;
   Phi_K_2_pt = 0;

   Jpsi_1_Lxy = 0;
   Jpsi_2_Lxy = 0;
   Phi_Lxy = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Jpsi_1_mass", &Jpsi_1_mass, &b_Jpsi_1_mass);
   fChain->SetBranchAddress("Jpsi_1_massErr", &Jpsi_1_massErr, &b_Jpsi_1_massErr);
   fChain->SetBranchAddress("Jpsi_1_massDiff", &Jpsi_1_massDiff, &b_Jpsi_1_massDiff);
   fChain->SetBranchAddress("Jpsi_1_ctau", &Jpsi_1_ctau, &b_Jpsi_1_ctau);
   fChain->SetBranchAddress("Jpsi_1_ctauErr", &Jpsi_1_ctauErr, &b_Jpsi_1_ctauErr);
   fChain->SetBranchAddress("Jpsi_1_Chi2", &Jpsi_1_Chi2, &b_Jpsi_1_Chi2);
   fChain->SetBranchAddress("Jpsi_1_ndof", &Jpsi_1_ndof, &b_Jpsi_1_ndof);
   fChain->SetBranchAddress("Jpsi_1_VtxProb", &Jpsi_1_VtxProb, &b_Jpsi_1_VtxProb);
   fChain->SetBranchAddress("Jpsi_1_px", &Jpsi_1_px, &b_Jpsi_1_px);
   fChain->SetBranchAddress("Jpsi_1_py", &Jpsi_1_py, &b_Jpsi_1_py);
   fChain->SetBranchAddress("Jpsi_1_pz", &Jpsi_1_pz, &b_Jpsi_1_pz);
   fChain->SetBranchAddress("Jpsi_1_phi", &Jpsi_1_phi, &b_Jpsi_1_phi);
   fChain->SetBranchAddress("Jpsi_1_eta", &Jpsi_1_eta, &b_Jpsi_1_eta);
   fChain->SetBranchAddress("Jpsi_1_pt", &Jpsi_1_pt, &b_Jpsi_1_pt);
   fChain->SetBranchAddress("Jpsi_1_mu_1_Idx", &Jpsi_1_mu_1_Idx, &b_Jpsi_1_mu_1_Idx);
   fChain->SetBranchAddress("Jpsi_1_mu_2_Idx", &Jpsi_1_mu_2_Idx, &b_Jpsi_1_mu_2_Idx);
   fChain->SetBranchAddress("Jpsi_2_mass", &Jpsi_2_mass, &b_Jpsi_2_mass);
   fChain->SetBranchAddress("Jpsi_2_massErr", &Jpsi_2_massErr, &b_Jpsi_2_massErr);
   fChain->SetBranchAddress("Jpsi_2_massDiff", &Jpsi_2_massDiff, &b_Jpsi_2_massDiff);
   fChain->SetBranchAddress("Jpsi_2_ctau", &Jpsi_2_ctau, &b_Jpsi_2_ctau);
   fChain->SetBranchAddress("Jpsi_2_ctauErr", &Jpsi_2_ctauErr, &b_Jpsi_2_ctauErr);
   fChain->SetBranchAddress("Jpsi_2_Chi2", &Jpsi_2_Chi2, &b_Jpsi_2_Chi2);
   fChain->SetBranchAddress("Jpsi_2_ndof", &Jpsi_2_ndof, &b_Jpsi_2_ndof);
   fChain->SetBranchAddress("Jpsi_2_VtxProb", &Jpsi_2_VtxProb, &b_Jpsi_2_VtxProb);
   fChain->SetBranchAddress("Jpsi_2_px", &Jpsi_2_px, &b_Jpsi_2_px);
   fChain->SetBranchAddress("Jpsi_2_py", &Jpsi_2_py, &b_Jpsi_2_py);
   fChain->SetBranchAddress("Jpsi_2_pz", &Jpsi_2_pz, &b_Jpsi_2_pz);
   fChain->SetBranchAddress("Jpsi_2_phi", &Jpsi_2_phi, &b_Jpsi_2_phi);
   fChain->SetBranchAddress("Jpsi_2_eta", &Jpsi_2_eta, &b_Jpsi_2_eta);
   fChain->SetBranchAddress("Jpsi_2_pt", &Jpsi_2_pt, &b_Jpsi_2_pt);
   fChain->SetBranchAddress("Jpsi_2_mu_1_Idx", &Jpsi_2_mu_1_Idx, &b_Jpsi_2_mu_1_Idx);
   fChain->SetBranchAddress("Jpsi_2_mu_2_Idx", &Jpsi_2_mu_2_Idx, &b_Jpsi_2_mu_2_Idx);
   fChain->SetBranchAddress("Pri_mass", &Pri_mass, &b_Pri_mass);
   fChain->SetBranchAddress("Pri_massErr", &Pri_massErr, &b_Pri_massErr);
   fChain->SetBranchAddress("Pri_ctau", &Pri_ctau, &b_Pri_ctau);
   fChain->SetBranchAddress("Pri_ctauErr", &Pri_ctauErr, &b_Pri_ctauErr);
   fChain->SetBranchAddress("Pri_Chi2", &Pri_Chi2, &b_Pri_Chi2);
   fChain->SetBranchAddress("Pri_ndof", &Pri_ndof, &b_Pri_ndof);
   fChain->SetBranchAddress("Pri_VtxProb", &Pri_VtxProb, &b_Pri_VtxProb);
   fChain->SetBranchAddress("Pri_px", &Pri_px, &b_Pri_px);
   fChain->SetBranchAddress("Pri_py", &Pri_py, &b_Pri_py);
   fChain->SetBranchAddress("Pri_pz", &Pri_pz, &b_Pri_pz);
   fChain->SetBranchAddress("Pri_phi", &Pri_phi, &b_Pri_phi);
   fChain->SetBranchAddress("Pri_eta", &Pri_eta, &b_Pri_eta);
   fChain->SetBranchAddress("Pri_pt", &Pri_pt, &b_Pri_pt);
   fChain->SetBranchAddress("Phi_mass", &Phi_mass, &b_Phi_mass);
   fChain->SetBranchAddress("Phi_massErr", &Phi_massErr, &b_Phi_massErr);
   fChain->SetBranchAddress("Phi_massDiff", &Phi_massDiff, &b_Phi_massDiff);
   fChain->SetBranchAddress("Phi_ctau", &Phi_ctau, &b_Phi_ctau);
   fChain->SetBranchAddress("Phi_ctauErr", &Phi_ctauErr, &b_Phi_ctauErr);
   fChain->SetBranchAddress("Phi_Chi2", &Phi_Chi2, &b_Phi_Chi2);
   fChain->SetBranchAddress("Phi_ndof", &Phi_ndof, &b_Phi_ndof);
   fChain->SetBranchAddress("Phi_VtxProb", &Phi_VtxProb, &b_Phi_VtxProb);
   fChain->SetBranchAddress("Phi_px", &Phi_px, &b_Phi_px);
   fChain->SetBranchAddress("Phi_py", &Phi_py, &b_Phi_py);
   fChain->SetBranchAddress("Phi_pz", &Phi_pz, &b_Phi_pz);
   fChain->SetBranchAddress("Phi_phi", &Phi_phi, &b_Phi_phi);
   fChain->SetBranchAddress("Phi_eta", &Phi_eta, &b_Phi_eta);
   fChain->SetBranchAddress("Phi_pt", &Phi_pt, &b_Phi_pt);
   fChain->SetBranchAddress("Phi_K_1_Idx", &Phi_K_1_Idx, &b_Phi_K_1_Idx);
   fChain->SetBranchAddress("Phi_K_2_Idx", &Phi_K_2_Idx, &b_Phi_K_2_Idx);
   fChain->SetBranchAddress("Jpsi_1_mu_1_px", &Jpsi_1_mu_1_px, &b_Jpsi_1_mu_1_px);
   fChain->SetBranchAddress("Jpsi_1_mu_1_py", &Jpsi_1_mu_1_py, &b_Jpsi_1_mu_1_py);
   fChain->SetBranchAddress("Jpsi_1_mu_1_pz", &Jpsi_1_mu_1_pz, &b_Jpsi_1_mu_1_pz);
   fChain->SetBranchAddress("Jpsi_1_mu_1_phi", &Jpsi_1_mu_1_phi, &b_Jpsi_1_mu_1_phi);
   fChain->SetBranchAddress("Jpsi_1_mu_1_eta", &Jpsi_1_mu_1_eta, &b_Jpsi_1_mu_1_eta);
   fChain->SetBranchAddress("Jpsi_1_mu_1_pt", &Jpsi_1_mu_1_pt, &b_Jpsi_1_mu_1_pt);
   fChain->SetBranchAddress("Jpsi_1_mu_1_isPatLooseMuon", &Jpsi_1_mu_1_isPatLooseMuon, &b_Jpsi_1_mu_1_isPatLooseMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_1_isPatSoftMuon", &Jpsi_1_mu_1_isPatSoftMuon, &b_Jpsi_1_mu_1_isPatSoftMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_1_isPatMediumMuon", &Jpsi_1_mu_1_isPatMediumMuon, &b_Jpsi_1_mu_1_isPatMediumMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_1_isPatTightMuon", &Jpsi_1_mu_1_isPatTightMuon, &b_Jpsi_1_mu_1_isPatTightMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_2_px", &Jpsi_1_mu_2_px, &b_Jpsi_1_mu_2_px);
   fChain->SetBranchAddress("Jpsi_1_mu_2_py", &Jpsi_1_mu_2_py, &b_Jpsi_1_mu_2_py);
   fChain->SetBranchAddress("Jpsi_1_mu_2_pz", &Jpsi_1_mu_2_pz, &b_Jpsi_1_mu_2_pz);
   fChain->SetBranchAddress("Jpsi_1_mu_2_phi", &Jpsi_1_mu_2_phi, &b_Jpsi_1_mu_2_phi);
   fChain->SetBranchAddress("Jpsi_1_mu_2_eta", &Jpsi_1_mu_2_eta, &b_Jpsi_1_mu_2_eta);
   fChain->SetBranchAddress("Jpsi_1_mu_2_pt", &Jpsi_1_mu_2_pt, &b_Jpsi_1_mu_2_pt);
   fChain->SetBranchAddress("Jpsi_1_mu_2_isPatLooseMuon", &Jpsi_1_mu_2_isPatLooseMuon, &b_Jpsi_1_mu_2_isPatLooseMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_2_isPatSoftMuon", &Jpsi_1_mu_2_isPatSoftMuon, &b_Jpsi_1_mu_2_isPatSoftMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_2_isPatMediumMuon", &Jpsi_1_mu_2_isPatMediumMuon, &b_Jpsi_1_mu_2_isPatMediumMuon);
   fChain->SetBranchAddress("Jpsi_1_mu_2_isPatTightMuon", &Jpsi_1_mu_2_isPatTightMuon, &b_Jpsi_1_mu_2_isPatTightMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_1_px", &Jpsi_2_mu_1_px, &b_Jpsi_2_mu_1_px);
   fChain->SetBranchAddress("Jpsi_2_mu_1_py", &Jpsi_2_mu_1_py, &b_Jpsi_2_mu_1_py);
   fChain->SetBranchAddress("Jpsi_2_mu_1_pz", &Jpsi_2_mu_1_pz, &b_Jpsi_2_mu_1_pz);
   fChain->SetBranchAddress("Jpsi_2_mu_1_phi", &Jpsi_2_mu_1_phi, &b_Jpsi_2_mu_1_phi);
   fChain->SetBranchAddress("Jpsi_2_mu_1_eta", &Jpsi_2_mu_1_eta, &b_Jpsi_2_mu_1_eta);
   fChain->SetBranchAddress("Jpsi_2_mu_1_pt", &Jpsi_2_mu_1_pt, &b_Jpsi_2_mu_1_pt);
   fChain->SetBranchAddress("Jpsi_2_mu_1_isPatLooseMuon", &Jpsi_2_mu_1_isPatLooseMuon, &b_Jpsi_2_mu_1_isPatLooseMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_1_isPatSoftMuon", &Jpsi_2_mu_1_isPatSoftMuon, &b_Jpsi_2_mu_1_isPatSoftMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_1_isPatMediumMuon", &Jpsi_2_mu_1_isPatMediumMuon, &b_Jpsi_2_mu_1_isPatMediumMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_1_isPatTightMuon", &Jpsi_2_mu_1_isPatTightMuon, &b_Jpsi_2_mu_1_isPatTightMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_2_px", &Jpsi_2_mu_2_px, &b_Jpsi_2_mu_2_px);
   fChain->SetBranchAddress("Jpsi_2_mu_2_py", &Jpsi_2_mu_2_py, &b_Jpsi_2_mu_2_py);
   fChain->SetBranchAddress("Jpsi_2_mu_2_pz", &Jpsi_2_mu_2_pz, &b_Jpsi_2_mu_2_pz);
   fChain->SetBranchAddress("Jpsi_2_mu_2_phi", &Jpsi_2_mu_2_phi, &b_Jpsi_2_mu_2_phi);
   fChain->SetBranchAddress("Jpsi_2_mu_2_eta", &Jpsi_2_mu_2_eta, &b_Jpsi_2_mu_2_eta);
   fChain->SetBranchAddress("Jpsi_2_mu_2_pt", &Jpsi_2_mu_2_pt, &b_Jpsi_2_mu_2_pt);
   fChain->SetBranchAddress("Jpsi_2_mu_2_isPatLooseMuon", &Jpsi_2_mu_2_isPatLooseMuon, &b_Jpsi_2_mu_2_isPatLooseMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_2_isPatSoftMuon", &Jpsi_2_mu_2_isPatSoftMuon, &b_Jpsi_2_mu_2_isPatSoftMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_2_isPatMediumMuon", &Jpsi_2_mu_2_isPatMediumMuon, &b_Jpsi_2_mu_2_isPatMediumMuon);
   fChain->SetBranchAddress("Jpsi_2_mu_2_isPatTightMuon", &Jpsi_2_mu_2_isPatTightMuon, &b_Jpsi_2_mu_2_isPatTightMuon);
   fChain->SetBranchAddress("Phi_K_1_px", &Phi_K_1_px, &b_Phi_K_1_px);
   fChain->SetBranchAddress("Phi_K_1_py", &Phi_K_1_py, &b_Phi_K_1_py);
   fChain->SetBranchAddress("Phi_K_1_pz", &Phi_K_1_pz, &b_Phi_K_1_pz);
   fChain->SetBranchAddress("Phi_K_1_phi", &Phi_K_1_phi, &b_Phi_K_1_phi);
   fChain->SetBranchAddress("Phi_K_1_eta", &Phi_K_1_eta, &b_Phi_K_1_eta);
   fChain->SetBranchAddress("Phi_K_1_pt", &Phi_K_1_pt, &b_Phi_K_1_pt);
   fChain->SetBranchAddress("Phi_K_2_px", &Phi_K_2_px, &b_Phi_K_2_px);
   fChain->SetBranchAddress("Phi_K_2_py", &Phi_K_2_py, &b_Phi_K_2_py);
   fChain->SetBranchAddress("Phi_K_2_pz", &Phi_K_2_pz, &b_Phi_K_2_pz);
   fChain->SetBranchAddress("Phi_K_2_phi", &Phi_K_2_phi, &b_Phi_K_2_phi);
   fChain->SetBranchAddress("Phi_K_2_eta", &Phi_K_2_eta, &b_Phi_K_2_eta);
   fChain->SetBranchAddress("Phi_K_2_pt", &Phi_K_2_pt, &b_Phi_K_2_pt);
   Notify();

   // Initialize the output tree branches
   filtered_Jpsi_1_mass = 0;
   filtered_Jpsi_1_massErr = 0;
   filtered_Jpsi_1_massDiff = 0;
   filtered_Jpsi_1_ctau = 0;
   filtered_Jpsi_1_ctauErr = 0;
   filtered_Jpsi_1_Chi2 = 0;
   filtered_Jpsi_1_ndof = 0;
   filtered_Jpsi_1_VtxProb = 0;
   filtered_Jpsi_1_px = 0;
   filtered_Jpsi_1_py = 0;
   filtered_Jpsi_1_pz = 0;
   filtered_Jpsi_1_phi = 0;
   filtered_Jpsi_1_eta = 0;
   filtered_Jpsi_1_pt = 0;
   filtered_Jpsi_1_mu_1_Idx = 0;
   filtered_Jpsi_1_mu_2_Idx = 0;

   filtered_Jpsi_2_mass = 0;
   filtered_Jpsi_2_massErr = 0;
   filtered_Jpsi_2_massDiff = 0;
   filtered_Jpsi_2_ctau = 0;
   filtered_Jpsi_2_ctauErr = 0;
   filtered_Jpsi_2_Chi2 = 0;
   filtered_Jpsi_2_ndof = 0;
   filtered_Jpsi_2_VtxProb = 0;
   filtered_Jpsi_2_px = 0;
   filtered_Jpsi_2_py = 0;
   filtered_Jpsi_2_pz = 0;
   filtered_Jpsi_2_phi = 0;
   filtered_Jpsi_2_eta = 0;
   filtered_Jpsi_2_pt = 0;
   filtered_Jpsi_2_mu_1_Idx = 0;
   filtered_Jpsi_2_mu_2_Idx = 0;

   filtered_Pri_mass = 0;
   filtered_Pri_massErr = 0;
   filtered_Pri_ctau = 0;
   filtered_Pri_ctauErr = 0;
   filtered_Pri_Chi2 = 0;
   filtered_Pri_ndof = 0;
   filtered_Pri_VtxProb = 0;
   filtered_Pri_px = 0;
   filtered_Pri_py = 0;
   filtered_Pri_pz = 0;
   filtered_Pri_phi = 0;
   filtered_Pri_eta = 0;
   filtered_Pri_pt = 0;

   filtered_Phi_mass = 0;
   filtered_Phi_massErr = 0;
   filtered_Phi_massDiff = 0;
   filtered_Phi_Chi2 = 0;
   filtered_Phi_ndof = 0;
   filtered_Phi_VtxProb = 0;
   filtered_Phi_px = 0;
   filtered_Phi_py = 0;
   filtered_Phi_pz = 0;
   filtered_Phi_phi = 0;
   filtered_Phi_eta = 0;
   filtered_Phi_pt = 0;
   filtered_Phi_ctau = 0;
   filtered_Phi_ctauErr = 0;
   filtered_Phi_K_1_Idx = 0;
   filtered_Phi_K_2_Idx = 0;

      // Transferring the muon data to the output tree
   filtered_Jpsi_1_mu_1_px = 0;
   filtered_Jpsi_1_mu_1_py = 0;
   filtered_Jpsi_1_mu_1_pz = 0;
   filtered_Jpsi_1_mu_1_phi = 0;
   filtered_Jpsi_1_mu_1_eta = 0;
   filtered_Jpsi_1_mu_1_pt = 0;
   filtered_Jpsi_1_mu_1_isPatLooseMuon = 0;
   filtered_Jpsi_1_mu_1_isPatSoftMuon = 0;
   filtered_Jpsi_1_mu_1_isPatMediumMuon = 0;
   filtered_Jpsi_1_mu_1_isPatTightMuon = 0;

   filtered_Jpsi_1_mu_2_px = 0;
   filtered_Jpsi_1_mu_2_py = 0;
   filtered_Jpsi_1_mu_2_pz = 0;
   filtered_Jpsi_1_mu_2_phi = 0;
   filtered_Jpsi_1_mu_2_eta = 0;
   filtered_Jpsi_1_mu_2_pt = 0;
   filtered_Jpsi_1_mu_2_isPatLooseMuon = 0;
   filtered_Jpsi_1_mu_2_isPatSoftMuon = 0;
   filtered_Jpsi_1_mu_2_isPatMediumMuon = 0;
   filtered_Jpsi_1_mu_2_isPatTightMuon = 0;

   filtered_Jpsi_2_mu_1_px = 0;
   filtered_Jpsi_2_mu_1_py = 0;
   filtered_Jpsi_2_mu_1_pz = 0;
   filtered_Jpsi_2_mu_1_phi = 0;
   filtered_Jpsi_2_mu_1_eta = 0;
   filtered_Jpsi_2_mu_1_pt = 0;
   filtered_Jpsi_2_mu_1_isPatLooseMuon = 0;
   filtered_Jpsi_2_mu_1_isPatSoftMuon = 0;
   filtered_Jpsi_2_mu_1_isPatMediumMuon = 0;
   filtered_Jpsi_2_mu_1_isPatTightMuon = 0;

   filtered_Jpsi_2_mu_2_px = 0;
   filtered_Jpsi_2_mu_2_py = 0;
   filtered_Jpsi_2_mu_2_pz = 0;
   filtered_Jpsi_2_mu_2_phi = 0;
   filtered_Jpsi_2_mu_2_eta = 0;
   filtered_Jpsi_2_mu_2_pt = 0;
   filtered_Jpsi_2_mu_2_isPatLooseMuon = 0;
   filtered_Jpsi_2_mu_2_isPatSoftMuon = 0;
   filtered_Jpsi_2_mu_2_isPatMediumMuon = 0;
   filtered_Jpsi_2_mu_2_isPatTightMuon = 0;

   filtered_Phi_K_1_px = 0;
   filtered_Phi_K_1_py = 0;
   filtered_Phi_K_1_pz = 0;
   filtered_Phi_K_1_phi = 0;
   filtered_Phi_K_1_eta = 0;
   filtered_Phi_K_1_pt = 0;

   filtered_Phi_K_2_px = 0;
   filtered_Phi_K_2_py = 0;
   filtered_Phi_K_2_pz = 0;
   filtered_Phi_K_2_phi = 0;
   filtered_Phi_K_2_eta = 0;
   filtered_Phi_K_2_pt = 0;
}

Bool_t secCut::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void secCut::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t secCut::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void secCut::ClearBranches(){
   Jpsi_1_Lxy->clear();
   Jpsi_2_Lxy->clear();
   Phi_Lxy->clear();

   filtered_Jpsi_1_mass->clear();
   filtered_Jpsi_1_massErr->clear();
   filtered_Jpsi_1_massDiff->clear();
   filtered_Jpsi_1_ctau->clear();
   filtered_Jpsi_1_ctauErr->clear();
   filtered_Jpsi_1_Chi2->clear();
   filtered_Jpsi_1_ndof->clear();
   filtered_Jpsi_1_VtxProb->clear();
   filtered_Jpsi_1_px->clear();
   filtered_Jpsi_1_py->clear();
   filtered_Jpsi_1_pz->clear();
   filtered_Jpsi_1_phi->clear();
   filtered_Jpsi_1_eta->clear();
   filtered_Jpsi_1_pt->clear();
   filtered_Jpsi_1_mu_1_Idx->clear();
   filtered_Jpsi_1_mu_2_Idx->clear();

   filtered_Jpsi_2_mass->clear();
   filtered_Jpsi_2_massErr->clear();
   filtered_Jpsi_2_massDiff->clear();
   filtered_Jpsi_2_ctau->clear();
   filtered_Jpsi_2_ctauErr->clear();
   filtered_Jpsi_2_Chi2->clear();
   filtered_Jpsi_2_ndof->clear();
   filtered_Jpsi_2_VtxProb->clear();
   filtered_Jpsi_2_px->clear();
   filtered_Jpsi_2_py->clear();
   filtered_Jpsi_2_pz->clear();
   filtered_Jpsi_2_phi->clear();
   filtered_Jpsi_2_eta->clear();
   filtered_Jpsi_2_pt->clear();
   filtered_Jpsi_2_mu_1_Idx->clear();
   filtered_Jpsi_2_mu_2_Idx->clear();

   filtered_Pri_mass->clear();
   filtered_Pri_massErr->clear();
   filtered_Pri_ctau->clear();
   filtered_Pri_ctauErr->clear();
   filtered_Pri_Chi2->clear();
   filtered_Pri_ndof->clear();
   filtered_Pri_VtxProb->clear();
   filtered_Pri_px->clear();
   filtered_Pri_py->clear();
   filtered_Pri_pz->clear();
   filtered_Pri_phi->clear();
   filtered_Pri_eta->clear();
   filtered_Pri_pt->clear();

   filtered_Phi_mass->clear();
   filtered_Phi_massErr->clear();
   filtered_Phi_massDiff->clear();
   filtered_Phi_Chi2->clear();
   filtered_Phi_ndof->clear();
   filtered_Phi_VtxProb->clear();
   filtered_Phi_ctau->clear();
   filtered_Phi_ctauErr->clear();
   filtered_Phi_px->clear();
   filtered_Phi_py->clear();
   filtered_Phi_pz->clear();
   filtered_Phi_phi->clear();
   filtered_Phi_eta->clear();
   filtered_Phi_pt->clear();
   filtered_Phi_K_1_Idx->clear();
   filtered_Phi_K_2_Idx->clear();

      // Transferring the muon data to the output tree
   filtered_Jpsi_1_mu_1_px->clear();
   filtered_Jpsi_1_mu_1_py->clear();
   filtered_Jpsi_1_mu_1_pz->clear();
   filtered_Jpsi_1_mu_1_phi->clear();
   filtered_Jpsi_1_mu_1_eta->clear();
   filtered_Jpsi_1_mu_1_pt->clear();
   filtered_Jpsi_1_mu_1_isPatLooseMuon->clear();
   filtered_Jpsi_1_mu_1_isPatSoftMuon->clear();
   filtered_Jpsi_1_mu_1_isPatMediumMuon->clear();
   filtered_Jpsi_1_mu_1_isPatTightMuon->clear();

   filtered_Jpsi_1_mu_2_px->clear();
   filtered_Jpsi_1_mu_2_py->clear();
   filtered_Jpsi_1_mu_2_pz->clear();
   filtered_Jpsi_1_mu_2_phi->clear();
   filtered_Jpsi_1_mu_2_eta->clear();
   filtered_Jpsi_1_mu_2_pt->clear();
   filtered_Jpsi_1_mu_2_isPatLooseMuon->clear();
   filtered_Jpsi_1_mu_2_isPatSoftMuon->clear();
   filtered_Jpsi_1_mu_2_isPatMediumMuon->clear();
   filtered_Jpsi_1_mu_2_isPatTightMuon->clear();

   filtered_Jpsi_2_mu_1_px->clear();
   filtered_Jpsi_2_mu_1_py->clear();
   filtered_Jpsi_2_mu_1_pz->clear();
   filtered_Jpsi_2_mu_1_phi->clear();
   filtered_Jpsi_2_mu_1_eta->clear();
   filtered_Jpsi_2_mu_1_pt->clear();
   filtered_Jpsi_2_mu_1_isPatLooseMuon->clear();
   filtered_Jpsi_2_mu_1_isPatSoftMuon->clear();
   filtered_Jpsi_2_mu_1_isPatMediumMuon->clear();
   filtered_Jpsi_2_mu_1_isPatTightMuon->clear();

   filtered_Jpsi_2_mu_2_px->clear();
   filtered_Jpsi_2_mu_2_py->clear();
   filtered_Jpsi_2_mu_2_pz->clear();
   filtered_Jpsi_2_mu_2_phi->clear();
   filtered_Jpsi_2_mu_2_eta->clear();
   filtered_Jpsi_2_mu_2_pt->clear();
   filtered_Jpsi_2_mu_2_isPatLooseMuon->clear();
   filtered_Jpsi_2_mu_2_isPatSoftMuon->clear();
   filtered_Jpsi_2_mu_2_isPatMediumMuon->clear();
   filtered_Jpsi_2_mu_2_isPatTightMuon->clear();

   filtered_Phi_K_1_px->clear();
   filtered_Phi_K_1_py->clear();
   filtered_Phi_K_1_pz->clear();
   filtered_Phi_K_1_phi->clear();
   filtered_Phi_K_1_eta->clear();
   filtered_Phi_K_1_pt->clear();

   filtered_Phi_K_2_px->clear();
   filtered_Phi_K_2_py->clear();
   filtered_Phi_K_2_pz->clear();
   filtered_Phi_K_2_phi->clear();
   filtered_Phi_K_2_eta->clear();
   filtered_Phi_K_2_pt->clear();
}

void secCut::InitFilteredTree()
{
   filteredTree = new TTree("filteredTree", "filteredTree");
   // Define the branches for the filtered data.
   filteredTree->Branch("Jpsi_1_Lxy", &Jpsi_1_Lxy);
   filteredTree->Branch("Jpsi_2_Lxy", &Jpsi_2_Lxy);
   filteredTree->Branch("Phi_Lxy", &Phi_Lxy);

   filteredTree->Branch("Jpsi_1_mass", &filtered_Jpsi_1_mass);
   filteredTree->Branch("Jpsi_1_massErr", &filtered_Jpsi_1_massErr);
   filteredTree->Branch("Jpsi_1_massDiff", &filtered_Jpsi_1_massDiff);
   filteredTree->Branch("Jpsi_1_ctau", &filtered_Jpsi_1_ctau);
   filteredTree->Branch("Jpsi_1_ctauErr", &filtered_Jpsi_1_ctauErr);
   filteredTree->Branch("Jpsi_1_Chi2", &filtered_Jpsi_1_Chi2);
   filteredTree->Branch("Jpsi_1_ndof", &filtered_Jpsi_1_ndof);
   filteredTree->Branch("Jpsi_1_VtxProb", &filtered_Jpsi_1_VtxProb);
   filteredTree->Branch("Jpsi_1_px", &filtered_Jpsi_1_px);
   filteredTree->Branch("Jpsi_1_py", &filtered_Jpsi_1_py);
   filteredTree->Branch("Jpsi_1_pz", &filtered_Jpsi_1_pz);
   filteredTree->Branch("Jpsi_1_phi", &filtered_Jpsi_1_phi);
   filteredTree->Branch("Jpsi_1_eta", &filtered_Jpsi_1_eta);
   filteredTree->Branch("Jpsi_1_pt", &filtered_Jpsi_1_pt);
   filteredTree->Branch("Jpsi_1_mu_1_Idx", &filtered_Jpsi_1_mu_1_Idx);
   filteredTree->Branch("Jpsi_1_mu_2_Idx", &filtered_Jpsi_1_mu_2_Idx);

   filteredTree->Branch("Jpsi_2_mass", &filtered_Jpsi_2_mass);
   filteredTree->Branch("Jpsi_2_massErr", &filtered_Jpsi_2_massErr);
   filteredTree->Branch("Jpsi_2_massDiff", &filtered_Jpsi_2_massDiff);
   filteredTree->Branch("Jpsi_2_ctau", &filtered_Jpsi_2_ctau);
   filteredTree->Branch("Jpsi_2_ctauErr", &filtered_Jpsi_2_ctauErr);
   filteredTree->Branch("Jpsi_2_Chi2", &filtered_Jpsi_2_Chi2);
   filteredTree->Branch("Jpsi_2_ndof", &filtered_Jpsi_2_ndof);
   filteredTree->Branch("Jpsi_2_VtxProb", &filtered_Jpsi_2_VtxProb);
   filteredTree->Branch("Jpsi_2_px", &filtered_Jpsi_2_px);
   filteredTree->Branch("Jpsi_2_py", &filtered_Jpsi_2_py);
   filteredTree->Branch("Jpsi_2_pz", &filtered_Jpsi_2_pz);
   filteredTree->Branch("Jpsi_2_phi", &filtered_Jpsi_2_phi);
   filteredTree->Branch("Jpsi_2_eta", &filtered_Jpsi_2_eta);
   filteredTree->Branch("Jpsi_2_pt", &filtered_Jpsi_2_pt);
   filteredTree->Branch("Jpsi_2_mu_1_Idx", &filtered_Jpsi_2_mu_1_Idx);
   filteredTree->Branch("Jpsi_2_mu_2_Idx", &filtered_Jpsi_2_mu_2_Idx);

   filteredTree->Branch("Pri_mass", &filtered_Pri_mass);
   filteredTree->Branch("Pri_massErr", &filtered_Pri_massErr);
   filteredTree->Branch("Pri_ctau", &filtered_Pri_ctau);
   filteredTree->Branch("Pri_ctauErr", &filtered_Pri_ctauErr);
   filteredTree->Branch("Pri_Chi2", &filtered_Pri_Chi2);
   filteredTree->Branch("Pri_ndof", &filtered_Pri_ndof);
   filteredTree->Branch("Pri_VtxProb", &filtered_Pri_VtxProb);
   filteredTree->Branch("Pri_px", &filtered_Pri_px);
   filteredTree->Branch("Pri_py", &filtered_Pri_py);
   filteredTree->Branch("Pri_pz", &filtered_Pri_pz);
   filteredTree->Branch("Pri_phi", &filtered_Pri_phi);
   filteredTree->Branch("Pri_eta", &filtered_Pri_eta);
   filteredTree->Branch("Pri_pt", &filtered_Pri_pt);

   filteredTree->Branch("Phi_mass", &filtered_Phi_mass);
   filteredTree->Branch("Phi_massErr", &filtered_Phi_massErr);
   filteredTree->Branch("Phi_massDiff", &filtered_Phi_massDiff);
   filteredTree->Branch("Phi_ctau", &filtered_Phi_ctau);
   filteredTree->Branch("Phi_ctauErr", &filtered_Phi_ctauErr);
   filteredTree->Branch("Phi_Chi2", &filtered_Phi_Chi2);
   filteredTree->Branch("Phi_ndof", &filtered_Phi_ndof);
   filteredTree->Branch("Phi_VtxProb", &filtered_Phi_VtxProb);
   filteredTree->Branch("Phi_px", &filtered_Phi_px);
   filteredTree->Branch("Phi_py", &filtered_Phi_py);
   filteredTree->Branch("Phi_pz", &filtered_Phi_pz);
   filteredTree->Branch("Phi_phi", &filtered_Phi_phi);
   filteredTree->Branch("Phi_eta", &filtered_Phi_eta);
   filteredTree->Branch("Phi_pt", &filtered_Phi_pt);
   filteredTree->Branch("Phi_K_1_Idx", &filtered_Phi_K_1_Idx);
   filteredTree->Branch("Phi_K_2_Idx", &filtered_Phi_K_2_Idx);

   // Register all muon-related branches
   filteredTree->Branch("Jpsi_1_mu_1_px", &filtered_Jpsi_1_mu_1_px);
   filteredTree->Branch("Jpsi_1_mu_1_py", &filtered_Jpsi_1_mu_1_py);
   filteredTree->Branch("Jpsi_1_mu_1_pz", &filtered_Jpsi_1_mu_1_pz);
   filteredTree->Branch("Jpsi_1_mu_1_phi", &filtered_Jpsi_1_mu_1_phi);
   filteredTree->Branch("Jpsi_1_mu_1_eta", &filtered_Jpsi_1_mu_1_eta);
   filteredTree->Branch("Jpsi_1_mu_1_pt", &filtered_Jpsi_1_mu_1_pt);
   filteredTree->Branch("Jpsi_1_mu_1_isPatLooseMuon", &filtered_Jpsi_1_mu_1_isPatLooseMuon);
   filteredTree->Branch("Jpsi_1_mu_1_isPatSoftMuon", &filtered_Jpsi_1_mu_1_isPatSoftMuon);
   filteredTree->Branch("Jpsi_1_mu_1_isPatMediumMuon", &filtered_Jpsi_1_mu_1_isPatMediumMuon);
   filteredTree->Branch("Jpsi_1_mu_1_isPatTightMuon", &filtered_Jpsi_1_mu_1_isPatTightMuon);
   
   filteredTree->Branch("Jpsi_1_mu_2_px", &filtered_Jpsi_1_mu_2_px);
   filteredTree->Branch("Jpsi_1_mu_2_py", &filtered_Jpsi_1_mu_2_py);
   filteredTree->Branch("Jpsi_1_mu_2_pz", &filtered_Jpsi_1_mu_2_pz);
   filteredTree->Branch("Jpsi_1_mu_2_phi", &filtered_Jpsi_1_mu_2_phi);
   filteredTree->Branch("Jpsi_1_mu_2_eta", &filtered_Jpsi_1_mu_2_eta);
   filteredTree->Branch("Jpsi_1_mu_2_pt", &filtered_Jpsi_1_mu_2_pt);
   filteredTree->Branch("Jpsi_1_mu_2_isPatLooseMuon", &filtered_Jpsi_1_mu_2_isPatLooseMuon);
   filteredTree->Branch("Jpsi_1_mu_2_isPatSoftMuon", &filtered_Jpsi_1_mu_2_isPatSoftMuon);
   filteredTree->Branch("Jpsi_1_mu_2_isPatMediumMuon", &filtered_Jpsi_1_mu_2_isPatMediumMuon);
   filteredTree->Branch("Jpsi_1_mu_2_isPatTightMuon", &filtered_Jpsi_1_mu_2_isPatTightMuon);

   filteredTree->Branch("Jpsi_2_mu_1_px", &filtered_Jpsi_2_mu_1_px);
   filteredTree->Branch("Jpsi_2_mu_1_py", &filtered_Jpsi_2_mu_1_py);
   filteredTree->Branch("Jpsi_2_mu_1_pz", &filtered_Jpsi_2_mu_1_pz);
   filteredTree->Branch("Jpsi_2_mu_1_phi", &filtered_Jpsi_2_mu_1_phi);
   filteredTree->Branch("Jpsi_2_mu_1_eta", &filtered_Jpsi_2_mu_1_eta);
   filteredTree->Branch("Jpsi_2_mu_1_pt", &filtered_Jpsi_2_mu_1_pt);
   filteredTree->Branch("Jpsi_2_mu_1_isPatLooseMuon", &filtered_Jpsi_2_mu_1_isPatLooseMuon);
   filteredTree->Branch("Jpsi_2_mu_1_isPatSoftMuon", &filtered_Jpsi_2_mu_1_isPatSoftMuon);
   filteredTree->Branch("Jpsi_2_mu_1_isPatMediumMuon", &filtered_Jpsi_2_mu_1_isPatMediumMuon);
   filteredTree->Branch("Jpsi_2_mu_1_isPatTightMuon", &filtered_Jpsi_2_mu_1_isPatTightMuon);

   filteredTree->Branch("Jpsi_2_mu_2_px", &filtered_Jpsi_2_mu_2_px);
   filteredTree->Branch("Jpsi_2_mu_2_py", &filtered_Jpsi_2_mu_2_py);
   filteredTree->Branch("Jpsi_2_mu_2_pz", &filtered_Jpsi_2_mu_2_pz);
   filteredTree->Branch("Jpsi_2_mu_2_phi", &filtered_Jpsi_2_mu_2_phi);
   filteredTree->Branch("Jpsi_2_mu_2_eta", &filtered_Jpsi_2_mu_2_eta);
   filteredTree->Branch("Jpsi_2_mu_2_pt", &filtered_Jpsi_2_mu_2_pt);
   filteredTree->Branch("Jpsi_2_mu_2_isPatLooseMuon", &filtered_Jpsi_2_mu_2_isPatLooseMuon);
   filteredTree->Branch("Jpsi_2_mu_2_isPatSoftMuon", &filtered_Jpsi_2_mu_2_isPatSoftMuon);
   filteredTree->Branch("Jpsi_2_mu_2_isPatMediumMuon", &filtered_Jpsi_2_mu_2_isPatMediumMuon);
   filteredTree->Branch("Jpsi_2_mu_2_isPatTightMuon", &filtered_Jpsi_2_mu_2_isPatTightMuon);

   filteredTree->Branch("Phi_K_1_px", &filtered_Phi_K_1_px);
   filteredTree->Branch("Phi_K_1_py", &filtered_Phi_K_1_py);
   filteredTree->Branch("Phi_K_1_pz", &filtered_Phi_K_1_pz);
   filteredTree->Branch("Phi_K_1_phi", &filtered_Phi_K_1_phi);
   filteredTree->Branch("Phi_K_1_eta", &filtered_Phi_K_1_eta);
   filteredTree->Branch("Phi_K_1_pt", &filtered_Phi_K_1_pt);
   
   filteredTree->Branch("Phi_K_2_px", &filtered_Phi_K_2_px);
   filteredTree->Branch("Phi_K_2_py", &filtered_Phi_K_2_py);
   filteredTree->Branch("Phi_K_2_pz", &filtered_Phi_K_2_pz);
   filteredTree->Branch("Phi_K_2_phi", &filtered_Phi_K_2_phi);
   filteredTree->Branch("Phi_K_2_eta", &filtered_Phi_K_2_eta);
   filteredTree->Branch("Phi_K_2_pt", &filtered_Phi_K_2_pt);
}
#endif // #ifdef secCut_cxx
