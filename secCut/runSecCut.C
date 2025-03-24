#include "/home/storage2/users/xingcheng/CMSSW_14_0_18/src/JpsiJpsiPhi-workspace/includes/ParticleCand.C"
#include "secCut.C"

//#define RUN_JOB

void runSecCut(){
    TChain *chain = new TChain("outputTree","");
    #ifdef RUN_JOB
    chain->Add("JOB_DATA");
    #else
    chain->Add("/home/storage2/users/xingcheng/CMSSW_14_0_18/src/JpsiJpsiPhi-workspace/preCut/preCut_Run3all_first.root");
    #endif
    secCut mySecCut(chain);
    mySecCut.Loop();
}
