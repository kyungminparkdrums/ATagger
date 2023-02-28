#include <vector>
#include <TFile.h>
#include <iostream>
using namespace std;

void plotAll() {

  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov26_GJ_AncestorMatching/Nov26_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov26_2022_GJ_pT40toInf_All/221126_233246/0000/output_1-1.root");

//  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov20/Nov20_2022_NoDeltaR/HAHMHToAA_AToGG_MA-0p8GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov20_2022_NoDeltaR_AA_0p8GeV/221124_202418/0000/AA_NoDeltaR_0p8.root");
//  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_HLT_PhoCut/Nov11_2022/HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov11_2022_AA_HLT_PhoCut_0p6/221111_193428/0000/output_1.root");

  TDirectory *dir = (TDirectory*)f->Get("fevt");
  TTree *tree = (TTree*)dir->Get("RHTree");

  int nEntries = tree->GetEntries();

  cout << nEntries << endl;

  //tree->Print();
  ULong64_t eventId;
  std::vector<float> *pho_E = 0; 

  tree->SetBranchAddress("pho_E", &pho_E);

  TH1D* hPhoE = new TH1D("hPhoE", "Photon E", 100, 0, 200);

  for (int i=0; i < 1000; i++) {
  //for (int i=0; i < nEntries; i++) {
    tree->GetEntry(i);
    //cout << eventId << endl;
    //cout << pho_E[0] << endl;
    //break;

    int n_pho_E = pho_E->size();
    float *data_pho_E = pho_E->data();

    for (int j=0; j < n_pho_E; j++) {
      hPhoE->Fill(data_pho_E[j]);    
      //std::cout << "Entry: " << i << ", photon idx: " << j << ", E: " << data_pho_E[j] << std::endl;
    }
  }

  TCanvas c1;
  hPhoE->GetXaxis()->SetTitle("E [GeV]");
  hPhoE->GetYaxis()->SetTitle("Events / 2 GeV");
  hPhoE->Draw("colz");
  c1.SaveAs("E_colz.pdf", "pdf");

  //tree->ResetBranchAddresses();
  //delete pho_E;
}
