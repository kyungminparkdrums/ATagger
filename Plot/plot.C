#include <vector>
#include <TFile.h>
#include <iostream>
using namespace std;

void plot() {

  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/TMP/output_1-1.root");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov26_GJ_AncestorMatching/Nov26_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov26_2022_GJ_pT40toInf_All/221126_233246/0000/output_1-1.root");


  TDirectory *dir = (TDirectory*)f->Get("fevt");
  TTree *tree = (TTree*)dir->Get("RHTree");

  int nEntries = tree->GetEntries();

  cout << nEntries << endl;

  //tree->Print();
  ULong64_t eventId;
  std::vector<float> *pho_E = 0; 

  tree->SetBranchAddress("pho_E", &pho_E);

  TH1D* hPhoE = new TH1D("hPhoE", "Photon E", 100, 0, 200);

  //for (int i=0; i < 1000; i++) {
  for (int i=0; i < nEntries; i++) {
    tree->GetEntry(i);

    int n_pho_E = pho_E->size();
    float *data_pho_E = pho_E->data();

    for (int j=0; j < n_pho_E; j++) {
      hPhoE->Fill(data_pho_E[j]);    
    }
  }

  TFile* outFile = new TFile("outFile0.root", "RECREATE");

  outFile->cd();
  
  hPhoE->Write();  
    
  outFile->Close();

/*
  TCanvas c1;
  //hPhoE->GetXaxis()->SetTitle("E [GeV]");
  //hPhoE->GetYaxis()->SetTitle("Events / 2 GeV");
  //hPhoE->Draw("colz");
  c1.SaveAs("E_colz.pdf", "pdf");

  //tree->ResetBranchAddresses();
  //delete pho_E;
*/
}
