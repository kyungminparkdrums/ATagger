#include <vector>
#include <TFile.h>
#include <iostream>
using namespace std;

void cutflow() {

  TFile *f = TFile::Open("root://cmsdata.phys.cmu.edu//store/user/kypark/ATagger_GJ_40toInf_pThist/Dec4_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Dec4_2022_GJ_pT40toInf/221204_074523/0000/output_1-1.root");
  //TFile *f = TFile::Open("root://cmsdata.phys.cmu.edu//store/user/kypark/ATagger_GJ_40toInf_pThist/Dec4_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Dec4_2022_GJ_pT40toInf/221204_074523/0000/output_1-1.root");
  //TFile *f = TFile::Open("root://cmsdata.phys.cmu.edu//store/user/kypark/ATagger_GJ_40toInf/Nov28_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov28_2022_GJ_pT40toInf/221129_050015/0000/output_1-9.root ");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov20_GJ/Nov20_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov20_2022_GJ_pT40toInf_All/221124_212252/0000/output_1-70.root");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov24_NoDeltaR_AA_Full/Nov24_2022_NoDeltaR/HAHMHToAA_AToGG_MA-0p1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov24_2022_NoDeltaR_AA_0p1GeV/221125_002134/0000/output_3.root");
//  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov20/Nov20_2022_NoDeltaR/HAHMHToAA_AToGG_MA-0p8GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov20_2022_NoDeltaR_AA_0p8GeV/221124_202418/0000/AA_NoDeltaR_0p8.root");
//  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_HLT_PhoCut/Nov11_2022/HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov11_2022_AA_HLT_PhoCut_0p6/221111_193428/0000/output_1.root");

  TDirectory *dir = (TDirectory*)f->Get("fevt");

  TH1F *hNpassed_hlt = (TH1F*)dir->Get("hNpassed_hlt");
  TH1F *hNpassed_img = (TH1F*)dir->Get("hNpassed_img");


  TH1F *hNRecoPho_HLT = (TH1F*)dir->Get("hNRecoPho_HLT");
  TH1F *hNRecoPho_hits = (TH1F*)dir->Get("hNRecoPho_hits");
  TH1F *hNRecoPho_genMatching = (TH1F*)dir->Get("hNRecoPho_genMatching");

  

  int nAnalyzed = hNpassed_hlt->GetBinContent(1);
  int nPassTrigger = hNpassed_img->GetBinContent(1);
  int nPassPresel = hNpassed_img->GetBinContent(2);

  cout << "Number of files analyzed: " << nAnalyzed << endl;
  cout << "After Trigger: " << nPassTrigger << endl;
  cout << "After Photon Preselection: " << nPassPresel << endl;


  int nPho_Trigger = 0;
  int nPho_Presel = 0;

  for (int i=1; i <= hNRecoPho_hits->GetNbinsX(); i++) {
  //for (int i=1; i <= hNRecoPho_HLT->GetNbinsX(); i++) {
    //cout << i-1 << " photons : " << hNRecoPho_HLT->GetBinContent(i) << endl;
    cout << i-1 << " photons : " << hNRecoPho_hits->GetBinContent(i) << endl;
    nPho_Presel += (i-1)*hNRecoPho_hits->GetBinContent(i);
    //nPho_Trigger += (i-1)*hNRecoPho_HLT->GetBinContent(i);
  }

  cout << "nPhotons after pre-selection: " << nPho_Presel << endl;
  //cout << "nPhotons after trigger: " << nPho_Trigger << endl;

/*
  TCanvas c1;
  hNRecoPho_HLT->Draw();
  c1.SaveAs("hNrecoPho_HLT.pdf", "pdf");

  hNRecoPho_hits->Draw();
  c1.SaveAs("hNrecoPho_hits.pdf", "pdf");

  hNRecoPho_genMatching->Draw();
  c1.SaveAs("hNrecoPho_genMatching.pdf", "pdf");
*/
}
