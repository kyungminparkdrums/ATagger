#include <vector>
#include <TFile.h>
#include <string>
#include <iostream>
using namespace std;

void plotE() {

  TFile *f = new TFile("/eos/uscms/store/user/kyungmip/Ng2.root");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/output_1.root");
  //TFile *f = new TFile("/uscms/home/kyungmip/nobackup/CMSSW_10_6_18_GJ_PhoIso_Ngamma2/src/output_numEvent1000.root");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/output_1-2.root");
  //TFile *f = new TFile("/eos/uscms/store/user/kyungmip/ATagger_Nov27_AA/Nov27_2022/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov27_2022_AA_1GeV/221127_084716/0000/output_1.root");

  TDirectory *dir = (TDirectory*)f->Get("fevt");
  TTree *tree = (TTree*)dir->Get("RHTree");

  int nEntries = tree->GetEntries();

  cout << nEntries << endl;

  //tree->Print();
  std::vector<std::vector<float>> *RHimg_energy = 0;

  tree->SetBranchAddress("RHimg_energy", &RHimg_energy);

  //TH2D* hRechits[14];
  TH2D* hRechits = new TH2D("hRechits", "EB RecHits (32x32)", 32, 0, 32, 32, 0, 32);

  //for (int i=0; i < 1; i++) {
  for (int i=0; i < 100; i++) {
  //for (int i=0; i < nEntries; i++) {
    tree->GetEntry(i);

    //char name[20];
    //sprintf(name, "hRechits_%d", i);
    //hRechits[i] = new TH2D(name, "EB RecHits Image (32x32)", 32, 0, 32, 32, 0, 32);

   

    for (int k = 0; k < RHimg_energy->size(); k++) {
      for (int idx_ = 0; idx_ < (RHimg_energy->at(k)).size(); idx_++) {
        double hit = (RHimg_energy->at(k)).at(idx_);
        int ieta = idx_ / 32;
        int iphi = idx_ % 32;

        if (hit != 0) {
          //std::cout << hit << std::endl;
          //std::cout << "ieta = " << ieta << ", iphi = " << iphi << std::endl;
          //std::cout << "idx_: " << idx_ << std::endl;

          //hRechits[i]->Fill(iphi, ieta, hit);          
          hRechits->Fill(iphi, ieta, hit);          
        }
      }
    }
//    char saveName[20];
//    sprintf(saveName, "hRechits_%d_logZ.pdf", i);
//    TCanvas c1;
//    hRechits[i]->GetXaxis()->SetTitle("iPhi");
//    hRechits[i]->GetYaxis()->SetTitle("iEta");
//    hRechits[i]->Draw();
//    hRechits[i]->Draw("colz");

//    c1.SetLogz();
//    c1.SaveAs(saveName, "pdf");
    //break;
  }

    Double_t scale = 1./100;
    //Double_t scale = 1./(hRechits->Integral());
    hRechits->Scale(scale);
    
    TCanvas c1;
    gStyle->SetOptStat(0);
    c1.SetLogz();
    hRechits->GetXaxis()->SetTitle("iPhi");
    hRechits->GetYaxis()->SetTitle("iEta");
    hRechits->Draw("colz");
    c1.SaveAs("GJ_100_LogEvents_avg.pdf", "pdf");
}
