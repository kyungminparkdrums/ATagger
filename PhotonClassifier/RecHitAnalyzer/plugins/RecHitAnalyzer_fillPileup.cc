#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

Int_t            nPUInfo_;
std::vector<int>      nPU_;
std::vector<int>      puBX_;
std::vector<float>    puTrue_;
std::vector<int>      mPU_;
std::vector<float>    mPUTrue_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesPileup ( TTree* tree, edm::Service<TFileService> &fs )
{

  tree->Branch("nPUInfo",       &nPUInfo_);
  tree->Branch("nPU",           &nPU_);
  tree->Branch("puBX",          &puBX_);
  tree->Branch("puTrue",        &puTrue_);

  tree->Branch("mPU",       &mPU_);
  tree->Branch("mPUTrue",   &mPUTrue_);


} // branchesEvtWgt()

// Fill EvtWeight _________________________________________________________________//
void RecHitAnalyzer::fillPileup ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  // from: https://github.com/cmkuo/ggAnalysis/blob/94X/ggNtuplizer/plugins/ggNtuplizer_genParticles.cc#L293-L301

  nPUInfo_ = 0;
  nPU_.clear();
  puBX_.clear();
  puTrue_.clear();
  mPU_.clear();
  mPUTrue_.clear();

  edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
  
  if (genPileupHandle.isValid()) {
    for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        mPU_.push_back(pu->getPU_NumInteractions());
        mPUTrue_.push_back(pu->getTrueNumInteractions());
	std::cout << "Crossing: " << pu->getBunchCrossing() << " N int: " << pu->getPU_NumInteractions() << " true N int: " << pu->getTrueNumInteractions() << std::endl;
      }

      nPU_   .push_back(pu->getPU_NumInteractions());
      puTrue_.push_back(pu->getTrueNumInteractions());
      puBX_  .push_back(pu->getBunchCrossing());

      nPUInfo_++;
    }
  } else {
    std::cout << "no PileupSummaryInfo in event" << std::endl;
  }
  
} // fillPileup()
