#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

struct gen_obj {
  unsigned int idx;
  double pt;
};

std::vector<gen_obj> vAs;
std::vector<unsigned int> vGenAIdxs;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesPhoObjSel_AA ( TTree* tree, edm::Service<TFileService> &fs )
{
  
  //tree->Branch("PhoObj",                  &vPhoObj_);
  tree->Branch("PhoObj_recoIdx",          &vPhoObj_recoIdx_);
  tree->Branch("PhoObj_genMatchedIdx",    &vPhoObj_genMatchedIdx_);
  tree->Branch("PhoObj_genMatchedPdgId",  &vPhoObj_genMatchedPdgId_);
  tree->Branch("PhoObj_AncestorPdgId",    &vPhoObj_ancestorPdgId_);

  tree->Branch("vA_E", &vA_E_);  
  tree->Branch("vA_pT", &vA_pT_);  
  tree->Branch("vA_eta", &vA_eta_);  
  tree->Branch("vA_phi", &vA_phi_);  
  tree->Branch("vA_mass", &vA_mass_);  
  tree->Branch("vA_DR", &vA_DR_);  

}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runPhoObjSel_AA ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  std::vector<int> vGenPdgIds;
  std::vector<int> vAncestorIds;

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  ////////// Apply selection //////////

  vAs.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    if ( std::abs(iGen->pdgId()) != 35 ) continue;
    //if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 25 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    if ( iGen->daughter(0)->pdgId() != 22 || iGen->daughter(1)->pdgId() != 22 ) continue;

    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    vAs.push_back( Gen_obj );
    vH += iGen->p4();

  } // gen particles
  if ( vAs.size() != 2 ) return false;
  //mHgen_ = vH.mass();

  // Sort As by pT, for abitrary N
  //std::sort( vAs.begin(), vAs.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );
  vGenAIdxs.clear();

    for ( unsigned int iG = 0; iG < vAs.size(); iG++ ) {
      reco::GenParticleRef iGen( genParticles, vAs[iG].idx );
      if ( debug ) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;
    //vPreselPhoIdxs_.push_back( vAs[iG].idx );
      vGenAIdxs.push_back( vAs[iG].idx );
    }

  if ( debug ) std::cout << " Pho collection size:" << photons->size() << std::endl;

  std::cout << "vRegressPhoIdxs_ size: " << vRegressPhoIdxs_.size() << std::endl;
  std::cout << "vGenAIdxs size: " << vGenAIdxs.size() << std::endl;

  vPhoObj_recoIdx_.clear();
  vPhoObj_ancestorPdgId_.clear();
  vPhoObj_genMatchedPdgId_.clear();
  vPhoObj_genMatchedIdx_.clear();

  vA_E_.clear();
  vA_pT_.clear();
  vA_eta_.clear();
  vA_phi_.clear();
  vA_mass_.clear();
  vA_DR_.clear();

  std::cout << "CLEARED VECTORS" << std::endl;

  float dPhi, dEta, dR, recoDR;
  int recoDR_idx;

  int s23Idx = -1;
  int ancIdx = -1;
  int genId = 10000000;
 
  float  minDR = 10000000.;
  for ( unsigned int iG : vGenAIdxs ) {

    std::cout << "WILL GO OVER A GENS" << std::endl;
    reco::GenParticleRef iGen( genParticles, iG );

    vA_E_.push_back( std::abs(iGen->energy()) );
    vA_pT_.push_back( std::abs(iGen->pt()) );
    vA_eta_.push_back( iGen->eta() );
    vA_phi_.push_back( iGen->phi() );
    vA_mass_.push_back( iGen->mass() );
    vA_DR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );

    //dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );
    //dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    //hdPhidEta->Fill( dPhi, dEta );

    // Get index to dR-matched preselected photon
    recoDR = 2*0.04;
    recoDR_idx = -1;
  
    minDR = 10000000.;
  
    std::cout << "WILL LOOP THROUGH REGRESS PHO" << std::endl; 

    //for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
 
    for ( int iP : vRegressPhoIdxs_ ) {
      //PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
      PhotonRef iPho( photons, iP );
      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if (dR < minDR) {
        std::cout << "dR smaller than minDR!!" << std::endl;
        recoDR_idx = iP;
        minDR = dR;
        }
     std::cout << "LOOPED THROUGH REGRESS PHOTONS!!" << std::endl;
    }
    if (minDR < recoDR) {
      genId = iGen->pdgId();

      std::cout << "PUSH BACK PHOTONS" << std::endl;
      vPhoObj_recoIdx_.push_back( recoDR_idx );
      vPhoObj_ancestorPdgId_.push_back( genId );
      vPhoObj_genMatchedPdgId_.push_back( genId );
      vPhoObj_genMatchedIdx_.push_back( s23Idx );
    }
  } // reco pho

  // need to return bool
  return true;
}

