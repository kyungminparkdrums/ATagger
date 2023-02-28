#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

struct gen_obj {
  unsigned int idx;
  double pt;
};

std::vector<gen_obj> vPhoObj;
std::vector<unsigned int> vGenIdxs;
std::vector<int> vGenPdgIds;
std::vector<int> vAncestorIds;

// check if gen particle has s=23 ancestor: s=23 is initial state particle
bool hasStatus23Ancestor ( const reco::Candidate* part ) {
  
  if (part->status() == 23) return true;
  if (part->numberOfMothers() == 0) return false;

  bool hasAncestor = false;
  for ( unsigned int i = 0; i < part->numberOfMothers(); i++ ) {
    hasAncestor = (hasAncestor || hasStatus23Ancestor(part->mother(i)));
  }

  return hasAncestor;
}

// find s=23 ancestor
const reco::Candidate* status23Ancestor ( const reco::Candidate* part ) {

  if (part->status() == 23) return part;
  if (part->numberOfMothers() == 0) return NULL;

  const reco::Candidate* ancestor = NULL;
  for ( unsigned int i = 0; i < part->numberOfMothers(); i++ ) {
    if (!ancestor) {
      ancestor = status23Ancestor(part->mother(i));
    }
  }

  return ancestor;
}

// recursively print all daughter particles
void printProgeny ( const reco::GenParticleRef part ) {
  void printProgeny_rec(const reco::Candidate* part, unsigned int l);
  std::cout << "L0 - ID: " << part->pdgId() << " status: " << part->status() << " pT: " << part->pt() << std::endl;

  for (unsigned int i = 0; i <part->numberOfDaughters(); i++) {
    printProgeny_rec(part->daughter(i), 0);
  }
}

// function for recursion above (type needs to be Candidate)
void printProgeny_rec ( const reco::Candidate* part, unsigned int l ) {
  std::cout << "L" << l+1 << " - ID: " << part->pdgId() << " status: " << part->status() << " pT: " << part->pt() << std::endl;
  
  for ( unsigned int i = 0; i < part->numberOfDaughters(); i++) {
    printProgeny_rec(part->daughter(i),l+1);
  }
}

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesPhoObjSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  
  //tree->Branch("PhoObj",                  &vPhoObj_);
  tree->Branch("PhoObj_recoIdx",          &vPhoObj_recoIdx_);
  tree->Branch("PhoObj_genMatchedIdx",    &vPhoObj_genMatchedIdx_);
  tree->Branch("PhoObj_genMatchedPdgId",  &vPhoObj_genMatchedPdgId_);
  tree->Branch("PhoObj_AncestorPdgId",    &vPhoObj_ancestorPdgId_);

}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runPhoObjSel_fake ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  ////////// Apply selection //////////

  if ( debug ) std::cout << " Pho collection size:" << photons->size() << std::endl;
  
  vPhoObj.clear();

  // Gather all initial state particle (23 ancestors) in vPhoObj
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
  
    reco::GenParticleRef iGen( genParticles, iG );
    if ( iGen->status() != 23 ) continue;
    
    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    vPhoObj.push_back( Gen_obj );

  } 
  
  vPhoObj_recoIdx_.clear();
  vPhoObj_ancestorPdgId_.clear();
  vPhoObj_genMatchedPdgId_.clear();
  vPhoObj_genMatchedIdx_.clear();
  vGenIdxs.clear();

  std::cout << "vRegressPhoIdxs_ size: " << vRegressPhoIdxs_.size() << std::endl;
 
  float dR;
  // iterate over photons that passed pre-selection (including barrel and all that)
  for ( int iP : vRegressPhoIdxs_ ) {
    
    PhotonRef iPho( photons, iP );
 
    int s23Idx = -1;
    int ancIdx = -1;
    int genId = 10000000;
    float minDR = 10000000.;

    if ( debug ) std::cout << iP << " Photon ID: " << iPho->pdgId() << " pT " << iPho->pt() << std::endl;                                                                                                                    
    const reco::Candidate* ancestor;
    for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
      
      const reco::Candidate* tmpAncestor = NULL;
      reco::GenParticleRef iGen( genParticles, iG );
      if (iGen->status() != 1) continue;  // only look at final state particle

      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());

      if (dR < minDR) {

        // find s = 23 ancestor	
        for ( unsigned int i = 0; i < iGen->numberOfMothers(); i++ ) {
          if (!tmpAncestor) {
            tmpAncestor = status23Ancestor(iGen->mother(i));
          }
        }

        // once you have matched the s = 23 ancestor
        if (tmpAncestor) {
	  gen_obj goAnc;
	  for ( gen_obj go : vPhoObj ) {
	    if ( go.pt == tmpAncestor->pt() ) {
	      goAnc = go;
	      break;
	    }
	  }
	  ancestor = tmpAncestor;
	  ancIdx = goAnc.idx;
	  genId = iGen->pdgId();
          s23Idx = iG;
          minDR = dR;
	}
      }
    }

    if ( ancestor->pdgId() == 22 ) {
      printProgeny_rec( ancestor, 0 );
    }
    
    reco::GenParticleRef fGen( genParticles, s23Idx );
   
    if ( minDR > 0.08 ) {
      std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
    } else if ( ancestor->pdgId() == 22 ) {
      std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
      std::cout << "Ancestor is photon - ignoring for GJet" << std::endl; //for GJet
    //} else if ( ancestor->pdgId() != 22 ) {
    //  std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
    //  std::cout << "Ancestor isn't photon - ignoring for DiPhoton" << std::endl; //for DiPhoton
    } else {
      std::cout << "Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
     
      // If  
      vPhoObj_recoIdx_.push_back( iP );
      vPhoObj_ancestorPdgId_.push_back( ancestor->pdgId() );
      vPhoObj_genMatchedPdgId_.push_back( genId );
      vPhoObj_genMatchedIdx_.push_back( s23Idx );
      vGenIdxs.push_back( ancIdx );
    }
  }
  
  // need to return bool
  return true;
}


bool RecHitAnalyzer::runPhoObjSel_prompt ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  ////////// Apply selection //////////

  if ( debug ) std::cout << " Pho collection size:" << photons->size() << std::endl;
  
  vPhoObj.clear();

  // Gather all initial state particle (23 ancestors) in vPhoObj
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
  
    reco::GenParticleRef iGen( genParticles, iG );
    if ( iGen->status() != 23 ) continue;
    
    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    vPhoObj.push_back( Gen_obj );

  } 
  
  vPhoObj_recoIdx_.clear();
  vPhoObj_ancestorPdgId_.clear();
  vPhoObj_genMatchedPdgId_.clear();
  vPhoObj_genMatchedIdx_.clear();
  vGenIdxs.clear();

  std::cout << "vRegressPhoIdxs_ size: " << vRegressPhoIdxs_.size() << std::endl;
 
  float dR;
  // iterate over photons that passed pre-selection (including barrel and all that)
  for ( int iP : vRegressPhoIdxs_ ) {
    
    PhotonRef iPho( photons, iP );
 
    int s23Idx = -1;
    int ancIdx = -1;
    int genId = 10000000;
    float minDR = 10000000.;

    if ( debug ) std::cout << iP << " Photon ID: " << iPho->pdgId() << " pT " << iPho->pt() << std::endl;                                                                                                                    
    const reco::Candidate* ancestor;
    for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
      
      const reco::Candidate* tmpAncestor = NULL;
      reco::GenParticleRef iGen( genParticles, iG );
      if (iGen->status() != 1) continue;  // only look at final state particle

      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());

      if (dR < minDR) {

        // find s = 23 ancestor	
        for ( unsigned int i = 0; i < iGen->numberOfMothers(); i++ ) {
          if (!tmpAncestor) {
            tmpAncestor = status23Ancestor(iGen->mother(i));
          }
        }

        // once you have matched the s = 23 ancestor
        if (tmpAncestor) {
	  gen_obj goAnc;
	  for ( gen_obj go : vPhoObj ) {
	    if ( go.pt == tmpAncestor->pt() ) {
	      goAnc = go;
	      break;
	    }
	  }
	  ancestor = tmpAncestor;
	  ancIdx = goAnc.idx;
	  genId = iGen->pdgId();
          s23Idx = iG;
          minDR = dR;
	}
      }
    }

    if ( ancestor->pdgId() == 22 ) {
      printProgeny_rec( ancestor, 0 );
    }
    
    reco::GenParticleRef fGen( genParticles, s23Idx );
   
    if ( minDR > 0.08 ) {
      std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
    //} else if ( ancestor->pdgId() == 22 ) {
    //  std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
    //  std::cout << "Ancestor is photon - ignoring for GJet" << std::endl; //for GJet
    } else if ( ancestor->pdgId() != 22 ) {
      std::cout << "REJECTED Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
      std::cout << "Ancestor isn't photon - ignoring for DiPhoton" << std::endl; //for DiPhoton
    } else {
      std::cout << "Min dR: " << minDR << " reco candidate idx " << iP << ", gen pdgId " << genId << " idx " << s23Idx << " ancestor truth ID " << ancestor->pdgId() << std::endl;
     
      // If  
      vPhoObj_recoIdx_.push_back( iP );
      vPhoObj_ancestorPdgId_.push_back( ancestor->pdgId() );
      vPhoObj_genMatchedPdgId_.push_back( genId );
      vPhoObj_genMatchedIdx_.push_back( s23Idx );
      vGenIdxs.push_back( ancIdx );
    }
  }
  
  // need to return bool
  return true;
}




// Fill branches ___________________________________________________________________//
void RecHitAnalyzer::fillPhoObjSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // do something if needed
}

