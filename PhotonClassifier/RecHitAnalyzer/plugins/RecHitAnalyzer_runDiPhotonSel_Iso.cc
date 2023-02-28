#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

struct pho_obj {
  unsigned int idx;
  double pt;
};


// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runDiPhotonSel_Iso ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  ////////// Apply selection //////////

/* COMMENTED OUT BY KYUNGMIN
  // Count number of "reco" photons
  std::vector<unsigned int> vRecoPhoIdxs;
  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    PhotonRef iPho( photons, iP );
    //if ( std::abs(iPho->pt()) < 5. ) continue;
    if ( std::abs(iPho->pt()) < 10. ) continue;
    vRecoPhoIdxs.push_back( iP );
  }
  if ( debug ) std::cout << " Reco pho size:" << vRecoPhoIdxs.size() << std::endl;
  nRecoPho_ = vRecoPhoIdxs.size();

  // Ensure at least 2 kinematic trigger-like photons
  hNpassed_kin->Fill(0.);
  std::vector<unsigned int> vKinPhoIdxs;
  for ( unsigned int iP : vRecoPhoIdxs ) {
    PhotonRef iPho( photons, iP );
    if ( debug ) std::cout << "Photon: " << iP << " pT: " << std::abs(iPho->pt()) << " eta: " << std::abs(iPho->eta()) << std::endl;
    if ( std::abs(iPho->pt()) <= 18. ) continue;
    if ( std::abs(iPho->eta()) >= 1.442 ) continue;
    vKinPhoIdxs.push_back( iP );
  }
  if ( debug ) std::cout << " Kinetic pho size:" << vKinPhoIdxs.size() << std::endl;
  if ( vKinPhoIdxs.size() < 2 ) return false;
  hNpassed_kin->Fill(1.);

  // Ensure two presel photons
  hNpassed_presel->Fill(0.);
  std::vector<pho_obj> vPhos;
  for ( unsigned int iP : vKinPhoIdxs ) {

    PhotonRef iPho( photons, iP );
    if ( debug ) std::cout << "Photon: " << iP << " full5x5_r9: " << std::abs(iPho->full5x5_r9()) << " hadTowOverEm: " << std::abs(iPho->hadTowOverEm()) << " hasPixelSeed: " << iPho->hasPixelSeed() << std::endl;

    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;

    if ( debug ) std::cout << "sigmaIetaIeta: " << std::abs(iPho->full5x5_sigmaIetaIeta()) << " userFloat: " << std::abs(iPho->userFloat("phoPhotonIsolation")) << " trkSumPt: " << iPho->trkSumPtHollowConeDR03() << std::endl;

    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
    }
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
    
    pho_obj Pho_obj = { iP, std::abs(iPho->pt()) };
    vPhos.push_back( Pho_obj );

  } // kinematic photons
  if ( debug ) std::cout << " Presel pho size:" << vPhos.size() << std::endl;
  if ( vPhos.size() != 2 ) return false;
  hNpassed_presel->Fill(1.);

  // Sort photons by pT, for abitrary N
  std::sort( vPhos.begin(), vPhos.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );
  for ( unsigned int iP = 0; iP < vPhos.size(); iP++ ) {
    PhotonRef iPho( photons, vPhos[iP].idx );
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
  }

  // Check if any photon pairing passes invariant mass cut
  hNpassed_mGG->Fill(0.);
  std::vector<int> vPhoIdxs;
  bool passedMassCut = false;
  for ( unsigned int j = 0; j < vPhos.size()-1; j++ ) {

    PhotonRef jPho( photons, vPhos[j].idx );

    for ( unsigned int k = 1; k < vPhos.size(); k++ ) {

      if ( k <= j ) continue;
      PhotonRef kPho( photons, vPhos[k].idx );
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vDiPho = jPho->p4() + kPho->p4();
      if ( debug ) std::cout << " >> m0:" << vDiPho.mass() << std::endl;
      if ( debug ) std::cout << k << " " << j << " >> m0:" << vDiPho.mass() << std::endl;

      if ( vDiPho.mass() > 90. ) {
	vPhoIdxs.push_back( vPhos[j].idx );
        vPhoIdxs.push_back( vPhos[k].idx );
        m0_ = vDiPho.mass();
        passedMassCut = true;
        break;
      }

    } //k
    if ( passedMassCut ) break;
  } // j
  if ( !passedMassCut ) return false;
  if ( debug ) std::cout << " >> m0:" << m0_ << std::endl;


  // Apply diphoton pT cuts
  float ptCut[2]   = { 30., 18. };
  float ptOmCut[2] = {  3.,  4. };
  vPreselPhoIdxs_.clear();
  for ( unsigned int iP = 0; iP < vPhoIdxs.size(); iP++ ) {

    PhotonRef iPho( photons, vPhoIdxs[iP] );

    if ( std::abs(iPho->pt()) < ptCut[iP] ) continue;
    if ( std::abs(iPho->pt()) < m0_/ptOmCut[iP] ) continue;
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
    
    vPreselPhoIdxs_.push_back( vPhoIdxs[iP] );

  } // vPhoIdxs
  if ( vPreselPhoIdxs_.size() != 2 ) return false;
  if ( debug ) std::cout << " Reco pho size:" << vPhos.size() << std::endl;
  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  hNpassed_mGG->Fill(1.);

*/


  hNpassed_hlt->Fill(0.);  // KYUNGMIN

  // Check HLT trigger decision
  edm::Handle<edm::TriggerResults> trgs;
  iEvent.getByToken( trgResultsT_, trgs );

  const edm::TriggerNames &triggerNames = iEvent.triggerNames( *trgs );
  if ( debug ) std::cout << " N triggers:" << trgs->size() << std::endl;
  for ( unsigned int iT = 0; iT < trgs->size(); iT++ ) {
    if ( debug ) std::cout << " name["<<iT<<"]:"<<triggerNames.triggerName(iT)<< std::endl;
  }

  int hltAccept = -1;
  //std::string trgName = "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v*";
  std::string trgName = "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*";
  std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), trgName );
  if ( debug ) std::cout << " N matches: " << trgMatches.size() << std::endl;

  if ( !trgMatches.empty() ) {

    hltAccept = 0;
    for ( auto const& iT : trgMatches ) {
      if ( debug ) std::cout << " name["<<triggerNames.triggerIndex(*iT)<<"]:"<< *iT << " -> " << trgs->accept(triggerNames.triggerIndex(*iT)) << std::endl;
      if ( trgs->accept(triggerNames.triggerIndex(*iT)) ) hltAccept = 1;
      break;
    }
  }
  hltAccept_ = hltAccept;

  // ADDED BY KYUNGMIN
  //std::cout << "hltAccept: " << hltAccept << std::endl;

  // ADDED BY KYUNGMIN
  vTrigPhoIdxs_.clear();
  vPreselPhoIdxs_.clear();

  std::vector<pho_obj> vPhosNoCut;
  vPhoNoCutIdxs_.clear();
  vPhosNoCut.clear();

  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    PhotonRef iPho( photons, iP );
    vTrigPhoIdxs_.push_back( iP );

    vPhoNoCutIdxs_.push_back( iP );

    pho_obj Pho_obj = { iP, std::abs(iPho->pt()) };
    vPhosNoCut.push_back( Pho_obj );
  }
  std::cout << " vTrigPhoIdxs_ size: " << vTrigPhoIdxs_.size() << std::endl;
  std::cout << " vPhoNoCutIdxs_ size: " << vPhoNoCutIdxs_.size() << std::endl;

  //std::sort( vPhosNoCut.begin(), vPhosNoCut.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  //leadingPt = vPhosNoCut[0].pt;
  //subleadingPt = vPhosNoCut[1].pt;

  if ( hltAccept == 0 ) return false;  // KYUNGMIN
  //nRecoPho_ = vPreselPhoIdxs_.size(); // KYUNGMIN
  hNpassed_hlt->Fill(1.);  // KYUNGMIN

  // Ensure two presel photons
  hNpassed_presel->Fill(0.);
  std::vector<pho_obj> vPhos;
  for ( unsigned int iP : vTrigPhoIdxs_ ) {

    PhotonRef iPho( photons, iP );
    if ( debug ) std::cout << "Photon: " << iP << " full5x5_r9: " << std::abs(iPho->full5x5_r9()) << " hadTowOverEm: " << std::abs(iPho->hadTowOverEm()) << " hasPixelSeed: " << iPho->hasPixelSeed() << std::endl;

    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;

    if ( debug ) std::cout << "sigmaIetaIeta: " << std::abs(iPho->full5x5_sigmaIetaIeta()) << " userFloat: " << std::abs(iPho->userFloat("phoPhotonIsolation")) << " trkSumPt: " << iPho->trkSumPtHollowConeDR03() << std::endl;

    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
    }
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
    
    pho_obj Pho_obj = { iP, std::abs(iPho->pt()) };
    vPhos.push_back( Pho_obj );

    vPreselPhoIdxs_.push_back( iP );
  } // kinematic photons

  if ( debug ) std::cout << " Presel pho size:" << vPhos.size() << std::endl;
  //if ( vPhos.size() != 2 ) return false;
  hNpassed_presel->Fill(1.);

  std::sort( vPhos.begin(), vPhos.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  if (vPhos.size() >= 2) {
  leadingPt = vPhos[0].pt;
  subleadingPt = vPhos[1].pt;
  }
  nRecoPho_ = vPreselPhoIdxs_.size(); // KYUNGMIN

  return true;
}

