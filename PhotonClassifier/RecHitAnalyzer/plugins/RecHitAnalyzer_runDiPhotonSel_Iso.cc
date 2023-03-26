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

  // Require reco photon to be in barrel, and > 10 GeV pT
  vRecoPhoIdxs_.clear();

  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    PhotonRef iPho( photons, iP );
    if ( std::abs(iPho->pt()) < 10. ) continue;
    if ( std::abs(iPho->eta()) > 1.442 ) continue;
    vRecoPhoIdxs_.push_back( iP );
  }

  std::cout << " vRecoPhoIdxs_ size: " << vRecoPhoIdxs_.size() << std::endl;

  // Apply Trigger
  hNpassed_hlt->Fill(0.);

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

  vTrigPhoIdxs_.clear();
  vPreselPhoIdxs_.clear();

  std::vector<pho_obj> vPhos;
  vPhos.clear();

  if ( hltAccept == 0 ) return false;
  hNpassed_hlt->Fill(1.);

  for ( unsigned int iP : vRecoPhoIdxs_ ) {
    PhotonRef iPho( photons, iP );
    
    vTrigPhoIdxs_.push_back( iP );
  }
  std::cout << " vTrigPhoIdxs_ size: " << vTrigPhoIdxs_.size() << std::endl;

  if ( vTrigPhoIdxs_.size() != 2 ) return false;  // require two barrel reco photons that pass the trigger

  // Apply EGamma Iso cuts
  hNpassed_presel->Fill(0.);

  for ( unsigned int iP : vTrigPhoIdxs_ ) {
    std::cout << "iP in vTrigPhoIdxs_: " << iP << std::endl;
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
  if ( vPhos.size() != 2 ) return false;  // require two barrel reco photons that pass the iso cuts
  hNpassed_presel->Fill(1.);

  std::sort( vPhos.begin(), vPhos.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  leadingPt = vPhos[0].pt;
  subleadingPt = vPhos[1].pt;

  hLeadingPhoPt->Fill(leadingPt);
  hSubleadingPhoPt->Fill(subleadingPt);

  //nRecoPho_ = vPreselPhoIdxs_.size();

  return true;
}


