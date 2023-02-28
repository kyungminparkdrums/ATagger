#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

struct pho_obj {
  unsigned int idx;
  double pt;
};


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesDiPhotonSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("m0",        &m0_);
  tree->Branch("FC_inputs", &vFC_inputs_);
  tree->Branch("hltAccept", &hltAccept_);
  //tree->Branch("nPreselPho",  &nPreselPho_);
  tree->Branch("nRecoPho",  &nRecoPho_);
  tree->Branch("minDR",     &vMinDR_);

  hNpassed_kin      = fs->make<TH1F>("hNpassed_kin", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_presel   = fs->make<TH1F>("hNpassed_presel", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_mGG      = fs->make<TH1F>("hNpassed_mGG", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_nRecoPho = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_hlt      = fs->make<TH1F>("hNpassed_hlt", "isPassed;isPassed;N", 2, 0., 2);
}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runDiPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  ////////// Apply selection //////////

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
  vTrigPhoIdxs_.clear();
  vPreselPhoIdxs_.clear();

  std::vector<pho_obj> vPhosNoCut;
  vPhoNoCutIdxs_.clear();
  vPhosNoCut.clear();

  std::vector<pho_obj> vPhos;
  vPhos.clear();

  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    PhotonRef iPho( photons, iP );
    vTrigPhoIdxs_.push_back( iP );

    vPhoNoCutIdxs_.push_back( iP );

    pho_obj Pho_obj = { iP, std::abs(iPho->pt()) };
    vPhosNoCut.push_back( Pho_obj );

    vPhos.push_back( Pho_obj );
    vPreselPhoIdxs_.push_back( iP );
  }
  std::cout << " vTrigPhoIdxs_ size: " << vTrigPhoIdxs_.size() << std::endl;

  //std::sort( vPhosNoCut.begin(), vPhosNoCut.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  //leadingPt = vPhosNoCut[0].pt;
  //subleadingPt = vPhosNoCut[1].pt;

  if ( hltAccept == 0 ) return false;  // KYUNGMIN
  //nRecoPho_ = vPreselPhoIdxs_.size(); // KYUNGMIN
  hNpassed_hlt->Fill(1.);  // KYUNGMIN

  if ( debug ) std::cout << " Presel pho size:" << vPhos.size() << std::endl;
  if ( vPhos.size() != 2 ) return false;
  hNpassed_presel->Fill(1.);

  std::sort( vPhos.begin(), vPhos.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  leadingPt = vPhos[0].pt;
  subleadingPt = vPhos[1].pt;

  nRecoPho_ = vPreselPhoIdxs_.size(); // KYUNGMIN

  return true;
}

// Fill branches ___________________________________________________________________//
void RecHitAnalyzer::fillDiPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // do we need this???

  // Fill kinematic variables
  float dphi[2] = { 0., 0. };
  vFC_inputs_.clear();
  for ( unsigned int iP = 0; iP < vPhoObj_recoIdx_.size(); iP++ ) {
    PhotonRef iPho( photons, vPhoObj_recoIdx_[iP] );
    vFC_inputs_.push_back( iPho->pt()/m0_ );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );

  // Get dR of closest reco photon to presel photon
  float minDR, dR, ikP;
  vMinDR_.clear();
  for ( unsigned int jP = 0; jP < vPhoObj_recoIdx_.size(); jP++ ) {
    PhotonRef jPho( photons, vPhoObj_recoIdx_[jP] );
    ikP = -1;
    minDR = 100.;
    for ( unsigned int kP = 0; kP < photons->size(); kP++ ) {
      PhotonRef kPho( photons, kP );

      if ( std::find(vPhoObj_recoIdx_.begin(), vPhoObj_recoIdx_.end(), kP) == vPhoObj_recoIdx_.end() ) continue;
      if ( std::abs(kPho->pt()) < 10. ) continue;

      dR = reco::deltaR( jPho->eta(),jPho->phi(), kPho->eta(),kPho->phi() );
      if ( dR < minDR ) {
	if ( debug ) std::cout << jP << " " << ikP << " dR " << dR << std::endl;
	minDR = dR;
	ikP = kP;
      }
    } //k
    if ( debug ) std::cout << jP << " " << ikP << " minDR " << minDR << std::endl;
    vMinDR_.push_back( minDR );
  } //j

}
