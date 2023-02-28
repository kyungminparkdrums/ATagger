// --------------------------------------------
//
// Package:    PhotonClassifier
// Class:      RecHitAnalyzer
// 
//
// Author: Andrew C. Roberts
// Started 2022/5/18
// Last Updated 2022/5/18
//
// --------------------------------------------

#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{

  photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
  jetCollectionT_ = consumes<JetCollection>(iConfig.getParameter<edm::InputTag>("jetCollection"));
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  ESRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedESRecHitCollection"));
  AODEBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedAODEBRecHitCollection"));
  AODEERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedAODEERecHitCollection"));
  AODESRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedAODESRecHitCollection"));
  RECOEBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  RECOEERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  RECOESRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ESRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  genJetCollectionT_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetCollection"));
  trackCollectionT_ = consumes<pat::IsolatedTrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  rhoLabel_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));
  trgResultsT_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trgResults"));
  genInfoT_ = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"));
  lheEventT_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lhe"));
  pfCollectionT_ = consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCollection"));
  vertexCollectionT_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  puCollection_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupCollection"));

  z0PVCut_   = iConfig.getParameter<double>("z0PVCut");

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Output Tree
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("eventId", &eventId_);
  RHTree->Branch("runId",   &runId_);
  RHTree->Branch("lumiId",  &lumiId_);

  RHTree->Branch("SC_iphi", &vIphi_Emax_);
  RHTree->Branch("SC_ieta", &vIeta_Emax_);

  branchesDiPhotonSel ( RHTree, fs );
  branchesPhoObjSel   ( RHTree, fs );
  //branchesPhoObjSel_AA   ( RHTree, fs );
  branchesPhoVars     ( RHTree, fs );
  branchesEvtWgt      ( RHTree, fs );
  branchesPileup      ( RHTree, fs );
  branchesRecHits     ( RHTree, fs );
  branchesPFTracks    ( RHTree, fs );

  hNpassed_img = fs->make<TH1F>("hNpassed_img", "isPassed;isPassed;N", 2, 0., 2);
  hNRecoPho_noCut = fs->make<TH1F>("hNRecoPho_noCut", "NRecoPho_noCut;NRecoPho_noCut;N", 10, 0., 10);
  hNRecoPho_HLT = fs->make<TH1F>("hNRecoPho_HLT", "NRecoPho_HLT;NRecoPho_HLT;N", 10, 0., 10);
  hNRecoPho_presel = fs->make<TH1F>("hNRecoPho_presel", "NRecoPho_presel;NRecoPho_presel;N", 10, 0., 10);
  hNRecoPho_hits = fs->make<TH1F>("hNRecoPho_hits", "NRecoPho_hits;NRecoPho_hits;N", 10, 0., 10);
  hNRecoPho_genMatching = fs->make<TH1F>("hNRecoPho_genMatching", "NRecoPho_genMatching;NRecoPho_genMatching;N", 10, 0., 10);

  hLeadingPhoPt = fs->make<TH1F>("hLeadingPhoPt", "Leading_Photon_Pt;Leading_Photon_Pt;N", 1000, 0., 1000);
  hSubleadingPhoPt = fs->make<TH1F>("hSubleadingPhoPt", "Subleading_Photon_Pt;Subleading_Photon_Pt;N", 1000, 0., 1000);
}

RecHitAnalyzer::~RecHitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
//
// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  //edm::Handle<EcalRecHitCollection> EERecHitsH;
  //iEvent.getByToken(EERecHitCollectionT_, EERecHitsH);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  // Run explicit jet selection
  bool hasPassed;
  vPreselPhoIdxs_.clear();
  nTotal += nPhotons;
 
  std::cout << ">> nPhotons: " << nPhotons << std::endl;  // added by KYUNGMIN
  std::cout << ">> nTotal: " << nTotal << std::endl;      // added by KYUNGMIN
 
  hasPassed = runDiPhotonSel_Iso ( iEvent, iSetup ); // trigger + isolation preselections + N(reco pho) == 2
  //hasPassed = runDiPhotonSel ( iEvent, iSetup ); // trigger + N(reco pho) == 2
  std::cout << ">> vPreselPhoIdxs_.size(): " << vPreselPhoIdxs_.size() << std::endl;
  std::cout << ">> hasPassed: " << hasPassed << std::endl;      // added by KYUNGMIN

  hNRecoPho_noCut->Fill(vPhoNoCutIdxs_.size());
  hLeadingPhoPt->Fill(leadingPt);
  hSubleadingPhoPt->Fill(subleadingPt);

  if ( !hasPassed ) return;
 
  hNRecoPho_HLT->Fill(vTrigPhoIdxs_.size());
  hNRecoPho_presel->Fill(vPreselPhoIdxs_.size());

  nTrigPassed += vTrigPhoIdxs_.size();
  std::cout << ">> nTrigPassed: " << nTrigPassed << std::endl;
  nPreselPassed += vPreselPhoIdxs_.size();
  std::cout << ">> nPreselPassed: " << nPreselPassed << std::endl;
  
  // Get coordinates of photon supercluster seed
  // and check if in workable location in EB
  hNpassed_img->Fill(0.);
  nPho = 0;
  int iphi_Emax, ieta_Emax;
  float Emax;
  GlobalPoint pos_Emax;
  vIphi_Emax_.clear();
  vIeta_Emax_.clear();
  vRegressPhoIdxs_.clear();
  vPos_Emax_.clear();
  int iphi_, ieta_; // rows:ieta, cols:iphi

  for ( unsigned int iP : vPreselPhoIdxs_ ) {

    // KYUNGMIN
    std::cout << "For loop for vPreselPhoIdxs_: " << iP << std::endl;

    PhotonRef iPho( photons, iP );

    // Get underlying super cluster
    reco::SuperClusterRef const& iSC = iPho->superCluster();
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );
    if ( debug ) std::cout << " >> SChits.size: " << SCHits.size() << std::endl;

    // Get Emax crystal
    Emax = 0.;
    iphi_Emax = -1;
    ieta_Emax = -1;

    //KYUNGMIN
    std::cout << "SCHits.size: " << SCHits.size() << std::endl;

    // Loop over SC hits of photon
    for(unsigned iH(0); iH != SCHits.size(); ++iH) {

      // Get DetId
      //if ( SCHits[iH].first.subdetId() != EcalBarrel ) continue;
      if ( SCHits[iH].first.subdetId() != EcalBarrel ) {
        std::cout << iH << ", HITS NOT IN BARREL" << std::endl;
        continue;
      }
      else {
        std::cout << iH << ", HITS FOUND IN BARREL" << std::endl;
      }
      EcalRecHitCollection::const_iterator iRHit( EBRecHitsH->find(SCHits[iH].first) );
      if ( iRHit == EBRecHitsH->end() ) {
        std::cout << iH << ", EBRecHitsH->end(), PASS" << std::endl;
        continue;
      }
      //if ( iRHit == EBRecHitsH->end() ) continue;

      // Convert coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      //EBDetId ebId( iSC->seed()->seed() );
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,1,...,85]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]
      iphi_ = ebId.iphi()-1; // [0,...,359]


      // KYUNGMIN
      std::cout << "iRHit->energy: " << iRHit->energy() << ", Emax: " << Emax << std::endl;
      std::cout << "ieta: " << ieta_ << ", iphi: " << iphi_ << std::endl;

      // Keep coordinates of shower max
      if ( iRHit->energy() > Emax ) {
        Emax = iRHit->energy();
        iphi_Emax = iphi_;
        ieta_Emax = ieta_;
        pos_Emax = caloGeom->getPosition(ebId);
      }
      if ( debug ) std::cout << " >> " << iH << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << iRHit->energy() << std::endl;
    } // SC hits

    // Apply selection on position of shower seed
    if ( Emax <= zs ) {
      std::cout << iP << ", EMAX LESS THAN zs = " << zs << std::endl;
      continue;
    }
    else {
      std::cout << iP << ", EMAX GREATER THAN zs = " << zs << std::endl;
    }
    //if ( Emax <= zs ) continue;
    if ( ieta_Emax > 169 - 16 || ieta_Emax < 15 ) {
      std::cout << "SEED NOT CENTERED AROUND [15, 15]" << std::endl;
      continue; // seed centered on [15,15] so must be padded by 15 below and 16 above
    }
    //if ( ieta_Emax > 169 - 16 || ieta_Emax < 15 ) continue; // seed centered on [15,15] so must be padded by 15 below and 16 above
    vIphi_Emax_.push_back( iphi_Emax );
    vIeta_Emax_.push_back( ieta_Emax );
    vPos_Emax_.push_back( pos_Emax );
    vRegressPhoIdxs_.push_back( iP );
    if ( debug ) std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    nPho++;

    //KYUNGMIN
    std::cout << "nPho: " << nPho << std::endl;

  } // Photons
 
 
  // Enforce selection
  if ( debug ) std::cout << " >> nPho: " << nPho << std::endl;
  if ( nPho == 0 ) return; // pho obj selection - ignore event if no photon objects
  std::cout << " >> nPho: " << nPho << std::endl;
  if ( debug ) std::cout << " >> Passed cropping. " << std::endl;


  //nPassed += nPho;   // commented out by KYUNGMIN
  hNRecoPho_hits->Fill(nPho);

  //runPhoObjSel_AA ( iEvent, iSetup );
  runPhoObjSel_fake ( iEvent, iSetup );
  //runPhoObjSel_prompt ( iEvent, iSetup );
  nPhoObjPassed += vPhoObj_recoIdx_.size();

  std::cout << "vPhoObj_recoIdx_ Size: " << vPhoObj_recoIdx_.size() << std::endl;
  std::cout << "vRegressPhoIdxs_ Size: " << vRegressPhoIdxs_.size() << std::endl;
  

  if ( debug ) std::cout << "PhoObj Size: " << vPhoObj_recoIdx_.size() << std::endl;
  if ( vPhoObj_recoIdx_.size() == 0 ) return;
  
  hNRecoPho_genMatching->Fill(vPhoObj_recoIdx_.size());

  fillDiPhotonSel   ( iEvent, iSetup );
  fillPhoObjSel     ( iEvent, iSetup ); // right now does nothing
  fillPhoVars       ( iEvent, iSetup );
  fillEvtWgt        ( iEvent, iSetup );
  fillPileup        ( iEvent, iSetup );
  fillRecHits       ( iEvent, iSetup );
  fillPFTracks      ( iEvent, iSetup );
 
  nPassed += nPho;   // commented out by KYUNGMIN
  //nPassed += nPhoObjPassed;   // commented out by KYUNGMIN

  // KYUNGMIN
//  if ( vPreselPhoIdxs_.size() == 0 ) return;
//  hNpassed_img->Fill(0.);
//  nPassed = nPreselPassed;  // KYUNGMIN

  RHTree->Fill();
  hNpassed_img->Fill(1.);

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
RecHitAnalyzer::beginJob()
{
  nTotal = 0;
  nTrigPassed = 0;
  nPreselPassed = 0;
  nPhoObjPassed = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitAnalyzer::endJob()
{
  std::cout << ">> After Trigger: " << nTrigPassed << std::endl;
  std::cout << ">> After Egamma Iso Cuts: " << nPreselPassed << "/" << nTotal << std::endl;
  std::cout << ">> Gen-matching etc: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
