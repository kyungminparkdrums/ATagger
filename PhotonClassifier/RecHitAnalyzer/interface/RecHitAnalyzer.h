#ifndef RecHitAnalyzer_h
#define RecHitAnalyzer_h

//---------------------------------------------                                                                                                                                                          
//                                                                                                                                                                                                       
// Package: PhotonClassifier                                                                                                                                                                             
// Class: RecHitAnalyzer                                                                                                                                                                                 
//                                                                                                                                                                                                       
// Author: Andrew C. Roberts                                                                                                                                                                             
// Started 2022/5/18                                                                                                                                                                                     
// Last Updated 2022/5/18                                                                                                                                                                                
//                                                                                                                                                                                                       
//---------------------------------------------

#include <vector>
#include "TCanvas.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include "TProfile2D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//#ifdef __MAKECINT__ 
//#pragma link C++ class std::vector<std::vector<float>>+;
//#pragma link C++ class std::vector<std::vector<std::vector<float>>>+;
//#endif

//edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >
//
// class declaration
//
//using pat::MuonCollection;
//using pat::MuonRef;
using pat::ElectronCollection;
using pat::ElectronRef;
using pat::JetCollection;
using pat::JetRef;
using pat::PhotonCollection;
using pat::PhotonRef;

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit RecHitAnalyzer(const edm::ParameterSet&);
    ~RecHitAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<ElectronCollection> electronCollectionT_;
    //edm::EDGetTokenT<MuonCollection> muonCollectionT_;
    edm::EDGetTokenT<PhotonCollection> photonCollectionT_;
    edm::EDGetTokenT<JetCollection> jetCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> ESRecHitCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;

    edm::EDGetTokenT<EcalRecHitCollection> AODEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> AODEERecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> AODESRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOEERecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOESRecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
    edm::EDGetTokenT<pat::IsolatedTrackCollection> trackCollectionT_;
    edm::EDGetTokenT<double> rhoLabel_;
    edm::EDGetTokenT<edm::TriggerResults> trgResultsT_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoT_;
    edm::EDGetTokenT<LHEEventProduct> lheEventT_;

    edm::EDGetTokenT<reco::VertexCollection> vertexCollectionT_;
    typedef std::vector<pat::PackedCandidate>  PFCollection;
    edm::EDGetTokenT<PFCollection> pfCollectionT_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puCollection_;
    
    //static const int nPhotons = 2;
    //static const int nPhotons = 1;

    TProfile2D *hEB_energy;
    TProfile2D *hEB_time;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_time_;

    TProfile2D * hRHimg_energy;
    TProfile2D * hRHimg_time;
    TH1F * hSC_mass;
    TH1F * hDR;
    TH1F * hdEta;
    TH1F * hdPhi;
    TH3F * hdPhidEtaM;
    TH2F * hdPhidEta;
    TH1F * hSC_pT;

    TTree* RHTree;

    unsigned int nPho;
    unsigned long long eventId_;
    unsigned int runId_;
    unsigned int lumiId_;
    double z0PVCut_;

    // selection - option to add different selection scripts
    void branchesDiPhotonSel ( TTree*, edm::Service<TFileService>& );
    bool runDiPhotonSel ( const edm::Event&, const edm::EventSetup& );
    bool runDiPhotonSel_Iso ( const edm::Event&, const edm::EventSetup& );
    void fillDiPhotonSel ( const edm::Event&, const edm::EventSetup& );
    
    void branchesPhoObjSel_AA ( TTree*, edm::Service<TFileService>& );
    void branchesPhoObjSel ( TTree*, edm::Service<TFileService>& );
    bool runPhoObjSel_AA ( const edm::Event&, const edm::EventSetup& );
    bool runPhoObjSel_fake ( const edm::Event&, const edm::EventSetup& );
    bool runPhoObjSel_prompt ( const edm::Event&, const edm::EventSetup& );
    void fillPhoObjSel ( const edm::Event&, const edm::EventSetup& );

    // general photon object variables
    void branchesPhoVars ( TTree*, edm::Service<TFileService>& );
    void fillPhoVars ( const edm::Event&, const edm::EventSetup& );

    // event weight
    void branchesEvtWgt ( TTree*, edm::Service<TFileService>& );
    void fillEvtWgt ( const edm::Event&, const edm::EventSetup& );

    // pileup
    void branchesPileup ( TTree*, edm::Service<TFileService>& );
    void fillPileup ( const edm::Event&, const edm::EventSetup& );

    // ecal rec hits
    void branchesRecHits ( TTree*, edm::Service<TFileService>& );
    void fillRecHits ( const edm::Event&, const edm::EventSetup& );

    // particle flow tracks - tracks associated with particle flow candidates
    void branchesPFTracks ( TTree*, edm::Service<TFileService>& );
    void fillPFTracks ( const edm::Event&, const edm::EventSetup& );

    std::map<unsigned int, std::vector<unsigned int>> mGenPi0_RecoPho;

    std::vector<unsigned int> vRecoPhoIdxs_;
    std::vector<unsigned int> vTrigPhoIdxs_;
    std::vector<unsigned int> vPreselPhoIdxs_;
    std::vector<unsigned int> vRegressPhoIdxs_;
    std::vector<float> vIphi_Emax_;
    std::vector<float> vIeta_Emax_;
    std::vector<GlobalPoint> vPos_Emax_;

    std::vector<int> vPhoNoCutIdxs_;
    float leadingPt;
    float subleadingPt;

    // RecHitAnalyzer_fillRecHits.cc
    std::vector<std::vector<float>> vRHimg_energy_;
    std::vector<std::vector<float>> vRHimg_energyT_;
    std::vector<std::vector<float>> vRHimg_energyZ_;
    std::vector<std::vector<float>> vRHimg_time_;
    
    std::vector<std::vector<float>> vRHgraph_nodes_img_;
    std::vector<std::vector<float>> vRHgraph_nodes_r04_;

    // RecHitAnalyzer_fillPFTracks.cc
    std::vector<std::vector<float>> vPFTimg_pT_;
    std::vector<std::vector<float>> vPFTimg_QpT_;
    std::vector<std::vector<float>> vPFTimg_pT_PV_;
    std::vector<std::vector<float>> vPFTimg_pT_nPV_;

    std::vector<std::vector<float>> vPFTimg_d0_;
    std::vector<std::vector<float>> vPFTimg_z0_;
    
    std::vector<std::vector<float>> vPFTgraph_nodes_img_;
    std::vector<std::vector<float>> vPFTgraph_nodes_r04_;

    TH2F *hTracks_EB;
    TH2F *hTracksPt_EB;
    
    // RecHitAnalyzer_fillPhoVars.cc
    std::vector<float> vPho_pT_;
    std::vector<float> vPho_E_;
    std::vector<float> vPho_eta_;
    std::vector<float> vPho_phi_;
    std::vector<float> vPho_ecalEPostCorr_;

    std::vector<float> vPho_r9_;
    std::vector<float> vPho_sieie_;
    std::vector<float> vPho_phoIso_;
    std::vector<float> vPho_chgIso_;
    std::vector<float> vPho_chgIsoWrongVtx_;
    std::vector<float> vPho_Eraw_;
    std::vector<float> vPho_phiWidth_;
    std::vector<float> vPho_etaWidth_;
    std::vector<float> vPho_scEta_;
    std::vector<float> vPho_sieip_;
    std::vector<float> vPho_s4_;
    std::vector<float> vPho_rho_;

    std::vector<float> vPho_neuIso_;
    std::vector<float> vPho_ecalIso_;
    std::vector<float> vPho_trkIso_;
    std::vector<float> vPho_hasPxlSeed_;
    std::vector<float> vPho_passEleVeto_;
    std::vector<float> vPho_HoE_;
    std::vector<float> vPho_phoIsoCorr_;
    std::vector<float> vPho_ecalIsoCorr_;
    std::vector<float> vPho_neuIsoCorr_;
    std::vector<float> vPho_chgIsoCorr_;
    std::vector<float> vPho_bdt_;

    std::vector<float> vPhoObj_recoIdx_;
    std::vector<float> vPhoObj_genMatchedIdx_;
    std::vector<float> vPhoObj_genMatchedPdgId_;
    std::vector<float> vPhoObj_ancestorPdgId_;
   
    std::vector<int> vGenStatus;

    std::vector<float> vA_mass_;
    std::vector<float> vA_DR_;
    std::vector<float> vA_E_;
    std::vector<float> vA_pT_;
    std::vector<float> vA_eta_;
    std::vector<float> vA_phi_;
    double mHgen_;
 
    int nTotal, nTrigPassed, nPreselPassed, nPhoObjPassed, nPassed;
    TH1F * hNpassed_kin;
    TH1F * hNpassed_presel;
    TH1F * hNpassed_mGG;
    TH1F * hNpassed_nRecoPho;
    TH1F * hNpassed_hlt;
    TH1F * hNpassed_img;
    TH1F * hNRecoPho_reco;
    TH1F * hNRecoPho_HLT;
    TH1F * hNRecoPho_presel;
    TH1F * hNRecoPho_hits;
    TH1F * hNRecoPho_genMatching;
    TH1F * hLeadingPhoPt;
    TH1F * hSubleadingPhoPt;

    float m0_;
    std::vector<float> vFC_inputs_;
    int hltAccept_;
    unsigned int nRecoPho_;
    std::vector<float> vMinDR_;
    double evtWeight_;

};

//
// constants, enums and typedefs
//
static const float zs = 0.;
static const int nEE = 2;

static const int crop_size = 32;
//static const bool debug = true;
static const bool debug = false;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100

static const int HBHE_IETA_MAX_FINE = 20;
static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;                                                                                                                                
static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1                                                                                                                                  
static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;                                                                                                                                
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17                                                                                                                           
static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;                                                                                                                                      
static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;                                                                                                                                       
static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;                                                                                                                                       
static const int ECAL_IETA_MAX_EXT = 140;

//
// static data member definitions
//

#endif
