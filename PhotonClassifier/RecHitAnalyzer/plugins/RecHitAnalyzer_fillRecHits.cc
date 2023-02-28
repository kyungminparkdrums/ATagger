#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <math.h>

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesRecHits ( TTree* tree, edm::Service<TFileService> &fs )
{
  hRHimg_energy = fs->make<TProfile2D>("SC_energy", "E(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );
  hRHimg_time = fs->make<TProfile2D>("SC_time", "t(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );

  RHTree->Branch("RHimg_energy",   &vRHimg_energy_);
  RHTree->Branch("RHimg_energyT",  &vRHimg_energyT_);
  RHTree->Branch("RHimg_energyZ",  &vRHimg_energyZ_);
  RHTree->Branch("RHimg_time",     &vRHimg_time_);
  
  RHTree->Branch("RHgraph_nodes_img",  &vRHgraph_nodes_img_);
  RHTree->Branch("RHgraph_nodes_r04",  &vRHgraph_nodes_r04_);

} // branchesSC()

// Fill SC rechits _________________________________________________________________//
void RecHitAnalyzer::fillRecHits ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  int idx_; // rows:ieta, cols:iphi                                                                                                                                                                      
  float eta, phi;
  GlobalPoint pos;

  std::cout << "Regress: ";
  for ( int iP : vRegressPhoIdxs_ ) {
    std::cout << iP << " ";
  }
  std::cout << "Reco: ";
  for ( int iP : vPhoObj_recoIdx_ ) {
    std::cout << iP << " ";
  }
  std::cout << std::endl;

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();


  // %%% MAKE IMAGES %%%
  vRHimg_energy_.clear();
  vRHimg_energyT_.clear();
  vRHimg_energyZ_.clear();
  vRHimg_time_.clear();
  std::vector<float> RHimg_energy;
  std::vector<float> RHimg_energyT;
  std::vector<float> RHimg_energyZ;
  std::vector<float> RHimg_time;

  int iphi_, ieta_;
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  for ( unsigned int iP(0); iP < vRegressPhoIdxs_.size(); iP++ ) {
  //for ( int iP : vRegressPhoIdxs_ ) {
 
    std::cout << "COMPARE vRegressPhoIdxs and vRegressPhoIdxs" << std::endl;
    if ( std::find(vPhoObj_recoIdx_.begin(), vPhoObj_recoIdx_.end(), vRegressPhoIdxs_[iP]) == vPhoObj_recoIdx_.end() ) {
      continue;
    }
 
    std::cout << "GOING THROUGH vRegressPhoIdxs_" << std::endl;
 
    if ( debug ) std::cout << "Img Index " << vRegressPhoIdxs_[iP] << std::endl;
    RHimg_energy.assign(crop_size*crop_size,0.);
    RHimg_energyT.assign(crop_size*crop_size,0.);
    RHimg_energyZ.assign(crop_size*crop_size,0.);
    RHimg_time.assign(crop_size*crop_size,0.);

    iphi_shift = vIphi_Emax_[iP] - 15;
    ieta_shift = vIeta_Emax_[iP] - 15;
    if ( debug ) std::cout << " >> Doing pho img: iphi_Emax,ieta_Emax: " << vIphi_Emax_[iP] << ", " << vIeta_Emax_[iP] << std::endl;

    for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
        iRHit != EBRecHitsH->end();
        ++iRHit) {

      if ( iRHit->energy() < zs ) continue;

      // Convert detector coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      iphi_ = ebId.iphi()-1; // [0,...,359]
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]

      // Convert to [0,...,31][0,...,31]
      ieta_crop = ieta_ - ieta_shift;
      iphi_crop = iphi_ - iphi_shift;
      if ( iphi_crop >= EBDetId::MAX_IPHI ) iphi_crop = iphi_crop - EBDetId::MAX_IPHI; // get wrap-around hits
      if ( iphi_crop < 0 ) iphi_crop = iphi_crop + EBDetId::MAX_IPHI; // get wrap-around hits

      if ( ieta_crop < 0 || ieta_crop > crop_size-1 ) continue;
      if ( iphi_crop < 0 || iphi_crop > crop_size-1 ) continue;

      // Convert to [0,...,32*32-1]
      idx_ = ieta_crop*crop_size + iphi_crop;

      auto pos = caloGeom->getPosition(ebId);

      // Fill branch arrays
      RHimg_energy[idx_] = iRHit->energy();
      RHimg_energyT[idx_] = iRHit->energy()/TMath::CosH(pos.eta());
      RHimg_energyZ[idx_] = iRHit->energy()*std::abs(TMath::TanH(pos.eta()));
      RHimg_time[idx_] = iRHit->time();

      if ( debug ) std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 
      if ( debug ) std::cout << "idx,ieta,iphi,E:" <<idx_<<","<< ieta_crop << "," << iphi_crop << "," << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      std::cout << "iRHit->energy(): " << iRHit->energy() << " at iphi_crop: " << iphi_crop << std::endl;
      hRHimg_energy->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hRHimg_time->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // EB rechits

    vRHimg_energy_.push_back( RHimg_energy );
    vRHimg_energyT_.push_back( RHimg_energyT );
    vRHimg_energyZ_.push_back( RHimg_energyZ );
    vRHimg_time_.push_back( RHimg_time );

  } // photons image

  // %%% MAKE GRAPH %%%
  vRHgraph_nodes_img_.clear();
  vRHgraph_nodes_r04_.clear();
  std::vector<float> RHgraph_nodes_img;
  std::vector<float> RHgraph_nodes_r02;
  std::vector<float> RHgraph_nodes_r04;
  std::vector<float> RHgraph_nodes_r06;
  std::vector<float> rechit;
  float radius, phi_shift, phi_rel, eta_shift, eta_rel;
  float iphi_rel, ieta_rel;
  for ( unsigned int iP(0); iP < vRegressPhoIdxs_.size(); iP++ ) {
    if ( std::find(vPhoObj_recoIdx_.begin(), vPhoObj_recoIdx_.end(), vRegressPhoIdxs_[iP]) == vPhoObj_recoIdx_.end() ) {
      continue;
    }

    if ( debug ) std::cout << "Graph Index " << vRegressPhoIdxs_[iP] << std::endl;

    RHgraph_nodes_img.clear();
    RHgraph_nodes_r02.clear();
    RHgraph_nodes_r04.clear();
    RHgraph_nodes_r06.clear();
    
    iphi_shift = vIphi_Emax_[iP] - 15;
    ieta_shift = vIeta_Emax_[iP] - 15;
    phi_shift  = vPos_Emax_[iP].phi();
    eta_shift  = vPos_Emax_[iP].eta();
    if ( debug ) std::cout << " >> Doing pho graph: iphi_Emax,ieta_Emax: " << vIphi_Emax_[iP] << ", " << vIeta_Emax_[iP] << std::endl;

    for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
        iRHit != EBRecHitsH->end();
        ++iRHit) {

      if ( iRHit->energy() < zs ) continue;

      // Convert detector coordinates to ordinals                                                                                                                                                        
      EBDetId ebId( iRHit->id() );
      pos = caloGeom->getPosition( ebId ); 
      phi = pos.phi();
      eta = pos.eta();
      
      phi_rel = phi - phi_shift;
      eta_rel = eta - eta_shift;

      if (abs(phi_rel + 2*M_PI) < abs(phi_rel)) {
	phi_rel = phi_rel + 2*M_PI;
      } else if (abs(phi_rel - 2*M_PI) < abs(phi_rel)) {
	phi_rel = phi_rel - 2*M_PI;
      }

      radius = sqrt( pow(phi_rel, 2.) + pow(eta_rel, 2.) );
      if ( radius > 0.8 ) continue;

      iphi_ = ebId.iphi()-1; // [0,...,359]       
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]  
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]                                     

      iphi_rel = iphi_ - vIphi_Emax_[iP];
      ieta_rel = ieta_ - vIeta_Emax_[iP];

      if (abs(iphi_rel + EBDetId::MAX_IPHI) < abs(iphi_rel)) {
        iphi_rel = iphi_rel + EBDetId::MAX_IPHI;
      } else if (abs(iphi_rel - EBDetId::MAX_IPHI) < abs(iphi_rel)) {
        iphi_rel = iphi_rel - EBDetId::MAX_IPHI;
      }

      ieta_crop = ieta_rel + 15;
      iphi_crop = iphi_rel + 15;

      // Fill branch arrays                                                                                                                                                                              
      rechit.clear();
      rechit.push_back(iRHit->energy());
      rechit.push_back(float(phi_rel));
      rechit.push_back(float(eta_rel));

      if ( ieta_crop >= 0 && ieta_crop <= crop_size-1 && iphi_crop >= 0 && iphi_crop <= crop_size-1 ) {
	RHgraph_nodes_img.insert( RHgraph_nodes_img.end(), rechit.begin(), rechit.end() );
      }            
      if ( radius <= 0.2 ) RHgraph_nodes_r02.insert(  RHgraph_nodes_r02.end(), rechit.begin(), rechit.end() );
      if ( radius <= 0.4 ) RHgraph_nodes_r04.insert(  RHgraph_nodes_r04.end(), rechit.begin(), rechit.end() );
      if ( radius <= 0.6 ) RHgraph_nodes_r06.insert(  RHgraph_nodes_r06.end(), rechit.begin(), rechit.end() );
      if ( debug ) std::cout << "RecHit energy " << iRHit->energy() << " radius " << radius << " hit len " << rechit.size() 
			     << " phi rel, eta rel " << phi_rel << " " << eta_rel << std::endl; 

    } // rec hits
    
    std::cout << "rh sizes fct 3 img, r=0.2, r=0.4, r=0.6: " << RHgraph_nodes_img.size() <<", "<< RHgraph_nodes_r02.size() << ", "
	      << RHgraph_nodes_r04.size() << ", " << RHgraph_nodes_r06.size() << std::endl;
  
    vRHgraph_nodes_img_.push_back( RHgraph_nodes_img );
    vRHgraph_nodes_r04_.push_back( RHgraph_nodes_r04 );

  } // photons graph

} // fillRecHits()
