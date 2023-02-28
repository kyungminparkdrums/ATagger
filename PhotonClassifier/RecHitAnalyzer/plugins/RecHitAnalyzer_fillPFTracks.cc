#include "PhotonClassifier/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill PFTracks in EB 
// Store PFTracks in EB image

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesPFTracks ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("PFTimg_pT",      &vPFTimg_pT_);
  tree->Branch("PFTimg_QpT",     &vPFTimg_QpT_);
  tree->Branch("PFTimg_pT_PV",   &vPFTimg_pT_PV_);
  tree->Branch("PFTimg_pT_nPV",  &vPFTimg_pT_nPV_);

  tree->Branch("PFTimg_d0",      &vPFTimg_d0_);
  tree->Branch("PFTimg_z0",      &vPFTimg_z0_);

  tree->Branch("PFTgraph_nodes_img",      &vPFTgraph_nodes_img_); 
  tree->Branch("PFTgraph_nodes_r04",      &vPFTgraph_nodes_r04_);

} // branchesEB()

// Fill TRK rechits at EB/EE ______________________________________________________________//
void RecHitAnalyzer::fillPFTracks ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int idx_; // rows:ieta, cols:iphi
  float eta, phi;
  GlobalPoint pos;

  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::VertexCollection> vertexInfo;
  iEvent.getByToken(vertexCollectionT_, vertexInfo);
  const reco::VertexCollection& vtxs = *vertexInfo;

  vPFTimg_pT_.clear();
  vPFTimg_pT_PV_.clear();
  vPFTimg_pT_nPV_.clear();

  vPFTimg_d0_.clear();
  vPFTimg_z0_.clear();


  std::vector<float> PFTrack;
  std::vector<float> PFTrackPt;
  std::vector<float> PFTrackQPt;
  std::vector<float> PFTrackPt_PV;
  std::vector<float> PFTrackPt_nPV;
  
  std::vector<float> PFTrackd0;
  std::vector<float> PFTrackz0;

  int iphi_, ieta_;
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  //float x, y, z, r;
  
  // &&& MAKE IMAGE &&&
  for ( unsigned int iP(0); iP < vRegressPhoIdxs_.size(); iP++ ) {
    if ( std::find(vPhoObj_recoIdx_.begin(), vPhoObj_recoIdx_.end(), vRegressPhoIdxs_[iP]) == vPhoObj_recoIdx_.end() ) {
      continue;
    }

    PFTrack.assign(crop_size*crop_size,0.);
    PFTrackPt.assign(crop_size*crop_size,0.);
    PFTrackQPt.assign(crop_size*crop_size,0.);
    PFTrackPt_PV.assign(crop_size*crop_size,0.);
    PFTrackPt_nPV.assign(crop_size*crop_size,0.);

    PFTrackd0.assign(crop_size*crop_size,0.);
    PFTrackz0.assign(crop_size*crop_size,0.);

    iphi_shift = vIphi_Emax_[iP] - 15;
    ieta_shift = vIeta_Emax_[iP] - 15;
    if ( debug ) std::cout << " >> Doing pfcs: iphi_Emax,ieta_Emax: " << vIphi_Emax_[iP] << ", " << vIeta_Emax_[iP] << std::endl;

    for ( PFCollection::const_iterator iPFC = pfCandsH_->begin();
	  iPFC != pfCandsH_->end(); ++iPFC ) {

      const reco::Track* thisTrk = iPFC->bestTrack();
      if(!thisTrk) continue;

      eta = iPFC->eta();
      phi = iPFC->phi();
      if ( std::abs(eta) > 1.5 ) continue;
      
      float d0    =  ( !vtxs.empty() ? thisTrk->dxy(vtxs[0].position()) : thisTrk->dxy() );
      float z0    =  ( !vtxs.empty() ? thisTrk->dz(vtxs[0].position())  : thisTrk->dz() );

      DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
      if ( id.subdetId() != EcalBarrel ) continue;

      float thisTrkPt = thisTrk->pt();
      float thisTrkQPt = (thisTrk->pt()*thisTrk->charge());
      
      const reco::Track thisTrack = *thisTrk;
      
      if ( debug ) std::cout << "track pt: " << thisTrack.pt() << std::endl;
      if ( debug ) std::cout << "track n good hits: " << thisTrack.numberOfValidHits() << std::endl;
      if ( debug ) std::cout << "track dz: " << thisTrack.dz() << std::endl;
      if ( debug ) std::cout << "track dxy: " << thisTrack.dxy() << std::endl;
      if ( debug ) std::cout << "track vx: " << thisTrack.vx() << std::endl;
      if ( debug ) std::cout << "track vy: " << thisTrack.vy() << std::endl;
      if ( debug ) std::cout << "track vz: " << thisTrack.vz() << std::endl;
      
      // Convert detector coordinates to ordinals                                                                                                                                                         
      EBDetId ebId( id );
      iphi_ = ebId.iphi()-1; 
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]                                                                                                                     
      ieta_ += EBDetId::MAX_IETA; // [0,...,169}                                                                                                                                                      

      ieta_crop = ieta_ - ieta_shift;
      iphi_crop = iphi_ - iphi_shift;
      if ( iphi_crop >= EBDetId::MAX_IPHI ) iphi_crop = iphi_crop - EBDetId::MAX_IPHI; // get wrap-around hits                                                                                            
      if ( iphi_crop < 0 ) iphi_crop = iphi_crop + EBDetId::MAX_IPHI; // get wrap-around hits                                                                                                             

      if ( ieta_crop < 0 || ieta_crop > crop_size-1 ) continue;
      if ( iphi_crop < 0 || iphi_crop > crop_size-1 ) continue;

      idx_ = ieta_crop*crop_size + iphi_crop;
      
      PFTrack[idx_] += 1;
      PFTrackPt[idx_] += thisTrkPt;
      PFTrackQPt[idx_] += thisTrkQPt;

      PFTrackd0[idx_] += abs(d0);
      PFTrackz0[idx_] += abs(z0);

      if(fabs(z0) < z0PVCut_){
	PFTrackPt_PV[idx_] += thisTrkPt;
      }else{
	PFTrackPt_nPV[idx_] += thisTrkPt;
      }
      if ( debug ) std::cout << iP << " pfcs " << idx_ << " " << thisTrkPt << " " << d0 << " " << z0 << std::endl;
      
      if (isinf(thisTrkPt)){
	std::cout << "Track pt inf!!" << std::endl;
      }
      if (isnan(d0)){
	std::cout << "d0 nan!!" << std::endl;
      }
      if (isnan(z0)){
	std::cout << "z0 nan!!" << std::endl;
      }

    }

    for (int i = 0; i < crop_size*crop_size; i++) {
    
      if (PFTrackd0[i] > 0.) PFTrackd0[i] = PFTrackd0[i] / PFTrack[i];
      if (PFTrackz0[i] > 0.) PFTrackz0[i] = PFTrackz0[i] / PFTrack[i];
    }
    
    vPFTimg_pT_.push_back( PFTrackPt );
    vPFTimg_QpT_.push_back( PFTrackQPt );
    vPFTimg_pT_PV_.push_back( PFTrackPt_PV );
    vPFTimg_pT_nPV_.push_back( PFTrackPt_nPV );
    
    vPFTimg_d0_.push_back( PFTrackd0 );
    vPFTimg_z0_.push_back( PFTrackz0 );

  }

  // %%% MAKE GRAPH %%%
  vPFTgraph_nodes_img_.clear();
  vPFTgraph_nodes_r04_.clear();
  std::vector<float> PFTgraph_nodes_img;
  std::vector<float> PFTgraph_nodes_r02;
  std::vector<float> PFTgraph_nodes_r04;
  std::vector<float> PFTgraph_nodes_r06;
  std::vector<float> pftrack;
  float radius, phi_shift, eta_shift, phi_rel, eta_rel, iphi_rel, ieta_rel;
  for ( unsigned int iP(0); iP < vRegressPhoIdxs_.size(); iP++ ) {
    if ( std::find(vPhoObj_recoIdx_.begin(), vPhoObj_recoIdx_.end(), vRegressPhoIdxs_[iP]) == vPhoObj_recoIdx_.end() ) {
      continue;
    }
    
    PFTgraph_nodes_img.clear();
    PFTgraph_nodes_r02.clear();
    PFTgraph_nodes_r04.clear();
    PFTgraph_nodes_r06.clear();

    iphi_shift = vIphi_Emax_[iP] - 15;
    ieta_shift = vIeta_Emax_[iP] - 15;
    phi_shift  = vPos_Emax_[iP].phi();
    eta_shift  = vPos_Emax_[iP].eta();
    for ( PFCollection::const_iterator iPFC = pfCandsH_->begin();
          iPFC != pfCandsH_->end(); ++iPFC ) {

      const reco::Track* thisTrk = iPFC->bestTrack();
      if(!thisTrk) continue;

      eta = iPFC->eta();
      phi = iPFC->phi();
      if ( std::abs(eta) > 1.5 ) continue;

      float d0    =  ( !vtxs.empty() ? thisTrk->dxy(vtxs[0].position()) : thisTrk->dxy() );
      float z0    =  ( !vtxs.empty() ? thisTrk->dz(vtxs[0].position())  : thisTrk->dz() );

      DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
      if ( id.subdetId() != EcalBarrel ) continue;

      float thisTrkPt = thisTrk->pt();
      float thisTrkQPt = (thisTrk->pt()*thisTrk->charge());

      const reco::Track thisTrack = *thisTrk;

      phi_rel = phi - phi_shift;
      eta_rel = eta - eta_shift;

      if (abs(phi_rel + 2*M_PI) < abs(phi_rel)) {
        phi_rel = phi_rel + 2*M_PI;
      } else if (abs(phi_rel - 2*M_PI) < abs(phi_rel)) {
        phi_rel = phi_rel - 2*M_PI;
      }

      radius = sqrt( pow(phi_rel, 2.) + pow(eta_rel, 2.) );
      if ( radius > 0.8 ) continue;

      EBDetId ebId( id );
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
      pftrack.clear();
      pftrack.push_back(thisTrkPt);
      pftrack.push_back(float(phi_rel));
      pftrack.push_back(float(eta_rel));
      pftrack.push_back(float(d0));
      pftrack.push_back(float(z0));
      pftrack.push_back(thisTrk->charge());

      if ( ieta_crop >= 0 && ieta_crop <= crop_size-1 && iphi_crop >= 0 && iphi_crop <= crop_size-1 ) {
        PFTgraph_nodes_img.insert( PFTgraph_nodes_img.end(), pftrack.begin(), pftrack.end() ); 
      }
      if ( radius <= 0.2 ) PFTgraph_nodes_r02.insert( PFTgraph_nodes_r02.end(), pftrack.begin(), pftrack.end() );
      if ( radius <= 0.4 ) PFTgraph_nodes_r04.insert( PFTgraph_nodes_r04.end(), pftrack.begin(), pftrack.end() );
      if ( radius <= 0.6 ) PFTgraph_nodes_r06.insert( PFTgraph_nodes_r06.end(), pftrack.begin(), pftrack.end() );
      if ( debug ) std::cout << "Track pt " << thisTrkPt << " radius " << radius << " hit len " << pftrack.size()
			     << " phi_rel, eta_rel " << phi_rel << " " << eta_rel << std::endl;
    }

    std::cout << "pft sizes fct 6 img, r=0.2, r=0.4, r=0.6: " << PFTgraph_nodes_img.size() <<", "<< PFTgraph_nodes_r02.size() << ", "
              << PFTgraph_nodes_r04.size() << ", " << PFTgraph_nodes_r06.size() << std::endl;
    
    vPFTgraph_nodes_img_.push_back( PFTgraph_nodes_img );
    vPFTgraph_nodes_r04_.push_back( PFTgraph_nodes_r04 );

  } // photons graph

} // fillPFTracks()
