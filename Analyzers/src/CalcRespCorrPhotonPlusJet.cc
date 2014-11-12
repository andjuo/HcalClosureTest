//
#include "HcalClosureTest/Analyzers/interface/CalcRespCorrPhotonPlusJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include <iostream>
#include <vector>
#include <set>
#include <map>
using namespace std;

#include <boost/regex.hpp>

// -------------------------------------------------

inline
unsigned int helper_findTrigger(const std::vector<std::string>& list,
				const std::string& name)
{
  boost::regex re(std::string("^(")+name+"|"+name+"_v\\d*)$");
  for (unsigned int i = 0,n = list.size() ; i < n ; ++i) {
    if(boost::regex_match(list[i],re)) return i;
  }
  return list.size();
}

// -------------------------------------------------

inline void HERE(const char *msg) { std::cout << msg << std::endl; }

template <class T>
inline void HERE(const char *msg, const T& x) {
  std::cout << msg << std::endl;
  std::cout << x << std::endl;
}

// -------------------------------------------------

void printElementsInBlocks(const PFCandidate& cand,
                           std::ostream& out=std::cout)  {
  if(!out) return;
  PFBlockRef firstRef;
  assert(!cand.elementsInBlocks().empty() );
  out << cand << "\n";
  for(unsigned i=0; i<cand.elementsInBlocks().size(); i++) {
    HERE("iBlock=",i);
    PFBlockRef blockRef = cand.elementsInBlocks()[i].first;
    if(blockRef.isNull()) {
      cerr<<"ERROR! no block ref!";
      continue;
    }
    if (0) {
      // original version
      if(!i) {
	out<<(*blockRef);
	firstRef = blockRef;
      }
      else if( blockRef!=firstRef) {
	cerr<<"WARNING! This PFCandidate is not made from a single block"<<endl;
      }
      out<<"\t"<<cand.elementsInBlocks()[i].second<<endl;
    }
    else {
      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
      for (unsigned int iEle=0; iEle<elements.size(); ++iEle) {
	//std::cout << "test " << iEle << ", " << elements[iEle].index() << " vs " << cand.elementsInBlocks()[i].second << "\n";
	if (elements[iEle].index() == cand.elementsInBlocks()[i].second) {
	  out << (*blockRef) << "\n";
	  //out << "found iEle=" << iEle << "\n";
	}
      }
    }
  }
  HERE("\tprint done\n\n");
}


// -------------------------------------------------

CalcRespCorrPhotonPlusJet::CalcRespCorrPhotonPlusJet(const edm::ParameterSet& iConfig)
{
  // set parameters
  debug_               = iConfig.getUntrackedParameter<int>("debug", 0);
  rhoCollection_       = iConfig.getParameter<edm::InputTag>("rhoColl");
  photonCollName_      = iConfig.getParameter<std::string>("photonCollName");
  caloJetCollName_     = iConfig.getParameter<std::string>("caloJetCollName");
  caloJetCorrName_     = iConfig.getParameter<std::string>("caloJetCorrName");
  pfJetCollName_       = iConfig.getParameter<std::string>("pfJetCollName");
  pfJetCorrName_       = iConfig.getParameter<std::string>("pfJetCorrName");
  genJetCollName_      = iConfig.getParameter<std::string>("genJetCollName");
  genParticleCollName_ = iConfig.getParameter<std::string>("genParticleCollName");
  genEventInfoName_    = iConfig.getParameter<std::string>("genEventInfoName");
  hbheRecHitName_      = iConfig.getParameter<std::string>("hbheRecHitName");
  hfRecHitName_        = iConfig.getParameter<std::string>("hfRecHitName");
  hoRecHitName_        = iConfig.getParameter<std::string>("hoRecHitName");
  rootHistFilename_    = iConfig.getParameter<std::string>("rootHistFilename");

  allowNoPhoton_       = iConfig.getParameter<bool>("allowNoPhoton");
  photonPtMin_         = iConfig.getParameter<double>("photonPtMin");
  photonJetDPhiMin_    = iConfig.getParameter<double>("photonJetDPhiMin");
  jetEtMin_            = iConfig.getParameter<double>("jetEtMin");
  jet2EtMax_           = iConfig.getParameter<double>("jet2EtMax");
  jet3EtMax_           = iConfig.getParameter<double>("jet3EtMax");
  photonTrigNamesV_    = iConfig.getParameter<std::vector<std::string>>("photonTriggers");
  jetTrigNamesV_       = iConfig.getParameter<std::vector<std::string>>("jetTriggers");
  writeTriggerPrescale_= iConfig.getParameter<bool>("writeTriggerPrescale");

  doCaloJets_          = iConfig.getParameter<bool>("doCaloJets");
  doPFJets_            = iConfig.getParameter<bool>("doPFJets");
  doGenJets_           = iConfig.getParameter<bool>("doGenJets");

  // set it here to ensure the value is defined
  eventWeight_ = 1.0;
}

CalcRespCorrPhotonPlusJet::~CalcRespCorrPhotonPlusJet()
{
  if (0) std::cout << CLHEP::electron_charge; // get rid of compiler error
}
  
//
// member functions
//
  
// ------------ method called to for each event  ------------
void CalcRespCorrPhotonPlusJet::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup)
{ 
  // Make preliminary checks to see that it is worth analyzing the event
  // 1. At least one photon with a goot pT
  // 2. At least one jet
  // 3. Trigger fired

  // 1st. Get Photons //
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByLabel(photonCollName_, photons);

  if(!photons.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound)
      << " could not find PhotonCollection named " << photonCollName_ << ".\n";
    return;
  }

  if ((photons->size()==0) && !allowNoPhoton_) {
    if (debug_) std::cout << "no photons in the event\n";
    return;
  }
  nPhotons_= photons->size();

  ///////// Run Over Photons /////////

  // sort photons by Et //
  // counter is needed later to get the reference to the ptr
  std::set<PhotonPair, PhotonPairComp> photonpairset;
  int counter=0;
  for(reco::PhotonCollection::const_iterator it=photons->begin(); it!=photons->end(); ++it) {
    const reco::Photon* photon=&(*it);
    photonpairset.insert( PhotonPair(photon, photon->pt(), counter) );
    counter++;
  }

  ///////////////////////////////
  // TAG = Highest Et photon
  ///////////////////////////////

  // find highest Et photon //
  PhotonPair  photon_tag;
  PhotonPair  photon_2nd;
  counter=0;
  for(std::set<PhotonPair, PhotonPairComp>::const_iterator it=photonpairset.begin(); it!=photonpairset.end(); ++it) {
    PhotonPair photon=(*it);
    ++counter;
    if(counter==1) photon_tag = photon;
    else if (counter==2) photon_2nd = photon;
    else break;
  }

  if(!photon_tag.photon() && !allowNoPhoton_)return; // should be unreachable

  // cut on photon pt
  if (photon_tag.isValid() && ( photon_tag.pt() < photonPtMin_ )) {
    if (debug_) std::cout << "largest photonPt=" << photon_tag.pt()<<std::endl;
    return;
  }

  // 2nd. Get Jets
  edm::Handle<reco::CaloJetCollection> calojets;
  edm::Handle<reco::PFJetCollection> pfjets;
  nCaloJets_=0;
  nPFJets_=0;
  nGenJets_=0;

  unsigned int anyJetCount=0;
  if (doCaloJets_) {
    iEvent.getByLabel(caloJetCollName_,calojets);
    if(!calojets.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find CaloJetCollection named " << caloJetCollName_ << ".\n";
      return;
    }
    anyJetCount+= calojets->size();
    nCaloJets_= calojets->size();
  }

  if (doPFJets_) {
    iEvent.getByLabel(pfJetCollName_,pfjets);
    if(!pfjets.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find PFJetCollection named " << pfJetCollName_ << ".\n";
      return;
    }
    anyJetCount+= pfjets->size();
    nPFJets_ = pfjets->size();
  }

  if (anyJetCount==0) {
    if (debug_) std::cout << "event contains no jets\n";
    return;
  }

  if (debug_) std::cout << "nPhotons=" << nPhotons_ << ", nCaloJets=" << nCaloJets_ << ", nPFJets=" << nPFJets_ << std::endl;

  // 3rd. Check the trigger
  photonTrigFired_.clear();
  photonTrigPrescale_.clear();
  jetTrigFired_.clear();
  jetTrigPrescale_.clear();

  // HLT Trigger
  //HLTConfigProvider ;


  // assign "trig fired" if no triggers are specified
  bool photonTrigFlag= false;
  bool jetTrigFlag= false;
  if ((photonTrigNamesV_.size()==1) &&
      (photonTrigNamesV_[0].length()==0)) photonTrigFlag=true;
  if ((jetTrigNamesV_.size()==1) &&
      (jetTrigNamesV_[0].length()==0)) jetTrigFlag=true;

  // If needed, process trigger information
  if (!photonTrigFlag || !jetTrigFlag) {
    // check the triggers

    edm::Handle<edm::TriggerResults> triggerResults;
    if( !iEvent.getByLabel(edm::InputTag("TriggerResults::HLT"),triggerResults) ) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find TriggerResults::HLT\n";
      return;
    }
    const edm::TriggerNames &evTrigNames =iEvent.triggerNames(*triggerResults);
    size_t id = 0;
    for (size_t i=0; i<photonTrigNamesV_.size(); ++i) {
      const std::string trigName=photonTrigNamesV_.at(i);
    //for (size_t i=0; i<evTrigNames.triggerNames().size(); i++) {
    //  std::string trigName= evTrigNames.triggerNames().at(i);
    //  if (trigName.find("_Photon")==std::string::npos) continue;
      id= helper_findTrigger(evTrigNames.triggerNames(),trigName);
      if (id==evTrigNames.size()) {
	photonTrigFired_.push_back(0);
	photonTrigPrescale_.push_back(-1);
	// we may want to debug the names of the used triggers
	//throw edm::Exception(edm::errors::ProductNotFound)
	//  << " could not find trigger " << trigName << "\n";
	//return;
	continue;
      }
      int fired= triggerResults->accept(id);
      if (fired) photonTrigFlag=true;
      //std::cout << "trigger " << trigName << " fired=" << fired << "\n";
      //std::cout << "identified as " << evTrigNames.triggerNames().at(id) << "\n";
      photonTrigFired_.push_back(fired);
      if (!writeTriggerPrescale_) photonTrigPrescale_.push_back(-1);
      else {
	// for triggers with two L1 seeds this fails
	std::pair<int,int> prescaleVals= hltConfig_.prescaleValues(iEvent,evSetup, evTrigNames.triggerName(id));
	photonTrigPrescale_.push_back(prescaleVals.first * prescaleVals.second);
      }
    }
    for (size_t i=0; i<jetTrigNamesV_.size(); ++i) {
      const std::string trigName=jetTrigNamesV_.at(i);
      id= helper_findTrigger(evTrigNames.triggerNames(),trigName);
      if (id==evTrigNames.size()) {
	jetTrigFired_.push_back(0);
	jetTrigPrescale_.push_back(-1);
	continue;
      }
      int fired= triggerResults->accept(id);
      if (fired) jetTrigFlag=true;
      jetTrigFired_.push_back(fired);
      std::pair<int,int> prescaleVals= hltConfig_.prescaleValues(iEvent,evSetup,evTrigNames.triggerName(id));
      jetTrigPrescale_.push_back(prescaleVals.first * prescaleVals.second);
    }
  }

  if (!photonTrigFlag && !jetTrigFlag) {
    if (debug_) std::cout << "no trigger fired" << std::endl;
    return;
  }


  //
  // proceed with further checks
  //

  //  cout<<"in analyze method...."<<endl;
  tagPho_pfiso_mycharged03.clear();

  edm::Handle<std::vector<reco::GenJet>> genjets;
  edm::Handle<std::vector<reco::GenParticle> > genparticles;

  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel("particleFlow", pfHandle);

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", vtxHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  iEvent.getByLabel("gsfElectrons", gsfElectronHandle);

  edm::Handle<double> rhoHandle_2012;
  iEvent.getByLabel(rhoCollection_, rhoHandle_2012);
  rho2012_ = *(rhoHandle_2012.product());

  ///  std::cout << "getting convH" << std::endl;
  edm::Handle<reco::ConversionCollection> convH;
  iEvent.getByLabel("allConversions", convH);

  /////  std::cout << "getting beamSpotHandle" << std::endl;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);


  if(doGenJets_){
    // Get GenJets
    iEvent.getByLabel(genJetCollName_,genjets);
    if(!genjets.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find GenJet vector named " << genJetCollName_ << ".\n";
      return;
    }
    nGenJets_= genjets->size();

    // Get GenParticles
    iEvent.getByLabel(genParticleCollName_,genparticles);
    if(!genparticles.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find GenParticle vector named " << genParticleCollName_ << ".\n";
      return;
    }

    // Get weights
    edm::Handle<GenEventInfoProduct> genEventInfoProduct;
    iEvent.getByLabel(genEventInfoName_, genEventInfoProduct);
    if(!genEventInfoProduct.isValid()){
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find GenEventInfoProduct named " << genEventInfoName_ << " \n";
      return;
    }
    eventWeight_ = genEventInfoProduct->weight();
  }

  runNumber_ = iEvent.id().run();
  lumiBlock_ = iEvent.id().luminosityBlock();
  eventNumber_ = iEvent.id().event();

  //
  // Fill photon info
  //
  edm::Handle<edm::ValueMap<Bool_t> > loosePhotonQual, tightPhotonQual;
  iEvent.getByLabel("PhotonIDProd", "PhotonCutBasedIDLoose", loosePhotonQual);
  iEvent.getByLabel("PhotonIDProd", "PhotonCutBasedIDTight", tightPhotonQual);
  if (!loosePhotonQual.isValid() || !tightPhotonQual.isValid()) {
    std::cout << "failed to get photon qualifiers" << std::endl;
  }

  // fill tag photon variables
  if (!photon_tag.isValid()) {
    tagPho_pt_=-1;
    pho_2nd_pt_=-1;
    tagPho_energy_=-1;
    tagPho_eta_=0;
    tagPho_phi_=0;
    tagPho_sieie_=0;
    tagPho_HoE_=0;
    tagPho_r9_=0;
    tagPho_EcalIsoDR04_=0;
    tagPho_HcalIsoDR04_=0;
    tagPho_HcalIsoDR0412_=0;
    tagPho_TrkIsoHollowDR04_=0;
    tagPho_pfiso_myphoton03_=0;
    tagPho_pfiso_myneutral03_=0;
    tagPho_pfiso_mycharged03.clear();
    tagPho_pixelSeed_=0;
    tagPho_ConvSafeEleVeto_=0;
    tagPho_idTight_=0;
    tagPho_idLoose_=0;
    tagPho_genPt_=0;
    tagPho_genEnergy_=0;
    tagPho_genEta_=0;
    tagPho_genPhi_=0;
    tagPho_genDeltaR_=0;
  }
  else {
  tagPho_pt_    = photon_tag.photon()->pt();
  pho_2nd_pt_   = (photon_2nd.photon()) ? photon_2nd.photon()->pt() : -1.;
  tagPho_energy_     = photon_tag.photon()->energy();
  tagPho_eta_   = photon_tag.photon()->eta();
  tagPho_phi_   = photon_tag.photon()->phi();
  tagPho_sieie_ = photon_tag.photon()->sigmaIetaIeta();
  tagPho_HoE_   = photon_tag.photon()->hadTowOverEm();
  tagPho_r9_    = photon_tag.photon()->r9();
  tagPho_pixelSeed_ = photon_tag.photon()->hasPixelSeed();
  tagPho_TrkIsoHollowDR04_ =  photon_tag.photon()->trkSumPtHollowConeDR04();
  tagPho_EcalIsoDR04_ = photon_tag.photon()->ecalRecHitSumEtConeDR04();
  tagPho_HcalIsoDR04_ = photon_tag.photon()->hcalTowerSumEtConeDR04();
  tagPho_HcalIsoDR0412_ = photon_tag.photon()->hcalTowerSumEtConeDR04() + (photon_tag.photon()->hadronicOverEm() - photon_tag.photon()->hadTowOverEm())*(photon_tag.photon()->energy()/cosh((photon_tag.photon()->eta())));

  tagPho_pfiso_myphoton03_  = pfEcalIso(photon_tag.photon(), pfHandle, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
  tagPho_pfiso_myneutral03_ = pfHcalIso(photon_tag.photon(), pfHandle, 0.3, 0.0, reco::PFCandidate::h0);
  tagPho_pfiso_mycharged03.push_back(pfTkIsoWithVertex(photon_tag.photon(), pfHandle, vtxHandle, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h));

  tagPho_ConvSafeEleVeto_ = ((int)ConversionTools::hasMatchedPromptElectron(photon_tag.photon()->superCluster(), gsfElectronHandle, convH, beamSpotHandle->position()));

  edm::Ref<reco::PhotonCollection> photonRef(photons, photon_tag.idx());
  tagPho_idLoose_ = (loosePhotonQual.isValid()) ? (*loosePhotonQual)[photonRef] : -1;
  tagPho_idTight_ = (tightPhotonQual.isValid()) ? (*tightPhotonQual)[photonRef] : -1;
  //std::cout << "photon tag ID = " << tagPho_idLoose_ << " and " << tagPho_idTight_ << std::endl;

  tagPho_genPt_=0;
  tagPho_genEnergy_=0;
  tagPho_genEta_=0;
  tagPho_genPhi_=0;
  tagPho_genDeltaR_=0;
  if (doGenJets_) {
    tagPho_genDeltaR_=9999.;
    for (std::vector<reco::GenParticle>::const_iterator itmc=genparticles->begin();
	 itmc!=genparticles->end(); itmc++) {
      if (itmc->status() == 1 && itmc->pdgId()==22) {
	float dR= deltaR(tagPho_eta_,tagPho_phi_,
			 itmc->eta(),itmc->phi());
	if (dR < tagPho_genDeltaR_) {
	  tagPho_genPt_     = itmc->pt();
	  tagPho_genEnergy_ = itmc->energy();
	  tagPho_genEta_    = itmc->eta();
	  tagPho_genPhi_    = itmc->phi();
	  tagPho_genDeltaR_ = dR;
	}
      }
    }
  }
  }

  // Run over caloJets //

  if (doCaloJets_ && (nCaloJets_>0)) {
        // Get jet corrections
    const JetCorrector* correctorCalo = JetCorrector::getJetCorrector(caloJetCorrName_, evSetup);

    //////////////////////////////
    // Event Selection
    //////////////////////////////

    // sort jets by corrected et
    std::set<CaloJetCorretPair, CaloJetCorretPairComp> calojetcorretpairset;
    for(reco::CaloJetCollection::const_iterator it=calojets->begin(); it!=calojets->end(); ++it) {
      const reco::CaloJet* jet=&(*it);
      calojetcorretpairset.insert( CaloJetCorretPair(jet, correctorCalo->correction(jet->p4())) );
    }
    
    // highest two (corrected) et jets
    CaloJetCorretPair calo_tag, calo_probe, calo_third;
    int cntr=0;
    for(std::set<CaloJetCorretPair, CaloJetCorretPairComp>::const_iterator it=calojetcorretpairset.begin(); it!=calojetcorretpairset.end(); ++it) {
      CaloJetCorretPair jet=(*it);
      ++cntr;
      if(cntr==1) calo_tag=jet;
      else if(cntr==2) calo_probe=jet;
      else if (cntr==3) calo_third=jet;
      else break;
    }

    // determine which cut results in failure
    int failSelCalo=0;
    
    // Selection cuts
    if (calo_tag.scaledEt() < jetEtMin_) failSelCalo |= 1;
    if (calc_dPhi(photon_tag,calo_tag) < photonJetDPhiMin_) failSelCalo |= 2;
    if (deltaR(photon_tag,calo_tag.jet())<0.5) failSelCalo |= 4;
    if (calo_probe.isValid() && (calo_probe.scaledEt() > jet2EtMax_))
      failSelCalo |= 8;
    if (calo_third.isValid() && (calo_third.scaledEt() > jet3EtMax_))
      failSelCalo |= 16;

    if (!failSelCalo) {
      // a good event

      // dump
      if(debug_) {
	std::cout << "Run: " << iEvent.id().run() << "; Event: " << iEvent.id().event() << std::endl;
	for(reco::CaloJetCollection::const_iterator it=calojets->begin(); it!=calojets->end(); ++it) {
	  const reco::CaloJet *jet=&(*it);
	  std::cout << "istag=" << (jet==calo_tag.jet()) << "; isprobe=" << (jet==calo_probe.jet()) << "; et=" << jet->et() << "; eta=" << jet->eta() << std::endl;
	}
      }

      // fill tag jet variables
      tcalojet_pt_    = calo_tag.jet()->pt();
      tcalojet_p_     = calo_tag.jet()->p();
      tcalojet_eta_   = calo_tag.jet()->eta();
      tcalojet_phi_   = calo_tag.jet()->phi();
      tcalojet_emf_   = calo_tag.jet()->emEnergyFraction();
      tcalojet_scale_ = calo_tag.scale();
      tcalojet_EBE_   = calo_tag.jet()->emEnergyInEB();
      tcalojet_EEE_   = calo_tag.jet()->emEnergyInEE();
      tcalojet_HBE_   = calo_tag.jet()->hadEnergyInHB();
      tcalojet_HEE_   = calo_tag.jet()->hadEnergyInHE();
      tcalojet_HFEtot_  = calo_tag.jet()->emEnergyInHF() + calo_tag.jet()->hadEnergyInHF();
      tcalojet_HFEcalE_ = calo_tag.jet()->emEnergyInHF();
      tcalojet_HFHadE_  = calo_tag.jet()->hadEnergyInHF();
      tcalojet_ntwrs_=0;
      std::vector<CaloTowerPtr> tagconst= calo_tag.jet()->getCaloConstituents();
      for(std::vector<CaloTowerPtr>::const_iterator it=tagconst.begin(); it!=tagconst.end(); ++it) {
	int ieta=(*it)->id().ieta();
	//int ietaAbs=(*it)->id().ietaAbs();
	tcalojet_twr_ieta_[tcalojet_ntwrs_]=ieta;
	tcalojet_twr_totE_[tcalojet_ntwrs_] = (*it)->energy();
	tcalojet_twr_eme_ [tcalojet_ntwrs_] = (*it)->emEnergy();
	tcalojet_twr_hade_[tcalojet_ntwrs_] = (*it)->hadEnergy();
	tcalojet_twr_hadoe_[tcalojet_ntwrs_] = (*it)->outerEnergy();
	++tcalojet_ntwrs_;
      }

      // fill probe jet variables
      if (!calo_probe.isValid()) {
	pcalojet_pt_=0;
	pcalojet_p_=0;
	pcalojet_eta_=0;
	pcalojet_phi_=0;
	pcalojet_emf_=0;
	pcalojet_scale_=0;
	pcalojet_EBE_=0;
	pcalojet_EEE_=0;
	pcalojet_HBE_=0;
	pcalojet_HEE_=0;
	pcalojet_HFEtot_=0;
	pcalojet_ntwrs_=0;
      }
      else { // valid pointer
      pcalojet_pt_    = calo_probe.jet()->pt();
      pcalojet_p_     = calo_probe.jet()->p();
      pcalojet_eta_   = calo_probe.jet()->eta();
      pcalojet_phi_   = calo_probe.jet()->phi();
      pcalojet_emf_   = calo_probe.jet()->emEnergyFraction();
      pcalojet_scale_ = calo_probe.scale();
      pcalojet_EBE_   = calo_probe.jet()->emEnergyInEB();
      pcalojet_EEE_   = calo_probe.jet()->emEnergyInEE();
      pcalojet_HBE_   = calo_probe.jet()->hadEnergyInHB();
      pcalojet_HEE_   = calo_probe.jet()->hadEnergyInHE();
      pcalojet_HFEtot_ = calo_probe.jet()->emEnergyInHF() + calo_probe.jet()->hadEnergyInHF();
      pcalojet_ntwrs_=0;
      std::vector<CaloTowerPtr> probeconst=calo_probe.jet()->getCaloConstituents();
      for(std::vector<CaloTowerPtr>::const_iterator it=probeconst.begin(); it!=probeconst.end(); ++it) {
	int ieta=(*it)->id().ieta();
	//int ietaAbs=(*it)->id().ietaAbs();
	pcalojet_twr_ieta_[pcalojet_ntwrs_]=ieta;
	pcalojet_twr_eme_[pcalojet_ntwrs_] = (*it)->emEnergy();
	pcalojet_twr_hade_[pcalojet_ntwrs_] = (*it)->hadEnergy();
	++pcalojet_ntwrs_;
      }

      // save info about 3rd jet
      if (!calo_third.isValid()) {
	calo_thirdjet_et_=0;
	calo_thirdjet_pt_=0;
	calo_thirdjet_p_=0;
	calo_thirdjet_px_=0;
	calo_thirdjet_py_=0;
	calo_thirdjet_E_=0;
	calo_thirdjet_eta_=0;
	calo_thirdjet_phi_=0;
	calo_thirdjet_scale_=0;
      }
      else {
	calo_thirdjet_et_ = calo_third.jet()->et();
	calo_thirdjet_pt_ = calo_third.jet()->pt();
	calo_thirdjet_p_ = calo_third.jet()->p();
	calo_thirdjet_px_= calo_third.jet()->px();
	calo_thirdjet_py_ = calo_third.jet()->py();
	calo_thirdjet_E_ = calo_third.jet()->energy();
	calo_thirdjet_eta_ = calo_third.jet()->eta();
	calo_thirdjet_phi_ = calo_third.jet()->phi();
	calo_thirdjet_scale_ = calo_third.scale();
      }
      }

      if(doGenJets_){
	// fill genjet tag/probe variables
	tcalojet_gendr_ = 99999.;
	tcalojet_genpt_ = 0;
	tcalojet_genp_  = 0;
	pcalojet_gendr_ = 99999.;
	pcalojet_genpt_ = 0;
	pcalojet_genp_  = 0;
	for(std::vector<reco::GenJet>::const_iterator it=genjets->begin(); it!=genjets->end(); ++it){
	  const reco::GenJet* jet=&(*it);
	  if (!calo_probe.isValid()) {
	    pcalojet_gendr_=0;
	    pcalojet_genpt_=0;
	    pcalojet_genp_=0;
	  }
	  else {
	    double drProbe=deltaR(jet, calo_probe.jet());
	    if(drProbe<pcalojet_gendr_) {
	      pcalojet_gendr_ = drProbe;
	      pcalojet_genpt_ = jet->pt();
	      pcalojet_genp_ = jet->p();
	    }
	  }
	  double dr=deltaR(jet, calo_tag.jet());
	  if(dr<tcalojet_gendr_) {
	    tcalojet_gendr_ = dr;
	    tcalojet_genpt_ = jet->pt();
	    tcalojet_genp_ = jet->p();
	  }
	}
      }

      calo_tree_->Fill();
    }

  }

  // Run over PFJets //

  if(doPFJets_ && (nPFJets_>0)){
    unsigned int debugEvent = 0;


    // Get RecHits in HB and HE
    edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hbhereco;
    iEvent.getByLabel(hbheRecHitName_,hbhereco);
    if(!hbhereco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HBHERecHit named " << hbheRecHitName_ << ".\n";
      return;
    }
    
    // Get RecHits in HF
    edm::Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>> hfreco;
    iEvent.getByLabel(hfRecHitName_,hfreco);
    if(!hfreco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HFRecHit named " << hfRecHitName_ << ".\n";
      return;
    }

    // Get RecHits in HO
    edm::Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>> horeco;
    iEvent.getByLabel(hoRecHitName_,horeco);
    if(!horeco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HORecHit named " << hoRecHitName_ << ".\n";
      return;
    }

    // Get geometry
    edm::ESHandle<CaloGeometry> geoHandle;
    evSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *HBGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 1);
    const CaloSubdetectorGeometry *HEGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 2);
    const CaloSubdetectorGeometry *HOGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 3);
    const CaloSubdetectorGeometry *HFGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 4);
    
    int HBHE_n = 0;
    for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
      HBHE_n++;
      h_hbherecoieta_->Fill((*ith).id().ieta());
      if(iEvent.id().event() == debugEvent){
	std::cout << (*ith).id().ieta() << " " << (*ith).id().iphi() << std::endl;
	h_rechitspos_->Fill((*ith).id().ieta(), (*ith).id().iphi());
      }
    }
    int HF_n = 0;
    for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
      HF_n++;
      if(iEvent.id().event() == debugEvent){
	h_rechitspos_->Fill((*ith).id().ieta(), (*ith).id().iphi());
      }
    }
    int HO_n = 0;
    for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
      HO_n++;
      if(iEvent.id().event() == debugEvent){
	h_rechitspos_->Fill((*ith).id().ieta(), (*ith).id().iphi());
      }
    }
    h_HBHE_n_->Fill(HBHE_n);
    h_HF_n_->Fill(HF_n);
    h_HO_n_->Fill(HO_n);
    
    // Get jet corrections
    const JetCorrector* correctorPF = JetCorrector::getJetCorrector(pfJetCorrName_,evSetup);
    
    //////////////////////////////
    // Event Selection
    //////////////////////////////
    
    // sort jets by corrected et
    std::set<PFJetCorretPair, PFJetCorretPairComp> pfjetcorretpairset;
    for(reco::PFJetCollection::const_iterator it=pfjets->begin(); it!=pfjets->end(); ++it) {
      const reco::PFJet* jet=&(*it);
      //double minDr=99999;
      //double dR= deltaR(photon_tag,jet);
      //if (dR < minDr) minDr=dR;
      //if(minDr<0.5) continue;
      // shorter version
      if ((deltaR(photon_tag,jet)<0.5) &&
	  (it!=pfjets->begin())) // do not allow photon to be inside leading jet
	  continue;
      //int index = it-pfjets->begin();
      //edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfjets,index));
      //reco::PFJetRef jetRef = it->castTo<reco::PFJetRef>();
      double jec = correctorPF->correction(*it, iEvent, evSetup);
      //cout<<index<<'\t'<<jec<<'\t'<<it->et()<<'\t'<<it->pt()<<endl;
      pfjetcorretpairset.insert( PFJetCorretPair(jet, jec));//correctorPF->correction(jet->p4())) );
    }

    PFJetCorretPair pfjet_probe;
    PFJetCorretPair pf_2ndjet;
    PFJetCorretPair pf_3rdjet;
    int cntr=0;
    for(std::set<PFJetCorretPair, PFJetCorretPairComp>::const_iterator it=pfjetcorretpairset.begin(); it!=pfjetcorretpairset.end(); ++it) {
      PFJetCorretPair jet=(*it);
      ++cntr;
      if(cntr==1) pfjet_probe = jet;
      else if(cntr==2) pf_2ndjet = jet;
      else if(cntr==3) pf_3rdjet = jet;
      else break;
    }

    // Check selection
    int failSelPF = 0;
    
    if (pfjet_probe.scaledEt() < jetEtMin_) failSelPF |= 1;
    if (calc_dPhi(photon_tag,pfjet_probe) < photonJetDPhiMin_) failSelPF |= 2;
    if (deltaR(photon_tag,pfjet_probe.jet())<0.5) failSelPF |= 4;
    if (pf_2ndjet.isValid() && (pf_2ndjet.scaledEt() > jet2EtMax_))
      failSelPF |= 8;
    if (pf_3rdjet.isValid() && (pf_3rdjet.scaledEt() > jet3EtMax_))
      failSelPF |= 16;

    if (!failSelPF) {
      // a good event

      // prepare the container -- later in the loop
      // clear_leadingPfJetVars();

      // put values into 3rd jet quantities
      if (pf_3rdjet.isValid()) {
	pf_thirdjet_et_ = pf_3rdjet.jet()->et();
	pf_thirdjet_pt_ = pf_3rdjet.jet()->pt();
	pf_thirdjet_p_  = pf_3rdjet.jet()->p();
	pf_thirdjet_px_ = pf_3rdjet.jet()->px();
	pf_thirdjet_py_ = pf_3rdjet.jet()->py();
	pf_thirdjet_E_  = pf_3rdjet.jet()->energy();
	pf_thirdjet_eta_= pf_3rdjet.jet()->eta();
	pf_thirdjet_phi_= pf_3rdjet.jet()->phi();
	pf_thirdjet_scale_= pf_3rdjet.scale();
      }
      else {
	pf_thirdjet_et_ = 0;
	pf_thirdjet_pt_ = pf_thirdjet_p_ = 0;
	pf_thirdjet_px_ = pf_thirdjet_py_ = 0;
	pf_thirdjet_E_ = pf_thirdjet_eta_ = pf_thirdjet_phi_ = 0;
	pf_thirdjet_scale_=0;
      }

    int types = 0;
    int ntypes = 0;
    
    /////////////////////////////////////////////
    // Get PF constituents and fill HCAL towers
    /////////////////////////////////////////////
    
    // fill jet variables
    // First start from a second jet, then fill the first jet
    PFJetCorretPair pfjet_probe_store = pfjet_probe;
    for (int iJet=2; iJet>0; iJet--) {
      // prepare the container
      clear_leadingPfJetVars();

      if (iJet==2) pfjet_probe= pf_2ndjet;
      else pfjet_probe = pfjet_probe_store;

      if(!pfjet_probe.jet()) {
	if (iJet==2) {
	  // put zeros into 2nd jet quantities
	  copy_leadingPfJetVars_to_pfJet2();
	}
	else {
	  std::cerr << "error in the code: leading pf jet is null"<< std::endl;
	}
	continue;
      }

      // temporary variables
      std::map<int,std::pair<int,std::set<float>>> ppfjet_rechits;
      std::map<float,int> ppfjet_clusters;

      // fill the values
      ppfjet_pt_    = pfjet_probe.jet()->pt();
      ppfjet_p_     = pfjet_probe.jet()->p();
      ppfjet_E_     = pfjet_probe.jet()->energy();
      ppfjet_eta_   = pfjet_probe.jet()->eta();
      ppfjet_phi_   = pfjet_probe.jet()->phi();
      ppfjet_NeutralHadronFrac_  = pfjet_probe.jet()->neutralHadronEnergyFraction();
      ppfjet_NeutralEMFrac_      = pfjet_probe.jet()->neutralEmEnergyFraction();
      ppfjet_nConstituents_      = pfjet_probe.jet()->nConstituents();
      ppfjet_ChargedHadronFrac_  = pfjet_probe.jet()->chargedHadronEnergyFraction();
      ppfjet_ChargedMultiplicity_= pfjet_probe.jet()->chargedMultiplicity();
      ppfjet_ChargedEMFrac_      = pfjet_probe.jet()->chargedEmEnergyFraction();
      ppfjet_scale_ = pfjet_probe.scale();
      ppfjet_ntwrs_=0;
      ppfjet_cluster_n_=0;
      ppfjet_ncandtracks_=0;
      
      if(iEvent.id().event() == debugEvent){
	std::cout << "Probe eta: " << ppfjet_eta_ << " phi: " << ppfjet_phi_ << std::endl;
      }
      //std::cout << "Probe eta: " << ppfjet_eta_ << " phi: " << ppfjet_phi_ << std::endl; //debug
      
      // Get PF constituents and fill HCAL towers
      std::vector<reco::PFCandidatePtr> probeconst=pfjet_probe.jet()->getPFConstituents();
      int iPF=0;
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=probeconst.begin(); it!=probeconst.end(); ++it){
	bool hasTrack = false;
	reco::PFCandidate::ParticleType candidateType = (*it)->particleId();
	iPF++;

	// some debug print
	//if ((*it)->hoEnergy()>0.1) std::cout << " hoEn: " << (*it)->hoEnergy() << ": raw=" << (*it)->rawHoEnergy() << "\n";

	if (0 && (candidateType==reco::PFCandidate::h)) {
	  std::cout << iPF << "hadron PF info:\n";
	  std::cout << " energy: " << (*it)->energy() << "\n";
	  std::cout << " ecalEn: " << (*it)->ecalEnergy() << ": raw=" << (*it)->rawEcalEnergy() << "\n";
	  std::cout << " hcalEn: " << (*it)->hcalEnergy() << ": raw=" << (*it)->rawHcalEnergy() << "\n";
	  std::cout << " hoEn: " << (*it)->hoEnergy() << ": raw=" << (*it)->rawHoEnergy() << "\n";
	  printElementsInBlocks(*it);
	}

	// store information
	switch(candidateType){
	case reco::PFCandidate::X:
	  ppfjet_unkown_E_ += (*it)->energy();
	  ppfjet_unkown_px_ += (*it)->px();
	  ppfjet_unkown_py_ += (*it)->py();
	  ppfjet_unkown_pz_ += (*it)->pz();
	  ppfjet_unkown_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_unkown_n_++;
	  continue;
	case reco::PFCandidate::h:
	  {
	    ppfjet_had_E_.push_back((*it)->energy());
	    ppfjet_had_px_.push_back((*it)->px());
	    ppfjet_had_py_.push_back((*it)->py());
	    ppfjet_had_pz_.push_back((*it)->pz());
	    ppfjet_had_EcalE_.push_back((*it)->ecalEnergy());
	    ppfjet_had_rawHcalE_.push_back((*it)->rawHcalEnergy());
	    ppfjet_had_id_.push_back(0);
	    ppfjet_had_ntwrs_.push_back(0);
	    ppfjet_had_n_++;
	    
	    if(doGenJets_){
	      //cout<<"38%"<<endl;
	      float gendr = 99999;
	      float genE = 0;
	      int genpdgId = 0;
	      for(std::vector<reco::GenParticle>::const_iterator itmc = genparticles->begin(); itmc != genparticles->end(); itmc++){
		if(itmc->status() == 1 && itmc->pdgId() > 100){
		  double dr = deltaR((*it)->eta(),(*it)->phi(),itmc->eta(),itmc->phi());
		  if(dr < gendr){
		    gendr = dr;
		    genE = itmc->energy();
		    genpdgId = itmc->pdgId();
		  }
		}
	      }
	      ppfjet_had_E_mctruth_.push_back(genE);
	      ppfjet_had_mcpdgId_.push_back(genpdgId);
	    }
	    
	    reco::TrackRef trackRef = (*it)->trackRef();
	    if(trackRef.isNonnull()){
	      reco::Track track = *trackRef;
	      ppfjet_candtrack_px_.push_back(track.px());
	      ppfjet_candtrack_py_.push_back(track.py());
	      ppfjet_candtrack_pz_.push_back(track.pz());
	      ppfjet_candtrack_EcalE_.push_back((*it)->ecalEnergy());
	      ppfjet_had_candtrackind_.push_back(ppfjet_ncandtracks_);
	      hasTrack = true;
	      ppfjet_ncandtracks_++;
	    }
	    else{
	      ppfjet_had_candtrackind_.push_back(-2);
	    }
	  }
	  break;
	case reco::PFCandidate::e:
	  ppfjet_electron_E_ += (*it)->energy();
	  ppfjet_electron_px_ += (*it)->px();
	  ppfjet_electron_py_ += (*it)->py();
	  ppfjet_electron_pz_ += (*it)->pz();
	  ppfjet_electron_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_electron_n_++;
	  continue;
	case reco::PFCandidate::mu:
	  ppfjet_muon_E_ += (*it)->energy();
	  ppfjet_muon_px_ += (*it)->px();
	  ppfjet_muon_py_ += (*it)->py();
	  ppfjet_muon_pz_ += (*it)->pz();
	  ppfjet_muon_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_muon_n_++;
	  continue;
	case reco::PFCandidate::gamma:
	  ppfjet_photon_E_ += (*it)->energy();
	  ppfjet_photon_px_ += (*it)->px();
	  ppfjet_photon_py_ += (*it)->py();
	  ppfjet_photon_pz_ += (*it)->pz();
	  ppfjet_photon_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_photon_n_++;
	  continue;
	case reco::PFCandidate::h0:
	  {
	    ppfjet_had_E_.push_back((*it)->energy());
	    ppfjet_had_px_.push_back((*it)->px());
	    ppfjet_had_py_.push_back((*it)->py());
	    ppfjet_had_pz_.push_back((*it)->pz());
	    ppfjet_had_EcalE_.push_back((*it)->ecalEnergy());
	    ppfjet_had_rawHcalE_.push_back((*it)->rawHcalEnergy());
	    ppfjet_had_id_.push_back(1);
	    ppfjet_had_candtrackind_.push_back(-1);
	    ppfjet_had_ntwrs_.push_back(0);
	    ppfjet_had_n_++;

	    if(doGenJets_){
 //cout<<"44%"<<endl;
	      float gendr = 99999;
	      float genE = 0;
	      int genpdgId = 0;
	      for(std::vector<reco::GenParticle>::const_iterator itmc = genparticles->begin(); itmc != genparticles->end(); itmc++){
		if(itmc->status() == 1 && itmc->pdgId() > 100){
		  double dr = deltaR((*it)->eta(),(*it)->phi(),itmc->eta(),itmc->phi());
		  if(dr < gendr){
		    gendr = dr;
		    genE = itmc->energy();
		    genpdgId = itmc->pdgId();
		  }
		}
	      }
	      ppfjet_had_E_mctruth_.push_back(genE);
	      ppfjet_had_mcpdgId_.push_back(genpdgId);
	    }
	    
	    break;
	  }
	case reco::PFCandidate::h_HF:
	  {
	    ppfjet_had_E_.push_back((*it)->energy());
	    ppfjet_had_px_.push_back((*it)->px());
	    ppfjet_had_py_.push_back((*it)->py());
	    ppfjet_had_pz_.push_back((*it)->pz());
	    ppfjet_had_EcalE_.push_back((*it)->ecalEnergy());
	    ppfjet_had_rawHcalE_.push_back((*it)->rawHcalEnergy());
	    ppfjet_had_id_.push_back(2);
	    ppfjet_had_candtrackind_.push_back(-1);
	    ppfjet_had_ntwrs_.push_back(0);
	    ppfjet_had_n_++;
	    
	    if(doGenJets_){
 //cout<<"45%"<<endl;
	      float gendr = 99999;
	      float genE = 0;
	      int genpdgId = 0;
	      for(std::vector<reco::GenParticle>::const_iterator itmc = genparticles->begin(); itmc != genparticles->end(); itmc++){
		if(itmc->status() == 1 && itmc->pdgId() > 100){
		  double dr = deltaR((*it)->eta(),(*it)->phi(),itmc->eta(),itmc->phi());
		  if(dr < gendr){
		    gendr = dr;
		    genE = itmc->energy();
		    genpdgId = itmc->pdgId();
		  }
		}
	      }
	      ppfjet_had_E_mctruth_.push_back(genE);
	      ppfjet_had_mcpdgId_.push_back(genpdgId);
	    }
	    
	    break;
	  }
	case reco::PFCandidate::egamma_HF:
	  {
	    ppfjet_had_E_.push_back((*it)->energy());
	    ppfjet_had_px_.push_back((*it)->px());
	    ppfjet_had_py_.push_back((*it)->py());
	    ppfjet_had_pz_.push_back((*it)->pz());
	    ppfjet_had_EcalE_.push_back((*it)->ecalEnergy());
	    ppfjet_had_rawHcalE_.push_back((*it)->rawHcalEnergy());
	    ppfjet_had_id_.push_back(3);
	    ppfjet_had_candtrackind_.push_back(-1);
	    ppfjet_had_ntwrs_.push_back(0);
	    ppfjet_had_n_++;

	    if(doGenJets_){
 //cout<<"50%"<<endl;
	      float gendr = 99999;
	      float genE = 0;
	      int genpdgId = 0;
	      for(std::vector<reco::GenParticle>::const_iterator itmc = genparticles->begin(); itmc != genparticles->end(); itmc++){
		if(itmc->status() == 1 && itmc->pdgId() > 100){
		  double dr = deltaR((*it)->eta(),(*it)->phi(),itmc->eta(),itmc->phi());
		  if(dr < gendr){
		    gendr = dr;
		    genE = itmc->energy();
		    genpdgId = itmc->pdgId();
		  }
		}
	      }
	      ppfjet_had_E_mctruth_.push_back(genE);
	      ppfjet_had_mcpdgId_.push_back(genpdgId);
	    }
	    
	    break;
	  }
	}

	float HFHAD_E = 0;
	float HFEM_E = 0;
	int HFHAD_n_ = 0;
	int HFEM_n_ = 0;
	int HF_type_ = 0;
	int maxElement=(*it)->elementsInBlocks().size();
	for(int e=0; e<maxElement; ++e){
	  // Get elements from block
	  reco::PFBlockRef blockRef = (*it)->elementsInBlocks()[e].first;
	  const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	  for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){
	      if(elements[iEle].type() == reco::PFBlockElement::HCAL){ // Element is HB or HE
		HF_type_ |= 0x1;
		// Get cluster and hits
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;
		double cluster_dR = deltaR(ppfjet_eta_,ppfjet_phi_,cluster.eta(),cluster.phi());
		if(ppfjet_clusters.count(cluster_dR) == 0){
		  ppfjet_clusters[cluster_dR] = ppfjet_cluster_n_;
		  ppfjet_cluster_eta_.push_back(cluster.eta());
		  ppfjet_cluster_phi_.push_back(cluster.phi());
		  ppfjet_cluster_dR_.push_back(cluster_dR);
		  ppfjet_cluster_n_++;
		}
		int cluster_ind = ppfjet_clusters[cluster_dR];
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();

		// Run over hits and match
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  int etaPhiPF = getEtaPhi(hitsAndFracs[iHit].first);

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = getEtaPhi((*ith).id());
		    if(etaPhiPF == etaPhiRecHit){
		      ppfjet_had_ntwrs_.at(ppfjet_had_n_ - 1)++;
		      if(ppfjet_rechits.count((*ith).id()) == 0){
			ppfjet_twr_ieta_.push_back((*ith).id().ieta());
			ppfjet_twr_iphi_.push_back((*ith).id().iphi());
			ppfjet_twr_depth_.push_back((*ith).id().depth());
			ppfjet_twr_subdet_.push_back((*ith).id().subdet());
			ppfjet_twr_hade_.push_back((*ith).energy());
			ppfjet_twr_frac_.push_back(hitsAndFracs[iHit].second);
			ppfjet_rechits[(*ith).id()].second.insert(hitsAndFracs[iHit].second);
			ppfjet_twr_hadind_.push_back(ppfjet_had_n_ - 1);
			ppfjet_twr_elmttype_.push_back(0);
			ppfjet_twr_clusterind_.push_back(cluster_ind);
			if(hasTrack){
			  ppfjet_twr_candtrackind_.push_back(ppfjet_ncandtracks_ - 1);
			}
			else{
			  ppfjet_twr_candtrackind_.push_back(-1);
			}
			switch((*ith).id().subdet()){
			case HcalSubdetector::HcalBarrel:
			  {
			    const CaloCellGeometry *thisCell = HBGeom->getGeometry((*ith).id().rawId());
			    const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
			    float avgeta = (cv[0].eta() + cv[2].eta())/2.0;
			    float avgphi = (static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
			    if(cv[0].phi() < cv[2].phi()) std::cout << "pHB" << cv[0].phi() << " " << cv[2].phi() << std::endl;
			    if(cv[0].phi() < cv[2].phi()) avgphi = (2.0*3.141592653 + static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;

			    //if(pf_Event_ == 9413996) //debug
			    //printf("pHB ieta: %3d iphi: %2d eta0: %6f phi0: %6f eta2: %6f phi2: %6f dR: %f\n",(*ith).id().ieta(),(*ith).id().iphi(),static_cast<double>(cv[0].eta()),static_cast<double>(cv[0].phi()),static_cast<double>(cv[2].eta()),static_cast<double>(cv[2].phi()),static_cast<double>(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi))); //debug
			    ppfjet_twr_dR_.push_back(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi));
			    break;
			  }
			case HcalSubdetector::HcalEndcap:
			  {
			    const CaloCellGeometry *thisCell = HEGeom->getGeometry((*ith).id().rawId());
			    const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
			    float avgeta = (cv[0].eta() + cv[2].eta())/2.0;
			    float avgphi = (static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
			    if(cv[0].phi() < cv[2].phi()) std::cout << "pHE" << cv[0].phi() << " " << cv[2].phi() << std::endl;
			    if(cv[0].phi() < cv[2].phi()) avgphi = (2.0*3.141592653 + static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;

			    //printf("pHE ieta: %3d iphi: %2d eta0: %6f phi0: %6f eta2: %6f phi2: %6f dR: %f\n",(*ith).id().ieta(),(*ith).id().iphi(),static_cast<double>cv[0].eta(),static_cast<double>cv[0].phi(),static_cast<double>cv[2].eta(),static_cast<double>cv[2].phi(),static_cast<double>deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi)); //debug
			    /*printf("  cv0: %f cv2: %f sum: %f avg: %f\n",
				   static_cast<double>cv[0].phi(),
				   static_cast<double>cv[2].phi(),
				   (static_cast<double>cv[0].phi() + static_cast<double>cv[2].phi()),
				   static_cast<double>((cv[0].phi() + cv[2].phi())/2.0));*/
			    ppfjet_twr_dR_.push_back(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi));
			    break;
			  }
			default:
			  ppfjet_twr_dR_.push_back(-1);
			  break;
			}
			ppfjet_rechits[(*ith).id()].first = ppfjet_ntwrs_;
			++ppfjet_ntwrs_;
		      }
		      else if(ppfjet_rechits[(*ith).id()].second.count(hitsAndFracs[iHit].second) == 0){
			ppfjet_twr_frac_.at(ppfjet_rechits[(*ith).id()].first) += hitsAndFracs[iHit].second;
			if(cluster_dR < ppfjet_cluster_dR_.at(ppfjet_twr_clusterind_.at(ppfjet_rechits[(*ith).id()].first))){
			  ppfjet_twr_clusterind_.at(ppfjet_rechits[(*ith).id()].first) = cluster_ind;
			}
			ppfjet_rechits[(*ith).id()].second.insert(hitsAndFracs[iHit].second);
		      }
		    } // Test if ieta,iphi matches
		  } // Loop over rechits
		} // Loop over hits
	      } // Test if element is from HCAL
	      else if(elements[iEle].type() == reco::PFBlockElement::HFHAD){ // Element is HF
		types |= 0x2;
		ntypes++;
		HFHAD_n_++;
		HF_type_ |= 0x2;
		
		h_etaHFHAD_->Fill((*it)->eta());

		for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		  if((*ith).id().depth() == 1) continue; // Remove long fibers
		  const CaloCellGeometry *thisCell = HFGeom->getGeometry((*ith).id().rawId());
		  const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
		  
		  bool passMatch = false;
		  if((*it)->eta() < cv[0].eta() && (*it)->eta() > cv[2].eta()){
		    if((*it)->phi() < cv[0].phi() && (*it)->phi() > cv[2].phi()) passMatch = true;
		    else if(cv[0].phi() < cv[2].phi()){
		      std::cout << "HFHAD probe" << std::endl;
		      if((*it)->phi() < cv[0].phi()) passMatch = true;
		      else if((*it)->phi() > cv[2].phi()) passMatch = true;
		    }
		  }
		  
		  if(passMatch){
		    ppfjet_had_ntwrs_.at(ppfjet_had_n_ - 1)++;
		    ppfjet_twr_ieta_.push_back((*ith).id().ieta());
		    ppfjet_twr_iphi_.push_back((*ith).id().iphi());
		    ppfjet_twr_depth_.push_back((*ith).id().depth());
		    ppfjet_twr_subdet_.push_back((*ith).id().subdet());
		    ppfjet_twr_hade_.push_back((*ith).energy());
		    ppfjet_twr_frac_.push_back(1.0);
		    ppfjet_twr_hadind_.push_back(ppfjet_had_n_ - 1);
		    ppfjet_twr_elmttype_.push_back(1);
		    ppfjet_twr_clusterind_.push_back(-1);
		    ppfjet_twr_candtrackind_.push_back(-1);
		    float avgeta = (cv[0].eta() + cv[2].eta())/2.0;
		    float avgphi = (static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
		    if(cv[0].phi() < cv[2].phi()) std::cout << "pHFhad" << cv[0].phi() << " " << cv[2].phi() << std::endl;
		    if(cv[0].phi() < cv[2].phi()) avgphi = (2.0*3.141592653 + static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
		    ppfjet_twr_dR_.push_back(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi));
		    ++ppfjet_ntwrs_;
		    HFHAD_E += (*ith).energy();
		  }
		}		
	      }
	      else if(elements[iEle].type() == reco::PFBlockElement::HFEM){ // Element is HF
		types |= 0x4;
		ntypes++;
		HFEM_n_++;
		HF_type_ |= 0x4;

		h_etaHFEM_->Fill((*it)->eta());
		
		for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		  if((*ith).id().depth() == 2) continue; // Remove short fibers
		  const CaloCellGeometry *thisCell = HFGeom->getGeometry((*ith).id().rawId());
		  const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
		  
		  bool passMatch = false;
		  if((*it)->eta() < cv[0].eta() && (*it)->eta() > cv[2].eta()){
		    if((*it)->phi() < cv[0].phi() && (*it)->phi() > cv[2].phi()) passMatch = true;
		    else if(cv[0].phi() < cv[2].phi()){
		      std::cout << "HFEM probe" << std::endl;
		      if((*it)->phi() < cv[0].phi()) passMatch = true;
		      else if((*it)->phi() > cv[2].phi()) passMatch = true;
		    }
		  }
		  
		  if(passMatch){
		    ppfjet_had_ntwrs_.at(ppfjet_had_n_ - 1)++;
		    ppfjet_twr_ieta_.push_back((*ith).id().ieta());
		    ppfjet_twr_iphi_.push_back((*ith).id().iphi());
		    ppfjet_twr_depth_.push_back((*ith).id().depth());
		    ppfjet_twr_subdet_.push_back((*ith).id().subdet());
		    ppfjet_twr_hade_.push_back((*ith).energy());
		    ppfjet_twr_frac_.push_back(1.0);
		    ppfjet_twr_hadind_.push_back(ppfjet_had_n_ - 1);
		    ppfjet_twr_elmttype_.push_back(2);
		    ppfjet_twr_clusterind_.push_back(-1);
		    ppfjet_twr_candtrackind_.push_back(-1);
		    float avgeta = (cv[0].eta() + cv[2].eta())/2.0;
		    float avgphi = (static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
		    if(cv[0].phi() < cv[2].phi()) std::cout << "pHFem" << cv[0].phi() << " " << cv[2].phi() << std::endl;
		    if(cv[0].phi() < cv[2].phi()) avgphi = (2.0*3.141592653 + static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
		    ppfjet_twr_dR_.push_back(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi));
		    ++ppfjet_ntwrs_;
		    HFEM_E += (*ith).energy();
		  }
		}
	      }
	      else if(elements[iEle].type() == reco::PFBlockElement::HO){ // Element is HO
		types |= 0x8;
		ntypes++;
		HF_type_ |= 0x8;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;
		double cluster_dR = deltaR(ppfjet_eta_,ppfjet_phi_,cluster.eta(),cluster.phi());
		if(ppfjet_clusters.count(cluster_dR) == 0){
		  ppfjet_clusters[cluster_dR] = ppfjet_cluster_n_;
		  ppfjet_cluster_eta_.push_back(cluster.eta());
		  ppfjet_cluster_phi_.push_back(cluster.phi());
		  ppfjet_cluster_dR_.push_back(cluster_dR);
		  ppfjet_cluster_n_++;
		}
		int cluster_ind = ppfjet_clusters[cluster_dR];

		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  int etaPhiPF = getEtaPhi(hitsAndFracs[iHit].first);

		  int tmpzside = ((hitsAndFracs[iHit].first.rawId() >> 13) & 0x1) ? 1 : -1;
		  int tmpieta = ((hitsAndFracs[iHit].first.rawId() >> 7) & 0x3F);
		  h_ietaHO_->Fill(tmpzside*tmpieta);

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = getEtaPhi((*ith).id());
		    if(etaPhiPF == etaPhiRecHit){
		      ppfjet_had_ntwrs_.at(ppfjet_had_n_ - 1)++;
		      if(ppfjet_rechits.count((*ith).id()) == 0){
			ppfjet_twr_ieta_.push_back((*ith).id().ieta());
			ppfjet_twr_iphi_.push_back((*ith).id().iphi());
			ppfjet_twr_depth_.push_back((*ith).id().depth());
			ppfjet_twr_subdet_.push_back((*ith).id().subdet());
			ppfjet_twr_hade_.push_back((*ith).energy());
			ppfjet_twr_frac_.push_back(hitsAndFracs[iHit].second);
			ppfjet_rechits[(*ith).id()].second.insert(hitsAndFracs[iHit].second);
			ppfjet_twr_hadind_.push_back(ppfjet_had_n_ - 1);
			ppfjet_twr_elmttype_.push_back(3);
			ppfjet_twr_clusterind_.push_back(cluster_ind);
			if(hasTrack){
			  ppfjet_twr_candtrackind_.push_back(ppfjet_ncandtracks_ - 1);
			}
			else{
			  ppfjet_twr_candtrackind_.push_back(-1);
			}
			const CaloCellGeometry *thisCell = HOGeom->getGeometry((*ith).id().rawId());
			const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
			float avgeta = (cv[0].eta() + cv[2].eta())/2.0;
			float avgphi = (static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;
			if(cv[0].phi() < cv[2].phi()) std::cout << "pHO" << cv[0].phi() << " " << cv[2].phi() << std::endl;
			if(cv[0].phi() < cv[2].phi()) avgphi = (2.0*3.141592653 + static_cast<double>(cv[0].phi()) + static_cast<double>(cv[2].phi()))/2.0;

			ppfjet_twr_dR_.push_back(deltaR(ppfjet_eta_,ppfjet_phi_,avgeta,avgphi));
			ppfjet_rechits[(*ith).id()].first = ppfjet_ntwrs_;
			++ppfjet_ntwrs_;
		      }
		      else if(ppfjet_rechits[(*ith).id()].second.count(hitsAndFracs[iHit].second) == 0){
			ppfjet_twr_frac_.at(ppfjet_rechits[(*ith).id()].first) += hitsAndFracs[iHit].second;
			if(cluster_dR < ppfjet_cluster_dR_.at(ppfjet_twr_clusterind_.at(ppfjet_rechits[(*ith).id()].first))){
			  ppfjet_twr_clusterind_.at(ppfjet_rechits[(*ith).id()].first) = cluster_ind;
			}
			ppfjet_rechits[(*ith).id()].second.insert(hitsAndFracs[iHit].second);
		      }
		    } // Test if ieta,iphi match
		  } // Loop over rechits
		} // Loop over hits
	      } // Test if element is from HO
	    } // Test for right element index
	  } // Loop over elements
	} // Loop over elements in blocks
	h_HFHAD_n_->Fill(HFHAD_n_);
	h_HFEM_n_->Fill(HFEM_n_);
	switch(candidateType){
	case reco::PFCandidate::h_HF:
	  h_HFHAD_type_->Fill(HF_type_);
	  ppfjet_had_emf_.push_back(HFEM_E/(HFEM_E + HFHAD_E));
	  break;
	case reco::PFCandidate::egamma_HF:
	  h_HFEM_type_->Fill(HF_type_);
	  ppfjet_had_emf_.push_back(-1);
	  break;
	default:
	  ppfjet_had_emf_.push_back(-1);
	  break;
	}
      } // Loop over PF constitutents
 
    if(doGenJets_){
      //cout<<"78%"<<endl;
      // fill genjet variables
      ppfjet_gendr_ = 99999.;
      ppfjet_genpt_ = 0;
      ppfjet_genp_  = 0;
      for(std::vector<reco::GenJet>::const_iterator it=genjets->begin(); it!=genjets->end(); ++it){
	const reco::GenJet* jet=&(*it);
	double dr=deltaR(jet, pfjet_probe.jet());
	if(dr<ppfjet_gendr_) {
	  ppfjet_gendr_ = dr;
	  ppfjet_genpt_ = jet->pt();
	  ppfjet_genp_ = jet->p();
	  ppfjet_genE_ = jet->energy();
	}
      }
    } // doGenJets_
    if (iJet==2) {
      copy_leadingPfJetVars_to_pfJet2();
    }
    }

    h_types_->Fill(types);
    h_ntypes_->Fill(ntypes);



    // fill photon+jet variables
 
    pf_tree_->Fill();
    }
  }
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void CalcRespCorrPhotonPlusJet::beginJob()
{

  //std::cout << "Start beginJob()" << std::endl;

  // book histograms
  rootfile_ = new TFile(rootHistFilename_.c_str(), "RECREATE");

  if(doPFJets_){
    h_types_ = new TH1D("h_types","h_types",16,0,16);
    h_ntypes_ = new TH1D("h_ntypes","h_ntypes",50,0,50);
    h_ietaHCAL_ = new TH1D("h_ietaHCAL","h_ietaHCAL",83,-41.5,41.5);
    h_etaHFHAD_ = new TH1D("h_etaHFHAD","h_etaHFHAD",100,-5.5,5.5);
    h_etaHFEM_ = new TH1D("h_etaHFEM","h_etaHFEM",100,-5.5,5.5);
    h_ietaHO_ = new TH1D("h_ietaHO","h_ietaHO",83,-41.5,41.5);
    h_HFHAD_n_ = new TH1D("h_HFHAD_n","h_HFHAD_n",10,0,10);
    h_HFEM_n_ = new TH1D("h_HFEM_n","h_HFEM_n",10,0,10);
    h_HFHAD_type_ = new TH1D("h_HFHAD_type","h_HFHAD_type",16,0,16);
    h_HFEM_type_ = new TH1D("h_HFEM_type","h_HFEM_type",16,0,16);
    h_HBHE_n_ = new TH1D("h_HBHE_n","h_HBHE_n",200,0,200);
    h_HF_n_ = new TH1D("h_HF_n","h_HF_n",200,0,200);
    h_HO_n_ = new TH1D("h_HO_n","h_HO_n",200,0,200);
    h_twrietas_ = new TH1D("h_twrietas","h_twrietas",20,0,20);
    h_rechitspos_ = new TH2D("h_rechitspos","h_rechitspos",83,-41.5,41.5,72,-0.5,71.5);
    h_hbherecoieta_ = new TH1D("h_hbherecoieta","h_hbherecoieta",83,-41.5,41.5);
  }

  // Save info about the triggers and other misc items
  {
    rootfile_->mkdir("miscItems");
    rootfile_->cd("miscItems");
    misc_tree_= new TTree("misc_tree","tree for misc.info");
    misc_tree_->Branch("photonTriggerNames",&photonTrigNamesV_);
    misc_tree_->Branch("jetTriggerNames",&jetTrigNamesV_);
    // put time stamp
    time_t ltime;
    ltime=time(NULL);
    TString str = TString(asctime(localtime(&ltime)));
    if (str[str.Length()-1]=='\n') str.Remove(str.Length()-1,1);
    TObjString date(str);
    date.Write(str.Data());
    rootfile_->cd();
  }

  // create the trees for the calo/pf jets
  if (doCaloJets_) {
    calo_tree_ = new TTree("calo_gammajettree", "tree for gamma+jet balancing using CaloJets");
    assert(calo_tree_);
  }
  if (doPFJets_) {
    pf_tree_ = new TTree("pf_gammajettree", "tree for gamma+jet balancing using PFJets");
    assert(pf_tree_);
  }

  //
  // Photon and event info. Duplicate into the trees
  //
  for (int iJet=0; iJet<2; iJet++) {
    bool doJet=(iJet==0) ? doCaloJets_ : doPFJets_;
    if (!doJet) continue;
    TTree *tree= (iJet==0) ? calo_tree_ : pf_tree_;

    // Event triggers
    tree->Branch("photonTrig_fired", &photonTrigFired_);
    tree->Branch("photonTrig_prescale", &photonTrigPrescale_);
    tree->Branch("jetTrig_fired", &jetTrigFired_);
    tree->Branch("jetTrig_prescale", &jetTrigPrescale_);

    // Event info
    tree->Branch("RunNumber",&runNumber_, "RunNumber/I");
    tree->Branch("LumiBlock",&lumiBlock_, "LumiBlock/I");
    tree->Branch("EventNumber",&eventNumber_, "EventNumber/I");
    tree->Branch("EventWeight",&eventWeight_, "EventWeight/F");

    // Photon info
    tree->Branch("rho2012", &rho2012_, "rho2012/F");
    tree->Branch("tagPho_pt",    &tagPho_pt_,    "tagPho_pt/F");
    tree->Branch("pho_2nd_pt",   &pho_2nd_pt_,   "pho_2nd_pt/F");
    tree->Branch("tagPho_energy",     &tagPho_energy_,     "tagPho_energy/F");
    tree->Branch("tagPho_eta",   &tagPho_eta_,   "tagPho_eta/F");
    tree->Branch("tagPho_phi",   &tagPho_phi_,   "tagPho_phi/F");
    tree->Branch("tagPho_sieie", &tagPho_sieie_, "tagPho_sieie/F");
    tree->Branch("tagPho_HoE",   &tagPho_HoE_,   "tagPho_HoE/F");
    tree->Branch("tagPho_r9",    &tagPho_r9_,    "tagPho_r9/F");
    tree->Branch("tagPho_EcalIsoDR04",&tagPho_EcalIsoDR04_, "tagPho_EcalIsoDR04/F");
    tree->Branch("tagPho_HcalIsoDR04",&tagPho_HcalIsoDR04_, "tagPho_HcalIsoDR04/F");
    tree->Branch("tagPho_HcalIsoDR0412",&tagPho_HcalIsoDR0412_, "tagPho_HcalIsoDR0412/F");
    tree->Branch("tagPho_TrkIsoHollowDR04",&tagPho_TrkIsoHollowDR04_, "tagPho_TrkIsoHollowDR04/F");
    tree->Branch("tagPho_pfiso_myphoton03",&tagPho_pfiso_myphoton03_, "tagPho_pfiso_myphoton03/F");
    tree->Branch("tagPho_pfiso_myneutral03",&tagPho_pfiso_myneutral03_, "tagPho_pfiso_myneutral03/F");
    tree->Branch("tagPho_pfiso_mycharged03","std::vector<std::vector<float> >", &tagPho_pfiso_mycharged03);
    tree->Branch("tagPho_pixelSeed",    &tagPho_pixelSeed_,    "tagPho_pixelSeed/I");
    tree->Branch("tagPho_ConvSafeEleVeto", &tagPho_ConvSafeEleVeto_, "tagPho_ConvSafeEleVeto/I");
    tree->Branch("tagPho_idTight",&tagPho_idTight_, "tagPho_idTight/I");
    tree->Branch("tagPho_idLoose",&tagPho_idLoose_, "tagPho_idLoose/I");
    // gen.info on photon
    if(doGenJets_){
      tree->Branch("tagPho_genPt",&tagPho_genPt_, "tagPho_genPt/F");
      tree->Branch("tagPho_genEnergy",&tagPho_genEnergy_,"tagPho_genEnergy/F");
      tree->Branch("tagPho_genEta",&tagPho_genEta_, "tagPho_genEta/F");
      tree->Branch("tagPho_genPhi",&tagPho_genPhi_, "tagPho_genPhi/F");
      tree->Branch("tagPho_genDeltaR",&tagPho_genDeltaR_,"tagPho_genDeltaR/F");
    }
      // counters
    tree->Branch("nPhotons",&nPhotons_, "nPhotons/I");
    tree->Branch("nGenJets",&nGenJets_, "nGenJets/I");
  }

  //
  // Calo jets
  //

  if(doCaloJets_){
    //hPassSelCalo_ = new TH1D("hPassSelectionCalo", "Selection Pass Failures CaloJets",200,-0.5,199.5);

    calo_tree_->Branch("nCaloJets",&nCaloJets_, "nCaloJets/I");

    calo_tree_->Branch("tcalojet_et",&tcalojet_et_, "tcalojet_et/F");
    calo_tree_->Branch("tcalojet_pt",&tcalojet_pt_, "tcalojet_pt/F");
    calo_tree_->Branch("tcalojet_p",&tcalojet_p_, "tcalojet_p/F");
    calo_tree_->Branch("tcalojet_eta",&tcalojet_eta_, "tcalojet_eta/F");
    calo_tree_->Branch("tcalojet_phi",&tcalojet_phi_, "tcalojet_phi/F");
    calo_tree_->Branch("tcalojet_emf",&tcalojet_emf_, "tcalojet_emf/F");
    calo_tree_->Branch("tcalojet_scale",&tcalojet_scale_, "tcalojet_scale/F");
    if(doGenJets_){
      calo_tree_->Branch("tcalojet_genpt",&tcalojet_genpt_, "tcalojet_genpt/F");
      calo_tree_->Branch("tcalojet_genp",&tcalojet_genp_, "tcalojet_genp/F");
      calo_tree_->Branch("tcalojet_gendr",&tcalojet_gendr_, "tcalojet_gendr/F");
    }
    calo_tree_->Branch("tcalojet_EBE",&tcalojet_EBE_, "tcalojet_EBE/F");
    calo_tree_->Branch("tcalojet_EEE",&tcalojet_EEE_, "tcalojet_EEE/F");
    calo_tree_->Branch("tcalojet_HBE",&tcalojet_HBE_, "tcalojet_HBE/F");
    calo_tree_->Branch("tcalojet_HEE",&tcalojet_HEE_, "tcalojet_HEE/F");
    calo_tree_->Branch("tcalojet_HFEtot",&tcalojet_HFEtot_, "tcalojet_HFEtot/F");
    calo_tree_->Branch("tcalojet_HFEcalE",&tcalojet_HFEcalE_, "tcalojet_HFEcalE/F");
    calo_tree_->Branch("tcalojet_HFHadE",&tcalojet_HFHadE_, "tcalojet_HFHadE/F");
    calo_tree_->Branch("tcalojet_ntwrs",&tcalojet_ntwrs_, "tcalojet_ntwrs/I");
    calo_tree_->Branch("tcalojet_twr_ieta",tcalojet_twr_ieta_, "tcalojet_twr_ieta[tcalojet_ntwrs]/I");
    calo_tree_->Branch("tcalojet_twr_totE",tcalojet_twr_totE_, "tcalojet_twr_totE[tcalojet_ntwrs]/F");
    calo_tree_->Branch("tcalojet_twr_eme",tcalojet_twr_eme_, "tcalojet_twr_eme[tcalojet_ntwrs]/F");
    calo_tree_->Branch("tcalojet_twr_hade",tcalojet_twr_hade_, "tcalojet_twr_hade[tcalojet_ntwrs]/F");
    calo_tree_->Branch("tcalojet_twr_hadoe",tcalojet_twr_hadoe_, "tcalojet_twr_hadoe[tcalojet_ntwrs]/F");

    calo_tree_->Branch("pcalojet_et",&pcalojet_et_, "pcalojet_et/F");
    calo_tree_->Branch("pcalojet_pt",&pcalojet_pt_, "pcalojet_pt/F");
    calo_tree_->Branch("pcalojet_p",&pcalojet_p_, "pcalojet_p/F");
    calo_tree_->Branch("pcalojet_eta",&pcalojet_eta_, "pcalojet_eta/F");
    calo_tree_->Branch("pcalojet_phi",&pcalojet_phi_, "pcalojet_phi/F");
    calo_tree_->Branch("pcalojet_emf",&pcalojet_emf_, "pcalojet_emf/F");
    calo_tree_->Branch("pcalojet_scale",&pcalojet_scale_, "pcalojet_scale/F");
    if(doGenJets_){
      calo_tree_->Branch("pcalojet_genpt",&pcalojet_genpt_, "pcalojet_genpt/F");
      calo_tree_->Branch("pcalojet_genp",&pcalojet_genp_, "pcalojet_genp/F");
      calo_tree_->Branch("pcalojet_gendr",&pcalojet_gendr_, "pcalojet_gendr/F");
    }
    calo_tree_->Branch("pcalojet_EBE",&pcalojet_EBE_, "pcalojet_EBE/F");
    calo_tree_->Branch("pcalojet_EEE",&pcalojet_EEE_, "pcalojet_EEE/F");
    calo_tree_->Branch("pcalojet_HBE",&pcalojet_HBE_, "pcalojet_HBE/F");
    calo_tree_->Branch("pcalojet_HEE",&pcalojet_HEE_, "pcalojet_HEE/F");
    calo_tree_->Branch("pcalojet_HFEtot",&pcalojet_HFEtot_, "pcalojet_HFEtot/F");
    calo_tree_->Branch("pcalojet_ntwrs",&pcalojet_ntwrs_, "pcalojet_ntwrs/I");
    calo_tree_->Branch("pcalojet_twr_ieta",pcalojet_twr_ieta_, "pcalojet_twr_ieta[pcalojet_ntwrs]/I");
    calo_tree_->Branch("pcalojet_twr_eme",pcalojet_twr_eme_, "pcalojet_twr_eme[pcalojet_ntwrs]/F");
    calo_tree_->Branch("pcalojet_twr_hade",pcalojet_twr_hade_, "pcalojet_twr_hade[pcalojet_ntwrs]/F");

    calo_tree_->Branch("calo_thirdjet_et", &calo_thirdjet_et_, "calo_thirdjet_et/F");
    calo_tree_->Branch("calo_thirdjet_pt", &calo_thirdjet_pt_, "calo_thirdjet_pt/F");
    calo_tree_->Branch("calo_thirdjet_p", &calo_thirdjet_p_, "calo_thirdjet_p/F");
    calo_tree_->Branch("calo_thirdjet_px", &calo_thirdjet_px_, "calo_thirdjet_px/F");
    calo_tree_->Branch("calo_thirdjet_py", &calo_thirdjet_py_, "calo_thirdjet_py/F");
    calo_tree_->Branch("calo_thirdjet_E", &calo_thirdjet_E_, "calo_thirdjet_E/F");
    calo_tree_->Branch("calo_thirdjet_eta", &calo_thirdjet_eta_, "calo_thirdjet_eta/F");
    calo_tree_->Branch("calo_thirdjet_phi", &calo_thirdjet_phi_, "calo_thirdjet_phi/F");
    calo_tree_->Branch("calo_thirdjet_scale", &calo_thirdjet_scale_, "calo_thirdjet_scale/F");
  }

    //////// Particle Flow ////////

  if (doPFJets_) {

    pf_tree_->Branch("nPFJets",&nPFJets_, "nPFJets/I");

    // Leading jet info
    pf_tree_->Branch("ppfjet_pt",&ppfjet_pt_, "ppfjet_pt/F");
    pf_tree_->Branch("ppfjet_p",&ppfjet_p_, "ppfjet_p/F");
    pf_tree_->Branch("ppfjet_E",&ppfjet_E_, "ppfjet_E/F");
    pf_tree_->Branch("ppfjet_eta",&ppfjet_eta_, "ppfjet_eta/F");
    pf_tree_->Branch("ppfjet_phi",&ppfjet_phi_, "ppfjet_phi/F");
    pf_tree_->Branch("ppfjet_scale",&ppfjet_scale_, "ppfjet_scale/F");
    pf_tree_->Branch("ppfjet_NeutralHadronFrac", &ppfjet_NeutralHadronFrac_, "ppfjet_NeutralHadronFrac/F");
    pf_tree_->Branch("ppfjet_NeutralEMFrac", &ppfjet_NeutralEMFrac_, "ppfjet_NeutralEMFrac/F");
    pf_tree_->Branch("ppfjet_nConstituents", &ppfjet_nConstituents_, "ppfjet_nConstituents/I");
    pf_tree_->Branch("ppfjet_ChargedHadronFrac", &ppfjet_ChargedHadronFrac_, "ppfjet_ChargedHadronFrac/F");
    pf_tree_->Branch("ppfjet_ChargedMultiplicity", &ppfjet_ChargedMultiplicity_, "ppfjet_ChargedMultiplicity/F");
    pf_tree_->Branch("ppfjet_ChargedEMFrac", &ppfjet_ChargedEMFrac_, "ppfjet_ChargedEMFrac/F");
    if(doGenJets_){
      pf_tree_->Branch("ppfjet_genpt",&ppfjet_genpt_, "ppfjet_genpt/F");
      pf_tree_->Branch("ppfjet_genp",&ppfjet_genp_, "ppfjet_genp/F");
      pf_tree_->Branch("ppfjet_genE",&ppfjet_genE_, "ppfjet_genE/F");
      pf_tree_->Branch("ppfjet_gendr",&ppfjet_gendr_, "ppfjet_gendr/F");
    }
    pf_tree_->Branch("ppfjet_unkown_E",&ppfjet_unkown_E_, "ppfjet_unkown_E/F");
    pf_tree_->Branch("ppfjet_electron_E",&ppfjet_electron_E_, "ppfjet_electron_E/F");
    pf_tree_->Branch("ppfjet_muon_E",&ppfjet_muon_E_, "ppfjet_muon_E/F");
    pf_tree_->Branch("ppfjet_photon_E",&ppfjet_photon_E_, "ppfjet_photon_E/F");
    pf_tree_->Branch("ppfjet_unkown_px",&ppfjet_unkown_px_, "ppfjet_unkown_px/F");
    pf_tree_->Branch("ppfjet_electron_px",&ppfjet_electron_px_, "ppfjet_electron_px/F");
    pf_tree_->Branch("ppfjet_muon_px",&ppfjet_muon_px_, "ppfjet_muon_px/F");
    pf_tree_->Branch("ppfjet_photon_px",&ppfjet_photon_px_, "ppfjet_photon_px/F");
    pf_tree_->Branch("ppfjet_unkown_py",&ppfjet_unkown_py_, "ppfjet_unkown_py/F");
    pf_tree_->Branch("ppfjet_electron_py",&ppfjet_electron_py_, "ppfjet_electron_py/F");
    pf_tree_->Branch("ppfjet_muon_py",&ppfjet_muon_py_, "ppfjet_muon_py/F");
    pf_tree_->Branch("ppfjet_photon_py",&ppfjet_photon_py_, "ppfjet_photon_py/F");
    pf_tree_->Branch("ppfjet_unkown_pz",&ppfjet_unkown_pz_, "ppfjet_unkown_pz/F");
    pf_tree_->Branch("ppfjet_electron_pz",&ppfjet_electron_pz_, "ppfjet_electron_pz/F");
    pf_tree_->Branch("ppfjet_muon_pz",&ppfjet_muon_pz_, "ppfjet_muon_pz/F");
    pf_tree_->Branch("ppfjet_photon_pz",&ppfjet_photon_pz_, "ppfjet_photon_pz/F");
    pf_tree_->Branch("ppfjet_unkown_EcalE",&ppfjet_unkown_EcalE_, "ppfjet_unkown_EcalE/F");
    pf_tree_->Branch("ppfjet_electron_EcalE",&ppfjet_electron_EcalE_, "ppfjet_electron_EcalE/F");
    pf_tree_->Branch("ppfjet_muon_EcalE",&ppfjet_muon_EcalE_, "ppfjet_muon_EcalE/F");
    pf_tree_->Branch("ppfjet_photon_EcalE",&ppfjet_photon_EcalE_, "ppfjet_photon_EcalE/F");
    pf_tree_->Branch("ppfjet_unkown_n",&ppfjet_unkown_n_, "ppfjet_unkown_n/I");
    pf_tree_->Branch("ppfjet_electron_n",&ppfjet_electron_n_, "ppfjet_electron_n/I");
    pf_tree_->Branch("ppfjet_muon_n",&ppfjet_muon_n_, "ppfjet_muon_n/I");
    pf_tree_->Branch("ppfjet_photon_n",&ppfjet_photon_n_, "ppfjet_photon_n/I");
    pf_tree_->Branch("ppfjet_had_n",&ppfjet_had_n_, "ppfjet_had_n/I");
    pf_tree_->Branch("ppfjet_had_E",&ppfjet_had_E_);
    pf_tree_->Branch("ppfjet_had_px",&ppfjet_had_px_);
    pf_tree_->Branch("ppfjet_had_py",&ppfjet_had_py_);
    pf_tree_->Branch("ppfjet_had_pz",&ppfjet_had_pz_);
    pf_tree_->Branch("ppfjet_had_EcalE",&ppfjet_had_EcalE_);
    pf_tree_->Branch("ppfjet_had_rawHcalE",&ppfjet_had_rawHcalE_);
    pf_tree_->Branch("ppfjet_had_emf",&ppfjet_had_emf_);
    pf_tree_->Branch("ppfjet_had_id",&ppfjet_had_id_);
    pf_tree_->Branch("ppfjet_had_candtrackind",&ppfjet_had_candtrackind_);
    if(doGenJets_){
      pf_tree_->Branch("ppfjet_had_E_mctruth",&ppfjet_had_E_mctruth_);
      pf_tree_->Branch("ppfjet_had_mcpdgId",&ppfjet_had_mcpdgId_);
    }
    pf_tree_->Branch("ppfjet_had_ntwrs",&ppfjet_had_ntwrs_);
    pf_tree_->Branch("ppfjet_ntwrs",&ppfjet_ntwrs_, "ppfjet_ntwrs/I");
    pf_tree_->Branch("ppfjet_twr_ieta",&ppfjet_twr_ieta_);
    pf_tree_->Branch("ppfjet_twr_iphi",&ppfjet_twr_iphi_);
    pf_tree_->Branch("ppfjet_twr_depth",&ppfjet_twr_depth_);
    pf_tree_->Branch("ppfjet_twr_subdet",&ppfjet_twr_subdet_);
    pf_tree_->Branch("ppfjet_twr_hade",&ppfjet_twr_hade_);
    pf_tree_->Branch("ppfjet_twr_frac",&ppfjet_twr_frac_);
    pf_tree_->Branch("ppfjet_twr_candtrackind",&ppfjet_twr_candtrackind_);
    pf_tree_->Branch("ppfjet_twr_hadind",&ppfjet_twr_hadind_);
    pf_tree_->Branch("ppfjet_twr_elmttype",&ppfjet_twr_elmttype_);
    pf_tree_->Branch("ppfjet_twr_dR",&ppfjet_twr_dR_);
    pf_tree_->Branch("ppfjet_twr_clusterind",&ppfjet_twr_clusterind_);
    pf_tree_->Branch("ppfjet_cluster_n",&ppfjet_cluster_n_, "ppfjet_cluster_n/I");
    pf_tree_->Branch("ppfjet_cluster_eta",&ppfjet_cluster_eta_);
    pf_tree_->Branch("ppfjet_cluster_phi",&ppfjet_cluster_phi_);
    pf_tree_->Branch("ppfjet_cluster_dR",&ppfjet_cluster_dR_);
    pf_tree_->Branch("ppfjet_ncandtracks",&ppfjet_ncandtracks_, "ppfjet_ncandtracks/I");
    pf_tree_->Branch("ppfjet_candtrack_px",&ppfjet_candtrack_px_);
    pf_tree_->Branch("ppfjet_candtrack_py",&ppfjet_candtrack_py_);
    pf_tree_->Branch("ppfjet_candtrack_pz",&ppfjet_candtrack_pz_);
    pf_tree_->Branch("ppfjet_candtrack_EcalE",&ppfjet_candtrack_EcalE_);

    // Subleading jet info
    pf_tree_->Branch("pfjet2_pt",&pfjet2_pt_, "pfjet2_pt/F");
    pf_tree_->Branch("pfjet2_p",&pfjet2_p_, "pfjet2_p/F");
    pf_tree_->Branch("pfjet2_E",&pfjet2_E_, "pfjet2_E/F");
    pf_tree_->Branch("pfjet2_eta",&pfjet2_eta_, "pfjet2_eta/F");
    pf_tree_->Branch("pfjet2_phi",&pfjet2_phi_, "pfjet2_phi/F");
    pf_tree_->Branch("pfjet2_scale",&pfjet2_scale_, "pfjet2_scale/F");
    pf_tree_->Branch("pfjet2_NeutralHadronFrac", &pfjet2_NeutralHadronFrac_, "pfjet2_NeutralHadronFrac/F");
    pf_tree_->Branch("pfjet2_NeutralEMFrac", &pfjet2_NeutralEMFrac_, "pfjet2_NeutralEMFrac/F");
    pf_tree_->Branch("pfjet2_nConstituents", &pfjet2_nConstituents_, "pfjet2_nConstituents/I");
    pf_tree_->Branch("pfjet2_ChargedHadronFrac", &pfjet2_ChargedHadronFrac_, "pfjet2_ChargedHadronFrac/F");
    pf_tree_->Branch("pfjet2_ChargedMultiplicity", &pfjet2_ChargedMultiplicity_, "pfjet2_ChargedMultiplicity/F");
    pf_tree_->Branch("pfjet2_ChargedEMFrac", &pfjet2_ChargedEMFrac_, "pfjet2_ChargedEMFrac/F");
    if(doGenJets_){
      pf_tree_->Branch("pfjet2_genpt",&pfjet2_genpt_, "pfjet2_genpt/F");
      pf_tree_->Branch("pfjet2_genp",&pfjet2_genp_, "pfjet2_genp/F");
      pf_tree_->Branch("pfjet2_genE",&pfjet2_genE_, "pfjet2_genE/F");
      pf_tree_->Branch("pfjet2_gendr",&pfjet2_gendr_, "pfjet2_gendr/F");
    }
    pf_tree_->Branch("pfjet2_unkown_E",&pfjet2_unkown_E_, "pfjet2_unkown_E/F");
    pf_tree_->Branch("pfjet2_electron_E",&pfjet2_electron_E_, "pfjet2_electron_E/F");
    pf_tree_->Branch("pfjet2_muon_E",&pfjet2_muon_E_, "pfjet2_muon_E/F");
    pf_tree_->Branch("pfjet2_photon_E",&pfjet2_photon_E_, "pfjet2_photon_E/F");
    pf_tree_->Branch("pfjet2_unkown_px",&pfjet2_unkown_px_, "pfjet2_unkown_px/F");
    pf_tree_->Branch("pfjet2_electron_px",&pfjet2_electron_px_, "pfjet2_electron_px/F");
    pf_tree_->Branch("pfjet2_muon_px",&pfjet2_muon_px_, "pfjet2_muon_px/F");
    pf_tree_->Branch("pfjet2_photon_px",&pfjet2_photon_px_, "pfjet2_photon_px/F");
    pf_tree_->Branch("pfjet2_unkown_py",&pfjet2_unkown_py_, "pfjet2_unkown_py/F");
    pf_tree_->Branch("pfjet2_electron_py",&pfjet2_electron_py_, "pfjet2_electron_py/F");
    pf_tree_->Branch("pfjet2_muon_py",&pfjet2_muon_py_, "pfjet2_muon_py/F");
    pf_tree_->Branch("pfjet2_photon_py",&pfjet2_photon_py_, "pfjet2_photon_py/F");
    pf_tree_->Branch("pfjet2_unkown_pz",&pfjet2_unkown_pz_, "pfjet2_unkown_pz/F");
    pf_tree_->Branch("pfjet2_electron_pz",&pfjet2_electron_pz_, "pfjet2_electron_pz/F");
    pf_tree_->Branch("pfjet2_muon_pz",&pfjet2_muon_pz_, "pfjet2_muon_pz/F");
    pf_tree_->Branch("pfjet2_photon_pz",&pfjet2_photon_pz_, "pfjet2_photon_pz/F");
    pf_tree_->Branch("pfjet2_unkown_EcalE",&pfjet2_unkown_EcalE_, "pfjet2_unkown_EcalE/F");
    pf_tree_->Branch("pfjet2_electron_EcalE",&pfjet2_electron_EcalE_, "pfjet2_electron_EcalE/F");
    pf_tree_->Branch("pfjet2_muon_EcalE",&pfjet2_muon_EcalE_, "pfjet2_muon_EcalE/F");
    pf_tree_->Branch("pfjet2_photon_EcalE",&pfjet2_photon_EcalE_, "pfjet2_photon_EcalE/F");
    pf_tree_->Branch("pfjet2_unkown_n",&pfjet2_unkown_n_, "pfjet2_unkown_n/I");
    pf_tree_->Branch("pfjet2_electron_n",&pfjet2_electron_n_, "pfjet2_electron_n/I");
    pf_tree_->Branch("pfjet2_muon_n",&pfjet2_muon_n_, "pfjet2_muon_n/I");
    pf_tree_->Branch("pfjet2_photon_n",&pfjet2_photon_n_, "pfjet2_photon_n/I");
    pf_tree_->Branch("pfjet2_had_n",&pfjet2_had_n_, "pfjet2_had_n/I");
    pf_tree_->Branch("pfjet2_had_E",&pfjet2_had_E_);
    pf_tree_->Branch("pfjet2_had_px",&pfjet2_had_px_);
    pf_tree_->Branch("pfjet2_had_py",&pfjet2_had_py_);
    pf_tree_->Branch("pfjet2_had_pz",&pfjet2_had_pz_);
    pf_tree_->Branch("pfjet2_had_EcalE",&pfjet2_had_EcalE_);
    pf_tree_->Branch("pfjet2_had_rawHcalE",&pfjet2_had_rawHcalE_);
    pf_tree_->Branch("pfjet2_had_emf",&pfjet2_had_emf_);
    pf_tree_->Branch("pfjet2_had_id",&pfjet2_had_id_);
    pf_tree_->Branch("pfjet2_had_candtrackind",&pfjet2_had_candtrackind_);
    if(doGenJets_){
      pf_tree_->Branch("pfjet2_had_E_mctruth",&pfjet2_had_E_mctruth_);
      pf_tree_->Branch("pfjet2_had_mcpdgId",&pfjet2_had_mcpdgId_);
    }
    pf_tree_->Branch("pfjet2_had_ntwrs",&pfjet2_had_ntwrs_);
    pf_tree_->Branch("pfjet2_ntwrs",&pfjet2_ntwrs_, "pfjet2_ntwrs/I");
    pf_tree_->Branch("pfjet2_twr_ieta",&pfjet2_twr_ieta_);
    pf_tree_->Branch("pfjet2_twr_iphi",&pfjet2_twr_iphi_);
    pf_tree_->Branch("pfjet2_twr_depth",&pfjet2_twr_depth_);
    pf_tree_->Branch("pfjet2_twr_subdet",&pfjet2_twr_subdet_);
    pf_tree_->Branch("pfjet2_twr_hade",&pfjet2_twr_hade_);
    pf_tree_->Branch("pfjet2_twr_frac",&pfjet2_twr_frac_);
    pf_tree_->Branch("pfjet2_twr_candtrackind",&pfjet2_twr_candtrackind_);
    pf_tree_->Branch("pfjet2_twr_hadind",&pfjet2_twr_hadind_);
    pf_tree_->Branch("pfjet2_twr_elmttype",&pfjet2_twr_elmttype_);
    pf_tree_->Branch("pfjet2_twr_dR",&pfjet2_twr_dR_);
    pf_tree_->Branch("pfjet2_twr_clusterind",&pfjet2_twr_clusterind_);
    pf_tree_->Branch("pfjet2_cluster_n",&pfjet2_cluster_n_, "pfjet2_cluster_n/I");
    pf_tree_->Branch("pfjet2_cluster_eta",&pfjet2_cluster_eta_);
    pf_tree_->Branch("pfjet2_cluster_phi",&pfjet2_cluster_phi_);
    pf_tree_->Branch("pfjet2_cluster_dR",&pfjet2_cluster_dR_);
    pf_tree_->Branch("pfjet2_ncandtracks",&pfjet2_ncandtracks_, "pfjet2_ncandtracks/I");
    pf_tree_->Branch("pfjet2_candtrack_px",&pfjet2_candtrack_px_);
    pf_tree_->Branch("pfjet2_candtrack_py",&pfjet2_candtrack_py_);
    pf_tree_->Branch("pfjet2_candtrack_pz",&pfjet2_candtrack_pz_);
    pf_tree_->Branch("pfjet2_candtrack_EcalE",&pfjet2_candtrack_EcalE_);

    // third pf jet
    pf_tree_->Branch("pf_thirdjet_et", &pf_thirdjet_et_, "pf_thirdjet_et/F");
    pf_tree_->Branch("pf_thirdjet_pt", &pf_thirdjet_pt_, "pf_thirdjet_pt/F");
    pf_tree_->Branch("pf_thirdjet_p", &pf_thirdjet_p_, "pf_thirdjet_p/F");
    pf_tree_->Branch("pf_thirdjet_px", &pf_thirdjet_px_, "pf_thirdjet_px/F");
    pf_tree_->Branch("pf_thirdjet_py", &pf_thirdjet_py_, "pf_thirdjet_py/F");
    pf_tree_->Branch("pf_thirdjet_E", &pf_thirdjet_E_, "pf_thirdjet_E/F");
    pf_tree_->Branch("pf_thirdjet_eta", &pf_thirdjet_eta_, "pf_thirdjet_eta/F");
    pf_tree_->Branch("pf_thirdjet_phi", &pf_thirdjet_phi_, "pf_thirdjet_phi/F");
    pf_tree_->Branch("pf_thirdjet_scale", &pf_thirdjet_scale_, "pf_thirdjet_scale/F");
  }

  return;
  ////  std::cout << "End beginJob()" << std::endl;
}  

// ------------ method called once each job just after ending the event loop  ------------
void 
CalcRespCorrPhotonPlusJet::endJob() {
  ///  std::cout << "Start endJob()" << std::endl;

  // write miscItems
  rootfile_->cd();
  rootfile_->cd("miscItems");
  misc_tree_->Fill();
  misc_tree_->Write();

  rootfile_->cd();

  if (doCaloJets_) {
    calo_tree_->Write();
  }

  if(doPFJets_){
    // write histograms
    h_types_->Write();
    h_ntypes_->Write();
    h_ietaHCAL_->Write();
    h_etaHFHAD_->Write();
    h_etaHFEM_->Write();
    h_ietaHO_->Write();
    h_HFHAD_n_->Write();
    h_HFEM_n_->Write();
    h_HFHAD_type_->Write();
    h_HFEM_type_->Write();
    h_HBHE_n_->Write();
    h_HF_n_->Write();
    h_HO_n_->Write();
    h_twrietas_->Write();
    h_rechitspos_->Write();
    h_hbherecoieta_->Write();
    pf_tree_->Write();
  }
  rootfile_->Close();
  ////  std::cout << "End endJob()" << std::endl;
}


// ---------------------------------------------------------------------

void CalcRespCorrPhotonPlusJet::beginRun(const edm::Run &iRun,
					 const edm::EventSetup &setup)
{
  //std::cout << "beginRun()" << std::endl;
  if (debug_) std::cout <<"Initializing trigger information for individual run"<<std::endl;
  bool changed(true);
  std::string processName="HLT";
  if (hltConfig_.init(iRun,setup,processName,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
     // The HLT config has actually changed wrt the previous Run, hence rebook your
     // histograms or do anything else dependent on the revised HLT config
    }
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    throw edm::Exception(edm::errors::ProductNotFound)
      << " HLT config extraction failure with process name " << processName;
    // In this case, all access methods will return empty values!
  }
}

// ---------------------------------------------------------------------

// helper function

float CalcRespCorrPhotonPlusJet::pfEcalIso(const reco::Photon* localPho1, edm::Handle<reco::PFCandidateCollection> pfHandle, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel, float etaStripEndcap, float energyBarrel, float energyEndcap, reco::PFCandidate::ParticleType pfToUse) {
  ////std::cout << "Inside pfEcalIso" << std::endl;
  reco::Photon* localPho = localPho1->clone();
  float dRVeto;
  float etaStrip;

  if (localPho->isEB()) {
    dRVeto = dRVetoBarrel;
    etaStrip = etaStripBarrel;
  } else {
    dRVeto = dRVetoEndcap;
    etaStrip = etaStripEndcap;
  }
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();
  int nsize = forIsolation->size();
  float sum = 0;
  for (int i=0; i<nsize; i++) {
    const reco::PFCandidate& pfc = (*forIsolation)[i];
    if (pfc.particleId() ==  pfToUse) {
      // Do not include the PFCandidate associated by SC Ref to the reco::Photon                                                       
      if(pfc.superClusterRef().isNonnull() && localPho->superCluster().isNonnull()) {
        if (pfc.superClusterRef() == localPho->superCluster())
          continue;
      }

      if (localPho->isEB()) {
        if (fabs(pfc.pt()) < energyBarrel)
          continue;
      } else {
        if (fabs(pfc.energy()) < energyEndcap)
          continue;
      }
      // Shift the photon direction vector according to the PF vertex                                                                  
      math::XYZPoint pfvtx = pfc.vertex();
      math::XYZVector photon_directionWrtVtx(localPho->superCluster()->x() - pfvtx.x(),
                                             localPho->superCluster()->y() - pfvtx.y(),
                                             localPho->superCluster()->z() - pfvtx.z());

      float dEta = fabs(photon_directionWrtVtx.Eta() - pfc.momentum().Eta());
      float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(), pfc.momentum().Eta(), pfc.momentum().Phi());

      if (dEta < etaStrip)
        continue;

      if(dR > dRmax || dR < dRVeto)
        continue;

      sum += pfc.pt();
    }
  }
  return sum;
}

// ---------------------------------------------------------------------

float CalcRespCorrPhotonPlusJet::pfHcalIso(const reco::Photon* localPho,edm::Handle<reco::PFCandidateCollection> pfHandle,float dRmax, float dRveto,reco::PFCandidate::ParticleType pfToUse) {
  //// std::cout << "Inside pfHcalIso" << std::endl;
  return pfEcalIso(localPho, pfHandle, dRmax, dRveto, dRveto, 0.0, 0.0, 0.0, 0.0, pfToUse);

}

// ---------------------------------------------------------------------

std::vector<float> CalcRespCorrPhotonPlusJet::pfTkIsoWithVertex(const reco::Photon* localPho1, edm::Handle<reco::PFCandidateCollection> pfHandle, edm::Handle<reco::VertexCollection> vtxHandle, float dRmax, float dRvetoBarrel, float dRvetoEndcap, float ptMin, float dzMax, float dxyMax, reco::PFCandidate::ParticleType pfToUse) {

  //  std::cout << "Inside pfTkIsoWithVertex()" << std::endl;
  reco::Photon* localPho = localPho1->clone();

  float dRveto;
  if (localPho->isEB())
    dRveto = dRvetoBarrel;
  else
    dRveto = dRvetoEndcap;

  std::vector<float> result;
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();

  //Calculate isolation sum separately for each vertex
  //  std::cout << "vtxHandle->size() = " << vtxHandle->size() << std::endl;
  for(unsigned int ivtx=0; ivtx<(vtxHandle->size()); ++ivtx) {
    //std::cout << "Vtx " << ivtx << std::endl;
    // Shift the photon according to the vertex                                                                                        
    reco::VertexRef vtx(vtxHandle, ivtx);
    math::XYZVector photon_directionWrtVtx(localPho->superCluster()->x() - vtx->x(),
                                           localPho->superCluster()->y() - vtx->y(),
                                           localPho->superCluster()->z() - vtx->z());
    //std::cout << "pfTkIsoWithVertex :: Will Loop over the PFCandidates" << std::endl;
    float sum = 0;
    // Loop over the PFCandidates                                                                                                      
    for(unsigned i=0; i<forIsolation->size(); i++) {
      //    std::cout << "inside loop" << std::endl; 
      const reco::PFCandidate& pfc = (*forIsolation)[i];

      //require that PFCandidate is a charged hadron
      // std::cout << "pfToUse=" << pfToUse << std::endl;
      //  std::cout<< "pfc.particleId()=" << pfc.particleId() << std::endl;

      if (pfc.particleId() == pfToUse) {
	//std::cout << "\n ***** HERE pfc.particleId() == pfToUse " << std::endl;
	//std::cout << "pfc.pt()=" << pfc.pt() << std::endl;
        if (pfc.pt() < ptMin)
          continue;

        float dz = fabs(pfc.trackRef()->dz(vtx->position()));
        if (dz > dzMax) continue;

        float dxy = fabs(pfc.trackRef()->dxy(vtx->position()));
        if(fabs(dxy) > dxyMax) continue;
        float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
        if(dR > dRmax || dR < dRveto) continue;
        sum += pfc.pt();
	//	std::cout << "pt=" << pfc.pt() << std::endl;
      }
    }
    //    std::cout << "sum=" << sum << std::endl;
    sum = sum*1.0;
    result.push_back(sum);
  }
  //  std::cout << "Will return result" << std::endl;
  // std::cout << "result" << &result << std::endl;
  return result;
  //std::cout << "Result returned" << std::endl;
}

// ---------------------------------------------------------------------

void CalcRespCorrPhotonPlusJet::clear_leadingPfJetVars() {
  ppfjet_pt_ = ppfjet_p_ = ppfjet_E_ = 0;
  ppfjet_eta_ = ppfjet_phi_ = ppfjet_scale_ = 0.;
  ppfjet_NeutralHadronFrac_ = ppfjet_NeutralEMFrac_ = 0.;
  ppfjet_nConstituents_ = 0;
  ppfjet_ChargedHadronFrac_ = ppfjet_ChargedMultiplicity_ = 0;
  ppfjet_ChargedEMFrac_ = 0.;
  ppfjet_gendr_ = ppfjet_genpt_ = ppfjet_genp_ = ppfjet_genE_ = 0.;
  // Reset particle variables
  ppfjet_unkown_E_ = ppfjet_unkown_px_ = ppfjet_unkown_py_ = ppfjet_unkown_pz_ = ppfjet_unkown_EcalE_ = 0.0;
  ppfjet_electron_E_ = ppfjet_electron_px_ = ppfjet_electron_py_ = ppfjet_electron_pz_ = ppfjet_electron_EcalE_ = 0.0;
  ppfjet_muon_E_ = ppfjet_muon_px_ = ppfjet_muon_py_ = ppfjet_muon_pz_ = ppfjet_muon_EcalE_ = 0.0;
  ppfjet_photon_E_ = ppfjet_photon_px_ = ppfjet_photon_py_ = ppfjet_photon_pz_ = ppfjet_photon_EcalE_ = 0.0;
  ppfjet_unkown_n_ = ppfjet_electron_n_ = ppfjet_muon_n_ = ppfjet_photon_n_ = 0;
  ppfjet_had_n_ = 0;
  ppfjet_ntwrs_ = 0;
  ppfjet_cluster_n_ = 0;
  ppfjet_ncandtracks_ = 0;

  ppfjet_had_E_.clear();
  ppfjet_had_px_.clear();
  ppfjet_had_py_.clear();
  ppfjet_had_pz_.clear();
  ppfjet_had_EcalE_.clear();
  ppfjet_had_rawHcalE_.clear();
  ppfjet_had_emf_.clear();
  ppfjet_had_E_mctruth_.clear();
  ppfjet_had_id_.clear();
  ppfjet_had_candtrackind_.clear();
  ppfjet_had_mcpdgId_.clear();
  ppfjet_had_ntwrs_.clear();
  ppfjet_twr_ieta_.clear();
  ppfjet_twr_iphi_.clear();
  ppfjet_twr_depth_.clear();
  ppfjet_twr_subdet_.clear();
  ppfjet_twr_candtrackind_.clear();
  ppfjet_twr_hadind_.clear();
  ppfjet_twr_elmttype_.clear();
  ppfjet_twr_hade_.clear();
  ppfjet_twr_frac_.clear();
  ppfjet_twr_dR_.clear();
  ppfjet_twr_clusterind_.clear();
  ppfjet_cluster_eta_.clear();
  ppfjet_cluster_phi_.clear();
  ppfjet_cluster_dR_.clear();
  ppfjet_candtrack_px_.clear();
  ppfjet_candtrack_py_.clear();
  ppfjet_candtrack_pz_.clear();
  ppfjet_candtrack_EcalE_.clear();

}

// ---------------------------------------------------------------------

void CalcRespCorrPhotonPlusJet::copy_leadingPfJetVars_to_pfJet2() {
  pfjet2_pt_ = ppfjet_pt_;
  pfjet2_p_ = ppfjet_p_;
  pfjet2_E_ = ppfjet_E_;
  pfjet2_eta_ = ppfjet_eta_;
  pfjet2_phi_ = ppfjet_phi_;
  pfjet2_scale_ = ppfjet_scale_;
  pfjet2_NeutralHadronFrac_ = ppfjet_NeutralHadronFrac_;
  pfjet2_NeutralEMFrac_ = ppfjet_NeutralEMFrac_;
  pfjet2_nConstituents_ = ppfjet_nConstituents_;
  pfjet2_ChargedHadronFrac_ = ppfjet_ChargedHadronFrac_;
  pfjet2_ChargedMultiplicity_ = ppfjet_ChargedMultiplicity_;
  pfjet2_ChargedEMFrac_ = ppfjet_ChargedEMFrac_;

  pfjet2_gendr_ = ppfjet_gendr_;
  pfjet2_genpt_ = ppfjet_genpt_;
  pfjet2_genp_ = ppfjet_genp_;
  pfjet2_genE_ = ppfjet_genE_;

  pfjet2_unkown_E_ = ppfjet_unkown_E_;
  pfjet2_unkown_px_ = ppfjet_unkown_px_;
  pfjet2_unkown_py_ = ppfjet_unkown_py_;
  pfjet2_unkown_pz_ = ppfjet_unkown_pz_;
  pfjet2_unkown_EcalE_ = ppfjet_unkown_EcalE_;

  pfjet2_electron_E_ = ppfjet_electron_E_;
  pfjet2_electron_px_ = ppfjet_electron_px_;
  pfjet2_electron_py_ = ppfjet_electron_py_;
  pfjet2_electron_pz_ = ppfjet_electron_pz_;
  pfjet2_electron_EcalE_ = ppfjet_electron_EcalE_;

  pfjet2_muon_E_ = ppfjet_muon_E_;
  pfjet2_muon_px_ = ppfjet_muon_px_;
  pfjet2_muon_py_ = ppfjet_muon_py_;
  pfjet2_muon_pz_ = ppfjet_muon_pz_;
  pfjet2_muon_EcalE_ = ppfjet_muon_EcalE_;

  pfjet2_photon_E_ = ppfjet_photon_E_;
  pfjet2_photon_px_ = ppfjet_photon_px_;
  pfjet2_photon_py_ = ppfjet_photon_py_;
  pfjet2_photon_pz_ = ppfjet_photon_pz_;
  pfjet2_photon_EcalE_ = ppfjet_photon_EcalE_;

  pfjet2_unkown_n_ = ppfjet_unkown_n_;
  pfjet2_electron_n_ = ppfjet_electron_n_;
  pfjet2_muon_n_ = ppfjet_muon_n_;
  pfjet2_photon_n_ = ppfjet_photon_n_;
  pfjet2_had_n_ = ppfjet_had_n_;

  pfjet2_had_E_ = ppfjet_had_E_;
  pfjet2_had_px_ = ppfjet_had_px_;
  pfjet2_had_py_ = ppfjet_had_py_;
  pfjet2_had_pz_ = ppfjet_had_pz_;
  pfjet2_had_EcalE_ = ppfjet_had_EcalE_;
  pfjet2_had_rawHcalE_ = ppfjet_had_rawHcalE_;
  pfjet2_had_emf_ = ppfjet_had_emf_;
  pfjet2_had_E_mctruth_ = ppfjet_had_E_mctruth_;

  pfjet2_had_id_ = ppfjet_had_id_;
  pfjet2_had_candtrackind_ = ppfjet_had_candtrackind_;
  pfjet2_had_mcpdgId_ = ppfjet_had_mcpdgId_;
  pfjet2_had_ntwrs_ = ppfjet_had_ntwrs_;

  pfjet2_ntwrs_ = ppfjet_ntwrs_;
  pfjet2_twr_ieta_ = ppfjet_twr_ieta_;
  pfjet2_twr_iphi_ = ppfjet_twr_iphi_;
  pfjet2_twr_depth_ = ppfjet_twr_depth_;
  pfjet2_twr_subdet_ = ppfjet_twr_subdet_;
  pfjet2_twr_candtrackind_ = ppfjet_twr_candtrackind_;
  pfjet2_twr_hadind_ = ppfjet_twr_hadind_;
  pfjet2_twr_elmttype_ = ppfjet_twr_elmttype_;
  pfjet2_twr_clusterind_ = ppfjet_twr_clusterind_;

  pfjet2_twr_hade_ = ppfjet_twr_hade_;
  pfjet2_twr_frac_ = ppfjet_twr_frac_;
  pfjet2_twr_dR_ = ppfjet_twr_dR_;

  pfjet2_cluster_n_ = ppfjet_cluster_n_;
  pfjet2_cluster_eta_ = ppfjet_cluster_eta_;
  pfjet2_cluster_phi_ = ppfjet_cluster_phi_;
  pfjet2_cluster_dR_ = ppfjet_cluster_dR_;

  pfjet2_ncandtracks_ = ppfjet_ncandtracks_;
  pfjet2_candtrack_px_ = ppfjet_candtrack_px_;
  pfjet2_candtrack_py_ = ppfjet_candtrack_py_;
  pfjet2_candtrack_pz_ = ppfjet_candtrack_pz_;
  pfjet2_candtrack_EcalE_ = ppfjet_candtrack_EcalE_;
}

// ---------------------------------------------------------------------

double CalcRespCorrPhotonPlusJet::deltaR(const reco::Jet* j1, const reco::Jet* j2)
{
  double deta = j1->eta()-j2->eta();
  double dphi = std::fabs(j1->phi()-j2->phi());
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return std::sqrt(deta*deta + dphi*dphi);
}

// ---------------------------------------------------------------------

double CalcRespCorrPhotonPlusJet::deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  double deta = eta1 - eta2;
  double dphi = std::fabs(phi1 - phi2);
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return std::sqrt(deta*deta + dphi*dphi);
}

// ---------------------------------------------------------------------

/*
// DetId rawId bits xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//                  1111222      3333345555556666666
//   1 = detector
//   2 = subdetector
//   3 = depth
//   4 = zside: 0 = negative z, 1 = positive z \
//   5 = abs(ieta)                              | ieta,iphi
//   6 = abs(iphi)                             /
*/

// ---------------------------------------------------------------------

int CalcRespCorrPhotonPlusJet::getEtaPhi(const DetId id)
{
  return id.rawId() & 0x3FFF; // Get 14 least-significant digits
}

// ---------------------------------------------------------------------

int CalcRespCorrPhotonPlusJet::getEtaPhi(const HcalDetId id)
{
  return id.rawId() & 0x3FFF; // Get 14 least-significant digits
}

// ---------------------------------------------------------------------

//define this as a plug-in
DEFINE_FWK_MODULE(CalcRespCorrPhotonPlusJet);
