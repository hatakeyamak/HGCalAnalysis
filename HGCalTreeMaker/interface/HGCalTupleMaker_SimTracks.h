#ifndef HGCalTupleMaker_SimTracks_h
#define HGCalTupleMaker_SimTracks_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

//
// Class definition
// 
class HGCalTupleMaker_SimTracks : public edm::EDProducer {
 protected:

  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackContainer_;

  const edm::InputTag   m_inputTag;
  const std::string     m_prefix;
  const std::string     m_suffix;

  bool debug=false;

  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup ) { 

    if (debug) std::cout << "HGCalTupleMaker_SimTracks: produce starts" << std::endl;
    
    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------

    edm::Handle<edm::SimTrackContainer> simTracks;

    iEvent.getByToken(SimTrackContainer_, simTracks);

    //
    // Loop over SimTracks
    //
    const edm::SimTrackContainer& trks = *(simTracks.product());
    edm::SimTrackContainer::const_iterator trksiter;

    if (debug) std::cout << "SimTrack Size: " << trks.size() << std::endl;

    for (trksiter = trks.begin(); trksiter != trks.end(); trksiter++){ 

      double simPt=trksiter->momentum().pt();
      double simEta=trksiter->momentum().eta();
      double simPhi=trksiter->momentum().phi();

      if (debug) std::cout << "SimTrack pt, eta, phi: " << simPt << " " << simEta << " " << simPhi << std::endl;
    
      v_pt     -> push_back ( simPt );
      v_eta    -> push_back ( simEta );
      v_phi    -> push_back ( simPhi );
      /*
      v_posx   -> push_back ( gcoord.x() );
      v_posy   -> push_back ( gcoord.y() );
      v_posz   -> push_back ( gcoord.z() );
      */

    }

    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);

  }
  
 public:
  
 HGCalTupleMaker_SimTracks(const edm::ParameterSet& iConfig) :
    m_inputTag       (iConfig.getUntrackedParameter<edm::InputTag>("Source")),
    m_prefix         (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix         (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {

    SimTrackContainer_ = consumes<edm::SimTrackContainer>(m_inputTag);

    produces<std::vector<float> > ( m_prefix + "Pt"     + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    //produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    //produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    //produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );

  }

  /*
  std::unique_ptr<std::vector<float> > v_posx;
  std::unique_ptr<std::vector<float> > v_posy;
  std::unique_ptr<std::vector<float> > v_posz;
  */
  std::unique_ptr<std::vector<float> > v_pt;
  std::unique_ptr<std::vector<float> > v_eta;
  std::unique_ptr<std::vector<float> > v_phi;

 private:

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup){
  }
  
 protected:

  void loadAlgo(){
    v_pt     = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    /*
    v_posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    */

  }
  
  void dumpAlgo( edm::Event & iEvent ){

    iEvent.put( move(v_pt     ), m_prefix + "Pt"     + m_suffix );
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );
    /*
    iEvent.put( move(v_posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(v_posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(v_posz   ), m_prefix + "Posz"   + m_suffix );
    */

  }  
    
};

#endif
