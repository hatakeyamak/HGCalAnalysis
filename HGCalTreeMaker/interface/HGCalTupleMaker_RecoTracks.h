#ifndef HGCalTupleMaker_RecoTracks_h
#define HGCalTupleMaker_RecoTracks_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// muons and tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

class HGCalTupleMaker_RecoTracks : public edm::EDProducer {
 public:
  explicit HGCalTupleMaker_RecoTracks(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::EDGetTokenT<reco::TrackCollection> RecoTrackContainer_;
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  bool debug=false;
};

#endif
