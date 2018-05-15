#ifndef HGCalTupleMaker_Event_h
#define HGCalTupleMaker_Event_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HGCalTupleMaker_Event : public edm::EDProducer {
 public:
  explicit HGCalTupleMaker_Event(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );

  bool debug=false;
  
};

#endif
