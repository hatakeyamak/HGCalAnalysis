import FWCore.ParameterSet.Config as cms

hgcalTupleHBHERecHits = cms.EDProducer("HGCalTupleMaker_HBHERecHits",
  source = cms.untracked.InputTag("hbhereco"),
  Prefix = cms.untracked.string  ("HBHERecHit"),
  Suffix = cms.untracked.string  ("")
)
