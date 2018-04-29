import FWCore.ParameterSet.Config as cms

hgcalTupleSimTracks = cms.EDProducer("HGCalTupleMaker_SimTracks",
  Source = cms.untracked.InputTag("g4SimHits",""),
  Source_SimVtx = cms.untracked.InputTag("g4SimHits",""),
  Prefix = cms.untracked.string  ("SimTracks"),
  Suffix = cms.untracked.string  ("")
)

