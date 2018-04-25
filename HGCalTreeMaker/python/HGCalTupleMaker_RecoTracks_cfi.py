import FWCore.ParameterSet.Config as cms

hgcalTupleGeneralTracks = cms.EDProducer("HGCalTupleMaker_RecoTracks",
   source = cms.untracked.InputTag("generalTracks"),
   Prefix = cms.untracked.string("GeneralTracks"),
   Suffix = cms.untracked.string("")                                 
)                                 
