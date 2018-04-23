import FWCore.ParameterSet.Config as cms

hgcalTupleGenParticles = cms.EDProducer("HGCalTupleMaker_GenParticles",
  source    = cms.untracked.InputTag("genParticles"),
  Prefix    = cms.untracked.string  ("GenPar"),
  Suffix    = cms.untracked.string  ("")
)

