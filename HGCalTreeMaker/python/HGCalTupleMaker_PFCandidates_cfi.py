import FWCore.ParameterSet.Config as cms

hgcalTuplePFCandidates = cms.EDProducer("HGCalTupleMaker_PFCandidates",
  source    = cms.untracked.InputTag('particleFlow', ''),
  PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("PFPar"),
  Suffix    = cms.untracked.string  ("")
)

hgcalTuplePackedPFCandidates = hgcalTuplePFCandidates.clone(
  source    = cms.untracked.InputTag('packedPFCandidates', ''),
  PackedCandidate = cms.untracked.bool(True),
)
