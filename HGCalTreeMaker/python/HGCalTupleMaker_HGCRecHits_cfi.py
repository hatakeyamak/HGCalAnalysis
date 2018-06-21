import FWCore.ParameterSet.Config as cms

hgcalTupleHGCRecHits = cms.EDProducer("HGCalTupleMaker_HGCRecHits",
  source = cms.untracked.VInputTag(
        cms.untracked.InputTag("HGCalRecHit","HGCEERecHits"),
        cms.untracked.InputTag("HGCalRecHit","HGCHEFRecHits"),
        cms.untracked.InputTag("HGCalRecHit","HGCHEBRecHits")
        ),
  geometrySource = cms.untracked.vstring(
        'HGCalEESensitive',
        'HGCalHESiliconSensitive',
        'HGCalHEScintillatorSensitive'
  ),
  Prefix = cms.untracked.string  ("HGCRecHit"),
  Suffix = cms.untracked.string  ("")
)
