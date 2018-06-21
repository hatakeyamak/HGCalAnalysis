import FWCore.ParameterSet.Config as cms

hgcalTupleHGCDigis = cms.EDProducer("HGCalTupleMaker_HGCDigis",
  source = cms.untracked.VInputTag(
        cms.untracked.InputTag("mix","HGCDigisEE"),
        cms.untracked.InputTag("mix","HGCDigisHEfront")
      #cms.untracked.InputTag("mix","HGCDigisHEback")
        ),
  geometrySource = cms.untracked.vstring(
        'HGCalEESensitive',
        'HGCalHESiliconSensitive'
      #'HCal'
  ),
  Prefix = cms.untracked.string  ("HGCDigi"),
  Suffix = cms.untracked.string  ("")
)

