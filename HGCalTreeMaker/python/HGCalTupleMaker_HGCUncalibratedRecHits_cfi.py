import FWCore.ParameterSet.Config as cms

hgcalTupleHGCUncalibratedRecHits = cms.EDProducer("HGCalTupleMaker_HGCUncalibratedRecHits",
  source = cms.untracked.VInputTag(
        cms.untracked.InputTag('HGCalUncalibRecHit:HGCEEUncalibRecHits'),
        cms.untracked.InputTag('HGCalUncalibRecHit:HGCHEFUncalibRecHits'),
        cms.untracked.InputTag('HGCalUncalibRecHit:HGCHEBUncalibRecHits')
        ),
  geometrySource = cms.untracked.vstring(
        'HGCalEESensitive',
        'HGCalHESiliconSensitive',
        'HGCalHEScintillatorSensitive'
  ),
  Prefix = cms.untracked.string  ("HGCUncalibratedRecHit"),
  Suffix = cms.untracked.string  ("")
)

