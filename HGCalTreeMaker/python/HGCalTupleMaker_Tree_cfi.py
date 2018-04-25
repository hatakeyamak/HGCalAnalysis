import FWCore.ParameterSet.Config as cms

hgcalTupleTree = cms.EDAnalyzer("HGCalTupleMaker_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_hgcalTupleEvent_*_*',
        'keep *_hgcalTupleHBHERecHits_*_*',
        'keep *_hgcalTupleHORecHits_*_*',
        'keep *_hgcalTupleHFRecHits_*_*',
        'keep *_hgcalTupleHGCRecHits_*_*',
        'keep *_hgcalTupleHGC*RecHits_*_*',
        'keep *_hgcalTupleHGCSimHits_*_*',
        'keep *_hgcalTupleHGC*SimHits_*_*',
        'keep *_hgcalTupleGenParticles_*_*',
        'keep *_hgcalTupleSimTracks_*_*',
        'keep *_hgcalTupleGeneralTracks_*_*'
    )

)             
