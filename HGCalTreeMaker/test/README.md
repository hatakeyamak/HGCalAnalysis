
## HGCAL tuple maker (normal mode)

```
cmsRun run_HGCalTupleMaker_2023.py
```

## HGCAL sequence:

Kevin's March 1, 2018 talk:
Link: 

### SimHits:
HGCHitsEE
HGCHitsHEfront
HcalHits
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsEE"       "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEfront"   "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEback"   "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HcalHits"        "SIM"

### Digis:
https://github.com/cms-sw/cmssw/blob/master/SimCalorimetry/HGCalSimProducers/src/HGCDigitizer.cc
https://github.com/cms-sw/cmssw/blob/master/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py
HGCDigisEE
HGCDigisHEfront
HGCDigisHEback
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisEE"      "HLT"     
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisHEfront"   "HLT"     
edm::SortedCollection<HGCDataFrame<HcalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HcalDetId,HGCSample> > >    "mix"                       "HGCDigisHEback"   "HLT"     

### RecHits (uncalib):
https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.cc
https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCEEUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEBUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEFUncalibRecHits"   "RECO"    

### RedHits (calib):
https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitProducer.cc
https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCEERecHits"    "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEBRecHits"   "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEFRecHits"   "RECO"    


