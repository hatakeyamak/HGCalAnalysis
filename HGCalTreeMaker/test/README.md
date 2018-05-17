
## HGCAL tuple maker 

### Normal mode
```
cmsRun run_HGCalTupleMaker_2023.py
```

### Full(?) content mode (including uncalibrated rechits, digis)
```
cmsRun run_HGCalTupleMaker_2023_full.py
```

- - - -

## HGCAL sequence:

Kevin's March 1, 2018 talk:
Link: 

### SimHits:
```
HGCHitsEE
HGCHitsHEfront
HcalHits
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsEE"       "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEfront"   "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEback"   "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HcalHits"        "SIM"
```

### Digis:
* https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCDigitizer.cc
* https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py
```
HGCDigisEE
HGCDigisHEfront
HGCDigisHEback
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisEE"      "HLT"     
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisHEfront"   "HLT"     
edm::SortedCollection<HGCDataFrame<HcalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HcalDetId,HGCSample> > >    "mix"                       "HGCDigisHEback"   "HLT"     
```
ADC (16 bits)
charge (fC? TBC)

### RecHits (uncalib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitWorkerWeights.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h
```
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCEEUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEFUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEBUncalibRecHits"   "RECO"    
```
ADC * adcLSB = amplitude (unit?)  
adcLSB=0.10 for HGCEE & HGCHEF  
adcLSB=0.25 for HGCHEB  

### RedHits (calib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitProducer.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py
```
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCEERecHits"    "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEFRecHits"   "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEBRecHits"   "RECO"    
```
energy (GeV)

