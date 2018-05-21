
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
* Additional:
* HGCEE: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCEEDigitizer.cc
* HGCHEFront: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCHEfrontDigitizer.cc
* HGC Si digitizer base: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCDigitizerBase.cc
* HGC FE electronics for Si: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCFEElectronics.cc
* HGCHEBack: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCHEbackDigitizer.cc
```
HGCDigisEE
HGCDigisHEfront
HGCDigisHEback
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisEE"      "HLT"     
edm::SortedCollection<HGCDataFrame<HGCalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HGCalDetId,HGCSample> > >    "mix"                       "HGCDigisHEfront"   "HLT"     
edm::SortedCollection<HGCDataFrame<HcalDetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<HcalDetId,HGCSample> > >    "mix"                       "HGCDigisHEback"   "HLT"     
```
```
ADC (16 bits? 12 bits? 10 bits?)
HEFront:

HEBack:  
keV->MIP->#p.e.(poisson smeared)->(xtalk&saturation)->MIP->(noise added)->chargeColl[i](MIPs)->ADC 
```

### RecHits (uncalib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py (fCPerMIP)
* Additional:
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitWorkerWeights.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalRecHitSimpleAlgo.h
```
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCEEUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEFUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEBUncalibRecHits"   "RECO"    
```
```
CEE & HEF: ADC * adcLSB / fCPerMIP_[thickness-1] -> amplitude (# of MIP equivalent)  
BH:        ADC * adcLSB                          -> amplitude (# of MIP equivalent)  
adcLSB=0.10, fC/MIP=1.25 (or 2.57 or 3.88 for some layers) for HGCEE & HGCHEF  
adcLSB=0.25 for HGCHEB  
```

### RedHits (calib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitProducer.cc (weights_[layer],rcorr_[thickness] etc)
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py
* Additional:
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitWorkerSimple.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalRecHitSimpleAlgo.h
```
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCEERecHits"    "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEFRecHits"   "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEBRecHits"   "RECO"    
```
```
uncalibRH.amplitude() * weights_[layer](dEdx, MeV, corrected for energy loss in passive materials) * 0.001 (MeV->GeV)  
 rcorr_[thickness] (0.88, 0.92, or 1) * cce_correction --> energy (GeV)
```

