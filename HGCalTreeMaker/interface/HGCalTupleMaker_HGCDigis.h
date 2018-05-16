#ifndef HGCalTupleMaker_HGCDigis_h
#define HGCalTupleMaker_HGCDigis_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

// HGCAL objects
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

//
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// HGCAL & HCAL Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

#include "DataFormats/HcalDetId/interface/HcalTestNumbering.h"

#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"

// Used:
// http://cmslxr.fnal.gov/source/AnalysisAlgos/SiStripClusterInfoProducer/plugins/SiStripProcessedRawDigiProducer.cc
// as an example

class HGCalTupleMaker_HGCDigis : public edm::EDProducer {
 protected:

  std::vector<edm::InputTag> m_HGCDigisTags;

  edm::EDGetTokenT<HGCEEDigiCollection> m_HGCEEDigisToken;
  edm::EDGetTokenT<HGCHEDigiCollection> m_HGCHEDigisToken;
  edm::EDGetTokenT<HGCBHDigiCollection> m_HGCBHDigisToken;

  std::vector<std::string> m_geometrySource;

  const std::string     m_prefix;
  const std::string     m_suffix;

  int SampleIndx;
  
  bool debug=false;

  //HGC Geometry                                                                                                     
  std::vector<const HGCalDDDConstants*> hgcCons_;
  std::vector<const HGCalGeometry*>     hgcGeometry_;
  const HcalDDDSimConstants*            hcCons_;
  const HcalDDDRecConstants*            hcConr_;
  const CaloSubdetectorGeometry*        hcGeometry_;
  
  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup) { 

    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------
    
    edm::Handle<HGCEEDigiCollection> HGCEEDigis;
    edm::Handle<HGCHEDigiCollection> HGCHEDigis;
    edm::Handle<HGCBHDigiCollection> HGCBHDigis;

    // Loop over input tags
    for( typename std::vector<edm::InputTag>::const_iterator
	   tag = m_HGCDigisTags.begin(); tag != m_HGCDigisTags.end(); ++tag ) {
    
      unsigned index(tag - m_HGCDigisTags.begin());
      if (debug) std::cout << index << std::endl;

      std::string nameDetector_ = m_geometrySource[index];
      if (debug) std::cout << nameDetector_ << std::endl;

      //----------
      if (nameDetector_ == "HGCalEESensitive") {

    	edm::ESHandle<HGCalGeometry> geom;
    	iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom);
    	if (!geom.isValid()) {
    	  edm::LogWarning("HGCalTupleMaker_HGCDigis")
    	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
    	  return;
    	}
    	const HGCalGeometry* geom0 = geom.product();

	iEvent.getByToken(m_HGCEEDigisToken, HGCEEDigis);
	const HGCEEDigiCollection* eeDigis = HGCEEDigis.product(); // get a ptr to the product
	for(auto it = eeDigis->begin(); it != eeDigis->end(); ++it) {
	  //worker_->run1(evt, it, *eeUncalibRechits);
	  
	  HGCEEDetId detId     = (it->id());
	  int        layer     = detId.layer();
	  HGCSample  hgcSample = it->sample(SampleIndx);
	  uint16_t   gain      = hgcSample.toa();
	  uint16_t   adc       = hgcSample.data();
	  double     charge    = adc*gain;
	  fill(detId, geom0, index, layer, adc, charge);	  

	}

      }
      //----------
      else if (nameDetector_ == "HGCalHESiliconSensitive") {

    	edm::ESHandle<HGCalGeometry> geom;
    	iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom);
    	if (!geom.isValid()) {
    	  edm::LogWarning("HGCalTupleMaker_HGCDigis")
    	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
    	  return;
    	}
    	const HGCalGeometry* geom0 = geom.product();

	iEvent.getByToken(m_HGCHEDigisToken, HGCHEDigis);
	const HGCHEDigiCollection* heDigis = HGCHEDigis.product(); // get a ptr to the product
	for(auto it = heDigis->begin(); it != heDigis->end(); ++it) {
	  //worker_->run1(evt, it, *heUncalibRechits);

	   HGCHEDetId detId     = (it->id());
	   int        layer     = detId.layer();
	   HGCSample  hgcSample = it->sample(SampleIndx);
	   uint16_t   gain      = hgcSample.toa();
	   uint16_t   adc       = hgcSample.data();
	   double     charge    = adc*gain;
	   fill(detId, geom0, index, layer, adc, charge);
	}

      }
      //----------
      else if (nameDetector_ == "HCal") {

    	edm::ESHandle<CaloGeometry> geom;
    	iSetup.get<CaloGeometryRecord>().get(geom);
    	if (!geom.isValid()) {
    	  edm::LogWarning("HGCalTupleMaker_HGCDigis")
    	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
    	  return;
    	}
    	const CaloGeometry* geom0 = geom.product();

	iEvent.getByToken(m_HGCBHDigisToken, HGCBHDigis);
	const HGCBHDigiCollection* bhDigis = HGCBHDigis.product(); // get a ptr to tbh product
	for(auto it = bhDigis->begin(); it != bhDigis->end(); ++it) {
	  //worker_->run1(evt, it, *bhUncalibRechits);

	  HcalDetId  detId     = (it->id());
	  int        layer     = detId.depth();
	  HGCSample  hgcSample = it->sample(SampleIndx);
	  uint16_t   gain      = hgcSample.toa();
	  uint16_t   adc       = hgcSample.data();
	  double     charge    = adc*gain;
	  fill(detId, geom0, index, layer, adc, charge);
	}

      }
      //----------
      else {

	std::cout << "Warning!!!" << std::endl;

      }
      
    }
    
    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);

    }

 public:
  
 HGCalTupleMaker_HGCDigis(const edm::ParameterSet& iConfig) :
  m_HGCDigisTags (iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("source")),
    m_geometrySource (iConfig.getUntrackedParameter<std::vector<std::string> >("geometrySource")),
    m_prefix         (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix         (iConfig.getUntrackedParameter<std::string>  ("Suffix")),
    SampleIndx       (iConfig.getUntrackedParameter<int>("SampleIndx",5)) {
    
    m_HGCEEDigisToken = consumes<HGCEEDigiCollection>(m_HGCDigisTags[0]);
    m_HGCHEDigisToken = consumes<HGCHEDigiCollection>(m_HGCDigisTags[1]);
    m_HGCBHDigisToken = consumes<HGCBHDigiCollection>(m_HGCDigisTags[2]);
    
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Layer"  + m_suffix );
    produces<std::vector<uint16_t > > ( m_prefix + "ADC" + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Charge" + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Index"  + m_suffix );
    //produces<std::vector<int>   > ( m_prefix + "Flags"  + m_suffix );
    //produces<std::vector<int>   > ( m_prefix + "Aux"    + m_suffix );
    //produces<std::vector<float> > ( m_prefix + "Time"   + m_suffix );
  }

  std::unique_ptr<std::vector<float> > v_eta;
  std::unique_ptr<std::vector<float> > v_phi;
  std::unique_ptr<std::vector<int  > > v_layer;
  std::unique_ptr<std::vector<uint16_t > > v_adc;
  std::unique_ptr<std::vector<float> > v_charge;
  std::unique_ptr<std::vector<float> > v_posx;
  std::unique_ptr<std::vector<float> > v_posy;
  std::unique_ptr<std::vector<float> > v_posz;
  std::unique_ptr<std::vector<int>   > v_index;

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup){

    //initiating HGC Geometry
    for (size_t i=0; i<m_geometrySource.size(); i++) {
      
      // HCAL for BH/HEB
      if (m_geometrySource[i].find("HCal") != std::string::npos) {
        edm::ESHandle<HcalDDDSimConstants> pHSNDC;
        iSetup.get<HcalSimNumberingRecord>().get(pHSNDC);
        if (pHSNDC.isValid()) {
          hcCons_ = pHSNDC.product();
          hgcCons_.push_back(0);
        } else {
          edm::LogWarning("HGCalValid") << "Cannot initiate HcalDDDSimConstants: "
                                        << m_geometrySource[i] << std::endl;
        }
        edm::ESHandle<HcalDDDRecConstants> pHRNDC;
        iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
        if (pHRNDC.isValid()) {
          hcConr_ = pHRNDC.product();
        } else {
          edm::LogWarning("HGCalValid") << "Cannot initiate HcalDDDRecConstants: "
                                        << m_geometrySource[i] << std::endl;
        }
        edm::ESHandle<CaloGeometry> caloG;
        iSetup.get<CaloGeometryRecord>().get(caloG);
        if (caloG.isValid()) {
          const CaloGeometry* geo = caloG.product();
          hcGeometry_ = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
          hgcGeometry_.push_back(0);
        } else {
          edm::LogWarning("HGCalValid") << "Cannot initiate HcalGeometry for "
                                        << m_geometrySource[i] << std::endl;
        }
      }
      // HGC for EE & HEF
      else {
        edm::ESHandle<HGCalDDDConstants> hgcCons;
        iSetup.get<IdealGeometryRecord>().get(m_geometrySource[i],hgcCons);
        if (hgcCons.isValid()) {
          hgcCons_.push_back(hgcCons.product());
        } else {
          edm::LogWarning("HGCalValid") << "Cannot initiate HGCalDDDConstants for "
                                        << m_geometrySource[i] << std::endl;
        }
        edm::ESHandle<HGCalGeometry> hgcGeom;
        iSetup.get<IdealGeometryRecord>().get(m_geometrySource[i],hgcGeom);     
        if(hgcGeom.isValid()) {
          hgcGeometry_.push_back(hgcGeom.product());    
        } else {
          edm::LogWarning("HGCalValid") << "Cannot initiate HGCalGeometry for "
                                        << m_geometrySource[i] << std::endl;
        }
      }
    } // loop over geometry source
  } // beginRun


 protected:

  void loadAlgo(){
    v_eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_layer  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_adc    = std::unique_ptr<std::vector<uint16_t > > ( new std::vector<uint16_t > ());
    v_charge = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_index = std::unique_ptr<std::vector<int>   > ( new std::vector<int  > ());
  }
  
  void dumpAlgo( edm::Event & iEvent ){
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(v_layer  ), m_prefix + "Layer"  + m_suffix );
    iEvent.put( move(v_adc    ), m_prefix + "ADC"    + m_suffix );
    iEvent.put( move(v_charge ), m_prefix + "Charge" + m_suffix );
    iEvent.put( move(v_posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(v_posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(v_posz   ), m_prefix + "Posz"   + m_suffix );
    iEvent.put( move(v_index  ), m_prefix + "Index" + m_suffix );
  }  

  template<class T1, class T2>
  void fill(const T1& detId, const T2* geom, 
	    int index, int layer, uint16_t adc, double charge) {

    GlobalPoint global = geom->getPosition(detId);

    if (debug) std::cout << layer << " " << global.z() << std::endl;

    v_eta    -> push_back ( global.eta() );
    v_phi    -> push_back ( global.phi() );
    v_layer  -> push_back ( layer        );
    v_adc    -> push_back ( adc          );
    v_charge -> push_back ( charge       );
    v_posx   -> push_back ( global.x()   );
    v_posy   -> push_back ( global.y()   );
    v_posz   -> push_back ( global.z()   );
    v_index  -> push_back ( index        );

  }
    
};

#endif
