// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Ken H.
//  Written on June 21, 2018
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
//   You should have "Fig" directory for plots 
//
// Usage : 
//
//   $ root -b  
//   root> .L ana_GeomValidation_v2.C++ 
//   root> ana_GeomValidation_v2("ntuples_pt5_numEvent10.root","hgcal_histograms_pt5.root",-1)
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not

//
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="Events/bin")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetStats(0);
}


const char* GetDetName(int Subdet) 
{ 
    const char* DetName;
    if(Subdet==1) DetName = "HB"; 
    if(Subdet==2) DetName = "HE"; 
    if(Subdet==3) DetName = "HO"; 
    if(Subdet==4) DetName = "HF"; 
    return DetName;
}

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max);
void book2D(TList *v_hist, std::string name, int n, double min, double max, int n2, double min2, double max2);
void bookHistograms(TList *v_hist);

//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value);
void fill2D(TList *v_hist, std::string name, double value, double value2);

//
void HGCALResponseCheckRun(TString rootfile="../HgcalNtuples_239995_Viktor.root", TString outfile="hgcal_histograms.root", int maxevents=-1, int skipevents=0, int option=2) 
{ 

    cout << "[Pion response analyzer] Running option " << option << " for " << endl; 

    // fit pannel display option
    gStyle->SetOptFit(1011);

    //
    // Get the tree from the PFG ntuple 
    //
    TChain *ch = new TChain("hgcalTupleTree/tree");
    ch->Add(rootfile);
    printf("%d;\n",ch->GetNtrees());
    printf("%lld;\n",ch->GetEntries());

    TTreeReader     fReader(ch);  //!the tree reader

    //
    // Set up TTreeReader's
    // -- use MakeSelector of root
    //

    // Readers to access the data (delete the ones you do not need).
    TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
    TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
    TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
    TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
    TTreeReaderArray<double> GeneralTracksD0 = {fReader, "GeneralTracksD0"};
    TTreeReaderArray<double> GeneralTracksDZ = {fReader, "GeneralTracksDZ"};
    TTreeReaderArray<double> GeneralTracksEta = {fReader, "GeneralTracksEta"};
    TTreeReaderArray<double> GeneralTracksPhi = {fReader, "GeneralTracksPhi"};
    TTreeReaderArray<double> GeneralTracksPt = {fReader, "GeneralTracksPt"};
    TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
    TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
    TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
    TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
    TTreeReaderArray<float> HGCDigiCharge = {fReader, "HGCDigiCharge"};
    TTreeReaderArray<float> HGCDigiEta = {fReader, "HGCDigiEta"};
    TTreeReaderArray<float> HGCDigiPhi = {fReader, "HGCDigiPhi"};
    TTreeReaderArray<float> HGCDigiPosx = {fReader, "HGCDigiPosx"};
    TTreeReaderArray<float> HGCDigiPosy = {fReader, "HGCDigiPosy"};
    TTreeReaderArray<float> HGCDigiPosz = {fReader, "HGCDigiPosz"};
    TTreeReaderArray<float> HGCRecHitEnergy = {fReader, "HGCRecHitEnergy"};
    TTreeReaderArray<float> HGCRecHitEta = {fReader, "HGCRecHitEta"};
    TTreeReaderArray<float> HGCRecHitPhi = {fReader, "HGCRecHitPhi"};
    TTreeReaderArray<float> HGCRecHitPosx = {fReader, "HGCRecHitPosx"};
    TTreeReaderArray<float> HGCRecHitPosy = {fReader, "HGCRecHitPosy"};
    TTreeReaderArray<float> HGCRecHitPosz = {fReader, "HGCRecHitPosz"};
    TTreeReaderArray<float> HGCSimHitsEnergy = {fReader, "HGCSimHitsEnergy"};
    TTreeReaderArray<float> HGCSimHitsEta = {fReader, "HGCSimHitsEta"};
    TTreeReaderArray<float> HGCSimHitsPhi = {fReader, "HGCSimHitsPhi"};
    TTreeReaderArray<float> HGCSimHitsPosx = {fReader, "HGCSimHitsPosx"};
    TTreeReaderArray<float> HGCSimHitsPosy = {fReader, "HGCSimHitsPosy"};
    TTreeReaderArray<float> HGCSimHitsPosz = {fReader, "HGCSimHitsPosz"};
    TTreeReaderArray<float> HGCSimHitsTime = {fReader, "HGCSimHitsTime"};
    TTreeReaderArray<float> HGCUncalibratedRecHitAmplitude = {fReader, "HGCUncalibratedRecHitAmplitude"};
    TTreeReaderArray<float> HGCUncalibratedRecHitEta = {fReader, "HGCUncalibratedRecHitEta"};
    TTreeReaderArray<float> HGCUncalibratedRecHitPhi = {fReader, "HGCUncalibratedRecHitPhi"};
    TTreeReaderArray<float> HGCUncalibratedRecHitPosx = {fReader, "HGCUncalibratedRecHitPosx"};
    TTreeReaderArray<float> HGCUncalibratedRecHitPosy = {fReader, "HGCUncalibratedRecHitPosy"};
    TTreeReaderArray<float> HGCUncalibratedRecHitPosz = {fReader, "HGCUncalibratedRecHitPosz"};
    TTreeReaderArray<float> SimTracksEta = {fReader, "SimTracksEta"};
    TTreeReaderArray<float> SimTracksPhi = {fReader, "SimTracksPhi"};
    TTreeReaderArray<float> SimTracksPt = {fReader, "SimTracksPt"};
    TTreeReaderArray<float> SimTracksR = {fReader, "SimTracksR"};
    TTreeReaderArray<float> SimTracksZ = {fReader, "SimTracksZ"};
    TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
    TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
    TTreeReaderArray<int> GeneralTracksNValidHits = {fReader, "GeneralTracksNValidHits"};
    TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
    TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
    TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
    TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
    TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
    TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
    TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
    TTreeReaderArray<int> HGCDigiIndex = {fReader, "HGCDigiIndex"};
    TTreeReaderArray<int> HGCDigiLayer = {fReader, "HGCDigiLayer"};
    TTreeReaderArray<int> HGCRecHitIndex = {fReader, "HGCRecHitIndex"};
    TTreeReaderArray<int> HGCRecHitLayer = {fReader, "HGCRecHitLayer"};
    TTreeReaderArray<int> HGCSimHitsIndex = {fReader, "HGCSimHitsIndex"};
    TTreeReaderArray<int> HGCSimHitsLayer = {fReader, "HGCSimHitsLayer"};
    TTreeReaderArray<int> HGCSimHitsSubdet = {fReader, "HGCSimHitsSubdet"};
    TTreeReaderArray<int> HGCUncalibratedRecHitIndex = {fReader, "HGCUncalibratedRecHitIndex"};
    TTreeReaderArray<int> HGCUncalibratedRecHitLayer = {fReader, "HGCUncalibratedRecHitLayer"};
    TTreeReaderArray<int> SimTracksCharge = {fReader, "SimTracksCharge"};
    TTreeReaderArray<int> SimTracksPID = {fReader, "SimTracksPID"};
    TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
    TTreeReaderValue<UInt_t> event = {fReader, "event"};
    TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
    TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
    TTreeReaderValue<UInt_t> run = {fReader, "run"};
    TTreeReaderArray<unsigned short> HGCDigiADC = {fReader, "HGCDigiADC"};
    
    //
    TList *v_hist = new TList();
    bookHistograms(v_hist); // most of histograms booked here

    TH1F *h_RecHitEtGenPt = new TH1F("h_RecHitEtGenPt","h_RecHitEtGenPt",100,0.,2.);
    TH1F *h_RecHitEtGenPt_Eta1p6_2p0 = new TH1F("h_RecHitEtGenPt_Eta1p6_2p0","h_RecHitEtGenPt_Eta1p6_2p0",100,0.,2.);
    TH1F *h_RecHitEtGenPt_Eta2p0_2p4 = new TH1F("h_RecHitEtGenPt_Eta2p0_2p4","h_RecHitEtGenPt_Eta2p0_2p4",100,0.,2.);
    TH1F *h_RecHitEtGenPt_Eta2p4_2p8 = new TH1F("h_RecHitEtGenPt_Eta2p4_2p8","h_RecHitEtGenPt_Eta2p4_2p8",100,0.,2.);
    
    std::string strtmp;
   
    //
    // Loop over entries
    //
    unsigned int nentries = (Int_t)ch->GetEntries();
    //nentries = 20000; // FIXME 
    cout << "[Pion Response analyzer] The number of entries is: " << nentries << endl;

    // main event loop
    //for(unsigned int ievent = 0; ievent<nentries; ievent++) 
    //{
    //ch->GetEntry(ievent); 
    int ievent=0;
    while (fReader.Next()) {
  
      // Progress indicator 
      ievent++;
      if(ievent%10==0) cout << "[HGCAL Response analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
      if (maxevents>0 && ievent>maxevents) break;
      if (ievent<=skipevents) continue;
      
      //std::cout << "Simhit size: " << HGCSimHitsEnergy.GetSize() << std::endl;
      //std::cout << "Digi size:   " << HGCDigiADC.GetSize() << std::endl;
      //std::cout << "Rechit size: " << HGCRecHitEnergy.GetSize() << std::endl;

      //
      // --- Looping over digis
      //
      for (int irc = 0, nrc =  HGCDigiADC.GetSize(); irc < nrc; ++irc) {
	int ilayer = HGCDigiLayer[irc];
	double z = HGCDigiPosz[irc];
	double x = HGCDigiPosx[irc];
	double y = HGCDigiPosy[irc];
	double r = pow(pow(HGCDigiPosx[irc],2)+pow(HGCDigiPosy[irc],2),0.5);
	//printf("Digi: %8d %8.2f %8.2f %8.2f %8.2f\n",ilayer,z,r,x,y);

	strtmp="Digis_layer_Z";
	fill2D(v_hist, strtmp, ilayer, fabs(z));
	strtmp="Digis_layer_R";
	fill2D(v_hist, strtmp, ilayer, fabs(r));
	strtmp="Digis_Z_R";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));

	if      (HGCDigiIndex[irc]==0) strtmp="Digis_Z_R_CEE";
	else if (HGCDigiIndex[irc]==1) strtmp="Digis_Z_R_CEHF";
	else if (HGCDigiIndex[irc]==2) strtmp="Digis_Z_R_CEHB";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));

	if      (HGCDigiIndex[irc]==0) strtmp="CEE";
	else if (HGCDigiIndex[irc]==1) strtmp="CEHF";
	else if (HGCDigiIndex[irc]==2) strtmp="CEHB";	
	fill2D(v_hist, "Digis_XY_"+strtmp, x, y);
	fill2D(v_hist, "Digis_XY_"+strtmp+"_layer"+std::to_string(ilayer), x, y);
	
      }

      //
      // --- Looping over simhits
      //
      for (int irc = 0, nrc =  HGCSimHitsEnergy.GetSize(); irc < nrc; ++irc) {
	int ilayer = HGCSimHitsLayer[irc];
	double z = HGCSimHitsPosz[irc];
	double x = HGCSimHitsPosx[irc];
	double y = HGCSimHitsPosy[irc];
	double r = pow(pow(HGCSimHitsPosx[irc],2)+pow(HGCSimHitsPosy[irc],2),0.5);
	//printf("Simhits: %8d %8.2f %8.2f %8.2f %8.2f\n",ilayer,z,r,x,y);

	strtmp="Simhits_layer_Z";
	fill2D(v_hist, strtmp, ilayer, fabs(z));
	strtmp="Simhits_layer_R";
	fill2D(v_hist, strtmp, ilayer, fabs(r));
	strtmp="Simhits_Z_R";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));

	if      (HGCSimHitsIndex[irc]==0) strtmp="Simhits_Z_R_CEE";
	else if (HGCSimHitsIndex[irc]==1) strtmp="Simhits_Z_R_CEHF";
	else if (HGCSimHitsIndex[irc]==2) strtmp="Simhits_Z_R_CEHB";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));

	if      (HGCSimHitsIndex[irc]==0) strtmp="CEE";
	else if (HGCSimHitsIndex[irc]==1) strtmp="CEHF";
	else if (HGCSimHitsIndex[irc]==2) strtmp="CEHB";	
	fill2D(v_hist, "Simhits_XY_"+strtmp, x, y);
	fill2D(v_hist, "Simhits_XY_"+strtmp+"_layer"+std::to_string(ilayer), x, y);

      }

      //
      // --- Looping over rechits
      //
      for (int irc = 0, nrc =  HGCRecHitEnergy.GetSize(); irc < nrc; ++irc) {
	TLorentzVector TLVRecHit; 
	double RecHitPt=HGCRecHitEnergy[irc]/cosh(HGCRecHitEta[irc]); //p=E for rechits
	TLVRecHit.SetPtEtaPhiE(RecHitPt,HGCRecHitEta[irc],HGCRecHitPhi[irc],HGCRecHitEnergy[irc]);
	int ilayer = HGCRecHitLayer[irc];
	double z = HGCRecHitPosz[irc];
	double x = HGCRecHitPosx[irc];
	double y = HGCRecHitPosy[irc];
	double r = pow(pow(HGCRecHitPosx[irc],2)+pow(HGCRecHitPosy[irc],2),0.5);
	//printf("%8d %8.2f %8.2f\n",ilayer,z,r);

	strtmp="Rechits_layer_Z";
	fill2D(v_hist, strtmp, ilayer, fabs(z));
	strtmp="Rechits_layer_R";
	fill2D(v_hist, strtmp, ilayer, fabs(r));
	strtmp="Rechits_Z_R";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));

	if      (HGCRecHitIndex[irc]==0) strtmp="Rechits_Z_R_CEE";
	else if (HGCRecHitIndex[irc]==1) strtmp="Rechits_Z_R_CEHF";
	else if (HGCRecHitIndex[irc]==2) strtmp="Rechits_Z_R_CEHB";
	fill2D(v_hist, strtmp, fabs(z), fabs(r));
	
	if      (HGCRecHitIndex[irc]==0) strtmp="CEE";
	else if (HGCRecHitIndex[irc]==1) strtmp="CEHF";
	else if (HGCRecHitIndex[irc]==2) strtmp="CEHB";	
	fill2D(v_hist, "Rechits_XY_"+strtmp, x, y);
	fill2D(v_hist, "Rechits_XY_"+strtmp+"_layer"+std::to_string(ilayer), x, y);

      }

      //
      // --- Loop over pions
      //
      for (int iGenPar = 0, nGenPar =  GenParPt.GetSize(); iGenPar < nGenPar; ++iGenPar) {
        //std::cout << GenParPdgId[iGenPar] << std::endl;
        TLorentzVector TLVPion; TLVPion.SetPtEtaPhiM(GenParPt[iGenPar],GenParEta[iGenPar],GenParPhi[iGenPar],GenParM[iGenPar]);
        //TLVPion.Print();
        //if (fabs(TLVPion.Eta())<1.8 || fabs(TLVPion.Eta())>2.4) continue;

        // Loop over HGCRecHits
        double SumEt=0.;
        for (int irc = 0, nrc =  HGCRecHitEnergy.GetSize(); irc < nrc; ++irc) {
          TLorentzVector TLVRecHit; 
          double RecHitPt=HGCRecHitEnergy[irc]/cosh(HGCRecHitEta[irc]); //p=E for rechits
          TLVRecHit.SetPtEtaPhiE(RecHitPt,HGCRecHitEta[irc],HGCRecHitPhi[irc],HGCRecHitEnergy[irc]);      
          double dR=TLVRecHit.DeltaR(TLVPion);
          if (dR<0.3) SumEt+=TLVRecHit.Pt();      
        }
        //std::cout << SumEt/TLVPion.Pt() <<std::endl;
        h_RecHitEtGenPt->Fill( SumEt/TLVPion.Pt());
	if      (fabs(TLVPion.Eta())>1.6 && fabs(TLVPion.Eta())<=2.0) h_RecHitEtGenPt_Eta1p6_2p0->Fill( SumEt/TLVPion.Pt());
	else if (fabs(TLVPion.Eta())>2.0 && fabs(TLVPion.Eta())<=2.4) h_RecHitEtGenPt_Eta2p0_2p4->Fill( SumEt/TLVPion.Pt());
	else if (fabs(TLVPion.Eta())>2.4 && fabs(TLVPion.Eta())<=2.8) h_RecHitEtGenPt_Eta2p4_2p8->Fill( SumEt/TLVPion.Pt());
      
      }
	
    }   // Event loop ends

    // output file for histograms
    TFile file_out(outfile,"RECREATE");

    v_hist->Write();

    h_RecHitEtGenPt->Write();
    h_RecHitEtGenPt_Eta1p6_2p0->Write();
    h_RecHitEtGenPt_Eta2p0_2p4->Write();
    h_RecHitEtGenPt_Eta2p4_2p8->Write();
    
    file_out.ls();
    file_out.Close();

}

//
// Main function
//
void ana_GeomValidation_v2(TString rootfile="PFGtuples_local/*root",TString outfile="hgcal_histograms.root",int maxevents=-1,int skipevents=0)
{
  HGCALResponseCheckRun(rootfile, outfile, maxevents, skipevents, 0);
}

//
// Book histograms
//
void bookHistograms(TList *v_hist)
{

  Char_t histo[100];

  // Digis
  sprintf(histo, "Digis_layer_Z");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,250.,550.);
  sprintf(histo, "Digis_layer_R");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,0.,400.);
  sprintf(histo, "Digis_Z_R");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

  sprintf(histo, "Digis_Z_R_CEE");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Digis_Z_R_CEHF");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Digis_Z_R_CEHB");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);

  sprintf(histo, "Digis_XY_CEE");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Digis_XY_CEHF");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Digis_XY_CEHB");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);

  for (int ilayer=1; ilayer<=28; ilayer++){
    sprintf(histo, "Digis_XY_CEE_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Digis_XY_CEHF_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Digis_XY_CEHB_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  }

  // Simhits  
  sprintf(histo, "Simhits_layer_Z");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,250.,550.);
  sprintf(histo, "Simhits_layer_R");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,0.,400.);
  sprintf(histo, "Simhits_Z_R");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

  sprintf(histo, "Simhits_Z_R_CEE");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Simhits_Z_R_CEHF");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Simhits_Z_R_CEHB");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);

  sprintf(histo, "Simhits_XY_CEE");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Simhits_XY_CEHF");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Simhits_XY_CEHB");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);

  for (int ilayer=1; ilayer<=28; ilayer++){
    sprintf(histo, "Simhits_XY_CEE_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Simhits_XY_CEHF_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Simhits_XY_CEHB_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  }

  // Rechits
  sprintf(histo, "Rechits_layer_Z");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,250.,550.);
  sprintf(histo, "Rechits_layer_R");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,0.,400.);
  sprintf(histo, "Rechits_Z_R");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

  sprintf(histo, "Rechits_Z_R_CEE");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_Z_R_CEHF");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_Z_R_CEHB");
  book2D(v_hist, histo, 300,250.,550.,100,0.,400.);

  sprintf(histo, "Rechits_XY_CEE");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Rechits_XY_CEHF");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  sprintf(histo, "Rechits_XY_CEHB");
  book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);

  for (int ilayer=1; ilayer<=28; ilayer++){
    sprintf(histo, "Rechits_XY_CEE_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Rechits_XY_CEHF_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
    sprintf(histo, "Rechits_XY_CEHB_layer%d",ilayer);
    book2D(v_hist, histo, 400,-400.,400.,400,-400.,400.);
  }
  
}

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max)
{
  TH1D *htemp = new TH1D(name.c_str(), name.c_str(), n, min, max);
  v_hist->Add(htemp);
}
//
// Book 2D histograms
//
void book2D(TList *v_hist, std::string name, int n, double min, double max, int n2, double min2, double max2)
{
  TH2D *htemp = new TH2D(name.c_str(), name.c_str(), n, min, max, n2, min2, max2);
  v_hist->Add(htemp);
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1D* htemp = (TH1D*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}
//
// Fill 2D histograms
//
void fill2D(TList *v_hist, std::string name, double value, double value2)
{
  TH2D* htemp = (TH2D*) v_hist->FindObject(name.c_str());
  htemp->Fill(value,value2);
}
