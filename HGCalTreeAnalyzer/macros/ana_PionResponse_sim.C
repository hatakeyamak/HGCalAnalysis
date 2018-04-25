// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Ken H.
//  Written on Feb 28, 2018
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
//   root> .L ana_PionResponse_sim.C++ 
//   root> ana_PionResponse_sim("../../HGCalTreeMaker/test/results_pt25_v02.root","hgcal_histograms_pt25_sim.root",-1)
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
void HGCALResponseCheckRun(TString rootfile="../HGCalNtuples_239995_Viktor.root", TString outfile="hgcal_histograms.root", int maxevents=-1, int option=2) 
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
    TTreeReaderArray<float> HGCRecHitEnergy = {fReader, "HGCRecHitEnergy"};
    TTreeReaderArray<float> HGCRecHitEta = {fReader, "HGCRecHitEta"};
    TTreeReaderArray<float> HGCRecHitPhi = {fReader, "HGCRecHitPhi"};
    TTreeReaderArray<float> HGCRecHitPosx = {fReader, "HGCRecHitPosx"};
    TTreeReaderArray<float> HGCRecHitPosy = {fReader, "HGCRecHitPosy"};
    TTreeReaderArray<float> HGCRecHitPosz = {fReader, "HGCRecHitPosz"};
    // TTreeReaderArray<float> HGCSimHitsEnergy = {fReader, "HGCSimHitsEnergy"};
    // TTreeReaderArray<float> HGCSimHitsEta = {fReader, "HGCSimHitsEta"};
    // TTreeReaderArray<float> HGCSimHitsPhi = {fReader, "HGCSimHitsPhi"};
    // TTreeReaderArray<float> HGCSimHitsPosx = {fReader, "HGCSimHitsPosx"};
    // TTreeReaderArray<float> HGCSimHitsPosy = {fReader, "HGCSimHitsPosy"};
    // TTreeReaderArray<float> HGCSimHitsPosz = {fReader, "HGCSimHitsPosz"};
    // TTreeReaderArray<float> HGCSimHitsTime = {fReader, "HGCSimHitsTime"};
    TTreeReaderArray<float> SimTracksEta = {fReader, "SimTracksEta"};
    TTreeReaderArray<float> SimTracksPhi = {fReader, "SimTracksPhi"};
    TTreeReaderArray<float> SimTracksPt = {fReader, "SimTracksPt"};
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
    TTreeReaderArray<int> HGCRecHitIndex = {fReader, "HGCRecHitIndex"};
    TTreeReaderArray<int> HGCRecHitLayer = {fReader, "HGCRecHitLayer"};
    // TTreeReaderArray<int> HGCSimHitsIndex = {fReader, "HGCSimHitsIndex"};
    // TTreeReaderArray<int> HGCSimHitsLayer = {fReader, "HGCSimHitsLayer"};
    // TTreeReaderArray<int> HGCSimHitsSubdet = {fReader, "HGCSimHitsSubdet"};
    TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
    TTreeReaderValue<UInt_t> event = {fReader, "event"};
    TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
    TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
    TTreeReaderValue<UInt_t> run = {fReader, "run"};
    
    // 
    TH1F *h_RecHitEtGenPt = new TH1F("h_RecHitEtGenPt","h_RecHitEtGenPt",100,0.,2.);

    // only when # of simtracks = 2
    TH1F *h_RecHitEtGenPt_simtrk = new TH1F("h_RecHitEtGenPt_simtrk","h_RecHitEtGenPt_simtrk",100,0.,2.);

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
      if(ievent%100==0) cout << "[HGCAL Response analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
      if (maxevents>0 && ievent>maxevents) break;

      //std::cout << *event << std::endl;
	
      // Loop over pions
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
	if (SimTracksPt.GetSize()==2) h_RecHitEtGenPt_simtrk->Fill( SumEt/TLVPion.Pt());
      } // Loop over pions ends	

    }   // Event loop ends


    // output file for histograms
    TFile file_out(outfile,"RECREATE");

    h_RecHitEtGenPt->Fit("gaus");
    h_RecHitEtGenPt->Write();
    h_RecHitEtGenPt_simtrk->Fit("gaus");
    h_RecHitEtGenPt_simtrk->Write();

    file_out.ls();
    file_out.Close();

}

//
// Main function
//
void ana_PionResponse_sim(TString rootfile="PFGtuples_local/*root",TString outfile="hgcal_histograms.root",int maxevents=-1)
{
  HGCALResponseCheckRun(rootfile, outfile, maxevents, 0);
}
