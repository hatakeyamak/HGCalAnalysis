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
//   root> .L ana_GeomValidation.C++ 
//   root> ana_GeomValidation("ntuples_pt5_numEvent10.root","hgcal_histograms_pt5.root",-1)
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
void HGCALResponseCheckRun(TString rootfile="../HgcalNtuples_239995_Viktor.root", TString outfile="hgcal_histograms.root", int maxevents=-1, int option=2) 
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

    TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
    TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
    TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
    TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
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
    TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
    TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
    TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
    TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
    TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
    TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
    TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
    TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
    TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
    TTreeReaderArray<int> HGCRecHitLayer = {fReader, "HGCRecHitLayer"};
    TTreeReaderArray<int> HGCRecHitIndex = {fReader, "HGCRecHitIndex"};
    TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
    TTreeReaderValue<UInt_t> event = {fReader, "event"};
    TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
    TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
    TTreeReaderValue<UInt_t> run = {fReader, "run"};
    
    //
    TList *v_hist = new TList();
    bookHistograms(v_hist); // most of histograms booked here
    
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
      if(ievent%100==0) cout << "[HGCAL Response analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
      if (maxevents>0 && ievent>maxevents) break;
      
      for (int irc = 0, nrc =  HGCRecHitEnergy.GetSize(); irc < nrc; ++irc) {
	TLorentzVector TLVRecHit; 
	double RecHitPt=HGCRecHitEnergy[irc]/cosh(HGCRecHitEta[irc]); //p=E for rechits
	TLVRecHit.SetPtEtaPhiE(RecHitPt,HGCRecHitEta[irc],HGCRecHitPhi[irc],HGCRecHitEnergy[irc]);
	int ilayer = HGCRecHitLayer[irc];
	double z = HGCRecHitPosz[irc];
	double r = pow(pow(HGCRecHitPosx[irc],2)+pow(HGCRecHitPosy[irc],2),0.5);
	printf("%8d %8.2f %8.2f\n",ilayer,z,r);

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
	
      }
    
    }   // Event loop ends

    // output file for histograms
    TFile file_out(outfile,"RECREATE");

    v_hist->Write();
    
    file_out.ls();
    file_out.Close();

}

//
// Main function
//
void ana_GeomValidation(TString rootfile="PFGtuples_local/*root",TString outfile="hgcal_histograms.root",int maxevents=-1)
{
  HGCALResponseCheckRun(rootfile, outfile, maxevents, 0);
}

//
// Book histograms
//
void bookHistograms(TList *v_hist)
{

  Char_t histo[100];

  sprintf(histo, "Rechits_layer_Z");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,250.,550.);
  sprintf(histo, "Rechits_layer_R");
  book2D(v_hist, histo, 55, -0.5, 54.5,100,0.,400.);
  sprintf(histo, "Rechits_Z_R");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

  sprintf(histo, "Rechits_Z_R_CEE");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_Z_R_CEHF");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_Z_R_CEHB");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

  sprintf(histo, "Rechits_XY_CEE");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_XY_CEHF");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);
  sprintf(histo, "Rechits_Z_R_CEHB");
  book2D(v_hist, histo, 100,250.,550.,100,0.,400.);

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
