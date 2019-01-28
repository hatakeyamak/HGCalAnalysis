{

  gStyle->SetTitleSize(0.08,"x");
  gStyle->SetTitleSize(0.08,"y");
  gStyle->SetOptStat(0);
  
  double Nevt=50.*200.;  // num of events x PU
  double space=2.*3.14159*0.25;  // eta binning 0.25
   
  TCanvas *c1 = new TCanvas("c1","c1",900,600);
  
  TFile *_file0 = TFile::Open("PF_D35_histos.root");
  THStack *hs = new THStack("hs","");  
  
  PFTask_PFEta_pt_ChargedHadron->Scale(1./Nevt/space);
  PFTask_PFEta_pt_ChargedHadron_NoAssociation->Scale(1./Nevt/space);
  PFTask_PFEta_pt_NeutralHadron->Scale(1./Nevt/space);
  PFTask_PFEta_pt_Photon->Scale(1./Nevt/space);
  PFTask_PFEta_pt_Electron->Scale(1./Nevt/space);
  PFTask_PFEta_pt_Muon->Scale(1./Nevt/space);
  PFTask_PFEta_pt_HFPhoton->Scale(1./Nevt/space);
  PFTask_PFEta_pt_HFHadron->Scale(1./Nevt/space);
  
  PFTask_PFEta_pt_ChargedHadron->SetFillColor(2);
  PFTask_PFEta_pt_ChargedHadron_NoAssociation->SetFillColor(11);
  PFTask_PFEta_pt_NeutralHadron->SetFillColor(6);
  PFTask_PFEta_pt_Photon->SetFillColor(3);
  PFTask_PFEta_pt_Electron->SetFillColor(4);
  PFTask_PFEta_pt_Muon->SetFillColor(12);
  PFTask_PFEta_pt_HFPhoton->SetFillColor(8);
  PFTask_PFEta_pt_HFHadron->SetFillColor(7);
  
  PFTask_PFEta_pt_ChargedHadron->SetMaximum(3.0);
  PFTask_PFEta_pt_ChargedHadron->SetTitle("");
  //PFTask_PFEta_pt_ChargedHadron->SetLabelSize(0.04,"xyz");
  PFTask_PFEta_pt_ChargedHadron->SetTitleSize(0.06,"xyz");
  PFTask_PFEta_pt_ChargedHadron->SetTitleOffset(0.7,"xyz");
  PFTask_PFEta_pt_ChargedHadron->GetXaxis()->SetTitle("#eta");
  PFTask_PFEta_pt_ChargedHadron->GetYaxis()->SetTitle("#Sigma p_{T}(GeV) / #Delta#eta#Delta#phi");
  PFTask_PFEta_pt_ChargedHadron->Draw();

  PFTask_PFEta_pt_ChargedHadron->Draw();
  hs->Add(PFTask_PFEta_pt_ChargedHadron_NoAssociation);
  hs->Add(PFTask_PFEta_pt_ChargedHadron);
  hs->Add(PFTask_PFEta_pt_Photon);
  hs->Add(PFTask_PFEta_pt_Electron);
  hs->Add(PFTask_PFEta_pt_Muon);
  hs->Add(PFTask_PFEta_pt_NeutralHadron);
  hs->Add(PFTask_PFEta_pt_HFPhoton);
  hs->Add(PFTask_PFEta_pt_HFHadron);
  hs->Draw("hist,same");

  auto legend = new TLegend(0.4,0.4,0.65,0.89);
  legend->SetLineColor(0);
  legend->AddEntry(PFTask_PFEta_pt_ChargedHadron,"Charged hadron (Assoc.)","f");
  legend->AddEntry(PFTask_PFEta_pt_ChargedHadron_NoAssociation,"Charged hadron (Unassoc.)","f");
  legend->AddEntry(PFTask_PFEta_pt_NeutralHadron,"Neutral hadron","f");
  legend->AddEntry(PFTask_PFEta_pt_Photon,"Photon","f");
  legend->AddEntry(PFTask_PFEta_pt_Electron,"Electron","f");
  legend->AddEntry(PFTask_PFEta_pt_Muon,"Muon","f");
  legend->AddEntry(PFTask_PFEta_pt_HFPhoton,"HF photon","f");
  legend->AddEntry(PFTask_PFEta_pt_HFHadron,"HF hadron","f");
  legend->Draw();

  gPad->RedrawAxis();
  
  c1->SaveAs("plot_pfeta.pdf");
  
}
