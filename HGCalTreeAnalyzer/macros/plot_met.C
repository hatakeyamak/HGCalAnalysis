{

  std::vector<std::string> files;
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_1.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_2.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_3.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_14.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_5.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_6.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_7.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_8.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_19.root");
  files.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_Step3_v3/180718_215359/0000/step3_inMINIAODSIM_10.root");

  std::vector<std::string> files_D30;
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_10.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_100.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_102.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_103.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_104.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_108.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_109.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_11.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_112.root");
  files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_114.root");
  //files_D30.push_back("/cms/data/store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_2_0_D30_Step3_v3/180718_231532/0000/step3_inMINIAODSIM_119.root");
  
  TChain chain_D28("Events");
  for(vector<string>::const_iterator i = files.begin(); i != files.end(); ++i) {
    //cout << *i << " "; // this will print all the contents of *features*
    TString fname = *i;
    chain_D28.Add(fname);
  }

  TChain chain_D30("Events");
  for(vector<string>::const_iterator i = files_D30.begin(); i != files_D30.end(); ++i) {
    //cout << *i << " "; // this will print all the contents of *features*
    TString fname = *i;
    chain_D30.Add(fname);
  }

  //
  // c1
  //
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  gPad->SetLogy();
  
  chain_D28.Draw("patMETs_slimmedMETs__RECO.obj.pt()");
  TH1F *pfmet2 = (TH1F*) htemp->Clone("pfmet2");  
  htemp->SetLineColor(8);
  pfmet2->SetLineColor(8);
  pfmet2->Draw();
  pfmet2->Print();

  chain_D30.Draw("patMETs_slimmedMETs__RECO.obj.pt()","","sames");
  TH1F *pfmet0 = (TH1F*) htemp->Clone("pfmet0");
  pfmet0->Draw("same");
  pfmet0->Print();

  c1->Print("pfmet.pdf");

  //
  // c2
  //
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  gPad->SetLogy();
  
  chain_D28.Draw("recoGenMETs_genMetTrue__HLT.obj.pt()");
  TH1F *genmet2 = (TH1F*) htemp->Clone("genmet2");  
  htemp->SetLineColor(8);
  genmet2->SetLineColor(8);
  genmet2->Draw();
  genmet2->Print();

  chain_D30.Draw("recoGenMETs_genMetTrue__HLT.obj.pt()","","sames");
  TH1F *genmet0 = (TH1F*) htemp->Clone("genmet0");
  genmet0->Draw("same");
  genmet0->Print();
  
  c2->Print("genmet.pdf");

  //
  // c3
  //
  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  gPad->SetLogy();
  
  chain_D28.Draw("patMETs_slimmedMETs__RECO.obj.pt()","recoGenMETs_genMetTrue__HLT.obj.pt()<10.");
  TH1F *pfmet_lowgenmet2 = (TH1F*) htemp->Clone("pfmet_lowgenmet2 age");  
  htemp->SetLineColor(8);
  pfmet_lowgenmet2->SetLineColor(8);
  pfmet_lowgenmet2->Draw();
  pfmet_lowgenmet2->Print();

  chain_D30.Draw("patMETs_slimmedMETs__RECO.obj.pt()","recoGenMETs_genMetTrue__HLT.obj.pt()<10.","sames");
  TH1F *pfmet_lowgenmet0 = (TH1F*) htemp->Clone("pfmet_lowgenmet0 noage");
  pfmet_lowgenmet0->Draw("same");
  pfmet_lowgenmet0->Print();

  c3->Print("pfmet_lowgenmet.pdf");

  //
  // c4 - pfjets
  //
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->cd();
  gPad->SetLogy();

  chain_D28.Draw("patJets_slimmedJets__RECO.obj.eta()","");
  TH1F *pfjeteta2 = (TH1F*) htemp->Clone("pfjeteta2 age");  
  htemp->SetLineColor(8);
  pfjeteta2->SetLineColor(8);
  pfjeteta2->Draw();
  pfjeteta2->Print();

  chain_D30.Draw("patJets_slimmedJets__RECO.obj.eta()","","sames");
  TH1F *pfjeteta0 = (TH1F*) htemp->Clone("pfjeteta0 noage");
  pfjeteta0->Draw("same");
  pfjeteta0->Print();

  c4->Print("pfjet_eta.pdf");

  //
  // c5 - calojets
  //
  TCanvas *c5 = new TCanvas("c5","c5");
  c5->cd();
  gPad->SetLogy();

  chain_D28.Draw("recoCaloJets_slimmedCaloJets__RECO.obj.eta()","");
  TH1F *calojeteta2 = (TH1F*) htemp->Clone("calojeteta2 age");  
  htemp->SetLineColor(8);
  calojeteta2->SetLineColor(8);
  calojeteta2->Draw();
  calojeteta2->Print();

  chain_D30.Draw("recoCaloJets_slimmedCaloJets__RECO.obj.eta()","","sames");
  TH1F *calojeteta0 = (TH1F*) htemp->Clone("calojeteta0 noage");
  calojeteta0->Draw("same");
  calojeteta0->Print();

  c5->Print("calojet_eta.pdf");
  
}
