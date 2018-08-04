{

  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  TFile *_file0 = TFile::Open("ntuples_digi_pt25_numEvent1000.root");
  TTree *tree = (TTree*)_file0->Get("hgcalTupleTree/tree");

  TH2F *h9 = new TH2F("h9","h9",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h9","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==9","colz");
  c1->SaveAs("h9_ietaiphi.pdf");

  TH2F *h10 = new TH2F("h10","h10",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h10","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==10","colz");
  c1->SaveAs("h10_ietaiphi.pdf");

  TH2F *h11 = new TH2F("h11","h11",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h11","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==11","colz");
  c1->SaveAs("h11_ietaiphi.pdf");

  TH2F *h12 = new TH2F("h12","h12",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h12","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==12","colz");
  c1->SaveAs("h12_ietaiphi.pdf");

  TH2F *h13 = new TH2F("h13","h13",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h13","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==13","colz");
  c1->SaveAs("h13_ietaiphi.pdf");

  TH2F *h14 = new TH2F("h14","h14",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h14","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==14","colz");
  c1->SaveAs("h14_ietaiphi.pdf");

  TH2F *h15 = new TH2F("h15","h15",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h15","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==15","colz");
  c1->SaveAs("h15_ietaiphi.pdf");

  TH2F *h16 = new TH2F("h16","h16",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h16","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==16","colz");
  c1->SaveAs("h16_ietaiphi.pdf");

  TH2F *h17 = new TH2F("h17","h17",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h17","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==17","colz");
  c1->SaveAs("h17_ietaiphi.pdf");

  TH2F *h18 = new TH2F("h18","h18",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h18","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==18","colz");
  c1->SaveAs("h18_ietaiphi.pdf");

  TH2F *h19 = new TH2F("h19","h19",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h19","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==19","colz");
  c1->SaveAs("h19_ietaiphi.pdf");

  TH2F *h20 = new TH2F("h20","h20",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h20","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==20","colz");
  c1->SaveAs("h20_ietaiphi.pdf");

  TH2F *h21 = new TH2F("h21","h21",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h21","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==21","colz");
  c1->SaveAs("h21_ietaiphi.pdf");

  TH2F *h22 = new TH2F("h22","h22",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h22","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==22","colz");
  c1->SaveAs("h22_ietaiphi.pdf");

  TH2F *h23 = new TH2F("h23","h23",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h23","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==23","colz");
  c1->SaveAs("h23_ietaiphi.pdf");

  TH2F *h24 = new TH2F("h24","h24",366,-0.5,365.5,50,0.5,50.5);
  tree->Draw("abs(HGCSimHitsIEta):abs(HGCSimHitsIPhi)>>h24","HGCSimHitsIPhi>-50.&&HGCSimHitsLayer==24","colz");
  c1->SaveAs("h24_ietaiphi.pdf");

}
