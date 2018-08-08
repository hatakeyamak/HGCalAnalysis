{

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gPad->SetGridx();
  gPad->SetGridy();
  
  TFile *_file0 = TFile::Open("ntuples_digi_pt25_numEvent10000.root");
  TTree *tree = (TTree*)_file0->Get("hgcalTupleTree/tree");

  char tmp[100];
  char tmp2[100];

  TH2F *h_cee[28];
  for (int i=1; i<=28; i++){
    sprintf(tmp,"h%d_cee_waferuv",i);
    h_ceh[i] = new TH2F(tmp,tmp,31,-15.5,15.5,31,-15.5,15.5);
    sprintf(tmp,"HGCSimHitsWaferU:HGCSimHitsWaferV>>h%d_cee_waferuv",i);
    sprintf(tmp2,"HGCSimHitsIndex==0&&HGCSimHitsLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_simhits_CEE_waferuv.pdf(","pdf");
    else if (i==28) c1->Print("plot_simhits_CEE_waferuv.pdf)","pdf");
    else            c1->Print("plot_simhits_CEE_waferuv.pdf","pdf");
  }

  TH2F *h_cee_xy[28];
  for (int i=1; i<=28; i++){
    sprintf(tmp,"h%d_cee",i);
    h_cee[i] = new TH2F(tmp,tmp,61,-30.5,30.5,61,-30.5,30.5);
    sprintf(tmp,"2.*HGCSimHitsWaferV:-2*HGCSimHitsWaferU+HGCSimHitsWaferV>>h%d_cee",i);
    sprintf(tmp2,"HGCSimHitsIndex==0&&HGCSimHitsLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_simhits_CEE_waferxy.pdf(","pdf");
    else if (i==28) c1->Print("plot_simhits_CEE_waferxy.pdf)","pdf");
    else            c1->Print("plot_simhits_CEE_waferxy.pdf","pdf");

  }

  TH2F *h_ceh[24];
  for (int i=1; i<=24; i++){
    sprintf(tmp,"h%d_cee_waferuv",i);
    h_ceh[i] = new TH2F(tmp,tmp,31,-15.5,15.5,31,-15.5,15.5);
    sprintf(tmp,"HGCSimHitsWaferU:HGCSimHitsWaferV>>h%d_cee_waferuv",i);
    sprintf(tmp2,"HGCSimHitsIndex==1&&HGCSimHitsLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_simhits_CEH_waferuv.pdf(","pdf");
    else if (i==24) c1->Print("plot_simhits_CEH_waferuv.pdf)","pdf");
    else            c1->Print("plot_simhits_CEH_waferuv.pdf","pdf");
  }

  TH2F *h_ceh_xy[24];
  for (int i=1; i<=24; i++){
    sprintf(tmp,"h%d_ceh",i);
    h_ceh_xy[i] = new TH2F(tmp,tmp,61,-30.5,30.5,61,-30.5,30.5);
    sprintf(tmp,"2.*HGCSimHitsWaferV:-2*HGCSimHitsWaferU+HGCSimHitsWaferV>>h%d_ceh",i);
    sprintf(tmp2,"HGCSimHitsIndex==1&&HGCSimHitsLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_simhits_CEH_waferxy.pdf(","pdf");
    else if (i==24) c1->Print("plot_simhits_CEH_waferxy.pdf)","pdf");
    else            c1->Print("plot_simhits_CEH_waferxy.pdf","pdf");

  }

  /*
  TH2F *h9 = new TH2F("h9","h9",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h9","HGCSimHitsIndex==1&&HGCSimHitsLayer==9","colz");
  c1->SaveAs("h9_CEH_celluv.pdf");

  TH2F *h10 = new TH2F("h10","h10",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h10","HGCSimHitsIndex==1&&HGCSimHitsLayer==10","colz");
  c1->SaveAs("h10_CEH_celluv.pdf");

  TH2F *h11 = new TH2F("h11","h11",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h11","HGCSimHitsIndex==1&&HGCSimHitsLayer==11","colz");
  c1->SaveAs("h11_CEH_celluv.pdf");

  TH2F *h12 = new TH2F("h12","h12",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h12","HGCSimHitsIndex==1&&HGCSimHitsLayer==12","colz");
  c1->SaveAs("h12_CEH_celluv.pdf");

  TH2F *h13 = new TH2F("h13","h13",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h13","HGCSimHitsIndex==1&&HGCSimHitsLayer==13","colz");
  c1->SaveAs("h13_CEH_celluv.pdf");

  TH2F *h14 = new TH2F("h14","h14",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h14","HGCSimHitsIndex==1&&HGCSimHitsLayer==14","colz");
  c1->SaveAs("h14_CEH_celluv.pdf");

  TH2F *h15 = new TH2F("h15","h15",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h15","HGCSimHitsIndex==1&&HGCSimHitsLayer==15","colz");
  c1->SaveAs("h15_CEH_celluv.pdf");

  TH2F *h16 = new TH2F("h16","h16",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h16","HGCSimHitsIndex==1&&HGCSimHitsLayer==16","colz");
  c1->SaveAs("h16_CEH_celluv.pdf");

  TH2F *h17 = new TH2F("h17","h17",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h17","HGCSimHitsIndex==1&&HGCSimHitsLayer==17","colz");
  c1->SaveAs("h17_CEH_celluv.pdf");

  TH2F *h18 = new TH2F("h18","h18",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h18","HGCSimHitsIndex==1&&HGCSimHitsLayer==18","colz");
  c1->SaveAs("h18_CEH_celluv.pdf");

  TH2F *h19 = new TH2F("h19","h19",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h19","HGCSimHitsIndex==1&&HGCSimHitsLayer==19","colz");
  c1->SaveAs("h19_CEH_celluv.pdf");

  TH2F *h20 = new TH2F("h20","h20",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h20","HGCSimHitsIndex==1&&HGCSimHitsLayer==20","colz");
  c1->SaveAs("h20_CEH_celluv.pdf");

  TH2F *h21 = new TH2F("h21","h21",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h21","HGCSimHitsIndex==1&&HGCSimHitsLayer==21","colz");
  c1->SaveAs("h21_CEH_celluv.pdf");

  TH2F *h22 = new TH2F("h22","h22",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h22","HGCSimHitsIndex==1&&HGCSimHitsLayer==22","colz");
  c1->SaveAs("h22_CEH_celluv.pdf");

  TH2F *h23 = new TH2F("h23","h23",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h23","HGCSimHitsIndex==1&&HGCSimHitsLayer==23","colz");
  c1->SaveAs("h23_CEH_celluv.pdf");

  TH2F *h24 = new TH2F("h24","h24",31,-15.5,15.5,31,-15.5,15.5);
  tree->Draw("HGCSimHitsWaferU:HGCSimHitsWaferV>>h24","HGCSimHitsIndex==1&&HGCSimHitsLayer==24","colz");
  c1->SaveAs("h24_CEH_celluv.pdf");
  */

}
