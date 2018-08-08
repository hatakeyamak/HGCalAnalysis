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
    sprintf(tmp,"HGCDigiWaferU:HGCDigiWaferV>>h%d_cee_waferuv",i);
    sprintf(tmp2,"HGCDigiIndex==0&&HGCDigiLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_digis_CEE_waferuv.pdf(","pdf");
    else if (i==28) c1->Print("plot_digis_CEE_waferuv.pdf)","pdf");
    else            c1->Print("plot_digis_CEE_waferuv.pdf","pdf");
  }

  TH2F *h_cee_xy[28];
  for (int i=1; i<=28; i++){
    sprintf(tmp,"h%d_cee",i);
    h_cee[i] = new TH2F(tmp,tmp,61,-30.5,30.5,61,-30.5,30.5);
    sprintf(tmp,"2.*HGCDigiWaferV:-2*HGCDigiWaferU+HGCDigiWaferV>>h%d_cee",i);
    sprintf(tmp2,"HGCDigiIndex==0&&HGCDigiLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_digis_CEE_waferxy.pdf(","pdf");
    else if (i==28) c1->Print("plot_digis_CEE_waferxy.pdf)","pdf");
    else            c1->Print("plot_digis_CEE_waferxy.pdf","pdf");

  }

  TH2F *h_ceh[24];
  for (int i=1; i<=24; i++){
    sprintf(tmp,"h%d_cee_waferuv",i);
    h_ceh[i] = new TH2F(tmp,tmp,31,-15.5,15.5,31,-15.5,15.5);
    sprintf(tmp,"HGCDigiWaferU:HGCDigiWaferV>>h%d_cee_waferuv",i);
    sprintf(tmp2,"HGCDigiIndex==1&&HGCDigiLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_digis_CEH_waferuv.pdf(","pdf");
    else if (i==24) c1->Print("plot_digis_CEH_waferuv.pdf)","pdf");
    else            c1->Print("plot_digis_CEH_waferuv.pdf","pdf");
  }

  TH2F *h_ceh_xy[24];
  for (int i=1; i<=24; i++){
    sprintf(tmp,"h%d_ceh",i);
    h_ceh_xy[i] = new TH2F(tmp,tmp,61,-30.5,30.5,61,-30.5,30.5);
    sprintf(tmp,"2.*HGCDigiWaferV:-2*HGCDigiWaferU+HGCDigiWaferV>>h%d_ceh",i);
    sprintf(tmp2,"HGCDigiIndex==1&&HGCDigiLayer==%d",i);
    std::cout << tmp  << std::endl;
    std::cout << tmp2 << std::endl;
    tree->Draw(tmp,tmp2,"colz");
    if (i==1)       c1->Print("plot_digis_CEH_waferxy.pdf(","pdf");
    else if (i==24) c1->Print("plot_digis_CEH_waferxy.pdf)","pdf");
    else            c1->Print("plot_digis_CEH_waferxy.pdf","pdf");

  }

}
