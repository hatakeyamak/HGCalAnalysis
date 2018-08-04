{

  TFile *_file0 = TFile::Open("hgcal_histograms_pt25_v02.root");

  gStyle->SetOptStat(0);

  TH2F *Simhits_Z_R = (TH2F*)_file0->Get("Simhits_Z_R");
  TH2F *Simhits_Z_R_CEE = (TH2F*)_file0->Get("Simhits_Z_R_CEE");
  TH2F *Simhits_Z_R_CEHF = (TH2F*)_file0->Get("Simhits_Z_R_CEHF");
  TH2F *Simhits_Z_R_CEHB = (TH2F*)_file0->Get("Simhits_Z_R_CEHB");

  TH2F *Digis_Z_R = (TH2F*)_file0->Get("Digis_Z_R");
  TH2F *Digis_Z_R_CEE = (TH2F*)_file0->Get("Digis_Z_R_CEE");
  TH2F *Digis_Z_R_CEHF = (TH2F*)_file0->Get("Digis_Z_R_CEHF");
  TH2F *Digis_Z_R_CEHB = (TH2F*)_file0->Get("Digis_Z_R_CEHB");

  TH2F *Rechits_Z_R = (TH2F*)_file0->Get("Rechits_Z_R");
  TH2F *Rechits_Z_R_CEE = (TH2F*)_file0->Get("Rechits_Z_R_CEE");
  TH2F *Rechits_Z_R_CEHF = (TH2F*)_file0->Get("Rechits_Z_R_CEHF");
  TH2F *Rechits_Z_R_CEHB = (TH2F*)_file0->Get("Rechits_Z_R_CEHB");

  bool detail=false;
  
  //-----
  /*
  double z_layer[52]={
  3198.0,    3207.1,    3222.4,    3231.5,    3246.8,
  3255.9,    3271.2,    3280.3,    3295.6,    3304.7,
  3320.0,    3329.1,    3344.4,    3353.5,    3368.8,
  3377.9,    3393.2,    3402.3,    3417.6,    3426.7,
  3442.0,    3451.1,    3466.4,    3475.5,    3490.8,
  3499.9,    3515.2,    3524.3,    3577.4,    3626.4,
  3675.4,    3724.4,    3773.4,    3822.4,    3871.4,
  3920.4,    3969.4,    4020.3,    4071.2,    4122.1,
  4206.0,    4289.9,    4373.8,    4457.7,    4541.6,
  4625.5,    4709.4,    4793.3,    4877.2,    4961.1,
  5045.0,    5128.9};
  */

  //
  // Philip Bloch's maps
  // https://www.dropbox.com/sh/sw82xi3ayyraddq/AACvBvknyqF4bAbaYdOeJeaEa?dl=0
  //
  double z_layer[52]={
    319.8, 320.7, 322.1, 323.0, 324.4,
    325.3, 326.8, 327.6, 329.1, 330.0,
    331.4, 332.3, 333.8, 334.6, 336.1,
    337.0, 338.4, 339.3, 340.8, 341.6,
    343.1, 344.0, 345.4, 346.3, 347.7,
    348.6, 350.1, 350.9, 356.2, 361.0,
    365.8, 370.6, 375.4, 380.1, 384.9,
    389.7, 394.5, 399.3, 404.0, 408.8,
    417.5, 426.2, 434.9, 443.5, 452.2,    
    460.9, 469.6, 478.3, 486.9, 495.6,    
    504.3, 513.0
  };

  double rout_layer[52]={
    1567.5, 1567.5, 1575.6, 1575.6, 1583.7,
    1583.7, 1591.8, 1591.8, 1599.9, 1599.9,
    1608.0, 1608.0, 1616.1, 1616.1, 1624.2,
    1624.2, 1632.2, 1632.2, 1640.4, 1640.4,
    1648.5, 1648.5, 1656.6, 1656.6, 1664.7,
    1664.7, 1672.8, 1672.8, 1696.9, 1713.5,
    1730.1, 1746.7, 1763.3, 1779.9, 1796.4,
    1844.2, 1907.6, 1971.0, 2034.5, 2097.9,
    2184.6, 2299.8, 2415.0, 2530.2, 2645.3,
    2664.0, 2664.0, 2664.0, 2664.0, 2664.0,
    2664.0, 2664.0
  };

  double rmid_layer[52]={
    1567.5, 1567.5, 1575.6, 1575.6, 1583.7,
    1583.7, 1591.8, 1591.8, 1599.9, 1599.9,
    1608.0, 1608.0, 1616.1, 1616.1, 1624.2,
    1624.2, 1632.2, 1632.2, 1640.4, 1640.4,
    1648.5, 1648.5, 1656.6, 1656.6, 1664.7,
    1664.7, 1672.8, 1672.8, 1696.9, 1713.5,
    1730.1, 1746.7, 1763.3, 1779.9, 1796.4,
    1844.2, 1462.9, 1421.1, 1300.6, 1300.6,  // <=== boundary kicks in from this line
    1267.8, 1139.7,  981.1,  981.1,  981.1,
    981.1,  981.1,  981.1,  981.1,  981.1,
    1114.9, 1114.9
  };

  double rout_layer_cm[52];
  double rin_layer_cm[52];
  double rmid_layer_cm[52];
  for (int i=0; i<=51; i++){
    rout_layer_cm[i]=rout_layer[i]/10.;
    rmid_layer_cm[i]=rmid_layer[i]/10.;
    rin_layer_cm[i]=z_layer[i]/sinh(3.);
    std::cout << i << " " << z_layer[i] << " " << rout_layer_cm[i] << " " << rin_layer_cm[i] << std::endl;
  }
  
  TPolyLine *pline_out = new TPolyLine(52,z_layer,rout_layer_cm);
  TPolyLine *pline_in = new TPolyLine(52,z_layer,rin_layer_cm);
  TPolyLine *pline_mid = new TPolyLine(52,z_layer,rmid_layer_cm);
  TLine *line = new TLine((z_layer[27]+z_layer[28])/2.,(rout_layer_cm[27]+rout_layer_cm[28])/2.,
			  (z_layer[27]+z_layer[28])/2.,(rin_layer_cm[27]+rin_layer_cm[28])/2.);
  
  TEllipse *el_out[52];
  TEllipse *el_in[52];
  TEllipse *el_mid[52];
  for (int i=0; i<=51; i++){
    el_out[i] = new TEllipse(0.0,0.0,rout_layer_cm[i]);
    el_in[i]  = new TEllipse(0.0,0.0,rin_layer_cm[i]);
    el_mid[i] = new TEllipse(0.0,0.0,rmid_layer_cm[i]);
    el_out[i]->SetLineWidth(2);
    el_out[i]->SetLineColor(2);
    el_out[i]->SetFillStyle(4000);
    el_in[i]->SetLineWidth(2);
    el_in[i]->SetLineColor(2);
    el_in[i]->SetFillStyle(4000);
    el_mid[i]->SetLineWidth(2);
    el_mid[i]->SetLineColor(2);
    el_mid[i]->SetFillStyle(4000);
  }  

  TH2F *Simhits_XY_CEE[28];
  TH2F *Simhits_XY_CEHF[28];
  TH2F *Simhits_XY_CEHB[28];
  Char_t histo[100];
  Char_t pdffile[100];
  for (int i=0; i<=27; i++){
    sprintf(histo, "Simhits_XY_CEE_layer%d",i+1);
    Simhits_XY_CEE[i] = (TH2F*)_file0->Get(histo);
    Simhits_XY_CEE[i]->GetXaxis()->SetRange(76,325);
    Simhits_XY_CEE[i]->GetYaxis()->SetRange(76,325);
    sprintf(histo, "Simhits_XY_CEHF_layer%d",i+1);
    Simhits_XY_CEHF[i] = (TH2F*)_file0->Get(histo);
    Simhits_XY_CEHF[i]->GetXaxis()->SetRange(66,335);
    Simhits_XY_CEHF[i]->GetYaxis()->SetRange(66,335);
    sprintf(histo, "Simhits_XY_CEHB_layer%d",i+1);
    Simhits_XY_CEHB[i] = (TH2F*)_file0->Get(histo);
    Simhits_XY_CEHB[i]->GetXaxis()->SetRange(66,335);
    Simhits_XY_CEHB[i]->GetYaxis()->SetRange(66,335);
  }  
  
  TH2F *Digis_XY_CEE[28];
  TH2F *Digis_XY_CEHF[28];
  TH2F *Digis_XY_CEHB[28];
  for (int i=0; i<=27; i++){
    sprintf(histo, "Digis_XY_CEE_layer%d",i+1);
    Digis_XY_CEE[i] = (TH2F*)_file0->Get(histo);
    Digis_XY_CEE[i]->GetXaxis()->SetRange(76,325);
    Digis_XY_CEE[i]->GetYaxis()->SetRange(76,325);
    sprintf(histo, "Digis_XY_CEHF_layer%d",i+1);
    Digis_XY_CEHF[i] = (TH2F*)_file0->Get(histo);
    Digis_XY_CEHF[i]->GetXaxis()->SetRange(51,350);
    Digis_XY_CEHF[i]->GetYaxis()->SetRange(51,350);
    sprintf(histo, "Digis_XY_CEHB_layer%d",i+1);
    Digis_XY_CEHB[i] = (TH2F*)_file0->Get(histo);
    Digis_XY_CEHB[i]->GetXaxis()->SetRange(51,350);
    Digis_XY_CEHB[i]->GetYaxis()->SetRange(51,350);
  }  

  TH2F *Rechits_XY_CEE[28];
  TH2F *Rechits_XY_CEHF[28];
  TH2F *Rechits_XY_CEHB[28];
  for (int i=0; i<=27; i++){
    sprintf(histo, "Rechits_XY_CEE_layer%d",i+1);
    Rechits_XY_CEE[i] = (TH2F*)_file0->Get(histo);
    Rechits_XY_CEE[i]->GetXaxis()->SetRange(76,325);
    Rechits_XY_CEE[i]->GetYaxis()->SetRange(76,325);
    sprintf(histo, "Rechits_XY_CEHF_layer%d",i+1);
    Rechits_XY_CEHF[i] = (TH2F*)_file0->Get(histo);
    Rechits_XY_CEHF[i]->GetXaxis()->SetRange(66,335);
    Rechits_XY_CEHF[i]->GetYaxis()->SetRange(66,335);
    sprintf(histo, "Rechits_XY_CEHB_layer%d",i+1);
    Rechits_XY_CEHB[i] = (TH2F*)_file0->Get(histo);
    Rechits_XY_CEHB[i]->GetXaxis()->SetRange(66,335);
    Rechits_XY_CEHB[i]->GetYaxis()->SetRange(66,335);
  }  
  
  //-----
  TCanvas *c0 = new TCanvas("c0","c0",2800,1000);
  c0->SetGridx();
  c0->SetGridy();

  c0->Divide(2,1);
  c0->cd(1);
  Simhits_Z_R->Draw("colz");
  pline_out->SetLineColor(2);
  pline_out->SetLineWidth(2);
  pline_out->Draw();
  pline_in->SetLineColor(2);
  pline_in->SetLineWidth(2);
  pline_in->Draw();
  pline_mid->SetLineColor(2);
  pline_mid->SetLineWidth(2);
  //pline_mid->Draw();
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->Draw();
  
  c0->cd(2);
  Simhits_Z_R_CEE->Draw("colz");
  
  c0->Print("plot_geom.pdf(","pdf");

  c0->cd(1);
  Simhits_Z_R_CEHF->Draw("colz");
  
  c0->cd(2);
  Simhits_Z_R_CEHB->Draw("colz");

  //if (ilayer==36)      c1->Print("plot_geom.pdf(","pdf");
  //else if (ilayer==51) c1->Print("plot_geom.pdf)","pdf");
  //else                 c1->Print("plot_geom.pdf","pdf");

  c0->Print("plot_geom.pdf)","pdf");

  //
  // CE-E
  //----
  TCanvas *c1 = new TCanvas("c1","c1",2000,1000);
  c1->SetGridx();
  c1->SetGridy();

  c1->Divide(2,1);
  
  for (int i=0; i<=27; i++){

    double location = i%2+1;
    
    c1->cd(location);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Simhits_XY_CEE[i]->Draw("colz");
    el_in[i]->Draw();
    el_out[i]->Draw();

    sprintf(pdffile, "Simhits_XY_CEE_layer%d.png",i);
    if (location==2) c1->Print(pdffile,"png");
    
  }

  //
  // CE-H
  //----
  for (int i=0; i<=7; i++){

    double location = i%2+1;
    
    c1->cd(location);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Simhits_XY_CEHF[i]->Draw("colz");
    el_in[i+28]->Draw();
    el_out[i+28]->Draw();

    sprintf(pdffile, "Simhits_XY_CEHF_layer%d.png",i);
    if (location==2 && i!=27) c1->Print(pdffile,"png");
    
  }

  //
  // CE-H mixed
  //----
  for (int i=8; i<=23; i++){

    c1->cd(1);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Simhits_XY_CEHF[i]->Draw("colz");
    el_in[i+28]->Draw();
    el_mid[i+28]->Draw();

    c1->cd(2);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Simhits_XY_CEHB[i]->Draw("colz");
    el_mid[i+28]->Draw();
    el_out[i+28]->Draw();

    sprintf(pdffile, "Simhits_XY_CEH_layer%d.png",i);
    c1->Print(pdffile,"png");
    
  }

  //
  // Digis
  //-----
  c0->cd(1);
  Digis_Z_R->Draw("colz");
  pline_out->SetLineColor(2);
  pline_out->SetLineWidth(2);
  pline_out->Draw();
  pline_in->SetLineColor(2);
  pline_in->SetLineWidth(2);
  pline_in->Draw();
  pline_mid->SetLineColor(2);
  pline_mid->SetLineWidth(2);
  //pline_mid->Draw();
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->Draw();
  
  c0->cd(2);
  Digis_Z_R_CEE->Draw("colz");
  
  c0->Print("plot_geom_digis.pdf(","pdf");

  c0->cd(1);
  Digis_Z_R_CEHF->Draw("colz");
  
  c0->cd(2);
  Digis_Z_R_CEHB->Draw("colz");

  c0->Print("plot_geom_digis.pdf)","pdf");

  //
  // CE-E
  //----
  for (int i=0; i<=27; i++){

    double location = i%2+1;
    
    c1->cd(location);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Digis_XY_CEE[i]->Draw("colz");
    el_in[i]->Draw();
    el_out[i]->Draw();

    sprintf(pdffile, "Digis_XY_CEE_layer%d.png",i);
    if (location==2) c1->Print(pdffile,"png");
    
  }

  //
  // CE-H
  //----
  for (int i=0; i<=7; i++){

    double location = i%2+1;
    
    c1->cd(location);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Digis_XY_CEHF[i]->Draw("colz");
    el_in[i+28]->Draw();
    el_out[i+28]->Draw();

    sprintf(pdffile, "Digis_XY_CEHF_layer%d.png",i);
    if (location==2 && i!=27) c1->Print(pdffile,"png");
    
  }

  //
  // CE-H mixed
  //----
  for (int i=8; i<=23; i++){

    c1->cd(1);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Digis_XY_CEHF[i]->Draw("colz");
    el_in[i+28]->Draw();
    el_mid[i+28]->Draw();

    c1->cd(2);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    Digis_XY_CEHB[i]->Draw("colz");
    el_mid[i+28]->Draw();
    el_out[i+28]->Draw();

    sprintf(pdffile, "Digis_XY_CEH_layer%d.png",i);
    c1->Print(pdffile,"png");
    
  }

  //
  // Rechits
  //-----
  c0->cd(1);
  Rechits_Z_R->Draw("colz");
  pline_out->SetLineColor(2);
  pline_out->SetLineWidth(2);
  pline_out->Draw();
  pline_in->SetLineColor(2);
  pline_in->SetLineWidth(2);
  pline_in->Draw();
  pline_mid->SetLineColor(2);
  pline_mid->SetLineWidth(2);
  //pline_mid->Draw();
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->Draw();
  
  c0->cd(2);
  Rechits_Z_R_CEE->Draw("colz");
  
  c0->Print("plot_geom_rechits.pdf(","pdf");

  c0->cd(1);
  Rechits_Z_R_CEHF->Draw("colz");
  
  c0->cd(2);
  Rechits_Z_R_CEHB->Draw("colz");

  c0->Print("plot_geom_rechits.pdf)","pdf");
  
  /*

  TFile *_file0 = TFile::Open("et100_eta1p7_v4.root");
  TFile *_file1 = TFile::Open("et100_eta1p7_v6.root");

  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->SetGridx();
  c1->SetGridy();

  bool develop=false;
  
  double z_layer[69]={
  3198.0,    3207.1,    3222.4,    3231.5,    3246.8,
  3255.9,    3271.2,    3280.3,    3295.6,    3304.7,
  3320.0,    3329.1,    3344.4,    3353.5,    3368.8,
  3377.9,    3393.2,    3402.3,    3417.6,    3426.7,
  3442.0,    3451.1,    3466.4,    3475.5,    3490.8,
  3499.9,    3515.2,    3524.3,    3577.4,    3626.4,
  3675.4,    3724.4,    3773.4,    3822.4,    3871.4,
  3920.4,    3969.4,    4020.3,    4071.2,    4122.1,
  4206.0,    4289.9,    4373.8,    4457.7,    4541.6,
  4625.5,    4709.4,    4793.3,    4877.2,    4961.1,
  5045.0,    5128.9,       0.0,    3971.2,    4022.1,
  4073.0,    4123.9,    4207.8,    4291.7,    4375.6,
  4459.5,    4543.4,    4627.3,    4711.2,    4795.1,
  4879.0,    4962.9,    5046.8,    5130.7};

  int eta_index[69]={
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,19,24, 25,30,29,32,37,
    41,41,42,43,44, 45,45,46,47};
  
  TH2Poly *map_Si[100];
  TH2F    *map_BH[100];
  char title[100];
  for (int ilayer=36;ilayer<=51;ilayer++){
    //
    // Silicon segmentation
    //
    sprintf(title,"map_1_%d",ilayer);
    map_Si[ilayer] = (TH2Poly*)_file0->Get(title);
    map_Si[ilayer]->SetMinimum(0.5);
    sprintf(title,"Layer %d",ilayer);
    map_Si[ilayer]->SetTitle(title);
    //
    // Sintilator segmentation
    //
    if (ilayer>=36 && ilayer<=39) {
      sprintf(title,"map_TH2F_2_%d",ilayer);
      map_BH[ilayer] = (TH2F*)_file1->Get(title);
    }
    if (ilayer>=40 && ilayer<=51) {
      sprintf(title,"map_TH2F_3_%d",ilayer);
      map_BH[ilayer] = (TH2F*)_file1->Get(title);
    }
  }
  
  //----- Plotting

  TEllipse *el1;
  TEllipse *el2;
  TEllipse *el3;
  double z;
  double eta_low=1.4;
  double r_max;
  double eta_hi;
  double r_min;
  
  int rmin_index[69];

  for (int ilayer=36;ilayer<=51;ilayer++){

    z = z_layer[ilayer+17];
    r_max=z*tan(2.*atan(exp(-eta_low)));

    if (ilayer<=40){
      eta_hi=eta_low+double(91)*0.01745;
      r_min=z*tan(2.*atan(exp(-eta_hi)));    
    }
    else {
      eta_hi=eta_low+double(73)*0.02182;
      r_min=z*tan(2.*atan(exp(-eta_hi)));    
    }
      
    el1 = new TEllipse(0.,0.,r_min,0.);
    el2 = new TEllipse(0.,0.,r_max,0.);
    el1->SetFillStyle(4000);
    el2->SetFillStyle(4000);
    
    std::cout << "r_min and r_max of that layer: " << r_min << " " << r_max << " " << ilayer+17 << std::endl;
    
    int ilayer_org=ilayer+17; // 0--68 convention
    // BH fine part
    if (ilayer>=36 && ilayer<=39){
      double eta=1.4+double(eta_index[ilayer_org])*0.01745;
      double z = z_layer[ilayer_org];
      double r = z*tan(2.*atan(exp(-eta)));    
      el3 = new TEllipse(0.,0.,r,0.);
      el3->SetFillStyle(4000);
      std::cout << "r,eta,eta_index of the Sinti inner boundary: "  << r << " " << eta << " " << eta_index[ilayer_org] << std::endl;

      double eta2=1.4+double(eta_index[ilayer_org]-2)*0.01745;
      double r2 = z*tan(2.*atan(exp(-eta2)));    

      std::cout << "r-r2: " << r-r2 << std::endl;      

    }
    // BH coarse part
    if (ilayer>=40 && ilayer<=51){
      double eta=1.4+double(eta_index[ilayer_org])*0.02182;
      double z = z_layer[ilayer_org];
      double r = z*tan(2.*atan(exp(-eta)));    
      el3 = new TEllipse(0.,0.,r,0.);
      el3->SetFillStyle(4000);
      std::cout << "r,eta,eta_index of the Sinti inner boundary: "  << r << " " << eta << " " << eta_index[ilayer_org] << std::endl;

      double eta2=1.4+double(eta_index[ilayer_org]-2)*0.02182;
      double r2 = z*tan(2.*atan(exp(-eta2)));    

      std::cout << "r-r2: " << r-r2 << std::endl;      
    }    
    
    map_Si[ilayer]->Draw("colz");
    map_BH[ilayer]->Draw("colz,pol,same");

    el2->SetLineColor(2);
    el3->SetLineColor(2);
    el1->SetLineColor(2);
    el2->Draw();
    el3->Draw();
    el1->Draw();
    map_Si[ilayer]->Draw("colz,same");
    map_BH[ilayer]->Draw("colz,pol,same");

    //
    // only for development
    // 
    if (develop) {
      double rmin=3000.;
      int eta_index;
      TH1D* proj = map_BH[ilayer]->ProjectionY();
      for (int ibinx=0;ibinx<=proj->GetNbinsX();ibinx++){
	double r = proj->GetBinCenter(ibinx);
	if (r<rmin && proj->GetBinContent(ibinx)){ rmin=r; eta_index=proj->GetNbinsX()-ibinx+1;}
      }      
      printf("%3d %f %d\n",ilayer,rmin,eta_index);
    }
      
    if (ilayer==36)      c1->Print("plot_geom.pdf(","pdf");
    else if (ilayer==51) c1->Print("plot_geom.pdf)","pdf");
    else                 c1->Print("plot_geom.pdf","pdf");

    sprintf(title,"plot_geom_layer%d.png",ilayer);
    c1->SaveAs(title);
    
  }

  */
    
}
