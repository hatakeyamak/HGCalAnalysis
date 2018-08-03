{

  TFile *_file0 = TFile::Open("et100_eta1p7_v4.root");
  TFile *_file1 = TFile::Open("et100_eta1p7_v6.root");

  gStyle->SetOptStat(0);
  
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
    
}
