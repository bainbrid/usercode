{

  TFile outfile("efficiency.root");
  
  TCanvas c1;

  c1->SetGridy();

  int Xmin=-3;
  int Xmax =3;

  eff_eta_reco_M250_100pb_fast->SetMarkerSize(1);
  eff_eta_reco_M250_100pb_fast->SetMarkerStyle(20);
  eff_eta_reco_M250_100pb_fast->SetMarkerColor(kBlack);
  eff_eta_reco_M250_100pb_fast->GetXaxis()->SetRangeUser(Xmin,Xmax);
  eff_eta_reco_M250_100pb_fast->Draw("histep");

  eff_eta_reco_M250_100pb_fast->GetYaxis()->SetTitle("efficiency");

  eff_eta_reco_IDcut_M250_100pb_fast->SetMarkerSize(1);
  eff_eta_reco_IDcut_M250_100pb_fast->SetMarkerStyle(4);
  eff_eta_reco_IDcut_M250_100pb_fast->SetMarkerColor(kRed);
  eff_eta_reco_IDcut_M250_100pb_fast->Draw("samehistep");  

  TLegend *legend = new TLegend(0.614322,0.615385,0.874372,0.868881);
  legend->SetFillColor(kWhite);
  legend->AddEntry(eff_eta_reco_M250_100pb_fast,"gsf ele - no cut","p");
  legend->AddEntry(eff_eta_reco_IDcut_M250_100pb_fast,"gsf ele - ID cut","p");
  legend->Draw();

  //regions out of fiducial region

  Int_t n = 5;
  Double_t x[5] = {-3,-2,-2,-3,-3};
  Double_t y[5] = {0,0,1,1,0};
  TPolyLine * pl = new TPolyLine(n,x,y,"F");
  pl->SetFillStyle(3005);   
  pl->SetFillColor(kRed);   
  pl->Draw();

  Double_t x1[5] = {-1.6,-1.4,-1.4,-1.6,-1.6};
  Double_t y1[5] = {0,0,1,1,0};
  TPolyLine * pl1 = new TPolyLine(n,x1,y1,"F");
  pl1->SetFillStyle(3005);   
  pl1->SetFillColor(kRed);   
  pl1->Draw();

  Double_t x2[5] = {1.4,1.6,1.6,1.4,1.4};
  Double_t y2[5] = {0,0,1,1,0};
  TPolyLine * pl2 = new TPolyLine(n,x2,y2,"F");
  pl2->SetFillStyle(3005);   
  pl2->SetFillColor(kRed);   
  pl2->Draw();

  Double_t x3[5] = {2,3,3,2,2};
  Double_t y3[5] = {0,0,1,1,0};
  TPolyLine * pl3 = new TPolyLine(n,x3,y3,"F");
  pl3->SetFillStyle(3005);   
  pl3->SetFillColor(kRed);   
  pl3->Draw();



  //************************************************************

  TCanvas c2;
  eff_pt_reco_IDcut_M650_100pb_full->SetMarkerSize(1);
  eff_pt_reco_IDcut_M650_100pb_full->SetMarkerStyle(20);
  eff_pt_reco_IDcut_M650_100pb_full->SetMarkerColor(kBlack);
  eff_pt_reco_IDcut_M650_100pb_full->Draw("histep");

  eff_pt_reco_IDcut_M650_100pb_fast->SetFillColor(kBlue);
  eff_pt_reco_IDcut_M650_100pb_fast->SetFillStyle(3004);
  eff_pt_reco_IDcut_M650_100pb_fast->SetMarkerColor(kBlue);
  eff_pt_reco_IDcut_M650_100pb_fast->Draw("samehist");  

  TLegend *legend = new TLegend(0.614322,0.615385,0.874372,0.868881);
  legend->SetFillColor(kWhite);
  legend->AddEntry(eff_pt_reco_IDcut_M650_100pb_full,"full sim","p");
  legend->AddEntry(eff_pt_reco_IDcut_M650_100pb_fast,"fast sim","f");
  legend->Draw();


  TCanvas c3;
  eff_eta_reco_IDcut_M650_100pb_full->SetMarkerSize(1);
  eff_eta_reco_IDcut_M650_100pb_full->SetMarkerStyle(20);
  eff_eta_reco_IDcut_M650_100pb_full->SetMarkerColor(kBlack);
  eff_eta_reco_IDcut_M650_100pb_full->Draw("histep");

  eff_eta_reco_IDcut_M650_100pb_fast->SetFillColor(kBlue);
  eff_eta_reco_IDcut_M650_100pb_fast->SetFillStyle(3004);
  eff_eta_reco_IDcut_M650_100pb_fast->SetMarkerColor(kBlue);
  eff_eta_reco_IDcut_M650_100pb_fast->Draw("samehist");  

  TLegend *legend = new TLegend(0.614322,0.615385,0.874372,0.868881);
  legend->SetFillColor(kWhite);
  legend->AddEntry(eff_eta_reco_IDcut_M650_100pb_full,"full sim","p");
  legend->AddEntry(eff_eta_reco_IDcut_M650_100pb_fast,"fast sim","f");
  legend->Draw();


  
}
