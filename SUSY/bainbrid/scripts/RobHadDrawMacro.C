HadDrawMacro()
{
  TFile* myOutPut = new TFile("OverlayedPlots.root","recreate");
  // example long
  TString name("cosTall"); 
  TString dirmame("AfterHTPlots");
  int rebin = 3;
  bool normalized = true;
  bool logerithmic = false;


  aDrawBkgdPlots("MuoPlusEleMultiplicity_Gen_soft","countplots",1,normalized,logerithmic,myOutPut,true);


  myOutPut->Write();
  myOutPut->Close();
}


TCanvas* aDrawBkgdPlots(TString name, TString dirmame, int rebin, bool normalized, bool logerithmic,TDirectory* myOutPut,bool plotSumAll=false)
{
  SetSomeStyles();
  TLegend *myLegend = new TLegend(0.6084746,0.5070671,0.9186441,0.9293286,NULL,"brNDC");
 
  myLegend->SetFillColor(0);
  myLegend->SetLineColor(0); 
 

  TCanvas *aCanvas = getaCanvas(name,myOutPut);
    

  TH1D* lm1 = readHist(name,"lm1_MHT/lm1.root",dirmame,rebin);
 
  TH1D* qcd_100_250 = readHist(name,"qcd_MHT/qcd_100_250.root",dirmame,rebin);
  TH1D* qcd_250_500 = readHist(name,"qcd_fast_MHT/qcd_250_500_fast.root",dirmame,rebin);
  TH1D* qcd_500_1000 = readHist(name,"qcd_MHT/qcd_500_1000.root",dirmame,rebin);
  TH1D* qcd_1000_Inf = readHist(name,"qcd_MHT/qcd_1000_Inf.root",dirmame,rebin);
  TH1D* w = readHist(name,"W_MHT/wjets.root",dirmame,rebin);
  TH1D* Zinv = readHist(name,"Zinv_MHT/zinv.root",dirmame,rebin);
  TH1D* tt = readHist(name,"tt_MHT/tt.root",dirmame,rebin);
  TH1D* Zjets = readHist(name,"ZJets_MHT/zjets.root",dirmame,rebin);



  TH1D* allQCD = qcd_100_250->Clone();
  allQCD->Add(qcd_250_500,1);
  allQCD->Add(qcd_500_1000,1);
  allQCD->Add(qcd_1000_Inf,1); 


  TH1D* Z = Zinv->Clone();
  Z->Add(Zjets,1);

 
  allQCD->SetLineColor(kGreen+2);
  allQCD->SetFillColor(kGreen+2);
  allQCD->SetFillStyle(3005);
  
  lm1->SetLineColor(kRed);
  w->SetLineColor(kMagenta);
  Zinv->SetLineColor(kBlack);
  Zjets->SetLineColor(kOrange+7);
  tt->SetLineColor(kBlue);


  myLegend->AddEntry(lm1,"LM1", "L");
  myLegend->AddEntry(tt,"t#bar{t}", "L");
  myLegend->AddEntry(w,"W", "L");
  myLegend->AddEntry(Z,"Z", "L");
  myLegend->AddEntry(allQCD,"QCD", "f");

  double aMax =0.;
  if(lm1->GetMaximum()>aMax) aMax=lm1->GetMaximum();  
  if(w->GetMaximum()>aMax) aMax=w->GetMaximum();  
  if(Zinv->GetMaximum()>aMax) aMax=Zinv->GetMaximum();  
  if(tt->GetMaximum()>aMax) aMax=tt->GetMaximum();  
  if(allQCD->GetMaximum()>aMax) aMax=allQCD->GetMaximum();
 
  allQCD->GetYaxis()->SetTitleOffset(1.43);
  allQCD->GetYaxis()->SetTitleSize(0.06);
  allQCD->GetXaxis()->SetTitleSize(0.06);
  allQCD->GetXaxis()->SetTitleOffset(0.9);

  allQCD->SetMaximum(aMax*1.25);
  allQCD->SetMinimum(0.1);
  
  if(normalized == true){
    
   
    allQCD->DrawNormalized("Ehist");
  
    tt->DrawNormalized("hsame");
    Z->DrawNormalized("hsame");
    w->DrawNormalized("hsame");
   
    lm1->DrawNormalized("hsame");
   
   
  }
  else {
    allQCD->Draw("h");
    w->Draw("sameH");
    Z->Draw("sameH");
   
    tt->Draw("sameh");
    lm1->Draw("sameH");
   
  }
  myOutPut->cd();

  myLegend->Draw("same");

 
 
  aCanvas->Write();

  return aCanvas;
}







TH1* readHist(TString nameHist,TString nameFile,TString Dirname, int rebin)
{
  TFile* file =  new TFile(nameFile);
  // file->ls();
  TDirectory* dir = (TDirectory*)file->Get(Dirname);
  // dir->ls();
  // cout << " name hist " << nameHist << " name file " << nameFile << endl;
  TH1* hist = (TH1*)dir->Get(nameHist);
  hist->SetLineWidth(3);
  // if(rebin>0) hist->Rebin(rebin);
  hist->GetXaxis()->SetTitleSize(.055);
  hist->GetYaxis()->SetTitleSize(.055);
  hist->GetXaxis()->SetLabelSize(.05);
  hist->GetYaxis()->SetLabelSize(.05);
  hist->SetStats(kFALSE);

  return hist;
}

TCanvas* getaCanvas(TString name,TDirectory* afile)
{
  afile->cd();
  TCanvas* aCanvas = new TCanvas(name);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  aCanvas->Range(-288.2483,-2.138147,1344.235,6.918939);
  aCanvas->SetFillColor(0);
  aCanvas->SetBorderMode(0);
  aCanvas->SetBorderSize(2);
  if(logerithmic == true)aCanvas->SetLogy();
  aCanvas->SetLeftMargin(0.1765705);
  aCanvas->SetRightMargin(0.05772496);
  aCanvas->SetTopMargin(0.04778761);
  aCanvas->SetBottomMargin(0.1256637);
  aCanvas->SetFrameFillStyle(0);
  aCanvas->SetFrameLineWidth(2);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameFillStyle(0);
  aCanvas->SetFrameLineWidth(2);
  aCanvas->SetFrameBorderMode(0);
 
  
  return aCanvas;
}


