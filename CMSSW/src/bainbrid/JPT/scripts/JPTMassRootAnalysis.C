#define JPTRootAnalysis_cxx
#include "JPTMassRootAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

Float_t JPTRootAnalysis::deltaPhi(Float_t phi1, Float_t phi2)
{
  Float_t pi = 3.1415927;
  Float_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Float_t JPTRootAnalysis::deltaEta(Float_t eta1, Float_t eta2)
{
  Float_t deta = fabs(eta1-eta2);
  return deta;
}

Float_t JPTRootAnalysis::deltaR(Float_t eta1, Float_t eta2,
		                Float_t phi1, Float_t phi2)
{
  Float_t dr = sqrt( deltaEta(eta1, eta2) * deltaEta(eta1, eta2) +
		     deltaPhi(phi1, phi2) * deltaPhi(phi1, phi2) );
  return dr;
}

void JPTRootAnalysis::setTDRStyle(Int_t ylog) {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(2);

  tdrStyle->SetEndErrorSize(4);
  
  tdrStyle->SetMarkerStyle(20);

  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  tdrStyle->SetOptDate(0);

  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);


  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);


  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(ylog);
  tdrStyle->SetOptLogz(0);

  tdrStyle->SetPaperSize(15.,15.);

  tdrStyle->cd();
}

void JPTRootAnalysis::Loop()
{
   if (fChain == 0) return;

   const Int_t nbins = 10;
   const Int_t nh = 10;

   const Int_t nbx = nbins+1;
   const Float_t xbins[nbx]={20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,150.};
   TH1F* hResRaw = new TH1F("hResRaw", "ResRaw", nbins, xbins);
   TH1F* hResJPTInCone = new TH1F("hResJPTInCone", "ResJPTInCone", nbins, xbins);
   TH1F* hResZSP = new TH1F("hResZSP", "ResZSP", nbins, xbins);
   TH1F* hResJPT = new TH1F("hResJPT", "ResJPT", nbins, xbins);
   TH1F* hScaleRaw = new TH1F("hScaleRaw", "ScaleRaw", nbins, xbins);
   TH1F* hScaleJPTInCone = new TH1F("hScaleJPTInCone", "ScaleJPTInCone", nbins, xbins);
   TH1F* hScaleZSP = new TH1F("hScaleZSP", "ScaleZSP", nbins, xbins);
   TH1F* hScaleJPT = new TH1F("hScaleJPT", "ScaleJPT", nbins, xbins);

   TH1F* hEtRaw[nh];
   TH1F* hEtJPTInCone[nh]; 
   TH1F* hEtZSP[nh];
   TH1F* hEtJPT[nh];

   const char* namesEtRaw[nh] = {"hEtRaw1","hEtRaw2","hEtRaw3","hEtRaw4","hEtRaw5","hEtRaw6","hEtRaw7","hEtRaw8","hEtRaw9","hEtRaw10"};
   const char* titleEtRaw[nh] = {"EtRaw1","EtRaw2","EtRaw2","EtRaw4","EtRaw5","EtRaw6","EtRaw7","EtRaw8","EtRaw9","EtRaw10"};

   const char* namesEtJPTInCone[nh] = {"hEtJPTInCone1","hEtJPTInCone2","hEtJPTInCone3","hEtJPTInCone4","hEtJPTInCone5","hEtJPTInCone6","hEtJPTInCone7","hEtJPTInCone8","hEtJPTInCone9","hEtJPTInCone10"};
   const char* titleEtJPTInCone[nh] = {"EtJPTInCone1","EtJPTInCone2","EtJPTInCone3","EtJPTInCone4","EtJPTInCone5","EtJPTInCone6","EtJPTInCone7","EtJPTInCone8","EtJPTInCone9","EtJPTInCone10"};

   const char* namesEtZSP[nh] = {"hEtZSP1","hEtZSP2","hEtZSP3","hEtZSP4","hEtZSP5","hEtZSP6","hEtZSP7","hEtZSP8","hEtZSP9","hEtZSP10"};
   const char* titleEtZSP[nh] = {"EtZSP1","EtZSP2","EtZSP3","EtZSP4","EtZSP5","EtZSP6","EtZSP7","EtZSP8","EtZSP9","EtZSP10"};

   const char* namesEtJPT[nh] = {"hEtJPT1","hEtJPT2","hEtJPT3","hEtJPT4","hEtJPT5","hEtJPT6","hEtJPT7","hEtJPT8","hEtJPT9","hEtJPT10"};
   const char* titleEtJPT[nh] = {"EtJPT1","EtJPT2","EtJPT3","EtJPT4","EtJPT5","EtJPT6","EtJPT7","EtJPT8","EtJPT9","EtJPT10"};

   for(Int_t ih=0; ih < nh; ih++) { 
     hEtRaw[ih]  = new TH1F(namesEtRaw[ih], titleEtRaw[ih], 60, 0., 3.);
     hEtJPTInCone[ih]  = new TH1F(namesEtJPTInCone[ih], titleEtJPTInCone[ih], 60, 0., 3.);
     hEtZSP[ih]  = new TH1F(namesEtZSP[ih], titleEtZSP[ih], 60, 0., 3.);
     hEtJPT[ih]  = new TH1F(namesEtJPT[ih], titleEtJPT[ih], 60, 0., 3.);
   }

   TH1F * hEtGen  = new TH1F( "hEtGen", "EtGen", 20, 0., 200.);
   TH1F * hEtaGen = new TH1F( "hEtaGen", "EtaGen", 16, 0., 2.1);
   TH1F * hDR     = new TH1F( "hDR", "DR", 100, 0., 10.);

   Float_t DR;
   Float_t DRcut = 2.0;
   Float_t etaMin = 0.0;
   Float_t etaMax = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for(Int_t ih = 0; ih < nh; ih++) {
	if(EtGen1 >= xbins[ih] && EtGen1 < xbins[ih+1] 
	   && DRMAXgjet1 < 0.30
	   && fabs(EtaGen1) > etaMin && fabs(EtaGen1) <= etaMax) {
	  if(EtRaw1/EtGen1 > 0.1) {


	    if(EtGen2 < 20.) {
	      hEtRaw[ih]->Fill(MassRaw1/MassGen1);
	      hEtZSP[ih]->Fill(MassZSP1/MassGen1);
	      hEtJPT[ih]->Fill(MassJPT1/MassGen1);
	    } 

	    if(EtGen2 > 20.) {
	      DR = deltaR(EtaGen1, EtaGen2, PhiGen1, PhiGen2);
	      if(DR > DRcut) {
		hDR->Fill(DR);
		hEtRaw[ih]->Fill(MassRaw1/MassGen1);
		hEtZSP[ih]->Fill(MassZSP1/MassGen1);
		hEtJPT[ih]->Fill(MassJPT1/MassGen1);
		hEtGen->Fill(EtGen1);
		hEtaGen->Fill(EtaGen1);
	      }
	    }
	  }
	}
	if(EtGen2 >= xbins[ih] && EtGen2 < xbins[ih+1] 
	   && DRMAXgjet2 < 0.30
	   && fabs(EtaGen2) > etaMin && fabs(EtaGen2) <= etaMax) {
	  if(EtRaw2/EtGen2 > 0.1) {

	    DR = deltaR(EtaGen1, EtaGen2, PhiGen1, PhiGen2);
	    if(DR > DRcut) {
	      hDR->Fill(DR);
	      hEtRaw[ih]->Fill(MassRaw2/MassGen2);
	      hEtZSP[ih]->Fill(MassZSP2/MassGen2);
	      hEtJPT[ih]->Fill(MassJPT2/MassGen2);
	      hEtGen->Fill(EtGen2);
	      hEtaGen->Fill(EtaGen2);
	    }
	  }
	}
      }
   }

   setTDRStyle(0);
   gStyle->SetOptFit();

  tdrStyle->SetOptLogy(1);

   TCanvas* c1 = new TCanvas("X1","Y1",1);
   c1->Divide(2,2);
   char name[50];
   for(Int_t ih = 0; ih < nh; ++ih) {

     TAxis* xaxisRes = hResJPT->GetXaxis();
     Double_t EbinCenter = xaxisRes->GetBinCenter(ih+1);
     cout <<" bin center = " << EbinCenter << endl;

     Double_t mean = 1000.;
     Double_t meanErr = 1000.;
     Double_t sigma = 1000.;
     Double_t sigmaErr = 1000.;
     Float_t resolution = 1000.;
     Float_t resolutionErr = 1000.;
     
     c1->cd(1);
     Int_t binMax = hEtJPT[ih]->GetMaximumBin();
     TAxis* xaxis = hEtJPT[ih]->GetXaxis();
     Double_t binCenter = xaxis->GetBinCenter(binMax);
     Double_t rms = hEtJPT[ih]->GetRMS();
     Double_t rFitMin = binCenter - 2.0 * rms; 
     Double_t rFitMax = binCenter + 2.0 * rms;
     hEtJPT[ih]->Fit("gaus","","",rFitMin,rFitMax);
     TF1 *fit = hEtJPT[ih]->GetFunction("gaus"); 
     gStyle->SetOptFit();
     mean = 1000.;
     meanErr = 1000.;
     sigma = 1000.;
     sigmaErr = 1000.;
     resolution = 1000.;
     resolutionErr = 1000.;
     if ( fit ) {
       mean  = fit->GetParameter(1);
       meanErr  = fit->GetParError(1);
       sigma = fit->GetParameter(2);
       sigmaErr = fit->GetParError(2);
       if ( mean > 0. ) { resolution = sigma/mean; }
       if ( mean > 0. && sigma > 0. ) {
	 resolutionErr = resolution * sqrt((meanErr/mean)*(meanErr/mean) + (sigmaErr/sigma)*(sigmaErr/sigma));
       }
     }
     if ( resolution < 999. ) hResJPT->Fill(EbinCenter,resolution);
     if ( resolutionErr < 999. ) hResJPT->SetBinError(ih+1,resolutionErr);    
     if ( mean < 999. ) hScaleJPT->Fill(EbinCenter,mean);
     if ( meanErr < 999. ) hScaleJPT->SetBinError(ih+1,meanErr);    

     c1->cd(2);
     binMax = hEtZSP[ih]->GetMaximumBin();
     xaxis = hEtZSP[ih]->GetXaxis();
     binCenter = xaxis->GetBinCenter(binMax);
     rms = hEtZSP[ih]->GetRMS();
     rFitMin = binCenter - 2.0 * rms; 
     rFitMax = binCenter + 2.0 * rms;
     hEtZSP[ih]->Fit("gaus","","",rFitMin,rFitMax);
     fit = hEtZSP[ih]->GetFunction("gaus"); 
     mean = 1000.;
     meanErr = 1000.;
     sigma = 1000.;
     sigmaErr = 1000.;
     resolution = 1000.;
     resolutionErr = 1000.;
     if ( fit ) {
       mean  = fit->GetParameter(1);
       meanErr  = fit->GetParError(1);
       sigma = fit->GetParameter(2);
       sigmaErr = fit->GetParError(2);
       if ( mean > 0. ) { resolution = sigma/mean; }
       if ( mean > 0. && sigma > 0. ) {
	 resolutionErr = resolution * sqrt((meanErr/mean)*(meanErr/mean) + (sigmaErr/sigma)*(sigmaErr/sigma));
       }
     }
     if ( resolution < 999. ) hResZSP->Fill(EbinCenter,resolution);
     if ( resolutionErr < 999. ) hResZSP->SetBinError(ih+1,resolutionErr);    
     if ( mean < 999. ) hScaleZSP->Fill(EbinCenter,mean);
     if ( meanErr < 999. ) hScaleZSP->SetBinError(ih+1,meanErr);    

     c1->cd(4);
     binMax = hEtRaw[ih]->GetMaximumBin();
     xaxis = hEtRaw[ih]->GetXaxis();
     binCenter = xaxis->GetBinCenter(binMax);
     rms = hEtRaw[ih]->GetRMS();
     rFitMin = binCenter - 2.0 * rms; 
     rFitMax = binCenter + 2.0 * rms;
     hEtRaw[ih]->Fit("gaus","","",rFitMin,rFitMax);
     fit = hEtRaw[ih]->GetFunction("gaus"); 
     mean = 1000.;
     meanErr = 1000.;
     sigma = 1000.;
     sigmaErr = 1000.;
     resolution = 1000.;
     resolutionErr = 1000.;
     if ( fit ) {
       mean  = fit->GetParameter(1);
       meanErr  = fit->GetParError(1);
       sigma = fit->GetParameter(2);
       sigmaErr = fit->GetParError(2);
       if ( mean > 0. ) { resolution = sigma/mean; }
       if ( mean > 0. && sigma > 0. ) {
	 resolutionErr = resolution * sqrt((meanErr/mean)*(meanErr/mean) + (sigmaErr/sigma)*(sigmaErr/sigma));
       }
     }
     if ( resolution < 999. ) hResRaw->Fill(EbinCenter,resolution);
     if ( resolutionErr < 999. ) hResRaw->SetBinError(ih+1, resolutionErr);    
     if ( mean < 999. ) hScaleRaw->Fill(EbinCenter,mean);
     if ( meanErr < 999. ) hScaleRaw->SetBinError(ih+1,meanErr);
     sprintf(name,"hCalo1_%d.eps",ih);
     c1->SaveAs(name);
   }

   tdrStyle->SetOptLogy(0);

   TFile efile("test.root","recreate");
   hResRaw->Write();
   hResJPTInCone->Write();
   hResJPT->Write();
   hScaleRaw->Write();
   hScaleJPTInCone->Write();
   hScaleJPT->Write();
   efile.Close();

   TCanvas* c40 = new TCanvas("X2","Y2",1);

   hResJPT->GetXaxis()->SetTitle("E_{T} Gen, GeV");
   hResJPT->GetYaxis()->SetTitle("Energy resolution, % ");

   hResJPT->SetMaximum(0.50);
   hResJPT->SetMinimum(0.05);
   hResJPT->SetMarkerStyle(21);
   hResJPT->SetMarkerSize(1.2);
   hResJPT->Draw("histPE1");

   hResZSP->SetMarkerSize(1.0);
   hResZSP->SetMarkerStyle(24);
   hResZSP->Draw("samePE1");

   hResRaw->SetMarkerSize(1.5);
   hResRaw->SetMarkerStyle(22);
   hResRaw->Draw("samePE1");

   TLatex *t = new TLatex();
   t->SetTextSize(0.042);
   TLegend *leg = new TLegend(0.45,0.5,0.85,0.8,NULL,"brNDC");
   leg->SetFillColor(10);
   leg->AddEntry(hResRaw,"Raw calo jets","P");
   leg->AddEntry(hResZSP,"ZSP corr","P");
   leg->AddEntry(hResJPT,"ZSP+JPT corr","P");
   leg->Draw();  
   t->DrawLatex(25,0.42,"CMSSW219");
   t->DrawLatex(25,0.40,"RelVal QCD 80-120 GeV, |#eta ^{jet}|< 1.0");

   c40->SaveAs("resJPT219.eps");

   TCanvas* c20 = new TCanvas("X3","Y3",1);

   hScaleJPT->GetXaxis()->SetTitle("E_{T} Gen, GeV");
   hScaleJPT->GetYaxis()->SetTitle("E_{T}^{reco}/E_{T}^{gen}");

   hScaleJPT->SetMaximum(1.2);
   hScaleJPT->SetMinimum(0.2);
   hScaleJPT->SetMarkerStyle(21);
   hScaleJPT->SetMarkerSize(1.2);
   hScaleJPT->Draw("histPE1");

   hScaleZSP->SetMarkerSize(1.0);
   hScaleZSP->SetMarkerStyle(24);
   hScaleZSP->Draw("samePE1");

   hScaleRaw->SetMarkerSize(1.5);
   hScaleRaw->SetMarkerStyle(22);
   hScaleRaw->Draw("samePE1");

   TLatex *t = new TLatex();
   t->SetTextSize(0.042);
   TLegend *leg = new TLegend(0.5,0.15,0.9,0.35,NULL,"brNDC");
   leg->SetFillColor(10);
   leg->AddEntry(hScaleRaw,"Raw calo jets","P");
   leg->AddEntry(hScaleZSP,"ZSP corr","P");
   leg->AddEntry(hScaleJPT,"ZSP+JPT corr","P");
   leg->Draw();  
   t->DrawLatex(25,1.12,"CMSSW219");
   t->DrawLatex(25,1.06,"RelVal QCD 80-120 GeV, |#eta ^{jet}|< 1.0");

   c20->SaveAs("ScaleJPT219.eps");
}
