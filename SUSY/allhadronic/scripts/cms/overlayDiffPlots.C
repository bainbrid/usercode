#include "style.C"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TText.h>
#include <vector>
#include <iostream>
#include <sstream>

typedef unsigned int uint;

// -----------------------------------------------------------------------------
//
TCanvas* createCanvas(TString name,TDirectory* afile, bool log)
{
  afile->cd();
  TCanvas* aCanvas = new TCanvas(name);
  //gStyle->SetOptFit(1);
  //gStyle->SetOptStat("mr");
  aCanvas->Range(-288.2483,-2.138147,1344.235,6.918939);
  aCanvas->SetFillColor(0);
  aCanvas->SetBorderMode(0);
  aCanvas->SetBorderSize(2);
  if ( log == true)aCanvas->SetLogy();
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

// -----------------------------------------------------------------------------
//
TH1* getHisto( TString nameFile,
	       TString nameHist,
	       TString Dirname, 
	       int rebin ) {
  TString name = nameFile;
  TFile* file =  new TFile(name);
  if (file) { std::cout << "Opened file: " << file->GetName() << std::endl; }
  else { 
    std::cout << "Could not find file: " << name << std::endl; 
    return 0; 
  }
  
  TDirectory* dir = (TDirectory*)file->Get(Dirname);
  if (dir) { std::cout << "Opened dir: " << dir->GetName() << std::endl; }
  else { 
    std::cout << "Could not find dir: " << Dirname << std::endl; 
    return 0; 
  }
  
  TH1* hist = (TH1*)dir->Get(nameHist);
  if (hist) { std::cout << "Opened histo: " << hist->GetName() << std::endl; }
  else { 
    std::cout << "Could not find histo: " << nameHist << std::endl; 
    return 0; 
  }

  hist->SetLineWidth(3);
  if ( rebin > 0 ) { hist->Rebin(rebin); }
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->SetStats(kFALSE);
  return hist;
}

// -----------------------------------------------------------------------------
//
TCanvas* createPlot( std::vector<TString> paths, 
		     TString canvas_name, 
		     std::vector<TString> names, 
		     std::vector<TString> titles, 
		     TString dirname, 
		     double lumi,
		     int rebin, 
		     bool norm, 
		     bool log,
		     TDirectory* file,
		     double min = -1.,
		     double max = -1. )
{

  // SetSomeStyles();

  // Create legend
  TLegend* legend = new TLegend( 0.6, 0.5, 0.9, 0.7,  NULL, "brNDC" );
  legend->SetFillColor(0);
  legend->SetLineColor(0); 
  
  // Create canvas
  TCanvas* aCanvas = createCanvas( canvas_name, file, log );

  TPaveText* stats = new TPaveText( 0.6, 0.3, 0.9, 0.45, "brNDC" );
  stats->SetFillColor(0);
  stats->SetLineColor(0);

  TLatex* prelim = new TLatex( 0.55, 0.8, "#scale[0.8]{CMS preliminary 2010}" );
  prelim->SetTextSize(0.04);
  prelim->SetNDC();
  std::stringstream ssl; ssl << "#scale[0.8]{#int L dt = " << lumi << "pb^{-1}, #sqrt{s} = 7 TeV}";
  TLatex* lumis = new TLatex( 0.51, 0.87, ssl.str().c_str() );
  lumis->SetTextSize(0.04);
  lumis->SetNDC();
    
  // First histo be drawn
  bool first = true;
  
  // Loop through histogram names
  std::vector<TH1D*> his;
  for ( uint iname = 0; iname < names.size(); ++iname ) {
    his.push_back( (TH1D*)getHisto( paths[iname], names[iname], dirname, rebin ) );
  }
  
  // Loop through histograms
  uint colour = 0;
  double aMax = 0.;
  double aMin = 1.e12;
  for ( uint ihis = 0; ihis < his.size(); ++ihis ) {
    if ( !his[ihis] ) { continue; }
    
    // Line colour and fill
    colour++;
    if ( colour==3 ) { colour++; }
    if ( colour==5 ) { colour++; }
    his[ihis]->SetMarkerStyle(21);
    his[ihis]->SetMarkerColor(colour);
    his[ihis]->SetMarkerSize(0.8);
    his[ihis]->SetLineColor(colour);
    his[ihis]->SetLineStyle(1);
    his[ihis]->SetLineWidth(1);
    
    // Populate legend
    if ( titles.size() > ihis ) { legend->AddEntry( his[ihis], titles[ihis], "EPL" ); }
    else                        { legend->AddEntry( his[ihis], names[ihis], "EPL" ); }

    // Populate stats box
    std::stringstream ss;
    ss << "Mean=" << int(his[ihis]->GetMean()*100.)/100. << ", RMS=" << int(his[ihis]->GetRMS()*100.)/100.;
    TText* text = stats->AddText(ss.str().c_str());
    text->SetTextSize(0.03);
    text->SetTextColor(colour);
    
    // Calc min/max number of entries
    if ( his[ihis]->GetMaximum() > aMax ) { aMax = his[ihis]->GetMaximum(); }
    if ( his[ihis]->GetMinimum(1.e-12) < aMin ) { aMin = his[ihis]->GetMinimum(1.e-12); }

  }
  
  if ( !his.empty() ) {
    if ( his[0] ) his[0]->GetYaxis()->SetTitleOffset(1.43);
    if ( his[0] ) his[0]->GetYaxis()->SetTitleSize(0.06);
    if ( his[0] ) his[0]->GetXaxis()->SetTitleSize(0.06);
    if ( his[0] ) his[0]->GetXaxis()->SetTitleOffset(0.9);
  }

  for ( uint ihis = 0; ihis < his.size(); ++ihis ) {

    if ( !his[ihis] ) { continue; }
    
    //his[ihis]->GetYaxis()->SetTitle("a.u.");
    
    if ( log ) {
      his[ihis]->SetMaximum( aMax * 10. );
      his[ihis]->SetMinimum( aMin * 0.1 );
    } else {
      his[ihis]->SetMaximum( aMax * 1.1 );
      his[ihis]->SetMinimum( aMin * 0.9 );
    }

    if ( min > 0. ) his[ihis]->SetMinimum( min );
    if ( max > 0. ) his[ihis]->SetMaximum( max );
    
    if ( norm ) {
      TString options = "";
      if ( first ) { options = "Ehist"; first = false; }
      else { options = "hsame"; }
      if ( his[ihis]->GetEntries() > 0. ) { his[ihis]->DrawNormalized(options); }
    } else {
      TString options = "";
      if ( first ) { options = "h"; first = false; }
      else { options = "hsame"; }
      his[ihis]->Draw(options);
    }
    
  } // Loop through histos

  file->cd();
  legend->Draw("same");
  stats->Draw("same");
  prelim->Draw("same");
  lumis->Draw("same");
  aCanvas->Modified();
  aCanvas->SaveAs( std::string(canvas_name+".png").c_str() );
  //aCanvas->SaveAs( std::string(canvas_name+".C").c_str() );
  aCanvas->Write();
  return aCanvas;

}

// -----------------------------------------------------------------------------
//
TCanvas* createPlot( TString path, 
		     TString canvas_name, 
		     std::vector<TString> names, 
		     std::vector<TString> titles, 
		     TString dirname, 
		     double lumi,
		     int rebin, 
		     bool norm, 
		     bool log,
		     TDirectory* file,
		     double min = -1.,
		     double max = -1. )
{
  std::vector<TString> paths;
  paths.push_back(path);
  return createPlot( paths, 
		     canvas_name, 
		     names, 
		     titles, 
		     dirname, 
		     lumi,
		     rebin, 
		     norm, 
		     log,
		     file,
		     min,
		     max );
}

// -----------------------------------------------------------------------------
//
TCanvas* createPlot( TString path, 
		     TString canvas_name, 
		     TString name, 
		     TString title, 
		     TString dirname, 
		     double lumi,
		     int rebin, 
		     bool norm, 
		     bool log,
		     TDirectory* file,
		     double min = -1.,
		     double max = -1. )
{
  std::vector<TString> names;
  std::vector<TString> titles;
  names.push_back(name);
  titles.push_back(title);
  return createPlot( path, 
		     canvas_name, 
		     names, 
		     titles, 
		     dirname, 
		     lumi,
		     rebin, 
		     norm, 
		     log,
		     file,
		     min,
		     max );
}

// -----------------------------------------------------------------------------
/*
  Arguments to createPlot():
  - path to input files
  - canvas title
  - histogram title
  - histogram directory
  - rebin histograms?
  - normalize all histograms to unit area?
  - logarithmic y-scale
  - output file
  - min for y-axis (optiional)
  - max for y-axis (optiional)
*/  
int overlayDiffPlots() {

  setTDRStyle();

  double lumi = 15.1;

  // Path to input files
  TString path("../../../results/batch/101102_1/100_100_50/");

  // Output file
  TFile* output = new TFile( "Plots.root", "RECREATE" );
  if ( !output || output->IsZombie() ) { return -1; }

  // -------------------- Jet multiplicities --------------------
   
  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_QCD", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_QCD", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_Data", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_Data", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_SM", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_SM", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }


  // ---------- diff ht bins ----------

  //std::string histo = "BabyCaloOverMeffAfterDeadEcal";
  //   std::string histo = "BabyPfOverMeffAfterDeadEcal";
  std::string histo = "BabyMhtOverMetAfterDeadEcal";
  //    std::string histo = "GenBabyOverMeffNoMet";
  //   std::string histo = "GenBabyMhtOverMetNoMet";
  
  std::string bin0 = "HT250";
   std::string bin1 = "HT300";
   std::string bin2 = "HT350";

   //    std::string bin0 = "Meff400";
   //    std::string bin1 = "Meff450";
   //    std::string bin2 = "Meff500";

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.5_all");
    names.push_back(histo+"_"+bin0+"_aT0.51_all");
    names.push_back(histo+"_"+bin0+"_aT0.52_all");
    names.push_back(histo+"_"+bin0+"_aT0.53_all");
    names.push_back(histo+"_"+bin0+"_aT0.54_all");
    names.push_back(histo+"_"+bin0+"_aT0.55_all");
    createPlot( path, histo+"_HT250", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin1+"_aT0.5_all");
    names.push_back(histo+"_"+bin1+"_aT0.51_all");
    names.push_back(histo+"_"+bin1+"_aT0.52_all");
    names.push_back(histo+"_"+bin1+"_aT0.53_all");
    names.push_back(histo+"_"+bin1+"_aT0.54_all");
    names.push_back(histo+"_"+bin1+"_aT0.55_all");
    createPlot( path, histo+"_HT300", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin2+"_aT0.5_all");
    names.push_back(histo+"_"+bin2+"_aT0.51_all");
    names.push_back(histo+"_"+bin2+"_aT0.52_all");
    names.push_back(histo+"_"+bin2+"_aT0.53_all");
    names.push_back(histo+"_"+bin2+"_aT0.54_all");
    names.push_back(histo+"_"+bin2+"_aT0.55_all");
    createPlot( path, histo+"_HT350", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.51_all");
    names.push_back(histo+"_"+bin1+"_aT0.51_all");
    names.push_back(histo+"_"+bin2+"_aT0.51_all");
    createPlot( path, histo+"_Check0.51", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.52_all");
    names.push_back(histo+"_"+bin1+"_aT0.52_all");
    names.push_back(histo+"_"+bin2+"_aT0.52_all");
    createPlot( path, histo+"_Check0.52", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.53_all");
    names.push_back(histo+"_"+bin1+"_aT0.53_all");
    names.push_back(histo+"_"+bin2+"_aT0.53_all");
    createPlot( path, histo+"_Check0.53", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.54_all");
    names.push_back(histo+"_"+bin1+"_aT0.54_all");
    names.push_back(histo+"_"+bin2+"_aT0.54_all");
    createPlot( path, histo+"_Check0.54", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.55_all");
    names.push_back(histo+"_"+bin1+"_aT0.55_all");
    names.push_back(histo+"_"+bin2+"_aT0.55_all");
    createPlot( path, histo+"_Check0.55", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }


  // --------------------------------------------------------------------------------

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.5_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.51_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.52_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.53_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.54_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.55_all");
    createPlot( path, "GenBabyOverHtNoMet", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    //names.push_back("GenHtMultiplicity_HT0");
    names.push_back("GenHtMultiplicity_HT250");
    names.push_back("GenHtMultiplicity_HT300");
    names.push_back("GenHtMultiplicity_HT350");
    createPlot( path, "GenHtMultiplicity", names, titles, "QcdBkgdEst", lumi, 1, false, true, output );
  }
  
  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_Ht0_aT0_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.5_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.51_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.52_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.53_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.54_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.55_all");
    createPlot( path, "GenBabyOverHt_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_HT250_aT0_all");
    names.push_back("GenBabyOverHt_HT350_aT0_all");
    names.push_back("GenBabyOverHt_HT250_aT0.5_all");
    names.push_back("GenBabyOverHt_HT350_aT0.5_all");
    names.push_back("GenBabyOverHt_HT250_aT0.51_all");
    names.push_back("GenBabyOverHt_HT350_aT0.51_all");
    names.push_back("GenBabyOverHt_HT250_aT0.52_all");
    names.push_back("GenBabyOverHt_HT350_aT0.52_all");
    createPlot( path, "GenBabyOverHt_Vs_Ht_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0_all");
    names.push_back("BabyOverHt_Ht0_aT0.5_all");
    names.push_back("BabyOverHt_Ht0_aT0.51_all");
    names.push_back("BabyOverHt_Ht0_aT0.52_all");
    names.push_back("BabyOverHt_Ht0_aT0.53_all");
    names.push_back("BabyOverHt_Ht0_aT0.54_all");
    names.push_back("BabyOverHt_Ht0_aT0.55_all");
    createPlot( path, "BabyOverHt_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_HT250_aT0_all");
    names.push_back("BabyOverHt_HT350_aT0_all");
    names.push_back("BabyOverHt_HT250_aT0.5_all");
    names.push_back("BabyOverHt_HT350_aT0.5_all");
    names.push_back("BabyOverHt_HT250_aT0.51_all");
    names.push_back("BabyOverHt_HT350_aT0.51_all");
    names.push_back("BabyOverHt_HT250_aT0.52_all");
    names.push_back("BabyOverHt_HT350_aT0.52_all");
    createPlot( path, "BabyOverHt_Vs_Ht_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0_all");
    names.push_back("BabyOverHt_HT250_aT0_all");
    names.push_back("BabyOverHt_HT300_aT0_all");
    names.push_back("BabyOverHt_HT350_aT0_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0.51_all");
    names.push_back("BabyOverHt_HT250_aT0.51_all");
    names.push_back("BabyOverHt_HT300_aT0.51_all");
    names.push_back("BabyOverHt_HT350_aT0.51_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0.52_all");
    names.push_back("BabyOverHt_HT250_aT0.52_all");
    names.push_back("BabyOverHt_HT300_aT0.52_all");
    names.push_back("BabyOverHt_HT350_aT0.52_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0_all");
    names.push_back("BabyOverMeff_Meff400_aT0_all");
    names.push_back("BabyOverMeff_Meff450_aT0_all");
    names.push_back("BabyOverMeff_Meff500_aT0_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0.51_all");
    names.push_back("BabyOverMeff_Meff400_aT0.51_all");
    names.push_back("BabyOverMeff_Meff450_aT0.51_all");
    names.push_back("BabyOverMeff_Meff500_aT0.51_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0.52_all");
    names.push_back("BabyOverMeff_Meff400_aT0.52_all");
    names.push_back("BabyOverMeff_Meff450_aT0.52_all");
    names.push_back("BabyOverMeff_Meff500_aT0.52_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_Ht0_aT0_all");
    names.push_back("GenBabyOverHt_HT250_aT0_all");
    names.push_back("GenBabyOverHt_HT300_aT0_all");
    names.push_back("GenBabyOverHt_HT350_aT0_all");
    createPlot( path, "GenBabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverMeff_Meff0_aT0_all");
    names.push_back("GenBabyOverMeff_Meff400_aT0_all");
    names.push_back("GenBabyOverMeff_Meff450_aT0_all");
    names.push_back("GenBabyOverMeff_Meff500_aT0_all");
    createPlot( path, "GenBabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  output->Write();
  output->Close();
  delete output; 

  return 0;

}


