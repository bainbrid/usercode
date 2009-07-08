#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TKey.h"
#include "TClass.h"
#include "TTree.h"
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TKey.h"
#include "TAxis.h"
#include "Riostream.h"

using namespace std;

// -----------------------------------------------------------------------------

void PlotEfficiency( char* inputfile, 
		     char* namehisto, 
		     char* plotName,
		     char* xAxisTitle, 
		     int rebin, 
		     char* outputfile)
{

  TFile infile(inputfile);
  infile.cd();

  TCanvas c1;
  c1.Update();
  c1.Clear();

  // Cumulative

  TH1F* histo1 = ((TH1F*)gDirectory->Get(namehisto)->Clone("Cumulative"));

  for ( int ii = 1; ii <= histo1->GetNbinsX(); ii++ ) {
    float total = 0;
    for ( int jj = ii; jj <= histo1->GetNbinsX(); jj++ ) {
      total += histo1->GetBinContent(jj);
    }
    histo1->SetBinContent(ii,total);
  }
  histo1->Rebin(rebin);

  // Integral
  
  TH1F *histo2 = ((TH1F*)gDirectory->Get(namehisto)->Clone("Integral"));
  
  float total = 0;
  for ( int ii = 1; ii <= histo2->GetNbinsX(); ii++ ) {
    total += histo2->GetBinContent(ii);
  }

  for ( int jj = 1; jj <= histo2->GetNbinsX(); jj++ ) {
    histo2->SetBinContent(jj,total);
  }
  histo2->Rebin(rebin);

  // Efficiency

  TH1F *histo_eff = ((TH1F*)histo1->Clone()); //histo_eff->Sumw2();

  histo_eff->Divide(histo1,histo2,1,1);
  histo_eff->SetStats(0); 

  for(int bin=1;bin<=histo_eff->GetNbinsX();bin++)
    {
      float eff = histo_eff->GetBinContent(bin);
      float one_minus_eff = 1- histo_eff->GetBinContent(bin);
      float Ntot =  float(histo2->GetBinContent(bin));

      if(Ntot>0)
	{
	  float err_eff = sqrt( (eff * one_minus_eff) / Ntot);
	  histo_eff->SetBinError(bin,err_eff);
	}
    }

  histo_eff->SetLineColor(kBlack);

  histo_eff->SetMarkerSize(1);
  histo_eff->SetMarkerStyle(20);
  histo_eff->SetMarkerColor(kBlack);

  histo_eff->GetYaxis()->SetRangeUser(0.,1);
  histo_eff->GetYaxis()->SetTitle("efficiency");
  histo_eff->GetXaxis()->SetTitle(xAxisTitle);

  histo_eff->SetTitle("");
  histo_eff->SetName(plotName);

  histo_eff->Draw("histep");

  float eff_tot=float(histo1->GetEntries())/float(histo2->GetEntries());
  float Ntot = float(histo2->GetEntries());
  float err_eff_tot = sqrt( (eff_tot * (1 - eff_tot) ) / Ntot);

//   cout << "******************************" << endl;
//   cout << plotName << endl;
//   cout << "eff +\- err_eff: " << eff_tot << " +/- " << err_eff_tot <<endl;
//   cout << "******************************" << endl;

  TFile outfile(outputfile,"UPDATE");

  outfile.cd();

  histo1->Write();
  histo2->Write();
  histo_eff->Write();

  outfile.Close();

  infile.Close();

}

// -----------------------------------------------------------------------------

void PlotEffGraph( char* inputfile, 
		   char* namehisto1, 
		   char* namehisto2, 
		   char* plotName,
		   char* xAxisTitle, 
		   int rebin, 
		   char* outputfile)
{
  
  TFile infile(inputfile);

  infile.cd();

  TCanvas c1;
  c1.Update();
  c1.Clear();
  c1.cd();

  // Get efficiency histograms
  TH1F* histo1 = ((TH1F*)gDirectory->Get(namehisto2)->Clone()); // bkgd
  TH1F* histo2 = ((TH1F*)gDirectory->Get(namehisto1)->Clone()); // signal
  histo1->Rebin(rebin);
  histo2->Rebin(rebin);

  // Create 2D histo
  TAxis* xaxis = histo1->GetXaxis();
  TAxis* yaxis = histo2->GetYaxis();

  vector<float> s_value;
  vector<float> s_error;
  for ( int ii = 1; ii <= xaxis->GetNbins(); ii++ ) {
    s_value.push_back( histo1->GetBinContent(ii) );
    s_error.push_back( histo1->GetBinError(ii) );
  }
  
  vector<float> b_value;
  vector<float> b_error;
  for ( int jj = 1; jj <= xaxis->GetNbins(); jj++ ) {
    b_value.push_back( histo2->GetBinContent(jj) );
    b_error.push_back( histo2->GetBinError(jj) );
  }

//   vector<float> r_value;
//   vector<float> r_error;
//   for ( int jj = 1; jj <= xaxis->GetNbins(); jj++ ) {
//     r_value.push_back( s_value );
//     r_error.push_back( 0. );
//   }

  // Create TGraph
  TGraphErrors* graph_eff = new TGraphErrors( s_value.size(), 
					      &s_value.front(), &b_value.front(), 
					      &s_error.front(), &b_error.front() );
  
  graph_eff->GetXaxis()->SetTitle("Efficiency(bkgd)");
  graph_eff->GetYaxis()->SetTitle("Efficiency(signal)");
  graph_eff->GetXaxis()->SetRangeUser(0.,1.);
  graph_eff->GetYaxis()->SetRangeUser(0.,1.);
  graph_eff->SetMarkerColor(2);
  graph_eff->SetMarkerStyle(20);
  graph_eff->Draw("ALP");
  stringstream name;
  name << plotName << ".jpg";
  c1.SaveAs(name.str().c_str());

  TFile outfile(outputfile,"UPDATE");
  
  outfile.cd();
  
  graph_eff->Write();
  
  outfile.Close();

  infile.Close();

}

// -----------------------------------------------------------------------------

void efficiency()
{

  gROOT->SetStyle("Plain");

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  
  gStyle->SetTitleBorderSize(0);
  
  gROOT->ForceStyle();

  char* outputfile="efficiency.root";
  
  TFile outfile(outputfile,"RECREATE");
  outfile.Close();

  cout << "********** START *********" << endl;

  // Alpha_T

  PlotEfficiency( "gmsb.root", "test/AlphaT/AlphaT", "GMSB_Eff_AlphaT", "\\alpha_t", 1, outputfile );
  PlotEfficiency( "bkgd.root", "test/AlphaT/AlphaT", "BKGD_Eff_AlphaT", "\\alpha_t", 1, outputfile );
  PlotEffGraph( "efficiency.root", "GMSB_Eff_AlphaT", "BKGD_Eff_AlphaT", "Eff_AlphaT", "\\alpha_t", 1, outputfile );

  // Beta_T

  PlotEfficiency( "gmsb.root", "test/AlphaT/BetaT", "GMSB_Eff_BetaT", "\\beta_t", 1, outputfile );
  PlotEfficiency( "bkgd.root", "test/AlphaT/BetaT", "BKGD_Eff_BetaT", "\\beta_t", 1, outputfile );
  PlotEffGraph( "efficiency.root", "GMSB_Eff_BetaT", "BKGD_Eff_BetaT", "Eff_BetaT", "\\beta_t", 1, outputfile );
  
  cout << "********** FINISHED *********" << endl;
  
}

// -----------------------------------------------------------------------------

int main()
{
  efficiency();
}


