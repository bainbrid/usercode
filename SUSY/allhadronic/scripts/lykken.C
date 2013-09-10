#include "common/style.C"
#include <TCanvas.h>
#include <TColor.h>
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
int lykken() {
  
  setTDRStyle();
  
  TFile* file =  new TFile("../../results/v30/Ratio__sm_pythia.root");
  if (!file) { return 0; }
  
  TDirectory* dir = (TDirectory*)file->Get("QcdBkgdEst");
  if (!dir) { return 0; }
  
  TH1* his = (TH1*)dir->Get("HtAfterTrackless_aT0_all");
  if (!his) { return 0; }

  TCanvas* c1 = new TCanvas("Lykken");
  c1->SetLogy();
  
  his->Scale(602./100.);
  his->Draw("h");
  
  
}

