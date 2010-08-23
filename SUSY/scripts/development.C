#include <map>
#include <string>

// -----------------------------------------------------------------------------
/*
  Arguments to overlayPlots() method:
  Relative path to input files
  Canvas title
  Histogram title"
  Histogram directory 
  Rebin histograms
  Normalize all histograms
  Logarithmic y-scale
  Output file name
*/  
int overlay() {
  
  // SetSomeStyles();
  
  // Path to input files
  TString input_path = "../results/";
  
  // Jet algorithm to be used
  TString algo = "IC5Calo";
  
  // Define samples to be used
  std::map<TString,TString> samples;
  samples.push_back( std::pair<TString,TString>(TString("qcd"),TString("QCD_Pythia_Merged")) );
//   samples["qcd"]   = "QCD_Pythia_Merged";
//   samples["wjets"] = "WJets";
//   samples["tt"]    = "TTbarJets";
//   samples["zinv"]  = "Zinv";
//   samples["zjets"] = "ZJets";
//   samples["lm0"]   = "LM0";
//   samples["lm1"]   = "LM1";
  
  // Output file
  TFile* output_file = new TFile( "OverlayPlots.root", "RECREATE" );
  if ( !output_file || output_file->IsZombie() ) { return -1; }

  
  if (1) {
    overlayPlots( algo, samples, input_path, "Hadronic", "Eta_1", "LeadingJetEta", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "Eta_2", "SecondJetEta", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "DeltaEta_1", "DeltaEta", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "DeltaR_1", "DeltaR", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "Eta_1", "LeadingJetEtaNorm", 10, true, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "Eta_2", "SecondJetEtaNorm", 10, true, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "DeltaEta_1", "DeltaEtaNorm", 10, true, true, output_file );
    overlayPlots( algo, samples, input_path, "Hadronic", "DeltaR_1", "DeltaRNorm", 10, true, true, output_file );
  }
  
  if (1) {
    overlayPlots( algo, samples, input_path, "KinSuiteComPlot", "AlphaT_0", "AlphaT", 5, false, true, output_file );
    overlayPlots( algo, samples, input_path, "KinSuiteComPlot", "AlphaT_0", "AlphaTzoom", 1, false, true, output_file );
    overlayPlots( algo, samples, input_path, "CommonPlotsPost", "JetMultiplicity", "JetMultiplicity", 1, false, true, output_file );
    overlayPlots( algo, samples, input_path, "CommonPlotsPost", "SumEtall", "HT", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "CommonPlotsPost", "MHTall", "MHT", 10, false, true, output_file );
    overlayPlots( algo, samples, input_path, "CommonPlotsPost", "MPTMHTdPhi", "DeltaPhi_MPT_MHT", 2, false, true, output_file );
    overlayPlots( algo, samples, input_path, "CommonPlotsPost", "biasedJetMHTdPhiall", "BiasedDPhi", 2, false, true, output_file );
  }
  
  // Close output file
  output_file->Write();
  output_file->Close();
  delete output_file; 
  
}

// -----------------------------------------------------------------------------
//
bool overlayPlots( const TString& input_path,
		   const std::map<TString,TString>& samples,
		   const TString& input_path, 
		   const TString& histo_dir, 
		   const TString& histo_title, 
		   const TString& canvas_title, 
		   int rebin_histo, 
		   bool normalize_histo, 
		   bool log_scale,
		   TFile* output_file ) {
  
//   // Create legend
//   TLegend* legend = new TLegend( 0.75, 0.65, 0.92, 0.92, NULL, "brNDC" );
//   legend->SetFillColor(0);
//   legend->SetLineColor(0); 
  
//   // Create canvas
//   TCanvas* canvas = createCanvas( output_file, 
// 				  canvas_title, 
// 				  log_scale );
  
//   // Create histograms for different samples
//   std::map<TString,TH1D*> histos;
//   bool ok = retrieveHistos( input_path,
// 			    output_file,
// 			    histo_dir,
// 			    histo_title,
// 			    algo, 
// 			    samples, 
// 			    histos );
//   if ( !ok ) { return 0; }
  
//   // Combine Z+jets and Z->inv samples
//   bool combine = true;
//   TH1D* z_all = 0;
//   if ( combine ) { z_all = combineZ( histos ) }
  
//   // Set line attributes for different histograms
//   setLineAttributes( histos );

//   // Populate legend
//   legend->AddEntry( qcd, " QCD", "f" );
//   legend->AddEntry( lm0, " SUSY LM0", "L" );
//   legend->AddEntry( lm1, " SUSY LM1", "L" );
//   legend->AddEntry( tt_jets, " t#bar{t}+jets", "L" );
//   legend->AddEntry( w_jets, " W+jets", "L" );
//   if ( combine ) {
//     legend->AddEntry( z_all, " Z", "L" );
//   } else {
//     legend->AddEntry( z_jets, " Z+jets", "L" );
//     legend->AddEntry( z_inv, " Z#rightarrow#nu#nu", "L" );
//   }

//   // Calc maximum number of entries
//   double aMax = 0.;
//   if ( qcd->GetMaximum()     > aMax ) { aMax = qcd->GetMaximum(); }
//   if ( lm0->GetMaximum()     > aMax ) { aMax = lm0->GetMaximum(); }
//   if ( lm1->GetMaximum()     > aMax ) { aMax = lm1->GetMaximum(); }
//   if ( tt_jets->GetMaximum() > aMax ) { aMax = tt_jets->GetMaximum(); }  
//   if ( w_jets->GetMaximum()  > aMax ) { aMax = w_jets->GetMaximum(); }  
//   if ( combine ) {
//     if ( z_all->GetMaximum()  > aMax ) { aMax = z_all->GetMaximum(); }  
//   } else {
//     if ( z_inv->GetMaximum()   > aMax ) { aMax = z_inv->GetMaximum(); }  
//     if ( z_jets->GetMaximum()  > aMax ) { aMax = z_jets->GetMaximum(); }  
//   }

//   // Calc minimum number of entries
//   double aMin = 1.e12;
//   if ( qcd->GetMinimum(1.e-12)     < aMin ) { aMin = qcd->GetMinimum(1.e-12); }
//   if ( lm0->GetMinimum(1.e-12)     < aMin ) { aMin = lm0->GetMinimum(1.e-12); }
//   if ( lm1->GetMinimum(1.e-12)     < aMin ) { aMin = lm1->GetMinimum(1.e-12); }
//   if ( tt_jets->GetMinimum(1.e-12) < aMin ) { aMin = tt_jets->GetMinimum(1.e-12); }  
//   if ( w_jets->GetMinimum(1.e-12)  < aMin ) { aMin = w_jets->GetMinimum(1.e-12); }  
//   if ( combine ) {
//     if ( z_all->GetMinimum(1.e-12)   < aMin ) { aMin = z_all->GetMinimum(1.e-12); }  
//   } else {
//     if ( z_inv->GetMinimum(1.e-12)   < aMin ) { aMin = z_inv->GetMinimum(1.e-12); }  
//     if ( z_jets->GetMinimum(1.e-12)  < aMin ) { aMin = z_jets->GetMinimum(1.e-12); }  
//   }

//   if ( qcd ) qcd->GetYaxis()->SetTitleOffset(1.43);
//   if ( qcd ) qcd->GetYaxis()->SetTitleSize(0.06);
//   if ( qcd ) qcd->GetXaxis()->SetTitleSize(0.06);
//   if ( qcd ) qcd->GetXaxis()->SetTitleOffset(0.9);

//   if ( log_scale ) {
//     if ( qcd ) qcd->SetMaximum( aMax * 10. );
//     if ( qcd ) qcd->SetMinimum( aMin * 0.1 );
//   } else {
//     if ( qcd ) qcd->SetMaximum( aMax * 1.1 );
//     if ( qcd ) qcd->SetMinimum( aMin * 0.9 );
//   }

//   if ( normalize_histo ) {
//     if ( qcd ) qcd->DrawNormalized("Ehist");
//     if ( lm0->GetEntries() > 0. )     { lm0->DrawNormalized("hsame"); }
//     if ( lm1->GetEntries() > 0. )     { lm1->DrawNormalized("hsame"); }
//     if ( tt_jets->GetEntries() > 0. ) { tt_jets->DrawNormalized("hsame"); }
//     if ( w_jets->GetEntries() > 0. )  { w_jets->DrawNormalized("hsame"); }
//     if ( combine ) {
//       if ( z_all->GetEntries() > 0. )   { z_all->DrawNormalized("hsame"); }
//     } else {
//       if ( z_inv->GetEntries() > 0. )   { z_inv->DrawNormalized("hsame"); }
//       if ( z_jets->GetEntries() > 0. )  { z_jets->DrawNormalized("hsame"); }
//     }
//   } else {
//     if ( qcd ) qcd->Draw("h");
//     lm0->Draw("sameH");
//     lm1->Draw("sameH");
//     if ( tt_jets ) tt_jets->Draw("sameh");
//     w_jets->Draw("sameH");
//     if ( combine ) {
//       z_all->Draw("sameH");
//     } else {
//       z_inv->Draw("sameH");
//       z_jets->Draw("sameH");
//     }
//   }
  
//   output_file->cd();
//   legend->Draw("same");
//   canvas->Write();
//   return canvas;

}

// // -----------------------------------------------------------------------------
// //
// TCanvas* createCanvas( TFile* output_file, 
// 		       const TString& canvas_title,
// 		       bool log_scale ) {

//   output_file->cd();
//   TCanvas* canvas = new TCanvas(canvas_title);
//   gStyle->SetOptFit(1);
//   gStyle->SetOptStat(0);
//   canvas->Range(-288.2483,-2.138147,1344.235,6.918939);
//   canvas->SetFillColor(0);
//   canvas->SetBorderMode(0);
//   canvas->SetBorderSize(2);
//   canvas->SetLeftMargin(0.1765705);
//   canvas->SetRightMargin(0.05772496);
//   canvas->SetTopMargin(0.04778761);
//   canvas->SetBottomMargin(0.1256637);
//   canvas->SetFrameFillStyle(0);
//   canvas->SetFrameLineWidth(2);
//   canvas->SetFrameBorderMode(0);
//   canvas->SetFrameFillStyle(0);
//   canvas->SetFrameLineWidth(2);
//   canvas->SetFrameBorderMode(0);
  
//   if ( log_scale ) { canvas->SetLogy(); }
  
//   return canvas;
  
// }

// // -----------------------------------------------------------------------------
// //
// bool retrieveHistos( const TString& input_path, 
// 		     const TString& output_file, 
// 		     const TString& histo_dir, 
// 		     const TString& histo_title, 
// 		     const TString& algo, 
// 		     const std::map<TString,TString>& samples,
// 		     std::map<TString,TH1D*>& histos ) {
  
//   // Clear histos map
//   histos.clear();

//   // Loop through samples
//   std::map<TString,TString>::const_iterator i = samples.begin();
//   std::map<TString,TString>::const_iterator j = samples.end();
//   for ( ; i != j; ++i ) {

//     // Check if histogram already exists for given sample
//     TH1D* histo = histo(i->first,histos);
//     if ( histo ) { return false; }
    
//     // Extract filename stub for given sample
//     TString sample = sample(i->first,samples);
//     if ( sample.empty() ) { return false; }
    
//     // Construct complete filename for given algo and sample
//     TString input_file = a + "_" + sample + ".root";
    
//     // Create histogram for given algo and sample
//     TH1D* histo = retrieveHisto( input_path, 
// 				 input_file, 
// 				 histo_dir, 
// 				 histo_title, 
// 				 rebin_histo );
//     if ( histo ) { histos[i->first] = histo; }
    
//   }

//   return true;

// }

// // -----------------------------------------------------------------------------
// //
// TH1* retrieveHisto( TString input_path,
// 		    TString output_file,
// 		    TString histo_dir, 
// 		    TString histo_title,
// 		    int rebin_histo ) {
  
//   // Open output file and book histogram
//   TString name = input_path + output_file;
//   TFile* file =  new TFile( name, "RECREATE" );
//   if ( !file || file->IsZombie() ) { return 0; }
//   TDirectory* dir = (TDirectory*)file->Get( histo_dir );
//   if ( !dir ) { return 0; }
//   TH1* histo = (TH1*)dir->Get(histo_title);
//   if ( !histo ) {
//     std::cout << " Unable to book histogram with:"
// 	      << " Name: " << histo_title
// 	      << " Dir: " <<  histo_dir
// 	      << " File: " << output_file
// 	      << std::endl;
//     return 0;	
//   }

//   // Some histogram settings
//   histo->SetLineWidth(1);
//   histo->GetXaxis()->SetTitleSize(0.055);
//   histo->GetYaxis()->SetTitleSize(0.055);
//   histo->GetXaxis()->SetLabelSize(0.05);
//   histo->GetYaxis()->SetLabelSize(0.05);
//   histo->SetStats(kFALSE);
//   if ( rebin_histo > 0 ) { hist->Rebin(rebin_histo); }
  
//   return histo;

// }

// // -----------------------------------------------------------------------------
// //
// TString sample( const TString& s, 
// 		    const std::map<TString,TString>& m ) {
//   std::map<TString,TString>::const_iterator i = m.find(s);
//   if ( i != m.end() ) { return i->second; }
//   else { return ""; }
// }

// // -----------------------------------------------------------------------------
// //
// TH1D* histo( const TString& s, 
// 	     const std::map<TString,TH1D*>& m ) {
//   std::map<TString,TH1D*>::const_iterator i = m.find(s);
//   if ( i != m.end() && i->second ) { return i->second; }
//   else { return 0; }
// }

// // -----------------------------------------------------------------------------
// //
// TH1D* combineZ( std::map<TString,TH1D*>& m ) {
//   TH1D* z_all = 0;
//   TH1D* z_inv = histo("zinv");
//   TH1D* z_jets = histo("zjets");
//   if ( z_inv && z_jets ) { 
//     z_all = z_inv->Clone(); 
//     z_all->Add(z_jets,1);
//   } else if ( z_inv ) { 
//     z_all = z_inv->Clone(); 
//   } else if ( z_jets ) { 
//     z_all = z_jets->Clone(); 
//   }
//   if ( !histo("zall",m) && z_all ) { m["zall"] = z_all; }
//   return z_all;
// }

// // -----------------------------------------------------------------------------
// //
// bool setLineAttributes( const std::map<TString,TH1D*>& m ) {
  
//   TH1D* qcd = histo("qcd");
//   if ( qcd ) { 
//     qcd->SetLineColor(kGreen+2);
//     qcd->SetFillColor(kGreen+2);
//     qcd->SetFillStyle(3003);
//   }

//   TH1D* tt = histo("tt");
//   if ( tt ) { 
//     tt->SetLineColor(kBlue);
//     tt->SetLineStyle(1);
//     tt->SetLineWidth(1);
//   }

//   TH1D* w = histo("wjets");
//   if ( w ) { 
//     w_jets->SetLineColor(kBlue);
//     w_jets->SetLineStyle(3);
//     w_jets->SetLineWidth(1);
//   }

//   TH1D* z_all = histo("zall");
//   if ( z_all ) { 
//     z_all->SetLineColor(kBlack);
//     z_all->SetLineStyle(3);
//     z_all->SetLineWidth(1);
//   } else {
//     TH1D* z_inv = histo("zinv");
//     if ( z_inv ) { 
//       z_inv->SetLineColor(kBlack);
//       z_inv->SetLineStyle(1);
//       z_inv->SetLineWidth(1);
//     } 
//     TH1D* z_jets = histo("zjets");
//     if ( z_jets ) { 
//       z_jets->SetLineColor(kBlack);
//       z_jets->SetLineStyle(3);
//       z_jets->SetLineWidth(1);
//     }
//   }

//   TH1D* lm0 = histo("lm0");
//   if ( lm0 ) { 
//     lm0->SetLineColor(kRed);
//     lm0->SetLineStyle(1);
//     lm0->SetLineWidth(2);
//   }

//   TH1D* lm1 = histo("lm1");
//   if ( lm1 ) { 
//     lm1->SetLineColor(kRed);
//     lm1->SetLineStyle(3);
//     lm1->SetLineWidth(2);
//   }

//   return true;

// }
