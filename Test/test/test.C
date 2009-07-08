void test( void ) {

  TDirectory* signal_dir = (TDirectory*)_file0->Get("test/AlphaT");
  dir->cd();
  TH1D* histo = (TH1D*)dir->Get("AlphaT")->Clone();
  
  Double_t integral1 = integral( histo );

}


Double_t integral( TH1D* histo, Int_t xbin1, Int_t xbin2 ) {
  
  

  
}




void photonId( void ) {


// Efficiency for SignalPhotons

TH1D* signal_efficiency_loose = (TH1D*)signal_dir->Get("EfficiencyLoosePhotons")->Clone();
TH1D* signal_efficiency_em = (TH1D*)signal_dir->Get("EfficiencyLooseEMObjects")->Clone();
TH1D* signal_efficiency_all = (TH1D*)signal_dir->Get("EfficiencyEMObjects")->Clone();

TH1D* gen_signal = (TH1D*)signal_dir->Get("GenPhotons")->Clone();

TCanvas* gen_signal_canvas = new TCanvas;
gen_signal_canvas->cd();

TGraphAsymmErrors* gen_graph_signal_tight = new TGraphAsymmErrors(signal_efficiency_tight,gen_signal);
gen_graph_signal_tight->SetTitle("Efficiency for SignalPhotons");
gen_graph_signal_tight->SetMarkerColor(2);
gen_graph_signal_tight->SetMarkerStyle(20);
gen_graph_signal_tight->GetYaxis()->SetRangeUser(0.,1.);
gen_graph_signal_tight->Draw("ALP");

TGraphAsymmErrors* gen_graph_signal_loose = new TGraphAsymmErrors(signal_efficiency_loose,gen_signal);
gen_graph_signal_loose->SetMarkerColor(3);
gen_graph_signal_loose->SetMarkerStyle(21);
gen_graph_signal_loose->Draw("LP");

TGraphAsymmErrors* gen_graph_signal_em = new TGraphAsymmErrors(signal_efficiency_em,gen_signal);
gen_graph_signal_em->SetMarkerColor(4);
gen_graph_signal_em->SetMarkerStyle(22);
gen_graph_signal_em->Draw("LP");

TGraphAsymmErrors* gen_graph_signal_all = new TGraphAsymmErrors(signal_efficiency_all,gen_signal);
gen_graph_signal_all->SetMarkerColor(1);
gen_graph_signal_all->SetMarkerStyle(25);
gen_graph_signal_all->Draw("LP");

TLegend* gen_signal_legend = new TLegend(0.6,0.2,0.8,0.4);
gen_signal_legend->AddEntry(gen_graph_signal_tight,"TightPhoton","LPE");
gen_signal_legend->AddEntry(gen_graph_signal_loose,"LoosePhoton","LPE");
gen_signal_legend->AddEntry(gen_graph_signal_em,"LooseEM","LPE");
gen_signal_legend->AddEntry(gen_graph_signal_all,"AllEM","LPE");
gen_signal_legend->Draw();



}












// Double_t signal_multiplicity_truth_integral_double = 0;
// for ( Int_t i = 0; i < signal_multiplicity_truth->GetSize(); i++ ) { 
//   signal_multiplicity_truth_integral_double += 
//     signal_multiplicity_truth->GetBinCenter() * 
//     signal_multiplicity_truth->GetBinContent();
// }

// Double_t* signal_multiplicity_truth_integral_array = new Double_t[ signal_multiplicity_truth->GetSize() ];
// for ( Int_t i = 0; i < signal_multiplicity_truth->GetSize(); i++ ) { 
//   signal_multiplicity_truth_integral_array[i] = temp; 
// }

// Double_t* signal_matched_tight_integral_array = new Double_t[ signal_matched_tight->GetSize() ];
// for ( Int_t i = 0; i < signal_matched_tight->GetSize(); i++ ) { 
//   Double_t temp = 0;
//   for ( Int_t j = i; j < signal_matched_tight->GetSize(); j++ ) { 
//     temp += signal_matched_tight->GetBinContent();
//   }
//   signal_matched_tight_integral_array[i] = temp;
// }

// TH1D* signal_matched_tight = (TH1D*)signal_dir->Get("MatchedTightPhotons")->Clone();

// Set(Int_t n, const Double_t* array)

// signal_matched_tight 
