void photonId( void ) {

// Purity plot for SignalPhotons

TDirectory* signal_dir = (TDirectory*)_file0->Get("simplePhotonIDAnalysis/SignalPhotons");
signal_dir->cd();

TH1D* signal_matched_tight = (TH1D*)signal_dir->Get("MatchedTightPhotons")->Clone();
TH1D* signal_matched_loose = (TH1D*)signal_dir->Get("MatchedLoosePhotons")->Clone();
TH1D* signal_matched_em = (TH1D*)signal_dir->Get("MatchedLooseEMObjects")->Clone();
TH1D* signal_matched_all = (TH1D*)signal_dir->Get("MatchedEMObjects")->Clone();

TH1D* signal_total_tight = (TH1D*)signal_dir->Get("TotalTightPhotons")->Clone();
TH1D* signal_total_loose = (TH1D*)signal_dir->Get("TotalLoosePhotons")->Clone();
TH1D* signal_total_em = (TH1D*)signal_dir->Get("TotalLooseEMObjects")->Clone();
TH1D* signal_total_all = (TH1D*)signal_dir->Get("TotalEMObjects")->Clone();

TCanvas* signal_canvas = new TCanvas;
signal_canvas->cd();

TGraphAsymmErrors* graph_signal_tight = new TGraphAsymmErrors(signal_matched_tight,signal_total_tight);
graph_signal_tight->SetTitle("Purity for SignalPhotons");
graph_signal_tight->SetMarkerColor(2);
graph_signal_tight->SetMarkerStyle(20);
graph_signal_tight->GetYaxis()->SetRangeUser(0.,1.);
graph_signal_tight->Draw("ALP");

TGraphAsymmErrors* graph_signal_loose = new TGraphAsymmErrors(signal_matched_loose,signal_total_loose);
graph_signal_loose->SetMarkerColor(3);
graph_signal_loose->SetMarkerStyle(21);
graph_signal_loose->Draw("LP");

TGraphAsymmErrors* graph_signal_em = new TGraphAsymmErrors(signal_matched_em,signal_total_em);
graph_signal_em->SetMarkerColor(4);
graph_signal_em->SetMarkerStyle(22);
graph_signal_em->Draw("LP");

TGraphAsymmErrors* graph_signal_all = new TGraphAsymmErrors(signal_matched_all,signal_total_all);
graph_signal_all->SetMarkerColor(1);
graph_signal_all->SetMarkerStyle(25);
graph_signal_all->Draw("LP");

TLegend* signal_legend = new TLegend(0.6,0.2,0.8,0.4);
signal_legend->AddEntry(graph_signal_tight,"TightPhoton","LPE");
signal_legend->AddEntry(graph_signal_loose,"LoosePhoton","LPE");
signal_legend->AddEntry(graph_signal_em,"LooseEM","LPE");
signal_legend->AddEntry(graph_signal_all,"AllEM","LPE");
signal_legend->Draw();

// Purity plot for BkgdPhotons

TDirectory* bkgd_dir = (TDirectory*)_file0->Get("simplePhotonIDAnalysis/BkgdPhotons");
bkgd_dir->cd();

TH1D* bkgd_matched_tight = (TH1D*)bkgd_dir->Get("MatchedTightPhotons")->Clone();
TH1D* bkgd_matched_loose = (TH1D*)bkgd_dir->Get("MatchedLoosePhotons")->Clone();
TH1D* bkgd_matched_em = (TH1D*)bkgd_dir->Get("MatchedLooseEMObjects")->Clone();
TH1D* bkgd_matched_all = (TH1D*)bkgd_dir->Get("MatchedEMObjects")->Clone();

TH1D* bkgd_total_tight = (TH1D*)bkgd_dir->Get("TotalTightPhotons")->Clone();
TH1D* bkgd_total_loose = (TH1D*)bkgd_dir->Get("TotalLoosePhotons")->Clone();
TH1D* bkgd_total_em = (TH1D*)bkgd_dir->Get("TotalLooseEMObjects")->Clone();
TH1D* bkgd_total_all = (TH1D*)bkgd_dir->Get("TotalEMObjects")->Clone();

TCanvas* bkgd_canvas = new TCanvas;
bkgd_canvas->cd();

TGraphAsymmErrors* graph_bkgd_tight = new TGraphAsymmErrors(bkgd_matched_tight,bkgd_total_tight);
graph_bkgd_tight->SetTitle("Purity for BkgdPhotons");
graph_bkgd_tight->SetMarkerColor(2);
graph_bkgd_tight->SetMarkerStyle(20);
graph_bkgd_tight->GetYaxis()->SetRangeUser(0.,1.);
graph_bkgd_tight->Draw("ALP");

TGraphAsymmErrors* graph_bkgd_loose = new TGraphAsymmErrors(bkgd_matched_loose,bkgd_total_loose);
graph_bkgd_loose->SetMarkerColor(3);
graph_bkgd_loose->SetMarkerStyle(21);
graph_bkgd_loose->Draw("LP");

TGraphAsymmErrors* graph_bkgd_em = new TGraphAsymmErrors(bkgd_matched_em,bkgd_total_em);
graph_bkgd_em->SetMarkerColor(4);
graph_bkgd_em->SetMarkerStyle(22);
graph_bkgd_em->Draw("LP");

TGraphAsymmErrors* graph_bkgd_all = new TGraphAsymmErrors(bkgd_matched_all,bkgd_total_all);
graph_bkgd_all->SetMarkerColor(1);
graph_bkgd_all->SetMarkerStyle(25);
graph_bkgd_all->Draw("LP");

TLegend* bkgd_legend = new TLegend(0.6,0.2,0.8,0.4);
bkgd_legend->AddEntry(graph_bkgd_tight,"TightPhoton","LPE");
bkgd_legend->AddEntry(graph_bkgd_loose,"LoosePhoton","LPE");
bkgd_legend->AddEntry(graph_bkgd_em,"LooseEM","LPE");
bkgd_legend->AddEntry(graph_bkgd_all,"AllEM","LPE");
bkgd_legend->Draw();

// Purity plot for FakePhotons

TDirectory* fake_dir = (TDirectory*)_file0->Get("simplePhotonIDAnalysis/FakePhotons");
fake_dir->cd();

TH1D* fake_matched_tight = (TH1D*)fake_dir->Get("MatchedTightPhotons")->Clone();
TH1D* fake_matched_loose = (TH1D*)fake_dir->Get("MatchedLoosePhotons")->Clone();
TH1D* fake_matched_em = (TH1D*)fake_dir->Get("MatchedLooseEMObjects")->Clone();
TH1D* fake_matched_all = (TH1D*)fake_dir->Get("MatchedEMObjects")->Clone();

TH1D* fake_total_tight = (TH1D*)fake_dir->Get("TotalTightPhotons")->Clone();
TH1D* fake_total_loose = (TH1D*)fake_dir->Get("TotalLoosePhotons")->Clone();
TH1D* fake_total_em = (TH1D*)fake_dir->Get("TotalLooseEMObjects")->Clone();
TH1D* fake_total_all = (TH1D*)fake_dir->Get("TotalEMObjects")->Clone();

TCanvas* fake_canvas = new TCanvas;
fake_canvas->cd();

TGraphAsymmErrors* graph_fake_tight = new TGraphAsymmErrors(fake_matched_tight,fake_total_tight);
graph_fake_tight->SetTitle("Fake Rate for Photons");
graph_fake_tight->SetMarkerColor(2);
graph_fake_tight->SetMarkerStyle(20);
graph_fake_tight->GetYaxis()->SetRangeUser(0.,1.);
graph_fake_tight->Draw("ALP");

TGraphAsymmErrors* graph_fake_loose = new TGraphAsymmErrors(fake_matched_loose,fake_total_loose);
graph_fake_loose->SetMarkerColor(3);
graph_fake_loose->SetMarkerStyle(21);
graph_fake_loose->Draw("LP");

TGraphAsymmErrors* graph_fake_em = new TGraphAsymmErrors(fake_matched_em,fake_total_em);
graph_fake_em->SetMarkerColor(4);
graph_fake_em->SetMarkerStyle(22);
graph_fake_em->Draw("LP");

TGraphAsymmErrors* graph_fake_all = new TGraphAsymmErrors(fake_matched_all,fake_total_all);
graph_fake_all->SetMarkerColor(1);
graph_fake_all->SetMarkerStyle(25);
graph_fake_all->Draw("LP");

TLegend* fake_legend = new TLegend(0.6,0.6,0.8,0.8);
fake_legend->AddEntry(graph_fake_tight,"TightPhoton","LPE");
fake_legend->AddEntry(graph_fake_loose,"LoosePhoton","LPE");
fake_legend->AddEntry(graph_fake_em,"LooseEM","LPE");
fake_legend->AddEntry(graph_fake_all,"AllEM","LPE");
fake_legend->Draw();

// Efficiency for SignalPhotons

TDirectory* signal_dir = (TDirectory*)_file0->Get("simplePhotonIDAnalysis/SignalPhotons");
signal_dir->cd();

TH1D* signal_efficiency_tight = (TH1D*)signal_dir->Get("EfficiencyTightPhotons")->Clone();
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

// Efficiency for BkgdPhotons

TDirectory* bkgd_dir = (TDirectory*)_file0->Get("simplePhotonIDAnalysis/BkgdPhotons");
bkgd_dir->cd();

TH1D* bkgd_efficiency_tight = (TH1D*)bkgd_dir->Get("EfficiencyTightPhotons")->Clone();
TH1D* bkgd_efficiency_loose = (TH1D*)bkgd_dir->Get("EfficiencyLoosePhotons")->Clone();
TH1D* bkgd_efficiency_em = (TH1D*)bkgd_dir->Get("EfficiencyLooseEMObjects")->Clone();
TH1D* bkgd_efficiency_all = (TH1D*)bkgd_dir->Get("EfficiencyEMObjects")->Clone();

TH1D* gen_bkgd = (TH1D*)bkgd_dir->Get("GenPhotons")->Clone();

TCanvas* gen_bkgd_canvas = new TCanvas;
gen_bkgd_canvas->cd();

TGraphAsymmErrors* gen_graph_bkgd_tight = new TGraphAsymmErrors(bkgd_efficiency_tight,gen_bkgd);
gen_graph_bkgd_tight->SetTitle("Efficiency for BkgdPhotons");
gen_graph_bkgd_tight->SetMarkerColor(2);
gen_graph_bkgd_tight->SetMarkerStyle(20);
gen_graph_bkgd_tight->GetYaxis()->SetRangeUser(0.,1.);
gen_graph_bkgd_tight->Draw("ALP");

TGraphAsymmErrors* gen_graph_bkgd_loose = new TGraphAsymmErrors(bkgd_efficiency_loose,gen_bkgd);
gen_graph_bkgd_loose->SetMarkerColor(3);
gen_graph_bkgd_loose->SetMarkerStyle(21);
gen_graph_bkgd_loose->Draw("LP");

TGraphAsymmErrors* gen_graph_bkgd_em = new TGraphAsymmErrors(bkgd_efficiency_em,gen_bkgd);
gen_graph_bkgd_em->SetMarkerColor(4);
gen_graph_bkgd_em->SetMarkerStyle(22);
gen_graph_bkgd_em->Draw("LP");

TGraphAsymmErrors* gen_graph_bkgd_all = new TGraphAsymmErrors(bkgd_efficiency_all,gen_bkgd);
gen_graph_bkgd_all->SetMarkerColor(1);
gen_graph_bkgd_all->SetMarkerStyle(25);
gen_graph_bkgd_all->Draw("LP");

TLegend* gen_bkgd_legend = new TLegend(0.6,0.2,0.8,0.4);
gen_bkgd_legend->AddEntry(gen_graph_bkgd_tight,"TightPhoton","LPE");
gen_bkgd_legend->AddEntry(gen_graph_bkgd_loose,"LoosePhoton","LPE");
gen_bkgd_legend->AddEntry(gen_graph_bkgd_em,"LooseEM","LPE");
gen_bkgd_legend->AddEntry(gen_graph_bkgd_all,"AllEM","LPE");
gen_bkgd_legend->Draw();

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
