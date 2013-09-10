{
//=========Macro generated from canvas: pValueVersusAlphaT/
//=========  (Mon Jul  4 14:22:50 2011) by ROOT version5.27/06b
   TCanvas *pValueVersusAlphaT = new TCanvas("pValueVersusAlphaT", "",0,22,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   pValueVersusAlphaT->Range(0.482561,-0.1585366,0.6167073,1.060976);
   pValueVersusAlphaT->SetFillColor(0);
   pValueVersusAlphaT->SetBorderMode(0);
   pValueVersusAlphaT->SetBorderSize(2);
   pValueVersusAlphaT->SetTickx(1);
   pValueVersusAlphaT->SetTicky(1);
   pValueVersusAlphaT->SetLeftMargin(0.13);
   pValueVersusAlphaT->SetRightMargin(0.05);
   pValueVersusAlphaT->SetTopMargin(0.05);
   pValueVersusAlphaT->SetBottomMargin(0.13);
   pValueVersusAlphaT->SetFrameFillStyle(0);
   pValueVersusAlphaT->SetFrameBorderMode(0);
   pValueVersusAlphaT->SetFrameFillStyle(0);
   pValueVersusAlphaT->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(1);
   grae->SetName("Graph");
   grae->SetTitle("pValue");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.55,0.9273498);
   grae->SetPointError(0,0,0,0,0);
   
   TH1F *Graph10 = new TH1F("Graph10","pValue",100,0.45,1.65);
   Graph10->SetMinimum(0.8273498);
   Graph10->SetMaximum(2.02735);
   Graph10->SetDirectory(0);
   Graph10->SetStats(0);
   Graph10->SetFillColor(63);
   Graph10->SetLineStyle(0);
   Graph10->SetMarkerStyle(20);
   Graph10->GetXaxis()->SetNdivisions(1005);
   Graph10->GetXaxis()->SetLabelFont(42);
   Graph10->GetXaxis()->SetLabelOffset(0.007);
   Graph10->GetXaxis()->SetLabelSize(0.05);
   Graph10->GetXaxis()->SetTitleSize(0.06);
   Graph10->GetXaxis()->SetTitleOffset(0.9);
   Graph10->GetXaxis()->SetTitleFont(42);
   Graph10->GetYaxis()->SetLabelFont(42);
   Graph10->GetYaxis()->SetLabelOffset(0.007);
   Graph10->GetYaxis()->SetLabelSize(0.05);
   Graph10->GetYaxis()->SetTitleSize(0.06);
   Graph10->GetYaxis()->SetTitleOffset(1.05);
   Graph10->GetYaxis()->SetTitleFont(42);
   Graph10->GetZaxis()->SetLabelFont(42);
   Graph10->GetZaxis()->SetLabelOffset(0.007);
   Graph10->GetZaxis()->SetLabelSize(0.05);
   Graph10->GetZaxis()->SetTitleSize(0.06);
   Graph10->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph10);
   
   multigraph->Add(grae,"p");
   multigraph->Draw("al");
   multigraph->GetXaxis()->SetTitle("#alpha_{T} cut value");
   multigraph->GetXaxis()->SetRange(1,10);
   multigraph->GetXaxis()->SetNdivisions(1005);
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelOffset(0.007);
   multigraph->GetXaxis()->SetLabelSize(0.05);
   multigraph->GetXaxis()->SetTitleSize(0.06);
   multigraph->GetXaxis()->SetTitleOffset(0.9);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("p-value");
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelOffset(0.007);
   multigraph->GetYaxis()->SetLabelSize(0.05);
   multigraph->GetYaxis()->SetTitleSize(0.06);
   multigraph->GetYaxis()->SetTitleOffset(1.05);
   multigraph->GetYaxis()->SetTitleFont(42);
   pValueVersusAlphaT->Modified();
   pValueVersusAlphaT->cd();
   pValueVersusAlphaT->SetSelected(pValueVersusAlphaT);
}
