// -----------------------------------------------------------------------------
/*
*/  
int ratio2() {

  TCanvas* canvas = new TCanvas("Ratio");
  canvas->SetLogy();
  canvas->SetFillColor(0);
  gStyle->SetOptStat(0);
  
  TLegend* legend = new TLegend( 0.6, 0.6, 0.8, 0.8, NULL, "brNDC" );
  legend->SetFillColor(0);
  legend->SetLineColor(0); 
  
  const Int_t n = 6;
  float pt[n] = { 50., 30., 20., 50., 30., 20. };
  TString sample[n] = { "Gen50", "Gen30", "Gen20", "Reco50", "Reco30", "Reco20" };
  Int_t style[n]  = { 20, 21, 22, 24, 25, 26 };
  Int_t colour[n] = { 2, 3, 4, 2, 3, 4 };
//   float pt[n] = { 20., 20. };
//   TString sample[n] = { "Gen20", "Reco20" };
//   Int_t style[n]  = { 20, 24 };
//   Int_t colour[n] = { 2, 3 };

  const Int_t ngr = 1000;
  double x3[ngr];
  double r[ngr];
  int count = 0;

  for ( Int_t i = 0; i < n; ++i ) {
    TString name = "results/" + sample[i] + "_QCDPythia6.root";
    TFile* file = new TFile(name);
    file->cd();
    TH1* numerator = (TH1*)file->Get("Ratio50/HtPostAlphaT0.5_2");
    TH1* denominator = (TH1*)file->Get("Ratio50/HtPreAlphaT0.5_2");
    int rebin = 40;
    numerator->Rebin(rebin);
    denominator->Rebin(rebin);
    TH1* ratio = numerator->Clone();
    ratio->GetXaxis()->SetRangeUser(0.,1400.);
    ratio->Divide(denominator);
    ratio->SetMarkerStyle(style[i]);
    ratio->SetMarkerColor(colour[i]);

    if ( i < 3 ) {
      for ( Int_t ii = 1; ii < ratio->GetNbinsX(); ++ii ) {
	double val = ratio->GetBinContent(ii);
	if ( val == 0. ) { continue; }
	if ( count < ngr ) { 
	  double ht = ratio->GetBinLowEdge(ii);
	  if ( ht < 150. || ht > 650. ) { continue; }
	  double temp = ( 2. * pt[i] ) / ( ht + pt[i] );
	  std::cout << " ht: " << ht 
		    << " pt: " << pt[i]
		    << " x3: " << temp
		    << std::endl;
	  x3[count] = temp;
	  r[count] = val;
	  count++;
	}
      }
    }

    canvas->cd();
    if ( i == 0 ) ratio->Draw("");
    else ratio->Draw("same");
    legend->AddEntry( ratio, sample[i], "lep" );

    if (0) {

      int nbins = ratio->GetNbinsX();
      int bin_width = ratio->GetBinWidth(1);
      
      double lower = 0.;
      double upper = 1400.;
      
      int bin_lower = int( ( lower - ratio->GetBinLowEdge(1) ) / bin_width );
      for ( Int_t ii = bin_lower; ii < ratio->GetNbinsX()-1; ++ii ) {
	if ( ratio->GetBinContent(ii) > 0. ) { 
	  lower = ratio->GetBinCenter(ii);
	  break;
	}
      }
      int bin_upper = int( ( upper - ratio->GetBinLowEdge(1) ) / bin_width );
      for ( Int_t ii = bin_upper; ii > 0; --ii ) {
	if ( ratio->GetBinContent(ii) > 0. ) { 
	  upper = ratio->GetBinCenter(ii);
	  break;
	}
      }
      if (0) {
	std::cout << " bin_width: " << bin_width
		  << " bin_lower: " << bin_lower
		  << " bin_upper: " << bin_upper
		  << " lower: " << lower
		  << " upper: " << upper
		  << std::endl;
      }

      TF1* fit = new TF1(sample[i],"expo",lower,upper); 
      fit->SetLineColor(colour[i]);
      fit->SetLineWidth(1);
      ratio->Fit(sample[i],"QR","same");

    }

  }
  
  canvas->cd();
  legend->Draw("same");
  canvas->Update();

  TCanvas* c2 = new TCanvas("C2");
  c2->SetLogy();
  c2->SetFillColor(0);
  gStyle->SetOptStat(0);
  if ( count > 0 ) {
    TGraph* graph = new TGraph(count,x3,r); 
    graph->Draw("a*");
  }

  
}
