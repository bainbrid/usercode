#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>

// -----------------------------------------------------------------------------
//
class CompareEnergy {
public:
  CompareEnergy( std::string&,
		 std::string&,
		 std::string& );
  ~CompareEnergy() {;}
  void setTDRStyle();
};

// -----------------------------------------------------------------------------
//
CompareEnergy::CompareEnergy( std::string& vectorial,
			      std::string& tracks ) {
  
  setTDRStyle();

  const Int_t nbins = 10;
  const Int_t nbx = nbins + 1;
  const Float_t xbins[nbx] = { 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 150. };

  // Scale

  {
  
    TCanvas* scale = new TCanvas("ScaleX015","ScaleY015",1);
  
    TFile* vectorial_file = TFile::Open( vectorial.c_str() );
    TH1F* vectorial_histo = 0;
    if ( vectorial_file ) { vectorial_histo = (TH1F*)vectorial_file->Get("hScaleJPT"); }
    else { std::cout << "file is null" << std::endl; }
    if ( vectorial_histo ) { 
      vectorial_histo->GetXaxis()->SetTitle("E_{T}^{gen} [GeV]");
      vectorial_histo->GetYaxis()->SetTitle("E_{T}^{reco}/E_{T}^{gen}");
      vectorial_histo->GetXaxis()->SetTitleOffset(1.2);
      vectorial_histo->GetYaxis()->SetTitleOffset(1.9);
      vectorial_histo->SetMaximum(1.01);
      vectorial_histo->SetMinimum(0.95);
      vectorial_histo->SetMarkerStyle(24);
      vectorial_histo->SetMarkerSize(1.0);
      vectorial_histo->Draw("histPE1");
    }
    else { std::cout << "histo is null" << std::endl; }

    TFile* tracks_file = TFile::Open( tracks.c_str() );
    TH1F* tracks_histo = 0;
    if ( tracks_file ) { tracks_histo = (TH1F*)tracks_file->Get("hScaleJPT"); }
    else { std::cout << "file is null" << std::endl; }
    if ( tracks_histo ) { 
      tracks_histo->SetMarkerStyle(20);
      tracks_histo->SetMarkerSize(1.0);
      tracks_histo->Draw("samePE1"); 
    }
    else { std::cout << "histo is null" << std::endl; }

    TLatex* scale_latex = new TLatex();
    if ( scale_latex ) { 
      scale_latex->SetTextSize(0.03);
      scale_latex->DrawLatex(35,1.004,"CMSSW_3_4_0, RelVal QCD 80-120 GeV, |#eta^{jet}| < 1.0");
    }
  
    TLegend* scale_legend = new TLegend(0.5,0.6,0.9,0.7,NULL,"brNDC");
    if ( scale_legend ) {
      scale_legend->SetFillColor(0);
      scale_legend->AddEntry(vectorial_histo,"Vectorial, tracks only (#DeltaR = 0.30)","P");
      scale_legend->AddEntry(tracks_histo,"Vectorial, tracks only (#DeltaR = 0.15)","P");
      scale_legend->Draw();  
    }

    scale->SaveAs("EnergyScaleComparison-0.15.png");

  }

  // Resolution

  {
  
    TCanvas* resolution = new TCanvas("ResX015","ResY015",1);
  
    TFile* vectorial_file = TFile::Open( vectorial.c_str() );
    TH1F* vectorial_histo = 0;
    if ( vectorial_file ) { vectorial_histo = (TH1F*)vectorial_file->Get("hResJPT"); }
    else { std::cout << "file is null" << std::endl; }
    if ( vectorial_histo ) { 
      vectorial_histo->GetXaxis()->SetTitle("E_{T}^{gen} [GeV]");
      vectorial_histo->GetYaxis()->SetTitle("Energy resolution [%]");
      vectorial_histo->GetXaxis()->SetTitleOffset(1.2);
      vectorial_histo->GetYaxis()->SetTitleOffset(1.9);
      vectorial_histo->SetMaximum(0.20);
      vectorial_histo->SetMinimum(0.05);
      vectorial_histo->SetMarkerStyle(24);
      vectorial_histo->SetMarkerSize(1.0);
      vectorial_histo->Draw("histPE1"); 
    }
    else { std::cout << "histo is null" << std::endl; }

    TFile* tracks_file = TFile::Open( tracks.c_str() );
    TH1F* tracks_histo = 0;
    if ( tracks_file ) { tracks_histo = (TH1F*)tracks_file->Get("hResJPT"); }
    else { std::cout << "file is null" << std::endl; }
    if ( tracks_histo ) { 
      tracks_histo->SetMarkerStyle(20);
      tracks_histo->SetMarkerSize(1.0);
      tracks_histo->Draw("samePE1"); 
    }
    else { std::cout << "histo is null" << std::endl; }

    TLatex* resolution_latex = new TLatex();
    if ( resolution_latex ) { 
      resolution_latex->SetTextSize(0.03);
      resolution_latex->DrawLatex(35,0.185,"CMSSW_3_4_0, RelVal QCD 80-120 GeV, |#eta^{jet}| < 1.0");
    }
  
    TLegend* resolution_legend = new TLegend(0.5,0.6,0.9,0.7,NULL,"brNDC");
    if ( resolution_legend ) {
      resolution_legend->SetFillColor(0);
      resolution_legend->AddEntry(vectorial_histo,"Vectorial, tracks only (#DeltaR = 0.30)","P");
      resolution_legend->AddEntry(tracks_histo,"Vectorial, tracks only (#DeltaR = 0.15)","P");
      resolution_legend->Draw();  
    }

    resolution->SaveAs("EnergyResolutionComparison-0.15.png");

  }

}

// -----------------------------------------------------------------------------
//
void CompareEnergy::setTDRStyle() {
  
  TStyle* tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); // Height of canvas
  tdrStyle->SetCanvasDefW(600); // Width of canvas
  tdrStyle->SetCanvasDefX(0);   // Position on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
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
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  //tdrStyle->SetTitleOffset(1.2, "y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_

  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);

  tdrStyle->cd();

}
