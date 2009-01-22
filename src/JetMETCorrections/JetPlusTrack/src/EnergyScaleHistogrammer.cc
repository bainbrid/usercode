#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistogrammer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistograms.h"
#include "TFile.h"
//#include "Rtypes.h"
#include "TStyle.h"

// -----------------------------------------------------------------------------
// 
EnergyScaleHistogrammer::EnergyScaleHistogrammer( const edm::ParameterSet& pset,
						  const edm::ActivityRegistry& reg ) 
  : histos_(),
    file_( new TFile( pset.getParameter<std::string>("RootFileName").c_str(), "RECREATE" ) ),
    style_( new TStyle( "TDRstyle" , "Style for P-TDR" ) )
{
  LogTrace("EnergyScale")
    << "[EnergyScaleHistogrammer::" << __func__ << "]"
    << " Constructing...";  
  if ( file_ ) { file_->cd(); }
  style();
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistogrammer::~EnergyScaleHistogrammer() {
  LogTrace("EnergyScale")
    << "[EnergyScaleHistogrammer::" << __func__ << "]"
    << " Destructing...";  
  std::vector<EnergyScaleHistograms*>::iterator ii = histos_.begin();
  std::vector<EnergyScaleHistograms*>::iterator jj = histos_.end();
  for ( ; ii != jj; ++ii ) { 
    if ( *ii ) { 
      (*ii)->write( file_ );
      (*ii)->end( style_ );
      delete *ii; 
    } 
  }
  histos_.clear();
  if ( file_ ) {
    file_->Close();
    file_->Write();
    delete file_;
  }
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistogrammer* EnergyScaleHistogrammer::Instance() {
  return edm::Service<EnergyScaleHistogrammer>().operator->(); 
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistograms* EnergyScaleHistogrammer::Histograms( const ObjectTags& tags ) {
  EnergyScaleHistogrammer* service = EnergyScaleHistogrammer::Instance();
  if ( service ) { return service->histograms( tags ); }
  else { return 0; }
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistograms* EnergyScaleHistogrammer::histograms( const ObjectTags& tags ) {
  std::vector<EnergyScaleHistograms*>::iterator ii = histos_.begin();
  std::vector<EnergyScaleHistograms*>::iterator jj = histos_.end();
  for ( ; ii != jj; ++ii ) { if ( *ii && (*ii)->tags() == tags ) { return *ii; } }
  EnergyScaleHistograms* histos = new EnergyScaleHistograms( tags );
  histos->book( file_ );
  histos_.push_back( histos );
  LogTrace("EnergyScale")
    << "[EnergyScaleHistogrammer::" << __func__ << "]"
    << " Added new EnergyScaleHistograms object!" 
    << " (size=" << histos_.size() << ")";
  return histos;
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleHistogrammer::style() {
  
  if ( !style_ ) {
    edm::LogError("EnergyScale")
      << "[EnergyScaleHistogrammer::" << __func__ << "]"
      << " NULL pointer to TStyle object!";
  }
  
  // Canvas
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); // Height of canvas
  style_->SetCanvasDefW(600); // Width of canvas
  style_->SetCanvasDefX(0);   // Position on screen
  style_->SetCanvasDefY(0);

  // Pad
  style_->SetPadBorderMode(0);
  //style_->SetPadBorderSize(Width_t size = 1);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // Frame
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);
  
  // Histo
  // style_->SetHistFillColor(1);
  // style_->SetHistFillStyle(0);
  style_->SetHistLineColor(1);
  style_->SetHistLineStyle(0);
  style_->SetHistLineWidth(2);
  // style_->SetLegoInnerR(Float_t rad = 0.5);
  // style_->SetNumberContours(Int_t number = 20);

  style_->SetEndErrorSize(4);
  // style_->SetErrorMarker(20);
  // style_->SetErrorX(0.);
  style_->SetMarkerStyle(20);

  // Fit/function
  style_->SetOptFit(1);
  style_->SetFitFormat("5.4g");
  style_->SetFuncColor(1);
  style_->SetFuncStyle(1);
  style_->SetFuncWidth(1);
  
  // Date
  style_->SetOptDate(0);
  // style_->SetDateX(Float_t x = 0.01);
  // style_->SetDateY(Float_t y = 0.01);

  // Stats box
  style_->SetOptFile(0);
  style_->SetOptStat(0); 
  // style_->SetOptStat("mr"); // display the mean and RMS
  style_->SetStatColor(kWhite);
  style_->SetStatFont(42);
  style_->SetStatFontSize(0.025);
  style_->SetStatTextColor(1);
  style_->SetStatFormat("6.4g");
  style_->SetStatBorderSize(1);
  style_->SetStatH(0.1);
  style_->SetStatW(0.15);
  // style_->SetStatStyle(Style_t style = 1001);
  // style_->SetStatX(Float_t x = 0);
  // style_->SetStatY(Float_t y = 0);

  // Margins
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.13);
  style_->SetPadLeftMargin(0.13);
  style_->SetPadRightMargin(0.05);

  // Global title
  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);
  // style_->SetTitleH(0); // Set the height of the title box
  // style_->SetTitleW(0); // Set the width of the title box
  // style_->SetTitleX(0); // Set the position of the title box
  // style_->SetTitleY(0.985); // Set the position of the title box
  // style_->SetTitleStyle(Style_t style = 1001);
  // style_->SetTitleBorderSize(2);

  // Axis titles
  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.06, "XYZ");
  // style_->SetTitleXSize(Float_t size = 0.02);
  // style_->SetTitleYSize(Float_t size = 0.02);
  style_->SetTitleXOffset(0.9);
  style_->SetTitleYOffset(1.05);
  // style_->SetTitleOffset(1.1, "Y");

  // Axis labels
  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.05, "XYZ");

  // Axis
  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1); // Tick marks on the opposite side of the frame
  style_->SetPadTickY(1);

  // Log plots
  style_->SetOptLogx(0);
  style_->SetOptLogy(0);
  style_->SetOptLogz(0);
  
  // Postscript options
  style_->SetPaperSize(15.,15.);
  // style_->SetPaperSize(20.,20.);
  // style_->SetPaperSize(7.5,7.5);
  
  // style_->SetLineScalePS(Float_t scale = 3);
  // style_->SetLineStyleString(Int_t i, const char* text);
  // style_->SetHeaderPS(const char* header);
  // style_->SetTitlePS(const char* pstitle);
  
  // style_->SetBarOffset(Float_t baroff = 0.5);
  // style_->SetBarWidth(Float_t barwidth = 0.5);
  // style_->SetPaintTextFormat(const char* format = "g");
  // style_->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // style_->SetTimeOffset(Double_t toffset);
  // style_->SetHistMinimumZero(kTRUE);
  
  style_->cd();
  
}
