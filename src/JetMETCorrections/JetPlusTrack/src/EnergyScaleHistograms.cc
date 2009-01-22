#include "JetMETCorrections/JetPlusTrack/interface/EnergyScaleHistograms.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <cmath>

// -----------------------------------------------------------------------------
// 
const float EnergyScaleHistograms::eBins_[nBins_+1] = { 
  20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,150. 
};

// -----------------------------------------------------------------------------
// 
EnergyScaleHistograms::EnergyScaleHistograms( const ObjectTags& tags ) 
  : tags_(tags),
    hScale_(),
    hEt_(0),
    hDR_(0),
    hScaleVsE_(0),
    hResVsE_(0)
{
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]"
    << " Constructing...";
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistograms::EnergyScaleHistograms() 
  : tags_(),
    hScale_(),
    hEt_(0),
    hDR_(0),
    hScaleVsE_(0),
    hResVsE_(0)
{
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]"
    << " Constructing...";
}

// -----------------------------------------------------------------------------
// 
EnergyScaleHistograms::~EnergyScaleHistograms() {
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]"
    << " Destructing...";  
  
  std::vector<TH1F*>::iterator ii = hScale_.begin();
  std::vector<TH1F*>::iterator jj = hScale_.end();
  for ( ; ii != jj; ++ii ) { if ( *ii ) { delete (*ii); } } 
  hScale_.clear();
  if ( hEt_ ) { delete hEt_; }
  if ( hDR_ ) { delete hDR_; }
  if ( hScaleVsE_ ) { delete hScaleVsE_; }
  if ( hResVsE_ ) { delete hResVsE_; }
  
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleHistograms::analyze( const std::vector<LorentzVectorPair>& input ) {
  
  // Iterate through LorentzVectorPair objects
  std::vector<LorentzVectorPair>::const_iterator ii = input.begin();
  std::vector<LorentzVectorPair>::const_iterator jj = input.end();
  for ( ; ii != jj; ++ii ) {
    
    if ( !ii->both() ) { continue; }
    
    uint16_t nhist = EnergyScaleHistograms::nBins_;
    for ( uint16_t ihist = 0; ihist < nhist; ++ihist ) {
      
      // Determine energy bin
      if( ii->gen().perp() >= eBins_[ihist] && 
	  ii->gen().perp() < eBins_[ihist+1]  ) { 
	
	float scale = 
	  ( ii->gen().perp() > 0. ) ? 
	  ( ii->reco().perp() / ii->gen().perp() ) :
	  ( -1. );
	
	// Check scale
	if ( scale < 0.1 ) { continue; }
	
	// Check separation
	if( ii->dR() > 0.3 ) { continue; }
	
	// Check eta
	if ( fabs( ii->gen().eta() ) < 0. ||
	     fabs( ii->gen().eta() ) > 1. ) { continue; }
	
	// Build pair of gen jets
	LorentzVectorPair other( ii->gen() );
	if ( input.size() < 2 ) { other.reco( 0., 0., 1., 1. ); }
	else {
	  if ( ii < input.end() - 1 ) { other.reco( (ii+1)->gen() ); }
	  else { other.reco( (ii-1)->gen() ); }
	}
	
	// Check dR for gen jets
	if ( other.dR() < 2. ) { continue; }
	
	// Fill histograms
	if ( ihist < hScale_.size() && hScale_[ihist] ) { hScale_[ihist]->Fill( scale ); }
	if ( hEt_ ) { hEt_->Fill( ii->gen().perp() ); }
	if ( hDR_ ) { hDR_->Fill( ii->dR() ); }
	
      }	
    }
  }
  
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleHistograms::book( TFile* const file ) {
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]";
  
  if ( file ) { file->cd(); }
  
  hScale_.resize( nBins_, 0 );
  for ( uint16_t ii = 0; ii < nBins_; ii++ ) { 
    std::stringstream scale_ss;
    scale_ss << "Scale_Bin" << ii << "_" << tags_.reco().str();
    hScale_[ii] = new TH1F( scale_ss.str().c_str(), 
			    scale_ss.str().c_str(), 
			    60, 0., 3. );
  }
  
  std::stringstream et_ss;
  et_ss << "Et_" << tags_.gen().str(); 
  hEt_ = new TH1F( et_ss.str().c_str(), 
		   et_ss.str().c_str(), 
		   20, 0., 200. );
  
  std::stringstream dr_ss;
  dr_ss << "DR_" << "_" << tags_.reco().str(); 
  hDR_ = new TH1F( dr_ss.str().c_str(), 
		   dr_ss.str().c_str(), 
		   100, 0., 10. );
  
  std::stringstream scalevse_ss;
  scalevse_ss << "ScaleVersusE_" << tags_.reco().str();
  hScaleVsE_ = new TH1F( scalevse_ss.str().c_str(), 
			 scalevse_ss.str().c_str(), 
			 nBins_, eBins_ );
  
  std::stringstream resvse_ss;
  resvse_ss << "ResVersusE_" << tags_.reco().str();
  hResVsE_   = new TH1F( resvse_ss.str().c_str(), 
			 resvse_ss.str().c_str(),
			 nBins_, eBins_ );
  
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleHistograms::write( TFile* const file ) {
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]";
  
  if ( file ) { file->cd(); }
  
  std::vector<TH1F*>::const_iterator ii = hScale_.begin();
  std::vector<TH1F*>::const_iterator jj = hScale_.end();
  for ( ; ii != jj; ++ii ) { if ( *ii ) { (*ii)->Write(); } } 
  if ( hEt_ )       { hEt_->Write(); } 
  if ( hDR_ )       { hDR_->Write(); }
  if ( hScaleVsE_ ) { hScaleVsE_->Write(); }
  if ( hResVsE_ )   { hResVsE_->Write(); }
  
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleHistograms::end( TStyle* const style ) {
  LogTrace("EnergyScale")
    << "[EnergyScaleHistograms::"<<__func__<<"]";
  
  // Set style
  if ( style ) { style->SetOptFit(); }
  
  // Create canvas 
  TCanvas* canvas = new TCanvas( "X","Y",1 );
  canvas->Divide(2,2);
  
  // Fit Et histograms to get scale and resolution vs E
  for ( uint32_t ihist = 0; ihist < nBins_; ++ihist ) {

    // Check histo pointers
    if ( !hScaleVsE_ || !hResVsE_ || !hScale_[ihist] ) { continue; }
    
    // Get center of appropriate bin for "scale vs E" histo 
    double scale_centre = hScaleVsE_->GetXaxis()->GetBinCenter(ihist+1);
    
    // Get center of appropriate bin for "resolution vs E" histo 
    double res_centre = hResVsE_->GetXaxis()->GetBinCenter(ihist+1);
    
    // Use canvas
    canvas->cd(1);
    
    // Find mid-value for bin with max content and histogram rms
    double  num = hScale_[ihist]->GetEntries();
    int32_t max = hScale_[ihist]->GetMaximumBin();
    double  cen = hScale_[ihist]->GetXaxis()->GetBinCenter(max);
    double  rms = hScale_[ihist]->GetRMS();
    
    // Fit histo
    hScale_[ihist]->Fit("gaus", "", "",
			cen - 2.0 * rms,
			cen + 2.0 * rms );
    TF1* fit = hScale_[ihist]->GetFunction("gaus"); 
    if ( style ) { style->SetOptFit(); }
    
    // Find scale and resolutions per bin
    double mean      = fit->GetParameter(1);
    double mean_err  = fit->GetParError(1);
    double sigma     = fit->GetParameter(2);
    double sigma_err = fit->GetParError(2);
    float  res       = 0.;
    float  res_err   = 0.;
    if ( mean < -1.e-9 || mean > 1.e-9 ) { 
      res = sigma / mean; 
      if ( sigma < -1.e-9 || sigma > 1.e-9 ) {
	res_err = res * sqrt( ( mean_err / mean ) * 
			      ( mean_err / mean ) + 
			      ( sigma_err / sigma ) *
			      ( sigma_err / sigma ) );
      }
    }
    
    std::stringstream ss;
    ss << "[EnergyScaleHistograms::" << __func__ << "]" 
       << " Scale and resolution for bin #" << ihist
       << " and energy " << eBins_[ihist] << std::endl
       << " Scale      : " << mean << " +/- " << mean_err << std::endl
       << " Resolution : " << res << " +/- " << res_err << std::endl
       << " Stats      : " << num << " +/- " << sqrt(num) << " (";
    float stat = 0;
    if ( num ) { stat = static_cast<uint32_t>( 1000. * sqrt(num) / num ) / 10.; }
    ss << stat << "%)";
    LogTrace("EnergyScale") << ss.str();
    
    // Fill histograms "vs E"
    hScaleVsE_->Fill( scale_centre, mean );
    hScaleVsE_->SetBinError( ihist+1, mean_err );    
    hResVsE_->Fill( res_centre, res );
    hResVsE_->SetBinError( ihist+1, res_err );    

  }
  
  // Scale histo
  {
    TCanvas* c20 = new TCanvas("X","Y",1);
  
    hScaleVsE_->GetXaxis()->SetTitle("E_{T} Gen [GeV]");
    hScaleVsE_->GetYaxis()->SetTitle("E_{T}^{reco}/E_{T}^{gen}");
  
    hScaleVsE_->SetMinimum(0.3);
    hScaleVsE_->SetMaximum(1.2);
    hScaleVsE_->SetMarkerStyle(21);
    hScaleVsE_->SetMarkerSize(1.2);
    hScaleVsE_->Draw("histPE1");
  
    TLegend* leg = new TLegend( 0.5, 0.15, 0.9, 0.35, NULL, "brNDC" );
    leg->SetFillColor(10);
    leg->AddEntry( hScaleVsE_, tags_.reco().str().c_str(),"P" );
    leg->Draw();  
  
    TLatex* t = new TLatex();
    t->SetTextSize( 0.042 );
    t->DrawLatex( 25, 1.12, "CMSSW220" );
    t->DrawLatex( 25, 1.06, "RelVal QCD 80-120 GeV, |#eta ^{jet}|< 1.0" );

    std::string f = "Scale_" + tags_.reco().str() + ".gif";
    c20->SaveAs( f.c_str() );
  }

  // Resolution histo
  {
    TCanvas* c40 = new TCanvas("X","Y",1);
    
    hResVsE_->GetXaxis()->SetTitle("E_{T} Gen [GeV]");
    hResVsE_->GetYaxis()->SetTitle("Energy resolution [%]");
    
    hResVsE_->SetMinimum(0.05);
    hResVsE_->SetMaximum(0.45);
    hResVsE_->SetMarkerStyle(21);
    hResVsE_->SetMarkerSize(1.2);
    hResVsE_->Draw("histPE1");

    TLegend* leg = new TLegend( 0.45, 0.5, 0.85, 0.8, NULL, "brNDC" );
    leg->SetFillColor(10);
    leg->AddEntry( hScaleVsE_, tags_.reco().str().c_str(),"P" );
    leg->Draw();  

    TLatex *t = new TLatex();
    t->SetTextSize(0.042);
    t->DrawLatex( 25, 0.42, "CMSSW220" );
    t->DrawLatex( 25, 0.40, "RelVal QCD 80-120 GeV, |#eta ^{jet}|< 1.0" );

    std::string f = "Resolution_" + tags_.reco().str() + ".gif";
    c40->SaveAs( f.c_str() );
  }
  
}
