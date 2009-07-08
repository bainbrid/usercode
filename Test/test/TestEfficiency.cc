#include "TestEfficiency.h"
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cstdio>
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

// -----------------------------------------------------------------------------
//
TestEfficiency::TestEfficiency( const edm::ParameterSet& pset ) 
  : append_( pset.getUntrackedParameter<std::string>("AppendString") ),
    file_( pset.getUntrackedParameter<std::string>("OutputFile") ),
    signal_( pset.getUntrackedParameter< std::vector<std::string> >("SignalFiles") ),
    bkgd_( pset.getUntrackedParameter< std::vector<std::string> >("BkgdFiles") ),
    histos_( pset.getUntrackedParameter< std::vector<std::string> >("Histograms") )
{
}

// -----------------------------------------------------------------------------
//
void TestEfficiency::beginJob( const edm::EventSetup& ) {

  TFile output_file( file_.c_str(), "RECREATE" );

  efficiencyPlots( output_file, signal_ );
  efficiencyPlots( output_file, bkgd_ );

  efficiencyCurves( output_file );
  significanceCurves( output_file );

  output_file.Close();

}

// -----------------------------------------------------------------------------
//
void TestEfficiency::efficiencyPlots( TFile& output_file,
				      const vstring& input_files ) {
  
  output_file.cd();
  TDirectory* output_dir = (TDirectory*)gDirectory;
  
  // Iterate through files
  vstring::const_iterator ifile = input_files.begin();
  vstring::const_iterator jfile = input_files.end();
  for ( ; ifile != jfile; ++ifile ) {

    // Open file
    TFile input_file( ifile->c_str() );
    input_file.cd();
    TDirectory* input_dir = (TDirectory*)gDirectory;

    // Create sub-directory and cd into it
    output_dir->mkdir( ifile->c_str() )->cd();
    
    // Iterate through different histograms
    vstring::const_iterator ihisto = histos_.begin();
    vstring::const_iterator jhisto = histos_.end();
    for ( ; ihisto != jhisto; ++ihisto ) {
      
      TH1F* orig = (TH1F*)input_dir->Get( ihisto->c_str() );
      if ( orig ) { 
	
	// Create cumulative histogram
	TH1F* histo_cumulative = (TH1F*)orig->Clone( std::string( std::string( orig->GetName() ) + "_Cumulative" ).c_str() );
	if ( histo_cumulative ) {
	  for ( int ii = 1; ii < histo_cumulative->GetNbinsX()+1; ii++ ) {
	    float total = 0;
	    for ( int jj = ii; jj < histo_cumulative->GetNbinsX()+1; jj++ ) { 
	      total += orig->GetBinContent(jj);
	    }
	    histo_cumulative->SetBinContent(ii,total); 
	  }
	  histo_cumulative->SetEntries( orig->GetEntries() ); 
	} else { cout << "No clone histo!" << endl; }

	// Create scaling histogram
	TH1F* histo_scale = (TH1F*)orig->Clone( std::string( std::string( orig->GetName() ) + "_Scale" ).c_str() );
	if ( histo_scale ) {
	  for ( int ii = 1; ii < histo_scale->GetNbinsX()+1; ii++ ) {
	    float total = 0;
	    for ( int jj = 1; jj < histo_scale->GetNbinsX()+1; jj++ ) { 
	      total += orig->GetBinContent(jj);
	    }
	    histo_scale->SetBinContent(ii,total); 
	  }
	  histo_scale->SetEntries( orig->GetEntries() ); 
	} else { cout << "No clone histo!" << endl; }
	
	// Create efficiency histogram
	TH1F* histo_eff = (TH1F*)orig->Clone( std::string( std::string( orig->GetName() ) + "_Efficiency" ).c_str() );
	if ( histo_eff && histo_cumulative && histo_scale ) { 
	  histo_eff->Divide( histo_cumulative, histo_scale, 1., 1., "B" ); 
	}
	
	orig->Write();
	histo_cumulative->Write();
	histo_scale->Write();
	histo_eff->Write();

      } else { cout << "No original histo!" << endl; }
      
    } // loop over histos
    
    input_file.Close();

  } // loop over files
  
}

// -----------------------------------------------------------------------------
//
void TestEfficiency::efficiencyCurves( TFile& output_file ) {

  output_file.cd();
  TDirectory* output_dir = (TDirectory*)gDirectory;

  for ( int iplot = 0; iplot < 2; ++iplot ) {
  
    // Iterate through signal files
    vstring::const_iterator isignal = signal_.begin();
    vstring::const_iterator jsignal = signal_.end();
    for ( ; isignal != jsignal; ++isignal ) {
    
      // Iterate through bkgd files
      vstring::const_iterator ibkgd = bkgd_.begin();
      vstring::const_iterator jbkgd = bkgd_.end();
      for ( ; ibkgd != jbkgd; ++ibkgd ) {

	// Cntr
	int cntr = 0;
      
	// Create canvas
	TCanvas canvas;
	canvas.Update();
	canvas.Clear();
      
	TLegend* legend = new TLegend(0.6,0.2,0.8,0.4);

	// Iterate through histograms
	vstring::const_iterator ihisto = histos_.begin();
	vstring::const_iterator jhisto = histos_.end();
	for ( ; ihisto != jhisto; ++ihisto ) {

	  std::stringstream signal_ss;
	  signal_ss << *isignal << ihisto->substr( ihisto->find_last_of("/"), std::string::npos ) << "_Efficiency";
	  TH1F* signal_histo = (TH1F*)output_dir->Get( signal_ss.str().c_str() );

	  std::stringstream bkgd_ss;
	  bkgd_ss << *ibkgd << ihisto->substr( ihisto->find_last_of("/"), std::string::npos ) << "_Efficiency";
	  TH1F* bkgd_histo = (TH1F*)output_dir->Get( bkgd_ss.str().c_str() );

	  vector<float> signal_value;
	  vector<float> signal_error;
	  if ( signal_histo ) {
	    for ( int ii = 0; ii < signal_histo->GetNbinsX(); ii++ ) {
	      float value = signal_histo->GetBinContent(ii+1) > 0. ? signal_histo->GetBinContent(ii+1) : 1.e-9;
	      signal_value.push_back( value );
	      signal_error.push_back( signal_histo->GetBinError(ii+1) );
	    }
	  } else { cout << "No signal histo!" << endl; }
	  
	  vector<float> bkgd_value;
	  vector<float> bkgd_error;
	  if ( bkgd_histo ) {
	    for ( int ii = 0; ii < bkgd_histo->GetNbinsX(); ii++ ) {
	      float value = bkgd_histo->GetBinContent(ii+1) > 0. ? bkgd_histo->GetBinContent(ii+1) : 1.e-9;
	      bkgd_value.push_back( value );
	      bkgd_error.push_back( bkgd_histo->GetBinError(ii+1) );
	    }
	  } else { cout << "No bkgd histo!" << endl; }
	  
	  if ( !signal_value.empty() ) {
	    
	    TGraphErrors* graph = new TGraphErrors( signal_value.size(), 
						    &bkgd_value.front(), &signal_value.front(), 
						    &bkgd_error.front(), &signal_error.front() );
	    
	    graph->GetXaxis()->SetTitle("Eff_{bkgd}");
	    graph->GetYaxis()->SetTitle("Eff_{signal}");
	    if ( iplot == 0 ) {
	      graph->GetXaxis()->SetRangeUser(0.,1.);
	      graph->GetYaxis()->SetRangeUser(0.,1.); 
	    } else {
	      canvas.cd();
	      canvas.SetLogx();
	      canvas.SetLogy();
	      graph->GetXaxis()->SetLimits(1.e-10,0.1);
	      graph->GetYaxis()->SetRangeUser(1.e-10,0.1); 
	    }
	    graph->SetMarkerStyle(24+cntr%10);
	    graph->SetMarkerColor(2+cntr%7);

	    std::string his = ihisto->substr( ihisto->find_last_of("/")+1, std::string::npos );
	    legend->AddEntry( graph, his.c_str() ,"LPE" );

	    canvas.cd();
	    if ( !cntr ) { graph->Draw("ALP"); }
	    else { graph->Draw("LP"); }
	    legend->Draw();
	    
	    output_dir->cd();
	    graph->Write();

	    
	  }

	  cntr++;
	
	} // loop through histograms
      
	// Save canvas
	stringstream ss;
	ss << "Eff_" << append_ << "_" << *isignal << "_" << *ibkgd << "_" << iplot << ".png";
	canvas.SaveAs( ss.str().c_str() );
	
      } // loop through bkgd samples
    
    } // loop through signal samples
  
  }
  
}

// -----------------------------------------------------------------------------
//
void TestEfficiency::significanceCurves( TFile& output_file ) {

  output_file.cd();
  TDirectory* output_dir = (TDirectory*)gDirectory;

  // Iterate through signal files
  vstring::const_iterator isignal = signal_.begin();
  vstring::const_iterator jsignal = signal_.end();
  for ( ; isignal != jsignal; ++isignal ) {
    
    // Iterate through bkgd files
    vstring::const_iterator ibkgd = bkgd_.begin();
    vstring::const_iterator jbkgd = bkgd_.end();
    for ( ; ibkgd != jbkgd; ++ibkgd ) {

      // Cntr
      int cntr = 0;
      
      // Create canvas
      TCanvas canvas;
      canvas.Update();
      canvas.Clear();
      
      TLegend* legend = new TLegend(0.2,0.2,0.4,0.4);

      // Iterate through histograms
      vstring::const_iterator ihisto = histos_.begin();
      vstring::const_iterator jhisto = histos_.end();
      for ( ; ihisto != jhisto; ++ihisto ) {

	std::stringstream signal_ss;
	signal_ss << *isignal << ihisto->substr( ihisto->find_last_of("/"), std::string::npos ) << "_Cumulative";
	TH1F* signal_histo = (TH1F*)output_dir->Get( signal_ss.str().c_str() );
	  
	std::stringstream bkgd_ss;
	bkgd_ss << *ibkgd << ihisto->substr( ihisto->find_last_of("/"), std::string::npos ) << "_Cumulative";
	TH1F* bkgd_histo = (TH1F*)output_dir->Get( bkgd_ss.str().c_str() );

	vector<double> signal_value;
	vector<double> signal_error;
	if ( signal_histo ) {
	  for ( int ii = 0; ii < signal_histo->GetNbinsX(); ii++ ) {
	    signal_value.push_back( signal_histo->GetBinContent(ii+1) );
	    signal_error.push_back( signal_histo->GetBinError(ii+1) );
	  }
	} else { cout << "No signal histo!" << endl; }
	  
	vector<double> bkgd_value;
	vector<double> bkgd_error;
	if ( bkgd_histo ) {
	  for ( int ii = 0; ii < bkgd_histo->GetNbinsX(); ii++ ) {
	    bkgd_value.push_back( bkgd_histo->GetBinContent(ii+1) );
	    bkgd_error.push_back( bkgd_histo->GetBinError(ii+1) );
	  }
	} else { cout << "No bkgd histo!" << endl; }

	vector<double> bin_value;
	vector<double> bin_error;
	if ( signal_histo ) {
	  for ( int ii = 0; ii < signal_histo->GetNbinsX(); ii++ ) {
	    bin_value.push_back( signal_histo->GetXaxis()->GetBinCenter(ii+1) );
	    bin_error.push_back( 0. ); //signal_histo->GetXaxis()->GetBinWidth(ii+1)/2. );
	  }
	} else { cout << "No signal histo!" << endl; }

	vector<double> ratio_value;
	vector<double> ratio_error;
	for ( uint32_t ii = 0; ii < signal_value.size(); ii++ ) {
	  double ratio = 1.e-6;
	  double error = 0.;
	  if ( signal_value[ii] > 0. ) { 
	    if ( bkgd_value[ii] > 0. ) { 
	      ratio = signal_value[ii] > 0. ? signal_value[ii] / sqrt( bkgd_value[ii] ) : 0.; 
	      double tmp1 = signal_value[ii] > 0. ? signal_error[ii] / signal_value[ii] : 0.;
	      double tmp2 = bkgd_value[ii] > 0. ? bkgd_error[ii] / bkgd_value[ii] : 0.;
	      error = ratio * sqrt( tmp1*tmp1 + tmp2*tmp2 );
	    } else {
	      ratio = 10. + 1.5 * cntr;
	      error = 0.;
	    }
	  }
	  ratio_value.push_back( ratio );
	  ratio_error.push_back( error );
	}
	
	TH1F* histo = (TH1F*)output_dir->Get( signal_ss.str().c_str() );
	if ( histo ) {
	  std::string name = signal_ss.str() + "_" + bkgd_ss.str() + "_S/B";
	  TH1F* ratio_histo = new TH1F( name.c_str(), "", 
					histo->GetXaxis()->GetNbins(), 
					histo->GetXaxis()->GetXmin(), 
					histo->GetXaxis()->GetXmax() );
	  if ( ratio_histo ) {
	    float entries = histo->GetEntries();
	    for ( uint16_t ii = 0; ii < ratio_value.size(); ii++ ) {
	      ratio_histo->SetBinContent(ii+1,ratio_value[ii]);
	      ratio_histo->SetBinError(ii+1,ratio_error[ii]);
	    }
	    ratio_histo->SetEntries( entries );
	      
	    TGraph* graph = new TGraph( signal_value.size(), 
					&bin_value.front(), &ratio_value.front() );
	    
	    graph->GetXaxis()->SetTitle("#alpha_{T}");
	    graph->GetYaxis()->SetTitle("S / #surd B");
	    graph->GetXaxis()->SetLimits(-0.5,2.5);
	    graph->GetYaxis()->SetRangeUser(1.e-7,20.); 
	    graph->SetMarkerStyle(24+cntr%10);
	    graph->SetMarkerColor(2+cntr%7);
	    
	    std::string his = ihisto->substr( ihisto->find_last_of("/")+1, std::string::npos );
	    legend->AddEntry( graph, his.c_str() ,"LPE" );
	    
	    canvas.cd();
	    if ( !cntr ) { graph->Draw("AP"); }
	    else { graph->Draw("P"); }
	    canvas.SetLogy();
	    legend->Draw();
	      
	    output_dir->cd();
	    ratio_histo->Write();
	    graph->Write();
	      
	  } else { cout << "No ratio histo!" << endl; }
	} else { cout << "No base histo!" << endl; }
	  
	cntr++;
	    
      } // loop through histograms
	
	// Save canvas
      stringstream ss;
      ss << "Sig_" << append_ << "_" << *isignal << "_" << *ibkgd << ".png";
      canvas.SaveAs( ss.str().c_str() );

    } // loop through bkgd samples
    
  } // loop through signal samples
  
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestEfficiency);











	//       // Some cosmetics
	//       histo_eff->SetLineColor(kBlack);
	//       histo_eff->SetMarkerSize(1);
	//       histo_eff->SetMarkerStyle(20);
	//       histo_eff->SetMarkerColor(kBlack);
	//       histo_eff->GetYaxis()->SetRangeUser(1.e-9,1);
	//       histo_eff->GetYaxis()->SetTitle("efficiency");
	//       histo_eff->GetXaxis()->SetTitle(xAxisTitle);
	//       histo_eff->SetTitle("");
	//       histo_eff->SetName(plotName);
	//       histo_eff->Draw("histep");
      
	
// 	// Calculate bin errors
// 	for( int bin = 1; bin < histo_eff->GetNbinsX()+1; bin++ ) {
// 	  float eff = histo_eff->GetBinContent(bin);
// 	  float tot = histo->GetBinContent(0);
// 	  float err = 0.;
// 	  if ( tot > 0. ) { err = sqrt( ( eff * ( 1. - eff ) ) / tot ); }
// 	  histo_eff->SetBinError(bin,err);
// 	  cout << " bin "
// 	       << bin << " eff " 
// 	       << eff << " 1-eff " 
// 	       << 1 - eff << " tot " 
// 	       << tot << " err "
// 	       << err << " his-err "
// 	       << histo_eff->GetBinError(bin) << " eff+err " 
// 	       << histo_eff->GetBinContent(bin) + err
// 	       << endl;
// 	}


