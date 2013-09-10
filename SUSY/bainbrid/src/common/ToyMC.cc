#include "ToyMC.hh"
#include "Types.hh"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TTimeStamp.h"

#include <iomanip>
#include <sstream>
#include <string>

using namespace Operation;

// -----------------------------------------------------------------------------
//
ToyMC::ToyMC( const Utils::ParameterSet& ps ) {

  // Time the method
  time_t start = TTimeStamp().GetSec();

  output_ = new TFile( ps.Get<std::string>("OutputFile").c_str(),"recreate" );
  initDir( output_, ps.Get<std::string>("Directory").c_str() );

  // Depth and binning
  std::vector<Int_t> tnb( ps.Get< std::vector<Int_t> >("nBins") );
  const Int_t ndepth = tnb.size();
  Int_t nbins[ndepth];
  for ( Int_t ii = 0; ii < ndepth; ++ii ) { nbins[ii] = tnb[ii]; } 
  
  // Arrays for nested loops
  std::vector<Double_t> tat( ps.Get< std::vector<Double_t> >("nAT") );
  const Int_t nat = tat.size();
  Double_t at[nat];
  for ( Int_t ii = 0; ii < nat; ++ii ) { at[ii] = tat[ii]; } 

  std::vector<Double_t> tpt( ps.Get< std::vector<Double_t> >("nPT") );
  const Int_t npt = tpt.size();
  Double_t pt[npt];
  for ( Int_t ii = 0; ii < npt; ++ii ) { pt[ii] = tpt[ii]; } 

  std::vector<Double_t> tht( ps.Get< std::vector<Double_t> >("nHT") );
  const Int_t nht = tht.size();
  Double_t ht[nht];
  for ( Int_t ii = 0; ii < nht; ++ii ) { ht[ii] = tht[ii]; } 

  Double_t ratio[nat][npt][nht];
  Int_t length[nat][npt];

  Int_t loop = 0;
  Int_t nloops = nat * npt * nht;

  // Loop through bins of AlphaT, MinPt and HT 
  for ( Int_t iat = 0; iat < nat; ++iat ) {
    for ( Int_t ipt = 0; ipt < npt; ++ipt ) {
      length[iat][ipt] = 0;
      for ( Int_t iht = 0; iht < nht; ++iht ) {
	std::cout << "Completed " 
		  << 100.*float(loop)/float(nloops) 
		  << "%..." 
		  << std::endl; loop++;

	Double_t x3 = ( 2. * pt[ipt] ) / ( ht[iht] + pt[ipt] );
	
	Double_t numerator = 0.;
	Double_t denominator = 0.;
	
	Int_t depth = 0;
	integrate( depth, ndepth, nbins, x3, at[iat], numerator, denominator, 0., 1., 0., 1. );
	
	// Ratio
	if ( denominator > 0. ) { 
	  ratio[iat][ipt][iht] = numerator/denominator;
	  if ( ratio[iat][ipt][iht] > 0. ) { length[iat][ipt]++; }
	}

      }
    }
  }

  // Debug
  if (1) {
    for ( Int_t iat = 0; iat < nat; ++iat ) {
      std::cout << " AlphaT: " << at[iat] << std::endl;
      for ( Int_t ipt = 0; ipt < npt; ++ipt ) {
	std::cout << "  Pt: " << pt[ipt] << std::endl; 
	std::cout << "  Length: " << length[iat][ipt] << std::endl; 
	for ( Int_t iht = 0; iht < nht; ++iht ) {
	  Double_t x3 = ( 2. * pt[ipt] ) / ( ht[iht] + pt[ipt] );
	  std::cout << "   HT: " << ht[iht]
		    << " x3: " << x3
		    << " ratio: " << ratio[iat][ipt][iht]
		    << std::endl;
	}
      }
    }
  }

  if (1) {

    // Canvas for ratios from theory 
    TCanvas* c2 = new TCanvas( "c2", "" );
    c2->SetLogy();
    c2->cd();
    TMultiGraph* mg2 = new TMultiGraph();
    for ( Int_t iat = 0; iat < nat; iat++ ) {
      for ( Int_t ipt = 0; ipt < npt; ipt++ ) {
	if ( length[iat][ipt] == 0 ) { continue; }
	TGraph* gr = new TGraph(length[iat][ipt],ht,ratio[iat][ipt]);
	std::stringstream ss;
	ss << "AlphaT=" << at[iat] << ", Pt=" << pt[ipt];
	mg2->Add(gr,"lp");
	gr->SetTitle(TString(ss.str()));
	gr->SetLineColor(2+iat);
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(20+ipt);
	gr->SetMarkerColor(2+iat);
	gr->SetMarkerSize(1.5);
      }
    }
    mg2->Draw("a");
    //mg2->GetYaxis()->SetRangeUser(1.e-5,1.e-1);
    c2->Update();
    c2->BuildLegend(0.8,0.1,0.99,0.9);
  
    // Save canvases
    c2->cd();
    c2->SaveAs("c2.C");
    c2->SaveAs("c2.png");

    // Canvas for ratios from theory 
    TCanvas* c3 = new TCanvas( "c3", "" );
    c3->SetLogy();
    c3->cd();
    TMultiGraph* mg3 = new TMultiGraph();
    for ( Int_t iat = 0; iat < nat; iat++ ) {
      for ( Int_t ipt = 0; ipt < npt; ipt++ ) {
	if ( length[iat][ipt] == 0 ) { continue; }
	Double_t x3_[nht];
	for ( Int_t iht = 0; iht < nht; iht++ ) {
	  x3_[iht] = ( 2. * pt[ipt] ) / ( ht[iht] + pt[ipt] );
	}
	TGraph* gr = new TGraph(length[iat][ipt],x3_,ratio[iat][ipt]);
	std::stringstream ss;
	ss << "AlphaT=" << at[iat] << ", Pt=" << pt[ipt];
	mg3->Add(gr,"lp");
	gr->SetTitle(TString(ss.str()));
	gr->SetLineColor(2+iat);
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(20+ipt);
	gr->SetMarkerColor(2+iat);
	gr->SetMarkerSize(1.5);
      }
    }
    mg3->Draw("a");
    //mg3->GetYaxis()->SetRangeUser(1.e-5,1.e-1);
    c3->Update();
    c3->BuildLegend(0.8,0.1,0.99,0.9);

  }

  time_t stop = TTimeStamp().GetSec();
  std::cout << " Time taken: " << stop - start <<  " seconds" << std::endl;

  output_->Write();

}
// -----------------------------------------------------------------------------
//
ToyMC::~ToyMC() {
  if ( output_ ) {
    output_->cd();
    output_->Write();
    output_->Close();
    delete output_;
  }
}

// -----------------------------------------------------------------------------
/**
   Returns false if (x2,x1) fails any of the kinematic constraints
*/
bool ToyMC::constrain( Double_t x1, Double_t x2, Double_t x3 ) {
  return ( x1 >= x2 &&                   // jet ordering by Pt
	   ( x1 + x2 ) <= 2. &&          // from relation "x1 + x2 + x3 = 2"
	   x1 <= 1. &&                   // from "lost jet" and relation "xmiss = -x1 -x2"
	   ( x1 + x2 ) >= 1. &&          // from "lost jet" and relation "xmiss = -x1 -x2"
	   ( x1 + x2 ) >= ( 2. - x3 ) ); // mimimum value of x1 + x2 given x3
}

// -----------------------------------------------------------------------------
//
void ToyMC::integrate( Int_t depth, 
		       Int_t ndepth, 
		       Int_t* nbins, 
		       Double_t x3,
		       Double_t alpha_t,
		       Double_t& numerator,
		       Double_t& denominator, 
		       Double_t xlow, 
		       Double_t xhigh, 
		       Double_t ylow, 
		       Double_t yhigh  ) {
  
  bool debug1 = false;
  bool debug2 = false;
  
  // Count number of iterations
  depth++;

  Double_t exponent = 0.;
  for ( Int_t ii = 0; ii < ndepth; ++ii ) { exponent += log10(nbins[depth-1]); }
  if ( exponent > 12. ) {
    std::cout << " Too many bins required! " << std::endl;
    return;
  }
  
  // Bin widths
  Double_t xwidth = ( xhigh - xlow ) / nbins[depth-1];
  Double_t ywidth = ( yhigh - ylow ) / nbins[depth-1];

  if (debug1) std::cout << "depth: " << depth
			<< " ndepth: " << ndepth
			<< " nbins[depth-1]: " << nbins[depth-1]
			<< " x3: " << x3
			<< " at: " << alpha_t
			<< " xlow: " << xlow
			<< " xhigh: " << xhigh
			<< " ylow: " << ylow
			<< " yhigh: " << yhigh
			<< " xwidth: " << xwidth
			<< " ywidth: " << ywidth
			<< std::endl;

  // Loop through bins
  Int_t npass = 0;
  Int_t nfail = 0;
  Int_t nedge = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t nstraddles = 0;

  for ( Int_t xbin = 0; xbin < nbins[depth-1]; ++xbin ) { 
    for ( Int_t ybin = 0; ybin < nbins[depth-1]; ++ybin ) { 
      
      // (x,y) for bottom-left corner of bin
      Double_t xbl = xlow + xwidth * xbin;
      Double_t ybl = ylow + ywidth * ybin;

      // (x,y) for bottom-right corner of bin
      Double_t xbr = xlow + xwidth * (xbin+1);
      Double_t ybr = ylow + ywidth * ybin;

      // (x,y) for top-left corner of bin
      Double_t xtl = xlow + xwidth * xbin;
      Double_t ytl = ylow + ywidth * (ybin+1);

      // (x,y) for top-right corner of bin
      Double_t xtr = xlow + xwidth * (xbin+1);
      Double_t ytr = ylow + ywidth * (ybin+1);
      
      // Check if all corners of bin pass or fail the kinematic constraints
      bool pass = constrain(ybl,xbl,x3) && constrain(ybr,xbr,x3) && constrain(ytl,xtl,x3) && constrain(ytr,xtr,x3);
      bool fail = !constrain(ybl,xbl,x3) && !constrain(ybr,xbr,x3) && !constrain(ytl,xtl,x3) && !constrain(ytr,xtr,x3);
      bool edge = !pass && !fail;
      
      if (pass) npass++;
      if (fail) nfail++;
      if (edge) nedge++;
      
      // If all fail, do not consider bin
      if ( fail ) { continue; }

      // AlphaT values at corners of bin
      Double_t nbl = ybl + xbl - 1;
      Double_t nbr = ybr + xbr - 1;
      Double_t ntl = ytl + xtl - 1;
      Double_t ntr = ytr + xtr - 1;
      Double_t abl = nbl < 0. ? 1.e9 : xbl / ( 2. * sqrt( nbl ) ); 
      Double_t abr = nbr < 0. ? 1.e9 : xbr / ( 2. * sqrt( nbr ) ); 
      Double_t atl = ntl < 0. ? 1.e9 : xtl / ( 2. * sqrt( ntl ) ); 
      Double_t atr = ntr < 0. ? 1.e9 : xtr / ( 2. * sqrt( ntr ) ); 
      
      // Difference b/w values at corners and cut value
      Double_t dbl = abl - alpha_t;
      Double_t dbr = abr - alpha_t;
      Double_t dtl = atl - alpha_t;
      Double_t dtr = atr - alpha_t;

      // Check if bin falls below, above or straddles cut AlphaT value
      bool below = dbl<0. && dbr<0. && dtl<0. && dtr<0.;
      bool above = dbl>0. && dbr>0. && dtl>0. && dtr>0.;
      bool straddles = !above && !below;
      
      if (above) nabove++;
      if (below) nbelow++;
      if (straddles) nstraddles++;
      
      // Cross-section values at corners
      bool bbl = xbl < 1. && ybl < 1.;
      bool bbr = xbr < 1. && ybr < 1.;
      bool btl = xtl < 1. && ytl < 1.;
      bool btr = xtr < 1. && ytr < 1.;
      Double_t cbl = !bbl ? 0. : ( xbl*xbl + ybl*ybl ) / ( ( 1. - xbl ) * ( 1. - ybl ) ); 
      Double_t cbr = !bbr ? 0. : ( xbr*xbr + ybr*ybr ) / ( ( 1. - xbr ) * ( 1. - ybr ) ); 
      Double_t ctl = !btl ? 0. : ( xtl*xtl + ytl*ytl ) / ( ( 1. - xtl ) * ( 1. - ytl ) ); 
      Double_t ctr = !btr ? 0. : ( xtr*xtr + ytr*ytr ) / ( ( 1. - xtr ) * ( 1. - ytr ) ); 
      
      if ( ( edge || straddles ) && depth < ndepth ) {

	if (debug2) std::cout << " depth: " << depth << std::endl;

	if (debug2) std::cout << " xbin: " << xbin
			      << " ybin " << ybin
			      << " bl: " << xbl
			      << "," << ybl
			      << " br: " << xbr
			      << "," << ybr
			      << " tl: " << xtl
			      << "," << ytl
			      << " tr: " << xtr
			      << "," << ytr
			      << std::endl;

	if (debug2) std::cout << " pbl: " << constrain(ybl,xbl,x3)
			      << " pbr: " << constrain(ybr,xbr,x3)
			      << " ptl: " << constrain(ytl,xtl,x3)
			      << " ptr: " << constrain(ytr,xtr,x3)
			      << " pass: " << pass
			      << " fail: " << fail
			      << " edge: " << edge
			      << std::endl;

	if (debug2) std::cout << " abl: " << abl
			      << " abr: " << abr
			      << " atl: " << atl
			      << " atr: " << atr
			      << std::endl;

	if (debug2) std::cout << " dbl: " << dbl
			      << " dbr: " << dbr
			      << " dtl: " << dtl
			      << " dtr: " << dtr
			      << std::endl;
	
	if (debug2) std::cout << " above: " << above
			      << " below: " << below
			      << " straddles: " << straddles
			      << std::endl;
	
	if (debug2) std::cout << " cbl: " << cbl
			      << " cbr: " << cbr
			      << " ctl: " << ctl
			      << " ctr: " << ctr
			      << std::endl;
	
	// If bin at edge of kinematically allowed region or alpha_t cut is straddles bin, iterate again
	integrate( depth, ndepth, nbins, x3, alpha_t, numerator, denominator, xbl, xbr, ybl, ytl );
	
      } else {
	
	// Check bin widths
	if ( xwidth == 0. || ywidth == 0. ) { continue; }
	
 	// Add value at bottom-left corner to numerator and donominator 
 	//std::cout << " n: " << numerator << " d: " << denominator << std::endl;
	denominator += (cbl*xwidth*ywidth)/(nbins[depth-1]*nbins[depth-1]);
 	//std::cout << " d: depth/above/below/straddles: " << depth << "/" << above << "/" << below << "/" << straddles << std::endl;
	if ( abl > alpha_t ) {
	  numerator += (cbl*xwidth*ywidth)/(nbins[depth-1]*nbins[depth-1]); 
 	  //std::cout << " n: depth/above/below/straddles: " << depth << "/" << above << "/" << below << "/" << straddles << std::endl;
	}
 	//std::cout << " n: " << *numerator << " d: " << *denominator << std::endl;
      }
      
    } // loop
  } // loop
  
  if (debug1) std::cout << " depth: " << depth
			<< " npass: " << npass
			<< " nfail: " << nfail
			<< " nedge: " << nedge
			<< " nabove: " << nabove
			<< " nbelow: " << nbelow
			<< " nstraddles: " << nstraddles
			<< " numerator: " << numerator
 			<< " denominator: " << denominator
			<< std::endl;

}


