#include "common/style.C"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

typedef unsigned int uint;

// -----------------------------------------------------------------------------

typedef std::vector<std::string> StringV;
typedef std::vector<StringV> StringVV;

typedef std::vector<double> DoubleV;
typedef std::vector<DoubleV> DoubleVV;
typedef std::vector<DoubleVV> DoubleVVV;
typedef std::vector<DoubleVVV> DoubleVVVV;
typedef std::vector<DoubleVVVV> DoubleVVVVV;

typedef std::vector<int> IntV;
typedef std::vector<IntV> IntVV;
typedef std::vector<IntVV> IntVVV;
typedef std::vector<IntVVV> IntVVVV;
typedef std::vector<IntVVVV> IntVVVVV;

void resize( DoubleVVVVV& v, int nfile, int nmulti, int nat, int npt, int nht ) {
  v.clear();
  v.resize( nfile, DoubleVVVV( nmulti, DoubleVVV( nat, DoubleVV( npt, DoubleV( nht, 0. ) ) ) ) );
}

void resize( IntVVVV& v, int nfile, int nmulti, int nat, int npt ) {
  v.clear();
  v.resize( nfile, IntVVV( nmulti, IntVV( nat, IntV( npt, 0 ) ) ) );
}

// -----------------------------------------------------------------------------
// Poisson errors for n<20
void poissonErr( double x, double& errh, double& errl ) {
  double poisson_eh[6] = { 1.84, 2.30, 2.63, 2.92, 3.16, 3.38 };
  double poisson_el[6] = { 0.00, 0.83, 1.92, 1.63, 1.91, 2.16 };
  if ( x < 5. ) {
    int n = int(x);
    double f = x - double(int(x));
    errh = poisson_eh[n] + f*( poisson_eh[n+1] - poisson_eh[n] );
    errl = poisson_el[n] + f*( poisson_el[n+1] - poisson_el[n] );
  } else if ( x < 20. ) {
    // Gehrels (1986 ApJ, 303, 336) is an approx, with accuracy better than 2% 
    errh = 1.0 + sqrt(x+0.75); // equ 7, S=1
    errl = 0.0 + sqrt(x-0.25); // equ 11, S=1
  } else {
    // Normal approximation
    errh = sqrt(x);
    errl = sqrt(x);
  }
}

double poissonErrH( double x ) {
  double errh = 0., errl = 0.;
  poissonErr(x,errh,errl);
  return errh;
}

double poissonErrL( double x ) {
  double errh = 0., errl = 0.;
  poissonErr(x,errh,errl);
  return errl;
}

// -----------------------------------------------------------------------------
//  
int qcd( const std::vector<float>& at, const std::vector<float>& pval ) {
  
  if ( at.empty() || pval.empty() || at.size() != pval.size() ) { return -1; }
  
  setTDRStyle();
  
  TCanvas* c1 = new TCanvas("pValueVersusAlphaT","");
  c1->SetFillColor(0);
  c1->SetLineColor(0); 
  //c1->SetRightMargin(0.2);
  
//   TLegend* legend = new TLegend( 0.82, 0.5, 0.98, 0.9, NULL, "brNDC" );
//   legend->SetFillColor(0);
//   legend->SetLineColor(0); 
//   legend->SetShadowColor(0); 
  
  TMultiGraph* mg = new TMultiGraph();
  for ( uint ii = 0; ii < 1; ++ii ) {
    
    TGraphAsymmErrors* gr = 0;
    
    std::vector<float> at_e(at.size(),0.);
    std::vector<float> pval_e(pval.size(),0.);
    
    gr = new TGraphAsymmErrors(at.size(),
			       &at.front(),
			       &pval.front(),
			       &at_e.front(),
			       &at_e.front(),
			       &pval_e.front(),
			       &pval_e.front()
			       );

    mg->Add(gr,"p");
    gr->SetTitle("pValue");
    gr->SetLineColor(1+ii);
    gr->SetMarkerColor(1+ii);
    gr->SetMarkerStyle(20+ii);
    //legend->AddEntry( gr, TString("pValue"), "lep" );
  }
    
  c1->cd();
  mg->Draw("al");
  mg->GetXaxis()->SetTitle("#alpha_{T} cut value");
  mg->GetYaxis()->SetTitle("p-value");
  mg->GetXaxis()->SetRangeUser(0.5,0.6);
  mg->GetYaxis()->SetRangeUser(0.,1.);
  //legend->Draw("same");
  c1->Update();
  c1->SaveAs(TString("pValueVersusAlphaT.pdf")); 
  c1->SaveAs(TString("pValueVersusAlphaT.C")); 
  
}

// -----------------------------------------------------------------------------
//  
double poissonf( double* x, double* par ) {                                                                              
  return par[0] * TMath::Poisson(x[0],par[1]);
}                                                                              

// -----------------------------------------------------------------------------
// 
void checkStats( std::vector<std::string>& type,
		 std::vector<int>& multi,
		 std::vector<int>& at, 
		 uint ifile, uint imulti, uint iat, 
		 double observed, double estimate, double errh, double errl, 
		 double b0_pass, double b0_fail, 
		 double b1_pass, double b1_fail,
		 double b2_fail,
		 double& b2_mean, 
		 double& b2_rms, 
		 double& b2_err ) {

  TDirectory* dir = gDirectory;
  TStyle* style = gStyle;
  gStyle->SetOptStat(1111111);
  
  std::stringstream ss;
  ss << "_" << type[ifile] << "_multi" << multi[imulti] << "_aT" << at[iat]/1000.;
  TCanvas* c2 = new TCanvas(TString("CheckStats"+ss.str()),
			    TString("CheckStats"+ss.str()));
  
  uint number_of_trials = 10000;
  
  // Poisson error
  double b0_pass_err = sqrt(b0_pass);
  double b1_pass_err = sqrt(b1_pass);
  
  // Gaussian error
  double b0_fail_err = sqrt(b0_fail);
  double b1_fail_err = sqrt(b1_fail);

  double xmin = 0.;
  double xmax = 0.;
  
  xmin = 1. * int( b0_fail - 5. * b0_fail_err );
  xmax = 1. * int( b0_fail + 5. * b0_fail_err );
  TH1D* b0_fail_his = new TH1D(TString("b0_fail"+ss.str()),"",int(xmax-xmin),xmin,xmax);

  xmin = 1. * int( b1_fail - 5. * b1_fail_err );
  xmax = 1. * int( b1_fail + 5. * b1_fail_err );
  TH1D* b1_fail_his = new TH1D(TString("b1_fail"+ss.str()),"",int(xmax-xmin),xmin,xmax);
  
  xmin = -1.;
  xmax = 1. * int( b0_pass + 10. * b0_pass_err );
  TH1D* b0_pass_his = new TH1D(TString("b0_pass"+ss.str()),"",int(xmax-xmin),xmin,xmax);
  
  xmin = -1.;
  xmax = 1. * int( b1_pass + 10. * b1_pass_err );
  TH1D* b1_pass_his = new TH1D(TString("b1_pass"+ss.str()),"",int(xmax-xmin),xmin,xmax);
  
  xmin = -1.;
  xmax = 1. * int( b1_pass + 10. * b1_pass_err );
  TH1D* b2_pred_his = new TH1D(TString("b2_pred"+ss.str()),"",int(xmax-xmin),xmin,xmax);
  
  TRandom3 rand;
  for ( uint ii = 0; ii < number_of_trials; ++ii ) {

    double b0_pass_rand = rand.PoissonD(b0_pass);
    double b1_pass_rand = rand.PoissonD(b1_pass);
    double b0_fail_rand = rand.Gaus(b0_fail,b0_fail_err);
    double b1_fail_rand = rand.Gaus(b1_fail,b1_fail_err);
    
    b0_pass_his->Fill(b0_pass_rand);
    b1_pass_his->Fill(b1_pass_rand);
    b0_fail_his->Fill(b0_fail_rand);
    b1_fail_his->Fill(b1_fail_rand);

    double b0_ratio = ( b0_fail_rand > 0. ? b0_pass_rand/b0_fail_rand : 0. );
    double b1_ratio = ( b1_fail_rand > 0. ? b1_pass_rand/b1_fail_rand : 0. );
    
    double ratio = ( b0_ratio > 0. && b1_ratio > 0. ? b1_ratio/b0_ratio : 0. );
    double b2_ratio = b1_ratio * ratio;
    
    double b2_pred_rand = b2_ratio * b2_fail;
    if ( ratio > 0. ) { b2_pred_his->Fill(b2_pred_rand); }
    
  }

  // Pass/fail distributions
  c2->Divide(2,3);
  c2->cd(1);
  b0_pass_his->Draw();
  c2->cd(2);
  b0_fail_his->Draw();
  c2->cd(3);
  b1_pass_his->Draw();
  c2->cd(4);
  b1_fail_his->Draw();
  
  // Predicted value (with fit)
  c2->cd(5);
  TF1 pois("pois",poissonf,xmin,xmax,2);
  pois.SetParName(0,"Const");
  pois.SetParName(1,"Mean");
  pois.SetParameter(0,b2_pred_his->GetMaximum());
  pois.SetParameter(1,b2_pred_his->GetMean());
  gStyle->SetOptFit(1111111);
  b2_pred_his->Fit("pois","Q");
  b2_pred_his->Draw();

  // 68% integral
  double denominator = b2_pred_his->Integral();
  b2_pred_his->SetAxisRange(estimate-errl,estimate+errh);
  double numerator = b2_pred_his->Integral();
  b2_pred_his->SetAxisRange(xmin,xmax);
  
  c2->SaveAs(TString("CheckStats"+ss.str()+".png"));
  
  b2_mean = pois.GetParameter(1);
  b2_rms = sqrt(pois.GetParameter(1));
  b2_err = pois.GetParError(1);
  
  dir->cd();
  gStyle = style;
  
  std::cout << " STAT CHECK:"
	    << " sample:\"" << type[ifile] << "\""
	    << " multi: " << multi[imulti]
	    << " aT:" << at[iat]/1000.
	    << " obs: " << observed
	    << "+/-" << sqrt(observed)
	    << " pred: " << estimate
	    << "+" << errh
	    << "-" << errl
	    << " stat: " << b2_mean
	    << "+/-" << b2_rms
	    << " (" << (int(( denominator > 0. ? numerator/denominator : 0. )*1000.))/10. << "%)"
	    << std::endl;
  
}

// -----------------------------------------------------------------------------
// 
double dr( double x, int decimal_places ) {
  if ( decimal_places < 0 ) { return x; }
  if ( x == 0. ) { return 0.; }
  return floor(x*pow(10.,decimal_places)+0.5)/pow(10.,decimal_places); 
}

// -----------------------------------------------------------------------------
// 
double sr( double x, int significant_figures ) {
  if ( significant_figures < 0 ) { return x; }
  if ( x == 0. ) { return 0.; }
  double s = floor(log10(x))-significant_figures+1;
  double f = pow(10.,fabs(s));
//   if ( f == 0. || 1/f == 0. ) {
//     std::cout << " x: " << x
// 	      << " n: " << significant_figures
// 	      << " s: " << s
// 	      << " abs(s): " << abs(s)
// 	      << " pow(10.,abs(s)): " << pow(10.,abs(s)) 
// 	      << std::endl;
//   }
  if (s<0) f = 1/f;
  return int(x/f)*f;
}

// -----------------------------------------------------------------------------
// 
void calcErr( double val, double& errh, double& errl, 
	      double val1, double errh1, double errl1,
	      double val2, double errh2, double errl2,
	      bool print, std::stringstream& ss ) {
  errh = 0.;
  errl = 0.;
  if ( val1 > 0. ) {
    errh = (errh1/val1)*(errh1/val1);
    errl = (errl1/val1)*(errl1/val1);
//    } else {
//      errh = errh1*errh1; 
//      errl = errl1*errl1;
  }
  if ( val2 > 0. ) {
    errh += (errh2/val2)*(errh2/val2);
    errl += (errl2/val2)*(errl2/val2);
//    } else {
//      errh += errh2*errh2;
//      errl += errl2*errl2;
  }
  if ( val > 0 ) {
    errh = val * sqrt(errh);
    errl = val * sqrt(errl);
//    } else {
//      errh = sqrt(errh);
//      errl = sqrt(errl);
  }

//   if ( !( val1 > 0. && val2 > 0. ) ) {
//     errh = 0.;
//     errl = 0.;
//   }

  if ( print ) {
    ss << " out: " << val
       << "+" << errh
       << "-" << errl
       << " (" 
       << dr( (val>0.?errh/val:0.), 2 ) 
       << "," 
       << dr( (val>0.?errl/val:-1.), 2 )
       << ") " 
       << " in1: " << val1
       << "+" << errh1
       << "-" << errl1
       << " (" 
       << dr( (val1>0.?errh1/val1:0.), 2 ) 
       << "," 
       << dr( (val1>0.?errl1/val1:0.), 2 ) 
       << ") " 
       << " in2: " << val2
       << "+" << errh2
       << "-" << errl2
       << " (" 
       << dr( (val2>0.?errh2/val2:0.), 2 )
       << "," 
       << dr( (val2>0.?errl2/val2:0.), 2 )
       << ") ";
  }

}

// -----------------------------------------------------------------------------
// 
void calcErr( double val, double& errh, double& errl, 
	      double val1, double errh1, double errl1,
	      double val2, double errh2, double errl2 ) {
  std::stringstream null;
  calcErr( val, errh, errl, 
	   val1, errh1, errl1,
	   val2, errh2, errl2, 
	   false, null );
}

// -----------------------------------------------------------------------------
//  
TH1D* rebin( TFile* file, const char* name, int multi, double width, double xlow, double xhigh ) {

  // Check for null pointer
  if ( !file ) { 
    std::cout << "Unable to retrieve file!" << std::endl;
    return 0; 
  }
  file->cd();

  // Retreive original histogram
  TH1D* input = 0;
  
  if ( multi == 0 ) {
    std::stringstream tmp; tmp << name << "_all";
    //std::cout << "name: " << tmp.str() << std::endl;
    input = (TH1D*)file->Get( tmp.str().c_str() )->Clone();
  } else if ( multi < 0 ) {
    for ( uint ii = abs(multi); ii < 10; ++ii ) {
      std::stringstream tmp; tmp << name << "_" << ii;
      //std::cout << "name: " << tmp.str() << std::endl;
      if ( int(ii) == abs(multi) ) { input = (TH1D*)file->Get( tmp.str().c_str() )->Clone(); }
      else { input->Add( (TH1D*)file->Get( tmp.str().c_str() ) ); }
    }
  } else {
    std::stringstream tmp; tmp << name << "_" << multi;
    input = (TH1D*)file->Get( tmp.str().c_str() )->Clone();
  }
  
  if ( !input ) { 
    std::cout << "Unable to retrieve histogram with name \"" << name << "\"" << std::endl;
    return 0; 
  }

  if (true) {
    if ( input ) {
      TH1D* his = input;
      std::stringstream ss;
      ss << "INPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }
  
  // Rebin original histogram
  input->Rebin( int(width/(input->GetBinWidth(1))) );
  
  // Create new rebinned histogram
  TH1D* output = new TH1D( TString(input->GetName())+"_Rebinned", "", int((xhigh+width-xlow)/width), xlow, xhigh+width );
  
  // Populate new histogram
  for ( int bin = 1; bin <= input->GetNbinsX(); ++bin ) {
    Double_t centre = input->GetBinLowEdge(bin) + input->GetBinWidth(bin)/2.;
    output->Fill( centre, input->GetBinContent(bin) );
  }
  // Modify last HT bin to be inclusive
  if ( true ) { 
    Double_t centre = output->GetBinLowEdge(output->GetNbinsX()) + output->GetBinWidth(output->GetNbinsX())/2.;
    output->Fill( centre, input->GetBinContent(output->GetNbinsX()+1) );
  }
  
  output->Sumw2();
  
  if (true) {
    if ( input ) {
      TH1D* his = input;
      std::stringstream ss;
      ss << "REBINNED HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }
  
  if (true) {
    if ( output ) {
      TH1D* his = output;
      std::stringstream ss;
      ss << "OUTPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }

  return output;

}

// // -----------------------------------------------------------------------------
// //  
// TH1D* rebinNew( TFile* file, const char* name, int multi, int nbins = -1, double xlow = -1., double xhigh = -1. ) {

//   bool debug = true;

//   // Check for null pointer
//   if ( file ) { 
//     //std::cout << "Opened file with name \"" << file->GetName() << "\"" << std::endl;
//   } else {
//     std::cout << "Unable to retrieve file!" << std::endl;
//     return 0; 
//   }
//   file->cd();

//   // Retreive original histogram
//   TH1D* input = 0;
  
//   if ( multi == 0 ) {
//     std::stringstream tmp; tmp << name << "_all";
//     input = (TH1D*)file->Get( tmp.str().c_str() )->Clone();
//   } else if ( multi < 0 ) {
//     for ( uint ii = abs(multi); ii < 9; ++ii ) {
//       std::stringstream tmp; tmp << name << "_" << ii;
//       if ( int(ii) == abs(multi) ) { input = (TH1D*)file->Get( tmp.str().c_str() )->Clone(); }
//       else { input->Add( (TH1D*)file->Get( tmp.str().c_str() ) ); }
//     }
//   } else {
//     std::stringstream tmp; tmp << name << "_" << multi;
//     input = (TH1D*)file->Get( tmp.str().c_str() )->Clone();
//   }
  
//   if ( input ) { 
//     //std::cout << "Retrieved histogram with name \"" << input->GetName() << "\"" << std::endl;
//   } else {
//     std::cout << "Unable to retrieve histogram with name \"" << name << "\"" << std::endl;
//     return 0; 
//   }

//   if (debug) {
//     if ( input ) {
//       TH1D* his = input;
//       std::stringstream ss;
//       ss << "INPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
// 	 << std::endl
// 	 << " underflow: " << his->GetBinContent(0)
// 	 << std::endl
// 	 << " underflow_err: " << his->GetBinError(0)
// 	 << std::endl
// 	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " entries: " << his->GetEntries()
// 	 << std::endl
// 	 << " integral: " << his->Integral()
// 	 << std::endl
// 	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl;
//       for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
// 	ss << " bin: " << bin
// 	   << " val: " << his->GetBinContent(bin)
// 	   << " err: " << his->GetBinError(bin)
// 	   << std::endl;
//       }
//       std::cout << ss.str();
//     }
//   }
  
//   // Rebin original histogram
//   if (debug) std::cout << " INPUT: " << input->GetName()
// 		       << " Nbins: " << input->GetXaxis()->GetNbins()
// 		       << " Xmin: " << input->GetXaxis()->GetXmin()
// 		       << " Xmax: " << input->GetXaxis()->GetXmax()
// 		       << std::endl;
//   if ( nbins > 0. ) { input->Rebin( int(input->GetXaxis()->GetNbins()/nbins) ); }
//   if (debug) std::cout << " REBINNED: " << input->GetName()
// 		       << " Nbins: " << input->GetXaxis()->GetNbins()
// 		       << " Xmin: " << input->GetXaxis()->GetXmin()
// 		       << " Xmax: " << input->GetXaxis()->GetXmax()
// 		       << std::endl;
  
//   // Create temp histogram
//   TH1D* temp = 0;
//   if ( nbins > 1 ) { 
//     int bins = nbins - 1;
//     double xmin = xlow;
//     double xmax = xhigh;
//     double width = ( xmax - xmin ) / ( bins * 1. );
//     xmax -= width;
//     temp = new TH1D( TString(input->GetName())+"_Temp", "", 
// 		     bins, 
// 		     xmin, 
// 		     xmax );
//     if (debug) std::cout << " TEMP: " << temp->GetName()
// 			 << " Nbins: " << temp->GetXaxis()->GetNbins()
// 			 << " Xmin: " << temp->GetXaxis()->GetXmin()
// 			 << " Xmax: " << temp->GetXaxis()->GetXmax()
// 			 << std::endl;
//   } else {
//     int bins = temp->GetXaxis()->GetNbins() - 1;
//     double xmin = temp->GetXaxis()->GetXmin();
//     double xmax = temp->GetXaxis()->GetXmax();
//     double width = ( xmax - xmin ) / ( bins * 1. );
//     xmax -= width;
//     temp = new TH1D( TString(input->GetName())+"_Temp", "", 
// 		     bins, 
// 		     xmin, 
// 		     xmax );
//     if (debug) std::cout << " TEMP: " << temp->GetName()
// 			 << " Nbins: " << temp->GetXaxis()->GetNbins()
// 			 << " Xmin: " << temp->GetXaxis()->GetXmin()
// 			 << " Xmax: " << temp->GetXaxis()->GetXmax()
// 			 << std::endl;
//   }
//   temp->Sumw2();
  
//   // Populate temp histogram
//   for ( int bin = 1; bin <= temp->GetNbinsX(); ++bin ) {
//     Double_t centre = temp->GetBinLowEdge(bin) + temp->GetBinWidth(bin)/2.;
//     int find_bin = input->FindBin(centre);
//     temp->SetBinContent( bin, input->GetBinContent( find_bin ) );
//     temp->SetBinError( bin, input->GetBinError( find_bin ) );
//   }
//   temp->SetBinContent( bin, input->GetBinContent( input->GetXaxis()->GetNbins()+1 ) );
//   temp->SetBinError( bin, input->GetBinError( input->GetXaxis()->GetNbins()+1 ) );
//   temp->SetEntries( input->GetEntries() );

//   // Create new rebinned histogram
//   TH1D* output = 0;
//   if ( nbins > 0 ) { 
//     output = new TH1D( TString(input->GetName())+"_Output", "", 
// 		       nbins, 
// 		       xlow, 
// 		       xhigh );
//     if (debug) std::cout << " OUTPUT: " << output->GetName()
// 			 << " Nbins: " << output->GetXaxis()->GetNbins()
// 			 << " Xmin: " << output->GetXaxis()->GetXmin()
// 			 << " Xmax: " << output->GetXaxis()->GetXmax()
// 			 << std::endl;
//   } else {
//     output = new TH1D( TString(input->GetName())+"_Output", "", 
// 		       input->GetXaxis()->GetNbins(), 
// 		       input->GetXaxis()->GetXmin(), 
// 		       input->GetXaxis()->GetXmax() );
//     if (debug) std::cout << " OUTPUT: " << output->GetName()
// 			 << " Nbins: " << output->GetXaxis()->GetNbins()
// 			 << " Xmin: " << output->GetXaxis()->GetXmin()
// 			 << " Xmax: " << output->GetXaxis()->GetXmax()
// 			 << std::endl;
//   }
//   output->Sumw2();

//   // Populate new histogram
//   for ( int bin = 1; bin <= output->GetNbinsX(); ++bin ) {
//     Double_t centre = output->GetBinLowEdge(bin) + output->GetBinWidth(bin)/2.;
//     int find_bin = temp->FindBin(centre);
//     output->SetBinContent( bin, temp->GetBinContent( find_bin ) );
//     output->SetBinError( bin, temp->GetBinError( find_bin ) );
//   }
//   output->SetEntries( temp->GetEntries() );

//   if (debug) {
//     if ( temp ) {
//       TH1D* his = temp;
//       std::stringstream ss;
//       ss << "TEMP HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
// 	 << std::endl
// 	 << " underflow: " << his->GetBinContent(0)
// 	 << std::endl
// 	 << " underflow_err: " << his->GetBinError(0)
// 	 << std::endl
// 	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " entries: " << his->GetEntries()
// 	 << std::endl
// 	 << " integral: " << his->Integral()
// 	 << std::endl
// 	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl;
//       for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
// 	ss << " bin: " << bin
// 	   << " val: " << his->GetBinContent(bin)
// 	   << " err: " << his->GetBinError(bin)
// 	   << std::endl;
//       }
//       std::cout << ss.str();
//     }
//   }
  
// //   // Modify last bin to be inclusive
// //   if ( true ) { 
// //     // Define value of centre of last bin of output histo
// //     Double_t centre = output->GetBinLowEdge(output->GetNbinsX()) + output->GetBinWidth(output->GetNbinsX())/2.;
// //     // Add to the last bin of output histo with overflow from output histo
// //     output->SetBinContent( output->GetNbinsX(), 
// // 			   output->GetBinContent(output->GetNbinsX()) +
// // 			   output->GetBinContent(output->GetNbinsX()+1) );
// //     output->SetBinError( output->GetNbinsX(), 
// // 			 sqrt(output->GetBinError(output->GetNbinsX()) *
// // 			      output->GetBinError(output->GetNbinsX()) +
// // 			      output->GetBinError(output->GetNbinsX()+1) *
// // 			      output->GetBinError(output->GetNbinsX()+1)) );
// //     // Add to the last bin of output histo with overflow from input histo
// //     output->SetBinContent( output->GetNbinsX(), 
// // 			   output->GetBinContent(output->GetNbinsX()) +
// // 			   input->GetBinContent(output->GetNbinsX()+1) );
// //     output->SetBinError( output->GetNbinsX(), 
// // 			 sqrt(output->GetBinError(output->GetNbinsX()) *
// // 			      output->GetBinError(output->GetNbinsX()) +
// // 			      output->GetBinError(input->GetNbinsX()+1) *
// // 			      output->GetBinError(input->GetNbinsX()+1)) );
// //   }
  
//   if (debug) {
//     if ( output ) {
//       TH1D* his = output;
//       std::stringstream ss;
//       ss << "OUTPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
// 	 << std::endl
// 	 << " underflow: " << his->GetBinContent(0)
// 	 << std::endl
// 	 << " underflow_err: " << his->GetBinError(0)
// 	 << std::endl
// 	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
// 	 << std::endl
// 	 << " entries: " << his->GetEntries()
// 	 << std::endl
// 	 << " integral: " << his->Integral()
// 	 << std::endl
// 	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
// 	 << std::endl;
//       for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
// 	ss << " bin: " << bin
// 	   << " val: " << his->GetBinContent(bin)
// 	   << " err: " << his->GetBinError(bin)
// 	   << std::endl;
//       }
//       std::cout << ss.str();
//     }
//   }

//   if (temp) { delete temp; }
//   return output;

// }

// -----------------------------------------------------------------------------
//  
//TH1D* rebinNew( TFile* file, const char* name, int multi, int nbins = -1, double xlow = -1., double xhigh = -1. ) {
TH1D* rebinNew( TFile* file, const char* name, int multi, int nbins = -1, double* xarray = 0 ) {

  bool debug = true;

  // Check for null pointer
  if ( file ) { 
    //std::cout << "Opened file with name \"" << file->GetName() << "\"" << std::endl;
  } else {
    std::cout << "Unable to retrieve file!" << std::endl;
    return 0; 
  }
  file->cd();

  // Retreive original histogram
  TH1D* input = 0;
  
  if ( multi == 0 ) {
    std::stringstream tmp; tmp << name << "_all";
    TH1D* temp = (TH1D*)file->Get( tmp.str().c_str() );
    if (temp) { input = (TH1D*)temp->Clone(); } 
    //else { std::cout << "Unable to retrieve histo with name " << tmp.str() << std::endl; }
  } else if ( multi < 0 ) {
    for ( uint ii = abs(multi); ii <= 20; ++ii ) {
      std::stringstream tmp; tmp << name << "_" << ii;
      if ( int(ii) == abs(multi) ) { 
	TH1D* temp = (TH1D*)file->Get( tmp.str().c_str() );
	if (temp) { input = (TH1D*)temp->Clone(); } 
	//else { std::cout << "Unable to retrieve histo with name " << tmp.str() << std::endl; }
      }
      else { 
	TH1D* temp = (TH1D*)file->Get( tmp.str().c_str() );
	if (temp) { input->Add( (TH1D*)temp ); } 
	//else { std::cout << "Unable to retrieve histo with name " << tmp.str() << std::endl; }
      }
    }
  } else {
    std::stringstream tmp; tmp << name << "_" << multi;
    TH1D* temp = (TH1D*)file->Get( tmp.str().c_str() );
    if (temp) { input = (TH1D*)temp->Clone(); } 
    //else { std::cout << "Unable to retrieve histo with name " << tmp.str() << std::endl; }
  }
  
  if ( input ) { 
    //std::cout << "Retrieved histogram with name \"" << input->GetName() << "\"" << std::endl;
  } else {
    //std::cout << "Unable to retrieve histogram with name \"" << name << "\"" << std::endl;
    return 0; 
  }

  if (debug) {
    if ( input ) {
      TH1D* his = input;
      std::stringstream ss;
      ss << "INPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " underflow_err: " << his->GetBinError(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }
  
  // Rebin original histogram
//   std::cout << " INPUT: Nbins: " << input->GetXaxis()->GetNbins()
// 	    << " Xmin: " << input->GetXaxis()->GetXmin()
// 	    << " Xmax: " << input->GetXaxis()->GetXmax()
// 	    << std::endl;
//   if ( nbins > 0. ) {
//     input->Rebin( int(input->GetXaxis()->GetNbins()/nbins) );
//   }
//   std::cout << " INPUT: Nbins: " << input->GetXaxis()->GetNbins()
// 	    << " Xmin: " << input->GetXaxis()->GetXmin()
// 	    << " Xmax: " << input->GetXaxis()->GetXmax()
// 	    << std::endl;

  if (debug) {
    if ( input ) {
      TH1D* his = input;
      std::stringstream ss;
      ss << "REBINNED HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " underflow_err: " << his->GetBinError(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }
  
  // Create new rebinned histogram
  TH1D* output = 0;
//   if ( nbins > 0. ) { 
//     output = new TH1D( TString(input->GetName())+"_Rebinned", "", 
// 		       nbins, 
// 		       xlow, 
// 		       xhigh );
// //     std::cout << " OUTPUT: Nbins: " << output->GetXaxis()->GetNbins()
// // 	      << " Xmin: " << output->GetXaxis()->GetXmin()
// // 	      << " Xmax: " << output->GetXaxis()->GetXmax()
// // 	      << std::endl;
//   } else {
//     output = new TH1D( TString(input->GetName())+"_Rebinned", "", 
// 		       input->GetXaxis()->GetNbins(), 
// 		       input->GetXaxis()->GetXmin(), 
// 		       input->GetXaxis()->GetXmax() );
// //     std::cout << " OUTPUT: Nbins: " << output->GetXaxis()->GetNbins()
// // 	      << " Xmin: " << output->GetXaxis()->GetXmin()
// // 	      << " Xmax: " << output->GetXaxis()->GetXmax()
// // 	      << std::endl;
//   }
//   if (output) output->Sumw2();

  // Populate new histogram
  if (!output) output = input;
  for ( int bin = 1; bin <= output->GetNbinsX(); ++bin ) {
    Double_t centre = output->GetBinLowEdge(bin) + output->GetBinWidth(bin)/2.;
    int find_bin = input->FindBin(centre);
    output->SetBinContent( bin, input->GetBinContent( find_bin ) );
    output->SetBinError( bin, input->GetBinError( find_bin ) );
  }
  output->SetEntries( input->GetEntries() );
  
  // Modify last bin to be inclusive
  if ( true ) { 
    // Add to the last bin of output histo with overflow from output histo
//      std::cout << " content: " << output->GetBinContent( output->GetNbinsX() )
//  	      << " error: " << output->GetBinError( output->GetNbinsX() )  
//  	      << std::endl;
    output->SetBinContent( output->GetNbinsX(), 
			   output->GetBinContent(output->GetNbinsX()) +
			   output->GetBinContent(output->GetNbinsX()+1) );
    output->SetBinError( output->GetNbinsX(), 
			 sqrt(output->GetBinError(output->GetNbinsX()) *
			      output->GetBinError(output->GetNbinsX()) +
			      output->GetBinError(output->GetNbinsX()+1) *
			      output->GetBinError(output->GetNbinsX()+1)) );
    output->SetBinContent( output->GetNbinsX()+1, 0. );
    output->SetBinError( output->GetNbinsX()+1, 0. );
//      std::cout << " content: " << output->GetBinContent( output->GetNbinsX() )
//  	      << " error: " << output->GetBinError( output->GetNbinsX() )  
//  	      << std::endl;
//     // Add to the last bin of output histo with overflow from input histo
//     output->SetBinContent( output->GetNbinsX(), 
// 			   output->GetBinContent(output->GetNbinsX()) +
// 			   input->GetBinContent(input->GetNbinsX()+1) );
//     output->SetBinError( output->GetNbinsX(), 
// 			 sqrt(output->GetBinError(output->GetNbinsX()) *
// 			      output->GetBinError(output->GetNbinsX()) +
// 			      input->GetBinError(input->GetNbinsX()+1) *
// 			      input->GetBinError(input->GetNbinsX()+1)) );
//      std::cout << " content: " << output->GetBinContent( output->GetNbinsX() )
// 	       << " error: " << output->GetBinError( output->GetNbinsX() )  
//  	      << std::endl;
  }
  
  if (debug) {
    if ( output ) {
      TH1D* his = output;
      std::stringstream ss;
      ss << "OUTPUT HISTO: " << file->GetName() << "/" << name << " multi: " << (multi<0?">=":"==") << abs(multi) 
	 << std::endl
	 << " underflow: " << his->GetBinContent(0)
	 << std::endl
	 << " underflow_err: " << his->GetBinError(0)
	 << std::endl
	 << " overflow: " << his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl
	 << " overflow_err: " << his->GetBinError(his->GetNbinsX()+1)
	 << std::endl
	 << " entries: " << his->GetEntries()
	 << std::endl
	 << " integral: " << his->Integral()
	 << std::endl
	 << " all: " << his->Integral() + his->GetBinContent(0) + his->GetBinContent(his->GetNbinsX()+1)
	 << std::endl;
      for ( int bin = 1; bin <= his->GetNbinsX(); ++bin ) {
	ss << " bin: " << bin
	   << " val: " << his->GetBinContent(bin)
	   << " err: " << his->GetBinError(bin)
	   << std::endl;
      }
      std::cout << ss.str();
    }
  }

  return output;

}

// -----------------------------------------------------------------------------
/*
 */  
TH1D* integral( TH1D* input ) {
  
  // Check input
  if ( !input ) { return 0; }

  // Create clone
  TString name( TString( input->GetName() ) + "_Integral" );
  TH1D* output = new TH1D( name, "", 
			   input->GetNbinsX(), 
			   input->GetXaxis()->GetXmin(), 
			   input->GetXaxis()->GetXmin() );
  
  // Integral
  double total = 0.;
  
  // Overflow
  total += input->GetBinContent( input->GetNbinsX() );
  
  // Loop through histo bins
  for ( int bin = 1; bin < input->GetNbinsX()+1; ++bin ) {
    Double_t centre = input->GetBinLowEdge(bin) + input->GetBinWidth(bin)/2.;
    total += input->GetBinContent( bin );
    output->Fill( centre, total );
  }
  output->Sumw2();

  std::cout << "total: " << total 
	    << " integral: " << input->Integral()
	    << " integral: " << output->Integral()
	    << std::endl;

  return output;
}

// -----------------------------------------------------------------------------
// 
void calcRatio( const uint nfile, 
		const uint nmulti, 
		const uint nat, 
		const uint npt, 
		const uint nht, 
 		StringV his, 
 		StringVV files,
 		DoubleV lumis,
 		DoubleV weights,
 		IntV multi, 
 		IntV at, 
 		DoubleV pt, 
 		DoubleV ht, 
 		int    ht_nbin,
 		double ht_min, 
 		double ht_max,
 		DoubleVVVVV& numer, 
 		DoubleVVVVV& numer_errh, 
 		DoubleVVVVV& numer_errl, 
 		DoubleVVVVV& denom, 
 		DoubleVVVVV& denom_errh, 
 		DoubleVVVVV& denom_errl, 
 		DoubleVVVVV& ratio, 
 		DoubleVVVVV& errh, 
 		DoubleVVVVV& errl,
 		IntVVVV& length,
 		int data_file,
 		double lumi,
		bool efficiency = false
		) {
  
  // Reset
  for ( uint ifile = 0; ifile < nfile; ++ifile ) {
    for ( uint imulti = 0; imulti < nmulti; ++imulti ) {
      for ( uint iat = 0; iat < nat; ++iat ) {	
	for ( uint ipt = 0; ipt < npt; ++ipt ) {
	  for ( uint iht = 0; iht < nht; ++iht ) {
	    numer[ifile][imulti][iat][ipt][iht] = 0.;
	    numer_errh[ifile][imulti][iat][ipt][iht] = 0.;
	    numer_errl[ifile][imulti][iat][ipt][iht] = 0.;
	    denom[ifile][imulti][iat][ipt][iht] = 0.;
	    denom_errh[ifile][imulti][iat][ipt][iht] = 0.; 
	    denom_errl[ifile][imulti][iat][ipt][iht] = 0.; 
	    ratio[ifile][imulti][iat][ipt][iht] = 0.;
	    length[ifile][imulti][iat][ipt] = 0; 
	    errh[ifile][imulti][iat][ipt][iht] = 0.;
	    errl[ifile][imulti][iat][ipt][iht] = 0.;
	  }
	}
      }
    }
  }

  // Loop through collections of data/MC files 
  for ( uint ifile = 0; ifile < nfile; ++ifile ) {
    
    // Loop through data/MC files with each collection
    for ( uint jfile = 0; jfile < files[ifile].size(); ++jfile ) {
      
      // Open file
      std::cout << " ifile: " << ifile
		<< " jfile: " << jfile
		<< " file: " << files[ifile][jfile]
		<< std::endl;
      TString filename(files[ifile][jfile]);
      TFile* file = new TFile(filename);
      if ( file->IsZombie() || !(file->IsOpen()) ) { 
	if (file) { delete file; }
	continue; 
      }
      file->cd();
      
      // Loop through bins of multiplicity
      for ( uint imulti = 0; imulti < nmulti; ++imulti ) {

	// Loop through bins of AlphaT
	for ( uint iat = 0; iat < nat; ++iat ) {
      
	  // Loop through bins of pT
	  for ( uint ipt = 0; ipt < npt; ++ipt ) {
	    
	    length[ifile][imulti][iat][ipt] = 0;

	    // Directory and histo names
	    std::stringstream dir;
	    std::stringstream name;
	    std::stringstream pre;
	    std::stringstream post;
	    dir << "QcdBkgdEst/";
	    name << his[ifile] << "_aT"; 
	    pre << dir.str() << name.str() << "0";
	    post << dir.str() << name.str() << at[iat]/1000.;
	
	    // Create ratio histo
	    TH1D* above = rebinNew( file, post.str().c_str(), multi[imulti] );
	    TH1D* below = rebinNew( file, pre.str().c_str(), multi[imulti] );
	  
	    file->cd();
	    
	    if ( above && below ) {

	      above->Scale(lumis[ifile]/100.);
	      below->Scale(lumis[ifile]/100.);

	      if ( efficiency == false ) { below->Add(above,-1.); }
	      
	      if ( false ) { 
		std::cout << " CONTENTS: file: " << files[ifile][jfile] 
			  << " multi: " << multi[imulti]
			  << " aT:" << at[iat]
			  << "  HT:" << ht[ipt]
			  << " n: " 
			  << above->GetBinContent(ipt+1) 
			  << " ne: " 
			  << above->GetBinError(ipt+1) 
			  << " d: " 
			  << below->GetBinContent(ipt+1) 
			  << " de: " 
			  << below->GetBinError(ipt+1) 
			  << std::endl;
	      }

	      for ( uint iht = 0; iht < nht; ++iht ) {
		
		double a  = above->GetBinContent(iht+1);
		double aeh = a > 20. ? above->GetBinError(iht+1) : poissonErrH(a);
		double ael = a > 20. ? above->GetBinError(iht+1) : poissonErrL(a);;
		double b  = below->GetBinContent(iht+1);
		double beh = b > 20. ? below->GetBinError(iht+1) : poissonErrH(b);
		double bel = b > 20. ? below->GetBinError(iht+1) : poissonErrL(b);
		
		numer[ifile][imulti][iat][ipt][iht]      += a;
		numer_errh[ifile][imulti][iat][ipt][iht] += aeh*aeh;
		numer_errl[ifile][imulti][iat][ipt][iht] += ael*ael;
		denom[ifile][imulti][iat][ipt][iht]      += b;
		denom_errh[ifile][imulti][iat][ipt][iht] += beh*beh; 
		denom_errl[ifile][imulti][iat][ipt][iht] += bel*bel; 
	      
	      } 
	    }
	    
	    if (above) { delete above; }
	    if (below) { delete below; }
	    
	  } // Loop through pT bins
	} // Loop through aT bins
      } // Loop through multiplicity bins
      
      // Close file that was opened
      if (file) {
	file->cd();
	file->Close();
	delete file;
      } 
      
    } // Loop through files (jfile)
    
    for ( uint imulti = 0; imulti < nmulti; ++imulti ) {
      for ( uint iat = 0; iat < nat; ++iat ) {
	for ( uint ipt = 0; ipt < npt; ++ipt ) {
	  for ( uint iht = 0; iht < nht; ++iht ) {
	    
	    double a   = numer[ifile][imulti][iat][ipt][iht];
	    double aeh = sqrt(numer_errh[ifile][imulti][iat][ipt][iht]);
	    double ael = sqrt(numer_errl[ifile][imulti][iat][ipt][iht]);
	    numer_errh[ifile][imulti][iat][ipt][iht] = aeh;
	    numer_errl[ifile][imulti][iat][ipt][iht] = ael;
	      
	    double b   = denom[ifile][imulti][iat][ipt][iht];
	    double beh = sqrt(denom_errh[ifile][imulti][iat][ipt][iht]);
	    double bel = sqrt(denom_errl[ifile][imulti][iat][ipt][iht]);
	    denom_errh[ifile][imulti][iat][ipt][iht] = beh;
	    denom_errl[ifile][imulti][iat][ipt][iht] = bel;
	      
	    double r  = b > 0. ? a/b : 0.;
	    double rerrl = 0.;
	    double rerrh = 0.;
	    calcErr( r, rerrh, rerrl, a, aeh, ael, b, beh, bel );
	      
	    ratio[ifile][imulti][iat][ipt][iht] = r;
	    if ( r > 0. ) { length[ifile][imulti][iat][ipt]++; }
	    errh[ifile][imulti][iat][ipt][iht] = rerrh;
	    errl[ifile][imulti][iat][ipt][iht] = rerrl;
	      
	  }
	}
      }
    }
      
  } // Loop through files (ifile)
  
  bool print = false;
  if (print) {
    std::stringstream ss;
    for ( uint ifile = 0; ifile < nfile; ++ifile ) {
      for ( uint imulti = 0; imulti < nmulti; ++imulti ) {
	for ( uint iat = 0; iat < nat; ++iat ) {
	  for ( uint iht = 0; iht < nht; ++iht ) {
	    ss << " CHECK:"
	       << " ifile: " << ifile
	       << " multi: " << multi[imulti]
	       << " AlphaT: " << at[iat]
	       << " HT :" << ht[iht]
	       << " n: " << numer[ifile][imulti][iat][iht][iht]
	       << " n_errh: " << numer_errh[ifile][imulti][iat][iht][iht]
	       << " d: " << denom[ifile][imulti][iat][iht][iht]
	       << " d_errh: " << denom_errh[ifile][imulti][iat][iht][iht]
	       << " r: " << ratio[ifile][imulti][iat][iht][iht]
	       << " eh: " << errh[ifile][imulti][iat][iht][iht]
	       << " el: " << errl[ifile][imulti][iat][iht][iht];
	    if ( ratio[ifile][imulti][iat][iht][iht] <= 0. ) { ss << " (0,0)"; }
	    else { 
	      ss << " (" 
		 << dr((errh[ifile][imulti][iat][iht][iht]/ratio[ifile][imulti][iat][iht][iht]),2) 
		 << ","
		 << dr((errl[ifile][imulti][iat][iht][iht]/ratio[ifile][imulti][iat][iht][iht]),2)
		 << ")"; 
	    }
	    ss << std::endl;
	  }
	}
      }
    }
    std::cout << ss.str();
  }
  
}

// -----------------------------------------------------------------------------
// 
void calcRatio( const uint nfile, 
		const uint nmulti, 
		const uint nat, 
		const uint npt, 
		const uint nht, 
 		StringV his, 
 		StringVV files,
 		DoubleV lumis,
 		DoubleV weights,
 		IntV multi, 
 		IntV at, 
 		DoubleV pt, 
 		DoubleV ht, 
 		double ht_min, 
 		double ht_max,
 		double ht_step,
 		DoubleVVVVV& numer, 
 		DoubleVVVVV& numer_errh, 
 		DoubleVVVVV& numer_errl, 
 		DoubleVVVVV& denom, 
 		DoubleVVVVV& denom_errh, 
 		DoubleVVVVV& denom_errl, 
 		DoubleVVVVV& ratio, 
 		DoubleVVVVV& errh, 
 		DoubleVVVVV& errl,
 		IntVVVV& length,
 		int data_file,
 		double lumi,
		bool efficiency = false
		) {

  calcRatio( nfile, nmulti, nat, npt, nht, 
	     his, files, lumis, weights, multi, at, pt, ht, 
	     (int)nht, ht_min, ht_max, 
	     numer, numer_errh, numer_errl, 
	     denom, denom_errh, denom_errl, 
	     ratio, errh, errl, length, 
	     data_file, lumi,
	     efficiency );
  
}
