#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <vector>


void calcIntegral( TH1D* cum, bool greater_than ) {

  double val = 0.;
  double err = 0.;

  // Loop through bins 
  int nbins = cum->GetNbinsX();
  for ( int ibin = 0; ibin <= nbins+1; ibin++ ) {
    
    // Check direction
    int bin = ibin;
    if ( greater_than ) { bin = nbins+1 - ibin; }
    
    // Running totals
    val += cum->GetBinContent(bin);
    err += cum->GetBinError(bin) * cum->GetBinError(bin);
    
    // Set bin content and error
    cum->SetBinContent(bin,val);
    cum->SetBinError(bin,sqrt(err));
    
  }
  
  // Number of entries
  cum->SetEntries(nbins+1);

}


void eff() {

  gStyle->SetOptStat(1111111);
  
  std::vector<std::string> var;
  //var.push_back("MHTovHT");
  var.push_back("HadronicAlphaT");
  
  std::vector<std::string> his;
  his.push_back("eq2");
  his.push_back("ge3");
  //his.push_back("ge2");

  std::vector<std::string> pt;
  //pt.push_back("40_20");
  //pt.push_back("40_40");
  //pt.push_back("40_80");
  //pt.push_back("40_93");
   pt.push_back("30");
    pt.push_back("34");
    pt.push_back("38");
    pt.push_back("42");
//    pt.push_back("53");
//    pt.push_back("66");
//    pt.push_back("80");
//   pt.push_back("93");
  //    pt.push_back("50");
  //    pt.push_back("57");
  //    pt.push_back("64");
  //    pt.push_back("71");
  //    pt.push_back("50");
  //    pt.push_back("54");
  //    pt.push_back("59");
  //    pt.push_back("63");
   
   std::vector<std::string> ht;
//     ht.push_back("150");
//      ht.push_back("200");
//      ht.push_back("250");
//      ht.push_back("300");
   ht.push_back("350");
   ht.push_back("400");
   ht.push_back("450");
   ht.push_back("500");
   //    ht.push_back("550");
   //    ht.push_back("600");
   //    ht.push_back("650");
   //    ht.push_back("700");
  
  for ( int kk = 0; kk < pt.size(); ++kk ) {
    for ( int ll = 0; ll < ht.size(); ++ll ) {
      
      TFile* f =  new TFile(std::string("AK5Calo_LM1_"+pt[kk]+"_"+ht[ll]+".root").c_str(),"READ");
      if ( f && !f->IsZombie() ) { 
	std::cout << "Opened file: " << f->GetName() << std::endl; 
      } else { 
	std::cout << "Could not find file " << std::endl; 
	//return; 
      }
	  
      TDirectory* d1 = (f==0?0:(TDirectory*)f->Get(std::string("alphaTnumbers"+ht[ll]).c_str()));
      if (d1) { std::cout << "Opened dir: " << d1->GetName() << std::endl; }
      else { 
	std::cout << "Could not find dir 1 " << std::endl; 
	//return; 
      }
      
      TDirectory* d2 = (f==0?0:(TDirectory*)f->Get(std::string("alphaTnumbers"+ht[ll]+"NoTrigger").c_str()));
      if (d2) { std::cout << "Opened dir: " << d2->GetName() << std::endl; }
      else { 
	std::cout << "Could not find dir 2 " << std::endl; 
	//return; 
      }
      
      for ( int jj = 0; jj < his.size(); ++jj ) {
	for ( int ii = 0; ii < var.size(); ++ii ) {

	  std::cout << " ii " << ii 
		    << " jj " << jj
		    << " kk " << kk
		    << " ll " << ll 
		    << std::endl;

	  std::string str;
	  if      ( his[jj] == "eq2" ) { str = var[ii]+"_2"; }
	  else if ( his[jj] == "ge3" ) { str = var[ii]+"Geq3_all"; }
	  else                         { str = var[ii]+"_all"; }

	  TH1D* h1 = (d1==0?0:(TH1D*)d1->Get(str.c_str()));
	  if (h1) { std::cout << "Opened histo: " << h1->GetName() << std::endl; }
	  else { 
	    std::cout << "Could not find histo 1 " << std::endl; 
	    //continue; 
	  }
	  
	  TH1D* t1 = (h1==0?0:(TH1D*)h1->Clone());
	  if (t1) calcIntegral(t1,true);
	  //t1->ComputeIntegral();
	  //Double_t* i1 = t1->GetIntegral();
	  //t1->SetContent(i1);
	  
	  TH1D* h2 = (d2==0?0:(TH1D*)d2->Get(str.c_str()));
	  if (h2) { std::cout << "Opened histo: " << h2->GetName() << std::endl; }
	  else { 
	    std::cout << "Could not find histo 2 " << std::endl; 
	    //continue; 
	  }

	  TH1D* t2 = (h2==0?0:(TH1D*)h2->Clone());
	  if (t2) calcIntegral(t2,true);
	  //t2->ComputeIntegral();
	  //Double_t* i2 = t2->GetIntegral();
	  //t2->SetContent(i2);

	  TGraphAsymmErrors* gr1 = new TGraphAsymmErrors();
	  if ( h1 && h2 ) gr1->BayesDivide(h1,h2);
	  
	  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors();
	  if ( t1 && t2 ) gr2->BayesDivide(t1,t2);
	  
	  TCanvas* c1 = new TCanvas(std::string("eff_"+var[ii]+"_"+his[jj]+"_"+pt[kk]+"_"+ht[ll]).c_str(),
				    std::string("eff_"+var[ii]+"_"+his[jj]+"_"+pt[kk]+"_"+ht[ll]).c_str(),
				    600,600);
	  c1->Divide(2,3);
	  
	  c1->cd(1);
	  if ( h1 ) { 
	    h1->GetXaxis()->SetRangeUser(0.,1.);
	    h1->Draw();
	  }

	  c1->cd(2);
	  if ( t1 ) {
  	    t1->GetXaxis()->SetTitle("#alpha_{T} > cut value");
	    t1->GetXaxis()->SetRangeUser(0.,1.);
	    t1->Draw();
	  }

	  c1->cd(3);
	  if ( h2 ) { 
	    h2->GetXaxis()->SetRangeUser(0.,1.);
	    h2->Draw();
	  }

	  c1->cd(4);
	  if ( t2 ) {
  	    t2->GetXaxis()->SetTitle("#alpha_{T} > cut value");
	    t2->GetXaxis()->SetRangeUser(0.,1.);
	    t2->Draw();
	  }

	  c1->cd(5);
	  if ( h1 && h2 ) { 
  	    gr1->GetXaxis()->SetTitle("#alpha_{T}");
  	    gr1->GetYaxis()->SetTitle("Efficiency");
	    gr1->GetXaxis()->SetRangeUser(0.,1.);
	    gr1->Draw("alp");
	  }
	  
 	  c1->cd(6);
 	  if ( t1 && t2 ) {
  	    gr2->GetXaxis()->SetTitle("#alpha_{T} > cut value");
  	    gr2->GetYaxis()->SetTitle("Efficiency");
  	    gr2->GetXaxis()->SetRangeUser(0.4,1.0);
  	    gr2->GetYaxis()->SetRangeUser(0.9,1.);
 	    gr2->Draw("alp");
 	  }
	  
	  c1->Update();
	  c1->SaveAs(std::string("eff_"+var[ii]+"_"+his[jj]+"_"+pt[kk]+"_"+ht[ll]+".pdf").c_str());
	    
	} // ll
      } // kk
    } // jj
  } // ii

}
