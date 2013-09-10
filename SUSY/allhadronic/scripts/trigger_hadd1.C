#include "./common/style.C"
#include <string.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include <map>
#include <iostream>
#include <sstream>
#include <string>

void trigger_hadd() {

  // in an interactive ROOT session, edit the file names
  // Target and FileList, then
  // root > .L hadd.C
  // root > hadd()

  int versus = 2;

  double lumi = 353.;
  
  //std::string file_name("Trigger_HT_Run2011A_PromptReco_v1");
  //std::string file_name("Trigger_HT_Run2011_promptReco_DCS");
  std::string file_name("Trigger_HT42_incomplete");
  std::string dir_name("Triggers");
  std::string his_name;
  if ( versus == 1 || versus == 2 ) his_name = "TriggersVsRunNumber";
  else if ( versus == 0 ) his_name = "TriggersVsPrescales";
  
  typedef std::map<int,float> Runs;
  typedef std::map<std::string,Runs> Triggers;
  Triggers trigs;
  std::vector<int> runs;
  
  for ( int ii = 0; ii < 1000; ++ii ) {
    
    // Construct file name
    std::stringstream ss;
    ss << "./" << file_name << "_" << ii << ".root";
    
    // Attempt to open file
    TFile* file = new TFile(ss.str().c_str());
    if ( file->IsZombie() || !(file->IsOpen()) ) { 
      if (file) { delete file; }
      std::cout << "Could not find file: " << ss.str() << std::endl; 
      continue; 
    } else {
      file->cd();
      std::cout << "Opened file: " << file->GetName() << std::endl; 
    }
    
    TDirectory* dir = (TDirectory*)file->Get(dir_name.c_str());
    if (dir) { std::cout << "Opened dir: " << dir->GetName() << std::endl; }
    else { 
      std::cout << "Could not find dir: " << dir_name << std::endl; 
      return; 
    }
    
    TH2D* his = (TH2D*)dir->Get(his_name.c_str());
    if (his) { std::cout << "Opened histo: " << his->GetName() << std::endl; }
    else { 
      std::cout << "Could not find histo: " << his_name << std::endl; 
      return; 
    }

    // Loop through x-bins 
    for ( int xbin = 0; xbin < his->GetNbinsX(); ++xbin ) {
      
      // Retrieve trigger name (and its prescale) from x-bin label
      std::string trigger = his->GetXaxis()->GetBinLabel(xbin+1);
      
      // Loop through y-bins (prescale values)
      for ( int ybin = 0; ybin < his->GetNbinsY(); ++ybin ) {
	
	// Check for non-zero content
	double weight = his->GetBinContent(xbin+1,ybin+1);
	if ( weight == 0. ) { continue; }
	
	// Retrieve prescale value (or run number from y-bin label)
	int yval = (int)his->GetYaxis()->GetBinLowEdge(ybin+1);
	if ( versus > 0 ) yval = atoi(his->GetYaxis()->GetBinLabel(ybin+1));

	// Cache runs used	
	if ( versus > 0 && weight > 0 ) { 
	  if ( find( runs.begin(), runs.end(), yval ) == runs.end() ) {
	    runs.push_back( yval );
	  }
	}
	
	// Populate map with histogram contents
	Triggers::const_iterator iter = trigs.find( trigger );
	if ( iter == trigs.end() ) {
	  // If trigger name not found
	  trigs[trigger][yval] = weight;
	} else {
	  if ( iter->second.find(yval) == iter->second.end() ) {
	    // If trigger name found but prescale not found
	    trigs[trigger][yval] = weight;
	  } else {
	    // If both trigger name and prescale found
	    if ( versus != 2 ) trigs[trigger][yval] += weight;
	  }
	}
	
      } // y-bins
    } // x-bins

    if (file) {
      file->Close();
      delete file;
    }

  } // files

  // Open new file to save triggers and prescales
  std::stringstream sss;
  sss << file_name << ".root";
  TFile* file = new TFile(sss.str().c_str(),"RECREATE");
  if ( file->IsZombie() || !(file->IsOpen()) ) { 
    if (file) { delete file; }
    std::cout << "Could not create file: " << sss.str() << std::endl; 
    return; 
  } else {
    file->cd();
    std::cout << "Opened file: " << file->GetName() << std::endl; 
  }

  // Create new dir
  TDirectory* dir = file->mkdir("Triggers");
  dir->cd();

  // Create histogram
  TH2D* his = new TH2D( "Triggers", 
			( versus > 0 ? ";;Run number" : ";;Prescales" ), 
			trigs.size(), 0., trigs.size(),
			10000, 0., 10000. );
  if ( versus > 0 ) his->GetYaxis()->Set( runs.size(), 0., runs.size() );
  his->Sumw2();


  // Cache info
  std::string last_trigger = "None";
  std::string last_version = "None";
  std::vector<int> indices;
  
  // Print some debug  
  std::stringstream ss;
  ss << " List of " << trigs.size()
     << " triggers that fired"
     << std::endl;
  Triggers::const_iterator ii = trigs.begin();
  Triggers::const_iterator jj = trigs.end();
  for ( ; ii != jj; ++ii ) {
    ss << " Trigger \"" << ii->first << "\""
       << std::endl;
    
    std::string curr_trigger = ii->first.substr( 0, ii->first.find("_v") );
    std::string curr_version = ii->first.substr( ii->first.find("_v") );
    if ( last_trigger == "None" || 
	 ( curr_trigger == last_trigger && 
	   curr_version != last_version ) ) { 
      int tmp = ii->second.empty() ? 0 : ii->second.begin()->first;
      if ( tmp && find( indices.begin(), indices.end(), tmp ) == indices.end() ) { 
	indices.push_back(tmp); 
	std::cout << " run " << tmp << std::endl;
      }
    }
    last_trigger = curr_trigger;
    last_version = curr_version;
    
    Runs::const_iterator iii = ii->second.begin();
    Runs::const_iterator jjj = ii->second.end();
    for ( ; iii != jjj; ++iii ) {
      std::stringstream run; run << iii->first << std::endl;
      his->Fill(ii->first.substr(4).c_str(),run.str().c_str(),iii->second*1.);
      ss << ( versus > 0 ? "  Run number: " : "  Prescale: " )
	 << iii->first
	 << " Events: " << iii->second
	 << std::endl;
    }    
  }
  //std::cout << ss.str() << std::endl;

  for ( int index = 0; index < indices.size(); ++index ) {
    std::cout << "found "<< indices[index] << std::endl;
  }
  
  // Set labels 
  if ( versus > 0 ) his->GetYaxis()->LabelsOption("a");
  for ( int index = 2; index < his->GetNbinsY(); ++index ) {
    int tmp = atoi(his->GetYaxis()->GetBinLabel( index ));
    std::cout << " tmp " << tmp << std::endl; 
    if ( find( indices.begin(), indices.end(), tmp ) == indices.end() ) { 
      his->GetYaxis()->SetBinLabel( index, "" );
      std::cout << " match " << std::endl;
    } else {
      std::cout << " NO match " << std::endl;
    }
  }
  
  // Print Trigger list for PSet
//   std::stringstream ssss;
//   ssss << " Triggers = [" << std::endl;
//   for ( ii = trigs.begin(); ii != trigs.end(); ++ii ) {
//     ssss << "\"" << ii->first << "\"," << std::endl;
//   }
//   ssss << "],";
//   std::cout << ssss.str() << std::endl;

  // Draw with some cosmetics
  setTDRStyle();
  set_plot_style();
  int size = ( versus > 0 ? 900 : 600 );
  TCanvas* canvas = new TCanvas("Triggers","Triggers",size,size);
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.2);
  canvas->SetBottomMargin(0.3);
  his->Draw("COLZ");
  his->GetXaxis()->LabelsOption("v");
  his->GetZaxis()->SetRangeUser(1.,10000.);
  his->GetYaxis()->SetTitleOffset(2.);
  if ( versus > 0 ) his->GetYaxis()->SetNoExponent();
  
  double max = pow( 10., int( log10(his->GetMaximum())*10. + 0.5 )/10. );
  his->GetZaxis()->SetTitle("Prescale");
  his->GetZaxis()->SetTitleOffset(1.4);
  his->GetZaxis()->SetRangeUser(0.,max);
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(1);
  if ( versus == 0 ) canvas->SetLogy();
  if ( versus > 0 ) canvas->SetLogz();
  
  TLatex* prelim = new TLatex( 0.16, 0.96, "#scale[0.8]{CMS 2011 preliminary}" );
  prelim->SetTextSize(0.03);
  prelim->SetNDC();
  
  std::stringstream ss1; ss1 << "#int L dt = " << lumi << " pb^{-1}, #sqrt{s} = 7 TeV";
  TLatex* lumis = new TLatex( 0.8, 0.96, ss1.str().c_str() );
  lumis->SetTextAlign(32);
  lumis->SetTextSize(0.024);
  lumis->SetNDC();

  std::stringstream ss3; ss3 << "N_{events} = " << his->Integral();
  TLatex* stats2 = new TLatex( 0.75, 0.80, ss3.str().c_str() );
  stats2->SetTextAlign(32);
  stats2->SetTextSize(0.03);
  stats2->SetNDC();

  std::stringstream ss2; ss2 << "N_{triggers} = " << his->GetEntries();
  TLatex* stats1 = new TLatex( 0.75, 0.75, ss2.str().c_str() );
  stats1->SetTextAlign(32);
  stats1->SetTextSize(0.03);
  stats1->SetNDC();
  
  prelim->Draw("same");
  //if ( versus == 0 ) 
  lumis->Draw("same");
  if ( versus == 0 ) stats1->Draw("same");
  if ( versus == 0 ) stats2->Draw("same");
  canvas->Modified();
  canvas->SaveAs( std::string(file_name+".pdf").c_str() );
  canvas->SaveAs( std::string(file_name+".C").c_str() );
  canvas->Write();

//   // Write histo to file and close
//   his->Write();
//   if (file) {
//     file->Close();
//     delete file;
//   }
  
}
