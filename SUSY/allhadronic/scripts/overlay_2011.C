#include "common/overlay.C"

// -----------------------------------------------------------------------------
int overlay_2011() {
  
  setTDRStyle();
  
  // Output file
  TFile* output_file = new TFile( "Plots.root", "RECREATE" );
  if ( !output_file || output_file->IsZombie() ) { return -1; }

  // -----------------------------------------------------------------------------
  // HT distributions
   
  if (1) {
    
    std::string file("../../results/v60/Ratio__data.root");

    std::vector<std::string> files, names, titles;
    std::vector<int> marker_style; 
    std::vector<int> marker_colour; 
    std::vector<float> marker_size; 
    std::vector<float> lumis; 

    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT275_aT0_2");
    titles.push_back("275 < H_{T} #leq 325");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT325_aT0_2");
    titles.push_back("325 < H_{T} #leq 375");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT375_aT0_2");
    titles.push_back("375 < H_{T} #leq 475");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT475_aT0_2");
    titles.push_back("475 < H_{T} #leq 575");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT575_aT0_2");
    titles.push_back("575 < H_{T} #leq 675");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT675_aT0_2");
    titles.push_back("675 < H_{T} #leq 775");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT775_aT0_2");
    titles.push_back("775 < H_{T} #leq 875");
    files.push_back(file);
    names.push_back("HtDistrPassTrackless_HT875_aT0_2");
    titles.push_back("875 < H_{T}");
    
    marker_style.resize(files.size(),21);
    for ( int ii = 0; ii < files.size(); ++ii ) { marker_colour.push_back(kRed-9+ii); }
    marker_size.resize(files.size(),1.0);
    lumis.resize(files.size(),100.);

//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT275_aT0_2");
//     titles.push_back("275 < H_{T} #leq 325");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT325_aT0_2");
//     titles.push_back("325 < H_{T} #leq 375");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT375_aT0_2");
//     titles.push_back("375 < H_{T} #leq 475");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT475_aT0_2");
//     titles.push_back("475 < H_{T} #leq 575");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT575_aT0_2");
//     titles.push_back("575 < H_{T} #leq 675");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT675_aT0_2");
//     titles.push_back("675 < H_{T} #leq 775");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT775_aT0_2");
//     titles.push_back("775 < H_{T} #leq 875");
//     files.push_back("../../results/v24/Ratio_sm_pythia.root");
//     names.push_back("HtDistrPassTrackless_HT875_aT0_2");
//     titles.push_back("875 < H_{T}");

//     marker_style.resize(files.size(),25);
//     for ( int ii = 0; ii < files.size(); ++ii ) { marker_colour.push_back(kBlue-9+ii); }
//     marker_size.resize(files.size(),1.0);
//     lumis.resize(files.size(),602.);

    createPlot( output_file,
		"Data_HtDistrPassTrackless_aT0", 
		files, 
		names, 
		titles, 
		marker_style, 
		marker_colour, 
		marker_size, 
		lumis,
		769., 
		1,     // rebin 
		false, // normalise
		true  // log
		);
    
  }
   
  // -----------------------------------------------------------------------------
  // Multiplicity plots
   
  if (0) {

    bool norm = true;
    bool log = false;
    
    std::string file("../../results/v40/Ratio__data.root");
    
    std::vector<std::string> files, names, titles;
    std::vector<int> marker_style; 
    std::vector<int> marker_colour; 
    std::vector<float> marker_size; 
    std::vector<float> lumis; 

    files.push_back(file);
    names.push_back("HtMultiplicity_HT275_aT0");
    titles.push_back("275 < H_{T} #leq 325");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT325_aT0");
    titles.push_back("325 < H_{T} #leq 375");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT375_aT0");
    //titles.push_back("375 < H_{T} #leq 475");
    titles.push_back("375 < H_{T}");

    for ( int ii = 0; ii < files.size(); ++ii ) { marker_style.push_back(24+ii); }
    for ( int ii = 0; ii < files.size(); ++ii ) { marker_colour.push_back(kBlue-9+ii*2); }
    marker_size.resize(files.size(),2.0);
    lumis.resize(files.size(),100.);
    
    createPlot( output_file,
		"HtMultiplicity_LowHT", 
		files, 
		names, 
		titles, 
		marker_style, 
		marker_colour, 
		marker_size, 
		lumis,
		602., 
		1,     // rebin 
		norm, // normalise
		log  // log
		);


    files.clear();
    names.clear();
    titles.clear();
    marker_style.clear();
    marker_size.clear();
    marker_colour.clear();

    files.push_back(file);
    names.push_back("HtMultiplicity_HT475_aT0");
    titles.push_back("475 < H_{T} #leq 575");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT575_aT0");
    titles.push_back("575 < H_{T} #leq 675");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT675_aT0");
    titles.push_back("675 < H_{T} #leq 775");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT775_aT0");
    titles.push_back("775 < H_{T} #leq 875");
    files.push_back(file);
    names.push_back("HtMultiplicity_HT875_aT0");
    titles.push_back("875 < H_{T}");

    for ( int ii = 0; ii < files.size(); ++ii ) { marker_style.push_back(24+ii); }
    for ( int ii = 0; ii < files.size(); ++ii ) { marker_colour.push_back(kBlue-9+ii*2); }
    marker_size.resize(files.size(),2.0);
    lumis.resize(files.size(),100.);
    
    createPlot( output_file,
		"HtMultiplicity_HighHT", 
		files, 
		names, 
		titles, 
		marker_style, 
		marker_colour, 
		marker_size, 
		lumis,
		602., 
		1,     // rebin 
		norm, // normalise
		log  // log
		);
    
  }

  // -----------------------------------------------------------------------------
  // Write to disk

  output_file->Write();
  output_file->Close();
  delete output_file; 

  return 0;

}


