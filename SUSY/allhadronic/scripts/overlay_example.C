#include "overlay.C"

// -----------------------------------------------------------------------------
/*
  Arguments to createPlot():
  - path to input files
  - canvas title
  - histogram title
  - histogram directory
  - rebin histograms?
  - normalize all histograms to unit area?
  - logarithmic y-scale
  - output file
  - min for y-axis (optiional)
  - max for y-axis (optiional)
*/  
int overlay_example() {

  setTDRStyle();

  double lumi = 15.1;

  // Path to input files
  TString path("../../../results/batch/101102_1/100_100_50/");

  // Output file
  TFile* output = new TFile( "Plots.root", "RECREATE" );
  if ( !output || output->IsZombie() ) { return -1; }

  // -------------------- Jet multiplicities --------------------
   
  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_QCD", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_qcd.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_QCD", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_Data", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_data.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_Data", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_1/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Scaled_SM", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (1) {
    std::vector<TString> files, names, titles;
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT250");
    titles.push_back("250 < H_{T} #leq 300");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT300");
    titles.push_back("300 < H_{T} #leq 350");
    files.push_back("../../../results/batch/101028e_2/100_100_50/QcdBkgdEst_sm.root");
    names.push_back("HtMultiplicity_HT350");
    titles.push_back("H_{T} > 350");
    createPlot( files, "Multiplicity_Fixed_SM", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }


  // ---------- diff ht bins ----------

  //std::string histo = "BabyCaloOverMeffAfterDeadEcal";
  //   std::string histo = "BabyPfOverMeffAfterDeadEcal";
  std::string histo = "BabyMhtOverMetAfterDeadEcal";
  //    std::string histo = "GenBabyOverMeffNoMet";
  //   std::string histo = "GenBabyMhtOverMetNoMet";
  
  std::string bin0 = "HT250";
   std::string bin1 = "HT300";
   std::string bin2 = "HT350";

   //    std::string bin0 = "Meff400";
   //    std::string bin1 = "Meff450";
   //    std::string bin2 = "Meff500";

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.5_all");
    names.push_back(histo+"_"+bin0+"_aT0.51_all");
    names.push_back(histo+"_"+bin0+"_aT0.52_all");
    names.push_back(histo+"_"+bin0+"_aT0.53_all");
    names.push_back(histo+"_"+bin0+"_aT0.54_all");
    names.push_back(histo+"_"+bin0+"_aT0.55_all");
    createPlot( path, histo+"_HT250", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin1+"_aT0.5_all");
    names.push_back(histo+"_"+bin1+"_aT0.51_all");
    names.push_back(histo+"_"+bin1+"_aT0.52_all");
    names.push_back(histo+"_"+bin1+"_aT0.53_all");
    names.push_back(histo+"_"+bin1+"_aT0.54_all");
    names.push_back(histo+"_"+bin1+"_aT0.55_all");
    createPlot( path, histo+"_HT300", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin2+"_aT0.5_all");
    names.push_back(histo+"_"+bin2+"_aT0.51_all");
    names.push_back(histo+"_"+bin2+"_aT0.52_all");
    names.push_back(histo+"_"+bin2+"_aT0.53_all");
    names.push_back(histo+"_"+bin2+"_aT0.54_all");
    names.push_back(histo+"_"+bin2+"_aT0.55_all");
    createPlot( path, histo+"_HT350", names, titles, "QcdBkgdEst", lumi, 2, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.51_all");
    names.push_back(histo+"_"+bin1+"_aT0.51_all");
    names.push_back(histo+"_"+bin2+"_aT0.51_all");
    createPlot( path, histo+"_Check0.51", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.52_all");
    names.push_back(histo+"_"+bin1+"_aT0.52_all");
    names.push_back(histo+"_"+bin2+"_aT0.52_all");
    createPlot( path, histo+"_Check0.52", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.53_all");
    names.push_back(histo+"_"+bin1+"_aT0.53_all");
    names.push_back(histo+"_"+bin2+"_aT0.53_all");
    createPlot( path, histo+"_Check0.53", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.54_all");
    names.push_back(histo+"_"+bin1+"_aT0.54_all");
    names.push_back(histo+"_"+bin2+"_aT0.54_all");
    createPlot( path, histo+"_Check0.54", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back(histo+"_"+bin0+"_aT0.55_all");
    names.push_back(histo+"_"+bin1+"_aT0.55_all");
    names.push_back(histo+"_"+bin2+"_aT0.55_all");
    createPlot( path, histo+"_Check0.55", names, titles, "QcdBkgdEst", lumi, 2, true, true, output );
  }


  // --------------------------------------------------------------------------------

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.5_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.51_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.52_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.53_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.54_all");
    names.push_back("GenBabyOverHtNoMet_HT250_aT0.55_all");
    createPlot( path, "GenBabyOverHtNoMet", names, titles, "QcdBkgdEst", lumi, 1, true, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    //names.push_back("GenHtMultiplicity_HT0");
    names.push_back("GenHtMultiplicity_HT250");
    names.push_back("GenHtMultiplicity_HT300");
    names.push_back("GenHtMultiplicity_HT350");
    createPlot( path, "GenHtMultiplicity", names, titles, "QcdBkgdEst", lumi, 1, false, true, output );
  }
  
  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_Ht0_aT0_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.5_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.51_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.52_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.53_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.54_all");
    names.push_back("GenBabyOverHt_Ht0_aT0.55_all");
    createPlot( path, "GenBabyOverHt_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_HT250_aT0_all");
    names.push_back("GenBabyOverHt_HT350_aT0_all");
    names.push_back("GenBabyOverHt_HT250_aT0.5_all");
    names.push_back("GenBabyOverHt_HT350_aT0.5_all");
    names.push_back("GenBabyOverHt_HT250_aT0.51_all");
    names.push_back("GenBabyOverHt_HT350_aT0.51_all");
    names.push_back("GenBabyOverHt_HT250_aT0.52_all");
    names.push_back("GenBabyOverHt_HT350_aT0.52_all");
    createPlot( path, "GenBabyOverHt_Vs_Ht_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0_all");
    names.push_back("BabyOverHt_Ht0_aT0.5_all");
    names.push_back("BabyOverHt_Ht0_aT0.51_all");
    names.push_back("BabyOverHt_Ht0_aT0.52_all");
    names.push_back("BabyOverHt_Ht0_aT0.53_all");
    names.push_back("BabyOverHt_Ht0_aT0.54_all");
    names.push_back("BabyOverHt_Ht0_aT0.55_all");
    createPlot( path, "BabyOverHt_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_HT250_aT0_all");
    names.push_back("BabyOverHt_HT350_aT0_all");
    names.push_back("BabyOverHt_HT250_aT0.5_all");
    names.push_back("BabyOverHt_HT350_aT0.5_all");
    names.push_back("BabyOverHt_HT250_aT0.51_all");
    names.push_back("BabyOverHt_HT350_aT0.51_all");
    names.push_back("BabyOverHt_HT250_aT0.52_all");
    names.push_back("BabyOverHt_HT350_aT0.52_all");
    createPlot( path, "BabyOverHt_Vs_Ht_VsAlphaT", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0_all");
    names.push_back("BabyOverHt_HT250_aT0_all");
    names.push_back("BabyOverHt_HT300_aT0_all");
    names.push_back("BabyOverHt_HT350_aT0_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0.51_all");
    names.push_back("BabyOverHt_HT250_aT0.51_all");
    names.push_back("BabyOverHt_HT300_aT0.51_all");
    names.push_back("BabyOverHt_HT350_aT0.51_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverHt_Ht0_aT0.52_all");
    names.push_back("BabyOverHt_HT250_aT0.52_all");
    names.push_back("BabyOverHt_HT300_aT0.52_all");
    names.push_back("BabyOverHt_HT350_aT0.52_all");
    createPlot( path, "BabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0_all");
    names.push_back("BabyOverMeff_Meff400_aT0_all");
    names.push_back("BabyOverMeff_Meff450_aT0_all");
    names.push_back("BabyOverMeff_Meff500_aT0_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0.51_all");
    names.push_back("BabyOverMeff_Meff400_aT0.51_all");
    names.push_back("BabyOverMeff_Meff450_aT0.51_all");
    names.push_back("BabyOverMeff_Meff500_aT0.51_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("BabyOverMeff_Meff0_aT0.52_all");
    names.push_back("BabyOverMeff_Meff400_aT0.52_all");
    names.push_back("BabyOverMeff_Meff450_aT0.52_all");
    names.push_back("BabyOverMeff_Meff500_aT0.52_all");
    createPlot( path, "BabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverHt_Ht0_aT0_all");
    names.push_back("GenBabyOverHt_HT250_aT0_all");
    names.push_back("GenBabyOverHt_HT300_aT0_all");
    names.push_back("GenBabyOverHt_HT350_aT0_all");
    createPlot( path, "GenBabyJets_Ht", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  if (0) {
    std::vector<TString> files, names, titles;
    names.push_back("GenBabyOverMeff_Meff0_aT0_all");
    names.push_back("GenBabyOverMeff_Meff400_aT0_all");
    names.push_back("GenBabyOverMeff_Meff450_aT0_all");
    names.push_back("GenBabyOverMeff_Meff500_aT0_all");
    createPlot( path, "GenBabyJets_Meff", names, titles, "QcdBkgdEst", lumi, 5, false, true, output );
  }

  output->Write();
  output->Close();
  delete output; 

  return 0;

}


