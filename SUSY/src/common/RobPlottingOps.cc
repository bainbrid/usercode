#include "RobPlottingOps.hh"
#include "AlphaT.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "GenMatrixBin.hh"
#include "RobOps.hh"
#include "Jet.hh"
#include "KinSuite.hh"
#include "Lepton.hh"
#include "Math/VectorUtil.h"
#include "Photon.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "Types.hh"
#include <sstream>
#include <string>
#include <iomanip>

#include "JetData.hh"
#include "CommonOps.hh"

using namespace Operation;

// -----------------------------------------------------------------------------
//
RobPlottingOps::RobPlottingOps( const Utils::ParameterSet& ps ) :
  // Misc
  dirName_("Plots"),
  nMin_(1),
  nMax_(10),
  verbose_(false),
  // Dalitz
  dalitz_(false),
  hDalitzRhoARhoB_(),
  hDalitzRhoARhoB2_(),
  hDalitzRhoARhoB3_(),
  hDalitzRhoARhoBn_(),
  hDalitzMTj1Mjj_(),
  hDalitzMTj2Mjj_(),
  hDalitzMTj1MTjj_(),
  hDalitzMTj2MTjj_(),
  hDalitzMT2j1M2jj_(),
  hDalitzMT2j2M2jj_(),
  hDalitzMT2j1MT2jj_(),
  hDalitzMT2j2MT2jj_(),
  // AlphaT
  alphaT_(false),
  hAlphaT_(),
  hAlphaTVal_(),
  hAlphaTnJets_(),
  // PtHat
  ptHat_(false),
  hPtHat_(),
  // MET
  met_(false),
  hMetDiff1_(),
  hMetDiff2_(),
  hMetCC_(),
  hMetRaw_(),
  hMetCalo_(),
  // Event summary
  summary_(false),
  // Cross-cleaning
  cc_(false),
  // Kinematic plots
  kine_(false),
  hEta_(),
  hDEta_(),
  hDR_(),
  hHt_(),
  ht_(false),
  hGenHt_(),
  genHt_(false),
  hMht_(),
  mht_(false),
  hMeff_(),
  meff_(false),
  hGenMeff_(),
  genMeff_(false),
  hMt2_(),
  mt2_(false),
  // Jet response
  response_(false),
  // AlphaT ratio 
  ratio_(false),
  alphaTcut_(0.55),
  hHtEqPre_(),
  hHtEqPost_(),
  hGenHtEqPre_(),
  hGenHtEqPost_(),
  hHtGtPre_(),
  hHtGtPost_(),
  hGenHtGtPre_(),
  hGenHtGtPost_(),
  hHtLtPre_(),
  hHtLtPost_(),
  hGenHtLtPre_(),
  hGenHtLtPost_()
{

  if ( ps.Contains("DirName") ) dirName_ = ps.Get<std::string>("DirName");
  if ( ps.Contains("MinObjects") ) nMin_ = ps.Get<int>("MinObjects");
  if ( ps.Contains("MaxObjects") ) nMax_ = ps.Get<int>("MaxObjects");
  if ( ps.Contains("Verbose") ) verbose_ = ps.Get<bool>("Verbose");
  if ( ps.Contains("Dalitz") ) dalitz_ = ps.Get<bool>("Dalitz");
  if ( ps.Contains("AlphaT") ) alphaT_ = ps.Get<bool>("AlphaT");
  if ( ps.Contains("PtHat") ) ptHat_ = ps.Get<bool>("PtHat");
  if ( ps.Contains("MET") ) met_ = ps.Get<bool>("MET");
  if ( ps.Contains("Summary") ) summary_ = ps.Get<bool>("Summary");
  if ( ps.Contains("CC") ) cc_ = ps.Get<bool>("CC");
  if ( ps.Contains("Kine") ) kine_ = ps.Get<bool>("Kine"); 
  if ( ps.Contains("HT") ) ht_ = ps.Get<bool>("HT"); 
  if ( ps.Contains("GenHT") ) genHt_ = ps.Get<bool>("GenHT"); 
  if ( ps.Contains("MHT") ) mht_ = ps.Get<bool>("MHT"); 
  if ( ps.Contains("Meff") ) meff_ = ps.Get<bool>("Meff"); 
  if ( ps.Contains("GenMeff") ) genMeff_ = ps.Get<bool>("GenMeff"); 
  if ( ps.Contains("Mt2") ) mt2_ = ps.Get<bool>("Mt2"); 
  if ( ps.Contains("Response") ) response_ = ps.Get<bool>("Response");
  if ( ps.Contains("Ratio") ) ratio_ = ps.Get<bool>("Ratio"); 
  if ( ps.Contains("AlphaTcut") ) alphaTcut_ = ps.Get<double>("AlphaTcut"); 
  
}

// -----------------------------------------------------------------------------
//
RobPlottingOps::~RobPlottingOps() {}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::Start( Event::Data& ev ) {
  initDir( ev.OutputFile(), dirName_.c_str() );
  BookHistos();
}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::BookHistos() {
  if ( dalitz_ )   { dalitz(); }
  if ( alphaT_ )   { alphaT(); }
  if ( ptHat_ )    { ptHat(); }
  if ( met_ )      { met(); }
  if ( kine_ )     { kine(); }
  if ( ht_ )       { ht(); }
  if ( genHt_ )    { genHt(); }
  if ( mht_ )      { mht(); }
  if ( meff_ )     { meff(); }
  if ( genMeff_ )  { genMeff(); }
  if ( mt2_ )      { mt2(); }
  if ( response_ ) { response(); }
  if ( ratio_ )    { ratio(); }
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::Process( Event::Data& ev ) {
  if ( summary_ )  { summary(ev); }
  if ( dalitz_ )   { dalitz(ev); }
  if ( alphaT_ )   { alphaT(ev); }
  if ( ptHat_ )    { ptHat(ev); }
  if ( met_ )      { met(ev); }
  if ( cc_ )       { cc(ev); }
  if ( kine_ )     { kine(ev); }
  if ( ht_ )       { ht(ev); }
  if ( genHt_ )    { genHt(ev); }
  if ( mht_ )      { mht(ev); }
  if ( meff_ )     { meff(ev); }
  if ( genMeff_ )  { genMeff(ev); }
  if ( mt2_ )      { mt2(ev); }
  if ( response_ ) { response(ev); }
  if ( ratio_ )    { ratio(ev); }
  return true;
}

// -----------------------------------------------------------------------------
//
std::ostream& RobPlottingOps::Description( std::ostream& ostrm ) {
  ostrm << "Rob plots ";
  ostrm << "(bins " << nMin_ << " to " << nMax_ << ") ";
  return ostrm;
}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::dalitz() {
  
  BookHistArray( hDalitzRhoARhoB_, 
		 "DalitzRho", 
		 ";#rho_{A};#rho_{B}",
		 50, 0., 0.5, 
		 50, 0., 0.5,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzRhoARhoB2_, 
		 "DalitzRho2", 
		 ";#rho_{A};#rho_{B}",
		 50, 0., 0.5, 
		 50, 0., 0.5,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzRhoARhoB3_, 
		 "DalitzRho3", 
		 ";#rho_{A};#rho_{B}",
		 50, 0., 0.5, 
		 50, 0., 0.5,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzRhoARhoBn_, 
		 "DalitzRhoN", 
		 ";#rho_{A};#rho_{B}",
		 50, 0., 0.5, 
		 50, 0., 0.5,
		 nMax_+1, 0, 1, true );
  
  BookHistArray( hDalitzMTj1Mjj_,
		 "hDalitzMTj1Mjj",
		 ";M_{T}(j_{1},#slash{E}_{T});M(j_{1},j_{2})", 
		 100, 0., 1000., 
		 100, 0., 1000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMTj2Mjj_,
		 "hDalitzMTj2Mjj",
		 ";M_{T}(j_{2},#slash{E}_{T});M(j_{1},j_{2})", 
		 100, 0., 1000., 
		 100, 0., 1000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMTj1MTjj_,
		 "hDalitzMTj1MTjj",
		 ";M_{T}(j_{1},#slash{E}_{T});M_{T}(j_{1},j_{2})", 
		 100, 0., 1000., 
		 100, 0., 1000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMTj2MTjj_,
		 "hDalitzMTj2MTjj",
		 ";M_{T}(j_{2},#slash{E}_{T});M_{T}(j_{1},j_{2})", 
		 100, 0., 1000., 
		 100, 0., 1000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMT2j1M2jj_,
		 "hDalitzMT2j1M2jj",
		 ";M_{T}^{2}(j_{1},#slash{E}_{T});M^{2}(j_{1},j_{2})", 
		 100, 0., 1000000., 
		 100, 0., 1000000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMT2j2M2jj_,
		 "hDalitzMT2j2M2jj",
		 ";M_{T}^{2}(j_{2},#slash{E}_{T});M^{2}(j_{1},j_{2})", 
		 100, 0., 1000000., 
		 100, 0., 1000000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMT2j1MT2jj_,
		 "hDalitzMT2j1MT2jj",
		 ";M_{T}^{2}(j_{1},#slash{E}_{T});M_{T}^{2}(j_{1},j_{2})", 
		 100, 0., 1000000., 
		 100, 0., 1000000.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hDalitzMT2j2MT2jj_,
		 "hDalitzMT2j2MT2jj",
		 ";M_{T}^{2}(j_{2},#slash{E}_{T});M_{T}^{2}(j_{1},j_{2})",
		 100, 0., 1000000., 
		 100, 0., 1000000.,
		 nMax_+1, 0, 1, true );

}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::ht() {
  BookHistArray( hHt_, 
		 "HT", 
		 ";HT [GeV];",
		 100, 0., 2000., 
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::genHt() {
  BookHistArray( hGenHt_, 
		 "GenHT", 
		 ";GenHT [GeV];",
		 100, 0., 2000., 
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::mht() {
  BookHistArray( hMht_, 
		 "MHT", 
		 ";MHT [GeV];",
		 100, 0., 1000., 
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::meff() {
  BookHistArray( hMeff_, 
		 "Meff", 
		 ";M_{eff} [GeV];",
		 100, 0., 2000.,
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::genMeff() {
  BookHistArray( hGenMeff_, 
		 "GenMeff", 
		 ";M_{eff} [GeV];",
		 100, 0., 2000.,
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::mt2() {
  BookHistArray( hMt2_, 
		 "Mt2", 
		 ";Mt2 [GeV];",
		 100, 0., 1000.,
		 nMax_+1, 0, 1, true );
}  

// -----------------------------------------------------------------------------
//
void RobPlottingOps::alphaT() {

  BookHistArray( hAlphaT_, 
		 "AlphaT", 
		 ";#alpha_{T}",
		 300, 0., 3.,
		 nMax_+1, 0, 1, true );

  BookHistArray( hAlphaTVal_, 
		 "AlphaTVal", 
		 ";old;new",
		 800, 0.3, 0.7, 
		   800, 0.3, 0.7,
		 nMax_+1, 0, 1, true );

  BookHistArray( hAlphaTnJets_, 
		 "AlphaTnJets", 
		 ";N_{jet};#alpha_{T}",
		 10, 0.5, 10.5, 
		 100, 0.0, 1.0,
		 1, 0, 1, true );

}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::ptHat() {
  BookHistArray( hPtHat_,
		 "PtHat",
		 "",
		 100, 0., 2500.,
		 1, 0, 1, true );
}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::met() {

  BookHistArray( hMetDiff1_, 
		 "CC1", 
		 "Ratio(RecoilMET,RawMET)",
		 1000, 0., 10.,
		 1, 0, 1, true );

  BookHistArray( hMetDiff2_, 
		 "CC2", 
		 "Ratio(RecoilMET,CaloMET)",
		 1000, 0., 10.,
		 1, 0, 1, true );

  BookHistArray( hMetCC_, 
		 "MET(CC)", 
		 "",
		 1000, 0., 1000.,
		 1, 0, 1, true );

  BookHistArray( hMetRaw_, 
		 "MET(Raw)", 
		 "",
		 1000, 0., 1000.,
		 1, 0, 1, true );

  BookHistArray( hMetCalo_, 
		 "MET(Calo)", 
		 "",
		 1000, 0., 1000.,
		 1, 0, 1, true );

//   BookHistArray( hEventDisplay_, 
// 		 "EventDisplay", 
// 		 "",
// 		 400, -2000., 2000.,
// 		 400, -2000., 2000.,
// 		 101, 0, 1, true );

}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::kine() {

  BookHistArray( hDEta_,
		 "DeltaEta",
		 ";#Delta#eta(j_{i},j_{i+1});", 
		 120, 
		 -12., 
		 12., 
		 nMax_+1, 0, 1, true);

  BookHistArray( hEta_,
		 "Eta",
		 ";#eta(j_{i});", 
		 240, 
		 -12., 
		 12., 
		 nMax_+1, 0, 1, true);
  
  BookHistArray( hDR_,
		 "DeltaR",
		 ";#DeltaR(j_{i},j_{i+1});", 
		 100, 
		 0., 
		 10., 
		 nMax_+1, 0, 1, true);

}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::response() {

  BookHistArray( hRespCorr_,
		 "ResponseCorr",
		 ";p_{T}^{gen};p_{T}^{corr}/p_{T}^{gen}",
		 100, 
		 0., 
		 2000., 
		 150,
		 0.,
		 1.5,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespEtaCorr_,
		 "ResponseEtaCorr",
		 ";#eta^{gen};p_{T}^{corr}/p_{T}^{gen}",
		 100, 
		 -5, 
		 5., 
		 150,
		 0.,
		 1.5,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespRaw_,
		 "ResponseRaw",
		 ";p_{T}^{gen};p_{T}^{raw}/p_{T}^{gen}",
		 100, 
		 0., 
		 2000., 
		 150,
		 0.,
		 1.5,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespEtaRaw_,
		 "ResponseEtaRaw",
		 ";#eta^{gen};p_{T}^{raw}/p_{T}^{gen}",
		 100, 
		 -5., 
		 5., 
		 150,
		 0.,
		 1.5,
		 nMax_+1, 0, 1, true);
  
  BookHistArray( hCorr_,
		 "Correction",
		 ";p_{T}^{corr}/p_{T}^{gen};f_{correction}",
		 150,
		 0.,
		 1.5,
		 100,
		 0.,
		 10.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hJetID_,
		 "JetID",
		 ";JetID",
		 4,
		 -0.5,
		 3.5,
		 150,
		 0.,
		 1.5,
		 1, 0, 1, true);
  
}

// -----------------------------------------------------------------------------
//
void RobPlottingOps::ratio() {

  std::stringstream ss; 
  ss << "AlphaT" << alphaTcut_;
  
  BookHistArray( hHtEqPre_, 
		 TString("HtEqPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hHtEqPost_, 
		 TString("HtEqPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtEqPre_, 
		 TString("GenHtEqPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtEqPost_, 
		 TString("GenHtEqPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );
  
  BookHistArray( hHtGtPre_, 
		 TString("HtGtPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hHtGtPost_, 
		 TString("HtGtPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtGtPre_, 
		 TString("GenHtGtPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtGtPost_, 
		 TString("GenHtGtPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );
  
  BookHistArray( hHtLtPre_, 
		 TString("HtLtPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hHtLtPost_, 
		 TString("HtLtPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtLtPre_, 
		 TString("GenHtLtPre"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

  BookHistArray( hGenHtLtPost_, 
		 TString("GenHtLtPost"+ss.str()), 
		 ";HT^{reco} [GeV];R(#alpha_{T})",
		 200, 0., 2000., 
		 nMax_+1, 0, 1, true );

}  

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::dalitz( Event::Data& ev ) {

  // Event weight
  Double_t weight = ev.GetEventWeight();

  UInt_t n = 0;
  double rho_a = -1.;
  double rho_b = -1.;

  if (0) { 

    // Common objects
    std::vector<Event::Jet const *> jets = ev.JD_CommonJets().accepted;
    n = jets.size();
    
    // MHT ( from "ev.JD_CommonJets().accepted")
    LorentzV mht = ev.CommonMHT();
    
    // Alpha_t variable
    std::vector<bool> pseudo;
    AlphaT()( jets, pseudo );
    if ( pseudo.size() != jets.size() ) { abort(); }
    
    // Variables for Dalitz plots
    LorentzV lv1(0.,0.,0.,0.);
    LorentzV lv2(0.,0.,0.,0.);
    if ( n == 2 ) {
      if ( jets[0] ) lv1 = *jets[0];
      if ( jets[1] ) lv2 = *jets[1];
    } else if ( n > 2 ) {
      for ( unsigned int i = 0; i < jets.size(); ++i ) {
	if ( jets[i] ) {
	  if ( pseudo[i] ) { lv1 += *jets[i]; }
	  else             { lv2 += *jets[i]; }
	}
      }
      if ( lv2.Pt() > lv1.Pt() ) { LorentzV tmp = lv1; lv1 = lv2; lv2 = tmp; }
    }
    
    rho_a = lv2.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );
    rho_b = lv1.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );

  } else if (0) { 

    std::vector<LorentzV> jets = ev.CommonObjects();
    n = jets.size();

    // MHT for all objects
    LorentzV mht = ev.CommonRecoilMET();
    
    // Alpha_t variable
    std::vector<bool> pseudo;
    AlphaT()( jets, pseudo );
    if ( pseudo.size() != jets.size() ) { abort(); }

    // Variables for Dalitz plots
    LorentzV lv1(0.,0.,0.,0.);
    LorentzV lv2(0.,0.,0.,0.);
    if ( n == 2 ) {
      lv1 = ev.CommonObjects()[0];
      lv2 = ev.CommonObjects()[1];
    } else if ( n > 2 ) {
      for ( unsigned int i = 0; i < ev.CommonObjects().size(); ++i ) {
	if ( pseudo[i] ) { lv1 += ev.CommonObjects()[i]; }
	else             { lv2 += ev.CommonObjects()[i]; }
      }
      if ( lv2.Pt() > lv1.Pt() ) { LorentzV tmp = lv1; lv1 = lv2; lv2 = tmp; }
    }

    rho_a = lv2.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );
    rho_b = lv1.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );

  } else {

  // Retrieve dalitz
    LorentzV mht;
    LorentzV lv1;
    LorentzV lv2;
    std::pair<double,double> rho = RobOps::dalitz( ev, mht, lv1, lv2, false );
    rho_a = rho.first;
    rho_b = rho.second;

  }

//   // Number of GenJets 
//   int ngen = 0;
//   ICF_LorentzVs::const_iterator igen = ev.genJetP4()->begin();
//   ICF_LorentzVs::const_iterator jgen = ev.genJetP4()->end();
//   for ( ; igen != jgen; ++igen ) {
//     if ( igen->Et() > 50 && fabs(igen->Eta()) < 3.0 ) { ngen++; }
//   }
  
  // DalitzRho plots
  if ( !hDalitzRhoARhoB_.empty() ) { 
    hDalitzRhoARhoB_[0]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen == 2 ) hDalitzRhoARhoB2_[0]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen == 3 ) hDalitzRhoARhoB3_[0]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen >= 4 ) hDalitzRhoARhoBn_[0]->Fill( rho_a, rho_b, weight ); 
  }
  if ( n >= nMin_ && n <= nMax_ && n < hDalitzRhoARhoB_.size() ) {
    hDalitzRhoARhoB_[n]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen == 2 ) hDalitzRhoARhoB2_[n]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen == 3 ) hDalitzRhoARhoB3_[n]->Fill( rho_a, rho_b, weight ); 
//     if ( ngen >= 4 ) hDalitzRhoARhoBn_[n]->Fill( rho_a, rho_b, weight ); 
  }

//   // DalitzMt plots
//   double m_jj  = ( lv1 + lv2 ).M();
//   double m2_jj  = ( lv1 + lv2 ).M2();
//   double mt_jj  = ( lv1 + lv2 ).Mt();
//   double mt2_jj  = ( lv1 + lv2 ).Mt2();

//   double mt_j1_mht = ( lv1 + mht ).Mt();
//   double mt2_j1_mht = ( lv1 + mht ).Mt2();

//   double mt_j2_mht = ( lv2 + mht ).Mt();
//   double mt2_j2_mht = ( lv2 + mht ).Mt2();


//   // M_{T}(j1,MHT) vs M(j1,j2)
//   if ( !hDalitzMTj1Mjj_.empty() ) { hDalitzMTj1Mjj_[0]->Fill( mt_j1_mht, m_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMTj1Mjj_.size() ) {
//     hDalitzMTj1Mjj_[n]->Fill( mt_j1_mht, m_jj, weight );
//   }

//   // M_{T}(j2,MHT) vs M(j1,j2)
//   if ( !hDalitzMTj2Mjj_.empty() ) { hDalitzMTj2Mjj_[0]->Fill( mt_j2_mht, m_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMTj2Mjj_.size() ) {
//     hDalitzMTj2Mjj_[n]->Fill( mt_j2_mht, m_jj, weight );
//   }

//   // M_{T}(j1,MHT) vs M_{T}(j1,j2)
//   if ( !hDalitzMTj1MTjj_.empty() ) { hDalitzMTj1MTjj_[0]->Fill( mt_j1_mht, mt_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMTj1MTjj_.size() ) {
//     hDalitzMTj1MTjj_[n]->Fill( mt_j1_mht, mt_jj, weight );
//   }

//   // M_{T}(j2,MHT) vs M_{T}(j1,j2)
//   if ( !hDalitzMTj2MTjj_.empty() ) { hDalitzMTj2MTjj_[0]->Fill( mt_j2_mht, mt_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMTj2MTjj_.size() ) {
//     hDalitzMTj2MTjj_[n]->Fill( mt_j2_mht, mt_jj, weight );
//   }

//   // M_{T}^{2}(j1,MHT) vs M^{2}(j1,j2)
//   if ( !hDalitzMT2j1M2jj_.empty() ) { hDalitzMT2j1M2jj_[0]->Fill( mt2_j1_mht, m2_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMT2j1M2jj_.size() ) {
//     hDalitzMT2j1M2jj_[n]->Fill( mt2_j1_mht, m2_jj, weight );
//   }

//   // M_{T}^{2}(j2,MHT) vs M^{2}(j1,j2)
//   if ( !hDalitzMT2j2M2jj_.empty() ) { hDalitzMT2j2M2jj_[0]->Fill( mt2_j2_mht, m2_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMT2j2M2jj_.size() ) {
//     hDalitzMT2j2M2jj_[n]->Fill( mt2_j2_mht, m2_jj, weight );
//   }

//   // M_{T}^{2}(j1,MHT) vs M_{T}^{2}(j1,j2)
//   if ( !hDalitzMT2j1MT2jj_.empty() ) { hDalitzMT2j1MT2jj_[0]->Fill( mt2_j1_mht, mt2_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMT2j1MT2jj_.size() ) {
//     hDalitzMT2j1MT2jj_[n]->Fill( mt2_j1_mht, mt2_jj, weight );
//   }

//   // M_{T}^{2}(j2,MHT) vs M_{T}^{2}(j1,j2)
//   if ( !hDalitzMT2j2MT2jj_.empty() ) { hDalitzMT2j2MT2jj_[0]->Fill( mt2_j2_mht, mt2_jj, weight ); }
//   if ( n >= nMin_ && n <= nMax_ && n < hDalitzMT2j2MT2jj_.size() ) {
//     hDalitzMT2j2MT2jj_[n]->Fill( mt2_j2_mht, mt2_jj, weight );
//   }
  
  return true;
  
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::alphaT( Event::Data& ev ) {

  // Event weight
  Double_t weight = ev.GetEventWeight();

  // Common objects
  UInt_t n = ev.CommonObjects().size();
  
  // Alpha_t variable
  std::vector<bool> pseudo;
  double alpha_t = AlphaT()( ev.CommonObjects(), pseudo );
  if ( pseudo.size() != n ) { abort(); }
  
  // AlphaT Plots
  if ( !hAlphaT_.empty() ) { hAlphaT_[0]->Fill( alpha_t, weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hAlphaT_.size() ) {
    hAlphaT_[n]->Fill( alpha_t, weight );
  }
  
  // Validation
  if ( n >= nMin_ && n <= nMax_ && n < hAlphaTVal_.size() && n <= Event::KS_MinDEt_MaxN() ) {
    hAlphaTVal_[0]->Fill( ev.CommonAlphaT(), alpha_t, weight ); 
    hAlphaTVal_[n]->Fill( ev.CommonAlphaT(), alpha_t, weight );
  }
  
  // AlphaT vs nJets
  if ( !hAlphaTnJets_.empty() ) { hAlphaTnJets_[0]->Fill( n*1., alpha_t, weight ); }

  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::ptHat( Event::Data& ev ) {

  // Event weight
  Double_t weight = ev.GetEventWeight();

  // PtHat plot
  if ( ptHat_ ) { hPtHat_[0]->Fill( ev.pthat(), weight ); }

  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::met( Event::Data& ev ) {
  
  // Event weight
  Double_t weight = ev.GetEventWeight();
  
  // Recoil MET from common objects
  LorentzV common_met = ev.CommonRecoilMET();

  // Recoil MET from raw ntuple objects
  LorentzV raw_met = RecoilMET(ev);
  
  // Calo MET from ntuple
  LorentzV calo_met = ev.CaloMET();
  
  double diff1 = ( fabs(common_met.Pt()) - fabs(raw_met.Pt()) ) / fabs(raw_met.Pt());

  double diff2 = ( fabs(common_met.Pt()) - fabs(calo_met.Pt()) ) / fabs(calo_met.Pt());

  // Histo
  if ( hMetDiff1_[0] ) { hMetDiff1_[0]->Fill( diff1, weight ); }
  if ( hMetDiff2_[0] ) { hMetDiff2_[0]->Fill( diff2, weight ); }
  
  // Histo
  if ( hMetCC_[0] ) { hMetCC_[0]->Fill( common_met.Pt(), weight ); }
  if ( hMetRaw_[0] ) { hMetRaw_[0]->Fill( raw_met.Pt(), weight ); }
  if ( hMetCalo_[0] ) { hMetCalo_[0]->Fill( calo_met.Pt(), weight ); }

  return true;
  
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::kine( Event::Data& ev ) {
  
  // Event weight
  Double_t weight = ev.GetEventWeight();
  
  // Common objects
  UInt_t n = ev.CommonObjects().size();
						
  // Loop through common objects 
  for ( unsigned int ii = 0; ii < n; ++ii ) {
    
    // Eta
    if ( !hEta_.empty() ) { hEta_[0]->Fill( ev.CommonObjects()[ii].Eta(), weight ); }
    if ( n >= nMin_ && n <= nMax_ && n < hEta_.size() ) {
      hEta_[ii]->Fill( ev.CommonObjects()[ii].Eta(), weight ); 
    }
    
    if ( ( n - ii ) > 1 ) {
      double deta = ev.CommonObjects()[ii].Eta() - ev.CommonObjects()[ii+1].Eta();
      double dphi = ROOT::Math::VectorUtil::DeltaPhi( ev.CommonObjects()[ii], ev.CommonObjects()[ii+1] );
      double dr = sqrt( deta * deta + dphi *dphi );
      // DeltaEta
      if ( !hDEta_.empty() ) { hDEta_[0]->Fill( deta, weight ); }
      if ( n >= nMin_ && n <= nMax_ && n < hDEta_.size() ) {
	hDEta_[ii]->Fill( deta, weight ); 
      }
      // DeltaR
      if ( !hDR_.empty() ) { hDR_[0]->Fill( dr, weight ); }
      if ( n >= nMin_ && n <= nMax_ && n < hDR_.size() ) {
	hDR_[ii]->Fill( dr, weight ); 
      }
    }

  }

  return true; 

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::response( Event::Data& ev ) {
  
  // Event weight
  Double_t weight = ev.GetEventWeight();
  
  std::stringstream ss;
  if ( verbose_ ) { 
    ss << "  Found " << ev.JD_CommonJets().accepted.size() << " jets!"
       << std::endl;
  }

  bool fail = false;
  unsigned int n = 0;
  std::vector<Event::Jet const*>::const_iterator ijet = ev.JD_CommonJets().accepted.begin();
  std::vector<Event::Jet const*>::const_iterator jjet = ev.JD_CommonJets().accepted.end();
  for ( ; ijet != jjet; ++ijet ) {
    
    n++;
    
    if ( !(*ijet) ) { continue; }
    
    // Index of original jet, as read from ntuple
    int index = (*ijet)->GetIndex(); 
    
    // Retrieve match index
    if ( index < int( ev.genJetMatchIndex()->size() ) ) {
      index = ev.genJetMatchIndex()->at(index);
    } else { index = -1; }
    
    // Check if index is valid for GenJet collection
    if ( index < 0 || index >= int(ev.genJetP4()->size()) ) { continue; }
    
    // Check if index is valid for CorrFactor collection
    if ( index < 0 || index >= int(ev.jetCorrFactor()->size()) ) { continue; }
    
//     std::cout << " int(ev.jetCorrFactor()->size()): " << int(ev.jetCorrFactor()->size())
// 	      << " int(ev.genJetP4()->size()): " << int(ev.genJetP4()->size())
// 	      << " (*ijet)->GetIndex(): " << (*ijet)->GetIndex()
// 	      << " index: " << index
// 	      << std::endl;

    // Retrieve LV for matched GenJet
    LorentzV gen_p4 = ev.genJetP4()->at(index);
    
    // Jet P4 and corrections
    LorentzV corr_p4 = **ijet;
    double corr = ev.jetCorrFactor()->at(index);
    LorentzV raw_p4 = corr != 0. ? ( corr_p4 / corr ) : LorentzV(0.,0.,0.,0.);
    
    double resp_corr = gen_p4.Pt() > 0. ? corr_p4.Pt() / gen_p4.Pt() : -1.;
    double resp_raw  = gen_p4.Pt() > 0. ? raw_p4.Pt()  / gen_p4.Pt() : -1.;

    if ( resp_corr < 0.5 ) { fail = true; }

    if ( !hJetID_.empty() ) { 
      if ( (*ijet)->GetTightId() ) { hJetID_[0]->Fill( 3., resp_corr, weight ); }
      else if ( (*ijet)->GetLooseId() ) { hJetID_[0]->Fill( 2., resp_corr, weight ); }
      else if ( (*ijet)->GetMinimalId() ) { hJetID_[0]->Fill( 1., resp_corr, weight ); }
      else { hJetID_[0]->Fill( 0., resp_corr, weight ); }
    }
    
    // Response plots

    if ( !hRespCorr_.empty() ) { hRespCorr_[0]->Fill( gen_p4.Pt(), resp_corr, weight ); }
    if ( index+1 < int(hRespCorr_.size()) ) { hRespCorr_[index+1]->Fill( gen_p4.Pt(), resp_corr, weight ); }
    
    if ( !hRespEtaCorr_.empty() ) { hRespEtaCorr_[0]->Fill( gen_p4.Eta(), resp_corr, weight ); }
    if ( index+1 < int(hRespEtaCorr_.size()) ) { hRespEtaCorr_[index+1]->Fill( gen_p4.Eta(), resp_corr, weight ); }

    if ( !hRespRaw_.empty() ) { hRespRaw_[0]->Fill( gen_p4.Pt(), resp_raw, weight ); }
    if ( index+1 < int(hRespRaw_.size()) ) { hRespRaw_[index+1]->Fill( gen_p4.Pt(), resp_raw, weight ); }

    if ( !hRespEtaRaw_.empty() ) { hRespEtaRaw_[0]->Fill( gen_p4.Eta(), resp_raw, weight ); }
    if ( index+1 < int(hRespEtaRaw_.size()) ) { hRespEtaRaw_[index+1]->Fill( gen_p4.Eta(), resp_raw, weight ); }
    
    if ( !hCorr_.empty() ) { hCorr_[0]->Fill( resp_corr, corr, weight ); }
    if ( index+1 < int(hCorr_.size()) ) { hCorr_[index+1]->Fill( resp_corr, corr, weight ); }
    
    if ( verbose_ ) {
      ss << " GenJetIndex  : " << std::setw(2) << index 
	 << std::endl
	 << " CorrJetIndex : " << std::setw(2) << int( ijet - ev.JD_CommonJets().accepted.begin() )
	 << std::endl
	 << " RawJetIndex  : " << std::setw(2) << (*ijet)->GetIndex() 
	 << std::endl
	 << " GenJet       : "
	 << "  Et: " << std::setw(8) << std::setprecision(4) << gen_p4.Et() 
	 << "  Eta: " << std::setw(8) << std::setprecision(4) << gen_p4.Eta() 
	 << "  Phi: " << std::setw(8) << std::setprecision(4) << gen_p4.Phi() 
	 << "  Px: " << std::setw(8) << std::setprecision(4) << gen_p4.Px() 
	 << "  Py: " << std::setw(8) << std::setprecision(4) << gen_p4.Py() 
	 << std::endl
	 << " CorrJet      : "
	 << "  Et: " << std::setw(8) << std::setprecision(4) << (*ijet)->Et() 
	 << "  Eta: " << std::setw(8) << std::setprecision(4) << (*ijet)->Eta() 
	 << "  Phi: " << std::setw(8) << std::setprecision(4) << (*ijet)->Phi() 
	 << "  Px: " << std::setw(8) << std::setprecision(4) << (*ijet)->Px() 
	 << "  Py: " << std::setw(8) << std::setprecision(4) << (*ijet)->Py() 
	 << std::endl
	 << " RawJet       : "
	 << "  Et: " << std::setw(8) << std::setprecision(4) << raw_p4.Et() 
	 << "  Eta: " << std::setw(8) << std::setprecision(4) << raw_p4.Eta() 
	 << "  Phi: " << std::setw(8) << std::setprecision(4) << raw_p4.Phi() 
	 << "  Px: " << std::setw(8) << std::setprecision(4) << raw_p4.Px() 
	 << "  Py: " << std::setw(8) << std::setprecision(4) << raw_p4.Py() 
	 << std::endl
	 << " CorrFactor   : " << std::setw(6) << std::setprecision(4) << corr
	 << std::endl
	 << " ResponseCorr : " << std::setw(6) << std::setprecision(4) << resp_corr
	 << std::endl;
// 	 << " ResponseRaw  : " << std::setw(6) << std::setprecision(4) << resp_raw
// 	 << std::endl;
    }

   }
  if ( verbose_ ) { std::cout << ss.str(); }

//   Compute::Jet::CommonJets RB_CommonJets(ev);
//   std::cout << " HERE1: has odd jet: " 
// 	    << RB_CommonJets().hasOdd
// 	    << std::endl;
  
//   Operation::OddJet RB_OddJet;
//   std::cout << " HERE2: pas odd jet veto " 
// 	    << RB_OddJet.Process(ev)
// 	    << std::endl;

  if ( verbose_ ) { std::cout << " RESPONSE FAIL: " << fail << std::endl; }
  return !fail;

//   return true;
  
}

// -----------------------------------------------------------------------------
//
LorentzV RobPlottingOps::RecoilMET( Event::Data& ev ) {
  
  LorentzV mData(0.,0.,0.,0.);

  std::vector<Event::Jet>::const_iterator ijet = ev.JD_Jets().begin();
  std::vector<Event::Jet>::const_iterator jjet = ev.JD_Jets().end();
  for ( ; ijet != jjet; ++ijet ) { mData -= (*ijet); }

  std::vector<Event::Lepton>::const_iterator imuon = ev.LD_Muons().begin();
  std::vector<Event::Lepton>::const_iterator jmuon = ev.LD_Muons().end();
  for ( ; imuon != jmuon; ++imuon ) { mData -= (*imuon); }
  
  std::vector<Event::Lepton>::const_iterator ielectron = ev.LD_Electrons().begin();
  std::vector<Event::Lepton>::const_iterator jelectron = ev.LD_Electrons().end();
  for ( ; ielectron != jelectron; ++ielectron ) { mData -= (*ielectron); }
    
  std::vector<Event::Photon>::const_iterator iphoton = ev.PD_Photons().begin();
  std::vector<Event::Photon>::const_iterator jphoton = ev.PD_Photons().end();
  for ( ; iphoton != jphoton; ++iphoton ) { mData -= (*iphoton); }
  
  return mData;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::ratio( Event::Data& ev ) {

  Double_t weight = ev.GetEventWeight();

  // RECO jets and variables
  std::vector<LorentzV> jets = ev.HadronicObjects(); 
  uint    n_reco = jets.size();
  double at_reco = AlphaT()( jets );
  double ht_reco = RobOps::ht( jets );

  // GEN jets
  std::vector<LorentzV> gen = RobOps::genJets( ev );
  uint    n_gen = gen.size();
  double at_gen = AlphaT()( gen );
  //double ht_gen = RobOps::ht( gen );
  
  for ( int ii = 0; ii < 3; ++ii ) {

    
    std::vector<TH1D*>* rPre;
    std::vector<TH1D*>* gPre;
    std::vector<TH1D*>* rPost;
    std::vector<TH1D*>* gPost;
    if ( ii == 0 && n_reco == 2 && n_gen == 2 ) { 
      rPre  = &hHtEqPre_; 
      rPost = &hHtEqPost_; 
      gPre  = &hGenHtEqPre_; 
      gPost = &hGenHtEqPost_; 
    } else if ( ii == 1 && n_reco == 2 && n_gen > 2 ) { 
      rPre  = &hHtGtPre_; 
      rPost = &hHtGtPost_; 
      gPre  = &hGenHtGtPre_; 
      gPost = &hGenHtGtPost_; 
    } else if ( ii == 2 && n_reco == 2 && n_gen < 2 ) { 
      rPre  = &hHtLtPre_; 
      rPost = &hHtLtPost_; 
      gPre  = &hGenHtLtPre_; 
      gPost = &hGenHtLtPost_; 
    } else { continue; }

    // RECO plots
    if ( !(*rPre).empty() && n_reco > 3 ) { 
      (*rPre)[0]->Fill( ht_reco, weight ); 
    }
    if ( n_reco >= nMin_ && n_reco <= nMax_ && n_reco < (*rPre).size() ) {
      (*rPre)[n_reco]->Fill( ht_reco, weight ); 
    }
    if ( at_reco >= alphaTcut_ ) {
      if ( !(*rPost).empty() && n_reco > 3 ) { 
	(*rPost)[0]->Fill( ht_reco, weight ); 
      }
      if ( n_reco >= nMin_ && n_reco <= nMax_ && n_reco < (*rPost).size() ) {
	(*rPost)[n_reco]->Fill( ht_reco, weight ); 
      }
    }
    
    // GEN plots
    if ( !(*gPre).empty() && n_reco > 3 ) { 
      (*gPre)[0]->Fill( ht_reco, weight ); 
    }
    if ( n_reco >= nMin_ && n_reco <= nMax_ && n_reco < (*gPre).size() ) {
      (*gPre)[n_reco]->Fill( ht_reco, weight ); 
    }
    if ( at_gen >= alphaTcut_ ) {
      if ( !(*gPost).empty() && n_reco > 3 ) { 
	(*gPost)[0]->Fill( ht_reco, weight ); 
      }
      if ( n_reco >= nMin_ && n_reco <= nMax_ && n_reco < (*gPost).size() ) {
	(*gPost)[n_reco]->Fill( ht_reco, weight ); 
      }
    }
    
  }

  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::ht( Event::Data& ev ) {
  Double_t weight = ev.GetEventWeight();
  uint n = 0;
  float ht = 0.;
  if (1) {
    std::vector<LorentzV> jets = ev.HadronicObjects(); 
    n = jets.size();
    ht = RobOps::ht( jets );
  } else {
    n = ev.CommonObjects().size();
    ht = ev.CommonHT();
  }
  if ( !hHt_.empty() ) { hHt_[0]->Fill( ht, weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hHt_.size() ) {
    hHt_[n]->Fill( ht, weight ); 
  }

  return true;
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::genHt( Event::Data& ev ) {
  
  std::vector<LorentzV> gen_jets = RobOps::genJets( ev );
  float gen_ht = RobOps::ht( gen_jets );
  
  Double_t weight = ev.GetEventWeight();
  uint n = ev.HadronicObjects().size();
  //uint n = gen_jets.size();
  if ( !hGenHt_.empty() ) { hGenHt_[0]->Fill( gen_ht, weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hGenHt_.size() ) {
    hGenHt_[n]->Fill( gen_ht, weight ); 
  }

  if ( 0 ) {

    {
      std::vector<LorentzV> all = *(ev.genJetP4());
      std::stringstream ss;
      ss << "GEN: " << all.size() << std::endl;
      std::vector<LorentzV>::const_iterator i = all.begin();
      std::vector<LorentzV>::const_iterator j = all.end();
      for ( ; i != j; ++i ) {
	ss << " pt/eta/phi: " << i->Pt() 
	   << "/" << i->Eta()
	   << "/" << i->Phi()
	   << std::endl;
      }    
      std::cout << ss.str();
    }
    
    {
      std::stringstream ss;
      ss << "SELECTED: " << gen_jets.size() << std::endl;
      std::vector<LorentzV>::const_iterator ij = gen_jets.begin();
      std::vector<LorentzV>::const_iterator jj = gen_jets.end();
      for ( ; ij != jj; ++ij ) {
	ss << " pt/eta/phi: " << ij->Pt() 
	   << "/" << ij->Eta()
	   << "/" << ij->Phi()
	   << std::endl;
      }
      std::cout << ss.str();
    }

    {
      std::vector<LorentzV> jets = ev.CommonObjects();
      std::stringstream ss;
      ss << "COMMON: " << jets.size() << std::endl;
      std::vector<LorentzV>::const_iterator ijet = jets.begin();
      std::vector<LorentzV>::const_iterator jjet = jets.end();
      for ( ; ijet != jjet; ++ijet ) {
	ss << " pt/eta/phi: " << ijet->Pt() 
	   << "/" << ijet->Eta()
	   << "/" << ijet->Phi()
	   << std::endl;
      }
      std::cout << ss.str();
    }
    
    {
      std::vector<LorentzV> hjets = ev.HadronicObjects();
      std::stringstream ss;
      ss << "HAD: " << hjets.size() << std::endl;
      std::vector<LorentzV>::const_iterator ijet = hjets.begin();
      std::vector<LorentzV>::const_iterator jjet = hjets.end();
      for ( ; ijet != jjet; ++ijet ) {
	ss << " pt/eta/phi: " << ijet->Pt() 
	   << "/" << ijet->Eta()
	   << "/" << ijet->Phi()
	   << std::endl;
      }
      std::cout << ss.str();
    }

    std::cout << std::endl;

  }

  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::mht( Event::Data& ev ) {
  Double_t weight = ev.GetEventWeight();
  uint n = ev.CommonObjects().size();
  if ( !hMht_.empty() ) { hMht_[0]->Fill( ev.CommonMHT().Pt(), weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hMht_.size() ) {
    hMht_[n]->Fill( ev.CommonMHT().Pt(), weight ); 
  }
  return true;
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::meff( Event::Data& ev ) {
  Double_t weight = ev.GetEventWeight();
  uint n = ev.CommonObjects().size();
  if ( !hMeff_.empty() ) { hMeff_[0]->Fill( ev.CommonMeff(), weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hMeff_.size() ) {
    hMeff_[n]->Fill( ev.CommonMeff(), weight ); 
  }
  return true;
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::genMeff( Event::Data& ev ) {
  
  std::vector<LorentzV> gen_jets = RobOps::genJets( ev );
  float gen_meff = RobOps::meff( gen_jets );
  
  Double_t weight = ev.GetEventWeight();
  uint n = ev.CommonObjects().size();
  if ( !hGenMeff_.empty() ) { hGenMeff_[0]->Fill( gen_meff, weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hGenMeff_.size() ) {
    hGenMeff_[n]->Fill( gen_meff, weight ); 
  }
  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::mt2( Event::Data& ev ) {
  Double_t weight = ev.GetEventWeight();
  uint n = ev.CommonObjects().size();
  if ( !hMt2_.empty() ) { hMt2_[0]->Fill( ev.CommonMT2(), weight ); }
  if ( n >= nMin_ && n <= nMax_ && n < hMt2_.size() ) {
    hMt2_[n]->Fill( ev.CommonMT2(), weight ); 
  }
  return true;
}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::summary( Event::Data& ev ) {

  std::stringstream ss;
  
  // Event info
  ss << "[RobPlottingOps::summary]" 
     << " Run: " << ev.RunNumber() 
     << " Lumi: " << ev.LumiSection() 
     << " Event: " << ev.EventNumber() 
     << " Weight: " << ev.GetEventWeight() << std::endl;

  // Number of raw objects
  ss << " Raw:   "
     << " Jets: " << ev.JD_Jets().size()
     << " Muons: " << ev.LD_Muons().size()
     << " Electrons: " << ev.LD_Electrons().size()
     << " Photons: " << ev.PD_Photons().size()  << std::endl;

  // Number of objects surviving CC 
  int nj = 0;
  {
    int n = 0;
    std::vector<Event::Jet const*>::const_iterator i = ev.CC().Jets.begin();
    std::vector<Event::Jet const*>::const_iterator j = ev.CC().Jets.end();
    for ( ; i != j; ++i ) { if ( *i ) { n++; } }
    nj = n;
  }

  int nm = 0;
  {
    int n = 0;
    std::vector<Event::Lepton const*>::const_iterator i = ev.CC().Muons.begin();
    std::vector<Event::Lepton const*>::const_iterator j = ev.CC().Muons.end();
    for ( ; i != j; ++i ) { if ( *i ) { n++; } }
    nm = n;
  }
  
  int ne = 0;
  {
    int n = 0;
    std::vector<Event::Lepton const*>::const_iterator i = ev.CC().Electrons.begin();
    std::vector<Event::Lepton const*>::const_iterator j = ev.CC().Electrons.end();
    for ( ; i != j; ++i ) { if ( *i ) { n++; } }
    ne = n;
  }
  
  int np = 0;
  {
    int n = 0;
    std::vector<Event::Photon const*>::const_iterator i = ev.CC().Photons.begin();
    std::vector<Event::Photon const*>::const_iterator j = ev.CC().Photons.end();
    for ( ; i != j; ++i ) { if ( *i ) { n++; } }
    np = n;
  }

  ss << " CC:    "
     << " Jets: " << nj
     << " Muons: " << nm
     << " Electrons: " << ne
     << " Photons: " << np << std::endl;
  
  // Number of obects accepted as being "Common"
  ss << " Common:"
     << " Jets: " << ev.JD_CommonJets().accepted.size()
     << " Muons: " << ev.LD_CommonMuons().accepted.size()
     << " Electrons: " << ev.LD_CommonElectrons().accepted.size()
     << " Photons: " << ev.PD_CommonPhotons().accepted.size()
     << std::endl;

  // Jet Pt
  if ( ev.JD_CommonJets().accepted.size() > 0 &&
       ev.JD_CommonJets().accepted[0] ) {
    ss << " LeadingJetPt: " << ev.JD_CommonJets().accepted[0]->Pt();
  }
  if ( ev.JD_CommonJets().accepted.size() > 1 &&
       ev.JD_CommonJets().accepted[1] ) {
    ss << " SecondJetPt: " << ev.JD_CommonJets().accepted[1]->Pt();
  }
  if ( ev.JD_CommonJets().accepted.size() > 2 &&
       ev.JD_CommonJets().accepted[2] ) {
    ss << " ThirdJetPt: " << ev.JD_CommonJets().accepted[2]->Pt();
  }
  ss << std::endl;

  // Odd obects and BadMuonInJet
  ss << " OddObjects:"
     << " Jets: " << ev.JD_CommonJets().hasOdd
     << " Muons: " << ev.LD_CommonMuons().hasOdd
     << " Electrons: " << ev.LD_CommonElectrons().hasOdd
     << " Photons: " << ev.PD_CommonPhotons().hasOdd
     << " BadMuonInJet: " << !badMuonInJet().Process(ev) << std::endl;
  
  // HT 
  double denominator = ( ev.CommonRecoilMET() + ev.JD_CommonJets().killedHT ).Pt();
  double missed_ht = denominator > 0. ? ev.CommonRecoilMET().Pt() / denominator : 0.;
  ss << " HT: " << ev.CommonHT()
     << " MHTratio: " << missed_ht << std::endl;
  
  // MET
  LorentzV gen_met(0.,0.,0.,0.);
  ICF_LorentzVs::const_iterator igen = ev.genJetP4()->begin();
  ICF_LorentzVs::const_iterator jgen = ev.genJetP4()->end();
  for ( ; igen != jgen; ++igen ) { gen_met -= (*igen); }
  ss << " CaloMET: " << ev.CaloMET().Pt()
     << " RecoilMET(GenObjects): " << gen_met.Pt()
     << " RecoilMET(RawObjects): " << RecoilMET(ev).Pt()
     << " RecoilMET(CommonObjects): " << ev.CommonRecoilMET().Pt() << std::endl;
  
  // AlphaT and biased DPhi
  ss << " AlphaT: " << AlphaT()( ev.CommonObjects() ) 
     << " BiasedDPhi: " << ev.BiasedCommonRecoilMETJetDPhi() << std::endl;
  
  // 
  
  // Output
  std::cout << ss.str();

  return true;

}

// -----------------------------------------------------------------------------
//
bool RobPlottingOps::cc( Event::Data& ev ) {
  
  // Recoil MET from common objects
  LorentzV common_met = ev.CommonRecoilMET();
  
  // Recoil MET from raw ntuple objects
  LorentzV raw_met = RecoilMET(ev);

  // Calo MET from ntuple
  LorentzV calo_met = ev.CaloMET();

  // Recoil MET from gen jets
  LorentzV gen_met(0.,0.,0.,0.);
  ICF_LorentzVs::const_iterator igen = ev.genJetP4()->begin();
  ICF_LorentzVs::const_iterator jgen = ev.genJetP4()->end();
  for ( ; igen != jgen; ++igen ) { gen_met -= (*igen); }

  // Checks
  double diff1 = ( fabs(common_met.Pt()) - fabs(raw_met.Pt()) ) / fabs(raw_met.Pt());
  double diff2 = ( fabs(common_met.Pt()) - fabs(calo_met.Pt()) ) / fabs(calo_met.Pt());
  bool fail = ( diff1 > 1. || diff2 > 1. ) ? true : false;
  
  // Debug
  std::stringstream ss;
  ss << "[RobPlottingOps::cc]"
     << "  FAIL: " << fail << std::endl
     << " HT                      : " << ev.CommonHT() << std::endl
     << " CaloMET                 : " << calo_met.Pt() << std::endl
     << " RecoilMHT(GenJets)      : " << gen_met.Pt() << std::endl
     << " RecoilMET(RawObjects)   : " << raw_met.Pt() << std::endl
     << " RecoilMET(CommonObjects): " << common_met.Pt() 
     << std::endl;
  //     std::cout << ss.str() << std::endl;
  //     ss.flush();
  
//   {
    
//     ss << " GEN OBJECTS FROM NTUPLE:" 
//        << std::endl;
    
//     GenMatrixBin gen(&ev);
//     //gen.PrintStableVis();

//     // Jets
//     ss << "  Found " << gen.the_GenJets.size() << " jets!"
//        << std::endl;
//     for ( unsigned int ijet = 0; ijet < gen.the_GenJets.size(); ++ijet ) {
//       //if ( gen.the_GenJets.at(ijet)->Et() < 50. || fabs( gen.the_GenJets.at(ijet)->Eta() ) > 3.0 ) { continue; }
//       ss << "  "
// 	 << " Index: " << std::setw(2) << ijet
// 	 << " Et: " << std::setw(8) << std::setprecision(4) << gen.the_GenJets.at(ijet)->Et() 
// 	 << " Eta: " << std::setw(8) << std::setprecision(4) << gen.the_GenJets.at(ijet)->Eta() 
// 	 << " Phi: " << std::setw(8) << std::setprecision(4) << gen.the_GenJets.at(ijet)->Phi() 
// 	 << " Px: " << std::setw(8) << std::setprecision(4) << gen.the_GenJets.at(ijet)->Px() 
// 	 << " Py: " << std::setw(8) << std::setprecision(4) << gen.the_GenJets.at(ijet)->Py() 
// 	 << "  --- "
// 	 << std::endl;
//     }

//     // Photons
//     ss << "  Found " << gen.the_GenPhot.size() << " photons!"
//        << std::endl;
//     for ( unsigned int ipho = 0; ipho < gen.the_GenPhot.size(); ++ipho ) {
//       //if ( gen.the_GenPhot.at(ijet)->Et() < 25. || fabs( gen.the_GenPhot.at(ijet)->Eta() ) > 2.5 ) { continue; }
//       ss << "  "
// 	 << " Index: " << std::setw(2) << ipho
// 	 << " Et: " << std::setw(8) << std::setprecision(4) << gen.the_GenPhot.at(ipho)->Et() 
// 	 << " Eta: " << std::setw(8) << std::setprecision(4) << gen.the_GenPhot.at(ipho)->Eta() 
// 	 << " Phi: " << std::setw(8) << std::setprecision(4) << gen.the_GenPhot.at(ipho)->Phi() 
// 	 << "  --- "
// 	 << std::endl;
//     }

//   }

  {
    
    ss << " GEN JETS FROM NTUPLE:" 
       << std::endl;
    
    // Jets
    ss << "  Found " << ev.genJetP4()->size() << " gen jets!"
       << std::endl;
    ICF_LorentzVs::const_iterator ijet = ev.genJetP4()->begin();
    ICF_LorentzVs::const_iterator jjet = ev.genJetP4()->end();
    for ( ; ijet != jjet; ++ijet ) {
      ss << "  "
	 << " Index: " << std::setw(2) << int( ijet - ev.genJetP4()->begin() )
	 << " Et: " << std::setw(8) << std::setprecision(4) << ijet->Et() 
	 << " Eta: " << std::setw(8) << std::setprecision(4) << ijet->Eta() 
	 << " Phi: " << std::setw(8) << std::setprecision(4) << ijet->Phi() 
	 << " Px: " << std::setw(8) << std::setprecision(4) << ijet->Px() 
	 << " Py: " << std::setw(8) << std::setprecision(4) << ijet->Py() 
	 << "  --- "
	 << std::endl;
    }
    
  }

  {
    
    ss << " RAW OBJECTS FROM NTUPLE:" 
       << std::endl;

    // Jets
    ss << "  Found " << ev.JD_Jets().size() << " jets!"
       << std::endl;
    std::vector<Event::Jet>::const_iterator ijet = ev.JD_Jets().begin();
    std::vector<Event::Jet>::const_iterator jjet = ev.JD_Jets().end();
    for ( ; ijet != jjet; ++ijet ) {

      // Match index
      int index = int( ijet - ev.JD_Jets().begin() );
      if ( index < int( ev.genJetMatchIndex()->size() ) ) {
	index = ev.genJetMatchIndex()->at(index);
      } else { index = -1; }
      
      ss << "  "
	 << " Index: " << std::setw(2) << ijet->GetIndex() 
	//<< " CC: " << std::setw(1) << ijet->WasItcc() 
	 << " Et: " << std::setw(8) << std::setprecision(4) << ijet->Et() 
	 << " Eta: " << std::setw(8) << std::setprecision(4) << ijet->Eta() 
	 << " Phi: " << std::setw(8) << std::setprecision(4) << ijet->Phi() 
	 << " Px: " << std::setw(8) << std::setprecision(4) << ijet->Px() 
	 << " Py: " << std::setw(8) << std::setprecision(4) << ijet->Py() 
	 << "  --- "
	 << " MatchIndex: " << std::setw(2) << index
	 << " MinimalId: " << std::setw(1) << ijet->GetMinimalId()
	 << " LooseId: " << std::setw(1) << ijet->GetLooseId()
	 << " TightId: " << std::setw(1) << ijet->GetTightId()
	 << " EMF: " << std::setw(8) << std::setprecision(4) << ijet->GetEmFrac() 
	 << std::endl;
    }

    // Muons
    ss << "  Found " << ev.LD_Muons().size() << " muons!"
       << std::endl;
    std::vector<Event::Lepton>::const_iterator imuon = ev.LD_Muons().begin();
    std::vector<Event::Lepton>::const_iterator jmuon = ev.LD_Muons().end();
    for ( ; imuon != jmuon; ++imuon ) {
      ss << "  "
     << " Index: " << std::setw(2) << imuon->GetIndex()
    //<< " CC: " << std::setw(1) << imuon->WasItcc()
     << " Et: " << std::setw(8) << std::setprecision(4) << imuon->Et()
     << " Eta: " << std::setw(8) << std::setprecision(4) << imuon->Eta()
     << " Phi: " << std::setw(8) << std::setprecision(4) << imuon->Phi()
     << "  --- "
     << " LooseId: " << std::setw(1) << imuon->GetLooseId()
     << " TightId: " << std::setw(1) << imuon->GetTightId()
     << " TrkIso: " << std::setw(8) << std::setprecision(4) << imuon->GetTrkIsolation()
     << " EcalIso: " << std::setw(8) << std::setprecision(4) << imuon->GetEcalIsolation()
     << " HcalIso: " << std::setw(8) << std::setprecision(4) << imuon->GetHcalIsolation()
     << " CombIso: " << std::setw(8) << std::setprecision(4) << imuon->GetCombIsolation()
     << std::endl;
    }

    // Electrons
    ss << "  Found " << ev.LD_Electrons().size() << " electrons!"
       << std::endl;
    std::vector<Event::Lepton>::const_iterator ielectron = ev.LD_Electrons().begin();
    std::vector<Event::Lepton>::const_iterator jelectron = ev.LD_Electrons().end();
    for ( ; ielectron != jelectron; ++ielectron ) {
      ss << "  "
     << " Index: " << std::setw(2) << ielectron->GetIndex()
    //<< " CC: " << std::setw(1) << ielectron->WasItcc()
     << " Et: " << std::setw(8) << std::setprecision(4) << ielectron->Et()
     << " Eta: " << std::setw(8) << std::setprecision(4) << ielectron->Eta()
     << " Phi: " << std::setw(8) << std::setprecision(4) << ielectron->Phi()
     << "  --- "
     << " LooseId: " << std::setw(1) << ielectron->GetLooseId()
     << " TightId: " << std::setw(1) << ielectron->GetTightId()
     << " TrkIso: " << std::setw(8) << std::setprecision(4) << ielectron->GetTrkIsolation()
     << " EcalIso: " << std::setw(8) << std::setprecision(4) << ielectron->GetEcalIsolation()
     << " HcalIso: " << std::setw(8) << std::setprecision(4) << ielectron->GetHcalIsolation()
     << " CombIso: " << std::setw(8) << std::setprecision(4) << ielectron->GetCombIsolation()
     << std::endl;
    }

    // Photons
    ss << "  Found " << ev.PD_Photons().size() << " photons!"
       << std::endl;
    std::vector<Event::Photon>::const_iterator iphoton = ev.PD_Photons().begin();
    std::vector<Event::Photon>::const_iterator jphoton = ev.PD_Photons().end();
    for ( ; iphoton != jphoton; ++iphoton ) {
      ss << "  "
     << " Index: " << std::setw(2) << iphoton->GetIndex()
    //<< " CC: " << std::setw(1) << iphoton->WasItcc()
     << " Et: " << std::setw(8) << std::setprecision(4) << iphoton->Et()
     << " Eta: " << std::setw(8) << std::setprecision(4) << iphoton->Eta()
     << " Phi: " << std::setw(8) << std::setprecision(4) << iphoton->Phi()
     << "  --- "
     << " LooseId: " << std::setw(1) << iphoton->IsItLoose()
     << " TightId: " << std::setw(1) << iphoton->IsItTight()
     << " TrkIso: " << std::setw(8) << std::setprecision(4) << iphoton->trkSumPtHollowConeDR04()
     << " EcalIso: " << std::setw(8) << std::setprecision(4) << iphoton->ecalRecHitSumEtConeDR04()
     << " HcalIso: " << std::setw(8) << std::setprecision(4) << iphoton->hcalTowerSumEtConeDR04()
     << " HadOverEm: " << std::setw(8) << std::setprecision(4) << iphoton->hadronicOverEm()
     << " r9: " << std::setw(8) << std::setprecision(4) << iphoton->r9()
     << " sigmaIetaIeta: " << std::setw(8) << std::setprecision(4) << iphoton->sigmaIetaIeta()
     << std::endl;
    }

  }
  
  {
    
    ss << " CC OBJECTS:"
       << std::endl;

    // Jets
    int njets = 0;
    std::stringstream jets;
    std::vector<Event::Jet const*>::const_iterator ijet = ev.CC().Jets.begin();
    std::vector<Event::Jet const*>::const_iterator jjet = ev.CC().Jets.end();
    for ( ; ijet != jjet; ++ijet ) {
      if ( !(*ijet) ) { continue; }
      jets << "  "
	   << " Index: " << std::setw(2) << (*ijet)->GetIndex() 
	//<< " CC: " << std::setw(1) << (*ijet)->WasItcc() 
	   << " Et: " << std::setw(8) << std::setprecision(4) << (*ijet)->Et() 
	   << " Eta: " << std::setw(8) << std::setprecision(4) << (*ijet)->Eta() 
	   << " Phi: " << std::setw(8) << std::setprecision(4) << (*ijet)->Phi() 
	   << " Px: " << std::setw(8) << std::setprecision(4) << (*ijet)->Px() 
	   << " Py: " << std::setw(8) << std::setprecision(4) << (*ijet)->Py() 
	   << "  --- "
	   << " EMF: " << std::setw(8) << std::setprecision(4) << (*ijet)->GetEmFrac() 
	   << std::endl;
      njets++;
    }
    ss << "  Found " << njets << " jets with details:" << std::endl;
    ss << jets.str();

    // Muons
    int nmuons = 0;
    std::stringstream muons;
    std::vector<Event::Lepton const*>::const_iterator imuon = ev.CC().Muons.begin();
    std::vector<Event::Lepton const*>::const_iterator jmuon = ev.CC().Muons.end();
    for ( ; imuon != jmuon; ++imuon ) {
      if ( !(*imuon) ) { continue; }
      nmuons++;
      muons << std::setw(2) << (*imuon)->GetIndex() << " ";
    }
    ss << "  Found " << nmuons << " muons with indices:";
    ss << muons.str() << std::endl;

    // Electrons
    int nelectrons = 0;
    std::stringstream electrons;
    std::vector<Event::Lepton const*>::const_iterator ielectron = ev.CC().Electrons.begin();
    std::vector<Event::Lepton const*>::const_iterator jelectron = ev.CC().Electrons.end();
    for ( ; ielectron != jelectron; ++ielectron ) {
      if ( !(*ielectron) ) { continue; }
      nelectrons++;
      electrons << std::setw(2) << (*ielectron)->GetIndex() << " ";
    }
    ss << "  Found " << nelectrons << " electrons with indices:";
    ss << electrons.str() << std::endl;

    // Photons
    int nphotons = 0;
    std::stringstream photons;
    std::vector<Event::Photon const*>::const_iterator iphoton = ev.CC().Photons.begin();
    std::vector<Event::Photon const*>::const_iterator jphoton = ev.CC().Photons.end();
    for ( ; iphoton != jphoton; ++iphoton ) {
      if ( !(*iphoton) ) { continue; }
      nphotons++;
      photons << std::setw(2) << (*iphoton)->GetIndex() << " ";
    }
    ss << "  Found " << nphotons << " photons with indices:";
    ss << photons.str() << std::endl;

  }

  {
    
    ss << " CC KILL MAP:"
       << std::endl;
    if ( ev.CC().KillMap.empty() ) {
      ss << "  Empty" << std::endl;
    } else {
      std::stringstream kills;
      std::vector<Compute::KillOp>::const_iterator ii = ev.CC().KillMap.begin();
      std::vector<Compute::KillOp>::const_iterator jj = ev.CC().KillMap.end();
      for ( ; ii != jj; ++ii ) { 
	kills << "  ";
	ii->terse(kills); 
	kills << std::endl; 
      }
      ss << kills.str();
    }

  }
  
  {

    ss << " ACCEPTED OBJECTS:"
       << std::endl;

    // Jets
    ss << "  Found " << ev.JD_CommonJets().accepted.size() << " jets with indices:";
    std::vector<Event::Jet const*>::const_iterator ijet = ev.JD_CommonJets().accepted.begin();
    std::vector<Event::Jet const*>::const_iterator jjet = ev.JD_CommonJets().accepted.end();
    for ( ; ijet != jjet; ++ijet ) {
      if ( !(*ijet) ) { continue; }
      ss << std::setw(2) << (*ijet)->GetIndex() << " ";
    }
    ss << std::endl;

    // Muons
    ss << "  Found " << ev.LD_CommonMuons().accepted.size() << " muons with indices: ";
    std::vector<Event::Lepton const*>::const_iterator imuon = ev.LD_CommonMuons().accepted.begin();
    std::vector<Event::Lepton const*>::const_iterator jmuon = ev.LD_CommonMuons().accepted.end();
    for ( ; imuon != jmuon; ++imuon ) {
      if ( !(*imuon) ) { continue; }
      ss << std::setw(2) << (*imuon)->GetIndex() << " ";
    }
    ss << std::endl;

    // Electrons
    ss << "  Found " << ev.LD_CommonElectrons().accepted.size() << " electrons with indices: ";
    std::vector<Event::Lepton const*>::const_iterator ielectron = ev.LD_CommonElectrons().accepted.begin();
    std::vector<Event::Lepton const*>::const_iterator jelectron = ev.LD_CommonElectrons().accepted.end();
    for ( ; ielectron != jelectron; ++ielectron ) {
      if ( !(*ielectron) ) { continue; }
      ss << std::setw(2) << (*ielectron)->GetIndex() << " ";
    }
    ss << std::endl;

    // Photons
    ss << "  Found " << ev.PD_CommonPhotons().accepted.size() << " photons with indices: ";
    std::vector<Event::Photon const*>::const_iterator iphoton = ev.PD_CommonPhotons().accepted.begin();
    std::vector<Event::Photon const*>::const_iterator jphoton = ev.PD_CommonPhotons().accepted.end();
    for ( ; iphoton != jphoton; ++iphoton ) {
      if ( !(*iphoton) ) { continue; }
      ss << std::setw(2) << (*iphoton)->GetIndex() << " ";
    }
    ss << std::endl;

  }

  {
    
    ss << " OBJECTS KILLED:"
       << std::endl;

    // Jets
    ss << "  Found " << ev.JD_CommonJets().killed.size() << " jets with indices:";
    std::vector<Event::Jet const*>::const_iterator ijet = ev.JD_CommonJets().killed.begin();
    std::vector<Event::Jet const*>::const_iterator jjet = ev.JD_CommonJets().killed.end();
    for ( ; ijet != jjet; ++ijet ) {
      if ( !(*ijet) ) { continue; }
      ss << " " /* << std::setw(2) */ << (*ijet)->GetIndex() ;
    }
    ss << std::endl;

    // Muons
    ss << "  Found " << ev.LD_CommonMuons().killed.size() << " muons with indices:";
    std::vector<Event::Lepton const*>::const_iterator imuon = ev.LD_CommonMuons().killed.begin();
    std::vector<Event::Lepton const*>::const_iterator jmuon = ev.LD_CommonMuons().killed.end();
    for ( ; imuon != jmuon; ++imuon ) {
      if ( !(*imuon) ) { continue; }
      ss << " " /* << std::setw(2) */ << (*imuon)->GetIndex();
    }
    ss << std::endl;

    // Electrons
    ss << "  Found " << ev.LD_CommonElectrons().killed.size() << " electrons with indices:";
    std::vector<Event::Lepton const*>::const_iterator ielectron = ev.LD_CommonElectrons().killed.begin();
    std::vector<Event::Lepton const*>::const_iterator jelectron = ev.LD_CommonElectrons().killed.end();
    for ( ; ielectron != jelectron; ++ielectron ) {
      if ( !(*ielectron) ) { continue; }
      ss << " " /* << std::setw(2) */ << (*ielectron)->GetIndex();
    }
    ss << std::endl;

    // Photons
    ss << "  Found " << ev.PD_CommonPhotons().killed.size() << " photons with indices:";
    std::vector<Event::Photon const*>::const_iterator iphoton = ev.PD_CommonPhotons().killed.begin();
    std::vector<Event::Photon const*>::const_iterator jphoton = ev.PD_CommonPhotons().killed.end();
    for ( ; iphoton != jphoton; ++iphoton ) {
      if ( !(*iphoton) ) { continue; }
      ss << " " /* << std::setw(2) */ << (*iphoton)->GetIndex();
    }
    ss << std::endl;

    std::cout << ss.str() << std::endl; 

  }

  return true;

}
