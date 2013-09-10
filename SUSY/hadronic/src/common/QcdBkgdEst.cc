#include "QcdBkgdEst.hh"
#include "AlphaT.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "HadronicOps.hh"
#include "HadronicPlottingOps.hh"
#include "Jet.hh"
#include "JetData.hh"
#include "KinSuite.hh"
#include "Math/VectorUtil.h"
#include "BkgdEstPlottingOps.hh"
#include "TH1D.h"
#include "TH2D.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

/*
  - pthat
  - mht cleaning
  - aT/HT scatter, different binning?
  - truth!!
*/

using namespace Operation;

// -----------------------------------------------------------------------------
//
bool SortByPt( const LorentzV& first, const LorentzV& second ) {
  return ( first.Pt() > second.Pt() );
}

// -----------------------------------------------------------------------------
//
QcdBkgdEst::QcdBkgdEst( const Utils::ParameterSet& ps ) :
  // Misc
  dirName_("QcdBkgdEst"),
  nMin_(2),
  nMax_(10),
  verbose_(false),
  filter_(-1),
  noQcd_(false),
  noCleaningInDenom_(false),
  norm_(2.),
  scale_(true),
  scaleSignal_(false),
  scaleSingle_(true),
  maxLeadingJetEta_(0.),
  babyCut_(0.),
  addMuon_(false),
  // Dead ECAL
  deadEcalRegionDrCut_(0.),
  deadEcalRegionEtaCut_(0.),
  minJetPtCut_(0.),
  nBadCells_(0),
  useFile_(false),
  deadECALfile_(""),
  minBiasDeltaPhiCut_(0.),
  dead_(0),
  //RecHit cleaning
  rechit_(0),
  // Analysis defaults
  alphaT_(),
  aT_(0.55),
  ht_(0),
  pt1_(0.),
  pt2_(0.),
  pt3_(0.),
  // Binning
  inclusive_(true),
  htNbins_(0),
  htXlow_(0.),
  htXhigh_(0.),
  htBins_(),
  htPt1Bins_(),
  htPt2Bins_(),
  htPt3Bins_(),
  offset_(0.),
  // Histograms
  hPtHat_(),
  htHistos_(false),
  hHtDistr_(),
  mhtHistos_(false),
  hMhtDistr_(),
  multiHistos_(false),
  hMultiplicity_(),
  minBiasDeltaPhiHistos_(false),
  hMinBiasDeltaPhi_(),
  babyJetsHistos_(false),
  hBabyJets_(),
  hPassAlphaT_(),
  hPassDeadEcal_(),
  hPassBaby_(),
  hPassTrackless_(),
  hPassRecHit_(),
  hPassMinBiasDPhi_(),
  hFailAlphaT_(),
  hFailDeadEcal_(),
  hFailBaby_(),
  hFailTrackless_(),
  hFailMinBiasDPhi_(),
  // RECO cntrs
  cntr_baseline(0.),
  cntr_leading_jet_eta(0.),
  cntr_odd_jet_veto(0.),
  cntr_bad_muon_in_jet(0.),
  cntr_leading_jet_pt(0.),
  cntr_second_jet_pt(0.),
  cntr_correct_reco_bin(0.),
  cntr_jet_multiplicity(0.),
  cntr_at(0.),
  cntr_pass_dead_ecal(0.),
  cntr_pass_baby_cut(0.),
  cntr_pass_pt_over_ht(0.),
  signal_(),
  monitor_(),
  names_()
{
  
  vstring signal;
  vstring monitor;
  
  // Optional parameters
  if ( ps.Contains("DirName") ) dirName_ = ps.Get<std::string>("DirName");
  if ( ps.Contains("MinObjects") ) nMin_ = ps.Get<int>("MinObjects");
  if ( ps.Contains("MaxObjects") ) nMax_ = ps.Get<int>("MaxObjects");
  if ( ps.Contains("Verbose") ) verbose_ = ps.Get<bool>("Verbose");
  if ( ps.Contains("NoQCD") ) noQcd_ = ps.Get<bool>("NoQCD");
  if ( ps.Contains("NoCleaningInDenom") ) noCleaningInDenom_ = ps.Get<bool>("NoCleaningInDenom");
  if ( ps.Contains("Filter") ) filter_ = ps.Get<int>("Filter"); 
  if ( ps.Contains("ScaleInSignalRegion") ) scaleSignal_ = ps.Get<bool>("ScaleInSignalRegion"); 
  if ( ps.Contains("ScaleSingleJetThreshold") ) scaleSingle_ = ps.Get<bool>("ScaleSingleJetThreshold"); 
  if ( ps.Contains("SignalTriggers") ) signal = ps.Get< vstring >("SignalTriggers"); 
  if ( ps.Contains("MonitorTriggers") ) monitor = ps.Get< vstring >("MonitorTriggers"); 
  if ( ps.Contains("Inclusive") ) inclusive_ = ps.Get<bool>("Inclusive"); 
  if ( ps.Contains("AlphaTDefault") ) aT_ = ps.Get<double>("AlphaTDefault"); 
  if ( ps.Contains("HtOffset") ) offset_ = ps.Get<double>("HtOffset"); 
  
  if ( ps.Contains("HtHistos") ) htHistos_ = ps.Get<bool>("HtHistos"); 
  if ( ps.Contains("MhtHistos") ) mhtHistos_ = ps.Get<bool>("MhtHistos"); 
  if ( ps.Contains("MultiplicityHistos") ) multiHistos_ = ps.Get<bool>("MultiplicityHistos"); 
  if ( ps.Contains("BabyJetsHistos") ) babyJetsHistos_ = ps.Get<bool>("BabyJetsHistos"); 
  if ( ps.Contains("MinBiasDPhiHistos") ) minBiasDeltaPhiHistos_ = ps.Get<bool>("MinBiasDPhiHistos"); 

  if ( ps.Contains("AddMuon") ) addMuon_ = ps.Get<bool>("AddMuon"); 
  
  // Required parameters
  maxLeadingJetEta_ = ps.Get<double>("MaxLeadingJetEta");
  babyCut_ = ps.Get<double>("BabyJetsCut"); 
  alphaT_ = ps.Get< std::vector<double> >("AlphaT"); 
  scale_ = ps.Get<bool>("Scale"); 
  ht_ = ps.Get<int>("HtBin"); 
  pt1_ = ps.Get<double>("MinJetPt1");
  pt2_ = ps.Get<double>("MinJetPt2");
  pt3_ = ps.Get<double>("MinJetPt3");
  htNbins_ = (uint)ps.Get<int>("HtNbins"); 
  htXlow_ = ps.Get<double>("HtXlow"); 
  htXhigh_ = ps.Get<double>("HtXhigh"); 

  // Override binning if defined
  if ( ps.Contains("HtBinning") ) {
    htBins_ = ps.Get< std::vector<double> >("HtBinning"); 
    if ( !htBins_.empty() ) {
      htNbins_ = htBins_.size() - 1;
      htXlow_ = htBins_.front();
      htXhigh_ = htBins_.back();
    }
  }
  
  useFile_ = ps.Get<bool>("UseDeadEcalFile"); 
  deadEcalRegionDrCut_ = ps.Get<double>("DeadEcalRegionDrCut"); 
  deadEcalRegionEtaCut_ = ps.Get<double>("DeadEcalRegionEtaCut"); 
  minJetPtCut_ = ps.Get<double>("MinJetPtCut"); 
  nBadCells_ = ps.Get<int>("NBadCells"); 
  deadECALfile_ = ps.Get<std::string>("DeadEcalFile"); 
  minBiasDeltaPhiCut_ = ps.Get<double>("MinBiasDeltaPhiCut"); 
  dead_ = new DeadECALCut(deadEcalRegionEtaCut_,
			  deadEcalRegionDrCut_,
			  minBiasDeltaPhiCut_,
			  minJetPtCut_,
			  nBadCells_,
			  useFile_,
			  deadECALfile_);

  rechit_ = new sumRecHitPtCut( ps.Get<double>("SumRecHitPtCut") );

  // Parse signal trigger strings
  signal_.clear();
  signal_.resize( htNbins_, vstring(100) );
  unsigned int ii = 0;
  unsigned int jj = 0;
  vstring::const_iterator iter = signal.begin();
  vstring::const_iterator jter = signal.end();
  for ( ; iter != jter; ++iter ) { 
    if ( *iter == "" ) { signal_[ii].resize(jj); ii++; jj = 0; }
    else { 
      if ( ii < signal_.size() && jj < signal_[ii].size() ) { signal_[ii][jj] = *iter; } 
      else { break; }
      jj++; 
    }
  }
  if ( ii < htNbins_ ) { signal_.resize(ii); signal_.resize(htNbins_,vstring(signal_[ii-1])); }

  for ( vvstring::const_iterator i = signal_.begin(); i != signal_.end(); ++i ) {
    std::cout << i->size() << ": ";
    for ( vstring::const_iterator j = i->begin(); j != i->end(); ++j ) {
      std::cout << "\"" << *j << "\" ";
    }
    std::cout << std::endl;
  }

  // Parse monitor trigger strings
  monitor_.clear();
  monitor_.resize( htNbins_, vstring(100) );
  ii = 0;
  jj = 0;
  iter = monitor.begin();
  jter = monitor.end();
  for ( ; iter != jter; ++iter ) { 
    if ( *iter == "" ) { monitor_[ii].resize(jj); ii++; jj = 0; }
    else { 
      if ( ii < monitor_.size() && jj < monitor_[ii].size() ) { monitor_[ii][jj] = *iter; } 
      else { break; }
      jj++; 
    }
  }
  if ( ii < htNbins_ ) { monitor_.resize(ii); monitor_.resize(htNbins_,vstring(monitor_[ii-1])); }

  for ( vvstring::const_iterator i = monitor_.begin(); i != monitor_.end(); ++i ) {
    std::cout << i->size() << ": ";
    for ( vstring::const_iterator j = i->begin(); j != i->end(); ++j ) {
      std::cout << "\"" << *j << "\" ";
    }
    std::cout << std::endl;
  }

  // Define binning in terms of HT, Meff and pT thresholds
  binning();
  
  // Add 0. cut value 
  alphaT_.insert( alphaT_.begin(), 0. );

}

// -----------------------------------------------------------------------------
//
QcdBkgdEst::~QcdBkgdEst() {
  std::cout << std::fixed 
	    << "RECO cut flow:" << std::endl
    //<< "cntr_signal_region: " << cntr_baseline << std::endl
	    << "cntr_correct_reco_bin: " << cntr_correct_reco_bin << std::endl
	    << "cntr_leading_jet_eta: " << cntr_leading_jet_eta << std::endl
	    << "cntr_bad_muon_in_jet: " << cntr_bad_muon_in_jet << std::endl
	    << "cntr_odd_jet_veto: " << cntr_odd_jet_veto << std::endl
	    << "cntr_leading_jet_pt: " << cntr_leading_jet_pt << std::endl
	    << "cntr_second_jet_pt: " << cntr_second_jet_pt << std::endl
	    << "cntr_jet_multiplicity: " << cntr_jet_multiplicity << std::endl
	    << "cntr_at: " << cntr_at << std::endl
	    << "cntr_pass_dead_ecal: " << cntr_pass_dead_ecal << std::endl 
	    << "cntr_pass_baby_cut: " << cntr_pass_baby_cut << std::endl
	    << "cntr_pass_pt_over_ht: " << cntr_pass_pt_over_ht << std::endl;

  std::cout << " HT triggers found:  " << std::endl;
  vstring::const_iterator ii = names_.begin();
  vstring::const_iterator jj = names_.end();
  for ( ; ii != jj; ++ii ) {
    std::cout << *ii << std::endl;
  }

}

// -----------------------------------------------------------------------------
//
std::ostream& QcdBkgdEst::Description( std::ostream& ss ) {
  ss << "[QcdBkgdEst::Description]" 
     << std::endl
     << " Writing histograms to directory \"" << dirName_ << "\""
     << "for events with " << nMin_ << "<objects<" << nMax_ << "."
     << std::endl
     << " Verbose=" << verbose_ 
     << std::endl;

  ss << " Thresholds scaled according to HT bins:" << std::endl
     << " --------------------------------------------------------" << std::endl
     << " |  HT  |    (pT1,pT2,pT3) | Meff |          (x1,x2,x3) |" << std::endl
     << " --------------------------------------------------------" << std::endl;
  for ( uint ii = 0; ii < htNbins_; ++ii ) {
    double ht = htBins_[ii];
    double pt1 = htPt1Bins_[ii];
    double pt2 = htPt2Bins_[ii];
    double pt3 = htPt3Bins_[ii];
    double meff = htBins_[ii] + htPt3Bins_[ii];
    ss << " | " << std::setw(4) << std::setprecision(3) << ht
       << " | (" << std::setw(4) << std::setprecision(3) << pt1
       << "," << std::setw(4) << std::setprecision(3) << pt2
       << "," << std::setw(4) << std::setprecision(3) << pt3
       << ") | " << std::setw(4) << std::setprecision(3) << meff
       << " | (" << std::setw(5) << std::setprecision(3) << (norm_*pt1/meff)
       << "," << std::setw(5) << std::setprecision(3) << (norm_*pt2/meff) 
       << "," << std::setw(5) << std::setprecision(3) << (norm_*pt3/meff) 
       << ") |" << std::endl;
  }
  ss << " --------------------------------------------------------" << std::endl;
  
  ss << std::fixed << std::setprecision(8);
  return ss;

}

// -----------------------------------------------------------------------------
//
void QcdBkgdEst::binning() {

  // HT binning
  if ( htBins_.empty() ) {
    for ( uint ii = 0; ii < htNbins_; ++ii ) { 
      double width = ( htXhigh_ - htXlow_ ) / htNbins_;
      htBins_.push_back(htXlow_+width*ii);
    }
    htBins_.push_back(htXhigh_);
  }

  // Normalisation: x_{1} + x_{2} + x_{3} = norm

  // pT thresholds that scale with HT
  double meff_default = htBins_[ht_] + pt3_;
  double x3_default = ( norm_ * pt3_ ) / meff_default; 
  double x3_factor = ( norm_ - x3_default ) / x3_default;
  for ( unsigned int iht = 0; iht < htNbins_; ++iht ) {
    double ht = htBins_[iht];
    double pt3 = ht / x3_factor;
    double meff = pt3 + ht;
    double pt1 = pt1_ * meff / meff_default;
    double pt2 = pt2_ * meff / meff_default;
    htPt1Bins_.push_back( scale_ && ( ht <= htBins_[ht_] || scaleSignal_ ) ? pt1 : pt1_ );
    htPt2Bins_.push_back( scale_ && ( ht <= htBins_[ht_] || scaleSignal_ ) ? pt2 : pt2_ );
    htPt3Bins_.push_back( scale_ && ( ht <= htBins_[ht_] || ( scaleSignal_ && scaleSingle_ ) ) ? pt3 : pt3_ );
  }

  // Apply offset to HT binning
  for ( uint ii = 0; ii < htNbins_; ++ii ) { htBins_[ii] += offset_; }

  std::cout << (*this);
  
}

// -----------------------------------------------------------------------------
//
bool QcdBkgdEst::Process( Event::Data& ev ) {

  //std::cout << "[QcdBkgdEst::Process] " << ev.GetEventNumber() << std::endl;
  
  // Use as filter for interesting events
  bool keep_event = false;

  // Event weight
  Double_t weight = ev.GetEventWeight();

  // Histogram of pT hat 
  try {
    double pt_hat = ev.pthat(); 
    if ( hPtHat_[0] ) hPtHat_[0]->Fill( pt_hat, weight ); 
  } catch (...){}
  
  // MET
  LorentzV calo_met = *(ev.metP4caloTypeII());

  // Retrieve "common" jets from ntuple
  std::vector<Event::Jet const*> common = ev.JD_CommonJets().accepted;

  // Retrieve "odd" jets from ntuple
  std::vector<Event::Jet const*> odd = ev.JD_CommonJets().odd;

  // Number of HT/Meff bins 
  int nbins = htNbins_; 

  // MHT in signal region
  double baby_val = 1000.;
  if ( ht_ < htNbins_ ) { 
    
    // Undo muon correction in CaloMET (to be used only when muon vetoes removed!)
    LorentzV met = calo_met;
    if (addMuon_) {
      std::vector<Event::Lepton const*> muons = ev.LD_CommonMuons().accepted;
      std::vector<Event::Lepton const*>::const_iterator imuon = muons.begin();
      std::vector<Event::Lepton const*>::const_iterator jmuon = muons.end();
      for ( ; imuon != jmuon; ++imuon ) { if ( *imuon ) { met += **imuon; } }
    }

    LorentzV mht_reco_signal(0.,0.,0.,0.);
    std::vector<Event::Jet const*>::const_iterator ireco = common.begin();
    std::vector<Event::Jet const*>::const_iterator jreco = common.end();
    for ( ; ireco != jreco; ++ireco ) {
      if ( (*ireco)->Pt() > htPt3Bins_[ht_] ) { mht_reco_signal -= **ireco; }
    }
    baby_val = met.Pt() > 0. ? ( mht_reco_signal.Pt() / met.Pt() ) : 1000.; 
  }

  
  // (Reverse!) Loop through HT/Meff bins 
  for ( int jbin = 0; jbin < nbins; ++jbin ) {
    int ibin = (nbins-1) - jbin;

    // Define jet pT thresholds 
    double pt1_threshold = htPt1Bins_[ibin];
    double pt2_threshold = htPt2Bins_[ibin];
    double pt3_threshold = htPt3Bins_[ibin];
    
    // -------------------- RECO INFO --------------------

    // Calculate aT, HT, Meff and MHT from common jets
    std::vector<LorentzV> reco; 
    reco.reserve(10);
    LorentzV mht_reco(0.,0.,0.,0.);
    double at_reco = 0.;
    double ht_reco = 0.;
    double meff_reco = 0.;
    double x1_reco = 0.;
    double x2_reco = 0.;
    
    // Baby jets
    std::vector<LorentzV> baby_jets; 
    baby_jets.reserve(10);
    LorentzV baby_mht(0.,0.,0.,0.);
    double baby_ht = 0.;
    
    // Loop through common jets and calc "reco" variables
    std::vector<Event::Jet const*>::const_iterator ireco = common.begin();
    std::vector<Event::Jet const*>::const_iterator jreco = common.end();
    for ( ; ireco != jreco; ++ireco ) {
      if ( (*ireco)->Pt() > pt3_threshold ) { 
	reco.push_back(**ireco); 
	ht_reco += (*ireco)->Et();
	mht_reco -= **ireco;
      } else if ( (*ireco)->Pt() > 10. ) {
	baby_jets.push_back(**ireco); 
	baby_ht += (*ireco)->Et();
	baby_mht -= **ireco;
      }
    }

    // Meff
    meff_reco = ht_reco + mht_reco.Pt();

    // Calc AlphaT and Dalitz variables
    std::vector<bool> pseudo_reco;
    at_reco = AlphaT()( reco, pseudo_reco );
    calcDalitzVars( reco, pseudo_reco, mht_reco, x1_reco, x2_reco );
      
    // Odd jet veto
    bool odd_jet_veto = false;
    std::vector<Event::Jet const*>::const_iterator iodd = odd.begin();
    std::vector<Event::Jet const*>::const_iterator jodd = odd.end();
    for ( ; iodd != jodd; ++iodd ) {
      if ( (*iodd)->Pt() >= pt3_threshold ) { odd_jet_veto = true; }
    }
      
    // Bad muon in jets
    bool bad_muon_in_jet = false;
    std::vector<Event::Lepton>::const_iterator imuon = ev.LD_Muons().begin();
    std::vector<Event::Lepton>::const_iterator jmuon = ev.LD_Muons().end();
    for ( ; imuon != jmuon; ++imuon ) {
      if( imuon->Pt() < 10. ) { continue; }
      bool found_match = false;
      std::vector<LorentzV>::const_iterator ijet = reco.begin();
      std::vector<LorentzV>::const_iterator jjet = reco.end();
      for ( ; ijet != jjet; ++ijet ) {
	if ( fabs( ROOT::Math::VectorUtil::DeltaR(*ijet,*imuon) ) < 0.5 ) {
	  if ( imuon->Pt() > 0.5 * ijet->Pt() ) { bad_muon_in_jet = true; }
	  found_match = true;
	}
      }
      //if ( !found_match && !imuon->WasItcc() ) { bad_muon_in_jet = true; }
    }

    // Order jets by pT
    std::sort( reco.begin(), reco.end(), SortByPt );
      
    // Check multiplicity and pT of leading and second jets
    bool jet_multiplicity = ( reco.size() >= 2 );
    bool leading_jet_pt   = ( reco.size() > 0 && reco[0].Pt() > pt1_threshold );
    bool second_jet_pt    = ( reco.size() > 1 && reco[1].Pt() > pt2_threshold );
    bool leading_jet_eta  = ( reco.size() > 0 && fabs(reco[0].Eta()) < maxLeadingJetEta_ );
    
    // Check if MHT/MET (baby jet cleaning) passes cut
    //if ( true ) baby_val = meff_reco > 0. ? ( mht_reco.Pt() - calo_met.Pt() ) / meff_reco  : 1000.; 
    bool pass_baby_cut = ( baby_val < babyCut_ );       
    
    // Check if event falls within appropriate HT or Meff bin (last bin is inclusive)
    bool correct_reco_bin = false;
    double ht_lower = htBins_[ibin];
    double ht_upper = ( (nbins-ibin) != 1 ? htBins_[ibin+1] : 10000. ); 
    if ( inclusive_ ) { if ( ht_reco > ht_lower ) { correct_reco_bin = true; } }
    else { if ( ht_reco > ht_lower && ht_reco < ht_upper ) { correct_reco_bin = true; } }
    
//     std::cout << " event " << ev.GetEventNumber()
// 	      << " ibin " << ibin
// 	      << " htBins_[ibin] " << htBins_[ibin]
// 	      << " ht_lower " << ht_lower
// 	      << " ht_upper " << ht_upper
// 	      << " ht_reco " << ht_reco
// 	      << " reco.size() " << reco.size()
// 	      << " nodd " << odd.size()
// 	      << " correct_reco_bin " << correct_reco_bin
// 	      << std::endl;
    
    // Check if event passes "dead ECAL" cut
    bool pass_dead_ecal = dead_ ? dead_->Process( ev, reco, baby_jets ) : true;
      
    // Track Sum Pt over HT cut
    bool pass_pt_over_ht = passVertexSumPtOverHt( ev, ht_reco );
    
    // Min biased delta phi
    double min_bias_dphi = minBiasDeltaPhi( reco );
    bool pass_min_bias_dphi = ( min_bias_dphi > 0.5 );
    
    // Sum RecHit Pt cleaning cut
    bool pass_rechit = rechit_ ? rechit_->Process(ev) : true; 
    
    // -------------------- Debug --------------------

//     //if ( true && correct_reco_bin ) {
//     if ( true ) {
//       std::cout << "[QcdBkgdEst::Process] RECO info:" << std::endl
// 		<< " Entry: " << ev.GetCurrentEntry() << std::endl
// 		<< " Event: " << ev.GetEventNumber() << std::endl
// 		<< " Njets: " << reco.size() << std::endl
// 		<< " HT: " << ht_reco << std::endl
// 		<< " (HT bin: " << ibin << ")" << std::endl
// 		<< " MHT: " << mht_reco.Pt() << std::endl
// 		<< " AlphaT: " << at_reco << std::endl
// 		<< " Meff: " << meff_reco << std::endl
// 		<< " MET consistency: " << baby_val << std::endl
// 		<< " x1: " << x1_reco << std::endl
// 		<< " x2: " << x2_reco << std::endl
// 		<< " Njets>=2: " << jet_multiplicity << std::endl
// 		<< " pT1>threshold: " << leading_jet_pt << std::endl
// 		<< " pT2>threshold: " << second_jet_pt << std::endl
// 		<< " eta1<threshold: " << leading_jet_eta << std::endl
// 		<< " OddJet: " << odd_jet_veto << std::endl
// 		<< " BadMuonInJet: " << bad_muon_in_jet << std::endl
// 		<< " PassDeadECAL: " << pass_dead_ecal << std::endl
// 		<< " PassMetConsistency: " << pass_baby_cut << std::endl;
//     }
    
//     if ( correct_reco_bin ) {
//       if ( !jet_multiplicity ||
// 	   !leading_jet_pt ||
// 	   !second_jet_pt ||
// 	   !leading_jet_eta ||
// 	   odd_jet_veto || 
// 	   bad_muon_in_jet ||
// 	   at_reco < aT_ || 
// 	   !pass_dead_ecal ||
// 	   !pass_baby_cut ||
// 	   !pass_pt_over_ht ) { 
// 	int signal = trigger( ev, signal_[ibin] );
// 	int monitor = trigger( ev, monitor_[ibin] );
// 	std::cout << "[QcdBkgdEst::Process] Event Summary:" << std::endl
// 		  << " Entry: " << ev.GetCurrentEntry() << std::endl
// 		  << " Event: " << ev.GetEventNumber() << std::endl
// 		  << " Signal trigger prescale: " << signal
// 		  << " Monitor trigger prescale: " << monitor
// 		  << " Njets: " << reco.size() << std::endl
// 		  << " HT: " << ht_reco << std::endl
// 		  << " (HT bin: " << ibin << ")" << std::endl
// 		  << " MHT: " << mht_reco.Pt() << std::endl
// 		  << " AlphaT: " << at_reco << std::endl
// 		  << " Meff: " << meff_reco << std::endl
// 		  << " MET consistency: " << baby_val << std::endl
// 		  << " x1: " << x1_reco << std::endl
// 		  << " x2: " << x2_reco << std::endl
// 		  << " njets: " << reco.size() << std::endl;
// 	std::vector<LorentzV>::const_iterator ireco = reco.begin();
// 	std::vector<LorentzV>::const_iterator jreco = reco.end();
// 	for ( ; ireco != jreco; ++ireco ) {
// 	  std::cout << "  ET: " << ireco->Et() 
// 		    << "  E: " << ireco->E() 
// 		    << "  pT: " << ireco->Pt() 
// 		    << " eta: " << ireco->eta()
// 		    << " phi: " << ireco->phi() 
// 		    << std::endl;
// 	}
// 	std::cout << jet_multiplicity << " "
// 		  << leading_jet_pt << " "
// 		  << second_jet_pt << " "
// 		  << leading_jet_eta << " "
// 		  << !odd_jet_veto << " " 
// 		  << !bad_muon_in_jet << " "
// 		  << ( at_reco > aT_ ) << " " 
// 		  << pass_dead_ecal << " "
// 		  << pass_baby_cut << " "
// 		  << pass_pt_over_ht << " "
// 		  << std::endl;
//       }
//     }

    if ( //mht_reco_signal.Pt() > 100 && //@@ ADDITONAL CUT REMOVE REMOVE !!!!!!!!
	 correct_reco_bin && 
	 jet_multiplicity &&
	 leading_jet_pt &&
	 second_jet_pt &&
	 leading_jet_eta &&
	 !odd_jet_veto && 
	 !bad_muon_in_jet ) {
      
      // Retrieve trigger and prescale info from event
      int signal = trigger( ev, signal_[ibin] );
      int monitor = trigger( ev, monitor_[ibin] );
      
      // Decide if event is interesting or not (-1 = no filter; >nbins = all bins; other = particular bin)
      if ( signal > 0 && 
	   at_reco > aT_ && 
	   pass_dead_ecal &&
	   pass_baby_cut &&
	   pass_pt_over_ht ) { 
	
	//fill( reco.size(), hMinBiasDeltaPhi_[0], min_bias_dphi, weight*signal ); 
	//fill( reco.size(), hMinBiasDeltaPhi_[ibin+1], min_bias_dphi, weight*signal ); 

	if ( filter_ >= nbins || ibin == filter_ ) { 
	  keep_event = true; 
	  if ( verbose_ ) { 
	    std::cout << "[QcdBkgdEst::Process] Event Summary:" << std::endl
		      << " Entry: " << ev.GetCurrentEntry() << std::endl
		      << " Event: " << ev.GetEventNumber() << std::endl
		      << " Signal trigger prescale: " << signal
		      << " Monitor trigger prescale: " << monitor
		      << " Njets: " << reco.size() << std::endl
		      << " HT: " << ht_reco << std::endl
		      << " (HT bin: " << ibin << ")" << std::endl
		      << " MHT: " << mht_reco.Pt() << std::endl
		      << " AlphaT: " << at_reco << std::endl
		      << " Meff: " << meff_reco << std::endl
		      << " MET consistency: " << baby_val << std::endl
		      << " x1: " << x1_reco << std::endl
		      << " x2: " << x2_reco << std::endl
		      << " njets: " << reco.size() << std::endl;
	    std::vector<LorentzV>::const_iterator ireco = reco.begin();
	    std::vector<LorentzV>::const_iterator jreco = reco.end();
	    for ( ; ireco != jreco; ++ireco ) {
	      std::cout << "  ET: " << ireco->Et() 
			<< "  E: " << ireco->E() 
			<< "  pT: " << ireco->Pt() 
			<< " eta: " << ireco->eta()
			<< " phi: " << ireco->phi() 
			<< std::endl;
	    }
	  }
	}
      }
      
//       // Determine HT used to fill histograms
//       double ht_binned = inclusive_ ? htBins_[ibin] + 1.e-3 : ht_reco;
      
//       // HT distributions
//       for ( uint iat = 0; iat < ( noQcd_? 1 : alphaT_.size() ); ++iat ) {
	
// 	// Check AlphaT is above threshold
// 	bool pass_at = ( at_reco > alphaT_[iat] );
	
// 	// Check if appropriate trigger fired
// 	int prescale = -1;
// 	if      ( iat == 0 && monitor > 0 ) { prescale = monitor; }
// 	else if (  pass_at &&  signal > 0 ) { prescale =  signal; }
// 	else if ( !pass_at && monitor > 0 ) { prescale = monitor; }
// 	else                                { continue; }
	
// 	// Adjust event weight by prescale
// 	double w = weight;//*prescale;
	
// 	if ( pass_at ) {
// 	  fill( reco.size(), hPassAlphaT_[iat], ht_binned, w ); 
// 	  bool ok = ( noCleaningInDenom_ && ( iat == 0 ) );
// 	  if ( babyJetsHistos_ ) {
// 	    fill( reco.size(), hBabyJets_[0][iat], baby_val, w ); 
// 	    fill( reco.size(), hBabyJets_[ibin+1][iat], baby_val, w ); 
// 	  }
// 	  if ( pass_dead_ecal || ok ) { 
// 	    fill( reco.size(), hPassDeadEcal_[iat], ht_binned, w ); 
// 	    if ( pass_baby_cut || ok ) { 
// 	      fill( reco.size(), hPassBaby_[iat], ht_binned, w ); 
// 	      if ( pass_pt_over_ht || ok ) { 
// 		fill( reco.size(), hPassTrackless_[iat], ht_binned, w ); 
// 		if ( pass_min_bias_dphi || ok ) { 
// 		  fill( reco.size(), hPassMinBiasDPhi_[iat], ht_binned, w ); 
// 		} else {
// 		  fill( reco.size(), hFailMinBiasDPhi_[iat], ht_binned, w ); 
// 		}
// 	      } else {
// 		fill( reco.size(), hFailTrackless_[iat], ht_binned, w ); 
// 		fill( reco.size(), hFailMinBiasDPhi_[iat], ht_binned, w ); 
// 	      }
// 	    } else {
// 	      fill( reco.size(), hFailBaby_[iat], ht_binned, w ); 
// 	      fill( reco.size(), hFailTrackless_[iat], ht_binned, w ); 
// 	      fill( reco.size(), hFailMinBiasDPhi_[iat], ht_binned, w ); 
// 	    }
// 	  } else {
// 	    fill( reco.size(), hFailDeadEcal_[iat], ht_binned, w ); 
// 	    fill( reco.size(), hFailBaby_[iat], ht_binned, w ); 
// 	    fill( reco.size(), hFailTrackless_[iat], ht_binned, w ); 
// 	    fill( reco.size(), hFailMinBiasDPhi_[iat], ht_binned, w ); 
// 	  }
// 	} else { 
// 	  fill( reco.size(), hFailAlphaT_[iat], ht_binned, w ); 
// 	  fill( reco.size(), hFailDeadEcal_[iat], ht_binned, w ); 
// 	  fill( reco.size(), hFailBaby_[iat], ht_binned, w ); 
// 	  fill( reco.size(), hFailTrackless_[iat], ht_binned, w ); 
// 	  fill( reco.size(), hFailMinBiasDPhi_[iat], ht_binned, w ); 
// 	}
	
//       } // iat

      // HT distributions
      for ( uint iat = 0; iat < ( noQcd_? 1 : alphaT_.size() ); ++iat ) {
	
	// Check AlphaT is above threshold
	if ( at_reco < alphaT_[iat] ) { continue; }
	
	// Extract prescale to weight histograms
	double w = weight;
	if      ( iat  > 0 &&  signal > 0 ) { w *=  signal*1.; }
	else if ( iat == 0 && monitor > 0 ) { w *= monitor*1.; }
	else                                { continue; }
	
	// Defines HT bin to be used in case of inclusive binning
	double ht_binned = inclusive_ ? htBins_[ibin] + 1.e-3 : ht_reco;
	
	// Histograms after "bare" AlphaT cut
	fill( reco.size(), hPassAlphaT_[iat], ht_binned, w ); 
	
	if ( mhtHistos_ ) {
	  fill( reco.size(), hMhtDistr_[0][iat], mht_reco.Pt(), w ); 
	  fill( reco.size(), hMhtDistr_[ibin+1][iat], mht_reco.Pt(), w ); 
	}
		
	// Histograms filled after cleaning cuts
	bool ok = ( noCleaningInDenom_ && ( iat == 0 ) );
	if ( pass_dead_ecal || ok ) { 
	  fill( reco.size(), hPassDeadEcal_[iat], ht_binned, w ); 
	  
	  // Additional histograms
	  if ( babyJetsHistos_ ) {
	    fill( reco.size(), hBabyJets_[0][iat], baby_val, w ); 
	    fill( reco.size(), hBabyJets_[ibin+1][iat], baby_val, w ); 
	  }
	  
	  if ( pass_baby_cut || ok ) { 
	    fill( reco.size(), hPassBaby_[iat], ht_binned, w ); 

	    if ( pass_pt_over_ht || ok ) { 
	      fill( reco.size(), hPassTrackless_[iat], ht_binned, w ); 
	      
	      if ( pass_rechit || ok ) { 
		fill( reco.size(), hPassRecHit_[iat], ht_binned, w ); 
		
		// Additional histograms 
		if ( htHistos_ ) {
		  fill( reco.size(), hHtDistr_[0][iat], ht_reco, w ); 
		  fill( reco.size(), hHtDistr_[ibin+1][iat], ht_reco, w ); 
		}
		if ( multiHistos_ ) {
		  fill( reco.size(), hMultiplicity_[0][iat], reco.size(), w ); 
		  fill( reco.size(), hMultiplicity_[ibin+1][iat], reco.size(), w ); 
		}
		if ( minBiasDeltaPhiHistos_ ) {
		  fill( reco.size(), hMinBiasDeltaPhi_[0][iat], min_bias_dphi, w ); 
		  fill( reco.size(), hMinBiasDeltaPhi_[ibin+1][iat], min_bias_dphi, w ); 
		}
		
		if ( pass_min_bias_dphi || ok ) { 
		  fill( reco.size(), hPassMinBiasDPhi_[iat], ht_binned, w ); 
		}
	      }
	    }
	  }
	}
      }

    } // pre-selection
	  
    // Check cut flow efficiencies
    if ( (nbins-ibin) == 1 ) {
      //cntr_baseline += weight;
      if ( correct_reco_bin ) {
	cntr_correct_reco_bin += weight;
	if ( leading_jet_eta ) {
	  cntr_leading_jet_eta += weight;
	  if ( !bad_muon_in_jet ) {
	    cntr_bad_muon_in_jet += weight;
	    if ( !odd_jet_veto ) { 
	      cntr_odd_jet_veto += weight;
	      if ( leading_jet_pt ) {
		cntr_leading_jet_pt += weight;
		if ( second_jet_pt ) {
		  cntr_second_jet_pt += weight;
		  if ( jet_multiplicity ) {
		    cntr_jet_multiplicity += weight;
		    if ( at_reco > aT_ ) {
		      cntr_at += weight;
		      if ( pass_dead_ecal ) {
			cntr_pass_dead_ecal += weight;
			if ( pass_baby_cut ) {
			  cntr_pass_baby_cut += weight;
			  if ( pass_pt_over_ht ) {
			    cntr_pass_pt_over_ht += weight;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
	  
    // Only consider each event once only
    if ( inclusive_ && correct_reco_bin ) { break; }
    
  } // nBins

  return keep_event;
  
}

// -----------------------------------------------------------------------------
//
void QcdBkgdEst::calcDalitzVars( const std::vector<LorentzV>& jets, 
				 const std::vector<bool>& pseudo,
				 const LorentzV& mht, 
				 double& x1, 
				 double& x2 ) {
  
  // Construct pseudo di-jets
  std::vector<LorentzV> dijets(2,LorentzV(0.,0.,0.,0.));
  for ( unsigned int i = 0; i < jets.size(); ++i ) {
    if ( pseudo[i] ) { dijets[0] += jets[i]; }
    else { dijets[1] += jets[i]; }
  }
  std::sort( dijets.begin(), dijets.end(), SortByPt );
  
  // Calculate x1, x2
  x1 = ( norm_ * dijets[0].Pt() ) / ( dijets[0].Pt() + dijets[1].Pt() + mht.Pt() );
  x2 = ( norm_ * dijets[1].Pt() ) / ( dijets[0].Pt() + dijets[1].Pt() + mht.Pt() );

}

// // -----------------------------------------------------------------------------
// //
// void QcdBkgdEst::fill( uint multiplicity,
// 		       double alpha_t,
// 		       double at_cut,
// 		       std::vector<TH1D*>& his,
// 		       double val,
// 		       double weight ) {
  
//   // Events passing AlphaT
//   if ( alpha_t >= at_cut ) {
//     if ( !his.empty() && 
// 	 multiplicity >= nMin_ ) { 
//       if ( his[0] ) his[0]->Fill( val, weight ); 
//     }
//     if ( multiplicity >= nMin_ && 
// 	 multiplicity <= nMax_ &&
// 	 multiplicity < his.size() ) {
//       if ( his[multiplicity] ) his[multiplicity]->Fill( val, weight ); 
//     }
//   }

// }

// -----------------------------------------------------------------------------
//
void QcdBkgdEst::fill( uint multiplicity,
		       std::vector<TH1D*>& his,
		       double valx,
		       double weight ) {
  
  if ( his.empty() ) { return; }
  if ( multiplicity < nMin_ ) { return; }
  
  if ( his[0] ) his[0]->Fill( valx, weight ); 
  
  if ( his.size() == 1 ) { return; }
  //if ( multiplicity >= his.size() ) { return; }

  if ( multiplicity < nMax_ ) { 
    if ( his[multiplicity] ) his[multiplicity]->Fill( valx, weight ); 
  } else { 
    if ( his[nMax_] ) his[nMax_]->Fill( valx, weight ); 
  }
  
}

// -----------------------------------------------------------------------------
//
void QcdBkgdEst::fill( uint multiplicity,
		       std::vector<TH2D*>& his,
		       double valx,
		       double valy,
		       double weight ) {
  
  if ( his.empty() ) { return; }
  if ( multiplicity < nMin_ ) { return; }
  
  if ( his.size() == 1 ) { return; }
  
  if ( his[0] ) his[0]->Fill( valx, valy, weight ); 
  if ( multiplicity < nMax_ ) { 
    if ( his[multiplicity] ) his[multiplicity]->Fill( valx, valy, weight ); 
  } else  { 
    if ( his[nMax_] ) his[nMax_]->Fill( valx, valy, weight ); 
  }

}

// -----------------------------------------------------------------------------
//
void QcdBkgdEst::fill2d( uint multiplicity,
			 std::vector<TH2D*>& his,
			 double x, 
			 double y,
			 double weight ) {
  
  if ( !his.empty() && 
       multiplicity >= nMin_ ) { 
    if ( his[0] ) { his[0]->Fill( x, y, weight ); }
  }
  if ( multiplicity >= nMin_ && 
       multiplicity <= nMax_ &&
       multiplicity < his.size() ) { 
    if ( his[multiplicity] ) { his[multiplicity]->Fill( x, y, weight ); }
  }

}


// -----------------------------------------------------------------------------
//
void QcdBkgdEst::Start( Event::Data& ev ) {
  
  // Top-level directory
  initDir( ev.OutputFile(), dirName_.c_str() );

  try {
    ev.pthat(); 
    BookHistArray( hPtHat_,
		   "GenPtHat",
		   ";#hat{p_{T}} [GeV];", 
		   400,0.,2000.,
		   1, 0, 1, true );
  } catch (...){}

  if ( multiHistos_ ) {
    hMultiplicity_.resize( htNbins_+1, vvTH1D( alphaT_.size(), vTH1D() ) );
    for ( uint ii = 0; ii <= htNbins_; ++ii ) { 
      for ( uint jj = 0; jj < alphaT_.size(); ++jj ) { 
	std::stringstream ss;
	ss << "HtMultiplicity_" << ht(ii) << "_" << at(jj);
	BookHistArray( hMultiplicity_[ii][jj], 
		       ss.str(),
		       ";N_{jets};", 
		       21, -0.5, 20.5, 
		       1, 0, 1, true );
      }
    }
  }
  
  if ( htHistos_ ) {
    hHtDistr_.resize( htNbins_+1, vvTH1D( alphaT_.size(), vTH1D() ) );
    for ( uint ii = 0; ii <= htNbins_; ++ii ) { 
      for ( uint jj = 0; jj < alphaT_.size(); ++jj ) { 
	std::stringstream ss;
	ss << "HtDistrPassTrackless_" << ht(ii) << "_" << at(jj);
	BookHistArray( hHtDistr_[ii][jj],
		       ss.str(),
		       ";HT (GeV);", 
		       int((htXhigh_-htXlow_)/10.),  htXlow_, htXhigh_, 
		       nMax_+1, 0, 1, true );
      }
    }
  }

  if ( minBiasDeltaPhiHistos_ ) {
    hMinBiasDeltaPhi_.resize( htNbins_+1, vvTH1D( alphaT_.size(), vTH1D() ) );
    for ( uint ii = 0; ii <= htNbins_; ++ii ) { 
      for ( uint jj = 0; jj < alphaT_.size(); ++jj ) { 
	std::stringstream ss;
	ss << "MinBiasDeltaPhi_" << ht(ii) << "_" << at(jj);
	BookHistArray( hMinBiasDeltaPhi_[ii][jj], 
		       ss.str(),
		       ";#Delta#phi^{*}_{min} [radians];", 
		       32, 0., 3.2,
		       nMax_+1, 0, 1, true );
      }
    }
  }
  
  if ( babyJetsHistos_ ) {
    hBabyJets_.resize( htNbins_+1, vvTH1D( alphaT_.size(), vTH1D() ) );
    for ( uint ii = 0; ii <= htNbins_; ++ii ) { 
      for ( uint jj = 0; jj < alphaT_.size(); ++jj ) { 
	std::stringstream ss;
	ss << "BabyJets_" << ht(ii) << "_" << at(jj);
	BookHistArray( hBabyJets_[ii][jj],
		       ss.str(),
		       ";Baby jet cleaning variable;", 
		       200,0.,10.,
		       nMax_+1, 0, 1, true );
      }
    }
  }
  
  hPassAlphaT_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterAlphaT_" << at(ii);
    BookHistArray( hPassAlphaT_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }
  
  hPassDeadEcal_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterDeadEcal_" << at(ii);
    BookHistArray( hPassDeadEcal_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }
  
  hPassBaby_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterBaby_" << at(ii);
    BookHistArray( hPassBaby_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }
  
  hPassTrackless_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterTrackless_" << at(ii);
    BookHistArray( hPassTrackless_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }

  hPassRecHit_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterRecHit_" << at(ii);
    BookHistArray( hPassRecHit_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }
  
  hPassMinBiasDPhi_.resize( alphaT_.size(), vTH1D() );
  for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
    std::stringstream ss;
    ss << "HtAfterMinBiasDPhi_" << at(ii);
    BookHistArray( hPassMinBiasDPhi_[ii], 
		   ss.str(),
		   ";HT (GeV);", 
		   htNbins_, &htBins_.front(), 
		   nMax_+1, 0, 1, true );
  }

//   hFailAlphaT_.resize( alphaT_.size(), vTH1D() );
//   for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
//     std::stringstream ss;
//     ss << "HtFailingAlphaT_" << at(ii);
//     BookHistArray( hFailAlphaT_[ii], 
// 		   ss.str(),
// 		   ";HT (GeV);", 
// 		   htNbins_, &htBins_.front(), 
// 		   nMax_+1, 0, 1, true );
//   }
  
//   hFailDeadEcal_.resize( alphaT_.size(), vTH1D() );
//   for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
//     std::stringstream ss;
//     ss << "HtFailingDeadEcal_" << at(ii);
//     BookHistArray( hFailDeadEcal_[ii], 
// 		   ss.str(),
// 		   ";HT (GeV);", 
// 		   htNbins_, &htBins_.front(), 
// 		   nMax_+1, 0, 1, true );
//   }
  
//   hFailBaby_.resize( alphaT_.size(), vTH1D() );
//   for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
//     std::stringstream ss;
//     ss << "HtFailingBaby_" << at(ii);
//     BookHistArray( hFailBaby_[ii], 
// 		   ss.str(),
// 		   ";HT (GeV);", 
// 		   htNbins_, &htBins_.front(), 
// 		   nMax_+1, 0, 1, true );
//   }
  
//   hFailTrackless_.resize( alphaT_.size(), vTH1D() );
//   for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
//     std::stringstream ss;
//     ss << "HtFailingTrackless_" << at(ii);
//     BookHistArray( hFailTrackless_[ii], 
// 		   ss.str(),
// 		   ";HT (GeV);", 
// 		   htNbins_, &htBins_.front(), 
// 		   nMax_+1, 0, 1, true );
//   }
  
//   hFailMinBiasDPhi_.resize( alphaT_.size(), vTH1D() );
//   for ( uint ii = 0; ii < alphaT_.size(); ++ii ) { 
//     std::stringstream ss;
//     ss << "HtFailingMinBiasDPhi_" << at(ii);
//     BookHistArray( hFailMinBiasDPhi_[ii], 
// 		   ss.str(),
// 		   ";HT (GeV);", 
// 		   htNbins_, &htBins_.front(), 
// 		   nMax_+1, 0, 1, true );
//   }
  
  if ( mhtHistos_ ) {
    hMhtDistr_.resize( htNbins_+1, vvTH1D( alphaT_.size(), vTH1D() ) );
    for ( uint ii = 0; ii <= htNbins_; ++ii ) { 
      for ( uint jj = 0; jj < alphaT_.size(); ++jj ) { 
	std::stringstream ss;
	ss << "MhtDistrPassAlphaT_" << ht(ii) << "_" << at(jj);
	BookHistArray( hMhtDistr_[ii][jj],
		       ss.str(),
		       ";MHT (GeV);", 
		       1000, 0., 1000., 
		       nMax_+1, 0, 1, true );
      }
    }
  }
  
}

// -----------------------------------------------------------------------------
//
std::string QcdBkgdEst::ht( int ii ) {
  std::stringstream ss;
  if ( ii == 0 ) { ss << "HT0"; }
  else { ss << "HT" << (int(10.*htBins_[ii-1])/10.); }
  return ss.str();
}

// -----------------------------------------------------------------------------
//
std::string QcdBkgdEst::at( int ii ) {
  std::stringstream ss;
  ss << "aT" << alphaT_[ii]; 
  return ss.str();
}

// -----------------------------------------------------------------------------
//
Double_t QcdBkgdEst::DeadECALDR(TString deadECALfile,float jetphi,float jeteta, int nBad){
  // cout << " hi ho " << deadECALfile << endl;
  TLorentzVector testjetVec(0,0,0,0);
  double DRmin = 100;
  testjetVec.SetPtEtaPhiE(100,jeteta,jetphi,100);

  //cout << " in dead ecal " << endl;
  ifstream myfile(deadECALfile);
  if(!myfile.good()) cout << "file is not there " << endl;
  if(myfile.is_open()){

    while(!myfile.eof()){

      float phi,eta;
      int iEta,iPhi,nBadCells;
      myfile>>eta;
      myfile>>phi;
      myfile>>iEta;
      myfile>>iPhi;
      myfile>>nBadCells;
      // cout <<eta<< " "<< phi<<endl;
      if(nBadCells > nBad){
	TLorentzVector testVec(0,0,0,0);
	testVec.SetPtEtaPhiE(100,eta,phi,100);
	float DR = fabs(ROOT::Math::VectorUtil::DeltaR(testVec,testjetVec));
	if(DR < DRmin) DRmin = DR;
      }
    }
   
    myfile.close();
  } 
  return DRmin;

}

// -----------------------------------------------------------------------------
//
bool QcdBkgdEst::deadECALCut( Event::Data& ev,
			      std::vector<LorentzV>& jets,
			      LorentzV met,
			      std::vector<LorentzV>& baby_jets ) {
  
  double biasedDPhi = 100;

  float deadECALDR = 10;

  LorentzV loweredMHT = met;
  for(std::vector<LorentzV>::const_iterator iM = baby_jets.begin();iM != baby_jets.end();++iM){
    if( (iM)->Pt() > minJetPtCut_)
      loweredMHT -= (*iM);
  }

  for( std::vector<LorentzV>::const_iterator i = jets.begin();
       i != jets.end();
       ++i ){

    float newBiasDPhi = fabs(ROOT::Math::VectorUtil::DeltaPhi(*i,loweredMHT + (*i))) ;

    if(newBiasDPhi < biasedDPhi){
      biasedDPhi = newBiasDPhi;
      if (useFile_) { deadECALDR = DeadECALDR_File(deadECALfile_,(i)->Phi(),(i)->Eta(),nBadCells_); }
      else {
	try{deadECALDR = DeadECALDR( ev,(i)->Phi(),(i)->Eta(),nBadCells_);}
	catch(...){
	  deadECALDR = DeadECALDR_File(deadECALfile_,(i)->Phi(),(i)->Eta(),nBadCells_);
        }
      }
    }
  }

  for( std::vector<LorentzV>::const_iterator iI = baby_jets.begin(); iI != baby_jets.end();
       ++iI) {
    if((iI)->Pt() > minJetPtCut_){
      float newBiasDPhi_2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(*iI, loweredMHT + (*iI) )) ;
      if(newBiasDPhi_2 < biasedDPhi){
	biasedDPhi = newBiasDPhi_2;
	if (useFile_) { deadECALDR = DeadECALDR_File(deadECALfile_,(iI)->Phi(),(iI)->Eta(),nBadCells_); }
	else {
	  try{deadECALDR = DeadECALDR(ev,(iI)->Phi(),(iI)->Eta(),nBadCells_);}
	  catch(...){
	    deadECALDR = DeadECALDR_File(deadECALfile_,(iI)->Phi(),(iI)->Eta(),nBadCells_);
	  } 
	} 
      }
    }
  }
  
  if( biasedDPhi < minBiasDeltaPhiCut_ ){
    if(deadECALDR < deadEcalRegionDrCut_)return false;
  }
  
  return true;
  
}

  // -----------------------------------------------------------------------------
  //
  Double_t QcdBkgdEst::DeadECALDR( Event::Data& ev,float jetphi,float jeteta, int nBad){
    double DRmin = 100;
    TLorentzVector testjetVec(0,0,0,0);
    TLorentzVector testVec(0,0,0,0);
    int idx = 0;
 
    testjetVec.SetPtEtaPhiE(100,jeteta,jetphi,100);
    for(std::vector<PolarLorentzV>::const_iterator ii=ev.ecalDeadTowerTrigPrimP4()->begin();ii!=ev.ecalDeadTowerTrigPrimP4()->end();++ii,++idx){
      if((int)ev.ecalDeadTowerNBadXtals()->at(idx) >= nBad){//hard coded I know
	testVec.SetPtEtaPhiE(100,(*ii).Eta(),(*ii).Phi(),100);
	float DR = fabs(ROOT::Math::VectorUtil::DeltaR(*ii,testjetVec));
	if(DR < DRmin) DRmin = DR;
      }
    
    }

    return DRmin;
  }

  //-----------------------------------------------------------------------------
  //
  Double_t QcdBkgdEst::DeadECALDR_File(TString deadECALfile,float jetphi,float jeteta, int nBad){
    // cout << " hi ho " << deadECALfile << endl;
    TLorentzVector testjetVec(0,0,0,0);
    double DRmin = 100;
    testjetVec.SetPtEtaPhiE(100,jeteta,jetphi,100);

    //cout << " in dead ecal " << endl;
    ifstream myfile(deadECALfile);
    if(!myfile.good()) cout << "file is not there " << endl;
    if(myfile.is_open()){

      while(!myfile.eof()){

	float phi,eta;
	int iEta,iPhi,nBadCells;
	myfile>>eta;
	myfile>>phi;
	myfile>>iEta;
	myfile>>iPhi;
	myfile>>nBadCells;
	// cout <<eta<< " "<< phi<<endl;
	if(nBadCells >= nBad){
	  TLorentzVector testVec(0,0,0,0);
	  testVec.SetPtEtaPhiE(100,eta,phi,100);
	  float DR = fabs(ROOT::Math::VectorUtil::DeltaR(testVec,testjetVec));
	  if(DR < DRmin) DRmin = DR;
	}
      }
   
      myfile.close();
    } 
    return DRmin;

  } 

// -----------------------------------------------------------------------------
//
bool QcdBkgdEst::passVertexSumPtOverHt( Event::Data& ev, double ht ) {

  if( !ev.vertexSumPt.enabled() ) { return true; }
  
  double vertex_pt = 0.;
  std::vector<floatle>::const_iterator ivtx = ev.vertexSumPt()->begin();
  std::vector<floatle>::const_iterator jvtx = ev.vertexSumPt()->end();
  for ( ; ivtx != jvtx; ++ivtx ) { 

    int pos = static_cast<int>( ivtx - ev.vertexSumPt()->begin() );
    
    if( !ev.vertexIsFake()->at( pos ) && 
	ev.vertexNdof()->at( pos ) > 4 &&
	fabs( ev.vertexPosition()->at( pos ).Z() ) < 24.0 && 
	ev.vertexPosition()->at( pos ).Rho() < 2.0 ) { vertex_pt += *ivtx; }
    
  }
  
  if ( vertex_pt / ( ht > 0. ? ht : ev.CommonHT() ) > 0.1 ) { return true; } 
  else { return false; }

}

// -----------------------------------------------------------------------------
//
int QcdBkgdEst::trigger( const Event::Data& ev, const vstring& triggers ) { 

  if ( triggers.empty() ) { return 1; }

  int prescale = -1;

  vstring::const_iterator iname = triggers.begin();
  vstring::const_iterator jname = triggers.end();
  for ( ; iname != jname; ++iname ) {

    // Find trigger 
    const std::map<std::string,bool>* triggers = ev.hlt();
    std::map<std::string,bool>::const_iterator itrigger = triggers->find(*iname); 
    if ( itrigger == triggers->end() ) { continue; }

    // Check if trigger fired
    if ( !itrigger->second ) { continue; }
    
    // Find prescale
    std::map<std::string,int>::const_iterator iprescale = ev.hlt_prescaled()->find(itrigger->first);
    if ( iprescale == ev.hlt_prescaled()->end() ) { continue; }
    
    // Store lowest prescale of triggers that fire
    if ( prescale < 0 || iprescale->second < prescale ) { prescale = iprescale->second; }
    
  }
  
  return prescale;

}

// -----------------------------------------------------------------------------
//
double QcdBkgdEst::minBiasDeltaPhi( const std::vector<LorentzV>& objects ) { 

  // Calculate recoil
  LorentzV recoil = LorentzV(0.,0.,0.,0.);
  std::vector<LorentzV>::const_iterator ii = objects.begin();
  std::vector<LorentzV>::const_iterator jj = objects.end();
  for ( ; ii != jj; ++ii ) { recoil += *ii; }
  recoil.SetPxPyPzE( -recoil.Px(), -recoil.Py(), -recoil.Pz(), recoil.E() );

  // Test each configuration
  double min = 10.;
  ii = objects.begin();
  jj = objects.end();
  for ( ; ii != jj; ++ii ) { 
    recoil += *ii; 
    double test = fabs( ROOT::Math::VectorUtil::DeltaPhi( *ii, recoil ) );
    if ( test < min ) { min = test; }
    recoil -= *ii;
  }

  return min;

}




// 	  // Retrieve trigger and prescale info from event
// 	  int prescale = -1; 
// 	  if ( !triggers_.empty() ) {
// 	    vstring::const_iterator iname = triggers_.begin();
// 	    vstring::const_iterator jname = triggers_.end();
// 	    for ( ; iname != jname; ++iname ) {

// 	      // Find trigger 
// 	      const std::map<std::string,bool>* triggers = ev.hlt();
// 	      std::map<std::string,bool>::const_iterator itrigger = triggers->find(*iname); 
// 	      if ( itrigger == triggers->end() ) { 
// 		//std::cout << " missing name " << *iname << std::endl;
// 		continue; 
// 	      }

// 	      // Find prescale
// 	      std::map<std::string,int>::const_iterator iprescale = ev.hlt_prescaled()->find(itrigger->first);
// 	      if ( iprescale == ev.hlt_prescaled()->end() ) { 
// 		//std::cout << " missing prescale " << *iname << std::endl;
// 		continue; 
// 	      }
	      
// 	      // Extract HT value
// 	      double trigger_ht = -1.;
// 	      size_t pos = itrigger->first.find("HLT_HT");
// 	      std::string substr = itrigger->first.substr(pos+6,3);
// 	      if ( pos != std::string::npos ) { trigger_ht = 1.*atoi( substr.c_str() ); }
// // 	      std::cout << " trigger " << itrigger->first
// // 			<< " pos " << pos 
// // 			<< " substr " << substr
// // 			<< " trigger_ht " << trigger_ht
// // 			<< std::endl;

// 	      // HT monitor trigger names comprise 12 characters
// 	      bool monitor_trigger = ( itrigger->first.size() == 12 );
	      
// 	      // Use only HT monitor (signal) triggers for left (right) side; use unprescaled HT monitor triggers always
// 	      bool use_trigger = ( ( ( at_reco > at_ ) && !monitor_trigger ) || 
// 				   ( ( at_reco < at_ ) &&  monitor_trigger ) );//||
// 	      //( monitor_trigger && iprescale->second == 1 ) );
	      
// // 	      if ( trigger_ht != 250 &&
// // 		   trigger_ht != 260 &&
// // 		   trigger_ht != 300 &&
// // 		   trigger_ht != 350 ) {
// //  	      if ( monitor_trigger ) {
// //  		std::cout << "TEST  "
// // 			  << itrigger->first << " "
// // 			  << itrigger->second << " "
// //  			  << at_reco << " "
// //  			  << monitor_trigger << " "
// //  			  << trigger_ht << " "
// //  			  << use_trigger << " "
// //  			  << iprescale->second << " "
// //  			  << ht_reco << " "
// //  			  << htBins_[ibin] << " "
// //  			  << ibin
// //  			  << std::endl;
// //  	      }

// 	      // HT offset to guarantee plateau
// 	      double offset = 1.e-3;
	      
// 	      // Store lowest prescale of triggers that fire
// 	      if ( itrigger->second && 
// 		   (htBins_[ibin]+offset) > trigger_ht && 
// 		   //(ht_reco+offset) > trigger_ht && 
// 		   use_trigger &&
// 		   ( prescale < 0 || iprescale->second < prescale ) ) { 
// 		prescale = iprescale->second; 
// 	      }

// // 	      std::cout << " trigger " << itrigger->first
// // 			<< " fired " << itrigger->second
// // 			<< " trigger_ht " << trigger_ht
// // 			<< " ht_reco " << ht_reco
// // 			<< " prescale " << iprescale->second
// // 			<< std::endl;

// 	    }
// 	  } else { prescale = 1; }
