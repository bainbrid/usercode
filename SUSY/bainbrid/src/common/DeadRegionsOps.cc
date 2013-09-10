#include "DeadRegionsOps.hh"
#include "EventData.hh"
#include "GenMatrixBin.hh"
#include "Math/PtEtaPhiE4D.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"
#include <sstream>
#include <string>
#include <iomanip>

using namespace Operation;

// -----------------------------------------------------------------------------
//
DeadRegionsOps::DeadRegionsOps( const Utils::ParameterSet& ps ) :
  dead_(),
  minDPhi_(0.),
  respMax_(0.),
  veto_(false),
  problem_(false),
  verbose_(false),
  histos_(false),
  dirName_("DeadRegions"),
  nMin_(2),
  nMax_(6),
  hDPhiNewVsOld_(),
  hRespGen_(),
  hRespProjVsGen_(),
  hRespRecoilVsGen_(),
  hRespProj_(),
  hRespRecoil_(),
  hRespRecoilVsProj_(),
  hRespProjVsDPhi_(),
  hRespRecoilVsDPhi_(),
  hDEta_(),
  hDPhi_(),
  hDR_(),
  hGenEtaVsGenPhi_(),
  hGenEtaVsGenPhiAdj_()
{

  std::vector<double> eta = ps.Get< std::vector<double> >("Eta");
  std::vector<double> phi = ps.Get< std::vector<double> >("Phi");
  if ( eta.size() != phi.size() ) {
    std::cout << "[DeadRegionsOps::DeadRegionsOps]"
	      << " WARNING! Eta and Phi vectors differ in size!"
	      << std::endl;
    problem_ = true;
  } else {
    dead_.clear();
    for ( uint ii = 0; ii < eta.size(); ++ii ) {
      dead_.push_back( std::make_pair(eta[ii],phi[ii]) );
    }
  }

  if ( ps.Contains("MinDPhi") ) minDPhi_ = ps.Get<double>("MinDPhi");
  if ( ps.Contains("ResponseMax") ) respMax_ = ps.Get<double>("ResponseMax");
  if ( ps.Contains("Veto") ) veto_ = ps.Get<bool>("Veto");
  if ( ps.Contains("Verbose") ) verbose_ = ps.Get<bool>("Verbose");
  if ( ps.Contains("Histos") ) histos_ = ps.Get<bool>("Histos");
  if ( ps.Contains("DirName") ) dirName_ = ps.Get<std::string>("DirName");
  if ( ps.Contains("MinObjects") ) nMin_ = ps.Get<int>("MinObjects");
  if ( ps.Contains("MaxObjects") ) nMax_ = ps.Get<int>("MaxObjects");
  
}

// -----------------------------------------------------------------------------
//
DeadRegionsOps::~DeadRegionsOps() {}

// -----------------------------------------------------------------------------
//
void DeadRegionsOps::Start( Event::Data& ev ) {
  initDir( ev.OutputFile(), dirName_.c_str() );
  BookHistos();
}

// -----------------------------------------------------------------------------
//
void DeadRegionsOps::BookHistos() {

  if ( !histos_ ) { return; }

  BookHistArray( hDPhiNewVsOld_,
		 "hDPhiNewVsOld",
		 ";Old;New",
		 320, 
		 0., 
		 3.2, 
		 320, 
		 0., 
		 3.2, 
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespGen_,
		 "hRespGen",
		 ";p_{T}^{jet}/p_{T}^{gen}",
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespProjVsGen_,
		 "hRespProjVsGen",
		 ";p_{T}^{jet}/p_{T}^{proj};p_{T}^{jet}/p_{T}^{gen}",
		 500,
		 -5.,
		 5.,
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespRecoilVsGen_,
		 "hRespRecoilVsGen",
		 ";p_{T}^{jet}/p_{T}^{recoil};p_{T}^{jet}/p_{T}^{gen}",
		 500,
		 -5.,
		 5.,
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespProj_,
		 "hRespProj",
		 ";p_{T}^{jet}/p_{T}^{proj}",
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespRecoil_,
		 "hRespRecoil",
		 ";p_{T}^{jet}/p_{T}^{recoil}",
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hRespRecoilVsProj_,
		 "hRespRecoilVsProj",
		 ";p_{T}^{jet}/p_{T}^{proj};p_{T}^{jet}/p_{T}^{recoil}",
		 500,
		 -5.,
		 5.,
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);
  
  BookHistArray( hRespProjVsDPhi_,
		 "hRespProjVsDPhi",
		 ";min_{i}#Delta#phi;p_{T}^{jet}/p_{T}^{proj}",
		 320, 
		 0., 
		 3.2, 
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);
  
  BookHistArray( hRespRecoilVsDPhi_,
		 "hRespRecoilVsDPhi",
		 ";min_{i}#Delta#phi;p_{T}^{jet}/p_{T}^{recoil}",
		 320, 
		 0., 
		 3.2, 
		 500,
		 -5.,
		 5.,
		 nMax_+1, 0, 1, true);

  BookHistArray( hDEta_,
		 "DEta",
		 ";#Delta#eta_{min}(j_{i},dead);", 
		 100, 
		 0., 
		 1., 
		 nMax_+1, 0, 1, true);

  BookHistArray( hDPhi_,
		 "DPhi",
		 ";#Delta#phi_{min}(j_{i},dead);", 
		 100, 
		 0., 
		 1., 
		 nMax_+1, 0, 1, true);

  BookHistArray( hDR_,
		 "DR",
		 ";#DeltaR_{min}(j_{i},dead);", 
		 100, 
		 0., 
		 1., 
		 nMax_+1, 0, 1, true);

  BookHistArray( hGenEtaVsGenPhi_,
		 "GenEtaVsGenPhi",
		 ";#phi;#eta",
		 64,
		 -3.2,
		 3.2,
		 75, 
		 -1.5, 
		 1.5, 
		 nMax_+1, 0, 1, true);

  BookHistArray( hGenEtaVsGenPhiAdj_,
		 "GenEtaVsGenPhiAdjusted",
		 ";#phi (adjusted);#eta",
		 64,
		 0.,
		 6.4,
		 75, 
		 -1.5, 
		 1.5, 
		 nMax_+1, 0, 1, true);
  
}

// -----------------------------------------------------------------------------
//
bool DeadRegionsOps::Process( Event::Data& ev ) {
  
  // Check map
  if ( problem_ ) { return false; }

  badlyMeasuredGenJets(ev);

  // Retrieve jets
  Jets jets = ev.JD_CommonJets().accepted;

  // MinBiasedDPhi and Response
  double biased_dphi = -1.;
  double response = 1.; 
  double response1 = 1.; 
  Jets::const_iterator jet = minBiasedDeltaPhi( jets, biased_dphi, response, response1 );

  // Check valid jet is returned
  if ( jet == jets.end() ) { return true; }
  
  // Retrieve GenJet P4
  LorentzV gen = matchedGenJet( ev, jet );

  // Find closest dead region for given jet
  double deta = -1.;
  double dphi = -1.;
  double dr = closestDeadRegion( jet, deta, dphi );
    
  // Plots

  if ( histos_ && jet != jets.end() ) {
    
    // Event weight
    Double_t weight = ev.GetEventWeight();
    
    // Retrieve Jet index in ntuple
    int index = (*jet)->GetIndex(); 
    
    // BiasedDPhi, new against old
    double check = ev.BiasedCommonRecoilMETJetDPhi();
    if ( !hDPhiNewVsOld_.empty() ) { hDPhiNewVsOld_[0]->Fill( check, biased_dphi, weight ); }
    if ( index+1 < int(hDPhiNewVsOld_.size()) ) { hDPhiNewVsOld_[index+1]->Fill( check, biased_dphi, weight ); }
    
    if ( gen.Pt() > 0. ) {
      
      // Responses
      
      double resp_corr = gen.Pt() > 0. ? (*jet)->Pt() / gen.Pt() : -1.;
      
      if ( !hRespGen_.empty() ) { hRespGen_[0]->Fill( resp_corr, weight ); }
      if ( index+1 < int(hRespGen_.size()) ) { hRespGen_[index+1]->Fill( resp_corr, weight ); }
      
      if ( !hRespProjVsGen_.empty() ) { hRespProjVsGen_[0]->Fill( resp_corr, response, weight ); }
      if ( index+1 < int(hRespProjVsGen_.size()) ) { hRespProjVsGen_[index+1]->Fill( resp_corr, response, weight ); }
      
      if ( !hRespRecoilVsGen_.empty() ) { hRespRecoilVsGen_[0]->Fill( resp_corr, response1, weight ); }
      if ( index+1 < int(hRespRecoilVsGen_.size()) ) { hRespRecoilVsGen_[index+1]->Fill( resp_corr, response1, weight ); }
      
    }
    
    // Responses
    
    if ( !hRespProj_.empty() ) { hRespProj_[0]->Fill( response, weight ); }
    if ( index+1 < int(hRespProj_.size()) ) { hRespProj_[index+1]->Fill( response, weight ); }
    
    if ( !hRespRecoil_.empty() ) { hRespRecoil_[0]->Fill( response1, weight ); }
    if ( index+1 < int(hRespRecoil_.size()) ) { hRespRecoil_[index+1]->Fill( response1, weight ); }
    
    if ( !hRespRecoilVsProj_.empty() ) { hRespRecoilVsProj_[0]->Fill( response, response1, weight ); }
    if ( index+1 < int(hRespRecoilVsProj_.size()) ) { hRespRecoilVsProj_[index+1]->Fill( response, response1, weight ); }

    // Response vs biased DPhi
    
    if ( !hRespProjVsDPhi_.empty() ) { hRespProjVsDPhi_[0]->Fill( biased_dphi, response, weight ); }
    if ( index+1 < int(hRespProjVsDPhi_.size()) ) { hRespProjVsDPhi_[index+1]->Fill( biased_dphi, response, weight ); }
    
    if ( !hRespRecoilVsDPhi_.empty() ) { hRespRecoilVsDPhi_[0]->Fill( biased_dphi, response1, weight ); }
    if ( index+1 < int(hRespRecoilVsDPhi_.size()) ) { hRespRecoilVsDPhi_[index+1]->Fill( biased_dphi, response1, weight ); }
    
    // Closest dead region
    
    if ( dphi < minDPhi_ ) { 
      
      if ( dr >= 0. ) {
	
	if ( !hDEta_.empty() ) { hDEta_[0]->Fill( deta, weight ); }
	if ( index+1 < int(hDEta_.size()) ) { hDEta_[index+1]->Fill( deta, weight ); }
	
	if ( !hDPhi_.empty() ) { hDPhi_[0]->Fill( dphi, weight ); }
	if ( index+1 < int(hDPhi_.size()) ) { hDPhi_[index+1]->Fill( dphi, weight ); }
	
	if ( !hDR_.empty() ) { hDR_[0]->Fill( dr, weight ); }
	if ( index+1 < int(hDR_.size()) ) { hDR_[index+1]->Fill( dr, weight ); }
	
      }
      
    }
    
  }
  
  // Veto event 
  if ( !problem_ && veto_ && 
       jet != jets.end() && dphi < minDPhi_ && 
       response > 0. && response < respMax_ ) { return false; }
  else { return true; }
  
}

// -----------------------------------------------------------------------------
//
std::ostream& DeadRegionsOps::Description( std::ostream& ostrm ) {
  ostrm << " Dead regions ";
  return ostrm;
}

// -----------------------------------------------------------------------------
//
void DeadRegionsOps::badlyMeasuredGenJets( Event::Data& ev ) {

  // Event weight
  Double_t weight = ev.GetEventWeight();

  // Retrieve jets
  Jets jets = ev.JD_CommonJets().accepted;

  // "Badly measured" GenJets 
  bool found = false;
  ICF_LorentzVs::const_iterator igen = ev.genJetP4()->begin();
  ICF_LorentzVs::const_iterator jgen = ev.genJetP4()->end();
  for ( ; igen != jgen; ++igen ) {
    if ( igen->Pt() < 50 ) { continue; }
    found = false;
    Jets::const_iterator ijet = jets.begin();
    Jets::const_iterator jjet = jets.end();
    for ( ; ijet != jjet; ++ijet ) {
      if ( !(*ijet) ) { continue; }
      double deta = (*ijet)->Eta() - igen->Eta();
      double dphi = ROOT::Math::VectorUtil::DeltaPhi( *igen, **ijet );
      double dr = sqrt( deta * deta + dphi * dphi );
      if ( dr < 0.7 && (*ijet)->Pt() / igen->Pt() > 0.5 ) { 
	found = true; 
	break;
      }
    }
    if ( !found ) { break; }
  }

  // Plots of eta/phi of "badly measured" GenJets
  if ( histos_ && !found && igen != ev.genJetP4()->end() ) {
    int index = int( igen - ev.genJetP4()->begin() );
    if ( !hGenEtaVsGenPhi_.empty() ) { 
      hGenEtaVsGenPhi_[0]->Fill( igen->Phi(), igen->Eta(), weight ); 
    }
    if ( index+1 < int(hGenEtaVsGenPhi_.size()) ) { 
      hGenEtaVsGenPhi_[index+1]->Fill( igen->Phi(), igen->Eta(), weight ); 
    }
    double phi = igen->Phi() < 0. ? igen->Phi() + 2 * 3.141592 : igen->Phi();
    if ( !hGenEtaVsGenPhiAdj_.empty() ) { 
      hGenEtaVsGenPhiAdj_[0]->Fill( phi, igen->Eta(), weight ); 
    }
    if ( index+1 < int(hGenEtaVsGenPhiAdj_.size()) ) { 
      hGenEtaVsGenPhiAdj_[index+1]->Fill( phi, igen->Eta(), weight ); 
    }
  }

}
  
// -----------------------------------------------------------------------------
//
DeadRegionsOps::Jets::const_iterator DeadRegionsOps::minBiasedDeltaPhi( Jets& jets,
									double& min_biased_dphi, 
									double& response,
									double& response1 ) {
  
  // Calculate MHT for jet system
  LorentzV mht(0.,0.,0.,0.);
  Jets::const_iterator ijet = jets.begin();
  Jets::const_iterator jjet = jets.end();
  for ( ; ijet != jjet; ++ijet ) {
    if ( !(*ijet) ) { continue; }
    mht += **ijet;
  }
  mht.SetPx( -1.*mht.Px() );
  mht.SetPy( -1.*mht.Py() );
  mht.SetPz( -1.*mht.Pz() );
  
  // Counters
  Jets::const_iterator jet = jets.end();
  double dphi = 10.;
  double resp = -1.;
  double resp1 = -1.;
  
  // Calculate recoil and response for given jet
  Jets::const_iterator ii = jets.begin();
  Jets::const_iterator jj = jets.end();
  for ( ; ii != jj; ++ii ) {
    if ( !(*ii) ) { continue; }
    LorentzV recoil = mht + **ii; 
    LorentzV met( recoil.Px(), recoil.Py(), 0., 0. );
    LorentzV dir( (*ii)->Px(), (*ii)->Py(), 0., 0. );
    double proj_jet = dir.Dot(met) / dir.Dot(dir);
    double resp_pt = (*ii)->Pt() / met.Pt(); 
    double resp_proj = 1. / proj_jet;
    double diff = fabs( ROOT::Math::VectorUtil::DeltaPhi(dir,met) );
    if ( diff < dphi ) { 
      jet = ii;
      dphi = diff;
      resp = resp_proj;
      resp1 = resp_pt;
    }
    
    // Debug
    if ( false && dphi < 0.0001 ) {
      std::cout << " Jet: pt: " << dir.Pt()
		<< " phi: " << dir.Phi()
		<< " Met: pt: " << met.Pt()
		<< " phi: " << met.Phi()
		<< " dphi: " << dphi
		<< " proj_jet: " << proj_jet
		<< " resp_pt: " << resp_pt
		<< " resp_proj: " << resp_proj
		<< std::endl;
    }

  }
  
  // Return values
  if ( jet != jets.end() ) { 
    min_biased_dphi = dphi;
    response = resp;
    response1 = resp1;
    return jet;
  } else {
    return jets.end();
  }

}

// -----------------------------------------------------------------------------
//
LorentzV DeadRegionsOps::matchedGenJet( Event::Data& ev,
					Jets::const_iterator jet ) {
  
  // GenJet P4
  LorentzV p4(0.,0.,0.,0.);
  
  // Retrieve Jet index in ntuple
  int index = (*jet)->GetIndex(); 
  
  // Retrieve GenJet index
  if ( index < int( ev.genJetMatchIndex()->size() ) ) {
    index = ev.genJetMatchIndex()->at(index);
  } else { index = -1; }
  
  // Check if index is valid for GenJet collection
  if ( index >= 0 && index < int(ev.genJetP4()->size()) ) { 
    
    // Retrieve LV for matched GenJet
    p4 = ev.genJetP4()->at(index);
    
  }
  
  return p4;
  
}

// -----------------------------------------------------------------------------
//
double DeadRegionsOps::closestDeadRegion( Jets::const_iterator jet,
					  double& deta,
					  double& dphi ) { 
  
  // Counters
  double dr_min = 1.e6;
  double deta_min = 0.;
  double dphi_min = 0.;
  
  // Find nearest dead region for given jet
  Regions::const_iterator ii = dead_.begin();
  Regions::const_iterator jj = dead_.end();
  for ( ; ii != jj; ++ii ) {
    ROOT::Math::PtEtaPhiE4D<double> lv( 1., ii->first, ii->second, 1. );
    double deta = (*jet)->Eta() - lv.Eta();
    double dphi = ROOT::Math::VectorUtil::DeltaPhi( **jet, lv );
    double dr = sqrt( deta * deta + dphi *dphi );
    if ( dr < dr_min ) { 
      dr_min = dr; 
      deta_min = deta;
      dphi_min = dphi;
    }
  }

  // Return values
  if ( dr_min < 1.e5 ) {
    deta = deta_min;
    dphi = dphi_min;
    return dr_min;
  } else {
    return -1.;
  }

}
