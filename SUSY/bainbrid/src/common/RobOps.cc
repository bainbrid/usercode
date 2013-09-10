#include "RobOps.hh"
#include "AlphaT.hh"
#include "Jet.hh"
#include "KinSuite.hh"
#include "Types.hh"
#include <iostream>
#include <string>

using namespace Operation;

// -----------------------------------------------------------------------------
//
RobAlphaT::RobAlphaT( float cut ) 
  : cut_(cut) 
{;}

// -----------------------------------------------------------------------------
//
bool RobAlphaT::Process( Event::Data& ev ) {
  if ( ev.CommonObjects().size() < 2 ||
       ev.CommonObjects().size() > 50 ) { return false; }
  if ( AlphaT()( ev.CommonObjects() ) > cut_ ) { return true; }
  return false;
}

// -----------------------------------------------------------------------------
//
std::ostream& RobAlphaT::Description( std::ostream &ostrm ) {
  ostrm << "RobAlphaT " << cut_ << " " ;
  return ostrm;
}

// -----------------------------------------------------------------------------
//
RobOps::RobOps( const Utils::ParameterSet& ps ) 
  : algo_(""),
    cut_(0.),
    jets_()
{
  if ( ps.Contains("Algo") ) { algo_ = ps.Get<std::string>("Algo"); }
  if ( ps.Contains("Cut") )  { cut_ = ps.Get<double>("Cut"); }
}

// -----------------------------------------------------------------------------
//
RobOps::~RobOps() {}

// -----------------------------------------------------------------------------
//
void RobOps::Start( Event::Data& ev ) {
  jets_.clear();
}

// -----------------------------------------------------------------------------
//
bool RobOps::Process( Event::Data& ev ) {

  uint n = ev.CommonObjects().size();

  if ( algo_ == "PreScale" ) {

    static int count = 0;
    static int freq = int( cut_ + 1.e-6 );
    count++;
    return !(count%freq);
    
  } else if ( algo_ == "GenMet" ) {
    
    if ( ev.genMetP4AK5()->Pt() < cut_ ) { return true; }
    else { return false; }
    
  } else if ( algo_ == "AlphaT" ) {
    
    // AlphaT
    if ( n < 2 || n > 50 ) { return false; }
    if ( AlphaT()( ev.CommonObjects() ) > cut_ ) { 
      if ( n >= jets_.size() ) { jets_.resize( n+1, 0. ); }
      jets_[n] += ev.GetEventWeight();
      return true;
    } else { return false; }
    
  } else if ( algo_ == "HadronicAlphaT" ) {

    std::vector<LorentzV> jets = ev.HadronicObjects();
    uint njets = jets.size();
    if ( njets < 2 || njets > 50 ) { return false; }
    if ( AlphaT()( jets ) > cut_ ) { return true; }
    else { return false; }
    
  } else if ( algo_ == "InvariantMass" ) {
    
    // Invariant mass
    
    LorentzV mht = ev.CommonRecoilMET();

    std::vector<bool> pseudo;
    AlphaT()( ev.CommonObjects(), pseudo );
    if ( pseudo.size() != n ) { abort(); }
  
    LorentzV lv1(0.,0.,0.,0.);
    LorentzV lv2(0.,0.,0.,0.);
    if ( n == 2 ) {
      lv1 = ev.CommonObjects()[0];
      lv2 = ev.CommonObjects()[1];
    } else if ( n > 2 ) {
      for ( unsigned int i = 0; i < n; ++i ) {
	if ( pseudo[i] ) { lv1 += ev.CommonObjects()[i]; }
	else             { lv2 += ev.CommonObjects()[i]; }
      }
      if ( lv2.Pt() > lv1.Pt() ) { LorentzV tmp = lv1; lv1 = lv2; lv2 = tmp; }
    }
    
    return ( lv1 + lv2 ).M() < cut_;

  } else if ( algo_ == "Dalitz" ) {
    
    // Dalitz 
    LorentzV mht;
    LorentzV lv1;
    LorentzV lv2;
    std::pair<double,double> rho = RobOps::dalitz( ev, mht, lv1, lv2, false );
    double rho_a = rho.first;
    double rho_b = rho.second;
    if ( rho_a < 0. || rho_b < 0. ) { return false; }
    double sigma = ( rho_a * rho_a + rho_b * rho_b ) / ( ( 1 - 2 * rho_a ) * ( 1 - 2 * rho_b ) );
    if ( sigma < cut_ ) { 
      if ( n >= jets_.size() ) { jets_.resize( n+1, 0. ); }
      jets_[n] += ev.GetEventWeight();
      
//       std::cout<< " event: " << ev.EventNumber()
// 	       << " mht: " << mht.Pt()
// 		<< " pt1: " << lv1.Pt()
// 		<< " pt2: " << lv2.Pt()
// 		<< " rho_a: " << rho_a
// 		<< " rho_b: " << rho_b
// 		<< " sigma: " << sigma
// 		<< std::endl;
      
      return true; 
    } else { return false; }

  } else if ( algo_ == "BiasedDPhi" ) {
    
    return ( ev.BiasedCommonRecoilMETJetDPhi() > cut_ );
    
  } else {
    std::cout << "[RobOps::Process]"
	      << " Unknown Algorithm!..."
	      << std::endl;
  }

  return false; 
  
}

// -----------------------------------------------------------------------------
//
std::pair<double,double> RobOps::dalitz( const Event::Data& ev,
					      LorentzV& mht,
					      LorentzV& lv1,
					      LorentzV& lv2,
					      bool thrust ) {

  std::vector<LorentzV> common = ev.CommonObjects();
  //sort( common.begin(), common.end(), KinSuite::Compare );

  //mht = ev.CommonMHT();
  mht = ev.CommonRecoilMET();
  lv1 = LorentzV(0.,0.,0.,0.);
  lv2 = LorentzV(0.,0.,0.,0.);

  if ( common.size() < 2 ) { 
    
    mht = LorentzV(0.,0.,0.,0.);
    lv1 = LorentzV(0.,0.,0.,0.);
    lv2 = LorentzV(0.,0.,0.,0.);
    return std::make_pair( -1., -1. ); 

//   } else if ( common.size() == 1 ) { 

//     lv1 = common[0];
//     lv2 = LorentzV(0.,0.,0.,0.);
    
//   } else if ( common.size() == 2 ) { 
    
//     lv1 = common[0];
//     lv2 = common[1];
    
//   } else if ( common.size() == 3 ) { 

//     if (0) {
      
//       lv1 = common[0];
//       lv2 = common[1] + common[2];
      
//     } else { 
      
//       std::vector<bool> pseudo;
//       AlphaT()( common, pseudo );
//       if ( pseudo.size() != common.size() ) { abort(); }
      
//       for ( unsigned int i = 0; i < common.size(); ++i ) {
// 	if ( pseudo[i] ) { lv1 += common[i]; }
// 	else             { lv2 += common[i]; }
//       }
      
//     }
    
  } else {
      
    if ( thrust ) {

      ThrustStuff t( common );
      lv1 = t.Pjet1;
      lv2 = t.Pjet2;
	
    } else {
	
      std::vector<bool> pseudo;
      /* double alpha_t = */ AlphaT()( common, pseudo );
      if ( pseudo.size() != common.size() ) { abort(); }

//       std::stringstream ss;
//        ss << " event: " << ev.EventNumber()
// 	 << " alpha_t: " << alpha_t
// 	 << " common.size(): " << common.size()
// 	 << " pseudo.size(): " << pseudo.size() << std::endl;
      for ( unsigned int i = 0; i < common.size(); ++i ) {
	if ( pseudo[i] ) { lv1 = lv1 + common[i]; }
	else             { lv2 = lv2 + common[i]; }
// 	ss << " i: " << i 
// 	   << " pseudo[i]: " << pseudo[i]
// 	   << " common[i].Phi: " << common[i].Phi()
// 	   << " common[i].Pt: " << common[i].Pt()
// 	   << " lv1.Pt(): " << lv1.Pt()
// 	   << " lv2.Pt(): " << lv2.Pt() << std::endl;
      }
      //std::cout << ss.str();
    }
      
  }

  // Order by pT
  if ( lv2.Pt() > lv1.Pt() ) { 
    LorentzV tmp = lv1; 
    lv1 = lv2; 
    lv2 = tmp; 
  }
  
  // Dalitz variables
  double rho_a = lv2.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );
  double rho_b = lv1.Pt() / ( lv1.Pt() + lv2.Pt() + mht.Pt() );
  
  return std::make_pair( rho_a, rho_b );
    
}

// -----------------------------------------------------------------------------
//
std::vector<LorentzV> RobOps::genJets( const Event::Data& ev ) {
  
  static double min_pt = Utils::GetConfig<double>("Common.Jets.PtCut");
  static double max_eta = Utils::GetConfig<double>("Common.Jets.EtaCut");
  
  static bool once = true;
  if (once) {
    once = false;
    std::cout << "[RobOps::genJets]"
  	      << " GenJets.MinPt: " << min_pt 
  	      << " GenJets.MaxEta: " <<max_eta
  	      << std::endl;
  }
  
  std::vector<LorentzV> jets;
  std::vector<LorentzV>::const_iterator ijet = ev.genJetP4()->begin();
  std::vector<LorentzV>::const_iterator jjet = ev.genJetP4()->end();
  for ( ; ijet != jjet; ++ijet ) {
    if ( ijet->Pt() > min_pt && fabs(ijet->Eta()) < max_eta ) { 
      jets.push_back( *ijet ); 
    }
  }    
  
  return jets;
  
}

// -----------------------------------------------------------------------------
//
double RobOps::ht( const std::vector<LorentzV>& jets ) {
  double ht = 0.;
  std::vector<LorentzV>::const_iterator ijet = jets.begin();
  std::vector<LorentzV>::const_iterator jjet = jets.end();
  for ( ; ijet != jjet; ++ijet ) { ht += ijet->Et(); } 
  return ht;
}

// -----------------------------------------------------------------------------
//
LorentzV RobOps::mht( const std::vector<LorentzV>& jets ) {
  LorentzV recoil(0.,0.,0.,0.);
  std::vector<LorentzV>::const_iterator ijet = jets.begin();
  std::vector<LorentzV>::const_iterator jjet = jets.end();
  for ( ; ijet != jjet; ++ijet ) { recoil -= *ijet; } 
  return recoil;
}

// -----------------------------------------------------------------------------
//
double RobOps::meff( const std::vector<LorentzV>& jets ) {
  LorentzV recoil(0.,0.,0.,0.);
  double sum_et = 0.; //@@ should be HT???
  std::vector<LorentzV>::const_iterator ijet = jets.begin();
  std::vector<LorentzV>::const_iterator jjet = jets.end();
  for ( ; ijet != jjet; ++ijet ) { 
    recoil -= *ijet;
    sum_et += ijet->Et(); //@@ should Pt???
  }
  return sum_et + recoil.Pt();
}

// -----------------------------------------------------------------------------
//
std::ostream& RobOps::Description( std::ostream& ss ) {
  ss << "RobOps is \"" << algo_
     << "\" with cut value " << cut_ 
     << " " << std::endl;
//   double total = 0.;
//   std::vector<double>::const_iterator ii = jets_.begin();
//   std::vector<double>::const_iterator jj = jets_.end();
//   for ( ; ii != jj; ++ii ) { total += *ii; }
//   ss << " Total=" << total << std::endl;
  std::vector<double>::const_iterator iii = jets_.begin();
  std::vector<double>::const_iterator jjj = jets_.end();
  for ( ; iii != jjj; ++iii ) { 
    ss << " nJets=" << int( iii - jets_.begin() ) << ": " << *iii << std::endl;
  }
  ss << " Total   ";
  return ss;
}
 
