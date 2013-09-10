#include "Ops.hh"
#include "AlphaT.hh"
#include <iostream>
#include <string>
#include <fstream>
#include "vector"
#include <cmath>
#include <iomanip>
#include "Types.hh"

using namespace Operation;

// ----------------------------------------------------------------
<<<<<<< .mine
AlphatTriggerCut::AlphatTriggerCut( float cut, float HTthresh, float JetThreshold, float MhtJetThreshold = -1. )
=======
AlphatTriggerCut::AlphatTriggerCut( float cut, float HTthresh, float LowerBound)
>>>>>>> .r2245
  : cut_(cut),
  HTthresh_(HTthresh),
  LowerBound_(LowerBound)
<<<<<<< .mine
{
  if ( MhtThreshold_ < 0. ) { MhtThreshold_ = JetThreshold_; }
}
=======
  {;}
// ----------------------------------------------------------------
bool AlphatTriggerCut::Process( Event::Data& ev){

  double MHTx = 0.;
  double MHTy = 0.;
  double HT = 0.;
  for(unsigned int i = 0; i <  ev.JD_Jets().size(); i++){
    if(ev.JD_Jets()[i].Et() > LowerBound_ && fabs(ev.JD_Jets()[i].Eta()) < 3.){
      double  MHT = 0.;
      HT+=ev.JD_Jets()[i].Et();
      MHTx-=ev.JD_Jets()[i].Et()*cos(ev.JD_Jets()[i].Phi());
      MHTy-=ev.JD_Jets()[i].Et()*sin(ev.JD_Jets()[i].Phi());
      MHT = sqrt(MHTy*MHTy + MHTx*MHTx);
      if(i == 1){ //dijet case check alpha T
        double alpha = (HT - (ev.JD_Jets()[0].Et() - ev.JD_Jets()[1].Et() )) / (2.0*sqrt(HT*HT - MHT*MHT));
        if(HT > HTthresh_ && alpha > cut_){ return true;}
      } // end Dijets
      if(i >1 ){ // nJet case check MHT/HT (BetaT)
        double alpha = (HT) / (2.0*sqrt(HT*HT - MHT*MHT));
        if(HT > HTthresh_ && alpha > cut_){return true;}
      } // end NJets
    }//end jet selection
  }//end jet loop
>>>>>>> .r2245

<<<<<<< .mine
// ----------------------------------------------------------------
bool AlphatTriggerCut::Process( Event::Data& ev ) {

//   bool pass = false;
  
//   int nJets = 0;
//   int mJets = 0;

//   double aT = 0.;
//   double bT = 0.;
//   double DHT = 0.;
//   double MHT = 0.;
//   double MHTx = 0.;
//   double MHTy = 0.;
//   double HT = 0.;
//   double HT1 = 0.;
  
//   std::stringstream ss;
//   ss << "[AlphatTriggerCut::Process]" << std::endl;

//   ss << " Number of LD muons: " << ev.LD_Muons().size() << std::endl;
//   for ( unsigned int i = 0; i < ev.LD_Muons().size(); ++i ) {
//     Event::Lepton mu( ev.LD_Muons()[i] );
//     if ( mu.Et() > 10. ) {
//       ss << " i: " << i 
// 	 << " ET: " << mu.Et()
// 	 << " eta: " << mu.Eta()
// 	 << " loose: " << mu.GetLooseId()
// 	 << " tight: " << mu.GetTightId()
// 	 << std::endl;
//     } else { continue; }
//   }

//   ss << " Number of LD electrons: " << ev.LD_Electrons().size() << std::endl;
//   for ( unsigned int i = 0; i < ev.LD_Electrons().size(); ++i ) {
//     Event::Lepton ele( ev.LD_Electrons()[i] );
//     if ( ele.Et() > 10. ) {
//       ss << " i: " << i 
// 	 << " ET: " << ele.Et()
// 	 << " eta: " << ele.Eta()
// 	 << " loose: " << ele.GetLooseId()
// 	 << " tight: " << ele.GetTightId()
// 	 << std::endl;
//     } else { continue; }
//   }

//   ss << " Number of LD photons: " << ev.PD_Photons().size() << std::endl;
//   for ( unsigned int i = 0; i < ev.PD_Photons().size(); ++i ) {
//     Event::Photon phot( ev.PD_Photons()[i] );
//     if ( phot.Et() > 10. ) {
//       ss << " i: " << i 
// 	 << " ET: " << phot.Et()
// 	 << " eta: " << phot.Eta()
// 	 << " loose: " << phot.IsItLoose()
// 	 << " tight: " << phot.IsItTight()
// 	 << std::endl;
//     } else { continue; }
//   }

//   ss << " Number of JD jets: " << ev.JD_Jets().size() << std::endl;
//   for ( unsigned int i = 0; i < ev.JD_Jets().size(); ++i ) {
    
//     Event::Jet jet( ev.JD_Jets()[i] );

//     if ( jet.Et() > 10. ) {
//       ss << " i: " << i 
// 	 << " ET: " << jet.Et()
// 	 << " eta: " << jet.Eta()
// 	 << " corr: " << ev.jetCorrFactor()->at(i);
//     } else { continue; }

//     // HT calc
//     if ( jet.Et() > JetThreshold_ ) {
//       HT1 += jet.Et();
//       mJets++;
//       ss << " SELECTED FOR HT!";
//     }
    
//     // MHT and HT calc
//     if ( jet.Et() > MhtThreshold_ && fabs(jet.Eta()) < 3. ) {
//       HT += jet.Et();
//       MHTx -= jet.Px();
//       MHTy -= jet.Py();
//       DHT += pow(-1.,(int)i) * jet.Et(); //@@ only to be used in dijet case!
//       nJets++;
//       ss << " SELECTED FOR MHT!";
//     }

//     ss << std::endl;
    
//   }
//   MHT = sqrt( MHTx*MHTx + MHTy*MHTy );

//   if ( nJets == 2 ) {
//     aT = ( HT - DHT ) / ( 2. * sqrt( ( HT*HT ) - ( MHT*MHT ) ) );
//     if ( ( HT1 >= HTthresh_ ) && ( aT >= cut_ ) ) { pass = true; }
//   }
//   else if ( nJets > 2 ) {
//     bT = sqrt( 1. / ( 4. * ( 1. - ( (MHT/HT) * (MHT/HT) ) ) ) );
//     if ( ( HT1 >= HTthresh_ ) && ( bT >= cut_ ) ) { pass = true; }
//   }
  
//   ss << " Number of jets selected for HT: " << mJets
//      << std::endl
//      << " Number of jets selected for MHT: " << nJets
//      << std::endl
//      << " HT: " << HT1
//      << std::endl
//      << " HT: " << HT << " MHT: " << MHT << " DHT: " << DHT << " MHT/HT: " << MHT/HT
//      << std::endl
//      << " aT: " << aT << " bT: " << bT
//      << " pass: " << std::boolalpha << pass
//      << std::endl;
  
  //std::cout << ss.str(); 

//   return pass;

  // Some initialization
  int njets = 0;
  LorentzV mht;
  double ht = 0.;
  double dht = 0.;
  for ( unsigned int i = 0; i < ev.JD_Jets().size(); ++i ) {
    
    // Check jet is above minimum pT threshold
    LorentzV jet( ev.JD_Jets()[i] );
    if ( jet.Et() < JetThreshold_ ) { continue; }
    
    // variables used in AlphaT calculation
    njets++;
    ht += jet.Et();
    mht -= jet;
    dht += ( njets < 2 ? jet.Et() : -1.* jet.Et() ); //@@ only use for njets < 4
    
    // Reject monojet events
    if ( njets < 2 ) { continue; }
    
    // Calc AlphaT value
    double aT = 0.;
    if ( njets == 2 || njets == 3 ) {
      aT = ( ht - fabs(dht) ) / ( 2. * sqrt( ( ht*ht ) - ( mht.Pt()*mht.Pt() ) ) );
    } else if ( njets > 3 ) {
      aT = 1. / ( 2. * sqrt( ( 1. - ( (mht.Pt()/ht) * (mht.Pt()/ht) ) ) ) );
    }
    
    // Check if above HT and AlphaT thresholds
    if ( ( ht >= HTthresh_ ) && ( aT >= cut_ ) ) { return true; } 

  }

=======
>>>>>>> .r2245
  return false;

}

std::ostream& AlphatTriggerCut::Description( std::ostream &ostrm ) {
  ostrm << "AlphatTriggerCut " << cut_ << " " << " Jet Scale: " <<LowerBound_ ;
  return ostrm;
}



MeffTriggerCut::MeffTriggerCut( float cut, float JetThreshold, float MhtJetThreshold)
  : cut_(cut),
  JetThreshold_(JetThreshold),
  MhtThreshold_(MhtJetThreshold)
  {;}
bool MeffTriggerCut::Process( Event::Data& ev){
  double MHTx = 0.;
  double MHTy = 0.;
  double  HT =0.;

  vector<Event::Jet> newJets;
  newJets.clear();
  for(unsigned int i = 0; i <  ev.JD_Jets().size(); i++){
    if(ev.JD_Jets()[i].Et() > JetThreshold_ && fabs(ev.JD_Jets()[i].Eta()) < 3.){
      newJets.push_back(ev.JD_Jets()[i]);
    }
  }
  for(unsigned int j = 0; j < newJets.size(); j++){
    if(newJets[j].Et() >= JetThreshold_){HT += newJets[j].Et();}
    if(newJets[j].Et() >= MhtThreshold_){
      MHTx-=newJets[j].Et()*cos(newJets[j].Phi());
      MHTy-=newJets[j].Et()*sin(newJets[j].Phi());
    }
  }
  double  MHT = sqrt(MHTx*MHTx + MHTy*MHTy);

  if(MHT + HT  > cut_){return true;}

  return false;

}


std::ostream& MeffTriggerCut::Description( std::ostream &ostrm ) {
  ostrm << "MeffTriggerCut " << cut_ << " " << " Jet Scale: " <<JetThreshold_ ;
  return ostrm;
}



HtTriggerCut::HtTriggerCut( float cut, float JetThreshold)
  : cut_(cut),
  JetThreshold_(JetThreshold)
  {;}
bool HtTriggerCut::Process( Event::Data& ev){
  double  HT =0.;

  vector<Event::Jet> newJets;
  newJets.clear();
  for(unsigned int i = 0; i <  ev.JD_Jets().size(); i++){
    if(ev.JD_Jets()[i].Et() > JetThreshold_ && fabs(ev.JD_Jets()[i].Eta()) < 3.){
      newJets.push_back(ev.JD_Jets()[i]);
    }
  }
  for(unsigned int j = 0; j < newJets.size(); j++){
    if(newJets[j].Et() >= JetThreshold_){HT += newJets[j].Et();}
  }

  if(HT  > cut_){return true;}

  return false;

}


std::ostream& HtTriggerCut::Description( std::ostream &ostrm ) {
  ostrm << "HtTriggerCut " << cut_ << " " << " Jet Scale: " <<JetThreshold_ ;
  return ostrm;
}





