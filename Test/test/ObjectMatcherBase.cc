#include "bainbrid/Test/test/ObjectMatcherBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "bainbrid/Test/test/EnergyScaleHistogrammer.h"
#include "bainbrid/Test/test/EnergyScaleHistograms.h"
#include "bainbrid/Test/test/ObjectTags.h"
#include "bainbrid/Test/test/ObjectMatcher.cc"
#include <sstream>
#include <string>

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::ObjectMatcherBase( const edm::ParameterSet& pset ) 
  : tags_( ObjectTags(pset) ),
    verbose_(false)
{
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Constructing...";  
  
  std::stringstream ss;
  ss << "[ObjectMatcherBase::" << __func__ << "] GenObject:" << std::endl
     << " Type      : \"" << tags_.gen().type_ << "\"" << std::endl
     << " Label     : \"" << tags_.gen().label_ << "\"" << std::endl
     << " Instance  : \"" << tags_.gen().instance_ << "\"" << std::endl
     << " Process   : \"" << tags_.gen().process_ << "\"" << std::endl;
  ss << "[ObjectMatcherBase::" << __func__ << "] RecoObject" << std::endl
     << " Type     : \"" << tags_.reco().type_ << "\"" << std::endl
     << " Label    : \"" << tags_.reco().label_ << "\"" << std::endl
     << " Instance : \"" << tags_.reco().instance_ << "\"" << std::endl
     << " Process  : \"" << tags_.reco().process_ << "\"" << std::endl;
  edm::LogVerbatim("EnergyScale") << ss.str();
  
  std::string name = tags_.reco().str() + ".root";
  hOutputFile = new TFile( name.c_str(), "RECREATE" ) ;
  hOutputFile->cd();
  t1 = new TTree("t1","analysis tree");
  t1->Branch("EtaGen1",&EtaGen1,"EtaGen1/D");
  t1->Branch("PhiGen1",&PhiGen1,"PhiGen1/D");
  t1->Branch("EtaRaw1",&EtaRaw1,"EtaRaw1/D");
  t1->Branch("PhiRaw1",&PhiRaw1,"PhiRaw1/D");
  t1->Branch("EtGen1",&EtGen1,"EtGen1/D");
  t1->Branch("EtRaw1",&EtRaw1,"EtRaw1/D");
  t1->Branch("EtMCJ1",&EtMCJ1,"EtMCJ1/D");
  t1->Branch("EtZSP1",&EtZSP1,"EtZSP1/D");
  t1->Branch("EtJPT1",&EtJPT1,"EtJPT1/D");
  t1->Branch("DRMAXgjet1",&DRMAXgjet1,"DRMAXgjet1/D");

  t1->Branch("EtaGen2",&EtaGen2,"EtaGen2/D");
  t1->Branch("PhiGen2",&PhiGen2,"PhiGen2/D");
  t1->Branch("EtaRaw2",&EtaRaw2,"EtaRaw2/D");
  t1->Branch("PhiRaw2",&PhiRaw2,"PhiRaw2/D");
  t1->Branch("EtGen2",&EtGen2,"EtGen2/D");
  t1->Branch("EtRaw2",&EtRaw2,"EtRaw2/D");
  t1->Branch("EtMCJ2",&EtMCJ2,"EtMCJ2/D");
  t1->Branch("EtZSP2",&EtZSP2,"EtZSP2/D");
  t1->Branch("EtJPT2",&EtJPT2,"EtJPT2/D");
  t1->Branch("DRMAXgjet2",&DRMAXgjet2,"DRMAXgjet2/D");

}

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::ObjectMatcherBase() 
  : tags_(),
    verbose_(false)
{;}

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase::~ObjectMatcherBase() {
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Destructing...";  

  hOutputFile->cd() ;
  if ( t1 ) t1->Write();
  if ( hOutputFile ) {
    hOutputFile->Write() ;
    hOutputFile->Close() ;
    delete hOutputFile;
  }
}

// -----------------------------------------------------------------------------
// 
ObjectMatcherBase* ObjectMatcherBase::Instance( const edm::ParameterSet& pset ) {
  LogTrace("EnergyScale")
    << "[ObjectMatcherBase::"<<__func__<<"]"
    << " Creating instance...";  
  
  std::string g = pset.getParameter<std::string>( "GenObjectType" );
  std::string r = pset.getParameter<std::string>( "RecoObjectType" );
  
  if ( g == "GenJet" && 
       r == "CaloJet" ) { 
    return new ObjectMatcher<reco::GenJet,reco::CaloJet>(pset);
  } else if ( g == "GenJet" && 
	      r == "PatJet" ) { 
    return new ObjectMatcher<reco::GenJet,pat::Jet>(pset);
  } else {
    edm::LogError("EnergyScale")
      << "[ObjectMatcherBase::"<<__func__<<"]"
      << " Unexpected string value for ObjectType! "
      << " GenObjectType: \"" << g << "\""
      << " RecoObjectType: \"" << r << "\"";
    return 0;
  } 
  
}

// -----------------------------------------------------------------------------
//
void ObjectMatcherBase::analyze( const edm::Event& event, 
				 const edm::EventSetup& setup ) {
  
  // Build vector of suitable gen objects:
  std::vector<math::XYZTLorentzVector> gen_objects;
  gen( event, setup, gen_objects );
  
  // Check if gen objects were found
  if ( gen_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable GenObjects found to generate 4-vectors!";
    return; 
  } else {
    if ( verbose_ ) {
      LogTrace("EnergyScale")
	<< "[ObjectMatcherBase::" << __func__ << "]"
	<< " Found " << gen_objects.size() << " GenObjects!";
    }
  }

//   {
//     std::vector<math::XYZTLorentzVector>::const_iterator ig = gen_objects.begin();
//     std::vector<math::XYZTLorentzVector>::const_iterator jg = gen_objects.end();
//     for ( ; ig != jg; ++ig ) {
//       std::cout << "TWO GenJets:"
// 		<< " Event: " << event.id().event()
// 		<< " GenJet# " << int( ig - gen_objects.begin() )
// 		<< " e= " << ig->e()
// 		<< " et= " << ig->et()
// 		<< " pt= " << ig->pt()
// 		<< " px= " << ig->px()
// 		<< " py= " << ig->py()
// 		<< " pz= " << ig->pz()
// 		<< std::endl;
//     }
//   }
  
  // Build vector of suitable reco objects:
  std::vector<math::XYZTLorentzVector> reco_objects;
  reco( event, setup, reco_objects );
  
  // Check if reco objects were found
  if ( reco_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable RecoObjects found to generate 4-vectors!";
    return; 
  } else {
    if ( verbose_ ) {
      LogTrace("EnergyScale")
	<< "[ObjectMatcherBase::" << __func__ << "]"
	<< " Found " << reco_objects.size() << " RecoObjects!";
    }
  }

  // Create container for 4-vector pairs
  std::vector<LorentzVectorPair> pairs;
  pairs.clear();
  
  // Iterate through 4-vectors of gen objects
  std::vector<math::XYZTLorentzVector>::const_iterator gg = gen_objects.begin();
  std::vector<math::XYZTLorentzVector>::const_iterator ig = gen_objects.begin();
  std::vector<math::XYZTLorentzVector>::const_iterator jg = gen_objects.end();
  for ( ; ig != jg; ++ig ) {
    
    // Store 4-vector of gen object 
    pairs.push_back( LorentzVectorPair(*ig) );
    
    // Iterate through 4-vectors of reco objects and find match
    std::vector<math::XYZTLorentzVector>::const_iterator rr = reco_objects.end(); 
    std::vector<math::XYZTLorentzVector>::const_iterator ir = reco_objects.begin(); 
    std::vector<math::XYZTLorentzVector>::const_iterator jr = reco_objects.end(); 
    for ( ; ir != jr; ++ir ) {
      if ( !pairs.back().both() || 
	   reco::deltaR( ig->eta(), ig->phi(), ir->eta(), ir->phi() ) < pairs.back().dR() ) {
	pairs.back().reco( *ir ); 
      }
    }
  }
  
//   // Iterate through 4-vectors of gen objects
//   std::vector<math::XYZTLorentzVector>::const_iterator gg = gen_objects.begin();
//   std::vector<math::XYZTLorentzVector>::const_iterator ig = gen_objects.begin();
//   std::vector<math::XYZTLorentzVector>::const_iterator jg = gen_objects.end();
//   for ( ; ig != jg; ++ig ) {

//     // Store 4-vector of gen object 
//     pairs.push_back( LorentzVectorPair(*ig) );
    
//     // Iterate through 4-vectors of reco objects and find match
//     std::vector<math::XYZTLorentzVector>::const_iterator rr = reco_objects.end(); 
//     std::vector<math::XYZTLorentzVector>::const_iterator ir = reco_objects.begin(); 
//     std::vector<math::XYZTLorentzVector>::const_iterator jr = reco_objects.end(); 
//     for ( ; ir != jr; ++ir ) {
//       float dr = pairs.back().both() ? pairs.back().dR() : 1000.;
//       if ( ig->deltaR(*ir) < dr ) { pairs.back().reco( *ir ); }
//     }
//   }

  // Debug
  if ( edm::isDebugEnabled() && verbose_ ) {
    std::vector<LorentzVectorPair>::const_iterator pp = pairs.begin();
    std::vector<LorentzVectorPair>::const_iterator ip = pairs.begin();
    std::vector<LorentzVectorPair>::const_iterator jp = pairs.end();
    for ( ; ip != jp; ++ip ) {
      std::stringstream ss;
      ss << "[ObjectMatcherBase::" << __func__ << "]"
	 << " Matched pair #" << static_cast<uint32_t>( ip - pp );
      if ( ip->both() ) {
	ss << " D(R): " << ip->dR()
	   << " D(eta): " << ip->dEta()
	   << " D(phi): " << ip->dPhi();
      }
      ss << std::endl
	 << "  GenObject:  " 
	 << " e: " << ip->gen().e() 
	 << " pt: " << ip->gen().pt() 
	 << " px: " << ip->gen().px() 
	 << " py: " << ip->gen().py() 
	 << " pz: " << ip->gen().pz() 
	 << " eta: " << ip->gen().eta() 
	 << " phi: " << ip->gen().phi() 
	 << std::endl;
      if ( ip->both() ) {
	ss << "  RecoObject: " 
	   << " e: " << ip->reco().e() 
	   << " pt: " << ip->reco().pt() 
	   << " px: " << ip->reco().px() 
	   << " py: " << ip->reco().py() 
	   << " pz: " << ip->reco().pz() 
	   << " eta: " << ip->reco().eta() 
	   << " phi: " << ip->reco().phi();
      } else {
	ss << "  RecoObject: no match!";
      }
      LogTrace("EnergyScale") << ss.str();
    }
  }
  
  // "Old skool"
  if (0) { analyze( event, gen_objects, reco_objects ); } //@@ Sasha's
  else { analyze( event, pairs ); } //@@ Hybrid
  
  // Access to histos object
  EnergyScaleHistograms* histos = EnergyScaleHistogrammer::Histograms( tags_ ); 
  if ( histos ) { histos->analyze( pairs ); }
  else {
    edm::LogError("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " NULL pointer to EnergyScaleHistograms object!";
    return;
  }
  
}

// -----------------------------------------------------------------------------
//
std::string ObjectMatcherBase::id() {
  return std::string( "<" + tags_.gen().type_ + "," + tags_.reco().type_ + ">" );
}

// -----------------------------------------------------------------------------
//
void ObjectMatcherBase::analyze( const edm::Event& iEvent, 
				 const std::vector<math::XYZTLorentzVector>& gen_objects,
				 const std::vector<math::XYZTLorentzVector>& reco_objects ) {

  // Check if gen objects were found
  if ( gen_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable GenObjects found to generate 4-vectors!";
    return; 
  } else {
    if ( verbose_ ) {
      LogTrace("EnergyScale")
	<< "[ObjectMatcherBase::" << __func__ << "]"
	<< " Found " << gen_objects.size() << " GenObjects!";
    }
  }
  
  // Check if reco objects were found
  if ( reco_objects.empty() ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " No suitable RecoObjects found to generate 4-vectors!";
    return; 
  } else {
    if ( verbose_ ) {
      LogTrace("EnergyScale")
	<< "[ObjectMatcherBase::" << __func__ << "]"
	<< " Found " << reco_objects.size() << " RecoObjects!";
    }
  }

  std::vector<math::XYZTLorentzVector> gjets;
  gjets.clear();

  // initialize tree variables
  EtaGen1 = 0.;
  PhiGen1 = 0.;
  EtaRaw1 = 0.;
  PhiRaw1 = 0.;
  EtGen1  = 0.;
  EtRaw1  = 0.;
  EtMCJ1  = 0.;
  EtZSP1  = 0.;
  EtJPT1  = 0.;
  DRMAXgjet1 = 1000.;

  EtaGen2 = 0.;
  PhiGen2 = 0.;
  EtaRaw2 = 0.;
  PhiRaw2 = 0.;
  EtGen2  = 0.;
  EtRaw2  = 0.;
  EtMCJ2  = 0.;
  EtZSP2  = 0.;
  EtJPT2  = 0.;
  DRMAXgjet2 = 1000.;

  // l1 and l2 are leptons from Z->ll to be checked they are not gen jets (DR match)
  math::XYZTLorentzVector l1(0.,0.,1.,1.);
  math::XYZTLorentzVector l2(0.,0.,1.,1.);

  int jg = 0;
  //    for(GenJetCollection::const_iterator gjet = genjets->begin(); 
  //        gjet != genjets->end(); ++gjet ) {
  for(std::vector<math::XYZTLorentzVector>::const_iterator gjet = gen_objects.begin(); 
      gjet != gen_objects.end(); ++gjet ) {

      math::XYZTLorentzVector jet(gjet->px(), gjet->py(), gjet->pz(), gjet->e());

//       std::cout << "ALL GenJets:"
// 		<< " Event: " << iEvent.id().event()
// 		<< " GenJet# " << int( gjet - gen_objects.begin() )
// 		<< " e= " << jet.e()
// 		<< " et= " << jet.et()
// 		<< " pt= " << jet.pt()
// 		<< " px= " << jet.px()
// 		<< " py= " << jet.py()
// 		<< " pz= " << jet.pz()
// 		<< std::endl;

    if(gjet->pt() >= 20.) {
      double drjl1 = reco::deltaR( l1.eta(), l1.phi(), jet.eta(), jet.phi() );
      double drjl2 = reco::deltaR( l2.eta(), l2.phi(), jet.eta(), jet.phi() );

//       std::cout <<" Gen Jet " << jg
// 		<<" pt = " << gjet->pt()
// 		<<" px = " << gjet->px()
// 		<<" py = " << gjet->py()
// 		<<" pz = " << gjet->pz()
// 		<<" energy = " << gjet->energy()
// 		<<" j eta = " << gjet->eta()
// 		<<" j phi = " << gjet->phi() 
// 		<<" l1 eta = " << l1.eta() 
// 		<<" l1 phi = " << l1.phi() 
// 		<<" l2 eta = " << l2.eta() 
// 		<<" l2 phi = " << l2.phi() 
// 		<<" dr1 = " << drjl1 
// 		<<" dr2 = " << drjl2 << std::endl;

      if(drjl1 > 1.0 && drjl2 > 1.0) 
	{
	  jg++;
	  if(jg <= 2) {
	    gjets.push_back(jet);
	  }
	}

    }
  }

//   std::vector<math::XYZTLorentzVector>::const_iterator ii = gjets.begin();
//   std::vector<math::XYZTLorentzVector>::const_iterator jj = gjets.end();
//   for ( ; ii != jj; ++ii ) {
//     std::cout << "TWO GenJets:"
// 	      << " Event: " << iEvent.id().event()
// 	      << " GenJet# " << int( ii - gjets.begin() )
// 	      << " e= " << ii->e()
// 	      << " et= " << ii->et()
// 	      << " pt= " << ii->pt()
// 	      << " px= " << ii->px()
// 	      << " py= " << ii->py()
// 	      << " pz= " << ii->pz()
// 	      << std::endl;
//   }
  
  if(gjets.size() > 0) {
       
    // loop over jets and do matching with gen jets
    int jc = 0;
       
    //       for( CaloJetCollection::const_iterator cjet = calojets->begin(); 
    // 	   cjet != calojets->end(); ++cjet ){ 
    for(std::vector<math::XYZTLorentzVector>::const_iterator cjet = reco_objects.begin(); 
	cjet != reco_objects.end(); ++cjet ){ 
      //
      math::XYZTLorentzVector cjetv(cjet->px(), cjet->py(), cjet->pz(), cjet->e());

//       std::cout << "ALL RecoJets:"
// 		<< " Event: " << iEvent.id().event()
// 		<< " RecoJet# " << int( cjet - reco_objects.begin() ) 
// 		<< " e= " << cjetv.e()
// 		<< " et= " << cjetv.et()
// 		<< " pt= " << cjetv.pt()
// 		<< " px= " << cjetv.px()
// 		<< " py= " << cjetv.py()
// 		<< " pz= " << cjetv.pz()
// 		<< std::endl;
    
      double DRgjet1 = reco::deltaR( gjets[0].eta(), gjets[0].phi(), cjetv.eta(), cjetv.phi() );
      
      if(DRgjet1 < DRMAXgjet1) {
	DRMAXgjet1 = DRgjet1;
	
	EtaGen1 = gjets[0].eta();
	PhiGen1 = gjets[0].phi();
	EtGen1  = gjets[0].pt();
	
	EtaRaw1 = cjet->eta(); 
	PhiRaw1 = cjet->phi();
	EtRaw1  = cjet->pt();
	EtJPT1  = cjetv.pt(); 
      }
      if(gjets.size() == 2) {
	double DRgjet2 = reco::deltaR( gjets[1].eta(), gjets[1].phi(), cjetv.eta(), cjetv.phi() );
	if(DRgjet2 < DRMAXgjet2) { 
	  DRMAXgjet2 = DRgjet2;
	  
	  EtaGen2 = gjets[1].eta();
	  PhiGen2 = gjets[1].phi();
	  EtGen2  = gjets[1].pt();

	  EtaRaw2 = cjet->eta(); 
	  PhiRaw2 = cjet->phi();
	  EtRaw2  = cjet->pt();
	  EtJPT2  = cjetv.pt(); 
	}
      }
      jc++;
    }
  }
  // fill tree
  t1->Fill();

//   std::cout << " Event: " << iEvent.id().event() << std::endl
// 	    << " Scale1: " << ( EtGen1 > 0. ? EtJPT1/EtGen1 : -1. )
// 	    << "  EtJPT1: " << EtJPT1
// 	    << " EtGen1: " << EtGen1
// 	    << " DRMAXgjet1: " << DRMAXgjet1
// 	    << " EtaGen1: " << EtaGen1
// 	    << " PhiGen1: " << PhiGen1
// 	    << std::endl
// 	    << " Scale2: " << ( EtGen2 > 0. ? EtJPT2/EtGen2 : -1. )
// 	    << "  EtJPT2: " << EtJPT2
// 	    << " EtGen2: " << EtGen2
// 	    << " DRMAXgjet2: " << DRMAXgjet2
// 	    << " EtaGen2: " << EtaGen2 
// 	    << " PhiGen2: " << PhiGen2 
// 	    << std::endl;
  
}

// -----------------------------------------------------------------------------
//
void ObjectMatcherBase::analyze( const edm::Event& iEvent, 
				 const std::vector<LorentzVectorPair>& pairs ) {
  
  // Check if gen objects were found
  if ( pairs.size() != 2 ) { 
    edm::LogWarning("EnergyScale")
      << "[ObjectMatcherBase::" << __func__ << "]"
      << " Size of LorentzVectorPair vector is not equal to two!";
  } else {
    if ( verbose_ ) {
      LogTrace("EnergyScale")
	<< "[ObjectMatcherBase::" << __func__ << "]"
	<< " Found " << pairs.size() << " LorentzVectorPair objects!";
    }
  }

  // initialize tree variables
  EtaGen1 = 0.;
  PhiGen1 = 0.;
  EtaRaw1 = 0.;
  PhiRaw1 = 0.;
  EtGen1  = 0.;
  EtRaw1  = 0.;
  EtMCJ1  = 0.;
  EtZSP1  = 0.;
  EtJPT1  = 0.;
  DRMAXgjet1 = 1000.;

  EtaGen2 = 0.;
  PhiGen2 = 0.;
  EtaRaw2 = 0.;
  PhiRaw2 = 0.;
  EtGen2  = 0.;
  EtRaw2  = 0.;
  EtMCJ2  = 0.;
  EtZSP2  = 0.;
  EtJPT2  = 0.;
  DRMAXgjet2 = 1000.;
  
  if ( pairs.size() > 0 && pairs[0].both() ) {
    DRMAXgjet1 = pairs[0].dR();
    EtaGen1 = pairs[0].gen().eta();
    PhiGen1 = pairs[0].gen().phi();
    EtGen1  = pairs[0].gen().pt();
    //EtRaw1  = pairs[0].reco().pt(); //@@
    EtJPT1  = pairs[0].reco().pt(); 
  }
  
  if ( pairs.size() > 1 && pairs[1].both() ) {
    DRMAXgjet2 = pairs[1].dR();
    EtaGen2 = pairs[1].gen().eta();
    PhiGen2 = pairs[1].gen().phi();
    EtGen2  = pairs[1].gen().pt();
    //EtRaw2  = pairs[1].reco().pt(); //@@
    EtJPT2  = pairs[1].reco().pt(); 
  }
  
  t1->Fill();

//   std::cout << " Event: " << iEvent.id().event() << std::endl
// 	    << " Scale1: " << ( EtGen1 > 0. ? EtJPT1/EtGen1 : -1. )
// 	    << "  EtJPT1: " << EtJPT1
// 	    << " EtGen1: " << EtGen1
// 	    << " DRMAXgjet1: " << DRMAXgjet1
// 	    << " EtaGen1: " << EtaGen1
// 	    << " PhiGen1: " << PhiGen1
// 	    << std::endl
// 	    << " Scale2: " << ( EtGen2 > 0. ? EtJPT2/EtGen2 : -1. )
// 	    << "  EtJPT2: " << EtJPT2
// 	    << " EtGen2: " << EtGen2
// 	    << " DRMAXgjet2: " << DRMAXgjet2
// 	    << " EtaGen2: " << EtaGen2 
// 	    << " PhiGen2: " << PhiGen2 
// 	    << std::endl;

}





//   // Create 4-vector pairs
//   std::vector<LorentzVectorPair> pairs;
//   pairs.resize( gen_objects.size() );
  
//   // Iterate through 4-vectors of reco objects and find match
//   std::vector<math::XYZTLorentzVector>::const_iterator rr = reco_objects.end(); 
//   std::vector<math::XYZTLorentzVector>::const_iterator ir = reco_objects.begin(); 
//   std::vector<math::XYZTLorentzVector>::const_iterator jr = reco_objects.end(); 
//   for ( ; ir != jr; ++ir ) {
//     // Iterate through 4-vectors of gen objects
//     std::vector<math::XYZTLorentzVector>::const_iterator gg = gen_objects.begin();
//     std::vector<math::XYZTLorentzVector>::const_iterator ig = gen_objects.begin();
//     std::vector<math::XYZTLorentzVector>::const_iterator jg = gen_objects.end();
//     for ( ; ig != jg; ++ig ) {
//       LorentzVectorPair& tmp = pairs[ static_cast<uint32_t>( ig - gg ) ];
//       tmp.gen( *ig );
//       if ( !tmp.both() || ig->deltaR(*ir) < tmp.dR() ) {
// 	tmp.reco( *ir ); 
//       }
//     }
//   }





//   // Create 4-vector pairs
//   std::vector<LorentzVectorPair> pairs;

//   // Record of matched reco objects
//   std::vector< std::vector<math::XYZTLorentzVector>::const_iterator > matched;
  
//   // Iterate through 4-vectors of gen objects
//   std::vector<math::XYZTLorentzVector>::const_iterator ig = gen_objects.begin();
//   std::vector<math::XYZTLorentzVector>::const_iterator jg = gen_objects.end();
//   for ( ; ig != jg; ++ig ) {
    
//     // Create 4-vector pair
//     LorentzVectorPair match( *ig );
    
//     // Iterate through 4-vectors of reco objects and find match
//     double dr = 1.e6;
//     std::vector<math::XYZTLorentzVector>::const_iterator iter = reco_objects.end(); 
//     std::vector<math::XYZTLorentzVector>::const_iterator ir = reco_objects.begin(); 
//     std::vector<math::XYZTLorentzVector>::const_iterator jr = reco_objects.end(); 
//     for ( ; ir != jr; ++ir ) {
//       //       if ( std::find( matched.begin(), 
//       // 		      matched.end(), 
//       // 		      ir ) == matched.end() ) { // Check not already matched
// 	double idr = ig->deltaR(*ir);
// 	if ( idr < dr ) { // Check if closer in dR
// 	  iter = ir;
// 	  dr = idr;
// 	}
// 	//       }
//     }

//     // Check if match is found 
//     if ( iter != jr ) { 
      
//       matched.push_back(iter);
//       match.reco(*iter);
      
//       // Debug
//       if ( edm::isDebugEnabled() ) {
// 	std::stringstream ss;
// 	ss << "[ObjectMatcherBase::" << __func__ << "]"
// 	   << " Matched GenObject #" << static_cast<uint32_t>( ig - gen_objects.begin() )
// 	   << " to RecoObject #" << static_cast<uint32_t>( iter - reco_objects.begin() ) << std::endl
// 	   << "  GenObject:  pt/eta/phi: " 
// 	   << ig->pt() << "/" 
// 	   << ig->eta() << "/"
// 	   << ig->phi() << std::endl
// 	   << "  RecoObject: pt/eta/phi: " 
// 	   << iter->pt() << "/" 
// 	   << iter->eta() << "/"
// 	   << iter->phi() << std::endl
// 	   << "  DR: " << dr << std::endl;
// 	LogTrace("EnergyScale") << ss.str();
//       }

//     }
    
//     pairs.push_back(match);
    
//   }
