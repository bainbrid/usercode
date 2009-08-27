#include "TestNtuplizer.h"
#include "CommonTools/Utils/interface/EtComparator.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

TestNtuplizer::TestNtuplizer( const edm::ParameterSet& pset ) 
  : tree_(0),
    n_(0)
{
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Data","data");
  tree_->SetAutoSave(10);
  tree_->Branch( "Test", &n_, "Test/int" );
  tree_->Branch( "Array", e_, "Array[Test]/double" );
}

// -----------------------------------------------------------------------------
//
void TestNtuplizer::analyze( const edm::Event& event, 
			 const edm::EventSetup& setup ) {
  n_ = event.id().event();
  e_[n_] = event.id().event();
  tree_->Fill();
}

// -----------------------------------------------------------------------------
//
void TestNtuplizer::beginJob( const edm::EventSetup& ) {}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestNtuplizer);
