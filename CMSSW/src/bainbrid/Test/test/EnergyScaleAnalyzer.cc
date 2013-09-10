#include "bainbrid/Test/test/EnergyScaleAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "bainbrid/Test/test/EnergyScaleHistogrammer.h" 
#include "bainbrid/Test/test/EnergyScaleHistograms.h" 
#include "bainbrid/Test/test/ObjectMatcherBase.h" 

// -----------------------------------------------------------------------------
// 
EnergyScaleAnalyzer::EnergyScaleAnalyzer( const edm::ParameterSet& pset ) 
  : matcher_( ObjectMatcherBase::Instance(pset) )
{
  LogTrace("EnergyScale")
    << "[EnergyScaleAnalyzer::"<<__func__<<"]"
    << " Constructing...";  
  if ( matcher_ ) { matcher_->verbose( pset.getUntrackedParameter<bool>("verbose",false) ); }
}

// -----------------------------------------------------------------------------
// 
EnergyScaleAnalyzer::~EnergyScaleAnalyzer() {
  LogTrace("EnergyScale")
    << "[EnergyScaleAnalyzer::"<<__func__<<"]"
    << " Destructing...";  

  // Do stuff here
  if ( matcher_ ) { delete matcher_; }
}

// -----------------------------------------------------------------------------
// 
void EnergyScaleAnalyzer::analyze( const edm::Event& event, 
				   const edm::EventSetup& setup ) {
  if ( matcher_ ) { matcher_->analyze( event, setup ); } 
}
