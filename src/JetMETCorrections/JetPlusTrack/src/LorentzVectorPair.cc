#include "JetMETCorrections/JetPlusTrack/interface/LorentzVectorPair.h"

// -----------------------------------------------------------------------------
// 
const double LorentzVectorPair::pi_ = 3.14159265359;

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair() 
  : gen_(),
    reco_(),
    both_(false)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair( const HepLorentzVector& gen_object ) 
  : gen_( gen_object ),
    reco_(),
    both_(false)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair( const HepLorentzVector& gen_object,
			      const HepLorentzVector& reco_object ) 
  : gen_( gen_object ),
    reco_( reco_object ),
    both_(true)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::~LorentzVectorPair() {;}
