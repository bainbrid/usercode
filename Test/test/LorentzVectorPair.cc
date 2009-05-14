#include "bainbrid/Test/test/LorentzVectorPair.h"

// -----------------------------------------------------------------------------
// 
const double LorentzVectorPair::pi_ = 3.14159265359;

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair( const HepLorentzVector& gen_object,
				      const HepLorentzVector& reco_object ) 
  : gen_( gen_object ),
    reco_( reco_object ),
    g_(true),
    r_(true)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair( const HepLorentzVector& gen_object ) 
  : gen_( gen_object ),
    reco_(),
    g_(true),
    r_(false)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair() 
  : gen_(),
    reco_(),
    g_(false),
    r_(false)
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::LorentzVectorPair( const LorentzVectorPair& input ) 
  : gen_( input.gen_ ),
    reco_( input.reco_ ),
    g_( input.g_ ),
    r_( input.r_ )
{;}

// -----------------------------------------------------------------------------
// 
LorentzVectorPair::~LorentzVectorPair() {;}
