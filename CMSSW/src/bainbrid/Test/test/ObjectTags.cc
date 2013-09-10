#include "bainbrid/Test/test/ObjectTags.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// -----------------------------------------------------------------------------
// 
ObjectTags::ObjectTags( const edm::ParameterSet& pset ) 
  : gen_( pset.getParameter<std::string>("GenObjectType"),
	  pset.getParameter<edm::InputTag>( "GenObjectTag" ) ),
    reco_( pset.getParameter<std::string>("RecoObjectType"),
	   pset.getParameter<edm::InputTag>( "RecoObjectTag" ) ) 
{;}

// -----------------------------------------------------------------------------
// 
ObjectTags::ObjectTags( const ObjectTags& input ) 
  : gen_( input.gen() ),
    reco_( input.reco() )
{;}

// -----------------------------------------------------------------------------
// 
ObjectTags::ObjectTags() 
  : gen_(),
    reco_()
{;}

// -----------------------------------------------------------------------------
// 
bool ObjectTags::operator== ( const ObjectTags& input ) const {
  return ( gen_ == input.gen() &&
	   reco_ == input.reco() );
}

// -----------------------------------------------------------------------------
// 
ObjectTags::Tag::Tag( const std::string& type, 
		      const edm::InputTag& tag ) 
  : type_( type ),
    label_( tag.label() ),
    instance_( tag.instance() ),
    process_( tag.process() )
{;}

// -----------------------------------------------------------------------------
// 
ObjectTags::Tag::Tag( const Tag& input ) 
  : type_( input.type_ ),
    label_( input.label_ ),
    instance_( input.instance_ ),
    process_( input.process_ )
{;}

// -----------------------------------------------------------------------------
// 
ObjectTags::Tag::Tag() 
  : type_(""),
    label_(""),
    instance_(""),
    process_("")
{;}

// -----------------------------------------------------------------------------
// 
bool ObjectTags::Tag::operator== ( const Tag& input ) const {
  return ( type_ == input.type_ &&
	   label_ == input.label_ &&
	   instance_ == input.instance_ &&
	   process_ == input.process_ );
}

// -----------------------------------------------------------------------------
// 
std::string ObjectTags::Tag::str() const { 
  std::stringstream ss;
  ss << type_ << "_" 
     << label_ << "_"
     << instance_ << "_"
     << process_;
  return ss.str();
}

// -----------------------------------------------------------------------------
// 
edm::InputTag ObjectTags::Tag::tag() const { 
  std::stringstream ss;
  ss << label_ << ":"
     << instance_ << ":"
     << process_;
  return edm::InputTag( ss.str() );
}
