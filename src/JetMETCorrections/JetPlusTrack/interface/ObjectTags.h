#ifndef JetMETCorrections_JetPlusTrack_ObjectTags_H
#define JetMETCorrections_JetPlusTrack_ObjectTags_H

#include "FWCore/ParameterSet/interface/InputTag.h"
#include <string>

namespace edm { class ParameterSet; }

/**
   Container for sim and reco object types and tags
*/
class ObjectTags {
  
 public:

  ObjectTags( const edm::ParameterSet& );
  
  ObjectTags( const ObjectTags& );

  ObjectTags();

  bool operator== ( const ObjectTags& ) const;
  
  class Tag {
  public:
    Tag( const std::string& type, 
	 const edm::InputTag& tag );
    Tag( const Tag& );
    Tag();
    bool operator== ( const Tag& ) const;
    std::string str() const;
    edm::InputTag tag() const;
  public:
    std::string type_;
    std::string label_;
    std::string instance_;
    std::string process_;
  };
  
  Tag gen() const;
  
  Tag reco() const;
  
 private:
  
  Tag gen_;
  
  Tag reco_;
  
};

inline ObjectTags::Tag ObjectTags::gen() const { return gen_; }
inline ObjectTags::Tag ObjectTags::reco() const { return reco_; }

#endif // JetMETCorrections_JetPlusTrack_ObjectTags_H
