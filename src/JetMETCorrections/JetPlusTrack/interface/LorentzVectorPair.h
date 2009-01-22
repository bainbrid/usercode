#ifndef JetMETCorrections_JetPlusTrack_LorentzVectorPair_H
#define JetMETCorrections_JetPlusTrack_LorentzVectorPair_H

#include "CLHEP/Vector/LorentzVector.h"

class LorentzVectorPair {

 public:

  // ---------- Con(de)structors ----------
  
  /// 4-vector of just gen object
  explicit LorentzVectorPair( const HepLorentzVector& gen_object );

  /// 4-vectors of gen and reco object
  explicit LorentzVectorPair( const HepLorentzVector& gen_object,
			      const HepLorentzVector& reco_object );
  
  /// Destructor
  ~LorentzVectorPair();

  // ---------- Setter methods ----------
  
  /// Add 4-vector for reco object 
  void reco( const HepLorentzVector& reco_object );
  
  /// Add reco object 
  void reco( const double& px,
	     const double& py,
	     const double& pz,
	     const double& e );

  // ---------- Getter methods ----------
  
  /// 4-vector of gen object
  HepLorentzVector gen() const;

  /// 4-vector of reco object
  HepLorentzVector reco() const;
  
  /// Indicates if both gen and reco objects are defined
  bool both() const;
  
  /// Eta between gen and reco objects
  double dEta() const;
  
  /// Phi between gen and reco objects
  double dPhi() const;
  
  /// DR between gen and reco objects
  double dR() const;

 private:
  
  // ---------- Utility methods ----------
  
  static const double pi_;
  double dEta( double, double ) const;
  double dPhi( double, double ) const;
  double dR( double, double ) const;
  
  /// Private default constructor
  LorentzVectorPair();

 private:
  
  /// 4-vector of gen object
  HepLorentzVector gen_;
  
  /// 4-vector of reco object
  HepLorentzVector reco_;

  /// Indicates if both gen and reco objects are defined
  bool both_;

};

// ---------- Inline methods ----------

inline void LorentzVectorPair::reco( const HepLorentzVector& object ) { 
  reco_ = object; 
  both_ = true;
}

inline void LorentzVectorPair::reco( const double& px,
				 const double& py,
				 const double& pz,
				 const double& e ) { 
  reco_ = HepLorentzVector(px,py,pz,e); 
  both_ = true;
}

inline HepLorentzVector LorentzVectorPair::gen() const { return gen_; }

inline HepLorentzVector LorentzVectorPair::reco() const {
  if ( both_ ) { return reco_; }
  else { return HepLorentzVector(); }
}

inline double LorentzVectorPair::dEta() const { 
  if ( both_ ) { return dEta( reco_.eta(), gen_.eta() ); }
  else { return 0.; }
}

inline bool LorentzVectorPair::both() const { return both_; }

inline double LorentzVectorPair::dPhi() const { 
  if ( both_ ) { return dPhi( reco_.phi(), gen_.phi() ); }
  else { return 0.; }
}

inline double LorentzVectorPair::dR() const { 
  if ( both_ ) { return dR( dEta(), dPhi() ); }
  else { return 0.; }
}

inline double LorentzVectorPair::dEta( double phi1, double phi2 ) const {
  double dphi = fabs( phi1 - phi2 );
  return ( dphi > pi_ ?  2. * pi_- dphi : dphi );
}

inline double LorentzVectorPair::dPhi( double eta1, double eta2 ) const {
  return fabs( eta1 - eta2 ); 
}

inline double LorentzVectorPair::dR( double deta, double dphi ) const {
  return sqrt( deta* deta + dphi * dphi );
}

#endif // ObjectMETCorrections_ObjectPlusTrack_LorentzVectorPair_H
