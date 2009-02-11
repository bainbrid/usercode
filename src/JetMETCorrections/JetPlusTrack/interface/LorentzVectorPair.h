#ifndef JetMETCorrections_JetPlusTrack_LorentzVectorPair_H
#define JetMETCorrections_JetPlusTrack_LorentzVectorPair_H

#include "CLHEP/Vector/LorentzVector.h"

class LorentzVectorPair {

 public:

  // ---------- Con(de)structors ----------
  
  /// 4-vectors of gen and reco object
  explicit LorentzVectorPair( const HepLorentzVector& gen_object,
			      const HepLorentzVector& reco_object );
  
  /// 4-vector of just gen object
  explicit LorentzVectorPair( const HepLorentzVector& gen_object );
  
  /// Default constructor
  LorentzVectorPair();

  /// Copy constructor
  LorentzVectorPair( const LorentzVectorPair& );

  /// Destructor
  ~LorentzVectorPair();

  // ---------- Setter methods ----------
  
  /// Add 4-vector for gen object 
  void gen( const HepLorentzVector& gen_object );
  
  /// Add gen object 
  void gen( const double& px,
	    const double& py,
	    const double& pz,
	    const double& e );
  
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
  
  /// DR between gen and reco objects
  double deltaR() const;
  
  // ---------- Utility methods ----------
  
  static const double pi_;
  static double dEta( double, double );
  static double dPhi( double, double );
  static double dR( double, double );
  
 private:
  
  /// 4-vector of gen object
  HepLorentzVector gen_;
  
  /// 4-vector of reco object
  HepLorentzVector reco_;
  
  /// Indicates if gen object is defined
  bool g_;
  
  /// Indicates if reco object is defined
  bool r_;
  
};

// ---------- Inline methods ----------

inline void LorentzVectorPair::gen( const HepLorentzVector& object ) { 
  gen_ = object; 
  g_ = true;
}

inline void LorentzVectorPair::gen( const double& px,
				    const double& py,
				    const double& pz,
				    const double& e ) { 
  gen_ = HepLorentzVector(px,py,pz,e); 
  g_ = true;
}

inline void LorentzVectorPair::reco( const HepLorentzVector& object ) { 
  reco_ = object; 
  r_ = true;
}

inline void LorentzVectorPair::reco( const double& px,
				     const double& py,
				     const double& pz,
				     const double& e ) { 
  reco_ = HepLorentzVector(px,py,pz,e); 
  r_ = true;
}

inline HepLorentzVector LorentzVectorPair::gen() const { 
  if ( g_ ) { return gen_; }
  else { return HepLorentzVector(); }
}

inline HepLorentzVector LorentzVectorPair::reco() const {
  if ( r_ ) { return reco_; }
  else { return HepLorentzVector(); }
}

inline bool LorentzVectorPair::both() const { return g_ && r_; }

inline double LorentzVectorPair::dEta() const { 
  if ( both() ) { return dEta( reco_.eta(), gen_.eta() ); }
  else { return -1.; }
}

inline double LorentzVectorPair::dPhi() const { 
  if ( both() ) { return dPhi( reco_.phi(), gen_.phi() ); }
  else { return -1.; }
}

inline double LorentzVectorPair::deltaR() const { 
  if ( both() ) { return gen_.deltaR(reco_); }
  else { return -1.; }
}

inline double LorentzVectorPair::dR() const { 
  if ( both() ) { return dR( dEta(), dPhi() ); }
  else { return -1.; }
}

inline double LorentzVectorPair::dEta( double eta1, double eta2 ) {
  return fabs( eta1 - eta2 ); 
}

inline double LorentzVectorPair::dPhi( double phi1, double phi2 ) {
  double dphi = fabs( phi1 - phi2 );
  return ( dphi > pi_ ?  2. * pi_- dphi : dphi );
}

/* inline double LorentzVectorPair::dEta( double phi1, double phi2 ) { */
/*   double dphi = fabs( phi1 - phi2 ); */
/*   return ( dphi > pi_ ?  2. * pi_- dphi : dphi ); */
/* } */

/* inline double LorentzVectorPair::dPhi( double eta1, double eta2 ) { */
/*   return fabs( eta1 - eta2 );  */
/* } */

inline double LorentzVectorPair::dR( double deta, double dphi ) {
  return sqrt( deta * deta + dphi * dphi );
}

#endif // ObjectMETCorrections_ObjectPlusTrack_LorentzVectorPair_H
