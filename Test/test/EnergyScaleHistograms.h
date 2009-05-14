#ifndef JetMETCorrections_JetPlusTrack_EnergyScaleHistograms_H
#define JetMETCorrections_JetPlusTrack_EnergyScaleHistograms_H

#include "bainbrid/Test/test/LorentzVectorPair.h"
#include "bainbrid/Test/test/ObjectTags.h" 
#include <vector>

class EnergyScaleHistogrammer;
class TFile;
class TH1F;
class TStyle;

/**
   @brief Service 
*/
class EnergyScaleHistograms {
  
 public:
  
  /// Compares 4-vectors from generator and reco objects
  void analyze( const std::vector<LorentzVectorPair>& );
  
  /// Tags used to identify gen and reco objects
  ObjectTags tags() const;

 private:

  // Only service has access to private constructors
  friend class EnergyScaleHistogrammer;

  /// Constructor called by service
  explicit EnergyScaleHistograms( const ObjectTags& );
  
  /// Private default constructor
  EnergyScaleHistograms();

  /// Destructor
  virtual ~EnergyScaleHistograms();

  /// Books histogram objects
  void book( TFile* const );
  
  /// Writes histogram objects to root file
  void write( TFile* const );

  /// Creates scale and resolution histograms
  void end( TStyle* const );

  void cd( TFile* const );

 private:
  
  /// Tags used to identify objects
  ObjectTags tags_;

  /// Number of energy bins
  static const uint32_t nBins_ = 10;
  
  /// Energy bins
  static const float eBins_[nBins_+1];

  /// Scale for each energy bin
  std::vector<TH1F*> hScale_; 
  
  /// Transverse energy of generator object
  TH1F* hEt_;
  
  /// DR b/w generator and reco objects
  TH1F* hDR_;
  
  /// Scale versus energy bin
  TH1F* hScaleVsE_;
  
  /// Resolution versus energy bin
  TH1F* hResVsE_;
  
};

inline ObjectTags EnergyScaleHistograms::tags() const { return tags_; }

#endif // JetMETCorrections_JetPlusTrack_EnergyScaleHistograms_H
