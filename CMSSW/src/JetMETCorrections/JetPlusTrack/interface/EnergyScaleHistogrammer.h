#ifndef JetMETCorrections_JetPlusTrack_EnergyScaleHistogrammer_H
#define JetMETCorrections_JetPlusTrack_EnergyScaleHistogrammer_H

#include <vector>

namespace edm { 
  class ActivityRegistry; 
  class ParameterSet; 
}
class EnergyScaleHistograms;
class ObjectTags;
class TStyle;
class TFile;

/**
   @brief Service 
*/
class EnergyScaleHistogrammer {
  
 public:
  
  /// Constructor called by framework
  explicit EnergyScaleHistogrammer( const edm::ParameterSet&,
				    const edm::ActivityRegistry& );
  
  /// Destructor called by framework
  virtual ~EnergyScaleHistogrammer();
  
  /// Static method to return instance of this class
  static EnergyScaleHistogrammer* Instance();
  
  /// Static method to return instance of histogram class
  static EnergyScaleHistograms* Histograms( const ObjectTags& );
  
  /// Returns histograms for object defined by tags
  EnergyScaleHistograms* histograms( const ObjectTags& );
  
 private:
  
  /// Private default constructor
  EnergyScaleHistogrammer();
  
  /// Defines root style
  void style();

  /// Tidy up
  void end();
  
 private:
  
  std::vector<EnergyScaleHistograms*> histos_;
  
  TFile* file_;
  
  TStyle* style_;
  
};

#endif // JetMETCorrections_JetPlusTrack_EnergyScaleHistogrammer_H
