#ifndef PhysicsTools_PatAlgos_PATSingleVertexSelector_h
#define PhysicsTools_PatAlgos_PATSingleVertexSelector_h
//
// $Id: PATSingleVertexSelector.h,v 1.1 2008/07/22 12:47:02 gpetrucc Exp $
//

/**
  \class    pat::PATSingleVertexSelector PATSingleVertexSelector.h "PhysicsTools/PatAlgos/plugins/PATSingleVertexSelector.h"
  \brief    Produces a list containing a single vertex selected by some criteria


  \author   Giovanni Petrucciani
  \version  $Id: PATSingleVertexSelector.h,v 1.1 2008/07/22 12:47:02 gpetrucc Exp $
*/

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "PhysicsTools/Utilities/interface/StringCutObjectSelector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"

namespace pat {

  class PATSingleVertexSelector : public edm::EDFilter {

    public:

      explicit PATSingleVertexSelector(const edm::ParameterSet & iConfig);
      ~PATSingleVertexSelector();

      virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      enum Mode { First, NearestToCand, FromCand, FromBeamSpot };
      typedef StringCutObjectSelector<reco::Vertex>    VtxSel;
      typedef StringCutObjectSelector<reco::Candidate> CandSel;

      static Mode parseMode(const std::string &name) ;
      bool filter_(Mode mode, edm::Event & iEvent, const edm::EventSetup & iSetup) ;
      bool hasMode_(Mode mode) const ;
      // configurables
      std::vector<Mode> modes_; // mode + optional fallbacks
      edm::InputTag vertices_;
      std::vector<edm::InputTag> candidates_;
      std::auto_ptr<VtxSel > vtxPreselection_;
      std::auto_ptr<CandSel> candPreselection_;
      edm::InputTag beamSpot_;
      // transient data. meaningful while 'filter()' is on the stack
      std::vector<const reco::Vertex *> selVtxs_;
      const reco::Candidate *           bestCand_;
      
  };

}

#endif

