#include "TriggerTurnOns.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "KinSuite.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "Types.hh"
#include "mt2_bisect.hh"
#include "AlphaT.hh"
#include "Jet.hh"
#include "Math/VectorUtil.h"
#include "JetData.hh"
#include "TMath.h"

using namespace Operation;

// -----------------------------------------------------------------------------
TriggerTurnOns::TriggerTurnOns( const Utils::ParameterSet& ps ) :
// Misc
dirName_( ps.Get<std::string>("DirName") ),
	nMin_( ps.Get<int>("MinObjects") ),
	nMax_( ps.Get<int>("MaxObjects") ),
// MT2
	Plots_( ps.Get<bool>("Plots") )

	{}

// -----------------------------------------------------------------------------
//
TriggerTurnOns::~TriggerTurnOns() {}

// -----------------------------------------------------------------------------
//
void TriggerTurnOns::Start( Event::Data& ev ) {
	initDir( ev.OutputFile(), dirName_.c_str() );
	BookHistos();
}

// -----------------------------------------------------------------------------
//
void TriggerTurnOns::BookHistos() {
	if ( Plots_ )           { Plots(); }
}

// -----------------------------------------------------------------------------
//
bool TriggerTurnOns::Process( Event::Data& ev ) {
	if ( Plots_ )               { Plots(ev); }
	return true;
}

// -----------------------------------------------------------------------------
//
std::ostream& TriggerTurnOns::Description( std::ostream& ostrm ) {
	ostrm << "Hadronic Common Plots ";
	ostrm << "(bins " << nMin_ << " to " << nMax_ << ") ";
	return ostrm;
}

// -----------------------------------------------------------------------------
//
void TriggerTurnOns::Plots() {
	BookHistArray( HT_,
		"HT",
		";n",
		3000,0.,3000.,
		nMax_+1, 0, 1, true );

	BookHistArray( AlphaT_,
		"AlphaT",
		";n",
		500,0.,5.,
		nMax_+1, 0, 1, true );


	BookHistArray( MHT_,
		"MHT",
		";n",
		100,0.,1000.,
		nMax_+1, 0, 1, true );

}

bool TriggerTurnOns::Plots( Event::Data& ev ) {
	UInt_t n = ev.CommonObjects().size();
	Double_t weight = ev.GetEventWeight();

	if ( n >= nMin_ && n <= nMax_ && n < HT_.size()) {
		HT_[0]->Fill(ev.CommonHT(),weight);
		HT_[n]->Fill(ev.CommonHT(),weight);
	}
	if ( n >= nMin_ && n <= nMax_ && n < AlphaT_.size()) {
		AlphaT_[0]->Fill(ev.CommonAlphaT(),weight);
		AlphaT_[n]->Fill(ev.CommonAlphaT(),weight);
	}
	if ( n >= nMin_ && n <= nMax_ && n < MHT_.size()) {
		MHT_[0]->Fill(ev.CommonMHT().Pt(),weight);
		MHT_[n]->Fill(ev.CommonMHT().Pt(),weight);
	}

}
