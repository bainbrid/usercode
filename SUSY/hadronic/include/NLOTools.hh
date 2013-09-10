#ifndef NLOTOOLS
#define NLOTOOLS

#include "Types.hh"
#include "Math/VectorUtil.h"
#include "ThrustStuff.hh"
#include "TH3.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <strstream>
#include <iostream>
#include <fstream>

class TH1D;
class TH2D;

namespace NLO {

  Double_t GetNLOCross(Double_t m0,Double_t m12, TString process, TString kFactorFile);
  TString GetProcess(Event::Data& ev);

  double ISRProducer(Event::Data& ev);
}

#endif // hadronic_include_BkgdEstPlottingOps_hh
