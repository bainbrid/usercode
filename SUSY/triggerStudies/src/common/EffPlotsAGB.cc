#include "EffPlotsAGB.hh"

using namespace AGBTrig;

AGBTrigEffPlots::AGBTrigEffPlots(const Utils::ParameterSet& ps):
  _dirName(ps.Get<std::string>("DirName")),
  _TrigMHTcut(ps.Get<double>("MHTcut")),
  _JetPtThres(ps.Get<double>("JetPt"))
{}


AGBTrigEffPlots::~AGBTrigEffPlots() {}

void AGBTrigEffPlots::Start(Event::Data& ev){
  initDir(ev.OutputFile(), _dirName.c_str());
  BookHistos();
}

void AGBTrigEffPlots::BookHistos() {

  pfCommonMHT = new TH1D("pfCommonMHT","pfCommonMHT",200,0,100);
  pfCalcMHT = new TH1D("pfCalcMHT","pfCalcMHT",200,0,100);
  pfCommonMHT_postCut = new TH1D("pfCommonMHT_postCut","pfCommonMHT_postCut",200,0,100);
  pfCalcMHT_postCut = new TH1D("pfCalcMHT_postCut","pfCalcMHT_postCut",200,0,100);

}


std::ostream& AGBTrigEffPlots::Description(std::ostream& ostm) {
  ostm << "AGB Trigger efficiencies";
  return ostm;
}



bool AGBTrigEffPlots::Process(Event::Data& ev){


  double w = ev.GetEventWeight();
  double MHTx = 0.;
  double MHTy = 0.;
  for (std::vector<Event::Jet>::const_iterator j = ev.JD_Jets().begin(); j != ev.JD_Jets().end(); ++j){
    if(j->Pt()>=5.){
      MHTx -= j->Px();
      MHTy -= j->Py();
    }
  }
  
  double MHT = sqrt(MHTx*MHTx + MHTy*MHTy);
  
  pfCalcMHT->Fill(MHT,w);

  double commonMHT = ev.CommonMHT().Pt();

  pfCommonMHT->Fill(commonMHT,w);

  if(commonMHT > _TrigMHTcut){
    pfCommonMHT_postCut->Fill(commonMHT,w);
  }

  if(MHT > _TrigMHTcut){
    pfCalcMHT_postCut->Fill(MHT,w);
  }
  
  return true;

}
