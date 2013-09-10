#include "NLOTools.hh"
#include "AlphaT.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "GenMatrixBin.hh"
#include "Jet.hh"
#include "KinSuite.hh"
#include "Lepton.hh"
#include "Math/VectorUtil.h"
#include "Photon.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "Types.hh"
#include <sstream>
#include <string>
#include <iomanip>

#include <strstream>
#include <iostream>
#include <fstream>


#include "JetData.hh"
#include "CommonOps.hh"


using namespace Operation;

Double_t NLO::GetNLOCross(Double_t m0,Double_t m12, TString process, TString kFactorFile){

  // cout << " Trying to open NLO file" << endl;
  ifstream myfile(kFactorFile);

  if(myfile.is_open()){
    // cout << "File is open " << endl;
    std::string title;
    getline(myfile,title,'\n');
       // cout << " new " << endl;

    while(!myfile.eof()){

      getline(myfile,title,'\n');
   // cout << "title " << title << endl;

      int posm0 = title.find("m0=",0);
      int posm0_end = title.find(",",posm0);
      int length = posm0_end-posm0-3;
   // cout << "lentgth " << length << endl;
      if(length > 0){
        TString thism0 = title.substr(posm0+3,length);
    // cout << " thismo " << thism0 << endl;

        int posm12 = title.find("m1/2=",posm0_end);
        int posm12_end = title.find(",",posm12);
        length = posm12_end - posm12-5;
        TString thism12 = title.substr(posm12+5,length);
     // cout  << " thism12 " << thism12 << endl;

        TString thisprocess[10];

        int at = 0;
        int it = 0;
    // cout << " search m0 " << m0 << " search m12 " << m12 << endl;
        if(thism0.Atof() == m0 && thism12.Atof() == m12){
          for(int i = title.find("|",posm12_end); i != std::string::npos; i = title.find("|",i))
          {
            if(at > 0) thisprocess[it-1] = title.substr(at+2,i - at-3);
            at = i;
            i++;
            if(at > 0){
       // cout << " thism0 " << thism0.Atof() << " thism12 " << thism12.Atof();
      // cout << " searchprocess " << process << " process " << thisprocess[it-1].Atof() << " it " <<  endl;
            }
            it++;
          }
          if(process == "ng") return thisprocess[0].Atof();
          if(process == "ns") return thisprocess[1].Atof();
          if(process == "nn") return thisprocess[2].Atof();
          if(process == "ll") return thisprocess[3].Atof();
          if(process == "sb") return thisprocess[4].Atof();
          if(process == "ss") return thisprocess[5].Atof();
          if(process == "tb") return thisprocess[6].Atof();
          if(process == "bb") return thisprocess[7].Atof();
          if(process == "gg") return thisprocess[8].Atof();
          if(process == "sg") return thisprocess[9].Atof();

        }
      }

    }

  }
  return -1;

}

double NLO::ISRProducer(Event::Data& ev){

  // double x = 0;

  LorentzV newPart(0.,0.,0.,0.);

  // cout << " new  " << endl;
  for (std::vector<Event::GenObject>::const_iterator j = ev.GenParticles().begin();  j != ev.GenParticles().end(); ++j) {
    
    
    if((fabs((*j).GetMotherID()) == 21 || fabs((*j).GetMotherID()) <=6 ) && (*j).GetStatus() == 3 ){

        
      if( (fabs((*j).GetID()) >= 1000001 && fabs((*j).GetID()) <= 1000045) ||
	  (fabs((*j).GetID()) >= 2000001 && fabs((*j).GetID()) <= 2000016) || 
	  fabs((*j).GetID()) == 1000021 ){
      
	
	newPart = newPart+(*j);
	//	cout << " id " << (*j).GetID() << " status " << (*j).GetStatus() << " single pt " << (*j).Pt() << endl;
	
      }
    }
    
  }
  // cout << " newpart " << newPart.Pt() << endl;
  return newPart.Pt();
}


TString NLO::GetProcess(Event::Data& ev){

  bool verbose = false;

  int squarks = 0;
  int antisquarks = 0;
  int gluinos = 0;
  int sleptons = 0;
  int neutralinos = 0;
  int charginos = 0;
  int sbottoms = 0;
  int stops = 0;

  TString process = "ss";


  for (std::vector<Event::GenObject>::const_iterator j = ev.GenParticles().begin();  j != ev.GenParticles().end(); ++j) {


    if(fabs((*j).GetMotherID()) == 21 || fabs((*j).GetMotherID()) <=6  ){

      //select squarks
      if( ((*j).GetID() >= 1000001 && (*j).GetID() <= 1000004) ||
      ((*j).GetID() >= 2000001 && (*j).GetID() <= 2000004) ){
        squarks++;
      }
      //select antisquarks
      if( ((*j).GetID() <= -1000001 && (*j).GetID() >= -1000004) ||
      ((*j).GetID() <= -2000001 && (*j).GetID() >= -2000004) ){
        antisquarks++;
      }
      if( fabs((*j).GetID()) == 1000005 || fabs((*j).GetID()) == 2000005)sbottoms++;
      if( fabs((*j).GetID()) == 1000006 || fabs((*j).GetID()) == 2000006)stops++;


      //select gluinos
      if( fabs((*j).GetID()) == 1000021 )gluinos++;

      //select sleptons
      if( (fabs((*j).GetID()) >= 1000011 && fabs((*j).GetID()) <= 1000016) ||
        (fabs((*j).GetID()) >= 2000011 && fabs((*j).GetID()) <= 2000016))sleptons++;
      //select neutralinos
      if( fabs((*j).GetID()) == 1000022 || fabs((*j).GetID()) == 1000023 ||  fabs((*j).GetID()) == 1000025 || fabs((*j).GetID()) == 1000035 ||  fabs((*j).GetID()) == 1000045  )neutralinos++;

      if( fabs((*j).GetID()) == 1000024 || fabs((*j).GetID()) == 1000037  )charginos++;

    }
  }

  

  if(verbose == true)cout << " neutralinos " << neutralinos << " charginos " << charginos << " gluinos " << gluinos << " squarks " << squarks << " antisquarks " << antisquarks << " sleptons " << sleptons << " stops " << stops << " sbottoms " << sbottoms << endl;

  if((neutralinos + charginos) ==1 && gluinos == 1)process = "ng";
  if((neutralinos + charginos)==1 && (squarks + antisquarks) == 1) process = "ns";
  if(neutralinos + charginos == 2)process = "nn";
  if(sleptons == 2)process = "ll";
  if(squarks ==1 && antisquarks == 1)process = "sb";
  if(squarks == 2)process = "ss";
  if(stops ==1 && sbottoms == 1)process = "tb";
  if(sbottoms == 2)process = "bb";
  if(gluinos == 2)process = "gg";
  if(squarks + antisquarks == 1 && gluinos == 1)process = "sg";
 

  return process;
}



