#include "TestMerge.h"
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TClass.h"
#include "TTree.h"
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

// -----------------------------------------------------------------------------
//
TestMerge::TestMerge( const edm::ParameterSet& pset ) 
  : files_(),
    file_( pset.getUntrackedParameter<std::string>("OutputFile") ),
    histos_(),
    histos2d_(),
    lumi_( pset.getUntrackedParameter<double>("NormalisedLumi") )
{
  files();
  files( pset );
}

// -----------------------------------------------------------------------------
//
TestMerge::Data::Data( std::string name, float xsec ) 
  : name_(name),
    xSec_(xsec)
{;}

// -----------------------------------------------------------------------------
//
TestMerge::Data::Data() 
  : name_(""),
    xSec_(0.)
{;}

// -----------------------------------------------------------------------------
//
void TestMerge::files() {
  map<TFile*,Data>::iterator ii = files_.begin();
  map<TFile*,Data>::iterator jj = files_.end();
  for ( ; ii != jj; ++ii ) { 
    if ( ii->first ) { 
      if ( ii->first->IsOpen() ) { ii->first->Close(); }
      delete ii->first; 
    }
  }
  files_.clear();
}

// -----------------------------------------------------------------------------
//
void TestMerge::files( const PSet& pset ) {
  VPSet psets = pset.getUntrackedParameter<VPSet>("InputFiles");
  VPSet::const_iterator ii = psets.begin();
  VPSet::const_iterator jj = psets.end();
  for ( ; ii != jj; ++ii ) {
    std::string name = ii->getUntrackedParameter<std::string>("FileName");
    // Search for "name"
    map<TFile*,Data>::iterator iii = files_.begin();
    map<TFile*,Data>::iterator jjj = files_.end();
    bool found = false;
    while ( !found && iii != jjj ) { 
      if ( name == iii->second.name_ ) { found = true; } 
      ++iii;
    }
    // If not found, create entry
    if ( !found ) {
      TFile* file = new TFile( name.c_str() );
      files_[ file ] = Data( name, ii->getUntrackedParameter<double>("XSection") );
    }
  }
}

// -----------------------------------------------------------------------------
//
void TestMerge::beginJob( const edm::EventSetup& ) {

  // Create output file
  TFile* output_file = new TFile( file_.c_str(), "RECREATE" );
  
  // Merge input files
  TList* input_files = new TList();
  std::map<TFile*,Data>::const_iterator ii = files_.begin();
  std::map<TFile*,Data>::const_iterator jj = files_.end();
  for ( ; ii != jj; ++ii ) { input_files->Add( ii->first ); }
  merge( output_file, input_files );
  
  // Close input files
  files();
  delete input_files;

  // Close output files
  output_file->Close();
  delete output_file;
  
}

// -----------------------------------------------------------------------------
//
void TestMerge::merge( TDirectory* target, TList* sourcelist ) {
  
  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  
  TFile *first_source = (TFile*)sourcelist->First();
  cout << "Target name: " << first_source->GetName() << endl;
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it

      //      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;

      { // Calculate event weight
	std::map<TFile*,Data>::const_iterator iter = files_.find( first_source );
	if ( iter != files_.end() ) {
	  float events = -1.;
	  TH1F* histo = (TH1F*)first_source->Get("test/Common/CutFlow_Efficiency");
	  if ( histo ) { events = histo->GetBinContent(1); }
	  else { cout << "No histo!" << endl; }
	  if ( events > 0. && iter->second.xSec_ > 0. ) {
	    h1->Sumw2();
	    h1->Scale( lumi_ / ( (float)events / iter->second.xSec_ ) );
	  } else { cout << "Null values!" << endl; }
	}
      }
      
      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
        
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
	  TH1 *h2 = (TH1*)key2->ReadObj();
	  
	  { // Calculate event weight
	    std::map<TFile*,Data>::const_iterator iter = files_.find( nextsource );
	    if ( iter != files_.end() ) {
	      float events = -1.;
	      TH1F* histo = (TH1F*)nextsource->Get("test/Common/CutFlow_Efficiency");
	      if ( histo ) { events = histo->GetBinContent(1); }
	      else { cout << "No histo!" << endl; }
	      if ( events > 0. && iter->second.xSec_ > 0. ) {
		h2->Sumw2();
		h2->Scale( lumi_ / ( (float)events / iter->second.xSec_ ) );
	      } else { cout << "Null values!" << endl; }
	    }
	  }
	  
	  h1->Add( h2 );
	  delete h2;
        }

        nextsource = (TFile*)sourcelist->After( nextsource );
      }
    }
    else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
      
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();

      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
      while ( nextsource ) {
     	  
	globChain->Add(nextsource->GetName());
	nextsource = (TFile*)sourcelist->After( nextsource );
      }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      merge( newdir, sourcelist );

    } else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " 
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( "TTree" ))
	globChain->Merge(target->GetFile(),0,"keep");
      else
	obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);

}
  
// -----------------------------------------------------------------------------
//
TH1D* TestMerge::histo( const std::string& histogram_name ) {
  std::map<std::string,TH1D*>::const_iterator ii = histos_.find(histogram_name);
  if ( ii != histos_.end() ) { return ii->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}
  
// -----------------------------------------------------------------------------
//
TH2D* TestMerge::histo2d( const std::string& histogram_name ) {
  std::map<std::string,TH2D*>::const_iterator jj = histos2d_.find(histogram_name);
  if ( jj != histos2d_.end() ) { return jj->second; }
  edm::LogWarning("TEST") << "Cannot find string: " << histogram_name;
  return 0;
}

// -----------------------------------------------------------------------------
// 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestMerge);
