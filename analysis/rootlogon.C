#include "/Users/admin/Work/STAR/util/defs.h"

//R__ADD_INCLUDE_PATH(/Users/admin/Work/STAR/util)
//#include "drawHistos.h"

void rootlogon(){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  const TString dir = "/Users/admin/Work/STAR/util";

#if defined(__CINT__)
  gROOT->LoadMacro(Form("%s/drawHistos.C",dir.Data()));
  gROOT->LoadMacro(Form("%s/HistUtility.C",dir.Data()));
#endif
} 
