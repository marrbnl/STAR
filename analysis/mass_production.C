#include "TFile.h"

const TString wdir = "/Users/admin/Work/STAR/analysis";

//================================================
void mass_production()
{
  EventQA();
}

//================================================
void EventQA()
{
  //vertex();
  //track();
  //mtdhits();
  matching();
}

//================================================
void vertex()
{
  gROOT->LoadMacro(Form("%s/qa_vertex.C++",wdir.Data()));
  qa_vertex(1);
}

//================================================
void track()
{
  gROOT->LoadMacro(Form("%s/qa_track.C",wdir.Data()));
  qa_track(1);
}

//================================================
void mtdhits()
{
  gROOT->LoadMacro(Form("%s/qa_MtdHit.C",wdir.Data()));
  qa_MtdHit(1);
}

//================================================
void matching()
{
  gROOT->LoadMacro(Form("%s/ana_Match.C",wdir.Data()));
  ana_Match(1);
}
