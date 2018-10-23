#include <cstring>
#include "TStyle.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include <cmath>
using namespace std;

TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(double xl, double xh, double yl, double yh, double size = 0.04, const Int_t font = 42);
void SetPadMargin(TVirtualPad *pad, const Double_t bottomMargin = 0.12, const Double_t leftMargin = 0.12, const Double_t rightMargin = 0.05, const Double_t topMargin = 0.10);

TRandom3 *myRandom;

TH2F *hTpcTrackRes;
TH1F *hTrkResBin[400];
TF1 *hMuonPtEff[3];
TH1F *hMcJpsiPt;
TH1F *hInJpsiPt;
TH1F *hInJpsiCent;
TH1F *hOutJpsiPt[3];
TH1F *hOutJpsiCent[3];

//================================================
void toyMC_Cent()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  myRandom = new TRandom3();
  TDatime *clock = new TDatime();
  myRandom->SetSeed(clock->GetTime());


  // Get true refMult vs Ncoll distribution
  TH2F *hRefMultVsNcoll = 0x0;

  // Get different efficiency curve vs. refMult
  const int nEff = 2;
  const int nMult = 2;
  TF1 *funcEff[nEff];
  funcEff[0] = new TF1("funcEff_0", "pol1",0,1000);
  funcEff[0]->SetParameters(0.9, -2e-4);
  funcEff[1] = new TF1("funcEff_1", "pol0",0,1000);
  funcEff[1]->SetParameter(0,0.6);
  TCanvas *c = new TCanvas("cEff","cEff",800,600);
  TH1F *hplot = new TH1F("hplot",";Mult;",100,0,1000);
  hplot->SetYTitle("Efficiency");
  hplot->DrawCopy();
  for(int i=0; i<2; i++)
    {
      funcEff[i]->SetLineColor(i+1);
      funcEff[i]->SetLineWidth(1.5);
      funcEff[i]->Draw("sames");
    }
  TPaveText *title = GetTitleText("TPC tracking efficiency");
  title->Draw();
  TLegend *leg = new TLegend(0.4,0.2,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(funcEff[0], "Efficiency 1", "L");
  leg->AddEntry(funcEff[1], "Efficiency 2", "L");
  leg->Draw();

  // Get ngTrack distribution to mimic efficiency distribution
  
  TH1F *hMcMult[nMult];
  
}


//-----------------------------------------
TPaveText *GetTitleText(TString title, const Float_t size, const Int_t font)
{
  TPaveText* t1=new TPaveText(0.3530151,0.8968531,0.6532663,0.9965035,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,title.Data());
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}

//-----------------------------------------
TPaveText *GetPaveText(double xl, double xh, double yl, double yh, double size, const Int_t font)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}

//-----------------------------------------
void SetPadMargin(TVirtualPad *pad, const Double_t bottomMargin, const Double_t leftMargin, const Double_t rightMargin, const Double_t topMargin)
{
  if(!pad) return;
  pad->SetLeftMargin(leftMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetTopMargin(topMargin);
}
