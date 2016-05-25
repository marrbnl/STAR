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
#include "TLegend.h"

#include <cmath>
using namespace std;
const double pi = 3.1415926;

TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);
TRandom3 *myRandom;

//================================================
void ana_Smear()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  // get the random number generator 
  myRandom = new TRandom3();
  TDatime *clock = new TDatime();
  myRandom->SetSeed(clock->GetTime());

  TFile *fin = TFile::Open(Form("Xinjie.Run14.AuAu200.UpsilonSmear.root"),"read");
  // pt resolution from embedding
  TH2F *hTrkResVsPt = (TH2F*)fin->Get("PrimTrkRes_vs_TruePt_cent0060");
  const int nHistos = hTrkResVsPt->GetNbinsX();
  TH1F *hTrkResBin[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTrkResVsPt->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  // single muon efficiency that includes tracking, matching, trigger efficiencies
  // For this smearing exercise, the only thing that matters is the shape of the efficiency
  // McTrkPtEff_cent0060_scaled = 5 * McTrkPtEff_cent0060
  TF1 *funcEff = (TF1*)fin->Get("Fit_McTrkPtEff_cent0060_scaled"); // this function does not fit the efficiency, but it is ok enough for now

  // Upslion pT distribution from pp
  TH1F *hMcUpsilonPt = (TH1F*)fin->Get("Upsilon_pp200_midRap");

  // Toy Monte Carlo
  const double upsilonMass = 9.46; // Y(1S) mass
  const double muMass = 0.1057; // muon mass
  const double pt1_cut = 1.5, pt2_cut = 1.5; 
  const double shift = 0; // for Run14 data, no shift is needed
  const double sigma = 0.0089; // additional smearing try to change to 0.0069 and 0.0108 for uncertainty

  TH2F *hRcMassVsPt = new TH2F("hRcMassVsPt","Mass distribution of reconstructed #Upsilon;p_{T} (GeV/c);M_{#mu#mu} (GeV/c^2)",20,0,10,600,8,14);
  const int nExpr = 1e4; // # of events for MC
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  = myRandom->Uniform(0,10);
      double weight = hMcUpsilonPt->GetBinContent(hMcUpsilonPt->FindFixBin(mc_pt));
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.8, 0.8);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+upsilonMass*upsilonMass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,upsilonMass);
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;

      double pt1  = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      
      double pt2  = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();

      if(pt1<0.5 || pt2<0.5) continue; // apply basic cuts to speed things up

      // tracking efficiency from embedding
      double eff1 = funcEff->Eval(pt1);
      double eff2 = funcEff->Eval(pt2);
      double pro1 = myRandom->Uniform(0., 1.);
      double pro2 = myRandom->Uniform(0., 1.);
      if(pro1>eff1 || pro2>eff2) continue;

      // momentum resolution from ebemdding
      int mom_index1 = hTrkResVsPt->GetXaxis()->FindBin(pt1)-1;
      int mom_index2 = hTrkResVsPt->GetXaxis()->FindBin(pt2)-1;
      if(mom_index1>=nHistos) mom_index1 = nHistos-1;
      if(mom_index2>=nHistos) mom_index2 = nHistos-1;
      double dpt1 = hTrkResBin[mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[mom_index2]->GetRandom();
      double emb_pt1 = (1-dpt1) * pt1;
      double emb_pt2 = (1-dpt2) * pt2;

      // additional smearing & shifting
      double rc_pt1 = emb_pt1 * myRandom->Gaus(1+shift/sqrt(emb_pt1),sqrt(emb_pt1)*sigma);
      double rc_pt2 = emb_pt2 * myRandom->Gaus(1+shift/sqrt(emb_pt2),sqrt(emb_pt2)*sigma);

      double leadPt = rc_pt1 > rc_pt2 ? rc_pt1 : rc_pt2;
      double subPt  = rc_pt1 > rc_pt2 ? rc_pt2 : rc_pt1;
      if(leadPt<pt1_cut || subPt<pt2_cut) continue;

      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hRcMassVsPt->Fill(rc_parent.Pt(),rc_parent.M(),weight);
    }
  TCanvas *c = new TCanvas();
  hRcMassVsPt->Draw("colz");
}

//-------------------------------------------------------
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter)
{
  float betax = parent.Px()/parent.E();
  float betay = parent.Py()/parent.E();
  float betaz = parent.Pz()/parent.E();
  daughter.Boost(betax,betay,betaz);
  return daughter;
}

//-------------------------------------------------------
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass) 
{
  Double_t e = parent.M()/2.;
  Double_t p = sqrt(e*e-dmass*dmass);
  Double_t costheta = myRandom->Uniform(-1.,1.);
  Double_t phi = myRandom->Uniform(0,TMath::Pi()*2);
  Double_t pz = p*costheta;
  Double_t px = p*sqrt(1.-costheta*costheta)*cos(phi);
  Double_t py = p*sqrt(1.-costheta*costheta)*sin(phi);
  TLorentzVector daughter(px,py,pz,e);
  return myBoost(parent,daughter);
}
