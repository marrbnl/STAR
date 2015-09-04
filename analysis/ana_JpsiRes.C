#include "TStyle.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TPaveText.h"
using namespace std;

const double pi = 3.1415926;
const double pt1_cut = 1.5, pt2_cut = 1.0;
const double low_mass = 2.9;
const double high_mass = 3.3;
const char *charge_name[3] = {"","_pos","_neg"};
const double jpsiMass = 3.097;
const double muMass = 0.1057;
const int nCentBins = 4;
const int centBins_low[nCentBins]  = {5,13,9,5};
const int centBins_high[nCentBins] = {16,16,12,8};
const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60"};
const char *cent_Title[nCentBins] = {"0060","0020","2040","4060"};

//const char *run_type = "Run13_pp500";

const char *run_type = "Run14_AuAu200";
const int year = 2014;
TRandom3 *myRandom;

void smear(const int icent);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);

TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TCanvas *draw1D(TH1F *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);
TCanvas *draw2D(TH2F *h, const TString hTitle = "" , const Float_t size = 0.04, const Bool_t logz = kTRUE, const char *drawOption = "colz");


//================================================
void ana_JpsiRes()
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

  smear(0);
}

//================================================
void smear(const int icent)
{
  // check with the embedded data
  TFile *fEmbed = TFile::Open("Rootfiles/Run14.AuAu200.TrkEff.all.root","read");

  TString name[3] = {"MCinput","TPCreco","MTDreco"};
  TH1F *hJpsiPt[3], *hJpsiMass[3];
  for(int i=0; i<3; i++)
    {
      hJpsiPt[i]   = (TH1F*)fEmbed->Get(Form("%s_Jpsi_pT_%s",name[i].Data(),cent_Title[icent]));
      hJpsiMass[i] = (TH1F*)fEmbed->Get(Form("%s_Jpsi_mass_%s",name[i].Data(),cent_Title[icent]));   
      hJpsiMass[i]->GetXaxis()->SetRangeUser(2.6,3.6);
      hJpsiMass[i]->SetMarkerStyle(21);
    }

  TCanvas *c = draw1D(hJpsiPt[0],Form("%s: p_{T} distribution of J/#psi",name[0].Data()),kFALSE,kTRUE);


  // tracking efficiency
  TFile *fEff = TFile::Open("Rootfiles/Run14.AuAu200.TrkEff.root","read");
  TH3F *hTpcEff[2];
  TH2F *hMtdEff[2][2];
  for(int i=0; i<2; i++)
    {
      hTpcEff[i] = (TH3F*)fEff->Get(Form("McTrkPtEtaPhiTpc%s_%s_Eff",charge_name[i+1],cent_Title[icent]));
      hMtdEff[0][i] = (TH2F*)fEff->Get(Form("McTrkPtEtaMtdEff%s_3tray_%s",charge_name[i+1],cent_Title[icent]));
      hMtdEff[1][i] = (TH2F*)fEff->Get(Form("McTrkPtEtaMtdEff%s_5tray_%s",charge_name[i+1],cent_Title[icent]));
    }

  // track resolution
  TH2F *hTrkResVsPt = (TH2F*)fEff->Get(Form("PrimTrkRes_vs_TruePt_%s",cent_Title[icent]));
  c = draw2D(hTrkResVsPt);
  TF1 *funcTrkRes = (TF1*)fEff->Get(Form("FuncPrimTrkRes_vs_TruePt_%s",cent_Title[icent]));
  TH2F *hTrkResVsPtScale = (TH2F*)hTrkResVsPt->Clone(Form("%s_scaled",hTrkResVsPt->GetName()));
  hTrkResVsPtScale->Reset();
  for(int binx=1; binx<=hTrkResVsPt->GetNbinsX(); binx++)
    {
      double bincenter = hTrkResVsPt->GetXaxis()->GetBinCenter(binx);
      double resolution = funcTrkRes->Eval(bincenter);
      double scale = 0.01/resolution;
      for(int biny=1; biny<=hTrkResVsPt->GetNbinsY(); biny++)
	{
	  int newbiny = hTrkResVsPt->GetYaxis()->FindFixBin(hTrkResVsPt->GetYaxis()->GetBinCenter(biny)*scale);
	  hTrkResVsPtScale->SetBinContent(binx,newbiny,hTrkResVsPt->GetBinContent(binx,biny));
	  hTrkResVsPtScale->SetBinError(binx,newbiny,hTrkResVsPt->GetBinError(binx,biny));
	}
    }
  c = draw2D(hTrkResVsPtScale);
  //TH1F *hTrkRes = (TH1F*)hTrkResVsPtScale->ProjectionY(Form("hTrkRes_cent%s",cent_Title[icent]),hTrkResVsPtScale->GetXaxis()->FindFixBin(0.5),-1);
  const double sample_pt = 9.9;
  int selectbin = hTrkResVsPt->GetXaxis()->FindFixBin(sample_pt);
  TH1F *hTrkRes = (TH1F*)hTrkResVsPt->ProjectionY(Form("hTrkRes_cent%s",cent_Title[icent]),selectbin,selectbin);
  c = draw1D(hTrkRes,"");

  const int nHistos = hTrkResVsPt->GetNbinsX();
  TH1F *hTrkResBin[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTrkResVsPt->ProjectionY(Form("hTrkRes_Bin%d_cent%s",i+1,cent_Title[icent]),i+1,i+1);
    }


  TH1F *hMcJpsiPt   = new TH1F(Form("hMcJpsiPt_%s",cent_Title[icent]),Form("p_{T} distribution of input MC J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
  TH1F *hRcJpsiPt[2];
  TH1F *hRcJpsiMass[2];
  for(int i=0; i<2; i++)
    {
      hRcJpsiPt[i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i+1].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
      hRcJpsiMass[i] = new TH1F(Form("hRcJpsiMass_%s_%s",name[i+1].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),1400,0,14);
    }

  const bool debug = 0;
  const int nExpr = 1e6;
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("\n+++ %d +++\n",i);
      double mc_pt =  myRandom->Uniform(0., 10.);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.8, 0.8);
      hMcJpsiPt->Fill(mc_pt);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+jpsiMass*jpsiMass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,jpsiMass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;
      
      int charge1 = myRandom->Uniform(-0.5,0.5)>0? 1 : -1;
      int charge2 = charge1 * -1;
      
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      if(debug) printf("daugther 1: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt1,eta1,phi1);
      
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);
      
      // tracking efficiency
      if(fabs(eta1)>0.8 || fabs(eta2)>0.8) continue;
      if(fabs(pt1)<0.5  || fabs(pt2)<0.5) continue;

      int index1 = (charge1 > 0) ? 0 : 1;
      double phi1_new = phi1;
      if(phi1_new<0) phi1_new += 2*pi;
      phi1_new += pi/12;
      double probility1 = hTpcEff[index1]->GetBinContent(hTpcEff[index1]->GetXaxis()->FindFixBin(pt1),
							 hTpcEff[index1]->GetYaxis()->FindFixBin(eta1),
							 hTpcEff[index1]->GetZaxis()->FindFixBin(phi1_new));
      if(debug) printf("Efficiency 1 = %3.2f\n",probility1);
      if( myRandom->Uniform(0., 1.) >  probility1) continue;

      int index2 = (charge2 > 0) ? 0 : 1;
      double phi2_new = phi2;
      if(phi2_new<0) phi2_new += 2*pi;
      phi2_new += pi/12;
      double probility2 = hTpcEff[index2]->GetBinContent(hTpcEff[index2]->GetXaxis()->FindFixBin(pt2),
							 hTpcEff[index2]->GetYaxis()->FindFixBin(eta2),
							 hTpcEff[index2]->GetZaxis()->FindFixBin(phi2_new));
      if(debug) printf("Efficiency 2 = %3.2f\n",probility2);
      if( myRandom->Uniform(0., 1.) >  probility2) continue;
      

      // momentum resolution
      // double rc_pt1 = myRandom->Gaus(pt1,funcTrkRes->Eval(pt1)*pt1);
      // double rc_pt2 = myRandom->Gaus(pt2,funcTrkRes->Eval(pt2)*pt2);
      // double rc_pt1 = (1 + hTrkRes->GetRandom() * funcTrkRes->Eval(pt1) / funcTrkRes->Eval(sample_pt)) * pt1;
      // double rc_pt2 = (1 + hTrkRes->GetRandom() * funcTrkRes->Eval(pt2) / funcTrkRes->Eval(sample_pt)) * pt2;
      
      double rc_pt1 = (1-hTrkResBin[hTrkResVsPt->GetXaxis()->FindFixBin(pt1)-1]->GetRandom()) * pt1;
      double rc_pt2 = (1-hTrkResBin[hTrkResVsPt->GetXaxis()->FindFixBin(pt2)-1]->GetRandom()) * pt2;

      if(debug) printf("rc daug 1: pt = %3.2f\n",rc_pt1);
      if(debug) printf("rc daug 2: pt = %3.2f\n",rc_pt2);
      if(rc_pt1<1.0 || rc_pt2<1.0) continue;
      if(rc_pt1<1.5 && rc_pt2<1.5) continue;
      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hRcJpsiPt[0]->Fill(rc_parent.Pt());
      hRcJpsiMass[0]->Fill(rc_parent.M());

      // MTD matching
      int tray_index = 1;
      if(phi1<3*pi/15 && phi1>12*pi/15) tray_index = 0;
      double mtd_eff_1 = hMtdEff[tray_index][index1]->GetBinContent(hMtdEff[tray_index][index1]->GetXaxis()->FindFixBin(pt1),
								    hMtdEff[tray_index][index1]->GetYaxis()->FindFixBin(eta1));
      if( myRandom->Uniform(0., 1.) >  mtd_eff_1) continue;
      
      double mtd_eff_2 = hMtdEff[tray_index][index2]->GetBinContent(hMtdEff[tray_index][index2]->GetXaxis()->FindFixBin(pt2),
								    hMtdEff[tray_index][index2]->GetYaxis()->FindFixBin(eta2));
      if( myRandom->Uniform(0., 1.) >  mtd_eff_2) continue;
      hRcJpsiPt[1]->Fill(rc_parent.Pt());
      hRcJpsiMass[1]->Fill(rc_parent.M());
      
    }

  for(int i=0; i<2; i++)
    {
      hJpsiPt[i+1]->Sumw2();
      hJpsiPt[i+1]->Scale(1./hJpsiPt[i+1]->Integral());
      hJpsiPt[i+1]->SetMarkerStyle(21);
      c = draw1D(hJpsiPt[i+1],Form("%s: p_{T} distribution of J/#psi",name[i+1].Data()),kFALSE,kTRUE);
      hRcJpsiPt[i]->Sumw2();
      hRcJpsiPt[i]->Scale(1./hRcJpsiPt[i]->Integral());
      hRcJpsiPt[i]->SetLineColor(2);
      hRcJpsiPt[i]->Draw("sames HIST");

      hJpsiMass[i+1]->Sumw2();
      hJpsiMass[i+1]->Scale(1./hJpsiMass[i+1]->Integral());
      c = draw1D(hJpsiMass[i+1],Form("%s: mass distribution of J/#psi",name[i+1].Data()),kTRUE,kTRUE);
      hRcJpsiMass[i]->Sumw2();
      hRcJpsiMass[i]->Scale(1./hRcJpsiMass[i]->Integral());
      hRcJpsiMass[i]->SetLineColor(2);
      hRcJpsiMass[i]->Draw("sames HIST");

      TH1F *hMassRatio = (TH1F*)hRcJpsiMass[i]->Clone(Form("hMassRatio_%s_%s",name[i+1].Data(),cent_Title[icent]));
      hMassRatio->Divide(hJpsiMass[i+1]);
      hMassRatio->GetXaxis()->SetRangeUser(2.8,3.4);
      hMassRatio->GetYaxis()->SetRangeUser(0,2);
      hMassRatio->SetMarkerStyle(20);
      hMassRatio->SetMarkerColor(2);
      c = draw1D(hMassRatio,Form("%s: ratio of mass distribution of J/#psi;M_{#mu#mu} (GeV/c^{2});ToyModel/Embed",name[i+1].Data()),kFALSE,kTRUE);
    }
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

//-----------------------------------------
TCanvas *draw1D(TH1F *h, TString hTitle, Bool_t setLog, Bool_t drawP, const Float_t size, const TString drawOpt, const Int_t titleFont,
		const Int_t wh, const Int_t ww)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),wh,ww);
  if(setLog) gPad->SetLogy();
  if(hTitle.Length()>0) h->SetTitle(hTitle);
  TPaveText *t1 = GetTitleText(h->GetTitle(),size,titleFont);
  h->SetTitle("");
  if(drawOpt.Length()>0)
    h->DrawCopy(drawOpt.Data());
  else
    {
      if(drawP)
	h->DrawCopy("PE");
      else
	h->DrawCopy("HIST");
    }
  t1->Draw();
  return c;
}

//-----------------------------------------
TCanvas *draw2D(TH2F *h, const TString hTitle, const Float_t size, const Bool_t logz, const char *drawOption)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),800,600);
  if(hTitle.Length()>0)
    h->SetTitle(hTitle);
  if(logz) gPad->SetLogz();
  TPaveText *t1 = GetTitleText(h->GetTitle(),size);
  h->SetTitle("");
  h->DrawCopy(drawOption);
  t1->Draw();
  return c;
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
