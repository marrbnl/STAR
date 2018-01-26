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

void getJpsiEff(const int savePlot, const int saveHisto);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
void getJpsiEff(const double mass, const int nExpr, TH1F *hInputPt, TF1 *funcTpcTrkEff, TH1F *hMcJpsiPt, TH1F *hRcJpsiPt, int debug = 0);
TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size = 0.04, const Int_t font = 42);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);

const int year = 2014;
const char *run_type = "Run14_AuAu200";
const double jpsi_rapidity = 0.5;
const double pt1_cut = 1.5, pt2_cut = 1.3;
const int nPtBins_pt = 9;
const double xPtBins_pt[nPtBins_pt+1]  = {0,1,2,3,4,5,6,8,10,15};
const int nCentBins_pt = 5;
const char *cent_Name_pt[nCentBins_pt] = {"0-80","0-20","20-40","40-60","60-80"};
const char *cent_Title_pt[nCentBins_pt] = {"0080","0020","2040","4060","6080"};
const int nPtBins_npart = 2;
const double ptBins_low_npart[nPtBins_npart]  = {0,5};
const double ptBins_high_npart[nPtBins_npart] = {15,15};
const char *pt_Name_npart[nPtBins_npart] = {"0-15","5-15"};
const int nCentBins_npart[nPtBins_npart] = {9,8};
const char *cent_Name_npart[18] = {"0-80","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-80","0-10","10-20","20-30","30-40","40-50","50-60","60-80","60-70"};
const char *cent_Title_npart[18] = {"0080","0010","1020","2030","3040","4050","5060","6070","7080","0080","0010","1020","2030","3040","4050","5060","6080","6070"};
const TString sysType = "TpcTracking";
// const int nSys = 12;
// const TString sysName[nSys] = {"default",
// 			       "dcaUp","dcaDown","NHitsUp","NDedxUp",
// 			       "dcaUp_NHitsUp", "dcaUp_NDedxUp", "dcaUp_NHitsUp_NDedxUp",
// 			       "dcaDown_NHitsUp", "dcaDown_NDedxUp", "dcaDown_NHitsUp_NDedxUp",
// 			       "NHitsUp_NDedxUp"};
const int nSys = 5;
const TString sysName[nSys] = {"default", "dcaUp","dcaDown","NHitsUp","NDedxUp",};

// const TString sysType = "MuonPid";
// const int nSys = 27;
// const TString sysName[nSys] = {"default",
// 			       "dzUp","dzDown","dyUp","dyDown","nSigPiUp","nSigPiDown",
// 			       "dzUp_dyUp","dzUp_dyDown","dzUp_nSigPiUp","dzUp_nSigPiDown","dzUp_dyUp_nSigPiUp","dzUp_dyUp_nSigPiDown","dzUp_dyDown_nSigPiUp","dzUp_dyDown_nSigPiDown",
// 			       "dzDown_dyUp","dzDown_dyDown","dzDown_nSigPiUp","dzDown_nSigPiDown","dzDown_dyUp_nSigPiUp","dzDown_dyUp_nSigPiDown","dzDown_dyDown_nSigPiUp","dzDown_dyDown_nSigPiDown",
// 			       "dyUp_nSigPiUp","dyUp_nSigPiDown",
// 			       "dyDown_nSigPiUp","dyDown_nSigPiDown"};

const double muMass = 0.1057;
const double jpsiMass = 3.906;
const double pi = 3.1415926;
TFile *f;
TCanvas *c;
TRandom3 *myRandom;

//================================================
void make_SysTpcAndPid()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  myRandom = new TRandom3();
  myRandom->SetSeed(0);
  getJpsiEff(1,0);
}
//================================================
void getJpsiEff(const int savePlot, const int saveHisto)
{
  const int nCentBins      = 13; 
  const char* cent_Name[nCentBins]    = {"00-80","00-20","20-40","40-60","60-80","00-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
  const char* cent_Title[nCentBins]   = {"0080","0020","2040","4060","6080","0010","1020","2030","3040","4050","5060","6070","7080"};

  if(saveHisto) f = TFile::Open(Form("Rootfiles/%s.Sys.%s.root",run_type,sysType.Data()),"update");
  else          f = TFile::Open(Form("Rootfiles/%s.Sys.%s.root",run_type,sysType.Data()),"read");

  const int nTrkPtBins = 13;
  const double xTrkPtBins[14] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0,2.5,3.0,5.0,10,20};
  TString embed_name[2] = {"Mc","Tpc"};
  if(sysType=="MuonPid") embed_name[1] = "Pid";
  TH1F *hTrkPt[nCentBins][nSys][2];
  TH1F *hTrkEffVsPt[nCentBins][nSys];
  for(int k=0; k<nCentBins; k++)
    {
      for(int s=0; s<nSys; s++)
	{
	  for(int i=0; i<2; i++)
	    {
	      TH1F *h1tmp = (TH1F*)f->Get(Form("Sys%s_TrkEffVsPtEmb_%s_cent%s",sysName[s].Data(),embed_name[i].Data(),cent_Title[k]));
	      hTrkPt[k][s][i] = (TH1F*)h1tmp->Rebin(nTrkPtBins, Form("%s_rebin",h1tmp->GetName()), xTrkPtBins);
	    }
	  hTrkEffVsPt[k][s] = (TH1F*)hTrkPt[k][s][1]->Clone(Form("Sys%s_TrkEffVsPtEmb_cent%s",sysName[s].Data(),cent_Title[k]));
	  hTrkEffVsPt[k][s]->Divide(hTrkPt[k][s][0]);
	}
    }

  TF1 *funcTrkEff[nCentBins][nSys];
  for(int k=0; k<nCentBins; k++)
    {
      c = new TCanvas(Form("TpcTrkEff_cent%s",cent_Title[k]),Form("TpcTrkEff_cent%s",cent_Title[k]),1100,700);
      if(sysType=="TpcTracking") c->Divide(3,2);
      if(sysType=="MuonPid") c->Divide(7,4);
      for(int s=0; s<nSys; s++)
	{
	  funcTrkEff[k][s] = new TF1(Form("Sys%s_TrkEffVsPtEmb_cent%s_FitFunc",sysName[s].Data(),cent_Title[k]),"[0]-exp([1]+[2]*x+[3]*x*x)",1,20);
	  if(k==0) funcTrkEff[k][s]->SetParameters(0.7,-2,0.2);
	  else
	    {
	      funcTrkEff[k][s]->FixParameter(1, funcTrkEff[0][0]->GetParameter(1));
	      funcTrkEff[k][s]->FixParameter(2, funcTrkEff[0][0]->GetParameter(2));
	    }
	  hTrkEffVsPt[k][s]->Fit(funcTrkEff[k][s], "R0Q");
	  hTrkEffVsPt[k][s]->SetMarkerStyle(20);
	  hTrkEffVsPt[k][s]->GetXaxis()->SetRangeUser(0.2, 10);
	  if(sysType=="TpcTracking")
	    {
	      if(k==4 || k>=10) hTrkEffVsPt[k][s]->GetYaxis()->SetRangeUser(0.3,1.2);
	      else hTrkEffVsPt[k][s]->GetYaxis()->SetRangeUser(0.5,0.8);
	    }
	  else if(sysType=="MuonPid")
	    {
	      hTrkEffVsPt[k][s]->GetYaxis()->SetRangeUser(0,0.5);
	    }
	  hTrkEffVsPt[k][s]->SetTitle(";p_{T} (GeV/c);TPC efficiency");
	  c->cd(s+1);
	  hTrkEffVsPt[k][s]->Draw();
	  funcTrkEff[k][s]->SetLineColor(4);
	  funcTrkEff[k][s]->SetLineStyle(2);
	  funcTrkEff[k][s]->SetLineWidth(2);
	  funcTrkEff[k][s]->Draw("sames");
	  TPaveText *t1 = GetTitleText(sysName[s].Data(),0.06);
	  t1->Draw();
	  t1 = GetPaveText(0.6,0.7,0.75,0.8,0.07);
	  t1->AddText(Form("p0 = %2.2f",funcTrkEff[k][s]->GetParameter(0)));
	  t1->SetTextColor(4);
	  t1->Draw();
	}
      c->cd(1);
      TLegend *leg0 = new TLegend(0.2,0.2,0.7,0.4);
      leg0->SetBorderSize(0);
      leg0->SetFillColor(0);
      leg0->SetTextSize(0.06);
      leg0->SetHeader(Form("%s: %s%%",run_type,cent_Name[k]));
      leg0->AddEntry(hTrkEffVsPt[0][k], "Embedding", "P");
      leg0->AddEntry(funcTrkEff[0][k], "Fit function: p0-e^{p1*(x-p2)}", "L");
      leg0->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_TpcAndPid/%s_FitTrkEffVsPtEmb_cent%s.pdf",run_type,sysType.Data(),cent_Title[k]));
    }
  return;

  // get J/psi efficiency 
  const int nExpr = 1e6;
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  TH1F *hInputJpsi  = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hMcJpsiPt[nCentBins][nSys];
  TH1F *hRcJpsiPt[nCentBins][nSys];
  for(int k=0; k<nCentBins; k++)
    {
      for(int s=0; s<nSys; s++)
	{
	  hMcJpsiPt[k][s] = new TH1F(Form("Sys%s_JpsiPtMc_cent%s",sysName[s].Data(),cent_Title[k]),";p_{T} (GeV/c)",nPtBins_pt,xPtBins_pt);
	  hMcJpsiPt[k][s]->Sumw2();
	  hRcJpsiPt[k][s] = new TH1F(Form("Sys%s_JpsiPtRc_cent%s",sysName[s].Data(),cent_Title[k]),";p_{T} (GeV/c)",nPtBins_pt,xPtBins_pt);
	  hRcJpsiPt[k][s]->Sumw2();
	}
    }
  for(int k=0; k<nCentBins; k++)
    {
      for(int s=0; s<nSys; s++)
	{
	  printf("[i] %s, centrality %s%%,\n",sysName[s].Data(), cent_Name[k]);
	  getJpsiEff(jpsiMass, nExpr, hInputJpsi, funcTrkEff[k][s], hMcJpsiPt[k][s], hRcJpsiPt[k][s], 0);
	}
    }

  // efficeincy vs. pT
  TH1F *hJpsiEffVsPt[nCentBins_pt][nSys];
  for(int k=0; k<nCentBins_pt; k++)
    {
      for(int s=0; s<nSys; s++)
	{
	  hJpsiEffVsPt[k][s] = (TH1F*)hRcJpsiPt[k][s]->Clone(Form("Sys%s_JpsiEffVsPtEmb_cent%s",sysName[s].Data(),cent_Title_pt[k]));
	  hJpsiEffVsPt[k][s]->Divide(hMcJpsiPt[k][s]);
	}
    }
  c = draw1D(hJpsiEffVsPt[0][0]);

  // efficicy vs. centrality
  TH1F *hJpsiEffVsCent[nPtBins_npart][nSys];
  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int s=0; s<nSys; s++)
	{
	  hJpsiEffVsCent[i][s] = new TH1F(Form("Sys%s_JpsiEffVsCentEmb_Pt%1.0f",sysName[s].Data(),ptBins_low_npart[i]),"",nCentBins_npart[i],0,nCentBins_npart[i]);
	  for(int k=0; k<nCentBins_npart[i]; k++)
	    {
	      int index;
	      if(k==0) index = 0;
	      else index = k+nCentBins_pt -1; 
	      if(i==1 && k==nCentBins_npart[i]-1) index = 4;
	      int low_bin = hMcJpsiPt[index][s]->FindFixBin(ptBins_low_npart[i]+1e-4);
	      int up_bin  = hMcJpsiPt[index][s]->FindFixBin(ptBins_high_npart[i]-1e-4);
	      double nMcJpsi_err, nRcJpsi_err;
	      double nMcJpsi = hMcJpsiPt[index][s]->IntegralAndError(low_bin, up_bin, nMcJpsi_err);
	      double nRcJpsi = hRcJpsiPt[index][s]->IntegralAndError(low_bin, up_bin, nRcJpsi_err);
	      double eff = nRcJpsi/nMcJpsi;
	      hJpsiEffVsCent[i][s]->SetBinContent(k+1, eff);
	      hJpsiEffVsCent[i][s]->SetBinError(k+1, eff*nRcJpsi_err/nRcJpsi);
	    }
	}
    }
  c = draw1D(hJpsiEffVsCent[1][0]);


  if(saveHisto)
    {
      f->cd();
      for(int k=0; k<nCentBins_pt; k++)
	{
	  for(int s=0; s<nSys; s++)
	    {
	      hJpsiEffVsPt[k][s]->Write("",TObject::kOverwrite);
	    }
	}

      for(int i=0; i<nPtBins_npart; i++)
	{
	  for(int s=0; s<nSys; s++)
	    {
	      hJpsiEffVsCent[i][s]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//-------------------------------------------------------
void getJpsiEff(const double mass, const int nExpr, TH1F *hInputPt, TF1 *funcTpcTrkEff, TH1F *hMcJpsiPt, TH1F *hRcJpsiPt, int debug)
{
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  = myRandom->Uniform(0,15);
      double weight = hInputPt->GetBinContent(hInputPt->FindFixBin(mc_pt));
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;

      double pt1  = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      if(debug) printf("daugther 1: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt1,eta1,phi1);
      
      double pt2  = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);

      hMcJpsiPt->Fill(parent.Pt(), weight);
      if(fabs(eta1) > 0.5 || fabs(eta2) > 0.5) continue;
      if(pt1<1.3 || pt2<1.3) continue;
      if(pt1<1.5 && pt2<1.5) continue;
      double eff1 = funcTpcTrkEff->Eval(pt1);
      double eff2 = funcTpcTrkEff->Eval(pt2);
      hRcJpsiPt->Fill(parent.Pt(), weight*eff1*eff2);
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
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size, const Int_t font)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
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
TCanvas *draw1D(TH1 *h, TString hTitle, Bool_t setLog, Bool_t drawP, const Float_t size, const TString drawOpt, const Int_t titleFont,
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
