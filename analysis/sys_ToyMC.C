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

#define YEAR 2014
#if (YEAR==2014)
const char *run_type = "Run14_AuAu200";
#elif (YEAR==2013)
char *run_type = "Run13_pp500";
#elif (YEAR==2015)
char *run_type = "Run15_pp200";
#endif

const double pi = 3.1415926;
const double muMass = 0.1057;
const int year = YEAR;
const int nPart = 4;
const char *part_name[nPart]  = {"Jpsi","Ups1S","Ups2S","Ups3S"};
const char *part_title[nPart] = {"Jpsi","Y(1S)","Y(2S)","Y(3S)"};
const double part_mass[nPart] = {3.097, 9.46, 10.023, 10.355};

TRandom3 *myRandom;
void anaSys(const int saveHisto = 0);
void plotSys_MtdTrigEff(const int savePlot = 0, const int saveHisto = 0);
void makeHisto(TString name, const double mass, const int nExpr = 1e4);
void toyMC(const double mass, const int nExpr, const int debug = 0);

TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);

void ScaleHistoTitle(const TH1 *h, 
                     const Double_t xTitleSize = 28, const Double_t xTitleOffset = 0.9, const Double_t xLabelSize = 20,
                     const Double_t yTitleSize = 28, const Double_t yTitleOffset = 0.9, const Double_t yLabelSize = 20,
                     const Int_t font = 42);
TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(double xl, double xh, double yl, double yh, double size = 0.04, const Int_t font = 42);
TLine *GetLine(double xl, double yl, double xh, double yh, Color_t color=2, Width_t width=2, Style_t style=2);
void SetPadMargin(TVirtualPad *pad, const Double_t bottomMargin = 0.12, const Double_t leftMargin = 0.12, const Double_t rightMargin = 0.05, const Double_t topMargin = 0.10);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);
TCanvas *draw2D(TH2 *h, const TString hTitle = "" , const Float_t size = 0.04, const Bool_t logz = kTRUE, const char *drawOption = "colz");
TCanvas *drawGraph(TGraph *h, const TString hTitle = "",Bool_t setLog = kFALSE, const Float_t size = 0.04, const char *drawOption = "AP", const Int_t font = 62);
void offset_x_with_asym_error(TGraphAsymmErrors* g, Double_t xoff);

TH2F *hTpcTrackRes;
TH1F *hTrkResBin[400];
TF1 *hMuonPtEff[3];
TH1F *hMcJpsiPt;
TH1F *hInJpsiPt;
TH1F *hInJpsiCent;
TH1F *hOutJpsiPt[3];
TH1F *hOutJpsiCent[3];

//================================================
void sys_ToyMC()
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

  anaSys(1);
  //plotSys_MtdTrigEff(1, 1);
}

//================================================
void plotSys_MtdTrigEff(const int savePlot, const int saveHisto)
{
  const int nLumi = 2;
  const int nSys = 2;
  const char *name_lumi[nLumi] = {"prod_high","prod_low"};
  const char *sys_legend[nSys] = {"M1: uncorr. stat. err.","M2: corr. stat. err."};
  TGraphAsymmErrors *gEffSysVsPt[nPart][nLumi][nSys];
  TGraphAsymmErrors *gEffSysVsCent[nPart][nLumi][nSys];  
  TH1F *hEffVsPt[nPart][nLumi][nSys][3];
  TH1F *hEffVsCent[nPart][nLumi][nSys][3];

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.JpsiMuon.sys.root",run_type),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s.JpsiMuon.sys.root",run_type),"read");
  for(int i=0; i<nPart; i++)
    {
      for(int k=0; k<nLumi; k++)
	{
	  for(int s=0; s<nSys; s++)
	    {
	      hEffVsPt[i][k][s][0] = (TH1F*)fin->Get(Form("%sEffVsPt_MtdTrigEff_BinCount_%s_FitFunc",part_name[i],name_lumi[k]));
	      hEffVsPt[i][k][s][1] = (TH1F*)fin->Get(Form("%sEffVsPt_MtdTrigEff_BinCount_%s_Sys%d_up",part_name[i],name_lumi[k],s+1));
	      hEffVsPt[i][k][s][2] = (TH1F*)fin->Get(Form("%sEffVsPt_MtdTrigEff_BinCount_%s_Sys%d_low",part_name[i],name_lumi[k],s+1));
	      int nBins = hEffVsPt[i][k][s][0]->GetNbinsX();
	      gEffSysVsPt[i][k][s] = new TGraphAsymmErrors(nBins);
	      gEffSysVsPt[i][k][s]->SetName(Form("%sEffVsPt_MtdTrigEff_BinCount_%s_Sys%d",part_name[i],name_lumi[k],s+1));
	      for(int bin=1; bin<=nBins; bin++)
		{
		  double eh = 0, el = 0;
		  for(int e=0; e<2; e++)
		    {
		      double diff = hEffVsPt[i][k][s][e+1]->GetBinContent(bin)/hEffVsPt[i][k][s][0]->GetBinContent(bin) - 1;
		      if(e==0) eh = fabs(diff);
		      else     el = fabs(diff);
		    }
		  gEffSysVsPt[i][k][s]->SetPoint(bin-1, hEffVsPt[i][k][s][0]->GetBinCenter(bin), 1);
		  gEffSysVsPt[i][k][s]->SetPointError(bin-1, 0, 0, el, eh);
		}

	      hEffVsCent[i][k][s][0] = (TH1F*)fin->Get(Form("%sEffVsCent_MtdTrigEff_BinCount_%s_FitFunc",part_name[i],name_lumi[k]));
	      hEffVsCent[i][k][s][1] = (TH1F*)fin->Get(Form("%sEffVsCent_MtdTrigEff_BinCount_%s_Sys%d_up",part_name[i],name_lumi[k],s+1));
	      hEffVsCent[i][k][s][2] = (TH1F*)fin->Get(Form("%sEffVsCent_MtdTrigEff_BinCount_%s_Sys%d_low",part_name[i],name_lumi[k],s+1));
	      nBins = hEffVsCent[i][k][s][0]->GetNbinsX();
	      gEffSysVsCent[i][k][s] = new TGraphAsymmErrors(nBins);
	      gEffSysVsCent[i][k][s]->SetName(Form("%sEffVsCent_MtdTrigEff_BinCount_%s_Sys%d",part_name[i],name_lumi[k],s+1));
	      for(int bin=1; bin<=nBins; bin++)
		{
		  double eh = 0, el = 0;
		  for(int e=0; e<2; e++)
		    {
		      double diff = hEffVsCent[i][k][s][e+1]->GetBinContent(bin)/hEffVsCent[i][k][s][0]->GetBinContent(bin) - 1;
		      if(e==0) eh = fabs(diff);
		      else     el = fabs(diff);
		    }
		  gEffSysVsCent[i][k][s]->SetPoint(bin-1, hEffVsCent[i][k][s][0]->GetBinCenter(bin), 1);
		  gEffSysVsCent[i][k][s]->SetPointError(bin-1, 0, 0, el, eh);
		}
	    }
	}
    }

  TH1F *hplot = new TH1F("hplot",";p_{T} (GeV/c);Sys. Uncert.",100,0,15);
  TH1F *hplot2 = new TH1F("hplo2t",";;Sys. Uncert.",2,0,2);
  hplot2->GetXaxis()->SetBinLabel(1, "p_{T} > 0 GeV/c");
  hplot2->GetXaxis()->SetBinLabel(2, "p_{T} > 5 GeV/c");
  TCanvas *c1[nPart];
  for(int i=0; i<nPart; i++)
    {
      c1[i] = new TCanvas(Form("%s_c1",part_name[i]),Form("%s_c1",part_name[i]),1100,700);
      c1[i]->Divide(2,2);
      TLegend *leg = new TLegend(0.2,0.7,0.4,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(run_type);
      for(int k=0; k<nLumi; k++)
	{
	  c1[i]->cd(k+1);
	  SetPadMargin(gPad, 0.13, 0.13);
	  ScaleHistoTitle(hplot,0.055,1,0.05,0.055,1,0.05,62);
	  if(i==0)
	    {
	      hplot->GetYaxis()->SetRangeUser(0.85,1.3);
	      hplot->GetXaxis()->SetRangeUser(0,15);
	    }
	  else
	    {
	      hplot->GetYaxis()->SetRangeUser(0.9,1.15);
	      hplot->GetXaxis()->SetRangeUser(0,10);
	    }
	  hplot->DrawCopy();
	  TPaveText *t1 = GetTitleText(Form("%s: uncertainty due to online trigger (%s)",part_title[i],name_lumi[k]),0.055);
	  t1->Draw();
	  for(int s=0; s<nSys; s++)
	    {
	      gEffSysVsPt[i][k][s]->SetMarkerStyle(20+s);
	      gEffSysVsPt[i][k][s]->SetMarkerColor(1+s);
	      gEffSysVsPt[i][k][s]->SetLineColor(1+s);
	      offset_x_with_asym_error(gEffSysVsPt[i][k][s], 0.2*s);
	      gEffSysVsPt[i][k][s]->Draw("samesPE");
	      if(k==0) leg->AddEntry(gEffSysVsPt[i][k][s],sys_legend[s],"P");
	    }

	  c1[i]->cd(k+3);
	  SetPadMargin(gPad, 0.13, 0.13);
	  ScaleHistoTitle(hplot2,0.055,1,0.05,0.055,1,0.05,62);
	  hplot2->GetYaxis()->SetRangeUser(0.9,1.1);
	  hplot2->GetXaxis()->SetLabelFont(62);
	  hplot2->GetXaxis()->SetLabelOffset(0.01);
	  hplot2->GetXaxis()->SetLabelSize(0.08);
	  hplot2->DrawCopy();
	  t1->Draw();
	  for(int s=0; s<nSys; s++)
	    {
	      gEffSysVsCent[i][k][s]->SetMarkerStyle(20+s);
	      gEffSysVsCent[i][k][s]->SetMarkerColor(1+s);
	      gEffSysVsCent[i][k][s]->SetLineColor(1+s);
	      offset_x_with_asym_error(gEffSysVsCent[i][k][s], 0.2*s);
	      gEffSysVsCent[i][k][s]->Draw("samesPE");
	    }
	}
      c1[i]->cd(1);
      leg->Draw();
      if(savePlot) c1[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdTrigEffSys_Comp.pdf",run_type,part_name[i]));
    }

  // Final efficiency
  // Method 2 as the default and symmetrize the error
  // Statistical weight: prod_mid/high (73.4%); prod_/low (26.6%)
  TCanvas *c = 0x0;
  const int sys_index = 0;
  const double weight[nLumi] = {0.734, 0.266};
  TH1F *hEffSysVsPtFinal[nPart];
  TH1F *hEffSysVsCentFinal[nPart];
  for(int i=0; i<nPart; i++)
    {
      hEffSysVsPtFinal[i] = (TH1F*)hEffVsPt[i][0][sys_index][0]->Clone(Form("%sEffVsPt_MtdTrigEff_FinalSys",part_name[i]));
      hEffSysVsPtFinal[i]->Reset();
      int nbins = hEffSysVsPtFinal[i]->GetNbinsX();
      for(int bin=1; bin<=nbins; bin++)
	{
	  double eh = gEffSysVsPt[i][0][sys_index]->GetErrorYhigh(bin-1) * weight[0] + gEffSysVsPt[i][1][sys_index]->GetErrorYhigh(bin-1) * weight[1];
	  double el = gEffSysVsPt[i][0][sys_index]->GetErrorYlow(bin-1)  * weight[0] + gEffSysVsPt[i][1][sys_index]->GetErrorYlow(bin-1)  * weight[1];
	  double err = (eh>el) ? eh : el;
	  hEffSysVsPtFinal[i]->SetBinContent(bin, 1);
	  hEffSysVsPtFinal[i]->SetBinError(bin, err);
	}
      hEffSysVsPtFinal[i]->GetYaxis()->SetRangeUser(0.85,1.15);
      c = draw1D(hEffSysVsPtFinal[i],Form("%s: uncertainty due to online trigger (%s);p_{T} (GeV/c)",part_title[i],run_type));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdTrigEff_FinalSysVsPt.pdf",run_type,part_name[i]));

      hEffSysVsCentFinal[i] = (TH1F*)hEffVsCent[i][0][sys_index][0]->Clone(Form("%sEffVsCent_MtdTrigEff_FinalSys",part_name[i]));
      hEffSysVsCentFinal[i]->Reset();
      hEffSysVsCentFinal[i]->GetXaxis()->SetBinLabel(1, "p_{T} > 0 GeV/c");
      hEffSysVsCentFinal[i]->GetXaxis()->SetBinLabel(2, "p_{T} > 5 GeV/c");
      hEffSysVsCentFinal[i]->GetXaxis()->SetLabelOffset(0.01);
      hEffSysVsCentFinal[i]->GetXaxis()->SetLabelSize(0.05);
      nbins = hEffSysVsCentFinal[i]->GetNbinsX();
      for(int bin=1; bin<=nbins; bin++)
	{
	  double eh = gEffSysVsCent[i][0][sys_index]->GetErrorYhigh(bin-1) * weight[0] + gEffSysVsCent[i][1][sys_index]->GetErrorYhigh(bin-1) * weight[1];
	  double el = gEffSysVsCent[i][0][sys_index]->GetErrorYlow(bin-1)  * weight[0] + gEffSysVsCent[i][1][sys_index]->GetErrorYlow(bin-1)  * weight[1];
	  double err = (eh>el) ? eh : el;
	  hEffSysVsCentFinal[i]->SetBinContent(bin, 1);
	  hEffSysVsCentFinal[i]->SetBinError(bin, err);
	}
      hEffSysVsCentFinal[i]->GetYaxis()->SetRangeUser(0.85,1.15);
      c = draw1D(hEffSysVsCentFinal[i],Form("%s: uncertainty due to online trigger (%s);;",part_title[i],run_type));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdTrigEff_FinalSysVsCent.pdf",run_type,part_name[i]));
    }

  if(saveHisto)
    {
      for(int i=0; i<nPart; i++)
	{
	  hEffSysVsPtFinal[i]->Write("",TObject::kOverwrite);
	  hEffSysVsCentFinal[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void anaSys(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  // single muon efficiency
  TFile *fMuonEff = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type),"read");
  hMuonPtEff[0] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff75_BinCount_FitFunc");
  hMuonPtEff[1] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff75_BinCount_FitFunc_Sysup");
  hMuonPtEff[2] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff75_BinCount_FitFunc_Sysdown");

  // save histogram
  TFile *fout = 0x0;
  if(saveHisto)
    {
      TString fName = fMuonEff->GetName();
      fName.ReplaceAll(".root",".sys.root");
      fout = TFile::Open(fName.Data(),"update");
    }

  for(int i=0; i<4; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7);
      
      if(saveHisto)
	{
	  fout->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void makeHisto(TString name, const double mass, const int nExpr)
{
  int nbins_tmp = 0;
  double xbins_tmp[10] = {0};
  if(name.Contains("Jpsi"))
    {
      TFile *fin = 0x0;
      if(year==2014)
	{
	  fin = TFile::Open("Rootfiles/models.root","read");
	  hMcJpsiPt = (TH1F*)fin->Get(Form("TBW_JpsiYield_AuAu200_cent0060"));
	}
      nbins_tmp = 9;
      double xbins_tmp_tmp[10] = {0,1,2,3,4,5,6,8,10,15};
      std::copy(std::begin(xbins_tmp_tmp), std::end(xbins_tmp_tmp), std::begin(xbins_tmp));
    }
  else
    {
      TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
      fBol->SetParameter(0,1.11);
      fBol->SetNpx(1000);
      hMcJpsiPt = (TH1F*)fBol->GetHistogram();
      nbins_tmp = 3;
      double xbins_tmp_tmp[4] = {0,2,4,10};
      std::copy(std::begin(xbins_tmp_tmp), std::end(xbins_tmp_tmp), std::begin(xbins_tmp));
    }
  hMcJpsiPt->Scale(1./hMcJpsiPt->Integral());
  
  // book histograms
  const int nbinsPt = nbins_tmp;
  double xbinsPt[nbinsPt+1];
  for(int i=0; i<nbinsPt+1; i++)
    {
      xbinsPt[i] = xbins_tmp[i];
    }
  TString hName;
  hInJpsiPt = new TH1F(Form("hInPartPt_%s",name.Data()), "", nbinsPt, xbinsPt);
  hInJpsiPt->Sumw2();
  hInJpsiCent = new TH1F(Form("hInPartCent_%s",name.Data()), "", 2, 0, 2);
  hInJpsiCent->Sumw2();
  for(int i=0; i<3; i++)
    {
      hName =  hMuonPtEff[i]->GetName();
      hName.ReplaceAll("DataJpsiMuon",Form("%sEffVsPt",name.Data()));
      hOutJpsiPt[i] = new TH1F(hName.Data(), "", nbinsPt, xbinsPt);
      hOutJpsiPt[i]->Sumw2();

      hName.ReplaceAll("VsPt","VsCent");
      hOutJpsiCent[i] = new TH1F(hName.Data(), "", 2, 0, 2);
      hOutJpsiCent[i]->Sumw2();
    }
  toyMC(mass, nExpr);

  TCanvas *c = 0x0;
  for(int e=0; e<3; e++)
    {
      hOutJpsiPt[e]->Divide(hInJpsiPt);
      hOutJpsiPt[e]->SetMarkerStyle(21);
      hOutJpsiPt[e]->SetMarkerColor(TMath::Power(2, e));
      hOutJpsiPt[e]->SetLineColor(TMath::Power(2, e));
      hOutJpsiPt[e]->GetYaxis()->SetRangeUser(0.4,1);
    }
  c = draw1D(hOutJpsiPt[0],"");
  hOutJpsiPt[1]->Draw("sames");
  hOutJpsiPt[2]->Draw("sames");

  for(int e=0; e<3; e++)
    {
      hOutJpsiCent[e]->Divide(hInJpsiCent);
      hOutJpsiCent[e]->SetMarkerStyle(21);
      hOutJpsiCent[e]->SetMarkerColor(TMath::Power(2, e));
      hOutJpsiCent[e]->SetLineColor(TMath::Power(2, e));
      hOutJpsiCent[e]->GetYaxis()->SetRangeUser(0.4,1);
    }
  c = draw1D(hOutJpsiCent[0],"");
  hOutJpsiCent[1]->Draw("sames");
  hOutJpsiCent[2]->Draw("sames");
}

//================================================
void toyMC(const double mass, const int nExpr, const int debug)
{
  double pt1_cut = 1.5;
  double pt2_cut = 1.3;
  if(mass>8) //Upsilon
    {
      pt1_cut = 4.0;
      pt2_cut = 1.5;
    }
  int nHisto = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("+++ Event %d +++\n",i+1);
      double mc_pt  = myRandom->Uniform(0,20);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      double weight = hMcJpsiPt->GetBinContent(hMcJpsiPt->FindFixBin(mc_pt));

      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      
      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);

      // acceptance cut
      if(fabs(eta1)>0.5 || fabs(eta2)>0.5) continue;
      if(fabs(pt1)<0.5  || fabs(pt2)<0.5)  continue;

      // momentum resolution 
      int mom_index1 = hTpcTrackRes->GetXaxis()->FindBin(pt1)-1;
      if(mom_index1>=nHisto) mom_index1=nHisto-1;
      int mom_index2 = hTpcTrackRes->GetXaxis()->FindBin(pt2)-1;
      if(mom_index2>=nHisto) mom_index2=nHisto-1;
      double dpt1 = hTrkResBin[mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[mom_index2]->GetRandom();
      double rc_pt1 = (1-dpt1) * pt1;
      double rc_pt2 = (1-dpt2) * pt2;
      double leadpt = rc_pt1 > rc_pt2 ? rc_pt1 : rc_pt2;
      double subpt  = rc_pt1 < rc_pt2 ? rc_pt1 : rc_pt2;
      if(leadpt<pt1_cut || subpt<pt2_cut) continue;

      hInJpsiPt->Fill(mc_pt,weight);
      if(mc_pt>0) hInJpsiCent->Fill(0.5,weight);
      if(mc_pt>5) hInJpsiCent->Fill(1.5,weight);

      // single muon efficiency
      for(int e=0; e<3; e++)
	{
	  double eff1 = hMuonPtEff[e]->Eval(rc_pt1);
	  double eff2 = hMuonPtEff[e]->Eval(rc_pt2);
	  if(myRandom->Uniform(0., 1.)<eff1 && myRandom->Uniform(0,1)<eff2)
	    {
	      hOutJpsiPt[e]->Fill(mc_pt,weight);
	      if(mc_pt>0) hOutJpsiCent[e]->Fill(0.5,weight);
	      if(mc_pt>5) hOutJpsiCent[e]->Fill(1.5,weight);
	    }
	}
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
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass) 
{
  double e = parent.M()/2.;
  double p = sqrt(e*e-dmass*dmass);
  double costheta = myRandom->Uniform(-1.0,1.0);
  double phi = myRandom->Uniform(0,TMath::Pi()*2);
  double pz = p*costheta;
  double px = p*sqrt(1.-costheta*costheta)*cos(phi);
  double py = p*sqrt(1.-costheta*costheta)*sin(phi);
  TLorentzVector daughter(px,py,pz,e);
  return myBoost(parent,daughter);
}


//-----------------------------------------
void ScaleHistoTitle(const TH1 *h, 
                     const Double_t xTitleSize, const Double_t xTitleOffset, const Double_t xLabelSize,
                     const Double_t yTitleSize, const Double_t yTitleOffset, const Double_t yLabelSize,
                     const Int_t font)
{
  if(!h) return;
  h->GetXaxis()->SetTitleFont(font);
  h->GetXaxis()->SetLabelFont(font);
  h->GetYaxis()->SetTitleFont(font);
  h->GetYaxis()->SetLabelFont(font);

  h->GetXaxis()->SetTitleSize(xTitleSize);
  h->GetXaxis()->SetTitleOffset(xTitleOffset);
  h->GetXaxis()->SetLabelSize(xLabelSize);
  h->GetYaxis()->SetTitleSize(yTitleSize);
  h->GetYaxis()->SetTitleOffset(yTitleOffset);
  h->GetYaxis()->SetLabelSize(yLabelSize);
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

//-----------------------------------------
TCanvas *draw2D(TH2 *h, const TString hTitle, const Float_t size, const Bool_t logz, const char *drawOption)
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
TCanvas *drawGraph(TGraph *h, const TString hTitle,Bool_t setLog, const Float_t size, const char *drawOption, const Int_t font)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),800,600);
  if(setLog)
    gPad->SetLogy();
  if(hTitle.Length()>0)
    h->SetTitle(hTitle);
  TPaveText *t1 = GetTitleText(h->GetHistogram()->GetTitle(),size,font);
  h->SetTitle("");
  h->Draw(drawOption);
  t1->Draw("sames");
  return c;
}

//--------------------------------------------
void offset_x_with_asym_error(TGraphAsymmErrors* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  Double_t* exl = g->GetEXlow();
  Double_t* exh = g->GetEXhigh();
  Double_t* eyl = g->GetEYlow();
  Double_t* eyh = g->GetEYhigh();

  for (Int_t j=0; j<npoints; j++)
    {
      g->SetPoint(j, x[j] + xoff, y[j]);
      g->SetPointError(j, exl[j], exh[j], eyl[j], eyh[j]);
    }
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

//--------------------------------------------
TLine *GetLine(double xl, double yl, double xh, double yh, Color_t color, Width_t width, Style_t style)
{
  TLine *line = new TLine(xl,yl,xh,yh);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  return line;
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

