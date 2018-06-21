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

#define YEAR 2014
#if (YEAR==2014)
const char *run_type = "Run14_AuAu200";
#elif (YEAR==2013)
char *run_type = "Run13_pp500";
#elif (YEAR==2015)
char *run_type = "Run15_pAu200";
#elif (YEAR==2016)
char *run_type = "Run16_AuAu200";
#endif

const double pi = 3.1415926;
const double muMass = 0.1057;
const int year = YEAR;
const int nPart = 4;
const char *part_name[nPart]  = {"Jpsi","Ups1S","Ups2S","Ups3S"};
const char *part_title[nPart] = {"Jpsi","Y(1S)","Y(2S)","Y(3S)"};
const double part_mass[nPart] = {3.097, 9.46, 10.023, 10.355};

TRandom3 *myRandom;
void anaSys(const char* type, const int saveHisto = 0);
void plotSys(const char* type, const int savePlot = 0, const int saveHisto = 0, const int saveAN = 0);
void mtdTrigEff(const int saveHisto = 0);
void mtdTrigEffLS(const int saveHisto = 0);
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

  const char* type = "MtdTrigEff";
  //const char* type = "Dtof";
  //const char* type = "MtdMthEff";
  //anaSys(type, 1);
  //plotSys(type, 1, 1, 0);
  //mtdTrigEff(0);
  mtdTrigEffLS(1);
}

//================================================
void plotSys(const char* type, const int savePlot, const int saveHisto, const int saveAN)
{
  char name_lumi[256], title_lumi[256], name_eff[256], title_eff[256];
  if(strcmp(type,"MtdTrigEff")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      sprintf(name_eff, "TacDiffEff");
      sprintf(title_eff, "MTD trigger efficiency");
    }

  else if(strcmp(type,"Dtof")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      if(year==2014) sprintf(name_eff, "Dtof0.75Eff");
      if(year==2015) sprintf(name_eff, "Dtof1.00Eff");
      sprintf(title_eff, "#Deltatof cut efficiency");
    }
  else if(strcmp(type,"MtdMthEff")==0)
    {
      sprintf(name_lumi, "");
      sprintf(title_lumi, "");
      sprintf(name_eff, "MtdMthEff");
      sprintf(title_eff, "MTD matching efficiency");
    }

  // const char* name_lumi = "";
  // const char* title_lumi = "";
  // const char* name_eff  = "DtofEff46";
  // const char* title_eff = "#Deltatof cut efficiency";

  const char *sys_name[3] = {"", "_Sysup", "_Sysdown"};
  TH1F *hEffVsPt[nPart][3];
  TH1F *hEffVsCent[nPart][3];

  TH1F *hEffSysVsPt[nPart];
  TH1F *hEffSysVsCent[nPart];  

  TFile *fin = 0x0;
  if(strcmp(type,"MtdTrigEff")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"read");
    }
  if(strcmp(type,"Dtof")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"read");
    }
  if(strcmp(type,"MtdMthEff")==0)
    {
      if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type),"update");
      else          fin = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type),"read");
    }

  for(int i=0; i<nPart; i++)
    {
      for(int s=0; s<3; s++)
	{
	  if(strcmp(type,"MtdTrigEff")==0)
	    {
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("%s_%sEffVsPt_%s%s",run_type,part_name[i],name_eff,sys_name[s]));
	    }
	  if(strcmp(type,"Dtof")==0)
	    {
	      char *tmp_name = Form("%s",sys_name[s]);
	      if(s==0) 
		{
		  if(year==2014) tmp_name = Form("_FitFunc");
		  if(year==2015) tmp_name = Form("_Syscenter");
		}
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("TagAndProbe_%sEffVsPt_%s%s",part_name[i],name_eff,tmp_name));
	    }
	  if(strcmp(type,"MtdMthEff")==0)
	    {
	      cout << Form("%sEffVsPt_%s%s",part_name[i],name_eff,sys_name[s]) << endl;
	      hEffVsPt[i][s] = (TH1F*)fin->Get(Form("%sEffVsPt_%s%s",part_name[i],name_eff,sys_name[s]));
	    }
	}
      int nBins = hEffVsPt[i][0]->GetNbinsX();
      hEffSysVsPt[i] = (TH1F*)hEffVsPt[i][0]->Clone(Form("%s_%sEffVsPt_Sys_%s",run_type,part_name[i],name_eff));
      hEffSysVsPt[i]->Reset();
      for(int bin=1; bin<=nBins; bin++)
	{
	  double diff1 = fabs(hEffVsPt[i][1]->GetBinContent(bin)/hEffVsPt[i][0]->GetBinContent(bin) - 1);
	  double diff2 = fabs(hEffVsPt[i][2]->GetBinContent(bin)/hEffVsPt[i][0]->GetBinContent(bin) - 1);
	  double error = diff1 > diff2 ? diff1 : diff2;
	  hEffSysVsPt[i]->SetBinContent(bin, 1);
	  hEffSysVsPt[i]->SetBinError(bin, error);
	}

      for(int s=0; s<3; s++)
	{
	  if(strcmp(type,"MtdTrigEff")==0)
	    {
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("%s_%sEffVsCent_%s%s",run_type,part_name[i],name_eff,sys_name[s]));
	    }
	  if(strcmp(type,"Dtof")==0)
	    {
	      char *tmp_name = Form("%s",sys_name[s]);
	      if(s==0) 
		{
		  if(year==2014) tmp_name = Form("_FitFunc");
		  if(year==2015) tmp_name = Form("_Syscenter");
		}
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("TagAndProbe_%sEffVsCent_%s%s",part_name[i],name_eff,tmp_name));
	    }
	  if(strcmp(type,"MtdMthEff")==0)
	    {
	      hEffVsCent[i][s] = (TH1F*)fin->Get(Form("%sEffVsCent_%s%s",part_name[i],name_eff,sys_name[s]));
	    }
	}
      nBins = hEffVsCent[i][0]->GetNbinsX();
      hEffSysVsCent[i] = (TH1F*)hEffVsCent[i][0]->Clone(Form("%s_%sEffVsCent_Sys_%s",run_type,part_name[i],name_eff));
      hEffSysVsCent[i]->Reset();
      for(int bin=1; bin<=nBins; bin++)
	{
	  double diff1 = fabs(hEffVsCent[i][1]->GetBinContent(bin)/hEffVsCent[i][0]->GetBinContent(bin) - 1);
	  double diff2 = fabs(hEffVsCent[i][2]->GetBinContent(bin)/hEffVsCent[i][0]->GetBinContent(bin) - 1);
	  double error = diff1 > diff2 ? diff1 : diff2;
	  hEffSysVsCent[i]->SetBinContent(bin, 1);
	  hEffSysVsCent[i]->SetBinError(bin, error);
	}
    }

  TCanvas *c1[nPart]; 
  TCanvas *c2[nPart]; 
  char outName[256];
  for(int i=0; i<nPart; i++)
    {
      c1[i]= new TCanvas(Form("%s_%s_SysVsPt",part_name[i],name_eff),Form("%s_%s_SysVsPt",part_name[i],name_eff),800,600);
      SetPadMargin(gPad, 0.13, 0.13);
      ScaleHistoTitle(hEffSysVsPt[i],0.05,1,0.04,0.05,1.2,0.04,42);
      hEffSysVsPt[i]->SetTitle(";p_{T} (GeV/c);Sys. Uncert.");
      if(year==2014) 
	{
	  if(strcmp(type,"MtdMthEff")==0)
	    hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.8,1.2);
	  else 
	    hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.9,1.1);
	}
      if(year==2015) hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(year==2016) hEffSysVsPt[i]->GetYaxis()->SetRangeUser(0.9,1.1);
      hEffSysVsPt[i]->SetMarkerStyle(20);
      hEffSysVsPt[i]->Draw();
      TPaveText *t1 = GetTitleText(Form("%s: uncertainty for %s",part_title[i],title_eff),0.05,42);
      t1->Draw();
      if(savePlot) 
	{
	  sprintf(outName, "~/Work/STAR/analysis/Plots/%s/ana_%s/Sys.%s_%sVsPt.pdf",run_type,type,part_name[i],name_eff);
	  c1[i]->SaveAs(outName);
	}
      if(saveAN && i<2) 
	{
	  if(strcmp(type,"Dtof")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_%sDtofEffSys.pdf",part_name[i]);
	    }
	  if(strcmp(type,"MtdTrigEff")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_%sMtdTrigEffSys.pdf",part_name[i]);
	    }
	  if(strcmp(type,"MtdMthEff")==0) 
	    {
	      sprintf(outName, "~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_%sMtdMthEffSys.pdf",part_name[i]);
	    }
	  c1[i]->SaveAs(outName);
	}

      c2[i]= new TCanvas(Form("%s_%s_SysVsCent",part_name[i],name_eff),Form("%s_%s_SysVsCent",part_name[i],name_eff),800,600);
      SetPadMargin(gPad, 0.13, 0.13);
      ScaleHistoTitle(hEffSysVsCent[i],0.05,1,0.04,0.05,1.2,0.04,42);
      hEffSysVsCent[i]->GetXaxis()->SetBinLabel(1, "p_{T} > 0 GeV/c");
      hEffSysVsCent[i]->GetXaxis()->SetBinLabel(2, "p_{T} > 5 GeV/c");
      if(strcmp(type,"MtdMthEff")==0)
	hEffSysVsCent[i]->GetYaxis()->SetRangeUser(0.8,1.2);
      else
	hEffSysVsCent[i]->GetYaxis()->SetRangeUser(0.9,1.1);
      hEffSysVsCent[i]->GetXaxis()->SetLabelFont(42);
      hEffSysVsCent[i]->GetXaxis()->SetLabelOffset(0.01);
      hEffSysVsCent[i]->GetXaxis()->SetLabelSize(0.07);
      hEffSysVsCent[i]->SetMarkerStyle(20);
      hEffSysVsCent[i]->Draw();
      t1->Draw();
      if(savePlot) 
	{
	  sprintf(outName, "~/Work/STAR/analysis/Plots/%s/ana_%s/Sys.%s_%sVsCent.pdf",run_type,type,part_name[i],name_eff);
	  c2[i]->SaveAs(outName);
	}
    }

  if(saveHisto)
    {
      for(int i=0; i<nPart; i++)
	{
	  hEffSysVsPt[i]->Write("",TObject::kOverwrite);
	  hEffSysVsCent[i]->Write("",TObject::kOverwrite);
	}
    }

}

//================================================
void anaSys(const char* type, const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  // single muon efficiency
  // TFile *fMuonEff = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type),"read");
  // hMuonPtEff[0] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc");
  // hMuonPtEff[1] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc_Sysup");
  // hMuonPtEff[2] = (TF1*)fMuonEff->Get("DataJpsiMuon_DtofEff46_BinCount_FitFunc_Sysdown");

  TFile *fdata = 0x0;
  if(strcmp(type,"MtdTrigEff")==0)
    {
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"read");
      hMuonPtEff[0] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff",run_type));
      hMuonPtEff[1] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff_Sysup",run_type));
      hMuonPtEff[2] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff_Sysdown",run_type));
    }
  else if(strcmp(type,"Dtof")==0)
    {
      double dtofCut = 0.75;
      if(year==2015)
	{
	  dtofCut = 1.0;
	}
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"read");
      if(year==2014) hMuonPtEff[0] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_FitFunc",dtofCut));
      if(year==2015) hMuonPtEff[0] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Syscenter",dtofCut));
      hMuonPtEff[1] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Sysup",dtofCut));
      hMuonPtEff[2] = (TF1*)fdata->Get(Form("TagAndProbe_Muon_Dtof%2.2fEff_Sysdown",dtofCut));
    }
  else if(strcmp(type,"MtdMthEff")==0)
    {
      if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type),"update");
      else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type),"read");
      hMuonPtEff[0] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type0");
      hMuonPtEff[0]->SetName("Muon_MtdMthEff");
      hMuonPtEff[1] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type1");
      hMuonPtEff[1]->SetName("Muon_MtdMthEff_Sysup");
      hMuonPtEff[2] = (TF1*)fdata->Get("funcTrkPtBtmBlEff_Type1");
      hMuonPtEff[2]->SetName("Muon_MtdMthEff_Sysdown");
    }

  for(int i=0; i<4; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e8);
      
      if(saveHisto)
	{
	  fdata->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void mtdTrigEff(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"read");
  hMuonPtEff[0] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff",run_type));

  TFile *fEff = TFile::Open(Form("Rootfiles/%s.MtdTrigEff.root",run_type),"read");
  hMuonPtEff[1] = (TF1*)fEff->Get(Form("%s_gTacDiffEffFinalFit_BinCount_prod_mid_Run15_pp200",run_type));
  hMuonPtEff[2] = (TF1*)fEff->Get(Form("%s_gTacDiffEffFinalFit_BinCount_prod_high_Run15_pp200",run_type));
  hMuonPtEff[1]->SetName(Form("%s_Muon_TacDiffEff_prod_mid",run_type));
  hMuonPtEff[2]->SetName(Form("%s_Muon_TacDiffEff_prod_high",run_type));

  for(int e=0; e<3; e++)
    {
      hMuonPtEff[e]->SetName(Form("%s_TrigStudy",hMuonPtEff[e]->GetName()));
      cout <<  hMuonPtEff[e]->GetName() << endl;
    }

  for(int i=0; i<1; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7);
      
      if(saveHisto)
	{
	  fdata->cd();
	  for(int e=0; e<3; e++)
	    {
	      hOutJpsiPt[e]->Write("",TObject::kOverwrite);
	      hOutJpsiCent[e]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void mtdTrigEffLS(const int saveHisto)
{
  // track momentum resolution
  TFile *fRes = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root"),"read");
  hTpcTrackRes = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
  int nHistos = hTpcTrackRes->GetNbinsX();
  for(int i=0; i<nHistos; i++)
    {
      hTrkResBin[i] = (TH1F*)hTpcTrackRes->ProjectionY(Form("hTrkRes_Bin%d",i+1),i+1,i+1);
    }

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdTrigEff.root",run_type),"read");
  hMuonPtEff[0] = (TF1*)fdata->Get(Form("%s_Muon_TacDiffEff",run_type));

  TFile *fEff = 0x0;
  if(saveHisto) fEff = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type),"update");
  else          fEff = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type),"read");

  hMuonPtEff[1] = (TF1*)fEff->Get(Form("%s_gTacDiffEff_LS_prod_mid_func",run_type));
  hMuonPtEff[2] = (TF1*)fEff->Get(Form("%s_gTacDiffEff_LS_prod_high_func",run_type));
  hMuonPtEff[1]->SetName(Form("%s_Muon_TacDiffEff_prod_mid",run_type));
  hMuonPtEff[2]->SetName(Form("%s_Muon_TacDiffEff_prod_high",run_type));

  for(int e=0; e<3; e++)
    {
      hMuonPtEff[e]->SetName(Form("%s_TrigStudy",hMuonPtEff[e]->GetName()));
    }

  for(int i=0; i<1; i++)
    {
      makeHisto(part_name[i],part_mass[i],1e7);
      
      if(saveHisto)
	{
	  fEff->cd();
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
      if(year==2014 || year==2016)
	{
	  fin = TFile::Open("Rootfiles/models.root","read");
	  hMcJpsiPt = (TH1F*)fin->Get(Form("TBW_JpsiYield_AuAu200_cent0060"));
	}
      else if(year==2015)
        {
          fin = TFile::Open("Rootfiles/JpsiSpectraShapepp200.root","read");
          TF1 *funcJpsi = (TF1*)fin->Get("TsallisPowerLawFitJpsipp200");
          funcJpsi->SetNpx(1000);
          hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
          hMcJpsiPt->SetName(Form("Jpsi_Yield_cent00100"));
          for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
            {
              hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
            }
        }

      if(year==2014)
	{
	  nbins_tmp = 9;
	  double xbins_tmp_tmp[10] = {0.15,1,2,3,4,5,6,8,10,15};
	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
	}
      else if(year==2015)
        {
          nbins_tmp = 8;
          double xbins_tmp_tmp[9] = {0,1,2,3,4,5,6,8,10};
    	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
        }
      else if(year==2016)
	{
	  nbins_tmp = 6;
	  double xbins_tmp_tmp[7] = {0,1,2,3,4,6,10};
	  for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
	}
    }
  else
    {
      TF1 *fBol = new TF1("Boltzmann","x/(exp(x/[0]+1))",0,10);
      fBol->SetParameter(0,1.11);
      fBol->SetNpx(1000);
      hMcJpsiPt = (TH1F*)fBol->GetHistogram();
      nbins_tmp = 3;
      double xbins_tmp_tmp[4] = {0,2,4,10};
      for(int i=0; i<nbins_tmp+1; i++) xbins_tmp[i] = xbins_tmp_tmp[i];
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
      hName.ReplaceAll("Muon",Form("%sEffVsPt",name.Data()));
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

