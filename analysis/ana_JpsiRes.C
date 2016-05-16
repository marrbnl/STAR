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

#define YEAR 2013

#if (YEAR==2014)
// +++ Run14 analysis +++
const char *run_type = "Run14_AuAu200";
const double pt1_cut = 1.0, pt2_cut = 1.0;
const double low_mass = 2.9;
const double high_mass = 3.3;
const int nCentBins = 4;
const int centBins_low[nCentBins]  = {5,13,9,5};
const int centBins_high[nCentBins] = {16,16,12,8};
const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60"};
const char *cent_Title[nCentBins] = {"0060","0020","2040","4060"};
const int nPtBins = 7;
const double ptBins_low[nPtBins]  = {0,1,2,3,4,6,8};
const double ptBins_high[nPtBins] = {10,2,3,4,6,8,10};
const char *pt_Name[nPtBins] = {"0-10","1-2","2-3","3-4","4-6","6-8","8-10"};
#elif (YEAR==2013)
// +++ Run13 analysis +++
const char *run_config = "VtxCut.";
//const char *run_config = "";
char *run_type = "Run13_pp500";
const double pt1_cut = 1.5, pt2_cut = 1.0;
const double low_mass = 2.9, high_mass = 3.25;
const int nPtBins = 7;
const double ptBins_low[nPtBins]  = {0,0,1,2,3,4,6};
const double ptBins_high[nPtBins] = {8,1,2,3,4,6,8};
const char *pt_Name[nPtBins] = {"0-8","0-1","1-2","2-3","3-4","4-6","6-8"};
const int nCentBins = 1;
const int centBins_low[nCentBins]  = {0};
const int centBins_high[nCentBins] = {-1};
const char *cent_Name[nCentBins] = {"0-100"};
const char *cent_Title[nCentBins] = {"00100"};
#endif

const double pi = 3.1415926;
const char *charge_name[3] = {"","_pos","_neg"};
const int gHistos = 6;
const char *det_name[gHistos] = {"Mc","TPCreco","TPCsmear","MTDmth","MTDresp","MTDtrig"};
const double jpsiMass = 3.096;
const double muMass = 0.1057;
const int year = YEAR;
const int kNPtBins = 200;
const double kLowPtBound = 0;
const double kHighPtBound = 20;
const int kNMassBins = 1400;
const double kLowMassBound = 0;
const double kHighMassBound = 14;

TRandom3 *myRandom;
TH3F *hTpcEff[nCentBins][2];
TH2F *hMtdEff[nCentBins][2][2];
TH2F *hTrkResVsPt[nCentBins];
TH1F *hTrkResBin[nCentBins][200];
TF1 *funcTrkRes[nCentBins];

TF1 *funcTrigEff[nCentBins];
TF1 *funcTrigEffCorr[nCentBins];

TH1F *hMtdRespEffCosmic;
TH1F *hMtdRespEffEmbed;

void tuneResolution(const int icent, const bool savePlot);
void upsilon(const int icent, const bool savePlot);
void trigEff();
void smear(const double mass, const int icent, const int nExpr, const double shift, const double sigma, TH1F *hInputPt, const bool debug, TH2F *hJpsiMassVsPt[gHistos]);
void tuneSmear(const double masss, const int icent, const int nExpr, const double shift, const double sigma, TH1F *hInputPt, TH2F *hRcJpsiMassVsPt);

TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);
double combinedFit(double *x, double *par);

TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size = 0.04, const Int_t font = 42);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);
TCanvas *draw2D(TH2 *h, const TString hTitle = "" , const Float_t size = 0.04, const Bool_t logz = kTRUE, const char *drawOption = "colz");
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color=2, Width_t width=2, Style_t style=2);


TH1F *hJpsiLineShape;

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


  // tracking efficiency
  TFile *fEff = 0x0;
  if(year==2014) fEff = TFile::Open("Rootfiles/Run14.AuAu200.TrkEff.root","read");
  if(year==2013) fEff = TFile::Open("Rootfiles/Run13.pp500.TrkEff.root","read");
  for(int icent=0; icent<nCentBins; icent++)
    {
      for(int i=0; i<2; i++)
	{
	  hTpcEff[icent][i] = (TH3F*)fEff->Get(Form("McTrkPtEtaPhiTpc%s_%s_Eff",charge_name[i+1],cent_Title[icent]));
	  hMtdEff[icent][0][i] = (TH2F*)fEff->Get(Form("McTrkPtEtaMtdEff%s_3tray_%s",charge_name[i+1],cent_Title[icent]));
	  hMtdEff[icent][1][i] = (TH2F*)fEff->Get(Form("McTrkPtEtaMtdEff%s_5tray_%s",charge_name[i+1],cent_Title[icent]));
	}

      // track resolution
      hTrkResVsPt[icent] = (TH2F*)fEff->Get(Form("PrimTrkRes_vs_TruePt_%s",cent_Title[icent]));
      funcTrkRes[icent] = (TF1*)fEff->Get(Form("FuncPrimTrkRes_vs_TruePt_%s",cent_Title[icent]));
      int nHistos = hTrkResVsPt[icent]->GetNbinsX();
      for(int i=0; i<nHistos; i++)
	{
	  hTrkResBin[icent][i] = (TH1F*)hTrkResVsPt[icent]->ProjectionY(Form("hTrkRes_Bin%d_cent%s",i+1,cent_Title[icent]),i+1,i+1);
	}
    }

  // single muon trigger efficiency
  if(year==2014)
    {
      TFile *fTrig =  TFile::Open(Form("Rootfiles/Run14.AuAu200.MuonTrigEff.root"),"read");
      for(int icent=0; icent<nCentBins; icent++)
	{
	  funcTrigEff[icent] = (TF1*)fTrig->Get(Form("MuonTrigEff_cent%s",cent_Title[icent]));
	  funcTrigEffCorr[icent] = (TF1*)fTrig->Get(Form("MuonTrigEffCorr_cent%s",cent_Title[icent]));
	}
    }
  if(year==2013)
    {
      for(int icent=0; icent<nCentBins; icent++)
	{
	  funcTrigEff[icent] = 0x0;
	  funcTrigEffCorr[icent] = 0x0;
	}
    }

  // MTD response efficiency
  TFile *fCosmic =  0x0;
  if(year==2014) fCosmic = TFile::Open(Form("Rootfiles/Run14.AuAu200.MtdResponseEff.root"),"read");
  if(year==2013) fCosmic = TFile::Open(Form("Rootfiles/Run13.pp500.MtdResponseEff.root"),"read");

  hMtdRespEffCosmic = (TH1F*)fCosmic->Get("MtdResponseEff_cosmic");
  hMtdRespEffEmbed = (TH1F*)fCosmic->Get("MtdResponseEff_embed");

  //tuneResolution(0,1);
  //upsilon(0,1);
  //trigEff();
}

//================================================
void upsilon(const int icent, const bool savePlot)
{
  TString name[3] = {"TPCemb","TPCreco","MTDreco"};
  TFile *fin = TFile::Open("Rootfiles/Upsion.pp200.root","read");
  TH1F *hMcUpsilonPt = (TH1F*)fin->Get("Upsilon_pp200_midRap");
  double upsilonMass[3] = {9.46,10.023,10.355};

  TH2F *hRcMassVsPt[3][gHistos];
  for(int j=0; j<3; j++)
    {
      for(int i=0; i<gHistos; i++)
	{
	  hRcMassVsPt[j][i] = new TH2F(Form("hRcMassVsPt_Upsilon%dS_%s_%s",j+1,det_name[i],cent_Title[icent]),Form("Mass distribution of reconstructed #Upsilon%dS (%s%%);p_{T} (GeV/c)",j+1,cent_Name[icent]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
	  hRcMassVsPt[j][i]->Sumw2();
	}

      smear(upsilonMass[j], icent, 1e8, 0, 0.017, hMcUpsilonPt, 0, hRcMassVsPt[j]);
    }


  for(int j=0; j<3; j++)
    {
      for(int i=1; i<2; i++)
	{
	  TH1F *hRcPt = (TH1F*)hRcMassVsPt[j][i]->ProjectionX(Form("hRcPt_Upsilon%dS_%d",j+1,i));
	  hRcPt->SetLineColor(2);
	  TCanvas *c = draw1D(hRcPt,Form("%s: p_{T} distribution of #Upsilon%dS (%s%%)",det_name[i],j+1,cent_Name[icent]),kTRUE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_Pt_cent%s.pdf",run_type,det_name[i],j+1,cent_Title[icent]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_Pt_cent%s.png",run_type,det_name[i],j+1,cent_Title[icent]));
	    }

	  TH1F *hMass = (TH1F*)hRcMassVsPt[j][i]->ProjectionY(Form("%s_proj",hRcMassVsPt[j][i]->GetName()));
	  hMass->SetLineColor(2);
	  c = draw1D(hMass,Form("%s: mass distribution of #Upsilon%dS (%s%%)",det_name[i],j+1,cent_Name[icent]),kFALSE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_mass_cent%s.pdf",run_type,det_name[i],j+1,cent_Title[icent]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_mass_cent%s.png",run_type,det_name[i],j+1,cent_Title[icent]));
	    }
	}
    }

  if(1)
    {
      cout << "Save histogram" << endl;
      TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.Upsilon.LineShape.root","recreate");
      for(int j=0; j<3; j++)
	{
	  for(int i=0; i<gHistos; i++)
	    {
	      TH1F *hMass = (TH1F*)hRcMassVsPt[j][i]->ProjectionY(Form("hRcMass_Upsilon%dS_%s_%s",j+1,det_name[i],cent_Title[icent]));
	      hMass->Write();
	    }
	} 
    }
}

//================================================
void trigEff()
{
  TFile *fin = 0x0;
  if(year==2013) fin = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
  if(year==2014) fin = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");

  TF1 *funcMcJpsiPt[nCentBins];
  TH1F *hMcJpsiPt[nCentBins];
  TH2F *hRcJpsiMass[nCentBins][gHistos];

  TH1F *hJpsiEff[nCentBins][gHistos-1][2];
  const char *eff_name[gHistos-1] = {"Tpc","Smear","Match","Resp","Trig"};

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  double sigma = 0, shift = 0;
  if(year==2013) 
    {
      sigma = 0.021;
      shift = -0.01;
    }
  if(year==2014) sigma = 0.016;

  for(int icent=0; icent<nCentBins; icent++)
    {
      cout << "Cent " << icent << endl;
      if(year==2014)
	hMcJpsiPt[icent] = (TH1F*)fin->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[icent]));

      if(year==2013)
	{
	  TF1 *funcJpsi = (TF1*)fin->Get("ffpt");
	  funcJpsi->SetNpx(1000);
	  hMcJpsiPt[icent] = (TH1F*)funcJpsi->GetHistogram();
	  hMcJpsiPt[icent]->SetName(Form("GlobalFit_Jpsi_Yield_cent%s",cent_Title[icent]));
	  for(int bin=1; bin<=hMcJpsiPt[icent]->GetNbinsX(); bin++)
	    {
	      hMcJpsiPt[icent]->SetBinContent(bin,hMcJpsiPt[icent]->GetBinCenter(bin)*hMcJpsiPt[icent]->GetBinContent(bin));
	    }
	}

      for(int i=0; i<gHistos; i++)
	{
	  hRcJpsiMass[icent][i] = new TH2F(Form("hRcJpsiMass_%s_%s",det_name[i],cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
	  hRcJpsiMass[icent][i]->Sumw2();
	}
      smear(jpsiMass, 0, 4e8, shift, sigma, hMcJpsiPt[icent], 0, hRcJpsiMass[icent]);

      for(int i=0; i<gHistos-1; i++)
	{
	  hJpsiEff[icent][i][0] = (TH1F*)hRcJpsiMass[icent][i+1]->ProjectionX(Form("Jpsi%sEff_cent%s",eff_name[i],cent_Title[icent]));
	  hJpsiEff[icent][i][1] = (TH1F*)hJpsiEff[icent][i][0]->Rebin(nbins,Form("%s_rebin",hJpsiEff[icent][i][0]->GetName()),xbins);
	  TH1F *hAll = (TH1F*)hRcJpsiMass[icent][i]->ProjectionX(Form("JpsiPt_%s_cent%s",det_name[i],cent_Title[icent]));
	  TH1F *hAllrebin = (TH1F*)hAll->Rebin(nbins,Form("%s_rebin",hAll->GetName()),xbins);
	  hJpsiEff[icent][i][0]->Divide(hAll);
	  hJpsiEff[icent][i][1]->Divide(hAllrebin);
	  TCanvas *c = draw1D(hJpsiEff[icent][i][0]);
	  hJpsiEff[icent][i][1]->Draw("sames");
	}
    }
 
  TString outName = "";
  if(year==2014) outName = Form("Run14.AuAu200.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  if(year==2013) outName = Form("Run13.pp500.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  TFile *fEff = TFile::Open(Form("Rootfiles/%s",outName.Data()),"update");

  for(int icent=0; icent<nCentBins; icent++)
    {
      for(int i=0; i<gHistos-1; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hJpsiEff[icent][i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
  
}

//================================================
void tuneResolution(const int icent, const bool savePlot)
{
  const int anaType = 3; // 0 - scan; 1 - determine smear & shift; 2 - generate final smear
  const int nShiftScan = 10;
  const double shiftStep = -0.005;
  const int nSmearScan = 10;
  const double smearStep = 0.002;
  const int index = 5;

  // input spectrum shape
  TH1F *hMcJpsiPt;
  if(year==2013)
    {
      TFile *fWeight = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
      TF1 *funcJpsi = (TF1*)fWeight->Get("ffpt");
      funcJpsi->SetNpx(1000);
      hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
      hMcJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield_cent%s",cent_Title[0]));
      for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
	{
	  hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
	}
    }
  if(year==2014)
    {
      TFile *fWeight = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");
      hMcJpsiPt = (TH1F*)fWeight->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[icent]));
    } 

  // scan smear and shift
  TH2F *hRcJpsiMassScan[nShiftScan][nSmearScan][gHistos];
  if(anaType == 0)
    {
      for(int i=0; i<nShiftScan; i++)
	{
	  for(int k=0; k<nSmearScan; k++)
	    {
	      for(int j=0; j<gHistos; j++)
		{
		  hRcJpsiMassScan[i][k][j] = new TH2F(Form("hRcJpsiMassVsPt_%s_%s_scan%d_%d",det_name[j],cent_Title[icent],i,k),Form("Mass distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c);mass (GeV/c^{2})",cent_Name[icent]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
		}

	      if(i>0 && k>0) continue;
	      printf("[i] Scan (%d,%d)\n",i,k);
	      double sigma = k*smearStep;
	      double shift = i*shiftStep;

	      //smear(jpsiMass, icent, 1e7, shift, sigma, hMcJpsiPt, 0, hRcJpsiMassScan[i][k]);
	      tuneSmear(jpsiMass, icent, 1e6, shift, sigma, hMcJpsiPt, hRcJpsiMassScan[i][k][index]);
	    }
	}
      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","recreate");
      if(year==2013) fout = TFile::Open(Form("Rootfiles/Run13.pp500.TrkResScan.%sroot",run_config),"recreate");
      for(int i=0; i<nShiftScan; i++){
	for(int k=0; k<nSmearScan; k++){
	  for(int j=0; j<gHistos; j++){
	    if(i>0 && k>0) continue;
	    hRcJpsiMassScan[i][k][j]->Write();
	  }
	}
      }
    }

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // final smear and shift
  TH2F *hRcJpsiMass[3][gHistos];
  const char *sysName_scan[3] = {"def","min","max"};
  if(anaType == 2)
    {
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<gHistos; j++)
	    {
	      hRcJpsiMass[i][j] = new TH2F(Form("hRcJpsiMass_%s_%s_final_%s",det_name[j],cent_Title[icent],sysName_scan[i]),Form("Mass distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c);mass (GeV/c^{2})",cent_Name[icent]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
	      hRcJpsiMass[i][j]->Sumw2();
	    }
	}

      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","update");
      if(year==2013) fout = TFile::Open(Form("Rootfiles/Run13.pp500.TrkResScan.%sroot",run_config),"update");
      TH1F *hFinalShift = (TH1F*)fout->Get(Form("hFinalShift_cent%s",cent_Title[icent]));
      TH1F *hFinalSmear = (TH1F*)fout->Get(Form("hFinalSmear_cent%s",cent_Title[icent]));

      for(int i=0; i<3; i++)
	{
	  double shift = hFinalShift->GetBinContent(i+1);
	  double smearvalue = hFinalSmear->GetBinContent(i+1);
	  printf("[i] Start final scanning: %s = (%4.4f,%4.4f)\n",sysName_scan[i],shift,smearvalue);
	  smear(jpsiMass, icent, 1e8, shift, smearvalue, hMcJpsiPt, 0, hRcJpsiMass[i]);
	}

      TH1F *hJpsiEff[3][2];
      for(int i=0; i<3; i++)
	{
	  hJpsiEff[i][0] = (TH1F*)hRcJpsiMass[i][2]->ProjectionX(Form("JpsiSmearEff_%s_cent%s",sysName_scan[i],cent_Title[icent]));
	  hJpsiEff[i][1] = (TH1F*)hJpsiEff[i][0]->Rebin(nbins,Form("JpsiSmearEff_%s_cent%s_rebin",sysName_scan[i],cent_Title[icent]),xbins);

	  TH1F *hAll = (TH1F*)hRcJpsiMass[i][1]->ProjectionX(Form("JpsiPt_all_%s_cent%s",sysName_scan[i],cent_Title[icent]));
	  TH1F *hAllrebin = (TH1F*)hAll->Rebin(nbins,Form("%s_rebin",hAll->GetName()),xbins);
	  hJpsiEff[i][0]->Divide(hAll);
	  hJpsiEff[i][1]->Divide(hAllrebin);
	}

      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<gHistos; j++){
	    hRcJpsiMass[i][j]->Write("",TObject::kOverwrite);
	  }
	  hJpsiEff[i][0]->Write("",TObject::kOverwrite);
	  hJpsiEff[i][1]->Write("",TObject::kOverwrite);
	}
    }

  // check final smear and shift
  TFile *fscan = 0x0;
  TFile *fdata = 0x0;
  if(anaType == 3)
    {
      if(year==2014) fscan = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","update");
      if(year==2013) fscan = TFile::Open(Form("Rootfiles/Run13.pp500.TrkResScan.%sroot",run_config),"update");
      
      TH2F *hMassVsPtEmbed = (TH2F*)fscan->Get(Form("hRcJpsiMass_%s_%s_final_def",det_name[index],cent_Title[icent]));
      
      TObjArray embedSlice;
      hMassVsPtEmbed->RebinX(5);
      hMassVsPtEmbed->FitSlicesY(0, 0, -1, 0, "QNR", &embedSlice);
      TH1F *hEmbSmearMean  =  (TH1F*)((TH1F*)embedSlice[1])->Clone("hEmbSmearMean");
      TH1F *hEmbSmearSigma =  (TH1F*)((TH1F*)embedSlice[2])->Clone("hEmbSmearSigma");      
  
      if(year==2014) fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
      if(year==2013) fdata = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut),"read");
  
      // mean
      TH1F *hDataMean = (TH1F*)fdata->Get(Form("Jpsi_FitMean_cent%s",cent_Title[icent]));
      hDataMean->SetMarkerStyle(21);
      hDataMean->SetMarkerColor(2);
      hDataMean->SetLineColor(2);
      hDataMean->SetTitle("Mean of J/#Psi peak");
      hDataMean->GetYaxis()->SetRangeUser(3.04,3.16);
      TCanvas *c = draw1D(hDataMean);
      hEmbSmearMean->SetMarkerStyle(25);
      hEmbSmearMean->Draw("sames");
      TLegend *leg = new TLegend(0.18,0.68,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hDataMean,"Real data","PL");
      leg->AddEntry(hEmbSmearMean,"Smeared embedding","PL");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_SmearEmbedVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_SmearEmbedVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      // sigma
      TH1F *hDataSigma = (TH1F*)fdata->Get(Form("Jpsi_FitSigma_cent%s",cent_Title[icent]));
      hDataSigma->SetMarkerStyle(21);
      hDataSigma->SetMarkerColor(2);
      hDataSigma->SetLineColor(2);
      hDataSigma->SetTitle("Width of J/#Psi peak");
      c = draw1D(hDataSigma);
      hEmbSmearSigma->SetMarkerStyle(25);
      hEmbSmearSigma->Draw("sames");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_SmearEmbedVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_SmearEmbedVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      // pt distribution
      TH1F *hRawCounts = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]));
      hRawCounts->Scale(1./hRawCounts->Integral("width"));
      hRawCounts->SetMarkerStyle(21);
      hRawCounts->SetMarkerColor(2);
      hRawCounts->SetLineColor(2);
      hRawCounts->SetTitle("Raw p_{T} distribution of J/#Psi");
      hRawCounts->GetYaxis()->SetRangeUser(0,0.5);
      c = draw1D(hRawCounts);
      TH1F *hCountsEmb = (TH1F*)hMassVsPtEmbed->ProjectionX("hCountsEmb");
      hCountsEmb->Rebin(2);
      hCountsEmb->Scale(1./hCountsEmb->Integral("width"));
      hCountsEmb->SetMarkerStyle(25);
      hCountsEmb->Draw("sames");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_SmearEmbedVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_SmearEmbedVsData_cent%s.png",run_type,cent_Title[icent]));
	}
      

      // signal shape
      TH1F *hSignal = (TH1F*)fdata->Get(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));
      TF1 *fitSignal = (TF1*)fdata->Get(Form("FitJpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));
      TF1 *funcBkg = new TF1("funcBkg","pol3",2.5,4);
      for(int i=0; i<4; i++)
	{
	  funcBkg->SetParameter(i, fitSignal->GetParameter(i+3));
	}
      TAxis *axis = hSignal->GetXaxis();
      TH1F *hSignalSub = new TH1F("hSignalSub","Compare J/#psi signal shape;M_{#mu#mu} (GeV/c^{2})",axis->GetNbins(),axis->GetXmin(), axis->GetXmax());
      for(int ibin=1; ibin<=hSignalSub->GetNbinsX(); ibin++)
	{
	  double value = hSignal->GetBinContent(ibin);
	  double error = hSignal->GetBinError(ibin);
	  double bkg = funcBkg->Eval(hSignal->GetBinCenter(ibin));
	  hSignalSub->SetBinContent(ibin, value - bkg);
	  hSignalSub->SetBinError(ibin, sqrt(error*error+fabs(bkg)));
	}
      hSignalSub->SetMarkerStyle(21);
      hSignalSub->GetXaxis()->SetRangeUser(2.5,4);
      hSignalSub->SetMaximum(1.6*hSignalSub->GetMaximum());
      hSignalSub->Scale(1./hSignalSub->Integral("width"));
      hSignalSub->GetYaxis()->SetRangeUser(-2,15);
      hSignalSub->SetMarkerColor(2);
      hSignalSub->SetLineColor(2);
      c = draw1D(hSignalSub);
      TH1F *hSignalEmb = (TH1F*)hMassVsPtEmbed->ProjectionY("hSignalEmb");
      hSignalEmb->SetMarkerStyle(25);
      hSignalEmb->Scale(1./hSignalEmb->Integral("width"));
      hSignalEmb->Draw("sames");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiShape_SmearEmbedVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiShape_SmearEmbedVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      // smearing efficiency
      TH1F *hJpsiEff[3];
      TGraphErrors *gEffSys = new TGraphErrors();
      for(int i=0; i<3; i++)
	{
	  hJpsiEff[i] = (TH1F*)fscan->Get(Form("JpsiSmearEff_%s_cent%s_rebin",sysName_scan[i],cent_Title[icent]));
	}
      gEffSys->Set(hJpsiEff[0]->GetNbinsX());
      gEffSys->SetName(Form("JpsiSmearEff_def_cent%s_sys_rebin",cent_Title[icent]));
      for(int i=0; i<gEffSys->GetN(); i++)
	{
	  double x = hJpsiEff[0]->GetBinCenter(i+1);
	  double y = hJpsiEff[0]->GetBinContent(i+1);
	  double ex = hJpsiEff[0]->GetBinWidth(i+1)/2;
	  double ey = fabs(hJpsiEff[0]->GetBinContent(i+1)-hJpsiEff[1]->GetBinContent(i+1));
	  double check = fabs(hJpsiEff[0]->GetBinContent(i+1)-hJpsiEff[2]->GetBinContent(i+1));
	  if(ey<check) ey = check;
	  gEffSys->SetPoint(i,x,y);
	  gEffSys->SetPointError(i,ex,ey);
	}
      hJpsiEff[0]->GetYaxis()->SetRangeUser(0.96,1.04);
      hJpsiEff[0]->SetTitle("Correction factor for smearing;p_{T} (GeV/c)");
      c = draw1D(hJpsiEff[0]);
      gEffSys->SetFillStyle(1001);
      gEffSys->SetFillColor(kGray);
      gEffSys->Draw("sames e2");
      hJpsiEff[0]->SetMarkerStyle(21);
      hJpsiEff[0]->Draw("samesP");
      leg = new TLegend(0.48,0.65,0.62,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hJpsiEff[0],"Correction factor","PL");
      leg->AddEntry(gEffSys,"Uncertainty","F");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/SmearEff_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/SmearEff_cent%s.png",run_type,cent_Title[icent]));
	}
      fscan->cd();
      gEffSys->Write("",TObject::kOverwrite);
    }

  if(anaType == 1)
    {
      if(year==2014) fscan = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","update");
      if(year==2013) fscan = TFile::Open(Form("Rootfiles/Run13.pp500.TrkResScan.%sroot",run_config),"update");

      // check with embedding
      TFile *fembed = 0x0;
      if(year==2014) fembed = TFile::Open(Form("Rootfiles/Run14.AuAu200.JpsiEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut),"read");
      if(year==2013) fembed = TFile::Open(Form("Rootfiles/Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut),"read");

      TH2F *hEmbedJpsi[2];
      hEmbedJpsi[0] = (TH2F*)fembed->Get(Form("TPCreco_Jpsi_MassVsPt_%s",cent_Title[icent]));
      hEmbedJpsi[1] = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_%s_scan0_0",det_name[index],cent_Title[icent]));

      TObjArray embSlices[2];
      TH1F *hEmbMean[2];
      TH1F *hEmbSigma[2];
      for(int i=0; i<2; i++)
	{
	  hEmbedJpsi[i]->RebinX(5);
	  hEmbedJpsi[i]->FitSlicesY(0, 0, -1, 0, "QNR", &embSlices[i]);
	  hEmbMean[i] = (TH1F*)((TH1F*)embSlices[i][1])->Clone(Form("hEmbMean_%d",i));
	  hEmbMean[i]->SetMarkerStyle(21+i*4);
	  hEmbMean[i]->SetMarkerColor(i+1);
	  hEmbMean[i]->SetLineColor(i+1);
	  hEmbMean[i]->GetXaxis()->SetRangeUser(0,8);
	  hEmbMean[i]->GetYaxis()->SetRangeUser(3.085,3.115);
	  hEmbMean[i]->SetTitle(";p_{T} (GeV/c);Mean");

	  hEmbSigma[i] = (TH1F*)((TH1F*)embSlices[i][2])->Clone(Form("hEmbSigma_%d",i));
	  hEmbSigma[i]->SetMarkerStyle(21+i*4);
	  hEmbSigma[i]->SetMarkerColor(i+1);
	  hEmbSigma[i]->SetLineColor(i+1);
	  hEmbSigma[i]->GetXaxis()->SetRangeUser(0,8);
	  hEmbSigma[i]->GetYaxis()->SetRangeUser(0.02,0.06);
	  hEmbSigma[i]->SetTitle(";p_{T} (GeV/c);#sigma");
	}
      TCanvas *c1 = new TCanvas("embed_mean", "embed_mean", 800, 600);
      hEmbMean[0]->Draw();
      hEmbMean[1]->Draw("sames");
      TPaveText *t1 = GetTitleText("Mean of J/#Psi mass peak");
      t1->Draw();
      TLegend *leg = new TLegend(0.18,0.65,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hEmbMean[0],"Embedding","PL");
      leg->AddEntry(hEmbMean[1],"ToyMC","PL");
      leg->Draw();
      if(savePlot)
	{
	  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ToyMcVsEmbed_cent%s.pdf",run_type,cent_Title[icent]));
	  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ToyMcVsEmbed_cent%s.png",run_type,cent_Title[icent]));
	}

      c1 = new TCanvas("embed_sigma", "embed_sigma", 800, 600);
      hEmbSigma[0]->Draw();
      hEmbSigma[1]->Draw("sames");
      t1 = GetTitleText("Width of J/#Psi mass peak");
      t1->Draw();
      leg->Draw();
      if(savePlot)
	{
	  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ToyMcVsEmbed_cent%s.pdf",run_type,cent_Title[icent]));
	  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ToyMcVsEmbed_cent%s.png",run_type,cent_Title[icent]));
	}

      // check with real data
      gStyle->SetOptFit(0);
      const int makePlot = 0;
      const char *sysName[3] = {"def","min","max"};
      const char *sysTitle[3] = {"default values","1#sigma lower limit","1#sigma upper limit"};
      if(year==2014) fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
      if(year==2013) fdata = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut),"read");

      TH1F *hFinalShift = new TH1F(Form("hFinalShift_cent%s",cent_Title[icent]),"",3,0,3);
      TH1F *hFinalSmear = new TH1F(Form("hFinalSmear_cent%s",cent_Title[icent]),"",3,0,3);
      for(int i=0; i<3; i++)
	{
	  hFinalShift->GetXaxis()->SetBinLabel(i+1,sysName[i]);
	  hFinalSmear->GetXaxis()->SetBinLabel(i+1,sysName[i]);
	}
      
      // determine shift
      TCanvas *cShift = new TCanvas("scan_shift", "scan_shift", 800, 600);
      TObjArray shiftSlices[nShiftScan];
      TH1F *hShiftMean[nShiftScan];
      TF1 *hFitShiftMean[nShiftScan];
      for(int i=0; i<nShiftScan; i++)
	{
	  TH2F *htmp = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_%s_scan%d_%d",det_name[index],cent_Title[icent],i,0));
	  TH2F *h2 = (TH2F*)htmp->Clone(Form("%s_clone",htmp->GetName()));
	  if(i>0) h2->RebinX(5);
	  hShiftMean[i] = (TH1F*)h2->ProjectionX(Form("hShiftMean_%d",i));
	  hShiftMean[i]->Reset();
	  TCanvas *cShiftFit = 0x0;
	  if(makePlot)
	    {
	      cShiftFit = new TCanvas(Form("fit_shift_%d",i), Form("fit_shift_%d",i), 1200, 800);
	      cShiftFit->Divide(5,4);
	    }
	  for(int j=0; j<20; j++)
	    {
	      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hShiftMass_%d_%d",i,j),j+1,j+1);
	      h1->GetXaxis()->SetRangeUser(2.6,3.8);
	      TF1 *funcTmp = new TF1(Form("funcTemp_shift_%d_%d",i,j),"gaus",2.7,3.5);
	      funcTmp->SetLineColor(2);
	      h1->Fit(funcTmp,"IRQ0");
	      hShiftMean[i]->SetBinContent(j+1,funcTmp->GetParameter(1));
	      hShiftMean[i]->SetBinError(j+1,funcTmp->GetParError(1));
	      if(makePlot)
		{
		  cShiftFit->cd(j+1);
		  h1->Draw();
		  funcTmp->Draw("sames");
		}
	    }
	  hShiftMean[i]->GetXaxis()->SetRangeUser(0,8);
	  hShiftMean[i]->GetYaxis()->SetRangeUser(2.9,3.15);
	  hShiftMean[i]->SetMarkerStyle(20+i);
	  hShiftMean[i]->SetTitle(";p_{T} (GeV/c);Mean");
	  cShift->cd();
	  if(i==0) hShiftMean[i]->Draw();
	  else     hShiftMean[i]->Draw("sames");
	  hFitShiftMean[i] = new TF1(Form("hFitShiftMean_%d",i),"pol4",0,8);
	  hShiftMean[i]->Fit(hFitShiftMean[i],"IRQ0");
	  hFitShiftMean[i]->SetLineColor(4);
	  hFitShiftMean[i]->SetLineStyle(2);
	  hFitShiftMean[i]->Draw("sames");
	}
      TH1F *hDataMean = (TH1F*)fdata->Get(Form("Jpsi_FitMean_cent%s",cent_Title[icent]));
      hDataMean->SetMarkerStyle(21);
      hDataMean->SetMarkerColor(2);
      hDataMean->SetLineColor(2);
      hDataMean->Draw("sames");
      t1 = GetTitleText("Mean of J/#Psi mass peak");
      t1->Draw();
      t1 = GetPaveText(0.45,0.65,0.3,0.35);
      t1->AddText("Real data");
      t1->SetTextAlign(11);
      t1->SetTextColor(2);
      t1->Draw();
      t1 = GetPaveText(0.45,0.65,0.25,0.3);
      t1->AddText("Scan shifting");
      t1->SetTextAlign(11);
      t1->Draw();
      t1 = GetPaveText(0.45,0.65,0.2,0.25);
      t1->AddText("Fit to scan data");
      t1->SetTextAlign(11);
      t1->SetTextColor(4);
      t1->Draw();
      if(savePlot)
	{
	  cShift->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ScanToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  cShift->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ScanToyMcVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      
      printf("[i] Scan shift\n");
      TH1F *hShiftChi2[3];
      TF1 *hFitShiftChi2[3];
      for(int i=0; i<3; i++)
	{
	  hShiftChi2[i] = new TH1F(Form("hShiftChi2_%s",sysName[i]),";shift (%); #chi^{2}", nShiftScan, (nShiftScan-1)*shiftStep+shiftStep/2, 0-shiftStep/2);
	  for(int j=0; j<nShiftScan; j++)
	    {
	      int bin = nShiftScan - j;
	      double chi2 = 0;
	      for(int ibin=1; ibin<=hDataMean->GetNbinsX(); ibin++)
		{
		  double value = hDataMean->GetBinContent(ibin);
		  double error = hDataMean->GetBinError(ibin);
		  if(i==1) value -= error;
		  if(i==2) value += error;
		  chi2 += TMath::Power((value-hFitShiftMean[j]->Eval(hDataMean->GetBinCenter(ibin)))/error,2);
		}
	      hShiftChi2[i]->SetBinContent(bin,chi2);
	      hShiftChi2[i]->SetBinError(bin,1e-10);
	    }
	  c1 = new TCanvas(Form("shift_chi2_%s",sysName[i]), Form("shift_chi2_%s",sysName[i]), 800, 600);
	  hShiftChi2[i]->SetMarkerStyle(20);
	  hShiftChi2[i]->Draw("P");
	  TAxis *axis = hShiftChi2[i]->GetXaxis();
	  hFitShiftChi2[i] = new TF1(Form("hFitShiftChi2_%d",i),"pol4",axis->GetXmin(),axis->GetXmax());
	  hShiftChi2[i]->Fit(hFitShiftChi2[i],"IRQ0");
	  hFitShiftChi2[i]->SetLineColor(4);
	  hFitShiftChi2[i]->SetLineStyle(2);
	  hFitShiftChi2[i]->Draw("sames");
	  t1 = GetTitleText(Form("Scan shift: use %s of data points",sysTitle[i]));
	  t1->Draw();
	  t1 = GetPaveText(0.6,0.7,0.7,0.8);
	  t1->AddText(Form("Minimum #chi^{2} at %4.2f%%",hFitShiftChi2[i]->GetMinimumX()*100));
	  t1->Draw();
	  hFinalShift->SetBinContent(i+1,hFitShiftChi2[i]->GetMinimumX());
	  if(savePlot)
	    {
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_Chi2_%s_cent%s.pdf",run_type,sysName[i],cent_Title[icent]));
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_Chi2_%s_cent%s.png",run_type,sysName[i],cent_Title[icent]));
	    }
	}


      // determine smear
      TCanvas *cSmear = new TCanvas("scan_smear", "scan_smear", 800, 600);
      TObjArray smearSlices[nSmearScan];
      TH1F *hSmearSigma[nSmearScan];
      TF1 *hFitSmearSigma[nSmearScan];
      for(int i=0; i<nSmearScan; i++)
	{
	  TH2F *htmp = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_%s_scan%d_%d",det_name[index],cent_Title[icent],0,i));
	  TH2F *h2 = (TH2F*)htmp->Clone(Form("%s_clone2",htmp->GetName()));
	  if(i>0) h2->RebinX(5);
	  hSmearSigma[i] = (TH1F*)h2->ProjectionX(Form("hSmearSigma_%d",i));
	  hSmearSigma[i]->Reset();
	  TCanvas *cSmearFit = 0x0;
	  if(makePlot)
	    {
	      cSmearFit = new TCanvas(Form("fit_smear_%d",i), Form("fit_smear_%d",i), 1200, 800);
	      cSmearFit->Divide(5,4);
	    }
	  for(int j=0; j<20; j++)
	    {
	      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hSmearMass_%d_%d",i,j),j+1,j+1);
	      h1->GetXaxis()->SetRangeUser(2.6,3.8);
	      TF1 *funcTmp = new TF1(Form("funcTemp_smear_%d_%d",i,j),"gaus",2.7,3.5);
	      funcTmp->SetLineColor(2);
	      h1->Fit(funcTmp,"IRQ0");
	      hSmearSigma[i]->SetBinContent(j+1,funcTmp->GetParameter(2));
	      hSmearSigma[i]->SetBinError(j+1,funcTmp->GetParError(2));
	      if(makePlot)
		{
		  cSmearFit->cd(j+1);
		  h1->Draw();
		  funcTmp->Draw("sames");
		}
	    }
	  hSmearSigma[i]->GetXaxis()->SetRangeUser(0,8);
	  hSmearSigma[i]->GetYaxis()->SetRangeUser(0.02,0.2);
	  hSmearSigma[i]->SetMarkerStyle(20+i);
	  hSmearSigma[i]->SetTitle(";p_{T} (GeV/c);#sigma");
	  cSmear->cd();
	  if(i==0) hSmearSigma[i]->Draw();
	  else     hSmearSigma[i]->Draw("sames");
	  hFitSmearSigma[i] = new TF1(Form("hFitSmearSigma_%d",i),"pol4",0,8);
	  hSmearSigma[i]->Fit(hFitSmearSigma[i],"IRQ0");
	  hFitSmearSigma[i]->SetLineColor(4);
	  hFitSmearSigma[i]->SetLineStyle(2);
	  hFitSmearSigma[i]->Draw("sames");
	}
      TH1F *hDataSigma = (TH1F*)fdata->Get(Form("Jpsi_FitSigma_cent%s",cent_Title[icent]));
      hDataSigma->SetMarkerStyle(21);
      hDataSigma->SetMarkerColor(2);
      hDataSigma->SetLineColor(2);
      hDataSigma->Draw("sames");
      t1 = GetTitleText("Width of J/#Psi mass peak");
      t1->Draw();
      t1 = GetPaveText(0.25,0.35,0.8,0.85);
      t1->AddText("Real data");
      t1->SetTextAlign(11);
      t1->SetTextColor(2);
      t1->Draw();
      t1 = GetPaveText(0.25,0.35,0.75,0.8);
      t1->AddText("Scan smearing");
      t1->SetTextAlign(11);
      t1->Draw();
      t1 = GetPaveText(0.25,0.35,0.7,0.75);
      t1->AddText("Fit to scan data");
      t1->SetTextAlign(11);
      t1->SetTextColor(4);
      t1->Draw();
      if(savePlot)
	{
	  cSmear->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ScanToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  cSmear->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ScanToyMcVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      printf("[i] Scan smear\n");
      TH1F *hSmearChi2[3];
      TF1 *hFitSmearChi2[3];
      for(int i=0; i<3; i++)
	{
	  hSmearChi2[i] = new TH1F(Form("hSmearChi2_%s",sysName[i]),";smear (%); #chi^{2}", nSmearScan, 0-smearStep/2, (nSmearScan-1)*smearStep+smearStep/2);
	  for(int j=0; j<nSmearScan; j++)
	    {
	      int bin = j + 1;
	      double chi2 = 0;
	      for(int ibin=1; ibin<=hDataSigma->GetNbinsX(); ibin++)
		{
		  double value = hDataSigma->GetBinContent(ibin);
		  double error = hDataSigma->GetBinError(ibin);
		  if(i==1) value -= error;
		  if(i==2) value += error;
		  chi2 += TMath::Power((value-hFitSmearSigma[j]->Eval(hDataSigma->GetBinCenter(ibin)))/error,2);
		}
	      hSmearChi2[i]->SetBinContent(bin,chi2);
	      hSmearChi2[i]->SetBinError(bin,1e-10);
	    }
	  c1 = new TCanvas(Form("smear_chi2_%s",sysName[i]), Form("smear_chi2_%s",sysName[i]), 800, 600);
	  hSmearChi2[i]->SetMarkerStyle(20);
	  hSmearChi2[i]->Draw("P");
	  TAxis *axis = hSmearChi2[i]->GetXaxis();
	  hFitSmearChi2[i] = new TF1(Form("hFitShiftChi2_%d",i),"pol4",axis->GetXmin(),axis->GetXmax());
	  hSmearChi2[i]->Fit(hFitSmearChi2[i],"IRQ0");
	  hFitSmearChi2[i]->SetLineColor(4);
	  hFitSmearChi2[i]->SetLineStyle(2);
	  hFitSmearChi2[i]->Draw("sames");
	  t1 = GetTitleText(Form("Scan smear: use %s of data points",sysTitle[i]));
	  t1->Draw();
	  t1 = GetPaveText(0.25,0.35,0.7,0.8);
	  t1->AddText(Form("Minimum #chi^{2} at %4.2f%%",hFitSmearChi2[i]->GetMinimumX()*100));
	  t1->Draw();
	  hFinalSmear->SetBinContent(i+1,hFitSmearChi2[i]->GetMinimumX());
	  if(savePlot)
	    {
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_Chi2_%s_cent%s.pdf",run_type,sysName[i],cent_Title[icent]));
	      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_Chi2_%s_cent%s.png",run_type,sysName[i],cent_Title[icent]));
	    }
	}

      fscan->cd();
      hFinalShift->Write("",TObject::kOverwrite);
      hFinalSmear->Write("",TObject::kOverwrite);
    }
  
  // real data
  /*
  TFile *fdata = 0x0;
  if(year==2014) fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  if(year==2013) fdata = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.VtxCut.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  TH1F *hRawCounts = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]));
  TH1F *hSignal = (TH1F*)fdata->Get(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));
  TF1 *fitSignal = (TF1*)fdata->Get(Form("FitJpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));

  const int nHistos = 5;
  TString name[nHistos] = {"TPCemb","TPCreco","MTDreco","MTDtrig","MTDresp"};

  const int shift_index = 5;
  const int smear_index = 10;
  
  if(anaType == 2)
    {


      // real data
      TH1F *hSignalFit = (TH1F*)hSignal->Clone("hSignalFit");
      hSignalFit->GetYaxis()->SetRangeUser(0,500);
      TCanvas *c = draw1D(hSignalFit);
      TF1 *funcSig = new TF1("funcSig","gausn(0)+pol3(3)",2.5,4);
      funcSig->SetParameter(0,20);
      funcSig->SetParameter(1,3.09);
      funcSig->SetParameter(2,0.06);
      hSignalFit->Fit(funcSig,"IR");
      funcSig->SetLineColor(2);
      funcSig->Draw("sames");
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Data_JpsiMass_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Data_JpsiMass_cent%s.png",run_type,cent_Title[icent]));
	}

      TGraphErrors *gMean = new TGraphErrors(1);
      gMean->SetPoint(0,0,funcSig->GetParameter(1));
      gMean->SetPointError(0,0,funcSig->GetParError(1));
      gMean->SetMarkerStyle(20);
      gMean->SetMarkerColor(2);
      gMean->SetLineColor(2);
      TGraphErrors *gSmear = new TGraphErrors(1);
      gSmear->SetPoint(0,0,funcSig->GetParameter(2));
      gSmear->SetPointError(0,0,funcSig->GetParError(2));
      gSmear->SetMarkerStyle(20);
      gSmear->SetMarkerColor(2);
      gSmear->SetLineColor(2);
      

      gStyle->SetStatFont(62);
      gStyle->SetStatFontSize(0.06);
      // determine shift
      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","read");
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.TrkResScan.root","read");
      TCanvas *cShift = new TCanvas("Fit_shift", "Fit_shift", 1200, 800);
      cShift->Divide(4,3);
      TH1F *hMean = new TH1F(Form("hMean_shift"),"Mean of J/#psi lineshape from MC;shift (%); mean",10, -1.9, 0.1);
      for(int i=0; i<nShiftScan; i++)
	{
	  TH2F *h2 = (TH2F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[index].Data(),cent_Title[icent],i,0));
	  TH1F *hRcJpsiMass = (TH1F*)h2->ProjectionY(Form("%s_proj",h2->GetName()));
	  cShift->cd(i+1);
	  hRcJpsiMass->GetXaxis()->SetRangeUser(2.8,3.6);
	  hRcJpsiMass->Draw();
	  TF1 *funcTmp = new TF1(Form("funcTemp_shift_%d",i),"gaus",3.0,3.2);
	  funcTmp->SetLineColor(2);
	  hRcJpsiMass->Fit(funcTmp,"IRQ");
	  hMean->SetBinContent(nShiftScan-i,funcTmp->GetParameter(1));
	  hMean->SetBinError(nShiftScan-i,funcTmp->GetParError(1));
	}
      hMean->SetMarkerStyle(21);
      hMean->GetYaxis()->SetRangeUser(3.03,3.1);
      c = draw1D(hMean);
      gMean->Draw("samesP");
      TF1 *funcMean = new TF1(Form("funcMean"),"pol1",hMean->GetXaxis()->GetXmin(),hMean->GetXaxis()->GetXmax());
      hMean->Fit(funcMean,"IR0");
      funcMean->SetLineColor(4);
      funcMean->Draw("sames");
      TLegend *leg = new TLegend(0.18,0.65,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gMean,"Real data","PE");
      leg->AddEntry(hMean,"Scan toyMC","PL");
      leg->AddEntry(funcMean,"Fit to scan","L");
      leg->Draw();
      double shift_y[3] = { funcSig->GetParameter(1), 
			    funcSig->GetParameter(1) - funcSig->GetParError(1), 
			    funcSig->GetParameter(1) + funcSig->GetParError(1)};
      for(int i=0; i<3; i++)
	{
	  TLine *line = GetLine(hMean->GetXaxis()->GetXmin(), shift_y[i],
				hMean->GetXaxis()->GetXmax(), shift_y[i],
				2,2,2);
	  line->Draw();
	}
      if(savePlot)
	{
	  cShift->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_FitScanShift_cent%s.pdf",run_type,cent_Title[icent]));
	  cShift->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_FitScanShift_cent%s.png",run_type,cent_Title[icent]));

	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_ScanShift_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_ScanShift_cent%s.png",run_type,cent_Title[icent]));
	}

      // determine smear
      TCanvas *cSmear = new TCanvas("Fit_smear", "Fit_smear", 1200, 800);
      cSmear->Divide(5,4);
      TH1F *hSmear = new TH1F(Form("hSmear_shift"),"Width of J/#psi lineshape from MC;smear (%); wdith", 20, -0.1, 3.9);
      for(int k=0; k<nSmearScan; k++)
	{
	  TH2F *h2 = (TH2F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[index].Data(),cent_Title[icent],shift_index,k));
	  TH1F *hRcJpsiMass = (TH1F*)h2->ProjectionY(Form("%s_proy",h2->GetName()));
	  cSmear->cd(k+1);
	  hRcJpsiMass->GetXaxis()->SetRangeUser(2.8,3.6);
	  hRcJpsiMass->Draw();
	  TF1 *funcTmp = new TF1(Form("funcTemp_smear_%d",k),"gaus",2.9,3.3);
	  funcTmp->SetLineColor(2);
	  hRcJpsiMass->Fit(funcTmp,"IRQ");
	  hSmear->SetBinContent(k+1,funcTmp->GetParameter(2));
	  hSmear->SetBinError(k+1,funcTmp->GetParError(2));
	}
      hSmear->SetMarkerStyle(21);
      c = draw1D(hSmear);
      gSmear->Draw("samesP");
      TF1 *funcSmear = new TF1(Form("funcSmear"),"pol3",hSmear->GetXaxis()->GetXmin(),hSmear->GetXaxis()->GetXmax());
      hSmear->Fit(funcSmear,"IR0");
      funcSmear->SetLineColor(4);
      funcSmear->Draw("sames");
      leg = new TLegend(0.18,0.65,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gSmear,"Real data","PE");
      leg->AddEntry(hSmear,"Scan toyMC","PL");
      leg->AddEntry(funcSmear,"Fit to scan","L");
      leg->Draw();
      double smear_y[3] = { funcSig->GetParameter(2), 
			    funcSig->GetParameter(2) - funcSig->GetParError(2), 
			    funcSig->GetParameter(2) + funcSig->GetParError(2)};
      for(int i=0; i<3; i++)
	{
	  TLine *line = GetLine(hSmear->GetXaxis()->GetXmin(), smear_y[i],
				hSmear->GetXaxis()->GetXmax(), smear_y[i],
				2,2,2);
	  line->Draw();
	}
      if(savePlot)
	{
	  cSmear->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_FitScanSmear_cent%s.pdf",run_type,cent_Title[icent]));
	  cSmear->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_FitScanSmear_cent%s.png",run_type,cent_Title[icent]));

	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_ScanSmear_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ToyMc_ScanSmear_cent%s.png",run_type,cent_Title[icent]));
	}

      // compare with real data
      TF1 *funcBkg = new TF1("funcBkg","pol3",2.5,4);
      for(int i=0; i<4; i++)
	{
	  funcBkg->SetParameter(i, funcSig->GetParameter(i+3));
	}
      TAxis *axis = hSignal->GetXaxis();
      TH1F *hSignalSub = new TH1F("hSignalSub","Compare J/#psi signal shape;M_{#mu#mu} (GeV/c^{2})",axis->GetNbins(),axis->GetXmin(), axis->GetXmax());
      for(int ibin=1; ibin<=hSignalSub->GetNbinsX(); ibin++)
	{
	  double value = hSignal->GetBinContent(ibin);
	  double error = hSignal->GetBinError(ibin);
	  double bkg = funcBkg->Eval(hSignal->GetBinCenter(ibin));
	  hSignalSub->SetBinContent(ibin, value - bkg);
	  hSignalSub->SetBinError(ibin, sqrt(error*error+fabs(bkg)));
	}
      hSignalSub->SetMarkerStyle(21);
      hSignalSub->GetXaxis()->SetRangeUser(2.5,4);
      hSignalSub->SetMaximum(1.6*hSignalSub->GetMaximum());
      hSignalSub->SetMinimum(-20);
      c = draw1D(hSignalSub);
      leg = new TLegend(0.18,0.68,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hSignalSub,"Real data","PE");
      double scale = hSignalSub->GetBinContent(hSignalSub->FindFixBin(3.07));
      TH1F *hRcJpsiMass = (TH1F*)((TH2F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[index].Data(),cent_Title[icent],0,0)))->ProjectionY("scan_0");
      hRcJpsiMass->Rebin(20);
      hRcJpsiMass->Scale(scale/hRcJpsiMass->GetBinContent(hRcJpsiMass->FindFixBin(3.1)));
      hRcJpsiMass->SetLineColor(4);
      hRcJpsiMass->SetLineWidth(1.5);
      hRcJpsiMass->Draw("sameHIST");
      leg->AddEntry(hRcJpsiMass,"Default embedding","L");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Compare_JpsiShape_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Compare_JpsiShape_cent%s.png",run_type,cent_Title[icent]));
	}

      hRcJpsiMass = (TH1F*)((TH2F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[index].Data(),cent_Title[icent],shift_index,smear_index)))->ProjectionY("scan_final");
      hRcJpsiMass->Rebin(20);
      hRcJpsiMass->Scale(scale/hRcJpsiMass->GetBinContent(hRcJpsiMass->FindFixBin(3.07)));
      hRcJpsiMass->SetLineColor(2);
      hRcJpsiMass->Draw("sameHIST");
      leg->AddEntry(hRcJpsiMass,"Embedding with smear & shift","L");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Check_JpsiShape_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Check_JpsiShape_cent%s.png",run_type,cent_Title[icent]));
	}
    }

  if(anaType == 3)
    {
      gStyle->SetStatFont(62);
      gStyle->SetStatFontSize(0.06);
      // determine shift
      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","read");
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.TrkResScan.root","read");
      TCanvas *cMC = new TCanvas("Fit_MC", "Fit_MC", 1200, 800);
      cMC->Divide(3,2);
      TH2F *hJpsiMassVspt = (TH2F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[index].Data(),cent_Title[icent],shift_index,smear_index));
      const int nbins = nPtBins -1;
      double xbins[nbins+1];
      for(int i=0; i<nbins; i++)
	xbins[i] = ptBins_low[i+1];
      xbins[nbins] = ptBins_high[nbins];
      TH1F *hMean  = new TH1F(Form("hMean_MC"),"Mean of J/#psi lineshape;p_{T} (GeV/c); mean (GeV/c^{2})",nbins,xbins);
      TH1F *hSigma = new TH1F(Form("hSigma_MC"),"Width of J/#psi lineshape;p_{T} (GeV/c); #sigma (GeV/c^{2})",nbins,xbins);
      for(int i=0; i<nbins; i++)
	{
	  int low_bin = hJpsiMassVspt->GetXaxis()->FindFixBin(xbins[i]+0.001);
	  int up_bin  = hJpsiMassVspt->GetXaxis()->FindFixBin(xbins[i+1]-0.001);
	  TH1F *hRcJpsiMass = (TH1F*)hJpsiMassVspt->ProjectionY(Form("hJpsiMass_%d",i),low_bin,up_bin);
	  cMC->cd(i+1);
	  hRcJpsiMass->Rebin(10);
	  hRcJpsiMass->GetXaxis()->SetRangeUser(2.8,3.6);
	  hRcJpsiMass->Draw();
	  TF1 *funcTmp = new TF1(Form("funcTemp_%d",i),"gaus",2.9,3.2);
	  funcTmp->SetLineColor(2);
	  hRcJpsiMass->Fit(funcTmp,"IRQ");
	  hMean->SetBinContent(i+1,funcTmp->GetParameter(1));
	  hMean->SetBinError(i+1,funcTmp->GetParError(1));
	  hSigma->SetBinContent(i+1,funcTmp->GetParameter(2));
	  hSigma->SetBinError(i+1,funcTmp->GetParError(2));
	}
      hMean->SetMarkerStyle(20);
      hMean->SetMarkerColor(2);
      hMean->SetLineColor(2);
      hMean->GetYaxis()->SetRangeUser(3.04,3.14);
      TCanvas *c = draw1D(hMean);
      TLegend *leg = new TLegend(0.18,0.68,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hMean,"Simulation","PL");
      TH1F *hMeanData = (TH1F*)fdata->Get("Jpsi_FitMean_cent00100");
      hMeanData->Draw("sames");
      leg->AddEntry(hMeanData,"Real data","PL");
      leg->Draw();

      hSigma->SetMarkerStyle(20);
      hSigma->SetMarkerColor(2);
      hSigma->SetLineColor(2);
      hSigma->GetYaxis()->SetRangeUser(0,0.12);
      c = draw1D(hSigma);
      leg = new TLegend(0.18,0.68,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hSigma,"Simulation","PL");
      TH1F *hSigmaData = (TH1F*)fdata->Get("Jpsi_FitSigma_cent00100");
      hSigmaData->Draw("sames");
      leg->AddEntry(hSigmaData,"Real data","PL");
      leg->Draw();
    }
  */
}

//================================================
double combinedFit(double *x, double *par)
{
  double xx = x[0];
  int bin = hJpsiLineShape->GetXaxis()->FindFixBin(xx);
  //int bin = TMath::Nint((xx-2)/hJpsiLineShape->GetBinWidth(1)) + 1;
  //cout << xx << "  " << bin << "  " << hJpsiLineShape->GetBinContent(bin) << endl;
  return par[0]*hJpsiLineShape->GetBinContent(bin) + par[1] + par[2]*xx + par[3]*xx*xx + par[4]*xx*xx*xx;
  //return par[0] + par[1] + par[2]*xx + par[3]*xx*xx + par[4]*xx*xx*xx;
  //return par[0]*hJpsiLineShape->GetBinContent(bin);
}

//================================================
void smear(const double mass, const int icent, const int nExpr, const double shift, const double sigma, TH1F *hInputPt, const bool debug, 
	   TH2F *hJpsiMassVsPt[gHistos])
{
  TF1 *funcsmear = (TF1*)funcTrkRes[icent]->Clone(Form("%s_smear",funcTrkRes[icent]->GetName()));
  funcsmear->SetParameter(0,funcTrkRes[icent]->GetParameter(0)*(1+sigma));

  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("\n+++ %d +++\n",i);
      double mc_pt =  hInputPt->GetRandom();
      //double mc_pt =  myRandom->Uniform(0,20);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.8, 0.8);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;
      hJpsiMassVsPt[0]->Fill(parent.Pt(),parent.M());
      
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
      int binx1 = hTpcEff[icent][index1]->GetXaxis()->FindBin(pt1);
      if(pt1>10) binx1 = hTpcEff[icent][index1]->GetXaxis()->FindBin(10);
      double probility1 = hTpcEff[icent][index1]->GetBinContent(binx1,
								hTpcEff[icent][index1]->GetYaxis()->FindBin(eta1),
								hTpcEff[icent][index1]->GetZaxis()->FindBin(phi1_new));
      if(debug) printf("Efficiency 1 = %3.2f\n",probility1);
      if( myRandom->Uniform(0., 1.) >  probility1) continue;

      int index2 = (charge2 > 0) ? 0 : 1;
      double phi2_new = phi2;
      if(phi2_new<0) phi2_new += 2*pi;
      phi2_new += pi/12;
      int binx2 = hTpcEff[icent][index2]->GetXaxis()->FindBin(pt2);
      if(pt2>10) binx2 = hTpcEff[icent][index2]->GetXaxis()->FindBin(10);
      double probility2 = hTpcEff[icent][index2]->GetBinContent(binx2,
								hTpcEff[icent][index2]->GetYaxis()->FindBin(eta2),
								hTpcEff[icent][index2]->GetZaxis()->FindBin(phi2_new));
      if(debug) printf("Efficiency 2 = %3.2f\n",probility2);
      if( myRandom->Uniform(0., 1.) >  probility2) continue;

      // momentum resolution & shift
      int mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt1)-1;
      if(pt1>10) mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;
      int mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt2)-1;
      if(pt2>10) mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;

      double dpt1 = hTrkResBin[icent][mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[icent][mom_index2]->GetRandom();

      double rc_emb_pt1 = (1-dpt1) * pt1;
      double rc_emb_pt2 = (1-dpt2) * pt2;

      double rc_pt1 = rc_emb_pt1 * myRandom->Gaus(1+shift/sqrt(rc_emb_pt1),rc_emb_pt1*sigma);
      double rc_pt2 = rc_emb_pt2 * myRandom->Gaus(1+shift/sqrt(rc_emb_pt2),rc_emb_pt2*sigma);
      

      /*
      double res1 = funcTrkRes[icent]->Eval(pt1);
      double res2 = funcTrkRes[icent]->Eval(pt2);
      double res1_smear = funcsmear->Eval(pt1);
      double res2_smear = funcsmear->Eval(pt2);
      double rc_pt1 = (1-dpt1/res1*res1_smear) * pt1 * (1+shift);
      double rc_pt2 = (1-dpt2/res2*res2_smear) * pt2 * (1+shift);
      */

      if(debug) printf("rc daug 1: pt = %3.2f\n",rc_pt1);
      if(debug) printf("rc daug 2: pt = %3.2f\n",rc_pt2);
      if(rc_pt1<pt1_cut || rc_pt2<pt2_cut) continue;
      TLorentzVector rc_emb_daughter1, rc_emb_daughter2;
      rc_emb_daughter1.SetPtEtaPhiM(rc_emb_pt1,eta1,phi1,muMass);
      rc_emb_daughter2.SetPtEtaPhiM(rc_emb_pt2,eta2,phi2,muMass);
      TLorentzVector rc_emb_parent = rc_emb_daughter1 + rc_emb_daughter2;
      hJpsiMassVsPt[1]->Fill(rc_emb_parent.Pt(),rc_emb_parent.M());

      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hJpsiMassVsPt[2]->Fill(rc_parent.Pt(),rc_parent.M());

      // MTD matching
      int tray_index = 1;
      if(phi1<3*pi/15 && phi1>12*pi/15) tray_index = 0;
      binx1 = hMtdEff[icent][tray_index][index1]->GetXaxis()->FindBin(rc_pt1);
      if(rc_pt1>10) binx1 = hMtdEff[icent][tray_index][index1]->GetXaxis()->FindBin(10);
      double mtd_eff_1 = hMtdEff[icent][tray_index][index1]->GetBinContent(binx1,
									   hMtdEff[icent][tray_index][index1]->GetYaxis()->FindFixBin(eta1));
      if(debug) printf("MTD eff 1 = %3.2f at bin = (%d,%d)\n",mtd_eff_1,binx1,hMtdEff[icent][tray_index][index1]->GetYaxis()->FindFixBin(eta1));
      if( myRandom->Uniform(0., 1.) >  mtd_eff_1) continue;
      
      binx2 = hMtdEff[icent][tray_index][index2]->GetXaxis()->FindBin(rc_pt2);
      if(rc_pt2>10) binx2 = hMtdEff[icent][tray_index][index2]->GetXaxis()->FindBin(10);
      double mtd_eff_2 = hMtdEff[icent][tray_index][index2]->GetBinContent(binx2,
									   hMtdEff[icent][tray_index][index2]->GetYaxis()->FindFixBin(eta2));
      if(debug) printf("MTD eff 2 = %3.2f\n",mtd_eff_2);
      if( myRandom->Uniform(0., 1.) >  mtd_eff_2) continue;
      hJpsiMassVsPt[3]->Fill(rc_parent.Pt(),rc_parent.M());

      // MTD efficiency
      if(funcTrigEff[icent] && funcTrigEffCorr[icent])
	{
	  double trig_eff_1 = funcTrigEff[icent]->Eval(rc_pt1) * funcTrigEffCorr[icent]->Eval(rc_pt1);
	  if(debug) printf("MTD trig eff 1 = %3.2f\n",trig_eff_1);
	  if( myRandom->Uniform(0., 1.) >  trig_eff_1) continue;

	  double trig_eff_2 = funcTrigEff[icent]->Eval(rc_pt2) * funcTrigEffCorr[icent]->Eval(rc_pt2);
	  if(debug) printf("MTD trig eff 2 = %3.2f\n",trig_eff_2);
	  if( myRandom->Uniform(0., 1.) >  trig_eff_2) continue;
	}
      hJpsiMassVsPt[4]->Fill(rc_parent.Pt(),rc_parent.M());

      // MTD response efficiency
      int resp_index1 = hMtdRespEffCosmic->FindBin(rc_pt1);
      if(rc_pt1>10) resp_index1 = hMtdRespEffCosmic->FindBin(9.5);
      double mtd_resp_eff_1 = hMtdRespEffCosmic->GetBinContent(resp_index1)/hMtdRespEffEmbed->GetBinContent(resp_index1);
      if(debug) printf("MTD resp eff 1 = %3.2f\n",mtd_resp_eff_1);
      if( myRandom->Uniform(0., 1.) >  mtd_resp_eff_1) continue;
	  
      int resp_index2 = hMtdRespEffCosmic->FindBin(rc_pt2);
      if(rc_pt2>10) resp_index1 = hMtdRespEffCosmic->FindBin(9.5);
      double mtd_resp_eff_2 = hMtdRespEffCosmic->GetBinContent(resp_index2)/hMtdRespEffEmbed->GetBinContent(resp_index2);
      if(debug) printf("MTD resp eff 2 = %3.2f\n",mtd_resp_eff_2);
      if( myRandom->Uniform(0., 1.) >  mtd_resp_eff_2) continue;
    
      hJpsiMassVsPt[5]->Fill(rc_parent.Pt(),rc_parent.M());
    }
}

//-------------------------------------------------------
void tuneSmear(const double masss, const int icent, const int nExpr, 
               const double shift, const double sigma, 
               TH1F *hInputPt, TH2F *hRcJpsiMassVsPt)
{
  TF1 *funcsmear = (TF1*)funcTrkRes[icent]->Clone(Form("%s_smear",funcTrkRes[icent]->GetName()));
  funcsmear->SetParameter(0,funcTrkRes[icent]->GetParameter(0)*(1+sigma));

  for(int i=0; i<nExpr; i++)
    {
      //double mc_pt =  hInputPt->GetRandom();
      double mc_pt  = myRandom->Uniform(0,20);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.8, 0.8);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+masss*masss) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,masss);
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;

      double pt1  = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      
      double pt2  = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();

      if(pt1<0.5 || pt2<0.5) continue;

      // momentum resolution & shift
      int mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt1)-1;
      if(pt1>10) mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;
      int mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt2)-1;
      if(pt2>10) mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;

      double dpt1 = hTrkResBin[icent][mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[icent][mom_index2]->GetRandom();
      double emb_pt1 = (1-dpt1) * pt1;
      double emb_pt2 = (1-dpt2) * pt2;

      double rc_pt1 = emb_pt1 * myRandom->Gaus(1+shift/sqrt(emb_pt1),emb_pt1*sigma);
      double rc_pt2 = emb_pt2 * myRandom->Gaus(1+shift/sqrt(emb_pt2),emb_pt2*sigma);

      // double res1 = funcTrkRes[icent]->Eval(pt1);
      // double res2 = funcTrkRes[icent]->Eval(pt2);
      // double res1_smear = funcsmear->Eval(pt1);
      // double res2_smear = funcsmear->Eval(pt2);
      // double rc_pt1 = (1-dpt1/res1*res1_smear) * pt1 * (1+shift/sqrt(pt1));
      // double rc_pt2 = (1-dpt2/res2*res2_smear) * pt2 * (1+shift/sqrt(pt2));


      if(rc_pt1<pt1_cut || rc_pt2<pt2_cut) continue;
      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hRcJpsiMassVsPt->Fill(rc_parent.Pt(),rc_parent.M());
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
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size, const Int_t font)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}

//--------------------------------------------
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color, Width_t width, Style_t style)
{
  TLine *line = new TLine(xl,yl,xh,yh);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  return line;
}
