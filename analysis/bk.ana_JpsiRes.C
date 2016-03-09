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
#else if (YEAR==2013)
// +++ Run13 analysis +++
const char *run_config = "";
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
const double jpsiMass = 3.097;
const double muMass = 0.1057;
const int year = YEAR;

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

void flatPt(const int icent, const bool savePlot);
void tuneResolution(const int icent, const bool savePlot);
void upsilon(const int icent, const bool savePlot);
void trigEff();
void smear(const double masss, const int icent, const int nExpr, const double shift, const double sigma, TH1F *hInputPt, 
	   TH1F *hRcJpsiPtEmbed, TH2F *hRcJpsiMassEmbed,
	   TH1F *hRcJpsiPtTpc, TH2F *hRcJpsiMassTpc, 
	   TH1F *hRcJpsiPtMtd, TH2F *hRcJpsiMassMtd, 
	   TH1F *hRcJpsiPtMtdTrig, TH2F *hRcJpsiMassMtdTrig,
	   TF1 *funcTrigEff, TF1 *funcTrigEffCorr,
	   const bool debug,
	   TH1F *hRcJpsiPtMtdResp=0, TH2F *hRcJpsiMassMtdResp=0);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);
double combinedFit(double *x, double *par);

TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
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

  //flatPt(0, 0);
  tuneResolution(0,1);
  //upsilon(0,1);
  //trigEff();
}

//================================================
void upsilon(const int icent, const bool savePlot)
{
  TString name[3] = {"TPCemb","TPCreco","MTDreco"};
  TFile *fin = TFile::Open("Rootfiles/Upsion.pp200.root","read");
  TH1F *hMcJpsiPt = (TH1F*)fin->Get("Upsilon_pp200_midRap");
  double upsilonMass[3] = {9.46,10.023,10.355};

  TH1F *hRcJpsiPt[3][3];
  TH2F *hRcJpsiMass[3][3];
  for(int j=0; j<3; j++)
    {
      for(int i=0; i<3; i++)
	{
	  hRcJpsiPt[j][i] = new TH1F(Form("hRcPt_Upsilon%dS_%s_%s",j+1,name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed #Upsilon%dS (%s%%);p_{T} (GeV/c)",j+1,cent_Name[icent]),110,0,11);
	  hRcJpsiMass[j][i] = new TH2F(Form("hRcMass_Upsilon%dS_%s_%s",j+1,name[i].Data(),cent_Title[icent]),Form("mass distribution of reconstructed #Upsilon%dS (%s%%);p_{T} (GeV/c)",j+1,cent_Name[icent]),110,0,11,400,8,12);
	}

      smear(upsilonMass[j], icent, 1e8, 0, 0.017, hMcJpsiPt, hRcJpsiPt[j][0], hRcJpsiMass[j][0],hRcJpsiPt[j][1], hRcJpsiMass[j][1], hRcJpsiPt[j][2], hRcJpsiMass[j][2], 0, 0, 0, 0, 0, 0, 0);
    }

  for(int j=0; j<3; j++)
    {
      for(int i=1; i<2; i++)
	{
	  hRcJpsiPt[j][i]->Sumw2();
	  hRcJpsiPt[j][i]->SetLineColor(2);
	  TCanvas *c = draw1D(hRcJpsiPt[j][i],Form("%s: p_{T} distribution of #Upsilon%dS (%s%%)",name[i].Data(),j+1,cent_Name[icent]),kTRUE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_Pt_cent%s.pdf",run_type,name[i].Data(),j+1,cent_Title[icent]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_Pt_cent%s.png",run_type,name[i].Data(),j+1,cent_Title[icent]));
	    }

	  hRcJpsiMass[j][i]->Sumw2();
	  TH1F *hMass = (TH1F*)hRcJpsiMass[j][i]->ProjectionY(Form("%s_proj",hRcJpsiMass[j][i]->GetName()));
	  hMass->SetLineColor(2);
	  c = draw1D(hMass,Form("%s: mass distribution of #Upsilon%dS (%s%%)",name[i].Data(),j+1,cent_Name[icent]),kFALSE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_mass_cent%s.pdf",run_type,name[i].Data(),j+1,cent_Title[icent]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/%s_Upsilon%dS_mass_cent%s.png",run_type,name[i].Data(),j+1,cent_Title[icent]));
	    }
	}
    }

  if(1)
    {
      cout << "Save histogram" << endl;
      TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.Upsilon.LineShape.root","recreate");
      for(int j=0; j<3; j++)
	{
	  hRcJpsiMass[j][1]->Write();
	}
    } 
}

//================================================
void trigEff()
{
  const int nHistos = 5;
  TString name[nHistos] = {"TPCemb","TPCreco","MTDreco","MTDtrig","MTDresp"};
  TFile *fin = 0x0;
  if(year==2013) fin = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
  if(year==2014) fin = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");

  TF1 *funcMcJpsiPt[nCentBins];
  TH1F *hMcJpsiPt[nCentBins];
  TH1F *hRcJpsiPt[nCentBins][nHistos];
  TH2F *hRcJpsiMass[nCentBins][nHistos];

  TH1F *hJpsiSmearEff[nCentBins][2];
  TH1F *hJpsiTrigEff[nCentBins][2];
  TH1F *hJpsiRespEff[nCentBins][2];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  // const int nbins = 10;
  // double xbins[nbins+1] = {0,1,2,3,4,5,6,7,8,9,10};

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

      for(int i=0; i<nHistos; i++)
	{
	  hRcJpsiPt[icent][i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
	  hRcJpsiPt[icent][i]->Sumw2();
	  hRcJpsiMass[icent][i] = new TH2F(Form("hRcJpsiMass_%s_%s",name[i].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),110,0,11,4000,2,4);
	  hRcJpsiMass[icent][i]->Sumw2();
	}
      smear(jpsiMass, 0, 4e8, shift, sigma, hMcJpsiPt[icent], hRcJpsiPt[icent][0], hRcJpsiMass[icent][0], hRcJpsiPt[icent][1], hRcJpsiMass[icent][1], hRcJpsiPt[icent][2], hRcJpsiMass[icent][2], hRcJpsiPt[icent][3], hRcJpsiMass[icent][3],funcTrigEff[0], funcTrigEffCorr[0], 0, hRcJpsiPt[icent][4], hRcJpsiMass[icent][4]);

      hJpsiSmearEff[icent][0] = (TH1F*)hRcJpsiPt[icent][1]->Clone(Form("JpsiSmearEff_cent%s",cent_Title[icent]));
      hJpsiSmearEff[icent][0]->Divide(hRcJpsiPt[icent][0]);
      hJpsiSmearEff[icent][1] = (TH1F*)hRcJpsiPt[icent][1]->Rebin(nbins,Form("JpsiSmearEff_cent%s_rebin",cent_Title[icent]),xbins);
      TH1F *htmp = (TH1F*)hRcJpsiPt[icent][0]->Rebin(nbins,Form("JpsiPtSmear_cent%s_rebin",cent_Title[icent]),xbins);
      hJpsiSmearEff[icent][1]->Divide(htmp);
      hJpsiSmearEff[icent][1]->SetMarkerStyle(20);
      hJpsiSmearEff[icent][1]->SetMarkerColor(2);
      TCanvas *c = draw1D(hJpsiSmearEff[icent][0]);
      hJpsiSmearEff[icent][1]->Draw("sames");

      hJpsiTrigEff[icent][0] = (TH1F*)hRcJpsiPt[icent][3]->Clone(Form("JpsiTrigEff_cent%s",cent_Title[icent]));
      hJpsiTrigEff[icent][0]->Divide(hRcJpsiPt[icent][2]);
      hJpsiTrigEff[icent][1] = (TH1F*)hRcJpsiPt[icent][3]->Rebin(nbins,Form("JpsiTrigEff_cent%s_rebin",cent_Title[icent]),xbins);
      htmp = (TH1F*)hRcJpsiPt[icent][2]->Rebin(nbins,Form("JpsiPtMtd_cent%s_rebin",cent_Title[icent]),xbins);
      hJpsiTrigEff[icent][1]->Divide(htmp);
      hJpsiTrigEff[icent][1]->SetMarkerStyle(20);
      hJpsiTrigEff[icent][1]->SetMarkerColor(2);
      c = draw1D(hJpsiTrigEff[icent][0]);
      hJpsiTrigEff[icent][1]->Draw("sames");

      hJpsiRespEff[icent][0] = (TH1F*)hRcJpsiPt[icent][4]->Clone(Form("JpsiRespEff_cent%s",cent_Title[icent]));
      hJpsiRespEff[icent][0]->Divide(hRcJpsiPt[icent][3]);
      hJpsiRespEff[icent][1] = (TH1F*)hRcJpsiPt[icent][4]->Rebin(nbins,Form("JpsiRespEff_cent%s_rebin",cent_Title[icent]),xbins);
      htmp = (TH1F*)hRcJpsiPt[icent][3]->Rebin(nbins,Form("JpsiRespEff_cent%s_rebin",cent_Title[icent]),xbins);
      hJpsiRespEff[icent][1]->Divide(htmp);
      hJpsiRespEff[icent][1]->SetMarkerStyle(20);
      hJpsiRespEff[icent][1]->SetMarkerColor(2);
      c = draw1D(hJpsiRespEff[icent][0]);
      hJpsiRespEff[icent][1]->Draw("sames");
    }
 
  TString outName = "";
  if(year==2014) outName = Form("Run14.AuAu200.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  if(year==2013) outName = Form("Run13.pp500.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  TFile *fEff = TFile::Open(Form("Rootfiles/%s",outName.Data()),"update");

  for(int icent=0; icent<nCentBins; icent++)
    {
      for(int i=0; i<2; i++)
	{
	  hJpsiSmearEff[icent][i]->Write("",TObject::kOverwrite);
	  hJpsiTrigEff[icent][i]->Write("",TObject::kOverwrite);
	  hJpsiRespEff[icent][i]->Write("",TObject::kOverwrite);
	}
    }
  
}

//================================================
void tuneResolution(const int icent, const bool savePlot)
{
  const int anaType = 3; // 0 - check raw jpsi distribution; 1 - scan; 2 - check scan; 3 - pt dependence
  const int nShiftScan = 10;
  const int nSmearScan = 20;
  const int index = 4;

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

  // real data
  TFile *fdata = 0x0;
  if(year==2014) fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  if(year==2013) fdata = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.VtxCut.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  TH1F *hRawCounts = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]));
  TH1F *hSignal = (TH1F*)fdata->Get(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));
  TF1 *fitSignal = (TF1*)fdata->Get(Form("FitJpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[0]));


  const int nHistos = 5;
  TString name[nHistos] = {"TPCemb","TPCreco","MTDreco","MTDtrig","MTDresp"};
  // check raw jpsi distribution
  if(anaType == 0)
    {
      TH1F *hRcJpsiPt[nHistos];
      TH2F *hRcJpsiMass[nHistos];
      for(int i=0; i<nHistos; i++)
	{
	  hRcJpsiPt[i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
	  hRcJpsiMass[i] = new TH2F(Form("hRcJpsiMass_%s_%s",name[i].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),110,0,11,4000,2,4);
	}
      smear(jpsiMass, icent, 5e7, 0, 0, hMcJpsiPt, hRcJpsiPt[0], hRcJpsiMass[0], hRcJpsiPt[1], hRcJpsiMass[1], hRcJpsiPt[2], hRcJpsiMass[2], hRcJpsiPt[3], hRcJpsiMass[3], funcTrigEff[icent], funcTrigEffCorr[icent], 0, hRcJpsiPt[4], hRcJpsiMass[4]);

      hRawCounts->SetMarkerStyle(21);
      hRawCounts->Scale(1./hRawCounts->Integral());
      TCanvas *c = draw1D(hRawCounts,Form("p_{T} distribution of J/#psi in MTD"),kFALSE,kTRUE);
      hRcJpsiPt[index]->Sumw2();
      hRcJpsiPt[index]->SetLineColor(2);
      TH1F *hRebin = (TH1F*)hRcJpsiPt[index]->Rebin(hRawCounts->GetXaxis()->GetNbins(),Form("%s_rebin",hRawCounts->GetName()),hRawCounts->GetXaxis()->GetXbins()->GetArray());
      hRebin->Scale(1./hRebin->Integral());
      hRebin->Draw("sameHIST");
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_ToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_ToyMcVsData_cent%s.png",run_type,cent_Title[icent]));
	}
    }  

  TH1F *hRcJpsiPtScan[nShiftScan][nSmearScan][nHistos];
  TH2F *hRcJpsiMassScan[nShiftScan][nSmearScan][nHistos];
  if(anaType == 1)
    {
      for(int i=0; i<nShiftScan; i++)
	{
	  for(int k=0; k<nSmearScan; k++)
	    {
	      printf("[i] Scan (%d,%d)\n",i+1,k+1);
	      double sigma = k*0.002;
	      double shift = i*(-0.002);
	      for(int j=0; j<nHistos; j++)
		{
		  hRcJpsiPtScan[i][k][j] = new TH1F(Form("hRcJpsiPt_%s_%s_scan%d_%d",name[j].Data(),cent_Title[icent],i,k),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
		  hRcJpsiMassScan[i][k][j] = new TH2F(Form("hRcJpsiMass_%s_%s_scan%d_%d",name[j].Data(),cent_Title[icent],i,k),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),110,0,11,4000,2,4);
		}
	      smear(jpsiMass, icent, 1e7, shift, sigma, hMcJpsiPt, 
		    hRcJpsiPtScan[i][k][0], hRcJpsiMassScan[i][k][0], hRcJpsiPtScan[i][k][1], hRcJpsiMassScan[i][k][1],
		    hRcJpsiPtScan[i][k][2], hRcJpsiMassScan[i][k][2], hRcJpsiPtScan[i][k][3], hRcJpsiMassScan[i][k][3],
		    funcTrigEff[icent], funcTrigEffCorr[icent], 0, hRcJpsiPtScan[i][k][4], hRcJpsiMassScan[i][k][4]);
	    }
	}

      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","recreate");
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.TrkResScan.root","recreate");
      for(int i=0; i<nShiftScan; i++)
	{
	  for(int k=0; k<nSmearScan; k++)
	    {
	      for(int j=0; j<nHistos; j++)
		{
		  hRcJpsiPtScan[i][k][j]->Write();
		  hRcJpsiMassScan[i][k][j]->Write();
		}
	    }
	}
    }

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
void smear(const double masss, const int icent, const int nExpr, const double shift, const double sigma, TH1F *hInputPt, 
	   TH1F *hRcJpsiPtEmbed, TH2F *hRcJpsiMassEmbed, TH1F *hRcJpsiPtTpc, TH2F *hRcJpsiMassTpc, 
	   TH1F *hRcJpsiPtMtd, TH2F *hRcJpsiMassMtd, TH1F *hRcJpsiPtMtdTrig, TH2F *hRcJpsiMassMtdTrig,
	   TF1 *funcTrigEff, TF1 *funcTrigEffCorr,
	   const bool debug,
	   TH1F *hRcJpsiPtMtdResp, TH2F *hRcJpsiMassMtdResp)
{
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("\n+++ %d +++\n",i);
      double mc_pt =  hInputPt->GetRandom();
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.8, 0.8);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+masss*masss) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,masss);
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

      double res1 = funcTrkRes[icent]->Eval(pt1);
      double res2 = funcTrkRes[icent]->Eval(pt2);
      
      TF1 *functmp = (TF1*)funcTrkRes[icent]->Clone(Form("%s_%d",funcTrkRes[icent]->GetName(),i));
      functmp->SetParameter(0,funcTrkRes[icent]->GetParameter(0)*(1+sigma));
      double res1_smear = functmp->Eval(pt1);
      double res2_semar = functmp->Eval(pt2);

      double rc_pt1 = (1-dpt1/res1*res1_smear) * pt1 * (1+shift);
      double rc_pt2 = (1-dpt2/res2*res2_smear) * pt2 * (1+shift);

      if(debug) printf("rc daug 1: pt = %3.2f\n",rc_pt1);
      if(debug) printf("rc daug 2: pt = %3.2f\n",rc_pt2);
      if(rc_pt1<pt1_cut || rc_pt2<pt2_cut) continue;
      TLorentzVector rc_emb_daughter1, rc_emb_daughter2;
      rc_emb_daughter1.SetPtEtaPhiM(rc_emb_pt1,eta1,phi1,muMass);
      rc_emb_daughter2.SetPtEtaPhiM(rc_emb_pt2,eta2,phi2,muMass);
      TLorentzVector rc_emb_parent = rc_emb_daughter1 + rc_emb_daughter2;
      hRcJpsiPtEmbed->Fill(rc_emb_parent.Pt());
      hRcJpsiMassEmbed->Fill(rc_emb_parent.Pt(),rc_emb_parent.M());

      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hRcJpsiPtTpc->Fill(rc_parent.Pt());
      hRcJpsiMassTpc->Fill(rc_parent.Pt(),rc_parent.M());

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
      hRcJpsiPtMtd->Fill(rc_parent.Pt());
      hRcJpsiMassMtd->Fill(rc_parent.Pt(),rc_parent.M());

      // MTD efficiency
      if(funcTrigEff && funcTrigEffCorr)
	{
	  double trig_eff_1 = funcTrigEff->Eval(rc_pt1) * funcTrigEffCorr->Eval(rc_pt1);
	  if(debug) printf("MTD trig eff 1 = %3.2f\n",trig_eff_1);
	  if( myRandom->Uniform(0., 1.) >  trig_eff_1) continue;

	  double trig_eff_2 = funcTrigEff->Eval(rc_pt2) * funcTrigEffCorr->Eval(rc_pt2);
	  if(debug) printf("MTD trig eff 2 = %3.2f\n",trig_eff_2);
	  if( myRandom->Uniform(0., 1.) >  trig_eff_2) continue;
	}
      if(hRcJpsiPtMtdTrig) hRcJpsiPtMtdTrig->Fill(rc_parent.Pt());
      if(hRcJpsiMassMtdTrig) hRcJpsiMassMtdTrig->Fill(rc_parent.Pt(),rc_parent.M());

      // MTD response efficiency
      if(hMtdRespEffCosmic && hMtdRespEffEmbed)
	{
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
	}
      if(hRcJpsiPtMtdResp) hRcJpsiPtMtdResp->Fill(rc_parent.Pt());
      if(hRcJpsiMassMtdResp) hRcJpsiMassMtdResp->Fill(rc_parent.Pt(),rc_parent.M());
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

//--------------------------------------------
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color, Width_t width, Style_t style)
{
  TLine *line = new TLine(xl,yl,xh,yh);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  return line;
}
