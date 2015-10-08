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

#include <cmath>
using namespace std;

const double pi = 3.1415926;
const double pt1_cut = 1.0, pt2_cut = 1.0;
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
const int nPtBins = 7;
const double ptBins_low[nPtBins]  = {0,1,2,3,4,6,8};
const double ptBins_high[nPtBins] = {10,2,3,4,6,8,10};

//const char *run_type = "Run13_pp500";

const char *run_type = "Run14_AuAu200";
const int year = 2014;
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
void smear(const double masss, const int icent, const int nExpr, const double sigma, TH1F *hInputPt, 
	   TH1F *hRcJpsiPtTpc, TH1F *hRcJpsiMassTpc, TH1F *hRcJpsiPtMtd, TH1F *hRcJpsiMassMtd, 
	   TH1F *hRcJpsiPtMtdTrig, TH1F *hRcJpsiMassMtdTrig,
	   TF1 *funcTrigEff, TF1 *funcTrigEffCorr,
	   const bool debug,
	   TH1F *hRcJpsiPtMtdResp=0, TH1F *hRcJpsiMassMtdResp=0);
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
  TFile *fEff = TFile::Open("Rootfiles/Run14.AuAu200.TrkEff.root","read");
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
  TFile *fTrig =  TFile::Open(Form("Rootfiles/Run14.AuAu200.MuonTrigEff.root"),"read");
  for(int icent=0; icent<nCentBins; icent++)
    {
      funcTrigEff[icent] = (TF1*)fTrig->Get(Form("MuonTrigEff_cent%s",cent_Title[icent]));
      funcTrigEffCorr[icent] = (TF1*)fTrig->Get(Form("MuonTrigEffCorr_cent%s",cent_Title[icent]));
    }

  // MTD response efficiency
  TFile *fCosmic =  TFile::Open(Form("Rootfiles/Run14.AuAu200.MtdResponseEff.root"),"read");
  hMtdRespEffCosmic = (TH1F*)fCosmic->Get("MtdResponseEff_cosmic");
  hMtdRespEffEmbed = (TH1F*)fCosmic->Get("MtdResponseEff_embed");

  //flatPt(0, 0);
  //tuneResolution(0,1);
  //upsilon(0,1);
  trigEff();
}

//================================================
void upsilon(const int icent, const bool savePlot)
{
  TString name[2] = {"TPCreco","MTDreco"};
  TFile *fin = TFile::Open("Rootfiles/Upsion.pp200.root","read");
  TH1F *hMcJpsiPt = (TH1F*)fin->Get("Upsilon_pp200_midRap");
  double upsilonMass[3] = {9.46,10.023,10.355};

  TH1F *hRcJpsiPt[3][2];
  TH1F *hRcJpsiMass[3][2];
  for(int j=0; j<3; j++)
    {
      for(int i=0; i<2; i++)
	{
	  hRcJpsiPt[j][i] = new TH1F(Form("hRcPt_Upsilon%dS_%s_%s",j+1,name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed #Upsilon%dS (%s%%);p_{T} (GeV/c)",j+1,cent_Name[icent]),110,0,11);
	  hRcJpsiMass[j][i] = new TH1F(Form("hRcMass_Upsilon%dS_%s_%s",j+1,name[i].Data(),cent_Title[icent]),Form("mass distribution of reconstructed #Upsilon%dS (%s%%);p_{T} (GeV/c)",j+1,cent_Name[icent]),400,8,12);
	}

      smear(upsilonMass[j], icent, 1e8, 0.017, hMcJpsiPt, hRcJpsiPt[j][0], hRcJpsiMass[j][0], hRcJpsiPt[j][1], hRcJpsiMass[j][1], 0, 0, 0, 0, 0);
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
	  hRcJpsiMass[j][i]->SetLineColor(2);
	  //hRcJpsiMass[j][i]->Scale(1./hRcJpsiMass[j][i]->Integral());
	  //hRcJpsiMass[j][i]->GetXaxis()->SetRangeUser(2.8,3.4);
	  c = draw1D(hRcJpsiMass[j][i],Form("%s: mass distribution of #Upsilon%dS (%s%%)",name[i].Data(),j+1,cent_Name[icent]),kFALSE,kFALSE);
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
  TString name[4] = {"TPCreco","MTDreco","MTDtrig","MTDresp"};
  TFile *fin = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");

  TF1 *funcMcJpsiPt[nCentBins];
  TH1F *hMcJpsiPt[nCentBins];
  TH1F *hRcJpsiPt[nCentBins][4];
  TH1F *hRcJpsiMass[nCentBins][4];
  TH1F *hJpsiTrigEff[nCentBins][2];
  TH1F *hJpsiRespEff[nCentBins][2];

  // const int nbins = nPtBins -1;
  // double xbins[nbins+1];
  // for(int i=0; i<nbins; i++)
  //   xbins[i] = ptBins_low[i+1];
  // xbins[nbins] = ptBins_high[nbins];
  const int nbins = 10;
  double xbins[nbins+1] = {0,1,2,3,4,5,6,7,8,9,10};

  for(int icent=0; icent<nCentBins; icent++)
    {
      cout << "Cent " << icent << endl;
      hMcJpsiPt[icent] = (TH1F*)fin->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[icent]));

      for(int i=0; i<4; i++)
	{
	  hRcJpsiPt[icent][i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
	  hRcJpsiPt[icent][i]->Sumw2();
	  hRcJpsiMass[icent][i] = new TH1F(Form("hRcJpsiMass_%s_%s",name[i].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),4000,2,4);
	  hRcJpsiMass[icent][i]->Sumw2();
	}
      smear(jpsiMass, 0, 4e8, 0.016, hMcJpsiPt[icent], hRcJpsiPt[icent][0], hRcJpsiMass[icent][0], hRcJpsiPt[icent][1], hRcJpsiMass[icent][1], hRcJpsiPt[icent][2], hRcJpsiMass[icent][2], funcTrigEff[0], funcTrigEffCorr[0], 0, hRcJpsiPt[icent][3], hRcJpsiMass[icent][3]);

      hJpsiTrigEff[icent][0] = (TH1F*)hRcJpsiPt[icent][2]->Clone(Form("JpsiTrigEff_cent%s",cent_Title[icent]));
      hJpsiTrigEff[icent][0]->Divide(hRcJpsiPt[icent][1]);
      TCanvas *c = draw1D(hJpsiTrigEff[icent][0]);

      hJpsiTrigEff[icent][1] = (TH1F*)hRcJpsiPt[icent][2]->Rebin(nbins,Form("JpsiTrigEff_cent%s_rebin",cent_Title[icent]),xbins);
      TH1F *htmp = (TH1F*)hRcJpsiPt[icent][1]->Rebin(nbins,Form("JpsiPtMtd_cent%s_rebin",cent_Title[icent]),xbins);
      hJpsiTrigEff[icent][1]->Divide(htmp);
      hJpsiTrigEff[icent][1]->SetMarkerStyle(20);
      hJpsiTrigEff[icent][1]->SetMarkerColor(2);
      hJpsiTrigEff[icent][1]->Draw("sames");

      hJpsiRespEff[icent][0] = (TH1F*)hRcJpsiPt[icent][3]->Clone(Form("JpsiRespEff_cent%s",cent_Title[icent]));
      hJpsiRespEff[icent][0]->Divide(hRcJpsiPt[icent][2]);
      c = draw1D(hJpsiRespEff[icent][0]);

      hJpsiRespEff[icent][1] = (TH1F*)hRcJpsiPt[icent][3]->Rebin(nbins,Form("JpsiRespEff_cent%s_rebin",cent_Title[icent]),xbins);
      htmp = (TH1F*)hRcJpsiPt[icent][2]->Rebin(nbins,Form("JpsiRespEff_cent%s_rebin",cent_Title[icent]),xbins);
      hJpsiRespEff[icent][1]->Divide(htmp);
      hJpsiRespEff[icent][1]->SetMarkerStyle(20);
      hJpsiRespEff[icent][1]->SetMarkerColor(2);
      hJpsiRespEff[icent][1]->Draw("sames");
    }
 
  TString outName = Form("Run14.AuAu200.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  TFile *fEff = TFile::Open(Form("Rootfiles/%s",outName.Data()),"update");

  for(int icent=0; icent<nCentBins; icent++)
    {
      for(int i=0; i<2; i++)
	{
	  hJpsiTrigEff[icent][i]->Write("",TObject::kOverwrite);
	  hJpsiRespEff[icent][i]->Write("",TObject::kOverwrite);
	}
    }
  
}

//================================================
void tuneResolution(const int icent, const bool savePlot)
{
  bool scan = 0;
  bool fit = 0;
  const int nScan = 40;

  TF1 *funcTrigEff = 0;
  TF1 *funcTrigEffCorr = 0;
  TFile *fEff = TFile::Open("Rootfiles/Run14.AuAu200.JpsiEff.root","update");
  funcTrigEff = (TF1*)fEff->Get(Form("MuonTrigEff_cent%s",cent_Title[icent]));
  funcTrigEffCorr = (TF1*)fEff->Get(Form("MuonTrigEffCorr_cent%s",cent_Title[icent]));

  TString name[3] = {"TPCreco","MTDreco","MTDtrig"};
  TFile *fin = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");
  //TH1F *hMcJpsiPt = (TH1F*)fin->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[icent]));
  TF1 *funcMcJpsiPt = (TF1*)fin->Get(Form("Fit_Jpsi_InvYield_AuAu200_Combined_cent%s",cent_Title[icent]));
  funcMcJpsiPt->SetNpx(1000);
  TH1F *hMcJpsiPt = (TH1F*)funcMcJpsiPt->GetHistogram();
  TH1F *hRcJpsiPt[3];
  TH1F *hRcJpsiMass[3];
  for(int i=0; i<3; i++)
    {
      hRcJpsiPt[i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
      hRcJpsiMass[i] = new TH1F(Form("hRcJpsiMass_%s_%s",name[i].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),4000,2,4);
    }
  if(!fit) smear(jpsiMass, icent, 5e7, 0, hMcJpsiPt, hRcJpsiPt[0], hRcJpsiMass[0], hRcJpsiPt[1], hRcJpsiMass[1], hRcJpsiPt[2], hRcJpsiMass[2], funcTrigEff, funcTrigEffCorr, 0);

  // real data
  TFile *fdata = TFile::Open("Rootfiles/Pico.Run14.AuAu200.jpsi.yield.root","read");
  TH1F *hRawCounts = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]));
  TH1F *hSignal = (TH1F*)fdata->Get(Form("Jpsi_Signal_cent%s_pt0-10",cent_Title[icent]));
  TF1 *fitSignal = (TF1*)fdata->Get(Form("FitJpsi_Signal_cent%s_pt0-10",cent_Title[icent]));

  // compare
  if(!fit)
    {
      hRawCounts->SetMarkerStyle(21);
      hRawCounts->Scale(1./hRawCounts->Integral());
      TCanvas *c = draw1D(hRawCounts,Form("p_{T} distribution of J/#psi in MTD"),kFALSE,kTRUE);
      for(int i=2; i<3; i++)
	{
	  hRcJpsiPt[i]->Sumw2();
	  hRcJpsiPt[i]->SetLineColor(2);
	  TH1F *hRebin = (TH1F*)hRcJpsiPt[i]->Rebin(hRawCounts->GetXaxis()->GetNbins(),Form("%s_rebin",hRawCounts->GetName()),hRawCounts->GetXaxis()->GetXbins()->GetArray());
	  hRebin->Scale(1./hRebin->Integral());
	  hRebin->Draw("sameHIST");
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_ToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiPt_ToyMcVsData_cent%s.png",run_type,cent_Title[icent]));
	}

      for(int i=2; i<3; i++)
	{
	  hRcJpsiMass[i]->Sumw2();
	  hRcJpsiMass[i]->SetLineColor(2);
	  hRcJpsiMass[i]->Scale(1./hRcJpsiMass[i]->Integral());
	  hRcJpsiMass[i]->GetXaxis()->SetRangeUser(2.8,3.4);
	  c = draw1D(hRcJpsiMass[i],Form("%s: mass distribution of J/#psi",name[i].Data()),kFALSE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/TBW_%s_JpsiMass_cent%s.pdf",run_type,name[i].Data(),cent_Title[icent]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/TBW_%s_JpsiMass_cent%s.png",run_type,name[i].Data(),cent_Title[icent]));
	    }
	}
    }  

  TH1F *hRcJpsiPtScan[nScan][3];
  TH1F *hRcJpsiMassScan[nScan][3];
  if(scan)
    {
      for(int i=0; i<nScan; i++)
	{
	  double sigma = i*0.001;
	  for(int j=0; j<3; j++)
	    {
	      hRcJpsiPtScan[i][j] = new TH1F(Form("hRcJpsiPt_%s_%s_scan%d",name[j].Data(),cent_Title[icent],i),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
	      hRcJpsiMassScan[i][j] = new TH1F(Form("hRcJpsiMass_%s_%s_scan%d",name[j].Data(),cent_Title[icent],i),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),4000,2,4);
	    }
	  smear(jpsiMass, icent, 1e7, sigma, hMcJpsiPt, hRcJpsiPtScan[i][0], hRcJpsiMassScan[i][0], hRcJpsiPtScan[i][1], hRcJpsiMassScan[i][1],hRcJpsiPtScan[i][2], hRcJpsiMassScan[i][2], funcTrigEff, funcTrigEffCorr, 0);
	}

      TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","recreate");
      for(int i=0; i<nScan; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hRcJpsiPtScan[i][j]->Write();
	      hRcJpsiMassScan[i][j]->Write();
	    }
	}
    }

  if(fit)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.TrkResScan.root","read");
      TH1F *hChi2 = new TH1F("hChi2","#Chi^{2}/ndf vs smear factor;#sigma (%);Chi^{2}/ndf",40,-0.05,3.95);
      TF1 *func[nScan];
      for(int i=0; i<nScan; i++)
	{
	  // fit
	  hRcJpsiMassScan[i][2] = (TH1F*)fout->Get(Form("hRcJpsiMass_%s_%s_scan%d",name[2].Data(),cent_Title[icent],i));
	  hJpsiLineShape = (TH1F*)hRcJpsiMassScan[i][2]->Clone(Form("hJpsiLineShape_scan%d",i));
	  hJpsiLineShape->Rebin(10);
	  //draw1D(hJpsiLineShape);
	  func[i] = new TF1(Form("CombinedFit_cent%s_scan%d",cent_Title[icent],i),combinedFit,2.5,4,5);
	  for(int j=1; j<5; j++)
	    func[i]->FixParameter(j,fitSignal->GetParameter(j+2));
	  
	  TH1F *hdata = (TH1F*)hSignal->Clone(Form("%s_scan%d",hSignal->GetName(),i));
	  hdata->Fit(func[i],"IR0");
	  TCanvas *c = draw1D(hdata);
	  func[i]->SetNpx(1000);
	  TH1F *hfit = (TH1F*)func[i]->GetHistogram();
	  hfit->SetName(Form("h%s",func[i]->GetName()));
	  hfit->SetLineColor(2);
	  hfit->Draw("sames");
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Scan%d_FitData_cent%s.pdf",run_type,i,cent_Title[icent]));
	    }

	  hChi2->SetBinContent(i+1, func[i]->GetChisquare()/func[i]->GetNDF());
	  hChi2->SetMarkerStyle(21);
	  hChi2->SetBinError(i+1, 0);
	}
      TCanvas *c = draw1D(hChi2);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ScanRes_Chi2_cent%s.pdf",run_type,cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/ScanRes_Chi2_cent%s.png",run_type,cent_Title[icent]));
	}
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
void flatPt(const int icent, const bool savePlot)
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

  TH1F *hMcJpsiPt   = new TH1F(Form("hMcJpsiPt_%s",cent_Title[icent]),Form("p_{T} distribution of input MC J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
  TH1F *hRcJpsiPt[2];
  TH1F *hRcJpsiMass[2];
  for(int i=0; i<2; i++)
    {
      hRcJpsiPt[i] = new TH1F(Form("hRcJpsiPt_%s_%s",name[i+1].Data(),cent_Title[icent]),Form("p_{T} distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c)",cent_Name[icent]),110,0,11);
      hRcJpsiMass[i] = new TH1F(Form("hRcJpsiMass_%s_%s",name[i+1].Data(),cent_Title[icent]),Form("Mass distribution of reconstructed J/#psi (%s%%);mass (GeV/c^{2})",cent_Name[icent]),1400,0,14);
    }

  smear(jpsiMass, icent, 1e7, 0, hJpsiPt[0], hRcJpsiPt[0], hRcJpsiMass[0], hRcJpsiPt[1], hRcJpsiMass[1], 0, 0, 0, 0, 0);

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
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiPt_cent%s.pdf",run_type,name[i+1].Data(),cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiPt_cent%s.png",run_type,name[i+1].Data(),cent_Title[icent]));
	}

      hJpsiMass[i+1]->Sumw2();
      hJpsiMass[i+1]->Scale(1./hJpsiMass[i+1]->Integral());
      c = draw1D(hJpsiMass[i+1],Form("%s: mass distribution of J/#psi",name[i+1].Data()),kTRUE,kTRUE);
      hRcJpsiMass[i]->Sumw2();
      hRcJpsiMass[i]->Scale(1./hRcJpsiMass[i]->Integral());
      hRcJpsiMass[i]->SetLineColor(2);
      hRcJpsiMass[i]->Draw("sames HIST");
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMass_cent%s_log.pdf",run_type,name[i+1].Data(),cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMass_cent%s_log.png",run_type,name[i+1].Data(),cent_Title[icent]));
	}
      TH1F *htmp = (TH1F*)hJpsiMass[i+1]->Clone(Form("%s_linear",hJpsiMass[i+1]->GetName()));
      c = draw1D(htmp,Form("%s: mass distribution of J/#psi",name[i+1].Data()),kFALSE,kTRUE);
      hRcJpsiMass[i]->Draw("sames HIST");
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMass_cent%s.pdf",run_type,name[i+1].Data(),cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMass_cent%s.png",run_type,name[i+1].Data(),cent_Title[icent]));
	}

      TH1F *hMassRatio = (TH1F*)hRcJpsiMass[i]->Clone(Form("hMassRatio_%s_%s",name[i+1].Data(),cent_Title[icent]));
      hMassRatio->Divide(hJpsiMass[i+1]);
      hMassRatio->GetXaxis()->SetRangeUser(2.8,3.4);
      hMassRatio->GetYaxis()->SetRangeUser(0,2);
      hMassRatio->SetMarkerStyle(20);
      hMassRatio->SetMarkerColor(2);
      c = draw1D(hMassRatio,Form("%s: ratio of mass distribution of J/#psi;M_{#mu#mu} (GeV/c^{2});ToyModel/Embed",name[i+1].Data()),kFALSE,kTRUE);
      TLine *line = GetLine(2.8,1,3.4,1,1);
      line->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMassRatio_cent%s.pdf",run_type,name[i+1].Data(),cent_Title[icent]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/FlatPt_%s_JpsiMassRatio_cent%s.png",run_type,name[i+1].Data(),cent_Title[icent]));
	}
    }
}

//================================================
void smear(const double masss, const int icent, const int nExpr, const double sigma, TH1F *hInputPt, 
	   TH1F *hRcJpsiPtTpc, TH1F *hRcJpsiMassTpc, TH1F *hRcJpsiPtMtd, TH1F *hRcJpsiMassMtd, 
	   TH1F *hRcJpsiPtMtdTrig, TH1F *hRcJpsiMassMtdTrig,
	   TF1 *funcTrigEff, TF1 *funcTrigEffCorr,
	   const bool debug,
	   TH1F *hRcJpsiPtMtdResp, TH1F *hRcJpsiMassMtdResp)
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

      // momentum resolution      
      int mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt1)-1;
      if(pt1>10) mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;
      int mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(pt2)-1;
      if(pt2>10) mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindBin(9)-1;
      double rc_pt1 = (1-hTrkResBin[icent][mom_index1]->GetRandom() + myRandom->Gaus(0,sigma)) * pt1;
      double rc_pt2 = (1-hTrkResBin[icent][mom_index2]->GetRandom() + myRandom->Gaus(0,sigma)) * pt2;

      if(debug) printf("rc daug 1: pt = %3.2f\n",rc_pt1);
      if(debug) printf("rc daug 2: pt = %3.2f\n",rc_pt2);
      if(rc_pt1<pt1_cut || rc_pt2<pt2_cut) continue;
      TLorentzVector rc_daughter1, rc_daughter2;
      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;
      hRcJpsiPtTpc->Fill(rc_parent.Pt());
      hRcJpsiMassTpc->Fill(rc_parent.M());

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
      hRcJpsiMassMtd->Fill(rc_parent.M());

      // MTD efficiency
      if(!funcTrigEff || !funcTrigEffCorr) continue;
      double trig_eff_1 = funcTrigEff->Eval(rc_pt1) * funcTrigEffCorr->Eval(rc_pt1);
      if(debug) printf("MTD trig eff 1 = %3.2f\n",trig_eff_1);
      if( myRandom->Uniform(0., 1.) >  trig_eff_1) continue;

      double trig_eff_2 = funcTrigEff->Eval(rc_pt2) * funcTrigEffCorr->Eval(rc_pt2);
      if(debug) printf("MTD trig eff 2 = %3.2f\n",trig_eff_2);
      if( myRandom->Uniform(0., 1.) >  trig_eff_2) continue;
      if(hRcJpsiPtMtdTrig) hRcJpsiPtMtdTrig->Fill(rc_parent.Pt());
      if(hRcJpsiMassMtdTrig) hRcJpsiMassMtdTrig->Fill(rc_parent.M());

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
      if(hRcJpsiPtMtdResp) hRcJpsiPtMtdResp->Fill(rc_parent.Pt());
      if(hRcJpsiMassMtdResp) hRcJpsiMassMtdResp->Fill(rc_parent.M());
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
