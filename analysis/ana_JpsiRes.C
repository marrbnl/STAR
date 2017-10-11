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

#define YEAR 2014

#if (YEAR==2013)
const char *run_config = "VtxCut.";
const char *run_type = "Run13_pp500";
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
#elif (YEAR==2014)
const char *run_config = "";
const char *run_type = "Run14_AuAu200";
const double pt1_cut = 1.5, pt2_cut = 1.3;
const double low_mass = 2.9;
const double high_mass = 3.3;
const int nCentBins = 11;
const char *cent_Name[nCentBins] = {"0-80","0-10","10-20","20-30","30-40","40-50","50-60","60-80","0-20","20-40","40-60"};
const char *cent_Title[nCentBins] = {"0080","0010","1020","2030","3040","4050","5060","6080","0020","2040","4060"};
const int nPtBins = 10;
const double ptBins_low[nPtBins]  = {0,0,1,2,3,4,5,6,8,10};
const double ptBins_high[nPtBins] = {15,1,2,3,4,5,6,8,10,15};
const char *pt_Name[nPtBins] = {"0-15","0-1","1-2","2-3","3-4","4-5","5-6","6-8","8-10","10-15"};
#elif (YEAR==2015)
const char *run_config = "";
const char *run_type = "Run15_pAu200";
const double pt1_cut = 1.0, pt2_cut = 1.0;
const double low_mass = 2.9, high_mass = 3.3;
const int nPtBins = 9;
const double ptBins_low[nPtBins]  = {0,0,1,2,3,4,5,6,8};
const double ptBins_high[nPtBins] = {10,1,2,3,4,5,6,8,10};
const char *pt_Name[nPtBins] = {"0-10","0-1","1-2","2-3","3-4","4-5","5-6","6-8","6-10"};
const int nCentBins = 1;
const int centBins_low[nCentBins]  = {0};
const int centBins_high[nCentBins] = {-1};
const char *cent_Name[nCentBins] = {"0-80"};
const char *cent_Title[nCentBins] = {"0080"};
#elif (YEAR==2016)
const char *run_config = "";
const char *run_type = "Run16_AuAu200";
const double pt1_cut = 1.5, pt2_cut = 1.3;
const double low_mass = 2.9, high_mass = 3.3;
const int nPtBins = 5;
const double ptBins_low[nPtBins]  = {0,0,2,4,6};
const double ptBins_high[nPtBins] = {10,2,4,6,10};
const char *pt_Name[nPtBins] = {"0-10","0-2","2-4","4-6","6-10"};
const int nCentBins = 1;
const int centBins_low[nCentBins]  = {0};
const int centBins_high[nCentBins] = {-1};
const char *cent_Name[nCentBins] = {"0-80"};
const char *cent_Title[nCentBins] = {"0080"};
#endif

const double pi = 3.1415926;
const double jpsiMass = 3.096;
const double muMass = 0.1057;
const int year = YEAR;
const int kNPtBins = 40;
const double kLowPtBound = 0;
const double kHighPtBound = 20;
const int kNMassBins = 1400;
const double kLowMassBound = 0;
const double kHighMassBound = 14;

TRandom3 *myRandom;
TH2F *hTrkResVsPt[nCentBins];
TH1F *hTrkResBin[nCentBins][400];
TF1 *funcTrkRes[nCentBins];
TH1F *hDeltaPt;

void tuneResolution(const int icent, const bool savePlot);
void smear(const double mass, const int icent, const int nExpr, const double shift, const double sigma, TF1 *hInputPt, TH2F *hRcJpsiMassVsPt,  const bool isWeight = 0, const bool debug = 0);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);

TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size = 0.04, const Int_t font = 42);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);
TCanvas *draw2D(TH2 *h, const TString hTitle = "" , const Float_t size = 0.04, const Bool_t logz = kTRUE, const char *drawOption = "colz");
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color=2, Width_t width=2, Style_t style=2);

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
    
  // tracking resolution for single muons
  TFile *fRes = 0x0;
  if(year==2014) fRes = TFile::Open(Form("Rootfiles/%s.EmbTrkEff.root",run_type),"read");
  else if(year==2015) fRes = TFile::Open(Form("Rootfiles/%s.PtResolution.root",run_type),"read");
  else if(year==2016) fRes = TFile::Open("Rootfiles/Upsilon1S.ptrk.Embed.root","read");
  if(!fRes)
    {
      printf("[e] No track resolution file avaialbe!\n");
    }
  else
    {
      if(year==2014)
	{
	  hDeltaPt = (TH1F*)fRes->Get("hdpTOverPt_0080");
	}
      for(int icent=0; icent<nCentBins; icent++)
	{
	  if(year==2014)
	    {
	      hTrkResVsPt[icent] = (TH2F*)fRes->Get(Form("pTrkRes_vs_TruePt_%s",cent_Title[icent]));
	      funcTrkRes[icent] = (TF1*)fRes->Get(Form("pTrkResFit_%s",cent_Title[icent]));
	    }
	  else if(year==2015)
	    {
	      hTrkResVsPt[icent] = (TH2F*)fRes->Get(Form("pTrkRes_vs_TruePt_%s",cent_Title[icent]));
	    }
	  else if(year==2016)
	    {
	      hTrkResVsPt[icent] = (TH2F*)fRes->Get(Form("McTrkPtReso_cent%s",cent_Title[icent]));
	    }
	  int nHistos = hTrkResVsPt[icent]->GetNbinsX();
	  for(int i=0; i<nHistos; i++)
	    {
	      hTrkResBin[icent][i] = (TH1F*)hTrkResVsPt[icent]->ProjectionY(Form("hTrkRes_Bin%d_cent%s",i+1,cent_Title[icent]),i+1,i+1);
	    }
	}
    }

  tuneResolution(0,1);
}

//================================================
void tuneResolution(const int icent, const bool savePlot)
{
  const int anaType = 1; // 0 - scan; 1 - determine smear & shift; 2 - generate final smear
  const int nShiftScan = 1;
  const double shiftStep = 0.0005;
  const int nSmearScan = 20;
  const double smearStep = 0.0005;
  const TString outName = Form("Rootfiles/%s.TrkResScan.%sroot",run_type,run_config);
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // initilize pointers
  TCanvas *c = 0x0, *c1 = 0x0;
  TPaveText *t1 = 0x0;
  TFile *fscan = 0x0;
  TLegend *leg = 0x0;
  
  printf("+++++ Enter tuneResolution() +++\n");

  // input spectrum shape from data
  TF1 *funcInputJpsi = 0x0;
  TFile *fdata = TFile::Open(Form("Rootfiles/%s.TrkResScan.input.root",run_type),"read");
  if(year==2014)
    {
      funcInputJpsi = (TF1*)fdata->Get("Fit_Jpsi_FitYield_cent0080_weight");
    }

  // scan smear and shift
  TH2F *hRcJpsiMassScan[nShiftScan][nSmearScan];
  if(anaType == 0)
    {
      for(int i=0; i<nShiftScan; i++)
	{
	  for(int k=0; k<nSmearScan; k++)
	    {
	      hRcJpsiMassScan[i][k] = new TH2F(Form("hRcJpsiMassVsPt_%s_scan%d_%d",cent_Title[icent],i,k),Form("Mass distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c);mass (GeV/c^{2})",cent_Name[icent]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
	      if(i>0 && k>0) continue;
	      double sigma = k*smearStep;
	      double shift = 0;
	      if(i==0)        shift = 0;
	      else if(i%2==1) shift = (i+1)/2*shiftStep;
	      else            shift = (i+1)/2*shiftStep*(-1);
	      printf("[i] Scan (%d,%d) = (%5.4f, %5.4f)\n",i,k,shift,sigma);
	      smear(jpsiMass, icent, 5e6, shift, sigma, funcInputJpsi, hRcJpsiMassScan[i][k], 0, 0);
	    }
	}
      TFile *fout = TFile::Open(outName.Data(),"recreate");
      for(int i=0; i<nShiftScan; i++){
	for(int k=0; k<nSmearScan; k++){
	    if(i>0 && k>0) continue;
	    hRcJpsiMassScan[i][k]->Write();
	}
      }
    }


  // determine smear and shift by comparing to real data
  if(anaType == 1)
    {
      fscan = TFile::Open(outName.Data(),"update");

      // check with embedding
      if(year==2013 || year==2014)
	{
	  TFile *fembed = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");

	  TH2F *hEmbedJpsi[2];
	  hEmbedJpsi[0] = (TH2F*)fembed->Get(Form("hJpsiMassVsPt_Tpc_cent%s",cent_Title[icent]));
	  hEmbedJpsi[1] = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_scan0_8",cent_Title[icent]));

	  TObjArray embSlices[2];
	  TH1F *hEmbMean[2];
	  TH1F *hEmbSigma[2];
	  TF1 *f1 = new TF1("f1","gaus",3.0,3.2);
	  for(int i=0; i<2; i++)
	    {
	      c = draw2D(hEmbedJpsi[i]);
	      hEmbedJpsi[i]->FitSlicesY(f1, 0, -1, 0, "QNR", &embSlices[i]);
	      hEmbMean[i] = (TH1F*)((TH1F*)embSlices[i][1])->Clone(Form("hEmbMean_%d",i));
	      hEmbMean[i]->SetMarkerStyle(21+i*4);
	      hEmbMean[i]->SetMarkerColor(i+1);
	      hEmbMean[i]->SetLineColor(i+1);
	      hEmbMean[i]->GetXaxis()->SetRangeUser(0,15);
	      hEmbMean[i]->GetYaxis()->SetRangeUser(3.08,3.12);
	      hEmbMean[i]->SetTitle(";p_{T} (GeV/c);Mean");

	      hEmbSigma[i] = (TH1F*)((TH1F*)embSlices[i][2])->Clone(Form("hEmbSigma_%d",i));
	      hEmbSigma[i]->SetMarkerStyle(21+i*4);
	      hEmbSigma[i]->SetMarkerColor(i+1);
	      hEmbSigma[i]->SetLineColor(i+1);
	      hEmbSigma[i]->GetXaxis()->SetRangeUser(0,15);
	      if(year==2013) hEmbSigma[i]->GetYaxis()->SetRangeUser(0.02,0.06);
	      if(year==2014) hEmbSigma[i]->GetYaxis()->SetRangeUser(0.02,0.1);
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
	  leg->SetHeader(run_type);
	  leg->AddEntry(hEmbMean[0],"Embedding","PL");
	  leg->AddEntry(hEmbMean[1],"ToyMC","PL");
	  leg->Draw();
	  if(savePlot)
	    c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ToyMcVsEmbed_cent%s.pdf",run_type,cent_Title[icent]));
	  c1 = new TCanvas("embed_sigma", "embed_sigma", 800, 600);
	  hEmbSigma[0]->Draw();
	  hEmbSigma[1]->Draw("sames");
	  t1 = GetTitleText("Width of J/#Psi mass peak");
	  t1->Draw();
	  leg->Draw();
	  if(savePlot)
	    c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ToyMcVsEmbed_cent%s.pdf",run_type,cent_Title[icent]));
	}

      // check with real data
      gStyle->SetOptFit(0);
      const int makePlot = 0;
      const char *sysName[3] = {"def","min","max"};
      const char *sysTitle[3] = {"default values","1#sigma lower limit","1#sigma upper limit"};
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
	  TH2F *htmp = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_scan%d_%d",cent_Title[icent],i,0));
	  TH2F *h2 = (TH2F*)htmp->Clone(Form("%s_clone",htmp->GetName()));
	  hShiftMean[i] = (TH1F*)h2->ProjectionX(Form("hShiftMean_%d",i));
	  hShiftMean[i]->Reset();
	  TCanvas *cShiftFit = 0x0;
	  if(makePlot)
	    {
	      cShiftFit = new TCanvas(Form("fit_shift_%d",i), Form("fit_shift_%d",i), 1200, 800);
	      cShiftFit->Divide(7,6);
	    }
	  for(int j=0; j<hShiftMean[i]->GetNbinsX(); j++)
	    {
	      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hShiftMass_%d_%d",i,j),j+1,j+1);
	      if(h1->GetEntries()<100) continue;
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
	  hShiftMean[i]->GetXaxis()->SetRangeUser(ptBins_low[0],ptBins_high[0]);
	  hShiftMean[i]->GetYaxis()->SetRangeUser(2.9,3.15);
	  if(year==2015) hShiftMean[i]->GetYaxis()->SetRangeUser(3.05,3.15);
	  if(year==2016) hShiftMean[i]->GetYaxis()->SetRangeUser(3,3.2);
	  hShiftMean[i]->SetMarkerStyle(20+i);
	  hShiftMean[i]->SetTitle(";p_{T} (GeV/c);Mean");
	  cShift->cd();
	  if(i==0) hShiftMean[i]->Draw();
	  else     hShiftMean[i]->Draw("sames");
	  hFitShiftMean[i] = new TF1(Form("hFitShiftMean_%d",i),"pol4",ptBins_low[0],ptBins_high[0]);
	  hShiftMean[i]->Fit(hFitShiftMean[i],"IRQ0");
	  hFitShiftMean[i]->SetLineColor(4);
	  hFitShiftMean[i]->SetLineStyle(2);
	  hFitShiftMean[i]->Draw("sames");
	}
      TH1F *hDataMean = (TH1F*)fdata->Get(Form("Jpsi_FitMean_cent%s_weight",cent_Title[icent]));
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
	cShift->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_ScanToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));

      
      printf("[i] Scan shift\n");
      TH1F *hShiftChi2[3];
      TF1 *hFitShiftChi2[3];
      for(int i=0; i<3; i++)
	{
	  hShiftChi2[i] = new TH1F(Form("hShiftChi2_%s",sysName[i]),";shift; #chi^{2}", nShiftScan, -1*(nShiftScan/2)*shiftStep-shiftStep/2, (nShiftScan/2)*shiftStep+shiftStep/2);
	  for(int j=0; j<nShiftScan; j++)
	    {
	      double shift = 0;
	      if(j==0)        shift = 0;
	      else if(j%2==1) shift = (j+1)/2*shiftStep;
	      else            shift = (j+1)/2*shiftStep*(-1);
	      int bin = hShiftChi2[i]->FindFixBin(shift);
	      double chi2 = 0;
	      for(int ibin=1; ibin<=hDataMean->GetNbinsX(); ibin++)
		{
		  double value = hDataMean->GetBinContent(ibin);
		  double error = hDataMean->GetBinError(ibin);
		  if(i==1) value -= 1.5*error;
		  if(i==2) value += 1.5*error;
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
	    c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiMean_Chi2_%s_cent%s.pdf",run_type,sysName[i],cent_Title[icent]));
	}

      // determine smear
      TCanvas *cSmear = new TCanvas("scan_smear", "scan_smear", 800, 600);
      TObjArray smearSlices[nSmearScan];
      TH1F *hSmearSigma[nSmearScan];
      TF1 *hFitSmearSigma[nSmearScan];
      for(int i=0; i<nSmearScan; i++)
	{
	  TH2F *htmp = (TH2F*)fscan->Get(Form("hRcJpsiMassVsPt_%s_scan%d_%d",cent_Title[icent],0,i));
	  TH2F *h2 = (TH2F*)htmp->Clone(Form("%s_clone2",htmp->GetName()));
	  hSmearSigma[i] = (TH1F*)h2->ProjectionX(Form("hSmearSigma_%d",i));
	  hSmearSigma[i]->Reset();
	  TCanvas *cSmearFit = 0x0;
	  if(makePlot)
	    {
	      cSmearFit = new TCanvas(Form("fit_smear_%d",i), Form("fit_smear_%d",i), 1200, 800);
	      cSmearFit->Divide(7,6);
	    }
	  for(int j=0; j<hSmearSigma[i]->GetNbinsX(); j++)
	    {
	      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hSmearMass_%d_%d",i,j),j+1,j+1);
	      if(h1->GetEntries()<100) continue;
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
	  hSmearSigma[i]->GetXaxis()->SetRangeUser(ptBins_low[0],ptBins_high[0]);
	  hSmearSigma[i]->GetYaxis()->SetRangeUser(0.02,0.2);
	  hSmearSigma[i]->SetMarkerStyle(20+i);
	  hSmearSigma[i]->SetTitle(";p_{T} (GeV/c);#sigma");
	  cSmear->cd();
	  if(i==0) hSmearSigma[i]->Draw();
	  else     hSmearSigma[i]->Draw("sames");
	  hFitSmearSigma[i] = new TF1(Form("hFitSmearSigma_%d",i),"pol4",ptBins_low[0],ptBins_high[0]);
	  hSmearSigma[i]->Fit(hFitSmearSigma[i],"IRQ0");
	  hFitSmearSigma[i]->SetLineColor(4);
	  hFitSmearSigma[i]->SetLineStyle(2);
	  hFitSmearSigma[i]->Draw("sames");
	}
      TH1F *hDataSigma = (TH1F*)fdata->Get(Form("Jpsi_FitSigma_cent%s_weight",cent_Title[icent]));
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
	cSmear->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_ScanToyMcVsData_cent%s.pdf",run_type,cent_Title[icent]));

      printf("[i] Scan smear\n");
      TH1F *hSmearChi2[3];
      TF1 *hFitSmearChi2[3];
      for(int i=0; i<3; i++)
	{
	  hSmearChi2[i] = new TH1F(Form("hSmearChi2_%s",sysName[i]),";smear; #chi^{2}", nSmearScan, 0-smearStep/2, (nSmearScan-1)*smearStep+smearStep/2);
	  for(int j=0; j<nSmearScan; j++)
	    {
	      int bin = j + 1;
	      double chi2 = 0;
	      for(int ibin=1; ibin<=hDataSigma->GetNbinsX(); ibin++)
		{
		  double value = hDataSigma->GetBinContent(ibin);
		  double error = hDataSigma->GetBinError(ibin);
		  if(i==1) value -= 1.0*error;
		  if(i==2) value += 1.0*error;
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
	    c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/JpsiSigma_Chi2_%s_cent%s.pdf",run_type,sysName[i],cent_Title[icent]));
	}
      fscan->cd();
      hFinalShift->Write("",TObject::kOverwrite);
      hFinalSmear->Write("",TObject::kOverwrite);
    }

  // final smear and shift
  const char *sysName_scan[3] = {"def","min","max"};
  if(anaType == 2)
    {
      TH2F *hRcJpsiMass[nCentBins][3];
      TFile *fout = TFile::Open(outName.Data(),"update");
      TH1F *hFinalShift = (TH1F*)fout->Get(Form("hFinalShift_cent%s",cent_Title[icent]));
      TH1F *hFinalSmear = (TH1F*)fout->Get(Form("hFinalSmear_cent%s",cent_Title[icent]));

      TH1F *hJpsiEff[nCentBins][3][2];
      for(int k=0; k<nCentBins; k++)
	{
	  for(int i=0; i<3; i++)
	    {
	      hRcJpsiMass[k][i] = new TH2F(Form("hRcJpsiMass_%s_final_%s",cent_Title[k],sysName_scan[i]),Form("Mass distribution of reconstructed J/#psi (%s%%);p_{T} (GeV/c);mass (GeV/c^{2})",cent_Name[k]),kNPtBins, kLowPtBound, kHighPtBound, kNMassBins, kLowMassBound, kHighMassBound);
	      double shiftvalue = hFinalShift->GetBinContent(i+1);
	      double smearvalue = hFinalSmear->GetBinContent(i+1);
	      if(year==2014 || year==2015 || year==2016) shiftvalue = 0;
	      printf("[i] Start final scanning: %s = (%4.4f,%4.4f)\n",sysName_scan[i],shiftvalue,smearvalue);
	      smear(jpsiMass, k, 1e7, shiftvalue, smearvalue, funcInputJpsi, hRcJpsiMass[k][i], 1, 0);
	    }

      
	  TH1F *hAll = (TH1F*)hRcJpsiMass[k][0]->ProjectionX(Form("JpsiPt_all_cent%s",cent_Title[k]));
	  TH1F *hAllrebin = (TH1F*)hAll->Rebin(nbins,Form("%s_rebin",hAll->GetName()),xbins);
	  for(int i=0; i<3; i++)
	    {
	      hJpsiEff[k][i][0] = (TH1F*)hRcJpsiMass[k][i]->ProjectionX(Form("JpsiSmearEff_%s_cent%s",sysName_scan[i],cent_Title[k]));
	      hJpsiEff[k][i][1] = (TH1F*)hJpsiEff[k][i][0]->Rebin(nbins,Form("JpsiSmearEff_%s_cent%s_rebin",sysName_scan[i],cent_Title[k]),xbins);
	      hJpsiEff[k][i][0]->Divide(hAll);
	      hJpsiEff[k][i][1]->Divide(hAllrebin);
	    }

	  for(int i=0; i<3; i++)
	    {
	      hRcJpsiMass[k][i]->Write("",TObject::kOverwrite);
	      hJpsiEff[k][i][0]->Write("",TObject::kOverwrite);
	      hJpsiEff[k][i][1]->Write("",TObject::kOverwrite);
	    }
	}
    }

  // check final smear and shift
  if(anaType == 3)
    {
      gStyle->SetOptFit(0);
      fscan = TFile::Open(outName.Data(),"update");
      TH2F *hMassVsPtEmbed[3];
      TObjArray embedSlice[3];
      TH1F *hEmbSmearMean[3];
      TH1F *hEmbSmearSigma[3];

      for(int i=0; i<3; i++)
	{
	  hMassVsPtEmbed[i] = (TH2F*)fscan->Get(Form("hRcJpsiMass_%s_final_%s",cent_Title[icent],sysName_scan[i]));
	  hMassVsPtEmbed[i]->FitSlicesY(0, 0, -1, 0, "QNR", &embedSlice[i]);
	  hEmbSmearMean[i]  =  (TH1F*)((TH1F*)embedSlice[i][1])->Clone(Form("hEmbSmearMean_%s",sysName_scan[i]));
	  hEmbSmearSigma[i] =  (TH1F*)((TH1F*)embedSlice[i][2])->Clone(Form("hEmbSmearSigma_%s",sysName_scan[i])); 
	}     
  
      // mean
      TH1F *hDataMean = (TH1F*)fdata->Get(Form("Jpsi_FitMean_cent%s_weight",cent_Title[icent]));
      hDataMean->SetMarkerStyle(21);
      hDataMean->SetMarkerColor(2);
      hDataMean->SetLineColor(2);
      hDataMean->SetTitle("Mean of J/#Psi peak");
      hDataMean->GetYaxis()->SetRangeUser(3.04,3.16);
      TCanvas *c = draw1D(hDataMean);
      hEmbSmearMean[0]->SetMarkerStyle(25);
      hEmbSmearMean[0]->Draw("sames");
      TLegend *leg = new TLegend(0.18,0.72,0.42,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hDataMean,"Real data","PL");
      leg->AddEntry(hEmbSmearMean[0],"Smeared embedding","PL");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/SmearEmbedVsData_JpsiMean_cent%s.pdf",run_type,cent_Title[icent]));

      // sigma
      TH1F *hDataSigma = (TH1F*)fdata->Get(Form("Jpsi_FitSigma_cent%s_weight",cent_Title[icent]));
      hDataSigma->SetMarkerStyle(21);
      hDataSigma->SetMarkerColor(2);
      hDataSigma->SetLineColor(2);
      hDataSigma->SetTitle("Width of J/#Psi peak");
      hDataSigma->GetYaxis()->SetRangeUser(0.02,0.14);
      c = draw1D(hDataSigma);
      hEmbSmearSigma[0]->SetMarkerStyle(25);
      hEmbSmearSigma[0]->Draw("sames");
      TF1* hFitSmearTmp = 0x0;
      for(int i=1; i<3; i++)
	{
	  hFitSmearTmp = new TF1(Form("hFitSmearTmp_%d",i),"pol4",ptBins_low[0],ptBins_high[0]);
	  hEmbSmearSigma[i]->Fit(hFitSmearTmp,"IRQ0");
	  hFitSmearTmp->SetLineColor(6);
	  hFitSmearTmp->SetLineStyle(2);
	  hFitSmearTmp->Draw("sames");
	}
      leg->Draw();
      TLegend *leg1 = new TLegend(0.18,0.68,0.42,0.72);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.04);
      leg1->AddEntry(hFitSmearTmp,"Uncertainty","L");
      leg1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/SmearEmbedVsData_JpsiSigma_cent%s.pdf",run_type,cent_Title[icent]));

      // signal shape
      TH1F *hJpsiShape[2][nPtBins-1];
      c = new TCanvas("CompareSigShape","CompareSigShape",1100,700);
      if(year==2014) c->Divide(5,2);
      for(int bin=0; bin<nPtBins-1; bin++)
	{
	  // data
	  TH1F *hSignal = (TH1F*)fdata->Get(Form("Jpsi_Signal_pt%s_cent%s_weight",pt_Name[bin+1],cent_Title[icent]));
	  TF1 *fitSignal = (TF1*)fdata->Get(Form("Jpsi_FitSig_pt%s_cent%s_weight",pt_Name[bin+1],cent_Title[icent]));
	  TF1 *funcBkg = new TF1(Form("funcBkg_pt%s",pt_Name[bin+1]),"pol3",2.5,4);
	  TF1 *funcSignal = new TF1(Form("funcSignal_pt%s",pt_Name[bin+1]), "gaus", 2.7, 3.5);
	  for(int i=0; i<4; i++)
	    funcBkg->SetParameter(i, fitSignal->GetParameter(i+3));
	  for(int i=0; i<3; i++)
	    funcSignal->SetParameter(i, fitSignal->GetParameter(i));

	  if(year==2015)
	    {
	      funcBkg    = (TF1*)fdata->Get("mfBkgCent0_Pt0");
	      funcSignal = (TF1*)fdata->Get("mfSigCent0_Pt0");
	    }

	  TAxis *axis = hSignal->GetXaxis();
	  hJpsiShape[0][bin] = (TH1F*)hSignal->Clone(Form("hSignalSub_pt%s",pt_Name[bin+1]));
	  hJpsiShape[0][bin]->Reset();
	  hJpsiShape[0][bin]->SetTitle("Compare J/#psi signal shape;M_{#mu#mu} (GeV/c^{2})");
	  for(int ibin=1; ibin<=hJpsiShape[0][bin]->GetNbinsX(); ibin++)
	    {
	      double value = hSignal->GetBinContent(ibin);
	      double error = hSignal->GetBinError(ibin);
	      double bkg = funcBkg->Integral(hSignal->GetXaxis()->GetBinLowEdge(ibin),hSignal->GetXaxis()->GetBinUpEdge(ibin))/hSignal->GetBinWidth(ibin);
	      //bkg = 0;
	      hJpsiShape[0][bin]->SetBinContent(ibin, value - bkg);
	      hJpsiShape[0][bin]->SetBinError(ibin, sqrt(error*error+fabs(bkg)));
	    }
	  hJpsiShape[0][bin]->SetMarkerStyle(21);
	  hJpsiShape[0][bin]->GetXaxis()->SetRangeUser(2.5,4);
	  if(year==2016) hJpsiShape[0][bin]->GetXaxis()->SetRangeUser(2.7,3.5);
	  double datascale = hJpsiShape[0][bin]->Integral(hJpsiShape[0][bin]->FindFixBin(2.9+1e-4),hJpsiShape[0][bin]->FindFixBin(3.3-1e-4),"width");
	  hJpsiShape[0][bin]->Scale(1./datascale);
	  hJpsiShape[0][bin]->SetLineColor(2);
	  hJpsiShape[0][bin]->SetMarkerColor(2);
	  hJpsiShape[0][bin]->SetMaximum(1.5*hJpsiShape[0][bin]->GetMaximum());
	  hJpsiShape[0][bin]->GetXaxis()->SetRangeUser(2.7,3.5);
	  hJpsiShape[0][bin]->SetTitle("");
	  c->cd(bin+1);
	  hJpsiShape[0][bin]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[bin+1],ptBins_high[bin+1]),0.07);
	  t1->Draw();

	  // embedding
	  hJpsiShape[1][bin] = (TH1F*)hMassVsPtEmbed[0]->ProjectionY(Form("Embed_JpsiShape_Pt%d",bin));
	  hJpsiShape[1][bin]->Reset();
	  int low_bin = hMassVsPtEmbed[0]->GetXaxis()->FindFixBin(ptBins_low[bin+1]+1e-4);
	  int high_bin = hMassVsPtEmbed[0]->GetXaxis()->FindFixBin(ptBins_high[bin+1]-1e-4);
	  for(int ibin=low_bin; ibin<=high_bin; ibin++)
	    {
	      TH1F *htmp = (TH1F*)hMassVsPtEmbed[0]->ProjectionY(Form("Embed_tmp_bin%d",ibin),ibin,ibin);
	      double pt = hMassVsPtEmbed[0]->GetXaxis()->GetBinCenter(ibin);
	      hJpsiShape[1][bin]->Add(htmp, funcInputJpsi->Eval(pt));
	    }
	  double scale_factor = hJpsiShape[0][bin]->GetBinContent(hJpsiShape[0][bin]->FindFixBin(3.095))/hJpsiShape[1][bin]->GetBinContent(hJpsiShape[1][bin]->FindFixBin(3.095));
	  hJpsiShape[1][bin]->SetMarkerStyle(25);
	  hJpsiShape[1][bin]->Scale(scale_factor);
	  hJpsiShape[1][bin]->Draw("samesHIST");
	}
      c->cd(10);
      leg = new TLegend(0.1,0.5,0.42,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.07);
      leg->SetHeader(run_type);
      leg->AddEntry(hJpsiShape[0][0],"Data","PL");
      leg->AddEntry(hJpsiShape[1][0],"Smeared embedding","L");
      leg->Draw();

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
      TH1F *hSys = (TH1F*)hJpsiEff[0]->Clone("hSysPlot");
      hSys->Reset();
      hSys->GetYaxis()->SetRangeUser(0.96,1.04);
      hSys->SetTitle("Systematic uncertainty for smearing;p_{T} (GeV/c)");
      c = draw1D(hSys);
      gEffSys->SetFillStyle(1001);
      gEffSys->SetFillColor(kGray);
      gEffSys->Draw("sames e2");
      TLine *line = GetLine(0,1,10,1,1);
      line->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/SmearSys_cent%s.pdf",run_type,cent_Title[icent]));
      fscan->cd();
      gEffSys->Write("",TObject::kOverwrite);
    }

  
}

//-------------------------------------------------------
void smear(const double mass, const int icent, const int nExpr, const double shift, const double sigma, 
	   TF1 *hInputPt, TH2F *hRcJpsiMassVsPt, const bool isWeight, const bool debug)
{
  const int nHisto = hTrkResVsPt[icent]->GetNbinsX();
  hRcJpsiMassVsPt->Sumw2();
  TF1 *hTrkResTmp = (TF1*)funcTrkRes[icent]->Clone(Form("hTrkRes"));
  hTrkResTmp->SetParameter(0, sigma);
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  = myRandom->Uniform(0,20);
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

      if(pt1<0.5 || fabs(eta1) > 0.5) continue;
      if(pt2<0.5 || fabs(eta2) > 0.5) continue;

      double rc_pt1 = pt1 - (hDeltaPt->GetRandom() * hTrkResTmp->Eval(pt1)/0.01) * pt1;
      double rc_pt2 = pt2 - (hDeltaPt->GetRandom() * hTrkResTmp->Eval(pt2)/0.01) * pt2;

      /*
      //if(fabs(eta1) > 0.5 || fabs(eta2) > 0.5) continue;

      // momentum resolution & shift
      int mom_index1 = hTrkResVsPt[icent]->GetXaxis()->FindFixBin(pt1)-1;
      if(mom_index1>=nHisto) mom_index1=nHisto-1;
      int mom_index2 = hTrkResVsPt[icent]->GetXaxis()->FindFixBin(pt2)-1;
      if(mom_index2>=nHisto) mom_index2=nHisto-1;

      double dpt1 = hTrkResBin[icent][mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[icent][mom_index2]->GetRandom();
      double emb_pt1 = (1-dpt1) * pt1;
      double emb_pt2 = (1-dpt2) * pt2;

      double rc_pt1 = 0, rc_pt2 = 0;
      if(year==2013 || year ==2014 || year==2015)
	{
	  rc_pt1 = emb_pt1 * myRandom->Gaus(1+shift/sqrt(emb_pt1),TMath::Power(emb_pt1, 0.95)*sigma);
	  rc_pt2 = emb_pt2 * myRandom->Gaus(1+shift/sqrt(emb_pt2),TMath::Power(emb_pt2, 0.95)*sigma);
	}
      else if(year==2016)
	{
	  // double boundary = 2;
	  // int sign1 = emb_pt1 > boundary ? 1 : -1;
	  // int sign2 = emb_pt2 > boundary ? 1 : -1;
	  // rc_pt1 = emb_pt1 * myRandom->Gaus(1+sign1 * TMath::Power(fabs(emb_pt1-boundary),1)*shift,sqrt(emb_pt1)*sigma);
	  // rc_pt2 = emb_pt2 * myRandom->Gaus(1+sign2 * TMath::Power(fabs(emb_pt2-boundary),1)*shift,sqrt(emb_pt2)*sigma);
	  rc_pt1 = emb_pt1 * myRandom->Gaus(1+shift,TMath::Power(emb_pt1,1)*sigma);
	  rc_pt2 = emb_pt2 * myRandom->Gaus(1+shift,TMath::Power(emb_pt2,1)*sigma);
	}
      */
      double leadpt = rc_pt1 > rc_pt2 ? rc_pt1 : rc_pt2;
      double subpt  = rc_pt1 < rc_pt2 ? rc_pt1 : rc_pt2;
      if(leadpt<pt1_cut || subpt<pt2_cut) continue;
      TLorentzVector rc_daughter1, rc_daughter2;
      // rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      // rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);
      //cout << pt1 << "  " << emb_pt1 << "  " << rc_pt1 << endl;
      if(debug)
	{
	  printf("\n+++++ Before +++++\n");
	  printf("Jpsi:   px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f, pt = %4.2f, M = %4.2f\n",
		 parent.Px(), parent.Py(), parent.Pz(),parent.Phi(), parent.Eta(), parent.E(),parent.Pt(),parent.M());
	  printf("Muon1: px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f\n",
		 daughter1.Px(), daughter1.Py(), daughter1.Pz(), daughter1.Phi(), daughter1.Eta(), daughter1.E());
	  printf("Muon2: px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f\n",
		 daughter2.Px(), daughter2.Py(), daughter2.Pz(), daughter2.Phi(), daughter2.Eta(), daughter2.E());
	}

      rc_daughter1.SetPtEtaPhiM(rc_pt1,eta1,phi1,muMass);
      rc_daughter2.SetPtEtaPhiM(rc_pt2,eta2,phi2,muMass);

      // rc_daughter1.SetXYZM(rc_pt1*cos(phi1), rc_pt1*sin(phi1), daughter1.Pz(), muMass);  
      // rc_daughter2.SetXYZM(rc_pt2*cos(phi2), rc_pt2*sin(phi2), daughter2.Pz(), muMass);  
      TLorentzVector rc_parent = rc_daughter1 + rc_daughter2;

      if(debug)
	{
	  printf("=== after\n");
	  printf("Jpsi:   px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f, pt = %4.2f, M = %4.2f\n",
		 rc_parent.Px(), rc_parent.Py(), rc_parent.Pz(), rc_parent.Phi(), rc_parent.Eta(), rc_parent.E(),rc_parent.Pt(),rc_parent.M());
	  printf("After1:  px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f\n",
		 rc_daughter1.Px(), rc_daughter1.Py(), rc_daughter1.Pz(), rc_daughter1.Phi(), rc_daughter1.Eta(), rc_daughter1.E());
	  printf("After2:  px = %4.2f, py = %4.2f, pz = %4.2f, phi = %4.2f, eta = %4.2f, e = %4.2f\n",
		 rc_daughter2.Px(), rc_daughter2.Py(), rc_daughter2.Pz(), rc_daughter2.Phi(), rc_daughter2.Eta(), rc_daughter2.E());
	}

      double weight = 1; 
      if(isWeight) weight = hInputPt->Eval(mc_pt);
      hRcJpsiMassVsPt->Fill(rc_parent.Pt(),rc_parent.M(),weight);
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
