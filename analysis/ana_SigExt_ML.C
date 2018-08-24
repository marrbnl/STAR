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

// +++ Run14 analysis +++
const double filter_factor = 0.485; 
const char *run_config = "";
//const char *run_config = "Dca3cm.";
//const char* run_config = "Dca3cm.prod_high.";
char *run_type = "Run14_AuAu200";
const int gApplyWeight = 1;
const char *gWeightName[2] = {"","_weight"};
const int gNZdcRate = 11;
const int gNTrgSetup = 5;
const char *gTrgSetupName[5] = {"","_P0","_P1","_P2","_P3"};
const char *gTrgSetupTitle[5] = {"","_prod","_prod_low","_prod_mid","_prod_high"};
const int start_run = 15074000;
const int end_run   = 15168000;
const double jpsi_rapidity = 0.5;
const double pt1_cut = 1.5, pt2_cut = 1.3;
const double low_mass = 2.92, high_mass = 3.32;
const double rej_low_mass = 2.9, rej_high_mass = 3.3;
// analysis of pt dependence
const int nPtBins_pt = 10;
const double ptBins_low_pt[nPtBins_pt]  = {0,0,1,2,3,4,5,6,8,10};
const double ptBins_high_pt[nPtBins_pt] = {15,1,2,3,4,5,6,8,10,15};
const char *pt_Name_pt[nPtBins_pt] = {"0-15","0-1","1-2","2-3","3-4","4-5","5-6","6-8","8-10","10-15"};
const char *pt_Title_pt[nPtBins_pt] = {"0.15 < p_{T} < 15 GeV/c","0.15 < p_{T} < 1 GeV/c","1 < p_{T} < 2 GeV/c","2 < p_{T} < 3 GeV/c","3 < p_{T} < 4 GeV/c","4 < p_{T} < 5 GeV/c","5 < p_{T} < 6 GeV/c","6 < p_{T} < 8 GeV/c","8 < p_{T} < 10 GeV/c","10 < p_{T} < 15 GeV/c"};
const int nCentBins_pt = 5;
const int centBins_low_pt[nCentBins_pt]  = {1,13,9,5,1};
const int centBins_high_pt[nCentBins_pt] = {16,16,12,8,4};
const char *cent_Name_pt[nCentBins_pt] = {"0-80","0-20","20-40","40-60","60-80"};
const char *cent_Title_pt[nCentBins_pt] = {"0080","0020","2040","4060","6080"};
/* const int nCentBins_pt = 1; */
/* const int centBins_low_pt[nCentBins_pt]  = {7}; */
/* const int centBins_high_pt[nCentBins_pt] = {8}; */
/* const char *cent_Name_pt[nCentBins_pt] = {"40-50"}; */
/* const char *cent_Title_pt[nCentBins_pt] = {"4050"}; */
// analysis of npart dependence
const int nPtBins_npart = 2;
const double ptBins_low_npart[nPtBins_npart]  = {0,5};
const double ptBins_high_npart[nPtBins_npart] = {15,15};
const char *pt_Name_npart[nPtBins_npart] = {"0-15","5-15"};
const char *pt_Title_npart[nPtBins_npart] = {"p_{T} > 0.15 GeV/c","p_{T} > 5 GeV/c"};
const int nCentBins_npart[nPtBins_npart] = {8,7};
const int centBins_low_npart[16]  = { 15,13,11,9,7,5,3,1, 15,13,11,9,7,5,1,3 };
const int centBins_high_npart[16] = { 16,14,12,10,8,6,4,2, 16,14,12,10,8,6,4,4 };
const char *cent_Name_npart[16] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80", "0-10","10-20","20-30","30-40","40-50","50-60","60-80","60-70"};
const char *cent_Title_npart[16] = {"0010","1020","2030","3040","4050","5060","6070","7080", "0010","1020","2030","3040","4050","5060","6080","6070"};

double funcTotalPol1(double *x, double *par);
double funcTotalPol2(double *x, double *par);
double funcTotalPol3(double *x, double *par);
double funcTotalPol4(double *x, double *par);

void ana_SigExt();
void prod_pt();


TFile *f = 0x0;
const int year = YEAR;
TString run_cfg_name;
const int useEmbWidth = 1;
TH1F *hSignalTemplate;
TF1  *funcBkgTemplate;
TH1F *hMixBkgTemplate;
//================================================
void ana_SigExt()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  f = TFile::Open(Form("./output/%s.jpsi.%sroot",run_type,run_config),"read");
  run_cfg_name = Form("%s",run_config);

  if(f)
    {
      TH1F *hStat = (TH1F*)f->Get("hEventStat");
      printf("%s\n",f->GetName());
      printf("all         events: %4.4e\n",hStat->GetBinContent(1));
      printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
      printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
    }

  
  //prod_pt();
  MLfit();
}


//================================================
void prod_pt()
{
  for(int t=0; t<5; t++)
    {
      for(int i=0; i<nCentBins_pt; i++)
	{
	  //for(int j=0; j<13; j++)
	  for(int j=0; j<1; j++)
	    {
	      if(t>0 && j>0) continue;
	      //Run14_signal(t,i,j,1,1,0);
	      MLfit(t,i,1,1);
	    }
	}
    }
}

//================================================
double funcBkgTmpPol3(double *x, double *par)
{
  double xx = x[0];
  int bin = hMixBkgTemplate->GetXaxis()->FindFixBin(xx);
  double xmin = hMixBkgTemplate->GetXaxis()->GetBinLowEdge(bin);
  double xmax = hMixBkgTemplate->GetXaxis()->GetBinUpEdge(bin);
  double bkg = hMixBkgTemplate->GetBinContent(bin);
  double sig = 0, res = 0;
  const int nstep = 10;
  const double size = (xmax-xmin)/nstep;
  for(int i=0; i<nstep; i++)
    {
      double xtemp = xmin+(i+0.5)*size;
      sig += par[0]*TMath::Gaus(xtemp, par[1], par[2], true)*size;
      res += (par[3] + par[4]*xtemp + par[5]*xtemp*xtemp + par[6]*xtemp*xtemp*xtemp)*size;
    }
  sig = sig/(xmax-xmin);
  res = res/(xmax-xmin);
  return sig+bkg+res;
}


//===============================================
void MLfit(const int isetup = 1, const int icent = 2, int savePlot = 0, int saveHisto = 0)
{
  gStyle->SetOptFit(0);

  // re-assign global constants
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  xbins[0] = 0.15;

  // prepare name and title
  const TString suffix = Form("cent%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]);
  const TString suf_title = Form("cent%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]);

  // global setup
  double  g_mix_scale_low = 2.7;
  double  g_mix_scale_high = 3.8;
  TString g_mix_func = "pol1";
  double  g_bin_width_1 = 0.02;
  double  g_bin_width_2 = 0.02;
  double  g_bin_width_3 = 0.02;
  TString g_func1 = "pol3";
  int     g_func1_npar = 5;
  TString g_func2 = "pol1";
  int     g_func2_npar = 2;
  double  g_sig_fit_min = 2.5;
  double  g_sig_fit_max = 4.0;

  // get the Jpsi width from smeared simulation
  TFile *finput = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidth = (TH1F*)finput->Get(Form("SmearEmb_JpsiWidth_cent%s_def",cent_Title[icent]));

  // get histograms
  TFile* fin = TFile::Open(Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config),"read");
  printf("\n===== Running configuration =====\n");
  printf("Input filename: %s\n",fin->GetName());
  printf("Histogram name: %s\n",suffix.Data());
  printf("Centrality: %s\n",cent_Name[icent]);
  printf("==================================\n\n");

  TH1F *hSeUL[nPtBins];
  TH1F *hSeLS[nPtBins];
  TH1F *hMixUL[nPtBins];
  TH1F *hMixLS[nPtBins];
  TH1F *hMixBkg[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      hSeUL[i] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s%s%s",pt_Name[i],cent_Name[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]));
      hSeLS[i] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s%s%s",pt_Name[i],cent_Name[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]));
      hMixUL[i] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixLS[i] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
    }
  const int binWidthScale = int(hSeUL[0]->GetBinWidth(1)/hMixUL[0]->GetBinWidth(1));

  // mixed event scaling
  TH1F *hMixScale[nPtBins];
  TF1 *funcScale[nPtBins];
  TCanvas *cScaling = new TCanvas(Form("mix_scale_%s",cent_Name[icent]),Form("mix_scale_%s",cent_Name[icent]),1100,650);
  int nxpad, nypad;
  if(nPtBins-1<=4) { nxpad = 2; nypad = 2; }
  else             { nxpad = (nPtBins-1)/2 + (nPtBins-1)%2; nypad = 2; }
  cScaling->Divide(nxpad,nypad);
  for(int i=0; i<nPtBins; i++)
    {
      hMixScale[i] = (TH1F*)hSeLS[i]->Clone(Form("LS_ratio_SE_to_ME_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixScale[i]->Rebin(5);
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone",hMixLS[i]->GetName()));
      htmp->Rebin(int(hMixScale[i]->GetBinWidth(1)/hMixLS[i]->GetBinWidth(1)));
      hMixScale[i]->Divide(htmp);

      double g_mix_scale_low_tmp = g_mix_scale_low;
      double g_mix_scale_high_tmp = g_mix_scale_high;

      // fitting method
      funcScale[i] = new TF1(Form("Fit_%s",hMixScale[i]->GetName()),g_mix_func,g_mix_scale_low_tmp,g_mix_scale_high_tmp);
      hMixScale[i]->Fit(funcScale[i],"IR0Q");

      if(i>0)
	{
	  // plotting
	  cScaling->cd(i);
	  hMixScale[i]->SetTitle("");
	  hMixScale[i]->SetMarkerStyle(21);
	  if(i==1) hMixScale[i]->GetXaxis()->SetRangeUser(2.7,4);
	  else hMixScale[i]->GetXaxis()->SetRangeUser(2.5,4);
	  hMixScale[i]->SetMaximum(1.1*hMixScale[i]->GetMaximum());
	  hMixScale[i]->Draw();
	  funcScale[i]->SetLineColor(2);
	  funcScale[i]->Draw("sames");
	  TPaveText *t = GetTitleText(Form("%s (%s%%)",pt_Title_pt[i],cent_Name[icent]),0.06);
	  t->Draw();
	}
    }
  t = GetPaveText(0.2,0.35,0.5,0.6,0.06);
  t->SetTextFont(62);
  t->AddText("LS: SE/ME");
  t->SetTextColor(4);
  cScaling->cd(1);
  t->Draw();

  // Acceptance
  TCanvas *cAcc = new TCanvas(Form("Acceptance_%s",cent_Name[icent]),Form("Acceptance_%s",cent_Name[icent]),1100,650);
  cAcc->Divide(nxpad,nypad);
  TH1F *hAcc[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      double g_bin_width_tmp = g_bin_width_1;
      if(i>=2) g_bin_width_tmp = g_bin_width_2;
      if(icent==4 && i==6) g_bin_width_tmp = g_bin_width_3;
      if((icent==0 || icent==1) && i==nPtBins-1) g_bin_width_tmp = g_bin_width_3;
      if((icent==2 || icent==3) && i==nPtBins-2) g_bin_width_tmp = g_bin_width_3;
      hAcc[i] = (TH1F*)hMixUL[i]->Clone(Form("Mix_Acceptance_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone2",hMixLS[i]->GetName()));
      hAcc[i]->Rebin(g_bin_width_tmp/hAcc[i]->GetBinWidth(1));
      htmp->Rebin(g_bin_width_tmp/htmp->GetBinWidth(1));
      hAcc[i]->Divide(htmp);
      if(i>0)
	{
	  cAcc->cd(i);
	  hAcc[i]->SetTitle(";M [GeV/c^{2}]");
	  hAcc[i]->SetMarkerStyle(21);
	  hAcc[i]->GetXaxis()->SetRangeUser(2,4);
	  hAcc[i]->GetYaxis()->SetRangeUser(0.9,1.1);
	  hAcc[i]->Draw("PE");
	  TPaveText *t = GetTitleText(Form("%s (%s%%)",pt_Title_pt[i],cent_Name[icent]),0.06);
	  t->Draw();
	  TLine *line = GetLine(2,1,4,1,1);
	  line->Draw();
	}
    }
  t = GetPaveText(0.3,0.4,0.2,0.3,0.07);
  t->SetTextFont(62);
  t->AddText("Mix: UL/LS");
  t->SetTextColor(4);
  cAcc->cd(1);
  t->Draw();

  // jpsi signal
  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",cent_Name[icent]),Form("InvMass_%s",cent_Name[icent]),1100,650);
  cSignal->Divide(nxpad,nypad);
  for(int i=0; i<nPtBins; i++)
    {
      hMixBkg[i] = (TH1F*)hMixUL[i]->Clone(Form("mix_bkg_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      
      // scale mix event background
      int low_bin = hMixBkg[i]->FindFixBin(2.3);
      int up_bin  = hMixBkg[i]->FindFixBin(4.2);
      for(int ibin=low_bin; ibin<=up_bin; ibin++)
	{
	  double mass = hMixBkg[i]->GetBinCenter(ibin);
	  double scale = funcScale[i]->Eval(mass);
	  hMixBkg[i]->SetBinContent(ibin, scale * hMixBkg[i]->GetBinContent(ibin));
	  hMixBkg[i]->SetBinError(ibin, scale * hMixBkg[i]->GetBinError(ibin));
	}
      
      double g_bin_width_tmp = g_bin_width_1;
      if(i>=2) g_bin_width_tmp = g_bin_width_2;
      if(icent==4 && i==6) g_bin_width_tmp = g_bin_width_3;
      if((icent==0 || icent==1) && i==nPtBins-1) g_bin_width_tmp = g_bin_width_3;
      if((icent==2 || icent==3) && i==nPtBins-2) g_bin_width_tmp = g_bin_width_3;
      hSeUL[i]->Rebin(g_bin_width_tmp/hSeUL[i]->GetBinWidth(1));
      hSeUL[i]->SetTitle("");
      hSeUL[i]->SetMarkerStyle(21);
      hSeUL[i]->SetMarkerColor(2);
      hSeUL[i]->SetLineColor(2);
      hSeUL[i]->GetXaxis()->SetRangeUser(2.5,4);
      hSeLS[i]->Rebin(hSeUL[i]->GetBinWidth(1)/hSeLS[i]->GetBinWidth(1));
      hMixBkg[i]->SetLineColor(4);
      hMixBkg[i]->Rebin(hSeUL[i]->GetBinWidth(1)/hMixBkg[i]->GetBinWidth(1));
      if(i>0)
	{
	  cSignal->cd(i);
	  hSeUL[i]->Draw();
	  hSeLS[i]->Draw("sames HIST");
	  hMixBkg[i]->Draw("sames HIST");
	  TPaveText *t = GetTitleText(Form("%s (%s%%)",pt_Title_pt[i],cent_Name[icent]),0.06);
	  t->Draw();
	}
    }
  cSignal->cd(10);
  leg = new TLegend(0.2,0.6,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.07);
  leg->AddEntry(hSeUL[1],"Unlike sign","P");
  leg->AddEntry(hSeLS[1],"Like sign (++)+(--)","L");
  leg->AddEntry(hMixBkg[1],"Mix UL","L");
  leg->Draw();
  

  // try three different methods
  // 1 - fit UL with likelihood method
  // 2 - fit UL with chi2 method
  // 3 - fit UL-ME with chi2 method
  const int nMethods = 3;
  TH1F *hSignal[nPtBins][nMethods];
  TF1 *funcSignal[nPtBins][nMethods];
  TF1 *funcBkg[nPtBins][nMethods];
  TH1F *hChi2[nMethods], *hMean[nMethods], *hSigma[nMethods], *hFitYield[nMethods], *hBinCountYield[nMethods], *hSigToBkg[nMethods];
  double bkg_params[5];
  double bkg_matrix[5*5];
  TCanvas *cFit1[nMethods], *cFit2[nMethods];
  TCanvas *c = 0x0;
  for(int m=0; m<3; m++)
    {
      if(m>0) continue;
      hChi2[m] = new TH1F(Form("Jpsi_Chi2_%s_M%d",suffix.Data(),m),"Chi2/NDF of Jpsi fit",nbins,xbins);
      hMean[m] = new TH1F(Form("Jpsi_FitMean_%s_M%d",suffix.Data(),m),"Mean of Jpsi peak",nbins,xbins);
      hSigma[m]= new TH1F(Form("Jpsi_FitSigma_%s_M%d",suffix.Data(),m),"Sigma of Jpsi peak",nbins,xbins);
      hFitYield[m] = new TH1F(Form("Jpsi_FitYield_%s_M%d",suffix.Data(),m),"Jpsi yield from fitting",nbins,xbins);
      hBinCountYield[m] = new TH1F(Form("Jpsi_BinCountYield_%s_M%d",suffix.Data(),m),"Jpsi yield from bin counting",nbins,xbins);
      hSigToBkg[m] = new TH1F(Form("Jpsi_SigToBkg_%s_M%d",suffix.Data(),m),"Signal-to-background ratio for Jpsi",nbins,xbins);

      cFit2[m] = new TCanvas(Form("Fit2_Jpsi_%s_M%d",cent_Name[icent],m),Form("Fit2_Jpsi_%s_M%d",cent_Name[icent],m),1400,650);
      cFit2[m]->Divide(5,2);

      cFit1[m] = new TCanvas(Form("Fit1_Jpsi_%s_M%d",cent_Name[icent],m),Form("Fit1_Jpsi_%s_M%d",cent_Name[icent],m),1400,650);
      cFit1[m]->Divide(5,2);

      TCanvas *ctmp = new TCanvas("ctmp","ctmp",800,600);

      TString funcForm;
      int nPar;
      for(int i=0; i<nPtBins; i++)
	{
	  if(icent==3 && i==nPtBins-1) continue;
	  if(icent==4 && i>=nPtBins-3) continue;

	  printf("+++++ %1.0f < pT < %1.0f +++++\n",ptBins_low[i],ptBins_high[i]);
	  hSignal[i][m] = (TH1F*)hSeUL[i]->Clone(Form("Jpsi_Signal_cent%s_pt%s_M%d",cent_Title[icent],pt_Name[i],m));
	  hSignal[i][m]->SetMarkerStyle(24);
	  hSignal[i][m]->SetMarkerSize(0.8);
	  hSignal[i][m]->SetLineColor(1);
	  hSignal[i][m]->SetMarkerColor(1);
	  TFitResultPtr ptr;
	  if(m==0)
	    {
	      nPar = 4;
	      funcForm = "pol3";
	      hMixBkgTemplate = hMixBkg[i];
	      funcSignal[i][m] = new TF1(Form("Jpsi_FitSig_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),funcBkgTmpPol3,g_sig_fit_min,g_sig_fit_max,7);

	      //funcSignal[i][m] = new TF1(Form("Jpsi_FitSig_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),Form("gausn(0)+pol3(3)"),g_sig_fit_min,g_sig_fit_max);
	      funcSignal[i][m]->SetParameters(hSignal[i][m]->GetMaximum(), 3.09, 0.06, 100, -10, 1, 1);
	      if(i>0) funcSignal[i][m]->FixParameter(2,hEmbJpsiWidth->GetBinContent(i));
	      if(i==1) funcSignal[i][m]->SetRange(2.75,g_sig_fit_max);
	      ptr = hSignal[i][m]->Fit(funcSignal[i][m],"IRS0L");
	    }
	  else if(m==1)
	    {
	      nPar = 4;
	      funcForm = "pol3";

	      hMixBkgTemplate = hMixBkg[i];
	      funcSignal[i][m] = new TF1(Form("Jpsi_FitSig_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),funcBkgTmpPol3,g_sig_fit_min,g_sig_fit_max,7);
	      //funcSignal[i][m] = new TF1(Form("Jpsi_FitSig_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),Form("gausn(0)+pol3(3)"),g_sig_fit_min,g_sig_fit_max);
	      funcSignal[i][m]->SetParameters(hSignal[i][m]->GetMaximum(), 3.09, 0.06, 100, -10, 1, 1);
	      if(i>0) funcSignal[i][m]->FixParameter(2,hEmbJpsiWidth->GetBinContent(i));
	      if(i==1) funcSignal[i][m]->SetRange(2.75,g_sig_fit_max);
	      ptr = hSignal[i][m]->Fit(funcSignal[i][m],"IRS0Q");
	    }
	  else if(m==2)
	    {
	      for(int bin=1; bin<=hSignal[i][m]->GetNbinsX(); bin++)
		{
		  if(hSignal[i][m]->GetBinContent(bin)==0)
		    {
		      hSignal[i][m]->SetBinContent(bin,0);
		      hSignal[i][m]->SetBinError(bin,1.4);
		    }
		}
	      hSignal[i][m]->Add(hMixBkg[i],-1);
	      hSignal[i][m]->SetLineColor(1);
	      hSignal[i][m]->SetMarkerColor(1);

	      // Fit signal
	      if(i<1)
		{
		  funcForm = g_func1; 
		  nPar = g_func1_npar;
		}
	      else
		{
		  funcForm = g_func2; 
		  nPar = g_func2_npar;
		}

	      funcSignal[i][m] = new TF1(Form("Jpsi_FitSig_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	      funcSignal[i][m]->SetParameters(hSignal[i][m]->GetMaximum(), 3.09, 0.06);
	      if(i>0) funcSignal[i][m]->FixParameter(2,hEmbJpsiWidth->GetBinContent(i));
	      if(i==1) funcSignal[i][m]->SetRange(2.75,g_sig_fit_max);
	      ptr = hSignal[i][m]->Fit(funcSignal[i][m],"IRS0Q");
	    }

	  double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
	  double bin_width = hSignal[i][m]->GetBinWidth(1);
	  double fit_yield = funcSignal[i][m]->GetParameter(0)/bin_width;
	  double fit_yield_err = funcSignal[i][m]->GetParError(0)/bin_width;
	  double fit_mean = funcSignal[i][m]->GetParameter(1);
	  double fit_sigma = fabs(funcSignal[i][m]->GetParameter(2));
	  printf("Fitting = %1.1f +/- %1.1f (%1.1fsigma)\n",fit_yield,fit_yield_err,fit_yield/fit_yield_err);
	  hFitYield[m]->SetBinContent(i,fit_yield);
	  hFitYield[m]->SetBinError(i,fit_yield_err);
	  hChi2[m]->SetBinContent(i, funcSignal[i][m]->GetChisquare()/funcSignal[i][m]->GetNDF());   
	  hChi2[m]->SetBinError(i, 1e-10);
	  hMean[m]->SetBinContent(i,funcSignal[i][m]->GetParameter(1));
	  hMean[m]->SetBinError(i,funcSignal[i][m]->GetParError(1));
	  hSigma[m]->SetBinContent(i,fit_sigma);
	  hSigma[m]->SetBinError(i,funcSignal[i][m]->GetParError(2));

	  // Extract background matrix
	  funcBkg[i][m] = new TF1(Form("Jpsi_FitBkg_pt%s_%s_M%d",pt_Name[i],suffix.Data(),m),Form("%s",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	  int nSigPar = 3;
	  for(int j=0; j<nPar; j++)
	    {
	      funcBkg[i][m]->SetParameter(j,funcSignal[i][m]->GetParameter(nSigPar+j));
	      funcBkg[i][m]->SetParError(j,funcSignal[i][m]->GetParError(nSigPar+j));
	      bkg_params[j] = funcSignal[i][m]->GetParameter(nSigPar+j);
	    }
 
	  for(int j=nSigPar; j<nSigPar+nPar; j++)
	    {
	      for(int k=nSigPar; k<nSigPar+nPar; k++)
		{
		  bkg_matrix[(j-nSigPar)*nPar+k-nSigPar] = matrix[j*(nSigPar+nPar)+k];
		}
	    }

	  // bin counting
	  double low_mass_tmp = fit_mean - 3.5 * fit_sigma;
	  double high_mass_tmp = fit_mean + 3.5 * fit_sigma;
	  int low_bin = hSignal[i][m]->FindFixBin(low_mass_tmp+1e-4) + 1;
	  int high_bin = hSignal[i][m]->FindFixBin(high_mass_tmp-1e-4) - 1;
	  double low_bin_mass = hSignal[i][m]->GetXaxis()->GetBinLowEdge(low_bin);
	  double high_bin_mass = hSignal[i][m]->GetXaxis()->GetBinUpEdge(high_bin);
	  double yield_all_err;
	  double yield_all = hSignal[i][m]->IntegralAndError(low_bin,high_bin,yield_all_err);
	  double yield_bkg = funcBkg[i][m]->Integral(low_bin_mass,high_bin_mass) * 1./hSignal[i][m]->GetBinWidth(1);
	  double yield_bkg_err = funcBkg[i][m]->IntegralError(low_bin_mass,high_bin_mass,bkg_params,bkg_matrix)* 1./hSignal[i][m]->GetBinWidth(1);
	  if(m<2)
	    {
	      // account for mixed event background
	      double mix_bkg_err;
	      double mix_bkg = hMixBkg[i]->IntegralAndError(low_bin, high_bin, mix_bkg_err);
	      yield_bkg = yield_bkg + mix_bkg;
	      yield_bkg_err = sqrt(yield_bkg_err*yield_bkg_err + mix_bkg_err*mix_bkg_err);
	    }

	  TF1 *funcJpsi_temp = new TF1(Form("funcJpsi_temp_%d",i),"gausn");
	  for(int par=0; par<3; par++)
	    {
	      funcJpsi_temp->SetParameter(par, funcSignal[i][m]->GetParameter(par));
	      funcJpsi_temp->SetParError(par, funcSignal[i][m]->GetParError(par));
	    }
	  double efficiency = funcJpsi_temp->Integral(low_bin_mass,high_bin_mass)/funcJpsi_temp->GetParameter(0);
	  double yield_sig = (yield_all - yield_bkg)/efficiency;
	  double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err)/efficiency;
	  double all = hSeUL[i]->Integral(hSeUL[i]->FindBin(low_mass_tmp+1e-4),hSeUL[i]->FindBin(high_mass_tmp-1e-4));
	  hBinCountYield[m]->SetBinContent(i,yield_sig);
	  hBinCountYield[m]->SetBinError(i,yield_sig_err);
	  hSigToBkg[m]->SetBinContent(i,yield_sig/(all-yield_sig));
	  printf("Count = %1.1f +/- %1.1f (%1.1fsigma)\n",yield_sig,yield_sig_err,yield_sig/yield_sig_err);
	  printf("Background = %1.1f +/- %1.1f\n",yield_bkg,yield_bkg_err);
	  printf("all = %1.2f, signal = %1.2f\n",all,yield_sig);

	  // plotting
	  funcSignal[i][m]->SetLineStyle(2);
	  funcSignal[i][m]->SetLineColor(2);
	  TH1F *funcSignalHist = (TH1F*)funcSignal[i][m]->GetHistogram();
	  TH1F *fitRatio = (TH1F*)hSignal[i][m]->Clone(Form("%s_ratio",hSignal[i][m]->GetName()));
	  int minbin = hSignal[i][m]->FindBin(funcSignal[i][m]->GetXmin());
	  int maxbin = hSignal[i][m]->FindBin(funcSignal[i][m]->GetXmax())-1;
	  for(int bin=minbin; bin<=maxbin; bin++)
	    {
	      if(m<2)
		{
		  double scale = funcSignalHist->GetBinContent(funcSignalHist->FindBin(hSignal[i][m]->GetBinCenter(bin)));
		  if(scale>0)
		    {
		      fitRatio->SetBinContent(bin, fitRatio->GetBinContent(bin)/scale);
		      fitRatio->SetBinError(bin, fitRatio->GetBinError(bin)/scale);
		    }
		  else
		    {
		      fitRatio->SetBinContent(bin, -1);
		      fitRatio->SetBinError(bin, 1e-10);
		    }
		}
	      else
		{
		  double fityield = funcSignal[i][m]->Integral(fitRatio->GetXaxis()->GetBinLowEdge(bin), fitRatio->GetXaxis()->GetBinUpEdge(bin))/fitRatio->GetBinWidth(bin);
		  fitRatio->SetBinContent(bin, fitRatio->GetBinContent(bin)-fityield);
		}
	    }

	  if(m<2)
	    {
	      fitRatio->SetYTitle("Data/Fit");
	      if(i<nPtBins/2) fitRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
	      else if(i<8) fitRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
	      else    fitRatio->GetYaxis()->SetRangeUser(0, 2);
	    }
	  else
	    {
	      fitRatio->SetYTitle("Data-Fit");
	      fitRatio->SetMaximum(1.8*fitRatio->GetMaximum());
	    }
	  int index;
	  if(i<nPtBins/2) 
	    {
	      c = cFit1[m];
	      index = i;
	    }
	  else
	    {
	      c = cFit2[m];
	      index = i-nPtBins/2;
	    }
	  
	  c->cd(index+6);
	  fitRatio->GetYaxis()->SetTitleOffset(1.3);
	  fitRatio->Draw();

	  c->cd(index+1);
	  hSignal[i][m]->SetMaximum(1.5*hSignal[i][m]->GetMaximum());
	  hSignal[i][m]->DrawCopy();
	  if(m<2) hMixBkg[i]->Draw("samesHIST");
	  funcSignalHist->Draw("sames");
	  if(m==2)
	    {
	      funcBkg[i][m]->SetLineColor(4);
	      funcBkg[i][m]->SetLineStyle(2);
	      funcBkg[i][m]->Draw("sames");
	    }

	  TPaveText *t = GetTitleText(Form("%s (%s%%)",pt_Title_pt[i],cent_Name[icent]),0.05);
	  t->Draw();
	  TLine *line = GetLine(low_bin_mass,hSignal[i][m]->GetMinimum(),low_bin_mass,hSignal[i][m]->GetMaximum()*0.6,1);
	  line->Draw();
	  TLine *line = GetLine(high_bin_mass,hSignal[i][m]->GetMinimum(),high_bin_mass,hSignal[i][m]->GetMaximum()*0.6,1);
	  line->Draw();

	  t = GetPaveText(0.16,0.3,0.75,0.88,0.045);
	  t->SetTextFont(62);
	  t->AddText(Form("Fitting: %1.0f #pm %1.0f (%1.1f#sigma)",fit_yield,fit_yield_err,fit_yield/fit_yield_err));
	  t->AddText(Form("Counting: %1.0f #pm %1.0f (%1.1f#sigma)",yield_sig,yield_sig_err,yield_sig/yield_sig_err));
	  t->SetTextAlign(11);
	  t->Draw();
	}
      t = GetPaveText(0.7,0.8,0.3,0.35,0.06);
      t->SetTextFont(62);
      if(m<2) t->AddText("UL");
      else t->AddText("UL-MIX(UL)");
      t->SetTextColor(4);
      cFit1[m]->cd(1);
      t->Draw();
      cFit2[m]->cd(1);
      t->Draw();

      if(savePlot) 
	{
	  cFit1[m]->SaveAs(Form("./Plots/%s/ana_SigExt/ML_M%d_SigFit1_%s%s.pdf",run_type,m,cent_Title[icent],gTrgSetupTitle[isetup]));
	  cFit2[m]->SaveAs(Form("./Plots/%s/ana_SigExt/ML_M%d_SigFit2_%s%s.pdf",run_type,m,cent_Title[icent],gTrgSetupTitle[isetup]));
	}
    }
  return;

  TList *list = new TList();
  TString legName[nMethods] = {"UL: ML fit", "UL: Chi2 fit", "UL-ME: Chi2 fit"};
  double pt_high = 15;
  if(icent==2 || icent==3) pt_high = 10;
  if(icent==4) pt_high = 6;

  // compare signal mean
  for(int m=0; m<nMethods; m++)
    {
      list->Add(hMean[m]);
    }
  c = drawHistos(list,"FitMean",Form("%s: mean of J/psi peak;p_{T} (GeV/c);Mean",run_type),true,0.15,pt_high,kTRUE,3.04,3.16,false,true,legName,true,Form("%s%%",cent_Name[icent]),0.15,0.35,0.7,0.88,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompMean_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));

  // chi2/ndf
  for(int m=0; m<nMethods; m++)
    {
      list->Add(hChi2[m]);
    }
  c = drawHistos(list,"FitChi2",Form("%s: #chi^{2}/NDF of fits;p_{T} (GeV/c);#chi^{2}/NDF",run_type),true,0.15,pt_high,kTRUE,0,2.5,false,true,legName,true,Form("%s%%",cent_Name[icent]),0.15,0.35,0.7,0.88,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompChi2_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));

  // fit yield
  for(int m=0; m<nMethods; m++)
    {
      list->Add(hFitYield[m]);
    }
  c = drawHistos(list,"FitYield",Form("%s: J/psi yield from fitting;p_{T} (GeV/c);counts",run_type),true,0.15,pt_high,true,1,hFitYield[2]->GetBinContent(1)*1.5,false,true,legName,true,Form("%s%%",cent_Name[icent]),0.55,0.75,0.65,0.85,kTRUE);
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompFitYield_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));

  c = drawHistos(list,"FitYieldLog",Form("%s: J/psi yield from fitting;p_{T} (GeV/c);counts",run_type),true,0.15,pt_high,true,1,hFitYield[2]->GetBinContent(1)*10,true,true,legName,true,Form("%s%%",cent_Name[icent]),0.55,0.75,0.65,0.85,kTRUE);
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompFitYieldLog_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));
  list->Clear();
  
  // signal significance & difference with bin counting
  TH1F *hYieldDiff[nMethods];
  TH1F *hYieldSig[nMethods];
  for(int m=0; m<nMethods; m++)
    {
      hYieldDiff[m] = (TH1F*)hFitYield[m]->Clone(Form("Jpsi_YieldDiff_%s_M%d",suffix.Data(),m));
      hYieldSig[m] = (TH1F*)hFitYield[m]->Clone(Form("Jpsi_YieldSig_%s_M%d",suffix.Data(),m));
      for(int bin=1; bin<=nbins; bin++)
	{
	  double fityield = hFitYield[m]->GetBinContent(bin);
	  double fiterr = hFitYield[m]->GetBinError(bin);
	  double countyield = hBinCountYield[m]->GetBinContent(bin);

	  if(fityield<=0) continue;

	  hYieldSig[m]->SetBinContent(bin, fityield/fiterr);
	  hYieldSig[m]->SetBinError(bin, 1e-10);
	  hYieldDiff[m]->SetBinContent(bin, countyield/fityield - m*0.5);
	  hYieldDiff[m]->SetBinError(bin, fiterr/fityield);
	  cout << m << "  " << bin << "  " << countyield/fityield - m*0.5 << endl;
	}
    }

  for(int m=0; m<nMethods; m++)
    {
      list->Add(hYieldSig[m]);
    }
  c = drawHistos(list,"FitSig",Form("%s: signal significance from fits;p_{T} (GeV/c);sig.",run_type),true,0.15,pt_high,true,0,30,false,true,legName,true,Form("%s%%",cent_Name[icent]),0.55,0.75,0.65,0.85,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompSignalSig_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));


  for(int m=0; m<nMethods; m++)
    {
      list->Add(hYieldDiff[m]);
    }
  c = drawHistos(list,"YieldDiff",Form("%s: J/psi yield difference between fitting and bin-counting;p_{T} (GeV/c);Counting/Fitting",run_type),true,0.15,pt_high,true,-0.5,2,false,true,legName,true,Form("%s%%",cent_Name[icent]),0.55,0.75,0.65,0.85,kTRUE);
  for(int m=0; m<nMethods; m++)
    {
      TLine *line = GetLine(0.15, 1-m*0.5, pt_high, 1-m*0.5, hYieldDiff[m]->GetLineColor());
      line->Draw();
    }
  list->Clear();
  if(savePlot) c->SaveAs(Form("./Plots/%s/ana_SigExt/ML_CompYieldDiff_%s%s.pdf",run_type,cent_Title[icent],gTrgSetupTitle[isetup]));

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config),"update");
      for(int m=0; m<3; m++)
	{
	  hFitYield[m]->Write("",TObject::kOverwrite);
	}
      fout->Close();
    }

}
