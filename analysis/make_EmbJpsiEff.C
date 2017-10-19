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

void makeJpsiTpcEffVsZdc(const int savePlot, const int saveHisto);
void makeDataJpsiWeight(const int saveHisto);
void makeJpsi(const bool saveHisto);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
void getJpsiEff(const double mass, const int nExpr, TH1F *hInputPt, TF1 *funcTpcTrkEff, TH1F *hMcJpsiPt, TH1F *hRcJpsiPt, int debug = 0);
TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);

const int year = 2014;
const char *run_type = "Run14_AuAu200";
const int gNZdcRate = 11;
const int gNTrgSetup = 5;
const char *gTrgSetupName[5] = {"","_P0","_P1","_P2","_P3"};
const char *gTrgSetupTitle[5] = {"","_prod","_prod_low","_prod_mid","_prod_high"};
const double jpsi_rapidity = 0.5;
const double pt1_cut = 1.5, pt2_cut = 1.3;
const int nPtBins_pt = 10;
const double ptBins_low_pt[nPtBins_pt]  = {0,0,1,2,3,4,5,6,8,10};
const double ptBins_high_pt[nPtBins_pt] = {15,1,2,3,4,5,6,8,10,15};
const char *pt_Name_pt[nPtBins_pt] = {"0-15","0-1","1-2","2-3","3-4","4-5","5-6","6-8","8-10","10-15"};
const int nCentBins_pt = 5;
const int centBins_low_pt[nCentBins_pt]  = {1,13,9,5,1};
const int centBins_high_pt[nCentBins_pt] = {16,16,12,8,4};
const char *cent_Name_pt[nCentBins_pt] = {"0-80","0-20","20-40","40-60","60-80"};
const char *cent_Title_pt[nCentBins_pt] = {"0080","0020","2040","4060","6080"};
const int nPtBins_npart = 2;
const double ptBins_low_npart[nPtBins_npart]  = {0,5};
const double ptBins_high_npart[nPtBins_npart] = {15,15};
const char *pt_Name_npart[nPtBins_npart] = {"0-15","5-15"};
const int nCentBins_npart[nPtBins_npart] = {8,7};
const int centBins_low_npart[16]  = { 15,13,11,9,7,5,3,1, 15,13,11,9,7,5,1,3 };
const int centBins_high_npart[16] = { 16,14,12,10,8,6,4,2, 16,14,12,10,8,6,4,4 };
const char *cent_Name_npart[16] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80", "0-10","10-20","20-30","30-40","40-50","50-60","60-80","60-70"};
const char *cent_Title_npart[16] = {"0010","1020","2030","3040","4050","5060","6070","7080", "0010","1020","2030","3040","4050","5060","6080","6070"};

const int nEffType = 8;
const char *trkEffType[nEffType] = {"MC","Tpc","MtdMth","Fake","MuonPid","MtdTrig","TrigUnit","Embed"};
const char *weight_name[2] = {"","_w"};
const double muMass = 0.1057;
const double pi = 3.1415926;

TFile *f;
TCanvas *c;
TRandom3 *myRandom;

//================================================
void make_EmbJpsiEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.Embed.root"),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/%s.Embed.Jpsi.root",run_type),"read");
    }
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("[i] # of events: %4.4e\n",hStat->GetBinContent(3));

  makeJpsi(1);
  //makeDataJpsiWeight(1);
  // makeJpsiTpcEffVsZdc(1, 1);
}

//================================================
void makeJpsiTpcEffVsZdc(const int savePlot, const int saveHisto)
{
  const int* nCentBins      = nCentBins_npart; 
  const int* centBins_low   = centBins_low_npart;
  const int* centBins_high  = centBins_high_npart;
  const char** cent_Name    = cent_Name_npart;
  const char** cent_Title   = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  const int nbins = nPtBins_pt -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low_pt[i+1];
  xbins[nbins] = ptBins_high_pt[nbins];

  myRandom = new TRandom3();
  myRandom->SetSeed(0);
  const int nExpr = 1e6;
  const double mass = 3.096;
  TFile *ftrk = TFile::Open(Form("Rootfiles/%s.EmbTrkEff.root",run_type),"read");
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");

  // calculate Jpsi efficiency due to inclusive TPC tracking
  TF1 *funcTrkEffAll = (TF1*)ftrk->Get("func_TpcTrkEff");
  TH1F *hMcJpsiPtAll = new TH1F(Form("JpsiPtMc_All"),";p_{T} (GeV/c)",15, 0, 15);
  TH1F *hRcJpsiPtAll = new TH1F(Form("JpsiPtRc_All"),";p_{T} (GeV/c)",15, 0, 15);
  TH1F *hTpcJpsiEffAll;
  TH1F *hMcJpsiPtRebinAll;
  TH1F *hRcJpsiPtRebinAll;
  TH1F *hTpcJpsiEffRebinAll;

  TH1F *hInputJpsiAll = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  getJpsiEff(mass, nExpr, hInputJpsiAll, funcTrkEffAll, hMcJpsiPtAll, hRcJpsiPtAll, 0);
  hTpcJpsiEffAll = (TH1F*)hRcJpsiPtAll->Clone(Form("JpsiTpcEff_All"));
  hTpcJpsiEffAll->Divide(hMcJpsiPtAll);

  hMcJpsiPtRebinAll = (TH1F*)hMcJpsiPtAll->Rebin(nbins, Form("%s_rebin",hMcJpsiPtAll->GetName()), xbins);
  hRcJpsiPtRebinAll = (TH1F*)hRcJpsiPtAll->Rebin(nbins, Form("%s_rebin",hRcJpsiPtAll->GetName()), xbins);
  hTpcJpsiEffRebinAll = (TH1F*)hRcJpsiPtRebinAll->Clone(Form("%s_rebin",hTpcJpsiEffAll->GetName()));
  hTpcJpsiEffRebinAll->Divide(hMcJpsiPtRebinAll);

  TH1F *hTpcJpisEffVsCentAll[nPtBins_npart];
  for(int i=0; i<nPtBins_npart; i++)
    {
      hTpcJpisEffVsCentAll[i] = new TH1F(Form("JpsiTpcEffVsCent_Pt%1.0f_All",ptBins_low_npart[i]),"",1,0,1);
      int low_bin = hMcJpsiPtAll->FindFixBin(ptBins_low_npart[i]+1e-4);
      int up_bin = hMcJpsiPtAll->FindFixBin(ptBins_high_npart[i]-1e-4);
      double nMcJpsi = hMcJpsiPtAll->Integral(low_bin, up_bin);
      double nRcJpsi = hRcJpsiPtAll->Integral(low_bin, up_bin);
      hTpcJpisEffVsCentAll[i]->SetBinContent(1, nRcJpsi/nMcJpsi);
      hTpcJpisEffVsCentAll[i]->SetBinError(1, 1e-5);
    }


  // calculate Jpsi efficiency due to TPC tracking vs. centrality vs. ZDCrate
  TH1F *hInputJpsi[kNCent];
  for(int k=0; k<kNCent; k++)
    {
      if(k==0 || k==1) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
      if(k==2 || k==3) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
      if(k>=4)         hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");
      hInputJpsi[k]->SetName(Form("hInputJpsi_cent%s",cent_Name[k]));
    }
  TF1 *funcTrkEff[gNZdcRate][kNCent];
  for(int k=0; k<kNCent; k++)
    {
      for(int j=0; j<gNZdcRate; j++)
	{
	  funcTrkEff[j][k] = (TF1*)ftrk->Get(Form("func_TpcTrkEff_cent%s_Zdc%d-%d",cent_Title[k],j*10,j*10+10));
	}
    }
  TH1F *hJpsiTpcEffVsCent[nPtBins_npart][gNZdcRate][3];
  TH1F *hMcJpsiPt[gNZdcRate][kNCent][3];
  TH1F *hRcJpsiPt[gNZdcRate][kNCent][3];
  TH1F *hTpcJpsiEff[gNZdcRate][kNCent][3];
  TH1F *hMcJpsiPtRebin[gNZdcRate][kNCent][3];
  TH1F *hRcJpsiPtRebin[gNZdcRate][kNCent][3];
  TH1F *hTpcJpsiEffRebin[gNZdcRate][kNCent][3];
  const char *type_name[3] = {"","_up","_down"};
  for(int j=0; j<gNZdcRate; j++)
    {
      for(int i=0; i<nPtBins_npart; i++)
	{
	  for(int t=0; t<3; t++)
	    {
	      hJpsiTpcEffVsCent[i][j][t] = new TH1F(Form("JpsiTpcEffVsCent_Pt%1.0f_Zdc%d-%d%s",ptBins_low_npart[i],j*10,j*10+10,type_name[t]),"",nCentBins[i],0,nCentBins[i]);
	    }
	}

      for(int k=0; k<kNCent; k++)
	{
	  for(int t=0; t<3; t++)
	    {
	      hMcJpsiPt[j][k][t] = new TH1F(Form("JpsiPtMc_cent%s_Zdc%d-%d%s",cent_Title[k],j*10,j*10+10,type_name[t]),";p_{T} (GeV/c)",15, 0, 15);
	      hMcJpsiPt[j][k][t]->Sumw2();
	      hRcJpsiPt[j][k][t] = new TH1F(Form("JpsiPtRc_cent%s_Zdc%d-%d%s",cent_Title[k],j*10,j*10+10,type_name[t]),";p_{T} (GeV/c)",15, 0, 15);
	      hRcJpsiPt[j][k][t]->Sumw2();
	    }
	}
    }

  for(int j=0; j<gNZdcRate; j++)
    {
      for(int k=0; k<kNCent; k++)
	{
	  for(int t=0; t<3; t++)
	    {
	      printf("[i] Centrality %s%%%s, %d < ZDC < %d\n",cent_Name[k],type_name[t],j*10,j*10+10);
	      TF1 *trkEff = (TF1*)funcTrkEff[j][k]->Clone(Form("%s%s",funcTrkEff[j][k]->GetName(),type_name[t]));
	      if(t==1) trkEff->SetParameter(0, funcTrkEff[j][k]->GetParameter(0)+funcTrkEff[j][k]->GetParError(0));
	      if(t==2) trkEff->SetParameter(0, funcTrkEff[j][k]->GetParameter(0)-funcTrkEff[j][k]->GetParError(0));

	      getJpsiEff(mass, nExpr, hInputJpsi[k], trkEff, hMcJpsiPt[j][k][t], hRcJpsiPt[j][k][t], 0);
	      hTpcJpsiEff[j][k][t] = (TH1F*)hRcJpsiPt[j][k][t]->Clone(Form("JpsiTpcEff_cent%s_Zdc%d-%d%s",cent_Title[k],j*10,j*10+10,type_name[t]));
	      hTpcJpsiEff[j][k][t]->Divide(hMcJpsiPt[j][k][t]);

	      hMcJpsiPtRebin[j][k][t] = (TH1F*)hMcJpsiPt[j][k][t]->Rebin(nbins, Form("%s_rebin",hMcJpsiPt[j][k][t]->GetName()), xbins);
	      hRcJpsiPtRebin[j][k][t] = (TH1F*)hRcJpsiPt[j][k][t]->Rebin(nbins, Form("%s_rebin",hRcJpsiPt[j][k][t]->GetName()), xbins);
	      hTpcJpsiEffRebin[j][k][t] = (TH1F*)hRcJpsiPtRebin[j][k][t]->Clone(Form("%s_rebin",hTpcJpsiEff[j][k][t]->GetName()));
	      hTpcJpsiEffRebin[j][k][t]->Divide(hMcJpsiPtRebin[j][k][t]);
	    }
	}
    }

  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int j=0; j<gNZdcRate; j++)
	{
	  for(int t=0; t<3; t++)
	    {
	      for(int k=0; k<nCentBins[i]; k++)
		{
		  int low_bin = hMcJpsiPt[j][k][t]->FindFixBin(ptBins_low_npart[i]+1e-4);
		  int up_bin = hMcJpsiPt[j][k][t]->FindFixBin(ptBins_high_npart[i]-1e-4);
		  double nMcJpsi = hMcJpsiPt[j][k][t]->Integral(low_bin, up_bin);
		  double nRcJpsi = hRcJpsiPt[j][k][t]->Integral(low_bin, up_bin);
		  if(i==1 && k==nCentBins[i]-1)
		    {
		      nMcJpsi += hMcJpsiPt[j][k+1][t]->Integral(low_bin, up_bin);
		      nRcJpsi += hRcJpsiPt[j][k+1][t]->Integral(low_bin, up_bin);
		    }
		  double eff = nRcJpsi/nMcJpsi;
		  hJpsiTpcEffVsCent[i][j][t]->SetBinContent(k+1, eff);
		  hJpsiTpcEffVsCent[i][j][t]->SetBinError(k+1, 1e-5);
		}
	    }
	}
    }

  // plot jpsi efficiency vs. pt
  for(int k=0; k<kNCent; k++)
    {
      c = new TCanvas(Form("TpcJpsiEff_cent%s",cent_Title[k]),Form("TpcJpsiEff_cent%s",cent_Title[k]),1100,700);
      c->Divide(4,3);
      for(int j=0; j<gNZdcRate; j++)
	{
	  if(j==0) continue;
	  c->cd(j+1);
	  hTpcJpsiEff[j][k][0]->SetMarkerStyle(20);
	  hTpcJpsiEff[j][k][0]->GetXaxis()->SetRangeUser(0.2, 10);
	  hTpcJpsiEff[j][k][0]->GetYaxis()->SetRangeUser(0, 0.3);
	  hTpcJpsiEff[j][k][0]->SetTitle(";p_{T} (GeV/c);TPC efficiency");
	  hTpcJpsiEff[j][k][0]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%d < ZDCrate < %d kHz",j*10,10+j*10),0.06);
	  t1->Draw();
	}
    }

  // compare inclusive jpsi efficiency
  double xmin = 0.15, xmax = 0.45, ymin = 0.15, ymax = 0.4;
  TLegend *leg[2];
  for(int l=0; l<2; l++)
    {
      leg[l] = new TLegend(xmin+0.3*l,ymin,xmax+0.3*l,ymax);
      leg[l]->SetBorderSize(0);
      leg[l]->SetFillColor(0);
      leg[l]->SetTextSize(0.03);
    }
  const int color[6] = {1, 2, 3, 4, 6, 7};
  TH1F *hTpcTrkEffVsCent[gNZdcRate];
  for(int j=0; j<gNZdcRate; j++)
    {
      hJpsiTpcEffVsCent[0][j][0]->GetYaxis()->SetRangeUser(0,0.1);
      hJpsiTpcEffVsCent[0][j][0]->GetXaxis()->SetLabelSize(0.05);
      hJpsiTpcEffVsCent[0][j][0]->SetMarkerStyle(20+j);
      hJpsiTpcEffVsCent[0][j][0]->SetMarkerColor(color[j%6]);
      hJpsiTpcEffVsCent[0][j][0]->SetLineColor(color[j%6]);
      leg[j/6]->AddEntry(hJpsiTpcEffVsCent[0][j][0],Form("%d < ZDC < %d kHz",j*10,j*10+10),"P");
      if(j==0) 
	{
	  for(int k=0; k<nCentBins[0]; k++)
	    {
	      hJpsiTpcEffVsCent[0][j][0]->GetXaxis()->SetBinLabel(k+1, Form("%s%%",cent_Name_npart[k]));
	    }
	  c = draw1D(hJpsiTpcEffVsCent[0][j][0]);
	  TPaveText *t1 = GetTitleText(Form("%s: TPC tracking efficiency for J/psi",run_type),0.04);
	  t1->Draw();
	}
      else     hJpsiTpcEffVsCent[0][j][0]->Draw("sames");
    }
  for(int l=0; l<2; l++) leg[l]->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Embed_JpsiTpcEffVsCentVsZdc.pdf",run_type));
  
  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int i=0; i<nPtBins_npart; i++)
	{
	  hTpcJpisEffVsCentAll[i]->Write("", TObject::kOverwrite);
	  for(int j=0; j<gNZdcRate; j++)
	    {
	      for(int t=0; t<3; t++)
		{
		  hJpsiTpcEffVsCent[i][j][t]->Write("", TObject::kOverwrite);
		}
	    }
	}

      for(int j=0; j<gNZdcRate; j++)
	{
	  for(int k=0; k<kNCent; k++)
	    {
	      for(int t=0; t<3; t++)
		{
		  hTpcJpsiEff[j][k][t]->Write("", TObject::kOverwrite);
		  hTpcJpsiEffRebin[j][k][t]->Write("", TObject::kOverwrite);
		}
	    }
	}
      hTpcJpsiEffAll->Write("", TObject::kOverwrite);
      hTpcJpsiEffRebinAll->Write("", TObject::kOverwrite);
    }
}

//================================================
void makeDataJpsiWeight(const int saveHisto)
{
  const int* nCentBins      = nCentBins_npart; 
  const int* centBinsLow    = centBins_low_npart;
  const int* centBinsHigh   = centBins_high_npart;
  const char** centName     = cent_Name_npart;
  const char** centTitle    = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  // unlike pairs in data
  TFile *fdata = TFile::Open(Form("output/%s.jpsi.root",run_type),"read");
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnJpsiInfo[3];
  for(int i=0; i<3; i++)
    {
      hnJpsiInfo[i] = (THnSparseF*)fdata->Get(Form("m%sWeight_di_mu",hName[i]));
      hnJpsiInfo[i]->Sumw2();
    }
  hnJpsiInfo[1]->Add(hnJpsiInfo[2]);

  const char* pName[2] = {"UL","LS"};
  TH1F *hJpsiInvMass[2][kNCent][gNTrgSetup][gNZdcRate];
  for(int i=0; i<2; i++)
    {
      hnJpsiInfo[i]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
      hnJpsiInfo[i]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
      for(int j=0; j<kNCent; j++)
	{
	  hnJpsiInfo[i]->GetAxis(4)->SetRange(centBins_low_npart[j],centBins_high_npart[j]);
	  for(int t=0; t<gNTrgSetup; t++)
	    {
	      if(t>0) hnJpsiInfo[i]->GetAxis(5)->SetRange(t,t);
	      for(int p=0; p<gNZdcRate; p++)
		{
		  hnJpsiInfo[i]->GetAxis(6)->SetRange(p+1,p+1);
		  hJpsiInvMass[i][j][t][p] = (TH1F*)hnJpsiInfo[i]->Projection(0);
		  hJpsiInvMass[i][j][t][p]->SetName(Form("Data_InvMass_%s_cent%s_zdc%d-%d%s",pName[i],cent_Title_npart[j],p*10,p*10+10,gTrgSetupName[t]));
		  hJpsiInvMass[i][j][t][p]->SetTitle("");
		  hnJpsiInfo[i]->GetAxis(6)->SetRange(0,-1);
		}
	      hnJpsiInfo[i]->GetAxis(5)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i]->GetAxis(4)->SetRange(0,-1);
	}
    }  

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<kNCent; j++)
	    {
	      for(int t=0; t<gNTrgSetup; t++)
		{
		  for(int p=0; p<gNZdcRate; p++)
		    {
		      hJpsiInvMass[i][j][t][p]->Write("", TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}


//================================================
void makeJpsi(const bool saveHisto)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  THnSparseF *hnJpsiInfo[nEffType][2];
  TH1F *hJpsiInvMass[nEffType][nCentBins][2];
  TH1F *hJpsiPt[nEffType][nCentBins][2];
  TH2F *hJpsiMassVsPt[nEffType][nCentBins][2];
  TH1F *hJpsiRapdity[nEffType][nCentBins][2];

  for(int i=0; i<nEffType; i++)
    {
      for(int w=0; w<1; w++)
	{
	  hnJpsiInfo[i][w] = (THnSparseF*)f->Get(Form("hJpsiInfo_%s_di_mu%s",trkEffType[i],weight_name[w]));
	  hnJpsiInfo[i][w]->GetAxis(0)->SetRangeUser(2.8, 3.4); // cut on jpsi mass
	  hnJpsiInfo[i][w]->GetAxis(2)->SetRangeUser(-1*jpsi_rapidity+0.01, jpsi_rapidity-0.01); // cut on jpsi rapidity
	  if(i>0)
	    {
	      hnJpsiInfo[i][w]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	      hnJpsiInfo[i][w]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  for(int k=0; k<nCentBins; k++)
	    {
	      hnJpsiInfo[i][w]->GetAxis(5)->SetRange(centBins_low[k],centBins_high[k]);

	      hJpsiInvMass[i][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(0);
	      hJpsiInvMass[i][k][w]->SetName(Form("hJpsiInvMass_%s_cent%s%s",trkEffType[i],cent_Title[k],weight_name[w]));
	      hJpsiInvMass[i][k][w]->SetTitle("");
	      hJpsiInvMass[i][k][w]->Sumw2();

	      hJpsiPt[i][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(1);
	      hJpsiPt[i][k][w]->SetName(Form("hJpsiPt_%s_cent%s%s",trkEffType[i],cent_Title[k],weight_name[w]));
	      hJpsiPt[i][k][w]->SetTitle("");
	      hJpsiPt[i][k][w]->SetBinContent(hJpsiPt[i][k][w]->GetNbinsX()+1,0); // reset overflow bin
	      hJpsiPt[i][k][w]->Sumw2();

	      hJpsiMassVsPt[i][k][w] = (TH2F*)hnJpsiInfo[i][w]->Projection(0,1);
	      hJpsiMassVsPt[i][k][w]->SetName(Form("hJpsiMassVsPt_%s_cent%s%s",trkEffType[i],cent_Title[k],weight_name[w]));
	      hJpsiMassVsPt[i][k][w]->SetTitle("");
	      hJpsiMassVsPt[i][k][w]->Sumw2();

	      hJpsiRapdity[i][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(2);
	      hJpsiRapdity[i][k][w]->SetName(Form("hJpsiRapdity_%s_cent%s%s",trkEffType[i],cent_Title[k],weight_name[w]));
	      hJpsiRapdity[i][k][w]->SetTitle("");
	      hJpsiRapdity[i][k][w]->Sumw2();

	      hnJpsiInfo[i][w]->GetAxis(5)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i][w]->GetAxis(0)->SetRange(0,-1);
	  hnJpsiInfo[i][w]->GetAxis(2)->SetRange(0,-1);
	  hnJpsiInfo[i][w]->GetAxis(3)->SetRange(0,-1);
	  hnJpsiInfo[i][w]->GetAxis(4)->SetRange(0,-1);
	}
    }

  // Get response matrix
  const int kNCent = nCentBins_npart[0];
  TH2F *hJpsiPtMcVsRc[kNCent][gNZdcRate];
  THnSparseF *hnJpsiMatch = (THnSparseF*)f->Get("mhJpsiMatch_di_mu");
  hnJpsiMatch->GetAxis(2)->SetRangeUser(pt2_cut+0.01,100);
  TH2F *hJpsiPtMcVsRcAll = (TH2F*)hnJpsiMatch->Projection(0,1);
  hJpsiPtMcVsRcAll->SetName(Form("Embed_JpsiPtMcVsRc_All"));
  for(int j=0; j<kNCent; j++)
    {
      hnJpsiMatch->GetAxis(3)->SetRange(centBins_low_npart[j],centBins_high_npart[j]);
      for(int p=0; p<gNZdcRate; p++)
	{
	  hnJpsiMatch->GetAxis(4)->SetRange(p+1,p+1);
	  hJpsiPtMcVsRc[j][p] = (TH2F*)hnJpsiMatch->Projection(0,1);
	  hJpsiPtMcVsRc[j][p]->SetName(Form("Embed_JpsiPtMcVsRc_cent%s_zdc%d-%d",cent_Title_npart[j],p*10,p*10+10));
	  hnJpsiMatch->GetAxis(4)->SetRange(0,-1);
	}
      hnJpsiMatch->GetAxis(3)->SetRange(0,-1);
    }
  hnJpsiInfo[0][0]->GetAxis(0)->SetRangeUser(2.8, 3.4); // cut on jpsi mass
  hnJpsiInfo[0][0]->GetAxis(2)->SetRangeUser(-1*jpsi_rapidity+0.01, jpsi_rapidity-0.01); // cut on jpsi rapidity
  TH1F *hJpsiPtMcAll = (TH1F*)hnJpsiInfo[0][0]->Projection(1);
  hJpsiPtMcAll->SetName(Form("Embed_JpsiPtMc_All"));

  TH1F *hJpsiPtMc[kNCent][gNZdcRate];
  for(int j=0; j<kNCent; j++)
    {
      hnJpsiInfo[0][0]->GetAxis(5)->SetRange(centBins_low_npart[j],centBins_high_npart[j]);
      for(int p=0; p<gNZdcRate; p++)
	{
	  hnJpsiInfo[0][0]->GetAxis(6)->SetRange(p+1,p+1);
	  hJpsiPtMc[j][p] = (TH1F*)hnJpsiInfo[0][0]->Projection(1);
	  hJpsiPtMc[j][p]->SetName(Form("Embed_JpsiPtMc_cent%s_zdc%d-%d",cent_Title_npart[j],p*10,p*10+10));
	  hnJpsiInfo[0][0]->GetAxis(6)->SetRange(0, -1);
	}
      hnJpsiInfo[0][0]->GetAxis(5)->SetRange(0, -1);
    }
  hnJpsiInfo[0][0]->GetAxis(0)->SetRange(0, -1);
  hnJpsiInfo[0][0]->GetAxis(2)->SetRange(0, -1);

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int i=0; i<nEffType; i++)
	{
	  for(int w=0; w<1; w++)
	    {
	      for(int k=0; k<nCentBins; k++)
		{
		  hJpsiMassVsPt[i][k][w]->Write("", TObject::kOverwrite);
		  hJpsiInvMass[i][k][w]->Write("", TObject::kOverwrite);
		  hJpsiPt[i][k][w]->Write("", TObject::kOverwrite);
		  hJpsiRapdity[i][k][w]->Write("", TObject::kOverwrite);
		}
	    }
	}

      hJpsiPtMcAll->SetTitle("");
      hJpsiPtMcAll->Write("", TObject::kOverwrite);
      hJpsiPtMcVsRcAll->SetTitle("");
      hJpsiPtMcVsRcAll->Write("", TObject::kOverwrite);
      for(int j=0; j<kNCent; j++)
	{
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      hJpsiPtMc[j][p]->SetTitle("");
	      hJpsiPtMc[j][p]->Write("", TObject::kOverwrite);
	      hJpsiPtMcVsRc[j][p]->SetTitle("");
	      hJpsiPtMcVsRc[j][p]->Write("", TObject::kOverwrite);
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
