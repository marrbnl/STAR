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
double getJpsiEff(const double mass, const int nExpr, const double trkEff, TH1F *hInputPt, int debug = 0);
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
const int nCentBins_npart[nPtBins_npart] = {8,7};
const int centBins_low_npart[16]  = { 15,13,11,9,7,5,3,1, 15,13,11,9,7,5,1,3 };
const int centBins_high_npart[16] = { 16,14,12,10,8,6,4,2, 16,14,12,10,8,6,4,4 };
const char *cent_Name_npart[16] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80", "0-10","10-20","20-30","30-40","40-50","50-60","60-80","60-70"};
const char *cent_Title_npart[16] = {"0010","1020","2030","3040","4050","5060","6070","7080", "0010","1020","2030","3040","4050","5060","6080","6070"};

const char *trkEffType[7] = {"MC","Tpc","MtdMth","Fake","MuonPid","MtdTrig","TrigUnit"};
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

  //makeJpsi(1);
  makeDataJpsiWeight(1);
  //makeJpsiTpcEffVsZdc(1, 1);
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

  // Get TPC efficiency vs. centrality vs. ZDCrate
  TFile *ftrk = TFile::Open(Form("Rootfiles/%s.EmbTrkEff.root",run_type),"read");
  TH1F *hMcTrkPtInZdc[2][gNZdcRate][kNCent];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<gNZdcRate; j++)
	{
	  for(int k=0; k<kNCent; k++)
	    {
	      if(i==0) hMcTrkPtInZdc[i][j][k] = (TH1F*)ftrk->Get(Form("McTrkPt_MC_cent%s_Zdc%d-%d",cent_Title[k],j*10,j*10+10));
	      if(i==1) hMcTrkPtInZdc[i][j][k] = (TH1F*)ftrk->Get(Form("McTrkPt_Tpc_cent%s_Zdc%d-%d",cent_Title[k],j*10,j*10+10));
	    }
	}
    }

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
  TH1F *hMcTrkEff[gNZdcRate];
  for(int j=0; j<gNZdcRate; j++)
    {
      hMcTrkEff[j] = new TH1F(Form("McTrkPtEff_Zdc%d",j),";;Efficiency",kNCent,0,kNCent);
      for(int k=0; k<kNCent; k++)
	{
	  hMcTrkEff[j]->GetXaxis()->SetBinLabel(k+1, Form("%s%%",cent_Name[k]));
	  double low_bin = hMcTrkPtInZdc[0][j][k]->FindFixBin(pt2_cut+1e-4);
	  double all = hMcTrkPtInZdc[0][j][k]->Integral(low_bin, -1);
	  double acc = hMcTrkPtInZdc[1][j][k]->Integral(low_bin, -1);
	  if(all>0 && acc>0)
	    {
	      hMcTrkEff[j]->SetBinContent(k+1, acc/all);
	      hMcTrkEff[j]->SetBinError(k+1, acc/all*sqrt(1/acc+1/all));
	    }
	  else
	    {
	      printf("[w] Missing %s%% ZDC%d-%d\n",cent_Name[k],j*10,j*10+10);
	    }
	}
      hMcTrkEff[j]->GetXaxis()->SetLabelSize(0.05);
      hMcTrkEff[j]->SetMarkerStyle(20+j);
      hMcTrkEff[j]->SetMarkerColor(color[j%6]);
      hMcTrkEff[j]->SetLineColor(color[j%6]);
      leg[j/6]->AddEntry(hMcTrkEff[j],Form("%d < ZDC < %d kHz",j*10,j*10+10),"P");
      if(j==0) 
	{
	  c = draw1D(hMcTrkEff[j]);
	  TPaveText *t1 = GetTitleText(Form("%s: TPC tracking efficiency",run_type),0.04);
	  t1->Draw();
	}
      else     hMcTrkEff[j]->Draw("sames");
    }
  for(int l=0; l<2; l++)
    {
      leg[l]->Draw();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Embed_TrkTpcEffVsCentVsZdc.pdf",run_type));
    }

  // calculate Jpsi efficiency due to TPC tracking
  TH1F *hMcJpsiEff[gNZdcRate];
  for(int j=0; j<gNZdcRate; j++)
    {
      hMcJpsiEff[j] = new TH1F(Form("McJpsiTpcEff_Zdc%d-%d",j*10,j*10+10),";;Efficiency",kNCent,0,kNCent);
    }
  myRandom = new TRandom3();
  myRandom->SetSeed(0);
  const int nExpr = 1e6;
  const double mass = 3.096;
  TH1F *hInputJpsi[kNCent];
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  for(int k=0; k<kNCent; k++)
    {
      if(k==0 || k==1) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
      if(k==2 || k==3) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
      if(k>=4)         hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");
      hInputJpsi[k]->SetName(Form("hInputJpsi_cent%s",cent_Name[k]));
    }

  for(int j=0; j<gNZdcRate; j++)
    {
      for(int k=0; k<kNCent; k++)
	{
	  double trkEff = hMcTrkEff[j]->GetBinContent(k+1);
	  double jpsiEff = getJpsiEff(mass, nExpr, trkEff, hInputJpsi[k], 0);
	  printf("[i] Ceentrality %s%%, %d < ZDC < %d, trk_eff = %4.2f, jpsi_eff = %4.2f\n",cent_Name[k],j*10,j*10+10,trkEff,jpsiEff);

	  hMcJpsiEff[j]->SetBinContent(k+1, jpsiEff);
	  hMcJpsiEff[j]->SetBinError(k+1, 1e-5);
	  hMcJpsiEff[j]->GetXaxis()->SetBinLabel(k+1, Form("%s%%",cent_Name[k]));
	}
    }

  for(int j=0; j<gNZdcRate; j++)
    {
      hMcJpsiEff[j]->GetXaxis()->SetLabelSize(0.05);
      hMcJpsiEff[j]->SetMarkerStyle(20+j);
      hMcJpsiEff[j]->SetMarkerColor(color[j%6]);
      hMcJpsiEff[j]->SetLineColor(color[j%6]);
      hMcJpsiEff[j]->GetYaxis()->SetRangeUser(-0.02,0.1);
      if(j==0) 
	{
	  c = draw1D(hMcJpsiEff[j]);
	  TPaveText *t1 = GetTitleText(Form("%s: TPC tracking efficiency for Jpsi",run_type),0.04);
	  t1->Draw();
	}
      else     hMcJpsiEff[j]->Draw("sames");
    }
  for(int l=0; l<2; l++)
    {
      leg[l]->Draw();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Embed_JpsiTpcEffVsCentVsZdc.pdf",run_type));
    }
  
  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int j=0; j<gNZdcRate; j++)
	{
	  hMcJpsiEff[j]->Write("", TObject::kOverwrite);
	}
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
  TH1F *hJpsiInvMass[2][kNCent][gNZdcRate];
  for(int i=0; i<2; i++)
    {
      hnJpsiInfo[i]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
      hnJpsiInfo[i]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
      for(int j=0; j<kNCent; j++)
	{
	  hnJpsiInfo[i]->GetAxis(5)->SetRange(centBins_low_npart[j],centBins_high_npart[j]);
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      hnJpsiInfo[i]->GetAxis(7)->SetRange(p+1,p+1);
	      hJpsiInvMass[i][j][p] = (TH1F*)hnJpsiInfo[i]->Projection(0);
	      hJpsiInvMass[i][j][p]->SetName(Form("Data_InvMass_%s_cent%s_zdc%d-%d",pName[i],cent_Title_npart[j],p*10,p*10+10));
	      hJpsiInvMass[i][j][p]->SetTitle("");
	      hnJpsiInfo[i]->GetAxis(7)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i]->GetAxis(5)->SetRange(0,-1);
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
	      for(int p=0; p<gNZdcRate; p++)
		{
		  hJpsiInvMass[i][j][p]->Write("", TObject::kOverwrite);
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

  THnSparseF *hnJpsiInfo[7][2];
  TH1F *hJpsiInvMass[7][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPt[7][gNTrgSetup][nCentBins][2];
  TH2F *hJpsiMassVsPt[7][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiRapdity[7][gNTrgSetup][nCentBins][2];

  for(int i=0; i<7; i++)
    {
      for(int w=0; w<2; w++)
	{
	  hnJpsiInfo[i][w] = (THnSparseF*)f->Get(Form("hJpsiInfo_%s_di_mu%s",trkEffType[i],weight_name[w]));
	  hnJpsiInfo[i][w]->GetAxis(2)->SetRangeUser(-1*jpsi_rapidity+0.01, jpsi_rapidity-0.01); // cut on jpsi rapidity
	  if(i>0)
	    {
	      hnJpsiInfo[i][w]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	      hnJpsiInfo[i][w]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      if(j>0) hnJpsiInfo[i][w]->GetAxis(6)->SetRange(j,j);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnJpsiInfo[i][w]->GetAxis(5)->SetRange(centBins_low[k],centBins_high[k]);

		  hJpsiInvMass[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(0);
		  hJpsiInvMass[i][j][k][w]->SetName(Form("hJpsiInvMass_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiInvMass[i][j][k][w]->SetTitle("");
		  hJpsiInvMass[i][j][k][w]->Sumw2();

		  hJpsiPt[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(1);
		  hJpsiPt[i][j][k][w]->SetName(Form("hJpsiPt_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiPt[i][j][k][w]->SetTitle("");
		  hJpsiPt[i][j][k][w]->SetBinContent(hJpsiPt[i][j][k][w]->GetNbinsX()+1,0); // reset overflow bin
		  hJpsiPt[i][j][k][w]->Sumw2();

		  hJpsiMassVsPt[i][j][k][w] = (TH2F*)hnJpsiInfo[i][w]->Projection(0,1);
		  hJpsiMassVsPt[i][j][k][w]->SetName(Form("hJpsiMassVsPt_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiMassVsPt[i][j][k][w]->SetTitle("");
		  hJpsiMassVsPt[i][j][k][w]->Sumw2();

		  hJpsiRapdity[i][j][k][w] = (TH1F*)hnJpsiInfo[i][w]->Projection(2);
		  hJpsiRapdity[i][j][k][w]->SetName(Form("hJpsiRapdity_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiRapdity[i][j][k][w]->SetTitle("");
		  hJpsiRapdity[i][j][k][w]->Sumw2();

		  hnJpsiInfo[i][w]->GetAxis(5)->SetRange(0,-1);
		}
	      hnJpsiInfo[i][w]->GetAxis(6)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i][w]->GetAxis(3)->SetRange(0,-1);
	  hnJpsiInfo[i][w]->GetAxis(4)->SetRange(0,-1);
	}
    }
  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int i=0; i<7; i++)
	{
	  for(int w=0; w<2; w++)
	    {
	      for(int j=0; j<gNTrgSetup; j++)
		{
		  for(int k=0; k<nCentBins; k++)
		    {
		      hJpsiMassVsPt[i][j][k][w]->Write("", TObject::kOverwrite);
		      hJpsiInvMass[i][j][k][w]->Write("", TObject::kOverwrite);
		      hJpsiPt[i][j][k][w]->Write("", TObject::kOverwrite);
		      hJpsiRapdity[i][j][k][w]->Write("", TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}


//-------------------------------------------------------
double getJpsiEff(const double mass, const int nExpr, const double trkEff, TH1F *hInputPt, int debug)
{
  double nOutJpsi = 0;
  double nInJpsi = 0;
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

      nInJpsi += 1 * weight;
      if(fabs(eta1) > 0.5 || fabs(eta2) > 0.5) continue;
      if(pt1<1.3 || pt2<1.3) continue;
      if(pt1<1.5 && pt2<1.5) continue;
      if(myRandom->Uniform(0,1)>trkEff) continue;
      if(myRandom->Uniform(0,1)>trkEff) continue;
      nOutJpsi += 1 * weight;
    }
      
  return nOutJpsi/nInJpsi;
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
