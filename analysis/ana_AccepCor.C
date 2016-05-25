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
#if (YEAR==2014)
const char *run_type = "Run14_AuAu200";
const int nPtBins = 10;
const double ptBins_low[nPtBins]  = {0,0,1,2,3,4,5,6,8,10};
const double ptBins_high[nPtBins] = {15,1,2,3,4,5,6,8,10,15};
const int gNTrgSetup = 5;
const char *gTrgSetupName[5] = {"","_P0","_P1","_P2","_P3"};
const char *gTrgSetupTitle[5] = {"","_prod","_prod_low","_prod_mid","_prod_high"};
const int nCentBins = 4;
const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60"};
const char *cent_Title[nCentBins] = {"0060","0020","2040","4060"};
#elif (YEAR==2013)
char *run_type = "Run13_pp500";
#elif (YEAR==2015)
char *run_type = "Run15_pp200";
#endif

const double pt1_cut = 1.5;
const double pt2_cut = 1.2;
const double pi = 3.1415926;
const double jpsiMass = 3.096;
const double muMass = 0.1057;
const int year = YEAR;
const int    gMtdNChannels             = 24;   // Total number of MTD channels per module. One cell has two channels.
const double gMtdCellLength            = 87.0; // Length of a MTD cell (cm)
const double gMtdCellWidth             = 3.8;  // Width of a MTD cell (cm)
const double gMtdCellGap               = 0.6;  // Gap between MTD cells (cm)
const double gMtdBacklegPhiWidth       = 8.*pi/180.;   // Width of backleg in phi (rad)
const double gMtdBacklegPhiGap         = 4.*pi/180.;   // Gap between backleg in phi (rad)
const double gMtdFirstBacklegPhiCenter = 90.*pi/180.;  // Center of backleg 1 at phi = 90 degree (rad)
const double gMtdRadius                = 403; // Minimum radius of MTD system extracted from geometry file (cm)

TRandom3 *myRandom;
void MtdAccepLoss(const int savePlot = 0, const int saveHisto = 0);
void TrkEffCentDep(const int savePlot = 0, const int saveHisto = 0);

TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
void getMtdPos(const double phi, const double z, int &backleg, int &module, int &cell);
double rotatePhi(double phi);
bool isInMtd(const int backleg, const int module, int cell);

void ScaleHistoTitle(const TH1 *h, 
                     const Double_t xTitleSize = 28, const Double_t xTitleOffset = 0.9, const Double_t xLabelSize = 20,
                     const Double_t yTitleSize = 28, const Double_t yTitleOffset = 0.9, const Double_t yLabelSize = 20,
                     const Int_t font = 42);
TPaveText *GetTitleText(TString title, const Float_t size = 0.04, const Int_t font = 62);
TPaveText *GetPaveText(double xl, double xh, double yl, double yh, double size = 0.04, const Int_t font = 42);
TCanvas *draw1D(TH1 *h, TString hTitle = "",Bool_t setLog = kFALSE, Bool_t drawP = kTRUE, const Float_t size = 0.04, const TString drawOpt = "", const Int_t titleFont = 62,
		const Int_t wh = 800, const Int_t ww = 600);
TCanvas *draw2D(TH2 *h, const TString hTitle = "" , const Float_t size = 0.04, const Bool_t logz = kTRUE, const char *drawOption = "colz");
TLine *GetLine(double xl, double yl, double xh, double yh, Color_t color=2, Width_t width=2, Style_t style=2);

TH1F *hMcJpsiPt;

//================================================
void ana_AccepCor()
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

  if(year==2013)
    {
      TFile *fWeight = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
      TF1 *funcJpsi = (TF1*)fWeight->Get("ffpt");
      funcJpsi->SetNpx(1000);
      hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
      hMcJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield_cent0100"));
      for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
	{
	  hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
	}
    }
  if(year==2014)
    {
      TFile *fWeight = TFile::Open("Rootfiles/Run14_AuAu200.Input.root","read");
      hMcJpsiPt = (TH1F*)fWeight->Get(Form("hInputJpsiShape_Cent0"));
    } 


  //MtdAccepLoss(1,1);
  TrkEffCentDep(1,1);
}

//================================================
void TrkEffCentDep(const int savePlot, const int saveHisto)
{
  gStyle->SetStatY(0.6);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  // extract single muon efficiency for each centrality
  // Jpsi coutns as weights
  TFile *fYield = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiCounts[gNTrgSetup-1];
  double nJpsi[nCentBins-1][gNTrgSetup-1] = {0};
  double nJpsiCent[nCentBins-1] = {0};
  double nJpsiAll = 0;
  for(int i=0; i<gNTrgSetup-1; i++)
    {
      hJpsiCounts[i] = (TH1F*)fYield->Get(Form("NJpsiInCent_weight%s",gTrgSetupName[i+1]));
      for(int bin=1; bin<=hJpsiCounts[i]->GetNbinsX(); bin++)
	{
	  nJpsi[bin-1][i] = hJpsiCounts[i]->GetBinContent(bin);
	  nJpsiCent[bin-1] += hJpsiCounts[i]->GetBinContent(bin);
	  nJpsiAll += hJpsiCounts[i]->GetBinContent(bin);
	  printf("[i] %s%% %s: %4.2f Jpsi\n",cent_Title[bin],gTrgSetupTitle[i+1],nJpsi[bin-1][i]);
	}
    }


  // tracking efficiency
  const int nHistos = nCentBins + 6;
  const char *centTitle[nHistos] = {"0060","0020","2040","4060","010","1020","2030","3040","4050","5060"};
  const char *trkEffType[2] = {"MC","Tpc"};
  TFile *fTrkEff = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"read");
  TH1F *hMcTrkPt[2][gNTrgSetup-1][nHistos];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  for(int k=0; k<nHistos; k++)
	    {
	      TH1F *htmp= (TH1F*)fTrkEff->Get(Form("hMcTrkPt_%s_cent%s%s",trkEffType[i],centTitle[k],gTrgSetupTitle[j+1]));
	      hMcTrkPt[i][j][k] = new TH1F(Form("hMcTrkPt_%s_cent%s%s_clone",trkEffType[i],centTitle[k],gTrgSetupTitle[j+1]),"",(int)(20/htmp->GetBinWidth(1)),0,20);
	      for(int bin=1; bin<=hMcTrkPt[i][j][k]->GetNbinsX(); bin++)
		{
		  double jbin = htmp->FindFixBin(hMcTrkPt[i][j][k]->GetBinCenter(bin));
		  hMcTrkPt[i][j][k]->SetBinContent(bin,htmp->GetBinContent(jbin));
		  hMcTrkPt[i][j][k]->SetBinError(bin,TMath::Sqrt(htmp->GetBinContent(jbin)));
		}
	      hMcTrkPt[i][j][k]->Rebin(2);
	    }
	}
    }
  
  // trigger efficiency
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/Run14.AuAu200.Input.root"),"read");
  TH1F *hTrigEff[gNTrgSetup-1][nCentBins-1];
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      for(int k=0; k<nCentBins-1; k++)
	{ 
	  TH1F *htmp = (TH1F*)fTrigEff->Get(Form("MuonTrigEff_Cent%d%s",k,gTrgSetupName[j+1]));
	  cout << htmp->GetName() << endl;
	  hTrigEff[j][k] = new TH1F(Form("MuonTrigEff_Cent%s%s",cent_Title[k+1],gTrgSetupName[j+1]),"",(int)(30/htmp->GetBinWidth(1)),0,30);
	  for(int bin=1; bin<=hTrigEff[j][k]->GetNbinsX(); bin++)
	    {
	      double x = hTrigEff[j][k]->GetBinCenter(bin);
	      if(x>9) x = htmp->GetBinCenter(htmp->GetNbinsX());
	      double jbin = htmp->FindFixBin(x);
	      hTrigEff[j][k]->SetBinContent(bin,htmp->GetBinContent(jbin));
	      hTrigEff[j][k]->SetBinError(bin, 1e-10);
	    }
	}
    }

  // 1/eff efficiency
  TH1F *hTrkOneOverEff[gNTrgSetup-1][nHistos];
  for(int j=0; j<gNTrgSetup-1; j++)
    {
      for(int k=0; k<nHistos; k++)
	{
	  hTrkOneOverEff[j][k] = (TH1F*)hMcTrkPt[0][j][k]->Clone(Form("hTrkOneOverEff_cent%s%s",centTitle[k],gTrgSetupTitle[j+1]));
	  hTrkOneOverEff[j][k]->Divide(hMcTrkPt[1][j][k]);
	  int index = k-1;
	  if(k==0) index = 0;
	  if(k>=nCentBins) index = (k-4)/2;
	  for(int bin=1; bin<=hTrkOneOverEff[j][k]->GetNbinsX(); bin++)
	    {
	      double x = hTrkOneOverEff[j][k]->GetBinCenter(bin);
	      double jbin = hTrigEff[j][index]->FindFixBin(x);
	      double eff = hTrigEff[j][index]->GetBinContent(jbin);
	      hTrkOneOverEff[j][k]->SetBinContent(bin,hTrkOneOverEff[j][k]->GetBinContent(bin)/eff);
	      hTrkOneOverEff[j][k]->SetBinError(bin,hTrkOneOverEff[j][k]->GetBinError(bin)/eff);
	    }
	}
    }
  
  // final avergae efficiency
  TH1F *hTrkEff[nHistos];
  for(int k=0; k<nHistos; k++)
    {
      hTrkEff[k] = (TH1F*)hMcTrkPt[0][0][k]->Clone(Form("hTrkEff_cent%s",centTitle[k]));
      hTrkEff[k]->Reset();
    }
  //0-60%
  for(int k=1; k<nCentBins; k++)
    { 
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  hTrkEff[0]->Add(hTrkOneOverEff[j][k], nJpsi[k-1][j]/nJpsiAll);
	}
    }
  // other centralities
  for(int k=1; k<nHistos; k++)
    { 
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  int index = k-1;
	  if(k>=nCentBins) index = (k-4)/2;
	  hTrkEff[k]->Add(hTrkOneOverEff[j][k], nJpsi[index][j]/nJpsiCent[index]);
	}
    }

  for(int k=0; k<nHistos; k++)
    { 
      for(int bin=1; bin<=hTrkEff[k]->GetNbinsX(); bin++)
	{
	  double value = hTrkEff[k]->GetBinContent(bin);
	  double error = hTrkEff[k]->GetBinError(bin);
	  hTrkEff[k]->SetBinContent(bin, 1./value);
	  hTrkEff[k]->SetBinError(bin, 1./value * error/value);
	}
    }

  TLegend *leg = new TLegend(0.3,0.2,0.5,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  for(int k=0; k<nCentBins; k++)
    {
      hTrkEff[k]->SetMarkerStyle(20+k);
      hTrkEff[k]->SetMarkerColor(TMath::Power(2,k));
      hTrkEff[k]->SetTitle("TPC tracking efficiency");
      hTrkEff[k]->GetYaxis()->SetRangeUser(0,0.8);
      leg->AddEntry(hTrkEff[k],Form("%s%%",cent_Name[k]),"PL");
    }
  TCanvas *c = draw1D(hTrkEff[0]);
  for(int k=1; k<nCentBins; k++)
    {
      hTrkEff[k]->Draw("samesP");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompTrkEff_InCent.pdf",run_type));

  c = new TCanvas("fit_trkeff","fit_trkeff",1200,600);
  c->Divide(4,3);
  TF1 *func[nHistos];
  for(int k=0; k<nHistos; k++)
    {
      TH1F *hfit = (TH1F*)hTrkEff[k]->Clone(Form("Fit_%s",hTrkEff[k]->GetName()));
      hfit->GetXaxis()->SetRangeUser(0.5,20);
      hfit->GetYaxis()->SetRangeUser(0,0.8);
      hfit->SetMarkerStyle(25);
      hfit->SetMarkerColor(1);
      func[k] = new TF1(Form("func_%d",k),"[0]*exp(-pow([1]/x,[2])-pow([3]/x/x,[4])-pow([5]/x/x/x,[6]))",1.2,20);
      func[k]->SetParameters(1,0.1,1,0.1,1,0.1,1);
      hfit->Fit(func[k],"RQ0");
      c->cd(k+1);
      hfit->SetTitle("");
      hfit->Draw();
      func[k]->SetLineColor(4);
      func[k]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency (%s%%)",centTitle[k]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitTrkEff_InCent.pdf",run_type));

  // Toy MC
  TH1F *hInJpsiPt;
  TH1F *hOutJpsiPt[nHistos];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  hInJpsiPt = new TH1F(Form("hInJpsiPt"),"p_{T} distribution of input J/#Psi;p_{T} (GeV/c)",nbins,xbins);
  hInJpsiPt->Sumw2();
  for(int i=0; i<nHistos; i++)
    {
      hOutJpsiPt[i] = new TH1F(Form("hOutJpsiPt_%d",i),"p_{T} distribution of reconstructed J/#Psi;p_{T} (GeV/c)",nbins,xbins);
      hOutJpsiPt[i]->Sumw2();
    }


  const int nExpr = 5e7;
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  =  myRandom->Uniform(0,20);
      double weight = hMcJpsiPt->GetBinContent(hMcJpsiPt->FindFixBin(mc_pt));
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+jpsiMass*jpsiMass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,jpsiMass);
      hInJpsiPt->Fill(mc_pt,weight);

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      double z1 = TMath::Tan(pi/2-daughter1.Theta()) * gMtdRadius;

      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      double z2 = TMath::Tan(pi/2-daughter2.Theta()) * gMtdRadius;

      double leadpt = pt1 > pt2 ? pt1 : pt2;
      double subpt  = pt1 < pt2 ? pt1 : pt2;
      if(leadpt<pt1_cut || subpt<pt2_cut) continue;
      
      for(int k=0; k<nHistos; k++)
	{
	  double eff1 = func[k]->Eval(pt1);
	  double eff2 = func[k]->Eval(pt2);
	  double pro1 = myRandom->Uniform(0., 1.);
	  double pro2 = myRandom->Uniform(0., 1.);
	  if(pro1<eff1 && pro2<eff2)
	    {
	      hOutJpsiPt[k]->Fill(mc_pt,weight);
	    }
	}
    }

  TH1F *hJpsiEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEff[k] = (TH1F*)hOutJpsiPt[k]->Clone(Form("hJpsiEffCorr_cent%s",centTitle[k]));
      hJpsiEff[k]->Divide(hOutJpsiPt[1]);
      hJpsiEff[k]->SetMarkerStyle(20+k);
      hJpsiEff[k]->SetMarkerColor(TMath::Power(2,k));
      hJpsiEff[k]->SetTitle("Ratio of J/#psi efficiency");
      hJpsiEff[k]->GetYaxis()->SetRangeUser(0.8,1.2);
    }
  c = draw1D(hJpsiEff[0]);
  hJpsiEff[2]->Draw("sames");
  hJpsiEff[3]->Draw("sames");
  leg = new TLegend(0.5,0.65,0.7,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  leg->AddEntry(hTrkEff[0],Form("%s%%/%s%%",cent_Name[0],cent_Name[1]),"PL");
  leg->AddEntry(hTrkEff[2],Form("%s%%/%s%%",cent_Name[2],cent_Name[1]),"PL");
  leg->AddEntry(hTrkEff[3],Form("%s%%/%s%%",cent_Name[3],cent_Name[1]),"PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffCorr_ToCentral.pdf",run_type));

  // Npart analysis
  int ptCuts[2] = {0,5};
  TH1F *hJpsiEffCent[2];
  for(int i=0; i<2; i++)
    {
      hJpsiEffCent[i] = new TH1F(Form("hJpsiEffCor_pt%d-15",ptCuts[i]),"",6,0,6);
      for(int bin=1; bin<=6; bin++)
	{
	  int low_bin = hOutJpsiPt[4]->FindFixBin(ptCuts[i]+1e-4);
	  int high_bin = hOutJpsiPt[4]->FindFixBin(15-1e-4);
	  double denom_err;
	  double denom = hOutJpsiPt[4]->IntegralAndError(low_bin,high_bin,denom_err);
	  double numer_err;
	  double numer =  hOutJpsiPt[10-bin]->IntegralAndError(low_bin,high_bin,numer_err);
	  double value = numer/denom;
	  double err = value*TMath::Sqrt(numer_err*numer_err/numer/numer+denom_err*denom_err/denom/denom);
	  hJpsiEffCent[i]->SetBinContent(bin,value);
	  hJpsiEffCent[i]->SetBinError(bin,err);
	}
      hJpsiEffCent[i]->SetMarkerStyle(20+i);
      hJpsiEffCent[i]->SetMarkerColor(TMath::Power(2,i));
      hJpsiEffCent[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    }
  c = draw1D(hJpsiEffCent[0],"Ratio of J/#psi efficiency to 0-10%");
  hJpsiEffCent[1]->Draw("samesP");

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiEff[k]->SetTitle("");
	  hJpsiEff[k]->Write("",TObject::kOverwrite);
	}
      fout->Close();

      fout =  TFile::Open(Form("Rootfiles/%s.Npart.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update"); 
      for(int i=0; i<2; i++)
	{
	  hJpsiEffCent[i]->Write("",TObject::kOverwrite);
	}
    }
  
}

//================================================
void MtdAccepLoss(const int savePlot, const int saveHisto)
{
  const int nHistos = 4;
  TH1F *hInJpsiPt;
  TH1F *hOutJpsiPt[nHistos];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  hInJpsiPt = new TH1F(Form("hInJpsiPt"),"p_{T} distribution of input J/#Psi;p_{T} (GeV/c)",nbins,xbins);
  hInJpsiPt->Sumw2();
  for(int i=0; i<nHistos; i++)
    {
      hOutJpsiPt[i] = new TH1F(Form("hOutJpsiPt_%d",i),"p_{T} distribution of reconstructed J/#Psi;p_{T} (GeV/c)",nbins,xbins);
      hOutJpsiPt[i]->Sumw2();
    }

  TH1F *hRespEff[150];
  TFile *fRespEff = TFile::Open(Form("Rootfiles/Run%dResponseEffViaPtTemplate.root",year-2000),"read");
  for(int i=0; i<150; i++)
    {
      int backleg = i/5 + 1;
      int module = i%5 + 1;
      TF1 *func  = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",i/5,i%5));
      hRespEff[i] = new TH1F(Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),1500,0,30);
      for(int bin=1; bin<=hRespEff[i]->GetNbinsX(); bin++)
        {
          double x = hRespEff[i]->GetBinCenter(bin);
	  hRespEff[i]->SetBinContent(bin,func->Eval(x));
	  hRespEff[i]->SetBinError(bin,1e-10);
        }
    }

  const int nExpr = 5e8;
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  =  myRandom->Uniform(0,20);
      double weight = hMcJpsiPt->GetBinContent(hMcJpsiPt->FindFixBin(mc_pt));
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+jpsiMass*jpsiMass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,jpsiMass);
      hInJpsiPt->Fill(mc_pt,weight);

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      double z1 = TMath::Tan(pi/2-daughter1.Theta()) * gMtdRadius;

      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      double z2 = TMath::Tan(pi/2-daughter2.Theta()) * gMtdRadius;

      double leadpt = pt1 > pt2 ? pt1 : pt2;
      double subpt  = pt1 < pt2 ? pt1 : pt2;
      if(leadpt<pt1_cut || subpt<pt2_cut) continue;

      int backleg1, module1, cell1;
      getMtdPos(phi1, z1, backleg1, module1, cell1);
      bool isMtd1 = isInMtd(backleg1, module1, cell1);

      int backleg2, module2, cell2;
      getMtdPos(phi2, z2, backleg2, module2, cell2);
      bool isMtd2 = isInMtd(backleg2, module2, cell2);

      if(!isMtd1 || !isMtd2) continue;

      bool isAcc1 = false;
      int index1 = (backleg1-1)*5+module1-1;
      int bin1 = hRespEff[index1]->FindBin(pt1);
      double eff1 = hRespEff[index1]->GetBinContent(bin1);
      double prob1 = myRandom->Uniform(0,1);
      if(prob1<eff1) isAcc1 = true;

      bool isAcc2 = false;
      int index2 = (backleg2-1)*5 + module2-1;
      int bin2 = hRespEff[index2]->FindFixBin(pt2);
      double eff2 = hRespEff[index2]->GetBinContent(bin2);
      double prob2 = myRandom->Uniform(0,1);
      if(prob2<eff2) isAcc2 = true;
      if(!isAcc1 || !isAcc2) continue;

      hOutJpsiPt[0]->Fill(mc_pt,weight);

      if( !(backleg1==8 || backleg1==19 || (backleg1==15&&module1==4)) &&
	  !(backleg2==8 || backleg2==19 || (backleg2==15&&module2==4)))
	{
	  hOutJpsiPt[1]->Fill(mc_pt,weight);
	}


      if( !(backleg1==8 || (backleg1==15&&module1==4)) &&
	  !(backleg2==8 || (backleg2==15&&module2==4)))
	{
	  hOutJpsiPt[2]->Fill(mc_pt,weight);
	}

      if( !( (backleg1==8&&module1<=2) || (backleg1==15&&module1==4)) &&
	  !( (backleg2==8&&module2<=2) || (backleg2==15&&module2==4)))
	{
	  hOutJpsiPt[3]->Fill(mc_pt,weight);
	}
    }


  const char *runName[nHistos-1] = {"Run<15098067","15098068<=Run<15106131","Run>=15106131"};
  TH1F *hAcc[nHistos-1];
  for(int i=0; i<nHistos-1; i++)
    {
      hAcc[i] = (TH1F*)hOutJpsiPt[i+1]->Clone(Form("hAccepLoss_%d",i));
      hAcc[i]->Divide(hOutJpsiPt[0]);
      hAcc[i]->SetMarkerStyle(21);
      hAcc[i]->GetYaxis()->SetRangeUser(0.85,1.05);
      TCanvas *c = draw1D(hAcc[i],Form("%s: MTD acceptance loss for J/#Psi (%s)",run_type,runName[i]));
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_AccepCor/JpsiEff_Acc%d.pdf",run_type,i));
    }
 
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"recreate");
      for(int i=0; i<nHistos-1; i++)
	{
	  hAcc[i]->Write("",TObject::kOverwrite);
	}
    }
}


//_____________________________________________________________________________
bool isInMtd(const int backleg, const int module, int cell)
{
  if(backleg<1 || backleg>30) return false;
  if(module<1 || module>5) return false;
  if(cell<0 || cell>11) return false;
  
  if(year==2013)
    {
      if((backleg>=8 && backleg<=9) || 
	 (backleg>=11 && backleg<=21) ||
	 (backleg>=23 && backleg<=24) ||
	 (backleg==7 && module==5))
	return false;
    }
  else if(year==2014)
    {
      if((backleg==9 || backleg==23) ||
	 (backleg>=12 && backleg<=20 && (module==1 || module==5)))
	return false;
    }
  else if(year==2015)
    {
      if((backleg==9 || backleg==23) ||
	 (backleg>=12 && backleg<=20 && (module==1 || module==5)))
	return false;
    }

 
  return true;
}

//_____________________________________________________________________________
void getMtdPos(const double phi, const double z, int &backleg, int &module, int &cell)
{
  backleg = -1; module = -1; cell = -1;

  double phiTemp = rotatePhi(phi);
  backleg = (Int_t)(phiTemp/(gMtdBacklegPhiWidth+gMtdBacklegPhiGap));
  backleg += 24;
  if(backleg>30) backleg -= 30;

  double temp = (z+2.5*gMtdCellLength)/gMtdCellLength;
  if(temp>0) module = (Int_t)temp + 1;
  if(module<1 || module>5) module = -1;

  double lowEdge = gMtdFirstBacklegPhiCenter + (backleg-1)*(gMtdBacklegPhiWidth+gMtdBacklegPhiGap) - (gMtdNChannels/4.)*(gMtdCellWidth+gMtdCellGap)/gMtdRadius; //approximation
  lowEdge = rotatePhi(lowEdge);
  double cellPhi = phi - lowEdge;
  cellPhi = rotatePhi(cellPhi);
  cell = (Int_t) ( cellPhi/((gMtdCellWidth+gMtdCellGap)/gMtdRadius));
  if(module>3) cell = 11 - cell;
  if(cell<0 || cell>12) cell = -1;
}

//_____________________________________________________________________________
double rotatePhi(double phi) 
{
  double outPhi = phi;
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
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
