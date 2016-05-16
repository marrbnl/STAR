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
void makeHisto(const int savePlot = 0, const int saveHisto = 0);
void anaSys(const int savePlot = 0);
void singleMuon(const int savePlot = 0);
void toyMC(const double mass, const int nExpr, const int debug,
	   const int inputType, TH1F *hJpsiTruth, TH1F *hRespEff[150], 
	   TH1F *hInputJpsiPt, TH1F *hMtdJpsiPt, TH1F *hAccJpsiPt,
	   TH1F *hMuonMap = 0x0, TH1F *hMuonMapTriggered = 0x0);

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

TH1F *hInputJpsiEta;
TH1F *hInputJpsiPhi;
TH1F *hInputJpsiY;
TH2F *hJpsiPtVsOpenAngle;
TH1F *hAccMuonMap;

//================================================
void ana_MtdRespEff()
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
  hAccMuonMap = 0x0;


  //anaSys(0);
  makeHisto(1,1);
  //singleMuon(0);
}

//================================================
void makeHisto(const int savePlot, const int saveHisto)
{
  const int nPtBins = 20;
  const double minPt = 0;
  const double maxPt = 10;

  // general QA plots
  hInputJpsiEta = new TH1F(Form("hInputJpsiEta"),"#eta distribution of input J/psi;#eta",100,-1.5,1.5);
  hInputJpsiPhi = new TH1F(Form("hInputJpsiPhi"),"#varphi distribution of input J/psi;#varphi",100,-1*pi,pi);
  hInputJpsiY = new TH1F(Form("hInputJpsiY"),"y distribution of input J/psi;y",100,-1.5,1.5);
  hAccMuonMap = new TH1F(Form("hAccMuonMap"),"Acceptance muon multiplicity;(backleg-1)*5+module",150,0.5,150.5);
  hJpsiPtVsOpenAngle= new TH2F(Form("hJpsiPtVsOpenAngle"),"Opening angle between muon daughters vs J/psi p_{T};p_{T} (GeV/c);#Delta#varphi",20,0,10,100,0,pi);

  const int mode = 0; // 0 - module wise; 1 - module average; 2 - stat. err.; 3 - systematics
  const char *mode_name[4] = {"PerMod","AvgMod","StatPerMod","SysPerMod"};
  printf("[i] Process mode %d: %s\n",mode,mode_name[mode]);
  TH1F *hInputJpsiPt = new TH1F(Form("hInputJpsiPt_%s",mode_name[mode]),"p_{T} distribution of input J/psi;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
  TH1F *hMtdJpsiPt = new TH1F(Form("hMtdJpsiPt_%s",mode_name[mode]),"p_{T} distribution of J/psi in MTD acceptance;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
  TH1F *hAccJpsiPt = new TH1F(Form("hAccJpsiPt_%s",mode_name[mode]),"p_{T} distribution of J/psi with MTD responsce;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
  TH1F *hMuonMap = new TH1F(Form("hMuonMap_%s",mode_name[mode]),"Muon multiplicity;(backleg-1)*5+module",150,0.5,150.5);
  TH1F *hMuonMapTriggered = new TH1F(Form("hMuonMapTriggered_%s",mode_name[mode]),"Muon multiplicity (Dimuon trigger);(backleg-1)*5+module",150,0.5,150.5);
  hInputJpsiPt->Sumw2();
  hMtdJpsiPt->Sumw2();
  hAccJpsiPt->Sumw2();
  hMuonMap->Sumw2();
  hMuonMapTriggered->Sumw2();
  TH1F *hJpsiEff = 0x0;

  // input spectrum shape
  TH1F *hMcJpsiPt;
  if(year==2013)
    {
      TFile *fWeight = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
      TF1 *funcJpsi = (TF1*)fWeight->Get("ffpt");
      funcJpsi->SetNpx(1000);
      hMcJpsiPt = (TH1F*)funcJpsi->GetHistogram();
      hMcJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield_cent00100"));
      for(int bin=1; bin<=hMcJpsiPt->GetNbinsX(); bin++)
	{
	  hMcJpsiPt->SetBinContent(bin,hMcJpsiPt->GetBinCenter(bin)*hMcJpsiPt->GetBinContent(bin));
	}
    }
  if(year==2014)
    {
      TFile *fWeight = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");
      hMcJpsiPt = (TH1F*)fWeight->Get(Form("TBW_Jpsi_Yield_cent0060"));
    } 


  // module wise response efficiency
  TH1F *hRespEff[150];
  TFile *fRespEff = TFile::Open(Form("Rootfiles/Run%dResponseEffViaPtTemplate.root",year-2000),"read");
  for(int i=0; i<150; i++)
    {
      int backleg = i/5 + 1;
      int module = i%5 + 1;
      TF1 *func  = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",i/5,i%5));
      TF1 *func1 = (TF1*)fRespEff->Get(Form("fPtMtdEffBkl%d_Mod%d",i/5,i%5));
      hRespEff[i] = new TH1F(Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),1000,0,20);
      for(int bin=1; bin<=hRespEff[i]->GetNbinsX(); bin++)
	{
	  double x = hRespEff[i]->GetBinCenter(bin);
	  if(year==2013)
	    {
	      hRespEff[i]->SetBinContent(bin,func->Eval(x));
	      hRespEff[i]->SetBinError(bin,1e-10);
	    }
	  else
	    {		      
	      if(i/5>=9 && i/5<=21) hRespEff[i]->SetBinContent(bin,func1->Eval(x));
	      else                  hRespEff[i]->SetBinContent(bin,func->Eval(x));
	    }
	}
    }

  
  if(mode==0)
    {
      toyMC(jpsiMass, 1e8, 0, 0, hMcJpsiPt, hRespEff, hInputJpsiPt, hMtdJpsiPt, hAccJpsiPt, hMuonMap, hMuonMapTriggered);
      
      hInputJpsiPt->SetMarkerStyle(20);
      TCanvas *c = new TCanvas("Input_Jpsi","Input_Jpsi",1100,700);
      c->Divide(2,2);
      c->cd(1); hInputJpsiPt->Draw();
      c->cd(2); hInputJpsiEta->Draw();
      c->cd(3); hInputJpsiPhi->SetMinimum(0); hInputJpsiPhi->Draw();
      c->cd(4); hInputJpsiY->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_InputJpsi.pdf",run_type,mode_name[mode]));

      c = draw2D(hJpsiPtVsOpenAngle);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_OpenAngleVsJpsiPt.pdf",run_type,mode_name[mode]));

      hJpsiEff = (TH1F*)hAccJpsiPt->Clone(Form("JpsiEffVsPt_%s",mode_name[mode]));
      hJpsiEff->Divide(hMtdJpsiPt);
      hJpsiEff->SetMarkerStyle(21);
      hJpsiEff->GetYaxis()->SetRangeUser(0,1);
      c = draw1D(hJpsiEff,"MTD response efficiency for J/psi");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_JpsiEffVsPt.pdf",run_type,mode_name[mode]));

      TH1F *hMuonEff = (TH1F*)hMuonMap->Clone("hMuonEff");
      hMuonEff->Sumw2();
      hMuonEff->Divide(hAccMuonMap);
      hMuonEff->SetMarkerStyle(21);
      hMuonEff->GetYaxis()->SetRangeUser(0,1);
      c = draw1D(hMuonEff,"MTD response efficiency for single muon");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_MuonEffVsPt.pdf",run_type,mode_name[mode]));

      hMuonMapTriggered->Scale(1./hMuonMapTriggered->Integral());
      hMuonMapTriggered->SetLineColor(2);
      hMuonMapTriggered->SetMaximum(1.5*hMuonMapTriggered->GetMaximum());
      c = draw1D(hMuonMapTriggered,"Normalized muon multiplicity in each module",false,false);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_CheckMuonEffInMod.pdf",run_type,mode_name[mode]));

      hMuonMap->Scale(1./hMuonMap->Integral());
      hMuonMap->Draw("sames HIST");
      TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hMuonMap,"Single-muon trigger","L");
      leg->AddEntry(hMuonMapTriggered,"Di-muon trigger","L");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_MuonMultiplicity.pdf",run_type,mode_name[mode]));

      if(saveHisto)
	{
	  TFile *fout = 0x0;
	  if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","update");
	  if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root","update");
	  if(year==2015) fout = TFile::Open("Rootfiles/Run15.pp200.MtdResponseEff.root","update");
	  for(int i=0; i<150; i++)
	    {
	      hRespEff[i]->Write("",TObject::kOverwrite);
	    }
	  hMuonMap->Write("",TObject::kOverwrite);
	  hMuonMapTriggered->Write("",TObject::kOverwrite);
	  hJpsiEff->Write("",TObject::kOverwrite);
	}
    }
  else if(mode==1)
    {
      TFile *fin = 0x0;
      if(year==2013) fin = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","read");

      TH1F *hWeight = (TH1F*)fin->Get("hMuonMapTriggered_PerMod");
      TH1F *hRespEffAvg[150];
      for(int i=0; i<150; i++)
	{
	  int backleg = i/5 + 1;
	  int module = i%5 + 1;
	  hRespEffAvg[i] = (TH1F*)hRespEff[i]->Clone(Form("MtdRespAvgEffvsPt_Bkl%d_Mod%d",backleg,module));
	  hRespEffAvg[i]->Reset();
	}

      for(int bin=1; bin<=hRespEffAvg[0]->GetNbinsX(); bin++)
	{
	  double sumn = 0;
	  double sumN = 0;
	  for(int i=0; i<150; i++)
	    {
	      double weight = hWeight->GetBinContent(i+1);
	      double eff = hRespEff[i]->GetBinContent(bin);
	      sumn += weight;
	      sumN += eff<1e-5? 0 : weight / eff;
	    }
	  double value = sumN<1e-5 ? 0 : sumn/sumN;
	  
	  for(int i=0; i<150; i++)
	    {
	      hRespEffAvg[i]->SetBinContent(bin,value);
	      hRespEffAvg[i]->SetBinError(bin,1e-10);
	    }
	}
      toyMC(jpsiMass, 1e6, 0, 0, hMcJpsiPt, hRespEffAvg, hInputJpsiPt, hMtdJpsiPt, hAccJpsiPt, hMuonMap, hMuonMapTriggered);

      TCanvas *c = draw1D(hRespEffAvg[0]);

      hJpsiEff = (TH1F*)hAccJpsiPt->Clone(Form("JpsiEffVsPt_%s",mode_name[mode]));
      hJpsiEff->Divide(hMtdJpsiPt);
      hJpsiEff->SetMarkerStyle(21);
      hJpsiEff->GetYaxis()->SetRangeUser(0,1);
      c = draw1D(hJpsiEff,"MTD response efficiency for J/psi");
      TH1F *hJpsiEffPerMod = (TH1F*)fin->Get("JpsiEffVsPt_PerMod");
      hJpsiEffPerMod->SetLineColor(2);
      hJpsiEffPerMod->SetMarkerColor(2);
      hJpsiEffPerMod->Draw("sames");
      TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hJpsiEff,"Module-average efficiency","L");
      leg->AddEntry(hJpsiEffPerMod,"Module-wise efficiency","L");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_CompJpsiEffVsPt.pdf",run_type,mode_name[mode]));

      TH1F *hJpsiEffRatio = (TH1F*)hJpsiEff->Clone("hJpsiEffRatio");
      hJpsiEffRatio->GetYaxis()->SetRangeUser(0.5,1.5);
      hJpsiEffRatio->Divide(hJpsiEffPerMod);
      c = draw1D(hJpsiEffRatio,"Ratio of response efficiency: Module-average/Module-wise");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_RatioJpsiEff.pdf",run_type,mode_name[mode]));
    }
  else if(mode==2)
    {
      // statistical errors of scale factors
      
      gStyle->SetOptFit(0);
      TH1F *hSysMtdRespEff = new TH1F(Form("hSysMtdRespEff_1"),"Systematic uncertainty of MTD response efficiency;J/#Psi p_{T} (GeV/c)",nPtBins,minPt,maxPt);

      // statistical error on efficiency
      TH1F *hSclFacErr = (TH1F*)fRespEff->Get("hSclFacErr");

      const int nIter = 50;
      TH1F *hSysInputJpsi[nIter];
      TH1F *hSysMtdJpsiPt[nIter];
      TH1F *hSysAccJpsiPt[nIter];
      TH1F *hSysJpsiPtEff[nIter];

      TH1F *hRespEffTmp[150];
      TCanvas *c = 0x0;
      for(int ih = 0; ih<nIter; ih++)
	{
	  hSysInputJpsi[ih] = new TH1F(Form("hSysInputJpsi_%d",ih),"p_{T} distribution of input J/psi;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  hSysMtdJpsiPt[ih] = new TH1F(Form("hSysMtdJpsiPt_%d",ih),"p_{T} distribution of J/psi in MTD;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  hSysAccJpsiPt[ih] = new TH1F(Form("hSysAccJpsiPt_%d",ih),"p_{T} distribution of triggered J/psi;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  for(int i=0; i<150; i++)
	    {
	      hRespEffTmp[i] = 0x0;
	      hRespEffTmp[i] = (TH1F*)hRespEff[i]->Clone(Form("hRespEffTmp_%d_%d",i,ih));
	      double error = hSclFacErr->GetBinContent(i+1);
	      double scale = myRandom->Gaus(1, error);
	      hRespEffTmp[i]->Scale(scale);
	      //printf("BL = %d, Mod = %d, err = %4.4f, scale = %4.4f\n",i/5+1,i%5+1,error,scale);
	    }
	  toyMC(jpsiMass, 1e7, 0, 1, hMcJpsiPt, hRespEffTmp, hSysInputJpsi[ih], hSysMtdJpsiPt[ih], hSysAccJpsiPt[ih]);
	  hSysJpsiPtEff[ih] = (TH1F*)hSysAccJpsiPt[ih]->Clone(Form("hSysJpsiPtEff_%d",ih));
	  hSysJpsiPtEff[ih]->Sumw2();
	  hSysJpsiPtEff[ih]->Divide(hSysMtdJpsiPt[ih]);
	  hSysJpsiPtEff[ih]->SetMarkerStyle(20);
	  hSysJpsiPtEff[ih]->SetMarkerColor(ih+1);
	  hSysJpsiPtEff[ih]->SetLineColor(ih+1);
	  hSysJpsiPtEff[ih]->GetYaxis()->SetRangeUser(0,1);
	  if(ih==0) c = draw1D(hSysJpsiPtEff[ih],"MTD response efficiency for J/psi");
	  else hSysJpsiPtEff[ih]->Draw("sames");
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_JpsiEff_Iterations1.pdf",run_type));

      TH1F *hJpsiEffPtBin[nPtBins];
      TF1 *funcJpsiEffPtBin[nPtBins];
      TCanvas *cFit = new TCanvas("cFit","cFit",1100,700);
      cFit->Divide(5,4);
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiEffPtBin[i] = new TH1F(Form("hJpsiEffPtBin%d",i+1),"",200,0,1);
	  hJpsiEffPtBin[i]->Sumw2();
	  for(int ih = 0; ih<nIter; ih++)
	    {
	      hJpsiEffPtBin[i]->Fill(hSysJpsiPtEff[ih]->GetBinContent(i+1));
	    }
	  double init_mean = hJpsiEffPtBin[i]->GetMean();
	  funcJpsiEffPtBin[i] = new TF1(Form("funcJpsiEffPtBin%d",i+1),"gaus",init_mean-0.4,init_mean+0.4);
	  hJpsiEffPtBin[i]->Fit(funcJpsiEffPtBin[i],"0IR");
	  double mean = funcJpsiEffPtBin[i]->GetParameter(1);
	  double sigma = funcJpsiEffPtBin[i]->GetParameter(2);
	  cFit->cd(i+1);
	  hJpsiEffPtBin[i]->SetMarkerStyle(25);
	  hJpsiEffPtBin[i]->GetXaxis()->SetRangeUser(mean-5*sigma,mean+5*sigma);
	  hJpsiEffPtBin[i]->SetMaximum(1.5*hJpsiEffPtBin[i]->GetMaximum());
	  hJpsiEffPtBin[i]->Draw("P");
	  funcJpsiEffPtBin[i]->SetLineColor(2);
	  funcJpsiEffPtBin[i]->Draw("sames");
	  hSysMtdRespEff->SetBinContent(i+1,1);
	  hSysMtdRespEff->SetBinError(i+1,sigma/mean);
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",i*0.5,i*0.5+0.5),0.07);
	  t1->Draw();
	}
      if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_FitProfile.pdf",run_type));
      hSysMtdRespEff->SetMarkerStyle(21);
      c = draw1D(hSysMtdRespEff);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/SysMtdRespEff1.pdf",run_type));

      if(saveHisto)
	{
	  TFile *fout = 0x0;
	  if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","update");
	  if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root","update");
	  if(year==2015) fout = TFile::Open("Rootfiles/Run15.pp200.MtdResponseEff.root","update");
	  hSysMtdRespEff->Write("",TObject::kOverwrite);
	}
    }
  else if(mode==3)
    {
      // if(year==2013)
      // 	{
      // 	  printf("[e] this mode does not apply to 2013!\n");
      // 	  return;
      // 	}
      TH1F *hSysMtdRespEff = new TH1F(Form("hSysMtdRespEff_2"),"Systematic uncertainty of MTD response efficiency;J/#Psi p_{T} (GeV/c)",nPtBins,minPt,maxPt);

      // average the difference for bottom backleg
      TH1F *hRatio = new TH1F("hRatio","",1000,0,20);
      TH1F *hRatioMod[150];
      int nBins = hRatio->GetNbinsX();
      int nMods = 0;
      for(int i=0; i<150; i++)
	{
	  int bl = i/5 + 1;
	  int mod = i%5 +1;
	  hRatioMod[i] = new TH1F(Form("hRatio_%d",i),"",1000,0,20);
	  if(bl<=9 || bl>=23) continue;
	  if(!isInMtd(bl,mod,0)) continue;
	  nMods ++;
	  TF1 *func1 = (TF1*)fRespEff->Get(Form("fPtMtdEffBkl%d_Mod%d",i/5,i%5));
	  TF1 *func2 = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",i/5,i%5));
	  for(int bin=1; bin<=hRatioMod[i]->GetNbinsX(); bin++)
	    {
	      double x = hRatioMod[i]->GetBinCenter(bin);
	      double ratio = 1;
	      if(x>1) ratio = func1->Eval(x)/func2->Eval(x);
	      hRatioMod[i]->SetBinContent(bin,ratio);
	      hRatioMod[i]->SetBinError(bin,1e-10);
	    }
	}
      printf("[i] %d modules are used for average\n",nMods);
      for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
	{
	  double value = 0;
	  for(int i=0; i<150; i++)
	    {
	      if(hRatioMod[i]->GetEntries()<1) continue;
	      value += hRatioMod[i]->GetBinContent(bin);
	    }
	  hRatio->SetBinContent(bin,value/nMods);
	  hRatio->SetBinError(bin,1e-10);
	}
      hRatio->SetMarkerStyle(25);
      hRatio->SetLineStyle(2);
      hRatio->GetXaxis()->SetRangeUser(0,10);
      hRatio->GetYaxis()->SetRangeUser(0.6,1.3);
      TCanvas *c1 = draw1D(hRatio,"Average difference to Template for bottom backlegs;p_{T} (GeV/c);Cosmic/Template");
      if(savePlot) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/RatioToTemplate_BottomBL.pdf",run_type));

      const int nIter = 2;
      TH1F *hSysInputJpsi[nIter];
      TH1F *hSysMtdJpsiPt[nIter];
      TH1F *hSysAccJpsiPt[nIter];
      TH1F *hSysJpsiPtEff[nIter];
      TH1F *hRespEffTmp[150];
      for(int ih = 0; ih<nIter; ih++)
	{
	  hSysInputJpsi[ih] = new TH1F(Form("hSysInputJpsi_%d",ih),"p_{T} distribution of input J/psi;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  hSysMtdJpsiPt[ih] = new TH1F(Form("hSysMtdJpsiPt_%d",ih),"p_{T} distribution of J/psi in MTD;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  hSysAccJpsiPt[ih] = new TH1F(Form("hSysAccJpsiPt_%d",ih),"p_{T} distribution of triggered J/psi;p_{T} (GeV/c)",nPtBins,minPt,maxPt);
	  
	  for(int i=0; i<150; i++)
	    {
	      hRespEffTmp[i] = new TH1F(Form("hRespEffTmp_%d_%d",i,ih),"",1000,0,20);
	      TF1 *func1 = (TF1*)fRespEff->Get(Form("fPtMtdEffBkl%d_Mod%d",i/5,i%5));
	      TF1 *func2 = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",i/5,i%5));
	      for(int bin=1; bin<=hRespEffTmp[i]->GetNbinsX(); bin++)
		{
		  double x = hRespEffTmp[i]->GetBinCenter(bin);
		  double ratio = 1;
		  if(ih==1) ratio = hRatio->GetBinContent(hRatio->FindFixBin(x));
		  if(year==2013)
		    {
		      hRespEffTmp[i]->SetBinContent(bin,func2->Eval(x)*ratio);
		    }
		  else
		    {
		      if(i/5>=9 && i/5<=21) hRespEffTmp[i]->SetBinContent(bin,func1->Eval(x));
		      else                  hRespEffTmp[i]->SetBinContent(bin,func2->Eval(x)*ratio);
		    }
		  hRespEffTmp[i]->SetBinError(bin,1e-10);
		}
	      //printf("BL = %d, Mod = %d, err = %4.4f, scale = %4.4f\n",i/5+1,i%5+1,error,scale);
	    }
	  toyMC(jpsiMass, 1e9, 0, 1, hMcJpsiPt, hRespEffTmp, hSysInputJpsi[ih], hSysMtdJpsiPt[ih], hSysAccJpsiPt[ih]);
	  hSysJpsiPtEff[ih] = (TH1F*)hSysAccJpsiPt[ih]->Clone(Form("hSysJpsiPtEff_%d",ih));
	  hSysJpsiPtEff[ih]->Sumw2();
	  hSysJpsiPtEff[ih]->Divide(hSysMtdJpsiPt[ih]);
	  hSysJpsiPtEff[ih]->SetMarkerStyle(20);
	  hSysJpsiPtEff[ih]->SetMarkerColor(ih+1);
	  hSysJpsiPtEff[ih]->SetLineColor(ih+1);
	  hSysJpsiPtEff[ih]->GetYaxis()->SetRangeUser(0,1);
	}
      TCanvas *c = draw1D(hSysJpsiPtEff[0],"MTD response efficiency for J/psi");
      hSysJpsiPtEff[1]->Draw("sames");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_JpsiEff_Iterations2.pdf",run_type));

      TH1F *hCheck = (TH1F*)hSysJpsiPtEff[1]->Clone("hCheck");
      hCheck->Divide(hSysJpsiPtEff[0]);
      hCheck->GetYaxis()->SetRangeUser(0.8,1.2);
      c = draw1D(hCheck);
      TF1 *funcCheck = new TF1("funcCheck","pol0",0,10);
      hCheck->Fit(funcCheck,"IR");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_JpsiEffRatio.pdf",run_type));
      
      for(int i=0; i<nPtBins; i++)
	{
	  hSysMtdRespEff->SetBinContent(i+1,1);
	  double error = fabs(1-hSysJpsiPtEff[1]->GetBinContent(i+1)/hSysJpsiPtEff[0]->GetBinContent(i+1));
	  hSysMtdRespEff->SetBinError(i+1,error);
	  printf("[i] pT = %1.1f: stat %6.2f%%, diff %6.2f%%\n",hSysJpsiPtEff[1]->GetBinCenter(i+1),
		 hSysJpsiPtEff[0]->GetBinError(i+1)/hSysJpsiPtEff[0]->GetBinContent(i+1)*100,
		 error * 100);
	}

      hSysMtdRespEff->SetMarkerStyle(21);
      c = draw1D(hSysMtdRespEff);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/SysMtdRespEff2.pdf",run_type));

      if(saveHisto)
	{
	  TFile *fout = 0x0;
	  if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","update");
	  if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root","update");
	  if(year==2015) fout = TFile::Open("Rootfiles/Run15.pp200.MtdResponseEff.root","update");
	  hSysMtdRespEff->Write("",TObject::kOverwrite);
	}
    }
  cout << "Done :) " << endl;
}


//================================================
void toyMC(const double mass, const int nExpr, const int debug, const int inputType, TH1F *hJpsiTruth,  TH1F *hRespEff[150],
	   TH1F *hInputJpsiPt, TH1F *hMtdJpsiPt, TH1F *hAccJpsiPt,
	   TH1F *hMuonMap, TH1F *hMuonMapTriggered)
{
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("+++ Event %d +++\n",i+1);
      double mc_pt;
      if(inputType==0)  mc_pt =  hJpsiTruth->GetRandom();
      if(inputType==1)  mc_pt = myRandom->Uniform(0,12);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      hInputJpsiPt->Fill(parent.Pt());
      hInputJpsiEta->Fill(parent.Eta());
      hInputJpsiPhi->Fill(parent.Phi());
      hInputJpsiY->Fill(parent.Rapidity());
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());

      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      double pt1 = daughter1.Pt();
      double eta1 = daughter1.Eta();
      double phi1 = daughter1.Phi();
      double z1 = TMath::Tan(pi/2-daughter1.Theta()) * gMtdRadius;
      int backleg1, module1, cell1;
      getMtdPos(phi1, z1, backleg1, module1, cell1);
      if(debug) printf("daugther 1: pt = %3.2f eta = %3.2f phi = %3.2f z = %3.2f\n",pt1,eta1,phi1/pi*180,z1);
      if(debug) printf("            backleg = %d, module = %d, cell = %d\n",backleg1,module1,cell1);
      bool isMtd1 = isInMtd(backleg1, module1, cell1);
      int isAcc1 = 0;
      if(isMtd1)
	{
	  int index1 = (backleg1-1)*5+module1-1;
	  int bin1 = hRespEff[index1]->FindBin(pt1);
	  double eff1 = hRespEff[index1]->GetBinContent(bin1);
	  double prob1 = myRandom->Uniform(0,1);
	  if(prob1<eff1)       isAcc1 = 1;
	  else           isAcc1 = 0;

	  //if(backleg1==4 && module1==1)
	    //printf("1: prob1 is less than eff1: %f <? %f = %d, pt1 = %f, bl = %d, mod = %d\n",prob1,eff1,isAcc1, pt1,backleg1,module1);
	}
      
      TLorentzVector daughter2 = parent - daughter1;
      double pt2 = daughter2.Pt();
      double eta2 = daughter2.Eta();
      double phi2 = daughter2.Phi();
      double z2 = TMath::Tan(pi/2-daughter2.Theta()) * gMtdRadius;
      int backleg2, module2, cell2;
      getMtdPos(phi2, z2, backleg2, module2, cell2);
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);
      bool isMtd2 = isInMtd(backleg2, module2, cell2);
      bool isAcc2 = false;
      if(isMtd2)
	{
	  int index2 = (backleg2-1)*5 + module2-1;
	  int bin2 = hRespEff[index2]->FindFixBin(pt2);
	  double eff2 = hRespEff[index2]->GetBinContent(bin2);
	  double prob2 = myRandom->Uniform(0,1);
	  if(prob2<eff2) isAcc2 = true;
	  else           isAcc2 = false;
	}

      double dphi = rotatePhi(daughter2.Phi()-daughter1.Phi());
      if(dphi>pi) dphi = 2*pi - dphi;
      hJpsiPtVsOpenAngle->Fill(parent.Pt(), dphi);

      double leadpt = pt1 > pt2 ? pt1 : pt2;

      if(leadpt<pt1_cut || pt1<pt2_cut || pt2<pt2_cut) continue;

      if(hAccMuonMap)
	{
	  if(isMtd1) hAccMuonMap->Fill((backleg1-1)*5+module1);
	  if(isMtd2) hAccMuonMap->Fill((backleg2-1)*5+module2);
	}

      if(hMuonMap)
	{
	  if(isAcc1) hMuonMap->Fill((backleg1-1)*5+module1);
	  if(isAcc2) hMuonMap->Fill((backleg2-1)*5+module2);
	}

      if(!isMtd1 || !isMtd2) continue;
      hMtdJpsiPt->Fill(parent.Pt());
      
      if(!isAcc1 || !isAcc2) continue;
      hAccJpsiPt->Fill(parent.Pt());

      if(hMuonMapTriggered)
	{
	  hMuonMapTriggered->Fill((backleg1-1)*5+module1);
	  hMuonMapTriggered->Fill((backleg2-1)*5+module2);
	}
    }
}


//================================================
void anaSys(const int savePlot)
{
  const int nBL = 13;
  TCanvas *cBL[nBL];
  TH1F *hCosmic[nBL][5];
  TF1 *hFitCosmic[nBL][5];
  TF1 *hFitTemplate[nBL][5];
  TH1F *hRatio[nBL][5];

  TFile *fRespEff = TFile::Open(Form("Rootfiles/Run%dResponseEffViaPtTemplate.root",year-2000),"read");
  for(int i=0; i<nBL; i++)
    {
      int nMod = 5;
      int bl = i+10;
      if(bl>11 && bl<21) 
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1100,700);
	  cBL[i]->Divide(3,2);
	  nMod = 3;
	}
      else
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1500,600);
	  cBL[i]->Divide(5,2);
	}

      
      for(int j=0; j<5; j++)
	{
	  hRatio[i][j] = new TH1F(Form("hRatio_%d_%d",bl,j+1),"",500,0,10);
	  hFitCosmic[i][j] = (TF1*)fRespEff->Get(Form("fPtMtdEffBkl%d_Mod%d",bl-1,j));
	  hFitTemplate[i][j] = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",bl-1,j));
	  hCosmic[i][j] = (TH1F*)fRespEff->Get(Form("hPtMtdEffBkl%d_Mod%d",bl-1,j));
	  int index = j;
	  if(bl>11 && bl<21)
	    {
	      if(j==0 || j==4) continue;
	      index = j-1;
	    }
	  
	  cBL[i]->cd(index+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  ScaleHistoTitle(hCosmic[i][j],0.05,1,0.04,0.05,1,0.04,62);
	  hCosmic[i][j]->SetTitle(";p_{T} (GeV/c);MTD response efficiency");
	  hCosmic[i][j]->GetXaxis()->SetRangeUser(0,4);
	  hCosmic[i][j]->Draw();
	  hFitCosmic[i][j]->Draw("sames");
	  hFitTemplate[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("Bkl=%d, Mod=%d",bl,j+1),0.06);
	  t1->Draw();
	  TLegend *leg = new TLegend(0.45,0.2,0.7,0.42);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->AddEntry(hCosmic[i][j],"Cosmic data","P");
	  leg->AddEntry(hFitCosmic[i][j],"Fit to cosmic","L");
	  leg->AddEntry(hFitTemplate[i][j],"Template","L");
	  leg->Draw();

	  for(int bin=1; bin<=hRatio[i][j]->GetNbinsX(); bin++)
	    {
	      double x = hRatio[i][j]->GetBinCenter(bin);
	      if(x<1) continue;
	      hRatio[i][j]->SetBinContent(bin,hFitCosmic[i][j]->Eval(x)/hFitTemplate[i][j]->Eval(x));
	      hRatio[i][j]->SetBinError(bin,1e-10);
	    }
	  cBL[i]->cd(index+nMod+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  hRatio[i][j]->SetMarkerStyle(25);
	  hRatio[i][j]->SetLineStyle(2);
	  hRatio[i][j]->SetTitle(";p_{T} (GeV/c);Template/Cosmic");
	  ScaleHistoTitle(hRatio[i][j],0.05,1,0.04,0.05,1,0.04,62);
	  hRatio[i][j]->GetXaxis()->SetRangeUser(0,4);
	  hRatio[i][j]->GetYaxis()->SetRangeUser(0,2);
	  hRatio[i][j]->Draw("P");
	  TLine *line = new TLine(0,1,4,1);
	  line->SetLineColor(2);
	  line->SetLineStyle(2);
	  line->Draw();
	}
      if(savePlot) cBL[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/BL%d_TemplateVsCosmic.pdf",run_type,bl));
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
	 (backleg>=12 && backleg<=20 && (module==1 || module==5)) ||
	 (backleg==15 && module==4))
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

//================================================
void singleMuon(const int savePlot)
{
  TH1F *hMuonPtTruth;
  TH1F *hMuonPtMtd;
  TH1F *hMuonPtResp;
  TH1F *hMuonMapResp;
  TH1F *hMuonEffPerMod;
  TH1F *hMuonEffAvgMod;

  hMuonPtTruth = new TH1F("hMuonPtTruth","p_{T} distribution of input muons;p_{T,#mu} (GeV/c)",100,0,10);
  hMuonPtMtd = new TH1F("hMuonPtMtd","p_{T} distribution of muons in MTD acceptance;p_{T,#mu} (GeV/c)",100,0,10);
  hMuonPtResp = new TH1F("hMuonPtResp","p_{T} distribution of muons with MTD response;p_{T,#mu} (GeV/c)",100,0,10);
  hMuonMapResp = new TH1F("hMuonMapResp","Muon map with MTD response;(backleg-1)*5+module",150,0.5,150.5);

  TH1F *hRespEff[150];
  TFile *fRespEff = 0x0;
  if(year==2013) fRespEff = TFile::Open("Rootfiles/Run13ResponseEfficiency2016Mar24.root","read");
  for(int i=0; i<150; i++)
    {
      int backleg = i/5 + 1;
      int module = i%5 + 1;
      TF1 *func = (TF1*)fRespEff->Get(Form("ResPonseEffvsPt_Bkl%d_Mod%d",backleg,module));
      hRespEff[i] = new TH1F(Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),Form("MtdRespEffvsPt_Bkl%d_Mod%d",backleg,module),1000,0,20);
      for(int bin=1; bin<=hRespEff[i]->GetNbinsX(); bin++)
	{
	  hRespEff[i]->SetBinContent(bin,func->Eval(hRespEff[i]->GetBinCenter(bin)));
	  hRespEff[i]->SetBinError(bin,1e-10);
	}
    }

  const int nExpr = 1e6;
  for(int i=0; i<nExpr; i++)
    {
      double mc_pt  =  myRandom->Uniform(0,15);
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_eta = myRandom->Uniform(-0.8, 0.8);
      TLorentzVector muon;
      muon.SetPtEtaPhiM (mc_pt, mc_eta, mc_phi, muMass);
      hMuonPtTruth->Fill(mc_pt);

      double z = TMath::Tan(pi/2-muon.Theta()) * gMtdRadius;
      int backleg, module, cell;
      getMtdPos(mc_phi, z, backleg, module, cell);
      bool isMtd = isInMtd(backleg, module, cell);
      if(!isMtd) continue;
      hMuonPtMtd->Fill(mc_pt);
      hMuonMapResp->Fill((backleg-1)*5+module);

      int isAcc = 0;
      int index = (backleg-1)*5+module-1;
      int bin = hRespEff[index]->FindBin(mc_pt);
      double eff = hRespEff[index]->GetBinContent(bin);
      double prob = myRandom->Uniform(0,1);
      if(prob<eff) isAcc = 1;
      else         isAcc = 0;
      if(!isAcc) continue;
      hMuonPtResp->Fill(mc_pt);

    }

  TCanvas *c = draw1D(hMuonPtTruth);
  hMuonEffPerMod = (TH1F*)hMuonPtResp->Clone("hMuonEffPerMod");
  hMuonEffPerMod->Sumw2();
  hMuonEffPerMod->Divide(hMuonPtMtd);
  hMuonEffPerMod->SetMarkerStyle(21);
  hMuonEffPerMod->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hMuonEffPerMod,"MTD response efficiency for single muon");

  // calculate average efficiency
  TH1F *hWeight = (TH1F*)hMuonMapResp->Clone("hWeight");
  hWeight->Scale(1./hWeight->Integral());
  hMuonEffAvgMod = (TH1F*)hRespEff[0]->Clone("hMuonEffAvgMod");
  hMuonEffAvgMod->Reset();
  for(int bin=1; bin<=hMuonEffAvgMod->GetNbinsX(); bin++)
    {
      double value = 0;
      for(int i=0; i<150; i++)
	{
	  double weight = hWeight->GetBinContent(i+1);
	  double eff = hRespEff[i]->GetBinContent(bin);
	  value += weight * eff;
	}
      hMuonEffAvgMod->SetBinContent(bin,value);
      hMuonEffAvgMod->SetBinError(bin,1e-10);
    }
  hMuonEffAvgMod->SetMarkerStyle(25);
  hMuonEffAvgMod->SetMarkerColor(2);
  hMuonEffAvgMod->SetLineColor(2);
  hMuonEffAvgMod->Draw("sames");

  TLegend *leg = new TLegend(0.4,0.3,0.6,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hMuonEffPerMod,"Module-wise efficiency","PL");
  leg->AddEntry(hMuonEffAvgMod,"Module-average efficiency","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_Check_SingleMuonEff.pdf",run_type));
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
