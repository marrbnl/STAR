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
void makeMtdAcc(const int savePlot = 0, const int saveAN = 0);

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
      TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
      hMcJpsiPt = (TH1F*)fWeight->Get(Form("TBW_JpsiYield_AuAu200_cent0060"));
    } 


  MtdAccepLoss(1,1);
  //makeMtdAcc(0, 1);
}

//================================================
void makeMtdAcc(const int savePlot, const int saveAN)
{
  TFile *fdata = TFile::Open(Form("output/%s.jpsi.root",run_type));
  THnSparseF *mhMtdRunStatus = (THnSparseF*)fdata->Get("mhMtdRunStatus_di_mu");
  int nRange = 0;
  TH2F *hMtdMap[10];
  int runRange[10];
  if(year==2014)
    {
      nRange = 7;
      runRange[0] = 15074103;
      runRange[1] = 15078034;
      runRange[2] = 15098066;
      runRange[3] = 15099003;
      runRange[4] = 15106130;
      runRange[5] = 15131036;
      runRange[6] = 15132019;
      runRange[7] = 15167014;
    }
  
  for(int i=0; i<nRange; i++)
    {
      mhMtdRunStatus->GetAxis(2)->SetRange(i+2, i+2);
      hMtdMap[i] = (TH2F*)mhMtdRunStatus->Projection(1,0);
      hMtdMap[i]->SetName(Form("hMtdMap_RunRange%d",i+1));
      hMtdMap[i]->GetYaxis()->SetRangeUser(0,59);
      TCanvas *c = draw2D(hMtdMap[i],Form("%s: run %d - %d",run_type,runRange[i]+1,runRange[i+1])); 
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_AccepCor/MtdAcceptanceLoss_RunRange%d.pdf",run_type,i));
      if(saveAN)
	c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_MtdAccLoss%d.pdf",i+2));
    }

}

//================================================
void MtdAccepLoss(const int savePlot, const int saveHisto)
{
  const int nHistos = 9;
  TH1F *hInJpsiPt;
  TH1F *hOutJpsiPt[nHistos];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  xbins[0] = 0.15;

  hInJpsiPt = new TH1F(Form("hInJpsiPt"),"p_{T} distribution of input J/#Psi;p_{T} (GeV/c)",nbins,xbins);
  hInJpsiPt->Sumw2();
  for(int i=0; i<nHistos; i++)
    {
      hOutJpsiPt[i] = new TH1F(Form("hOutJpsiPt_%d",i),"p_{T} distribution of reconstructed J/#Psi;p_{T} (GeV/c)",nbins,xbins);
      hOutJpsiPt[i]->Sumw2();
    }

  TH1F *hRespEff[150];
  TFile *fRespEff = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");
  TF1 *func  = 0x0;
  for(int i=0; i<150; i++)
    {
      int bl = i/5 + 1;
      int mod = i%5 + 1;
      if(bl>9 && bl<23)
	{
	  func = (TF1*)fRespEff->Get(Form("Cosmic_FitRespEff_BL%d_Mod%d",bl,mod));
	}
      else
	{
	  func = (TF1*)fRespEff->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",bl,mod));
	}

      hRespEff[i] = new TH1F(Form("MtdRespEffvsPt_Bkl%d_Mod%d",bl,mod),Form("MtdRespEffvsPt_Bkl%d_Mod%d",bl,mod),1500,0,30);
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
      double mc_pt  =  myRandom->Uniform(0.15,20);
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

      if( (backleg1==24 && ((module1<=3&&cell1<=1) || (module1>3&&cell1>=10)) ) || 
	  (backleg2==24 && ((module2<=3&&cell2<=1) || (module2>3&&cell1>=10)) ) ) continue;

      if( (backleg1==8 && ((module1<=3&&cell1>=9) || (module1>3&&cell1<=2)) ) || 
	  (backleg2==8 && ((module2<=3&&cell2>=9) || (module2>3&&cell1<=2)) ) ) continue;

      if ( (backleg1==15&&module1==4) || (backleg2==15&&module2==4) ) continue;

      hOutJpsiPt[0]->Fill(mc_pt,weight);

      if( !(backleg1==8 || backleg1==19 || (backleg1==17&&module1==3&&cell1==7)) &&
	  !(backleg2==8 || backleg2==19 || (backleg2==17&&module2==3&&cell2==7)))
	{
	  hOutJpsiPt[1]->Fill(mc_pt,weight);
	  hOutJpsiPt[3]->Fill(mc_pt,weight);
	}

      if( !(backleg1==8 || backleg1==19 || backleg1==17) &&
	  !(backleg2==8 || backleg2==19 || backleg2==17))
	{
	  hOutJpsiPt[2]->Fill(mc_pt,weight);
	}


      if( !(backleg1==8 || (backleg1>=11&&backleg1<=14) || (backleg1==17&&module1==3&&cell1==7)) &&
	  !(backleg2==8 || (backleg2>=11&&backleg2<=14) || (backleg2==17&&module2==3&&cell2==7)))
	{
	  hOutJpsiPt[4]->Fill(mc_pt,weight);
	}


      if( !(backleg1==8 || (backleg1==17&&module1==3&&cell1==7) ) &&
	  !(backleg2==8 || (backleg2==17&&module2==3&&cell2==7)) )
	{
	  hOutJpsiPt[5]->Fill(mc_pt,weight);
	}

      if( !(backleg1==17&&module1==3&&cell1==7) &&
	  !(backleg2==17&&module2==3&&cell2==7))
	{
	  hOutJpsiPt[6]->Fill(mc_pt,weight);
	}

      if( !(backleg1>=15&&backleg1<=17) &&
	  !(backleg2>=15&&backleg2<=17))
	{
	  hOutJpsiPt[7]->Fill(mc_pt,weight);
	}
      hOutJpsiPt[8]->Fill(mc_pt,weight);
    }


  int runRange[nHistos] = {15074104, 15077035, 15078021, 15098066, 15099002, 15106130, 15131038, 15132019, 15167014};
  TH1F *hAcc[nHistos-1];
  for(int i=0; i<nHistos-1; i++)
    {
      hAcc[i] = (TH1F*)hOutJpsiPt[i+1]->Clone(Form("hAccepLoss_RunRange%d",i));
      hAcc[i]->Divide(hOutJpsiPt[0]);
      hAcc[i]->SetMarkerStyle(21);
      hAcc[i]->GetYaxis()->SetRangeUser(0.85,1.05);
      TCanvas *c = draw1D(hAcc[i],Form("%s: MTD acceptance loss for J/#Psi (%d-%d)",run_type,runRange[i]+1,runRange[i+1]));
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
