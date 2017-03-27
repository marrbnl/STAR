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

const int year = 2015;
const char *run_type = "Run15_pAu200";
const double pi = 3.1415926;
const double muMass = 0.1057;
const double smear[2] = {0.0133, 0.0161};
const int    gMtdNChannels             = 24;   // Total number of MTD channels per module. One cell has two channels.
const double gMtdCellLength            = 87.0; // Length of a MTD cell (cm)
const double gMtdCellWidth             = 3.8;  // Width of a MTD cell (cm)
const double gMtdCellGap               = 0.6;  // Gap between MTD cells (cm)
const double gMtdBacklegPhiWidth       = 8.*pi/180.;   // Width of backleg in phi (rad)
const double gMtdBacklegPhiGap         = 4.*pi/180.;   // Gap between backleg in phi (rad)
const double gMtdFirstBacklegPhiCenter = 90.*pi/180.;  // Center of backleg 1 at phi = 90 degree (rad)
const double gMtdRadius                = 403; // Minimum radius of MTD system extracted from geometry file (cm)

TRandom3 *myRandom;
TH2F *hTrkResVsPt[2];
TH1F *hTrkResBin[2][400];

void anaEff(const int savePlot);
double toyMC(const int index, const double mass, const TH1F *hInput, const int nExpr, const TF1 *hEmbEff, const TH1F *hPidEff, TH1F *hRecMuonPt, const int debug=0);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
void getMtdPos(const double phi, const double z, int &backleg, int &module, int &cell);
double rotatePhi(double phi);
bool isInMtd(const int backleg, const int module, int cell);

//================================================
void ana_Psi2S()
{
  gStyle->SetOptStat(0);

  //makePlots();
  anaEff(0);
}

//================================================
void anaEff(const int savePlot)
{
  const char *data_name[2] = {"pp","pAu"};
  const char *part_name[2] = {"Jpsi","Psi2S"};
  const double part_mass[2] = {3.097, 3.686};

  TFile *fin = TFile::Open(Form("Rootfiles/%s.Psi2S.root",run_type),"read");
  TH1F *hJpsi = (TH1F*)fin->Get("hJpsi_default");
  TH1F *hPsi2SInput[3];
  hPsi2SInput[0] = (TH1F*)fin->Get("hPsi2S_default");
  hPsi2SInput[1] = (TH1F*)fin->Get("hPsi2S_low");
  hPsi2SInput[2] = (TH1F*)fin->Get("hPsi2S_high");
  TF1 *funcEmbEff[2];
  TH1F *hPidEffInput[2][2];
  for(int i=0; i<2; i++)
    {
      funcEmbEff[i] = (TF1*)fin->Get(Form("Run15_%s200_EmbEffFit",data_name[i]));
      hPidEffInput[i][0] = (TH1F*)fin->Get(Form("Run15_%s200_PidEff_fit",data_name[i]));
      hPidEffInput[i][1] = (TH1F*)fin->Get(Form("Run15_%s200_PidEff_data",data_name[i]));

      TFile* fRes = TFile::Open(Form("Rootfiles/Run15_%s200.PtResolution.root",data_name[i]),"read");
      hTrkResVsPt[i] = (TH2F*)fRes->Get("PrimTrkRes_vs_TruePt_cent0080");
      hTrkResVsPt[i]->SetName(Form("Run15_%s200_%s",data_name[i],hTrkResVsPt[i]->GetName()));
      int nHistos = hTrkResVsPt[i]->GetNbinsX();
      for(int bin=0; bin<nHistos; bin++)
	{
	  hTrkResBin[i][bin] = (TH1F*)hTrkResVsPt[i]->ProjectionY(Form("Run15_%s200_hTrkRes_Bin%d",data_name[i],bin+1),bin+1,bin+1);
	}
    }


  const int nExpr = 1e7;
  myRandom = new TRandom3();
  TDatime *clock = new TDatime();
  myRandom->SetSeed(clock->GetTime());

  double eff[2][2];
  TH1F *hPidEff[2] = {0x0, 0x0};
  TH1F *hInput[2] = {0x0, 0x0};
  const char *sys_name[4] = {"Default", "Psi2S_low","Psi2S_high","PidEffDataPoint"};
  TH1F *hRecMuonPt[4][2][2];
  for(int k=0; k<4; k++)
    {
      // default setup
      hInput[0] = hJpsi; 
      hInput[1] = hPsi2SInput[0];
      hPidEff[0] = hPidEffInput[0][0]; 
      hPidEff[1] = hPidEffInput[1][0]; 

      printf("\n[i] %s\n",sys_name[k]);
      if(k==1 || k==2) hInput[1] = hPsi2SInput[k];
      if(k==3) 
	{
	  hPidEff[0] = hPidEffInput[0][1];
	  hPidEff[1] = hPidEffInput[1][1]; 
	}
      for(int i=0; i<2; i++) // 0 - pp; 1 - pAu
	{
	  printf("[i] %s\n",data_name[i]);
	  for(int j=0; j<2; j++) // 0 - Jpsi; 1 - Psi(2S)
	    {
	      hRecMuonPt[k][i][j] = new TH1F(Form("hMuonPt_%s_%s_Sys%d",data_name[i],part_name[j],k),"",50,0,10);
	      eff[i][j] = toyMC(i, part_mass[j], hInput[j], nExpr, funcEmbEff[i], hPidEff[i], hRecMuonPt[k][i][j], 0);
	      printf("[i]%s: %4.2f%%\n", part_name[j],eff[i][j]*100);
	    }
	  printf("[i] Rel. eff. = %4.2f\n",eff[i][1]/eff[i][0]);
	}
    }

  TCanvas *c = new TCanvas("RecMuonPt","RecMuonPt",800,600);
  gPad->SetLogy();
  for(int j=0; j<2; j++) 
    {
      hRecMuonPt[0][0][j]->SetTitle("p_{T} distribution of decayed muons;p_{T}");
      hRecMuonPt[0][0][j]->Sumw2();
      hRecMuonPt[0][0][j]->Scale(1./hRecMuonPt[0][0][j]->Integral());
      hRecMuonPt[0][0][j]->SetMarkerStyle(21);
      hRecMuonPt[0][0][j]->SetMarkerColor(j+1);
      hRecMuonPt[0][0][j]->SetLineColor(j+1);
    }
  hRecMuonPt[0][0][0]->Draw();
  hRecMuonPt[0][0][1]->Draw("sames");
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hRecMuonPt[0][0][0],"J/psi","PL");
  leg->AddEntry(hRecMuonPt[0][0][1],"Psi(2S)","PL");
  leg->Draw();
  if(1)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Psi2S/RecMuonPt.pdf",run_type));
}

//================================================
double toyMC(const int index, const double mass, const TH1F *hInput, const int nExpr, const TF1 *hEmbEff, const TH1F *hPidEff, TH1F *hRecMuonPt, const int debug)
{
  const int nHisto = hTrkResVsPt[index]->GetNbinsX();
  double pt1_cut = 1.5;
  double pt2_cut = 1.3;
  int nInput = 0, nOutput = 0;
  for(int i=0; i<nExpr; i++)
    {
      if(debug) printf("+++ Event %d +++\n",i+1);
      double mc_pt  = hInput->GetRandom();
      double mc_phi = myRandom->Uniform(-1*pi, pi);
      double mc_y   = myRandom->Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+mass*mass) * TMath::SinH(mc_y);
      nInput ++;

      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,mass);
      if(debug) printf("parent:     pt = %3.2f eta = %3.2f phi = %3.2f\n",parent.Pt(),parent.Eta(),parent.Phi());

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
      if(debug) printf("daugther 2: pt = %3.2f eta = %3.2f phi = %3.2f\n",pt2,eta2,phi2);

      // acceptance cut
      if(fabs(eta1)>0.5 || fabs(eta2)>0.5) continue;
      if(pt1<0.5 || pt2<0.5) continue;

      // momentum resolution
      int mom_index1 = hTrkResVsPt[index]->GetXaxis()->FindBin(pt1)-1;
      if(mom_index1>=nHisto) mom_index1=nHisto-1;
      int mom_index2 = hTrkResVsPt[index]->GetXaxis()->FindBin(pt2)-1;
      if(mom_index2>=nHisto) mom_index2=nHisto-1;
      double dpt1 = hTrkResBin[index][mom_index1]->GetRandom();
      double dpt2 = hTrkResBin[index][mom_index2]->GetRandom();
      double emb_pt1 = (1-dpt1) * pt1;
      double emb_pt2 = (1-dpt2) * pt2;
      double rc_pt1 = emb_pt1 * myRandom->Gaus(1,sqrt(emb_pt1)*smear[index]);
      double rc_pt2 = emb_pt2 * myRandom->Gaus(1,sqrt(emb_pt2)*smear[index]);
      if(debug) cout << pt1 << "  " << emb_pt1 << "  " << rc_pt1 << endl;

      hRecMuonPt->Fill(rc_pt1);
      hRecMuonPt->Fill(rc_pt2);
      if(rc_pt1<pt1_cut && rc_pt2<pt1_cut) continue;
      if(rc_pt1<pt2_cut || rc_pt2<pt2_cut) continue;


      double tmp_pt1 = rc_pt1 > 10 ? 10 : rc_pt1;
      double tmp_pt2 = rc_pt2 > 10 ? 10 : rc_pt2;
      
      // TPC+MTD efficiency
      double emb_eff1 = hEmbEff->Eval(tmp_pt1);
      double emb_eff2 = hEmbEff->Eval(tmp_pt2);
      if(myRandom->Uniform(0., 1.)>emb_eff1 || myRandom->Uniform(0,1)>emb_eff2) continue;

      // PID efficiency
      double pid_eff1 = hPidEff->GetBinContent(hPidEff->FindFixBin(tmp_pt1));
      double pid_eff2 = hPidEff->GetBinContent(hPidEff->FindFixBin(tmp_pt2));
      if(myRandom->Uniform(0., 1.)>pid_eff1 || myRandom->Uniform(0,1)>pid_eff2) continue;

      // MTD acceptance
      int backleg1, module1, cell1;
      getMtdPos(phi1, z1, backleg1, module1, cell1);
      bool isMtd1 = isInMtd(backleg1, module1, cell1);
      int backleg2, module2, cell2;
      getMtdPos(phi2, z2, backleg2, module2, cell2);
      bool isMtd2 = isInMtd(backleg2, module2, cell2);
      if(!isMtd1 || !isMtd2) continue;

      nOutput ++;
    }
  double eff = nOutput*1./nInput; 
  double err = eff * sqrt(nOutput)/nOutput;
  printf("++++++++++++++++++++++++++++++++++++++\n");
  printf("Input %d, output %d, eff = %4.2f+/-%4.3f%%\n",nInput,nOutput,eff*100,err*100);
  return eff;
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
double rotatePhi(double phi) 
{
  double outPhi = phi;
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
}

/*
//================================================
void makePlots(const int savePlot = 1, const int saveHisto = 1)
{
  // spectrum shape
  TFile *fshape = TFile::Open("Rootfiles/JpsiSpectraShapepp200.root","read");
  TF1 *funcJpsi = (TF1*)fshape->Get("TsallisPowerLawFitJpsipp200");
  TH1F *hJpsi = new TH1F("hJpsi_default", "hJpsi_default", 100, 0, 10);
  for(int bin=1; bin<=hJpsi->GetNbinsX(); bin++)
    {
      double pt = hJpsi->GetBinCenter(bin);
      hJpsi->SetBinContent(bin, funcJpsi->Eval(pt)*pt);
      hJpsi->SetBinError(bin, 1e-10);
    }
  hJpsi->SetMarkerStyle(21);
  TCanvas *c = draw1D(hJpsi,"p_{T} distriubtion of J/psi in p+p collisions at 200 GeV;p_{T} (GeV/c);dN/dp_{T}");
  gPad->SetLogy();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Psi2S/Jpsi_shape.pdf",run_type));
  
  TH1F *hPsi2S[3];
  TF1 *funcRatio[3];
  const double p0[3] = {0.01341, 0.01341, 0.01341};
  const double p1[3] = {0.002513, 0.001013, 0.004013};
  const char* name[3] = {"default","low","high"};
  for(int i=0; i<3; i++)
    {
      funcRatio[i] = new TF1(Form("Psi2SToJpsiRatio_%d",i),"pol1",0,10);
      funcRatio[i]->SetParameters(p0[i], p1[i]);

      hPsi2S[i] = new TH1F(Form("hPsi2S_%s",name[i]), Form("hPsi2S_%s",name[i]), 100, 0, 10);
      for(int bin=1; bin<=hPsi2S[i]->GetNbinsX(); bin++)
	{
	  double pt = hPsi2S[i]->GetBinCenter(bin);
	  hPsi2S[i]->SetBinContent(bin, funcJpsi->Eval(pt)*funcRatio[i]->Eval(pt)*pt);
	  hPsi2S[i]->SetBinError(bin, 1e-10);
	}
    }
  TList *list = new TList;
  for(int i=0; i<3; i++)
    {
      list->Add(hPsi2S[i]);
    }
  TString legName[3] = {"Default shape","Lower limit","Upper limit"};
  c = drawHistos(list,"Psi2S_Shape","p_{T} distribution of #Psi(2S) in p+p collisions at 200 GeV;p_{T} (GeV/c);dN/dp_{T}",kFALSE,0,30,false,0.5,1.1,true,kTRUE,legName,kTRUE,"",0.3,0.5,0.16,0.36,true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Psi2S/Psi2S_shape.pdf",run_type));


  // TPC+MTD efficiency
  TFile *fEmbEff[2];
  TH1F *hInputEmbEff[2];
  TH1F *hEmbEff[2];
  TF1 *funcEmbEff[2];
  const char *data_name[2] = {"pp","pAu"};
  for(int i=0; i<2; i++)
    {
      fEmbEff[i] = TFile::Open(Form("Rootfiles/Run15%s200SingleMuonEffMtdMth.root",data_name[i]));
      hInputEmbEff[i] = (TH1F*)fEmbEff[i]->Get("mhTpcEff");
      hInputEmbEff[i]->SetName(Form("hInputEmbEff_%d",i));
      hEmbEff[i] = new TH1F(Form("Run15_%s200_EmbEff",data_name[i]),"",120,0,12);
      for(int bin=1; bin<=hEmbEff[i]->GetNbinsX(); bin++)
	{
	  hEmbEff[i]->SetBinContent(bin, hInputEmbEff[i]->GetBinContent(bin));
	  hEmbEff[i]->SetBinError(bin, hInputEmbEff[i]->GetBinError(bin));
	}
      funcEmbEff[i] = new TF1(Form("Run15_%s200_EmbEffFit",data_name[i]),"[0]-exp([1]+[2]*x**0.5+[3]*x+[4]*x**2)",1.2,12);
      funcEmbEff[i]->SetParameters(0.35, 92, -180, 97, -9);
      hEmbEff[i]->Fit(funcEmbEff[i],"0R");
      hEmbEff[i]->SetMarkerStyle(21);
      hEmbEff[i]->GetXaxis()->SetRangeUser(0,12);
      hEmbEff[i]->GetYaxis()->SetRangeUser(0,0.5);
      c = draw1D(hEmbEff[i],Form("Run15 %s200: TPC+MTD efficiency for single muons;p_{T} (GeV/c);Eff",data_name[i]));
      funcEmbEff[i]->SetLineColor(2);
      funcEmbEff[i]->Draw("sames");
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Psi2S/EmbedEff_%s200.pdf",run_type,data_name[i]));
    }

  // PID efficiency
  TFile *fPidEff[2];
  TGraphAsymmErrors *gInputPidEff[2];
  TH1F *hPidEff[2][2];
  TF1 *funcPidEff[2];
  const int nbins = 7;
  const double xbins[nbins+1] = {1.3,2,2.5,3,3.5,4,5,10};
  double x,y;
  for(int i=0; i<2; i++)
    {
      fPidEff[i] = TFile::Open(Form("Rootfiles/Run15%s200_PIDEff_SingleMuon.root",data_name[i]));
      gInputPidEff[i] = (TGraphAsymmErrors*)fPidEff[i]->Get("Eff_LR");
      gInputPidEff[i]->SetName(Form("hInputPidEff_%d",i));

      // data points
      hPidEff[i][0] = new TH1F(Form("Run15_%s200_PidEff_data",data_name[i]),"",nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  gInputPidEff[i]->GetPoint(bin-1, x, y);
	  hPidEff[i][0]->SetBinContent(bin, y);
	  hPidEff[i][0]->SetBinError(bin, gInputPidEff[i]->GetErrorYlow(bin-1));
	}

      // fitting
      funcPidEff[i] = new TF1(Form("Run15_%s200_PidEffFit",data_name[i]),"[0]-exp([1]+[2]*x)",1.3,10);
      funcPidEff[i]->SetParLimits(0,0,1);
      funcPidEff[i]->SetParameters(0.9, 3.5,-2.6);
      hPidEff[i][0]->Fit(funcPidEff[i],"0R");
      gInputPidEff[i]->SetMarkerStyle(21);
      gInputPidEff[i]->GetXaxis()->SetRangeUser(0,10);
      gInputPidEff[i]->GetXaxis()->SetLabelSize(0.04);
      gInputPidEff[i]->GetXaxis()->SetTitleSize(0.045);
      gInputPidEff[i]->GetXaxis()->SetTitleOffset(1);
      gInputPidEff[i]->GetYaxis()->SetRangeUser(0,1.2);
      gInputPidEff[i]->GetYaxis()->SetLabelSize(0.04);
      gInputPidEff[i]->GetYaxis()->SetTitleSize(0.045);
      gInputPidEff[i]->GetYaxis()->SetTitleOffset(1);
      c = drawGraph(gInputPidEff[i],Form("Run15 %s200: PID efficiency of likelihood ratio method;p_{T} (GeV/c);Eff",data_name[i]));
      funcPidEff[i]->SetLineColor(2);
      funcPidEff[i]->Draw("sames");
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Psi2S/PidEff_%s200.pdf",run_type,data_name[i]));
      funcPidEff[i]->SetNpx(100);
      hPidEff[i][1] = (TH1F*)funcPidEff[i]->GetHistogram();
      hPidEff[i][1]->SetName(Form("Run15_%s200_PidEff_fit",data_name[i]));
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Psi2S.root",run_type),"recreate");
      hJpsi->Write();
      for(int i=0; i<3; i++)
	{
	  hPsi2S[i]->Write();
	}
      for(int i=0; i<2; i++)
	{
	  hEmbEff[i]->Write();
	  funcEmbEff[i]->Write();
	  hPidEff[i][0]->Write();
	  funcPidEff[i]->Write();
	  hPidEff[i][1]->Write();
	}
    }
}
*/
