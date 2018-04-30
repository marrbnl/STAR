#include "TStyle.h"
#include "TSystem.h"
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
#include "TEfficiency.h"
#include "TFitResult.h"

#include <cmath>
using namespace std;
const int year = YEAR;

TFile *f;
TF1 *fResDyVsPt;
TF1 *fResDzVsPt;

//================================================
void ana_JpsiMuon()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  fResDyVsPt = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)");
  fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)");
  fResDyVsPt->SetParameters(-17.6867, 18.4528, 0.637142);
  fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.44421);

  //makeHisto();
  //makeHistoTn();
  muonPID();
  //pidSys();
  //TagAndProbe();
  //Run16nSigmaPiEff();
}

//================================================
void pidSys(const int savePlot = 1, const int saveHisto = 1)
{
  char* run_type = "Run15_pAu200";
  
  TFile *fdtof = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type), "read");
  TH1F *hdtofSys = (TH1F*)fdtof->Get(Form("%s_JpsiEffVsPt_Sys_Dtof1.00Eff",run_type));
  double sysOther = 0;
  if(run_type=="Run15_pp200")  sysOther = sqrt(0.005*0.005*3)*2;
  if(run_type=="Run15_pAu200") sysOther = sqrt(0.015*0.015+0.01*0.01*2)*2;

  for(int bin=1; bin<=hdtofSys->GetNbinsX(); bin++)
    {
      hdtofSys->SetBinError(bin, sqrt(pow(sysOther,2)+pow(hdtofSys->GetBinError(bin), 2)) );
    }
  ScaleHistoTitle(hdtofSys,0.045,1,0.035,0.045,1,0.035,42);
  hdtofSys->GetYaxis()->SetRangeUser(0.9, 1.1);
  c = draw1D(hdtofSys, Form("%s: muon PID uncertainties for J/psi",run_type));
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Final_PidEffSys.pdf",run_type));
	  
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "update");
      hdtofSys->Write(Form("%s_JpsiEffVsPt_PidSys",run_type),TObject::kOverwrite);
    }
}

//================================================
void muonPID(const int savePlot = 1, const int saveHisto = 1)
{
  char* run_type = "Run15_pp200";
  const int iSys = 1;

  const int nData = 2;
  const char* dataType[nData] = {"Data","Embed"};

  const int nVar = 4;
  const char* varName[nVar] = {"dEdx", "dY", "dZ", "dTof"};
  const char* varTitle[nVar] = {"n#sigma_{#pi}", "#Deltay", "#Deltaz", "#Deltatof"};
  const char* varUnit[nVar]  = {"", "[cm]", "[cm]", "[ns]"};

  const int nbins = 16;
  const double xbins[nbins+1] = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
  
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};

  const char* pairType[2] = {"UL", "LS"};
  const char* pidType[2] = {"noPID", "wPID"};

  TFile *fin = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "read");

  //==============================================
  // compare invariant mass
  //==============================================
  TH1F *hPairMass[nVar][nbins][2][2];
  TH1F *hPairMassCut[nVar][nbins][2][2];
  TH1F *hScale[nVar][2];
  const int pidIndex = 1;
  for(int i=0; i<nVar; i++)
    {
      for(int k=0; k<2; k++)
	{
	  hScale[i][k] = (TH1F*)fin->Get(Form("%s_%s_scaleFactor_%s",dataType[0],varName[i],pidType[k]));
	}

      TCanvas *c = new TCanvas(Form("InvMass_%s",varName[i]), Form("InvMass_%s",varName[i]), 1100, 700);
      c->Divide(2,2);
      for(int j=0; j<nbins; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      for(int t=0; t<2; t++)
		{
		  hPairMassCut[i][j][k][t] = (TH1F*)fin->Get(Form("%s_%smass_%sCut_Pt%d_%s",dataType[0],pairType[t],varName[i],j,pidType[k]));
		  hPairMassCut[i][j][k][t]->Rebin(4);

		  hPairMass[i][j][k][t] = (TH1F*)fin->Get(Form("%s_%smass_%s_Pt%d_%s",dataType[0],pairType[t],varName[i],j,pidType[k]));
		  hPairMass[i][j][k][t]->Rebin(4);
		  if(t==0)
		    {
		      hPairMass[i][j][k][t]->SetMarkerStyle(20);
		      hPairMass[i][j][k][t]->SetMarkerColor(2);
		      hPairMass[i][j][k][t]->SetLineColor(2);
		    }
		  hPairMass[i][j][k][t]->SetMaximum(1.3*hPairMass[i][j][k][t]->GetMaximum());
		}
	    }
	  if(j%4==0)
	    {
	      c->cd(j/4+1);
	      hPairMass[i][j][pidIndex][0]->DrawCopy("PE");
	      
	      if(j/4==2)
		{
		  TLegend *leg2 = new TLegend(0.12,0.65,0.32,0.85);
		  leg2->SetBorderSize(0);
		  leg2->SetFillColor(0);
		  leg2->SetTextSize(0.045);
		}
	      // draw signal region
	      double binwidth = hPairMass[i][j][pidIndex][0]->GetBinWidth(1);
	      for(int itmp=0; itmp<3; itmp++)
		{
		  TH1F *htmp = new TH1F(Form("%s_tmp%d",hPairMass[i][j][pidIndex][0]->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
		  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		    {
		      int jbin = hPairMass[i][j][pidIndex][0]->FindFixBin(htmp->GetBinCenter(ibin));
		      htmp->SetBinContent(ibin,hPairMass[i][j][pidIndex][0]->GetBinContent(jbin));
		    }
		  if(itmp==0) htmp->SetFillColor(kGreen-6);
		  else htmp->SetFillColor(kGray);
		  htmp->SetLineColor(htmp->GetFillColor());
		  htmp->Draw("sames");
		  if(j/4==2)
		    {
		      if(itmp==0) leg2->AddEntry(htmp,"Signal region","F");
		      if(itmp==1) leg2->AddEntry(htmp,"Side-band","F");
		    }
		}
	      
	      hPairMass[i][j][pidIndex][0]->DrawCopy("samesPE");
	      hPairMass[i][j][pidIndex][1]->DrawCopy("samesHIST");
	      TPaveText *t1 = GetTitleText(Form("%s: %1.1f < p_{T}^{#mu} < %1.1f GeV/c",run_type,xbins[j],xbins[j+1]),0.06);
	      t1->Draw();
	      if(j/4==2) leg2->Draw();

	      if(j/4==0)
		{
		  TPaveText *t1 = GetPaveText(0.15,0.5,0.8,0.85,0.045,62);
		  t1->SetTextColor(4);
		  if(pidIndex==0) t1->AddText("No muon PID cuts");
		  if(pidIndex==1) t1->AddText(Form("muon PID cuts except for %s",varTitle[i]));
		  t1->Draw();
		  TLegend *leg1 = new TLegend(0.12,0.15,0.32,0.3);
		  leg1->SetBorderSize(0);
		  leg1->SetFillColor(0);
		  leg1->SetTextSize(0.045);
		  leg1->AddEntry(hPairMass[i][j][pidIndex][0],"Unlike-sign","PL");
		  leg1->AddEntry(hPairMass[i][j][pidIndex][1],"Like-sign","L");
		  leg1->Draw();
		}
	      t1 = GetPaveText(0.7,0.8,0.65,0.7,0.045,62);
	      t1->AddText(Form("scale = %4.2f",hScale[i][pidIndex]->GetBinContent(j+1)));
	      t1->Draw();
	    }
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_massULvsLS_%s.pdf",run_type,dataType[0],varName[i],pidType[pidIndex]));
    }

  double all[nbins], all_err[nbins], acc[nbins], acc_err[nbins];
  double all_1[nbins], all_err_1[nbins], acc_1[nbins], acc_err_1[nbins];
  
  //==============================================
  // Method 1: use Jpsi counts
  //==============================================
  TH1F *hULmass, *hLSmass;
  TH1F *hMuonEffInvMassFit[nVar][2];
  for(int i=0; i<nVar; i++)
    {
      for(int k=0; k<2; k++)
	{
	  hMuonEffInvMassFit[i][k] = new TH1F(Form("%s_JpsiMuon_%s_Eff_FitInvMass_%s",dataType[0],varName[i],pidType[k]), "", nbins, xbins);
	  for(int m=0; m<2; m++)
	    {
	      TCanvas *c = new TCanvas(Form("FitInvMass_%s_%s_%d",varName[i],pidType[k],m), Form("FitInvMass_%s_%s_Cut%d",varName[i],pidType[k],m), 1100, 700);
	      c->Divide(4,4);
	      for(int j=0; j<nbins; j++)
		{
		  if(m==0)
		    {
		      hULmass = (TH1F*)hPairMass[i][j][k][0]->Clone(Form("%s_clone",hPairMass[i][j][k][0]->GetName()));
		      hLSmass = (TH1F*)hPairMass[i][j][k][1]->Clone(Form("%s_clone",hPairMass[i][j][k][1]->GetName()));
		    }
		  else
		    {
		      hULmass = (TH1F*)hPairMassCut[i][j][k][0]->Clone(Form("%s_clone",hPairMassCut[i][j][k][0]->GetName()));
		      hLSmass = (TH1F*)hPairMassCut[i][j][k][1]->Clone(Form("%s_clone",hPairMassCut[i][j][k][1]->GetName()));
		    }
		  hULmass->Add(hLSmass, -1);
		  c->cd(j+1);
		  hULmass->SetMarkerStyle(24);
		  hULmass->SetMarkerColor(1);
		  hULmass->SetLineColor(1);
		  TF1 *func = new TF1(Form("func_%s",hULmass->GetName()),"gausn(0)+pol1(3)",2.5, 4.0);
		  func->SetParameter(1, 3.09);
		  func->SetParameter(2, 0.05);
		  if(run_type=="Run15_pp200" && j==7)
		    func->SetParameter(2, 0.1);
		  hULmass->Fit(func, "IR0Q");
		  hULmass->GetXaxis()->SetRangeUser(2.5, 4);
		  hULmass->SetMaximum(1.5*hULmass->GetMaximum());
		  if(j==0 || j==4) hULmass->SetMaximum(1.2*hULmass->GetMaximum());
		  hULmass->Draw();
		  func->SetLineColor(2);
		  func->SetLineWidth(2);
		  func->Draw("sames");
		  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",xbins[j],xbins[j+1]), 0.08);
		  t1->Draw();
		  if(m==0)
		    {
		      all[j] = func->GetParameter(0);
		      all_err[j] = func->GetParError(0);
		      if(all[j]<0) all[j] = 0;
		    }
		  else
		    {
		      acc[j] = func->GetParameter(0);
		      acc_err[j] = func->GetParError(0);
		      if(acc[j]<0) acc[j] = 0;
		    }
		}
	      c->cd(1);
	      t1 = GetPaveText(0.15,0.45,0.75,0.85,0.08);
	      t1->AddText(run_type);
	      t1->SetTextColor(6);
	      t1->Draw();
	      c->cd(5);
	      t1 = GetPaveText(0.15,0.4,0.7,0.85,0.08);
	      if(m==0) t1->AddText(Form("w/o %s cut",varName[i]));
	      if(m==1) t1->AddText(Form("w/  %s cut",varName[i]));
	      t1->AddText(Form("%s",pidType[k]));
	      t1->SetTextColor(6);
	      t1->Draw();
	      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_FitInvMass_cut%d_%s.pdf",run_type,dataType[0],varName[i],m,pidType[k]));
	    }
	  GetEfficiencyCurve(hMuonEffInvMassFit[i][k], all, all_err, acc, acc_err);
	}
    }

  //==============================================
  // Method 2: use single muon distribution
  //==============================================
  TH1F *hPairDis[nVar][nbins][2][2];
  TH1F *hMuonDis[nData][nVar][nbins][2];
  for(int i=0; i<nVar; i++)
    {
      TCanvas *c = new TCanvas(Form("JpsiMuon_%s",varName[i]), Form("JpsiMuon_%s",varName[i]), 1100, 700);
      c->Divide(2,2);
      for(int j=0; j<nbins; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      for(int t=0; t<2; t++)
		{
		  hPairDis[i][j][k][t] = (TH1F*)fin->Get(Form("%s_%s_%s_Pt%d_%s",dataType[0],pairType[t],varName[i],j,pidType[k]));
		  hPairDis[i][j][k][t]->SetMarkerStyle(20+t*4);
		  hPairDis[i][j][k][t]->SetMarkerColor(1+t*3);
		  hPairDis[i][j][k][t]->SetLineColor(1+t*3);
		  hPairDis[i][j][k][t]->SetMaximum(1.3*hPairDis[i][j][k][t]->GetMaximum());
		}
	      hMuonDis[0][i][j][k] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_Pt%d_%s",dataType[0],varName[i],j,pidType[k]));
	      hMuonDis[0][i][j][k]->SetMarkerStyle(21);
	      hMuonDis[0][i][j][k]->SetMarkerColor(2);
	      hMuonDis[0][i][j][k]->SetLineColor(2);
	    }
	  if(j%4==0)
	    {
	      c->cd(j/4+1);
	      hPairDis[i][j][pidIndex][0]->DrawCopy("PE");
	      hPairDis[i][j][pidIndex][1]->DrawCopy("samesPE");
	      hMuonDis[0][i][j][pidIndex]->DrawCopy("samesPE");
	      TPaveText *t1 = GetTitleText(Form("%s: %1.1f < p_{T}^{#mu} < %1.1f GeV/c",run_type,xbins[j],xbins[j+1]),0.06);
	      t1->Draw();
	      if(j/4==0)
		{
		  TPaveText *t1 = GetPaveText(0.15,0.5,0.8,0.85,0.045,62);
		  t1->SetTextColor(4);
		  if(pidIndex==0) t1->AddText("No muon PID cuts");
		  if(pidIndex==1) t1->AddText(Form("muon PID cuts except for %s",varTitle[i]));
		  t1->Draw();
		}
	      if(j/4==1)
		{
		  TLegend *leg1 = new TLegend(0.65,0.6,0.85,0.85);
		  leg1->SetBorderSize(0);
		  leg1->SetFillColor(0);
		  leg1->SetTextSize(0.045);
		  leg1->AddEntry(hPairDis[i][j][pidIndex][0],"Unlike-sign","PL");
		  leg1->AddEntry(hPairDis[i][j][pidIndex][1],"Like-sign","PL");
		  leg1->AddEntry(hMuonDis[0][i][j][pidIndex],"UL-LS","PL");
		  leg1->Draw();
		}
	    }
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_muonULvsLS_%s.pdf",run_type,dataType[0],varName[i],pidType[pidIndex]));
    }

  //==============================================
  // Jpsi muon
  //==============================================
  TF1 *funcMuonDis[nData][nVar][nbins][2];
  TFitResultPtr ptrMuonDis[nData][nVar][nbins][2];
  TH1F *hMuonDisMean[nData][nVar][2];
  TH1F *hMuonDisSigma[nData][nVar][2];
  TH1F *hMuonEffFit[nData][nVar][2];
  TH1F *hMuonEffCount[nData][nVar][2];
  for(int d=0; d<nData; d++)
    {
      for(int i=0; i<nVar; i++)
	{
	  if(d==1 && i==nVar-1) continue;
	  for(int k=0; k<2; k++)
	    {
	      hMuonDisMean[d][i][k] = new TH1F(Form("%s_JpsiMuon_%s_FitMean_%s",dataType[d],varName[i],pidType[k]), "", nbins, xbins);
	      hMuonDisSigma[d][i][k] = new TH1F(Form("%s_JpsiMuon_%s_FitSigma_%s",dataType[d],varName[i],pidType[k]), "", nbins, xbins);
	      hMuonEffFit[d][i][k] = new TH1F(Form("%s_JpsiMuon_%s_Eff_Fitting_%s",dataType[d],varName[i],pidType[k]), "", nbins, xbins);
	      hMuonEffCount[d][i][k] = new TH1F(Form("%s_JpsiMuon_%s_Eff_BinCount_%s",dataType[d],varName[i],pidType[k]), "", nbins, xbins);

	      TCanvas *c = new TCanvas(Form("c_%d%d%d",d,i,k), Form("%s_Muon_Fit%s_%s",dataType[d],varName[i],pidType[k]), 1100, 700);
	      c->Divide(4,4);
	      for(int j=0; j<nbins; j++)
		{
		  all[j] = 0; all_err[j] = 0; acc[j] = 0; acc_err[j] = 0;
		  all_1[j] = 0; all_err_1[j] = 0; acc_1[j] = 0; acc_err_1[j] = 0;
		  if(d==1) hMuonDis[d][i][j][k] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_Pt%d_%s",dataType[d],varName[i],j,pidType[k]));
		  hMuonDis[d][i][j][k]->SetMarkerColor(1);
		  hMuonDis[d][i][j][k]->SetMarkerStyle(20);
		  hMuonDis[d][i][j][k]->SetLineColor(1);
		  if(d==0) 
		    {
		      if(i<nVar-1)hMuonDis[d][i][j][k]->Rebin(2);
		      hMuonDis[d][i][j][k]->SetMaximum(1.5*hMuonDis[d][i][j][k]->GetMaximum());
		    }
		  else
		    {
		      hMuonDis[d][i][j][k]->Rebin(2);
		      hMuonDis[d][i][j][k]->SetMaximum(1.2*hMuonDis[d][i][j][k]->GetMaximum());
		    }

		  double pt = (xbins[j]+xbins[j+1]) * 0.5;
		  if(i==0) hMuonDis[d][i][j][k]->GetXaxis()->SetRangeUser(-5, 5);
		  else if(i<3) hMuonDis[d][i][j][k]->GetXaxis()->SetRangeUser(-10*fResDzVsPt->Eval(pt), 10*fResDzVsPt->Eval(pt));
		  else  hMuonDis[d][i][j][k]->GetXaxis()->SetRangeUser(-1, 1.5);

		  if(i<3) 
		    {
		      funcMuonDis[d][i][j][k] = new TF1(Form("%s_JpsiMuon_Fit%s_Pt%d_%s",dataType[d],varName[i],j,pidType[k]), "gaus");
		      funcMuonDis[d][i][j][k]->SetParameters(hMuonDis[d][i][j][k]->GetMaximum(), 0, 10);
		      if(i==0) funcMuonDis[d][i][j][k]->SetRange(-4, 4);
		      else     funcMuonDis[d][i][j][k]->SetRange(-4*fResDzVsPt->Eval(pt), 4*fResDzVsPt->Eval(pt));
		    }
		  else    
		    {
		      if(savePlot || saveHisto)
			{
			  funcMuonDis[d][i][j][k] = new TF1(Form("%s_JpsiMuon_Fit%s_Pt%d_%s",dataType[d],varName[i],j,pidType[k]), CrystalBall, -0.5,1.0,5);
			  funcMuonDis[d][i][j][k]->SetParameters(-1.5, 0, 0.2, 1, hMuonDis[d][i][j][k]->GetMaximum()*0.8);
			}
		      else
			{
			  funcMuonDis[d][i][j][k] = new TF1(Form("%s_JpsiMuon_Fit%s_Pt%d_%s",dataType[d],varName[i],j,pidType[k]), "gaus", -0.5,1.0);
			}
		      if(j==0) funcMuonDis[d][i][j][k]->SetRange(-1.0, 1.5);
		    }
		  if(d==0 && run_type=="Run15_pp200" && j>=14) continue;
		  ptrMuonDis[d][i][j][k] = hMuonDis[d][i][j][k]->Fit(funcMuonDis[d][i][j][k], "I0QLRS");
		  c->cd(j+1);
		  hMuonDis[d][i][j][k]->Draw();
		  funcMuonDis[d][i][j][k]->SetLineStyle(2);
		  funcMuonDis[d][i][j][k]->SetLineColor(2);
		  funcMuonDis[d][i][j][k]->Draw("sames");
		  TPaveText *t1 = GetTitleText(Form("%s: %1.1f < p_{T}^{#mu} < %1.1f GeV/c",dataType[d],xbins[j],xbins[j+1]),0.08);
		  t1->Draw();
		  hMuonDisMean[d][i][k]->SetBinContent(j+1, funcMuonDis[d][i][j][k]->GetParameter(1));
		  hMuonDisMean[d][i][k]->SetBinError(j+1, funcMuonDis[d][i][j][k]->GetParError(1));
		  hMuonDisSigma[d][i][k]->SetBinContent(j+1, funcMuonDis[d][i][j][k]->GetParameter(2));
		  hMuonDisSigma[d][i][k]->SetBinError(j+1, funcMuonDis[d][i][j][k]->GetParError(2));

		  double all_min = 0, all_max = 0;
		  double cut_min = 0, cut_max = 0;
		  if(i==0)
		    {
		      all_min = -5; all_max = 5;
		      cut_min = -2; cut_max = 3;
		    }
		  else if(i==1)

		    {
		      all_min = -40; all_max = 40;
		      double nsigma = pt > 3 ? 3.5 : 3;
		      cut_min = nsigma * fResDyVsPt->Eval(pt) * -1;
		      cut_max = nsigma * fResDyVsPt->Eval(pt);
		    }
		  else if(i==2)
		    {
		      all_min = -40; all_max = 40;
		      double nsigma = pt > 3 ? 3.5 : 3;
		      cut_min = nsigma * fResDzVsPt->Eval(pt) * -1;
		      cut_max = nsigma * fResDzVsPt->Eval(pt);
		    }
		  else
		    {
		      all_min = -5; all_max = 5;
		      cut_min = -1;
		      cut_max = 1;
		    }
		  
		  // bin counting
		  all[j] = hMuonDis[d][i][j][k]->IntegralAndError(hMuonDis[d][i][j][k]->FindBin(all_min),
								  hMuonDis[d][i][j][k]->FindBin(all_max),
								  all_err[j]);
		  acc[j] = hMuonDis[d][i][j][k]->IntegralAndError(hMuonDis[d][i][j][k]->FindBin(cut_min),
								  hMuonDis[d][i][j][k]->FindBin(cut_max),
								  acc_err[j]);

		  // fitting
		  all_1[j] = funcMuonDis[d][i][j][k]->Integral(all_min, all_max);
		  all_err_1[j] = funcMuonDis[d][i][j][k]->IntegralError(all_min, all_max,
									funcMuonDis[d][i][j][k]->GetParameters(),
									ptrMuonDis[d][i][j][k]->GetCovarianceMatrix().GetMatrixArray());
		  acc_1[j] = funcMuonDis[d][i][j][k]->Integral(cut_min, cut_max);
		  acc_err_1[j] = funcMuonDis[d][i][j][k]->IntegralError(cut_min, cut_max,
									funcMuonDis[d][i][j][k]->GetParameters(),
									ptrMuonDis[d][i][j][k]->GetCovarianceMatrix().GetMatrixArray());

		  if(j==12) cout << all_1[j] << " > " << acc_1[j] << endl;

		}
	      GetEfficiencyCurve(hMuonEffCount[d][i][k], all, all_err, acc, acc_err);
	      GetEfficiencyCurve(hMuonEffFit[d][i][k], all_1, all_err_1, acc_1, acc_err_1);
	      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_FitMuonDis_%s.pdf",run_type,dataType[d],varName[i],pidType[k]));
	    }
	}
    }

  //==============================================
  // Compare mean & sigma
  //==============================================
  for(int i=0; i<3; i++)
    {
      TCanvas *c = new TCanvas(Form("cMean_%d",i), Form("Comp_%sMean",varName[i]), 1100, 500);
      c->Divide(2,1);
      TLegend *leg1 = new TLegend(0.15,0.7,0.35,0.88);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.035);
      for(int d=0; d<nData; d++)
	{
	  for(int k=0; k<2; k++)
	    {
	      c->cd(1);
	      if(i==0) hMuonDisMean[d][i][k]->GetYaxis()->SetRangeUser(-0.5, 1.5);
	      else     hMuonDisMean[d][i][k]->GetYaxis()->SetRangeUser(-4, 4);
	      hMuonDisMean[d][i][k]->SetTitle(Form(";p_{T}^{#mu} [GeV/c];<%s> %s",varTitle[i],varUnit[i]));
	      hMuonDisMean[d][i][k]->SetMarkerStyle(24-k*4+d);
	      hMuonDisMean[d][i][k]->SetMarkerColor(d+1);
	      hMuonDisMean[d][i][k]->SetLineColor(d+1);
	      hMuonDisMean[d][i][k]->Draw("samesPE");
	      TPaveText *t1 = GetTitleText(Form("%s: mean of %s distribution",run_type,varTitle[i]),0.04);
	      t1->Draw();
	      

	      c->cd(2);
	      if(i==0) hMuonDisSigma[d][i][k]->GetYaxis()->SetRangeUser(0, 3);
	      else     hMuonDisSigma[d][i][k]->GetYaxis()->SetRangeUser(0, 20);
	      hMuonDisSigma[d][i][k]->SetTitle(Form(";p_{T}^{#mu} [GeV/c];#sigma(%s) %s",varTitle[i],varUnit[i]));
	      hMuonDisSigma[d][i][k]->SetMarkerStyle(24-k*4+d);
	      hMuonDisSigma[d][i][k]->SetMarkerColor(d+1);
	      hMuonDisSigma[d][i][k]->SetLineColor(d+1);
	      hMuonDisSigma[d][i][k]->Draw("samesPE");
	      TPaveText *t1 = GetTitleText(Form("%s: width of %s distribution",run_type,varTitle[i]),0.04);
	      t1->Draw();

	      if(k==0) leg1->AddEntry(hMuonDisMean[d][i][k],Form("%s: no muon PID",dataType[d]),"P");
	      if(k==1) leg1->AddEntry(hMuonDisMean[d][i][k],Form("%s: muon PID except for %s",dataType[d],varTitle[i]),"P");
	    }
	}
      c->cd(2);
      leg1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataVsEmbed_%s_CompMean.pdf",run_type,varName[i]));
    }

  //==============================================
  // check correlation
  //==============================================
  for(int i=0; i<3; i++)
    {
      TCanvas *c = new TCanvas(Form("cEff_Cor_%d",i), Form("Correlation_%sEff",varName[i]), 800, 600);
      TLegend *leg1 = new TLegend(0.55,0.2,0.7,0.4);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.04);
      leg1->SetHeader(dataType[1]);
      for(int k=0; k<2; k++)
	{
	  hMuonEffFit[1][i][k]->SetMarkerStyle(20+k*4);
	  hMuonEffFit[1][i][k]->SetMarkerColor(k+1);
	  hMuonEffFit[1][i][k]->SetLineColor(k+1);
	  hMuonEffFit[1][i][k]->GetYaxis()->SetRangeUser(0.95, 1.02);
	  hMuonEffFit[1][i][k]->Draw("samesPE");
	  leg1->AddEntry(hMuonEffFit[1][i][k], Form("Fit: %s",pidType[k]), "P");	  
	}
      TPaveText *t1 = GetTitleText(Form("%s: %s efficiency for muon",run_type,varTitle[i]));
      t1->Draw();
      leg1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Embed_%s_Eff_CheckCor.pdf",run_type,varName[i]));
    }

  //==============================================
  // Compare efficiency
  //==============================================
  TFile *ftap = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "read");
  TFile *ftof = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type), "read");
  TGraphErrors *gDtofEff[4];
  for(int i=0; i<nVar; i++)
    {
      if(i<nVar-1) gDtofEff[i] = (TGraphErrors*)ftap->Get(Form("Data_TagAndProbe_%sEff",varName[i]));
      else         gDtofEff[i] = (TGraphErrors*)ftof->Get(Form("TagAndProbe_Muon_Dtof1.00Eff"));

      TCanvas *c = new TCanvas(Form("cEff_%d",i), Form("Comp_%sEff",varName[i]), 800, 600);
      TLegend *leg1 = new TLegend(0.15,0.15,0.35,0.3);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.035);
      TLegend *leg2 = new TLegend(0.5,0.2,0.75,0.3);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.035);
      for(int d=0; d<nData; d++)
	{
	  if(i==nVar-1 && d==nData-1) continue;

	  hMuonEffCount[d][i][1]->SetTitle(Form(";p_{T}^{#mu} [GeV/c];Eff.",varTitle[i]));
	  hMuonEffCount[d][i][1]->SetMarkerStyle(24+d);
	  hMuonEffCount[d][i][1]->SetMarkerColor(d+1);
	  hMuonEffCount[d][i][1]->SetLineColor(d+1);
	  hMuonEffCount[d][i][1]->GetYaxis()->SetRangeUser(0.9, 1.05);
	  hMuonEffCount[d][i][1]->Draw("samesPE");
	  if(d==0) leg1->AddEntry(hMuonEffCount[d][i][1], Form("%s: bin-count %s dis",dataType[d],varTitle[i]), "P");
	  else     leg2->AddEntry(hMuonEffCount[d][i][1], Form("%s: bin-count %s dis",dataType[d],varTitle[i]), "P");

	  hMuonEffFit[d][i][1]->SetMarkerStyle(20+d);
	  hMuonEffFit[d][i][1]->SetMarkerColor(d+1);
	  hMuonEffFit[d][i][1]->SetLineColor(d+1);
	  hMuonEffFit[d][i][1]->DrawCopy("samesPE");
	  if(d==0) leg1->AddEntry(hMuonEffFit[d][i][1], Form("%s: fit %s dis",dataType[d],varTitle[i]), "P");	  
	  else     leg2->AddEntry(hMuonEffFit[d][i][1], Form("%s: fit %s dis",dataType[d],varTitle[i]), "P");

	  if(d==0)
	    {
	      hMuonEffInvMassFit[i][1]->SetMarkerStyle(33);
	      hMuonEffInvMassFit[i][1]->SetMarkerSize(1.5);
	      hMuonEffInvMassFit[i][1]->SetMarkerColor(4);
	      hMuonEffInvMassFit[i][1]->SetLineColor(4);
	      hMuonEffInvMassFit[i][1]->Draw("samesPE");
	      leg1->AddEntry(hMuonEffInvMassFit[i][1], Form("%s: fit InvMass",dataType[0]), "P");

	      gDtofEff[i]->SetMarkerStyle(34);
	      gDtofEff[i]->SetMarkerSize(1.5);
	      gDtofEff[i]->SetMarkerColor(6);
	      gDtofEff[i]->SetLineColor(6);
	      gDtofEff[i]->Draw("samesPEZ");
	      leg1->AddEntry(gDtofEff[i], Form("%s: tag-and-probe",dataType[0]), "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s: %s efficiency for muon",run_type,varTitle[i]));
      t1->Draw();
      leg1->Draw();
      leg2->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Embed_%s_Eff_Comp.pdf",run_type,varName[i]));

      if(iSys && i<nVar-1)
	{
	  TH1F *hSys[2];
	  const char* sysName[2] = {"Up","Down"};
	  double sysValue = 0;
	  if(run_type=="Run15_pp200")
	    {
	      if(i==0) sysValue = 0.005;
	      if(i==1) sysValue = 0.005;
	      if(i==2) sysValue = 0.005;
	    }
	  if(run_type=="Run15_pAu200")
	    {
	      if(i==0) sysValue = 0.015;
	      if(i==1) sysValue = 0.01;
	      if(i==2) sysValue = 0.01;
	    }
	  for(int s=0; s<2; s++)
	    {
	      hSys[s] = (TH1F*)hMuonEffFit[1][i][1]->Clone(Form("%s_Sys%s",hMuonEffFit[1][i][1]->GetName(),sysName[s]));
	      for(int bin=1; bin<=hSys[s]->GetNbinsX(); bin++)
		{
		  hSys[s]->SetBinContent(bin, ((1-2*s)*sysValue+1)*hSys[s]->GetBinContent(bin));
		}
	      hSys[s]->SetLineColor(2);
	      hSys[s]->SetLineStyle(2);
	      hSys[s]->SetLineWidth(2);
	      hSys[s]->Draw("samesHIST");
	    }
	  TLegend *leg3 = new TLegend(0.5,0.15,0.75,0.2);
	  leg3->SetBorderSize(0);
	  leg3->SetFillColor(0);
	  leg3->SetTextSize(0.035);
	  leg3->AddEntry(hSys[0], "Uncertainty", "L");
	  leg3->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Embed_%s_EffSys.pdf",run_type,varName[i]));
	}
    }
	  
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "update");
      for(int d=0; d<nData; d++)
	{
	  for(int i=0; i<nVar; i++)
	    {
	      if(i==nVar-1 && d==nData-1) continue;
	      for(int k=0; k<2; k++)
		{
		  hMuonEffCount[d][i][k]->Write("",TObject::kOverwrite);
		  hMuonEffFit[d][i][k]->Write("",TObject::kOverwrite);
		  if(d==0)
		    {
		      hMuonEffInvMassFit[i][k]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}

//================================================
void TagAndProbe(const int savePlot = 1, const int saveHisto = 1)
{
 // const int year = 2014;
  // const char* run_type = "Run14_AuAu200";
  // const double dtofCut = 0.75;

  gStyle->SetOptFit(0);

  const int year = 2015;
  const char* run_type = "Run15_pp200";
  const double dtofCut = 1;

  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};

  if(year==2014) 
    {
      printf("[i] Please check !\n");
      return;
    }
  if(year==2015) f = TFile::Open(Form("Output/%s.jpsi.root",run_type),"read");
  THnSparseF *hnDtof = (THnSparseF*)f->Get("mhTaP_dtof_di_mu");
    
  TH2F *hInvMassVsPtInv[4][2][2];
  const char* varName[4] = {"dTof","dY","dZ","dEdx"};

  for(int i=0; i<2; i++)
    {
      hnDtof->GetAxis(3)->SetRange(i+1,i+1);

      for(int j=0; j<4; j++)
	{
	  // cut on other variables
	  for(int l=0; l<4; l++)
	    {
	      if(l==j) continue;
	      if(l==0)
		{
		  hnDtof->GetAxis(2)->SetRangeUser(-1, 0.2);
		}
	      else
		{
		  hnDtof->GetAxis(l+3)->SetRange(2,2);
		}
	    }

	  hInvMassVsPtInv[j][0][i] = (TH2F*)hnDtof->Projection(0,1);
	  hInvMassVsPtInv[j][0][i]->SetName(Form("mhTaP_Type%d_%s_All",i,varName[j]));

	  if(j==0)
	    {
	      hnDtof->GetAxis(2)->SetRangeUser(-1*dtofCut+1e-4, dtofCut-1e-4);
	    }
	  else
	    {
	      hnDtof->GetAxis(j+3)->SetRange(2,2);
	    }
	  hInvMassVsPtInv[j][1][i] = (TH2F*)hnDtof->Projection(0,1);
	  hInvMassVsPtInv[j][1][i]->SetName(Form("mhTaP_Type%d_%s_Acc",i,varName[j]));

	  for(int l=0; l<4; l++)
	    {
	      if(l==0)
		{
		  hnDtof->GetAxis(2)->SetRange(0,-1);
		}
	      else
		{
		  hnDtof->GetAxis(l+3)->SetRange(0,-1);
		}
	    }
	}
      hnDtof->GetAxis(3)->SetRange(0,-1);
    }

  // Get the mean pT
  TH1F *hMeanPt = new TH1F("hMeanPt", "Mean p_{T} of muons", nbins, 0, nbins);
  TH1F *hMuonPt = (TH1F*)hInvMassVsPtInv[0][0][0]->ProjectionX("hMuonPt");
  for(int bin=1; bin<=nbins; bin++)
    {
      hMuonPt->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
      hMeanPt->SetBinContent(bin, hMuonPt->GetMean());
      hMeanPt->SetBinError(bin, hMuonPt->GetMeanError());
    }

  //==============================================
  // Fit to invariant mass
  //==============================================
  TH1F *hInvMass[4][2][nbins];
  TF1 *funcInvMass[4][2][nbins];
  TH1F *hJpsiCounts[4][2];
  double xmin, xmax;
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<2; k++)
	{
	  hJpsiCounts[i][k] = new TH1F(Form("hJpsiCounts_%s_Cut%d",varName[i],k), "# of J/psi", nbins, xbins);
	  TCanvas *c = new TCanvas(Form("FitInvMass_%s_Cut%d",varName[i],k), Form("FitInvMass_%s_Cut%d",varName[i],k), 1100, 700);
	  c->Divide(3,2);
	  for(int j=0; j<nbins; j++)
	    {
	      int low_bin = hInvMassVsPtInv[i][k][0]->GetXaxis()->FindFixBin(xbins[j]+1e-4);
	      int up_bin  = hInvMassVsPtInv[i][k][0]->GetXaxis()->FindFixBin(xbins[j+1]-1e-4);
	      hInvMass[i][k][j] = (TH1F*)hInvMassVsPtInv[i][k][0]->ProjectionY(Form("hInvMass_bin%d_%s_Cut%d",j+1,varName[i],k),low_bin, up_bin);
	      hInvMass[i][k][j]->SetTitle("");
	      hInvMass[i][k][j]->Sumw2();
	      TH1F *hLS = (TH1F*)hInvMassVsPtInv[i][k][1]->ProjectionY(Form("hInvMass_bin%d_%s_Cut%d_UL",j+1,varName[i],k),low_bin, up_bin);
	      hLS->Sumw2();
	      hInvMass[i][k][j]->Add(hLS, -1);
	      if(j<3) hInvMass[i][k][j]->Rebin(5);
	      else if(j<5)  hInvMass[i][k][j]->Rebin(5);
	      else hInvMass[i][k][j]->Rebin(10);

	      funcInvMass[i][k][j] = new TF1(Form("func_%s",hInvMass[i][k][j]->GetName()),"gausn(0)+pol1(3)",2.5,4.0);
	      funcInvMass[i][k][j]->SetParameter(0, hInvMass[i][k][j]->GetMaximum());
	      funcInvMass[i][k][j]->SetParameter(1, 3.096);
	      funcInvMass[i][k][j]->SetParameter(2, 0.1);
	      if(i>0)
		{
		  funcInvMass[i][k][j]->SetParameters(funcInvMass[0][k][j]->GetParameters());
		  funcInvMass[i][k][j]->SetParameter(1, funcInvMass[0][k][j]->GetParameter(1));
		  funcInvMass[i][k][j]->SetParameter(2, funcInvMass[0][k][j]->GetParameter(2));
		}
	      funcInvMass[i][k][j]->SetParName(0, "N");
	      funcInvMass[i][k][j]->SetParName(1, "#mu");
	      funcInvMass[i][k][j]->SetParName(2, "#sigma");

	      hInvMass[i][k][j]->GetXaxis()->SetRangeUser(2.5, 4.0);
	      hInvMass[i][k][j]->Fit(funcInvMass[i][k][j], "IR0S");
	      c->cd(j+1);
	      hInvMass[i][k][j]->SetMaximum(1.5*hInvMass[i][k][j]->GetMaximum());
	      hInvMass[i][k][j]->SetMarkerStyle(24);
	      hInvMass[i][k][j]->Draw();
	      funcInvMass[i][k][j]->SetLineStyle(2);
	      funcInvMass[i][k][j]->SetLineColor(2);
	      funcInvMass[i][k][j]->GetRange(xmin, xmax);
	      TF1 *tmpBkg = new TF1(Form("TmpBkg_%d_%d",i,j), "pol1", xmin, xmax);
	      for(int ipar=0; ipar<tmpBkg->GetNpar(); ipar++)
		{
		  tmpBkg->SetParameter(ipar, funcInvMass[i][k][j]->GetParameter(ipar+3));
		}
	      tmpBkg->SetLineStyle(2);
	      tmpBkg->SetLineColor(4);
	      tmpBkg->Draw("sames");
	      funcInvMass[i][k][j]->Draw("sames");
	      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c (%s)",xbins[j],xbins[j+1],varName[i]), 0.05);
	      t1->Draw();
	      double bin_width = hInvMass[i][k][j]->GetBinWidth(1);
	      hJpsiCounts[i][k]->SetBinContent(j+1, fabs(funcInvMass[i][k][j]->GetParameter(0))/bin_width);
	      hJpsiCounts[i][k]->SetBinError(j+1, fabs(funcInvMass[i][k][j]->GetParError(0))/bin_width);
	      t1 = GetPaveText(0.6,0.8,0.7,0.85,0.05);
	      t1->AddText(Form("#chi^{2}/NDF = %3.1f/%d",funcInvMass[i][k][j]->GetChisquare(), funcInvMass[i][k][j]->GetNDF()));
	      t1->AddText(Form("# of J/psi = %1.0f#pm%1.0f",hJpsiCounts[i][k]->GetBinContent(j+1),hJpsiCounts[i][k]->GetBinError(j+1)));
	      t1->Draw();
	    }
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_FitInvMass_%s_cut%d.pdf",run_type,varName[i],k));
	}
    }
 

  //==============================================
  // extract efficiency
  //==============================================
  TH1F *hJpsiCountsClone[4][2];
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<2; k++)
	{
	  hJpsiCountsClone[i][k] = (TH1F*)hJpsiCounts[i][k]->Clone(Form("%s_clone",hJpsiCounts[i][k]->GetName()));
	  if(k==1)
	    {
	      for(int bin=1; bin<=nbins; bin++)
		{
		  if(hJpsiCounts[i][k]->GetBinContent(bin)>hJpsiCounts[i][0]->GetBinContent(bin))
		    {
		      hJpsiCountsClone[i][k]->SetBinContent(bin, hJpsiCounts[i][0]->GetBinContent(bin));
		    }
		}
	    }
	}
    }

  //return;
  double x, y;
  TGraphErrors *gDtofEff[4];
  for(int i=0; i<4; i++)
    {
      gDtofEff[i] = new TGraphErrors(nbins);
      gDtofEff[i]->SetName(Form("Data_TagAndProbe_%sEff",varName[i]));
      TGraphAsymmErrors *gEff = new TGraphAsymmErrors(hJpsiCountsClone[i][1], hJpsiCountsClone[i][0],"cl=0.683 b(1,1) mode");
      for(int ipoint=0; ipoint<nbins; ipoint++)
	{
	  gEff->GetPoint(ipoint,x,y);
	  double err_h = gEff->GetErrorYhigh(ipoint);
	  double err_l = gEff->GetErrorYlow(ipoint);
	  double err = err_h > err_l ? err_h : err_l;
	  if(hJpsiCounts[i][1]->GetBinContent(ipoint+1)>hJpsiCountsClone[i][0]->GetBinContent(ipoint+1))
	    {
	      y = getEffAboveOne( hJpsiCounts[i][1]->GetBinContent(ipoint+1), hJpsiCounts[i][1]->GetBinError(ipoint+1),
				  hJpsiCounts[i][0]->GetBinContent(ipoint+1), hJpsiCounts[i][0]->GetBinError(ipoint+1),
				  err);
	    }
	  cout << i << "  " << ipoint << "  " << y << endl;
	  gDtofEff[i]->SetPoint(ipoint, hMeanPt->GetBinContent(ipoint+1), y);
	  gDtofEff[i]->SetPointError(ipoint, hMeanPt->GetBinError(ipoint+1), err);
	}
      gDtofEff[i]->SetMarkerStyle(20+i/2*4+i%2);
      gDtofEff[i]->SetMarkerColor(color[i]);
      gDtofEff[i]->SetLineColor(color[i]);
      gDtofEff[i]->SetMarkerSize(1.3);
    }
 
  TH1F *hplot = new TH1F("hplot",Form("%s: muon PID efficiencies;p_{T,#mu} (GeV/c);Eff.",run_type), 100, 0, 10);
  hplot->GetXaxis()->SetRangeUser(1,8);
  hplot->GetYaxis()->SetRangeUser(0.6, 1.1);
  c = draw1D(hplot);
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("With other PID cuts applied");
  TString legName[4] = {"#Deltatof cut","#Deltay cut","#Deltaz cut","n#sigma_{#pi} cut"};
  for(int i=3; i>-1; i--)
    {
      gDtofEff[i]->Draw("samesPEZ");
      leg->AddEntry(gDtofEff[i], legName[i].Data(), "PL");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_JpsiEffInv.pdf",run_type));

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "update");
      for(int i=0; i<4; i++)
	{
	  gDtofEff[i]->SetTitle("");
	  gDtofEff[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void makeHistoTn(const int saveHisto = 1)
{
  const int dataType = 1;
  const char* typeName[2] = {"Data","Embed"};
  
  char* run_type = "Run15_pAu200";
  
  TFile *fin = NULL;
  if(dataType==0) fin = TFile::Open(Form("output/%s.jpsi.root",run_type), "read");
  if(dataType==1) fin = TFile::Open("Rootfiles/Run15_ppAu.embed.muonPid.root", "read");

  const int nVar = 4;
  const char* varName[nVar] = {"dEdx", "dY", "dZ", "dTof"};

  const int nbins = 16;
  const double xbins[nbins+1] = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
  
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};

  const char* pairType[2] = {"UL", "LS"};
  const char* pidType[2] = {"noPID", "wPID"};

  TH1F *hPairDis[nVar][nbins][2][2];
  TH1F *hPairMass[nVar][nbins][2][2];
  TH1F *hPairMassCut[nVar][nbins][2][2];
  TH1F *hScale[nVar][2];
  TH1F *hMuonDis[nVar][nbins][2];

  if(dataType==0)
    {
      THnSparseF *hn = (THnSparseF*)fin->Get("mhJpsiMuonPid_di_mu");
      for(int j=0; j<nbins; j++)
	{
	  double pt_min = xbins[j];
	  double pt_max = xbins[j+1];
	  double pt = (pt_min+pt_max)/2;
	  hn->GetAxis(1)->SetRangeUser(pt_min, pt_max);
	  for(int t=0; t<2; t++)
	    {
	      hn->GetAxis(6)->SetRange(t+1, t+1);
	      for(int i=0; i<nVar; i++)
		{
		  hPairMass[i][j][0][t] = (TH1F*)hn->Projection(0);
		  hPairMass[i][j][0][t]->Sumw2();
		  hPairMass[i][j][0][t]->SetName(Form("%s_%smass_%s_Pt%d_noPID",typeName[dataType],pairType[t],varName[i],j));
		  cout << Form("%s_%smass_%s_Pt%d_noPID",typeName[dataType],pairType[t],varName[i],j) << endl;

		  double cut_min_def = 0, cut_max_def = 0;
		  if(i==0) { cut_min_def = -2; cut_max_def = 3; }
		  else if(i==1)
		    {
		      double nsigma = pt > 3 ? 3.5 : 3;
		      cut_min_def = nsigma * fResDyVsPt->Eval(pt) * -1;
		      cut_max_def = nsigma * fResDyVsPt->Eval(pt);
		    }
		  else if(i==2)
		    {
		      double nsigma = pt > 3 ? 3.5 : 3;
		      cut_min_def = nsigma * fResDzVsPt->Eval(pt) * -1;
		      cut_max_def = nsigma * fResDzVsPt->Eval(pt);
		    }
		  else { cut_min_def = -1; cut_max_def = 1; }

		  hn->GetAxis(i+2)->SetRangeUser(cut_min_def, cut_max_def);
		  hPairMassCut[i][j][0][t] = (TH1F*)hn->Projection(0);
		  hPairMassCut[i][j][0][t]->Sumw2();
		  hPairMassCut[i][j][0][t]->SetName(Form("%s_%smass_%sCut_Pt%d_noPID",typeName[dataType],pairType[t],varName[i],j));
		  hn->GetAxis(i+2)->SetRange(0, -1);

		  hn->GetAxis(0)->SetRangeUser(min_mass[0], max_mass[0]);
		  hPairDis[i][j][0][t] = (TH1F*)hn->Projection(i+2);
		  hPairDis[i][j][0][t]->Sumw2();
		  hPairDis[i][j][0][t]->SetName(Form("%s_%s_%s_Pt%d_noPID",typeName[dataType],pairType[t],varName[i],j));
		  hn->GetAxis(0)->SetRange(0, -1);
		  double cut_min = 0, cut_max = 0;
		  for(int index=2; index<6; index++)
		    {
		      if(index==i+2) continue;
		      if(index==2) { cut_min = 0; cut_max = 2; }
		      if(index==3) { cut_min = -2*fResDyVsPt->Eval(pt); cut_max = 2*fResDyVsPt->Eval(pt); }
		      if(index==4) { cut_min = -2*fResDzVsPt->Eval(pt); cut_max = 2*fResDzVsPt->Eval(pt); }
		      if(index==5) { cut_min = -1; cut_max = 0.2; }
		      hn->GetAxis(index)->SetRangeUser(cut_min, cut_max);
		    }
		  hPairMass[i][j][1][t] = (TH1F*)hn->Projection(0);
		  hPairMass[i][j][1][t]->Sumw2();
		  hPairMass[i][j][1][t]->SetName(Form("%s_%smass_%s_Pt%d_wPID",typeName[dataType],pairType[t],varName[i],j));
		  hn->GetAxis(i+2)->SetRangeUser(cut_min_def, cut_max_def);
		  hPairMassCut[i][j][1][t] = (TH1F*)hn->Projection(0);
		  hPairMassCut[i][j][1][t]->Sumw2();
		  hPairMassCut[i][j][1][t]->SetName(Form("%s_%smass_%sCut_Pt%d_wPID",typeName[dataType],pairType[t],varName[i],j));
		  hn->GetAxis(i+2)->SetRange(0, -1);

		  hn->GetAxis(0)->SetRangeUser(min_mass[0], max_mass[0]);
		  hPairDis[i][j][1][t] = (TH1F*)hn->Projection(i+2);
		  hPairDis[i][j][1][t]->Sumw2();
		  hPairDis[i][j][1][t]->SetName(Form("%s_%s_%s_Pt%d_wPID",typeName[dataType],pairType[t], varName[i],j));
		  hn->GetAxis(0)->SetRange(0, -1);
		  for(int index=2; index<6; index++) hn->GetAxis(index)->SetRange(0, -1);
		}
	    }
	}

      // Get scale factor
      for(int i=0; i<nVar; i++)
	{
	  for(int k=0; k<2; k++)
	    {
	      hScale[i][k] = new TH1F(Form("%s_%s_scaleFactor_%s",typeName[dataType],varName[i],pidType[k]), "", nbins, xbins); 
	      for(int j=0; j<nbins; j++)
		{
		  int bin11 = hPairMass[i][j][k][0]->FindBin(min_mass[1]);
		  int bin12 = hPairMass[i][j][k][0]->FindBin(max_mass[1]);
		  int bin21 = hPairMass[i][j][k][0]->FindBin(min_mass[2]);
		  int bin22 = hPairMass[i][j][k][0]->FindBin(max_mass[2]);

		  double ul = hPairMass[i][j][k][0]->Integral(bin11, bin12)+hPairMass[i][j][k][0]->Integral(bin21, bin22);
		  double ls = hPairMass[i][j][k][1]->Integral(bin11, bin12)+hPairMass[i][j][k][1]->Integral(bin21, bin22);
		  double scale = ul/ls;
		  hScale[i][k]->SetBinContent(j+1, scale);
		  hMuonDis[i][j][k] = (TH1F*)hPairDis[i][j][k][0]->Clone(Form("%s_JpsiMuon_%s_Pt%d_%s",typeName[dataType],varName[i],j,pidType[k]));
		  hMuonDis[i][j][k]->Add(hPairDis[i][j][k][1], -1*scale);
		}
	    }
	}
    }

  if(dataType==1)
    {
      THnSparseF *hn = NULL;
      if(run_type=="Run15_pp200") hn = (THnSparseF*)fin->Get("hnppMuonInfo");
      if(run_type=="Run15_pAu200") hn = (THnSparseF*)fin->Get("hnpAuMuonInfo");
      
      for(int j=0; j<nbins; j++)
	{
	  double pt_min = xbins[j];
	  double pt_max = xbins[j+1];
	  hn->GetAxis(0)->SetRangeUser(pt_min, pt_max);
	  for(int i=0; i<3; i++)
	    {
	      int axis = i;
	      if(i==0) axis = 3;
	      hMuonDis[i][j][0] = (TH1F*)hn->Projection(axis);
	      hMuonDis[i][j][0]->SetName(Form("%s_JpsiMuon_%s_Pt%d_%s",typeName[dataType],varName[i],j,pidType[0]));

	      for(int index=1; index<4; index++)
		{
		  if(index==axis) continue;
		  if(index==1) hn->GetAxis(4)->SetRangeUser(-2, 2);
		  if(index==2) hn->GetAxis(5)->SetRangeUser(-2, 2);
		  if(index==3) hn->GetAxis(3)->SetRangeUser(0,  2);
		}
	      hMuonDis[i][j][1] = (TH1F*)hn->Projection(axis);
	      hMuonDis[i][j][1]->SetName(Form("%s_JpsiMuon_%s_Pt%d_%s",typeName[dataType],varName[i],j,pidType[1]));
	      for(int index=3; index<6; index++) hn->GetAxis(index)->SetRange(0, -1);
	    }
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type), "update");
      if(dataType==0)
	{
	  for(int i=0; i<nVar; i++)
	    {
	      for(int j=0; j<nbins; j++)
		{
		  for(int k=0; k<2; k++)
		    {
		      for(int t=0; t<2; t++)
			{
			  hPairDis[i][j][k][t]->SetTitle("");
			  hPairDis[i][j][k][t]->Write("",TObject::kOverwrite);
			  hPairMass[i][j][k][t]->SetTitle("");
			  hPairMass[i][j][k][t]->Write("",TObject::kOverwrite);
			  hPairMassCut[i][j][k][t]->SetTitle("");
			  hPairMassCut[i][j][k][t]->Write("",TObject::kOverwrite);
			}
		    }
		}
	      for(int k=0; k<2; k++)
		{
		  hScale[i][k]->SetTitle("");
		  hScale[i][k]->Write("",TObject::kOverwrite);
		}

	      for(int j=0; j<nbins; j++)
		{
		  for(int k=0; k<2; k++)
		    {
		      hMuonDis[i][j][k]->SetTitle("");
		      hMuonDis[i][j][k]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
      if(dataType==1)
	{
	  for(int i=0; i<3; i++)
	    {
	      for(int j=0; j<nbins; j++)
		{
		  for(int k=0; k<2; k++)
		    {
		      hMuonDis[i][j][k]->SetTitle("");
		      hMuonDis[i][j][k]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}


//================================================
void makeHisto(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fdata = 0x0;
  //if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
  if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.JpsiMuon.dtof0.4.root","read");
  if(year==2015) fdata = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");


  //==============================================
  // compare invariant mass
  //==============================================
  const int nHistos = 6;
  const int nbins = 5;
  const double xbins[nbins+1] = {0,1.0,2.0,3.0,5.0,10.0};
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  double counts[nHistos][2][3] = {0};
  TH2F *hInvMassVsPtUL[nHistos];
  TH2F *hInvMassVsPtLS[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hInvMassVsPtUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_UL_di_mu",name[i]));
      hInvMassVsPtUL[i]->Sumw2();
      hInvMassVsPtLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_LS_di_mu",name[i]));
      hInvMassVsPtLS[i]->Sumw2();
    }
  
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("InvMass_%s",name[i]),Form("InvMass_%s",name[i]),1100,700);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(run_type);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hInvMassVsPtUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hInvMassVsPtUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *hInvMassUL = (TH1F*)hInvMassVsPtUL[i]->ProjectionY(Form("Data_%s_InvMassUL_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassUL->SetMarkerStyle(20);
	  hInvMassUL->SetMarkerStyle(20);
	  hInvMassUL->Rebin(5);
	  hInvMassUL->SetMaximum(1.5*hInvMassUL->GetMaximum());
	  hInvMassUL->SetMinimum(0);
	  hInvMassUL->GetXaxis()->SetRangeUser(2.2,4);
	  hInvMassUL->SetTitle("");
	  if(bin==1) leg->AddEntry(hInvMassUL,"Unlike-sign","PL");

	  TH1F *hInvMassLS = (TH1F*)hInvMassVsPtLS[i]->ProjectionY(Form("Data_%s_InvMassLS_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassLS->SetMarkerStyle(24);
	  hInvMassLS->SetMarkerColor(2);
	  hInvMassLS->SetLineColor(2);
	  hInvMassLS->Rebin(5);
	  if(bin==1) leg->AddEntry(hInvMassLS,"Like-sign","PL");

	  c->cd(bin);
	  hInvMassUL->Draw("P");
	  // draw signal region
	  double binwidth = hInvMassUL->GetBinWidth(1);
	  for(int itmp=0; itmp<3; itmp++)
	    {
	      counts[i][0][itmp] = 0;
	      counts[i][1][itmp] = 0;
	      TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMassUL->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	      for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		{
		  int jbin = hInvMassUL->FindFixBin(htmp->GetBinCenter(ibin));
		  htmp->SetBinContent(ibin,hInvMassUL->GetBinContent(jbin));
		  counts[i][0][itmp] += hInvMassUL->GetBinContent(jbin);
		  counts[i][1][itmp] += hInvMassLS->GetBinContent(jbin);
		}
	      if(itmp==0) htmp->SetFillColor(7);
	      else htmp->SetFillColor(5);
	      htmp->Draw("sames");
	      if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	      if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	    }
	  hInvMassLS->Draw("samesP");
	  TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbins[bin-1],xbins[bin]),0.05);
	  t1->Draw();
	  double count_ul = counts[i][0][1] + counts[i][0][2];
	  double count_ls = counts[i][1][1] + counts[i][1][2];
	  double scale = count_ul/count_ls;
	  double error = scale * sqrt(1/count_ul + 1/count_ls);
	  TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("S/B = 1:%4.2f",counts[i][1][0]/(counts[i][0][0]-counts[i][1][0])));
	  t1->AddText(Form("UL/LS = %4.2f#pm%4.3f",scale,error));
	  t1->Draw();
	}
      c->cd(6);
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%s_InvMassLSvsUL.pdf",run_type,name[i]));
    }

  //==============================================
  // Muon distribution UL vs LS
  //==============================================
  const int rebins[nHistos] = {8,8,10,1,1,1};
  const double minimum[nHistos] = {-50, -50, -5, 0, 0, 0};
  const double maximum[nHistos] = {50, 50, 5, 2, 50, 50};

  TH2F *hJpsiMassVsMuonPtUL[nHistos];
  TH2F *hJpsiMassVsMuonPtLS[nHistos];
  TH1F *hMuonPtSideUL[nHistos];
  TH1F *hMuonPtSideLS[nHistos];
  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];
  TCanvas *c = new TCanvas(Form("MuonSB_UL_vs_LS"),Form("MuonSB_UL_vs_LS"),1100,700);
  c->Divide(3,2);
  for(int i=0; i<nHistos; i++)
    {
      // Get the muon distribution in the side-band region
      hJpsiMassVsMuonPtUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_%s_UL_di_mu",name[i]));
      hJpsiMassVsMuonPtUL[i]->Sumw2();
      hMuonPtSideUL[i] = (TH1F*)hJpsiMassVsMuonPtUL[i]->ProjectionX(Form("hMuonPt_%s_SideBand_UL",name[i]));
      hMuonPtSideUL[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPtUL[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPtUL[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  TH1F *htmp = (TH1F*)hJpsiMassVsMuonPtUL[i]->ProjectionX(Form("hMuonPt_%s_SideBand%d_UL",name[i],j),low_bin,up_bin);
	  hMuonPtSideUL[i]->Add(htmp);
	}

      hJpsiMassVsMuonPtLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_%s_LS_di_mu",name[i]));
      hJpsiMassVsMuonPtLS[i]->Sumw2();
      hMuonPtSideLS[i] = (TH1F*)hJpsiMassVsMuonPtLS[i]->ProjectionX(Form("hMuonPt_%s_SideBand_LS",name[i]));
      hMuonPtSideLS[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPtLS[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPtLS[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  TH1F *htmp = (TH1F*)hJpsiMassVsMuonPtLS[i]->ProjectionX(Form("hMuonPt_%s_SideBand%d_LS",name[i],j),low_bin,up_bin);
	  hMuonPtSideLS[i]->Add(htmp);
	}
      c->cd(i+1);
      gPad->SetLogy();
      hMuonPtSideUL[i]->GetXaxis()->SetRangeUser(0,10);
      hMuonPtSideUL[i]->SetMarkerStyle(20);
      hMuonPtSideUL[i]->Draw("PE");
      hMuonPtSideLS[i]->SetMarkerStyle(24);
      hMuonPtSideLS[i]->SetMarkerColor(2);
      hMuonPtSideLS[i]->SetLineColor(2);
      hMuonPtSideLS[i]->Draw("samesPE");

      hDataUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMuon_%s_UL_di_mu",name[i]));
      hDataUL[i]->Sumw2();
      hDataLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMuon_%s_LS_di_mu",name[i]));
      hDataLS[i]->Sumw2();
      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("mhJpsiMuon%sVsPt",name[i]));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);
    }

  TH1F *hUL[nHistos][nMuonPtBin];
  TH1F *hLS[nHistos][nMuonPtBin];
  TH1F *hMuonFineBin[nHistos][nMuonPtBin];
  TH1F *hMuon[nHistos][nMuonPtBin];
  double mean_pt[nHistos][nMuonPtBin];
  double mean_pt_err[nHistos][nMuonPtBin];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      // calculate mean pt of each bin
      TH1F *htmp = (TH1F*)hDataDisVsPt[i]->ProjectionX(Form("%s_tmp",hDataDisVsPt[i]->GetName()),1,hDataDisVsPt[i]->GetNbinsX());
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  htmp->GetXaxis()->SetRangeUser(muonPtBins[bin-1]+1e-4, muonPtBins[bin]-1e-4);
	  mean_pt[i][bin-1] = htmp->GetMean();
	  mean_pt_err[i][bin-1] = htmp->GetMeanError();

	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(muonPtBins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(muonPtBins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("Data_JpsiMuon_%s_UL_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetTitle("");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("Data_JpsiMuon_%s_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);

	  start_bin = hMuonPtSideUL[i]->GetXaxis()->FindBin(muonPtBins[bin-1]+1e-4);
	  end_bin   = hMuonPtSideUL[i]->GetXaxis()->FindBin(muonPtBins[bin]-1e-4);
	  double scale = hMuonPtSideUL[i]->Integral(start_bin,end_bin)/hMuonPtSideLS[i]->Integral(start_bin,end_bin);
	  hMuonFineBin[i][bin-1] = (TH1F*)hUL[i][bin-1]->Clone(Form("Data_JpsiMuon_%s_bin%d",name[i],bin));
	  hMuonFineBin[i][bin-1]->Add(hLS[i][bin-1],-1.0*scale);

	  hMuon[i][bin-1] = (TH1F*)hMuonFineBin[i][bin-1]->Clone(Form("Data_JpsiMuon_%s_bin%d_Rebin",name[i],bin));
	  if(i==0 && bin>=5) hMuon[i][bin-1]->Rebin(4);
	  else hMuon[i][bin-1]->Rebin(rebins[i]);

	  c->cd(bin);
	  hUL[i][bin-1]->Rebin(rebins[i]);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
      	  hUL[i][bin-1]->Draw("P");
	  hLS[i][bin-1]->Rebin(rebins[i]);
	  hLS[i][bin-1]->Draw("samesP");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",muonPtBins[bin-1],muonPtBins[bin]),0.06);
          t1->Draw();
	  t1 = GetPaveText(0.6,0.8,0.75,0.8,0.05);
	  t1->AddText(Form("Scale = %4.3f",scale));
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.15,0.62,0.35,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(run_type);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%s_ULvsLS.pdf",run_type,name[i]));
    }

  for(int i=0; i<nHistos; i++)
    {
      hDataDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      hDataDisVsPt[i]->GetYaxis()->SetRangeUser(minimum[i], maximum[i]);
      c = draw2D(hDataDisVsPt[i],Form("%s: %s distribution for J/#psi muon;p_{T,#mu} (GeV/c);%s%s",run_type,title[i],title[i],unit[i]));
      if(i<2)
	{
	  TF1 *func11 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_upBound",fResVsPt[i]->GetName()));
	  func11->SetRange(0.5,ptbound);
	  func11->SetParameter(0,nsigma1*func11->GetParameter(0));
	  func11->SetParameter(1,nsigma1*func11->GetParameter(1));
	  func11->SetLineColor(2);
	  func11->Draw("sames");
	  TF1 *func12 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_lowBound",fResVsPt[i]->GetName()));
	  func12->SetRange(0.5,ptbound);
	  func12->SetParameter(0,-nsigma1*func12->GetParameter(0));
	  func12->SetParameter(1,-nsigma1*func12->GetParameter(1));
	  func12->SetLineColor(2);
	  func12->Draw("sames");
	  if(ptbound<10)
	    {
	      TF1 *func21 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_upBound",fResVsPt[i]->GetName()));
	      func21->SetRange(ptbound,10);
	      func21->SetParameter(0,nsigma2*func21->GetParameter(0));
	      func21->SetParameter(1,nsigma2*func21->GetParameter(1));
	      func21->SetLineColor(2);
	      func21->Draw("sames");
	      TF1 *func22 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_lowBound",fResVsPt[i]->GetName()));
	      func22->SetRange(ptbound,10);
	      func22->SetParameter(0,-nsigma2*func22->GetParameter(0));
	      func22->SetParameter(1,-nsigma2*func22->GetParameter(1));
	      func22->SetLineColor(2);
	      func22->Draw("sames");
	    }
	}
      else if(i==2)
	{
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      else
	{
	  TLine *line = GetLine(0,cuts[i],10,cuts[i]);
	  line->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sVsPt.pdf",run_type,name[i]));
    }
 
  //==============================================
  // Fit data distributions to extract efficiency
  //==============================================
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  TF1 *func[nHistos][nMuonPtBin];
  TFitResultPtr ptr[nHistos][nMuonPtBin];
  for(int i=0; i<3; i++)
    {
      hFitDataMean[i] = new TH1F(Form("Data_JpsiMuon_%s_FitMean",name[i]),Form("Mean of %s;p_{T} (GeV/c)",title[i]),nMuonPtBin,muonPtBins);
      hFitDataSigma[i] = new TH1F(Form("Data_JpsiMuon_%s_FitSigma",name[i]),Form("Width of %s;p_{T} (GeV/c)",title[i]),nMuonPtBin,muonPtBins);

      TCanvas *c = new TCanvas(Form("Fit_%s",name[i]),Form("Fit_%s",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  func[i][bin-1] = new TF1(Form("Data_JpsiMuon_%sFit_bin%d",name[i],bin),"gaus",minimum[i], maximum[i]);
	  if(i==1 && bin==1) ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"R0QS");
	  else ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"IR0QS");
	  hFitDataMean[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(1));
	  hFitDataMean[i]->SetBinError(bin,func[i][bin-1]->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,func[i][bin-1]->GetParError(2));

	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	  hFit->Draw();
	  func[i][bin-1]->SetLineColor(4);
	  func[i][bin-1]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",muonPtBins[bin-1],muonPtBins[bin]),0.06);
	  t1->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sFit.pdf",run_type,name[i]));
    }
  
  // efficiency
  TGraphAsymmErrors *gCountDataEff[nHistos];
  TGraphAsymmErrors *gFitDataEff[nHistos];
  for(int i=0; i<3; i++)
    {
      TH1F* hBase1  = new TH1F(Form("hBase1_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hMatch1 = new TH1F(Form("hMatch1_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hBase2  = new TH1F(Form("hBase2_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hMatch2 = new TH1F(Form("hMatch2_%d",i), "", nMuonPtBin, muonPtBins);
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  double min = -999, max = -999;
	  double pt = mean_pt[i][bin-1];;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else
	    { min = -1; max = 3; }

	  // bin counting
	  double all_err;
	  int low_bin = hMuonFineBin[i][bin-1]->FindFixBin(minimum[i]+1e-4);
	  int high_bin = hMuonFineBin[i][bin-1]->FindFixBin(maximum[i]-1e-4);
	  double all = hMuonFineBin[i][bin-1]->IntegralAndError(low_bin, high_bin,all_err);
	  double acc_err;
	  low_bin = hMuonFineBin[i][bin-1]->FindFixBin(min+1e-4);
	  high_bin = hMuonFineBin[i][bin-1]->FindFixBin(max-1e-4);
	  double acc = hMuonFineBin[i][bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
	  double acc_corr = acc - (hMuonFineBin[i][bin-1]->GetBinContent(low_bin)/hMuonFineBin[i][bin-1]->GetBinWidth(low_bin)*(min-hMuonFineBin[i][bin-1]->GetXaxis()->GetBinLowEdge(low_bin)) + hMuonFineBin[i][bin-1]->GetBinContent(high_bin)/hMuonFineBin[i][bin-1]->GetBinWidth(high_bin)*(hMuonFineBin[i][bin-1]->GetXaxis()->GetBinUpEdge(high_bin)-max) );
	  //double acc_corr = acc;
	  cout << acc_corr/all << "  " << acc/all << endl;
	  hBase1->SetBinContent(bin,all);
	  hBase1->SetBinError(bin,all_err);
	  hMatch1->SetBinContent(bin,acc_corr);
	  hMatch1->SetBinError(bin,acc_corr*acc_err/acc);

	  // fitting
	  all = func[i][bin-1]->Integral(minimum[i],maximum[i]);
	  all_err = func[i][bin-1]->IntegralError(minimum[i],maximum[i],func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  acc = func[i][bin-1]->Integral(min,max);
	  acc_err = func[i][bin-1]->IntegralError(min,max,func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  hBase2->SetBinContent(bin,all);
	  hBase2->SetBinError(bin,all_err);
	  hMatch2->SetBinContent(bin,acc);
	  hMatch2->SetBinError(bin,acc_err);
	}

      // bin counting    
      TH1F *hMatch1_corr = (TH1F*)hMatch1->Clone(Form("%_tmp",hMatch1->GetName()));
      for(int ibin=1; ibin<=hMatch1->GetNbinsX(); ibin++)
	{
	  if(hMatch1->GetBinContent(ibin)>hBase1->GetBinContent(ibin))
	    {
	      hMatch1_corr->SetBinContent(ibin,hBase1->GetBinContent(ibin));
	    }
	}
      gCountDataEff[i] = new TGraphAsymmErrors(hMatch1_corr, hBase1,"cl=0.683 b(1,1) mode");
      gCountDataEff[i]->SetName(Form("Data_JpsiMuon_%sEff_BinCounting",name[i]));
      double x,y;
      hMatch1->Divide(hBase1);
      for(int ipoint=0; ipoint<gCountDataEff[i]->GetN(); ipoint++)
	{
	  gCountDataEff[i]->GetPoint(ipoint,x,y);
	  gCountDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	  gCountDataEff[i]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	  gCountDataEff[i]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);

	  if(hMatch1->GetBinContent(ipoint+1)>1)
	    {
	      gCountDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],hMatch1->GetBinContent(ipoint+1));	 
	      gCountDataEff[i]->SetPointEYhigh(ipoint,hMatch1->GetBinError(ipoint+1));
	      gCountDataEff[i]->SetPointEYlow(ipoint,hMatch1->GetBinError(ipoint+1));   
	    }
	}

      // fitting 
      gFitDataEff[i] = new TGraphAsymmErrors(hMatch2, hBase2,"cl=0.683 b(1,1) mode");
      gFitDataEff[i]->SetName(Form("Data_JpsiMuon_%sEff_Fitting",name[i]));
      for(int ipoint=0; ipoint<gFitDataEff[i]->GetN(); ipoint++)
	{
	  gFitDataEff[i]->GetPoint(ipoint,x,y);
	  gFitDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	  gFitDataEff[i]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	  gFitDataEff[i]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
	}

      // Counting vs fitting
      TCanvas* cEff = new TCanvas(Form("%sEff_CountVsFit",name[i]),Form("%sEff_CountVsFit",name[i]),800,600);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      gPad->SetGridy();
      gCountDataEff[i]->SetMarkerStyle(20);
      gCountDataEff[i]->GetYaxis()->SetRangeUser(0.5,1.2);
      gCountDataEff[i]->GetXaxis()->SetRangeUser(1.3,10);
      gCountDataEff[i]->SetMarkerSize(1.2);
      gCountDataEff[i]->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
      ScaleHistoTitle(gCountDataEff[i]->GetHistogram(),0.05,1.0,0.045,0.05,1.0,0.045,62);
      gCountDataEff[i]->Draw("APZ");
      gFitDataEff[i]->SetMarkerStyle(25);
      gFitDataEff[i]->SetMarkerColor(2);
      gFitDataEff[i]->SetLineColor(2);
      gFitDataEff[i]->SetMarkerSize(1.2);
      gFitDataEff[i]->Draw("PZsames");
      TPaveText *t1 = GetTitleText(Form("%s: %s efficiency",run_type,title[i]),0.04);
      t1->Draw();
      TLegend *leg = new TLegend(0.3,0.3,0.5,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(run_type);
      leg->AddEntry(gCountDataEff[i],Form("Bin counting"),"p");
      leg->AddEntry(gFitDataEff[i],Form("Fitting"),"p");
      leg->Draw();
      if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sEff_FitVsCount.pdf",run_type,name[i]));
    }

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
      for(int i=0; i<nHistos; i++)
	{
	  hDataDisVsPt[i]->Write(Form("Data_JpsiMuon_%sVsPt",name[i]),TObject::kOverwrite);
	  for(int bin=1; bin<=nMuonPtBin; bin++)
	    {
	      hMuonFineBin[i][bin-1]->Write("",TObject::kOverwrite);
	      hMuon[i][bin-1]->Write("",TObject::kOverwrite);
	      if(i<3)
		func[i][bin-1]->Write("",TObject::kOverwrite);
	    }
	  if(i<3)
	    {
	      hFitDataMean[i]->Write("",TObject::kOverwrite);
	      hFitDataSigma[i]->Write("",TObject::kOverwrite);
	      gCountDataEff[i]->Write("",TObject::kOverwrite);
	      gFitDataEff[i]->Write("",TObject::kOverwrite);
	    }
	}
    }
}


//================================================
void Run16nSigmaPiEff(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  const int nbins = 6;
  const double xbins[nbins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5, 10};
  TF1* funcSigma = new TF1("func_MuonSigma","pol0",1.3,10);
  funcSigma->SetParameter(0, 0.9985);
  funcSigma->SetParError(0, 0.02662);
  TF1* funcMean[2];
  funcMean[0] = new TF1("func_MuonMean","[0]+[1]*exp([2]/x)",1.3,10);
  funcMean[0]->SetParameters(-4.841, 5.1, 0.07498);
  double errors[3] = {41.72, 41.65, 0.5859};
  funcMean[0]->SetParErrors(errors);
  funcMean[1] = new TF1("func_MuonMean_sys","pol0",1.3,10);
  funcMean[1]->SetParameter(0, 0.4309);
  funcMean[1]->SetParError(0, 0.03729);
  TH1F *hEff[2];
  TF1 *funcEff[2];
  for(int i=0; i<1; i++)
    {
      hEff[i] = new TH1F(Form("hEff_%d",i),Form("hEff_%d",i),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  double pt = hEff[i]->GetBinCenter(bin);
	  TF1* functmp = new TF1(Form("functmp_%d_%d",i,bin),"gaus",-5,5);
	  functmp->SetParameters(1, funcMean[i]->Eval(pt), funcSigma->Eval(pt));
	  double eff = functmp->Integral(-1,3)/functmp->Integral(-5,5);
	  hEff[i]->SetBinContent(bin, eff);
	  hEff[i]->SetBinError(bin, 3e-3);
	}
      funcEff[i] = new TF1(Form("funcEff_%d",i),"[0]+[1]*exp([2]/x)",1.3,10);
      funcEff[i]->SetParameters(0.4588, 0.4374, 0.1118);
      hEff[i]->Fit(funcEff[i], "R0Q");
      hEff[i]->GetYaxis()->SetRangeUser(0.8,1.0);
      hEff[i]->SetMarkerStyle(21);
    }
  TCanvas *c = draw1D(hEff[0],"Run16_AuAu200: efficiency of -1 < n#sigma_{#pi} <3 cut for single muons;p_{T} (GeV/c)");
  funcEff[0]->SetLineColor(2);
  funcEff[0]->SetLineStyle(2);
  funcEff[0]->Draw("sames");
  TLegend *leg = new TLegend(0.45,0.25,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEff[0],"Efficiency via fitting","P");
  leg->AddEntry(funcEff[0],"Fit to efficiency","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run16_AuAu200/ana_JpsiMuon/DataJpsiMuon_nSigmaPiEff.pdf"));
  
}




