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
//#include "/Users/admin/Work/ALICE/Utility/drawHistos.h"
//#include "/Users/admin/Work/STAR/util/defs.h"

#include <cmath>
using namespace std;
const int year = YEAR;

TFile *f;

//================================================
void ana_JpsiMuon()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }

  //DeltaTof();
  MtdVpdTacDiff();
  //kink();
}

//================================================
void MtdVpdTacDiff(const Int_t savePlot = 0, const int saveHisto = 0)
{
  TList *list = new TList;

  const char *name = "AuAu200";
  const char *title = "Run14_AuAu200";
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum = 760;
  const double maximum = 840;
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  TH1F *htmp = 0x0;

  TFile *fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");

  //==============================================
  // single muon UL vs LS
  //==============================================
  TH2F *hJpsiMassVsMuonPt[2];
  TH1F *hMuonPtSide[2];
  hJpsiMassVsMuonPt[0] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_MtdVpdTacDiff_UL_di_mu"));
  hJpsiMassVsMuonPt[1] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_MtdVpdTacDiff_LS_di_mu"));
  for(int i=0; i<2; i++)
    {
      // Get the muon distribution in the side-band region
      hJpsiMassVsMuonPt[i]->Sumw2();
      hMuonPtSide[i] = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_MtdVpdTacDiff_SideBand_%d",i));
      hMuonPtSide[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  htmp = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_MtdVpdTacDiff_SideBand%d_%d",j,i),low_bin,up_bin);
	  hMuonPtSide[i]->Add(htmp);
	}
    }


  TH2F *hDataUL = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_UL_di_mu");
  hDataUL->Sumw2();
  TH2F *hDataLS = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_LS_di_mu");
  hDataLS->Sumw2();
  TH2F *hDataDisVsPt = (TH2F*)hDataUL->Clone(Form("JpsiMuonMtdVpdTacDiff"));
  hDataDisVsPt->Add(hDataLS, -1);
  

  const int rebin = 2;
  TH1F *hUL[nbins];
  TH1F *hLS[nbins];
  TH1F *hMuon[nbins];
  double mean_pt[nbins];
  double mean_pt_err[nbins];
  TCanvas *c = new TCanvas(Form("MtdVpdTacDiff_UL_vs_LS"),Form("MtdVpdTacDiff_UL_vs_LS"),1100,700);
  c->Divide(3,2);

  // calculate mean pt of each bin
  int bin1 = hDataDisVsPt->GetYaxis()->FindBin(minimum);
  int bin2 = hDataDisVsPt->GetYaxis()->FindBin(maximum);
  htmp = (TH1F*)hDataDisVsPt->ProjectionX(Form("%s_tmp",hDataDisVsPt->GetName()),bin1,bin2);
  for(int bin=1; bin<=nbins; bin++)
    {
      htmp->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
      mean_pt[bin-1] = htmp->GetMean();
      mean_pt_err[bin-1] = htmp->GetMeanError();

      c->cd(bin);
      int start_bin = hDataUL->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hDataUL->GetXaxis()->FindBin(xbins[bin]-1e-4);

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",name,bin),start_bin,end_bin);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->Rebin(rebin);
      hUL[bin-1]->SetMaximum(1.5*hUL[bin-1]->GetMaximum());
      hUL[bin-1]->GetXaxis()->SetRangeUser(minimum,maximum);
      hUL[bin-1]->SetTitle("");
      hUL[bin-1]->Draw("P");

      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",name,bin),start_bin,end_bin);
      hLS[bin-1]->SetMarkerStyle(24);
      hLS[bin-1]->SetMarkerColor(2);
      hLS[bin-1]->SetLineColor(2);
      hLS[bin-1]->Rebin(rebin);
      hLS[bin-1]->Draw("samesP");

      start_bin = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      end_bin   = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin]-1e-4);
      double scale = hMuonPtSide[0]->Integral(start_bin,end_bin)/hMuonPtSide[1]->Integral(start_bin,end_bin);

      hMuon[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("%s_DataMtdVpdTacDiff_bin%d",name,bin));
      hMuon[bin-1]->Add(hLS[bin-1],-1);

      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
      t1 = GetPaveText(0.6,0.8,0.75,0.8,0.05);
      t1->AddText(Form("Scale = %4.3f",scale));
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.15,0.63,0.35,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(title);
  leg->AddEntry(hUL[0],"Unlike-sign","PL");
  leg->AddEntry(hLS[0],"Like-sign","PL");
  leg->Draw();

  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdVpdTacDiff_ULvsLS.pdf",run_type,name));
  return;
  // Fit data to extract efficiency
  const int nLumi = 2;
  const char *name_lumi[nLumi] = {"prod_high","prod_low"};
  const double min[nHistos][nLumi] = { {789,786}, {910,910} };
  const double max[nHistos][nLumi] = { {837,837}, {960,960} };
  const double shift[nLumi] = {0, 130};
  TF1 *funcFitData[nHistos][nbins];
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  TF1 *funcBay[nHistos][nLumi][nbins];
  TCanvas *cFuncBay[nHistos][nLumi];
  TH1F *hFitDataEff[nHistos][nLumi];
  TGraphAsymmErrors *gFitDataEff[nHistos][nLumi];
  double x,y;
  for(int i=0; i<nHistos; i++)
    {
      hFitDataMean[i] = new TH1F(Form("%s_JpsiMuon_MtdVpdTacDiff_FitMean",name[i]),Form("%s: mean of MtdVpdTacDiff;p_{T} (GeV/c)",title[i]),nbins,xbins);
      hFitDataSigma[i] = new TH1F(Form("%s_JpsiMuon_MtdVpdTacDiff_FitSigma",name[i]),Form("%s: sigma of MtdVpdTacDiff;p_{T} (GeV/c)",title[i]),nbins,xbins);

      TH1F *hBase = (TH1F*)hFitDataMean[i]->Clone(Form("hBase_%d",i));
      TH1F *hMatch[nLumi];
      for(int k=0; k<nLumi; k++)
	{
	  hMatch[k] = (TH1F*)hFitDataMean[i]->Clone(Form("hMatch_%d_%d",i,k));
	}

      TCanvas *c = new TCanvas(Form("%s_FitMtdVpdTacDiff",name[i]),Form("%s_FitMtdVpdTacDiff",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  funcFitData[i][bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",min[i][0],maximum[i]);
	  funcFitData[i][bin-1]->SetParameter(2,5);
	  hFit->Fit(funcFitData[i][bin-1],"IR0Q");
	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->Draw();
	  funcFitData[i][bin-1]->SetLineColor(4);
	  funcFitData[i][bin-1]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  double all = funcFitData[i][bin-1]->Integral(minimum[i],maximum[i]);
	  double all_err = funcFitData[i][bin-1]->IntegralError(minimum[i],maximum[i]);
	  hBase->SetBinContent(bin,all);
	  hBase->SetBinError(bin,all_err);
	  
	  for(int k=0; k<nLumi; k++)
	    {
	      double acc = funcFitData[i][bin-1]->Integral(min[i][k],max[i][k]);
	      double acc_err = funcFitData[i][bin-1]->IntegralError(min[i][k],max[i][k]);
	      hMatch[k]->SetBinContent(bin,acc);
	      hMatch[k]->SetBinError(bin,acc_err);
	    }
	  hFitDataMean[i]->SetBinContent(bin,funcFitData[i][bin-1]->GetParameter(1)-shift[i]);
	  hFitDataMean[i]->SetBinError(bin,funcFitData[i][bin-1]->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,funcFitData[i][bin-1]->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,funcFitData[i][bin-1]->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.8,0.85,0.055);
      t1->AddText(title[i]);
      t1->SetTextColor(2);
      t1->SetTextFont(62);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_FitMtdVpdTacDiff.pdf",run_type,name[i]));

      for(int k=0; k<nLumi; k++)
	{
	  hFitDataEff[i][k] = DivideTH1ForEff(hMatch[k],hBase,Form("%s_JpsiMuon_MtdVpdTacDiff_FitEff_%s",name[i],name_lumi[k]));
	  gFitDataEff[i][k] = new TGraphAsymmErrors(hMatch[k], hBase,"cl=0.683 b(1,1) mode");
	  gFitDataEff[i][k]->SetName(Form("graph_%s_JpsiMuon_MtdVpdTacDiff_FitEff_%s",name[i],name_lumi[k]));
	  cFuncBay[i][k] = new TCanvas(Form("cFuncBay_%s_%s",name[i],name_lumi[k]),Form("cFuncBay_%s_%s",name[i],name_lumi[k]),800,600);
	  TH1F *htmp = new TH1F(Form("htmp_%s_%s",name[i],name_lumi[k]),";#epsilon;P(#epsilon)",1000,0.7,1.0);
	  if(i==0) htmp->GetYaxis()->SetRangeUser(0,0.5);
	  if(i==1) htmp->GetYaxis()->SetRangeUser(0,1.5);
	  htmp->Draw();
	  TLegend *leg = new TLegend(0.15,0.5,0.5,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  leg->SetHeader(Form("%s: %s",name[i],name_lumi[k]));
	  for(int ipoint=0; ipoint<nbins; ipoint++)
	    {
	      gFitDataEff[i][k]->GetPoint(ipoint,x,y);
	      gFitDataEff[i][k]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	      gFitDataEff[i][k]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	      gFitDataEff[i][k]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
	  
	      double pw  = hMatch[k]->GetBinContent(ipoint+1);
	      double pw2 = hMatch[k]->GetSumw2()->At(ipoint+1);
	      double tw  = hBase->GetBinContent(ipoint+1);
	      double tw2 = hBase->GetSumw2()->At(ipoint+1);
	      double alpha = 1, beta = 1;
	      double norm = tw/tw2;
	      double N = tw * norm;
	      double p = pw * norm;
	      double aa = pw * norm + alpha;
	      double bb = (tw-pw) * norm + beta;
	      double eff = TEfficiency::BetaMode(aa,bb);
	      double lower, upper;
	      double conf = 0.683;
	      TEfficiency::BetaShortestInterval(conf, aa, bb, lower, upper);
	      double gamma = 1./TMath::Beta(p+1,N-p+1);
	      funcBay[i][k][ipoint] = new TF1(Form("funcBay_%s_%s_%d",name[i],name_lumi[k],ipoint),"[0]*x**[1]*(1-x)**([2]-[1])",0,1);
	      double binwidth = (funcBay[i][k][ipoint]->GetXmax()-funcBay[i][k][ipoint]->GetXmin())/funcBay[i][k][ipoint]->GetNpx();
	      funcBay[i][k][ipoint]->SetParameters(gamma*binwidth,p,N);	

	      cFuncBay[i][k]->cd();
	      funcBay[i][k][ipoint]->SetLineColor(color[ipoint]);
	      funcBay[i][k][ipoint]->Draw("same");
	      leg->AddEntry(funcBay[i][k][ipoint],Form("p_{T} = %1.2f GeV/c",x),"l");
	    }
	  leg->Draw();
	  if(savePlot)
	    cFuncBay[i][k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_ProbDis.pdf",run_type,name[i],name_lumi[k]));
	}
    }

  list->Clear();
  TString legName[2];
  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataMean[i]);
      if(shift[i]==0) legName[i] = title[i];
      else            legName[i] = Form("%s: - %1.0f ch",title[i],shift[i]);
    }
  TCanvas *c = drawHistos(list,"MtdVpdTacDiff_Mean","Mean of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);Mean",false,0,0,true,780,810,false,true,legName,true,"",0.15,0.4,0.2,0.4,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_MtdVpdTacDiffMean.pdf",run_type));
  list->Clear();

  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataSigma[i]);
      legName[i] = title[i];
    }
  c = drawHistos(list,"MtdVpdTacDiff_Sigma","Sigma of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);#sigma",false,0,0,true,0,20,false,true,legName,true,"",0.15,0.3,0.73,0.88,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_MtdVpdTacDiffSigma.pdf",run_type));
  list->Clear();

  // compare pp vs AuAu
  c = new TCanvas(Form("Muon_pp_vs_AuAu"),Form("Muon_pp_vs_AuAu"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      c->cd(bin);
      TLegend *leg = new TLegend(0.15,0.7,0.6,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      for(int i=0; i<nHistos; i++)
	{
	  TH1F *htmp = (TH1F*)hMuon[i][bin-1]->Clone(Form("%s_clone",hMuon[i][bin-1]->GetName()));
	  htmp->Rebin(2);
	  htmp->SetMarkerStyle(21+i*4);
	  htmp->SetMarkerColor(i+1);
	  htmp->SetLineColor(i+1);
	  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	    {
	      int shift_bin = int(shift[i]/htmp->GetBinWidth(1));
	      if(shift_bin==0) continue;
	      int new_bin = ibin - shift_bin;
	      if(new_bin<=0) continue;
	      htmp->SetBinContent(new_bin, htmp->GetBinContent(ibin));
	      htmp->SetBinError(new_bin, htmp->GetBinError(ibin));
	      htmp->SetBinContent(ibin, 0);
	      htmp->SetBinError(ibin, 0);
	    }
	  htmp->Scale(1./htmp->GetBinContent(htmp->FindFixBin(hFitDataMean[0]->GetBinContent(bin))));
	  htmp->GetXaxis()->SetRangeUser(minimum[0],maximum[0]);
	  htmp->SetMaximum(1.6);
	  if(i==0) htmp->Draw();
	  else     htmp->Draw("samesP");
	  if(i==0) leg->AddEntry(htmp,title[i],"PL");
	  else     leg->AddEntry(htmp,Form("%s: - %1.0f ch",title[i],shift[i]),"PL");
	}
      if(bin==2)       leg->Draw();
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/CompMtdVpdTac_ppVsAuAu.pdf",run_type));

  
  // trigger efficiency
  TFile *fTrigEff = TFile::Open("Rootfiles/Run14.AuAu200.MuonTrigEff.root","read");
  TCanvas *cTrigEff[2];
  TF1 *funcdata[2];
  TLegend *legdata[2];
  TH1F *gFit_err[2][2];
  TF1 *fExtrap[2];
  for(int k=0; k<2; k++)
    {
      gFitDataEff[0][k]->SetMarkerStyle(21);
      gFitDataEff[0][k]->SetMarkerSize(1.2);
      gFitDataEff[0][k]->GetYaxis()->SetRangeUser(0.7,1.05);
      gFitDataEff[0][k]->GetXaxis()->SetRangeUser(1.2,8);
      cTrigEff[k] = drawGraph(gFitDataEff[0][k],"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);

      TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gFitDataEff[0][k]->Clone(Form("Fit_%s",gFitDataEff[0][k]->GetName()));
      funcdata[k] = new TF1(Form("FitFunc_%s",gFitDataEff[0][k]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
      gfit->Fit(funcdata[k],"R0");
      funcdata[k]->SetLineColor(2);
      funcdata[k]->Draw("sames");

      fExtrap[k] = (TF1*)fTrigEff->Get(Form("MuonTrigEff_cent0060_P%d",3-k*2));
      fExtrap[k]->SetLineStyle(5);
      fExtrap[k]->SetLineColor(4);
      fExtrap[k]->DrawCopy("sames");

      legdata[k] = new TLegend(0.4,0.25,0.6,0.45);
      legdata[k]->SetBorderSize(0);
      legdata[k]->SetFillColor(0);
      legdata[k]->SetTextSize(0.04);
      legdata[k]->SetHeader(Form("Run14_AuAu200, %s",name_lumi[k]));
      legdata[k]->AddEntry(gFitDataEff[0][k],"Data-driven: J/#Psi muons","P");
      legdata[k]->AddEntry(funcdata[k],"Fit to data","L");
      legdata[k]->AddEntry(fExtrap[k],"Extrapolation method","L");
      legdata[k]->Draw();
      if(savePlot) cTrigEff[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_MtdTrigEff.pdf",run_type,name[0],name_lumi[k]));
    }

  // systematic uncertainty
  const int mode = 0; // 0 - randomize efficiency; 1 - randomize MtdVpdTacDiff
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  const int nexpr = 1000;
  TGraphAsymmErrors *gEffRndm[nHistos][nLumi][nexpr];
  double rndm_all[nbins];
  double rndm_all_err[nbins];
  double rndm_acc[nbins];
  double rndm_acc_err[nbins];
  double rndm_acc_2[nbins];
  double rndm_acc_err_2[nbins];
  for(int i=0; i<1; i++)
    {
      for(int j = 0; j < nexpr; j++)
	{
	  if(mode==0)
	    {
	      TGraphAsymmErrors *gBaseEff = 0;
	      for(int k=0; k<2; k++)
		{
		  TGraphAsymmErrors *gBaseEff = gFitDataEff[i][k];
		  int npoint = gBaseEff->GetN();
		  gEffRndm[i][k][j] = new TGraphAsymmErrors(npoint);
		  gEffRndm[i][k][j]->SetName(Form("Fit%d_%s_%s",j,name[i],name_lumi[k]));
		  double *exh = gBaseEff->GetEXhigh();
		  double *exl = gBaseEff->GetEXlow();
		  double *eyh = gBaseEff->GetEYhigh();
		  double *eyl = gBaseEff->GetEYlow();
		  for(int ipoint=0; ipoint<nbins; ipoint++)
		    {
		      gBaseEff->GetPoint(ipoint,x,y);
		      y = funcBay[i][k][ipoint]->GetRandom(0,1);
		      //y = y - eyl[ipoint];
		      gEffRndm[i][k][j]->SetPoint(ipoint,x,y);
		      gEffRndm[i][k][j]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
		    }
		}
	    }
	  if(mode==1)
	    {
	      for(int bin=1; bin<=nbins; bin++)
		{
		  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit%d_%s",j,hMuon[i][bin-1]->GetName()));
		  for(int ibin=1; ibin<=hFit->GetNbinsX(); ibin++)
		    {
		      double value = hFit->GetBinContent(ibin);
		      double error = hFit->GetBinError(ibin);
		      value = rndm->Gaus(value, error);
		      hFit->SetBinContent(ibin, value);
		      hFit->SetBinError(ibin, error);
		    }
		  TF1 *func_rndm = new TF1(Form("func_%d_%d_%d",i,bin,nexpr),"gaus",min[i][0],maximum[i]);
		  func_rndm->SetParameter(2,5);
		  hFit->Fit(func_rndm,"IR0Q");
		  rndm_all[bin-1] = func_rndm->Integral(minimum[i],maximum[i]);
		  rndm_all_err[bin-1] = func_rndm->IntegralError(minimum[i],maximum[i]);
		  rndm_acc[bin-1] = func_rndm->Integral(min[i][0],max[i][0]);
		  rndm_acc_err[bin-1] = func_rndm->IntegralError(min[i][0],max[i][0]);
		  rndm_acc_2[bin-1] = func_rndm->Integral(min[i][1],max[i][1]);
		  rndm_acc_err_2[bin-1] = func_rndm->IntegralError(min[i][1],max[i][1]);
		}
	      for(int k=0; k<2; k++)
		{
		  TH1F *hBase = new TH1F(Form("gBase_%s_%s_%d",name[i],name_lumi[k],j),Form("gBase_%s_%s_%d",name[i],name_lumi[k],j),nbins,xbins);
		  TH1F *hMatch2 = new TH1F(Form("gMatch_%s_%s_%d",name[i],name_lumi[k],j),Form("gMatch_%s_%s_%d",name[i],name_lumi[k],j),nbins,xbins);
		  for(int bin=1; bin<=nbins; bin++)
		    {
		      hBase->SetBinContent(bin,rndm_all[bin-1]);
		      hBase->SetBinError(bin,rndm_all_err[bin-1]);
		      if(k==0)
			{
			  hMatch2->SetBinContent(bin,rndm_acc[bin-1]);
			  hMatch2->SetBinError(bin,rndm_acc_err[bin-1]);
			}
		      else
			{
			  hMatch2->SetBinContent(bin,rndm_acc_2[bin-1]);
			  hMatch2->SetBinError(bin,rndm_acc_err_2[bin-1]);
			}
		    }
		 gEffRndm[i][k][j] = new TGraphAsymmErrors(hMatch2, hBase,"cl=0.683 b(1,1) mode");
		 gEffRndm[i][k][j]->SetName(Form("Fit%d_%s_%s",j,name[i],name_lumi[k]));
		 for(int ipoint=0; ipoint<nbins; ipoint++)
		   {
		     gEffRndm[i][k][j]->GetPoint(ipoint,x,y);
		     gEffRndm[i][k][j]->SetPoint(ipoint,mean_pt[i][ipoint],y);
		     gEffRndm[i][k][j]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
		     gEffRndm[i][k][j]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
		   }
		}
	    }
	}
    }

  TH1F *hplot = new TH1F("hplot",";p_{T} (GeV/c);Efficiency",100,0,10);
  hplot->GetYaxis()->SetRangeUser(0.5,1.1);
  hplot->GetXaxis()->SetRangeUser(1,7);
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1,1.2,1.4,1.6,1.8,2,2.5,3,4,6,8};
  TH1F *hSpread[nHistos][2][nbins_rndm];
  TGraphErrors *gLimits[nHistos][nLumi][2];
  TF1 *funcLimits[nHistos][nLumi][2];
  for(int i=0; i<1; i++)
    {
      for(int k=0; k<2; k++)
	{
	  TCanvas *c= new TCanvas(Form("cEff_%s_%s",name[i],name_lumi[k]),Form("cEff_%s_%s",name[i],name_lumi[k]),800,600);
	  hplot->DrawCopy();
	  TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  leg->SetHeader(Form("%s, %s",title[i],name_lumi[k]));

	  TCanvas *cFit = new TCanvas(Form("cEffFit_%s_%s",name[i],name_lumi[k]),Form("cEffFit_%s_%s",name[i],name_lumi[k]),1200,800);
	  cFit->Divide(4,3);

	  for(int l=0; l<nbins_rndm; l++)
	    {
	      hSpread[i][k][l] = new TH1F(Form("hSpread_%d_%s_%d",i,name_lumi[k],l),Form("hSpread_%d_%s_%d",i,name_lumi[k],l),600,0,1.2);
	    }

	  for(int j = 0; j < nexpr; j++)
	    {
	      TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[i][k][j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
	      funcTmp->SetParLimits(0,0,1);
	      funcTmp->SetParLimits(1,0,5);
	      funcTmp->SetParLimits(2,0,1);
	      gEffRndm[i][k][j]->Fit(funcTmp,"R0Q");
	      c->cd();
	      funcTmp->SetLineStyle(2);
	      funcTmp->Draw("sames");
	      if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
	      for(int l=0; l<nbins_rndm; l++)
		{
		  hSpread[i][k][l]->Fill(funcTmp->Eval(xbins_rndm[l]));
		}
	    }

	  gLimits[i][k][0] = new TGraphErrors(nbins_rndm);
	  gLimits[i][k][1] = new TGraphErrors(nbins_rndm);
	  for(int l=0; l<nbins_rndm; l++)
	    {
	      cFit->cd(l+1);
	      if(k==0)
		{
		  hSpread[i][k][l]->Rebin(2);
		  if(l<2) hSpread[i][k][l]->GetXaxis()->SetRangeUser(0.5,1.1);
		  else if(l<4) hSpread[i][k][l]->GetXaxis()->SetRangeUser(0.8,1.1);
		  else hSpread[i][k][l]->GetXaxis()->SetRangeUser(0.85,1.1);
		}
	      else
		{
		  if(l<4) hSpread[i][k][l]->GetXaxis()->SetRangeUser(0.5,1.1);
		  else hSpread[i][k][l]->GetXaxis()->SetRangeUser(0.9,1.1);
		}
	      hSpread[i][k][l]->SetTitle(";Efficiency;");
	      hSpread[i][k][l]->Draw();
	      TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[l]),0.06);
	      t1->Draw();
	      double boundary = hSpread[i][k][l]->GetBinCenter(hSpread[i][k][l]->GetMaximumBin());
	      cout << boundary << endl;

	      TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[i][k][l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",boundary,1);
	      funcTmp2->FixParameter(1,boundary);
	      funcTmp2->SetParameter(2,0.01);
	      hSpread[i][k][l]->Fit(funcTmp2,"R0Q");
	      funcTmp2->SetLineColor(4);
	      funcTmp2->SetLineStyle(2);
	      funcTmp2->Draw("sames");
	      gLimits[i][k][1]->SetPoint(l,xbins_rndm[l],funcTmp2->GetParameter(1)+TMath::Abs(funcTmp2->GetParameter(2)));
	      gLimits[i][k][1]->SetPointError(l,0,funcTmp2->GetParError(2));


	      TF1 *funcTmp = new TF1(Form("FitFunc_%s",hSpread[i][k][l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,boundary);
	      funcTmp->FixParameter(1,funcTmp2->GetParameter(1));
	      funcTmp->SetParameter(2,0.01);
	      hSpread[i][k][l]->Fit(funcTmp,"R0Q");
	      funcTmp->SetLineColor(2);
	      funcTmp->SetLineStyle(2);
	      funcTmp->Draw("sames");
	      gLimits[i][k][0]->SetPoint(l,xbins_rndm[l],funcTmp->GetParameter(1)-TMath::Abs(funcTmp->GetParameter(2)));
	      gLimits[i][k][0]->SetPointError(l,1e-5,funcTmp->GetParError(2));
	    }
	  cFit->cd(12);
	  TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
	  t1->AddText(name[i]);
	  t1->AddText(name_lumi[k]);
	  t1->Draw();
	  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_RndmFit.pdf",run_type,name[i],name_lumi[k]));

	  c->cd();
	  leg->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_RndmCurve.pdf",run_type,name[i],name_lumi[k]));
	  for(int t=0; t<2; t++)
	    {
	      gLimits[i][k][t]->SetMarkerStyle(20);
	      gLimits[i][k][t]->SetMarkerColor(2+2*t);
	      gLimits[i][k][t]->SetLineColor(2+2*t);
	      gLimits[i][k][t]->Draw("samesP");

	      funcLimits[i][k][t] = new TF1(Form("Fit_%s",gLimits[i][k][t]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	      funcLimits[i][k][t]->SetParameters(0.95,2,0.2);
	      gLimits[i][k][t]->Fit(funcLimits[i][k][t],"R0");
	      funcLimits[i][k][t]->SetLineColor(gLimits[i][k][t]->GetMarkerColor());
	      funcLimits[i][k][t]->DrawCopy("sames");
	    }
	  leg->AddEntry(gLimits[i][k][0],"Lower limit","PL");
	  leg->AddEntry(gLimits[i][k][1],"Upper limit","PL");
	  leg->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_RndmLimits.pdf",run_type,name[i],name_lumi[k]));
	}
    }

  TCanvas *cSys[2];
  for(int k=0; k<2; k++)
    {
      TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gFitDataEff[0][k]->Clone(Form("%s_sys",gFitDataEff[0][k]->GetName()));
      cSys[k] = drawGraph(gtmp,"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);
      funcdata[k]->Draw("sames");
      fExtrap[k]->DrawCopy("sames");
      legdata[k]->Draw();
      for(int t=0; t<2; t++)
	{
	  funcLimits[0][k][t]->SetLineColor(2);
	  funcLimits[0][k][t]->SetLineStyle(2);
	  funcLimits[0][k][t]->Draw("sames");
	}
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.25);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(funcLimits[0][k][0],"Uncertainty","L");
      leg->Draw();
      if(savePlot) cSys[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_MtdTrigEffSys.pdf",run_type,name[0],name_lumi[k]));
    }

  // another way to estimate the uncertainty
  TCanvas *cSys2[2];
  TF1 *funcLimits2[nHistos][nLumi][2];
  for(int k=0; k<2; k++)
    {
      TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gFitDataEff[0][k]->Clone(Form("%s_sys2",gFitDataEff[0][k]->GetName()));
      cSys2[k] = drawGraph(gtmp,"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);
      funcdata[k]->Draw("sames");
      fExtrap[k]->DrawCopy("sames");
      legdata[k]->Draw();

      TGraphAsymmErrors *gBaseEff = gFitDataEff[0][k];
      int npoint = gBaseEff->GetN();
      for(int t=0; t<2; t++)
	{
	  TGraphAsymmErrors *gSys2 = new TGraphAsymmErrors(npoint);
	  gSys2->SetName(Form("Fit2_%s_%s_%d",name[0],name_lumi[k],t));
	  double *exh = gBaseEff->GetEXhigh();
	  double *exl = gBaseEff->GetEXlow();
	  double *eyh = gBaseEff->GetEYhigh();
	  double *eyl = gBaseEff->GetEYlow();
	  for(int ipoint=0; ipoint<nbins; ipoint++)
	    {
	      gBaseEff->GetPoint(ipoint,x,y);
	      if(t==0) y = y + eyh[ipoint];
	      if(t==1) y = y - eyl[ipoint];
	      gSys2->SetPoint(ipoint,x,y);
	      gSys2->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	    }
	  funcLimits2[0][k][t] = new TF1(Form("Fit2_%s",gSys2->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  funcLimits2[0][k][t]->SetParameters(0.95,2,0.2);
	  gSys2->Fit(funcLimits2[0][k][t],"R0Q");
	  funcLimits2[0][k][t]->SetLineColor(kMagenta-4);
	  funcLimits2[0][k][t]->SetLineStyle(2);
	  funcLimits2[0][k][t]->DrawCopy("sames");
	}
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.25);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(funcLimits2[0][k][0],"Uncertainty2","L");
      leg->Draw();
      if(savePlot) cSys2[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_MtdTrigEffSys2.pdf",run_type,name[0],name_lumi[k]));
    }
  
  // Bin counting method
  TGraphAsymmErrors *gCountDataEff[nHistos][nLumi];
  for(int i=0; i<1; i++)
    {
      TH1F *hBase = new TH1F(Form("%s_hBaseCount",name[i]),"",nbins,xbins);
      TH1F *hMatch[nLumi];
      for(int k=0; k<nLumi; k++)
	{
	  hMatch[k] = new TH1F(Form("%s_hMatchCount_%s",name[i],name_lumi[k]),"",nbins,xbins);
	}
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);

	  TH1F *hCount = (TH1F*)hDataUL[i]->ProjectionY(Form("%s_Count_UL_bin%d",name[i],bin),start_bin,end_bin);
	  TH1F *htmp   = (TH1F*)hDataLS[i]->ProjectionY(Form("%s_Count_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hCount->Add(htmp,-1);

	  double mean = funcFitData[i][bin-1]->GetParameter(1);
	  double value = 0, error = 0, fraction = 0;

	  // all
	  start_bin = hCount->FindFixBin(minimum[i]);
	  end_bin   = hCount->FindFixBin(maximum[i]);
	  int new_bin = -1;
	  for(int ibin=start_bin; ibin<=end_bin; ibin++)
	    {
	      if(hCount->GetXaxis()->GetBinLowEdge(ibin)<min[i][0])
		new_bin = hCount->FindFixBin(2*mean-hCount->GetBinCenter(ibin));
	      else
		new_bin = ibin;
	      //printf("[i] new bin = %d, old bin = %d, center = %f\n",new_bin,ibin,hCount->GetBinCenter(ibin));

	      value += hCount->GetBinContent(new_bin);
	      error += TMath::Power(hCount->GetBinError(new_bin),2);
	      //printf("[i] all bin %d: frac = %f, value = %f\n",ibin,fraction,value);
	    }
	  error = TMath::Sqrt(error);
	  hBase->SetBinContent(bin,value);
	  hBase->SetBinError(bin,error);
	  printf("[i] all: %f\n",value);

	  for(int k=0; k<nLumi; k++)
	    {
	      // accepted
	      start_bin = hCount->FindFixBin(min[i][k]);
	      end_bin   = hCount->FindFixBin(max[i][k]);
	      value = 0, error = 0, fraction = 0;
	      for(int ibin=start_bin; ibin<=end_bin; ibin++)
		{
		  if(hCount->GetXaxis()->GetBinLowEdge(ibin)<min[i][0])
		    new_bin = hCount->FindFixBin(2*mean-hCount->GetBinCenter(ibin));
		  else
		    new_bin = ibin;
		  value += hCount->GetBinContent(new_bin);
		  error += TMath::Power(hCount->GetBinError(new_bin),2);
		  //printf("[i] acc bin %d: frac = %f, value = %f\n",ibin,fraction,value);
		}
	      error = TMath::Sqrt(error);
	      hMatch[k]->SetBinContent(bin,value);
	      hMatch[k]->SetBinError(bin,error);
	      printf("[i] acc: %f\n",value);
	    }
	}
      double x,y;
      for(int k=0; k<nLumi; k++)
	{
	  gCountDataEff[i][k] = new TGraphAsymmErrors(hMatch[k], hBase,"cl=0.683 b(1,1) mode");
	  gCountDataEff[i][k]->SetName(Form("graph_%s_JpsiMuon_MtdVpdTacDiff_CountEff_%s",name[i],name_lumi[k]));
	  for(int ipoint=0; ipoint<nbins; ipoint++)
	    {
	      gCountDataEff[i][k]->GetPoint(ipoint,x,y);
	      gCountDataEff[i][k]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	      gCountDataEff[i][k]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	      gCountDataEff[i][k]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
	    }
	}
    }

  for(int k=0; k<nLumi; k++)
    {
      cSys[k]->cd();
      gCountDataEff[0][k]->SetMarkerStyle(24);
      gCountDataEff[0][k]->Draw("samesPE");
      TLegend *leg = new TLegend(0.4,0.15,0.6,0.2);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gCountDataEff[0][k],"Bin Counting","P");
      leg->Draw();
      if(savePlot) cSys[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_MtdTrigEffSys_BinCount.pdf",run_type,name[0],name_lumi[k]));
    }

  for(int k=0; k<nLumi; k++)
    {
      cSys2[k]->cd();
      gCountDataEff[0][k]->SetMarkerStyle(24);
      gCountDataEff[0][k]->Draw("samesPE");
      TLegend *leg = new TLegend(0.4,0.15,0.6,0.2);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gCountDataEff[0][k],"Bin Counting","P");
      leg->Draw();
      if(savePlot) cSys2[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_%s_MtdTrigEffSys2_BinCount.pdf",run_type,name[0],name_lumi[k]));
    }
  
}

//================================================
void DeltaTof(const Int_t savePlot = 0, const int saveHisto = 0)
{
  TList *list = new TList;

  const char *name = "AuAu200";
  const char *title = "Run14_AuAu_200";
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};

  TH2F *hInvMassVpPtUL[nHistos];
  TH2F *hInvMassVpPtLS[nHistos];
  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];
  TH2F *hDtofVsMod[nHistos];

  TFile *fdata[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      if(i==0) fdata[i] = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
      if(i==1) fdata[i] = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");
      hDataUL[i]       = (TH2F*)fdata[i]->Get("mhJpsiMuonDtof_UL_di_mu");
      hDataUL[i]->SetName(Form("%s_%d",hDataUL[i]->GetName(),i));
      hDataUL[i]->Sumw2();

      hDataLS[i]       = (TH2F*)fdata[i]->Get(Form("mhJpsiMuonDtof_LS_di_mu"));
      hDataLS[i]->SetName(Form("%s_%d",hDataLS[i]->GetName(),i));
      hDataLS[i]->Sumw2();

      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("JpsiMuonDtofVsPt_%d",i));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);

      hDtofVsMod[i] = (TH2F*)fdata[i]->Get(Form("mhDeltaTof_di_mu"));
      hDtofVsMod[i]->SetName(Form("%s_%d",hDtofVsMod[i]->GetName(),i));

      hInvMassVpPtUL[i]       = (TH2F*)fdata[i]->Get("mhJpsiMassVsPt_Dtof_UL_di_mu");
      hInvMassVpPtUL[i]->SetName(Form("%s_%s",hInvMassVpPtUL[i]->GetName(),name[i]));
      hInvMassVpPtUL[i]->Sumw2();
      hInvMassVpPtLS[i]       = (TH2F*)fdata[i]->Get("mhJpsiMassVsPt_Dtof_LS_di_mu");
      hInvMassVpPtLS[i]->SetName(Form("%s_%s",hInvMassVpPtLS[i]->GetName(),name[i]));
      hInvMassVpPtLS[i]->Sumw2();
      
    }

  for(int i=0; i<nHistos; i++)
    {
      hDtofVsMod[i]->GetYaxis()->SetRangeUser(-10,10);
      TCanvas *c = draw2D(hDtofVsMod[i],Form("%s: #Deltatof of tracks matched to MTD;Module;#Deltatof (ns)",title[i]));
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_DtofvsMod_MthTrk.pdf",run_type,name[i]));
    }

  //==============================================
  // compare invariant mass
  //==============================================
  TH1F *hInvMassUL[nHistos][nbins];
  TH1F *hInvMassLS[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_InvMass",name[i]),Form("%s_InvMass",name[i]),1100,700);
      c->Divide(3,2);
      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hInvMassVpPtUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hInvMassVpPtUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hInvMassUL[i][bin-1] = (TH1F*)hInvMassVpPtUL[i]->ProjectionY(Form("Data_InvMass_UL_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassUL[i][bin-1]->SetMarkerStyle(20);
	  hInvMassUL[i][bin-1]->SetMarkerStyle(20);
	  hInvMassUL[i][bin-1]->Rebin(10);
	  hInvMassUL[i][bin-1]->SetMaximum(1.5*hInvMassUL[i][bin-1]->GetMaximum());
	  hInvMassUL[i][bin-1]->GetXaxis()->SetRangeUser(2,4);
	  hInvMassUL[i][bin-1]->SetTitle("");
	  hInvMassUL[i][bin-1]->Draw("P");

	  hInvMassLS[i][bin-1] = (TH1F*)hInvMassVpPtLS[i]->ProjectionY(Form("Data_InvMass_LS_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassLS[i][bin-1]->SetMarkerStyle(24);
	  hInvMassLS[i][bin-1]->SetMarkerColor(2);
	  hInvMassLS[i][bin-1]->SetLineColor(2);
	  hInvMassLS[i][bin-1]->Rebin(10);
	  hInvMassLS[i][bin-1]->Draw("samesP");

	  TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hInvMassUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hInvMassLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_Dtof_InvMass_LSvsUL_InPtBins.pdf",run_type,name[i]));

      c = new TCanvas(Form("%s_InvMass_Ratio",name[i]),Form("%s_InvMass_Ratio",name[i]),1100,700);
      c->Divide(3,2);
      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  TH1F *hRatio = (TH1F*)hInvMassLS[i][bin-1]->Clone(Form("Data_InvMassRatio_%s_bin%d",name[i],bin));
	  hRatio->Divide(hInvMassUL[i][bin-1]);
	  if(i==0) hRatio->GetYaxis()->SetRangeUser(0.9,1.1);
	  if(i==1) hRatio->GetYaxis()->SetRangeUser(0.2,1.8);
	  hRatio->SetTitle("");
	  hRatio->Draw();
	  TPaveText *t1 = GetTitleText(Form("Ratio: LS/UL (%1.1f < p_{T}^{J/#psi} < %1.1f)",xbins[bin-1],xbins[bin]),0.05);
          t1->Draw();
	  TLine *line = GetLine(2,1,4,1,1);
	  line->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_Dtof_InvMass_LStoUL_InPtBins.pdf",run_type,name[i]));
    }

  //==============================================
  // UL vs LS
  //==============================================
  const int rebin = 5;
  TH1F *hUL[nHistos][nbins];
  TH1F *hLS[nHistos][nbins];
  TH1F *hMuon[nHistos][nbins];
  double mean_pt[nHistos][nbins];
  double mean_pt_err[nHistos][nbins];

  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      // calculate mean pt of each bin
      int bin1 = hDataDisVsPt[i]->GetYaxis()->FindBin(-5);
      int bin2 = hDataDisVsPt[i]->GetYaxis()->FindBin(5);
      TH1F *htmp = (TH1F*)hDataDisVsPt[i]->ProjectionX(Form("%s_tmp",name[i]),bin1,bin2);
      for(int bin=1; bin<=nbins; bin++)
	{
	  htmp->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
	  mean_pt[i][bin-1] = htmp->GetMean();
	  mean_pt_err[i][bin-1] = htmp->GetMeanError();

	  c->cd(bin);
	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("%s_DataDtof_UL_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->Rebin(rebin);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(-2,4);
	  hUL[i][bin-1]->SetTitle("");
	  hUL[i][bin-1]->Draw("P");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("%s_DataDtof_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);
	  hLS[i][bin-1]->Rebin(rebin);
	  hLS[i][bin-1]->Draw("samesP");

	  hMuon[i][bin-1] = (TH1F*)hUL[i][bin-1]->Clone(Form("%s_DataDtof_bin%d",name[i],bin));
	  hMuon[i][bin-1]->Add(hLS[i][bin-1],-1);

	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.5,0.6,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(title[i]);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_Dtof_ULvsLS_InPtBins.pdf",run_type,name[i]));
    }

  //==============================================
  // compare pp vs AuAu
  //==============================================
  TCanvas *c = new TCanvas("Muon_pp_vs_AuAu","Muon_pp_vs_AuAu",1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      c->cd(bin);
      TLegend *leg = new TLegend(0.15,0.7,0.6,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      for(int i=0; i<nHistos; i++)
	{
	  TH1F *htmp = (TH1F*)hMuon[i][bin-1]->Clone(Form("%s_clone",hMuon[i][bin-1]->GetName()));
	  htmp->Rebin(2);
	  htmp->SetMarkerStyle(21+i*4);
	  htmp->SetMarkerColor(i+1);
	  htmp->SetLineColor(i+1);
	  htmp->Scale(1./htmp->Integral());
	  htmp->GetXaxis()->SetRangeUser(-2,2);
	  htmp->SetMaximum(2*htmp->GetMaximum());
	  if(i==0) htmp->Draw();
	  else     htmp->Draw("samesP");
	  leg->AddEntry(htmp,title[i],"PL");
	}
      if(bin==1)       leg->Draw();
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/CompDtof_ppVsAuAu.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  TF1 *func[nHistos][nbins];
  TFitResultPtr ptr[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      hFitDataMean[i] = new TH1F(Form("%s_JpsiMuon_Dtof_FitMean",name[i]),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",title[i]),nbins,xbins);
      hFitDataSigma[i] = new TH1F(Form("%s_JpsiMuon_Dtof_FitSigma",name[i]),Form("%s: sigma of #Deltatof (AuAu-pp);p_{T} (GeV/c)",title[i]),nbins,xbins);
      TCanvas *c = new TCanvas(Form("%s_FitDtof",name[i]),Form("%s_FitDtof",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  cout << i << "  " << bin << endl;
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  hFit->GetXaxis()->SetRangeUser(-3,3);

	  // func[i][bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,3);
	  // func[i][bin-1]->SetParameter(2,0.1);
	  // ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"IR0QS");

          func[i][bin-1] = new TF1(Form("func_%d_%d",i,bin),CrystalBall,-3,3,5);
          func[i][bin-1]->SetParameters(-1, 0, 0.2, 50, hFit->GetMaximum());
	  if(i==1 && bin==3) func[i][bin-1]->SetParameters(-0.5, 0.5, 1.5, 1, hFit->GetMaximum());
	  if(i==0) func[i][bin-1]->SetParameters(-1.5, -0.05, 0.2, 50, hFit->GetMaximum());
	  if(i==0 && bin==6) func[i][bin-1]->SetParameters(-0.5, 0, 0.5, 100, hFit->GetMaximum());
          ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"IR0S");

	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->Draw();
	  func[i][bin-1]->SetLineColor(4);
	  func[i][bin-1]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  hFitDataMean[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(1));
	  hFitDataMean[i]->SetBinError(bin,func[i][bin-1]->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,func[i][bin-1]->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.8,0.85,0.055);
      t1->AddText(title[i]);
      t1->SetTextColor(2);
      t1->SetTextFont(62);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_JpsiMuon_FitDtof.pdf",run_type,name[i]));
    }
  

  list->Clear();
  TString legName1[2];
  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataMean[i]);
      legName1[i] = title[i];
    }
  c = drawHistos(list,"Dtof_Mean","Mean of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);Mean",false,0,0,true,-0.8,0.3,false,true,legName1,true,"",0.15,0.4,0.2,0.4,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_Dtof_FitMean.pdf",run_type));
  list->Clear();

  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataSigma[i]);
    }
  c = drawHistos(list,"Dtof_Sigma","Sigma of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);#sigma",false,0,0,true,0,0.5,false,true,legName1,true,"",0.15,0.3,0.73,0.88,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_Dtof_FitSigma.pdf",run_type));
  list->Clear();


  TH1F *hFitDataMeanDiff = new TH1F(Form("JpsiMuon_Dtof_FitMeanDiff"),Form("Difference in mean of #Deltatof;p_{T} (GeV/c)"),nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      hFitDataMeanDiff->SetBinContent(bin,hFitDataMean[0]->GetBinContent(bin)-hFitDataMean[1]->GetBinContent(bin));
      hFitDataMeanDiff->SetBinError(bin,TMath::Sqrt(TMath::Power(hFitDataMean[0]->GetBinError(bin),2)+TMath::Power(hFitDataMean[1]->GetBinError(bin),2)));
    }
  hFitDataMeanDiff->SetMarkerStyle(20);
  TF1 *funcDataMeanDiff = new TF1("func_DataMeanDiff","pol0",0,10);
  hFitDataMeanDiff->Fit(funcDataMeanDiff,"0RQ");
  hFitDataMeanDiff->GetYaxis()->SetRangeUser(-0.1,0.1);
  c = draw1D(hFitDataMeanDiff,"");
  funcDataMeanDiff->SetLineColor(4);
  funcDataMeanDiff->Draw("same");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_Dtof_FitMeanDiff.pdf",run_type));

  //==============================================
  // extract efficiency
  //==============================================
  const double minimum = -3;
  const double maximum = 3;
  const int nValue = 4;
  const double min[nValue] = {minimum,minimum,minimum,minimum};
  const double max[nValue] = {1, 0.75, 0.5, 0.25};

  // bin counting method for data
  TGraphAsymmErrors *gCountDataEff[nHistos][nValue];
  double x,y;
  for(int i=0; i<nHistos; i++)
    {
      for(int j=0; j<nValue; j++)
	{
	  // Counting method
	  TH1F *hBase = new TH1F(Form("hBase_%d_%d",i,j),Form("hBase_%d_%d",i,j),nbins,xbins);
	  TH1F *hMatch = new TH1F(Form("hMatch_%d_%d",i,j),Form("hMatch_%d_%d",i,j),nbins,xbins);
	  for(int bin=1; bin<=nbins; bin++)
	    {
	      double all_err;
	      int low_bin = hMuon[i][bin-1]->FindFixBin(minimum+1e-4);
	      int high_bin = hMuon[i][bin-1]->FindFixBin(maximum-1e-4);
	      double all = hMuon[i][bin-1]->IntegralAndError(low_bin, high_bin,all_err);

	      double acc_err;
	      low_bin = hMuon[i][bin-1]->FindFixBin(min[j]+1e-4);
	      high_bin = hMuon[i][bin-1]->FindFixBin(max[j]-1e-4);
	      double acc = hMuon[i][bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
	      hBase->SetBinContent(bin,all);
	      hBase->SetBinError(bin,all_err);
	      hMatch->SetBinContent(bin,acc);
	      hMatch->SetBinError(bin,acc_err);
	    }
	  gCountDataEff[i][j] = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
	  gCountDataEff[i][j]->SetName(Form("%s_JpsiMuon_CountEff_Dtof%1.1f",name[i],max[j]));
	  for(int ipoint=0; ipoint<nbins; ipoint++)
	    {
	      gCountDataEff[i][j]->GetPoint(ipoint,x,y);
	      gCountDataEff[i][j]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	      gCountDataEff[i][j]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	      gCountDataEff[i][j]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
	    }
	}
    }

  // fitting method
  double shift[2] = {0,funcDataMeanDiff->GetParameter(0)};
  TGraphAsymmErrors *gFitDataEff[nHistos][nValue][2]; // 0 - no shift; 1 - with shift
  for(int i=0; i<nHistos; i++)
    {
      for(int j=0; j<nValue; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      TH1F *hBase = new TH1F(Form("hBase_%d_%d_%d",i,j,k),Form("hBase_%d_%d",i,j),nbins,xbins);
	      TH1F *hMatch = new TH1F(Form("hMatch_%d_%d_%d",i,j,k),Form("hMatch_%d_%d",i,j),nbins,xbins);
	      for(int bin=1; bin<=nbins; bin++)
		{
		  double all = func[i][bin-1]->Integral(minimum,maximum);
		  double all_err = func[i][bin-1]->IntegralError(minimum,maximum,func[i][bin-1]->GetParameters(), ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
		  double acc = func[i][bin-1]->Integral(min[j],max[j]-shift[k]);
		  double acc_err = func[i][bin-1]->IntegralError(min[j],max[j]-shift[k],func[i][bin-1]->GetParameters(),  ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
		  hBase->SetBinContent(bin,all);
		  hBase->SetBinError(bin,all_err);
		  hMatch->SetBinContent(bin,acc);
		  hMatch->SetBinError(bin,acc_err);
		  //if(k==0) cout << i << " " << max[j] << "  " << bin << " -> " << acc << "/" << all << endl;
		}
	      gFitDataEff[i][j][k] = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
	      gFitDataEff[i][j][k]->SetName(Form("%s_JpsiMuon_FitEff_Dtof%1.1f_shift%1.1f",name[i],max[j],shift[k]));
	      for(int ipoint=0; ipoint<nbins; ipoint++)
		{
		  gFitDataEff[i][j][k]->GetPoint(ipoint,x,y);
		  gFitDataEff[i][j][k]->SetPoint(ipoint,mean_pt[i][ipoint],y);
		  gFitDataEff[i][j][k]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
		  gFitDataEff[i][j][k]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
		}
	    }
	}
    }


  // check pp 200: counting vs fitting
  c = new TCanvas(Form("%s_CountVsFit",name[1]),Form("%s_CountVsFit",name[1]),1100,700);
  c->Divide(2,2);
  for(int j=0; j<nValue; j++)
    {
      c->cd(j+1);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      gCountDataEff[1][j]->SetMarkerStyle(21);
      gCountDataEff[1][j]->GetYaxis()->SetRangeUser(0.6,1.05);
      gCountDataEff[1][j]->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
      ScaleHistoTitle(gCountDataEff[1][j]->GetHistogram(),0.05,1.0,0.045,0.05,1.0,0.045,62);
      gCountDataEff[1][j]->Draw("APZ");
      gFitDataEff[1][j][0]->SetMarkerStyle(24);
      gFitDataEff[1][j][0]->SetMarkerColor(2);
      gFitDataEff[1][j][0]->SetLineColor(2);
      gFitDataEff[1][j][0]->Draw("Psames");
      TPaveText *t1 = GetTitleText(Form("Efficiency for %1.0f < #DeltaTof < %1.2f ns",min[j],max[j]),0.06);
      t1->Draw();
      if(j==0)
	{
	  TLegend *leg = new TLegend(0.4,0.25,0.6,0.5);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->SetHeader(title[1]);
	  leg->AddEntry(gCountDataEff[1][j],"Bin counting","p");
	  leg->AddEntry(gFitDataEff[1][j][0],"Fitting","p");
	  leg->Draw();
	}
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/pp200_DtofEff_FitVsCount.pdf",run_type));

  // check AuAu 200: counting vs fitting
  c = new TCanvas(Form("%s_CountVsFit",name[0]),Form("%s_CountVsFit",name[0]),1100,700);
  c->Divide(2,2);
  for(int j=0; j<nValue; j++)
    {
      c->cd(j+1);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      gCountDataEff[0][j]->SetMarkerStyle(21);
      gCountDataEff[0][j]->GetYaxis()->SetRangeUser(0.6,1.05);
      gCountDataEff[0][j]->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
      ScaleHistoTitle(gCountDataEff[0][j]->GetHistogram(),0.05,1.0,0.045,0.05,1.0,0.045,62);
      gCountDataEff[0][j]->Draw("APZ");
      gFitDataEff[1][j][1]->SetMarkerStyle(24);
      gFitDataEff[1][j][1]->SetMarkerColor(2);
      gFitDataEff[1][j][1]->SetLineColor(2);
      gFitDataEff[1][j][1]->Draw("PZsames");
      gFitDataEff[0][j][0]->SetMarkerStyle(25);
      gFitDataEff[0][j][0]->SetMarkerColor(4);
      gFitDataEff[0][j][0]->SetLineColor(4);
      gFitDataEff[0][j][0]->Draw("PZsames");
      TPaveText *t1 = GetTitleText(Form("Efficiency for %1.0f < #DeltaTof < %1.2f ns",min[j],max[j]),0.06);
      t1->Draw();
      if(j==0)
	{
	  TLegend *leg = new TLegend(0.3,0.2,0.5,0.45);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->AddEntry(gCountDataEff[0][j],Form("%s: bin counting",name[0]),"p");
	  leg->AddEntry(gFitDataEff[0][j][0],Form("%s: fitting",name[0]),"p");
	  leg->AddEntry(gFitDataEff[1][j][1],Form("%s: fitting",name[1]),"p");
	  leg->Draw();
	}
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/AuAu200_DtofEff_FitVsCount.pdf",run_type));

}

//================================================
void kink(const Int_t save = 0)
{
  THnSparseF *hnKink = (THnSparseF*)f->Get(Form("mhMuonKink_%s",trigName[kTrigType]));
  hnKink->GetAxis(2)->SetRangeUser(1.6,100);

  // check invariant mass
  const char *name[3] = {"UL","LS","UL-LS"};
  TH1F *hInvMass[2];
  TH1F *hDr[2], *hDz[2];
  for(int i=0; i<2; i++)
    {
      hnKink->GetAxis(1)->SetRange(i+1,i+1);

      hnKink->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
      hDr[i] = (TH1F*)hnKink->Projection(3);
      hDr[i]->Sumw2();
      hDr[i]->SetName(Form("hDr_%s",name[i]));
      TH1F *htmp = (TH1F*)hnKink->Projection(6);
      htmp->Sumw2();
      htmp->SetName(Form("hDr_%s_2",name[i]));
      hDr[i]->Add(htmp);
      cout << hDr[i]->GetEntries() << endl;

      hDz[i] = (TH1F*)hnKink->Projection(4);
      hDz[i]->SetName(Form("hDz_%s",name[i]));
      hDz[i]->Sumw2();
      htmp = (TH1F*)hnKink->Projection(7);
      htmp->SetName(Form("hDz_%s_2",name[i]));
      htmp->Sumw2();
      hDz[i]->Add(htmp);
      hnKink->GetAxis(0)->SetRange(0,-1);

      hInvMass[i] = (TH1F*)hnKink->Projection(0);
      hInvMass[i]->SetName(Form("hInvMass_%s",name[i]));
      hInvMass[i]->Rebin(4);
      hInvMass[i]->SetMarkerStyle(21);
      hInvMass[i]->SetMarkerColor(color[1-i]);
      hInvMass[i]->SetLineColor(hInvMass[i]->GetMarkerColor());
    }
  hnKink->GetAxis(1)->SetRange(0,-1);
  TCanvas *c = draw1D(hInvMass[0],"Invariant mass distribution of dimuon pairs;M_{#mu#mu} (GeV/c^{2});Counts");
  hInvMass[1]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p_{T,1} > 1.5 GeV/c");
  leg->AddEntry(hInvMass[0],"Unlike-sign","P");
  leg->AddEntry(hInvMass[1],"Like-sign","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_InvMass.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_InvMass.pdf",run_type));
    }

  // analyze dr and dz distribution
  TH1F *hKinkDr[2], *hKinkDz[2];
  hKinkDr[0] = (TH1F*)hDr[0]->Clone("hKinkDr_US-LS");
  hKinkDr[0]->Add(hDr[1],-1);
  hKinkDz[0] = (TH1F*)hDz[0]->Clone("hKinkDz_US-LS");
  hKinkDz[0]->Add(hDz[1],-1);
  hnKink->GetAxis(1)->SetRange(1,1);
  hKinkDr[1] = (TH1F*)hnKink->Projection(3);
  hKinkDr[1]->Sumw2();
  hKinkDr[1]->SetName("hKinkDr_LS");
  TH1F *htmp = (TH1F*)hnKink->Projection(6);
  htmp->Sumw2();
  htmp->SetName("hKinkDr_LS_2");
  hKinkDr[1]->Add(htmp);
  hKinkDz[1] = (TH1F*)hnKink->Projection(4);
  hKinkDz[1]->Sumw2();
  hKinkDz[1]->SetName("hKinkDz_LS");
  htmp = (TH1F*)hnKink->Projection(7);
  htmp->Sumw2();
  htmp->SetName("hKinkDz_LS_2");
  hKinkDz[1]->Add(htmp);
  hnKink->GetAxis(1)->SetRange(0,-1);
  for(int i=0; i<2; i++)
    {
      hKinkDr[i]->Rebin(2);
      hKinkDr[i]->Scale(1./hKinkDr[i]->GetBinContent(hKinkDr[i]->FindFixBin(0)));
      hKinkDr[i]->SetMarkerStyle(21);
      hKinkDr[i]->SetMarkerColor(color[1-i]);
      hKinkDr[i]->SetLineColor(hKinkDr[i]->GetMarkerColor());

      hKinkDz[i]->Scale(1./hKinkDz[i]->GetBinContent(hKinkDz[i]->FindFixBin(0)));
      hKinkDz[i]->SetMarkerStyle(21);
      hKinkDz[i]->SetMarkerColor(color[1-i]);
      hKinkDz[i]->SetLineColor(hKinkDz[i]->GetMarkerColor());

      cout << "Fraction = " << hKinkDz[i]->Integral(hKinkDz[i]->FindFixBin(-0.15),hKinkDz[i]->FindFixBin(0.15))/hKinkDz[i]->Integral() << endl;
    }
  c = draw1D(hKinkDr[0],"Distance in transverse plane for muon candidates;#sqrt{(#Deltax)^{2}+(#Deltay)^{2}} (cm)");
  hKinkDr[1]->Draw("sames HIST");
  leg = new TLegend(0.6,0.65,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hKinkDr[0],"J/psi muon","P");
  leg->AddEntry(hKinkDr[1],"All muon","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dr.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dr.pdf",run_type));
    }

  hKinkDz[0]->GetXaxis()->SetRangeUser(-2,2);
  c = draw1D(hKinkDz[0],"Distance in z coordinate for muon candidates;#Deltaz (cm)");
  hKinkDz[1]->Draw("sames HIST");
  leg = new TLegend(0.6,0.65,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hKinkDz[0],"J/psi muon","P");
  leg->AddEntry(hKinkDz[1],"All muon","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dz.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dz.pdf",run_type));
    }
}



