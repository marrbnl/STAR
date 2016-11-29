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
  //getDtofCut();
  //kink();
}

//================================================
void MtdVpdTacDiff(const Int_t savePlot = 0, const int saveHisto = 0)
{
  TList *list = new TList;
  TH1F *htmp = 0x0;
  TCanvas *c = 0x0;

  const char *name = "AuAu200";
  const char *title = "Run14_AuAu200";
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum = 760;
  const double maximum = 840;
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};

  TFile *fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
  //TFile *fdata = TFile::Open("./output/Pico.Run14.AuAu200.JpsiMuon.dtof0.4.root","read");

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

  TH1F *hScaleFactor = new TH1F("hScaleFactor","Scale factor for unlike-sign background (online trigger efficiency);p_{T} (GeV/c)",nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      int start_bin = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin]-1e-4);
      double ul = hMuonPtSide[0]->Integral(start_bin,end_bin);
      double ls = hMuonPtSide[1]->Integral(start_bin,end_bin);
      double scale = ul/ls;
      double error = scale * sqrt(1/ul+1/ls);
      hScaleFactor->SetBinContent(bin, scale);
      hScaleFactor->SetBinError(bin, error);
    }
  hScaleFactor->SetMarkerStyle(21);
  hScaleFactor->GetYaxis()->SetRangeUser(0.95,1.15);
  c = draw1D(hScaleFactor);
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdVpdTacDiff_ScaleFactor.pdf",run_type));


  TH2F *hDataUL = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_UL_di_mu");
  hDataUL->Sumw2();
  TH2F *hDataLS = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_LS_di_mu");
  hDataLS->Sumw2();
  TH2F *hDataDisVsPt = (TH2F*)hDataUL->Clone(Form("JpsiMuonMtdVpdTacDiff"));
  hDataDisVsPt->Add(hDataLS, -1);
  const int nbins1 = 22;
  const double xbins1[nbins1+1] = {760,765,770,775,780,782,784,786,788,789,793,797,801,805,809,813,817,821,825,829,833,837,841};
  TH1F *hUL[nbins];
  TH1F *hLS[nbins];
  TH1F *hMuon[nbins];
  TH1F *hMuonFineBin[nbins];
  double mean_pt[nbins];
  double mean_pt_err[nbins];
  c = new TCanvas(Form("MtdVpdTacDiff_UL_vs_LS"),Form("MtdVpdTacDiff_UL_vs_LS"),1100,700);
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

      int start_bin = hDataUL->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hDataUL->GetXaxis()->FindBin(xbins[bin]-1e-4);

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",name,bin),start_bin,end_bin);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMarkerStyle(20);

      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",name,bin),start_bin,end_bin);
      hLS[bin-1]->SetMarkerStyle(24);
      hLS[bin-1]->SetMarkerColor(2);
      hLS[bin-1]->SetLineColor(2);

      double scale = hScaleFactor->GetBinContent(bin);
      hMuonFineBin[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("Data_JpsiMuon_MtdVpdTacDiff_bin%d",bin));
      hMuonFineBin[bin-1]->Add(hLS[bin-1],-1.0*scale);

      // rebin
      TH1F *htmp1 = (TH1F*)hUL[bin-1]->Clone(Form("%s_rebin_tmp",hUL[bin-1]->GetName()));
      hUL[bin-1] = (TH1F*)htmp1->Rebin(nbins1,Form("%s_Rebin",hUL[bin-1]->GetName()),xbins1);
      scaleHisto(hUL[bin-1], 1, 1, true);
      htmp1 = (TH1F*)hLS[bin-1]->Clone(Form("%s_rebin_tmp",hLS[bin-1]->GetName()));
      hLS[bin-1] = (TH1F*)htmp1->Rebin(nbins1,Form("%s_Rebin",hLS[bin-1]->GetName()),xbins1);
      scaleHisto(hLS[bin-1], 1, 1, true);
      hMuon[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("Data_JpsiMuon_MtdVpdTacDiff_bin%d_Rebin",bin));
      hMuon[bin-1]->Add(hLS[bin-1],-1.0*scale);

      c->cd(bin);
      hUL[bin-1]->SetMaximum(1.5*hUL[bin-1]->GetMaximum());
      hUL[bin-1]->GetXaxis()->SetRangeUser(minimum,maximum);
      hUL[bin-1]->SetTitle("");
      hUL[bin-1]->Draw("P");
      hLS[bin-1]->Draw("samesP");
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdVpdTacDiff_ULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("Data_JpsiMuon_MtdVpdTacDiff_FitMean"),Form("%s: mean of MtdVpdTacDiff;p_{T} (GeV/c)",name),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("Data_JpsiMuon_MtdVpdTacDiff_FitSigma"),Form("%s: sigma of MtdVpdTacDiff;p_{T} (GeV/c)",name),nbins,xbins);
  c = new TCanvas(Form("Fit_MtdVpdTacDiff"),Form("Fit_MtdVpdTacDiff"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));
      funcFitData[bin-1] = new TF1(Form("Data_JpsiMuon_MtdVpdTacDiffFit_bin%d",bin),"gaus",789,maximum);
      funcFitData[bin-1]->SetParameter(2,5);
      ptr[bin-1] = hFit->Fit(funcFitData[bin-1],"IR0QS");
      hFitDataMean->SetBinContent(bin,funcFitData[bin-1]->GetParameter(1));
      hFitDataMean->SetBinError(bin,funcFitData[bin-1]->GetParError(1));
      hFitDataSigma->SetBinContent(bin,funcFitData[bin-1]->GetParameter(2));
      hFitDataSigma->SetBinError(bin,funcFitData[bin-1]->GetParError(2));

      c->cd(bin);
      hFit->SetMaximum(1.5*hFit->GetMaximum());
      hFit->SetTitle("");
      hFit->GetXaxis()->SetRangeUser(minimum,maximum);
      hFit->Draw();
      funcFitData[bin-1]->SetLineColor(4);
      funcFitData[bin-1]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdVpdTacDiffFit.pdf",run_type));

  hFitDataMean->SetMarkerStyle(21);
  hFitDataMean->GetYaxis()->SetRangeUser(780,810);
  c = draw1D(hFitDataMean,"Mean of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);Mean");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdVpdTacDiff_FitMean.pdf",run_type));

  hFitDataSigma->SetMarkerStyle(21);
  hFitDataSigma->GetYaxis()->SetRangeUser(0,14);
  c = draw1D(hFitDataSigma,"Sigma of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);#sigma");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdVpdTacDiff_FitSigma.pdf",run_type));

  //==============================================
  // calculate efficiency
  //==============================================
  const int nLumi = 2;
  const char *name_lumi[nLumi] = {"prod_high","prod_low"};
  const double min[nLumi] = {789,786};
  const double max[nLumi] = {837,837};
  TF1 *funcBay[nLumi][nbins];
  TCanvas *cFuncBay[nLumi];
  double x,y;
  TGraphAsymmErrors *gDataEff[nLumi][2];
  TH1F *hBase[nLumi][2];
  TH1F *hMatch[nLumi][2];
  for(int k=0; k<nLumi; k++)
    {
      for(int i=0; i<2; i++)
	{
	  hBase[k][i]  = new TH1F(Form("hBase_%d_M%d",k,i),"",nbins,xbins);
	  hMatch[k][i] = new TH1F(Form("hMatch_%d_M%d",k,i),"",nbins,xbins);
	}

      // fitting method
      for(int bin=1; bin<=nbins; bin++)
	{
	  double all = funcFitData[bin-1]->Integral(minimum,maximum);
	  double all_err = funcFitData[bin-1]->IntegralError(minimum,maximum,funcFitData[bin-1]->GetParameters(),ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  hBase[k][0]->SetBinContent(bin,all);
	  hBase[k][0]->SetBinError(bin,all_err);
	  double acc = funcFitData[bin-1]->Integral(min[k],max[k]);
	  double acc_err = funcFitData[bin-1]->IntegralError(min[k],max[k],funcFitData[bin-1]->GetParameters(),ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  hMatch[k][0]->SetBinContent(bin,acc);
	  hMatch[k][0]->SetBinError(bin,acc_err);
	}
      gDataEff[k][0] = new TGraphAsymmErrors(hMatch[k][0], hBase[k][0],"cl=0.683 b(1,1) mode");
      gDataEff[k][0]->SetName(Form("Data_JpsiMuon_MtdTrigEff_Fitting_%s",name_lumi[k]));
      for(int ipoint=0; ipoint<nbins; ipoint++)
	{
	  gDataEff[k][0]->GetPoint(ipoint,x,y);
	  gDataEff[k][0]->SetPoint(ipoint,mean_pt[ipoint],y);
	  gDataEff[k][0]->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
	  gDataEff[k][0]->SetPointEXlow(ipoint,mean_pt_err[ipoint]);
	}

      // Bin counting method
      for(int bin=1; bin<=nbins; bin++)
	{
	  double mean = funcFitData[bin-1]->GetParameter(1);
	  double value = 0, error = 0;
	  int start_bin = hMuonFineBin[bin-1]->FindFixBin(minimum+1e-4);
	  int end_bin   = hMuonFineBin[bin-1]->FindFixBin(maximum-1e-4);
	  int new_bin = -1;
	  for(int ibin=start_bin; ibin<=end_bin; ibin++)
	    {
	      //cout << hMuonFineBin[bin-1] ->GetBinWidth(ibin) << endl;
	      if(hMuonFineBin[bin-1]->GetXaxis()->GetBinLowEdge(ibin)<min[0])
		new_bin = hMuonFineBin[bin-1]->FindFixBin(2*mean-hMuonFineBin[bin-1]->GetBinCenter(ibin));
	      else
		new_bin = ibin;
	      //printf("[i] new bin = %d, old bin = %d, center = %f\n",new_bin,ibin,hMuonFineBin[bin-1]->GetBinCenter(ibin));

	      value += hMuonFineBin[bin-1]->GetBinContent(new_bin);
	      error += TMath::Power(hMuonFineBin[bin-1]->GetBinError(new_bin),2);
	    }
	  error = TMath::Sqrt(error);
	  hBase[k][1]->SetBinContent(bin,value);
	  hBase[k][1]->SetBinError(bin,error);
	  //printf("[i] all: %f\n",value);

	  // accepted
	  start_bin = hMuonFineBin[bin-1]->FindFixBin(min[k]+1e-4);
	  end_bin   = hMuonFineBin[bin-1]->FindFixBin(max[k]-1e-4);
	  value = 0, error = 0;
	  for(int ibin=start_bin; ibin<=end_bin; ibin++)
	    {
	      if(hMuonFineBin[bin-1]->GetXaxis()->GetBinLowEdge(ibin)<min[0])
		new_bin = hMuonFineBin[bin-1]->FindFixBin(2*mean-hMuonFineBin[bin-1]->GetBinCenter(ibin));
	      else
		new_bin = ibin;
	      value += hMuonFineBin[bin-1]->GetBinContent(new_bin);
	      error += TMath::Power(hMuonFineBin[bin-1]->GetBinError(new_bin),2);
	    }
	  error = TMath::Sqrt(error);
	  hMatch[k][1]->SetBinContent(bin,value);
	  hMatch[k][1]->SetBinError(bin,error);
	  //printf("[i] acc: %f\n",value);
	}
      gDataEff[k][1] = new TGraphAsymmErrors(hMatch[k][1], hBase[k][1],"cl=0.683 b(1,1) mode");
      gDataEff[k][1]->SetName(Form("Data_JpsiMuon_MtdTrigEff_BinCount_%s",name_lumi[k]));
      for(int ipoint=0; ipoint<nbins; ipoint++)
	{
	  gDataEff[k][1]->GetPoint(ipoint,x,y);
	  gDataEff[k][1]->SetPoint(ipoint,mean_pt[ipoint],y);
	  gDataEff[k][1]->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
	  gDataEff[k][1]->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
	}

      // probability distribution
      cFuncBay[k] = new TCanvas(Form("cFuncBay_%s",name_lumi[k]),Form("cFuncBay_%s",name_lumi[k]),800,600);
      TH1F *htmp = new TH1F(Form("htmp_%s",name_lumi[k]),";#epsilon;P(#epsilon)",1000,0.5,1.0);
      htmp->GetYaxis()->SetRangeUser(0,0.5);
      htmp->Draw();
      TLegend *leg = new TLegend(0.15,0.5,0.5,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("%s: %s",name,name_lumi[k]));
      for(int ipoint=0; ipoint<nbins; ipoint++)
	{	  
	  gDataEff[k][1]->GetPoint(ipoint,x,y);
	  double pw  = hMatch[k][1]->GetBinContent(ipoint+1);
	  double pw2 = hMatch[k][1]->GetSumw2()->At(ipoint+1);
	  double tw  = hBase[k][1]->GetBinContent(ipoint+1);
	  double tw2 = hBase[k][1]->GetSumw2()->At(ipoint+1);
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
	  funcBay[k][ipoint] = new TF1(Form("funcBay_%s_%d",name_lumi[k],ipoint),"[0]*x**[1]*(1-x)**([2]-[1])",0,1);
	  double binwidth = (funcBay[k][ipoint]->GetXmax()-funcBay[k][ipoint]->GetXmin())/funcBay[k][ipoint]->GetNpx();
	  funcBay[k][ipoint]->SetParameters(gamma*binwidth,p,N);	

	  cFuncBay[k]->cd();
	  funcBay[k][ipoint]->SetLineColor(color[ipoint]);
	  funcBay[k][ipoint]->Draw("same");
	  leg->AddEntry(funcBay[k][ipoint],Form("p_{T} = %1.2f GeV/c",x),"l");
	}
      leg->Draw();
      if(savePlot)
	cFuncBay[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMUuon_MtdTrigEff_ProbDis_%s.pdf",run_type,name_lumi[k]));
    }
  
  // trigger efficiency
  TFile *fTrigEff = TFile::Open("Rootfiles/Run14.AuAu200.MuonTrigEff.root","read");
  TCanvas *cTrigEff[nLumi];
  TF1 *funcEff[nLumi];
  TLegend *legdata[nLumi];
  TH1F *gFit_err[nLumi][2];
  TF1 *fExtrap[nLumi];
  for(int k=0; k<nLumi; k++)
    {
      gDataEff[k][0]->SetMarkerStyle(24);
      gDataEff[k][0]->SetMarkerSize(1.2);
      gDataEff[k][0]->GetYaxis()->SetRangeUser(0.5,1);
      gDataEff[k][0]->GetXaxis()->SetRangeUser(1.2,8);
      cTrigEff[k] = drawGraph(gDataEff[k][0],"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);

      gDataEff[k][1]->SetMarkerStyle(21);
      gDataEff[k][1]->SetMarkerSize(1.2);
      gDataEff[k][1]->GetYaxis()->SetRangeUser(0.5,1);
      gDataEff[k][1]->GetXaxis()->SetRangeUser(1.2,8);
      gDataEff[k][1]->Draw("samesPE");

      TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gDataEff[k][1]->Clone(Form("Fit_%s",gDataEff[k][1]->GetName()));
      funcEff[k] = new TF1(Form("%s_FitFunc",gDataEff[k][1]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      gfit->Fit(funcEff[k],"R0");
      funcEff[k]->SetLineColor(2);
      funcEff[k]->Draw("sames");

      fExtrap[k] = (TF1*)fTrigEff->Get(Form("MuonTrigEff_cent0060_P%d",3-k*2));
      fExtrap[k]->SetLineStyle(5);
      fExtrap[k]->SetLineColor(4);
      fExtrap[k]->DrawCopy("sames");

      legdata[k] = new TLegend(0.4,0.2,0.6,0.45);
      legdata[k]->SetBorderSize(0);
      legdata[k]->SetFillColor(0);
      legdata[k]->SetTextSize(0.04);
      legdata[k]->SetHeader(Form("Run14_AuAu200, %s",name_lumi[k]));
      legdata[k]->AddEntry(gDataEff[k][0],"J/#Psi muons: fit method","P");
      legdata[k]->AddEntry(gDataEff[k][1],"J/#Psi muons: bin counting","P");
      legdata[k]->AddEntry(funcEff[k],"Fit to bin counting","L");
      legdata[k]->AddEntry(fExtrap[k],"Extrapolation method","L");
      legdata[k]->Draw();
      if(savePlot) cTrigEff[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdTrigEff_%s.pdf",run_type,name_lumi[k]));
    }

  //==============================================
  // systematic uncertainty
  //==============================================
  //Method 1 - randomize efficiency
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int defMed = 1;
  const int nexpr = 400;
  TGraphAsymmErrors *gEffRndm[nLumi][nexpr];
  for(int k=0; k<nLumi; k++)
    {
      for(int j = 0; j < nexpr; j++)
	{
	  TGraphAsymmErrors *gBaseEff = gDataEff[k][defMed]; // bin counting method 
	  int npoint = gBaseEff->GetN();
	  gEffRndm[k][j] = new TGraphAsymmErrors(npoint);
	  gEffRndm[k][j]->SetName(Form("%s_Rndm%d",gBaseEff->GetName(),j));
	  double *exh = gBaseEff->GetEXhigh();
	  double *exl = gBaseEff->GetEXlow();
	  double *eyh = gBaseEff->GetEYhigh();
	  double *eyl = gBaseEff->GetEYlow();
	  for(int ipoint=0; ipoint<npoint; ipoint++)
	    {
	      gBaseEff->GetPoint(ipoint,x,y);
	      y = funcBay[k][ipoint]->GetRandom(0,1);
	      gEffRndm[k][j]->SetPoint(ipoint,x,y);
	      gEffRndm[k][j]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	    }
	}
    }

  TH1F *hplot = new TH1F("hplot",";p_{T} (GeV/c);Efficiency",100,0,10);
  hplot->GetYaxis()->SetRangeUser(0.6,1);
  hplot->GetXaxis()->SetRangeUser(1.3,7);
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1.4,1.5,1.75,2,2.25,2.5,3,4,5,6,8};
  TH1F *hSpread[nLumi][nbins_rndm];
  TGraphErrors *gLimits[nLumi][2];
  TF1 *funcLimits[nLumi][2];
  const char *limit_name[2] = {"low","up"};
  for(int k=0; k<nLumi; k++)
    {
      c = new TCanvas(Form("cEff_%s",name_lumi[k]),Form("cEff_%s",name_lumi[k]),800,600);
      hplot->DrawCopy();
      TPaveText *t1 = GetTitleText("Estimation of trigger efficiency uncertainty");
      t1->Draw();
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("%s, %s",title,name_lumi[k]));

      TCanvas *cFit = new TCanvas(Form("cEffFit_%s",name_lumi[k]),Form("cEffFit_%s",name_lumi[k]),1200,800);
      cFit->Divide(4,3);
      for(int l=0; l<nbins_rndm; l++)
	{
	  hSpread[k][l] = new TH1F(Form("hSpread_%s_%d",name_lumi[k],l),Form("hSpread_%s_%d",name_lumi[k],l),600,0,1.2);
	}

      for(int j = 0; j < nexpr; j++)
	{
	  TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[k][j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
	  funcTmp->SetParLimits(0,0,1);
	  funcTmp->SetParLimits(1,0,5);
	  funcTmp->SetParLimits(2,0,1);
	  gEffRndm[k][j]->Fit(funcTmp,"R0Q");
	  c->cd();
	  funcTmp->SetLineStyle(2);
	  funcTmp->Draw("sames");
	  if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
	  for(int l=0; l<nbins_rndm; l++)
	    {
	      hSpread[k][l]->Fill(funcTmp->Eval(xbins_rndm[l]));
	    }
	}

      gLimits[k][0] = new TGraphErrors(nbins_rndm);
      gLimits[k][1] = new TGraphErrors(nbins_rndm);
      for(int l=0; l<nbins_rndm; l++)
	{
	  cFit->cd(l+1);
	  if(k==0)
	    {
	      hSpread[k][l]->Rebin(2);
	      if(l<2) hSpread[k][l]->GetXaxis()->SetRangeUser(0.3,1);
	      else if(l<4) hSpread[k][l]->GetXaxis()->SetRangeUser(0.65,1);
	      else hSpread[k][l]->GetXaxis()->SetRangeUser(0.75,1);
	    }
	  else
	    {
	      if(l<3) hSpread[k][l]->GetXaxis()->SetRangeUser(0.5,1);
	      else hSpread[k][l]->GetXaxis()->SetRangeUser(0.88,1);
	    }
	  hSpread[k][l]->SetTitle(";Efficiency;");
	  hSpread[k][l]->SetMaximum(1.3*hSpread[k][l]->GetMaximum());
	  hSpread[k][l]->Draw();
	  TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[l]),0.06);
	  t1->Draw();
	  double boundary = hSpread[k][l]->GetBinCenter(hSpread[k][l]->GetMaximumBin());
	  cout << boundary << endl;

	  TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[k][l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
	  funcTmp2->SetParameter(0,hSpread[k][l]->GetMaximum());
	  funcTmp2->SetParameter(1,hSpread[k][l]->GetMean());
	  funcTmp2->SetParameter(2,hSpread[k][l]->GetRMS());
	  hSpread[k][l]->Fit(funcTmp2,"R0Q");
	  funcTmp2->SetLineColor(4);
	  funcTmp2->SetLineStyle(2);
	  funcTmp2->Draw("sames");

	  double pt = xbins_rndm[l];
	  //double mean = funcTmp2->GetParameter(1);
	  double mean = funcEff[k]->Eval(pt);
	  //double width = TMath::Abs(funcTmp2->GetParameter(2));
	  //double error = funcTmp2->GetParError(2);
	  //double mean  = hSpread[k][l]->GetMean();
	  double nSigma = 1;
	  if(l<3) nSigma = 3;
	  double width = hSpread[k][l]->GetRMS() * nSigma;
	  double error = hSpread[k][l]->GetRMSError() * nSigma;
	  // if(l<2) 
	  //   {
	  //     mean  = hSpread[k][l]->GetMean();
	  //   }
	  gLimits[k][1]->SetPoint(l,pt,mean+width);
	  gLimits[k][1]->SetPointError(l,0,error);
	  gLimits[k][0]->SetPoint(l,pt,mean-width);
	  gLimits[k][0]->SetPointError(l,0,error);

	  TPaveText *t1 = GetPaveText(0.6,0.7,0.7,0.85,0.06);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("RMS=%4.3f",hSpread[k][l]->GetRMS()));
	  t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcTmp2->GetParameter(2))));
	  t1->Draw();
	}
      cFit->cd(12);
      TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
      t1->AddText(name);
      t1->AddText(name_lumi[k]);
      t1->Draw();
      if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_MuonJpsi_MtdTrigEff_RndmFit_%s.pdf",run_type,name_lumi[k]));

      c->cd();
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_MuonJpsi_MtdTrigEff_RndmCurve_%s.pdf",run_type,name_lumi[k]));
      for(int t=0; t<2; t++)
	{
	  gLimits[k][t]->SetMarkerStyle(20);
	  gLimits[k][t]->SetMarkerColor(2+2*t);
	  gLimits[k][t]->SetLineColor(2+2*t);
	  gLimits[k][t]->Draw("samesP");

	  funcLimits[k][t] = new TF1(Form("%s_Sys1_%s",gDataEff[k][defMed]->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  funcLimits[k][t]->SetParameters(0.95,2,0.2);
	  gLimits[k][t]->Fit(funcLimits[k][t],"R0");
	  funcLimits[k][t]->SetLineColor(gLimits[k][t]->GetMarkerColor());
	  funcLimits[k][t]->DrawCopy("sames");
	}
      leg->AddEntry(gLimits[k][0],"Lower limit","PL");
      leg->AddEntry(gLimits[k][1],"Upper limit","PL");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_MuonJpsi_MtdTrigEff_RndmLimits_%s.pdf",run_type,name_lumi[k]));
    }
  return;

  TCanvas *cSys[2];
  for(int k=0; k<2; k++)
    {
      TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gDataEff[k][defMed]->Clone(Form("%s_sys",gDataEff[k][defMed]->GetName()));
      gtmp->GetXaxis()->SetRangeUser(1.3,8);
      gtmp->GetYaxis()->SetRangeUser(0.6,1);
      cSys[k] = drawGraph(gtmp,"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);
      funcEff[k]->Draw("sames");
      gDataEff[k][0]->Draw("samesPE");
      for(int t=0; t<2; t++)
	{
	  funcLimits[k][t]->SetLineColor(4);
	  funcLimits[k][t]->SetLineStyle(4);
	  funcLimits[k][t]->Draw("sames");
	  cout << k << "  " << funcLimits[k][t]->Eval(1.4) << endl;
	}
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.5);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("Run14_AuAu200, %s",name_lumi[k]));
      leg->AddEntry(gDataEff[k][0],"J/#Psi muons: fit method","P");
      leg->AddEntry(gDataEff[k][1],"J/#Psi muons: bin counting","P");
      leg->AddEntry(funcEff[k],"Fit to bin counting","L");
      leg->AddEntry(funcLimits[k][0],"Uncertainty","L");
      leg->Draw();
      if(savePlot) cSys[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdTrigEffSys_%s.pdf",run_type,name_lumi[k]));
    }  
return;
  //Method 2 - assume correlated statistical errors
  TCanvas *cSys2[2];
  TF1 *funcLimits2[nLumi][2];
  for(int k=0; k<2; k++)
    {
      TGraphAsymmErrors *gtmp = (TGraphAsymmErrors*)gDataEff[k][defMed]->Clone(Form("%s_sys2",gDataEff[k][defMed]->GetName()));
      gtmp->GetXaxis()->SetRangeUser(1.3,8);
      gtmp->GetYaxis()->SetRangeUser(0.6,1);
      cSys2[k] = drawGraph(gtmp,"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      gPad->SetGrid(0,1);
      funcEff[k]->Draw("sames");
      gDataEff[k][0]->Draw("samesPE");

      TGraphAsymmErrors *gBaseEff = gDataEff[k][defMed];
      int npoint = gBaseEff->GetN();
      for(int t=0; t<2; t++)
	{
	  TGraphAsymmErrors *gSys2 = new TGraphAsymmErrors(npoint);
	  gSys2->SetName(Form("Fit2_%s_%s_%d",name,name_lumi[k],t));
	  double *exh = gBaseEff->GetEXhigh();
	  double *exl = gBaseEff->GetEXlow();
	  double *eyh = gBaseEff->GetEYhigh();
	  double *eyl = gBaseEff->GetEYlow();
	  double nsigma = 1;
	  for(int ipoint=0; ipoint<nbins; ipoint++)
	    {
	      gBaseEff->GetPoint(ipoint,x,y);
	      if(ipoint==0) nsigma = 1.5;
	      else          nsigma = 1;
	      if(t==0) y = y + nsigma*eyh[ipoint];
	      if(t==1) y = y - nsigma*eyl[ipoint];
	      gSys2->SetPoint(ipoint,x,y);
	      gSys2->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	    }
	  funcLimits2[k][t] = new TF1(Form("%s_Sys2_%s",gDataEff[k][defMed]->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  funcLimits2[k][t]->SetParameters(0.95,2,0.2);
	  gSys2->Fit(funcLimits2[k][t],"R0Q");
	  funcLimits2[k][t]->SetLineColor(kMagenta-4);
	  funcLimits2[k][t]->SetLineStyle(2);
	  funcLimits2[k][t]->DrawCopy("sames");
	}
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.5);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("Run14_AuAu200, %s",name_lumi[k]));
      leg->AddEntry(gDataEff[k][0],"J/#Psi muons: fit method","P");
      leg->AddEntry(gDataEff[k][1],"J/#Psi muons: bin counting","P");
      leg->AddEntry(funcEff[k],"Fit to bin counting","L");
      leg->AddEntry(funcLimits2[k][0],"Uncertainty","L");
      leg->Draw();
      if(savePlot) cSys2[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_MtdTrigEffSys2_%s.pdf",run_type,name_lumi[k]));
    }
 
  
  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type),"update");
      for(int bin=1; bin<=nbins; bin++)
	{
	  hMuonFineBin[bin-1]->Write("",TObject::kOverwrite);
	  hMuon[bin-1]->Write("",TObject::kOverwrite);
	}
      hFitDataMean->Write("",TObject::kOverwrite);
      hFitDataSigma->Write("",TObject::kOverwrite);

      for(int k=0; k<2; k++)
	{
	  gDataEff[k][0]->Write("",TObject::kOverwrite);
	  gDataEff[k][1]->Write("",TObject::kOverwrite);
	  funcEff[k]->Write("",TObject::kOverwrite);
	  for(int t=0; t<2; t++)
	    {
	      funcLimits[k][t]->Write("",TObject::kOverwrite);
	      funcLimits2[k][t]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void DeltaTof(const Int_t savePlot = 0, const int saveHisto = 0)
{
  const char *name = "AuAu200";
  const char *title = "Run14_AuAu_200";
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum = -3;
  const double maximum = 3;
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};

  TFile *fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");

  //==============================================
  // compare invariant mass
  //==============================================
  const int nbinsJpsi = 5;
  const double xbinsJpsi[nbinsJpsi+1] = {0,1.0,2.0,3.0,5.0,10.0};
  TH2F *hInvMassVsPtUL = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_Dtof_UL_di_mu"));
  TH2F *hInvMassVsPtLS = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_Dtof_LS_di_mu"));

  TCanvas *c = new TCanvas("InvMass_Dtof","InvMass_Dtof",1100,700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetHeader(run_type);
  for(int bin=1; bin<=nbinsJpsi; bin++)
    {
      int start_bin = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin-1]+1e-4);
      int end_bin   = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin]-1e-4);
      TH1F *hInvMassUL = (TH1F*)hInvMassVsPtUL->ProjectionY(Form("Data_Dtof_InvMassUL_bin%d",bin),start_bin,end_bin);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->Rebin(5);
      hInvMassUL->SetMaximum(1.5*hInvMassUL->GetMaximum());
      hInvMassUL->SetMinimum(0);
      hInvMassUL->GetXaxis()->SetRangeUser(2.2,4);
      hInvMassUL->SetTitle("");
      if(bin==1) leg->AddEntry(hInvMassUL,"Unlike-sign","PL");

      TH1F *hInvMassLS = (TH1F*)hInvMassVsPtLS->ProjectionY(Form("Data_Dtof_InvMassLS_bin%d",bin),start_bin,end_bin);
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
	  TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMassUL->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	    {
	      int jbin = hInvMassUL->FindFixBin(htmp->GetBinCenter(ibin));
	      htmp->SetBinContent(ibin,hInvMassUL->GetBinContent(jbin));
	    }
	  if(itmp==0) htmp->SetFillColor(7);
	  else htmp->SetFillColor(5);
	  htmp->Draw("sames");
	  if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	  if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	}
      hInvMassLS->Draw("samesP");
      TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbinsJpsi[bin-1],xbinsJpsi[bin]),0.05);
      t1->Draw();
      int low_bin  = hInvMassUL->FindFixBin(min_mass[0]+1e-4);
      int high_bin = hInvMassUL->FindFixBin(max_mass[0]-1e-4);
      double all = hInvMassUL->Integral(low_bin, high_bin);
      double bkg = hInvMassLS->Integral(low_bin, high_bin);
      TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
      t1->SetTextAlign(11);
      t1->AddText(Form("S/B = 1:%4.2f",bkg/(all-bkg)));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_Dtof_InvMassLSvsUL.pdf",run_type));

  //==============================================
  // single muon UL vs LS
  //==============================================
  TH2F *hJpsiMassVsMuonPt[2];
  TH1F *hMuonPtSide[2];
  hJpsiMassVsMuonPt[0] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_Dtof_UL_di_mu"));
  hJpsiMassVsMuonPt[1] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_Dtof_LS_di_mu"));
  for(int i=0; i<2; i++)
    {
      // Get the muon distribution in the side-band region
      hJpsiMassVsMuonPt[i]->Sumw2();
      hMuonPtSide[i] = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_Dtof_SideBand_%d",i));
      hMuonPtSide[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  htmp = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_MtdVpdTacDiff_SideBand%d_%d",j,i),low_bin,up_bin);
	  hMuonPtSide[i]->Add(htmp);
	}
    }

  TH1F *hScaleFactor = new TH1F("hScaleFactor","Scale factor for unlike-sign background (#Deltatof cut);p_{T} (GeV/c)",nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      int start_bin = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin]-1e-4);
      double ul = hMuonPtSide[0]->Integral(start_bin,end_bin);
      double ls = hMuonPtSide[1]->Integral(start_bin,end_bin);
      double scale = ul/ls;
      double error = scale * sqrt(1/ul+1/ls);
      hScaleFactor->SetBinContent(bin, scale);
      hScaleFactor->SetBinError(bin, error);
    }
  hScaleFactor->SetMarkerStyle(21);
  hScaleFactor->GetYaxis()->SetRangeUser(0.95,1.15);
  c = draw1D(hScaleFactor);
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_Dtof_ScaleFactor.pdf",run_type));


  TH2F *hDataUL = (TH2F*)fdata->Get("mhJpsiMuon_Dtof_UL_di_mu");
  hDataUL->Sumw2();
  TH2F *hDataLS = (TH2F*)fdata->Get("mhJpsiMuon_Dtof_LS_di_mu");
  hDataLS->Sumw2();
  TH2F *hDataDisVsPt = (TH2F*)hDataUL->Clone("JpsiMuonDtof");
  hDataDisVsPt->Add(hDataLS, -1);
  TH1F *hUL[nbins];
  TH1F *hLS[nbins];
  TH1F *hMuon[nbins];
  TH1F *hMuonFineBin[nbins];
  double mean_pt[nbins];
  double mean_pt_err[nbins];
  const int rebin = 10;
  c = new TCanvas(Form("MtdVpdTacDiff_UL_vs_LS"),Form("MtdVpdTacDiff_UL_vs_LS"),1100,700);
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

      int start_bin = hDataUL->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hDataUL->GetXaxis()->FindBin(xbins[bin]-1e-4);

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",name,bin),start_bin,end_bin);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMarkerStyle(20);

      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",name,bin),start_bin,end_bin);
      hLS[bin-1]->SetMarkerStyle(24);
      hLS[bin-1]->SetMarkerColor(2);
      hLS[bin-1]->SetLineColor(2);

      double scale = hScaleFactor->GetBinContent(bin);
      hMuonFineBin[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("Data_JpsiMuon_Dtof_bin%d",bin));
      hMuonFineBin[bin-1]->Add(hLS[bin-1],-1.0*scale);

      // rebin
      hMuon[bin-1] = (TH1F*)hMuonFineBin[bin-1]->Clone(Form("Data_JpsiMuon_Dtof_bin%d_Rebin",bin));
      hMuon[bin-1]->Rebin(rebin);
      
      c->cd(bin);
      hUL[bin-1]->Rebin(rebin);
      hUL[bin-1]->SetMaximum(1.5*hUL[bin-1]->GetMaximum());
      hUL[bin-1]->GetXaxis()->SetRangeUser(minimum,maximum);
      hUL[bin-1]->SetTitle("");
      hUL[bin-1]->Draw("P");
      hLS[bin-1]->Rebin(rebin);
      hLS[bin-1]->Draw("samesP");
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_Dtof_ULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("Data_JpsiMuon_Dtof_FitMean"),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",name),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("Data_JpsiMuon_Dtof_FitSigma"),Form("%s: sigma of #Deltatof;p_{T} (GeV/c)",name),nbins,xbins);
  c = new TCanvas(Form("Fit_Dtof"),Form("Fit_Dtof"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));

      //funcFitData[bin-1] = new TF1(Form("Data_JpsiMuon_Dtof_bin%d_Fit",bin),CrystalBall,-3,3,5);
      //funcFitData[bin-1]->SetParameters(-1.5, -0.05, 0.2, 50, hFit->GetMaximum());
      funcFitData[bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,3);
      ptr[bin-1] = hFit->Fit(funcFitData[bin-1],"IR0QS");

      hFitDataMean->SetBinContent(bin,funcFitData[bin-1]->GetParameter(1));
      hFitDataMean->SetBinError(bin,funcFitData[bin-1]->GetParError(1));
      hFitDataSigma->SetBinContent(bin,funcFitData[bin-1]->GetParameter(2));
      hFitDataSigma->SetBinError(bin,funcFitData[bin-1]->GetParError(2));

      c->cd(bin);
      hFit->SetMaximum(1.5*hFit->GetMaximum());
      hFit->SetTitle("");
      hFit->GetXaxis()->SetRangeUser(minimum,maximum);
      hFit->Draw();
      funcFitData[bin-1]->SetLineColor(4);
      funcFitData[bin-1]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_DtofFit.pdf",run_type));

  hFitDataMean->SetMarkerStyle(21);
  hFitDataMean->GetYaxis()->SetRangeUser(-1,1);
  c = draw1D(hFitDataMean,"Mean of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);Mean");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_Dtof_FitMean.pdf",run_type));

  hFitDataSigma->SetMarkerStyle(21);
  hFitDataSigma->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hFitDataSigma,"Width of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);#sigma");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_Dtof_FitSigma.pdf",run_type));

  //==============================================
  // extract efficiency
  //==============================================
  const double dtof_cut_min = -3;
  const double dtof_cut_max = 0.46;

  // bin counting method for data
  double x,y;
  TH1F *hBase = new TH1F("hBase","hBase",nbins,xbins);
  TH1F *hMatch = new TH1F("hMatch","hMatch",nbins,xbins);
  TH1F *hMatch2 = new TH1F("hMatch2","hMatch2",nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      double all_err;
      int low_bin = hMuonFineBin[bin-1]->FindFixBin(minimum+1e-4);
      int high_bin = hMuonFineBin[bin-1]->FindFixBin(maximum-1e-4);
      double all = hMuonFineBin[bin-1]->IntegralAndError(low_bin, high_bin,all_err);

      double acc_err;
      low_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_min+1e-4);
      high_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_max-1e-4);
      double acc = hMuonFineBin[bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
      hBase->SetBinContent(bin,all);
      hBase->SetBinError(bin,all_err);
      hMatch->SetBinContent(bin,acc);
      hMatch->SetBinError(bin,acc_err);
      if(acc>all) { acc = all; acc_err = all_err; }
      hMatch2->SetBinContent(bin,acc);
      hMatch2->SetBinError(bin,acc_err);
    }
  
  TGraphAsymmErrors *gCountDataEff = new TGraphAsymmErrors(hMatch2, hBase,"cl=0.683 b(1,1) mode");
  gCountDataEff->SetName(Form("Data_JpsiMuon_DtofEff%1.0f_BinCount",dtof_cut_max*100));
  hMatch->Divide(hBase);
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint,x,y);
      gCountDataEff->SetPoint(ipoint,mean_pt[ipoint],y);
      gCountDataEff->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
      gCountDataEff->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
      cout << hMatch->GetBinContent(ipoint+1) << endl;
      if(hMatch->GetBinContent(ipoint+1)>1)
	{
	  gCountDataEff->SetPoint(ipoint,mean_pt[ipoint],hMatch->GetBinContent(ipoint+1));
	  gCountDataEff->SetPointEYhigh(ipoint, hMatch->GetBinError(ipoint+1));
	  gCountDataEff->SetPointEYlow(ipoint, hMatch->GetBinError(ipoint+1));
	}
    }

  // fitting method
  hBase->Reset();
  hMatch->Reset();
  for(int bin=1; bin<=nbins; bin++)
    {
      double all = funcFitData[bin-1]->Integral(minimum,maximum);
      double all_err = funcFitData[bin-1]->IntegralError(minimum,maximum,funcFitData[bin-1]->GetParameters(), ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
      double acc = funcFitData[bin-1]->Integral(dtof_cut_min, dtof_cut_max);
      double acc_err = funcFitData[bin-1]->IntegralError(dtof_cut_min, dtof_cut_max,funcFitData[bin-1]->GetParameters(), ptr[i]->GetCovarianceMatrix().GetMatrixArray());
      hBase->SetBinContent(bin,all);
      hBase->SetBinError(bin,all_err);
      hMatch->SetBinContent(bin,acc);
      hMatch->SetBinError(bin,acc_err);
    }
  TGraphAsymmErrors *gFitDataEff = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
  gFitDataEff->SetName(Form("Data_JpsiMuon_DtofEff%1.0f_Fitting",dtof_cut_max*100));
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gFitDataEff->GetPoint(ipoint,x,y);
      gFitDataEff->SetPoint(ipoint,mean_pt[ipoint],y);
      gFitDataEff->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
      gFitDataEff->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
    }

   // trigger efficiency
  TH1F *hplot = new TH1F("hplot","",100,0,10);
  hplot->GetYaxis()->SetRangeUser(0.6,1.05);
  hplot->GetXaxis()->SetRangeUser(1.2,7);
  hplot->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
  c = new TCanvas("DtofEff_Comp","DtofEff_Comp",800,600);
  gPad->SetGridy();
  hplot->DrawCopy();
  gCountDataEff->SetMarkerStyle(21);
  gCountDataEff->SetMarkerSize(1.3);
  gCountDataEff->Draw("PZsames");
  gFitDataEff->SetMarkerSize(1.3);
  gFitDataEff->SetMarkerStyle(24);
  gFitDataEff->SetMarkerColor(2);
  gFitDataEff->SetLineColor(2);
  gFitDataEff->Draw("PZsames");
  TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("Fit_%s",gCountDataEff->GetName()));
  TF1 *funcEff = new TF1(Form("%s_FitFunc",gCountDataEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  gfit->Fit(funcEff,"R0");
  funcEff->SetLineColor(2);
  funcEff->Draw("sames");

  TPaveText *t1 = GetTitleText(Form("Efficiency for %1.0f < #DeltaTof < %1.2f ns",dtof_cut_min,dtof_cut_max),0.04);
  t1->Draw();
  TLegend *leg = new TLegend(0.3,0.2,0.5,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  leg->AddEntry(gCountDataEff,Form("Data: bin counting"),"p");
  leg->AddEntry(gFitDataEff,Form("Data: fitting"),"p");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Data_JpsiMuon_DtofEff_FitVsCount.pdf",run_type));

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type),"update");
      for(int bin=1; bin<=nbins; bin++)
	{
	  hMuonFineBin[bin-1]->Write("",TObject::kOverwrite);
	  hMuon[bin-1]->Write("",TObject::kOverwrite);
	  funcFitData[bin-1]->Write("",TObject::kOverwrite);
	  ptr[bin-1]->Write("",TObject::kOverwrite);
	}
      hFitDataMean->Write("",TObject::kOverwrite);
      hFitDataSigma->Write("",TObject::kOverwrite);
      gFitDataEff->Write("",TObject::kOverwrite);
      gCountDataEff->Write("",TObject::kOverwrite);
      funcEff->Write("",TObject::kOverwrite);
    }
}

//================================================
void getDtofCut(const Int_t savePlot = 0, const Int_t saveHisto = 0)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s.JpsiMuon.root",run_type));

  TGraphAsymmErrors *gDtofCut = (TGraphAsymmErrors*)fin->Get("Data_JpsiMuon_DtofEff46_Fitting"); // get mean pt
  gDtofCut->SetName("Data_JpsiMuon_DtofCut_Fitting");

  const double goal = 0.95;
  const int nbins = 6;
  TF1 *func[nbins];
  TFitResultPtr *ptr[nbins];
  double x,y;
  for(int i=0; i<nbins; i++)
    {
      func[i] = (TF1*)fin->Get(Form("Data_JpsiMuon_Dtof_bin%d_Fit",i+1));
      ptr[i] = (TFitResultPtr*)fin->Get(Form("Data_JpsiMuon_Dtof_bin%d_Fit",i+1));
      double all = func[i]->Integral(-3,3);
      double cut = 3;
      for(int j=0; j<1000; j++)
	{
	  double val = 3 - j*1e-3*6;
	  double eff = func[i]->Integral(-3,val)/all;
	  if(eff<goal)
	    {
	      cut = val;
	      break;
	    }
	}
      gDtofCut->GetPoint(i,x,y);
      gDtofCut->SetPoint(i,x,cut);
      gDtofCut->SetPointError(i,0,0,0,0);
    }
  c = drawGraph(gDtofCut);
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



