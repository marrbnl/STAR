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
}

//================================================
void MtdVpdTacDiff(const Int_t savePlot = 1, const int saveHisto = 1)
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
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdVpdTacDiff_ScaleFactor.pdf",run_type));


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
      hMuonFineBin[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("DataJpsiMuon_MtdVpdTacDiff_bin%d",bin));
      hMuonFineBin[bin-1]->Add(hLS[bin-1],-1.0*scale);

      // rebin
      TH1F *htmp1 = (TH1F*)hUL[bin-1]->Clone(Form("%s_rebin_tmp",hUL[bin-1]->GetName()));
      hUL[bin-1] = (TH1F*)htmp1->Rebin(nbins1,Form("%s_Rebin",hUL[bin-1]->GetName()),xbins1);
      scaleHisto(hUL[bin-1], 1, 1, true);
      htmp1 = (TH1F*)hLS[bin-1]->Clone(Form("%s_rebin_tmp",hLS[bin-1]->GetName()));
      hLS[bin-1] = (TH1F*)htmp1->Rebin(nbins1,Form("%s_Rebin",hLS[bin-1]->GetName()),xbins1);
      scaleHisto(hLS[bin-1], 1, 1, true);
      hMuon[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("DataJpsiMuon_MtdVpdTacDiff_bin%d_Rebin",bin));
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdVpdTacDiff_ULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("DataJpsiMuon_MtdVpdTacDiff_FitMean"),Form("%s: mean of MtdVpdTacDiff;p_{T} (GeV/c)",name),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("DataJpsiMuon_MtdVpdTacDiff_FitSigma"),Form("%s: sigma of MtdVpdTacDiff;p_{T} (GeV/c)",name),nbins,xbins);
  c = new TCanvas(Form("Fit_MtdVpdTacDiff"),Form("Fit_MtdVpdTacDiff"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));
      funcFitData[bin-1] = new TF1(Form("DataJpsiMuon_MtdVpdTacDiffFit_bin%d",bin),"gaus",789,maximum);
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdVpdTacDiffFit.pdf",run_type));

  hFitDataMean->SetMarkerStyle(21);
  hFitDataMean->GetYaxis()->SetRangeUser(780,810);
  c = draw1D(hFitDataMean,"Mean of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);Mean");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdVpdTacDiff_FitMean.pdf",run_type));

  hFitDataSigma->SetMarkerStyle(21);
  hFitDataSigma->GetYaxis()->SetRangeUser(0,14);
  c = draw1D(hFitDataSigma,"Sigma of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);#sigma");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdVpdTacDiff_FitSigma.pdf",run_type));

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
      gDataEff[k][0]->SetName(Form("DataJpsiMuon_MtdTrigEff_Fitting_%s",name_lumi[k]));
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
      gDataEff[k][1]->SetName(Form("DataJpsiMuon_MtdTrigEff_BinCount_%s",name_lumi[k]));
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
	cFuncBay[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEff_ProbDis_%s.pdf",run_type,name_lumi[k]));
    }
  
  // trigger efficiency
  const int defMed = 1;
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

      TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gDataEff[k][defMed]->Clone(Form("Fit_%s",gDataEff[k][defMed]->GetName()));
      funcEff[k] = new TF1(Form("%s_FitFunc",gDataEff[k][defMed]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
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
      if(savePlot) cTrigEff[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEff_%s.pdf",run_type,name_lumi[k]));
    }

  //==============================================
  // systematic uncertainty
  //==============================================
  // three sources:
  // i)   statistical error: randomization method
  // ii)  difference between default and central value of randomization
  // iii) level of asymmetry using Run15 data
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
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
  TGraphErrors *gLimits[nLumi][3];
  TF1 *funcLimits[nLumi][3];
  const char *limit_name[3] = {"center","low","up"};
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
      gLimits[k][2] = new TGraphErrors(nbins_rndm);
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

	  TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[k][l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
	  funcTmp2->SetParameter(0,hSpread[k][l]->GetMaximum());
	  funcTmp2->SetParameter(1,hSpread[k][l]->GetMean());
	  funcTmp2->SetParameter(2,hSpread[k][l]->GetRMS());
	  hSpread[k][l]->Fit(funcTmp2,"R0Q");
	  funcTmp2->SetLineColor(4);
	  funcTmp2->SetLineStyle(2);
	  funcTmp2->Draw("sames");

	  double pt    = xbins_rndm[l];
	  double mean  = hSpread[k][l]->GetMean();
	  double width = hSpread[k][l]->GetRMS();
	  double error = hSpread[k][l]->GetRMSError();

	  gLimits[k][0]->SetPoint(l,pt,mean);
	  gLimits[k][0]->SetPointError(l,0,error);
	  gLimits[k][1]->SetPoint(l,pt,mean-width);
	  gLimits[k][1]->SetPointError(l,0,error);
	  gLimits[k][2]->SetPoint(l,pt,mean+width);
	  gLimits[k][2]->SetPointError(l,0,error);

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
      if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEff_RndmFit_%s.pdf",run_type,name_lumi[k]));

      c->cd();
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEff_RndmCurve_%s.pdf",run_type,name_lumi[k]));
      for(int t=0; t<3; t++)
	{
	  gLimits[k][t]->SetMarkerStyle(20);
	  gLimits[k][t]->SetMarkerColor(2+2*t);
	  gLimits[k][t]->SetLineColor(2+2*t);
	  gLimits[k][t]->Draw("samesP");

	  funcLimits[k][t] = new TF1(Form("%s_Sys1_%s",gDataEff[k][defMed]->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
	  funcLimits[k][t]->SetParameters(0.95,2,0.2);
	  gLimits[k][t]->Fit(funcLimits[k][t],"R0Q");
	  funcLimits[k][t]->SetLineColor(gLimits[k][t]->GetMarkerColor());
	  if(t>0) funcLimits[k][t]->SetLineStyle(2);
	  funcLimits[k][t]->DrawCopy("sames");
	}
      leg->AddEntry(gLimits[k][0],"Central value","PL");
      leg->AddEntry(gLimits[k][1],"Lower limit","PL");
      leg->AddEntry(gLimits[k][2],"Upper limit","PL");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEff_RndmLimits_%s.pdf",run_type,name_lumi[k]));
    }

  TCanvas *cSys[nLumi];
  TGraphErrors *gSys[nLumi][3]; // three sources for prod_low/high
  for(int k=0; k<nLumi; k++)
    {
      int npoints = gDataEff[k][defMed]->GetN();
      for(int s=0; s<3; s++)
	{
	  gSys[k][s] = new TGraphErrors(npoints);
	  gSys[k][s]->SetName(Form("%s_S%d",gDataEff[k][defMed]->GetName(),s));
	  gSys[k][s]->SetTitle(Form("Systematic uncertainty %d (%s);p_{T} (GeV/c)",s+1,name_lumi[k]));
	  gSys[k][s]->SetMarkerStyle(20);
	  gSys[k][s]->SetMarkerSize(1.5);
	  gSys[k][s]->SetMarkerColor(s+1);
	  gSys[k][s]->SetLineColor(s+1);
	}

      // source 1
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  gDataEff[k][defMed]->GetPoint(ipoint, x, y);
	  double def   = funcLimits[k][0]->Eval(x);
	  double lower = funcLimits[k][1]->Eval(x);
	  double upper = funcLimits[k][2]->Eval(x);
	  double error = (upper-def) > (def-lower) ? (upper-def) : (def-lower);
	  gSys[k][0]->SetPoint(ipoint, x, 1);
	  gSys[k][0]->SetPointError(ipoint, 0, error/def);
	}

      // source 2
      for(int ipoint=0; ipoint<npoints; ipoint++)
	{
	  gDataEff[k][defMed]->GetPoint(ipoint, x, y);
	  double def  = funcEff[k]->Eval(x);
	  double diff = funcLimits[k][0]->Eval(x);
	  double error = fabs(diff/def-1);
	  gSys[k][1]->SetPoint(ipoint, x, 1);
	  gSys[k][1]->SetPointError(ipoint, 0, error);
	}

      gSys[k][0]->GetXaxis()->SetRangeUser(1.3,8);
      gSys[k][0]->GetYaxis()->SetRangeUser(0.95,1.05);
      cSys[k] = drawGraph(gSys[k][0],"Systematic uncertainty for MTD trigger efficiency");
      gPad->SetGridy();
      TGraphErrors *gSysShift = (TGraphErrors*)gSys[k][1]->Clone(Form("%s_shift",gSys[k][1]->GetName()));
      offset_x(gSysShift, 0.05);
      gSysShift->Draw("samesPE");
      TLegend *leg = new TLegend(0.4,0.15,0.65,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(name_lumi[k]);
      leg->AddEntry(gSys[k][0], "Source 1", "PE");
      leg->AddEntry(gSys[k][1], "Source 2", "PE");  
      leg->Draw();
      if(savePlot) cSys[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEffSys_%s.pdf",run_type,name_lumi[k]));
    }

  // combine the systematic uncertainty
  const char *sys_name[2] = {"down", "up"};
  TGraphAsymmErrors *gDataEffSys[nLumi][2];
  TF1 *funcEffSys[nLumi][2];
  for(int k=0; k<nLumi; k++)
    {
      double *exh = gDataEff[k][defMed]->GetEXhigh();
      double *exl = gDataEff[k][defMed]->GetEXlow();
      double *eyh = gDataEff[k][defMed]->GetEYhigh();
      double *eyl = gDataEff[k][defMed]->GetEYlow();
      for(int i=0; i<2; i++)
	{
	  int npoint = gDataEff[k][defMed]->GetN();
	  gDataEffSys[k][i] = new TGraphAsymmErrors(npoint);
	  gDataEffSys[k][i]->SetName(Form("%s_FinalSys%d",gDataEff[k][defMed]->GetName(),i));
	  for(int ipoint=0; ipoint<npoint; ipoint++)
	    {
	      gDataEff[k][defMed]->GetPoint(ipoint,x,y);
	      double sys_all = 0;
	      for(int s=0; s<2; s++)
		{
		  double y_sys = gSys[k][s]->GetErrorYhigh(ipoint);
		  sys_all += y_sys * y_sys;
		}
	      sys_all = sqrt(sys_all);
	      double new_y = y + sys_all*(i*2-1);
	      gDataEffSys[k][i]->SetPoint(ipoint,x,new_y);
	      gDataEffSys[k][i]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint]/y*new_y,eyh[ipoint]/y*new_y);	 
	    }
	  gDataEffSys[k][i]->GetYaxis()->SetRangeUser(0.1,1.05);
	  c = drawGraph(gDataEffSys[k][i]);
	  gPad->SetGridy();
	  funcEffSys[k][i] = new TF1(Form("%s_FitFunc_Sys%s",gDataEff[k][defMed]->GetName(),sys_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
	  gDataEffSys[k][i]->Fit(funcEffSys[k][i],"R0Q");
	  funcEffSys[k][i]->SetLineColor(2);
	  funcEffSys[k][i]->Draw("sames");
	}

      cTrigEff[k]->cd();
      for(int i=0; i<2; i++)
	{
	  funcEffSys[k][i]->SetLineStyle(2);
	  funcEffSys[k][i]->SetLineColor(6);
	  funcEffSys[k][i]->Draw("sames");
	}
      TLegend *leg = new TLegend(0.4,0.15,0.6,0.2);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(funcEffSys[k][0], "Systematic uncertainty", "L");
      leg->Draw();
      if(savePlot) cTrigEff[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_MtdTrigEffWithSys_%s.pdf",run_type,name_lumi[k]));
    }

  // Get the average efficiency
  // Statistical weight: prod_mid/high (73.4%); prod_/low (26.6%)
  const double weight[nLumi] = {0.734, 0.266};
  TF1 *funcEffFinal = new TF1("DataJpsiMuon_MtdTrigEff_BinCount_FitFunc","[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  TF1 *funcEffSysFinal[2];
  for(int i=0; i<2; i++) funcEffSysFinal[i] = new TF1(Form("DataJpsiMuon_MtdTrigEff_BinCount_FitFunc_Sys%s",sys_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  TCanvas *c = new TCanvas("Avg_Eff","Avg_eff",800,600);
  c->Divide(2,2);
  TF1 *hprodlow = 0x0, *hprodhigh = 0x0, *havg = 0x0;
  for(int i=0; i<3; i++)
    {
      if(i==0)
	{
	  hprodhigh = funcEff[0];
	  hprodlow = funcEff[1];
	  havg = funcEffFinal;
	}
      else
	{
	  hprodhigh = funcEffSys[0][i-1];
	  hprodlow = funcEffSys[1][i-1];
	  havg = funcEffSysFinal[i-1];
	}
      TH1F *htmp = new TH1F(Form("htmp_%d",i),"",70,1,8);
      for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	{
	  double pt = htmp->GetBinCenter(bin);
	  htmp->SetBinContent(bin, weight[0]*hprodhigh->Eval(pt)+weight[1]*hprodlow->Eval(pt));
	  htmp->SetBinError(bin, 1e-10);
	}
      TF1 *funcTmp = new TF1(Form("funcTmp_%d",i),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      funcTmp->SetParameters(hprodlow->GetParameters());
      c->cd(i+1);
      htmp->GetYaxis()->SetRangeUser(0.5,1);
      htmp->SetMarkerStyle(24);
      htmp->Fit(funcTmp,"IRQ");
      havg->SetParameters(funcTmp->GetParameters());
    }
  c = new TCanvas("MtdTrigEff_Avg","MtdTrigEff_Avg",800,600);
  gPad->SetGridy();
  hplot->DrawCopy();
  TPaveText *t1 = GetTitleText("Avergae MTD trigger efficiency");
  t1->Draw();
  funcEffFinal->Draw("sames");
  for(int i=0; i<2; i++) 
    {
      funcEffSysFinal[i]->SetLineColor(6);
      funcEffSysFinal[i]->SetLineStyle(2);
      funcEffSysFinal[i]->Draw("sames");
    }
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(funcEffFinal, "MTD trigger efficiency", "L");
  leg->AddEntry(funcEffSysFinal[0], "Systematic uncertainty", "L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_AvgMtdTrigEffWithSys.pdf",run_type));

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

      for(int k=0; k<2; k++)
	{
	  gDataEff[k][0]->Write("",TObject::kOverwrite);
	  gDataEff[k][1]->Write("",TObject::kOverwrite);
	  funcEff[k]->Write("",TObject::kOverwrite);
	  for(int t=0; t<2; t++)
	    {
	      funcEffSys[k][t]->Write("",TObject::kOverwrite);
	    }
	}
      funcEffFinal->Write("",TObject::kOverwrite);
      for(int t=0; t<2; t++)
	{
	  funcEffSysFinal[t]->Write("",TObject::kOverwrite);
	}
    }

  /*
  // source 3 
  const int nData = 2;
  const char *data_name[nData] = {"Run15_pp200", "Run15_pAu200"};
  TFile *fdata_Run15[nData];
  fdata_Run15[0] = TFile::Open("output/Pico.Run15.pp200.JpsiMuon.Dz3cm.root","read");
  fdata_Run15[1] = TFile::Open("output/Pico.Run15.pAu200.JpsiMuon.root","read");
  double Run15_mean_pt[nData][nbins];
  double Run15_mean_pt_err[nData][nbins];
  TH1F *hMuonTacDiffFineBin_Run15[nData][nbins];
  TH1F *hMuonTacDiff_Run15[nData][nbins];
  TF1 *funcFitData_Run15[nData][nbins][2];
  TFitResultPtr ptr_Run15[nData][nbins][2];
  double Run15_cut_min[nData][nLumi][nbins];
  double Run15_cut_max[nData][nLumi][nbins];
  const double Run15_fit_min[nData] = {880, 860};
  const double Run15_fit_max[nData] = {960, 960};
  const int Run15_rebin[nData] = {4, 8};
  const int nDataUsed = 1;
  // vertex distribution
  TH1F* hVzDiff[nData];
  for(int i=0; i<nData; i++)
    {
      TH2F* h2 = (TH2F*)fdata_Run15[i]->Get("mhVzDiffVsTpcVz_di_mu");
      h2->SetName(Form("%s_%s",data_name[i],h2->GetName()));
      hVzDiff[i] = (TH1F*)h2->ProjectionY(Form("%s_hVzDiff",data_name[i]));
      hVzDiff[i]->Scale(1./hVzDiff[i]->GetBinContent(hVzDiff[i]->FindFixBin(0)));
      hVzDiff[i]->SetLineWidth(2);
      // TF1 *func = new TF1(Form("func_%s",h1->GetName()),"gaus(0)+pol2(3)",-20,20);
      // func->SetParameter(0,h1->GetMaximum());
      // func->SetParameter(1,0);
      // func->SetParameter(2,2);
      // h1->Fit(func,"IR0");
      // func->Draw("sames");
    }
  draw1D(hVzDiff[1],"distribution of vz difference",false,false);
  gPad->SetLogy();
  hVzDiff[0]->SetLineStyle(2);
  hVzDiff[0]->SetLineColor(2);
  hVzDiff[0]->Draw("sames");
  TLegend *leg = new TLegend(0.15,0.65,0.35,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hVzDiff[0],data_name[0],"L");
  leg->AddEntry(hVzDiff[1],data_name[1],"L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Run15_CompareDz.pdf",run_type));

  for(int i=0; i<nDataUsed; i++)
    {
      getMuonTacDiff(savePlot, data_name[i], nbins, xbins, min_mass, max_mass, fdata_Run15[i], Run15_mean_pt[i], Run15_mean_pt_err[i], hMuonTacDiffFineBin_Run15[i]);

      //+++++ Fit muon distributions
      c = new TCanvas(Form("%s_Fit_MtdVpdTacDiff",data_name[i]),Form("%s_Fit_MtdVpdTacDiff",data_name[i]), 1100,700);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.15,0.6,0.5,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(data_name[i]);
      TLegend *leg1 = new TLegend(0.15,0.45,0.5,0.6);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.045);
      for(int bin=1; bin<=nbins; bin++)
	{
	  hMuonTacDiff_Run15[i][bin-1] = (TH1F*)hMuonTacDiffFineBin_Run15[i][bin-1]->Clone(Form("%s_Rebin",hMuonTacDiffFineBin_Run15[i][bin-1]->GetName()));
	  hMuonTacDiff_Run15[i][bin-1]->Rebin(Run15_rebin[i]);
	  if(bin==nbins) continue;

	  TH1F *hFit = (TH1F*)hMuonTacDiff_Run15[i][bin-1]->Clone(Form("Fit_%s",hMuonTacDiff_Run15[i][bin-1]->GetName()));
	  // Gaussian function
	  funcFitData_Run15[i][bin-1][0] = new TF1(Form("%s_MtdVpdTacDiffFit0_bin%d",data_name[i],bin),"gaus",Run15_fit_min[i]+30+(1-i)*5,Run15_fit_max[i]);
	  funcFitData_Run15[i][bin-1][0]->SetParameter(2,5);
	  ptr_Run15[i][bin-1][0] = hFit->Fit(funcFitData_Run15[i][bin-1][0],"IR0S");

	  // Crystal-ball function
	  funcFitData_Run15[i][bin-1][1] = new TF1(Form("%s_MtdVpdTacDiffFit1_bin%d",data_name[i],bin),CrystalBall,Run15_fit_min[i],Run15_fit_max[i],5);
	  if(i==0) 
	    {
	      funcFitData_Run15[i][bin-1][1]->SetParameters(1, 928, 8, 1, funcFitData_Run15[i][bin-1][0]->GetParameter(0));
	      funcFitData_Run15[i][bin-1][1]->SetParLimits(3, 0, 20);
	    }
	  if(i==1) 
	    {
	      if(bin==1) funcFitData_Run15[i][bin-1][1]->SetParameters(1, 907, 8, 0.01, funcFitData_Run15[i][bin-1][0]->GetParameter(0));
	      else       funcFitData_Run15[i][bin-1][1]->SetParameters(1.5, 907, 7, 0.1, funcFitData_Run15[i][bin-1][0]->GetParameter(0));
	    }
	  ptr_Run15[i][bin-1][1] = hFit->Fit(funcFitData_Run15[i][bin-1][1],"IR0S");

	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->SetTitle("");
	  hFit->SetMarkerStyle(20);
	  hFit->GetXaxis()->SetRangeUser(Run15_fit_min[i],Run15_fit_max[i]);
	  hFit->Draw();
	  funcFitData_Run15[i][bin-1][1]->SetLineColor(6);
	  funcFitData_Run15[i][bin-1][1]->Draw("sames");
	  TF1 *funcplot = (TF1*)funcFitData_Run15[i][bin-1][0]->Clone(Form("clone_%s",funcFitData_Run15[i][bin-1][0]->GetName()));
	  funcplot->SetRange(900,960);
	  funcplot->SetLineColor(4);
	  funcplot->SetLineStyle(2);
	  funcplot->Draw("sames");
	  funcFitData_Run15[i][bin-1][0]->SetLineColor(4);
	  funcFitData_Run15[i][bin-1][0]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();

	  if(bin==1)
	    {
	      leg->AddEntry(hFit,"Data","P");
	      leg->AddEntry(funcFitData_Run15[i][bin-1][0],"Gaussian fit","L");
	      leg->AddEntry(funcFitData_Run15[i][bin-1][1],"Crystal-ball fit","L");
	    }

	  // determine equivalent cut ranges
	  for(int k=0; k<nLumi; k++)
	    {
	      Run15_cut_min[i][k][bin-1] = funcFitData_Run15[i][bin-1][0]->GetParameter(1) - (funcFitData[bin-1]->GetParameter(1)-min[k])/funcFitData[bin-1]->GetParameter(2) * funcFitData_Run15[i][bin-1][0]->GetParameter(2);
	      Run15_cut_max[i][k][bin-1] = funcFitData_Run15[i][bin-1][0]->GetParameter(1) + (max[k]-funcFitData[bin-1]->GetParameter(1))/funcFitData[bin-1]->GetParameter(2) * funcFitData_Run15[i][bin-1][0]->GetParameter(2);
	      
	      TLine *line = GetLine(Run15_cut_min[i][k][bin-1], 0, Run15_cut_min[i][k][bin-1], 0.7*hFit->GetMaximum(), k+1);
	      line->Draw();
	      if(bin==2) leg1->AddEntry(line,Form("Cut for %s",name_lumi[k]),"L");
	    }
	}
      c->cd(nbins);
      leg->Draw();
      leg1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MuonTacDiffFit.pdf",run_type,data_name[i]));
    }

  //+++++ calculate efficiency
  TGraphAsymmErrors *gDataEff_Run15[nData][nLumi][2];
  TF1 *funcEff_Run15[nData][nLumi][2];
  for(int i=0; i<nDataUsed; i++)
    {
      c = new TCanvas(Form("%s_MtdTrigEff",data_name[i]),Form("%s_MtdTrigEff",data_name[i]),1100,500);
      c->Divide(2,1);
      for(int k=0; k<nLumi; k++)
	{
	  for(int j=0; j<2; j++)
	    {
	      TH1F *hBaseTmp  = new TH1F(Form("hBaseTmp_%d_%d_M%d",i,k,j),"",nbins,xbins);
	      TH1F *hMatchTmp = new TH1F(Form("hMatchTmp_%d_%d_M%d",i,k,j),"",nbins,xbins);
	      for(int bin=1; bin<=nbins-1; bin++)
		{
		  double all = funcFitData_Run15[i][bin-1][j]->Integral(860,980);
		  double all_err = funcFitData_Run15[i][bin-1][j]->IntegralError(860,980,funcFitData_Run15[i][bin-1][j]->GetParameters(),ptr_Run15[i][bin-1][j]->GetCovarianceMatrix().GetMatrixArray());
		  hBaseTmp->SetBinContent(bin,all);
		  hBaseTmp->SetBinError(bin,all_err);
		  double acc = funcFitData_Run15[i][bin-1][j]->Integral(Run15_cut_min[i][k][bin-1],980);
		  double acc_err = funcFitData_Run15[i][bin-1][j]->IntegralError(Run15_cut_min[i][k][bin-1],980,funcFitData_Run15[i][bin-1][j]->GetParameters(),ptr_Run15[i][bin-1][j]->GetCovarianceMatrix().GetMatrixArray());
		  hMatchTmp->SetBinContent(bin,acc);
		  hMatchTmp->SetBinError(bin,acc_err);
		}
	      gDataEff_Run15[i][k][j] = new TGraphAsymmErrors(hMatchTmp, hBaseTmp,"cl=0.683 b(1,1) mode");
	      gDataEff_Run15[i][k][j]->SetName(Form("%s_MtdTrigEff_Fit%d_%s",data_name[i],j,name_lumi[k]));
	      gDataEff_Run15[i][k][j]->SetMarkerStyle(20+j);
	      gDataEff_Run15[i][k][j]->SetMarkerColor(j+1);
	      gDataEff_Run15[i][k][j]->SetLineColor(j+1);
	      for(int ipoint=0; ipoint<nbins-1; ipoint++)
		{
		  gDataEff_Run15[i][k][j]->GetPoint(ipoint,x,y);
		  gDataEff_Run15[i][k][j]->SetPoint(ipoint,Run15_mean_pt[i][ipoint],y);
		  gDataEff_Run15[i][k][j]->SetPointEXhigh(ipoint,Run15_mean_pt_err[i][ipoint]);
		  gDataEff_Run15[i][k][j]->SetPointEXlow(ipoint,Run15_mean_pt_err[i][ipoint]);
		}
	      funcEff_Run15[i][k][j] = new TF1(Form("%s_FitFunc",gDataEff_Run15[i][k][j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.3,5);
	      funcEff_Run15[i][k][j]->SetParameters(0.9, 5, 1);
	      gDataEff_Run15[i][k][j]->Fit(funcEff_Run15[i][k][j],"R0");
	    }
	  c->cd(k+1);
	  gDataEff_Run15[i][k][0]->GetYaxis()->SetRangeUser(0.55,1);
	  gDataEff_Run15[i][k][0]->SetTitle(Form(";p_{T}^{#mu} (GeV/c);"));
	  gDataEff_Run15[i][k][0]->Draw("AZP");
	  gDataEff_Run15[i][k][1]->Draw("samesPEZ");
	  for(int j=0; j<2; j++)
	    {
	      funcEff_Run15[i][k][j]->SetLineStyle(2);
	      funcEff_Run15[i][k][j]->SetLineColor(j+1);
	      funcEff_Run15[i][k][j]->Draw("sames");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s: MTD trigger efficiency for %s cut",data_name[i],name_lumi[k]),0.04);
	  t1->Draw();
	  if(k==0)
	    {
	      TLegend *leg = new TLegend(0.5,0.2,0.7,0.35);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.04);
	      leg->AddEntry(gDataEff_Run15[i][k][0], "Gaussian fit", "P");
	      leg->AddEntry(gDataEff_Run15[i][k][1], "Crystal-ball fit", "P");
	      leg->Draw();
	    }
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdTrigEff.pdf",run_type,data_name[i]));
    }
*/

}


//================================================
void getMuonTacDiff(const int savePlot, const char *name, const int nbins, const double *xbins, 
		    const double *min_mass, const double *max_mass,
		    TFile *fdata, double *mean_pt, double *mean_pt_err, TH1F **hMuonDis)
{
  TH1F *htmp = 0x0;

  //==============================================
  // compare invariant mass
  //==============================================
  const int nbinsJpsi = 5;
  const double xbinsJpsi[nbinsJpsi+1] = {0,1.0,2.0,3.0,5.0,10.0};
  TH2F *hInvMassVsPtUL = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_MtdVpdTacDiff_UL_di_mu"));
  TH2F *hInvMassVsPtLS = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_MtdVpdTacDiff_LS_di_mu"));
  hInvMassVsPtUL->Sumw2();
  hInvMassVsPtLS->Sumw2();

  TCanvas *c = new TCanvas(Form("%s_InvMass_MtdVpdTacDiff",name),Form("%s_InvMass_MtdVpdTacDiff",name),1100,700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetHeader(name);
  for(int bin=1; bin<=nbinsJpsi; bin++)
    {
      int start_bin = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin-1]+1e-4);
      int end_bin   = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin]-1e-4);
      TH1F *hInvMassUL = (TH1F*)hInvMassVsPtUL->ProjectionY(Form("%s_InvMassUL_bin%d",name,bin),start_bin,end_bin);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->Rebin(5);
      hInvMassUL->SetMaximum(1.5*hInvMassUL->GetMaximum());
      hInvMassUL->SetMinimum(0);
      hInvMassUL->GetXaxis()->SetRangeUser(2.2,4);
      hInvMassUL->SetTitle("");
      if(bin==1) leg->AddEntry(hInvMassUL,"Unlike-sign","PL");

      TH1F *hInvMassLS = (TH1F*)hInvMassVsPtLS->ProjectionY(Form("%s_InvMassLS_bin%d",name,bin),start_bin,end_bin);
      hInvMassLS->SetMarkerStyle(24);
      hInvMassLS->SetMarkerColor(2);
      hInvMassLS->SetLineColor(2);
      hInvMassLS->Rebin(5);
      if(bin==1) leg->AddEntry(hInvMassLS,"Like-sign","PL");

      c->cd(bin);
      hInvMassUL->Draw("P");
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
      t1->AddText(Form("S/B = %4.2f:1",(all-bkg)/bkg));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MuonTacDiff_InvMassLSvsUL.pdf",run_type,name));

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

  TH1F *hScaleFactor = new TH1F(Form("%s_hScaleFactor",name),"Scale factor for unlike-sign background (online trigger efficiency);p_{T} (GeV/c)",nbins,xbins);
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


  TH2F *hDataUL = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_UL_di_mu");
  hDataUL->Sumw2();
  TH2F *hDataLS = (TH2F*)fdata->Get("mhJpsiMuon_MtdVpdTacDiff_LS_di_mu");
  hDataLS->Sumw2();
  TH2F *hDataDisVsPt = (TH2F*)hDataUL->Clone(Form("JpsiMuonMtdVpdTacDiff"));
  hDataDisVsPt->Add(hDataLS, -1);
  TH1F *hUL[nbins];
  TH1F *hLS[nbins];

  // calculate mean pt of each bin
  int bin1 = hDataDisVsPt->GetYaxis()->FindBin(500);
  int bin2 = hDataDisVsPt->GetYaxis()->FindBin(1000);
  htmp = (TH1F*)hDataDisVsPt->ProjectionX(Form("%s_tmp",hDataDisVsPt->GetName()),bin1,bin2);
  TCanvas *c = new TCanvas(Form("%s_MtdVpdTacDiff_USvsLS",name),Form("%s_MtdVpdTacDiff_USvsLS",name),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      htmp->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
      mean_pt[bin-1] = htmp->GetMean();
      mean_pt_err[bin-1] = htmp->GetMeanError();

      int start_bin = hDataUL->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hDataUL->GetXaxis()->FindBin(xbins[bin]-1e-4);

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",name,bin),start_bin,end_bin);
      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",name,bin),start_bin,end_bin);
      double scale = hScaleFactor->GetBinContent(bin);
      hMuonDis[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("%s_DataJpsiMuon_MtdVpdTacDiff_bin%d",name,bin));
      hMuonDis[bin-1]->Add(hLS[bin-1],-1.0*scale);
      //hMuonDis[bin-1]->Add(hLS[bin-1],-1.0);

      c->cd(bin);
      hUL[bin-1]->Rebin(4);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMaximum(1.5*hUL[bin-1]->GetMaximum());
      hUL[bin-1]->GetXaxis()->SetRangeUser(870,960);
      hUL[bin-1]->SetTitle("");
      hUL[bin-1]->Draw("P");
      hLS[bin-1]->Rebin(4);
      hLS[bin-1]->SetMarkerStyle(24);
      hLS[bin-1]->SetMarkerColor(2);
      hLS[bin-1]->SetLineColor(2);
      hLS[bin-1]->Draw("samesP");
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
      t1 = GetPaveText(0.6,0.8,0.8,0.85,0.05);
      t1->AddText(Form("Scale = %4.3f",scale));
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.15,0.63,0.35,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(name);
  leg->AddEntry(hUL[0],"Unlike-sign","PL");
  leg->AddEntry(hLS[0],"Like-sign","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MuonTacDiff_ULvsLS.pdf",run_type,name));
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
  const double dtof_cut_min = -3;
  const double dtof_cut_max = 0.75;

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
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_InvMassLSvsUL.pdf",run_type));

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
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_ScaleFactor.pdf",run_type));


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
      hMuonFineBin[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("DataJpsiMuon_Dtof_bin%d",bin));
      hMuonFineBin[bin-1]->Add(hLS[bin-1],-1.0*scale);

      // rebin
      hMuon[bin-1] = (TH1F*)hMuonFineBin[bin-1]->Clone(Form("DataJpsiMuon_Dtof_bin%d_Rebin",bin));
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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_MuonULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("DataJpsiMuon_Dtof_FitMean"),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",name),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("DataJpsiMuon_Dtof_FitSigma"),Form("%s: sigma of #Deltatof;p_{T} (GeV/c)",name),nbins,xbins);
  c = new TCanvas(Form("Fit_Dtof"),Form("Fit_Dtof"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));

      funcFitData[bin-1] = new TF1(Form("DataJpsiMuon_Dtof_bin%d_Fit",bin),CrystalBall,-3,3,5);
      funcFitData[bin-1]->SetParameters(-1.5, -0.05, 0.2, 50, hFit->GetMaximum());
      //funcFitData[bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,3);
      ptr[bin-1] = hFit->Fit(funcFitData[bin-1],"IR0QS");
      ptr[bin-1]->SetName(Form("DataJpsiMuon_Dtof_bin%d_FitResult",bin));

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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_Fit.pdf",run_type));

  hFitDataMean->SetMarkerStyle(21);
  hFitDataMean->GetYaxis()->SetRangeUser(-1,1);
  c = draw1D(hFitDataMean,"Mean of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);Mean");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_FitMean.pdf",run_type));

  hFitDataSigma->SetMarkerStyle(21);
  hFitDataSigma->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hFitDataSigma,"Width of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);#sigma");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_FitSigma.pdf",run_type));

  //==============================================
  // extract efficiency
  //==============================================
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

  // recalculate statistical errors
  hMatch->Divide(hBase);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int nRndm = 1e3;
  for(int bin=1; bin<=nbins; bin++)
    {
      if(hMatch->GetBinContent(bin)<1) continue;
      TH1F *hRndmDis = new TH1F(Form("RndmDis_Bin%d",bin),"",50,0.8,1.2);
      double all, all_err, acc, acc_err;
      int all_low_bin = hMuonFineBin[bin-1]->FindFixBin(minimum+1e-4);
      int all_high_bin = hMuonFineBin[bin-1]->FindFixBin(maximum-1e-4);
      int acc_low_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_min+1e-4);
      int acc_high_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_max-1e-4);
      int nbinsx = hMuonFineBin[bin-1]->GetNbinsX();
      for(int i=0; i<nRndm; i++)
	{
	  TH1F *htmp = (TH1F*)hMuonFineBin[bin-1]->Clone(Form("%s_rndm%d", hMuonFineBin[bin-1]->GetName(), i));
	  htmp->Reset();
	  for(int ibin=1; ibin<=nbinsx; ibin++)
	    {
	      double mean = hMuonFineBin[bin-1]->GetBinContent(ibin);
	      double error = hMuonFineBin[bin-1]->GetBinError(ibin);
	      htmp->SetBinContent(ibin, rndm->Gaus(mean, error));
	      htmp->SetBinError(ibin, error);
	    }
	  all = htmp->IntegralAndError(all_low_bin, all_high_bin, all_err);
	  acc = htmp->IntegralAndError(acc_low_bin, acc_high_bin, acc_err);
	  hRndmDis->Fill(acc/all);
	}
      hRndmDis->SetMarkerStyle(21);
      c = draw1D(hRndmDis,Form("Randomized efficiency for %1.1f < p_{T}^{#mu} < %1.1f GeV/c;Efficiency",xbins[bin-1],xbins[bin]));
      TF1 *funcRndm = new TF1(Form("Fit_%s",hRndmDis->GetName()),"gaus",0.8,1.2);
      hRndmDis->Fit(funcRndm,"IR0Q");
      funcRndm->SetLineColor(2);
      funcRndm->Draw("sames");
      TPaveText *t1 = GetPaveText(0.2,0.35,0.7,0.85,0.04);
      t1->SetTextAlign(11);
      t1->AddText(Form("RMS=%4.3f",hRndmDis->GetRMS()));
      t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcRndm->GetParameter(2))));
      t1->Draw();
      printf("[i] Check the mean for bin %d: %4.2f =? %4.2f\n",bin, hMatch->GetBinContent(bin), funcRndm->GetParameter(1));
      hMatch->SetBinError(bin, hRndmDis->GetRMS());
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_RndmBin%d.pdf",run_type,dtof_cut_max*100,bin));
    }
  
  TGraphAsymmErrors *gCountDataEff = new TGraphAsymmErrors(hMatch2, hBase,"cl=0.683 b(1,1) mode");
  gCountDataEff->SetName(Form("DataJpsiMuon_DtofEff%1.0f_BinCount",dtof_cut_max*100));
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint,x,y);
      gCountDataEff->SetPoint(ipoint,mean_pt[ipoint],y);
      gCountDataEff->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
      gCountDataEff->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
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
      double acc_err = funcFitData[bin-1]->IntegralError(dtof_cut_min, dtof_cut_max,funcFitData[bin-1]->GetParameters(), ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
      hBase->SetBinContent(bin,all);
      hBase->SetBinError(bin,all_err);
      hMatch->SetBinContent(bin,acc);
      hMatch->SetBinError(bin,acc_err);
    }
  TGraphAsymmErrors *gFitDataEff = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
  gFitDataEff->SetName(Form("DataJpsiMuon_DtofEff%1.0f_Fitting",dtof_cut_max*100));
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
  TCanvas *cEff = new TCanvas("DtofEff_Comp","DtofEff_Comp",800,600);
  gPad->SetGridy();
  hplot->DrawCopy();
  gFitDataEff->SetMarkerSize(1.3);
  gFitDataEff->SetMarkerStyle(24);
  gFitDataEff->SetMarkerColor(2);
  gFitDataEff->SetLineColor(2);
  gFitDataEff->Draw("PZsames");
  gCountDataEff->SetMarkerStyle(21);
  gCountDataEff->SetMarkerSize(1.3);
  gCountDataEff->Draw("PZsames");
  TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("Fit_%s",gCountDataEff->GetName()));
  TF1 *funcEff = new TF1(Form("%s_FitFunc",gCountDataEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  gfit->Fit(funcEff,"R0");
  funcEff->SetLineColor(2);
  funcEff->Draw("sames");

  TPaveText *t1 = GetTitleText(Form("Efficiency for %1.0f < #DeltaTof < %1.2f ns cut",dtof_cut_min,dtof_cut_max),0.04);
  t1->Draw();
  TLegend *leg = new TLegend(0.45,0.25,0.6,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  leg->AddEntry(gFitDataEff,Form("Data: fitting"),"p");
  leg->AddEntry(gCountDataEff,Form("Data: bin counting"),"p");
  leg->AddEntry(funcEff,"Fit to bin counting","L");
  leg->Draw();
  if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_FitVsCount.pdf",run_type,dtof_cut_max*100));


  //==============================================
  // systematic uncertainty
  //==============================================
  // three sources:
  // i)   statistical error: randomization method
  // ii)  difference between default and central value of randomization
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int nexpr = 400;
  TGraphAsymmErrors *gEffRndm[nexpr];
  double *exh = gCountDataEff->GetEXhigh();
  double *exl = gCountDataEff->GetEXlow();
  double *eyh = gCountDataEff->GetEYhigh();
  double *eyl = gCountDataEff->GetEYlow();
  for(int j = 0; j < nexpr; j++)
    {
      int npoint = gCountDataEff->GetN();
      gEffRndm[j] = new TGraphAsymmErrors(npoint);
      gEffRndm[j]->SetName(Form("%s_Rndm%d",gCountDataEff->GetName(),j));
      for(int ipoint=0; ipoint<npoint; ipoint++)
	{
	  gCountDataEff->GetPoint(ipoint,x,y);
	  double sigma = eyl[ipoint] > eyh[ipoint] ? eyl[ipoint] : eyl[ipoint];
	  y = rndm->Gaus(y, sigma);
	  gEffRndm[j]->SetPoint(ipoint,x,y);
	  gEffRndm[j]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	}
    }

  hplot->GetYaxis()->SetRangeUser(0.6,1.1);
  hplot->GetXaxis()->SetRangeUser(1.3,7);
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1.4,1.5,1.75,2,2.25,2.5,3,4,5,6,8};
  TH1F *hSpread[nbins_rndm];
  TGraphErrors *gLimits[3];
  TF1 *funcLimits[3];
  const char *limit_name[3] = {"center","low","up"};

  c = new TCanvas(Form("cEff"),Form("cEff"),800,600);
  hplot->DrawCopy();
  TPaveText *t1 = GetTitleText("Estimation of trigger efficiency uncertainty");
  t1->Draw();
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%s",title));

  for(int l=0; l<nbins_rndm; l++)
    {
      hSpread[l] = new TH1F(Form("hSpread_%d",l),Form("hSpread_%d",l),600,0,1.2);
    }
  for(int j = 0; j < nexpr; j++)
    {
      TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
      funcTmp->SetParLimits(0,0,1);
      funcTmp->SetParLimits(1,0,5);
      funcTmp->SetParLimits(2,0,1);
      gEffRndm[j]->Fit(funcTmp,"R0Q");
      c->cd();
      funcTmp->SetLineStyle(2);
      funcTmp->Draw("sames");
      if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
      for(int l=0; l<nbins_rndm; l++)
	{
	  hSpread[l]->Fill(funcTmp->Eval(xbins_rndm[l]));
	}
    }

  TCanvas *cFit = new TCanvas(Form("cEffFit"),Form("cEffFit"),1200,800);
  cFit->Divide(4,3);
  for(int t=0; t<3; t++) gLimits[t] = new TGraphErrors(nbins_rndm);
  for(int l=0; l<nbins_rndm; l++)
    {
      cFit->cd(l+1);
      hSpread[l]->Rebin(2);
      if(l<2) hSpread[l]->GetXaxis()->SetRangeUser(0.3,1);
      else if(l<4) hSpread[l]->GetXaxis()->SetRangeUser(0.65,1);
      else hSpread[l]->GetXaxis()->SetRangeUser(0.75,1);
      hSpread[l]->SetTitle(";Efficiency;");
      hSpread[l]->SetMaximum(1.3*hSpread[l]->GetMaximum());
      hSpread[l]->Draw();
      TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[l]),0.06);
      t1->Draw();
      TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
      funcTmp2->SetParameter(0,hSpread[l]->GetMaximum());
      funcTmp2->SetParameter(1,hSpread[l]->GetMean());
      funcTmp2->SetParameter(2,hSpread[l]->GetRMS());
      hSpread[l]->Fit(funcTmp2,"R0Q");
      funcTmp2->SetLineColor(4);
      funcTmp2->SetLineStyle(2);
      funcTmp2->Draw("sames");

      double pt    = xbins_rndm[l];
      double mean  = hSpread[l]->GetMean();
      double width = hSpread[l]->GetRMS();
      double error = hSpread[l]->GetRMSError();
      gLimits[0]->SetPoint(l,pt,mean);
      gLimits[0]->SetPointError(l,0,error);
      gLimits[1]->SetPoint(l,pt,mean-width);
      gLimits[1]->SetPointError(l,0,error);
      gLimits[2]->SetPoint(l,pt,mean+width);
      gLimits[2]->SetPointError(l,0,error);

      TPaveText *t1 = GetPaveText(0.6,0.7,0.7,0.85,0.06);
      t1->SetTextAlign(11);
      t1->AddText(Form("RMS=%4.3f",width));
      t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcTmp2->GetParameter(2))));
      t1->Draw();
    }
  cFit->cd(12);
  TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
  t1->AddText(title);
  t1->Draw();
  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_RndmFit.pdf",run_type,dtof_cut_max*100));
  c->cd();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_RndmCurve.pdf",run_type,dtof_cut_max*100));
  for(int t=0; t<3; t++)
    {
      gLimits[t]->SetMarkerStyle(20);
      gLimits[t]->SetMarkerColor(2+2*t);
      gLimits[t]->SetLineColor(2+2*t);
      gLimits[t]->Draw("samesP");

      funcLimits[t] = new TF1(Form("%s_Sys1_%s",gCountDataEff->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
      funcLimits[t]->SetParameters(0.95,2,0.2);
      gLimits[t]->Fit(funcLimits[t],"R0Q");
      funcLimits[t]->SetLineColor(gLimits[t]->GetMarkerColor());
      if(t>0) funcLimits[t]->SetLineStyle(2);
      funcLimits[t]->DrawCopy("sames");
    }
  leg->AddEntry(gLimits[0],"Central value","PL");
  leg->AddEntry(gLimits[1],"Lower limit","PL");
  leg->AddEntry(gLimits[2],"Upper limit","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_RndmLimits.pdf",run_type,dtof_cut_max*100));

  int npoints = gCountDataEff->GetN();
  TGraphErrors *gSys[2];
  for(int s=0; s<2; s++)
    {
      gSys[s] = new TGraphErrors(npoints);
      gSys[s]->SetName(Form("%s_S%d",gCountDataEff->GetName(),s));
      gSys[s]->SetTitle(Form("Systematic uncertainty;p_{T} (GeV/c)"));
      gSys[s]->SetMarkerStyle(20);
      gSys[s]->SetMarkerColor(s+1);
      gSys[s]->SetMarkerSize(1.3);
      gSys[s]->SetLineColor(s+1);
    }
  // source 1
  for(int ipoint=0; ipoint<npoints; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint, x, y);
      double def   = funcLimits[0]->Eval(x);
      double lower = funcLimits[1]->Eval(x);
      double upper = funcLimits[2]->Eval(x);
      double error = (upper-def) > (def-lower) ? (upper-def) : (def-lower);
      gSys[0]->SetPoint(ipoint, x, 1);
      gSys[0]->SetPointError(ipoint, 0, error/def);
    }
  // source 2
  for(int ipoint=0; ipoint<npoints; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint, x, y);
      double def  = funcEff->Eval(x);
      double diff = funcLimits[0]->Eval(x);
      double error = fabs(diff/def-1);
      gSys[1]->SetPoint(ipoint, x, 1);
      gSys[1]->SetPointError(ipoint, 0, error);
    }
  gSys[0]->GetXaxis()->SetRangeUser(1.3,8);
  gSys[0]->GetYaxis()->SetRangeUser(0.94,1.06);
  c = drawGraph(gSys[0],"Systematic uncertainty for #Deltatof cut efficiency");
  gPad->SetGridy();
  TGraphErrors *gSysShift = (TGraphErrors*)gSys[1]->Clone(Form("%s_shift",gSys[1]->GetName()));
  offset_x(gSysShift, 0.05);
  gSysShift->Draw("samesPE");
  TLegend *leg = new TLegend(0.4,0.15,0.6,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gSys[0], "Source 1", "PE");
  leg->AddEntry(gSys[1], "Source 2", "PE");  
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_Sys.pdf",run_type,dtof_cut_max*100));

  // combine the systematic uncertainty
  const char *sys_name[2] = {"down", "up"};
  TGraphAsymmErrors *gDataEffSys[2];
  TF1 *funcEffSys[2];
  for(int i=0; i<2; i++)
    {
      int npoint = gCountDataEff->GetN();
      gDataEffSys[i] = new TGraphAsymmErrors(npoint);
      gDataEffSys[i]->SetName(Form("%s_FinalSys%d",gCountDataEff->GetName(),i));
      for(int ipoint=0; ipoint<npoint; ipoint++)
	{
	  gCountDataEff->GetPoint(ipoint,x,y);
	  double sys_all = 0;
	  for(int s=0; s<2; s++)
	    {
	      double y_sys = gSys[s]->GetErrorYhigh(ipoint);
	      sys_all += y_sys * y_sys;
	    }
	  sys_all = sqrt(sys_all);
	  double new_y = y+sys_all*(i*2-1);
	  gDataEffSys[i]->SetPoint(ipoint,x,new_y);
	  gDataEffSys[i]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint]/y*new_y,eyh[ipoint]/y*new_y);	 
	}
      gDataEffSys[i]->GetYaxis()->SetRangeUser(0.6,1.05);
      c = drawGraph(gDataEffSys[i]);
      gPad->SetGridy();
      funcEffSys[i] = new TF1(Form("%s_FitFunc_Sys%s",gCountDataEff->GetName(),sys_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
      funcEffSys[i]->SetParLimits(0, 0, 1);
      funcEffSys[i]->SetParameter(1, 3.8);
      funcEffSys[i]->SetParameter(2, 0.85);
      gDataEffSys[i]->Fit(funcEffSys[i],"R0");
      funcEffSys[i]->SetLineColor(2);
      funcEffSys[i]->Draw("sames");
    }

  cEff->cd();
  for(int i=0; i<2; i++)
    {
      funcEffSys[i]->SetLineStyle(2);
      funcEffSys[i]->SetLineColor(6);
      funcEffSys[i]->Draw("sames");
    }
  TLegend *leg = new TLegend(0.45,0.2,0.6,0.25);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(funcEffSys[0], "Systematic uncertainty", "L");
  leg->Draw();
  if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_DtofEff%1.0f_FitVsCountWithSys.pdf",run_type,dtof_cut_max*100));
  

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
      for(int i=0; i<2; i++)
	{
	  funcEffSys[i]->Write("",TObject::kOverwrite);
	}
    }
}



