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

  //DeltaTof();
}


//================================================
void DeltaTof(const Int_t savePlot = 0, const int saveHisto = 0)
{
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum = -3;
  const double maximum = 3;
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  const double dtof_cut_min = -3;
  const double dtof_cut_max = 0.25;

  TFile *fdata = 0x0;
  if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.JpsiMuon.root","read");
  if(year==2016) fdata = TFile::Open("./output/Pico.Run16.AuAu200.JpsiMuon.root","read");

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

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",run_type,bin),start_bin,end_bin);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMarkerStyle(20);

      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",run_type,bin),start_bin,end_bin);
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
  leg->SetHeader(run_type);
  leg->AddEntry(hUL[0],"Unlike-sign","PL");
  leg->AddEntry(hLS[0],"Like-sign","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataJpsiMuon_Dtof_MuonULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("DataJpsiMuon_Dtof_FitMean"),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",run_type),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("DataJpsiMuon_Dtof_FitSigma"),Form("%s: sigma of #Deltatof;p_{T} (GeV/c)",run_type),nbins,xbins);
  c = new TCanvas(Form("Fit_Dtof"),Form("Fit_Dtof"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));

      funcFitData[bin-1] = new TF1(Form("DataJpsiMuon_Dtof_bin%d_Fit",bin),CrystalBall,-3,3,5);
      funcFitData[bin-1]->SetParameters(-1.5, -0.05, 0.2, 50, hFit->GetMaximum());
      // funcFitData[bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,0.5);
      // funcFitData[bin-1]->SetParameters(hFit->GetMaximum(),-1,0.2);
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
      if(hMatch->GetBinContent(bin)<0.98) continue;
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
      if(hMatch->GetBinContent(ipoint+1)>0.98)
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
  if(year==2016)
    hplot->GetYaxis()->SetRangeUser(0.5,1.16);
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
  TGraphAsymmErrors *gplot = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("plot_%s",gCountDataEff->GetName()));
  gplot->Draw("PZsames");
  if(year==2016)
    {
      // for the lowest pT bin, use the efficiency from fitting method
      // gFitDataEff->GetPoint(0,x,y);
      // gCountDataEff->SetPoint(0,x,y);
      // gCountDataEff->SetPointError(0, gFitDataEff->GetErrorXlow(0), gFitDataEff->GetErrorXhigh(0), 
      // 				   gFitDataEff->GetErrorYlow(0), gFitDataEff->GetErrorYhigh(0));
      // gCountDataEff->SetMarkerStyle(25);
      // gCountDataEff->SetMarkerSize(1.3);
      // gCountDataEff->Draw("PZsames");
    }

  TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("Fit_%s",gCountDataEff->GetName()));
  TF1 *funcEff = new TF1(Form("%s_FitFunc",gCountDataEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  if(year==2016) funcEff->SetRange(1.5,7);
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
  leg->AddEntry(gplot,Form("Data: bin counting"),"p");
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
  leg->SetHeader(Form("%s",run_type));

  for(int l=0; l<nbins_rndm; l++)
    {
      hSpread[l] = new TH1F(Form("hSpread_%d",l),Form("hSpread_%d",l),600,0,1.2);
    }
  for(int j = 0; j < nexpr; j++)
    {
      TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
      if(year==2016) funcTmp->SetRange(1.5,7);
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
  t1->AddText(run_type);
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
      if(year==2016) funcLimits[t]->SetRange(1.5,7);
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
	  double new_y = funcEff->Eval(x)+sys_all*(i*2-1);
	  gDataEffSys[i]->SetPoint(ipoint,x,new_y);
	  gDataEffSys[i]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint]/y*new_y,eyh[ipoint]/y*new_y);	 
	}
      gDataEffSys[i]->GetYaxis()->SetRangeUser(0.6,1.05);
      c = drawGraph(gDataEffSys[i]);
      gPad->SetGridy();
      funcEffSys[i] = new TF1(Form("%s_FitFunc_Sys%s",gCountDataEff->GetName(),sys_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
      if(year==2016) funcEffSys[i]->SetRange(1.5,7);
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



