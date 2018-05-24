const int year = YEAR;

const char* dataType[2] = {"Embed","Cosmic"};
const int nVz = 5;
const double vzMin[nVz] = {-100, -50, -4, 40, 90};
const double vzMax[nVz] = {-90,  -40, 4,  50, 100};
const int nCellBin = 5;
const int start_cell[nCellBin] = {1, 3, 4, 16, 17};
const int end_cell[nCellBin]   = {2, 3, 15, 16, 18};
const int nEta = 8;
const int nPhiDiffCut = 6;
const double phiDiffCuts[nPhiDiffCut] = {0.025, 0.05, 0.1, 0.15, 0.2, 1};
const int nYears = 3;
const int years[4] = {2014, 2015, 2016, 2017};
const int nBtmBL = 13;
const double ptCut = 5;
const double vzCut = 100;
const double etaCut = 0.4;
const double phiDiffCut = 0.05;
const double pi = TMath::Pi();
const int nPtBins = 20;
const double xPtBins[nPtBins+1] = {0,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0,8.0,10,15};
const int colors[10] = {1, 2, 4, 6, 8, 1, 2, 4, 6, 8};

//================================================
void ana_MtdMthEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  studyMthEff();
  //makeHisto();
  //makeHisto2();
  //makeHisto3();
  //respEff();

}

//================================================
void studyMthEff(const int savePlot = 1, const int saveHisto = 1)
{
  gStyle->SetOptFit(0);
  TFile *fin = TFile::Open(Form("Rootfiles/%s.MtdMthEff.root",run_type), "read");
  TString legName[2] = {"Embedding","Cosmic ray"};
  TH1F *h1tmp = NULL;

  // phi correlation
  TH2F *hPosVsMomPhi[2];
  TH1F *hPhiDiff[2];
  TCanvas *c = new TCanvas("cPosVsMomPhi", "cPosVsMomPhi", 1100, 700);
  c->Divide(2,2);
  for(int i=0; i<2; i++)
    {
      c->cd(i+1);
      gPad->SetLogz();
      hPosVsMomPhi[i] = (TH2F*)fin->Get(Form("mhPosVsMomPhi_Type%d",i));
      hPosVsMomPhi[i]->SetTitle("");
      ScaleHistoTitle(hPosVsMomPhi[i],0.05,0.8,0.045,0.05,1,0.045);
      hPosVsMomPhi[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%s: #varphi correlation (p_{T} > 5 GeV/c)",dataType[i]),0.055);
      t1->Draw();
      
      hPhiDiff[i] = (TH1F*)fin->Get(Form("mhPosMomPhiDiff_Type%d",i));
      c->cd(i+3);
      ScaleHistoTitle(hPhiDiff[i],0.05,0.9,0.045,0.05,1,0.045);
      hPhiDiff[i]->GetXaxis()->SetRangeUser(-1,1);
      hPhiDiff[i]->SetTitle("");
      hPhiDiff[i]->Draw();
      t1 = GetTitleText(Form("%s: #varphi difference distribution (p_{T} > 5 GeV/c)",dataType[i]),0.055);
      t1->Draw();
    }
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/PosVsMomPhi.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_PosVsMomPhi.pdf"));
    }

  // eff vs. track eta and vz
  printf("+++++ histograms vs. eta and vz +++++\n");
  TH1F *hTrkPtVsVzAll[2][nEta];
  TH1F *hTrkPtVsVzAcc[2][nEta];
  TH1F *hTrkPtVsVzEff[2][nEta];
  TCanvas *c = new TCanvas("EffVsVz_InEtaBin","EffVsVz_InEtaBin",1100,500);
  c->Divide(2,1);
  TLegend *leg1 = new TLegend(0.2,0.15,0.4,0.35);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  TLegend *leg2 = new TLegend(0.5,0.15,0.7,0.35);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.04);
  for(int i=0; i<2; i++)
    {
      c->cd(i+1);
      for(int j=0; j<nEta; j++)
	{
	  hTrkPtVsVzAll[i][j] = (TH1F*)fin->Get(Form("hTrkPtVsVzAll_pt%1.0f_eta%d_Type%d",ptCut,j,i));
	  hTrkPtVsVzAcc[i][j] = (TH1F*)fin->Get(Form("hTrkPtVsVzAcc_pt%1.0f_eta%d_Type%d",ptCut,j,i));
	  hTrkPtVsVzAll[i][j]->Rebin(10);
	  hTrkPtVsVzAcc[i][j]->Rebin(10);
	  hTrkPtVsVzEff[i][j] = (TH1F*)hTrkPtVsVzAcc[i][j]->Clone(Form("hTrkPtVsVzEff_pt%1.0f_eta%d_Type%d",ptCut,j,i));
	  hTrkPtVsVzEff[i][j]->Divide(hTrkPtVsVzAll[i][j]);
	  for(int bin=1; bin<=hTrkPtVsVzEff[i][j]->GetNbinsX(); bin++)
	    {
	      if(hTrkPtVsVzEff[i][j]->GetBinError(bin)/hTrkPtVsVzEff[i][j]->GetBinContent(bin)>0.5)
		{
		  hTrkPtVsVzEff[i][j]->SetBinContent(bin, 0);
		  hTrkPtVsVzEff[i][j]->SetBinError(bin, 0);
		}
	    }
	  hTrkPtVsVzEff[i][j]->SetMarkerStyle(20+j);
	  hTrkPtVsVzEff[i][j]->SetMarkerColor(colors[j%4]);
	  hTrkPtVsVzEff[i][j]->SetLineColor(colors[j%4]);
	  hTrkPtVsVzEff[i][j]->GetYaxis()->SetRangeUser(0.7,1.1);
	  hTrkPtVsVzEff[i][j]->GetYaxis()->SetTitle("Resp. Eff.");
	  hTrkPtVsVzEff[i][j]->SetTitle("");
	  if(j==0) hTrkPtVsVzEff[i][j]->Draw();
	  else     hTrkPtVsVzEff[i][j]->Draw("sames");
	  if(i==0)
	    {
	      if(j<nEta/2) leg1->AddEntry(hTrkPtVsVzEff[i][j],Form("%1.1f < #eta < %1.1f",-0.8+j*0.2,-0.6+j*0.2),"P");
	      else         leg2->AddEntry(hTrkPtVsVzEff[i][j],Form("%1.1f < #eta < %1.1f",-0.8+j*0.2,-0.6+j*0.2),"P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%d: %s (p_{T} > %1.0f GeV/c)",year,legName[i].Data(),ptCut),0.045);
      t1->Draw();
    }
  c->cd(1);
  leg1->Draw();
  leg2->Draw();
  c->cd(2);
  TPaveText *t1 = GetPaveText(0.2,0.5,0.8,0.85,0.045);
  t1->AddText(Form("|#varphi_{mom} - #varphi_{pos}| < %2.2f",phiDiffCut));
  t1->SetTextColor(4);
  t1->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthEff_VsVzInEtaBin.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_VsVzInEtaBin.pdf"));
    }
  printf("+++++ Done +++++\n\n");

  TCanvas *c = new TCanvas("MtdRespMthEff_InVzBin","MtdRespMthEff_InVzBin",1100,700);
  c->Divide(2,2);
  TH1F *hTrkPtRespInVzAll[2][nVz];
  TH1F *hTrkPtRespInVzAcc[2][nVz];
  TH1F *hTrkPtRespInVzEff[2][nVz];
  TH1F *hTrkPtMthInVzAll[2][nVz];
  TH1F *hTrkPtMthInVzAcc[2][nVz];
  TH1F *hTrkPtMthInVzEff[2][nVz];
  for(int i=0; i<2; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad,0.13,0.1,0.1,0.1);
      TLegend *leg = new TLegend(0.5,0.2,0.7,0.6);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      if(i==0) leg->SetHeader(Form("%d, |#eta| < %1.1f",year,etaCut));
      if(i==1) leg->SetHeader(Form("%d, |#eta| < %1.1f",year,etaCut));
      for(int k=0; k<nVz; k++)
	{
	  hTrkPtMthInVzAll[i][k] = (TH1F*)fin->Get(Form("hTrkPtMthAll_Vz%d_Type%d",k,i));
	  hTrkPtMthInVzAcc[i][k] = (TH1F*)fin->Get(Form("hTrkPtMthAcc_Vz%d_Type%d",k,i));
	  hTrkPtMthInVzEff[i][k] = (TH1F*)hTrkPtMthInVzAcc[i][k]->Clone(Form("hTrkPtMthEff_Vz%d_Type%d",k,i));
	  hTrkPtMthInVzEff[i][k]->Divide(hTrkPtMthInVzAll[i][k]);
	  hTrkPtMthInVzEff[i][k]->SetMarkerStyle(20+k);
	  hTrkPtMthInVzEff[i][k]->SetMarkerColor(colors[k]);
	  hTrkPtMthInVzEff[i][k]->SetLineColor(colors[k]);
	  hTrkPtMthInVzEff[i][k]->GetYaxis()->SetRangeUser(0, 0.7);
	  hTrkPtMthInVzEff[i][k]->SetTitle(";p_{T} [GeV/c];Match eff.");
	  ScaleHistoTitle(hTrkPtMthInVzEff[i][k],0.05,1,0.045,0.05,1,0.045);
	  if(k==0) hTrkPtMthInVzEff[i][k]->Draw();
	  else     hTrkPtMthInVzEff[i][k]->Draw("samesPE");
	  if(i==0) leg->AddEntry(hTrkPtMthInVzEff[i][k], Form("%1.0f < v_{z} < %1.0f cm",vzMin[k],vzMax[k]), "PL");
	  if(i==1) leg->AddEntry(hTrkPtMthInVzEff[i][k], Form("%1.0f < DCA_{z} < %1.0f cm",vzMin[k],vzMax[k]), "PL");
	}
      TPaveText *t1 = GetTitleText(Form("%s: MTD matching efficiency vs. p_{T}",legName[i].Data()),0.05);
      t1->Draw();
      leg->Draw();

      c->cd(i+3);
      SetPadMargin(gPad,0.13,0.1,0.1,0.1);
      for(int k=0; k<nVz; k++)
	{
	  hTrkPtRespInVzAll[i][k] = (TH1F*)fin->Get(Form("hTrkPtRespAll_Vz%d_Type%d",k,i));
	  hTrkPtRespInVzAcc[i][k] = (TH1F*)fin->Get(Form("hTrkPtRespAcc_Vz%d_Type%d",k,i));
	  hTrkPtRespInVzEff[i][k] = (TH1F*)hTrkPtRespInVzAcc[i][k]->Clone(Form("hTrkPtRespEff_Vz%d_Type%d",k,i));
	  hTrkPtRespInVzEff[i][k]->Divide(hTrkPtRespInVzAll[i][k]);
	  hTrkPtRespInVzEff[i][k]->SetMarkerStyle(20+k);
	  hTrkPtRespInVzEff[i][k]->SetMarkerColor(colors[k]);
	  hTrkPtRespInVzEff[i][k]->SetLineColor(colors[k]);
	  hTrkPtRespInVzEff[i][k]->GetYaxis()->SetRangeUser(0, 1);
	  hTrkPtRespInVzEff[i][k]->SetTitle(";p_{T} [GeV/c];Resp. eff.");
	  ScaleHistoTitle(hTrkPtRespInVzEff[i][k],0.05,1,0.045,0.05,1,0.045);
	  if(k==0) hTrkPtRespInVzEff[i][k]->Draw();
	  else     hTrkPtRespInVzEff[i][k]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s: MTD response efficiency vs. p_{T}",legName[i].Data()),0.05);
      t1->Draw();
    }
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthRespEffVsPt_InVzBin.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_InVzBin.pdf"));
    }

  TCanvas *c = new TCanvas("MtdRespMthEff_InEtaBin","MtdRespMthEff_InEtaBin",1100,700);
  c->Divide(2,2);
  TH1F *hTrkPtRespInEtaAll[2][nEta];
  TH1F *hTrkPtRespInEtaAcc[2][nEta];
  TH1F *hTrkPtRespInEtaEff[2][nEta];
  TH1F *hTrkPtMthInEtaAll[2][nEta];
  TH1F *hTrkPtMthInEtaAcc[2][nEta];
  TH1F *hTrkPtMthInEtaEff[2][nEta];
  for(int i=0; i<2; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad,0.13,0.1,0.1,0.1);
      TLegend *leg1 = new TLegend(0.3,0.2,0.5,0.45);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.045);
      TLegend *leg2 = new TLegend(0.6,0.2,0.8,0.45);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.045);
      if(i==0) leg1->SetHeader(Form("%d, |v_{z}| < 50 cm",year));
      if(i==1) leg1->SetHeader(Form("%d, |DCA_{z}| < 50 cm",year));
      for(int k=0; k<nEta; k++)
	{
	  hTrkPtMthInEtaAll[i][k] = (TH1F*)fin->Get(Form("hTrkPtMthAll_Eta%d_Type%d",k,i));
	  hTrkPtMthInEtaAcc[i][k] = (TH1F*)fin->Get(Form("hTrkPtMthAcc_Eta%d_Type%d",k,i));
	  hTrkPtMthInEtaEff[i][k] = (TH1F*)hTrkPtMthInEtaAcc[i][k]->Clone(Form("hTrkPtMthEff_Eta%d_Type%d",k,i));
	  hTrkPtMthInEtaEff[i][k]->Divide(hTrkPtMthInEtaAll[i][k]);
	  hTrkPtMthInEtaEff[i][k]->SetMarkerStyle(20+k);
	  hTrkPtMthInEtaEff[i][k]->SetMarkerColor(colors[k]);
	  hTrkPtMthInEtaEff[i][k]->SetLineColor(colors[k]);
	  hTrkPtMthInEtaEff[i][k]->GetYaxis()->SetRangeUser(0, 0.7);
	  hTrkPtMthInEtaEff[i][k]->SetTitle(";p_{T} [GeV/c];Match eff.");
	  ScaleHistoTitle(hTrkPtMthInEtaEff[i][k],0.05,1,0.045,0.05,1,0.045);
	  if(k==0) hTrkPtMthInEtaEff[i][k]->Draw();
	  else     hTrkPtMthInEtaEff[i][k]->Draw("samesPE");
	  if(k<nEta/2) leg1->AddEntry(hTrkPtMthInEtaEff[i][k], Form("%1.1f < #eta < %1.1f",-0.8+k*0.2,-0.6+k*0.2), "PL");
	  else         leg2->AddEntry(hTrkPtMthInEtaEff[i][k], Form("%1.1f < #eta < %1.1f",-0.8+k*0.2,-0.6+k*0.2), "PL");
	}
      TPaveText *t1 = GetTitleText(Form("%s: MTD matching efficiency vs. p_{T}",legName[i].Data()),0.05);
      t1->Draw();
      leg1->Draw();
      leg2->Draw();

      c->cd(i+3);
      SetPadMargin(gPad,0.13,0.1,0.1,0.1);
      for(int k=1; k<nEta-1; k++)
	{
	  hTrkPtRespInEtaAll[i][k] = (TH1F*)fin->Get(Form("hTrkPtRespAll_Eta%d_Type%d",k,i));
	  hTrkPtRespInEtaAcc[i][k] = (TH1F*)fin->Get(Form("hTrkPtRespAcc_Eta%d_Type%d",k,i));
	  hTrkPtRespInEtaEff[i][k] = (TH1F*)hTrkPtRespInEtaAcc[i][k]->Clone(Form("hTrkPtRespEff_Eta%d_Type%d",k,i));
	  hTrkPtRespInEtaEff[i][k]->Divide(hTrkPtRespInEtaAll[i][k]);
	  hTrkPtRespInEtaEff[i][k]->SetMarkerStyle(20+k);
	  hTrkPtRespInEtaEff[i][k]->SetMarkerColor(colors[k]);
	  hTrkPtRespInEtaEff[i][k]->SetLineColor(colors[k]);
	  hTrkPtRespInEtaEff[i][k]->GetYaxis()->SetRangeUser(0.4, 1);
	  hTrkPtRespInEtaEff[i][k]->SetTitle(";p_{T} [GeV/c];Resp. eff.");
	  ScaleHistoTitle(hTrkPtRespInEtaEff[i][k],0.05,1,0.045,0.05,1,0.045);
	  if(k==1) hTrkPtRespInEtaEff[i][k]->Draw();
	  else     hTrkPtRespInEtaEff[i][k]->Draw("samesPE");
	}
      TPaveText *t1 = GetTitleText(Form("%s: MTD response efficiency vs. p_{T}",legName[i].Data()),0.05);
      t1->Draw();
    }
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthRespEffVsPt_InEtaBin.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_InEtaBin.pdf"));
    }

  // check phiDiffCut
  TCanvas *c = new TCanvas("MtdRespMthEff_InPhiCut","MtdRespMthEff_InPhiCut",1100,700);
  c->Divide(2,2);
  TH1F *hTrkPtRespInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtRespInPhiDiffAcc[nPhiDiffCut];
  TH1F *hTrkPtRespInPhiDiffEff[nPhiDiffCut];
  TH1F *hTrkPtMthInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtMthInPhiDiffAcc[nPhiDiffCut];
  TH1F *hTrkPtMthInPhiDiffEff[nPhiDiffCut];
  c->cd(1);
  SetPadMargin(gPad,0.13,0.1,0.1,0.1);
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader(Form("%d, |#eta|<%1.1f, |DCA_{z}|<%1.0f cm",year,etaCut,vzCut));
  for(int k=0; k<nPhiDiffCut; k++)
    {
      hTrkPtRespInPhiDiffAll[k] = (TH1F*)fin->Get(Form("hTrkPtRespAll_PhiDiff%d_Type%d",k,1));
      hTrkPtRespInPhiDiffAcc[k] = (TH1F*)fin->Get(Form("hTrkPtRespAcc_PhiDiff%d_Type%d",k,1));
      hTrkPtRespInPhiDiffEff[k] = (TH1F*)hTrkPtRespInPhiDiffAcc[k]->Clone(Form("hTrkPtRespEff_PhiDiff%d_Type%d",k,1));
      hTrkPtRespInPhiDiffEff[k]->Divide(hTrkPtRespInPhiDiffAll[k]);
      hTrkPtRespInPhiDiffEff[k]->SetMarkerStyle(20+k);
      hTrkPtRespInPhiDiffEff[k]->SetMarkerColor(colors[k]);
      hTrkPtRespInPhiDiffEff[k]->SetLineColor(colors[k]);
      hTrkPtRespInPhiDiffEff[k]->GetYaxis()->SetRangeUser(0, 0.7);
      hTrkPtRespInPhiDiffEff[k]->SetTitle(";p_{T} [GeV/c];Resp. eff.");
      ScaleHistoTitle(hTrkPtRespInPhiDiffEff[k],0.045,1,0.04,0.045,1,0.04);
      hTrkPtRespInPhiDiffEff[k]->GetYaxis()->SetRangeUser(0,1);
      if(k==0) hTrkPtRespInPhiDiffEff[k]->Draw();
      else     hTrkPtRespInPhiDiffEff[k]->Draw("samesPE");
      if(k<nPhiDiffCut-1) leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("|#varphi_{mom}-#varphi_{pos}|<%4.3f",phiDiffCuts[k]), "PL");
      else                leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("no |#varphi_{mom}-#varphi_{pos}| cut"), "PL");
    }
  TPaveText *t1 = GetTitleText(Form("%s: MTD response efficiency vs. p_{T}",legName[1].Data()),0.045);
  t1->Draw();
  leg->Draw();

  c->cd(2);
  SetPadMargin(gPad,0.13,0.1,0.1,0.1);
  for(int k=0; k<nPhiDiffCut; k++)
    {
      hTrkPtMthInPhiDiffAll[k] = (TH1F*)fin->Get(Form("hTrkPtMthAll_PhiDiff%d_Type%d",k,1));
      hTrkPtMthInPhiDiffAcc[k] = (TH1F*)fin->Get(Form("hTrkPtMthAcc_PhiDiff%d_Type%d",k,1));
      hTrkPtMthInPhiDiffEff[k] = (TH1F*)hTrkPtMthInPhiDiffAcc[k]->Clone(Form("hTrkPtMthEff_PhiDiff%d_Type%d",k,1));
      hTrkPtMthInPhiDiffEff[k]->Divide(hTrkPtMthInPhiDiffAll[k]);
      hTrkPtMthInPhiDiffEff[k]->SetMarkerStyle(20+k);
      hTrkPtMthInPhiDiffEff[k]->SetMarkerColor(colors[k]);
      hTrkPtMthInPhiDiffEff[k]->SetLineColor(colors[k]);
      hTrkPtMthInPhiDiffEff[k]->GetYaxis()->SetRangeUser(0, 0.7);
      hTrkPtMthInPhiDiffEff[k]->SetTitle(";p_{T} [GeV/c];Match eff.");
      ScaleHistoTitle(hTrkPtRespInPhiDiffEff[k],0.045,1,0.04,0.045,1,0.04);
      if(k==0) hTrkPtMthInPhiDiffEff[k]->Draw();
      else     hTrkPtMthInPhiDiffEff[k]->Draw("samesPE");
    }
  TPaveText *t1 = GetTitleText(Form("%s: MTD matching efficiency vs. p_{T}",legName[1].Data()),0.045);
  t1->Draw();
  leg->Draw();

  c->cd(3);
  for(int k=0; k<nPhiDiffCut-1; k++)
    {
      TH1F *hratio = (TH1F*)hTrkPtRespInPhiDiffEff[k]->Clone(Form("%s_ratio",hTrkPtRespInPhiDiffEff[k]->GetName()));
      hratio->Divide(hTrkPtRespInPhiDiffEff[nPhiDiffCut-1]);
      hratio->GetYaxis()->SetRangeUser(0.8,1.2);
      hratio->SetYTitle("Ratio to inclusive");
      if(k==0) hratio->Draw();
      else     hratio->Draw("samesPE");
    }
  TPaveText *t1 = GetTitleText(Form("%s: ratio of MTD response efficiency",legName[1].Data()),0.045);
  t1->Draw();

  c->cd(4);
  for(int k=0; k<nPhiDiffCut-1; k++)
    {
      TH1F *hratio = (TH1F*)hTrkPtMthInPhiDiffEff[k]->Clone(Form("%s_ratio",hTrkPtMthInPhiDiffEff[k]->GetName()));
      hratio->Divide(hTrkPtMthInPhiDiffEff[nPhiDiffCut-1]);
      hratio->GetYaxis()->SetRangeUser(0.5,1.2);
      hratio->SetYTitle("Ratio to inclusive");
      if(k==0) hratio->Draw();
      else     hratio->Draw("samesPE");
    }
  TPaveText *t1 = GetTitleText(Form("%s: ratio of MTD matching efficiency",legName[1].Data()),0.045);
  t1->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthRespEffVsPt_InPhiDiffCutBin.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_InPhiDiffBin.pdf"));
    }

  TH2F *hCellVsPhiDiff = (TH2F*)fin->Get("hCellVsPhiDiff_BtmBL_Type1");
  c = draw2D(hCellVsPhiDiff,Form("%d: correlation between cell and #Delta#varphi (BL10-22, 2 < p_{T} < 4 GeV/c)",year));
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/PhiDiffVsCell_BtmBL.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_PhiDiffVsCell_BtmBL.pdf"));
    }

  TH1F *hTrkPtVsCellInPhiDiff[nPhiDiffCut];
  TCanvas *c2 = new TCanvas("cTrkPtVsCellInPhiDiff","cTrkPtVsCellInPhiDiff",800, 600);
  TLegend *leg = new TLegend(0.2,0.55,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader(Form("|#eta|<%1.1f, |DCA_{z}|<%1.0f cm, 2 < p_{T} < 4 GeV/c",etaCut,vzCut));
  for(int k=nPhiDiffCut-1; k>-1; k--)
    {
      hTrkPtVsCellInPhiDiff[k] = (TH1F*)fin->Get(Form("hTrkPtVsCell_BtmBL_PhiDiff%d_Type1",k));
      hTrkPtVsCellInPhiDiff[k]->Sumw2();
      hTrkPtVsCellInPhiDiff[k]->SetMarkerStyle(20+k);
      hTrkPtVsCellInPhiDiff[k]->SetMarkerColor(colors[k]);
      hTrkPtVsCellInPhiDiff[k]->SetLineColor(colors[k]);
      hTrkPtVsCellInPhiDiff[k]->Scale(1./hTrkPtVsCellInPhiDiff[k]->Integral());
      hTrkPtVsCellInPhiDiff[k]->GetYaxis()->SetRangeUser(0.03, 0.11);
      hTrkPtVsCellInPhiDiff[k]->SetTitle();
      if(k==nPhiDiffCut-1) 
	{
	  hTrkPtVsCellInPhiDiff[k]->Draw("P");
	  leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("no |#varphi_{mom}-#varphi_{pos}| cut"), "PL");
	}
      else     
	{
	  hTrkPtVsCellInPhiDiff[k]->Draw("samesP");
	  leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("|#varphi_{mom}-#varphi_{pos}|<%4.3f",phiDiffCuts[k]), "PL");
	}         
    }
  leg->Draw();
  TPaveText *t1 = GetTitleText(Form("%d: cell distribution in BL10-22",year),0.04);
  t1->Draw();
  if(savePlot) 
    c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/PhiDiff_CompCell_BtmBL.pdf",run_type));
  if(gSaveAN)
    {
      c2->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_CompCell_BtmBL.pdf"));
    }


  TCanvas *c3 = new TCanvas("cTrkPtEffInCell","cTrkPtEffInCell",1100, 600);
  c3->Divide(3,2);
  TH1F *hTrkPtInCellInPhiDiffAll[nPhiDiffCut][nCellBin];
  TH1F *hTrkPtInCellInPhiDiffAcc[nPhiDiffCut][nCellBin];
  TH1F *hTrkPtInCellInPhiDiffEff[nPhiDiffCut][nCellBin];
  TLegend *leg = new TLegend(0.2,0.4,0.6,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(Form("%d: |#eta|<%1.1f, |DCA_{z}|<%1.0f cm",year,etaCut,vzCut));
  for(int j=0; j<nCellBin; j++)
    {
      c3->cd(j+1);
      for(int k=0; k<nPhiDiffCut; k++)
	{
	  hTrkPtInCellInPhiDiffAll[k][j] = (TH1F*)fin->Get(Form("hTrkPtAll_BtmBL_Cell%d_PhiDiff%d_Type1",j,k));
	  hTrkPtInCellInPhiDiffAcc[k][j] = (TH1F*)fin->Get(Form("hTrkPtAcc_BtmBL_Cell%d_PhiDiff%d_Type1",j,k));
	  hTrkPtInCellInPhiDiffEff[k][j] = (TH1F*)hTrkPtInCellInPhiDiffAcc[k][j]->Clone(Form("hTrkPtEff_BtmBL_Cell%d_PhiDiff%d_Type1",j,k));
	  hTrkPtInCellInPhiDiffEff[k][j]->Divide(hTrkPtInCellInPhiDiffAll[k][j]);
	  hTrkPtInCellInPhiDiffEff[k][j]->SetMarkerStyle(20+k);
	  hTrkPtInCellInPhiDiffEff[k][j]->SetMarkerColor(colors[k]);
	  hTrkPtInCellInPhiDiffEff[k][j]->SetLineColor(colors[k]);
	  if(j!=2) hTrkPtInCellInPhiDiffEff[k][j]->SetMaximum(1.5*hTrkPtInCellInPhiDiffEff[k][j]->GetMaximum());
	  hTrkPtInCellInPhiDiffEff[k][j]->SetTitle("");
	  if(k==0) hTrkPtInCellInPhiDiffEff[k][j]->Draw();
	  else     hTrkPtInCellInPhiDiffEff[k][j]->Draw("sames");
	  if(j==0)
	    {
	      if(k==nPhiDiffCut-1) 
		{
		  leg->AddEntry(hTrkPtInCellInPhiDiffEff[k][j], Form("no |#varphi_{mom}-#varphi_{pos}| cut"), "PL");
		}
	      else     
		{
		  leg->AddEntry(hTrkPtInCellInPhiDiffEff[k][j], Form("|#varphi_{mom}-#varphi_{pos}|<%4.3f",phiDiffCuts[k]), "PL");
		}
	    }
	}
      TPaveText *t1 = GetTitleText(Form("BL 10-22, Cell: %d -> %d",start_cell[j]-4,end_cell[j]-4),0.06);
      t1->Draw();
    }
  c3->cd(6);
  leg->Draw();
  if(savePlot) 
    c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/PhiDiff_CompEffInCell_BtmBL.pdf",run_type));
  if(gSaveAN)
    {
      c3->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_CompEffInCell_BtmBL.pdf"));
    }

  TH1F *hTrkPtInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtInPhiDiffAcc[nPhiDiffCut];
  TH1F *hTrkPtInPhiDiffEff[nPhiDiffCut];
  TCanvas *c4 = new TCanvas("cTrkPtEffInCellWeighted","cTrkPtEffInCellWeighted",1100, 500);
  c4->Divide(2,1);
  c4->cd(1);
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.55);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("|#eta|<%1.1f, |DCA_{z}|<%1.0f cm",etaCut,vzCut));
  for(int k=nPhiDiffCut-1; k>-1; k--)
    {
      hTrkPtInPhiDiffAll[k] = (TH1F*)fin->Get(Form("hTrkPtAll_BtmBL_PhiDiff%d_Type1",k));
      TH1F *hAll = (TH1F*)hTrkPtInPhiDiffAll[k]->Rebin(nPtBins, Form("%s_rebin",hTrkPtInPhiDiffAll[k]->GetName()), xPtBins);
      hTrkPtInPhiDiffAcc[k] = (TH1F*)fin->Get(Form("hTrkPtAcc_BtmBL_PhiDiff%d_Type1",k));
      hTrkPtInPhiDiffEff[k] = (TH1F*)hTrkPtInPhiDiffAcc[k]->Rebin(nPtBins, Form("hTrkPtEff_BtmBL_PhiDiff%d_Type1",k), xPtBins);
      hTrkPtInPhiDiffEff[k]->Divide(hAll);
      hTrkPtInPhiDiffEff[k]->SetMarkerStyle(20+k);
      hTrkPtInPhiDiffEff[k]->SetMarkerColor(colors[k]);
      hTrkPtInPhiDiffEff[k]->SetLineColor(colors[k]);
      hTrkPtInPhiDiffEff[k]->SetTitle("");
      if(k==nPhiDiffCut-1) 
	{
	  hTrkPtInPhiDiffEff[k]->Draw();
	  leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("no |#varphi_{mom}-#varphi_{pos}| cut"), "PL");
	}
      else     
	{
	  hTrkPtInPhiDiffEff[k]->Draw("sames");
	  leg->AddEntry(hTrkPtRespInPhiDiffEff[k], Form("|#varphi_{mom}-#varphi_{pos}|<%4.3f",phiDiffCuts[k]), "PL");
	}
    }
  TPaveText *t1 = GetTitleText(Form("Cosmic ray: matching efficiency for BL10-22"),0.04);
  t1->Draw();
  leg->Draw();
  c4->cd(2);
  for(int k=0; k<nPhiDiffCut; k++)
    {
      TH1F *hRatio = (TH1F*)hTrkPtInPhiDiffEff[k]->Clone(Form("hRatio_PhiDiff%d",k));
      hRatio->Divide(hTrkPtInPhiDiffEff[nPhiDiffCut-3]);
      hRatio->SetYTitle("ratio");
      if(k==0) hRatio->Draw();
      else     hRatio->Draw("sames");
    }
  if(savePlot) 
    c4->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/PhiDiff_CompEff_BtmBL.pdf",run_type));
  if(gSaveAN)
    {
      c4->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_CompEff_BtmBL.pdf"));
    }

  // compare for track phi and eta
  TH1F *hTrkPhiAll[2];
  TH1F *hTrkPhiAcc[2];
  TH1F *hTrkPhiEff[2];
  TH1F *hTrkEtaAll[2];
  TH1F *hTrkEtaAcc[2];
  TH1F *hTrkEtaEff[2];
  TH1F *hTrkDCAzAll[2];
  TH1F *hTrkDCAzAcc[2];
  TH1F *hTrkDCAzEff[2];
  for(int i=0; i<2; i++)
    {
      hTrkPhiAll[i] = (TH1F*)fin->Get(Form("hTrkPhiAll_pt%1.0f_Type%d",ptCut,i));
      hTrkPhiAll[i]->SetMarkerStyle(20+i*4);
      hTrkPhiAll[i]->SetMarkerColor(i+1);
      hTrkPhiAll[i]->SetLineColor(i+1);
      hTrkPhiAcc[i] = (TH1F*)fin->Get(Form("hTrkPhiAcc_pt%1.0f_Type%d",ptCut,i));
      hTrkPhiEff[i] = (TH1F*)hTrkPhiAcc[i]->Clone(Form("hTrkPhiEff_pt%1.0f_Type%d",ptCut,i));
      hTrkPhiEff[i]->Divide(hTrkPhiAll[i]);
      hTrkPhiEff[i]->SetMarkerStyle(20+i*4);
      hTrkPhiEff[i]->SetMarkerColor(i+1);
      hTrkPhiEff[i]->SetLineColor(i+1);
      hTrkPhiEff[i]->GetYaxis()->SetRangeUser(0,1);

      hTrkEtaAll[i] = (TH1F*)fin->Get(Form("hTrkEtaAll_pt%1.0f_Type%d",ptCut,i));
      hTrkEtaAll[i]->SetMarkerStyle(20+i*4);
      hTrkEtaAll[i]->SetMarkerColor(i+1);
      hTrkEtaAll[i]->SetLineColor(i+1);
      hTrkEtaAcc[i] = (TH1F*)fin->Get(Form("hTrkEtaAcc_pt%1.0f_Type%d",ptCut,i));
      hTrkEtaEff[i] = (TH1F*)hTrkEtaAcc[i]->Clone(Form("hTrkEtaEff_pt%1.0f_Type%d",ptCut,i));
      hTrkEtaEff[i]->Divide(hTrkEtaAll[i]);
      hTrkEtaEff[i]->SetMarkerStyle(20+i*4);
      hTrkEtaEff[i]->SetMarkerColor(i+1);
      hTrkEtaEff[i]->SetLineColor(i+1);

      hTrkDCAzAll[i] = (TH1F*)fin->Get(Form("hTrkDCAzAll_pt%1.0f_Type%d",ptCut,i));
      hTrkDCAzAll[i]->SetMarkerStyle(20+i*4);
      hTrkDCAzAll[i]->SetMarkerColor(i+1);
      hTrkDCAzAll[i]->SetLineColor(i+1);
      hTrkDCAzAcc[i] = (TH1F*)fin->Get(Form("hTrkDCAzAcc_pt%1.0f_Type%d",ptCut,i));
      hTrkDCAzEff[i] = (TH1F*)hTrkDCAzAcc[i]->Clone(Form("hTrkDCAzEff_pt%1.0f_Type%d",ptCut,i));
      hTrkDCAzEff[i]->Divide(hTrkDCAzAll[i]);
      hTrkDCAzEff[i]->SetMarkerStyle(20+i*4);
      hTrkDCAzEff[i]->SetMarkerColor(i+1);
      hTrkDCAzEff[i]->SetLineColor(i+1);
    }

  TCanvas *cPhi = new TCanvas("cMthVsPhi", "cMthVsPhi", 1200, 700);
  cPhi->Divide(3, 2);
  for(int i=0; i<2; i++)
    {
      cPhi->cd(1);
      hTrkDCAzAll[i]->Scale(1./hTrkDCAzAll[i]->Integral());
      hTrkDCAzAll[i]->SetTitle(";DCA_{z} (v_{z}) [cm];");
      hTrkDCAzAll[i]->SetMaximum(1.2*hTrkDCAzAll[i]->GetMaximum());
      if(i==0)
	{
	  hTrkDCAzAll[i]->Draw();
	  t1 = GetTitleText(Form("DCA_{z} (v_{z}) distribution of projected tracks"),0.05);
	  t1->Draw();
	}
      else     hTrkDCAzAll[i]->Draw("sames");

      cPhi->cd(2);
      hTrkPhiAll[i]->Scale(1./hTrkPhiAll[i]->Integral());
      hTrkPhiAll[i]->SetTitle(";#varphi;");
      hTrkPhiAll[i]->SetMaximum(3*hTrkPhiAll[i]->GetMaximum());
      if(i==0)
	{
	  hTrkPhiAll[i]->Draw();
	  t1 = GetTitleText(Form("#varphi distribution of projected tracks"),0.05);
	  t1->Draw();
	}
      else     hTrkPhiAll[i]->Draw("sames");

      cPhi->cd(3);
      hTrkEtaAll[i]->Scale(1./hTrkEtaAll[i]->Integral());
      hTrkEtaAll[i]->SetTitle(";#eta;");
      hTrkEtaAll[i]->SetMaximum(2*hTrkEtaAll[i]->GetMaximum());
      if(i==0) 
	{
	  hTrkEtaAll[i]->Draw();
	  t1 = GetTitleText(Form(" #eta distribution of projected tracks"),0.05);
	  t1->Draw();
	}
      else     hTrkEtaAll[i]->Draw("sames");

      cPhi->cd(4);
      hTrkDCAzEff[i]->SetTitle(";DCA_{z} (v_{z}) [cm];");
      hTrkDCAzEff[i]->GetYaxis()->SetRangeUser(0.5,0.65);
      if(i==0)
	{
	  hTrkDCAzEff[i]->Draw();
	  t1 = GetTitleText(Form("Matching efficiency vs. DCA_{z} (v_{z})"),0.05);
	  t1->Draw();
	}
      else     hTrkDCAzEff[i]->Draw("sames");

      cPhi->cd(5);
      hTrkPhiEff[i]->SetTitle(";#varphi;");
      hTrkPhiEff[i]->GetYaxis()->SetRangeUser(0,1.2);
      if(i==0)
	{
	  hTrkPhiEff[i]->Draw();
	  t1 = GetTitleText(Form("Matching efficiency vs. #varphi"),0.05);
	  t1->Draw();
	}
      else     hTrkPhiEff[i]->Draw("sames");

      cPhi->cd(6);
      hTrkEtaEff[i]->SetTitle(";#eta;");
      hTrkEtaEff[i]->GetYaxis()->SetRangeUser(0,0.8);
      if(i==0)
	{
	  hTrkEtaEff[i]->Draw();
	  t1 = GetTitleText(Form("Matching efficiency vs. #eta"),0.05);
	  t1->Draw();
	}
      else     hTrkEtaEff[i]->Draw("sames");
    }
  if(savePlot) 
    cPhi->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthEff_TrkDCAzPhiEtaDis.pdf",run_type));
  
  // compare matching efficiency for a given vz range
  TH1F *hTrkPtInVzAll[2][nVz];
  TH1F *hTrkPtInVzAcc[2][nVz];
  TH1F *hTrkPtInVzEff[2][nVz];
  TCanvas *c = new TCanvas("MthEffVsPtInVz","MthEffVsPtInVz",1100,700);
  c->Divide(3,2);
  for(int k=0; k<nVz; k++)
    {
      c->cd(k+1);
      for(int i=0; i<2; i++)
	{
	  hTrkPtInVzAll[i][k] = (TH1F*)fin->Get(Form("hTrkPtInVzAll_Vz%d_Type%d",k,i));
	  hTrkPtInVzAcc[i][k] = (TH1F*)fin->Get(Form("hTrkPtInVzAcc_Vz%d_Type%d",k,i));
	  hTrkPtInVzEff[i][k] = (TH1F*)hTrkPtInVzAcc[i][k]->Clone(Form("hTrkPtInVzEff_Vz%d_Type%d",k,i));
	  hTrkPtInVzEff[i][k]->Divide(hTrkPtInVzAll[i][k]);
	  hTrkPtInVzEff[i][k]->GetXaxis()->SetRangeUser(0,15);
	  hTrkPtInVzEff[i][k]->GetYaxis()->SetRangeUser(0,0.8);
	  hTrkPtInVzEff[i][k]->SetTitle("");
	  hTrkPtInVzEff[i][k]->SetMarkerStyle(20+i*4);
	  hTrkPtInVzEff[i][k]->SetMarkerColor(i+1);
	  hTrkPtInVzEff[i][k]->SetLineColor(i+1);
	  if(i==0) hTrkPtInVzEff[i][k]->Draw("PE");
	  else     hTrkPtInVzEff[i][k]->Draw("samesPE");

	  TPaveText *t1 = GetTitleText(Form("MTD matching efficiency"), 0.045);
	  t1->Draw();
	  t1 = GetPaveText(0.6,0.8,0.2,0.4,0.05);
	  t1->AddText(Form("%1.0f < v_{z} < %1.0f cm",vzMin[k],vzMax[k]));
	  t1->AddText(Form("|#eta| < %1.1f",etaCut));
	  t1->Draw();
	}
    }
  c->cd(6);
  TLegend *leg = new TLegend(0.15,0.5,0.5,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(Form("Run%d",year-2000));
  for(int i=0; i<2; i++)
    {
      leg->AddEntry(hTrkPtInVzEff[i][0], legName[i].Data(), "PL");
    }
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CompMthEffVsPt_InVzBin.pdf",run_type));

  // compare eff vs. cell
  TH1F *hTrkPtVsCellAll[2][30]; // 0 - embed; 1 - cosmic; 2 - weighted cosmic
  TH1F *hTrkPtVsCellAcc[2][30];
  TH1F *hTrkPtVsCellEff[2][30];
  TCanvas *c = new TCanvas("TrkPtVsCellInBL","TrkPtVsCellInBL",1100,700);
  c->Divide(6,5);
  for(int j=0; j<30; j++)
    {
      c->cd(j+1);
      for(int i=0; i<2; i++)
	{
	  hTrkPtVsCellAll[i][j] = (TH1F*)fin->Get(Form("hTrkPtVsCellAll_BL%d_Type%d",j+1,i));
	  hTrkPtVsCellAcc[i][j] = (TH1F*)fin->Get(Form("hTrkPtVsCellAcc_BL%d_Type%d",j+1,i));
	  hTrkPtVsCellEff[i][j] = (TH1F*)hTrkPtVsCellAcc[i][j]->Clone(Form("hTrkPtVsCellEff_BL%d_Type%d",j+1,i));
	  hTrkPtVsCellEff[i][j]->Divide(hTrkPtVsCellAll[i][j]);
	  if(j==8 || j==22) continue;
	  hTrkPtVsCellAll[i][j]->Scale(1./hTrkPtVsCellAll[i][j]->Integral());
	  hTrkPtVsCellAll[i][j]->SetTitle("");
	  hTrkPtVsCellAll[i][j]->SetMarkerStyle(20+i*4);
	  hTrkPtVsCellAll[i][j]->SetMarkerColor(i+1);
	  hTrkPtVsCellAll[i][j]->SetLineColor(i+1);
	  hTrkPtVsCellAll[i][j]->GetYaxis()->SetRangeUser(0,0.2);
	  if(i==0) hTrkPtVsCellAll[i][j]->Draw();
	  else     hTrkPtVsCellAll[i][j]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.1);
      t1->Draw();
    }
  c->cd(1);
  t1 = GetPaveText(0.45,0.65,0.45,0.8,0.1);
  t1->AddText(Form("p_{T} > %1.0f GeV/c",ptCut));
  t1->AddText(Form("|v_{z}| < %1.0f cm",vzCut));
  t1->AddText(Form("|#eta| < %1.1f",etaCut));
  t1->SetTextAlign(11);
  t1->Draw();
  c->cd(9);
  TLegend *leg = new TLegend(0.15,0.15,0.7,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.1);
  leg->SetHeader(Form("Run%d",year-2000));
  leg->AddEntry(hTrkPtVsCellAll[0][0], "Embed", "P");
  leg->AddEntry(hTrkPtVsCellAll[1][0], "Cosmic ray", "P");
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CompTrkPtVsCellInBL.pdf",run_type));

  TCanvas *c = new TCanvas("MthEffVsCellInBL","MthEffVsCellInBL",1100,700);
  c->Divide(6,5);
  for(int j=0; j<30; j++)
    {
      c->cd(j+1);
      for(int i=0; i<2; i++)
	{
	  hTrkPtVsCellEff[i][j]->SetTitle("");
	  hTrkPtVsCellEff[i][j]->SetMarkerStyle(20+i*4);
	  hTrkPtVsCellEff[i][j]->SetMarkerColor(i+1);
	  hTrkPtVsCellEff[i][j]->SetLineColor(i+1);
	  hTrkPtVsCellEff[i][j]->GetYaxis()->SetRangeUser(0,1.2);
	  if(i==0) hTrkPtVsCellEff[i][j]->Draw();
	  else     hTrkPtVsCellEff[i][j]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.1);
      t1->Draw();
    }
  c->cd(1);
  t1 = GetPaveText(0.3,0.6,0.1,0.4,0.1);
  t1->AddText(Form("p_{T} > %1.0f GeV/c",ptCut));
  t1->AddText(Form("|v_{z}| < %1.0f cm",vzCut));
  t1->AddText(Form("|#eta| < %1.1f",etaCut));
  t1->SetTextAlign(11);
  t1->Draw();
  c->cd(9);
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CompMthEffVsCellInBL.pdf",run_type));

  // Compare for each BL
  TH1F *hTrkPtInBlAll[3][30];
  TH1F *hTrkPtInBlAcc[3][30];
  TH1F *hTrkPtInBlEff[3][30];
  TCanvas *c = new TCanvas("MthEffVsPtInBL","MthEffVsPtInBL",1100,700);
  c->Divide(6,5);
  for(int j=0; j<30; j++)
    {
      c->cd(j+1);
      for(int i=0; i<3; i++)
	{
	  if(i<2)
	    {
	      hTrkPtInBlAll[i][j] = (TH1F*)fin->Get(Form("hTrkPtInBlAll_BL%d_Type%d",j+1,i));
	      hTrkPtInBlAcc[i][j] = (TH1F*)fin->Get(Form("hTrkPtInBlAcc_BL%d_Type%d",j+1,i));
	      int bin_index = hTrkPtInBlAll[i][j]->FindFixBin(ptCut);
	    }
	  else
	    {
	      hTrkPtInBlAll[i][j] = (TH1F*)fin->Get(Form("hTrkPtInBlAll_BL%d_Type1_Weight",j+1));
	      hTrkPtInBlAcc[i][j] = (TH1F*)fin->Get(Form("hTrkPtInBlAcc_BL%d_Type1_Weight",j+1));
	    }
	  hTrkPtInBlEff[i][j] = (TH1F*)hTrkPtInBlAcc[i][j]->Clone(Form("hTrkPtInBlEff_BL%d_Type%d",j+1,i));
	  hTrkPtInBlEff[i][j]->Divide(hTrkPtInBlAll[i][j]);
	  hTrkPtInBlEff[i][j]->SetTitle("");
	  hTrkPtInBlEff[i][j]->SetMarkerStyle(20+i*4);
	  hTrkPtInBlEff[i][j]->SetMarkerColor(pow(2,i));
	  hTrkPtInBlEff[i][j]->SetLineColor(pow(2,i));
	  hTrkPtInBlEff[i][j]->GetYaxis()->SetRangeUser(0,1);
	  if(i==0) hTrkPtInBlEff[i][j]->Draw();
	  else     hTrkPtInBlEff[i][j]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.1);
      t1->Draw();
    }
  c->cd(1);
  t1 = GetPaveText(0.5,0.7,0.15,0.4,0.1);
  t1->SetTextAlign(11);
  t1->AddText(Form("|v_{z}| < %1.0f cm",vzCut));
  t1->AddText(Form("|#eta| < %1.1f",etaCut));
  t1->SetTextColor(6);
  t1->Draw();
  c->cd(9);
  TLegend *leg = new TLegend(0.15,0.15,0.7,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.08);
  leg->SetHeader(Form("Run%d",year-2000));
  leg->AddEntry(hTrkPtInBlEff[0][0], "Embed", "P");
  leg->AddEntry(hTrkPtInBlEff[1][0], "Cosmic ray", "P");
  leg->AddEntry(hTrkPtInBlEff[2][0], "Cosmic ray (weighted)", "P");
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CompMthEffVsPtInBL.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_CompEffInBL.pdf"));
    }

  TH1F *hTrkPtInBlEffRatio[30];
  TCanvas *c = new TCanvas("MthEffVsPtInBLRatio","MthEffVsPtInBLRatio",1100,700);
  c->Divide(6,5);
  for(int j=0; j<30; j++)
    {
      c->cd(j+1);
      gPad->SetGrid(0,1);
      hTrkPtInBlEffRatio[j] = (TH1F*)hTrkPtInBlEff[0][j]->Clone(Form("hTrkPtInBlEffRatio_BL%d",j+1));
      hTrkPtInBlEffRatio[j]->Divide((TH1F*)hTrkPtInBlEff[2][j]);
      hTrkPtInBlEffRatio[j]->GetYaxis()->SetRangeUser(0.5,1.3);
      hTrkPtInBlEffRatio[j]->Draw();
      TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.1);
      t1->Draw();
      if(j==8 || j==22) continue;
      TLine *line = GetLine(0,1,15,1);
      line->Draw();
    }
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CompMthEffVsPtInBLRatio.pdf",run_type));

  // compare efficiency vs. pt in different cells
  TH1F *hTrkPtInCellAll[2][30][5];
  TH1F *hTrkPtInCellAcc[2][30][5];
  TH1F *hTrkPtInCellEff[2][30][5];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<30; j++)
	{
	  double nAll = 0, nAcc = 0;
	  for(int k=0; k<5; k++)
	    {
	      hTrkPtInCellAll[i][j][k] = (TH1F*)fin->Get(Form("hTrkPtAll_BL%d_Cell%d_Type%d",j+1,k,i));
	      hTrkPtInCellAcc[i][j][k] = (TH1F*)fin->Get(Form("hTrkPtAcc_BL%d_Cell%d_Type%d",j+1,k,i));
	      hTrkPtInCellEff[i][j][k] = (TH1F*)hTrkPtInCellAcc[i][j][k]->Clone(Form("hTrkPtEff_BL%d_Cell%d_Type%d",j+1,k,i));
	      hTrkPtInCellEff[i][j][k]->Divide(hTrkPtInCellAll[i][j][k]);
	    }
	}
    }
  for(int j=9; j<22; j++)
    {
      int bl = j+1;
      TCanvas *c = new TCanvas(Form("MthEffVsPt_InCell_BL%d",bl),Form("MthEffVsPt_InCell_BL%d",bl),900,500);
      c->Divide(3,2);
      for(int k=0; k<5; k++)
	{
	  c->cd(k+1);
	  for(int i=0; i<2; i++)
	    {
	      hTrkPtInCellEff[i][j][k]->SetMarkerStyle(20+i*4);
	      hTrkPtInCellEff[i][j][k]->SetMarkerColor(i+1);
	      hTrkPtInCellEff[i][j][k]->SetLineColor(i+1);
	      hTrkPtInCellEff[i][j][k]->GetYaxis()->SetRangeUser(0,1);
	      hTrkPtInCellEff[i][j][k]->SetTitle(";p_{T} [GeV/c];Match Eff.");
	      if(i==0) hTrkPtInCellEff[i][j][k]->Draw("PE");
	      else     hTrkPtInCellEff[i][j][k]->Draw("samesPE");
	    }
	  TPaveText *t1 = GetTitleText(Form("BL = %d, Cell: %d -> %d",bl,start_cell[k]-4,end_cell[k]-4),0.06);
	  t1->Draw();
	}
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthEffVsPt_InCell_BL%d.pdf",run_type,bl));
      if(bl==11 || bl==17)
	{
	  if(gSaveAN)
	    {
	      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_CompInCell_BL%d.pdf",bl));
	    }
	}
    }

  // avergae difference between cosmic ray and embedding
  TH1F *hTrkPtBtmBlAll[2];
  TH1F *hTrkPtBtmBlAcc[2];
  TH1F *hTrkPtBtmBlEff[2];
  TH1F *hWeightEmbedBL = new TH1F("hWeightEmbedBL","hWeightEmbedBL",13,10,22);
  for(int j=9; j<=21; j++)
    {
      hWeightEmbedBL->SetBinContent(j-8, hTrkPtInBlAll[0][j]->Integral());
    }
  hWeightEmbedBL->Scale(1./hWeightEmbedBL->Integral());

  for(int i=0; i<2; i++)
    {
      hTrkPtBtmBlAll[i] = (TH1F*)hTrkPtInBlAll[i*2][0]->Clone(Form("hTrkPtBtmBlAll_Type%d",i));
      hTrkPtBtmBlAll[i]->Reset("AC");
      hTrkPtBtmBlAcc[i] = (TH1F*)hTrkPtInBlAcc[i*2][0]->Clone(Form("hTrkPtBtmBlAcc_Type%d",i));
      hTrkPtBtmBlAcc[i]->Reset("AC");

      for(int j=9; j<=21; j++)
	{
	  double weight = 1;
	  if(i==1) weight = hWeightEmbedBL->GetBinContent(j-8)/hTrkPtInBlAll[i*2][j]->Integral();
	  hTrkPtBtmBlAll[i]->Add(hTrkPtInBlAll[i*2][j], weight);
	  hTrkPtBtmBlAcc[i]->Add(hTrkPtInBlAcc[i*2][j], weight); 
	}
      hTrkPtBtmBlEff[i] = (TH1F*)hTrkPtBtmBlAcc[i]->Clone(Form("hTrkPtBtmBlEff_Type%d",i));
      hTrkPtBtmBlEff[i]->Divide(hTrkPtBtmBlAll[i]);
    }

  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%d",year));
  TH1F *hMtdMthEffComp[2];
  TF1 *funcTrkPtBtmBlEff[2];
  for(int i=0; i<2; i++)
    {
      hMtdMthEffComp[i] = (TH1F*)hTrkPtBtmBlEff[i]->Clone(Form("hMtdMthEffComp_Type%d",i));
      hMtdMthEffComp[i]->SetTitle("");
      hMtdMthEffComp[i]->SetMarkerStyle(20+i*4);
      hMtdMthEffComp[i]->SetMarkerColor(i+1);
      hMtdMthEffComp[i]->SetLineColor(i+1);
      hMtdMthEffComp[i]->GetYaxis()->SetRangeUser(0.2,0.8);

      //funcTrkPtBtmBlEff[i] = new TF1(Form("funcTrkPtBtmBlEff_Type%d",i),"[0]/(x-[1])+[2]/(x-[3])+[4]",1.3,20);
      funcTrkPtBtmBlEff[i] = new TF1(Form("funcTrkPtBtmBlEff_Type%d",i),"[0]/([1]*x+[2]*x*x+[5]*x*x*x-[3])+[4]",1.3,20);
      funcTrkPtBtmBlEff[i]->SetParameters(-0.5, 1.2, -0.5, 1.2, 0.8,0);
      hMtdMthEffComp[i]->Fit(funcTrkPtBtmBlEff[i], "R0QM");
      funcTrkPtBtmBlEff[i]->SetLineColor(i+1);
      funcTrkPtBtmBlEff[i]->SetLineStyle(2);
      leg->AddEntry(hMtdMthEffComp[i],legName[i].Data(),"P");
    }
  c = draw1D(hMtdMthEffComp[0],"MTD matching efficiency");
  hMtdMthEffComp[1]->Draw("sames");
  for(int i=0; i<2; i++)
    {
      funcTrkPtBtmBlEff[i]->Draw("samesHIST");
    }
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/FitMthEffVsPtInBtmBL.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_BtmBL_CosmicVsEmbed.pdf"));
    }

  TH1F *hRatio = (TH1F*)hMtdMthEffComp[0]->Clone("hratio");
  hRatio->Divide(hMtdMthEffComp[1]);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.2);
  c = draw1D(hRatio,"Embed/Cosmic: MTD matching efficiency ratio for BL9-21");
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MthEffVsPtInBtmBLRatio.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_MthEff_BtmBL_CosmicOverEmbed.pdf"));
    }
  
  if(saveHisto == 1)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type), "recreate");
      for(int i=0; i<2; i++)
	{
	  funcTrkPtBtmBlEff[i]->Write("",TObject::kOverwrite);
	}
    }

  return;
  // day-by-day
  TH1F *hDayAll[nYears];
  TH1F *hDayAcc[nYears];
  TH1F *hDayEff[nYears];
  for(int i=0; i<nYears; i++)
    {
      hDayAll[i] = (TH1F*)fin->Get(Form("hTrkDayAll_pt%1.0f_Year%d",ptCut,years[i]));
      hDayAcc[i] = (TH1F*)fin->Get(Form("hTrkDayAcc_pt%1.0f_Year%d",ptCut,years[i]));

      hDayEff[i] = (TH1F*)hDayAcc[i]->Clone(Form("hDayEff_pt%1.0f_Year%d",ptCut,years[i]));
      hDayEff[i]->Divide(hDayAll[i]);
      hDayEff[i]->GetYaxis()->SetRangeUser(0, 0.8);
      hDayEff[i]->SetMarkerStyle(20+i);
      hDayEff[i]->SetMarkerSize(1.5);
      hDayEff[i]->SetMarkerColor(pow(2,i));
      hDayEff[i]->SetLineColor(pow(2,i));
    }
  c = draw1D(hDayEff[0],Form("MTD matching efficiency (|v_{z}| < %1.0f cm, p_{T} > %1.0f GeV/c);day;Efficiency",vzCut,ptCut));
  for(int i=1; i<nYears; i++) hDayEff[i]->Draw("sames");
  TLegend *leg = new TLegend(0.6,0.2,0.8,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Cosmic ray");
  for(int i=0; i<nYears; i++)
    {
      leg->AddEntry(hDayEff[i], Form("Year %d",years[i]), "PL");
    }
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/CosmicMthEffVsDay.pdf",run_type));
}


//================================================
void makeHisto(const int saveHisto = 1)
{
  const TString legName[2] = {"Embedding","Cosmic ray"}
  TFile *fdata[2];
  fdata[0] = TFile::Open(Form("output/%s.Embed.Jpsi.root",run_type),"read");
  fdata[1] = TFile::Open(Form("output/Run%d.cosmic.root",year-2000),"read");
  THnSparseF *hnMthEff[2];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i] = (THnSparseF*)fdata[i]->Get("mhMtdMatchEff");
      hnMthEff[i]->SetName(Form("%s_Type%d",hnMthEff[i]->GetName(),i));
      printf("[i] Input file %s\n",fdata[i]->GetName());
    }

  TH1F *h1tmp = 0x0, *h2tmp = 0x0;

  // phi correlation between mom and pos
  TH2F *hPosVsMomPhi[2];
  TH1F *hPhiDiff[2];
  for(int i=0; i<2; i++)
    {
      hPosVsMomPhi[i] = (TH2F*)fdata[i]->Get("mhPosVsMomPhi");
      hPosVsMomPhi[i]->SetName(Form("%s_Type%d",hPosVsMomPhi[i]->GetName(),i));
      hPhiDiff[i] = (TH1F*)fdata[i]->Get("mhPosMomPhiDiff");
      hPhiDiff[i]->SetName(Form("%s_Type%d",hPhiDiff[i]->GetName(),i));
    }

  // resp. eff vs. track eta and vz
  printf("+++++ histograms vs. eta and vz +++++\n");
  TH1F *hTrkPtVsVzAll[2][nEta];
  TH1F *hTrkPtVsVzAcc[2][nEta];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(0)->SetRangeUser(ptCut, 20);
      hnMthEff[i]->GetAxis(8)->SetRangeUser(-35,35); // localz
      hnMthEff[i]->GetAxis(3)->SetRangeUser(3,8); // cell
      if(i==1) hnMthEff[i]->GetAxis(9)->SetRangeUser(-1*phiDiffCut, phiDiffCut); // phiDiff cut
      for(int j=0; j<nEta; j++)
	{
	  hnMthEff[i]->GetAxis(7)->SetRange(j+2, j+2);
	  hTrkPtVsVzAll[i][j] = (TH1F*)hnMthEff[i]->Projection(4);
	  hTrkPtVsVzAll[i][j]->Sumw2();
	  hTrkPtVsVzAll[i][j]->SetName(Form("hTrkPtVsVzAll_pt%1.0f_eta%d_Type%d",ptCut,j,i));

	  hnMthEff[i]->GetAxis(5)->SetRange(2, 2);
	  hTrkPtVsVzAcc[i][j] = (TH1F*)hnMthEff[i]->Projection(4);
	  hTrkPtVsVzAcc[i][j]->Sumw2();
	  hTrkPtVsVzAcc[i][j]->SetName(Form("hTrkPtVsVzAcc_pt%1.0f_eta%d_Type%d",ptCut,j,i));
	  hnMthEff[i]->GetAxis(5)->SetRange(0, -1);
	}
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(0)->SetRange(0,-1);
      hnMthEff[i]->GetAxis(3)->SetRange(0,-1);
      hnMthEff[i]->GetAxis(8)->SetRange(0,-1);
      if(i==1) hnMthEff[i]->GetAxis(9)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  // compare vs. track eta and phi
  printf("+++++ histograms vs. eta and phi +++++\n");
  TH1F *hTrkPhiAll[2];
  TH1F *hTrkPhiAcc[2];
  TH1F *hTrkEtaAll[2];
  TH1F *hTrkEtaAcc[2];
  TH1F *hTrkDCAzAll[2];
  TH1F *hTrkDCAzAcc[2];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(0)->SetRangeUser(ptCut, 20);
      hnMthEff[i]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5);
      hnMthEff[i]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut);

      hTrkDCAzAll[i] = (TH1F*)hnMthEff[i]->Projection(4);
      hTrkDCAzAll[i]->Sumw2();
      hTrkDCAzAll[i]->SetName(Form("hTrkDCAzAll_pt%1.0f_Type%d",ptCut,i));
      hTrkPhiAll[i] = (TH1F*)hnMthEff[i]->Projection(6);
      hTrkPhiAll[i]->Sumw2();
      hTrkPhiAll[i]->SetName(Form("hTrkPhiAll_pt%1.0f_Type%d",ptCut,i));
      hTrkEtaAll[i] = (TH1F*)hnMthEff[i]->Projection(7);
      hTrkEtaAll[i]->Sumw2();
      hTrkEtaAll[i]->SetName(Form("hTrkEtaAll_pt%1.0f_Type%d",ptCut,i));

      hnMthEff[i]->GetAxis(5)->SetRange(2, 2);
      hTrkDCAzAcc[i] = (TH1F*)hnMthEff[i]->Projection(4);
      hTrkDCAzAcc[i]->Sumw2();
      hTrkDCAzAcc[i]->SetName(Form("hTrkDCAzAcc_pt%1.0f_Type%d",ptCut,i));
      hTrkPhiAcc[i] = (TH1F*)hnMthEff[i]->Projection(6);
      hTrkPhiAcc[i]->Sumw2();
      hTrkPhiAcc[i]->SetName(Form("hTrkPhiAcc_pt%1.0f_Type%d",ptCut,i));
      hTrkEtaAcc[i] = (TH1F*)hnMthEff[i]->Projection(7);
      hTrkEtaAcc[i]->Sumw2();
      hTrkEtaAcc[i]->SetName(Form("hTrkEtaAcc_pt%1.0f_Type%d",ptCut,i));

      hnMthEff[i]->GetAxis(5)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(0)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(4)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  TH1F *h1tmp = 0x0;
  printf("+++++ histograms vs. vz +++++\n");
  TH1F *hTrkPtInVzAll[2][nVz];
  TH1F *hTrkPtInVzAcc[2][nVz];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut); // eta cut
      for(int k=0; k<nVz; k++)
	{
	  hnMthEff[i]->GetAxis(4)->SetRangeUser(vzMin[k]+0.5, vzMax[k]-0.5);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_%d_Vz%d_All",h1tmp->GetName(),i,k));
	  hTrkPtInVzAll[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInVzAll_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtInVzAll[i][k]->Sumw2();

	  hnMthEff[i]->GetAxis(5)->SetRange(2,2);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_%d_Vz%d_Acc",h1tmp->GetName(),i,k));
	  hTrkPtInVzAcc[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInVzAcc_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtInVzAcc[i][k]->Sumw2();
	  hnMthEff[i]->GetAxis(5)->SetRange(0,-1);
	}
      hnMthEff[i]->GetAxis(4)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  // apply phiDiff cut
  printf("+++++ From now on, phiDiffCut = %4.3f is applied in cosmic ray +++\n\n",phiDiffCut);
  hnMthEff[1]->GetAxis(9)->SetRangeUser(-1*phiDiffCut,1*phiDiffCut);

  // compare restricted phi & eta region
  printf("+++++ histograms in center +++++\n");
  TH1F *hTrkPtInCenterAll[2];
  TH1F *hTrkPtInCenterAcc[2];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5);
      hnMthEff[i]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut);
      hnMthEff[i]->GetAxis(1)->SetRange(15,17);
     
      h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
      h1tmp->SetName(Form("h1center_all_%d",i));
      hTrkPtInCenterAll[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInCenterAll_pt%1.0f_Type%d",ptCut,i), xPtBins);
      hTrkPtInCenterAll[i]->Sumw2();

      hnMthEff[i]->GetAxis(5)->SetRange(2, 2);
      h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
      h1tmp->SetName(Form("h1center_acc_%d",i));
      hTrkPtInCenterAcc[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInCenterAcc_pt%1.0f_Type%d",ptCut,i), xPtBins);
      hTrkPtInCenterAcc[i]->Sumw2();
      hnMthEff[i]->GetAxis(5)->SetRange(0, -1);

      hnMthEff[i]->GetAxis(4)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(1)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");


  printf("+++++ histograms vs. BL and cell +++++\n");
  TH3F *hTrkPtVsBlVsCellAll[2];
  TH3F *hTrkPtVsBlVsCellAcc[2];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5);
      hnMthEff[i]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut);

      hnMthEff[i]->GetAxis(2)->SetRange(1,3);
      hTrkPtVsBlVsCellAll[i] = (TH3F*)hnMthEff[i]->Projection(0,1,3);
      hTrkPtVsBlVsCellAll[i]->SetName(Form("hTrkPtVsBlVsCellAll_Type%d",i));
      hTrkPtVsBlVsCellAll[i]->Sumw2();
      hnMthEff[i]->GetAxis(2)->SetRange(4,5);
      TH3F *h3tmp = (TH3F*)hnMthEff[i]->Projection(0,1,3);
      h3tmp->SetName(Form("%s_2",hTrkPtVsBlVsCellAll[i]->GetName()));
      h3tmp->Sumw2();
      TH3F *h3clone = (TH3F*)h3tmp->Clone(Form("%s_clone",h3tmp->GetName()));
      h3clone->Reset("AC");
      int nbinsx = h3tmp->GetNbinsX();
      int nbinsy = h3tmp->GetNbinsY();
      int nbinsz = h3tmp->GetNbinsZ();
      for(int ibin = 1; ibin<=nbinsx; ibin++)
	{
	  for(int jbin=1; jbin<=nbinsy; jbin++)
	    {
	      for(int kbin=1; kbin<=nbinsz; kbin++)
		{
		  h3clone->SetBinContent(ibin, jbin, kbin, h3tmp->GetBinContent(ibin, jbin, nbinsz-kbin+1));
		  h3clone->SetBinError(ibin, jbin, kbin, h3tmp->GetBinError(ibin, jbin, nbinsz-kbin+1));
		}
	    }
	}
      hTrkPtVsBlVsCellAll[i]->Add(h3clone);
      hnMthEff[i]->GetAxis(2)->SetRange(0, -1);


      hnMthEff[i]->GetAxis(5)->SetRange(2, 2);
      hnMthEff[i]->GetAxis(2)->SetRange(1,3);
      hTrkPtVsBlVsCellAcc[i] = (TH3F*)hnMthEff[i]->Projection(0,1,3);
      hTrkPtVsBlVsCellAcc[i]->SetName(Form("hTrkPtVsBlVsCellAcc_Type%d",i));
      hTrkPtVsBlVsCellAcc[i]->Sumw2();
      hnMthEff[i]->GetAxis(2)->SetRange(4,5);
      h3tmp = (TH3F*)hnMthEff[i]->Projection(0,1,3);
      h3tmp->SetName(Form("%s_2",hTrkPtVsBlVsCellAcc[i]->GetName()));
      h3tmp->Sumw2();
      h3clone = (TH3F*)h3tmp->Clone(Form("%s_clone",h3tmp->GetName()));
      h3clone->Reset("AC");
      for(int ibin = 1; ibin<=nbinsx; ibin++)
	{
	  for(int jbin=1; jbin<=nbinsy; jbin++)
	    {
	      for(int kbin=1; kbin<=nbinsz; kbin++)
		{
		  h3clone->SetBinContent(ibin, jbin, kbin, h3tmp->GetBinContent(ibin, jbin, nbinsz-kbin+1));
		  h3clone->SetBinError(ibin, jbin, kbin, h3tmp->GetBinError(ibin, jbin, nbinsz-kbin+1));
		}
	    }
	}
      hTrkPtVsBlVsCellAcc[i]->Add(h3clone);
      hnMthEff[i]->GetAxis(2)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(5)->SetRange(0, -1);

      hnMthEff[i]->GetAxis(8)->SetRange(0, -1); // localz
      hnMthEff[i]->GetAxis(3)->SetRange(0, -1); // cell
      hnMthEff[i]->GetAxis(4)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  printf("+++++ histograms vs. BL +++++\n");
  TH1F *hTrkPtInBlAll[3][30];
  TH1F *hTrkPtInBlAcc[3][30];
  TH1F *hWeightVsCell[30];
  int lowPtBin = hTrkPtVsBlVsCellAll[0]->GetXaxis()->FindFixBin(ptCut);

  for(int i=0; i<2; i++)
    {
      for(int j=0; j<30; j++)
	{
	  h1tmp = (TH1F*)hTrkPtVsBlVsCellAll[i]->ProjectionX(Form("hTrkPtInBlAll_BL%d_Type%d_tmp",j+1,i),j+1,j+1,0,-1);
	  hTrkPtInBlAll[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInBlAll_BL%d_Type%d",j+1,i), xPtBins);

	  h1tmp = (TH1F*)hTrkPtVsBlVsCellAcc[i]->ProjectionX(Form("hTrkPtInBlAcc_BL%d_Type%d_tmp",j+1,i),j+1,j+1,0,-1);
	  hTrkPtInBlAcc[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtInBlAcc_BL%d_Type%d",j+1,i), xPtBins);

	  if(i==0)
	    {
	      hWeightVsCell[j] = (TH1F*)hTrkPtVsBlVsCellAll[i]->ProjectionZ(Form("hWeightVsCell_BL%d_Type%d",j+1,i),lowPtBin,-1,j+1,j+1);
	      if(hWeightVsCell[j]->Integral()>0) hWeightVsCell[j]->Scale(1./hWeightVsCell[j]->Integral());
	    }
	  else
	    {
	      TH1F *h1All = (TH1F*)h1tmp->Clone( Form("hTrkPtInBlAll_BL%d_Type%d_Weight_tmp",j+1,i));
	      h1All->Reset("AC");
	      TH1F *h1Acc = (TH1F*)h1tmp->Clone( Form("hTrkPtInBlAcc_BL%d_Type%d_Weight_tmp",j+1,i));
	      h1Acc->Reset("AC");
	      for(int k=0; k<hWeightVsCell[j]->GetNbinsX(); k++)
		{
		  h1tmp = (TH1F*)hTrkPtVsBlVsCellAll[i]->ProjectionX(Form("hTrkPtInBlAll_BL%d_Cell%d_Type%d_tmp",j+1,k,i),j+1,j+1,k+1,k+1);
		  double integral = h1tmp->Integral();
		  if(integral>0)
		    {
		      h1All->Add(h1tmp, hWeightVsCell[j]->GetBinContent(k+1)/integral);
		    }
		  h1tmp = (TH1F*)hTrkPtVsBlVsCellAcc[i]->ProjectionX(Form("hTrkPtInBlAcc_BL%d_Cell%d_Type%d_tmp",j+1,k,i),j+1,j+1,k+1,k+1);
		  if(integral>0)
		    {
		      h1Acc->Add(h1tmp, hWeightVsCell[j]->GetBinContent(k+1)/integral);
		    }
		}
	      hTrkPtInBlAll[2][j] = (TH1F*)h1All->Rebin(nPtBins, Form("hTrkPtInBlAll_BL%d_Type%d_Weight",j+1,i), xPtBins);
	      hTrkPtInBlAcc[2][j] = (TH1F*)h1Acc->Rebin(nPtBins, Form("hTrkPtInBlAcc_BL%d_Type%d_Weight",j+1,i), xPtBins);
	    }
	} 
    }
  printf("+++++ Done! +++++\n\n");

  printf("+++++ histograms vs. cell +++++\n");
  // efficiency vs. cell
  TH1F *hTrkPtVsCellAll[2][30];
  TH1F *hTrkPtVsCellAcc[2][30];
  TH1F *hTrkPtInCellAll[2][30][nCellBin];
  TH1F *hTrkPtInCellAcc[2][30][nCellBin];
  int bin_index = hTrkPtInBlAll[0][0]->FindFixBin(ptCut);
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<30; j++)
	{
	  hTrkPtVsCellAll[i][j] = (TH1F*)hTrkPtVsBlVsCellAll[i]->ProjectionZ(Form("hTrkPtVsCellAll_BL%d_Type%d",j+1,i), lowPtBin, hTrkPtVsBlVsCellAll[i]->GetNbinsX()+1, j+1, j+1);
	  hTrkPtVsCellAcc[i][j] = (TH1F*)hTrkPtVsBlVsCellAcc[i]->ProjectionZ(Form("hTrkPtVsCellAcc_BL%d_Type%d",j+1,i), lowPtBin, hTrkPtVsBlVsCellAcc[i]->GetNbinsX()+1, j+1, j+1);

	  for(int k=0; k<nCellBin; k++)
	    {
	      h1tmp = (TH1F*)hTrkPtVsBlVsCellAll[i]->ProjectionX(Form("hTrkPtAll_BL%d_Cell%d_Type%d_tmp",j+1,k,i), j+1, j+1, start_cell[k], end_cell[k]);
	      hTrkPtInCellAll[i][j][k] =  (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAll_BL%d_Cell%d_Type%d",j+1,k,i), xPtBins);

	      h1tmp = (TH1F*)hTrkPtVsBlVsCellAcc[i]->ProjectionX(Form("hTrkPtAcc_BL%d_Cell%d_Type%d_tmp",j+1,k,i), j+1, j+1, start_cell[k], end_cell[k]);
	      hTrkPtInCellAcc[i][j][k] =  (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAcc_BL%d_Cell%d_Type%d",j+1,k,i), xPtBins);
	    }
	}
    }
  printf("+++++ Done! +++++\n\n");

  if(saveHisto)
    {
      printf("+++++ save histograms +++++\n");
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdMthEff.root",run_type),"update");
      for(int i=0; i<2; i++)
	{
	  hPosVsMomPhi[i]->Write("",TObject::kOverwrite);
	  hPhiDiff[i]->Write("",TObject::kOverwrite);

	  hTrkPtInCenterAll[i]->Write("",TObject::kOverwrite);
	  hTrkPtInCenterAcc[i]->Write("",TObject::kOverwrite);

	  for(int j=0; j<nEta; j++)
	    {
	      hTrkPtVsVzAll[i][j]->Write("",TObject::kOverwrite);
	      hTrkPtVsVzAcc[i][j]->Write("",TObject::kOverwrite);
	    }
	  hTrkDCAzAll[i]->Write("",TObject::kOverwrite);
	  hTrkDCAzAcc[i]->Write("",TObject::kOverwrite);
	  hTrkPhiAll[i]->Write("",TObject::kOverwrite);
	  hTrkPhiAcc[i]->Write("",TObject::kOverwrite);
	  hTrkEtaAll[i]->Write("",TObject::kOverwrite);
	  hTrkEtaAcc[i]->Write("",TObject::kOverwrite);
	  for(int k=0; k<nVz; k++)
	    {
	      hTrkPtInVzAll[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtInVzAcc[i][k]->Write("",TObject::kOverwrite);
	    }

	  for(int j=0; j<30; j++)
	    {
	      hTrkPtInBlAll[i][j]->Write("",TObject::kOverwrite);
	      hTrkPtInBlAcc[i][j]->Write("",TObject::kOverwrite);
	      if(i==1)
		{
		  hTrkPtInBlAll[2][j]->Write("",TObject::kOverwrite);
		  hTrkPtInBlAcc[2][j]->Write("",TObject::kOverwrite);
		}

	      hTrkPtVsCellAll[i][j]->Write("",TObject::kOverwrite);
	      hTrkPtVsCellAcc[i][j]->Write("",TObject::kOverwrite);

	      for(int k=0; k<5; k++)
		{
		  hTrkPtInCellAll[i][j][k]->Write("",TObject::kOverwrite);
		  hTrkPtInCellAcc[i][j][k]->Write("",TObject::kOverwrite);
		}
	    }
	}
    }
}


//================================================
void makeHisto3(const int saveHisto = 1)
{
  const TString legName[2] = {"Embedding","Cosmic ray"}
  TFile *fdata[2];
  fdata[0] = TFile::Open(Form("output/%s.Embed.Jpsi.root",run_type),"read");
  fdata[1] = TFile::Open(Form("output/Run%d.cosmic.root",year-2000),"read");
  THnSparseF *hnMthEff[2];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i] = (THnSparseF*)fdata[i]->Get("mhMtdMatchEff");
      hnMthEff[i]->SetName(Form("%s_Type%d",hnMthEff[i]->GetName(),i));
      printf("[i] Input file %s\n",fdata[i]->GetName());
    }

  TH1F *h1tmp = 0x0;

  /*
  // response & match eff vs. track eta and vz
  printf("+++++ histograms vs. vz +++++\n");
  TH1F *hTrkPtRespInVzAll[2][nVz];
  TH1F *hTrkPtRespInVzAcc[2][nVz];
  TH1F *hTrkPtMthInVzAll[2][nVz];
  TH1F *hTrkPtMthInVzAcc[2][nVz];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(1)->SetRange(9, 23);
      hnMthEff[i]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut); // eta cut
      for(int k=0; k<nVz; k++)
	{
	  hnMthEff[i]->GetAxis(4)->SetRangeUser(vzMin[k]+0.5, vzMax[k]-0.5);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Mth_%d_Vz%d_All",h1tmp->GetName(),i,k));
	  hTrkPtMthInVzAll[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAll_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtMthInVzAll[i][k]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRangeUser(-35, 35); // localz
	  hnMthEff[i]->GetAxis(3)->SetRangeUser(3,8); // cell
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Resp_%d_Vz%d_All",h1tmp->GetName(),i,k));
	  hTrkPtRespInVzAll[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAll_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtRespInVzAll[i][k]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(3)->SetRange(0,-1);


	  hnMthEff[i]->GetAxis(5)->SetRange(2,2);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Mth_%d_Vz%d_Acc",h1tmp->GetName(),i,k));
	  hTrkPtMthInVzAcc[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAcc_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtMthInVzAcc[i][k]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRangeUser(-35, 35); // localz
	  hnMthEff[i]->GetAxis(3)->SetRangeUser(3,8); // cell
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Resp_%d_Vz%d_Acc",h1tmp->GetName(),i,k));
	  hTrkPtRespInVzAcc[i][k] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAcc_Vz%d_Type%d",k,i), xPtBins);
	  hTrkPtRespInVzAcc[i][k]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(3)->SetRange(0,-1);

	  hnMthEff[i]->GetAxis(4)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(5)->SetRange(0, -1);
	}
      hnMthEff[i]->GetAxis(1)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(7)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  printf("+++++ histograms vs. eta +++++\n");
  TH1F *hTrkPtRespInEtaAll[2][nEta];
  TH1F *hTrkPtRespInEtaAcc[2][nEta];
  TH1F *hTrkPtMthInEtaAll[2][nEta];
  TH1F *hTrkPtMthInEtaAcc[2][nEta];
  for(int i=0; i<2; i++)
    {
      hnMthEff[i]->GetAxis(1)->SetRange(9, 23);
      hnMthEff[i]->GetAxis(4)->SetRangeUser(-50,50); // vz cut
      for(int j=0; j<nEta; j++)
	{
	  hnMthEff[i]->GetAxis(7)->SetRange(j+2, j+2);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Mth_%d_Eta%d_All",h1tmp->GetName(),i, j));
	  hTrkPtRespInEtaAll[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAll_Eta%d_Type%d",j,i), xPtBins);
	  hTrkPtRespInEtaAll[i][j]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRangeUser(-35, 35); // localz
	  hnMthEff[i]->GetAxis(3)->SetRangeUser(3,8); // cell
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Resp_%d_Eta%d_All",h1tmp->GetName(),i,j));
	  hTrkPtMthInEtaAll[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAll_Eta%d_Type%d",j,i), xPtBins);
	  hTrkPtMthInEtaAll[i][j]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(3)->SetRange(0,-1);

	  hnMthEff[i]->GetAxis(5)->SetRange(2,2);
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Mth_%d_Eta%d_Acc",h1tmp->GetName(),i,j));
	  hTrkPtRespInEtaAcc[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAcc_Eta%d_Type%d",j,i), xPtBins);
	  hTrkPtRespInEtaAcc[i][j]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRangeUser(-35, 35); // localz
	  hnMthEff[i]->GetAxis(3)->SetRangeUser(3,8); // cell
	  h1tmp = (TH1F*)hnMthEff[i]->Projection(0);
	  h1tmp->SetName(Form("%s_Resp_%d_Eta%d_Acc",h1tmp->GetName(),i,j));
	  hTrkPtMthInEtaAcc[i][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAcc_Eta%d_Type%d",j,i), xPtBins);
	  hTrkPtMthInEtaAcc[i][j]->Sumw2();
	  hnMthEff[i]->GetAxis(8)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(3)->SetRange(0,-1);
	  hnMthEff[i]->GetAxis(5)->SetRange(0, -1);
	}
      hnMthEff[i]->GetAxis(1)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(4)->SetRange(0, -1);
    }
  printf("+++++ Done! +++++\n\n");

  printf("+++++ histograms vs. phiDiffCut +++++\n");
  // check the effect of phiDiff cut
  TH1F *hTrkPtRespInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtRespInPhiDiffAcc[nPhiDiffCut];
  TH1F *hTrkPtMthInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtMthInPhiDiffAcc[nPhiDiffCut];

  hnMthEff[1]->GetAxis(1)->SetRange(9, 23); // bottom backleg
  hnMthEff[1]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5); // vz cut
  hnMthEff[1]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut); // eta cut
  for(int i=0; i<nPhiDiffCut; i++)
    {
      hnMthEff[1]->GetAxis(9)->SetRangeUser(-1*phiDiffCuts[i]+1e-4,1*phiDiffCuts[i]-1e-4); // phiDiff cut
      
      h1tmp = (TH1F*)hnMthEff[1]->Projection(0);
      h1tmp->SetName(Form("%s_Mth_%d_All",h1tmp->GetName(),i));
      hTrkPtMthInPhiDiffAll[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAll_PhiDiff%d_Type%d",i,1), xPtBins);
      hTrkPtMthInPhiDiffAll[i]->Sumw2();
      hnMthEff[1]->GetAxis(8)->SetRangeUser(-35,35); // localz
      hnMthEff[1]->GetAxis(3)->SetRangeUser(3,8); // cell
      h1tmp = (TH1F*)hnMthEff[1]->Projection(0);
      h1tmp->SetName(Form("%s_Resp_%d_All",h1tmp->GetName(),i));
      hTrkPtRespInPhiDiffAll[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAll_PhiDiff%d_Type%d",i,1), xPtBins);
      hTrkPtRespInPhiDiffAll[i]->Sumw2();
      hnMthEff[1]->GetAxis(8)->SetRange(0,-1);
      hnMthEff[1]->GetAxis(3)->SetRange(0,-1);

      hnMthEff[1]->GetAxis(5)->SetRange(2,2);
      h1tmp = (TH1F*)hnMthEff[1]->Projection(0);
      h1tmp->SetName(Form("%s_Mth_%d_Acc",h1tmp->GetName(),i));
      hTrkPtMthInPhiDiffAcc[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtMthAcc_PhiDiff%d_Type%d",i,1), xPtBins);
      hTrkPtMthInPhiDiffAcc[i]->Sumw2();
      hnMthEff[1]->GetAxis(8)->SetRangeUser(-35,35); // localz
      hnMthEff[1]->GetAxis(3)->SetRangeUser(3,8); // cell
      h1tmp = (TH1F*)hnMthEff[1]->Projection(0);
      h1tmp->SetName(Form("%s_Resp_%d_Acc",h1tmp->GetName(),i));
      hTrkPtRespInPhiDiffAcc[i] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtRespAcc_PhiDiff%d_Type%d",i,1), xPtBins);
      hTrkPtRespInPhiDiffAcc[i]->Sumw2();
      hnMthEff[1]->GetAxis(8)->SetRange(0,-1);
      hnMthEff[1]->GetAxis(3)->SetRange(0,-1);
      hnMthEff[1]->GetAxis(5)->SetRange(0, -1);

      hnMthEff[1]->GetAxis(9)->SetRange(0,-1);
    }
  hnMthEff[1]->GetAxis(4)->SetRange(0,-1);
  hnMthEff[1]->GetAxis(1)->SetRange(0, -1);
  hnMthEff[1]->GetAxis(7)->SetRange(0, -1);
  printf("+++++ Done! +++++\n\n");
  */

  // check the effect of phiDiffCut   
  hnMthEff[1]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5); // vz cut
  hnMthEff[1]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut); // eta cut
  TH2F *hCellVsPhiDiff = (TH2F*)hnMthEff[1]->Projection(3, 9);
  hCellVsPhiDiff->SetName(Form("hCellVsPhiDiff_BtmBL_Type%d",1));
  TH1F *hTrkPtVsCellInPhiDiff[nPhiDiffCut];
  TH3F *hTrkPtVsBlVsCellAll[nPhiDiffCut];
  TH3F *hTrkPtVsBlVsCellAcc[nPhiDiffCut];
  for(int k=0; k<nPhiDiffCut; k++)
    {	      
      hnMthEff[1]->GetAxis(9)->SetRangeUser(-1*phiDiffCuts[k],1*phiDiffCuts[k]);
      hTrkPtVsCellInPhiDiff[k] = (TH1F*)hnMthEff[1]->Projection(3);
      hTrkPtVsCellInPhiDiff[k]->SetName(Form("hTrkPtVsCell_BtmBL_PhiDiff%d_Type1",k));
	      
      hTrkPtVsBlVsCellAll[k] = (TH3F*)hnMthEff[1]->Projection(0,1,3);
      hTrkPtVsBlVsCellAll[k]->SetName(Form("hTrkPtVsBlVsCellAll_PhiDiff%d_Type1",k));
      hTrkPtVsBlVsCellAll[k]->Sumw2();
      hnMthEff[1]->GetAxis(5)->SetRange(2, 2);
      hTrkPtVsBlVsCellAcc[k] = (TH3F*)hnMthEff[1]->Projection(0,1,3);
      hTrkPtVsBlVsCellAcc[k]->SetName(Form("hTrkPtVsBlVsCellAcc_PhiDiff%d_Type1",k));
      hTrkPtVsBlVsCellAcc[k]->Sumw2();
      hnMthEff[1]->GetAxis(5)->SetRange(0, -1);
    }
  hnMthEff[1]->GetAxis(9)->SetRange(0,-1);

  TH1F *hTrkPtInCellInPhiDiffAll[nPhiDiffCut][nCellBin];
  TH1F *hTrkPtInCellInPhiDiffAcc[nPhiDiffCut][nCellBin];
  for(int k=0; k<nPhiDiffCut; k++)
    {
      for(int j=0; j<nCellBin; j++)
	{
	  h1tmp = (TH1F*)hTrkPtVsBlVsCellAll[k]->ProjectionX(Form("hTrkPtAll_BtmBL_Cell%d_PhiDiff%d_Type1_tmp",j,k),1,hTrkPtVsBlVsCellAll[k]->GetNbinsY(),start_cell[j],end_cell[j]);
	  hTrkPtInCellInPhiDiffAll[k][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAll_BtmBL_Cell%d_PhiDiff%d_Type1",j,k), xPtBins);

	  h1tmp = (TH1F*)hTrkPtVsBlVsCellAcc[k]->ProjectionX(Form("hTrkPtAcc_BtmBL_Cell%d_PhiDiff%d_Type1_tmp",j,k),1,hTrkPtVsBlVsCellAcc[k]->GetNbinsY(),start_cell[j],end_cell[j]);
	  hTrkPtInCellInPhiDiffAcc[k][j] = (TH1F*)h1tmp->Rebin(nPtBins, Form("hTrkPtAcc_BtmBL_Cell%d_PhiDiff%d_Type1",j,k), xPtBins);
	}
    }
  

  //get the weights
  hnMthEff[0]->GetAxis(4)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5); // vz cut
  hnMthEff[0]->GetAxis(7)->SetRangeUser(-1*etaCut,1*etaCut); // eta cut
  TH2F *hTrkAllBLvsCell = (TH2F*)hnMthEff[0]->Projection(3,1);
  hTrkAllBLvsCell->Scale(1./hTrkAllBLvsCell->Integral()); 
  hnMthEff[0]->GetAxis(4)->SetRange(0,-1);
  hnMthEff[0]->GetAxis(7)->SetRange(0,-1);

  TH1F *hTrkPtInPhiDiffAll[nPhiDiffCut];
  TH1F *hTrkPtInPhiDiffAcc[nPhiDiffCut];
  TH1F *h1tmpAll, *h1tmpAcc;
  for(int k=0; k<nPhiDiffCut; k++)
    {	
      hTrkPtInPhiDiffAll[k] = (TH1F*)hTrkPtVsBlVsCellAll[k]->ProjectionX(Form("hTrkPtAll_BtmBL_PhiDiff%d_Type1",k));
      hTrkPtInPhiDiffAll[k]->Reset("AC");
      hTrkPtInPhiDiffAcc[k] = (TH1F*)hTrkPtVsBlVsCellAll[k]->ProjectionX(Form("hTrkPtAcc_BtmBL_PhiDiff%d_Type1",k));
      hTrkPtInPhiDiffAcc[k]->Reset("AC");
      for(int bl=10; bl<=22; bl++) 
	{
	  for(int cell=0; cell<18; cell++)
	    {
	      h1tmpAll = (TH1F*)hTrkPtVsBlVsCellAll[k]->ProjectionX(Form("hTrkPtAll_BL%d_Cell%d_PhiDiff%d_Type1",bl,cell,k), bl, bl, cell+1, cell+1);
	      h1tmpAcc = (TH1F*)hTrkPtVsBlVsCellAcc[k]->ProjectionX(Form("hTrkPtAcc_BL%d_Cell%d_PhiDiff%d_Type1",bl,cell,k), bl, bl, cell+1, cell+1);
	      double weight = hTrkAllBLvsCell->GetBinContent(bl,cell+1)/h1tmpAll->Integral();
	      hTrkPtInPhiDiffAll[k]->Add(h1tmpAll, weight);
	      hTrkPtInPhiDiffAcc[k]->Add(h1tmpAcc, weight);
	    }
	}
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdMthEff.root",run_type),"update");
      printf("+++++ save histograms in %s +++++\n",fout->GetName());

      hCellVsPhiDiff->Write("",TObject::kOverwrite);
      for(int k=0; k<nPhiDiffCut; k++)
	{
	  hTrkPtVsCellInPhiDiff[k]->Write("",TObject::kOverwrite);
	  hTrkPtInPhiDiffAll[k]->Write("",TObject::kOverwrite);
	  hTrkPtInPhiDiffAcc[k]->Write("",TObject::kOverwrite);
	  for(int j=0; j<nCellBin; j++)
	    {
	      hTrkPtInCellInPhiDiffAll[k][j]->Write("",TObject::kOverwrite);
	      hTrkPtInCellInPhiDiffAcc[k][j]->Write("",TObject::kOverwrite);	      
	    }
	}

      /*
      for(int i=0; i<2; i++)
	{
	  for(int k=0; k<nVz; k++)
	    {
	      hTrkPtRespInVzAll[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtRespInVzAcc[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtMthInVzAll[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtMthInVzAcc[i][k]->Write("",TObject::kOverwrite);
	    }

	  for(int k=0; k<nEta; k++)
	    {
	      hTrkPtRespInEtaAll[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtRespInEtaAcc[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtMthInEtaAll[i][k]->Write("",TObject::kOverwrite);
	      hTrkPtMthInEtaAcc[i][k]->Write("",TObject::kOverwrite);
	    }
	}

      for(int i=0; i<nPhiDiffCut; i++)
	{
	  hTrkPtRespInPhiDiffAll[i]->Write("",TObject::kOverwrite);
	  hTrkPtRespInPhiDiffAcc[i]->Write("",TObject::kOverwrite);
	  hTrkPtMthInPhiDiffAll[i]->Write("",TObject::kOverwrite);
	  hTrkPtMthInPhiDiffAcc[i]->Write("",TObject::kOverwrite);
	}
      */
      fout->Close();
    }
}

//================================================
TH1F *reverseCell1D(const TH1F *h1)
{
  TH1F *h1clone = (TH1F*)h1->Clone(Form("%s_clone",h1->GetName()));
  h1clone->Reset("AC");
  int nbinsx = h1->GetNbinsX();
  for(int ibin=1; ibin<=nbinsx; ibin++)
    {
      h1clone->SetBinContent(ibin, h1->GetBinContent(nbinsx-ibin+1));
      h1clone->SetBinError(ibin, h1->GetBinError(nbinsx-ibin+1));
    }
  return h1clone;
}


//================================================
TH2F *reverseCell2D(const TH2F *h2, const int axis)
{
  TH2F *h2clone = (TH2F*)h2->Clone(Form("%s_clone",h2->GetName()));
  h2clone->Reset("AC");
  int nbinsx = h2->GetNbinsX();
  int nbinsy = h2->GetNbinsY();
  for(int ibin=1; ibin<=nbinsx; ibin++)
    {
      for(int jbin=1; jbin<=nbinsy; jbin++)
	{
	  if(axis==0)
	    {
	      h2clone->SetBinContent(ibin, jbin, h2->GetBinContent(nbinsx-ibin+1, jbin));
	      h2clone->SetBinError(ibin, jbin, h2->GetBinError(nbinsx-ibin+1, jbin));
	    }
	  else if(axis==1)
	    {
	      h2clone->SetBinContent(ibin, jbin, h2->GetBinContent(ibin, nbinsy-jbin+1));
	      h2clone->SetBinError(ibin, jbin, h2->GetBinError(ibin, nbinsy-jbin+1));
	    }
	}
    }
  return h2clone;
}

//================================================
TH3F *reverseCell3D(const TH3F *h3, const int axis)
{
  TH3F *h3clone = (TH3F*)h3->Clone(Form("%s_clone",h3->GetName()));
  h3clone->Reset("AC");
  int nbinsx = h3->GetNbinsX();
  int nbinsy = h3->GetNbinsY();
  int nbinsz = h3->GetNbinsZ();
  for(int ibin=1; ibin<=nbinsx; ibin++)
    {
      for(int jbin=1; jbin<=nbinsy; jbin++)
	{
	  for(int kbin=1; kbin<=nbinsz; kbin++)
	    {
	      if(axis==0)
		{
		  h3clone->SetBinContent(ibin, jbin, kbin, h3->GetBinContent(nbinsx-ibin+1, jbin, kbin));
		  h3clone->SetBinError(ibin, jbin, kbin, h3->GetBinError(nbinsx-ibin+1, jbin, kbin));
		}
	      else if(axis==1)
		{
		  h3clone->SetBinContent(ibin, jbin, kbin, h3->GetBinContent(ibin, nbinsy-jbin+1, kbin));
		  h3clone->SetBinError(ibin, jbin, kbin, h3->GetBinError(ibin, nbinsy-jbin+1, kbin));
		}
	      else if(axis==2)
		{
		  h3clone->SetBinContent(ibin, jbin, kbin, h3->GetBinContent(ibin, jbin, nbinsz-kbin+1));
		  h3clone->SetBinError(ibin, jbin, kbin, h3->GetBinError(ibin, jbin, nbinsz-kbin+1));
		}
	    }
	}
    }
  return h3clone;
}

//================================================
void makeHisto2(const int saveHisto = 1)
{
  TFile *fdata[nYears];
  THnSparseF *hnMthEff[nYears];
  TH1F *hDayAll[nYears];
  TH1F *hDayAcc[nYears];

  for(int i=0; i<nYears; i++)
    {
      cout << "Processing " << years[i] << endl;
      fdata[i] = TFile::Open(Form("output/Run%d.cosmic.root",years[i]-2000),"read");
      hnMthEff[i] = (THnSparseF*)fdata[i]->Get("mhMtdMatchEff");
      hnMthEff[i]->SetName(Form("%s_Year%d",hnMthEff[i]->GetName(),years[i]));
      hnMthEff[i]->GetAxis(2)->SetRangeUser(-1*vzCut+0.5, vzCut-0.5);
      hnMthEff[i]->GetAxis(0)->SetRangeUser(ptCut, 15);

      hDayAll[i] = (TH1F*)hnMthEff[i]->Projection(6);
      hDayAll[i]->Sumw2();
      hDayAll[i]->SetName(Form("hTrkDayAll_pt%1.0f_Year%d",ptCut,years[i]));

      hnMthEff[i]->GetAxis(3)->SetRange(2, 2);
      hDayAcc[i] = (TH1F*)hnMthEff[i]->Projection(6);
      hDayAcc[i]->Sumw2();
      hDayAcc[i]->SetName(Form("hTrkDayAcc_pt%1.0f_Year%d",ptCut,years[i]));

      hnMthEff[i]->GetAxis(3)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(0)->SetRange(0, -1);
      hnMthEff[i]->GetAxis(2)->SetRange(0, -1);
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdMthEff.root",run_type),"update");
      for(int i=0; i<nYears; i++)
	{
	  hDayAll[i]->Write("",TObject::kOverwrite);
	  hDayAcc[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void respEff(const int savePlot = 0)
{
  TFile *fcosmic  = TFile::Open("output/Run14.cosmic.root","read");
  TFile *fcosmic2 = TFile::Open("output/Run14.cosmic.MtdAcc.root","read");
 
  TH2F *hProjTrkPtVsBL = (TH2F*)fcosmic->Get("mhProjTrkPtVsBL");
  TH2F *hMthTrkPtVsBL = (TH2F*)fcosmic->Get("mhMthTrkPtVsBL");
  TH2F *hProjTrkPtVsBL2 = (TH2F*)fcosmic2->Get("mhProjTrkPtVsBL");
  TH2F *hMthTrkPtVsBL2 = (TH2F*)fcosmic2->Get("mhMthTrkPtVsBL");

  const char *anaName[2] = {"R","T"};
  TH1F *hProjTrkPt[2][30][5];
  TH1F *hMthTrkPt[2][30][5];
  TH1F *hMthEff[2][30][5];
  TF1  *funcMthEff[2][30][5];
  TCanvas *cFit[2][6];
  TH1F *hMthEffVsBL[3];
  for(int k=0; k<2; k++)
    {
      hMthEffVsBL[k] = new TH1F(Form("hMthEffVsBL_%s",anaName[k]),Form("hMthEffVsBL_%s",anaName[k]),150,0,150);
      int nAll = 0, nAcc = 0;
      for(int i=0; i<30; i++)
	{
	  if(i%5==0) 
	    {
	      cFit[k][i/5] = new TCanvas(Form("Fit_BL%d-%d_%s",i+1,i+6,anaName[k]), Form("Fit_BL%d-%d_%s",i+1,i+6,anaName[k]), 1100, 700);
	      cFit[k][i/5]->Divide(5, 5);
	    }
	  for(int j=0; j<5; j++)
	    {
	      if(k==0)
		{
		  hProjTrkPt[k][i][j] = (TH1F*)hProjTrkPtVsBL->ProjectionY(Form("hProjTrkPt_BL%d_Mod%d_%s",i+1,j+1,anaName[k]),i*5+j+1,i*5+j+1);
		  hMthTrkPt[k][i][j] = (TH1F*)hMthTrkPtVsBL->ProjectionY(Form("hMthTrkPt_BL%d_Mod%d_%s",i+1,j+1,anaName[k]),i*5+j+1,i*5+j+1);
		}
	      else
		{
		  hProjTrkPt[k][i][j] = (TH1F*)hProjTrkPtVsBL2->ProjectionY(Form("hProjTrkPt_BL%d_Mod%d_%s",i+1,j+1,anaName[k]),i*5+j+1,i*5+j+1);
		  hMthTrkPt[k][i][j] = (TH1F*)hMthTrkPtVsBL2->ProjectionY(Form("hMthTrkPt_BL%d_Mod%d_%s",i+1,j+1,anaName[k]),i*5+j+1,i*5+j+1);
		}
	      // else
	      // 	{
	      // 	  hProjTrkPt[k][i][j] = (TH1F*)fcosmic2->Get(Form("hPtMtdProBkl%d_Mod%d_2",i,j));
	      // 	  hMthTrkPt[k][i][j] = (TH1F*)fcosmic2->Get(Form("hPtMtdMatBkl%d_Mod%d_2",i,j));
	      // 	}
	      nAll += hProjTrkPt[k][i][j]->Integral(hProjTrkPt[k][i][j]->FindBin(0),hProjTrkPt[k][i][j]->FindBin(20));
	      nAcc += hMthTrkPt[k][i][j]->Integral(hMthTrkPt[k][i][j]->FindBin(0),hMthTrkPt[k][i][j]->FindBin(20));
	      hProjTrkPt[k][i][j]->Sumw2();
	      hProjTrkPt[k][i][j]->Rebin(10);
	      hMthTrkPt[k][i][j]->Sumw2();
	      hMthTrkPt[k][i][j]->Rebin(10);
	      hMthEff[k][i][j] = (TH1F*)hMthTrkPt[k][i][j]->Clone(Form("hMthEff_BL%d_Mod%d_%s",i+1,j+1,anaName[k]));
	      hMthEff[k][i][j]->Divide(hProjTrkPt[k][i][j]);
	      funcMthEff[k][i][j] = new TF1(Form("funcMthEff_BL%d_Mod%d_%s",i+1,j+1,anaName[k]), "pol0", 5, 20);
	      hMthEff[k][i][j]->Fit(funcMthEff[k][i][j], "R0Q");

	      cFit[k][i/5]->cd(i%5*5+j+1);
	      hMthEff[k][i][j]->SetTitle("");
	      hMthEff[k][i][j]->SetMarkerStyle(24);
	      hMthEff[k][i][j]->GetYaxis()->SetRangeUser(0, 1);
	      hMthEff[k][i][j]->Draw();
	      funcMthEff[k][i][j]->SetLineStyle(2);
	      funcMthEff[k][i][j]->SetLineColor(2);
	      funcMthEff[k][i][j]->Draw("sames");
	      hMthEffVsBL[k]->SetBinContent(i*5+j+1, funcMthEff[k][i][j]->GetParameter(0));
	      hMthEffVsBL[k]->SetBinError(i*5+j+1,   funcMthEff[k][i][j]->GetParError(0));

	      TPaveText *t1 = GetTitleText(Form("BL = %d, Mod = %d",i+1, j+1), 0.09);
	      t1->Draw();
	    }
	}
      printf("[i] %s: %d/%d = %4.2f%%\n",anaName[k],nAcc,nAll,nAcc*1./nAll*100);
      if(savePlot) 
	{
	  for(int i=0; i<6; i++)
	    {
	      cFit[k][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MtdRespEff_FitBL%d-%d_%s.pdf",run_type,i*5+1,i*5+5,anaName[k]));
	    }
	}
    }

  // Takahito's results
  TFile *fTak = TFile::Open("Rootfiles/Run14ResponseEffViaPtTemplate.root","read");
  hMthEffVsBL[2] = (TH1F*)fTak->Get("hSclFac");
  
  for(int i=0; i<3; i++)
    {
      hMthEffVsBL[i]->SetMarkerStyle(20+2*i);
      hMthEffVsBL[i]->SetMarkerColor(pow(2,i));
      hMthEffVsBL[i]->SetLineColor(pow(2,i));
      hMthEffVsBL[i]->GetYaxis()->SetRangeUser(0, 1);
    }
  c = draw1D(hMthEffVsBL[0],"MTD response efficiency for p_{T} > 5 GeV/c;(BL-1)*5+Mod;Eff");
  hMthEffVsBL[1]->Draw("sames");
  hMthEffVsBL[2]->Draw("sames");
  TLegend *leg = new TLegend(0.5,0.15,0.8,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Run14 cosmic ray");
  // leg->AddEntry(hMthEffVsBL[0], "Rongrong", "P");
  // leg->AddEntry(hMthEffVsBL[1], "Takahito by Rongrong", "P");
  // leg->AddEntry(hMthEffVsBL[2], "Takahito by Takahito", "P");
  leg->AddEntry(hMthEffVsBL[2], "Takahito: default", "P");
  leg->AddEntry(hMthEffVsBL[0], "Rongrong: default", "P");
  leg->AddEntry(hMthEffVsBL[1], "Rongrong: MTD acceptance", "P");
  leg->Draw();
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdMthEff/MtdRespEff_Comp.pdf",run_type));
    }
}
