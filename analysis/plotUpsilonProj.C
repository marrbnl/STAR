#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TPDF.h"

TString run_cfg_name = "2015QM";

TH2D* histo(TString name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle);
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex);
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle , Int_t lineColor);
void drawLines(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle, Int_t lineColor);
void setpad(Double_t left, Double_t right, Double_t top, Double_t bottom);
void clearPad(TCanvas *c, Int_t nPads);
void setHisto(TH1D *h,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void setGraph(TGraphErrors *gr,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void setLegend(TLegend *leg,Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t textSize);
void pdfAction(TCanvas *c, TPDF *ps);

void plotUpsilonProj()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);


  TGraphErrors *grPPStat = new TGraphErrors(1);
  grPPStat->SetPoint(0,100,0.421);
  grPPStat->SetPointError(0,0,0.009);
  setGraph(grPPStat,20,1.5,1,1);

  TGraphErrors *grPPSys = new TGraphErrors(1);
  grPPSys->SetPoint(0,100,0.421);
  grPPSys->SetPointError(0,2,0.017);
  setGraph(grPPSys,20,1.5,1,1);
  grPPSys->SetFillStyle(1001);
  grPPSys->SetFillColor(kGray+1);

  TGraphErrors *grPbPbStat = new TGraphErrors(1);
  grPbPbStat->SetPoint(0,220,0.069);
  grPbPbStat->SetPointError(0,0,0.0248);
  setGraph(grPbPbStat,22,2,4,4);

  TGraphErrors *grPbPbSys = new TGraphErrors(1);
  grPbPbSys->SetPoint(0,220,0.069);
  grPbPbSys->SetPointError(0,2,0.0137);
  setGraph(grPbPbSys,22,1.3,4,4);
  grPbPbSys->SetFillStyle(1001);
  grPbPbSys->SetFillColor(kGray+1);

  TGraphAsymmErrors *grAuAuee0 = new TGraphAsymmErrors(1);
  grAuAuee0->SetPoint(0,180,0.2785);
  grAuAuee0->SetPointError(0,0,0,0.15,0);
  grAuAuee0->SetMarkerStyle(21);
  grAuAuee0->SetMarkerSize(1.5);
  grAuAuee0->SetMarkerColor(6);
  grAuAuee0->SetLineColor(6);
  grAuAuee0->SetFillColor(6);
  grAuAuee0->SetLineWidth(2);

  TGraphAsymmErrors *grAuAuee1 = new TGraphAsymmErrors(1);
  grAuAuee1->SetPoint(0,180,0.2785);
  grAuAuee1->SetPointError(0,10,10,0,0);
  grAuAuee1->SetMarkerStyle(29);
  grAuAuee1->SetMarkerSize(1.6);
  grAuAuee1->SetMarkerColor(6);
  grAuAuee1->SetLineColor(6);
  grAuAuee1->SetLineWidth(2);

  TGraphErrors *grAuAumumu = new TGraphErrors(1);
  grAuAumumu->SetPoint(0,125,0.3506);
  grAuAumumu->SetPointError(0,0,0.3382);
  setGraph(grAuAumumu,29,2.5,2,2);

  TGraphErrors *grAuAumumu2 = new TGraphErrors(1);
  grAuAumumu2->SetPoint(0,130,0.069);
  grAuAumumu2->SetPointError(0,0,0.069* 0.3382/sqrt(14)/0.3506 * TMath::Sqrt(0.3506/0.069));
  setGraph(grAuAumumu2,29,2.5,1,1);

  TGraphErrors *grAuAumumu3 = new TGraphErrors(1);
  grAuAumumu3->SetPoint(0,130,0.3506);
  grAuAumumu3->SetPointError(0,0,0.3382/sqrt(14));
  setGraph(grAuAumumu3,29,2.5,1,1);


  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  TLegend *leg = new TLegend();
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->SetTextFont(22);
  leg->SetTextSize(0.04);

  TH2D *dd = (TH2D *)histo("dd",80,240,0,1.2,"Collision System","#varUpsilon(2S+3S)/#varUpsilon(1S)");
  c1->cd();
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  //gPad->SetGridy(1);
  //gStyle->SetGridColor(1);
  dd->DrawCopy("c");
  dd->GetXaxis()->SetLabelSize(0);
  dd->GetXaxis()->SetNdivisions(205);
  setLegend(leg,0.2,0.65,0.6,0.93,0.042);
  leg->AddEntry(grPPStat,"p+p (world-wide)","pl");
  leg->AddEntry(grPbPbStat,"CMS Pb+Pb@2.76 TeV (0-100%)","pl");
  leg->AddEntry(grAuAuee0,"STAR Au+Au@200 GeV (ee) (0-80%)","pl");
  leg->AddEntry(grAuAumumu,"STAR Au+Au@200 GeV (#mu#mu) (0-80%)","pl");
  leg->AddEntry(grAuAumumu2,"STAR #Upsilon#rightarrow#mu#mu Run14+16 projection","pl");
  grAuAumumu->Draw("pesamez");
  grAuAumumu2->Draw("pesamez");
  grPPSys->Draw("e2same");
  grPPStat->Draw("psame");
  grPbPbSys->Draw("e2same");
  grPbPbStat->Draw("psame");
  grAuAuee0->Draw("|>pesame");
  grAuAuee1->Draw("pzesame");
  leg->Draw("same");
  TPaveText *star = GetPaveText(0.7,0.85,0.5,0.55,0.05);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  TLine *line = GetLine(120,0.069,240,0.069,1);
  line->Draw();

  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtCMS.pdf",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtCMS.png",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtCMS.jpg",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtCMS.eps",run_type,run_cfg_name.Data()));

  TCanvas *c1 = new TCanvas("c3","c3",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  dd->DrawCopy("c");
  grAuAumumu->Draw("pesamez");
  grAuAumumu3->Draw("pesamez");
  grPPSys->Draw("e2same");
  grPPStat->Draw("psame");
  grPbPbSys->Draw("e2same");
  grPbPbStat->Draw("psame");
  grAuAuee0->Draw("|>pesame");
  grAuAuee1->Draw("pzesame");
  leg->Draw("same");
  star->Draw();

  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtStar.pdf",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtStar.png",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtStar.jpg",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_AtStar.eps",run_type,run_cfg_name.Data()));


  TCanvas *c1 = new TCanvas("c2","c2",800,600);
  c1->SetLogy();
  TH2D *dd = (TH2D *)histo("dd2",80,240,1e-2,50,"Collision System","#varUpsilon(2S+3S)/#varUpsilon(1S)");
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  dd->GetXaxis()->SetLabelSize(0);
  dd->GetXaxis()->SetNdivisions(205);
  dd->DrawCopy("c");
  grAuAumumu->Draw("pesamez");
  grAuAumumu2->Draw("pesamez");
  grPPSys->Draw("e2same");
  grPPStat->Draw("psame");
  grPbPbSys->Draw("e2same");
  grPbPbStat->Draw("psame");
  grAuAuee0->Draw("|>pesame");
  grAuAuee1->Draw("pzesame");
  leg->Draw("same");
  star->Draw();
  line->Draw();

  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_log.pdf",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_log.png",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_log.jpg",run_type,run_cfg_name.Data()));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14_Upsilon_Ratio_log.eps",run_type,run_cfg_name.Data()));

  cout<<"The program has completed!!"<<endl;
}
//__________________________________________________________________________
TH2D* histo(TString name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle)
{
  TH2D *dd = new TH2D(name.Data(),"",100,xlow,xup,100,ylow,yup);
  dd->GetXaxis()->SetTitle(xTitle.Data());
  dd->GetYaxis()->SetTitle(yTitle.Data());

  dd->GetXaxis()->SetTitleSize(0.055);
  dd->GetXaxis()->SetTitleOffset(0.9);
  dd->GetXaxis()->SetLabelSize(0.045);
  dd->GetYaxis()->SetTitleSize(0.055);
  dd->GetYaxis()->SetTitleOffset(1);
  dd->GetYaxis()->SetLabelSize(0.045);
  dd->GetXaxis()->CenterTitle(kTRUE);
  dd->GetYaxis()->CenterTitle(kTRUE);
  dd->GetXaxis()->SetNdivisions(512);
  return dd;
}
//__________________________________________________________________________
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
  TLatex *latex = new TLatex(x,y,text.Data());
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");
  return latex;
}
//__________________________________________________________________________
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
  TLine *l1 = new TLine(xlow,ylow,xup,yup);
  l1->SetLineWidth(lineWidth);
  l1->SetLineColor(lineColor);
  l1->SetLineStyle(lineStyle);
  l1->Draw("same");
  return l1;
}
//__________________________________________________________________________
void drawLines(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
  drawLine(xlow,ylow,xup,ylow,lineWidth,lineStyle,lineColor);
  drawLine(xlow,yup,xup,yup,lineWidth,lineStyle,lineColor);
  drawLine(xlow,ylow,xlow,yup,lineWidth,lineStyle,lineColor);
  drawLine(xup,ylow,xup,yup,lineWidth,lineStyle,lineColor);
}
//__________________________________________________________________________
void setpad(Double_t left, Double_t right, Double_t top, Double_t bottom)
{
  gPad->SetFillColor(10);
  gPad->SetBorderMode(0);
  gPad->SetBorderSize(0);
  gPad->SetFrameFillColor(10);
  gPad->SetFrameBorderMode(0);
  gPad->SetFrameBorderSize(0);
  gPad->SetLeftMargin(left);
  gPad->SetRightMargin(right);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
}
//__________________________________________________________________________
void clearPad(TCanvas *c, Int_t nPads)
{
  for(Int_t i=0;i<nPads;i++){
    c->cd(i+1);
    gPad->Clear();
  }
}
//__________________________________________________________________________
void setHisto(TH1D *h,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor){
  h->SetMarkerStyle(MarkerStyle);
  h->SetMarkerSize(MarkerSize);
  h->SetMarkerColor(MarkerColor);
  h->SetLineColor(LineColor);
};
//__________________________________________________________________________
void setGraph(TGraphErrors *gr,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor){
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerSize(MarkerSize);
  gr->SetMarkerColor(MarkerColor);
  gr->SetLineColor(LineColor);
};
//__________________________________________________________________________
void setLegend(TLegend *leg, Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t textSize){
  leg->Clear();
  leg->SetX1NDC(xlow);
  leg->SetY1NDC(ylow);
  leg->SetX2NDC(xup);
  leg->SetY2NDC(yup);
  leg->SetTextSize(textSize);
}
//__________________________________________________________________________
void pdfAction(TCanvas *c, TPDF *ps){
  ps->On();
  c->Update();
  c->cd();
  ps->NewPage();
  ps->Off();
};
