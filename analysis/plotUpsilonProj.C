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
#include "TPaveText.h"
TString run_cfg_name = "2016sQM";

TH2D* histo(TString name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle);
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex);
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle , Int_t lineColor);
void drawLines(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle, Int_t lineColor);
void SetPad(Double_t left, Double_t right, Double_t top, Double_t bottom);
void clearPad(TCanvas *c, Int_t nPads);
void setHisto(TH1D *h,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void setGraph(TGraphErrors *gr,Int_t MarkerStyle, Double_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void setLegend(TLegend *leg,Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t textSize);
void pdfAction(TCanvas *c, TPDF *ps);

void plotUpsilonProj()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptDate(0);

  TGraphErrors *grPPStat = new TGraphErrors(1);
  grPPStat->SetPoint(0,100,0.421);
  grPPStat->SetPointError(0,0,0.009);
  setGraph(grPPStat,22,2,1,1);

  TGraphErrors *grPPSys = new TGraphErrors(1);
  grPPSys->SetPoint(0,100,0.421);
  grPPSys->SetPointError(0,2.2,0.017);
  setGraph(grPPSys,22,1.6,1,1);
  grPPSys->SetFillStyle(1001);
  grPPSys->SetFillColor(kGray+1);

  TGraphErrors *grPbPbStat = new TGraphErrors(1);
  grPbPbStat->SetPoint(0,215,0.09);
  grPbPbStat->SetPointError(0,0,0.02);
  setGraph(grPbPbStat,20,1.3,4,4);
  grPbPbStat->SetMarkerSize(1.5);

  TGraphErrors *grPbPbSys = new TGraphErrors(1);
  grPbPbSys->SetPoint(0,215,0.09);
  grPbPbSys->SetPointError(0,2.1,0.02);
  setGraph(grPbPbSys,20,1.3,4,4);
  grPbPbSys->SetFillStyle(1001);
  grPbPbSys->SetFillColor(kGray+1);

  TGraphAsymmErrors *grPbPb3s0 = new TGraphAsymmErrors(1);
  grPbPb3s0->SetPoint(0,235,0.04);
  grPbPb3s0->SetPointError(0,0,0,0.04,0);
  grPbPb3s0->SetMarkerStyle(20);
  grPbPb3s0->SetMarkerSize(1.5);
  grPbPb3s0->SetMarkerColor(4);
  grPbPb3s0->SetLineColor(4);
  grPbPb3s0->SetFillColor(4);
  grPbPb3s0->SetLineWidth(2);

  TGraphAsymmErrors *grPbPb3s1 = new TGraphAsymmErrors(1);
  grPbPb3s1->SetPoint(0,235,0.04);
  grPbPb3s1->SetPointError(0,10,10,0,0);
  grPbPb3s1->SetMarkerStyle(29);
  grPbPb3s1->SetMarkerSize(0);
  grPbPb3s1->SetMarkerColor(4);
  grPbPb3s1->SetLineColor(4);
  grPbPb3s1->SetLineWidth(2);

  TArrow *PbPb3sArr = new TArrow(235,0,235,0.04,0.05,"<|");
  PbPb3sArr->SetAngle(40);
  PbPb3sArr->SetArrowSize(0.02);
  PbPb3sArr->SetLineWidth(2);
  PbPb3sArr->SetLineColor(4);
  PbPb3sArr->SetFillColor(4);


  TGraphAsymmErrors *grAuAuee0 = new TGraphAsymmErrors(1);
  grAuAuee0->SetPoint(0,180,0.217);
  grAuAuee0->SetPointError(0,0,0,0.14,0);
  grAuAuee0->SetMarkerStyle(21);
  grAuAuee0->SetMarkerSize(1.6);
  grAuAuee0->SetMarkerColor(6);
  grAuAuee0->SetLineColor(6);
  grAuAuee0->SetFillColor(6);
  grAuAuee0->SetLineWidth(2);

  TArrow *AuAuee = new TArrow(180,0.1,180,0.217,0.05,"<|");
  AuAuee->SetAngle(40);
  AuAuee->SetArrowSize(0.03);
  AuAuee->SetLineWidth(2);
  AuAuee->SetLineColor(6);
  AuAuee->SetFillColor(6);

  TGraphAsymmErrors *grAuAuee1 = new TGraphAsymmErrors(1);
  grAuAuee1->SetPoint(0,180,0.217);
  grAuAuee1->SetPointError(0,10,10,0,0);
  grAuAuee1->SetMarkerStyle(29);
  grAuAuee1->SetMarkerSize(0);
  grAuAuee1->SetMarkerColor(6);
  grAuAuee1->SetLineColor(6);
  grAuAuee1->SetLineWidth(2);

  TGraphErrors *grAuAumumu = new TGraphErrors(1);
  grAuAumumu->SetPoint(0,155,0.375);
  grAuAumumu->SetPointError(0,0,0.132);
  setGraph(grAuAumumu,29,2.5,2,2);

  TGraphErrors *grAuAuSys = new TGraphErrors(1);
  grAuAuSys->SetPoint(0,155,0.375);
  grAuAuSys->SetPointError(0,2.3,0.083);
  setGraph(grAuAuSys,29,2.5,2,2);
  grAuAuSys->SetFillStyle(1001);
  grAuAuSys->SetFillColor(kGray+1);


  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  TLegend *leg = new TLegend();
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(22);
  leg->SetTextSize(0.04);

  TH2D *dd = (TH2D *)histo("dd",80,260,0,0.8,"","#varUpsilon(2S+3S)/#varUpsilon(1S)");
  c1->cd();
  SetPad(0.13,0.05,0.05,0.13);
  gPad->SetGrid(0,0);
  dd->GetXaxis()->SetLabelSize(0);
  dd->GetXaxis()->SetNdivisions(205);
  dd->GetXaxis()->SetTickLength(0);
  dd->DrawCopy("c");

  TLine *line1 = new TLine(120,0,120,0.8);
  line1->SetLineStyle(2);
  line1->Draw();

  TLine *line2 = new TLine(190,0,190,0.8);
  line2->SetLineStyle(2);
  line2->Draw();

  setLegend(leg,0.2,0.75,0.7,0.949,0.04);
  leg->AddEntry(grPPStat,"p+p (world-wide)","pl");
  leg->AddEntry(grPbPbStat,"CMS Pb+Pb@2.76 TeV (0-100%)","pl");
  //leg->AddEntry(grAuAuee0,"STAR Au+Au@200 GeV (ee) (0-60%)","l");
  leg->AddEntry(grAuAumumu,"STAR Au+Au@200 GeV (#mu#mu) (0-80%)","pl");

  grAuAuSys->Draw("e2same");
  grAuAumumu->Draw("pesamez");
  grPPSys->Draw("e2same");
  grPPStat->Draw("psame");
  grPbPbSys->Draw("e2same");
  grPbPbStat->Draw("psame");
  //grAuAuee0->Draw("|>same");
  //AuAuee->Draw();
  //grAuAuee1->Draw("pzesame");
  //grPbPb3s0->Draw("|>pesame");
  PbPb3sArr->Draw();
  grPbPb3s1->Draw("pzesame");
  leg->Draw("same");
  TPaveText *star = new TPaveText(0.75,0.6,0.85,0.7,"brNDC");
  star->SetFillStyle(0);
  star->SetBorderSize(0);
  star->AddText("STAR preliminary");
  star->SetTextFont(22);
  star->SetTextSize(0.045);
  star->SetTextColor(2);
  star->Draw("");
  
  TPaveText *cms = new TPaveText(0.7,0.24,0.8,0.29,"brNDC");
  cms->AddText("#varUpsilon(2S)/#varUpsilon(1S)");
  cms->SetTextFont(22);
  cms->SetTextColor(4);
  cms->SetTextSize(.034);
  cms->SetFillColor(0);
  cms->SetFillStyle(0);
  cms->SetBorderSize(0);
  cms->Draw("");

  TPaveText *cms2 = new TPaveText(0.79,0.17,0.89,0.22,"brNDC");
  cms2->AddText("#varUpsilon(3S)/#varUpsilon(1S)");
  cms2->SetTextFont(22);
  cms2->SetTextColor(4);
  cms2->SetTextSize(.034);
  cms2->SetFillColor(0);
  cms2->SetFillStyle(0);
  cms2->SetBorderSize(0);
  cms2->Draw("");

  TPaveText *pp = new TPaveText(0.18,0.08,0.23,0.1,"brNDC");
  pp->AddText("pp");
  pp->SetTextFont(22);
  pp->SetTextSize(.045);
  pp->SetFillColor(0);
  pp->SetFillStyle(0);
  pp->SetBorderSize(0);
  pp->Draw("");

  TPaveText *rhic = new TPaveText(0.36,0.08,0.56,0.1,"brNDC");
  rhic->AddText("AuAu @ RHIC");
  rhic->SetTextFont(22);
  rhic->SetTextSize(.045);
  rhic->SetFillColor(0);
  rhic->SetFillStyle(0);
  rhic->SetBorderSize(0);
  rhic->Draw("");

  TPaveText *lhc = new TPaveText(0.75,0.08,0.85,0.1,"brNDC");
  lhc->AddText("PbPb @ LHC");
  lhc->SetTextFont(22);
  lhc->SetTextSize(.045);
  lhc->SetFillColor(0);
  lhc->SetFillStyle(0);
  lhc->SetBorderSize(0);
  lhc->Draw("");

  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/2016sQM/Run14_Upsilon_Ratio_AtStar.pdf"));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/2016sQM/Run14_Upsilon_Ratio_AtStar.png"));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/2016sQM/Run14_Upsilon_Ratio_AtStar.gif"));
  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/2016sQM/Run14_Upsilon_Ratio_AtStar.eps"));

  return;

  TCanvas *c1 = new TCanvas("c2","c2",800,600);
  c1->SetLogy();
  TH2D *dd = (TH2D *)histo("dd2",80,260,1e-2,50,"Collision System","#varUpsilon(2S+3S)/#varUpsilon(1S)");
  SetPad(0.13,0.13,0.05,0.05);
  gPad->SetGrid(0,0);
  dd->GetXaxis()->SetLabelSize(0);
  dd->GetXaxis()->SetNdivisions(205);
  dd->GetXaxis()->SetTickLength(0);
  dd->DrawCopy("c");

  grAuAuSys->Draw("e2same");
  grAuAumumu->Draw("pesamez");
  //grAuAumumu2->Draw("pesamez");
  grPPSys->Draw("e2same");
  grPPStat->Draw("psame");
  grPbPbSys->Draw("e2same");
  grPbPbStat->Draw("psame");
  grAuAuee0->Draw("|>pesame");
  grAuAuee1->Draw("pzesame");
  grPbPb3s0->Draw("|>pesame");
  grPbPb3s1->Draw("pzesame");
  leg->Draw("same");
  star->Draw();
  cms->Draw();
  cms2->Draw();
  //line->Draw();
/*
  c1->SaveAs(Form("Run14_Upsilon_Ratio_log.pdf"));
  c1->SaveAs(Form("Run14_Upsilon_Ratio_log.png"));
  c1->SaveAs(Form("Run14_Upsilon_Ratio_log.jpg"));
  c1->SaveAs(Form("Run14_Upsilon_Ratio_log.eps"));
*/
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
void SetPad(Double_t left, Double_t right, Double_t top, Double_t bottom)
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
