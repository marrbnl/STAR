//================================================
void forZebo()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c = new TCanvas("RHIC_vs_ALICE","RHIC_vs_ALICE",850,750);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  TH1F *h = new TH1F("histogram",";Event activity;#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,4.5);
  h->GetYaxis()->SetRangeUser(0,13.2);
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetLabelFont(62);
  h->GetYaxis()->SetLabelSize(0.035);
  h->DrawCopy();

  TFile *fin = TFile::Open(Form("HP2015.pp.Jpsi.root"),"read");

  // STAR
  TGraphAsymmErrors *ratio = (TGraphAsymmErrors*)fin->Get("STAR_Run13_pp500_Jpsi_mumu_pT0");
  TGraphAsymmErrors *sys = (TGraphAsymmErrors*)fin->Get("STAR_Run13_pp500_Jpsi_mumu_pT0_sys");
  ratio->Draw("sames PEZ");
  sys->Draw("sameE5");

  TGraphErrors *gr = (TGraphErrors*)fin->Get("STAR_Run11_pp500_Jpsi_ee_pT4");
  TGraphAsymmErrors *grSys = (TGraphAsymmErrors*)fin->Get("STAR_Run11_pp500_Jpsi_ee_pT4_sys");
  gr->Draw("sames PEZ");
  grSys->Draw("sameE5");


  // ALICE
  TGraphAsymmErrors *gRatioJ = (TGraphAsymmErrors*)fin->Get("ALICE_pp7000_Jpsi_ee_pT0");
  TGraphAsymmErrors *gSysJ = (TGraphAsymmErrors*)fin->Get("ALICE_pp7000_Jpsi_ee_pT0_sys");
  gRatioJ->Draw("sames PEZ");
  gSysJ->Draw("sameE5");

  TGraphAsymmErrors *gRatioD = (TGraphAsymmErrors*)fin->Get("ALICE_pp7000_D_pT2_4");
  TGraphAsymmErrors *gSysD = (TGraphAsymmErrors*)fin->Get("ALICE_pp7000_D_pT2_4_sys");
  gRatioD->Draw("sames PEZ");
  gSysD->Draw("sameE5");

  leg = new TLegend(0.15,0.65,0.34,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.031);
  leg->SetHeader("p+p collisions");
  leg->AddEntry(ratio,"STAR 500 GeV: J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5, p_{T}>0 GeV/c","P");
  leg->AddEntry(gr,"STAR 500 GeV: J/#psi#rightarrowe^{+}e^{-}, |y|<1, p_{T}>4 GeV/c","P");
  leg->AddEntry(gRatioJ,"ALICE 7 TeV: J/#psi#rightarrowe^{+}e^{-}, |y|<0.9, p_{T}>0 GeV/c","P");
  leg->AddEntry(gRatioD,"ALICE 7 TeV: D meason, |y|<0.5, 2<p_{T}<4 GeV/c","P");
  leg->Draw();

  TPaveText *syserr = new TPaveText(0.38,0.14,0.85,0.2,"brNDC");
  syserr->SetFillStyle(0);
  syserr->SetBorderSize(0);
  syserr->SetTextSize(0.03);
  syserr->SetTextAlign(11);
  syserr->AddText("STAR data points");
  syserr->AddText("+15% one-sided error along both x- and y- direction");
  syserr->SetTextFont(32);
  syserr->SetTextColor(1);
  syserr->Draw();

  TLine *line = new TLine(0,0,4.5,4.5);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  TPaveText *star = new TPaveText(0.28,0.5,0.38,0.55,"brNDC");
  star->SetFillStyle(0);
  star->SetBorderSize(0);
  star->SetTextSize(0.045); 
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
}
