// const int nMultBins = 4;
// const int low_tofMult[nMultBins] = {0,0,4,7};
// const int high_tofMult[nMultBins] = {30,3,6,30};
const int nMultBins = 4;
const int low_tofMult[nMultBins] = {0,1,3,6};
const int high_tofMult[nMultBins] = {30,2,5,30};
const Double_t low_mass = 2.8;
const Double_t high_mass = 3.3;
TString run_cfg_name = "2015HP";

//================================================
void plot_pp500()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  run_type = "Run13_pp500";

  //Run11_xsec();
  evtAct();
  //Run13_signal();
  //compareToModel();
  //Run13_ptDependence();
  //Run11_TofMult();
}

//================================================
void Run11_TofMult(const bool save = 0)
{
  TFile *f = TFile::Open("Rootfiles/BBCvsTOFM.root","read");
  TH2F *hTofMult = (TH2F*)f->Get("hBBCvsTofRefM");
  hTofMult->GetXaxis()->SetRange(150,2000);
  hTofMult->GetYaxis()->SetRangeUser(0,60);
  c = draw2D(hTofMult,";Beam-Beam Counter coincidence rate;TofMult");
  gPad->SetRightMargin(0.16);
  TProfile *pro = (TProfile*)hTofMult->ProfileX("pro");
  pro->SetMarkerStyle(21);
  pro->Draw("sames");
  leg = new TLegend(0.15,0.75,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Run11 p+p @ 500 GeV");
  leg->AddEntry(pro,"Mean of each BBCrate slice","P");
  leg->Draw();
  TPaveText *star = GetPaveText(0.6,0.8,0.8,0.9,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_BBCvsTOFM.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_BBCvsTOFM.png",run_type,run_cfg_name.Data()));
    }
}

//================================================
void compareToModel(const bool save = 1)
{
  // Run11 cross section

  // model
  const double br = 0.0594;
  const double phase = 4*pi;
  ifstream fmodel;
  fmodel.open("Model/nlo_promp_jpsi_yield_500_0_1.dat");
  const int npoints = 36;
  double x[npoints];
  double y[npoints];
  double xe[npoints];
  double yl[npoints];
  double yh[npoints];
  double yel[npoints];
  double yeh[npoints];
  for(int i=0; i<npoints; i++)
    {
      fmodel >> x[i] >> y[i] >> yh[i] >> yl[i];
      double scale = br * 1./phase * 1./x[i];
      y[i] *= scale;
      yh[i] *= scale;
      yl[i] *= scale;
      yel[i] = TMath::Abs(yl[i]-y[i]);
      yeh[i] = TMath::Abs(yh[i]-y[i]);
      xe[i] = 0;
    }

  TCanvas *c2 = new TCanvas("c2","c2", 700, 700);
  c2->SetFillColor(10);
  c2->SetBorderMode(0);
  c2->SetBorderSize(0);
  c2->SetFrameFillColor(10);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderSize(0);
  c2->SetTopMargin(0.03);
  c2->SetRightMargin(0.03);
  c2->SetLeftMargin(0.14);
  c2->SetBottomMargin(0.12);
  gPad->SetLogy();
  TH1F *h = new TH1F("h2","",10,0,25);
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetYaxis()->CenterTitle(1);
  h->GetYaxis()->SetLabelFont(62);
  h->GetYaxis()->SetTitle("B#times1/(2#pip_{T})#timesd^{2}#sigma/(dp_{T}dy)   (nb/GeV/c)^{2}");   
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetRangeUser(1e-6,1);
  h->GetXaxis()->SetRangeUser(0,22.5);
  h->Draw();
  // data
  TFile *fdata = TFile::Open(Form("./Rootfiles/sptrum.root"),"read");	
  TGraphErrors *gData = (TGraphErrors*)fdata->Get("Jpsi_pp500");
  //gData->SetMarkerStyle(29);
  //gData->SetMarkerSize(2);
  gData->Draw("sames PEZ");
  for(int i=0; i<19; i++)
    {
      TBox *box = (TBox*)fdata->Get(Form("sys_uncert_%d",i));
      box->Draw();
    }

  leg = new TLegend(0.5,0.8,0.7,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader("p+p @ 500 GeV");
  leg->AddEntry(gData,"STAR  J/#psi#rightarrowe^{+}e^{-}, |y|<1","P");
  leg->Draw();

  TPaveText *t1 = GetPaveText(0.25,0.5,0.2,0.3,0.04);
  t1->AddText("STAR Preliminary");
  t1->SetTextFont(20);
  t1->SetTextColor(2);
  t1->Draw();

  if(save) 
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data.eps",run_type,run_cfg_name.Data()));
    }

  const int modelColor = kGray+1;
  TGraphAsymmErrors *gSysModel = new TGraphAsymmErrors(npoints,x,y,xe,xe,yel,yeh);
  gSysModel->SetFillStyle(3000);
  gSysModel->SetLineColor(modelColor);
  gSysModel->SetFillColor(modelColor);
  gSysModel->Draw("sames E3");

  gData->Draw("sames PEZ");
  for(int i=0; i<19; i++)
    {
      TBox *box = (TBox*)fdata->Get(Form("sys_uncert_%d",i));
      box->Draw();
    }

  leg = new TLegend(0.5,0.75,0.7,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(gSysModel,"NLO NRQCD: prompt J/#psi","F");
  leg->Draw();



  if(save) 
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data_vs_NRQCD.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data_vs_NRQCD.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data_vs_NRQCD.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_data_vs_NRQCD.eps",run_type,run_cfg_name.Data()));
    }

}

//================================================
void Run13_signal(const bool save = 1)
{
  TFile *f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.EvtAct.root"),"read");
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhJpsiCutStudy_%s",trigName[kTrigType]));
  hn->GetAxis(1)->SetRangeUser(1.5+0.001,100); // pt1 cut
  hn->GetAxis(2)->SetRange(1,1); // unlike-sign
  hn->GetAxis(6)->SetRange(1,350);

  TH1F *hInvMass = (TH1F*)hn->Projection(0);
  hInvMass->SetName(Form("hInvMass"));
  hInvMass->SetMarkerStyle(20);
  hInvMass->SetMarkerColor(2);
  hInvMass->Rebin(2);
  hInvMass->GetXaxis()->SetRangeUser(2.5,3.5);
  TF1 *func = new TF1(Form("func"),"gaus(0)+expo(3)+pol0(5)",2.5,3.5);
  //TF1 *func = new TF1(Form("func"),"gaus(0)+pol3(3)",2,4);
  func->SetParLimits(1,3.05,3.15);
  func->SetParLimits(2,0.04,0.1);
  hInvMass->Fit(func,"IR0Q");
  TF1 *func1 = new TF1(Form("signal"),"gaus",func->GetXmin(),func->GetXmax());
  func1->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
  hInvMass->SetMaximum(1.3*hInvMass->GetMaximum());
  c = draw1D(hInvMass,Form(";M_{#mu#mu} (GeV/c^{2});Counts"));
  func->SetLineColor(4);
  func->Draw("sames");

  // Singal counts
  int low_bin = hInvMass->FindFixBin(low_mass+0.001);
  int high_bin = hInvMass->FindFixBin(high_mass-0.001);
  double nAll = hInvMass->Integral(low_bin,high_bin);
  nAll -= hInvMass->GetBinContent(low_bin) * (low_mass-hInvMass->GetXaxis()->GetBinLowEdge(low_bin))/hInvMass->GetBinWidth(low_bin);
  nAll -= hInvMass->GetBinContent(high_bin) * (hInvMass->GetXaxis()->GetBinUpEdge(high_bin)-high_mass)/hInvMass->GetBinWidth(high_bin);
  double allErr = TMath::Sqrt(nAll);
  double nBkg = (func->Integral(low_mass,high_mass) - func1->Integral(low_mass,high_mass)) * 1./hInvMass->GetBinWidth(1);
  double bkgErr = TMath::Sqrt(nBkg);
  double nSignal = nAll - nBkg;
  double sigErr = TMath::Sqrt(nAll+nBkg);
  TPaveText *signif = GetPaveText(0.13,0.3,0.68,0.88);
  signif->SetTextAlign(11);
  signif->SetTextFont(62);
  signif->AddText("Run13 p+p @ 500 GeV");
  signif->AddText("p_{T,J/#psi} > 0 GeV/c");
  signif->AddText("J/#psi#rightarrow#mu^{+}#mu^{-}");
  //signif->AddText(Form("Signal = %3.1f #pm %3.1f",nSignal,sigErr));
  signif->Draw();

  leg = new TLegend(0.6,0.7,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hInvMass,"Unlike-sign pairs","PL");
  leg->AddEntry(func,"Fit signal+bkg","L");
  leg->Draw();

  TPaveText *star = GetPaveText(0.2,0.3,0.15,0.2,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Fit_Jpsi.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Fit_Jpsi.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Fit_Jpsi.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Fit_Jpsi.eps",run_type,run_cfg_name.Data()));
    }
}

//================================================
void evtAct(const bool save = 1)
{
  TPaveText *syserr = GetPaveText(0.5,0.85,0.15,0.18,0.03);
  syserr->AddText("+15% one-sided error along both x- and y- direction");
  syserr->SetTextFont(32);
  syserr->SetTextColor(1);

  TPaveText *syserr_2 = GetPaveText(0.38,0.85,0.14,0.2,0.03);
  syserr_2->SetTextAlign(11);
  syserr_2->AddText("STAR data points");
  syserr_2->AddText("+15% one-sided error along both x- and y- direction");
  syserr_2->SetTextFont(32);
  syserr_2->SetTextColor(1);

  gStyle->SetEndErrorSize(6);
  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run13.pp500.jpsi.EvtAct.Rank1.VtxDz0.pt0-10.root"),"read");

  const char *draw_style = "sameE5";

  TCanvas *c = new TCanvas("Yield_vs_activity","Yield_vs_activity",850,750);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  TH1F *h = new TH1F("histogram",";#frac{TofMult}{<TofMult>};#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,4.5);
  h->GetYaxis()->SetRangeUser(0,13.2);
  ScaleHistoTitle(h,0.04,1.2,0.035,0.04,1.2,0.035,62);
  h->DrawCopy();

  TGraphAsymmErrors *ratio = (TGraphAsymmErrors*)fin->Get("Run13_Jpsi_vs_evtAct");
  ratio->SetTitle();
  ratio->SetMarkerStyle(29);
  ratio->SetMarkerColor(2);
  ratio->SetLineColor(2);
  ratio->SetMarkerSize(3);
  ratio->GetHistogram()->GetXaxis()->SetRangeUser(0,10);
  ratio->Draw("sames PEZ");
  TGraphAsymmErrors *sys = (TGraphAsymmErrors*)ratio->Clone("gSys");
  for(int i=0; i<sys->GetN(); i++)
    {
      double x,y;
      sys->GetPoint(i,x,y);
      double exl = 0.1;
      //double exl = 0;
      double exh = 0.1;
      double eyl = y * TMath::Sqrt(10*10+5*5+1.5*1.5)/100;
      double eyh = y * TMath::Sqrt(10*10+5*5+1.5*1.5)/100;
      sys->SetPoint(i,x,y);
      sys->SetPointError(i,exl,exh,eyl,eyh);
    }
  sys->SetFillStyle(0);
  sys->SetLineColor(ratio->GetLineColor());
  sys->SetLineWidth(2);
  sys->Draw(draw_style);
  
  // Compare to Qian
  TFile *fHT = TFile::Open("Rootfiles/Spectrum_in_bin_Rank.root");
  TGraphErrors *gr = (TGraphErrors*)fHT->Get("event_Act");
  gr->SetMarkerStyle(29);
  gr->SetMarkerColor(4);
  gr->SetLineColor(4);
  gr->SetMarkerSize(3);
  gr->Draw("sames PEZ");

  TGraphAsymmErrors *grSys = new TGraphAsymmErrors(gr->GetN());
  for(int i=0; i<gr->GetN(); i++)
    {
      double x,y;
      gr->GetPoint(i,x,y);
      //double exl = 0.15*x;
      double exl = 0.1;
      double exh = 0.1;
      double eyl = y * TMath::Sqrt(1.3*1.3*2)/100;
      double eyh = y * TMath::Sqrt(1.3*1.3*2)/100;
      grSys->SetPoint(i,x,y);
      grSys->SetPointError(i,exl,exh,eyl,eyh);
    }
  grSys->SetFillStyle(0);
  grSys->SetLineColor(gr->GetLineColor());
  grSys->SetLineWidth(2);
  grSys->Draw(draw_style);


  leg = new TLegend(0.18,0.7,0.34,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(ratio,"J/#psi#rightarrow#mu^{+}#mu^{-}, p_{T,J/#psi} > 0 GeV/c","PL");
  leg->AddEntry(gr,"J/#psi#rightarrowe^{+}e^{-}, p_{T,J/#psi} > 4 GeV/c","PL");
  leg->Draw();

  TPaveText *t1 = GetPaveText(0.15,0.5,0.88,0.93);
  t1->SetTextFont(62);
  t1->AddText("STAR p+p @ 500 GeV");
  t1->Draw();

  TLine *line = GetLine(0,0,4.5,4.5,1);
  line->Draw();

  TPaveText *star = GetPaveText(0.28,0.38,0.5,0.55,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  syserr->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct.eps",run_type,run_cfg_name.Data()));
    }
 
  //return;
  // Compare to PYTHIA
  TCanvas *c = new TCanvas("RHIC_vs_PYTHIA","RHIC_vs_PYTHIA",850,750);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  h->SetXTitle("Event activity");
  h->DrawCopy();

  TGraphErrors *gPythia[2];
  TFile *fpythia = TFile::Open("Rootfiles/jpsi-pt.default.20150629.root");
  for(int i=0; i<2; i++)
    {
      gPythia[i] = (TGraphErrors*)fpythia->Get(Form("Graph;%d",i+1));
      gPythia[i]->SetMarkerStyle(20);
      gPythia[i]->SetMarkerColor(4-i*2);
      gPythia[i]->SetLineColor(4-i*2);
      gPythia[i]->SetMarkerSize(2);
      gPythia[i]->SetFillStyle(3003);
      gPythia[i]->SetFillColor(4-i*2);
      //if(i==0) gPythia[i]->SetFillColor(kRed-2);
      gPythia[i]->SetLineStyle(2);
      gPythia[i]->SetLineWidth(2);
      gPythia[i]->Draw("samesE3");

      TGraphErrors* copy = (TGraphErrors*)gPythia[i]->Clone(Form("%s_clone",gPythia[i]->GetName()));
      for(int ipoint=0; ipoint<copy->GetN(); ipoint++)
	{
	  copy->SetPointError(ipoint,0,0);
	}
      copy->Draw("sames C");
    }
  leg = new TLegend(0.15,0.65,0.34,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p collisions @ 500 GeV");
  leg->AddEntry(ratio,"STAR: J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5, p_{T}>0 GeV/c","P");
  leg->AddEntry(gr,"STAR: J/#psi#rightarrowe^{+}e^{-}, |y|<1, p_{T}>4 GeV/c","P");
  leg->AddEntry(gPythia[1],"PYTHIA8.183 default: p_{T}>0 GeV/c","FL");
  leg->AddEntry(gPythia[0],"PYTHIA8.183 default: p_{T}>4 GeV/c","FL");
  leg->Draw();
  line->Draw();
  star->Draw();
  syserr_2->Draw();

  ratio->Draw("sames PEZ");
  sys->Draw(draw_style);
  gr->Draw("sames PEZ");
  grSys->Draw(draw_style);

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA.eps",run_type,run_cfg_name.Data()));
    }
  //return;

  TCanvas *c = new TCanvas("RHIC_vs_PYTHIA_Model","RHIC_vs_PYTHIA_Model",850,750);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  h->DrawCopy();
  for(int i=0; i<2; i++)
    {
      gPythia[i]->Draw("sames E3");
      TGraphErrors* copy = (TGraphErrors*)gPythia[i]->Clone(Form("%s_clone",gPythia[i]->GetName()));
      for(int ipoint=0; ipoint<copy->GetN(); ipoint++)
	{
	  copy->SetPointError(ipoint,0,0);
	}
      copy->Draw("sames C");
    }
  leg->Draw();
  line->Draw();
  star->Draw();
  TGraphErrors *hModel = new TGraphErrors(9);
  double model_prediction[9] = {1,2.268,4.140,6.868,10.7,15.87,22.876,31.856,43.632};
  for(int ipoint=0; ipoint<9; ipoint++)
    {
      hModel->SetPoint(ipoint,ipoint+1,model_prediction[ipoint]);
    }
  hModel->SetMarkerStyle(20);
  hModel->SetLineColor(2);
  hModel->SetLineWidth(2);
  hModel->Draw("sames C");
  leg = new TLegend(0.15,0.6,0.34,0.65);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hModel,"Percolation model: p_{T}>0 GeV/c","L");
  leg->Draw();
  syserr_2->Draw();

  ratio->Draw("sames PEZ");
  sys->Draw(draw_style);
  gr->Draw("sames PEZ");
  grSys->Draw(draw_style);

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA_vs_model.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA_vs_model.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA_vs_model.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Data_vs_PYTHIA_vs_model.eps",run_type,run_cfg_name.Data()));
    }
  //return;
  // compare to ALICE data
  TCanvas *c = new TCanvas("RHIC_vs_ALICE","RHIC_vs_ALICE",850,750);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  // TH1F *h = new TH1F("histogram_2",";#frac{Multiplicity}{<Multiplicity>};#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,4.5);
  // h->GetYaxis()->SetRangeUser(0,13);
  // ScaleHistoTitle(h,0.04,1.2,0.035,0.04,1.2,0.035,62);
  h->DrawCopy();
  TGraphAsymmErrors *ratio_clone = (TGraphAsymmErrors*)ratio->Clone("ratio_clone");
  TGraphAsymmErrors *sys_clone = (TGraphAsymmErrors*)sys->Clone("sys_clone");
  TGraphAsymmErrors *grSys_clone = (TGraphAsymmErrors*)grSys->Clone("grSys_clone");
  TGraphErrors *gr_clone = (TGraphErrors*)gr->Clone("grSys_clone");
  //ratio_clone->SetMarkerSize(1.8);
  //gr_clone->SetMarkerSize(1.8);
  ratio_clone->Draw("sames PEZ");
  sys_clone->Draw(draw_style);
  gr_clone->Draw("sames PEZ");
  grSys_clone->Draw(draw_style);

  // J/psi
  const int nj = 5;
  double xj[nj] = {0.44,1.18,1.78,2.63,4.00};
  double exlj[nj] = {0,0,0,0,0};
  double sxlj[nj] = {0.03,0.07,0.11,0.17,0.25};
  double yj[nj] = {0.274, 0.991, 2.272, 3.086, 7.901};
  double eylj[nj] = {0.043,0.166,0.32,0.487,1.397};
  double sylj[nj] = {0.007,0.052,0.275,0.312,0.198};

  TGraphAsymmErrors *gRatioJ = new TGraphAsymmErrors(nj,xj,yj,exlj,exlj,eylj,eylj);
  gRatioJ->SetMarkerStyle(24);
  gRatioJ->SetMarkerSize(2);
  gRatioJ->SetMarkerColor(kGreen+2);
  gRatioJ->SetLineColor(gRatioJ->GetMarkerColor());
  //gRatioJ->SetLineWidth(2);
  gRatioJ->Draw("sames PEZ");
  
  TGraphAsymmErrors *gSysJ = new TGraphAsymmErrors(nj,xj,yj,sxlj,sxlj,sylj,sylj);
  gSysJ->SetFillStyle(0);
  gSysJ->SetLineColor(gRatioJ->GetMarkerColor());
  gSysJ->Draw("sameE5");

  // D-meason
  const int nd = 6;
  double xd[nd] = {0.45,1.18,1.78,2.63,4.01,6.11};
  double exld[nd] = {0,0,0,0,0,0};
  double exhd[nd] = {0,0,0,0,0,0};
  double sxld[nd] = {0.03,0.07,0.11,0.17,0.25,0.39};
  double sxhd[nd] = {0.03,0.07,0.10,0.15,0.23,0.35};
  double yd[nd] = {0.21,1.08,2.28,4.28,9.02,16.51};
  double eyd[nd] = {0.02,0.06,0.12,0.19,0.57,2.60};
  double syd[nd] = {0.01,0.06,0.13,0.17,0.47,2.11};

  TGraphAsymmErrors *gRatioD = new TGraphAsymmErrors(nd,xd,yd,exld,exld,eyd,eyd);
  gRatioD->SetMarkerStyle(20);
  gRatioD->SetMarkerSize(2);
  gRatioD->Draw("sames PEZ");
  
  TGraphAsymmErrors *gSysD = new TGraphAsymmErrors(nd,xd,yd,sxld,sxhd,syd,syd);
  gSysD->SetFillStyle(0);
  gSysD->SetLineColor(gRatioD->GetLineColor());
  gSysD->Draw("sameE5");

  leg = new TLegend(0.15,0.65,0.34,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.031);
  leg->SetHeader("p+p collisions");
  leg->AddEntry(ratio,"STAR 500 GeV: J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5, p_{T}>0 GeV/c","P");
  leg->AddEntry(gr,"STAR 500 GeV: J/#psi#rightarrowe^{+}e^{-}, |y|<1, p_{T}>4 GeV/c","P");
  leg->AddEntry(gRatioJ,"ALICE 7 TeV: J/#psi#rightarrowe^{+}e^{-}, |y|<0.9, p_{T}>0 GeV/c","P");
  leg->AddEntry(gRatioD,"ALICE 7 TeV: D meson, |y|<0.5, 2<p_{T}<4 GeV/c","P");
  leg->Draw();

  TLine *line = GetLine(0,0,4.5,4.5,1);
  line->Draw();
  star->Draw();
  syserr_2->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/RHIC_vs_LHC.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/RHIC_vs_LHC.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/RHIC_vs_LHC.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/RHIC_vs_LHC.eps",run_type,run_cfg_name.Data()));
    }

  //return;

  // compare to more pT bins
  TCanvas *c = new TCanvas("Yield_vs_activity_2","Yield_vs_activity_2",850,750);
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  TH1F *h = new TH1F("histogram_3",";#frac{TofMult}{<TofMult>};#frac{N_{J/#psi}}{<N_{J/#psi}>}",7,0,4);
  h->GetYaxis()->SetRangeUser(0,10);
  ScaleHistoTitle(h,0.04,1.2,0.035,0.04,1.2,0.035,62);
  h->DrawCopy();

  ratio->Draw("sames PEZ");
  sys->Draw(draw_style);
  TFile *fHT2 = TFile::Open("Rootfiles/Spectrum_in_bin_Rank.root");
  TGraphErrors *graph[3];
  graph[0] = (TGraphErrors*)fHT2->Get("event_Act");
  graph[1] = (TGraphErrors*)fHT2->Get("evAct_lowpT");
  graph[2] = (TGraphErrors*)fHT2->Get("evAct_higpT");
  TGraphErrors *graphSys[3];
  for(int i=0; i<3; i++)
    {
      graphSys[i] = (TGraphErrors*)graph[i]->Clone(Form("grSys_%d",i));
      for(int j=0; j<graphSys[i]->GetN(); j++)
	{
	  double x,y;
	  graphSys[i]->GetPoint(j,x,y);
	  double exl = 0.1;
	  double exh = 0.1;
	  double eyl = y * TMath::Sqrt(1.3*1.3*2)/100;
	  double eyh = eyl;
	  graphSys[i]->SetPoint(j,x,y);
	  graphSys[i]->SetPointError(j,exl,eyl);
	}
      graph[i]->SetLineWidth(1);
      if(i==0)
	{
	  graph[i]->SetMarkerStyle(29);
	  graph[i]->SetMarkerColor(4);
	  graph[i]->SetLineColor(4);
	  graph[i]->SetMarkerSize(3);
	}
      if(i==1)
	{
	  graph[i]->SetMarkerStyle(20);
	  graph[i]->SetMarkerColor(kGreen+2);
	  graph[i]->SetLineColor(kGreen+2);
	  graph[i]->SetMarkerSize(2);
	}
      if(i==2)
	{
	  graph[i]->SetMarkerStyle(21);
	  graph[i]->SetMarkerColor(6);
	  graph[i]->SetLineColor(6);
	  graph[i]->SetMarkerSize(2);
	}
      graphSys[i]->SetFillStyle(0);
      graphSys[i]->SetLineColor(graph[i]->GetLineColor());
    }

  for(int i=1; i<3; i++)
    {
      graph[i]->Draw("same PEZ");
      graphSys[i]->Draw(draw_style);
    }
  leg = new TLegend(0.18,0.67,0.34,0.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(ratio,"J/#psi#rightarrow#mu^{+}#mu^{-}, p_{T,J/#psi} > 0 GeV/c","PL");
  leg->AddEntry(graph[1],"J/#psi#rightarrowe^{+}e^{-}, 4 < p_{T,J/#psi} < 8 GeV/c","PL");
  leg->AddEntry(graph[2],"J/#psi#rightarrowe^{+}e^{-}, p_{T,J/#psi} > 8 GeV/c","PL");
  leg->Draw();
  t1->Draw();
  TLine *line = GetLine(0,0,4,4,1);
  line->Draw();
  star->Draw();
  syserr->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct_PtBin.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct_PtBin.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct_PtBin.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Jpsi_vs_EvtAct_PtBin.eps",run_type,run_cfg_name.Data()));
    }

  if(0)
    {
      TFile *fout = TFile::Open("Rootfiles/HP2015.pp.Jpsi.root","recreate");
      ratio->Write("STAR_Run13_pp500_Jpsi_mumu_pT0");
      sys->Write("STAR_Run13_pp500_Jpsi_mumu_pT0_sys");
      gr->Write("STAR_Run11_pp500_Jpsi_ee_pT4");
      grSys->Write("STAR_Run11_pp500_Jpsi_ee_pT4_sys");
      graph[1]->Write("STAR_Run11_pp500_Jpsi_ee_pT4_8");
      graphSys[1]->Write("STAR_Run11_pp500_Jpsi_ee_pT4_8_sys");
      graph[2]->Write("STAR_Run11_pp500_Jpsi_ee_pT8");
      graphSys[2]->Write("STAR_Run11_pp500_Jpsi_ee_pT8_sys");
      gPythia[1]->Write("PYTHIA8.183_pp500_Jpsi_pT0");
      gPythia[0]->Write("PYTHIA8.183_pp500_Jpsi_pT4");
      hModel->Write("PercolationModel_pp500_Jpsi_pT0");
      gRatioJ->Write("ALICE_pp7000_Jpsi_ee_pT0");
      gSysJ->Write("ALICE_pp7000_Jpsi_ee_pT0_sys");
      gRatioD->Write("ALICE_pp7000_D_pT2_4");
      gSysD->Write("ALICE_pp7000_D_pT2_4_sys");
    }
}

//================================================
void Run11_xsec(const bool save = 1)
{
  TFile *fxsec = TFile::Open("./Rootfiles/Spectrum_in_bin_Rank.root","read");
  TFile *fsys  = TFile::Open("./Rootfiles/HTppSpectrum.root","read"); 
  TGraphErrors *gSysAll = (TGraphErrors*)fsys->Get("gSpectrum_sysErr");
  TGraphErrors *gSpectrum[5];
  TGraphErrors *gSys[5];
  const char *name_temp [5] = {"gall","gbin1","gbin2","gbin3","gbin4"};
  const int marker_style[5] = {20,21,24,25,29};
  const int marker_color[5] = {1,2,4,6,kGreen+2};
  const double marker_size[5] = {1.5,1.5,1.5,1.5,2};
  for(int i=0; i<5; i++)
    {
      TGraphErrors *gtmp = (TGraphErrors*)fxsec->Get(name_temp[i]);
      gSpectrum[i] = new TGraphErrors(gtmp->GetN());
      double s,t;
      for(int ipoint=0; ipoint<gSpectrum[i]->GetN(); ipoint++)
	{
	  gtmp->GetPoint(ipoint,s,t);
	  gSpectrum[i]->SetPoint(ipoint,s,t);
	  gSpectrum[i]->SetPointError(ipoint,gtmp->GetErrorX(ipoint),gtmp->GetErrorY(ipoint));
	}
      gSpectrum[i]->SetMarkerStyle(marker_style[i]);
      gSpectrum[i]->SetMarkerColor(marker_color[i]);
      gSpectrum[i]->SetMarkerSize(marker_size[i]);
      gSpectrum[i]->GetYaxis()->SetRangeUser(5e-6,10);
      gSpectrum[i]->SetLineColor(marker_color[i]);

      int nPoints = gSpectrum[i]->GetN();
      gSys[i] = new TGraphErrors(nPoints);
      double x,y;
      double x1,y1,x2,y2;
      double scale;
      for(int ipoint=0; ipoint<nPoints; ipoint++)
	{
	  gSpectrum[i]->GetPoint(ipoint,x,y);
	  gSys[i]->SetPoint(ipoint,x,y);
	  for(int jpoint=0; jpoint<gSysAll->GetN(); jpoint++)
	    {
	      gSysAll->GetPoint(jpoint,x1,y1);
	      gSysAll->GetPoint(jpoint+1, x2, y2);
	      if(x<=x2*mean_pt && x>=x1*mean_pt)
		{
		  scale = (gSysAll->GetErrorY(jpoint)/y1+gSysAll->GetErrorY(jpoint+1)/y2)/2;
		  //cout << x << "  " << x1 << "  " << x2 << endl;
		  //cout << scale << " = " << gSysAll->GetErrorY(jpoint)/y1 << "  " << gSysAll->GetErrorY(jpoint+1)/y2 << endl;
		  break;
		}
	    }
	  //additional 15% uncertainty
	  if(i>0) scale = TMath::Sqrt(scale*scale+0.15*0.15);
	  gSys[i]->SetPointError(ipoint,0.4,y*scale);
	}
      gSys[i]->SetFillStyle(0);
      gSys[i]->SetLineColor(marker_color[i]);
    }
  TCanvas *c = new TCanvas("xsec","xsec",850,750);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.05);
  gSpectrum[0]->SetTitle(";p_{T} [GeV/c];#sigma_{pp}/#sigma_{i} Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]");
  ScaleHistoTitle(gSpectrum[0]->GetHistogram(),0.045,1.1,0.035,0.045,1.3,0.035,62);
  gSpectrum[0]->Draw("AP");
  for(int i=0; i<5; i++)
    {
      gSpectrum[i]->Draw("sames PE");
      gSys[i]->Draw("sames E5");
    }
  leg = new TLegend(0.57,0.64,0.73,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gSpectrum[0],"Inclusive","P");
  leg->AddEntry(gSpectrum[1],"TofMult: 1 - 4","P");
  leg->AddEntry(gSpectrum[2],"TofMult: 5 - 8","P");
  leg->AddEntry(gSpectrum[3],"TofMult: 9 - 12","P");
  leg->AddEntry(gSpectrum[4],"TofMult: >= 13","P");
  leg->Draw();

  TPaveText *t1 = GetPaveText(0.15,0.5,0.88,0.93);
  t1->SetTextFont(62);
  t1->AddText("Run11 p+p @ 500 GeV");
  t1->Draw();
  TPaveText *star = GetPaveText(0.28,0.38,0.18,0.23,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Xsec_TofMultBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Xsec_TofMultBins.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Xsec_TofMultBins.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Xsec_TofMultBins.eps",run_type,run_cfg_name.Data()));
    }

  // Rpp
  TGraphErrors *gRpp[4];
  TGraphErrors *gSysRpp[4];
  for(int i=0; i<4; i++)
    {
      gRpp[i] = (TGraphErrors*)fxsec->Get(Form("Rppbin%d",i+1));
      gRpp[i]->SetMarkerStyle(marker_style[i+1]);
      gRpp[i]->SetMarkerColor(marker_color[i+1]);
      gRpp[i]->SetMarkerSize(marker_size[i+1]);
      gRpp[i]->GetYaxis()->SetRangeUser(0,6);
      gRpp[i]->SetLineColor(marker_color[i+1]);

      int nPoints = gRpp[i]->GetN();
      gSysRpp[i] = new TGraphErrors(nPoints);
      double x,y;
      double error = sqrt(1.3*1.3*2+15*15)/100;
      for(int ipoint=0; ipoint<nPoints; ipoint++)
	{
	  gRpp[i]->GetPoint(ipoint,x,y);
	  gSysRpp[i]->SetPoint(ipoint,x,y);
	  gSysRpp[i]->SetPointError(ipoint,0.4,y*error);
	}
      gSysRpp[i]->SetFillStyle(0);
      gSysRpp[i]->SetLineColor(marker_color[i+1]);
    }
  TCanvas *c = new TCanvas("Rpp","Rpp",850,750);
  SetPadMargin(gPad,0.13,0.14,0.05,0.05);
  gRpp[0]->SetTitle(";p_{T} [GeV/c];#frac{d#sigma/dp_{T}}{<d#sigma/dp_{T}>}");
  ScaleHistoTitle(gRpp[0]->GetHistogram(),0.045,1.1,0.035,0.045,1.1,0.035,62);
  gRpp[0]->Draw("AP");
  for(int i=0; i<4; i++)
    {
      gRpp[i]->Draw("sames PE");
      gSysRpp[i]->Draw("sames E5");
    }
  leg = new TLegend(0.57,0.7,0.73,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gRpp[0],"TofMult: 1 - 4","P");
  leg->AddEntry(gRpp[1],"TofMult: 5 - 8","P");
  leg->AddEntry(gRpp[2],"TofMult: 9 - 12 * 1/2","P");
  leg->AddEntry(gRpp[3],"TofMult: >= 13 * 1/3","P");
  leg->Draw();

  TPaveText *t1 = GetPaveText(0.15,0.5,0.88,0.93);
  t1->SetTextFont(62);
  t1->AddText("Run11 p+p @ 500 GeV");
  t1->Draw();

  TLine *line = GetLine(4,1,13.8,1,1);
  line->Draw();

  TPaveText *star = GetPaveText(0.28,0.38,0.68,0.73,0.045);
  star->AddText("STAR preliminary");
  star->SetTextFont(20);
  star->SetTextColor(2);
  star->Draw();
  
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Rpp_TofMultBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Rpp_TofMultBins.png",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Rpp_TofMultBins.jpg",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run11_Rpp_TofMultBins.eps",run_type,run_cfg_name.Data()));
    }
  
}

//================================================
void Run13_ptDependence(const bool save = 1)
{
  const int marker_style[3] = {21,24,29};
  const int marker_color[3] = {2,4,kGreen+2};
  const double marker_size[3] = {1.5,1.5,2};

  TFile *file[2];
  file[0] = TFile::Open("Rootfiles/Pico.Run13.pp500.jpsi.EvtAct.Rank0.pt0-2.root","read");
  file[1] = TFile::Open("Rootfiles/Pico.Run13.pp500.jpsi.EvtAct.Rank0.pt2-4.root","read");

  TGraphErrors *graph[3];
  for(int i=0; i<3; i++) graph[i] = new TGraphErrors(2);
  double x,y;
  double pt[2] = {1,3};
  double scale[3] = {1,2,3};
  for(int i=0; i<2; i++)
    {
      TGraphAsymmErrors	*gr = (TGraphAsymmErrors*)file[i]->Get("Run13_Jpsi_vs_evtAct");
      for(int ipoint=0; ipoint<gr->GetN(); ipoint++)
	{
	  gr->GetPoint(ipoint,x,y);
	  graph[ipoint]->SetPoint(i,pt[i],y/scale[ipoint]);
	  graph[ipoint]->SetPointError(i,0,gr->GetErrorYlow(ipoint)/scale[ipoint]);
	}
    }
  for(int i=0; i<3; i++)
    {
      graph[i]->SetMarkerStyle(marker_style[i]);
      graph[i]->SetMarkerColor(marker_color[i]);
      graph[i]->SetMarkerSize(marker_size[i]);
      graph[i]->GetYaxis()->SetRangeUser(0,6);
      graph[i]->SetLineColor(marker_color[i]);
    }

  TCanvas *c = new TCanvas("Rpp","Rpp",850,750);
  SetPadMargin(gPad,0.13,0.14,0.05,0.05);
  graph[0]->SetTitle(";p_{T} [GeV/c];#frac{d#sigma/dp_{T}}{<d#sigma/dp_{T}>}");
  graph[0]->GetYaxis()->SetRangeUser(0,6);
  ScaleHistoTitle(graph[0]->GetHistogram(),0.045,1.1,0.035,0.045,1.1,0.035,62);
  graph[0]->Draw("AP");
  graph[1]->Draw("sames PE");
  graph[2]->Draw("sames PE");

  leg = new TLegend(0.57,0.7,0.73,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Run13 dimuon channel");
  leg->AddEntry(graph[0],"TofMult: 1 - 2","P");
  leg->AddEntry(graph[1],"TofMult: 3 - 5 * 1/2","P");
  leg->AddEntry(graph[2],"TofMult: >= 6 * 1/3","P");
  leg->Draw();

  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Rpp_TofMultBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run13_Rpp_TofMultBins.png",run_type,run_cfg_name.Data()));
    }

}