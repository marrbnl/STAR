const char *run_type = "Run14_AuAu200";
TString run_cfg_name = "Projection";
const double luminosity = 3.8;
const double ppInelastic = 42.; // mb
const double ppInelasticErr = 3.; // mb
const int nCentBins = 4;
const double ncoll[nCentBins] = {393., 785., 300., 95.}; 
const double ncollErr[nCentBins] = {27., 29., 31., 21.};
const double npart[nCentBins] = {161, 280, 142, 62};
const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60"};
const char *cent_Title[nCentBins] = {"0060","0020","2040","4060"};

//================================================
void projection()
{  
  //jpsi();
  //upsilon();
  efficiency();
}

//================================================
void efficiency(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  TFile *fmuon = TFile::Open("Rootfiles/Run14.AuAu200.JpsiEff.pt1.0.pt1.0.root","read");
  TH1F *hJpsiPt[2], *hJpsiPtRebin[2];
  hJpsiPt[0] = (TH1F*)fmuon->Get("MCinput_Jpsi_pT_0020_WeightPt");
  hJpsiPt[1] = (TH1F*)fmuon->Get("MTDreco_Jpsi_pT_0020_WeightPt");
  const int nbins = 10;
  double xbins[nbins+1] = {0,1,2,3,4,5,6,7,8,9,10};
  for(int i=0; i<2; i++)
    {
      hJpsiPtRebin[i] = (TH1F*)hJpsiPt[i]->Rebin(nbins,Form("%s_rebin",hJpsiPt[i]->GetName()),xbins);
    }
  TH1F *hJpsiEffMuon = (TH1F*)hJpsiPtRebin[1]->Clone("hJpsiEffMuon");
  hJpsiEffMuon->Divide(hJpsiPtRebin[0]);
  hJpsiEffMuon->Scale(0.8);

  TFile *ftrig = TFile::Open("Rootfiles/Run14.AuAu200.JpsiTrigEff.pt1.0.pt1.0.root","read");
  TH1F *htrig = (TH1F*)ftrig->Get("JpsiTrigEff_cent0020_rebin");
  TH1F *hresp = (TH1F*)ftrig->Get("JpsiRespEff_cent0020_rebin");
  hJpsiEffMuon->Multiply(htrig);
  hJpsiEffMuon->Multiply(hresp);

  const int nbins = 6;
  const double xbins[nbins+1] = {3,4,5,6,7,8,10};
  const double eff[nbins] = {1.05e-3, 9e-3, 2.4e-2, 4.2e-2, 6.2e-2, 9e-2};
  TH1F *hJpsiEffElec = new TH1F("JpsiEff_ee","",nbins,xbins);
  for(int ibin=1; ibin<=nbins; ibin++)
    {
      hJpsiEffElec->SetBinContent(ibin,eff[ibin-1]);
      hJpsiEffElec->SetBinError(ibin,1e-10);
    }
  hJpsiEffElec->SetMarkerStyle(20);
  hJpsiEffElec->SetMarkerColor(4);
  hJpsiEffElec->SetLineColor(4);
  hJpsiEffElec->SetMarkerSize(1.5);
  

  TCanvas *c1 = new TCanvas("Jpsi_eff","Jpsi_eff",800,650);
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hJpsiEffMuon,0.05,1,0.04,0.05,1.1,0.04,62);
  hJpsiEffMuon->GetYaxis()->SetRangeUser(1e-3,1);
  hJpsiEffMuon->GetYaxis()->SetNdivisions(505);
  hJpsiEffMuon->SetTitle(";p_{T} (GeV/c);J/#psi efficiency#timesacceptance");
  hJpsiEffMuon->SetMarkerStyle(21);
  hJpsiEffMuon->SetMarkerColor(2);
  hJpsiEffMuon->SetLineColor(2);
  hJpsiEffMuon->SetMarkerSize(1.5);
  hJpsiEffMuon->Draw();
  hJpsiEffElec->Draw("samesP");

  TLegend *leg = new TLegend(0.2,0.7,0.4,0.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(26);
  leg->SetHeader("AuAu @ 200 GeV");
  leg->AddEntry(hJpsiEffElec,"Electron: BEMC trigger with E_{T} > 4.3 GeV","PL");
  leg->AddEntry(hJpsiEffMuon,"Muon: MTD trigger","PL");
  leg->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/CompJpsiEff.pdf",run_type,run_cfg_name.Data()));
    }
}

//================================================
void jpsi(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  const int marker_style[nCentBins] = {kFullCircle, kFullStar, kFullSquare, kFullCross};
  const double marker_size[nCentBins] = {1.5,2,1.5,2};
  const int marker_color[nCentBins] = {1,2,4,6};
  const int pub_marker_style[nCentBins] = {kOpenCircle, kOpenStar, kOpenSquare, kOpenCross};
  const int pub_marker_color_low[nCentBins] = {kGray+2, kRed+2, kBlue+2, kMagenta+2};
  const int pub_marker_color_high[nCentBins] = {1,2,4,6};
  const double scale_factor[nCentBins] = {10,1,0.2,0.1};
  const int lowpt_color = 9;
  const int highpt_color = 12;

  TFile *fpub = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");
  TH1F *hJpsipp = (TH1F*)fpub->Get("Jpsi_xsec_pp200_combined");

  // Publication
  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  TH1F *hTBW[nCentBins];
  TGraphAsymmErrors *gRaaLowPt[nCentBins];
  TGraphAsymmErrors *gRaaLowPtSys[nCentBins];
  TGraphAsymmErrors *gRaaHighPt[nCentBins];
  TGraphAsymmErrors *gRaaHighPtSys[nCentBins];
  TH1F *hJpsipp = (TH1F*)fpub->Get("Jpsi_xsec_pp200_combined");
  double x,y;
  double x1,y1;
  int new_color = kGreen+2;
  for(int k=0; k<nCentBins; k++)
    {
      gAuAuLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[k]));
      gAuAuLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[k]));
      gAuAuHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[k]));
      gAuAuHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[k]));
      hTBW[k] = (TH1F*)fpub->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[k]));
      // scaleGraph(gAuAuLowPt[k], scale_factor[k]);
      // scaleGraph(gAuAuLowPtSys[k], scale_factor[k]);
      // scaleGraph(gAuAuHighPt[k], scale_factor[k]);
      // scaleGraph(gAuAuHighPtSys[k], scale_factor[k]);
      // hTBW[k]->Scale(scale_factor[k]);

      gRaaLowPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",cent_Title[k]));
      gRaaLowPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",cent_Title[k]));
      gRaaHighPt[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",cent_Title[k]));
      gRaaHighPtSys[k] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",cent_Title[k]));

      // make AuAu error only
      Int_t npoints = gRaaLowPt[k]->GetN();
      for (Int_t j=0; j<npoints; j++)
	{ 
	  gRaaLowPt[k]->GetPoint(j,x,y);
	  gAuAuLowPt[k]->GetPoint(j,x1,y1);
	  gRaaLowPt[k]->SetPointError(j, 0, 0, gAuAuLowPt[k]->GetErrorYlow(j)/y1*y, gAuAuLowPt[k]->GetErrorYhigh(j)/y1*y);
	  //cout << gRaaLowPt[k]->GetErrorYlow(j)/y << " =? " << gAuAuLowPt[k]->GetErrorYlow(j)/y1 << endl;
	}

      Int_t npoints = gRaaHighPt[k]->GetN();
      for (Int_t j=0; j<npoints; j++)
	{ 
	  gRaaHighPt[k]->GetPoint(j,x,y);
	  gAuAuHighPt[k]->GetPoint(j,x1,y1);
	  gRaaHighPt[k]->SetPointError(j, 0, 0, gAuAuHighPt[k]->GetErrorYlow(j)/y1*y, gAuAuHighPt[k]->GetErrorYhigh(j)/y1*y);
	}

      
      gRaaLowPt[k]->SetMarkerStyle(kOpenStar);
      gRaaLowPt[k]->SetMarkerSize(2);
      gRaaLowPt[k]->SetMarkerColor(new_color);
      gRaaLowPt[k]->SetLineColor(new_color);
      gRaaLowPtSys[k]->SetMarkerColor(new_color);
      gRaaLowPtSys[k]->SetLineColor(new_color);

      gRaaHighPt[k]->SetMarkerStyle(29);
      gRaaHighPt[k]->SetMarkerSize(2);
      gRaaHighPt[k]->SetMarkerColor(new_color);
      gRaaHighPt[k]->SetLineColor(new_color);
      gRaaHighPtSys[k]->SetMarkerColor(new_color);
      gRaaHighPtSys[k]->SetLineColor(new_color);
      
      offset_x_with_asym_error(gRaaLowPt[k], -0.4);
      offset_x_with_asym_error(gRaaLowPtSys[k], -0.4);
      offset_x_with_asym_error(gRaaHighPt[k], -0.4);
      offset_x_with_asym_error(gRaaHighPtSys[k], -0.4);

      
      Int_t npoints = gRaaLowPtSys[k]->GetN();
      for (Int_t j=0; j<npoints; j++)
	{
	  gRaaLowPtSys[k]->SetPointError(j, 0.1, 0.1, gRaaLowPtSys[k]->GetErrorYlow(j), gRaaLowPtSys[k]->GetErrorYhigh(j));
	}

      Int_t npoints = gRaaHighPtSys[k]->GetN();
      for (Int_t j=0; j<npoints; j++)
	{
	  gRaaHighPtSys[k]->SetPointError(j, 0.1, 0.1, gRaaHighPtSys[k]->GetErrorYlow(j), gRaaHighPtSys[k]->GetErrorYhigh(j));
	}
      
      
    }

  // data
  TFile *fSys = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hSys = (TH1F*)fSys->Get(Form("Sys_all_%s",cent_Title[0]));
  TFile *fdata = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiInvYield[nCentBins];
  TGraphErrors *hJpsiXsec[nCentBins];
  TGraphErrors *hJpsiXsecSys[nCentBins];
  TGraphErrors *hJpsiRaa[nCentBins];
  TGraphErrors *hJpsiRaaSys[nCentBins];
  double x, y;
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)fdata->Get(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsec[k] = new TGraphErrors(hJpsiInvYield[k]);
      hJpsiXsec[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiXsecSys[k] = new TGraphErrors(hJpsiInvYield[k]);
      hJpsiXsecSys[k]->SetName(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[k]));
      for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	{
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  hJpsiXsec[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsec[k]->SetPointError(i,0,hJpsiXsec[k]->GetErrorY(i)*scale_factor[k]);

	  hJpsiXsecSys[k]->SetPoint(i,x,y*scale_factor[k]);
	  hJpsiXsecSys[k]->SetPointError(i,0.2,hSys->GetBinContent(i+1)*y*scale_factor[k]);
	}

      hJpsiRaa[k] = new TGraphErrors();
      hJpsiRaa[k]->SetName(Form("Graph_Jpsi_Raa_cent%s",cent_Title[k]));
      hJpsiRaaSys[k] = new TGraphErrors();
      hJpsiRaaSys[k]->SetName(Form("Graph_Jpsi_Raa_cent%s_sys",cent_Title[k]));
      for(int i=0; i<hJpsiXsec[k]->GetN(); i++)
	{
	  hJpsiXsec[k]->GetPoint(i,x,y);
	  double AuAu_val = y/scale_factor[k];
	  double AuAu_err = hJpsiXsec[k]->GetErrorY(i)/scale_factor[k];
	  double pp_val = hJpsipp->GetBinContent(i+1);
	  double pp_err = hJpsipp->GetBinError(i+1);
	  double prefix = ppInelastic/ncoll[k] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  //double err = val * sqrt(AuAu_err*AuAu_err/AuAu_val/AuAu_val + pp_err*pp_err/pp_val/pp_val);
	  double err = prefix * AuAu_err / pp_val; // AuAu error only
	  hJpsiRaa[k]->SetPoint(i,x,val);
	  hJpsiRaa[k]->SetPointError(i,0,err);
	  hJpsiRaa[k]->SetMarkerStyle(21);
	  hJpsiRaa[k]->SetMarkerColor(1);
	  hJpsiRaa[k]->SetLineColor(1);
	  hJpsiRaa[k]->SetMarkerSize(1.5);
	  printf("%s: raa -> %2.2f +/- %2.2f%%, AuAu -> %2.2e +/- %2.2f%%, pp -> %2.2e +/-%2.2f%%\n",pt_Name[i+1],val,err/val*100,AuAu_val,AuAu_err/AuAu_val*100,pp_val,pp_err/pp_val*100);

	  hJpsiRaaSys[k]->SetPoint(i,x,val);
	  hJpsiRaaSys[k]->SetPointError(i,0.1,val*hSys->GetBinContent(i+1));
	  hJpsiRaaSys[k]->SetFillStyle(0);
	}
    }

 
  // RAA plot
  TCanvas *c3 = new TCanvas("Raa_Jpsi_vs_pub_0-10","Raa_Jpsi_vs_pub_0-10",800,600);
  SetPadMargin(gPad,0.13,0.13,0.05,0.02);
  TH1F *hRaa = new TH1F("Raa_Jpsi",";p_{T} (GeV/c);R_{AA}",10,-0.1,9.8);
  hRaa->GetYaxis()->SetRangeUser(0.05,1.6);
  hRaa->GetYaxis()->CenterTitle();
  ScaleHistoTitle(hRaa,28,1,24,28,1,24,63);
  hRaa->SetMaximum(1.7);
  hRaa->DrawCopy();
  TLine *line = GetLine(hRaa->GetXaxis()->GetXmin(),1,hRaa->GetXaxis()->GetXmax(),1,1);
  line->Draw();

  gRaaLowPtSys[1]->Draw("sameE5");
  gRaaHighPtSys[1]->Draw("sameE5");
  gRaaLowPt[1]->Draw("sames PEZ");
  gRaaHighPt[1]->Draw("sames PEZ");

  hJpsiRaaSys[1]->Draw("sameE5");
  hJpsiRaa[1]->Draw("samesPEZ");

  // projection
  TGraphErrors *hJpsiRaaProj = new TGraphErrors();
  for(int i=0; i<hJpsiRaa[1]->GetN(); i++)
    {
      hJpsiRaa[1]->GetPoint(i,x,y);
      hJpsiRaaProj->SetPoint(i,x+0.25,y);
      hJpsiRaaProj->SetPointError(i,0,hJpsiRaa[1]->GetErrorY(i)/sqrt(2./0.3));
    }
  hJpsiRaaProj->SetMarkerStyle(21);
  hJpsiRaaProj->SetMarkerColor(2);
  hJpsiRaaProj->SetLineColor(2);
  hJpsiRaaProj->SetMarkerSize(1.5);
  hJpsiRaaProj->Draw("samesPEZ");

 TGraphErrors *hJpsippError = new TGraphErrors();
  for(int i=0; i<hJpsiRaa[1]->GetN(); i++)
    {
      hJpsiRaa[1]->GetPoint(i,x,y);
      hJpsippError->SetPoint(i,x-0.2,y);
      hJpsippError->SetPointError(i,0,y*hJpsipp->GetBinError(i+1)/hJpsipp->GetBinContent(i+1));
    }
  hJpsippError->SetMarkerStyle(20);
  hJpsippError->SetMarkerColor(4);
  hJpsippError->SetLineColor(4);
  hJpsippError->SetMarkerSize(1.5);
  hJpsippError->Draw("samesPEZ");

  double gerr = sqrt(ppInelasticErr*ppInelasticErr/ppInelastic/ppInelastic+ncollErr[1]*ncollErr[1]/ncoll[1]/ncoll[1]);
  TBox *box = new TBox(0.1,1-gerr,0.2,1+gerr);
  box->SetFillStyle(1001);
  box->SetFillColor(kGray);
  box->Draw();

  TPaveText *t1 = GetPaveText(0.6,0.8,0.2,0.25,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV 0-20%%"));
  t1->Draw();
  TLegend *leg = new TLegend(0.2,0.7,0.46,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(63);
  leg->SetTextSize(23);
  leg->SetHeader("J/#psi#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hJpsippError,"Error for published pp reference: stat+sys","PE");
  leg->AddEntry(hJpsiRaa[1],"QM2015 (0.3*10nb^{-1}, Au+Au stat+sys error)","PE");
  leg->AddEntry(hJpsiRaaProj,"Run14+Run16 (2*10nb^{-1}, Au+Au stat error)","PE");
  leg->AddEntry(gRaaLowPt[1],"Publication J/#psi#rightarrowe^{+}e^{-} (Au+Au stat+sys error)","PE");
  leg->Draw();

  if(savePlot)
    {
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14+16_JpsiRaa_0020.pdf",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14+16_JpsiRaa_0020.png",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14+16_JpsiRaa_0020.jpg",run_type,run_cfg_name.Data()));
      c3->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Run14+16_JpsiRaa_0020.eps",run_type,run_cfg_name.Data()));
    }

}

//================================================
void upsilon(const bool savePlot = 1)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("Rootfiles/UpsFinal.root","read");
  TH1F *hUL = (TH1F*)f->Get("SE_UL");
  TH1F *hLS = (TH1F*)f->Get("SE_LS");
  TH1F *hSub = (TH1F*)f->Get("SE_Ups");
  TF1 *func[4];
  func[0] = (TF1*)f->Get("fun");
  func[1] = (TF1*)f->Get("r1s");
  func[2] = (TF1*)f->Get("r2s");
  func[3] = (TF1*)f->Get("r3s");

  TCanvas *c1 = new TCanvas("Fit_Jpsi_All","Fit_Jpsi_All",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  ScaleHistoTitle(hUL,0.05,1,0.04,0.05,1.1,0.04,62);
  hUL->SetMaximum(120);
  hUL->SetTitle("");
  hUL->SetLineColor(2);
  hUL->Draw();
  hLS->Draw("sames HIST");

  int lowbin = hUL->FindBin(9);
  int highbin = hUL->FindBin(12);
  printf("[i] %d counts in unlike-sign\n",hUL->Integral(lowbin,highbin));
  printf("[i] %d counts in like-sign\n",hLS->Integral(lowbin,highbin));
  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.5,0.55,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hUL,"Unlike-sign pairs","P");
  leg->AddEntry(hLS,"Like-sign pairs","L");
  leg->Draw();
  TPaveText *star = GetPaveText(0.15,0.4,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.pdf",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.png",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.jpg",run_type,run_cfg_name.Data()));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_counts.eps",run_type,run_cfg_name.Data()));
    }

  TCanvas *c2 = new TCanvas("Fit_Upsilon","Fit_Upsilon",800,650);
  SetPadMargin(gPad,0.12,0.12,0.05,0.05);
  //gPad->SetGrid(1,1);
  ScaleHistoTitle(hSub,0.05,1,0.04,0.05,1.1,0.04,62);
  hSub->SetTitle("");
  hSub->SetMarkerColor(2);
  hSub->SetLineColor(1);
  hSub->GetYaxis()->SetRangeUser(-20,40);
  hSub->Draw();
  for(int i=0; i<4; i++)
    {
      func[i]->Draw("sames");
    }
  TF1 *funcBkg = new TF1("funcBkg","[0]*expo(1)",8,12);
  for(int i=0; i<3; i++)
    {
      funcBkg->SetParameter(i,func[0]->GetParameter(9+i));
    }
  funcBkg->SetLineColor(kGreen+1);
  funcBkg->SetLineStyle(2);
  funcBkg->Draw("sames");
  printf("Fit: chi2/ndf = %1.1f/%1.0f\n",func[0]->GetChisquare(), func[0]->GetNDF());

  TPaveText *t1 = GetPaveText(0.15,0.25,0.85,0.9,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("Au+Au @ 200 GeV #it{L} ~ %1.1f nb^{-1}",luminosity));
  t1->Draw();
  TLegend *leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetHeader("#Upsilon#rightarrow#mu^{+}#mu^{-}, |y| < 0.5");
  leg->AddEntry(hSub,"UL-LS","P");
  leg->AddEntry(func[0],"Combined fit","L");
  leg->Draw();
  TPaveText *t1 = GetPaveText(0.65,0.85,0.58,0.63,0.04,62);
  t1->SetTextAlign(11);
  t1->AddText(Form("#chi^{2}/NDF = %1.1f/%d\n",func[0]->GetChisquare(), func[0]->GetNDF()));
  t1->Draw();
  TPaveText *star = GetPaveText(0.7,0.85,0.15,0.2,0.04,20);
  star->AddText("STAR preliminary");
  star->SetTextColor(2);
  star->Draw();

  if(savePlot)
    {
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.pdf",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.png",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.jpg",run_type,run_cfg_name.Data()));
      c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/%s/Upsilon_UL-LS.eps",run_type,run_cfg_name.Data()));
    }
 
}


void scaleGraph(TGraphAsymmErrors *gr, double scale)
{
  int npoints = gr->GetN();
  Double_t x,y;
  for(int i=0; i<npoints; i++)
    {
      gr->GetPoint(i,x,y);
      double err_xlow = gr->GetErrorXlow(i);
      double err_xhigh = gr->GetErrorXhigh(i);
      double err_ylow = gr->GetErrorYlow(i);
      double err_yhigh = gr->GetErrorYhigh(i);
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,err_xlow,err_xhigh,err_ylow*scale,err_yhigh*scale);
    }
}

void drawBoxes(Int_t n, Double_t* x, Double_t dx, Double_t* y, Double_t* dy, Int_t lineWidth, Int_t lineColor){
  for(int i=0;i<n;i++){
    TBox *box = new TBox(x[i]-dx,y[i]-dy[i],x[i]+dx,y[i]+dy[i]);
    box->SetLineWidth(lineWidth);
    box->SetLineColor(lineColor);
    box->SetFillStyle(0);
    box->Draw("lsame");
  }
}

