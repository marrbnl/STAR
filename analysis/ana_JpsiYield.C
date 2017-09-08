TFile *f;
const int year = YEAR;
TString run_cfg_name;

//================================================
void ana_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2013)
    {
      f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
  printf("HFT di-muon events: %4.2f%%\n",hStat->GetBinContent(15)/hStat->GetBinContent(10)*100);

  yieldVsPt();
  //yieldVsNpart();
  //yieldVsLumi();
  //fitYield();
  //pt2scan();
  //HftTracking();
}


//================================================
void fitYield(int icent = 0, int savePlot = 0, int saveHisto = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiYield = (TH1F*)fdata->Get(Form("Jpsi_FitYield_cent%s_weight",cent_Title[icent]));
  TH1F *hInvJsiYield = (TH1F*)hJpsiYield->Clone("hInvJsiYield");
  for(int bin=1; bin<=hJpsiYield->GetNbinsX(); bin++)
    {
      double scale = hJpsiYield->GetBinWidth(bin);
      double pt = hJpsiYield->GetBinCenter(bin);
      hJpsiYield->SetBinContent(bin, hJpsiYield->GetBinContent(bin)/scale);
      hJpsiYield->SetBinError(bin, hJpsiYield->GetBinError(bin)/scale);
      hInvJsiYield->SetBinContent(bin, hInvJsiYield->GetBinContent(bin)/scale/pt);
      hInvJsiYield->SetBinError(bin, hInvJsiYield->GetBinError(bin)/scale/pt);
    }
  TF1 *funcJpsi = new TF1(Form("Fit_%s_tmp",hJpsiYield->GetName()),"exp([0]+[1]*x)",0,15);
  hInvJsiYield->Fit(funcJpsi,"IR0");
  hInvJsiYield->SetTitle(Form("%s: raw invariant J/psi yield;p_{T} (GeV/c);dN/p_{T}dp_{T}",run_type));
  hInvJsiYield->SetMarkerStyle(21);
  c = draw1D(hInvJsiYield);
  gPad->SetLogy();
  funcJpsi->SetLineColor(2);
  funcJpsi->Draw("sames");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiRes/Fit_JpsiFitYield_cent%s.pdf",run_type,cent_Title[icent]));
    }
  
  hJpsiYield->SetTitle(Form("%s: raw J/psi yield;p_{T} (GeV/c);dN/dp_{T}",run_type));
  hJpsiYield->SetMarkerStyle(21);
  c = draw1D(hJpsiYield);
  TF1 *funcInputJpsi = new TF1(Form("Fit_%s",hJpsiYield->GetName()),"exp([0]+[1]*x)*x",0,20);
  funcInputJpsi->SetParameters(funcJpsi->GetParameters());
  gPad->SetLogy();
  funcInputJpsi->SetLineColor(2);
  funcInputJpsi->Draw("sames");

  if(saveHisto)
    {
      fdata->cd();
      funcInputJpsi->Write("",TObject::kOverwrite);
    }
}

//================================================
void yieldVsPt(int savePlot = 0)
{
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TFile *fin = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config),"read");

  bool drawLegend = true;
  if(nCentBins==1) drawLegend = false;

  TList *list = new TList;
  TString legName[nCentBins];
  TH1F *hFitYield[nCentBins];
  TH1F *hBinCountYield[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hFitYield[i] = (TH1F*)fin->Get(Form("Jpsi_FitYield_cent%s%s",cent_Title[i],gWeightName[gApplyWeight]));
      hBinCountYield[i] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[i],gWeightName[gApplyWeight]));
      if(i==2 || i==3)
	{
	  hFitYield[i]->SetBinContent(nPtBins-1,0);
	  hFitYield[i]->SetBinError(nPtBins-1,1e-10);
	  hBinCountYield[i]->SetBinContent(nPtBins-1,0);
	  hBinCountYield[i]->SetBinError(nPtBins-1,1e-10);
	}
      else if(i==4)
	{
	  for(int bin=7; bin<=9; bin++)
	    {
	      hFitYield[i]->SetBinContent(bin,0);
	      hFitYield[i]->SetBinError(bin,1e-10);
	      hBinCountYield[i]->SetBinContent(bin,0);
	      hBinCountYield[i]->SetBinError(bin,1e-10);
	    }
	}
      double count_1, error_1, count_2, error_2;
      count_1 = hFitYield[i]->IntegralAndError(1,nPtBins-1,error_1);
      count_2 = hBinCountYield[i]->IntegralAndError(1,nPtBins-1,error_2);
      printf("Total # of J/psi in %s%%: fit - %1.0f (%2.2f#sigma), bin count - %1.0f (%2.2f#sigma)\n",cent_Name[i],count_1,count_1/error_1,count_2,count_2/error_2);
      scaleHisto(hFitYield[i],1,1,kTRUE,kFALSE,kFALSE);
      list->Add(hFitYield[i]);
      legName[i] = Form("%s%%",cent_Name[i]);
    }
  if(year==2013) c = drawHistos(list,"Jpsi_FitYield","J/psi yield from bin counting;p_{T} (GeV/c);Counts/(1 GeV/c)",kFALSE,0,30,kTRUE,1e-6,1.2*hBinCountYield[0]->GetMaximum(),kFALSE,drawLegend,legName,drawLegend,"Centrality",0.45,0.75,0.55,0.8,kTRUE);
  else c = drawHistos(list,"Jpsi_FitYield","Raw J/psi yield in each p_{T} bin;p_{T} (GeV/c);Counts/(1 GeV/c)",kFALSE,0,30,kTRUE,1,5e4,true,kTRUE,legName,kTRUE,run_type,0.45,0.75,0.65,0.88,kTRUE,0.05);
  
  for(int i=0; i<nCentBins; i++)
    {
      scaleHisto(hBinCountYield[i],1,1,kTRUE,kFALSE,kFALSE);
      hBinCountYield[i]->SetMarkerStyle(24);
      hBinCountYield[i]->SetMarkerColor(color[i]);
      hBinCountYield[i]->SetLineColor(color[i]);
      TGraphErrors *gr = new TGraphErrors(hBinCountYield[i]);
      offset_x(gr,0.2);
      gr->Draw("samesPEZ");
    }
  leg = new TLegend(0.65,0.7,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hFitYield[0],"Fitting","P");
  leg->AddEntry(hBinCountYield[0],"Bin counting","P");
  leg->Draw();
  

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiYieldVsPt.pdf",run_type,run_cfg_name.Data()));
    }

  TH1F *hSignif[nCentBins][2];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hSignif[i][j] = (TH1F*)hFitYield[0]->Clone(Form("hSignif_%s_%d",cent_Name[i],j));
	  hSignif[i][j]->Reset();
	}
      for(int bin=1; bin<=hSignif[i][0]->GetNbinsX(); bin++)
	{
	  hSignif[i][0]->SetBinContent(bin,hFitYield[i]->GetBinContent(bin)/hFitYield[i]->GetBinError(bin));
	  hSignif[i][1]->SetBinContent(bin,hBinCountYield[i]->GetBinContent(bin)/hBinCountYield[i]->GetBinError(bin));
	}
    }

  const char *name[2] = {"Fitting","BinCount"};
  for(int j=0; j<2; j++)
    {
      list->Clear();
      for(int i=0; i<nCentBins; i++)
	{
	  list->Add(hSignif[i][j]);
	}
      c = drawHistos(list,Form("%s_Signif",name[j]),Form("%s: J/psi significance;p_{T} (GeV/c);significance",name[j]),kFALSE,0,20,kTRUE,0.1,20,kFALSE,drawLegend,legName,drawLegend,run_type,0.55,0.7,0.6,0.85,kTRUE);
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiSignifVsPt_%s.pdf",run_type,run_cfg_name.Data(),name[j]));
	}
    }

  // compare fit parameters
  TH1F *hMean[nCentBins];
  TH1F *hSigma[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hMean[i] = (TH1F*)fin->Get(Form("Jpsi_FitMean_cent%s_weight",cent_Title[i]));
      hSigma[i] = (TH1F*)fin->Get(Form("Jpsi_FitSigma_cent%s_weight",cent_Title[i]));
    }
  list->Clear();
  for(int i=0; i<nCentBins; i++)
    {
      if(i==2 || i==3)
	{
	  hMean[i]->SetBinContent(nPtBins-1,0);
	  hMean[i]->SetBinError(nPtBins-1,1e-10);
	}
      else if(i==4)
	{
	  for(int bin=7; bin<=9; bin++)
	    {
	      hMean[i]->SetBinContent(bin,0);
	      hMean[i]->SetBinError(bin,1e-10);
	    }
	}
      list->Add(hMean[i]);
    }
  c = drawHistos(list,"JpsiMeanVsPt","Mean of J/psi peak;p_{T} (GeV/c);Mean",kFALSE,0,20,kTRUE,3.04,3.2,kFALSE,drawLegend,legName,drawLegend,run_type,0.15,0.35,0.6,0.88,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiMeanVsPt.pdf",run_type,run_cfg_name.Data()));
    }

  list->Clear();
  for(int i=0; i<nCentBins; i++)
    {
      if(i==2 || i==3)
	{
	  hSigma[i]->SetBinContent(nPtBins-1,-1);
	  hSigma[i]->SetBinError(nPtBins-1,1e-10);
	}
      else if(i==4)
	{
	  for(int bin=7; bin<=9; bin++)
	    {
	      hSigma[i]->SetBinContent(bin,-1);
	      hSigma[i]->SetBinError(bin,1e-10);
	    }
	}
      list->Add(hSigma[i]);
    }
  c = drawHistos(list,"JpsiSigmaVsPt","Width of J/psi peak;p_{T} (GeV/c);#sigma",kFALSE,0,20,kTRUE,0,0.15,kFALSE,drawLegend,legName,drawLegend,run_type,0.15,0.35,0.6,0.88,kTRUE);
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiSigmaVsPt.pdf",run_type,run_cfg_name.Data()));
    }
}

//================================================
void yieldVsNpart(int savePlot = 1)
{
  const int nPtBins         = nPtBins_npart;
  const double* ptBins_low  = ptBins_low_npart;
  const double* ptBins_high = ptBins_high_npart;
  const char** pt_Name      = pt_Name_npart;

  TFile *fin = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config),"read");
  TH1F *hFitYield[nPtBins];
  TH1F *hBinCountYield[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hFitYield[i] = (TH1F*)fin->Get(Form("Jpsi_FitYield_pt%s%s",pt_Name[i],gWeightName[gApplyWeight]));
      hBinCountYield[i] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_pt%s%s",pt_Name[i],gWeightName[gApplyWeight]));
      
      hFitYield[i]->SetMarkerStyle(20);
      hFitYield[i]->SetMarkerColor(i+1);
      hFitYield[i]->SetLineColor(i+1);
      hFitYield[i]->SetMarkerSize(1.2);
      hFitYield[i]->GetXaxis()->SetLabelSize(0.045);
      hFitYield[i]->GetYaxis()->SetLabelSize(0.035);
      c = draw1D(hFitYield[i],Form("%s: raw J/psi yield in each centrality bin",run_type));
      gPad->SetLogy();
      hBinCountYield[i]->SetMarkerStyle(24);
      hBinCountYield[i]->SetMarkerColor(i+1);
      hBinCountYield[i]->SetLineColor(i+1);
      TGraphErrors *gr = new TGraphErrors(hBinCountYield[i]);
      offset_x(gr,0.2);
      gr->Draw("samesPEZ");
      leg = new TLegend(0.2,0.25,0.4,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetHeader(Form("p_{T} > %1.0f GeV/c",ptBins_low[i]));
      leg->AddEntry(hFitYield[i],"Fitting","P");
      leg->AddEntry(hBinCountYield[i],"Bin counting","P");
      leg->Draw();
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiYieldVsCent_pt%s.pdf",run_type,run_cfg_name.Data(),pt_Name[i]));
	}
    }

  TH1F *hSignif[nPtBins][2];
  for(int i=0; i<nPtBins; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hSignif[i][j] = (TH1F*)hFitYield[i]->Clone(Form("hSignif_%s_%d",pt_Name[i],j));
	  hSignif[i][j]->Reset();
	}
      for(int bin=1; bin<=hSignif[i][0]->GetNbinsX(); bin++)
	{
	  hSignif[i][0]->SetBinContent(bin,hFitYield[i]->GetBinContent(bin)/hFitYield[i]->GetBinError(bin));
	  hSignif[i][1]->SetBinContent(bin,hBinCountYield[i]->GetBinContent(bin)/hBinCountYield[i]->GetBinError(bin));
	}
      hSignif[i][0]->GetYaxis()->SetRangeUser(0,20);
      hSignif[i][0]->SetMarkerColor(i+1);
      hSignif[i][0]->SetMarkerSize(1.2);
      c = draw1D(hSignif[i][0],Form("%s: significance of J/psi signal from fitting (p_{T} > %1.0f GeV/c)",run_type,ptBins_low[i]));
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiSignifVsCent_pt%s_Fitting.pdf",run_type,run_cfg_name.Data(),pt_Name[i]));
	}
    }
}


//================================================
void HftTracking(int savePlot = 1)
{
  const int ipt = 0, icent = 0;
  // Events with HFT
  const int hftindex = 1;
  const char *hftname[2] = {"Hft",""};
  const char *hfttitle[2] = {" with HFT",""};

  TFile *fin = TFile::Open("output/Pico.Run14.AuAu200.jpsi.Vz6cm.root","read");
  TH1F *hStat = (TH1F*)fin->Get("hEventStat");
  printf("+++ check this +++\n");
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
  printf("HFT di-muon events: %4.4e, %4.2f%%\n",hStat->GetBinContent(15),hStat->GetBinContent(15)/hStat->GetBinContent(10)*100);
  const double events[2] = {hStat->GetBinContent(15),hStat->GetBinContent(10)};

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3] = {0x0};
  TH1F *hInvMass[3] = {0x0};
  
  // same event
  char name[512];
  for(Int_t j=0; j<3; j++) // pair type
    { 
      hnInvMass[j] = (THnSparseF*)fin->Get(Form("m%s%s_%s",hName[j],hftname[hftindex],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRange(centBins_low[icent],centBins_high[icent]);
      hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins_low[ipt]+0.01,ptBins_high[ipt]-0.01);

      hInvMass[j] = (TH1F*)hnInvMass[j]->Projection(0);
      hInvMass[j]->SetName(Form("InvMassHft_jpsi_%d",j));
      hInvMass[j]->Sumw2();
    }
  hInvMass[1]->Add(hInvMass[2]);

  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.histo.root",pt1_cut,pt2_cut));
  TFile *fin = TFile::Open(fileName,"read");

  // mix event
  TH1F *hSeUL  = (TH1F*)hInvMass[0]->Clone("hSeUL");
  TH1F *hSeLS  = (TH1F*)hInvMass[1]->Clone("hSeLS");
  TH1F *hMixUL = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[ipt],cent_Name[icent]));
  TH1F *hMixLS = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[ipt],cent_Name[icent]));

  double g_mix_scale_low = 2.5;
  double g_mix_scale_high = 4;
  double g_bin_width = 0.04; // 40 MeV
  TString g_func1 = "pol3";
  TString g_func2 = "pol0";
  int g_func1_npar = 4;
  int g_func2_npar = 1;
  double g_sig_fit_min = 2.6;
  double g_sig_fit_max = 4.0;

  double se = 0, se_err = 0, me = 0, me_err = 0;
  int low_bin = hSeLS->FindFixBin(g_mix_scale_low+1e-4);
  int high_bin = hSeLS->FindFixBin(g_mix_scale_high-1e-4);
  se = hSeLS->IntegralAndError(low_bin,high_bin,se_err);
  
  int low_bin_me = hMixLS->FindFixBin(g_mix_scale_low+1e-4);
  int high_bin_me = hMixLS->FindFixBin(g_mix_scale_high-1e-4);
  me = hMixLS->IntegralAndError(low_bin,high_bin,me_err);

  double scale = se/me;
  double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);

  // jpsi signal
  TH1F *hMixBkg = (TH1F*)hMixUL->Clone(Form("mix_bkg_pt%s_cent%s",pt_Name[ipt],cent_Name[icent]));
  hMixBkg->Scale(scale);
  hSeUL->Rebin(g_bin_width/hSeUL->GetBinWidth(1));
  hSeUL->SetTitle("");
  hSeUL->SetMarkerStyle(21);
  hSeUL->SetMarkerColor(2);
  hSeUL->SetLineColor(2);
  hSeUL->GetXaxis()->SetRangeUser(2.5,4);
  hSeLS->Rebin(g_bin_width/hSeLS->GetBinWidth(1));
  hMixBkg->SetLineColor(4);
  hMixBkg->Rebin(g_bin_width/hMixBkg->GetBinWidth(1));
  hSeUL->SetMaximum(1.2*hSeUL->GetMaximum());
  c = draw1D(hSeUL,Form("Dimuon events%s: %1.1f < p_{T} < %1.1f GeV/c (%s%%)",hfttitle[hftindex],ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]));
  hSeLS->Draw("sames HIST");
  hMixBkg->Draw("sames HIST");
  leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hSeUL,"Unlike sign","P");
  leg->AddEntry(hSeLS,"Like sign (++)+(--)","L");
  leg->AddEntry(hMixBkg,"Mix UL","L");
  leg->Draw();

  // Fit residual
  TH1F *hSignal = (TH1F*)hSeUL->Clone(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[ipt]));
  //hSignal->Add(hMixBkg,-1);
  hSignal->Add(hSeLS,-1);
  hSignal->SetLineColor(1);
  hSignal->SetMarkerColor(1);
  TF1 *funcSignal = new TF1(Form("Jpsi_FitSig_pt%s",pt_Name[ipt]),Form("gausn(0)+pol3(3)"),g_sig_fit_min,g_sig_fit_max);
  funcSignal->SetParameter(0,100);
  funcSignal->SetParameter(1,3.09);
  funcSignal->SetParameter(2,0.1);
  funcSignal->SetLineColor(2);
  TFitResultPtr ptr = hSignal->Fit(funcSignal,"IRS0Q");
  double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
  double fit_yield = funcSignal->GetParameter(0)/hSignal->GetBinWidth(1);
  double fit_yield_err = funcSignal->GetParError(0)/hSignal->GetBinWidth(1);

  // Extract background matrix
  const int nParameter = 4;
  double bkg_params[nParameter];
  double bkg_matrix[nParameter*nParameter];
  TF1 *funcBkg = new TF1(Form("Jpsi_FitBkg_pt%s_%s",pt_Name[ipt]),"pol3",g_sig_fit_min,g_sig_fit_max);
  for(int j=0; j<nParameter; j++)
    {
      funcBkg->SetParameter(j,funcSignal->GetParameter(3+j));
      funcBkg->SetParError(j,funcSignal->GetParError(3+j));
      bkg_params[j] = funcSignal->GetParameter(3+j);
    }
 
  for(int j=3; j<3+nParameter; j++)
    {
      for(int k=3; k<3+nParameter; k++)
	{
	  bkg_matrix[(j-3)*nParameter+k-3] = matrix[j*(3+nParameter)+k];
	}
    }

  // bin counting
  double low_mass_tmp = 2.96;;
  double high_mass_tmp = 3.24;
  int low_bin = hSignal->FindFixBin(low_mass_tmp+1e-4);
  int high_bin = hSignal->FindFixBin(high_mass_tmp-1e-4);
  double yield_all_err;
  double yield_all = hSignal->IntegralAndError(low_bin,high_bin,yield_all_err);
  double yield_bkg = funcBkg->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal->GetBinWidth(1);
  double yield_bkg_err = funcBkg->IntegralError(low_mass_tmp,high_mass_tmp,bkg_params,bkg_matrix)* 1./hSignal->GetBinWidth(1);
  double yield_sig = yield_all - yield_bkg;
  double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err);
  
  hSignal->SetMaximum(2.5*hSignal->GetMaximum());
  c = draw1D(hSignal,Form("Dimuon events%s: %1.1f < p_{T} < %1.1f GeV/c (%s%%)",hfttitle[hftindex],ptBins_low[ipt],ptBins_high[ipt],cent_Name[icent]));
  funcSignal->Draw("sames");
  funcBkg->SetLineColor(4);
  funcBkg->Draw("sames");
  TLine *line = GetLine(low_mass_tmp,hSignal->GetMinimum()*1.5,low_mass_tmp,hSignal->GetMaximum()*0.3,1);
  line->Draw();
  TLine *line = GetLine(high_mass_tmp,hSignal->GetMinimum()*1.5,high_mass_tmp,hSignal->GetMaximum()*0.3,1);
  line->Draw();

  t = GetPaveText(0.16,0.3,0.63,0.88,0.04);
  t->SetTextFont(62);
  t->AddText(Form("N_{evt} = %2.1e",events[hftindex]));
  t->AddText("Fitting");
  t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",fit_yield,fit_yield_err,fit_yield/fit_yield_err));
  t->AddText("Counting");
  t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",yield_sig,yield_sig_err,yield_sig/yield_sig_err));
  t->SetTextAlign(11);
  t->Draw();

  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sDimuon%s_JpsiInvMass.pdf",run_type,run_cfg_name.Data(),hftname[hftindex]));

}

//================================================
void yieldVsLumi(int savePlot = 1, int saveHisto = 1)
{
  const int gApplyWeight = 1;
  const int nSetup = 4;
  const char *setup_name[nSetup] = {"production_2014","production_low_2014","production_mid_2014","production_high_2014"};
  TList *list = new TList;
  TString legName[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      legName[i] = setup_name[i];
    }

  // luminosity
  TH2F *hTrgSetupVsZdcRate = (TH2F*)f->Get("mhTrgSetupVsZdcRate_di_mu");
  TH1F *hZdcRate[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      hZdcRate[i] = (TH1F*)hTrgSetupVsZdcRate->ProjectionX(Form("hZdcRate_%d",i),i+1,i+1);
    }
  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hZdcRate[i]);
    }
  c = drawHistos(list,"ZdcRate","Distribution of ZDC rate",true,0,120,true,1,5e7,false,true,legName,true,run_type,0.5,0.8,0.62,0.88,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sZdcRate_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  fileName.ReplaceAll(".root",".yield.root");
  TFile *fin;
  if(saveHisto) fin = TFile::Open(fileName,"update");
  else          fin = TFile::Open(fileName,"read");

  TH1F *hJpsiYield[nSetup][nCentBins];
  for(int i=0; i<nSetup; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hJpsiYield[i][j] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_cent%s%s%s",cent_Title[j],gWeightName[gApplyWeight],gTrgSetupName[i+1]));
	}
    }

  // # of events in each setup/centrality
  TH1F *hCent = (TH1F*)f->Get("mhCentrality_di_mu");
  double centFraction[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      centFraction[i] = hCent->Integral(centBins_low[i],centBins_high[i])/hCent->Integral(0,-1);
    }


  // # of Jpsi signal per event
  double nEvents[nSetup][nCentBins];
  TH1F *hEvtRunAcc = (TH1F*)f->Get("mhEvtRunAcc_di_mu");
  for(int i=0; i<nSetup; i++)
    {
      for(int j=0; j<nCentBins; j++) nEvents[i][j] = 0;
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s.list",run_type,setup_name[i]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  if(runnumber>0)
	    {
	      for(int j=0; j<nCentBins; j++)
		{
		  nEvents[i][j] += hEvtRunAcc->GetBinContent(hEvtRunAcc->FindFixBin(runnumber)) * centFraction[j];
		}
	    }
	}
    }

  // # of Jpsi in each setup/centrality
  TH1F *hJpsiInCent[nSetup];
  TH1F *hJpsiPerEvent[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      hJpsiPerEvent[i] = new TH1F(Form("hJpsiPerEvent%s%s",gWeightName[gApplyWeight],gTrgSetupName[i+1]),"Jpsi counts per event in each centrality bin (p_{T} > 0)",nCentBins-1,0,nCentBins-1);
      hJpsiInCent[i] = new TH1F(Form("NJpsiInCent%s%s",gWeightName[gApplyWeight],gTrgSetupName[i+1]),"Jpsi counts in each centrality bin (p_{T} > 0)",nCentBins-1,0,nCentBins-1);

      for(int j=0; j<nCentBins; j++)
	{
	  printf("# of events in %s %s: %1.0f\n",gTrgSetupName[i+1],cent_Title[j],nEvents[i][j]);
	  hJpsiInCent[i]->SetBinContent(j, hJpsiYield[i][j]->GetBinContent(0));
	  hJpsiInCent[i]->SetBinError(j, hJpsiYield[i][j]->GetBinError(0));
	  if(j>0) hJpsiInCent[i]->GetXaxis()->SetBinLabel(j,Form("%s%%",cent_Name[j]));

	  hJpsiPerEvent[i]->SetBinContent(j, hJpsiYield[i][j]->GetBinContent(0)/nEvents[i][j]);
	  hJpsiPerEvent[i]->SetBinError(j, hJpsiYield[i][j]->GetBinError(0)/nEvents[i][j]);
	  if(j>0) hJpsiPerEvent[i]->GetXaxis()->SetBinLabel(j,Form("%s%%",cent_Name[j]));
	}
      hJpsiInCent[i]->GetXaxis()->SetLabelSize(0.05);
      hJpsiPerEvent[i]->GetXaxis()->SetLabelSize(0.05);
    }

  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hJpsiInCent[i]);
    }
  c = drawHistos(list,"Jpsi_Counts","",kFALSE,0,30,true,0,1e4,kFALSE,true,legName,true,run_type,0.45,0.75,0.62,0.88,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiCounts_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hJpsiPerEvent[i]);
    }
  c = drawHistos(list,"Jpsi_Yield","",kFALSE,0,30,true,1e-5,5e-4,true,true,legName,true,run_type,0.15,0.3,0.62,0.88,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiYield_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nSetup; i++)
	{
	  hJpsiInCent[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void pt2scan(int savePlot = 0)
{
  const int ncuts = 2;
  double pt2_cuts[ncuts] = {1.0,1.2};
  TFile *fin[ncuts];
  TH1F *hJpsiYield[ncuts];
}
