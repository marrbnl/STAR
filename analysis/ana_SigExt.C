TFile *f = 0x0;
const int year = YEAR;
TString run_cfg_name;

//================================================
void ana_SigExt()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  f = TFile::Open(Form("./output/%s.jpsi.%sroot",run_type,run_config),"read");
  run_cfg_name = Form("%s",run_config);

  if(f)
    {
      TH1F *hStat = (TH1F*)f->Get("hEventStat");
      printf("all         events: %4.4e\n",hStat->GetBinContent(1));
      printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
      printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
    }

  
  //prod_pt();
  Run14_signal();
  //Run14_npart();
}

//================================================
void prod_pt()
{
  for(int t=0; t<5; t++)
    {
      for(int i=0; i<4; i++)
	{
	  for(int j=0; j<8; j++)
	    {
	      if(t>0 && j>0) continue;
	      Run14_signal(t,i,j,1,1);
	    }
	}
    }
}

//===============================================
void Run14_signal(const int isetup = 0, const int icent = 0, const int isys = 0, int savePlot = 0, int saveHisto = 0)
{
  // re-assign global constants
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  // prepare name and title
  const int nSys = 10;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc","_LargeFit","_SmallFit","_SigFunc","_FixSig"};
  const char *sys_title[nSys] = {"","Sys.LargeScale.","Sys.SmallScale.","Sys.ScaleFit.","Sys.Rebin.","Sys.BkgFunc.","Sys.LargeFit.","Sys.SmallFit.","Sys.SigFunc.","Sys.FixSig."};
  const TString suffix = Form("cent%s%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup],sys_name[isys]);
  const TString suf_title = Form("cent%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]);

  // open input and output files
  TString fileName   = Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);
  TFile *fin = TFile::Open(fileName,"read");

  printf("\n===== Running configuration =====\n");
  printf("Input filename: %s\n",fileName.Data());
  printf("Histogram name: %s\n",suffix.Data());
  printf("Centrality: %s\n",cent_Name[icent]);
  printf("==================================\n\n");

  // global setup
  double  g_mix_scale_low = 2.7;
  double  g_mix_scale_high = 3.5;
  TString g_mix_func = "pol0";
  double  g_bin_width = 0.05;
  TString g_func1 = "pol3";
  int     g_func1_npar = 4;
  TString g_func2 = "pol1";
  int     g_func2_npar = 2;
  double  g_sig_fit_min = 2.5;
  double  g_sig_fit_max = 4.0;
  if(isys==1) { g_mix_scale_low = 2.5; g_mix_scale_high = 3.7; }
  if(isys==2) { g_mix_scale_low = 2.9; g_mix_scale_high = 3.3; }
  if(isys==3) { g_mix_func = "pol1"; g_mix_scale_low = 2.5; g_mix_scale_high = 4.0; g_func1 = "pol3"; g_func1_npar = 4; }
  if(isys==4) { g_bin_width = 0.02; }
  if(isys==5) { g_func1 = "pol4"; g_func1_npar = 5; g_func2 = "pol0"; g_func2_npar = 1; }
  if(isys==6) { g_sig_fit_min = 2.3; g_sig_fit_max = 4.2; }
  if(isys==7) { g_sig_fit_min = 2.6; g_sig_fit_max = 3.8; }
  if(isys==8) { g_sig_fit_min = 2.5; g_sig_fit_max = 4.0; g_bin_width = 0.04;}


  // get histograms
  TH1F *hSeUL[nPtBins];
  TH1F *hSeLS[nPtBins];
  TH1F *hMixUL[nPtBins];
  TH1F *hMixLS[nPtBins];
  TH1F *hMixBkg[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      hSeUL[i] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s%s%s",pt_Name[i],cent_Name[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]));
      hSeLS[i] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s%s%s",pt_Name[i],cent_Name[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]));
      hMixUL[i] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixLS[i] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
    }
  const int binWidthScale = int(hSeUL[0]->GetBinWidth(1)/hMixUL[0]->GetBinWidth(1));
  
  int nxpad, nypad;
  if(nPtBins-1<=4) { nxpad = 2; nypad = 2; }
  else             { nxpad = (nPtBins-1)/2 + (nPtBins-1)%2; nypad = 2; }

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // mixed event scaling
  TH1F *hMixScale[nPtBins];
  TF1 *funcScale[nPtBins];
  TH1F *hFitScaleFactor = new TH1F(Form("FitScaleFactor_%s",suffix.Data()),"Mixed-event scale factor from fitting",nbins,xbins);
  TH1F *hBinCountScaleFactor = new TH1F(Form("BinCountScaleFactor_%s",suffix.Data()),"Mixed-event scale factor from bin counting",nbins,xbins);
  TCanvas *cScaling = new TCanvas(Form("mix_scale_%s",cent_Name[icent]),Form("mix_scale_%s",cent_Name[icent]),1100,650);
  cScaling->Divide(nxpad,nypad);
  for(int i=0; i<nPtBins; i++)
    {
      hMixScale[i] = (TH1F*)hSeLS[i]->Clone(Form("LS_ratio_SE_to_ME_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixScale[i]->Rebin(5);
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone",hMixLS[i]->GetName()));
      htmp->Rebin(int(hMixScale[i]->GetBinWidth(1)/hMixLS[i]->GetBinWidth(1)));
      hMixScale[i]->Divide(htmp);

      double g_mix_scale_low_tmp = g_mix_scale_low;
      double g_mix_scale_high_tmp = g_mix_scale_high;

      // fitting method
      funcScale[i] = new TF1(Form("Fit_%s",hMixScale[i]->GetName()),g_mix_func,g_mix_scale_low_tmp,g_mix_scale_high_tmp);
      hMixScale[i]->Fit(funcScale[i],"IR0Q");
      hFitScaleFactor->SetBinContent(i,funcScale[i]->GetParameter(0));
      hFitScaleFactor->SetBinError(i,funcScale[i]->GetParError(0));

      // bin counting method
      double se = 0, se_err = 0, me = 0, me_err = 0;
      int low_bin = hSeLS[i]->FindFixBin(g_mix_scale_low_tmp+1e-4);
      int high_bin = hSeLS[i]->FindFixBin(g_mix_scale_high_tmp-1e-4);
      for(int bin=low_bin; bin<=high_bin; bin++)
	{
	  se += hSeLS[i]->GetBinContent(bin);
	  se_err += TMath::Power(hSeLS[i]->GetBinError(bin),2);
	}

      int low_bin_me = hMixLS[i]->FindFixBin(g_mix_scale_low_tmp+1e-4);
      int high_bin_me = hMixLS[i]->FindFixBin(g_mix_scale_high_tmp-1e-4);
      for(int bin=low_bin_me; bin<=high_bin_me; bin++)
	{
	  me += hMixLS[i]->GetBinContent(bin);
	  me_err += TMath::Power(hMixLS[i]->GetBinError(bin),2);
	}
      se_err = TMath::Sqrt(se_err);
      me_err = TMath::Sqrt(me_err);
      double scale = se/me;
      double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);
      hBinCountScaleFactor->SetBinContent(i,scale);
      hBinCountScaleFactor->SetBinError(i,scale_error);

      if(i>-1)
	{
	  // plotting
	  cScaling->cd(i+1);
	  hMixScale[i]->SetTitle("");
	  hMixScale[i]->SetMarkerStyle(21);
	  hMixScale[i]->GetXaxis()->SetRangeUser(2.5,4);
	  hMixScale[i]->SetMaximum(1.1*hMixScale[i]->GetMaximum());
	  if(i<5) hMixScale[i]->GetYaxis()->SetRangeUser(0.6e-3, 0.75e-3);
	  hMixScale[i]->Draw();
	  funcScale[i]->SetLineColor(2);
	  funcScale[i]->Draw("sames");
	  TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
	  t->Draw();
	  t = GetPaveText(0.16,0.3,0.7,0.85,0.05);
	  t->SetTextFont(62);
	  t->AddText(Form("Fit: %4.3e",funcScale[i]->GetParameter(0)));
	  t->AddText(Form("Count: %4.3e",scale));
	  t->SetTextAlign(11);
	  t->Draw();
	}
    }
  t = GetPaveText(0.15,0.35,0.5,0.6,0.06);
  t->SetTextFont(62);
  t->AddText("LS: SE/ME");
  t->SetTextColor(4);
  cScaling->cd(1);
  t->Draw();
  if(savePlot && isys==0) 
    {
      cScaling->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sBkg.FitSEtoME_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  TList *list = new TList;
  TString legName[2] = {"Fitting","Bin Counting"};
  list->Add(hFitScaleFactor);
  list->Add(hBinCountScaleFactor);
  c = drawHistos(list,"Mix_scaleFactor",Form("Scale factor for mixed events (%s%%)",cent_Name[icent]),kFALSE,0,30,kTRUE,hFitScaleFactor->GetMinimum()*0.95,hFitScaleFactor->GetMaximum()*1.08,kFALSE,kTRUE,legName,kTRUE,"Scale factor",0.5,0.7,0.68,0.85,kTRUE);
  if(savePlot && isys==0) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sBkg.ScaleFactor_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  // Acceptance
  TCanvas *cAcc = new TCanvas(Form("Acceptance_%s",cent_Name[icent]),Form("Acceptance_%s",cent_Name[icent]),1100,650);
  cAcc->Divide(nxpad,nypad);
  TH1F *hAcc[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hAcc[i] = (TH1F*)hMixUL[i]->Clone(Form("Mix_Acceptance_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone2",hMixLS[i]->GetName()));
      hAcc[i]->Rebin(g_bin_width/hAcc[i]->GetBinWidth(1));
      htmp->Rebin(g_bin_width/htmp->GetBinWidth(1));
      hAcc[i]->Divide(htmp);
      if(i>0)
	{
	  cAcc->cd(i);
	  hAcc[i]->SetTitle("");
	  hAcc[i]->SetMarkerStyle(21);
	  hAcc[i]->GetXaxis()->SetRangeUser(2,4);
	  hAcc[i]->GetYaxis()->SetRangeUser(0.9,1.1);
	  hAcc[i]->Draw("PE");
	  TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
	  t->Draw();
	  TLine *line = GetLine(2,1,4,1,1);
	  line->Draw();
	}
    }
  t = GetPaveText(0.3,0.4,0.2,0.3,0.07);
  t->SetTextFont(62);
  t->AddText("Mix: UL/LS");
  t->SetTextColor(4);
  cAcc->cd(1);
  t->Draw();
  if(savePlot && isys==0) 
    {
      cAcc->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sBkg.Mix_Acceptance_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  // jpsi signal
  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",cent_Name[icent]),Form("InvMass_%s",cent_Name[icent]),1100,650);
  cSignal->Divide(nxpad,nypad);

  for(int i=0; i<nPtBins; i++)
    {
      hMixBkg[i] = (TH1F*)hMixUL[i]->Clone(Form("mix_bkg_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      
      // scale mix event background
      int low_bin = hMixBkg[i]->FindFixBin(2.3);
      int up_bin  = hMixBkg[i]->FindFixBin(4.2);
      for(int ibin=low_bin; ibin<=up_bin; ibin++)
	{
	  double mass = hMixBkg[i]->GetBinCenter(ibin);
	  double scale = funcScale[i]->Eval(mass);
	  hMixBkg[i]->SetBinContent(ibin, scale * hMixBkg[i]->GetBinContent(ibin));
	  hMixBkg[i]->SetBinError(ibin, scale * hMixBkg[i]->GetBinError(ibin));
	}
      
      hSeUL[i]->Rebin(g_bin_width/hSeUL[i]->GetBinWidth(1));
      hSeUL[i]->SetTitle("");
      hSeUL[i]->SetMarkerStyle(21);
      hSeUL[i]->SetMarkerColor(2);
      hSeUL[i]->SetLineColor(2);
      hSeUL[i]->GetXaxis()->SetRangeUser(2.5,4);
      hSeLS[i]->Rebin(g_bin_width/hSeLS[i]->GetBinWidth(1));
      hMixBkg[i]->SetLineColor(4);
      hMixBkg[i]->Rebin(g_bin_width/hMixBkg[i]->GetBinWidth(1));
      if(i>0)
	{
	  cSignal->cd(i);
	  hSeUL[i]->Draw();
	  hSeLS[i]->Draw("sames HIST");
	  hMixBkg[i]->Draw("sames HIST");
	  TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
	  t->Draw();
	}
    }
  cSignal->cd(10);
  leg = new TLegend(0.2,0.6,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.07);
  leg->AddEntry(hSeUL[1],"Unlike sign","P");
  leg->AddEntry(hSeLS[1],"Like sign (++)+(--)","L");
  leg->AddEntry(hMixBkg[1],"Mix UL","L");
  leg->Draw();
  if(savePlot && isys==0) 
    {
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sSig.ULvsLS_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  // Fit residual
  double fix_mean[nPtBins] = {0};
  double fix_sigma[nPtBins];
  if( (icent>0 || isetup > 0) && isys<9)
    {
      TFile *fFit = 0;
      if(isys==0)  fFit = TFile::Open(outName,"read");
      else         fFit = TFile::Open(outNameSys,"read");
      TH1F *h1 = (TH1F*)fFit->Get(Form("Jpsi_FitMean_cent%s%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0],sys_name[isys]));
      TH1F *h2 = (TH1F*)fFit->Get(Form("Jpsi_FitSigma_cent%s%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0],sys_name[isys]));
      for(int bin=0; bin<=h1->GetNbinsX(); bin++)
	{
	  fix_mean[bin] = h1->GetBinContent(bin);
	  fix_sigma[bin] = h2->GetBinContent(bin);
	}
      fFit->Close();
    }
  if(isys==9)
    {
      TFile *fFit = TFile::Open(outName,"read");
      TF1 *func1 = (TF1*)fFit->Get(Form("Jpsi_FitMeanFit_cent%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0]));
      TF1 *func2 = (TF1*)fFit->Get(Form("Jpsi_FitSigmaFit_cent%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0]));
      for(int bin=0; bin<nPtBins; bin++)
	{
	  double pt = (ptBins_low[bin] + ptBins_high[bin]) * 0.5;
	  cout << pt << endl;
	  fix_mean[bin] = func1->Eval(pt);
	  fix_sigma[bin] = func2->Eval(pt);
	}
      fFit->Close();
    }

  
  TCanvas *cFit1 = new TCanvas(Form("Fit_Jpsi_%s_All",cent_Name[icent]),Form("Fit_Jpsi_%s_All",cent_Name[icent]),800,650);
  TCanvas *cFit = new TCanvas(Form("Fit_Jpsi_%s",cent_Name[icent]),Form("Fit_Jpsi_%s",cent_Name[icent]),1400,650);
  cFit->Divide(nxpad,nypad);
  TH1F *hSignal[nPtBins];
  TH1F *hSignalSave[nPtBins];
  TF1 *funcSignal[nPtBins];
  TF1 *funcBkg[nPtBins];
  TH1F *hMean = new TH1F(Form("Jpsi_FitMean_%s",suffix.Data()),"Mean of Jpsi peak",nbins,xbins);
  TH1F *hSigma = new TH1F(Form("Jpsi_FitSigma_%s",suffix.Data()),"Sigma of Jpsi peak",nbins,xbins);
  TH1F *hFitYield = new TH1F(Form("Jpsi_FitYield_%s",suffix.Data()),"Jpsi yield from fitting",nbins,xbins);
  TH1F *hBinCountYield = new TH1F(Form("Jpsi_BinCountYield_%s",suffix.Data()),"Jpsi yield from bin counting",nbins,xbins);
  TH1F *hSigToBkg = new TH1F(Form("Jpsi_SigToBkg_%s",suffix.Data()),"Signal-to-background ratio for Jpsi",nbins,xbins);

  TString funcForm;
  int nPar;
  for(int i=0; i<nPtBins; i++)
    {
      printf("+++++ %1.0f < pT < %1.0f +++++\n",ptBins_low[i],ptBins_high[i]);

      hSignal[i] = (TH1F*)hSeUL[i]->Clone(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      for(int bin=1; bin<=hSignal[i]->GetNbinsX(); bin++)
	{
	  if(hSignal[i]->GetBinContent(bin)==0)
	    {
	      hSignal[i]->SetBinContent(bin,0);
	      hSignal[i]->SetBinError(bin,1.4);
	    }
	}
      if(i==-1)
	{
	  hSeLS[i]->Multiply(hAcc[i]);
	  hSignal[i]->Add(hSeLS[i],-1);
	}
      else
	{
	  hSignal[i]->Add(hMixBkg[i],-1);
	}
      hSignal[i]->SetLineColor(1);
      hSignal[i]->SetMarkerColor(1);
      hSignalSave[i] = (TH1F*)hSignal[i]->Clone(Form("Jpsi_Signal_pt%s_%s",pt_Name[i],suffix.Data()));

      // Fit signal
      if(i<5)
      	{
	  funcForm = g_func1; 
	  nPar = g_func1_npar;
      	}
      else
      	{
	  funcForm = g_func2; 
	  nPar = g_func2_npar;
      	}

      if(isys!=8)
	{
	  funcSignal[i] = new TF1(Form("Jpsi_FitSig_pt%s_%s",pt_Name[i],suffix.Data()),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	  funcSignal[i]->SetParameter(0,hSignal[i]->GetMaximum());
	  funcSignal[i]->SetParameter(1,3.09);
	  funcSignal[i]->SetParameter(2,0.05);
	}
      else
	{
	  // use crystal-ball function to fit
	  if(funcForm=="pol1") funcSignal[i] = new TF1(Form("Jpsi_FitSig_pt%s_%s",pt_Name[i],suffix.Data()),CrystalBallPlusPol1,g_sig_fit_min,g_sig_fit_max,7);
	  if(funcForm=="pol3") funcSignal[i] = new TF1(Form("Jpsi_FitSig_pt%s_%s",pt_Name[i],suffix.Data()),CrystalBallPlusPol3,g_sig_fit_min,g_sig_fit_max,9);
	  if(i>7) funcSignal[i]->SetRange(2.7,3.9);
	  funcSignal[i]->SetParameter(0,hSignal[i]->GetMaximum());
	  funcSignal[i]->SetParameter(1,3.09);
	  funcSignal[i]->SetParameter(2,0.1);
	  funcSignal[i]->SetParameter(3,0.1);
	  funcSignal[i]->SetParameter(4,0.2);
	}
      if(i==1)
	{
	  funcSignal[i]->SetParameter(2,0.1);
	  if(isys==6) funcSignal[i]->SetRange(2.7,g_sig_fit_max);
	  else if(isys==7) funcSignal[i]->SetRange(2.8,g_sig_fit_max);
	  else funcSignal[i]->SetRange(2.75,g_sig_fit_max);
	}

      if(fix_mean[0]>0)
	{
	  if(!(isys==9 && i==0))
	    {
	      funcSignal[i]->FixParameter(1,fix_mean[i]);
	      funcSignal[i]->FixParameter(2,fix_sigma[i]);
	    }
	}
      funcSignal[i]->SetLineColor(2);
      TFitResultPtr ptr = hSignal[i]->Fit(funcSignal[i],"IRS0Q");
      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
      double fit_yield = funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1);
      double fit_yield_err = funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1);
      double fit_mean = funcSignal[i]->GetParameter(1);
      double fit_sigma = fabs(funcSignal[i]->GetParameter(2));
      if(isys==8)
	{
	  TF1 *funcSignalTmp = new TF1(Form("tmp_%s",funcSignal[i]->GetName()),CrystalBall,g_sig_fit_min,g_sig_fit_max,5);
	  for(int p=0; p<5; p++)
	    {
	      funcSignalTmp->SetParameter(p, funcSignal[i]->GetParameter(p));
	      funcSignalTmp->SetParError(p, funcSignal[i]->GetParError(p));
	    }
	  fit_yield = funcSignalTmp->Integral(2.8,3.3)/hSignal[i]->GetBinWidth(1);
	  fit_yield_err = sqrt(fit_yield);
	}
      printf("Fitting = %1.1f +/- %1.1f (%1.1fsigma)\n",fit_yield,fit_yield_err,fit_yield/fit_yield_err);
      hFitYield->SetBinContent(i,fit_yield);
      hFitYield->SetBinError(i,fit_yield_err);
      hMean->SetBinContent(i,funcSignal[i]->GetParameter(1));
      hMean->SetBinError(i,funcSignal[i]->GetParError(1));
      hSigma->SetBinContent(i,fit_sigma);
      hSigma->SetBinError(i,funcSignal[i]->GetParError(2));

      // Extract background matrix
      const int nParameter = nPar;
      double bkg_params[nParameter];
      double bkg_matrix[nParameter*nParameter];
      funcBkg[i] = new TF1(Form("Jpsi_FitBkg_pt%s_%s",pt_Name[i],suffix.Data()),Form("%s",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
      int nSigPar = 3;
      if(isys==8) nSigPar = 5;
      for(int j=0; j<nPar; j++)
	{
	  funcBkg[i]->SetParameter(j,funcSignal[i]->GetParameter(nSigPar+j));
	  funcBkg[i]->SetParError(j,funcSignal[i]->GetParError(nSigPar+j));
	  bkg_params[j] = funcSignal[i]->GetParameter(nSigPar+j);
	}
 
      for(int j=nSigPar; j<nSigPar+nPar; j++)
	{
	  for(int k=nSigPar; k<nSigPar+nPar; k++)
	    {
	      bkg_matrix[(j-nSigPar)*nPar+k-nSigPar] = matrix[j*(nSigPar+nPar)+k];
	    }
	}

      if(0)
	{
	  ptr->GetCovarianceMatrix().Print();
	  for(int im=0; im<(3+nPar)*(3+nPar); im++) cout << matrix[im] << "  ";
	  cout << endl;

	  for(int im=0; im<(3+nPar); im++) cout << parameters[im] << "  ";
	  cout << endl;
	}

      // bin counting
      double low_mass_tmp = fit_mean - 3.5 * fit_sigma;
      double high_mass_tmp = fit_mean + 3.5 * fit_sigma;
      int low_bin = hSignal[i]->FindFixBin(low_mass_tmp+1e-4);
      int high_bin = hSignal[i]->FindFixBin(high_mass_tmp-1e-4);
      double low_bin_mass = hSignal[i]->GetXaxis()->GetBinLowEdge(low_bin);
      double high_bin_mass = hSignal[i]->GetXaxis()->GetBinUpEdge(high_bin);
      double yield_all_err;
      double yield_all = hSignal[i]->IntegralAndError(low_bin,high_bin,yield_all_err);
      double yield_bkg = funcBkg[i]->Integral(low_bin_mass,high_bin_mass) * 1./hSignal[i]->GetBinWidth(1);
      double yield_bkg_err = funcBkg[i]->IntegralError(low_bin_mass,high_bin_mass,bkg_params,bkg_matrix)* 1./hSignal[i]->GetBinWidth(1);
      double yield_sig = yield_all - yield_bkg;
      double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err);
      double all = hSeUL[i]->Integral(hSeUL[i]->FindBin(low_mass_tmp+1e-4),hSeUL[i]->FindBin(high_mass_tmp-1e-4));
      hBinCountYield->SetBinContent(i,yield_sig);
      hBinCountYield->SetBinError(i,yield_sig_err);
      hSigToBkg->SetBinContent(i,yield_sig/(all-yield_sig));
      printf("Count = %1.1f +/- %1.1f (%1.1fsigma)\n",yield_sig,yield_sig_err,yield_sig/yield_sig_err);
      printf("Background = %1.1f +/- %1.1f\n",yield_bkg,yield_bkg_err);
      printf("all = %1.2f, signal = %1.2f\n",all,yield_sig);

      // plotting
      if(i==0) cFit1->cd();
      else     cFit->cd(i);	  
      hSignal[i]->SetMaximum(2*hSignal[i]->GetMaximum());
      if(i<4)  hSignal[i]->SetMaximum(1.2*hSignal[i]->GetMaximum());
      hSignal[i]->Draw();
      funcSignal[i]->Draw("sames");
      funcBkg[i]->SetLineColor(4);
      funcBkg[i]->SetLineStyle(2);
      funcBkg[i]->Draw("sames");
      TPaveText *t = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.05);
      t->Draw();
      TLine *line = GetLine(low_mass_tmp,hSignal[i]->GetMinimum()*1.5,low_mass_tmp,hSignal[i]->GetMaximum()*0.3,1);
      line->Draw();
      TLine *line = GetLine(high_mass_tmp,hSignal[i]->GetMinimum()*1.5,high_mass_tmp,hSignal[i]->GetMaximum()*0.3,1);
      line->Draw();

      t = GetPaveText(0.16,0.3,0.63,0.88,0.045);
      t->SetTextFont(62);
      t->AddText("Fitting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",fit_yield,fit_yield_err,fit_yield/fit_yield_err));
      t->AddText("Counting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",yield_sig,yield_sig_err,yield_sig/yield_sig_err));
      t->SetTextAlign(11);
      t->Draw();
    }
  t = GetPaveText(0.7,0.8,0.3,0.35,0.06);
  t->SetTextFont(62);
  t->AddText("UL-MIX(UL)");
  t->SetTextColor(4);
  cFit->cd(1);
  t->Draw();

  if(savePlot) 
    {
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sSig.Fit_All_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sSig.Fit_Jpsi_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  TF1 *funcMean = new TF1(Form("Jpsi_FitMeanFit_%s",suffix.Data()), "pol0", 0, 15);
  hMean->Fit(funcMean,"IR0Q");
  hMean->SetMarkerStyle(21);
  hMean->GetYaxis()->SetRangeUser(3.06,3.12);
  c = draw1D(hMean,Form("Mean of Gaussian fit to Jpsi signal (%s%%);p_{T} (GeV/c);mean",cent_Name[icent]));
  funcMean->SetLineColor(4);
  funcMean->SetLineStyle(2);
  funcMean->Draw("sames");
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sSig.FitMeanVsPt_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  TF1 *funcSigma = new TF1(Form("Jpsi_FitSigmaFit_%s",suffix.Data()), "pol1", 0, 15);
  hSigma->Fit(funcSigma,"IR0Q");
  hSigma->SetMarkerStyle(21);
  hSigma->GetYaxis()->SetRangeUser(0,0.15);
  c = draw1D(hSigma,Form("Sigma of Gaussian fit to Jpsi signal (%s%%);p_{T} (GeV/c);#sigma",cent_Name[icent]));
  funcSigma->SetLineColor(4);
  funcSigma->SetLineStyle(2);
  funcSigma->Draw("sames");
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/%s%sSig.FitSigmaVsPt_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],suf_title.Data()));
    }

  if(saveHisto)
    {
      cout << outName.Data() << endl;
      TFile *fout = 0;
      if(isetup==0 && icent==0 && gApplyWeight==0) 
	{
	if(isys==0) fout = TFile::Open(outName,"update");
	else        fout = TFile::Open(outNameSys,"update");
	}
      else
	{
	  if(isys==0) fout = TFile::Open(outName,"update");
	  else        fout = TFile::Open(outNameSys,"update");
	}
      hFitYield->Write("",TObject::kOverwrite);
      hBinCountYield->Write("",TObject::kOverwrite);
      hMean->Write("",TObject::kOverwrite);
      funcMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      funcSigma->Write("",TObject::kOverwrite);
      if(isys==0)
	{
	  hFitScaleFactor->Write("",TObject::kOverwrite);
	  hBinCountScaleFactor->Write("",TObject::kOverwrite);
	  hSigToBkg->Write("",TObject::kOverwrite);
	  for(int i=0; i<nPtBins; i++)
	    {
	      hSignalSave[i]->Write("",TObject::kOverwrite);
	      funcSignal[i]->Write("",TObject::kOverwrite);
	      funcBkg[i]->Write("",TObject::kOverwrite);
	      hMixBkg[i]->Write(Form("MEbkg_pt%s_%s",pt_Name[i],suffix.Data()),TObject::kOverwrite);
	      hSeUL[i]->Write(Form("DataUL_pt%s_%s",pt_Name[i],suffix.Data()),TObject::kOverwrite);
	      hSeLS[i]->SetTitle("");
	      hSeLS[i]->Write(Form("DataLS_pt%s_%s",pt_Name[i],suffix.Data()),TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

//================================================
void makeYieldRun13(const int icent = 0, const int isys = 0, int savePlot = 0, int saveHisto = 0)
{
  const char *sys_name[7] = {"","_LargeScale","_SmallScale","_pol1","_LargeFit","_SmallFit","_Rebin"};
  const char *sys_title[7] = {"","Sys.LargeScale.","Sys.SmallScale.","Sys.pol1.","Sys.LargeFit.","Sys.SmallFit.","Sys.Rebin"};

  printf("Process centrality %s%%\n",cent_Name[icent]);
  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  TFile *fin = TFile::Open(fileName,"read");
  cout << fileName << endl;

  TString outName = fileName;
  outName.ReplaceAll(".root",".yield.root");
  TString outNameSys = outName;
  outNameSys.ReplaceAll(".yield.root",".sys.signal.root");

  TH1F *hSeUL[nPtBins];
  TH1F *hSeLS[nPtBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      hSeUL[i] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hSeLS[i] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
    }
  
  int nxpad, nypad;
  if(nPtBins-1==6 || nPtBins-1==5) { nxpad = 3; nypad = 2; }
  if(nPtBins-1<=4) { nxpad = 2; nypad = 2; }

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // jpsi signal
  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",cent_Name[icent]),Form("InvMass_%s",cent_Name[icent]),1100,650);
  cSignal->Divide(nxpad+1,nypad);
  int rebin_signal = 10;
  if(isys==6) rebin_signal = 5;
  for(int i=0; i<nPtBins; i++)
    {
      hSeUL[i]->Rebin(rebin_signal);
      hSeUL[i]->SetTitle("");
      hSeUL[i]->SetMarkerStyle(21);
      hSeUL[i]->SetMarkerColor(2);
      hSeUL[i]->SetLineColor(2);
      hSeUL[i]->GetXaxis()->SetRangeUser(2,4);
      //hSeUL[i]->SetMaximum(1.2*hSeUL[i]->GetMaximum());
      hSeLS[i]->Rebin(rebin_signal);
      cSignal->cd(i+1);
      hSeUL[i]->Draw();
      hSeLS[i]->Draw("sames HIST");
      TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]),0.06);
      t->Draw();
    }
  cSignal->cd(nPtBins+1);
  leg = new TLegend(0.2,0.5,0.5,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->AddEntry(hSeUL[1],"Unlike sign","P");
  leg->AddEntry(hSeLS[1],"Like sign (++)+(--)","L");
  leg->Draw();
  if(savePlot) 
    {
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sJpsiSignal_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sJpsiSignal_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  // Fit signal
  TCanvas *cFit1 = new TCanvas(Form("Fit_Jpsi_%s_All",cent_Name[icent]),Form("Fit_Jpsi_%s_All",cent_Name[icent]),800,650);
  TCanvas *cFit = new TCanvas(Form("Fit_Jpsi_%s",cent_Name[icent]),Form("Fit_Jpsi_%s",cent_Name[icent]),1100,650);
  cFit->Divide(nxpad,nypad);
  TH1F *hSignal[nPtBins];
  TH1F *hSignalSave[nPtBins];
  TF1 *funcSignal[nPtBins];
  TH1F *hMean = new TH1F(Form("Jpsi_FitMean_cent%s%s",cent_Title[icent],sys_name[isys]),Form("Jpsi_FitMean_cent%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TH1F *hSigma = new TH1F(Form("Jpsi_FitSigma_cent%s%s",cent_Title[icent],sys_name[isys]),Form("Jpsi_FitSigma_cent%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TH1F *hFitYield = new TH1F(Form("Jpsi_FitYield_cent%s%s",cent_Title[icent],sys_name[isys]),Form("Jpsi_FitYield_cent%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TH1F *hBinCountYield = new TH1F(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[icent],sys_name[isys]),Form("Jpsi_BinCountYield_cent%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TH1F *hSigToBkg = new TH1F(Form("Jpsi_SigToBkg_cent%s%s",cent_Title[icent],sys_name[isys]),Form("Jpsi_SigToBkg_cent%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TString funcForm1 = "pol3";
  int nPar1 = 4;
  TString funcForm2 = "pol3";
  int nPar2 = 4;
  double fit_min = 2.5, fit_max = 4.0;
  if(isys==3)
    {
      funcForm2 = "pol1"; 
      const int nPar2 = 2;
    }
  if(isys==4)
    {
      fit_min = 2.5; fit_max = 3.9;
    }
  if(isys==5)
    {
      fit_min = 2.7; fit_max = 3.7;
    }

  TString funcForm;
  int nPar;
  for(int i=0; i<nPtBins; i++)
    {
      printf("+++++ %1.0f < pT < %1.0f +++++\n",ptBins_low[i],ptBins_high[i]);

      // Fitting
      hSignal[i] = (TH1F*)hSeUL[i]->Clone(Form("Jpsi_Signal_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      hSignal[i]->GetXaxis()->SetRangeUser(2.5,4);
      for(int bin=1; bin<=hSignal[i]->GetNbinsX(); bin++)
	{
	  if(hSignal[i]->GetBinContent(bin)==0)
	    {
	      hSignal[i]->SetBinContent(bin,0);
	      hSignal[i]->SetBinError(bin,1.4);
	    }
	}
      hSignal[i]->SetLineColor(1);
      hSignal[i]->SetMarkerColor(1);
      hSignalSave[i] = (TH1F*)hSignal[i]->Clone(Form("%s_save",hSignal[i]->GetName()));
      if(i==0)
      	{
	  funcForm = funcForm1; 
	  nPar = nPar1;
      	}
      else
      	{
	  funcForm = funcForm2; 
	  nPar = nPar2;
      	}

      
      TF1 *functmp = new TF1(Form("Fit%s_temp",hSignal[i]->GetName()),Form("gausn(0)+%s(3)",funcForm.Data()),fit_min,fit_max);
      functmp->SetParameter(0,20);
      functmp->FixParameter(1,3.09);
      functmp->FixParameter(2,0.05);
      hSignal[i]->Fit(functmp,"IRS0Q");
      
      funcSignal[i] = new TF1(Form("Fit%s",hSignal[i]->GetName()),Form("gausn(0)+%s(3)",funcForm.Data()),fit_min,fit_max);
      for(int ip=0; ip<functmp->GetNpar(); ip++)
      	{
      	  funcSignal[i]->SetParameter(ip,functmp->GetParameter(ip));
      	}
      funcSignal[i]->SetLineColor(2);
      TFitResultPtr ptr = hSignal[i]->Fit(funcSignal[i],"IRS0");
      printf("Fit signal = %1.1f +/- %1.1f (%1.1fsigma)\n",funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParameter(0)/funcSignal[i]->GetParError(0));
      hFitYield->SetBinContent(i,funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1));
      hFitYield->SetBinError(i,funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1));
      hMean->SetBinContent(i,funcSignal[i]->GetParameter(1));
      hMean->SetBinError(i,funcSignal[i]->GetParError(1));
      hSigma->SetBinContent(i,funcSignal[i]->GetParameter(2));
      hSigma->SetBinError(i,funcSignal[i]->GetParError(2));
      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();

      // bin counting
      double low_mass_tmp = low_mass;
      double high_mass_tmp = high_mass;
      if(i>4) 
	{
	  low_mass_tmp = 2.8;
	  high_mass_tmp = 3.3;
	}
      //double low_mass_tmp = funcSignal[i]->GetParameter(1) - 3*funcSignal[i]->GetParameter(2);
      //double high_mass_tmp = funcSignal[i]->GetParameter(1) + 3*funcSignal[i]->GetParameter(2);
      //cout << low_mass_tmp << " -> " << high_mass_tmp << endl;
      int low_bin = hSignal[i]->FindFixBin(low_mass_tmp+1e-4);
      int high_bin =hSignal[i]->FindFixBin(high_mass_tmp-1e-4);
      //cout << low_bin << " -> " << high_bin << endl;
      double error;
      double count = hSignal[i]->IntegralAndError(low_bin,high_bin,error);
      printf("All = %1.0f +/- %1.1f\n",count,error);
      TF1 *functmp = new TF1(Form("bkg_%s",hSignal[i]->GetName()),Form("%s",funcForm.Data()),fit_min,fit_max);
      functmp->SetLineColor(4);
      const int nParameter = nPar;
      double bkg_params[nParameter];
      for(int j=0; j<nPar; j++)
	{
	  functmp->SetParameter(j,funcSignal[i]->GetParameter(3+j));
	  functmp->SetParError(j,funcSignal[i]->GetParError(3+j));
	  bkg_params[j] = funcSignal[i]->GetParameter(3+j);
	}
      double bkg_matrix[nParameter*nParameter];
      for(int j=3; j<3+nPar; j++)
	{
	  for(int k=3; k<3+nPar; k++)
	    {
	      bkg_matrix[(j-3)*nPar+k-3] = matrix[j*(3+nPar)+k];
	    }
	}
      if(0)
	{
	  for(int im=0; im<nPar*nPar; im++) cout << bkg_matrix[im] << "  ";
	  cout << endl;
	  double *parameters = ptr->GetParams();
	  for(int im=0; im<nPar; im++) cout << bkg_params[im] << "  ";
	  cout << endl;
	}
      double bkg = functmp->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal[i]->GetBinWidth(1);
      double bkg_err = functmp->IntegralError(low_mass_tmp,high_mass_tmp,bkg_params,bkg_matrix) * 1./hSignal[i]->GetBinWidth(1);
      double signal = count - bkg;
      double sig_err = TMath::Sqrt(error*error+bkg_err*bkg_err);
      printf("Fit combined = %1.1f +/- %1.1f -> %1.1f +/- %1.1f (%1.1fsigma)\n",bkg,bkg_err,signal,sig_err,signal/sig_err);

      // fit background only
      TF1 *func_tmp = 0;
      if(funcForm=="pol1") func_tmp = new TF1(Form("func_tmp_%d",i),Polynomial1,fit_min,fit_max,nPar);
      if(funcForm=="pol3") func_tmp = new TF1(Form("func_tmp_%d",i),Polynomial3,fit_min,fit_max,nPar);
      for(int ip=0; ip<func_tmp->GetNpar(); ip++)
      	{
      	  //func_tmp->SetParameter(ip,functmp->GetParameter(ip+3));
      	}
      TH1F *hBkgSig = (TH1F*)hSignal[i]->Clone(Form("%s_fit_bkg",hSignal[i]->GetName()));
      TFitResultPtr ptr = hBkgSig->Fit(func_tmp,"IR0QS");
      TF1 *bkgfunc_tmp = new TF1(Form("bkgfunc_tmp_%d",i),Form("%s",funcForm.Data()),fit_min,fit_max);
      for(int j=0; j<nPar; j++)
      	{
      	  bkgfunc_tmp->SetParameter(j,func_tmp->GetParameter(j));
      	  bkgfunc_tmp->SetParError(j,func_tmp->GetParError(j));
      	}
      double bkg_2 = bkgfunc_tmp->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal[i]->GetBinWidth(1);
      double bkg_err_2 = bkgfunc_tmp->IntegralError(low_mass_tmp,high_mass_tmp,ptr->GetParams(), ptr->GetCovarianceMatrix().GetMatrixArray())* 1./hSignal[i]->GetBinWidth(1);
      double signal_2 = count-bkg_2;
      double sig_err_2 = sqrt(error*error+bkg_err_2*bkg_err_2);
      printf("Fit background = %1.1f +/- %1.1f -> %1.1f +/- %1.1f (%1.1fsigma)\n",bkg_2,bkg_err_2,signal_2,sig_err_2,signal_2/sig_err_2);
      hBinCountYield->SetBinContent(i,signal_2);
      hBinCountYield->SetBinError(i,sig_err_2);
      double all = hSeUL[i]->Integral(hSeUL[i]->FindBin(low_mass_tmp+1e-4),hSeUL[i]->FindBin(high_mass_tmp-1e-4));
      hSigToBkg->SetBinContent(i,signal_2/(all-signal_2));


      // plotting
      if(i==0) cFit1->cd();
      else     cFit->cd(i);	  
      hSignal[i]->SetMaximum(2*hSignal[i]->GetMaximum());
      if(i<4)  hSignal[i]->SetMaximum(1.2*hSignal[i]->GetMaximum());
      hSignal[i]->Draw();
      funcSignal[i]->Draw("sames");
      functmp->Draw("sames");
      TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]),0.06);
      t->Draw();
      TLine *line = GetLine(low_mass_tmp,hSignal[i]->GetMinimum()*1,low_mass_tmp,hSignal[i]->GetMaximum()*0.5,1);
      line->Draw();
      TLine *line = GetLine(high_mass_tmp,hSignal[i]->GetMinimum()*1,high_mass_tmp,hSignal[i]->GetMaximum()*0.5,1);
      line->Draw();

      t = GetPaveText(0.16,0.3,0.63,0.88,0.045);
      t->SetTextFont(62);
      t->AddText("Fitting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParameter(0)/funcSignal[i]->GetParError(0)));
      t->AddText("Counting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",signal_2,sig_err_2,signal_2/sig_err_2));
      t->SetTextAlign(11);
      t->Draw();

      func_tmp->SetLineColor(6);
      func_tmp->SetLineStyle(2);
      func_tmp->Draw("sames");
    }



  if(savePlot) 
    {
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_All_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_All_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  hMean->SetMarkerStyle(21);
  hMean->GetYaxis()->SetRangeUser(3.0,3.12);
  c = draw1D(hMean,Form("Mean of Gaussian fit to Jpsi signal;p_{T} (GeV/c);mean"));
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

    }

  hSigma->SetMarkerStyle(21);
  hSigma->GetYaxis()->SetRangeUser(0,0.12);
  c = draw1D(hSigma,Form("Sigma of Gaussian fit to Jpsi signal;p_{T} (GeV/c);#sigma"));

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  if(saveHisto)
    {
      cout << outName.Data() << endl;
      TFile *fout = 0;
      if(icent==0) 
	{
	if(isys==0) fout = TFile::Open(outName,"recreate");
	else if(isys==1) fout = TFile::Open(outNameSys,"recreate");
	else fout = TFile::Open(outNameSys,"update");
	}
      else
	{
	  if(isys==0) fout = TFile::Open(outName,"update");
	  else        fout = TFile::Open(outNameSys,"update");
	}
      //hFitScaleFactor->Write("",TObject::kOverwrite);
      //hBinCountScaleFactor->Write("",TObject::kOverwrite);
      hMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      hFitYield->Write("",TObject::kOverwrite);
      hBinCountYield->Write("",TObject::kOverwrite);
      hSigToBkg->Write("",TObject::kOverwrite);
      if(isys==0)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      //hAcc[i]->Write("",TObject::kOverwrite);
	      hSignal[i]->Write("",TObject::kOverwrite);
	      hSignalSave[i]->Write("",TObject::kOverwrite);
	      funcSignal[i]->Write("",TObject::kOverwrite);
	      //hMixBkg[i]->Write("",TObject::kOverwrite);
	      hSeUL[i]->Write("",TObject::kOverwrite);
	      hSeLS[i]->Write("",TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

//================================================
void makeYieldRun15(const int icent = 0, int savePlot = 0, int saveHisto = 0)
{
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.JpsiYield.root",run_type),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s.JpsiYield.root",run_type),"read");

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // jpsi signal
  TH1F *hSignal[nPtBins];
  TF1 *funcSignal[nPtBins];
  TH1F *hMean = new TH1F(Form("Jpsi_FitMean_cent%s",cent_Title[icent]),"",nbins,xbins);
  TH1F *hSigma = new TH1F(Form("Jpsi_FitSigma_cent%s",cent_Title[icent]),"",nbins,xbins);
  TH1F *hBinCountYield = new TH1F(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]),"",nbins,xbins);

  TCanvas *c = new TCanvas("Fit_JpsiSignal","Fit_JpsiSignal",1100,700);
  c->Divide(3,3);
  for(int i=0; i<nPtBins; i++)
    {
      printf("+++++ %1.0f < pT < %1.0f +++++\n",ptBins_low[i],ptBins_high[i]);

      // Fitting
      hSignal[i] = (TH1F*)fin->Get(Form("mhMassCent0_Pt%d",i));
      hSignal[i]->GetXaxis()->SetTitleOffset(0.85);
      funcSignal[i] = (TF1*)fin->Get(Form("mfAllCent0_Pt%d",i));
      if(i>0)
	{
	  c->cd(i);
	  hSignal[i]->SetMaximum(1.3*hSignal[i]->GetMaximum());
	  hSignal[i]->SetTitle("");
	  hSignal[i]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]),0.08);
	  t1->Draw();
	}

      //funcSignal[i]->SetParErrors(func1->GetParErrors());
      hBinCountYield->SetBinContent(i,funcSignal[i]->GetParameter(0));
      hBinCountYield->SetBinError(i,funcSignal[i]->GetParError(0));
      hMean->SetBinContent(i,funcSignal[i]->GetParameter(1));
      hMean->SetBinError(i,funcSignal[i]->GetParError(1));
      hSigma->SetBinContent(i,funcSignal[i]->GetParameter(2));
      hSigma->SetBinError(i,funcSignal[i]->GetParError(2));
    }  
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiSignal_cent%s.pdf",run_type,cent_Title[icent]));   
  c = draw1D(hSignal[0]);
  TPaveText *t1 = GetPaveText(0.15,0.3,0.75,0.8,0.04);
  t1->AddText(run_type);
  t1->Draw();
  t1->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiSignalAll_cent%s.pdf",run_type,cent_Title[icent]));   

  TH1F *hInputTmp = (TH1F*)hBinCountYield->Clone(Form("%s_tmp",hBinCountYield->GetName()));
  for(int bin=1; bin<=hInputTmp->GetNbinsX(); bin++)
    {
      double dpt = hInputTmp->GetBinWidth(bin);
      double pt = hInputTmp->GetBinCenter(bin);
      double yield = hInputTmp->GetBinContent(bin);
      double error = hInputTmp->GetBinError(bin);
      hInputTmp->SetBinContent(bin, yield/dpt);
      hInputTmp->SetBinError(bin, error/dpt);
    }
  TF1 *funcInputJpsi = new TF1(Form("fitfunc_%s",hBinCountYield->GetName()),"[0]*exp([1]*x+[2]*x*x)*x*x",0.1,10);
  hInputTmp->Fit(funcInputJpsi,"R0Q");
  hInputTmp->SetMarkerStyle(21);
  c = draw1D(hInputTmp, Form("%s: raw J/#Psi yield as a function of p_{T};p_{T} (GeV/c);dN/dp_{T}",run_type), true, true);
  funcInputJpsi->SetLineStyle(2);
  funcInputJpsi->SetLineColor(2);
  funcInputJpsi->Draw("sames");
  leg = new TLegend(0.15,0.2,0.4,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hInputTmp,"Data","PL");
  leg->AddEntry(funcInputJpsi,"Fit to data: [0]*e^{[1]*x+[2]*x^{2}}*x^{2}","L");
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiYield_cent%s.pdf",run_type,cent_Title[icent]));

  // re-generate the resolution plot
  TFile* fRes = 0x0;
  if(saveHisto) fRes = TFile::Open(Form("Rootfiles/%s.PtResolution.root",run_type),"update");
  else          fRes = TFile::Open(Form("Rootfiles/%s.PtResolution.root",run_type),"read");
  TH1F *mhRelPtResVsPt = (TH1F*)fRes->Get("mhRelPtResVsPt");
  TH1F *mhRelPtShiftVsPt = (TH1F*)fRes->Get("mhRelPtShiftVsPt");
  TH2F *hRelRes = new TH2F(Form("PrimTrkRes_vs_TruePt_cent%s",cent_Title[icent]),"",60,0,12,200,-1,1);
  for(int ibin=1; ibin<=60; ibin++)
    {
      TF1 *funcTmp = new TF1(Form("funcTmp_%d",ibin),"gaus",-1,1);
      double pt = mhRelPtResVsPt->GetBinCenter(ibin);
      double weight = funcInputJpsi->Eval(pt);
      funcTmp->SetParameters(1, mhRelPtShiftVsPt->GetBinContent(ibin), mhRelPtResVsPt->GetBinContent(ibin));
      for(int i=0; i<1e6; i++)
	{
	  hRelRes->Fill(pt, funcTmp->GetRandom(), weight);
	}
    }
  c = draw2D(hRelRes);

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSignal[i]->Write(Form("Jpsi_Signal_pt%s_cent%s",pt_Name[i],cent_Title[icent]),TObject::kOverwrite);
	  funcSignal[i]->Write(Form("Jpsi_FitSig_pt%s_cent%s",pt_Name[i],cent_Title[icent]),TObject::kOverwrite);
	}
      hMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      hBinCountYield->Write("",TObject::kOverwrite);
      funcInputJpsi->Write("",TObject::kOverwrite);

      fRes->cd();
      hRelRes->Write("",TObject::kOverwrite);
    }
}


//================================================
void makeYieldRun16(const int icent = 0, int savePlot = 0, int saveHisto = 0)
{
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open("Rootfiles/Run16_AuAu200.JpsiYield.root","update");
  else          fin = TFile::Open("Rootfiles/Run16_AuAu200.JpsiYield.root","read");

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // jpsi signal
  TH1F *hSignal[nPtBins];
  TF1 *funcSignal[nPtBins];
  TH1F *hMean = new TH1F(Form("Jpsi_FitMean_cent%s",cent_Title[icent]),"",nbins,xbins);
  TH1F *hSigma = new TH1F(Form("Jpsi_FitSigma_cent%s",cent_Title[icent]),"",nbins,xbins);
  TH1F *hBinCountYield = new TH1F(Form("Jpsi_BinCountYield_cent%s",cent_Title[icent]),"",nbins,xbins);

  TCanvas *c = new TCanvas("Fit_JpsiSignal","Fit_JpsiSignal",1100,700);
  c->Divide(3,2);
  for(int i=0; i<nPtBins; i++)
    {
      printf("+++++ %1.0f < pT < %1.0f +++++\n",ptBins_low[i],ptBins_high[i]);

      // Fitting
      if(i==0)
	{
	  hSignal[i] = (TH1F*)fin->Get("hFitbgRawSignal");
	  funcSignal[i] = (TF1*)fin->Get("fSig");
	}
      else
	{
	  hSignal[i] = (TH1F*)fin->Get(Form("hFitbgRawSignalPtBin%d",i-1));
	  funcSignal[i] = (TF1*)fin->Get(Form("fSig_PtBin%d",i-1));
	  c->cd(i);
	  hSignal[i]->SetMaximum(1.3*hSignal[i]->GetMaximum());
	  hSignal[i]->Draw();
	}

      TList *list = (TList*)hSignal[i]->GetListOfFunctions();
      TF1 *func1 = (TF1*)list->At(0);
      funcSignal[i]->SetParErrors(func1->GetParErrors());
      hBinCountYield->SetBinContent(i,funcSignal[i]->GetParameter(0));
      hBinCountYield->SetBinError(i,funcSignal[i]->GetParError(0));
      hMean->SetBinContent(i,funcSignal[i]->GetParameter(1));
      hMean->SetBinError(i,funcSignal[i]->GetParError(1));
      hSigma->SetBinContent(i,funcSignal[i]->GetParameter(2));
      hSigma->SetBinError(i,funcSignal[i]->GetParError(2));
    }  
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiSignal_cent%s.pdf",run_type,cent_Title[icent]));   
  c = draw1D(hSignal[0]);
  TPaveText *t1 = GetPaveText(0.15,0.3,0.75,0.8);
  t1->AddText("Run16 AuAu200");
  t1->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiSignalAll_cent%s.pdf",run_type,cent_Title[icent]));   

  TH1F *hInputTmp = (TH1F*)hBinCountYield->Clone(Form("%s_tmp",hBinCountYield->GetName()));
  for(int bin=1; bin<=hInputTmp->GetNbinsX(); bin++)
    {
      double dpt = hInputTmp->GetBinWidth(bin);
      double pt = hInputTmp->GetBinCenter(bin);
      double yield = hInputTmp->GetBinContent(bin);
      double error = hInputTmp->GetBinError(bin);
      hInputTmp->SetBinContent(bin, yield/dpt);
      hInputTmp->SetBinError(bin, error/dpt);
    }
  TF1 *funcInputJpsi = new TF1(Form("fitfunc_%s",hBinCountYield->GetName()),"[0]*exp([1]*x+[2]*x*x)*x*x",0.1,10);
  hInputTmp->Fit(funcInputJpsi,"R0Q");
  hInputTmp->SetMarkerStyle(21);
  c = draw1D(hInputTmp, Form("%s: raw J/#Psi yield as a function of p_{T};p_{T} (GeV/c);dN/dp_{T}",run_type), true, true);
  funcInputJpsi->SetLineStyle(2);
  funcInputJpsi->SetLineColor(2);
  funcInputJpsi->Draw("sames");
  leg = new TLegend(0.15,0.2,0.4,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hInputTmp,"Data","P");
  leg->AddEntry(funcInputJpsi,"Fit to data: [0]*e^{[1]*x+[2]*x^{2}}*x","L");
  leg->Draw();
  if(savePlot) 
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/make_JpsiYield/FitJpsiYield_cent%s.pdf",run_type,cent_Title[icent]));

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSignal[i]->Write(Form("Jpsi_Signal_pt%s_cent%s",pt_Name[i],cent_Title[icent]),TObject::kOverwrite);
	  funcSignal[i]->Write(Form("Jpsi_FitSig_pt%s_cent%s",pt_Name[i],cent_Title[icent]),TObject::kOverwrite);
	}
      hMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      hBinCountYield->Write("",TObject::kOverwrite);
      funcInputJpsi->Write("",TObject::kOverwrite);
    }
}

//================================================
double Polynomial3(double *x, double *par)
{
  if(x[0]>rej_low_mass && x[0]<rej_high_mass)
    {
      TF1::RejectPoint();
      return 0;
    }
  
  return par[0]+par[1]*x[0]+par[2]*x[0]**2+par[3]*x[0]**3;
}

//================================================
double Polynomial1(double *x, double *par)
{
  if(x[0]>rej_low_mass && x[0]<rej_high_mass)
    {
      TF1::RejectPoint();
      return 0;
    }
  
  return par[0]+par[1]*x[0];
}

