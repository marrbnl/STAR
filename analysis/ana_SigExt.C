TFile *f = 0x0;
const int year = YEAR;
TString run_cfg_name;
const int useEmbWidth = 1;

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
  //Run14_signal();
  //prod_npart();
  Run14_npart();
  //getJpsiWidthEmbed();
}


//================================================
void getJpsiWidthEmbed(int icent = 0, int savePlot = 1, int saveHisto = 1)
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

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TFile *fdata = TFile::Open(Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TF1 *funcInputJpsi = (TF1*)fdata->Get(Form("Fit_Jpsi_FitYield_cent%s_weight",cent_Title[icent]));

  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"update");
  const char *sysName_scan[3] = {"def","min","max"};
  TH2F *hMassVsPtEmbed[3];

  TH1F *hInvMass[3][nbins];
  for(int i=0; i<3; i++)
    {
      hMassVsPtEmbed[i] = (TH2F*)fscan->Get(Form("hRcJpsiMass_%s_final_%s",cent_Title[icent],sysName_scan[i]));
      for(int j=0; j<nbins; j++)
	{
	  hInvMass[i][j] = (TH1F*)hMassVsPtEmbed[i]->ProjectionY(Form("InvMass_Pt%d_%s",j,sysName_scan[i]));
	  hInvMass[i][j]->Reset("AC");
	  int low_bin = hMassVsPtEmbed[i]->GetXaxis()->FindFixBin(xbins[j]+1e-4);
	  int up_bin = hMassVsPtEmbed[i]->GetXaxis()->FindFixBin(xbins[j+1]-1e-4);
	  for(int bin=low_bin; bin<=up_bin; bin++)
	    {
	      TH1F *htmp = (TH1F*)hMassVsPtEmbed[i]->ProjectionY(Form("Projy_bin%d_%s",bin,sysName_scan[i]),bin,bin);
	      double scale = funcInputJpsi->Integral(hMassVsPtEmbed[i]->GetXaxis()->GetBinLowEdge(bin),
						     hMassVsPtEmbed[i]->GetXaxis()->GetBinUpEdge(bin));
	      hInvMass[i][j]->Add(htmp, scale); 
	    }
	}
    }

  TCanvas *cFit[3];
  TF1 *funcInvMass[3][nbins];
  TH1F *hInvMassMean[3];
  TH1F *hInvMassSigma[3];
  for(int i=0; i<3; i++)
    {
      cFit[i] = new TCanvas(Form("FitInvMass_%s",sysName_scan[i]), Form("FitInvMass_%s",sysName_scan[i]), 1100, 700);
      cFit[i]->Divide(3,3);
      hInvMassMean[i] = new TH1F(Form("SmearEmb_JpsiMean_%s",sysName_scan[i]), "Mean of J/psi peak;p_{T} (GeV/c)", nbins, xbins);
      hInvMassSigma[i] = new TH1F(Form("SmearEmb_JpsiWidth_%s",sysName_scan[i]), "Width of J/psi peak;p_{T} (GeV/c)", nbins, xbins);     
      for(int j=0; j<nbins; j++)
	{
	  funcInvMass[i][j] = new TF1(Form("func_%s",hInvMass[i][j]->GetName()), "gaus", 2.9, 3.3);
	  hInvMass[i][j]->Fit(funcInvMass[i][j], "R0Q");
	  cFit[i]->cd(j+1);
	  hInvMass[i][j]->GetXaxis()->SetRangeUser(2.8,3.4);
	  hInvMass[i][j]->SetMarkerStyle(24);
	  hInvMass[i][j]->Draw();
	  funcInvMass[i][j]->SetLineColor(4);
	  funcInvMass[i][j]->SetLineStyle(2);
	  funcInvMass[i][j]->Draw("sames");	  
	  hInvMassMean[i]->SetBinContent(j+1, funcInvMass[i][j]->GetParameter(1));
	  hInvMassMean[i]->SetBinError(j+1, funcInvMass[i][j]->GetParError(1));
	  hInvMassSigma[i]->SetBinContent(j+1, funcInvMass[i][j]->GetParameter(2));
	  hInvMassSigma[i]->SetBinError(j+1, funcInvMass[i][j]->GetParError(2));
	}
    }

  TList *list = new TList;
  TString legName[3] = {"Default","Lower limit", "Upper limit"};
  list->Add(hInvMassSigma[0]);
  list->Add(hInvMassSigma[1]);
  list->Add(hInvMassSigma[2]);
  TCanvas *c = drawHistos(list,"Emb_JpsiWidth","Jpsi width extracted from smeared embedding",kFALSE,0,30,kTRUE,0,0.14,kFALSE,kTRUE,legName,kTRUE,"Embedding",0.5,0.7,0.68,0.85,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Emb_JpsiWidthVsPt.pdf",run_type));

  // Jpsi width vs. npart
  TCanvas *cFitIntegr = new TCanvas("cFitIntegr", "cFitIntegr", 1100, 700);
  cFitIntegr->Divide(3,2);
  TH1F* hInvMassIntegr[3][nPtBins_npart];
  TH1F *hSigmaIntegr[3];
  TF1 *funcInvMassIntegr[3][nPtBins_npart];
  for(int i=0; i<3; i++)
    {
      hSigmaIntegr[i] = new TH1F(Form("SmearEmb_JpsiWidthIntegr_%s",sysName_scan[i]), "Width of J/psi peak", nPtBins_npart, 0, nPtBins_npart);
       for(int j=0; j<nPtBins_npart; j++)
	{
	  hSigmaIntegr[i]->GetXaxis()->SetBinLabel(j+1, Form("p_{T} > %1.0f GeV/c",ptBins_low_npart[j]));
	  hInvMassIntegr[i][j] = (TH1F*)hMassVsPtEmbed[i]->ProjectionY(Form("InvMass_Pt%s_%s",pt_Name_npart[j],sysName_scan[i]));
	  hInvMassIntegr[i][j]->Reset("AC");
	  int low_bin = hMassVsPtEmbed[i]->GetXaxis()->FindFixBin(ptBins_low_npart[j]+1e-4);
	  int up_bin = hMassVsPtEmbed[i]->GetXaxis()->FindFixBin(ptBins_high_npart[j]-1e-4);
	  for(int bin=low_bin; bin<=up_bin; bin++)
	    {
	      TH1F *htmp = (TH1F*)hMassVsPtEmbed[i]->ProjectionY(Form("Projy2_bin%d_%s",bin,sysName_scan[i]),bin,bin);
	      double scale = funcInputJpsi->Integral(hMassVsPtEmbed[i]->GetXaxis()->GetBinLowEdge(bin),
						     hMassVsPtEmbed[i]->GetXaxis()->GetBinUpEdge(bin));
	      hInvMassIntegr[i][j]->Add(htmp, scale); 
	    }

	  funcInvMassIntegr[i][j] = new TF1(Form("func_%s",hInvMassIntegr[i][j]->GetName()), "gaus", 3.0, 3.2);
	  hInvMassIntegr[i][j]->Fit(funcInvMassIntegr[i][j], "R0Q");
	  hSigmaIntegr[i]->SetBinContent(j+1, funcInvMassIntegr[i][j]->GetParameter(2));
	  hSigmaIntegr[i]->SetBinError(j+1, funcInvMassIntegr[i][j]->GetParError(2));
	  cFitIntegr->cd(j*3+i+1);
	  hInvMassIntegr[i][j]->SetMarkerStyle(24);
	  hInvMassIntegr[i][j]->GetXaxis()->SetRangeUser(2.8, 3.4);
	  hInvMassIntegr[i][j]->SetTitle("");
	  hInvMassIntegr[i][j]->Draw();
	  TPaveText *t1 = GetTitleText(Form("p_{T} > %1.0f GeV/c (%s)",ptBins_low_npart[j],sysName_scan[i]));
	  t1->Draw();
	  funcInvMassIntegr[i][j]->SetLineColor(4);
	  funcInvMassIntegr[i][j]->SetLineStyle(2);
	  funcInvMassIntegr[i][j]->Draw("sames");	 
	}
    }
  list->Clear();
  list->Add(hSigmaIntegr[0]);
  list->Add(hSigmaIntegr[1]);
  list->Add(hSigmaIntegr[2]);
  TCanvas *c = drawHistos(list,"Emb_JpsiWidthVsNpart","Jpsi width extracted from smeared embedding",kFALSE,0,30,kTRUE,0,0.14,kFALSE,kTRUE,legName,kTRUE,"Embedding",0.5,0.7,0.68,0.85,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Emb_JpsiWidthVsCent.pdf",run_type));


  if(saveHisto)
    {
      fscan->cd();
      for(int i=0; i<3; i++)
	{
	  hInvMassMean[i]->Write("",TObject::kOverwrite);
	  hInvMassSigma[i]->Write("",TObject::kOverwrite);
	  hSigmaIntegr[i]->Write("",TObject::kOverwrite);
	}
    }	  
}

//================================================
void prod_pt()
{
  for(int t=0; t<5; t++)
    {
      for(int i=0; i<nCentBins_pt; i++)
	{
	  for(int j=0; j<12; j++)
	    {
	      if(t>0 && j>0) continue;
	      Run14_signal(t,i,j,1,1);
	    }
	}
    }
}

//===============================================
void Run14_signal(const int isetup = 0, const int icent = 0, const int isys = 11, int savePlot = 0, int saveHisto = 0)
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
  const int nSys = 12;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc1","_BkgFunc2","_LargeFit","_SmallFit","_SigFunc","_FixSig","_FixSigUp"};
  const char *sys_title[nSys] = {"","Sys.LargeScale.","Sys.SmallScale.","Sys.ScaleFit.","Sys.Rebin.","Sys.BkgFunc1.","Sys.BkgFunc2.","Sys.LargeFit.","Sys.SmallFit.","Sys.SigFunc.","Sys.FixSig.","Sys.FixSigUp."};
  const TString suffix = Form("cent%s%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup],sys_name[isys]);
  const TString suf_title = Form("cent%s%s%s",cent_Title[icent],gWeightName[gApplyWeight],gTrgSetupName[isetup]);

  // open input and output files
  TString fileName   = Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);
  TFile *fin = TFile::Open(fileName,"read");

  // get the Jpsi width from embedding
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type),"read");
  TH1F *hEmbJpsiWidth[3];
  hEmbJpsiWidth[0] = (TH1F*)fscan->Get("SmearEmb_JpsiWidth_def");
  hEmbJpsiWidth[1] = (TH1F*)fscan->Get("SmearEmb_JpsiWidth_min");
  hEmbJpsiWidth[2] = (TH1F*)fscan->Get("SmearEmb_JpsiWidth_max");

  printf("\n===== Running configuration =====\n");
  printf("Input filename: %s\n",fileName.Data());
  printf("Histogram name: %s\n",suffix.Data());
  printf("Centrality: %s\n",cent_Name[icent]);
  printf("==================================\n\n");

  // global setup
  double  g_mix_scale_low = 2.7;
  double  g_mix_scale_high = 3.8;
  TString g_mix_func = "pol1";
  double  g_bin_width_1 = 0.04;
  double  g_bin_width_2 = 0.05;
  TString g_func1 = "pol3";
  int     g_func1_npar = 5;
  TString g_func2 = "pol1";
  int     g_func2_npar = 2;
  double  g_sig_fit_min = 2.5;
  double  g_sig_fit_max = 4.0;
  if(isys==1) { g_mix_scale_low = 2.5; g_mix_scale_high = 4.0; }
  if(isys==2) { g_mix_scale_low = 2.8; g_mix_scale_high = 3.7; }
  if(isys==3) { g_mix_func = "pol0"; g_mix_scale_low = 2.7; g_mix_scale_high = 3.8;}
  if(isys==4) { g_bin_width_1 = 0.02; g_bin_width_4 = 0.02;}
  if(isys==5) { g_func1 = "pol2"; g_func1_npar = 3; g_func2 = "pol0"; g_func2_npar = 1; }
  if(isys==6) { g_func1 = "pol4"; g_func1_npar = 5; g_func2 = "pol2"; g_func2_npar = 3; g_sig_fit_min = 2.55; }
  if(isys==7) { g_sig_fit_min = 2.3; g_sig_fit_max = 4.2; }
  if(isys==8) { g_sig_fit_min = 2.6; g_sig_fit_max = 3.8; }
  if(isys==9) { g_sig_fit_min = 2.5; g_sig_fit_max = 4.0; }


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

      if(i>0)
	{
	  // plotting
	  cScaling->cd(i);
	  hMixScale[i]->SetTitle("");
	  hMixScale[i]->SetMarkerStyle(21);
	  if(i==1) hMixScale[i]->GetXaxis()->SetRangeUser(2.7,4);
	  else hMixScale[i]->GetXaxis()->SetRangeUser(2.5,4);
	  hMixScale[i]->SetMaximum(1.1*hMixScale[i]->GetMaximum());
	  //if(i<5) hMixScale[i]->GetYaxis()->SetRangeUser(0.6e-3, 0.75e-3);
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
  t = GetPaveText(0.2,0.35,0.5,0.6,0.06);
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
      double g_bin_width_tmp = g_bin_width_1;
      if(i>=2) g_bin_width_tmp = g_bin_width_2;
      hAcc[i] = (TH1F*)hMixUL[i]->Clone(Form("Mix_Acceptance_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone2",hMixLS[i]->GetName()));
      hAcc[i]->Rebin(g_bin_width_tmp/hAcc[i]->GetBinWidth(1));
      htmp->Rebin(g_bin_width_tmp/htmp->GetBinWidth(1));
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
      
      double g_bin_width_tmp = g_bin_width_1;
      if(i>=2) g_bin_width_tmp = g_bin_width_2;
      hSeUL[i]->Rebin(g_bin_width_tmp/hSeUL[i]->GetBinWidth(1));
      hSeUL[i]->SetTitle("");
      hSeUL[i]->SetMarkerStyle(21);
      hSeUL[i]->SetMarkerColor(2);
      hSeUL[i]->SetLineColor(2);
      hSeUL[i]->GetXaxis()->SetRangeUser(2.5,4);
      hSeLS[i]->Rebin(hSeUL[i]->GetBinWidth(1)/hSeLS[i]->GetBinWidth(1));
      hMixBkg[i]->SetLineColor(4);
      hMixBkg[i]->Rebin(hSeUL[i]->GetBinWidth(1)/hMixBkg[i]->GetBinWidth(1));
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
  if(isys==9 || isys==10)
    {
      TFile *fFit = TFile::Open(outName,"read");
      TH1F *h1 = (TH1F*)fFit->Get(Form("Jpsi_FitMean_cent%s%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0],sys_name[0]));
      TH1F *h2 = (TH1F*)fFit->Get(Form("Jpsi_FitSigma_cent%s%s%s%s",cent_Title[0],gWeightName[gApplyWeight],gTrgSetupName[0],sys_name[0]));
      for(int bin=0; bin<=h1->GetNbinsX(); bin++)
  	{
  	  fix_mean[bin] = h1->GetBinContent(bin);
  	  fix_sigma[bin] = h2->GetBinContent(bin);
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
      if(i<1)
      	{
	  funcForm = g_func1; 
	  nPar = g_func1_npar;
      	}
      else
      	{
	  funcForm = g_func2; 
	  nPar = g_func2_npar;
      	}
      if(isys==3)
	{
	  if(icent<=1 && i>=2 && i<=4) 
	    {
	      funcForm = g_func1; 
	      nPar = g_func1_npar;
	    }
	}

      if(isys!=9)
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
	  if(i>=9) funcSignal[i]->SetRange(2.6,3.7);
	  funcSignal[i]->SetParameter(0,hSignal[i]->GetMaximum());
	  funcSignal[i]->SetParameter(1,fix_mean[i]);
	  funcSignal[i]->SetParameter(2,fix_sigma[i]);
	  funcSignal[i]->SetParameter(3,0.1);
	  funcSignal[i]->SetParameter(4,0.2);
	  if(i==9 && icent==0)
	    {
	      funcSignal[i]->SetParameter(1,3.09);
	      funcSignal[i]->SetParameter(2,0.1);
	    }
	  if(i==8 && icent==2)
	    {
	      funcSignal[i]->SetParameter(0, 10);
	      funcSignal[i]->SetParameter(1,3.09);
	      funcSignal[i]->SetParameter(2,0.05);
	    }
	}
      if(i==1)
	{
	  funcSignal[i]->SetParameter(2,0.15);
	  if(isys==7) funcSignal[i]->SetRange(2.7,g_sig_fit_max);
	  else if(isys==8) funcSignal[i]->SetRange(2.8,g_sig_fit_max);
	  else funcSignal[i]->SetRange(2.75,g_sig_fit_max);
	  if(icent==4 && isys==8) funcSignal[i]->SetParameter(2,0.01);
	}
      if(useEmbWidth && i>0)
	{
	  funcSignal[i]->FixParameter(2,hEmbJpsiWidth[0]->GetBinContent(i));
	}

      if(isys==10)
	{
	  if(useEmbWidth)
	    {
	      if(i>0) funcSignal[i]->FixParameter(2,hEmbJpsiWidth[1]->GetBinContent(i));
	    }
	  else
	    {
	      funcSignal[i]->FixParameter(1,fix_mean[i]);
	      funcSignal[i]->FixParameter(2,fix_sigma[i]);
	    }
	}
      if(isys==11)
	{
	  if(i>0) funcSignal[i]->FixParameter(2,hEmbJpsiWidth[2]->GetBinContent(i));
	}
      funcSignal[i]->SetLineColor(2);
      TFitResultPtr ptr = hSignal[i]->Fit(funcSignal[i],"IRS0Q");
      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
      double fit_yield = funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1);
      double fit_yield_err = funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1);
      double fit_mean = funcSignal[i]->GetParameter(1);
      double fit_sigma = fabs(funcSignal[i]->GetParameter(2));
      if(isys==9)
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
      if(isys==9) nSigPar = 5;
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
      int low_bin = hSignal[i]->FindFixBin(low_mass_tmp+1e-4) + 1;
      int high_bin = hSignal[i]->FindFixBin(high_mass_tmp-1e-4) - 1;
      double low_bin_mass = hSignal[i]->GetXaxis()->GetBinLowEdge(low_bin);
      double high_bin_mass = hSignal[i]->GetXaxis()->GetBinUpEdge(high_bin);
      double yield_all_err;
      double yield_all = hSignal[i]->IntegralAndError(low_bin,high_bin,yield_all_err);
      double yield_bkg = funcBkg[i]->Integral(low_bin_mass,high_bin_mass) * 1./hSignal[i]->GetBinWidth(1);
      double yield_bkg_err = funcBkg[i]->IntegralError(low_bin_mass,high_bin_mass,bkg_params,bkg_matrix)* 1./hSignal[i]->GetBinWidth(1);
      TF1 *funcJpsi_temp = new TF1(Form("funcJpsi_temp_%d",i),"gausn");
      for(int par=0; par<3; par++)
	{
	  funcJpsi_temp->SetParameter(par, funcSignal[i]->GetParameter(par));
	  funcJpsi_temp->SetParError(par, funcSignal[i]->GetParError(par));
	}
      double efficiency = funcJpsi_temp->Integral(low_bin_mass,high_bin_mass)/funcJpsi_temp->GetParameter(0);
      double yield_sig = (yield_all - yield_bkg)/efficiency;
      double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err)/efficiency;
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
      TLine *line = GetLine(low_bin_mass,hSignal[i]->GetMinimum()*1.5,low_bin_mass,hSignal[i]->GetMaximum()*0.3,1);
      line->Draw();
      TLine *line = GetLine(high_bin_mass,hSignal[i]->GetMinimum()*1.5,high_bin_mass,hSignal[i]->GetMaximum()*0.3,1);
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
  double meanFit_max = 15;
  if(icent==2) meanFit_max = 10;
  if(icent==3) meanFit_max = 10;
  if(icent==4) meanFit_max = 6;
  TF1 *funcMean = new TF1(Form("Jpsi_FitMeanFit_%s",suffix.Data()), "pol0", 0, meanFit_max);
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

  TF1 *funcSigma = new TF1(Form("Jpsi_FitSigmaFit_%s",suffix.Data()), "pol1", 0, meanFit_max);
  if(icent==4) funcSigma->FixParameter(1, 0.003683);
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
void prod_npart()
{
  for(int i=0; i<11; i++)
    {
      Run14_npart(i,1,1);
    }
}

//===============================================
void Run14_npart(const int isys = 0, int savePlot = 0, int saveHisto = 0)
{
  // re-assign global constants
  const int nPtBins         = nPtBins_npart;
  const double* ptBins_low  = ptBins_low_npart;
  const double* ptBins_high = ptBins_high_npart;
  const char** pt_Name      = pt_Name_npart;
  const int* nCentBins      = nCentBins_npart; 
  const int* centBins_low   = centBins_low_npart;
  const int* centBins_high  = centBins_high_npart;
  const char** cent_Name    = cent_Name_npart;
  const char** cent_Title   = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  // prepare name and title
  const int nSys = 11;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc1","_BkgFunc2","_LargeFit","_SmallFit","_SigFunc","_FixSig"};
  const char *sys_title[nSys] = {"","Sys.LargeScale.","Sys.SmallScale.","Sys.ScaleFit.","Sys.Rebin.","Sys.BkgFunc1.","Sys.BkgFunc2.","Sys.LargeFit.","Sys.SmallFit.","Sys.SigFunc.","Sys.FixSig."};
  const TString suffix = Form("%s%s",gWeightName[gApplyWeight],sys_name[isys]);

  // open input and output files
  TString fileName   = Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);
  TFile *fin = TFile::Open(fileName,"read");
  printf("\n===== Running configuration =====\n");
  printf("Input filename: %s\n",fileName.Data());
  printf("Histogram name: %s\n",suffix.Data());
  printf("==================================\n\n");

  // global setup
  double  g_mix_scale_low = 2.7;
  double  g_mix_scale_high = 3.8;
  TString g_mix_func = "pol1";
  double  g_bin_width = 0.04;
  TString g_func1 = "pol3";
  int     g_func1_npar = 4;
  TString g_func2 = "pol1";
  int     g_func2_npar = 2;
  double  g_sig_fit_min = 2.5;
  double  g_sig_fit_max = 4.0;
  if(isys==1) { g_mix_scale_low = 2.5; g_mix_scale_high = 4.0; }
  if(isys==2) { g_mix_scale_low = 2.8; g_mix_scale_high = 3.7; }
  if(isys==3) { g_mix_func = "pol0";}
  if(isys==4) { g_bin_width = 0.02; }
  if(isys==5) { g_func1 = "pol2"; g_func1_npar = 3; g_func2 = "pol0"; g_func2_npar = 1; }
  if(isys==6) { g_func1 = "pol4"; g_func1_npar = 5; g_func2 = "pol2"; g_func2_npar = 3; }
  if(isys==7) { g_sig_fit_min = 2.3; g_sig_fit_max = 4.2; }
  if(isys==8) { g_sig_fit_min = 2.6; g_sig_fit_max = 3.8; }
  if(isys==9) { g_sig_fit_min = 2.5; g_sig_fit_max = 4.0; g_bin_width = 0.04;}

  // get histograms
  if(isys==0)
    {
      const int nSetup = gNTrgSetup;
    }
  else
    {
      const int nSetup = 1;
    }
  TH1F *hSeUL[nPtBins][kNCent][nSetup];
  TH1F *hSeLS[nPtBins][kNCent][nSetup];
  TH1F *hMixUL[nPtBins][kNCent];
  TH1F *hMixLS[nPtBins][kNCent];
  for(Int_t i=0; i<nPtBins; i++)
    {
      for(int k=0; k<nCentBins[i]; k++)
	{
	  for(int t=0; t<nSetup; t++)
	    {
	      hSeUL[i][k][t] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s%s%s",pt_Name[i],cent_Name[i*kNCent+k],gWeightName[gApplyWeight],gTrgSetupName[t]));
	      hSeLS[i][k][t] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s%s%s",pt_Name[i],cent_Name[i*kNCent+k],gWeightName[gApplyWeight],gTrgSetupName[t]));
	    }
	  hMixUL[i][k] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[i*kNCent+k]));
	  hMixLS[i][k] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[i*kNCent+k]));
	}
    }  
  int nxpad = 4, nypad = 2;

  // mixed event scaling
  TList *list = new TList;
  TString legName[2] = {"Fitting","Bin Counting"};
  TH1F *hMixScale[nPtBins][kNCent][nSetup];
  TF1 *funcScale[nPtBins][kNCent][nSetup];
  TH1F *hFitScaleFactor[nPtBins][nSetup];
  TH1F *hBinCountScaleFactor[nPtBins][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<nSetup; t++)
	{
	  TString tmpName = Form("pt%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t]);
	  hFitScaleFactor[i][t] = new TH1F(Form("FitScaleFactor_%s",tmpName.Data()),"Mixed-event scale factor from fitting",nCentBins[i],0,nCentBins[i]);
	  hBinCountScaleFactor[i][t] = new TH1F(Form("BinCountScaleFactor_%s",tmpName.Data()),"Mixed-event scale factor from bin counting",nCentBins[i],0,nCentBins[i]);
	  for(int bin=1; bin<=nCentBins[i]; bin++)
	    {
	      hFitScaleFactor[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	      hBinCountScaleFactor[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	    }
	  TCanvas *cScaling = new TCanvas(Form("mix_scale_%s",tmpName.Data()),Form("mix_scale_%s",tmpName.Data()),1100,650);
	  cScaling->Divide(nxpad,nypad);

	  for(int k=0; k<nCentBins[i]; k++)
	    {
	      hMixScale[i][k][t] = (TH1F*)hSeLS[i][k][t]->Clone(Form("%s_MixScale",hSeLS[i][k][t]->GetName()));
	      hMixScale[i][k][t]->Rebin(2);
	      TH1F *htmp = (TH1F*)hMixLS[i][k]->Clone(Form("%s_%d_clone",hMixLS[i][k]->GetName(),t));
	      htmp->Rebin(int(hMixScale[i][k][t]->GetBinWidth(1)/hMixLS[i][k]->GetBinWidth(1)));
	      hMixScale[i][k][t]->Divide(htmp);

	      // fitting method
	      funcScale[i][k][t] = new TF1(Form("Fit_%s",hMixScale[i][k][t]->GetName()),g_mix_func,g_mix_scale_low,g_mix_scale_high);
	      hMixScale[i][k][t]->Fit(funcScale[i][k][t],"IR0Q");
	      hFitScaleFactor[i][t]->SetBinContent(k+1,funcScale[i][k][t]->GetParameter(0));
	      hFitScaleFactor[i][t]->SetBinError(k+1,funcScale[i][k][t]->GetParError(0));

	      // bin counting method
	      double se = 0, se_err = 0, me = 0, me_err = 0;
	      int low_bin = hSeLS[i][k][t]->FindFixBin(g_mix_scale_low+1e-4);
	      int high_bin = hSeLS[i][k][t]->FindFixBin(g_mix_scale_high-1e-4);
	      for(int bin=low_bin; bin<=high_bin; bin++)
		{
		  se += hSeLS[i][k][t]->GetBinContent(bin);
		  se_err += TMath::Power(hSeLS[i][k][t]->GetBinError(bin),2);
		}
	      
	      int low_bin_me = hMixLS[i][k]->FindFixBin(g_mix_scale_low+1e-4);
	      int high_bin_me = hMixLS[i][k]->FindFixBin(g_mix_scale_high-1e-4);
	      for(int bin=low_bin_me; bin<=high_bin_me; bin++)
		{
		  me += hMixLS[i][k]->GetBinContent(bin);
		  me_err += TMath::Power(hMixLS[i][k]->GetBinError(bin),2);
		}
	      se_err = TMath::Sqrt(se_err);
	      me_err = TMath::Sqrt(me_err);
	      double scale = se/me;
	      double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);
	      hBinCountScaleFactor[i][t]->SetBinContent(k+1,scale);
	      hBinCountScaleFactor[i][t]->SetBinError(k+1,scale_error);
	      
	      // plotting
	      cScaling->cd(k+1);
	      hMixScale[i][k][t]->SetTitle("");
	      hMixScale[i][k][t]->SetMarkerStyle(21);
	      hMixScale[i][k][t]->GetXaxis()->SetRangeUser(2.5,4);
	      hMixScale[i][k][t]->SetMaximum(1.5*hMixScale[i][k][t]->GetMaximum());
	      hMixScale[i][k][t]->SetMinimum(0.5*hMixScale[i][k][t]->GetMinimum());
	      hMixScale[i][k][t]->Draw();
	      funcScale[i][k][t]->SetLineColor(2);
	      funcScale[i][k][t]->Draw("sames");
	      TPaveText *text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],cent_Name[i*kNCent+k]),0.05);
	      text->Draw();
	      text = GetPaveText(0.16,0.3,0.7,0.85,0.05);
	      text->SetTextFont(62);
	      text->AddText(Form("Fit: %4.3e",funcScale[i][k][t]->GetParameter(0)));
	      text->AddText(Form("Count: %4.3e",scale));
	      text->SetTextAlign(11);
	      text->Draw();
	    }
	  text = GetPaveText(0.15,0.7,0.15,0.3,0.05);
	  text->SetTextFont(62);
	  text->AddText(Form("Run14_AuAu_200%s",gTrgSetupTitle[t]));
	  text->AddText("LS: SE/ME");
	  text->SetTextAlign(11);
	  text->SetTextColor(4);
	  cScaling->cd(1);
	  text->Draw();
	  if(savePlot && isys==0)
	    { 	 
	      cScaling->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sBkg.FitSEtoME_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	    }

	  list->Add(hFitScaleFactor[i][t]);
	  list->Add(hBinCountScaleFactor[i][t]);
	  c = drawHistos(list,Form("MixScale_%s",tmpName.Data()),Form("Scale factor for mixed events (p_{T} > %1.1f GeV/c)",ptBins_low[i]),kFALSE,0,30,kTRUE,0,hBinCountScaleFactor[i][t]->GetMaximum()*1.2,kFALSE,kTRUE,legName,kTRUE,Form("AuAu_200%s",gTrgSetupTitle[t]),0.5,0.7,0.68,0.85,kTRUE);
	  list->Clear();
	  if(savePlot && isys==0)
	    { 
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sBkg.ScaleFactor_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	    }
	}
    }

  // jpsi signal
  TH1F *hMixBkg[nPtBins][kNCent][nSetup];
  TH1F *hSignal[nPtBins][kNCent][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<nSetup; t++)
	{
	  TString tmpName = Form("pt%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t]);
	  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",tmpName.Data()),Form("InvMass_%s",tmpName.Data()),1100,650);
	  cSignal->Divide(nxpad,nypad);

	  for(int k=0; k<nCentBins[i]; k++)
	    {
	      hMixBkg[i][k][t] = (TH1F*)hMixUL[i][k]->Clone(Form("MixBkg_cent%s%s",cent_Name[i*kNCent+k],tmpName.Data()));
	      hMixBkg[i][k][t]->Scale(hBinCountScaleFactor[i][t]->GetBinContent(k+1));
	      hSeUL[i][k][t]->Rebin(g_bin_width/hSeUL[i][k][t]->GetBinWidth(1));
	      hSeUL[i][k][t]->SetTitle("");
	      hSeUL[i][k][t]->SetMarkerStyle(21);
	      hSeUL[i][k][t]->SetMarkerColor(2);
	      hSeUL[i][k][t]->SetLineColor(2);
	      hSeUL[i][k][t]->GetXaxis()->SetRangeUser(2.5,4);
	      hSeLS[i][k][t]->Rebin(g_bin_width/hSeLS[i][k][t]->GetBinWidth(1));
	      hMixBkg[i][k][t]->SetLineColor(4);
	      hMixBkg[i][k][t]->Rebin(g_bin_width/hMixBkg[i][k][t]->GetBinWidth(1));
	      cSignal->cd(nCentBins[i]-k);
	      hSeUL[i][k][t]->Draw("PE");
	      hSeLS[i][k][t]->Draw("sames HIST");
	      hMixBkg[i][k][t]->Draw("sames HIST");
	      text = GetTitleText(Form("p_{T} > %1.0f GeV/c (%s%%)",ptBins_low[i],cent_Name[i*kNCent+k]),0.06);
	      text->Draw();
	      hSignal[i][k][t] = (TH1F*)hSeUL[i][k][t]->Clone(Form("JpsiSignal_cent%s%s",cent_Name[i*kNCent+k],tmpName.Data()));
	      for(int bin=1; bin<=hSignal[i][k][t]->GetNbinsX(); bin++)
		{
		  if(hSignal[i][k][t]->GetBinContent(bin)==0)
		    {
		      hSignal[i][k][t]->SetBinContent(bin,0);
		      hSignal[i][k][t]->SetBinError(bin,1.4);
		    }
		}
	      hSignal[i][k][t]->Add(hMixBkg[i][k][t],-1);
	      hSignal[i][k][t]->SetLineColor(1);
	      hSignal[i][k][t]->SetMarkerColor(1);
	    }
	  cSignal->cd(1);
	  leg = new TLegend(0.55,0.65,0.75,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->AddEntry(hSeUL[i][1][t],"Unlike sign","P");
	  leg->AddEntry(hSeLS[i][1][t],"Like sign","L");
	  leg->AddEntry(hMixBkg[i][1][t],"Mix UL","L");
	  leg->Draw();
	  if(savePlot && isys==0) 
	    {
	      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sSig.ULvsLS_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	    }
	}
    }

  // extract yields from 0-80%
  TH1F *hJpsiWeight[nPtBins];
  for(int i=0; i<nPtBins; i++)
    hJpsiWeight[i] = 0x0;

  if(isys==0)
    {
      for(int i=0; i<1; i++)
	{
	  c = new TCanvas(Form("Fit_InvMass_All_pt%d",i),Form("Fit_InvMass_All_pt%d",i),1100,650);
	  c->Divide(2,2);
	  hJpsiWeight[i] = new TH1F(Form("hJpsiYieldVsProd_pt%s_0080",pt_Name[i]),Form("hJpsiYieldVsProd_pt%s_0080",pt_Name[i]),4,0,4);
	  for(int t=1; t<nSetup; t++)
	    {
	      c->cd(t);
	      TH1F *htmp = (TH1F*)hSignal[i][0][t]->Clone(Form("hSignal_All_pt%s%s",pt_Name[i],gTrgSetupName[t]));
	      htmp->Reset();
	      for(int k=0; k<nCentBins[i]; k++)
		{
		  htmp->Add(hSignal[i][k][t]);
		}
	      TF1 *functmp = new TF1(Form("FitAll_pt%s_P%d",pt_Name[i],t-1),Form("gausn(0)+pol3(3)"),g_sig_fit_min,g_sig_fit_max);
	      functmp->SetParameter(0,100);
	      functmp->SetParameter(1,3.09);
	      functmp->SetParameter(2,0.1);
	      functmp->SetLineColor(2);
	      htmp->Fit(functmp,"IRS0Q");
	      htmp->SetMaximum(1.3*htmp->GetMaximum());
	      htmp->Draw("P");
	      functmp->SetLineColor(2);
	      functmp->Draw("sames");
	      double fit_yield = functmp->GetParameter(0)/htmp->GetBinWidth(1);
	      double fit_yield_err = functmp->GetParError(0)/htmp->GetBinWidth(1);
	      hJpsiWeight[i]->SetBinContent(t,fit_yield);
	      hJpsiWeight[i]->SetBinError(t,fit_yield_err);
	      TPaveText *t1 = GetTitleText(Form("%s%s",run_type,gTrgSetupTitle[t]),0.06);
	      t1->Draw();
	      t1 = GetPaveText(0.15,0.35,0.75,0.85,0.05);
	      t1->AddText(Form("N_{J/#psi} = %1.0f #pm %1.0f",fit_yield,fit_yield_err));
	      t1->Draw();
	    }
	}  
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%sFitJpsiVsProd_cent0080_pt0-15.pdf",run_type,run_cfg_name.Data()));
	}
    }

  double fix_mean[nPtBins] = {0};
  double fix_sigma[nPtBins];
  if(isys==10)
    {
      TFile *fFit = TFile::Open(outName,"read");
      for(int i=0; i<nPtBins; i++)
	{
	  TF1 *func1 = fFit->Get(Form("Jpsi_FitMeanFit_pt%s_weight",pt_Name[i]));
	  fix_mean[i] = func1->GetParameter(0);

	  TF1 *func2 = fFit->Get(Form("Jpsi_FitSigmaFit_pt%s_weight",pt_Name[i]));
	  fix_sigma[i] = func2->GetParameter(0);
	}
      fFit->Close();
    }

  TH1F *hMean[nPtBins][nSetup];
  TF1  *funcMean[nPtBins][nSetup];
  TH1F *hSigma[nPtBins][nSetup];
  TF1  *funcSigma[nPtBins][nSetup];
  TH1F *hChi2[nPtBins][nSetup];
  TH1F *hFitYield[nPtBins][nSetup];
  TH1F *hBinCountYield[nPtBins][nSetup];
  TF1 *funcSignal[nPtBins][kNCent][nSetup];
  TF1 *funcBkg[nPtBins][kNCent][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<1; t++)
	{
	  TString tmpName = Form("pt%s%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t],sys_name[isys]);
	  hMean[i][t] = new TH1F(Form("Jpsi_FitMean_%s",tmpName.Data()),"Mean of Jpsi peak",nCentBins[i],0,nCentBins[i]);
	  hSigma[i][t] = new TH1F(Form("Jpsi_FitSigma_%s",tmpName.Data()),"Sigma of Jpsi peak",nCentBins[i],0,nCentBins[i]);
	  hChi2[i][t] = new TH1F(Form("Jpsi_FitChi2_%s",tmpName.Data()),"Chi2/NDF of fitting Jpsi peak",nCentBins[i],0,nCentBins[i]);
	  hFitYield[i][t] = new TH1F(Form("Jpsi_FitYield_%s",tmpName.Data()),"Jpsi yield from fitting",nCentBins[i],0,nCentBins[i]);
	  hBinCountYield[i][t] = new TH1F(Form("Jpsi_BinCountYield_%s",tmpName.Data()),"Jpsi yield from bin counting",nCentBins[i],0,nCentBins[i]);
	  funcMean[i][t] = new TF1(Form("Jpsi_FitMeanFit_%s",tmpName.Data()),"pol0",0,nCentBins[i]);
	  funcSigma[i][t] = new TF1(Form("Jpsi_FitSigmaFit_%s",tmpName.Data()),"pol0",0,nCentBins[i]);
	  
	  for(int bin=1; bin<=nCentBins[i]; bin++)
	    {
	      hMean[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	      hSigma[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	      hChi2[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	      hFitYield[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	      hBinCountYield[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",cent_Name[i*kNCent+bin-1]));
	    }
	  
	  TCanvas *cFit = new TCanvas(Form("Fit_Jpsi_%s",tmpName.Data()),Form("Fit_Jpsi_%s",tmpName.Data()),1100,650);
	  cFit->Divide(nxpad,nypad);

	  for(int k=0; k<nCentBins[i]; k++)
	    {
	      TString funcForm;
	      int nPar;
	      if(i==0 && k<5)
		{
		  funcForm = g_func1;
		  nPar = g_func1_npar;
		}
	      else
		{
		  funcForm = g_func2;
		  nPar = g_func2_npar;
		}

	      if(isys!=9)
		{
		  funcSignal[i][k][t] = new TF1(Form("Jpsi_FitSig_cent%s_%s",cent_Name[i*kNCent+k],tmpName.Data()),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
		  funcSignal[i][k][t]->SetParameter(0,funcSignal[i][k][t]->GetMaximum());
		  funcSignal[i][k][t]->SetParameter(1,3.1);
		  funcSignal[i][k][t]->SetParameter(2,0.06);
		}
	      else
		{
		  if(funcForm=="pol1") funcSignal[i][k][t] = new TF1(Form("Jpsi_FitSig_cent%s_%s",cent_Name[i*kNCent+k],tmpName.Data()),CrystalBallPlusPol1,g_sig_fit_min,g_sig_fit_max,7);
		  if(funcForm=="pol3") funcSignal[i][k][t] = new TF1(Form("Jpsi_FitSig_cent%s_%s",cent_Name[i*kNCent+k],tmpName.Data()),CrystalBallPlusPol3,g_sig_fit_min,g_sig_fit_max,9);
		  funcSignal[i][k][t]->SetParameter(0,hSignal[i][k][t]->GetMaximum());
		  funcSignal[i][k][t]->SetParameter(1,3.09);
		  funcSignal[i][k][t]->SetParameter(2,0.05);
		  funcSignal[i][k][t]->SetParameter(3,0.1);
		  funcSignal[i][k][t]->SetParameter(4,0.2);
		}
	      if(isys==10)
		{
		  funcSignal[i][k][t]->FixParameter(1,fix_mean[i]);
		  funcSignal[i][k][t]->FixParameter(2,fix_sigma[i]);
		}
	      if(t>0)
		{
		  funcSignal[i][k][t]->FixParameter(1,funcSignal[i][k][0]->GetParameter(1));
		  funcSignal[i][k][t]->FixParameter(2,funcSignal[i][k][0]->GetParameter(2));
		}

	      funcSignal[i][k][t]->SetLineColor(2);
	      TFitResultPtr ptr = hSignal[i][k][t]->Fit(funcSignal[i][k][t],"IRS0Q");
	      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
	      double fit_yield = funcSignal[i][k][t]->GetParameter(0)/hSignal[i][k][t]->GetBinWidth(1);
	      double fit_yield_err = funcSignal[i][k][t]->GetParError(0)/hSignal[i][k][t]->GetBinWidth(1);
	      double fit_mean = funcSignal[i][k][t]->GetParameter(1);
	      double fit_sigma = funcSignal[i][k][t]->GetParameter(2);
	      if(isys==9)
		{
		  TF1 *funcSignalTmp = new TF1(Form("tmp_%s",funcSignal[i][k][t]->GetName()),CrystalBall,g_sig_fit_min,g_sig_fit_max,5);
		  for(int p=0; p<5; p++)
		    {
		      funcSignalTmp->SetParameter(p, funcSignal[i][k][t]->GetParameter(p));
		      funcSignalTmp->SetParError(p, funcSignal[i][k][t]->GetParError(p));
		    }
		  fit_yield = funcSignalTmp->Integral(2.7,3.4)/hSignal[i][k][t]->GetBinWidth(1);
		  fit_yield_err = sqrt(fit_yield);
		}
	      hMean[i][t]->SetBinContent(k+1,funcSignal[i][k][t]->GetParameter(1));
	      hMean[i][t]->SetBinError(k+1,funcSignal[i][k][t]->GetParError(1));
	      hSigma[i][t]->SetBinContent(k+1,funcSignal[i][k][t]->GetParameter(2));
	      hSigma[i][t]->SetBinError(k+1,funcSignal[i][k][t]->GetParError(2));
	      hChi2[i][t]->SetBinContent(k+1,funcSignal[i][k][t]->GetChisquare()/funcSignal[i][k][t]->GetNDF());
	      hChi2[i][t]->SetBinError(k+1, 1e-5);
	      hFitYield[i][t]->SetBinContent(k+1,fit_yield);
	      hFitYield[i][t]->SetBinError(k+1,fit_yield_err);
	      if(t==0)
		{
		  printf("+++++ %1.0f < pT < %1.0f, %s +++++\n",ptBins_low[i],ptBins_high[i],cent_Name[i*kNCent+k]);
		  printf("# of Jpsi = %1.0f #pm %1.0f\n\n",fit_yield,fit_yield_err);
		}

	      // Extract background matrix
	      const int nParameter = nPar;
	      double bkg_params[nParameter];
	      double bkg_matrix[nParameter*nParameter];
	      funcBkg[i][k][t] = new TF1(Form("Jpsi_FitBkg_cent%s_%s",cent_Name[i*kNCent+k],tmpName.Data()),funcForm.Data(),g_sig_fit_min,g_sig_fit_max);
	      int nSigPar = 3;
	      if(isys==9) nSigPar = 5;
	      for(int j=0; j<nPar; j++)
		{
		  funcBkg[i][k][t]->SetParameter(j,funcSignal[i][k][t]->GetParameter(nSigPar+j));
		  funcBkg[i][k][t]->SetParError(j,funcSignal[i][k][t]->GetParError(nSigPar+j));
		  bkg_params[j] = funcSignal[i][k][t]->GetParameter(nSigPar+j);
		}

	      for(int j=nSigPar; j<nSigPar+nPar; j++)
		{
		  for(int m=nSigPar; m<nSigPar+nPar; m++)
		    {
		      bkg_matrix[(j-nSigPar)*nPar+m-nSigPar] = matrix[j*(nSigPar+nPar)+m];
		    }
		}

	      // bin counting
	      double low_mass_tmp = fit_mean - 3.5 * fit_sigma;
	      double high_mass_tmp = fit_mean + 3.5 * fit_sigma;
	      int low_bin = hSignal[i][k][t]->FindFixBin(low_mass_tmp+1e-4) + 1;
	      int high_bin = hSignal[i][k][t]->FindFixBin(high_mass_tmp-1e-4) - 1;
	      double low_bin_mass = hSignal[i][k][t]->GetXaxis()->GetBinLowEdge(low_bin);
	      double high_bin_mass = hSignal[i][k][t]->GetXaxis()->GetBinUpEdge(high_bin);
	      double yield_all_err;
	      double yield_all = hSignal[i][k][t]->IntegralAndError(low_bin,high_bin,yield_all_err);
	      double yield_bkg = funcBkg[i][k][t]->Integral(low_bin_mass,high_bin_mass) * 1./hSignal[i][k][t]->GetBinWidth(1);
	      double yield_bkg_err = funcBkg[i][k][t]->IntegralError(low_bin_mass,high_bin_mass,bkg_params,bkg_matrix)* 1./hSignal[i][k][t]->GetBinWidth(1);
	      double yield_sig = yield_all - yield_bkg;
	      double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err);

	      hBinCountYield[i][t]->SetBinContent(k+1,yield_sig);
	      hBinCountYield[i][t]->SetBinError(k+1,yield_sig_err);

	      // plotting
	      cFit->cd(k+1);	  
	      hSignal[i][k][t]->SetMaximum(2*hSignal[i][k][t]->GetMaximum());
	      hSignal[i][k][t]->Draw();
	      funcSignal[i][k][t]->Draw("sames");
	      funcBkg[i][k][t]->SetLineColor(4);
	      funcBkg[i][k][t]->Draw("sames");
	      text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],cent_Name[i*kNCent+k]),0.06);
	      text->Draw();
	      TLine *line = GetLine(low_bin_mass,hSignal[i][k][t]->GetMinimum()*1.5,low_bin_mass,hSignal[i][k][t]->GetMaximum()*0.3,1);
	      line->Draw();
	      TLine *line = GetLine(high_bin_mass,hSignal[i][k][t]->GetMinimum()*1.5,high_bin_mass,hSignal[i][k][t]->GetMaximum()*0.3,1);
	      line->Draw();

	      text = GetPaveText(0.16,0.3,0.63,0.88,0.045);
	      text->SetTextFont(62);
	      text->AddText("Fitting");
	      text->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",fit_yield,fit_yield_err,fit_yield/fit_yield_err));
	      text->AddText("Counting");
	      text->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",yield_sig,yield_sig_err,yield_sig/yield_sig_err));
	      text->SetTextAlign(11);
	      text->Draw();
	    }
	  text = GetPaveText(0.7,0.8,0.35,0.4,0.06);
	  text->SetTextFont(62);
	  text->AddText("UL-MIX(UL)");
	  text->SetTextColor(4);
	  cFit->cd(1);
	  text->Draw();
	  if(savePlot) 
	    cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sSig.Fit_Jpsi_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	}
      hMean[i][0]->SetMarkerStyle(21);
      hMean[i][0]->GetYaxis()->SetRangeUser(3.08-0.02*i,3.12+0.02*i);
      hMean[i][0]->Fit(funcMean[i][0],"IR0Q");
      c = draw1D(hMean[i][0],Form("Mean of Gaussian fit to Jpsi signal (p_{T} > %1.0f GeV/c);;mean",ptBins_low[i]));
      funcMean[i][0]->SetLineColor(4);
      funcMean[i][0]->SetLineStyle(2);
      funcMean[i][0]->Draw("sames");
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sSig.FitMeanVsCent_pt%s%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],pt_Name[i],gWeightName[gApplyWeight]));
	}

      hSigma[i][0]->SetMarkerStyle(21);
      hSigma[i][0]->GetYaxis()->SetRangeUser(0,0.1);
      hSigma[i][0]->Fit(funcSigma[i][0],"IR0Q");
      c = draw1D(hSigma[i][0],Form("Sigma of Gaussian fit to Jpsi signal (p_{T} > %1.0f GeV/c);;#sigma",ptBins_low[i]));
      funcSigma[i][0]->SetLineColor(4);
      funcSigma[i][0]->SetLineStyle(2);
      funcSigma[i][0]->Draw("sames");
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_SigExt/Npart.%s%sSig.FitSigmaVsCent_pt%s%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],pt_Name[i],gWeightName[gApplyWeight]));
	}
    }

  if(saveHisto)
    {
      cout << outName.Data() << endl;
      TFile *fout = 0;
      if(isys==0) fout = TFile::Open(outName,"update");
      else        fout = TFile::Open(outNameSys,"update");
      for(int i=0; i<nPtBins; i++)
	{
	  for(int t=0; t<1; t++)
	    {
	      hFitYield[i][t]->Write("",TObject::kOverwrite);
	      hBinCountYield[i][t]->Write("",TObject::kOverwrite);
	      hMean[i][t]->Write("",TObject::kOverwrite);
	      funcMean[i][t]->Write("",TObject::kOverwrite);
	      hSigma[i][t]->Write("",TObject::kOverwrite);
	      funcSigma[i][t]->Write("",TObject::kOverwrite);
	      hChi2[i][t]->Write("",TObject::kOverwrite);
	    }
	  if(hJpsiWeight[i]) hJpsiWeight[i]->Write("",TObject::kOverwrite);
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

