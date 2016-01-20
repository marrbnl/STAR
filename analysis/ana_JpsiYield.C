TFile *f;
const Bool_t iPico = 1;
const int year = YEAR;
const int mix_type = 0; // 0 - Shuai; 1 - Rongrong
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
      run_type = "Run13_pp500";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";

      if(iPico)
	f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else
	f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  //makeHistos();
  //makeYield();
  makeYieldRun13();
  //prod_makeYield();
  //anaYield();
}

//================================================
void anaYield(int savePlot = 1)
{
  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  fileName.ReplaceAll(".root",".yield.root");

  TFile *fin = TFile::Open(fileName,"read");

  bool drawLegend = true;
  if(nCentBins==1) drawLegend = false;

  TList *list = new TList;
  TString legName[nCentBins];
  TH1F *hFitYield[nCentBins];
  TH1F *hBinCountYield[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hFitYield[i] = (TH1F*)fin->Get(Form("Jpsi_FitYield_cent%s",cent_Title[i]));
      hBinCountYield[i] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[i]));
      double count_1, error_1, count_2, error_2;
      count_1 = hFitYield[i]->IntegralAndError(1,6,error_1);
      count_2 = hBinCountYield[i]->IntegralAndError(1,6,error_2);
      printf("Total # of J/psi in %s%%: fit - %1.0f (%2.2f), bin count - %1.0f (%2.2f)\n",cent_Name[i],count_1,count_1/error_1,count_2,count_2/error_2);
      list->Add(hFitYield[i]);
      legName[i] = Form("%s%%",cent_Name[i]);
    }
  if(year==2013) c = drawHistos(list,"Jpsi_FitYield","J/psi yield from bin counting;p_{T} (GeV/c);Counts",kFALSE,0,30,kTRUE,1e-6,1.2*hBinCountYield[0]->GetMaximum(),kFALSE,drawLegend,legName,drawLegend,"Centrality",0.45,0.75,0.55,0.8,kTRUE);
  else c = drawHistos(list,"Jpsi_FitYield","J/psi yield from bin counting;p_{T} (GeV/c);Counts",kFALSE,0,30,kTRUE,1e-6,1600,kFALSE,kTRUE,legName,kTRUE,"Centrality",0.45,0.75,0.55,0.8,kTRUE);
  
  for(int i=0; i<nCentBins; i++)
    {
      hBinCountYield[i]->SetMarkerStyle(24);
      hBinCountYield[i]->SetMarkerColor(color[i]);
      hBinCountYield[i]->SetLineColor(color[i]);
      TGraphErrors *gr = new TGraphErrors(hBinCountYield[i]);
      offset_x(gr,0.2);
      gr->Draw("samesPEZ");
    }
  leg = new TLegend(0.65,0.65,0.85,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hFitYield[0],"Fitting","P");
  leg->AddEntry(hBinCountYield[0],"Bin counting","P");
  leg->Draw();
  

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Yield_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Yield_CentBins.png",run_type,run_cfg_name.Data()));
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

  const char *name[2] = {"Fitting","BinCounting"};
  for(int j=0; j<2; j++)
    {
      list->Clear();
      for(int i=0; i<nCentBins; i++)
	{
	  list->Add(hSignif[i][j]);
	}
      c = drawHistos(list,Form("%s_Signif",name[j]),Form("%s: J/psi significance;p_{T} (GeV/c)",name[j]),kFALSE,0,30,kTRUE,1e-6,15,kFALSE,drawLegend,legName,drawLegend,"Centrality",0.15,0.3,0.6,0.85,kTRUE);
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Significance_%s.pdf",run_type,run_cfg_name.Data(),name[j]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Significance_%s.png",run_type,run_cfg_name.Data(),name[j]));
	}
    }
}

//================================================
void prod_makeYield()
{
  for(int i=0; i<4; i++)
    for(int j=0; j<7; j++)
      makeYield(i,j,1,1);
}

//===============================================
void makeYield(const int icent = 0, const int isys = 0, int savePlot = 0, int saveHisto = 0)
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
  TH1F *hMixUL[nPtBins];
  TH1F *hMixLS[nPtBins];
  TH1F *hMixBkg[nPtBins];

  for(Int_t i=0; i<nPtBins; i++)
    {
      hSeUL[i] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hSeLS[i] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixUL[i] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixLS[i] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
    }
  
  int nxpad, nypad;
  if(nPtBins-1==6 || nPtBins-1==5) { nxpad = 3; nypad = 2; }
  if(nPtBins-1<=4) { nxpad = 2; nypad = 2; }

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // mixed event scaling
  TH1F *hMixScale[nPtBins];
  TF1 *funcScale[nPtBins];
  TH1F *hFitScaleFactor = new TH1F(Form("hFitScaleFactor_%s%s",cent_Title[icent],sys_name[isys]),Form("hFitScaleFactor_%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TH1F *hBinCountScaleFactor = new TH1F(Form("hBinCountScaleFactor_%s%s",cent_Title[icent],sys_name[isys]),Form("hBinCountScaleFactor_%s%s",cent_Title[icent],sys_name[isys]),nbins,xbins);
  TCanvas *cScaling = new TCanvas(Form("mix_scale_%s",cent_Name[icent]),Form("mix_scale_%s",cent_Name[icent]),1100,650);
  cScaling->Divide(nxpad,nypad);
  double scale_low = 2.7;
  double scale_high = 3.5;
  if(isys==1)
    {
      scale_low = 2.5; scale_high = 3.7;
    }
  if(isys==2)
    {
      scale_low = 2.9; scale_high = 3.3;
    }
  for(int i=0; i<nPtBins; i++)
    {
      hMixScale[i] = (TH1F*)hSeLS[i]->Clone(Form("LS_ratio_SE_to_ME_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone",hMixLS[i]->GetName()));
      hMixScale[i]->Rebin(4);
      htmp->Rebin(4);
      hMixScale[i]->Divide(htmp);

      // fitting method
      funcScale[i] = new TF1(Form("Fit_%s",hMixScale[i]->GetName()),"pol0",scale_low,scale_high);
      hMixScale[i]->Fit(funcScale[i],"IR0Q");
      hFitScaleFactor->SetBinContent(i,funcScale[i]->GetParameter(0));
      hFitScaleFactor->SetBinError(i,funcScale[i]->GetParError(0));

      // bin counting method
      double se = 0, se_err = 0, me = 0, me_err = 0;
      int low_bin = hSeLS[i]->FindFixBin(scale_low+1e-4);
      int high_bin = hSeLS[i]->FindFixBin(scale_high-1e-4);
      for(int bin=low_bin; bin<=high_bin; bin++)
	{
	  se += hSeLS[i]->GetBinContent(bin);
	  se_err += TMath::Power(hSeLS[i]->GetBinError(bin),2);
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
	  hMixScale[i]->GetXaxis()->SetRangeUser(2.5,4);
	  hMixScale[i]->SetMaximum(1.2*hMixScale[i]->GetMaximum());
	  hMixScale[i]->Draw();
	  funcScale[i]->SetLineColor(2);
	  funcScale[i]->Draw("sames");
	  TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
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
  TList *list = new TList;
  TString legName[2] = {"Fitting","Bin Counting"};
  list->Add(hFitScaleFactor);
  list->Add(hBinCountScaleFactor);
  c = drawHistos(list,"Mix_scaleFactor",Form("Scale factor for mixed events (%s%%)",cent_Name[icent]),kFALSE,0,30,kTRUE,hFitScaleFactor->GetMinimum()*0.95,hFitScaleFactor->GetMaximum()*1.08,kFALSE,kTRUE,legName,kTRUE,"Scale factor",0.5,0.7,0.68,0.85,kTRUE);

  if(savePlot) 
    {
      cScaling->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_LS_SE_to_ME_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cScaling->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_LS_SE_to_ME_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sScaleFactor_SE_to_ME_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sScaleFactor_SE_to_ME_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  // Acceptance
  TCanvas *cAcc = new TCanvas(Form("Acceptance_%s",cent_Name[icent]),Form("Acceptance_%s",cent_Name[icent]),1100,650);
  cAcc->Divide(nxpad,nypad);
  TH1F *hAcc[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hAcc[i] = (TH1F*)hMixUL[i]->Clone(Form("Mix_Acceptance_cent%s_pt%s",cent_Title[icent],pt_Name[i]));
      TH1F *htmp = (TH1F*)hMixLS[i]->Clone(Form("%s_clone2",hMixLS[i]->GetName()));
      hAcc[i]->Rebin(10);
      htmp->Rebin(10);
      hAcc[i]->Divide(htmp);
      if(i>0)
	{
	  cAcc->cd(i);
	  hAcc[i]->SetTitle("");
	  hAcc[i]->GetXaxis()->SetRangeUser(2,4);
	  hAcc[i]->GetYaxis()->SetRangeUser(0.9,1.1);
	  hAcc[i]->Draw("HIST");
	  TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
	  t->Draw();
	  TLine *line = GetLine(2,1,4,1,1);
	  line->Draw();
	}
    }
  t = GetPaveText(0.2,0.35,0.7,0.8,0.06);
  t->SetTextFont(62);
  t->AddText("Mix: UL/LS");
  t->SetTextColor(4);
  cAcc->cd(1);
  t->Draw();

  if(savePlot) 
    {
      cAcc->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sMix_Acceptance_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cAcc->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sMix_Acceptance_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    } 

  // jpsi signal
  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",cent_Name[icent]),Form("InvMass_%s",cent_Name[icent]),1100,650);
  cSignal->Divide(nxpad,nypad);
  int rebin_signal = 10;
  if(isys==6)
    {
      rebin_signal = 5;
    }

  for(int i=0; i<nPtBins; i++)
    {
      hMixBkg[i] = (TH1F*)hMixUL[i]->Clone(Form("mix_bkg_pt%s_cent%s",pt_Name[i],cent_Name[icent]));
      hMixBkg[i]->Scale(hBinCountScaleFactor->GetBinContent(i));
      hSeUL[i]->Rebin(rebin_signal);
      hSeUL[i]->SetTitle("");
      hSeUL[i]->SetMarkerStyle(21);
      hSeUL[i]->SetMarkerColor(2);
      hSeUL[i]->SetLineColor(2);
      hSeUL[i]->GetXaxis()->SetRangeUser(2.5,4);
      hSeLS[i]->Rebin(rebin_signal);
      hMixBkg[i]->SetLineColor(4);
      hMixBkg[i]->Rebin(rebin_signal);
      if(i>0)
	{
	  cSignal->cd(i);
	  hSeUL[i]->Draw();
	  hSeLS[i]->Draw("sames HIST");
	  hMixBkg[i]->Draw("sames HIST");
	  TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
	  t->Draw();
	}
    }
  cSignal->cd(1);
  leg = new TLegend(0.45,0.63,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hSeUL[1],"Unlike sign","P");
  leg->AddEntry(hSeLS[1],"Like sign (++)+(--)","L");
  leg->AddEntry(hMixBkg[1],"Mix UL","L");
  leg->Draw();
  if(savePlot) 
    {
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sJpsiSignal_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sJpsiSignal_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  // Fit residual
  double fix_mean[nPtBins] = {0};
  double fix_sigma[nPtBins];
  if(icent>0)
    {
      TFile *fFit = 0;
      if(isys==0)  fFit = TFile::Open(outName,"read");
      else         fFit = TFile::Open(outNameSys,"read");
      TH1F *h1 = (TH1F*)fFit->Get(Form("Jpsi_FitMean_cent%s%s",cent_Title[0],sys_name[isys]));
      TH1F *h2 = (TH1F*)fFit->Get(Form("Jpsi_FitSigma_cent%s%s",cent_Title[0],sys_name[isys]));
      for(int bin=0; bin<=h1->GetNbinsX(); bin++)
	{
	  fix_mean[bin] = h1->GetBinContent(bin);
	  fix_sigma[bin] = h2->GetBinContent(bin);
	}
      fFit->Close();
    }

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
  TString funcForm2 = "pol1";
  int nPar2 = 2;
  double fit_min = 2.6, fit_max = 3.8;
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
      for(int bin=1; bin<=hSignal[i]->GetNbinsX(); bin++)
	{
	  if(hSignal[i]->GetBinContent(bin)==0)
	    {
	      hSignal[i]->SetBinContent(bin,0);
	      hSignal[i]->SetBinError(bin,1.4);
	    }
	}
      hSignal[i]->Add(hMixBkg[i],-1);
      //hSignal[i]->Add(hSeLS[i],-1);
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
      funcSignal[i] = new TF1(Form("Fit%s",hSignal[i]->GetName()),Form("gausn(0)+%s(3)",funcForm.Data()),fit_min,fit_max);
      funcSignal[i]->SetParameter(1,3.09);
      funcSignal[i]->SetParameter(2,0.06);
      //funcSignal[i]->SetParameter(0,300);
      //funcSignal[i]->SetParLimits(1,3.08,3.1);
      //funcSignal[i]->SetParLimits(2,0.04,0.07);
      if(icent>0 && fix_mean[0]>0)
	{
	  funcSignal[i]->FixParameter(1,fix_mean[i]);
	  funcSignal[i]->FixParameter(2,fix_sigma[i]);
	}
      funcSignal[i]->SetLineColor(2);
      TFitResultPtr ptr = hSignal[i]->Fit(funcSignal[i],"IRS0Q");
      printf("Fit signal = %1.1f +/- %1.1f (%1.1fsigma)\n",funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParameter(0)/funcSignal[i]->GetParError(0));
      hFitYield->SetBinContent(i,funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1));
      hFitYield->SetBinError(i,funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1));
      hMean->SetBinContent(i,funcSignal[i]->GetParameter(1));
      hMean->SetBinError(i,funcSignal[i]->GetParError(1));
      hSigma->SetBinContent(i,funcSignal[i]->GetParameter(2));
      hSigma->SetBinError(i,funcSignal[i]->GetParError(2));
      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();
      if(0)
	{
	  ptr->GetCovarianceMatrix().Print();
	  for(int im=0; im<(3+nPar)*(3+nPar); im++) cout << matrix[im] << "  ";
	  cout << endl;
	  double *parameters = ptr->GetParams();
	  for(int im=0; im<(3+nPar); im++) cout << parameters[im] << "  ";
	  cout << endl;
	}

      // bin counting
      double low_mass_tmp = low_mass;
      double high_mass_tmp = high_mass;
      if(i<4) 
	{
	  low_mass_tmp = 2.95;
	  high_mass_tmp = 3.25;
	}

      int low_bin = hSignal[i]->FindFixBin(low_mass_tmp+1e-4);
      int high_bin =hSignal[i]->FindFixBin(high_mass_tmp-1e-4);
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
      hBinCountYield->SetBinContent(i,signal);
      hBinCountYield->SetBinError(i,sig_err);
      double all = hSeUL[i]->Integral(hSeUL[i]->FindBin(low_mass_tmp+1e-4),hSeUL[i]->FindBin(high_mass_tmp-1e-4));
      hSigToBkg->SetBinContent(i,signal/(all-signal));
      printf("all = %1.2f, signal = %1.2f\n",all,signal);

      // plotting
      if(i==0) cFit1->cd();
      else     cFit->cd(i);	  
      hSignal[i]->SetMaximum(2*hSignal[i]->GetMaximum());
      if(i<4)  hSignal[i]->SetMaximum(1.2*hSignal[i]->GetMaximum());
      hSignal[i]->Draw();
      funcSignal[i]->Draw("sames");
      functmp->Draw("sames");
      TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[icent]),0.06);
      t->Draw();
      TLine *line = GetLine(low_mass_tmp,hSignal[i]->GetMinimum()*1.5,low_mass_tmp,hSignal[i]->GetMaximum()*0.3,1);
      line->Draw();
      TLine *line = GetLine(high_mass_tmp,hSignal[i]->GetMinimum()*1.5,high_mass_tmp,hSignal[i]->GetMaximum()*0.3,1);
      line->Draw();

      t = GetPaveText(0.16,0.3,0.63,0.88,0.045);
      t->SetTextFont(62);
      t->AddText("Fitting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",funcSignal[i]->GetParameter(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParError(0)/hSignal[i]->GetBinWidth(1),funcSignal[i]->GetParameter(0)/funcSignal[i]->GetParError(0)));
      t->AddText("Counting");
      t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",signal,sig_err,signal/sig_err));
      t->SetTextAlign(11);
      t->Draw();
    }
  t = GetPaveText(0.7,0.8,0.16,0.22,0.06);
  t->SetTextFont(62);
  t->AddText("LS-MIX(LS)");
  t->SetTextColor(4);
  cFit->cd(1);
  t->Draw();

  if(savePlot) 
    {
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_All_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_All_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  hMean->SetMarkerStyle(21);
  hMean->GetYaxis()->SetRangeUser(3.06,3.12);
  c = draw1D(hMean,Form("Mean of Gaussian fit to Jpsi signal (%s%%);p_{T} (GeV/c);mean",cent_Name[icent]));
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

    }

  hSigma->SetMarkerStyle(21);
  hSigma->GetYaxis()->SetRangeUser(0,0.12);
  c = draw1D(hSigma,Form("Sigma of Gaussian fit to Jpsi signal (%s%%);p_{T} (GeV/c);#sigma",cent_Name[icent]));

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
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
      hFitScaleFactor->Write("",TObject::kOverwrite);
      hBinCountScaleFactor->Write("",TObject::kOverwrite);
      hMean->Write("",TObject::kOverwrite);
      hSigma->Write("",TObject::kOverwrite);
      hFitYield->Write("",TObject::kOverwrite);
      hBinCountYield->Write("",TObject::kOverwrite);
      hSigToBkg->Write("",TObject::kOverwrite);
      if(isys==0)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      hAcc[i]->Write("",TObject::kOverwrite);
	      hSignal[i]->Write("",TObject::kOverwrite);
	      hSignalSave[i]->Write("",TObject::kOverwrite);
	      funcSignal[i]->Write("",TObject::kOverwrite);
	      hMixBkg[i]->Write("",TObject::kOverwrite);
	      hSeUL[i]->Write("",TObject::kOverwrite);
	      hSeLS[i]->Write("",TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

//================================================
void makeYieldRun13(const int icent = 0, const int isys = 0, int savePlot = 1, int saveHisto = 1)
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
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sJpsiSignal_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sJpsiSignal_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
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
      functmp->FixParameter(2,0.06);
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
      TFitResultPtr ptr = hSignal[i]->Fit(func_tmp,"IR0QS");
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
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_All_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_All_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
    }

  hMean->SetMarkerStyle(21);
  hMean->GetYaxis()->SetRangeUser(3.0,3.12);
  c = draw1D(hMean,Form("Mean of Gaussian fit to Jpsi signal;p_{T} (GeV/c);mean"));
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Mean_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));

    }

  hSigma->SetMarkerStyle(21);
  hSigma->GetYaxis()->SetRangeUser(0,0.12);
  c = draw1D(hSigma,Form("Sigma of Gaussian fit to Jpsi signal;p_{T} (GeV/c);#sigma"));

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%s%sFit_Jpsi_Sigma_cent%s.png",run_type,run_cfg_name.Data(),sys_title[isys],cent_Title[icent]));
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
void makeHistos()
{
  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[3];
  TH1F *hInvMass[nCentBins][nPtBins][3];

  // same event
  for(Int_t j=0; j<3; j++)
    {
      hnInvMass[j] = (THnSparseF*)f->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
      hnInvMass[j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

      for(Int_t i=0; i<nPtBins; i++)
	{
	  hnInvMass[j]->GetAxis(1)->SetRangeUser(ptBins_low[i]+0.01,ptBins_high[i]-0.01);
	  for(int k=0; k<nCentBins; k++)
	    {
	      hnInvMass[j]->GetAxis(7)->SetRange(centBins_low[k],centBins_high[k]);
	      hInvMass[k][i][j] = (TH1F*)hnInvMass[j]->Projection(0);
	      hInvMass[k][i][j]->SetName(Form("%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",hName[j],trigName[kTrigType],i,k));
	      hInvMass[k][i][j]->Sumw2();
	      hnInvMass[j]->GetAxis(7)->SetRange(0,-1);
	    }
	  hnInvMass[j]->GetAxis(1)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
    }

  for(int k=0; k<nCentBins; k++)
    {
      for(Int_t i=0; i<nPtBins; i++)
	{
	  hInvMass[k][i][1]->Add(hInvMass[k][i][2]);
	}
    }

  // mixed event
  TFile *fmix = 0;
  if(iPico)
    {
      if(year==2014) 
	{
	  char *mixName = Form("Pico.Run14.AuAu200.jpsi.mix.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
	  fmix = TFile::Open(Form("Output/%s",mixName),"read");
	}
    }
  if(fmix) 
    {
      cout << "Mix file: " << fmix->GetName() << endl;
      TH1F *hMixInvMass[nCentBins][nPtBins][3];
      if(mix_type==0)
	{
	  printf("INFO: using Shuai's mixed events\n");
	  TH3D *hMixMmumuvsPtCen[3];
	  hMixMmumuvsPtCen[0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
	  hMixMmumuvsPtCen[1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
	  hMixMmumuvsPtCen[2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");
	  for(Int_t j=0; j<3; j++)
	    {
	      for(int i=0; i<nPtBins; i++)
		{
		  int ybin_min = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBins_low[i]+1e-4);
		  int ybin_max = hMixMmumuvsPtCen[j]->GetYaxis()->FindFixBin(ptBins_high[i]-1e-4);
		  for(int k=0; k<nCentBins; k++)
		    {
		      TH1F *htmp = (TH1F*)hMixMmumuvsPtCen[j]->ProjectionZ(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",hName[j],trigName[kTrigType],i,k),centBins_low[k],centBins_high[k],ybin_min,ybin_max);
		      hMixInvMass[k][i][j] = new TH1F(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),2800,0,14);
		      for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
			{
			  hMixInvMass[k][i][j]->SetBinContent(bin,htmp->GetBinContent(bin));
			  hMixInvMass[k][i][j]->SetBinError(bin,htmp->GetBinError(bin));
			}
		    }
		}
	    }
	}
      else if(type==1)
	{
	  printf("INFO: using Rongrong's mixed events\n");
	  THnSparseF *hnMixInvMass[3];
	  for(Int_t j=0; j<3; j++)
	    {
	      hnMixInvMass[j] = (THnSparseF*)fmix->Get(Form("%s_%s",hName[j],trigName[kTrigType]));
	      hnMixInvMass[j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
	      hnMixInvMass[j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
	  
	      for(Int_t i=0; i<nPtBins; i++)
		{
		  hnMixInvMass[j]->GetAxis(1)->SetRangeUser(ptBins_low[i]+0.01,ptBins_high[i]-0.01);
		  for(int k=0; k<nCentBins; k++)
		    {
		      hnMixInvMass[j]->GetAxis(4)->SetRange(centBins_low[k],centBins_high[k]);
		      hMixInvMass[k][i][j] = (TH1F*)hnMixInvMass[j]->Projection(0);
		      hMixInvMass[k][i][j]->SetName(Form("mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",hName[j],trigName[kTrigType],i,k));
		      hMixInvMass[k][i][j]->Sumw2();
		      hnMixInvMass[j]->GetAxis(4)->SetRange(0,-1);
		    }
		  hnMixInvMass[j]->GetAxis(1)->SetRange(0,-1);
		}
	      hnMixInvMass[j]->GetAxis(2)->SetRange(0,-1);
	      hnMixInvMass[j]->GetAxis(3)->SetRange(0,-1);
	    }
	}

      for(int k=0; k<nCentBins; k++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hMixInvMass[k][i][1]->Add(hMixInvMass[k][i][2]);
	    }
	}
    }

  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  TFile *fout = TFile::Open(fileName,"recreate");
  for(int k=0; k<nCentBins; k++)
    {
      for(Int_t i=0; i<nPtBins; i++)
	{
	  hInvMass[k][i][0]->Write(Form("InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[k]));
	  hInvMass[k][i][1]->Write(Form("InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[k]));
	  if(fmix)
	    {
	      hMixInvMass[k][i][0]->Write(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],cent_Name[k]));
	      hMixInvMass[k][i][1]->Write(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],cent_Name[k]));
	    }
	}
    }
  fout->Close();
}

//================================================
double Polynomial3(double *x, double *par)
{
  if(x[0]>low_mass && x[0]<high_mass)
    {
      TF1::RejectPoint();
      return 0;
    }
  
  return par[0]+par[1]*x[0]+par[2]*x[0]**2+par[3]*x[0]**3;
}

//================================================
double Polynomial1(double *x, double *par)
{
  if(x[0]>low_mass && x[0]<high_mass)
    {
      TF1::RejectPoint();
      return 0;
    }
  
  return par[0]+par[1]*x[0];
}

//================================================
void checkStatError(TH1F *hSignal, TF1 *funcSignal, const int N = 1e3)
{
  TRandom3 *rndm = new TRandom3();

  TF1 *func = new TF1("func",Polynomial1,2.7,4,2);
  TH1F *h = (TH1F*)hSignal->Clone("signal");
  TFitResultPtr ptr = h->Fit(func,"IR0QS");
  draw1D(h);
  //func->Draw("sames");
  // TF1 *sig = new TF1("sig","gausn");
  // for(int j=0; j<3; j++)
  //   {
  //     sig->SetParameter(j,func->GetParameter(j));
  //     sig->SetParError(j,func->GetParError(j));
  //   }
  // double bkg = (func->Integral(low_mass,high_mass)-sig->GetParameter(0)) * 1./hSignal->GetBinWidth(1);
  TF1 *bkgfunc = new TF1("bkgfunc","pol1",2.7,4);
  for(int j=0; j<2; j++)
    {
      bkgfunc->SetParameter(j,func->GetParameter(j));
      bkgfunc->SetParError(j,func->GetParError(j));
    }
  double bkg = bkgfunc->Integral(low_mass,high_mass) * 1./hSignal->GetBinWidth(1);
  double bkg_err = TMath::Sqrt(TMath::Abs(bkg));
  double bkg_err_2 = bkgfunc->IntegralError(low_mass,high_mass,ptr.Get()->GetParams(), ptr->GetCovarianceMatrix()->GetMatrixArray())* 1./hSignal->GetBinWidth(1);
  bkgfunc->SetLineColor(2);
  bkgfunc->SetLineStyle(2);
  bkgfunc->Draw("sames");
  cout << bkg << "  " << bkg_err_2 << endl;
  return;
  TH1F *hMean = new TH1F("hMean","Distribution background counts;N",50,-1000,1200);
  int low_bin = h->GetXaxis()->FindFixBin(2.5);
  int high_bin = h->GetXaxis()->FindFixBin(4.5);
  for(int i=0; i<N; i++)
    {
      TH1F *htmp = (TH1F*)h->Clone(Form("htmp_%d",i));
      for(int bin=low_bin; bin<=high_bin; bin++)
	{
	  htmp->SetBinContent(bin,rndm->Gaus(h->GetBinContent(bin),h->GetBinError(bin)));
	  htmp->SetBinError(bin,h->GetBinError(bin));
	}
      htmp->Fit(func,"IRQ0");
      //draw1D(htmp);
      //func->Draw("sames");
      for(int j=0; j<2; j++)
	{
	  bkgfunc->SetParameter(j,func->GetParameter(j));
	  bkgfunc->SetParError(j,func->GetParError(j));
	}
      double new_bkg = bkgfunc->Integral(low_mass,high_mass) * 1./hSignal->GetBinWidth(1);
      //cout << new_bkg << endl;
      hMean->Fill(new_bkg);
    }
  TF1 *func2 = new TF1("fun2","gaus",-500,1200);
  fun2->SetParameter(1,bkg);
  fun2->SetParameter(2,bkg_err_2);
  hMean->Fit(func2,"IR0");
  draw1D(hMean);
  func2->Draw("sames");
}


      /*
      // fit background only
      TF1 *func_tmp = 0;
      if(funcForm=="pol1") func_tmp = new TF1(Form("func_tmp_%d",i),Polynomial1,fit_min,fit_max,nPar);
      if(funcForm=="pol3") func_tmp = new TF1(Form("func_tmp_%d",i),Polynomial3,fit_min,fit_max,nPar);
      TFitResultPtr ptr = hSignal[i]->Fit(func_tmp,"IR0QSN");
      TF1 *bkgfunc_tmp = new TF1(Form("bkgfunc_tmp_%d",i),Form("%s",funcForm.Data()),fit_min,fit_max);
      for(int j=0; j<nPar; j++)
      	{
      	  bkgfunc_tmp->SetParameter(j,func_tmp->GetParameter(j));
      	  bkgfunc_tmp->SetParError(j,func_tmp->GetParError(j));
      	}
      double bkg_2 = bkgfunc_tmp->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal[i]->GetBinWidth(1);
      double bkg_err_2 = bkgfunc_tmp->IntegralError(low_mass_tmp,high_mass_tmp,ptr->GetParams(), ptr->GetCovarianceMatrix().GetMatrixArray())* 1./hSignal[i]->GetBinWidth(1);
      bkgfunc_tmp->SetLineColor(4);
      bkgfunc_tmp->SetLineStyle(2);
      bkgfunc_tmp->Draw("sames");
      double signal_2 = count-bkg_2;
      double sig_err_2 = sqrt(error*error+bkg_err_2*bkg_err_2);
      printf("Fit background = %1.1f +/- %1.1f -> %1.1f +/- %1.1f (%1.1fsigma)\n",bkg_2,bkg_err_2,signal_2,sig_err_2,signal_2/sig_err_2);

      // Side band
      int side_low_bin = low_bin - (high_bin-low_bin + 1);
      int side_high_bin = high_bin + (high_bin-low_bin + 1);
      double bkg_3 = 0, bkg_err_3 = 0;
      for(int bin=side_low_bin; bin<=side_high_bin; bin++)
	{
	  if(bin>=low_bin && bin<=high_bin) continue;
	  bkg_3 += hSignal[i]->GetBinContent(bin);
	  bkg_err_3 += TMath::Power(hSignal[i]->GetBinError(bin),2);
	}
      bkg_err_3 = sqrt(bkg_err_3);
      //bkg_3 /= 2.;
      //bkg_err_3 /= 2.;
      double signal_3 = count-bkg_3;
      double sig_err_3 = sqrt(error*error+bkg_err_3*bkg_err_3);
      printf("Side-band = %1.1f +/- %1.1f -> %1.1f +/- %1.1f (%1.1fsigma)\n",bkg_3,bkg_err_3,signal_3,sig_err_3,signal_3/sig_err_3);
      */

