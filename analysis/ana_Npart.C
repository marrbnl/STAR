TFile *f;
const int year = YEAR;
TString run_cfg_name;
TString outFileName;
const int gApplyWeight = 1;
const int nPtBins = 2;
const double ptBins_low[nPtBins]  = {0,5};
const double ptBins_high[nPtBins] = {15,15};
const char *pt_Name[nPtBins] = {"0-15","5-15"};
const int nCentBins = 6;
const int CentBins_low[nCentBins]  = {5,7,9, 11,13,15};
const int CentBins_high[nCentBins] = {6,8,10,12,14,16};
const char *Cent_Name[nCentBins] = {"50-60","40-50","30-40","20-30","10-20","0-10"};

//================================================
void ana_Npart()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2014)
    {
      f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      outFileName = Form("Rootfiles/%s.Npart.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut);
    }
  run_cfg_name = Form("%s",run_config);

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("\n+++++ Event statistics +++++\n");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  //makeDataHistos();
  //makeEmbed();
  //makeJpsiYield();
  //xsec();
  //prodYieldSys();
}

//================================================
void prodYieldSys()
{
  for(int i=0; i<7; i++)
    {
      makeJpsiYield(i,1,1);
    }
}

//================================================
void xsec(const bool savePlot = 1, const bool saveHisto = 1)
{
  // Get correction factor for othe centralities
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(outFileName.Data(),"update");
  else          fin = TFile::Open(outFileName.Data(),"read");

  // =============================================
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  double mb_events_tmp[4];
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEvtRun = (TH1F*)f->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)f->Get("mhEvtRunAcc_di_mu");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents[4];
  for(int i=0; i<4; i++)
    {
      hEqMbEvents[i] = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
      mb_events_tmp[i] = 0;
      for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
	{
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	  double run = hEvtRunAcc->GetBinCenter(bin);
	  double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(run));
	  if(nEventsTaken==0) 
	    {
	      if(i==0) printf("[w] check run %1.0f\n",run);
	      continue;
	    }
	  double nEventsRun = hEvtRun->GetBinContent(bin);
	  double rf = hRF->GetBinContent(hRF->FindFixBin(run));
	  if(rf==0)
	    {
	      printf("[w] rf = 0 for run %1.0f\n",run);
	      rf = 0.49;
	    }
	  double eq_mb = hEqMbEvents[i]->GetBinContent(hEqMbEvents[i]->FindFixBin(run));
	  mb_events_tmp[i] += nEventsRun/rf/nEventsTaken * eq_mb;

	}
      if(i==0) printf("Effective # of MB events for %s%%: %4.4e\n",cent_Name[i],mb_events_tmp[i]);
    }
  double mb_events[nCentBins] = {mb_events_tmp[3]/2,mb_events_tmp[3]/2,mb_events_tmp[2]/2,mb_events_tmp[2]/2,mb_events_tmp[1]/2,mb_events_tmp[1]/2};
  for(int i=0; i<nCentBins; i++)
    {
      printf("Effective # of MB events for %s%%: %4.4e\n",Cent_Name[i],mb_events[i]);
    }
  
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================
  //
  //

  // =============================================
  // MTD acceptance loss
  double evtCount[4] = {0,0,0,0};
  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
    {
      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
      double run = hEvtRunAcc->GetBinCenter(bin);
      double evt = hEvtRunAcc->GetBinContent(bin);
      if(run<15098067) evtCount[0] += evt;
      else if(run<15106131) evtCount[1] += evt;
      else evtCount[2] += evt;
      evtCount[3] += evt;
    }
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"read");
  TH1F *hAccLoss[3];
  for(int i=0; i<3; i++)
    {
      hAccLoss[i] = (TH1F*)fAcc->Get(Form("hAccepLoss_%d",i));
    }
  TH1F *hAccCorr = (TH1F*)hAccLoss[0]->Clone("hAccCorr");
  hAccCorr->Reset();
  for(int i=0; i<3; i++)
    {
      hAccCorr->Add(hAccLoss[i],evtCount[i]/evtCount[3]);
    }

  TFile *fWeight = TFile::Open("Rootfiles/Run14_AuAu200.Input.root","read");
  TH1F *hMcJpsiPt = (TH1F*)fWeight->Get(Form("hInputJpsiShape_Cent0"));
  double mtd_acc[nPtBins] = {0,0};
  for(int i=0; i<nPtBins; i++)
    {
      int low_bin = hAccCorr->FindFixBin(ptBins_low[i]);
      for(int bin=low_bin; bin<=hAccCorr->GetNbinsX(); bin++)
	{
	  double bin1 = hMcJpsiPt->FindFixBin(hAccCorr->GetXaxis()->GetBinLowEdge(bin));
	  double bin2 = hMcJpsiPt->FindFixBin(hAccCorr->GetXaxis()->GetBinUpEdge(bin));
	  mtd_acc[i] += hAccCorr->GetBinContent(bin) * hMcJpsiPt->Integral(bin1,bin2);
	}
      mtd_acc[i] /= hMcJpsiPt->Integral(hMcJpsiPt->FindFixBin(ptBins_low[i]),hMcJpsiPt->FindFixBin(ptBins_high[i]));
      printf("[i] MTD acceptance correction is %4.2f for pt > %1.0f\n",mtd_acc[i],ptBins_low[i]);
    }
  
  // =============================================
  //
  //
  TList *list = new TList;
  TString legName[nPtBins];

  // Jpsi raw counts &  efficiency
  TH1F *hJpsiCounts[nPtBins];
  TH1F *hJpsiEff[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiEff[i] = (TH1F*)fin->Get(Form("hJpsiEffFinal_pt%s",pt_Name[i]));
      hJpsiCounts[i] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_pt%s_weight",pt_Name[i]));
      hJpsiCounts[i]->GetXaxis()->SetLabelSize(0.05);
      hJpsiCounts[i]->SetMarkerSize(1.5);
      list->Add(hJpsiCounts[i]);
      legName[i] = Form("p_{T,J/#Psi} > %1.0f GeV/c",ptBins_low[i]);
    }
  c = drawHistos(list,"hJpsiCounts","Number of J/#Psi in each centrality bin;;N",kFALSE,0,10,true,10,1e6,true,kTRUE,legName,kTRUE,"Run14_AuAu200",0.45,0.65,0.65,0.85,kTRUE);
  if(savePlot)      
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/JpsiRawCounts.pdf",run_type));
  list->Clear();

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiInvYield[i] = (TH1F*)hJpsiCounts[i]->Clone(Form("Jpsi_InvYield_pt%s",pt_Name[i]));
      hJpsiInvYield[i]->Divide(hJpsiEff[i]);
      hJpsiInvYield[i]->Scale(1./mtd_acc[i]);
      hJpsiInvYield[i]->Scale(1./mb_events[i]); // N_evt
      hJpsiInvYield[i]->Scale(1./(2*pi)); // 2pi
      hJpsiInvYield[i]->Scale(1./1.6); // dy
      hJpsiInvYield[i]->SetMarkerStyle(21);
      hJpsiInvYield[i]->SetMarkerColor(2);
      hJpsiInvYield[i]->SetLineColor(2);
      list->Add(hJpsiInvYield[i]);
    }
  c = drawHistos(list,"hJpsiXsec","J/#Psi invariant yield per event;;1/N_{MB} dN/(2#pidy)",kFALSE,0,10,true,1e-8,1e-2,true,kTRUE,legName,kTRUE,"Run14_AuAu200",0.45,0.65,0.65,0.85,kTRUE);
  if(savePlot)      
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/JpsiXsec.pdf",run_type));
  list->Clear();

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiInvYield[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void makeEmbed(const bool savePlot = 1, const bool saveHisto = 1)
{
  THnSparseF *hnJpsiInfo[2];
  TFile *femb = TFile::Open(Form("./output/Run14.AuAu200.Jpsi.Embed.%sroot",run_config),"read");
  hnJpsiInfo[0] = (THnSparseF*)femb->Get("hJpsiInfo_MC_di_mu_w");
  hnJpsiInfo[1] = (THnSparseF*)femb->Get("hJpsiInfo_MtdTrig_di_mu_w");


  TH1F *hJpsiCounts[2][nPtBins][gNTrgSetup];
  TH1F *hJpsiEff[nPtBins][gNTrgSetup];
  for(int j=0; j<2; j++)
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      hJpsiCounts[j][i][t] = new TH1F(Form("JpsiCounts_pt%s%s_%d",pt_Name[i],gTrgSetupName[t],j),"",nCentBins,0,nCentBins);
	      for(int bin=1; bin<=nCentBins; bin++)
		{
		  hJpsiCounts[j][i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
		  hJpsiCounts[j][i][t]->GetXaxis()->SetLabelSize(0.05);
		}
	    }
	}
    }
  

  for(int j=0; j<2; j++)
    {
      if(j>0)
	{
	  hnJpsiInfo[j]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	  hnJpsiInfo[j]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	}
      for(int t=0; t<gNTrgSetup; t++)
	{
	  if(t>0) hnJpsiInfo[j]->GetAxis(6)->SetRange(t,t);
	  for(int k=0; k<nCentBins; k++)
	    {
	      hnJpsiInfo[j]->GetAxis(5)->SetRange(CentBins_low[k],CentBins_high[k]);
	      TH1F *htmp = (TH1F*)hnJpsiInfo[j]->Projection(1);
	      htmp->SetName(Form("hJpsiPt_%d%d%d",j,t,k));
	      htmp->SetTitle("");
	      for(int i=0; i<nPtBins; i++)
		{
		  int low_bin = htmp->FindFixBin(ptBins_low[i]+1e-4);
		  int hi_bin  = htmp->FindFixBin(ptBins_high[i]-1e-4);
		  double err;
		  double yield = htmp->IntegralAndError(low_bin,hi_bin,err);
		  hJpsiCounts[j][i][t]->SetBinContent(k+1,yield);
		  hJpsiCounts[j][i][t]->SetBinError(k+1,err);
		}
	      hnJpsiInfo[j]->GetAxis(5)->SetRange(0,-1);
	    }
	  hnJpsiInfo[j]->GetAxis(6)->SetRange(0,-1);
	}
      hnJpsiInfo[j]->GetAxis(3)->SetRange(0,-1);
      hnJpsiInfo[j]->GetAxis(4)->SetRange(0,-1);
    }

  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  hJpsiEff[i][t] = (TH1F*)hJpsiCounts[1][i][t]->Clone(Form("JpsiEff_pt%s%s",pt_Name[i],gTrgSetupName[t]));
	  hJpsiEff[i][t]->Divide(hJpsiCounts[0][i][t]);
	}
    }

  // Get Jpsi counts as weights
  TFile *fYield = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiWeight[gNTrgSetup-1];
  double nJpsi[nCentBins-1][gNTrgSetup-1];
  double nJpsiCent[nCentBins-1];
  for(int k=0; k<nCentBins-1; k++)
    {
      nJpsiCent[k] = 0;
    }
  for(int i=0; i<gNTrgSetup-1; i++)
    {
      hJpsiWeight[i] = (TH1F*)fYield->Get(Form("NJpsiInCent_weight%s",gTrgSetupName[i+1]));
      for(int bin=1; bin<=hJpsiWeight[i]->GetNbinsX(); bin++)
	{
	  nJpsi[bin-1][i] = hJpsiWeight[i]->GetBinContent(bin);
	  nJpsiCent[bin-1] += hJpsiWeight[i]->GetBinContent(bin);
	}
    }
  for(int i=0; i<gNTrgSetup-1; i++)
    for(int bin=1; bin<=hJpsiWeight[i]->GetNbinsX(); bin++)
      printf("[i] %s%% %s: %4.2f Jpsi\n",cent_Title[bin],gTrgSetupTitle[i+1],nJpsi[bin-1][i]/nJpsiCent[bin-1]);

  // Get correction factor for othe centralities
  TList *list = new TList;
  TString legName[nPtBins];
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(outFileName.Data(),"update");
  else          fin = TFile::Open(outFileName.Data(),"read");
  TH1F *hJpsiEffCorr[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiEffCorr[i] = (TH1F*)fin->Get(Form("hJpsiEffCor_pt%s",pt_Name[i]));
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  hJpsiEffCorr[i]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	  hJpsiEffCorr[i]->GetXaxis()->SetLabelSize(0.05);
	}
      list->Add(hJpsiEffCorr[i]);
      legName[i] = Form("p_{T,J/#Psi} > %1.0f GeV/c",ptBins_low[i]);
    }
  c = drawHistos(list,"hJpsiEffCor","Ratio of J/#psi efficiency to 0-10%",kFALSE,0,10,kFALSE,0.5,1,kFALSE,kTRUE,legName,kTRUE,"Run14_AuAu200",0.45,0.65,0.65,0.85,kTRUE);
  if(savePlot)      
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/JpsiEffCorr.pdf",run_type));
  list->Clear();

  // Get average efficiency
  TH1F *hJpsiEffFinal[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiEffFinal[i] = (TH1F*)hJpsiEff[i][0]->Clone(Form("hJpsiEffFinal_pt%s",pt_Name[i]));
      hJpsiEffFinal[i]->Reset();
      // work on 0-10% first
      double avg_eff = 0, avg_err = 0;
      for(int t=0; t<gNTrgSetup-1; t++)
	{
	  double eff = hJpsiEff[i][t+1]->GetBinContent(nCentBins);
	  double err = hJpsiEff[i][t+1]->GetBinError(nCentBins);
	  avg_eff += 1/eff * nJpsi[0][t]/nJpsiCent[0];
	  avg_err += TMath::Power(nJpsi[0][t]/nJpsiCent[0]*err/eff,2);
	}
      avg_eff = 1/avg_eff;
      avg_err = avg_eff * TMath::Sqrt(avg_err);
      hJpsiEffFinal[i]->SetBinContent(nCentBins,avg_eff);
      hJpsiEffFinal[i]->SetBinError(nCentBins,avg_err);

      // corrected for other centralities
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  cout << Cent_Name[bin-1] << " -> " << hJpsiEffCorr[i]->GetBinContent(bin) << endl;
	  hJpsiEffFinal[i]->SetBinContent(bin,hJpsiEffFinal[i]->GetBinContent(nCentBins)*hJpsiEffCorr[i]->GetBinContent(bin));
	  hJpsiEffFinal[i]->SetBinError(bin,hJpsiEffFinal[i]->GetBinError(nCentBins)*hJpsiEffCorr[i]->GetBinContent(bin));
	}
      list->Add(hJpsiEffFinal[i]);
    }
  c = drawHistos(list,"hJpsiEffFinal","J/#psi efficiency in different centrality",kFALSE,0,10,true,0,0.015,kFALSE,kTRUE,legName,kTRUE,"Run14_AuAu200",0.45,0.65,0.65,0.85,kTRUE);
  if(savePlot)      
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/JpsiEffFinal.pdf",run_type));
  
  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      fin->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiEffFinal[i]->Write("",TObject::kOverwrite);
	}
    }
}

//===============================================
void makeJpsiYield(const int isys = 4, int savePlot = 0, int saveHisto = 1)
{
  const char *sys_name[7] = {"","_LargeScale","_SmallScale","_pol1","_LargeFit","_SmallFit","_Rebin"};
  const char *sys_title[7] = {"","Sys.LargeScale.","Sys.SmallScale.","Sys.pol1.","Sys.LargeFit.","Sys.SmallFit.","Sys.Rebin"};
  const TString suffix = Form("%s%s",gWeightName[gApplyWeight],sys_name[isys]);

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(outFileName.Data(),"update");
  else          fin = TFile::Open(outFileName.Data(),"read");
  printf("\n+++++ %s +++++\n",suffix.Data());
  cout << fin->GetName() << endl;

  // global setup
  double g_mix_scale_low = 2.6;
  double g_mix_scale_high = 3.6;
  double g_bin_width = 0.04; // 40 MeV
  TString g_func1 = "pol1", g_func2 = "pol3";
  int g_func1_npar = 2, g_func2_npar = 4;
  double g_sig_fit_min = 2.6;
  double g_sig_fit_max = 3.8;
  if(isys==1) { g_mix_scale_low = 2.5; g_mix_scale_high = 3.7; }
  if(isys==2) { g_mix_scale_low = 2.7; g_mix_scale_high = 3.5; }
  if(isys==3) { g_func1 = "pol0"; g_func1_npar = 1; g_func2 = "pol2"; g_func2_npar = 3;}
  if(isys==4) { g_sig_fit_min = 2.5; g_sig_fit_max = 3.9; }
  if(isys==5) { g_sig_fit_min = 2.7; g_sig_fit_max = 3.7; }
  if(isys==6) { g_bin_width = 0.02; }

  const int nSetup = 1;

  // get histograms
  TH1F *hSeUL[nPtBins][nCentBins][nSetup];
  TH1F *hSeLS[nPtBins][nCentBins][nSetup];
  TH1F *hMixUL[nPtBins][nCentBins];
  TH1F *hMixLS[nPtBins][nCentBins];
  for(Int_t i=0; i<nPtBins; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  for(int t=0; t<nSetup; t++)
	    {
	      hSeUL[i][k][t] = (TH1F*)fin->Get(Form("InvMass_UL_pt%s_cent%s%s%s",pt_Name[i],Cent_Name[k],gWeightName[gApplyWeight],gTrgSetupName[t]));
	      hSeLS[i][k][t] = (TH1F*)fin->Get(Form("InvMass_LS_pt%s_cent%s%s%s",pt_Name[i],Cent_Name[k],gWeightName[gApplyWeight],gTrgSetupName[t]));
	    }
	  hMixUL[i][k] = (TH1F*)fin->Get(Form("Mix_InvMass_UL_pt%s_cent%s",pt_Name[i],Cent_Name[k]));
	  hMixLS[i][k] = (TH1F*)fin->Get(Form("Mix_InvMass_LS_pt%s_cent%s",pt_Name[i],Cent_Name[k]));
	}
    }  
  int nxpad = 3, nypad = 2;

  // mixed event scaling
  TList *list = new TList;
  TString legName[2] = {"Fitting","Bin Counting"};
  TH1F *hMixScale[nPtBins][nCentBins][nSetup];
  TF1 *funcScale[nPtBins][nCentBins][nSetup];
  TH1F *hFitScaleFactor[nPtBins][nSetup];
  TH1F *hBinCountScaleFactor[nPtBins][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<nSetup; t++)
	{
	  TString tmpName = Form("pt%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t]);
	  hFitScaleFactor[i][t] = new TH1F(Form("FitScaleFactor_%s",tmpName.Data()),"Mixed-event scale factor from fitting",nCentBins,0,nCentBins);
	  hBinCountScaleFactor[i][t] = new TH1F(Form("BinCountScaleFactor_%s",tmpName.Data()),"Mixed-event scale factor from bin counting",nCentBins,0,nCentBins);
	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      hFitScaleFactor[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	      hBinCountScaleFactor[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	    }
	  TCanvas *cScaling = new TCanvas(Form("mix_scale_%s",tmpName.Data()),Form("mix_scale_%s",tmpName.Data()),1100,650);
	  cScaling->Divide(nxpad,nypad);

	  for(int k=0; k<nCentBins; k++)
	    {
	      hMixScale[i][k][t] = (TH1F*)hSeLS[i][k][t]->Clone(Form("%s_MixScale",hSeLS[i][k][t]->GetName()));
	      hMixScale[i][k][t]->Rebin(2);
	      TH1F *htmp = (TH1F*)hMixLS[i][k]->Clone(Form("%s_%d_clone",hMixLS[i][k]->GetName(),t));
	      htmp->Rebin(int(hMixScale[i][k][t]->GetBinWidth(1)/hMixLS[i][k]->GetBinWidth(1)));
	      hMixScale[i][k][t]->Divide(htmp);

	      // fitting method
	      funcScale[i][k][t] = new TF1(Form("Fit_%s",hMixScale[i][k][t]->GetName()),"pol0",g_mix_scale_low,g_mix_scale_high);
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
	      cScaling->cd(nCentBins-k);
	      hMixScale[i][k][t]->SetTitle("");
	      hMixScale[i][k][t]->SetMarkerStyle(21);
	      hMixScale[i][k][t]->GetXaxis()->SetRangeUser(2.5,4);
	      hMixScale[i][k][t]->SetMaximum(1.5*hMixScale[i][k][t]->GetMaximum());
	      hMixScale[i][k][t]->SetMinimum(0.5*hMixScale[i][k][t]->GetMinimum());
	      hMixScale[i][k][t]->Draw();
	      funcScale[i][k][t]->SetLineColor(2);
	      funcScale[i][k][t]->Draw("sames");
	      TPaveText *text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],Cent_Name[k]),0.05);
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
	  if(savePlot) 	 
	    cScaling->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/%s%sFitSeToMe_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));

	  list->Add(hFitScaleFactor[i][t]);
	  list->Add(hBinCountScaleFactor[i][t]);
	  c = drawHistos(list,Form("MixScale_%s",tmpName.Data()),Form("Scale factor for mixed events (p_{T} > %1.1f GeV/c)",ptBins_low[i]),kFALSE,0,30,kTRUE,0,hBinCountScaleFactor[i][t]->GetMaximum()*1.2,kFALSE,kTRUE,legName,kTRUE,Form("AuAu_200%s",gTrgSetupTitle[t]),0.5,0.7,0.68,0.85,kTRUE);
	  list->Clear();
	  if(savePlot) 
	    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/%s%sMixScaleFactor_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	}
    }

  // jpsi signal
  TH1F *hMixBkg[nPtBins][nCentBins][nSetup];
  TH1F *hSignal[nPtBins][nCentBins][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<nSetup; t++)
	{
	  TString tmpName = Form("pt%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t]);
	  TCanvas *cSignal = new TCanvas(Form("InvMass_%s",tmpName.Data()),Form("InvMass_%s",tmpName.Data()),1100,650);
	  cSignal->Divide(nxpad,nypad);

	  for(int k=0; k<nCentBins; k++)
	    {
	      hMixBkg[i][k][t] = (TH1F*)hMixUL[i][k]->Clone(Form("MixBkg_cent%s%s",Cent_Name[k],tmpName.Data()));
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
	      cSignal->cd(nCentBins-k);
	      hSeUL[i][k][t]->Draw("PE");
	      hSeLS[i][k][t]->Draw("sames HIST");
	      hMixBkg[i][k][t]->Draw("sames HIST");
	      text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],Cent_Name[k]),0.06);
	      text->Draw();
	      hSignal[i][k][t] = (TH1F*)hSeUL[i][k][t]->Clone(Form("JpsiSignal_cent%s%s",Cent_Name[k],tmpName.Data()));
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
	  leg = new TLegend(0.5,0.6,0.7,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.06);
	  leg->AddEntry(hSeUL[i][1][t],"Unlike sign","P");
	  leg->AddEntry(hSeLS[i][1][t],"Like sign (++)+(--)","L");
	  leg->AddEntry(hMixBkg[i][1][t],"Mix UL","L");
	  leg->Draw();
	  if(savePlot) 
	    cSignal->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/%s%sJpsiSignal_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	}
    }

  // Fit residual

  // extract yields from 0-60%
  TH1F *hJpsiWeight[nPtBins];
  for(int i=0; i<1; i++)
    {
      c = new TCanvas(Form("Fit_InvMass_All_pt%d",i),Form("Fit_InvMass_All_pt%d",i),1100,650);
      c->Divide(2,2);
      hJpsiWeight[i] = new TH1F(Form("hJpsiWeight_pt%s",pt_Name[i]),Form("hJpsiWeight_pt%s",pt_Name[i]),4,0,4);
      for(int t=1; t<nSetup; t++)
  	{
	  c->cd(t);
	  TH1F *htmp = (TH1F*)hSignal[i][0][t]->Clone(Form("hSignal_All_pt%s_P%d",pt_Name[i],t-1));
	  htmp->Reset();
  	  for(int k=0; k<nCentBins; k++)
  	    {
  	      htmp->Add(hSignal[i][k][t]);
  	    }
	  TF1 *functmp = new TF1(Form("FitAll_pt%s_P%d",pt_Name[i],t-1),Form("gausn(0)+pol3(3)"),g_sig_fit_min,g_sig_fit_max);
	  functmp->SetParameter(0,100);
	  functmp->SetParameter(1,3.09);
	  functmp->SetParameter(2,0.1);
	  functmp->SetLineColor(2);
	  htmp->Fit(functmp,"IRS0Q");
	  htmp->Draw("P");
	  functmp->SetLineColor(2);
	  functmp->Draw("sames");
	  double fit_yield = functmp->GetParameter(0)/htmp->GetBinWidth(1);
	  double fit_yield_err = functmp->GetParError(0)/htmp->GetBinWidth(1);
	  hJpsiWeight[i]->SetBinContent(t,fit_yield);
	  hJpsiWeight[i]->SetBinError(t,fit_yield_err);
  	}
    }  

  TH1F *hMean[nPtBins][nSetup];
  TH1F *hSigma[nPtBins][nSetup];
  TH1F *hFitYield[nPtBins][nSetup];
  TH1F *hBinCountYield[nPtBins][nSetup];
  TF1 *funcSignal[nPtBins][nCentBins][nSetup];
  TF1 *funcBkg[nPtBins][nCentBins][nSetup];
  for(int i=0; i<nPtBins; i++)
    {
      for(int t=0; t<1; t++)
	{
	  TString tmpName = Form("pt%s%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[t],sys_name[isys]);
	  hMean[i][t] = new TH1F(Form("Jpsi_FitMean_%s",tmpName.Data()),"Mean of Jpsi peak",nCentBins,0,nCentBins);
	  hSigma[i][t] = new TH1F(Form("Jpsi_FitSigma_%s",tmpName.Data()),"Sigma of Jpsi peak",nCentBins,0,nCentBins);
	  hFitYield[i][t] = new TH1F(Form("Jpsi_FitYield_%s",tmpName.Data()),"Jpsi yield from fitting",nCentBins,0,nCentBins);
	  hBinCountYield[i][t] = new TH1F(Form("Jpsi_BinCountYield_%s",tmpName.Data()),"Jpsi yield from bin counting",nCentBins,0,nCentBins);
	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      hMean[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	      hSigma[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	      hFitYield[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	      hBinCountYield[i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	    }
	  
	  TCanvas *cFit = new TCanvas(Form("Fit_Jpsi_%s",tmpName.Data()),Form("Fit_Jpsi_%s",tmpName.Data()),1100,650);
	  cFit->Divide(nxpad,nypad);

	  for(int k=nCentBins-1; k>-1; k--)
	    {
	      TString funcForm;
	      int nPar;
	      if(i==0 && k>2)
		{
		  funcForm = g_func2;
		  nPar = g_func2_npar;
		}
	      else
		{
		  funcForm = g_func1;
		  nPar = g_func1_npar;
		}

	      funcSignal[i][k][t] = new TF1(Form("Jpsi_FitSig_cent%s_%s",Cent_Name[k],tmpName.Data()),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	      funcSignal[i][k][t]->SetParameter(0,100);
	      funcSignal[i][k][t]->SetParameter(1,3.09);
	      funcSignal[i][k][t]->SetParameter(2,0.1);
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
	      hMean[i][t]->SetBinContent(k+1,funcSignal[i][k][t]->GetParameter(1));
	      hMean[i][t]->SetBinError(k+1,funcSignal[i][k][t]->GetParError(1));
	      hSigma[i][t]->SetBinContent(k+1,funcSignal[i][k][t]->GetParameter(2));
	      hSigma[i][t]->SetBinError(k+1,funcSignal[i][k][t]->GetParError(2));
	      hFitYield[i][t]->SetBinContent(k+1,fit_yield);
	      hFitYield[i][t]->SetBinError(k+1,fit_yield_err);

	      // Extract background matrix
	      const int nParameter = nPar;
	      double bkg_params[nParameter];
	      double bkg_matrix[nParameter*nParameter];
	      funcBkg[i][k][t] = new TF1(Form("Jpsi_FitBkg_cent%s_%s",Cent_Name[k],tmpName.Data()),funcForm.Data(),g_sig_fit_min,g_sig_fit_max);

	      // from combined fit
	      for(int j=0; j<nPar; j++)
		{
		  funcBkg[i][k][t]->SetParameter(j,funcSignal[i][k][t]->GetParameter(3+j));
		  funcBkg[i][k][t]->SetParError(j,funcSignal[i][k][t]->GetParError(3+j));
		  bkg_params[j] = funcSignal[i][k][t]->GetParameter(3+j);
		}

	      for(int j=3; j<3+nPar; j++)
		{
		  for(int m=3; m<3+nPar; m++)
		    {
		      bkg_matrix[(j-3)*nPar+m-3] = matrix[j*(3+nPar)+m];
		    }
		}

	      // bin counting
	      double low_mass_tmp = low_mass;
	      double high_mass_tmp = high_mass;
	      if(i==0) 
		{
		  low_mass_tmp = 2.96;
		  high_mass_tmp = 3.24;
		}

	      int low_bin = hSignal[i][k][t]->FindFixBin(low_mass_tmp+1e-4);
	      int high_bin = hSignal[i][k][t]->FindFixBin(high_mass_tmp-1e-4);
	      double yield_all_err;
	      double yield_all = hSignal[i][k][t]->IntegralAndError(low_bin,high_bin,yield_all_err);
	      double yield_bkg = funcBkg[i][k][t]->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal[i][k][t]->GetBinWidth(1);
	      double yield_bkg_err = funcBkg[i][k][t]->IntegralError(low_mass_tmp,high_mass_tmp,bkg_params,bkg_matrix)* 1./hSignal[i][k][t]->GetBinWidth(1);
	      double yield_sig = yield_all - yield_bkg;
	      double yield_sig_err = TMath::Sqrt(yield_all_err*yield_all_err+yield_bkg_err*yield_bkg_err);

	      hBinCountYield[i][t]->SetBinContent(k+1,yield_sig);
	      hBinCountYield[i][t]->SetBinError(k+1,yield_sig_err);

	      // plotting
	      cFit->cd(nCentBins-k);	  
	      hSignal[i][k][t]->SetMaximum(2*hSignal[i][k][t]->GetMaximum());
	      hSignal[i][k][t]->Draw();
	      funcSignal[i][k][t]->Draw("sames");
	      funcBkg[i][k][t]->SetLineColor(4);
	      funcBkg[i][k][t]->Draw("sames");
	      text = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],Cent_Name[k]),0.06);
	      text->Draw();
	      TLine *line = GetLine(low_mass_tmp,hSignal[i][k][t]->GetMinimum()*1.5,low_mass_tmp,hSignal[i][k][t]->GetMaximum()*0.3,1);
	      line->Draw();
	      TLine *line = GetLine(high_mass_tmp,hSignal[i][k][t]->GetMinimum()*1.5,high_mass_tmp,hSignal[i][k][t]->GetMaximum()*0.3,1);
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
	    cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Npart/%s%sFit_Jpsi_%s.pdf",run_type,run_cfg_name.Data(),sys_title[isys],tmpName.Data()));
	}
    }

  if(saveHisto)
    {
      if(isys==0)
	{
	  fin->cd();
	  for(int i=0; i<nPtBins; i++)
	    {
	      for(int t=0; t<1; t++)
		{
		  hFitYield[i][t]->Write("",TObject::kOverwrite);
		  hBinCountYield[i][t]->Write("",TObject::kOverwrite);
		}
	    }
	  hJpsiWeight[0]->Write("",TObject::kOverwrite);
	}
      else
	{
	  TString fileName = fin->GetName();
	  fileName.ReplaceAll(".root",".sys.root");
	  TFile *fout = TFile::Open(fileName,"update");
	  for(int i=0; i<nPtBins; i++)
	    {
	      for(int t=0; t<1; t++)
		{
		  cout << hFitYield[i][t]->GetName() << endl;
		  hFitYield[i][t]->Write("",TObject::kOverwrite);
		  hBinCountYield[i][t]->Write("",TObject::kOverwrite);
		}
	    }
	  fout->Close();
	}
    }
}

