const int year = YEAR;

//================================================
void ana_JpsiXsec()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);


  xsec_Run14();
  //compare();
  //trgSetup();
  //shiftDataPoint();
}

//================================================
void trgSetup(const bool savePlot = 0)
{
  //==============================================
  // Cross section vs. pT
  //==============================================
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  // =============================================
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  double mb_events[nCentBins][gNTrgSetup];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  mb_events[i][j] = 0;
	}
    }
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hEqMbEvents[i] = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
      mb_events[i][0] = 0;
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
	  mb_events[i][0] += nEventsRun/rf/nEventsTaken * eq_mb;
	}
      printf("Effective # of MB events for %s%%: %4.4e\n",cent_Name[i],mb_events[i][0]);
    }

  const char *trgSetupName[4] = {"","_low","_mid","_high"};
  for(int j=1; j<gNTrgSetup; j++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type,trgSetupName[j-1]));
      int runnumber;
      while(fruns >> runnumber)
	{
	  int bin = hEvtRunAcc->FindBin(runnumber);
	  if(bin<1 || bin>hEvtRunAcc->GetNbinsX()) continue;
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	  double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(runnumber));
	  if(nEventsTaken==0) 
	    {
	      if(i==0) printf("[w] check run %1.0f\n",runnumber);
	      continue;
	    }
	  double nEventsRun = hEvtRun->GetBinContent(bin);
	  double rf = hRF->GetBinContent(hRF->FindFixBin(runnumber));
	  if(rf==0)
	    {
	      printf("[w] rf = 0 for run %1.0f\n",run);
	      rf = 0.49;
	    }
	  for(int i=0; i<nCentBins; i++)
	    {
	      double eq_mb = hEqMbEvents[i]->GetBinContent(hEqMbEvents[i]->FindFixBin(runnumber));
	      mb_events[i][j] += nEventsRun/rf/nEventsTaken * eq_mb;
	    }
	}
      printf("[i] # of events for %d: %1.0f, %4.2f%%\n",j,mb_events[0][j],mb_events[0][j]/mb_events[0][0]*100);
    }
  printf("[i]Check %1.0f =? %1.0f\n",mb_events[0][0],mb_events[0][1]+mb_events[0][2]+mb_events[0][3]+mb_events[0][4]);
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================

  // =============================================
  // MTD acceptance loss
  const int nRunRange = 8;
  double evtCount[nRunRange+1];
  for(int i=0; i<nRunRange+1; i++) evtCount[i] = 0;
  int runRange[nRunRange+1] = {15074104, 15077035, 15078021, 15098066, 15099002, 15106130, 15131038, 15132019, 15167014};
  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
    {
      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
      double run = hEvtRunAcc->GetBinCenter(bin);
      double evt = hEvtRunAcc->GetBinContent(bin);
      for(int i=0; i<nRunRange; i++)
	{
	  if(run>runRange[i] && run<=runRange[i+1]) evtCount[i] += evt;
	}
      evtCount[nRunRange] += evt;
    }
  TList *list = new TList;
  TString legName[nRunRange+1];
  legName[nRunRange] = "Average correction factor";
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"read");
  TH1F *hAccLoss[nRunRange];
  for(int i=0; i<nRunRange; i++)
    {
      hAccLoss[i] = (TH1F*)fAcc->Get(Form("hAccepLoss_RunRange%d",i));
    }
  TH1F *hAccCorr = (TH1F*)hAccLoss[0]->Clone("hAccCorr");
  hAccCorr->Reset();
  int colors[10] = {1, 2, 4, 6, 1, 2, 4, 6, 1, 2};
  for(int i=0; i<nRunRange; i++)
    {
      hAccCorr->Add(hAccLoss[i],evtCount[i]/evtCount[nRunRange]);
      list->Add(hAccLoss[i]);
      hAccLoss[i]->SetMarkerStyle(20+i);
      hAccLoss[i]->SetMarkerColor(colors[i]);
      hAccLoss[i]->SetLineColor(colors[i]);
      legName[i] = Form("Run (%d,%d]: %2.1f%%",runRange[i],runRange[i+1],evtCount[i]/evtCount[nRunRange]*100);
    }
  list->Add(hAccCorr);
  c = drawHistos(list,"MTD_Acceptance","Correction for MTD acceptance loss",kFALSE,0,30,true,0.5,1.02,kFALSE,kTRUE,legName,kTRUE,"",0.35,0.75,0.15,0.45,kTRUE,0.04,0.04,false,1,false,false);
  // =============================================
  //
  //

  // Jpsi efficiency
  // additional correction for trigger efficiency is needed
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.MtdTrigEff.root"),"read");
  TH1F *hJpsiTrigEff[3];
  hJpsiTrigEff[0] = (TH1F*)fTrigEff->Get("Run14_AuAu200_JpsiEffVsPt_TacDiffEff_TrigStudy");
  hJpsiTrigEff[1] = (TH1F*)fTrigEff->Get("Run14_AuAu200_JpsiEffVsPt_TacDiffEff_prod_low_TrigStudy");
  hJpsiTrigEff[2] = (TH1F*)fTrigEff->Get("Run14_AuAu200_JpsiEffVsPt_TacDiffEff_prod_high_TrigStudy");
  hJpsiTrigEff[1]->Divide(hJpsiTrigEff[0]);
  hJpsiTrigEff[2]->Divide(hJpsiTrigEff[0]);


  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiEff[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiEff[k][j] = (TH1F*)fEff->Get(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[j]));
	  if(j>=1 && j<=3) hJpsiEff[k][j]->Multiply(hJpsiTrigEff[1]);
	  if(j==4) hJpsiEff[k][j]->Multiply(hJpsiTrigEff[2]);
	}
    }
  

  // Jpsi raw counts
  char * yieldName = Form("Run14_AuAu200.JpsiYield.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiCounts[k][j] = (TH1F*)fYield->Get(Form("Jpsi_FitYield_cent%s_weight%s",cent_Title[k],gTrgSetupName[j]));
	}
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins][gNTrgSetup];
  TH1F *hJpsiInvRatio[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[k][j] = (TH1F*)hJpsiCounts[k][j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[k],gTrgSetupTitle[j]));
	  hJpsiInvYield[k][j]->Divide(hJpsiEff[k][j]);
	  hJpsiInvYield[k][j]->Divide(hAccCorr);
	  cout << hJpsiEff[k][j]->GetBinContent(1) << endl;
	  for(int bin=1; bin<=hJpsiInvYield[k][j]->GetNbinsX(); bin++)
	    {
	      double bin_width = hJpsiInvYield[k][j]->GetBinWidth(bin); // dpT
	      double bin_center = hJpsiInvYield[k][j]->GetBinCenter(bin); // pT 
	      if(bin==1)
		{
		  bin_width = 0.85;
		  bin_center = 0.575;
		}
	      hJpsiInvYield[k][j]->SetBinContent(bin,hJpsiInvYield[k][j]->GetBinContent(bin)/bin_width/bin_center);
	      hJpsiInvYield[k][j]->SetBinError(bin,hJpsiInvYield[k][j]->GetBinError(bin)/bin_width/bin_center);
	    }
	  hJpsiInvYield[k][j]->Scale(1./mb_events[k][j]); // N_evt
	  hJpsiInvYield[k][j]->Scale(1./(2*pi)); // 2pi
	  hJpsiInvYield[k][j]->SetMarkerStyle(21);
	  hJpsiInvYield[k][j]->SetMarkerColor(color[j]);
	  hJpsiInvYield[k][j]->SetLineColor(color[j]);
	  hJpsiInvYield[k][j]->SetMarkerSize(1.2);

	  hJpsiInvRatio[k][j] = (TH1F*)hJpsiInvYield[k][j]->Clone(Form("Jpsi_InvYieldRatio_cent%s%s",cent_Title[k],gTrgSetupTitle[j]));
	  hJpsiInvRatio[k][j]->Divide(hJpsiInvYield[k][0]);
	  for(int bin=1; bin<=hJpsiInvRatio[k][j]->GetNbinsX(); bin++)
	    {
	      double err = hJpsiInvYield[k][j]->GetBinError(bin)/hJpsiInvYield[k][j]->GetBinContent(bin);
	      hJpsiInvRatio[k][j]->SetBinError(bin,hJpsiInvRatio[k][j]->GetBinContent(bin)*err);
	    }
	}
    }

  TCanvas *c = new TCanvas("Comp_JetXsec","Comp_JetXsec",1100,700);
  c->Divide(3,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",15,0,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-11,1e-4);
  ScaleHistoTitle(hAuAu,0.055,1,0.045,0.055,1.4,0.045,62);
  TLegend *leg = new TLegend(0.1,0.6,0.4,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05); 
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      if(k==2 || k==3) hAuAu->GetXaxis()->SetRangeUser(0,10);
      if(k==4) hAuAu->GetXaxis()->SetRangeUser(0,6);
      hAuAu->DrawCopy();
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[k][j]->Draw("sames");
	  if(k==0) leg->AddEntry(hJpsiInvYield[k][j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
	}
      TPaveText *t1 = GetPaveText(0.25,0.35,0.2,0.3,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_TrgSetupComp.pdf",run_type,run_config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_InvYield_CompLumi.pdf"));
    }


  TCanvas *c = new TCanvas("Ratio_JetXsec","Ratio_JetXsec",1100,700);
  c->Divide(3,2);
  hAuAu->SetYTitle("Ratio to combined");
  hAuAu->GetYaxis()->SetRangeUser(0,2.5);
  TLegend *leg = new TLegend(0.1,0.6,0.4,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05); 
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      if(k==0 || k==1) hAuAu->GetXaxis()->SetRangeUser(0.15,15);
      if(k==2 || k==3) hAuAu->GetXaxis()->SetRangeUser(0.15,10);
      if(k==4) hAuAu->GetXaxis()->SetRangeUser(0.15,6);
      hAuAu->DrawCopy();
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hJpsiInvRatio[k][j]->Draw("sames");
	  if(k==0) leg->AddEntry(hJpsiInvRatio[k][j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
	}
      TPaveText *t1 = GetPaveText(0.25,0.35,0.8,0.9,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_TrgSetupRatio.pdf",run_type,run_config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_InvYield_RatioLumi.pdf"));
    }

}

//================================================
void xsec_Run14(const bool savePlot = 0, const bool saveHisto = 0)
{

  //==============================================
  // Cross section vs. pT
  //==============================================
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;


  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
  const double dimuon_events = hStat->GetBinContent(10);

  // =============================================
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  double mb_events[nCentBins];
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hEqMbEvents[i] = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
      mb_events[i] = 0;
      for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
	{
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	  double run = hEvtRunAcc->GetBinCenter(bin);
	  double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(run));
	  if(nEventsTaken==0) 
	    {
	      if(i==0) printf("[w] check run %1.0f, # of analyzed = %1.0f\n",run,hEvtRun->GetBinContent(bin));
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
	  mb_events[i] += nEventsRun/rf/nEventsTaken * eq_mb;

	}
      printf("Effective # of MB events for %s%%: %4.4e\n",cent_Name[i],mb_events[i]);
    }
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================
  //
  //

  // =============================================
  // MTD acceptance loss
  const int nRunRange = 7;
  double evtCount[nRunRange+1];
  for(int i=0; i<nRunRange+1; i++) evtCount[i] = 0;
  int runRange[nRunRange+1] = {15074104, 15078034, 15098066, 15099003, 15106130, 15131036, 15132019, 15167014};
  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
    {
      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
      double run = hEvtRunAcc->GetBinCenter(bin);
      double evt = hEvtRunAcc->GetBinContent(bin);
      for(int i=0; i<nRunRange; i++)
	{
	  if(run>runRange[i] && run<=runRange[i+1]) evtCount[i] += evt;
	}
      evtCount[nRunRange] += evt;
    }
  TList *list = new TList;
  TString legName[nRunRange+1];
  legName[nRunRange] = "Average correction factor";
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"read");
  TH1F *hAccLoss[nRunRange];
  for(int i=0; i<nRunRange; i++)
    {
      hAccLoss[i] = (TH1F*)fAcc->Get(Form("hAccepLoss_RunRange%d",i+1));
    }
  TH1F *hAccCorr = (TH1F*)hAccLoss[0]->Clone("hAccCorr");
  hAccCorr->Reset();
  int colors[10] = {1, 2, 4, 6, 1, 2, 4, 6, 1, 2};
  for(int i=0; i<nRunRange; i++)
    {
      hAccCorr->Add(hAccLoss[i],evtCount[i]/evtCount[nRunRange]);
      list->Add(hAccLoss[i]);
      hAccLoss[i]->SetMarkerStyle(20+i);
      hAccLoss[i]->SetMarkerColor(colors[i]);
      hAccLoss[i]->SetLineColor(colors[i]);
      legName[i] = Form("Run (%d,%d]: %2.1f%%",runRange[i],runRange[i+1],evtCount[i]/evtCount[nRunRange]*100);
    }
  list->Add(hAccCorr);
  c = drawHistos(list,"MTD_Acceptance","Correction for MTD acceptance loss",kFALSE,0,30,true,0.5,1.02,kFALSE,kTRUE,legName,kTRUE,"",0.35,0.75,0.15,0.45,kTRUE,0.04,0.04,false,1,false,false);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MTD_AccepLoss.pdf",run_type));
    }
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_MtdAccLossAvg.pdf"));
    }
  // =============================================
  //
  //

  // Jpsi efficiency
  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEff[k] = (TH1F*)fEff->Get(Form("JpsiEffVsPt_cent%s_final",cent_Title[k]));
      printf("cent %s\n",cent_Title[k]);
      for(int bin=1; bin<=hJpsiEff[k]->GetNbinsX(); bin++)
	{
	  printf("[i] pt = %4.2f, eff = %4.5f\n",hJpsiEff[k]->GetBinCenter(bin),hJpsiEff[k]->GetBinContent(bin));
	}
    }
  

  // Jpsi raw counts
  char * yieldName = Form("Run14_AuAu200.JpsiYield.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiCounts[k] = (TH1F*)fYield->Get(Form("Jpsi_FitYield_cent%s_weight",cent_Title[k]));
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)hJpsiCounts[k]->Clone(Form("Jpsi_InvYieldVsPt_cent%s",cent_Title[k]));
      hJpsiInvYield[k]->Divide(hJpsiEff[k]);
      hJpsiInvYield[k]->Divide(hAccCorr);
      cout << hJpsiEff[k]->GetBinContent(1) << endl;
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  double bin_width = hJpsiInvYield[k]->GetBinWidth(bin); // dpT
	  double bin_center = hJpsiInvYield[k]->GetBinCenter(bin); // pT 
	  if(bin==1)
	    {
	      bin_width = 0.85;
	      bin_center = 0.575;
	    }
	  hJpsiInvYield[k]->SetBinContent(bin,hJpsiInvYield[k]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[k]->SetBinError(bin,hJpsiInvYield[k]->GetBinError(bin)/bin_width/bin_center);
	}

      hJpsiInvYield[k]->Scale(1./mb_events[k]); // N_evt
      hJpsiInvYield[k]->Scale(1./(2*pi)); // 2pi
      hJpsiInvYield[k]->SetMarkerStyle(21);
      hJpsiInvYield[k]->SetMarkerColor(2);
      hJpsiInvYield[k]->SetLineColor(2);
      hJpsiInvYield[k]->SetMarkerSize(1.5);
    }


  // published results
  TFile *fpub = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");
  TGraphAsymmErrors *gAuAuLowPt[3];
  TGraphAsymmErrors *gAuAuLowPtSys[3];
  TGraphAsymmErrors *gAuAuHighPt[3];
  TGraphAsymmErrors *gAuAuHighPtSys[3];
  TH1F *hTBW[3];
  for(int i=0; i<3; i++)
    {
      gAuAuLowPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[i+1]));
      gAuAuLowPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[i+1]));
      gAuAuHighPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[i+1]));
      gAuAuHighPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[i+1]));
      hTBW[i] = (TH1F*)fpub->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[i+1]));
    }

  // sQM2016
  TFile *fsQM = TFile::Open("Rootfiles/2016sQM.root","read");
  TGraphErrors *gAuAuYield_sQM[3];
  TGraphErrors *gAuAuYieldSys_sQM[3];
  for(int i=0; i<3; i++)
    {
      gAuAuYield_sQM[i] = (TGraphErrors*)fsQM->Get(Form("Graph_Jpsi_InvYield_cent%s",cent_Title[i+1]));
      gAuAuYieldSys_sQM[i] = (TGraphErrors*)fsQM->Get(Form("Graph_Jpsi_InvYield_cent%s_sys",cent_Title[i+1]));
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(3,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",15,0.15,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-11,1e-4);
  hAuAu->GetYaxis()->SetNdivisions(505);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1.2,0.05,62);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.16,0.02,0.02);
      hAuAu->Draw();

      if(k>0 && k<4)
	{
	  gAuAuLowPt[k-1]->SetMarkerStyle(24);
	  gAuAuLowPt[k-1]->SetMarkerColor(1);
	  gAuAuLowPt[k-1]->SetLineColor(1);
	  gAuAuLowPt[k-1]->Draw("sames PE");
	  gAuAuLowPtSys[k-1]->SetMarkerColor(1);
	  gAuAuLowPtSys[k-1]->SetLineColor(1);
	  gAuAuLowPtSys[k-1]->Draw("sameE5");
	  gAuAuHighPt[k-1]->SetMarkerStyle(26);
	  gAuAuHighPt[k-1]->SetMarkerColor(1);
	  gAuAuHighPt[k-1]->SetLineColor(1);
	  gAuAuHighPt[k-1]->Draw("sames PE");
	  gAuAuHighPtSys[k-1]->SetMarkerColor(1);
	  gAuAuHighPtSys[k-1]->SetLineColor(1);
	  gAuAuHighPtSys[k-1]->Draw("sameE5");
	  hTBW[k-1]->Draw("sames");

	  gAuAuYield_sQM[k-1]->SetMarkerStyle(29);
	  gAuAuYield_sQM[k-1]->SetMarkerColor(kGreen+2);
	  gAuAuYield_sQM[k-1]->SetLineColor(kGreen+2);
	  //gAuAuYield_sQM[k-1]->Draw("sames PE");
	}

      hJpsiInvYield[k]->Draw("sames");
      TPaveText *t1 = GetPaveText(0.7,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }

  c->cd(6);
  TLegend *leg = new TLegend(0.1,0.3,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->AddEntry(gAuAuLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (Low p_{T})","P");
  leg->AddEntry(gAuAuHighPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (High p_{T})","P");
  leg->AddEntry(hTBW[0],"TBW fit (#beta=0) to J/#psi#rightarrowe^{+}e^{-}","L");
  leg->AddEntry(hJpsiInvYield[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");
  //leg->AddEntry(gAuAuYield_sQM[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5 (sQM2016)","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_compareToPub.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_compareToPub.png",run_type,run_config));
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi_ratio","AuAu200_Jpsi_ratio",1100,700);
  c->Divide(2,2);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  double x,y;
  for(int k=0; k<3; k++)
    {
      TH1F *hRatio = (TH1F*)hJpsiInvYield[k+1]->Clone(Form("%s_ratio",hJpsiInvYield[k]->GetName()));
      int endbin = hRatio->FindFixBin(9);

      TH1F *hRatio_sQM = (TH1F*)hRatio->Clone(Form("%s_clone",hRatio->GetName()));
      hRatio_sQM->Reset();
      for(int bin=1; bin<=endbin; bin++)
	{
	  int start_bin = hTBW[k]->FindBin(hJpsiInvYield[k]->GetXaxis()->GetBinLowEdge(bin)+1e-6);
	  int end_bin   = hTBW[k]->FindBin(hJpsiInvYield[k]->GetXaxis()->GetBinUpEdge(bin)-1e-6);
	  double scale = 0;
	  for(int ibin=start_bin; ibin<=end_bin; ibin++)
	    {
	      scale += hTBW[k]->GetBinContent(ibin) * hTBW[k]->GetBinCenter(ibin) * hTBW[k]->GetBinWidth(ibin);
	    }
	  scale = scale / hRatio->GetBinCenter(bin) / hRatio->GetBinWidth(bin);
	  hRatio->SetBinContent(bin,hRatio->GetBinContent(bin)/scale);
	  hRatio->SetBinError(bin,hRatio->GetBinError(bin)/scale);

	  //gAuAuYield_sQM[k]->GetPoint(bin-1, x, y);
	  //hRatio_sQM->SetBinContent(bin,y/scale);
	  //hRatio_sQM->SetBinError(bin,gAuAuYield_sQM[k]->GetErrorYhigh(bin-1)/scale);
	  
	}
      c->cd(k+1);
      gPad->SetGridy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
 
      // final
      hRatio->SetTitle(";p_{T} (GeV/c);Ratio to TBW");
      hRatio->GetXaxis()->SetRangeUser(0.15,10);
      hRatio->GetYaxis()->SetRangeUser(0,2);
      ScaleHistoTitle(hRatio,0.06,1,0.05,0.06,1,0.05,62);
      hRatio->Draw();

      /*
      // sQM2016
      hRatio_sQM->SetMarkerStyle(29);
      hRatio_sQM->SetMarkerColor(kGreen+2);
      hRatio_sQM->SetLineColor(kGreen+2);
      hRatio_sQM->SetMarkerSize(2);
      hRatio_sQM->Draw("samesP");
      */
      // published
      TGraphAsymmErrors *gRatioLowPt = (TGraphAsymmErrors*)gAuAuLowPt[k]->Clone(Form("%s_ratio",gAuAuLowPt[k]->GetName()));
      for(int ipoint=0; ipoint<gRatioLowPt->GetN(); ipoint++)
	{
	  gAuAuLowPt[k]->GetPoint(ipoint, x, y);
	  double scale = hTBW[k]->GetBinContent(hTBW[k]->FindFixBin(x));
	  gRatioLowPt->SetPoint(ipoint, x, y/scale);
	  gRatioLowPt->SetPointError(ipoint, gAuAuLowPt[k]->GetErrorXlow(ipoint), gAuAuLowPt[k]->GetErrorXhigh(ipoint),
				     gAuAuLowPt[k]->GetErrorYlow(ipoint)/scale, gAuAuLowPt[k]->GetErrorYhigh(ipoint)/scale);
	}
      gRatioLowPt->Draw("samesPEZ");

      TGraphAsymmErrors *gRatioHighPt = (TGraphAsymmErrors*)gAuAuHighPt[k]->Clone(Form("%s_ratio",gAuAuHighPt[k]->GetName()));
      for(int ipoint=0; ipoint<gRatioHighPt->GetN(); ipoint++)
	{
	  gAuAuHighPt[k]->GetPoint(ipoint, x, y);
	  double scale = hTBW[k]->GetBinContent(hTBW[k]->FindFixBin(x));
	  gRatioHighPt->SetPoint(ipoint, x, y/scale);
	  gRatioHighPt->SetPointError(ipoint, gAuAuHighPt[k]->GetErrorXlow(ipoint), gAuAuHighPt[k]->GetErrorXhigh(ipoint),
				     gAuAuHighPt[k]->GetErrorYlow(ipoint)/scale, gAuAuHighPt[k]->GetErrorYhigh(ipoint)/scale);
	}
      gRatioHighPt->SetMarkerStyle(26);
      gRatioHighPt->Draw("samesPEZ");      

      TPaveText *t1 = GetPaveText(0.2,0.3,0.85,0.9,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k+1]));
      t1->Draw();
    }
  c->cd(4);
  TLegend *leg = new TLegend(0.1,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->AddEntry(gAuAuLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (Low p_{T})","P");
  leg->AddEntry(gAuAuHighPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (High p_{T})","P");
  leg->AddEntry(hJpsiInvYield[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");
  //leg->AddEntry(gAuAuYield_sQM[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5 (sQM2016)","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_RatioToTBW.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_RatioToTBW.png",run_type,run_config));
    }
  return;

  // Raa distribution
  TFile *fppRef = TFile::Open("Rootfiles/Paper.Run14_AuAu200.Jpsi.root","read");
  const double ppInelastic = 42.; // mb
  double ncoll[nCentBins] = {291.9, 766.47, 290.87, 91.33, 21.57};
  TH1F *hJpsipp = (TH1F*)fppRef->Get("hpp200JpsiVsPtFinal");
  TH1F *hJpsiRaa[nCentBins];
  TCanvas *c = new TCanvas("AuAu200_Jpsi_Raa","AuAu200_Jpsi_Raa",1100,700);
  c->Divide(3,2);
  for(int k=0; k<nCentBins; k++)
    {
      int nbins = hJpsiInvYield[k]->GetNbinsX();
      hJpsiRaa[k] = (TH1F*)hJpsiInvYield[k]->Clone(Form("Jpsi_RaaVsPt_cent%s",cent_Title[k]));
      hJpsiRaa[k]->Reset();
      for(int bin=1; bin<=nbins; bin++)
	{
	  double AuAu_val = hJpsiInvYield[k]->GetBinContent(bin);
	  double AuAu_err = hJpsiInvYield[k]->GetBinError(bin);
	  double pp_val = hJpsipp->GetBinContent(bin);
	  double pp_err = hJpsipp->GetBinError(bin);
	  double prefix = ppInelastic/ncoll[k] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  double err = prefix * AuAu_err / pp_val;
	  hJpsiRaa[k]->SetBinContent(bin, val);
	  hJpsiRaa[k]->SetBinError(bin, err);
	}

      c->cd(k+1);
      if(k==2 || k==3) hJpsiRaa[k]->GetXaxis()->SetRangeUser(0.15,10);
      if(k==4) hJpsiRaa[k]->GetXaxis()->SetRangeUser(0.15,6);
      hJpsiRaa[k]->GetYaxis()->SetRangeUser(0,1.5);
      hJpsiRaa[k]->SetTitle(";p_{T} (GeV/c);R_{AA}");
      hJpsiRaa[k]->Draw();
      TPaveText *t1 = GetTitleText(Form("%s%%",cent_Name[k]),0.06);
      t1->Draw();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsPt.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsPt.png",run_type,run_config));
    }


  // compare Raa distribtion using different references
  const int kcent = 1;
  TH1F *hJpsiRaaVsPt[3];
  TGraphAsymmErrors *gJPsiRaaVsPtSys[3];
  TGraphAsymmErrors *gppRefSys[3];
  double x,y;
  for(int i=0; i<3; i++)
    {
      int nbins = hJpsiInvYield[kcent]->GetNbinsX();
      hJpsiRaaVsPt[i] = (TH1F*)hJpsiInvYield[kcent]->Clone(Form("Jpsi_RaaVsPt_cent%s_ppRef%d",cent_Title[kcent],i+1));
      hJpsiRaaVsPt[i]->Reset();
      hJpsiRaaVsPt[i]->SetTitle(Form(";p_{T} (GeV/c);R_{AA}"));
      gJPsiRaaVsPtSys[i] = new TGraphAsymmErrors(nbins);

      // pp reference
      if(i==0) gppRefSys[i] = (TGraphAsymmErrors*)fppRef->Get("hpp200JpsiVsPtFinalSys");
      if(i==1) gppRefSys[i] = (TGraphAsymmErrors*)fppRef->Get("hpp200JpsiVsPtFinalSysTot_Run15_dimuon");
      if(i==2) gppRefSys[i] = (TGraphAsymmErrors*)fppRef->Get("hpp200JpsiVsPtFinalSysTot_Run12_dielectron");
      for(int bin=1; bin<=nbins; bin++)
	{
	  if(i==1 && bin==nbins) continue;
	  gppRefSys[i]->GetPoint(bin-1, x, y);
	  double AuAu_val = hJpsiInvYield[kcent]->GetBinContent(bin);
	  double AuAu_err = hJpsiInvYield[kcent]->GetBinError(bin);
	  double dpt = hJpsiInvYield[kcent]->GetBinWidth(bin);
	  double pp_val = y;
	  double pp_err_h = gppRefSys[i]->GetErrorYhigh(bin-1);
	  double pp_err_l = gppRefSys[i]->GetErrorYlow(bin-1);
	  double prefix = ppInelastic/ncoll[kcent] * 1e6;
	  double val = prefix * AuAu_val / pp_val;
	  double err = prefix * AuAu_err / pp_val;
	  hJpsiRaaVsPt[i]->SetBinContent(bin, val);
	  hJpsiRaaVsPt[i]->SetBinError(bin, err);
	  gJPsiRaaVsPtSys[i]->SetPoint(bin-1, x, val);
	  gJPsiRaaVsPtSys[i]->SetPointError(bin-1, dpt/2, dpt/2, pp_err_h/pp_val*val, pp_err_l/pp_val*val);
	}
      int clr = color[i];
      hJpsiRaaVsPt[i]->SetMarkerStyle(20+i);
      hJpsiRaaVsPt[i]->SetMarkerSize(1.2);
      hJpsiRaaVsPt[i]->SetMarkerColor(clr);
      hJpsiRaaVsPt[i]->SetLineColor(clr);
      hJpsiRaaVsPt[i]->GetYaxis()->SetRangeUser(0,1);
      gJPsiRaaVsPtSys[i]->SetLineColor(clr);
      gJPsiRaaVsPtSys[i]->SetLineWidth(1);
    }
  
  TCanvas *c = new TCanvas("JpsiRaa_CompRef","JpsiRaa_CompRef", 1200, 500);
  c->Divide(3, 1);
  const char *refType[3] = {"Published STAR+PHENIX", "Run15_dimuon","Run12_dielectron"};
  for(int i=0; i<3; i++)
    {
      c->cd(i+1);
      hJpsiRaaVsPt[i]->Draw("");
      gJPsiRaaVsPtSys[i]->SetFillStyle(0);
      gJPsiRaaVsPtSys[i]->Draw("sameE5");
      TLegend *leg = new TLegend(0.15,0.15,0.35,0.25);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(62);
      leg->SetTextSize(0.035);
      leg->AddEntry(hJpsiRaaVsPt[i],Form("pp reference: %s",refType[i]),"PL");
      leg->AddEntry(gJPsiRaaVsPtSys[i],"Total stat.+sys. error from reference","F");
      leg->Draw();
      TPaveText *t1 = GetTitleText(Form("J/psi R_{AA} in %s%%",cent_Name[kcent]),0.06);
      t1->Draw();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsPt_CompRef.pdf",run_type,run_config));
    }

  //==============================================
  // Cross section vs. centrality
  //==============================================

  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  const int kNCent = nCentBins_npart[0];
  double mb_events_npart[nPtBins_npart][kNCent];
  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int k=0; k<nCentBins_npart[i]; k++)
	{
	  TH1F *hEqMbEvts = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title_npart[i*kNCent+k]));
	  mb_events_npart[i][k] = 0;
	  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
	    {
	      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	      double run = hEvtRunAcc->GetBinCenter(bin);
	      double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(run));
	      if(nEventsTaken==0) 
		{
		  if(i==0 && k==0) printf("[w] check run %1.0f\n",run);
		  continue;
		}
	      double nEventsRun = hEvtRun->GetBinContent(bin);
	      double rf = hRF->GetBinContent(hRF->FindFixBin(run));
	      if(rf==0)
		{
		  printf("[w] rf = 0 for run %1.0f\n",run);
		  rf = 0.49;
		}
	      double eq_mb = hEqMbEvts->GetBinContent(hEqMbEvts->FindFixBin(run));
	      mb_events_npart[i][k] += nEventsRun/rf/nEventsTaken * eq_mb;
	    }
	  if(i==0) printf("Effective # of MB events for %s%%: %4.4e\n",cent_Name_npart[i*kNCent+k],mb_events_npart[i][k]);
	}
    }
  
  printf("++++++++++MTD acceptance loss+++++++++++++\n");
  // =============================================
  // MTD acceptance loss
  TFile *fWeight = TFile::Open("Rootfiles/Run14_AuAu200.Input.root","read");
  TH1F *hMcJpsiPt = (TH1F*)fWeight->Get(Form("hInputJpsiShape_Cent0"));
  double mtd_acc[nPtBins_npart] = {0,0};
  for(int i=0; i<nPtBins_npart; i++)
    {
      int low_bin = hAccCorr->FindFixBin(ptBins_low_npart[i]+1e-4);
      int high_bin = hAccCorr->FindFixBin(ptBins_high_npart[i]-1e-4);
      for(int bin=low_bin; bin<=high_bin; bin++)
	{
	  double bin1 = hMcJpsiPt->FindFixBin(hAccCorr->GetXaxis()->GetBinLowEdge(bin)+1e-4);
	  double bin2 = hMcJpsiPt->FindFixBin(hAccCorr->GetXaxis()->GetBinUpEdge(bin)-1e-4);
	  mtd_acc[i] += hAccCorr->GetBinContent(bin) * hMcJpsiPt->Integral(bin1,bin2);
	}
      mtd_acc[i] /= hMcJpsiPt->Integral(hMcJpsiPt->FindFixBin(ptBins_low_npart[i]),hMcJpsiPt->FindFixBin(ptBins_high_npart[i]));
      printf("[i] MTD acceptance correction is %4.3f for pt > %1.0f\n",mtd_acc[i],ptBins_low_npart[i]);
    }
  
  // =============================================
  // cross-section vs. centrality
  TH1F *hJpsiCounts_npart[nPtBins_npart];
  TH1F *hJpsiEff_npart[nPtBins_npart];
  TH1F *hJpsiInvYield_npart[nPtBins_npart];
  for(int i=0; i<nPtBins_npart; i++)
    {
      int numCentBins = nCentBins_npart[i];
      hJpsiEff_npart[i] = (TH1F*)fEff->Get(Form("JpsiEffVsCent_Pt%1.0f_final",ptBins_low_npart[i]));
      for(int bin=1; bin<=numCentBins; bin++)
	{
	  printf("[i] Eff for pt > %1.0f, cent %s: %8.4f%%\n",ptBins_low_npart[i],cent_Name_npart[bin-1],hJpsiEff_npart[i]->GetBinContent(bin)*100);
	}
      hJpsiCounts_npart[i] = (TH1F*)fYield->Get(Form("Jpsi_FitYield_pt%s_weight",pt_Name_npart[i]));
      TH1F *htmpYield = (TH1F*)hJpsiCounts_npart[i]->Clone(Form("%s_clone",hJpsiCounts_npart[i]->GetName()));
      htmpYield->Divide(hJpsiEff_npart[i]);
      htmpYield->Scale(1./mtd_acc[i]);
      htmpYield->Scale(1./(2*pi)); // 2pi

      hJpsiInvYield_npart[i] = new TH1F(Form("Jpsi_InvYieldVsCent_pt%s",pt_Name_npart[i]),"",numCentBins,0,numCentBins);
      for(int bin=1; bin<=numCentBins; bin++)
	{
	  hJpsiInvYield_npart[i]->SetBinContent(bin, htmpYield->GetBinContent(numCentBins-bin+1)/mb_events_npart[i][numCentBins-bin]);
	  hJpsiInvYield_npart[i]->SetBinError(bin, htmpYield->GetBinError(numCentBins-bin+1)/mb_events_npart[i][numCentBins-bin]);
	  hJpsiInvYield_npart[i]->GetXaxis()->SetBinLabel(bin, Form("%s%%",cent_Name_npart[i*kNCent+numCentBins-bin]));
	}
      c = draw1D( hJpsiInvYield[i]);
    }


  TGraphErrors *gJpsiRaaVsCent_sQM[nPtBins_npart];
  TH1F *hJpsiRaaVsCent_sQM[nPtBins_npart];
  double ncoll_part[nPtBins_npart][kNCent] = {941.2, 593.7, 366.4, 216.5, 120.1, 62.2, 29.7, 13.1,
					      941.2, 593.7, 366.4, 216.5, 120.1, 62.2, 21.6, 29.7};
  double ncoll_sQM[6]     = {964, 609, 377, 224, 124, 64};
  TH1F *hppNpart = (TH1F*)fppRef->Get("hpp200JpsiVsCentFinal");
  TH1F *hJpsiRaaVsCent[nPtBins_npart];
  double x,y;
  for(int i=0; i<nPtBins_npart; i++)
    {
      int numCentBins = nCentBins_npart[i];
      gJpsiRaaVsCent_sQM[i] = (TGraphErrors*)fsQM->Get(Form("MTD_Run14AuAu_JpsiRaaVsNpart%1.0fGeV",ptBins_low_npart[i]));
      hJpsiRaaVsCent_sQM[i] = new TH1F(Form("hJpsiRaaVsCent_sQM_pt%s",pt_Name_npart[i]),"",numCentBins,0,numCentBins);
      for(int ipoint=0; ipoint<gJpsiRaaVsCent_sQM[i]->GetN(); ipoint++)
	{
	  gJpsiRaaVsCent_sQM[i]->GetPoint(ipoint, x, y);
	  hJpsiRaaVsCent_sQM[i]->SetBinContent(ipoint+3-i, y/ncoll_sQM[ipoint]*ncoll_part[i][ipoint]);
	  hJpsiRaaVsCent_sQM[i]->SetBinError(ipoint+3-i, gJpsiRaaVsCent_sQM[i]->GetErrorY(ipoint)/ncoll_sQM[ipoint]*ncoll_part[i][ipoint]);
	}
      hJpsiRaaVsCent_sQM[i]->SetMarkerStyle(29);
      hJpsiRaaVsCent_sQM[i]->SetMarkerColor(kGreen+2);
      hJpsiRaaVsCent_sQM[i]->SetLineColor(kGreen+2);
      hJpsiRaaVsCent_sQM[i]->SetMarkerSize(2);

      hJpsiRaaVsCent[i] = (TH1F*)hJpsiInvYield_npart[i]->Clone(Form("Jpsi_RaaVsCent_pt%s",pt_Name_npart[i]));
      double pp_yield = hppNpart->GetBinContent(i+1);
      for(int bin=1; bin<=numCentBins; bin++)
	{
	  double auau_yield = hJpsiInvYield_npart[i]->GetBinContent(bin) * 2*pi;
	  double auau_err = hJpsiInvYield_npart[i]->GetBinError(bin) * 2*pi;
	  double raa = ppInelastic/ncoll_part[i][numCentBins-bin] * 1e6 * auau_yield/pp_yield;
	  double raa_err = auau_err/auau_yield * raa;
	  hJpsiRaaVsCent[i]->SetBinContent(bin, raa);
	  hJpsiRaaVsCent[i]->SetBinError(bin, raa_err);
	}
      hJpsiRaaVsCent[i]->SetMarkerStyle(21);
      hJpsiRaaVsCent[i]->SetMarkerColor(2);
      hJpsiRaaVsCent[i]->SetLineColor(2);
      hJpsiRaaVsCent[i]->SetMarkerSize(1.5);
      hJpsiRaaVsCent[i]->GetYaxis()->SetRangeUser(0,1.2);
      hJpsiRaaVsCent[i]->GetXaxis()->SetLabelSize(0.05);
      TCanvas *c = draw1D(hJpsiRaaVsCent[i],Form("%s: J/#Psi R_{AA} vs. centrality for %s;;R_{AA}",run_type,pt_Title_npart[i]));
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsCent_Pt%1.0f.pdf",run_type,run_config,ptBins_low_npart[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsCent_Pt%1.0f.png",run_type,run_config,ptBins_low_npart[i]));
	}

      hJpsiRaaVsCent_sQM[i]->Draw("samesP");
      TLegend *leg = new TLegend(0.3,0.2,0.5,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(62);
      leg->SetTextSize(0.04);
      leg->AddEntry(hJpsiRaaVsCent[i],"Latest analysis","PL");
      leg->AddEntry(hJpsiRaaVsCent_sQM[i],"2016sQM","PL");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsCent_CompSQM2016_Pt%1.0f.pdf",run_type,run_config,ptBins_low_npart[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiRaaVsCent_CompSQM2016_Pt%1.0f.png",run_type,run_config,ptBins_low_npart[i]));
	}
	  
    }

  if(saveHisto)
    {
      char *outname = Form("Run14_AuAu200.JpsiXsec.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outname),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiInvYield[k]->Write();
	  hJpsiRaa[k]->Write();
	}
      for(int i=0; i<nPtBins_npart; i++)
	{
	  hJpsiInvYield_npart[i]->Write();
	  hJpsiRaaVsCent[i]->Write();
	}
    }
}

//===============================================
void shiftDataPoint(const bool savePlot = 1, const bool saveHisto = 1)
{
  //==============================================
  // Cross section vs. pT
  //==============================================
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  char *filename = Form("Run14_AuAu200.JpsiXsec.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s",filename),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s",filename),"read");
 
  TH1F *hJpsiYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiYield[k] = (TH1F*)fin->Get(Form("Jpsi_InvYieldVsPt_cent%s",cent_Title[k]));
    }

  // add uncertainty from signal extraction to the statistical error
  TFile *ferr = TFile::Open(Form("Rootfiles/%s.Sys.JpsiYield.root",run_type));
  TH1F *hJpsiSigExtErr[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiSigExtErr[k] = (TH1F*)ferr->Get(Form("Sys_signalExt_%s",cent_Title[k]));
      for(int bin=1; bin<=hJpsiYield[k]->GetNbinsX(); bin++)
	{
	  double err1 = hJpsiYield[k]->GetBinError(bin)/hJpsiYield[k]->GetBinContent(bin);
	  double err2 = fabs(1-hJpsiSigExtErr[k]->GetBinContent(bin));
	  double err = sqrt(err1*err1 + err2*err2);
	  hJpsiYield[k]->SetBinContent(bin, hJpsiYield[k]->GetBinContent(bin)*hJpsiYield[k]->GetBinCenter(bin));
	  hJpsiYield[k]->SetBinError(bin, err*hJpsiYield[k]->GetBinContent(bin));
	}
    }

  TF1 *funcYield[nCentBins];
  TH1F *hJpsiPtPos[nCentBins];
  TCanvas *cFit = new TCanvas("FitYield","FitYield",1100,700);
  cFit->Divide(3,2);
  for(int k=0; k<nCentBins; k++)
    {
      double max_pt = 15;
      if(k==2 || k==3) max_pt = 10;
      if(k==4) max_pt = 6;
      funcYield[k] = new TF1(Form("Fit_%s",hJpsiYield[k]->GetName()), "[0]*pow(1+(x+(x*x))/([1]*[2]),-[1])*x", 0,max_pt);
      funcYield[k]->SetParameter(1, 6);
      funcYield[k]->SetParameter(2, 2);
      //if(k==2) funcYield[k]->SetParameter(1, 10);
      hJpsiYield[k]->Fit(funcYield[k], "R0Q");
      cFit->cd(k+1);
      gPad->SetLogy();
      hJpsiYield[k]->SetMarkerStyle(24);
      hJpsiYield[k]->GetYaxis()->SetRangeUser(1e-11, 1e-4);
      hJpsiYield[k]->GetXaxis()->SetRangeUser(0.15, max_pt);
      hJpsiYield[k]->SetTitle(";p_{T} (GeV/c);d^{2}N/dp_{T}d#eta");
      hJpsiYield[k]->Draw("PE");
      funcYield[k]->SetLineStyle(2);
      funcYield[k]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s%%: J/psi yield distribution",cent_Name[k]),0.045);
      t1->Draw();      
	  
      hJpsiPtPos[k] = (TH1F*)hJpsiYield[k]->Clone(Form("hJpsiPtPos_cent%s",cent_Title[k]));
      hJpsiPtPos[k]->Reset();
      for(int bin=1; bin<=hJpsiYield[k]->GetNbinsX(); bin++)
	{
	  double min_pt = hJpsiYield[k]->GetXaxis()->GetBinLowEdge(bin);
	  double max_pt = hJpsiYield[k]->GetXaxis()->GetBinUpEdge(bin);
	  double y_data = hJpsiYield[k]->GetBinContent(bin);
	  double y_fit = funcYield[k]->Integral(min_pt, max_pt)/(max_pt-min_pt);
	  double x_est = funcYield[k]->GetX(y_fit,min_pt,max_pt);
	  //printf("[i] bin %d, yield = %4.4e, func x = %4.4f y = %4.4e y1 = %4.4e\n",bin,y,x_est,funcYield->Eval(x_est),funcInvYield[k][i]->Eval(x_est));
	  hJpsiPtPos[k]->SetBinContent(bin, x_est);
	  hJpsiPtPos[k]->SetBinError(bin, 1e-10);
	  printf("[i] centrality %s%%, bin [%1.0f,%1.0f], y_data = %4.2e, y_fit = %4.2e, pt = %4.2f\n",cent_Name[k],min_pt,max_pt,y_data,y_fit,x_est);
	}
    }
  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sFitJpsiPtPos.pdf",run_type,run_config));
  if(gSaveAN)
    {
      cFit->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_FitJpsiInvYield.pdf"));
    }

  TList *list = new TList;
  TString legName[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      list->Add(hJpsiPtPos[k]);
      legName[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,"PtPos_cent",Form("%s: estimated p_{T} position;p_{T} (GeV/c);p_{T} position",run_type),false,2.0,3.8,true,0,15,kFALSE,kTRUE,legName,true,"",0.15,0.3,0.6,0.85,kTRUE,0.04,0.04); 
  TLine *line = GetLine(0,0,15,15,1);
  line->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiPtPos.pdf",run_type,run_config));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch6_JpsiPtPos.pdf"));
    }

  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  funcYield[k]->Write(Form("func_JpsiYieldVsPt_cent%s",cent_Title[k]),TObject::kOverwrite);
	  hJpsiPtPos[k]->Write("",TObject::kOverwrite);
	}
    }
  

}

//================================================
void compare(const bool savePlot = 1)
{
  // const int nFile = 2;
  // const char *name[nFile] = {"Pico.Run14.AuAu200.jpsi.dtof1.root","Pico.Run14.AuAu200.jpsi.LooseCut.root"};
  // const TString legName[nFile] = {"dtof < 1 ns","Loose cut"};
  // const char *save_name = "dtof1VsLooseCut";

  const int nFile = 2;
  const char *name[nFile] = {"Pico.Run14.AuAu200.jpsi.dtof1.Xsec.root","Pico.Run14.AuAu200.jpsi.dtof1.Xsec.pt1.5.pt1.5.root"};
  const TString legName[nFile] = {"p_{T,1} > 1.5, p_{T,2} > 1 GeV/c","p_{T,1}, p_{T,2} > 1.5 GeV/c"};
  const char *save_name = "1GeVvs1.5GeV";

  TFile *f[nFile];
  TH1F *hYield[nCentBins][nFile];
  for(int i=0; i<nFile; i++)
    {
      f[i] = TFile::Open(Form("Rootfiles/%s",name[i]),"read");
      for(int k=0; k<nCentBins; k++)
	{
	  hYield[k][i] = (TH1F*)f[i]->Get(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
	  hYield[k][i]->SetName(Form("%s_%d",hYield[k][i]->GetName(),i));
	  hYield[k][i]->SetMarkerStyle(21+i*4);
	  hYield[k][i]->SetMarkerColor(1+i);
	  hYield[k][i]->SetLineColor(1+i);
	}
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-10,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      hAuAu->Draw();
      hYield[k][0]->Draw("sames");
      hYield[k][1]->Draw("sames");

      TPaveText *t1 = GetPaveText(0.7,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }

  c->cd(1);
  TLegend *leg = new TLegend(0.18,0.2,0.42,0.48);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int i=0; i<nFile; i++)
    {
      leg->AddEntry(hYield[0][i],legName[i].Data(),"P");
    }
  leg->Draw();

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/CompareXsec_%s.pdf",run_type,save_name));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/CompareXsec_%s.png",run_type,save_name));
    }
}


/*
//===============================================
void shiftDataPoint(const bool savePlot = 0, const bool saveHisto = 0)
{
  //==============================================
  // Cross section vs. pT
  //==============================================
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  char *filename = Form("Run14_AuAu200.JpsiXsec.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s",filename),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s",filename),"read");
 
  TH1F *hJpsiInvYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)fin->Get(Form("Jpsi_InvYieldVsPt_cent%s",cent_Title[k]));
    }

  // add uncertainty from signal extraction to the statistical error
  TFile *ferr = TFile::Open(Form("Rootfiles/%s.Sys.JpsiYield.root",run_type));
  TH1F *hJpsiSigExtErr[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiSigExtErr[k] = (TH1F*)ferr->Get(Form("Sys_signalExt_%s",cent_Title[k]));
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  double err1 = hJpsiInvYield[k]->GetBinError(bin)/hJpsiInvYield[k]->GetBinContent(bin);
	  double err2 = fabs(1-hJpsiSigExtErr[k]->GetBinContent(bin));
	  double err = sqrt(err1*err1 + err2*err2);
	  hJpsiInvYield[k]->SetBinError(bin, err*hJpsiInvYield[k]->GetBinContent(bin));
	}
    }

  // use iterative procedure to obtain bin position 
  const int nIter = 1;
  TH1F *hJpsiInvYield[nCentBins]; 
  TF1 *funcInvYield[nCentBins][nIter];
  double ptPos[nCentBins][nPtBins][nIter+1];
  for(int k=0; k<nCentBins; k++)
    {
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  for(int i=0; i<nIter+1; i++)
	    {
	      ptPos[k][bin-1][i] = hJpsiInvYield[k]->GetBinCenter(bin);
	    }
	}
    }

  double x,y;
  TH1F *hJpsiPtPos[nCentBins][nIter];
  TList *list = new TList;
  TString legName[nIter];
  TF1 *funcYield = new TF1("funcYield", "[0]*pow(1+(x+(x*x))/([1]*[2]),-[1])*x", 0,15);
  TCanvas *cFit = new TCanvas("FitYield","FitYield",1100,700);
  cFit->Divide(3,2);
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<nIter; i++)
	{
	  gJpsiInvYield[k][i] = new TGraphErrors(hJpsiInvYield[k]);
	  for(int ipoint=0; ipoint<gJpsiInvYield[k][i]->GetN(); ipoint++)
	    {
	      x = ptPos[k][ipoint][i];
	      y = hJpsiInvYield[k]->GetBinContent(ipoint+1) * hJpsiInvYield[k]->GetBinCenter(ipoint+1)/x;
	      double x_err = 0;
	      double y_err = hJpsiInvYield[k]->GetBinError(ipoint+1) * hJpsiInvYield[k]->GetBinCenter(ipoint+1)/x;
	      gJpsiInvYield[k][i]->SetPoint(ipoint, x, y);
	      gJpsiInvYield[k][i]->SetPointError(ipoint, x_err, y_err);
	    }
	  double max_pt = 15;
	  if(k==2 || k==3) max_pt = 10;
	  if(k==4) max_pt = 6;
	  gJpsiInvYield[k][i]->SetName(Form("%s_iter%d",hJpsiInvYield[k]->GetName(),i+1));
	  funcInvYield[k][i] = new TF1(Form("Fit_%s",gJpsiInvYield[k][i]->GetName()), "[0]*pow(1+(x+(x*x))/([1]*[2]),-[1])", 0,max_pt);
	  if(k==4) funcInvYield[k][i]->SetRange(1,max_pt);
	  funcInvYield[k][i]->SetParameter(1, 5);
	  funcInvYield[k][i]->SetParameter(2, 1);
	  gJpsiInvYield[k][i]->Fit(funcInvYield[k][i], "IR0Q");
	  cFit->cd(k+1);
	  gPad->SetLogy();
	  gJpsiInvYield[k][i]->SetMarkerStyle(24);
	  gJpsiInvYield[k][i]->GetYaxis()->SetRangeUser(1e-11, 1e-4);
	  gJpsiInvYield[k][i]->GetXaxis()->SetRangeUser(0, max_pt);
	  gJpsiInvYield[k][i]->SetTitle(";p_{T} (GeV/c);d^{2}N/p_{T}dp_{T}d#eta");
	  gJpsiInvYield[k][i]->Draw("AZPE");
	  funcInvYield[k][i]->SetLineStyle(2);
	  funcInvYield[k][i]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%s%%: J/psi invariant yield distribution",cent_Name[k]),0.045);
	  t1->Draw();

	  for(int ipar=0; ipar<funcYield->GetNpar(); ipar++)
	    {
	      funcYield->SetParameter(ipar, funcInvYield[k][i]->GetParameter(ipar));
	      funcYield->SetParError(ipar, funcInvYield[k][i]->GetParError(ipar));
	    }	      
	  
	  hJpsiPtPos[k][i] = (TH1F*)hJpsiInvYield[k]->Clone(Form("hJpsiPtPos_cent%s_iter%d",cent_Title[k],i+1));
	  hJpsiPtPos[k][i]->Reset();
	  for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	    {
	      y = hJpsiInvYield[k]->GetBinContent(bin) * hJpsiInvYield[k]->GetBinCenter(bin);
	      double min_pt = hJpsiInvYield[k]->GetXaxis()->GetBinLowEdge(bin);
	      double max_pt = hJpsiInvYield[k]->GetXaxis()->GetBinUpEdge(bin);
	      double x_est = funcYield->GetX(y,min_pt,max_pt);
	      if(x_est<min_pt) y_est = min_pt;
	      if(x_est>max_pt) y_est = max_pt;
	      //printf("[i] bin %d, yield = %4.4e, func x = %4.4f y = %4.4e y1 = %4.4e\n",bin,y,x_est,funcYield->Eval(x_est),funcInvYield[k][i]->Eval(x_est));

	      ptPos[k][bin-1][i+1] = x_est;
	      hJpsiPtPos[k][i]->SetBinContent(bin, x_est-ptPos[k][bin-1][i]);
	      hJpsiPtPos[k][i]->SetBinError(bin, 1e-10);
	      if(i==nIter-1)
		printf("[i] centrality %s%%, bin [%1.0f,%1.0f], pt = %4.2f\n",cent_Name[k],min_pt,max_pt,x_est);
	    }
	}
      if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sFitJpsiPtPos_cent%s.pdf",run_type,run_config,cent_Title[k]));

      list->Clear();
      for(int i=0; i<nIter; i++)
	{
	  list->Add(hJpsiPtPos[k][i]);
	  legName[i] = Form("Iter %d",i+1);
	}
      c = drawHistos(list,Form("PtPos_cent%s",cent_Title[k]),Form("%s: estimated p_{T} position w.r.t. previous iteration",run_type),false,2.0,3.8,true,-3,1,kFALSE,kTRUE,legName,true,Form("%s%%",cent_Name[k]),0.15,0.3,0.15,0.45,kTRUE,0.04,0.04); 
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiPtPos_cent%s.pdf",run_type,run_config,cent_Title[k]));
    }

}
 */

