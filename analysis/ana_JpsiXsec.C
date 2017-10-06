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
  //sfToMb();
  //xsec_Run13();
  //trgSetup();
}

//================================================
void trgSetup(const bool savePlot = 0)
{
  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  // =============================================
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  double mb_events[nCentBins][gNTrgSetup];
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
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hEvtRunAcc->FindBin(runnumber);
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
      printf("[i] # of events for %d: %1.0f\n",j,mb_events[0][j]);
    }
  printf("[i]Check %1.0f =? %1.0f\n",mb_events[0][0],mb_events[0][1]+mb_events[0][2]+mb_events[0][3]+mb_events[0][4]);
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================
  

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
  TList *list = new TList;
  TString legName[4] = {"Run < 15098067","15098068 <= Run < 15106131","Run >= 15106131","Average correction factor"};
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
      list->Add(hAccLoss[i]);
      legName[i] = Form("%s: %2.1f%%",legName[i].Data(),evtCount[i]/evtCount[3]*100);
    }
  list->Add(hAccCorr);
  c = drawHistos(list,"MTD_Acceptance","Correction for MTD acceptance loss",kFALSE,0,30,true,0.8,1.02,kFALSE,kTRUE,legName,kTRUE,"",0.35,0.75,0.18,0.45,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MTD_AccepLoss.pdf",run_type));
    }
  // =============================================
  //
  //

  // Jpsi efficiency
  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.JpsiEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiEff[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  if(j==0) hJpsiEff[k][j] = (TH1F*)fEff->Get(Form("JpsiPtEff_cent0020_corr"));
	  else     hJpsiEff[k][j] = (TH1F*)fEff->Get(Form("JpsiPtEff_Final_cent0020%s_w_rebin",gTrgSetupName[j]));
	  hJpsiEff[k][j]->SetName(Form("JpsiPtEff_cent%s%s",cent_Title[k],gTrgSetupTitle[j]));
	  TH1F *hCorr = (TH1F*)fEff->Get(Form("hJpsiEffCorr_cent%s%s",cent_Title[k],gTrgSetupTitle[j]));
	  hJpsiEff[k][j]->Multiply(hCorr);
	}
    }
  

  // Jpsi raw counts
  char * yieldName = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiCounts[k][j] = (TH1F*)fYield->Get(Form("Jpsi_BinCountYield_cent%s_weight%s",cent_Title[k],gTrgSetupName[j]));
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
	      hJpsiInvYield[k][j]->SetBinContent(bin,hJpsiInvYield[k][j]->GetBinContent(bin)/bin_width/bin_center);
	      hJpsiInvYield[k][j]->SetBinError(bin,hJpsiInvYield[k][j]->GetBinError(bin)/bin_width/bin_center);
	    }
	  hJpsiInvYield[k][j]->Scale(1./mb_events[k][j]); // N_evt
	  hJpsiInvYield[k][j]->Scale(1./(2*pi)); // 2pi
	  hJpsiInvYield[k][j]->Scale(1./1.6); // dy
	  hJpsiInvYield[k][j]->SetMarkerStyle(21);
	  hJpsiInvYield[k][j]->SetMarkerColor(color[j]);
	  hJpsiInvYield[k][j]->SetLineColor(color[j]);
	  hJpsiInvYield[k][j]->SetMarkerSize(1.5);

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
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",15,0,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-11,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  TLegend *leg = new TLegend(0.45,0.7,0.8,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05); 
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      hAuAu->DrawCopy();
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[k][j]->Draw("sames");
	  if(k==0) leg->AddEntry(hJpsiInvYield[k][j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
	}
      TPaveText *t1 = GetPaveText(0.2,0.3,0.2,0.3,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_TrgSetupComp.pdf",run_type,run_config));

  TCanvas *c = new TCanvas("Ratio_JetXsec","Ratio_JetXsec",1100,700);
  c->Divide(2,2);
  hAuAu->SetYTitle("Ratio to combined");
  hAuAu->GetYaxis()->SetRangeUser(0,3.5);
  TLegend *leg = new TLegend(0.45,0.7,0.8,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05); 
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      hAuAu->DrawCopy();
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hJpsiInvRatio[k][j]->Draw("sames");
	  if(k==0) leg->AddEntry(hJpsiInvRatio[k][j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
	}
      TPaveText *t1 = GetPaveText(0.2,0.3,0.8,0.9,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_TrgSetupRatio.pdf",run_type,run_config));

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
  TList *list = new TList;
  TString legName[4] = {"Run < 15098067","15098068 <= Run < 15106131","Run >= 15106131","Average correction factor"};
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
      list->Add(hAccLoss[i]);
      legName[i] = Form("%s: %2.1f%%",legName[i].Data(),evtCount[i]/evtCount[3]*100);
    }
  list->Add(hAccCorr);
  c = drawHistos(list,"MTD_Acceptance","Correction for MTD acceptance loss",kFALSE,0,30,true,0.8,1.02,kFALSE,kTRUE,legName,kTRUE,"",0.35,0.75,0.18,0.45,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MTD_AccepLoss.pdf",run_type));
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
	  hJpsiInvYield[k]->SetBinContent(bin,hJpsiInvYield[k]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[k]->SetBinError(bin,hJpsiInvYield[k]->GetBinError(bin)/bin_width/bin_center);
	}

      hJpsiInvYield[k]->Scale(1./mb_events[k]); // N_evt
      hJpsiInvYield[k]->Scale(1./(2*pi)); // 2pi
      //hJpsiInvYield[k]->Scale(1./1.6); // dy
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
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",15,0,15);
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
	  gAuAuYield_sQM[k-1]->Draw("sames PE");
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
  leg->AddEntry(gAuAuYield_sQM[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5 (sQM2016)","P");
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

	  gAuAuYield_sQM[k]->GetPoint(bin-1, x, y);
	  hRatio_sQM->SetBinContent(bin,y/scale);
	  hRatio_sQM->SetBinError(bin,gAuAuYield_sQM[k]->GetErrorYhigh(bin-1)/scale);
	  
	}
      c->cd(k+1);
      gPad->SetGridy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
 
      // final
      hRatio->SetTitle(";p_{T} (GeV/c);Ratio to TBW");
      hRatio->GetXaxis()->SetRangeUser(0,10);
      hRatio->GetYaxis()->SetRangeUser(0,2);
      ScaleHistoTitle(hRatio,0.06,1,0.05,0.06,1,0.05,62);
      hRatio->Draw();

      // sQM2016
      hRatio_sQM->SetMarkerStyle(29);
      hRatio_sQM->SetMarkerColor(kGreen+2);
      hRatio_sQM->SetLineColor(kGreen+2);
      hRatio_sQM->SetMarkerSize(2);
      hRatio_sQM->Draw("samesP");

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
  leg->AddEntry(gAuAuYield_sQM[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5 (sQM2016)","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_RatioToTBW.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_RatioToTBW.png",run_type,run_config));
    }

  // Raa distribution
  const double ppInelastic = 42.; // mb
  double ncoll[nCentBins] = {291.9, 766.47, 290.87, 91.33, 21.57};
  TH1F *hJpsipp = (TH1F*)fpub->Get("hPPJpsiFinal");
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
      if(k==2 || k==3) hJpsiRaa[k]->GetXaxis()->SetRangeUser(0,10);
      if(k==4) hJpsiRaa[k]->GetXaxis()->SetRangeUser(0,6);
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
  TH1F *hppNpart = (TH1F*)fpub->Get("pp200_Jpsi_Integrated");
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
      TCanvas *c = draw1D(hJpsiRaaVsCent[i],Form("%s: J/#Psi R_{AA} vs. centrality for p_{T} > %1.0f GeV/c;;R_{AA}",run_type,ptBins_low_npart[i]));
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



//================================================
void xsec_Run13(const bool savePlot = 0, const bool saveHisto = 1)
{
  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %d\n",hStat->GetBinContent(1));
  printf("all di-muon events: %d\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %d\n",hStat->GetBinContent(10));
  const double dimuon_events = hStat->GetBinContent(3);

  // Luminosity
  const double sample_Lum = 28.272; // pb-1
  const double sample_Nev = 118476807;
  const double nsd_xsec = 34 * 1e6; // nb
  TFile *fLumi = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  TH1F *hMBevents = (TH1F*)fLumi->Get("hMBevents");
  const double process_MTD = hMBevents->GetBinContent(1);
  const double process_VPD = hMBevents->GetBinContent(2);

  /*
  TFile *fVtxEff = TFile::Open("Rootfiles/Run13.pp500.VPDMB.VtxCutEff.root","read");
  TH1F *hVtxEff = (TH1F*)fVtxEff->Get("hVtxEff");
  const double mb_vtx_eff = hVtxEff->GetBinContent(1);
  */
  if(run_config=="VtxCut.") const double ana_Lum = dimuon_events * process_VPD/process_MTD * 1./nsd_xsec; // nb-1
  else                      const double ana_Lum = dimuon_events * sample_Lum / sample_Nev * 1e3;

  printf("+++++++++++++++++++++++++++++++++\n");
  printf("# of dimuon events is: %d\n",dimuon_events);
  printf("# of VPD events is: %e\n",dimuon_events * process_VPD/process_MTD);
  printf("Effective luminosity %5.3f nb-1\n",ana_Lum);
  printf("+++++++++++++++++++++++++++++++++\n");

  // trigger bias
  TFile *fTrgBias = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  TH1F *hTrgBias = (TH1F*)fTrgBias->Get("VPDMB_TrigBias");

  // Jpsi efficiency
  char *embedEffName = Form("Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  char *trigEffName = Form("Run13.pp500.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  printf("Embed   eff: %s\n",embedEffName);
  printf("Trigger eff: %s\n",trigEffName);
  TFile *fEmbedEff = TFile::Open(Form("Rootfiles/%s",embedEffName),"read");
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/%s",trigEffName),"read");
  TH1F *hJpsiEffEmbed[nCentBins];
  TH1F *hJpsiEffTrig[nCentBins];
  TH1F *hJpsiRespEff[nCentBins];
  TH1F *hJpsiSmearEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEffEmbed[k] = (TH1F*)fEmbedEff->Get(Form("MTDreco_Jpsi_pT_%s_WeightPt_Eff_rebin",cent_Title[k]));
      hJpsiEffTrig[k] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s_rebin",cent_Title[k]));
      hJpsiRespEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiRespEff_cent%s_rebin",cent_Title[k]));
      hJpsiSmearEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiSmearEff_cent%s_rebin",cent_Title[k]));
    }

  // Jpsi raw counts
  char * yieldName = Form("Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiCounts[k] = (TH1F*)fYield->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[k]));
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)hJpsiCounts[k]->Clone(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiInvYield[k]->Divide(hJpsiEffEmbed[k]);
      //hJpsiInvYield[k]->Divide(hJpsiEffTrig[k]); // 100% trigger efficiency
      hJpsiInvYield[k]->Divide(hJpsiRespEff[k]);
      hJpsiInvYield[k]->Divide(hJpsiSmearEff[k]);
      if(run_config=="VtxCut.")  hJpsiInvYield[k]->Multiply(hTrgBias);
      cout << hJpsiEffEmbed[k]->GetBinContent(1) << "  " << hJpsiEffTrig[k]->GetBinContent(1) << endl;
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  double bin_width = hJpsiInvYield[k]->GetBinWidth(bin); // dpT
	  double bin_center = hJpsiInvYield[k]->GetBinCenter(bin); // pT 
	  hJpsiInvYield[k]->SetBinContent(bin,hJpsiInvYield[k]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[k]->SetBinError(bin,hJpsiInvYield[k]->GetBinError(bin)/bin_width/bin_center);
	}

      hJpsiInvYield[k]->Scale(1./ana_Lum); // N_evt
      hJpsiInvYield[k]->Scale(1./(2*pi)); // 2pi
      hJpsiInvYield[k]->Scale(1./1.6); // dy
      hJpsiInvYield[k]->SetMarkerStyle(21);
      hJpsiInvYield[k]->SetMarkerColor(1);
      hJpsiInvYield[k]->SetLineColor(1);
      hJpsiInvYield[k]->SetMarkerSize(1.5);
    }

  TCanvas *c = new TCanvas("c2","c2", 700, 700);
  SetPadMargin(gPad,0.12, 0.14, 0.03,0.03);
  gPad->SetLogy();
  TH1F *h = new TH1F("h2",";p_{T} (GeV/c);B#times1/(2#pip_{T})#timesd^{2}#sigma/(dp_{T}dy)   (nb/GeV/c)^{2}",1000,0,25);
  ScaleHistoTitle(h,0.045,1,0.035,0.045,1.4,0.035,62);
  h->GetYaxis()->SetRangeUser(1e-6,100);
  h->GetXaxis()->SetRangeUser(0,22.5);
  h->GetYaxis()->CenterTitle(1);
  h->DrawCopy();

  leg = new TLegend(0.45,0.65,0.65,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader("p+p @ 500 GeV");

  // CGC
  TFile *fCGC = TFile::Open("Rootfiles/Published/Model/pp500_Jpsi_xsec/CGC.root","read");
  TGraphAsymmErrors *gCgcLowPt = (TGraphAsymmErrors*)fCGC->Get("CGC_lowPt");
  TGraphAsymmErrors *gCgcHighPt = (TGraphAsymmErrors*)fCGC->Get("CGC_highPt");
  gCgcLowPt->Draw("sames E3");
  gCgcHighPt->Draw("sames E3");
  leg->AddEntry(gCgcLowPt,"CGC+NRQCD","F");
  leg->AddEntry(gCgcHighPt,"NLO+NRQCD","F");

  /*      
  //TFile *fFit = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
  TFile *fFit = TFile::Open("Rootfiles/pt500GeVfit_new.root","read");
  TF1 *funcJpsi  = (TF1*)fFit->Get("ffpt");
  TF1 *funcJpsi1 = (TF1*)fFit->Get("ffpt1");
  TF1 *funcJpsi2 = (TF1*)fFit->Get("ffpt2");
  const int nPoints = 1000;
  TGraphAsymmErrors *gGJpsiSys = new TGraphAsymmErrors(nPoints);
  for(int i=0; i<nPoints; i++)
    {
      double x = 20./nPoints * (i+1);
      double y = funcJpsi->Eval(x);
      gGJpsiSys->SetPoint(i,x,y);
      gGJpsiSys->SetPointError(i, 0, 0, funcJpsi1->Eval(x)-y, y-funcJpsi2->Eval(x));
    }
  gGJpsiSys->SetFillStyle(1001);
  gGJpsiSys->SetLineColor(kGray+1);
  gGJpsiSys->SetFillColor(kGray+1);  
  gGJpsiSys->Draw("sames E3");
  funcJpsi->Draw("sames");
  leg->AddEntry(gGJpsiSys,"Global fit","F");
  */

  // Run11 e+e-
  TFile *fee = TFile::Open(Form("Rootfiles/2015HP/sptrum.root"),"read");
  TGraphErrors *gData = (TGraphErrors*)fee->Get("Jpsi_pp500");
  gData->Draw("sames PEZ");
  for(int i=0; i<19; i++)
    {
      TBox *box = (TBox*)fee->Get(Form("sys_uncert_%d",i));
      box->Draw();
    }
  leg->AddEntry(gData,"Run11: J/#psi#rightarrowe^{+}e^{-}, |y|<1","P");

  hJpsiInvYield[0]->Draw("samesP");
  leg->AddEntry(hJpsiInvYield[0],"Run13: J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");

  // Combined fit
  TCanvas *cFit = new TCanvas("cFit","cFit", 700, 700);
  SetPadMargin(gPad,0.12, 0.14, 0.03,0.03);
  gPad->SetLogy();
  h->DrawCopy();
  TH1F *hmm = (TH1F*)hJpsiInvYield[0]->Clone("hmm");
  hmm->Draw("sames");
  TH1F *hee = gData->GetHistogram();
  
  double x,y;
  for(int ipoint=0; ipoint<=gData->GetN(); ipoint++)
    {
      gData->GetPoint(ipoint, x, y);
      int bin = hee->FindFixBin(x);
      hee->SetBinContent(bin, y);
      hee->SetBinError(bin, gData->GetErrorY(ipoint));
      //cout << bin << "  " << x << "  " << y << endl;
    }
  hee->SetMarkerStyle(20);
  hee->SetMarkerColor(2);
  hee->SetLineColor(2);
  hee->Draw("sames");

  ROOT::Fit::BinData data; 
  //ROOT::Fit::FillData(data, hmm); 
  ROOT::Fit::FillData(data, hee); 

  cout << "data size is " << data.Size() << endl;

  TF1 * funcFit = new TF1("funcFit","[0]*exp([1]+[2]*x+[3]*x*x)",3.5,14);
  ROOT::Math::WrappedTF1 wf(*funcFit);
  ROOT::Fit::Fitter fitter;
  fitter.SetFunction(wf);
  fitter.Fit(data);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  funcFit->Draw("sames");

  c->cd();
  //funcFit->Draw("sames");
  //leg->AddEntry(funcFit,"Exp fit to Run11","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_Compare.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_Compare.png",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_Compare.eps",run_type,run_config));
    }
  //return;


  // ratio to Global fit
  TCanvas *c = new TCanvas("pp500_Jpsi_ratio","pp500_Jpsi_ratio",800,600);
  //gPad->SetLogy();
  h->SetTitle(";p_{T} (GeV/c);Run13/Run11Fit");
  ScaleHistoTitle(h,0.045,1,0.035,0.045,1,0.035,62);
  h->GetYaxis()->SetRangeUser(0.001,2);
  h->GetXaxis()->SetRangeUser(4,18);
  h->GetYaxis()->CenterTitle(1);
  h->Draw();

  /*
  TGraphAsymmErrors *gRJpsiSys = new TGraphAsymmErrors(nPoints);
  for(int i=0; i<nPoints; i++)
    {
      double x = 20./nPoints * (i+1);
      double y = funcJpsi->Eval(x);
      gRJpsiSys->SetPoint(i,x,1);
      gRJpsiSys->SetPointError(i, 0, 0, funcJpsi1->Eval(x)/y-1, 1-funcJpsi2->Eval(x)/y);
    }
  gRJpsiSys->SetFillStyle(1001);
  gRJpsiSys->SetLineColor(kGray);
  gRJpsiSys->SetFillColor(kGray);  
  gRJpsiSys->Draw("sames E3");
  */

  TF1 *funcJpsi = (TF1*)funcFit->Clone("funcJpsi");
  funcJpsi->SetRange(0,20);
  funcJpsi->SetNpx(1000);
  TH1F *hFitJpsiPt = (TH1F*)funcJpsi->GetHistogram();
  hFitJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield"));
  for(int bin=1; bin<=hFitJpsiPt->GetNbinsX(); bin++)
    {
      hFitJpsiPt->SetBinContent(bin,hFitJpsiPt->GetBinContent(bin));
    }

  TH1F *hRatio = (TH1F*)hJpsiInvYield[0]->Clone(Form("%s_ratio",hJpsiInvYield[0]->GetName()));
  for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
    {
      int start_bin = hFitJpsiPt->FindBin(hRatio->GetXaxis()->GetBinLowEdge(bin)+1e-6);
      int end_bin   = hFitJpsiPt->FindBin(hRatio->GetXaxis()->GetBinUpEdge(bin)-1e-6);
      double scale = 0;
      for(int ibin=start_bin; ibin<=end_bin; ibin++)
	{
	  scale += hFitJpsiPt->GetBinContent(ibin) * hFitJpsiPt->GetBinCenter(ibin) * hFitJpsiPt->GetBinWidth(ibin);
	}
      scale = scale / hRatio->GetBinCenter(bin) / hRatio->GetBinWidth(bin);
      cout << scale << " -> " << hRatio->GetBinContent(bin) << endl;
      hRatio->SetBinContent(bin,hRatio->GetBinContent(bin)/scale);
      hRatio->SetBinError(bin,hRatio->GetBinError(bin)/scale);
    }
  hRatio->Draw("sames");

  TGraphErrors *gRatio = (TGraphErrors*)gData->Clone("gRatio");
  double x,y;
  double *xarr = gData->GetX();
  double *yarr = gData->GetY();
  for(int ip=0; ip<gData->GetN(); ip++)
    {
      double ex = gData->GetErrorX(ip);
      double ey = gData->GetErrorY(ip);
      double scale = funcJpsi->Eval(xarr[ip]);
      //cout << xarr[ip] << "  " << yarr[ip]/scale << "  " << scale << "  " << ey << endl;
      gRatio->SetPoint(ip,xarr[ip],yarr[ip]/scale);
      gRatio->SetPointError(ip,ex,ey/scale);
     
    }
  gRatio->Draw("sames PEZ");

  TLine *line = GetLine(4,1,16,1,1);
  line->Draw();

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToCFit.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToCFit.png",run_type,run_config));
    }

  if(saveHisto)
    {
      char *outname = "";
      if(year==2013) outname = Form("Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut);
      if(year==2014) outname = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut);
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outname),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiInvYield[k]->Write();
	}
    }

}
