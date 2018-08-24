const int nPtBins        = 9;
const double xPtBins[10] = {0,1,2,3,4,5,6,8,10,15}
const int nCentBins      = nCentBins_pt; 
const int* centBinsLow   = centBins_low_pt;
const int* centBinsHigh  = centBins_high_pt;
const char** cent_Name    = cent_Name_pt;
const char** cent_Title   = cent_Title_pt;
const int icent = 4;
const char* gPairName[2] = {"UL", "LS"};
const int gPairType = 1;

//================================================
void check_Cuts()
{
  gStyle->SetOptStat(0);

  //makeHistos();
  xsec();
}

//================================================
void xsec(const int savePlot = 0)
{
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  double mb_events[gNTrgSetup];
  double mtd_events[gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      mb_events[j] = 0;
      mtd_events[j] = 0;
    }

  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  printf("[i] Process %s\n",fdata->GetName());
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");

  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[icent]));

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
	  int lumiBin = hNeventsTake->FindFixBin(runnumber);
	  double nEventsTaken = hNeventsTake->GetBinContent(lumiBin);
	  if(nEventsTaken==0) 
	    {
	      printf("[w] check run %1.0f\n",runnumber);
	      continue;
	    }
	  double nEventsRun = hEvtRun->GetBinContent(bin);
	  double rf = hRF->GetBinContent(hRF->FindFixBin(runnumber));
	  if(rf==0)
	    {
	      printf("[w] rf = 0 for run %1.0f\n",runnumber);
	      rf = 0.49;
	    }

	  double eq_mb = hEqMbEvents->GetBinContent(hEqMbEvents->FindFixBin(runnumber));
	  //mb_events[j] += nEventsRun/rf/nEventsTaken * eq_mb;
	  //mb_events[0] += nEventsRun/rf/nEventsTaken * eq_mb;
	  mb_events[j] += eq_mb;
	  mb_events[0] += eq_mb;

	  mtd_events[j] += nEventsRun;
	  mtd_events[0] += nEventsRun;
	}
    }
  for(int j=0; j<gNTrgSetup; j++)
    {
      printf("[i] MB events for %s%s is: %1.0f (%4.2f%%)\n",run_type,gTrgSetupTitle[j],mb_events[j],mb_events[j]/mb_events[0]*100);
    }
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================


  // =============================================
  // MTD acceptance loss
  const int nRunRange = 8;
  double evtCount[nRunRange+1][gNTrgSetup];
  for(int i=0; i<nRunRange+1; i++) 
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  evtCount[i][j] = 0;
	}
    }
  int runRange[nRunRange+1] = {15074104, 15077035, 15078021, 15098066, 15099002, 15106130, 15131038, 15132019, 15167014};
  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
    {
      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
      double run = hEvtRunAcc->GetBinCenter(bin);
      double evt = hEvtRunAcc->GetBinContent(bin);
      for(int i=0; i<nRunRange; i++)
	{
	  if(run>runRange[i] && run<=runRange[i+1]) evtCount[i][0] += evt;
	}
      evtCount[nRunRange][0] += evt;
    }    

  for(int j=1; j<gNTrgSetup; j++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type,trgSetupName[j-1]));
      int runnumber;
      while(fruns >> runnumber)
	{
	  int bin = hEvtRunAcc->FindBin(runnumber);
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	  double evt = hEvtRunAcc->GetBinContent(bin);
	  for(int i=0; i<nRunRange; i++)
	    {
	      if(runnumber>runRange[i] && runnumber<=runRange[i+1]) evtCount[i][j] += evt;
	    }
	  evtCount[nRunRange][j] += evt;
	}
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
  int colors[10] = {1, 2, 4, 6, 1, 2, 4, 6, 1, 2};
  TH1F *hAccCorr[gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      hAccCorr[j] = (TH1F*)hAccLoss[0]->Clone(Form("hAccCorr%s",gTrgSetupTitle[j]));
      hAccCorr[j]->Reset();
      for(int i=0; i<nRunRange; i++)
	{
	  hAccCorr[j]->Add(hAccLoss[i],evtCount[i][j]/evtCount[nRunRange][j]);
	  if(j==0)
	    {
	      list->Add(hAccLoss[i]);
	      hAccLoss[i]->SetMarkerStyle(20+i);
	      hAccLoss[i]->SetMarkerColor(colors[i]);
	      hAccLoss[i]->SetLineColor(colors[i]);
	      legName[i] = Form("Run (%d,%d]: %2.1f%%",runRange[i],runRange[i+1],evtCount[i][j]/evtCount[nRunRange][j]*100);
	    }
	}
      if(j==0)
	{
	  list->Add(hAccCorr[j]);
	}
    }
  c = drawHistos(list,"MTD_Acceptance","Correction for MTD acceptance loss",kFALSE,0,30,true,0.5,1.02,kFALSE,kTRUE,legName,kTRUE,"",0.35,0.75,0.15,0.45,kTRUE,0.04,0.04,false,1,false,false);

  list->Clear();
  const TString legName2[5] = {"Inclusive","prod","prod_low","prod_mid","prod_high"};
  for(int j=0; j<gNTrgSetup; j++)
    {
      list->Add(hAccCorr[j]);
    }
  c = drawHistos(list,"MTDAcc_Lumi","Correction for MTD acceptance loss",kFALSE,0,30,true,0.5,1.02,kFALSE,kTRUE,legName2,kTRUE,"",0.35,0.75,0.15,0.45,kTRUE,0.04,0.04);
  // =============================================
  //
  //

  // Jpsi efficiency
  // take into account efficiency vs. luminosity
  // trigger efficiency; muonPid efficiency
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.Sys.MtdTrigEff.root"),"read");
  TFile *fTrigEffLS = TFile::Open(Form("Rootfiles/Run14_AuAu200.StudyLumiDep.root"),"read");
  TH1F *hJpsiTrigEff[5];
  TH1F *hJpsiTrigEffVsCent[5];
  TH1F *hJpsiTrigEffLS[5];
  TH1F *hJpsiTrigEffVsCentLS[5];
  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiTrigEff[j]       = (TH1F*)fTrigEff->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      if(j>0) hJpsiTrigEff[j]->Divide(hJpsiTrigEff[0]);

      hJpsiTrigEffVsCent[j] = (TH1F*)fTrigEff->Get(Form("Run14_AuAu200_JpsiEffVsCent_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      //printf("[i] Jpsi trigger efficiency for %s: %4.3f\n",gTrgSetupTitle[j], hJpsiTrigEffVsCent[j]->GetBinContent(1));

      hJpsiTrigEffLS[j]       = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      hJpsiTrigEffLS[j]->SetName(Form("%s_LS",hJpsiTrigEffLS[j]->GetName()));
      if(j>0) hJpsiTrigEffLS[j]->Divide(hJpsiTrigEffLS[0]);

      hJpsiTrigEffVsCentLS[j] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsCent_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      hJpsiTrigEffVsCentLS[j]->SetName(Form("%s_LS",hJpsiTrigEffVsCentLS[j]->GetName()));
      if(j>0) printf("[i] LS trigger efficiency for %s: %4.3f\n",gTrgSetupTitle[j],hJpsiTrigEffVsCentLS[j]->GetBinContent(1));      
    }

  const double muonPidScale[gNTrgSetup] = {1, 1.008, 1.006, 1.004, 0.99};
  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",pt1_cut,pt2_cut,run_config),"read");
  TH1F *hJpsiEff[gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiEff[j] = (TH1F*)fEff->Get(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[icent],gTrgSetupTitle[j]));
      if(j>0) hJpsiEff[j]->Multiply(hJpsiTrigEffLS[j]);
      hJpsiEff[j]->Scale(muonPidScale[j]);
      hJpsiEff[j]->Multiply(hAccCorr[j]);
    }

  // inclusive efficiency
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  TH1F *hModel = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hWeight = (TH1F*)hJpsiEff[0]->Clone("hWeight");
  hWeight->Reset("AC");
  for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
    {
      double low_pt  = hWeight->GetXaxis()->GetBinLowEdge(bin);
      double high_pt = hWeight->GetXaxis()->GetBinUpEdge(bin);
      double jpsi_yield = hModel->Integral(hModel->FindFixBin(low_pt+1e-5),hModel->FindFixBin(high_pt-1e-5));
      hWeight->SetBinContent(bin, jpsi_yield);
    }
  hWeight->Scale(1./hWeight->Integral());
  double total_eff[gNTrgSetup];
  for(int j=1; j<gNTrgSetup; j++)
    {
      total_eff[j] = 0;
      for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
	{
	  total_eff[j] += hWeight->GetBinContent(bin) * hJpsiEff[j]->GetBinContent(bin);
	}
      printf("[i] Total efficiency for %s%s is %4.3f%%\n",run_type,gTrgSetupTitle[j],total_eff[j]*100);
    }

  TCanvas *c = new TCanvas("CompEff", "CompEff", 800, 600);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  for(int j=1; j<gNTrgSetup; j++)
    {
      TH1F *hRatio = (TH1F*)hJpsiEff[j]->Clone(Form("%s_clone",hJpsiEff[j]->GetName()));
      hRatio->Divide(hJpsiEff[2]);
      hRatio->SetMarkerStyle(20+j);
      hRatio->SetMarkerColor(color[j]);
      hRatio->SetLineColor(color[j]);
      hRatio->GetYaxis()->SetRangeUser(0.1,1.1);
      if(j==1) hRatio->Draw();
      else     hRatio->Draw("sames");
      leg->AddEntry(hRatio, gTrgSetupTitle[j], "P");
    }
  leg->Draw();

  // Jpsi raw counts
  TH1F *hJpsiCounts[gNTrgSetup];
  TFile *fYield = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiCounts[j] = (TH1F*)fYield->Get(Form("%s_PairPt_%s_cent%s%s",run_config,gPairName[gPairType],cent_Title[icent],gTrgSetupTitle[j]));
      printf("[i] Total counts of %s for %s%s: %1.0f%\n",gPairName[gPairType],run_type,gTrgSetupTitle[j],hJpsiCounts[j]->GetEntries());
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[gNTrgSetup];
  TH1F *hJpsiInvRatio[gNTrgSetup];
  for(int j=1; j<gNTrgSetup; j++)
    {
      hJpsiInvYield[j] = (TH1F*)hJpsiCounts[j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[icent],gTrgSetupTitle[j]));
      hJpsiInvYield[j]->Divide(hJpsiEff[j]);
      for(int bin=1; bin<=hJpsiInvYield[j]->GetNbinsX(); bin++)
	{
	  double bin_width = hJpsiInvYield[j]->GetBinWidth(bin); // dpT
	  double bin_center = hJpsiInvYield[j]->GetBinCenter(bin); // pT 
	  hJpsiInvYield[j]->SetBinContent(bin,hJpsiInvYield[j]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[j]->SetBinError(bin,hJpsiInvYield[j]->GetBinError(bin)/bin_width/bin_center);
	}
      hJpsiInvYield[j]->Scale(1./(2*pi)); // 2pi
      if(j==1) hJpsiInvYield[0] = (TH1F*)hJpsiInvYield[j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[icent],gTrgSetupTitle[0]));
      else     hJpsiInvYield[0]->Add((TH1F*)hJpsiInvYield[j]);
    }

  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiInvYield[j]->Scale(1./mb_events[j]); // N_evt
      hJpsiInvYield[j]->SetMarkerStyle(21);
      hJpsiInvYield[j]->SetMarkerColor(color[j]);
      hJpsiInvYield[j]->SetLineColor(color[j]);
      hJpsiInvYield[j]->SetMarkerSize(1.2);

      hJpsiInvRatio[j] = (TH1F*)hJpsiInvYield[j]->Clone(Form("Jpsi_InvYieldRatio_cent%s%s",cent_Title[icent],gTrgSetupTitle[j]));
      hJpsiInvRatio[j]->Divide(hJpsiInvYield[0]);
      for(int bin=1; bin<=hJpsiInvRatio[j]->GetNbinsX(); bin++)
	{
	  double err = hJpsiInvYield[j]->GetBinError(bin)/hJpsiInvYield[j]->GetBinContent(bin);
	  hJpsiInvRatio[j]->SetBinError(bin,hJpsiInvRatio[j]->GetBinContent(bin)*err);
	}
    }

  TCanvas *c = new TCanvas("Comp_JetXsec","Comp_JetXsec",1100,500);
  c->Divide(2,1);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",15,0,15);
  hAuAu->GetYaxis()->SetRangeUser(1e-11,1e-2);
  ScaleHistoTitle(hAuAu,0.045,1,0.04,0.045,1.4,0.04,62);
  TLegend *leg = new TLegend(0.2,0.15,0.45,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  leg->SetHeader(Form("%s%%",cent_Name[icent]));
  c->cd(1);
  gPad->SetLogy();
  SetPadMargin(gPad,0.13,0.13,0.05,0.1);
  hAuAu->GetXaxis()->SetRangeUser(0,10);
  hAuAu->DrawCopy();
  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiInvYield[j]->Draw("sames");
      leg->AddEntry(hJpsiInvYield[j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
    }
  leg->Draw();
  TPaveText *t1 = GetTitleText(Form("%s pairs",gPairName[gPairType]),0.05);
  t1->Draw();

  c->cd(2);
  hAuAu->SetYTitle("Ratio to inclusive");
  hAuAu->GetYaxis()->SetRangeUser(0.8,1.2);
  TLegend *leg = new TLegend(0.1,0.6,0.4,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05); 
  SetPadMargin(gPad,0.13,0.13,0.05,0.1);
  hAuAu->DrawCopy();
  for(int j=0; j<gNTrgSetup; j++)
    {
      hJpsiInvRatio[j]->Draw("sames");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/%s%sYieldRatio.pdf",run_type,run_config,gPairName[gPairType]));
}

//================================================
void makeHistos(const int saveHisto = 1)
{
  // data
  TFile *fin = TFile::Open(Form("output/Run14_AuAu200.jpsi.%sroot",run_config), "read");
  THnSparseF *hnInvMass[2] = {0x0};
  hnInvMass[0] = (THnSparseF*)fin->Get("mhJpsiInfo_di_mu");
  hnInvMass[1] = (THnSparseF*)fin->Get("mhBkgLSPos_di_mu");
  hnInvMass[1]->Add((THnSparseF*)fin->Get("mhBkgLSNeg_di_mu"));

  TH1F *hPairPt[5][2] = {0x0};
  for(Int_t j=0; j<2; j++) // pair type
    { 
      hnInvMass[j]->GetAxis(0)->SetRangeUser(3.0, 3.2); // mass cut
      hnInvMass[j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);
      hnInvMass[j]->GetAxis(4)->SetRange(centBinsLow[icent],centBinsHigh[icent]);

      for(int t=0; t<gNTrgSetup; t++) // trigger setup
	{
	  if(t>0) hnInvMass[j]->GetAxis(5)->SetRange(t,t);
	  TH1F *htmp = (TH1F*)hnInvMass[j]->Projection(1);
	  htmp->SetName(Form("htmp_%d%d",j,t));
	  htmp->Sumw2();
	  hPairPt[t][j] = (TH1F*)htmp->Rebin(nPtBins, Form("%s_PairPt_%s_cent%s%s",run_config,gPairName[j],cent_Title[icent],gTrgSetupTitle[t]), xPtBins);
	  hPairPt[t][j]->SetTitle();
	  hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(0)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(2)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(0)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(0)->SetRange(0,-1);
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root", "update");
      for(Int_t j=0; j<2; j++) // pair type
	{ 
	  for(int t=0; t<gNTrgSetup; t++) // trigger setup
	    {
	      hPairPt[t][j]->Write("", TObject::kOverwrite);
	    }
	}
    }
}
