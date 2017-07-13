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

  //makeEmbed();
  //xsec();
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


