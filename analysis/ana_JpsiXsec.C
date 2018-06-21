const int year = YEAR;
#include <fstream>      // std::ofstream

//================================================
void ana_JpsiXsec()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);


  //xsec_Run14();
  //compare();
  trgSetup();
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
  TH1F *hPSdimuon = (TH1F*)fLumi->Get("hPreScale_dimuon");
  TH1F *hLTdimuon = (TH1F*)fLumi->Get("hLiveTime_dimuon");

  // TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_BHT2-VPDMB-30");
  // TH1F *hPSdimuon = (TH1F*)fLumi->Get("hPreScale_BHT2-VPDMB-30");
  // TH1F *hLTdimuon = (TH1F*)fLumi->Get("hLiveTime_BHT2-VPDMB-30");

  TH1F *hPSmb = (TH1F*)fLumi->Get("hPreScale_VPD-ZDC-novtx-mon");
  TH1F *hLTmb = (TH1F*)fLumi->Get("hLiveTime_VPD-ZDC-novtx-mon");
  TH1F *hNeventsTakeMb = (TH1F*)fLumi->Get("hNevents_VPD-ZDC-novtx-mon");
  TH1F *hEqMbEvents[nCentBins];

  int nRun = 0, nGoodRun = 0;
  for(int bin=1; bin<=hEvtRun->GetNbinsX(); bin++)
    {
      if(hEvtRun->GetBinContent(bin)>0) nRun++;
      if(hEvtRunAcc->GetBinContent(bin)>0) nGoodRun++;
    }
  cout << nRun << "  " << nGoodRun << endl;
  
  for(int i=0; i<nCentBins; i++)
    {
      hEqMbEvents[i] = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
      mb_events[i][0] = 0;
      for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
	{
	  double run = hEvtRunAcc->GetBinCenter(bin);
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
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
  TH1F *hPSdimuonInTrig[4];
  TH1F *hLTdimuonInTrig[4];
  TH1F *hNdimuonInTrig[4];
  TH1F *hPSmbInTrig[4];
  TH1F *hLTmbInTrig[4];
  TH1F *hNmbInTrig[4];
  TH1F *hEqvRatioInTrig[4];

  long int dimuonevent[gNTrgSetup];
  long int dimuoneventacc[gNTrgSetup];
  dimuonevent[0] = 0;
  dimuoneventacc[0] = 0;
  for(int j=1; j<gNTrgSetup; j++)
    {
      hPSdimuonInTrig[j-1] = (TH1F*)hPSdimuon->Clone(Form("%s_%s",hPSdimuon->GetName(),trgSetupName[j-1]));
      hPSdimuonInTrig[j-1]->Reset("C");
      hLTdimuonInTrig[j-1] = (TH1F*)hLTdimuon->Clone(Form("%s_%s",hLTdimuon->GetName(),trgSetupName[j-1]));
      hLTdimuonInTrig[j-1]->Reset("C");
      hNdimuonInTrig[j-1] = (TH1F*)hNeventsTake->Clone(Form("%s_%s",hNeventsTake->GetName(),trgSetupName[j-1]));
      hNdimuonInTrig[j-1]->Reset("C");
      hPSmbInTrig[j-1] = (TH1F*)hPSmb->Clone(Form("%s_%s",hPSmb->GetName(),trgSetupName[j-1]));
      hPSmbInTrig[j-1]->Reset("C");
      hLTmbInTrig[j-1] = (TH1F*)hLTmb->Clone(Form("%s_%s",hLTmb->GetName(),trgSetupName[j-1]));
      hLTmbInTrig[j-1]->Reset("C");
      hNmbInTrig[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("%s_%s",hNeventsTakeMb->GetName(),trgSetupName[j-1]));
      hNmbInTrig[j-1]->Reset("C");
      hEqvRatioInTrig[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEqvRatioInTrig_%s",trgSetupName[j-1]));
      hEqvRatioInTrig[j-1]->Reset("AC");

      dimuonevent[j] = 0;
      dimuoneventacc[j] = 0;
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type,trgSetupName[j-1]));
      int runnumber;
      while(fruns >> runnumber)
	{
	  int bin = hEvtRunAcc->FindBin(runnumber);
	  if(bin<1 || bin>hEvtRunAcc->GetNbinsX()) continue;
	  if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
	  int lumiBin = hNeventsTake->FindFixBin(runnumber);

	  dimuoneventacc[j] += hEvtRunAcc->GetBinContent(bin);
	  dimuoneventacc[0] += hEvtRunAcc->GetBinContent(bin);
	  dimuonevent[j] += hNeventsTake->GetBinContent(lumiBin);
	  dimuonevent[0] += hNeventsTake->GetBinContent(lumiBin);

	  double nEventsTaken = hNeventsTake->GetBinContent(lumiBin);
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
	  hPSdimuonInTrig[j-1]->SetBinContent(lumiBin, hPSdimuon->GetBinContent(lumiBin));
	  hPSdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hLTdimuonInTrig[j-1]->SetBinContent(lumiBin, hLTdimuon->GetBinContent(lumiBin));
	  hLTdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hNdimuonInTrig[j-1]->SetBinContent(lumiBin, hNeventsTake->GetBinContent(lumiBin));
	  hNdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hPSmbInTrig[j-1]->SetBinContent(lumiBin, hPSmb->GetBinContent(lumiBin));
	  hPSmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hLTmbInTrig[j-1]->SetBinContent(lumiBin, hLTmb->GetBinContent(lumiBin));
	  hLTmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hNmbInTrig[j-1]->SetBinContent(lumiBin, hNeventsTakeMb->GetBinContent(lumiBin));
	  hNmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  //hEqvRatioInTrig[j-1]->SetBinContent(lumiBin, nEventsRun/rf/nEventsTaken * hEqMbEvents[0]->GetBinContent(hEqMbEvents[0]->FindFixBin(runnumber))/hNeventsTake->GetBinContent(lumiBin)/hPSdimuon->GetBinContent(lumiBin));
	  double scalefactor = 1;
	  //if(j==1 || j==2) scalefactor = 1.5;
	  hEqvRatioInTrig[j-1]->SetBinContent(lumiBin, scalefactor*hNeventsTakeMb->GetBinContent(lumiBin)*hPSmb->GetBinContent(lumiBin)/hNeventsTake->GetBinContent(lumiBin)/hPSdimuon->GetBinContent(lumiBin));
	  hEqvRatioInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	}
      for(int i=0; i<nCentBins; i++)
	{
	  printf("[i] # of events for prod%s in %s: %1.0f, %4.2f%%\n",trgSetupName[j-1],cent_Name[i],mb_events[i][j],mb_events[i][j]/mb_events[i][0]*100);
	}
    }
  for(int j=1; j<gNTrgSetup; j++)
    {
      cout << "All: " << j << " -> " <<  dimuonevent[j] << " = " << dimuonevent[j]*1./dimuonevent[0]*100 << endl;
      cout << "Acc: " << j << " -> " <<  dimuoneventacc[j] << " = " << dimuoneventacc[j]*1./dimuoneventacc[0]*100 << endl;
    }

  for(int i=0; i<nCentBins; i++)
    {
      printf("[i]Cent %s, Check %1.0f =? %1.0f\n",cent_Name[i],mb_events[i][0],mb_events[i][1]+mb_events[i][2]+mb_events[i][3]+mb_events[i][4]);
    }
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================

  // check the pre-scale and live-time
  TCanvas *c = new TCanvas("check_ps_lt", "check_ps_lt", 1200, 700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  const int firstrun = 15103010;
  const int lastrun  = 15103063;

  for(int j=1; j<gNTrgSetup; j++)
    {
      c->cd(1);
      hPSdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hPSdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hPSdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hPSdimuonInTrig[j-1]->SetTitle("Dimuon: pre-scale;runId;");
      hPSdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5);
      hPSdimuonInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) hPSdimuonInTrig[j-1]->Draw("P");
      else     hPSdimuonInTrig[j-1]->Draw("samesP");
      leg->AddEntry(hPSdimuonInTrig[j-1], Form("prod%s",trgSetupName[j-1]), "P");

      c->cd(2);
      hNdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hNdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hNdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hNdimuonInTrig[j-1]->SetTitle("# of dimuon events on tape;runId;");
      hNdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e6);
      hNdimuonInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) hNdimuonInTrig[j-1]->Draw("P");
      else     hNdimuonInTrig[j-1]->Draw("samesP");

      /*
      c->cd(3);
      hLTdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hLTdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hLTdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hLTdimuonInTrig[j-1]->SetTitle("Dimuon live-time");
      hLTdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTdimuonInTrig[j-1]->Draw("P");
      else     hLTdimuonInTrig[j-1]->Draw("samesP");
      */
      c->cd(3);
      hEqvRatioInTrig[j-1]->SetMarkerStyle(19+j);
      hEqvRatioInTrig[j-1]->SetMarkerColor(color[j-1]);
      hEqvRatioInTrig[j-1]->SetLineColor(color[j-1]);
      hEqvRatioInTrig[j-1]->SetTitle(";runId;");
      hEqvRatioInTrig[j-1]->GetYaxis()->SetRangeUser(0,50);
      hEqvRatioInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) 
	{
	  hEqvRatioInTrig[j-1]->Draw("P");
	  //TPaveText *t1 = GetTitleText("Equivalent # of MB events after vertex cuts per dimuon event (0-80%)",0.035);
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events on tape per dimuon event",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrig[j-1]->Draw("samesP");

      c->cd(4);
      hPSmbInTrig[j-1]->SetMarkerStyle(19+j);
      hPSmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hPSmbInTrig[j-1]->SetLineColor(color[j-1]);
      hPSmbInTrig[j-1]->SetTitle("MB: pre-scale;runId;");
      hPSmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e4);
      hPSmbInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) hPSmbInTrig[j-1]->Draw("P");
      else     hPSmbInTrig[j-1]->Draw("samesP");

      c->cd(5);
      hNmbInTrig[j-1]->SetMarkerStyle(19+j);
      hNmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hNmbInTrig[j-1]->SetLineColor(color[j-1]);
      hNmbInTrig[j-1]->SetTitle("# of MB events on tape;runId;");
      hNmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5e4);
      hNmbInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) hNmbInTrig[j-1]->Draw("P");
      else     hNmbInTrig[j-1]->Draw("samesP");

      /*
      c->cd(6);
      hLTmbInTrig[j-1]->SetMarkerStyle(19+j);
      hLTmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hLTmbInTrig[j-1]->SetLineColor(color[j-1]);
      hLTmbInTrig[j-1]->SetTitle("MB live-time");
      hLTmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTmbInTrig[j-1]->Draw("P");
      else     hLTmbInTrig[j-1]->Draw("samesP");
      */
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sTrgSetupComp_PSandLT.pdf",run_type,run_config));

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
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/TrgSetupComp_MtdAcc.pdf",run_type));
  // =============================================
  //
  //
  //return;

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
      printf("[i] Jpsi trigger efficiency for %s: %4.3f\n",gTrgSetupTitle[j], hJpsiTrigEffVsCent[j]->GetBinContent(1));

      hJpsiTrigEffLS[j]       = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      hJpsiTrigEffLS[j]->SetName(Form("%s_LS",hJpsiTrigEffLS[j]->GetName()));
      if(j>0) hJpsiTrigEffLS[j]->Divide(hJpsiTrigEffLS[0]);

      hJpsiTrigEffVsCentLS[j] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsCent_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
      hJpsiTrigEffVsCentLS[j]->SetName(Form("%s_LS",hJpsiTrigEffVsCentLS[j]->GetName()));
      printf("[i] LS trigger efficiency for %s: %4.3f\n",gTrgSetupTitle[j],hJpsiTrigEffVsCentLS[j]->GetBinContent(1));      
    }

  const double muonPidScale[gNTrgSetup] = {1, 1.008, 1.006, 1.004, 0.99};
  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.EmbJpsiEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiEff[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiEff[k][j] = (TH1F*)fEff->Get(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[j]));
	  //if(j>0) hJpsiEff[k][j]->Multiply(hJpsiTrigEff[j]);
	  if(j>0) hJpsiEff[k][j]->Multiply(hJpsiTrigEffLS[j]);
	  hJpsiEff[k][j]->Scale(muonPidScale[j]);
	  hJpsiEff[k][j]->Multiply(hAccCorr[j]);
	}
    }

  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  TH1F *hModel = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hWeight = (TH1F*)hJpsiEff[0][0]->Clone("hWeight");
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
  for(int j=0; j<gNTrgSetup; j++)
    {
      total_eff[j] = 0;
      for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
	{
	  total_eff[j] += hWeight->GetBinContent(bin) * hJpsiEff[0][j]->GetBinContent(bin);
	  //total_eff[j] += hWeight->GetBinContent(bin) * hAccCorr[j]->GetBinContent(bin);
	}
      printf("[i] Total efficiency for %s is %4.3f%%\n",gTrgSetupTitle[j],total_eff[j]*100);
    }

  TCanvas *c = new TCanvas("CompEff", "CompEff", 1100, 700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  TH1F *hRatio = (TH1F*)hJpsiEff[k][j]->Clone(Form("%s_clone",hJpsiEff[k][j]->GetName()));
	  hRatio->Divide(hJpsiEff[k][0]);
	  hRatio->SetMarkerStyle(20+j);
	  hRatio->SetMarkerColor(color[j]);
	  hRatio->SetLineColor(color[j]);
	  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
	  if(j==1) hRatio->Draw();
	  else     hRatio->Draw("sames");
	  if(k==0)leg->AddEntry(hRatio, gTrgSetupTitle[j], "P");
	}
    }
  c->cd(1);
  leg->Draw();
  return;

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

  // systematic uncertianty on J/psi signal extraction
  const int checkSys = 1;
  if(checkSys)
    {
      TFile *fsys = TFile::Open(Form("Rootfiles/%s.StudyLumiDep.root",run_type), "read");
      TH1F *hSys[nCentBins][gNTrgSetup];
      for(int k=0; k<nCentBins; k++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      hSys[k][j]    = (TH1F*)fsys->Get(Form("Sys_signalExt_cent%s%s",cent_Title[k],gTrgSetupName[j]));
	      int nbins = hJpsiCounts[k][j]->GetNbinsX();
	      for(int bin=1; bin<=nbins; bin++)
		{
		  double err = sqrt( pow(hSys[k][j]->GetBinContent(bin), 2) + pow(hJpsiCounts[k][j]->GetBinError(bin)/hJpsiCounts[k][j]->GetBinContent(bin), 2) );
		  hJpsiCounts[k][j]->SetBinError(bin, hJpsiCounts[k][j]->GetBinContent(bin)*err);
		}
	    }
	}
    }

  // calcualte statistical power
  double jpsiCounts[gNTrgSetup];
  double jpsiCountsErr[gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {
      jpsiCounts[j] = 0;
      jpsiCountsErr[j] = 0;
      for(int bin=1; bin<=hJpsiCounts[0][j]->GetNbinsX(); bin++)
	{
	  jpsiCounts[j] += hJpsiCounts[0][j]->GetBinContent(bin);
	  jpsiCountsErr[j] += pow(hJpsiCounts[0][j]->GetBinError(bin), 2);
	}
      jpsiCountsErr[j] = sqrt(jpsiCountsErr[j]);
      printf("[i] Cent %s, %s: jpsi counts = %1.0f #pm %1.0f (%4.2f%%)\n",cent_Name[0],gTrgSetupName[j],jpsiCounts[j],jpsiCountsErr[j],jpsiCountsErr[j]/jpsiCounts[j]*100);
    }
  
  return;

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins][gNTrgSetup];
  TH1F *hJpsiInvRatio[nCentBins][gNTrgSetup];
  for(int k=0; k<nCentBins; k++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[k][j] = (TH1F*)hJpsiCounts[k][j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[k],gTrgSetupTitle[j]));
	  hJpsiInvYield[k][j]->Divide(hJpsiEff[k][j]);
	  //cout << hJpsiEff[k][j]->GetBinContent(1) << endl;
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
  TH1F *hspread[nCentBins][gNTrgSetup];
  TRandom3 *gRand = new TRandom3();
  gRand->SetSeed(0);
  TCanvas *cspread = new TCanvas("cSpread","cSpread",1100,700);
  cspread->Divide(2, 2);
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
	  c->cd(k+1);
	  hJpsiInvRatio[k][j]->Draw("sames");
	  TF1 *func = new TF1(Form("Fit_%s",hJpsiInvRatio[k][j]->GetName()), "[0]", 0, 15);
	  if(k==2 || k==3) func->SetRange(0,10);
	  if(k==3 && j==3) func->SetRange(0,6);
	  if(k==4) func->SetRange(0,6);
	  hJpsiInvRatio[k][j]->Fit(func,"0QR");
	  func->SetLineColor(hJpsiInvRatio[k][j]->GetMarkerColor());
	  func->SetLineStyle(2);
	  func->Draw("sames");
	  

	  if(k==0)
	    {
	      hspread[k][j] = new TH1F(Form("hspread_%d_%d",k,j),"",150,0,1.5);
	      for(int iexpr=0; iexpr<4000; iexpr++)
		{
		  TH1F *htmp = (TH1F*)hJpsiInvRatio[k][j]->Clone(Form("%s_tmp%d",hJpsiInvRatio[k][j]->GetName(),iexpr));
		  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		    {
		      htmp->SetBinContent(bin, gRand->Gaus(hJpsiInvRatio[k][j]->GetBinContent(bin), hJpsiInvRatio[k][j]->GetBinError(bin)));
		      htmp->SetBinError(bin, hJpsiInvRatio[k][j]->GetBinError(bin)/hJpsiInvRatio[k][j]->GetBinContent(bin)*htmp->GetBinContent(bin));
		    }
		  TF1 *functmp = new TF1(Form("Fit_%s",htmp->GetName()), "[0]", 0, 10);
		  if(k==2 || k==3) functmp->SetRange(0,10);
		  if(k==3 && j==3) functmp->SetRange(0,6);
		  if(k==4) functmp->SetRange(0,6);
		  htmp->Fit(functmp,"0QR");
		  hspread[k][j]->Fill(functmp->GetParameter(0));
		}
	      cspread->cd(k*gNTrgSetup+j);
	      hspread[k][j]->SetMarkerStyle(20);
	      hspread[k][j]->Draw();
	      TF1 *funcSpread = new TF1(Form("Fit_%s",hspread[k][j]->GetName()), "gaus", 0.8, 1.5);
	      if(j==2) funcSpread->SetRange(1.0, 1.5);
	      funcSpread->SetParameter(1, func->GetParameter(0));
	      funcSpread->SetParameter(2, func->GetParError(0));
	      hspread[k][j]->Fit(funcSpread, "0QR");
	      funcSpread->SetLineColor(2);
	      funcSpread->Draw("sames");
	      printf("[i] cent %s, %10s: gaus = %4.3f#pm%4.3f\n",cent_Title[k],gTrgSetupTitle[j], funcSpread->GetParameter(1), funcSpread->GetParameter(2));
	      TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.8,0.055);
	      t1->AddText(Form("mean = %4.3f",funcSpread->GetParameter(1)));
	      t1->AddText(Form("width = %4.3f",funcSpread->GetParameter(2)));
	      t1->Draw();
	      TPaveText *t1 = GetTitleText(Form("%s%%: randomoization exercise (2014%s)",cent_Name[k],gTrgSetupTitle[j]),0.055);
	      t1->Draw();
	    }
	  
	  
	  double chi2 = 0;
	  int ndf = 0;
	  for(int bin=1; bin<=hJpsiInvRatio[k][j]->GetNbinsX(); bin++)
	    {
	      if(hJpsiInvRatio[k][j]->GetBinContent(bin)<=0) continue;
	      ndf++;
	      chi2 += pow( (hJpsiInvRatio[k][j]->GetBinContent(bin)-1)/hJpsiInvRatio[k][j]->GetBinError(bin), 2);
	    }
	  printf("[i] cent %s, %10s: p0 = %4.3f #pm %4.3f, %4.3fsigma\n",cent_Title[k],gTrgSetupTitle[j],func->GetParameter(0),func->GetParError(0),fabs(func->GetParameter(0)-1)/func->GetParError(0));
	  //printf("chi2 = %4.2f/%d, p-value = %4.2f\n",chi2,ndf,TMath::Prob(chi2, ndf));
	  if(k==0) leg->AddEntry(hJpsiInvRatio[k][j],Form("Run14_AuAu200%s",gTrgSetupTitle[j]),"P");
	}
      c->cd(k+1);
      TPaveText *t1 = GetPaveText(0.25,0.35,0.8,0.9,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  if(savePlot) cSpread->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sTrgSetupComp_RandFitError.pdf",run_type,run_config));
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
void xsec_Run14(const bool savePlot = 1, const bool saveHisto = 0)
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
  TGraphErrors *gAuAuRun11[3];
  TGraphErrors *gAuAuRun11Sys[3];
  TH1F *hTBW[3];
  for(int i=0; i<3; i++)
    {
      gAuAuLowPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[i+1]));
      gAuAuLowPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[i+1]));
      gAuAuHighPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[i+1]));
      gAuAuHighPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[i+1]));
      gAuAuRun11[i] = (TGraphErrors*)fpub->Get(Form("Run11_JpsInvYield_cent%s",cent_Title[i+1]));
      gAuAuRun11Sys[i] = (TGraphErrors*)fpub->Get(Form("Run11_JpsInvYieldSys_cent%s",cent_Title[i+1]));
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
	  gAuAuRun11[k-1]->SetMarkerStyle(20);
	  gAuAuRun11[k-1]->SetMarkerColor(9);
	  gAuAuRun11[k-1]->SetMarkerSize(1.5);
	  gAuAuRun11[k-1]->SetLineColor(9);
	  gAuAuRun11[k-1]->Draw("sames PE");
	  gAuAuRun11Sys[k-1]->SetMarkerColor(9);
	  gAuAuRun11Sys[k-1]->SetLineColor(9);
	  gAuAuRun11Sys[k-1]->SetFillStyle(0);
	  gAuAuRun11Sys[k-1]->Draw("sameE5");
	  hTBW[k-1]->Draw("sames");

	  gAuAuYield_sQM[k-1]->SetMarkerStyle(29);
	  gAuAuYield_sQM[k-1]->SetMarkerColor(kGreen+2);
	  gAuAuYield_sQM[k-1]->SetLineColor(kGreen+2);
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
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_compareToPub.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYieldVsPt_compareToPub.png",run_type,run_config));
    }


  TFile *fSys = TFile::Open(Form("Rootfiles/%s.Sys.JpsiXsec.root",run_type),"read");
  TH1F *hAuAuJpsiSys[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hAuAuJpsiSys[k]= (TH1F*)fSys->Get(Form("JpsiSysVsPt_All_cent%s",cent_Title[k]));
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi_ratio","AuAu200_Jpsi_ratio",1100,700);
  c->Divide(2,2);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  double x,y;
  for(int k=0; k<3; k++)
    {
      TH1F *hRatio = (TH1F*)hJpsiInvYield[k+1]->Clone(Form("%s_ratio",hJpsiInvYield[k+1]->GetName()));
      int endbin = hRatio->FindFixBin(9);
      TGraphErrors *gRatioSys = new TGraphErrors(hJpsiInvYield[k+1]->GetNbinsX());

      TGraphErrors *gRatioRun11   = (TGraphErrors*)gAuAuRun11[k]->Clone(Form("%s_clone",gAuAuRun11[k]->GetName()));
      TGraphErrors *gRatioRun11Sys = (TGraphErrors*)gAuAuRun11Sys[k]->Clone(Form("%s_clone",gAuAuRun11Sys[k]->GetName()));
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
	  gRatioSys->SetPoint(bin-1, hRatio->GetBinCenter(bin), hRatio->GetBinContent(bin));
	  gRatioSys->SetPointError(bin-1, 0.2, hAuAuJpsiSys[k+1]->GetBinContent(bin)*hRatio->GetBinContent(bin));

	  // Run11 results
	  if(bin<=5)
	    {
	      gRatioRun11->GetPoint(bin-1, x, y);
	      gRatioRun11->SetPoint(bin-1, x, y/scale);
	      gRatioRun11->SetPointError(bin-1, 0, gRatioRun11->GetErrorY(bin-1)/scale);
	      gRatioRun11Sys->SetPoint(bin-1, x, y/scale);
	      gRatioRun11Sys->SetPointError(bin-1, 0.2, gRatioRun11Sys->GetErrorY(bin-1)/scale);
	    }
	}
      c->cd(k+1);
      gPad->SetGridy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
 
      // final
      hRatio->SetTitle(";p_{T} (GeV/c);Ratio to TBW");
      hRatio->GetXaxis()->SetRangeUser(0.15,10);
      hRatio->GetYaxis()->SetRangeUser(0,2.5);
      ScaleHistoTitle(hRatio,0.06,1,0.05,0.06,1,0.05,62);
      hRatio->Draw();
      gRatioSys->SetFillStyle(0);
      gRatioSys->SetLineColor(2);
      gRatioSys->Draw("samesE5");

      // published
      gRatioRun11Sys->Draw("samesE5");
      gRatioRun11->Draw("samesPE");

      TGraphAsymmErrors *gRatioLowPt    = (TGraphAsymmErrors*)gAuAuLowPt[k]->Clone(Form("%s_ratio",gAuAuLowPt[k]->GetName()));
      TGraphAsymmErrors *gRatioLowPtSys = (TGraphAsymmErrors*)gAuAuLowPtSys[k]->Clone(Form("%s_ratio",gAuAuLowPtSys[k]->GetName()));
      for(int ipoint=0; ipoint<gRatioLowPt->GetN(); ipoint++)
	{
	  gAuAuLowPt[k]->GetPoint(ipoint, x, y);
	  double scale = hTBW[k]->GetBinContent(hTBW[k]->FindFixBin(x));
	  gRatioLowPt->SetPoint(ipoint, x, y/scale);
	  gRatioLowPt->SetPointError(ipoint, gAuAuLowPt[k]->GetErrorXlow(ipoint), gAuAuLowPt[k]->GetErrorXhigh(ipoint),
				     gAuAuLowPt[k]->GetErrorYlow(ipoint)/scale, gAuAuLowPt[k]->GetErrorYhigh(ipoint)/scale);
	  gRatioLowPtSys->SetPoint(ipoint, x, y/scale);
	  gRatioLowPtSys->SetPointError(ipoint, gAuAuLowPtSys[k]->GetErrorXlow(ipoint), gAuAuLowPtSys[k]->GetErrorXhigh(ipoint),
					gAuAuLowPtSys[k]->GetErrorYlow(ipoint)/scale, gAuAuLowPtSys[k]->GetErrorYhigh(ipoint)/scale);
	}
      gRatioLowPtSys->Draw("samesE5");  
      gRatioLowPt->Draw("samesPEZ");

      TGraphAsymmErrors *gRatioHighPt    = (TGraphAsymmErrors*)gAuAuHighPt[k]->Clone(Form("%s_ratio",gAuAuHighPt[k]->GetName()));
      TGraphAsymmErrors *gRatioHighPtSys = (TGraphAsymmErrors*)gAuAuHighPtSys[k]->Clone(Form("%s_ratio",gAuAuHighPtSys[k]->GetName()));
      for(int ipoint=0; ipoint<gRatioHighPt->GetN(); ipoint++)
	{
	  gAuAuHighPt[k]->GetPoint(ipoint, x, y);
	  double scale = hTBW[k]->GetBinContent(hTBW[k]->FindFixBin(x));
	  gRatioHighPt->SetPoint(ipoint, x, y/scale);
	  gRatioHighPt->SetPointError(ipoint, gAuAuHighPt[k]->GetErrorXlow(ipoint), gAuAuHighPt[k]->GetErrorXhigh(ipoint),
				     gAuAuHighPt[k]->GetErrorYlow(ipoint)/scale, gAuAuHighPt[k]->GetErrorYhigh(ipoint)/scale);
	  gRatioHighPtSys->SetPoint(ipoint, x, y/scale);
	  gRatioHighPtSys->SetPointError(ipoint, gAuAuHighPtSys[k]->GetErrorXlow(ipoint), gAuAuHighPtSys[k]->GetErrorXhigh(ipoint),
				     gAuAuHighPtSys[k]->GetErrorYlow(ipoint)/scale, gAuAuHighPtSys[k]->GetErrorYhigh(ipoint)/scale);
	}
      gRatioHighPtSys->Draw("samesE5");  
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
  leg->AddEntry(gAuAuLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (Low p_{T}, 2010)","P");
  leg->AddEntry(gAuAuHighPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (High p_{T})","P");
  leg->AddEntry(gAuAuRun11[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1 (Low p_{T}, 2011)","P");
  leg->AddEntry(hJpsiInvYield[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");
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

