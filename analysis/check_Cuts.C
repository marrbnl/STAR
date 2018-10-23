#include <map>
const int nPtBins        = 9;
const double xPtBins[10] = {0,1,2,3,4,5,6,8,10,15};
const int nCentBins      = nCentBins_pt; 
const int* centBinsLow   = centBins_low_pt;
const int* centBinsHigh  = centBins_high_pt;
const char** cent_Name    = cent_Name_pt;
const char** cent_Title   = cent_Title_pt;
const char* gPairName[2] = {"UL", "LS"};
const int gPairType = 1;
const TString legNameLumi[gNTrgSetup] = {"inclusive", "prod", "prod_low", "prod_mid", "prod_high"};

//================================================
void check_Cuts()
{
  gStyle->SetOptStat(0);

  //makeHistos();
  xsec();
  //LSpairInMB();
  //MbVsDm();
  //JpsiTree();
  //dataQa();
}

//================================================
void xsec(const int savePlot = 0, const int saveHisto = 0)
{
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

  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  printf("[i] Process %s\n",fdata->GetName());
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");

  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");

  const char *trgSetupName[4] = {"","_low","_mid","_high"};
  const int eqMbMethod = 1;
  if(eqMbMethod==0)
    {
      TH1F *hEqMbEvents[nCentBins];
      for(int i=0; i<nCentBins; i++)
	{
	  hEqMbEvents[i] = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
	}
      for(int j=1; j<gNTrgSetup; j++)
	{
	  ifstream fruns;
	  fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type.Data(),trgSetupName[j-1]));
	  int runnumber;
	  int nRuns = 0;
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
	      double nEventsAcc = hEvtRunAcc->GetBinContent(bin);
	      double rf = hRF->GetBinContent(hRF->FindFixBin(runnumber));
	      if(rf==0)
		{
		  printf("[w] rf = 0 for run %1.0f\n",runnumber);
		  rf = 0.49;
		}

	      if(runnumber==15149073 || runnumber==15166014 || runnumber==15166015) continue;

	      for(int i=0; i<nCentBins; i++)
		{
		  double eq_mb = hEqMbEvents[i]->GetBinContent(hEqMbEvents[i]->FindFixBin(runnumber));
		  mb_events[i][j] += nEventsRun/rf/nEventsTaken * eq_mb;
		  mb_events[i][0] += nEventsRun/rf/nEventsTaken * eq_mb;
		}
	      nRuns ++;
	    }
	  printf("[i] Total # of runs for %s = %d\n",trgSetupName[j-1],nRuns);
	}
    }
  if(eqMbMethod==1)
    {
      TH1F *hEqMb[gNTrgSetup];
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hEqMb[j] = (TH1F*)fLumi->Get(Form("hEqMbEvt_ScaleMbData%s",gTrgSetupTitle[j]));
	  for(int i=0; i<nCentBins; i++)
	    {
	      mb_events[i][j] = hEqMb[j]->Integral((centBinsLow[i]+1)/2, centBinsLow[i]/2);
	    }
	}
    }
  for(int j=0; j<gNTrgSetup; j++)
    {
      printf("[i] MB events for %s%s is: %1.0f (%4.2f%%)\n",run_type.Data(),gTrgSetupTitle[j],mb_events[0][j],mb_events[0][j]/mb_events[0][0]*100);
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
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type.Data(),trgSetupName[j-1]));
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
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type.Data()),"read");
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
  TFile *fTrigEffLS = TFile::Open(Form("Rootfiles/Run14_AuAu200.StudyTrigEff.root"),"read");
  TH1F *hJpsiTrigEffLS[gNTrgSetup];
  TH1F *hJpsiTrigEffVsCentLS[gNTrgSetup];
  for(int j=0; j<gNTrgSetup; j++)
    {

      if(j==0)
	{
	  hJpsiTrigEffLS[j]       = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
	  hJpsiTrigEffVsCentLS[j] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsCent_TacDiffEff%s_TrigStudy",gTrgSetupTitle[j]));
	}
      else
	{
	  hJpsiTrigEffLS[j]       = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_M1_TrigStudy",gTrgSetupTitle[j]));
	  hJpsiTrigEffVsCentLS[j] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsCent_TacDiffEff%s_M1_TrigStudy",gTrgSetupTitle[j]));
	}
      hJpsiTrigEffLS[j]->SetName(Form("%s_LS",hJpsiTrigEffLS[j]->GetName()));
      if(j>0) 
	{
	  hJpsiTrigEffLS[j]->Divide(hJpsiTrigEffLS[0]);
	  printf("[i] LS trigger efficiency for %s: %4.3f\n",gTrgSetupTitle[j],hJpsiTrigEffVsCentLS[j]->GetBinContent(1));      
	}
    }

  TFile *fEff = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"read");
  printf("[w] efficiency name: %s\n",fEff->GetName());
  TH1F *hJpsiEff[nCentBins][gNTrgSetup];
  TH1F *hJpsiEffTpc[nCentBins][gNTrgSetup];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiEff[i][j] = (TH1F*)fEff->Get(Form("LSEffVsPt_cent%s%s_final",cent_Title[i],gTrgSetupTitle[j]));
	  hJpsiEffTpc[i][j] = (TH1F*)hJpsiEff[i][j]->Clone(Form("LSTpcEffVsPt_cent%s%s_final",cent_Title[i],gTrgSetupTitle[j]));
	  if(j>0) hJpsiEff[i][j]->Multiply(hJpsiTrigEffLS[j]);
	  hJpsiEff[i][j]->Multiply(hAccCorr[j]);
	}
    }

  // additional TPC tracking inefficiency
  TFile *fLumiDepEff = TFile::Open(Form("Rootfiles/%s.Sys.LumiDepEff.root",run_type.Data()), "read");
  TH2F *hJpsiCorrInPtBin = (TH2F*)fLumiDepEff->Get("hJpsiCorrInPtBin_Sys0");
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  double eff = hJpsiCorrInPtBin->GetBinContent(i+1, gNTrgSetup-j);
	  hJpsiEff[i][j]->Scale(eff);
	}
    }

  // TPC efficiency
  TCanvas *c = new TCanvas("CompEff", "CompEff", 800, 600);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  for(int j=1; j<gNTrgSetup; j++)
    {
      TH1F *hRatio = (TH1F*)hJpsiEffTpc[0][j]->Clone(Form("%s_clone",hJpsiEff[0][j]->GetName()));
      hRatio->Divide(hJpsiEffTpc[0][2]);
      hRatio->SetMarkerStyle(20+j);
      hRatio->SetMarkerColor(gColor[j]);
      hRatio->SetLineColor(gColor[j]);
      hRatio->GetYaxis()->SetRangeUser(0.1,1.1);
      if(j==1) hRatio->Draw();
      else     hRatio->Draw("sames");
      leg->AddEntry(hRatio, gTrgSetupTitle[j], "P");
    }
  leg->Draw();

  // Jpsi raw counts
  TH1F *hJpsiCounts[nCentBins][gNTrgSetup];
  TFile *fYield = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiCounts[i][j] = (TH1F*)fYield->Get(Form("%s_PairPt_%s_cent%s%s",run_config,gPairName[gPairType],cent_Title[i],gTrgSetupTitle[j]));
	  printf("[i] %s for %s in %s%%: %1.0f%\n",gPairName[gPairType],gTrgSetupTitle[j],cent_Name[i],hJpsiCounts[i][j]->Integral());
	}
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins][gNTrgSetup];
  TH1F *hJpsiInvRatio[nCentBins][gNTrgSetup];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[i][j] = (TH1F*)hJpsiCounts[i][j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[i],gTrgSetupTitle[j]));
	  hJpsiInvYield[i][j]->Divide(hJpsiEff[i][j]);
	  for(int bin=1; bin<=hJpsiInvYield[i][j]->GetNbinsX(); bin++)
	    {
	      double bin_width = hJpsiInvYield[i][j]->GetBinWidth(bin); // dpT
	      double bin_center = hJpsiInvYield[i][j]->GetBinCenter(bin); // pT 
	      hJpsiInvYield[i][j]->SetBinContent(bin,hJpsiInvYield[i][j]->GetBinContent(bin)/bin_width/bin_center);
	      hJpsiInvYield[i][j]->SetBinError(bin,hJpsiInvYield[i][j]->GetBinError(bin)/bin_width/bin_center);
	    }
	  hJpsiInvYield[i][j]->Scale(1./(2*pi)); // 2pi
	  if(j==1) hJpsiInvYield[i][0] = (TH1F*)hJpsiInvYield[i][j]->Clone(Form("Jpsi_InvYield_cent%s%s",cent_Title[i],gTrgSetupTitle[0]));
	  else     hJpsiInvYield[i][0]->Add((TH1F*)hJpsiInvYield[i][j]);
	}
    }

  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiInvYield[i][j]->Scale(1./mb_events[i][j]); // N_evt
	  hJpsiInvYield[i][j]->SetMarkerStyle(21);
	  hJpsiInvYield[i][j]->SetMarkerColor(gColor[j]);
	  hJpsiInvYield[i][j]->SetLineColor(gColor[j]);
	  hJpsiInvYield[i][j]->SetMarkerSize(1.2);

	  hJpsiInvRatio[i][j] = (TH1F*)hJpsiInvYield[i][j]->Clone(Form("Jpsi_InvYieldRatio_cent%s%s",cent_Title[i],gTrgSetupTitle[j]));
	  hJpsiInvRatio[i][j]->Divide(hJpsiInvYield[i][0]);
	  for(int bin=1; bin<=hJpsiInvRatio[i][j]->GetNbinsX(); bin++)
	    {
	      double err = hJpsiInvYield[i][j]->GetBinError(bin)/hJpsiInvYield[i][j]->GetBinContent(bin);
	      hJpsiInvRatio[i][j]->SetBinError(bin,hJpsiInvRatio[i][j]->GetBinContent(bin)*err);
	    }
	}
    }

  TCanvas *c = new TCanvas("Ratio_JetXsec","Ratio_JetXsec",1100,700);
  c->Divide(3,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);Ratio to inclusive",15,0,15);
  ScaleHistoTitle(hAuAu,0.045,1,0.04,0.045,1.4,0.04,62);
  TLegend *leg2 = new TLegend(0.3,0.4,0.5,0.8);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.055);
  leg2->SetHeader(run_type.Data());
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      if(i>=3) hAuAu->GetYaxis()->SetRangeUser(0.3,1.6);
      else  hAuAu->GetYaxis()->SetRangeUser(0.6,1.3);
      hAuAu->GetXaxis()->SetRangeUser(0,8);
      hAuAu->DrawCopy();
      TLegend *leg = new TLegend(0.18,0.15,0.9,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.035); 
      leg->SetNColumns(2);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hJpsiInvRatio[i][j]->Draw("sames");
	  TF1 *func = new TF1(Form("Fit_%s",hJpsiInvRatio[i][j]->GetName()), "[0]", 0, 8);
	  hJpsiInvRatio[i][j]->Fit(func,"0QR");
	  func->SetLineColor(hJpsiInvRatio[i][j]->GetMarkerColor());
	  func->SetLineStyle(2);
	  func->Draw("sames");
	  leg->AddEntry(func,Form("%4.3f +/- %4.3f",func->GetParameter(0),func->GetParError(0)),"L");
	  if(i==0) leg2->AddEntry(hJpsiInvRatio[i][j], gLegNameTrg[j].Data(), "PL");
	}
      leg->Draw();
      TPaveText *t1 = GetTitleText(Form("%s%%: %s pairs",cent_Name[i],gPairName[gPairType]),0.055);
      t1->Draw();
    }
  c->cd(6);
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%sYieldRatio.pdf",run_type.Data(),gPairName[gPairType]));

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiMb.root","update");
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hAccCorr[j]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void MbVsDm(const int savePlot = 0)
{
  const int nFile = 2;
  const char* trig_name[3] = {"di_mu", "mb", "mb_ps"};
  TFile *fstudy[nFile];
  fstudy[0] = TFile::Open("output/Run14_AuAu200.Study.root", "read");
  fstudy[1] = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");

  // check # of events
  double nMbEvtEq[nFile+1][gNTrgSetup];

  // mb
  TH2F* hMtdQaCentEqMb = (TH2F*)fstudy[1]->Get("mhMtdQaCentEqMb_mb");
  double prescale[4] = {1.0766, 1.118, 1.0438, 0.9827};
  nMbEvtEq[1][0] = 0;
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMtdQaCentEqMb->ProjectionY("",k,k);
      nMbEvtEq[1][k] = htmp->Integral(1,16);
      nMbEvtEq[1][0] += nMbEvtEq[1][k];
    }
  // dimuon
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEqMb[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      hEqMb[k] = (TH1F*)fLumi->Get(Form("hEqMbEvt_ScaleMbData%s",gTrgSetupTitle[k]));
      nMbEvtEq[0][k] = hEqMb[k]->Integral(1,8);
    }

  TH2F *hMbEvt[2];
  TH2F *hMbEvtDm[2];
  for(int i=0; i<nFile; i++)
    {
      hMbEvt[i] = (TH2F*)fstudy[i]->Get(Form("mhMtdQaCent_%s",trig_name[i]));
      hMbEvtDm[i] = (TH2F*)fstudy[i]->Get(Form("mhMtdQaCentDm_%s",trig_name[i]));

      for(int k=1; k<gNTrgSetup; k++)
	{
	  TH1F *htmp  = (TH1F*)hMbEvt[i]->ProjectionY("",k,k);
	  TH1F *htmp2 = (TH1F*)hMbEvtDm[i]->ProjectionY("",k,k);
	  printf("[i] %s%s has %1.0f dm events and %1.0f eq. mb events - %4.3f%%\n",trig_name[i],gTrgSetupTitle[k],htmp2->GetEntries(),nMbEvtEq[i][k],htmp2->GetEntries()/nMbEvtEq[i][k]*100);
	}
    }

  // MB + prescale
  nMbEvtEq[2][0] = 0;
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMtdQaCentEqMb->ProjectionY("",k,k);
      nMbEvtEq[2][k] = htmp->Integral(1,16);
      nMbEvtEq[2][0] += nMbEvtEq[2][k];
    }


  // # of matched track pair
  // UL+LS M > 0 

  // mass dependence
  THnSparseF *hnInvMass[nFile];
  TH1F *hPairMass[nFile+1][gNTrgSetup];
  for(int i=0; i<nFile; i++)
    {
      hnInvMass[i] = (THnSparseF*)fstudy[i]->Get(Form("mhMtdQaInvMass_%s",trig_name[i])); 
      hnInvMass[i]->GetAxis(7)->SetRange(2,2); // fire trigger
      //hnInvMass[i]->GetAxis(5)->SetRange(2,2); // muon candidate
      hnInvMass[i]->GetAxis(6)->SetRange(2,2); // tagged as dimuon event
      //hnInvMass[i]->GetAxis(8)->SetRange(2,2); // tagged as mb event
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hnInvMass[i]->GetAxis(3)->SetRange(k, k);
	  hPairMass[i][k] = (TH1F*)hnInvMass[i]->Projection(1);
	  hPairMass[i][k]->SetName(Form("hPirMass_%s%s",trig_name[i],gTrgSetupTitle[k]));
	  hPairMass[i][k]->Sumw2();
	  hPairMass[i][k]->Rebin(4);
	  printf("[i] %s has %1.0f LS pairs\n",hPairMass[i][k]->GetName(),hPairMass[i][k]->Integral());
	  hPairMass[i][k]->Scale(1./nMbEvtEq[i][k]);
	  hnInvMass[i]->GetAxis(3)->SetRange(0,-1);
	}
    }

  TCanvas *c = new TCanvas("cPairMassComp", "cPairMassComp", 1100, 700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.3,0.15,0.5,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int k=2; k<gNTrgSetup; k++)
    {
      c->cd(k-1);
      gPad->SetLogy();
      for(int i=0; i<nFile; i++)
	{
	  hPairMass[i][k]->SetMarkerStyle(18+k+4*i);
	  hPairMass[i][k]->SetMarkerColor(gColor[k]);
	  hPairMass[i][k]->SetLineColor(gColor[k]);
	  hPairMass[i][k]->SetTitle("");
	  if(i==0) hPairMass[i][k]->Draw();
	  else     hPairMass[i][k]->Draw("sames");
	  if(k==2)
	    {
	      if(i==0) leg->AddEntry(hPairMass[i][k], "Dimuon", "P");
	      if(i==1) leg->AddEntry(hPairMass[i][k], "MB+Dimuon", "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s",gLegNameTrg[k].Data()), 0.055);
      t1->Draw();

      c->cd(k+2);
      TH1F *hRatio = (TH1F*)hPairMass[0][k]->Clone(Form("%s_ratio",hPairMass[0][k]->GetName()));
      hRatio->Divide(hPairMass[1][k]);
      hRatio->GetYaxis()->SetRangeUser(0.5,1.2);
      hRatio->Draw();
      TF1 *funcRatio = new TF1(Form("%s_func",hRatio->GetName()), "pol0", 0, 4);
      hRatio->Fit(funcRatio, "IR0Q");
      funcRatio->SetLineColor(hRatio->GetMarkerColor());
      funcRatio->SetLineStyle(2);
      funcRatio->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: DM/(MB+DM)",gLegNameTrg[k].Data()), 0.055);
      t1->Draw();
      leg2 = new TLegend(0.2,0.15,0.5,0.35);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.05);
      leg2->AddEntry(funcRatio, Form("%4.3f +/- %4.3f",funcRatio->GetParameter(0),funcRatio->GetParError(0)), "L");
      leg2->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_CompMthPairM_MbVsDm.pdf",run_type.Data()));
  return;

  // pt dependence
  const int nbins = 7;
  double xbins[8] = {0, 1, 2, 3, 4, 5, 6, 8};
  TH1F *hPairPt[nFile+1][gNTrgSetup];
  for(int i=0; i<nFile; i++)
    {
      //hnInvMass[i]->GetAxis(1)->SetRangeUser(3.0,3.2);
      //hnInvMass[i]->GetAxis(0)->SetRange(1,1);
      hnInvMass[i]->GetAxis(7)->SetRange(2,2);
      hnInvMass[i]->GetAxis(6)->SetRange(2,2); // tagged as dimuon event
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k>0) hnInvMass[i]->GetAxis(3)->SetRange(k, k);
	  TH1F *h1tmp = (TH1F*)hnInvMass[i]->Projection(2);
	  h1tmp->SetName(Form("hPirPt_%s%s_tmp",trig_name[i],gTrgSetupTitle[k]));
	  hPairPt[i][k] = (TH1F*)h1tmp->Rebin(nbins, Form("hPirPt_%s%s",trig_name[i],gTrgSetupTitle[k]), xbins);
	  hPairPt[i][k]->Sumw2();
	  hPairPt[i][k]->Scale(1./nMbEvtEq[i][k]);
	  hnInvMass[i]->GetAxis(3)->SetRange(0,-1);
	}
    }
  TH2F *hMtdQaPairPtEq = (TH2F*)fstudy[1]->Get("mhMtdQaPairPtEq_mb");
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k==0) h1tmp = (TH1F*)hMtdQaPairPtEq->ProjectionY(Form("hPirPt_%s%s_tmp",trig_name[2],gTrgSetupTitle[k]));
      else     h1tmp = (TH1F*)hMtdQaPairPtEq->ProjectionY(Form("hPirPt_%s%s_tmp",trig_name[2],gTrgSetupTitle[k]), k, k);
      hPairPt[2][k] = (TH1F*)h1tmp->Rebin(nbins, Form("hPirPt_%s%s",trig_name[2],gTrgSetupTitle[k]), xbins);
      hPairPt[2][k]->Sumw2();
      hPairPt[2][k]->Scale(1./nMbEvtEq[2][k]);
    }
  
  TCanvas *c = new TCanvas("cPairPtComp", "cPairPtComp", 1100, 700);
  c->Divide(3,2);
  leg = new TLegend(0.3,0.15,0.5,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int k=2; k<gNTrgSetup; k++)
    {
      c->cd(k-1);
      gPad->SetLogy();
      for(int i=0; i<3; i++)
	{
	  hPairPt[i][k]->SetMarkerStyle(18+k+4*i);
	  hPairPt[i][k]->SetMarkerColor(gColor[k]);
	  hPairPt[i][k]->SetLineColor(gColor[k]);
	  hPairPt[i][k]->SetTitle("");
	  if(i==0) hPairPt[i][k]->Draw();
	  else     hPairPt[i][k]->Draw("sames");
	  if(k==2)
	    {
	      if(i==0) leg->AddEntry(hPairPt[i][k], "Dimuon", "P");
	      if(i==1) leg->AddEntry(hPairPt[i][k], "MB+Dimuon", "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s",gLegNameTrg[k].Data()), 0.055);
      t1->Draw();

      c->cd(k+2);
      TH1F *hRatio = (TH1F*)hPairPt[0][k]->Clone(Form("%s_ratio",hPairPt[2][k]->GetName()));
      hRatio->Divide(hPairPt[1][k]);
      hRatio->GetYaxis()->SetRangeUser(0.5,1.2);
      hRatio->Draw();
      TF1 *funcRatio = new TF1(Form("%s_func",hRatio->GetName()), "pol0", 0, 8);
      hRatio->Fit(funcRatio, "IR0Q");
      funcRatio->SetLineColor(hRatio->GetMarkerColor());
      funcRatio->SetLineStyle(2);
      funcRatio->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: DM/(MB+DM)",gLegNameTrg[k].Data()), 0.055);
      t1->Draw();
      leg2 = new TLegend(0.2,0.15,0.5,0.35);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.05);
      leg2->AddEntry(funcRatio, Form("%4.3f +/- %4.3f",funcRatio->GetParameter(0),funcRatio->GetParError(0)), "L");
      leg2->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_CompMthPairPt_MbVsDm.pdf",run_type.Data()));
  
}

//================================================
void LSpairInMB(const int savePlot = 0)
{
  //**************** QA ******************
  const int iMB = 1;
  const char* trig_name[2] = {"di_mu", "mb"};
  TFile *fstudy = 0x0;
  if(iMB==0) fstudy = TFile::Open("output/Run14_AuAu200.Study.MthTrk.root", "read");
  else       fstudy = TFile::Open("output/Run14_AuAu200.MB.P15ic.Study.Dca3cm.root", "read");
  //else       fstudy = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");

  // check # of events
  double nMbEvt[gNTrgSetup];
  double nMbEvtEq[gNTrgSetup];
  double nMbEvtDm[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      nMbEvt[k]   = 0;
      nMbEvtEq[k] = 0;
      nMbEvtDm[k] = 0;
    }

  TH2F *hMbEvt = (TH2F*)fstudy->Get(Form("mhMtdQaCent_%s",trig_name[iMB]));
  nMbEvt[0] = hMbEvt->GetEntries();
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMbEvt->ProjectionY("",k,k);
      nMbEvt[k] = htmp->Integral(1,16);
    }
  TH2F *hMbEvtDm = (TH2F*)fstudy->Get(Form("mhMtdQaCentDm_%s",trig_name[iMB]));
  nMbEvtDm[0] = hMbEvtDm->GetEntries();
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMbEvtDm->ProjectionY("",k,k);
      nMbEvtDm[k] = htmp->GetEntries();
    }

  TH2F* hMtdQaCentEqMb = (TH2F*)fstudy->Get(Form("mhMtdQaCentEqMb_%s",trig_name[iMB]));
  if(iMB==0)
    {
      TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
      TH1F *hEqMb[gNTrgSetup];
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hEqMb[k] = (TH1F*)fLumi->Get(Form("hEqMbEvt_ScaleMbData%s",gTrgSetupTitle[k]));
	  nMbEvtEq[k]= hEqMb[k]->Integral(1,8);
	}
    }
  else
    {
      for(int k=1; k<gNTrgSetup; k++)
	{
	  TH1F *htmp = (TH1F*)hMtdQaCentEqMb->ProjectionY("",k,k);
	  nMbEvtEq[k] = htmp->Integral(1,16);
	  nMbEvtEq[0] += nMbEvtEq[k];
	}
    }


  // check single tracks
  const int nHistos = 4;
  const char* histoName[nHistos] = {"mhMtdQaTrkPt", "mhMtdQaTofTrkPt", "mhMtdQaMthTrkPt", "mhMtdQaMuonTrkPt"};
  TH1F *hTrkPt[nHistos][gNTrgSetup];
  for(int i=0; i<nHistos; i++)
    {
      TH2F *h2 = (TH2F*)fstudy->Get(Form("%s_mb",histoName[i]));
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k==0) hTrkPt[i][k] = (TH1F*)h2->ProjectionY(Form("%s%s",histoName[i],gTrgSetupTitle[k]), 1, 4);
	  else     hTrkPt[i][k] = (TH1F*)h2->ProjectionY(Form("%s%s",histoName[i],gTrgSetupTitle[k]), k, k);
	  hTrkPt[i][k]->Rebin(10);
	  hTrkPt[i][k]->Sumw2();
	  printf("[i] %s has %1.0f counts\n",hTrkPt[i][k]->GetName(),hTrkPt[i][k]->GetEntries());
	  hTrkPt[i][k]->Scale(1./nMbEvt[k]);
	}
    }
  TCanvas *c = new TCanvas("cTrkPt", "cTrkPt", 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2,0.15,0.8,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  leg->SetNColumns(3);
  for(int i=0; i<nHistos; i++)
    {
      c->cd(i+1);
      for(int k=0; k<gNTrgSetup; k++)
	{
	  TH1F *hRatio = (TH1F*)hTrkPt[i][k]->Clone(Form("%s_Ratio",hTrkPt[i][k]->GetName()));
	  hRatio->Divide(hTrkPt[i][2]);
	  hRatio->SetMarkerStyle(20+k);
	  hRatio->SetMarkerColor(gColor[k]);
	  hRatio->SetLineColor(gColor[k]);
	  hRatio->GetYaxis()->SetRangeUser(0.5, 1.05);
	  hRatio->SetTitle("");
	  if(k==0) hRatio->Draw("PE");
	  else hRatio->Draw("samesPE");
	  if(i==0) leg->AddEntry(hRatio, legNameLumi[k].Data(), "PL");
	}
      TPaveText *t1 = GetTitleText(Form("Run14_AuAu200: %s (0-80%%)",histoName[i]),0.045);
      t1->Draw();
      leg->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%s_CompTrkPt.pdf",run_type.Data(),trig_name[iMB]));
  
  // obtain # of pairs
  const int nbins = 7;
  double xbins[8] = {0, 1, 2, 3, 4, 5, 6, 8};

  //Type vs. M_{#mu#mu} vs. p_{T,#mumu} vs. trigSetup vs. centrality vs. isMuon vs. isDm vs. isTrig vs. isMB vs. NHitsFit_{min}
  THnSparseF *hnInvMass = (THnSparseF*)fstudy->Get(Form("mhMtdQaInvMass_%s",trig_name[iMB]));
  hnInvMass->GetAxis(4)->SetRange(1,16); // centrality
  hnInvMass->GetAxis(1)->SetRangeUser(2.6,3.6);
  if(iMB==0)
    {
      hnInvMass->GetAxis(1)->SetRangeUser(3.0,3.2);
      hnInvMass->GetAxis(0)->SetRange(1,1);
    }
  TH1F *hPirPt[4][gNTrgSetup];
  const char* pairType[4] = {"Matched", "MthTrig", "Muon", "MuonTrig"};
  for(int k=0; k<gNTrgSetup; k++)
    {
      if(k>0) hnInvMass->GetAxis(3)->SetRange(k, k);
      for(int j=0; j<4; j++)
	{
	  if(j==1) 
	    {
	      hnInvMass->GetAxis(7)->SetRange(2,2);
	    }
	  if(j==2)
	    {
	      hnInvMass->GetAxis(5)->SetRange(2,2);
	    }
	  if(j==3)
	    {
	      hnInvMass->GetAxis(5)->SetRange(2,2);
	      hnInvMass->GetAxis(7)->SetRange(2,2);
	    }
	  TH1F *h1tmp = (TH1F*)hnInvMass->Projection(2);
	  h1tmp->SetName(Form("h%sPirPt_%s_tmp",pairType[j],gTrgSetupTitle[k]));
	  hPirPt[j][k] = (TH1F*)h1tmp->Rebin(nbins, Form("h%sPirPt_%s",pairType[j],gTrgSetupTitle[k]), xbins);
	  hnInvMass->GetAxis(5)->SetRange(0,-1);
	  hnInvMass->GetAxis(6)->SetRange(0,-1);
	  hnInvMass->GetAxis(7)->SetRange(0,-1);
	}
      hnInvMass->GetAxis(3)->SetRange(0,-1);
    }

   for(int k=0; k<gNTrgSetup; k++)
     {
       printf("[i] %s: %1.0f counts in %1.0f MB events or %1.0f DM events (%4.3f%%) or %1.0f eq. MB events\n",gTrgSetupTitle[k],hPirPt[1][k]->GetEntries(),nMbEvt[k],nMbEvtDm[k],nMbEvtDm[k]/nMbEvtEq[k]*100,nMbEvtEq[k]);
     }

  // PID efficiency on LS pair
  TCanvas *c = new TCanvas("cPidEff_Pair", "cPidEff_Pair", 800, 600);
  TLegend *leg = new TLegend(0.2,0.15,0.45,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  for(int k=0; k<gNTrgSetup; k++)
    {
      TH1F *hPidEff = (TH1F*)hPirPt[3][k]->Clone(Form("hPidEff_%s",gTrgSetupTitle[k]));
      hPidEff->Sumw2();
      hPidEff->Divide(hPirPt[1][k]);
      hPidEff->SetMarkerStyle(20+k);
      hPidEff->SetMarkerColor(gColor[k]);
      hPidEff->SetLineColor(gColor[k]);
      hPidEff->GetYaxis()->SetRangeUser(0, 1);
      if(iMB==1) hPidEff->GetYaxis()->SetRangeUser(0, 1.2);
      hPidEff->SetTitle("");
      if(k==0) hPidEff->Draw("PE");
      else hPidEff->Draw("samesPE");
      leg->AddEntry(hPidEff, legNameLumi[k].Data(), "PL");
    }
  TPaveText *t1 = GetTitleText(Form("Run14_AuAu200_%s: PID efficiency (0-80%%)",trig_name[iMB]),0.04);
  t1->Draw();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%s_PidEff.pdf",run_type.Data(),trig_name[iMB]));


   // TPC tracking efficiency
  TFile *fEff = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"read");
  printf("[w] efficiency name: %s\n",fEff->GetName());
  TH1F *hJpsiEff[gNTrgSetup];
  for(int k=1; k<gNTrgSetup; k++)
    {
      hJpsiEff[k] = (TH1F*)fEff->Get(Form("LSEffVsPt_cent%s%s_final",cent_Title[0],gTrgSetupTitle[k]));
      if(k==1) 
	{
	  hJpsiEff[0] = (TH1F*)hJpsiEff[k]->Clone(Form("LSEffVsPt_cent%s%s_final",cent_Title[0],gTrgSetupTitle[0]));
	  hJpsiEff[0]->Reset("AC");
	}
      hJpsiEff[0]->Add(hJpsiEff[k], nMbEvtEq[k]/nMbEvtEq[0]);
    }
  TCanvas *c = new TCanvas("cTpcEff_Pair", "cTpcEff_Pair", 800, 600);
  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 
  for(int k=0; k<gNTrgSetup; k++)
    {
      TH1F *hRatio = (TH1F*)hJpsiEff[k]->Clone(Form("hTpcEffRatio_%s",gTrgSetupTitle[k]));
      hRatio->Divide(hJpsiEff[2]);
      hRatio->SetMarkerStyle(20+k);
      hRatio->SetMarkerColor(gColor[k]);
      hRatio->SetLineColor(gColor[k]);
      hRatio->GetYaxis()->SetRangeUser(0, 1);
      hRatio->GetYaxis()->SetRangeUser(0.5, 1.1);
      hRatio->SetTitle("");
      if(k==0) hRatio->Draw("PE");
      else hRatio->Draw("samesPE");
      leg->AddEntry(hRatio, legNameLumi[k].Data(), "PL");
    }
  TPaveText *t1 = GetTitleText(Form("Run14_AuAu200: TPC efficiency (0-80%%)"),0.045);
  t1->Draw();
  leg->Draw();

  // Acceptance correction
  TFile *fstudy = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiMb.root", "read");
  TH1F *hAccCorr[gNTrgSetup];
  for(int k=0; k<gNTrgSetup; k++)
    {
      hAccCorr[k] = (TH1F*)fstudy->Get(Form("hAccCorr%s",gTrgSetupTitle[k]));
      if(k==1) 
	{
	  hAccCorr[0] = (TH1F*)hAccCorr[k]->Clone(Form("hAccCorr%s",gTrgSetupTitle[0]));
	  hAccCorr[0]->Reset("AC");
	}
      if(iMB==0)  hAccCorr[0]->Add(hAccCorr[k], nMbEvtEq[k]/nMbEvtEq[0]);
      else        hAccCorr[0]->Add(hAccCorr[k], nMbEvt[k]/nMbEvt[0]);
    }

  // LS Trigger efficiency
  TFile *fTrigEffLS = TFile::Open(Form("Rootfiles/Run14_AuAu200.StudyTrigEff.root"),"read");
  TH1F *hJpsiTrigEffLS[gNTrgSetup];
  for(int k=1; k<gNTrgSetup; k++)
    {
      hJpsiTrigEffLS[k] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_M1_TrigStudy",gTrgSetupTitle[k]));
      if(k==1) 
	{
	  hJpsiTrigEffLS[0] = (TH1F*)hJpsiTrigEffLS[k]->Clone(Form("Run14_AuAu200_JpsiEffVsPt_TacDiffEff%s_M1_TrigStudy",gTrgSetupTitle[0]));
	  hJpsiTrigEffLS[0]->Reset("AC");
	}
      if(iMB==0) hJpsiTrigEffLS[0]->Add(hJpsiTrigEffLS[k], nMbEvtEq[k]/nMbEvtEq[0]);
      else       hJpsiTrigEffLS[0]->Add(hJpsiTrigEffLS[k], nMbEvt[k]/nMbEvt[0]);
    }

  // Matched track trigger efficiency
  TH1F *hJpsiTrigEffMth[gNTrgSetup];
  for(int k=1; k<gNTrgSetup; k++)
    {
      hJpsiTrigEffMth[k] = (TH1F*)fTrigEffLS->Get(Form("Run14_AuAu200_MtdTrig_JpsiEffVsPt%s_M2",gTrgSetupTitle[k]));
      if(k==1) 
	{
	  hJpsiTrigEffMth[0] = (TH1F*)hJpsiTrigEffMth[k]->Clone(Form("Run14_AuAu200_MtdTrig_JpsiEffVsPt%s_M2",gTrgSetupTitle[0]));
	  hJpsiTrigEffMth[0]->Reset("AC");
	}
      if(iMB==0) hJpsiTrigEffMth[0]->Add(hJpsiTrigEffMth[k], nMbEvtEq[k]/nMbEvtEq[0]);
      else       hJpsiTrigEffMth[0]->Add(hJpsiTrigEffMth[k], nMbEvtEq[k]/nMbEvtEq[0]);
    }


  // plot invariant yield
  for(int j=0; j<4; j++)
     {
       printf("+++++ %s +++++\n",pairType[j]);
       TCanvas *c = new TCanvas(Form("c%sPair",pairType[j]), Form("c%sPair",pairType[j]), 1100, 500);
       c->Divide(2,1);

       c->cd(1);
       gPad->SetLogy();
       SetPadMargin(gPad,0.13,0.13,0.05,0.1);
       TLegend *leg = new TLegend(0.2,0.15,0.45,0.4);
       leg->SetBorderSize(0);
       leg->SetFillColor(0);
       leg->SetTextSize(0.04); 
       for(int k=0; k<gNTrgSetup; k++)
	 {
	   hPirPt[j][k]->Sumw2();

	   if(iMB==0) 
	     hPirPt[j][k]->Scale(1./nMbEvtEq[k]);
	   else
	     hPirPt[j][k]->Scale(1./nMbEvt[k]);
	   /*
	   hPirPt[j][k]->Scale(1./nMbEvtDm[k]);
	   */	   
	   double acc = 1, trig_eff = 1, tpc_eff = 1;
	   for(int bin=1; bin<=hPirPt[j][k]->GetNbinsX(); bin++)
	     {
	       //if(iMB==0 || (iMB==1 && (j==1 || j==3))) 
		 acc = hAccCorr[k]->GetBinContent(bin);
		 //else 
		 //acc = 1;

	       if(j==1)
		 trig_eff = hJpsiTrigEffMth[k]->GetBinContent(bin);
	       else if(j==3)
		 trig_eff = hJpsiTrigEffLS[k]->GetBinContent(bin);
	       else
		 trig_eff = 1;

	       tpc_eff = hJpsiEff[k]->GetBinContent(bin);

	       hPirPt[j][k]->SetBinContent(bin, hPirPt[j][k]->GetBinContent(bin)/tpc_eff/trig_eff/acc);
	       hPirPt[j][k]->SetBinError(bin, hPirPt[j][k]->GetBinError(bin)/tpc_eff/trig_eff/acc);
	     }
	   hPirPt[j][k]->SetMarkerStyle(21);
	   hPirPt[j][k]->SetMarkerColor(gColor[k]);
	   hPirPt[j][k]->SetLineColor(gColor[k]);
	   hPirPt[j][k]->SetMarkerSize(1.2);
	   hPirPt[j][k]->SetTitle(";p_{T} (GeV/c);1/N dN/dp_{T}");
	   ScaleHistoTitle(hPirPt[j][k],0.045,1,0.04,0.045,1.4,0.04,62);
	   if(k==0) hPirPt[j][k]->Draw();
	   else     hPirPt[j][k]->Draw("sames");
	   leg->AddEntry(hPirPt[j][k],legNameLumi[k].Data(),"P");
	 }
       leg->Draw();
       TPaveText *t1 = GetTitleText(Form("Run14_AuAu200_%s: %s pairs (0-80%%)",trig_name[iMB],pairType[j]),0.04);
       t1->Draw();

       c->cd(2);
       SetPadMargin(gPad,0.13,0.13,0.05,0.1);
       leg = new TLegend(0.2,0.15,0.8,0.26);
       leg->SetBorderSize(0);
       leg->SetFillColor(0);
       leg->SetTextSize(0.04); 
       leg->SetNColumns(2);
       for(int k=0; k<gNTrgSetup; k++)
	 {
	   TH1F *hRatio = (TH1F*)hPirPt[j][k]->Clone(Form("%s_clone",hPirPt[j][k]->GetName()));
	   for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
	     {
	       if(hPirPt[j][0]->GetBinContent(bin)<=0) continue;
	       hRatio->SetBinContent(bin, hPirPt[j][k]->GetBinContent(bin)/hPirPt[j][0]->GetBinContent(bin));
	       hRatio->SetBinError(bin, hPirPt[j][k]->GetBinError(bin)/hPirPt[j][0]->GetBinContent(bin));
	     }
	   hRatio->GetYaxis()->SetRangeUser(0.75, 1.25);
	   if(iMB==1 && j>1) hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
	   hRatio->SetTitle(";p_{T} (GeV/c);Ratio");
	   if(k==0) hRatio->Draw();
	   else     hRatio->Draw("sames");
	   TF1 *func = new TF1(Form("%s_fit",hRatio->GetName()), "pol0", 0, 8);
	   hRatio->Fit(func, "IR0Q");
	   func->SetLineColor(gColor[k]);
	   func->SetLineStyle(2);
	   func->Draw("sames");
	   printf("[i] %s: ratio = %4.2f +/- %4.2f\n",legNameLumi[k].Data(),func->GetParameter(0),func->GetParError(0));
	   if(k>0) leg->AddEntry(func, Form("%2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)), "L");
	 }
       t1 = GetTitleText(Form("Run14_AuAu200_%s: ratio to inclusive",trig_name[iMB]),0.04);
       t1->Draw();
       leg->Draw();
       if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%s_%sYield.pdf",run_type.Data(),trig_name[iMB],pairType[j]));
     }
}

//================================================
void JpsiTree(const int savePlot = 0)
{
  TFile *ftree[2];
  TFile *fstudy[2];
  
  ftree[0] = TFile::Open("output/Run14_AuAu200.MB.JpsiTree.root", "read");
  ftree[1] = TFile::Open("output/Run14_AuAu200.JpsiTree.root", "read");

  fstudy[0] = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");
  fstudy[1] = TFile::Open("output/Run14_AuAu200.Study.root", "read");

  // event statistics
  const char* trig_name[2] = {"mb", "di_mu"};
  TH2F *hMtdQaCent[2];
  TH2F *hMtdQaCentDm[2];
  TH2F *hMtdQaCentEqMb[2];
  for(int i=0; i<2; i++)
    {
      hMtdQaCent[i] = (TH2F*)fstudy[i]->Get(Form("mhMtdQaCent_%s",trig_name[i]));
      hMtdQaCentDm[i] = (TH2F*)fstudy[i]->Get(Form("mhMtdQaCentDm_%s",trig_name[i]));
      hMtdQaCentEqMb[i] = (TH2F*)fstudy[i]->Get(Form("mhMtdQaCentEqMb_%s",trig_name[i]));
    }

  double nDmEvtInMb[gNTrgSetup];
  nDmEvtInMb[0] = 0;
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMtdQaCentDm[0]->ProjectionY("",k,k);
      nDmEvtInMb[k] = htmp->Integral(1,16);
      nDmEvtInMb[0] += nDmEvtInMb[k];
    }

  double nEqMbForDmInMb[gNTrgSetup];
  nEqMbForDmInMb[0] = 0;
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMtdQaCentEqMb[0]->ProjectionY("",k,k);
      nEqMbForDmInMb[k] = htmp->Integral(1,16);
      nEqMbForDmInMb[0] += nEqMbForDmInMb[k];
    }

  double nDmEvtInDm[gNTrgSetup];
  nDmEvtInDm[0] = 0;
  for(int k=1; k<gNTrgSetup; k++)
    {
      TH1F *htmp = (TH1F*)hMtdQaCentDm[1]->ProjectionY("",k,k);
      nDmEvtInDm[k] = htmp->Integral(1,16);
      nDmEvtInDm[0] += nDmEvtInDm[k];
    }


  TH1F *hPairPtAll[2][gNTrgSetup];
  TH1F *hPairPtAcc[2][gNTrgSetup];
  for(int i=0; i<2; i++)
    {
      for(int k=0; k<gNTrgSetup; k++)
	{
	  hPairPtAll[i][k] = new TH1F(Form("hPairPtAll_%d%s",i,gTrgSetupTitle[k]), ";p_{T} [GeV/c];counts", 10, 0, 10);
	  hPairPtAcc[i][k] = new TH1F(Form("hPairPtAcc_%d%s",i,gTrgSetupTitle[k]), ";p_{T} [GeV/c];counts", 10, 0, 10);
	}
    }
	
  TH1F *hNMtdHitsAll[2];
  TH1F *hNMtdHitsOlp[2];
  TH1F *hNMthTrksAll[2];
  TH1F *hNMthTrksOlp[2];
  for(int i=0; i<2; i++)
    {
      hNMtdHitsAll[i] = new TH1F(Form("hNMtdHitsAll_%s",trig_name[i]), ";N;", 20, 0, 20);
      hNMtdHitsOlp[i] = new TH1F(Form("hNMtdHitsOlp_%s",trig_name[i]), ";N;", 20, 0, 20);
      hNMthTrksAll[i] = new TH1F(Form("hNMthTrksAll_%s",trig_name[i]), ";N;", 20, 0, 20);
      hNMthTrksOlp[i] = new TH1F(Form("hNMthTrksOlp_%s",trig_name[i]), ";N;", 20, 0, 20); 
    }
  double nDmEvtInMbDm[gNTrgSetup] = {0, 0, 0, 0, 0};

  int runId[2], evtId[2], trigSetup[2], centrality[2]; float zdcRate[2];
  int nMthTrks[2]; float trackPt[2][100];
  int nMtdHits[2], backleg[2][100], module[2][100], cell[2][100], isTrig[2][100];
  int nPairs[2], type[2][100], isMtdTrig[2][100], isMuonPid[2][100];
  float pairM[2][100], pairPt[2][100], pairEta[2][100], pairPhi[2][100];
  TTree *jpsiTree[2];
  for(int i=0; i<2; i++)
    {
      jpsiTree[i] = (TTree*)ftree[i]->Get("mOutJpsiTree");
      jpsiTree[i]->SetName(Form("%s_%s",jpsiTree[i]->GetName(),trig_name[i]));

      jpsiTree[i]->SetBranchAddress("runId", &runId[i]);
      jpsiTree[i]->SetBranchAddress("evtId", &evtId[i]);
      jpsiTree[i]->SetBranchAddress("trigSetup", &trigSetup[i]);
      jpsiTree[i]->SetBranchAddress("centrality", &centrality[i]);
      jpsiTree[i]->SetBranchAddress("zdcRate", &zdcRate[i]);

      jpsiTree[i]->SetBranchAddress("nMthTrks", &nMthTrks[i]);
      jpsiTree[i]->SetBranchAddress("trackPt", &trackPt[i]);

      jpsiTree[i]->SetBranchAddress("nMtdHits", &nMtdHits[i]);
      jpsiTree[i]->SetBranchAddress("backleg", &backleg[i]);
      jpsiTree[i]->SetBranchAddress("module", &module[i]);
      jpsiTree[i]->SetBranchAddress("cell", &cell[i]);
      jpsiTree[i]->SetBranchAddress("isTrig", &isTrig[i]);

      jpsiTree[i]->SetBranchAddress("nPairs", &nPairs[i]);
      jpsiTree[i]->SetBranchAddress("type", &type[i]);
      jpsiTree[i]->SetBranchAddress("isMtdTrig", &isMtdTrig[i]);
      jpsiTree[i]->SetBranchAddress("isMuonPid", &isMuonPid[i]);
      jpsiTree[i]->SetBranchAddress("pairM", &pairM[i]);
      jpsiTree[i]->SetBranchAddress("pairPt", &pairPt[i]);
      jpsiTree[i]->SetBranchAddress("pairEta", &pairEta[i]);
      jpsiTree[i]->SetBranchAddress("pairPhi", &pairPhi[i]);
    }
  const int nMbEvts = jpsiTree[0]->GetEntries();
  const int nDmEvts = jpsiTree[1]->GetEntries();
  cout << nMbEvts << "  " << nDmEvts << endl;

  // build a map for accepted event Id
  TH1F *hRunId = new TH1F("hRunId", "", 94000, 15074000, 15168000);
  jpsiTree[1]->Draw("runId >> hRunId", "", "");
  map<int,int> runIndex;
  int counter = 0;
  for(int bin=1; bin<=94000; bin++)
    {
      if(hRunId->GetBinContent(bin)<=0) continue;
      runIndex[15074000+bin-1] = counter;
      counter++;
    }
  TFile *fdm = TFile::Open("output/Run14_AuAu200.jpsi.root");	
  TH1F *hDmEvtAcc = (TH1F*)fdm->Get("mhEvtRunAcc_di_mu");
  hRunId->Sumw2();
  for(int bin=1; bin<=hRunId->GetNbinsX(); bin++)
    {
      int jbin = hDmEvtAcc->FindBin(hRunId->GetBinCenter(bin)-0.4);
      double dmevt = hDmEvtAcc->GetBinContent(jbin);
      if(hRunId->GetBinContent(bin)<=0) continue;
      hRunId->SetBinContent(bin, hRunId->GetBinContent(bin)/dmevt);
      hRunId->SetBinError(bin, hRunId->GetBinError(bin)/dmevt);
    }
  hRunId->SetMarkerStyle(24);
  draw1D(hRunId);
  cout << hRunId->GetEntries() << endl;
  return;

  int evtIds[3000][250];
  int evtCounter[3000];
  for(int i=0; i<3000; i++) evtCounter[i] = 0;
  for(int e=0; e<nDmEvts; e++)
    {
      jpsiTree[1]->GetEntry(e);
      evtIds[runIndex[runId[1]]][evtCounter[runIndex[runId[1]]]] = evtId[1];
      evtCounter[runIndex[runId[1]]]++;
      hNMtdHitsAll[1]->Fill(nMtdHits[1]);
      hNMthTrksAll[1]->Fill(nMthTrks[1]);
      hNMtdHitsOlp[1]->Fill(nMtdHits[1]);
      hNMthTrksOlp[1]->Fill(nMthTrks[1]);
    }

  for(int e=0; e<nMbEvts; e++)
    {
      jpsiTree[0]->GetEntry(e);
      if(runId[0]<15077036) continue;

      hNMtdHitsAll[0]->Fill(nMtdHits[0]);
      hNMthTrksAll[0]->Fill(nMthTrks[0]);
      for(int p=0; p<nPairs[0]; p++)
	{
	  if(isMtdTrig[0][p]) 
	    {
	      hPairPtAll[0][0]->Fill(pairPt[0][p]);
	      hPairPtAll[0][trigSetup[0]+1]->Fill(pairM[0][p]);
	    }
	  if(isMtdTrig[0][p]&&isMuonPid[0][p]) 
	    {
	      hPairPtAll[1][0]->Fill(pairPt[0][p]);
	      hPairPtAll[1][trigSetup[0]+1]->Fill(pairM[0][p]);
	    }
	} 
      bool isEvtInDm = false;
      int index = runIndex[runId[0]];
      for(int j=0; j<evtCounter[index]; j++)
	{
	  if(evtId[0]==evtIds[index][j])
	    {
	      isEvtInDm = true;
	      break;
	    }
	}

      if(!isEvtInDm) continue;
      nDmEvtInMbDm[0]++;
      nDmEvtInMbDm[trigSetup[0]+1]++;

      hNMtdHitsOlp[0]->Fill(nMtdHits[0]);
      hNMthTrksOlp[0]->Fill(nMthTrks[0]);
      for(int p=0; p<nPairs[0]; p++)
	{
	  if(isMtdTrig[0][p]) 
	    {
	      hPairPtAcc[0][0]->Fill(pairPt[0][p]);
	      hPairPtAcc[0][trigSetup[0]+1]->Fill(pairM[0][p]);
	    }
	  if(isMtdTrig[0][p]&&isMuonPid[0][p]) 
	    {
	      hPairPtAcc[1][0]->Fill(pairPt[0][p]);
	      hPairPtAcc[1][trigSetup[0]+1]->Fill(pairM[0][p]);
	    }
	}
    }

  // check NMtdHits and NMtdTrks
  TCanvas *cMtdHits = new TCanvas("cMtdHits", "cMtdHits", 800, 600);
  gPad->SetLogy();
  for(int i=0; i<2; i++)
    {
      hNMtdHitsAll[i]->Sumw2();
      hNMtdHitsAll[i]->SetMarkerStyle(20+i*4);
      hNMtdHitsAll[i]->SetMarkerColor(gColor[i]);
      hNMtdHitsAll[i]->SetLineColor(gColor[i]);
      if(i==0) hNMtdHitsAll[i]->Draw();
      else hNMtdHitsAll[i]->Draw("samesPE");

      if(i==0)
	{
	  hNMtdHitsOlp[i]->Sumw2();
	  hNMtdHitsOlp[i]->SetMarkerStyle(21+i*4);
	  hNMtdHitsOlp[i]->SetMarkerColor(gColor[i+2]);
	  hNMtdHitsOlp[i]->SetLineColor(gColor[i+2]);
	  hNMtdHitsOlp[i]->Draw("samesPE");
	}
    }
  leg = new TLegend(0.5,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hNMtdHitsAll[0], "MB", "PE");
  leg->AddEntry(hNMtdHitsAll[1], "DM", "PE");
  leg->AddEntry(hNMtdHitsOlp[0], "MB+DM from MB", "PE");
  leg->Draw();
     
  // equivalent # of mb events
  for(int k=0; k<gNTrgSetup; k++)
    {
      printf("[i] %9s: dm_in_mb = %1.0f, dm_in_mb&dm = %1.0f (%4.4f); eq_mb_for_dm_in_mb = %4.3e, dm_in_dm = %4.3e, eq_mb_for_dm_in_dm = %4.3e\n",
	     gLegNameTrg[k].Data(),nDmEvtInMb[k],nDmEvtInMbDm[k],nDmEvtInMbDm[k]/nDmEvtInMb[k],nEqMbForDmInMb[k],nDmEvtInDm[k],nDmEvtInDm[k]/nDmEvtInMbDm[k]*nEqMbForDmInMb[k]);
    }
  
  // two-pass iteration efficiency 
  const char* typeName[2] = {"TrigMth", "TrigMuon"};
  for(int i=0; i<2; i++)
    {
      TCanvas *c = new TCanvas(Form("c%sEff",typeName[i]), Form("c%sEff",typeName[i]), 1100, 700);
      c->Divide(4,2);
      TLegend *leg = new TLegend(0.3,0.15,0.5,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      for(int k=1; k<gNTrgSetup; k++)
	{
	  c->cd(k);
	  gPad->SetLogy();
	  hPairPtAll[i][k]->Sumw2();
	  hPairPtAll[i][k]->SetMarkerStyle(25);
	  hPairPtAll[i][k]->Draw("PE");
	  hPairPtAcc[i][k]->SetMarkerStyle(21);
	  hPairPtAcc[i][k]->SetMarkerColor(gColor[k]);
	  hPairPtAcc[i][k]->SetLineColor(gColor[k]);
	  hPairPtAcc[i][k]->Draw("samesPE");
	  if(i==0) printf("[i] %s: %1.0f/%1.0f LS pairs\n",hPairPtAll[i][k]->GetName(),hPairPtAcc[i][k]->Integral(1,4),hPairPtAll[i][k]->Integral(1,4));

	  if(k==1)
	    {
	      leg->AddEntry(hPairPtAll[i][k], "MB", "PL");
	      leg->AddEntry(hPairPtAcc[i][k], "MB+DM", "PL");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s: %s",gLegNameTrg[k].Data(),typeName[i]), 0.055);
	  t1->Draw();

	  c->cd(k+4);
	  TGraphAsymmErrors *gRatio = new TGraphAsymmErrors(hPairPtAcc[i][k], hPairPtAll[i][k],"cl=0.683 b(1,1) mode");
	  gRatio->SetMarkerStyle(hPairPtAcc[i][k]->GetMarkerStyle());
	  gRatio->SetMarkerColor(hPairPtAcc[i][k]->GetMarkerColor());
	  gRatio->SetLineColor(hPairPtAcc[i][k]->GetLineColor());
	  gRatio->SetName(Form("%s_ratio",hPairPtAcc[i][k]->GetName()));
	  gRatio->GetYaxis()->SetRangeUser(0.8,1.2);
	  gRatio->Draw("AP");
	  TF1 *funcRatio = new TF1(Form("%s_func",gRatio->GetName()), "pol0", 0, 10);
	  gRatio->Fit(funcRatio, "IR0Q");
	  funcRatio->SetLineColor(gRatio->GetMarkerColor());
	  funcRatio->SetLineStyle(2);
	  funcRatio->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%s: (MB+DM)/MB",gLegNameTrg[k].Data()), 0.055);
	  t1->Draw();
	  leg2 = new TLegend(0.2,0.15,0.5,0.35);
	  leg2->SetBorderSize(0);
	  leg2->SetFillColor(0);
	  leg2->SetTextSize(0.05);
	  leg2->AddEntry(funcRatio, Form("%4.3f +/- %4.3f",funcRatio->GetParameter(0),funcRatio->GetParError(0)), "L");
	  leg2->Draw();
	}
      c->cd(1);
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_ParTrkEff_%s.pdf",run_type.Data(),typeName[i]));
    }
}


//================================================
void dataQa(const int savePlot = 1)
{
  TF1 *fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)");
  TF1 *fResDyVsPt = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)");
  fResDzVsPt->SetParameters(-21.04, 21.09, 0.693);
  fResDyVsPt->SetParameters(-12.61, 13.43, 0.889);

  //**************** QA ******************
  const int iMB = 1;
  const char* trig_name[2] = {"di_mu", "mb"};
  TFile *fstudy = 0x0;
  if(iMB==0) fstudy = TFile::Open("output/Run14_AuAu200.Study.MthTrk.root", "read");
  else       fstudy = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");
  const char* hName[13] = {"mhMtdQaMthTrkDca", "mhMtdQaMthTrkNHitsFit", "mhMtdQaNMthTrk", "mhMtdQaMthTrkPt",
			   "mhMtdQaMthTrkEtaPhi", "mhMtdQaMthTrkNHitsDedx", "mhMtdQaMthTrkNHitsFrac", "mhMtdQaMthTrkNSigmaPi",
			   "mhMtdQaNMtdHit", "mhMtdQaHitMap", "mhMtdQaMthDy", "mhMtdQaMthDz", "mhMtdQaMthDtof"};
  const int nPtCuts = 3;
  const double ptCuts[nPtCuts+1] = {1, 2, 3, 5};
  for(int i=0; i<13; i++)
    {
      TString saveName = hName[i];
      saveName.ReplaceAll("mhMtdQa","");
      cout << saveName.Data() << endl;
      if(i==2 || i==3 || i==8)
	{
	  TH2F *h2 = 0x0;
	  if(i==2 || i==3)
	    {
	      h2 = (TH2F*)fstudy->Get(Form("%s_%s",hName[i],trig_name[iMB]));
	    }
	  else
	    {
	      if(iMB==0) h2 = (TH2F*)fstudy->Get(Form("%s_%s",hName[i],trig_name[iMB]));
	      else
		{
		  TH3F *h3 = (TH3F*)fstudy->Get(Form("%s_%s",hName[i],trig_name[iMB]));
		  h2 = (TH2F*)h3->Project3D("yx");
		}
	    }

	  TCanvas *c = new TCanvas(Form("c%s",saveName.Data()), Form("c%s",saveName.Data()), 800, 600);
	  gPad->SetLogy();
	  TLegend *leg = new TLegend(0.5,0.6,0.65,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  leg->SetHeader(run_type.Data());
	  for(int k=0; k<4; k++)
	    {
	      TH1F *h1proj = (TH1F*)h2->ProjectionY(Form("%s_proj%d",h2->GetName(),k),k+1,k+1);
	      if(h1proj->GetEntries()>0) h1proj->Scale(1./h1proj->GetEntries());
	      h1proj->SetMarkerStyle(20+k);
	      h1proj->SetMarkerColor(gColor[k+1]);
	      h1proj->SetLineColor(gColor[k+1]);
	      h1proj->SetTitle("");
	      if(k==0) h1proj->Draw();
	      else     h1proj->Draw("sames");
	      leg->AddEntry(h1proj, Form("%s: mean = %4.2f",gTrgSetupTitle[k+1],h1proj->GetMean()), "L");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s",h2->GetTitle()),0.045);
	  t1->Draw();
	  leg->Draw();
	}
      else 
	{
	  TH3F *h3 = (TH3F*)fstudy->Get(Form("%s_%s",hName[i],trig_name[iMB]));
	  TCanvas *c = new TCanvas(Form("c%s",saveName.Data()), Form("c%s",saveName.Data()), 1100, 700);
	  c->Divide(2,2);
	  if(i==4 || i==9)
	    {
	      for(int k=0; k<4; k++)
		{
		  h3->GetXaxis()->SetRange(k+1,k+1);
		  TH2F *h2proj = (TH2F*)h3->Project3D(Form("proj%d_zy",k));
		  c->cd(k+1);
		  h2proj->SetTitle("");
		  ScaleHistoTitle(h2proj, 0.05, 0.9, 0.045, 0.05, 0.9, 0.045);
		  h2proj->Draw("colz");
		  TPaveText *t1 = GetTitleText(Form("%s (%s)",h3->GetTitle(),gTrgSetupTitle[k+1]),0.055);
		  t1->Draw();
		}
	    }
	  else
	    {
	      TLegend *leg = new TLegend(0.25,0.4,0.5,0.85);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.06);
	      leg->SetHeader(run_type.Data());
	      for(int j=0; j<nPtCuts; j++)
		{
		  TLegend *leg2 = new TLegend(0.7,0.5,0.85,0.85);
		  leg2->SetBorderSize(0);
		  leg2->SetFillColor(0);
		  leg2->SetTextSize(0.055);
		  leg2->SetHeader("efficiency");
		  int low_pt_bin = h3->GetZaxis()->FindBin(ptCuts[j]+1e-4);
		  int high_pt_bin = h3->GetZaxis()->FindBin(ptCuts[j+1]-1e-4);
		  c->cd(j+1);
		  if(i==0) gPad->SetLogy();
		  double min, max;
		  if(i==7 || i==10 || i==11 || i==12)
		    {
		      if(i==7)  { min = -1; max = 3; }
		      if(i==10) 
			{
			  if(j<2) 
			    { 
			      min = -2 * fResDyVsPt->Eval((ptCuts[j]+ptCuts[j+1])/2); 
			      max = -1 * min;
			    }
			  else
			    {
			      min = -2.5 * fResDyVsPt->Eval((ptCuts[j]+ptCuts[j+1])/2); 
			      max = -1 * min;
			    }
			}
		      if(i==11)
			{
			  if(j<2) 
			    { 
			      min = -2 * fResDzVsPt->Eval((ptCuts[j]+ptCuts[j+1])/2); 
			      max = -1 * min;
			    }
			  else
			    {
			      min = -2.5 * fResDzVsPt->Eval((ptCuts[j]+ptCuts[j+1])/2); 
			      max = -1 * min;
			    }
			}
		      if(i==12)
			{
			  min = -1e4;
			  max = 0.75;
			}
		    }
		  for(int k=0; k<4; k++)
		    {
		      TH1F *h1proj = (TH1F*)h3->ProjectionY(Form("%s_pt%d_%d",h3->GetName(),j,k),k+1,k+1,low_pt_bin,high_pt_bin);
		      h1proj->Sumw2();
		      if(h1proj->GetEntries()>0) h1proj->Scale(1./h1proj->Integral(0,-1));
		      h1proj->SetMaximum(1.2*h1proj->GetMaximum());
		      h1proj->SetMarkerStyle(20+k);
		      h1proj->SetMarkerColor(gColor[k+1]);
		      h1proj->SetLineColor(gColor[k+1]);
		      h1proj->SetTitle("");
		      if(i==0) h1proj->GetXaxis()->SetRangeUser(0,3);
		      if(k==0) h1proj->Draw("PE");
		      else     h1proj->Draw("samesPE");
		      if(j==0) leg->AddEntry(h1proj, gTrgSetupTitle[k+1], "PL");
		      if(i==7 || i==10 || i==11 || i==12)
			{
			  double frac = h1proj->Integral(h1proj->FindBin(min), h1proj->FindBin(max));
			  leg2->AddEntry(h1proj, Form("%4.2f%%",frac*100), "P");
			}
		    }
		  TPaveText *t1 = GetTitleText(Form("%s (%1.0f < p_{T} < %1.0f GeV/c)",h3->GetTitle(),ptCuts[j],ptCuts[j+1]),0.055);
		  t1->Draw();
		  if(i==7 || i==10 || i==11 || i==12)
		    {
		      TLine *line1 = GetLine(min, 0, min, h1proj->GetMaximum(), 2, 1);
		      line1->Draw();
		      TLine *line2 = GetLine(max, 0, max, h1proj->GetMaximum(), 2, 1);
		      line2->Draw();
		      leg2->Draw();
		    }
		}
	      c->cd(4);
	      leg->Draw();
	    }
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%s_Qa%s.pdf",run_type.Data(),trig_name[iMB],saveName.Data()));
    }

  if(iMB==1)
    {
      const int nHistos2 = 2;
      const char* histoName2[nHistos2] = {"mhMtdQaNgTrk", "mhMtdQaNMtdHit"};
      TProfile *hPro[nHistos2][gNTrgSetup];
      TH2F *h2 = NULL;
      for(int i=0; i<nHistos2; i++)
	{
	  TH3F *h3 = (TH3F*)fstudy->Get(Form("%s_%s",histoName2[i],trig_name[iMB]));
	  for(int k=0; k<gNTrgSetup; k++)
	    {
	      if(k>0) h3->GetXaxis()->SetRange(k,k);
	      h2 = (TH2F*)h3->Project3D("zy");
	      h2->SetName(Form("%s_%d",h2->GetName(),k));
	      hPro[i][k] = (TProfile*)h2->ProfileY(Form("%s_pro",h2->GetName()));
	    }

	  TCanvas *c = new TCanvas(Form("c%s",histoName2[i]), Form("c%s",histoName2[i]), 800, 600);
	  TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04); 
	  for(int k=2; k<gNTrgSetup; k++)
	    {
	      hPro[i][k]->SetMarkerStyle(20+k);
	      hPro[i][k]->SetMarkerColor(gColor[k]);
	      hPro[i][k]->SetLineColor(gColor[k]);
	      hPro[i][k]->SetTitle("");
	      if(i==0) hPro[i][k]->GetYaxis()->SetRangeUser(0,2500);
	      if(i==1) hPro[i][k]->GetYaxis()->SetRangeUser(1,2.2);
	      if(k==2) hPro[i][k]->Draw("PE");
	      else hPro[i][k]->Draw("samesPE");
	      leg->AddEntry(hPro[i][k], legNameLumi[k].Data(), "PL");
	    }
	  TPaveText *t1 = NULL;
	  if(i==0) t1 = GetTitleText(Form("Run14_AuAu200_MB: # of global tracks (0-80%%)"),0.04);
	  if(i==1) t1 = GetTitleText(Form("Run14_AuAu200_MB: # of good MTD hits (0-80%%)"),0.04);
	  t1->Draw();
	  leg->Draw();
	  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Xsec_%s_%s.pdf",run_type.Data(),trig_name[iMB],histoName2[i]));
	}
    }
}

//================================================
void makeHistos(const int saveHisto = 1)
{
  // data
  const int nCentBins      = nCentBins_pt; 
  TFile *fin = TFile::Open(Form("output/Run14_AuAu200.jpsi.%sroot",run_config), "read");
  THnSparseF *hnInvMass[2] = {0x0};
  hnInvMass[0] = (THnSparseF*)fin->Get("mhJpsiInfo_di_mu");
  hnInvMass[1] = (THnSparseF*)fin->Get("mhBkgLSPos_di_mu");
  hnInvMass[1]->Add((THnSparseF*)fin->Get("mhBkgLSNeg_di_mu"));

  TH1F *hPairPt[nCentBins][gNTrgSetup][2] = {0x0};
  for(Int_t j=0; j<2; j++) // pair type
    { 
      hnInvMass[j]->GetAxis(0)->SetRangeUser(3.0, 3.2); // mass cut
      hnInvMass[j]->GetAxis(2)->SetRangeUser(pt1_cut+0.01,100);
      hnInvMass[j]->GetAxis(3)->SetRangeUser(pt2_cut+0.01,100);

      for(int i=0; i<nCentBins; i++)
	{
	  hnInvMass[j]->GetAxis(4)->SetRange(centBinsLow[i],centBinsHigh[i]);
	  for(int t=0; t<gNTrgSetup; t++) // trigger setup
	    {
	      if(t>0) hnInvMass[j]->GetAxis(5)->SetRange(t,t);
	      TH1F *htmp = (TH1F*)hnInvMass[j]->Projection(1);
	      htmp->SetName(Form("htmp_%d%d",j,t));
	      htmp->Sumw2();
	      hPairPt[i][t][j] = (TH1F*)htmp->Rebin(nPtBins, Form("%s_PairPt_%s_cent%s%s",run_config,gPairName[j],cent_Title[i],gTrgSetupTitle[t]), xPtBins);
	      hPairPt[i][t][j]->SetTitle();
	      hnInvMass[j]->GetAxis(5)->SetRange(0,-1);
	      printf("[i] %s has %1.0f counts\n", hPairPt[i][t][j]->GetName(), hPairPt[i][t][j]->Integral());
	    }
	  hnInvMass[j]->GetAxis(4)->SetRange(0,-1);
	}
      hnInvMass[j]->GetAxis(0)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(2)->SetRange(0,-1);
      hnInvMass[j]->GetAxis(3)->SetRange(0,-1);
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root", "update");
      for(Int_t j=0; j<2; j++) // pair type
	{ 
	  for(int i=0; i<nCentBins; i++)
	    {
	      for(int t=0; t<gNTrgSetup; t++) // trigger setup
		{
		  hPairPt[i][t][j]->Write("", TObject::kOverwrite);
		}
	    }
	}
    }
}
