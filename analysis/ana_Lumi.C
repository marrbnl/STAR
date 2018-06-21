const int year = YEAR;
TFile *f;

//================================================
void ana_Lumi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //makeHistoLumi();
  //makeHistoData();
  //mergeHisto();

  //Lumi2013();
  Lumi2014();
  //Lumi2015();
  //pAuCent();
  //Nevents();
}

//================================================
void pAuCent(const int savePlot = 0)
{
  run_type = "Run15_pAu200";
  f = TFile::Open("Rootfiles/Run15_pAu200.Centrality.root", "read");

  TH1F *hEff = (TH1F*)f->Get("hEffi_1");
  TF1 *funcEff = new TF1("FuncEff","[0]-exp(-1*[1]*(x-[2]))",0,20);
  hEff->Fit(funcEff,"R0");
  funcEff->SetNpx(20);
  //funcEff->SetParameter(1, 0.2);
  hEff->SetXTitle("N_{trk}");
  hEff->SetYTitle("Efficiency");
  hEff->SetMarkerStyle(21);
  hEff->GetXaxis()->SetRangeUser(0,40);
  c = draw1D(hEff);
  funcEff->SetLineColor(2);
  funcEff->SetLineStyle(2);
  funcEff->Draw("sames");
  TLegend *leg = new TLegend(0.15, 0.6, 0.35, 0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEff,"HIJING","P");
  leg->AddEntry(funcEff,"Fit: p0*e^{-p1*(x-p2)}","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_VpdEffFit.pdf",run_type));

  TH1F *hNPrimTrk = (TH1F*)f->Get("NGoodPrimaryTracks_8");
  hNPrimTrk->SetXTitle("N_{trk}");
  hNPrimTrk->SetYTitle("a.u.");
  hNPrimTrk->GetXaxis()->SetRangeUser(0,80);
  hNPrimTrk->Scale(1./hNPrimTrk->Integral());
  c = draw1D(hNPrimTrk,"",true,false);
  
  TH1F *hPrimTrkWeight = (TH1F*)hNPrimTrk->Clone("hPrimTrkWeight");
  TH1F *hWeight = (TH1F*)funcEff->GetHistogram();
  for(int bin=1; bin<=hPrimTrkWeight->GetNbinsX(); bin++)
    {
      double weight = hWeight->GetBinContent(bin)/hWeight->GetBinContent(20);
      if(bin>20) weight = 1;
      hPrimTrkWeight->SetBinContent(bin, hPrimTrkWeight->GetBinContent(bin)/weight);
      hPrimTrkWeight->SetBinError(bin, hPrimTrkWeight->GetBinError(bin)/weight);
    }
  hPrimTrkWeight->SetLineColor(2);
  hPrimTrkWeight->Draw("sames HIST");
  printf("[i] Scale factor = %4.3f\n",hPrimTrkWeight->Integral()/hNPrimTrk->Integral());
  TLegend *leg = new TLegend(0.15, 0.4, 0.35, 0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hNPrimTrk,"Raw","L");
  leg->AddEntry(hPrimTrkWeight,"VPD eff. weighted","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_VpdEffCorr.pdf",run_type));

  if(0)
    {
      TFile *fout = TFile::Open("Rootfiles/Run15_pAu200.VpdEff.root","recreate");
      funcEff->Write();
      fout->Close();
    }
}

//================================================
void Nevents(const char *icent = "6080")
{
  if(year!=2014)
    {
      printf("[e] Not suitable for %d\n",year);
      return;
    }

  // Get the dimuon events number

  /*
  TFile *fdata = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");
  */

   
  TFile *fdata = TFile::Open(Form("./Rootfiles/Run14AuAu200MTD_RunNumbers_before_cuts.root"),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("hRunNumber");
  TH1F *hEvtRunAcc = 0x0;
  

  // =============================================
  // Effective number of MB events
  printf("+++++++++++++++++++++++++++++++++\n");
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents =  (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",icent));
  
  double nAnaEventsAll = 0;
  double nAnaEventNoEqMb = 0;

  double mb_events[5];
  for(int i=0; i<5; i++) mb_events[i] = 0;
  for(int bin=1; bin<=hEvtRun->GetNbinsX(); bin++)
    {
      if(hEvtRun->GetBinContent(bin)<=0) continue;
      if(hEvtRunAcc && hEvtRunAcc->GetBinContent(bin)<=0) continue;
      int run = TMath::FloorNint(hEvtRun->GetBinCenter(bin));
      double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(run));
      if(nEventsTaken==0) 
	{
	  if(i==0) printf("[w] check run %d\n",run);
	  continue;
	}
      double nEventsRun = hEvtRun->GetBinContent(bin);
      nAnaEventsAll += nEventsRun;
      double rf = hRF->GetBinContent(hRF->FindFixBin(run));
      if(rf==0)
	{
	  printf("[w] rf = 0 for run %d\n",run);
	  rf = 0.49;
	}
      double eq_mb = hEqMbEvents->GetBinContent(hEqMbEvents->FindFixBin(run));
      if(eq_mb==0)
	{
	  nAnaEventNoEqMb += nEventsRun;
	  printf("[w] eq_mb = 0 for run %d, nAnaEvent = %1.0f\n",run, nEventsRun);
	}
      mb_events[0] += nEventsRun/rf/nEventsTaken * eq_mb;
    }
  printf("[i] # of analyzed dimuon events: %4.4e\n",nAnaEventsAll);
  printf("[i] Missing equivalent MB events: %4.4e\n",nAnaEventNoEqMb);
  printf("[i] Effective # of MB events for %s%%: %4.4e\n",icent,mb_events[0]/(nAnaEventsAll-nAnaEventNoEqMb)*nAnaEventsAll);

  const char *trgSetupName[4] = {"","_low","_mid","_high"};
  for(int j=1; j<5; j++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type,trgSetupName[j-1]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hEvtRun->FindBin(runnumber+0.1);
	  if(hEvtRun->GetBinContent(bin)<=0) continue;
	  if(hEvtRunAcc && hEvtRunAcc->GetBinContent(bin)<=0) continue;
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
	      printf("[w] rf = 0 for run %1.0f (%s)\n",runnumber,trgSetupName[j-1]);
	      rf = 0.49;
	    }
	  double eq_mb = hEqMbEvents->GetBinContent(hEqMbEvents->FindFixBin(runnumber));
	  mb_events[j] += nEventsRun/rf/nEventsTaken * eq_mb;
	}
      printf("[i] # of events for %s: %4.4e\n",trgSetupName[j-1],mb_events[j]);
    }
  printf("[i]Check %1.0f =? %1.0f\n",mb_events[0],mb_events[1]+mb_events[2]+mb_events[3]+mb_events[4]);
  printf("+++++++++++++++++++++++++++++++++\n");

}

//================================================
void Lumi2013(const int savePlot = 0, const int saveHisto = 0)
{
  f = TFile::Open(Form("./output/Run13.MB.VtxCut.root"),"read");

  THnSparseF *hn = (THnSparseF*)f->Get("hEventVtx");
  hn->GetAxis(0)->SetRangeUser(start_run,end_run);

  // VPDMB efficiency
  hn->GetAxis(1)->SetRange(2,2);
  hn->GetAxis(3)->SetRange(3,3);
  TH1F *hBbcTpcVz = (TH1F*)hn->Projection(4);
  hBbcTpcVz->SetName("hBbcTpcVz");
  c = draw1D(hBbcTpcVz,"BBCMB: distribution of TPC vz with Ranking>0;vz_{TPC} (cm)",false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/BBCMB_TpcVz_PosRank.pdf",run_type));

  TH1F *hGoodBbc = (TH1F*)hn->Projection(0);
  hGoodBbc->SetName("hGoodBbc");
  hn->GetAxis(2)->SetRange(2,2);
  TH1F *hGoodBbcVpd = (TH1F*)hn->Projection(0);
  hGoodBbcVpd->SetName("hGoodBbcVpd");
  hGoodBbcVpd->Sumw2();
  hGoodBbc->Sumw2();
  hGoodBbcVpd->Divide(hGoodBbc);
  hGoodBbcVpd->SetMarkerStyle(21);
  hGoodBbcVpd->SetMarkerSize(1.0);
  TF1 *func1 = new TF1("func1","pol0",start_run,end_run);
  hGoodBbcVpd->Fit(func1,"IR0");
  c = draw1D(hGoodBbcVpd,"Efficiency of VPDMB trigger;RunId;(BBC&PosRank&VPD)/(BBC&PosRank)",false,true);
  func1->SetLineColor(2);
  func1->Draw("sames");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/BBCMB_VpdEff.pdf",run_type));
  
  hn->GetAxis(1)->SetRange(0,-1);
  hn->GetAxis(2)->SetRange(0,-1);
  hn->GetAxis(3)->SetRange(0,-1);

  // Obtain vertex cut efficiency for VPDMB
  hn->GetAxis(2)->SetRange(2,2);
  TH1F *hVpdAll = (TH1F*)hn->Projection(0);
  hVpdAll->SetName("hVpdAll");
  TH2F *hVpdVsZdcAll = (TH2F*)hn->Projection(0,6);
  hVpdVsZdcAll->SetName("hVpdVsZdcAll");

  TH1F *hTpcVpdDz[3];
  hTpcVpdDz[0]= (TH1F*)hn->Projection(5);
  hTpcVpdDz[0]->SetName("hTpcVpdDz_All");  

  hn->GetAxis(3)->SetRange(3,3);
  hTpcVpdDz[1]= (TH1F*)hn->Projection(5);
  hTpcVpdDz[1]->SetName("hTpcVpdDz_PosRank");

  hn->GetAxis(4)->SetRangeUser(-50,50);
  hTpcVpdDz[2]= (TH1F*)hn->Projection(5);
  hTpcVpdDz[2]->SetName("hTpcVpdDz_vz50cm");
  TH2F *hDzVsRunRaw = (TH2F*)hn->Projection(5,0);
  hDzVsRunRaw->SetName("hDzVsRunRaw");

  for(int i=0; i<3; i++)
    {
      hTpcVpdDz[i]->SetLineColor(color[i]);
    }
  c = draw1D(hTpcVpdDz[0],"VPDMB: distribution of vz difference between TPC and VPD;vz_{TPC} - vz_{VPD} (cm)",true,false);
  hTpcVpdDz[1]->Draw("samesHIST");
  hTpcVpdDz[2]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.15, 0.65, 0.3, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hTpcVpdDz[0],"All","L");
  leg->AddEntry(hTpcVpdDz[1],"+Ranking > 0","L");
  leg->AddEntry(hTpcVpdDz[2],"+|vz_{TPC}| < 50 cm","L");
  leg->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_TpcVpdDz.pdf",run_type));

  hn->GetAxis(5)->SetRangeUser(-5,5);
  TH1F *hVpdAcc = (TH1F*)hn->Projection(0);
  hVpdAcc->SetName("hVpdAcc");
  TH2F *hVpdVsZdcAcc = (TH2F*)hn->Projection(0,6);
  hVpdVsZdcAcc->SetName("hVpdVsZdcAcc");

  TH1F *hVpdFrac = (TH1F*)hVpdAcc->Clone("hVpdFrac");
  hVpdAll->Sumw2();
  hVpdFrac->Divide(hVpdAll);
  hVpdFrac->SetMarkerStyle(21);
  c = draw1D(hVpdFrac,"Fraction of VPDMB events within vertxing cuts;RunId;fraction");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VPDMB_VtxCutEff.pdf",run_type));
  

  // Get equivalent # of MB events
  TFile *fin = 0;
  if(!saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  else           fin = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
  TH1F *hNeventMtd = (TH1F*)fin->Get("hNevents_MTD-dimuon");
  TH1F *hPrescaleMtd = (TH1F*)fin->Get("hPreScale_MTD-dimuon");
  TH1F *hLivetimeMtd = (TH1F*)fin->Get("hLiveTime_MTD-dimuon");
  TH1F *hNeventVpd = (TH1F*)fin->Get("hNevents_VPDMB");
  TH1F *hPrescaleVpd = (TH1F*)fin->Get("hPreScale_VPDMB");
  TH1F *hLivetimeVpd = (TH1F*)fin->Get("hLiveTime_VPDMB");

  int nbins = hNeventMtd->GetNbinsX();
  double nMtdEvents = 0;
  double nVpdEvents = 0;
  double nVpdEventsAcc = 0;

  TH1F *hDzWeight = 0x0;
  //hDzVsRunRaw->SetName("hDzVsRunRaw");
  c = draw2D(hDzVsRunRaw);
  for(int bin=1; bin<=nbins; bin++)
    {
      double mtdEvts =  hNeventMtd->GetBinContent(bin);
      if(mtdEvts<=0) continue;
      if(hPrescaleMtd->GetBinContent(bin)==0) continue;
      if(hLivetimeVpd->GetBinContent(bin)==0) continue;
      nMtdEvents += mtdEvts;
      
      int jbin = hVpdAcc->FindFixBin(hNeventMtd->GetBinCenter(bin));
      double scale = hPrescaleVpd->GetBinContent(bin)/hPrescaleMtd->GetBinContent(bin) * hLivetimeMtd->GetBinContent(bin)/hLivetimeVpd->GetBinContent(bin);
      nVpdEventsAcc += hNeventVpd->GetBinContent(bin) * hVpdFrac->GetBinContent(jbin) * scale;
      //printf("[i] Run %d for VPD/MTD = %4.2f\n",14130000+bin-1,hNeventVpd->GetBinContent(bin)*scale/mtdEvts);


      int kbin = hDzVsRunRaw->GetXaxis()->FindFixBin(hNeventMtd->GetBinCenter(bin));
      cout << bin << "  " << kbin << endl;
      
      TH1F *htmp = (TH1F*)hDzVsRunRaw->ProjectionY(Form("hDz_%d",bin),kbin,kbin);
      
      if(!hDzWeight)
	{
	  hDzWeight = (TH1F*)htmp->Clone("hDzWeight");
	  hDzWeight->Scale(scale);
	}
      else
	{
	  hDzWeight->Add(htmp,scale);
	}
      
      //cout << nVpdEventsAcc << endl;
      //cout << bin << "  " << hNeventMtd->GetBinCenter(bin) << "  " << jbin << "  " << hVpdFrac->GetBinContent(jbin)
      //	   << "  " << hPrescaleVpd->GetBinContent(bin) << "  " << hLivetimeVpd->GetBinContent(bin) 
      //	   << "  " << hPrescaleMtd->GetBinContent(bin) << "  " << hLivetimeMtd->GetBinContent(bin) << endl;
    }
  printf("[i] Processed MTD = %1.0f --> VPD = %1.0f\n",nMtdEvents,nVpdEventsAcc);

  c = draw1D(hDzWeight);

  TH1F *hMBevents = new TH1F("hMBevents","",2,0,2);
  hMBevents->SetBinContent(1,nMtdEvents);
  hMBevents->SetBinContent(2,nVpdEventsAcc);
  if(saveHisto)
    {
      hMBevents->Write("",TObject::kOverwrite);
    }
}

//================================================
void Lumi2014(const int savePlot = 0, const int saveHisto = 0)
{
  TList *list = new TList;
  TString legName[3] = {"|vr_{TPC}| < 2 cm", "+|vz_{TPC}| < 100 cm", "+ |vz_{TPC}-vz_{VPD}| < 3 cm"};
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  const char *setupName[4] = {"prod","prod_low","prod_mid","prod_high"};
  const int nCentBins = 16;
  const char *cent_Name[nCentBins] = {"0-80","0-20","20-40","40-60","60-80","0-60","0-10","10-30","30-60","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
  const char *cent_Title[nCentBins] = {"0080","0020","2040","4060","6080","0060","0010","1030","3060","1020","2030","3040","4050","5060","6070","7080"};

  TString legName2[4];
  for(int i=0; i<4; i++)
    {
      legName2[i] = Form("AuAu_200_%s_2014",trgSetupName[i]);
    }
  TString legName3[2] = {"w/o weights","w/ weights"};
  TString legName4[5];
  for(int i=0; i<5; i++) legName4[i] = Form("%s%%",cent_Name[i]);
  const int max_vz = 100;
  const int max_dz = 3;
  const int max_vr = 2;

  // MTD triggers
  TFile *fMtd = TFile::Open(Form("./output/Run14_AuAu200.jpsi.root"),"read");
  TH2F *hCentWeight = (TH2F*)fMtd->Get("mhCentWeight_di_mu");
  c = draw2D(hCentWeight,"Event weights vs. multiplicity");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/EventWeight.pdf",run_type));
  if(gSaveAN)   c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_EventWeight.pdf"));

  /// vertex finding efficiency
  TH1F *hEvtRun = (TH1F*)fMtd->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fMtd->Get("mhEvtRunAcc_di_mu");
  TH1F *hMtdVtxEff = (TH1F*)hEvtRunAcc->Clone("hMtdVtxEff");
  hMtdVtxEff->Sumw2();
  TGraphAsymmErrors *graph = new TGraphAsymmErrors(hEvtRunAcc,hEvtRun,"cl=0.683 b(1,1) mode");
  graph->SetName("gMtdVtxFindEff");
  hMtdVtxEff->Divide(hEvtRun);
  double x,y;
  for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
    {
      graph->GetPoint(ipoint,x,y);
      int bin = hMtdVtxEff->FindFixBin(x);
      double err1 = graph->GetErrorYhigh(ipoint);
      double err2 = graph->GetErrorYlow(ipoint);
      double err = (err1>err2) ? err1 : err2;
      hMtdVtxEff->SetBinError(bin,err);
      if(y==0)
	{
	  // mainly due to tof not included in the run
	  // consequently, vpd information is not read out
	  printf("[w] 0 efficiency for run %d: %1.0f/%1.0f\n",bin+start_run-1,hEvtRunAcc->GetBinContent(bin),hEvtRun->GetBinContent(bin));
	}
    }

  TH1F *hMtdVtxEffLumi[4];
  for(int i=0; i<4; i++)
    {
      hMtdVtxEffLumi[i] = (TH1F*)hMtdVtxEff->Clone(Form("MtdEffPerRun_%s_cent0100",trgSetupName[i]));
      hMtdVtxEffLumi[i]->Reset();

      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[i]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hMtdVtxEff->FindFixBin(runnumber);
	  hMtdVtxEffLumi[i]->SetBinContent(bin,hMtdVtxEff->GetBinContent(bin));
	  hMtdVtxEffLumi[i]->SetBinError(bin,hMtdVtxEff->GetBinError(bin));
	}
      list->Add(hMtdVtxEffLumi[i]);
    }
  c = drawHistos(list,"MtdVtxEffLumi",Form("Dimuon: vertex cut efficiency (centrality integrated)"),kFALSE,0,30,kTRUE,0.5,1.1,kFALSE,kTRUE,legName2,kTRUE,"trgsetupname",0.3,0.5,0.16,0.36,true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_VtxCutEffLumi.pdf",run_type));
  list->Clear();

  // VPD-ZDC-novtx-mon
  TFile *fMB = TFile::Open(Form("./Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
 
  // centrality
  TH1F *hMbCent[2];
  hMbCent[0] = (TH1F*)fMB->Get("Cent_all");
  hMbCent[1] = (TH1F*)fMB->Get("Cent_all_w");
  for(int i=0; i<2; i++)
    {
      hMbCent[i]->SetLineWidth(2);
      list->Add(hMbCent[i]);
      if(i==1)
	{
	  for(int bin=1; bin<=hMbCent[i]->GetNbinsX()/2; bin++)
	    {
	      cout << "bin = " << bin << ": " << (hMbCent[i]->GetBinContent(bin*2-1)+hMbCent[i]->GetBinContent(bin*2)) << endl;
	    }
	}
    }
  c = drawHistos(list,"MbCent","VPD-ZDC-novtx-mon: distribution of centrality bins;cent",true,0,16,true,0,1.2*hMbCent[1]->GetMaximum(),kFALSE,kTRUE,legName3,kTRUE,"",0.5,0.7,0.2,0.35,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_CentWeight.pdf",run_type));
  if(gSaveAN)   c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_Centrality.pdf"));
  list->Clear();

  TH1F *hMbCentLumi[4];
  for(int i=0; i<4; i++)
    {
      hMbCentLumi[i] = (TH1F*)fMB->Get(Form("Cent_%s_w",setupName[i]));
      //hMbCentLumi[i]->Rebin(2);
      //hMbCentLumi[i]->Scale(1./hMbCentLumi[i]->Integral(15,16)*2);
      hMbCentLumi[i]->Scale(1./hMbCentLumi[i]->Integral()*16);
      hMbCentLumi[i]->SetLineWidth(2);
      list->Add(hMbCentLumi[i]);
    }
  c = drawHistos(list,"MbCentLumi","VPD-ZDC-novtx-mon: distribution of centrality bins;cent",true,0,16,kTRUE,0.5,1.2,kFALSE,kTRUE,legName2,kTRUE,"trgsetupname",0.3,0.5,0.16,0.36,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_CentVsLumi.pdf",run_type));
  list->Clear();

  // Vertex distribution
  TH1F *hTpcVz[4];
  for(int i=0; i<4; i++)
    {
      hTpcVz[i] = (TH1F*)fMB->Get(Form("TpcVz_cent0080_%s",setupName[i]));
      hTpcVz[i]->Scale(1./hTpcVz[i]->Integral());
      list->Add(hTpcVz[i]);
    }
  c = drawHistos(list,"TpcVz","VPD-ZDC-novtx-mon: TPC vz distribution (0-80%);vz_{TPC} (cm)",false,0,30,true,0,0.02,kFALSE,kTRUE,legName2,kTRUE,"0-80%",0.18,0.38,0.68,0.88,false);
  TLine *line = GetLine(-1*max_vz,0,-1*max_vz,0.8*hTpcVz[3]->GetMaximum());
  line->Draw();
  line = GetLine(max_vz,0,max_vz,0.8*hTpcVz[3]->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVz.pdf",run_type));
  list->Clear();

  TH1F *hTpcVr[4];
  for(int i=0; i<4; i++)
    {
      hTpcVr[i] = (TH1F*)fMB->Get(Form("TpcVr_cent0080_%s",setupName[i]));
      hTpcVr[i]->Scale(1./hTpcVr[i]->Integral());
      list->Add(hTpcVr[i]);
    }
  c = drawHistos(list,"TpcVr","VPD-ZDC-novtx-mon: TPC vr distribution (0-80%);vr_{TPC} (cm)",false,0,30,true,0.00001,1000,true,kTRUE,legName2,kTRUE,"0-80%",0.18,0.38,0.68,0.88,false);
  line = GetLine(max_vr,0,max_vr,0.8*hTpcVr[3]->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVr.pdf",run_type));
  list->Clear();

  TH1F *hTpcVpdDz[4];
  for(int i=0; i<4; i++)
    {
      hTpcVpdDz[i] = (TH1F*)fMB->Get(Form("TpcVpdDz_cent0080_%s",setupName[i]));
      hTpcVpdDz[i]->Scale(1./hTpcVpdDz[i]->Integral());
      list->Add(hTpcVpdDz[i]);
    }
  c = drawHistos(list,"TpcVpdDz","VPD-ZDC-novtx-mon:#Deltaz distribution for vertices (0-80%);vz_{TPC}-vz_{VPD} (cm)",true,-10,10,true,1e-6,1e3,true,kTRUE,legName2,kTRUE,"0-80%",0.18,0.38,0.68,0.88,false);
  line = GetLine(-1*max_dz,0,-1*max_dz,0.8*hTpcVpdDz[3]->GetMaximum());
  line->Draw();
  line = GetLine(max_dz,0,max_dz,0.8*hTpcVpdDz[3]->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVpdDz.pdf",run_type));
  list->Clear();

  // vertex cut efficiency
  TH1F *hNEvents[nCentBins][4];
  for(int i=0; i<nCentBins; i++)
    {
      hNEvents[i][0] = (TH1F*)fMB->Get(Form("NEvents_cent%s_all",cent_Title[i]));
      hNEvents[i][1] = (TH1F*)fMB->Get(Form("NEvents_VrCut_cent%s_all",cent_Title[i]));
      hNEvents[i][2] = (TH1F*)fMB->Get(Form("NEvents_VrVzCut_cent%s_all",cent_Title[i]));
      hNEvents[i][3] = (TH1F*)fMB->Get(Form("NEvents_VrVzDzCut_cent%s_all",cent_Title[i]));
    }

  TH1F *hVtxEff[nCentBins][3];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hVtxEff[i][j] = (TH1F*)hNEvents[i][j+1]->Clone(Form("%s_eff",hNEvents[i][j+1]->GetName()));
	  TGraphAsymmErrors *graph = new TGraphAsymmErrors(hVtxEff[i][j],hNEvents[i][0],"cl=0.683 b(1,1) mode");
	  graph->SetName(Form("graph_%d_%d",i,j));
	  hVtxEff[i][j]->Divide(hNEvents[i][0]);
	  double x,y;
	  for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
	    {
	      graph->GetPoint(ipoint,x,y);
	      int bin = hVtxEff[i][j]->FindFixBin(x);
	      double err1 = graph->GetErrorYhigh(ipoint);
	      double err2 = graph->GetErrorYlow(ipoint);
	      double err = (err1>err2) ? err1 : err2;
	      hVtxEff[i][j]->SetBinError(bin,err);
	      if(i==0 && y==0)
		{
		  printf("[w] MB: 0 efficiency for run %d: %1.0f/%1.0f\n",bin+start_run-1,hNEvents[i][j+1]->GetBinContent(bin),hNEvents[i][0]->GetBinContent(bin));
		}
	      if(y!=hVtxEff[i][j]->GetBinContent(bin))
		{
		  //printf("[e] Mismatch eff: %4.2f != %4.2f\n",y,hVtxEff[i][j]->GetBinContent(bin));
		}
	    }
	}
    }

  TCanvas *c = new TCanvas("MB_VtxCutEff","MB_VtxCutEff",1100,650);
  c->Divide(3,2);
  for(int i=0; i<5; i++)
    {
      c->cd(i+1);
      for(int j=0; j<3; j++)
	{
	  hVtxEff[i][j]->SetMarkerStyle(20);
	  hVtxEff[i][j]->SetMarkerColor(color[j]);
	  hVtxEff[i][j]->SetLineColor(color[j]);
	  hVtxEff[i][j]->GetYaxis()->SetRangeUser(0.5,1.1);
	  hVtxEff[i][j]->SetTitle("");
	  if(j==0) hVtxEff[i][j]->Draw();
	  else     hVtxEff[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("VPD-ZDC-novtx-mon: vertex cut efficiency (%s%%)",cent_Name[i]),0.05);
	  t1->Draw();
	}
      if(i==0)
	{
	  TLegend *leg = new TLegend(0.5,0.15,0.7,0.4);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  for(int j=0; j<3; j++)
	    {
	      leg->AddEntry(hVtxEff[i][j],legName[j].Data(),"P");
	    }
	  leg->Draw();
	}
    }
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_VtxCutEff.pdf",run_type));
  if(gSaveAN)   c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_VtxCutEff.png"));

  // different luminosity
  TH1F *hVtxEffLumi[nCentBins][4];
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hVtxEffLumi[j][i] = (TH1F*)hVtxEff[j][2]->Clone(Form("EffPerRun_TpcVpdDzCut_%s_cent%s",trgSetupName[i],cent_Title[j]));
	  hVtxEffLumi[j][i]->Reset();
	}
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[i]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  for(int j=0; j<nCentBins; j++)
	    {
	      int bin = hVtxEff[j][1]->FindFixBin(runnumber);
	      hVtxEffLumi[j][i]->SetBinContent(bin,hVtxEff[j][1]->GetBinContent(bin));
	      hVtxEffLumi[j][i]->SetBinError(bin,hVtxEff[j][1]->GetBinError(bin));
	    }
	}
    }
  c = new TCanvas("MB_VtxCutEffLumi","MB_VtxCutEffLumi",1100,650);
  c->Divide(2,2);
  for(int j=0; j<4; j++)
    {
      c->cd(j+1);
      for(int i=0; i<4; i++)
	{
	  hVtxEffLumi[j][i]->SetMarkerStyle(20);
	  hVtxEffLumi[j][i]->SetMarkerColor(color[i]);
	  hVtxEffLumi[j][i]->SetLineColor(color[i]);
	  hVtxEffLumi[j][i]->GetYaxis()->SetRangeUser(0.5,1.1);
	  hVtxEffLumi[j][i]->SetTitle("");
	  if(i==0) hVtxEffLumi[j][i]->Draw();
	  else     hVtxEffLumi[j][i]->Draw("sames");
	  if(i==0)
	    {
	      TPaveText *t1 = GetTitleText(Form("VPD-ZDC-novtx-mon: vertex cut efficiency (%s%%)",cent_Name[j]),0.05);
	      t1->Draw();
	    }
	}
      if(j==0)
	{
	  TLegend *leg = new TLegend(0.3,0.15,0.5,0.4);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  for(int i=0; i<4; i++)
	    {
	      leg->AddEntry(hVtxEffLumi[j][i],legName2[i].Data(),"P");
	    }
	  leg->Draw();
	}
    }
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_VtxCutEffLumi.pdf",run_type));

  // rejection factor during production
  const int nRuns = 38;
  const int runIDs[nRuns] = {15078103, 15078104, 15078107, 15078108, 15079059, 
			     15079061, 15084002, 15084022, 15084052, 15088003, 
			     15090006, 15097032, 15097034, 15102021, 15104018, 
			     15104039, 15104059, 15106008, 15106009, 15106010, 
			     15106011, 15107077, 15110032, 15110038, 15119021, 
			     15132026, 15142054, 15151035, 15151036, 15151037, 
			     15151038, 15151039, 15151040, 15151041, 15151042, 
			     15151043, 15162019, 15166023};
  TFile *frf = TFile::Open("./output/Run14_AuAu200.RejectFactor.root","read");
  TH1F *hEvtAll = (TH1F*)frf->Get("hEvtAll");
  TH1F *hEvtAcc = (TH1F*)frf->Get("hEvtAcc");
  TH1F *hRF = (TH1F*)hEvtAcc->Clone("hRejectFactor_dimuon");
  hRF->Sumw2();
  hRF->Divide(hEvtAll);

  list->Clear();
  TH1F *hRFLumi[4];
  TF1 *funcRF[4];
  TCanvas *c = new TCanvas("Fit_RF_Lumi","Fit_RF_Lumi",800,600);
  c->Divide(2,2);
  for(int i=0; i<4; i++)
    {
      hRFLumi[i] = (TH1F*)hRF->Clone(Form("RejectFactor_%s_cent01000",trgSetupName[i]));
      hRFLumi[i]->Reset();

      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[i]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hRF->FindFixBin(runnumber);
	  bool isGoodRun = true;
	  for(int irun=0; irun<nRuns; irun++)
	    {
	      if(runnumber==runIDs[irun])
		{
		  isGoodRun = false;
		  break;
		}
	    }
	  if(!isGoodRun) continue;
	  double rf = hRF->GetBinContent(bin);
	  hRFLumi[i]->SetBinContent(bin,rf);
	  hRFLumi[i]->SetBinError(bin,hRF->GetBinError(bin));
	  if((rf>0.51 || rf<0.466) && rf > 0)
	    {
	      cout << "Bad run for RF: " << runnumber << " with rf = " << rf << endl;
	    }
	}

      TH1F *hRFtmp = (TH1F*)hRFLumi[i]->Clone(Form("Fit_%s",hRFLumi[i]->GetName()));
      funcRF[i] = new TF1(Form("func_%s",hRFLumi[i]->GetName()),"pol0",start_run,end_run);
      hRFtmp->Fit(funcRF[i],"IR0Q");
      c->cd(i+1);
      hRFtmp->SetMarkerStyle(25);
      hRFtmp->GetYaxis()->SetRangeUser(0.4,0.6);
      hRFtmp->SetTitle("");
      hRFtmp->Draw();
      funcRF[i]->SetLineColor(2);
      funcRF[i]->Draw("sames");
      TPaveText *t1 = GetTitleText("Dimuon: fraction of accepted events",0.05);
      t1->Draw();
      t1 = GetPaveText(0.2,0.6,0.15,0.2,0.05);
      t1->AddText(Form("AuAu_200_%s_2014",trgSetupName[i]));
      t1->Draw();
      list->Add(hRFLumi[i]);
    }
  c = drawHistos(list,"MtdRFLumi",Form("Dimuon: fraction of events accepted during two-pass reconstruction"),kFALSE,0,30,kTRUE,0.35,0.55,kFALSE,kTRUE,legName2,kTRUE,"trgsetupname",0.3,0.5,0.16,0.36,true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_RejectFactorLumi.pdf",run_type));
  if(gSaveAN)   c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_RejectFactor.pdf"));
  list->Clear();

  // check missing runs
  printf("+++ check RF +++\n");
  int nRF = 0;
  double nRFevts = 0;
  int nbins = hEvtRunAcc->GetNbinsX();
  for(int ibin=1; ibin<=nbins; ibin++)
    {
      double run = hEvtRunAcc->GetBinCenter(ibin);
      if(hEvtRunAcc->GetBinContent(ibin)>0)
	{
	  int bin = hRF->FindFixBin(run);
	  if(hRF->GetBinContent(bin)<=0)
	    {
	      nRF++;
	      int index = -1;
	      for(int i=0; i<4; i++)
		{
		  if(index>-1) break;
		  for(int jbin=1; jbin<=hMtdVtxEffLumi[i]->GetNbinsX(); jbin++)
		    {
		      double eff = hMtdVtxEffLumi[i]->GetBinContent(hMtdVtxEffLumi[i]->FindFixBin(run));
		      if(eff>0)
			{
			  index = i;
			  break;
			}
		    }
		}
	      if(index<0) 
		{
		  printf("[w] Can't find run %1.0f in any trigger setup\n",run);
		  continue;
		}
	      hRF->SetBinContent(bin,funcRF[index]->GetParameter(0));
	      hRF->SetBinError(bin,funcRF[index]->GetParError(0));
	      nRFevts += hEvtRunAcc->GetBinContent(ibin);
	      printf("[w] Miss RF for run %1.0f with %4.1fK (%s)\n",run,hEvtRunAcc->GetBinContent(ibin)/1e3,trgSetupName[index]);
	    }
	}
    }
  printf("+++ %d runs (%4.2fM) are missing +++\n\n",nRF,nRFevts/1e6);
  
  // get equivalent # of MB events
  /// Raw # of events
  TH1F *hNEventsRawAll[5];
  TH1F *hNEventsRaw[nCentBins][5];
  TH1F *hNevtVtxCutWeight[nCentBins][5];

  hNEventsRawAll[0] = (TH1F*)fMB->Get("NEvents_all");
  for(int j=0; j<4; j++)
    {
      hNEventsRawAll[j+1] = (TH1F*)fMB->Get(Form("NEvents_%s",setupName[j]));
    }
  for(int i=0; i<nCentBins; i++)
    {
     hNEventsRaw[i][0] = (TH1F*)fMB->Get(Form("NEvents_cent%s_all",cent_Title[i]));
     hNevtVtxCutWeight[i][0] = (TH1F*)fMB->Get(Form("NEvents_VrVzDzCut_cent%s_all_w",cent_Title[i]));
     for(int j=0; j<4; j++)
       {
	 hNevtVtxCutWeight[i][j+1] = (TH1F*)fMB->Get(Form("NEvents_VrVzDzCut_cent%s_%s_w",cent_Title[i],setupName[j]));
       }
    }

  TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
  TH1F *hEqMbEvt = (TH1F*)fout->Get("hEqMbEvt_dimuon");
  TH1F *hEqMbEvtVtxCutWeight[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hEqMbEvtVtxCutWeight[i] = (TH1F*)hEqMbEvt->Clone(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[i]));
      hEqMbEvtVtxCutWeight[i]->Reset();
      for(int ibin=1; ibin<=hEqMbEvt->GetNbinsX(); ibin++)
	{
	  double run = hEqMbEvt->GetBinCenter(ibin);
	  double eqevt = hEqMbEvt->GetBinContent(ibin);
	  if(eqevt<=0) continue;
	  int bin = hNEventsRawAll[0]->FindFixBin(run);
	  double rawevt = hNEventsRawAll[0]->GetBinContent(bin);
	  double evtw = hNevtVtxCutWeight[i][0]->GetBinContent(bin);
	  if(rawevt<=0) continue;
	  double value = eqevt*evtw/rawevt;
	  double error = evtw<=0? 0 : value * 1./sqrt(evtw);
	  hEqMbEvtVtxCutWeight[i]->SetBinContent(ibin, value);
	  hEqMbEvtVtxCutWeight[i]->SetBinError(ibin, error);
	}
    }

  // Fit event fraction after vertex and centrality selections
  TH1F *hNevtRatio[nCentBins][4];
  TF1 *funcMbVtxEff[nCentBins][4];
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hNevtRatio[j][i] = (TH1F*)hNevtVtxCutWeight[j][i+1]->Clone(Form("Fit_%s",hNevtVtxCutWeight[j][i+1]->GetName()));
	  hNevtRatio[j][i]->Divide(hNEventsRawAll[i+1]);
	  if(i==3)
	    {
	      double run = 15119021;
	      int bin = hNevtRatio[j][i]->FindFixBin(run);
	      hNevtRatio[j][i]->SetBinContent(bin,0);
	      hNevtRatio[j][i]->SetBinError(bin,0);
	    }
	  funcMbVtxEff[j][i] = new TF1(Form("Func_%s",hNevtVtxCutWeight[j][i+1]->GetName()),"pol0",start_run,end_run);
	  hNevtRatio[j][i]->Fit(funcMbVtxEff[j][i],"IR0Q");
	  printf("[i] %s%% for %s: %4.2f%%\n",cent_Name[j],setupName[i],funcMbVtxEff[j][i]->GetParameter(0)*100);
	}
    }

  c = new TCanvas("Fit_MB_VtxEff","Fit_MB_VtxEff",1200,700);
  c->Divide(5,4);
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<5; j++)
	{
	  c->cd(i*5+j+1);
	  hNevtRatio[j][i]->SetMarkerStyle(29);
	  hNevtRatio[j][i]->GetYaxis()->SetRangeUser(0,1);
	  hNevtRatio[j][i]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%s: %s%%",setupName[i],cent_Name[j]),0.08);
	  t1->Draw();
	  funcMbVtxEff[j][i]->SetLineColor(2);
	  funcMbVtxEff[j][i]->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.5,0.7,0.4,0.6,0.08);
	  t1->AddText(Form("p0 = %4.2f",funcMbVtxEff[j][i]->GetParameter(0)));
	  t1->Draw();
	}
    }
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/MB_EvtFractionWithVtxCutCent.pdf",run_type));
  if(gSaveAN)   c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch1_FitVtxEff.pdf"));

  // check vertex efficiency
  printf("+++ check vertex efficiency +++\n");
  int nVE = 0;
  double nVEevts = 0;
  for(int ibin=1; ibin<=hEvtRunAcc->GetNbinsX(); ibin++)
    {
      double run = hEvtRunAcc->GetBinCenter(ibin);
      if(hEvtRunAcc->GetBinContent(ibin)>0)
	{
	  for(int i=0; i<nCentBins; i++)
	    {
	      int bin = hEqMbEvtVtxCutWeight[i]->FindFixBin(run);
	      if(hEqMbEvtVtxCutWeight[i]->GetBinContent(bin)<=0)
		{
		  if(i==0) nVE++;
		  int index = -1;
		  for(int j=0; j<4; j++)
		    {
		      if(index>-1) break;
		      for(int jbin=1; jbin<=hMtdVtxEffLumi[j]->GetNbinsX(); jbin++)
			{
			  double eff = hMtdVtxEffLumi[j]->GetBinContent(hMtdVtxEffLumi[j]->FindFixBin(run));
			  if(eff>0)
			    {
			      index = j;
			      break;
			    }
			}
		    }
		  if(index<0 && i==0) 
		    {
		      printf("[w] Can't find run %1.0f in any trigger setup\n",run);
		      continue;
		    }
		  if(i==0) 
		    {
		      nVEevts += hEvtRunAcc->GetBinContent(ibin);
		      printf("[w] Miss vtx eff for run %1.0f with %4.1fK (%s)\n",run,hEvtRunAcc->GetBinContent(ibin)/1e3,trgSetupName[index]);
		    }
		  double eqevt = hEqMbEvt->GetBinContent(hEqMbEvt->FindFixBin(run));
		  hEqMbEvtVtxCutWeight[i]->SetBinContent(bin,funcMbVtxEff[i][index]->GetParameter(0)*eqevt);
		  hEqMbEvtVtxCutWeight[i]->SetBinError(bin,funcMbVtxEff[i][index]->GetParError(0)*eqevt);
		}
	    }
	}
    }
  printf("+++ %d runs (%4.2fM) are missing +++\n\n",nVE,nVEevts/1e6);

  list->Clear();
  for(int i=0; i<5; i++)
    {
      list->Add(hEqMbEvtVtxCutWeight[i]);
    }
  c = drawHistos(list,"MtdEqMbEvt",Form("Dimuon: equivalent MB events after vertex cuts with event weights"),kFALSE,0,30,false,0.35,0.55,kFALSE,kTRUE,legName4,kTRUE,"Centrality",0.3,0.5,0.68,0.86,true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_EquMbEvtVtxCutWeight.pdf",run_type));
  list->Clear();


  if(saveHisto)
    {
      fout->cd();
      for(int i=0; i<nCentBins; i++)
	{
	  hEqMbEvtVtxCutWeight[i]->Write("",TObject::kOverwrite);
	}
      hRF->Write("",TObject::kOverwrite);
    }
}


//================================================
void Lumi2015(const int savePlot = 0, const int saveHisto = 0)
{
  const char *run_type = "Run15_pAu200";
  const char *trigName = "VPDMBnovtx";
  f = TFile::Open(Form("./output/%s.MB.root",run_type),"read");
  cout << f->GetName() << endl;

  // luminosity distributions
  TProfile *hBbcRate = (TProfile*)f->Get("hBbcRate");
  hBbcRate->SetMarkerStyle(24);
  c = draw1D(hBbcRate,Form("%s: BBC rate vs. runId;runId",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_BBCrate.pdf",run_type,trigName));

  TProfile *hZdcRate = (TProfile*)f->Get("hZdcRate");
  hZdcRate->SetMarkerStyle(24);
  c = draw1D(hZdcRate,Form("%s: ZDC rate vs. runId;runId",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_ZDCrate.pdf",run_type,trigName));

  TH2F *hBbcVsZdc = new TH2F("hBbcVsZdc",Form("%s: ZDC rate vs. BBC rate;BBC (kHz);ZDC (kHz)",run_type),300,0,3e3,450,0,45);
  for(int bin=1; bin<=hBbcRate->GetNbinsX(); bin++)
    {
      if(hBbcRate->GetBinContent(bin)<=0) continue;
      hBbcVsZdc->Fill(hBbcRate->GetBinContent(bin)/1e3,hZdcRate->GetBinContent(bin)/1e3);
    }
  c = draw2D(hBbcVsZdc);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_ZDCvsBBC.pdf",run_type,trigName));

  // vertex distributions
  TH1F* hTpcVz = (TH1F*)f->Get("hTpcVz");
  c = draw1D(hTpcVz,Form("%s: TPC vertex z distribution;vz_{TPC} (cm)",run_type),false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVz.pdf",run_type,trigName));

  TH1F* hVzDiff = (TH1F*)f->Get("hVzDiff");
  c = draw1D(hVzDiff,Form("%s: TPC-VPD vertex z distribution;vz_{TPC}-vz_{VPD} (cm)",run_type),false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDiff.pdf",run_type,trigName));

  const int nTrigId = 6;
  const int trigId[nTrigId] = {470004, 470014, 480004, 490004, 500004, 510004};
  TH1F *hTpcVzTrigId[nTrigId];
  TH2F *hTpcVzVsTrigId = (TH2F*)f->Get("hTpcVzVsTrigId");
  TList *list = new TList;
  TString legName[nTrigId];
  double max = 0;
  for(int i=0; i<nTrigId; i++)
    {
      hTpcVzTrigId[i] = (TH1F*)hTpcVzVsTrigId->ProjectionY(Form("hTpcVz_%d",trigId[i]),7+i, 7+i);
      if(hTpcVzTrigId[i]->GetEntries()<=0) continue;
      list->Add(hTpcVzTrigId[i]);
      legName[list->GetEntries()-1] = Form("Trigger Id: %d",trigId[i]);
      if(hTpcVzTrigId[i]->GetMaximum()>max) max = hTpcVzTrigId[i]->GetMaximum();
    }
  c = drawHistos(list,Form("%s_TpcVzTrigId",run_type),Form("%s: TPC vertex z distribution",run_type),false,0,10,true,0,1.3*max,kFALSE,kTRUE,legName,kTRUE,run_type,0.15,0.3,0.6,0.85,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVzTrigId.pdf",run_type,trigName));

  const int nRanking = 2;
  TH1F *hTpcVzRank[nRanking];
  TH2F *hTpcVzVsRanking = (TH2F*)f->Get("hTpcVzVsRanking");
  const char* rank_name[nRanking] = {"Ranking > 0", "Ranking < 0"};
  TString legName1[nRanking];
  double nEvents[nRanking];
  for(int i=0; i<nRanking; i++)
    {
      hTpcVzRank[i] = (TH1F*)hTpcVzVsRanking->ProjectionY(Form("hTpcVz_%d",i),3-2*i,3-2*i);
      list->Add(hTpcVzRank[i]);
      nEvents[i] = hTpcVzRank[i]->Integral(hTpcVzRank[i]->FindFixBin(-100+0.01), hTpcVzRank[i]->FindFixBin(100-0.01));
    }
  for(int i=0; i<nRanking; i++)
    {
      legName1[i] = Form("%s: %4.1f%%",rank_name[i],nEvents[i]/(nEvents[0]+nEvents[1])*100);
    }
  c = drawHistos(list,Form("%s_TpcVzRanking",run_type),Form("%s: TPC vertex z distribution;vz_{TPC} (cm)",run_type),false,0,10,false,0,1,kFALSE,kTRUE,legName1,kTRUE,"|vz_{TPC}| < 100 cm",0.13,0.3,0.65,0.85,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVzRank.pdf",run_type,trigName));

  TH1F *hVzDiffRank[nRanking];
  TH2F *hDiffVzVsRanking = (TH2F*)f->Get("hDiffVzVsRanking");
  for(int i=0; i<nRanking; i++)
    {
      hVzDiffRank[i] = (TH1F*)hDiffVzVsRanking->ProjectionY(Form("hVzDiff_%d",i),3-2*i,3-2*i);
      list->Add(hVzDiffRank[i]);
      nEvents[i] = hVzDiffRank[i]->Integral(hVzDiffRank[i]->FindFixBin(-6+1e-4), hVzDiffRank[i]->FindFixBin(6-1e-4));
    }
  for(int i=0; i<nRanking; i++)
    {
      legName1[i] = Form("%s: %4.1f%%",rank_name[i],nEvents[i]/(nEvents[0]+nEvents[1])*100);
    }
  c = drawHistos(list,Form("%s_VzDiffRanking",run_type),Form("%s: TPC-VPD vertex z distribution;vz_{TPC}-vz_{VPD} (cm)",run_type),false,0,10,false,0,1,true,kTRUE,legName1,kTRUE,"|vz_{TPC}-vz_{VPD}| < 6 cm",0.13,0.3,0.65,0.85,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDiffRank.pdf",run_type,trigName));

  // Tpc vs Vpd vz
  TH2F *hTpcVsVpdVz[4];
  const char *hname[4] = {"All","Rank","Dz","DzRank"};
  const char *htitle[4] = {"All","Ranking>0","|#Deltavz|<6cm","Ranking>0,|#Deltavz|<6cm"};
  TString legName2[4];
  for(int i=0; i<4; i++)
    {
      hTpcVsVpdVz[i]= (TH2F*)f->Get(Form("hTpcVsVpdVz%s",hname[i]));
      TH1F *htmp = (TH1F*)hTpcVsVpdVz[i]->ProjectionY();
      if(i==0) htmp->Rebin(10);
      list->Add(htmp);
      legName2[i] = htitle[i];
    }
  c = drawHistos(list,Form("%s_CompTpcVz",run_type),Form("%s: TPC vertex z distribution;vz_{TPC} (cm)",run_type),false,0,10,false,0,1,false,kTRUE,legName2,kTRUE,"",0.13,0.3,0.65,0.85,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_CompTpcVz.pdf",run_type,trigName));

  c = draw2D(hTpcVsVpdVz[0],Form("%s: vz_{TPC} vs. vz_{VPD} (%s);vz_{VPD} (cm); vz_{TPC} (cm)",run_type,htitle[0]));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVsVpdVz%s.png",run_type,trigName,hname[0]));

  TH1F *hTpcVzBad = (TH1F*)hTpcVsVpdVz[0]->ProjectionY("hTpcVzBad",1,100);
  hTpcVzBad->GetXaxis()->SetRangeUser(-5,5);
  c = draw1D(hTpcVzBad,Form("%s: TPC vertex z distribution (vz_{VPD} < -900 cm);vz_{TPC} (cm)",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVzBad.pdf",run_type,trigName));

  // Obtain vertex cut efficiency for VPDMB
  TH1F* hNeventAll    = (TH1F*)f->Get("hNeventAll");
  hNeventAll->Sumw2();
  TH1F* hNeventTpcVz  = (TH1F*)f->Get("hNeventTpcVz");
  TH1F* hNeventVzDiff = (TH1F*)f->Get("hNeventVzDiff");
  TH1F* hNeventVzDzRank = (TH1F*)f->Get("hNeventVzDzRank");

  TH1F *hTpcVzEff = (TH1F*)hNeventTpcVz->Clone("hTpcVzEff");
  hTpcVzEff->Divide(hNeventAll);
  hTpcVzEff->SetMarkerStyle(20);
  c = draw1D(hTpcVzEff,Form("%s: efficiency of TPC vz cut;RunId",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVzEff.pdf",run_type,trigName));

  TH1F *hVzDiffEff = (TH1F*)hNeventVzDiff->Clone("hVzDiffEff");
  hVzDiffEff->Divide(hNeventAll);
  hVzDiffEff->SetMarkerStyle(20);
  c = draw1D(hVzDiffEff,Form("%s: efficiency of vertex vz cut;RunId",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDiffEff.pdf",run_type,trigName));

  TH1F *hVzDzRankEff = (TH1F*)hNeventVzDzRank->Clone("hVzDzRankEff");
  hVzDzRankEff->Divide(hNeventAll);
  hVzDzRankEff->SetMarkerStyle(20);
  c = draw1D(hVzDzRankEff,Form("%s: efficiency of good vertex;RunId",run_type));
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDzRankEff.pdf",run_type,trigName));

  TH1F *hFinalEff = hVzDzRankEff;
  //TH1F *hFinalEff = hVzDiffEff;

  TH2F *hVtxEffVsBbcRate = new TH2F("hVtxEffVsBbcRate",Form("%s: vertex cut efficiency vs. BBC rate;BBC (kHz);Efficiency",run_type),300,0,3e3,100,0,1);
  for(int bin=1; bin<=hBbcRate->GetNbinsX(); bin++)
    {
      if(hBbcRate->GetBinContent(bin)<=0) continue;
      if(hBbcRate->GetBinCenter(bin)<=16050000) continue;
      hVtxEffVsBbcRate->Fill(hBbcRate->GetBinContent(bin)/1e3,hFinalEff->GetBinContent(bin));
    }
  c = draw2D(hVtxEffVsBbcRate);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDiffEffVsBBC.pdf",run_type,trigName));
  TF1 *funcEff = new TF1("funcEff","pol2",200,2800);
  if(run_type=="Run15_pp200")
    {
      TProfile *hpro = (TProfile*)hVtxEffVsBbcRate->ProfileX();
      hpro->SetMarkerStyle(24);
      hpro->Draw("sames");
      hpro->Fit(funcEff,"0");
      funcEff->SetLineColor(4);
      funcEff->SetLineStyle(2);
      c->cd();
      funcEff->Draw("sames");
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_VzDiffEffVsBBCFit.pdf",run_type,trigName));
    }

  // # of MTD events analyzed
  TFile *fdata = TFile::Open("Rootfiles/Run15DimuonStatisitcs.root","read");
  TH1F *hdata = 0x0;
  TProfile* hBbcData = 0x0;
  if(run_type=="Run15_pp200")
    {
      hdata = (TH1F*)fdata->Get("mhRunStat_Run15pp200");
      hBbcData = (TProfile*)fdata->Get("mhRunBbcRate_Run15pp200");
    }
  if(run_type=="Run15_pAu200") 
    {
      hdata = (TH1F*)fdata->Get("mhRunStat_Run15pAu200");
      hBbcData = (TProfile*)fdata->Get("mhRunBbcRate_Run15pAu200");
    }

  // Get equivalent # of MB events
  TFile *fin = 0;
  if(!saveHisto) fin = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  else           fin = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
  TH1F *hNeventMtd = (TH1F*)fin->Get("hNevents_dimuon");
  TH1F *hPrescaleMtd = (TH1F*)fin->Get("hPreScale_dimuon");
  TH1F *hLivetimeMtd = (TH1F*)fin->Get("hLiveTime_dimuon");
  TH1F *hNeventVpd = (TH1F*)fin->Get("hNevents_VPDMB-novtx");
  TH1F *hPrescaleVpd = (TH1F*)fin->Get("hPreScale_VPDMB-novtx");
  TH1F *hLivetimeVpd = (TH1F*)fin->Get("hLiveTime_VPDMB-novtx");

  double nMissEff = 0;
  double nMissVpdLive = 0;
  double nMtdEvents = 0;
  double nVpdEventsAcc = 0;
  for(int bin=1; bin<=hdata->GetNbinsX(); bin++)
    {
      int runnumber = atoi(hdata->GetXaxis()->GetBinLabel(bin));
      if(runnumber<16050000) continue;
      double nMtdEventAna = hdata->GetBinContent(bin);
      nMtdEvents += nMtdEventAna;

      int jbin = hNeventMtd->FindFixBin(runnumber);
      double nMtdEventAll = hNeventMtd->GetBinContent(jbin);
      double psMtd        = hPrescaleMtd->GetBinContent(jbin);
      double liveMtd      = hLivetimeMtd->GetBinContent(jbin);
      double nVpdEventAll = hNeventVpd->GetBinContent(jbin);
      double psVpd        = hPrescaleVpd->GetBinContent(jbin);
      double liveVpd      = hLivetimeVpd->GetBinContent(jbin);
      double vtxEff       = hFinalEff->GetBinContent(hFinalEff->FindFixBin(runnumber));
      if(vtxEff<=0)
	{
	  printf("[w] Missing efficiency for run %d\n",runnumber);
	  nMissEff += nMtdEventAna;
	  double bbcrate = hBbcData->GetBinContent(bin);
	  vtxEff = funcEff->Eval(bbcrate);
	  cout << bbcrate << "  " << vtxEff << endl;
	}
      if(liveVpd<=0)
	{
	  printf("[w] Missing VPDMB livetime for %d\n",runnumber);
	  nMissVpdLive += nMtdEventAna;
	  continue;
	}
      nVpdEventsAcc += nMtdEventAna/nMtdEventAll * nVpdEventAll * psVpd/liveVpd * liveMtd/psMtd * vtxEff;
      //printf("%1.0f = %1.0f, %1.0f = %1.0f, eff = %4.2f\n",runnumber,hNeventMtd->GetBinCenter(jbin),nMtdEventAna,nMtdEventAll,vtxEff); 
    }
  printf("[i] # of MTD event analyzed: %4.4e\n", nMtdEvents);
  printf("[i] # of MTD event without eff: %4.2e (%4.2f%%)\n", nMissEff, nMissEff/nMtdEvents*100);
  printf("[i] # of MTD event missing liveVPD: %4.2e (%4.2f%%)\n", nMissVpdLive, nMissVpdLive/nMtdEvents*100);
  printf("[i] # of equivalent MB events: %4.6e\n", nVpdEventsAcc);
}

//================================================
void makeHistoLumi(const int savePlot = 1, const int saveHisto = 1)
{
  //const char *trigName = "VPD-ZDC-novtx-mon";
  //const char *trigName = "VPDMB-novtx";
  //const char *trigName = "dimuon";
  const char *trigName = "BHT2-VPDMB-30";
  
  TH1F *hPreScale = new TH1F(Form("hPreScale_%s",trigName), Form("%s: per-scale per run;runId;Pre-scale",trigName), end_run-start_run+1, start_run-0.5, end_run+0.5);
  TH1F *hLiveTime = new TH1F(Form("hLiveTime_%s",trigName), Form("%s: live-time per run;runId;Live-time",trigName), end_run-start_run+1, start_run-0.5, end_run+0.5);
  TH1F *hNevents  = new TH1F(Form("hNevents_%s",trigName), Form("%s: # of events per run;runId;Nevts",trigName), end_run-start_run+1, start_run-0.5, end_run+0.5);
  TH1F *hEqMbEvt  = new TH1F(Form("hEqMbEvt_%s",trigName), Form("%s: # of equivalent MB events per run;runId;Nevts",trigName), end_run-start_run+1, start_run-0.5, end_run+0.5);


  ifstream flumi, fnevt;
  flumi.open(Form("Rootfiles/Luminosity/%s/lum_perrun_%s.txt",run_type,trigName));
  fnevt.open(Form("Rootfiles/Luminosity/%s/nev_perday_unix_%s.txt",run_type,trigName));
  char lumi[512], nevt[512];
  int nevt_old = 0;
  int nevt_new = 0;
  int totalevents = 0;
  double crosssection = 1;
  if(year==2014) crosssection = 6e6; // ub
  while(!flumi.eof() && !fnevt.eof())
    {
      flumi.getline(lumi,512);
      fnevt.getline(nevt,512);
      string slumi = lumi;
      string snevt = nevt;
      if(slumi.length()<=0) break;
      int counter = 0;
      int run, nmbevts, starttime;
      double prescale, livetime;
      while(slumi.find(" ")!=string::npos)
	{
	  string sub = slumi.substr(0,slumi.find(" "));
	  slumi.erase(0,slumi.find(" ")+1);
	  if(counter==0) run = atoi(sub.c_str());
	  if(counter==1) starttime = atoi(sub.c_str());
	  if(counter==4) nmbevts   = floor(atof(sub.c_str())*crosssection+0.5);
	  if(counter==5) prescale = atof(sub.c_str());
	  if(counter==6) livetime = atof(sub.c_str());
	  counter++;
	}
      int starttime2 = atoi(snevt.substr(0,12).c_str());
      nevt_new = atoi(snevt.erase(0,12).c_str());
      int nevents = nevt_new - nevt_old;
      nevt_old = nevt_new;
      totalevents += nevents;
      if(run>=start_run && run<=end_run)
	{
	  int bin = run - start_run + 1;
	  hPreScale->SetBinContent(bin, prescale);
	  hPreScale->SetBinError(bin, 1e-10);
	  hLiveTime->SetBinContent(bin, livetime);
	  hLiveTime->SetBinError(bin, 1e-10);
	  hEqMbEvt->SetBinContent(bin, nmbevts);
	  hEqMbEvt->SetBinError(bin, 1e-10);
	  if(starttime2==starttime)
	    {
	      if(year==2013 && trigName=="MTD-dimuon")
		{
		  if(run==14143081) nevents = 7384;
		  if(run==14143090) nevents = 7670;
		}
	      hNevents->SetBinContent(bin, nevents);
	      hNevents->SetBinError(bin, 1e-10);
	    }
	  else
	    {
	      printf("[e] mismatch %d != %d\n",starttime2,starttime);
	    }
	}
      //printf("[i] Run %d starts at %d has %d events\n",run,starttime,nevents);
    }
  printf("+++ Check event counts: %d =? %d +++\n",totalevents,nevt_new);

  hPreScale->SetMarkerStyle(21);
  c = draw1D(hPreScale);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_prescale.pdf",run_type,trigName));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_prescale.png",run_type,trigName));
    }

  hLiveTime->SetMarkerStyle(21);
  c = draw1D(hLiveTime);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_livetime.pdf",run_type,trigName));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_livetime.png",run_type,trigName));
    }

  hNevents->SetMarkerStyle(21);
  c = draw1D(hNevents);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_nevents.pdf",run_type,trigName));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_nevents.png",run_type,trigName));
    }

  hEqMbEvt->SetMarkerStyle(21);
  c = draw1D(hEqMbEvt);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_EqMbEvt.pdf",run_type,trigName));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_EqMbEvt.png",run_type,trigName));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
      hPreScale->Write("",TObject::kOverwrite);
      hLiveTime->Write("",TObject::kOverwrite);
      hNevents->Write("",TObject::kOverwrite);
      hEqMbEvt->Write("",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
void makeHistoData(const int savePlot = 1, const int saveHisto = 1)
{
  if(year!=2014)
    {
      printf("[e] Not suitable for %d\n",year);
      return;
    }
  const char *trigName = "VPD-ZDC-novtx-mon";
  const int max_vz = 100;
  const int max_dz = 3;
  const int max_vr = 2;

  const char *wName[2] = {"","_w"};
  const char *setupName[5] = {"all","prod","prod_low","prod_mid","prod_high"};
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  const int nCentBins = 17;
  const int centBins_low[nCentBins]  = {5,  13, 9,  5, 15, 11, 5,  13, 11, 9,  7, 5, 1,  1, 3, 1, 1};
  const int centBins_high[nCentBins] = {16, 16, 12, 8, 16, 14, 10, 14, 12, 10, 8, 6, 16, 4, 4, 2, 1};
  const char *cent_Name[nCentBins] = {"0-60","0-20","20-40","40-60","0-10","10-30","30-60","10-20","20-30","30-40","40-50","50-60","0-80","60-80","60-70","70-80","75-80"};
  const char *cent_Title[nCentBins] = {"0060","0020","2040","4060","0010","1030","3060","1020","2030","3040","4050","5060","0080","6080","6070","7080","7580"};
  TH1F *hNEvents[5][2];
  TH1F *hTpcVz[5][nCentBins];
  TH1F *hDiffVz[5][nCentBins];
  TH1F *hTpcVr[5][nCentBins];
  TH1F *hNEventsAll[5][nCentBins][2];
  TH1F *hNEventsVr[5][nCentBins][2];
  TH1F *hNEventsVz[5][nCentBins][2];
  TH1F *hNEventsAcc[5][nCentBins][2];
  TH1F *hCent[5][2];

  TH2F *hTpcVzVsRun[nCentBins];
  TH2F *hTpcVrVsRun[nCentBins];
  TH2F *hDiffVzVsRun[nCentBins];
  TH2F *hCentVsRun[2];


  const char *lumiType = "prod_high";
  TFile *fMB = TFile::Open(Form("./output/Run14_AuAu200.MB.VtxEff.%s.root",lumiType),"read");
  THnSparseF *hn[2];
  hn[0] = (THnSparseF*)fMB->Get("mhMbEvtEff");
  hn[1] = (THnSparseF*)fMB->Get("mhMbEvtEffWeight");

  for(int i=0; i<2; i++)
    {
      hNEvents[0][i] = (TH1F*)hn[i]->Projection(0);
      hNEvents[0][i]->SetName(Form("NEvents_%s%s",setupName[0],wName[i]));

      for(int j=0; j<nCentBins; j++)
	{
	  hn[i]->GetAxis(5)->SetRange(centBins_low[j],centBins_high[j]);
	  hNEventsAll[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsAll[0][j][i]->SetName(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	  if(i==0)
	    {
	      hTpcVz[0][j] = (TH1F*)hn[i]->Projection(1);
	      hTpcVz[0][j]->SetName(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hTpcVzVsRun[j] = (TH2F*)hn[i]->Projection(1,0);
	      hTpcVzVsRun[j]->SetName(Form("TpcVzVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      //draw2D(hTpcVzVsRun[j]);

	      hTpcVr[0][j] = (TH1F*)hn[i]->Projection(3);
	      hTpcVr[0][j]->SetName(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hTpcVrVsRun[j] = (TH2F*)hn[i]->Projection(3,0);
	      hTpcVrVsRun[j]->SetName(Form("TpcVrVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	      hDiffVz[0][j] = (TH1F*)hn[i]->Projection(2);
	      hDiffVz[0][j]->SetName(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	      hDiffVzVsRun[j] = (TH2F*)hn[i]->Projection(2,0);
	      hDiffVzVsRun[j]->SetName(Form("DiffVzVsRun_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));
	    }

	  hn[i]->GetAxis(3)->SetRangeUser(-1,max_vr-1e-4);
	  hNEventsVr[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsVr[0][j][i]->SetName(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(1)->SetRangeUser(-1*max_vz+1e-4,max_vz-1e-4);
	  hNEventsVz[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsVz[0][j][i]->SetName(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(2)->SetRangeUser(-1*max_dz+1e-4,max_dz-1e-4);
	  hNEventsAcc[0][j][i] = (TH1F*)hn[i]->Projection(0);
	  hNEventsAcc[0][j][i]->SetName(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[0],wName[i]));

	  hn[i]->GetAxis(5)->SetRange(0,-1);
	  if(j==0)
	    {
	      hCent[0][i] = (TH1F*)hn[i]->Projection(5);
	      hCent[0][i]->SetName(Form("Cent_%s%s",setupName[0],wName[i]));
	      hCentVsRun[i] = (TH2F*)hn[i]->Projection(5,0);
	      hCentVsRun[i]->SetName(Form("CentVsRun_%s%s",setupName[0],wName[i]));
	    }

	  hn[i]->GetAxis(1)->SetRange(0,-1);
	  hn[i]->GetAxis(2)->SetRange(0,-1);
	  hn[i]->GetAxis(3)->SetRange(0,-1);

	}
    }


  for(int k=0; k<4; k++)
    {
      for(int i=0; i<2; i++)
	{
	  hNEvents[k+1][i] = (TH1F*)hNEvents[0][i]->Clone(Form("NEvents_%s%s",setupName[k+1],wName[i]));
	  hNEvents[k+1][i]->Reset();
	  for(int j=0; j<nCentBins; j++)
	    {
	      if(i==0)
		{
		  hTpcVz[k+1][j] = (TH1F*)hTpcVz[0][j]->Clone(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hTpcVz[k+1][j]->Reset();

		  hTpcVr[k+1][j] = (TH1F*)hTpcVr[0][j]->Clone(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hTpcVr[k+1][j]->Reset();

		  hDiffVz[k+1][j] = (TH1F*)hDiffVz[0][j]->Clone(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  hDiffVz[k+1][j]->Reset();
		}
	      hNEventsAll[k+1][j][i] = (TH1F*)hNEventsAll[0][j][i]->Clone(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsAll[k+1][j][i]->Reset();

	      hNEventsVr[k+1][j][i] = (TH1F*)hNEventsVr[0][j][i]->Clone(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsVr[k+1][j][i]->Reset();

	      hNEventsVz[k+1][j][i] = (TH1F*)hNEventsVz[0][j][i]->Clone(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsVz[k+1][j][i]->Reset();

	      hNEventsAcc[k+1][j][i] = (TH1F*)hNEventsAcc[0][j][i]->Clone(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      hNEventsAcc[k+1][j][i]->Reset();

	      if(j==0)
		{
		  hCent[k+1][i] = (TH1F*)hCent[0][i]->Clone(Form("Cent_%s%s",setupName[k+1],wName[i]));
		  hCent[k+1][i]->Reset();
		}
	    }
	}
    }

  for(int k=0; k<4; k++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s_2014.list",run_type,trgSetupName[k]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int bin = hNEvents[0][0]->FindFixBin(runnumber);

	  for(int i=0; i<2; i++)
	    {
	      hNEvents[k+1][i]->SetBinContent(bin,hNEvents[0][i]->GetBinContent(bin));
	      hNEvents[k+1][i]->SetBinError(bin,hNEvents[0][i]->GetBinError(bin));
	      for(int j=0; j<nCentBins; j++)
		{
		  if(i==0)
		    {
		      TH1F *htmp = (TH1F*)hTpcVzVsRun[j]->ProjectionY(Form("%s_%d",hTpcVzVsRun[j]->GetName(),bin),bin,bin);
		      hTpcVz[k+1][j]->Add(htmp);

		      htmp = (TH1F*)hTpcVrVsRun[j]->ProjectionY(Form("%s_%d",hTpcVrVsRun[j]->GetName(),bin),bin,bin);
		      hTpcVr[k+1][j]->Add(htmp);

		      htmp = (TH1F*)hDiffVzVsRun[j]->ProjectionY(Form("%s_%d",hDiffVzVsRun[j]->GetName(),bin),bin,bin);
		      hDiffVz[k+1][j]->Add(htmp);
		    }
	  
		  hNEventsAll[k+1][j][i]->SetBinContent(bin,hNEventsAll[0][j][i]->GetBinContent(bin));
		  hNEventsAll[k+1][j][i]->SetBinError(bin,hNEventsAll[0][j][i]->GetBinError(bin));

		  hNEventsVr[k+1][j][i]->SetBinContent(bin,hNEventsVr[0][j][i]->GetBinContent(bin));
		  hNEventsVr[k+1][j][i]->SetBinError(bin,hNEventsVr[0][j][i]->GetBinError(bin));

		  hNEventsVz[k+1][j][i]->SetBinContent(bin,hNEventsVz[0][j][i]->GetBinContent(bin));
		  hNEventsVz[k+1][j][i]->SetBinError(bin,hNEventsVz[0][j][i]->GetBinError(bin));

		  hNEventsAcc[k+1][j][i]->SetBinContent(bin,hNEventsAcc[0][j][i]->GetBinContent(bin));
		  hNEventsAcc[k+1][j][i]->SetBinError(bin,hNEventsAcc[0][j][i]->GetBinError(bin));
		  if(j==0)
		    {
		      htmp = (TH1F*)hCentVsRun[i]->ProjectionY(Form("%s_%d",hCentVsRun[i]->GetName(),bin),bin,bin);
		      hCent[k+1][i]->Add(htmp);
		    }
		}
	    }
	}
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.%s.root",run_type,lumiType),"update");
      for(int k=0; k<5; k++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hNEvents[k][i]->SetTitle("");
	      hNEvents[k][i]->Write("",TObject::kOverwrite);
	      for(int j=0; j<nCentBins; j++)
		{
		  hNEventsAll[k][j][i]->SetTitle("");
		  hNEventsVr[k][j][i]->SetTitle("");
		  hNEventsVz[k][j][i]->SetTitle("");
		  hNEventsAcc[k][j][i]->SetTitle("");
		  hNEventsAll[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVr[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVz[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsAcc[k][j][i]->Write("",TObject::kOverwrite);
		  if(j==0)
		    {
		      hCent[k][i]->SetTitle("");
		      hCent[k][i]->Write("",TObject::kOverwrite);
		    }
		  if(i==0)
		    {
		      hTpcVz[k][j]->SetTitle("");
		      hTpcVr[k][j]->SetTitle("");
		      hDiffVz[k][j]->SetTitle("");
		      hTpcVz[k][j]->Write("",TObject::kOverwrite);
		      hTpcVr[k][j]->Write("",TObject::kOverwrite);
		      hDiffVz[k][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}


//================================================
void mergeHisto(const int saveHisto = 1)
{
  if(year!=2014)
    {
      printf("[e] Not suitable for %d\n",year);
      return;
    }

  const char *wName[2] = {"","_w"};
  const char *setupName[5] = {"all","prod","prod_low","prod_mid","prod_high"};
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  const int nCentBins = 17;
  const char *cent_Title[nCentBins] = {"0060","0020","2040","4060","0010","1030","3060","1020","2030","3040","4050","5060","0080","6080","6070","7080","7580"};
  TH1F *hNEvents[5][2];
  TH1F *hTpcVz[5][nCentBins];
  TH1F *hDiffVz[5][nCentBins];
  TH1F *hTpcVr[5][nCentBins];
  TH1F *hNEventsAll[5][nCentBins][2];
  TH1F *hNEventsVr[5][nCentBins][2];
  TH1F *hNEventsVz[5][nCentBins][2];
  TH1F *hNEventsAcc[5][nCentBins][2];
  TH1F *hCent[5][2];

  TFile *f[2];
  f[0] = TFile::Open(Form("Rootfiles/%s.Luminosity.prod_low.root",run_type),"read");
  f[1] = TFile::Open(Form("Rootfiles/%s.Luminosity.prod_high.root",run_type),"read");
  TFile *fCurr = 0x0;
  for(int k=0; k<4; k++)
    {
    if(k<3) fCurr = f[0];
    else    fCurr = f[1];
      for(int i=0; i<2; i++)
	{
	  hNEvents[k+1][i] = (TH1F*)fCurr->Get(Form("NEvents_%s%s",setupName[k+1],wName[i]));
	  if(k==0) hNEvents[0][i] = (TH1F*)hNEvents[k+1][i]->Clone(Form("NEvents_%s%s",setupName[k],wName[i]));
	  else     hNEvents[0][i]->Add(hNEvents[k+1][i]);

	  for(int j=0; j<nCentBins; j++)
	    {
	      if(i==0)
		{
		  hTpcVz[k+1][j]  = (TH1F*)fCurr->Get(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  if(k==0) hTpcVz[0][j] = (TH1F*)hTpcVz[k+1][j]->Clone(Form("TpcVz_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
		  else     hTpcVz[0][j]->Add(hTpcVz[k+1][j]);

		  hTpcVr[k+1][j]  = (TH1F*)fCurr->Get(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  if(k==0) hTpcVr[0][j] = (TH1F*)hTpcVr[k+1][j]->Clone(Form("TpcVr_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
		  else     hTpcVr[0][j]->Add(hTpcVr[k+1][j]);

		  hDiffVz[k+1][j] = (TH1F*)fCurr->Get(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
		  if(k==0) hDiffVz[0][j] = (TH1F*)hDiffVz[k+1][j]->Clone(Form("TpcVpdDz_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
		  else     hDiffVz[0][j]->Add(hDiffVz[k+1][j]);
		}
	      hNEventsAll[k+1][j][i] = (TH1F*)fCurr->Get(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      if(k==0) hNEventsAll[0][j][i] = (TH1F*)hNEventsAll[k+1][j][i]->Clone(Form("NEvents_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
	      else     hNEventsAll[0][j][i]->Add(hNEventsAll[k+1][j][i]);

	      hNEventsVr[k+1][j][i]  = (TH1F*)fCurr->Get(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      if(k==0) hNEventsVr[0][j][i] = (TH1F*)hNEventsVr[k+1][j][i]->Clone(Form("NEvents_VrCut_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
	      else     hNEventsVr[0][j][i]->Add(hNEventsVr[k+1][j][i]);

	      hNEventsVz[k+1][j][i]  = (TH1F*)fCurr->Get(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      if(k==0) hNEventsVz[0][j][i] = (TH1F*)hNEventsVz[k+1][j][i]->Clone(Form("NEvents_VrVzCut_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
	      else     hNEventsVz[0][j][i]->Add(hNEventsVz[k+1][j][i]);

	      hNEventsAcc[k+1][j][i] = (TH1F*)fCurr->Get(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[k+1],wName[i]));
	      if(k==0) hNEventsAcc[0][j][i] = (TH1F*)hNEventsAcc[k+1][j][i]->Clone(Form("NEvents_VrVzDzCut_cent%s_%s%s",cent_Title[j],setupName[k],wName[i]));
	      else     hNEventsAcc[0][j][i]->Add(hNEventsAcc[k+1][j][i]);

	      if(j==0)
		{
		  hCent[k+1][i] = (TH1F*)fCurr->Get(Form("Cent_%s%s",setupName[k+1],wName[i]));
		  if(k==0) hCent[0][i] = (TH1F*)hCent[k+1][i]->Clone(Form("Cent_%s%s",setupName[k],wName[i]));
		  else     hCent[0][i]->Add(hCent[k+1][i]);
		}
	    }
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"recreate");
      for(int k=0; k<5; k++)
	{
	  for(int i=0; i<2; i++)
	    {
	      hNEvents[k][i]->SetTitle("");
	      hNEvents[k][i]->Write("",TObject::kOverwrite);
	      for(int j=0; j<nCentBins; j++)
		{
		  hNEventsAll[k][j][i]->SetTitle("");
		  hNEventsVr[k][j][i]->SetTitle("");
		  hNEventsVz[k][j][i]->SetTitle("");
		  hNEventsAcc[k][j][i]->SetTitle("");
		  hNEventsAll[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVr[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsVz[k][j][i]->Write("",TObject::kOverwrite);
		  hNEventsAcc[k][j][i]->Write("",TObject::kOverwrite);
		  if(j==0)
		    {
		      hCent[k][i]->SetTitle("");
		      hCent[k][i]->Write("",TObject::kOverwrite);
		    }
		  if(i==0)
		    {
		      hTpcVz[k][j]->SetTitle("");
		      hTpcVr[k][j]->SetTitle("");
		      hDiffVz[k][j]->SetTitle("");
		      hTpcVz[k][j]->Write("",TObject::kOverwrite);
		      hTpcVr[k][j]->Write("",TObject::kOverwrite);
		      hDiffVz[k][j]->Write("",TObject::kOverwrite);
		    }
		}
	    }
	}
    }
}

