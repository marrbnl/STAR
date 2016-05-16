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

  //Lumi2013();
  Lumi2014();
  //MB();
  //MTD();
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
  TString legName2[4];
  for(int i=0; i<4; i++)
    {
      legName2[i] = Form("AuAu_200_%s_2014",trgSetupName[i]);
    }
  TString legName3[2] = {"w/o weights","w/ weights"};
  TString legName4[nCentBins];
  for(int i=0; i<nCentBins; i++) legName4[i] = Form("%s%%",cent_Name[i]);
  const int max_vz = 100;
  const int max_dz = 3;
  const int max_vr = 2;


  // MTD triggers
  TFile *fMtd = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.root"),"read");
  TH2F *hCentWeight = (TH2F*)fMtd->Get("mhCentWeight_di_mu");
  c = draw2D(hCentWeight,"Event weights vs. multiplicity");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/EventWeight.pdf",run_type));

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
    }
  c = drawHistos(list,"MbCent","VPD-ZDC-novtx-mon: distribution of centrality bins;cent",true,0,17,true,0,1.2*hMbCent[1]->GetMaximum(),kFALSE,kTRUE,legName3,kTRUE,"",0.5,0.7,0.2,0.35,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_CentWeight.pdf",run_type));
  list->Clear();

  // Vertex distribution
  TH1F *hTpcVz[4];
  for(int i=0; i<4; i++)
    {
      hTpcVz[i] = (TH1F*)fMB->Get(Form("TpcVz_cent0060_%s",setupName[i]));
      hTpcVz[i]->Scale(1./hTpcVz[i]->Integral());
      list->Add(hTpcVz[i]);
    }
  c = drawHistos(list,"TpcVz","VPD-ZDC-novtx-mon: TPC vz distribution (0-60%);vz_{TPC} (cm)",false,0,30,true,0,0.02,kFALSE,kTRUE,legName2,kTRUE,"0-60%",0.18,0.38,0.68,0.88,false);
  TLine *line = GetLine(-1*max_vz,0,-1*max_vz,0.8*hTpcVz[3]->GetMaximum());
  line->Draw();
  line = GetLine(max_vz,0,max_vz,0.8*hTpcVz[3]->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVz.pdf",run_type));
  list->Clear();

  TH1F *hTpcVr[4];
  for(int i=0; i<4; i++)
    {
      hTpcVr[i] = (TH1F*)fMB->Get(Form("TpcVr_cent0060_%s",setupName[i]));
      hTpcVr[i]->Scale(1./hTpcVr[i]->Integral());
      list->Add(hTpcVr[i]);
    }
  c = drawHistos(list,"TpcVr","VPD-ZDC-novtx-mon: TPC vr distribution (0-60%);vr_{TPC} (cm)",false,0,30,true,0.00001,1000,true,kTRUE,legName2,kTRUE,"0-60%",0.18,0.38,0.68,0.88,false);
  line = GetLine(max_vr,0,max_vr,0.8*hTpcVr[3]->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVr.pdf",run_type));
  list->Clear();

  TH1F *hTpcVpdDz[4];
  for(int i=0; i<4; i++)
    {
      hTpcVpdDz[i] = (TH1F*)fMB->Get(Form("TpcVpdDz_cent0060_%s",setupName[i]));
      hTpcVpdDz[i]->Scale(1./hTpcVpdDz[i]->Integral());
      list->Add(hTpcVpdDz[i]);
    }
  c = drawHistos(list,"TpcVpdDz","VPD-ZDC-novtx-mon:#Deltaz distribution for vertices (0-60%);vz_{TPC}-vz_{VPD} (cm)",true,-10,10,true,1e-6,1e3,true,kTRUE,legName2,kTRUE,"0-60%",0.18,0.38,0.68,0.88,false);
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

  list->Clear();
  TH1F *hVtxEff[nCentBins][3];
  TCanvas *c = new TCanvas("MB_VtxCutEff","MB_VtxCutEff",1100,650);
  c->Divide(2,2);
  for(int i=0; i<nCentBins; i++)
    {
      list->Clear();
      c->cd(i+1);
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
  for(int j=0; j<nCentBins; j++)
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
  TFile *frf = TFile::Open("./output/Run14.AuAu200.RejectFactor.root","read");
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
	  hRFLumi[i]->SetBinContent(bin,hRF->GetBinContent(bin));
	  hRFLumi[i]->SetBinError(bin,hRF->GetBinError(bin));
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
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_FitRejectFactorLumi.pdf",run_type));
  c = drawHistos(list,"MtdRFLumi",Form("Dimuon: fraction of events accepted during two-pass reconstruction"),kFALSE,0,30,kTRUE,0.35,0.55,kFALSE,kTRUE,legName2,kTRUE,"trgsetupname",0.3,0.5,0.16,0.36,true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_RejectFactorLumi.pdf",run_type));
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
  for(int i=0; i<nCentBins; i++)
    {
     hNEventsRaw[i][0] = (TH1F*)fMB->Get(Form("NEvents_cent%s_all",cent_Title[i]));
     hNevtVtxCutWeight[i][0] = (TH1F*)fMB->Get(Form("NEvents_VrVzDzCut_cent%s_all_w",cent_Title[i]));
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
	  double evtcent = hNEventsRaw[i][0]->GetBinContent(bin);
	  double evtw = hNevtVtxCutWeight[i][0]->GetBinContent(bin);
	  if(rawevt<=0) continue;
	  double value = eqevt/rawevt*evtw;
	  double error = evtw<=0? 0 : value * 1./sqrt(evtw);
	  hEqMbEvtVtxCutWeight[i]->SetBinContent(ibin, value);
	  hEqMbEvtVtxCutWeight[i]->SetBinError(ibin, error);
	}
    }

  // Fit event fraction after vertex and centrality selections
  TF1 *funcMbVtxEff[4][4];
  c = new TCanvas("Fit_MB_VtxEff","Fit_MB_VtxEff",1200,700);
  c->Divide(4,4);
  for(int i=0; i<4; i++)
    {
      hNEventsRawAll[i+1] = (TH1F*)fMB->Get(Form("NEvents_%s",setupName[i]));
      for(int j=0; j<nCentBins; j++)
	{
	  hNevtVtxCutWeight[j][i+1] = (TH1F*)fMB->Get(Form("NEvents_VrVzDzCut_cent%s_%s_w",cent_Title[j],setupName[i]));
	  TH1F *htmp = (TH1F*)hNevtVtxCutWeight[j][i+1]->Clone(Form("Fit_%s",hNevtVtxCutWeight[j][i+1]->GetName()));
	  htmp->Divide(hNEventsRawAll[i+1]);
	  if(i==3)
	    {
	      double run = 15119021;
	      int bin = htmp->FindFixBin(run);
	      htmp->SetBinContent(bin,0);
	      htmp->SetBinError(bin,0);
	    }
	  funcMbVtxEff[j][i] = new TF1(Form("Func_%s",hNevtVtxCutWeight[j][i+1]->GetName()),"pol0",start_run,end_run);
	  htmp->Fit(funcMbVtxEff[j][i],"IR0Q");
	  c->cd(i*4+j+1);
	  htmp->SetMarkerStyle(29);
	  htmp->GetYaxis()->SetRangeUser(0,1);
	  htmp->Draw();
	  TPaveText *t1 = GetTitleText(Form("%s: %s%%",setupName[i],cent_Name[j]),0.08);
	  t1->Draw();
	  funcMbVtxEff[j][i]->SetLineColor(2);
	  funcMbVtxEff[j][i]->Draw("sames");
	}
    }
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/MB_EvtFractionWithVtxCutCent.pdf",run_type));
  //setupName
  
  // check vertex efficiency
  printf("+++ check vertex efficiency +++\n");
  int nVE = 0;
  double nVEevts = 0;
  for(int ibin=1; ibin<=hEvtRunAcc->GetNbinsX(); ibin++)
    {
      double run = hEvtRunAcc->GetBinCenter(ibin);
      if(hEvtRunAcc->GetBinContent(ibin)>0)
	{
	  for(int i=0; i<4; i++)
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
  for(int i=0; i<nCentBins; i++)
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
void makeHistoLumi(const int savePlot = 0, const int saveHisto = 0)
{
  const char *trigName = "VPD-ZDC-novtx-mon";
  //const char *trigName = "dimuon";
  
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
void makeHistoData(const int savePlot = 0, const int saveHisto = 1)
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

  TFile *fMB = TFile::Open(Form("./output/Run14.AuAu200.MB.VtxEff.root"),"read");
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
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
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
void MTD(const int savePlot = 1)
{
  TH1F *hRunId[2];
  TH1F *hZdcRate[2];
  for(int j=0; j<3; j++)
    {
      TFile *f = TFile::Open(Form("output/temp/%d.root",10+j),"read");
      THnSparseF *hnQA = (THnSparseF*)f->Get("mhEvtWithTwoMuons_di_mu");
      hnQA->SetName(Form("%s_%d",hnQA->GetName(),j));
      for(int i=0; i<2; i++)
	{
	  if(i==1) hnQA->GetAxis(5)->SetRange(2,2);
	  TH1F *htmp = (TH1F*)hnQA->Projection(0);
	  htmp->SetName(Form("RunId_%d_%d",j,i));
	  if(j==0) hRunId[i] = (TH1F*)htmp->Clone(Form("RunId_%d",i));
	  else     hRunId[i]->Add(htmp);
	  
	  htmp = (TH1F*)hnQA->Projection(3);
	  htmp->SetName(Form("hZdcRate_%d_%d",j,i));
	  if(j==0) hZdcRate[i] = (TH1F*)htmp->Clone(Form("hZdcRate_%d",i));
	  else     hZdcRate[i]->Add(htmp);
	}
    }
  for(int i=0; i<2; i++)
    {
      hRunId[i]->Sumw2();
      hZdcRate[i]->Sumw2();
    }
  TH1F *hRunIdEff = (TH1F*)hRunId[1]->Clone("RunId_Eff");
  hRunIdEff->Divide(hRunId[0]);
  c = draw1D(hRunIdEff);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_RunId.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_RunId.png",run_type));
    }

  TH1F *hZdcRateEff = (TH1F*)hZdcRate[1]->Clone("ZdcRate_Eff");
  hZdcRateEff->Divide(hZdcRate[0]);
  hZdcRateEff->GetXaxis()->SetRangeUser(0,100);
  hZdcRateEff->SetMarkerStyle(21);
  c = draw1D(hZdcRateEff,"Fraction of events with two muon candidates;ZdcRate (kHz);Fraction");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_ZdcRate.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_ZdcRate.png",run_type));
    }
}

//================================================
void MB(const int savePlot = 1)
{
  //gStyle->SetOptStat(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *fData = TFile::Open("output/Pico.Run14.AuAu200.jpsi.TightTof.root","read");
  TH1F *hVtxZMtd = (TH1F*)fData->Get("mhTpcVz_di_mu");
  hVtxZMtd->Scale(1./hVtxZMtd->GetBinContent(hVtxZMtd->FindBin(0)));
  hVtxZMtd->SetLineColor(4);
  TH2F *hVzDiffVsTpcVz = (TH2F*)fData->Get("mhVzDiffVsTpcVz_di_mu");
  TH1F *hDzMtd = (TH1F*)hVzDiffVsTpcVz->ProjectionY(Form("hDzMtd"),hVzDiffVsTpcVz->GetXaxis()->FindBin(-100),hVzDiffVsTpcVz->GetXaxis()->FindBin(100));


  TFile *fMB = TFile::Open("output/Pico.Run14.AuAu200.MB.root","read");
  THnSparseF *hnCentQA = (THnSparseF*)fMB->Get("mhEventCent_qa_mb");
  TH1F *hVertexMB = (TH1F*)hnCentQA->Projection(4);
  hVertexMB->Scale(1./hVertexMB->GetBinContent(hVertexMB->FindBin(0)));
  hVertexMB->SetLineColor(2);
  //hVertexMB->Draw("sames HIST");
  //printf("MB 2: %2.2f%%\n",hVertexMB->Integral(hVertexMB->FindBin(-100),hVertexMB->FindBin(100))/hVertexMB->Integral(0,-1)*100);

  const char *name[2] = {"Mid","Low"};
  TFile *f[2];
  TH1F *hVtxZDimuonNoCut[2];
  TH1F *hDzDimuonVzCut[2];
  TH1F *hVtxZNoCut[2];
  TH1F *hVtxDzVzCut[2];
  TH1F *hZdcRateAllWithCut[2];
  TH1F *hZdcRateDimuonWithCut[2];
  TH1F *hRFvsRate[2];

  for(int i=0; i<1; i++)
    {
      printf("+++ process production_%s\n",name[i]);
      f[i] = TFile::Open(Form("output/Run14.AuAu200.MB.Prod%s.root",name[i]),"read");
      THnSparseF *hnLumiQA = (THnSparseF*)f[i]->Get("mhLumiQA_mb");
      hnLumiQA->SetName(Form("%s_%s",hnLumiQA->GetName(),name[i]));
      if(i==0) hnLumiQA->GetAxis(0)->SetRangeUser(15119017,15163026);
      //if(i==0) hnLumiQA->GetAxis(0)->SetRangeUser(15163026,15167014);

      hnLumiQA->GetAxis(10)->SetRange(2,2);
      hVtxZDimuonNoCut[i] = (TH1F*)hnLumiQA->Projection(6);
      hVtxZDimuonNoCut[i]->SetName(Form("hVtxZDimuonNoCut_%s",name[i]));
      cout << "Dimuon before cuts: " << hVtxZDimuonNoCut[i]->Integral(0,-1) << endl;
      hnLumiQA->GetAxis(6)->SetRangeUser(-100,100);
      hDzDimuonVzCut[i] = (TH1F*)hnLumiQA->Projection(7);
      hDzDimuonVzCut[i]->SetName(Form("hDzDimuonVzCut_%s",name[i]));
      hnLumiQA->GetAxis(6)->SetRange(0,-1);
      hnLumiQA->GetAxis(10)->SetRange(0,-1);

      int eventAll = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventAll << endl;

      hnLumiQA->GetAxis(4)->SetRange(2,2);
      int eventVtx = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventVtx << endl;
      
      hnLumiQA->GetAxis(5)->SetRange(2,2);
      int eventGoodVtx = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventGoodVtx << endl;

      /*
	hnLumiQA->GetAxis(8)->SetRange(2,2);
	int eventVptTacHigh = hnLumiQA->Projection(3)->Integral(0,-1);
	cout << eventVptTacHigh << endl;

	hnLumiQA->GetAxis(9)->SetRange(2,2);
	int eventVptTacLow = hnLumiQA->Projection(3)->Integral(0,-1);
	cout << eventVptTacLow << endl;
      */

      hVtxZNoCut[i] = (TH1F*)hnLumiQA->Projection(6);
      cout << hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100)) << endl;
      hVtxZNoCut[i]->SetName(Form("hVtxZNoCut_%s",name[i]));
      hnLumiQA->GetAxis(6)->SetRangeUser(-100,100);
      int eventTpcVz = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventTpcVz << endl;

      hVtxDzVzCut[i] = (TH1F*)hnLumiQA->Projection(7);
      hVtxDzVzCut[i]->SetName(Form("hVtxDzVzCut_%s",name[i]));
      hnLumiQA->GetAxis(7)->SetRangeUser(-3,3);

      hZdcRateAllWithCut[i] = (TH1F*)hnLumiQA->Projection(3);
      hZdcRateAllWithCut[i]->SetName(Form("hZdcRateAllWithCut_%s",name[i]));
      int eventVzDiff = hZdcRateAllWithCut[i]->Integral(0,-1);
      cout << eventVzDiff << endl;

      hnLumiQA->GetAxis(10)->SetRange(2,2);
      hZdcRateDimuonWithCut[i] = (TH1F*)hnLumiQA->Projection(3);
      hZdcRateDimuonWithCut[i]->SetName(Form("hZdcRateDimuonWithCut_%s",name[i]));
      int eventDimuon = hZdcRateDimuonWithCut[i]->Integral(0,-1);
      cout << eventDimuon << endl;

      // vz distribution
      hVtxZNoCut[i]->Scale(1./hVtxZNoCut[i]->GetBinContent(hVtxZNoCut[i]->FindBin(0)));
      hVtxZNoCut[i]->GetYaxis()->SetRangeUser(0.01,10);
      c = draw1D(hVtxZNoCut[i],Form("TPC vz distribution (production_%s)",name[i]),kTRUE,kFALSE);
      printf("MB 1: %2.2f%%\n",hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100))/hVtxZNoCut[i]->Integral(0,-1)*100);

      hVtxZDimuonNoCut[i]->SetMarkerStyle(20);
      hVtxZDimuonNoCut[i]->Sumw2();
      hVtxZDimuonNoCut[i]->Scale(1./hVtxZDimuonNoCut[i]->GetBinContent(hVtxZDimuonNoCut[i]->FindBin(0)));
      hVtxZDimuonNoCut[i]->SetLineColor(6);
      hVtxZDimuonNoCut[i]->SetMarkerColor(6);
      hVtxZDimuonNoCut[i]->Draw("sames P");
  
      hVtxZMtd->Draw("sames HIST");
      printf("Dimuon 1: %2.2f%%\n",hVtxZMtd->Integral(hVtxZMtd->FindBin(-100),hVtxZMtd->FindBin(100))/hVtxZMtd->Integral(0,-1)*100);

      TLegend *leg = new TLegend(0.15,0.68,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hVtxZNoCut[i],Form("VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100))/hVtxZNoCut[i]->Integral(0,-1)*100),"L");
      leg->AddEntry(hVtxZDimuonNoCut[i],Form("dimuon in VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxZDimuonNoCut[i]->Integral(hVtxZDimuonNoCut[i]->FindBin(-100),hVtxZDimuonNoCut[i]->FindBin(100))/hVtxZDimuonNoCut[i]->Integral(0,-1)*100),"P");
      leg->AddEntry(hVtxZMtd,Form("Production_high dimuon trigger: %2.1f%%",hVtxZMtd->Integral(hVtxZMtd->FindBin(-100),hVtxZMtd->FindBin(100))/hVtxZMtd->Integral(0,-1)*100),"L");
      leg->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVtxZ.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVtxZ.png",run_type,name[i]));
	}

      // dz distribution
      hVtxDzVzCut[i]->Scale(1./hVtxDzVzCut[i]->GetBinContent(hVtxDzVzCut[i]->FindBin(0)));
      hVtxDzVzCut[i]->GetYaxis()->SetRangeUser(1e-5,1e2);
      c = draw1D(hVtxDzVzCut[i],Form("vz_{TPC}-vz_{VPD} distribution (production_%s)",name[i]),kTRUE,kFALSE);
      hDzMtd->Scale(1./hDzMtd->GetBinContent(hDzMtd->FindBin(0)));
      hDzMtd->SetLineColor(4);
      hDzMtd->Draw("sames HIST");
      hDzDimuonVzCut[i]->Sumw2();
      hDzDimuonVzCut[i]->Scale(1./hDzDimuonVzCut[i]->GetBinContent(hDzDimuonVzCut[i]->FindBin(0)));
      hDzDimuonVzCut[i]->SetLineColor(6);
      hDzDimuonVzCut[i]->SetMarkerStyle(20);
      hDzDimuonVzCut[i]->SetMarkerColor(6);
      hDzDimuonVzCut[i]->Draw("sames P");
      TLegend *leg = new TLegend(0.15,0.68,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hVtxDzVzCut[i],Form("VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxDzVzCut[i]->Integral(hVtxDzVzCut[i]->FindBin(-3),hVtxDzVzCut[i]->FindBin(3))/hVtxDzVzCut[i]->Integral(0,-1)*100),"L");
      leg->AddEntry(hDzDimuonVzCut[i],Form("dimuon in VPD-ZDC-NoVtx-MB: %2.1f%%",hDzDimuonVzCut[i]->Integral(hDzDimuonVzCut[i]->FindBin(-3),hDzDimuonVzCut[i]->FindBin(3))/hDzDimuonVzCut[i]->Integral(0,-1)*100),"P");
      leg->AddEntry(hDzMtd,Form("Production_high dimuon trigger: %2.1f%%",hDzMtd->Integral(hDzMtd->FindBin(-3),hDzMtd->FindBin(3))/hDzMtd->Integral(0,-1)*100),"L");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVpdDz.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVpdDz.png",run_type,name[i]));
	}


      // vs zdc rate
      hZdcRateAllWithCut[i]->Sumw2();
      hZdcRateDimuonWithCut[i]->Sumw2();
      hRFvsRate[i] = (TH1F*)hZdcRateAllWithCut[i]->Clone(Form("hRF_vs_zdc_%s",name[i]));
      hRFvsRate[i]->Divide(hZdcRateDimuonWithCut[i]);
      hRFvsRate[i]->SetMarkerStyle(21+i*4);
      hRFvsRate[i]->SetMarkerColor(i+1);
      hRFvsRate[i]->SetLineColor(i+1);
      hRFvsRate[i]->GetXaxis()->SetRangeUser(0,65);
      hRFvsRate[i]->GetYaxis()->SetRangeUser(0,50);
      TH1F *htmp = (TH1F*)hRFvsRate[i]->Clone(Form("hRF_vs_zdc_%s_clone",name[i]));
      htmp->GetXaxis()->SetRangeUser(25,65);
      htmp->GetYaxis()->SetRangeUser(5,45);
      htmp->SetMarkerStyle(21);
      htmp->SetMarkerColor(1);
      htmp->SetLineColor(1);
      TF1 *func = new TF1(Form("func_%s",name[i]),"pol0",25,65);
      htmp->Fit(func,"IR0Q");
      c = draw1D(htmp,Form("Rejection factor (production_%s);zdcRate (kHz)",name[i]));
      func->SetLineColor(2);
      func->Draw("sames");

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_SF_vs_Lumi.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_SF_vs_Lumi.png",run_type,name[i]));
	}
    }
  return;
  c = draw1D(hRFvsRate[1],Form("Rejection factor;zdcRate (kHz)"));
  hRFvsRate[0]->Draw("sames");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/SF_vs_Lumi_MidLow.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/SF_vs_Lumi_MidLow.png",run_type));
    }

}
