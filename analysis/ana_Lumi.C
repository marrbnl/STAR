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

  //makeHisto();

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
  f = TFile::Open(Form("./output/Run14.AuAu200.MB.VtxEff.root"),"read");

  THnSparseF *hn = (THnSparseF*)f->Get("mhMbEvtEff");

  // Vertex distribution
  TH1F *hTpcVz = (TH1F*)hn->Projection(1);
  c = draw1D(hTpcVz,"VPD-ZDC-novtx-mon: vertex z distribution;vz_{TPC} (cm)",false,false);
  TLine *line = GetLine(-100,0,-100,0.8*hTpcVz->GetMaximum());
  line->Draw();
  line = GetLine(100,0,100,0.8*hTpcVz->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVz.pdf",run_type));

  TH1F *hVtxDz = (TH1F*)hn->Projection(2);
  c = draw1D(hVtxDz,"VPD-ZDC-novtx-mon: #Deltaz distribution for vertices;vz_{TPC}-vz_{VPD} (cm)",true,false);
  TLine *line = GetLine(-5,0,-5,0.8*hVtxDz->GetMaximum());
  line->Draw();
  line = GetLine(5,0,5,0.8*hVtxDz->GetMaximum());
  line->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_TpcVpdDz.pdf",run_type));

  // vertex cut efficiency
  TH1F *hNEvents[nCentBins][3];
  for(int i=0; i<nCentBins; i++)
    {
      hn->GetAxis(4)->SetRange(centBins_low[i],centBins_high[i]);
     
      hNEvents[i][0] = (TH1F*)hn->Projection(0);
      hNEvents[i][0]->SetName(Form("NEventsPerRun_All_cent%s",cent_Title[i]));

      hn->GetAxis(1)->SetRangeUser(-100+1e-4,100-1e-4);
      hNEvents[i][1] = (TH1F*)hn->Projection(0);
      hNEvents[i][1]->SetName(Form("NEventsPerRun_TpcVzCut_cent%s",cent_Title[i]));

      hn->GetAxis(2)->SetRangeUser(-5+1e-4,5-1e-4);
      hNEvents[i][2] = (TH1F*)hn->Projection(0);
      hNEvents[i][2]->SetName(Form("NEventsPerRun_TpcVpdDzCut_cent%s",cent_Title[i]));

      for(int j=0; j<3; j++) hNEvents[i][j]->Sumw2();

      hn->GetAxis(1)->SetRange(0,-1);
      hn->GetAxis(2)->SetRange(0,-1);
    }
  hn->GetAxis(4)->SetRange(0,-1);

  TList *list = new TList;
  TString legName[2] = {"|vz_{TPC}| < 100 cm", "+ |vz_{TPC}-vz_{VPD}| < 3 cm"};
  TH1F *hVtxEff[nCentBins][2];
  for(int i=0; i<nCentBins; i++)
    {
      hVtxEff[i][0] = (TH1F*)hNEvents[i][1]->Clone(Form("EffPerRun_TpcVzCut_cent%s",cent_Title[i]));
      hVtxEff[i][1] = (TH1F*)hNEvents[i][2]->Clone(Form("EffPerRun_TpcVpdDzCut_cent%s",cent_Title[i]));
      list->Clear();
      for(int j=0; j<2; j++)
	{
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
	      if(y==0)
		{
		  printf("[w] 0 efficiency for run %d: %1.0f/%1.0f\n",bin+start_run-1,hNEvents[i][j+1]->GetBinContent(bin),hNEvents[i][0]->GetBinContent(bin));
		}
	      if(y!=hVtxEff[i][j]->GetBinContent(bin))
		{
		  printf("[e] Mismatch eff: %4.2f != %4.2f\n",y,hVtxEff[i][j]->GetBinContent(bin));
		}
	    }
	  list->Add(hVtxEff[i][j]);
	}
      c = drawHistos(list,Form("VtxCutEff_cent%s",cent_Title[i]),Form("VPD-ZDC-novtx-mon: vertex cut efficiency (%s%%)",cent_Name[i]),kFALSE,0,30,kTRUE,-0.1,1.5,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.2,0.35,true);
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_VtxCutEff_cent%s.pdf",run_type,cent_Title[i]));
    }

  // different luminosity
  TH1F *hVtxEffLumi[nCentBins][4];
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hVtxEffLumi[j][i] = (TH1F*)hVtxEff[j][1]->Clone(Form("EffPerRun_TpcVpdDzCut_%s_cent%s",trgSetupName[i],cent_Title[j]));
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
	      int bin = hVtxEff[j][1]->FindFixBin(runnumber+0.5);
	      hVtxEffLumi[j][i]->SetBinContent(bin,hVtxEff[j][1]->GetBinContent(bin));
	      hVtxEffLumi[j][i]->SetBinError(bin,hVtxEff[j][1]->GetBinError(bin));
	    }
	}
    }
  draw1D(hVtxEffLumi[0][3]);
  TString legName2[4];
  for(int j=0; j<nCentBins; j++)
    {
      list->Clear();
      for(int i=0; i<4; i++)
	{
	  legName2[i] = Form("AuAu_200_%s_2014",trgSetupName[i]);
	  list->Add(hVtxEffLumi[j][i]);
	}
      c = drawHistos(list,Form("VtxCutEff_cent%s_%s",cent_Title[j],trgSetupName[i]),Form("VPD-ZDC-novtx-mon: vertex cut efficiency (%s%%)",cent_Name[j]),kFALSE,0,30,kTRUE,-0.1,1.5,kFALSE,kTRUE,legName2,kTRUE,"trgsetupname",0.5,0.7,0.2,0.35,true);
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_VtxCutEff_%s_cent%s.pdf",run_type,trgSetupName[i],cent_Title[j]));
    }
      
  return;

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
void makeHisto(const int savePlot = 1, const int saveHisto = 1)
{
  const char *trigName = "VPD-ZDC-novtx-mon";
  
  TH1F *hPreScale = new TH1F(Form("hPreScale_%s",trigName), Form("%s: per-scale per run;runId;Pre-scale",trigName), end_run-start_run, start_run, end_run);
  TH1F *hLiveTime = new TH1F(Form("hLiveTime_%s",trigName), Form("%s: live-time per run;runId;Live-time",trigName), end_run-start_run, start_run, end_run);
  TH1F *hNevents = new TH1F(Form("hNevents_%s",trigName), Form("%s: # of events per run;runId;Nevts",trigName), end_run-start_run, start_run, end_run);

  if(year==2013)
    {
      ifstream flumi, fnevt;
      flumi.open(Form("Rootfiles/Luminosity/%s/lum_perrun_%s.txt",run_type,trigName));
      fnevt.open(Form("Rootfiles/Luminosity/%s/nev_perday_unix_%s.txt",run_type,trigName));
      char lumi[512], nevt[512];
      int nevt_old = 0;
      int nevt_new = 0;
      int totalevents = 0;
      while(!flumi.eof() && !fnevt.eof())
	{
	  flumi.getline(lumi,512);
	  fnevt.getline(nevt,512);
	  string slumi = lumi;
	  string snevt = nevt;
	  if(slumi.length()<=0) break;
	  int counter = 0;
	  int run, starttime;
	  double prescale, livetime;
	  while(slumi.find(" ")!=string::npos)
	    {
	      string sub = slumi.substr(0,slumi.find(" "));
	      slumi.erase(0,slumi.find(" ")+1);
	      if(counter==0) run = atoi(sub.c_str());
	      if(counter==1) starttime = atoi(sub.c_str());
	      if(counter==5) prescale = atof(sub.c_str());
	      if(counter==6) livetime = atof(sub.c_str());
	      counter++;
	    }
	  int starttime2 = atoi(snevt.substr(0,12).c_str());
	  nevt_new = atoi(snevt.erase(0,12).c_str());
	  int nevents = nevt_new - nevt_old;
	  nevt_old = nevt_new;
	  totalevents += nevents;
	  if(run>=start_run && run<end_run)
	    {
	      int bin = run - start_run + 1;
	      hPreScale->SetBinContent(bin, prescale);
	      hPreScale->SetBinError(bin, 1e-10);
	      hLiveTime->SetBinContent(bin, livetime);
	      hLiveTime->SetBinError(bin, 1e-10);
	      if(starttime2==starttime)
		{
		  if(trigName=="MTD-dimuon")
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
	  printf("[i] Run %d starts at %d has %d events\n",run,starttime,nevents);
	}
      if(totalevents == nevt_new) printf("+++ Check is successful +++\n");
    }
  else if(year==2014)
    {
      ifstream flumi;
      flumi.open(Form("Rootfiles/Luminosity/%s/lum_perrun_%s.txt",run_type,trigName));
      char lumi[512], nevt[512];
      int nevt_old = 0;
      int nevt_new = 0;
      int totalevents = 0;
      while(!flumi.eof())
	{
	  flumi.getline(lumi,512);
	  string slumi = lumi;
	  if(slumi.length()<=0) break;
	  int counter = 0;
	  int run, nevents;
	  double prescale, livetime;
	  while(slumi.find(" ")!=string::npos)
	    {
	      string sub = slumi.substr(0,slumi.find(" "));
	      slumi.erase(0,slumi.find(" ")+1);
	      if(counter==0) run = atoi(sub.c_str());
	      if(counter==4) nevents  = floor(atof(sub.c_str())*6e6+0.5);
	      if(counter==5) prescale = atof(sub.c_str());
	      if(counter==6) livetime = atof(sub.c_str());
	      counter++;
	    }
	  if(run>=start_run && run<end_run)
	    {
	      int bin = run - start_run + 1;
	      hPreScale->SetBinContent(bin, prescale);
	      hPreScale->SetBinError(bin, 1e-10);
	      hLiveTime->SetBinContent(bin, livetime);
	      hLiveTime->SetBinError(bin, 1e-10);
	      hNevents->SetBinContent(bin, nevents);
	      hNevents->SetBinError(bin, 1e-10);
	    }
	  printf("[i] Run %d has %d events\n",run,nevents);
	}
    }

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

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"update");
      hPreScale->Write("",TObject::kOverwrite);
      hLiveTime->Write("",TObject::kOverwrite);
      hNevents->Write("",TObject::kOverwrite);
      fout->Close();
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
