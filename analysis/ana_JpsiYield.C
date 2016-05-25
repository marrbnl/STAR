TFile *f;
const int year = YEAR;
TString run_cfg_name;

//================================================
void ana_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2013)
    {
      f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);

  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));

  //anaYield();
  yieldVsLumi();
  //pt2scan();
}

//================================================
void yieldVsLumi(int savePlot = 1, int saveHisto = 1)
{
  const int gApplyWeight = 1;
  const int nSetup = 4;
  const char *setup_name[nSetup] = {"production_2014","production_low_2014","production_mid_2014","production_high_2014"};
  TList *list = new TList;
  TString legName[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      legName[i] = setup_name[i];
    }

  // luminosity
  TH2F *hTrgSetupVsZdcRate = (TH2F*)f->Get("mhTrgSetupVsZdcRate_di_mu");
  TH1F *hZdcRate[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      hZdcRate[i] = (TH1F*)hTrgSetupVsZdcRate->ProjectionX(Form("hZdcRate_%d",i),i+1,i+1);
    }
  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hZdcRate[i]);
    }
  c = drawHistos(list,"ZdcRate","Distribution of ZDC rate",true,0,120,true,1,5e7,false,true,legName,true,run_type,0.5,0.8,0.62,0.88,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sZdcRate_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  fileName.ReplaceAll(".root",".yield.root");
  TFile *fin;
  if(saveHisto) fin = TFile::Open(fileName,"update");
  else          fin = TFile::Open(fileName,"read");

  TH1F *hJpsiYield[nSetup][nCentBins];
  for(int i=0; i<nSetup; i++)
    {
      for(int j=0; j<nCentBins; j++)
	{
	  hJpsiYield[i][j] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_cent%s%s%s",cent_Title[j],gWeightName[gApplyWeight],gTrgSetupName[i+1]));
	}
    }

  // # of events in each setup/centrality
  TH1F *hCent = (TH1F*)f->Get("mhCentrality_di_mu");
  double centFraction[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      centFraction[i] = hCent->Integral(centBins_low[i],centBins_high[i])/hCent->Integral(0,-1);
    }


  // # of Jpsi signal per event
  double nEvents[nSetup][nCentBins];
  TH1F *hEvtRunAcc = (TH1F*)f->Get("mhEvtRunAcc_di_mu");
  for(int i=0; i<nSetup; i++)
    {
      for(int j=0; j<nCentBins; j++) nEvents[i][j] = 0;
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_%s.list",run_type,setup_name[i]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  if(runnumber>0)
	    {
	      for(int j=0; j<nCentBins; j++)
		{
		  nEvents[i][j] += hEvtRunAcc->GetBinContent(hEvtRunAcc->FindFixBin(runnumber)) * centFraction[j];
		}
	    }
	}
    }

  // # of Jpsi in each setup/centrality
  TH1F *hJpsiInCent[nSetup];
  TH1F *hJpsiPerEvent[nSetup];
  for(int i=0; i<nSetup; i++)
    {
      hJpsiPerEvent[i] = new TH1F(Form("hJpsiPerEvent%s%s",gWeightName[gApplyWeight],gTrgSetupName[i+1]),"Jpsi counts per event in each centrality bin (p_{T} > 0)",nCentBins-1,0,nCentBins-1);
      hJpsiInCent[i] = new TH1F(Form("NJpsiInCent%s%s",gWeightName[gApplyWeight],gTrgSetupName[i+1]),"Jpsi counts in each centrality bin (p_{T} > 0)",nCentBins-1,0,nCentBins-1);

      for(int j=0; j<nCentBins; j++)
	{
	  printf("# of events in %s %s: %1.0f\n",gTrgSetupName[i+1],cent_Title[j],nEvents[i][j]);
	  hJpsiInCent[i]->SetBinContent(j, hJpsiYield[i][j]->GetBinContent(0));
	  hJpsiInCent[i]->SetBinError(j, hJpsiYield[i][j]->GetBinError(0));
	  if(j>0) hJpsiInCent[i]->GetXaxis()->SetBinLabel(j,Form("%s%%",cent_Name[j]));

	  hJpsiPerEvent[i]->SetBinContent(j, hJpsiYield[i][j]->GetBinContent(0)/nEvents[i][j]);
	  hJpsiPerEvent[i]->SetBinError(j, hJpsiYield[i][j]->GetBinError(0)/nEvents[i][j]);
	  if(j>0) hJpsiPerEvent[i]->GetXaxis()->SetBinLabel(j,Form("%s%%",cent_Name[j]));
	}
      hJpsiInCent[i]->GetXaxis()->SetLabelSize(0.05);
      hJpsiPerEvent[i]->GetXaxis()->SetLabelSize(0.05);
    }

  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hJpsiInCent[i]);
    }
  c = drawHistos(list,"Jpsi_Counts","",kFALSE,0,30,true,0,1e4,kFALSE,true,legName,true,run_type,0.45,0.75,0.62,0.88,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiCounts_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  list->Clear();
  for(int i=0; i<nSetup; i++)
    {
      list->Add(hJpsiPerEvent[i]);
    }
  c = drawHistos(list,"Jpsi_Yield","",kFALSE,0,30,true,1e-5,5e-4,true,true,legName,true,run_type,0.15,0.3,0.62,0.88,kTRUE);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsiYield_TrgSetup.pdf",run_type,run_cfg_name.Data()));

  if(saveHisto)
    {
      fin->cd();
      for(int i=0; i<nSetup; i++)
	{
	  hJpsiInCent[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void anaYield(int savePlot = 1)
{
  TString fileName = f->GetName();
  fileName.ReplaceAll("output","Rootfiles");
  fileName.ReplaceAll(".root",Form(".pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut));
  fileName.ReplaceAll(".root",".yield.root");

  TFile *fin = TFile::Open(fileName,"read");

  bool drawLegend = true;
  if(nCentBins==1) drawLegend = false;

  TList *list = new TList;
  TString legName[nCentBins];
  TH1F *hFitYield[nCentBins];
  TH1F *hBinCountYield[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hFitYield[i] = (TH1F*)fin->Get(Form("Jpsi_FitYield_cent%s",cent_Title[i]));
      hBinCountYield[i] = (TH1F*)fin->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[i]));
      double count_1, error_1, count_2, error_2;
      count_1 = hFitYield[i]->IntegralAndError(1,nPtBins-1,error_1);
      count_2 = hBinCountYield[i]->IntegralAndError(1,nPtBins-1,error_2);
      printf("Total # of J/psi in %s%%: fit - %1.0f (%2.2f#sigma), bin count - %1.0f (%2.2f#sigma)\n",cent_Name[i],count_1,count_1/error_1,count_2,count_2/error_2);
      scaleHisto(hFitYield[i],1,1,kTRUE,kFALSE,kFALSE);
      list->Add(hFitYield[i]);
      legName[i] = Form("%s%%",cent_Name[i]);
    }
  if(year==2013) c = drawHistos(list,"Jpsi_FitYield","J/psi yield from bin counting;p_{T} (GeV/c);Counts/(1 GeV/c)",kFALSE,0,30,kTRUE,1e-6,1.2*hBinCountYield[0]->GetMaximum(),kFALSE,drawLegend,legName,drawLegend,"Centrality",0.45,0.75,0.55,0.8,kTRUE);
  else c = drawHistos(list,"Jpsi_FitYield","Raw J/psi yield in each p_{T} bin;p_{T} (GeV/c);Counts/(1 GeV/c)",kFALSE,0,30,kTRUE,1,5e4,true,kTRUE,legName,kTRUE,run_type,0.45,0.75,0.65,0.88,kTRUE,0.05);
  
  for(int i=0; i<nCentBins; i++)
    {
      scaleHisto(hBinCountYield[i],1,1,kTRUE,kFALSE,kFALSE);
      hBinCountYield[i]->SetMarkerStyle(24);
      hBinCountYield[i]->SetMarkerColor(color[i]);
      hBinCountYield[i]->SetLineColor(color[i]);
      TGraphErrors *gr = new TGraphErrors(hBinCountYield[i]);
      offset_x(gr,0.2);
      gr->Draw("samesPEZ");
    }
  leg = new TLegend(0.65,0.7,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hFitYield[0],"Fitting","P");
  leg->AddEntry(hBinCountYield[0],"Bin counting","P");
  leg->Draw();
  

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Yield_CentBins.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Yield_CentBins.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hSignif[nCentBins][2];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hSignif[i][j] = (TH1F*)hFitYield[0]->Clone(Form("hSignif_%s_%d",cent_Name[i],j));
	  hSignif[i][j]->Reset();
	}
      for(int bin=1; bin<=hSignif[i][0]->GetNbinsX(); bin++)
	{
	  hSignif[i][0]->SetBinContent(bin,hFitYield[i]->GetBinContent(bin)/hFitYield[i]->GetBinError(bin));
	  hSignif[i][1]->SetBinContent(bin,hBinCountYield[i]->GetBinContent(bin)/hBinCountYield[i]->GetBinError(bin));
	}
    }

  const char *name[2] = {"Fitting","BinCounting"};
  for(int j=0; j<2; j++)
    {
      list->Clear();
      for(int i=0; i<nCentBins; i++)
	{
	  list->Add(hSignif[i][j]);
	}
      c = drawHistos(list,Form("%s_Signif",name[j]),Form("%s: J/psi significance;p_{T} (GeV/c);significance",name[j]),kFALSE,0,20,kTRUE,0.1,20,kFALSE,drawLegend,legName,drawLegend,run_type,0.55,0.7,0.6,0.85,kTRUE);
      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Significance_%s.pdf",run_type,run_cfg_name.Data(),name[j]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiYield/%sJpsi_Significance_%s.png",run_type,run_cfg_name.Data(),name[j]));
	}
    }
}

//================================================
void pt2scan(int savePlot = 0)
{
  const int ncuts = 2;
  double pt2_cuts[ncuts] = {1.0,1.2};
  TFile *fin[ncuts];
  TH1F *hJpsiYield[ncuts];
}
