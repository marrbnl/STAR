const int year = YEAR;

//================================================
void ana_RefMult()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //correctRefMultRun14();
  run14CentCalib();
  //Run16vs14();
}

//================================================
void correctRefMultRun14(const int savePlot = 1)
{
  const char* trigger = "mb";
  //const char* trigger = "di_mu";
  const char* lumi[2] = {"prod_mid/low","prod_high"};

  TFile *fin = NULL;
  if(trigger=="mb") fin = TFile::Open("output/Run14_AuAu200.MB.RunDepQA.root", "read");
  if(trigger=="di_mu") fin = TFile::Open("output/Run14_AuAu200.RunDepQA.root", "read");
  THnSparseF *hnRefMult = (THnSparseF*)fin->Get(Form("mhQagRefMultCorr_%s",trigger));
  
  TH1F *hRefMult[2][2];
  for(int i=0; i<2; i++)
    {
      if(i==0) hnRefMult->GetAxis(2)->SetRange(1,3);
      else     hnRefMult->GetAxis(2)->SetRange(4,4);
      for(int j=0; j<2; j++)
	{
	  hnRefMult->GetAxis(1)->SetRange(j+1, j+1);
	  hRefMult[i][j] = (TH1F*)hnRefMult->Projection(0);
	  hRefMult[i][j]->SetName(Form("hRefMult_Lumi%d_Run%d",i+1,j+1));
	  hnRefMult->GetAxis(1)->SetRange(0, -1);
	}
      hnRefMult->GetAxis(2)->SetRange(0, -1);
    }

  TString legName[2] = {"Run < 15161050", "Run > 15161050"};
  TLegend *leg;
  for(int i=0; i<2; i++)
    {    
      if(trigger=="mb")    leg = new TLegend(0.15,0.3,0.3,0.5);
      if(trigger=="di_mu") leg = new TLegend(0.15,0.7,0.3,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      for(int j=0; j<2; j++)
	{
	  if(trigger=="mb") hRefMult[i][j]->Rebin(10);
	  hRefMult[i][j]->Sumw2();
	  hRefMult[i][j]->Scale(1./hRefMult[i][j]->GetEntries());
	  hRefMult[i][j]->SetMarkerStyle(21+j*3);
	  hRefMult[i][j]->SetMarkerColor(j+1);
	  hRefMult[i][j]->SetLineColor(j+1);
	  leg->AddEntry(hRefMult[i][j], legName[j], "P");
	}
      c = draw1D(hRefMult[i][0],Form("%s: gRefMultCorr distribution (%s, %s)",run_type,trigger,lumi[i]));
      if(trigger=="mb") gPad->SetLogy();
      hRefMult[i][1]->Draw("sames");
      leg->Draw();
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Comp_gRefMult_Lumi%d_%s.pdf",run_type,i,trigger));
    }

}

//================================================
void run14CentCalib(const int savePlot = 0)
{
  // centrality calibration
  TFile *fcalib = TFile::Open("Rootfiles/run14.PicoTree.May20_ReDo_VPDNoVtx_MTD.root","read");
  TH3F *gRefMultVsZdcVsVz = (TH3F*)fcalib->Get("mgRefMultVsZdcVsVzTriggerCorrWg_0");
  TH1F *hData = (TH1F*)gRefMultVsZdcVsVz->ProjectionZ("hgRefMult_CorrWg");
  TCanvas *c = draw1D(hData,"VPD-ZDC-novtx-mon: gRefMult distribution after calibration w/ weight",true);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Guannan_VpdZdcNoVtx_gRefMultCorrWg.pdf",run_type));

  const int nCentBins = 16;
  const double xCentBins[nCentBins+1] = {11, 16, 23, 32, 44, 59, 77, 99, 126, 157, 193, 235, 283, 338, 401, 472, 700}; // SL15c
  TH1F *hCalibCent = new TH1F("hCalibCent","",nCentBins,0,nCentBins);
  for(int bin=1; bin<=nCentBins; bin++)
    {
      double err;
      double val = hData->IntegralAndError(hData->FindFixBin(xCentBins[bin-1]+0.1),
					   hData->FindFixBin(xCentBins[bin]-0.1),
					   err);
      hCalibCent->GetXaxis()->SetBinLabel(bin, Form("%d-%d%%",80-bin*5, 85-bin*5));
      hCalibCent->SetBinContent(bin, val);
      hCalibCent->SetBinError(bin, err);
    }
  hCalibCent->Scale(1./hCalibCent->Integral(15,16)*2);
  hCalibCent->GetXaxis()->SetLabelSize(0.045);
  hCalibCent->GetYaxis()->SetRangeUser(0.95,1.05);
  c = draw1D(hCalibCent,"VPD-ZDC-novtx-mon: centrality distribution after calibration w/ weight",false,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Guannan_VpdZdcNoVtx_Centrality.pdf",run_type));

  // gRefMult in prod_high
  TFile *fdata = TFile::Open("output/Run14_AuAu200.MB.VtxEff.prod_high.root","read");
  THnSparseF *mhEventMult = (THnSparseF*)fdata->Get("mhEventMult_qa_mb");
  TH1F *hgRefMultHigh = (TH1F*)mhEventMult->Projection(1);
  hgRefMultHigh->Sumw2();
  hgRefMultHigh->Scale(1./hgRefMultHigh->Integral(hgRefMultHigh->FindFixBin(100),hgRefMultHigh->FindFixBin(300)));
  hgRefMultHigh->GetXaxis()->SetRangeUser(0,700);
  hgRefMultHigh->GetYaxis()->SetRangeUser(1e-6,1e-1);
  c = draw1D(hgRefMultHigh,Form("%s: gRefMultCorr distribution w/ weight",run_type),true);
  TH1F *hgRefMultMid = (TH1F*)hData->Clone("hgRefMultMid");
  hgRefMultMid->Scale(1./hgRefMultMid->Integral(hgRefMultMid->FindFixBin(100),hgRefMultMid->FindFixBin(300)));
  hgRefMultMid->SetMarkerColor(2);
  hgRefMultMid->SetLineColor(2);
  hgRefMultMid->Draw("sames");
  for(int i=0; i<nCentBins; i++)
    {
      double max = hgRefMultHigh->GetBinContent(hgRefMultHigh->FindBin(xCentBins[i]));
      TLine *line = GetLine(xCentBins[i], 1e-6, xCentBins[i], max);
      line->Draw();
    }
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader("VPD-ZDC-novtx-mon");
  leg->AddEntry(hgRefMultMid,"prod_low/mid","L");
  leg->AddEntry(hgRefMultHigh,"prod_high","L");
  leg->Draw();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/VpdZdcNoVtx_gRefMultCorrWg_HighToMid.pdf",run_type));
}


//================================================
void Run16vs14(const int savePlot = 0)
{
  const int nData = 2;
  const char* data_name[nData] = {"Run14_AuAu200","Run16_AuAu200"};
  TFile *fdata[nData];
  fdata[0] = TFile::Open("output/Run14_AuAu200.MB.RefMult.root","read");
  fdata[1] = TFile::Open("output/Run16_AuAu200.MB.RefMult.root","read");
  THnSparseF *mhRefMult[nData];
  TH1F *hTpcVz[nData];
  TH1F *hZdcRate[nData];
  for(int i=0; i<nData; i++)
    {
      mhRefMult[i] = (THnSparseF*)fdata[i]->Get("mhRefMult");
      mhRefMult[i]->SetName(Form("mhRefMult_%s",data_name[i]));
      printf("[i] %s: %1.0f events\n",data_name[i],mhRefMult[i]->GetEntries());

      hTpcVz[i] = (TH1F*)mhRefMult[i]->Projection(0);
      hTpcVz[i]->SetName(Form("hTpcVz_%s",data_name[i]));

      hZdcRate[i] = (TH1F*)mhRefMult[i]->Projection(1);
      hZdcRate[i]->SetName(Form("hZdcRate_%s",data_name[i]));
    }

  TList *list = new TList;
  TString legName[nData];
  for(int i=0; i<nData; i++)
    {
      hTpcVz[i]->Scale(1./hTpcVz[i]->GetEntries());
      list->Add(hTpcVz[i]);
      legName[i] = data_name[i];
    }
  c = drawHistos(list,"TpcVz","Distribution of TPC vz;vz (cm)",false,0,120,true,0,0.016,false,true,legName,true,"VPDMB-novtx",0.6,0.8,0.65,0.88,false);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run16_AuAu200/ana_RefMult/Compare14vs16_TpcVz.pdf"));

  for(int i=0; i<nData; i++)
    {
      hZdcRate[i]->Scale(1./hZdcRate[i]->GetEntries());
      list->Add(hZdcRate[i]);
    }
  c = drawHistos(list,"ZdcRate","Distribution of ZDC rate;ZDCcorr (kHz)",false,0,120,true,0,0.3,false,true,legName,true,"VPDMB-novtx",0.6,0.8,0.65,0.88,false);
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run16_AuAu200/ana_RefMult/Compare14vs16_ZDCrate.pdf"));

  // gRefMult vs. TPCvz & ZDCrate
  const int nZdcRate = 3;
  const double zdc_cut[nZdcRate+1] = {40,45,50,55};
  const int nTpcVz = 5;
  const double vz_cut[nTpcVz+1] = {-100,-30,-10,10,30,100};
  TH1F *hgRefMult[nZdcRate][nTpcVz][nData];
  for(int j=0; j<nZdcRate; j++)
    {
      c = new TCanvas(Form("gRefMult_Zdc%1.0f-%1.0f",zdc_cut[j], zdc_cut[j+1]),Form("gRefMult_Zdc%1.0f-%1.0f",zdc_cut[j], zdc_cut[j+1]),1100,700);
      c->Divide(3, 2);
      for(int k=0; k<nTpcVz; k++)
	{
	  for(int i=0; i<nData; i++)
	    {
	      mhRefMult[i]->GetAxis(0)->SetRangeUser(vz_cut[k]+1e-4, vz_cut[k+1]-1e-4);
	      mhRefMult[i]->GetAxis(1)->SetRangeUser(zdc_cut[j]+1e-4, zdc_cut[j+1]-1e-4);
	      hgRefMult[j][k][i] = (TH1F*)mhRefMult[i]->Projection(3);
	      hgRefMult[j][k][i]->SetName(Form("hgRefMult_%d_%d_%d",i,j,k));
	      hgRefMult[j][k][i]->Sumw2();
	      hgRefMult[j][k][i]->Rebin(5);
	      hgRefMult[j][k][i]->Scale(1./hgRefMult[j][k][i]->GetEntries());
	      hgRefMult[j][k][i]->SetLineColor(i+1);
	      hgRefMult[j][k][i]->SetMarkerStyle(21+i*4);
	      hgRefMult[j][k][i]->SetMarkerColor(i+1);
	      hgRefMult[j][k][i]->GetXaxis()->SetRangeUser(0,700);
	      hgRefMult[j][k][i]->SetTitle(";gRefMult");
	      mhRefMult[i]->GetAxis(0)->SetRange(0,-1);
	      mhRefMult[i]->GetAxis(1)->SetRange(0,-1);
	    }
	  c->cd(k+1);
	  gPad->SetLogy();
	  hgRefMult[j][k][0]->Draw("PE");
	  hgRefMult[j][k][1]->Draw("samesPE");
	  TPaveText *t1 = GetTitleText(Form("%1.0f < vz_{TPC} < %1.0f cm",vz_cut[k], vz_cut[k+1]),0.055);
	  t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.15,0.2,0.35,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader("VPDMB-novtx");
      for(int i=0; i<nData; i++)
	{
	  leg->AddEntry(hgRefMult[j][0][i],data_name[i],"P");
	}
      leg->Draw();
      c->cd(6);
      TPaveText *t1 = GetPaveText(0.2,0.8,0.5,0.7,0.07,62);
      t1->AddText(Form("%1.0f < ZDCcorr < %1.0f kHz",zdc_cut[j], zdc_cut[j+1]));
      t1->Draw(); 
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run16_AuAu200/ana_RefMult/Compare14vs16_gRefMult_ZDC%1.0f_%1.0f.pdf",zdc_cut[j], zdc_cut[j+1]));
    }
}
