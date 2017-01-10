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

  Run16vs14();
}


//================================================
void Run16vs14(const int savePlot = 0)
{
  const int nData = 2;
  const char* data_name[nData] = {"Run14_AuAu200","Run16_AuAu200"};
  TFile *fdata[nData];
  fdata[0] = TFile::Open("output/Run14.AuAu200.MB.RefMult.root","read");
  fdata[1] = TFile::Open("output/Run16.AuAu200.MB.RefMult.root","read");
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
