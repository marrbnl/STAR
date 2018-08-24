//================================================
void check_Tracking()
{
  gStyle->SetOptStat(0);

  trackQA();
  //tightCutsEmb();
  //trackYield();
}

//================================================
void trackYield(const int savePlot = 0)
{
  const char* trigTitle = "vpd30";
  TFile *fin = 0x0;
  char *trigName = "";
  char *libVer = "";

  if(trigTitle=="vpd30")
    {
      fin = TFile::Open("output/Run14_AuAu200.MB.P17id.Study.TpcTracking.root", "read");
      trigName = "mb";
      libVer = "P17id";
    }
  else if(trigTitle=="dimuon")
    {
      fin = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
      trigName = "di_mu";
      libVer = "P15ie";
    }

  // vertex and centrality
  THnSparseF *hnCent = (THnSparse*)fin->Get(Form("mhnCent_%s",trigName));
  hnCent->GetAxis(1)->SetRange(1,16);
  const int nHisto = 4;
  const char* histoName[nHisto] = {"gRefMultCorr", "cent", "ZDC", "TpcVz"};
  TH1F *histo[nHisto][3];
  TH1F *hCent[3];
  for(int k=0; k<3; k++)
    {
      hnCent->GetAxis(4)->SetRange(k+2, k+2);
      for(int i=0; i<nHisto; i++)
	{
	  histo[i][k] = (TH1F*)hnCent->Projection(i);
	  //histo[i][k]->Sumw2();
	  histo[i][k]->SetName(Form("%s%s_0080",histoName[i], gTrgSetupTitle[k+2]));
	  if(i==1)
	    {
	      hCent[k] = (TH1F*)histo[i][k]->Clone(Form("%s_clone",histo[i][k]->GetName()));
	    }
	  if(i!=2) histo[i][k]->Scale(1./histo[i][k]->Integral());
	}
    }
  hnCent->GetAxis(4)->SetRange(0, -1);

  for(int i=0; i<nHisto; i++)
    {
      TCanvas *c = new TCanvas(Form("c%s",histoName[i]), Form("c%s",histoName[i]), 800, 600);
      TLegend *leg =  new TLegend(0.65, 0.68, 0.75, 0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.035);
      for(int k=0; k<3; k++)
	{
	  histo[i][k]->SetMarkerStyle(21+k);
	  histo[i][k]->SetMarkerColor(color[k+1]);
	  histo[i][k]->SetLineColor(color[k+1]);
	  histo[i][k]->SetMarkerSize(1.2);
	  histo[i][k]->SetTitle("");
	  if(i==0)
	    {
	      gPad->SetLogy();
	      histo[i][k]->GetXaxis()->SetRangeUser(0, 700);
	      histo[i][k]->GetYaxis()->SetRangeUser(1e-7, 1e-1);
	      histo[i][k]->SetXTitle("gRefMultCorr");
	    }
	  else if(i==1)
	    {
	      histo[i][k]->GetYaxis()->SetRangeUser(0.045, 0.08);
	    }
	  else if(i==2)
	    {
	      histo[i][k]->GetYaxis()->SetRangeUser(0, 7e5);
	    }
	  if(k==0) histo[i][k]->Draw();
	  else     histo[i][k]->Draw("sames");
	  leg->AddEntry(histo[i][k], Form("%s",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s_%s (%s, 0-80%%)",run_type,trigTitle,libVer), 0.04);
      t1->Draw();
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_%sInLumi.pdf",run_type,trigName,histoName[i]));
    }

  // track yield 
  const int nCentBins = 4; 
  int centBins_low[nCentBins]  = {15,11,7,1};
  int centBins_high[nCentBins] = {16,12,8,4};
  TString cent_Name[nCentBins];
  TString cent_Title[nCentBins];

  double nEvents[nCentBins][3];
  TH1F *hTrkPt[nCentBins][3];
  THnSparseF *hnTrk = (THnSparseF*)fin->Get(Form("mhnTrkPt_%s",trigName));
  for(int j=0; j<nCentBins; j++)
    {
      cent_Name[j] = Form("%d%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
      cent_Title[j] = Form("%d-%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);

      hnTrk->GetAxis(3)->SetRange(centBins_low[j], centBins_high[j]);
      for(int k=0; k<3; k++)
	{
	  hnTrk->GetAxis(4)->SetRange(k+2, k+2);
	  nEvents[j][k] = hCent[k]->Integral(centBins_low[j], centBins_high[j]);
	  printf("[i] %s, %s%%, N = %1.0f\n",gTrgSetupTitle[k+2], cent_Title[j].Data(), nEvents[j][k]);

	  hTrkPt[j][k] = (TH1F*)hnTrk->Projection(0);
	  hTrkPt[j][k]->Scale(1./nEvents[j][k]);
	  hTrkPt[j][k]->SetName(Form("hTrkPt_%d%d%s",cent_Name[j].Data(),gTrgSetupTitle[k+2]));
	}
      hnTrk->GetAxis(4)->SetRange(0, -1);
    }
  hnTrk->GetAxis(3)->SetRange(0, -1);

  // raw distribution
  TCanvas *c = new TCanvas(Form("cTrkPtRaw_%s",trigName), Form("cTrkPtRaw_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.6, 0.6, 0.8, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      gPad->SetLogy();
      for(int k=0; k<3; k++)
	{
	  hTrkPt[j][k]->SetMarkerStyle(21+k);
	  hTrkPt[j][k]->SetMarkerColor(color[k+1]);
	  hTrkPt[j][k]->SetLineColor(color[k+1]);
	  hTrkPt[j][k]->SetMarkerSize(1.2);
	  hTrkPt[j][k]->SetTitle("");
	  if(k==0) hTrkPt[j][k]->Draw();
	  else     hTrkPt[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkPt[j][k], Form("%s",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s, %s (%s%%)",run_type,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkPtRaw.pdf",run_type,trigName));

  TCanvas *c = new TCanvas(Form("cTrkPtRawRatio_%s",trigName), Form("cTrkPtRawRatio_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2, 0.15, 0.4, 0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=1; k<3; k++)
	{
	  TH1F *hRatio = (TH1F*)hTrkPt[j][k]->Clone(Form("%s_ratio", hTrkPt[j][k]->GetName()));
	  hRatio->Divide(hTrkPt[j][0]);
	  hRatio->GetYaxis()->SetRangeUser(0.3, 1.4);
	  if(k==0) hRatio->Draw();
	  else     hRatio->Draw("sames");
	  if(j==0) leg->AddEntry(hRatio, Form("%s/prod_low",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s, %s (%s%%)",run_type,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkPtRawRatio.pdf",run_type,trigName));

  // tracking efficiency from embedding
  TFile *fEmb = 0x0;
  if(trigTitle=="vpd30")
    {
      fEmb = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
    }
  else if(trigTitle=="dimuon")
    {
      fEmb = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
    }
  TH1F *hTrkEff[nCentBins][3];
  THnSparseF *hnEmbTrk[2];
  hnEmbTrk[0] = (THnSparseF*)fEmb->Get(Form("mhMcTrkPtEff_Tpc_%s",trigName));
  hnEmbTrk[1] = (THnSparseF*)fEmb->Get(Form("mhMcTrkPtEff_MC_%s",trigName));
  for(int i=0; i<2; i++)
    {
      hnEmbTrk[i]->GetAxis(1)->SetRangeUser(-0.5, 0.5);
      for(int j=0; j<nCentBins; j++)
	{
	  hnEmbTrk[i]->GetAxis(2)->SetRange(centBins_low[j], centBins_high[j]);
	  for(int k=0; k<3; k++)
	    {
	      hnEmbTrk[i]->GetAxis(4)->SetRange(k+2, k+2);
	      TH1F *h1tmp = (TH1F*)hnEmbTrk[i]->Projection(0);
	      h1tmp->SetName(Form("hEmbTrkPt_%s%s",cent_Name[j].Data(),gTrgSetupTitle[k+2]));
	      h1tmp->Sumw2();
	      h1tmp->Rebin(10);
	      if(i==0) hTrkEff[j][k] = (TH1F*)h1tmp->Clone(Form("hEmbTrkEff_%s%s",cent_Name[j].Data(),gTrgSetupTitle[k+2]));
	      else     hTrkEff[j][k]->Divide(h1tmp);
	    }
	  hnEmbTrk[i]->GetAxis(4)->SetRange(0, -1);
	}
      hnEmbTrk[i]->GetAxis(2)->SetRange(0, -1);
      hnEmbTrk[i]->GetAxis(1)->SetRange(0, -1);
    }

  TCanvas *c = new TCanvas(Form("cTrkEff_%s",trigName), Form("cTrkEff_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.6, 0.25, 0.8, 0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=0; k<3; k++)
	{
	  hTrkEff[j][k]->SetMarkerStyle(21+k);
	  hTrkEff[j][k]->SetMarkerColor(color[k+1]);
	  hTrkEff[j][k]->SetLineColor(color[k+1]);
	  hTrkEff[j][k]->SetMarkerSize(1.2);
	  hTrkEff[j][k]->GetXaxis()->SetRangeUser(0,5);
	  hTrkEff[j][k]->GetYaxis()->SetRangeUser(0,1);
	  hTrkEff[j][k]->SetTitle(";p_{T} (GeV/c);Efficiency");
	  if(k==0) hTrkEff[j][k]->Draw();
	  else     hTrkEff[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkEff[j][k], Form("%s",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s embedding, %s (%s%%)",run_type,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkEffVsPt.pdf",run_type,trigName));

  // corrected yield
  TCanvas *c = new TCanvas(Form("cTrkPtCorr_%s",trigName), Form("cTrkPtCorr_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.6, 0.6, 0.8, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  TH1F *hTrkPtCorr[nCentBins][3];
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      gPad->SetLogy();
      for(int k=0; k<3; k++)
	{
	  hTrkPtCorr[j][k] = (TH1F*)hTrkPt[j][k]->Clone(Form("hTrkPtCorr_%d%d%s",cent_Name[j].Data(),gTrgSetupTitle[k+2]));
	  for(int bin=1; bin<=hTrkPtCorr[j][k]->GetNbinsX(); bin++)
	    {
	      double eff = hTrkEff[j][k]->GetBinContent(bin);
	      if(eff<=0) continue;
	      hTrkPtCorr[j][k]->SetBinContent(bin, hTrkPtCorr[j][k]->GetBinContent(bin)/eff);
	      hTrkPtCorr[j][k]->SetBinError(bin, hTrkPtCorr[j][k]->GetBinError(bin)/eff);
	    }
	  if(trigTitle=="vpd30") hTrkPtCorr[j][k]->GetXaxis()->SetRangeUser(0,5);
	  hTrkPtCorr[j][k]->SetTitle(";p_{T} (GeV/c);1/N dN/dp_{T}");
	  if(k==0) hTrkPtCorr[j][k]->Draw();
	  else     hTrkPtCorr[j][k]->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkPtCorr[j][k], Form("%s",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s, %s (%s%%)",run_type,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkPtCorr.pdf",run_type,trigName));

  TCanvas *c = new TCanvas(Form("cTrkPtCorrRatio_%s",trigName), Form("cTrkPtCorrRatio_%s",trigName), 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2, 0.15, 0.4, 0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int j=0; j<nCentBins; j++)
    {
      c->cd(j+1);
      for(int k=1; k<3; k++)
	{
	  TH1F *hRatio = (TH1F*)hTrkPtCorr[j][k]->Clone(Form("%s_ratio", hTrkPtCorr[j][k]->GetName()));
	  hRatio->Divide(hTrkPtCorr[j][0]);
	  hRatio->GetYaxis()->SetRangeUser(0.5, 1.2);
	  hRatio->SetYTitle("Ratio");
	  if(k==0) hRatio->Draw();
	  else     hRatio->Draw("sames");
	  if(j==0) leg->AddEntry(hRatio, Form("%s/prod_low",gTrgSetupTitle[k+2]), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%s, %s (%s%%)",run_type,trigTitle,cent_Title[j].Data()), 0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TrkPtCorrRatio.pdf",run_type,trigName));
}

//================================================
void trackQA(const int savePlot = 0)
{
  const int nCentBins       = 8; 
  TString cent_Name[nCentBins];
  TString cent_Title[nCentBins];
  const int centBins_low[nCentBins]  = {1,5,7,9,11,13,15,16};
  const int centBins_high[nCentBins] = {4,6,8,10,12,14,15,16}; 
  for(int j=0; j<nCentBins; j++)
    {
      cent_Name[j]  = Form("%d-%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
      cent_Title[j] = Form("%d%d",(16-centBins_high[j])*5,(17-centBins_low[j])*5);
    }
  const int nType = 2;
  const char* typeName[nType] = {"NHitsFit", "Dca"};
  const double min_pt = 2;
  const double max_pt = 10;

  //----------------------------------------------
  const int nData = 3;
  const char* dataTitle[3] = {"Data_dimu", "Data_mb", "Embed_dimu"};
  const char* histoName[3] = {"di_mu", "mb", "di_mu"};
  const char* saveName = histoName[0];
  TString legNameData[3] = {"Data dimuon: #mu cand.", "Data mb: charged h", "Embed: #mu"};
  TFile *fin[3];
  fin[0] = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.MB.Study.TpcTracking.root", "read");
  fin[2] = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.Dca3cm.root", "read");

  /*
  const int nData = 2;
  const char* dataTitle[3] = {"Data", "Embed_MC", "Embed_data"};
  const char* histoName = "mb";
  const char* legTitle = "vpd30";
  TString legNameData[3] = {"Data: hadron", "Embed: MC #pi", "Embed: data #pi"};
  TFile *fin[3];
  fin[0] = TFile::Open("output/Run14_AuAu200.MB.P17id.Study.TpcTracking.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
  fin[2] = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
  */

  THnSparseF *hnQa[nData][nType];
  TH1F *hTrkDis[nData][nType][nCentBins][gNTrgSetup-1];
  THnSparseF* hnTmp;
  for(int d=0; d<nData; d++)
    {
      for(int i=0; i<nType; i++)
	{
	  if(d==0)      hnQa[d][i] = (THnSparseF*)fin[d]->Get(Form("mhnTrkPt_%s",histoName[d]));
	  else if(d==1) hnQa[d][i] = (THnSparseF*)fin[d]->Get(Form("mhnTrkPt_%s",histoName[d]));
	  else if(d==2) hnQa[d][i] = (THnSparseF*)fin[d]->Get(Form("hTrk%s_MCreco_%s",typeName[i],histoName[d]));
	  hnQa[d][i]->SetName(Form("%s_%s",dataTitle[d], typeName[i]));

	  hnQa[d][i]->GetAxis(0)->SetRangeUser(min_pt, max_pt);
	  for(int j=0; j<nCentBins; j++)
	    {
	      hnQa[d][i]->GetAxis(3)->SetRange(centBins_low[j], centBins_high[j]);
	      for(int k=0; k<gNTrgSetup-1; k++)
		{
		  hnQa[d][i]->GetAxis(4)->SetRange(k+1, k+1);
		  if(d==0) hTrkDis[d][i][j][k] = (TH1F*)hnQa[d][i]->Projection(i+1);
		  else if(d==1) hTrkDis[d][i][j][k] = (TH1F*)hnQa[d][i]->Projection(i+1);
		  else if(d==2) hTrkDis[d][i][j][k] = (TH1F*)hnQa[d][i]->Projection(1);
		  hTrkDis[d][i][j][k]->SetName(Form("%s_Trk%s_cent%s%s",dataTitle[d], typeName[i],cent_Title[j].Data(), gTrgSetupTitle[k+1]));
		  if(d>1) hTrkDis[d][i][j][k]->Sumw2();
		  if(hTrkDis[d][i][j][k]->Integral()>0)
		    hTrkDis[d][i][j][k]->Scale(1./hTrkDis[d][i][j][k]->Integral());
		}
	      hnQa[d][i]->GetAxis(4)->SetRange(0, -1);
	    }
	  hnQa[d][i]->GetAxis(3)->SetRange(0, -1);
	  hnQa[d][i]->GetAxis(0)->SetRange(0, -1);
	}
    }

  // different data sets
  TCanvas *cData[nType][nCentBins];
  TH1F *hCutAll[nType][nData][gNTrgSetup-1];
  TH1F *hCutAcc[nType][nData][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      for(int d=0; d<nData; d++)
	{
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hCutAll[i][d][k] = new TH1F(Form("h%s_All_%s%s",typeName[i], dataTitle[d], gTrgSetupTitle[k+1]), "", nCentBins, 0, nCentBins);
	      hCutAcc[i][d][k] = new TH1F(Form("h%s_Acc_%s%s",typeName[i], dataTitle[d], gTrgSetupTitle[k+1]), "", nCentBins, 0, nCentBins);
	    }
	}
      for(int j=0; j<nCentBins; j++)
	{
	  cData[i][j] = new TCanvas(Form("cData_%s_%s",typeName[i], cent_Name[j].Data()), Form("cData_%s_%s",typeName[i], cent_Name[j].Data()), 1100, 700);
	  cData[i][j]->Divide(2, 2);
	  TLegend *leg =  0x0;
	  if(i==0) leg = new TLegend(0.15, 0.6, 0.3, 0.85);
	  if(i==1) leg = new TLegend(0.5, 0.6, 0.7, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      cData[i][j]->cd(k+1);
	      if(i==1) gPad->SetLogy();
	      double all, all_err, acc, acc_err;
	      for(int d=0; d<nData; d++)
		{
		  TH1F *hplot = (TH1F*)hTrkDis[d][i][j][k]->Clone(Form("%s_clone",hTrkDis[d][i][j][k]->GetName()));
		  hplot->SetMarkerStyle(20+2*d);
		  hplot->SetMarkerColor(color[d]);
		  hplot->SetLineColor(color[d]);
		  hplot->SetMarkerSize(1.2);
		  hplot->SetTitle("");
		  if(i==0) 
		    {
		      hplot->GetYaxis()->SetRangeUser(0,0.12);
		      all = hTrkDis[d][i][j][k]->IntegralAndError(hTrkDis[d][i][j][k]->FindFixBin(10),
								  hTrkDis[d][i][j][k]->GetNbinsX(),
								  all_err);
		      acc = hTrkDis[d][i][j][k]->IntegralAndError(hTrkDis[d][i][j][k]->FindFixBin(20),
								  hTrkDis[d][i][j][k]->GetNbinsX(),
								  acc_err);
		    }
		  if(i==1) 
		    {
		      hplot->GetXaxis()->SetRangeUser(0,3.5);
		      hplot->GetYaxis()->SetRangeUser(5e-4, 1);
		      all = hTrkDis[d][i][j][k]->IntegralAndError(1, hTrkDis[d][i][j][k]->FindFixBin(3), all_err);
		      acc = hTrkDis[d][i][j][k]->IntegralAndError(1, hTrkDis[d][i][j][k]->FindFixBin(1), acc_err);
		    }
		  hCutAll[i][d][k]->SetBinContent(j+1, all);
		  hCutAll[i][d][k]->SetBinError(j+1, all_err);
		  hCutAcc[i][d][k]->SetBinContent(j+1, acc);
		  hCutAcc[i][d][k]->SetBinError(j+1, acc_err);

		  if(d==0) hplot->DrawCopy("P");
		  else     hplot->DrawCopy("samesP");
		  if(k==0) leg->AddEntry(hplot, legNameData[d].Data(), "P");
		}
	      TPaveText *t1 = GetTitleText(Form("%s%s: %s%%",run_type,gTrgSetupTitle[k+1],cent_Name[j].Data()), 0.055);
	      t1->Draw();
	    }
	  cData[i][j]->cd(1);
	  leg->Draw();
	  if(savePlot) cData[i][j]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sVsData_%s_%s.pdf",run_type,typeName[i],cent_Title[j].Data(),saveName));
	}
    }

  TCanvas *cEffVsCent[nData];
  TCanvas *cEffVsData[nData];
  TGraphAsymmErrors *gCutEff[nType][nCentBins-1][nData];
  for(int i=0; i<nType; i++)
    {
      cEffVsCent[i] = new TCanvas(Form("cEffVsCent_%s",typeName[i]), Form("cEffVsCent_%s",typeName[i]), 1100, 700);
      cEffVsCent[i]->Divide(2,2);
      TLegend *leg = new TLegend(0.15, 0.15, 0.3, 0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cEffVsCent[i]->cd(k+1);
	  for(int d=0; d<nData; d++)
	    {
	      gCutEff[i][d][k] = new TGraphAsymmErrors(hCutAcc[i][d][k], hCutAll[i][d][k],"cl=0.683 b(1,1) mode");
	      if(d==0)
		{
		  TH1F *hplot = (TH1F*)hCutAcc[i][d][k]->Clone(Form("plot_%s",hCutAcc[i][d][k]->GetName()));
		  for(int bin=1; bin<=hplot->GetNbinsX(); bin++)
		    {
		      hplot->GetXaxis()->SetBinLabel(bin, Form("%s%%",cent_Name[bin-1].Data()));
		    }
		  hplot->GetXaxis()->SetLabelSize(0.06);
		  hplot->Reset();
		  if(i==0) hplot->GetYaxis()->SetRangeUser(0.85, 1.0);
		  if(i==1) hplot->GetYaxis()->SetRangeUser(0.7, 1.0);
		  hplot->Draw();
		}
	      gCutEff[i][d][k]->SetMarkerStyle(20+2*d);
	      gCutEff[i][d][k]->SetMarkerColor(color[d]);
	      gCutEff[i][d][k]->SetLineColor(color[d]);
	      gCutEff[i][d][k]->SetMarkerSize(1.2);
	      gCutEff[i][d][k]->SetTitle("");
	      gCutEff[i][d][k]->Draw("samesPE");
	      if(j==0) leg->AddEntry(gCutEff[i][d][k], legNameData[d].Data(), "P");
	    }
	  if(i==0)
	    {
	      TPaveText *t1 = GetTitleText(Form("eff = (NHits>=20)/(NHits>=10) (%s)",gTrgSetupTitle[k+1]),0.055);
	      t1->Draw();
	    }
	  if(i==1)
	    {
	      TPaveText *t1 = GetTitleText(Form("eff = (DCA<=3)/(DCA<=1) (%s)",gTrgSetupTitle[k+1]),0.055);
	      t1->Draw();
	    }
	}
      cEffVsCent[i]->cd(3);
      leg->Draw();
      if(savePlot) cEffVsCent[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sEffComp_%s.pdf",run_type,typeName[i],saveName));
    }
  return;

  // in luminosity bins
  TCanvas *cCent[nData][nType];
  for(int d=0; d<2; d++)
    {
      for(int i=0; i<nType; i++)
	{
	  cCent[d][i] = new TCanvas(Form("cCent_%s_%s",dataTitle[d], typeName[i]), Form("cCent_%s_%s",dataTitle[d], typeName[i]), 1100, 700);
	  cCent[d][i]->Divide(2, 2);
	  TLegend *leg =  0x0;
	  if(i==0) leg = new TLegend(0.2, 0.45, 0.4, 0.8);
	  if(i==1) leg = new TLegend(0.6, 0.45, 0.8, 0.8);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(Form("%1.0f < p_{T} < %1.0f GeV/c, |#eta| < 0.5",min_pt,max_pt));
	  for(int j=0; j<nCentBins-1; j++)
	    {
	      if(j%2==1) continue;
	      cCent[d][i]->cd(j/2+1);
	      if(i==1) gPad->SetLogy();
	      for(int k=0; k<gNTrgSetup-1; k++)
		{
		  TH1F *hplot = (TH1F*)hTrkDis[d][i][j][k]->Clone(Form("%s_clone",hTrkDis[d][i][j][k]->GetName()));
		  hplot->SetMarkerStyle(20+k);
		  hplot->SetMarkerColor(color[k]);
		  hplot->SetLineColor(color[k]);
		  hplot->SetMarkerSize(1.2);
		  hplot->SetTitle("");
		  if(i==0) 
		    {
		      hplot->GetYaxis()->SetRangeUser(0,0.1);
		    }
		  if(i==1) 
		    {
		      hplot->GetXaxis()->SetRangeUser(0,3.5);
		      if(d==1) hplot->GetYaxis()->SetRangeUser(1e-4, 1);
		      if(d==0) hplot->GetYaxis()->SetRangeUser(1e-3, 1);
		    }
		  if(k==0) hplot->DrawCopy("P");
		  else     hplot->DrawCopy("samesP");
		  if(j==0) leg->AddEntry(hplot, Form("%s",gTrgSetupTitle[k+1]), "P");
		}
	      TPaveText *t1 = GetTitleText(Form("%s_%s (%d-%d%%, %s)",run_type,dataTitle[d],(16-centBins_high[j+1])*5,(17-centBins_low[j+1])*5,legTitle), 0.055);
	      t1->Draw();
	    }
	  cCent[d][i]->cd(1);
	  leg->Draw();
	  if(savePlot) cCent[d][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%sVsCent_%s.pdf",run_type,typeName[i],dataTitle[d]));
	}
    }
}

//================================================
void tightCutsEmb(const int savePlot = 0)
{
  const char* trigTitle = "dimuon";
  TFile *fin = 0x0;
  char *trigName = "";
  char *libVer = "";

  if(trigTitle=="vpd30")
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed.Pion.P17id.root", "read");
      trigName = "mb";
      libVer = "P17id";
    }
  else if(trigTitle=="dimuon")
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
      trigName = "di_mu";
      libVer = "P15ie";
    }

  const int nType = 3;
  const char* typeName[nType] = {"MC", "Tpc", "ElecCheck"};
  THnSparseF *hnTpc[nType];
  TH1F *hTrkPtVsCent[nType][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      hnTpc[i] = (THnSparseF*)fin->Get(Form("mhMcTrkPtEff_%s_%s",typeName[i],trigName));

      hnTpc[i]->GetAxis(0)->SetRangeUser(2, 20);
      hnTpc[i]->GetAxis(1)->SetRangeUser(-0.5, 0.5);
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hnTpc[i]->GetAxis(4)->SetRange(k+1, k+1);
	  hTrkPtVsCent[i][k] = (TH1F*)hnTpc[i]->Projection(2);
	  hTrkPtVsCent[i][k]->SetName(Form("hTrkPtVsCent_%s%s",typeName[i],gTrgSetupTitle[k+1]));
	  hTrkPtVsCent[i][k]->Sumw2();
	  hTrkPtVsCent[i][k]->Rebin(2);
	}
      hnTpc[i]->GetAxis(0)->SetRange(0,-1);
      hnTpc[i]->GetAxis(1)->SetRange(0,-1);
      hnTpc[i]->GetAxis(4)->SetRange(0,-1);
    }

  TList *list = new TList;
  TString legName_trg[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      legName_trg[k] = Form("AuAu_200%s",gTrgSetupTitle[k+1]);
    }
  TH1F *hTrkEffVsCent[2][gNTrgSetup-1];
  const char *typeTitle[2] = {"standard cuts", "tighter cuts"};
  for(int i=0; i<2; i++)
    {    
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hTrkEffVsCent[i][k] = (TH1F*)hTrkPtVsCent[i+1][k]->Clone(Form("hTrkEffVsCent_%s%s",typeName[i+1],gTrgSetupTitle[k+1]));
	  hTrkEffVsCent[i][k]->Divide(hTrkPtVsCent[0][k]);
	  hTrkEffVsCent[i][k]->SetMarkerStyle(20+k);
	  hTrkEffVsCent[i][k]->SetMarkerColor(color[k]);
	  hTrkEffVsCent[i][k]->SetLineColor(color[k]);
	  for(int bin=1; bin<=hTrkEffVsCent[i][k]->GetNbinsX(); bin++)
	    {
	      hTrkEffVsCent[i][k]->GetXaxis()->SetBinLabel(bin,Form("%d-%d%%",80-bin*10,90-bin*10));
	    }
	  list->Add(hTrkEffVsCent[i][k]);
	}
      c = drawHistos(list,Form("TpcEff_vs_Cent_%s",typeName[i+1]),Form("%s: TPC tracking efficiency above 2 GeV/c (%s);;Efficiency",trigTitle,typeTitle[i]),false,2.0,3.8,true,0,1,kFALSE,kTRUE,legName_trg,true,"|#eta| < 0.5",0.5,0.7,0.2,0.45,kTRUE,0.04,0.04); 
      c->SetGrid(0,1);
      list->Clear();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Tracking_%s_TpcEffVsCent_Cut%d.pdf",run_type,trigName,i));
    }
}
