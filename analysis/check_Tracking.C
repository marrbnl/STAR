//================================================
void check_Tracking()
{
  gStyle->SetOptStat(0);

  embQA();
  //dataQA();
  //tightCutsEmb();
}

//================================================
void dataQA(const int savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  //----------------------------------------------
  TFile *fin = TFile::Open("output/Run14_AuAu200.Study.TpcTracking.root", "read");
  const char* dataTitle = "Run14_AuAu200_MTD";
  const char* dataName = "DataMtd";

  THnSparseF *hnQa = (THnSparseF*)fin->Get("mhnTrkPt_di_mu");
  const int nType = 2;
  const char* typeName[nType] = {"NHitsFit", "Dca"};

  hnQa->GetAxis(0)->SetRangeUser(2, 10);
  TH1F *hTrkDis[nType][nCentBins-1][gNTrgSetup-1];
  for(int j=0; j<nCentBins-1; j++)
    {
      hnQa->GetAxis(3)->SetRange(centBins_low[j+1], centBins_high[j+1]);
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hnQa->GetAxis(4)->SetRange(k+1, k+1);
	  for(int i=0; i<nType; i++)
	    {
	      hTrkDis[i][j][k] = (TH1F*)hnQa->Projection(i+1);
	      hTrkDis[i][j][k]->SetName(Form("Trk%s_cent%s%s",typeName[i],cent_Title[j+1], gTrgSetupTitle[k+1]));
	      hTrkDis[i][j][k]->Scale(1./hTrkDis[i][j][k]->Integral());
	    }
	}
      hnQa->GetAxis(4)->SetRange(0, -1);
    }
  hnQa->GetAxis(3)->SetRange(0, -1);

  // in centrality bins
  TCanvas *cCent[nType];
  for(int i=0; i<nType; i++)
    {
      cCent[i] = new TCanvas(Form("cCent_%s",typeName[i]), Form("cCent_%s",typeName[i]), 1100, 700);
      cCent[i]->Divide(2, 2);
      TLegend *leg =  0x0;
      if(i==0) leg = new TLegend(0.2, 0.45, 0.4, 0.8);
      else     leg = new TLegend(0.6, 0.45, 0.8, 0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader("p_{T} > 2 GeV/c");
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cCent[i]->cd(k+1);
	  if(i==1) gPad->SetLogy();
	  for(int j=0; j<nCentBins-1; j++)
	    {
	      TH1F *hplot = (TH1F*)hTrkDis[i][j][k]->Clone(Form("%s_clone",hTrkDis[i][j][k]->GetName()));
	      hplot->SetMarkerStyle(22+j);
	      hplot->SetMarkerColor(color[j]);
	      hplot->SetLineColor(color[j]);
	      hplot->SetMarkerSize(1.2);
	      hplot->SetTitle("");
	      if(i==0) 
		{
		  hplot->GetYaxis()->SetRangeUser(0,0.12);
		}
	      if(i==1) 
		{
		  hplot->GetXaxis()->SetRangeUser(0,1.5);
		  hplot->GetYaxis()->SetRangeUser(5e-3, 1);
		}
	      if(j==0) hplot->DrawCopy("P");
	      else     hplot->DrawCopy("samesP");
	      if(k==0) leg->AddEntry(hplot, Form("%s%%",cent_Name[j+1]), "P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s%s",dataTitle,gTrgSetupTitle[k+1]), 0.05);
	  t1->Draw();
	}
      cCent[i]->cd(1);
      leg->Draw();
      if(savePlot) cCent[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/%s_%sVsCent.pdf",run_type,dataName,typeName[i]));
    }
  return;

  // in luminosity bins
  TCanvas *cLumi[nType];
  for(int i=0; i<nType; i++)
    {
      cLumi[i] = new TCanvas(Form("cLumi_%s",typeName[i]), Form("cLumi_%s",typeName[i]), 1100, 700);
      cLumi[i]->Divide(2, 2);
      TLegend *leg =  0x0;
      if(i<3) leg = new TLegend(0.25, 0.5, 0.45, 0.8);
      else    leg = new TLegend(0.6, 0.5, 0.8, 0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      for(int j=0; j<nCentBins-1; j++)
	{
	  cLumi[i]->cd(j+1);
	  if(i==1 || i==3) gPad->SetLogy();
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hTrkDis[i][j][k]->SetMarkerStyle(20+k);
	      hTrkDis[i][j][k]->SetMarkerColor(color[k]);
	      hTrkDis[i][j][k]->SetLineColor(color[k]);
	      hTrkDis[i][j][k]->SetMarkerSize(1.2);
	      hTrkDis[i][j][k]->SetTitle("");
	      if(i==3) hTrkDis[i][j][k]->GetXaxis()->SetRangeUser(0,1.5);
	      if(i==0 || i==2) 
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(0,0.08);
		}
	      else if(i==1)
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(1e-4, 10);
		}
	      else
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(1e-3, 1);
		}
	      if(k==0) hTrkDis[i][j][k]->DrawCopy("P");
	      else     hTrkDis[i][j][k]->DrawCopy("samesP");
	      if(j==0) leg->AddEntry(hTrkDis[i][j][k], Form("%s",gTrgSetupTitle[k+1]), "P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s %s%%",run_type,cent_Name[j+1]), 0.05);
	  t1->Draw();
	}
      cLumi[i]->cd(1);
      leg->Draw();
      if(savePlot) cLumi[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Embed_%sVsLumi.pdf",run_type,typeName[i]));
    }
}

//================================================
void trackQA(const int savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  //----------------------------------------------
  const int nData = 3;
  const int nType = 2;
  const char* typeName[nType] = {"NHitsFit", "Dca"};

  TFile *fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");

  THnSparseF *hnQa[nType];
  TH1F *hTrkDis[nType][nCentBins-1][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      hnQa[i] = (THnSparseF*)fin->Get(Form("mhQa%s_di_mu",typeName[i]));
      hnQa[i]->GetAxis(0)->SetRangeUser(2, 20);
      for(int j=0; j<nCentBins-1; j++)
	{
	  hnQa[i]->GetAxis(3)->SetRange(centBins_low[j+1], centBins_high[j+1]);
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hnQa[i]->GetAxis(4)->SetRange(k+1, k+1);
	      hTrkDis[i][j][k] = (TH1F*)hnQa[i]->Projection(1);
	      hTrkDis[i][j][k]->SetName(Form("Trk%s_cent%s%s",typeName[i],cent_Title[j+1], gTrgSetupTitle[k+1]));
	      hTrkDis[i][j][k]->Sumw2();
	      hTrkDis[i][j][k]->Scale(1./hTrkDis[i][j][k]->Integral());
	    }
	  hnQa[i]->GetAxis(4)->SetRange(0, -1);
	}
      hnQa[i]->GetAxis(3)->SetRange(0, -1);
    }

  // in centrality bins
  TCanvas *cCent[nType];
  for(int i=0; i<nType; i++)
    {
      cCent[i] = new TCanvas(Form("cCent_%s",typeName[i]), Form("cCent_%s",typeName[i]), 1100, 700);
      cCent[i]->Divide(2, 2);
      TLegend *leg =  0x0;
      if(i<3) leg = new TLegend(0.2, 0.45, 0.4, 0.8);
      else    leg = new TLegend(0.6, 0.45, 0.8, 0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader("p_{T} > 2 GeV/c, |#eta| < 0.5");
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  cCent[i]->cd(k+1);
	  if(i==1 || i==3) gPad->SetLogy();
	  for(int j=0; j<nCentBins-1; j++)
	    {
	      TH1F *hplot = (TH1F*)hTrkDis[i][j][k]->Clone(Form("%s_clone",hTrkDis[i][j][k]->GetName()));
	      hplot->SetMarkerStyle(22+j);
	      hplot->SetMarkerColor(color[j]);
	      hplot->SetLineColor(color[j]);
	      hplot->SetMarkerSize(1.2);
	      hplot->SetTitle("");
	      if(i==3) hplot->GetXaxis()->SetRangeUser(0,1.5);
	      if(i==0 || i==2) 
		{
		  hplot->GetYaxis()->SetRangeUser(0,0.08);
		}
	      else if(i==1)
		{
		  hplot->GetYaxis()->SetRangeUser(1e-4, 10);
		}
	      else
		{
		  hplot->GetYaxis()->SetRangeUser(1e-3, 1);
		}
	      if(j==0) hplot->DrawCopy("P");
	      else     hplot->DrawCopy("samesP");
	      if(k==0) leg->AddEntry(hplot, Form("%s%%",cent_Name[j+1]), "P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s%s",run_type,gTrgSetupTitle[k+1]), 0.05);
	  t1->Draw();
	}
      cCent[i]->cd(1);
      leg->Draw();
      if(savePlot) cCent[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Embed_%sVsCent.pdf",run_type,typeName[i]));
    }

  // in luminosity bins
  TCanvas *cLumi[nType];
  for(int i=0; i<nType; i++)
    {
      cLumi[i] = new TCanvas(Form("cLumi_%s",typeName[i]), Form("cLumi_%s",typeName[i]), 1100, 700);
      cLumi[i]->Divide(2, 2);
      TLegend *leg =  0x0;
      if(i<3) leg = new TLegend(0.2, 0.45, 0.4, 0.8);
      else    leg = new TLegend(0.6, 0.45, 0.8, 0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.045);
      leg->SetHeader("p_{T} > 2 GeV/c, |#eta| < 0.5");
      for(int j=0; j<nCentBins-1; j++)
	{
	  cLumi[i]->cd(j+1);
	  if(i==1 || i==3) gPad->SetLogy();
	  for(int k=0; k<gNTrgSetup-1; k++)
	    {
	      hTrkDis[i][j][k]->SetMarkerStyle(20+k);
	      hTrkDis[i][j][k]->SetMarkerColor(color[k]);
	      hTrkDis[i][j][k]->SetLineColor(color[k]);
	      hTrkDis[i][j][k]->SetMarkerSize(1.2);
	      hTrkDis[i][j][k]->SetTitle("");
	      if(i==3) hTrkDis[i][j][k]->GetXaxis()->SetRangeUser(0,1.5);
	      if(i==0 || i==2) 
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(0,0.08);
		}
	      else if(i==1)
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(1e-4, 10);
		}
	      else
		{
		  hTrkDis[i][j][k]->GetYaxis()->SetRangeUser(1e-3, 1);
		}
	      if(k==0) hTrkDis[i][j][k]->DrawCopy("P");
	      else     hTrkDis[i][j][k]->DrawCopy("samesP");
	      if(j==0) leg->AddEntry(hTrkDis[i][j][k], Form("%s",gTrgSetupTitle[k+1]), "P");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s %s%%",run_type,cent_Name[j+1]), 0.05);
	  t1->Draw();
	}
      cLumi[i]->cd(1);
      leg->Draw();
      if(savePlot) cLumi[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Embed_%sVsLumi.pdf",run_type,typeName[i]));
    }
}

//================================================
void tightCutsEmb(const int savePlot = 1)
{
  TFile *fin = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root", "read");
  const int nType = 3;
  const char* typeName[nType] = {"MC", "Tpc", "ElecCheck"};
  THnSparseF *hnTpc[nType];
  TH1F *hTrkPtVsCent[nType][gNTrgSetup-1];
  for(int i=0; i<nType; i++)
    {
      hnTpc[i] = (THnSparseF*)fin->Get(Form("mhMcTrkPtEff_%s_di_mu",typeName[i]));

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
      c = drawHistos(list,Form("TpcEff_vs_Cent_%s",typeName[i+1]),Form("TPC tracking efficiency above 2 GeV/c (%s);;Efficiency",typeTitle[i]),false,2.0,3.8,true,0,1,kFALSE,kTRUE,legName_trg,true,"|#eta| < 0.5",0.5,0.7,0.2,0.45,kTRUE,0.04,0.04); 
      c->SetGrid(0,1);
      list->Clear();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Embed_TpcEffVsCent_Cut%d.pdf",run_type,i));
    }
}
