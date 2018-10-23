const int year = YEAR;

//================================================
void check_EmbedLumi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);


  tofMult();
}

//================================================
void tofMult(const int savePlot = 0)
{
  TFile *fin = NULL;
  char *trig_name = "";
  if(iMbEmb)
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed_MB.Jpsi.root", "read");
      trig_name = "mb";
    }
  else
    {
      fin = TFile::Open("output/Run14_AuAu200.Embed.StudyCent.root", "read");
      trig_name = "di_mu";
    }

  // global properties
  THnSparseF *hn = (THnSparseF*)fin->Get(Form("mhTofMult_%s",trig_name));
  TH2F *hTofMultVsRefMult[4];
  TH2F *hTofMultVsgTrack[4];
  TH2F *hgTrackVsRefMult[4];
  for(int i=0; i<4; i++)
    {
      hn->GetAxis(2)->SetRangeUser(20+i*20+1e-4,40+i*20-1e-4);
      hTofMultVsRefMult[i] = (TH2F*)hn->Projection(4,3);
      hTofMultVsRefMult[i]->SetName(Form("hTofMultVsRefMult_Zdc%d-%d",20+i*20,40+i*20));
      hTofMultVsgTrack[i] = (TH2F*)hn->Projection(4,5);
      hTofMultVsgTrack[i]->SetName(Form("hTofMultVsgTrack_Zdc%d-%d",20+i*20,40+i*20));
      hgTrackVsRefMult[i] = (TH2F*)hn->Projection(3,5);
      hgTrackVsRefMult[i]->SetName(Form("hgTrackVsRefMult_Zdc%d-%d",20+i*20,40+i*20));
    }

  TCanvas *c = new TCanvas("cTofMultVsRefMult", "cTofMultVsRefMult", 1100, 700);
  c->Divide(2,2);
  for(int i=0; i<4; i++)
    {
      c->cd(i+1);
      hTofMultVsRefMult[i]->SetYTitle("# of TOF hits");
      hTofMultVsRefMult[i]->GetYaxis()->SetTitleOffset(1.2);
      hTofMultVsRefMult[i]->SetTitle("");
      hTofMultVsRefMult[i]->GetXaxis()->SetRangeUser(0,700);
      hTofMultVsRefMult[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d < ZDCrate < %d kHz",20+i*20,40+i*20),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Cent_Embed_TofMultVsRefMult.png",run_type.Data()));
  

  c = new TCanvas("cTofMultVsgTrack", "cTofMultVsgTrack", 1100, 700);
  c->Divide(2,2);
  for(int i=0; i<4; i++)
    {
      c->cd(i+1);
      hTofMultVsgTrack[i]->SetYTitle("# of TOF hits");
      hTofMultVsgTrack[i]->GetYaxis()->SetTitleOffset(1.2);
      hTofMultVsgTrack[i]->SetXTitle("# of global tracks");
      hTofMultVsgTrack[i]->GetXaxis()->SetTitleOffset(1.2);
      hTofMultVsgTrack[i]->SetTitle("");
      hTofMultVsgTrack[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d < ZDCrate < %d kHz",20+i*20,40+i*20),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Cent_%sEmbed_TofMultVsgTrack.png",run_type.Data(),gEmbTypeName[iMbEmb]));


  TH1F *hgTrackInTofMult[4][4];
  const double lowTofMult[4] = {0, 400, 1000, 2800};
  const double highTofMult[4] = {200, 1000, 2000, 3000};
  c = new TCanvas("cgTrackInTofMult", "cgTrackInTofMult", 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.55,0.55,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int j=0; j<4; j++)
    {
      c->cd(j+1);
      gPad->SetLogy();
      for(int i=0; i<4; i++)
	{
	  int low_bin = hTofMultVsgTrack[i]->GetYaxis()->FindBin(lowTofMult[j]);
	  int high_bin =  hTofMultVsgTrack[i]->GetYaxis()->FindBin(highTofMult[j]);
	  hgTrackInTofMult[i][j] = (TH1F*)hTofMultVsgTrack[i]->ProjectionX(Form("hgTrackInTofMult%d_Zdc%d-%d",j,20+i*20,40+i*20),low_bin,high_bin);
	  hgTrackInTofMult[i][j]->SetLineColor(gColor[i]);
	  hgTrackInTofMult[i][j]->Sumw2();
	  hgTrackInTofMult[i][j]->Rebin(20);
	  hgTrackInTofMult[i][j]->Scale(1./hgTrackInTofMult[i][j]->Integral());
	  hgTrackInTofMult[i][j]->SetMarkerStyle(20+i);
	  hgTrackInTofMult[i][j]->SetMarkerColor(gColor[i]);
	  if(i==0) hgTrackInTofMult[i][j]->Draw("P");
	  else     hgTrackInTofMult[i][j]->Draw("PEsames");
	  if(j==0) leg->AddEntry(hgTrackInTofMult[i][j], Form("%d < ZDC < %d kHz",20+i*20,40+i*20), "P");
	}
      TPaveText *t1 = GetTitleText(Form("%1.0f < TofMult < %1.0f",lowTofMult[j],highTofMult[j]),0.06);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();

  c = new TCanvas("cgTrackVsRefMult", "cgTrackVsRefMult", 1100, 700);
  c->Divide(2,2);
  for(int i=0; i<4; i++)
    {
      c->cd(i+1);
      hgTrackVsRefMult[i]->SetXTitle("# of global tracks");
      hgTrackVsRefMult[i]->GetYaxis()->SetTitleOffset(1.2);
      hgTrackVsRefMult[i]->SetTitle("");
      hgTrackVsRefMult[i]->Draw("colz");
      TPaveText *t1 = GetTitleText(Form("%d < ZDCrate < %d kHz",20+i*20,40+i*20),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Cent_%sEmbed_gTrackVsRefMult.png",run_type.Data(),gEmbTypeName[iMbEmb]));

  // tracking efficiency
  THnSparseF *hnTrkEff = (THnSparseF*)fin->Get(Form("mhTrkEffVsTofMult_%s",trig_name));
  const int nZdcBins = 4;
  const int lowZdcBins[nZdcBins] = {20, 40, 60, 80};
  const int highZdcBins[nZdcBins] = {40, 60, 80, 100};
  const double minVz = -100, maxVz = 100;
  const double ptCut = 1.5;
  const char* multName[3] = {"gRefMultCorr", "TofMult", "ngTrack"};
  const int rebin = 20;
  TH1F *hMcTrkInMult[3][nZdcBins];
  TH1F *hRcTrkInMult[3][nZdcBins];
  TH1F *hTrkEffInMult[3][nZdcBins];
  
  hnTrkEff->GetAxis(1)->SetRangeUser(ptCut, 100);
  hnTrkEff->GetAxis(2)->SetRangeUser(minVz, maxVz);
  for(int i=0; i<nZdcBins; i++)
    {
      hnTrkEff->GetAxis(3)->SetRangeUser(lowZdcBins[i]+1e-4, highZdcBins[i]-1e-4);
      for(int j=0; j<3; j++)
	{
	  hMcTrkInMult[j][i] = (TH1F*)hnTrkEff->Projection(j+4);
	  hMcTrkInMult[j][i]->SetName(Form("hMcTrk_%s_Zdc%d-%d",multName[j],lowZdcBins[i],highZdcBins[i]));
	  hMcTrkInMult[j][i]->Sumw2();
	  if(j==0) hMcTrkInMult[j][i]->Rebin(rebin/2);
	  else     hMcTrkInMult[j][i]->Rebin(rebin);

	  hnTrkEff->GetAxis(0)->SetRange(2,2);
	  hRcTrkInMult[j][i] = (TH1F*)hnTrkEff->Projection(j+4);
	  hRcTrkInMult[j][i]->SetName(Form("hRcTrk_%s_Zdc%d-%d",multName[j],lowZdcBins[i],highZdcBins[i]));
	  hRcTrkInMult[j][i]->Sumw2();
	  if(j==0) hRcTrkInMult[j][i]->Rebin(rebin/2);
	  else     hRcTrkInMult[j][i]->Rebin(rebin);
	  hnTrkEff->GetAxis(0)->SetRange(0,-1);

	  hTrkEffInMult[j][i] = (TH1F*)hRcTrkInMult[j][i]->Clone(Form("hTrkEff_%s_Zdc%d-%d",multName[j],lowZdcBins[i],highZdcBins[i]));
	  hTrkEffInMult[j][i]->Divide(hMcTrkInMult[j][i]);
	}
    }
  c = new TCanvas("cTrkEffInMult", "cTrkEffInMult", 1100, 700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2,0.35,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  leg->SetHeader(Form("p_{T}>%1.1f GeV/c, |#eta|<0.8, |v_{z}|<%1.0f cm",ptCut,maxVz));
  for(int j=0; j<3; j++)
    {
      c->cd(j+1);
      for(int i=0; i<nZdcBins; i++)
	{
	  hTrkEffInMult[j][i]->SetMarkerStyle(20+i);
	  hTrkEffInMult[j][i]->SetMarkerColor(gColor[i]);
	  hTrkEffInMult[j][i]->SetLineColor(gColor[i]);
	  hTrkEffInMult[j][i]->GetYaxis()->SetRangeUser(0,1);
	  hTrkEffInMult[j][i]->SetTitle("");
	  if(j==0) hTrkEffInMult[j][i]->GetXaxis()->SetRangeUser(0,800);

	  TH1F *htmp = (TH1F*)hTrkEffInMult[j][i]->Clone(Form("%s_clone",hTrkEffInMult[j][i]->GetName()));
	  htmp->Divide(hTrkEffInMult[j][0]);
	  if(i==0) htmp->Draw();
	  else     htmp->Draw("sames");
	  if(j==0) leg->AddEntry(hTrkEffInMult[j][i], Form("%d < ZDCrate < %d kHz",20+i*20,40+i*20), "P");
	}
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency vs. %s",multName[j]),0.06);
      t1->Draw();
    }
  c->cd(4);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/Cent_%sEmbed_TpcEffVsMult_Vz%1.0f.pdf",run_type.Data(),gEmbTypeName[iMbEmb],maxVz));
}
