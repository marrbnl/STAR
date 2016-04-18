const int year = YEAR;

//================================================
void ana_JpsiEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }

  //ploEff();
  //plotEmbedEff();
  makeMtdRespEff();
  //makeTrigEff();
  //makeEmbedEff();
  //hadronEmbed();
  //compTrigEff();
}

//================================================
void ploEff(const int savePlot = 0)
{
  TList *list = new TList;
  TString legName[nCentBins];

  char *trigEffName = "";
  char *embedEffName = "";
  if(year==2013) 
    {
      embedEffName = Form("Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
      trigEffName = Form("Run13.pp500.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
    }
  if(year==2014)
    {
      embedEffName = Form("Run14.AuAu200.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
      trigEffName = Form("Run14.AuAu200.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
    }

  TFile *fEmbedEff = TFile::Open(Form("Rootfiles/%s",embedEffName),"read");
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/%s",trigEffName),"read");
  TH1F *hJpsiEffEmbed[nCentBins];
  TH1F *hJpsiSmearEff[nCentBins];
  TH1F *hJpsiEffTrig[nCentBins];
  TH1F *hJpsiRespEff[nCentBins];
  TH1F *hJpsiEffTotal[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEffEmbed[k] = (TH1F*)fEmbedEff->Get(Form("MTDreco_Jpsi_pT_%s_WeightPt_Eff_rebin",cent_Title[k]));
      hJpsiEffEmbed[k]->Scale(0.8/0.5);
      hJpsiSmearEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiSmearEff_cent%s_rebin",cent_Title[k]));
      hJpsiEffTrig[k] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s_rebin",cent_Title[k]));
      hJpsiRespEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiRespEff_cent%s_rebin",cent_Title[k]));
      hJpsiEffTotal[k] = (TH1F*)hJpsiEffEmbed[k]->Clone(Form("JpsiEffTotal_cent%s_rebin",cent_Title[k]));
      if(hJpsiSmearEff[k]) hJpsiEffTotal[k]->Multiply(hJpsiSmearEff[k]);
      hJpsiEffTotal[k]->Multiply(hJpsiEffTrig[k]);
      hJpsiEffTotal[k]->Multiply(hJpsiRespEff[k]);
    }

  TCanvas *c = new TCanvas("JpsiTrigEff","JpsiTrigEff",1100,700);
  c->Divide(2,2);
  TH1F *hJpsiTrigEff[nCentBins][2];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiTrigEff[k][0] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s_rebin",cent_Title[k]));
      hJpsiTrigEff[k][1] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s",cent_Title[k]));
      c->cd(k+1);
      hJpsiTrigEff[k][0]->GetYaxis()->SetRangeUser(0.5,0.9);
      hJpsiTrigEff[k][0]->DrawCopy();
      hJpsiTrigEff[k][1]->DrawCopy("sames");
    }

  TList *list = new TList;
  TString legName[nCentBins];
  bool drawLegend = true;
  if(nCentBins==1) drawLegend = false;
  // smearing efficiency
  if(hJpsiSmearEff[0])
    {
      list->Clear();
      for(int k=0; k<nCentBins; k++)
	{
	  TH1F *htmp = (TH1F*)hJpsiSmearEff[k]->Clone(Form("%s_clone",hJpsiSmearEff[k]->GetName()));
	  list->Add(htmp);
	  legName[k] = Form("%s%%",cent_Name[k]);
	}
      c = drawHistos(list,"SmearEff","Correction factor due to smearing;p_{T,J/#psi} (GeV/c);Eff",kFALSE,0,30,kTRUE,0.9,1.1,kFALSE,drawLegend,legName,drawLegend,"Centrality",0.15,0.3,0.6,0.85,kTRUE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffSmear_CentBins.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffSmear_CentBins.png",run_type));
	}
    }

  // response efficiency
  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      TH1F *htmp = (TH1F*)hJpsiRespEff[k]->Clone(Form("%s_clone",hJpsiRespEff[k]->GetName()));
      list->Add(htmp);
      legName[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,"RespEff","Efficiency for MTD response;p_{T,J/#psi} (GeV/c);Eff",kFALSE,0,30,kTRUE,0.5,0.9,kFALSE,drawLegend,legName,drawLegend,"Centrality",0.15,0.3,0.6,0.85,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffResp_CentBins.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffResp_CentBins.png",run_type));
    }

  int dx = 1100, dy = 700;
  if(nCentBins==1) { dx = 800; dy = 600; }
  TCanvas *c = new TCanvas("JpsiEffTotal","JpsiEffTotal",dx,dy);
  if(nCentBins>1) c->Divide(nCentBins/2,2);
  TLegend *leg = new TLegend(0.15,0.5,0.25,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      hJpsiEffTotal[k]->GetYaxis()->SetRangeUser(0.0001,1.2);
      hJpsiEffTotal[k]->SetTitle(";p_{T,J/#psi} (GeV/c);Efficiency");
      hJpsiEffTotal[k]->SetMarkerStyle(20);
      hJpsiEffTotal[k]->SetMarkerColor(1);
      hJpsiEffTotal[k]->SetLineColor(1);
      hJpsiEffTotal[k]->SetMarkerSize(1.5);
      hJpsiEffTotal[k]->Draw();

      hJpsiEffEmbed[k]->SetMarkerStyle(21);
      hJpsiEffEmbed[k]->SetMarkerColor(2);
      hJpsiEffEmbed[k]->SetLineColor(2);
      hJpsiEffEmbed[k]->SetMarkerSize(1.5);
      hJpsiEffEmbed[k]->Draw("sames");

      if(hJpsiSmearEff[k])
	{
	  hJpsiSmearEff[k]->SetMarkerStyle(25);
	  hJpsiSmearEff[k]->SetMarkerColor(kGreen+2);
	  hJpsiSmearEff[k]->SetLineColor(kGreen+2);
	  hJpsiSmearEff[k]->SetMarkerSize(1.5);
	  hJpsiSmearEff[k]->Draw("sames");
	}

      hJpsiRespEff[k]->SetMarkerStyle(22);
      hJpsiRespEff[k]->SetMarkerColor(4);
      hJpsiRespEff[k]->SetLineColor(4);
      hJpsiRespEff[k]->SetMarkerSize(1.5);
      hJpsiRespEff[k]->Draw("sames");

      if(year!=2013)
	{
	  hJpsiEffTrig[k]->SetMarkerStyle(23);
	  hJpsiEffTrig[k]->SetMarkerColor(6);
	  hJpsiEffTrig[k]->SetLineColor(6);
	  hJpsiEffTrig[k]->SetMarkerSize(1.5);
	  hJpsiEffTrig[k]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("J/#psi efficiency in %s%%",cent_Name[k]),0.05);
      if(nCentBins==1) t1 = GetTitleText(Form("J/#psi efficiency"));
      t1->Draw();

      if(k==0)
	{
	  leg->AddEntry(hJpsiEffTotal[k],"Total efficiency","PL");	    
	  leg->AddEntry(hJpsiEffEmbed[k],"Tracking, PID, acc (|#eta_{mc}^{J/#psi}|<0.5)","PL");
	  if(hJpsiSmearEff[k]) leg->AddEntry(hJpsiSmearEff[k],"Smear embedding","PL");
	  leg->AddEntry(hJpsiRespEff[k],"MTD response","PL");
	  if(year!=2013) leg->AddEntry(hJpsiEffTrig[k],"MTD trigger efficiency","PL");
	}
    }
  c->cd(1);
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffTotal_CentBins.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffTotal_CentBins.png",run_type));
    }
}

//================================================
void plotEmbedEff(const int savePlot = 1, const int saveHisto = 1)
{
  TFile *fEmbed = 0x0, *fData = 0x0;

  if(year==2013)
    {
      fEmbed = TFile::Open(Form("output/Run13.pp500.jpsi.Embed.%sroot",run_config),"read");
      fData = TFile::Open(Form("output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.TightTof.root","read");
      fData = TFile::Open("output/Pico.Run14.AuAu200.jpsi.TightTof.root","read");
    }

  // compare global event characteristic between data and embedding
  TList *list = new TList;
  const TString legName[2] = {"Embedding","Data"};
  if(year==2014)
    {
      TH1F *hCent[2];
      hCent[0] = (TH1F*)fEmbed->Get("mhCentrality_di_mu");
      hCent[0]->SetName("mhCentrality_embed");
      hCent[1] = (TH1F*)fData->Get("mhCentrality_di_mu");
      list->Clear();
      for(int i=0; i<2; i++)
	{
	  hCent[i]->Sumw2();
	  hCent[i]->Scale(1./hCent[i]->GetBinContent(16));
	  list->Add(hCent[i]);
	}
      c = drawHistos(list,"Compare_centrality","Centrality distribution",kFALSE,0,0,kFALSE,0,0,kTRUE,kTRUE,legName,kTRUE,"",0.2,0.4,0.65,0.8,kFALSE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.png",run_type));
	}

      TH1F *hgRefMult[2];
      hgRefMult[0] = (TH1F*)fEmbed->Get("mhgRefMult_di_mu");
      hgRefMult[0]->SetName("mhgRefMult_embed");
      hgRefMult[1] = (TH1F*)fData->Get("mhgRefMult_di_mu");
      list->Clear();
      for(int i=0; i<2; i++)
	{
	  hgRefMult[i]->Sumw2();
	  hgRefMult[i]->Scale(1./hgRefMult[i]->Integral());
	  list->Add(hgRefMult[i]);
	}
      c = drawHistos(list,"Compare_gRefMult","gRefMult distribution",kTRUE,0,800,kFALSE,0,0,kFALSE,kTRUE,legName,kTRUE,"",0.65,0.8,0.7,0.85,kFALSE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.png",run_type));
	}

      TH1F *hgRefMultCorr[2];
      hgRefMultCorr[0] = (TH1F*)fEmbed->Get("mhgRefMultCorr_di_mu");
      hgRefMultCorr[0]->SetName("mhgRefMultCorr_embed");
      hgRefMultCorr[1] = (TH1F*)fData->Get("mhgRefMultCorr_di_mu");
      list->Clear();
      for(int i=0; i<2; i++)
	{
	  hgRefMultCorr[i]->Sumw2();
	  hgRefMultCorr[i]->Scale(1./hgRefMultCorr[i]->Integral());
	  list->Add(hgRefMultCorr[i]);
	}
      c = drawHistos(list,"Compare_gRefMultCorr","Corrected gRefMult distribution",kTRUE,0,800,kFALSE,0,0,kFALSE,kTRUE,legName,kTRUE,"",0.65,0.8,0.7,0.85,kFALSE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.png",run_type));
	}
    }

  // compare vertex
  TH1F *hVertex[2];
  if(year==2013) hVertex[0] = (TH1F*)fEmbed->Get("mhDefaultVtxZ_di_mu");
  if(year==2014) hVertex[0] = (TH1F*)fEmbed->Get("mhChosenVtxZ_di_mu");
  hVertex[0]->Rebin(2);
  hVertex[1] = (TH1F*)fData->Get("mhTpcVz_di_mu");
  list->Clear();
  for(int i=0; i<2; i++)
    {
      hVertex[i]->Sumw2();
      hVertex[i]->Scale(1./hVertex[i]->Integral());
      list->Add(hVertex[i]);
    }
  c = drawHistos(list,"Compare_TpcVz","Z distribution of TPC vertex",kFALSE,0,800,kFALSE,0,0,kFALSE,kTRUE,legName,kTRUE,"",0.65,0.8,0.7,0.85,kFALSE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.png",run_type));
    }

  // Jpsi distribution
  char *inName = "";
  if(year==2013) inName = Form("Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  if(year==2014) inName = Form("Run14.AuAu200.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);

  TFile *fin = 0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s",inName),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s",inName),"read");
  TString name[3] = {"MCinput","TPCreco","MTDreco"};
  TString variable[4] = {"mass","pT","rapidity","phi"};
  const char *mc_name[4] = {"invariant mass","p_{T}","y","#varphi"};
  TH1F *hJpsi[nCentBins][3][4];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<4; j++)
	    {
	      hJpsi[k][i][j] = (TH1F*)fin->Get(Form("%s_Jpsi_%s_%s",name[i].Data(),variable[j].Data(),cent_Title[k]));
	      hJpsi[k][i][j]->Sumw2();
	      hJpsi[k][i][j]->SetMinimum(0.1);
	      hJpsi[k][i][j]->SetMarkerStyle(20+k);
	      hJpsi[k][i][j]->SetMarkerColor(color[k]);
	      hJpsi[k][i][j]->SetLineColor(color[k]);
	      hJpsi[k][i][j]->SetTitle("");
	    }
	}
    }

  for(int i=0; i<3; i++)
    {
      TCanvas *c = new TCanvas(Form("JpsiDis_%s",name[i].Data()),Form("JpsiDis_%s",name[i].Data()),1100,700);
      TLegend *leg = new TLegend(0.65,0.55,0.85,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      c->Divide(2,2);
      for(int j=0; j<4; j++)
	{
	  c->cd(j+1);
	  for(int k=0; k<nCentBins; k++)
	    {
	      if(j==0) 
		{
		  leg->AddEntry(hJpsi[k][i][j], Form("%s%%",cent_Name[k]),"P");
		  hJpsi[k][i][j]->GetXaxis()->SetRangeUser(2.9,3.3);
		  hJpsi[k][i][j]->SetXTitle("M_{#mu#mu} (GeV/c^{2})");
		}
	      if(k==0) hJpsi[k][i][j]->DrawCopy();
	      else     hJpsi[k][i][j]->DrawCopy("samesP");
	    }
	  TPaveText *t1 = GetTitleText(Form("%s: %s distribution of J/#psi",name[i].Data(),variable[j].Data()),0.05);
	  t1->Draw();
	}
      c->cd(1);
      if(nCentBins>1) leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%s%s_JpsiDis_CentBins.pdf",run_type,run_config,name[i].Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%s%s_JpsiDis_CentBins.png",run_type,run_config,name[i].Data()));
	}
    }

  // matched Jpsi for response
  TCanvas *c = new TCanvas("Check_JpsiMatch","Check_JpsiMatch",1100,700);
  c->Divide(2,2);
  TH2F *hJpsiPtMatch[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiPtMatch[k] = (TH2F*)fin->Get(Form("JpsiPt_TrueVsReco_%s",cent_Title[k]));
      c->cd(k+1);
      hJpsi[k][2][1]->Draw();
      TH1F *htmp = (TH1F*)hJpsiPtMatch[k]->ProjectionX(Form("%s_tmp",hJpsiPtMatch[k]->GetName()));
      htmp->Draw("samesHIST");
    }

 // weight input spectrum
  printf("+++ Weight jpsi pt distribution +++\n");
  TH1F *funcWeight[nCentBins];
  if(year==2013)
    {
      TFile *fWeight = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
      TF1 *funcJpsi = (TF1*)fWeight->Get("ffpt");
      funcJpsi->SetNpx(1000);
      for(int k=0; k<nCentBins; k++)
	{
	  funcWeight[k] = (TH1F*)funcJpsi->GetHistogram();
	  funcWeight[k]->SetName(Form("GlobalFit_Jpsi_Yield_cent%s",cent_Title[k]));
	  for(int bin=1; bin<=funcWeight[k]->GetNbinsX(); bin++)
	    {
	      funcWeight[k]->SetBinContent(bin,funcWeight[k]->GetBinCenter(bin)*funcWeight[k]->GetBinContent(bin));
	    }
	}  
    }
  if(year==2014)
    {
      TFile *fWeight = TFile::Open("Rootfiles/Published/Jpsi_Raa_200/Publication.Jpsi.200GeV.root","read");
      for(int k=0; k<nCentBins; k++)
	{
	  funcWeight[k] = (TH1F*)fWeight->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[k]));
	}  
    }

  TH1F *hJpsiPt[nCentBins][3][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  if(k<3) hJpsiPt[k][i][0] = (TH1F*)hJpsi[k][i][1]->Clone(Form("%s_FlatPt",hJpsi[k][i][1]->GetName()));
	  else    hJpsiPt[k][i][0] = (TH1F*)hJpsi[2][i][1]->Clone(Form("%s_FlatPt",hJpsi[k][i][1]->GetName()));
	}
    }

  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  // use 20-40% embedding data for 40-60% centrality bin
	  if(i<2) 
	    {
	      if(k<3) hJpsiPt[k][i][1] = (TH1F*)hJpsi[k][i][1]->Clone(Form("%s_WeightPt",hJpsi[k][i][1]->GetName()));
	      else    hJpsiPt[k][i][1] = (TH1F*)hJpsi[2][i][1]->Clone(Form("%s_WeightPt",hJpsi[k][i][1]->GetName()));
	    }

	  if(i==0)
	    {
	      for(int bin=1; bin<=hJpsiPt[k][i][1]->GetNbinsX(); bin++)
		{
		  double weight = funcWeight[k]->GetBinContent(funcWeight[k]->FindBin(hJpsiPt[k][i][1]->GetBinCenter(bin)));
		  hJpsiPt[k][i][1]->SetBinContent(bin,hJpsiPt[k][i][1]->GetBinContent(bin)*weight);
		  hJpsiPt[k][i][1]->SetBinError(bin,hJpsiPt[k][i][1]->GetBinError(bin)*weight);
		}
	    }
	  else if(i==2)
	    {
	      if(k<3) TH2F *h2tmp = (TH2F*)hJpsiPtMatch[k]->Clone(Form("%s_tmp",hJpsiPtMatch[k]->GetName()));
	      else    TH2F *h2tmp = (TH2F*)hJpsiPtMatch[2]->Clone(Form("%s_tmp",hJpsiPtMatch[k]->GetName()));

	      for(int biny=1; biny<=h2tmp->GetNbinsY(); biny++)
		{
		  double weight = funcWeight[k]->GetBinContent(funcWeight[k]->FindBin(h2tmp->GetYaxis()->GetBinCenter(biny)));
		  for(int binx=1; binx<=h2tmp->GetNbinsX(); binx++)
		    {
		      h2tmp->SetBinContent(binx,biny,weight*h2tmp->GetBinContent(binx,biny));
		      h2tmp->SetBinError(binx,biny,weight*h2tmp->GetBinError(binx,biny));
		    }
		}
	      hJpsiPt[k][i][1] = (TH1F*)h2tmp->ProjectionX(Form("%s_WeightPt",hJpsi[k][i][1]->GetName()));
	    }
	  hJpsiPt[k][i][1]->SetTitle("");
	}
    }

  printf("+++ Rebin jpsi distribution +++\n");
  // rebin input spectrum
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  TH1F *hJpsiPtRebin[nCentBins][3][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hJpsiPtRebin[k][i][j] = (TH1F*)hJpsiPt[k][i][j]->Rebin(nbins,Form("%s_rebin",hJpsiPt[k][i][j]->GetName()),xbins);
	      cout << hJpsiPtRebin[k][i][j]->GetName() << "  "  << i << "  " << j << "  " << hJpsiPtRebin[k][i][j]->GetBinContent(4) << endl;
	      scaleHisto(hJpsiPtRebin[k][i][j], 1, 1, kTRUE, kFALSE, kFALSE);
	      scaleHisto(hJpsiPt[k][i][j], 1, 1, kTRUE, kFALSE, kFALSE);
	    }
	}
    }

  // j/psi efficiency
  TH1F *hJpsiEff[nCentBins][2][2], *hJpsiEffRebin[nCentBins][2][2];
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hJpsiEff[k][i][j] = (TH1F*)hJpsiPt[k][i+1][j]->Clone(Form("%s_Eff",hJpsiPt[k][i+1][j]->GetName()));
	      hJpsiEff[k][i][j]->Divide(hJpsiPt[k][0][j]);
	      
	      hJpsiEffRebin[k][i][j] = (TH1F*)hJpsiPtRebin[k][i+1][j]->Clone(Form("%s_Eff_rebin",hJpsiPt[k][i+1][j]->GetName()));
	      hJpsiEffRebin[k][i][j]->Divide(hJpsiPtRebin[k][0][j]);
	      cout << hJpsiEffRebin[k][i][j]->GetName() << "  " << i << "  " << j << "  " << hJpsiEffRebin[k][i][j]->GetBinContent(4) << endl;
	    }
	}
    }

  TString legName2[2] = {"Flat pT","Weighted pT"};
  const char *legNameTitle[2] = {"FlatPt","WeightedPt"};
  int cx = 1100, cy = 700;
  if(nCentBins==1) { cx = 800; cy = 600; }
  for(int j=0; j<2; j++)
    {
      TCanvas *c = new TCanvas(Form("McInput_JpsiPt_%s",legName2[j].Data()),Form("McInput_JpsiPt_%s",legName2[j].Data()),cx,cy);
      if(nCentBins>1) c->Divide(2,2);
      for(int k=0; k<nCentBins; k++)
	{
	  c->cd(k+1);
	  if(nCentBins>1)
	    {
	      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	      ScaleHistoTitle(hJpsiPt[k][0][j],0.06,1,0.05,0.06,1,0.05,62);
	    }
	  else
	    {
	      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
	      ScaleHistoTitle(hJpsiPt[k][0][j],0.045,1,0.035,0.045,1,0.035,62);
	    }
	  hJpsiPt[k][0][j]->SetMarkerStyle(24);
	  hJpsiPt[k][0][j]->SetMarkerColor(color[k]);
	  hJpsiPt[k][0][j]->SetLineColor(color[k]);
	  hJpsiPt[k][0][j]->Draw();
	  hJpsiPtRebin[k][0][j]->SetMarkerStyle(21);
	  hJpsiPtRebin[k][0][j]->SetMarkerColor(color[k]);
	  hJpsiPtRebin[k][0][j]->SetLineColor(color[k]);
	  hJpsiPtRebin[k][0][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText("MC input p_{T} disribution of J/#psi");
	  t1->Draw();
	  TLegend *leg = new TLegend(0.5,0.5,0.6,0.8);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.04);
	  leg->SetHeader(Form("%s%%",cent_Name[k]));
	  leg->AddEntry(hJpsiPt[k][0][j],Form("%s",legName2[j].Data()),"P");
	  leg->AddEntry(hJpsiPtRebin[k][0][j],Form("%s rebinned",legName2[j].Data()),"P");
	  leg->Draw();
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MCJpsiPt_%s_CentBins.pdf",run_type,legNameTitle[j]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MCJpsiPt_%s_CentBins.png",run_type,legNameTitle[j]));
	}
    }

  TCanvas *c = new TCanvas(Form("JpsiPtEff"),Form("JpsiPtEff"),cx,cy);
  if(nCentBins>1) c->Divide(2,2);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      if(nCentBins>1) SetPadMargin(gPad,0.15,0.15,0.05,0.05);
      else            SetPadMargin(gPad,0.13,0.13,0.05,0.05);
      for(int j=0; j<2; j++)
	{
	  hJpsiEff[k][1][j]->SetMarkerStyle(24+j);
	  hJpsiEff[k][1][j]->SetMarkerColor(j+1);
	  hJpsiEff[k][1][j]->SetLineColor(j+1);
	  hJpsiEffRebin[k][1][j]->SetMarkerStyle(20+j);
	  hJpsiEffRebin[k][1][j]->SetMarkerColor(j+1);
	  hJpsiEffRebin[k][1][j]->SetLineColor(j+1);

	  if(nCentBins>1) ScaleHistoTitle(hJpsiEff[k][1][j],0.06,1,0.05,0.06,1,0.05,62);
	  else            ScaleHistoTitle(hJpsiEff[k][1][j],0.045,1,0.035,0.045,1,0.035,62);
	  hJpsiEff[k][1][j]->SetTitle(";p_{T} (GeV/c);Efficiency");
	  if(j==0) hJpsiEff[k][1][j]->DrawCopy();
	  else  hJpsiEff[k][1][j]->DrawCopy("sames");
	  hJpsiEffRebin[k][1][j]->DrawCopy("sames");
	}
      if(nCentBins>1)
	{
	  TPaveText *t1 = GetPaveText(0.7,0.8,0.2,0.3,0.06,62);
	  t1->AddText(Form("%s%%",cent_Name[k]));
	  t1->Draw();
	}
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.18,0.6,0.42,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hJpsiEff[0][1][0],"Flat pT","P");
  leg->AddEntry(hJpsiEff[0][1][1],"Weighted pT","P");
  leg->AddEntry(hJpsiEffRebin[0][1][0],"Flat pT rebinned","P");
  leg->AddEntry(hJpsiEffRebin[0][1][1],"Weighted pT rebinned","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiPtEff_CentBins.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiPtEff_CentBins.png",run_type,run_config));
    }

  if(nCentBins>1)
    {
      TList *list = new TList;
      TString legName3[nCentBins];
      for(int k=0; k<nCentBins; k++)
	{
	  TH1F *htmp = (TH1F*)hJpsiEffRebin[k][1][1]->Clone(Form("%s_tmp",hJpsiEffRebin[k][1][1]->GetName()));
	  list->Add(htmp);
	  legName3[k] = Form("%s%%",cent_Name[k]);
	}
      c = drawHistos(list,"JpsiEff_CentBins","J/#psi efficiency from embedding;p_{T} (GeV/c);Efficiency",kTRUE,0,10,kTRUE,0,0.03,kFALSE,kTRUE,legName3,kTRUE,"",0.2,0.4,0.6,0.85,kTRUE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiEffEmbed_CentBins.pdf",run_type,run_config));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiEffEmbed_CentBins.png",run_type,run_config));
	}
    }

  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  for(int i=0; i<3; i++)
	    {
	      for(int j=0; j<2; j++)
		{
		  hJpsiPt[k][i][j]->Write("",TObject::kOverwrite);
		  if(i>0) hJpsiEff[k][i-1][j]->Write("",TObject::kOverwrite);

		  hJpsiPtRebin[k][i][j]->Write("",TObject::kOverwrite);
		  if(i>0) hJpsiEffRebin[k][i-1][j]->Write("",TObject::kOverwrite);
		}
	    }
	}
    }
}

//================================================
void makeEmbedEff()
{
  gStyle->SetOptStat(1);
  TString fileName ;
  if(year==2013)
    {
      fileName = Form("Run13.pp500.jpsi.Embed.%sroot",run_config);
    }
  else if(year==2014)
    {
      fileName = Form("Run14.AuAu200.Jpsi.Embed.%sroot",run_config);
    }

  printf("Embedding file: %s\n",fileName.Data());

  TFile *f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  THnSparseF *hJpsiUS_mc = (THnSparseF*)f->Get("hJpsiInfo_di_mu");
  printf("+++ Jpsi distribution +++\n");
  // j/psi distribution
  TH1F *hJpsi[nCentBins][3][4];
  TH2F *hJpsiMassVsPt[nCentBins][3];
  TString name[3] = {"MCinput","TPCreco","MTDreco"};
  TString variable[4] = {"mass","pT","rapidity","phi"};
  const char *mc_name[4] = {"invariant mass","p_{T}","y","#varphi"};
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hJpsiUS_mc->GetAxis(10)->SetRange(centBins_low[k],centBins_high[k]);
      for(int i=0; i<3; i++)
	{
	  hJpsiUS_mc->GetAxis(7)->SetRange(i+4,i+4);
	  if(i>0)
	    {
	      hJpsiUS_mc->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
	      hJpsiUS_mc->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  hJpsiMassVsPt[k][i] = (TH2F*)hJpsiUS_mc->Projection(0,1);
	  hJpsiMassVsPt[k][i]->SetName(Form("%s_Jpsi_MassVsPt_%s",name[i].Data(),cent_Title[k]));
	  
	  for(int j=0; j<4; j++)
	    {
	      if(i>0 && j>0) hJpsiUS_mc->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
	      hJpsi[k][i][j] = (TH1F*)hJpsiUS_mc->Projection(j);
	      hJpsi[k][i][j]->SetName(Form("%s_Jpsi_%s_%s",name[i].Data(),variable[j].Data(),cent_Title[k]));
	      hJpsiUS_mc->GetAxis(0)->SetRange(0,-1);
	    }
	  hJpsiUS_mc->GetAxis(4)->SetRange(0,-1);
	  hJpsiUS_mc->GetAxis(5)->SetRange(0,-1);
	  hJpsiUS_mc->GetAxis(7)->SetRange(i+4,i+4);
	}
      hJpsiUS_mc->GetAxis(10)->SetRange(0,-1);
    }

  // matched Jpsi for response
  THnSparseF *hnJpsiMatch = (THnSparseF*)f->Get("mhJpsiMatch_di_mu");
  hnJpsiMatch->GetAxis(2)->SetRangeUser(low_mass+0.001, high_mass-0.001);
  hnJpsiMatch->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
  hnJpsiMatch->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
  TH2F *hJpsiPtMatch[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      if(nCentBins>1) hnJpsiMatch->GetAxis(6)->SetRange(centBins_low[k],centBins_high[k]);
      hJpsiPtMatch[k] = (TH2F*)hnJpsiMatch->Projection(1,3);
      hJpsiPtMatch[k]->Sumw2();
      hJpsiPtMatch[k]->SetName(Form("JpsiPt_TrueVsReco_%s",cent_Title[k]));
    }
  hnJpsiMatch->GetAxis(6)->SetRange(0,-1);
  hnJpsiMatch->GetAxis(2)->SetRange(0,-1);
  hnJpsiMatch->GetAxis(4)->SetRange(0,-1);
  hnJpsiMatch->GetAxis(5)->SetRange(0,-1);

  
  char *outName = "";
  if(year==2013) outName = Form("Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  if(year==2014) outName = Form("Run14.AuAu200.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);

  TFile *fout = TFile::Open(Form("Rootfiles/%s",outName),"recreate");
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiPtMatch[k]->Write();
      for(int i=0; i<3; i++)
	{
	  hJpsiMassVsPt[k][i]->SetTitle("");
	  hJpsiMassVsPt[k][i]->Write();
	  for(int j=0; j<4; j++)
	    {
	      hJpsi[k][i][j]->SetTitle("");
	      hJpsi[k][i][j]->Write();
	    }
	}
    }
  fout->Close();
}


//================================================
void makeMtdRespEff(bool savePlot = 0, bool saveHisto = 0)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.15);
  const char *name[2] = {"cosmic","embed"};
  TGraphAsymmErrors *hMtdResponseEff[2];  
  TFile *fin = 0x0, *fEmbed = 0x0;
  if(year==2014)
    {
      //cosmic
      fin = TFile::Open("Rootfiles/Run14.cosmic.MtdEff.root","read");
      hMtdResponseEff[0] = (TGraphAsymmErrors*)fin->Get("RateDown");
      fEmbed =  TFile::Open("output/Run14.AuAu200.Jpsi.Embed.LooseCut.root","read");
    }
  else if(year==2013)
    {
      fin = TFile::Open("Rootfiles/Run1314CosmicRayMTDResponseEfficiency.root","read");
      hMtdResponseEff[0] = (TGraphAsymmErrors*)fin->Get("Run13BOTVsPt");
      fEmbed =  TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
    }

  int npoints = hMtdResponseEff[0]->GetN();
  const int nbins = npoints - 1;
  double xbins[nbins+1];
  double x,y;
  for(int ip=0; ip<nbins; ip++)
    {
      hMtdResponseEff[0]->GetPoint(ip,x,y);
      xbins[ip] = x - hMtdResponseEff[0]->GetErrorXlow(ip);
      if(ip==nbins-1) xbins[ip+1] = x + hMtdResponseEff[0]->GetErrorXhigh(ip);
    }
  if(xbins[nbins]>10) xbins[nbins] = 10;

  // embedding
  TH2F *hProjVsPt = (TH2F*)fEmbed->Get("mhProjTrack");
  draw2D(hProjVsPt);
  TH2F *hMatchVsPt = (TH2F*)fEmbed->Get("mhMatchTrack");
  draw2D(hMatchVsPt);
  TH1F *htmp = (TH1F*)hProjVsPt->ProjectionX("hProj_tmp");
  TH1F *hProj = (TH1F*)htmp->Rebin(nbins,"hProj_1",xbins);
  htmp  = (TH1F*)hMatchVsPt->ProjectionX("hMatch_tmp");
  TH1F *hMatch = (TH1F*)htmp->Rebin(nbins,"hMatch_1",xbins);
  hMtdResponseEff[1] = new TGraphAsymmErrors();
  hMtdResponseEff[1]->Divide(hMatch,hProj,"cl=0.683 b(1,1) mode");

  TH1F *hEff[2];
  for(int i=0; i<2; i++)
    {
      hEff[i] = new TH1F(Form("MtdResponseEff_%s",name[i]), "MTD response efficiency;p_{T} (GeV/c);Efficiency", nbins, xbins);
      for(int ip=0; ip<nbins; ip++)
	{
	  hMtdResponseEff[i]->GetPoint(ip,x,y);
	  double error = hMtdResponseEff[i]->GetErrorYhigh(ip);
	  if(error<hMtdResponseEff[i]->GetErrorYlow(ip)) error = hMtdResponseEff[i]->GetErrorYlow(ip);
	  hEff[i]->SetBinContent(ip+1, y);
	  hEff[i]->SetBinError(ip+1, error);
	}
      hEff[i]->SetMarkerStyle(20+i);
      hEff[i]->SetLineColor(2-i);
      hEff[i]->SetMarkerColor(2-i);
      hEff[i]->GetXaxis()->SetRangeUser(0,10);
      hEff[i]->GetYaxis()->SetRangeUser(0,1.2);
    }
  ScaleHistoTitle(hEff[0],0.045,1,0.035,0.045,1,0.035,62);


  // Fit embedding data
  TH1F *hEmbedFit = (TH1F*)hEff[1]->Clone("hEmbedFit");
  TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2])-pow([3]/x/x,[4]))",1,20);
  func->SetParameters(1,0.1,1,0.01,1);
  hEmbedFit->Fit(func,"IR0");
  hEmbedFit->GetYaxis()->SetRangeUser(0,1.5);
  c = draw1D(hEmbedFit,"MTD response efficiency from embedding;p_{T} (GeV/c);Efficiency");
  func->SetLineColor(4);
  func->Draw("sames");
  TLegend *leg = new TLegend(0.3,0.2,0.5,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEmbedFit,"Embedding","PE");
  leg->AddEntry(func,"Fit function: [0]*exp{-(#frac{[1]}{x})^{[2]}-(#frac{[3]}{x^{2}})^{[4]}}","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitEmbedMtdRespEff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitEmbedMtdRespEff.png",run_type));
    }

  // compare
  c = draw1D(hEff[0],"MTD response efficiency;p_{T} (GeV/c);Efficiency");
  hEff[1]->Draw("sames P");

  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEff[0],"Cosmic ray","PL");
  leg->AddEntry(hEff[1],"Embedding","PL");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMtdRespEff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMtdRespEff.png",run_type));
    }

  // make ratio
  TH1F *hMtdRespCorr = (TH1F*)hEff[0]->Clone("MtdRespEffCorr");
  hMtdRespCorr->Divide(hEff[1]);
  c = draw1D(hMtdRespCorr,"Correction for MTD response efficiency: cosmic/embed;p_{T} (GeV/c);Cosmic/Embed");


  if(saveHisto)
    {
      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root","recreate");
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","recreate");
      
      for(int i=0; i<2; i++)
	{
	  hEff[i]->Write();
	}
      func->Write(Form("Fit_%s",hEff[1]->GetName()));
    }
}

//================================================
Double_t epsilonfit_1(Double_t *x, Double_t *par) 
{
  return (par[0] - par[1] * TMath::Erf(-x[0]+par[2]));
}




//================================================
void hadronEmbed()
{
  TFile *f = TFile::Open("output/Run14.AuAu200.Hadron.embed.root","read");
  TFile *f2 = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");

  const int nPart = 6;
  const int geantId[nPart] = {8,9,11,12,5,6};
  const char *partName[nPart] = {"Piplus","Piminus","Kplus","Kminus","Muplus","Muminus"};
  TH3F *hHadronPtVsPid[2][2];
  TH1F *hHadronPt[nCentBins][nPart][2];
  TH1F *hHadronEff[nCentBins][nPart];
  hHadronPtVsPid[0][0] = (TH3F*)f->Get("mhHadronMcPt");
  hHadronPtVsPid[0][1] = (TH3F*)f->Get("mhHadronRcPtMtd");
  hHadronPtVsPid[1][0] = (TH3F*)f2->Get("mhHadronMcPt");
  hHadronPtVsPid[1][0]->SetName("mhHadronMcPtJpsi");
  hHadronPtVsPid[1][1] = (TH3F*)f2->Get("mhHadronRcPtMtd");
  hHadronPtVsPid[1][1]->SetName("mhHadronRcPtMtdJpsi");
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hHadronPtVsPid[i][j]->Sumw2();
	}
    }

  TList *list = new TList;
  TString legName[nPart];
  for(int k=0; k<nCentBins; k++)
    {
      list->Clear();
      for(int i=0; i<nPart; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      if(i<4) hHadronPt[k][i][j] = (TH1F*)hHadronPtVsPid[0][j]->ProjectionX(Form("%s_%s_Cent%s",hHadronPtVsPid[0][j]->GetName(),partName[i],cent_Title[k]),geantId[i]+1,geantId[i]+1,centBins_low[k],centBins_high[k]);
	      else hHadronPt[k][i][j] = (TH1F*)hHadronPtVsPid[1][j]->ProjectionX(Form("%s_%s_Cent%s",hHadronPtVsPid[1][j]->GetName(),partName[i],cent_Title[k]),geantId[i]+1,geantId[i]+1,centBins_low[k],centBins_high[k]);
	    }
	  hHadronEff[k][i] = (TH1F*)hHadronPt[k][i][1]->Clone(Form("mhHadronMtdEff_%s_Cent%s",partName[i],cent_Title[k]));
	  hHadronEff[k][i]->Divide(hHadronPt[k][i][0]);
	  TH1F *htmp = (TH1F*)hHadronEff[k][i]->Clone(Form("%s_clone",hHadronEff[k][i]->GetName()));
	  list->Add(htmp);
	  legName[i] = partName[i];
	}
      c = drawHistos(list,Form("McTrkEff_cent%s",cent_Title[k]),Form("Efficiency of embedded tracks (%s%%);p_{T} (GeV/c);Efficiency",cent_Name[k]),kTRUE,0,13,kTRUE,0,0.8,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.3,0.45,0.6,0.85,kTRUE);
    }

  TString legName2[nCentBins];
  for(int i=0; i<nPart; i++)
    {
      list->Clear();
      for(int k=0; k<nCentBins; k++)
	{
	  legName[k] = cent_Name[k];
	  TH1F *htmp = (TH1F*)hHadronEff[k][i]->Clone(Form("%s_clone2",hHadronEff[k][i]->GetName()));
	  list->Add(htmp);
	}
      c = drawHistos(list,Form("McTrkEff_%s",partName[i]),Form("Efficiency of embedded %s;p_{T} (GeV/c);Efficiency",partName[i]),kTRUE,0,13,kTRUE,0,0.8,kFALSE,kTRUE,legName,kTRUE,"|#eta_{mc}|<0.5",0.3,0.45,0.6,0.85,kTRUE);
    }


  TFile *fout = TFile::Open("Rootfiles/Run14.AuAu200.HadronEff.root","recreate");
  for(int k=0; k<nCentBins; k++)
    {
      for(int i=0; i<nPart; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hHadronPt[k][i][j]->Write();
	    }
	  hHadronEff[k][i]->Write();
	}
    }
  fout->Close();
}


//================================================
void makeTrigEff(const int savePlot = 0, const int saveHisto = 0)
{
  TF1 *funcTrigEff[nCentBins];
  TF1 *funcTrigEffCorr[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      funcTrigEff[k] = new TF1(Form("MuonTrigEff_cent%s",cent_Title[k]),"[0]-exp(-1*[1]*(x-[2]))",0,12);
      if(k==0) funcTrigEff[k]->SetParameters(0.908, 2.09256, 1.03406e-14);      
      if(k==1) funcTrigEff[k]->SetParameters(0.905465, 2.05892, 9.27036e-15);
      if(k==2) funcTrigEff[k]->SetParameters(0.920748, 1.89851, 1.57097e-14);
      if(k==3) funcTrigEff[k]->SetParameters(0.883885, 2.1119, 3.33068e-15);

      funcTrigEffCorr[k] = new TF1(Form("MuonTrigEffCorr_cent%s",cent_Title[k]),"[0]-exp(-1*[1]*(x-[2]))",0,12);
      if(k==0) funcTrigEffCorr[k]->SetParameters(0.936584, 17.1834, 0.893653);      
      if(k==1) funcTrigEffCorr[k]->SetParameters(0.945798, 11.6541, 0.838267);
      if(k==2) funcTrigEffCorr[k]->SetParameters(0.935776, 5.75038, 0.499548);
      if(k==3) funcTrigEffCorr[k]->SetParameters(0.938118, 13.9462, 0.875862);
    }

  TCanvas *c = new TCanvas("MuonTrigEff","MuonTrigEff",1100,700);
  c->Divide(2,2);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
      gPad->SetGrid(1,1);
      ScaleHistoTitle(funcTrigEff[k]->GetHistogram(),0.06,1,0.05,0.06,0.9,0.05,62);
      funcTrigEff[k]->SetTitle(";p_{T} (GeV/c);Trigger efficiency");
      funcTrigEff[k]->SetMaximum(1.1);
      funcTrigEff[k]->SetMinimum(0.5);
      funcTrigEff[k]->SetLineColor(1);
      funcTrigEff[k]->Draw();
      funcTrigEffCorr[k]->SetLineColor(4);
      funcTrigEffCorr[k]->Draw("sames");
      TPaveText *t1 = GetPaveText(0.7,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.45,0.25,0.65,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(funcTrigEff[0],"VPDMB5","L");
  leg->AddEntry(funcTrigEffCorr[0],"NoVtx/VPDMB5","L");
  leg->Draw();
  
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins.png",run_type));
    }

  TH1F *hMuonTrigEff[nCentBins];
  TString legName[nCentBins];
  TList *list = new TList;
  for(int k=0; k<nCentBins; k++)
    {
      hMuonTrigEff[k] = new TH1F(Form("MuonTrigEff_cent%s",cent_Title[k]),"Muon efficiency",1000,0,10);
      for(int bin=1; bin<=hMuonTrigEff[k]->GetNbinsX(); bin++)
	{
	  double x = hMuonTrigEff[k]->GetXaxis()->GetBinCenter(bin);
	  hMuonTrigEff[k]->SetBinContent(bin,funcTrigEff[k]->Eval(x)*funcTrigEffCorr[k]->Eval(x));
	  hMuonTrigEff[k]->SetBinError(bin,0);
	}
      list->Add(hMuonTrigEff[k]);
      legName[k] = Form("%s%%",cent_Name[k]);
    }
  c = drawHistos(list,"MuonTrigEff_CentBins","Single muon trigger efficiency;p_{T} (GeV/c);Efficiency",kTRUE,0.5,10,kTRUE,0.5,1.0,kFALSE,kTRUE,legName,kTRUE,"",0.4,0.6,0.3,0.55,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins.png",run_type));
    }
  
  if(saveHisto)
    {
      char *outName = "Run14.AuAu200.MuonTrigEff.root";
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  funcTrigEff[k]->Write("",TObject::kOverwrite);
	  funcTrigEffCorr[k]->Write("",TObject::kOverwrite);
	}
      fout->Close();
    }
}


//================================================
void compTrigEff(const int savePlot = 0)
{
  TF1 *funcTrigEff[2];
  TF1 *funcTrigEffCorr[2];
  funcTrigEff[0] = new TF1(Form("MuonTrigEff_mid"),"[0]-exp(-1*[1]*(x-[2]))",0,12);
  funcTrigEff[0]->SetParameters(0.908, 2.09256, 1.03406e-14);      
  funcTrigEffCorr[0] = new TF1(Form("MuonTrigEffCorr_mid"),"[0]-exp(-1*[1]*(x-[2]))",0,12);
  funcTrigEffCorr[0]->SetParameters(0.936584, 17.1834, 0.893653); 

  funcTrigEff[1] = new TF1(Form("MuonTrigEff_low"),"[0]-exp(-1*[1]*(x-[2]))",0,12);
  funcTrigEff[1]->SetParameters(0.951456, 2.32223, 6.61693e-14);      
  funcTrigEffCorr[1] = new TF1(Form("MuonTrigEffCorr_low"),"[0]-exp(-1*[1]*(x-[2]))",0,12);
  funcTrigEffCorr[1]->SetParameters(0.96023, 10.5511, 0.803); 

  TH1F *hMuonTrigEff[2];
  TString legName[2] = {"Production_mid","Production_low"};
  TList *list = new TList;
  for(int k=0; k<2; k++)
    {
      hMuonTrigEff[k] = new TH1F(Form("MuonTrigEff_%d",k),"Single muon efficiency",1000,0,10);
      for(int bin=1; bin<=hMuonTrigEff[k]->GetNbinsX(); bin++)
	{
	  double x = hMuonTrigEff[k]->GetXaxis()->GetBinCenter(bin);
	  hMuonTrigEff[k]->SetBinContent(bin,funcTrigEff[k]->Eval(x)*funcTrigEffCorr[k]->Eval(x));
	  hMuonTrigEff[k]->SetBinError(bin,0);
	}
      list->Add(hMuonTrigEff[k]);
    }
  c = drawHistos(list,"MuonTrigEff_CentBins","Single muon trigger efficiency;p_{T} (GeV/c);Efficiency",kTRUE,0.5,10,kTRUE,0.5,1.0,kFALSE,kTRUE,legName,kTRUE,"",0.4,0.6,0.3,0.55,kTRUE);

  TH1F *hRatio = (TH1F*)hMuonTrigEff[1]->Clone("hRatio");
  hRatio->GetXaxis()->SetRangeUser(1,3);
  hRatio->Divide(hMuonTrigEff[0]);
  c = draw1D(hRatio,";p_{T,#mu} (GeV/c);Prod_low/Prod_mid");
  /*
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMuonTrigEffCombined.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMuonTrigEffCombined.png",run_type));
    }
  */
}

