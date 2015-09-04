const double pt1_cut = 1.5, pt2_cut = 1.0;
const double low_mass = 2.9;
const double high_mass = 3.3;

const int year = 2014;

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


  embedding();
  //hadronEmbed();
}

//================================================
void embedding(const int savePlot = 0)
{
  TFile *fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
  TFile *fData = TFile::Open("output/Pico.Run14.AuAu200.jpsi.root","read");
  TH1F *hCent[2];
  hCent[0] = (TH1F*)fEmbed->Get("mhCentrality_di_mu");
  hCent[0]->SetName("mhCentrality_embed");
  hCent[1] = (TH1F*)fData->Get("mhCentrality_di_mu");
  for(int i=0; i<2; i++)
    {
      hCent[i]->Sumw2();
      hCent[i]->Scale(1./hCent[i]->GetBinContent(16));
      hCent[i]->SetMarkerStyle(20);
      hCent[i]->SetMarkerColor(2-i);
      hCent[i]->SetLineColor(2-i);
    }
  c = draw1D(hCent[0],"Centrality distribution",kTRUE,kTRUE);
  hCent[1]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.2,0.65,0.4,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hCent[0],"Embedding","P");
  leg->AddEntry(hCent[1],"Data","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.png",run_type));
    }

  TH1F *hgRefMult[2];
  hgRefMult[0] = (TH1F*)fEmbed->Get("mhgRefMult_di_mu");
  hgRefMult[0]->SetName("mhgRefMult_embed");
  hgRefMult[1] = (TH1F*)fData->Get("mhgRefMult_di_mu");
  for(int i=0; i<2; i++)
    {
      hgRefMult[i]->Sumw2();
      hgRefMult[i]->Scale(1./hgRefMult[i]->Integral());
      hgRefMult[i]->SetMarkerStyle(20);
      hgRefMult[i]->SetMarkerColor(2-i);
      hgRefMult[i]->SetLineColor(2-i);
    }
  draw1D(hgRefMult[0],"gRefMult distribution",kFALSE,kTRUE);
  hgRefMult[1]->Draw("samesHIST");
  TLegend *leg = new TLegend(0.6,0.65,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hgRefMult[0],"Embedding","P");
  leg->AddEntry(hgRefMult[1],"Data","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.png",run_type));
    }

  TH1F *hgRefMultCorr[2];
  hgRefMultCorr[0] = (TH1F*)fEmbed->Get("mhgRefMultCorr_di_mu");
  hgRefMultCorr[0]->SetName("mhgRefMultCorr_embed");
  hgRefMultCorr[1] = (TH1F*)fData->Get("mhgRefMultCorr_di_mu");
  for(int i=0; i<2; i++)
    {
      hgRefMultCorr[i]->Sumw2();
      hgRefMultCorr[i]->Scale(1./hgRefMultCorr[i]->Integral());
      hgRefMultCorr[i]->SetMarkerStyle(20);
      hgRefMultCorr[i]->SetMarkerColor(2-i);
      hgRefMultCorr[i]->SetLineColor(2-i);
    }
  draw1D(hgRefMultCorr[0],"Corrected gRefMult distribution",kFALSE,kTRUE);
  hgRefMultCorr[1]->Draw("samesHIST");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.png",run_type));
    }

  TH1F *hVertex[2];
  hVertex[0] = (TH1F*)fEmbed->Get("mhChosenVtxZ_di_mu");
  hVertex[1] = (TH1F*)fData->Get("mhTpcVz_di_mu");
  for(int i=0; i<2; i++)
    {
      hVertex[i]->Sumw2();
      hVertex[i]->Scale(1./hVertex[i]->Integral());
      hVertex[i]->SetMarkerStyle(20);
      hVertex[i]->SetMarkerColor(2-i);
      hVertex[i]->SetLineColor(2-i);
    }
  draw1D(hVertex[0],"Z distribution of TPC vertex",kFALSE,kTRUE);
  hVertex[1]->Draw("samesHIST");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.png",run_type));
    }

  /*
  // Jpsi efficiency
  TString name[3] = {"MCinput","TPCreco","MTDreco"};
  TString variable[4] = {"mass","pT","rapidity","phi"};
  const char *mc_name[4] = {"invariant mass","p_{T}","y","#varphi"};
  TH1F *hJpsi[3][4];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<4; j++)
	{
	  hJpsi[i][j] = (TH1F*)fin->Get(Form("%s_Jpsi_%s_%s",name[i].Data(),variable[j].Data(),cent_Title[0]));
	  hJpsi[i][j]->SetMinimum(0);
	  if(j==0) 
	    {
	      hJpsi[i][j]->GetXaxis()->SetRangeUser(2,4);
	      hJpsi[i][j]->SetXTitle("M_{#mu#mu} (GeV/c^{2})");
	    }
	  c = draw1D(hJpsi[i][j],Form("%s: %s distribution of J/psi",name[i].Data(),mc_name[j]),kFALSE,kFALSE);
	  if(savePlot)
	    {
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%s%s_Jpsi_%s_cent%s.pdf",run_type,run_cfg_name.Data(),name[i].Data(),variable[j].Data(),cent_Title[0]));
	      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Embed/%s%s_Jpsi_%s_cent%s.png",run_type,run_cfg_name.Data(),name[i].Data(),variable[j].Data(),cent_Title[0]));
	    }
	}
    }

  TString legName2[2] = {"Flat pT","Weighted pT"};
  TH1F *hJpsiPt[nCentBins][2];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiPt[k][0] = (TH1F*)fin->Get(Form("MCinput_Jpsi_pT_%s_FlatPt",cent_Title[k]));
      hJpsiPt[k][1] = (TH1F*)fin->Get(Form("MCinput_Jpsi_pT_%s_WeightPt",cent_Title[k]));
      list->Clear();
      for(int i=0; i<2; i++) list->Add(hJpsiPt[k][i]);
      char *centrality = Form(" (%s%%)",cent_Name[k]);
      if(year==2013) centrality = "";
      c = drawHistos(list,Form("Mcinput_Jpsi_pt_Cent%d",k),Form("p_{T} distribution of MC input J/psi %s",centrality),kTRUE,0,10,kTRUE,0.1,1e6,kTRUE,kTRUE,legName2,kTRUE,"",0.2,0.4,0.2,0.35,kTRUE);
    }
  */
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

