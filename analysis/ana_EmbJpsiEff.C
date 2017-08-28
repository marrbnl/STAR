const int year = YEAR;
const char *trkEffType[6] = {"MC","Tpc","MtdMth","MuonPid","MtdTrig","TrigUnit"};
const char *weight_name[2] = {"","_w"};

//================================================
void ana_EmbJpsiEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //getJpsiWeight();
  embJpsiEff();

  //makeTrigEff();

  //ploEff();
  //plotEmbedEff();
  //hadronEmbed();
  //compTrigEff();
}



//================================================
void embJpsiEff(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");

  // Jpsi efficiency vs. pT
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** ptName       = pt_Name_pt;
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TList *list = new TList;
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }

  const int nHistos = 6;
  TH1F *hJpsiPt[nHistos][nCentBins];
  TH1F *hJpsiPtEff[nHistos][nCentBins];
  for(int i=0; i<nHistos; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiPt[i][k] = (TH1F*)fin->Get(Form("hJpsiPt_%s_cent%s",trkEffType[i],cent_Title[k]));
	  hJpsiPt[i][k]->Rebin(2);
	  int index = i-1;
	  if(i==0) index = 0;
	  hJpsiPtEff[i][k] = DivideTH1ForEff(hJpsiPt[i][k],hJpsiPt[index][k],Form("hJpsiPtEff_%s_cent%s",trkEffType[i],cent_Title[k]));
	}
    }

  // various efficiency
  const int kcent = 0;
  for(int i=1; i<nHistos; i++)
    {
      list->Add(hJpsiPtEff[i][kcent]);
    }
  TString legName2[5] = {"TPC tracking + p_{T,#mu} cut","MTD acceptance & response","Muon PID","MTD triggering","Trigger unit"};
  c = drawHistos(list,"JpsiEff_AllEffs",Form("%s: efficiencies for J/#psi ;p_{T} (GeV/c);Efficiency",run_type),true,0,15,true,1e-5,1.7,false,kTRUE,legName2,true,Form("%s%%",cent_Name[kcent]),0.2,0.4,0.63,0.88,kTRUE,0.04,0.035);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/JpsiEff_AllTypes.pdf",run_type));
  list->Clear();
  return;

  //==============================================
  // weight input spectrum shape
  //==============================================

  // check if the histograms are consistent
  TH2F *hJpsiPtMcVsRcAll = (TH2F*)fin->Get("Embed_JpsiPtMcVsRc_All");
  TH1F *hJpsiPtRcCheck[2];
  hJpsiPtRcCheck[0] = (TH1F*)hJpsiPtMcVsRcAll->ProjectionX("hJpsiPtRc_0");
  hJpsiPtRcCheck[1] = (TH1F*)hJpsiPt[5][0]->Clone("hJpsiPtRc_1");
  c = draw1D(hJpsiPtRcCheck[0]);
  hJpsiPtRcCheck[1]->SetLineColor(2);
  hJpsiPtRcCheck[1]->Draw("sames");

  // weight with Jpsi spectrum shape
  const int kNCent = nCentBins_npart[0];
  TH1F *hInputJpsi[kNCent];
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  for(int k=0; k<kNCent; k++)
    {
      if(k==0 || k==1) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
      if(k==2 || k==3) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
      if(k>=4)         hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");
      hInputJpsi[k]->SetName(Form("hInputJpsi_cent%s",cent_Name_npart[k]));
      hInputJpsi[k]->Scale(1./hInputJpsi[k]->Integral());
    }

  // weight with relative abundance
  TH1F *hJpsiCorr[gNTrgSetup][gNZdcRate];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorr[t][p] = (TH1F*)fin->Get(Form("Data_JpsiCountWeight_Zdc%d-%d%s",p*10,p*10+10,gTrgSetupName[t]));
	}
    }

  // matching between MC vs. RC
  TH2F *hJpsiPtMcVsRc[kNCent][gNZdcRate];
  TH1F *hJpsiPtMc[kNCent][gNZdcRate];
  TH1F *hJpsiPtRc[kNCent][gNZdcRate];
  TH1F *hJpsiEff[kNCent][gNZdcRate];
  TH1F *hJpsiPtMcWeight[kNCent][gNZdcRate];
  TH1F *hJpsiPtRcWeight[kNCent][gNZdcRate];
  TH1F *hJpsiEffWeight[kNCent][gNZdcRate];
  TH1F *hJpsiPtMcRebin[kNCent][gNZdcRate];
  TH1F *hJpsiPtRcRebin[kNCent][gNZdcRate];
  TH1F *hJpsiEffRebin[kNCent][gNZdcRate];
  for(int k=0; k<kNCent; k++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiPtMc[k][p]     = (TH1F*)fin->Get(Form("Embed_JpsiPtMc_cent%s_zdc%d-%d",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtMcVsRc[k][p] = (TH2F*)fin->Get(Form("Embed_JpsiPtMcVsRc_cent%s_zdc%d-%d",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtRc[k][p]     = (TH1F*)hJpsiPtMcVsRc[k][p]->ProjectionX(Form("Embed_JpsiPtRc_cent%s_zdc%d-%d",cent_Title_npart[k],p*10,p*10+10));

	  // Mc input
	  hJpsiPtMcWeight[k][p] = (TH1F*)hJpsiPtMc[k][p]->Clone(Form("%s_weight",hJpsiPtMc[k][p]->GetName()));
	  int nBinsX = hJpsiPtMcWeight[k][p]->GetNbinsX();
	  for(int xbin = 1; xbin<= nBinsX; xbin++)
	    {
	      double low_pt  = hJpsiPtMcWeight[k][p]->GetXaxis()->GetBinLowEdge(xbin);
	      double high_pt = hJpsiPtMcWeight[k][p]->GetXaxis()->GetBinUpEdge(xbin);
	      double jpsi_yield = hInputJpsi[k]->Integral(hInputJpsi[k]->FindFixBin(low_pt+1e-5),hInputJpsi[k]->FindFixBin(high_pt-1e-5));
	      hJpsiPtMcWeight[k][p]->SetBinContent(xbin, hJpsiPtMcWeight[k][p]->GetBinContent(xbin) * jpsi_yield);
	      hJpsiPtMcWeight[k][p]->SetBinError(xbin, hJpsiPtMcWeight[k][p]->GetBinError(xbin) * jpsi_yield);
	    }
	  hJpsiPtMcWeight[k][p]->Scale(hJpsiPtMcVsRc[k][p]->GetYaxis()->GetBinWidth(1)/hJpsiPtMcWeight[k][p]->GetBinWidth(1));
	  hJpsiPtMcRebin[k][p] = (TH1F*)hJpsiPtMcWeight[k][p]->Rebin(nbins, Form("%s_weight_rebin",hJpsiPtMc[k][p]->GetName()), xbins);

	  // Rc output
	  TH2F *h2Tmp = (TH2F*)hJpsiPtMcVsRc[k][p]->Clone(Form("%s_tmp",hJpsiPtMcVsRc[k][p]->GetName()));
	  int nBinsY = h2Tmp->GetNbinsY();
	  int nBinsX = h2Tmp->GetNbinsX();
	  for(int ybin = 1; ybin<= nBinsY; ybin++)
	    {
	      double low_pt = h2Tmp->GetYaxis()->GetBinLowEdge(ybin);
	      double high_pt = h2Tmp->GetYaxis()->GetBinUpEdge(ybin);
	      double jpsi_yield = hInputJpsi[k]->Integral(hInputJpsi[k]->FindFixBin(low_pt+1e-5),hInputJpsi[k]->FindFixBin(high_pt-1e-5));
	      for(int xbin = 1; xbin <= nBinsX; xbin++)
		{
		  h2Tmp->SetBinContent(xbin, ybin, h2Tmp->GetBinContent(xbin, ybin) * jpsi_yield);
		  h2Tmp->SetBinError(xbin, ybin, h2Tmp->GetBinError(xbin, ybin) * jpsi_yield);
		}
	    }
	  hJpsiPtRcWeight[k][p] = (TH1F*)h2Tmp->ProjectionX(Form("Embed_JpsiPtRc_cent%s_zdc%d-%d_weight",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtRcRebin[k][p] = (TH1F*)hJpsiPtRcWeight[k][p]->Rebin(nbins, Form("%s_rebin",hJpsiPtRcWeight[k][p]->GetName()), xbins);

	  // efficiencies
	  int rebin = int(hJpsiPtRc[k][p]->GetBinWidth(1)/hJpsiPtMc[k][p]->GetBinWidth(1));
	  hJpsiPtMc[k][p]->Rebin(rebin);
	  hJpsiEff[k][p] = (TH1F*)hJpsiPtRc[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiEff[k][p]->Divide(hJpsiPtMc[k][p]);

	  hJpsiPtMcWeight[k][p]->Rebin(rebin);
	  hJpsiEffWeight[k][p] = (TH1F*)hJpsiPtRcWeight[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d_weight",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiEffWeight[k][p]->Divide(hJpsiPtMcWeight[k][p]);

	  hJpsiEffRebin[k][p] = (TH1F*)hJpsiPtRcRebin[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d_weight_rebin",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiEffRebin[k][p]->Divide(hJpsiPtMcRebin[k][p]);
	}
    }

  TH1F *hJpsiEffVsPt[gNTrgSetup][nCentBins];
  TH1F *hJpsiEffVsPtWeight[gNTrgSetup][nCentBins];
  TH1F *hJpsiEffVsPtFinal[gNTrgSetup][nCentBins];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  // reconstructed 
	  hJpsiEffVsPt[t][k] = (TH1F*)hJpsiEff[0][0]->Clone(Form("JpsiEffVsPt_cent%s%s",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPt[t][k]->Reset();

	  hJpsiEffVsPtWeight[t][k] = (TH1F*)hJpsiEffWeight[0][0]->Clone(Form("JpsiEffVsPt_cent%s%s_weight",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPtWeight[t][k]->Reset();

	  hJpsiEffVsPtFinal[t][k] = (TH1F*)hJpsiEffRebin[0][0]->Clone(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPtFinal[t][k]->Reset();

	  int low_CentBin = 8 - centBins_high[k]/2;
	  int high_CentBin = 8 - (centBins_low[k]+1)/2;
	  double scaleAll = 0;
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      for(int centBin = low_CentBin; centBin <= high_CentBin; centBin++)
		{
		  double scale = hJpsiCorr[t][p]->GetBinContent(centBin+1);
		  hJpsiEffVsPt[t][k]->Add(hJpsiEff[centBin][p], scale);
		  hJpsiEffVsPtWeight[t][k]->Add(hJpsiEffWeight[centBin][p], scale);
		  hJpsiEffVsPtFinal[t][k]->Add(hJpsiEffRebin[centBin][p], scale);
		  scaleAll+= scale;
		}
	    }
	  hJpsiEffVsPt[t][k]->Scale(1./scaleAll);
	  hJpsiEffVsPtWeight[t][k]->Scale(1./scaleAll);
	  hJpsiEffVsPtFinal[t][k]->Scale(1./scaleAll);
	  printf("[i] Total scale = %4.2f\n",scaleAll);
	}
    }
  
  // check the effect of weighting & binning
  TH1F *h11[3];
  h11[0] = (TH1F*)hJpsiEffVsPt[0][0]->Clone("h11_clone_0");
  h11[1] = (TH1F*)hJpsiEffVsPtWeight[0][0]->Clone("h11_clone_1");
  h11[2] = (TH1F*)hJpsiEffVsPtFinal[0][0]->Clone("h11_clone_2");
  // h11[0] = (TH1F*)hJpsiEff[0][10]->Clone("h11_clone_0");
  // h11[1] = (TH1F*)hJpsiEffWeight[0][10]->Clone("h11_clone_1");
  // h11[2] = (TH1F*)hJpsiEffRebin[0][10]->Clone("h11_clone_2");
  for(int i=0; i<3; i++)
    {
      list->Add(h11[i]);
    }
  TString legName1[3] = {"Unweighted","Weighted","Weighted && rebinned"};
  c = drawHistos(list,"JpsiEff_Check_Weight_Rebin",Form("%s: combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type),true,0,15,true,5e-4,5e-2,true,kTRUE,legName1,true,Form("%s%%",cent_Name[0]),0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal_Check_Weight_Rebin.pdf",run_type));
  list->Clear();

  // centrality dependence
  for(int k=0; k<3; k++)
    {
      list->Add(hJpsiEffVsPtFinal[0][k]);
    }
  c = drawHistos(list,"JpsiEffFinal_vs_cent",Form("%s: combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type),true,0,15,true,5e-4,5e-2,true,kTRUE,legName_cent,true,"",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal_vs_Cent.pdf",run_type));
  list->Clear();

  return;

  TLegend *leg1 = new TLegend(0.2,0.15,0.4,0.35);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  TLegend *leg2 = new TLegend(0.55,0.15,0.7,0.35);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.04);
  for(int j=1; j<gNTrgSetup; j++)
    {
      for(int k=1; k<nCentBins; k++)
	{
	  hJpsiPtEffFinalRebin[j][k][1]->GetYaxis()->SetRangeUser(1e-5,0.1);
	  hJpsiPtEffFinalRebin[j][k][1]->SetMarkerStyle(19+k);
	  hJpsiPtEffFinalRebin[j][k][1]->SetMarkerColor(color[j-1]);
	  hJpsiPtEffFinalRebin[j][k][1]->SetLineColor(color[j-1]);
	  if(j==1 && k==1) c = draw1D(hJpsiPtEffFinalRebin[j][k][1],Form("Combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency"),true,true);
	  else hJpsiPtEffFinalRebin[j][k][1]->Draw("sames");
	  if(k==1) leg1->AddEntry(hJpsiPtEffFinalRebin[j][k][1],Form("AuAu_200%s",gTrgSetupTitle[j]),"P");
	  if(j==1) leg2->AddEntry(hJpsiPtEffFinalRebin[j][k][1],Form("%s%%",cent_Name[k]),"P");
	}
    }
  leg1->Draw();
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal.pdf",run_type));
  return;

  TString legName3[4] = {"Run14_AuAu200_production","Run14_AuAu200_production_low","Run14_AuAu200_production_mid","Run14_AuAu200_production_high"};
  for(int j=1; j<gNTrgSetup; j++)
    {
      TH1F *hEff = (TH1F*)hJpsiPtEffFinalRebin[j][1][1]->Clone(Form("hJpseEff_%d",j));
      list->Add(hEff);
    }
  c = drawHistos(list,"JpsiEff_CompLumi",Form("Efficiency for J/#psi;p_{T} (GeV/c);Efficiency",gTrgSetupTitle[jsetup]),false,2.0,3.8,true,0,0.025,false,kTRUE,legName3,true,"0-20%",0.15,0.35,0.65,0.88,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEff_CompLumi.pdf",run_type));
  list->Clear();

  return;





  // final average efficiency
  TFile *fYield = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiCounts[gNTrgSetup-1];
  double nJpsi[nCentBins-1][gNTrgSetup-1];
  double nJpsiCent[nCentBins-1];
  double nJpsiAll = 0;
  for(int k=0; k<nCentBins-1; k++)
    {
      nJpsiCent[k] = 0;
    }
  for(int i=0; i<gNTrgSetup-1; i++)
    {
      hJpsiCounts[i] = (TH1F*)fYield->Get(Form("NJpsiInCent_weight%s",gTrgSetupName[i+1]));
      for(int bin=1; bin<=hJpsiCounts[i]->GetNbinsX(); bin++)
	{
	  nJpsi[bin-1][i] = hJpsiCounts[i]->GetBinContent(bin);
	  nJpsiAll += hJpsiCounts[i]->GetBinContent(bin);
	  nJpsiCent[bin-1] += hJpsiCounts[i]->GetBinContent(bin);
	  printf("[i] %s%% %s: %4.2f Jpsi\n",cent_Title[bin],gTrgSetupTitle[i+1],nJpsi[bin-1][i]);
	}
    }

  TH1F *hJpsiPtEffAvg[nCentBins];
  TH1F *hJpsiOneOverEff[gNTrgSetup][nCentBins];
  for(int j=0; j<gNTrgSetup; j++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiOneOverEff[j][k] = (TH1F*)hJpsiPtEffFinalRebin[j][k][1]->Clone(Form("JpsiOneOverEff_cent%s%s",cent_Title[k],gTrgSetupName[j]));
	  hJpsiOneOverEff[j][k]->Reset();
	  for(int bin=1; bin<=hJpsiOneOverEff[j][k]->GetNbinsX(); bin++)
	    {
	      double value = hJpsiPtEffFinalRebin[j][k][1]->GetBinContent(bin);
	      double error = hJpsiPtEffFinalRebin[j][k][1]->GetBinError(bin);
	      hJpsiOneOverEff[j][k]->SetBinContent(bin, 1./value);
	      hJpsiOneOverEff[j][k]->SetBinError(bin, 1./value * error/value);
	    }
	}
    }


  for(int k=0; k<nCentBins; k++)
    {
      printf("+++++ %s%% +++++\n",cent_Name[k]);
      hJpsiPtEffAvg[k] = (TH1F*)hJpsiPtEffFinalRebin[0][0][1]->Clone(Form("JpsiPtEff_cent%s",cent_Title[k]));
      hJpsiPtEffAvg[k]->Reset();
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  if(k==0)
	    {
	      for(int icent=0; icent<nCentBins-1; icent++)
		{
		  hJpsiPtEffAvg[k]->Add(hJpsiOneOverEff[j+1][icent+1], nJpsi[icent][j]/nJpsiAll);
		}
	    }
	  else
	    {
	      hJpsiPtEffAvg[k]->Add(hJpsiOneOverEff[j+1][k], nJpsi[k-1][j]/nJpsiCent[k-1]);
	    }
	}
      for(int bin=1; bin<=hJpsiPtEffAvg[k]->GetNbinsX(); bin++)
	{
	  double value = hJpsiPtEffAvg[k]->GetBinContent(bin);
	  double error = hJpsiPtEffAvg[k]->GetBinError(bin);
	  hJpsiPtEffAvg[k]->SetBinContent(bin, 1./value);
	  hJpsiPtEffAvg[k]->SetBinError(bin, 1./value * error/value);
	  printf("[i] pt = %2.1f with eff = %4.2e, err = %4.2f%%\n",hJpsiPtEffAvg[k]->GetBinCenter(bin),hJpsiPtEffAvg[k]->GetBinContent(bin),error/value*100);
	}
      list->Add(hJpsiPtEffAvg[k]);
    }
  c = drawHistos(list,"JpsiEffCent",Form("Combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency"),false,2.0,3.8,true,1e-4,0.1,true,kTRUE,legName_cent,true,run_type,0.5,0.7,0.15,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal_Cent.pdf",run_type));
  list->Clear();

  /*
  // correct efficiency w.r.t. 0-20%
  TH1F *hJpsiPtEffAvgCorr[nCentBins];
  TH1F *hEffCorr[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiPtEffAvgCorr[k] = (TH1F*)hJpsiPtEffAvg[1]->Clone(Form("JpsiPtEff_cent%s_corr",cent_Title[k]));
      hEffCorr[k] = (TH1F*)fin->Get(Form("hJpsiEffCorr_cent%s",cent_Title[k]));
      hJpsiPtEffAvgCorr[k]->Multiply(hEffCorr[k]);
      list->Add(hJpsiPtEffAvgCorr[k]);
    }
  c = drawHistos(list,"JpsiEffCentCorr",Form("Combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency"),false,2.0,3.8,true,1e-4,0.1,true,kTRUE,legName_cent,true,run_type,0.5,0.7,0.15,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal_Cent_corr.pdf",run_type));
  list->Clear();
  */
 
  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiPtEffAvg[k]->Write("",TObject::kOverwrite);
	  //hJpsiPtEffAvgCorr[k]->Write("",TObject::kOverwrite);
	}
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hJpsiPtEffFinalRebin[j][1][1]->Write("",TObject::kOverwrite);
	  hJpsiPtEffFinal[j][1][1]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void getJpsiWeight(const int savePlot = 1, const int saveHisto = 0)
{
  const int* nCentBins      = nCentBins_npart; 
  const int* centBins_low   = centBins_low_npart;
  const int* centBins_high  = centBins_high_npart;
  const char** cent_Name    = cent_Name_npart;
  const char** cent_Title   = cent_Title_npart;
  const int kNCent          = nCentBins[0];

  // check the acceptance effect using mixed event
  TFile *fmix = TFile::Open(Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  const char* pName[2] = {"UL","LS"};
  TH1F *hMixInvMass[2][nCentBins_pt];
  TH1F *hAcc[nCentBins_pt];
  TList *list = new TList;
  TString legName[nCentBins_pt];
  for(int j=0; j<nCentBins_pt; j++)
    {
      for(int i=0; i<2; i++)
	{
	  hMixInvMass[i][j] = (TH1F*)fmix->Get(Form("Mix_InvMass_%s_pt0-15_cent%s",pName[i],cent_Name_pt[j]));
	  if(j==nCentBins_pt-1) hMixInvMass[i][j]->Rebin(5);
	}
      hAcc[j] = (TH1F*)hMixInvMass[0][j]->Clone(Form("Mix_InvMass_Acc_pt0-15_cent%s",cent_Name_pt[j]));
      hAcc[j]->Divide(hMixInvMass[1][j]);
      legName[j] = Form("%s%%",cent_Name_pt[j]);
      list->Add(hAcc[j]);
    }
  TCanvas *c = drawHistos(list,"Mix_Acceptance",Form("%s: acceptance effect from mixed event;M_{#mu#mu} (GeV/c^{2});UL/LS",run_type),true,1,5,true,0.85,1.2,false,kTRUE,legName,kTRUE,"p_{T} > 0 GeV/c",0.45,0.7,0.65,0.88,kTRUE,0.04);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/MixEvt_AccEff.pdf",run_type));
    }

  // Get the # of UL-LS pairs from data
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");

  TH1F *hJpsiInvMass[2][kNCent][gNTrgSetup][gNZdcRate];
  TH1F *hJpsiCount[gNTrgSetup][gNZdcRate];
  double min_mass[2] = {2.6, 3.2};
  double max_mass[2] = {2.9, 3.5};
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCount[t][p]  = new TH1F(Form("Data_JpsiCount_Zdc%d-%d%s",p*10,p*10+10,gTrgSetupName[t]),";;Count",kNCent,0,kNCent);
	  if(t==0 && (p==1 || p==4) )
	    {
	      c = new TCanvas(Form("ULvsLS_ZDC%d",p),Form("ULvsLS_ZDC%d-%d",p*10,p*10+10),1100,700);
	      c->Divide(4,2);
	    }
	  for(int j=0; j<kNCent; j++)
	    {
	      hJpsiCount[t][p]->GetXaxis()->SetBinLabel(j+1, Form("%s%%",cent_Name[j]));
	      hJpsiCount[t][p]->GetXaxis()->SetLabelSize(0.045);

	      for(int i=0; i<2; i++)
		{
		  hJpsiInvMass[i][j][t][p] = (TH1F*)fin->Get(Form("Data_InvMass_%s_cent%s_zdc%d-%d%s",pName[i],cent_Title_npart[j],p*10,p*10+10,gTrgSetupName[t]));
		  int nbins = hJpsiInvMass[i][j][t][p]->GetNbinsX();
		  for(int bin=1; bin<=nbins; bin++)
		    {
		      hJpsiInvMass[i][j][t][p]->SetBinError(bin, sqrt(hJpsiInvMass[i][j][t][p]->GetBinContent(bin)));
		    }
		}

	      // normalize in side-band region
	      double nCount[2] = {0,0};
	      for(int i=0; i<2; i++)
		{
		  for(int m=0; m<2; m++)
		    {
		      int low_bin = hJpsiInvMass[i][j][t][p]->FindFixBin(min_mass[m]);
		      int high_bin = hJpsiInvMass[i][j][t][p]->FindFixBin(max_mass[m]);
		      nCount[i] += hJpsiInvMass[i][j][t][p]->Integral(low_bin, high_bin);
		    }
		}
	      double scale = nCount[0]/nCount[1];
	      hJpsiInvMass[1][j][t][p]->Scale(scale);
	      low_bin = hJpsiInvMass[0][j][t][p]->FindFixBin(3.0);
	      high_bin = hJpsiInvMass[0][j][t][p]->FindFixBin(3.2);
	      double nUL, nULerr, nLS, nLSerr;
	      nUL = hJpsiInvMass[0][j][t][p]->Integral(low_bin,high_bin);
	      nLS = hJpsiInvMass[1][j][t][p]->Integral(low_bin,high_bin);
	      nULerr = sqrt(nUL);
	      nLSerr = sqrt(nLS);
	      if(nUL>0 && nUL >= nLS)
		{
		  hJpsiCount[t][p]->SetBinContent(j+1, nUL-nLS);
		  hJpsiCount[t][p]->SetBinError(j+1, sqrt(nULerr*nULerr+nLSerr*nLSerr));
		}
	      else
		{
		  printf("[w] Missing %s%% ZDC%d-%d%s\n",cent_Name[j],p*10,p*10+10,gTrgSetupName[t]);
		}

	      if(t==0 && (p==1 || p==4))
		{
		  c->cd(j+1);
		  hJpsiInvMass[0][j][t][p]->Rebin(5);
		  hJpsiInvMass[0][j][t][p]->GetXaxis()->SetRangeUser(2.5,4);
		  hJpsiInvMass[0][j][t][p]->SetMarkerStyle(21);
		  hJpsiInvMass[0][j][t][p]->SetMaximum(1.3*hJpsiInvMass[0][j][t][p]->GetMaximum());
		  hJpsiInvMass[0][j][t][p]->Draw("P");
		  hJpsiInvMass[1][j][t][p]->Rebin(5);
		  hJpsiInvMass[1][j][t][p]->SetLineColor(4);
		  hJpsiInvMass[1][j][t][p]->Draw("samesHIST");
		  TPaveText *t1 = GetTitleText(Form("%s%%, %d < ZDC < %d kHz",cent_Name[j],p*10,p*10+10),0.06);
		  t1->Draw();
		  if(j==0)
		    {
		      t1 = GetPaveText(0.2,0.5,0.8,0.85,0.06);
		      t1->SetTextColor(6);
		      t1->AddText(run_type);
		      t1->Draw();
		    }
		}
	    }
	  if(t==0 && (p==1 || p==4) )
	    {
	      if(savePlot)
		c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_ULvsLS_Zdc%d-%d.pdf",run_type,p*10,p*10+10));
	    }
	}
    }


  // plot nJpsi vs. centrality vs. ZDC
  double xmin = 0.15, xmax = 0.45, ymin = 0.65, ymax = 0.88;
  TLegend *leg[2];
  for(int l=0; l<2; l++)
    {
      leg[l] = new TLegend(xmin+0.3*l,ymin,xmax+0.3*l,ymax);
      leg[l]->SetBorderSize(0);
      leg[l]->SetFillColor(0);
      leg[l]->SetTextSize(0.03);
    }
  const int color[6] = {1, 2, 3, 4, 6, 7};
  TH1F *hMcTrkEff[gNTrgSetup][gNZdcRate];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCount[t][p]->GetXaxis()->SetLabelSize(0.05);
	  hJpsiCount[t][p]->SetMarkerStyle(20+p);
	  hJpsiCount[t][p]->SetMarkerColor(color[p%6]);
	  hJpsiCount[t][p]->SetLineColor(color[p%6]);
	  hJpsiCount[t][p]->GetYaxis()->SetRangeUser(1,1e6);
	}
    }

  for(int p=0; p<gNZdcRate; p++)
    {
      leg[p/6]->AddEntry(hJpsiCount[0][p],Form("%d < ZDC < %d kHz",p*10,p*10+10),"P");
      if(p==0)
	{
	  c = draw1D(hJpsiCount[0][p]);
	  gPad->SetLogy();
	  TPaveText *t1 = GetTitleText(Form("%s: # of UL-LS counts from data",run_type),0.04);
	  t1->Draw();
	}
      else     hJpsiCount[0][p]->Draw("sames");
    }
  for(int l=0; l<2; l++)
    leg[l]->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_NJpsi_US-LS.pdf",run_type));
    }

  // Apply TPC efficiency
  TH1F *hMcJpsiEff[gNZdcRate]; 
  TH1F *hJpsiCorr[gNTrgSetup][gNZdcRate];
  double nJpsiAll[gNTrgSetup];
  for(int t=0; t<gNTrgSetup; t++)
    {
      nJpsiAll[t] = 0;
      for(int p=0; p<gNZdcRate; p++)
	{
	  hMcJpsiEff[p] = (TH1F*)fin->Get(Form("McJpsiTpcEff_Zdc%d-%d",p*10,p*10+10));
	  hJpsiCorr[t][p] = (TH1F*)hJpsiCount[t][p]->Clone(Form("Data_JpsiCountWeight_Zdc%d-%d%s",p*10,p*10+10,gTrgSetupName[t]));
	  int nbins = hMcJpsiEff[p]->GetNbinsX();
	  for(int bin=1; bin<=nbins; bin++)
	    {
	      double nJpsi = hJpsiCorr[t][p]->GetBinContent(bin);
	      double eff = hMcJpsiEff[p]->GetBinContent(bin);
	      if(eff>0)
		{
		  hJpsiCorr[t][p]->SetBinContent(bin, nJpsi/eff);
		  hJpsiCorr[t][p]->SetBinError(bin, hJpsiCorr[t][p]->GetBinError(bin)/eff);
		}
	      else
		{
		  hJpsiCorr[t][p]->SetBinContent(bin, 0);
		  hJpsiCorr[t][p]->SetBinError(bin, 0);
		}
	      if(p==0)
		{
		  hJpsiCorr[t][p]->SetBinContent(bin, 0);
		  hJpsiCorr[t][p]->SetBinError(bin, 0);
		}
	      nJpsiAll[t] += hJpsiCorr[t][p]->GetBinContent(bin);
	    }
	}
      printf("[i] Total # of corrected Jpsi: %1.0f for %s\n",nJpsiAll[t],gTrgSetupName[t]);
      
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorr[t][p]->Scale(1./nJpsiAll[t]);
	  hJpsiCorr[t][p]->GetYaxis()->SetRangeUser(1e-4,10);
	}
    }

  for(int p=0; p<gNZdcRate; p++)
    {
      if(p==0) 
	{
	  c = draw1D(hJpsiCorr[0][p]);
	  gPad->SetLogy();
	  TPaveText *t1 = GetTitleText(Form("%s: corrected (UL-LS) counts as weights",run_type),0.04);
	  t1->Draw();
	}
      else     
	hJpsiCorr[0][p]->Draw("sames");
    }
  for(int l=0; l<2; l++)
    leg[l]->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_Jpsi_Weights.pdf",run_type));
    }

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      fin->cd();
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      hJpsiCorr[t][p]->Write("", TObject::kOverwrite);
	    }
	}
    }
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
void plotEmbedEff(const int savePlot = 0, const int saveHisto = 0)
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


