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
  //compare();

  //ploEff();
  //plotEmbedEff();
  //hadronEmbed();
  //compTrigEff();
}


//================================================
void compare(const int savePlot = 0)
{
  TFile *fin[2];
  fin[0] = TFile::Open(Form("Rootfiles/%s.MbEmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");
  fin[1] = TFile::Open(Form("Rootfiles/%s.EmbJpsiEff.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");

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
  const int kNCent = nCentBins_npart[0];

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TList *list = new TList;
  const int nHistos = 6;
  TH1F *hJpsiPt[2][nHistos][nCentBins];
  TH1F *hJpsiPtEffs[2][nHistos][nCentBins];
  for(int j=0; j<2; j++)
    {
      for(int i=0; i<nHistos; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      hJpsiPt[j][i][k] = (TH1F*)fin[j]->Get(Form("hJpsiPt_%s_cent%s",trkEffType[i],cent_Title[k]));
	      hJpsiPt[j][i][k]->SetName(Form("%s_file%d",hJpsiPt[j][i][k]->GetName(),j));
	      if(j==0) hJpsiPt[j][i][k]->Rebin(4);
	      if(j==1) hJpsiPt[j][i][k]->Rebin(4);
	      int index = i-1;
	      if(i==0) index = 0;
	      hJpsiPtEffs[j][i][k] = DivideTH1ForEff(hJpsiPt[j][i][k],hJpsiPt[j][index][k],Form("hJpsiPtEff_%s_cent%s_file%d",trkEffType[i],cent_Title[k],j));
	    }
	}
    }

  // various efficiency
  const int kcent = 0;
  for(int i=1; i<nHistos; i++)
    {
      hJpsiPtEffs[0][i][kcent]->Divide(hJpsiPtEffs[1][i][kcent]);
      list->Add(hJpsiPtEffs[0][i][kcent]);
    }
  TString legName2[5] = {"TPC tracking + p_{T,#mu} cut","MTD acceptance & response","Muon PID","MTD triggering","Trigger unit"};
  c = drawHistos(list,"JpsiEff_AllEffs",Form("%s: efficiencies for J/#psi ;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0,15,true,0.8,1.2,false,kTRUE,legName2,true,Form("%s%%",cent_Name[kcent]),0.2,0.4,0.63,0.88,kTRUE,0.04,0.035);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Compare_JpsiEff_AllTypes.pdf",run_type.Data()));
  list->Clear();
  //return;

  // weighted and rebinned
  TH1F *hJpsiEffAll[2][3];
  for(int j=0; j<2; j++)
    {
      hJpsiEffAll[j][0] = (TH1F*)fin[j]->Get("Embed_JpsiEff_All");
      hJpsiEffAll[j][1] = (TH1F*)fin[j]->Get("Embed_JpsiEff_All_weight");
      hJpsiEffAll[j][2] = (TH1F*)fin[j]->Get("Embed_JpsiEff_All_weight_rebin");
      for(int i=0; i<3; i++)
	{
	  hJpsiEffAll[j][i]->SetName(Form("%s_file%d",hJpsiEffAll[j][i]->GetName(),j));
	}
    }	  
  for(int i=2; i<3; i++)
    {
      hJpsiEffAll[0][i]->Divide(hJpsiEffAll[1][i]);
      list->Add(hJpsiEffAll[0][i]);
    }
  TString legName1[3] = {"Unweighted","Weighted","Weighted && rebinned"};
  c = drawHistos(list,"JpsiEff_Check_Weight_Rebin",Form("%s: total efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0,15,true,0, 2,false,kTRUE,legName1,true,"0-80%",0.5,0.7,0.2,0.45,kTRUE);
  list->Clear();
  

  // final efficiency
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }
  TH1F *hJpsiEffFinal[2][nCentBins][gNTrgSetup];
  for(int j=0; j<2; j++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  for(int t=0; t<gNTrgSetup; t++)
	    {
	      hJpsiEffFinal[j][k][t] = (TH1F*)fin[j]->Get(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[t]));
	      hJpsiEffFinal[j][k][t]->SetName(Form("%s_file%d",hJpsiEffFinal[j][k][t]->GetName(),j));
	    }
	}
    }
  for(int k=0; k<nCentBins; k++)
    {
      for(int t=0; t<gNTrgSetup; t++)
	{
	  hJpsiEffFinal[0][k][t]->Divide(hJpsiEffFinal[1][k][t]);
	}
      list->Add(hJpsiEffFinal[0][k][0]);
    }
  TCanvas *c = new TCanvas("JpsiEff_lumi", "JpsiEff_lumi", 1100, 700);
  c->Divide(2,2);
  leg = new TLegend(0.6,0.6,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.055);
  for(int k=1; k<nCentBins; k++)
    {
      c->cd(k);
      for(int t=1; t<gNTrgSetup; t++)
	{
	  hJpsiEffFinal[0][k][t]->SetMarkerStyle(20);
	  hJpsiEffFinal[0][k][t]->SetMarkerColor(gColor[t-1]);
	  hJpsiEffFinal[0][k][t]->SetLineColor(gColor[t-1]);
	  hJpsiEffFinal[0][k][t]->GetYaxis()->SetRangeUser(0.2, 2.2);
	  if(k==0) hJpsiEffFinal[0][k][t]->Draw();
	  else hJpsiEffFinal[0][k][t]->Draw("sames");
	  if(k==1) leg->AddEntry(hJpsiEffFinal[0][k][t], gTrgSetupTitle[t], "P");
	}
      TPaveText *t1 = GetTitleText(Form("Ratio of total J/psi efficiency: MB/DM (%s%%)",cent_Name[k]),0.055);
      t1->Draw();
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Compare_JpsiEffInLumi.pdf",run_type.Data()));

  c = drawHistos(list,"JpsiEff_Final",Form("%s: total efficiency for J/#psi ;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0,15,true,0,2,false,kTRUE,legName_cent,true,"",0.2,0.4,0.63,0.88,kTRUE,0.04,0.035);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Compare_JpsiEff_Final.pdf",run_type.Data()));
  
}



//================================================
void embJpsiEff(const int savePlot = 1, const int saveHisto = 1)
{
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"read");

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
  xbins[0] = 0.15;

  TList *list = new TList;
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }

  const int nHistos = 6;
  TH1F *hJpsiPt[nHistos][nCentBins];
  TH1F *hJpsiPtEffs[nHistos][nCentBins];
  for(int i=0; i<nHistos; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiPt[i][k] = (TH1F*)fin->Get(Form("hJpsiPt_%s_cent%s",trkEffType[i],cent_Title[k]));
	  hJpsiPt[i][k]->Rebin(4);
	  int index = i-1;
	  if(i==0) index = 0;
	  hJpsiPtEffs[i][k] = DivideTH1ForEff(hJpsiPt[i][k],hJpsiPt[index][k],Form("hJpsiPtEff_%s_cent%s",trkEffType[i],cent_Title[k]));
	}
    }

  // various efficiency
  const int kcent = 0;
  for(int i=1; i<nHistos; i++)
    {
      list->Add(hJpsiPtEffs[i][kcent]);
    }
  TString legName2[5] = {"TPC tracking + p_{T,#mu} cut","MTD acceptance & response","Muon PID","MTD triggering","Trigger unit"};
  c = drawHistos(list,"JpsiEff_AllEffs",Form("%s: efficiencies for J/#psi ;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0,15,true,1e-5,1.7,false,kTRUE,legName2,true,Form("%s%%",cent_Name[kcent]),0.2,0.4,0.63,0.88,kTRUE,0.04,0.035);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEff_AllTypes.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_JpsiEffVsPt_AllTypes.pdf"));
    }
  list->Clear();

  //==============================================
  // weight input spectrum shape
  //==============================================

  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");

  // check if the histograms are consistent
  TH2F *hJpsiPtMcVsRcAll = (TH2F*)fin->Get("Embed_JpsiPtMcVsRc_All");
  TH1F *hJpsiPtRcCheck[2];
  hJpsiPtRcCheck[0] = (TH1F*)hJpsiPtMcVsRcAll->ProjectionX("hJpsiPtRc_0");
  hJpsiPtRcCheck[0]->Rebin(4);
  hJpsiPtRcCheck[1] = (TH1F*)hJpsiPt[5][0]->Clone("hJpsiPtRc_1");
  c = draw1D(hJpsiPtRcCheck[0]);
  hJpsiPtRcCheck[1]->SetLineColor(2);
  hJpsiPtRcCheck[1]->Draw("sames");

  //==============================================
  // inclusive efficiency
  TH1F *hInputJpsiAll = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hJpsiPtMcAll = (TH1F*)fin->Get("Embed_JpsiPtMc_All");
  TH1F *hJpsiPtRcAll = (TH1F*)hJpsiPtMcVsRcAll->ProjectionX("Embed_JpsiPtRc_All");

  TH1F *hJpsiPtMcWeightAll = (TH1F*)hJpsiPtMcAll->Clone(Form("%s_weight",hJpsiPtMcAll->GetName()));
  int nBinsX = hJpsiPtMcWeightAll->GetNbinsX();
  for(int xbin = 1; xbin<= nBinsX; xbin++)
    {
      double low_pt  = hJpsiPtMcWeightAll->GetXaxis()->GetBinLowEdge(xbin);
      double high_pt = hJpsiPtMcWeightAll->GetXaxis()->GetBinUpEdge(xbin);
      double jpsi_yield = hInputJpsiAll->Integral(hInputJpsiAll->FindFixBin(low_pt+1e-5),hInputJpsiAll->FindFixBin(high_pt-1e-5));
      hJpsiPtMcWeightAll->SetBinContent(xbin, hJpsiPtMcWeightAll->GetBinContent(xbin) * jpsi_yield);
      hJpsiPtMcWeightAll->SetBinError(xbin, hJpsiPtMcWeightAll->GetBinError(xbin) * jpsi_yield);
    }
  hJpsiPtMcWeightAll->Scale(hJpsiPtMcVsRcAll->GetYaxis()->GetBinWidth(1)/hJpsiPtMcWeightAll->GetBinWidth(1));

  TH2F *h2Tmp = (TH2F*)hJpsiPtMcVsRcAll->Clone(Form("%s_tmp",hJpsiPtMcVsRcAll->GetName()));
  int nBinsY = h2Tmp->GetNbinsY();
  int nBinsX = h2Tmp->GetNbinsX();
  for(int ybin = 1; ybin<= nBinsY; ybin++)
    {
      double low_pt = h2Tmp->GetYaxis()->GetBinLowEdge(ybin);
      double high_pt = h2Tmp->GetYaxis()->GetBinUpEdge(ybin);
      double jpsi_yield = hInputJpsiAll->Integral(hInputJpsiAll->FindFixBin(low_pt+1e-5),hInputJpsiAll->FindFixBin(high_pt-1e-5));
      for(int xbin = 1; xbin <= nBinsX; xbin++)
	{
	  h2Tmp->SetBinContent(xbin, ybin, h2Tmp->GetBinContent(xbin, ybin) * jpsi_yield);
	  h2Tmp->SetBinError(xbin, ybin, h2Tmp->GetBinError(xbin, ybin) * jpsi_yield);
	}
    }
  TH1F *hJpsiPtRcWeightAll = (TH1F*)h2Tmp->ProjectionX("Embed_JpsiPtRcAll");

  // efficiencies
  int rebin = int(hJpsiPtRcAll->GetBinWidth(1)/hJpsiPtMcAll->GetBinWidth(1));
  hJpsiPtMcAll->Rebin(rebin);
  TH1F *hJpsiPtEffAll = (TH1F*)hJpsiPtRcAll->Clone(Form("Embed_JpsiEff_All"));
  hJpsiPtEffAll->Divide(hJpsiPtMcAll);

  hJpsiPtMcWeightAll->Rebin(rebin);
  TH1F *hJpsiPtEffWeightAll = (TH1F*)hJpsiPtRcWeightAll->Clone(Form("Embed_JpsiEff_All_weight"));
  hJpsiPtEffWeightAll->Divide(hJpsiPtMcWeightAll);

  TH1F *hJpsiPtMcRebinAll = (TH1F*)hJpsiPtMcWeightAll->Rebin(nbins, Form("%s_weight_rebin",hJpsiPtMcAll->GetName()), xbins);
  TH1F *hJpsiPtRcRebinAll = (TH1F*)hJpsiPtRcWeightAll->Rebin(nbins, Form("%s_rebin",hJpsiPtRcWeightAll->GetName()), xbins);
  hJpsiPtEffRebinAll = (TH1F*)hJpsiPtRcRebinAll->Clone(Form("Embed_JpsiEff_All_weight_rebin"));
  hJpsiPtEffRebinAll->Divide(hJpsiPtMcRebinAll);

  // check the effect of weighting & binning
  TH1F *h11[3];
  h11[0] = (TH1F*)hJpsiPtEffAll->Clone("h11_clone_0");
  h11[1] = (TH1F*)hJpsiPtEffWeightAll->Clone("h11_clone_1");
  h11[2] = (TH1F*)hJpsiPtEffRebinAll->Clone("h11_clone_2");
  for(int i=0; i<3; i++)
    {
      list->Add(h11[i]);
    }
  TString legName1[3] = {"Unweighted","Weighted","Weighted && rebinned"};
  c = drawHistos(list,"JpsiEff_Check_Weight_Rebin",Form("%s: total efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0,15,true,5e-4,5e-2,true,kTRUE,legName1,true,"0-80%",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEff_Check_Weight_Rebin.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_CheckWeightRebin.pdf"));
    }
  list->Clear();
  //return;

  //==============================================
  // individual efficiency
  TH1F *hInputJpsi[nCentBinsEff];
  for(int k=0; k<nCentBinsEff; k++)
    {
      if(nCentBinsEff==8)
	{
	  if(k==0 || k==1) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
	  if(k==2 || k==3) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
	  if(k>=4)         hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");
	}
      if(nCentBinsEff==16)
	{
	  if(k>=0 && k<=3) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0020");
	  if(k>=4 && k<=7) hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent2040");
	  if(k>=8)         hInputJpsi[k] = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent4060");
	}
      hInputJpsi[k]->SetName(Form("hInputJpsi_cent%s",centNameEff[k]));
      hInputJpsi[k]->Scale(1./hInputJpsi[k]->Integral());
    }

  // matching between MC vs. RC
  TH2F *hJpsiPtMcVsRc[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtMc[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtRc[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtEff[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtMcWeight[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtRcWeight[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtEffWeight[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtMcRebin[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtRcRebin[nCentBinsEff][gNZdcRate];
  TH1F *hJpsiPtEffRebin[nCentBinsEff][gNZdcRate];
  for(int k=0; k<nCentBinsEff; k++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiPtMc[k][p]     = (TH1F*)fin->Get(Form("Embed_JpsiPtMc_cent%s_zdc%d-%d",centTitleEff[k],p*10,p*10+10));
	  hJpsiPtMcVsRc[k][p] = (TH2F*)fin->Get(Form("Embed_JpsiPtMcVsRc_cent%s_zdc%d-%d",centTitleEff[k],p*10,p*10+10));
	  hJpsiPtRc[k][p]     = (TH1F*)hJpsiPtMcVsRc[k][p]->ProjectionX(Form("Embed_JpsiPtRc_cent%s_zdc%d-%d",centTitleEff[k],p*10,p*10+10));

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
	  hJpsiPtEff[k][p] = (TH1F*)hJpsiPtRc[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtEff[k][p]->Divide(hJpsiPtMc[k][p]);

	  hJpsiPtMcWeight[k][p]->Rebin(rebin);
	  hJpsiPtEffWeight[k][p] = (TH1F*)hJpsiPtRcWeight[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d_weight",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtEffWeight[k][p]->Divide(hJpsiPtMcWeight[k][p]);

	  hJpsiPtEffRebin[k][p] = (TH1F*)hJpsiPtRcRebin[k][p]->Clone(Form("Embed_JpsiEff_cent%s_zdc%d-%d_weight_rebin",cent_Title_npart[k],p*10,p*10+10));
	  hJpsiPtEffRebin[k][p]->Divide(hJpsiPtMcRebin[k][p]);
	}
    }

  //==============================================
  // Apply efficiency scale factors w.r.t. inclusive efficiency
  //==============================================
  TH1F *hJpsiTpcEffAll = (TH1F*)fin->Get("JpsiTpcEff_All_rebin");
  TH1F *hJpsiTpcEff[nCentBinsEff][gNZdcRate];
  for(int k=0; k<nCentBinsEff; k++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiTpcEff[k][p] = (TH1F*)fin->Get(Form("JpsiTpcEff_cent%s_Zdc%d-%d_rebin",centTitleEff[k],p*10,p*10+10));
	}
    }

  TH1F *hJpsiPtEffCorrRebin[nCentBinsEff][gNZdcRate];
  for(int k=0; k<nCentBinsEff; k++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiPtEffCorrRebin[k][p] = (TH1F*)hJpsiPtEffRebinAll->Clone(Form("Embed_JpsiEffCorr_cent%s_zdc%d-%d_weight_rebin",centTitleEff[k],p*10,p*10+10));
	  int nXbins = hJpsiPtEffCorrRebin[k][p]->GetNbinsX();
	  for(int bin=1; bin<=nXbins; bin++)
	    {
	      double corr =  hJpsiTpcEff[k][p]->GetBinContent(bin)/hJpsiTpcEffAll->GetBinContent(bin);
	      hJpsiPtEffCorrRebin[k][p]->SetBinContent(bin, hJpsiPtEffCorrRebin[k][p]->GetBinContent(bin) * corr);
	      hJpsiPtEffCorrRebin[k][p]->SetBinError(bin, hJpsiPtEffCorrRebin[k][p]->GetBinError(bin) * corr);
	    }
	}
    }

  TCanvas *cEffCheck[nCentBinsEff];
  for(int k=0; k<nCentBinsEff; k++)
    {
      if(k%2==1) continue;
      cEffCheck[k] = new TCanvas(Form("CheckEffScale_%d",k),Form("CheckEffScale_%s",centTitleEff[k]),1100,700);
      cEffCheck[k]->Divide(4,3);
      for(int p=0; p<gNZdcRate; p++)
	{
	  cEffCheck[k]->cd(p+1);
	  gPad->SetLogy();
	  hJpsiPtEffRebin[k][p]->GetYaxis()->SetRangeUser(1e-3,0.1);
	  hJpsiPtEffRebin[k][p]->SetMarkerStyle(24);
	  hJpsiPtEffRebin[k][p]->Draw();
	  hJpsiPtEffCorrRebin[k][p]->SetMarkerStyle(20);
	  hJpsiPtEffCorrRebin[k][p]->SetMarkerColor(2);
	  hJpsiPtEffCorrRebin[k][p]->SetLineColor(2);
	  hJpsiPtEffCorrRebin[k][p]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%d<ZDCx<%d kHz",p*10,p*10+10),0.07);
	  t1->Draw();
	}
      cEffCheck[k]->cd(12);
      TPaveText *t1 = GetPaveText(0.2,0.6,0.4,0.8,0.07);
      t1->AddText(run_type.Data());
      t1->AddText(Form("%s%%",centNameEff[k]));
      t1->Draw();
      if(savePlot) cEffCheck[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEff_CheckCorr_%s.pdf",run_type.Data(),gEmbTypeName[iMbEmb],centTitleEff[k]));
      if(gSaveAN)
	{
	  cEffCheck[k]->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_CheckCorr_cent%s.pdf",centTitleEff[k]));
	}
    }

  //==============================================
  // weight with relative abundance from data
  //==============================================
  // function to fit single track efficiency
  TFile *fTrkEff = TFile::Open(Form("Rootfiles/%s.%sEmbTrkEff.root",run_type.Data(),gEmbTypeName[iMbEmb]),"read");
  TF1 *funcTrkEff[nCentBinsEff][gNZdcRate];
  for(int k=0; k<nCentBinsEff; k++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  funcTrkEff[k][p] = (TF1*)fTrkEff->Get(Form("func_TpcTrkEff_cent%s_Zdc%d-%d",centTitleEff[k],p*10,p*10+10));
	}
    }

  TH1F *hJpsiCorr[gNTrgSetup][gNZdcRate];  // 0 - fit US-LS; 1 - count US-LS; 2 - LS; 3 - MB events
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorr[t][p] = (TH1F*)fin->Get(Form("Data_JpsiCountWeight_Zdc%d-%d%s_M3",p*10,p*10+10,gTrgSetupName[t]));
	}
    }

  TH1F *hJpsiEffVsPtUnCorr[gNTrgSetup][nCentBins];
  TH1F *hJpsiEffVsPtFinal[gNTrgSetup][nCentBins];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiEffVsPtFinal[t][k] = (TH1F*)hJpsiPtEffCorrRebin[0][0]->Clone(Form("JpsiEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPtFinal[t][k]->Reset();

	  hJpsiEffVsPtUnCorr[t][k] = (TH1F*)hJpsiPtEffRebin[0][0]->Clone(Form("JpsiEffVsPt_cent%s%s_Uncorr",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPtUnCorr[t][k]->Reset();
	    
	  int low_CentBin = 0, high_CentBin = 0;
	  if(nCentBinsEff==8)
	    {
	      low_CentBin = 8 - centBins_high[k]/2;
	      high_CentBin = 8 - (centBins_low[k]+1)/2;
	    }
	  if(nCentBinsEff==16)
	    {
	      low_CentBin = 16 - centBins_high[k];
	      high_CentBin = 16 - centBins_low[k];
	    }
	  double scaleAll = 0;
	  double trkEffErr = 0; // error due to TPC tracking efficiency
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      for(int centBin = low_CentBin; centBin <= high_CentBin; centBin++)
		{
		  double scale = hJpsiCorr[t][p]->GetBinContent(centBin+1);
		  hJpsiEffVsPtFinal[t][k]->Add(hJpsiPtEffCorrRebin[centBin][p], scale);
		  hJpsiEffVsPtUnCorr[t][k]->Add(hJpsiPtEffRebin[centBin][p], scale);
		  scaleAll+= scale;
		  if(scale>0)
		    {
		      trkEffErr += pow(scale*funcTrkEff[centBin][p]->GetParError(0)/funcTrkEff[centBin][p]->GetParameter(0)*2*hJpsiPtEffCorrRebin[centBin][p]->GetBinContent(1), 2);
		    }
		}
	    }
	  hJpsiEffVsPtFinal[t][k]->Scale(1./scaleAll);
	  hJpsiEffVsPtUnCorr[t][k]->Scale(1./scaleAll);
	  trkEffErr = sqrt(trkEffErr)/scaleAll;
	  printf("[i] %s, %s: error from trk eff = %4.3e +/- %4.3e (%4.3f%%)\n",cent_Title[k],gTrgSetupTitle[t],hJpsiEffVsPtFinal[t][k]->GetBinContent(1),trkEffErr,trkEffErr/hJpsiEffVsPtFinal[t][k]->GetBinContent(1)*100);

	  // calculate error on efficiency due to weights
	  for(int bin=1; bin<=hJpsiEffVsPtFinal[t][k]->GetNbinsX(); bin++)
	    {
	      double eff = hJpsiEffVsPtFinal[t][k]->GetBinContent(bin);
	      double eff_err = 0;
	      for(int p=0; p<gNZdcRate; p++)
		{
		  for(int centBin = low_CentBin; centBin <= high_CentBin; centBin++)
		    {
		      double eff_i = hJpsiPtEffCorrRebin[centBin][p]->GetBinContent(bin);
		      double w_i     = hJpsiCorr[t][p]->GetBinContent(centBin+1);
		      double w_i_err = hJpsiCorr[t][p]->GetBinError(centBin+1);
		      eff_err += pow( (eff_i-eff)/scaleAll*w_i_err , 2);
		    }
		}
	      eff_err = sqrt(eff_err);
	      if(bin==1)
		{
		  //printf("[i] pt bin %d: eff = %4.3e +/- %4.3e (%4.3f%%)\n",bin,eff,eff_err,eff_err/eff*100);
		}
	    }
	}
    }

  // Like-sign
  TH1F *hJpsiCorrLS[gNTrgSetup][gNZdcRate];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorrLS[t][p] = (TH1F*)fin->Get(Form("Data_JpsiCountWeight_Zdc%d-%d%s_M2",p*10,p*10+10,gTrgSetupName[t]));
	}
    }
  TH1F *hJpsiEffVsPtFinalLS[gNTrgSetup][nCentBins];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiEffVsPtFinalLS[t][k] = (TH1F*)hJpsiPtEffCorrRebin[0][0]->Clone(Form("LSEffVsPt_cent%s%s_final",cent_Title[k],gTrgSetupTitle[t]));
	  hJpsiEffVsPtFinalLS[t][k]->Reset();

	  double scaleAll = 0;
	  int low_CentBin = 0, high_CentBin = 0;
	  if(nCentBinsEff==8)
	    {
	      low_CentBin = 8 - centBins_high[k]/2;
	      high_CentBin = 8 - (centBins_low[k]+1)/2;
	    }
	  if(nCentBinsEff==16)
	    {
	      low_CentBin = 16 - centBins_high[k];
	      high_CentBin = 16 - centBins_low[k];
	    }
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      for(int centBin = low_CentBin; centBin <= high_CentBin; centBin++)
		{
		  double scale = hJpsiCorrLS[t][p]->GetBinContent(centBin+1);
		  hJpsiEffVsPtFinalLS[t][k]->Add(hJpsiPtEffCorrRebin[centBin][p], scale);
		  scaleAll+= scale;
		}
	    }
	  hJpsiEffVsPtFinalLS[t][k]->Scale(1./scaleAll);
	}
    }

  // demonstrate low statistics
  for(int k=0; k<nCentBins; k++)
    {
      list->Add(hJpsiEffVsPtUnCorr[0][k]);
    }
  c = drawHistos(list,"JpsiEffFinal_vs_pt_InCentBin_uncorr",Form("%s: combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0.5,15,true,5e-4,5e-2,true,kTRUE,legName_cent,true,"",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEffFinal_vs_Pt_Uncorr.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  list->Clear();

  // luminosity dependence
  TString legName_lumi[gNTrgSetup];
  for(int t=0; t<gNTrgSetup; t++)
    {
      legName_lumi[t] = Form("%s%s",run_type.Data(),gTrgSetupTitle[t]);
      list->Add(hJpsiEffVsPtFinal[t][0]);
    }
  c = drawHistos(list,"JpsiEffFinal_vs_pt_InLumiBin",Form("%s: combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0.15,15,true,5e-4,5e-2,true,kTRUE,legName_lumi,true,"",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEffFinal_vs_Pt_Lumi.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  list->Clear();

  TH1F *hJpsiEffInLumiRatio[gNTrgSetup];
  for(int t=0; t<gNTrgSetup; t++)
    {
      hJpsiEffInLumiRatio[t] = (TH1F*)hJpsiEffVsPtFinal[t][0]->Clone(Form("%s_ratio",hJpsiEffVsPtFinal[t][0]->GetName()));
      hJpsiEffInLumiRatio[t]->Divide(hJpsiEffVsPtFinal[0][0]);
      list->Add(hJpsiEffInLumiRatio[t]);
    }
  c = drawHistos(list,"JpsiEffFinal_vs_pt_InLumiBin_ratio",Form("%s: ratio of J/#psi to that of inclusive;p_{T} (GeV/c);Ratio",run_type.Data()),true,0.15,15,true,0.5,1.5,false,kTRUE,legName_lumi,true,"",0.5,0.7,0.15,0.35,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEffFinal_vs_Pt_Lumi_ratio.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  list->Clear();

  // centrality dependence
  for(int k=0; k<nCentBins; k++)
    {
      list->Add(hJpsiEffVsPtFinal[0][k]);
    }
  c = drawHistos(list,"JpsiEffFinal_vs_pt_InCentBin",Form("%s: combined efficiency for J/#psi;p_{T} (GeV/c);Efficiency",run_type.Data()),true,0.15,15,true,5e-4,5e-2,true,kTRUE,legName_cent,true,"",0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEffFinal_vs_Pt.pdf",run_type.Data(),gEmbTypeName[iMbEmb]));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_JpsiEffFinalVsPt.pdf"));
    }
  list->Clear();


  //==============================================
  // npart analysis
  //==============================================
  TH1F *hJpsiCentEffAll[nPtBins_npart];
  for(int i=0; i<nPtBins_npart; i++)
    {
      hJpsiCentEffAll[i] = new TH1F(Form("Embed_JpsiCentEff_All_Pt%1.0f",ptBins_low_npart[i]), "", 1, 0, 1);
      int low_bin = hJpsiPtMcWeightAll->FindFixBin(ptBins_low_npart[i]+1e-4);
      int up_bin =  hJpsiPtMcWeightAll->FindFixBin(ptBins_high_npart[i]-1e-4);
      double nMcJpsi_err, nRcJpsi_err;
      double nMcJpsi = hJpsiPtMcWeightAll->IntegralAndError(low_bin, up_bin, nMcJpsi_err);
      double nRcJpsi = hJpsiPtRcWeightAll->IntegralAndError(low_bin, up_bin, nRcJpsi_err);
      double eff = nRcJpsi/nMcJpsi;
      double err = eff * TMath::Sqrt(TMath::Power(nMcJpsi_err/nMcJpsi, 2) + TMath::Power(nRcJpsi_err/nRcJpsi, 2));
      hJpsiCentEffAll[i]->SetBinContent(1, eff);
      hJpsiCentEffAll[i]->SetBinError(1, err);
    }
  TH1F *hJpsiTpcEffVsCentAll[nPtBins_npart];
  TH1F *hJpsiTpcEffVsCent[nPtBins_npart][gNZdcRate];
  for(int i=0; i<nPtBins_npart; i++)
    {
      hJpsiTpcEffVsCentAll[i] = (TH1F*)fin->Get(Form("JpsiTpcEffVsCent_Pt%1.0f_All",ptBins_low_npart[i]));
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiTpcEffVsCent[i][p] = (TH1F*)fin->Get(Form("JpsiTpcEffVsCent_Pt%1.0f_Zdc%d-%d",ptBins_low_npart[i],p*10,p*10+10));
	}
    }
  TH1F *hJpsiEffVsCentCorr[nPtBins_npart][gNZdcRate];
  for(int i=0; i<nPtBins_npart; i++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiEffVsCentCorr[i][p] = (TH1F*)hJpsiTpcEffVsCent[i][p]->Clone(Form("Embed_JpsiEffVsCentCorr_Pt%1.0f_zdc%d-%d",ptBins_low_npart[i],p*10,p*10+10));
	  hJpsiEffVsCentCorr[i][p]->Reset();
	  int nbinsX = hJpsiEffVsCentCorr[i][p]->GetNbinsX();
	  for(int bin=1; bin<=nbinsX; bin++)
	    {
	      double corr = hJpsiTpcEffVsCent[i][p]->GetBinContent(bin)/hJpsiTpcEffVsCentAll[i]->GetBinContent(1);
	      hJpsiEffVsCentCorr[i][p]->SetBinContent(bin, hJpsiCentEffAll[i]->GetBinContent(1) * corr);
	      hJpsiEffVsCentCorr[i][p]->SetBinError(bin, hJpsiCentEffAll[i]->GetBinError(1) * corr);
	    }
	}
    }

  TH1F *hJpsiEffVsCentFinal[gNTrgSetup][nPtBins_npart];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int i=0; i<nPtBins_npart; i++)
	{
	  hJpsiEffVsCentFinal[t][i] = (TH1F*)hJpsiEffVsCentCorr[i][0]->Clone(Form("JpsiEffVsCent_Pt%1.0f%s_final",ptBins_low_npart[i],gTrgSetupTitle[t]));
	  hJpsiEffVsCentFinal[t][i]->Reset();

	  int nbinsX = hJpsiEffVsCentFinal[t][i]->GetNbinsX();
	  for(int bin=1; bin<=nbinsX; bin++)
	    {
	      double scaleAll = 0;
	      double eff = 0;
	      double err = 0;
	      int low_CentBin = 0, high_CentBin = 0;
	      if(nCentBinsEff==8)
		{
		  low_CentBin = 8 - centBins_high_npart[bin+i*8-1]/2;
		  high_CentBin = 8 - (centBins_low_npart[bin+i*8-1]+1)/2;
		}
	      if(nCentBinsEff==16)
		{
		  low_CentBin = 16 - centBins_high_npart[bin+i*8-1];
		  high_CentBin = 16 - centBins_low_npart[bin+i*8-1];
		}
	      for(int p=0; p<gNZdcRate; p++)
		{
		  double scale = 0;
		  for(int ibin=low_CentBin; ibin<=high_CentBin; ibin++)
		    {
		      scale += hJpsiCorr[t][p]->GetBinContent(ibin+1);
		    }
		  scaleAll+= scale;
		  eff += hJpsiEffVsCentCorr[i][p]->GetBinContent(bin)*scale;
		  err += TMath::Power(hJpsiEffVsCentCorr[i][p]->GetBinError(bin)*scale, 2);
		}
	      hJpsiEffVsCentFinal[t][i]->SetBinContent(bin, eff/scaleAll);
	      hJpsiEffVsCentFinal[t][i]->SetBinError(bin, sqrt(err)/scaleAll);   
	      hJpsiEffVsCentFinal[t][i]->GetXaxis()->SetBinLabel(bin, Form("%s%%",cent_Name_npart[bin+i*8-1]));
	    }
	}
    }

  for(int i=0; i<nPtBins_npart; i++)
    {
      hJpsiEffVsCentFinal[0][i]->SetMarkerStyle(20);
      hJpsiEffVsCentFinal[0][i]->SetMarkerSize(1.5);
      hJpsiEffVsCentFinal[0][i]->GetXaxis()->SetLabelSize(0.05);

      if(i==0) hJpsiEffVsCentFinal[0][i]->GetYaxis()->SetRangeUser(0, 0.01);
      else     hJpsiEffVsCentFinal[0][i]->GetYaxis()->SetRangeUser(0, 0.02);
      c = draw1D(hJpsiEffVsCentFinal[0][i], Form("%s: combined efficiency for J/#psi (%s);;Efficiency",run_type.Data(),pt_Title_npart[i]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/%sJpsiEffFinal_vs_Cent_Pt%1.0f.pdf",run_type.Data(),gEmbTypeName[iMbEmb],ptBins_low_npart[i]));
      if(gSaveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_JpsiEffFinalVsCent_Pt%1.0f.pdf",ptBins_low_npart[i]));
	}
    }
      
  if(saveHisto)
    {
      fin->cd();
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      hJpsiEffVsPtFinal[t][k]->Write("",TObject::kOverwrite);
	      hJpsiEffVsPtFinalLS[t][k]->Write("",TObject::kOverwrite);
	    }
	}

      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int i=0; i<nPtBins_npart; i++)
	    {
	      hJpsiEffVsCentFinal[t][i]->Write("",TObject::kOverwrite);
	    }
	}
      hJpsiPtEffAll->Write("",TObject::kOverwrite);
      hJpsiPtEffWeightAll->Write("",TObject::kOverwrite);
      hJpsiPtEffRebinAll->Write("",TObject::kOverwrite);
    }
}

//================================================
void getJpsiWeight(const int savePlot = 1, const int saveHisto = 1)
{
  // check the acceptance effect using mixed event
  TFile *fmix = TFile::Open(Form("Rootfiles/%s.Jpsi.pt%1.1f.pt%1.1f.root",run_type.Data(),pt1_cut,pt2_cut),"read");
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
  TCanvas *c = drawHistos(list,"Mix_Acceptance",Form("%s: acceptance effect from mixed event;M_{#mu#mu} (GeV/c^{2});UL/LS",run_type.Data()),true,1,5,true,0.85,1.2,false,kTRUE,legName,kTRUE,"p_{T} > 0 GeV/c",0.45,0.7,0.65,0.88,kTRUE,0.04);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/MixEvt_AccEff.pdf",run_type.Data()));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_MixEvtAcc.pdf"));

  // Get the # of UL-LS pairs from data
  TFile *fscan = TFile::Open(Form("Rootfiles/%s.TrkResScan.root",run_type.Data()),"read");
  TH1F *hEmbJpsiWidth = (TH1F*)fscan->Get(Form("SmearEmb_JpsiWidthIntegr_Pt0-15_def"));

  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.%sEmbJpsiEff.pt%1.1f.pt%1.1f.%sroot",run_type.Data(),gEmbTypeName[iMbEmb],pt1_cut,pt2_cut,run_config),"read");

  TH1F *hJpsiInvMass[2][nCentBinsEff][gNTrgSetup][gNZdcRate];
  TH1F *hJpsiCount[gNTrgSetup][gNZdcRate][3]; // 0 - fit US-LS; 1 - count US-LS; 2 - LS; 3 - MB events
  double min_mass[2] = {2.6, 3.2};
  double max_mass[2] = {2.9, 3.5};
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int p=0; p<gNZdcRate; p++)
	{
	  for(int k=0; k<3; k++)
	    {
	      hJpsiCount[t][p][k]  = new TH1F(Form("Data_JpsiCount_Zdc%d-%d%s_M%d",p*10,p*10+10,gTrgSetupName[t],k),";;Count",nCentBinsEff,0,nCentBinsEff);
	      for(int j=0; j<nCentBinsEff; j++)
		{
		  hJpsiCount[t][p][k]->GetXaxis()->SetBinLabel(j+1, Form("%s%%",centNameEff[j]));
		  hJpsiCount[t][p][k]->GetXaxis()->SetLabelSize(0.045);
		}
	    }

	  if(t==0 && (p==1 || p==4) )
	    {
	      c = new TCanvas(Form("ULvsLS_ZDC%d%s",p,gTrgSetupName[t]),Form("ULvsLS_ZDC%d-%d%s",p*10,p*10+10,gTrgSetupName[t]),1100,700);
	      c->Divide(4,4);
	    }
	  for(int j=0; j<nCentBinsEff; j++)
	    {
	      for(int i=0; i<2; i++)
		{
		  hJpsiInvMass[i][j][t][p] = (TH1F*)fin->Get(Form("Data_InvMass_%s_cent%s_zdc%d-%d%s",pName[i],centTitleEff[j],p*10,p*10+10,gTrgSetupName[t]));
		  int nbins = hJpsiInvMass[i][j][t][p]->GetNbinsX();
		  for(int bin=1; bin<=nbins; bin++)
		    {
		      hJpsiInvMass[i][j][t][p]->SetBinError(bin, sqrt(hJpsiInvMass[i][j][t][p]->GetBinContent(bin)));
		      if(i==1)
			{
			  int jbin = hAcc[0]->FindFixBin(hJpsiInvMass[i][j][t][p]->GetBinCenter(bin));
			  double acc = hAcc[0]->GetBinContent(jbin);
			  hJpsiInvMass[i][j][t][p]->SetBinContent(bin, acc*hJpsiInvMass[i][j][t][p]->GetBinContent(bin));
			  hJpsiInvMass[i][j][t][p]->SetBinError(bin, acc*hJpsiInvMass[i][j][t][p]->GetBinError(bin));
			}
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
	      if(nCount[1]<=0 || nCount[0]<=0) scale = 1;
	      double yield, error;

	      // count US-LS
	      hJpsiInvMass[0][j][t][p]->Add(hJpsiInvMass[1][j][t][p], -1*scale);
	      yield = hJpsiInvMass[0][j][t][p]->IntegralAndError(hJpsiInvMass[0][j][t][p]->FindFixBin(3.0), hJpsiInvMass[0][j][t][p]->FindFixBin(3.2), error);
	      if(yield<0) { yield = 0; error = 0; }
	      //cout << hJpsiCount[t][p][1]->GetName() << "  " << yield << "  " << scale << endl;
	      hJpsiCount[t][p][1]->SetBinContent(j+1, yield);
	      hJpsiCount[t][p][1]->SetBinError(j+1, error);

	      // count LS
	      yield = hJpsiInvMass[1][j][t][p]->IntegralAndError(hJpsiInvMass[1][j][t][p]->FindFixBin(3.0), hJpsiInvMass[1][j][t][p]->FindFixBin(3.2), error);
	      if(yield<0) { yield = 0; error = 0; }
	      hJpsiCount[t][p][2]->SetBinContent(j+1, yield);
	      hJpsiCount[t][p][2]->SetBinError(j+1, error);

	      // fit US-LS
	      hJpsiInvMass[0][j][t][p]->Rebin(5);
	      TF1 *func = new TF1(Form("func_%s",hJpsiInvMass[0][j][t][p]->GetName()), "gausn(0)+pol1(3)", 1.5, 4);
	      func->FixParameter(1, 3.096);
	      if(nCentBinsEff==8)  func->FixParameter(2, hEmbJpsiWidth->GetBinContent(j+1));
	      if(nCentBinsEff==16) func->FixParameter(2, hEmbJpsiWidth->GetBinContent(j/2+1));
	      hJpsiInvMass[0][j][t][p]->Fit(func, "IR0Q");
	      if(func->GetParameter(0)>0 && func->GetParameter(0)<1e3)
		{
		  hJpsiCount[t][p][0]->SetBinContent(j+1, func->GetParameter(0)/hJpsiInvMass[0][j][t][p]->GetBinWidth(1));
		  hJpsiCount[t][p][0]->SetBinError(j+1, func->GetParError(0)/hJpsiInvMass[0][j][t][p]->GetBinWidth(1));
		}
	      else
		{
		  hJpsiCount[t][p][0]->SetBinContent(j+1, 0);
		  hJpsiCount[t][p][0]->SetBinError(j+1, 0);
		}

	      if(t==0 && (p==1 || p==4))
		{
		  c->cd(j+1);
		  hJpsiInvMass[0][j][t][p]->GetXaxis()->SetRangeUser(2.5,4);
		  hJpsiInvMass[0][j][t][p]->SetMarkerStyle(21);
		  hJpsiInvMass[0][j][t][p]->SetMaximum(3*hJpsiInvMass[0][j][t][p]->GetMaximum());
		  hJpsiInvMass[0][j][t][p]->Draw("P");
		  func->SetLineColor(2);
		  func->SetLineStyle(2);
		  func->Draw("sames");

		  TPaveText *t1 = GetTitleText(Form("%s%%, %d < ZDC < %d kHz",centNameEff[j],p*10,p*10+10),0.06);
		  t1->Draw();
		}

	    }
	  if(t==0 && (p==1 || p==4) )
	    {
	      if(savePlot)
		c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_ULvsLS_Zdc%d-%d.pdf",run_type.Data(),p*10,p*10+10));
	      if(gSaveAN)
		{
		  c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_DataULvsLS_ZDC%d.pdf",p));
		}
	    }
	}
    }

  // Apply TPC efficiency
  TH1F *hMcJpsiEff[gNZdcRate]; 
  for(int p=0; p<gNZdcRate; p++)
    {
      hMcJpsiEff[p] = (TH1F*)fin->Get(Form("JpsiTpcEffVsCent_Pt0_Zdc%d-%d",p*10,p*10+10));
    }
  TH1F *hJpsiCorr[gNTrgSetup][gNZdcRate][4];
  double nJpsiAll[gNTrgSetup][4];
  for(int t=0; t<gNTrgSetup; t++)
    {
      for(int k=0; k<3; k++)
	{
	  nJpsiAll[t][k] = 0;

	  for(int p=0; p<gNZdcRate; p++)
	    {
	      hJpsiCorr[t][p][k] = (TH1F*)hJpsiCount[t][p][k]->Clone(Form("Data_JpsiCountWeight_Zdc%d-%d%s_M%d",p*10,p*10+10,gTrgSetupName[t],k));
	      int nbins = hJpsiCorr[t][p][k]->GetNbinsX();
	      for(int bin=1; bin<=nbins; bin++)
		{
		  double nJpsi = hJpsiCorr[t][p][k]->GetBinContent(bin);
		  double eff = 0;
		  if(nCentBinsEff==8)  eff = hMcJpsiEff[p]->GetBinContent(bin);
		  if(nCentBinsEff==16) eff = hMcJpsiEff[p]->GetBinContent((bin-1)/2+1);
		  
		  if(eff>0)
		    {
		      hJpsiCorr[t][p][k]->SetBinContent(bin, nJpsi/eff);
		      hJpsiCorr[t][p][k]->SetBinError(bin, hJpsiCorr[t][p][k]->GetBinError(bin)/eff);
		    }
		  else
		    {
		      hJpsiCorr[t][p][k]->SetBinContent(bin, 0);
		      hJpsiCorr[t][p][k]->SetBinError(bin, 0);
		    }
		  if(p==0)
		    {
		      hJpsiCorr[t][p][k]->SetBinContent(bin, 0);
		      hJpsiCorr[t][p][k]->SetBinError(bin, 0);
		    }
		  nJpsiAll[t][k] += hJpsiCorr[t][p][k]->GetBinContent(bin);
		}
	    }

	  if(k==0) printf("[i] Total # of corrected Jpsi: %1.0f for %s\n",nJpsiAll[t][k],gTrgSetupName[t]);
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      hJpsiCorr[t][p][k]->Scale(1./nJpsiAll[t][k]);
	    }
	}
    }

  // Get weights from MB events
  TFile *fmb = TFile::Open("output/Run14_AuAu200.MB.P16id.VtxEff.root", "read");
  TH3F	*hEqMbEvtWeight = (TH3F*)fmb->Get("mhEqMbEvtWeight"); //Centrality vs. ZdcRate vs. TrgSetup
  if(nCentBinsEff==8)
    {
      const double ncoll[8] = {13.09, 29.66, 62.18, 120.06, 216.45, 366.44, 593.67, 941.24}; 
      const double raa[8] = {1, 0.87, 0.75, 0.63, 0.58, 0.55, 0.45, 0.41};
    }
  if(nCentBinsEff==16)
    {
      const double ncoll[nCentBinsEff] = {10.54, 15.98, 23.96, 35.64, 51.62, 72.79, 101.51, 138.71, 186.69, 246.95, 320.78, 411.86, 524.31, 663.03, 838.41, 1048.11}; 
      const double raa[nCentBinsEff] = {1, 1, 1, 0.87, 0.75, 0.68, 0.63, 0.60, 0.58, 0.56, 0.55, 0.51, 0.48, 0.44, 0.42, 0.41};
    }

  if(nCentBinsEff==8)  hEqMbEvtWeight->Rebin3D(2, 2, 1);
  if(nCentBinsEff==16) hEqMbEvtWeight->Rebin3D(1, 2, 1);
  TH2F *hEqMbEvt[gNTrgSetup];
  for(int t=0; t<gNTrgSetup; t++)
    {
      if(t>0) hEqMbEvtWeight->GetZaxis()->SetRange(t,t);
      hEqMbEvt[t] = (TH2F*)hEqMbEvtWeight->Project3D("yx");
      hEqMbEvt[t]->SetName(Form("hEqMbEvt_%s",gTrgSetupTitle[t]));
    }

  for(int t=0; t<gNTrgSetup; t++)
    {
      nJpsiAll[t][3] = 0;
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorr[t][p][3] = (TH1F*)hJpsiCount[t][p][0]->Clone(Form("Data_JpsiCountWeight_Zdc%d-%d%s_M%d",p*10,p*10+10,gTrgSetupName[t],3));
	  hJpsiCorr[t][p][3]->Reset("AC");
	  int nbins = hJpsiCorr[t][p][3]->GetNbinsX();
	  for(int bin=1; bin<=nbins; bin++)
	    {
	      hJpsiCorr[t][p][3]->SetBinContent(bin, hEqMbEvt[t]->GetBinContent(nCentBinsEff-bin+1, p+1)*ncoll[nCentBinsEff-bin]*raa[nCentBinsEff-bin]);
	      hJpsiCorr[t][p][3]->SetBinError(bin, hEqMbEvt[t]->GetBinError(nCentBinsEff-bin+1, p+1)*ncoll[nCentBinsEff-bin]*raa[nCentBinsEff-bin]);
	      nJpsiAll[t][3] += hJpsiCorr[t][p][3]->GetBinContent(bin);
	    }
	}
      for(int p=0; p<gNZdcRate; p++)
	{
	  hJpsiCorr[t][p][3]->Scale(1./nJpsiAll[t][3]);
	}
    }

  // compare different methods
  const TString legMethod[4] = {"Fit US-LS", "Count US-LS", "Count LS", "Eq. MB events"};
  for(int t=0; t<gNTrgSetup; t++)
    {
      TCanvas *c = new TCanvas(Form("cJpsiWeights_%s",gTrgSetupTitle[t]), Form("cJpsiWeights_%s",gTrgSetupTitle[t]), 1100, 700);
      c->Divide(4,3);
      TLegend *leg = new TLegend(0.2,0.4,0.8,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06); 
      leg->SetHeader(gLegNameTrg[t].Data());
      for(int p=0; p<gNZdcRate; p++)
	{
	  c->cd(p+1);
	  gPad->SetLogy();
	  for(int k=0; k<4; k++)
	    {
	      TH1F *hclone = (TH1F*)hJpsiCorr[t][p][k]->Clone(Form("%s_clone",hJpsiCorr[t][p][k]->GetName()));
	      hclone->SetMarkerStyle(20+k);
	      hclone->SetMarkerColor(gColor[k]);
	      hclone->SetLineColor(gColor[k]);
	      hclone->GetYaxis()->SetRangeUser(1e-5, 5e-1);
	      if(k==0) hclone->Draw();
	      else     hclone->Draw("samesPE");
	      if(p==0) leg->AddEntry(hclone, legMethod[k].Data(), "PL");
	    }
	  TPaveText *t1 = GetTitleText(Form("%d < ZDC < %d kHz",p*10,p*10+10),0.065);
	  t1->Draw();
	}
      c->cd(12);
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_CompJpisWeight%s.pdf",run_type.Data(),gTrgSetupTitle[t]));
    }

  // plot nJpsi vs. centrality vs. ZDC
  const int methodIndex = 3;
  TLegend *leg = new TLegend(0.15, 0.65, 0.75, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  leg->SetNColumns(2);
  const int color[6] = {1, 2, 3, 4, 6, 7};
  int nJpsiCounts[gNTrgSetup];
  for(int t=0; t<gNTrgSetup; t++)
    {
      nJpsiCounts[t] = 0;
      for(int p=0; p<gNZdcRate; p++)
	{
	  for(int k=0; k<4; k++)
	    {
	      if(k<3)
		{
		  hJpsiCount[t][p][k]->GetXaxis()->SetLabelSize(0.05);
		  hJpsiCount[t][p][k]->SetMarkerStyle(20+p);
		  hJpsiCount[t][p][k]->SetMarkerColor(color[p%6]);
		  hJpsiCount[t][p][k]->SetLineColor(color[p%6]);
		  hJpsiCount[t][p][k]->GetYaxis()->SetRangeUser(1,1e6);
		  if(k==0) nJpsiCounts[t] += hJpsiCount[t][p][k]->Integral();
		}
	      hJpsiCorr[t][p][k]->GetXaxis()->SetLabelSize(0.05);
	      hJpsiCorr[t][p][k]->SetMarkerStyle(20+p);
	      hJpsiCorr[t][p][k]->SetMarkerColor(color[p%6]);
	      hJpsiCorr[t][p][k]->SetLineColor(color[p%6]);
	      hJpsiCorr[t][p][k]->GetYaxis()->SetRangeUser(1e-6,100);
	    }
	}
      printf("[i] %s: # of Jpsi = %4.2f\n",gTrgSetupName[t],nJpsiCounts[t]);
    }

  const int trgSetupIndex = 0;
  for(int p=0; p<gNZdcRate; p++)
    {
      leg->AddEntry(hJpsiCount[trgSetupIndex][p][0],Form("%d < ZDC < %d kHz",p*10,p*10+10),"P");
      if(p==0)
	{
	  c = draw1D(hJpsiCount[trgSetupIndex][p][0]);
	  gPad->SetLogy();
	  TPaveText *t1 = GetTitleText(Form("%s: # of J/psi counts from data",run_type.Data()),0.04);
	  t1->Draw();
	}
      else     hJpsiCount[trgSetupIndex][p][0]->Draw("sames");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_NJpsi_US-LS.pdf",run_type.Data()));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_RawNJpsi.pdf"));

  // plot weight vs. centrality vs. ZDC
  for(int p=0; p<gNZdcRate; p++)
    {
      if(p==0) 
	{
	  c = draw1D(hJpsiCorr[trgSetupIndex][p][methodIndex]);
	  gPad->SetLogy();
	  TPaveText *t1 = GetTitleText(Form("%s: weights for efficiency calculation",run_type.Data()),0.04);
	  t1->Draw();
	}
      else     
	hJpsiCorr[trgSetupIndex][p][methodIndex]->Draw("sames");
    }
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_EmbJpsiEff/Data_Jpsi_Weights.pdf",run_type.Data()));
  if(gSaveAN) c->SaveAs(Form("~/Dropbox/STAR_Quarkonium/Run14_Jpsi/Analysis_note/Figures/Ch4_EffTotal_CorrNJpsi.pdf"));

  if(saveHisto)
    {
      printf("+++ Save histograms +++\n");
      fin->cd();
      for(int t=0; t<gNTrgSetup; t++)
	{
	  for(int p=0; p<gNZdcRate; p++)
	    {
	      for(int k=0; k<4; k++)
		{
		  hJpsiCorr[t][p][k]->Write("", TObject::kOverwrite);
		}
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffSmear_CentBins.pdf",run_type.Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffSmear_CentBins.png",run_type.Data()));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffResp_CentBins.pdf",run_type.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffResp_CentBins.png",run_type.Data()));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffTotal_CentBins.pdf",run_type.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffTotal_CentBins.png",run_type.Data()));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.pdf",run_type.Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_centrality.png",run_type.Data()));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.pdf",run_type.Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMult.png",run_type.Data()));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.pdf",run_type.Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_gRefMultCorr.png",run_type.Data()));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.pdf",run_type.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/EmbedVsData_TpcVz.png",run_type.Data()));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%s%s_JpsiDis_CentBins.pdf",run_type.Data(),run_config,name[i].Data()));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%s%s_JpsiDis_CentBins.png",run_type.Data(),run_config,name[i].Data()));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MCJpsiPt_%s_CentBins.pdf",run_type.Data(),legNameTitle[j]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MCJpsiPt_%s_CentBins.png",run_type.Data(),legNameTitle[j]));
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
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiPtEff_CentBins.pdf",run_type.Data(),run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiPtEff_CentBins.png",run_type.Data(),run_config));
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiEffEmbed_CentBins.pdf",run_type.Data(),run_config));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/%sJpsiEffEmbed_CentBins.png",run_type.Data(),run_config));
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


