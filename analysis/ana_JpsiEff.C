const int year = YEAR;
const char *trkEffType[5] = {"MC","Tpc","MtdMth","MuonPid","MtdTrig"};
const char *weight_name[2] = {"","_w"};

//================================================
void ana_JpsiEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //trackingEff();
  trackingRes();
  //JpsiEff();
  //makeTrigEff();

  //ploEff();
  //plotEmbedEff();
  //makeMtdRespEff();
  //makeTrigEff();
  //makeEmbedEff();
  //hadronEmbed();
  //compTrigEff();
}


//================================================
void JpsiEff(const int savePlot = 1, const int saveHisto = 1)
{
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.JpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.JpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");

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
  TString legName_trg[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      legName_trg[k] = Form("AuAu_200%s",gTrgSetupTitle[k+1]);
    }

  TH1F *hJpsiPt[5][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPtRebin[5][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPtEff[5][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPtEffRebin[5][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPtEffFinal[gNTrgSetup][nCentBins][2];
  TH1F *hJpsiPtEffFinalRebin[gNTrgSetup][nCentBins][2];

  for(int i=0; i<5; i++)
    {
      for(int w=0; w<2; w++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      for(int k=0; k<nCentBins; k++)
		{
		  hJpsiPt[i][j][k][w] = (TH1F*)fin->Get(Form("hJpsiPt_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j],weight_name[w]));
		  hJpsiPtRebin[i][j][k][w] = (TH1F*)hJpsiPt[i][j][k][w]->Rebin(nbins,Form("%s_rebin",hJpsiPt[i][j][k][w]->GetName()),xbins);

		  if(i>0)
		    {
		      hJpsiPtEff[i][j][k][w] = DivideTH1ForEff(hJpsiPt[i][j][k][w],hJpsiPt[i-1][j][k][w],Form("hJpsiPtEff_%s_cent%s%s%s",trkEffType[i],cent_Title[k],gTrgSetupName[j],weight_name[w]));
		      hJpsiPtEffRebin[i][j][k][w] = DivideTH1ForEff(hJpsiPtRebin[i][j][k][w],hJpsiPtRebin[i-1][j][k][w],Form("hJpsiPtEff_%s_cent%s%s%s_rebin",trkEffType[i],cent_Title[k],gTrgSetupName[j],weight_name[w]));
		    }

		  if(i==4)
		    {
		      hJpsiPtEffFinal[j][k][w] = DivideTH1ForEff(hJpsiPt[i][j][k][w],hJpsiPt[0][j][k][w],Form("hJpsiPtEff_Final_cent%s%s%s",cent_Title[k],gTrgSetupName[j],weight_name[w]));
		      hJpsiPtEffFinalRebin[j][k][w] = DivideTH1ForEff(hJpsiPtRebin[i][j][k][w],hJpsiPtRebin[0][j][k][w],Form("JpsiPtEff_Final_cent%s%s%s_rebin",cent_Title[k],gTrgSetupName[j],weight_name[w]));
		    }
		}
	    }
	}
    }

  // various efficiency
  const int jsetup = 4, kcent = 1;
  list->Add(hJpsiPtEffFinal[jsetup][kcent][0]);
  list->Add(hJpsiPtEffFinalRebin[jsetup][kcent][0]);
  list->Add(hJpsiPtEffFinalRebin[jsetup][kcent][1]);
  TString legName1[3] = {"Unweighted","Unweighted && rebinned","Weighted && rebinned"};
  c = drawHistos(list,"JpsiEff_Check_Weight_Rebin",Form("Combined efficiency for J/#psi (AuAu_200%s);p_{T} (GeV/c);Efficiency",gTrgSetupTitle[jsetup]),false,2.0,3.8,true,1e-4,0.1,true,kTRUE,legName1,true,Form("%s%%",cent_Name[kcent]),0.5,0.7,0.2,0.45,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEffFinal_Check_Weight_Rebin.pdf",run_type));
  list->Clear();

  for(int i=1; i<5; i++)
    {
      list->Add(hJpsiPtEffRebin[i][jsetup][kcent][1]);
    }
  TString legName2[4] = {"TPC tracking + p_{T,#mu} cut","MTD matching","Muon PID","MTD triggering"};
  c = drawHistos(list,"JpsiEff_AllEffs",Form("Efficiency for J/#psi (AuAu_200%s);p_{T} (GeV/c);Efficiency",gTrgSetupTitle[jsetup]),false,2.0,3.8,true,1e-5,1.5,false,kTRUE,legName2,true,Form("%s%%",cent_Name[kcent]),0.2,0.4,0.6,0.88,kTRUE);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/JpsiEff_Various.pdf",run_type));
  list->Clear();

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

  /*
  for(int k=0; k<nCentBins; k++)
    {
      printf("+++++ %s%% +++++\n",cent_Name[k]);
      TH1F *hBase = (TH1F*)hJpsiPtRebin[0][0][0][0]->Clone(Form("hBase_%d",k));
      hBase->Reset();
      TH1F *hMatch = (TH1F*)hJpsiPtRebin[0][0][0][0]->Clone(Form("hMatch_%d",k));
      hMatch->Reset();
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  if(k==0)
	    {
	      for(int icent=0; icent<nCentBins-1; icent++)
		{
		  hBase->Add(hJpsiPtRebin[0][j+1][icent+1][1], nJpsi[icent][j]/nJpsiAll);
		  hMatch->Add(hJpsiPtRebin[4][j+1][icent+1][1], nJpsi[icent][j]/nJpsiAll);
		}
	    }
	  else
	    {
	      hBase->Add(hJpsiPtRebin[0][j+1][k][1], nJpsi[k-1][j]/nJpsiCent[k-1]);
	      hMatch->Add(hJpsiPtRebin[4][j+1][k][1], nJpsi[k-1][j]/nJpsiCent[k-1]);
	    }
	}
      hJpsiPtEffAvg[k] = DivideTH1ForEff(hMatch,hBase,Form("JpsiPtEff_cent%s",cent_Title[k]));
      for(int bin=1; bin<=hJpsiPtEffAvg[k]->GetNbinsX(); bin++)
	{
	  printf("[i] pt = %2.1f with eff = %4.2e, err = %4.2f%%\n",hJpsiPtEffAvg[k]->GetBinCenter(bin),hJpsiPtEffAvg[k]->GetBinContent(bin),hJpsiPtEffAvg[k]->GetBinError(bin)/hJpsiPtEffAvg[k]->GetBinContent(bin)*100);
	}
      list->Add(hJpsiPtEffAvg[k]);
    }
  */
  


  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiPtEffAvg[k]->Write("",TObject::kOverwrite);
	  hJpsiPtEffAvgCorr[k]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void trackingRes(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"read");
  TList *list = new TList;
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }

  // track momentum resolution
  TH2F *hResVsTruePt[nCentBins];
  TH1F *hTrkPtRes[nCentBins];
  TH1F *hTrkPtShift[nCentBins];
  TF1 *funcRes[nCentBins];
  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      hResVsTruePt[k] = (TH2F*) fin->Get(Form("PrimTrkRes_vs_TruePt_cent%s",cent_Title[k]));
      hResVsTruePt[k]->RebinX(2);
      hResVsTruePt[k]->GetXaxis()->SetRangeUser(0,20);
      char *title = Form("Transverse momentum resolution of primary muon tracks (%s%%)",cent_Name[k]);
      c = draw2D(hResVsTruePt[k],title);
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/DeltaPt_vs_Pt_cent%s.pdf",run_type,cent_Title[k]));

      hTrkPtRes[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("PrimTrkRes_cent%s",cent_Title[k]));
      hTrkPtRes[k]->Reset();
      hTrkPtShift[k] = (TH1F*)hResVsTruePt[k]->ProjectionX(Form("PrimTrkShift_cent%s",cent_Title[k]));
      hTrkPtShift[k]->Reset();

      TCanvas *c = new TCanvas(Form("FitTrkRes_cent%d",k),Form("FitTrkRes_cent%d",k),800,600);
      c->Divide(10,10);
      for(int ibin=1; ibin<=hTrkPtRes[k]->GetNbinsX(); ibin++)
	{
	  TH1F *htmp = (TH1F*)hResVsTruePt[k]->ProjectionY(Form("TrkPt_bin%d_cent%d",ibin,k),ibin,ibin);
	  htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
	  htmp->SetMarkerStyle(20);
	  TF1 *func = new TF1(Form("func_bin%d",ibin),"gaus",-0.15,0.15);
	  htmp->Fit(func,"IR0Q");
	  c->cd(ibin);
	  htmp->Draw("P");
	  func->SetLineColor(2);
	  func->Draw("same");
	  hTrkPtRes[k]->SetBinContent(ibin,func->GetParameter(2));
	  hTrkPtRes[k]->SetBinError(ibin,func->GetParError(2));
	  hTrkPtShift[k]->SetBinContent(ibin,func->GetParameter(1));
	  hTrkPtShift[k]->SetBinError(ibin,func->GetParError(1));
	}

      TH1F *htmp = (TH1F*)hTrkPtRes[k]->Clone(Form("%s_clone",hTrkPtRes[k]->GetName()));
      list->Add(htmp);

      funcRes[k] = new TF1(Form("FitPrimTrkRes_cent%s",cent_Title[k]),"sqrt([0]^2*x^2+[1]^2)",1,10);
      funcRes[k]->SetParNames("a","b");
      funcRes[k]->SetParameter(0.0004,0.0001);
      hTrkPtRes[k]->Fit(funcRes[k],"IR0");
      hTrkPtRes[k]->SetMarkerStyle(21);
      hTrkPtRes[k]->GetYaxis()->SetRangeUser(0,0.15);
      c = draw1D(hTrkPtRes[k],Form("Transverse momentum resolution of primary tracks (%s%%);p_{T,true} (GeV/c);#sigma(p_{T})/p_{T}",cent_Name[k]));
      funcRes[k]->SetLineColor(2);
      funcRes[k]->Draw("same");
      leg = new TLegend(0.2,0.5,0.4,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hTrkPtRes[k],"Embedding data","P");
      leg->AddEntry(funcRes[k],"Fit: #sqrt{(a*p_{T})^{2}+b^{2}}","L");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitTrkPtRes_cent%s.pdf",run_type,cent_Title[k]));
    }
  bool drawLegend = true;
  if(year==2013) drawLegend = false;
  c = drawHistos(list,Form("TrkRes"),Form("p_{T} resolution of muon tracks;p_{T,MC} (GeV/c);#sigma(p_{T})/p_{T}"),kTRUE,0,20,kTRUE,0,0.1,kFALSE,drawLegend,legName_cent,drawLegend,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTrkPtRes_CentBins.pdf",run_type));

  list->Clear();
  for(int k=0; k<nCentBins; k++)
    {
      TH1F *htmp = (TH1F*)hTrkPtShift[k]->Clone(Form("%s_clone",hTrkPtShift[k]->GetName()));
      htmp->Scale(100);
      list->Add(htmp);
    }
  c = drawHistos(list,Form("TrkPtShift"),Form("p_{T} shift of muon tracks;p_{T,MC} (GeV/c);<#Deltap_{T}/p_{T}> (%%)"),kTRUE,0,20,kTRUE,-1,1,kFALSE,drawLegend,legName_cent,drawLegend,"",0.3,0.5,0.6,0.85,kTRUE);
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTrkPtShift_CentBins.pdf",run_type));

  if(saveHisto)
    {
      fin->cd();
      for(int k=0; k<nCentBins; k++)
	{
	  hTrkPtRes[k]->Write("",TObject::kOverwrite);
	  hTrkPtShift[k]->Write("",TObject::kOverwrite);
	  funcRes[k]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void trackingEff(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fin;
  if(saveHisto)   fin = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"update");
  else            fin = TFile::Open(Form("Rootfiles/%s.TrkEff.root",run_type),"read");
  TList *list = new TList;
  TString legName_cent[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      legName_cent[k] = Form("%s%%",cent_Name[k]);
    }
  TString legName_trg[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      legName_trg[k] = Form("AuAu_200%s",gTrgSetupTitle[k+1]);
    }

  // tracking efficiency
  TH1F *hMcTrkPt[5][gNTrgSetup][nCentBins];
  TH1F *hMcTrkPtEff[5][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsZdc[5][gNTrgSetup][nCentBins];
  TH2F *hMcTrkPtVsCent[5][gNTrgSetup];

  for(int i=0; i<5; i++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  hMcTrkPtVsCent[i][j] = (TH2F*)fin->Get(Form("hMcTrkPtVsCent_%s%s",trkEffType[i],gTrgSetupTitle[j]));
	  hMcTrkPtVsCent[i][j]->RebinX(2);
	  for(int k=0; k<nCentBins; k++)
	    {
	      hMcTrkPtVsZdc[i][j][k] = (TH2F*)fin->Get(Form("hMcTrkPtVsZdcRate_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));

	      hMcTrkPt[i][j][k] = (TH1F*)fin->Get(Form("hMcTrkPt_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
	      hMcTrkPt[i][j][k]->Rebin(5);
	      if(i>0) 
		{
		  hMcTrkPtEff[i][j][k] = DivideTH1ForEff(hMcTrkPt[i][j][k], hMcTrkPt[i-1][j][k], Form("hMcTrkPtEff_%s_cent%s%s",trkEffType[i],cent_Title[k],gTrgSetupTitle[j]));
		}	
	    }
	}
    }

  TCanvas *c = new TCanvas("TpcEff_Lumi","TpcEff_Lumi",1100,700);
  c->Divide(2,2);
  TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  for(int k=1; k<nCentBins; k++)
    {
      c->cd(k);
      SetPadMargin(gPad,0.13, 0.13, 0.03,0.1);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  TH1F *htmp = (TH1F*)hMcTrkPtEff[1][j][k]->Clone(Form("%s_clone",hMcTrkPtEff[1][j][k]->GetName()));
	  htmp->SetMarkerStyle(20+j-1);
	  htmp->SetLineColor(color[j-1]);
	  htmp->SetMarkerColor(color[j-1]);
	  htmp->GetYaxis()->SetRangeUser(0,1);
	  htmp->GetYaxis()->SetTitle("Efficiency");
	  ScaleHistoTitle(htmp,0.06,1,0.05,0.06,0.9,0.05,62);
	  if(j==1) htmp->Draw();
	  else     htmp->Draw("sames");
	  if(k==1)
	    {
	      leg->AddEntry(htmp,Form("AuAu_200%s",gTrgSetupTitle[j]),"PL");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency (%s%%)",cent_Name[k]),0.06);
      t1->Draw();
      TLine *line = GetLine(1.2,0,1.2,0.8,1);
      line->Draw();
    }
  c->cd(4);
  leg->SetHeader("|#eta_{MC}| < 0.5");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTpcEff_Lumi.pdf",run_type));

  const double pt_cut = 2;
  // efficiency vs luminosity
  TH1F *hTpcEffVsLumi[gNTrgSetup][nCentBins];
  for(int j=0; j<gNTrgSetup; j++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  int bin_cut = hMcTrkPtVsZdc[1][j][k]->GetYaxis()->FindFixBin(pt_cut);
	  hTpcEffVsLumi[j][k] = (TH1F*)hMcTrkPtVsZdc[1][j][k]->ProjectionX(Form("hTpcEffVsLumi_cent%s%s",cent_Title[k],gTrgSetupTitle[j]),bin_cut,-1);
	  TH1F *htmp = (TH1F*)hMcTrkPtVsZdc[0][j][k]->ProjectionX(Form("htmp_cent%s%s",cent_Title[k],gTrgSetupTitle[j]),bin_cut,-1);
	  hTpcEffVsLumi[j][k]->Divide(htmp);
	}
    }
  TCanvas *c = new TCanvas("TpcEff_vs_Lumi","TpcEff_vs_Lumi",1100,700);
  c->Divide(2,2);
  for(int k=1; k<nCentBins; k++)
    {
      c->cd(k);
      SetPadMargin(gPad,0.13, 0.13, 0.03,0.1);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  hTpcEffVsLumi[j][k]->SetMarkerStyle(20+j-1);
	  hTpcEffVsLumi[j][k]->SetLineColor(color[j-1]);
	  hTpcEffVsLumi[j][k]->SetMarkerColor(color[j-1]);
	  hTpcEffVsLumi[j][k]->GetYaxis()->SetRangeUser(0,1);
	  hTpcEffVsLumi[j][k]->GetXaxis()->SetRangeUser(0,120);
	  hTpcEffVsLumi[j][k]->GetYaxis()->SetTitle("Efficiency");
	  ScaleHistoTitle(hTpcEffVsLumi[j][k],0.06,1,0.05,0.06,0.9,0.05,62);
	  if(j==1) hTpcEffVsLumi[j][k]->Draw();
	  else     hTpcEffVsLumi[j][k]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("TPC tracking efficiency above %1.0f GeV/c (%s%%)",pt_cut,cent_Name[k]),0.06);
      t1->Draw();
    }
  c->cd(4);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTpcEff_vs_Lumi.pdf",run_type));

  // efficiency vs centrality
  list->Clear();
  TH1F *hTpcEffVsCent[gNTrgSetup];
  for(int j=1; j<gNTrgSetup; j++)
    {
      int bin_cut = hMcTrkPtVsCent[1][j]->GetYaxis()->FindFixBin(pt_cut);
      hTpcEffVsCent[j] = (TH1F*)hMcTrkPtVsCent[1][j]->ProjectionX(Form("hTpcEffVsCent%s",gTrgSetupTitle[j]),bin_cut,-1);
      TH1F *htmp = (TH1F*)hMcTrkPtVsCent[0][j]->ProjectionX(Form("htmp%s",gTrgSetupTitle[j]),bin_cut,-1);
      hTpcEffVsCent[j]->Divide(htmp);
      hTpcEffVsCent[j]->SetMarkerStyle(20+j-1);
      hTpcEffVsCent[j]->SetLineColor(color[j-1]);
      list->Add(hTpcEffVsCent[j]);
      for(int bin=1; bin<=hTpcEffVsCent[j]->GetNbinsX(); bin++)
	{
	  hTpcEffVsCent[j]->GetXaxis()->SetBinLabel(bin,Form("%d-%d%%",80-bin*10,90-bin*10));
	}
      hTpcEffVsCent[j]->GetXaxis()->SetLabelSize(0.05);
    }
  c = drawHistos(list,"TpcEff_vs_Cent",Form("TPC tracking efficiency for tracks above %1.0f GeV/c;;Efficiency",pt_cut),false,2.0,3.8,true,0,1,kFALSE,kTRUE,legName_trg,true,"|#eta_{MC}| < 0.5",0.5,0.7,0.2,0.45,kTRUE,0.04,0.04,false,1,false,true); 
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTpcEff_vs_Cent.pdf",run_type));
  
  // other efficiencies
  const char *trkEffTitle[3] = {"MTD matching","muon PID","MTD trigger"};
  for(int i=2; i<5; i++)
    {
      TLegend *leg1 = new TLegend(0.2,0.2,0.4,0.4);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.04);
      TLegend *leg2 = new TLegend(0.55,0.2,0.7,0.4);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.04);
      for(int j=1; j<gNTrgSetup; j++)
	{
	  for(int k=1; k<nCentBins; k++)
	    {
	      if(i==2)  hMcTrkPtEff[i][j][k]->GetYaxis()->SetRangeUser(0,0.5);
	      else      hMcTrkPtEff[i][j][k]->GetYaxis()->SetRangeUser(0,1.1);
	      hMcTrkPtEff[i][j][k]->SetMarkerStyle(19+k);
	      hMcTrkPtEff[i][j][k]->SetMarkerColor(color[j-1]);
	      if(j==1 && k==1) c = draw1D(hMcTrkPtEff[i][j][k],Form("Single muon efficiency of %s (|#eta_{#mu}| < 0.5);p_{T}^{mc} (GeV/c);Efficiency",trkEffTitle[i-2]));
	      else hMcTrkPtEff[i][j][k]->Draw("sames");
	      if(k==1) leg1->AddEntry(hMcTrkPtEff[i][j][k],Form("AuAu_200%s",gTrgSetupTitle[j]),"P");
	      if(j==1) leg2->AddEntry(hMcTrkPtEff[i][j][k],Form("%s%%",cent_Name[k]),"P");
	    }
	}
      leg1->Draw();
      leg2->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/McTrkEff_%s.pdf",run_type,trkEffType[i]));
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
  TFile *fin = 0x0, *fEmbed = 0x0;
  if(year==2014)
    {
      //cosmic
      fin = TFile::Open("Rootfiles/Run14ResponseEffViaPtTemplate.root","read");
      fEmbed =  TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
    }
  else if(year==2013)
    {
      fin = TFile::Open("Rootfiles/Run13ResponseEffViaPtTemplate.root","read");
      fEmbed =  TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
    }

  // embedding
  TH2F *hProjVsPt = (TH2F*)fEmbed->Get("mhProjTrack");
  draw2D(hProjVsPt);
  TH2F *hMatchVsPt = (TH2F*)fEmbed->Get("mhMatchTrack");
  draw2D(hMatchVsPt);
  TH1F *hProj = (TH1F*)hProjVsPt->ProjectionX("hProj");
  TH1F *hMatch = (TH1F*)hMatchVsPt->ProjectionX("hMatch");
  TH1F *hMtdRespEffEmbed = (TH1F*)hMatch->Clone(Form("MtdResponseEff_%s",name[1]));
  hMtdRespEffEmbed->Sumw2();
  hMtdRespEffEmbed->Divide(hProj);
  hMtdRespEffEmbed->SetMarkerStyle(20);
  hMtdRespEffEmbed->SetLineColor(1);
  hMtdRespEffEmbed->SetMarkerColor(1);
  hMtdRespEffEmbed->GetXaxis()->SetRangeUser(0,10);
  hMtdRespEffEmbed->GetYaxis()->SetRangeUser(0,1.2);
  ScaleHistoTitle(hMtdRespEffEmbed,0.045,1,0.035,0.045,1,0.035,62);

  // Fit embedding data
  TH1F *hEmbedFit = (TH1F*)hMtdRespEffEmbed->Clone("hEmbedFit");
  TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2])-pow([3]/x/x,[4])-pow([5]/x/x/x,[6]))",1,20);
  func->SetParameters(1,0.1,1,0.1,1,0.1,1);
  hEmbedFit->Fit(func,"IR0");
  hEmbedFit->GetYaxis()->SetRangeUser(0,1.8);
  c = draw1D(hEmbedFit,"MTD response efficiency from embedding;p_{T} (GeV/c);Efficiency");
  func->SetLineColor(4);
  func->Draw("sames");
  TLegend *leg = new TLegend(0.3,0.2,0.5,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEmbedFit,"Embedding","PE");
  leg->AddEntry(func,"Fit function: [0]*exp{-(#frac{[1]}{x})^{[2]}-(#frac{[3]}{x^{2}})^{[4]}-(#frac{[5]}{x^{3}})^{[6]}}","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitEmbedMtdRespEff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/FitEmbedMtdRespEff.png",run_type));
    }

  // cosmic
  TF1 *funcCosmic[150];
  for(int i=0; i<150; i++)
    {
      funcCosmic[i] = (TF1*)fin->Get(Form("fSclPtTmpBkl%d_Mod%d",i/5,i%5));
      funcCosmic[i]->SetLineColor(4);
      funcCosmic[i]->SetLineStyle(1);
    }
  hMtdRespEffEmbed->GetYaxis()->SetRangeUser(0,1.3);
  c = draw1D(hMtdRespEffEmbed,"MTD response efficiency;p_{T} (GeV/c);Efficiency");
  for(int i=0; i<150; i++)
    {
      funcCosmic[i]->Draw("sames");
    }
  TLegend *leg = new TLegend(0.2,0.75,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hMtdRespEffEmbed,"Embedding","PE");
  leg->AddEntry(funcCosmic[0],"Cosmic ray (module-wise)","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMtdRespEff.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/CompareMtdRespEff.png",run_type));
    }

  if(saveHisto)
    {
      TFile *fout = 0x0;
      if(year==2014) fout = TFile::Open("Rootfiles/Run14.AuAu200.MtdResponseEff.root","recreate");
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.MtdResponseEff.root","recreate");
      
      hMtdRespEffEmbed->Write();
      func->Write(Form("Fit_%s",hMtdRespEffEmbed->GetName()));
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
void makeTrigEff(const int savePlot = 1, const int saveHisto = 1)
{
  const char *trgSetupName[4] = {"AuAu_200_production_2014","AuAu_200_production_low_2014","AuAu_200_production_mid_2014","AuAu_200_production_high_2014"};
  TF1 *funcTrigEff[4][nCentBins];
  TF1 *funcTrigEffCorr[4][nCentBins];
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  funcTrigEff[i][k] = new TF1(Form("MuonTrigEff_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"[0]-exp(-1*[1]*(x-[2]))",0,12);
	  funcTrigEffCorr[i][k] = new TF1(Form("MuonTrigEffCorr_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"[0]-exp(-1*[1]*(x-[2]))",0,12);
	  if(i==0 || i==1)
	    {
	      if(k==0) funcTrigEff[i][k]->SetParameters(0.951456, 2.32223, 6.61693e-14);      
	      if(k==1) funcTrigEff[i][k]->SetParameters(0.948947, 2.55808, 0.0659758);
	      if(k==2) funcTrigEff[i][k]->SetParameters(0.950789, 2.51311, 7.49401e-14);
	      if(k==3) funcTrigEff[i][k]->SetParameters(0.968246, 2.45772, 8.27116e-15);

	      if(k==0) funcTrigEffCorr[i][k]->SetParameters(0.96023,  10.5511, 0.802917);      
	      if(k==1) funcTrigEffCorr[i][k]->SetParameters(0.96918,  9.40215, 0.778563);
	      if(k==2) funcTrigEffCorr[i][k]->SetParameters(0.964368, 3.65518, 0.298387);
	      if(k==3) funcTrigEffCorr[i][k]->SetParameters(0.966162, 2.77122, 3.6e-15);
	    }
	  else
	    {
	      if(k==0) funcTrigEff[i][k]->SetParameters(0.908, 2.09256, 1.03406e-14);      
	      if(k==1) funcTrigEff[i][k]->SetParameters(0.905465, 2.05892, 9.27036e-15);
	      if(k==2) funcTrigEff[i][k]->SetParameters(0.920748, 1.89851, 1.57097e-14);
	      if(k==3) funcTrigEff[i][k]->SetParameters(0.883885, 2.1119, 3.33068e-15);

	      if(k==0) funcTrigEffCorr[i][k]->SetParameters(0.936584, 17.1834, 0.893653);      
	      if(k==1) funcTrigEffCorr[i][k]->SetParameters(0.945798, 11.6541, 0.838267);
	      if(k==2) funcTrigEffCorr[i][k]->SetParameters(0.935776, 5.75038, 0.499548);
	      if(k==3) funcTrigEffCorr[i][k]->SetParameters(0.938118, 13.9462, 0.875862);
	    }
	}
    }

  for(int i=0; i<4; i++)
    {
      TCanvas *c = new TCanvas(Form("MuonTrigEff%s",gTrgSetupName[i+1]),Form("MuonTrigEff%s",gTrgSetupName[i+1]),1100,700);
      c->Divide(2,2);
      for(int k=0; k<nCentBins; k++)
	{
	  c->cd(k+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	  gPad->SetGrid(1,1);
	  ScaleHistoTitle(funcTrigEff[i][k]->GetHistogram(),0.06,1,0.05,0.06,0.9,0.05,62);
	  funcTrigEff[i][k]->SetTitle(Form("%s: %s%%;p_{T} (GeV/c);Trigger efficiency",trgSetupName[i],cent_Name[k]));
	  funcTrigEff[i][k]->SetMaximum(1.1);
	  funcTrigEff[i][k]->SetMinimum(0.5);
	  funcTrigEff[i][k]->SetLineColor(1);
	  funcTrigEff[i][k]->Draw();
	  funcTrigEffCorr[i][k]->SetLineColor(4);
	  funcTrigEffCorr[i][k]->Draw("sames");
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.45,0.25,0.65,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(funcTrigEff[i][0],"VPDMB5","L");
      leg->AddEntry(funcTrigEffCorr[i][0],"NoVtx/VPDMB5","L");
      leg->Draw();
  
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins%s.pdf",run_type,gTrgSetupName[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEff_CentBins%s.png",run_type,gTrgSetupName[i+1]));
	}
    }

  TH1F *hMuonTrigEff[4][nCentBins];
  TString legName[nCentBins];
  TList *list = new TList;
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hMuonTrigEff[i][k] = new TH1F(Form("CombinedMuonTrigEff_cent%s%s",cent_Title[k],gTrgSetupName[i+1]),"Muon efficiency",1000,0,10);
	  for(int bin=1; bin<=hMuonTrigEff[i][k]->GetNbinsX(); bin++)
	    {
	      double x = hMuonTrigEff[i][k]->GetXaxis()->GetBinCenter(bin);
	      hMuonTrigEff[i][k]->SetBinContent(bin,funcTrigEff[i][k]->Eval(x)*funcTrigEffCorr[i][k]->Eval(x));
	      hMuonTrigEff[i][k]->SetBinError(bin,0);
	    }
	  list->Add(hMuonTrigEff[i][k]);
	  legName[k] = Form("%s%%",cent_Name[k]);
	}
      c = drawHistos(list,Form("MuonTrigEff_CentBins%s",gTrgSetupName[i+1]),Form("Single muon trigger efficiency (%s);p_{T} (GeV/c);Efficiency",trgSetupName[i]),kTRUE,0.5,10,kTRUE,0.5,1.0,kFALSE,kTRUE,legName,kTRUE,"",0.4,0.6,0.3,0.55,kTRUE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins%s.pdf",run_type,gTrgSetupName[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiEff/MuonTrigEffCombined_CentBins%s.png",run_type,gTrgSetupName[i+1]));
	}
      list->Clear();
    }
  
  if(saveHisto)
    {
      char *outName = "Run14.AuAu200.MuonTrigEff.root";
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outName),"recreate");
      for(int i=0; i<4; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      funcTrigEff[i][k]->Write("",TObject::kOverwrite);
	      funcTrigEffCorr[i][k]->Write("",TObject::kOverwrite);
	      hMuonTrigEff[i][k]->Write("",TObject::kOverwrite);
	    }
	}
      fout->Close();
    }
}

