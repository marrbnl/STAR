TFile *f;
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;
TString file_name;

//================================================
void sys_JpsiYield()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2013)
    {
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //qa();
  //signalExtraction();
  //makeHistoSysCuts();
  //processSysCuts();
  //trackCuts();
  //triggerEfficiency();
  mergeSystematics();
}


//================================================
void mergeSystematics(int savePlot = 1, int saveHisto = 1)
{
  const int icent = 0;
  TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"update");
  const int nSys = 4;
  const char *name[nSys] = {"signalExt","TrackCuts","PidCuts","trigEff"};
  TH1F *hSysAll;
  TH1F *hSys[nSys];
  TH1F *hSysSub[nSys];
  for(int i=0; i<nSys; i++)
    {
      hSys[i] = (TH1F*)fout->Get(Form("Sys_%s_%s",name[i],cent_Title[icent]));
      hSysSub[i] = (TH1F*)hSys[i]->Clone(Form("Sys_%s_%s_sub",name[i],cent_Title[icent]));
      for(int bin=1; bin<=hSysSub[i]->GetNbinsX(); bin++)
	{
	  hSysSub[i]->SetBinContent(bin,hSysSub[i]->GetBinContent(bin)-1);
	}
    }
  hSysAll = (TH1F*)hSys[0]->Clone(Form("Sys_all_%s",cent_Title[icent]));
  hSysAll->Reset();
  for(int bin=1; bin<=hSysAll->GetNbinsX(); bin++)
    {
      double error = 0;
      for(int i=0; i<nSys; i++)
	{
	  error += TMath::Power(hSysSub[i]->GetBinContent(bin),2);
	}
      error = sqrt(error);
      hSysAll->SetBinContent(bin,error);
    }

  TList *list = new TList;
  TString legName[nSys+1] = {"Total","Signal extraction","Tracking efficiency","PID efficiency","Trigger efficiency"};
  list->Add(hSysAll);
  for(int i=0; i<nSys; i++)
    {
      list->Add(hSysSub[i]);
    }
  c = drawHistos(list,"JpsiYield_Sys","Systematic uncertainty for J/#psi spectrum;p_{T} (GeV/c);Uncertainties",kFALSE,0.5,10,kTRUE,0,0.35,kFALSE,kTRUE,legName,kTRUE,Form("%s%%",cent_Name[icent]),0.45,0.65,0.55,0.85,kFALSE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_All.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_All.png",run_type,run_cfg_name.Data()));
    }

  if(saveHisto)
    {
      fout->cd();
      hSysAll->Write("",TObject::kOverwrite);
    }
}

//================================================
void qa()
{
  const int icent = 0;
  const int ipt = 0;
  TCanvas *c = new TCanvas("qa","qa",800,600);
  const int nSys = 12;
  const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp",
			    "nSigmaPiUp","nSigmaPiDown",
			    "dzUp","dzDown","dyUp","dyDown","dtofUp"};

  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.sys.Cuts.root",run_config,pt1_cut,pt2_cut),"read");
  TH1F *histo[nSys];
  for(int s=0; s<nSys; s++)
    {
      histo[s] = (TH1F*)fin->Get((Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[ipt],cent_Title[icent])));
      histo[s]->SetMarkerColor(color[s]);
      histo[s]->SetLineColor(color[s]);
      if(s==0) histo[s]->Draw();
      else     histo[s]->Draw("sames");
    }
}

//================================================
void signalExtraction(int savePlot = 1, int saveHisto = 1)
{
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TH1F *hSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_signalExt_%s",cent_Title[i]),Form("Systematic uncertainty for signal extraction (%s%%)",cent_Name[i]),nbins,xbins);
    }

  const int nSys = 7;
  const char *sys_name[nSys] = {"","_LargeScale","_SmallScale","_LargeFit","_SmallFit","_Rebin","_pol1"};
  TString outName = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut);
  TString outNameSys = outName;
  outNameSys.ReplaceAll(".yield.root",".sys.signal.root");

  f = TFile::Open(Form("Rootfiles/%s",outName.Data()),"read");
  TFile *fSys = TFile::Open(Form("Rootfiles/%s",outNameSys.Data()),"read");
  TH1F *hSignal[nCentBins][nSys];

  const double sys_value = 0.06;
  TCanvas *c = new TCanvas(Form("Sys_signalExt"),Form("Sys_signalExt"),1100,700);
  c->Divide(2,2);
  double max[nCentBins][nPtBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int k=0; k<nPtBins; k++)
	{
	  max[i][k] = 0;
	}
    }

  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<6; j++)
	{
	  if(j==0) hSignal[i][j] = (TH1F*)f->Get(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[i],sys_name[j]));
	  else     hSignal[i][j] = (TH1F*)fSys->Get(Form("Jpsi_BinCountYield_cent%s%s",cent_Title[i],sys_name[j]));
	  hSignal[i][j]->SetMarkerSize(1.5);
	  hSignal[i][j]->SetMarkerStyle(21);
	  hSignal[i][j]->SetMarkerColor(color[j]);
	  
	  TH1F *htmp = (TH1F*)hSignal[i][j]->Clone(Form("%s_tmp",hSignal[i][j]->GetName()));
	  htmp->Divide(hSignal[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hSignal[i][j]->GetBinError(bin)/hSignal[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.4,1.7);
	  htmp->SetTitle(";p_{T} (GeV/c);Relative difference");
	  c->cd(i+1);
	  SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	  ScaleHistoTitle(htmp,0.06,1,0.05,0.06,1,0.05,62);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("Systematic uncertainty of signal extraction"),0.06);
	  t1->Draw();
	  t1 = GetPaveText(0.25,0.35,0.2,0.3,0.06,62);
	  t1->AddText(Form("%s%%",cent_Name[i]));
	  t1->Draw();
	  for(int k=0; k<nPtBins; k++)
	    {
	      double value = fabs(htmp->GetBinContent(k+1)-1);
	      if(max[i][k]<value) max[i][k] = value;
	    }
	}

      if(i==0)
	{
	  double xmin = 0.2, xmax = 0.5, ymin = 0.65, ymax = 0.88;
	  leg1 = new TLegend(xmin,ymin,xmax,ymax);
	  leg1->SetBorderSize(0);
	  leg1->SetFillColor(0);
	  leg1->SetTextSize(0.045);
	  leg1->AddEntry(hSignal[i][0],"Default with stat. err.","P");
	  leg1->AddEntry(hSignal[i][1],"Larger scale range","P");
	  leg1->AddEntry(hSignal[i][2],"Smaller scale range","P");
	  leg1->Draw();

	  leg2 = new TLegend(xmin+0.4,ymin,xmax+0.4,ymax);
	  leg2->SetBorderSize(0);
	  leg2->SetFillColor(0);
	  leg2->SetTextSize(0.045);
	  leg2->AddEntry(hSignal[i][3],"Larger fit range","P");
	  leg2->AddEntry(hSignal[i][4],"Smaller fit range","P");
	  leg2->AddEntry(hSignal[i][5],"Different binning","P");
	  leg2->Draw();
	}
      
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,max[0][bin-1]+1);
	  htmp->SetBinContent(bin,1-max[0][bin-1]);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_signalExt.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_signalExt.png",run_type,run_cfg_name.Data()));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"update");
      for(int i=0; i<nCentBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
      fout->Close();
    }
  
}


//================================================
void triggerEfficiency(int saveHisto = 1)
{
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TH1F *hSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_trigEff_%s",cent_Title[i]),Form("Systematic uncertainty for trigger efficiency (%s%%)",cent_Name[i]),nbins,xbins);
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,0.13+1);
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"update");
      for(int i=0; i<nCentBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
      fout->Close();
    }
}

//================================================
void trackCuts(int savePlot = 1, int saveHisto = 1)
{
  const int type = 0; // 0 - quality cuts; 1 - PID cuts
  const char *type_name[2] = {"TrackCuts","PidCuts"};
  const char *type_title[2] = {"track quality cuts","track PID cuts"};

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TH1F *hSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_%s_%s",type_name[type],cent_Title[i]),Form("Systematic uncertainty for %s (%s%%)",type_title[type],cent_Name[i]),nbins,xbins);
    }

  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.sys.Cuts.root",run_config,pt1_cut,pt2_cut),"read");

  if(type==0)
    {
      const int nSys = 5;
      const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp"};
    }
  else if(type==1)
    {
      const int nSys = 8;
      const char *name[nSys] = {"default","nSigmaPiUp","nSigmaPiDown", "dzUp","dzDown","dyUp","dyDown","dtofUp"};
    }

  TH1F *hYield[nCentBins][nSys];
  TH1F *hEff[nCentBins][nSys];
  TCanvas *c = new TCanvas(Form("Sys_signalExt"),Form("Sys_signalExt"),1100,700);
  c->Divide(2,2);
  double max[nCentBins][nPtBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int k=0; k<nPtBins; k++)
	{
	  max[i][k] = 0;
	}
    }
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      for(int j=0; j<nSys; j++)
	{
	  hYield[i][j] = (TH1F*)fin->Get(Form("Sys_%s_JpsiYield_cent%s",name[j], cent_Title[i]));
	  hEff[i][j] = (TH1F*)fin->Get(Form("Sys_%s_JpsiEff_cent%s",name[j], cent_Title[i]));
	  hYield[i][j]->Divide(hEff[i][j]);

	  hYield[i][j]->SetMarkerSize(1.5);
	  hYield[i][j]->SetMarkerStyle(21);
	  hYield[i][j]->SetMarkerColor(color[j]);
	  
	  TH1F *htmp = (TH1F*)hYield[i][j]->Clone(Form("%s_tmp",hYield[i][j]->GetName()));
	  htmp->Divide(hYield[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hYield[i][j]->GetBinError(bin)/hYield[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.4,1.7);
	  htmp->SetTitle(";p_{T} (GeV/c);Relative difference");
	  SetPadMargin(gPad,0.15,0.15,0.05,0.1);
	  ScaleHistoTitle(htmp,0.06,1,0.05,0.06,1,0.05,62);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("Systematic uncertainty of %s",type_title[type]),0.06);
	  t1->Draw();
	  t1 = GetPaveText(0.4,0.5,0.2,0.3,0.06,62);
	  t1->AddText(Form("%s%%",cent_Name[i]));
	  t1->Draw();
	  for(int k=0; k<nPtBins; k++)
	    {
	      double value = fabs(htmp->GetBinContent(k+1)-1);
	      if(max[i][k]<value) max[i][k] = value;
	    }
	}

      if(i==0)
	{
	  double xmin = 0.2, xmax = 0.5, ymin = 0.65, ymax = 0.88;
	  leg1 = new TLegend(xmin,ymin,xmax,ymax);
	  leg1->SetBorderSize(0);
	  leg1->SetFillColor(0);
	  leg1->SetTextSize(0.045);
	  leg1->AddEntry(hYield[i][0],"Default with stat. err.","P");
	  if(type==0)
	    {
	      leg1->AddEntry(hYield[i][1],"dca < 1.2 cm","P");
	      leg1->AddEntry(hYield[i][2],"dca < 0.8 cm","P");
	    }
	  else if(type==1)
	    {
	      leg1->AddEntry(hYield[i][1],"-1.25 < n#sigma_{#pi} < 3.25","P");
	      leg1->AddEntry(hYield[i][2],"-0.75 < n#sigma_{#pi} < 2.75","P");
	      leg1->AddEntry(hYield[i][3],"|dz| < 2.75#sigma","P");
	    }
	  leg1->Draw();

	  leg2 = new TLegend(xmin+0.4,ymin,xmax+0.4,ymax);
	  leg2->SetBorderSize(0);
	  leg2->SetFillColor(0);
	  leg2->SetTextSize(0.045);
	  if(type==0)
	    {
	      leg2->AddEntry(hYield[i][3],"NHitsFit >= 15","P");
	      leg2->AddEntry(hYield[i][4],"NHitsDedx >= 20","P");
	    }
	  else if(type==1)
	    {
	      leg2->AddEntry(hYield[i][4],"|dz| < 2.25#sigma","P");
	      leg2->AddEntry(hYield[i][5],"|dy| < 2.75#sigma","P");
	      leg2->AddEntry(hYield[i][6],"|dy| < 2.25#sigma","P");
	      leg2->AddEntry(hYield[i][7],"dtof < 1.5 ns","P");
	    }
	  
	  leg2->Draw();
	}
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,max[0][bin-1]+1);
	  htmp->SetBinContent(bin,1-max[0][bin-1]);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_%s.pdf",run_type,run_cfg_name.Data(),type_name[type]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/%sSys_%s.png",run_type,run_cfg_name.Data(),type_name[type]));
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"update");
      for(int i=0; i<nCentBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
      fout->Close();
    }
  
}

//================================================
void processSysCuts(int saveHisto = 0)
{
  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.sys.Cuts.root",run_config,pt1_cut,pt2_cut),"update");
  const int nSys = 12;
  const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp",
			    "nSigmaPiUp","nSigmaPiDown",
			    "dzUp","dzDown","dyUp","dyDown","dtofUp"};
  TString embed_name[2] = {"MCinput","MTDreco"};

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  // weight input distribution
  TFile *fWeight = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");
  TH1F *funcWeight[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      funcWeight[k] = (TH1F*)fWeight->Get(Form("TBW_Jpsi_Yield_cent%s",cent_Title[k]));
    }  


  TH1F *hSeUL[nSys][nCentBins][nPtBins];
  TH1F *hSeLS[nSys][nCentBins][nPtBins];
  TH1F *hMixUL[nSys][nCentBins][nPtBins];
  TH1F *hMixLS[nSys][nCentBins][nPtBins];
  TH1F *hMean[nSys][nCentBins];
  TH1F *hSigma[nSys][nCentBins];
  TH1F *hJpsiYield[nSys][nCentBins];

  TH1F *hJpsi[nSys][nCentBins][2];
  TH1F *hJpsiEff[nSys][nCentBins];

  for(int s=0; s<nSys; s++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hMean[s][k] = 0;
	  hSigma[s][k] = 0;
	  hJpsiYield[s][k] = 0;
	  hJpsiEff[s][k] = 0;
	}
    }



  // Extract raw yield
  double scale_low = 2.7, scale_high = 3.5;
  double fit_min = 2.6, fit_max = 3.8;
  int rebin_signal = 5;

  for(int s=0; s<nSys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);
      for(int k=0; k<nCentBins; k++)
	{
	  hMean[s][k] = new TH1F(Form("Sys_%s_Mean_cent%s",name[s],cent_Title[k]),Form("%s: fitted mean vs p_{T} (%s%%);p_{T} (GeV/c);mean",name[s],cent_Name[k]),nbins,xbins);
	  hSigma[s][k] = new TH1F(Form("Sys_%s_Sigma_cent%s",name[s],cent_Title[k]),Form("%s: fitted sigma vs p_{T} (%s%%);p_{T} (GeV/c);#sigma",name[s],cent_Name[k]),nbins,xbins);
	  hJpsiYield[s][k] = new TH1F(Form("Sys_%s_JpsiYield_cent%s",name[s],cent_Title[k]),Form("%s: Jpsi counts vs p_{T} (%s%%);p_{T} (GeV/c);Counts",name[s],cent_Name[k]),nbins,xbins);

	  // Get the mean and sigma from 0-10% centrality bin
	  TH1F *hFixMean, *hFixSigma;
	  if(k>0)
	    {
	      hFixMean = (TH1F*)fin->Get(Form("Sys_%s_Mean_cent%s",name[s],cent_Title[0]));
	      hFixSigma = (TH1F*)fin->Get(Form("Sys_%s_Sigma_cent%s",name[s],cent_Title[0]));
	    }

	  if(k==0)
	    {
	      TCanvas *c = new TCanvas(Form("Sys_%s",name[s]),Form("Sys_%s",name[s]),1100,700);
	      c->Divide(3,2);

	      TCanvas *cMix = new TCanvas(Form("Sys_%s_mix",name[s]),Form("Sys_%s_mix",name[s]),1100,700);
	      cMix->Divide(3,2);
	    }
	  for(Int_t i=1; i<nPtBins; i++)
	    {
	      hSeUL[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hSeLS[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hMixUL[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_Mix_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hMixLS[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_Mix_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));

	      // mix event background
	      double se = 0, se_err = 0, me = 0, me_err = 0;
	      int low_bin = hSeLS[s][k][i]->FindFixBin(scale_low+1e-4);
	      int high_bin = hSeLS[s][k][i]->FindFixBin(scale_high-1e-4);
	      for(int bin=low_bin; bin<=high_bin; bin++)
		{
		  se += hSeLS[s][k][i]->GetBinContent(bin);
		  se_err += TMath::Power(hSeLS[s][k][i]->GetBinError(bin),2);
		  me += hMixLS[s][k][i]->GetBinContent(bin);
		  me_err += TMath::Power(hMixLS[s][k][i]->GetBinError(bin),2);
		}
	      se_err = TMath::Sqrt(se_err);
	      me_err = TMath::Sqrt(me_err);
	      double scale = se/me;
	      double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);
	      TH1F *hMixBkg = (TH1F*)hMixUL[s][k][i]->Clone(Form("Sys_%s_mix_bkg_pt%s_cent%s",name[s],pt_Name[i],cent_Name[k]));
	      hMixBkg->Scale(scale);
	      hMixBkg->Rebin(rebin_signal);
	      cout << scale << endl;
	      if(k==0)
		{
		  cMix->cd(i);
		  TH1F *hScale = (TH1F*)hSeLS[s][k][i]->Clone(Form("%s_scale",hSeLS[s][k][i]->GetName()));
		  hScale->Divide(hMixLS[s][k][i]);
		  hScale->GetXaxis()->SetRangeUser(2.5,4);
		  hScale->SetMarkerStyle(21);
		  hScale->Draw();
		}

	      // signal
	      TH1F *hSignal = (TH1F*)hSeUL[s][k][i]->Clone(Form("Sys_%s_Jpsi_Signal_cent%s_pt%s",name[s],cent_Title[k],pt_Name[i]));
	      hSignal->Rebin(rebin_signal);
	      for(int bin=1; bin<=hSignal->GetNbinsX(); bin++)
		{
		  if(hSignal->GetBinContent(bin)==0)
		    {
		      hSignal->SetBinContent(bin,0);
		      hSignal->SetBinError(bin,1.4);
		    }
		}
	      hSignal->Add(hMixBkg,-1);
	      hSignal->SetLineColor(1);
	      hSignal->SetMarkerColor(1);

	      TF1 *funcSignal = new TF1(Form("Fit%s",hSignal->GetName()),"gausn(0)+pol1(3)",fit_min,fit_max);
	      funcSignal->SetParameter(1,3.09);
	      funcSignal->SetParameter(2,0.04);
	      if(k>0)
		{
		  funcSignal->FixParameter(1,hFixMean->GetBinContent(i));
		  funcSignal->FixParameter(2,hFixSigma->GetBinContent(i));
		}
	      funcSignal->SetLineColor(2);
	      TFitResultPtr ptr = hSignal->Fit(funcSignal,"IRS0Q");
	      hMean[s][k]->SetBinContent(i,funcSignal->GetParameter(1));
	      hMean[s][k]->SetBinError(i,funcSignal->GetParError(1));
	      hSigma[s][k]->SetBinContent(i,funcSignal->GetParameter(2));
	      hSigma[s][k]->SetBinError(i,funcSignal->GetParError(2));
	      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();

	      // bin counting
	      double low_mass_tmp = low_mass;
	      double high_mass_tmp = high_mass;
	      if(i<4) 
		{
		  low_mass_tmp = 2.95;
		  high_mass_tmp = 3.25;
		}
	      int low_bin = hSignal->FindFixBin(low_mass_tmp+1e-4);
	      int high_bin =hSignal->FindFixBin(high_mass_tmp-1e-4);
	      double error;
	      double count = hSignal->IntegralAndError(low_bin,high_bin,error);
	      TF1 *funcbkg = new TF1(Form("bkg_%s",hSignal->GetName()),"pol1",fit_min,fit_max);
	      const int nPar = 2;
	      double bkg_params[nPar];
	      for(int j=0; j<nPar; j++)
		{
		  funcbkg->SetParameter(j,funcSignal->GetParameter(3+j));
		  funcbkg->SetParError(j,funcSignal->GetParError(3+j));
		  bkg_params[j] = funcSignal->GetParameter(3+j);
		}
	      double bkg_matrix[nPar*nPar];
	      for(int j=3; j<3+nPar; j++)
		{
		  for(int l=3; l<3+nPar; l++)
		    {
		      bkg_matrix[(j-3)*nPar+l-3] = matrix[j*(3+nPar)+l];
		    }
		}
	      double bkg = funcbkg->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal->GetBinWidth(1);
	      double bkg_err = funcbkg->IntegralError(low_mass_tmp,high_mass_tmp,bkg_params,bkg_matrix) * 1./hSignal->GetBinWidth(1);
	      double signal = count - bkg;
	      double sig_err = TMath::Sqrt(error*error+bkg_err*bkg_err);
	      hJpsiYield[s][k]->SetBinContent(i,signal);
	      hJpsiYield[s][k]->SetBinError(i,sig_err);

	      // plotting
	      if(k==0)
		{
		  c->cd(i);	  
		  hSignal->SetTitle("");
		  hSignal->GetXaxis()->SetRangeUser(2.5,4);
		  hSignal->SetMarkerStyle(21);
		  hSignal->SetMaximum(2*hSignal->GetMaximum());
		  hSignal->Draw();
		  funcSignal->Draw("sames");
		  funcbkg->SetLineColor(4);
		  funcbkg->Draw("sames");
		  TPaveText *t = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c (%s%%)",ptBins_low[i],ptBins_high[i],cent_Name[k]),0.06);
		  t->Draw();
		  TLine *line = GetLine(low_mass_tmp,hSignal->GetMinimum()*1.5,low_mass_tmp,hSignal->GetMaximum()*0.3,1);
		  line->Draw();
		  TLine *line = GetLine(high_mass_tmp,hSignal->GetMinimum()*1.5,high_mass_tmp,hSignal->GetMaximum()*0.3,1);
		  line->Draw();

		  t = GetPaveText(0.16,0.3,0.7,0.88,0.045);
		  t->SetTextFont(62);
		  t->AddText("Counting");
		  t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",signal,sig_err,signal/sig_err));
		  t->SetTextAlign(11);
		  t->Draw();
		}
	    }
	  if(k==0)
	    {
	      TCanvas *c = new TCanvas(Form("Sys_%s_result",name[s]),Form("Sys_%s_result",name[s]),1100,700);
	      c->Divide(2,2);
	      c->cd(1);
	      hJpsiYield[s][k]->Draw();
	      c->cd(2);
	      hMean[s][k]->Draw();
	      c->cd(3);
	      hSigma[s][k]->Draw();
	    }

	  // +++++ embedding efficiency +++++

	  // Jpsi distribution
	  TH1F *hJpsiTmp = 0;
	  for(int i=0; i<2; i++)
	    {
	      if(k<3) hJpsiTmp = (TH1F*)fin->Get(Form("Sys_%s_Jpsi_pt_%s_cent%s",name[s],embed_name[i].Data(),cent_Title[k]));
	      else    hJpsiTmp = (TH1F*)fin->Get(Form("Sys_%s_Jpsi_pt_%s_cent%s",name[s],embed_name[i].Data(),cent_Title[2]));
	      hJpsiTmp->SetName(Form("Sys_%s_Jpsi_pt_%s_cent%s_tmp",name[s],embed_name[i].Data(),cent_Title[k]));
	      hJpsiTmp->Sumw2();
	      if(i==0)
		{
		  for(int bin=1; bin<=hJpsiTmp->GetNbinsX(); bin++)
		    {
		      double weight = funcWeight[k]->GetBinContent(funcWeight[k]->FindBin(hJpsiTmp->GetBinCenter(bin)));
		      hJpsiTmp->SetBinContent(bin,hJpsiTmp->GetBinContent(bin)*weight);
		      hJpsiTmp->SetBinError(bin,hJpsiTmp->GetBinError(bin)*weight);
		    }
		}
	      else if(i==1)
		{
		  if(k<3) TH2F *h2tmp = (TH2F*)fin->Get(Form("Sys_%s_JpsiPt_TrueVsReco_cent%s",name[s],cent_Title[k]));
		  else    TH2F *h2tmp = (TH2F*)fin->Get(Form("Sys_%s_JpsiPt_TrueVsReco_cent%s",name[s],cent_Title[2]));
		  h2tmp->SetName(Form("Sys_%s_JpsiPt_TrueVsReco_cent%s_tmp",name[s],cent_Title[k]));

		  for(int biny=1; biny<=h2tmp->GetNbinsY(); biny++)
		    {
		      double weight = funcWeight[k]->GetBinContent(funcWeight[k]->FindBin(h2tmp->GetYaxis()->GetBinCenter(biny)));
		      for(int binx=1; binx<=h2tmp->GetNbinsX(); binx++)
			{
			  h2tmp->SetBinContent(binx,biny,weight*h2tmp->GetBinContent(binx,biny));
			  h2tmp->SetBinError(binx,biny,weight*h2tmp->GetBinError(binx,biny));
			}
		    }
		  hJpsiTmp = (TH1F*)h2tmp->ProjectionX(Form("%s_WeightPt",hJpsiTmp->GetName()));
		}

	      hJpsi[s][k][i] = (TH1F*)hJpsiTmp->Rebin(nbins,Form("Sys_%s_Jpsi_pt_%s_cent%s_WeightPt_Rebin",name[s],embed_name[i].Data(),cent_Title[k]),xbins);
	    }

	  // Jpsi efficiency
	  hJpsiEff[s][k] = (TH1F*)hJpsi[s][k][1]->Clone(Form("Sys_%s_JpsiEff_cent%s",name[s],cent_Title[k]));
	  hJpsiEff[s][k]->Divide(hJpsi[s][k][0]);
	  if(k==0)
	    {
	      c->cd(4);
	      hJpsiEff[s][k]->Draw();
	    }
	}
    }

  if(saveHisto)
    {
      fin->cd();
      for(int s=0; s<nSys; s++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      if(hMean[s][k]) hMean[s][k]->Write("",TObject::kOverwrite);
	      if(hSigma[s][k]) hSigma[s][k]->Write("",TObject::kOverwrite);
	      if(hJpsiYield[s][k]) hJpsiYield[s][k]->Write("",TObject::kOverwrite);
	      if(hJpsiEff[s][k]) hJpsiEff[s][k]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void makeHistoSysCuts()
{
  const int nSys = 12;
  const char *name[nSys] = {"dcaUp","dcaDown","default","NHitsUp","NDedxUp",
			    "nSigmaPiUp","nSigmaPiDown",
			    "dzUp","dzDown","dyUp","dyDown","dtofUp"};

  const char *hName[3] = {"hJpsiInfo","hBkgLSPos","hBkgLSNeg"};
  THnSparseF *hnInvMass[nSys][3];
  TH1F *hInvMass[nSys][nCentBins][nPtBins][3];
  TH1F *hMixInvMass[nSys][nCentBins][nPtBins][3];
  TH3D *hMixMmumuvsPtCen[nSys][3];

  TH1F *hJpsi[nSys][nCentBins][2];
  TH2F *hJpsiPtMatch[nSys][nCentBins];
  TString embed_name[2] = {"MCinput","MTDreco"};

  int nsys = nSys;
  for(int s=0; s<nsys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);

      // same event
      TFile *fdata = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.root",name[s]),"read");
      for(Int_t j=0; j<3; j++)
	{
	  hnInvMass[s][j] = (THnSparseF*)fdata->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
	  hnInvMass[s][j]->SetName(Form("Sys_%s_%s",name[s],hnInvMass[s][j]->GetName()));
	  hnInvMass[s][j]->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[s][j]->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hnInvMass[s][j]->GetAxis(1)->SetRangeUser(ptBins_low[i]+0.01,ptBins_high[i]-0.01);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnInvMass[s][j]->GetAxis(7)->SetRange(centBins_low[k],centBins_high[k]);
		  hInvMass[s][k][i][j] = (TH1F*)hnInvMass[s][j]->Projection(0);
		  hInvMass[s][k][i][j]->SetName(Form("Sys%d_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k));
		  hInvMass[s][k][i][j]->Sumw2();
		  hnInvMass[s][j]->GetAxis(7)->SetRange(0,-1);
		}
	      hnInvMass[s][j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hnInvMass[s][j]->GetAxis(4)->SetRange(0,-1);
	  hnInvMass[s][j]->GetAxis(5)->SetRange(0,-1);
	}

      for(int k=0; k<nCentBins; k++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hInvMass[s][k][i][1]->Add(hInvMass[s][k][i][2]);
	    }
	}

      // mixed event
      TFile *fmix = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.mix.root",name[s]),"read");
      hMixMmumuvsPtCen[s][0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
      hMixMmumuvsPtCen[s][1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
      hMixMmumuvsPtCen[s][2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");
      for(Int_t j=0; j<3; j++)
	{
	  hMixMmumuvsPtCen[s][j]->SetName(Form("Sys%d_%s",s,hMixMmumuvsPtCen[s][j]->GetName()));
	  for(int i=0; i<nPtBins; i++)
	    {
	      int ybin_min = hMixMmumuvsPtCen[s][j]->GetYaxis()->FindFixBin(ptBins_low[i]+1e-4);
	      int ybin_max = hMixMmumuvsPtCen[s][j]->GetYaxis()->FindFixBin(ptBins_high[i]-1e-4);
	      for(int k=0; k<nCentBins; k++)
		{
		  TH1F *htmp = (TH1F*)hMixMmumuvsPtCen[s][j]->ProjectionZ(Form("Sys%d_mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",s,hName[j],trigName[kTrigType],i,k),centBins_low[k],centBins_high[k],ybin_min,ybin_max);
		  hMixInvMass[s][k][i][j] = new TH1F(Form("Sys%d_mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),2800,0,14);
		  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		    {
		      hMixInvMass[s][k][i][j]->SetBinContent(bin,htmp->GetBinContent(bin));
		      hMixInvMass[s][k][i][j]->SetBinError(bin,htmp->GetBinError(bin));
		    }
		}
	    }
	}

      for(int k=0; k<nCentBins; k++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hMixInvMass[s][k][i][1]->Add(hMixInvMass[s][k][i][2]);
	    }
	}

      // embeding
      TFile *fEmbed = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.embed.root",name[s]),"read");
      THnSparseF *hJpsiUS_mc = (THnSparseF*)fEmbed->Get("hJpsiInfo_di_mu");
      hJpsiUS_mc->SetName(Form("Sys%d_%s",s,hJpsiUS_mc->GetName()));
      
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiUS_mc->GetAxis(10)->SetRange(centBins_low[k],centBins_high[k]);
	  for(int i=0; i<2; i++)
	    {
	      hJpsiUS_mc->GetAxis(7)->SetRange(4+i*2,4+i*2);
	      if(i>0)
		{
		  hJpsiUS_mc->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
		  hJpsiUS_mc->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);
		  hJpsiUS_mc->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
		}
	  
	      hJpsi[s][k][i] = (TH1F*)hJpsiUS_mc->Projection(1);
	      hJpsi[s][k][i]->SetName(Form("Sys_%s_Jpsi_pt_%s_cent%s",name[s],embed_name[i].Data(),cent_Title[k]));
	      hJpsiUS_mc->GetAxis(0)->SetRange(0,-1);
	      hJpsiUS_mc->GetAxis(4)->SetRange(0,-1);
	      hJpsiUS_mc->GetAxis(5)->SetRange(0,-1);
	      hJpsiUS_mc->GetAxis(7)->SetRange(0,-1);
	    }
	  hJpsiUS_mc->GetAxis(10)->SetRange(0,-1);
	}

      // matched Jpsi for response
      THnSparseF *hnJpsiMatch = (THnSparseF*)fEmbed->Get("mhJpsiMatch_di_mu");
      hnJpsiMatch->SetName(Form("Sys%d_%s",s,hnJpsiMatch->GetName()));
      hnJpsiMatch->GetAxis(2)->SetRangeUser(low_mass+0.001, high_mass-0.001);
      hnJpsiMatch->GetAxis(4)->SetRangeUser(pt1_cut+0.01,100);
      hnJpsiMatch->GetAxis(5)->SetRangeUser(pt2_cut+0.01,100);

      for(int k=0; k<nCentBins; k++)
	{
	  hnJpsiMatch->GetAxis(6)->SetRange(centBins_low[k],centBins_high[k]);
	  hJpsiPtMatch[s][k] = (TH2F*)hnJpsiMatch->Projection(1,3);
	  hJpsiPtMatch[s][k]->Sumw2();
	  hJpsiPtMatch[s][k]->SetName(Form("Sys_%s_JpsiPt_TrueVsReco_cent%s",name[s],cent_Title[k]));
	}
      hnJpsiMatch->GetAxis(6)->SetRange(0,-1);
      hnJpsiMatch->GetAxis(2)->SetRange(0,-1);
      hnJpsiMatch->GetAxis(4)->SetRange(0,-1);
      hnJpsiMatch->GetAxis(5)->SetRange(0,-1);
    }

  TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.sys.Cuts.root",run_config,pt1_cut,pt2_cut),"recreate");
  for(int s=0; s<nsys; s++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hInvMass[s][k][i][0]->Write(Form("Sys_%s_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]),TObject::kOverwrite);
	      hInvMass[s][k][i][1]->Write(Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]),TObject::kOverwrite);
	      hMixInvMass[s][k][i][0]->Write(Form("Sys_%s_Mix_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]),TObject::kOverwrite);
	      hMixInvMass[s][k][i][1]->Write(Form("Sys_%s_Mix_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]),TObject::kOverwrite);
	    }
	}
    }

  for(int s=0; s<nsys; s++)
    {
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiPtMatch[s][k]->SetTitle("");
	  hJpsiPtMatch[s][k]->Write("",TObject::kOverwrite);
	  for(int i=0; i<2; i++)
	    {
	      hJpsi[s][k][i]->SetTitle("");
	      hJpsi[s][k][i]->Write("",TObject::kOverwrite);
	    }
	}
    }
  fout->Close();
}
