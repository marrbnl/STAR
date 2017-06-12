TFile *f;
const int year = YEAR;
TString file_name;

//================================================
void sys_JpsiXsec()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //qa();
  signalExtraction();
  //makeHistoSysCuts();
  //processSysCuts();
  //trackCuts();
  //triggerEfficiency();
  //MtdRespEff();
  //mergeSystematics();
}


//================================================
void mergeSystematics(const int icent = 0, int savePlot = 1, int saveHisto = 1)
{
  TFile *fout = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.systematics.root",run_config,pt1_cut,pt2_cut),"update");
  const int nSys = 5;
  const char *name[nSys] = {"signalExt","TrackCuts","PidCuts","trigEff","respEff"};
  TList *list = new TList;
  TString legName[nSys+1] = {"Total","Signal extraction","Tracking efficiency","PID efficiency","Trigger efficiency","Response efficiency"};
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
void signalExtraction(int savePlot = 0, int saveHisto = 0)
{
 // re-assign global constants
  const int nPtBins         = nPtBins_pt;
  const double* ptBins_low  = ptBins_low_pt;
  const double* ptBins_high = ptBins_high_pt;
  const char** pt_Name      = pt_Name_pt;
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

  TH1F *hSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_signalExt_%s",cent_Title[i]),Form("Systematic uncertainty for signal extraction (%s%%)",cent_Name[i]),nbins,xbins);
    }

  const int nSys = 10;
  const char *sys_name[nSys]  = {"","_LargeScale","_SmallScale","_ScaleFit","_Binning","_BkgFunc","_LargeFit","_SmallFit","_SigFunc","_FixSig"};
  const char *sys_leg[nSys] = {"Default","Larger bkg norm range","Smaller bkg norm range","Fit ME/SE w/ pol1","Binning","Res. bkg function","Larger sig fit range","Smaller sig fit range","Crystal-ball","Fix sig. mean&sigma"};

  TString outName    = Form("Rootfiles/%s.JpsiYield.pt%1.1f.pt%1.1f.%sroot",run_type,pt1_cut,pt2_cut,run_config);
  TString outNameSys = Form("Rootfiles/%s.Sys.JpsiYield.root",run_type);

  TFile *fdata = TFile::Open(Form("%s",outName.Data()),"read");
  TFile *fsys = 0x0;
  if(saveHisto) fsys = TFile::Open(Form("%s",outNameSys.Data()),"update");
  else          fsys = TFile::Open(Form("%s",outNameSys.Data()),"read");

  TH1F *hSignal[nCentBins][nSys];

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

  const char* method = "Fit";
  const int color_sys[nSys] = {1, 2, 3, 4, 6, 7, 8, 9, 28, 40};
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<9; j++)
	{
	  if(j==0) hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_%sYield_cent%s_weight%s",method,cent_Title[i],sys_name[j]));
	  else     hSignal[i][j] = (TH1F*)fsys->Get(Form("Jpsi_%sYield_cent%s_weight%s",method,cent_Title[i],sys_name[j]));
	  hSignal[i][j]->SetMarkerSize(1);
	  hSignal[i][j]->SetMarkerStyle(21);
	  hSignal[i][j]->SetMarkerColor(color_sys[j]);
	  
	  TH1F *htmp = (TH1F*)hSignal[i][j]->Clone(Form("%s_tmp",hSignal[i][j]->GetName()));
	  htmp->Divide(hSignal[i][0]);
	  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
	    {
	      if(j==0) htmp->SetBinError(bin,hSignal[i][j]->GetBinError(bin)/hSignal[i][j]->GetBinContent(bin));
	      else     htmp->SetBinError(bin,0);
	    }
	  htmp->GetYaxis()->SetRangeUser(0.6,1.7);
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
	  TLegend *leg[2];
	  for(int l=0; l<2; l++)
	    {
	      leg[l] = new TLegend(xmin+0.4*l,ymin,xmax+0.4*l,ymax);
	      leg[l]->SetBorderSize(0);
	      leg[l]->SetFillColor(0);
	      leg[l]->SetTextSize(0.045);
	    }
	  for(int j=0; j<9; j++)
	    leg[j/5]->AddEntry(hSignal[i][j],sys_leg[j],"P");
	  leg[0]->Draw();
	  leg[1]->Draw();
	}
      
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  double sys = max[0][bin-1];
	  if(bin>=2) sys = 0.06;
	  hSys[i]->SetBinContent(bin,sys+1);
	  htmp->SetBinContent(bin,1-sys);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_JpsiYield/Sys_signalExt.pdf",run_type));
    }

  if(saveHisto)
    {
      fsys->cd();
      for(int i=0; i<nCentBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
  
}

//================================================
void MtdRespEff(int saveHisto = 1)
{
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TH1F *hSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_respEff_%s",cent_Title[i]),Form("Systematic uncertainty for trigger efficiency (%s%%)",cent_Name[i]),nbins,xbins);
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,0.03+1);
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
	  hSys[i]->SetBinContent(bin,0.1+1);
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
  const int type = 1; // 0 - quality cuts; 1 - PID cuts
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
      //const int nSys = 8;
      //const char *name[nSys] = {"default","nSigmaPiUp","nSigmaPiDown", "dzUp","dzDown","dyUp","dyDown","dtofUp"};
      const int nSys = 7;
      const char *name[nSys] = {"default","nSigmaPiUp","nSigmaPiDown", "dzUp","dzDown","dyUp","dyDown"};
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
	      //leg2->AddEntry(hYield[i][7],"dtof < 1.5 ns","P");
	    }
	  
	  leg2->Draw();
	}
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      double avg_sys = 0;
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  double sys = max[0][bin-1];
	  if(type==0 && bin>=4) sys = 0.078;
	  if(type==1 && bin>=5) sys = 0.08;
	  hSys[i]->SetBinContent(bin,sys+1);
	  htmp->SetBinContent(bin,1-sys);
	  //avg_sys += max[0][bin-1];
	}
      //avg_sys /= hSys[i]->GetNbinsX();
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
      //TLine *line = GetLine(0,avg_sys+1,ptBins_high[nPtBins-1],avg_sys+1,1);
      //line->Draw();
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
void processSysCuts(int saveHisto = 1)
{
  TFile *fin = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.sys.Cuts.root",run_config,pt1_cut,pt2_cut),"update");
  const int nSys = 12;
  const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp",
			    "nSigmaPiUp","nSigmaPiDown",
			    "dzUp","dzDown","dyUp","dyDown","dtofUp"};

  // global setup
  double g_mix_scale_low = 2.6;
  double g_mix_scale_high = 4;
  double g_bin_width = 0.04; // 40 MeV
  TString g_func1 = "pol3";
  TString g_func2 = "pol0";
  int g_func1_npar = 4;
  int g_func2_npar = 1;
  double g_sig_fit_min = 2.5;
  double g_sig_fit_max = 4.0;

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];


  TH1F *hSeUL[nSys][nCentBins][nPtBins];
  TH1F *hSeLS[nSys][nCentBins][nPtBins];
  TH1F *hMixUL[nSys][nCentBins][nPtBins];
  TH1F *hMixLS[nSys][nCentBins][nPtBins];
  TH1F *hMean[nSys][nCentBins];
  TH1F *hSigma[nSys][nCentBins];
  TH1F *hJpsiYield[nSys][nCentBins];

  // Extract raw yield

  for(int s=0; s<nSys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);
      for(int k=0; k<nCentBins; k++)
	{
	  hMean[s][k] = new TH1F(Form("Sys_%s_Mean_cent%s",name[s],cent_Title[k]),Form("%s: fitted mean vs p_{T} (%s%%);p_{T} (GeV/c);mean",name[s],cent_Name[k]),nbins,xbins);
	  hSigma[s][k] = new TH1F(Form("Sys_%s_Sigma_cent%s",name[s],cent_Title[k]),Form("%s: fitted sigma vs p_{T} (%s%%);p_{T} (GeV/c);#sigma",name[s],cent_Name[k]),nbins,xbins);
	  hJpsiYield[s][k] = new TH1F(Form("Sys_%s_JpsiYield_cent%s",name[s],cent_Title[k]),Form("%s: Jpsi counts vs p_{T} (%s%%);p_{T} (GeV/c);Counts",name[s],cent_Name[k]),nbins,xbins);

	  TCanvas *c = new TCanvas(Form("Sys_%s_%s",name[s],cent_Title[k]),Form("Sys_%s_%s",name[s],cent_Title[k]),1100,700);
	  c->Divide(5,2);

	  //TCanvas *cmix = new TCanvas(Form("Sys_mix_%s_%s",name[s],cent_Title[k]),Form("Sys_mix_%s_%s",name[s],cent_Title[k]),1100,700);
	  //cmix->Divide(5,2);

	  TString funcForm;
	  int nPar;
	  for(Int_t i=1; i<nPtBins; i++)
	    {
	      hSeUL[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hSeLS[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hMixUL[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_Mix_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));
	      hMixLS[s][k][i] = (TH1F*)fin->Get(Form("Sys_%s_Mix_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],cent_Title[k]));

	      // mix event background
	      double g_mix_scale_low_tmp = g_mix_scale_low;
	      double g_mix_scale_high_tmp = g_mix_scale_high;
	      if(i==1)
		{
		  //g_mix_scale_low_tmp = 2.5;
		  //g_mix_scale_high_tmp = 2.7;
		}

	      double se = 0, se_err = 0, me = 0, me_err = 0;
	      int low_bin = hSeLS[s][k][i]->FindFixBin(g_mix_scale_low_tmp+1e-4);
	      int high_bin = hSeLS[s][k][i]->FindFixBin(g_mix_scale_high_tmp-1e-4);
	      for(int bin=low_bin; bin<=high_bin; bin++)
		{
		  se += hSeLS[s][k][i]->GetBinContent(bin);
		  se_err += TMath::Power(hSeLS[s][k][i]->GetBinError(bin),2);
		}

	      int low_bin_me = hMixLS[s][k][i]->FindFixBin(g_mix_scale_low_tmp+1e-4);
	      int high_bin_me = hMixLS[s][k][i]->FindFixBin(g_mix_scale_high_tmp-1e-4);
	      for(int bin=low_bin_me; bin<=high_bin_me; bin++)
		{
		  me += hMixLS[s][k][i]->GetBinContent(bin);
		  me_err += TMath::Power(hMixLS[s][k][i]->GetBinError(bin),2);
		}
	      se_err = TMath::Sqrt(se_err);
	      me_err = TMath::Sqrt(me_err);
	      double scale = se/me;
	      double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);
	      cout << scale << endl;
	      TH1F *hMixBkg = (TH1F*)hMixUL[s][k][i]->Clone(Form("Sys_%s_mix_bkg_pt%s_cent%s",name[s],pt_Name[i],cent_Name[k]));
	      hMixBkg->Scale(scale);
	      // cmix->cd(i);
	      // gPad->SetLogy();
	      // hMixBkg->GetXaxis()->SetRangeUser(2,5);
	      // hMixBkg->DrawCopy();
	      // hSeLS[s][k][i]->SetMarkerStyle(20);
	      // hSeLS[s][k][i]->Draw("sames");

	      hMixBkg->Rebin(g_bin_width/hMixBkg->GetBinWidth(1));

	      // signal
	      TH1F *hSignal = (TH1F*)hSeUL[s][k][i]->Clone(Form("Sys_%s_Jpsi_Signal_cent%s_pt%s",name[s],cent_Title[k],pt_Name[i]));
	      hSignal->Rebin(g_bin_width/hSignal->GetBinWidth(1));
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

	      if(i<5)
		{
		  funcForm = g_func1; 
		  nPar = g_func1_npar;
		}
	      else
		{
		  funcForm = g_func2; 
		  nPar = g_func2_npar;
		}

	      if(i==1) g_sig_fit_min = 2.7;
	      else     g_sig_fit_min = 2.5;
	      TF1 *funcSignal  = new TF1(Form("Fit_%s",hSignal->GetName()),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	      funcSignal->SetParameter(1,3.09);
	      funcSignal->SetParameter(2,0.1);
	      if(i==5) funcSignal->SetParameter(2,0.05);
	      // Get the mean and sigma from 0-10% centrality bin
	      if(k>0)
		{
		  funcSignal->FixParameter(1,hMean[s][0]->GetBinContent(i));
		  funcSignal->FixParameter(2,hSigma[s][0]->GetBinContent(i));
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
	      if(i<5) 
		{
		  low_mass_tmp = 2.96;
		  high_mass_tmp = 3.24;
		}
	      int low_bin = hSignal->FindFixBin(low_mass_tmp+1e-4);
	      int high_bin =hSignal->FindFixBin(high_mass_tmp-1e-4);
	      double error;
	      double count = hSignal->IntegralAndError(low_bin,high_bin,error);

	      // from combined fit
	      const int nParameter = nPar;
	      double bkg_params[nParameter];
	      double bkg_matrix[nParameter*nParameter];
	      TF1 *funcbkg = new TF1(Form("Bkg_%s",hSignal->GetName()),Form("%s",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	      for(int j=0; j<nPar; j++)
		{
		  funcbkg->SetParameter(j,funcSignal->GetParameter(3+j));
		  funcbkg->SetParError(j,funcSignal->GetParError(3+j));
		  bkg_params[j] = funcSignal->GetParameter(3+j);
		}
 
	      for(int j=3; j<3+nPar; j++)
		{
		  for(int m=3; m<3+nPar; m++)
		    {
		      bkg_matrix[(j-3)*nPar+m-3] = matrix[j*(3+nPar)+m];
		    }
		}

	      double bkg = funcbkg->Integral(low_mass_tmp,high_mass_tmp) * 1./hSignal->GetBinWidth(1);
	      double bkg_err = funcbkg->IntegralError(low_mass_tmp,high_mass_tmp,bkg_params,bkg_matrix) * 1./hSignal->GetBinWidth(1);
	      double signal = count - bkg;
	      double sig_err = TMath::Sqrt(error*error+bkg_err*bkg_err);
	      hJpsiYield[s][k]->SetBinContent(i,signal);
	      hJpsiYield[s][k]->SetBinError(i,sig_err);
	      // plotting
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
    }

  
  // +++++ embedding efficiency +++++
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

  TFile *fCorr = TFile::Open(Form("Rootfiles/%s.JpsiEff.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TH1F *hEffCorr[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hEffCorr[k] = (TH1F*)fCorr->Get(Form("hJpsiEffCorr_cent%s",cent_Title[k]));
    }

  TString embed_name[2] = {"MC","MtdTrig"};
  TH1F *hJpsi[nSys][gNTrgSetup][nCentBins][2];
  TH1F *hJpsiOneOverEff[nSys][gNTrgSetup][nCentBins];
  TH1F *hJpsiPtEffAvg0020[nSys];
  TH1F *hJpsiEff[nSys][nCentBins];
 
  for(int s=0; s<nSys; s++)
    {
      for(int j=0; j<gNTrgSetup; j++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      for(int i=0; i<2; i++)
		{
		  //cout << Form("Sys_%s_Jpsi_pt_%s_cent%s_%s",name[s],embed_name[i].Data(),cent_Title[k],gTrgSetupTitle[j]) << endl;
		  TH1F *hJpsiTmp = (TH1F*)fin->Get(Form("Sys_%s_Jpsi_pt_%s_cent%s_%s",name[s],embed_name[i].Data(),cent_Title[k],gTrgSetupTitle[j]));
		  //cout << hJpsiTmp->GetName() << endl;
		  hJpsi[s][j][k][i] = (TH1F*)hJpsiTmp->Rebin(nbins,Form("%s_rebin",hJpsiTmp->GetName()),xbins);
		}

	      // Jpsi efficiency
	      hJpsiOneOverEff[s][j][k] = (TH1F*)hJpsi[s][j][k][0]->Clone(Form("Sys_%s_JpsiEff_cent%s%s",name[s],cent_Title[k],gTrgSetupTitle[j]));
	      hJpsiOneOverEff[s][j][k]->Divide(hJpsi[s][j][k][1]);
	    }
	}

      
      // average efficiency for 0-20%
      hJpsiPtEffAvg0020[s] = (TH1F*)hJpsiOneOverEff[s][0][0]->Clone(Form("Sys_%s_JpsiPtEff_cent%s",name[s],cent_Title[1]));
      hJpsiPtEffAvg0020[s]->Reset();
      for(int j=0; j<gNTrgSetup-1; j++)
	{
	  hJpsiPtEffAvg0020[s]->Add(hJpsiOneOverEff[s][j+1][1], nJpsi[0][j]/nJpsiCent[0]);
	}
      for(int bin=1; bin<=hJpsiPtEffAvg0020[s]->GetNbinsX(); bin++)
	{
	  double value = hJpsiPtEffAvg0020[s]->GetBinContent(bin);
	  double error = hJpsiPtEffAvg0020[s]->GetBinError(bin);
	  hJpsiPtEffAvg0020[s]->SetBinContent(bin, 1./value);
	  hJpsiPtEffAvg0020[s]->SetBinError(bin, 1./value * error/value);
	  if(bin==8) hJpsiPtEffAvg0020[s]->SetBinError(bin, 1./value * 0.01);
	  //printf("[i] pt = %2.1f with eff = %4.2e, err = %4.2f%%\n",hJpsiPtEffAvg0020[s]->GetBinCenter(bin),hJpsiPtEffAvg0020[s]->GetBinContent(bin),error/value*100);
	}
      //c = draw1D(hJpsiPtEffAvg0020[s]);

      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiEff[s][k] = (TH1F*)hJpsiPtEffAvg0020[s]->Clone(Form("Sys_%s_JpsiEff_cent%s",name[s], cent_Title[k]));
	  hJpsiEff[s][k]->Multiply(hEffCorr[k]);
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
	  hnInvMass[s][j]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[s][j]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);

	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hnInvMass[s][j]->GetAxis(1)->SetRangeUser(ptBins_low[i]+0.01,ptBins_high[i]-0.01);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnInvMass[s][j]->GetAxis(5)->SetRange(centBins_low[k],centBins_high[k]);
		  hInvMass[s][k][i][j] = (TH1F*)hnInvMass[s][j]->Projection(0);
		  hInvMass[s][k][i][j]->SetName(Form("Sys%d_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k));
		  hInvMass[s][k][i][j]->Sumw2();
		  hnInvMass[s][j]->GetAxis(5)->SetRange(0,-1);
		}
	      hnInvMass[s][j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hnInvMass[s][j]->GetAxis(3)->SetRange(0,-1);
	  hnInvMass[s][j]->GetAxis(4)->SetRange(0,-1);
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
		  hMixInvMass[s][k][i][j] = new TH1F(Form("Sys%d_mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),1400,0,14);
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
    }


  THnSparseF *hnJpsiInfo[2];
  TH1F *hJpsi[nSys][gNTrgSetup][nCentBins][2];
  TH2F *hJpsiPtMatch[nSys][nCentBins];
  TString embed_name[2] = {"MC","MtdTrig"};
  for(int s=0; s<nsys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);

      // embeding
      TFile *fEmbed = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.embed.root",name[s]),"read");
      for(int i=0; i<2; i++)
	{
	  hnJpsiInfo[i] =  (THnSparseF*)fEmbed->Get(Form("hJpsiInfo_%s_di_mu_w",embed_name[i].Data()));
	  hnJpsiInfo[i]->SetName(Form("Sys%d_%s",s,hnJpsiInfo[i]->GetName()));
	  if(i>0)
	    {
	      hnJpsiInfo[i]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	      hnJpsiInfo[i]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	    }

	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      if(j>0) hnJpsiInfo[i]->GetAxis(6)->SetRange(j,j);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnJpsiInfo[i]->GetAxis(5)->SetRange(centBins_low[k],centBins_high[k]);

		  hJpsi[s][j][k][i] = (TH1F*)hnJpsiInfo[i]->Projection(1);
		  hJpsi[s][j][k][i]->SetName(Form("Sys_%s_Jpsi_pt_%s_cent%s_%s",name[s],embed_name[i].Data(),cent_Title[k],gTrgSetupTitle[j]));
		  hJpsi[s][j][k][i]->SetTitle("");
		  hJpsi[s][j][k][i]->SetBinContent(hJpsi[s][j][k][i]->GetNbinsX()+1,0); // reset overflow bin
		  hJpsi[s][j][k][i]->Sumw2();

		  hnJpsiInfo[i]->GetAxis(5)->SetRange(0,-1); 
		}
	      hnJpsiInfo[i]->GetAxis(6)->SetRange(0,-1);
	    }
	  hnJpsiInfo[i]->GetAxis(3)->SetRange(0,-1);
	  hnJpsiInfo[i]->GetAxis(4)->SetRange(0,-1);
	}
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
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<gNTrgSetup; j++)
	    {
	      for(int k=0; k<nCentBins; k++)
		{
		  hJpsi[s][j][k][i]->Write("",TObject::kOverwrite);
		}
	    }
	}
    }
  fout->Close();
}
