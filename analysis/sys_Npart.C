TFile *f;
const int year = YEAR;
TString run_cfg_name;
TString outFileName;
const int gApplyWeight = 1;
const int nPtBins = 2;
const double ptBins_low[nPtBins]  = {0,5};
const double ptBins_high[nPtBins] = {15,15};
const char *pt_Name[nPtBins] = {"0-15","5-15"};
const int nCentBins = 6;
const int CentBins_low[nCentBins]  = {5,7,9, 11,13,15};
const int CentBins_high[nCentBins] = {6,8,10,12,14,16};
const char *Cent_Name[nCentBins] = {"50-60","40-50","30-40","20-30","10-20","0-10"};

//================================================
void sys_Npart()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  f = TFile::Open(Form("Rootfiles/%s.Npart.pt%1.1f.pt%1.1f.sys.root",run_type,pt1_cut,pt2_cut),"update");

  //MtdRespEff();
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
  const int nSys = 5;
  const char *name[nSys] = {"signalExt","TrackCuts","PidCuts","trigEff","respEff"};
  TList *list = new TList;
  TString legName[nSys+1] = {"Total","Signal extraction","Tracking efficiency","PID efficiency","Trigger efficiency","Response efficiency"};
  TH1F *hSysAll[nPtBins];
  TH1F *hSys[nPtBins][nSys];
  TH1F *hSysSub[nPtBins][nSys];
  for(int i=0; i<nPtBins; i++)
    {
      for(int s=0; s<nSys; s++)
	{
	  hSys[i][s] = (TH1F*)f->Get(Form("Sys_%s_pt%s",name[s],pt_Name[i]));
	  hSys[i][s]->SetLineWidth(2);
	  hSys[i][s]->GetXaxis()->SetLabelSize(0.05);
	  hSysSub[i][s] = (TH1F*)hSys[i][s]->Clone(Form("Sys_%s_pt%s_sub",name[s],pt_Name[i]));
	  for(int bin=1; bin<=hSysSub[i][s]->GetNbinsX(); bin++)
	    {
	      hSysSub[i][s]->SetBinContent(bin,hSysSub[i][s]->GetBinContent(bin)-1);
	    }
	}
    
      hSysAll[i] = (TH1F*)hSys[0][i]->Clone(Form("Sys_all_pt%s",pt_Name[i]));
      hSysAll[i]->Reset();
      for(int bin=1; bin<=hSysAll[i]->GetNbinsX(); bin++)
	{
	  double error = 0;
	  for(int s=0; s<nSys; s++)
	    {
	      error += TMath::Power(hSysSub[i][s]->GetBinContent(bin),2);
	    }
	  error = sqrt(error);
	  hSysAll[i]->SetBinContent(bin,error);
	}
      
      list->Clear();
      list->Add(hSysAll[i]);
      for(int s=0; s<nSys; s++)
	{
	  list->Add(hSysSub[i][s]);
	}
      c = drawHistos(list,Form("Sys_Npart_pt%s",pt_Name[i]),Form("Systematic uncertainty for J/#psi yield (p_{T} > %1.0f GeV/c);;Uncertainties",ptBins_low[i]),kFALSE,0.5,10,kTRUE,0,0.4,kFALSE,kTRUE,legName,kTRUE,"",0.15,0.45,0.55,0.85,kFALSE);
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Npart/Sys_All_pt%s.pdf",run_type,pt_Name[i]));
    }

  if(saveHisto)
    {
      f->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSysAll[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void signalExtraction(int savePlot = 0, int saveHisto = 0)
{
  TH1F *hSys[nPtBins];
  TFile *fdata = TFile::Open(Form("Rootfiles/%s.Npart.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  const int nSys = 7;
  const char *sys_name[nSys] = {"","_LargeScale","_SmallScale","_LargeFit","_SmallFit","_Rebin","_pol1"};
  TH1F *hSignal[nPtBins][nSys];

  double max[nPtBins][nCentBins];
  for(int i=0; i<nPtBins; i++)
    {
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  max[i][bin-1] = 0;
	}
    }

  for(int i=0; i<nPtBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_signalExt_pt%s",pt_Name[i]),Form("Systematic uncertainty for signal extraction (%s GeV/c)",pt_Name[i]),nCentBins,0,nCentBins);
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  hSys[i]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	}

      TCanvas *c = new TCanvas(Form("Sys_signalExt_pt%s",pt_Name[i]),Form("Sys_signalExt_pt%s",pt_Name[i]),800,600);
      for(int j=0; j<6; j++)
	{
	  TString tmpName = Form("pt%s%s%s%s",pt_Name[i],gWeightName[gApplyWeight],gTrgSetupName[0],sys_name[j]);
	  if(j==0) hSignal[i][j] = (TH1F*)fdata->Get(Form("Jpsi_BinCountYield_%s",tmpName.Data()));
	  else     hSignal[i][j] = (TH1F*)f->Get(Form("Jpsi_BinCountYield_%s",tmpName.Data()));
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
	  htmp->SetTitle(";;Relative difference");
	  c->cd(i+1);
	  SetPadMargin(gPad,0.1,0.12,0.1,0.1);
	  htmp->GetXaxis()->SetLabelSize(0.05);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("Systematic uncertainty of signal extraction (p_{T,J/#Psi} > %1.0f GeV/c)",ptBins_low[i]),0.04);
	  t1->Draw();
	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      double value = fabs(htmp->GetBinContent(bin)-1);
	      if(max[i][bin-1]<value) max[i][bin-1] = value;
	    }
	}

      double xmin = 0.2, xmax = 0.5, ymin = 0.65, ymax = 0.88;
      leg1 = new TLegend(xmin,ymin,xmax,ymax);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.04);
      leg1->AddEntry(hSignal[i][0],"Default with stat. err.","P");
      leg1->AddEntry(hSignal[i][1],"Larger scale range","P");
      leg1->AddEntry(hSignal[i][2],"Smaller scale range","P");
      leg1->Draw();

      leg2 = new TLegend(xmin+0.4,ymin,xmax+0.4,ymax);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.04);
      leg2->AddEntry(hSignal[i][3],"Larger fit range","P");
      leg2->AddEntry(hSignal[i][4],"Smaller fit range","P");
      leg2->AddEntry(hSignal[i][5],"Different binning","P");
      //leg2->AddEntry(hSignal[i][6],"Background shape","P");
      leg2->Draw();
      
      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,max[i][bin-1]+1);
	  htmp->SetBinContent(bin,1-max[i][bin-1]);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
      if(savePlot)   c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Npart/Sys_signalExt_pt%s.pdf",run_type,pt_Name[i]));
    }


  if(saveHisto)
    {
      f->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
  
}


//================================================
void triggerEfficiency(int saveHisto = 1)
{
  TH1F *hSys[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_trigEff_pt%s",pt_Name[i]),Form("Systematic uncertainty for trigger efficiency (%s GeV/c)",pt_Name[i]),nCentBins,0,nCentBins);
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  hSys[i]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	  hSys[i]->SetBinContent(bin,0.1+1);
	}
    }

  if(saveHisto)
    {
      f->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void MtdRespEff(int saveHisto = 1)
{
  TH1F *hSys[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_respEff_pt%s",pt_Name[i]),Form("Systematic uncertainty for response efficiency (%s GeV/c)",pt_Name[i]),nCentBins,0,nCentBins);
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  hSys[i]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	  hSys[i]->SetBinContent(bin,0.03+1);
	}
    }

  if(saveHisto)
    {
      f->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void trackCuts(int savePlot = 0, int saveHisto = 0)
{
  const int type = 1; // 0 - quality cuts; 1 - PID cuts
  const char *type_name[2] = {"TrackCuts","PidCuts"};
  const char *type_title[2] = {"track quality cuts","track PID cuts"};

  if(type==0)
    {
      const int nSys = 5;
      const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp"};
    }
  else if(type==1)
    {
      const int nSys = 7;
      //const char *name[nSys] = {"default","nSigmaPiUp","nSigmaPiDown", "dzUp","dzDown","dyUp","dyDown","dtofUp"};
      const char *name[nSys] = {"default","nSigmaPiUp","nSigmaPiDown", "dzUp","dzDown","dyUp","dyDown"};
    }


  TH1F *hSys[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hSys[i] = new TH1F(Form("Sys_%s_pt%s",type_name[type],pt_Name[i]),Form("Systematic uncertainty for %s (%s)",type_title[type],pt_Name[i]),nCentBins,0,nCentBins);
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  hSys[i]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	}
    }

  TH1F *hYield[nPtBins][nSys];
  TH1F *hEff[nPtBins][nSys];

  double max[nPtBins][nCentBins];
  for(int i=0; i<nPtBins; i++)
    {
      for(int bin=1; bin<=nCentBins; bin++)
	{
	  max[i][bin-1] = 0;
	}
    }

  for(int i=0; i<nPtBins; i++)
    {
      TCanvas *c = new TCanvas(Form("Sys_trackCuts_pt%s",pt_Name[i]),Form("Sys_trackCuts_pt%s",pt_Name[i]),800,600);
      for(int j=0; j<nSys; j++)
	{
	  hYield[i][j] = (TH1F*)f->Get(Form("Sys_%s_JpsiYield_pt%s",name[j],pt_Name[i]));
	  hYield[i][j]->SetTitle(";;");
	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      hYield[i][j]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
	    }
	  hEff[i][j] = (TH1F*)f->Get(Form("Sys_%s_JpsiEff_pt%s",name[j],pt_Name[i]));
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
	  htmp->SetTitle(";;Relative difference");
	  SetPadMargin(gPad,0.1,0.12,0.1,0.1);
	  htmp->GetXaxis()->SetLabelSize(0.05);
	  if(j==0) htmp->Draw("P");
	  else htmp->Draw("P sames");
	  TPaveText *t1 = GetTitleText(Form("Systematic uncertainty of %s (p_{T,J/#Psi} > %1.0f GeV/c)",type_title[type],ptBins_low[i]),0.04);
	  t1->Draw();

	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      double value = fabs(htmp->GetBinContent(bin)-1);
	      if(max[i][bin-1]<value) max[i][bin-1] = value;
	    }
	}


      double xmin = 0.2, xmax = 0.5, ymin = 0.65, ymax = 0.88;
      leg1 = new TLegend(xmin,ymin,xmax,ymax);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.04);
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
      leg2->SetTextSize(0.04);
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

      TH1F *htmp = (TH1F*)hSys[i]->Clone(Form("%s_low",hSys[i]->GetName()));
      for(int bin=1; bin<=hSys[i]->GetNbinsX(); bin++)
	{
	  hSys[i]->SetBinContent(bin,max[i][bin-1]+1);
	  htmp->SetBinContent(bin,1-max[i][bin-1]);
	}
      hSys[i]->SetLineColor(2);
      hSys[i]->Draw("samesHIST");
      htmp->SetLineColor(2);
      htmp->Draw("samesHIST");
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/sys_Npart/Sys_%s_pt%s.pdf",run_type,type_name[type],pt_Name[i]));
    }


  if(saveHisto)
    {
      f->cd();
      for(int i=0; i<nPtBins; i++)
	{
	  hSys[i]->Write("",TObject::kOverwrite);
	}
    }
  
}

//================================================
void processSysCuts(int saveHisto = 0)
{
  const int nSys = 12;
  const char *name[nSys] = {"default","dcaUp","dcaDown","NHitsUp","NDedxUp",
			    "nSigmaPiUp","nSigmaPiDown",
			    "dzUp","dzDown","dyUp","dyDown","dtofUp"};

  // global setup
  double g_mix_scale_low = 2.6;
  double g_mix_scale_high = 3.6;
  double g_bin_width = 0.04; // 40 MeV
  TString g_func1 = "pol1", g_func2 = "pol3";
  int g_func1_npar = 2, g_func2_npar = 4;
  double g_sig_fit_min = 2.6;
  double g_sig_fit_max = 3.8;


  TH1F *hSeUL[nSys][nCentBins][nPtBins];
  TH1F *hSeLS[nSys][nCentBins][nPtBins];
  TH1F *hMixUL[nSys][nCentBins][nPtBins];
  TH1F *hMixLS[nSys][nCentBins][nPtBins];
  TH1F *hJpsiYield[nSys][nPtBins];
  TH1F *hJpsiEff[nSys][nPtBins];


  // Extract raw yield
  for(int s=0; s<nSys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiYield[s][i] = new TH1F(Form("Sys_%s_JpsiYield_pt%s",name[s],pt_Name[i]),"",nCentBins,0,nCentBins);

	  TCanvas *cSignal = new TCanvas(Form("Sys_%s_InvMass_pt%s",name[s],pt_Name[i]),Form("Sys_%s_InvMass_pt%s",name[s],pt_Name[i]),1100,650);
	  cSignal->Divide(3,2);

	  for(int k=nCentBins-1; k>-1; k--)
	    {
	      hSeUL[s][k][i] = (TH1F*)f->Get(Form("Sys_%s_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]));
	      hSeLS[s][k][i] = (TH1F*)f->Get(Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]));
	      hMixUL[s][k][i] = (TH1F*)f->Get(Form("Sys_%s_Mix_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]));
	      hMixLS[s][k][i] = (TH1F*)f->Get(Form("Sys_%s_Mix_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]));
	      for(int bin=1; bin<=hSeUL[s][k][i]->GetNbinsX(); bin++)
		{
		  double content = hSeUL[s][k][i]->GetBinContent(bin);
		  if(content==0)
		    {
		      hSeUL[s][k][i]->SetBinContent(bin,0);
		      hSeUL[s][k][i]->SetBinError(bin,1.4);
		    }
		  else
		    hSeUL[s][k][i]->SetBinError(bin,TMath::Sqrt(hSeUL[s][k][i]->GetBinContent(bin)));
		  hSeLS[s][k][i]->SetBinError(bin,TMath::Sqrt(hSeLS[s][k][i]->GetBinContent(bin)));
		}

	      // mix event background
	      double se = 0, se_err = 0, me = 0, me_err = 0;
	      int low_bin = hSeLS[s][k][i]->FindFixBin(g_mix_scale_low+1e-4);
	      int high_bin = hSeLS[s][k][i]->FindFixBin(g_mix_scale_high-1e-4);
	      for(int bin=low_bin; bin<=high_bin; bin++)
		{
		  se += hSeLS[s][k][i]->GetBinContent(bin);
		  se_err += TMath::Power(hSeLS[s][k][i]->GetBinError(bin),2);
		}
	      
	      int low_bin_me = hMixLS[s][k][i]->FindFixBin(g_mix_scale_low+1e-4);
	      int high_bin_me = hMixLS[s][k][i]->FindFixBin(g_mix_scale_high-1e-4);
	      for(int bin=low_bin_me; bin<=high_bin_me; bin++)
		{
		  me += hMixLS[s][k][i]->GetBinContent(bin);
		  me_err += TMath::Power(hMixLS[s][k][i]->GetBinError(bin),2);
		}
	      se_err = TMath::Sqrt(se_err);
	      me_err = TMath::Sqrt(me_err);
	      double scale = se/me;
	      double scale_error = scale * TMath::Sqrt(se_err*se_err/se/se+me_err*me_err/me/me);
	      TH1F *hMixBkg = (TH1F*)hMixUL[s][k][i]->Clone(Form("Sys_%s_mix_bkg_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]));
	      hMixBkg->Scale(scale);
	      hMixBkg->Rebin(g_bin_width/hMixBkg->GetBinWidth(1));

	      // signal
	      TH1F *hSignal = (TH1F*)hSeUL[s][k][i]->Clone(Form("Sys_%s_Jpsi_Signal_cent%s_pt%s",name[s],Cent_Name[k],pt_Name[i]));
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

	      TString funcForm;
	      int nPar;
	      if(i==0 && k>2)
		{
		  funcForm = g_func2;
		  nPar = g_func2_npar;
		}
	      else
		{
		  funcForm = g_func1;
		  nPar = g_func1_npar;
		}

	      TF1 *funcSignal  = new TF1(Form("Fit%s",hSignal->GetName()),Form("gausn(0)+%s(3)",funcForm.Data()),g_sig_fit_min,g_sig_fit_max);
	      funcSignal->SetParameter(0,100);
	      funcSignal->SetParameter(1,3.09);
	      funcSignal->SetParameter(2,0.1);
	      funcSignal->SetLineColor(2);
	      TFitResultPtr ptr = hSignal->Fit(funcSignal,"IRS0Q");
	      double *matrix = ptr->GetCovarianceMatrix().GetMatrixArray();

	      // bin counting
	      double low_mass_tmp = low_mass;
	      double high_mass_tmp = high_mass;
	      if(i==0) 
		{
		  low_mass_tmp = 2.96;
		  high_mass_tmp = 3.24;
		}
	      int low_bin = hSignal->FindFixBin(low_mass_tmp+1e-4);
	      int high_bin =hSignal->FindFixBin(high_mass_tmp-1e-4);
	      double error;
	      double count = hSignal->IntegralAndError(low_bin,high_bin,error);
	      TF1 *funcbkg = new TF1(Form("bkg_%s",hSignal->GetName()),funcForm.Data(),g_sig_fit_min,g_sig_fit_max);
	      const int nParameter = nPar;
	      double bkg_params[nParameter];
	      double bkg_matrix[nParameter*nParameter];
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
	      hJpsiYield[s][i]->SetBinContent(k+1,signal);
	      hJpsiYield[s][i]->SetBinError(k+1,sig_err);

	      // plotting
	      cSignal->cd(nCentBins-k);	  
	      hSignal->SetTitle("");
	      hSignal->GetXaxis()->SetRangeUser(2.5,4);
	      hSignal->SetMarkerStyle(21);
	      hSignal->SetMaximum(2*hSignal->GetMaximum());
	      hSignal->Draw();
	      funcSignal->Draw("sames");
	      funcbkg->SetLineColor(4);
	      funcbkg->Draw("sames");
	      TPaveText *t = GetTitleText(Form("p_{T} > %1.1f GeV/c (%s%%)",ptBins_low[i],Cent_Name[k]),0.06);;
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

  if(saveHisto)
    {
      f->cd();
      for(int s=0; s<nSys; s++)
	{
	  for(int i=0; i<nPtBins; i++)
	    {
	      hJpsiYield[s][i]->Write("",TObject::kOverwrite);
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
  TH1F *hInvMass[nSys][nPtBins][nCentBins][3];

  TH3D *hMixMmumuvsPtCen[nSys][3];
  TH1F *hMixInvMass[nSys][nPtBins][nCentBins][3];

  int nsys = nSys;
  for(int s=0; s<nsys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);

      // same event
      TFile *fdata = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.root",name[s]),"read");
      for(Int_t j=0; j<3; j++)
	{
	  hnInvMass[s][j] = (THnSparseF*)fdata->Get(Form("m%s_%s",hName[j],trigName[kTrigType]));
	  hnInvMass[s][j]->Sumw2();
	  hnInvMass[s][j]->SetName(Form("Sys_%s_%s",name[s],hnInvMass[s][j]->GetName()));
	  hnInvMass[s][j]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	  hnInvMass[s][j]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);

	  for(Int_t i=0; i<nPtBins; i++)
	    {
	      hnInvMass[s][j]->GetAxis(1)->SetRangeUser(ptBins_low[i]+0.01,ptBins_high[i]-0.01);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnInvMass[s][j]->GetAxis(5)->SetRange(CentBins_low[k],CentBins_high[k]);
		  hInvMass[s][i][k][j] = (TH1F*)hnInvMass[s][j]->Projection(0);
		  hInvMass[s][i][k][j]->SetName(Form("Sys%d_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k));
		  hInvMass[s][i][k][j]->Sumw2();
		  hnInvMass[s][j]->GetAxis(5)->SetRange(0,-1);
		}
	      hnInvMass[s][j]->GetAxis(1)->SetRange(0,-1);
	    }
	  hnInvMass[s][j]->GetAxis(3)->SetRange(0,-1);
	  hnInvMass[s][j]->GetAxis(4)->SetRange(0,-1);
	}

      for(Int_t i=0; i<nPtBins; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      hInvMass[s][i][k][1]->Add(hInvMass[s][i][k][2]);
	    }
	}

      // mixed event
      TFile *fmix = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.mix.root",name[s]),"read");
      hMixMmumuvsPtCen[s][0] = (TH3D*)fmix->Get("hMixULMmumuvsPtCen");
      hMixMmumuvsPtCen[s][1] = (TH3D*)fmix->Get("hMixLPosMmumuvsPtCen");
      hMixMmumuvsPtCen[s][2] = (TH3D*)fmix->Get("hMixLNegMmumuvsPtCen");
      for(Int_t j=0; j<3; j++)
	{
	  hMixMmumuvsPtCen[s][j]->Sumw2();
	  hMixMmumuvsPtCen[s][j]->SetName(Form("Sys%d_%s",s,hMixMmumuvsPtCen[s][j]->GetName()));
	  for(int i=0; i<nPtBins; i++)
	    {
	      int ybin_min = hMixMmumuvsPtCen[s][j]->GetYaxis()->FindFixBin(ptBins_low[i]+1e-4);
	      int ybin_max = hMixMmumuvsPtCen[s][j]->GetYaxis()->FindFixBin(ptBins_high[i]-1e-4);
	      for(int k=0; k<nCentBins; k++)
		{
		  TH1F *htmp = (TH1F*)hMixMmumuvsPtCen[s][j]->ProjectionZ(Form("Sys%d_mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d_tmp",s,hName[j],trigName[kTrigType],i,k),CentBins_low[k],CentBins_high[k],ybin_min,ybin_max);
		  hMixInvMass[s][i][k][j] = new TH1F(Form("Sys%d_mix_%s_%s_InvMass_jpsi_PtBin%d_CentBin%d",s,hName[j],trigName[kTrigType],i,k),htmp->GetTitle(),1400,0,14);
		  for(int bin=1; bin<=htmp->GetNbinsX(); bin++)
		    {
		      hMixInvMass[s][i][k][j]->SetBinContent(bin,htmp->GetBinContent(bin));
		      hMixInvMass[s][i][k][j]->SetBinError(bin,htmp->GetBinError(bin));
		    }
		}
	    }
	}

      for(Int_t i=0; i<nPtBins; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      hMixInvMass[s][i][k][1]->Add(hMixInvMass[s][i][k][2]);
	    }
	}
    }

  // embedding
  // Get Jpsi counts as weights
  TFile *fYield = TFile::Open(Form("Rootfiles/Pico.Run14.AuAu200.jpsi.pt%1.1f.pt%1.1f.yield.root",pt1_cut,pt2_cut),"read");
  TH1F *hJpsiWeight[gNTrgSetup-1];
  double nJpsi[nCentBins-1][gNTrgSetup-1];
  double nJpsiCent[nCentBins-1];
  for(int k=0; k<nCentBins-1; k++)
    {
      nJpsiCent[k] = 0;
    }
  for(int i=0; i<gNTrgSetup-1; i++)
    {
      hJpsiWeight[i] = (TH1F*)fYield->Get(Form("NJpsiInCent_weight%s",gTrgSetupName[i+1]));
      for(int bin=1; bin<=hJpsiWeight[i]->GetNbinsX(); bin++)
	{
	  nJpsi[bin-1][i] = hJpsiWeight[i]->GetBinContent(bin);
	  nJpsiCent[bin-1] += hJpsiWeight[i]->GetBinContent(bin);
	}
    }
  for(int i=0; i<gNTrgSetup-1; i++)
    for(int bin=1; bin<=hJpsiWeight[i]->GetNbinsX(); bin++)
      printf("[i] %s%% %s: %4.2f Jpsi\n",cent_Title[bin],gTrgSetupTitle[i+1],nJpsi[bin-1][i]/nJpsiCent[bin-1]);

  // Get correction factor for othe centralities
  TFile *fCorr = TFile::Open(Form("Rootfiles/%s.Npart.pt%1.1f.pt%1.1f.root",run_type,pt1_cut,pt2_cut),"read");
  TH1F *hJpsiEffCorr[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      hJpsiEffCorr[i] = (TH1F*)fCorr->Get(Form("hJpsiEffCor_pt%s",pt_Name[i]));
    }

  THnSparseF *hnJpsiInfo[2];
  TH1F *hJpsiCounts[nSys][2][nPtBins][gNTrgSetup];
  TH1F *hJpsiEff[nSys][nPtBins][gNTrgSetup];
  TH1F *hJpsiEffFinal[nSys][nPtBins];

  for(int s=0; s<nsys; s++)
    {
      printf("+++ Process %s +++\n",name[s]);
      // embeding
      TFile *fEmbed = TFile::Open(Form("output/Sys.%s.Run14.AuAu200.jpsi.embed.root",name[s]),"read");
      cout << fEmbed->GetName() << endl;
      hnJpsiInfo[0] = (THnSparseF*)fEmbed->Get("hJpsiInfo_MC_di_mu_w");
      hnJpsiInfo[1] = (THnSparseF*)fEmbed->Get("hJpsiInfo_MtdTrig_di_mu_w");

      for(int j=0; j<2; j++)
	{
	  for(int t=0; t<gNTrgSetup; t++)
	    {
	      for(int i=0; i<nPtBins; i++)
		{
		  hJpsiCounts[s][j][i][t] = new TH1F(Form("Sys%d_JpsiCounts_PtBin%d_TrgBin%d_%d",s,i,t,j),"",nCentBins,0,nCentBins);
		  for(int bin=1; bin<=nCentBins; bin++)
		    {
		      hJpsiCounts[s][j][i][t]->GetXaxis()->SetBinLabel(bin,Form("%s%%",Cent_Name[bin-1]));
		      hJpsiCounts[s][j][i][t]->GetXaxis()->SetLabelSize(0.05);
		    }
		}
	    }
	}

      for(int j=0; j<2; j++)
	{
	  hnJpsiInfo[j]->SetName(Form("Sys%d_%s",s,hnJpsiInfo[j]->GetName()));
	  if(j>0)
	    {
	      hnJpsiInfo[j]->GetAxis(3)->SetRangeUser(pt1_cut+0.01,100);
	      hnJpsiInfo[j]->GetAxis(4)->SetRangeUser(pt2_cut+0.01,100);
	    }
	  for(int t=0; t<gNTrgSetup; t++)
	    {
	      if(t>0) hnJpsiInfo[j]->GetAxis(6)->SetRange(t,t);
	      for(int k=0; k<nCentBins; k++)
		{
		  hnJpsiInfo[j]->GetAxis(5)->SetRange(CentBins_low[k],CentBins_high[k]);
		  TH1F *htmp = (TH1F*)hnJpsiInfo[j]->Projection(1);
		  htmp->SetName(Form("hJpsiPt_%d%d%d%d",s,j,t,k));
		  htmp->SetTitle("");
		  for(int i=0; i<nPtBins; i++)
		    {
		      int low_bin = htmp->FindFixBin(ptBins_low[i]+1e-4);
		      int hi_bin  = htmp->FindFixBin(ptBins_high[i]-1e-4);
		      double err;
		      double yield = htmp->IntegralAndError(low_bin,hi_bin,err);
		      hJpsiCounts[s][j][i][t]->SetBinContent(k+1,yield);
		      hJpsiCounts[s][j][i][t]->SetBinError(k+1,err);
		    }
		  hnJpsiInfo[j]->GetAxis(5)->SetRange(0,-1);
		}
	      hnJpsiInfo[j]->GetAxis(6)->SetRange(0,-1);
	    }
	  hnJpsiInfo[j]->GetAxis(3)->SetRange(0,-1);
	  hnJpsiInfo[j]->GetAxis(4)->SetRange(0,-1);
	}

      for(int i=0; i<nPtBins; i++)
	{
	  for(int t=0; t<gNTrgSetup; t++)
	    {
	      hJpsiEff[s][i][t] = (TH1F*)hJpsiCounts[s][1][i][t]->Clone(Form("Sys%d_JpsiEff_PtBin%d_TrgBin%d_%d",s,i,t,j));
	      hJpsiEff[s][i][t]->Divide(hJpsiCounts[s][0][i][t]);
	    }
	}
      // Get average efficiency
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiEffFinal[s][i] = (TH1F*)hJpsiEff[s][i][0]->Clone(Form("Sys_%s_JpsiEff_pt%s",name[s],pt_Name[i]));
	  hJpsiEffFinal[s][i]->Reset();
	  // work on 0-10% first
	  double avg_eff = 0, avg_err = 0;
	  for(int t=0; t<gNTrgSetup-1; t++)
	    {
	      double eff = hJpsiEff[s][i][t+1]->GetBinContent(nCentBins);
	      double err = hJpsiEff[s][i][t+1]->GetBinError(nCentBins);
	      avg_eff += 1/eff * nJpsi[0][t]/nJpsiCent[0];
	      avg_err += TMath::Power(nJpsi[0][t]/nJpsiCent[0]*err/eff,2);
	    }
	  avg_eff = 1/avg_eff;
	  avg_err = avg_eff * TMath::Sqrt(avg_err);
	  hJpsiEffFinal[s][i]->SetBinContent(nCentBins,avg_eff);
	  hJpsiEffFinal[s][i]->SetBinError(nCentBins,avg_err);

	  // corrected for other centralities
	  for(int bin=1; bin<=nCentBins; bin++)
	    {
	      hJpsiEffFinal[s][i]->SetBinContent(bin,hJpsiEffFinal[s][i]->GetBinContent(nCentBins)*hJpsiEffCorr[i]->GetBinContent(bin));
	      hJpsiEffFinal[s][i]->SetBinError(bin,hJpsiEffFinal[s][i]->GetBinError(nCentBins)*hJpsiEffCorr[i]->GetBinContent(bin));
	    }
	}
    }


  f->cd();
  for(int s=0; s<nsys; s++)
    {
      for(Int_t i=0; i<nPtBins; i++)
	{
	  for(int k=0; k<nCentBins; k++)
	    {
	      hInvMass[s][i][k][0]->Write(Form("Sys_%s_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]),TObject::kOverwrite);
	      hInvMass[s][i][k][1]->Write(Form("Sys_%s_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]),TObject::kOverwrite);
	      hMixInvMass[s][i][k][0]->Write(Form("Sys_%s_Mix_InvMass_UL_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]),TObject::kOverwrite);
	      hMixInvMass[s][i][k][1]->Write(Form("Sys_%s_Mix_InvMass_LS_pt%s_cent%s",name[s],pt_Name[i],Cent_Name[k]),TObject::kOverwrite);
	    }
	}
    }

  for(int s=0; s<nsys; s++)
    {
      for(Int_t i=0; i<nPtBins; i++)
	{
	  hJpsiEffFinal[s][i]->Write("",TObject::kOverwrite);
	}
    }
}
