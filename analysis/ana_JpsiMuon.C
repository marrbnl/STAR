const int year = YEAR;
TFile *f;

//================================================
void ana_JpsiMuon()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }

  //DeltaTof();
  MtdVpdTacDiff();
  //kink();
}

//================================================
void MtdVpdTacDiff(const Int_t savePlot = 1)
{
  TList *list = new TList;

  const int nHistos = 2;
  const char *name[2] = {"AuAu200","pp200"};
  const char *title[2] = {"Run14_AuAu_200","Run15_pp_200"};
  const int nbins = 6;
  const double xbins[nbins+1] = {1.2,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum[nHistos] = {760, 880};
  const double maximum[nHistos] = {840, 980};

  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];

  TFile *fdata[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      if(i==0) fdata[i] = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
      if(i==1) fdata[i] = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");
      hDataUL[i]       = (TH2F*)fdata[i]->Get("mhJpsiMuonMtdVpdTacDiff_UL_di_mu");
      hDataUL[i]->SetName(Form("%s_%d",hDataUL[i]->GetName(),i));
      hDataUL[i]->Sumw2();

      hDataLS[i]       = (TH2F*)fdata[i]->Get(Form("mhJpsiMuonMtdVpdTacDiff_LS_di_mu"));
      hDataLS[i]->SetName(Form("%s_%d",hDataLS[i]->GetName(),i));
      hDataLS[i]->Sumw2();

      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("JpsiMuonMtdVpdTacDiff_%d",i));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);
      
    }

  // UL vs LS
  const int rebin = 2;
  TH1F *hUL[nHistos][nbins];
  TH1F *hLS[nHistos][nbins];
  TH1F *hMuon[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->Rebin(rebin);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i],maximum[i]);
	  hUL[i][bin-1]->SetTitle("");
	  hUL[i][bin-1]->Draw("P");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);
	  hLS[i][bin-1]->Rebin(rebin);
	  hLS[i][bin-1]->Draw("samesP");

	  hMuon[i][bin-1] = (TH1F*)hUL[i][bin-1]->Clone(Form("%s_DataMtdVpdTacDiff_bin%d",name[i],bin));
	  hMuon[i][bin-1]->Add(hLS[i][bin-1],-1);

	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.5,0.6,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(title[i]);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_MtdVpdTacDiff_ULvsLS.pdf",run_type,name[i]));
    } 

  // Fit data to extract efficiency
  const double min[2] = {788,900};
  const double max[2] = {837,960};
  const double shift[2] = {0, 126};
  TH1F *hFitDataEff[nHistos];
  TH1F *hFitDataEfflow[nHistos];
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hFitDataMean[i] = new TH1F(Form("%s_JpsiMuon_MtdVpdTacDiff_FitMean",name[i]),Form("%s: mean of MtdVpdTacDiff;p_{T} (GeV/c)",title[i]),nbins,xbins);
      hFitDataSigma[i] = new TH1F(Form("%s_JpsiMuon_MtdVpdTacDiff_FitSigma",name[i]),Form("%s: sigma of MtdVpdTacDiff;p_{T} (GeV/c)",title[i]),nbins,xbins);

      TH1F *hBase = (TH1F*)hFitDataMean[i]->Clone(Form("hBase_%d",i));
      TH1F *hMatch = (TH1F*)hFitDataMean[i]->Clone(Form("hMatch_%d",i));
      TH1F *hMatch2 = (TH1F*)hBase->Clone(Form("hMatch2_%d",i));

      TCanvas *c = new TCanvas(Form("%s_FitMtdVpdTacDiff",name[i]),Form("%s_FitMtdVpdTacDiff",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  TF1 *func = new TF1(Form("func_%d_%d",i,bin),"gaus",min[i],maximum[i]);
	  func->SetParameter(2,5);
	  hFit->Fit(func,"IR0Q");
	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->Draw();
	  func->SetLineColor(4);
	  func->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  double all = func->Integral(minimum[i],maximum[i]);
	  double all_err = func->IntegralError(minimum[i],maximum[i]);
	  double acc = func->Integral(min[i],max[i]);
	  double acc_err = func->IntegralError(min[i],max[i]);
	  double acc2 = func->Integral(785,837);
	  double acc_err_2 = func->IntegralError(785,837);
	  hBase->SetBinContent(bin,all);
	  hBase->SetBinError(bin,all_err);
	  hMatch->SetBinContent(bin,acc);
	  hMatch->SetBinError(bin,acc_err);
	  hMatch2->SetBinContent(bin,acc2);
	  hMatch2->SetBinError(bin,acc_err_2);
	  hFitDataMean[i]->SetBinContent(bin,func->GetParameter(1)-shift[i]);
	  hFitDataMean[i]->SetBinError(bin,func->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,func->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,func->GetParError(2));
	}
      hFitDataEff[i] = DivideTH1ForEff(hMatch,hBase,Form("%s_JpsiMuon_MtdVpdTacDiff_FitEff_prodhigh",name[i]));
      hFitDataEfflow[i] = DivideTH1ForEff(hMatch2,hBase,Form("%s_JpsiMuon_MtdVpdTacDiff_FitEff_prodlow",name[i]));
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.8,0.85,0.055);
      t1->AddText(title[i]);
      t1->SetTextColor(2);
      t1->SetTextFont(62);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_FitMtdVpdTacDiff.pdf",run_type,name[i]));
    }

  TList *list = new TList;
  TString legName[2];
  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataMean[i]);
      if(shift[i]==0) legName[i] = title[i];
      else            legName[i] = Form("%s: - %1.0f ch",title[i],shift[i]);
    }
  c = drawHistos(list,"MtdVpdTacDiff_Mean","Mean of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);Mean",false,0,0,true,780,810,false,true,legName,true,"",0.15,0.4,0.2,0.4,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_MtdVpdTacDiffMean.pdf",run_type));
  list->Clear();

  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataSigma[i]);
      legName[i] = title[i];
    }
  c = drawHistos(list,"MtdVpdTacDiff_Sigma","Sigma of TAC_{MTD}-TAC_{VPD} distribution for J/#psi muons;p_{T} (GeV/c);#sigma",false,0,0,true,0,20,false,true,legName,true,"",0.15,0.3,0.73,0.88,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_MtdVpdTacDiffSigma.pdf",run_type));
  list->Clear();

  // compare pp vs AuAu
  TCanvas *c = new TCanvas(Form("Muon_pp_vs_AuAu"),Form("Muon_pp_vs_AuAu"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      c->cd(bin);
      TLegend *leg = new TLegend(0.15,0.7,0.6,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      for(int i=0; i<nHistos; i++)
	{
	  TH1F *htmp = (TH1F*)hMuon[i][bin-1]->Clone(Form("%s_clone",hMuon[i][bin-1]->GetName()));
	  htmp->Rebin(2);
	  htmp->SetMarkerStyle(21+i*4);
	  htmp->SetMarkerColor(i+1);
	  htmp->SetLineColor(i+1);
	  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	    {
	      int shift_bin = int(shift[i]/htmp->GetBinWidth(1));
	      if(shift_bin==0) continue;
	      int new_bin = ibin - shift_bin;
	      if(new_bin<=0) continue;
	      htmp->SetBinContent(new_bin, htmp->GetBinContent(ibin));
	      htmp->SetBinError(new_bin, htmp->GetBinError(ibin));
	      htmp->SetBinContent(ibin, 0);
	      htmp->SetBinError(ibin, 0);
	    }
	  htmp->Scale(1./htmp->GetBinContent(htmp->FindFixBin(hFitDataMean[0]->GetBinContent(bin))));
	  htmp->GetXaxis()->SetRangeUser(minimum[0],maximum[0]);
	  htmp->SetMaximum(1.6);
	  if(i==0) htmp->Draw();
	  else     htmp->Draw("samesP");
	  if(i==0) leg->AddEntry(htmp,title[i],"PL");
	  else     leg->AddEntry(htmp,Form("%s: - %1.0f ch",title[i],shift[i]),"PL");
	}
      if(bin==2)       leg->Draw();
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/CompMtdVpdTac_ppVsAuAu.pdf",run_type));

  // trigger efficiency
  TFile *fTrigEff = TFile::Open("Rootfiles/Run14.AuAu200.MuonTrigEff.root","read");
  const char *name[2] = {"prod_high","prod_low"};
  TH1F *hdata = 0x0;
  TF1 *fExt = 0x0;
  for(int i=0; i<2; i++)
    {
      if(i==0) 
	{
	  hdata = hFitDataEff[0];
	  fExt = (TF1*)fTrigEff->Get("MuonTrigEff_cent0060_P3");
	}
      else    
	{
	  hdata = hFitDataEfflow[0];
	  fExt = (TF1*)fTrigEff->Get("MuonTrigEff_cent0060_P1");
	}
      hdata->SetMarkerStyle(21);
      hdata->GetYaxis()->SetRangeUser(0.5,1.1);
      c = draw1D(hdata,"MTD trigger efficiency for single muons;p_{T} (GeV/c);efficiency");
      fExt->SetLineColor(4);
      fExt->Draw("sames");
      TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(Form("Run14_AuAu200, %s",name[i]));
      leg->AddEntry(hdata,"Data-driven: J/#Psi muons","P");
      leg->AddEntry(fExt,"Extrapolation method","L");
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/MtdTrigEff_%s.pdf",run_type,name[i]));
    }
}

//================================================
void DeltaTof(const Int_t savePlot = 0)
{
  TList *list = new TList;

  const int nHistos = 2;
  const char *name[2] = {"AuAu200","pp200"};
  const char *title[2] = {"Run14_AuAu_200","Run15_pp_200"};
  const int nbins = 6;
  const double xbins[nbins+1] = {1.2,1.5,2.0,2.5,3.0,5.0,10.0};

  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];
  TH2F *hDtofVsMod[nHistos];

  TFile *fdata[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      if(i==0) fdata[i] = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
      if(i==1) fdata[i] = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");
      hDataUL[i]       = (TH2F*)fdata[i]->Get("mhJpsiMuonDtof_UL_di_mu");
      hDataUL[i]->SetName(Form("%s_%d",hDataUL[i]->GetName(),i));
      hDataUL[i]->Sumw2();

      hDataLS[i]       = (TH2F*)fdata[i]->Get(Form("mhJpsiMuonDtof_LS_di_mu"));
      hDataLS[i]->SetName(Form("%s_%d",hDataLS[i]->GetName(),i));
      hDataLS[i]->Sumw2();

      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("JpsiMuonDtofVsPt_%d",i));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);

      hDtofVsMod[i] = (TH2F*)fdata[i]->Get(Form("mhDeltaTof_%s",trigName[kTrigType]));
      hDtofVsMod[i]->SetName(Form("%s_%d",hDtofVsMod[i]->GetName(),i));
      
    }

  for(int i=0; i<nHistos; i++)
    {
      hDtofVsMod[i]->GetYaxis()->SetRangeUser(-10,10);
      c = draw2D(hDtofVsMod[i],Form("%s: #Deltatof of tracks matched to MTD;Module;#Deltatof (ns)",title[i]));
      if(savePlot) 
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_DtofvsMod_MthTrk.pdf",run_type,name[i]));
    }

  // UL vs LS
  const int rebin = 5;
  TH1F *hUL[nHistos][nbins];
  TH1F *hLS[nHistos][nbins];
  TH1F *hMuon[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("%s_DataDtof_UL_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->Rebin(rebin);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(-2,4);
	  hUL[i][bin-1]->SetTitle("");
	  hUL[i][bin-1]->Draw("P");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("%s_DataDtof_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);
	  hLS[i][bin-1]->Rebin(rebin);
	  hLS[i][bin-1]->Draw("samesP");

	  hMuon[i][bin-1] = (TH1F*)hUL[i][bin-1]->Clone(Form("%s_DataDtof_bin%d",name[i],bin));
	  hMuon[i][bin-1]->Add(hLS[i][bin-1],-1);

	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.5,0.6,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(title[i]);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_Dtof_ULvsLS_InPtBins.pdf",run_type,name[i]));
    }

  // Fit data distributions
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  TF1 *func[nHistos][nbins];
  TFitResultPtr ptr[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      hFitDataMean[i] = new TH1F(Form("%s_JpsiMuon_Dtof_FitMean",name[i]),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",title[i]),nbins,xbins);
      hFitDataSigma[i] = new TH1F(Form("%s_JpsiMuon_Dtof_FitSigma",name[i]),Form("%s: sigma of #Deltatof;p_{T} (GeV/c)",title[i]),nbins,xbins);
      TCanvas *c = new TCanvas(Form("%s_FitDtof",name[i]),Form("%s_FitDtof",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  hFit->GetXaxis()->SetRangeUser(-3,3);
	  func[i][bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,3);
	  func[i][bin-1]->SetParameter(2,0.1);
	  ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"IR0QS");
	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->Draw();
	  func[i][bin-1]->SetLineColor(4);
	  func[i][bin-1]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  hFitDataMean[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(1));
	  hFitDataMean[i]->SetBinError(bin,func[i][bin-1]->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,func[i][bin-1]->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.8,0.85,0.055);
      t1->AddText(title[i]);
      t1->SetTextColor(2);
      t1->SetTextFont(62);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/%s_JpsiMuon_FitDtof.pdf",run_type,name[i]));
    }

  TList *list = new TList;
  TString legName1[2];
  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataMean[i]);
      legName1[i] = title[i];
    }
  c = drawHistos(list,"Dtof_Mean","Mean of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);Mean",false,0,0,true,-0.8,0.3,false,true,legName1,true,"",0.15,0.4,0.2,0.4,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_DtofMean.pdf",run_type));
  list->Clear();


  for(int i=0; i<nHistos; i++)
    {
      list->Add(hFitDataSigma[i]);
    }
  c = drawHistos(list,"Dtof_Sigma","Sigma of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);#sigma",false,0,0,true,0,0.5,false,true,legName1,true,"",0.15,0.3,0.73,0.88,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/Fit_DtofSigma.pdf",run_type));
  list->Clear();

  // extract efficiency
  const int nValue = 3;
  const double min[nValue] = {-4, -4, -4};
  const double max[nValue] = {1, 0.4, 0.2};
  TH1F *hFitDataEff[nHistos][nValue][2];
  for(int i=0; i<nHistos; i++)
    {
      for(int j=0; j<nValue; j++)
	{
	  // Fitting method
	  TH1F *hBase = new TH1F(Form("hBase_%d_%d",i,j),Form("hBase_%d_%d",i,j),nbins,xbins);
	  TH1F *hMatch = new TH1F(Form("hMatch_%d_%d",i,j),Form("hMatch_%d_%d",i,j),nbins,xbins);
	  for(int bin=1; bin<=nbins; bin++)
	    {
	      double all = func[i][bin-1]->Integral(-4,4);
	      double all_err = func[i][bin-1]->IntegralError(-4,4,func[i][bin-1]->GetParameters(),  ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	      double acc = func[i][bin-1]->Integral(min[j],max[j]);
	      double acc_err = func[i][bin-1]->IntegralError(min[j],max[j],func[i][bin-1]->GetParameters(),  ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	      hBase->SetBinContent(bin,all);
	      hBase->SetBinError(bin,all_err);
	      hMatch->SetBinContent(bin,acc);
	      hMatch->SetBinError(bin,acc_err);
	    }
	  hFitDataEff[i][j][0] = DivideTH1ForEff(hMatch,hBase,Form("%s_JpsiMuon_FitEff_Dtof%1.1f",name[i],max[j]));

	  // Counting method
	  hBase->Reset();
	  hMatch->Reset();
	  for(int bin=1; bin<=nbins; bin++)
	    {
	      double all_err;
	      int low_bin = hMuon[i][bin-1]->FindFixBin(-4);
	      int high_bin = hMuon[i][bin-1]->FindFixBin(4);
	      double all = hMuon[i][bin-1]->IntegralAndError(low_bin, high_bin,all_err);

	      double acc_err;
	      low_bin = hMuon[i][bin-1]->FindFixBin(min[j]);
	      high_bin = hMuon[i][bin-1]->FindFixBin(max[j]-1e-4);
	      double acc = hMuon[i][bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
	      hBase->SetBinContent(bin,all);
	      hBase->SetBinError(bin,all_err);
	      hMatch->SetBinContent(bin,acc);
	      hMatch->SetBinError(bin,acc_err);
	      cout << acc << " < " << all << endl;
	    }
	  hFitDataEff[i][j][1] = DivideTH1ForEff(hMatch,hBase,Form("%s_JpsiMuon_CountEff_Dtof%1.1f",name[i],max[j]));
	}
    }

  TString legName[4];
  for(int j=0; j<nValue; j++)
    {
      list->Clear();
      for(int i=0; i<nHistos; i++)
	{
	  for(int k=0; k<2; k++)
	    {
	      hFitDataEff[i][j][k]->SetMarkerStyle(21+k*4);
	      hFitDataEff[i][j][k]->SetMarkerColor(i+1);
	      hFitDataEff[i][j][k]->SetLineColor(i+1);
	      list->Add(hFitDataEff[i][j][k]);
	      if(k==0) legName[i*2+k] = Form("%s: fitting",title[i]);
	      else     legName[i*2+k] = Form("%s: counting",title[i]);
	    }
	}
      c = drawHistos(list,Form("Eff_Dtof%1.1f",max[j]),Form("Efficiency of #Deltatof < %1.1f ns cut for J/#Psi muons;p_{T} (GeV/c)",max[j]),kFALSE,0,10,true,0.2,1.2,kFALSE,kTRUE,legName,kTRUE,"",0.4,0.6,0.2,0.45,kTRUE,0.04,0.04,false,1,false,false);
      list->Clear();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/DataMuon_DtofEff_%1.1f.pdf",run_type,max[j]));
    }
}

//================================================
void kink(const Int_t save = 0)
{
  THnSparseF *hnKink = (THnSparseF*)f->Get(Form("mhMuonKink_%s",trigName[kTrigType]));
  hnKink->GetAxis(2)->SetRangeUser(1.6,100);

  // check invariant mass
  const char *name[3] = {"UL","LS","UL-LS"};
  TH1F *hInvMass[2];
  TH1F *hDr[2], *hDz[2];
  for(int i=0; i<2; i++)
    {
      hnKink->GetAxis(1)->SetRange(i+1,i+1);

      hnKink->GetAxis(0)->SetRangeUser(low_mass+0.001, high_mass-0.001);
      hDr[i] = (TH1F*)hnKink->Projection(3);
      hDr[i]->Sumw2();
      hDr[i]->SetName(Form("hDr_%s",name[i]));
      TH1F *htmp = (TH1F*)hnKink->Projection(6);
      htmp->Sumw2();
      htmp->SetName(Form("hDr_%s_2",name[i]));
      hDr[i]->Add(htmp);
      cout << hDr[i]->GetEntries() << endl;

      hDz[i] = (TH1F*)hnKink->Projection(4);
      hDz[i]->SetName(Form("hDz_%s",name[i]));
      hDz[i]->Sumw2();
      htmp = (TH1F*)hnKink->Projection(7);
      htmp->SetName(Form("hDz_%s_2",name[i]));
      htmp->Sumw2();
      hDz[i]->Add(htmp);
      hnKink->GetAxis(0)->SetRange(0,-1);

      hInvMass[i] = (TH1F*)hnKink->Projection(0);
      hInvMass[i]->SetName(Form("hInvMass_%s",name[i]));
      hInvMass[i]->Rebin(4);
      hInvMass[i]->SetMarkerStyle(21);
      hInvMass[i]->SetMarkerColor(color[1-i]);
      hInvMass[i]->SetLineColor(hInvMass[i]->GetMarkerColor());
    }
  hnKink->GetAxis(1)->SetRange(0,-1);
  c = draw1D(hInvMass[0],"Invariant mass distribution of dimuon pairs;M_{#mu#mu} (GeV/c^{2});Counts");
  hInvMass[1]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.15,0.62,0.3,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p_{T,1} > 1.5 GeV/c");
  leg->AddEntry(hInvMass[0],"Unlike-sign","P");
  leg->AddEntry(hInvMass[1],"Like-sign","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_InvMass.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_InvMass.pdf",run_type));
    }

  // analyze dr and dz distribution
  TH1F *hKinkDr[2], *hKinkDz[2];
  hKinkDr[0] = (TH1F*)hDr[0]->Clone("hKinkDr_US-LS");
  hKinkDr[0]->Add(hDr[1],-1);
  hKinkDz[0] = (TH1F*)hDz[0]->Clone("hKinkDz_US-LS");
  hKinkDz[0]->Add(hDz[1],-1);
  hnKink->GetAxis(1)->SetRange(1,1);
  hKinkDr[1] = (TH1F*)hnKink->Projection(3);
  hKinkDr[1]->Sumw2();
  hKinkDr[1]->SetName("hKinkDr_LS");
  htmp = (TH1F*)hnKink->Projection(6);
  htmp->Sumw2();
  htmp->SetName("hKinkDr_LS_2");
  hKinkDr[1]->Add(htmp);
  hKinkDz[1] = (TH1F*)hnKink->Projection(4);
  hKinkDz[1]->Sumw2();
  hKinkDz[1]->SetName("hKinkDz_LS");
  htmp = (TH1F*)hnKink->Projection(7);
  htmp->Sumw2();
  htmp->SetName("hKinkDz_LS_2");
  hKinkDz[1]->Add(htmp);
  hnKink->GetAxis(1)->SetRange(0,-1);
  for(int i=0; i<2; i++)
    {
      hKinkDr[i]->Rebin(2);
      hKinkDr[i]->Scale(1./hKinkDr[i]->GetBinContent(hKinkDr[i]->FindFixBin(0)));
      hKinkDr[i]->SetMarkerStyle(21);
      hKinkDr[i]->SetMarkerColor(color[1-i]);
      hKinkDr[i]->SetLineColor(hKinkDr[i]->GetMarkerColor());

      hKinkDz[i]->Scale(1./hKinkDz[i]->GetBinContent(hKinkDz[i]->FindFixBin(0)));
      hKinkDz[i]->SetMarkerStyle(21);
      hKinkDz[i]->SetMarkerColor(color[1-i]);
      hKinkDz[i]->SetLineColor(hKinkDz[i]->GetMarkerColor());

      cout << "Fraction = " << hKinkDz[i]->Integral(hKinkDz[i]->FindFixBin(-0.15),hKinkDz[i]->FindFixBin(0.15))/hKinkDz[i]->Integral() << endl;
    }
  c = draw1D(hKinkDr[0],"Distance in transverse plane for muon candidates;#sqrt{(#Deltax)^{2}+(#Deltay)^{2}} (cm)");
  hKinkDr[1]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.6,0.65,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hKinkDr[0],"J/psi muon","P");
  leg->AddEntry(hKinkDr[1],"All muon","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dr.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dr.pdf",run_type));
    }

  hKinkDz[0]->GetXaxis()->SetRangeUser(-2,2);
  c = draw1D(hKinkDz[0],"Distance in z coordinate for muon candidates;#Deltaz (cm)");
  hKinkDz[1]->Draw("sames HIST");
  TLegend *leg = new TLegend(0.6,0.65,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hKinkDz[0],"J/psi muon","P");
  leg->AddEntry(hKinkDz[1],"All muon","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dz.png",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiMuon/kink_Dz.pdf",run_type));
    }
}


