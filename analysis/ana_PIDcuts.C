const int year = YEAR;
TF1 *fResVsPt[3];
double ptbound = 100;
double nsigma1 = 2, nsigma2 = 2;

//================================================
void ana_PIDcuts()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
    }
  else if(year==2014)
    {
      ptbound = 3;
      nsigma2 = 2.5;
    }
  else if(year==2015)
    {
      run_type = "Run15_pp200";
      ptbound = 3;
      nsigma2 = 2.5;
    }

  fResVsPt[0]= new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[0]->SetParameters(-32.6793, 32.6034, 0.444217);
  fResVsPt[1] = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[1]->SetParameters(-17.6867, 18.4528, 0.637142);
  fResVsPt[2] = new TF1("fResDtofVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[2]->SetParameters(0.0817528, 0.0169419, 4.34897);


  //makeData();
  //anaEmbed();
  //anaData();
  nSigmaPi();
}

//================================================
void anaData(const bool savePlot = 1, const bool saveHisto = 0)
{
  const int year = 2014;
  TFile *fembed = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
  TFile *fdata = 0x0;
  if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
  if(year==2015) fdata = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");

  const int nHistos = 3;
  const char *name[nHistos] = {"Dz","Dy","NsigmaPi"};
  const char *title[nHistos] = {"#Deltaz","#Deltay","n#sigma_{#pi}"};
  const char *unit[nHistos] = {" (cm)"," (cm)",""};
  const int nbins = 6;
  const double xbins[nbins+1] = {1.2,1.5,2.0,2.5,3.0,5.0,10.0};
  const int rebin[nHistos] = {4,4,2};
  const double minimum[nHistos] = {-100, -50, -5};
  const double maximum[nHistos] = {100, 50, 5};

  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hInvMassVpPtUL[nHistos];
  TH2F *hInvMassVpPtLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];
  TH2F *hEmbedDisVsPt[nHistos];
  TH1F *hEmbedEff[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hDataUL[i]       = (TH2F*)fdata->Get(Form("mhJpsiMuon%s_UL_di_mu",name[i]));
      hDataUL[i]->Sumw2();
      hDataLS[i]       = (TH2F*)fdata->Get(Form("mhJpsiMuon%s_LS_di_mu",name[i]));
      hDataLS[i]->Sumw2();
      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("mhJpsiMuon%s_UL_di_mu_Clone",name[i]));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);

      hInvMassVpPtUL[i]       = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_UL_di_mu",name[i]));
      hInvMassVpPtUL[i]->Sumw2();
      hInvMassVpPtLS[i]       = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_LS_di_mu",name[i]));
      hInvMassVpPtLS[i]->Sumw2();

      hEmbedDisVsPt[i] = (TH2F*)fembed->Get(Form("Embed_JpsiMuon_%sVsPt",name[i]));
      hEmbedEff[i]     = (TH1F*)fembed->Get(Form("Embed_JpsiMuon_%s_Eff",name[i]));
    }

  /*
  // compare invariant mass
  TH1F *hInvMassUL[nHistos][nbins];
  TH1F *hInvMassLS[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_InvMass",name[i]),Form("%s_InvMass",name[i]),1100,700);
      c->Divide(3,2);
      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hInvMassVpPtUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hInvMassVpPtUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hInvMassUL[i][bin-1] = (TH1F*)hInvMassVpPtUL[i]->ProjectionY(Form("Data_InvMass_UL_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassUL[i][bin-1]->SetMarkerStyle(20);
	  hInvMassUL[i][bin-1]->SetMarkerStyle(20);
	  hInvMassUL[i][bin-1]->Rebin(10);
	  hInvMassUL[i][bin-1]->SetMaximum(1.5*hInvMassUL[i][bin-1]->GetMaximum());
	  hInvMassUL[i][bin-1]->GetXaxis()->SetRangeUser(2,4);
	  hInvMassUL[i][bin-1]->SetTitle("");
	  hInvMassUL[i][bin-1]->Draw("P");

	  hInvMassLS[i][bin-1] = (TH1F*)hInvMassVpPtLS[i]->ProjectionY(Form("Data_InvMass_LS_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassLS[i][bin-1]->SetMarkerStyle(24);
	  hInvMassLS[i][bin-1]->SetMarkerColor(2);
	  hInvMassLS[i][bin-1]->SetLineColor(2);
	  hInvMassLS[i][bin-1]->Rebin(10);
	  hInvMassLS[i][bin-1]->Draw("samesP");

	  TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hInvMassUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hInvMassLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_InvMass_ULvsLS_%s_InPtBins.pdf",run_type,name[i]));
    }
  */

  // UL vs LS
  TH1F *hUL[nHistos][nbins];
  TH1F *hLS[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("Data_UL_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->Rebin(rebin[i]);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	  hUL[i][bin-1]->SetTitle("");
	  hUL[i][bin-1]->Draw("P");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("Data_LS_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);
	  hLS[i][bin-1]->Rebin(rebin[i]);
	  hLS[i][bin-1]->Draw("samesP");

	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_ULvsLS_%s_InPtBins.pdf",run_type,name[i]));
    }

  for(int i=0; i<nHistos; i++)
    {
      hDataDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      hDataDisVsPt[i]->GetYaxis()->SetRangeUser(minimum[i], maximum[i]);
      c = draw2D(hDataDisVsPt[i],Form("Data: %s distribution for J/#psi muon;p_{T,#mu} (GeV/c);%s%s",title[i],title[i],unit[i]));
      if(i<2)
	{
	  TF1 *func11 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_upBound",fResVsPt[i]->GetName()));
	  func11->SetRange(0.5,ptbound);
	  func11->SetParameter(0,nsigma1*func11->GetParameter(0));
	  func11->SetParameter(1,nsigma1*func11->GetParameter(1));
	  func11->SetLineColor(2);
	  func11->Draw("sames");
	  TF1 *func12 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_lowBound",fResVsPt[i]->GetName()));
	  func12->SetRange(0.5,ptbound);
	  func12->SetParameter(0,-nsigma1*func12->GetParameter(0));
	  func12->SetParameter(1,-nsigma1*func12->GetParameter(1));
	  func12->SetLineColor(2);
	  func12->Draw("sames");
	  if(ptbound<10)
	    {
	      TF1 *func21 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_upBound",fResVsPt[i]->GetName()));
	      func21->SetRange(ptbound,10);
	      func21->SetParameter(0,nsigma2*func21->GetParameter(0));
	      func21->SetParameter(1,nsigma2*func21->GetParameter(1));
	      func21->SetLineColor(2);
	      func21->Draw("sames");
	      TF1 *func22 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_lowBound",fResVsPt[i]->GetName()));
	      func22->SetRange(ptbound,10);
	      func22->SetParameter(0,-nsigma2*func22->GetParameter(0));
	      func22->SetParameter(1,-nsigma2*func22->GetParameter(1));
	      func22->SetLineColor(2);
	      func22->Draw("sames");
	    }
	}
      else
	{
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_vs_pt.pdf",run_type,name[i]));
    }

  // Data vs embed
  TH1F *hData[nHistos][nbins];
  TH1F *hEmbed[nHistos][nbins];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_distribution",name[i]),Form("%s_distribution",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  c->cd(bin);
	  int start_bin = hDataDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hDataDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hData[i][bin-1] = (TH1F*)hDataDisVsPt[i]->ProjectionY(Form("Data_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hData[i][bin-1]->SetMarkerStyle(20);
	  hData[i][bin-1]->SetMarkerColor(1);
	  hData[i][bin-1]->SetLineColor(1);
	  hData[i][bin-1]->Rebin(rebin[i]);
	  hData[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	  hData[i][bin-1]->Scale(1./hData[i][bin-1]->Integral());
	  hData[i][bin-1]->SetTitle("");
	  hData[i][bin-1]->Draw("P");

	  hEmbed[i][bin-1] = (TH1F*)hEmbedDisVsPt[i]->ProjectionY(Form("Emebd_%s_bin%d",name[i],bin),start_bin,end_bin);
	  hEmbed[i][bin-1]->SetLineColor(2);
	  hEmbed[i][bin-1]->Rebin(rebin[i]);
	  hEmbed[i][bin-1]->Scale(1./hEmbed[i][bin-1]->Integral());
	  hEmbed[i][bin-1]->Draw("sames HIST");

	  double up = hData[i][bin-1]->GetMaximum();
	  double low = hData[i][bin-1]->GetMinimum();
	  double pt = (xbins[bin-1]+xbins[bin])/2;
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else
	    { min = -1; max = 3; }

	  TLine *line = GetLine(min,low,min,up,4);
	  line->Draw();
	  TLine *line = GetLine(max,low,max,up,4);
	  line->Draw();

          TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hData[i][0],"Data","P");
      leg->AddEntry(hEmbed[i][0],"Embedding","L");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataVsEmbed_%s_InPtBins.pdf",run_type,name[i]));
    }

  // Fit data to extract efficiency
  TH1F *hFitDataEff[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hFitDataEff[i] = new TH1F(Form("Data_JpsiMuon_%s_FitEff",name[i]),Form("%s efficiency from fitting data;p_{T} (GeV/c)",name[i]),nbins,xbins);
      TCanvas *c = new TCanvas(Form("FitData_%s",name[i]),Form("FitData_%s",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  TH1F *hFit = (TH1F*)hData[i][bin-1]->Clone(Form("Fit_%s",hData[i][bin-1]->GetName()));
	  TF1 *func = new TF1(Form("func_%d_%d",i,bin),"gaus",minimum[i], maximum[i]);
	  func->SetParameter(2,5);
	  hFit->Fit(func,"IR0Q");
	  c->cd(bin);
	  hFit->Draw();
	  func->SetLineColor(4);
	  func->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  double min = -999, max = -999;
	  double pt = (xbins[bin-1]+xbins[bin])/2;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else
	    { min = -1; max = 3; }
	  double all = func->Integral(minimum[i], maximum[i]);
	  double all_err = func->IntegralError(minimum[i], maximum[i]);
	  double acc = func->Integral(min,max);
	  double acc_err = func->IntegralError(min,max);
	  hFitDataEff[i]->SetBinContent(bin,acc/all);
	  hFitDataEff[i]->SetBinError(bin,acc_err/all);
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Fit_Data_muon_%s.pdf",run_type,name[i]));
    }
  

  TH1F *hDisAll[nHistos], *hDisAcc[nHistos], *hEff[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hDisAll[i] = new TH1F(Form("Data_JpsiMuon_%s_all",name[i]),Form("%s distribution of input muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      hDisAcc[i] = new TH1F(Form("Data_JpsiMuon_%s_Acc",name[i]),Form("%s distribution of accepted muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hDataDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hDataDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *htmp = (TH1F*)hDataDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  double all_err;
	  double all = htmp->IntegralAndError(0,-1,all_err);
	  hDisAll[i]->SetBinContent(bin,all);
	  hDisAll[i]->SetBinError(bin,all_err);

	  double pt = hDisAll[i]->GetBinCenter(bin);
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else
	    { min = -1; max = 3; }
	  double acc_err;
	  double acc = htmp->IntegralAndError(htmp->FindBin(min),htmp->FindBin(max),acc_err);
	  if(i==-1)
	    {
	      cout << all << "  " << all_err << endl;
	      cout << htmp->FindBin(min) << "  " << htmp->FindBin(max) << endl;
	      cout << acc << "  " << acc_err << endl;
	    }
	  hDisAcc[i]->SetBinContent(bin,acc);
	  hDisAcc[i]->SetBinError(bin,acc_err);
	}
      hEff[i] = (TH1F*)hDisAcc[i]->Clone(Form("Data_JpsiMuon_%s_Eff",name[i]));
      hEff[i]->Divide(hDisAll[i]);
      hEff[i]->SetMarkerStyle(21);
      hEff[i]->GetYaxis()->SetRangeUser(0.3,1.3);
      c = draw1D(hEff[i],Form("Single muon efficiency of %s cut;p_{T,#mu} (GeV/c);Efficiency",title[i]));
      gPad->SetGridy();
      hEmbedEff[i]->SetMarkerColor(2);
      hEmbedEff[i]->Draw("sames");
      hFitDataEff[i]->SetMarkerStyle(33);
      hFitDataEff[i]->SetMarkerSize(2);
      hFitDataEff[i]->SetMarkerColor(4);
      hFitDataEff[i]->SetLineColor(4);
      hFitDataEff[i]->Draw("sames");
      TH1F *hSys = 0x0;
      if(i==2)
	{
	  hSys = (TH1F*)hEmbedEff[i]->Clone(Form("%s_sys",hEmbedEff[i]->GetName()));
	  for(int ibin=1; ibin<=hSys->GetNbinsX(); ibin++)
	    {
	      hSys->SetBinContent(ibin,hSys->GetBinContent(ibin)-0.05);
	    }
	  hSys->SetLineColor(2);
	  hSys->SetLineWidth(2);
	  hSys->SetLineStyle(2);
	  hSys->Draw("HISTsames");
	}
      TLegend *leg = new TLegend(0.2,0.15,0.4,0.35);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hEff[i],"Data counting","P");
      leg->AddEntry(hFitDataEff[i],"Data fitting","P");
      leg->AddEntry(hEmbedEff[i],"Embedding","P");
      if(hSys) leg->AddEntry(hSys,"Uncertainty","L");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Data_muon_%s_eff.pdf",run_type,name[i]));
    }

  if(saveHisto)
    {
      fembed->cd();
      for(int i=0; i<nHistos; i++)
	{
	  hDisAll[i]->Write("",TObject::kOverwrite);
	  hDisAcc[i]->Write("",TObject::kOverwrite);
	  hEff[i]->Write("",TObject::kOverwrite);
	  hFitDataEff[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void anaEmbed(const bool savePlot = 0, const bool saveHisto = 0)
{
  TFile *fEmbed = 0x0;
  if(year==2013) fEmbed = TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
  if(year==2014) fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
  const char *name[4] = {"Dz","Dy","Dtof","NsigmaPi"};
  const char *title[4] = {"#Deltaz","#Deltay","#Deltatof","n#sigma_{#pi}"};
  const char *unit[4] = {" (cm)"," (cm)", " (ns)",""};
  TH2F *hMuonDisVsPt[4];
  hMuonDisVsPt[0] = (TH2F*)fEmbed->Get("hTrkDzVsPt_MCreco_di_mu");
  hMuonDisVsPt[1] = (TH2F*)fEmbed->Get("hTrkDyVsPt_MCreco_di_mu");
  hMuonDisVsPt[2] = (TH2F*)fEmbed->Get("hTrkDtofVsPt_MCreco_di_mu");
  hMuonDisVsPt[3] = (TH2F*)fEmbed->Get("hTrkNSigmaPi_MCreco_di_mu");
  const int nbins = 24;
  const double xbins[nbins+1] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0};

  // fit MC distribution
  TH1F *hFitSigma[2];
  TF1 *funcCut[2];
  for(int i=0; i<2; i++)
    {
      hFitSigma[i] = new TH1F(Form("Embed_JpsiMuon_Fit%s_Sigma",name[i]),Form("Embed: width of %s distribution from fitting;p_{T} (GeV/c);#sigma (cm)",name[i]),nbins,xbins);
      TCanvas *cfit = new TCanvas(Form("Fit_%s",name[i]),Form("Fit_%s",name[i]),1000,600);
      cfit->Divide(6,4);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *htmp = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  cfit->cd(bin);
	  htmp->SetTitle("");
	  htmp->GetXaxis()->SetRangeUser(-50,50);
	  htmp->Draw();
	  if(htmp->GetEntries()>0)
	    {
	      TF1 *functmp = new TF1(Form("func_%s_%d",name[i],bin),"gaus",-30,30);
	      htmp->Fit(functmp,"IR0Q");
	      functmp->SetLineColor(2);
	      functmp->Draw("same");
	      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xbins[bin-1],xbins[bin]),0.07);
	      t1->Draw();
	      hFitSigma[i]->SetBinContent(bin,functmp->GetParameter(2));
	      hFitSigma[i]->SetBinError(bin,functmp->GetParError(2));
	    }
	}
      hFitSigma[i]->SetMarkerStyle(21);
      hFitSigma[i]->SetMarkerColor(2);
      hFitSigma[i]->SetLineColor(2);
      hFitSigma[i]->GetYaxis()->SetRangeUser(0,30);
      c = draw1D(hFitSigma[i],"");
      funcCut[i]= new TF1(Form("Embed_JpsiMuon_Fit%s_sigma",name[i]),"[0]+[1]*exp([2]/x)",1,10);
      funcCut[i]->SetParameters(-10,10,1);
      hFitSigma[i]->Fit(funcCut[i],"IR0");
      funcCut[i]->SetLineColor(4);
      funcCut[i]->Draw("sames");
      fResVsPt[i]->Draw("sames");
    }


  // calculate pid efficiency in MC
  TH1F *hDisAll[4], *hDisAcc[4], *hEff[4];
  for(int i=0; i<4; i++)
    {
      hDisAll[i] = new TH1F(Form("Embed_JpsiMuon_%s_all",name[i]),Form("%s distribution of input muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      hDisAcc[i] = new TH1F(Form("Embed_JpsiMuon_%s_Acc",name[i]),Form("%s distribution of accepted muon tracks;p_{T} (GeV/c)",name[i]),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *htmp = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  double all = htmp->GetEntries();
	  hDisAll[i]->SetBinContent(bin,all);
	  hDisAll[i]->SetBinError(bin,sqrt(all));
	  double pt = hDisAll[i]->GetBinCenter(bin);
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else if (i==2)
	    {
	      min = -100;
	      max = fResVsPt[i]->Eval(pt);
	    }
	  else
	    { min = -1; max = 3; }
	  double acc = htmp->Integral(htmp->FindBin(min),htmp->FindBin(max));
	  hDisAcc[i]->SetBinContent(bin,acc);
	  hDisAcc[i]->SetBinError(bin,sqrt(acc));
	}
      hEff[i] = (TH1F*)hDisAcc[i]->Clone(Form("Embed_JpsiMuon_%s_Eff",name[i]));
      hEff[i]->Divide(hDisAll[i]);

      hMuonDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      c = draw2D(hMuonDisVsPt[i],Form("Embed: %s distribution for muon tracks;p_{T} (GeV/c);%s%s",title[i],title[i],unit[i]));
      if(i<2)
	{
	  TF1 *func11 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_upBound",fResVsPt[i]->GetName()));
	  func11->SetRange(0.5,ptbound);
	  func11->SetParameter(0,nsigma1*func11->GetParameter(0));
	  func11->SetParameter(1,nsigma1*func11->GetParameter(1));
	  func11->SetLineColor(2);
	  func11->Draw("sames");
	  TF1 *func12 = (TF1*)fResVsPt[i]->Clone(Form("%s_2sigma_lowBound",fResVsPt[i]->GetName()));
	  func12->SetRange(0.5,ptbound);
	  func12->SetParameter(0,-nsigma1*func12->GetParameter(0));
	  func12->SetParameter(1,-nsigma1*func12->GetParameter(1));
	  func12->SetLineColor(2);
	  func12->Draw("sames");
	  if(ptbound<10)
	    {
	      TF1 *func21 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_upBound",fResVsPt[i]->GetName()));
	      func21->SetRange(ptbound,10);
	      func21->SetParameter(0,nsigma2*func21->GetParameter(0));
	      func21->SetParameter(1,nsigma2*func21->GetParameter(1));
	      func21->SetLineColor(2);
	      func21->Draw("sames");
	      TF1 *func22 = (TF1*)fResVsPt[i]->Clone(Form("%s_2.5sigma_lowBound",fResVsPt[i]->GetName()));
	      func22->SetRange(ptbound,10);
	      func22->SetParameter(0,-nsigma2*func22->GetParameter(0));
	      func22->SetParameter(1,-nsigma2*func22->GetParameter(1));
	      func22->SetLineColor(2);
	      func22->Draw("sames");
	    }
	}
      else if(i==2)
	{
	  TF1 *func = (TF1*)fResVsPt[i]->Clone(Form("%s_upBound",fResVsPt[i]->GetName()));
	  func->SetRange(0.5,10);
	  func->Draw("sames");
	}
      else
	{
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_vs_pt.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_vs_pt.png",run_type,name[i]));
	}

      hEff[i]->SetMarkerStyle(21);
      hEff[i]->GetYaxis()->SetRangeUser(0.3,1.2);
      c = draw1D(hEff[i],Form("Embed: efficiency of %s cut;p_{T} (GeV/c);Efficiency",title[i]));
      gPad->SetGridy();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_eff.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Embed_muon_%s_eff.png",run_type,name[i]));
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
      for(int i=0; i<4; i++)
	{
	  hMuonDisVsPt[i]->Write(Form("Embed_JpsiMuon_%sVsPt",name[i]),TObject::kOverwrite);
	  hDisAll[i]->Write("",TObject::kOverwrite);
	  hDisAcc[i]->Write("",TObject::kOverwrite);
	  hEff[i]->Write("",TObject::kOverwrite);
	  if(i<2)  funcCut[i]->Write("",TObject::kOverwrite);
	}
    }
}


//================================================
void makeData()
{
  const int nHistos = 4;
  TH2F *hMuonDisVsPt[2][nHistos][2];
  const char *type[2] = {"US","LS"};
  const char *name[nHistos] = {"Dz","Dy","Dtof","NsigmaPi"};
  const char *charge[2] = {"neg","pos"};

  TFile *fin = 0x0;
  if(year==2013) fin = TFile::Open("output/Pico.Run13.pp500.jpsi.PIDcuts.root","read");

  THnSparseF *hnJpsiMuon[nHistos];
  for(int j=0; j<nHistos; j++)
    {
      hnJpsiMuon[j] = (THnSparseF*)fin->Get(Form("mhMuon%s_di_mu",name[j]));
      hnJpsiMuon[j]->GetAxis(0)->SetRangeUser(3.0,3.2);
      for(int i=0; i<2; i++)
	{
	  hnJpsiMuon[j]->GetAxis(1)->SetRange(i+1,i+1);
	  for(int k=0; k<2; k++)
	    {
	      for(int m=0; m<2; m++)
		{
		  hnJpsiMuon[j]->GetAxis(6+m)->SetRange(1+k*2, 1+k*2);
		  TH2F* h2 = (TH2F*)hnJpsiMuon[j]->Projection(4+m,2+m);
		  h2->SetName(Form("Data_JpsiMuon_%sVsPt_%s_%s_muon%d",name[j],charge[k],type[i],m+1));
		  h2->Sumw2();
		  //cout << h2->GetName() << ": " << h2->GetNbinsX() << ", " << h2->GetNbinsY() << endl;
		  if(m==0) hMuonDisVsPt[i][j][k] = (TH2F*)h2->Clone(Form("Data_JpsiMuon_%sVsPt_%s_%s",name[j],charge[k],type[i]));
		  else     hMuonDisVsPt[i][j][k]->Add(h2);
		  hnJpsiMuon[j]->GetAxis(6+m)->SetRange(0,-1);
		}
	    }
	  hnJpsiMuon[j]->GetAxis(1)->SetRange(0,-1);
	}   
    }

  TH2F *hDis[nHistos][3];
  for(int j=0; j<nHistos; j++)
    {
      for(int k=0; k<2; k++)
	{
	  hDis[j][k] = (TH2F*)hMuonDisVsPt[0][j][k]->Clone(Form("Data_JpsiMuon_%sVsPt_%s",name[j],charge[k]));
	  hDis[j][k]->Add(hMuonDisVsPt[1][j][k],-1);
	}
      hDis[j][2] = (TH2F*)hDis[j][0]->Clone(Form("Data_JpsiMuon_%sVsPt",name[j]));
      hDis[j][2]->Add((TH2F*)hDis[j][1]);
    }

  TFile *fout = 0x0;
  if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.PIDcuts.root","update");
  for(int j=0; j<nHistos; j++)
    {
      for(int k=2; k>-1; k--)
	{
	  hDis[j][k]->SetTitle("");
	  hDis[j][k]->Write("",TObject::kOverwrite);
	  //cout << hDis[j][k]->GetName() << ": " << hDis[j][k]->GetNbinsX() << ", " << hDis[j][k]->GetNbinsY() << endl;
	  if(k<2)
	    {
	      for(int i=0; i<2; i++)
		{
		  hMuonDisVsPt[i][j][k]->SetTitle("");
		  hMuonDisVsPt[i][j][k]->Write("",TObject::kOverwrite);
		}
	    }	    
	}
    }
}


//================================================
void nSigmaPi(const bool savePlot = 0, const bool saveHisto = 0)
{
  const int nbins = 6;
  const double xbins[nbins+1] = {1.2,1.5,2.0,2.5,3.0,5.0,10.0};
  TList *list = new TList;

  const char *trgName[3] = {"low","mid","high"};
  TFile *fMB = TFile::Open("output/Run14.AuAu200.Pion.root","read");
  TH2F *hDataPion[3];
  TH1F *hDataPionMean[3];
  TH1F *hDataPionSigma[3];
  TH1F *hDataPionEff[3];
  for(int i=0; i<3; i++)
    {
      hDataPion[i] = (TH2F*)fMB->Get(Form("hPtnSigmaPi_%s",trgName[i]));
      ScaleHistoTitle(hDataPion[i],0.04,1.1,0.035,0.04,0.8,0.035,62);
      c = draw2D(hDataPion[i],Form("Run14_AuAu200: n#sigma_{#pi} vs p_{T} for #pi from weak decays (prod_%s)",trgName[i]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataPion_nSigmaPiVsPt_%s.pdf",run_type,trgName[i]));

      hDataPionMean[i] = new TH1F(Form("hDataPionMean_%s",trgName[i]),Form("Data: mean of n#sigma_{#pi} distribution for #pi sample;p_{T} (GeV/c);<n#sigma_{#pi}>"),nbins,xbins);
      hDataPionSigma[i] = new TH1F(Form("hDataPionSigma_%s",trgName[i]),Form("Data: width of n#sigma_{#pi} distribution for #pi sample;p_{T} (GeV/c);#sigma(n#sigma_{#pi})"),nbins,xbins);

      hDataPionEff[i] = new TH1F(Form("hDataPionEff_%s",trgName[i]),Form("n#sigma_{#pi} efficiency from fitting data;p_{T} (GeV/c)"),nbins,xbins);
      TCanvas *c = new TCanvas(Form("FitDataPion_%s",trgName[i]),Form("FitDataPion_%s",trgName[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nbins; bin++)
	{
	  double pt1 = xbins[bin-1]; 
	  double pt2 = xbins[bin];
	  int start_bin = hDataPion[i]->GetXaxis()->FindBin(pt1+1e-4);
	  int end_bin = hDataPion[i]->GetXaxis()->FindBin(pt2-1e-4);
	  TH1F *hFit = (TH1F*)hDataPion[i]->ProjectionY(Form("%s_%d",hDataPion[i]->GetName(),bin),start_bin,end_bin);
	  hFit->Rebin(4);
	  TF1 *func = new TF1(Form("funcFitDataPion_%d_%d",i,bin),"gaus",-3,3);
	  func->SetParameter(2,5);
	  hFit->Fit(func,"IR0Q");
	  c->cd(bin);
	  hFit->SetMaximum(1.2*hFit->GetMaximum());
	  hFit->Draw();
	  func->SetLineColor(2);
	  func->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("#pi: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
	  t1->Draw();
	  double all = func->Integral(-3,3);
	  double all_err = func->IntegralError(-3,3);
	  double acc = func->Integral(-1.5,2.5);
	  double acc_err = func->IntegralError(-1.5,2.5);
	  hDataPionEff[i]->SetBinContent(bin,acc/all);
	  hDataPionEff[i]->SetBinError(bin,acc_err/all);
	  hDataPionMean[i]->SetBinContent(bin,func->GetParameter(1));
	  hDataPionMean[i]->SetBinError(bin,func->GetParError(1));
	  hDataPionSigma[i]->SetBinContent(bin,func->GetParameter(2));
	  hDataPionSigma[i]->SetBinError(bin,func->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.72,0.85,0.05);
      t1->AddText("Run14_AuAu200");
      t1->AddText(Form("prod_%s",trgName[i]));
      t1->SetTextColor(4);
      t1->SetTextFont(62);
      t1->SetTextAlign(11);
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataPion_FitnSigmaPi_%s.pdf",run_type,trgName[i]));
    }
  TString legName[3];
  for(int i=0; i<3; i++)
    {
      list->Add(hDataPionMean[i]);
      legName[i] = Form("prod_%s",trgName[i]);
    }
  c = drawHistos(list,"DataPion_nSigmaPi_Mean","",false,-0.5,0.1,true,-0.5,0.1,kFALSE,kTRUE,legName,kTRUE,run_type,0.5,0.7,0.15,0.4,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataPion_nSigmaPiMean.pdf",run_type));

  for(int i=0; i<3; i++)
    {
      list->Add(hDataPionSigma[i]);
    }
  c = drawHistos(list,"DataPion_nSigmaPi_Sigma","",false,-0.5,0.1,true,0.5,1.5,kFALSE,kTRUE,legName,kTRUE,run_type,0.5,0.7,0.15,0.4,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/DataPion_nSigmaPiSigma.pdf",run_type));

  // compare data with simulation
  const int nHistos = 4;
  const char *hName[nHistos] = {"EmbedMuon","EmbedPion","DataMuon","DataPion"};
  TH2F *hNsigmaPiVsPt[nHistos];

  TFile *fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
  hNsigmaPiVsPt[2] = (TH2F*)fdata->Get(Form("mhJpsiMuonNsigmaPi_UL_di_mu"));
  TH2F *hDataLS    = (TH2F*)fdata->Get(Form("mhJpsiMuonNsigmaPi_LS_di_mu"));
  hNsigmaPiVsPt[2]->Sumw2();
  hNsigmaPiVsPt[2]->Add(hDataLS, -1);

  TFile *fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");
  hNsigmaPiVsPt[0] = (TH2F*)fEmbed->Get("hTrkNSigmaPi_MCreco_di_mu");

  hNsigmaPiVsPt[3] = (TH2F*)hDataPion[0]->Clone("DataPion_hNsigmaPiVsPt");
  hNsigmaPiVsPt[3]->Add(hDataPion[1]);

  TFile *fEmbed2 = TFile::Open("output/Run14.AuAu200.Pion.Embed.root","read");
  hNsigmaPiVsPt[1] = (TH2F*)fEmbed2->Get("mhNsigmaPiVsPt");

  TH1F *hNsigmaPiMean[nHistos];
  TH1F *hNsigmaPiSigma[nHistos];
  TH1F *hNsigmaPiEff[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      ScaleHistoTitle(hNsigmaPiVsPt[i],0.04,1.1,0.035,0.04,0.8,0.035,62);
      c = draw2D(hNsigmaPiVsPt[i],Form("%s: n#sigma_{#pi} vs p_{T} distribution",hName[i]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/%s_nSigmaPiVsPt.pdf",run_type,hName[i]));

      if(i<2)
	{
	  const int nbins1 = 11;
	  const double xbins1[nbins1+1] = {1.2,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	}
      else
	{
	  const int nbins1 = 6;
	  const double xbins1[nbins1+1] = {1.2,1.5,2.0,2.5,3.0,5.0,10.0};
	}

      hNsigmaPiMean[i] = new TH1F(Form("hNsigmaPiMean_%s",hName[i]),Form("Mean of n#sigma_{#pi} distribution;p_{T} (GeV/c);<n#sigma_{#pi}>"),nbins1,xbins1);
      hNsigmaPiSigma[i] = new TH1F(Form("hNsigmaPiSigma_%s",hName[i]),Form("Width of n#sigma_{#pi} distribution;p_{T} (GeV/c);#sigma(n#sigma_{#pi})"),nbins1,xbins1);
      hNsigmaPiEff[i] = new TH1F(Form("hNsigmaPiEff_%s",hName[i]),Form("n#sigma_{#pi} efficiency from fitting;p_{T} (GeV/c)"),nbins1,xbins1);

      TCanvas *c = new TCanvas(Form("Fit_%s",hName[i]),Form("Fit_%s",hName[i]),1100,700);
      if(i<2) c->Divide(4,3);
      else    c->Divide(3,2);

      for(int bin=1; bin<=nbins1; bin++)
	{
	  double pt1 = xbins1[bin-1]; 
	  double pt2 = xbins1[bin];
	  int start_bin = hNsigmaPiVsPt[i]->GetXaxis()->FindBin(pt1+1e-4);
	  int end_bin = hNsigmaPiVsPt[i]->GetXaxis()->FindBin(pt2-1e-4);
	  TH1F *hFit = (TH1F*)hNsigmaPiVsPt[i]->ProjectionY(Form("%s_%d",hNsigmaPiVsPt[i]->GetName(),bin),start_bin,end_bin);
	  if(i>1) hFit->Rebin(2);
	  TF1 *func = new TF1(Form("funcFit_%s_%d",hName[i],bin),"gaus",-3,3);
	  func->SetParameter(2,5);
	  hFit->Fit(func,"IR0Q");
	  c->cd(bin);
	  hFit->SetMaximum(1.2*hFit->GetMaximum());
	  hFit->SetMarkerStyle(21);
	  hFit->GetXaxis()->SetRangeUser(-5,5);
	  hFit->Draw();
	  func->SetLineColor(2);
	  func->Draw("sames");
	  TPaveText *t1 = 0x0;
	  if(i%2==0) t1 = GetTitleText(Form("#mu: %1.1f < p_{T} < %1.1f",xbins1[bin-1],xbins1[bin]),0.06);
	  else       t1 = GetTitleText(Form("#pi: %1.1f < p_{T} < %1.1f",xbins1[bin-1],xbins1[bin]),0.06);
	  t1->Draw();
	  double all = func->Integral(-3,3);
	  double all_err = func->IntegralError(-3,3);
	  double acc, acc_err;
	  if(i%2==0)
	    {
	      acc = func->Integral(-1,3);
	      acc_err = func->IntegralError(-1,3);
	    }
	  else
	    {
	      acc = func->Integral(-1.5,2.5);
	      acc_err = func->IntegralError(-1.5,2.5);
	    }
	  hNsigmaPiEff[i]->SetBinContent(bin,acc/all);
	  hNsigmaPiEff[i]->SetBinError(bin,acc_err/all);
	  hNsigmaPiMean[i]->SetBinContent(bin,func->GetParameter(1));
	  hNsigmaPiMean[i]->SetBinError(bin,func->GetParError(1));
	  hNsigmaPiSigma[i]->SetBinContent(bin,func->GetParameter(2));
	  hNsigmaPiSigma[i]->SetBinError(bin,func->GetParError(2));
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.15,0.55,0.72,0.85,0.05);
      t1->AddText("Run14_AuAu200");
      t1->AddText(Form("%s",hName[i]));
      t1->SetTextColor(4);
      t1->SetTextFont(62);
      t1->SetTextAlign(11);
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/%s_FitnSigmaPi.pdf",run_type,hName[i]));
    }
  TString legName2[nHistos] = {"Embedding: muon","Embedding: pion","Data: muon","Data: pion"};
  for(int i=0; i<nHistos; i++)
    {
      list->Add(hNsigmaPiMean[i]);
    }
  c = drawHistos(list,"nSigmaPi_Mean","",false,-0.5,0.1,true,-0.5,0.5,kFALSE,kTRUE,legName2,kTRUE,run_type,0.5,0.7,0.15,0.4,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Compare_nSigmaPiMean.pdf",run_type));

  for(int i=0; i<nHistos; i++)
    {
      list->Add(hNsigmaPiSigma[i]);
    }
  c = drawHistos(list,"nSigmaPi_Sigma","",false,-0.5,0.1,true,0.2,1.2,kFALSE,kTRUE,legName2,kTRUE,run_type,0.5,0.7,0.15,0.4,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Compare_nSigmaPiSigma.pdf",run_type));

  for(int i=0; i<nHistos; i++)
    {
      list->Add(hNsigmaPiEff[i]);
    }
  c = drawHistos(list,"nSigmaPi_Eff","",false,-0.5,0.1,true,0.7,1,kFALSE,kTRUE,legName2,kTRUE,run_type,0.5,0.7,0.15,0.4,kTRUE);
  list->Clear();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_PIDcuts/Compare_nSigmaPiEff.pdf",run_type));
}
