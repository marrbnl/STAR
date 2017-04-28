const int year = YEAR;

TFile *f;

//================================================
void ana_MtdTrigEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //trigElecEff();
  //anaTrigEff();

  //make_ppAu();
  make_AuAu();
}
//================================================
void make_AuAu(const int savePlot = 0, const int saveHisto = 0)
{
  const int year = 2014;
  const char* config = "";
  const char *data_name = "Run14_AuAu200";
  const char *data_title  = "Run14 Au+Au @ 200 GeV";
  const char* type_name[3] = {"Muon","UL", "LS"};
  const int nbins = 7;
  const double lowbins[nbins] = {1.3, 1.3, 1.5, 2.0, 2.5, 3.0, 5.0};
  const double upbins[nbins]  = {10,  1.5, 2.0, 2.5, 3.0, 5.0, 10.0};
  const int nPtBins = nbins -1;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 
  if(year==2014)
    {
      const int nBinsTacDiff = 22;
      const double xBinsTacDiff[nBinsTacDiff+1] = {760,765,770,775,780,782,784,786,788,789,793,797,801,805,809,813,817,821,825,829,833,837,841};
    }
  else if(year==2016)
    {
      const int nBinsTacDiff = 19;
      const double xBinsTacDiff[nBinsTacDiff+1] = {920,930,935,940,945,950,951,955,960,965,970,975,980,985,990,995,1000,1005,1010,1015};
    }

  TFile *fdata = TFile::Open(Form("output/Run14.AuAu200.JpsiMuon%s.root",config),"read");

  //==============================================
  // compare invariant mass
  //==============================================
  const double min_mass[3] = {3.0, 2.6, 3.4};
  const double max_mass[3] = {3.2, 2.8, 3.6};
  const int nbinsJpsi = 5;
  const double xbinsJpsi[nbinsJpsi+1] = {0,1.0,2.0,3.0,5.0,10.0};
  TH2F *hInvMassVsPt[2];
  TH1F *hInvMass[2][nbinsJpsi];
  for(int i=0; i<2; i++)
    {
      hInvMassVsPt[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_MtdVpdTacDiff_%s_di_mu",type_name[i+1]));
      hInvMassVsPt[i]->Sumw2();
    }
      
  TCanvas *c = new TCanvas("InvMass", "InvMass",1100,700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetHeader(data_name);
  for(int bin=1; bin<=nbinsJpsi; bin++)
    {
      for(int i=0; i<2; i++)
	{
	  int start_bin = hInvMassVsPt[i]->GetXaxis()->FindBin(xbinsJpsi[bin-1]+1e-4);
	  int end_bin   = hInvMassVsPt[i]->GetXaxis()->FindBin(xbinsJpsi[bin]-1e-4);
	  hInvMass[i][bin-1] = (TH1F*)hInvMassVsPt[i]->ProjectionY(Form("%s_TacDiff_InvMass%s_bin%d",data_name,type_name[i+1],bin),start_bin,end_bin);
	  hInvMass[i][bin-1]->Rebin(5);
	  hInvMass[i][bin-1]->SetMarkerStyle(20+i*4);
	  hInvMass[i][bin-1]->SetMarkerColor(i+1);
	  hInvMass[i][bin-1]->SetLineColor(i+1);
	  hInvMass[i][bin-1]->SetMaximum(1.5*hInvMass[i][bin-1]->GetMaximum());
	  hInvMass[i][bin-1]->SetMinimum(0);
	  hInvMass[i][bin-1]->GetXaxis()->SetRangeUser(2.2,4);
	  hInvMass[i][bin-1]->SetTitle("");
	  if(bin==1) leg->AddEntry(hInvMass[i][bin-1],type_name[i+1],"PL");
	}

      c->cd(bin);
      hInvMass[0][bin-1]->Draw("P");
      double binwidth = hInvMass[0][bin-1]->GetBinWidth(1);
      for(int itmp=0; itmp<3; itmp++)
	{
	  TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMass[0][bin-1]->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	    {
	      int jbin = hInvMass[0][bin-1]->FindFixBin(htmp->GetBinCenter(ibin));
	      htmp->SetBinContent(ibin,hInvMass[0][bin-1]->GetBinContent(jbin));
	    }
	  if(itmp==0) htmp->SetFillColor(7);
	  else htmp->SetFillColor(5);
	  htmp->Draw("sames");
	  if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	  if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	}
      hInvMass[1][bin-1]->Draw("samesP");
      TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbinsJpsi[bin-1],xbinsJpsi[bin]),0.055);
      t1->Draw();
      int low_bin  = hInvMass[0][bin-1]->FindFixBin(min_mass[0]+1e-4);
      int high_bin = hInvMass[0][bin-1]->FindFixBin(max_mass[0]-1e-4);
      double all = hInvMass[0][bin-1]->Integral(low_bin, high_bin);
      double bkg = hInvMass[1][bin-1]->Integral(low_bin, high_bin);
      TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
      t1->SetTextAlign(11);
      t1->AddText(Form("UL-LS = %1.0f",all-bkg));
      t1->AddText(Form("S/B = %4.2f:1",(all-bkg)/bkg));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s.TacDiff_InvMassLSvsUL.pdf",data_name,config));

  //==============================================
  // scale factor
  //==============================================
  TH2F *hJpsiMassVsMuonPt[2];
  TH1F *hMuonPtSide[2];
  for(int i=0; i<2; i++)
    {
      hJpsiMassVsMuonPt[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_MtdVpdTacDiff_%s_di_mu",type_name[i+1]));
      hJpsiMassVsMuonPt[i]->SetName(Form("%s_%s",data_name,hJpsiMassVsMuonPt[i]->GetName()));
      hJpsiMassVsMuonPt[i]->Sumw2();
      hMuonPtSide[i] = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("%s_TacDiff_SideBand_%s",data_name,type_name[i+1]));
      hMuonPtSide[i]->Reset();
      for(int k=0; k<2; k++)
	{
	  int low_bin = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(min_mass[k+1]);
	  int up_bin  = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(max_mass[k+1]);
	  TH1F *htmp = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("%s_TacDiff_SideBand_%s_%d",data_name,type_name[i+1],k),low_bin,up_bin);
	  hMuonPtSide[i]->Add(htmp);
	}
    }
  TGraphErrors *gScaleFactor = new TGraphErrors(nbins);
  gScaleFactor->SetName(Form("%s_ScaleFactor",data_name));
  gScaleFactor->SetTitle(Form("%s: scale factor for like-sign background (TacDiff);p_{T} (GeV/c)",data_name));
  for(int bin=0; bin<nbins; bin++)
    {
      int start_bin = hMuonPtSide[0]->GetXaxis()->FindBin(lowbins[bin]+1e-4);
      int end_bin   = hMuonPtSide[0]->GetXaxis()->FindBin(upbins[bin]-1e-4);
      double ul = hMuonPtSide[0]->Integral(start_bin,end_bin);
      double ls = hMuonPtSide[1]->Integral(start_bin,end_bin);
      double scale = ul/ls;
      double error = scale * sqrt(1/ul+1/ls);
      double pt = 0.5*(lowbins[bin]+upbins[bin]);
      if(bin==0) pt = -0.5;
      gScaleFactor->SetPoint(bin, pt, scale);
      gScaleFactor->SetPointError(bin, 0, error);
    }
  gScaleFactor->SetMarkerStyle(21);
  c = drawGraph(gScaleFactor);
  TPaveText *t1 = GetPaveText(0.15,0.2,0.2,0.23,0.03);
  t1->AddText(Form("%1.1f < p_{T} < %1.1f",lowbins[0],upbins[0]));
  t1->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s.TacDiff_ScaleFactor.pdf",data_name,config));
  double *scales = gScaleFactor->GetY();

  //==============================================
  // MtdVpdTacDiff study
  //==============================================
  THnSparseF *hn = (THnSparseF*)fdata->Get("mhJpsiMuonTrigEff_di_mu");
  TH2F *hTacDiffVsTrigUnit[3][nbins];
  TH2F *hTacDiffVsTrigUnitBkg[2][nbins];
  for(int bin=0; bin<nbins; bin++)
    {
      hn->GetAxis(2)->SetRangeUser(lowbins[bin]+1e-4, upbins[bin]-1e-4);
      
      
	hn->GetAxis(0)->SetRangeUser(min_mass[0]+1e-4,max_mass[0]-1e-4);
      // signal
      for(int i=0; i<2; i++)
	{
	  hn->GetAxis(5)->SetRange(i+1,i+1);
	  hTacDiffVsTrigUnit[i+1][bin] = (TH2F*)hn->Projection(1,3);
	  hTacDiffVsTrigUnit[i+1][bin]->SetName(Form("%s_hTacDiffVsVzTrigUnit_%s_PtBin%d",data_name,type_name[i+1],bin));
	  hTacDiffVsTrigUnit[i+1][bin]->Sumw2();
	  hn->GetAxis(5)->SetRange(0,-1);
	}

      // background
      
      hn->GetAxis(2)->SetRange(0,-1);
      
      

      /*
      printf("[i] pt bin = %d\n",bin);
      // unlike-sign signal
      hn->GetAxis(5)->SetRange(1,1);
      hn->GetAxis(0)->SetRangeUser(min_mass[0]+1e-4,max_mass[0]-1e-4);
      hTacDiffVsTrigUnit[1][bin] = (TH2F*)hn->Projection(1,3);
      hTacDiffVsTrigUnit[1][bin]->SetName(Form("%s_hTacDiffVsVzTrigUnit_UL_PtBin%d",data_name,bin));
      hTacDiffVsTrigUnit[1][bin]->Sumw2();
      printf("UL [%2.1f, %2.1f] = %1.0f\n",min_mass[0],max_mass[0],hTacDiffVsTrigUnit[1][bin]->GetEntries());

      // unlike-sign side-band background
      hTacDiffVsTrigUnit[2][bin] = (TH2F*)hTacDiffVsTrigUnit[1][bin]->Clone(Form("%s_hTacDiffVsVzTrigUnit_Bkg_PtBin%d",data_name,bin));
      hTacDiffVsTrigUnit[2][bin]->Reset();
      for(int k=0; k<2; k++)
	{
	  hn->GetAxis(0)->SetRangeUser(min_mass[k+1]+1e-4,max_mass[k+1]-1e-4);
	  TH2F *h2tmp = (TH2F*)hn->Projection(1,3);
	  h2tmp->SetName(Form("%s_hTacDiffVsVzTrigUnit_UL_SB%d_PtBin%d",data_name,k,bin));
	  h2tmp->Sumw2();
	  hTacDiffVsTrigUnit[2][bin]->Add(h2tmp);
	  hn->GetAxis(0)->SetRange(0,-1);
	  printf("UL [%2.1f, %2.1f] = %1.0f\n",min_mass[k+1],max_mass[k+1],h2tmp->GetEntries());
	}
      hn->GetAxis(0)->SetRange(0,-1);
      
      // like-sign background
      hn->GetAxis(5)->SetRange(2,2);
      hn->GetAxis(0)->SetRangeUser(min_mass[1]+1e-4,max_mass[2]-1e-4);
      TH2F *h2tmp = (TH2F*)hn->Projection(1,3);
      h2tmp->SetName(Form("%s_hTacDiffVsVzTrigUnit_LS_PtBin%d",data_name,bin));
      h2tmp->Sumw2();
      hTacDiffVsTrigUnit[2][bin]->Add(h2tmp);
      printf("LS [%2.1f, %2.1f] = %1.0f\n",min_mass[1],max_mass[2],h2tmp->GetEntries());
      hn->GetAxis(5)->SetRange(0,-1);
      hn->GetAxis(0)->SetRange(0,-1);
      hn->GetAxis(2)->SetRange(0,-1);

      double factor = (2*max_mass[2]-2*min_mass[1]-min_mass[2]+max_mass[1])/(max_mass[0]-min_mass[0]);
      hTacDiffVsTrigUnit[2][bin]->Scale(1./factor);
      */

      // pure muon distribution
      hTacDiffVsTrigUnit[0][bin] = (TH2F*)hTacDiffVsTrigUnit[1][bin]->Clone(Form("%s_hTacDiffVsVzTrigUnit_Muon_PtBin%d",data_name,bin));
      hTacDiffVsTrigUnit[0][bin]->Add(hTacDiffVsTrigUnit[2][bin], -1*scales[bin]);
      cout << scales[bin] << endl;

      for(int k=0; k<3; k++)
	hTacDiffVsTrigUnit[k][bin]->GetYaxis()->SetRangeUser(760, 860);
    }
  c = draw2D(hTacDiffVsTrigUnit[0][0],Form("%s: J/#psi #mu %1.1f < p_{T} < %1.1f GeV/c",data_name,lowbins[0],upbins[0]));
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s.TacDiffVsTrigUnit_PtBin0.pdf",data_name,config));

  //==============================================
  // Show TacDiff distribution in each trigger unit
  int colors[5] = {1, 2, 4, 6, kGreen+2};
  const int nTrigUnit = 28;
  for(int bin=0; bin<nbins; bin++)
    {
      c = new TCanvas(Form("%s_TacDiffTrigUnit_ULvsLS_PtBin%d",data_name,bin),Form("%s_TacDiffTrigUnit_ULvsLS_PtBin%d",data_name,bin),1200,800);
      c->Divide(6,5);
      TLegend *leg = new TLegend(0.1,0.3,0.5,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.09);
      leg->SetHeader(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]));
      for(int j=0; j<nTrigUnit; j++)
	{
	  c->cd(j+2);
	  SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	  for(int k=2; k>-1; k--)
	    {
	      TH1F *htmp = (TH1F*)hTacDiffVsTrigUnit[k][bin]->ProjectionY(Form("%s_hTacDiff_%s_TrigUnit%d_PtBin%d",data_name,type_name[k],j+1,bin),j+2,j+2);
	      TH1F *hRebin = (TH1F*)htmp->Rebin(nBinsTacDiff, Form("%s_Rebin",htmp->GetName()), xBinsTacDiff);
	      hRebin->SetMarkerStyle(20+k);
	      hRebin->SetMaximum(1.5*hRebin->GetMaximum());
	      hRebin->SetMarkerColor(TMath::Power(2,k));
	      ScaleHistoTitle(hRebin, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	      hRebin->SetTitle("");
	      hRebin->Draw("samesPE");
	      if(j==0) leg->AddEntry(hRebin,type_name[k],"P");
	    }
	  TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	  t1->Draw();
	}
      c->cd(1);
      leg->Draw();
    }
  //return;


  //==============================================
  // Fit TacDiff distribution in each trigger unit
  TH1F *hTacDiffTrigUnit[nbins][nTrigUnit];
  TH1F *hTacDiffTrigUnitRebin[nbins][nTrigUnit];
  TH1F *hTacDiffMeanVsTrigUnit[nbins]; 
  TH1F *hTacDiffWidthVsTrigUnit[nbins]; 
  TH1F *hTacDiffMeanVsPt[nTrigUnit];
  TF1  *funcTacDiffTrigUnit[nbins][nTrigUnit];
  double fit_min = 789;
  double fit_max = 840;

  for(int j=0; j<nTrigUnit; j++)
    {
      hTacDiffMeanVsPt[j] = new TH1F(Form("%s_hTacDiffMeanVsPt_TrigUnit%d",data_name,j+1),Form("%s: mean of MtdVpdTacDiff vs. p_{T};p_{T} (GeV/c);Mean",data_name),nPtBins, xPtBins);
      for(int bin=0; bin<nbins; bin++)
	{
	  funcTacDiffTrigUnit[bin][j] = new TF1(Form("%s_FitTacDiff_TrigUnit%d_PtBin%d",data_name,j+1,bin),"gaus",fit_min,fit_max);
	}
    }
  for(int bin=0; bin<nbins; bin++)
    {
      hTacDiffMeanVsTrigUnit[bin] = new TH1F(Form("%s_hTacDiffMeanVsTrigUnit_PtBin%d",data_name,bin),Form("%s: mean of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);TrigUnit;Mean",data_name,lowbins[bin],upbins[bin]),nTrigUnit,1,nTrigUnit+1);
      hTacDiffWidthVsTrigUnit[bin] = new TH1F(Form("%s_hTacDiffWidthVsTrigUnit_PtBin%d",data_name,bin),Form("%s: width of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);TrigUnit;#sigma",data_name,lowbins[bin],upbins[bin]),nTrigUnit,1,nTrigUnit+1);
    }

  for(int bin=0; bin<nbins; bin++)
    {
      c = new TCanvas(Form("%s_TacDiff_vs_TrigUnit_PtBin%d",data_name,bin),Form("%s_TacDiff_vs_TrigUnit_PtBin%d",data_name,bin),1200,800);
      c->Divide(6,5);
      for(int j=0; j<nTrigUnit; j++)
	{
	  hTacDiffTrigUnit[bin][j] = (TH1F*)hTacDiffVsTrigUnit[0][bin]->ProjectionY(Form("%s_hTacDiff_TrigUnit%d_PtBin%d",data_name,j+1,bin),j+2,j+2);
	  hTacDiffTrigUnit[bin][j]->SetMarkerStyle(20);
	  hTacDiffTrigUnit[bin][j]->SetMarkerSize(0.8);
	  hTacDiffTrigUnitRebin[bin][j] = (TH1F*)hTacDiffTrigUnit[bin][j]->Rebin(nBinsTacDiff, Form("%s_hTacDiff_TrigUnit%d_PtBin%d_Rebin",data_name,j+1,bin), xBinsTacDiff);
	  c->cd(j+2);
	  SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	  TH1F *hFit = (TH1F*)hTacDiffTrigUnitRebin[bin][j]->Clone(Form("Fit_%s",hTacDiffTrigUnitRebin[bin][j]->GetName()));
	  hFit->Fit(funcTacDiffTrigUnit[bin][j],"IR0QL");
	  hTacDiffMeanVsTrigUnit[bin]->SetBinContent(j+1, funcTacDiffTrigUnit[bin][j]->GetParameter(1));
	  hTacDiffMeanVsTrigUnit[bin]->SetBinError(j+1, funcTacDiffTrigUnit[bin][j]->GetParError(1));
	  hTacDiffWidthVsTrigUnit[bin]->SetBinContent(j+1, funcTacDiffTrigUnit[bin][j]->GetParameter(2));
	  hTacDiffWidthVsTrigUnit[bin]->SetBinError(j+1, funcTacDiffTrigUnit[bin][j]->GetParError(2));
	  hTacDiffMeanVsPt[j]->SetBinContent(bin, funcTacDiffTrigUnit[bin][j]->GetParameter(1));
	  hTacDiffMeanVsPt[j]->SetBinError(bin, funcTacDiffTrigUnit[bin][j]->GetParError(1));
	  ScaleHistoTitle(hFit, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	  hFit->SetTitle("");
	  hFit->Draw("samesPE");
	  funcTacDiffTrigUnit[bin][j]->SetLineColor(4);
	  funcTacDiffTrigUnit[bin][j]->SetLineStyle(2);
	  funcTacDiffTrigUnit[bin][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	  t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.1,0.3,0.5,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.09);
      leg->SetHeader(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]));
      leg->AddEntry(hTacDiffTrigUnit[bin][0],Form("%s",data_title),"P");
      leg->AddEntry(funcTacDiffTrigUnit[bin][0],"Fit","L");
      leg->Draw();
    }
}

//================================================
void make_ppAu(const int savePlot = 1, const int saveHisto = 1)
{
  const char* config = "";
  const int nData = 2;
  const char *data_name[3] = {"Run15_pp200", "Run15_pAu200","Run15_ppAu200"};
  const char *data_leg[3]  = {"p+p @ 200 GeV", "p+Au @ 200 GeV", "p+p&p+Au @ 200 GeV"};
  const char* type_name[3] = {"Muon","UL", "LS"};
  const int nbins = 7;
  const double lowbins[nbins] = {1.3, 1.3, 1.5, 2.0, 2.5, 3.0, 5.0};
  const double upbins[nbins]  = {10,  1.5, 2.0, 2.5, 3.0, 5.0, 10.0};
  const int nPtBins = nbins -1;
  const double xPtBins[nPtBins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0}; 

  TFile *fdata[nData];
  for(int i=0; i<nData; i++)
    {
      fdata[i] = TFile::Open(Form("output/%s.JpsiMuon.%sroot",data_name[i],config),"read");
    }

  //==============================================
  // compare invariant mass
  //==============================================
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  const int nbinsJpsi = 5;
  const double xbinsJpsi[nbinsJpsi+1] = {0,1.0,2.0,3.0,5.0,10.0};
  TH2F *hInvMassVsPt[nData][2];
  TH1F *hInvMass[nData][2][nbinsJpsi];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hInvMassVsPt[i][j] = (TH2F*)fdata[i]->Get(Form("mhJpsiMassVsPt_MtdVpdTacDiff_%s_di_mu",type_name[j+1]));
	  hInvMassVsPt[i][j]->Sumw2();
	  hInvMassVsPt[i][j]->SetName(Form("%s_%s",data_name[i],hInvMassVsPt[i][j]->GetName()));
	}
      
      TCanvas *c = new TCanvas(Form("%s_InvMass",data_name[i]),Form("%s_InvMass",data_name[i]),1100,700);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(data_name[i]);
      for(int bin=1; bin<=nbinsJpsi; bin++)
	{
	  for(int j=0; j<2; j++)
	    {
	      int start_bin = hInvMassVsPt[i][j]->GetXaxis()->FindBin(xbinsJpsi[bin-1]+1e-4);
	      int end_bin   = hInvMassVsPt[i][j]->GetXaxis()->FindBin(xbinsJpsi[bin]-1e-4);
	      hInvMass[i][j][bin-1] = (TH1F*)hInvMassVsPt[i][j]->ProjectionY(Form("%s_TacDiff_InvMass%s_bin%d",data_name[i],type_name[j+1],bin),start_bin,end_bin);
	      hInvMass[i][j][bin-1]->Rebin(5);
	      hInvMass[i][j][bin-1]->SetMarkerStyle(20+j*4);
	      hInvMass[i][j][bin-1]->SetMarkerColor(j+1);
	      hInvMass[i][j][bin-1]->SetLineColor(j+1);
	      hInvMass[i][j][bin-1]->SetMaximum(1.5*hInvMass[i][j][bin-1]->GetMaximum());
	      hInvMass[i][j][bin-1]->SetMinimum(0);
	      hInvMass[i][j][bin-1]->GetXaxis()->SetRangeUser(2.2,4);
	      hInvMass[i][j][bin-1]->SetTitle("");
	      if(bin==1) leg->AddEntry(hInvMass[i][j][bin-1],type_name[j+1],"PL");
	    }

	  c->cd(bin);
	  hInvMass[i][0][bin-1]->Draw("P");
	  double binwidth = hInvMass[i][0][bin-1]->GetBinWidth(1);
	  for(int itmp=0; itmp<3; itmp++)
	    {
	      TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMass[i][0][bin-1]->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	      for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		{
		  int jbin = hInvMass[i][0][bin-1]->FindFixBin(htmp->GetBinCenter(ibin));
		  htmp->SetBinContent(ibin,hInvMass[i][0][bin-1]->GetBinContent(jbin));
		}
	      if(itmp==0) htmp->SetFillColor(7);
	      else htmp->SetFillColor(5);
	      htmp->Draw("sames");
	      if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	      if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	    }
	  hInvMass[i][1][bin-1]->Draw("samesP");
	  TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbinsJpsi[bin-1],xbinsJpsi[bin]),0.055);
	  t1->Draw();
	  int low_bin  = hInvMass[i][0][bin-1]->FindFixBin(min_mass[0]+1e-4);
	  int high_bin = hInvMass[i][0][bin-1]->FindFixBin(max_mass[0]-1e-4);
	  double all = hInvMass[i][0][bin-1]->Integral(low_bin, high_bin);
	  double bkg = hInvMass[i][1][bin-1]->Integral(low_bin, high_bin);
	  TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("UL-LS = %1.0f",all-bkg));
	  t1->AddText(Form("S/B = %4.2f:1",(all-bkg)/bkg));
	  t1->Draw();
	}
      c->cd(6);
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiff_InvMassLSvsUL.pdf",config,data_name[i]));
    }

  //==============================================
  // scale factor
  //==============================================
  TH2F *hJpsiMassVsMuonPt[nData][2];
  TH1F *hMuonPtSide[nData][2];
  TGraphErrors *gScaleFactor[nData];
  double *scales[2];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hJpsiMassVsMuonPt[i][j] = (TH2F*)fdata[i]->Get(Form("mhJpsiMassVsMuonPt_MtdVpdTacDiff_%s_di_mu",type_name[j+1]));
	  hJpsiMassVsMuonPt[i][j]->SetName(Form("%s_%s",data_name[i],hJpsiMassVsMuonPt[i][j]->GetName()));
	  hJpsiMassVsMuonPt[i][j]->Sumw2();
	  hMuonPtSide[i][j] = (TH1F*)hJpsiMassVsMuonPt[i][j]->ProjectionX(Form("%s_TacDiff_SideBand_%s",data_name[i],type_name[j+1]));
	  hMuonPtSide[i][j]->Reset();
	  for(int k=0; k<2; k++)
	    {
	      int low_bin = hJpsiMassVsMuonPt[i][j]->GetYaxis()->FindFixBin(min_mass[k+1]);
	      int up_bin  = hJpsiMassVsMuonPt[i][j]->GetYaxis()->FindFixBin(max_mass[k+1]);
	      TH1F *htmp = (TH1F*)hJpsiMassVsMuonPt[i][j]->ProjectionX(Form("%s_TacDiff_SideBand_%s_%d",data_name[i],type_name[j+1],k),low_bin,up_bin);
	      hMuonPtSide[i][j]->Add(htmp);
	    }
	}

      gScaleFactor[i] = new TGraphErrors(nbins);
      gScaleFactor[i]->SetName(Form("%s_ScaleFactor",data_name[i]));
      gScaleFactor[i]->SetTitle(Form("%s: scale factor for like-sign background (TacDiff);p_{T} (GeV/c)",data_name[i]));
      for(int bin=0; bin<nbins; bin++)
	{
	  int start_bin = hMuonPtSide[i][0]->GetXaxis()->FindBin(lowbins[bin]+1e-4);
	  int end_bin   = hMuonPtSide[i][0]->GetXaxis()->FindBin(upbins[bin]-1e-4);
	  double ul = hMuonPtSide[i][0]->Integral(start_bin,end_bin);
	  double ls = hMuonPtSide[i][1]->Integral(start_bin,end_bin);
	  double scale = ul/ls;
	  double error = scale * sqrt(1/ul+1/ls);
	  double pt = 0.5*(lowbins[bin]+upbins[bin]);
	  if(bin==0) pt = -0.5;
	  gScaleFactor[i]->SetPoint(bin, pt, scale);
	  gScaleFactor[i]->SetPointError(bin, 0, error);
	}
      gScaleFactor[i]->SetMarkerStyle(21);
      c = drawGraph(gScaleFactor[i]);
      TPaveText *t1 = GetPaveText(0.15,0.2,0.2,0.23,0.03);
      t1->AddText(Form("%1.1f < p_{T} < %1.1f",lowbins[0],upbins[0]));
      t1->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiff_ScaleFactor.pdf",config,data_name[i]));
      scales[i] = gScaleFactor[i]->GetY();
    }

  //==============================================
  // MtdVpdTacDiff study
  //==============================================
  THnSparseF *hn[nData]; 
  TH2F *hTacDiffVsTrigUnit[nData][3][nbins];
  TH2F *hTacDiffVsTrigUnitBkg[nData][3][nbins]; // 0 - Muon, 1 - UL, 2 - LS in [2.2,4.0];
  for(int i=0; i<2; i++)
    {
      hn[i] = (THnSparseF*)fdata[i]->Get("mhJpsiMuonTrigEff_di_mu");
      hn[i]->SetName(Form("%s_%s",data_name[i],hn[i]->GetName()));
      hn[i]->GetAxis(1)->SetTitle("#DeltaTacSum");

      for(int bin=0; bin<nbins; bin++)
	{
	  hn[i]->GetAxis(2)->SetRangeUser(lowbins[bin]+1e-4, upbins[bin]-1e-4);
      
	  // unlike-sign signal
	  hn[i]->GetAxis(5)->SetRange(1,1);
	  hTacDiffVsTrigUnitBkg[i][1][bin] = (TH2F*)hn[i]->Projection(1,3);
	  hTacDiffVsTrigUnitBkg[i][1][bin]->SetName(Form("%s_hTacDiffVsTrigUnit_BkgUL_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnitBkg[i][1][bin]->Sumw2();
	  
	  hn[i]->GetAxis(0)->SetRangeUser(min_mass[0]+1e-4,max_mass[0]-1e-4);
	  hTacDiffVsTrigUnit[i][1][bin] = (TH2F*)hn[i]->Projection(1,3);
	  hTacDiffVsTrigUnit[i][1][bin]->SetName(Form("%s_hTacDiffVsTrigUnit_UL_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnit[i][1][bin]->Sumw2();

	  // unlike-sign side-band background
	  hTacDiffVsTrigUnit[i][2][bin] = (TH2F*)hTacDiffVsTrigUnit[i][1][bin]->Clone(Form("%s_hTacDiffVsVzTrigUnit_Bkg_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnit[i][2][bin]->Reset();
	  for(int k=0; k<2; k++)
	    {
	      hn[i]->GetAxis(0)->SetRangeUser(min_mass[k+1]+1e-4,max_mass[k+1]-1e-4);
	      TH2F *h2tmp = (TH2F*)hn[i]->Projection(1,3);
	      h2tmp->SetName(Form("%s_hTacDiffVsVzTrigUnit_UL_SB%d_PtBin%d",data_name[i],k,bin));
	      h2tmp->Sumw2();
	      hTacDiffVsTrigUnit[i][2][bin]->Add(h2tmp);
	    }
	  hn[i]->GetAxis(0)->SetRange(0,-1);

	  // like-sign background
	  hn[i]->GetAxis(5)->SetRange(2,2);
	  hTacDiffVsTrigUnitBkg[i][2][bin] = (TH2F*)hn[i]->Projection(1,3);
	  hTacDiffVsTrigUnitBkg[i][2][bin]->SetName(Form("%s_hTacDiffVsTrigUnit_BkgLS_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnitBkg[i][2][bin]->Sumw2();

	  hn[i]->GetAxis(0)->SetRangeUser(min_mass[1]+1e-4,max_mass[2]-1e-4);
	  TH2F *h2tmp = (TH2F*)hn[i]->Projection(1,3);
	  h2tmp->SetName(Form("%s_hTacDiffVsVzTrigUnit_LSPtBin%d",data_name[i],bin));
	  h2tmp->Sumw2();
	  hTacDiffVsTrigUnit[i][2][bin]->Add(h2tmp);
	  hn[i]->GetAxis(5)->SetRange(0,-1);
	  hn[i]->GetAxis(0)->SetRange(0,-1);
	  double factor = (2*max_mass[2]-2*min_mass[1]-min_mass[2]+max_mass[1])/(max_mass[0]-min_mass[0]);
	  hTacDiffVsTrigUnit[i][2][bin]->Scale(1./factor);

	  // pure muon distribution
	  hTacDiffVsTrigUnit[i][0][bin] = (TH2F*)hTacDiffVsTrigUnit[i][1][bin]->Clone(Form("%s_hTacDiffVsTrigUnit_Muon_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnit[i][0][bin]->Add(hTacDiffVsTrigUnit[i][2][bin], -1*scales[i][bin]);

	  hTacDiffVsTrigUnitBkg[i][0][bin] = (TH2F*)hTacDiffVsTrigUnit[i][0][bin]->Clone(Form("%s_hTacDiffVsTrigUnit_BkgMuon_PtBin%d",data_name[i],bin));
	  //printf("[i] %s: %1.1f < pT < %1.1f, entry = %4.2f\n",data_name[i], lowbins[bin]+1e-4, upbins[bin]-1e-4, hTacDiffVsTrigUnit[i][1][bin]->Integral());
	}
    }

  //==============================================
  // MtdVpdTacDiff vs. trigger unit
  int colors[5] = {1, 2, 4, 6, kGreen+2};
  const int nDatas = 3;
  const int nTrigUnit = 28;
  for(int i=0; i<nData; i++)
    {
      hTacDiffVsTrigUnit[i][0][0]->GetYaxis()->SetRangeUser(850,1000);
      c = draw2D(hTacDiffVsTrigUnit[i][0][0],Form("%s: J/#psi #mu %1.1f < p_{T} < %1.1f GeV/c",data_name[i],lowbins[0],upbins[0]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffVsMod.pdf",config,data_name[i]));
    }

  // compare signal vs. background
  const char* bkg_leg[3] = {"J/#Psi Muon", "UL: [2.2,4.0] GeV/c^{2}", "LS: [2.2,4.0] GeV/c^{2}"};
  const char* bkg_name[3] = {"Muon","ULbkg","LSbkg"};
  TF1 *funcTacDiffTrigUnitBkg[nData][3][nbins][nTrigUnit];
  TH1F *hTacDiffMeanVsTrigUnitBkg[nData][3][nbins];
  double fit_bkg_min[nData] = {915, 895};
  double fit_bkg_max[nData] = {950, 935};
  double fit_min[nDatas] = {890, 870, -50};
  double fit_max[nDatas] = {960, 940, 50};
  for(int i=0; i<nData; i++)
    {
      for(int k=0; k<3; k++)
	{
	  for(int bin=0; bin<nbins; bin++)
	    {
	      hTacDiffMeanVsTrigUnitBkg[i][k][bin] = new TH1F(Form("%s_hTacDiffMeanVsTrigUnit_Bkg%s_PtBin%d",data_name[i],type_name[k], bin),Form("%s: mean of #DeltaTacSum (%1.1f < p_{T} < %1.1f GeV/c);TrigUnit;Mean",data_name[i],lowbins[bin],upbins[bin]),nTrigUnit,1,nTrigUnit+1);
	      c = new TCanvas(Form("%s_TacDiff_Bkg%s_PtBin%d",data_name[i],type_name[k],bin),Form("%s_TacDiff_Bkg%s_PtBin%d",data_name[i],type_name[k],bin),1200,800);
	      c->Divide(6,5);
	      TLegend *leg = new TLegend(0.1,0.2,0.5,0.7);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.09);
	      leg->SetHeader(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]));
	      for(int j=0; j<nTrigUnit; j++)
		{
		  c->cd(j+2);
		  SetPadMargin(gPad,0.15,0.08,0.02,0.1);
		  hTacDiffVsTrigUnitBkg[i][k][bin]->GetYaxis()->SetRangeUser(850,1000);
		  TH1F *htmp = (TH1F*)hTacDiffVsTrigUnitBkg[i][k][bin]->ProjectionY(Form("%s_hTacDiff_%s_TrigUnit%d_PtBin%d",data_name[i],type_name[k],j+1,bin),j+2,j+2);
		  if(htmp->Integral()!=0) htmp->Scale(1./htmp->Integral());
		  htmp->Rebin(5);
		  funcTacDiffTrigUnitBkg[i][k][bin][j] = new TF1(Form("%s_FitTacDiff_%s_TrigUnit%d_PtBin%d",data_name[i],type_name[k],j+1,bin), "gaus", fit_bkg_min[i], fit_bkg_max[i]);
		  if(k==0) funcTacDiffTrigUnitBkg[i][k][bin][j]->SetRange(fit_min[i], fit_max[i]);
		  htmp->Fit(funcTacDiffTrigUnitBkg[i][k][bin][j], "IR0Q");
		  hTacDiffMeanVsTrigUnitBkg[i][k][bin]->SetBinContent(j+1, funcTacDiffTrigUnitBkg[i][k][bin][j]->GetParameter(1));
		  hTacDiffMeanVsTrigUnitBkg[i][k][bin]->SetBinError(j+1, funcTacDiffTrigUnitBkg[i][k][bin][j]->GetParError(1));
		  htmp->SetMarkerStyle(24);
		  htmp->SetMarkerColor(TMath::Power(2,k));
		  htmp->SetLineColor(TMath::Power(2,k));
		  htmp->SetTitle("");
		  ScaleHistoTitle(htmp, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
		  htmp->Draw("");
		  funcTacDiffTrigUnitBkg[i][k][bin][j]->SetLineColor(6);
		  funcTacDiffTrigUnitBkg[i][k][bin][j]->SetLineStyle(2);
		  funcTacDiffTrigUnitBkg[i][k][bin][j]->Draw("sames");
		  if(j==0) 
		    {
		      leg->AddEntry(htmp, bkg_leg[k], "P");
		      leg->AddEntry(funcTacDiffTrigUnitBkg[i][k][bin][j], "Fit", "L");
		    }
		}
	      TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	      t1->Draw();
	      c->cd(1);
	      leg->Draw();
	      TPaveText *t1 = GetPaveText(0.18,0.58,0.75,0.85,0.09);
	      t1->AddText(Form("%s",data_name[i]));
	      t1->Draw();
	      if(savePlot)
		c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffFit_%s_PtBin%d.pdf",config,data_name[i],bkg_name[k],bin));
	    }
	}
    }

  for(int i=0; i<nData; i++)
    {
      c = new TCanvas(Form("%s_TacDiffMeanVsTrigUnit_BkgULvsLS",data_name[i]),Form("%s_TacDiffMeanVsTrigUnit_BkgULvsLS",data_name[i]),1200,600);
      c->Divide(4,2);
      for(int bin=0; bin<nbins; bin++)
	{
	  c->cd(bin+2);
	  SetPadMargin(gPad,0.13,0.13,0.02,0.1);
	  for(int k=0; k<3; k++)
	    {
	      hTacDiffMeanVsTrigUnitBkg[i][k][bin]->SetMarkerStyle(24+k);
	      hTacDiffMeanVsTrigUnitBkg[i][k][bin]->SetMarkerColor(TMath::Power(2,k));
	      hTacDiffMeanVsTrigUnitBkg[i][k][bin]->SetLineColor(TMath::Power(2,k));
	      TH1F *htmp = (TH1F*)hTacDiffMeanVsTrigUnitBkg[i][k][bin]->Clone(Form("%s_tmp",hTacDiffMeanVsTrigUnitBkg[i][k][bin]->GetName()));
	      htmp->GetYaxis()->SetRangeUser(fit_bkg_min[i], fit_bkg_max[i]-10);
	      htmp->SetTitle(";TrigUnit;<#DeltaTacSum>");
	      ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 1.2, 0.04);
	      if(k==0) htmp->Draw();
	      else     htmp->Draw("sames");
	    }
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
	  t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.25,0.4,0.6,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(data_name[i]);
      for(int k=0; k<3; k++) leg->AddEntry(hTacDiffMeanVsTrigUnitBkg[i][k][0], bkg_leg[k], "P");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffMeanVsTrigUnit_MuonVsBkg.pdf",config,data_name[i]));
    }

  // Take the difference between signal and background
  TH1F *hTacDiffMeanVsTrigUnitBkgDiff[nData][2][nbins];
  TH1F *hTacDiffSigBkgDiff[nData];
  for(int i=0; i<nData; i++)
    {
      hTacDiffSigBkgDiff[i] = new TH1F(Form("%s_hTacDiffSigBkgDiff",data_name[i]),"",nPtBins,xPtBins);
      c = new TCanvas(Form("%s_TacDiffMeanVsTrigUnit_BkgDiff",data_name[i]),Form("%s_TacDiffMeanVsTrigUnit_BkgDiff",data_name[i]),1200,600);
      c->Divide(4,2);
      for(int bin=0; bin<nbins; bin++)
	{
	  c->cd(bin+2);
	  SetPadMargin(gPad,0.13,0.13,0.02,0.1);
	  for(int k=0; k<2; k++)
	    {
	      hTacDiffMeanVsTrigUnitBkgDiff[i][k][bin] = (TH1F*)hTacDiffMeanVsTrigUnitBkg[i][k+1][bin]->Clone(Form("%s_hTacDiffMeanVsTrigUnit_BkgDiff%s_PtBin%d",data_name[i],type_name[k+1], bin));
	      hTacDiffMeanVsTrigUnitBkgDiff[i][k][bin]->Add(hTacDiffMeanVsTrigUnitBkg[i][0][bin], -1);
	      TH1F *htmp = (TH1F*)hTacDiffMeanVsTrigUnitBkgDiff[i][k][bin]->Clone(Form("%s_tmp",hTacDiffMeanVsTrigUnitBkgDiff[i][k][bin]->GetName()));
	      htmp->GetYaxis()->SetRangeUser(-10,10);
	      htmp->SetTitle(";TrigUnit;<#DeltaTacSum>^{bkg} - <#DeltaTacSum>^{#mu}");
	      ScaleHistoTitle(htmp, 0.05, 0.9, 0.045, 0.05, 1.2, 0.04);
	      if(k==0) 
		{
		  htmp->Draw();
		}
	      else     
		{
		  TF1 *func = new TF1(Form("funcTacDiff_BkgDiff_%d_%d_%d",i,bin,k),"pol0",1,28);
		  htmp->Fit(func,"IR0Q");
		  func->SetLineColor(htmp->GetLineColor());
		  func->SetLineStyle(2);
		  func->Draw("sames");
		  hTacDiffSigBkgDiff[i]->SetBinContent(bin, func->GetParameter(0));
		  hTacDiffSigBkgDiff[i]->SetBinError(bin, func->GetParError(0));
		  htmp->Draw("sames");
		}
	    }
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
	  t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.25,0.4,0.6,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(data_name[i]);
      for(int k=0; k<2; k++) leg->AddEntry(hTacDiffMeanVsTrigUnitBkgDiff[i][k][0], bkg_leg[k+1], "P");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffMeanVsTrigUnit_MuonMinusBkg.pdf",config,data_name[i]));
    }

  // Fit TacDiff distribution in each trigger unit
  TH1F *hTacDiffTrigUnit[nData][nbins][nTrigUnit];
  TH1F *hTacDiffMeanVsTrigUnit[nDatas][nbins]; 
  TH1F *hTacDiffMeanVsPt[nDatas][nTrigUnit];
  TH1F *hTacDiffWidthVsTrigUnit[nDatas][nbins]; 
  TF1  *funcTacDiffTrigUnit[nDatas][nbins][nTrigUnit];
  for(int i=0; i<nDatas; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  hTacDiffMeanVsTrigUnit[i][bin] = new TH1F(Form("%s_hTacDiffMeanVsTrigUnit_PtBin%d",data_name[i],bin),Form("%s: mean of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);TrigUnit;<#DeltaTacSum>",data_name[i],lowbins[bin],upbins[bin]),nTrigUnit,1,nTrigUnit+1);
	  hTacDiffWidthVsTrigUnit[i][bin] = new TH1F(Form("%s_hTacDiffWidthVsTrigUnit_PtBin%d",data_name[i],bin),Form("%s: width of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);TrigUnit;#sigma(#DeltaTacSum)",data_name[i],lowbins[bin],upbins[bin]),nTrigUnit,1,nTrigUnit+1);
	  for(int j=0; j<nTrigUnit; j++)
	    {
	      funcTacDiffTrigUnit[i][bin][j] = new TF1(Form("%s_FitTacDiff_TrigUnit%d_PtBin%d",data_name[i],j+1,bin),"gaus",fit_min[i],fit_max[i]);
	    }
	}
      for(int j=0; j<nTrigUnit; j++)
	{
	  hTacDiffMeanVsPt[i][j] = new TH1F(Form("%s_hTacDiffMeanVsPt_TrigUnit%d",data_name[i],j+1),Form("%s: mean of MtdVpdTacDiff vs. p_{T};p_{T} (GeV/c);Mean",data_name[i]),nPtBins, xPtBins);
	}
    }

  for(int i=0; i<nData; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  c = new TCanvas(Form("%s_MtdVdpTacDiff_vs_TrigUnit_PtBin%d",data_name[i],bin),Form("%s_MtdVdpTacDiff_vs_TrigUnit_PtBin%d",data_name[i],bin),1200,800);
	  c->Divide(6,5);
	  for(int j=0; j<nTrigUnit; j++)
	    {
	      hTacDiffVsTrigUnit[i][0][bin]->GetYaxis()->SetRangeUser(850,1000);
	      hTacDiffTrigUnit[i][bin][j] = (TH1F*)hTacDiffVsTrigUnit[i][0][bin]->ProjectionY(Form("%s_hTacDiff_TrigUnit%d_PtBin%d",data_name[i],j+1,bin),j+2,j+2);
	      hTacDiffTrigUnit[i][bin][j]->SetMarkerStyle(20);
	      hTacDiffTrigUnit[i][bin][j]->SetMarkerSize(0.8);
	      c->cd(j+2);
	      SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	      TH1F *hFit = (TH1F*)hTacDiffTrigUnit[i][bin][j]->Clone(Form("Fit_%s",hTacDiffTrigUnit[i][bin][j]->GetName()));
	      hFit->Rebin(5);
	      hFit->Fit(funcTacDiffTrigUnit[i][bin][j],"IR0QL");
	      hTacDiffMeanVsTrigUnit[i][bin]->SetBinContent(j+1, funcTacDiffTrigUnit[i][bin][j]->GetParameter(1));
	      hTacDiffMeanVsTrigUnit[i][bin]->SetBinError(j+1, funcTacDiffTrigUnit[i][bin][j]->GetParError(1));
	      hTacDiffWidthVsTrigUnit[i][bin]->SetBinContent(j+1, funcTacDiffTrigUnit[i][bin][j]->GetParameter(2));
	      hTacDiffWidthVsTrigUnit[i][bin]->SetBinError(j+1, funcTacDiffTrigUnit[i][bin][j]->GetParError(2));
	      hTacDiffMeanVsPt[i][j]->SetBinContent(bin, funcTacDiffTrigUnit[i][bin][j]->GetParameter(1));
	      hTacDiffMeanVsPt[i][j]->SetBinError(bin, funcTacDiffTrigUnit[i][bin][j]->GetParError(1));
	      ScaleHistoTitle(hFit, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	      hFit->SetTitle("");
	      hFit->Draw("samesPE");
	      funcTacDiffTrigUnit[i][bin][j]->SetLineColor(4);
	      funcTacDiffTrigUnit[i][bin][j]->SetLineStyle(2);
	      funcTacDiffTrigUnit[i][bin][j]->Draw("sames");
	      TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	      t1->Draw();
	    }
	  c->cd(1);
	  TLegend *leg = new TLegend(0.1,0.3,0.5,0.7);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.09);
	  leg->SetHeader(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]));
	  leg->AddEntry(hTacDiffTrigUnit[i][bin][0],Form("Run15 %s",data_leg[i]),"P");
	  leg->AddEntry(funcTacDiffTrigUnit[i][bin][0],"Fit","L");
	  leg->Draw();
	  if(savePlot)
	    c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffFit_PtBin%d.pdf",config,data_name[i],bin));
	}
    }

  // compare mean vs. trigUnit
  c = new TCanvas("TacDiffMeanVsTrigUnit_pp_vs_pAu", "TacDiffMeanVsTrigUnit_pp_vs_pAu", 800, 600);
  TLegend *leg = new TLegend(0.25,0.65,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(int i=0; i<nData; i++)
    {
      TH1F *htmp = (TH1F*)hTacDiffMeanVsTrigUnit[i][0]->Clone(Form("%s_tmp2",hTacDiffMeanVsTrigUnit[i][0]->GetName()));
      htmp->SetMarkerStyle(20+i);
      htmp->SetMarkerColor(TMath::Power(2,i));
      htmp->SetLineColor(TMath::Power(2,i));
      htmp->GetYaxis()->SetRangeUser(900,950);
      htmp->SetTitle(";TrigUnit;<#DeltaTacSum>");
      ScaleHistoTitle(htmp, 0.04, 0.9, 0.04, 0.04, 1.2, 0.04);
      if(i==0) htmp->Draw();
      else     htmp->Draw("sames");
      leg->AddEntry(htmp, data_name[i], "P");
    }
  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[0],upbins[0]),0.045);
  t1->Draw();
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15ppVspAu%s.TacDiffMeanVsTrigUnit_JpsiMuon.pdf",config));

  // mean vs. pT
  for(int i=0; i<nData; i++)
    {
      c = new TCanvas(Form("%s_MtdVdpTacDiff_vs_Pt",data_name[i]),Form("%s_MtdVdpTacDiff_vs_Pt",data_name[i]),1200,800);
      c->Divide(6,5);
      for(int j=0; j<nTrigUnit; j++)
	{
	  c->cd(j+2);
	  SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	  hTacDiffMeanVsPt[i][j]->SetMarkerStyle(20);
	  TH1F *htmp = (TH1F*)hTacDiffMeanVsPt[i][j]->Clone(Form("%s_clone",hTacDiffMeanVsPt[i][j]->GetName()));
	  ScaleHistoTitle(htmp, 0.06, 0.85, 0.04, 0.06, 1, 0.04);
	  htmp->GetXaxis()->SetRangeUser(1.3,4.5);
	  htmp->GetYaxis()->SetRangeUser(915-i*15,935-i*15);
	  htmp->SetTitle("");
	  htmp->Draw();
	  TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	  t1->Draw();
	}
      c->cd(1);
      TPaveText *t1 = GetPaveText(0.1,0.6,0.4,0.7,0.09);
      t1->AddText(Form("Run15"));
      t1->AddText(Form("%s",data_leg[i]));
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffMeanVsPt_InTrigUnit.pdf",config,data_name[i]));
    }

  // combine all trigger units 
  TH1F *hTacDiffTrigUnitAllRaw[nDatas][nbins];
  TH1F *hTacDiffTrigUnitAll[nDatas][nbins];
  TF1  *funcTacDiffTrigUnitAll[nDatas][nbins];
  TH1F *hTacDiffMeanVsPtAll[nDatas];
  TF1  *funcTacDiffMeanVsPtAll[nDatas];
  for(int i=0; i<nDatas; i++)
    {
      hTacDiffMeanVsPtAll[i] = new TH1F(Form("%s_hTacDiffMeanVsPt_Combined",data_name[i]),Form("%s: mean of MtdVpdTacDiff vs. p_{T};p_{T} (GeV/c);<#DeltaTacSum>",data_name[i]),nPtBins, xPtBins);
      funcTacDiffMeanVsPtAll[i] = new TF1(Form("%s_FitTacDiffVsPt_All",data_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.3,10);
      for(int bin=0; bin<nbins; bin++)
	{
	  hTacDiffTrigUnitAllRaw[i][bin] = new TH1F(Form("%s_hTacDiffAllRaw_PtBin%d",data_name[i],bin),";#DeltaTacSum;",100,-50,50);
	  hTacDiffTrigUnitAll[i][bin] = new TH1F(Form("%s_hTacDiffAll_PtBin%d",data_name[i],bin),";#DeltaTacSum;",100,-50,50);
	  funcTacDiffTrigUnitAll[i][bin] = new TF1(Form("%s_FitTacDiffAll_PtBin%d",data_name[i],bin),"gaus",-15,20);
	  //funcTacDiffTrigUnitAll[i][bin] = new TF1(Form("%s_FitTacDiffAll_PtBin%d",data_name[i],bin),CrystalBall,-50,50,5);
	}
    }
 
  for(int i=0; i<nDatas; i++)
    {
      if(i<nData)
	{
	  for(int bin=0; bin<nbins; bin++)
	    {
	      for(int j=0; j<nTrigUnit; j++)
		{
		  TH1F *htmp = (TH1F*)hTacDiffTrigUnitAll[i][bin]->Clone(Form("%s_%d_tmp",hTacDiffTrigUnitAll[i][bin]->GetName(),j+1));
		  htmp->Reset();
		  int nxbins = htmp->GetNbinsX();
		  double shift = hTacDiffMeanVsPt[i][j]->GetBinContent(0);
		  for(int ibin=1; ibin<=nxbins; ibin++)
		    {
		      int jbin = hTacDiffTrigUnit[i][bin][j]->FindFixBin(htmp->GetBinCenter(ibin) + shift);
		      htmp->SetBinContent(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinContent(jbin));
		      htmp->SetBinError(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinError(jbin));
		    }
		  hTacDiffTrigUnitAll[i][bin]->Add(htmp);

		  htmp = (TH1F*)hTacDiffTrigUnitAllRaw[i][bin]->Clone(Form("%s_%d_tmp",hTacDiffTrigUnitAllRaw[i][bin]->GetName(),j+1));
		  htmp->Reset();
		  nxbins = htmp->GetNbinsX();
		  shift = hTacDiffMeanVsPt[i][0]->GetBinContent(0)+2*(1-i);
		  for(int ibin=1; ibin<=nxbins; ibin++)
		    {
		      int jbin = hTacDiffTrigUnit[i][bin][j]->FindFixBin(htmp->GetBinCenter(ibin) + shift);
		      htmp->SetBinContent(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinContent(jbin));
		      htmp->SetBinError(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinError(jbin));
		    }
		  hTacDiffTrigUnitAllRaw[i][bin]->Add(htmp);
		}
	    }
	}
      else
	{
	  for(int bin=0; bin<nbins; bin++)
	    {
	      hTacDiffTrigUnitAllRaw[i][bin]->Add(hTacDiffTrigUnitAllRaw[0][bin]);
	      hTacDiffTrigUnitAllRaw[i][bin]->Add(hTacDiffTrigUnitAllRaw[1][bin]);

	      hTacDiffTrigUnitAll[i][bin]->Add(hTacDiffTrigUnitAll[0][bin]);
	      hTacDiffTrigUnitAll[i][bin]->Add(hTacDiffTrigUnitAll[1][bin]);
	    }
	}
    }

  // compare TacDiff before and after lineup
  for(int i=0; i<nData; i++)
    {
      TH1F *htmp1[2];
      htmp1[0] = (TH1F*)hTacDiffTrigUnitAllRaw[i][0]->Clone(Form("%s_plot",hTacDiffTrigUnitAllRaw[i][0]->GetName()));
      htmp1[1] = (TH1F*)hTacDiffTrigUnitAll[i][0]->Clone(Form("%s_plot",hTacDiffTrigUnitAll[i][0]->GetName()));

      for(j=0; j<2; j++)
	{
	  //htmp1[j]->Rebin(2);
	  htmp1[j]->SetMarkerStyle(20+j*4);
	  htmp1[j]->SetMarkerColor(j+1);
	  htmp1[j]->SetLineColor(j+1);
	}
      c = draw1D(htmp1[0],Form("%s: \DeltaTacSum distribution",data_name[i]));
      htmp1[1]->Draw("sames");
      TLegend *leg = new TLegend(0.2,0.65,0.4,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader("p_{T} > 1.3 GeV/c");
      leg->AddEntry(htmp1[0],"Raw","P");
      leg->AddEntry(htmp1[1],"Lined up","P");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffComp_Lineup.pdf",config,data_name[i]));
    }

  // Fit TacDiff after lineup vs. pT
  for(int i=0; i<nDatas; i++)
    {
      c = new TCanvas(Form("%s_MtdVdpTacDiffAll_vs_Pt",data_name[i]),Form("%s_MtdVdpTacDiffAll_vs_Pt",data_name[i]),1200,600);
      c->Divide(4,2);
 
      for(int bin=0; bin<nbins; bin++)
	{
	  c->cd(bin+2);
	  TH1F *hFit = (TH1F*)hTacDiffTrigUnitAll[i][bin]->Clone(Form("%s_fit",hTacDiffTrigUnitAll[i][bin]->GetName()));
	  hFit->Rebin(2);
	  hFit->SetMarkerStyle(25);
	  hFit->Fit(funcTacDiffTrigUnitAll[i][bin],"IR0QL");
	  hTacDiffMeanVsPtAll[i]->SetBinContent(bin, funcTacDiffTrigUnitAll[i][bin]->GetParameter(1));
	  hTacDiffMeanVsPtAll[i]->SetBinError(bin, funcTacDiffTrigUnitAll[i][bin]->GetParError(1));
	  hFit->SetTitle("");
	  hFit->Draw("PE");
	  funcTacDiffTrigUnitAll[i][bin]->SetLineColor(4);
	  funcTacDiffTrigUnitAll[i][bin]->SetLineStyle(2);
	  funcTacDiffTrigUnitAll[i][bin]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
	  t1->Draw();
	}
      c->cd(1);
      hTacDiffMeanVsPtAll[i]->SetMarkerStyle(21);
      TH1F *hFit = (TH1F*)hTacDiffMeanVsPtAll[i]->Clone(Form("%s_fit",hTacDiffMeanVsPtAll[i]->GetName()));
      hFit->SetTitle("");
      funcTacDiffMeanVsPtAll[i]->SetParameters(4,2,2);
      hFit->Fit(funcTacDiffMeanVsPtAll[i],"IR0QL");
      hFit->GetYaxis()->SetRangeUser(-3,10);
      hFit->Draw("");
      funcTacDiffMeanVsPtAll[i]->SetLineColor(2);
      funcTacDiffMeanVsPtAll[i]->SetLineStyle(2);
      funcTacDiffMeanVsPtAll[i]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("Run15 %s",data_leg[i]),0.06);
      t1->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffMeanVsPt_Combined.pdf",config,data_name[i]));
    }

  // Line up TacDiff distribution
  TH1F *hTacDiffTrigUnitLineup[nDatas][nbins][nTrigUnit];
  TH1F *hTacDiffPtLineup[nDatas][nbins];
  for(int i=0; i<nDatas; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  hTacDiffPtLineup[i][bin] = new TH1F(Form("%s_hTacDiffLineup_PtBin%d",data_name[i],bin),"",100,-50,50);
	  for(int j=0; j<nTrigUnit; j++)
	    {
	      hTacDiffTrigUnitLineup[i][bin][j] = new TH1F(Form("%s_hTacDiffLineup_TrigUnit%d_PtBin%d",data_name[i],j+1,bin),"",100,-50,50);
	      if(i<nData)
		{
		  int nxbins = hTacDiffTrigUnitLineup[i][bin][j]->GetNbinsX();
		  double shift = hTacDiffMeanVsTrigUnit[i][bin]->GetBinContent(j+1);
		  for(int ibin=1; ibin<=nxbins; ibin++)
		    {
		      int jbin = hTacDiffTrigUnit[i][bin][j]->FindFixBin(hTacDiffTrigUnitLineup[i][bin][j]->GetBinCenter(ibin) + shift);
		      hTacDiffTrigUnitLineup[i][bin][j]->SetBinContent(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinContent(jbin));
		      hTacDiffTrigUnitLineup[i][bin][j]->SetBinError(ibin, hTacDiffTrigUnit[i][bin][j]->GetBinError(jbin));
		    }
		}
	      else
		{
		  hTacDiffTrigUnitLineup[i][bin][j]->Add(hTacDiffTrigUnitLineup[0][bin][j]);
		  hTacDiffTrigUnitLineup[i][bin][j]->Add(hTacDiffTrigUnitLineup[1][bin][j]);
		}
	      hTacDiffPtLineup[i][bin]->Add(hTacDiffTrigUnitLineup[i][bin][j]);
	    }
	}
    }
	    
  // fit pp+pAu after line up
  for(int bin=0; bin<nbins; bin++)
    {
      c = new TCanvas(Form("%s_MtdVdpTacDiff_vs_TrigUnit_PtBin%d",data_name[2],bin),Form("%s_MtdVdpTacDiff_vs_TrigUnit_PtBin%d",data_name[2],bin),1200,800);
      c->Divide(6,5);
      for(int j=0; j<nTrigUnit; j++)
	{
	  c->cd(j+1);
	  SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	  TH1F *hFit = (TH1F*)hTacDiffTrigUnitLineup[2][bin][j]->Clone(Form("Fit_%s",hTacDiffTrigUnitLineup[2][bin][j]->GetName()));
	  hFit->Rebin(5);
	  hFit->SetMarkerStyle(20);
	  hFit->SetMarkerSize(0.8);
	  hFit->Fit(funcTacDiffTrigUnit[2][bin][j],"IR0QL");
	  hTacDiffMeanVsTrigUnit[2][bin]->SetBinContent(j+1, funcTacDiffTrigUnit[2][bin][j]->GetParameter(1));
	  hTacDiffMeanVsTrigUnit[2][bin]->SetBinError(j+1, funcTacDiffTrigUnit[2][bin][j]->GetParError(1));
	  hTacDiffWidthVsTrigUnit[2][bin]->SetBinContent(j+1, funcTacDiffTrigUnit[2][bin][j]->GetParameter(2));
	  hTacDiffWidthVsTrigUnit[2][bin]->SetBinError(j+1, funcTacDiffTrigUnit[2][bin][j]->GetParError(2));
	  ScaleHistoTitle(hFit, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	  hFit->SetTitle("");
	  hFit->Draw("samesPE");
	  funcTacDiffTrigUnit[2][bin][j]->SetLineColor(2);
	  funcTacDiffTrigUnit[2][bin][j]->SetLineStyle(2);
	  funcTacDiffTrigUnit[2][bin][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("TrigUnit = %d",j+1),0.09);
	  t1->Draw();

	  for(int k=0; k<2; k++)
	    {
	      TH1F *h1tmp = (TH1F*)hTacDiffTrigUnitLineup[k][bin][j]->Clone(Form("%s_clone",hTacDiffTrigUnitLineup[k][bin][j]->GetName()));
	      h1tmp->SetMarkerStyle(24);
	      h1tmp->SetMarkerColor(4+k*2);
	      h1tmp->Rebin(5);
	      h1tmp->Draw("samesPE");
	    }
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/%s%s_TacDiffFit_PtBin%d.pdf",config,data_name[2],bin));
    }
 
  // compare sigma bewteen pp and pAu
  TCanvas *c = new TCanvas("MtdVdpTacDiffSigma_pp_vs_pAu","MtdVdpTacDiffSigma_pp_vs_pAu",1200,600);
  c->Divide(3,2);
  int marker_style[nDatas] = {24,24,25};
  int marker_color[nDatas] = {2,4,1};
  for(int bin=0; bin<6; bin++)
    {
      c->cd(bin+1);
      for(int i=0; i<nDatas; i++)
	{
	  hTacDiffWidthVsTrigUnit[i][bin]->SetMarkerStyle(marker_style[i]);
	  hTacDiffWidthVsTrigUnit[i][bin]->SetMarkerColor(marker_color[i]);
	  hTacDiffWidthVsTrigUnit[i][bin]->SetLineColor(marker_color[i]);
	  hTacDiffWidthVsTrigUnit[i][bin]->GetYaxis()->SetRangeUser(5,18);
	  hTacDiffWidthVsTrigUnit[i][bin]->SetTitle("");
	  ScaleHistoTitle(hTacDiffWidthVsTrigUnit[i][bin], 0.05, 0.9, 0.045, 0.055, 0.8, 0.045);
	  if(i==0) hTacDiffWidthVsTrigUnit[i][bin]->Draw();
	  else     hTacDiffWidthVsTrigUnit[i][bin]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c",lowbins[bin],upbins[bin]),0.06);
      t1->Draw();
      if(bin==0)
	{
	  TLegend *leg = new TLegend(0.4,0.65,0.7,0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->SetHeader("Run15");
	  for(int i=0; i<nDatas; i++)
	    {
	      leg->AddEntry(hTacDiffWidthVsTrigUnit[i][bin],data_leg[i],"P");
	    }
	  leg->Draw();
	}
    }
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run14_AuAu200/ana_MtdTrigEff/Run15ppVspAu%s.TacDiffSigmaVsTrigUnit.pdf",config));
  
  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run15.ppAu200.MtdTrigEff.root","recreate");
      for(int i=0; i<nData; i++)
	{
	  for(int bin=0; bin<nbins; bin++)
	    {
	      hTacDiffVsTrigUnit[i][0][bin]->Write("",TObject::kOverwrite);
	      hTacDiffMeanVsTrigUnit[i][bin]->Write("",TObject::kOverwrite);
	      hTacDiffWidthVsTrigUnit[i][bin]->Write("",TObject::kOverwrite);
	      hTacDiffMeanVsTrigUnitBkgDiff[i][1][bin]->Write("",TObject::kOverwrite);
	    }
	  hTacDiffSigBkgDiff[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void trigElecEff(const int savePlot = 1)
{
  if(year==2014)
    {
      f = TFile::Open("output/Run14.AuAu200.MB.TrigElecEff.root","read");
    }
  else
    {
      printf("[e] No available input file!\n");
      return;
    }

  THnSparseF *hnTrigEff = (THnSparseF*)f->Get("mhMtdTrigElecEff_mb");
  const int nbins = 11;
  const double xbins[nbins+1] = {0,1,1.5,2,2.5,3,3.5,4,5,6,8,10};
  
  TList *list = new TList;
  // Efficiency vs. dTof
  const int nDtof = 4;
  const double dtof_cut[nDtof] = {1, 0.75, 0.5, 0.25};
  TH1F *hMuonPtDtof[nDtof][3];
  TH1F *hMuonEffDtof[nDtof][3];
  for(int i=0; i<nDtof; i++)
    {
      hnTrigEff->GetAxis(3)->SetRangeUser(-3,dtof_cut[i]);
      hMuonPtDtof[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][0]->SetName(Form("hMuonPtMatch_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtDtof[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][1]->SetName(Form("hMuonPtQT_DtofCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtDtof[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][2]->SetName(Form("hMuonPtHigh2_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(3)->SetRange(0,-1);
    }

  TString legName1[nDtof];
  for(int i=0; i<nDtof; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtDtof[i][j]->GetName()),xbins);
	  hMuonEffDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Clone(Form("hMuonEffDtof_%d_%d",i,j));
	  hMuonEffDtof[i][j]->Sumw2();
	  hMuonEffDtof[i][j]->Divide(hMuonPtDtof[i][0]);
	  hMuonEffDtof[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffDtof[i][j]); 
	}
      legName1[i] = Form("#Deltatof < %2.2f ns",dtof_cut[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Dtof",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName1,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareDtof.pdf",run_type));

 
  // Efficiency vs. centrality
  // use dtof < 1 ns cut
  TH1F *hMuonPtCent[nCentBins][3];
  TH1F *hMuonEffCent[nCentBins][3];
  for(int i=0; i<nCentBins; i++)
    {
      hnTrigEff->GetAxis(5)->SetRange(centBins_low[i], centBins_high[i]);
      hMuonPtCent[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][0]->SetName(Form("hMuonPtMatch_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtCent[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][1]->SetName(Form("hMuonPtQT_%s",cent_Title[i]));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtCent[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][2]->SetName(Form("hMuonPtHigh2_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(5)->SetRange(0,-1);
    }

  TString legName2[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtCent[i][j]->GetName()),xbins);
	  hMuonEffCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Clone(Form("hMuonEff_%s_%d",cent_Title[i],j));
	  hMuonEffCent[i][j]->Sumw2();
	  hMuonEffCent[i][j]->Divide(hMuonPtCent[i][0]);
	  hMuonEffCent[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffCent[i][j]); 
	}
      legName2[i] = Form("%s%%",cent_Name[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Cent",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName2,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareCent.pdf",run_type));

  // Efficiency vs. TPC vz
  // use dtof < 1 ns cut and 0-60%
  hnTrigEff->GetAxis(5)->SetRange(5,16);
  const int nTpcVz = 6;
  const double tpcvz_cut[7] = {-100,-30,-5,0,5,30,100};
  TH1F *hMuonPtTpcVz[nTpcVz][3];
  TH1F *hMuonEffTpcVz[nTpcVz][3];
  for(int i=0; i<nTpcVz; i++)
    {
      hnTrigEff->GetAxis(4)->SetRangeUser(tpcvz_cut[i]+1e-4,tpcvz_cut[i+1]-1e-4);
      hMuonPtTpcVz[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][0]->SetName(Form("hMuonPtMatch_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtTpcVz[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][1]->SetName(Form("hMuonPtQT_TpcVzCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtTpcVz[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][2]->SetName(Form("hMuonPtHigh2_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(4)->SetRange(0,-1);

      for(int j=0; j<3; j++)
	{
	  hMuonPtTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtTpcVz[i][j]->GetName()),xbins);
	}
    }

  TString legName3[nTpcVz];
  for(int i=3; i<6; i++)
    {
      TH1F *hRef = (TH1F*)hMuonPtTpcVz[i][0]->Clone(Form("%s_clone",hMuonPtTpcVz[i][0]->GetName()));
      hRef->Add(hMuonPtTpcVz[5-i][0]);
      for(int j=0; j<3; j++)
	{
	  hMuonEffTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Clone(Form("hMuonEff_%d_%d",i,j));
	  hMuonEffTpcVz[i][j]->Sumw2();
	  hMuonEffTpcVz[i][j]->Add(hMuonPtTpcVz[5-i][j]);
	  hMuonEffTpcVz[i][j]->Divide(hRef);
	  hMuonEffTpcVz[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffTpcVz[i][j]); 
	}
      legName3[i-3] = Form("%1.0f < |vz| < %1.0f cm",tpcvz_cut[i],tpcvz_cut[i+1]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_TpcVz",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName3,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareTpcVz.pdf",run_type));
  hnTrigEff->GetAxis(5)->SetRange(0,-1);

  // Final efficiency: dtof < 1 ns, 0-60%, |tpcVz| < 100 cm
  TH1F *hTrigElecEff = (TH1F*)hMuonEffCent[0][2]->Clone(Form("%s_TrigElecEff",run_type));
  TF1 *func = new TF1(Form("%s_FitFunc",hTrigElecEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,10);
  func->SetParLimits(0,0,1);
  hTrigElecEff->Fit(func,"R0Q");
  hTrigElecEff->GetYaxis()->SetRangeUser(0.85,1.1);
  c = draw1D(hTrigElecEff,Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type));
  func->SetLineColor(2);
  func->SetLineStyle(2);
  func->Draw("sames");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_Fit.pdf",run_type));
}
