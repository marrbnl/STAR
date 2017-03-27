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
  anaTrigEff();
}

//================================================
void anaTrigEff(const int savePlot = 0)
{
  const char* config = "";
  const int nData = 2;
  const char *data_name[3] = {"Run15_pp200", "Run15_pAu200","Run15_ppAu200"};
  const char* type_name[3] = {"Muon","UL", "LS"};
  const int nbins = 7;
  const double lowbins[nbins] = {1.3, 1.3, 1.5, 2.0, 2.5, 3.0, 5.0};
  const double upbins[nbins]  = {10,  1.5, 2.0, 2.5, 3.0, 5.0, 10.0};

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
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s_TacDiff_InvMassLSvsUL.pdf",config,data_name[i]));
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
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s_TacDiff_ScaleFactor.pdf",config,data_name[i]));
      scales[i] = gScaleFactor[i]->GetY();
    }

  //==============================================
  // MtdVpdTacDiff study
  //==============================================
  THnSparseF *hn[nData]; 
  TH2F *hTacDiffVsTrigUnit[nData][3][nbins];
  
  for(int i=0; i<2; i++)
    {
      hn[i] = (THnSparseF*)fdata[i]->Get("mhJpsiMuonTrigEff_di_mu");
      hn[i]->SetName(Form("%s_%s",data_name[i],hn[i]->GetName()));

      for(int bin=0; bin<nbins; bin++)
	{
	  hn[i]->GetAxis(2)->SetRangeUser(lowbins[bin]+1e-4, upbins[bin]-1e-4);
      
	  // unlike-sign signal
	  hn[i]->GetAxis(5)->SetRange(1,1);
	  hn[i]->GetAxis(0)->SetRangeUser(min_mass[0]+1e-4,max_mass[0]-1e-4);
	  hTacDiffVsTrigUnit[i][1][bin] = (TH2F*)hn[i]->Projection(1,3);
	  hTacDiffVsTrigUnit[i][1][bin]->SetName(Form("%s_hTacDiffVsVzTrigUnit_UL_PtBin%d",data_name[i],bin));
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
	      hn[i]->GetAxis(0)->SetRange(0,-1);
	    }
	  hn[i]->GetAxis(0)->SetRange(0,-1);

	  // like-sign background
	  hn[i]->GetAxis(5)->SetRange(2,2);
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
	  hTacDiffVsTrigUnit[i][0][bin] = (TH2F*)hTacDiffVsTrigUnit[i][1][bin]->Clone(Form("%s_hTacDiffVsVzTrigUnit_Muon_PtBin%d",data_name[i],bin));
	  hTacDiffVsTrigUnit[i][0][bin]->Add(hTacDiffVsTrigUnit[i][2][bin], -1*scales[i][bin]);
	  //printf("[i] %s: %1.1f < pT < %1.1f, entry = %4.2f\n",data_name[i], lowbins[bin]+1e-4, upbins[bin]-1e-4, hTacDiffVsTrigUnit[i][1][bin]->Integral());
	}
    }

  //==============================================
  // MtdVpdTacDiff vs. trig Module
  int colors[5] = {1, 2, 4, 6, kGreen+2};
  for(int i=0; i<nData; i++)
    {
      hTacDiffVsTrigUnit[i][0][0]->GetYaxis()->SetRangeUser(850,1000);
      c = draw2D(hTacDiffVsTrigUnit[i][0][0],Form("%s: J/#psi #mu %1.1f < p_{T} < %1.1f GeV/c",data_name[i],lowbins[0],upbins[0]));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s_TacDiffVsMod.pdf",config,data_name[i]));

      // c =  new TCanvas(Form("%s_TacDiffInMod",data_name[i]),Form("%s_TacDiffInMod",data_name[i]),1200,800);
      // c->Divide(6,5);
      // for(int j=0; j<30; j++)
      // 	{
      // 	  c->cd(j+1);
      // 	  for(int k=0; k<5; k++)
      // 	    {
      // 	      TH1F *h1tmp = (TH1F*)hTacDiffVsTrigUnit[i][0][0]->ProjectionY(Form("%s_hTacDiff_BL%d_Mod%d",data_name[i],j+1,k+1),j*5+k+1, j*5+k+1);
      // 	      h1tmp->Rebin(5);
      // 	      if(h1tmp->Integral()>0) h1tmp->Scale(1./h1tmp->Integral());
      // 	      h1tmp->SetLineColor(colors[k]);
      // 	      h1tmp->GetYaxis()->SetRangeUser(0,0.3);
      // 	      if(k==0) h1tmp->Draw("HIST");
      // 	      else     h1tmp->Draw("samesHIST");
      // 	    }
      // 	  TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.08);
      // 	  t1->Draw();
      // 	}
      // if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run15_pp200/ana_MtdTrigEff/%s%s_TacDiffInMod.pdf",config,data_name[i]));
    }
  return;

  // Fit TacDiff distribution in each BL
  const int nDatas = 3;
  TH1F *hTacDiffBL[nData][nbins][30];
  TH1F *hTacDiffBLMean[nDatas][nbins]; 
  TH1F *hTacDiffBLWidth[nDatas][nbins]; 
  TF1  *funcTacDiffBL[nDatas][nbins][30];
  double fit_min[nDatas] = {900, 860, -50};
  double fit_max[nDatas] = {960, 940, 50};
  for(int i=0; i<nDatas; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  hTacDiffBLMean[i][bin] = new TH1F(Form("%s_hTacDiffBLMean_PtBin%d",data_name[i],bin),Form("%s: mean of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);BL;Mean",data_name[i],lowbins[bin],upbins[bin]),30,0,30);
	  hTacDiffBLWidth[i][bin] = new TH1F(Form("%s_hTacDiffBLWidth_PtBin%d",data_name[i],bin),Form("%s: width of MtdVpdTacDiff (%1.1f < p_{T} < %1.1f GeV/c);BL;#sigma",data_name[i],lowbins[bin],upbins[bin]),30,0,30);
	  for(int j=0; j<30; j++)
	    {
	      funcTacDiffBL[i][bin][j] = new TF1(Form("%s_FitTacDif_BL%d_PtBin%d",data_name[i],j+1,bin),"gaus",fit_min[i],fit_max[i]);
	    }
	}
    }

  for(int i=0; i<nData; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  c = new TCanvas(Form("%s_MtdVdpTacDiff_vs_BL_PtBin%d",data_name[i],bin),Form("%s_MtdVdpTacDiff_vs_BL_PtBin%d",data_name[i],bin),1200,800);
	  c->Divide(6,5);
	  for(int j=0; j<30; j++)
	    {
	      hTacDiffVsTrigUnit[i][0][bin]->GetYaxis()->SetRangeUser(850,1000);
	      hTacDiffBL[i][bin][j] = (TH1F*)hTacDiffVsTrigUnit[i][0][bin]->ProjectionY(Form("%s_hTacDiff_BL%d_PtBin%d",data_name[i],j+1,bin),j*5+1, j*5+5);
	      c->cd(j+1);
	      SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	      TH1F *hFit = (TH1F*)hTacDiffBL[i][bin][j]->Clone(Form("Fit_%s",hTacDiffBL[i][bin][j]->GetName()));
	      hFit->Rebin(5);
	      hFit->SetMarkerStyle(20);
	      hFit->SetMarkerSize(0.8);
	      hFit->Fit(funcTacDiffBL[i][bin][j],"IR0QL");
	      hTacDiffBLMean[i][bin]->SetBinContent(j+1, funcTacDiffBL[i][bin][j]->GetParameter(1));
	      hTacDiffBLMean[i][bin]->SetBinError(j+1, funcTacDiffBL[i][bin][j]->GetParError(1));
	      hTacDiffBLWidth[i][bin]->SetBinContent(j+1, funcTacDiffBL[i][bin][j]->GetParameter(2));
	      hTacDiffBLWidth[i][bin]->SetBinError(j+1, funcTacDiffBL[i][bin][j]->GetParError(2));
	      ScaleHistoTitle(hFit, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	      hFit->SetTitle("");
	      hFit->Draw("samesPE");
	      funcTacDiffBL[i][bin][j]->SetLineColor(4);
	      funcTacDiffBL[i][bin][j]->SetLineStyle(2);
	      funcTacDiffBL[i][bin][j]->Draw("sames");
	      TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.09);
	      t1->Draw();
	    }
	}
    }

  // Line up TacDiff distribution
  TH1F *hTacDiffBlLineup[nDatas][nbins][30];
  TH1F *hTacDiffPtLineup[nDatas][nbins];
  for(int i=0; i<nDatas; i++)
    {
      for(int bin=0; bin<nbins; bin++)
	{
	  hTacDiffPtLineup[i][bin] = new TH1F(Form("%s_hTacDiffLineup_PtBin%d",data_name[i],bin),"",100,-50,50);
	  for(int j=0; j<30; j++)
	    {
	      hTacDiffBlLineup[i][bin][j] = new TH1F(Form("%s_hTacDiffLineup_BL%d_PtBin%d",data_name[i],j+1,bin),"",100,-50,50);
	      if(i<nData)
		{
		  int nxbins = hTacDiffBlLineup[i][bin][j]->GetNbinsX();
		  double shift = hTacDiffBLMean[i][bin]->GetBinContent(j+1);
		  for(int ibin=1; ibin<=nxbins; ibin++)
		    {
		      int jbin = hTacDiffBL[i][bin][j]->FindFixBin(hTacDiffBlLineup[i][bin][j]->GetBinCenter(ibin) + shift);
		      hTacDiffBlLineup[i][bin][j]->SetBinContent(ibin, hTacDiffBL[i][bin][j]->GetBinContent(jbin));
		      hTacDiffBlLineup[i][bin][j]->SetBinError(ibin, hTacDiffBL[i][bin][j]->GetBinError(jbin));
		    }
		}
	      else
		{
		  hTacDiffBlLineup[i][bin][j]->Add(hTacDiffBlLineup[0][bin][j]);
		  hTacDiffBlLineup[i][bin][j]->Add(hTacDiffBlLineup[1][bin][j]);
		}
	      hTacDiffPtLineup[i][bin]->Add(hTacDiffBlLineup[i][bin][j]);
	    }
	}
    }
	    
  // fit pp+pAu after line up
  for(int bin=0; bin<nbins; bin++)
    {
      c = new TCanvas(Form("%s_MtdVdpTacDiff_vs_BL_PtBin%d",data_name[2],bin),Form("%s_MtdVdpTacDiff_vs_BL_PtBin%d",data_name[2],bin),1200,800);
      c->Divide(6,5);
      for(int j=0; j<30; j++)
	{
	  c->cd(j+1);
	  SetPadMargin(gPad,0.15,0.05,0.05,0.1);
	  TH1F *hFit = (TH1F*)hTacDiffBlLineup[2][bin][j]->Clone(Form("Fit_%s",hTacDiffBlLineup[2][bin][j]->GetName()));
	  hFit->Rebin(5);
	  hFit->SetMarkerStyle(20);
	  hFit->SetMarkerSize(0.8);
	  hFit->Fit(funcTacDiffBL[2][bin][j],"IR0QL");
	  hTacDiffBLMean[2][bin]->SetBinContent(j+1, funcTacDiffBL[2][bin][j]->GetParameter(1));
	  hTacDiffBLMean[2][bin]->SetBinError(j+1, funcTacDiffBL[2][bin][j]->GetParError(1));
	  hTacDiffBLWidth[2][bin]->SetBinContent(j+1, funcTacDiffBL[2][bin][j]->GetParameter(2));
	  hTacDiffBLWidth[2][bin]->SetBinError(j+1, funcTacDiffBL[2][bin][j]->GetParError(2));
	  ScaleHistoTitle(hFit, 0.08, 0.85, 0.06, 0.08, 1, 0.06);
	  hFit->SetTitle("");
	  hFit->Draw("samesPE");
	  funcTacDiffBL[2][bin][j]->SetLineColor(4);
	  funcTacDiffBL[2][bin][j]->SetLineStyle(2);
	  funcTacDiffBL[2][bin][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("BL = %d",j+1),0.09);
	  t1->Draw();
	}
    }
 
  // Fit TacDiff vs. pT 
  
  
  return;
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/JpsiMuon_MtdVpdTacDiffBL_Fit.pdf",run_type));
  hTacDiffBLMean->SetMarkerStyle(21);
  hTacDiffBLMean->GetYaxis()->SetRangeUser(900,950);
  c = draw1D(hTacDiffBLMean);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/JpsiMuon_MtdVpdTacDiffBL_Mean.pdf",run_type));
  hTacDiffBLWidth->SetMarkerStyle(21);
  hTacDiffBLWidth->GetYaxis()->SetRangeUser(0,20);
  c = draw1D(hTacDiffBLWidth);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/JpsiMuon_MtdVpdTacDiffBL_Sigma.pdf",run_type));

  return;


  // for(int i=0; i<2; i++)
  //   {
  //     hMtdVdpTacDiffBL[i]->Scale(1./hMtdVdpTacDiffBL[i]->Integral());
  //     hMtdVdpTacDiffBL[i]->SetLineColor(i+1);
  //     hMtdVdpTacDiffBL[i]->SetMarkerStyle(21);
  //     if(i==0) c = draw1D(hMtdVdpTacDiffBL[i]);
  //     else     hMtdVdpTacDiffBL[i]->Draw("samesHIST");
  //   }
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
