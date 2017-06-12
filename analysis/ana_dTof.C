const int year = YEAR;
TFile *f = 0x0;

//================================================
void ana_Dtof()
{
  gStyle->SetOptStat(0);
  
  //scanDtofCut();
  DtofCutEff();
}

//================================================
void DtofCutEff(const int savePlot = 0)
{
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double dtofCut = 0.75;

  f = TFile::Open("Output/Run14_AuAu200.jpsi.root","read");
  TH2F *hInvMassVsPt[2];
  hInvMassVsPt[0] = (TH2F*)f->Get("mhTaP_dtof_di_mu_All");
  hInvMassVsPt[1] = (TH2F*)f->Get("mhTaP_dtof_di_mu_Acc");
  TH1F *hInvMass[2][nbins];
  for(int i=0; i<2; i++)
    {
      TCanvas *c = new TCanvas(Form("FitInvMass_%d",i), Form("FitInvMass_%d",i), 1100, 700);
      c->Divide(3,2);
      for(int j=0; j<nbins; j++)
	{
	  int low_bin = hInvMassVsPt[i]->GetXaxis()->FindFixBin(xbins[j]+1e-4);
	  int up_bin  = hInvMassVsPt[i]->GetXaxis()->FindFixBin(xbins[j+1]-1e-4);
	  hInvMass[i][j] = (TH1F*)hInvMassVsPt[i]->ProjectionY(Form("hInvMass_bin%d_%d",j+1,i),low_bin, up_bin);
	  c->cd(j+1);
	  hInvMass[i][j]->Sumw2();
	  hInvMass[i][j]->Rebin(3);
	  hInvMass[i][j]->SetMarkerStyle(24);
	  hInvMass[i][j]->Draw();
	}
    }
}

//================================================
void scanDtofCut(const int savePlot = 0)
{
  f = TFile::Open("Output/Run14_AuAu200.JpsiMuon.root","read");
  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];
  const int nScan = 5;
  const double dtofCut[nScan] = {1.0, 0.75, 0.5, 0.25, 0};
  TH1F *hJpsiMassUL[nPtBins][nScan];
  TH1F *hJpsiMassLS[nPtBins][nScan];
  TH1F *hJpsiMass[nPtBins][nScan];
  TH1F *hBestDtofCut[2]; // 0 - bin counting;1 - fitting
  TH1F *hJpsiSigVsDtof[nPtBins][2]; 
  TH1F *hJpsiSigVsPt[nScan][2];
  TH1F *hJpsiErrVsDtof[nPtBins][2];
  const char *method[2] = {"BinCount","Fit"};
  for(int k=0; k<2; k++)
    {
      hBestDtofCut[k] = new TH1F(Form("hBestDtofCut_%s",method[k]),";p_{T} (GeV/c);#Deltatof (ns)",nbins,xbins);
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiSigVsDtof[i][k] = new TH1F(Form("JpsiSigVsDtof_pt%1.0f-%1.0f_%s",ptBins_low[i],ptBins_high[i],method[k]),";#Deltatof cut (ns);Significance",nScan, -0.125,1.125);
	  hJpsiErrVsDtof[i][k] = new TH1F(Form("JpsiErrVsDtof_pt%1.0f-%1.0f_%s",ptBins_low[i],ptBins_high[i],method[k]),";#Deltatof cut (ns);Stat. Err.",nScan, -0.125,1.125);
	}
      for(int j=0; j<nScan; j++)
	{
	  hJpsiSigVsPt[j][k] = new TH1F(Form("hJpsiSigVsPt_dtof%2.2f_%s",dtofCut[j],method[k]),";p_{T} (GeV/c);Significance",nbins,xbins);
	}
    }

  
  THnSparseF* hn = (THnSparseF*)f->Get("mhJpsiMassVsDtof_di_mu");
  const int rebin = 5;
  for(int i=0; i<nPtBins; i++)
    {
      hn->GetAxis(1)->SetRangeUser(ptBins_low[i]+1e-4,ptBins_high[i]-1e-4); // pt cut
      for(int j=0; j<nScan; j++)
	{
	  hn->GetAxis(2)->SetRangeUser(-3, dtofCut[j]-1e-4); // dtof cut

	  hn->GetAxis(3)->SetRange(1,1); // unlike-sign
	  hJpsiMassUL[i][j] = (TH1F*)hn->Projection(0);
	  hJpsiMassUL[i][j]->SetName(Form("hJpsiMassUL_pt%1.0f-%1.0f_dtof%1.2f",ptBins_low[i],ptBins_high[i],dtofCut[j]));
	  hJpsiMassUL[i][j]->Sumw2();
	  hJpsiMassUL[i][j]->Rebin(rebin);

	  hn->GetAxis(3)->SetRange(2,2); // like-sign
	  hJpsiMassLS[i][j] = (TH1F*)hn->Projection(0);
	  hJpsiMassLS[i][j]->SetName(Form("hJpsiMassLS_pt%1.0f-%1.0f_dtof%1.2f",ptBins_low[i],ptBins_high[i],dtofCut[j]));
	  hJpsiMassLS[i][j]->Sumw2();
	  hJpsiMassLS[i][j]->Rebin(rebin);

	  hn->GetAxis(3)->SetRange(0,-1);
	  hn->GetAxis(2)->SetRange(0,-1);
	}
      hn->GetAxis(1)->SetRange(0,-1);
    }

  // bin counting method
  const double min_mass[3] = {3.0, 2.5, 3.3};
  const double max_mass[3] = {3.2, 2.8, 3.6};
  double counts[nPtBins][nScan][2][3] = {0};
  TCanvas *cBinCount[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      cBinCount[i] = new TCanvas(Form("cBinCount_%d",i),Form("cBinCount_pt%1.0f-%1.0f",ptBins_low[i],ptBins_high[i]),1100,700);
      cBinCount[i]->Divide(3,2);
      for(int j=0; j<nScan; j++)
	{
	  cBinCount[i]->cd(j+1);
	  hJpsiMassUL[i][j]->SetMarkerStyle(21);
	  hJpsiMassUL[i][j]->SetMarkerColor(2);
	  hJpsiMassUL[i][j]->SetLineColor(2);
	  hJpsiMassUL[i][j]->SetTitle("");
	  hJpsiMassUL[i][j]->GetXaxis()->SetRangeUser(2.5,4);
	  hJpsiMassUL[i][j]->SetMaximum(1.5*hJpsiMassUL[i][j]->GetMaximum());
	  hJpsiMassUL[i][j]->Draw("P");
	  TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T}^{J/#Psi} < %1.0f GeV/c, #Deltatof < %1.2f ns",ptBins_low[i],ptBins_high[i],dtofCut[j]),0.05);
	  t1->Draw();
	  double binwidth = hJpsiMassUL[i][j]->GetBinWidth(1);
	  for(int itmp=0; itmp<3; itmp++)
	    {
	      counts[i][j][0][itmp] = 0;
	      counts[i][j][1][itmp] = 0;
	      TH1F *htmp = new TH1F(Form("%s_tmp%d",hJpsiMassUL[i][j]->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	      for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		{
		  int jbin = hJpsiMassUL[i][j]->FindFixBin(htmp->GetBinCenter(ibin));
		  htmp->SetBinContent(ibin,hJpsiMassUL[i][j]->GetBinContent(jbin));
		  counts[i][j][0][itmp] += hJpsiMassUL[i][j]->GetBinContent(jbin);
		  counts[i][j][1][itmp] += hJpsiMassLS[i][j]->GetBinContent(jbin);
		}
	      if(itmp==0) 
		{
		  htmp->SetLineColor(7);
		  htmp->SetFillColor(7);
		}
	      else
		{
		  htmp->SetLineColor(5);
		  htmp->SetFillColor(5);
		}
	      htmp->Draw("sames");
	    }
	  hJpsiMassLS[i][j]->Draw("samesHIST");
	  double sig_ul = counts[i][j][0][0];
	  double sig_ls = counts[i][j][1][0];
	  double bkg_ul = counts[i][j][0][1] + counts[i][j][0][2];
	  double bkg_ls = counts[i][j][1][1] + counts[i][j][1][2];
	  double scale = bkg_ul/bkg_ls;
	  double sig = sig_ul - sig_ls * scale;
	  double sig_err = sqrt(sig_ul+sig_ls*scale*scale);
	  t1 = GetPaveText(0.6,0.8,0.65,0.85,0.045);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("N = %2.0f#pm%2.0f",sig,sig_err));
	  t1->AddText(Form("S/B = 1:%2.2f",(sig_ul-sig)/sig));
	  t1->Draw();
	  hJpsiSigVsDtof[i][0]->SetBinContent(nScan-j, sig/sig_err);
	  hJpsiSigVsDtof[i][0]->SetBinError(nScan-j, 1e-4);
	  hJpsiErrVsDtof[i][0]->SetBinContent(nScan-j, sig_err/sig*100);
	  hJpsiErrVsDtof[i][0]->SetBinError(nScan-j, 1e-4);	  

	  hJpsiSigVsPt[j][0]->SetBinContent(i, sig/sig_err);
	  hJpsiSigVsPt[j][0]->SetBinError(i, 1e-4);
	}
      cBinCount[i]->cd(6);
      hJpsiSigVsDtof[i][0]->SetMaximum(1.5*hJpsiSigVsDtof[i][0]->GetMaximum());
      hJpsiSigVsDtof[i][0]->SetMinimum(1);
      hJpsiSigVsDtof[i][0]->SetMarkerStyle(20);
      hJpsiSigVsDtof[i][0]->SetMarkerSize(1.2);
      hJpsiSigVsDtof[i][0]->Draw("P");

      TF1 *func = new TF1(Form("Fit_%s",hJpsiSigVsDtof[i][0]->GetName()),"pol4",-0.125,1.125);
      hJpsiSigVsDtof[i][0]->Fit(func,"R0Q");
      func->SetLineColor(4);
      func->Draw("sames");
      hBestDtofCut[0]->SetBinContent(i, func->GetMaximumX(0,1));
      hBestDtofCut[0]->SetBinError(i, 1e-4);
      if(savePlot)  cBinCount[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiSig_DtofScan_BinCount_PtBin%d.pdf",run_type,i));
    }

  // fitting method
  TF1 *funcSignal = 0x0;
  TCanvas *cFit[nPtBins];
  for(int i=0; i<nPtBins; i++)
    {
      cFit[i] = new TCanvas(Form("JpsiMassFit_%d",i),Form("JpsiMassFit_pt%1.0f-%1.0f",ptBins_low[i],ptBins_high[i]),1100,700);
      cFit[i]->Divide(3,2);
      for(int j=0; j<nScan; j++)
	{
	  hJpsiMass[i][j] = (TH1F*)hJpsiMassUL[i][j]->Clone(Form("hJpsiMass_pt%1.0f-%1.0f_dtof%1.2f",ptBins_low[i],ptBins_high[i],dtofCut[j]));
	  double scale = (counts[i][j][0][1] + counts[i][j][0][2])/(counts[i][j][1][1] + counts[i][j][1][2]);
	  hJpsiMass[i][j]->Add(hJpsiMassLS[i][j], -1*scale);

	  cFit[i]->cd(j+1);
	  hJpsiMass[i][j]->SetMarkerColor(1);
	  hJpsiMass[i][j]->SetLineColor(1);
	  hJpsiMass[i][j]->SetMarkerStyle(20);
	  hJpsiMass[i][j]->GetXaxis()->SetRangeUser(2.5,4);
	  hJpsiMass[i][j]->SetTitle("");
	  hJpsiMass[i][j]->SetMaximum(1.5*hJpsiMass[i][j]->GetMaximum());
	  hJpsiMass[i][j]->Draw();
	  TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T}^{J/#Psi} < %1.0f GeV/c, #Deltatof < %1.2f ns",ptBins_low[i],ptBins_high[i],dtofCut[j]),0.05);
	  t1->Draw();

	  funcSignal = new TF1(Form("Fit_%s",hJpsiMass[i][j]->GetName()),Form("gausn(0)+pol3(3)"),2.5,4.0);
	  funcSignal->SetParameter(1,3.09);
	  funcSignal->SetParameter(2,0.01);
	  //if(i==1) funcSignal->FixParameter(2,0.05);
	  hJpsiMass[i][j]->Fit(funcSignal,"IR0QB");
	  funcSignal->SetLineColor(2);
	  funcSignal->Draw("sames");
	  double binwidth = hJpsiMass[i][j]->GetBinWidth(1);
	  t1 = GetPaveText(0.6,0.8,0.6,0.85,0.045);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("N = %2.0f#pm%2.0f",funcSignal->GetParameter(0)/binwidth,funcSignal->GetParError(0)/binwidth));
	  t1->AddText(Form("#mu = %2.2f#pm%2.2f",funcSignal->GetParameter(1),funcSignal->GetParError(1)));
	  t1->AddText(Form("#sigma = %2.2f#pm%2.2f",funcSignal->GetParameter(2),funcSignal->GetParError(2)));
	  t1->Draw();
	  hJpsiSigVsDtof[i][1]->SetBinContent(nScan-j, funcSignal->GetParameter(0)/funcSignal->GetParError(0));
	  hJpsiSigVsDtof[i][1]->SetBinError(nScan-j, 1e-4);
	  hJpsiErrVsDtof[i][1]->SetBinContent(nScan-j, funcSignal->GetParError(0)/funcSignal->GetParameter(0)*100);
	  hJpsiErrVsDtof[i][1]->SetBinError(nScan-j, 1e-4);

	  hJpsiSigVsPt[j][1]->SetBinContent(i, funcSignal->GetParameter(0)/funcSignal->GetParError(0));
	  hJpsiSigVsPt[j][1]->SetBinError(i, 1e-4);
	}
      cFit[i]->cd(6);
      hJpsiSigVsDtof[i][1]->SetMaximum(1.5*hJpsiSigVsDtof[i][1]->GetMaximum());
      hJpsiSigVsDtof[i][1]->SetMinimum(1);
      hJpsiSigVsDtof[i][1]->SetMarkerStyle(20);
      hJpsiSigVsDtof[i][1]->SetMarkerSize(1.2);
      hJpsiSigVsDtof[i][1]->Draw("P");
      TF1 *func = new TF1(Form("Fit_%s",hJpsiSigVsDtof[i][1]->GetName()),"pol4",-0.125,1.125);
      hJpsiSigVsDtof[i][1]->Fit(func,"R0Q");
      func->SetLineColor(4);
      func->Draw("sames");
      hBestDtofCut[1]->SetBinContent(i, func->GetMaximumX(0,1));
      hBestDtofCut[1]->SetBinError(i, 1e-4);
      if(savePlot)  cFit[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiSig_DtofScan_Fit_PtBin%d.pdf",run_type,i));
    }

  // compare different scan parameters
  TList *list = new TList();
  const char *method_title[2] = {"Bin-counting","Fitting"};
  for(int k=0; k<2; k++)
    {
      TLegend *leg1 = new TLegend(0.2,0.65,0.4,0.87);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(0);
      leg1->SetTextSize(0.035);
      TLegend *leg2 = new TLegend(0.5,0.65,0.7,0.87);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextSize(0.035);
      list->Clear();
      for(int i=0; i<nPtBins; i++)
	{
	  TH1F *htmp = (TH1F*)hJpsiSigVsDtof[i][k]->Clone(Form("%s_clone",hJpsiSigVsDtof[i][k]->GetName()));
	  htmp->SetMarkerStyle(20+i);
	  htmp->SetMarkerSize(1.3);
	  list->Add(htmp);
	  if(i<nPtBins/2) leg1->AddEntry(htmp,Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]),"P");
	  else            leg2->AddEntry(htmp,Form("%1.0f < p_{T} < %1.0f GeV/c",ptBins_low[i],ptBins_high[i]),"P");
	}
      c = drawHistos(list,Form("JpsiSigVsPt_%s",method[k]),Form("%s: significance of J/#psi signal (%s);#Deltatof (ns);Significance",run_type,method_title[k]),false,0,100,true,0,39,kFALSE,false,0x0,false,"",0.5,0.65,0.65,0.8,kTRUE,0.04,0.04,false,1,false);
      leg1->Draw();
      leg2->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiSigVsDtof_PtBins_%s.pdf",run_type,method[k]));

      list->Clear();
      for(int i=0; i<nPtBins; i++)
	{
	  hJpsiErrVsDtof[i][k]->SetMarkerStyle(20+i);
	  hJpsiErrVsDtof[i][k]->SetMarkerSize(1.3);
	  list->Add(hJpsiErrVsDtof[i][k]);
	}
      c = drawHistos(list,Form("JpsiErrVsPt_%s",method[k]),Form("%s: statistical error of J/#psi signal (%s);#Deltatof (ns);Stat. Err. (%%)",run_type,method_title[k]),false,0,100,true,0,39,kFALSE,false,0x0,false,"",0.5,0.65,0.65,0.8,kTRUE,0.04,0.04,false,1,false);
      leg1->Draw();
      leg2->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiErrVsDtof_PtBins_%s.pdf",run_type,method[k]));
    } 

  // list->Clear();
  // TString legName[2];
  // for(int k=0; k<2; k++)
  //   {
  //     legName[k] = method_title[k];
  //     list->Add(hBestDtofCut[k]);
  //   }
  // c = drawHistos(list,Form("JpsiSig_BestDtofCut"),Form("%s: #Deltatof cut with best significance;p_{T} (GeV/c);#Deltatof (ns)",run_type),false,0,100,true,0,1.5,kFALSE,true,legName,true,"",0.5,0.65,0.65,0.8,kTRUE);
}
