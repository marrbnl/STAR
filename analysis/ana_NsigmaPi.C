TFile *f;
const int year = 2014;

//================================================
void ana_NsigmaPi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if(year==2013)
    {
      //f = TFile::Open("./output/Run13.pp500.jpsi.CutRanking.root","read");
    }
  else if(year==2014)
    {
    }

  //NsigmaPi_V0();
  NsigmaPi_m2();
  //jpsiMuon();
  //NsigmaPiDiff();

}

//================================================
void NsigmaPi_V0(const Int_t savePlot = 0)
{
  const int kNcent = 3, kNlumi = 4, kNpart = 3;
  const char* cent_name[kNcent] = {"0020","2040","4060"};
  const char* lumi_name[kNlumi] = {"all","low","mid","high"};
  const char* V0_name[kNpart] = {"K0s","Lambda","Photon"};

  const char* cent_title[kNcent] = {"0-20%","20-40%","40-60%"};
  const char* lumi_title[kNlumi] = {"all","prod_low","prod_mid","prod_high"}; 
  const char* V0_title[kNpart] = {"K^{0}_{s}","#Lambda","photon"};

  f = TFile::Open("./output/Pico.Run14.AuAu200.V0.root","read");
  TH2F *hV0MassVsPtUL[kNpart][kNcent][kNlumi];
  TH2F *hV0MassVsPtLS[kNpart][kNcent][kNlumi];
  
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=1; k<kNlumi; k++)
	    {
	      hV0MassVsPtUL[i][j][k] = (TH2F*)f->Get(Form("mhInvMassVsPtUL_%s_%s_%s_di_mu",V0_name[i],cent_name[j],lumi_name[k]));
	      hV0MassVsPtLS[i][j][k] = (TH2F*)f->Get(Form("mhInvMassVsPtLS_%s_%s_%s_di_mu",V0_name[i],cent_name[j],lumi_name[k]));
	      hV0MassVsPtUL[i][j][k]->Sumw2();
	      hV0MassVsPtLS[i][j][k]->Sumw2();
	      if(i==0)
		{
		  hV0MassVsPtUL[i][j][k]->RebinY(5);
		  hV0MassVsPtLS[i][j][k]->RebinY(5);
		}
	      else if(i==2)
		{
		  hV0MassVsPtUL[i][j][k]->RebinY(5);
		  hV0MassVsPtLS[i][j][k]->RebinY(5);
		}
	      if(k==1)
		{
		  hV0MassVsPtUL[i][j][0] = (TH2F*)hV0MassVsPtUL[i][j][k]->Clone(Form("mhInvMassVsPtUL_%s_%s_%s_di_mu",V0_name[i],cent_name[j],lumi_name[0]));
		  hV0MassVsPtLS[i][j][0] = (TH2F*)hV0MassVsPtLS[i][j][k]->Clone(Form("mhInvMassVsPtLS_%s_%s_%s_di_mu",V0_name[i],cent_name[j],lumi_name[0]));
		}
	      else
		{
		  hV0MassVsPtUL[i][j][0]->Add(hV0MassVsPtUL[i][j][k]);
		  hV0MassVsPtLS[i][j][0]->Add(hV0MassVsPtLS[i][j][k]);
		}
	    }
	}
    }

  const Int_t nV0PtBin = 6;
  const Double_t V0PtBins[nV0PtBin+1] = {0,1,2,3,4,5,10};
  const double min_mass[kNpart] = {0.48, 1.110, 0};
  const double max_mass[kNpart] = {0.52, 1.123, 0.1};
  const double v0_mass[kNpart] = {0.4976, 1.12, 0};

  // side-band for Kshort
  const double min_side[2] = {0.40, 0.56};
  const double max_side[2] = {0.44, 0.60};
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  TCanvas *c = new TCanvas(Form("MassVsPt_%d_%d",i,j),Form("MassVsPt_%s_%s",V0_name[i],cent_name[j]),1200,700);
	  c->Divide(2,2);
	  for(int k=0; k<kNlumi; k++)
	    {
	      c->cd(k+1);
	      gPad->SetLogz();
	      hV0MassVsPtUL[i][j][k]->SetTitle();
	      hV0MassVsPtUL[i][j][k]->Draw("colz");
	      TPaveText *t1 = GetTitleText(Form("Invariant mass vs. p_{T} for %s candidates (%s, %s)",V0_title[i],cent_title[j],lumi_title[k]),0.05);
	      t1->Draw();

	      TCanvas *c1 = new TCanvas(Form("Mass_%d_%d_%d",i,j,k),Form("MassVsPt_%s_%s_%s",V0_name[i],cent_name[j],lumi_name[k]),1200,700);
	      c1->Divide(3,2);
	      TLegend *leg = new TLegend(0.15,0.6,0.4,0.88);
	      leg->SetBorderSize(0);
	      leg->SetFillColor(0);
	      leg->SetTextSize(0.04);
	      leg->SetHeader(Form("%s, %s",cent_title[j],lumi_title[k]));
	      for(int p=0; p<nV0PtBin; p++)
		{
		  int low_bin  = hV0MassVsPtUL[i][j][k]->GetXaxis()->FindFixBin(V0PtBins[p]+1e-4);
		  int high_bin = hV0MassVsPtUL[i][j][k]->GetXaxis()->FindFixBin(V0PtBins[p+1]-1e-4);
		  TH1F *hUL = (TH1F*)hV0MassVsPtUL[i][j][k]->ProjectionY(Form("%s_bin%d",hV0MassVsPtUL[i][j][k]->GetName(),p),low_bin,high_bin);
		  hUL->SetMarkerStyle(20);
		  if(i==0) hUL->GetXaxis()->SetRangeUser(0.38,0.62);
		  if(i==1) hUL->GetXaxis()->SetRangeUser(1.08,1.15);
		  if(i==2) hUL->GetXaxis()->SetRangeUser(0,0.2);
		  hUL->SetMaximum((1.5*hUL->GetMaximum()));
		  c1->cd(p+1);
		  hUL->Draw("P");
		  if(p==0) leg->AddEntry(hUL,"Unlike-sign","P");

		  // draw signal region
		  double binwidth = hUL->GetBinWidth(1);
		  TH1F *htmp = new TH1F(Form("%s_signal",hUL->GetName()),"",int((max_mass[i]-min_mass[i])/binwidth+0.5),min_mass[i],max_mass[i]);
		  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		    {
		      int jbin = hUL->FindFixBin(htmp->GetBinCenter(ibin));
		      htmp->SetBinContent(ibin,hUL->GetBinContent(jbin));
		    }
		  htmp->SetFillColor(7);
		  htmp->Draw("sames");
		  if(p==0) leg->AddEntry(htmp,"Signal","F");


		  // side-band only for Kshort
		  if(i==0)
		    {
		      for(int itmp=0; itmp<2; itmp++)
			{
			  TH1F *htmp = new TH1F(Form("%s_tmp_%d",hUL->GetName(),itmp),"",int((max_side[itmp]-min_side[itmp])/binwidth),min_side[itmp],max_side[itmp]);
			  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
			    {
			      int jbin = hUL->FindFixBin(htmp->GetBinCenter(ibin));
			      htmp->SetBinContent(ibin,hUL->GetBinContent(jbin));
			    }
			  htmp->SetFillColor(5);
			  htmp->Draw("sames");
			  if(p==0 && itmp==0) leg->AddEntry(htmp,"Side-band","F");
			}
		    }

		  TH1F *hLS = (TH1F*)hV0MassVsPtLS[i][j][k]->ProjectionY(Form("%s_bin%d",hV0MassVsPtLS[i][j][k]->GetName(),p),low_bin,high_bin);
		  hLS->SetLineColor(2);
		  hLS->Draw("samesHIST");
		  if(p==0 && i==0) leg->AddEntry(hLS,"Rotational","L");
		  if(p==0 && i>=1) leg->AddEntry(hLS,"Like-sign","L");
		  TPaveText *t1 = GetTitleText(Form("%1.1f < %s p_{T} < %1.1f GeV/c",V0PtBins[p],V0_title[i],V0PtBins[p+1]),0.05);
		  t1->Draw();
		}
	      c1->cd(1);
	      leg->Draw();
	      if(savePlot)  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sMassInPtBins_%s_%s.pdf",run_type,V0_name[i],cent_name[j],lumi_name[k]));
	    }
	  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sMassVsPt_%s.pdf",run_type,V0_name[i],cent_name[j]));
	}
    }

  // pure pion from K0s decay (side-band)
  // pure pronton from Lambda decay (like-sign distribution)
  // pure electron from gamma conversion
  const char* part_name[kNpart] = {"Pion","Proton","Electron"};
  const char* part_title[kNpart] = {"pion","proton","electron"};

  TH2F *hNsigmaPiVsPtUL[kNpart][kNcent][kNlumi];
  TH2F *hNsigmaPiVsPtLS[kNpart][kNcent][kNlumi];
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=1; k<kNlumi; k++)
	    {
	      hNsigmaPiVsPtUL[i][j][k] = (TH2F*)f->Get(Form("mhNsigmaPiVsPtUL_%s_%s_%s_di_mu",part_name[i],cent_name[j],lumi_name[k]));
	      hNsigmaPiVsPtLS[i][j][k] = (TH2F*)f->Get(Form("mhNsigmaPiVsPtLS_%s_%s_%s_di_mu",part_name[i],cent_name[j],lumi_name[k]));
	      hNsigmaPiVsPtUL[i][j][k]->Sumw2();
	      hNsigmaPiVsPtLS[i][j][k]->Sumw2();
	      if(k==1)
		{
		  hNsigmaPiVsPtUL[i][j][0] = (TH2F*)hNsigmaPiVsPtUL[i][j][k]->Clone(Form("mhNsigmaPiVsPtUL_%s_%s_%s_di_mu",part_name[i],cent_name[j],lumi_name[0]));
		  hNsigmaPiVsPtLS[i][j][0] = (TH2F*)hNsigmaPiVsPtLS[i][j][k]->Clone(Form("mhNsigmaPiVsPtLS_%s_%s_%s_di_mu",part_name[i],cent_name[j],lumi_name[0]));
		}
	      else
		{
		  hNsigmaPiVsPtUL[i][j][0]->Add(hNsigmaPiVsPtUL[i][j][k]);
		  hNsigmaPiVsPtLS[i][j][0]->Add(hNsigmaPiVsPtLS[i][j][k]);
		}
	    }
	  for(int k=0; k<kNlumi; k++)
	    {
	      hNsigmaPiVsPtUL[i][j][k]->Add(hNsigmaPiVsPtLS[i][j][k],-1);
	    }
	}
    }

  TH1F *hFitMean[kNpart][kNcent][kNlumi];
  TH1F *hFitSigma[kNpart][kNcent][kNlumi];
  const Int_t nPartPtBin = 11;
  const Double_t partPtBins[nPartPtBin+1] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=0; k<kNlumi; k++)
	    {
	      hFitMean[i][j][k]  = new TH1F(Form("hNsigmaPiMean_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),";p_{T} (GeV/c);<n#sigma_{#pi}>",nPartPtBin,partPtBins);
	      hFitSigma[i][j][k] = new TH1F(Form("hNsigmaPiSigma_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),";p_{T} (GeV/c);#sigma(n#sigma_{#pi})",nPartPtBin,partPtBins);

	      TCanvas *c1 = new TCanvas(Form("FitNsigmaPi_%d_%d_%d",i,j,k),Form("FitNsigmaPi_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),1200,800);
	      c1->Divide(3,4);
	      for(int p=0; p<nPartPtBin; p++)
		{
		  int low_bin  = hNsigmaPiVsPtUL[i][j][k]->GetXaxis()->FindFixBin(partPtBins[p]+1e-4);
		  int high_bin = hNsigmaPiVsPtUL[i][j][k]->GetXaxis()->FindFixBin(partPtBins[p+1]-1e-4);
		  TH1F *hUL = (TH1F*)hNsigmaPiVsPtUL[i][j][k]->ProjectionY(Form("%s_bin%d",hNsigmaPiVsPtUL[i][j][k]->GetName(),p),low_bin,high_bin);
		  hUL->SetLineColor(1);
		  hUL->SetMarkerColor(1);
		  hUL->SetMarkerStyle(24);
		  hUL->Rebin(2);
		  if(i<2)  hUL->GetXaxis()->SetRangeUser(-5,5);
		  if(i==2) hUL->GetXaxis()->SetRangeUser(0,10);
		  hUL->SetMaximum((1.5*hUL->GetMaximum()));
		  TF1 *func = new TF1(Form("Fit_%s",hUL->GetName()),"gaus",-3,3);
		  if(i==1) func->SetRange(-4,4);
		  if(i==2) func->SetRange(2,7);
		  hUL->Fit(func,"R0Q");
		  c1->cd(p+2);
		  hUL->SetTitle("");
		  hUL->Draw("P");
		  func->SetLineColor(2);
		  func->Draw("sames");
		  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",partPtBins[p],partPtBins[p+1]),0.09);
		  t1->Draw();
		  hFitMean[i][j][k]->SetBinContent(p+1, func->GetParameter(1));
		  hFitMean[i][j][k]->SetBinError(p+1, func->GetParError(1));
		  hFitSigma[i][j][k]->SetBinContent(p+1, func->GetParameter(2));
		  hFitSigma[i][j][k]->SetBinError(p+1, func->GetParError(2));
		}
	      c1->cd(1);
	      TPaveText *t1 = GetPaveText(0.3,0.6,0.3,0.8,0.12);
	      t1->AddText(run_type);
	      t1->AddText(Form("%s, %s",cent_title[j],lumi_title[k]));
	      t1->AddText(Form("%s candidates",part_title[i]));
	      t1->SetTextAlign(11);
	      t1->Draw();
	      if(savePlot)  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitNsigmaPi_%s_%s_%s.pdf",run_type,part_title[i],cent_name[j],lumi_name[k]));
	    }
	}
    }

  for(int i=0; i<kNpart; i++)
    {
      TCanvas *c1 = new TCanvas(Form("NsigmaPiMean_%d",i),Form("CompNsigmaPiMean_%s",part_name[i]),1100,700);
      c1->Divide(2,2);
      TLegend *leg11 = new TLegend(0.6,0.6,0.8,0.85);
      leg11->SetBorderSize(0);
      leg11->SetFillColor(0);
      leg11->SetTextSize(0.05);
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=0; k<kNlumi; k++)
	    {
	      c1->cd(j+1);
	      hFitMean[i][j][k]->SetMarkerStyle(20+j);
	      hFitMean[i][j][k]->SetMarkerColor(color[k]);
	      hFitMean[i][j][k]->SetLineColor(color[k]);
	      if(i==0) hFitMean[i][j][k]->GetYaxis()->SetRangeUser(-2,2);
	      if(i==1) hFitMean[i][j][k]->GetYaxis()->SetRangeUser(-4,2);
	      if(i==2) hFitMean[i][j][k]->GetYaxis()->SetRangeUser(0,8);
	      SetPadMargin(gPad,0.12,0.1,0.1,0.1);
	      ScaleHistoTitle(hFitMean[i][j][k],0.05,1.0,0.04,0.05,0.9,0.04,62);
	      if(k==0) 
		{
		  hFitMean[i][j][k]->Draw();
		  TPaveText *t1 = GetTitleText(Form("Mean of n#sigma_{#pi} for %s from %s decay (%s)",part_title[i],V0_title[i],cent_title[j]),0.05);
		  t1->Draw();
		}
	      else     
		{
		  hFitMean[i][j][k]->Draw("sames");
		}
	      if(j==0) leg11->AddEntry(hFitMean[i][j][k],lumi_title[k],"P");
	    }
	}
      c1->cd(1);
      leg11->Draw();
      c1->cd(4);
      SetPadMargin(gPad,0.12,0.1,0.1,0.1);
      TLegend *leg12 = new TLegend(0.6,0.6,0.8,0.85);
      leg12->SetBorderSize(0);
      leg12->SetFillColor(0);
      leg12->SetTextSize(0.05);
      leg12->SetHeader(run_type);
      for(int j=0; j<kNcent; j++)
	{
	  if(j==0) hFitMean[i][j][0]->Draw();
	  else     hFitMean[i][j][0]->Draw("sames");
	  leg12->AddEntry(hFitMean[i][j][0],cent_title[j],"P");
	}
      TPaveText *t1 = GetTitleText(Form("Mean of n#sigma_{#pi} for %s from %s decay",part_title[i],V0_title[i]),0.05);
      t1->Draw();
      leg12->Draw();
      if(savePlot)  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sNSigmaPiMean.pdf",run_type,part_name[i]));


      // width vs. pt
      TCanvas *c2 = new TCanvas(Form("NsigmaPiWidth_%d",i),Form("CompNsigmaPiWidth_%s",part_name[i]),1100,700);
      c2->Divide(2,2);
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=0; k<kNlumi; k++)
	    {
	      c2->cd(j+1);
	      hFitSigma[i][j][k]->SetMarkerStyle(20+j);
	      hFitSigma[i][j][k]->SetMarkerColor(color[k]);
	      hFitSigma[i][j][k]->SetLineColor(color[k]);
	      hFitSigma[i][j][k]->GetYaxis()->SetRangeUser(0,2);
	      SetPadMargin(gPad,0.12,0.1,0.1,0.1);
	      ScaleHistoTitle(hFitSigma[i][j][k],0.05,1.0,0.04,0.05,0.9,0.04,62);
	      if(k==0) 
		{
		  hFitSigma[i][j][k]->Draw();
		  TPaveText *t1 = GetTitleText(Form("Width of n#sigma_{#pi} for %s from %s decay (%s)",part_title[i],V0_title[i],cent_title[j]),0.05);
		  t1->Draw();
		}
	      else     
		{
		  hFitSigma[i][j][k]->Draw("sames");
		}
	    }
	}
      c2->cd(1);
      leg11->Draw();
      c2->cd(4);
      SetPadMargin(gPad,0.12,0.1,0.1,0.1);
      for(int j=0; j<kNcent; j++)
	{
	  if(j==0) hFitSigma[i][j][0]->Draw();
	  else     hFitSigma[i][j][0]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("Width of n#sigma_{#pi} for %s from %s decay",part_title[i],V0_title[i]),0.05);
      t1->Draw();
      leg12->Draw();
      if(savePlot)  c2->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sNSigmaPiWidth.pdf",run_type,part_name[i]));
    }
}

//================================================
void jpsiMuon(const Int_t savePlot = 0, const Int_t saveHisto = 0)
{
}


//================================================
void NsigmaPiDiff(const Int_t savePlot = 0, const int saveHisto = 0)
{
  f = TFile::Open("output/Pico.Run14.AuAu200.V0.root","read");
  const int nHistos = 2;
  const char* part_name[nHistos] = {"Proton","Electron"};
  const char* part_title[nHistos] = {"proton","electron"};
  const char *hName[nHistos] = {"P","e"};
  const char *track[2] = {"primary","global"};
  TH2F *hDiffVsPt[2][nHistos];
  TH1F *hTheoryMean[2][nHistos];
  TH1F *hTheorySigma[2][nHistos];
  for(Int_t i=0; i<nHistos; i++)
    {
      for(int j=0; j<2; j++)
	{
	  if(j==0) hDiffVsPt[j][i] = (TH2F*)f->Get(Form("mhPrimTrkDiffpi%s_di_mu",hName[i]));
	  if(j==1) hDiffVsPt[j][i] = (TH2F*)f->Get(Form("mhTrkDiffpi%s_di_mu",hName[i]));
	  hDiffVsPt[j][i]->SetName(Form("NsigmaPiDiffVsPt_%s_%s",part_name[i],track[j]));

	  if(i==0) hDiffVsPt[j][i]->GetYaxis()->SetRangeUser(-4,10);
	  if(i==1) hDiffVsPt[j][i]->GetYaxis()->SetRangeUser(0,10);
	  c = draw2D(hDiffVsPt[j][i],Form("%s: n#sigma_{#pi}-n#sigma_{%s} distribution for %s tracks",run_type,hName[i],track[j]));

	  TObjArray fitslice;
	  hDiffVsPt[j][i]->FitSlicesY(0, 0, -1, 0, "QNR", &fitslice);      

	  hTheoryMean[j][i]  =  (TH1F*)((TH1F*)fitslice[1])->Clone(Form("NsigmaPiDiffMean_%s_%s",part_name[i],track[j]));
	  hTheoryMean[j][i]->SetMarkerStyle(21);
	  hTheoryMean[j][i]->Draw("sames");

	  hTheorySigma[j][i]  =  (TH1F*)((TH1F*)fitslice[2])->Clone(Form("NsigmaPiDiffSigma_%s_%s",part_name[i],track[j]));
	  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NSigmaPiDiff_%s_%s.pdf",run_type,part_name[i],track[j]));
	}
    }

  TList *list = new TList;
  TString legName[2] = {"Primary tracks","Global tracks"};
  for(Int_t i=0; i<nHistos; i++)
    {
      list->Clear();
      for(int j=0; j<2; j++)
	{
	  TH1F *hClone = (TH1F*)hTheoryMean[j][i]->Clone(Form("%s_clone",hTheoryMean[j][i]->GetName()));
	  list->Add(hClone);
	}
      c = drawHistos(list,Form("Compare_NsigmaPiDiffMean_%s",part_name[i]),Form("%s: mean of n#sigma_{#pi}-n#sigma_{%s} distributions;p_{T} (GeV/c);mean",run_type,hName[i]),kFALSE,0,100,false,-3,2,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.65,0.65,0.8,kTRUE);
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NSigmaPiDiff_%s_CompMean.pdf",run_type,part_name[i]));

      list->Clear();
      for(int j=0; j<2; j++)
	{
	  TH1F *hClone = (TH1F*)hTheorySigma[j][i]->Clone(Form("%s_clone",hTheorySigma[j][i]->GetName()));
	  list->Add(hClone);
	}
      c = drawHistos(list,Form("Compare_NsigmaPiDiffSigma_%s",part_name[i]),Form("%s: width of n#sigma_{#pi}-n#sigma_{%s} distributions;p_{T} (GeV/c);width",run_type,hName[i]),kFALSE,0,100,false,-3,2,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.65,0.65,0.8,kTRUE);
      if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NSigmaPiDiff_%s_CompSigma.pdf",run_type,part_name[i]));
    }

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile* fout = TFile::Open(Form("Rootfiles/%s.NsigmaPi.root",run_type),"update");
      for(Int_t i=0; i<nHistos; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      hDiffVsPt[j][i]->Write("",TObject::kOverwrite);
	      hTheoryMean[j][i]->Write("",TObject::kOverwrite);
	      hTheorySigma[j][i]->Write("",TObject::kOverwrite);
	    }
	}
    }
}


//================================================
void NsigmaPi_m2(const Int_t savePlot = 0, const int saveHisto = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.98);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.25); 

  const Int_t nPartPtBin = 10;
  const Double_t partPtBins[nPartPtBin+1] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
  const int kNpart = 3;
  const char* part_name[kNpart] = {"Pion","Proton","Electron"};
  const char* part_title[kNpart] = {"pion","proton","electron"};
  
  f = TFile::Open("output/Pico.Run14.AuAu200.V0.root","read");


  //==============================================
  // m2 vs pt distribution
  //==============================================
  TH2F *hPrimM2VsPt = (TH2F*)f->Get("mhPrimM2VsPt_di_mu");
  c = draw2D(hPrimM2VsPt, Form("%s: m^{2} vs p_{T} of primary tracks;p_{T} (GeV/c);m^{2}",run_type));
  Double_t xmin = hPrimM2VsPt->GetXaxis()->GetXmin();
  Double_t xmax = hPrimM2VsPt->GetXaxis()->GetXmax();
  TLine *line = GetLine(xmin,pion_mass*pion_mass,xmax,pion_mass*pion_mass,1);
  line->Draw();
  line = GetLine(xmin,kaon_mass*kaon_mass,xmax,kaon_mass*kaon_mass,1);
  line->Draw();
  line = GetLine(xmin,proton_mass*proton_mass,xmax,proton_mass*proton_mass,1);
  line->Draw();
  const double m2_cuts[4] = {0.016,0.021,0.85,0.89};
  for(int i=0; i<4; i++)
    {
      line = GetLine(xmin,m2_cuts[i],xmax,m2_cuts[i],2);
      line->Draw();
    }
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/m2Prim_vs_pt.pdf",run_type));
    }
  return;

  //==============================================
  // NsigmaPi distribution 
  //==============================================
  TH2F *hNsigmaPiVsPt[kNpart];
  TH1F *hNsigmaPi[kNpart][nPartPtBin];
  TF1 *funcNsigmaPi[kNpart][nPartPtBin];
  TH1F *hMean[kNpart];
  TH1F *hSigma[kNpart];
  for(Int_t i=0; i<3; i++)
    {
      hMean[i] = new TH1F(Form("m2Prim_%s_NsigmaPi_FitMean",part_name[i]),Form("Mean of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);<n#sigma_{#pi}>",part_title[i]),nPartPtBin,partPtBins);
      hSigma[i] = new TH1F(Form("m2Prim_%s_NsigmaPi_FitSigma",part_name[i]),Form("Width of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);#sigma(n#sigma_{#pi})",part_title[i]),nPartPtBin,partPtBins);
      hNsigmaPiVsPt[i] = (TH2F*)f->Get(Form("mhPrimNsigmaPiVsPt_%s_di_mu",part_name[i]));
      TCanvas *c = new TCanvas(Form("cNsigmaPi_%s",part_title[i]),Form("cNsigmaPi_%s",part_title[i]),1200,700);
      c->Divide(4,3);
      for(Int_t j=0; j<nPartPtBin; j++)
	{
	  int low_bin  = hNsigmaPiVsPt[i]->GetXaxis()->FindFixBin(partPtBins[j]+1e-4);
	  int high_bin = hNsigmaPiVsPt[i]->GetXaxis()->FindFixBin(partPtBins[j+1]-1e-4);
	  hNsigmaPi[i][j] = (TH1F*)hNsigmaPiVsPt[i]->ProjectionY(Form("m2Prim_%s_NsigmaPi_bin%d",part_name[i],j+1),low_bin,high_bin);
	  //hNsigmaPi[i][j]->Rebin(2);
	  hNsigmaPi[i][j]->Sumw2();

	  double ymax = hNsigmaPi[i][j]->GetMaximum();
	  if(i==0)
	    { 
	      funcNsigmaPi[i][j] = new TF1(Form("m2Prim_%s_NsigmaPiFunc_bin%d",part_name[i],j+1),"gaus(0)+gaus(3)",-4,4);
	      funcNsigmaPi[i][j]->SetParameters(ymax,-1,1,ymax,1,1);
	      hNsigmaPi[i][j]->GetXaxis()->SetRangeUser(-4,4);
	    }
	  else if(i==1) 
	    {
	      funcNsigmaPi[i][j] = new TF1(Form("m2Prim_%s_NsigmaPiFunc_bin%d",part_name[i],j+1),"gaus(0)",-6,6);
	      funcNsigmaPi[i][j]->SetParameters(ymax,1,1);
	      hNsigmaPi[i][j]->GetXaxis()->SetRangeUser(-6,6);
	    }
	  else if(i==2)
	    {
	      funcNsigmaPi[i][j] = new TF1(Form("m2Prim_%s_NsigmaPiFunc_bin%d",part_name[i],j+1),"gaus(0)+gaus(3)",0,6);
	      funcNsigmaPi[i][j]->SetParameters(ymax,3,0.5,ymax,5,0.5);
	      hNsigmaPi[i][j]->GetXaxis()->SetRangeUser(0,7);
	    }
	  hNsigmaPi[i][j]->Fit(funcNsigmaPi[i][j],"IR0");

	  c->cd(j+1);
	  SetPadMargin(gPad,0.12,0.1,0.02,0.02);
	  hNsigmaPi[i][j]->SetTitle("");
	  ScaleHistoTitle(hNsigmaPi[i][j],0.06,0.8,0.04,0.05,1.1,0.04,62);
	  hNsigmaPi[i][j]->SetMarkerStyle(21);
	  hNsigmaPi[i][j]->SetMaximum(1.3*hNsigmaPi[i][j]->GetMaximum());
	  hNsigmaPi[i][j]->Draw("P");
	  funcNsigmaPi[i][j]->SetLineColor(2);
	  funcNsigmaPi[i][j]->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.3,0.35,0.86,0.93,0.06);
	  t1->AddText(Form("%1.1f < p_{T} < %1.1f GeV/c",partPtBins[j],partPtBins[j+1]));
	  t1->SetTextFont(62);
	  t1->Draw();

	  hMean[i]->SetBinContent(j+1,funcNsigmaPi[i][j]->GetParameter(1));
	  hMean[i]->SetBinError(j+1,funcNsigmaPi[i][j]->GetParError(1));

	  hSigma[i]->SetBinContent(j+1,funcNsigmaPi[i][j]->GetParameter(2));
	  hSigma[i]->SetBinError(j+1,funcNsigmaPi[i][j]->GetParError(2));
	}
      c->cd(1);
      t1 = GetPaveText(0.2,0.25,0.6,0.65,0.07);
      t1->AddText(Form("%s",part_name[i]));
      t1->SetTextFont(62);
      t1->SetTextColor(4);
      t1->Draw();

      if(savePlot) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/m2Prim_%s_NsigmaPiFit.pdf",run_type,part_name[i]));
	}
    }

  /*
  TString legName[3] = {"Pion sample","Kaon sample", "Proton sample"};

  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hMean[i]);
    }
  c = drawHistos(list,"Mean_NsigmaPi",Form("%s: mean of n#sigma_{#pi} vs p_{T};p_{T} (GeV/c);<n#sigma_{#pi}>",trigName[kTrigType]),kFALSE,0,100,kTRUE,-2,2,kFALSE,kTRUE,legName,kTRUE,"",0.6,0.75,0.6,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPi_mean_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPi_mean_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  list->Clear();
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hSigma[i]);
    }
  c = drawHistos(list,"Sigma_NsigmaPi",Form("%s: sigma of n#sigma_{#pi} vs p_{T};p_{T} (GeV/c);#sigma(n#sigma_{#pi})",trigName[kTrigType]),kFALSE,0,100,kTRUE,0.5,2.5,kFALSE,kTRUE,legName,kTRUE,"",0.2,0.35,0.65,0.85,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPi_sigma_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPi_sigma_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }
  */

  const char *hName[] = {"mhPrimTrkDiffpiP","mhPrimTrkDiffpie"};
  TH1F *hDiff[2];
  for(Int_t i=0; i<2; i++)
    {
      TH2 *h2 = (TH2*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      hDiff[i] = (TH1F*)h2->ProfileX(Form("%s_%s_pro",hName[i],trigName[kTrigType]));
      c = draw2D(h2,"");
    }
  return;
  list->Clear();
  for(Int_t i=0; i<2; i++)
    {
      TH1F *hClone = (TH1F*)hMean[i]->Clone(Form("%s_clone",hMean[i]->GetName()));
      list->Add(hClone);
    }
  list->Add(hDiff[1]);
  list->Add(hDiff[2]);
  TString legName2[4] = {"n#sigma_{#pi} of kaon sample","n#sigma_{#pi} of proton sample","n#sigma_{#pi}-n#sigma_{k} of tracks","n#sigma_{#pi}-n#sigma_{p} of tracks"}
  c = drawHistos(list,"Compare_NsigmaPiDiff",Form("%s: difference in <n#sigma_{#pi}> between different species;p_{T} (GeV/c);difference",trigName[kTrigType]),kFALSE,0,100,kTRUE,-3,2,kFALSE,kTRUE,legName2,kTRUE,"",0.5,0.65,0.6,0.85,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }
}
