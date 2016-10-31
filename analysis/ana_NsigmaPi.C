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
      f = TFile::Open("./output/Pico.Run14.AuAu200.V0.root","read");
    }

  fix_NsigmaPi_V0();

  //fit_NsigmaPiDiff();
  //fix_NsigmaPi_m2();
}

//================================================
void fix_NsigmaPi_V0(const Int_t savePlot = 0)
{
  const int kNcent = 3, kNlumi = 4, kNpart = 2;
  const char* cent_name[kNcent] = {"0020","2040","4060"};
  const char* lumi_name[kNlumi] = {"all","low","mid","high"};
  const char* V0_name[kNpart] = {"K0s","Lambda"};

  const char* cent_title[kNcent] = {"0-20%","20-40%","40-60%"};
  const char* lumi_title[kNlumi] = {"all","prod_low","prod_mid","prod_high"}; 
  const char* V0_title[kNpart] = {"K^{0}_{s}","#Lambda"};

  TH2F *hV0MassVsPtUL[kNpart][kNcent][kNlumi];
  TH2F *hV0MassVsPtLS[kNpart][kNcent][kNlumi];
  
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=1; k<kNlumi; k++)
	    {
	      hV0MassVsPtUL[i][j][k] = (TH2F*)f->Get(Form("mh%sMassVsPtUL_di_mu_%s_%s",V0_name[i],cent_name[j],lumi_name[k]));
	      hV0MassVsPtLS[i][j][k] = (TH2F*)f->Get(Form("mh%sMassVsPtLS_di_mu_%s_%s",V0_name[i],cent_name[j],lumi_name[k]));
	      hV0MassVsPtUL[i][j][k]->Sumw2();
	      hV0MassVsPtLS[i][j][k]->Sumw2();
	      if(k==1)
		{
		  hV0MassVsPtUL[i][j][0] = (TH2F*)hV0MassVsPtUL[i][j][k]->Clone(Form("mh%sMassVsPtUL_di_mu_%s_%s",V0_name[i],cent_name[j],lumi_name[0]));
		  hV0MassVsPtLS[i][j][0] = (TH2F*)hV0MassVsPtLS[i][j][k]->Clone(Form("mh%sMassVsPtLS_di_mu_%s_%s",V0_name[i],cent_name[j],lumi_name[0]));
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
  const double min_mass[2] = {0.48, 1.110};
  const double max_mass[2] = {0.52, 1.123};
  const double min_fit[2] = {0.44, 1.11};
  const double max_fit[2] = {0.55, 1.15};
  const double v0_mass[2] = {0.4976, 1.12};
  TF1 *funcV0[kNpart][kNcent][kNlumi][nV0PtBin];
  TF1 *funcV0Bkg[kNpart][kNcent][kNlumi][nV0PtBin];

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
	      for(int p=0; p<nV0PtBin; p++)
		{
		  int low_bin  = hV0MassVsPtUL[i][j][k]->GetXaxis()->FindFixBin(V0PtBins[p]+1e-4);
		  int high_bin = hV0MassVsPtUL[i][j][k]->GetXaxis()->FindFixBin(V0PtBins[p+1]-1e-4);
		  TH1F *hUL = (TH1F*)hV0MassVsPtUL[i][j][k]->ProjectionY(Form("%s_bin%d",hV0MassVsPtUL[i][j][k]->GetName(),p),low_bin,high_bin);

		  if(i==1 && p==0) continue;

		  // fit V0 signal+background
		  funcV0[i][j][k][p] = new TF1(Form("func_%s",hUL->GetName()),"gaus(0)+pol3(3)",min_fit[i],max_fit[i]);
		  funcV0[i][j][k][p]->SetParameter(0,hUL->GetBinContent(hUL->FindFixBin(v0_mass[i])));
		  funcV0[i][j][k][p]->SetParameter(1, v0_mass[i]);
		  funcV0[i][j][k][p]->SetParameter(2, 0.001);
		  hUL->Fit(funcV0[i][j][k][p],"IRQ0");

		  funcV0Bkg[i][j][k][p] = new TF1(Form("funcBkg_%s",hUL->GetName()),"pol3",min_fit[i],max_fit[i]);
		  for(int ipar=0; ipar<4; ipar++)
		    {
		      funcV0Bkg[i][j][k][p]->SetParameter(ipar, funcV0[i][j][k][p]->GetParameter(ipar+3));
		      funcV0Bkg[i][j][k][p]->SetParError(ipar, funcV0[i][j][k][p]->GetParError(ipar+3));
		    }

		  hUL->SetLineColor(2);
		  hUL->SetMarkerColor(2);
		  hUL->SetMarkerStyle(20);
		  if(i==0) hUL->GetXaxis()->SetRangeUser(0.44,0.55);
		  if(i==1) hUL->GetXaxis()->SetRangeUser(1.08,1.15);
		  hUL->SetMaximum((1.5*hUL->GetMaximum()));
		  c1->cd(p+1);
		  hUL->Draw("P");					
		  funcV0[i][j][k][p]->SetLineColor(4);
		  funcV0[i][j][k][p]->SetLineStyle(2);
		  funcV0[i][j][k][p]->Draw("sames");
		  funcV0Bkg[i][j][k][p]->SetLineColor(4);
		  funcV0Bkg[i][j][k][p]->SetLineStyle(2);
		  funcV0Bkg[i][j][k][p]->Draw("sames");
		  TH1F *hLS = (TH1F*)hV0MassVsPtLS[i][j][k]->ProjectionY(Form("%s_bin%d",hV0MassVsPtLS[i][j][k]->GetName(),p),low_bin,high_bin);
		  hLS->Draw("samesHIST");
		  TPaveText *t1 = GetTitleText(Form("%1.1f < %s p_{T} < %1.1f GeV/c",V0PtBins[p],V0_title[i],V0PtBins[p+1]),0.05);
		  t1->Draw();
		  double ymax = hUL->GetMaximum();
		  TLine *line = GetLine(min_mass[i],0,min_mass[i],0.8*ymax,1);
		  line->Draw();
		  line = GetLine(max_mass[i],0,max_mass[i],0.8*ymax,1);
		  line->Draw();
		}
	      c1->cd(1);
	      TPaveText *t1 = GetPaveText(0.15,0.3,0.65,0.8,0.045);
	      t1->AddText(Form("%s, %s",cent_title[j],lumi_title[k]));
	      t1->AddText(Form("%s candidates",V0_title[i]));
	      t1->SetTextAlign(11);
	      t1->Draw();
	      if(savePlot)  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sMassInPtBins_%s_%s.pdf",run_type,V0_name[i],cent_name[j],lumi_name[k]));
	    }
	  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%sMassVsPt_%s.pdf",run_type,V0_name[i],cent_name[j]));
	}
    }
  return;

  // pure pion from K0s decay
  // pure pronton from Lambda decay
  const char* part_name[kNpart] = {"K0sPion","LambdaProton"};
  const char* part_title[kNpart] = {"pion","proton"};
  TH2F *hNsigmaPiVsPtUL[kNpart][kNcent][kNlumi];
  TH2F *hNsigmaPiVsPtLS[kNpart][kNcent][kNlumi];
  
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=1; k<kNlumi; k++)
	    {
	      hNsigmaPiVsPtUL[i][j][k] = (TH2F*)f->Get(Form("mh%sUL_di_mu_%s_%s",part_name[i],cent_name[j],lumi_name[k]));
	      hNsigmaPiVsPtLS[i][j][k] = (TH2F*)f->Get(Form("mh%sLS_di_mu_%s_%s",part_name[i],cent_name[j],lumi_name[k]));
	      hNsigmaPiVsPtUL[i][j][k]->Sumw2();
	      hNsigmaPiVsPtLS[i][j][k]->Sumw2();
	      if(k==1)
		{
		  hNsigmaPiVsPtUL[i][j][0] = (TH2F*)hNsigmaPiVsPtUL[i][j][k]->Clone(Form("mh%sUL_di_mu_%s_%s",part_name[i],cent_name[j],lumi_name[0]));
		  hNsigmaPiVsPtLS[i][j][0] = (TH2F*)hNsigmaPiVsPtLS[i][j][k]->Clone(Form("mh%sLS_di_mu_%s_%s",part_name[i],cent_name[j],lumi_name[0]));
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
  const Int_t nPartPtBin = 21;
  const Double_t partPtBins[nPartPtBin+1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
  for(int i=0; i<kNpart; i++)
    {
      for(int j=0; j<kNcent; j++)
	{
	  for(int k=0; k<kNlumi; k++)
	    {
	      hFitMean[i][j][k]  = new TH1F(Form("hNsigmaPiMean_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),";p_{T} (GeV/c);<n#sigma_{#pi}>",nPartPtBin,partPtBins);
	      hFitSigma[i][j][k] = new TH1F(Form("hNsigmaPiSigma_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),";p_{T} (GeV/c);#sigma(n#sigma_{#pi})",nPartPtBin,partPtBins);

	      TCanvas *c1 = new TCanvas(Form("FitNsigmaPi_%d_%d_%d",i,j,k),Form("FitNsigmaPi_%s_%s_%s",part_name[i],cent_name[j],lumi_name[k]),1200,800);
	      c1->Divide(4,6);
	      for(int p=0; p<nPartPtBin; p++)
		{
		  int low_bin  = hNsigmaPiVsPtUL[i][j][k]->GetXaxis()->FindFixBin(partPtBins[p]+1e-4);
		  int high_bin = hNsigmaPiVsPtUL[i][j][k]->GetXaxis()->FindFixBin(partPtBins[p+1]-1e-4);
		  TH1F *hUL = (TH1F*)hNsigmaPiVsPtUL[i][j][k]->ProjectionY(Form("%s_bin%d",hNsigmaPiVsPtUL[i][j][k]->GetName(),p),low_bin,high_bin);
		  hUL->SetLineColor(1);
		  hUL->SetMarkerColor(1);
		  hUL->SetMarkerStyle(24);
		  hUL->GetXaxis()->SetRangeUser(-5,5);
		  hUL->SetMaximum((1.5*hUL->GetMaximum()));
		  TF1 *func = new TF1(Form("Fit_%s",hUL->GetName()),"gaus",-2.5,2.5);
		  if(i==1) func->SetRange(-3,1);
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
void fit_NsigmaPiDiff(const Int_t save = 0)
{
  const Int_t func_color[4] = {2,8,4,6};
  gStyle->SetOptFit(1);
  // gStyle->SetStatY(0.98);                
  // gStyle->SetStatX(0.98);  
  // gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2); 

  const Int_t nTrkPtBin = 20;
  const Double_t trkPtBins[nTrkPtBin+1] = {1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,6.0,8.0,20.0};
  TH2F *hTrkDiffRebin[3];
  TH1F *hDiffRebin[3];

  const char *hName[] = {"mhTrkDiffpie","mhTrkDiffpiK","mhTrkDiffpiP"};
  TList *list = new TList;
  TH1F *hDiff[3];
  for(Int_t i=0; i<3; i++)
    {
      TH2 *h2 = (TH2*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      hDiff[i] = (TH1F*)h2->ProfileX(Form("%s_%s_pro",hName[i],trigName[kTrigType]));
      c = draw2D(h2,"");

      TString outname = hName[i];
      outname.ReplaceAll("mhTrk","NsigmaPi");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
      list->Add(hDiff[i]);

      hTrkDiffRebin[i] = new TH2F(Form("%s_%s_rebin",hName[i],trigName[kTrigType]),h2->GetTitle(),nTrkPtBin,trkPtBins,h2->GetNbinsY(),h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax());
      for(Int_t ibin=1; ibin<=hTrkDiffRebin[i]->GetNbinsX(); ibin++)
	{
	  for(Int_t jbin=1; jbin<=hTrkDiffRebin[i]->GetNbinsY(); jbin++)
	    {
	      Int_t low_x_bin = h2->GetXaxis()->FindBin(hTrkDiffRebin[i]->GetXaxis()->GetBinLowEdge(ibin)+0.01);
	      Int_t up_x_bin  = h2->GetXaxis()->FindBin(hTrkDiffRebin[i]->GetXaxis()->GetBinUpEdge(ibin)-0.01);
	      Double_t value = 0, error = 0;
	      for(Int_t iibin=low_x_bin; iibin<=up_x_bin; iibin++)
		{
		  value += h2->GetBinContent(iibin,jbin);
		  error += TMath::Power(h2->GetBinError(iibin,jbin),2);
		}
	      error = TMath::Sqrt(error);
	      hTrkDiffRebin[i]->SetBinContent(ibin,jbin,value);
	      hTrkDiffRebin[i]->SetBinError(ibin,jbin,error);
	    }
	}
      hDiffRebin[i] = (TH1F*)hTrkDiffRebin[i]->ProfileX(Form("%s_%s_pro_rebin",hName[i],trigName[kTrigType]));
    }

  TString legName[3] = {"n#sigma_{#pi}-n#sigma_{e}","n#sigma_{#pi}-n#sigma_{k}","n#sigma_{#pi}-n#sigma_{p}"}
  c = drawHistos(list,"NsigmaPiDiff",Form("%s: difference in <n#sigma_{#pi}> between different species;p_{T} (GeV/c);difference",trigName[kTrigType]),kFALSE,0,100,kTRUE,-3,8,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.65,0.6,0.85,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }


  list->Clear();
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hDiffRebin[i]);
    }
  c = drawHistos(list,"NsigmaPiDiff_Rebin",Form("%s: difference in <n#sigma_{#pi}> between different species;p_{T} (GeV/c);difference",trigName[kTrigType]),kFALSE,0,100,kTRUE,-3,8,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.65,0.6,0.85,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_rebin_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/NsigmaPiDiff_vs_pt_rebin_%s.png",run_type,trigName[kTrigType]));
    }

  TH2F *hNsigmaPiVsPt = (TH2F*)f->Get(Form("mhTrkNsigmaPi_%s",trigName[kTrigType]));
  draw2D(hNsigmaPiVsPt);

  TString partName[4] = {"pion","electron","kaon","proton"};
  TString partName2[3] = {"electron","kaon","proton"};

  TH1F *hMean[4];
  TH1F *hSigma[4];
  TH1F *hMeanDiff[3];
  TH1F *hParticleRatio[3];

  for(Int_t j=0; j<4; j++)
    {
      hMean[j]  = new TH1F(Form("FitNsigmaPiMean_%s",partName[j].Data()),Form("%s: fitted mean of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);<n#sigma_{#pi}>",trigName[kTrigType],partName[j].Data()),nTrkPtBin,trkPtBins);
      hSigma[j] = new TH1F(Form("FitNsigmaPiSigma_%s",partName[j].Data()),Form("%s: fitted sigma of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);#sigma(n#sigma_{#pi})",trigName[kTrigType],partName[j].Data()),nTrkPtBin,trkPtBins);

      if(j>0)
	{
	  hMeanDiff[j-1]  = new TH1F(Form("FitNsigmaPiMeanDiff_%s_%s",partName[0].Data(),partName[j].Data()),Form("%s: difference of fitted mean of n#sigma_{#pi} for %s and %s;p_{T} (GeV/c);#Delta<n#sigma_{#pi}>",trigName[kTrigType],partName[0].Data(),partName[j].Data()),nTrkPtBin,trkPtBins);
	  hParticleRatio[j-1]  = new TH1F(Form("FitParticleRatio_%s_%s",partName[0].Data(),partName[j].Data()),Form("%s: particle ratio of %s and %s;p_{T} (GeV/c);Ratio = %s/%s",trigName[kTrigType],partName[0].Data(),partName[j].Data()),nTrkPtBin,trkPtBins);
	}
    }

  TH1F *hNsigmaPi[nTrkPtBin];
  TF1 *func[nTrkPtBin][5];
  TCanvas *canvas[nTrkPtBin/4];
  for(Int_t i=0; i<nTrkPtBin; i++)
    {
      Int_t lowBin = hNsigmaPiVsPt->GetXaxis()->FindBin(trkPtBins[i]+0.01);
      Int_t hiBin  = hNsigmaPiVsPt->GetXaxis()->FindBin(trkPtBins[i+1]-0.01);
      
      printf("Track bin %d: %1.1f < pt < %1.1f GeV/c, %d < ibin < %d\n",i+1,trkPtBins[i],trkPtBins[i+1],lowBin,hiBin);
      hNsigmaPi[i] = (TH1F*)hNsigmaPiVsPt->ProjectionY(Form("hTrkNsigmaPi_%d",i),lowBin,hiBin);
      if(i%4==0)
      	{
      	  canvas[i/4] = new TCanvas(Form("TrkNsigmaPi_%d",i/4),Form("TrkNsigmaPi_%d",i/4),1100,750);
      	  canvas[i/4]->Divide(2,2);
      	}
      canvas[i/4]->cd(i%4+1);
      SetPadMargin(gPad,0.12,0.1,0.02,0.05);
      ScaleHistoTitle(hNsigmaPi[i],0.06,0.8,0.04,0.05,1.1,0.04,62);
      hNsigmaPi[i]->GetXaxis()->SetRangeUser(-6,10);
      hNsigmaPi[i]->SetMaximum(1.2*hNsigmaPi[i]->GetMaximum());
      hNsigmaPi[i]->Draw();

      //cout << hDiffRebin[1]->GetBinContent(i+1)-0.5 << "  " << hDiffRebin[1]->GetBinContent(i+1)+0.5 << endl;
      //func[i][0] = new TF1(Form("Fit_NsigmaPi_%d",i),"gaus(0)+gaus(3)+gaus(6)+gaus(9)",-10,10);
      func[i][0] = new TF1(Form("Fit_NsigmaPi_%d",i),"[0] * ( exp(-0.5*((x-[1])/[2])**2)+gaus(3)+gaus(6)+gaus(9) )",-10,10);

      func[i][0]->SetParNames("#pi yield","#pi mean","#pi #sigma",
			      "e/#pi","e mean","e #sigma",
			      "K/#pi","K mean","K #sigma");
      func[i][0]->SetParName(9,"P/#pi");
      func[i][0]->SetParName(10,"P mean");
      func[i][0]->SetParName(11,"P #sigma");

      func[i][0]->SetParLimits(0,0,1e10); // pion yield
      // func[i][0]->SetParLimits(3,0,1e10); // electron yield
      // func[i][0]->SetParLimits(6,0,1e10); // kaon yield
      // func[i][0]->SetParLimits(9,0,1e10); // proton yield

      func[i][0]->SetParLimits(3,0,1); // electron yield
      func[i][0]->SetParLimits(6,0,1); // kaon yield
      func[i][0]->SetParLimits(9,0,1); // proton yield

      func[i][0]->SetParameter(2,1); // pion sigma
      func[i][0]->SetParameter(5,1); // electron sigma
      func[i][0]->SetParameter(8,1); // kaon sigma
      func[i][0]->SetParameter(11,1); // proton sigma

      // func[i][0]->SetParameter(1,0); // pion mean
      // func[i][0]->SetParameter(4,hDiffRebin[0]->GetBinContent(i+1)); // electron mean
      // func[i][0]->SetParameter(7,hDiffRebin[1]->GetBinContent(i+1)); // kaon mean
      // func[i][0]->SetParameter(10,hDiffRebin[2]->GetBinContent(i+1)); // proton mean

      Double_t margin = 0.5;
      func[i][0]->SetParLimits(1,0-0.5,0+0.5); // pion mean
      func[i][0]->SetParLimits(4,hDiffRebin[0]->GetBinContent(i+1)-margin,hDiffRebin[0]->GetBinContent(i+1)+margin); // electron mean
      func[i][0]->SetParLimits(7,hDiffRebin[1]->GetBinContent(i+1)-margin,hDiffRebin[1]->GetBinContent(i+1)+margin); // kaon mean
      func[i][0]->SetParLimits(10,hDiffRebin[2]->GetBinContent(i+1)-margin,hDiffRebin[2]->GetBinContent(i+1)+margin); // proton mean
      hNsigmaPi[i]->Fit(func[i][0],"IRB");

      TPaveText *t1 = GetPaveText(0.25,0.33,0.85,0.92,0.05);
      t1->AddText(Form("%1.1f < p_{T} < %1.1f GeV/c",trkPtBins[i],trkPtBins[i+1]));
      t1->SetTextFont(62);
      t1->Draw();

      TLegend *leg = new TLegend(0.15,0.4,0.3,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      for(Int_t j=0; j<4; j++)
	{
	  func[i][j+1] = new TF1(Form("NsigmaPi_Species_%d",j),"gaus",-10,10);
	  if(j==0) func[i][j+1]->SetParameter(0,func[i][0]->GetParameter(3*j));
	  else     func[i][j+1]->SetParameter(0,func[i][0]->GetParameter(3*j)*func[i][0]->GetParameter(0)); 
	  func[i][j+1]->SetParameter(1,func[i][0]->GetParameter(3*j+1));
	  func[i][j+1]->SetParameter(2,func[i][0]->GetParameter(3*j+2));
	  func[i][j+1]->SetLineColor(func_color[j]);
	  func[i][j+1]->SetLineStyle(2);
	  func[i][j+1]->Draw("sames");
	  leg->AddEntry(func[i][j+1],partName[j].Data(),"L");

	  hMean[j]->SetBinContent(i+1,func[i][0]->GetParameter(3*j+1));
	  hMean[j]->SetBinError(i+1,func[i][0]->GetParError(3*j+1));

	  hSigma[j]->SetBinContent(i+1,TMath::Abs(func[i][0]->GetParameter(3*j+2)));
	  hSigma[j]->SetBinError(i+1,func[i][0]->GetParError(3*j+2));

	  if(j>0)
	    {
	      hMeanDiff[j-1]->SetBinContent(i+1,func[i][0]->GetParameter(3*j+1)-func[i][0]->GetParameter(1));
	      hMeanDiff[j-1]->SetBinError(i+1,TMath::Sqrt(TMath::Power(func[i][0]->GetParError(3*j+1),2) + TMath::Power(func[i][0]->GetParError(1),2)));

	      hParticleRatio[j-1]->SetBinContent(i+1,func[i][0]->GetParameter(3*j));
	      hParticleRatio[j-1]->SetBinError(i+1,func[i][0]->GetParError(3*j));
	    }
	}
      if(i%4==0)
	leg->Draw();
    }

  if(save) 
    {
      for(Int_t i=0; i<nTrkPtBin/4; i++)
	{
	  canvas[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitNsigmaPi_pt%d_%s.pdf",run_type,i,trigName[kTrigType]));
	  canvas[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitNsigmaPi_pt%d_%s.png",run_type,i,trigName[kTrigType]));
	}
    }

  // Different of mean
  list->Clear();
  for(Int_t i=0; i<3; i++)
    {
      list->Add(hDiffRebin[i]);
    }
  c = drawHistos(list,"NsigmaPiDiff_Rebin_compare",Form("%s: difference in <n#sigma_{#pi}> between different species;p_{T} (GeV/c);difference",trigName[kTrigType]),kFALSE,0,100,kTRUE,-3,8,kFALSE,kTRUE,legName,kTRUE,"",0.25,0.4,0.6,0.85,kTRUE);
  TLegend *leg = new TLegend(0.45,0.6,0.65,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  for(Int_t j=0; j<3; j++)
    {
      hMeanDiff[j]->SetMarkerStyle(25);
      hMeanDiff[j]->SetMarkerColor(color[j]);
      leg->AddEntry(hMeanDiff[j],Form("<n#sigma_{#pi}>_{%s} - <n#sigma_{#pi}>_{%s}",partName[j+1].Data(),partName[0].Data()));
      hMeanDiff[j]->Draw("samesP");
    }
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/Compare_NsigmaPiDiff_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/Compare_NsigmaPiDiff_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  // sigma
  list->Clear();
  for(Int_t j=0; j<4; j++)
    {
      TH1F *htmp = (TH1F*)hSigma[j]->Clone(Form("%s_clone",hSigma[j]->GetName()));
      list->Add(htmp);
    }
  c = drawHistos(list,"FitNSigmaPiSigma",Form("%s: fitted sigma of n#sigma_{#pi};p_{T} (GeV/c);#sigma(n#sigma_{#pi})",trigName[kTrigType]),kFALSE,0,100,kTRUE,0.1,4,kFALSE,kTRUE,partName,kTRUE,"",0.5,0.7,0.55,0.8,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitNsigmaPiWidth_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitNsigmaPiWidth_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  // e,K,P to pi ratio
  list->Clear();
  for(Int_t j=0; j<3; j++)
    {
      list->Add(hParticleRatio[j]);
    }
  c = drawHistos(list,"ParticleRatio",Form("%s: particle ratio to #pi from fitting;p_{T} (GeV/c);Ratio",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-4,10,kTRUE,kTRUE,partName2,kTRUE,"",0.5,0.7,0.18,0.35,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitParticleRatio_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitParticleRatio_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  hMean[0]->SetMarkerStyle(21);
  hMean[0]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c = draw1D(hMean[0],"");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitPionNsigmaPiMean_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitPionNsigmaPiMean_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  hSigma[0]->SetMarkerStyle(21);
  hSigma[0]->GetYaxis()->SetRangeUser(0.7,1.3);
  c = draw1D(hSigma[0],"");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitPionNsigmaPiSigma_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/FitPionNsigmaPiSigma_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }
}


//================================================
void fix_NsigmaPi_m2(const Int_t save = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.98);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.25); 

  const Int_t nTrkPtBin = 6;
  const Double_t trkPtBins[nTrkPtBin+1] = {1,1.5,2.0,2.5,3.0,5.0,20.0};
  const char *particleType[3] = {"pion","kaon","proton"};

  THnSparseF *hnTrkPid = (THnSparseF*)f->Get(Form("mhTrkPid_%s",trigName[kTrigType]));
  TH2F *hM2VsNsigmaPi = (TH2F*)hnTrkPid->Projection(2,1);
  hM2VsNsigmaPi->SetName(Form("hM2VsNsigmaPi_%s",trigName[kTrigType]));
  hM2VsNsigmaPi->SetTitle(Form("%s: m^{2} vs n#sigma_{#pi} of primary tracks;n#sigma_{#pi};m^{2}",trigName[kTrigType]));
  c = draw2D(hM2VsNsigmaPi);
  Double_t xmin = hM2VsNsigmaPi->GetXaxis()->GetXmin();
  Double_t xmax = hM2VsNsigmaPi->GetXaxis()->GetXmax();
  TLine *line = GetLine(xmin,pion_mass*pion_mass,xmax,pion_mass*pion_mass,1);
  line->Draw();
  line = GetLine(xmin,kaon_mass*kaon_mass,xmax,kaon_mass*kaon_mass,1);
  line->Draw();
  line = GetLine(xmin,proton_mass*proton_mass,xmax,proton_mass*proton_mass,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/M2_vs_NsigmaPi_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/M2_vs_NsigmaPi_%s.png",run_type,trigName[kTrigType]));
    }

  TH1F *hM2 = (TH1F*)hM2VsNsigmaPi->ProjectionY(Form("hM2_%s",trigName[kTrigType]));
  hM2->SetMarkerStyle(20);
  c = draw1D(hM2,Form("%s: m^{2} distribution of primary tracks;m^{2} (GeV/c^{2})^{2}",trigName[kTrigType]));
  
  TBox *box = new TBox(pion_m2_min,0,pion_m2_max,hM2->GetBinContent(hM2->FindBin(pion_mass*pion_mass)));
  box->SetFillStyle(1);
  box->SetFillColor(kGray);
  box->Draw();
  box = new TBox(kaon_m2_min,0,kaon_m2_max,hM2->GetBinContent(hM2->FindBin(kaon_mass*kaon_mass)));
  box->SetFillStyle(1);
  box->SetFillColor(kGray);
  box->Draw();
  box = new TBox(proton_m2_min,0,proton_m2_max,hM2->GetBinContent(hM2->FindBin(proton_mass*proton_mass)));
  box->SetFillStyle(1);
  box->SetFillColor(kGray);
  box->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/M2_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/M2_%s.png",run_type,trigName[kTrigType]));
    }

  TH1F *hNsigmaPi[3][nTrkPtBin];
  TF1 *funcNsigmaPi[3][nTrkPtBin];
  TH1F *hMean[3];
  TH1F *hSigma[3];
  for(Int_t i=0; i<3; i++)
    {
      hMean[i] = new TH1F(Form("hNsigmaPi_mean_%s",particleType[i]),Form("Mean of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);<n#sigma_{#pi}>",particleType[i]),nTrkPtBin,trkPtBins);
      hSigma[i] = new TH1F(Form("hNsigmaPi_sigma_%s",particleType[i]),Form("Sigma of n#sigma_{#pi} for %s candidates;p_{T} (GeV/c);<n#sigma_{#pi}>",particleType[i]),nTrkPtBin,trkPtBins);
      
      if(i==0) hnTrkPid->GetAxis(2)->SetRangeUser(pion_m2_min,pion_m2_max);
      else if (i==1) hnTrkPid->GetAxis(2)->SetRangeUser(kaon_m2_min,kaon_m2_max);
      else if (i==2) hnTrkPid->GetAxis(2)->SetRangeUser(proton_m2_min,proton_m2_max);
      for(Int_t j=0; j<nTrkPtBin; j++)
	{
	  hnTrkPid->GetAxis(0)->SetRangeUser(trkPtBins[j]+0.01,trkPtBins[j+1]-0.01);
	  hNsigmaPi[i][j] = (TH1F*)hnTrkPid->Projection(1);
	  hNsigmaPi[i][j]->Sumw2();
	  hNsigmaPi[i][j]->SetName(Form("hNsigmaPi_%s_pt_%1.1f_%1.1f",particleType[i],trkPtBins[j],trkPtBins[j+1]));
	  hnTrkPid->GetAxis(0)->SetRange(0,-1);
	}
      hnTrkPid->GetAxis(2)->SetRange(0,-1);


      TCanvas *c = new TCanvas(Form("cNsigmaPi_%s",particleType[i]),Form("cNsigmaPi_%s",particleType[i]),1200,600);
      c->Divide(3,2);
      for(Int_t j=0; j<nTrkPtBin; j++)
	{
	  c->cd(j+1);
	  SetPadMargin(gPad,0.12,0.1,0.02,0.02);
	  hNsigmaPi[i][j]->SetTitle("");
	  ScaleHistoTitle(hNsigmaPi[i][j],0.06,0.8,0.04,0.05,1.1,0.04,62);

	  funcNsigmaPi[i][j] = new TF1(Form("funcNsigmaPi_%s_pt_%1.1f_%1.1f",particleType[i],trkPtBins[j],trkPtBins[j+1]),"gaus",-5,5);
	  if(i==2 && j==0) 
	    {
	      funcNsigmaPi[i][j]->SetRange(-2.5,5.5);
	      hNsigmaPi[i][j]->GetXaxis()->SetRangeUser(-5,10);
	    }
	  hNsigmaPi[i][j]->Fit(funcNsigmaPi[i][j],"IR0");
	  hNsigmaPi[i][j]->SetMarkerStyle(21);
	  hNsigmaPi[i][j]->SetMaximum(1.3*hNsigmaPi[i][j]->GetMaximum());
	  hNsigmaPi[i][j]->Draw("P");
	  funcNsigmaPi[i][j]->SetLineColor(2);
	  funcNsigmaPi[i][j]->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.3,0.35,0.86,0.93,0.06);
	  t1->AddText(Form("%1.1f < p_{T} < %1.1f GeV/c",trkPtBins[j],trkPtBins[j+1]));
	  t1->SetTextFont(62);
	  t1->Draw();

	  hMean[i]->SetBinContent(j+1,funcNsigmaPi[i][j]->GetParameter(1));
	  hMean[i]->SetBinError(j+1,funcNsigmaPi[i][j]->GetParError(1));

	  hSigma[i]->SetBinContent(j+1,funcNsigmaPi[i][j]->GetParameter(2));
	  hSigma[i]->SetBinError(j+1,funcNsigmaPi[i][j]->GetParError(2));
	}
      c->cd(1);
      t1 = GetPaveText(0.2,0.25,0.6,0.65,0.07);
      t1->AddText(Form("%s",particleType[i]));
      t1->SetTextFont(62);
      t1->SetTextColor(4);
      t1->Draw();

      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/Fit_NsigmaPi_%s_%s.pdf",run_type,particleType[i],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_NsigmaPi/Fit_NsigmaPi_%s_%s.png",run_type,particleType[i],trigName[kTrigType]));
	}
    }

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


  const char *hName[] = {"mhTrkDiffpiK","mhTrkDiffpiP"};
  TH1F *hDiff[2];
  for(Int_t i=0; i<2; i++)
    {
      TH2 *h2 = (TH2*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      hDiff[i] = (TH1F*)h2->ProfileX(Form("%s_%s_pro",hName[i],trigName[kTrigType]));
      c = draw2D(h2,"");
    }

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
