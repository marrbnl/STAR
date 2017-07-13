const int year = YEAR;

//================================================
void ana_MtdRespEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  //ana_cosmicRay();
  embedVsCosmic();
  //systematics();
}

//================================================
void systematics(const int savePlot = 1, const int saveHisto = 1)
{
  // const char* part_name = "Ups1S";
  // const char* part_title = "Y(1S)";
  const char* part_name = "Jpsi";
  const char* part_title = "J/psi";

  const int year = 2014;
  if(year==2013)
    {
      const int nPtBins = 6;
      double xPtBins[7] = {0,1,2,3,4,6,8};
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      const int nPtBins = 9;
      double xPtBins[10] = {0,1,2,3,4,5,6,8,10,15};
      run_type = "Run14_AuAu200";
    }
  else if(year==2015)
    {
      const int nPtBins = 8;
      double xPtBins[9] = {0,1,2,3,4,5,6,8,10};
      run_type = "Run15_pp200";
    }
  else if(year==2016)
    {
      const int nPtBins = 6;
      double xPtBins[7] = {0,1,2,3,4,6,10};
      run_type = "Run14_AuAu200";
    }

  // open the corresponding files
  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdRespEff.root",run_type),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdRespEff.root",run_type),"read");

  TH1F *hSysMtdRespEff[2][3];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<3; j++)
	{
	  if(i==0)
	    {
	      hSysMtdRespEff[i][j] = new TH1F(Form("%s_hMtdRespEffSys%d",part_name,j+1),Form("%s: uncertainty of MTD response efficiency for %s;p_{T} (GeV/c)",run_type,part_title),nPtBins,xPtBins);
	    }
	  else
	    {
	      hSysMtdRespEff[i][j] = new TH1F(Form("%s_Npart_hMtdRespEffSys%d",part_name,j+1),Form("%s: uncertainty of MTD response efficiency for %s;p_{T} (GeV/c)",run_type,part_title),2,0,2);
	      hSysMtdRespEff[i][j]->GetXaxis()->SetBinLabel(1, "p_{T} > 0 GeV/c");
	      hSysMtdRespEff[i][j]->GetXaxis()->SetBinLabel(2, "p_{T} > 5 GeV/c");
	    }
	  if(j==2)
	    {
	      if(i==0) hSysMtdRespEff[i][j]->SetName(Form("%s_hMtdRespEffSysAll",part_name));
	      if(i==1) hSysMtdRespEff[i][j]->SetName(Form("%s_Npart_hMtdRespEffSysAll",part_name));
	    }
	}
    }

  TH1F *hJpsiEffPtBin[nPtBins+2][2];
  for(int i=0; i<nPtBins+2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  if(i<nPtBins) hJpsiEffPtBin[i][j] = (TH1F*)fdata->Get(Form("hJpsiEffSpreadSys%d_Pt%1.0f-%1.0f",j+1,xPtBins[i],xPtBins[i+1]));
	  else          hJpsiEffPtBin[i][j] = (TH1F*)fdata->Get(Form("hJpsiEffSpreadSys%d_Pt%d-15",j+1,(i-nPtBins)*5));
	}
    }

  // Uncertainty I: fit errors
  TF1 *funcJpsiEffPtBin[nPtBins+2];
  TCanvas *cFit = new TCanvas("cFit","cFit",1100,700);
  cFit->Divide(4,3);
  for(int i=0; i<nPtBins+2; i++)
    {
      double init_mean = hJpsiEffPtBin[i][0]->GetMean();
      double init_rms  = hJpsiEffPtBin[i][0]->GetRMS();
      funcJpsiEffPtBin[i] = new TF1(Form("funcJpsiEffPtBin%d",i+1),"gaus",init_mean-3*init_rms,init_mean+3*init_rms);
      if(year==2014)
	{
	  if(i<2) hJpsiEffPtBin[i][0]->Rebin(4);
	  else if(i!=nPtBins-1) hJpsiEffPtBin[i][0]->Rebin(2);
	}
      hJpsiEffPtBin[i][0]->Fit(funcJpsiEffPtBin[i],"0IR");
      double mean = funcJpsiEffPtBin[i]->GetParameter(1);
      double sigma = funcJpsiEffPtBin[i]->GetParameter(2);
      cFit->cd(i+1);
      hJpsiEffPtBin[i][0]->SetMarkerStyle(25);
      hJpsiEffPtBin[i][0]->GetXaxis()->SetRangeUser(mean-5*sigma,mean+5*sigma);
      hJpsiEffPtBin[i][0]->SetMaximum(1.8*hJpsiEffPtBin[i][0]->GetMaximum());
      hJpsiEffPtBin[i][0]->Draw("P");
      funcJpsiEffPtBin[i]->SetLineColor(2);
      funcJpsiEffPtBin[i]->Draw("sames");
      if(i<nPtBins)
	{
	  hSysMtdRespEff[0][0]->SetBinContent(i+1,1);
	  hSysMtdRespEff[0][0]->SetBinError(i+1,sigma/mean);
	  TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T} < %1.0f GeV/c",xPtBins[i],xPtBins[i+1]),0.07);
	  t1->Draw();
	}
      else
	{
	  hSysMtdRespEff[1][0]->SetBinContent(i-nPtBins+1, 1);
	  hSysMtdRespEff[1][0]->SetBinError(i-nPtBins+1, sigma/mean);
	  TPaveText *t1 = GetTitleText(Form("p_{T} > %d GeV/c",(i-nPtBins)*5),0.07);
	  t1->Draw();
	}
    }
  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/ToyMC_%s_FitSpread_Sys1.pdf",run_type,part_name));

  // Uncertainty II: template
  for(int i=0; i<nPtBins+2; i++)
    {
      double eff1 = 0, eff2 = 0;
      int nbins = hJpsiEffPtBin[i][1]->GetNbinsX();
      for(int bin=1; bin<=nbins; bin++)
	{
	  double eff = hJpsiEffPtBin[i][1]->GetBinCenter(bin);
	  if(hJpsiEffPtBin[i][1]->GetBinContent(bin)>0)
	    {
	      cout << eff << endl;
	      if(eff1==0) eff1 = eff;
	      else        eff2 = eff;
	    }
	}
      if(i<nPtBins)
	{
	  hSysMtdRespEff[0][1]->SetBinContent(i+1,1);
	  hSysMtdRespEff[0][1]->SetBinError(i+1,fabs(1-eff1/eff2));
	}
      else
	{
	  hSysMtdRespEff[1][1]->SetBinContent(i-nPtBins+1, 1);
	  hSysMtdRespEff[1][1]->SetBinError(i-nPtBins+1, fabs(1-eff1/eff2));
	}
    }

  // Combine the uncertainties
  const char *legName[3] = {"Fit error","Template","All"};
  for(int i=0; i<2; i++)
    {
      for(int bin=1; bin<=hSysMtdRespEff[i][2]->GetNbinsX(); bin++)
	{
	  double e1 = hSysMtdRespEff[i][0]->GetBinError(bin);
	  double e2 = hSysMtdRespEff[i][1]->GetBinError(bin);
	  double error = sqrt(e1*e1 + e2*e2);
	  hSysMtdRespEff[i][2]->SetBinContent(bin, 1);
	  hSysMtdRespEff[i][2]->SetBinError(bin, error);
	}
      hSysMtdRespEff[i][2]->SetMarkerStyle(21);
      if(year==2013) hSysMtdRespEff[i][2]->GetYaxis()->SetRangeUser(0.85,1.15);
      if(year==2014) hSysMtdRespEff[i][2]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(year==2015) hSysMtdRespEff[i][2]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(i==1)
	{
	  hSysMtdRespEff[i][2]->SetXTitle("");
	  hSysMtdRespEff[i][2]->GetXaxis()->SetLabelSize(0.06);
	}
      c = draw1D(hSysMtdRespEff[i][2],Form("%s: %s uncertainty due to MTD response efficiency",run_type,part_title));
      TLegend *leg = new TLegend(0.6,0.7,0.8,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.035);
      leg->AddEntry(hSysMtdRespEff[i][2],legName[2],"P");
      for(int j=0; j<2; j++)
	{
	  hSysMtdRespEff[i][j]->SetMarkerStyle(24);
	  hSysMtdRespEff[i][j]->SetMarkerColor(4+j*2);
	  hSysMtdRespEff[i][j]->SetLineColor(4+j*2);
	  TGraphErrors *gr = new TGraphErrors(hSysMtdRespEff[i][j]);
	  offset_x(gr,0.15+j*0.15);
	  gr->Draw("samesPEZ");
	  leg->AddEntry(hSysMtdRespEff[i][j],legName[j],"P");
	}
      leg->Draw();
      if(savePlot)
	{
	  if(i==0) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Jpsi_MtdRespEffSysAll.pdf",run_type));
	  if(i==1) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Npart.Jpsi_MtdRespEffSysAll.pdf",run_type));
	}
    }

  if(saveHisto)
    {
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<3; j++)
	    {
	      hSysMtdRespEff[i][j]->Write("",TObject::kOverwrite);
	      hSysMtdRespEff[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void embedVsCosmic(const int savePlot = 1, const int saveHisto = 1)
{
  // fit response efficiency in embedding
  TFile *femb = TFile::Open(Form("output/%s.Embed.Jpsi.root",run_type),"read");
  TH2F *hProjTrack = (TH2F*)femb->Get("mhProjTrack");
  hProjTrack->Sumw2();
  TH2F *hMthTrack = (TH2F*)femb->Get("mhMatchTrack");
  hMthTrack->Sumw2();
  c = draw2D(hProjTrack);
  c = draw2D(hMthTrack);
  TH1F *hAll[30][5];
  TH1F *hAcc[30][5];
  TH1F *hEff[30][5];
  TF1  *funcEff[30][5];
  TCanvas *cEff[6];
  for(int i=0; i<6; i++)
    {
      cEff[i] = new TCanvas(Form("EmbedMtdRespEff_BL%d-%d",i*5+1,i*5+5),Form("EmbedMtdRespEff_BL%d-%d",i*5+1,i*5+5),1100,800);
      cEff[i]->Divide(5,5);
    }
  const int nPtBins = 21;
  const double xPtBins[nPtBins+1] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10};
  TH1F *htmp = 0x0;
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  funcEff[i][j] = 0x0;
	  int bin = i*5+j+1;
	  htmp = (TH1F*)hProjTrack->ProjectionX(Form("hProj_%d_%d",i+1,j+1),bin,bin);
	  hAll[i][j] = (TH1F*)htmp->Rebin(nPtBins, Form("hProj_BL%d_Mod%d",i+1,j+1), xPtBins);
	  htmp = (TH1F*)hMthTrack->ProjectionX(Form("hMth_%d_%d",i+1,j+1),bin,bin);
	  hAcc[i][j] = (TH1F*)htmp->Rebin(nPtBins, Form("hMth_BL%d_Mod%d",i+1,j+1),xPtBins);
	  hEff[i][j] = (TH1F*)hAcc[i][j]->Clone(Form("Embed_MtdRespEff_BL%d_Mod%d",i+1,j+1));
	  hEff[i][j]->Divide(hAll[i][j]);
	  cEff[i/5]->cd(i%5*5+j+1);
	  SetPadMargin(gPad,0.13,0.13,0.02,0.1);
	  ScaleHistoTitle(hEff[i][j], 0.065, 0.9, 0.05, 0.06, 1.2, 0.05);
	  hEff[i][j]->SetMarkerStyle(20);
	  funcEff[i][j] = new TF1(Form("Embed_funcMtdRespEff_BL%d_Mod%d",i+1,j+1),"[0]/([1]*x+[2]*x*x-[3])+[4]",1.1,10);
	  if(htmp->GetEntries()<1) continue;
	  funcEff[i][j]->SetParameters(-0.1,1,1,1.1,0.98);
	  hEff[i][j]->Fit(funcEff[i][j],"IR0Q");
	  hEff[i][j]->Draw();
	  funcEff[i][j]->SetLineColor(4);
	  funcEff[i][j]->SetLineStyle(2);
	  funcEff[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("BL = %d, Mod = %d",i+1,j+1),0.07);
	  t1->Draw();
	}
    }
  if(savePlot)
    {
      for(int i=0; i<6; i++)
	{
	  cEff[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Embed.MtdRespEff_BL%d-%d.pdf",run_type, i*5+1, i*5+5));
	}
    }

  TH1F *hplot = new TH1F("hplot",Form("%s: MTD response efficiency from embedding;p_{T} (GeV/c);Eff",run_type),100,0,10);
  c = draw1D(hplot);
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  if(funcEff[i][j]->GetParameter(0)!=0)
	    funcEff[i][j]->Draw("sames");
	}
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Embed.MtdRespEff_vs_BL.pdf",run_type));

  // embed vs. cosmic ray
  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");

  TF1 *funcData[30][5];
  TF1 *funcRatio[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  funcData[i][j] = 0x0;
	  funcRatio[i][j] = 0x0;
	  if(!funcEff[i][j]) continue;
	  funcData[i][j] = (TF1*)fdata->Get(Form("funcMtdRespEffvsPt_Bkl%d_Mod%d",i+1,j+1));
	  funcRatio[i][j] = new TF1(Form("MtdRespEff_Bkl%d_Mod%d",i+1,j+1),"([0]/(x-[1])+[2])/([3]/([4]*x+[5]*x*x-[6])+[7])",1.2,10);
	  for(int par=0; par<8; par++)
	    {
	      if(par<3) funcRatio[i][j]->SetParameter(par, funcData[i][j]->GetParameter(par));
	      else      funcRatio[i][j]->SetParameter(par, funcEff[i][j]->GetParameter(par-3));
	    }
	}
    }
  hplot->SetName("hRatio");
  hplot->GetYaxis()->SetRangeUser(0.1,1.2);
  hplot->GetXaxis()->SetRangeUser(1, 10);
  hplot->SetYTitle("Cosmic/Embedding");
  c = draw1D(hplot,"Cosmic/Embedding: MTD response efficiency");
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  if(funcEff[i][j]->GetParameter(0)!=0)
	    {
	      funcRatio[i][j]->SetLineStyle(2);
	      funcRatio[i][j]->Draw("sames");
	    }
	}
    }
  TLine *line = GetLine(1.3, 0.1, 1.3, 1.1, 2, 2, 1);
  line->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Embed.MtdRespEff_to_cosmic.pdf",run_type));
  if(saveHisto)
    {
      fdata->cd();
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      hEff[i][j]->Write("",TObject::kOverwrite);
	      funcEff[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void ana_cosmicRay(const int savePlot = 0)
{
  const int nBL = 13;
  TCanvas *cBL[nBL];
  TH1F *hCosmic[nBL][5];
  TF1 *hFitCosmic[nBL][5];
  TF1 *hFitTemplate[nBL][5];
  TH1F *hFitError[nBL][5];
  TF1 *hRatio[nBL][5];

  TFile *fRespEff = TFile::Open(Form("Rootfiles/Run%dResponseEffViaPtTemplate.root",year-2000),"read");
  TH1F *hplot = new TH1F("hplot","",100,0,10);
  for(int i=0; i<nBL; i++)
    {
      int nMod = 5;
      int bl = i+10;
      if(bl>11 && bl<21) 
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1100,700);
	  cBL[i]->Divide(3,2);
	  nMod = 3;
	}
      else
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1500,600);
	  cBL[i]->Divide(5,2);
	}

      
      for(int j=0; j<5; j++)
	{
	  hRatio[i][j] = new TF1(Form("hRatio_%d_%d",bl,j+1),"([0]/(x-[1])+[2])/(([3]/(x-[4])+[5])*[6])");
	  hFitCosmic[i][j] = (TF1*)fRespEff->Get(Form("fPtMtdEffBkl%d_Mod%d",bl-1,j));
	  hFitTemplate[i][j] = (TF1*)fRespEff->Get(Form("fSclPtTmpBkl%d_Mod%d",bl-1,j));
	  hFitError[i][j] = (TH1F*)fRespEff->Get(Form("MtdRespEffCosmic_BL%d_Mod%d",bl,j+1));
	  for(int bin=1; bin<=hFitError[i][j]->GetNbinsX(); bin++)
	    {
	      hFitError[i][j]->SetBinError(bin, 2*hFitError[i][j]->GetBinError(bin));
	    }
	  hCosmic[i][j] = (TH1F*)fRespEff->Get(Form("hPtMtdEffBkl%d_Mod%d",bl-1,j));
	  int index = j;
	  if(bl>11 && bl<21)
	    {
	      if(j==0 || j==4) continue;
	      if(bl==15 && j==3) continue;
	      index = j-1;
	    }
	  
	  cBL[i]->cd(index+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  ScaleHistoTitle(hFitError[i][j],0.05,1,0.04,0.05,1.1,0.04,62);
	  hFitError[i][j]->SetTitle(";p_{T} (GeV/c);MTD response efficiency");
	  hFitError[i][j]->GetXaxis()->SetRangeUser(1.3,4);
	  hFitError[i][j]->GetYaxis()->SetRangeUser(0.2, 1.1);
	  if(bl==19) hFitError[i][j]->GetYaxis()->SetRangeUser(0, 0.7);
	  hFitError[i][j]->SetFillStyle(3001);
	  hFitError[i][j]->SetFillColor(kBlue);
	  hFitError[i][j]->SetLineStyle(2);
	  hFitError[i][j]->Draw("e5");
	  hCosmic[i][j]->Draw("sames");
	  //hFitCosmic[i][j]->Draw("sames");
	  hFitTemplate[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("Bkl=%d, Mod=%d",bl,j+1),0.06);
	  t1->Draw();
	  TLegend *leg = new TLegend(0.45,0.2,0.7,0.42);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->AddEntry(hCosmic[i][j],"Cosmic data","P");
	  leg->AddEntry(hFitCosmic[i][j],"Fit to cosmic","L");
	  leg->AddEntry(hFitTemplate[i][j],"Template","L");
	  leg->Draw();
	  for(int ipar=0; ipar<7; ipar++)
	    {
	      if(ipar<3) hRatio[i][j]->SetParameter(ipar, hFitCosmic[i][j]->GetParameter(ipar));
	      else       hRatio[i][j]->SetParameter(ipar, hFitTemplate[i][j]->GetParameter(ipar-3));
	    }
	  cBL[i]->cd(index+nMod+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  hplot->SetTitle(";p_{T} (GeV/c);Cosmic/Template");
	  ScaleHistoTitle(hplot,0.05,1,0.04,0.05,1,0.04,62);
	  hplot->GetXaxis()->SetRangeUser(1.3,4);
	  hplot->GetYaxis()->SetRangeUser(0.8,1.2);
	  hplot->Draw("HIST");
	  hRatio[i][j]->SetRange(1.3,4);
	  hRatio[i][j]->Draw("sames");
	  TLine *line = new TLine(1.3,1,4,1);
	  line->SetLineColor(2);
	  line->SetLineStyle(2);
	  line->Draw();
	}
      if(savePlot) cBL[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/BL%d_TemplateVsCosmic.pdf",run_type,bl));
    }

}
