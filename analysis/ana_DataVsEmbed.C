const int year = YEAR;
TF1 *fResVsPt[3];
const int nHistos = 6;
const char *name[nHistos] = {"Dz","Dy","NsigmaPi","Dca","NHitsFit","NHitsDedx"};
const char *title[nHistos] = {"#Deltaz","#Deltay","n#sigma_{#pi}","DCA","NHitsFit","NHitsDedx"};
const char *unit[nHistos] = {" (cm)"," (cm)",""," (cm)","",""};
const double cuts[nHistos] = {0, 0, 0, 1, 15, 10};
const Int_t nMuonPtBin = 6;
const Double_t muonPtBins[nMuonPtBin+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0};
char* data_name[2] = {"Data","Embed"};
double ptbound = 100;
double nsigma1 = 2, nsigma2 = 2;

//================================================
void ana_DataVsEmbed()
{
  gStyle->SetOptStat(0);

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
  fResVsPt[0]->SetParameters(-21.04, 21.09, 0.693);
  fResVsPt[1] = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[1]->SetParameters(-12.61, 13.43, 0.889);
  fResVsPt[2] = new TF1("fResDtofVsPt","[0]+[1]*exp([2]/x)",1,10);
  fResVsPt[2]->SetParameters(0.0817528, 0.0169419, 4.34897);

  //makeData();
  //makeEmbed();
  dataVsEmbed();
  //DyDzCut();
  //Run16nSigmaPiEff();
}


//================================================
void Run16nSigmaPiEff(const bool savePlot = 0, const bool saveHisto = 0)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  const int nbins = 6;
  const double xbins[nbins+1] = {1.3, 1.5, 2.0, 2.5, 3.0, 5, 10};
  TF1* funcSigma = new TF1("func_MuonSigma","pol0",1.3,10);
  funcSigma->SetParameter(0, 0.9985);
  funcSigma->SetParError(0, 0.02662);
  TF1* funcMean[2];
  funcMean[0] = new TF1("func_MuonMean","[0]+[1]*exp([2]/x)",1.3,10);
  funcMean[0]->SetParameters(-4.841, 5.1, 0.07498);
  double errors[3] = {41.72, 41.65, 0.5859};
  funcMean[0]->SetParErrors(errors);
  funcMean[1] = new TF1("func_MuonMean_sys","pol0",1.3,10);
  funcMean[1]->SetParameter(0, 0.4309);
  funcMean[1]->SetParError(0, 0.03729);
  TH1F *hEff[2];
  TF1 *funcEff[2];
  for(int i=0; i<1; i++)
    {
      hEff[i] = new TH1F(Form("hEff_%d",i),Form("hEff_%d",i),nbins,xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
	  double pt = hEff[i]->GetBinCenter(bin);
	  TF1* functmp = new TF1(Form("functmp_%d_%d",i,bin),"gaus",-5,5);
	  functmp->SetParameters(1, funcMean[i]->Eval(pt), funcSigma->Eval(pt));
	  double eff = functmp->Integral(-1,3)/functmp->Integral(-5,5);
	  hEff[i]->SetBinContent(bin, eff);
	  hEff[i]->SetBinError(bin, 3e-3);
	}
      funcEff[i] = new TF1(Form("funcEff_%d",i),"[0]+[1]*exp([2]/x)",1.3,10);
      funcEff[i]->SetParameters(0.4588, 0.4374, 0.1118);
      hEff[i]->Fit(funcEff[i], "R0Q");
      hEff[i]->GetYaxis()->SetRangeUser(0.8,1.0);
      hEff[i]->SetMarkerStyle(21);
    }
  TCanvas *c = draw1D(hEff[0],"Run16_AuAu200: efficiency of -1 < n#sigma_{#pi} <3 cut for single muons;p_{T} (GeV/c)");
  funcEff[0]->SetLineColor(2);
  funcEff[0]->SetLineStyle(2);
  funcEff[0]->Draw("sames");
  TLegend *leg = new TLegend(0.45,0.25,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hEff[0],"Efficiency via fitting","P");
  leg->AddEntry(funcEff[0],"Fit to efficiency","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/Run16_AuAu200/ana_DataVsEmbe/DataJpsiMuon_nSigmaPiEff.pdf"));
  
}

//================================================
void DyDzCut(const bool savePlot = 0, const bool saveHisto = 0, const bool saveAN = 1)
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.15);

  TFile *fin = 0x0;
  if(saveHisto) fin = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
  else          fin = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"read");

  TH1F *hFitSigma[2][2];
  TF1 *func[2];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hFitSigma[j][i] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_FitSigma",data_name[j],name[i]));
	  hFitSigma[j][i]->SetMarkerStyle(21+j*3);
	  hFitSigma[j][i]->SetMarkerColor(j+1);
	  hFitSigma[j][i]->SetLineColor(j+1);
	}
      hFitSigma[1][i]->SetMarkerColor(1);
      hFitSigma[1][i]->SetLineColor(1);
      hFitSigma[0][i]->GetYaxis()->SetRangeUser(0,15);
      c = draw1D(hFitSigma[1][i],Form("Width of %s distribution for muons;p_{T} (GeV/c);#sigma",title[i]));
      //hFitSigma[1][i]->Draw("sames");

      func[i] = new TF1(Form("%s_JpsiMuon_%s_FitSigmaVsPt",data_name[1],name[i]),"[0]+[1]*exp([2]/x)",1,10);
      func[i]->SetParameters(-32.6793, 32.6034, 0.444217);
      hFitSigma[1][i]->Fit(func[i],"R0");
      func[i]->SetLineStyle(2);
      func[i]->SetLineColor(2);
      func[i]->Draw("sames");
      TLegend *leg = new TLegend(0.3,0.5,0.5,0.75);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(run_type);
      //leg->AddEntry(hFitSigma[0][i],"Data","P");
      leg->AddEntry(hFitSigma[1][i],"Embedding","P");
      //leg->AddEntry(func[i],"Fit to embedding: [0]+[1]*exp([2]/x)","L");
      leg->AddEntry(func[i],Form("Fit to embedding: %2.2f+%2.2f*exp(%2.2f/x)",func[i]->GetParameter(0), func[i]->GetParameter(1), func[i]->GetParameter(2)),"L");
      leg->Draw();
      if(savePlot)
 	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Embed_JpsiMuon_%s_FitSigmaVsPt.pdf",run_type,name[i]));
      if(saveAN)
	c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Chp2_Fit%sSigmaVsPt.pdf",name[i]));
    }

  if(saveHisto)
    {
      for(int i=0; i<2; i++)
	{
	  func[i]->Write("",TObject::kOverwrite);
	}
    }
}

//================================================
void dataVsEmbed(const bool savePlot = 0, const bool saveAN = 1)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"read");

  //==============================================
  // Compare single muon distribution
  //==============================================
  TH1F *hMuonDis[2][nHistos][nMuonPtBin];
  const int rebin[nHistos] = {4,4,5,4,5,5};
  const double minimum[nHistos] = {-50, -50, -5, 0, 0, 0};
  const double maximum[nHistos] = {50, 50, 5, 2, 50, 50};
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_distribution",name[i]),Form("%s_distribution",name[i]),1100,700);
      c->Divide(3,2);
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  c->cd(bin);
	  for(int j=0; j<2; j++)
	    {
	      hMuonDis[j][i][bin-1] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_bin%d",data_name[j],name[i],bin));
	      hMuonDis[j][i][bin-1]->SetMarkerStyle(20+j);
	      hMuonDis[j][i][bin-1]->SetMarkerColor(j+1);
	      hMuonDis[j][i][bin-1]->SetLineColor(j+1);
	      if(!(j==1&&i==3))hMuonDis[j][i][bin-1]->Rebin(rebin[i]);
	      if(i==2 && j==0) hMuonDis[j][i][bin-1]->Rebin(2);
	      hMuonDis[j][i][bin-1]->Scale(1./hMuonDis[j][i][bin-1]->Integral());
	      hMuonDis[j][i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	      hMuonDis[j][i][bin-1]->SetMaximum(1.5*hMuonDis[j][i][bin-1]->GetMaximum());
	      hMuonDis[j][i][bin-1]->SetTitle("");
	      hMuonDis[j][i][bin-1]->SetXTitle(Form("%s%s",title[i],unit[i]));
	      ScaleHistoTitle(hMuonDis[j][i][bin-1], 0.05, 0.9, 0.045, 0.05, 1.2, 0.04);
	      if(j==0) hMuonDis[j][i][bin-1]->Draw("P");
	      else     hMuonDis[j][i][bin-1]->Draw("sames HIST");
	    }

	  double up = hMuonDis[0][i][bin-1]->GetMaximum() * 0.7;
	  double low = hMuonDis[0][i][bin-1]->GetMinimum();
	  double pt = (muonPtBins[bin-1]+muonPtBins[bin])/2;
	  double min = -999, max = -999;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else if(i==2)
	    { min = -1; max = 3; }
	  else
	    { min = cuts[i]; max = 999; }

	  TLine *line = GetLine(min,low,min,up,4);
	  line->Draw();
	  TLine *line = GetLine(max,low,max,up,4);
	  line->Draw();

          TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f GeV/c",muonPtBins[bin-1],muonPtBins[bin]),0.06);
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = 0x0;
      if(i<4) leg = new TLegend(0.6,0.68,0.8,0.88);
      else    leg = new TLegend(0.2,0.68,0.4,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(hMuonDis[0][i][0],"Data","P");
      leg->AddEntry(hMuonDis[1][i][0],"Embedding","L");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/DataVsEmbed_%s_InPtBins.pdf",run_type,name[i]));
      if(saveAN)
	c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch2_Comp%s_InPtBins.pdf",name[i]));
    }
  return;

  //==============================================
  // Compare mean and sigma
  //==============================================
  TH1F *hFitSigma[2][nHistos];
  TH1F *hFitMean[2][nHistos];
  double mean_miny[3] = {-3, -2, -1};
  double mean_maxy[3] = {3, 5, 1.5};
  double sigma_miny[3] = {0, 0, 0.5};
  double sigma_maxy[3] = {15,15,1.5};
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<2; j++)
	{
	  hFitMean[j][i] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_FitMean",data_name[j],name[i]));
	  hFitMean[j][i]->SetMarkerStyle(21+j*3);
	  hFitMean[j][i]->SetMarkerColor(j+1);
	  hFitMean[j][i]->SetLineColor(j+1);

	  hFitSigma[j][i] = (TH1F*)fin->Get(Form("%s_JpsiMuon_%s_FitSigma",data_name[j],name[i]));
	  hFitSigma[j][i]->SetMarkerStyle(21+j*3);
	  hFitSigma[j][i]->SetMarkerColor(j+1);
	  hFitSigma[j][i]->SetLineColor(j+1);
	}
      hFitMean[0][i]->GetYaxis()->SetRangeUser(mean_miny[i],mean_maxy[i]);
      c = draw1D(hFitMean[0][i],Form("Mean of %s distribution for muons;p_{T} (GeV/c);Mean",title[i]));
      hFitMean[1][i]->Draw("sames");
      TLegend *leg = new TLegend(0.6,0.68,0.8,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(run_type);
      leg->AddEntry(hFitMean[0][i],"Data","P");
      leg->AddEntry(hFitMean[1][i],"Embedding","P");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/DataVsEmbed_%s_Mean.pdf",run_type,name[i]));

      hFitSigma[0][i]->GetYaxis()->SetRangeUser(sigma_miny[i],sigma_maxy[i]);
      c = draw1D(hFitSigma[0][i],Form("Width of %s distribution for muons;p_{T} (GeV/c);#sigma",title[i]));
      hFitSigma[1][i]->Draw("sames");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/DataVsEmbed_%s_Sigma.pdf",run_type,name[i]));
    }

  //==============================================
  // Compare efficiencies
  //==============================================
  TH1F *hMcEffCount[nHistos], *hMcEffFit[nHistos];
  TGraphAsymmErrors *gDataEffCount[nHistos], *gDataEffFit[nHistos];
  for(int i=0; i<3; i++)
    {
      hMcEffCount[i] = (TH1F*)fin->Get(Form("Embed_JpsiMuon_%sEff_BinCounting",name[i]));
      hMcEffFit[i] = (TH1F*)fin->Get(Form("Embed_JpsiMuon_%sEff_Fitting",name[i]));
      gDataEffCount[i] = (TGraphAsymmErrors*)fin->Get(Form("Data_JpsiMuon_%sEff_BinCounting",name[i]));
      gDataEffFit[i] = (TGraphAsymmErrors*)fin->Get(Form("Data_JpsiMuon_%sEff_Fitting",name[i]));

      hMcEffCount[i]->SetMarkerStyle(20);
      hMcEffCount[i]->SetMarkerColor(1);
      hMcEffCount[i]->SetLineColor(1);
      hMcEffCount[i]->GetYaxis()->SetRangeUser(0.5,1.1);
      ScaleHistoTitle(hMcEffCount[i],0.045,1.0,0.04,0.045,1.0,0.04,62);
      c = draw1D(hMcEffCount[i],Form("Efficency of %s cut for muons",title[i]));
      gPad->SetGridy();
      
      hMcEffFit[i]->SetMarkerStyle(25);
      hMcEffFit[i]->SetMarkerColor(1);
      hMcEffFit[i]->SetLineColor(1);
      hMcEffFit[i]->Draw("sames");

      gDataEffCount[i]->SetMarkerStyle(20);
      gDataEffCount[i]->SetMarkerColor(2);
      gDataEffCount[i]->SetLineColor(2);
      gDataEffCount[i]->Draw("samesPE");

      gDataEffFit[i]->SetMarkerStyle(25);
      gDataEffFit[i]->SetMarkerColor(2);
      gDataEffFit[i]->SetLineColor(2);
      gDataEffFit[i]->Draw("samesPE");

      TLegend *leg = new TLegend(0.5,0.15,0.8,0.4);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(run_type);
      leg->AddEntry(gDataEffCount[i],"Data: bin counting","P");
      leg->AddEntry(gDataEffFit[i],"Data: fitting","P");
      leg->AddEntry(hMcEffCount[i],"Embed: bin counting","P");
      leg->AddEntry(hMcEffFit[i],"Embed: fitting","P");
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/DataVsEmbed_%sEff.pdf",run_type,name[i]));
    }
}

//================================================
void makeEmbed(const bool savePlot = 0, const bool saveHisto = 1)
{
  TFile *fEmbed = 0x0;
  if(year==2013) fEmbed = TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
  if(year==2014) fEmbed = TFile::Open("output/Run14.AuAu200.Jpsi.Embed.root","read");

  TH2F *hMuonDisVsPt[nHistos];
  TH1F *hMuon[nHistos][nMuonPtBin];
  for(int i=0; i<nHistos; i++)
    {
      cout << Form("hTrk%sVsPt_MCreco_di_mu",name[i]) << endl;
      if(i<2) hMuonDisVsPt[i] = (TH2F*)fEmbed->Get(Form("hTrk%sVsPt_MCreco_di_mu",name[i]));
      else if(i==2) hMuonDisVsPt[i] = (TH2F*)fEmbed->Get(Form("hTrkNSigmaPi_MCreco_di_mu"));
      else hMuonDisVsPt[i] = (TH2F*)fEmbed->Get(Form("hTrk%s_MCreco_di_mu",name[i]));

      hMuonDisVsPt[i]->SetName(Form("Embed_JpsiMuon_%sVsPt",name[i]));
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(muonPtBins[bin-1]+1e-4);
	  int end_bin   = hMuonDisVsPt[i]->GetXaxis()->FindBin(muonPtBins[bin]-1e-4);
	  hMuon[i][bin-1] = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("Embed_JpsiMuon_%s_bin%d",name[i],bin),start_bin,end_bin);
	}
    }

  //==============================================
  // Display cuts
  //==============================================
  const double minimum[nHistos] = {-50, -50, -5, 0, 0, 0};
  const double maximum[nHistos] = {50, 50, 5, 2, 50, 50};
  for(int i=0; i<nHistos; i++)
    {
      hMuonDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      hMuonDisVsPt[i]->GetYaxis()->SetRangeUser(minimum[i], maximum[i]);
      c = draw2D(hMuonDisVsPt[i],Form("Embed: %s distribution for J/#psi muon;p_{T,#mu} (GeV/c);%s%s",title[i],title[i],unit[i]));
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
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      else
	{
	  TLine *line = GetLine(0,cuts[i],10,cuts[i]);
	  line->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Embed_JpsiMuon_%sVsPt.pdf",run_type,name[i]));
    }

  //==============================================
  // Fit MC distributions to extract efficiency
  //==============================================
  const int nbins = 19;
  const double xbins[nbins+1] = {1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0};
  TH1F *hFitSigma[nHistos];
  TH1F *hFitMean[nHistos];
  TH1F *hMuonFit[nHistos][nbins];
  TF1 *func[nHistos][nbins];
  TFitResultPtr ptr[nHistos][nbins];

  for(int i=0; i<3; i++)
    {
      hFitMean[i] = new TH1F(Form("Embed_JpsiMuon_%s_FitMean",name[i]),Form("Mean of %s for muons;p_{T} (GeV/c)",title[i]),nbins,xbins);
      hFitSigma[i] = new TH1F(Form("Embed_JpsiMuon_%s_FitSigma",name[i]),Form("Width of %s for muons;p_{T} (GeV/c)",title[i]),nbins,xbins);

      TCanvas *cfit = new TCanvas(Form("Fit_%s",name[i]),Form("Fit_%s",name[i]),1000,600);
      cfit->Divide(5,4);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin = hMuonDisVsPt[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  hMuonFit[i][bin-1] = (TH1F*)hMuonDisVsPt[i]->ProjectionY(Form("%s_%d",name[i],bin),start_bin,end_bin);
	  func[i][bin-1] = new TF1(Form("Embed_JpsiMuon_%s_bin%d",name[i],bin),"gaus",minimum[i], maximum[i]);
	  //func[i][bin-1] = new TF1(Form("Embed_JpsiMuon_%s_bin%d",name[i],bin),"gaus",-30,30);
	  ptr[i][bin-1] = hMuonFit[i][bin-1]->Fit(func[i][bin-1],"IR0QS");
	  hFitMean[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(1));
	  hFitMean[i]->SetBinError(bin,func[i][bin-1]->GetParError(1));
	  hFitSigma[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(2));
	  hFitSigma[i]->SetBinError(bin,func[i][bin-1]->GetParError(2));

	  cfit->cd(bin);
	  gPad->SetLogy();
	  hMuonFit[i][bin-1]->SetTitle("");
	  hMuonFit[i][bin-1]->SetMarkerStyle(20);
	  hMuonFit[i][bin-1]->SetMaximum(1.5*hMuonFit[i][bin-1]->GetMaximum());
	  hMuonFit[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	  hMuonFit[i][bin-1]->Draw("PE");
	  func[i][bin-1]->SetLineColor(4);
	  func[i][bin-1]->Draw("same");
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xbins[bin-1],xbins[bin]),0.07);
	  t1->Draw();
	}
      if(savePlot)
	cfit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Embed_JpsiMuon_%sFit.pdf",run_type,name[i]));
      hFitMean[i]->SetMarkerStyle(21);
      c = draw1D(hFitMean[i],"");

      hFitSigma[i]->SetMarkerStyle(21);
      c = draw1D(hFitSigma[i],"");
    }

  // efficiency
  TH1F *hMcEffCount[nHistos];
  TH1F *hMcEffFit[nHistos];
  for(int i=0; i<3; i++)
    {
      TH1F* hBase1  = new TH1F(Form("hBase1_%d",i), "", nbins, xbins);
      TH1F* hMatch1 = new TH1F(Form("hMatch1_%d",i), "", nbins, xbins);
      TH1F* hBase2  = new TH1F(Form("hBase2_%d",i), "", nbins, xbins);
      TH1F* hMatch2 = new TH1F(Form("hMatch2_%d",i), "", nbins, xbins);
      for(int bin=1; bin<=nbins; bin++)
	{
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

	  // bin counting
	  double all_err;
	  int low_bin = hMuonFit[i][bin-1]->FindFixBin(minimum[i]+1e-4);
	  int high_bin = hMuonFit[i][bin-1]->FindFixBin(maximum[i]-1e-4);
	  double all = hMuonFit[i][bin-1]->IntegralAndError(low_bin, high_bin,all_err);
	  double acc_err;
	  low_bin = hMuonFit[i][bin-1]->FindFixBin(min);
	  high_bin = hMuonFit[i][bin-1]->FindFixBin(max);
	  double acc = hMuonFit[i][bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
	  double acc_corr = acc - (hMuonFit[i][bin-1]->GetBinContent(low_bin)/hMuonFit[i][bin-1]->GetBinWidth(low_bin)*(min-hMuonFit[i][bin-1]->GetXaxis()->GetBinLowEdge(low_bin)) + hMuonFit[i][bin-1]->GetBinContent(high_bin)/hMuonFit[i][bin-1]->GetBinWidth(high_bin)*(hMuonFit[i][bin-1]->GetXaxis()->GetBinUpEdge(high_bin)-max) );
	  double acc_corr = acc;
	  hBase1->SetBinContent(bin,all);
	  hBase1->SetBinError(bin,all_err);
	  hMatch1->SetBinContent(bin,acc_corr);
	  hMatch1->SetBinError(bin,acc_corr*acc_err/acc);

	  // fitting
	  all = func[i][bin-1]->Integral(minimum[i],maximum[i]);
	  all_err = func[i][bin-1]->IntegralError(minimum[i],maximum[i],func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  acc = func[i][bin-1]->Integral(min,max);
	  acc_err = func[i][bin-1]->IntegralError(min,max,func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  hBase2->SetBinContent(bin,all);
	  hBase2->SetBinError(bin,all_err);
	  hMatch2->SetBinContent(bin,acc);
	  hMatch2->SetBinError(bin,acc_err);
	}

      // Counting vs fitting 
      hMcEffCount[i] = (TH1F*)hMatch1->Clone(Form("Embed_JpsiMuon_%sEff_BinCounting",name[i]));
      hMcEffCount[i]->Divide(hBase1);

      hMcEffFit[i] = (TH1F*)hMatch2->Clone(Form("Embed_JpsiMuon_%sEff_Fitting",name[i]));
      hMcEffFit[i]->Divide(hBase2);

      TCanvas* cEff = new TCanvas(Form("%sEff_CountVsFit",name[i]),Form("%sEff_CountVsFit",name[i]),800,600);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      gPad->SetGridy();
      hMcEffCount[i]->SetMarkerStyle(20);
      hMcEffCount[i]->GetYaxis()->SetRangeUser(0.5,1.05);
      hMcEffCount[i]->GetXaxis()->SetRangeUser(1.3,10);
      hMcEffCount[i]->SetMarkerSize(1.2);
      hMcEffCount[i]->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
      ScaleHistoTitle(hMcEffCount[i],0.05,1.0,0.045,0.05,1.0,0.045,62);
      hMcEffCount[i]->Draw("PE");
      hMcEffFit[i]->SetMarkerStyle(25);
      hMcEffFit[i]->SetMarkerColor(2);
      hMcEffFit[i]->SetLineColor(2);
      hMcEffFit[i]->SetMarkerSize(1.2);
      hMcEffFit[i]->Draw("PZsames");
      TPaveText *t1 = GetTitleText(Form("Embed: %s efficiency",title[i]),0.04);
      t1->Draw();
      TLegend *leg = new TLegend(0.3,0.3,0.5,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader("Embedding");
      leg->AddEntry(hMcEffCount[i],Form("Bin counting"),"p");
      leg->AddEntry(hMcEffFit[i],Form("Fitting"),"p");
      leg->Draw();
      if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Embed_JpsiMuon_%sEff_FitVsCount.pdf",run_type,name[i]));
    }

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
      for(int i=0; i<nHistos; i++)
	{
	  hMuonDisVsPt[i]->Write(Form("Embed_JpsiMuon_%sVsPt",name[i]),TObject::kOverwrite);
	  for(int bin=1; bin<=nMuonPtBin; bin++)
	    {
	      hMuon[i][bin-1]->Write("",TObject::kOverwrite);
	    }
	  if(i<3)
	    {
	      hFitMean[i]->Write("",TObject::kOverwrite);
	      hFitSigma[i]->Write("",TObject::kOverwrite);
	      hMcEffCount[i]->Write("",TObject::kOverwrite);
	      hMcEffFit[i]->Write("",TObject::kOverwrite);
	    }
	}
    }
}


//================================================
void makeData(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fdata = 0x0;
  //if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.jpsi.root","read");
  if(year==2014) fdata = TFile::Open("./output/Pico.Run14.AuAu200.JpsiMuon.dtof0.4.root","read");
  if(year==2015) fdata = TFile::Open("./output/Pico.Run15.pp200.jpsi.muon.root","read");


  //==============================================
  // compare invariant mass
  //==============================================
  const int nbins = 5;
  const double xbins[nbins+1] = {0,1.0,2.0,3.0,5.0,10.0};
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  double counts[nHistos][2][3] = {0};
  TH2F *hInvMassVsPtUL[nHistos];
  TH2F *hInvMassVsPtLS[nHistos];
  for(int i=0; i<nHistos; i++)
    {
      hInvMassVsPtUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_UL_di_mu",name[i]));
      hInvMassVsPtUL[i]->Sumw2();
      hInvMassVsPtLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_%s_LS_di_mu",name[i]));
      hInvMassVsPtLS[i]->Sumw2();
    }
  
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("InvMass_%s",name[i]),Form("InvMass_%s",name[i]),1100,700);
      c->Divide(3,2);
      TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.06);
      leg->SetHeader(run_type);
      for(int bin=1; bin<=nbins; bin++)
	{
	  int start_bin = hInvMassVsPtUL[i]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
	  int end_bin   = hInvMassVsPtUL[i]->GetXaxis()->FindBin(xbins[bin]-1e-4);
	  TH1F *hInvMassUL = (TH1F*)hInvMassVsPtUL[i]->ProjectionY(Form("Data_%s_InvMassUL_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassUL->SetMarkerStyle(20);
	  hInvMassUL->SetMarkerStyle(20);
	  hInvMassUL->Rebin(5);
	  hInvMassUL->SetMaximum(1.5*hInvMassUL->GetMaximum());
	  hInvMassUL->SetMinimum(0);
	  hInvMassUL->GetXaxis()->SetRangeUser(2.2,4);
	  hInvMassUL->SetTitle("");
	  if(bin==1) leg->AddEntry(hInvMassUL,"Unlike-sign","PL");

	  TH1F *hInvMassLS = (TH1F*)hInvMassVsPtLS[i]->ProjectionY(Form("Data_%s_InvMassLS_bin%d",name[i],bin),start_bin,end_bin);
	  hInvMassLS->SetMarkerStyle(24);
	  hInvMassLS->SetMarkerColor(2);
	  hInvMassLS->SetLineColor(2);
	  hInvMassLS->Rebin(5);
	  if(bin==1) leg->AddEntry(hInvMassLS,"Like-sign","PL");

	  c->cd(bin);
	  hInvMassUL->Draw("P");
	  // draw signal region
	  double binwidth = hInvMassUL->GetBinWidth(1);
	  for(int itmp=0; itmp<3; itmp++)
	    {
	      counts[i][0][itmp] = 0;
	      counts[i][1][itmp] = 0;
	      TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMassUL->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	      for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
		{
		  int jbin = hInvMassUL->FindFixBin(htmp->GetBinCenter(ibin));
		  htmp->SetBinContent(ibin,hInvMassUL->GetBinContent(jbin));
		  counts[i][0][itmp] += hInvMassUL->GetBinContent(jbin);
		  counts[i][1][itmp] += hInvMassLS->GetBinContent(jbin);
		}
	      if(itmp==0) htmp->SetFillColor(7);
	      else htmp->SetFillColor(5);
	      htmp->Draw("sames");
	      if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	      if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	    }
	  hInvMassLS->Draw("samesP");
	  TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbins[bin-1],xbins[bin]),0.05);
	  t1->Draw();
	  double count_ul = counts[i][0][1] + counts[i][0][2];
	  double count_ls = counts[i][1][1] + counts[i][1][2];
	  double scale = count_ul/count_ls;
	  double error = scale * sqrt(1/count_ul + 1/count_ls);
	  TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
	  t1->SetTextAlign(11);
	  t1->AddText(Form("S/B = 1:%4.2f",counts[i][1][0]/(counts[i][0][0]-counts[i][1][0])));
	  t1->AddText(Form("UL/LS = %4.2f#pm%4.3f",scale,error));
	  t1->Draw();
	}
      c->cd(6);
      leg->Draw();
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%s_InvMassLSvsUL.pdf",run_type,name[i]));
    }

  //==============================================
  // Muon distribution UL vs LS
  //==============================================
  const int rebin[nHistos] = {8,8,10,1,1,1};
  const double minimum[nHistos] = {-50, -50, -5, 0, 0, 0};
  const double maximum[nHistos] = {50, 50, 5, 2, 50, 50};

  TH2F *hJpsiMassVsMuonPtUL[nHistos];
  TH2F *hJpsiMassVsMuonPtLS[nHistos];
  TH1F *hMuonPtSideUL[nHistos];
  TH1F *hMuonPtSideLS[nHistos];
  TH2F *hDataUL[nHistos];
  TH2F *hDataLS[nHistos];
  TH2F *hDataDisVsPt[nHistos];
  TCanvas *c = new TCanvas(Form("MuonSB_UL_vs_LS"),Form("MuonSB_UL_vs_LS"),1100,700);
  c->Divide(3,2);
  for(int i=0; i<nHistos; i++)
    {
      // Get the muon distribution in the side-band region
      hJpsiMassVsMuonPtUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_%s_UL_di_mu",name[i]));
      hJpsiMassVsMuonPtUL[i]->Sumw2();
      hMuonPtSideUL[i] = (TH1F*)hJpsiMassVsMuonPtUL[i]->ProjectionX(Form("hMuonPt_%s_SideBand_UL",name[i]));
      hMuonPtSideUL[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPtUL[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPtUL[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  TH1F *htmp = (TH1F*)hJpsiMassVsMuonPtUL[i]->ProjectionX(Form("hMuonPt_%s_SideBand%d_UL",name[i],j),low_bin,up_bin);
	  hMuonPtSideUL[i]->Add(htmp);
	}

      hJpsiMassVsMuonPtLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_%s_LS_di_mu",name[i]));
      hJpsiMassVsMuonPtLS[i]->Sumw2();
      hMuonPtSideLS[i] = (TH1F*)hJpsiMassVsMuonPtLS[i]->ProjectionX(Form("hMuonPt_%s_SideBand_LS",name[i]));
      hMuonPtSideLS[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPtLS[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPtLS[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  TH1F *htmp = (TH1F*)hJpsiMassVsMuonPtLS[i]->ProjectionX(Form("hMuonPt_%s_SideBand%d_LS",name[i],j),low_bin,up_bin);
	  hMuonPtSideLS[i]->Add(htmp);
	}
      c->cd(i+1);
      gPad->SetLogy();
      hMuonPtSideUL[i]->GetXaxis()->SetRangeUser(0,10);
      hMuonPtSideUL[i]->SetMarkerStyle(20);
      hMuonPtSideUL[i]->Draw("PE");
      hMuonPtSideLS[i]->SetMarkerStyle(24);
      hMuonPtSideLS[i]->SetMarkerColor(2);
      hMuonPtSideLS[i]->SetLineColor(2);
      hMuonPtSideLS[i]->Draw("samesPE");

      hDataUL[i] = (TH2F*)fdata->Get(Form("mhJpsiMuon_%s_UL_di_mu",name[i]));
      hDataUL[i]->Sumw2();
      hDataLS[i] = (TH2F*)fdata->Get(Form("mhJpsiMuon_%s_LS_di_mu",name[i]));
      hDataLS[i]->Sumw2();
      hDataDisVsPt[i]  = (TH2F*)hDataUL[i]->Clone(Form("mhJpsiMuon%sVsPt",name[i]));
      hDataDisVsPt[i]->Add(hDataLS[i], -1);
    }

  TH1F *hUL[nHistos][nMuonPtBin];
  TH1F *hLS[nHistos][nMuonPtBin];
  TH1F *hMuonFineBin[nHistos][nMuonPtBin];
  TH1F *hMuon[nHistos][nMuonPtBin];
  double mean_pt[nHistos][nMuonPtBin];
  double mean_pt_err[nHistos][nMuonPtBin];
  for(int i=0; i<nHistos; i++)
    {
      TCanvas *c = new TCanvas(Form("%s_UL_vs_LS",name[i]),Form("%s_UL_vs_LS",name[i]),1100,700);
      c->Divide(3,2);

      // calculate mean pt of each bin
      TH1F *htmp = (TH1F*)hDataDisVsPt[i]->ProjectionX(Form("%s_tmp",hDataDisVsPt[i]->GetName()),1,hDataDisVsPt[i]->GetNbinsX());
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  htmp->GetXaxis()->SetRangeUser(muonPtBins[bin-1]+1e-4, muonPtBins[bin]-1e-4);
	  mean_pt[i][bin-1] = htmp->GetMean();
	  mean_pt_err[i][bin-1] = htmp->GetMeanError();

	  int start_bin = hDataUL[i]->GetXaxis()->FindBin(muonPtBins[bin-1]+1e-4);
	  int end_bin   = hDataUL[i]->GetXaxis()->FindBin(muonPtBins[bin]-1e-4);
	  hUL[i][bin-1] = (TH1F*)hDataUL[i]->ProjectionY(Form("Data_JpsiMuon_%s_UL_bin%d",name[i],bin),start_bin,end_bin);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetMarkerStyle(20);
	  hUL[i][bin-1]->SetTitle("");

	  hLS[i][bin-1] = (TH1F*)hDataLS[i]->ProjectionY(Form("Data_JpsiMuon_%s_LS_bin%d",name[i],bin),start_bin,end_bin);
	  hLS[i][bin-1]->SetMarkerStyle(24);
	  hLS[i][bin-1]->SetMarkerColor(2);
	  hLS[i][bin-1]->SetLineColor(2);

	  start_bin = hMuonPtSideUL[i]->GetXaxis()->FindBin(muonPtBins[bin-1]+1e-4);
	  end_bin   = hMuonPtSideUL[i]->GetXaxis()->FindBin(muonPtBins[bin]-1e-4);
	  double scale = hMuonPtSideUL[i]->Integral(start_bin,end_bin)/hMuonPtSideLS[i]->Integral(start_bin,end_bin);
	  hMuonFineBin[i][bin-1] = (TH1F*)hUL[i][bin-1]->Clone(Form("Data_JpsiMuon_%s_bin%d",name[i],bin));
	  hMuonFineBin[i][bin-1]->Add(hLS[i][bin-1],-1.0*scale);

	  hMuon[i][bin-1] = (TH1F*)hMuonFineBin[i][bin-1]->Clone(Form("Data_JpsiMuon_%s_bin%d_Rebin",name[i],bin));
	  if(i==0 && bin>=5) hMuon[i][bin-1]->Rebin(4);
	  else hMuon[i][bin-1]->Rebin(rebin[i]);

	  c->cd(bin);
	  hUL[i][bin-1]->Rebin(rebin[i]);
	  hUL[i][bin-1]->SetMaximum(1.5*hUL[i][bin-1]->GetMaximum());
	  hUL[i][bin-1]->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
      	  hUL[i][bin-1]->Draw("P");
	  hLS[i][bin-1]->Rebin(rebin[i]);
	  hLS[i][bin-1]->Draw("samesP");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",muonPtBins[bin-1],muonPtBins[bin]),0.06);
          t1->Draw();
	  t1 = GetPaveText(0.6,0.8,0.75,0.8,0.05);
	  t1->AddText(Form("Scale = %4.3f",scale));
          t1->Draw();
	}
      c->cd(1);
      TLegend *leg = new TLegend(0.15,0.62,0.35,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(run_type);
      leg->AddEntry(hUL[i][0],"Unlike-sign","PL");
      leg->AddEntry(hLS[i][0],"Like-sign","PL");
      leg->Draw();

      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%s_ULvsLS.pdf",run_type,name[i]));
    }

  for(int i=0; i<nHistos; i++)
    {
      hDataDisVsPt[i]->GetXaxis()->SetRangeUser(0,10);
      hDataDisVsPt[i]->GetYaxis()->SetRangeUser(minimum[i], maximum[i]);
      c = draw2D(hDataDisVsPt[i],Form("%s: %s distribution for J/#psi muon;p_{T,#mu} (GeV/c);%s%s",run_type,title[i],title[i],unit[i]));
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
	  TLine *line = GetLine(0,-1,10,-1);
	  line->Draw();
	  line = GetLine(0,3,10,3);
	  line->Draw(); 
	}
      else
	{
	  TLine *line = GetLine(0,cuts[i],10,cuts[i]);
	  line->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sVsPt.pdf",run_type,name[i]));
    }
 
  //==============================================
  // Fit data distributions to extract efficiency
  //==============================================
  TH1F *hFitDataMean[nHistos];
  TH1F *hFitDataSigma[nHistos];
  TF1 *func[nHistos][nMuonPtBin];
  TFitResultPtr ptr[nHistos][nMuonPtBin];
  for(int i=0; i<3; i++)
    {
      hFitDataMean[i] = new TH1F(Form("Data_JpsiMuon_%s_FitMean",name[i]),Form("Mean of %s;p_{T} (GeV/c)",title[i]),nMuonPtBin,muonPtBins);
      hFitDataSigma[i] = new TH1F(Form("Data_JpsiMuon_%s_FitSigma",name[i]),Form("Width of %s;p_{T} (GeV/c)",title[i]),nMuonPtBin,muonPtBins);

      TCanvas *c = new TCanvas(Form("Fit_%s",name[i]),Form("Fit_%s",name[i]),1100,700);
      c->Divide(3,2);

      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  TH1F *hFit = (TH1F*)hMuon[i][bin-1]->Clone(Form("Fit_%s",hMuon[i][bin-1]->GetName()));
	  func[i][bin-1] = new TF1(Form("Data_JpsiMuon_%sFit_bin%d",name[i],bin),"gaus",minimum[i], maximum[i]);
	  if(i==1 && bin==1) ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"R0QS");
	  else ptr[i][bin-1] = hFit->Fit(func[i][bin-1],"IR0QS");
	  hFitDataMean[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(1));
	  hFitDataMean[i]->SetBinError(bin,func[i][bin-1]->GetParError(1));
	  hFitDataSigma[i]->SetBinContent(bin,func[i][bin-1]->GetParameter(2));
	  hFitDataSigma[i]->SetBinError(bin,func[i][bin-1]->GetParError(2));

	  c->cd(bin);
	  hFit->SetMaximum(1.5*hFit->GetMaximum());
	  hFit->GetXaxis()->SetRangeUser(minimum[i], maximum[i]);
	  hFit->Draw();
	  func[i][bin-1]->SetLineColor(4);
	  func[i][bin-1]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",muonPtBins[bin-1],muonPtBins[bin]),0.06);
	  t1->Draw();
	}
      if(savePlot)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sFit.pdf",run_type,name[i]));
    }
  
  // efficiency
  TGraphAsymmErrors *gCountDataEff[nHistos];
  TGraphAsymmErrors *gFitDataEff[nHistos];
  for(int i=0; i<3; i++)
    {
      TH1F* hBase1  = new TH1F(Form("hBase1_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hMatch1 = new TH1F(Form("hMatch1_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hBase2  = new TH1F(Form("hBase2_%d",i), "", nMuonPtBin, muonPtBins);
      TH1F* hMatch2 = new TH1F(Form("hMatch2_%d",i), "", nMuonPtBin, muonPtBins);
      for(int bin=1; bin<=nMuonPtBin; bin++)
	{
	  double min = -999, max = -999;
	  double pt = mean_pt[i][bin-1];;
	  if(i<2)
	    {
	      double sigma = fResVsPt[i]->Eval(pt);
	      if(pt<ptbound) min = -1 * nsigma1 * sigma;
	      else           min = -1 * nsigma2 * sigma;
	      max = -1 * min;
	    }
	  else
	    { min = -1; max = 3; }

	  // bin counting
	  double all_err;
	  int low_bin = hMuonFineBin[i][bin-1]->FindFixBin(minimum[i]+1e-4);
	  int high_bin = hMuonFineBin[i][bin-1]->FindFixBin(maximum[i]-1e-4);
	  double all = hMuonFineBin[i][bin-1]->IntegralAndError(low_bin, high_bin,all_err);
	  double acc_err;
	  low_bin = hMuonFineBin[i][bin-1]->FindFixBin(min+1e-4);
	  high_bin = hMuonFineBin[i][bin-1]->FindFixBin(max-1e-4);
	  double acc = hMuonFineBin[i][bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
	  double acc_corr = acc - (hMuonFineBin[i][bin-1]->GetBinContent(low_bin)/hMuonFineBin[i][bin-1]->GetBinWidth(low_bin)*(min-hMuonFineBin[i][bin-1]->GetXaxis()->GetBinLowEdge(low_bin)) + hMuonFineBin[i][bin-1]->GetBinContent(high_bin)/hMuonFineBin[i][bin-1]->GetBinWidth(high_bin)*(hMuonFineBin[i][bin-1]->GetXaxis()->GetBinUpEdge(high_bin)-max) );
	  //double acc_corr = acc;
	  cout << acc_corr/all << "  " << acc/all << endl;
	  hBase1->SetBinContent(bin,all);
	  hBase1->SetBinError(bin,all_err);
	  hMatch1->SetBinContent(bin,acc_corr);
	  hMatch1->SetBinError(bin,acc_corr*acc_err/acc);

	  // fitting
	  all = func[i][bin-1]->Integral(minimum[i],maximum[i]);
	  all_err = func[i][bin-1]->IntegralError(minimum[i],maximum[i],func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  acc = func[i][bin-1]->Integral(min,max);
	  acc_err = func[i][bin-1]->IntegralError(min,max,func[i][bin-1]->GetParameters(),ptr[i][bin-1]->GetCovarianceMatrix().GetMatrixArray());
	  hBase2->SetBinContent(bin,all);
	  hBase2->SetBinError(bin,all_err);
	  hMatch2->SetBinContent(bin,acc);
	  hMatch2->SetBinError(bin,acc_err);
	}

      // bin counting    
      TH1F *hMatch1_corr = (TH1F*)hMatch1->Clone(Form("%_tmp",hMatch1->GetName()));
      for(int ibin=1; ibin<=hMatch1->GetNbinsX(); ibin++)
	{
	  if(hMatch1->GetBinContent(ibin)>hBase1->GetBinContent(ibin))
	    {
	      hMatch1_corr->SetBinContent(ibin,hBase1->GetBinContent(ibin));
	    }
	}
      gCountDataEff[i] = new TGraphAsymmErrors(hMatch1_corr, hBase1,"cl=0.683 b(1,1) mode");
      gCountDataEff[i]->SetName(Form("Data_JpsiMuon_%sEff_BinCounting",name[i]));
      double x,y;
      hMatch1->Divide(hBase1);
      for(int ipoint=0; ipoint<gCountDataEff[i]->GetN(); ipoint++)
	{
	  gCountDataEff[i]->GetPoint(ipoint,x,y);
	  gCountDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	  gCountDataEff[i]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	  gCountDataEff[i]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);

	  if(hMatch1->GetBinContent(ipoint+1)>1)
	    {
	      gCountDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],hMatch1->GetBinContent(ipoint+1));	 
	      gCountDataEff[i]->SetPointEYhigh(ipoint,hMatch1->GetBinError(ipoint+1));
	      gCountDataEff[i]->SetPointEYlow(ipoint,hMatch1->GetBinError(ipoint+1));   
	    }
	}

      // fitting 
      gFitDataEff[i] = new TGraphAsymmErrors(hMatch2, hBase2,"cl=0.683 b(1,1) mode");
      gFitDataEff[i]->SetName(Form("Data_JpsiMuon_%sEff_Fitting",name[i]));
      for(int ipoint=0; ipoint<gFitDataEff[i]->GetN(); ipoint++)
	{
	  gFitDataEff[i]->GetPoint(ipoint,x,y);
	  gFitDataEff[i]->SetPoint(ipoint,mean_pt[i][ipoint],y);
	  gFitDataEff[i]->SetPointEXhigh(ipoint,mean_pt_err[i][ipoint]);
	  gFitDataEff[i]->SetPointEXlow(ipoint,mean_pt_err[i][ipoint]);	 
	}

      // Counting vs fitting
      TCanvas* cEff = new TCanvas(Form("%sEff_CountVsFit",name[i]),Form("%sEff_CountVsFit",name[i]),800,600);
      SetPadMargin(gPad,0.13,0.13,0.05,0.1);
      gPad->SetGridy();
      gCountDataEff[i]->SetMarkerStyle(20);
      gCountDataEff[i]->GetYaxis()->SetRangeUser(0.5,1.2);
      gCountDataEff[i]->GetXaxis()->SetRangeUser(1.3,10);
      gCountDataEff[i]->SetMarkerSize(1.2);
      gCountDataEff[i]->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
      ScaleHistoTitle(gCountDataEff[i]->GetHistogram(),0.05,1.0,0.045,0.05,1.0,0.045,62);
      gCountDataEff[i]->Draw("APZ");
      gFitDataEff[i]->SetMarkerStyle(25);
      gFitDataEff[i]->SetMarkerColor(2);
      gFitDataEff[i]->SetLineColor(2);
      gFitDataEff[i]->SetMarkerSize(1.2);
      gFitDataEff[i]->Draw("PZsames");
      TPaveText *t1 = GetTitleText(Form("%s: %s efficiency",run_type,title[i]),0.04);
      t1->Draw();
      TLegend *leg = new TLegend(0.3,0.3,0.5,0.45);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader(run_type);
      leg->AddEntry(gCountDataEff[i],Form("Bin counting"),"p");
      leg->AddEntry(gFitDataEff[i],Form("Fitting"),"p");
      leg->Draw();
      if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbe/Data_JpsiMuon_%sEff_FitVsCount.pdf",run_type,name[i]));
    }

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.PIDcuts.root",run_type),"update");
      for(int i=0; i<nHistos; i++)
	{
	  hDataDisVsPt[i]->Write(Form("Data_JpsiMuon_%sVsPt",name[i]),TObject::kOverwrite);
	  for(int bin=1; bin<=nMuonPtBin; bin++)
	    {
	      hMuonFineBin[i][bin-1]->Write("",TObject::kOverwrite);
	      hMuon[i][bin-1]->Write("",TObject::kOverwrite);
	      if(i<3)
		func[i][bin-1]->Write("",TObject::kOverwrite);
	    }
	  if(i<3)
	    {
	      hFitDataMean[i]->Write("",TObject::kOverwrite);
	      hFitDataSigma[i]->Write("",TObject::kOverwrite);
	      gCountDataEff[i]->Write("",TObject::kOverwrite);
	      gFitDataEff[i]->Write("",TObject::kOverwrite);
	    }
	}
    }
}
