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
  makeEmbed();
  //dataVsEmbed();
  //DyDzCut();
  //Run16nSigmaPiEff();
}

//================================================
void DyDzCut(const bool savePlot = 0, const bool saveHisto = 0, const bool saveAN = 0)
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
void makeEmbed(const bool savePlot = 1, const bool saveHisto = 0)
{
  TFile *fEmbed = 0x0;
  if(year==2013) fEmbed = TFile::Open("output/Run13.pp500.jpsi.Embed.root","read");
  if(year==2014) fEmbed = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root","read");

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
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/Embed_JpsiMuon_%sVsPt.pdf",run_type,name[i]));
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
	cfit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/Embed_JpsiMuon_%sFit.pdf",run_type,name[i]));
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
      if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_DataVsEmbed/Embed_JpsiMuon_%sEff_FitVsCount.pdf",run_type,name[i]));
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

