const char *run_type = "Run14_AuAu200";

//================================================
void plot_AuAu200()
{  
  publication();
}

//================================================
void publication(const bool savePlot = 1, const bool saveHisto = 1)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *fout = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","update");

  const int nbins = nPtBins -1;
  double xbins[nbins+1];
  for(int i=0; i<nbins; i++)
    xbins[i] = ptBins_low[i+1];
  xbins[nbins] = ptBins_high[nbins];

  TFile *f = TFile::Open("Rootfiles/Publication.Jpsi.pp200GeV.root","read");
  gStyle->SetOptStat(0);

  // pp data
  TH1F *hpp = new TH1F("pp200_Jpsi",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",15,0,15);
  hpp->GetYaxis()->SetRangeUser(1e-6,10);
  TCanvas *c = draw1D(hpp,"",kTRUE);
  const int npp = 11;
  double xpp[npp] = {2.25, 2.75, 3.25, 3.75, 4.5, 5.5, 6.5, 7.5, 9, 11, 13};
  double exlpp[npp];
  double sxlpp[npp];
  for(int i=0; i<npp; i++) { exlpp[i] = 0; sxlpp[i] = 0.2; }
  double ypp[npp] = {0.68, 0.318, 0.187, 0.1032, 0.0334, 0.00905, 0.00154, 0.00084, 2.005e-4, 4.55e-5, 9.71e-6};
  double eylpp[npp] = {0.14, 0.050, 0.023, 0.0118, 0.0034, 0.00117, 0.00042, 0.00017, 2.42e-5, 7.2e-6, 2.41e-6};
  double sylpp[npp] = {0.06, 0.028, 0.018, 0.0089, 0.0028, 0.00077, 0.00013, 0.00025, 5.21e-5, 1.18e-5, 2.62e-6};
  double syhpp[npp] = {0.03, 0.007, 0.009, 0.0018, 0.0004, 0.00011, 0.00003, 0.00013, 1.6e-6, 4e-7, 1.6e-7};

  TGraphAsymmErrors *gHighPtPP = new TGraphAsymmErrors(npp, xpp, ypp, exlpp, exlpp, eylpp, eylpp);
  gHighPtPP->SetName("Jpsi_xsec_pp200_highPt");
  gHighPtPP->SetMarkerStyle(20);
  gHighPtPP->SetMarkerSize(1.5);
  gHighPtPP->Draw("sames PE");

  TGraphAsymmErrors *gHighPtPPSys = new TGraphAsymmErrors(npp, xpp, ypp, sxlpp, sxlpp, sylpp, syhpp);
  gHighPtPPSys->SetName("Jpsi_xsec_pp200_highPt_systematics");
  gHighPtPPSys->SetMarkerStyle(20);
  gHighPtPPSys->SetMarkerSize(0);
  gHighPtPPSys->SetFillStyle(0);
  gHighPtPPSys->Draw("sameE5");

  TH1F *hLowPtPP = (TH1F*)f->Get("hInvariantYield_pp");
  hLowPtPP->SetMarkerStyle(24);
  hLowPtPP->SetMarkerSize(1.5);
  hLowPtPP->Draw("sames");

  TH1F *hJpsiPP = new TH1F("Jpsi_xsec_pp200_combined",";p_{T} (GeV/c);Bd^{2}#sigma/(2#pip_{T}dp_{T}dy) [nb/(GeV/c)^{2}]",nbins,xbins);
  for(int ibin=1; ibin<=4; ibin++)
    {
      hJpsiPP->SetBinContent(ibin,hLowPtPP->GetBinContent(ibin));
      hJpsiPP->SetBinError(ibin,hLowPtPP->GetBinError(ibin));
    }
  double val = (ypp[4]*xpp[4]+ypp[5]*xpp[5])*1./(2.*5.);
  double err = sqrt(eylpp[4]*eylpp[4]*xpp[4]*xpp[4]+eylpp[5]*eylpp[5]*xpp[5]*xpp[5])*1./(2.*5.);
  double syl = (sylpp[4]*xpp[4]+sylpp[5]*xpp[5])*1./(2.*5.);
  double syh = (syhpp[4]*xpp[4]+syhpp[5]*xpp[5])*1./(2.*5.);
  double sys = (syl+syh)/2.;
  hJpsiPP->SetBinContent(5,val);
  hJpsiPP->SetBinError(5,sqrt(err*err+sys*sys));
  val = (ypp[6]*xpp[6]*1.+ypp[7]*xpp[7]*1.+ypp[8]*xpp[8]*2.)/(4.*8.);
  err = sqrt(eylpp[6]*eylpp[6]*xpp[6]*xpp[6]+eylpp[7]*eylpp[7]*xpp[7]*xpp[7]+eylpp[8]*eylpp[8]*xpp[8]*xpp[8]*2.*2.)/(4.*8.);
  syl = (sylpp[6]*xpp[6]*1.+sylpp[7]*xpp[7]*1.+sylpp[8]*xpp[8]*2.)/(4.*8.);
  syh = (syhpp[6]*xpp[6]*1.+syhpp[7]*xpp[7]*1.+syhpp[8]*xpp[8]*2.)/(4.*8.);
  sys = (syl+syh)/2.;
  hJpsiPP->SetBinContent(6,val);
  hJpsiPP->SetBinError(6,sqrt(err*err+sys*sys));

  hJpsiPP->SetMarkerStyle(21);
  hJpsiPP->SetMarkerSize(1.5);
  hJpsiPP->SetMarkerColor(2);
  hJpsiPP->SetLineColor(2);
  hJpsiPP->Draw("sames");
  leg = new TLegend(0.5,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("p+p @ 200 GeV");
  leg->AddEntry(gHighPtPP,"STAR high p_{T}","PL");
  leg->AddEntry(hJpsiPP,"PHENIX + STAR","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_xsec_pp.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_xsec_pp.png",run_type));
    }

  // AuAu
  const int nAuAuLowPt = 5;
  double xAuAuLowPt[nCentBins][nAuAuLowPt];
  double exlAuAuLowPt[nCentBins][nAuAuLowPt];
  double sxlAuAuLowPt[nCentBins][nAuAuLowPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuLowPt; j++)
	{
	  xAuAuLowPt[i][j] = 0.5 + j; 
	  exlAuAuLowPt[i][j] = 0; 
	  sxlAuAuLowPt[i][j] = 0.2; 
	}
    }
  double yAuAuLowPt[nCentBins][nAuAuLowPt];
  double eylAuAuLowPt[nCentBins][nAuAuLowPt];
  double sylAuAuLowPt[nCentBins][nAuAuLowPt];
  double syhAuAuLowPt[nCentBins][nAuAuLowPt];

  double xRaaLowPt[nCentBins][nAuAuLowPt];
  double exlRaaLowPt[nCentBins][nAuAuLowPt];
  double sxlRaaLowPt[nCentBins][nAuAuLowPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuLowPt; j++)
	{
	  xRaaLowPt[i][j] = 0.5 + j; 
	  exlRaaLowPt[i][j] = 0; 
	  sxlRaaLowPt[i][j] = 0.2; 
	}
    }
  double yRaaLowPt[nCentBins][nAuAuLowPt];
  double eylRaaLowPt[nCentBins][nAuAuLowPt];
  double sylRaaLowPt[nCentBins][nAuAuLowPt];
  double syhRaaLowPt[nCentBins][nAuAuLowPt];

  ifstream data_in;
  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.LowPt.txt");
  char data[10][256];
  double meanx, meany, staty, sysyhigh, sysylow, raa, statraa, sysraahigh, sysraalow, sysglobal;
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<10; i++)
	data_in >> data[i];

      int cent = counter/nAuAuLowPt;
      int ptbin = counter%nAuAuLowPt;
      yAuAuLowPt[cent][ptbin] = atof(data[1]);
      eylAuAuLowPt[cent][ptbin] = atof(data[2]);
      syhAuAuLowPt[cent][ptbin] = atof(data[3]);
      sylAuAuLowPt[cent][ptbin] = fabs(atof(data[4]));

      yRaaLowPt[cent][ptbin] = atof(data[5]);
      eylRaaLowPt[cent][ptbin] = atof(data[6]);
      syhRaaLowPt[cent][ptbin] = atof(data[7]);
      sylRaaLowPt[cent][ptbin] = fabs(atof(data[8]));
      //cout << sylAuAuLowPt[cent][ptbin] << endl;
      //cout << counter << "  " << cent << "  " << ptbin << endl;
      counter++;
      if(counter==nAuAuLowPt*nCentBins) break; 
    }
  data_in.close();

  const int nAuAuHighPt = 6;
  double xAuAuHighPt[nCentBins][nAuAuHighPt];
  double exlAuAuHighPt[nCentBins][nAuAuHighPt];
  double sxlAuAuHighPt[nCentBins][nAuAuHighPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nAuAuHighPt; j++)
	{
	  exlAuAuHighPt[i][j] = 0; 
	  sxlAuAuHighPt[i][j] = 0.2; 
	}
    }
  double yAuAuHighPt[nCentBins][nAuAuHighPt];
  double eylAuAuHighPt[nCentBins][nAuAuHighPt];
  double sylAuAuHighPt[nCentBins][nAuAuHighPt];
  double syhAuAuHighPt[nCentBins][nAuAuHighPt];

  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.HighPt.txt");
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<6; i++)
	data_in >> data[i];

      int cent = counter/nAuAuHighPt;
      int ptbin = counter%nAuAuHighPt;
      xAuAuHighPt[cent][ptbin] = atof(data[1]);
      yAuAuHighPt[cent][ptbin] = atof(data[2]);
      eylAuAuHighPt[cent][ptbin] = atof(data[3]);
      sylAuAuHighPt[cent][ptbin] = atof(data[4]);
      syhAuAuHighPt[cent][ptbin] = fabs(atof(data[5]));

      counter++;
      if(counter==nAuAuHighPt*nCentBins) break;
    }
  data_in.close();


  const int nRaaHighPt = 5;
  double xRaaHighPt[nCentBins][nRaaHighPt];
  double exlRaaHighPt[nCentBins][nRaaHighPt];
  double sxlRaaHighPt[nCentBins][nRaaHighPt];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<nRaaHighPt; j++)
	{
	  exlRaaHighPt[i][j] = 0; 
	  sxlRaaHighPt[i][j] = 0.2; 
	}
    }
  double yRaaHighPt[nCentBins][nRaaHighPt];
  double eylRaaHighPt[nCentBins][nRaaHighPt];
  double sylRaaHighPt[nCentBins][nRaaHighPt];
  double syhRaaHighPt[nCentBins][nRaaHighPt];
  data_in.open("/Users/admin/Work/STAR/analysis/Rootfiles/Publication.Jpsi.AuAu200.Raa.HighPt.txt");
  int counter = 0;
  while(!data_in.eof())
    {
      for(int i=0; i<6; i++)
	data_in >> data[i];

      int cent = counter/nRaaHighPt;
      int ptbin = counter%nRaaHighPt;
      xRaaHighPt[cent][ptbin] = atof(data[1]);
      yRaaHighPt[cent][ptbin] = atof(data[2]);
      eylRaaHighPt[cent][ptbin] = atof(data[3]);
      sylRaaHighPt[cent][ptbin] = atof(data[4]);
      syhRaaHighPt[cent][ptbin] = fabs(atof(data[5]));
      cout << xRaaHighPt[cent][ptbin] << endl;
      counter++;
      if(counter==nRaaHighPt*nCentBins) break;
    }
  data_in.close();
  
  // Jpsi Raa
  TCanvas *c = new TCanvas("AuAu200_Raa","AuAu200_Raa",1100,700);
  c->Divide(2,2);
  TH1F *hRaa = new TH1F("AuAu200_Raa",";p_{T} (GeV/c);R_{AA}",10,0,10);
  hRaa->GetYaxis()->SetRangeUser(0,1.8);
  ScaleHistoTitle(hRaa,0.06,1,0.05,0.06,1,0.05,62);

  TGraphAsymmErrors *gRaaLowPt[nCentBins];
  TGraphAsymmErrors *gRaaLowPtSys[nCentBins];
  TGraphAsymmErrors *gRaaHighPt[nCentBins];
  TGraphAsymmErrors *gRaaHighPtSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
      hRaa->Draw();
      
      gRaaLowPt[i] = new TGraphAsymmErrors(nAuAuLowPt, xRaaLowPt[i], yRaaLowPt[i], exlRaaLowPt[i], exlRaaLowPt[i], eylRaaLowPt[i], eylRaaLowPt[i]);
      gRaaLowPt[i]->SetName(Form("Jpsi_InvYield_Raa200_LowPt_cent%s",cent_Title[i]));
      gRaaLowPt[i]->SetMarkerStyle(20);
      gRaaLowPt[i]->SetMarkerSize(1.5);
      gRaaLowPt[i]->SetMarkerColor(2);
      gRaaLowPt[i]->SetLineColor(2);
      gRaaLowPt[i]->Draw("sames PE");

      gRaaLowPtSys[i] = new TGraphAsymmErrors(nAuAuLowPt, xRaaLowPt[i], yRaaLowPt[i], sxlRaaLowPt[i], sxlRaaLowPt[i], sylRaaLowPt[i], syhRaaLowPt[i]);
      gRaaLowPtSys[i]->SetName(Form("Jpsi_InvYield_Raa200_LowPt_systematics_cent%s",cent_Title[i]));
      gRaaLowPtSys[i]->SetMarkerStyle(20);
      gRaaLowPtSys[i]->SetMarkerSize(0);
      gRaaLowPtSys[i]->SetFillStyle(0);
      gRaaLowPtSys[i]->SetMarkerColor(2);
      gRaaLowPtSys[i]->SetLineColor(2);
      gRaaLowPtSys[i]->Draw("sameE5");

      gRaaHighPt[i] = new TGraphAsymmErrors(nRaaHighPt, xRaaHighPt[i], yRaaHighPt[i], exlRaaHighPt[i], exlRaaHighPt[i], eylRaaHighPt[i], eylRaaHighPt[i]);
      gRaaHighPt[i]->SetName(Form("Jpsi_InvYield_Raa200_HighPt_cent%s",cent_Title[i]));
      gRaaHighPt[i]->SetMarkerStyle(21);
      gRaaHighPt[i]->SetMarkerSize(1.5);
      gRaaHighPt[i]->SetMarkerColor(4);
      gRaaHighPt[i]->SetLineColor(4);
      gRaaHighPt[i]->Draw("sames PE");

      gRaaHighPtSys[i] = new TGraphAsymmErrors(nRaaHighPt, xRaaHighPt[i], yRaaHighPt[i], sxlRaaHighPt[i], sxlRaaHighPt[i], sylRaaHighPt[i], syhRaaHighPt[i]);
      gRaaHighPtSys[i]->SetName(Form("Jpsi_InvYield_Raa200_HighPt_systematics_cent%s",cent_Title[i]));
      gRaaHighPtSys[i]->SetMarkerStyle(21);
      gRaaHighPtSys[i]->SetMarkerSize(0);
      gRaaHighPtSys[i]->SetFillStyle(0);
      gRaaHighPtSys[i]->SetMarkerColor(4);
      gRaaHighPtSys[i]->SetLineColor(4);
      gRaaHighPtSys[i]->Draw("sameE5");

      TPaveText *t1 = GetPaveText(0.56,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[i]));
      t1->Draw();

      TLine *line = GetLine(0,1,10,1,1);
      line->Draw();
    }
  c->cd(1);
  leg = new TLegend(0.2,0.6,0.4,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Au+Au @ 200 GeV");
  leg->AddEntry(gRaaLowPt[0],"STAR: Run10","P");
  leg->AddEntry(gRaaHighPt[0],"STAR: Run10","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Raa.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Raa.png",run_type));
    }

  // combined data point for fit
  TCanvas *c = new TCanvas("AuAu200_Jpsi_fit","AuAu200_Jpsi_fit",1100,700);
  c->Divide(2,2);

  const int nbinsAuAu = 9;
  const double xbinsAuAu[nbinsAuAu+1] = {0,1,2,3,4,5,6,7,8,10};
  TH1F *hJpsiAuAu[nCentBins];
  TF1 *funcJpsiAuAu[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      hJpsiAuAu[i] = new TH1F(Form("Jpsi_InvYield_AuAu200_Combined_cent%s",cent_Title[i]),Form("Invariant yield of J/#psi (%s%%);p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",cent_Name[i]), nbinsAuAu, xbinsAuAu);
      for(int bin=1; bin<=nbinsAuAu; bin++)
	{
	  if(bin<=4)
	    {
	      hJpsiAuAu[i]->SetBinContent(bin, yAuAuLowPt[i][bin-1]);
	      hJpsiAuAu[i]->SetBinError(bin, eylAuAuLowPt[i][bin-1]);
	    }
	  else
	    {
	      hJpsiAuAu[i]->SetBinContent(bin,  yAuAuHighPt[i][bin-4]);
	      hJpsiAuAu[i]->SetBinError(bin, eylAuAuHighPt[i][bin-4]);
	    }
	}
      
      hJpsiAuAu[i]->SetMarkerStyle(21);
      hJpsiAuAu[i]->GetYaxis()->SetRangeUser(1e-10,1e-4);
      ScaleHistoTitle(hJpsiAuAu[i],0.06,1,0.05,0.06,1,0.05,62);
      c->cd(i+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);

      funcJpsiAuAu[i] = new TF1(Form("Fit_Jpsi_InvYield_AuAu200_Combined_cent%s",cent_Title[i]),"[0]*exp([1]*x+[2]*x*x+[3]*x*x*x)",0.1,10);
      if(i==3)funcJpsiAuAu[i]->SetParameter(0,1e-5);
      hJpsiAuAu[i]->Fit(funcJpsiAuAu[i],"IR0");
      hJpsiAuAu[i]->Draw();
      funcJpsiAuAu[i]->SetLineColor(2);
      funcJpsiAuAu[i]->Draw("sames");
    }


  // Jpsi yield
  TH1F *hTBW[4];
  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-10,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);

  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      c->cd(i+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.1);
      hAuAu->Draw();
      
      gAuAuLowPt[i] = new TGraphAsymmErrors(nAuAuLowPt, xAuAuLowPt[i], yAuAuLowPt[i], exlAuAuLowPt[i], exlAuAuLowPt[i], eylAuAuLowPt[i], eylAuAuLowPt[i]);
      gAuAuLowPt[i]->SetName(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[i]));
      gAuAuLowPt[i]->SetMarkerStyle(20);
      gAuAuLowPt[i]->SetMarkerSize(1.5);
      gAuAuLowPt[i]->SetMarkerColor(2);
      gAuAuLowPt[i]->SetLineColor(2);
      gAuAuLowPt[i]->Draw("sames PE");

      gAuAuLowPtSys[i] = new TGraphAsymmErrors(nAuAuLowPt, xAuAuLowPt[i], yAuAuLowPt[i], sxlAuAuLowPt[i], sxlAuAuLowPt[i], sylAuAuLowPt[i], syhAuAuLowPt[i]);
      gAuAuLowPtSys[i]->SetName(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[i]));
      gAuAuLowPtSys[i]->SetMarkerStyle(20);
      gAuAuLowPtSys[i]->SetMarkerSize(0);
      gAuAuLowPtSys[i]->SetFillStyle(0);
      gAuAuLowPtSys[i]->SetMarkerColor(2);
      gAuAuLowPtSys[i]->SetLineColor(2);
      gAuAuLowPtSys[i]->Draw("sameE5");

      gAuAuHighPt[i] = new TGraphAsymmErrors(nAuAuHighPt, xAuAuHighPt[i], yAuAuHighPt[i], exlAuAuHighPt[i], exlAuAuHighPt[i], eylAuAuHighPt[i], eylAuAuHighPt[i]);
      gAuAuHighPt[i]->SetName(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[i]));
      gAuAuHighPt[i]->SetMarkerStyle(21);
      gAuAuHighPt[i]->SetMarkerSize(1.5);
      gAuAuHighPt[i]->SetMarkerColor(4);
      gAuAuHighPt[i]->SetLineColor(4);
      gAuAuHighPt[i]->Draw("sames PE");

      gAuAuHighPtSys[i] = new TGraphAsymmErrors(nAuAuHighPt, xAuAuHighPt[i], yAuAuHighPt[i], sxlAuAuHighPt[i], sxlAuAuHighPt[i], sylAuAuHighPt[i], syhAuAuHighPt[i]);
      gAuAuHighPtSys[i]->SetName(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[i]));
      gAuAuHighPtSys[i]->SetMarkerStyle(21);
      gAuAuHighPtSys[i]->SetMarkerSize(0);
      gAuAuHighPtSys[i]->SetFillStyle(0);
      gAuAuHighPtSys[i]->SetMarkerColor(4);
      gAuAuHighPtSys[i]->SetLineColor(4);
      gAuAuHighPtSys[i]->Draw("sameE5");

      TPaveText *t1 = GetPaveText(0.56,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[i]));
      t1->Draw();

      hTBW[i] = (TH1F*)fout->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[i]));
      hTBW[i]->SetLineColor(1);
      hTBW[i]->SetLineStyle(2);
      hTBW[i]->SetLineWidth(2);
      hTBW[i]->Draw("sames");
    }
  c->cd(1);
  leg = new TLegend(0.2,0.2,0.4,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->SetHeader("Au+Au @ 200 GeV");
  leg->AddEntry(gAuAuLowPt[0],"STAR: Run10","P");
  leg->AddEntry(gAuAuHighPt[0],"STAR: Run10","P");
  leg->AddEntry(hTBW[0],"TBW (#beta=0)","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Yield.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/2015QM/Pub_Jpsi_Yield.png",run_type));
    }

  if(saveHisto)
    {
      fout->cd();

      gHighPtPP->Write("",TObject::kOverwrite);
      gHighPtPPSys->Write("",TObject::kOverwrite);
      hLowPtPP->Write("",TObject::kOverwrite);
      hJpsiPP->Write("",TObject::kOverwrite);

      for(int i=0; i<nCentBins; i++)
	{
	  gAuAuLowPt[i]->Write("",TObject::kOverwrite);
	  gAuAuLowPtSys[i]->Write("",TObject::kOverwrite);
	  gAuAuHighPt[i]->Write("",TObject::kOverwrite);
	  gAuAuHighPtSys[i]->Write("",TObject::kOverwrite);
	  hTBW[i]->Write("",TObject::kOverwrite);
	  hJpsiAuAu[i]->Write("",TObject::kOverwrite);
	  funcJpsiAuAu[i]->Write("",TObject::kOverwrite);
	}

      for(int i=0; i<nCentBins; i++)
	{
	  gRaaLowPt[i]->Write("",TObject::kOverwrite);
	  gRaaLowPtSys[i]->Write("",TObject::kOverwrite);
	  gRaaHighPt[i]->Write("",TObject::kOverwrite);
	  gRaaHighPtSys[i]->Write("",TObject::kOverwrite);
	}
    }
}

