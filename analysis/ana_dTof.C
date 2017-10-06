const int year = YEAR;
TFile *f = 0x0;

//================================================
void ana_Dtof()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  
  //scanDtofCut();
  TagAndProbe();
  //JpsiMuon();
}

//================================================
void TagAndProbe(const int savePlot = 1, const int saveHisto = 1)
{
  // const int year = 2014;
  // const char* run_type = "Run14_AuAu200";
  // const double dtofCut = 0.75;

  const int year = 2015;
  const char* run_type = "Run15_pp200";
  const double dtofCut = 1;

  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};

  if(year==2014) f = TFile::Open("Output/Run14_AuAu200.jpsi.root","read");
  else           f = TFile::Open(Form("Output/%s.TagAndProbe.root",run_type),"read");
  THnSparseF *hnDtof = (THnSparseF*)f->Get("mhTaP_dtof_di_mu");
  hnDtof->GetAxis(3)->SetRange(1,1);

  TH2F *hInvMassVsPt[2];
  hInvMassVsPt[0] = (TH2F*)hnDtof->Projection(0,1);
  hInvMassVsPt[0]->SetName("mhTaP_dtof_di_mu_All");
  if(year==2014) hnDtof->GetAxis(2)->SetRangeUser(-100, dtofCut-1e-4);
  if(year==2015) hnDtof->GetAxis(2)->SetRangeUser(-1*dtofCut+1e-4, dtofCut-1e-4);
  hInvMassVsPt[1] = (TH2F*)hnDtof->Projection(0,1);
  hInvMassVsPt[1]->SetName("mhTaP_dtof_di_mu_Acc");

  // Get the mean pT
  TH1F *hMeanPt = new TH1F("hMeanPt", "Mean p_{T} of muons", nbins, 0, nbins);
  TH1F *hMuonPt = (TH1F*)hInvMassVsPt[1]->ProjectionX("hMuonPt");
  for(int bin=1; bin<=nbins; bin++)
    {
      hMuonPt->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
      hMeanPt->SetBinContent(bin, hMuonPt->GetMean());
      hMeanPt->SetBinError(bin, hMuonPt->GetMeanError());
    }
  
  //==============================================
  // Fit to invariant mass
  //==============================================
  TString hTitle[2] = {"All", Form("#Deltatof<%2.2fns",dtofCut)};
  if(year==2015) hTitle[1] = Form("|#Deltatof|<%2.0fns",dtofCut);
  const TString hName[2] = {"All","Acc"};
  TH1F *hInvMass[2][nbins];
  TF1 *funcInvMass[2][nbins];
  TH1F *hJpsiCounts[2];
  double xmin, xmax;
  for(int i=0; i<2; i++)
    {
      hJpsiCounts[i] = new TH1F(Form("hJpsiCounts_%d",i), "# of J/psi", nbins, xbins);
      TCanvas *c = new TCanvas(Form("FitInvMass_%d",i), Form("FitInvMass_%d",i), 1100, 700);
      c->Divide(3,2);
      for(int j=0; j<nbins; j++)
	{
	  int low_bin = hInvMassVsPt[i]->GetXaxis()->FindFixBin(xbins[j]+1e-4);
	  int up_bin  = hInvMassVsPt[i]->GetXaxis()->FindFixBin(xbins[j+1]-1e-4);
	  hInvMass[i][j] = (TH1F*)hInvMassVsPt[i]->ProjectionY(Form("hInvMass_bin%d_%d",j+1,i),low_bin, up_bin);
	  hInvMass[i][j]->SetTitle("");
	  hInvMass[i][j]->Sumw2();
	  if(j<3) hInvMass[i][j]->Rebin(2);
	  else    hInvMass[i][j]->Rebin(5);
	  funcInvMass[i][j] = new TF1(Form("func_%s",hInvMass[i][j]->GetName()),"gausn(0)+pol3(3)",2.5,3.6);
	  funcInvMass[i][j]->SetParameter(0, 100);
	  funcInvMass[i][j]->FixParameter(1, 3.096);
	  funcInvMass[i][j]->SetParameter(2, 0.5);
	  if(j==0) 
	    {
	      funcInvMass[i][j]->SetParLimits(2, 0, 0.05);
	    }
	  if(j==2) 
	    {
	      funcInvMass[i][j]->SetParameter(0, 100);
	      funcInvMass[i][j]->SetParameter(2, 0.05);
	    }

	  if(year==2014)
	    {
	      if(j==0) funcInvMass[i][j]->SetRange(2.85, 3.4);
	      if(j==2) funcInvMass[i][j]->SetRange(2.5, 3.3);
	    }

	  if(i==1)
	    {
	      funcInvMass[i][j]->SetParameters(funcInvMass[0][j]->GetParameters());
	      funcInvMass[i][j]->FixParameter(2, funcInvMass[0][j]->GetParameter(2));
	    }

	  if(j<2) hInvMass[i][j]->GetXaxis()->SetRangeUser(2.8, 3.5);
	  hInvMass[i][j]->Fit(funcInvMass[i][j], "IR0S");
	  c->cd(j+1);
	  hInvMass[i][j]->SetMarkerStyle(24);
	  hInvMass[i][j]->Draw();
	  funcInvMass[i][j]->SetLineStyle(2);
	  funcInvMass[i][j]->SetLineColor(2);
	  funcInvMass[i][j]->GetRange(xmin, xmax);
	  TF1 *tmpBkg = new TF1(Form("TmpBkg_%d_%d",i,j), "pol3", xmin, xmax);
	  for(int ipar=0; ipar<4; ipar++)
	    {
	      tmpBkg->SetParameter(ipar, funcInvMass[i][j]->GetParameter(ipar+3));
	    }
	  tmpBkg->SetLineStyle(2);
	  tmpBkg->SetLineColor(4);
	  tmpBkg->Draw("sames");
	  funcInvMass[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T}^{#mu} < %1.1f GeV/c (%s)",xbins[j],xbins[j+1],hTitle[i].Data()), 0.05);
	  t1->Draw();
	  double bin_width = hInvMass[i][j]->GetBinWidth(1);
	  hJpsiCounts[i]->SetBinContent(j+1, fabs(funcInvMass[i][j]->GetParameter(0))/bin_width);
	  hJpsiCounts[i]->SetBinError(j+1, fabs(funcInvMass[i][j]->GetParError(0))/bin_width);
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_FitInvMass_%s.pdf",run_type,hName[i].Data()));
    }
  TList *list = new TList;
  for(int i=0; i<2; i++)
    {
      list->Add(hJpsiCounts[i]);
    }
  TCanvas *c = drawHistos(list,"JpsiCounts",Form("%s: # of Jpsi vs. muon p_{T};p_{T}^{#mu} (GeV/c);counts",run_type),kFALSE,0,30,kFALSE,0.5,1.1,kTRUE,kTRUE,hTitle,kTRUE,"",0.5,0.7,0.6,0.85,true);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_JpsiCounts.pdf",run_type));

  //==============================================
  // extract efficiency
  //==============================================
  TH1F *hJpsiCountsClone[2];
  for(int i=0; i<2; i++)
    {
      hJpsiCountsClone[i] = (TH1F*)hJpsiCounts[i]->Clone(Form("%s_clone",hJpsiCounts[i]->GetName()));
    }
  for(int bin=1; bin<nbins; bin++)
    {
      if(hJpsiCountsClone[1]->GetBinContent(bin)>hJpsiCountsClone[0]->GetBinContent(bin))
	{
	  hJpsiCountsClone[1]->SetBinContent(bin, hJpsiCountsClone[0]->GetBinContent(bin));
	}
    }
  
  TH1F *hRatio = hJpsiCounts[1]->Clone(Form("hRatio"));
  hRatio->Divide(hJpsiCounts[0]);
  for(int bin=1; bin<nbins; bin++)
    {
      if(hRatio->GetBinContent(bin)>1)
	{
	  double err = hJpsiCounts[0]->GetBinError(bin)/hJpsiCounts[0]->GetBinContent(bin);
	  hRatio->SetBinError(bin, err*hRatio->GetBinContent(bin));
	}
    }

  TGraphAsymmErrors *gEff = new TGraphAsymmErrors(hJpsiCountsClone[1], hJpsiCountsClone[0],"cl=0.683 b(1,1) mode");
  TGraphErrors *gDtofEff = new TGraphErrors(nbins);
  gDtofEff->SetName(Form("TagAndProbe_Muon_Dtof%2.2fEff",dtofCut));
  double x, y;
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gEff->GetPoint(ipoint,x,y);
      double err = gEff->GetErrorYhigh(ipoint) > gEff->GetErrorYlow(ipoint) ? gEff->GetErrorYhigh(ipoint) : gEff->GetErrorYlow(ipoint);
      if(y>=1)
	{
	  y = hRatio->GetBinContent(ipoint+1);
	  err = hRatio->GetBinError(ipoint+1);
	}
      gDtofEff->SetPoint(ipoint, hMeanPt->GetBinContent(ipoint+1), y);
      gDtofEff->SetPointError(ipoint, hMeanPt->GetBinError(ipoint+1), err);
    }
  TGraphErrors *gfit = (TGraphErrors*)gDtofEff->Clone(Form("Fit_%s",gDtofEff->GetName()));
  TF1 *funcDtofEff = new TF1(Form("%s_FitFunc",gDtofEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  if(year==2015) funcDtofEff->SetParameter(0, 0.985);
  gfit->Fit(funcDtofEff,"R0");
  funcDtofEff->SetLineColor(2);
  funcDtofEff->SetLineStyle(2);

  TH1F *hplot = new TH1F("hplot",Form(";p_{T,#mu} (GeV/c);Eff."), 100, 0, 10);
  //==============================================
  // systematic uncertainty
  //==============================================
  // i)   statistical error: randomization method
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  rndm->SetSeed(0);
  const int nexpr = 400;
  TGraphErrors *gEffRndm[nexpr];
  double *ex = gDtofEff->GetEX();
  double *ey = gDtofEff->GetEY();
  for(int j = 0; j < nexpr; j++)
    {
      int npoint = gDtofEff->GetN();
      gEffRndm[j] = new TGraphErrors(npoint);
      gEffRndm[j]->SetName(Form("%s_Rndm%d",gDtofEff->GetName(),j));
      for(int ipoint=0; ipoint<npoint; ipoint++)
	{
	  gDtofEff->GetPoint(ipoint,x,y);
	  y = rndm->Gaus(y, ey[ipoint]);
	  gEffRndm[j]->SetPoint(ipoint,x,y);
	  gEffRndm[j]->SetPointError(ipoint,ex[ipoint],ey[ipoint]);	 
	}
    }

  hplot->GetYaxis()->SetRangeUser(0.6,1.1);
  hplot->GetXaxis()->SetRangeUser(1.3,7);
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1.4,1.5,1.75,2,2.25,2.5,3,4,5,6,8};
  TH1F *hSpread[nbins_rndm];
  TGraphErrors *gLimits[3];
  TF1 *funcLimits[3];
  const char *limit_name[3] = {"center","down","up"};

  c = new TCanvas(Form("cEff"),Form("cEff"),800,600);
  hplot->DrawCopy();
  TPaveText *t1 = GetTitleText(Form("Estimation of uncertainty for #Deltatof < %2.2f ns cut",dtofCut));
  t1->Draw();
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%s",run_type));

  for(int l=0; l<nbins_rndm; l++)
    {
      hSpread[l] = new TH1F(Form("hSpread_%d",l),Form("hSpread_%d",l),600,0,1.2);
    }
  for(int j = 0; j < nexpr; j++)
    {
      TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
      funcTmp->SetParameters(funcDtofEff->GetParameters());
      gEffRndm[j]->Fit(funcTmp,"R0Q");
      c->cd();
      funcTmp->SetLineStyle(2);
      funcTmp->Draw("sames");
      if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
      for(int l=0; l<nbins_rndm; l++)
	{
	  hSpread[l]->Fill(funcTmp->Eval(xbins_rndm[l]));
	}
    }

  TCanvas *cFit = new TCanvas(Form("cEffFit"),Form("cEffFit"),1200,800);
  cFit->Divide(4,3);
  for(int t=0; t<3; t++) gLimits[t] = new TGraphErrors(nbins_rndm);
  for(int l=0; l<nbins_rndm; l++)
    {
      cFit->cd(l+1);
      hSpread[l]->Rebin(2);
      if(l<2) hSpread[l]->GetXaxis()->SetRangeUser(0.3,1);
      else if(l<4) hSpread[l]->GetXaxis()->SetRangeUser(0.65,1);
      else hSpread[l]->GetXaxis()->SetRangeUser(0.75,1);
      hSpread[l]->SetTitle(";Efficiency;");
      hSpread[l]->SetMaximum(1.3*hSpread[l]->GetMaximum());
      hSpread[l]->Draw();
      TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[l]),0.06);
      t1->Draw();
      TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
      funcTmp2->SetParameter(0,hSpread[l]->GetMaximum());
      funcTmp2->SetParameter(1,hSpread[l]->GetMean());
      funcTmp2->SetParameter(2,hSpread[l]->GetRMS());
      hSpread[l]->Fit(funcTmp2,"R0Q");
      funcTmp2->SetLineColor(4);
      funcTmp2->SetLineStyle(2);
      funcTmp2->Draw("sames");

      double pt    = xbins_rndm[l];
      double mean  = hSpread[l]->GetMean();
      double width = hSpread[l]->GetRMS();
      double error = hSpread[l]->GetRMSError();
      gLimits[0]->SetPoint(l,pt,mean);
      gLimits[0]->SetPointError(l,0,error);
      gLimits[1]->SetPoint(l,pt,mean-width);
      gLimits[1]->SetPointError(l,0,error);
      gLimits[2]->SetPoint(l,pt,mean+width);
      gLimits[2]->SetPointError(l,0,error);

      TPaveText *t1 = GetPaveText(0.6,0.7,0.7,0.85,0.06);
      t1->SetTextAlign(11);
      t1->AddText(Form("RMS=%4.3f",width));
      t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcTmp2->GetParameter(2))));
      t1->Draw();
    }
  cFit->cd(12);
  TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
  t1->AddText(run_type);
  t1->Draw();
  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f_RndmFit.pdf",run_type,dtofCut));
  c->cd();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f_RndmCurve.pdf",run_type,dtofCut));
  for(int t=0; t<3; t++)
    {
      gLimits[t]->SetMarkerStyle(20);
      gLimits[t]->SetMarkerColor(2+2*t);
      gLimits[t]->SetLineColor(2+2*t);
      gLimits[t]->Draw("samesP");

      funcLimits[t] = new TF1(Form("%s_Sys%s",gDtofEff->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
      if(year==2016) funcLimits[t]->SetRange(1.5,7);
      funcLimits[t]->SetParameters(0.95,2,0.2);
      gLimits[t]->Fit(funcLimits[t],"R0Q");
      funcLimits[t]->SetLineColor(gLimits[t]->GetMarkerColor());
      if(t>0) funcLimits[t]->SetLineStyle(2);
      funcLimits[t]->DrawCopy("sames");
    }
  leg->AddEntry(gLimits[0],"Central value","PL");
  leg->AddEntry(gLimits[1],"Lower limit","PL");
  leg->AddEntry(gLimits[2],"Upper limit","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f_RndmLimits.pdf",run_type,dtofCut));

  //==============================================
  // Final efficiency
  //==============================================
  hplot->GetXaxis()->SetRangeUser(1.3, 7);
  if(year==2014) 
    {
      hplot->GetYaxis()->SetRangeUser(0.6, 1.05);
      c = draw1D(hplot,Form("%s: efficiency for #Deltatof < %2.2f ns cut",run_type,dtofCut));
    }
  if(year==2015)
    {
      hplot->GetYaxis()->SetRangeUser(0.8, 1.05);
      c = draw1D(hplot,Form("%s: efficiency for |#Deltatof| < %2.0f ns cut",run_type,dtofCut));
    }
  gPad->SetGrid(0,1);
  gDtofEff->SetMarkerStyle(21);
  gDtofEff->SetMarkerSize(1.5);
  gDtofEff->Draw("samesPEZ");
  if(year==2014) funcDtofEff->Draw("sames");
  if(year==2015) funcLimits[0]->Draw("sames");
  TLegend *leg = new TLegend(0.2,0.3,0.35,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader("Tag-and-probe");
  leg->AddEntry(gDtofEff,"Data","p");
  if(year==2014) leg->AddEntry(funcDtofEff,"Fit to data","L");
  if(year==2015) leg->AddEntry(funcLimits[0],"Mean of randomization","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f.pdf",run_type,dtofCut));
  funcLimits[1]->DrawCopy("sames");
  funcLimits[2]->SetLineColor(4);
  funcLimits[2]->DrawCopy("sames");
  TLegend *leg2 = new TLegend(0.2,0.25,0.35,0.3);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(funcLimits[1],"Systematic uncertainty","L");
  leg2->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f_WithSys.pdf",run_type,dtofCut));

  //==============================================
  // Compare with Jpsi muon method
  //==============================================
  TFile *fout = 0x0;
  if(year==2014)
    {
      if(saveHisto) fout = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"update");
      else          fout = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"read");
      TGraphAsymmErrors *gJpsiMuon[2];
      gJpsiMuon[0] = (TGraphAsymmErrors*)fout->Get(Form("JpsiMuon_DtofEff%2.2f_Fitting",dtofCut));
      gJpsiMuon[1] = (TGraphAsymmErrors*)fout->Get(Form("JpsiMuon_DtofEff%2.2f_BinCount",dtofCut));
      for(int i=0; i<2; i++)
	{
	  gJpsiMuon[i]->SetMarkerStyle(20+i*4);
	  gJpsiMuon[i]->SetMarkerColor(kGreen+2);
	  gJpsiMuon[i]->SetLineColor(kGreen+2);
	  gJpsiMuon[i]->SetMarkerSize(1.5);
	  gJpsiMuon[i]->Draw("samesPEZ");
	}
      TLegend *leg3 = new TLegend(0.55,0.3,0.7,0.45);
      leg3->SetBorderSize(0);
      leg3->SetFillColor(0);
      leg3->SetTextSize(0.04);
      leg3->SetHeader("Jpsi muon method");
      leg3->AddEntry(gJpsiMuon[0],"Fitting","P");
      leg3->AddEntry(gJpsiMuon[1],"Bin-counting","P");
      leg3->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/TagAndProbe_DtofEff%2.2f_CompToJpsiMuon.pdf",run_type,dtofCut));
    }


  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      if(!fout) fout = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"update");
      gDtofEff->Write("",TObject::kOverwrite);
      funcDtofEff->Write("",TObject::kOverwrite);
      funcLimits[0]->Write("",TObject::kOverwrite);
      funcLimits[1]->Write("",TObject::kOverwrite);
      funcLimits[2]->Write("",TObject::kOverwrite);
      fout->Close();
    }
}

//================================================
void JpsiMuon(const Int_t savePlot = 1, const int saveHisto = 1)
{
  const int nbins = 6;
  const double xbins[nbins+1] = {1.3,1.5,2.0,2.5,3.0,5.0,10.0};
  const double minimum = -3;
  const double maximum = 3;
  const double min_mass[3] = {3.0, 2.6, 3.3};
  const double max_mass[3] = {3.2, 2.9, 3.6};
  const double dtof_cut_min = -3;
  const double dtof_cut_max = 0.75;

  TFile *fdata = 0x0;
  if(year==2014) fdata = TFile::Open("./output/Run14_AuAu200.JpsiMuon.root","read");
  if(year==2016) fdata = TFile::Open("./output/Pico.Run16.AuAu200.JpsiMuon.root","read");

  //==============================================
  // compare invariant mass
  //==============================================
  const int nbinsJpsi = 5;
  const double xbinsJpsi[nbinsJpsi+1] = {0,1.0,2.0,3.0,5.0,10.0};
  TH2F *hInvMassVsPtUL = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_Dtof_UL_di_mu"));
  TH2F *hInvMassVsPtLS = (TH2F*)fdata->Get(Form("mhJpsiMassVsPt_Dtof_LS_di_mu"));

  TCanvas *c = new TCanvas("InvMass_Dtof","InvMass_Dtof",1100,700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.2,0.4,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetHeader(run_type);
  for(int bin=1; bin<=nbinsJpsi; bin++)
    {
      int start_bin = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin-1]+1e-4);
      int end_bin   = hInvMassVsPtUL->GetXaxis()->FindBin(xbinsJpsi[bin]-1e-4);
      TH1F *hInvMassUL = (TH1F*)hInvMassVsPtUL->ProjectionY(Form("Data_Dtof_InvMassUL_bin%d",bin),start_bin,end_bin);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->SetMarkerStyle(20);
      hInvMassUL->Rebin(5);
      hInvMassUL->SetMaximum(1.5*hInvMassUL->GetMaximum());
      hInvMassUL->SetMinimum(0);
      hInvMassUL->GetXaxis()->SetRangeUser(2.2,4);
      hInvMassUL->SetTitle("");
      if(bin==1) leg->AddEntry(hInvMassUL,"Unlike-sign","PL");

      TH1F *hInvMassLS = (TH1F*)hInvMassVsPtLS->ProjectionY(Form("Data_Dtof_InvMassLS_bin%d",bin),start_bin,end_bin);
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
	  TH1F *htmp = new TH1F(Form("%s_tmp%d",hInvMassUL->GetName(),itmp),"",int((max_mass[itmp]-min_mass[itmp])/binwidth+0.5),min_mass[itmp],max_mass[itmp]);
	  for(int ibin=1; ibin<=htmp->GetNbinsX(); ibin++)
	    {
	      int jbin = hInvMassUL->FindFixBin(htmp->GetBinCenter(ibin));
	      htmp->SetBinContent(ibin,hInvMassUL->GetBinContent(jbin));
	    }
	  if(itmp==0) htmp->SetFillColor(7);
	  else htmp->SetFillColor(5);
	  htmp->Draw("sames");
	  if(bin==1 && itmp==0) leg->AddEntry(htmp,"Signal","F");
	  if(bin==1 && itmp==1) leg->AddEntry(htmp,"Side-band","F");
	}
      hInvMassLS->Draw("samesP");
      TPaveText *t1 = GetTitleText(Form("J/#psi: %1.1f < p_{T} < %1.1f GeV/c",xbinsJpsi[bin-1],xbinsJpsi[bin]),0.05);
      t1->Draw();
      int low_bin  = hInvMassUL->FindFixBin(min_mass[0]+1e-4);
      int high_bin = hInvMassUL->FindFixBin(max_mass[0]-1e-4);
      double all = hInvMassUL->Integral(low_bin, high_bin);
      double bkg = hInvMassLS->Integral(low_bin, high_bin);
      TPaveText *t1 = GetPaveText(0.15,0.5,0.72,0.85,0.05);
      t1->SetTextAlign(11);
      t1->AddText(Form("S/B = 1:%4.2f",bkg/(all-bkg)));
      t1->Draw();
    }
  c->cd(6);
  leg->Draw();
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_InvMassLSvsUL.pdf",run_type));

  //==============================================
  // single muon UL vs LS
  //==============================================
  TH2F *hJpsiMassVsMuonPt[2];
  TH1F *hMuonPtSide[2];
  hJpsiMassVsMuonPt[0] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_Dtof_UL_di_mu"));
  hJpsiMassVsMuonPt[1] = (TH2F*)fdata->Get(Form("mhJpsiMassVsMuonPt_Dtof_LS_di_mu"));
  for(int i=0; i<2; i++)
    {
      // Get the muon distribution in the side-band region
      hJpsiMassVsMuonPt[i]->Sumw2();
      hMuonPtSide[i] = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_Dtof_SideBand_%d",i));
      hMuonPtSide[i]->Reset();
      for(int j=0; j<2; j++)
	{
	  int low_bin = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(min_mass[j+1]);
	  int up_bin  = hJpsiMassVsMuonPt[i]->GetYaxis()->FindFixBin(max_mass[j+1]);
	  htmp = (TH1F*)hJpsiMassVsMuonPt[i]->ProjectionX(Form("hMuonPt_MtdVpdTacDiff_SideBand%d_%d",j,i),low_bin,up_bin);
	  hMuonPtSide[i]->Add(htmp);
	}
    }

  TH1F *hScaleFactor = new TH1F("hScaleFactor","Scale factor for unlike-sign background (#Deltatof cut);p_{T} (GeV/c)",nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      int start_bin = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hMuonPtSide[0]->GetXaxis()->FindBin(xbins[bin]-1e-4);
      double ul = hMuonPtSide[0]->Integral(start_bin,end_bin);
      double ls = hMuonPtSide[1]->Integral(start_bin,end_bin);
      double scale = ul/ls;
      double error = scale * sqrt(1/ul+1/ls);
      hScaleFactor->SetBinContent(bin, scale);
      hScaleFactor->SetBinError(bin, error);
    }
  hScaleFactor->SetMarkerStyle(21);
  hScaleFactor->GetYaxis()->SetRangeUser(0.95,1.15);
  c = draw1D(hScaleFactor);
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_ScaleFactor.pdf",run_type));


  TH2F *hDataUL = (TH2F*)fdata->Get("mhJpsiMuon_Dtof_UL_di_mu");
  hDataUL->Sumw2();
  TH2F *hDataLS = (TH2F*)fdata->Get("mhJpsiMuon_Dtof_LS_di_mu");
  hDataLS->Sumw2();
  TH2F *hDataDisVsPt = (TH2F*)hDataUL->Clone("JpsiMuonDtof");
  hDataDisVsPt->Add(hDataLS, -1);
  TH1F *hUL[nbins];
  TH1F *hLS[nbins];
  TH1F *hMuon[nbins];
  TH1F *hMuonFineBin[nbins];
  double mean_pt[nbins];
  double mean_pt_err[nbins];
  const int rebin = 10;
  c = new TCanvas(Form("MtdVpdTacDiff_UL_vs_LS"),Form("MtdVpdTacDiff_UL_vs_LS"),1100,700);
  c->Divide(3,2);

  // calculate mean pt of each bin
  int bin1 = hDataDisVsPt->GetYaxis()->FindBin(minimum);
  int bin2 = hDataDisVsPt->GetYaxis()->FindBin(maximum);
  htmp = (TH1F*)hDataDisVsPt->ProjectionX(Form("%s_tmp",hDataDisVsPt->GetName()),bin1,bin2);
  for(int bin=1; bin<=nbins; bin++)
    {
      htmp->GetXaxis()->SetRangeUser(xbins[bin-1]+1e-4, xbins[bin]-1e-4);
      mean_pt[bin-1] = htmp->GetMean();
      mean_pt_err[bin-1] = htmp->GetMeanError();

      int start_bin = hDataUL->GetXaxis()->FindBin(xbins[bin-1]+1e-4);
      int end_bin   = hDataUL->GetXaxis()->FindBin(xbins[bin]-1e-4);

      hUL[bin-1] = (TH1F*)hDataUL->ProjectionY(Form("%s_DataMtdVpdTacDiff_UL_bin%d",run_type,bin),start_bin,end_bin);
      hUL[bin-1]->SetMarkerStyle(20);
      hUL[bin-1]->SetMarkerStyle(20);

      hLS[bin-1] = (TH1F*)hDataLS->ProjectionY(Form("%s_DataMtdVpdTacDiff_LS_bin%d",run_type,bin),start_bin,end_bin);
      hLS[bin-1]->SetMarkerStyle(24);
      hLS[bin-1]->SetMarkerColor(2);
      hLS[bin-1]->SetLineColor(2);

      double scale = hScaleFactor->GetBinContent(bin);
      hMuonFineBin[bin-1] = (TH1F*)hUL[bin-1]->Clone(Form("DataJpsiMuon_Dtof_bin%d",bin));
      hMuonFineBin[bin-1]->Add(hLS[bin-1],-1.0*scale);

      // rebin
      hMuon[bin-1] = (TH1F*)hMuonFineBin[bin-1]->Clone(Form("DataJpsiMuon_Dtof_bin%d_Rebin",bin));
      hMuon[bin-1]->Rebin(rebin);
      
      c->cd(bin);
      hUL[bin-1]->Rebin(rebin);
      hUL[bin-1]->SetMaximum(1.5*hUL[bin-1]->GetMaximum());
      hUL[bin-1]->GetXaxis()->SetRangeUser(minimum,maximum);
      hUL[bin-1]->SetTitle("");
      hUL[bin-1]->Draw("P");
      hLS[bin-1]->Rebin(rebin);
      hLS[bin-1]->Draw("samesP");
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
      t1 = GetPaveText(0.6,0.8,0.75,0.8,0.05);
      t1->AddText(Form("Scale = %4.3f",scale));
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.15,0.63,0.35,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(run_type);
  leg->AddEntry(hUL[0],"Unlike-sign","PL");
  leg->AddEntry(hLS[0],"Like-sign","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_MuonULvsLS.pdf",run_type));

  //==============================================
  // Fit data distributions
  //==============================================
  TF1 *funcFitData[nbins];
  TFitResultPtr ptr[nbins];
  TH1F *hFitDataMean = new TH1F(Form("DataJpsiMuon_Dtof_FitMean"),Form("%s: mean of #Deltatof;p_{T} (GeV/c)",run_type),nbins,xbins);
  TH1F *hFitDataSigma = new TH1F(Form("DataJpsiMuon_Dtof_FitSigma"),Form("%s: sigma of #Deltatof;p_{T} (GeV/c)",run_type),nbins,xbins);
  c = new TCanvas(Form("Fit_Dtof"),Form("Fit_Dtof"),1100,700);
  c->Divide(3,2);
  for(int bin=1; bin<=nbins; bin++)
    {
      TH1F *hFit = (TH1F*)hMuon[bin-1]->Clone(Form("Fit_%s",hMuon[bin-1]->GetName()));

      funcFitData[bin-1] = new TF1(Form("DataJpsiMuon_Dtof_bin%d_Fit",bin),CrystalBall,-3,3,5);
      funcFitData[bin-1]->SetParameters(-1.5, -0.05, 0.2, 50, hFit->GetMaximum());
      // funcFitData[bin-1] = new TF1(Form("func_%d_%d",i,bin),"gaus",-3,0.5);
      // funcFitData[bin-1]->SetParameters(hFit->GetMaximum(),-1,0.2);
      ptr[bin-1] = hFit->Fit(funcFitData[bin-1],"IR0QS");
      ptr[bin-1]->SetName(Form("DataJpsiMuon_Dtof_bin%d_FitResult",bin));

      hFitDataMean->SetBinContent(bin,funcFitData[bin-1]->GetParameter(1));
      hFitDataMean->SetBinError(bin,funcFitData[bin-1]->GetParError(1));
      hFitDataSigma->SetBinContent(bin,funcFitData[bin-1]->GetParameter(2));
      hFitDataSigma->SetBinError(bin,funcFitData[bin-1]->GetParError(2));

      c->cd(bin);
      hFit->SetMaximum(1.5*hFit->GetMaximum());
      hFit->SetTitle("");
      hFit->GetXaxis()->SetRangeUser(minimum,maximum);
      hFit->Draw();
      funcFitData[bin-1]->SetLineColor(4);
      funcFitData[bin-1]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("J/#psi #mu: %1.1f < p_{T} < %1.1f",xbins[bin-1],xbins[bin]),0.06);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_Fit.pdf",run_type));

  hFitDataMean->SetMarkerStyle(21);
  hFitDataMean->GetYaxis()->SetRangeUser(-1,1);
  c = draw1D(hFitDataMean,"Mean of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);Mean");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_FitMean.pdf",run_type));

  hFitDataSigma->SetMarkerStyle(21);
  hFitDataSigma->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hFitDataSigma,"Width of #Deltatof distribution for J/#psi muons;p_{T} (GeV/c);#sigma");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_Dtof_FitSigma.pdf",run_type));


  //==============================================
  // extract efficiency
  //==============================================
  // bin counting method for data
  double x,y;
  TH1F *hBase = new TH1F("hBase","hBase",nbins,xbins);
  TH1F *hMatch = new TH1F("hMatch","hMatch",nbins,xbins);
  TH1F *hMatch2 = new TH1F("hMatch2","hMatch2",nbins,xbins);
  for(int bin=1; bin<=nbins; bin++)
    {
      double all_err;
      int low_bin = hMuonFineBin[bin-1]->FindFixBin(minimum+1e-4);
      int high_bin = hMuonFineBin[bin-1]->FindFixBin(maximum-1e-4);
      double all = hMuonFineBin[bin-1]->IntegralAndError(low_bin, high_bin,all_err);

      double acc_err;
      low_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_min+1e-4);
      high_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_max-1e-4);
      double acc = hMuonFineBin[bin-1]->IntegralAndError(low_bin, high_bin,acc_err);
      hBase->SetBinContent(bin,all);
      hBase->SetBinError(bin,all_err);
      hMatch->SetBinContent(bin,acc);
      hMatch->SetBinError(bin,acc_err);
      if(acc>all) { acc = all; acc_err = all_err; }
      hMatch2->SetBinContent(bin,acc);
      hMatch2->SetBinError(bin,acc_err);
    }

  // recalculate statistical errors
  hMatch->Divide(hBase);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int nRndm = 1e3;
  for(int bin=1; bin<=nbins; bin++)
    {
      if(hMatch->GetBinContent(bin)<0.98) continue;
      TH1F *hRndmDis = new TH1F(Form("RndmDis_Bin%d",bin),"",50,0.8,1.2);
      double all, all_err, acc, acc_err;
      int all_low_bin = hMuonFineBin[bin-1]->FindFixBin(minimum+1e-4);
      int all_high_bin = hMuonFineBin[bin-1]->FindFixBin(maximum-1e-4);
      int acc_low_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_min+1e-4);
      int acc_high_bin = hMuonFineBin[bin-1]->FindFixBin(dtof_cut_max-1e-4);
      int nbinsx = hMuonFineBin[bin-1]->GetNbinsX();
      for(int i=0; i<nRndm; i++)
	{
	  TH1F *htmp = (TH1F*)hMuonFineBin[bin-1]->Clone(Form("%s_rndm%d", hMuonFineBin[bin-1]->GetName(), i));
	  htmp->Reset();
	  for(int ibin=1; ibin<=nbinsx; ibin++)
	    {
	      double mean = hMuonFineBin[bin-1]->GetBinContent(ibin);
	      double error = hMuonFineBin[bin-1]->GetBinError(ibin);
	      htmp->SetBinContent(ibin, rndm->Gaus(mean, error));
	      htmp->SetBinError(ibin, error);
	    }
	  all = htmp->IntegralAndError(all_low_bin, all_high_bin, all_err);
	  acc = htmp->IntegralAndError(acc_low_bin, acc_high_bin, acc_err);
	  hRndmDis->Fill(acc/all);
	}
      hRndmDis->SetMarkerStyle(21);
      c = draw1D(hRndmDis,Form("Randomized efficiency for %1.1f < p_{T}^{#mu} < %1.1f GeV/c;Efficiency",xbins[bin-1],xbins[bin]));
      TF1 *funcRndm = new TF1(Form("Fit_%s",hRndmDis->GetName()),"gaus",0.8,1.2);
      hRndmDis->Fit(funcRndm,"IR0Q");
      funcRndm->SetLineColor(2);
      funcRndm->Draw("sames");
      TPaveText *t1 = GetPaveText(0.2,0.35,0.7,0.85,0.04);
      t1->SetTextAlign(11);
      t1->AddText(Form("RMS=%4.3f",hRndmDis->GetRMS()));
      t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcRndm->GetParameter(2))));
      t1->Draw();
      printf("[i] Check the mean for bin %d: %4.2f =? %4.2f\n",bin, hMatch->GetBinContent(bin), funcRndm->GetParameter(1));
      hMatch->SetBinError(bin, hRndmDis->GetRMS());
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_RndmBin%d.pdf",run_type,dtof_cut_max,bin));
    }
  
  TGraphAsymmErrors *gCountDataEff = new TGraphAsymmErrors(hMatch2, hBase,"cl=0.683 b(1,1) mode");
  gCountDataEff->SetName(Form("JpsiMuon_DtofEff%2.2f_BinCount",dtof_cut_max));
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint,x,y);
      gCountDataEff->SetPoint(ipoint,mean_pt[ipoint],y);
      gCountDataEff->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
      gCountDataEff->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
      if(hMatch->GetBinContent(ipoint+1)>0.98)
	{
	  gCountDataEff->SetPoint(ipoint,mean_pt[ipoint],hMatch->GetBinContent(ipoint+1));
	  gCountDataEff->SetPointEYhigh(ipoint, hMatch->GetBinError(ipoint+1));
	  gCountDataEff->SetPointEYlow(ipoint, hMatch->GetBinError(ipoint+1));
	}
    }

  // fitting method
  hBase->Reset();
  hMatch->Reset();
  for(int bin=1; bin<=nbins; bin++)
    {
      double all = funcFitData[bin-1]->Integral(minimum,maximum);
      double all_err = funcFitData[bin-1]->IntegralError(minimum,maximum,funcFitData[bin-1]->GetParameters(), ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
      double acc = funcFitData[bin-1]->Integral(dtof_cut_min, dtof_cut_max);
      double acc_err = funcFitData[bin-1]->IntegralError(dtof_cut_min, dtof_cut_max,funcFitData[bin-1]->GetParameters(), ptr[bin-1]->GetCovarianceMatrix().GetMatrixArray());
      hBase->SetBinContent(bin,all);
      hBase->SetBinError(bin,all_err);
      hMatch->SetBinContent(bin,acc);
      hMatch->SetBinError(bin,acc_err);
    }
  TGraphAsymmErrors *gFitDataEff = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
  gFitDataEff->SetName(Form("JpsiMuon_DtofEff%2.2f_Fitting",dtof_cut_max));
  for(int ipoint=0; ipoint<nbins; ipoint++)
    {
      gFitDataEff->GetPoint(ipoint,x,y);
      gFitDataEff->SetPoint(ipoint,mean_pt[ipoint],y);
      gFitDataEff->SetPointEXhigh(ipoint,mean_pt_err[ipoint]);
      gFitDataEff->SetPointEXlow(ipoint,mean_pt_err[ipoint]);	 
    }

   // trigger efficiency
  TH1F *hplot = new TH1F("hplot","",100,0,10);
  hplot->GetYaxis()->SetRangeUser(0.6,1.05);
  if(year==2016)
    hplot->GetYaxis()->SetRangeUser(0.5,1.16);
  hplot->GetXaxis()->SetRangeUser(1.2,7);
  hplot->SetTitle(";p_{T,#mu} (GeV/c);Efficiency");
  TCanvas *cEff = new TCanvas("DtofEff_Comp","DtofEff_Comp",800,600);
  gPad->SetGridy();
  hplot->DrawCopy();
  gFitDataEff->SetMarkerSize(1.3);
  gFitDataEff->SetMarkerStyle(24);
  gFitDataEff->SetMarkerColor(2);
  gFitDataEff->SetLineColor(2);
  gFitDataEff->Draw("PZsames");
  gCountDataEff->SetMarkerStyle(21);
  gCountDataEff->SetMarkerSize(1.3);
  TGraphAsymmErrors *gplot = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("plot_%s",gCountDataEff->GetName()));
  gplot->Draw("PZsames");
  if(year==2016)
    {
      // for the lowest pT bin, use the efficiency from fitting method
      // gFitDataEff->GetPoint(0,x,y);
      // gCountDataEff->SetPoint(0,x,y);
      // gCountDataEff->SetPointError(0, gFitDataEff->GetErrorXlow(0), gFitDataEff->GetErrorXhigh(0), 
      // 				   gFitDataEff->GetErrorYlow(0), gFitDataEff->GetErrorYhigh(0));
      // gCountDataEff->SetMarkerStyle(25);
      // gCountDataEff->SetMarkerSize(1.3);
      // gCountDataEff->Draw("PZsames");
    }

  TGraphAsymmErrors *gfit = (TGraphAsymmErrors*)gCountDataEff->Clone(Form("Fit_%s",gCountDataEff->GetName()));
  TF1 *funcEff = new TF1(Form("%s_FitFunc",gCountDataEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
  if(year==2016) funcEff->SetRange(1.5,7);
  gfit->Fit(funcEff,"R0");
  funcEff->SetLineColor(2);
  funcEff->Draw("sames");

  TPaveText *t1 = GetTitleText(Form("Efficiency for %1.0f < #DeltaTof < %1.2f ns cut",dtof_cut_min,dtof_cut_max),0.04);
  t1->Draw();
  TLegend *leg = new TLegend(0.45,0.25,0.6,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(run_type);
  leg->AddEntry(gFitDataEff,Form("Data: fitting"),"p");
  leg->AddEntry(gplot,Form("Data: bin counting"),"p");
  leg->AddEntry(funcEff,"Fit to bin counting","L");
  leg->Draw();
  if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_FitVsCount.pdf",run_type,dtof_cut_max));

  //==============================================
  // systematic uncertainty
  //==============================================
  // three sources:
  // i)   statistical error: randomization method
  // ii)  difference between default and central value of randomization
  gStyle->SetOptFit(0);
  TRandom3 *rndm = new TRandom3();
  TDatime *clock = new TDatime();
  rndm->SetSeed(clock->GetTime());
  const int nexpr = 400;
  TGraphAsymmErrors *gEffRndm[nexpr];
  double *exh = gCountDataEff->GetEXhigh();
  double *exl = gCountDataEff->GetEXlow();
  double *eyh = gCountDataEff->GetEYhigh();
  double *eyl = gCountDataEff->GetEYlow();
  for(int j = 0; j < nexpr; j++)
    {
      int npoint = gCountDataEff->GetN();
      gEffRndm[j] = new TGraphAsymmErrors(npoint);
      gEffRndm[j]->SetName(Form("%s_Rndm%d",gCountDataEff->GetName(),j));
      for(int ipoint=0; ipoint<npoint; ipoint++)
	{
	  gCountDataEff->GetPoint(ipoint,x,y);
	  double sigma = eyl[ipoint] > eyh[ipoint] ? eyl[ipoint] : eyl[ipoint];
	  y = rndm->Gaus(y, sigma);
	  gEffRndm[j]->SetPoint(ipoint,x,y);
	  gEffRndm[j]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint],eyh[ipoint]);	 
	}
    }

  hplot->GetYaxis()->SetRangeUser(0.6,1.1);
  hplot->GetXaxis()->SetRangeUser(1.3,7);
  const int nbins_rndm = 11;
  double xbins_rndm[nbins_rndm] = {1.4,1.5,1.75,2,2.25,2.5,3,4,5,6,8};
  TH1F *hSpread[nbins_rndm];
  TGraphErrors *gLimits[3];
  TF1 *funcLimits[3];
  const char *limit_name[3] = {"center","low","up"};

  c = new TCanvas(Form("cEff"),Form("cEff"),800,600);
  hplot->DrawCopy();
  TPaveText *t1 = GetTitleText("Estimation of trigger efficiency uncertainty");
  t1->Draw();
  TLegend *leg = new TLegend(0.4,0.2,0.6,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetHeader(Form("%s",run_type));

  for(int l=0; l<nbins_rndm; l++)
    {
      hSpread[l] = new TH1F(Form("hSpread_%d",l),Form("hSpread_%d",l),600,0,1.2);
    }
  for(int j = 0; j < nexpr; j++)
    {
      TF1 *funcTmp = new TF1(Form("FitFunc_%s",gEffRndm[j]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,7);
      if(year==2016) funcTmp->SetRange(1.5,7);
      funcTmp->SetParLimits(0,0,1);
      funcTmp->SetParLimits(1,0,5);
      funcTmp->SetParLimits(2,0,1);
      gEffRndm[j]->Fit(funcTmp,"R0Q");
      c->cd();
      funcTmp->SetLineStyle(2);
      funcTmp->Draw("sames");
      if(j==0) leg->AddEntry(funcTmp,"Randomization","L");
      for(int l=0; l<nbins_rndm; l++)
	{
	  hSpread[l]->Fill(funcTmp->Eval(xbins_rndm[l]));
	}
    }

  TCanvas *cFit = new TCanvas(Form("cEffFit"),Form("cEffFit"),1200,800);
  cFit->Divide(4,3);
  for(int t=0; t<3; t++) gLimits[t] = new TGraphErrors(nbins_rndm);
  for(int l=0; l<nbins_rndm; l++)
    {
      cFit->cd(l+1);
      hSpread[l]->Rebin(2);
      if(l<2) hSpread[l]->GetXaxis()->SetRangeUser(0.3,1);
      else if(l<4) hSpread[l]->GetXaxis()->SetRangeUser(0.65,1);
      else hSpread[l]->GetXaxis()->SetRangeUser(0.75,1);
      hSpread[l]->SetTitle(";Efficiency;");
      hSpread[l]->SetMaximum(1.3*hSpread[l]->GetMaximum());
      hSpread[l]->Draw();
      TPaveText *t1 = GetTitleText(Form("p_{T} = %1.1f GeV/c\n",xbins_rndm[l]),0.06);
      t1->Draw();
      TF1 *funcTmp2 = new TF1(Form("FitFunc2_%s",hSpread[l]->GetName()),"[0]*exp(-pow((x-[1])/sqrt(2)/[2],2))",0.5,1);
      funcTmp2->SetParameter(0,hSpread[l]->GetMaximum());
      funcTmp2->SetParameter(1,hSpread[l]->GetMean());
      funcTmp2->SetParameter(2,hSpread[l]->GetRMS());
      hSpread[l]->Fit(funcTmp2,"R0Q");
      funcTmp2->SetLineColor(4);
      funcTmp2->SetLineStyle(2);
      funcTmp2->Draw("sames");

      double pt    = xbins_rndm[l];
      double mean  = hSpread[l]->GetMean();
      double width = hSpread[l]->GetRMS();
      double error = hSpread[l]->GetRMSError();
      gLimits[0]->SetPoint(l,pt,mean);
      gLimits[0]->SetPointError(l,0,error);
      gLimits[1]->SetPoint(l,pt,mean-width);
      gLimits[1]->SetPointError(l,0,error);
      gLimits[2]->SetPoint(l,pt,mean+width);
      gLimits[2]->SetPointError(l,0,error);

      TPaveText *t1 = GetPaveText(0.6,0.7,0.7,0.85,0.06);
      t1->SetTextAlign(11);
      t1->AddText(Form("RMS=%4.3f",width));
      t1->AddText(Form("#sigma=%4.3f",TMath::Abs(funcTmp2->GetParameter(2))));
      t1->Draw();
    }
  cFit->cd(12);
  TPaveText *t1 = GetPaveText(0.3,0.4,0.3,0.6,0.08);
  t1->AddText(run_type);
  t1->Draw();
  if(savePlot) cFit->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_RndmFit.pdf",run_type,dtof_cut_max));
  c->cd();
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_RndmCurve.pdf",run_type,dtof_cut_max));
  for(int t=0; t<3; t++)
    {
      gLimits[t]->SetMarkerStyle(20);
      gLimits[t]->SetMarkerColor(2+2*t);
      gLimits[t]->SetLineColor(2+2*t);
      gLimits[t]->Draw("samesP");

      funcLimits[t] = new TF1(Form("%s_Sys1_%s",gCountDataEff->GetName(),limit_name[t]),"[0]-exp(-1*[1]*(x-[2]))",1,8);
      if(year==2016) funcLimits[t]->SetRange(1.5,7);
      funcLimits[t]->SetParameters(0.95,2,0.2);
      gLimits[t]->Fit(funcLimits[t],"R0Q");
      funcLimits[t]->SetLineColor(gLimits[t]->GetMarkerColor());
      if(t>0) funcLimits[t]->SetLineStyle(2);
      funcLimits[t]->DrawCopy("sames");
    }
  leg->AddEntry(gLimits[0],"Central value","PL");
  leg->AddEntry(gLimits[1],"Lower limit","PL");
  leg->AddEntry(gLimits[2],"Upper limit","PL");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_RndmLimits.pdf",run_type,dtof_cut_max));

  int npoints = gCountDataEff->GetN();
  TGraphErrors *gSys[2];
  for(int s=0; s<2; s++)
    {
      gSys[s] = new TGraphErrors(npoints);
      gSys[s]->SetName(Form("%s_S%d",gCountDataEff->GetName(),s));
      gSys[s]->SetTitle(Form("Systematic uncertainty;p_{T} (GeV/c)"));
      gSys[s]->SetMarkerStyle(20);
      gSys[s]->SetMarkerColor(s+1);
      gSys[s]->SetMarkerSize(1.3);
      gSys[s]->SetLineColor(s+1);
    }
  // source 1
  for(int ipoint=0; ipoint<npoints; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint, x, y);
      double def   = funcLimits[0]->Eval(x);
      double lower = funcLimits[1]->Eval(x);
      double upper = funcLimits[2]->Eval(x);
      double error = (upper-def) > (def-lower) ? (upper-def) : (def-lower);
      gSys[0]->SetPoint(ipoint, x, 1);
      gSys[0]->SetPointError(ipoint, 0, error/def);
    }
  // source 2
  for(int ipoint=0; ipoint<npoints; ipoint++)
    {
      gCountDataEff->GetPoint(ipoint, x, y);
      double def  = funcEff->Eval(x);
      double diff = funcLimits[0]->Eval(x);
      double error = fabs(diff/def-1);
      gSys[1]->SetPoint(ipoint, x, 1);
      gSys[1]->SetPointError(ipoint, 0, error);
    }
  gSys[0]->GetXaxis()->SetRangeUser(1.3,8);
  gSys[0]->GetYaxis()->SetRangeUser(0.94,1.06);
  c = drawGraph(gSys[0],"Systematic uncertainty for #Deltatof cut efficiency");
  gPad->SetGridy();
  TGraphErrors *gSysShift = (TGraphErrors*)gSys[1]->Clone(Form("%s_shift",gSys[1]->GetName()));
  offset_x(gSysShift, 0.05);
  gSysShift->Draw("samesPE");
  TLegend *leg = new TLegend(0.4,0.15,0.6,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gSys[0], "Source 1", "PE");
  leg->AddEntry(gSys[1], "Source 2", "PE");  
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_Sys.pdf",run_type,dtof_cut_max));

  // combine the systematic uncertainty
  const char *sys_name[2] = {"down", "up"};
  TGraphAsymmErrors *gDataEffSys[2];
  TF1 *funcEffSys[2];
  for(int i=0; i<2; i++)
    {
      int npoint = gCountDataEff->GetN();
      gDataEffSys[i] = new TGraphAsymmErrors(npoint);
      gDataEffSys[i]->SetName(Form("%s_FinalSys%d",gCountDataEff->GetName(),i));
      for(int ipoint=0; ipoint<npoint; ipoint++)
	{
	  gCountDataEff->GetPoint(ipoint,x,y);
	  double sys_all = 0;
	  for(int s=0; s<2; s++)
	    {
	      double y_sys = gSys[s]->GetErrorYhigh(ipoint);
	      sys_all += y_sys * y_sys;
	    }
	  sys_all = sqrt(sys_all);
	  double new_y = funcEff->Eval(x)+sys_all*(i*2-1);
	  gDataEffSys[i]->SetPoint(ipoint,x,new_y);
	  gDataEffSys[i]->SetPointError(ipoint,exl[ipoint],exh[ipoint],eyl[ipoint]/y*new_y,eyh[ipoint]/y*new_y);	 
	}
      gDataEffSys[i]->GetYaxis()->SetRangeUser(0.6,1.05);
      c = drawGraph(gDataEffSys[i]);
      gPad->SetGridy();
      funcEffSys[i] = new TF1(Form("%s_FitFunc_Sys%s",gCountDataEff->GetName(),sys_name[i]),"[0]-exp(-1*[1]*(x-[2]))",1.2,7);
      if(year==2016) funcEffSys[i]->SetRange(1.5,7);
      funcEffSys[i]->SetParLimits(0, 0, 1);
      funcEffSys[i]->SetParameter(1, 3.8);
      funcEffSys[i]->SetParameter(2, 0.85);
      gDataEffSys[i]->Fit(funcEffSys[i],"R0");
      funcEffSys[i]->SetLineColor(2);
      funcEffSys[i]->Draw("sames");
    }

  cEff->cd();
  for(int i=0; i<2; i++)
    {
      funcEffSys[i]->SetLineStyle(2);
      funcEffSys[i]->SetLineColor(6);
      funcEffSys[i]->Draw("sames");
    }
  TLegend *leg = new TLegend(0.45,0.2,0.6,0.25);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(funcEffSys[0], "Systematic uncertainty", "L");
  leg->Draw();
  if(savePlot) cEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Dtof/JpsiMuon_DtofEff%2.2f_FitVsCountWithSys.pdf",run_type,dtof_cut_max));
  

  //==============================================
  // save histograms
  //==============================================
  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.DtofEff.root",run_type),"update");
      // for(int bin=1; bin<=nbins; bin++)
      // 	{
      // 	  hMuonFineBin[bin-1]->Write("",TObject::kOverwrite);
      // 	  hMuon[bin-1]->Write("",TObject::kOverwrite);
      // 	  funcFitData[bin-1]->Write("",TObject::kOverwrite);
      // 	  ptr[bin-1]->Write("",TObject::kOverwrite);
      // 	}
      // hFitDataMean->Write("",TObject::kOverwrite);
      // hFitDataSigma->Write("",TObject::kOverwrite);
      gFitDataEff->Write("",TObject::kOverwrite);
      gCountDataEff->Write("",TObject::kOverwrite);
      funcEff->Write("",TObject::kOverwrite);
      for(int i=0; i<2; i++)
	{
	  funcEffSys[i]->Write("",TObject::kOverwrite);
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
}
