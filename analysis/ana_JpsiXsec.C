const int year = YEAR;

//================================================
void ana_JpsiXsec()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  if(year==2013)
    {
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }


  //xsec();
  //compare();
  //sfToMb();
  xsec_Run13();
}


//================================================
void xsec_Run13(const bool savePlot = 0, const bool saveHisto = 0)
{
  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %d\n",hStat->GetBinContent(1));
  printf("all di-muon events: %d\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %d\n",hStat->GetBinContent(10));
  const double dimuon_events = hStat->GetBinContent(3);

  // Luminosity
  const double sample_Lum = 28.272; // pb-1
  const double sample_Nev = 118476807;
  const double nsd_xsec = 34 * 1e6; // nb
  TFile *fLumi = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  TH1F *hMBevents = (TH1F*)fLumi->Get("hMBevents");
  const double process_MTD = hMBevents->GetBinContent(1);
  const double process_VPD = hMBevents->GetBinContent(2);

  /*
  TFile *fVtxEff = TFile::Open("Rootfiles/Run13.pp500.VPDMB.VtxCutEff.root","read");
  TH1F *hVtxEff = (TH1F*)fVtxEff->Get("hVtxEff");
  const double mb_vtx_eff = hVtxEff->GetBinContent(1);
  */
  if(run_config=="VtxCut.") const double ana_Lum = dimuon_events * process_VPD/process_MTD * 1./nsd_xsec; // nb-1
  else                      const double ana_Lum = dimuon_events * sample_Lum / sample_Nev * 1e3;

  printf("+++++++++++++++++++++++++++++++++\n");
  printf("# of dimuon events is: %d\n",dimuon_events);
  printf("# of VPD events is: %e\n",dimuon_events * process_VPD/process_MTD);
  printf("Effective luminosity %5.3f nb-1\n",ana_Lum);
  printf("+++++++++++++++++++++++++++++++++\n");

  // trigger bias
  TFile *fTrgBias = TFile::Open(Form("Rootfiles/%s.Luminosity.root",run_type),"read");
  TH1F *hTrgBias = (TH1F*)fTrgBias->Get("VPDMB_TrigBias");

  // Jpsi efficiency
  char *embedEffName = Form("Run13.pp500.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  char *trigEffName = Form("Run13.pp500.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  printf("Embed   eff: %s\n",embedEffName);
  printf("Trigger eff: %s\n",trigEffName);
  TFile *fEmbedEff = TFile::Open(Form("Rootfiles/%s",embedEffName),"read");
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/%s",trigEffName),"read");
  TH1F *hJpsiEffEmbed[nCentBins];
  TH1F *hJpsiEffTrig[nCentBins];
  TH1F *hJpsiRespEff[nCentBins];
  TH1F *hJpsiSmearEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEffEmbed[k] = (TH1F*)fEmbedEff->Get(Form("MTDreco_Jpsi_pT_%s_WeightPt_Eff_rebin",cent_Title[k]));
      hJpsiEffTrig[k] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s_rebin",cent_Title[k]));
      hJpsiRespEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiRespEff_cent%s_rebin",cent_Title[k]));
      hJpsiSmearEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiSmearEff_cent%s_rebin",cent_Title[k]));
    }

  // Jpsi raw counts
  char * yieldName = Form("Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiCounts[k] = (TH1F*)fYield->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[k]));
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)hJpsiCounts[k]->Clone(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiInvYield[k]->Divide(hJpsiEffEmbed[k]);
      //hJpsiInvYield[k]->Divide(hJpsiEffTrig[k]); // 100% trigger efficiency
      hJpsiInvYield[k]->Divide(hJpsiRespEff[k]);
      hJpsiInvYield[k]->Divide(hJpsiSmearEff[k]);
      if(run_config=="VtxCut.")  hJpsiInvYield[k]->Multiply(hTrgBias);
      cout << hJpsiEffEmbed[k]->GetBinContent(1) << "  " << hJpsiEffTrig[k]->GetBinContent(1) << endl;
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  double bin_width = hJpsiInvYield[k]->GetBinWidth(bin); // dpT
	  double bin_center = hJpsiInvYield[k]->GetBinCenter(bin); // pT 
	  hJpsiInvYield[k]->SetBinContent(bin,hJpsiInvYield[k]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[k]->SetBinError(bin,hJpsiInvYield[k]->GetBinError(bin)/bin_width/bin_center);
	}

      hJpsiInvYield[k]->Scale(1./ana_Lum); // N_evt
      hJpsiInvYield[k]->Scale(1./(2*pi)); // 2pi
      hJpsiInvYield[k]->Scale(1./1.6); // dy
      hJpsiInvYield[k]->SetMarkerStyle(21);
      hJpsiInvYield[k]->SetMarkerColor(1);
      hJpsiInvYield[k]->SetLineColor(1);
      hJpsiInvYield[k]->SetMarkerSize(1.5);
    }

  TCanvas *c = new TCanvas("c2","c2", 700, 700);
  SetPadMargin(gPad,0.12, 0.14, 0.03,0.03);
  gPad->SetLogy();
  TH1F *h = new TH1F("h2",";p_{T} (GeV/c);B#times1/(2#pip_{T})#timesd^{2}#sigma/(dp_{T}dy)   (nb/GeV/c)^{2}",1000,0,25);
  ScaleHistoTitle(h,0.045,1,0.035,0.045,1.4,0.035,62);
  h->GetYaxis()->SetRangeUser(1e-6,100);
  h->GetXaxis()->SetRangeUser(0,22.5);
  h->GetYaxis()->CenterTitle(1);
  h->DrawCopy();

  leg = new TLegend(0.45,0.65,0.65,0.95);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader("p+p @ 500 GeV");

  // CGC
  TFile *fCGC = TFile::Open("Rootfiles/Model/pp500_Jpsi_xsec/CGC.root","read");
  TGraphAsymmErrors *gCgcLowPt = (TGraphAsymmErrors*)fCGC->Get("CGC_lowPt");
  TGraphAsymmErrors *gCgcHighPt = (TGraphAsymmErrors*)fCGC->Get("CGC_highPt");
  gCgcLowPt->Draw("sames E3");
  gCgcHighPt->Draw("sames E3");
  leg->AddEntry(gCgcLowPt,"CGC+NRQCD","F");
  leg->AddEntry(gCgcHighPt,"NLO+NRQCD","F");

  /*      
  //TFile *fFit = TFile::Open("Rootfiles/GlobalFit.Jpsi.pp500.root","read");
  TFile *fFit = TFile::Open("Rootfiles/pt500GeVfit_new.root","read");
  TF1 *funcJpsi  = (TF1*)fFit->Get("ffpt");
  TF1 *funcJpsi1 = (TF1*)fFit->Get("ffpt1");
  TF1 *funcJpsi2 = (TF1*)fFit->Get("ffpt2");
  const int nPoints = 1000;
  TGraphAsymmErrors *gGJpsiSys = new TGraphAsymmErrors(nPoints);
  for(int i=0; i<nPoints; i++)
    {
      double x = 20./nPoints * (i+1);
      double y = funcJpsi->Eval(x);
      gGJpsiSys->SetPoint(i,x,y);
      gGJpsiSys->SetPointError(i, 0, 0, funcJpsi1->Eval(x)-y, y-funcJpsi2->Eval(x));
    }
  gGJpsiSys->SetFillStyle(1001);
  gGJpsiSys->SetLineColor(kGray+1);
  gGJpsiSys->SetFillColor(kGray+1);  
  gGJpsiSys->Draw("sames E3");
  funcJpsi->Draw("sames");
  leg->AddEntry(gGJpsiSys,"Global fit","F");
  */

  // Run11 e+e-
  TFile *fee = TFile::Open(Form("Rootfiles/2015HP/sptrum.root"),"read");
  TGraphErrors *gData = (TGraphErrors*)fee->Get("Jpsi_pp500");
  gData->Draw("sames PEZ");
  for(int i=0; i<19; i++)
    {
      TBox *box = (TBox*)fee->Get(Form("sys_uncert_%d",i));
      box->Draw();
    }
  leg->AddEntry(gData,"Run11: J/#psi#rightarrowe^{+}e^{-}, |y|<1","P");

  hJpsiInvYield[0]->Draw("samesP");
  leg->AddEntry(hJpsiInvYield[0],"Run13: J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");

  // Combined fit
  TCanvas *cFit = new TCanvas("cFit","cFit", 700, 700);
  SetPadMargin(gPad,0.12, 0.14, 0.03,0.03);
  gPad->SetLogy();
  h->DrawCopy();
  TH1F *hmm = (TH1F*)hJpsiInvYield[0]->Clone("hmm");
  hmm->Draw("sames");
  TH1F *hee = gData->GetHistogram();
  
  double x,y;
  for(int ipoint=0; ipoint<=gData->GetN(); ipoint++)
    {
      gData->GetPoint(ipoint, x, y);
      int bin = hee->FindFixBin(x);
      hee->SetBinContent(bin, y);
      hee->SetBinError(bin, gData->GetErrorY(ipoint));
      //cout << bin << "  " << x << "  " << y << endl;
    }
  hee->SetMarkerStyle(20);
  hee->SetMarkerColor(2);
  hee->SetLineColor(2);
  hee->Draw("sames");

  ROOT::Fit::BinData data; 
  //ROOT::Fit::FillData(data, hmm); 
  ROOT::Fit::FillData(data, hee); 

  cout << "data size is " << data.Size() << endl;

  TF1 * funcFit = new TF1("funcFit","[0]*exp([1]+[2]*x+[3]*x*x)",3.5,14);
  ROOT::Math::WrappedTF1 wf(*funcFit);
  ROOT::Fit::Fitter fitter;
  fitter.SetFunction(wf);
  fitter.Fit(data);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  funcFit->Draw("sames");

  c->cd();
  funcFit->Draw("sames");
  leg->AddEntry(funcFit,"Exp fit to Run11","L");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_Compare.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_Compare.png",run_type,run_config));
    }
  //return;


  // ratio to Global fit
  TCanvas *c = new TCanvas("pp500_Jpsi_ratio","pp500_Jpsi_ratio",800,600);
  //gPad->SetLogy();
  h->SetTitle(";p_{T} (GeV/c);Run13/Run11Fit");
  ScaleHistoTitle(h,0.045,1,0.035,0.045,1,0.035,62);
  h->GetYaxis()->SetRangeUser(0.001,2);
  h->GetXaxis()->SetRangeUser(4,18);
  h->GetYaxis()->CenterTitle(1);
  h->Draw();

  /*
  TGraphAsymmErrors *gRJpsiSys = new TGraphAsymmErrors(nPoints);
  for(int i=0; i<nPoints; i++)
    {
      double x = 20./nPoints * (i+1);
      double y = funcJpsi->Eval(x);
      gRJpsiSys->SetPoint(i,x,1);
      gRJpsiSys->SetPointError(i, 0, 0, funcJpsi1->Eval(x)/y-1, 1-funcJpsi2->Eval(x)/y);
    }
  gRJpsiSys->SetFillStyle(1001);
  gRJpsiSys->SetLineColor(kGray);
  gRJpsiSys->SetFillColor(kGray);  
  gRJpsiSys->Draw("sames E3");
  */

  TF1 *funcJpsi = (TF1*)funcFit->Clone("funcJpsi");
  funcJpsi->SetRange(0,20);
  funcJpsi->SetNpx(1000);
  TH1F *hFitJpsiPt = (TH1F*)funcJpsi->GetHistogram();
  hFitJpsiPt->SetName(Form("GlobalFit_Jpsi_Yield"));
  for(int bin=1; bin<=hFitJpsiPt->GetNbinsX(); bin++)
    {
      hFitJpsiPt->SetBinContent(bin,hFitJpsiPt->GetBinContent(bin));
    }

  TH1F *hRatio = (TH1F*)hJpsiInvYield[0]->Clone(Form("%s_ratio",hJpsiInvYield[0]->GetName()));
  for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
    {
      int start_bin = hFitJpsiPt->FindBin(hRatio->GetXaxis()->GetBinLowEdge(bin)+1e-6);
      int end_bin   = hFitJpsiPt->FindBin(hRatio->GetXaxis()->GetBinUpEdge(bin)-1e-6);
      double scale = 0;
      for(int ibin=start_bin; ibin<=end_bin; ibin++)
	{
	  scale += hFitJpsiPt->GetBinContent(ibin) * hFitJpsiPt->GetBinCenter(ibin) * hFitJpsiPt->GetBinWidth(ibin);
	}
      scale = scale / hRatio->GetBinCenter(bin) / hRatio->GetBinWidth(bin);
      cout << scale << " -> " << hRatio->GetBinContent(bin) << endl;
      hRatio->SetBinContent(bin,hRatio->GetBinContent(bin)/scale);
      hRatio->SetBinError(bin,hRatio->GetBinError(bin)/scale);
    }
  hRatio->Draw("sames");

  TGraphErrors *gRatio = (TGraphErrors*)gData->Clone("gRatio");
  double x,y;
  double *xarr = gData->GetX();
  double *yarr = gData->GetY();
  for(int ip=0; ip<gData->GetN(); ip++)
    {
      double ex = gData->GetErrorX(ip);
      double ey = gData->GetErrorY(ip);
      double scale = funcJpsi->Eval(xarr[ip]);
      //cout << xarr[ip] << "  " << yarr[ip]/scale << "  " << scale << "  " << ey << endl;
      gRatio->SetPoint(ip,xarr[ip],yarr[ip]/scale);
      gRatio->SetPointError(ip,ex,ey/scale);
     
    }
  gRatio->Draw("sames PEZ");

  TLine *line = GetLine(4,1,16,1,1);
  line->Draw();

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToCFit.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToCFit.png",run_type,run_config));
    }

  if(saveHisto)
    {
      char *outname = "";
      if(year==2013) outname = Form("Pico.Run13.pp500.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut);
      if(year==2014) outname = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut);
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outname),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiInvYield[k]->Write();
	}
    }

}


//================================================
void sfToMb(const bool savePlot = 0)
{
  const int nFiles = 3;
  const char *name[nFiles] = {"high","mid","low"};
  
  TH1F *hScale[nFiles];
  TList *list = new TList;
  TString legName[nFiles];
  for(int i=0; i<nFiles; i++)
    {
      //hScale[i] = new TH1F(Form("ScaleFactor_Prod_%s",name[i]),Form("Scale factor for production_%s;sf",name[i]),200,0,100);
      hScale[i] = new TH1F(Form("ScaleFactor_Prod_%s",name[i]),Form("Scale factor for production_%s;sf",name[i]),2000,0,2000);
      
      ifstream inFile;
      inFile.open(Form("Rootfiles/Run14.production.%s.sf.list",name[i]));
      double sf;
      int counter = 1;
      while(!inFile.eof())
	{
	  inFile >> sf;
	  //cout << sf << endl;
	  //hScale[i]->Fill(sf);
	  hScale[i]->SetBinContent(counter,sf);
	  counter++;
	}
      inFile.close();
      list->Add(hScale[i]);
      legName[i] = Form("production_%s",name[i]);
    }
  c = drawHistos(list,"ScaleFactor","Scale factor for dimuon trigger to MB;sf",kFALSE,0,30,kFALSE,1e-6,1600,kFALSE,kTRUE,legName,kTRUE,"",0.45,0.75,0.55,0.8,kTRUE);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/ScaleFactor_dimuon_to_MB.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/ScaleFactor_dimuon_to_MB.png",run_type));
    }
}

//================================================
void xsec(const bool savePlot = 1, const bool saveHisto = 1)
{
  // Get the dimuon events number
  TFile *fdata = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hStat = (TH1F*)fdata->Get("hEventStat");
  printf("all         events: %4.4e\n",hStat->GetBinContent(1));
  printf("all di-muon events: %4.4e\n",hStat->GetBinContent(3));
  printf("acc di-muon events: %4.4e\n",hStat->GetBinContent(10));
  const double dimuon_events = hStat->GetBinContent(10);

  // vertex cut efficiency
  TH1F *hTpcVz = (TH1F*)fdata->Get("mhTpcVz_di_mu");
  double vtx_eff = 1;
  if(hTpcVz)
    {
      hTpcVz->SetMaximum(1.2*hTpcVz->GetMaximum());
      c = draw1D(hTpcVz,"",kFALSE,kFALSE);
      TLine *line = GetLine(-100,0,-100,0.5*hTpcVz->GetMaximum());
      line->Draw();
      line = GetLine(100,0,100,0.5*hTpcVz->GetMaximum());
      line->Draw();
      vtx_eff = hTpcVz->Integral(hTpcVz->GetXaxis()->FindBin(-100),hTpcVz->GetXaxis()->FindBin(100))/hTpcVz->Integral(0,-1);
      TPaveText *t1 = GetPaveText(0.2,0.65,0.8,0.85,0.04,62);
      t1->AddText(Form("Fraction of events within |vz|<100cm = %4.3f%%",vtx_eff*100));
      t1->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/VtxEff.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/VtxEff.png",run_type));
	}
    }

  // MB centrality distribution
  TFile *fMB = TFile::Open("Rootfiles/VPDNovtx_Before_Ratio.root","read");
  TH1F *hRefMult = (TH1F*)fMB->Get("mgRefMultTriggerCorr_0");
  c = draw1D(hRefMult,"",kTRUE);
  double bounds[9] = {1000,401,283,193,126,77,44,23,11};
  double faction_0060 = 0;
  for(int i=1; i<9; i++)
    {
      double value = bounds[i];
      double height = hRefMult->GetBinContent(hRefMult->FindFixBin(value));
      TLine *line = GetLine(value,0,value,height,2,1,1);
      line->Draw();

      int high_bin = hRefMult->FindBin(bounds[i-1]-0.5);
      int low_bin = hRefMult->FindBin(bounds[i]+0.5);
      double nevents = hRefMult->Integral(low_bin,high_bin);
      printf("Centrality %d-%d%%: Nevents = %1.0f, fraction = %3.2f%%\n",i*10-10,i*10,nevents,nevents/hRefMult->GetEntries()*100);
      if(i<7) faction_0060 += nevents/hRefMult->GetEntries();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MB_gRefMultCorr.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MB_gRefMultCorr.png",run_type));
    }

  // MTD acceptance loss
  TH1F *hAccCounts = (TH1F*)fdata->Get("mhAccCounts_di_mu");
  double accept_corr = 1;
  if(hAccCounts)
    {
      c = draw1D(hAccCounts,"",kFALSE,kFALSE);
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MTD_acceptance_loss.pdf",run_type));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/MTD_acceptance_loss.png",run_type));
	}

      accept_corr = (hAccCounts->GetBinContent(1)*113./122 + 
		     hAccCounts->GetBinContent(2)*116./122 +
		     hAccCounts->GetBinContent(3)*120./122) / hAccCounts->GetEntries();
      accept_corr *= accept_corr;
    }

  // Effective number of MB events
  const double mb_events = dimuon_events / filter_factor * dimuon_to_mb * faction_0060 * accept_corr;

  printf("+++++++++++++++++++++++++++++++++\n");
  printf("# of dimuon events is: %4.4e\n",dimuon_events);
  printf("Vertex cut efficiency: %4.3f%%\n",vtx_eff*100);
  printf("Filter factor: %1.1f\n",filter_factor);
  printf("Scale factor: %3.4f\n",dimuon_to_mb);
  printf("Acceptance loss: %3.4f\n",sqrt(accept_corr)*100);
  printf("Effective # of MB events: %4.4e\n",mb_events);
  printf("+++++++++++++++++++++++++++++++++\n");


  // Jpsi efficiency
  char *embedEffName = Form("Run14.AuAu200.JpsiEff.%spt%1.1f.pt%1.1f.root",run_config,pt1_cut,pt2_cut);
  char *trigEffName = Form("Run14.AuAu200.JpsiTrigEff.pt%1.1f.pt%1.1f.root",pt1_cut,pt2_cut);
  printf("Embed   eff: %s\n",embedEffName);
  printf("Trigger eff: %s\n",trigEffName);
  TFile *fEmbedEff = TFile::Open(Form("Rootfiles/%s",embedEffName),"read");
  TFile *fTrigEff = TFile::Open(Form("Rootfiles/%s",trigEffName),"read");
  TH1F *hJpsiEffEmbed[nCentBins];
  TH1F *hJpsiEffTrig[nCentBins];
  TH1F *hJpsiRespEff[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiEffEmbed[k] = (TH1F*)fEmbedEff->Get(Form("MTDreco_Jpsi_pT_%s_WeightPt_Eff_rebin",cent_Title[k]));
      hJpsiEffTrig[k] = (TH1F*)fTrigEff->Get(Form("JpsiTrigEff_cent%s_rebin",cent_Title[k]));
      hJpsiRespEff[k] = (TH1F*)fTrigEff->Get(Form("JpsiRespEff_cent%s_rebin",cent_Title[k]));
    }

  // Jpsi raw counts
  char * yieldName = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.yield.root",run_config,pt1_cut,pt2_cut);
  TFile *fYield = TFile::Open(Form("Rootfiles/%s",yieldName),"read");
  cout << yieldName << endl;
  TH1F *hJpsiCounts[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiCounts[k] = (TH1F*)fYield->Get(Form("Jpsi_BinCountYield_cent%s",cent_Title[k]));
    }

  // Jpsi invariant yield
  TH1F *hJpsiInvYield[nCentBins];
  for(int k=0; k<nCentBins; k++)
    {
      hJpsiInvYield[k] = (TH1F*)hJpsiCounts[k]->Clone(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
      hJpsiInvYield[k]->Divide(hJpsiEffEmbed[k]);
      hJpsiInvYield[k]->Divide(hJpsiEffTrig[k]);
      hJpsiInvYield[k]->Divide(hJpsiRespEff[k]);
      cout << hJpsiEffEmbed[k]->GetBinContent(1) << "  " << hJpsiEffTrig[k]->GetBinContent(1) << endl;
      for(int bin=1; bin<=hJpsiInvYield[k]->GetNbinsX(); bin++)
	{
	  double bin_width = hJpsiInvYield[k]->GetBinWidth(bin); // dpT
	  double bin_center = hJpsiInvYield[k]->GetBinCenter(bin); // pT 
	  hJpsiInvYield[k]->SetBinContent(bin,hJpsiInvYield[k]->GetBinContent(bin)/bin_width/bin_center);
	  hJpsiInvYield[k]->SetBinError(bin,hJpsiInvYield[k]->GetBinError(bin)/bin_width/bin_center);
	}

      double event = mb_events * (centBins_high[k]-centBins_low[k]+1) * 1./12;
      hJpsiInvYield[k]->Scale(1./event); // N_evt
      hJpsiInvYield[k]->Scale(1./(2*pi)); // 2pi
      hJpsiInvYield[k]->Scale(1./1.6); // dy
      hJpsiInvYield[k]->SetMarkerStyle(21);
      hJpsiInvYield[k]->SetMarkerColor(2);
      hJpsiInvYield[k]->SetLineColor(2);
      hJpsiInvYield[k]->SetMarkerSize(1.5);
    }

  TFile *fpub = TFile::Open("Rootfiles/Publication.Jpsi.200GeV.root","read");
  TGraphAsymmErrors *gAuAuLowPt[nCentBins];
  TGraphAsymmErrors *gAuAuLowPtSys[nCentBins];
  TGraphAsymmErrors *gAuAuHighPt[nCentBins];
  TGraphAsymmErrors *gAuAuHighPtSys[nCentBins];
  TH1F *hTBW[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      gAuAuLowPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_cent%s",cent_Title[i]));
      gAuAuLowPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_LowPt_systematics_cent%s",cent_Title[i]));
      gAuAuHighPt[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_cent%s",cent_Title[i]));
      gAuAuHighPtSys[i] = (TGraphAsymmErrors*)fpub->Get(Form("Jpsi_InvYield_AuAu200_HighPt_systematics_cent%s",cent_Title[i]));
      hTBW[i] = (TH1F*)fpub->Get(Form("TBW_Jpsi_InvYield_cent%s",cent_Title[i]));
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-10,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      hAuAu->Draw();

      gAuAuLowPt[k]->SetMarkerStyle(24);
      gAuAuLowPt[k]->SetMarkerColor(1);
      gAuAuLowPt[k]->SetLineColor(1);
      gAuAuLowPt[k]->Draw("sames PE");
      gAuAuLowPtSys[k]->SetMarkerColor(1);
      gAuAuLowPtSys[k]->SetLineColor(1);
      gAuAuLowPtSys[k]->Draw("sameE5");
      gAuAuHighPt[k]->SetMarkerStyle(24);
      gAuAuHighPt[k]->SetMarkerColor(1);
      gAuAuHighPt[k]->SetLineColor(1);
      gAuAuHighPt[k]->Draw("sames PE");
      gAuAuHighPtSys[k]->SetMarkerColor(1);
      gAuAuHighPtSys[k]->SetLineColor(1);
      gAuAuHighPtSys[k]->Draw("sameE5");
      hTBW[k]->Draw("sames");

      hJpsiInvYield[k]->Draw("sames");
      TPaveText *t1 = GetPaveText(0.7,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }

  c->cd(1);
  TLegend *leg = new TLegend(0.18,0.2,0.42,0.48);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(gAuAuLowPt[0],"J/#psi#rightarrowe^{+}e^{-}, |y|<1","P");
  leg->AddEntry(hTBW[0],"TBW fit (#beta=0)","L");
  leg->AddEntry(hJpsiInvYield[0],"J/#psi#rightarrow#mu^{+}#mu^{-}, |y|<0.5","P");
  leg->Draw();
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_compareToPub.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_compareToPub.png",run_type,run_config));
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi_ratio","AuAu200_Jpsi_ratio",1100,700);
  c->Divide(2,2);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  for(int k=0; k<nCentBins; k++)
    {
      TH1F *hRatio = (TH1F*)hJpsiInvYield[k]->Clone(Form("%s_ratio",hJpsiInvYield[k]->GetName()));
      for(int bin=1; bin<=hRatio->GetNbinsX(); bin++)
	{
	  int start_bin = hTBW[k]->FindBin(hJpsiInvYield[k]->GetXaxis()->GetBinLowEdge(bin)+1e-6);
	  int end_bin   = hTBW[k]->FindBin(hJpsiInvYield[k]->GetXaxis()->GetBinUpEdge(bin)-1e-6);
	  double scale = 0;
	  for(int ibin=start_bin; ibin<=end_bin; ibin++)
	    {
	      scale += hTBW[k]->GetBinContent(ibin) * hTBW[k]->GetBinCenter(ibin) * hTBW[k]->GetBinWidth(ibin);
	    }
	  scale = scale / hJpsiInvYield[k]->GetBinCenter(bin) / hJpsiInvYield[k]->GetBinWidth(bin);
	  hRatio->SetBinContent(bin,hRatio->GetBinContent(bin)/scale);
	  hRatio->SetBinError(bin,hRatio->GetBinError(bin)/scale);
	}
      c->cd(k+1);
      gPad->SetGridy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
 
      hRatio->SetTitle(";p_{T} (GeV/c);Ratio to TBW");
      hRatio->GetYaxis()->SetRangeUser(0,2);
      ScaleHistoTitle(hRatio,0.06,1,0.05,0.06,1,0.05,62);
      hRatio->Draw();

      TPaveText *t1 = GetPaveText(0.2,0.3,0.85,0.9,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToTBW.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/%sJpsiInvYield_RatioToTBW.png",run_type,run_config));
    }


  if(saveHisto)
    {
      char *outname = Form("Pico.Run14.AuAu200.jpsi.%spt%1.1f.pt%1.1f.xsec.root",run_config,pt1_cut,pt2_cut);
      TFile *fout = TFile::Open(Form("Rootfiles/%s",outname),"recreate");
      for(int k=0; k<nCentBins; k++)
	{
	  hJpsiInvYield[k]->Write();
	}
    }

}

//================================================
void compare(const bool savePlot = 1)
{
  // const int nFile = 2;
  // const char *name[nFile] = {"Pico.Run14.AuAu200.jpsi.dtof1.root","Pico.Run14.AuAu200.jpsi.LooseCut.root"};
  // const TString legName[nFile] = {"dtof < 1 ns","Loose cut"};
  // const char *save_name = "dtof1VsLooseCut";

  const int nFile = 2;
  const char *name[nFile] = {"Pico.Run14.AuAu200.jpsi.dtof1.Xsec.root","Pico.Run14.AuAu200.jpsi.dtof1.Xsec.pt1.5.pt1.5.root"};
  const TString legName[nFile] = {"p_{T,1} > 1.5, p_{T,2} > 1 GeV/c","p_{T,1}, p_{T,2} > 1.5 GeV/c"};
  const char *save_name = "1GeVvs1.5GeV";

  TFile *f[nFile];
  TH1F *hYield[nCentBins][nFile];
  for(int i=0; i<nFile; i++)
    {
      f[i] = TFile::Open(Form("Rootfiles/%s",name[i]),"read");
      for(int k=0; k<nCentBins; k++)
	{
	  hYield[k][i] = (TH1F*)f[i]->Get(Form("Jpsi_InvYield_cent%s",cent_Title[k]));
	  hYield[k][i]->SetName(Form("%s_%d",hYield[k][i]->GetName(),i));
	  hYield[k][i]->SetMarkerStyle(21+i*4);
	  hYield[k][i]->SetMarkerColor(1+i);
	  hYield[k][i]->SetLineColor(1+i);
	}
    }

  TCanvas *c = new TCanvas("AuAu200_Jpsi","AuAu200_Jpsi",1100,700);
  c->Divide(2,2);
  TH1F *hAuAu = new TH1F("AuAu200_Jpsi",";p_{T} (GeV/c);d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{2}]",10,0,10);
  hAuAu->GetYaxis()->SetRangeUser(1e-10,1e-4);
  ScaleHistoTitle(hAuAu,0.06,1,0.05,0.06,1,0.05,62);
  for(int k=0; k<nCentBins; k++)
    {
      c->cd(k+1);
      gPad->SetLogy();
      SetPadMargin(gPad,0.15,0.15,0.05,0.02);
      hAuAu->Draw();
      hYield[k][0]->Draw("sames");
      hYield[k][1]->Draw("sames");

      TPaveText *t1 = GetPaveText(0.7,0.8,0.8,0.85,0.06,62);
      t1->AddText(Form("%s%%",cent_Name[k]));
      t1->Draw();
    }

  c->cd(1);
  TLegend *leg = new TLegend(0.18,0.2,0.42,0.48);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int i=0; i<nFile; i++)
    {
      leg->AddEntry(hYield[0][i],legName[i].Data(),"P");
    }
  leg->Draw();

  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/CompareXsec_%s.pdf",run_type,save_name));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_JpsiXsec/CompareXsec_%s.png",run_type,save_name));
    }
}
