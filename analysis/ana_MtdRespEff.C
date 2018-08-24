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
  //ana_embed();
  //embedVsCosmic();
  //systematics();


  compHV();
  //cosmic_QA();
}

//================================================
void cosmic_QA(const int savePlot = 0)
{						
  if(year==2014)
    {
      const int nDays = 33;
      const int days[nDays] = {50, 51, 75, 79, 86, 88, 89, 91, 92, 94, 96, 99, 100, 102, 103, 106, 108, 109, 115, 117, 125, 128, 134, 135, 137, 138, 139, 154, 155, 156, 157, 166, 174};
    }

  TFile *fin[nDays];
  TH2F *hProjTrkPtVsBL[nDays];
  TH2F *hMthTrkPtVsBL[nDays];
  TH1F *hProjTrkBtm[30][5];
  TH1F *hMthTrkBtm[30][5];
  TH1F *hMtdRespEff[30][5];
  for(int bl=1; bl<=30; bl++)
    {
      for(int mod=1; mod<=5; mod++)
	{
	  hProjTrkBtm[bl-1][mod-1] = new TH1F(Form("hProjTrkBtm_BL%d_Mod%d",bl,mod),"",nDays,0,nDays);
	  hMthTrkBtm[bl-1][mod-1]  = new TH1F(Form("hMthTrkBtm_BL%d_Mod%d",bl,mod),"",nDays,0,nDays);
	  for(int i=0; i<nDays; i++)
	    {
	      hProjTrkBtm[bl-1][mod-1]->GetXaxis()->SetBinLabel(i+1, Form("%d",days[i]));
	      hMthTrkBtm[bl-1][mod-1]->GetXaxis()->SetBinLabel(i+1, Form("%d",days[i]));
	    }
	}
    }

  TH1F *h1tmp = 0x0;
  const double ptCut = 5;
  TH1F *hStat[nDays];
  for(int i=0; i<nDays; i++)
    {
      fin[i] = TFile::Open(Form("output/cosmic/CosmicStudy.872D44E9F945FA1B27BC0A9576A401F1_%d.root",i),"read");
      hStat[i] = (TH1F*)fin[i]->Get("mhEventStat");
      hStat[i]->SetName(Form("%s_Day%d",hStat[i]->GetName(),days[i]));

      hProjTrkPtVsBL[i] = (TH2F*)fin[i]->Get("mhProjTrkPtVsBL");
      hProjTrkPtVsBL[i]->SetName(Form("%s_Day%d",hProjTrkPtVsBL[i]->GetName(),days[i]));

      hMthTrkPtVsBL[i]  = (TH2F*)fin[i]->Get("mhMthTrkPtVsBL");
      hMthTrkPtVsBL[i]->SetName(Form("%s_%d",hMthTrkPtVsBL[i]->GetName(), i));

      int nbinsx = hProjTrkPtVsBL[i]->GetNbinsX();
      int nbinsy = hProjTrkPtVsBL[i]->GetNbinsY();
      int ybinBound = hProjTrkPtVsBL[i]->GetYaxis()->FindFixBin(ptCut);
      double val, err;
      for(int bin=1; bin<=nbinsx; bin++)
	{
	  h1tmp = (TH1F*)hProjTrkPtVsBL[i]->ProjectionY(Form("hTrkPtTempProj_bin%d_day%d_tmp",bin,days[i]),bin,bin);
	  val = h1tmp->IntegralAndError(ybinBound, nbinsy, err);
	  hProjTrkBtm[(bin-1)/5][(bin-1)%5]->SetBinContent(i+1, val);
	  hProjTrkBtm[(bin-1)/5][(bin-1)%5]->SetBinError(i+1, err);

	  h1tmp = (TH1F*)hMthTrkPtVsBL[i]->ProjectionY(Form("hTrkPtTempMth_bin%d_day%d_tmp",bin,days[i]),bin,bin);
	  val = h1tmp->IntegralAndError(ybinBound, nbinsy, err);
	  hMthTrkBtm[(bin-1)/5][(bin-1)%5]->SetBinContent(i+1, val);
	  hMthTrkBtm[(bin-1)/5][(bin-1)%5]->SetBinError(i+1, err);
	}
    }
  
  TCanvas *cEff[30];
  for(int bl=1; bl<=30; bl++)
    {
      cEff[bl-1] = new TCanvas(Form("RespEff_vs_Day_BL%d",bl), Form("RespEff_vs_Day_BL%d",bl), 1100, 700);
      cEff[bl-1]->Divide(3, 2);
      for(int mod=1; mod<=5; mod++)
	{
	  hMtdRespEff[bl-1][mod-1] = GetEfficiencyCurve(hMthTrkBtm[bl-1][mod-1], hProjTrkBtm[bl-1][mod-1]);
	  hMtdRespEff[bl-1][mod-1]->SetName(Form("Cosmic_MtdRespEff_BL%d_Mod%d",bl,mod));
	  cEff[bl-1]->cd(mod);
	  hMtdRespEff[bl-1][mod-1]->GetYaxis()->SetRangeUser(0.2,1.1);
	  hMtdRespEff[bl-1][mod-1]->SetMarkerStyle(21);
	  hMtdRespEff[bl-1][mod-1]->Draw();
	  TPaveText *t1 = GetTitleText(Form("BL = %d, Mod = %d",bl,mod),0.06);
	  t1->Draw();
	}
      if(savePlot) cEff[bl-1]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_QA_RespEffVsDay_BL%d.pdf",run_type,bl));
    }

  for(int i=0; i<nDays; i++)
    {
      TH1F *hNMtdHits = (TH1F*)fin[i]->Get("mhNMtdHits");
      hNMtdHits->SetName(Form("%s_day%d",hNMtdHits->GetName(),days[i]));
      hNMtdHits->Sumw2();
      hNMtdHits->Scale(1./hStat[i]->GetBinContent(1));
      if(days[i]==96) 
	{
	  hNMtdHits->SetLineColor(2);
	  hNMtdHits->SetLineWidth(2);
	}
      if(i==0) 
	{
	  c = draw1D(hNMtdHits);
	  gPad->SetLogy();
	}
      else hNMtdHits->Draw("samesHIST");
    }
}

//================================================
void compHV(const int savePlot = 1, const int saveHisto = 1)
{
  const int nbins = 16; // Pt bins for efficiency
  const double xbins[]={0.0,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,15.0,20.};

  TFile *fin[3];
  fin[0]= TFile::Open(Form("output/Run14.cosmic.6300.root"),"read");
  fin[1]= TFile::Open(Form("output/Run14.cosmic.6400.root"),"read");
  fin[2]= TFile::Open(Form("output/Run14.cosmic.root"),"read");
  
  TH2F *hProjTrkPtVsBL[3];
  TH2F *hMthTrkPtVsBL[3];
  TH1F *hProjTrkPtBtm[3];
  TH1F *hMthTrkPtBtm[3];
  TH1F *hMtdRespEff[4];
  TF1  *funcMtdRespEff[4];

  TH1F *h1tmp = 0x0;
  for(int i=0; i<3; i++)
    {
      hProjTrkPtVsBL[i] = (TH2F*)fin[i]->Get("mhProjTrkPtVsBL");
      hProjTrkPtVsBL[i]->SetName(Form("%s_%d",hProjTrkPtVsBL[i]->GetName(), i));
      h1tmp = (TH1F*)hProjTrkPtVsBL[i]->ProjectionY(Form("hTrkPtTempProj_%d_tmp",i),45,110);
      hProjTrkPtBtm[i] = (TH1F*)h1tmp->Rebin(nbins, Form("hTrkPtTempProj_%d",i), xbins);

      hMthTrkPtVsBL[i]  = (TH2F*)fin[i]->Get("mhMthTrkPtVsBL");
      hMthTrkPtVsBL[i]->SetName(Form("%s_%d",hMthTrkPtVsBL[i]->GetName(), i));
      h1tmp = (TH1F*)hMthTrkPtVsBL[i]->ProjectionY(Form("hTrkPtTempMth_%d_tmp",i),45,110);
      hMthTrkPtBtm[i] = (TH1F*)h1tmp->Rebin(nbins, Form("hTrkPtTempMth_%d",i), xbins);

      hMtdRespEff[i] = GetEfficiencyCurve(hMthTrkPtBtm[i], hProjTrkPtBtm[i]);
      hMtdRespEff[i]->SetName(Form("Cosmic_RespEff_%d",i));
    }

  // get relative weights from data for 6300 vs. 6400
  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.root"),"read");
  printf("[i] Process %s\n", fdata->GetName());
  double mb_events[2] = {0, 0};
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hEqMbEvents  = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent0080_dimuon"));
  for(int bin=1; bin<=hEvtRunAcc->GetNbinsX(); bin++)
    {
      if(hEvtRunAcc->GetBinContent(bin)<=0) continue;
      double run = hEvtRunAcc->GetBinCenter(bin);
      double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindFixBin(run));
      if(nEventsTaken==0) continue;
      double nEventsRun = hEvtRun->GetBinContent(bin);
      double rf = hRF->GetBinContent(hRF->FindFixBin(run));
      if(rf==0) rf = 0.49;
      double eq_mb = nEventsRun/rf/nEventsTaken * hEqMbEvents->GetBinContent(hEqMbEvents->FindFixBin(run));
      if(run<15120072) mb_events[0] += eq_mb;
      else if(run<15125034) mb_events[1] += eq_mb;
      else if(run<15134021) mb_events[0] += eq_mb;
      else mb_events[1] += eq_mb;
    }
  double weights[2] = {mb_events[0]/(mb_events[0]+mb_events[1]), mb_events[1]/(mb_events[0]+mb_events[1])};
  printf("[i] weights for 6300V is %4.2f%%, for 6400V is %4.2f%%\n",weights[0]*100, weights[1]*100);
  hMtdRespEff[3] = (TH1F*)hMtdRespEff[0]->Clone(Form("Cosmic_RespEff_%d",3));
  hMtdRespEff[3]->Scale(weights[0]);
  hMtdRespEff[3]->Add(hMtdRespEff[1], weights[1]);

  TList *list = new TList;
  for(int i=0; i<4; i++)
    {
      hMtdRespEff[i]->SetMarkerStyle(20+i);
      hMtdRespEff[i]->SetLineColor(color[i]);
      hMtdRespEff[i]->SetMarkerColor(color[i]);
      list->Add(hMtdRespEff[i]);
    }
  TString legName[4] = {"HV = 6300 V", "HV = 6400 V", "Current efficiency", "Efficiency with correct weights"};
  TCanvas *c = drawHistos(list,"MtdRespEff","Run14 cosmic ray: MTD response efficiency;p_{T} (GeV/c);Eff",true,1.3,10,true,0,1,false,true,legName,true,"",0.4,0.6,0.3,0.5,true,0.04,0.04,false,0,false,false);
  for(int i=0; i<4; i++)
    {
      funcMtdRespEff[i] = new TF1(Form("Cosmic_FitRespEff_%d",i),"[0]/(x-[1])+[2]/(x-[3])+[4]",1.2,20);
      funcMtdRespEff[i]->SetParameters(-0.09,0.97,-0.03,1.12,0.9);
      hMtdRespEff[i]->Fit(funcMtdRespEff[i],"IR0");
      funcMtdRespEff[i]->SetLineColor(hMtdRespEff[i]->GetMarkerColor());
      funcMtdRespEff[i]->SetLineStyle(2);
      funcMtdRespEff[i]->Draw("sames");
    }
  if(savePlot)
    c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_CompHV.pdf",run_type));

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"update");
      hMtdRespEff[0]->Write("Cosmic_RespEff_Btm_6300",TObject::kOverwrite);
      funcMtdRespEff[0]->Write("Cosmic_RespEff_Btm_6300_fit",TObject::kOverwrite);
      hMtdRespEff[1]->Write("Cosmic_RespEff_Btm_6400",TObject::kOverwrite);
      funcMtdRespEff[1]->Write("Cosmic_RespEff_Btm_6400_fit",TObject::kOverwrite);
    }
}

//================================================
void systematics(const int savePlot = 1, const int saveHisto = 1)
{
  // const char* part_name = "Ups1S";
  // const char* part_title = "Y(1S)";
  const char* part_name = "Jpsi";
  const char* part_title = "J/psi";

  if(year==2013)
    {
      const int nPtBins = 5;
      double xPtBins[6] = {0,1.5,3,5,7,9};
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
      if(part_name=="Jpsi")
	{
	  const int nPtBins = 6;
	  double xPtBins[7] = {0,1,2,3,4,6,10};
	}
      else
	{
	  const int nPtBins = 3;
	  double xPtBins[7] = {0,2,4,10};
	}
      run_type = "Run16_AuAu200";
    }

  // open the corresponding files
  TFile *fdata = 0x0;
  if(saveHisto) fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdRespEff.root",run_type),"update");
  else          fdata = TFile::Open(Form("Rootfiles/%s.Sys.MtdRespEff.root",run_type),"read");

  TH1F *hSysMtdRespEff[2][4];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<4; j++)
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
	  if(j==3)
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
      if(year==2013)
	{
	  if(i<2) hJpsiEffPtBin[i][0]->Rebin(8);
	  else if(i!=nPtBins-1) hJpsiEffPtBin[i][0]->Rebin(2);
	}
      if(year==2014)
	{
	  if(i<2) hJpsiEffPtBin[i][0]->Rebin(4);
	  else if(i!=nPtBins-1) hJpsiEffPtBin[i][0]->Rebin(2);
	}
      if(year==2016)
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
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",xPtBins[i],xPtBins[i+1]),0.07);
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

  // Uncertainty III: matching efficiency
  TFile *fMth = TFile::Open(Form("Rootfiles/%s.Sys.MtdMthEff.root",run_type),"read");
  hSysMtdRespEff[0][2] = (TH1F*)fMth->Get(Form("%s_JpsiEffVsPt_Sys_MtdMthEff",run_type));
  hSysMtdRespEff[1][2] = (TH1F*)fMth->Get(Form("%s_JpsiEffVsCent_Sys_MtdMthEff",run_type));

  // Combine the uncertainties
  const char *legName[4] = {"Resp. eff: fit error","Resp. eff: template","Mth. eff: cosmic vs. embed (2014)","All"};
  for(int i=0; i<2; i++)
    {
      for(int bin=1; bin<=hSysMtdRespEff[i][3]->GetNbinsX(); bin++)
	{
	  double e1 = hSysMtdRespEff[i][0]->GetBinError(bin);
	  double e2 = hSysMtdRespEff[i][1]->GetBinError(bin);
	  double e3 = hSysMtdRespEff[i][2]->GetBinError(bin);
	  double error = sqrt(e1*e1 + e2*e2 + e3*e3);
	  hSysMtdRespEff[i][3]->SetBinContent(bin, 1);
	  hSysMtdRespEff[i][3]->SetBinError(bin, error);
	}
      hSysMtdRespEff[i][3]->SetMarkerStyle(21);
      if(year==2013) hSysMtdRespEff[i][3]->GetYaxis()->SetRangeUser(0.85,1.15);
      if(year==2014) hSysMtdRespEff[i][3]->GetYaxis()->SetRangeUser(0.9,1.1);
      if(year==2015) hSysMtdRespEff[i][3]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(year==2016) hSysMtdRespEff[i][3]->GetYaxis()->SetRangeUser(0.95,1.05);
      if(i==1)
	{
	  hSysMtdRespEff[i][3]->SetXTitle("");
	  hSysMtdRespEff[i][3]->GetXaxis()->SetLabelSize(0.06);
	}
      c = draw1D(hSysMtdRespEff[i][3],Form("%s: %s uncertainty due to MTD matching efficiency",run_type,part_title));
      TLegend *leg = new TLegend(0.4,0.7,0.6,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.035);
      leg->AddEntry(hSysMtdRespEff[i][3],legName[3],"P");
      for(int j=0; j<3; j++)
	{
	  hSysMtdRespEff[i][j]->SetMarkerStyle(24);
	  hSysMtdRespEff[i][j]->SetMarkerColor(color[j+1]);
	  hSysMtdRespEff[i][j]->SetLineColor(color[j+1]);
	  TGraphErrors *gr = new TGraphErrors(hSysMtdRespEff[i][j]);
	  if(i==0) offset_x(gr,0.15+j*0.15);
	  else     offset_x(gr,0.05+j*0.05);
	  gr->Draw("samesPEZ");
	  leg->AddEntry(hSysMtdRespEff[i][j],legName[j],"P");
	}
      leg->Draw();
      if(savePlot)
	{
	  if(i==0) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/%s_MtdRespEffSysVsPt.pdf",run_type,part_name));
	  if(i==1) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/%s_MtdRespEffSysVsCent.pdf",run_type,part_name));
	}
      if(gSaveAN)
	{
	  if(i==0) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_SysRespVsPt.pdf"));
	  if(i==1) c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch5_SysRespVsCent.pdf"));
	}
    }

  if(saveHisto)
    {
      fdata->cd();
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<4; j++)
	    {
	      hSysMtdRespEff[i][j]->Write("",TObject::kOverwrite);
	      hSysMtdRespEff[i][j]->Write("",TObject::kOverwrite);
	    }
	}
    }
}

//================================================
void embedVsCosmic(const int savePlot = 1)
{
  // embed vs. cosmic ray
  TFile *fin = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"read");

  TF1  *funcEmbed[30][5];
  TF1  *funcData[30][5];
  TH1F *hRatio[30][5];
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  funcEmbed[i][j] = (TF1*)fin->Get(Form("Embed_FitRespEff_BL%d_Mod%d",i+1,j+1));
	  if(i+1>9 && i+1<23) 
	    {
	      funcData[i][j] = (TF1*)fin->Get(Form("Cosmic_FitRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  else
	    {
	      funcData[i][j] = (TF1*)fin->Get(Form("Cosmic_TempRespEff_BL%d_Mod%d",i+1,j+1));
	    }
	  hRatio[i][j] = new TH1F(Form("RespEffRatio_BL%d_Mod%d",i+1,j+1), "", 100,0,20);

	  if(funcEmbed[i][j]->GetParameter(0)==0) continue;
	  for(int bin=1; bin<=hRatio[i][j]->GetNbinsX(); bin++)
	    {
	      double x = hRatio[i][j]->GetBinCenter(bin);
	      hRatio[i][j]->SetBinContent(bin, funcData[i][j]->Eval(x)/funcEmbed[i][j]->Eval(x));
	      hRatio[i][j]->SetBinError(bin, 1e-10);
	    }
	}
    }

  TH1F *hplot = new TH1F("hplot","Cosmic/Embedding: MTD response efficiency;p_{T} (GeV/c);Cosmic/Embedding",200,0,20);
  hplot->SetTitle("");
  if(year==2013) hplot->GetYaxis()->SetRangeUser(0,1.5);
  else           hplot->GetYaxis()->SetRangeUser(0.6,1);
  hplot->GetXaxis()->SetRangeUser(1.3, 10);
  TCanvas *c = new TCanvas("MtdResp_CosmicOverEmbed","MtdResp_CosmicOverEmbed",1100,700);
  c->Divide(6,5);
  int colors[5] = {1, 2, 4, 6, 8};
  for(int i=0; i<30; i++)
    {
      c->cd(i+1);
      hplot->Draw();
      for(int j=0; j<5; j++)
	{
	  if(funcEmbed[i][j]->GetParameter(0)!=0)
	    {
	      hRatio[i][j]->SetMarkerColor(colors[j]);
	      hRatio[i][j]->SetLineColor(colors[j]);
	      hRatio[i][j]->SetLineWidth(1);
	      hRatio[i][j]->Draw("sames");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("BL = %d",i+1),0.1);
      t1->Draw();
    }
  c->cd(9);
  TLegend *leg = new TLegend(0.15,0.15,0.7,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.08);
  leg->SetHeader(Form("Run%d: cosmic/embed",year-2000));
  for(int j=0; j<5; j++)
    {
      leg->AddEntry(hRatio[0][j], Form("Module %d",j+1), "L");
    }
  leg->Draw();

  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEff_toEmbed.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffResp_EmbOverCosmic.pdf"));
    }

  // compare efficiency at high pt
  TH1F *hRespEffScaleHighPt[2];
  hRespEffScaleHighPt[0] = (TH1F*)fin->Get("Embed_RespEffVsMod");
  hRespEffScaleHighPt[1] = (TH1F*)fin->Get("Cosmic_RespEffVsMod");
  TList *list = new TList();
  for(int i=0; i<2; i++)
    {
      list->Add(hRespEffScaleHighPt[i]);
    }
  TString legName[2] = {"Embedding","Cosmic"};
  if(year==2013) c = drawHistos(list,"CompRespEff",Form("Run%d: MTD response efficiency;module;Resp. Eff",year-2000),kFALSE,-100,100,kTRUE,0,1.1,kFALSE,kTRUE,legName,kTRUE,"",0.45,0.55,0.2,0.3);
  else           c = drawHistos(list,"CompRespEff",Form("Run%d: MTD response efficiency;module;Resp. Eff",year-2000),kFALSE,-100,100,kTRUE,0.4,1.1,kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.2,0.3);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEffHighPt_VsEmbed.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffResp_EmbVsCosmic.pdf"));
    }
}

//================================================
void ana_embed(const int savePlot = 0, const int saveHisto = 0)
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
  TH1F *hRespEffScaleHighPt = new TH1F("Embed_RespEffVsMod",Form("Run%d_embed: MTD response efficiency (p_{T} > 5 GeV/c);module;Resp. Eff",year-2000), 150, 0.5, 150.5);
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
	  hEff[i][j] = (TH1F*)hAcc[i][j]->Clone(Form("Embed_RespEff_BL%d_Mod%d",i+1,j+1));
	  hEff[i][j]->Divide(hAll[i][j]);
	  cEff[i/5]->cd(i%5*5+j+1);
	  SetPadMargin(gPad,0.13,0.13,0.02,0.1);
	  ScaleHistoTitle(hEff[i][j], 0.065, 0.9, 0.05, 0.06, 1.2, 0.05);
	  hEff[i][j]->SetMarkerStyle(20);
	  funcEff[i][j] = new TF1(Form("Embed_FitRespEff_BL%d_Mod%d",i+1,j+1),"[0]/([1]*x+[2]*x*x-[3])+[4]",1.1,10);
	  funcEff[i][j]->SetParameters(0,0,0,0,0);
	  if(htmp->GetEntries()<1) continue;
	  funcEff[i][j]->SetParameters(-0.1,1,1,1.1,0.98);
	  hEff[i][j]->Fit(funcEff[i][j],"IR0Q");
	  hEff[i][j]->Draw();
	  funcEff[i][j]->SetLineColor(4);
	  funcEff[i][j]->SetLineStyle(2);
	  funcEff[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("BL = %d, Mod = %d",i+1,j+1),0.07);
	  t1->Draw();
	  hRespEffScaleHighPt->SetBinContent(bin, funcEff[i][j]->Eval(10));
	  hRespEffScaleHighPt->SetBinError(bin, funcEff[i][j]->GetParError(4)/funcEff[i][j]->GetParameter(4)*funcEff[i][j]->Eval(10));
	}
    }
  if(savePlot)
    {
      for(int i=0; i<6; i++)
	{
	  cEff[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Embed_RespEff_BL%d-%d.pdf",run_type, i*5+1, i*5+5));
	}
    }
  if(gSaveAN)
    {
      for(int i=0; i<6; i++)
	{
	  cEff[i]->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffResp_EmbEff_BL%d_%d.pdf",i*5+1, i*5+5));
	}
    }

  hRespEffScaleHighPt->SetMarkerStyle(20);
  hRespEffScaleHighPt->GetYaxis()->SetRangeUser(0,1.2);
  c = draw1D(hRespEffScaleHighPt);
  if(savePlot)
   c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Embed_RespEffVsMod_HighPt.pdf",run_type));

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type),"update");
      for(int i=0; i<30; i++)
	{
	  for(int j=0; j<5; j++)
	    {
	      hEff[i][j]->Write("",TObject::kOverwrite);
	      funcEff[i][j]->Write("",TObject::kOverwrite);
	    }
	}
      hRespEffScaleHighPt->Write("",TObject::kOverwrite);
    }
}

//================================================
void ana_cosmicRay(const int savePlot = 0, const int saveHisto = 0)
{
  TFile *fin = 0x0;
  if(year==2013) fin = TFile::Open(Form("Rootfiles/Run%dResponseEffViaPtTemplate.root",year-2000),"read");
  else           fin = TFile::Open(Form("output/Run%d.cosmic.root",year-2000),"read");
  if(year>2013)
    {
      TH2F *hProjTrkPtVsBL = (TH2F*)fin->Get("mhProjTrkPtVsBL");
      TCanvas *c = draw2D(hProjTrkPtVsBL,Form("Run%d_cosmic: p_{T} of tracks projected to center of MTD module",year-2000));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_ProjTrkPtVsgMod.pdf",run_type));
      TH2F *hMthTrkPtVsBL  = (TH2F*)fin->Get("mhMthTrkPtVsBL");
      c = draw2D(hMthTrkPtVsBL,Form("Run%d_cosmic: p_{T} of tracks matched to hits in each MTD module",year-2000));
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_MatchTrkPtVsgMod.pdf",run_type));
    }

  //==============================================
  // get the template using the bottome backlegs
  //==============================================
  //const int nbins = 16; // Pt bins for efficiency
  //const double xbins[]={0.0,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,15.0,20.};
  const int nbins = 19; // Pt bins for efficiency
  const double xbins[]={0.0,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,15.0,20.};
  TH1F *hRespEffTemp;
  TF1 *funcRespEffTemp;
  TH1F *h1tmp = 0x0;
  if(year==2013)
    {
      TFile *f2014 = TFile::Open("Rootfiles/Run14_AuAu200.MtdRespEff.root", "read");
      hRespEffTemp = (TH1F*)f2014->Get("Cosmic_RespEff_Template");
      funcRespEffTemp = (TF1*)f2014->Get("Cosmic_FitRespEff_Template");  
      c = draw1D(hRespEffTemp,Form("Run%d_cosmic: efficiency template using bottom backlegs;p_{T} (GeV/c);Resp. Eff.",14));
    }
  else
    {
      h1tmp = (TH1F*)hProjTrkPtVsBL->ProjectionY("hTrkPtTempProj_tmp",45,110);
      TH1F *hTrkPtTempProj = (TH1F*)h1tmp->Rebin(nbins, "hTrkPtTempProj", xbins);
      h1tmp = (TH1F*)hMthTrkPtVsBL->ProjectionY("hTrkPtTempMth_tmp",45,110);
      TH1F *hTrkPtTempMth = (TH1F*)h1tmp->Rebin(nbins, "hTrkPtTempMth", xbins);
      TH1F *hRespEffTemp = GetEfficiencyCurve(hTrkPtTempMth, hTrkPtTempProj);
      hRespEffTemp->SetName("Cosmic_RespEff_Template");
      hRespEffTemp->SetMarkerStyle(20);
      TF1 *funcRespEffTemp = new TF1("Cosmic_FitRespEff_Template","[0]/(x-[1])+[2]/(x-[3])+[4]",1.2,20);
      //funcRespEffTemp->SetParameters(-0.09,0.97,-0.03,1.12,0.9);
      funcRespEffTemp->SetParameters(-0.09,1.04,-0.02,1.19,0.9);
      //TF1 *funcRespEffTemp = new TF1("Cosmic_FitRespEff_Template","[0]/(x-[1])+[2]",1.3,20);
      //funcRespEffTemp->SetParameters(-0.09,0.97,0.9);
      hRespEffTemp->Fit(funcRespEffTemp,"IR0");
      c = draw1D(hRespEffTemp,Form("Run%d_cosmic: efficiency template using bottom backlegs;p_{T} (GeV/c);Resp. Eff.",year-2000));
    }
  funcRespEffTemp->SetRange(1.15,20);
  funcRespEffTemp->SetLineColor(2);
  funcRespEffTemp->SetLineStyle(2);
  funcRespEffTemp->Draw("sames");
  TLegend *leg = new TLegend(0.45,0.25,0.7,0.42);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hRespEffTemp,"Cosmic ray data","P");
  leg->AddEntry(funcRespEffTemp,"Fit to cosmic","L");
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEffTemplateFit.pdf",run_type));
  if(gSaveAN)
    {
      c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffResp_TemplateFit.pdf"));
    }
  
  //==============================================
  // Fit individual module
  //==============================================
  const int nBL  = 30;
  const int nMod = 5;
  TCanvas *cBL[nBL];
  TH1F *hTrkPtProj[nBL][nMod];
  TH1F *hTrkPtMth[nBL][nMod];
  TH1F *hRespEffMod[nBL][nMod];
  TF1 *funcRespEffMod[nBL][nMod];
  TF1 *funcRespEffModHighPt[nBL][nMod];
  TF1 *funcRespEffTempScale[nBL][nMod];
  TH1F *hFitError[nBL][nMod];
  TH1F *hFitErrorHighPt[nBL][nMod];
  TH1F *hRatio[nBL][nMod];

  TH1F *hRespEffScaleHighPt = new TH1F("Cosmic_RespEffVsMod",Form("Run%d: MTD response efficiency (p_{T} > 5 GeV/c);module;Resp. Eff",year-2000), 150, 0.5, 150.5);
  const double xmin = 1.0, xmax = 20;
  TH1F *hplot = new TH1F("hplot","",200,0,20);
  const double CL = 0.68;
  TF1 *funcTmp = 0x0;
  TH1F *hErrTmp = 0x0;
  for(int i=0; i<nBL; i++)
    {
      int bl = i+1;
      //if(!(bl<=24 && bl>=8)) continue;
      //if(bl!=10) continue;

      // initilization
      for(int j=0; j<nMod; j++)
	{
	  int mod = j+1;
	  funcRespEffModHighPt[i][j] = new TF1(Form("Cosmic_FitRespEff_BL%d_Mod%d_HighPt",bl,mod),"pol0",5,20);
	  funcRespEffMod[i][j] = (TF1*)funcRespEffTemp->Clone(Form("Cosmic_FitRespEff_BL%d_Mod%d",bl,mod));
	  funcRespEffTempScale[i][j] = (TF1*)funcRespEffTemp->Clone(Form("Cosmic_TempRespEff_BL%d_Mod%d",bl,mod));
	  hRespEffMod[i][j] = (TH1F*)hRespEffTemp->Clone(Form("Cosmic_RespEff_BL%d_Mod%d",bl,mod));
	  hRespEffMod[i][j]->Reset("AC");
	  hFitError[i][j] = (TH1F*)hRespEffMod[i][j]->Clone(Form("Cosmic_FitRespEff_BL%d_Mod%d_Err",bl,mod));
	  hFitError[i][j]->Reset("AC");
	  hFitErrorHighPt[i][j] = (TH1F*)hRespEffMod[i][j]->Clone(Form("Cosmic_FitRespEff_BL%d_Mod%d_HighPt_Err",bl,mod));
	  hFitErrorHighPt[i][j]->Reset("AC");
	}

      if(bl==9 || bl==23) continue;
      if(year==2013 && ((bl>=11 && bl<=21) || bl==8 || bl==24)) continue;
      if(bl>11 && bl<21) 
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1100,700);
	  cBL[i]->Divide(3,2);
	}
      else
	{
	  cBL[i] = new TCanvas(Form("BL%d",bl),Form("BL%d",bl),1500,600);
	  cBL[i]->Divide(5,2);
	}

      
      for(int j=0; j<nMod; j++)
	{
	  int index = j;
	  if(bl>11 && bl<21)
	    {
	      if(j==0 || j==4) continue;
	      if(year==2014 && bl==15 && j==3) continue;
	      index = j-1;
	    }

	  int mod = j+1;
	  if(year>2013)
	    {
	      h1tmp = (TH1F*)hProjTrkPtVsBL->ProjectionY(Form("htmp_proj_%d_%d",bl,mod),(bl-1)*5+mod,(bl-1)*5+mod);
	      hTrkPtProj[i][j] = (TH1F*)h1tmp->Rebin(nbins, Form("hTrkPtProj_BL%d_Mod%d",bl,mod), xbins);
	      h1tmp = (TH1F*)hMthTrkPtVsBL->ProjectionY(Form("htmp_mth_%d_%d",bl,mod),(bl-1)*5+mod,(bl-1)*5+mod);
	      hTrkPtMth[i][j] = (TH1F*)h1tmp->Rebin(nbins, Form("hTrkPtMth_BL%d_Mod%d",bl,mod), xbins);
	      hRespEffMod[i][j] = GetEfficiencyCurve(hTrkPtMth[i][j], hTrkPtProj[i][j]);
	    }
	  else
	    {
	      hRespEffMod[i][j] = (TH1F*)fin->Get(Form("hPtMtdEffBkl%d_Mod%d",bl-1,mod-1));
	      hRespEffMod[i][j]->GetXaxis()->SetRangeUser(0,20);
	      hFitError[i][j] = (TH1F*)hRespEffMod[i][j]->Clone(Form("Cosmic_FitRespEff_BL%d_Mod%d_Err",bl,mod));
	      hFitError[i][j]->Reset("AC");
	      hFitErrorHighPt[i][j] = (TH1F*)hRespEffMod[i][j]->Clone(Form("Cosmic_FitRespEff_BL%d_Mod%d_HighPt_Err",bl,mod));
	      hFitErrorHighPt[i][j]->Reset("AC");
	    }
	  hRespEffMod[i][j]->SetName(Form("Cosmic_RespEff_BL%d_Mod%d",bl,mod));
	  hRespEffMod[i][j]->SetMarkerStyle(20);
	  

	  // fit high pt
	  hRespEffMod[i][j]->Fit(funcRespEffModHighPt[i][j], "IR0Q");
	  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hFitErrorHighPt[i][j], CL);
	  double chidf = TMath::Sqrt(funcRespEffModHighPt[i][j]->GetChisquare()/funcRespEffModHighPt[i][j]->GetNDF());
	  for(int ibin=1; ibin<=nbins; ibin++)
	    {
	      hFitErrorHighPt[i][j]->SetBinError(ibin, hFitErrorHighPt[i][j]->GetBinError(ibin)/chidf);
	    }
	  hRespEffScaleHighPt->SetBinContent((bl-1)*5+mod, funcRespEffModHighPt[i][j]->GetParameter(0));
	  hRespEffScaleHighPt->SetBinError  ((bl-1)*5+mod, hFitErrorHighPt[i][j]->GetBinError(1));

	  // entire pt
	  funcRespEffMod[i][j]->SetRange(1.2,20);
	  funcRespEffMod[i][j]->SetParameters(-0.1, 1.0, -0.01, 0.9, 0.93);
	  hRespEffMod[i][j]->Fit(funcRespEffMod[i][j], "R0Q");
	  funcRespEffMod[i][j]->SetParameter(4,funcRespEffTemp->GetParameter(4));
	  if(year==2014)
	    {
	      if(bl==10 && mod==3) funcRespEffMod[i][j]->SetParameters(-0.1, 1.0, -0.01, 0.9, 0.905);
	      if(bl==13 && mod==3) funcRespEffMod[i][j]->SetParameters(-0.1, 0.1, 0, 0.1, 0.85);
	      if(bl==19 && mod==3) funcRespEffMod[i][j]->SetParameters(-0.9, 0.9, 0, 0.9, 0.85);
	      if(bl==22 && mod==5) funcRespEffMod[i][j]->SetParameters(-0.1, 1.0, -0.01, 0.9, 0.905);
	    }
	  hRespEffMod[i][j]->Fit(funcRespEffMod[i][j], "R0");
	  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hFitError[i][j], CL);
	  chidf = TMath::Sqrt(funcRespEffMod[i][j]->GetChisquare()/funcRespEffMod[i][j]->GetNDF());
	  for(int ibin=1; ibin<=nbins; ibin++)
	    {
	      hFitError[i][j]->SetBinError(ibin, hFitError[i][j]->GetBinError(ibin)/chidf);
	    }

	  // scaled template
	  double scale = funcRespEffModHighPt[i][j]->Eval(10)/funcRespEffTempScale[i][j]->Eval(10);
	  funcRespEffTempScale[i][j]->SetParameter(0, funcRespEffTempScale[i][j]->GetParameter(0)*scale);
	  funcRespEffTempScale[i][j]->SetParameter(2, funcRespEffTempScale[i][j]->GetParameter(2)*scale);
	  funcRespEffTempScale[i][j]->SetParameter(4, funcRespEffTempScale[i][j]->GetParameter(4)*scale);
	  
	  // plotting
	  cBL[i]->cd(index+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  hplot->SetTitle(";p_{T} (GeV/c);Resp. Eff.");
	  ScaleHistoTitle(hplot,0.05,1,0.04,0.05,1.2,0.04,62);
	  hplot->GetXaxis()->SetRangeUser(xmin, xmax);
	  hplot->GetYaxis()->SetRangeUser(0.2, 1.1);
	  if(year==2013)
	    {
	      hplot->GetYaxis()->SetRangeUser(0, funcRespEffModHighPt[i][j]->GetParameter(0)+0.2);
	    }
	  hplot->DrawCopy();
	  if(year==2013)
	    {
	      funcTmp = funcRespEffModHighPt[i][j];
	      hErrTmp = hFitErrorHighPt[i][j];
	    }
	  else
	    {
	      if(bl>9 && bl<23)
		{
		  funcTmp = funcRespEffMod[i][j];
		  hErrTmp = hFitError[i][j];
		}
	      else
		{
		  funcTmp = funcRespEffModHighPt[i][j];
		  hErrTmp = hFitErrorHighPt[i][j];
		}
	    }
	  hErrTmp->SetFillStyle(3001);
	  hErrTmp->SetFillColor(kBlue);
	  hErrTmp->SetLineStyle(2);
	  hErrTmp->SetMarkerSize(0);
	  hErrTmp->Draw("sames e5");
	  hRespEffMod[i][j]->Draw("sames");
	  funcTmp->SetLineColor(4);
	  funcTmp->SetLineStyle(2);
	  funcTmp->SetLineWidth(2);
	  funcTmp->Draw("sames");
	  funcRespEffTempScale[i][j]->SetRange(1.15,20);
	  funcRespEffTempScale[i][j]->SetLineColor(6);
	  funcRespEffTempScale[i][j]->SetLineStyle(2);
	  funcRespEffTempScale[i][j]->SetLineWidth(2);
	  funcRespEffTempScale[i][j]->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("BL=%d, Mod=%d",bl,j+1),0.06);
	  t1->Draw();
	  TLegend *leg = new TLegend(0.45,0.2,0.7,0.42);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.05);
	  leg->AddEntry(hRespEffMod[i][j],"Cosmic data","P");
	  leg->AddEntry(funcTmp,"Fit to cosmic","L");
	  leg->AddEntry(funcRespEffTempScale[i][j],"Template","L");
	  leg->Draw();

	  // ratio plot
	  TH1F *hRelFitErr = (TH1F*)hErrTmp->Clone(Form("%s_rel",hErrTmp->GetName()));
	  for(int bin=1; bin<=hRelFitErr->GetNbinsX(); bin++)
	    {
	      double pt = hRelFitErr->GetBinCenter(bin);
	      hRelFitErr->SetBinContent(bin, 1);
	      if(pt<1.3) {
		hRelFitErr->SetBinError(bin,0);
	      }
	      else {
		hRelFitErr->SetBinError(bin, hErrTmp->GetBinError(bin)/hErrTmp->GetBinContent(bin));
	      }
	    }
	  TH1F *hDataRatio = (TH1F*)hRespEffMod[i][j]->Clone(Form("%s_ratio",hRespEffMod[i][j]->GetName()));
	  hDataRatio->Reset();
	  hDataRatio->SetMarkerStyle(20);
	  for(int bin=1; bin<=hDataRatio->GetNbinsX(); bin++)
	    {
	      double pt = hRespEffMod[i][j]->GetBinCenter(bin);
	      hDataRatio->SetBinContent(bin, hRespEffMod[i][j]->GetBinContent(bin)/funcTmp->Eval(pt));
	      hDataRatio->SetBinError(bin, hRespEffMod[i][j]->GetBinError(bin)/funcTmp->Eval(pt));
	    }

	  if(bl>11 && bl<21) cBL[i]->cd(index+3+1);
	  else cBL[i]->cd(index+5+1);
	  gPad->SetBottomMargin(0.12);
	  gPad->SetLeftMargin(0.12);
	  hplot->SetTitle(";p_{T} (GeV/c);Ratio to cosmic fit");
	  hplot->GetYaxis()->SetRangeUser(0.5,1.2);
	  //hplot->GetXaxis()->SetRangeUser(0.8,2);
	  hplot->DrawCopy("HIST");
	  TLine *line = GetLine(xmin, 1, xmax, 1, 4);
	  line->Draw();
	  hRelFitErr->Draw("sames e5");
	  hDataRatio->Draw("sames PE");
	  hRatio[i][j] = new TH1F(Form("hRatio_%d_%d",bl,j+1),"",1000,0,25);
	  for(int bin=1; bin<=hRatio[i][j]->GetNbinsX(); bin++)
	    {
	      double pt = hRatio[i][j]->GetBinCenter(bin);
	      if(pt<1.3)  hRatio[i][j]->SetBinContent(bin, 1);
	      else hRatio[i][j]->SetBinContent(bin, funcRespEffTempScale[i][j]->Eval(pt)/funcTmp->Eval(pt));
	    }
	  hRatio[i][j]->SetLineWidth(2);
	  hRatio[i][j]->SetLineStyle(2);
	  hRatio[i][j]->SetLineColor(6);
	  hRatio[i][j]->Draw("sames HIST");
	}
      if(savePlot)
	{
	  cBL[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEff_BL%d.pdf",run_type,bl));
	  cBL[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEff_BL%d.png",run_type,bl));
	}
      if(gSaveAN && (bl==1||bl==17))
	{
	  cBL[i]->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Ch4_EffResp_FitCosmic_BL%d.pdf",bl));
	}
    }

  // efficiency at high pt
  hRespEffScaleHighPt->SetMarkerStyle(20);
  hRespEffScaleHighPt->GetYaxis()->SetRangeUser(0,1);
  c = draw1D(hRespEffScaleHighPt);
  TF1 *funcRespEffVsMod = new TF1("Cosmic_FitRespEffVsMod","[0]",1,150);
  hRespEffScaleHighPt->Fit(funcRespEffVsMod,"R0");
  if(savePlot)
   c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdRespEff/Cosmic_RespEffVsMod_HighPt.pdf",run_type,bl));

  if(year==2014)
    {
      // use the average efficiency for BL8 and BL24-1
      // use the differece to the central value as the uncertainty
      int blk = 8;
      for(int j=0; j<nMod; j++)
	{
	  double scale = funcRespEffVsMod->Eval(10)/funcRespEffTempScale[blk-1][j]->Eval(10);
	  funcRespEffTempScale[blk-1][j]->SetParameter(0, funcRespEffTempScale[blk-1][j]->GetParameter(0)*scale);
	  funcRespEffTempScale[blk-1][j]->SetParameter(2, funcRespEffTempScale[blk-1][j]->GetParameter(2)*scale);
	  funcRespEffTempScale[blk-1][j]->SetParameter(4, funcRespEffTempScale[blk-1][j]->GetParameter(4)*scale);
	  for(int bin=1; bin<=hFitErrorHighPt[blk-1][j]->GetNbinsX(); bin++)
	    {
	      hFitErrorHighPt[blk-1][j]->SetBinContent(bin, funcRespEffVsMod->GetParameter(0));
	      hFitErrorHighPt[blk-1][j]->SetBinError(bin, fabs(funcRespEffVsMod->GetParameter(0)-funcRespEffModHighPt[blk-1][j]->GetParameter(0)));
	    }
	}
    }

  if(saveHisto)
    {
      TFile *fout = TFile::Open(Form("Rootfiles/%s.MtdRespEff.root",run_type), "update");
      hRespEffTemp->Write("",TObject::kOverwrite);
      funcRespEffTemp->Write("",TObject::kOverwrite);
      for(int i=0; i<nBL; i++)
	{
	  for(int j=0; j<nMod; j++)
	    {
	      hRespEffMod[i][j]->Write("",TObject::kOverwrite);
	      funcRespEffMod[i][j]->Write("",TObject::kOverwrite);
	      funcRespEffModHighPt[i][j]->Write("",TObject::kOverwrite);
	      funcRespEffTempScale[i][j]->Write("",TObject::kOverwrite);
	      hFitError[i][j]->Write("",TObject::kOverwrite);
	      hFitErrorHighPt[i][j]->Write("",TObject::kOverwrite);
	    }
	}
      hRespEffScaleHighPt->Write("",TObject::kOverwrite);
      funcRespEffVsMod->Write("",TObject::kOverwrite);
    }
}

