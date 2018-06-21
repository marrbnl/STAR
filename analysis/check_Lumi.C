const int year = YEAR;

//================================================
void check_Lumi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);


  //equiMB();
  //makeTacDiff();

  //makeVz();
  //anaVz();
}

//================================================
void anaVz(const int savePlot = 1)
{
  const int nCentBins       = nCentBins_pt; 
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  TFile *fin = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  const int nHistos = 3;
  const char* histoName[nHistos] = {"TpcVz", "DiffVz", "TpcVr"};
  const char* histoTitle[nHistos] = {"TPC v_{z}", "TPC-VPD v_{z}", "TPC v_{r}"};
  const char* legTitle[nHistos] = {"|v_{z}| < 100 cm", "|#Deltav_{z}| < 3 cm", "|v_{r}| < 2 cm"};
  const int nRunRanges = 8;
  const char* runRangeName[nRunRanges] = {"Day 074-083", "15084035-040", "Day 084-085", "Day 086-088", "Day 089-094", "Day 095-121", "Day 122-162", "Day 163-167"};

  TH1F *hVtxDis[4][nHistos][nRunRanges];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int r=0; r<nRunRanges; r++)
	    {
	      hVtxDis[k][j][r] = (TH1F*)fin->Get(Form("h%s_cent%s_%s_run%d_MB",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1],r));
	    }
	}
    }

  TCanvas *cVtx[3][2];
  for(int j=0; j<nHistos; j++)
    {
      for(int i=0; i<2; i++)
	{
	  cVtx[j][i] = new TCanvas(Form("c%s_%d",histoName[j],i), Form("c%s_%d",histoName[j],i), 1100, 700);
	  cVtx[j][i]->Divide(2,2);
	}
      for(int r=0; r<nRunRanges; r++)
	{
	  cVtx[j][r/4]->cd(r%4+1);
	  if(j>0) gPad->SetLogy();
	  TLegend *leg =  new TLegend(0.55, 0.6, 0.75, 0.85);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  leg->SetTextSize(0.045);
	  leg->SetHeader(legTitle[j]);
	  for(int k=0; k<4; k++)
	    {
	      if(hVtxDis[k][j][r]->Integral()>0) hVtxDis[k][j][r]->Scale(1./hVtxDis[k][j][r]->Integral());
	      hVtxDis[k][j][r]->SetMarkerStyle(20+k);
	      hVtxDis[k][j][r]->SetMarkerColor(color[k]);
	      hVtxDis[k][j][r]->SetLineColor(color[k]);
	      if(j==0) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(-106, -92);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(0, 0.006);
		}
	      else if(j==1) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(-10, 10);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-4, 10);
		}
	      else if(j==2) 
		{
		  hVtxDis[k][j][r]->GetXaxis()->SetRangeUser(0, 5);
		  hVtxDis[k][j][r]->GetYaxis()->SetRangeUser(1e-6, 10);
		}
	      if(k==0) hVtxDis[k][j][r]->Draw();
	      else     hVtxDis[k][j][r]->Draw("sames");
	      if(hVtxDis[k][j][r]->GetEntries()>0)
		{
		  int bin1, bin2;
		  if(j==0)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-100+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(100-1e-4);
		    }
		  else if(j==1)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(-5+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(5-1e-4);
		    }
		  else if(j==2)
		    {
		      bin1 = hVtxDis[k][j][r]->FindFixBin(0+1e-4);
		      bin2 = hVtxDis[k][j][r]->FindFixBin(2-1e-4);
		    }
		  double eff = hVtxDis[k][j][r]->Integral(bin1, bin2)/hVtxDis[k][j][r]->Integral(0,-1);
		  leg->AddEntry(hVtxDis[k][j][r], Form("%s: f = %4.1f%%",gTrgSetupTitle[k+1], eff*100), "P");
		}
	    }
	  TPaveText *t1 = GetTitleText(Form("%s MB trigger (%s%%, %s)",run_type,cent_Name_pt[0],runRangeName[r]),0.05);
	  t1->Draw();
	  leg->Draw();
	}

     for(int i=0; i<2; i++)
	{
	  if(savePlot) cVtx[j][i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/VtxEff_%sDis_%d.pdf",run_type,histoName[j],i));
	}
    }
	      
	  
  
}

//================================================
void makeVz(const int saveHisto = 1)
{
  const char *trgSetupName[4] = {"production","production_low","production_mid","production_high"};
  const int nHistos = 3;
  const char* histoName[nHistos] = {"TpcVz", "DiffVz", "TpcVr"};

  // MB data
  TFile *fMB = TFile::Open(Form("paper_code/Rootfiles/%s.data.jpsi.root",run_type),"read");
  const char* lumiName[2] = {"prod_low", "prod_high"};
  THnSparseF *hn[2];
  TH2F *hVtxDisVsRunAll[2][nHistos];
  for(int l=0; l<2; l++)
    {
      hn[l] = (THnSparseF*)fMB->Get(Form("mhMbEvtEffWeight_%s",lumiName[l]));
      hn[l]->GetAxis(5)->SetRange(1,16); // 0-80%
      for(int j=0; j<nHistos; j++)
	{
	  hVtxDisVsRunAll[l][j] = (TH2F*)hn[l]->Projection(j+1, 0);
	  hVtxDisVsRunAll[l][j]->SetName(Form("h%sVsRunAll_cent%s_%d",histoName[j],cent_Title_pt[0],l));
	}
      hn[l]->GetAxis(5)->SetRange(0,-1);
    }

  TH2F *hVtxDisVsRun[4][nHistos];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  hVtxDisVsRun[k][j] = (TH2F*)hVtxDisVsRunAll[0][j]->Clone(Form("h%sVsRun_cent%s_%s",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1]));
	  hVtxDisVsRun[k][j]->Reset();
	}
    }

  TH2F *h2tmp;
  for(int k=0; k<4; k++)
    {
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/Run14_AuAu200/AuAu_200_%s_2014.list",trgSetupName[k]));
      int runnumber;
      while(!fruns.eof())
	{
	  fruns >> runnumber;
	  int xbin = hVtxDisVsRunAll[0][0]->GetXaxis()->FindFixBin(runnumber);

	  for(int j=0; j<nHistos; j++)
	    {
	      if(k<3) h2tmp = hVtxDisVsRunAll[0][j];
	      else    h2tmp = hVtxDisVsRunAll[1][j];
	      int nybins = h2tmp->GetNbinsY();
	      for(int ybin=1; ybin<=nybins; ybin++)
		{
		  hVtxDisVsRun[k][j]->SetBinContent(xbin, ybin, h2tmp->GetBinContent(xbin, ybin));
		  hVtxDisVsRun[k][j]->SetBinError(xbin, ybin, h2tmp->GetBinError(xbin, ybin));
		}
	    }
	}
    }

  const int nRunRanges = 8;
  const int runRangeLow[nRunRanges]  = {15074000, 15084035, 15084500, 15085500, 15088500, 15094500, 15121500, 15162500};
  const int runRangeHigh[nRunRanges] = {15083500, 15084040, 15085500, 15088500, 15094500, 15121500, 15162500, 15167500};
  TH1F *hVtxDis[4][nHistos][nRunRanges];
  for(int k=0; k<4; k++)
    {
      for(int j=0; j<nHistos; j++)
	{
	  for(int r=0; r<nRunRanges; r++)
	    {
	      int bin1 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeLow[r]);
	      int bin2 = hVtxDisVsRun[k][j]->GetXaxis()->FindBin(runRangeHigh[r]);
	      hVtxDis[k][j][r] = (TH1F*)hVtxDisVsRun[k][j]->ProjectionY(Form("h%s_cent%s_%s_run%d_MB",histoName[j],cent_Title_pt[0],gTrgSetupTitle[k+1],r), bin1, bin2);
	    }
	}
    }


  if(saveHisto)
    {
      TFile *fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","update");
      for(int k=0; k<4; k++)
	{
	  for(int j=0; j<nHistos; j++)
	    {
	      for(int r=0; r<nRunRanges; r++)
		{
		  hVtxDis[k][j][r]->SetTitle("");
		  hVtxDis[k][j][r]->Write("",TObject::kOverwrite);
		}
	    }
	}
      fout->Close();
    }
}

//================================================
void makeTacDiff(const int savePlot = 1, const int saveHisto= 0)
{
  TFile *fout = 0x0;
  if(saveHisto) fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","update");
  else          fout = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  TF1 *funcppLS = (TF1*)fout->Get("fit_Run15_pp200_hTacDiff_LS_PtBin1");

  TFile *fTrig = TFile::Open("output/Run14_AuAu200.Study.MtdTrig.root","read");
  TH2F *hMT101TacDiffVsRun = (TH2F*)fTrig->Get("mhMT101TacDiffVsRun_di-muon");
  hMT101TacDiffVsRun->Sumw2();

  TH1F *hTacDiffMean = (TH1F*)hMT101TacDiffVsRun->ProjectionX("AuAu200_hMT101TacDiffMeanVsRun");
  hTacDiffMean->Reset();
  TH1F *hTacDiffSigma = (TH1F*)hMT101TacDiffVsRun->ProjectionX("AuAu200_hMT101TacDiffSigmaVsRun");
  hTacDiffSigma->Reset();
  TH1F *hTacDiffEff = (TH1F*)hMT101TacDiffVsRun->ProjectionX("AuAu200_hMT101TacDiffEffVsRun");
  hTacDiffEff->Reset();

  TCanvas *cFit[3];
  for(int i=0; i<3; i++)
    {
      cFit[i] = new TCanvas(Form("cFit_%d",i), Form("cFit_%d",i), 1100, 700);
      cFit[i]->Divide(3,3);
    }
  int nbins = hMT101TacDiffVsRun->GetNbinsX();
  int counter = 0;
  for(int bin=1; bin<=nbins; bin++)
    {
      double run = hMT101TacDiffVsRun->GetBinCenter(bin);
      if(bin%1000==0) printf("[i] Process run %1.0f\n",run);
      TH1F *h1tmp = (TH1F*)hMT101TacDiffVsRun->ProjectionY(Form("hMT101TacDiff_Run%d",bin), bin, bin);
      if(h1tmp->GetEntries()<=0) continue;
      bool isBadRun = false;
      for(int irun=0; irun<nBadRuns; irun++)
	{
	  if((int)run == badRunIDs[irun])
	    {
	      isBadRun = true;
	      break;
	    }
	}
      if(isBadRun) 
	{
	  printf("[w] Caught a bad run %1.0f\n",run);
	  continue;
	}
      TF1 *functmp = new TF1(Form("functmp_Run%d",bin),"gaus",789,800);
      if(h1tmp->GetBinContent(h1tmp->FindBin(788.5))/h1tmp->GetBinContent(h1tmp->FindBin(789.5))>0.5)
	{
	  functmp->SetRange(786,800);
	}
      if(run<=15093041 && run>=15089022 &&
	 h1tmp->GetBinContent(h1tmp->FindBin(789.5))/h1tmp->GetBinContent(h1tmp->FindBin(790.5))<0.5)
	{
	  functmp->SetRange(790,800);
	}
      h1tmp->Fit(functmp,"IR0Q");
      hTacDiffMean->SetBinContent(bin, functmp->GetParameter(1));
      hTacDiffMean->SetBinError(bin, functmp->GetParError(1));
      hTacDiffSigma->SetBinContent(bin, functmp->GetParameter(2));
      hTacDiffSigma->SetBinError(bin, functmp->GetParError(2)); 

      double tacSumCutMin = 786;
      double tacSumCutMax = 837;
      if(run<=15083031) tacSumCutMax = 823;
      if(h1tmp->GetBinContent(h1tmp->FindBin(788.5))/h1tmp->GetBinContent(h1tmp->FindBin(789.5))<0.5) tacSumCutMin = 789;
      if(run<=15093041 && run>=15089022 && h1tmp->GetBinContent(h1tmp->FindBin(789.5))/h1tmp->GetBinContent(h1tmp->FindBin(790.5))<0.5)
	{
	  tacSumCutMin = 790;
	}
      double tacSumCutMinNew = funcppLS->GetParameter(1) - (hTacDiffMean->GetBinContent(bin)-tacSumCutMin)/hTacDiffSigma->GetBinContent(bin)*funcppLS->GetParameter(2);
      double tacSumCutMaxNew = funcppLS->GetParameter(1) + (tacSumCutMax - hTacDiffMean->GetBinContent(bin))/hTacDiffSigma->GetBinContent(bin)*funcppLS->GetParameter(2);
      double eff = funcppLS->Integral(tacSumCutMinNew,tacSumCutMaxNew)/funcppLS->Integral(900,1200);
      hTacDiffEff->SetBinContent(bin, eff);
      hTacDiffEff->SetBinError(bin, 1e-10);
      
      //if(run<=15131055 && run>=15131025)
      if(run<=15164050 && run>=15164020)
	{
	  cFit[counter/9]->cd(counter%9+1);
	  h1tmp->GetXaxis()->SetRangeUser(780, 820);
	  h1tmp->SetMarkerStyle(20);
	  h1tmp->SetTitle("");
	  h1tmp->Draw();
	  functmp->SetLineColor(2);
	  functmp->SetLineStyle(2);
	  functmp->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("Run = %1.0f",hMT101TacDiffVsRun->GetBinCenter(bin)), 0.08);
	  t1->Draw();
	  TPaveText *t1 = GetPaveText(0.6,0.8,0.6,0.8,0.08);
	  t1->AddText(Form("#mu = %4.2f",functmp->GetParameter(1)));
	  t1->AddText(Form("#sigma = %4.2f",functmp->GetParameter(2)));
	  t1->Draw();
	  counter ++;
	}
    }
  if(savePlot)
    {
      for(int i=0; i<3; i++)
	{
	  cFit[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_FitTacSumDiffPerRun_%d.pdf",run_type,i));
	}
    }
  hTacDiffMean->SetMarkerStyle(20);
  c = draw1D(hTacDiffMean);

  hTacDiffSigma->SetMarkerStyle(20);
  c = draw1D(hTacDiffSigma);

  hTacDiffEff->SetMarkerStyle(20);
  c = draw1D(hTacDiffEff);
  
  if(saveHisto)
    {
      fout->cd();
      hTacDiffMean->Write("",TObject::kOverwrite);
      hTacDiffSigma->Write("",TObject::kOverwrite);    
      hTacDiffEff->Write("",TObject::kOverwrite);  
    }
}


//================================================
void equiMB(const bool savePlot = 1)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  // =============================================
  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hPSdimuon = (TH1F*)fLumi->Get("hPreScale_dimuon");
  TH1F *hLTdimuon = (TH1F*)fLumi->Get("hLiveTime_dimuon");

  // TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_BHT2-VPDMB-30");
  // TH1F *hPSdimuon = (TH1F*)fLumi->Get("hPreScale_BHT2-VPDMB-30");
  // TH1F *hLTdimuon = (TH1F*)fLumi->Get("hLiveTime_BHT2-VPDMB-30");

  TH1F *hPSmb = (TH1F*)fLumi->Get("hPreScale_VPD-ZDC-novtx-mon");
  TH1F *hLTmb = (TH1F*)fLumi->Get("hLiveTime_VPD-ZDC-novtx-mon");
  TH1F *hNeventsTakeMb = (TH1F*)fLumi->Get("hNevents_VPD-ZDC-novtx-mon");
  TH1F *hEqMbEventsGood  = (TH1F*)fLumi->Get(Form("EqMbEvtVtxCutWeight_cent%s_dimuon",cent_Title[0]));

  const char *trgSetupName[4] = {"","_low","_mid","_high"};
  TH1F *hPSdimuonInTrig[4];
  TH1F *hLTdimuonInTrig[4];
  TH1F *hNdimuonInTrig[4];
  TH1F *hPSmbInTrig[4];
  TH1F *hLTmbInTrig[4];
  TH1F *hNmbInTrig[4];
  TH1F *hEqvRatioInTrig[4];
  TH1F *hEqvRatioInTrigGood[4];

  double dimuonevent[gNTrgSetup];
  dimuonevent[0] = 0;
  for(int j=1; j<gNTrgSetup; j++)
    {
      hPSdimuonInTrig[j-1] = (TH1F*)hPSdimuon->Clone(Form("%s_%s",hPSdimuon->GetName(),trgSetupName[j-1]));
      hPSdimuonInTrig[j-1]->Reset("C");
      hLTdimuonInTrig[j-1] = (TH1F*)hLTdimuon->Clone(Form("%s_%s",hLTdimuon->GetName(),trgSetupName[j-1]));
      hLTdimuonInTrig[j-1]->Reset("C");
      hNdimuonInTrig[j-1] = (TH1F*)hNeventsTake->Clone(Form("%s_%s",hNeventsTake->GetName(),trgSetupName[j-1]));
      hNdimuonInTrig[j-1]->Reset("C");
      hPSmbInTrig[j-1] = (TH1F*)hPSmb->Clone(Form("%s_%s",hPSmb->GetName(),trgSetupName[j-1]));
      hPSmbInTrig[j-1]->Reset("C");
      hLTmbInTrig[j-1] = (TH1F*)hLTmb->Clone(Form("%s_%s",hLTmb->GetName(),trgSetupName[j-1]));
      hLTmbInTrig[j-1]->Reset("C");
      hNmbInTrig[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("%s_%s",hNeventsTakeMb->GetName(),trgSetupName[j-1]));
      hNmbInTrig[j-1]->Reset("C");
      hEqvRatioInTrig[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEqvRatioInTrig_%s",trgSetupName[j-1]));
      hEqvRatioInTrig[j-1]->Reset("AC");
      hEqvRatioInTrigGood[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEqvRatioInTrigGood_%s",trgSetupName[j-1]));
      hEqvRatioInTrigGood[j-1]->Reset("AC");

      dimuonevent[j] = 0;
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type,trgSetupName[j-1]));
      int runnumber;
      while(fruns >> runnumber)
	{
	  int bin = hEvtRunAcc->FindBin(runnumber);
	  if(bin<1 || bin>hEvtRunAcc->GetNbinsX()) continue;
	  if(hEvtRun->GetBinContent(bin)<=0) continue;
	  int lumiBin = hNeventsTake->FindFixBin(runnumber);

	  double nEventsTaken = hNeventsTake->GetBinContent(lumiBin);
	  dimuonevent[j] += nEventsTaken;
	  dimuonevent[0] += nEventsTaken;
	  if(nEventsTaken==0) 
	    {
	      printf("[w] check run %1.0f\n",runnumber);
	      continue;
	    }
	  double nEventsRun = hEvtRun->GetBinContent(bin);
	  double rf = hRF->GetBinContent(hRF->FindFixBin(runnumber));
	  if(rf==0)
	    {
	      printf("[w] rf = 0 for run %1.0f\n",run);
	      rf = 0.49;
	    }
	  double eq_mb_good = hEqMbEventsGood->GetBinContent(hEqMbEventsGood->FindFixBin(runnumber));
	  hPSdimuonInTrig[j-1]->SetBinContent(lumiBin, hPSdimuon->GetBinContent(lumiBin));
	  hPSdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hLTdimuonInTrig[j-1]->SetBinContent(lumiBin, hLTdimuon->GetBinContent(lumiBin));
	  hLTdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hNdimuonInTrig[j-1]->SetBinContent(lumiBin, hNeventsTake->GetBinContent(lumiBin));
	  hNdimuonInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hPSmbInTrig[j-1]->SetBinContent(lumiBin, hPSmb->GetBinContent(lumiBin));
	  hPSmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hLTmbInTrig[j-1]->SetBinContent(lumiBin, hLTmb->GetBinContent(lumiBin));
	  hLTmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hNmbInTrig[j-1]->SetBinContent(lumiBin, hNeventsTakeMb->GetBinContent(lumiBin));
	  hNmbInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hEqvRatioInTrig[j-1]->SetBinContent(lumiBin, hNeventsTakeMb->GetBinContent(lumiBin)*hPSmb->GetBinContent(lumiBin)/nEventsTaken/hPSdimuon->GetBinContent(lumiBin));
	  hEqvRatioInTrig[j-1]->SetBinError(lumiBin, 1e-10);
	  hEqvRatioInTrigGood[j-1]->SetBinContent(lumiBin, nEventsRun/rf * eq_mb_good/nEventsTaken/nEventsTaken);
	  hEqvRatioInTrigGood[j-1]->SetBinError(lumiBin, 1e-10);	  
	}
    }
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================

  // check the pre-scale and live-time
  TCanvas *c = new TCanvas("check_ps_lt", "check_ps_lt", 1200, 700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 

  for(int j=1; j<gNTrgSetup; j++)
    {
      c->cd(1);
      hPSdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hPSdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hPSdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hPSdimuonInTrig[j-1]->SetTitle("Dimuon: pre-scale;runId;");
      hPSdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5);
      if(j==1) hPSdimuonInTrig[j-1]->Draw("P");
      else     hPSdimuonInTrig[j-1]->Draw("samesP");
      leg->AddEntry(hPSdimuonInTrig[j-1], Form("prod%s",trgSetupName[j-1]), "P");

      c->cd(2);
      hNdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hNdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hNdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hNdimuonInTrig[j-1]->SetTitle("# of dimuon events on tape;runId;");
      hNdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e6);
      if(j==1) hNdimuonInTrig[j-1]->Draw("P");
      else     hNdimuonInTrig[j-1]->Draw("samesP");

      c->cd(3);
      hLTdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hLTdimuonInTrig[j-1]->SetMarkerColor(color[j-1]);
      hLTdimuonInTrig[j-1]->SetLineColor(color[j-1]);
      hLTdimuonInTrig[j-1]->SetTitle("Dimuon live-time");
      hLTdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTdimuonInTrig[j-1]->Draw("P");
      else     hLTdimuonInTrig[j-1]->Draw("samesP");

      c->cd(4);
      hPSmbInTrig[j-1]->SetMarkerStyle(19+j);
      hPSmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hPSmbInTrig[j-1]->SetLineColor(color[j-1]);
      hPSmbInTrig[j-1]->SetTitle("MB: pre-scale;runId;");
      hPSmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e4);
      if(j==1) hPSmbInTrig[j-1]->Draw("P");
      else     hPSmbInTrig[j-1]->Draw("samesP");

      c->cd(5);
      hNmbInTrig[j-1]->SetMarkerStyle(19+j);
      hNmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hNmbInTrig[j-1]->SetLineColor(color[j-1]);
      hNmbInTrig[j-1]->SetTitle("# of MB events on tape;runId;");
      hNmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5e4);
      if(j==1) hNmbInTrig[j-1]->Draw("P");
      else     hNmbInTrig[j-1]->Draw("samesP");

      c->cd(6);
      hLTmbInTrig[j-1]->SetMarkerStyle(19+j);
      hLTmbInTrig[j-1]->SetMarkerColor(color[j-1]);
      hLTmbInTrig[j-1]->SetLineColor(color[j-1]);
      hLTmbInTrig[j-1]->SetTitle("MB live-time");
      hLTmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTmbInTrig[j-1]->Draw("P");
      else     hLTmbInTrig[j-1]->Draw("samesP");
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_PSandLT.pdf",run_type));

  //const int firstrun = 15103035;
  //const int lastrun  = 15103063;
  const int firstrun = 0;
  const int lastrun = 1e8;
  TCanvas *c = new TCanvas("Comp_EqMbPerTrig","Comp_EqMbPerTrig",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrig[j-1]->SetMarkerStyle(19+j);
      hEqvRatioInTrig[j-1]->SetMarkerColor(color[j-1]);
      hEqvRatioInTrig[j-1]->SetLineColor(color[j-1]);
      hEqvRatioInTrig[j-1]->SetTitle(";runId;");
      hEqvRatioInTrig[j-1]->GetYaxis()->SetRangeUser(0,70);
      hEqvRatioInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) 
	{
	  hEqvRatioInTrig[j-1]->Draw("P");
	  //TPaveText *t1 = GetTitleText("Equivalent # of MB events after vertex cuts per dimuon event (0-80%)",0.035);
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events on tape per dimuon event",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrig[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBvsRun.pdf",run_type));

  // correct for MTD trigger efficiency
  TFile *flumiDep = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  TF1 *funcppLS = (TF1*)flumiDep->Get("fit_Run15_pp200_hTacDiff_LS_PtBin1");
  TH1F *hTacDiffMean = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffMeanVsRun");
  hTacDiffMean->GetYaxis()->SetRangeUser(785, 795);
  c = draw1D(hTacDiffMean,"Run14_AuAu200: mean of #DeltaTacSum distribution from MT101");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumMeanVsRun.pdf",run_type));

  TH1F *hTacDiffSigma = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffSigmaVsRun");
  hTacDiffSigma->GetYaxis()->SetRangeUser(4, 12);
  c = draw1D(hTacDiffSigma,"Run14_AuAu200: width of #DeltaTacSum distribution from MT101");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumWidthVsRun.pdf",run_type));

  TH1F *hTacDiffEff = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffEffVsRun");
  hTacDiffEff->GetYaxis()->SetRangeUser(0,1.0);
  c = draw1D(hTacDiffEff,"Run14_AuAu200: efficiency of #DeltaTacSum cut for single tracks");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumEffVsRun.pdf",run_type));

  TFile *fTrig = TFile::Open("output/Run14_AuAu200.Study.MtdTrig.root","read");
  TH2F *hMT101TacDiffVsRun = (TH2F*)fTrig->Get("mhMT101TacDiffVsRun_di-muon");
  hMT101TacDiffVsRun->Sumw2();

  // efficiency correction
  const int nRunRange = 8;
  int runRange[nRunRange+1] = {15074104, 15077035, 15078021, 15098066, 15099002, 15106130, 15131038, 15132019, 15167014};
  double accEff[nRunRange];  
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type),"read");
  TH1F *hAccLoss[nRunRange];
  for(int i=0; i<nRunRange; i++)
    {
      hAccLoss[i] = (TH1F*)fAcc->Get(Form("hAccepLoss_RunRange%d",i));
    }
  TFile *fWeight = TFile::Open("Rootfiles/models.root","read");
  TH1F *hModel = (TH1F*)fWeight->Get("TBW_JpsiYield_AuAu200_cent0060");
  TH1F *hWeight = (TH1F*)hAccLoss[0]->Clone("hWeight");
  hWeight->Reset("AC");
  for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
    {
      double low_pt  = hWeight->GetXaxis()->GetBinLowEdge(bin);
      double high_pt = hWeight->GetXaxis()->GetBinUpEdge(bin);
      double jpsi_yield = hModel->Integral(hModel->FindFixBin(low_pt+1e-5),hModel->FindFixBin(high_pt-1e-5));
      hWeight->SetBinContent(bin, jpsi_yield);
    }
  hWeight->Scale(1./hWeight->Integral());
  for(int i=0; i<nRunRange; i++)
    {
      accEff[i] = 0;
      for(int bin=1; bin<=hWeight->GetNbinsX(); bin++)
	{
	  accEff[i] += hWeight->GetBinContent(bin) * hAccLoss[i]->GetBinContent(bin);
	}
      printf("[i] Total efficiency for %d-%d is %4.3f%%\n",runRange[i],runRange[i+1],accEff[i]*100);
    }


  double eqMBevt[gNTrgSetup];
  double eqMBevtCorr[gNTrgSetup];
  double eqMBevtGood[gNTrgSetup];
  double eqMBevtGoodCorr[gNTrgSetup];
  eqMBevt[0] = 0;
  eqMBevtCorr[0] = 0;
  eqMBevtGood[0] = 0;
  eqMBevtGoodCorr[0] = 0;
  TH1F *hEqvRatioInTrigCorr[4];
  TH1F *hEqvRatioInTrigAccCorr[4];
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorr[j-1] = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_EffCorr",hEqvRatioInTrig[j-1]->GetName()));
      hEqvRatioInTrigCorr[j-1]->Reset("");
      hEqvRatioInTrigAccCorr[j-1] = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_EffAccCorr",hEqvRatioInTrig[j-1]->GetName()));
      hEqvRatioInTrigAccCorr[j-1]->Reset("");
      int nbins = hEqvRatioInTrig[j-1]->GetNbinsX();
      eqMBevt[j] = 0;
      eqMBevtCorr[j] = 0;
      eqMBevtGood[j] = 0;
      eqMBevtGoodCorr[j] = 0;
      for(int bin=1; bin<=nbins; bin++)
	{
	  double ratio = hEqvRatioInTrig[j-1]->GetBinContent(bin);
	  if(ratio<=0) continue;
	  double run = hEqvRatioInTrig[j-1]->GetBinCenter(bin);
	  int jbin = hTacDiffEff->GetXaxis()->FindBin(run);
	  double eff = hTacDiffEff->GetBinContent(jbin);
	  if(eff<=0) continue;
	  hEqvRatioInTrigCorr[j-1]->SetBinContent(bin, ratio*eff*eff);
	  hEqvRatioInTrigCorr[j-1]->SetBinError(bin, hEqvRatioInTrig[j-1]->GetBinError(bin)*eff*eff);
	  int index = -1;
	  for(int ir=0; ir<nRunRange; ir++)
	    {
	      if(run<=runRange[ir+1] && run>=runRange[ir])
		{
		  index = ir;
		  break;
		}
	    }
	  double accCorr = 1;
	  if(index>-1) accCorr = accEff[index];
	  double ratio_good = hEqvRatioInTrigGood[j-1]->GetBinContent(bin);
	  hEqvRatioInTrigAccCorr[j-1]->SetBinContent(bin, ratio*eff*eff*accCorr);
	  hEqvRatioInTrigAccCorr[j-1]->SetBinError(bin, 1e-10);
	  hEqvRatioInTrigGood[j-1]->SetBinContent(bin,  ratio_good*eff*eff*accCorr);
	  double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindBin(run));

	  eqMBevt[j] += ratio * nEventsTaken;
	  eqMBevt[0] += ratio * nEventsTaken;

	  eqMBevtCorr[j] += ratio * nEventsTaken * eff * eff * accCorr;
	  eqMBevtCorr[0] += ratio * nEventsTaken * eff * eff * accCorr;

	  eqMBevtGood[j] += ratio_good * nEventsTaken;
	  eqMBevtGood[0] += ratio_good * nEventsTaken;

	  eqMBevtGoodCorr[j] += hEqvRatioInTrigGood[j-1]->GetBinContent(bin) * nEventsTaken;
	  eqMBevtGoodCorr[0] += hEqvRatioInTrigGood[j-1]->GetBinContent(bin) * nEventsTaken;
	}
    }

  for(int j=1; j<gNTrgSetup; j++)
    {
      printf("[i] prod%10s: dm = %4.2e, eq_mb = %4.2e, eq_mb_corr = %4.2e, eq_mb_good = %4.2e, eq_mb_good_corr = %4.2e; eq_mb/dm = %2.2f, eq_mb_corr/dm = %2.2f, eq_mb_good_corr/dm = %2.2f\n",trgSetupName[j-1],dimuonevent[j],eqMBevt[j],eqMBevtCorr[j],eqMBevtGood[j],eqMBevtGoodCorr[j],eqMBevt[j]*1./dimuonevent[j],eqMBevtCorr[j]*1./dimuonevent[j],eqMBevtGoodCorr[j]*1./dimuonevent[j]);
    }

  TCanvas *c = new TCanvas("Comp_EqMbPerTrigCorr","Comp_EqMbPerTrigCorr",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorr[j-1]->GetYaxis()->SetRangeUser(5,15);
      if(j==1) 
	{
	  hEqvRatioInTrigCorr[j-1]->Draw("P");
	  //TPaveText *t1 = GetTitleText("Equivalent # of MB events after vertex cuts per dimuon event (0-80%)",0.035);
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events on tape per dimuon event (TrigEff Corr.)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigCorr[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBTrigCorrvsRun.pdf",run_type));


  TCanvas *c = new TCanvas("Comp_EqMbPerTrigAccCorr","Comp_EqMbPerTrigAccCorr",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigAccCorr[j-1]->GetYaxis()->SetRangeUser(5,15);
      if(j==1) 
	{
	  hEqvRatioInTrigAccCorr[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events on tape per dimuon event (TrigEff+Acc Corr.)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigAccCorr[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBTrigAccCorrvsRun.pdf",run_type));

  TCanvas *c = new TCanvas("Comp_EqMbPerTrigGood","Comp_EqMbPerTrigGood",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigGood[j-1]->SetMarkerStyle(19+j);
      hEqvRatioInTrigGood[j-1]->SetMarkerColor(color[j-1]);
      hEqvRatioInTrigGood[j-1]->SetLineColor(color[j-1]);
      hEqvRatioInTrigGood[j-1]->SetTitle(";runId;");
      hEqvRatioInTrigGood[j-1]->GetYaxis()->SetRangeUser(5,15);
      if(j==1) 
	{
	  hEqvRatioInTrigGood[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events in 0-80% per dimuon event (TrigEff+Acc Corr.)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigGood[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_GoodMBTrigAccCorrvsRun.pdf",run_type));
  

  return;

  // check eq. MB vs. ZDC rate
  TFile *fzdc = TFile::Open("output/Run14_AuAu200.RunDepQA.root", "read");
  TProfile *hZdc = (TProfile*)fzdc->Get("mhZDCrateVsRun_di_mu");
  TH2F *hEqvRatioVsZdc = new TH2F(Form("hEqvRatioVsZdc"),";ZDC (kHz);Eq. MB",200,0,100,60,0,1000);
  for(int j=1; j<gNTrgSetup; j++)
    {
      for(int bin=1; bin<=hEqvRatioInTrig[j-1]->GetNbinsX(); bin++)
	{
	  double run = hEqvRatioInTrig[j-1]->GetBinCenter(bin);
	  double eq = hEqvRatioInTrig[j-1]->GetBinContent(bin);
	  if(eq<=0) continue;
	  double zdc = hZdc->GetBinContent(hZdc->FindFixBin(run));
	  hEqvRatioVsZdc->Fill(zdc, eq);
	}
    }
  draw2D(hEqvRatioVsZdc);
}
