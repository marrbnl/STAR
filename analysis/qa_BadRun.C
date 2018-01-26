const int year = YEAR;
TFile *f;

//================================================
void qa_BadRun()
{
  gStyle->SetOptStat(0);
  if(year==2013)
    {
      f = TFile::Open("./output/Run13.pp500.jpsi.RunQA.root","read");
    }
  else if(year==2014)
    {
      f = TFile::Open("output/Run14_AuAu200.RunDepQA.root","read");
    }

  //findBadRuns();
  //badRunList();
  Run14Qa();
}

//================================================
void Run14Qa(const int savePlot = 1)
{
  // BBC vs. ZDC
  TH2F *hBBCvsZDC = (TH2F*)f->Get("mhBBCvsZDC_di_mu");
  hBBCvsZDC->GetYaxis()->SetRangeUser(0,0.12);
  c = draw2D(hBBCvsZDC,Form("%s: %s",run_type,hBBCvsZDC->GetTitle()));
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_BBCvsZDC.pdf",run_type));

  // gRefMultCorr vs Run
  TProfile *hgRefMult = (TProfile*)f->Get("mhgRefMultVsRun_di_mu");
  hgRefMult->SetMarkerStyle(24);
  hgRefMult->GetYaxis()->SetRangeUser(250, 420);
  c = draw1D(hgRefMult);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_RefMultVsRun.pdf",run_type));

  TProfile *hgRefMultCorr = (TProfile*)f->Get("mhgRefMultCorrVsRun_di_mu");
  hgRefMultCorr->SetMarkerStyle(24);
  hgRefMultCorr->GetYaxis()->SetRangeUser(380, 420);
  c = draw1D(hgRefMultCorr);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_RefMultCorrVsRun.pdf",run_type));

  // NHitsPoss
  TFile *femb = TFile::Open("output/Run14_AuAu200.Embed.Jpsi.root","read");
  THnSparseF *hnQaNHitsPoss = (THnSparseF*)femb->Get("mhQaNHitsPoss_di_mu");
  const int vz_min = 0;
  const int vz_max = 5;
  hnQaNHitsPoss->GetAxis(2)->SetRangeUser(vz_min+0.1, vz_max-0.1);
  TH1F *hNHitsPoss[3][4];
  c = new TCanvas("cNHitsPoss","cNHitsPoss",1100,700);
  c->Divide(2,2);
  const TString legName[3] = {"prod_low", "prod_mid", "prod_high"};
  for(int j=0; j<4; j++)
    {
      TLegend *leg = new TLegend(0.15,0.6,0.25,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetHeader(Form("%s",run_type));
      hnQaNHitsPoss->GetAxis(0)->SetRangeUser(j*2+0.1, j*2+1-0.1);
      for(int i=0; i<3; i++)
	{
	  hnQaNHitsPoss->GetAxis(3)->SetRange(i+2, i+2);
	  hNHitsPoss[i][j] = (TH1F*)hnQaNHitsPoss->Projection(1);
	  hNHitsPoss[i][j]->SetName(Form("hNHitsPoss_Pt%d_Lumi%d",j,i));
	  hNHitsPoss[i][j]->Sumw2();
	  hNHitsPoss[i][j]->Scale(1./hNHitsPoss[i][j]->Integral());
	  hNHitsPoss[i][j]->SetMaximum(5* hNHitsPoss[i][j]->GetMaximum());
	  hNHitsPoss[i][j]->SetMarkerStyle(20+i);
	  hNHitsPoss[i][j]->SetMarkerColor(color[i]);
	  hNHitsPoss[i][j]->SetLineColor(color[i]);
	  hNHitsPoss[i][j]->SetTitle();
	  leg->AddEntry(hNHitsPoss[i][j],Form("%s: <NHitsPoss> = %2.1f",legName[i].Data(),hNHitsPoss[i][j]->GetMean()),"p");
	  c->cd(j+1);
	  gPad->SetLogy();
	  if(i==0) hNHitsPoss[i][j]->Draw();
	  else hNHitsPoss[i][j]->Draw("sames");
	  hnQaNHitsPoss->GetAxis(3)->SetRange(0,-1);
	}
      hnQaNHitsPoss->GetAxis(0)->SetRange(0,-1);
      TPaveText *t1 = GetTitleText(Form("%d < p_{T} < %d GeV/c (%d < v_{z} < %d cm)",j*2, j*2+1,vz_min,vz_max),0.055);
      t1->Draw();
      leg->Draw();
    }
  hnQaNHitsPoss->GetAxis(2)->SetRange(0,-1);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_NHitsPoss_vz%d-%d.pdf",run_type,vz_min,vz_max));
  return;
  
  // matched track pT
  const int nRun1 = 11;
  const int badRuns1[nRun1] = {15089052, 15102011, 15110011, 15117063, 15123048, 15123050, 15128019, 15152049, 15157045, 15160040, 15163034};
  TH2F *mhQaMthTrkPtVsRun = (TH2F*)f->Get("mhQaMthTrkPtVsRun_di_mu");
  TProfile *hpro = mhQaMthTrkPtVsRun->ProfileX();
  draw2D(mhQaMthTrkPtVsRun);
  hpro->SetMarkerStyle(20);
  hpro->Draw("sames");
  TH1F *hMthTrkPt[nRun1+1];
  for(int i=0; i<nRun1; i++)
    {
      int run = badRuns1[i];
      int bin = mhQaMthTrkPtVsRun->GetXaxis()->FindBin(run);
      hMthTrkPt[i+1] = (TH1F*)mhQaMthTrkPtVsRun->ProjectionY(Form("mhQaMthTrkPt_Run%d",run),bin,bin);
      hMthTrkPt[i+1]->Sumw2();
      hMthTrkPt[i+1]->Scale(1./hMthTrkPt[i+1]->Integral());
      hMthTrkPt[i+1]->SetMarkerStyle(25);
      hMthTrkPt[i+1]->SetMarkerColor(2);
      hMthTrkPt[i+1]->SetLineColor(2);
      //printf("[i] mean = %4.2e+/-%4.2e, pro = %4.2e +/- %4.2e\n",hMthTrkPt[i+1]->GetMean(),hMthTrkPt[i+1]->GetMeanError(),hpro->GetBinContent(bin),hpro->GetBinError(bin));
    }
  hMthTrkPt[0] = (TH1F*)mhQaMthTrkPtVsRun->ProjectionY(Form("mhQaMthTrkPt_All"),0,-1);
  hMthTrkPt[0]->Sumw2();
  hMthTrkPt[0]->Scale(1./hMthTrkPt[0]->Integral());
  hMthTrkPt[0]->SetMarkerStyle(21);
  TCanvas *c1 = new TCanvas("QaMthTrkPt","QaMthTrkPt",1100,700);
  c1->Divide(4,3);
  for(int i=0; i<nRun1; i++)
    {
      c1->cd(i+1);
      gPad->SetLogy();
      hMthTrkPt[0]->Draw();
      hMthTrkPt[i+1]->Draw("sames");
    }


  // matched track phi
  TH2F *mhQaMthTrkPhiVsRun = (TH2F*)f->Get("mhQaMthTrkPhiVsRun_di_mu");
  TProfile *hpro = mhQaMthTrkPhiVsRun->ProfileX();
  draw2D(mhQaMthTrkPhiVsRun);
  hpro->SetMarkerStyle(20);
  hpro->Draw("sames");
  for(int bin=hpro->FindFixBin(15130000); bin<=hpro->FindFixBin(15134000); bin++)
    {
      if(hpro->GetBinContent(bin)<2.65 && hpro->GetBinContent(bin)>2)
	{
	  printf("[i] run %1.0f\n",hpro->GetBinCenter(bin));
	}
    }
  TH1F *hQaMthTrkPhi[2];
  hQaMthTrkPhi[0] = (TH1F*)mhQaMthTrkPhiVsRun->ProjectionY("hQaMthTrkPhiAll",0,-1);
  hQaMthTrkPhi[1] = (TH1F*)mhQaMthTrkPhiVsRun->ProjectionY("hQaMthTrkPhi",hpro->FindFixBin(15131039),hpro->FindFixBin(15132019));
  for(int i=0; i<2; i++)
    {
      hQaMthTrkPhi[i]->Sumw2();
      hQaMthTrkPhi[i]->Scale(1./hQaMthTrkPhi[i]->Integral(0,3));
      hQaMthTrkPhi[i]->SetMarkerStyle(21+i*4);
      hQaMthTrkPhi[i]->SetMarkerColor(i+1);
      hQaMthTrkPhi[i]->SetLineColor(i+1);
      hQaMthTrkPhi[i]->GetYaxis()->SetRangeUser(0, 0.6);
    }
  c = draw1D(hQaMthTrkPhi[0],"#varphi distribution of tracks matched to MTD");
  hQaMthTrkPhi[1]->Draw("sames");
  TLegend *leg = new TLegend(0.6,0.72,0.7,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hQaMthTrkPhi[0],"All runs","p");
  leg->AddEntry(hQaMthTrkPhi[1],"15131039-15132019","p");
  leg->Draw();
  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_MthTrkPhi.pdf",run_type));
    }

  // Dy distribution
  THnSparseF *hQaMthTrkDy = (THnSparseF*)f->Get("mhQaMthTrkDy_di_mu");
  hQaMthTrkDy->Sumw2();
  TH1F *hMtdDy[30][2];
  for(int i=0; i<30; i++)
    {
      hQaMthTrkDy->GetAxis(0)->SetRange(i+1, i+1);
      for(int j=0; j<2; j++)
	{
	  hQaMthTrkDy->GetAxis(3)->SetRange(j+1, j+1);
	  hMtdDy[i][j] = (TH1F*)hQaMthTrkDy->Projection(2);
	  hMtdDy[i][j]->SetName(Form("hMtdDy_BL%d_%d",i+1,j));
	  hMtdDy[i][j]->SetTitle("");
	  if(hMtdDy[i][j]->GetEntries()<=0) continue;
	  hMtdDy[i][j]->Scale(1./hMtdDy[i][j]->Integral());
	  hMtdDy[i][j]->SetMarkerStyle(21+j*4);
	  hMtdDy[i][j]->SetMarkerColor(j+1);
	  hMtdDy[i][j]->SetLineColor(j+1);
	  if(i==7 || i==23) cout << hMtdDy[i][j]->GetMean() << endl;
	}
    }
  TCanvas *c = new TCanvas("hMtdDy","hMtdDy",1100,700);
  c->Divide(6,5);
  for(int i=0; i<30; i++)
    {
      c->cd(i+1);
      hMtdDy[i][1]->Draw();
      hMtdDy[i][0]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("BL = %d",i+1),0.1);
      t1->Draw();
    }
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.08);
  leg->AddEntry(hMtdDy[0][0],"Run<15106050","L");
  leg->AddEntry(hMtdDy[0][1],"Run>15106050","L");
  c->cd(1);
  leg->Draw();

  if(savePlot) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/QA_MthDyInBl.pdf",run_type));
    }

  // NMtdHits vs. Run
  TH1F *hStat = (TH1F*)f->Get(Form("mhNeventVsRun_%s",trigName[kTrigType]));
  THnSparseF *hQaNMtdHits = (THnSparseF*)f->Get("mhQaNMtdHits_di_mu");
  TProfile *hNMtdHits[30];
  for(int i=0; i<30; i++)
    {
      hQaNMtdHits->GetAxis(0)->SetRange(i*5+1,i*5+5);
      TH2F *hNMtdHitsVsRun = (TH2F*)hQaNMtdHits->Projection(1,2);
      hNMtdHitsVsRun->SetName(Form("hNMtdHitsVsRun_BL%d",i+1));
      hNMtdHits[i] = (TProfile*)hNMtdHitsVsRun->ProfileX(Form("hNMtdHits_BL%d_Mod%d",i+1));
      hNMtdHits[i]->SetMarkerStyle(24);
    }
  const double nSigma = 0.4;
  TCanvas *c2[6];
  for(int i=0; i<6; i++)
    {
      c2[i] = new TCanvas(Form("NMtdHits_BL%d-%d",i*5+1,i*5+5),Form("NMtdHits_BL%d-%d",i*5+1,i*5+5),1100,700);
      c2[i]->Divide(3,2);
      for(int j=0; j<5; j++)
	{
	  c2[i]->cd(j+1);
	  int bl = i*5+j;
	  TH1F *hProj = new TH1F(Form("NMtdHitsDis_BL%d",bl+1),"",100,0,1);
	  int nruns = hNMtdHits[i*5+j]->GetNbinsX();
	  for(int bin=1; bin<=nruns; bin++)
	    {
	      if(hStat->GetBinContent(bin)<=0) continue;
	      hProj->Fill(hNMtdHits[i*5+j]->GetBinContent(bin));
	    }
	  double mean = hProj->GetMean();
	  double rms = hProj->GetRMS();
	  double xmin = hNMtdHits[i*5+j]->GetXaxis()->GetBinLowEdge(1);
	  double xmax = hNMtdHits[i*5+j]->GetXaxis()->GetBinUpEdge(nruns);
	  cout << mean << "  " << rms << endl;

	  hNMtdHits[i*5+j]->GetYaxis()->SetRangeUser(0.45, 0.6);
	  hNMtdHits[i*5+j]->Draw("PE");
	  double ymax = mean+nSigma*rms;
	  double ymin = mean-nSigma*rms;
	  TLine *line1 = GetLine(xmin, ymax, xmax, ymax);
	  TLine *line2 = GetLine(xmin, ymin, xmax, ymin);
	  line1->Draw();
	  line2->Draw();
	  if(bl+1!=8  && bl+1!=19)
	    {
	      for(int bin=1; bin<=nruns; bin++)
		{
		  if(hStat->GetBinContent(bin)<=0) continue;
		  double entry = hNMtdHits[i*5+j]->GetBinContent(bin);
		  if(entry > ymax || entry < ymin)
		    {
		      printf("[i] BL %d, bad run = %1.0f, entry = %4.4f\n",bl+1,hNMtdHits[i*5+j]->GetBinCenter(bin),entry);
		    }
		}
	    }
	}
    }

}

//================================================
void badRunList()
{
  TH1F *hStat = (TH1F*)f->Get(Form("mhNeventVsRun_%s",trigName[kTrigType]));
  int nEventAll = hStat->GetEntries();
  int nRunsAll = 0;
  printf("Total # of events = %e\n",nEventAll);

  // list of bad runs
  const int nFinalBadRuns = 38;
  const int finalBadRuns[nFinalBadRuns] = {15078103, 15078104, 15078107, 15078108, 15079059, 15079061, 15084002, 15084022, 15084052, 15088003, 15090006, 15097032, 15097034, 15102021, 15104018, 15104039, 15104059, 15106008, 15106009, 15106010, 15106011, 15107077, 15110032, 15110038, 15119021, 15132026, 15142054, 15151035, 15151036, 15151037, 15151038, 15151039, 15151040, 15151041, 15151042, 15151043, 15162019, 15166023};

  double nFinalBadEvents = 0;
  for(int i=0; i<nFinalBadRuns; i++)
    {
      int run = finalBadRuns[i];
      int bin = hStat->FindFixBin(run);
      double events = hStat->GetBinContent(bin);
      nFinalBadEvents += events;
      printf("[i] Run %d has %1.0f events\n",run,events);
    }
  printf("[i] Bad run list: %d runs, %1.0f events, %4.3f%%\n",nFinalBadRuns,nFinalBadEvents,nFinalBadEvents/nEventAll*100);
}

//================================================
void findBadRuns(const Bool_t save = 0)
{
  char outPDFName[256];
  sprintf(outPDFName,"PDF/%s.BadRunList.pdf",run_type);

  //----------------------------------------------------------------------------
  // Title page
  TDatime time;
  Int_t year  = time.GetYear();
  Int_t month = time.GetMonth();
  Int_t day   = time.GetDay();

  TCanvas *c1 = new TCanvas("1pad","1pad",800,600);
  SetPadMargin(gPad);


  TPaveText *t1 = new TPaveText(0.28,0.6,0.7,0.7,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.06);
  t1->AddText("MTD run-wise QA");
  t1->Draw();

  t1 = new TPaveText(0.28,0.5,0.7,0.6,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.06);
  if(year==2013) t1->AddText("for pp 500 GeV in Run13");
  if(year==2014) t1->AddText("for AuAu 200 GeV in Run14");
  t1->Draw();

  t1  = new TPaveText(0.28,0.3,0.7,0.4,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.05);
  t1->AddText(Form("%02d/%02d/%04d",month,day,year));
  t1->Draw();
  c1->Print(Form("%s(",outPDFName));

  // ----------------------------------------------------------------------------
  // Get TPDF
  TPDF *pdf = 0;
  TSeqCollection *col = gROOT->GetListOfSpecials();
  for(Int_t i=0; i<col->GetEntries(); i++)
    {
      TObject *obj = (TObject*)col->At(i);
      if( obj->IsA()->InheritsFrom("TPDF") && strcmp(obj->GetName(),outPDFName)==0 )
        pdf = dynamic_cast<TPDF*>obj;
    }
  if(!pdf) 
    {
      printf("No pointer to PDF file is available.\n");
      return;
    }
  pdf->Off();

  //----------------------------------------------------------------------------
  // analysis cuts
  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("Analysis cuts");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);

  t1->AddText("Trigger: di-muon");
  if(year==2014)
    {
      t1->AddText(Form("Vertex cut: |vz_{TPC}| < 100 cm, |vz_{TPC}-vz_{VPD}| < 3 cm, vr_{TPC} < 2 cm"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks"); 
  t1->AddText(Form("    NHitsFit >= 15"));     
  t1->AddText(Form("    Match to MTD"));
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);


  // ----------------------------------------------------------------------------
  // identify bad runs
  const int nHistos = 50;
  const int nBadRunsMax = 200;
  int nBadRuns[nHistos];
  double nBadEvents[nHistos];
  double nBadEventsNew[nHistos];
  double badRuns[nHistos][nBadRunsMax];
  double badRunStat[nHistos][nBadRunsMax];
  int isExist[nHistos][nBadRunsMax];

  for(int i=0; i<nHistos; i++)
    {
      nBadEvents[i] = 0;
      nBadEventsNew[i] = 0;
      nBadRuns[i] = 0;
      for(int j=0; j<nBadRunsMax; j++)
	{
	  badRuns[i][j] = 0;
	  badRunStat[i][j] = 0;
	  isExist[i][j] = 0;
	}
    }

  TH1F *hStat = (TH1F*)f->Get(Form("mhNeventVsRun_%s",trigName[kTrigType]));
  int nEventAll = hStat->GetEntries();
  int nRunsAll = 0;
  printf("Total # of events = %e\n",nEventAll);


  printf("\n+++++ %s +++++\n",hStat->GetName());
  TH1F *h1Clone = (TH1F*)hStat->Clone(Form("%s_clone",hStat->GetName()));
  h1Clone->Reset();
  int counter = 0;
  for(int bin=1; bin<=hStat->GetNbinsX(); bin++)
    {
      double entry = hStat->GetBinContent(bin);
      if(entry>0) nRunsAll ++;
      if(entry>0 && entry<1000) 
	{
	  badRuns[0][counter] = hStat->GetBinCenter(bin);
	  badRunStat[0][counter] = entry;
	  nBadEvents[0] += entry;
	  nBadRuns[0]++;
	  isExist[0][counter] = 0;
	  if(!isExist[0][counter])
	    {
	      nBadEventsNew[0] += entry;
	    }
	  h1Clone->SetBinContent(bin, entry);
	  h1Clone->SetBinError(bin, hStat->GetBinError(bin));
	  counter++;
	}
    }
  printf("Total # of runs = %d\n",nRunsAll);
  hStat->SetMarkerStyle(21);
  c = draw1D(hStat,"",0.04,kTRUE);
  h1Clone->SetMarkerStyle(21);
  h1Clone->SetMarkerColor(2);
  h1Clone->SetLineColor(2);
  h1Clone->Draw("sames");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/NeventVsRun_%s.pdf",run_type,trigName[kTrigType]));
    }
  printStat(nEventAll, nBadRuns[0], nBadEvents[0], badRuns[0], badRunStat[0], isExist[0], nBadEventsNew[0]);
  PaintCanvasToPDF(c, pdf);

  TList *list = (TList*)f->GetListOfKeys();
  Int_t nEntry = list->GetEntries();
  Int_t runCandidate[200], index[200];
  Int_t nCandidate = 0, nVarialbes = 0;
  const int nSigma = 4;
  const int nIter = 5;
  TCanvas *c1;
  for(Int_t i=0; i<nEntry; i++)
    {
      TKey *key = (TKey*)list->At(i);
      TObject *obj = key->ReadObj();
      if(!obj->IsA()->InheritsFrom("TProfile")) continue;
      TString objname = obj->GetName();
      TProfile *hpro = (TProfile*)obj;
      TString outname = objname;
      outname.Remove(outname.Index("_"));
      outname.ReplaceAll("mh","");
      printf("\n+++++ %s +++++\n",outname.Data());
      if(objname.Contains("BetaVsRun") || objname.Contains("MthTrkLocalyVsRun")
	 || objname.Contains("MthTrkLocalzVsRun")) continue;

      if(objname.Contains("mhNTrkVsRun") || 
	 objname.Contains("mhTrkPtVsRun") ||
	 objname.Contains("mhTrkEtaVsRun") ||
	 objname.Contains("mhTrkPhiVsRun") ||
	 objname.Contains("mhNsigmaPiVsRun") ||
	 objname.Contains("mhDedxVsRun") ||
	 objname.Contains("mhNsigmaElVsRun") ||
	 objname.Contains("mhNsigmaKaVsRun") ||
	 objname.Contains("mhNsigmaPrVsRun") ) continue;

      if(!objname.Contains("mhBBCrateVsRun") && !objname.Contains("mhZDCrateVsRun"))
	{
	  nVarialbes ++;
	  if(objname.Contains("mhMthTrkPtVsRun"))
	    {
	      TH2F *mhQaMthTrkPtVsRun = (TH2F*)f->Get("mhQaMthTrkPtVsRun_di_mu");
	      hpro = (TProfile*)mhQaMthTrkPtVsRun->ProfileX();
	    }
	  hpro->SetMarkerStyle(21);

	  double mean, rms;	  
	  TH1F *hprojection = new TH1F(Form("%s_px",hpro->GetName()),"",1000,hpro->GetYmin(),hpro->GetYmax());
	  double badRunTmp[3000];
	  int nBadRunTmp = 0;
	  for(int iter=0; iter<nIter; iter++)
	    {
	      hprojection->Reset();
	      Int_t nRun = 0;
	      Double_t variable[3000];
	      for(Int_t ibin=1; ibin<=hpro->GetNbinsX(); ibin++)
		{
		  if(hStat->GetBinContent(ibin)<=0) continue;
		  double run = hpro->GetBinCenter(ibin);
		  bool isGood = true;
		  if(iter>0)
		    {
		      for(int irun=0; irun<nBadRunTmp; irun++)
			{
			  if(run==badRunTmp[irun])
			    {
			      isGood = false;
			      break;
			    }
			}
		    }
		  if(!isGood) continue;
		  variable[nRun] = hpro->GetBinContent(ibin);
		  nRun ++;
		  hprojection->Fill(hpro->GetBinContent(ibin));
		}
	      mean = hprojection->GetMean();
	      rms = hprojection->GetRMS();
	      printf("[i] iter=%d, runs=%d, badrun=%d, mean = %4.3f, rms = %4.3f\n",iter,hprojection->GetEntries(),nBadRunTmp,mean,rms);
	      nBadRunTmp = 0;
	      for(Int_t ibin=1; ibin<=hpro->GetNbinsX(); ibin++)
		{
		  if(hStat->GetBinContent(ibin)<=0) continue;
		  double entry = hpro->GetBinContent(ibin);
		  double error = hpro->GetBinError(ibin);
		  double run = hpro->GetBinCenter(ibin);
		  if(fabs(entry-mean)>nSigma*rms && fabs(entry-mean)>nSigma*error) 
		    {
		      badRunTmp[nBadRunTmp] = run;
		      nBadRunTmp ++;
		    }
		}
	    }
	  

	  TH1F *hOutlier = (TH1F*)hStat->Clone(Form("%s_outlier",objname.Data()));
	  hOutlier->Reset("AC");
	  int counter = 0;
	  for(Int_t ibin=1; ibin<=hpro->GetNbinsX(); ibin++)
	    {
	      if(hStat->GetBinContent(ibin)<=0) continue;
	      double entry = hpro->GetBinContent(ibin);
	      double run = hpro->GetBinCenter(ibin);
	      double error = hpro->GetBinError(ibin);
	      double nevent = hStat->GetBinContent(ibin);
	      if(fabs(entry-mean)>nSigma*rms) 
		{
		  hOutlier->SetBinContent(ibin,entry);
		  hOutlier->SetBinError(ibin,error);
		  printf("run = %1.0f, entrty = %4.2f\n",run,entry);

		  if(objname.Contains("TpcVzVsRun") || objname.Contains("VpdVzVsRun"))
		    {
		      if(run>=15120011) continue;
		    }
		  nBadRuns[i+1]++;
		  badRuns[i+1][counter] = run;
		  nBadEvents[i+1] += nevent;
		  badRunStat[i+1][counter] = nevent;
		  isExist[i+1][counter] = 0;
		  for(int j=0; j<i+1; j++)
		    {
		      if(isExist[i+1][counter]) break;
		      for(int irun=0; irun<nBadRuns[j]; irun++)
			{
			  if(run==badRuns[j][irun])
			    {
			      isExist[i+1][counter] = 1;
			      break;
			    }
			}
		    }
		  if(!isExist[i+1][counter])
		    {
		      nBadEventsNew[i+1] += nevent;
		    }
		  counter++;
		}
	      else
		{
		  hOutlier->SetBinContent(ibin,-999);
		  hOutlier->SetBinError(ibin,1);
		}
	    }
	  printStat(nEventAll, nBadRuns[i+1], nBadEvents[i+1], badRuns[i+1], badRunStat[i+1], isExist[i+1], nBadEventsNew[i+1]);
	  hOutlier->SetMarkerStyle(21);
	  hOutlier->SetMarkerColor(2);
	  hOutlier->SetLineColor(2);
	  hpro->GetYaxis()->SetRangeUser(mean-15*rms, mean + 15*rms);
	  if(objname.Contains("mhMthTrkPtVsRun")) hpro->GetYaxis()->SetRangeUser(1.45, 1.75);
	  c1 = draw1D(hpro,"",false,true);
	  c1->cd();
	  hOutlier->Draw("sameP");
	}
      else
	{
	  hpro->SetMarkerStyle(21);
	  c1 = draw1D(hpro,"",false,true);
	  if(objname.Contains("ZDCrateVsRun"))
	    {
	      for(int bin=1; bin<=10000; bin++)
		{
		  if(hpro->GetBinContent(bin)>80)
		    {
		      printf("[i] run = %1.0f, ZDC = %4.2f, n = %f\n",hpro->GetBinCenter(bin),hpro->GetBinContent(bin),hStat->GetBinContent(bin));
		    }
		}
	    }
	}
      PaintCanvasToPDF(c1, pdf);
      if(save) 
	{
	  c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_BadRunList/%s_pro_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
  pdf->On();
  pdf->Close();
}

//================================================
void printStat(const double nEventAll, const int nBadRuns, const double nBadEvents, const double *badRuns, const double *badRunStat, const int *isExist, const double nBadEventsNew)
{
  printf("[i] # of bad runs: %d, events = %4.3e (%4.3e%%)\n",nBadRuns,nBadEvents,nBadEvents/nEventAll*100);
  printf("[i] New events: %4.3e (%4.3e%%)\n",nBadEventsNew,nBadEventsNew/nEventAll*100);
  for(int i=0; i<nBadRuns; i++)
    {
      printf("Bad run: %1.0f, %4.3e, %d\n",badRuns[i],badRunStat[i],isExist[i]);
    }
  printf("\n");
}

//================================================
TCanvas *makeSlide(TH1 *h, const char* title, int nBadRuns, double nBadStatTot, double *badruns, double *badrunstat)
{
  TString cName = Form("%s_slide",h->GetName());
  TCanvas *c = new TCanvas(cName.Data(), cName.Data(), 800, 600);
  TPad *pad1 = new TPad(Form("%s_pad1",h->GetName()),Form("%s_pad1",h->GetName()),0.01,0.25,0.6,0.95);
  pad1->SetBorderMode(0);
  pad1->SetFillColor(kWhite);
  pad1->SetFrameLineWidth(0);
  pad1->SetBorderSize(0);
  pad1->SetTopMargin(0.1);
  pad1->SetLeftMargin(0.1);
  pad1->SetRightMargin(0.1);
  pad1->SetBottomMargin(0.1);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameFillColor(0);
  pad1->Draw();
  pad1->cd();
  TPaveText *t1 = GetTitleText(title,0.04,62);
  h->SetTitle("");
  h->DrawCopy("PE");
  t1->Draw();
  c->Update();
  c->cd();
  TPad *pad2 = new TPad(Form("%s_pad2",h->GetName()),Form("%s_pad2",h->GetName()),0.61,0.25,0.99,0.95);
  pad2->SetBorderMode(0);
  pad2->SetFillColor(kWhite);
  pad2->SetFrameLineWidth(0);
  pad2->SetBorderSize(0);
  pad2->SetTopMargin(0.1);
  pad2->SetLeftMargin(0.1);
  pad2->SetRightMargin(0.1);
  pad2->SetBottomMargin(0.1);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameFillColor(0);
  pad2->Draw();
  pad2->cd();
  int half = nBadRuns/2;
  TPaveText *t2 = GetPaveText(0.01, 0.5, 0.01, 0.99);
  for(int i=0; i<half; i++)
    {
      t2->AddText(Form("%1.0f(%4.1e)",badruns[i],badrunstat[i]));
    }
  t2->Draw();
  TPaveText *t3 = GetPaveText(0.5, 0.99, 0.01, 0.99);
  for(int i=half; i<nBadRuns; i++)
    {
      t3->AddText(Form("%1.0f(%4.1e)",badruns[i],badrunstat[i]));
    }
  t3->Draw();

  return c;
}



/*

  t1->AddText("Trigger: di-muon");

  TH1F *h1 = (TH1F*)fmc->Get("hAnalysisCuts");
  Int_t counter = (Int_t)(h1->GetBinContent(3)/1e4);
  Double_t vtxz = h1->GetBinContent(1)/counter;
  if(1)
    {
      t1->AddText(Form("Vertex cut: |vz_{TPC}| < %d cm, |vz_{TPC}-vz_{VPD}| < 3 cm, vr < 2 cm",100));
    }
  else
    {
      t1->AddText(Form("Vertex cut: none"));
    }

  t1->AddText(Form("Track selection:"));
  t1->AddText("    Primary tracks");
  //t1->AddText(Form("    p_{T} > %1.1f GeV/c",h1->GetBinContent(2)/counter));
  t1->AddText(Form("    p_{T} > 1 GeV/c"));
  t1->AddText(Form("    |#eta| < %1.1f",h1->GetBinContent(4)/counter));  
  t1->AddText(Form("    NHitsFit >= %1.0f",h1->GetBinContent(5)/counter));     
  t1->AddText(Form("    NHitsDedx >= %1.0f",h1->GetBinContent(6)/counter));
  t1->AddText(Form("    NHitsFit/NHitsPoss >= %1.2f",h1->GetBinContent(10)/counter));
  if(TMath::Abs(h1->GetBinContent(7)/counter)<10)
    t1->AddText(Form("    global dca <= %1.1f cm",h1->GetBinContent(7)/counter));  
  c1->Clear();
  title->Draw();
  t1->Draw();
  PaintCanvasToPDF(c1,pdf);

  title  = new TPaveText(0.1,0.85,0.9,0.95,"brNDC");
  title->SetFillStyle(0);
  title->SetBorderSize(0);
  title->SetTextFont(62);
  title->SetTextSize(0.07);
  title->AddText("PID cuts in data");
  
  t1  = new TPaveText(0.05,0.1,0.5,0.8,"brNDC");
  t1->SetTextAlign(11);
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextFont(62);
  t1->SetTextSize(0.045);
  t1->AddText("Pion PID");
  if(year==2013) 
    t1->AddText("    0.016 < m^{2} < 0.021 (GeV/c^{2})^{2}");
  t1->AddText("    |n#sigma_{#pi}| < 2");
  t1->AddText("Muon PID");
  t1->AddText(Form("    %1.1f < n#sigma_{#pi} < %1.1f",h1->GetBinContent(8)/counter,h1->GetBinContent(9)/counter));
  t1->AddText(Form("    |#Deltaz| < %1.0f cm",h1->GetBinContent(14)/counter));
  t1->AddText(Form("    Under J/#psi mass peak of [3.0,3.2] GeV/c^{2}"));
  c1->Clear();
  title->Draw();
  t1->Draw();
*/
