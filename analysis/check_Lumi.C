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
  AuAuMbData();
  //makeTacDiff();
  //ppTacSumDiff();
  //MbAuAuTacDiff();
}

//================================================
void AuAuMbData(const int savePlot = 0)
{
  // Au+Au MB
  TFile* fstudy = TFile::Open("output/Run14_AuAu200.MB.Study.root", "read");
  //TFile* fstudy = TFile::Open("jpsi.study.test.histos.root", "read");
  double nMbEvt[4][gNTrgSetup];
  double nMbEvtDm[4][gNTrgSetup];
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<gNTrgSetup; k++)
	{
	  nMbEvt[i][k]   = 0;
	  nMbEvtDm[i][k] = 0;
	}
    }

  THnSparseF *hnMbEvt = (THnSparseF*)fstudy->Get("mhMtdQaCentAll_mb");
  THnSparseF *hnMbEvtDm = (THnSparseF*)fstudy->Get("mhMtdQaCentAllDm_mb");

  TH2F *hMbEvt[3];
  TH2F *hMbEvtDm[3];
  hMbEvt[0] = (TH2F*)hnMbEvt->Projection(1,0);
  hMbEvt[1] = (TH2F*)fstudy->Get("mhMtdQaCent_mb");
  hMbEvt[2] = (TH2F*)fstudy->Get("mhMtdQaCentWeight_mb");

  hMbEvtDm[0] = (TH2F*)hnMbEvtDm->Projection(1,0);
  hMbEvtDm[1] = (TH2F*)fstudy->Get("mhMtdQaCentDm_mb");
  hMbEvtDm[2] = (TH2F*)fstudy->Get("mhMtdQaCentWeightDm_mb");

  TH1F *htmp = 0x0;
  for(int i=0; i<3; i++)
    {
      for(int k=0; k<gNTrgSetup; k++)
	{
	  if(k==0) htmp = (TH1F*)hMbEvt[i]->ProjectionY("",1,4);
	  else     htmp = (TH1F*)hMbEvt[i]->ProjectionY("",k,k);
	  nMbEvt[i][k] = htmp->Integral(1,16);

	  if(k==0) htmp = (TH1F*)hMbEvtDm[i]->ProjectionY("",1,4);
	  else     htmp = (TH1F*)hMbEvtDm[i]->ProjectionY("",k,k);
	  nMbEvtDm[i][k] = htmp->Integral(1,16);
	}
    }

  nMbEvt[3][0] = 64091730216;
  nMbEvt[3][1] = 4501422084;
  nMbEvt[3][2] = 7579287027;
  nMbEvt[3][3] = 16920976928;
  nMbEvt[3][4] = 35090044177;

  nMbEvtDm[3][0] = 1098238107/0.49;
  nMbEvtDm[3][1] = 98906160/0.49;
  nMbEvtDm[3][2] = 189965145/0.49;
  nMbEvtDm[3][3] = 292502539/0.49;
  nMbEvtDm[3][4] = 516864263/0.49;

  const char* typeName[4] = {"MB_all","MB_acc","MB_acc_w","DM_acc"};
  for(int k=0; k<gNTrgSetup; k++)
    {
      nMbEvtDm[2][k] = nMbEvtDm[1][k];
      printf("\n\n+++++++++++ AuAu_200%s +++++++++++\n",gTrgSetupTitle[k]);
      for(int i=0; i<4; i++)
	{
	  printf("[i] %s: %1.0f MB events or %1.0f DM events  = %4.3f%%\n",typeName[i],nMbEvt[i][k],nMbEvtDm[i][k],nMbEvtDm[i][k]/nMbEvt[i][k]*100);
	}
    }

  // Get DM fraction vs. ZDC rate
  TH1F *hMbZdc[2][gNTrgSetup-1];
  TH1F *hDmZdc[2][gNTrgSetup-1];
  TH1F *hDmToMbRatio[2][gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      hnMbEvt->GetAxis(0)->SetRange(k+1, k+1);
      hMbZdc[0][k] = (TH1F*)hnMbEvt->Projection(3);
      hMbZdc[0][k]->SetName(Form("hMbZdcAll%s",gTrgSetupTitle[k+1]));
      hnMbEvt->GetAxis(1)->SetRange(1,16);
      hnMbEvt->GetAxis(2)->SetRange(2,2);
      hMbZdc[1][k] = (TH1F*)hnMbEvt->Projection(3);
      hMbZdc[1][k]->SetName(Form("hMbZdcVtx%s",gTrgSetupTitle[k+1]));
      hnMbEvt->GetAxis(2)->SetRange(0,-1);
      hnMbEvt->GetAxis(0)->SetRange(0,-1);
      hnMbEvt->GetAxis(1)->SetRange(0,-1);

      hnMbEvtDm->GetAxis(0)->SetRange(k+1, k+1);
      hDmZdc[0][k] = (TH1F*)hnMbEvtDm->Projection(3);
      hDmZdc[0][k]->SetName(Form("hDmZdcAll%s",gTrgSetupTitle[k+1]));
      hnMbEvtDm->GetAxis(1)->SetRange(1,16);
      hnMbEvtDm->GetAxis(2)->SetRange(2,2);
      hDmZdc[1][k] = (TH1F*)hnMbEvtDm->Projection(3);
      hDmZdc[1][k]->SetName(Form("hDmZdcVtx%s",gTrgSetupTitle[k+1]));
      hnMbEvtDm->GetAxis(0)->SetRange(0,-1);
      hnMbEvtDm->GetAxis(1)->SetRange(0,-1);
      hnMbEvtDm->GetAxis(2)->SetRange(0,-1);
    }

  const char* typeName[2] = {"All", "Vtx"};
  const char* typeTitle[2] = {"All", "0-80%, vertex cut"};
  for(int i=0; i<2; i++)
    {
      TCanvas *c = new TCanvas(Form("cRatio_%s",typeName[i]), Form("cRatio_%s",typeName[i]), 800, 600);
      TLegend *leg = new TLegend(0.65,0.65,0.85,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hDmToMbRatio[i][k] = (TH1F*)hDmZdc[i][k]->Clone(Form("hDmToMbRatio%s_%d",gTrgSetupTitle[k+1],i));
	  hDmToMbRatio[i][k]->Sumw2();
	  hDmToMbRatio[i][k]->Divide(hMbZdc[i][k]);
	  hDmToMbRatio[i][k]->SetMarkerStyle(20+k);
	  hDmToMbRatio[i][k]->SetMarkerColor(gColor[k+1]);
	  hDmToMbRatio[i][k]->SetLineColor(gColor[k+1]);
	  hDmToMbRatio[i][k]->SetTitle("");
	  hDmToMbRatio[i][k]->GetYaxis()->SetRangeUser(0, 0.07);
	  leg->AddEntry(hDmToMbRatio[i][k], gTrgSetupTitle[k+1], "P");
	  if(k==0) hDmToMbRatio[i][k]->Draw();
	  else     hDmToMbRatio[i][k]->Draw("sames");
	}
      TPaveText *t1 = GetTitleText(Form("DM/MB ratio (%s)",typeTitle[i]),0.04);
      t1->Draw();
      leg->Draw();
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_DmToMbVsZdc_%s.pdf",run_type.Data(),typeName[i]));
    }

  // MT101 distribution
  const double min_TacDiffCut[4] = {785+1,785+1,788+1,788+1};
  const double max_TacDiffCut[4] = {837, 837, 837, 837};
  TFile *fTrig = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
  TH3F *mhMT101TacDiff = (TH3F*)fTrig->Get("mhMT101TacDiff_di_mu");
  TH1F *hMT101TacDiff[gNTrgSetup-1][9];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      TCanvas *c = new TCanvas(Form("cMT101_%s",gTrgSetupTitle[k+1]), Form("cMT101_%s",gTrgSetupTitle[k+1]), 1100, 700);
      c->Divide(3,3);
      for(int i=0; i<9; i++)
	{
	  c->cd(i+1);
	  //gPad->SetLogy();
	  hMT101TacDiff[k][i] = (TH1F*)mhMT101TacDiff->ProjectionX(Form("hMT101TacDiff%s_Zdc%d",gTrgSetupTitle[k+1],i),k+1,k+1,3+i*2,4+i*2);
	  hMT101TacDiff[k][i]->SetTitle("");
	  if(hMT101TacDiff[k][i]->GetEntries()<=0) continue;
	  hMT101TacDiff[k][i]->Scale(1./hMT101TacDiff[k][i]->Integral(0,-1));
	  hMT101TacDiff[k][i]->GetXaxis()->SetRangeUser(750,820);
	  hMT101TacDiff[k][i]->Draw();
	  double eff = hMT101TacDiff[k][i]->Integral(hMT101TacDiff[k][i]->FindBin(min_TacDiffCut[k]),
						     hMT101TacDiff[k][i]->FindBin(max_TacDiffCut[k]));
	  TPaveText *t1 = GetPaveText(0.2,0.4,0.6,0.7,0.06);
	  t1->AddText(Form("eff = %4.2f%%",eff*100));
	  t1->Draw();
	  TPaveText *t1 = GetTitleText(Form("%d < ZDC < %d kHz (%s)",10+i*10,20+i*10,gTrgSetupTitle[k+1]),0.06);
	  t1->Draw();

	  TLine *line1 = GetLine(min_TacDiffCut[k], 0, min_TacDiffCut[k], hMT101TacDiff[k][i]->GetMaximum(), 2, 1);
	  line1->Draw();
	  TLine *line2 = GetLine(max_TacDiffCut[k], 0, max_TacDiffCut[k], hMT101TacDiff[k][i]->GetMaximum(), 2, 1);
	  line2->Draw();
	}
      if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuMb_TacSumDis%s.pdf",run_type.Data(),gTrgSetupTitle[k+1]));
    }
}


//================================================
void MbAuAuTacDiff(const int savePlot = 0)
{
  TFile *fin[2];
  fin[0] = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
  fin[1] = TFile::Open("output/Run14_AuAu200.Study.MtdTrig.root", "read");
  const char* trigName[2] = {"mb", "di_mu"};
  const char* trigTitle[2] = {"Vpd-NoVtx", "dimuon"};
  TH2F *hTacDiffVsLumi[2];
  TH3F *hTacDiffVsPtVsLumiMth[2];
  TH3F *hTacDiffVsPtVsLumiMuon[2];
  for(int i=0; i<2; i++)
    {
      hTacDiffVsLumi[i] = (TH2F*)fin[i]->Get(Form("mhMT101TacDiff_di_mu"));
      hTacDiffVsLumi[i]->SetName(Form("mhMT101TacDiff_%s",trigName[i]));
      hTacDiffVsLumi[i]->Sumw2();

      hTacDiffVsPtVsLumiMth[i] = (TH3F*)fin[i]->Get(Form("mhMT101TacDiffMatched_di_mu"));
      hTacDiffVsPtVsLumiMth[i]->SetName(Form("mhMT101TacDiffMatched_%s",trigName[i]));
      hTacDiffVsPtVsLumiMth[i]->Sumw2();

      hTacDiffVsPtVsLumiMuon[i] = (TH3F*)fin[i]->Get(Form("mhMT101TacDiffMuon_di_mu"));
      hTacDiffVsPtVsLumiMuon[i]->SetName(Form("mhMT101TacDiffMuon_%s",trigName[i]));
      hTacDiffVsPtVsLumiMuon[i]->Sumw2();
    }
  
  // MtdVpdTacDiff vs. luminosity
  TH1F *hTacDiffInLumi[2][gNTrgSetup-1];
  TCanvas *cTacVsLumi[2];
  TH1F *hplot = new TH1F("hplot", ";#DeltaTacSum;", 1024, 0, 1024);
  for(int i=0; i<2; i++)
    {
      cTacVsLumi[i] = new TCanvas(Form("cTacVsLumi_%s",trigName[i]),Form("cTacVsLumi_%s",trigName[i]),800,600);
      hplot->GetXaxis()->SetRangeUser(650,840);
      hplot->GetYaxis()->SetRangeUser(1e-5, 1.2);
      hplot->DrawCopy();
      TLegend *leg = new TLegend(0.2,0.55,0.4,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->SetHeader("All hits");
      for(int k=0; k<gNTrgSetup-1; k++)
	{
	  hTacDiffInLumi[i][k] = (TH1F*)hTacDiffVsLumi[i]->ProjectionX(Form("hTacDiff_%s%s",trigName[i],gTrgSetupTitle[k+1]),k+1,k+1);
	  hTacDiffInLumi[i][k]->Scale(1./hTacDiffInLumi[i][k]->GetBinContent(hTacDiffInLumi[i][k]->FindBin(790.5)));
	  hTacDiffInLumi[i][k]->SetMarkerStyle(20+k);
	  hTacDiffInLumi[i][k]->SetMarkerColor(gColor[k]);
	  hTacDiffInLumi[i][k]->SetLineColor(gColor[k]);
	  hTacDiffInLumi[i][k]->DrawCopy("sames");
	  leg->AddEntry(hTacDiffInLumi[i][k],Form("Run14%s",gTrgSetupTitle[k+1]),"P");
	}
      TPaveText *t1 = GetTitleText(Form("%s: %s",run_type.Data(),trigTitle[i]),0.045);
      t1->Draw();
      leg->Draw();
      if(savePlot) cTacVsLumi[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuTacDiffInLumi_%s.pdf",run_type.Data(),trigName[i]));
    }

  // MtdVpdTacDiff for different type of tracks
  TCanvas *cTacVsTrk = new TCanvas("cTacVsTrk", "cTacVsTrk", 1100, 700);
  cTacVsTrk->Divide(2,2);
  TLegend *leg = new TLegend(0.15,0.55,0.35,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(trigTitle[0]);
  TH1F *htmp = 0x0;
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      cTacVsTrk->cd(k+1);
      hplot->GetXaxis()->SetRangeUser(700,830);
      hplot->GetYaxis()->SetRangeUser(1e-5, 0.06);
      hplot->DrawCopy();
      for(int i=0; i<3; i++)
	{
	  if(i==0) htmp = (TH1F*)hTacDiffInLumi[0][k]->Clone(Form("hTacDiffAll_%s%s",trigName[0],gTrgSetupTitle[k+1]));
	  if(i==1) htmp = (TH1F*)hTacDiffVsPtVsLumiMth[0]->ProjectionX(Form("hTacDiffMth_%s%s",trigName[0],gTrgSetupTitle[k+1]),k+1,k+1,1,100);
	  if(i==2) htmp = (TH1F*)hTacDiffVsPtVsLumiMuon[0]->ProjectionX(Form("hTacDiffMuon_%s%s",trigName[0],gTrgSetupTitle[k+1]),k+1,k+1,1,100);
	  htmp->Scale(1./htmp->Integral());
	  htmp->SetMarkerStyle(22+i);
	  htmp->SetMarkerColor(gColor[i]);
	  htmp->SetLineColor(gColor[i]);
	  htmp->Draw("sames");
	  if(k==0) 
	    {
	      if(i==0) leg->AddEntry(htmp, "All hits", "P");
	      if(i==1) leg->AddEntry(htmp, "Mathed to MTD", "P");
	      if(i==2) leg->AddEntry(htmp, "Muon candidates", "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s%s",run_type.Data(),gTrgSetupTitle[k+1]),0.06);
      t1->Draw();
    }
  cTacVsTrk->cd(1);
  leg->Draw();
  if(savePlot) cTacVsTrk->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuTacDiff_MthVsMuon_%s.pdf",run_type.Data(),trigName[0]));
  
  // MtdVpdTacDiff vs. pt for matched tracks
  const int nPtBins = 6;
  const double lowPtBin[nPtBins] = {1.3, 1.5, 2.0, 2.5, 3.0, 5.0};
  const double upPtBin[nPtBins]  = {1.5, 2.0, 2.5, 3.0, 5.0, 8.0};
  TH1F *hTacDiffInPt[2][gNTrgSetup-1][nPtBins];
  TF1 *funcTacDiffInPt[2][gNTrgSetup-1][nPtBins];
  TGraphErrors *gTacDiffMeanVsPt[2][gNTrgSetup-1];
  TCanvas *cTacVsPt[gNTrgSetup-1];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      cTacVsPt[k] = new TCanvas(Form("cTacVsPt%s",gTrgSetupTitle[k+1]),Form("cTacVsPt%s",gTrgSetupTitle[k+1]),1100,700);
      cTacVsPt[k]->Divide(4,2);
      TLegend *leg = new TLegend(0.15,0.5,0.5,0.8);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.07);
      leg->SetHeader(Form("Run14%s",gTrgSetupTitle[k+1]));
      for(int i=0; i<2; i++)
	{
	  gTacDiffMeanVsPt[i][k] = new TGraphErrors(nPtBins);
	}
      for(int p=0; p<nPtBins; p++)
	{
	  cTacVsPt[k]->cd(p+2);
	  hplot->GetXaxis()->SetRangeUser(750,840);
	  hplot->GetYaxis()->SetRangeUser(1e-5, 1.2);
	  hplot->DrawCopy();
	  for(int i=0; i<2; i++)
	    {
	      int lowbin = hTacDiffVsPtVsLumiMuon[i]->GetZaxis()->FindFixBin(lowPtBin[p]+1e-4);
	      int upbin  = hTacDiffVsPtVsLumiMuon[i]->GetZaxis()->FindFixBin(upPtBin[p]-1e-4);
	      hTacDiffInPt[i][k][p] = (TH1F*)hTacDiffVsPtVsLumiMuon[i]->ProjectionX(Form("hTacDiff_pt%d_%s%s",p+1,trigName[i],gTrgSetupTitle[k+1]),k+1,k+1,lowbin,upbin);
	      if(i==0 && p>2) hTacDiffInPt[i][k][p]->Rebin(4);
	      hTacDiffInPt[i][k][p]->Scale(1./hTacDiffInPt[i][k][p]->GetBinContent(hTacDiffInPt[i][k][p]->FindBin(790.5)));
	      hTacDiffInPt[i][k][p]->Scale(1./hTacDiffInPt[i][k][p]->Integral());
	      funcTacDiffInPt[i][k][p] = new TF1(Form("fit_%s",hTacDiffInPt[i][k][p]->GetName()), "gaus", 790, 810);
	      if(k<2) funcTacDiffInPt[i][k][p]->SetRange(786,810);
	      else if(k<3) funcTacDiffInPt[i][k][p]->SetRange(789,810);
	      hTacDiffInPt[i][k][p]->Fit(funcTacDiffInPt[i][k][p], "IR0Q");
	      double scale = 1./hTacDiffInPt[i][k][p]->GetBinContent(hTacDiffInPt[i][k][p]->FindBin(funcTacDiffInPt[i][k][p]->GetParameter(1)));
	      hTacDiffInPt[i][k][p]->Scale(scale);
	      hTacDiffInPt[i][k][p]->SetMarkerStyle(20+5*i);
	      hTacDiffInPt[i][k][p]->SetMarkerColor(gColor[i]);
	      hTacDiffInPt[i][k][p]->SetLineColor(gColor[i]);
	      hTacDiffInPt[i][k][p]->DrawCopy("sames");
	      funcTacDiffInPt[i][k][p]->SetLineColor(gColor[i+2]);
	      funcTacDiffInPt[i][k][p]->SetLineStyle(2);
	      funcTacDiffInPt[i][k][p]->SetParameter(0, scale*funcTacDiffInPt[i][k][p]->GetParameter(0));
	      funcTacDiffInPt[i][k][p]->Draw("sames");
	      gTacDiffMeanVsPt[i][k]->SetPoint(p, (lowPtBin[p]+upPtBin[p])/2., funcTacDiffInPt[i][k][p]->GetParameter(1));
	      gTacDiffMeanVsPt[i][k]->SetPointError(p, (upPtBin[p]-lowPtBin[p])/2., funcTacDiffInPt[i][k][p]->GetParError(1));
	      if(p==0)
		{
		  leg->AddEntry(hTacDiffInPt[i][k][p], trigTitle[i], "P");
		}
	    }
	  TPaveText *t1 = GetTitleText(Form("%1.1f < p_{T} < %1.1f GeV/c",lowPtBin[p],upPtBin[p]),0.07);
	  t1->Draw();
	}
      cTacVsPt[k]->cd(8);
      htmp = (TH1F*)hplot->Clone("hplot_clone");
      htmp->GetXaxis()->SetRangeUser(0,9);
      htmp->GetYaxis()->SetRangeUser(790,802);
      htmp->SetTitle(";p_{T} (GeV/c);<#DeltaTacSum>");
      htmp->GetYaxis()->SetTitleOffset(1.3);
      htmp->DrawCopy();
      for(int i=0; i<2; i++)
	{
	  gTacDiffMeanVsPt[i][k]->SetMarkerStyle(20+5*i);
	  gTacDiffMeanVsPt[i][k]->SetMarkerColor(gColor[i]);
	  gTacDiffMeanVsPt[i][k]->SetLineColor(gColor[i]);
	  gTacDiffMeanVsPt[i][k]->Draw("samesPEZ");
	}
      cTacVsPt[k]->cd(1);
      leg->Draw();
      if(savePlot) cTacVsPt[k]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuTacDiffFit%s.pdf",run_type.Data(),gTrgSetupTitle[k+1]));
    }

  // compare efficiency to the one estimated using pp results
  const double min_TacDiffCut[gNTrgSetup-1] = {785+1,785+1,788+1,788+1};
  const double max_TacDiffCut[gNTrgSetup-1] = {837, 837, 837, 837};
  TGraphAsymmErrors *gTacDiffEff[gNTrgSetup-1][2];
  double pt_arr[nPtBins], all_arr[nPtBins], all_err_arr[nPtBins], acc_arr[nPtBins], acc_err_arr[nPtBins];
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      for(int p=0; p<nPtBins; p++)
	{
	  int lowbin = hTacDiffVsPtVsLumiMuon[0]->GetZaxis()->FindFixBin(lowPtBin[p]+1e-4);
	  int upbin  = hTacDiffVsPtVsLumiMuon[0]->GetZaxis()->FindFixBin(upPtBin[p]-1e-4);
	  htmp = (TH1F*)hTacDiffVsPtVsLumiMuon[0]->ProjectionX(Form("hTacDiff_pt%d_%s%s_tmp",p+1,trigName[0],gTrgSetupTitle[k+1]),k+1,k+1,lowbin,upbin);
	  pt_arr[p] = (lowPtBin[p]+upPtBin[p])/2.;
	  all_arr[p] = htmp->IntegralAndError(1, htmp->GetNbinsX(), all_err_arr[p]);
	  acc_arr[p] = htmp->IntegralAndError(htmp->FindFixBin(min_TacDiffCut[k]+0.1),
					      htmp->FindFixBin(max_TacDiffCut[k]-0.1),
					      acc_err_arr[p]);
	  //cout << pt_arr[p] << "  " << all_arr[p] << "  " << all_err_arr[p] << "  " << acc_arr[p] << "  " << acc_err_arr[p] << endl;
	  if(k==2 && p==-1)
	    {
	      htmp->Rebin(4);
	      htmp->Scale(1./htmp->GetBinContent(htmp->FindBin(funcTacDiffInPt[0][k][p]->GetParameter(1))));
	      htmp->GetXaxis()->SetRangeUser(750, 840);
	      htmp->SetMarkerStyle(20);
	      c = draw1D(htmp);
	      TF1 *functmp = new TF1(Form("fitCB_%s",htmp->GetName()), CrystalBall, 750, 840, 5);
	      functmp->SetParameters(1,funcTacDiffInPt[0][k][p]->GetParameter(1),funcTacDiffInPt[0][k][p]->GetParameter(2),0.2, 1);
	      htmp->Fit(functmp,"IR0L");
	      functmp->SetLineColor(2);
	      functmp->Draw("sames");
	      double rel_err_all = all_err_arr[p]/all_arr[p];
	      double rel_err_acc = acc_err_arr[p]/acc_arr[p];
	      all_arr[p] = functmp->Integral(600, 1000);
	      all_err_arr[p] = all_arr[p]*rel_err_all;
	      acc_arr[p] = functmp->Integral(min_TacDiffCut[k], max_TacDiffCut[k]);
	      acc_err_arr[p] = acc_arr[p]*rel_err_acc;
	    }
	}
      gTacDiffEff[k][0] = GetEfficiencyCurve(nPtBins, pt_arr, all_arr, all_err_arr, acc_arr, acc_err_arr);
      gTacDiffEff[k][0]->SetName(Form("%s_gTacDiffEff_LS%s_AuAuMB",run_type.Data(),gTrgSetupTitle[k+1]));
    }

  TFile *fdep = TFile::Open("Rootfiles/old.Run14_AuAu200.StudyLumiDep.root", "read");
  TF1 *funcTacEff[gNTrgSetup-1][2];
  double x, y, x1, y1;
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      gTacDiffEff[k][1] = (TGraphAsymmErrors*)fdep->Get(Form("%s_gTacDiffEff_LS%s",run_type.Data(),gTrgSetupTitle[k+1]));
      funcTacEff[k][1] = (TF1*)fdep->Get(Form("%s_gTacDiffEff_LS%s_func",run_type.Data(),gTrgSetupTitle[k+1]));

      for(int p=0; p<nPtBins; p++)
	{
	  gTacDiffEff[k][1]->GetPoint(p, x, y);

	  gTacDiffEff[k][0]->GetPoint(p, x1, y1);
	  gTacDiffEff[k][0]->SetPoint(p, x, y1);
	  gTacDiffEff[k][0]->SetPointEXhigh(p, gTacDiffEff[k][1]->GetErrorXhigh(p));
	  gTacDiffEff[k][0]->SetPointEXlow(p, gTacDiffEff[k][1]->GetErrorXlow(p));
	}
      funcTacEff[k][0] = new TF1(Form("%s_func",gTacDiffEff[k][0]->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1.3,7);
      funcTacEff[k][0]->SetParameters(0.9, 3.5, 0.8);
      TGraphAsymmErrors *gFitTmp = (TGraphAsymmErrors*)gTacDiffEff[k][0]->Clone(Form("%s_clone",gTacDiffEff[k][0]->GetName()));
      if(k==-1)
	{
	  gFitTmp->SetPointEYhigh(3, 0.2);
	  gFitTmp->SetPointEYlow(3, 0.2);
	}
      gFitTmp->Fit(funcTacEff[k][0],"IR0");
    }
  TCanvas *cTacEff = new TCanvas("cTacEff", "cTacEff", 1100, 700);
  cTacEff->Divide(2,2);
  TLegend *leg = new TLegend(0.5,0.3,0.8,0.5);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  for(int k=0; k<gNTrgSetup-1; k++)
    {
      cTacEff->cd(k+1);
      hplot->GetXaxis()->SetRangeUser(0,9);
      hplot->GetYaxis()->SetRangeUser(0.6,1.0);
      hplot->SetTitle(";p_{T} (GeV/c);Eff.");
      hplot->DrawCopy();
      for(int i=0; i<2; i++)
	{
	  gTacDiffEff[k][i]->SetMarkerStyle(20+i*4);
	  gTacDiffEff[k][i]->SetMarkerColor(2-i);
	  gTacDiffEff[k][i]->SetLineColor(2-i);
	  gTacDiffEff[k][i]->Draw("samesPEZ");
	  funcTacEff[k][i]->SetLineColor(2-i);
	  funcTacEff[k][i]->SetLineStyle(2);
	  funcTacEff[k][i]->Draw("sames");
	  if(k==0)
	    {
	      if(i==0) leg->AddEntry(gTacDiffEff[k][i], "2014 Au+Au MB", "P");
	      if(i==1) leg->AddEntry(gTacDiffEff[k][i], "2015 p+p dimuon", "P");
	    }
	}
      TPaveText *t1 = GetTitleText(Form("%s%s: trigger efficiency for muon candidates",run_type.Data(),gTrgSetupTitle[k+1]),0.05);
      t1->Draw();
    }
  cTacEff->cd(1);
  leg->Draw();
  if(savePlot) cTacEff->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_AuAuMBVspp.pdf",run_type.Data()));
}

//================================================
void ppTacSumDiff(const int savePlot = 1)
{
  // compare dTacSum from pp
  TH1F *hppTacDiff[3];

  // MT101
  TFile *fin = TFile::Open("output/Run15_pp200.Study.MtdTrig.root", "read");
  TH2F *hMT101TacDiff = (TH2F*)fin->Get("mhMT101TacDiff_di-muon");
  hppTacDiff[0] = (TH1F*)hMT101TacDiff->ProjectionY("hppTacDiff_MT101");
  hppTacDiff[0]->Sumw2();

  // matched tracks
  TFile *fpp = TFile::Open("Rootfiles/Run15_pp200.MtdTrigEff.root", "read");
  TH2F *hppTacDiffLS = (TH2F*)fpp->Get("Run15_pp200_hTacDiffVsTrigUnit_LS_PtBin0");
  hppTacDiff[1] = (TH1F*)hppTacDiffLS->ProjectionY("hppTacDiff_LS");

  TH2F *hppTacDiffMuon = (TH2F*)fpp->Get("Run15_pp200_hTacDiffVsTrigUnit_Muon_PtBin0");
  hppTacDiff[2] = (TH1F*)hppTacDiffMuon->ProjectionY("hppTacDiff_Muon");

  // compare
  const int peaks[3] = {920, 927, 930};
  TString legName[3] = {"MT101", "Like-sign", "Muon"};
  TList *list = new TList;
  for(int i=0; i<3; i++)
    {
      hppTacDiff[i]->Scale(1./hppTacDiff[i]->GetBinContent(hppTacDiff[i]->FindFixBin(peaks[i])));
      hppTacDiff[i]->SetMarkerStyle(20+i);
      hppTacDiff[i]->SetMarkerColor(gColor[i]);
      hppTacDiff[i]->SetLineColor(gColor[i]);
      hppTacDiff[i]->GetXaxis()->SetRangeUser(875, 1005);
      hppTacDiff[i]->GetYaxis()->SetRangeUser(0, 1.2);
      list->Add(hppTacDiff[i]);
    }
  TCanvas *c = drawHistos(list,"ppTacSumDiff","Run15_pp200: #DeltaTacSum distribution",true,875,1005,true,0,1.2,kFALSE,kTRUE,legName,true,"p_{T} > 1.3 GeV/c",0.6,0.8,0.65,0.85,true,0.04,0.04,false,0,false,false);
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_ppTacSumComp.pdf",run_type.Data()));

}

//================================================
void makeTacDiff(const int savePlot = 0, const int saveHisto= 0)
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
      //if(bin%1000==1) printf("[i] Process run %1.0f\n",run);
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
      //if(run<=15164050 && run>=15164020)
      if(run<=15109032 && run>=15109006)
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
	  cFit[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_FitTacSumDiffPerRun_%d.pdf",run_type.Data(),i));
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
void equiMB(const bool savePlot = 0)
{
  const int nCentBins       = nCentBins_pt; 
  const int* centBins_low   = centBins_low_pt;
  const int* centBins_high  = centBins_high_pt;
  const char** cent_Name    = cent_Name_pt;
  const char** cent_Title   = cent_Title_pt;

  // =============================================
  TFile *fdata = TFile::Open(Form("./output/Run14_AuAu200.jpsi.%sroot",run_config),"read");
  TH1F *hEvtRun = (TH1F*)fdata->Get("mhEvtRun_di_mu");
  TH1F *hEvtRunAcc = (TH1F*)fdata->Get("mhEvtRunAcc_di_mu");

  // =============================================
  TFile *fmb = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root", "read");
  TH1F *hNEventsAll = (TH1F*)fmb->Get("mhNEventsAll");
  TH1F *hNEventsDM  = (TH1F*)fmb->Get("mhNEventsDM");
  hNEventsAll->Sumw2();
  hNEventsDM->Sumw2();
  TH1F *hNEventsRatio = (TH1F*)hNEventsAll->Clone("hNEventsRatio");
  hNEventsRatio->Divide(hNEventsDM);

  // =============================================
  TFile *fLumi = TFile::Open(Form("Rootfiles/Run14_AuAu200.Luminosity.root"),"read");
  TH1F *hRF = (TH1F*)fLumi->Get("hRejectFactor_dimuon");
  TH1F *hNeventsTake = (TH1F*)fLumi->Get("hNevents_dimuon");
  TH1F *hPSdimuon = (TH1F*)fLumi->Get("hPreScale_dimuon");
  TH1F *hLTdimuon = (TH1F*)fLumi->Get("hLiveTime_dimuon");
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
  TH1F *hEqvRatioInTrigMB[4];
  TH1F *hEqvRatioInTrigGood[4];
  TH1F *hEvtRFCordimuon[4];

  double mb_event_auau[gNTrgSetup];
  double mtd_event_auau[gNTrgSetup];
  double dimuonevent[gNTrgSetup];

  mb_event_auau[0] = 0;
  mtd_event_auau[0] = 0;
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
      hEqvRatioInTrigMB[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEqvRatioInTrigMB_%s",trgSetupName[j-1]));
      hEqvRatioInTrigMB[j-1]->Reset("AC");
      hEqvRatioInTrigGood[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEqvRatioInTrigGood_%s",trgSetupName[j-1]));
      hEqvRatioInTrigGood[j-1]->Reset("AC");
      hEvtRFCordimuon[j-1] = (TH1F*)hNeventsTakeMb->Clone(Form("hEvtRFCordimuon_%s",trgSetupName[j-1]));
      hEvtRFCordimuon[j-1]->Reset("AC");

      mb_event_auau[j] = 0;
      mtd_event_auau[j] = 0;
      dimuonevent[j] = 0;
      ifstream fruns;
      fruns.open(Form("Rootfiles/Luminosity/%s/AuAu_200_production%s_2014.list",run_type.Data(),trgSetupName[j-1]));
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

	  mb_event_auau[j] += hNEventsAll->GetBinContent(hNEventsAll->FindBin(runnumber));
	  mb_event_auau[0] += hNEventsAll->GetBinContent(hNEventsAll->FindBin(runnumber));

	  mtd_event_auau[j]+= hNEventsDM->GetBinContent(hNEventsDM->FindBin(runnumber));
	  mtd_event_auau[0]+= hNEventsDM->GetBinContent(hNEventsDM->FindBin(runnumber));

	  if(nEventsTaken==0) 
	    {
	      printf("[w] check run %1.0f\n",runnumber);
	      continue;
	    }
	  double nEventsRun = hEvtRun->GetBinContent(bin);
	  double rf = hRF->GetBinContent(hRF->FindFixBin(runnumber));
	  if(rf==0)
	    {
	      printf("[w] rf = 0 for run %1.0f\n",runnumber);
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
	  hEqvRatioInTrigMB[j-1]->SetBinContent(lumiBin, hNEventsRatio->GetBinContent(hNEventsRatio->FindBin(runnumber))*1./hPSdimuon->GetBinContent(lumiBin));
	  hEqvRatioInTrigMB[j-1]->SetBinError(lumiBin, hNEventsRatio->GetBinError(hNEventsRatio->FindBin(runnumber))*1./hPSdimuon->GetBinContent(lumiBin));
	  hEqvRatioInTrigGood[j-1]->SetBinContent(lumiBin, nEventsRun/rf * eq_mb_good/nEventsTaken/nEventsTaken);
	  //hEqvRatioInTrigGood[j-1]->SetBinContent(lumiBin, eq_mb_good/nEventsTaken);
	  hEqvRatioInTrigGood[j-1]->SetBinError(lumiBin, 1e-10);
	  hEvtRFCordimuon[j-1]->SetBinContent(lumiBin, nEventsRun/rf/nEventsTaken);
	  hEvtRFCordimuon[j-1]->SetBinError(lumiBin, 1e-10);
	  if(j==2 && runnumber>=15165000)
	    {
	      printf("[i] run = %1.0f; dimuon: ps = %4.4f, lt = %4.4f, nEventsTaken = %1.1f; mb: ps = %4.4f, lt = %4.4f, nEventsTaken = %1.1f, nEventsAll = %1.1f, nEventsDM = %1.1f; ratio1 = %4.4f, ratio2 = %4.4f\n",
		     runnumber,hPSdimuon->GetBinContent(lumiBin), hLTdimuon->GetBinContent(lumiBin), nEventsTaken,
		     hPSmb->GetBinContent(lumiBin), hLTmb->GetBinContent(lumiBin), hNeventsTakeMb->GetBinContent(lumiBin), 
		     hNEventsAll->GetBinContent(hNEventsAll->FindBin(runnumber)), hNEventsDM->GetBinContent(hNEventsDM->FindBin(runnumber)),
		     hEqvRatioInTrig[j-1]->GetBinContent(lumiBin), hEqvRatioInTrigMB[j-1]->GetBinContent(lumiBin));
	      //printf("[i] run = %1.0f, nEventsRun = %1.1f, rf = %4.2f, nEventsTaken = %1.1f\n",runnumber,nEventsRun,rf,nEventsTaken);
	    }
	}
    }
  printf("+++++++++++++++++++++++++++++++++\n");

  
  printf("\n\n++++++++++++ Au+Au MB +++++++++++++\n");
  for(int j=0; j<gNTrgSetup; j++)
    {
      printf("[i] %s: %1.0f MB events, %1.0f dimuon events (%4.3f%%)\n",gTrgSetupTitle[j],mb_event_auau[j],mtd_event_auau[j],mtd_event_auau[j]/mb_event_auau[j]*100);
    }
  printf("+++++++++++++++++++++++++++++++++\n");
  // =============================================


  // check RF corrected dimuon events
  TCanvas *c = new TCanvas("Comp_EvtRFCordimuon","Comp_EvtRFCordimuon",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEvtRFCordimuon[j-1]->SetMarkerStyle(19+j);
      hEvtRFCordimuon[j-1]->SetMarkerColor(gColor[j-1]);
      hEvtRFCordimuon[j-1]->SetLineColor(gColor[j-1]);
      hEvtRFCordimuon[j-1]->GetYaxis()->SetRangeUser(0,1.5);
      if(j==1) 
	{
	  hEvtRFCordimuon[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("(# of analyzed dimuon events * RF)/(# of recorded dimuon events)",0.035);
	  t1->Draw();
	}
      else     hEvtRFCordimuon[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_DimuonRFcorrEvtvsRun.pdf",run_type.Data()));
  //return;

  // check the pre-scale and live-time
  TCanvas *c = new TCanvas("check_ps_lt", "check_ps_lt", 1200, 700);
  c->Divide(3,2);
  TLegend *leg = new TLegend(0.6,0.15,0.8,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04); 

  for(int j=1; j<gNTrgSetup; j++)
    {
      //if(j>=3) continue;
      c->cd(1);
      hPSdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hPSdimuonInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hPSdimuonInTrig[j-1]->SetLineColor(gColor[j-1]);
      hPSdimuonInTrig[j-1]->SetTitle("Dimuon: pre-scale;runId;");
      hPSdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5);
      if(j==1) hPSdimuonInTrig[j-1]->Draw("P");
      else     hPSdimuonInTrig[j-1]->Draw("samesP");
      leg->AddEntry(hPSdimuonInTrig[j-1], Form("prod%s",trgSetupName[j-1]), "P");

      c->cd(2);
      hNdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hNdimuonInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hNdimuonInTrig[j-1]->SetLineColor(gColor[j-1]);
      hNdimuonInTrig[j-1]->SetTitle("# of dimuon events on tape;runId;");
      hNdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e6);
      if(j==1) hNdimuonInTrig[j-1]->Draw("P");
      else     hNdimuonInTrig[j-1]->Draw("samesP");

      c->cd(3);
      hLTdimuonInTrig[j-1]->SetMarkerStyle(19+j);
      hLTdimuonInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hLTdimuonInTrig[j-1]->SetLineColor(gColor[j-1]);
      hLTdimuonInTrig[j-1]->SetTitle("Dimuon live-time");
      hLTdimuonInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTdimuonInTrig[j-1]->Draw("P");
      else     hLTdimuonInTrig[j-1]->Draw("samesP");

      c->cd(4);
      hPSmbInTrig[j-1]->SetMarkerStyle(19+j);
      hPSmbInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hPSmbInTrig[j-1]->SetLineColor(gColor[j-1]);
      hPSmbInTrig[j-1]->SetTitle("MB: pre-scale;runId;");
      hPSmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,3e4);
      if(j==1) hPSmbInTrig[j-1]->Draw("P");
      else     hPSmbInTrig[j-1]->Draw("samesP");

      c->cd(5);
      hNmbInTrig[j-1]->SetMarkerStyle(19+j);
      hNmbInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hNmbInTrig[j-1]->SetLineColor(gColor[j-1]);
      hNmbInTrig[j-1]->SetTitle("# of MB events on tape;runId;");
      hNmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1.5e4);
      if(j==1) hNmbInTrig[j-1]->Draw("P");
      else     hNmbInTrig[j-1]->Draw("samesP");

      c->cd(6);
      hLTmbInTrig[j-1]->SetMarkerStyle(19+j);
      hLTmbInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hLTmbInTrig[j-1]->SetLineColor(gColor[j-1]);
      hLTmbInTrig[j-1]->SetTitle("MB live-time");
      hLTmbInTrig[j-1]->GetYaxis()->SetRangeUser(0,1);
      if(j==1) hLTmbInTrig[j-1]->Draw("P");
      else     hLTmbInTrig[j-1]->Draw("samesP");
    }
  c->cd(1);
  leg->Draw();
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_PSandLT.pdf",run_type.Data()));
  return;

  // const int firstrun = 15109000;
  // const int lastrun  = 15109040;
  const int firstrun = 0;
  const int lastrun = 1e8;
  TCanvas *c = new TCanvas("Comp_EqMbPerTrig","Comp_EqMbPerTrig",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrig[j-1]->SetMarkerStyle(19+j);
      hEqvRatioInTrig[j-1]->SetMarkerColor(gColor[j-1]);
      hEqvRatioInTrig[j-1]->SetLineColor(gColor[j-1]);
      hEqvRatioInTrig[j-1]->SetTitle(";runId;");
      hEqvRatioInTrig[j-1]->GetYaxis()->SetRangeUser(0,70);
      hEqvRatioInTrig[j-1]->GetXaxis()->SetRangeUser(firstrun, lastrun);
      if(j==1) 
	{
	  hEqvRatioInTrig[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB events on tape per dimuon event",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrig[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBvsRun.pdf",run_type.Data()));

  TCanvas *c = new TCanvas("Comp_EqMbPerTrig_MB","Comp_EqMbPerTrig_MB",1100,700);
  c->Divide(2,2);
  for(int j=1; j<gNTrgSetup; j++)
    {
      c->cd(j);
      TH1F *hRatio = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_ratio",hEqvRatioInTrig[j-1]->GetName()));
      hRatio->Divide(hEqvRatioInTrigMB[j-1]);
      hRatio->GetYaxis()->SetRangeUser(0.1,2);
      hRatio->Draw("P");
      // hEqvRatioInTrig[j-1]->Draw("P");
      // hEqvRatioInTrigMB[j-1]->SetMarkerStyle(24);
      // hEqvRatioInTrigMB[j-1]->Draw("samesP");
      TF1 *func = new TF1(Form("fit_%s",hRatio->GetName()), "pol0", 15070000,15170000);
      hRatio->Fit(func,"R0");
      func->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("Equivalent # of MB events on tape per dimuon event (%s)",gTrgSetupTitle[j]),0.035);
      t1->Draw();
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBvsRun_UseAuAuMB.pdf",run_type.Data()));

  // correct for MTD trigger efficiency
  TFile *flumiDep = TFile::Open("Rootfiles/Run14_AuAu200.StudyLumiDep.root","read");
  TF1 *funcppLS = (TF1*)flumiDep->Get("fit_Run15_pp200_hTacDiff_LS_PtBin1");
  TH1F *hTacDiffMean = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffMeanVsRun");
  hTacDiffMean->GetYaxis()->SetRangeUser(785, 795);
  c = draw1D(hTacDiffMean,"Run14_AuAu200: mean of #DeltaTacSum distribution from MT101");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumMeanVsRun.pdf",run_type.Data()));

  TH1F *hTacDiffSigma = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffSigmaVsRun");
  hTacDiffSigma->GetYaxis()->SetRangeUser(4, 12);
  c = draw1D(hTacDiffSigma,"Run14_AuAu200: width of #DeltaTacSum distribution from MT101");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumWidthVsRun.pdf",run_type.Data()));

  TH1F *hTacDiffEff = (TH1F*)flumiDep->Get("AuAu200_hMT101TacDiffEffVsRun");
  hTacDiffEff->GetYaxis()->SetRangeUser(0,1.0);
  c = draw1D(hTacDiffEff,"Run14_AuAu200: efficiency of #DeltaTacSum cut for single tracks");
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/TrigEff_TacSumEffVsRun.pdf",run_type.Data()));

  /// get trigger effiiency from MB events
  TFile *fTrig = TFile::Open("output/Run14_AuAu200.MB.Study.MtdTrig.root","read");
  TH2F *hMT101TacDiffVsRun = (TH2F*)fTrig->Get("mhMT101TacDiffVsRun_di_mu");
  hMT101TacDiffVsRun->Sumw2();
  TH1F *hTacDiffEffAuAu = (TH1F*)hTacDiffEff->Clone("hTacDiffEffAuAu");
  hTacDiffEffAuAu->Reset("AC");
  gStyle->SetOptStat(1);
  TCanvas *c = new TCanvas("cAuAuMbTacSumDiff", "cAuAuMbTacSumDiff", 1100, 700);
  c->Divide(6,5);
  for(int j=1; j<gNTrgSetup; j++)
    {
      const int nbins = hEqvRatioInTrig[j-1]->GetNbinsX();
      for(int bin=1; bin<=nbins; bin++)
	{
	  if(hEqvRatioInTrig[j-1]->GetBinContent(bin)<=0) continue;
	  double run = hEqvRatioInTrig[j-1]->GetBinCenter(bin);
	  int ibin = hTacDiffEff->FindBin(run);
	  double eff1 = hTacDiffEff->GetBinContent(ibin);
	  int jbin = hMT101TacDiffVsRun->GetXaxis()->FindBin(run);
	  TH1F *htmp = (TH1F*)hMT101TacDiffVsRun->ProjectionY(Form("hproj_run%1.0f",run),jbin,jbin);
	  if(run>=15166001 && run<=15166030)
	    {
	      c->cd(run-15166001+1);
	      htmp->Draw();
	    }
	  double min = 786;
	  if(j>=3) 
	    {
	      min = 789;
	      if(run<=15093041 && run>=15089022) min = 790;
	    }
	  if(htmp->Integral()>0) 
	    {
	      double all = htmp->Integral(0,-1);
	      double acc = htmp->Integral(htmp->FindBin(min+1e-4), htmp->FindBin(837-1e-4));
	      double eff2 = acc/all;
	      hTacDiffEffAuAu->SetBinContent(ibin, eff2);
	      hTacDiffEffAuAu->SetBinError(ibin, sqrt(acc)/acc*eff2);
	    }
	}
    }
  draw1D(hTacDiffEffAuAu);

  // acceptance correction
  const int nRunRange = 8;
  int runRange[nRunRange+1] = {15074104, 15077035, 15078021, 15098066, 15099002, 15106130, 15131038, 15132019, 15167014};
  double accEff[nRunRange];  
  TFile *fAcc = TFile::Open(Form("Rootfiles/%s.AcceptanceLoss.root",run_type.Data()),"read");
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

  // MTD response efficiency correction
  TFile *fRespEff = TFile::Open("Rootfiles/Run14_AuAu200.MtdRespEff.root", "read");
  TH1F *hRespEff[2];
  hRespEff[0] = (TH1F*)fRespEff->Get("Cosmic_JpsiEffVsCent_Btm_6300_TrigStudy");
  hRespEff[1] = (TH1F*)fRespEff->Get("Cosmic_JpsiEffVsCent_Btm_6400_TrigStudy");

  // Apply efficiency correction
  double eqMBevt[gNTrgSetup];
  double eqMBevtCorr[gNTrgSetup];
  double eqMBevtGood[gNTrgSetup];
  double eqMBevtGoodCorr[gNTrgSetup];
  eqMBevt[0] = 0;
  eqMBevtCorr[0] = 0;
  eqMBevtGood[0] = 0;
  eqMBevtGoodCorr[0] = 0;
  TH1F *hEqvRatioInTrigCorr[4];
  TH1F *hEqvRatioInTrigCorrAcc[4];
  TH1F *hEqvRatioInTrigCorrAccResp[4];
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorr[j-1] = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_EffCorr",hEqvRatioInTrig[j-1]->GetName()));
      hEqvRatioInTrigCorr[j-1]->Reset("");
      hEqvRatioInTrigCorrAcc[j-1] = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_EffCorrAcc",hEqvRatioInTrig[j-1]->GetName()));
      hEqvRatioInTrigCorrAcc[j-1]->Reset("");
      hEqvRatioInTrigCorrAccResp[j-1] = (TH1F*)hEqvRatioInTrig[j-1]->Clone(Form("%s_EffCorrAccResp",hEqvRatioInTrig[j-1]->GetName()));
      hEqvRatioInTrigCorrAccResp[j-1]->Reset("");

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

	  // trigger efficiency
	  int jbin = hTacDiffEff->GetXaxis()->FindBin(run);
	  //double eff = hTacDiffEff->GetBinContent(jbin);
	  double eff = hTacDiffEffAuAu->GetBinContent(jbin);
	  if(eff<=0) continue;
	  hEqvRatioInTrigCorr[j-1]->SetBinContent(bin, ratio*eff*eff);
	  hEqvRatioInTrigCorr[j-1]->SetBinError(bin, hEqvRatioInTrig[j-1]->GetBinError(bin)*eff*eff);
	  if(run>=15109000 && run<=15109035)
	    {
	      printf("[i] run = %1.0f: ratio = %2.1f, eff = %4.2f, after = %2.1f\n",run,ratio,eff,ratio*eff*eff);
	    }
	  //if(run==15109008 || run==15109030) continue;

	  // acceptance 
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
	  hEqvRatioInTrigCorrAcc[j-1]->SetBinContent(bin, ratio*eff*eff*accCorr);
	  hEqvRatioInTrigCorrAcc[j-1]->SetBinError(bin, hEqvRatioInTrig[j-1]->GetBinError(bin)*eff*eff*accCorr);

	  // relaitve response efficiency w.r.t. 6400V
	  double respEff = 1;
	  if(run<15120072) respEff = 1;
	  else if(run<15125034) respEff = hRespEff[1]->GetBinContent(1)/hRespEff[0]->GetBinContent(1);
	  else if(run<15134021) respEff = 1;
	  else respEff = hRespEff[1]->GetBinContent(1)/hRespEff[0]->GetBinContent(1);
	  hEqvRatioInTrigCorrAccResp[j-1]->SetBinContent(bin, ratio*eff*eff*accCorr*respEff);
	  hEqvRatioInTrigCorrAccResp[j-1]->SetBinError(bin, hEqvRatioInTrig[j-1]->GetBinError(bin)*eff*eff*accCorr*respEff);

	  // good MB events within 0-80% with event weights
	  double ratio_good = hEqvRatioInTrigGood[j-1]->GetBinContent(bin);
	  hEqvRatioInTrigGood[j-1]->SetBinContent(bin,  ratio_good*eff*eff*accCorr*respEff);
	  double nEventsTaken = hNeventsTake->GetBinContent(hNeventsTake->FindBin(run));

	  eqMBevt[j] += ratio * nEventsTaken;
	  eqMBevt[0] += ratio * nEventsTaken;

	  eqMBevtCorr[j] += ratio * nEventsTaken * eff * eff * accCorr *respEff;
	  eqMBevtCorr[0] += ratio * nEventsTaken * eff * eff * accCorr *respEff;

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

  // apply trigger efficiency
  TCanvas *c = new TCanvas("Comp_EqMbPerTrigCorr","Comp_EqMbPerTrigCorr",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorr[j-1]->GetYaxis()->SetRangeUser(0,15);
      if(j==1) 
	{
	  hEqvRatioInTrigCorr[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB per dimuon (MtdTrigEff)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigCorr[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBTrigCorrvsRun.pdf",run_type.Data()));


  // + apply acceptance correction
  TCanvas *c = new TCanvas("Comp_EqMbPerTrigAccCorr","Comp_EqMbPerTrigAccCorr",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorrAcc[j-1]->GetYaxis()->SetRangeUser(0,15);
      if(j==1) 
	{
	  hEqvRatioInTrigCorrAcc[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB per dimuon (MtdTrigEff+Accept)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigCorrAcc[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBTrigAccCorrvsRun.pdf",run_type.Data()));


  // + apply response efficiency
  TCanvas *c = new TCanvas("Comp_EqMbPerTrigAccCorrResp","Comp_EqMbPerTrigAccCorrResp",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigCorrAccResp[j-1]->GetYaxis()->SetRangeUser(0,15);
      if(j==1) 
	{
	  hEqvRatioInTrigCorrAccResp[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB per dimuon (MtdTrigEff+Accept+RespEff)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigCorrAccResp[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_AllMBTrigAccCorrRespvsRun.pdf",run_type.Data()));

  // Good MB events with all corrections
  TCanvas *c = new TCanvas("Comp_EqMbPerTrigGood","Comp_EqMbPerTrigGood",800,600);
  for(int j=1; j<gNTrgSetup; j++)
    {
      hEqvRatioInTrigGood[j-1]->SetMarkerStyle(19+j);
      hEqvRatioInTrigGood[j-1]->SetMarkerColor(gColor[j-1]);
      hEqvRatioInTrigGood[j-1]->SetLineColor(gColor[j-1]);
      hEqvRatioInTrigGood[j-1]->SetTitle(";runId;");
      hEqvRatioInTrigGood[j-1]->GetYaxis()->SetRangeUser(0,15);
      if(j==1) 
	{
	  hEqvRatioInTrigGood[j-1]->Draw("P");
	  TPaveText *t1 = GetTitleText("Equivalent # of MB per dimuon (vtx, 0-80%, MtdTrigEff+Accept+RespEff)",0.035);
	  t1->Draw();
	}
      else     hEqvRatioInTrigGood[j-1]->Draw("samesP");
    }
  if(savePlot) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_lumiDep/EqMB_GoodMBTrigAccCorrvsRun.pdf",run_type.Data()));
  

  //return;

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
