const int year = 2014;

//================================================
void ana_Lumi()
{
  gStyle->SetOptStat(0);
  if(year==2013)
    {
      run_type = "Run13_pp500";
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }

  MB();
  //MTD();
}

//================================================
void MTD(const int savePlot = 1)
{
  TH1F *hRunId[2];
  TH1F *hZdcRate[2];
  for(int j=0; j<3; j++)
    {
      TFile *f = TFile::Open(Form("output/temp/%d.root",10+j),"read");
      THnSparseF *hnQA = (THnSparseF*)f->Get("mhEvtWithTwoMuons_di_mu");
      hnQA->SetName(Form("%s_%d",hnQA->GetName(),j));
      for(int i=0; i<2; i++)
	{
	  if(i==1) hnQA->GetAxis(5)->SetRange(2,2);
	  TH1F *htmp = (TH1F*)hnQA->Projection(0);
	  htmp->SetName(Form("RunId_%d_%d",j,i));
	  if(j==0) hRunId[i] = (TH1F*)htmp->Clone(Form("RunId_%d",i));
	  else     hRunId[i]->Add(htmp);
	  
	  htmp = (TH1F*)hnQA->Projection(3);
	  htmp->SetName(Form("hZdcRate_%d_%d",j,i));
	  if(j==0) hZdcRate[i] = (TH1F*)htmp->Clone(Form("hZdcRate_%d",i));
	  else     hZdcRate[i]->Add(htmp);
	}
    }
  for(int i=0; i<2; i++)
    {
      hRunId[i]->Sumw2();
      hZdcRate[i]->Sumw2();
    }
  TH1F *hRunIdEff = (TH1F*)hRunId[1]->Clone("RunId_Eff");
  hRunIdEff->Divide(hRunId[0]);
  c = draw1D(hRunIdEff);
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_RunId.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_RunId.png",run_type));
    }

  TH1F *hZdcRateEff = (TH1F*)hZdcRate[1]->Clone("ZdcRate_Eff");
  hZdcRateEff->Divide(hZdcRate[0]);
  hZdcRateEff->GetXaxis()->SetRangeUser(0,100);
  hZdcRateEff->SetMarkerStyle(21);
  c = draw1D(hZdcRateEff,"Fraction of events with two muon candidates;ZdcRate (kHz);Fraction");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_ZdcRate.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/Dimuon_TwoMuonEff_vs_ZdcRate.png",run_type));
    }
}

//================================================
void MB(const int savePlot = 1)
{
  //gStyle->SetOptStat(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *fData = TFile::Open("output/Pico.Run14.AuAu200.jpsi.TightTof.root","read");
  TH1F *hVtxZMtd = (TH1F*)fData->Get("mhTpcVz_di_mu");
  hVtxZMtd->Scale(1./hVtxZMtd->GetBinContent(hVtxZMtd->FindBin(0)));
  hVtxZMtd->SetLineColor(4);
  TH2F *hVzDiffVsTpcVz = (TH2F*)fData->Get("mhVzDiffVsTpcVz_di_mu");
  TH1F *hDzMtd = (TH1F*)hVzDiffVsTpcVz->ProjectionY(Form("hDzMtd"),hVzDiffVsTpcVz->GetXaxis()->FindBin(-100),hVzDiffVsTpcVz->GetXaxis()->FindBin(100));


  TFile *fMB = TFile::Open("output/Pico.Run14.AuAu200.MB.root","read");
  THnSparseF *hnCentQA = (THnSparseF*)fMB->Get("mhEventCent_qa_mb");
  TH1F *hVertexMB = (TH1F*)hnCentQA->Projection(4);
  hVertexMB->Scale(1./hVertexMB->GetBinContent(hVertexMB->FindBin(0)));
  hVertexMB->SetLineColor(2);
  //hVertexMB->Draw("sames HIST");
  //printf("MB 2: %2.2f%%\n",hVertexMB->Integral(hVertexMB->FindBin(-100),hVertexMB->FindBin(100))/hVertexMB->Integral(0,-1)*100);

  const char *name[2] = {"Mid","Low"};
  TFile *f[2];
  TH1F *hVtxZDimuonNoCut[2];
  TH1F *hDzDimuonVzCut[2];
  TH1F *hVtxZNoCut[2];
  TH1F *hVtxDzVzCut[2];
  TH1F *hZdcRateAllWithCut[2];
  TH1F *hZdcRateDimuonWithCut[2];
  TH1F *hRFvsRate[2];

  for(int i=0; i<1; i++)
    {
      printf("+++ process production_%s\n",name[i]);
      f[i] = TFile::Open(Form("output/Run14.AuAu200.MB.Prod%s.root",name[i]),"read");
      THnSparseF *hnLumiQA = (THnSparseF*)f[i]->Get("mhLumiQA_mb");
      hnLumiQA->SetName(Form("%s_%s",hnLumiQA->GetName(),name[i]));
      if(i==0) hnLumiQA->GetAxis(0)->SetRangeUser(15119017,15163026);
      //if(i==0) hnLumiQA->GetAxis(0)->SetRangeUser(15163026,15167014);

      hnLumiQA->GetAxis(10)->SetRange(2,2);
      hVtxZDimuonNoCut[i] = (TH1F*)hnLumiQA->Projection(6);
      hVtxZDimuonNoCut[i]->SetName(Form("hVtxZDimuonNoCut_%s",name[i]));
      cout << "Dimuon before cuts: " << hVtxZDimuonNoCut[i]->Integral(0,-1) << endl;
      hnLumiQA->GetAxis(6)->SetRangeUser(-100,100);
      hDzDimuonVzCut[i] = (TH1F*)hnLumiQA->Projection(7);
      hDzDimuonVzCut[i]->SetName(Form("hDzDimuonVzCut_%s",name[i]));
      hnLumiQA->GetAxis(6)->SetRange(0,-1);
      hnLumiQA->GetAxis(10)->SetRange(0,-1);

      int eventAll = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventAll << endl;

      hnLumiQA->GetAxis(4)->SetRange(2,2);
      int eventVtx = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventVtx << endl;
      
      hnLumiQA->GetAxis(5)->SetRange(2,2);
      int eventGoodVtx = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventGoodVtx << endl;

      /*
	hnLumiQA->GetAxis(8)->SetRange(2,2);
	int eventVptTacHigh = hnLumiQA->Projection(3)->Integral(0,-1);
	cout << eventVptTacHigh << endl;

	hnLumiQA->GetAxis(9)->SetRange(2,2);
	int eventVptTacLow = hnLumiQA->Projection(3)->Integral(0,-1);
	cout << eventVptTacLow << endl;
      */

      hVtxZNoCut[i] = (TH1F*)hnLumiQA->Projection(6);
      cout << hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100)) << endl;
      hVtxZNoCut[i]->SetName(Form("hVtxZNoCut_%s",name[i]));
      hnLumiQA->GetAxis(6)->SetRangeUser(-100,100);
      int eventTpcVz = hnLumiQA->Projection(3)->Integral(0,-1);
      cout << eventTpcVz << endl;

      hVtxDzVzCut[i] = (TH1F*)hnLumiQA->Projection(7);
      hVtxDzVzCut[i]->SetName(Form("hVtxDzVzCut_%s",name[i]));
      hnLumiQA->GetAxis(7)->SetRangeUser(-3,3);

      hZdcRateAllWithCut[i] = (TH1F*)hnLumiQA->Projection(3);
      hZdcRateAllWithCut[i]->SetName(Form("hZdcRateAllWithCut_%s",name[i]));
      int eventVzDiff = hZdcRateAllWithCut[i]->Integral(0,-1);
      cout << eventVzDiff << endl;

      hnLumiQA->GetAxis(10)->SetRange(2,2);
      hZdcRateDimuonWithCut[i] = (TH1F*)hnLumiQA->Projection(3);
      hZdcRateDimuonWithCut[i]->SetName(Form("hZdcRateDimuonWithCut_%s",name[i]));
      int eventDimuon = hZdcRateDimuonWithCut[i]->Integral(0,-1);
      cout << eventDimuon << endl;

      // vz distribution
      hVtxZNoCut[i]->Scale(1./hVtxZNoCut[i]->GetBinContent(hVtxZNoCut[i]->FindBin(0)));
      hVtxZNoCut[i]->GetYaxis()->SetRangeUser(0.01,10);
      c = draw1D(hVtxZNoCut[i],Form("TPC vz distribution (production_%s)",name[i]),kTRUE,kFALSE);
      printf("MB 1: %2.2f%%\n",hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100))/hVtxZNoCut[i]->Integral(0,-1)*100);

      hVtxZDimuonNoCut[i]->SetMarkerStyle(20);
      hVtxZDimuonNoCut[i]->Sumw2();
      hVtxZDimuonNoCut[i]->Scale(1./hVtxZDimuonNoCut[i]->GetBinContent(hVtxZDimuonNoCut[i]->FindBin(0)));
      hVtxZDimuonNoCut[i]->SetLineColor(6);
      hVtxZDimuonNoCut[i]->SetMarkerColor(6);
      hVtxZDimuonNoCut[i]->Draw("sames P");
  
      hVtxZMtd->Draw("sames HIST");
      printf("Dimuon 1: %2.2f%%\n",hVtxZMtd->Integral(hVtxZMtd->FindBin(-100),hVtxZMtd->FindBin(100))/hVtxZMtd->Integral(0,-1)*100);

      TLegend *leg = new TLegend(0.15,0.68,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hVtxZNoCut[i],Form("VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxZNoCut[i]->Integral(hVtxZNoCut[i]->FindBin(-100),hVtxZNoCut[i]->FindBin(100))/hVtxZNoCut[i]->Integral(0,-1)*100),"L");
      leg->AddEntry(hVtxZDimuonNoCut[i],Form("dimuon in VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxZDimuonNoCut[i]->Integral(hVtxZDimuonNoCut[i]->FindBin(-100),hVtxZDimuonNoCut[i]->FindBin(100))/hVtxZDimuonNoCut[i]->Integral(0,-1)*100),"P");
      leg->AddEntry(hVtxZMtd,Form("Production_high dimuon trigger: %2.1f%%",hVtxZMtd->Integral(hVtxZMtd->FindBin(-100),hVtxZMtd->FindBin(100))/hVtxZMtd->Integral(0,-1)*100),"L");
      leg->Draw();

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVtxZ.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVtxZ.png",run_type,name[i]));
	}

      // dz distribution
      hVtxDzVzCut[i]->Scale(1./hVtxDzVzCut[i]->GetBinContent(hVtxDzVzCut[i]->FindBin(0)));
      hVtxDzVzCut[i]->GetYaxis()->SetRangeUser(1e-5,1e2);
      c = draw1D(hVtxDzVzCut[i],Form("vz_{TPC}-vz_{VPD} distribution (production_%s)",name[i]),kTRUE,kFALSE);
      hDzMtd->Scale(1./hDzMtd->GetBinContent(hDzMtd->FindBin(0)));
      hDzMtd->SetLineColor(4);
      hDzMtd->Draw("sames HIST");
      hDzDimuonVzCut[i]->Sumw2();
      hDzDimuonVzCut[i]->Scale(1./hDzDimuonVzCut[i]->GetBinContent(hDzDimuonVzCut[i]->FindBin(0)));
      hDzDimuonVzCut[i]->SetLineColor(6);
      hDzDimuonVzCut[i]->SetMarkerStyle(20);
      hDzDimuonVzCut[i]->SetMarkerColor(6);
      hDzDimuonVzCut[i]->Draw("sames P");
      TLegend *leg = new TLegend(0.15,0.68,0.3,0.88);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hVtxDzVzCut[i],Form("VPD-ZDC-NoVtx-MB: %2.1f%%",hVtxDzVzCut[i]->Integral(hVtxDzVzCut[i]->FindBin(-3),hVtxDzVzCut[i]->FindBin(3))/hVtxDzVzCut[i]->Integral(0,-1)*100),"L");
      leg->AddEntry(hDzDimuonVzCut[i],Form("dimuon in VPD-ZDC-NoVtx-MB: %2.1f%%",hDzDimuonVzCut[i]->Integral(hDzDimuonVzCut[i]->FindBin(-3),hDzDimuonVzCut[i]->FindBin(3))/hDzDimuonVzCut[i]->Integral(0,-1)*100),"P");
      leg->AddEntry(hDzMtd,Form("Production_high dimuon trigger: %2.1f%%",hDzMtd->Integral(hDzMtd->FindBin(-3),hDzMtd->FindBin(3))/hDzMtd->Integral(0,-1)*100),"L");
      leg->Draw();
      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVpdDz.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_TpcVpdDz.png",run_type,name[i]));
	}


      // vs zdc rate
      hZdcRateAllWithCut[i]->Sumw2();
      hZdcRateDimuonWithCut[i]->Sumw2();
      hRFvsRate[i] = (TH1F*)hZdcRateAllWithCut[i]->Clone(Form("hRF_vs_zdc_%s",name[i]));
      hRFvsRate[i]->Divide(hZdcRateDimuonWithCut[i]);
      hRFvsRate[i]->SetMarkerStyle(21+i*4);
      hRFvsRate[i]->SetMarkerColor(i+1);
      hRFvsRate[i]->SetLineColor(i+1);
      hRFvsRate[i]->GetXaxis()->SetRangeUser(0,65);
      hRFvsRate[i]->GetYaxis()->SetRangeUser(0,50);
      TH1F *htmp = (TH1F*)hRFvsRate[i]->Clone(Form("hRF_vs_zdc_%s_clone",name[i]));
      htmp->GetXaxis()->SetRangeUser(25,65);
      htmp->GetYaxis()->SetRangeUser(5,45);
      htmp->SetMarkerStyle(21);
      htmp->SetMarkerColor(1);
      htmp->SetLineColor(1);
      TF1 *func = new TF1(Form("func_%s",name[i]),"pol0",25,65);
      htmp->Fit(func,"IR0Q");
      c = draw1D(htmp,Form("Rejection factor (production_%s);zdcRate (kHz)",name[i]));
      func->SetLineColor(2);
      func->Draw("sames");

      if(savePlot)
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_SF_vs_Lumi.pdf",run_type,name[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/%s_SF_vs_Lumi.png",run_type,name[i]));
	}
    }
  return;
  c = draw1D(hRFvsRate[1],Form("Rejection factor;zdcRate (kHz)"));
  hRFvsRate[0]->Draw("sames");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/SF_vs_Lumi_MidLow.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Lumi/SF_vs_Lumi_MidLow.png",run_type));
    }

}
