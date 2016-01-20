TFile *f;
const int year = YEAR;

//================================================
void ana_VtxEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TString cut_name = run_config;
  if(year==2013)
    {
      run_type = "Run13_pp500";
      f = TFile::Open(Form("./output/Run13.pp500.jpsi.VPDMB.VtxEff.root"),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
    }

  makeVtxEff();
}

//================================================
void makeVtxEff(int savePlot = 1, int saveHisto = 1)
{
  THnSparseF *hn = (THnSparseF*)f->Get("mhMbEvtEff");

  // ranking cut
  TH1F *hRankingAll = (TH1F*)hn->Projection(0);
  hRankingAll->SetName("hRankingAll");
  hRankingAll->Sumw2();
  printf("[i] All events: %1.0f\n",hRankingAll->GetEntries());
  hn->GetAxis(3)->SetRange(3,3);
  TH1F *hRankingCut = (TH1F*)hn->Projection(0);
  hRankingCut->Sumw2();
  hRankingCut->SetName("hRankingCut");  
  printf("[i] Ranking>0: %1.0f, %4.2f%%\n",hRankingCut->GetEntries(),hRankingCut->GetEntries()/hRankingAll->GetEntries()*100);
  TH1F *hRankingCutRatio = (TH1F*)hRankingCut->Clone("hRankingCutRatio");
  hRankingCutRatio->Divide(hRankingAll);
  hRankingCutRatio->SetMarkerStyle(21);
  c = draw1D(hRankingCutRatio,"Fraction of events with ranking > 0");
  /*
  TF1 *func1 = new TF1("func1","pol0");
  hRankingCutRatio->Fit(func1);
  return;
  */
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.RankingCut.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.RankingCut.png",run_type,run_config));
    }

  // TPC vz cut
  hn->GetAxis(1)->SetRangeUser(-50,50);
  TH1F *hTpcVzCut = (TH1F*)hn->Projection(0);
  hTpcVzCut->SetName("hTpcVzCut");  
  hTpcVzCut->Sumw2();
  printf("[i] |TpcVz| < 50 cm : %1.0f, %4.2f%%\n",hTpcVzCut->GetEntries(),hTpcVzCut->GetEntries()/hRankingCut->GetEntries()*100);
  TH1F *hTpcVzCutRatio = (TH1F*)hTpcVzCut->Clone("hTpcVzCutRatio");
  hTpcVzCutRatio->Divide(hRankingCut);
  hTpcVzCutRatio->SetMarkerStyle(21);
  c = draw1D(hTpcVzCutRatio,"Fraction of events with |TpcVz| < 50 cm");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.TpcVzCut.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.TpcVzCut.png",run_type,run_config));
    }

  // TPC-VPD vz cut
  hn->GetAxis(2)->SetRangeUser(-5,5);
  TH1F *hDVzCut = (TH1F*)hn->Projection(0);
  hDVzCut->SetName("hDVzCut");  
  hDVzCut->Sumw2();
  printf("[i] |TpcVz-VpdVz| < 5 cm : %1.0f, %4.2f%%\n",hDVzCut->GetEntries(),hDVzCut->GetEntries()/hTpcVzCut->GetEntries()*100);
  TH1F *hDVzCutRatio = (TH1F*)hDVzCut->Clone("hDVzCutRatio");
  hDVzCutRatio->Divide(hTpcVzCut);
  hDVzCutRatio->SetMarkerStyle(21);
  c = draw1D(hDVzCutRatio,"Fraction of events with |TpcVz-VpdVz| < 5 cm");
  if(savePlot)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.TpcVpdVzCut.pdf",run_type,run_config));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_VtxEff/VPDMB.TpcVpdVzCut.png",run_type,run_config));
    }
  
  TH1F *hVtxEff = new TH1F("hVtxEff","Vtx cut efficiency",1,0,1);
  hVtxEff->SetBinContent(1,hDVzCut->GetEntries()/hRankingAll->GetEntries());
  printf("[i] total efficiency = %4.2f%%\n",hDVzCut->GetEntries()/hRankingAll->GetEntries()*100);
  if(saveHisto)
    {
      TFile *fout = 0x0;
      if(year==2013) fout = TFile::Open("Rootfiles/Run13.pp500.VPDMB.VtxCutEff.root","recreate");
      hVtxEff->Write();
      fout->Close();
    } 
}
