const int year = YEAR;

TFile *f;

//================================================
void ana_MtdTrigEff()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.98);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  trigElecEff();
}

//================================================
void trigElecEff(const int savePlot = 1)
{
  if(year==2014)
    {
      f = TFile::Open("output/Run14.AuAu200.MB.TrigElecEff.root","read");
    }
  else
    {
      printf("[e] No available input file!\n");
      return;
    }

  THnSparseF *hnTrigEff = (THnSparseF*)f->Get("mhMtdTrigElecEff_mb");
  const int nbins = 11;
  const double xbins[nbins+1] = {0,1,1.5,2,2.5,3,3.5,4,5,6,8,10};
  
  TList *list = new TList;
  // Efficiency vs. dTof
  const int nDtof = 4;
  const double dtof_cut[nDtof] = {1, 0.75, 0.5, 0.25};
  TH1F *hMuonPtDtof[nDtof][3];
  TH1F *hMuonEffDtof[nDtof][3];
  for(int i=0; i<nDtof; i++)
    {
      hnTrigEff->GetAxis(3)->SetRangeUser(-3,dtof_cut[i]);
      hMuonPtDtof[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][0]->SetName(Form("hMuonPtMatch_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtDtof[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][1]->SetName(Form("hMuonPtQT_DtofCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtDtof[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtDtof[i][2]->SetName(Form("hMuonPtHigh2_DtofCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(3)->SetRange(0,-1);
    }

  TString legName1[nDtof];
  for(int i=0; i<nDtof; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtDtof[i][j]->GetName()),xbins);
	  hMuonEffDtof[i][j] = (TH1F*)hMuonPtDtof[i][j]->Clone(Form("hMuonEffDtof_%d_%d",i,j));
	  hMuonEffDtof[i][j]->Sumw2();
	  hMuonEffDtof[i][j]->Divide(hMuonPtDtof[i][0]);
	  hMuonEffDtof[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffDtof[i][j]); 
	}
      legName1[i] = Form("#Deltatof < %2.2f ns",dtof_cut[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Dtof",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName1,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareDtof.pdf",run_type));

 
  // Efficiency vs. centrality
  // use dtof < 1 ns cut
  TH1F *hMuonPtCent[nCentBins][3];
  TH1F *hMuonEffCent[nCentBins][3];
  for(int i=0; i<nCentBins; i++)
    {
      hnTrigEff->GetAxis(5)->SetRange(centBins_low[i], centBins_high[i]);
      hMuonPtCent[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][0]->SetName(Form("hMuonPtMatch_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtCent[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][1]->SetName(Form("hMuonPtQT_%s",cent_Title[i]));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtCent[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtCent[i][2]->SetName(Form("hMuonPtHigh2_%s",cent_Title[i]));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(5)->SetRange(0,-1);
    }

  TString legName2[nCentBins];
  for(int i=0; i<nCentBins; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hMuonPtCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtCent[i][j]->GetName()),xbins);
	  hMuonEffCent[i][j] = (TH1F*)hMuonPtCent[i][j]->Clone(Form("hMuonEff_%s_%d",cent_Title[i],j));
	  hMuonEffCent[i][j]->Sumw2();
	  hMuonEffCent[i][j]->Divide(hMuonPtCent[i][0]);
	  hMuonEffCent[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffCent[i][j]); 
	}
      legName2[i] = Form("%s%%",cent_Name[i]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_Cent",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName2,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareCent.pdf",run_type));

  // Efficiency vs. TPC vz
  // use dtof < 1 ns cut and 0-60%
  hnTrigEff->GetAxis(5)->SetRange(5,16);
  const int nTpcVz = 6;
  const double tpcvz_cut[7] = {-100,-30,-5,0,5,30,100};
  TH1F *hMuonPtTpcVz[nTpcVz][3];
  TH1F *hMuonEffTpcVz[nTpcVz][3];
  for(int i=0; i<nTpcVz; i++)
    {
      hnTrigEff->GetAxis(4)->SetRangeUser(tpcvz_cut[i]+1e-4,tpcvz_cut[i+1]-1e-4);
      hMuonPtTpcVz[i][0] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][0]->SetName(Form("hMuonPtMatch_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(2,2);
      hMuonPtTpcVz[i][1] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][1]->SetName(Form("hMuonPtQT_TpcVzCut%d",i));

      hnTrigEff->GetAxis(2)->SetRange(2,2);
      hMuonPtTpcVz[i][2] = (TH1F*)hnTrigEff->Projection(0);
      hMuonPtTpcVz[i][2]->SetName(Form("hMuonPtHigh2_TpcVzCut%d",i));

      hnTrigEff->GetAxis(1)->SetRange(0,-1);
      hnTrigEff->GetAxis(2)->SetRange(0,-1);
      hnTrigEff->GetAxis(4)->SetRange(0,-1);

      for(int j=0; j<3; j++)
	{
	  hMuonPtTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Rebin(nbins,Form("%s_rebin",hMuonPtTpcVz[i][j]->GetName()),xbins);
	}
    }

  TString legName3[nTpcVz];
  for(int i=3; i<6; i++)
    {
      TH1F *hRef = (TH1F*)hMuonPtTpcVz[i][0]->Clone(Form("%s_clone",hMuonPtTpcVz[i][0]->GetName()));
      hRef->Add(hMuonPtTpcVz[5-i][0]);
      for(int j=0; j<3; j++)
	{
	  hMuonEffTpcVz[i][j] = (TH1F*)hMuonPtTpcVz[i][j]->Clone(Form("hMuonEff_%d_%d",i,j));
	  hMuonEffTpcVz[i][j]->Sumw2();
	  hMuonEffTpcVz[i][j]->Add(hMuonPtTpcVz[5-i][j]);
	  hMuonEffTpcVz[i][j]->Divide(hRef);
	  hMuonEffTpcVz[i][j]->SetMarkerSize(1.5);
	  if(j==2) list->Add(hMuonEffTpcVz[i][j]); 
	}
      legName3[i-3] = Form("%1.0f < |vz| < %1.0f cm",tpcvz_cut[i],tpcvz_cut[i+1]);
    }
  TCanvas *c = drawHistos(list,"MuonPtEff_TpcVz",Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type),true,0,10,true,0.8,1.05,kFALSE,true,legName3,true,"",0.4,0.65,0.2,0.45,kTRUE);
  list->Clear();
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_CompareTpcVz.pdf",run_type));
  hnTrigEff->GetAxis(5)->SetRange(0,-1);

  // Final efficiency: dtof < 1 ns, 0-60%, |tpcVz| < 100 cm
  TH1F *hTrigElecEff = (TH1F*)hMuonEffCent[0][2]->Clone(Form("%s_TrigElecEff",run_type));
  TF1 *func = new TF1(Form("%s_FitFunc",hTrigElecEff->GetName()),"[0]-exp(-1*[1]*(x-[2]))",1,10);
  func->SetParLimits(0,0,1);
  hTrigElecEff->Fit(func,"R0Q");
  hTrigElecEff->GetYaxis()->SetRangeUser(0.85,1.1);
  c = draw1D(hTrigElecEff,Form("%s: trigger electronics efficiency;p_{T} (GeV/c)",run_type));
  func->SetLineColor(2);
  func->SetLineStyle(2);
  func->Draw("sames");
  if(savePlot)  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_MtdTrigEff/MtdTrigElecEff_Fit.pdf",run_type));
}
