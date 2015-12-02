TFile *f;
const char *run_config = "";
const int year = 2013;
TString run_cfg_name;

//================================================
void qa_FindVtx()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      fileName = Form("Run13.pp500.jpsi.FindVtx.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      fileName = Form("Run14.AuAu200.jpsi.FindVtx.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");
  run_cfg_name = run_config;

  compare();
}

//================================================
void compare(const int saveHisto = 0)
{
  const char *FVtxName[3] = {"default","closest","Mtd"};

  TH1F *hNPrimVtx  = (TH1F*)f->Get(Form("mhNPrimVtx_%s",trigName[kTrigType]));
  hNPrimVtx->GetXaxis()->SetRangeUser(0,20);
  c = draw1D(hNPrimVtx,"",kTRUE,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NPrimVtx.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NPrimVtx.png",run_type));
    }
  printf("[i] events = %d\n",hNPrimVtx->GetEntries());

  TH2F *hNVtxVsNGoodRank  = (TH2F*)f->Get(Form("mhNVtxVsNGoodRank_%s",trigName[kTrigType]));
  hNVtxVsNGoodRank->GetXaxis()->SetRangeUser(0,15);
  hNVtxVsNGoodRank->GetYaxis()->SetRangeUser(0,20);
  c = draw2D(hNVtxVsNGoodRank);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NPrimVtx_vs_NGoodRank.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NPrimVtx_vs_NGoodRank.png",run_type));
    }

  TH2F *hVtxRankVsNMtdTrk  = (TH2F*)f->Get(Form("mhVtxRankVsNMtdTrk_%s",trigName[kTrigType]));
  TH1F *hNMtdTrkPosRank = (TH1F*)hVtxRankVsNMtdTrk->ProjectionY("hNMtdTrkPosRank",3,3);
  TH1F *hNMtdTrkNegRank = (TH1F*)hVtxRankVsNMtdTrk->ProjectionY("hNMtdTrkNegRank",1,1);
  TList *list = new TList;
  list->Add(hNMtdTrkNegRank);
  list->Add(hNMtdTrkPosRank);
  TString legName[2] = {"Ranking < 0","Ranking > 0"};
  c = drawHistos(list,"VtxRankVsNMtdTrk",Form("%s: number of associated primary tracks matched to MTD; N",trigName[kTrigType]),kFALSE,0,100,kFALSE,-3,8,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.65,0.6,0.8,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NMtdTrk_vs_Ranking.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/NMtdTrk_vs_Ranking.png",run_type));
    }

  TH2F *hVtxClosestVsMtdMth  = (TH2F*)f->Get(Form("mhVtxClosestVsMtdMth_%s",trigName[kTrigType]));
  hVtxClosestVsMtdMth->GetXaxis()->SetRangeUser(0,15);
  hVtxClosestVsMtdMth->GetYaxis()->SetRangeUser(0,15);
  c = draw2D(hVtxClosestVsMtdMth);
  double counter = 0;
  for(int bin=1; bin<=hVtxClosestVsMtdMth->GetNbinsX(); bin++)
    {
      counter += hVtxClosestVsMtdMth->GetBinContent(bin,bin);
    }
  printf("[i] Fraction of VPD = MTD vertex is %4.2f%%\n",counter/hVtxClosestVsMtdMth->GetEntries()*100);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VpdVtxIndex_vs_MtdVtxIndex.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VpdVtxIndex_vs_MtdVtxIndex.png",run_type));
    }
  
  TH1F *hFVtxIndex[3];
  TString legName2[3] = {"Default","Closest to VPD","Matched to MTD"};
  list->Clear();
  for(int i=0; i<3; i++)
    {
      hFVtxIndex[i] = (TH1F*)f->Get(Form("mhFVtxIndex_%s_%s",FVtxName[i],trigName[kTrigType]));
      list->Add(hFVtxIndex[i]);
      printf("[i] %s: index = 0 --> %4.2f%%\n",FVtxName[i], hFVtxIndex[i]->GetBinContent(1)/hFVtxIndex[i]->GetEntries()*100);
    }
  c = drawHistos(list,"FVtxIndex",Form("%s: index of selected vertex; Index",trigName[kTrigType]),kTRUE,0,15,kTRUE,1,1e8,kTRUE,kTRUE,legName2,kTRUE,"",0.6,0.75,0.6,0.8,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxIndex.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxIndex.png",run_type));
    }

  TH1F *hFVtxRanking[3];
  list->Clear();
  for(int i=0; i<3; i++)
    {
      hFVtxRanking[i] = (TH1F*)f->Get(Form("mhFVtxRanking_%s_%s",FVtxName[i],trigName[kTrigType]));
      list->Add(hFVtxRanking[i]);
      printf("[i] %s: positive ranking %4.2f%%\n",FVtxName[i], hFVtxRanking[i]->GetBinContent(3)/hFVtxRanking[i]->GetEntries()*100);
    }
  c = drawHistos(list,"FVtxRanking",Form("%s: ranking of selected vertex",trigName[kTrigType]),kTRUE,-1.5,1.5,kTRUE,1,1e8,kFALSE,kTRUE,legName2,kTRUE,"",0.6,0.75,0.6,0.8,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxRanking.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxRanking.png",run_type));
    }

  TH1F *hFVtxNTrkUsed[3];
  list->Clear();
  for(int i=0; i<3; i++)
    {
      hFVtxNTrkUsed[i] = (TH1F*)f->Get(Form("mhFVtxNTrkUsed_%s_%s",FVtxName[i],trigName[kTrigType]));
      list->Add(hFVtxNTrkUsed[i]);
    }
  c = drawHistos(list,"FVtxNTrkUsed",Form("%s: # of tracks used for selected vertex",trigName[kTrigType]),kTRUE,0,100,kTRUE,1e2,5e7,kTRUE,kTRUE,legName2,kTRUE,"",0.6,0.75,0.6,0.8,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxNTrkUsed.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxNTrkUsed.png",run_type));
    }

  TH1F *hFVtxNMtdTrk[3];
  list->Clear();
  for(int i=0; i<3; i++)
    {
      hFVtxNMtdTrk[i] = (TH1F*)f->Get(Form("mhFVtxNMtdTrk_%s_%s",FVtxName[i],trigName[kTrigType]));
      list->Add(hFVtxNMtdTrk[i]);
    }
  c = drawHistos(list,"FVtxNMtdTrk",Form("%s: # of associated tracks matched to MTD for selected vertex",trigName[kTrigType]),kTRUE,0,10,kTRUE,1,1e8,kTRUE,kTRUE,legName2,kTRUE,"",0.6,0.75,0.6,0.8,kFALSE);
  if(saveHisto) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxNMtdTrk.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VtxNMtdTrk.png",run_type));
    }

  TH2F *hFVtxVzDiffVsVz[3];
  for(int i=0; i<3; i++)
    {
      hFVtxVzDiffVsVz[i] = (TH2F*)f->Get(Form("mhFVtxVzDiffVsVz_%s_%s",FVtxName[i],trigName[kTrigType]));
      c = draw2D(hFVtxVzDiffVsVz[i]);
      if(saveHisto) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VzDiffVsVtxVz_%s.pdf",run_type,FVtxName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/VzDiffVsVtxVz_%s.png",run_type,FVtxName[i]));
	}

      TH1F *hDz = (TH1F*)hFVtxVzDiffVsVz[i]->ProjectionY();
      printf("[i] %s: |dz| < 5 cm %4.2f%%\n",FVtxName[i], hDz->Integral(hDz->FindFixBin(-5),hDz->FindFixBin(5))/hDz->GetEntries()*100);
    }

  THnSparseF *hFVtxJpsi[3];
  TH1F *hJpsiM[3][4];
  const char *anaCuts[4] = {"all","ranking > 0","|dz| < 5 cm","ranking > 0 && |dz| < 5 cm"};
  for(int i=0; i<3; i++)
    {
      hFVtxJpsi[i] = (THnSparseF*)f->Get(Form("mhFVtxJpsi_%s_%s",FVtxName[i],trigName[kTrigType]));
      hFVtxJpsi[i]->GetAxis(1)->SetRangeUser(1.5,100);
      
      hJpsiM[i][0] = (TH1F*)hFVtxJpsi[i]->Projection(0);
      hJpsiM[i][0]->SetName(Form("FVtxJpsi_%s_all",FVtxName[i]));

      hFVtxJpsi[i]->GetAxis(2)->SetRange(3,3);
      hJpsiM[i][1] = (TH1F*)hFVtxJpsi[i]->Projection(0);
      hJpsiM[i][1]->SetName(Form("FVtxJpsi_%s_ranking",FVtxName[i]));
      hFVtxJpsi[i]->GetAxis(2)->SetRange(0,-1);

      hFVtxJpsi[i]->GetAxis(3)->SetRangeUser(-5,5);
      hJpsiM[i][2] = (TH1F*)hFVtxJpsi[i]->Projection(0);
      hJpsiM[i][2]->SetName(Form("FVtxJpsi_%s_vz",FVtxName[i]));

      hFVtxJpsi[i]->GetAxis(2)->SetRange(3,3);
      hJpsiM[i][3] = (TH1F*)hFVtxJpsi[i]->Projection(0);
      hJpsiM[i][3]->SetName(Form("FVtxJpsi_%s_vz_ranking",FVtxName[i]));
      hFVtxJpsi[i]->GetAxis(2)->SetRange(0,-1);
      hFVtxJpsi[i]->GetAxis(3)->SetRange(0,-1);

      TCanvas *c = new TCanvas(Form("JpsiM_%s",FVtxName[i]),Form("JpsiM_%s",FVtxName[i]),1100,650);
      c->Divide(2,2);
      for(int j=0; j<4; j++)
	{
	  c->cd(j+1);
	  hJpsiM[i][j]->Sumw2();
	  hJpsiM[i][j]->GetXaxis()->SetRangeUser(2.5,4);
	  hJpsiM[i][j]->SetMarkerStyle(21);
	  hJpsiM[i][j]->SetMaximum(hJpsiM[i][j]->GetMaximum()*1.5);
	  hJpsiM[i][j]->SetMinimum(0);
	  TF1 *func = new TF1(Form("func_%s",hJpsiM[i][j]->GetName()),"gaus(0)+pol3(3)",2.5,3.5);
	  func->SetParameter(1,3.09);
	  func->SetParameter(2,0.1);
	  hJpsiM[i][j]->Fit(func,"0IRQ");
	  hJpsiM[i][j]->SetTitle("");
	  hJpsiM[i][j]->Draw("P");
	  TF1 *funcbkg = new TF1(Form("bkg_%s",hJpsiM[i][j]->GetName()),"pol3",2.5,3.5);
	  funcbkg->SetLineColor(4);
	  for(int k=0; k<4; k++)
	    {
	      funcbkg->SetParameter(k,func->GetParameter(3+k));
	      funcbkg->SetParError(k,func->GetParError(3+k));
	    }
	  funcbkg->Draw("sames");
	  func->SetLineColor(2);
	  func->Draw("sames");
	  TPaveText *t1 = GetTitleText(Form("%s vertex: %s",legName2[i].Data(),anaCuts[j]),0.065);
	  t1->Draw();

	  double low_mass = 2.9;
	  double high_mass = 3.2;
	  int low_bin = hJpsiM[i][j]->FindFixBin(low_mass+1e-4);
	  int high_bin = hJpsiM[i][j]->FindFixBin(high_mass-1e-4);
	  double error;
	  double count = hJpsiM[i][j]->IntegralAndError(low_bin,high_bin,error);
	  double bkg = funcbkg->Integral(low_mass,high_mass) * 1./hJpsiM[i][j]->GetBinWidth(1);
	  double bkg_err = sqrt(bkg);
	  double signal = count - bkg;
	  double sig_err = TMath::Sqrt(error*error+bkg_err*bkg_err);
	  TPaveText *t = GetPaveText(0.16,0.3,0.7,0.88,0.05);
	  t->SetTextFont(62);
	  t->AddText("Counting");
	  t->AddText(Form("%1.0f #pm %1.0f (%1.1f#sigma)",signal,sig_err,signal/sig_err));
	  t->SetTextAlign(11);
	  t->Draw();
	  printf("[i] %s: %1.0f +/- %1.0f (%1.1f)\n",FVtxName[i], signal,sig_err,signal/sig_err);

	  double height = hJpsiM[i][j]->GetMaximum() * 0.6;
	  TLine *line = GetLine(low_mass, 0, low_mass, height, 1);
	  line->Draw();
	  TLine *line = GetLine(high_mass, 0, high_mass, height, 1);
	  line->Draw();
	}
      if(saveHisto) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/JpsiCounts_%s.pdf",run_type,FVtxName[i]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_FindVtx/JpsiCounts_%s.png",run_type,FVtxName[i]));
	}
    }
}
