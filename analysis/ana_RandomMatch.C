TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const Double_t pt_cuts[5] = {1,2,4,10,20};
const TString legName_charge[2] = {"negative","positive"};
//================================================
void ana_RandomMatch()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  //f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  f = TFile::Open("/Users/admin/Work/STAR/analysis/output/Run13.pp500.jpsi.CutRanking.root","read");

  MonteCarlo();
  //makeHisto();
}

//================================================
void MonteCarlo(const Int_t save = 1)
{
  const char *title[2] = {"z","y"};
  const Int_t type = 1;

  TFile *fin = TFile::Open("Rootfiles/Run13.pp500.RandomMatch.root","read");

  // MTD hit
  TH2F *hHitYVsModule = (TH2F*)fin->Get("mhMthMtdHitLocaly_di_mu");
  c = draw2D(hHitYVsModule);
  TH1F *hHitY = (TH1F*)hHitYVsModule->ProjectionY("hMtdHitLocaly");
  hHitY->Scale(1./hHitY->Integral());
  c = draw1D(hHitY,"Local y distribution of MTD hits",kFALSE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_MtdHit_localy.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_MtdHit_localy.png",run_type));
    }

  TList *list = new TList;
  // track y
  TH1F *hTrkYinPt[4][2];
  for(Int_t i=0; i<4; i++)
    {
      list->Clear();
      for(Int_t j=0; j<2; j++)
	{
	  hTrkYinPt[i][j] = (TH1F*)fin->Get(Form("Data_TrkYAtMtd_pt%1.0f_%1.0f_%s",pt_cuts[i],pt_cuts[i+1],legName_charge[j].Data()));
	  list->Add(hTrkYinPt[i][j]);
	}
      c = drawHistos(list,Form("Trk%s_pt%1.0f_%1.0f",title[type],pt_cuts[i],pt_cuts[i+1]),Form("%s distribution of tracks projected to MTD radius;%s (cm)",title[type],title[type]),kFALSE,0,100,kFALSE,1e-1,6e6,kFALSE,kTRUE,legName_charge,kTRUE,Form("%1.0f < p_{T} < %1.0f GeV/c",pt_cuts[i],pt_cuts[i+1]),0.15,0.25,0.65,0.85,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_Trk_localy_pt_%1.0f_%1.0f.pdf",run_type,pt_cuts[i],pt_cuts[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_Trk_localy_pt_%1.0f_%1.0f.png",run_type,pt_cuts[i],pt_cuts[i+1]));
	}
    }

  TH1F *hDyMC[4][2];

  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hDyMC[i][j]  = new TH1F(Form("MC_hD%s_pt%1.0f_%1.0f_%s",title[type],pt_cuts[i],pt_cuts[i+1],legName_charge[j].Data()),Form("Monte Carlo: #Delta%s distribution (%s, %1.0f < p_{T} < %1.0f); #Delta%s (cm)",title[type],legName_charge[j].Data(),pt_cuts[i],pt_cuts[i+1],title[type]),200,-100,100);
	}
    }
  //Double_t center[20] = {-24.2,-24.2,-24.2,-24.2,-19.8,-19.8,-15.4,-11,-6.6,-2.2,2.2,6.6,11,15.4,19.8,19.8,24.2,24.2,24.2,24.2};
  const Int_t nExpr = 1e5;

  gRandom->SetSeed(0);
  Double_t hit_y, trk_y;
  for(Int_t k=0; k<nExpr; k++)
    {
      for(Int_t i=0; i<4; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      hit_y = hHitY->GetRandom();
	      trk_y = hTrkYinPt[i][j]->GetRandom();
	      hDyMC[i][j]->Fill(trk_y-hit_y);
	    }
	}
    }

  
  for(Int_t i=0; i<4; i++)
    {
      list->Clear();
      for(Int_t j=0; j<2; j++)
	{
	  hDyMC[i][j]->Scale(1./hDyMC[i][j]->Integral());
	  hDyMC[i][j]->SetMaximum(1.2*hDyMC[i][j]->GetMaximum());
	  list->Add(hDyMC[i][j]);
	}

      c = drawHistos(list,Form("MC_hD%s_pt%1.0f_%1.0f",title[type],pt_cuts[i],pt_cuts[i+1]),Form("Monte Carlo: #Delta%s distribution; #Delta%s (cm)",title[type],title[type]),kFALSE,0,100,kFALSE,1e-1,6e6,kFALSE,kTRUE,legName_charge,kTRUE,Form("%1.0f < p_{T} < %1.0f GeV/c",pt_cuts[i],pt_cuts[i+1]),0.15,0.25,0.65,0.85,kFALSE);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_MC_dy_pt_%1.0f_%1.0f.pdf",run_type,pt_cuts[i],pt_cuts[i+1]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_MC_dy_pt_%1.0f_%1.0f.png",run_type,pt_cuts[i],pt_cuts[i+1]));
	}
    }

  // Compare to data
  c = new TCanvas("hDy_TrkPtBin","hDy_TrkPtBin",1200,650);
  c->Divide(2,2);

  TH1F *hDyData[4][2];
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hDyData[i][j] = (TH1F*)fin->Get(Form("Data_TrkDy_pt%1.0f_%1.0f_%s",pt_cuts[i],pt_cuts[i+1],legName_charge[j].Data()));
	  hDyData[i][j]->SetLineColor(j+1);
	  hDyData[i][j]->Scale(1./hDyData[i][j]->Integral());
	  hDyData[i][j]->SetMaximum(1.3*hDyData[i][j]->GetMaximum());
	  hDyData[i][j]->SetTitle("");
	  hDyMC[i][j]->SetLineColor(color[j+2]);
	  //hDyMC[i][j]->SetLineStyle(2);
	}
      c->cd(i+1);
      hDyData[i][0]->Draw();
      hDyData[i][1]->Draw("sames");
      hDyMC[i][0]->Draw("sames");
      hDyMC[i][1]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("#Deltay of matched track-hit pairs"),0.06);
      t1->Draw();
      t1 = GetPaveText(0.15,0.35,0.7,0.75,0.06);
      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",pt_cuts[i],pt_cuts[i+1]));
      t1->SetTextColor(4);
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.6,0.6,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hDyData[0][0],"Data: negative","L");
  leg->AddEntry(hDyData[0][1],"Data: positive","L");
  leg->AddEntry(hDyMC[0][0],"MC: negative","L");
  leg->AddEntry(hDyMC[0][1],"MC: positive","L");
  leg->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_dy_MC_vs_Data.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/RandomMatch_dy_MC_vs_Data.png",run_type));
    }

}

//================================================
void makeHisto(const Int_t save = 1)
{
  // dy distribution from real data
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  TH1F *hDy[4][2];
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hn->GetAxis(5)->SetRange(1+j*2,1+j*2);
	  hn->GetAxis(0)->SetRangeUser(pt_cuts[i]+0.1,pt_cuts[i+1]-0.1);
	  hDy[i][j] = (TH1F*)hn->Projection(2);
	  hDy[i][j]->SetName(Form("Data_TrkDy_pt%1.0f_%1.0f_%s",pt_cuts[i],pt_cuts[i+1],legName_charge[j].Data()));
	}
    }

  // Track y at MTD
  THnSparseF *hTrkProjYZ = (THnSparseF*)f->Get(Form("mhTrkYZatMtd_qa_%s",trigName[kTrigType]));
  hTrkProjYZ->GetAxis(5)->SetRange(2,2);
  TH1F *hTrkYinPt[4][2];
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hTrkProjYZ->GetAxis(4)->SetRange(1+j*2,1+j*2);
	  hTrkProjYZ->GetAxis(0)->SetRangeUser(pt_cuts[i]+0.1,pt_cuts[i+1]-0.1);
	  hTrkYinPt[i][j] = (TH1F*)hTrkProjYZ->Projection(1);
	  hTrkYinPt[i][j]->SetName(Form("Data_TrkYAtMtd_pt%1.0f_%1.0f_%s",pt_cuts[i],pt_cuts[i+1],legName_charge[j].Data()));
	}
    }

  // MTD hit y
  TH2F *hMtdLocaly = (TH2F*)f->Get(Form("mhMthMtdHitLocaly_%s",trigName[kTrigType]));

  if(save)
    {
      TFile *fout = TFile::Open("Rootfiles/Run13.pp500.RandomMatch.root","recreate");
      for(Int_t i=0; i<4; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      hDy[i][j]->Write();
	    }
	}

      for(Int_t i=0; i<4; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      hTrkYinPt[i][j]->Write();
	    }
	}
      hMtdLocaly->Write();
    }
}



