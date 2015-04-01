TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void ana_Match()
{						
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  TString cut_name = run_config;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  //f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  f = TFile::Open("./output/Run13.pp500.jpsi.CutRanking.root","read");

  //DeltaZVsPt();
  //DeltaZVsPos();
  DeltaY();

}

//================================================
void DeltaZVsPos(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  TH2F *hTrkDzVsBL = (TH2F*)hn->Projection(1,3);
  c = draw2D(hTrkDzVsBL,Form("%s: #Deltaz of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_BL_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_BL_%s.png",run_type,trigName[kTrigType]));
    }

  TList *list = new TList;
  TString legName[30];
  TH1F *hTrkDzInBL[30];
  Int_t counter = 0;
  for(Int_t i=0; i<30; i++)
    {
      hTrkDzInBL[i] = (TH1F*)hTrkDzVsBL->ProjectionY(Form("hDeltaZ_BL%d",i+1),i+1,i+1);
      if(hTrkDzInBL[i]->GetEntries()>0)
	{
	  legName[counter] = Form("Module %d",i+1);
	  hTrkDzInBL[i]->SetLineColor(counter+1);
	  list->Add(hTrkDzInBL[i]);
	  counter ++;
	}
    }
  c = drawHistos(list,"TrkDzInBL",Form("%s: #Deltaz of matched track-hit pairs in backleg;#Deltaz (cm)",trigName[kTrigType]),kTRUE,-100,100,kTRUE,0,2e4,kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.2,0.88,kFALSE,0.04,0.04,kFALSE,1,kTRUE,kFALSE);
  TLine *line = GetLine(0,0,0,1.8e4,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_in_BL_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_in_BL_%s.png",run_type,trigName[kTrigType]));
    }

  TH2F *hTrkDzVsMod = (TH2F*)hn->Projection(1,4);
  c = draw2D(hTrkDzVsMod,Form("%s: #Deltaz of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_Mod_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_Mod_%s.png",run_type,trigName[kTrigType]));
    }

  list->Clear();
  TString legName2[5];
  TH1F *hTrkDzInMod[5];
  for(Int_t i=0; i<5; i++)
    {
      hTrkDzInMod[i] = (TH1F*)hTrkDzVsMod->ProjectionY(Form("hDeltaZ_Mod%d",i+1),i+1,i+1);
      legName2[i] = Form("Module %d",i+1);
      list->Add(hTrkDzInMod[i]);
    }
  c = drawHistos(list,"TrkDzInMod",Form("%s: #Deltaz of matched track-hit pairs in module;#Deltaz (cm)",trigName[kTrigType]),kTRUE,-100,100,kTRUE,0,8e4,kFALSE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.6,0.88,kTRUE);
  TLine *line = GetLine(0,0,0,7e4,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_in_Mod_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_in_Mod_%s.png",run_type,trigName[kTrigType]));
    }
}

//================================================
void DeltaZVsPt(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  TH2F *hTrkDzVsPt = (TH2F*)hn->Projection(1,0);
  c = draw2D(hTrkDzVsPt,Form("%s: #Deltaz of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  Double_t pt_cut = 1;
  hTrkDzVsPt->GetXaxis()->SetRangeUser(pt_cut+0.1,100);
  TH1F *hMthDz = (TH1F*)hTrkDzVsPt->ProjectionY(Form("hTrkDzVsPt_%s_proj",trigName[kTrigType]));
  hMthDz->SetTitle(Form("%s: #Deltaz of matched track-hit pairs (p_{T}>%1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],pt_cut));
  TH1F *hClone = (TH1F*)hMthDz->Clone(Form("%s_clone",hMthDz->GetName()));
  c = draw1D(hClone,"",kFALSE,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaZ_%s.png",run_type,trigName[kTrigType]));
    }

  Double_t range = 50;
  TF1 *func = new TF1("func","gaus(0)+gaus(3)",-1*range,range);
  func->SetParameters(10000,0,10,1000,0,40);
  c = FitDeltaZ(hMthDz,func,range,20.);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/FitDz_Pt%1.0f_%s.pdf",run_type,pt_cut,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/FitDz_Pt%1.0f_%s.png",run_type,pt_cut,trigName[kTrigType]));
    }

  // pt dependence
  Double_t pt_cuts[5] = {1,2,3,5,20};
  for(Int_t i=0; i<4; i++)
    {
      hTrkDzVsPt->GetXaxis()->SetRangeUser(pt_cuts[i]+0.1,pt_cuts[i+1]-0.1);
      TH1F *htmp = (TH1F*)hTrkDzVsPt->ProjectionY(Form("hTrkDz_pt%1.0f-%1.0f_%s",pt_cuts[i],pt_cuts[i+1],trigName[kTrigType]));
      htmp->SetTitle(Form("%s: #Deltaz of matched track-hit pairs (%1.0f < p_{T} < %1.0f GeV/c);#Deltaz (cm)",trigName[kTrigType],pt_cuts[i],pt_cuts[i+1]));

      TF1 *func = new TF1(Form("func_pt%1.0f-%1.0f",pt_cuts[i],pt_cuts[i+1]),"gaus(0)+gaus(3)",-1*range,range);
      if(i==0) func->SetParameters(100,0,100,1000,0,10);
      if(i==1) func->SetParameters(1000,0,15,1000,0,60);
      if(i==2) func->SetParameters(1000,0,15,1000,0,60);
      if(i==3) func->SetParameters(1000,0,60,1000,0,15);
      c = FitDeltaZ(htmp,func,range,20.);
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/FitDz_Pt%1.0f_%1.0f_%s.pdf",run_type,pt_cuts[i],pt_cuts[i+1],trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/FitDz_Pt%1.0f_%1.0f_%s.png",run_type,pt_cuts[i],pt_cuts[i+1],trigName[kTrigType]));
	}
    }
 
}

//================================================
void DeltaY(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  TH2F *hTrkDyVsPt = (TH2F*)hn->Projection(2,0);
  c = draw2D(hTrkDyVsPt,Form("%s: #Deltay of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaY_vs_pt_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaY_vs_pt_%s.png",run_type,trigName[kTrigType]));
    }

  c = new TCanvas("hDy_TrkPtBin","hDy_TrkPtBin",1200,650);
  c->Divide(2,2);
  Double_t pt_cuts[5] = {1,2,4,10,20};
  TH1F *hDy[4][2];
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<2; j++)
	{
	  hn->GetAxis(5)->SetRange(1+j*2,1+j*2);
	  hn->GetAxis(0)->SetRangeUser(pt_cuts[i]+0.1,pt_cuts[i+1]-0.1);
	  hDy[i][j] = (TH1F*)hn->Projection(2);
	  hDy[i][j]->SetName(Form("hTrkDy_pt%1.0f_%1.0f_%d",pt_cuts[i],pt_cuts[i+1],j));
	  hDy[i][j]->SetLineColor(color[j]);
	  hDy[i][j]->SetMaximum(1.3*hDy[i][j]->GetMaximum());
	  hDy[i][j]->SetTitle("");
	}
	 
      c->cd(i+1);
      hDy[i][0]->Draw();
      hDy[i][1]->Draw("sames");
      TPaveText *t1 = GetTitleText(Form("%s: #Deltay of matched track-hit pairs",trigName[kTrigType]),0.06);
      t1->Draw();
      t1 = GetPaveText(0.15,0.35,0.7,0.75,0.06);
      t1->AddText(Form("%1.0f < p_{T} < %1.0f GeV/c",pt_cuts[i],pt_cuts[i+1]));
      t1->SetTextColor(4);
      t1->Draw();
    }
  c->cd(1);
  TLegend *leg = new TLegend(0.6,0.63,0.8,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hDy[0][0],"Negative","L");
  leg->AddEntry(hDy[0][1],"Positive","L");
  leg->Draw();

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaY_InPtBin_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/DeltaY_InPtBin_%s.png",run_type,trigName[kTrigType]));
    }
}

//================================================
TCanvas *FitDeltaZ(TH1 *histo, TF1 *func, const Double_t range1 = 50., Double_t range2 = 20.)
{
  histo->Fit(func,"R0");
  histo->GetYaxis()->SetNdivisions(505);
  c = draw1D(histo);
  TF1 *func1 = new TF1("func1","gaus",-1*range1,range1);
  func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
  func->SetLineColor(2);
  func->Draw("sames");
  func1->SetLineColor(4);
  func1->Draw("sames");
  TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.78,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range1,range1));
  t1->AddText(Form("S/B ~ %1.1f:1",((func->Integral(-1*range1,range1))-(func1->Integral(-1*range1,range1)))/(func1->Integral(-1*range1,range1))));
  t1->Draw();
  t1 = GetPaveText(0.2,0.3,0.47,0.6,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range2,range2));
  t1->AddText(Form("S/B ~ %1.1f:1",((func->Integral(-1*range2,range2))-(func1->Integral(-1*range2,range2)))/(func1->Integral(-1*range2,range2))));
  t1->Draw();
  return c;
}
