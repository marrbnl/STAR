TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "DeltaZ.new.";
//const char *run_config = "";
//const char *run_config = "CutRanking.";
const Bool_t iPico = 0;
const int year = 2013;
TString run_cfg_name;

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

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run13.pp500.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run13.pp500.jpsi.%sroot",run_config),"read");
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) f = TFile::Open(Form("./output/Pico.Run14.AuAu200.jpsi.%sroot",run_config),"read");
      else      f = TFile::Open(Form("./output/Run14.AuAu200.jpsi.%sroot",run_config),"read");
    }
  run_cfg_name = Form("%s",run_config);
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //DeltaZVsPt();
  DeltaZVsPos();
  //DeltaY();
  //makeHistos();
  //fitDzInMod();

}

//================================================
void fitDzInMod(const int save = 0)
{
  TFile *fin = TFile::Open(Form("Rootfiles/%s%s.TrkMthResidual.root",run_cfg_name.Data(),run_type));
  TH2F *hDzVsModBL = (TH2F*)fin->Get(Form("hTrkDzVsModBL_%s",trigName[kTrigType]));
  c = draw2D(hDzVsModBL,Form("%s: #Deltaz of matched track-hit pairs;module;#Deltaz (cm)",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_ModInBL_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_ModInBL_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  TCanvas *canvas[3];
  for(int i=0; i<3; i++)
    {
      canvas[i] = new TCanvas(Form("canvas_%d",i),Form("canvas_%d",i),1100,800);
      canvas[i]->Divide(5,5);
    }

  Int_t counter=0;
  TH1F *hDz[30][5];
  TH1F *hDzMean[5];
  TH1F *hDzSigma[5];
  for(int j=0; j<5; j++)
    {
      hDzMean[j] = new TH1F(Form("hTrkDzMean_Mod%d_%s",j+1,trigName[kTrigType]),Form("%s: mean of #Deltaz from fit;module;<#Deltaz> (cm)",trigName[kTrigType]),150,0,150);
      hDzSigma[j] = new TH1F(Form("hDzSigma_Mod%d_%s",j+1,trigName[kTrigType]),Form("%s: sigma of #Deltaz from fit;module;#sigma(#Deltaz) (cm)",trigName[kTrigType]),150,0,150);
    }

  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  hDz[i][j] = (TH1F*)fin->Get(Form("hTrkDz_BL%d_Mod%d_%s",i+1,j+1,trigName[kTrigType]));
	  if(i==6 && j==4) counter++;
	  if(hDz[i][j]->GetEntries()<=0) continue;
	  canvas[counter/25]->cd(counter%25+1);
	  hDz[i][j]->GetXaxis()->SetRangeUser(-60,100);
	  hDz[i][j]->SetTitle("");
	  hDz[i][j]->SetLineColor(1);
	  hDz[i][j]->SetLineWidth(2);
	  hDz[i][j]->Draw();
	  TF1 *func = new TF1(Form("Fit_%s",hDz[i][j]->GetName()),"gaus",-20,25);
	  hDz[i][j]->Fit(func,"IRQ0");
	  func->SetLineColor(2);
	  func->SetLineStyle(2);
	  func->Draw("sames");
	  TPaveText *t1 = GetPaveText(0.15,0.35,0.6,0.8,0.1);
	  t1->AddText(Form("BL %d",i+1));
	  t1->AddText(Form("Mod %d",j+1));
	  t1->Draw();
	  counter++;

	  hDzMean[j]->SetBinContent(i*5+j+1,func->GetParameter(1));
	  hDzMean[j]->SetBinError(i*5+j+1,func->GetParError(1));
	  hDzSigma[j]->SetBinContent(i*5+j+1,func->GetParameter(2));
	  hDzSigma[j]->SetBinError(i*5+j+1,func->GetParError(2));
	}
    }

  if(save) 
    {
      for(int i=0; i<3; i++)
	{
	  canvas[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDz_vs_Mod_%d_%s.pdf",run_type,run_cfg_name.Data(),i,trigName[kTrigType]));
	  canvas[i]->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDz_vs_Mod_%d_%s.png",run_type,run_cfg_name.Data(),i,trigName[kTrigType]));
	}
    }

  TList *list = new TList;
  list->Clear();
  TString legName[5];
  for(Int_t j=0; j<5; j++)
    {
      legName[j] = Form("Module %d",j+1);
      list->Add(hDzMean[j]);
    }
  c = drawHistos(list,"TrkDzMean","",kTRUE,0,150,kTRUE,-4,14,kFALSE,kTRUE,legName,kTRUE,"",0.45,0.6,0.6,0.88,kTRUE);
  TLine *line = GetLine(0,0,150,0,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDzMean_vs_Mod_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDzMean_vs_Mod_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  for(Int_t j=0; j<5; j++)
    {
      list->Add(hDzSigma[j]);
    }
  c = drawHistos(list,"TrkDzSigma","",kTRUE,0,150,kTRUE,6,19,kFALSE,kTRUE,legName,kTRUE,"",0.2,0.35,0.18,0.42,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDzSigma_vs_Mod_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sFitDzSigma_vs_Mod_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }
}



//================================================
void DeltaZVsPos(const Int_t save = 0)
{
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  TList *list = new TList;

  // dz vs BL
  TH2F *hTrkDzVsBL = (TH2F*)hn->Projection(1,3);
  c = draw2D(hTrkDzVsBL,Form("%s: #Deltaz of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_BL_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_BL_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  list->Clear();
  TString legName[30];
  TH1F *hTrkDzInBL[30];
  Int_t counter = 0;
  for(Int_t i=0; i<30; i++)
    {
      hTrkDzInBL[i] = (TH1F*)hTrkDzVsBL->ProjectionY(Form("hDeltaZ_BL%d",i+1),i+1,i+1);
      if(hTrkDzInBL[i]->GetEntries()>0)
	{
	  legName[counter] = Form("Module %d",i+1);
	  hTrkDzInBL[i]->SetLineColor(color[counter]);
	  list->Add(hTrkDzInBL[i]);
	  counter ++;
	}
    }
  c = drawHistos(list,"TrkDzInBL",Form("%s: #Deltaz of matched track-hit pairs in backleg;#Deltaz (cm)",trigName[kTrigType]),kTRUE,-100,100,kTRUE,0,1.2*hTrkDzInBL[1]->GetMaximum(),kFALSE,kTRUE,legName,kTRUE,"",0.15,0.25,0.2,0.88,kFALSE,0.04,0.04,kFALSE,1,kTRUE,kFALSE);
  TLine *line = GetLine(0,0,0,1.1*hTrkDzInBL[1]->GetMaximum(),1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_BL_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_BL_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  // dz vs Mod
  TH2F *hTrkDzVsMod = (TH2F*)hn->Projection(1,4);
  c = draw2D(hTrkDzVsMod,Form("%s: #Deltaz of matched track-hit pairs",trigName[kTrigType]));
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_Mod_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_vs_Mod_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  TH1F *hMthMod = (TH1F*)hTrkDzVsMod->ProjectionX("hMthMod");
  hMthMod->Sumw2();
  hMthMod->Scale(1./hMthMod->Integral());
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("mhMtdHitMap_%s",trigName[kTrigType]));
  TH1F *htmp = (TH1F*)hMtdHitMap->ProjectionY("hHitMod_finebin");
  htmp->Rebin(12);
  TH1F *hMtdHitMod = new TH1F(Form("hMtdHitMod_%s",trigName[kTrigType]),"# of MTD hits per module;module",5,1,6);
  for(int i=0; i<hMtdHitMod->GetNbinsX(); i++)
    {
      hMtdHitMod->SetBinContent(i+1,htmp->GetBinContent(i+1));
      hMtdHitMod->SetBinError(i+1,htmp->GetBinError(i+1));
    }
  hMtdHitMod->Scale(1./hMtdHitMod->Integral());
  list->Clear();
  list->Add(hMthMod);
  list->Add(hMtdHitMod);
  TString legName3[2] = {"Matched good hits","All good hits"};
  c = drawHistos(list,"MtdHitMod",Form("%s: MTD hits per module;module;probability",trigName[kTrigType]),kFALSE,0,5,kTRUE,0,0.5,kFALSE,kTRUE,legName3,kTRUE,"",0.15,0.25,0.6,0.88,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sCompMtdHitMod_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sCompMtdHitMod_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
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
  c = drawHistos(list,"TrkDzInMod",Form("%s: #Deltaz of matched track-hit pairs in module;#Deltaz (cm)",trigName[kTrigType]),kTRUE,-100,100,kTRUE,0,1.2*hTrkDzInMod[3]->GetMaximum(),kFALSE,kTRUE,legName2,kTRUE,"",0.15,0.25,0.6,0.88,kTRUE);
  TLine *line = GetLine(0,0,0,hTrkDzInMod[3]->GetMaximum()*1.05,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_Mod_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/ana_Match/%sDeltaZ_in_Mod_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
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
void makeHistos()
{
  TH2F *hDzVsModBL;
  TH1F *hDz[30][5];

  // dz vs mod in each BL
  THnSparseF *hn = (THnSparseF*)f->Get(Form("mhTrkDzDy_%s",trigName[kTrigType]));
  Int_t nbins = hn->GetAxis(1)->GetNbins();
  Double_t xmin = hn->GetAxis(1)->GetXmin();
  Double_t xmax = hn->GetAxis(1)->GetXmax();

  hDzVsModBL = new TH2F(Form("hTrkDzVsModBL_%s",trigName[kTrigType]),Form("%s: #Deltaz of matched track-hit pairs;module;#Deltaz (cm)",trigName[kTrigType]),150,0,150,nbins,xmin,xmax);
  for(int i=0; i<30; i++)
    {
      hn->GetAxis(3)->SetRange(i+1,i+1);
      for(int j=0; j<5; j++)
	{
	  hn->GetAxis(4)->SetRange(j+1,j+1);
	  hDz[i][j] = (TH1F*)hn->Projection(1);
	  hDz[i][j]->SetName(Form("hTrkDz_BL%d_Mod%d_%s",i+1,j+1,trigName[kTrigType]));
	  if(hDz[i][j]->GetEntries()<=0) continue;

	  for(int ibin=1; ibin<=hDz[i][j]->GetNbinsX(); ibin++)
	    {
	      hDzVsModBL->SetBinContent(i*5+j+1,ibin,hDz[i][j]->GetBinContent(ibin));
	      hDzVsModBL->SetBinError(i*5+j+1,ibin,hDz[i][j]->GetBinError(ibin));
	    }
	  hn->GetAxis(4)->SetRange(0,-1);
	}
      hn->GetAxis(3)->SetRange(0,-1);
    }
  TFile *fout = TFile::Open(Form("Rootfiles/%s%s.TrkMthResidual.root",run_cfg_name.Data(),run_type),"update");
  hDzVsModBL->Write();
  for(int i=0; i<30; i++)
    {
      for(int j=0; j<5; j++)
	{
	  hDz[i][j]->Write();
	}
    }
  fout->Close();
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


