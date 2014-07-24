TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "GlobalTrk.NoDCA.HLT";

//================================================
void qa_Match()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  //Track();
  //DeltaZ();
  //qualityCuts();
  randomMatch();
  //MonteCarlo();
  //makeHisto();
}

//================================================
void MonteCarlo(const Int_t save = 0)
{
  TH1F *hDy = new TH1F("hDy_MC","Monte Carlo: #Deltay distribution; #Deltay (cm)",200,-100,100);
  Double_t center[12] = {-24.2,-19.8,-15.4,-11,-6.6,-2.2,2.2,6.6,11,15.4,19.8,24.2};
  const Int_t nExpr = 1e4;

  for(Int_t i=0; i<nExpr; i++)
    {
      Double_t hit_y = center[(Int_t)(gRandom->Rndm()*12)];
      Double_t trk_y = (gRandom->Rndm()-0.5)*2*26.4;
      hDy->Fill(trk_y-hit_y);
    }
  draw1D(hDy);
}

//================================================
void makeHisto(const Int_t save = 1)
{
  const char *title[2] = {"z","y"};
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/AuAu200.Run14.jpsi.RotateMTD.HLT.root"),"read");
  THnSparseF *hn = (THnSparseF*)f1->Get("hTrkDzDy_di-muon");
  TH2F *hRes_All[2];
  TH2F *hRes_BL[2][30];
  TH2F *hRes_Mod[2][30][5];
  for(Int_t i=0; i<2; i++)
    {
      hRes_All[i] = (TH2F*)hn->Projection(i+1,0);
      hRes_All[i]->SetName(Form("hD%sVsPt_All",title[i]));
      for(Int_t j=0; j<30; j++)
	{
	  hn->GetAxis(3)->SetRange(j*60+1, j*60+60);
	  hRes_BL[i][j] = (TH2F*)hn->Projection(i+1,0);
	  hRes_BL[i][j]->SetName(Form("hD%sVsPt_BL%d",title[i],j+1));
	  for(Int_t k=0; k<5; k++)
	    {
	      hn->GetAxis(3)->SetRange(j*60+k*12+1, j*60+k*12+12);
	      hRes_Mod[i][j][k] = (TH2F*)hn->Projection(i+1,0);
	      hRes_Mod[i][j][k]->SetName(Form("hD%sVsPt_BL%d_Mod%d",title[i],j+1,k+1));
	      hn->GetAxis(3)->SetRange(0,-1);
	    }
	  hn->GetAxis(3)->SetRange(0,-1);
	}
    }

  if(save)
    {
      TFile *fout = TFile::Open("Rootfiles/AuAu200.RotateMTD.root","recreate");
      for(Int_t i=0; i<2; i++)
	{
	  hRes_All[i]->Write();
	}
      for(Int_t i=0; i<2; i++)
	{
	  for(Int_t j=0; j<30; j++)
	    {
	      hRes_BL[i][j]->Write();
	    }
	}

      for(Int_t i=0; i<2; i++)
	{
	  for(Int_t j=0; j<30; j++)
	    {
	      for(Int_t k=0; k<5; k++)
		{
		  hRes_Mod[i][j][k]->Write();
		}
	    }
	}
      
    }
}

//================================================
void randomMatch(const Int_t save = 0)
{
  const char *title[2] = {"z","y"};
  TString legName[2] = {"Standard","Rotated"};
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/AuAu200.Run14.jpsi.RotateMTD.HLT.root"),"read");

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f1->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au di-muon: channel vs backleg of good MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.Rotate_MtdGoodHitMap_%s.png",run_config,trigName[kTrigType]));

  // matched hit map
  TH2F *hMthMtdHitMap = (TH2F*)f1->Get(Form("hMthMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMthMtdHitMap,Form("Au+Au di-muon: channel vs backleg of matched MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.Rotate_MthMtdHitMap_%s.png",run_config,trigName[kTrigType]));

  TFile *f2 = TFile::Open("Rootfiles/AuAu200.RotateMTD.root","read");
  const Double_t pt_cut = 1;
  TCanvas *cRes[2][6];
  for(Int_t i=0; i<2; i++)
    {
      for(Int_t j=0; j<30; j++)
	{
	  if(j%5==0)
	    {
	      cRes[i][j/5] = new TCanvas(Form("D%s_Backleg%d-%d",title[i],j+1,j+5),Form("D%s_Backleg%d-%d",title[i],j+1,j+5),1350,800);
	      cRes[i][j/5]->Divide(5,5);
	    }
	  for(Int_t k=0; k<5; k++)
	    {
	      cRes[i][j/5]->cd(j%5*5+k+1);
	      TH2F *h2 = (TH2F*)f2->Get(Form("hD%sVsPt_BL%d_Mod%d",title[i],j+1,k+1));
	      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("%s_pro",h2->GetName()),h2->GetXaxis()->FindFixBin(pt_cut+0.1),-1);
	      h1->Draw();
	    }
	}
    }

  TH2F *hResi[2][2];
  hResi[0][0] = (TH2F*)f->Get("hTrkDz_di-muon");
  hResi[0][1] = (TH2F*)f2->Get("hDzVsPt_All");
  hResi[1][0] = (TH2F*)f->Get("hTrkDy_di-muon");
  hResi[1][1] = (TH2F*)f2->Get("hDyVsPt_All");
  TList *list[2];

  for(Int_t i=0; i<2; i++)
    {
      list[i] = new TList;
      for(Int_t j=0; j<2; j++)
	{
	  hResi[i][j]->SetName(Form("%s_%s",hResi[i][j]->GetName(),legName[j].Data()));
	  draw2D(hResi[i][j]);
	  TH1F *h = (TH1F*)hResi[i][j]->ProjectionY(Form("%s_pro",hResi[i][j]->GetName()),hResi[i][j]->GetXaxis()->FindFixBin(pt_cut+0.1),-1);
	  h->Sumw2();
	  h->Scale(1./h->Integral());
	  list[i]->Add(h);
	}
      c = drawHistos(list[i],list[i]->At(0)->GetName(),Form("Au+Au %s: #Delta%s distribution using %s tracks%s;#Delta%s (cm)",trigName[kTrigType],title[i],trk_name[trk_index],hlt_name[hlt_index],title[i]),kFALSE,0,100,kTRUE,0,0.02,kFALSE,kTRUE,legName,kTRUE,Form("p_{T} > %1.1f GeV/c",pt_cut),0.6,0.8,0.6,0.8,kFALSE);
      if(save)
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.Rotate_D%s_Pt%1.0fGeV_%s.png",run_config,title[i],pt_cut,trigName[kTrigType]));
    }
      
}

//================================================
void qualityCuts(const Int_t save = 0)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.13); 

  THnSparseF *hn = (THnSparseF*)f->Get(Form("hTrkHitDz_qa_%s",trigName[kTrigType]));
  Double_t pt_cut = 2;
  hn->GetAxis(0)->SetRangeUser(pt_cut+0.1,100);

  TH1F *hdz[3];
  Double_t nsigma_low[3] = {-100,-1,0};
  Double_t nsigma_hi[3]  = {100, 3, 3};
  Double_t range = 40, range2 = 20;
  for(Int_t i=0; i<3; i++)
    {
      hn->GetAxis(4)->SetRangeUser(nsigma_low[i]+0.1,nsigma_hi[i]-0.1);
      hdz[i] = (TH1F*)hn->Projection(5);
      hdz[i]->SetName(Form("hdz_nsigma_cut%d",i));
      if(i==0) hdz[i]->SetTitle(Form("Au+Au %s: #Deltaz distribution w/o n#sigma_{#pi} cut%s;#Deltaz (cm)",trigName[kTrigType],hlt_name[hlt_index]));
      else     hdz[i]->SetTitle(Form("Au+Au %s: #Deltaz distribution w/ %1.0f<n#sigma_{#pi}<%1.0f%s;#Deltaz (cm)",trigName[kTrigType],nsigma_low[i],nsigma_hi[i],hlt_name[hlt_index]));

      TF1 *func = new TF1(Form("func_%d",i),"gaus(0)+gaus(3)",-1*range,range);
      func->SetParameters(10000,0,10,100,0,100);
      //func->SetParameters(1000,0,10,100,0,50);
      hdz[i]->Fit(func,"R0");
      hdz[i]->GetYaxis()->SetNdivisions(505);
      c = draw1D(hdz[i]);
      TF1 *func1 = new TF1(Form("func1_%d",i),"gaus",-1*range,range);
      func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
      func->SetLineColor(2);
      func->Draw("sames");
      func1->SetLineColor(4);
      func1->Draw("sames");
      TPaveText *t1 = GetPaveText(0.2,0.3,0.65,0.85,0.04);
      t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range,range));
      t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*range,range))-(Int_t)(func1->Integral(-1*range,range))));
      t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*range,range))));
      t1->AddText(Form("S/B ~ %1.2f",(func->Integral(-1*range,range)-func1->Integral(-1*range,range))/(func1->Integral(-1*range,range))));
      t1->Draw();
      t1 = GetPaveText(0.2,0.3,0.4,0.6,0.04);
      t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range2,range2));
      t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*range2,range2))-(Int_t)(func1->Integral(-1*range2,range2))));
      t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*range2,range2))));
      t1->AddText(Form("S/B ~ %1.2f",(func->Integral(-1*range2,range2)-func1->Integral(-1*range2,range2))/(func1->Integral(-1*range2,range2))));
      t1->Draw();
      t1 = GetPaveText(0.72,0.8,0.4,0.55,0.04);
      t1->AddText(Form("%s tracks",trk_name[trk_index]));
      t1->AddText(Form("p_{T,trk} > %1.1f GeV/c",pt_cut));
      t1->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.FitDz_nsigmaCut%d_Pt%1.0fGeV_%s.png",run_config,i,pt_cut,trigName[kTrigType]));
    }
}

//================================================
void DeltaZ(const Int_t save = 0)
{
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.15);                
  gStyle->SetStatH(0.13); 

  TH2F *hTrkDzVsPt = (TH2F*)f->Get(Form("hTrkDz_%s",trigName[kTrigType]));
  hTrkDzVsPt->GetXaxis()->SetRangeUser(0,6);
  hTrkDzVsPt->GetYaxis()->SetRangeUser(-100,100);
  c = draw2D(hTrkDzVsPt,Form("Au+Au %s: #Deltaz of matched %s track-hit pairs%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.TrkPtDz_%s.png",run_config,trigName[kTrigType]));

  Double_t pt_cut = 2;
  hTrkDzVsPt->GetXaxis()->SetRangeUser(pt_cut,100);
  TH1F *hMthDz = (TH1F*)hTrkDzVsPt->ProjectionY(Form("hTrkDzVsPt_%s_proj",trigName[kTrigType]));
  Double_t range = 40;
  TF1 *func = new TF1("func","gaus(0)+gaus(3)",-1*range,range);
  func->SetParameters(1000,0,10,100,0,100);
  //func->SetParameters(1000,0,10,100,0,50);
  hMthDz->Fit(func,"R0");
  hMthDz->GetYaxis()->SetNdivisions(505);
  c = draw1D(hMthDz,Form("Au+Au %s: #Deltaz of matched %s track-hit pairs (p_{T}>%1.1f GeV/c);#Deltaz (cm)",trigName[kTrigType],trk_name[trk_index],pt_cut));
  TF1 *func1 = new TF1("func1","gaus",-1*range,range);
  func1->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
  func->SetLineColor(2);
  func->Draw("sames");
  func1->SetLineColor(4);
  func1->Draw("sames");
  TPaveText *t1 = GetPaveText(0.2,0.3,0.7,0.85,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range,range));
  t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*range,range))-(Int_t)(func1->Integral(-1*range,range))));
  t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*range,range))));
  t1->Draw();
  range = 20;
  t1 = GetPaveText(0.2,0.3,0.45,0.6,0.04);
  t1->AddText(Form("#Deltaz ~ [-%1.0f,%1.0f]",range,range));
  t1->AddText(Form("signal ~ %d",(Int_t)(func->Integral(-1*range,range))-(Int_t)(func1->Integral(-1*range,range))));
  t1->AddText(Form("background ~ %d",(Int_t)(func1->Integral(-1*range,range))));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.FitDz_Pt%1.0fGeV_%s.png",run_config,pt_cut,trigName[kTrigType]));
 
}

//================================================
void Track(const Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");

  // track multiplicity
  TH1F *hNTrk       = (TH1F*)f->Get(Form("hNTrk_%s",trigName[kTrigType]));
  TH1F *hMthMtdHitN = (TH1F*)f->Get(Form("hMthMtdHitN_%s",trigName[kTrigType]));
  scaleHisto( hNTrk,    hStat->GetBinContent(kTrigType+7), 1);
  scaleHisto( hMthMtdHitN, hStat->GetBinContent(kTrigType+7), 1);
  TList *list = new TList;
  list->Add(hNTrk);
  list->Add(hMthMtdHitN);
  TString legName[2];
  legName[0] = Form("Good tracks <N> = %2.2f",hNTrk->GetMean());
  legName[1] = Form("Matched tracks <N> = %2.2f",hMthMtdHitN->GetMean());
  c = drawHistos(list,"Track_multiplicity",Form("Au+Au %s: multiplicity of %s tracks%s;N_{trk};1/N_{evt} dN_{trk}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kFALSE,0,100,kTRUE,1e-8,10,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.6,0.8);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.NMthTrk_%s.png",run_config,trigName[kTrigType]));
}
