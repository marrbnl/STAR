TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "Match.Global.HLT";

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
  yzDistribution();
}


//================================================
void yzDistribution(const Int_t save = 0)
{
  const char *title[2] = {"z","y"};
  const char *trkname[2] = {"projected tracks","MTD hits"};
  const char *name[2] = {"track","hit"};
  
  THnSparseF *hYZ[2];
  hYZ[0] = (THnSparseF*)f->Get(Form("hTrkProjYZ_qa_%s",trigName[kTrigType]));
  hYZ[1] = (THnSparseF*)f->Get(Form("hHitYZ_qa_%s",trigName[kTrigType]));
  TH2F *hYBL[2], *hZBL[2];
  TH1F *hYAll[2], *hZAll[2];
  for(Int_t i=0; i<2; i++)
    {
      hYBL[i] = (TH2F*)hYZ[i]->Projection(0,2);
      hYBL[i]->SetName(Form("%s_YBL",trkname[i]));
      hYBL[i]->GetYaxis()->SetRangeUser(-50,50);

      hZBL[i] = (TH2F*)hYZ[i]->Projection(1,2);
      hZBL[i]->SetName(Form("%s_ZBL",trkname[i]));

      c1 = draw2D(hYBL[i],Form("y distribution of %s per module;backleg*5+module",trkname[i]));
      c2 = draw2D(hZBL[i],Form("z distribution of %s per module;backleg*5+module",trkname[i]));
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));


      hYAll[i] = (TH1F*)hYBL[i]->ProjectionY(Form("%s_Y_All",trkname[i]));
      hZAll[i] = (TH1F*)hZBL[i]->ProjectionY(Form("%s_Z_All",trkname[i]));     
      hYAll[i]->Scale(1./hYAll[i]->Integral());
      hZAll[i]->Scale(1./hZAll[i]->Integral());

      c1 = draw1D(hYAll[i],Form("Local y distribution of %s",trkname[i]),kFALSE,kFALSE);
      c2 = draw1D(hZAll[i],Form("Global z distribution of %s",trkname[i]),kFALSE,kFALSE);
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_%s.png",run_config,name[i],trigName[kTrigType]));
    }

  TH1F *hHitZ[5];
  TList *list1 = new TList;
  TString legName1[5];
  for(Int_t i=0; i<5; i++)
    {
      for(Int_t j=0; j<30; j++)
	{
	  TH1F *htmp = (TH1F*)hZBL[1]->ProjectionY(Form("hHitZ_BL%d_Mod%d",j+1,i+1),j*5+i+1,j*5+i+1);
	  if(j==0) hHitZ[i] = (TH1F*)htmp->Clone(Form("hHitZ_Mod%d_AllBL",i+1));
	  else     hHitZ[i]->Add(htmp);
	}
      hHitZ[i]->Sumw2();
      hHitZ[i]->Scale(1./hHitZ[i]->Integral());
      list1->Add(hHitZ[i]);
      legName1[i] = Form("Module %d",i+1);
    }
  c = drawHistos(list1,list1->At(0)->GetName(),Form("Au+Au %s: z distribution of MTD hits in each module;z (cm)",trigName[kTrigType]),kFALSE,0,100,kTRUE,1e-6,0.018,kFALSE,kTRUE,legName1,kTRUE,"",0.6,0.8,0.68,0.88,kFALSE);
  for(Int_t i=0; i<6; i++)
    {
      TLine *line = GetLine(-217.5+i*87, 0, -217.5+i*87, 0.0125, 2, 2, 2);
      line->Draw();
    }
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.hit_z_per_Module_%s.png",run_config,trigName[kTrigType]));

      
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
