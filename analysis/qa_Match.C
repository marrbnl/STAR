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
  //yzDistribution();
  //eLoss();
}




//================================================
void yzDistribution(const Int_t save = 1)
{
  const char *title[2] = {"z","y"};
  const char *trkname[2] = {"projected tracks","MTD hits"};
  const char *name[2] = {"track","hit"};

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au %s: channel vs backleg of good MTD hits%s",trigName[kTrigType],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.MtdGoodHitMap_%s.png",run_config,trigName[kTrigType]));

  // matched hit map
  TH2F *hMthMtdHitMap = (TH2F*)f->Get(Form("hMthMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMthMtdHitMap,Form("Au+Au %s: channel vs backleg of matched MTD hits%s",trigName[kTrigType],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.MthMtdHitMap_%s.png",run_config,trigName[kTrigType]));
  
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

      c1 = draw2D(hYBL[i],Form("Au+Au %s: local y distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]));
      c2 = draw2D(hZBL[i],Form("Au+Au %s: global z distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]));
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_vs_Module_%s.png",run_config,name[i],trigName[kTrigType]));

      if(i==1)
	{
	  TH2F *h2tmp = (TH2F*)hZBL[i]->Clone(Form("%s_zoomin",hZBL[i]->GetName()));
	  h2tmp->GetXaxis()->SetRangeUser(0,25);
	  c2 = draw2D(h2tmp,Form("Au+Au %s: global z distribution of %s per module;backleg*5+module",trigName[kTrigType],trkname[i]),0.04,kFALSE);
	  if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_vs_Module_zoomin_%s.png",run_config,name[i],trigName[kTrigType]));
	}


      hYAll[i] = (TH1F*)hYBL[i]->ProjectionY(Form("%s_Y_All",trkname[i]));
      hZAll[i] = (TH1F*)hZBL[i]->ProjectionY(Form("%s_Z_All",trkname[i]));     
      hYAll[i]->Scale(1./hYAll[i]->Integral());
      hZAll[i]->Scale(1./hZAll[i]->Integral());

      c1 = draw1D(hYAll[i],Form("Au+Au %s: local y distribution of %s",trigName[kTrigType],trkname[i]),kFALSE,kFALSE);
      c2 = draw1D(hZAll[i],Form("Au+Au %s: global z distribution of %s",trigName[kTrigType],trkname[i]),kFALSE,kFALSE);
      if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_y_%s.png",run_config,name[i],trigName[kTrigType]));
      if(save) c2->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.%s_z_%s.png",run_config,name[i],trigName[kTrigType]));
    }

  THnSparseF *hTrkYZ = (THnSparseF*)f->Get(Form("hTrkYZ_qa_%s",trigName[kTrigType]));
  const char *trkYZ_name[4] = {"start point","outer magnet","MTD (one module)","MTD (multiple modules)"};
  for(Int_t i=0; i<4; i++)
    {
      hTrkYZ->GetAxis(3)->SetRange(i+1,i+1);
      TH2F *h2 = (TH2F*)hTrkYZ->Projection(2,1);
      h2->SetName(Form("hTrkYZ_%s",trkYZ_name[i]));
      c = draw2D(h2,Form("Au+Au %s: #varphi vs global z of %s tracks at %s",trigName[kTrigType],trk_name[trk_index],trkYZ_name[i]),0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.track_phi_vs_z_%s_%s.png",run_config,trkYZ_name[i],trigName[kTrigType]));

      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("hTrkY_%s",trkYZ_name[i]));
      c = draw1D(h1,Form("Au+Au %s: #varphi distribution of %s tracks at %s",trigName[kTrigType],trk_name[trk_index],trkYZ_name[i]),kFALSE,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.track_phi_%s_%s.png",run_config,trkYZ_name[i],trigName[kTrigType]));
    }


  TH1F *hHitZ[5];
  TH2F *hHitZvsBL[5];
  TList *list1 = new TList;
  TString legName1[5];
  for(Int_t i=0; i<5; i++)
    {
      hHitZvsBL[i] = new TH2F(Form("hHitZvsBL_Mod%d",i+1),Form("Au+Au %s: global z distribution of MTD hits in module %d;backleg;z (cm)",trigName[kTrigType],i+1),30,0,30,hZBL[1]->GetNbinsY(),hZBL[1]->GetYaxis()->GetXmin(), hZBL[1]->GetYaxis()->GetXmax());
      for(Int_t j=0; j<30; j++)
	{
	  TH1F *htmp = (TH1F*)hZBL[1]->ProjectionY(Form("hHitZ_BL%d_Mod%d",j+1,i+1),j*5+i+1,j*5+i+1);
	  if(j==0) hHitZ[i] = (TH1F*)htmp->Clone(Form("hHitZ_Mod%d_AllBL",i+1));
	  else     hHitZ[i]->Add(htmp);
	  for(Int_t ibin=1; ibin<=hHitZvsBL[i]->GetNbinsY(); ibin++)
	    {
	      hHitZvsBL[i]->SetBinContent(j+1,ibin,htmp->GetBinContent(ibin));
	      hHitZvsBL[i]->SetBinError(j+1,ibin,htmp->GetBinError(ibin));
	    }
	}
      hHitZ[i]->Sumw2();
      hHitZ[i]->Scale(1./hHitZ[i]->Integral());
      list1->Add(hHitZ[i]);
      legName1[i] = Form("Module %d",i+1);
      c = draw2D(hHitZvsBL[i],"",0.04,kFALSE);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.hit_z_vs_BL_Mod%d_%s.png",run_config,i+1,trigName[kTrigType]));
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


//================================================
void eLoss(const Int_t save = 1)
{
  TCanvas *c = new TCanvas("energy_loss","energy_loss",800,600);
  TF1 *fEloss = new TF1("f2","[0]*exp(-pow([1]/x,[2]))",0.,20);
  fEloss->SetParameters(1.38147e+00,6.08655e-02,5.03337e-01);
  fEloss->SetMinimum(0);
  fEloss->SetTitle(";E_{track} (GeV);#DeltaE (GeV)");
  fEloss->Draw();
 
  TF1 *fEmc = new TF1("fEmc","pol0",0.,100.);
  fEmc->SetParameter(0,0.215);
  fEmc->SetLineColor(2);
  fEmc->Draw("sames");

  TF1 *fCoil = new TF1("fCoil","pol0",0.,100.);
  fCoil->SetParameter(0,0.176);
  fCoil->SetLineColor(4);
  fCoil->Draw("sames");

  TF1 *fMag = new TF1("fMag","[0]*exp(-pow([1]/x,[2]))-[3]",0.,100);
  fMag->SetParameters(1.38147e+00,6.08655e-02,5.03337e-01,fEmc->GetParameter(0),fCoil->GetParameter(0));
  fMag->SetLineColor(6);
  fMag->Draw("sames");

  TLegend *leg = new TLegend(0.6,0.3,0.8,0.6);
  leg->SetHeader("Energy loss");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(fEloss,"Sum","L");
  leg->AddEntry(fEmc,"EMC","L");
  leg->AddEntry(fCoil,"Coil","L");
  leg->AddEntry(fMag,"Magnet","L");
  leg->Draw();

  TPaveText *t1 = GetTitleText("Energy loss of tracks propagated to MTD",0.05);
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/EnergyLoss.png"));
}
