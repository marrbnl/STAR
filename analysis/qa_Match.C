TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

//================================================
void qa_Match()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  Track();
  //DeltaZ();
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

  Double_t pt_cut = 1;
  hTrkDzVsPt->GetXaxis()->SetRangeUser(pt_cut,100);
  TH1F *hMthDz = (TH1F*)hTrkDzVsPt->ProjectionY(Form("hTrkDzVsPt_%s_proj",trigName[kTrigType]));
  Double_t range = 40;
  TF1 *func = new TF1("func","gaus(0)+gaus(3)",-1*range,range);
  func->SetParameters(1000,0,10,100,0,50);
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
  legName[1] = Form("Matached tracks <N> = %2.2f",hMthMtdHitN->GetMean());
  c = drawHistos(list,"Track_multiplicity",Form("Au+Au %s: multiplicity of %s tracks%s;N_{trk};1/N_{evt} dN_{trk}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kFALSE,0,100,kTRUE,1e-8,10,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.6,0.8);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/%s.NMthTrk_%s.png",run_config,trigName[kTrigType]));
}
