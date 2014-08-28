TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void ana_DCA()
{
  gStyle->SetOptStat(0);
  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.%s.DCA1cm.root",run_config),"read");

  tracks();
}

//================================================
void tracks(const Int_t save = 1)
{
  TFile *f1 = TFile::Open(Form("~/Work/STAR/analysis/output/jpsi.AuAu200.Run14.PrimTrk.ClosePrimVtx.root"),"read");
  TH1F *hMthTrkPt[2];
  TList *list = new TList;
  TFile *ftmp = 0;
  for(Int_t i=0; i<2; i++)
    {
      if(i==0) ftmp = f1;
      if(i==1) ftmp = f;
      hMthTrkPt[i] = (TH1F*)ftmp->Get(Form("hMthTrkPt_%s",trigName[kTrigType]));
      hMthTrkPt[i]->SetName(Form("%s_%d",hMthTrkPt[i]->GetName(),i));

      TH1F *hStat = (TH1F*)ftmp->Get("hEventStat");
      scaleHisto( hMthTrkPt[i], hStat->GetBinContent(kTrigType+7), 1, kTRUE);
      list->Add(hMthTrkPt[i]);
    }
  TString legName[2] = {"gDCA < 3 cm","gDCA < 1 cm"};
  c = drawHistos(list,"MthTrkPt",Form("Au+Au %s: p_{T} distribution of matched %s tracks;p_{T} (GeV/c);1/N dN/dp_{T}",trigName[kTrigType],trk_name[trk_index]),kFALSE,0,100,kFALSE,1e-1,1.6e6,kTRUE,kTRUE,legName,kTRUE,"",0.6,0.8,0.65,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_DCA/%s.MthTrkPt_%s.png",run_config,trigName[kTrigType]));

  TH1F *hRatio = (TH1F*)hMthTrkPt[1]->Clone("hMthTrkPt_Ratio");
  hRatio->Divide(hMthTrkPt[0]);
  c = draw1D(hRatio,Form("Au+Au %s: ratio of p_{T} distribution of matched %s tracks;p_{T} (GeV/c)",trigName[kTrigType],trk_name[trk_index]));
  TPaveText *t1 = GetPaveText(0.3,0.4,0.2,0.4);
  t1->AddText("ratio = #frac{1/N dN/dp_{T} (gDCA<1cm)}{1/N dN/dp_{T} (gDCA<3cm)}");
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_DCA/%s.MthTrkPt_ratio_%s.png",run_config,trigName[kTrigType]));
}
