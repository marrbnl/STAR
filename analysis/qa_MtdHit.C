TFile *f;

//================================================
void qa_MtdHit(const Int_t save = 0)
{
  gStyle->SetOptStat(0);

  const char *cut = "GlobalTrk.HLT";
  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",cut),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");

  TString cut_name = cut;
  Int_t hlt_index = 0;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;

  // hit multiplicity
  TH1F *hMtdHitN    = (TH1F*)f->Get(Form("hMtdHitN_%s",trigName[kTrigType]));
  TH1F *hMthMtdHitN = (TH1F*)f->Get(Form("hMthMtdHitN_%s",trigName[kTrigType]));
  scaleHisto( hMtdHitN,    hStat->GetBinContent(kTrigType+7), 1);
  scaleHisto( hMthMtdHitN, hStat->GetBinContent(kTrigType+7), 1);
  TList *list = new TList;
  list->Add(hMtdHitN);
  list->Add(hMthMtdHitN);
  TString legName[2];
  legName[0] = Form("All hits <N> = %2.2f",hMtdHitN->GetMean());
  legName[1] = Form("Matached hits <N> = %2.2f",hMthMtdHitN->GetMean());
  c = drawHistos(list,"MTD_hit_multiplicity",Form("Au+Au di-muon: multiplicity of good MTD hits%s;N_{hit};1/N_{evt} dN_{hit}",hlt_name[hlt_index]),kFALSE,0,100,kTRUE,1e-8,10,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.7,0.6,0.8);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.MtdGoodHitN_%s.png",cut,trigName[kTrigType]));

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au di-muon: channel vs backleg of good MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.MtdGoodHitMap_%s.png",cut,trigName[kTrigType]));

  // matched hit map
  TH2F *hMthMtdHitMap = (TH2F*)f->Get(Form("hMtdMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMthMtdHitMap,Form("Au+Au di-muon: channel vs backleg of matched MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.MthMtdHitMap_%s.png",cut,trigName[kTrigType]));
}
