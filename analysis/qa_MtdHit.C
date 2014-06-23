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
  TH1F *hMtdHitN = (TH1F*)f->Get(Form("hMtdHitN_%s",trigName[kTrigType]));
  scaleHisto( hMtdHitN, hStat->GetBinContent(kTrigType+7), 1);
  hMtdHitN->SetMarkerStyle(21);
  c = draw1D(hMtdHitN,Form("Au+Au di-muon: multiplicity of good MTD hits%s;N_{hit};1/N_{evt} dN_{hit}",hlt_name[hlt_index]),kTRUE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.MtdGoodHitN_%s.png",cut,trigName[kTrigType]));

  // hit map
  TH2F *hMtdHitMap = (TH2F*)f->Get(Form("hMtdHitMap_%s",trigName[kTrigType]));
  c = draw2D(hMtdHitMap,Form("Au+Au di-muon: channel vs backleg of good MTD hits%s",hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.MtdGoodHitMap_%s.png",cut,trigName[kTrigType]));
}
