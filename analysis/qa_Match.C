TFile *f;
const char *cut = "GlobalTrk.HLT";
Int_t hlt_index = 0;
Int_t trk_index = 0;

//================================================
void qa_Match()
{						
  gStyle->SetOptStat(0);

  TString cut_name = cut;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",cut),"read");

  Track();
  //DeltaZ();
}

//================================================
void DeltaZ(const Int_t save = 0)
{
  TH2F *hTrkDzVsPt = (TH2F*)f->Get(Form("hTrkDz_%s",trigName[kTrigType]));
  hTrkDzVsPt->GetXaxis()->SetRangeUser(0,6);
  hTrkDzVsPt->GetYaxis()->SetRangeUser(-100,100);
  draw2D(hTrkDzVsPt,Form("Au+Au %s: #Deltaz of matched %s track-hit pairs%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));

  hTrkDzVsPt->GetXaxis()->SetRangeUser(1,100);
  TH1F *hMthDz = (TH1F*)hTrkDzVsPt->ProjectionY(Form("hTrkDzVsPt_%s_proj",trigName[kTrigType]));
  draw1D(hMthDz);
 
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
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_MtdHit/%s.NMthTrk_%s.png",cut,trigName[kTrigType]));
}
