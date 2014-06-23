
//================================================
void qa_Match()
{						
  gStyle->SetOptStat(0);

  //MtdHit();
  //Track();
  DeltaZ();
}

//================================================
void DeltaZ(const Int_t save = 0)
{
  TFile *f = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.NoDzCut.root","read");

  TH1F *hPrimTrkDz[3];
  for(Int_t i=0; i<3; i++)
    {
      TH2F *h2 = (TH2F*)f->Get(Form("hPrimTrkDz_%s",trigName[i]));
      draw2D(h2);
      h2->GetXaxis()->SetRangeUser(1,100);
      hPrimTrkDz[i] = (TH1F*)h2->ProjectionY(Form("hPrimTrkDz_%s_proj",trigName[i]));
      draw1D(hPrimTrkDz[i]);
    }
 
}


//================================================
void Track(const Int_t save = 1)
{
 TFile *f = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.DzCut.root","read");
 TH1F *hStat = (TH1F*)f->Get("mhEventStat");

  // pt distribution
  TH1F *hMthPrimTrkPt[3];
  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      hMthPrimTrkPt[i] = (TH1F*)f->Get(Form("hMthPrimTrkPt_%s",trigName[i]));
      scaleHisto( hMthPrimTrkPt[i], hStat->GetBinContent(i+4), 1, kTRUE);
      list->Add(hMthPrimTrkPt[i]);
    }
  c = drawHistos(list,"hMthPrimTrkPt","Au+Au: p_{T} distribution of matched primary tracks;p_{T} (GeV/c);1/N dN/dp_{T}",kFALSE,0,pi/2,kFALSE,-0.005,0.02,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.68,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/PrimTrkPt.png"));
}

//================================================
void MtdHit(const Int_t save = 1)
{
  TFile *f = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.DzCut.root","read");
  TH1F *hStat = (TH1F*)f->Get("mhEventStat");

  // hit multiplicity
  TH1F *hMthMtdHitN[3];
  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      hMthMtdHitN[i] = (TH1F*)f->Get(Form("hMthMtdHitN_%s",trigName[i]));
      scaleHisto( hMthMtdHitN[i], hStat->GetBinContent(i+4), 1, kTRUE);
      list->Add(hMthMtdHitN[i]);
    }
  c = drawHistos(list,"hMtdHitN","Au+Au: matched MTD hit multiplicity with trigger window cut;N_{hit};1/N_{evt} dN_{hit}",kFALSE,0,pi/2,kFALSE,-0.005,0.02,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.68,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/Mth_MtdGoodHitN.png"));

  // hit map
  TH2F *hMtdMtdHitMap[3];
  for(Int_t i=0; i<3; i++)
    {
      hMtdMtdHitMap[i] = (TH2F*)f->Get(Form("hMtdMtdHitMap_%s",trigName[i]));
      c = draw2D(hMtdMtdHitMap[i]);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_Match/Mth_MtdGoodHitMap_%s.png",trigName[i]));
    }

}
