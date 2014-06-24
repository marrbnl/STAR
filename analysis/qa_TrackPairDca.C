TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

//================================================
void qa_TrackPairDca()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  Dca();
}

//================================================
void Dca(const Int_t save = 0)
{
  TH1F *hTrkPairDcaDr = (TH1F*)f->Get(Form("hTrkPairDcaDr_%s",trigName[kTrigType]));
  hTrkPairDcaDr->SetMarkerStyle(20);
  hTrkPairDcaDr->SetMinimum(0.1);
  c = draw1D(hTrkPairDcaDr,Form("Au+Au %s: DCA_{R} of %s track pair to beam line%s;DCA_{R} (cm)",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  TLine *line = GetLine(5,0,5,hTrkPairDcaDr->GetMaximum()*0.8);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_TrackPairDca/%s.TrkPairDcaDr_%s.png",run_config,trigName[kTrigType]));

  TH1F *hTrkPairDcaDz = (TH1F*)f->Get(Form("hTrkPairDcaDz_%s",trigName[kTrigType]));
  hTrkPairDcaDz->SetMarkerStyle(20);
  //hTrkPairDcaDr->SetMinimum(0.1);
  c = draw1D(hTrkPairDcaDz,Form("Au+Au %s: DCA_{z} of %s track pair to VPD vz %s;DCA_{z} (cm)",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  line = GetLine(4,0,4,hTrkPairDcaDz->GetMaximum()*0.8);
  line->Draw();
  line = GetLine(-4,0,-4,hTrkPairDcaDz->GetMaximum()*0.8);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_TrackPairDca/%s.TrkPairDcaDz_%s.png",run_config,trigName[kTrigType]));
}
