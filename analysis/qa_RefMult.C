TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;

//================================================
void qa_RefMult()
{						
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");

  distribution();
}

//================================================
void distribution(const Int_t save = 1)
{
  TH1F *hRefMult = (TH1F*)f->Get(Form("hRefMult_%s",trigName[kTrigType]));
  hRefMult->SetMarkerStyle(20);
  hRefMult->GetXaxis()->SetRangeUser(0,500);
  c = draw1D(hRefMult,Form("Au+Au %s: reference multiplicity w.r.t. default primary vertex%s",trigName[kTrigType],hlt_name[hlt_index]),kTRUE,kTRUE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_RefMult/%s.RefMult_%s.png",run_config,trigName[kTrigType]));
}
