TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void qa_DCA()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  
  gDCA();
}

//================================================
void gDCA(const Int_t save = 1)
{
  TH2F *hTrkDca = (TH2F*)f->Get(Form("hTrkDca_qa_%s",trigName[kTrigType]));
  c = draw2D(hTrkDca,Form("Au+Au %s: dca vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_DCA/%s.%s_gDCA_vs_pt_%s.png",run_config,trk_name[trk_index],trigName[kTrigType]));
}
