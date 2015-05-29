TFile *f;
const Int_t trk_index = 0;
const Int_t hlt_index = 0;

const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void qa_muon()
{
  gStyle->SetOptStat(0);

  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  run_cfg_name = run_config;

  qa();
}

//================================================
void qa(const Int_t save = 0)
{
  TH2F *hMuonPhiVsPt = (TH2F*)f->Get(Form("mhMuonPhiVsPt_%s",trigName[kTrigType]));
  c = draw2D(hMuonPhiVsPt);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_phi_vs_pt.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_phi_vs_pt.png",run_type,run_cfg_name.Data()));
    }

  TH2F *hNMuonPosVsNeg = (TH2F*)f->Get(Form("mhNMuonPosVsNeg_%s",trigName[kTrigType]));
  c = draw2D(hNMuonPosVsNeg);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sNMuonPosVsNeg.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sNMuonPosVsNeg.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hNMuonPos = (TH1F*)hNMuonPosVsNeg->ProjectionY("hNMuonPos");
  TH1F *hNMuonNeg = (TH1F*)hNMuonPosVsNeg->ProjectionX("hNMuonNeg");
  TList *list = new TList;
  list->Add(hNMuonPos);
  list->Add(hNMuonNeg);
  TString legName[2];
  legName[0] = Form("N(#mu^{+}) = %4.3f",hNMuonPos->GetMean());
  legName[1] = Form("N(#mu^{-}) = %4.3f",hNMuonNeg->GetMean());
  c = drawHistos(list,"NMuon",Form("Number of muon candidates per event;N"),kFALSE,2.0,3.8,kFALSE,-20,300,kTRUE,kTRUE,legName,kTRUE,"",0.2,0.3,0.3,0.5,kTRUE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sNMuon_pos_neg.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sNMuon_pos_neg.png",run_type,run_cfg_name.Data()));
    }
}
