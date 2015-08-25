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
void qa(const Int_t save = 1)
{
  TList *list = new TList;

  TH2F *hMuonPhiVsPt = (TH2F*)f->Get(Form("mhMuonPhiVsPt_%s",trigName[kTrigType]));
  c = draw2D(hMuonPhiVsPt);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_phi_vs_pt.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_phi_vs_pt.png",run_type,run_cfg_name.Data()));
    }

  TH1F *hMuonPt = (TH1F*)hMuonPhiVsPt->ProjectionX("hMuonPt");
  TH1F *hPosPt = new TH1F("hPosPt","hPosPt",100,0,10);
  TH1F *hNegPt = new TH1F("hNegPt","hNegPt",100,0,10);
  for(int ibin=1; ibin<=hMuonPt->GetNbinsX(); ibin++)
    {
      double center = hMuonPt->GetBinCenter(ibin);
      double content = hMuonPt->GetBinContent(ibin);
      double error = hMuonPt->GetBinError(ibin);
      int bin = hPosPt->FindFixBin(TMath::Abs(center));
      if(center>0)
	{
	  hPosPt->SetBinContent(bin,content);
	  hPosPt->SetBinError(bin,error);
	}
      else
	{
	  hNegPt->SetBinContent(bin,content);
	  hNegPt->SetBinError(bin,error);
	}
    }
  list->Clear();
  list->Add(hPosPt);
  list->Add(hNegPt);
  TString legName2[2] = {"Positive","Negative"};
  c = sysCompare(list,"Muon_pt","p_{T} distribution of muon candidates;p_{T} (GeV/c);counts","Ratio: negative/positive",kTRUE,0,10,kFALSE,0.1,10,kTRUE,0.7,1.3,kTRUE,kTRUE,legName2,kTRUE,"Muon",0.2,0.4,0.2,0.4,kTRUE);
  TLine *line = GetLine(0,1,10,1,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_pt_neg_vs_pos.pdf",run_type,run_cfg_name.Data()));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_muon/%sMuon_pt_neg_vs_pos.png",run_type,run_cfg_name.Data()));
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
  list->Clear();
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
