void forXinjie()
{
  TFile *fin = TFile::Open("Rootfiles/Upsilon/Upsion.pp200.root","read");
  TH1F *hMcUpsilonPt = (TH1F*)fin->Get("Upsilon_pp200_midRap");

  // tracking efficiency & resolution
  TFile *fEff = TFile::Open(Form("Rootfiles/Run14_AuAu200.TrkEff.root",run_type),"read");
  TH2F *hTrkResVsPt = (TH2F*)fEff->Get(Form("PrimTrkRes_vs_TruePt_cent0060"));

  TH1F *hMcTrkPt[2];
  hMcTrkPt[0] = (TH1F*)fEff->Get("hMcTrkPt_MC_cent0060");
  hMcTrkPt[1] = (TH1F*)fEff->Get("hMcTrkPt_MtdTrig_cent0060");

  TH1F *hTrkEff = (TH1F*)hMcTrkPt[1]->Clone("McTrkPtEff_cent0060");
  hTrkEff->Divide(hMcTrkPt[0]);

  TH1F *hTrkEffScale = (TH1F*)hTrkEff->Clone("McTrkPtEff_cent0060_scaled");
  hTrkEffScale->Scale(5);
  
  hTrkEffScale->SetMarkerStyle(21);
  c = draw1D(hTrkEffScale);

  TF1 *func = new TF1("Fit_McTrkPtEff_cent0060_scaled","[0]*exp(-pow([1]/x,[2])-pow([3]/x/x,[4])-pow([5]/x/x/x,[6]))",1,20);
  func->SetParameters(1,0.1,1,0.1,1,0.1,1);
  hTrkEffScale->Fit(func,"R0");
  func->SetLineColor(4);
  func->Draw("sames");

  TFile *fout = TFile::Open("Rootfiles/Xinjie.Run14.AuAu200.UpsilonSmear.root","recreate");
  hMcUpsilonPt->Write("Upsilon_pp200_midRap");
  hTrkResVsPt->Write();
  hTrkEff->Write();
  hTrkEffScale->Write();
  func->Write();
  fout->Close();
}

 
