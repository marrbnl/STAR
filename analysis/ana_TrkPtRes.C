TFile *f;
const char *signal_name[2] = {"J/#psi","#Upsilon(1S)"};
Int_t trk_index = 0;
Int_t signal_index = 0;
const char *run_config = "Embed.Upsilon.PrimaryTrk";

//================================================
void ana_TrkPtRes()
{
  gStyle->SetOptStat(0);

  TString cut_name = run_config;
  if(cut_name.Contains("Global"))
    trk_index = 1;
  if(cut_name.Contains("Upsilon"))
    signal_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  printf("# of di-muon events: %d\n",hStat->GetBinContent(7));

  embedding();
}

//================================================
void embedding(const Int_t save = 1)
{
  TH2F *hpTrkRes = (TH2F*)f->Get("hpTrkPtRes_di-muon");
  draw2D(hpTrkRes);
  TObjArray pSlices;
  hpTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &pSlices);
  TH1D *hpRes = (TH1D*)pSlices[2];
  hpRes->SetYTitle("#sigma_{p_{T}}/p_{T}");
  hpRes->SetMarkerStyle(21);
  c = draw1D(hpRes,Form("Au+Au %s embedding: p_{T} resolution of primary muon tracks",trigName[kTrigType]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_TrkPtRes/%s.MC_primary_muon_ptRes.png",run_config));

  TH2F *hgTrkRes = (TH2F*)f->Get("hgTrkPtRes_di-muon");
  draw2D(hgTrkRes);
  TObjArray gSlices;
  hgTrkRes->FitSlicesY(0, 0, -1, 0, "QNR", &gSlices);
  TH1D *hgRes = (TH1D*)gSlices[2];
  hgRes->SetYTitle("#sigma_{p_{T}}/p_{T}");
  hgRes->SetMarkerStyle(21);
  c = draw1D(hgRes,Form("Au+Au %s embedding: p_{T} resolution of global muon tracks",trigName[kTrigType]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_TrkPtRes/%s.MC_global_muon_ptRes.png",run_config));

  // Compare with other measurement
  TH1F *hpMuon = (TH1F*)hpRes->Clone("hpTrkPtRes");
  hpMuon->GetYaxis()->SetRangeUser(0,0.28);
  c = draw1D(hpMuon,"p_{T} resolution of single tracks;p_{T} (GeV/c);#sigma_{p_{T}}/p_{T}   ");

  TH1F *hgMuon = (TH1F*)hgRes->Clone("hgTrkPtRes");
  hgMuon->SetMarkerStyle(20);
  hgMuon->SetMarkerColor(6);
  hgMuon->SetLineColor(6);
  hgMuon->Draw("sames");

  TFile *f1 = TFile::Open("Rootfiles/ElectronPtRes_Run13_Embedding.root");
  TH1F *hElec = (TH1F*)f1->Get("hJpsiElectronPtSigma");
  hElec->SetMarkerStyle(25);
  hElec->SetMarkerColor(2);
  hElec->SetLineColor(2);
  hElec->Draw("same");

  TFile *f2 = TFile::Open("Rootfiles/TrkPtRes_Run14_Cosmic.root");
  TH1F *hCosmic[3];
  hCosmic[0] = (TH1F*)f2->Get("h_ybin_80_2");
  hCosmic[1] = (TH1F*)f2->Get("h_ybin_100_2");
  hCosmic[2] = (TH1F*)f2->Get("h_ybin_120_2");
  for(Int_t i=0; i<1; i++)
    {
      hCosmic[i]->Scale(1./TMath::Sqrt(2));
      hCosmic[i]->SetStats(kFALSE);
      hCosmic[i]->SetMarkerStyle(29);
      hCosmic[i]->SetMarkerSize(1.2);
      hCosmic[i]->SetMarkerColor(4+2*i);
      hCosmic[i]->SetLineColor(4+2*i);
      hCosmic[i]->Draw("same");
    }

  TLegend *leg = new TLegend(0.12,0.65,0.5,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hpMuon,"Embed: Run14 Au+Au 200 GeV, #Upsilon(1S)#rightarrow#mu#mu, primary");
  leg->AddEntry(hgMuon,"Embed: Run14 Au+Au 200 GeV, #Upsilon(1S)#rightarrow#mu#mu, global");
  leg->AddEntry(hElec,"Embed: Run13 p+p 500 GeV, J/#psi#rightarrowee, primary");
  leg->AddEntry(hCosmic[0],"Cosmic ray: Run14 NbinsY=80");
  //leg->AddEntry(hCosmic[1],"Cosmic ray: Run14 NbinsY=100");
  //leg->AddEntry(hCosmic[2],"Cosmic ray: Run14 NbinsY=120");
  leg->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_TrkPtRes/compare_ptRes.png"));
}
