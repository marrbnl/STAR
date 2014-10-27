//================================================
void qa_trigger()
{
  //signal();
  QTcorrelation();
}


//================================================
void signal(const Int_t save = 0)
{
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

  TFile *f = TFile::Open("output/trigger.test.histos.root","read");
  
  TH1F *hNMtdHits = (TH1F*)f->Get("hNMtdHits_di-muon");
  c = draw1D(hNMtdHits,"",kTRUE,kFALSE);
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfMtdHits.png");

  TH1F *hNQTsignals = (TH1F*)f->Get("hNQTsignals_di-muon");
  c = draw1D(hNQTsignals,"",kTRUE,kFALSE);
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfQTsignals.png");

  TH1F *hNMIXsignals = (TH1F*)f->Get("hNMIXsignals_di-muon");
  c = draw1D(hNMIXsignals,"",kTRUE,kFALSE);
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfMT101signals.png");

  TH1F *hNMuons = (TH1F*)f->Get("hNMuons_di-muon");
  c = draw1D(hNMuons,"",kTRUE,kFALSE);
  cout << "2 muons: " << hNMuons->GetBinContent(3)/hNMuons->GetEntries() << endl;
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfTF201signals.png");

  TH1F *hNTrigMtdHits = (TH1F*)f->Get("hNTrigMtdHits_di-muon");
  c = draw1D(hNTrigMtdHits,"",kTRUE,kFALSE);
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfTrigMtdHits.png");

  TH1F *hNTrigGoodMtdHits = (TH1F*)f->Get("hNTrigGoodMtdHits_di-muon");
  c = draw1D(hNTrigGoodMtdHits,"",kTRUE,kFALSE);
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/NumberOfTrigGoodMtdHits.png");
}

//================================================
void QTcorrelation(const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("output/trigger.test.histos.root","read");
  
  TH2F *hMtdHitVsQtSignal = (TH2F*)f->Get("hMtdHitVsQtSignal_di-muon");
  hMtdHitVsQtSignal->SetXTitle("QT signal");
  hMtdHitVsQtSignal->SetYTitle("MTD hit");
  c = draw2D(hMtdHitVsQtSignal,"");
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/MtdHitVsQtSignal.png");
  cout << "MTD hit: " << hMtdHitVsQtSignal->Integral(1,1,2,33)/hMtdHitVsQtSignal->Integral(1,33,2,33) << endl;

  TH2F *hGoodMtdHitVsQtSignal = (TH2F*)f->Get("hGoodMtdHitVsQtSignal_di-muon");
  hGoodMtdHitVsQtSignal->SetXTitle("QT signal");
  hGoodMtdHitVsQtSignal->SetYTitle("MTD hit");
  c = draw2D(hGoodMtdHitVsQtSignal,"");
  if(save) c->SaveAs("~/Work/STAR/analysis/Plots/qa_trigger/GoodMtdHitVsQtSignal.png");
  cout << "MTD hit: " << hGoodMtdHitVsQtSignal->Integral(1,1,2,33)/hGoodMtdHitVsQtSignal->Integral(1,33,2,33) << endl;
}
