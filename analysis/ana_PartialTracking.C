TFile *f;

void ana_PartialTracking()
{
  gStyle->SetOptStat(0);
  f = TFile::Open("output/PartialTracking.histos.root","read");
  
  //tracks();
  //tracking();
  hits();
  //efficiency();
}

//================================================
void efficiency(const Int_t save = 0)
{
  TH1F *hEventCount = (TH1F*)f->Get("hEventCount");
  TH1F *hEventCountDaq10k = (TH1F*)f->Get("hEventCountDaq10k");
  for(Int_t ibin=1; ibin<=10; ibin++)
    {
      printf("%s: %d/%d = %1.2f%%\n",hEventCount->GetXaxis()->GetBinLabel(ibin),hEventCountDaq10k->GetBinContent(ibin),hEventCount->GetBinContent(ibin),hEventCountDaq10k->GetBinContent(ibin)*100./hEventCount->GetBinContent(ibin));
    }
}

//================================================
void hits(const Int_t save = 1)
{
  TH1F *hDeltaSector = (TH1F*)f->Get("hDeltaSector_di-muon");
  printf("Total number of matched hits: %d\n",hDeltaSector->GetEntries());
  printf("Found: %d = %1.2f%%\n", hDeltaSector->GetBinContent(1),hDeltaSector->GetBinContent(1)/hDeltaSector->GetEntries()*100);
  printf("Neighbor: %d = %1.2f%%\n", hDeltaSector->GetBinContent(2),hDeltaSector->GetBinContent(2)/hDeltaSector->GetEntries()*100);
  printf("Opposite eta: %d = %1.2f%%\n", hDeltaSector->GetBinContent(13),hDeltaSector->GetBinContent(13)/hDeltaSector->GetEntries()*100);
  printf("Opposite neighbor: %d = %1.2f%%\n", hDeltaSector->GetBinContent(14),hDeltaSector->GetBinContent(14)/hDeltaSector->GetEntries()*100);

  TH1F *hHitModuleOpposEta = (TH1F*)f->Get("hHitModuleOpposEta_di-muon");
  c = draw1D(hHitModuleOpposEta,"Module of MTD hits with matched tracks in opposite #eta",kFALSE,kFALSE,0.04,"TEXT45 HIST");
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_PartialTracking/HitModule_OpposEta.png"));
  cout << "Probility in Module 3 is: " << hHitModuleOpposEta->Integral(3,3)/hHitModuleOpposEta->Integral(1,5) << endl;

  TH2F *hMtdTrigHitMap = (TH2F*)f->Get("hMtdTrigHitMap_di-muon");
  c = draw2D(hMtdTrigHitMap);
  TH1F *hMtdModule = (TH1F*)hMtdTrigHitMap->ProjectionY("hMtdModule");
  hMtdModule->Rebin(12);
  c = draw1D(hMtdModule,"Multiplicity of MTD hits fired the trigger;Module",kFALSE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_PartialTracking/HitModule_Trigger.png"));
  cout << "Probility in Module 3 is: " << hMtdModule->Integral(3,3)/hMtdModule->Integral(1,5) << endl;

}

//================================================
void tracking(const Int_t save = 0)
{
  gStyle->SetOptStat(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);

 TH1F *hNFiredTpcSectors = (TH1F*)f->Get("hNFiredTpcSectors_di-muon");
  c = draw1D(hNFiredTpcSectors,"# of TPC sectors for partial tracking",kTRUE,kFALSE,0.05);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_PartialTracking/NFiredTpcSectors.png"));
}

//================================================
void tracks(const Int_t save = 0)
{
  TH1F *hNSectorsInTrk = (TH1F*)f->Get("hNSectorsInTrk_di-muon");
  c = draw1D(hNSectorsInTrk,"# of sectors contributing to matched global tracks",kTRUE,kFALSE,0.05);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_PartialTracking/NSectors_mth_gTrk.png"));

  TH2F *hMaxHitFractionVsPt = (TH2F*)f->Get("hMaxHitFractionVsPt_di-muon");
  hMaxHitFractionVsPt->GetYaxis()->SetRangeUser(0.4,1.1);
  printf("%1.2f%% of tracks have more than 90%% of hits in one sector\n",hMaxHitFractionVsPt->Integral(0,-1,hMaxHitFractionVsPt->GetYaxis()->FindFixBin(0.9),101)/hMaxHitFractionVsPt->Integral(1,101)*100);
  c = draw2D(hMaxHitFractionVsPt,"Largest fraction of hits in one sector vs p_{T};p_{T} (GeV/c);fraction",0.05);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/ana_PartialTracking/MaxHitFractionVsPt_mth_gTrk.png"));
}

