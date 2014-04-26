
//================================================
void hlt_rate(const char *day = "077", const Int_t save = 0)
{
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open(Form("Output/output.%s.root",day),"read");
  TH1F *hEventRate = (TH1F*)f->Get("hEventRate_di-muon");
  hEventRate->SetMaximum(10*hEventRate->GetMaximum());
  c = draw1D(hEventRate,"Event counts for different requirements (di-muon trigger);;Counts",kTRUE,kFALSE);
  hEventRate->SetMarkerSize(1.5);
  hEventRate->Draw("sames TEXT45");
  TPaveText *t1 = GetPaveText(0.3,0.7,0.5,0.88,0.03);
  t1->SetTextAlign(13);
  t1->AddText("1st bin: All events");
  t1->AddText("2nd bin: At least two p_{T,mth} > 1 GeV/c for all MTD hits");
  t1->AddText("3rd bin: At least two p_{T,mth} > 1 GeV/c for good MTD hits");
  t1->AddText("4th bin: At least two p_{T,mth} > 1 GeV/c && one p_{T,mth} > 1.5 GeV/c");
  t1->AddText("5th bin: At least two p_{T,mth} > 1 GeV/c && |#Deltaz|<50cm");
  t1->AddText("6th bin: At least two p_{T,mth} > 1 GeV/c && |#Deltaz|<30cm");
  t1->Draw();
  t1 = GetPaveText(0.2,0.3,0.15,0.24,0.035);
  t1->AddText("Au+Au");
  t1->AddText(" #sqrt{s_{NN}} = 200 GeV");
  t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/QA/Plots/HLT/%s_MtdEventRateForHLT.png",day));
    }
}
