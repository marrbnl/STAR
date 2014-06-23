void check()
{
  merge();
}

//================================================
void merge()
{
  TFile *f[20];
  f[0] = TFile::Open("Output/077.root");
  for(Int_t i=1; i<20; i++)
    f[i] = TFile::Open(Form("Output/jpsi.histos.077.%02d.root",i-1));

  TH1F *hStat[20];
  for(Int_t i=0; i<20; i++)
    {
      hStat[i] = (TH1F*)f[i]->Get("mhEventStat");
      hStat[i]->SetLineColor(i+1);
      if(i==0) draw1D(hStat[i],"",kFALSE,kFALSE);
      else     hStat[i]->Draw("HIST sames");
    }
}
