TFile *f;

//================================================
void qa_Event()
{
  gStyle->SetOptStat(0);
  f = TFile::Open("./output/Run13.pp500.jpsi.EventQA.root","read");
 
  ana();
}

//================================================
void ana(const Int_t save = 1)
{
  // statistics
 TH1F *hTriggerMix = (TH1F*)f->Get(Form("hTriggerMix"));
  c = draw1D(hTriggerMix,"",kTRUE,kFALSE);
  printStat(hTriggerMix);
  
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Event/TriggerMix.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Event/TriggerMix.png",run_type));
    }

  TH1F *hEventStat = (TH1F*)f->Get(Form("hEventStat"));
  c = draw1D(hEventStat,"",kTRUE,kFALSE);
  printStat(hEventStat);
  
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Event/EventStat.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_Event/EventStat.png",run_type));
    }
}


//================================================
void printStat(const TH1 *h)
{
  if(!h) return;
  printf("=== print stats for %s ===\n",h->GetName());
  for(Int_t ibin=1; ibin<=h->GetNbinsX(); ibin++)
    {
      printf("%15s: %5.4e\n",h->GetXaxis()->GetBinLabel(ibin),h->GetBinContent(ibin));
    }
}
