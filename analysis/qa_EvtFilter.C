//================================================
void qa_EvtFilter(const int save = 0)
{
  gStyle->SetOptStat(0);

  const char *config = "final";
  const char *track = "Global";
  const char *type[2] = {"Normal","Filtered"};

  TFile *f[2];
  //f[0] = TFile::Open(Form("output/jpsi.%s.%s.final.root",type[0],track),"read");
  //f[1] = TFile::Open(Form("output/jpsi.%s.%s.%s.root",type[1],track,config),"read");

  f[0] = TFile::Open(Form("output/jpsi.Normal.check.root"),"read");
  f[1] = TFile::Open(Form("output/jpsi.Filtered.DbV20150502.root"),"read");

  TH1F *hPtLead[2];

  for(int i=0; i<2; i++)
    {
      printf("==== Work on %s file with %s tracks =====\n",type[i],track);
      
      // trigger mixing
      TH1F *hTrigMix = (TH1F*)f[i]->Get("hTriggerMix");
      hTrigMix->SetName(Form("%s_%s",hTrigMix->GetName(),type[i]));
      for(int ibin=1; ibin<=hTrigMix->GetNbinsX(); ibin++)
	{
	  printf("%20s: %6d\n",hTrigMix->GetXaxis()->GetBinLabel(ibin),hTrigMix->GetBinContent(ibin));
	}

      // event counts
      TH1F *hEventCount = (TH1F*)f[i]->Get("mhEventCount_di_mu");
      hEventCount->SetName(Form("%s_%s",hEventCount->GetName(),type[i]));
      for(int ibin=1; ibin<=hEventCount->GetNbinsX(); ibin++)
	{
	  printf("%20s: %6d\n",hEventCount->GetXaxis()->GetBinLabel(ibin),hEventCount->GetBinContent(ibin));
	}
    }
  
  TList *list = new TList;
  TString legName[2] = {"Normal","Filtered"};

  for(int i=0; i<2; i++)
    {
      hPtLead[i] = (TH1F*)f[i]->Get("mhMthTrkLead_di_mu");
      hPtLead[i]->SetName(Form("%s_%s",hPtLead[i]->GetName(),type[i]));
      list->Add(hPtLead[i]);
    }
  c = drawHistos(list,"MthTrkLeadPt","",kTRUE,0,5,kFALSE,-20,300,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.7,0.6,0.8,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/EventFilter/LeadPt.pdf"));


  list->Clear();
  for(int i=0; i<2; i++)
    {
      TH2F *h2 = (TH2F*)f[i]->Get("mhDzVsPt_di_mu");
      h2->SetName(Form("%s_%s",h2,type[i]));
      TH1F *h1 = (TH1F*)h2->ProjectionY(Form("%s_proy",h2->GetName()));
      list->Add(h1);
    }
  c = drawHistos(list,"Dz","#Deltaz of matched tracks",kFALSE,0,5,kFALSE,-20,300,kFALSE,kTRUE,legName,kTRUE,"",0.2,0.4,0.6,0.8,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/EventFilter/DeltaZ.pdf"));
}
