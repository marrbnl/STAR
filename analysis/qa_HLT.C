TFile *f[2];
Int_t trk_index = -1;
const char *name[2] = {"Regular","HLT"};

//================================================
void qa_HLT()
{
  gStyle->SetOptStat(0);

  f[0] = TFile::Open(Form("Output/jpsi.AuAu200.Run14.GlobalTrk.All.root"),"read");
  f[1] = TFile::Open(Form("Output/jpsi.AuAu200.Run14.GlobalTrk.HLT.root"),"read");
  //f[0] = TFile::Open(Form("Output/jpsi.AuAu200.Run14.GlobalTrk.NoDCA.TestSample.root"),"read");
  //f[1] = TFile::Open(Form("Output/jpsi.AuAu200.Run14.GlobalTrk.NoDCA.HLT.root"),"read");
  trk_index = 1;

  match();
}


//================================================
void match(const Int_t save = 1)
{
  TH1F *hMthNhit[2];
  TList *list = new TList;
  TString legName[2];
  for(Int_t i=0; i<2; i++)
    {
      hMthNhit[i] = (TH1F*)f[i]->Get(Form("hMthMtdHitN_%s",trigName[kTrigType]));
      hMthNhit[i]->SetName(Form("%s_%s",name[i],hMthNhit[i]->GetName()));
      hMthNhit[i]->Sumw2();
      hMthNhit[i]->Scale(1./hMthNhit[i]->Integral());
      list->Add(hMthNhit[i]);
      legName[i] = Form("%s events: <N> = %2.2f",name[i],hMthNhit[i]->GetMean());
    }

  c = drawHistos(list,"NMthMtdHit",Form("Au+Au %s: # of MTD hits matched to %s tracks;N_{hit}",trigName[kTrigType],trk_name[trk_index]),kFALSE,0,100,kTRUE,1e-7,10,kTRUE,kTRUE,legName,kTRUE,"di-muon",0.5,0.7,0.65,0.85);

  if(save)
    {
      TString file_name = f[0]->GetName();
      if(file_name.Contains("NoDCA"))
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_HLT/compre_NMthHit_%s_NoDCA_%s.png",trk_name[trk_index],trigName[kTrigType]));
      else
	c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_HLT/compre_NMthHit_%s_%s.png",trk_name[trk_index],trigName[kTrigType]));
    }
}
