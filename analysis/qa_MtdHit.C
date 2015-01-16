TFile *f;

//================================================
void qa_MtdHit(const Int_t save = 0)
{
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  


  //f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  f = TFile::Open("./output/Run13.pp500.jpsi.EventQA.root","read");

  if(!save)
    {
      qa();
    }
  else
    {
      qa(1);
    }

}


//================================================
void qa(const Int_t save = 1)
{
  const char *hName[] = {"mhMtdRawHitN","mhMtdHitN","mhMthMtdHitN","mhMtdHitMap",
			 "mhMtdHitTrigTime","mhMtdRawHitMap","mhMthMtdHitMap","mhMthMtdHitLocaly","mhMthMtdHitGlobalz"};

  for(Int_t i=0; i<9; i++)
    {
      TH1 *h = (TH1*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      if(!h->IsA()->InheritsFrom("TH2")) 
	{
	  c = draw1D(h,"",kTRUE,kFALSE);
	  gStyle->SetOptStat(1);
	}
      else
	{
	  TH2F *h2 = (TH2F*)h;
	  c = draw2D(h2,"");
	  gStyle->SetOptStat(0);
	}

      TString outname = hName[i];
      outname.ReplaceAll("mh","");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_MtdHit/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_MtdHit/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }

  TH2F *h2 = (TH2F*)f->Get(Form("mhMtdHitMap_%s",trigName[kTrigType]));

  TH2F *hHitMapRot = new TH2F("hHitMapRot","di_muo: #varphi vs #eta distribution of MTD hits;#eta;#varphi",
			      h2->GetYaxis()->GetNbins(),-0.6,0.6,
			      h2->GetXaxis()->GetNbins(),0,2*pi);

  for(Int_t ibin=1; ibin<=h2->GetNbinsX(); ibin++)
    {
      for(Int_t jbin=1; jbin<=h2->GetNbinsY(); jbin++)
	{
	  Int_t iibin = jbin;
	  Int_t jjbin = ibin + 7;
	  if(ibin>23) jjbin = ibin-23;
	  hHitMapRot->SetBinContent(iibin,jjbin,h2->GetBinContent(ibin,jbin));
	  hHitMapRot->SetBinError(iibin,jjbin,h2->GetBinError(ibin,jbin));
	}
    }
  c = draw2D(hHitMapRot);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_MtdHit/MtdHitMap_Rotate_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_MtdHit/MtdHitMap_Rotate_%s.png",run_type,trigName[kTrigType]));
    }
}
