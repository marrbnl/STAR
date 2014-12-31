TFile *f;
Int_t vtx_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx";

//================================================
void qa_vertex(const Int_t save = 0)
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  f = TFile::Open("./output/Run13.pp500.jpsi.EventQA.root","read");

  TString cut_name = run_config;

  if(cut_name.Contains("ClosePrimVtx"))
    vtx_index = 1;
  else if (cut_name.Contains("MtdVtx"))
    vtx_index = 2;


  cuts();
  //qa();
}

//================================================
void cuts(const Int_t save = 0)
{
  gStyle->SetOptStat(0);

  // TPC vz cut
  TH1F *hTpcVz = (TH1F*)f->Get(Form("mhTpcVz_%s",trigName[kTrigType]));
  TF1 *func = new TF1("func_TpcVz","gaus",-30,35);
  hTpcVz->Fit(func,"IR0");
  func->SetLineColor(2);
  hTpcVz->SetMaximum(1.2*hTpcVz->GetMaximum());
  c = draw1D(hTpcVz,Form("v_{z} of primary vertex reconstructed using TPC"),kFALSE,kFALSE);
  TPaveText *t1 = GetJpsiPaveText(0.13,0.2,0.78,0.88, run_type, 0.035);
  t1->Draw();
  func->Draw("sames");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVz_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVz_%s.png",run_type,trigName[kTrigType]));
    }

  // TPC vy vs vx cut
  TH2F *hTpcVyVx = (TH2F*)f->Get(Form("mhTpcVyVx_%s",trigName[kTrigType]));
  hTpcVyVx->GetXaxis()->SetRangeUser(-0.5,0.5);
  hTpcVyVx->GetYaxis()->SetRangeUser(-1,0);
  c = draw2D(hTpcVyVx,Form("v_{y} vs v_{x} of primary vertex reconstructed using TPC"));
  t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVyVx_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVyVx_%s.png",run_type,trigName[kTrigType]));
    }


  // TPC-VPD cut
  //TH2F *hVzDiffVsTpcVz = (TH2F*)f->Get(Form("mhVzDiffVsTpcVz_%s",trigName[kTrigType]));
  //TH1F *hVzDiff = (TH1F*)hVzDiffVsTpcVz->ProjectionX(Form("hVzDiff_%s",trigName[kTrigType]));
  //hVzDiff->Rebin(2);
  TH1F *hVzDiff = (TH1F*)f->Get(Form("mhDiffVzWithCut_%s",trigName[kTrigType]));
  func = new TF1("func_VzDiff","gaus",-1.8,6.5);
  hVzDiff->Fit(func,"IR0");
  func->SetLineColor(2);
  hVzDiff->GetXaxis()->SetRangeUser(-30,50);
  hVzDiff->SetMaximum(1.2*hVzDiff->GetMaximum());
  c = draw1D(hVzDiff,Form("v_{z} difference of TPC-VPD"),kFALSE,kFALSE);
  t1->Draw();
  func->Draw("sames");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutDiffVz_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutDiffVz_%s.png",run_type,trigName[kTrigType]));
    }
}


//================================================
void qa(const Int_t save = 1)
{
  const char *hName[] = {"mhNPrimVtx","mhTpcVx","mhTpcVy","mhTpcVr",
			 "mhTpcVxVz","mhTpcVyVz","mhVpdVz","mhVzDiffVsTpcVz"};

  for(Int_t i=0; i<8; i++)
    {
      TH1 *h = (TH1*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      if(!h->IsA()->InheritsFrom("TH2")) 
	{
	  Bool_t logy = kTRUE;
	  if(run_type == "Run13_pp500")
	    {
	      if(hName[i]=="mhNPrimVtx")
		h->GetXaxis()->SetRangeUser(0,20);
	      if(hName[i]=="mhVpdVz")
		logy = kFALSE;
	    }
	  c = draw1D(h,"",logy,kFALSE);
	}
      else
	{
	  TH2F *h2 = (TH2F*)h;
	  c = draw2D(h2,"");
	}

      TString outname = hName[i];
      outname.ReplaceAll("mh","");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
}
