TFile *f;

const char *run_config = "";
const Bool_t iPico = 1;
const int year = 2014;
TString run_cfg_name;

//================================================
void qa_vertex()
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  TString fileName;

  if(year==2013)
    {
      run_type = "Run13_pp500";
      if(iPico) fileName = Form("Pico.Run13.pp500.jpsi.%sroot",run_config);
      else      fileName = Form("Run13.pp500.jpsi.%sroot",run_config);
    }
  else if(year==2014)
    {
      run_type = "Run14_AuAu200";
      if(iPico) fileName = Form("Pico.Run14.AuAu200.jpsi.%sroot",run_config);
      else      fileName = Form("Run14.AuAu200.jpsi.%sroot",run_config);
    }

  f = TFile::Open(Form("./output/%s",fileName.Data()),"read");

  run_cfg_name = run_config;

  //cuts();
  //qa();
}

//================================================
void cuts(const Int_t save = 0)
{
  gStyle->SetOptStat(0);

  // TPC vz cut
  TH1F *hTpcVz = (TH1F*)f->Get(Form("mhTpcVz_%s",trigName[kTrigType]));
  //TF1 *func = new TF1("func_TpcVz","gaus",-30,35);
  //hTpcVz->Fit(func,"IR0");
  //func->SetLineColor(2);
  hTpcVz->SetMaximum(1.2*hTpcVz->GetMaximum());
  c = draw1D(hTpcVz,Form("v_{z} of primary vertex reconstructed using TPC"),kFALSE,kFALSE);
  //TPaveText *t1 = GetJpsiPaveText(0.13,0.5,0.75,0.88, run_type, 0.035);
  //t1->Draw();
  //func->Draw("sames");
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVz_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVz_%s.png",run_type,trigName[kTrigType]));
    }

  // TPC vy vs vx cut
  TH2F *hTpcVyVx = (TH2F*)f->Get(Form("mhTpcVyVx_%s",trigName[kTrigType]));
  //hTpcVyVx->GetXaxis()->SetRangeUser(-0.5,0.5);
  //hTpcVyVx->GetYaxis()->SetRangeUser(-1,0);
  c = draw2D(hTpcVyVx,Form("v_{y} vs v_{x} of primary vertex reconstructed using TPC"));
  //t1->Draw();
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVyVx_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutTpcVyVx_%s.png",run_type,trigName[kTrigType]));
    }


  // TPC-VPD cut
  TH2F *hVzDiffVsTpcVz = (TH2F*)f->Get(Form("mhVzDiffVsTpcVz_%s",trigName[kTrigType]));
  TH1F *hVzDiff = (TH1F*)hVzDiffVsTpcVz->ProjectionY(Form("hVzDiff_%s",trigName[kTrigType]));
  //func = new TF1("func_VzDiff","gaus",-1.8,6.5);
  //hVzDiff->Fit(func,"IR0");
  //func->SetLineColor(2);
  //hVzDiff->GetXaxis()->SetRangeUser(-30,50);
  hVzDiff->SetMaximum(20*hVzDiff->GetMaximum());
  c = draw1D(hVzDiff,Form("v_{z} difference of TPC-VPD"),kTRUE,kFALSE);
  //t1->Draw();
  //func->Draw("sames");
  TPaveText *t1 = GetPaveText(0.15,0.6,0.8,0.85);
  t1->AddText(Form("Fraction of events within [-3,3] cm: %2.1f%%",hVzDiff->Integral(hVzDiff->FindFixBin(-3),hVzDiff->FindFixBin(3))/hVzDiff->GetEntries()*100));
  t1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutDiffVz_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CutDiffVz_%s.png",run_type,trigName[kTrigType]));
    }
}


//================================================
void qa(const Int_t save = 0)
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

//================================================
void ranking(const Int_t save = 1)
{
  gStyle->SetOptStat(0);

  TH2F *hNVtxVsGoodRank = (TH2F*)f->Get(Form("mhNVtxVsGoodRank_%s",trigName[kTrigType]));
  hNVtxVsGoodRank->GetXaxis()->SetRangeUser(0,15);
  hNVtxVsGoodRank->GetYaxis()->SetRangeUser(0,20);
  c = draw2D(hNVtxVsGoodRank);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/NVtxVsGoodRank_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/NVtxVsGoodRank_%s.png",run_type,trigName[kTrigType]));
    }

  TH2F *hNVtxVsNTrk = (TH2F*)f->Get(Form("mhNVtxVsNTrk_%s",trigName[kTrigType]));
  hNVtxVsNTrk->GetYaxis()->SetRangeUser(0,20);
  c = draw2D(hNVtxVsNTrk);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/NVtxVsNTrkUsed_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/NVtxVsNTrkUsed_%s.png",run_type,trigName[kTrigType]));
    }

  TH2F *hNTrkUsedVsRank = (TH2F*)f->Get(Form("mhNTrkUsedVsRank_%s",trigName[kTrigType]));
  draw2D(hNTrkUsedVsRank);
  TH1F *hNTrkUsed[2];
  TList *list = new TList;
  for(Int_t i=0; i<2; i++)
    {
      hNTrkUsed[i] = (TH1F*)hNTrkUsedVsRank->ProjectionY(Form("NTrkUsed_%d",i),i+1,i+1);
      hNTrkUsed[i]->SetMaximum(1.4*hNTrkUsed[i]->GetMaximum());
      list->Add(hNTrkUsed[i]);
    }
  TString legName[2] = {"Ranking > 0","Ranking < 0"};
  c = drawHistos(list,"hNTrkUsed",Form("%s: number of tracks used to reconstruct vertex;# of used tracks",trigName[kTrigType]),kFALSE,0,100,kFALSE,-3,8,kFALSE,kTRUE,legName,kTRUE,"",0.5,0.65,0.6,0.8,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CompareNTrkUsed_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/CompareNTrkUsed_%s.png",run_type,trigName[kTrigType]));
    }
}
