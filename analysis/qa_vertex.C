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
  if(iPico) run_cfg_name = Form("Pico.%s",run_cfg_name.Data());

  //cuts();
  qa();
}

//================================================
void cuts(const Int_t save = 1)
{
  gStyle->SetOptStat(0);

  // TPC vz cut
  TH1F *hTpcVz = (TH1F*)f->Get(Form("mhTpcVz_%s",trigName[kTrigType]));
  hTpcVz->SetMaximum(1.2*hTpcVz->GetMaximum());
  c = draw1D(hTpcVz,Form("v_{z} of primary vertex reconstructed using TPC"),kFALSE,kFALSE);
  double vz_cut;
  if(year==2014)
    {
      vz_cut = 100;
      TLine *line = GetLine(-1*vz_cut,0,-1*vz_cut,hTpcVz->GetMaximum()*0.8);
      line->Draw();
      TLine *line = GetLine(vz_cut,0,vz_cut,hTpcVz->GetMaximum()*0.8);
      line->Draw();
      TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85);
      t1->AddText(Form("Event fraction: %2.1f%%",hTpcVz->Integral(hTpcVz->FindFixBin(-1*vz_cut+1e-4),hTpcVz->FindFixBin(vz_cut-1e-4))/hTpcVz->GetEntries()*100));
      t1->SetTextFont(62);
      t1->Draw();
    }
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutTpcVz_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutTpcVz_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }

  // TPC vy vs vx cut
  TH2F *hTpcVyVx = (TH2F*)f->Get(Form("mhTpcVyVx_%s",trigName[kTrigType]));
  c = draw2D(hTpcVyVx,Form("v_{y} vs v_{x} of primary vertex reconstructed using TPC"));
  if(save)
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutTpcVyVx_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutTpcVyVx_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
    }


  // TPC-VPD cut
  bool determine_cut = false;
  TH2F *hVzDiffVsTpcVz = (TH2F*)f->Get(Form("mhVzDiffVsTpcVz_%s",trigName[kTrigType]));
  TH1F *hVzDiff = (TH1F*)hVzDiffVsTpcVz->ProjectionY(Form("hVzDiff_%s",trigName[kTrigType]));
  if(determine_cut)
    {
      TF1 *func = new TF1("func_VzDiff","gaus",-1.8,6.5);
      hVzDiff->Fit(func,"IR0");
      func->SetLineColor(2);
    }
  if(year==2014)
    {
      hVzDiff->GetXaxis()->SetRangeUser(-5,5);
    }
  hVzDiff->SetMaximum(10*hVzDiff->GetMaximum());
  c = draw1D(hVzDiff,Form("v_{z} difference of TPC-VPD"),kTRUE,kFALSE);
  double vzdiff_cut;
  if(determine_cut)
    {
      func->Draw("sames");
      double mean = func->GetParameter(1);
      double sigma = func->GetParameter(2);
      double nSigma = 3;
      TLine *line = GetLine(mean-nSigma*sigma,0,mean-nSigma*sigma,hVzDiff->GetMaximum());
      line->Draw();
      TLine *line = GetLine(mean+nSigma*sigma,0,mean+nSigma*sigma,hVzDiff->GetMaximum());
      line->Draw();
      printf("Cut on %2.1f - %2.1f\n",mean-nSigma*sigma,mean+nSigma*sigma);
    }
  if(year==2014)
    {
      vzdiff_cut = 3;
      TLine *line = GetLine(-1*vzdiff_cut,0,-1*vzdiff_cut,hVzDiff->GetMaximum()*0.1);
      line->Draw();
      TLine *line = GetLine(vzdiff_cut,0,vzdiff_cut,hVzDiff->GetMaximum()*0.1);
      line->Draw();
      TPaveText *t1 = GetPaveText(0.2,0.4,0.8,0.85);
      t1->AddText(Form("Event fraction: %2.1f%%",hVzDiff->Integral(hVzDiff->FindFixBin(-1*vzdiff_cut+1e-4),hVzDiff->FindFixBin(vzdiff_cut-1e-4))/hVzDiff->GetEntries()*100));
      t1->SetTextFont(62);
      t1->Draw();
    }

  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutDiffVz_%s.pdf",run_type,run_cfg_name.Data(),trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%sCutDiffVz_%s.png",run_type,run_cfg_name.Data(),trigName[kTrigType]));
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
	    }
	  if(hName[i]=="mhVpdVz")
	    logy = kFALSE;
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%s%s_%s.pdf",run_type,run_cfg_name.Data(),outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%s%s_%s.png",run_type,run_cfg_name.Data(),outname.Data(),trigName[kTrigType]));
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
