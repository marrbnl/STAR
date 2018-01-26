TFile *f;
const int year = YEAR;

//================================================
void qa_vertex()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  f = TFile::Open(Form("./output/%s.jpsi.root",run_type),"read");

  qa();
  //vtxCuts();
  //cuts();
  //vtxSelection();
}

//================================================
void qa(const Int_t save = 1, const int saveAN = 1)
{
  const char *hName[] = {"mhTpcVyVx_di_mu","mhTpcVyVxWithCut_di_mu","mhTpcVz_di_mu","mhTpcVzWithCut_di_mu",
			 "mhVzDiffVsTpcVz_di_mu","mhDiffVzWithCut_di_mu"};
  const char *hTitle[] = {"y vs. x of primary vertex (w/o vertex cut)","y vs. x of primary vertex (w/ vertex cut)",
			  "z distribution of primary vertex (w/o vertex cut)", "z distribution of primary vertex (w/ vertex cut)",
			  "correlation of vz between VPD and TPC (w/o vertex cut)", "#Deltavz distribution of primary vertex (w/ vertex cut)"};
  const char *xTitle[] = {"vx (cm)", "vx (cm)", "vz_{TPC} (cm)", "vz_{TPC} (cm)", "vz_{TPC} (cm)", "vz_{TPC}-vz_{VPD} (cm)"};
  const char *yTitle[] = {"vy (cm)", "vy (cm)", "", "", "vz_{TPC}-vz_{VPD} (cm)", ""};

  for(Int_t i=0; i<6; i++)
    {
      TH1 *h = (TH1*)f->Get(Form("%s",hName[i]));
      h->SetTitle(Form("%s: %s;%s;%s",run_type,hTitle[i],xTitle[i],yTitle[i]));
      if(!h->IsA()->InheritsFrom("TH2")) 
	{
	  if(i==5) h->GetXaxis()->SetRangeUser(-5,5);
	  Bool_t logy = kTRUE;
	  if(hName[i]=="mhVpdVz")
	    logy = kFALSE;
	  c = draw1D(h,"",logy,kFALSE);
	}
      else
	{
	  TH2F *h2 = (TH2F*)h;
	  if(i==0 || i==1 || i==4)
	    {
	      h2->RebinX(2);
	      h2->RebinY(2);
	    }
	  c = draw2D(h2,"");
	}

      TString outname = hName[i];
      outname.ReplaceAll("mh","");
      outname.ReplaceAll("_di_mu","");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/%s.pdf",run_type,outname.Data()));
	}
      if(saveAN)
	{
	  c->SaveAs(Form("~/Dropbox/STAR\ Quarkonium/Run14_Jpsi/Analysis\ note/Figures/Chp1_%s.pdf",outname.Data()));
	}
    }
}

//================================================
void vtxSelection(const Int_t save = 1)
{
  gStyle->SetOptStat(0);

  // TPC vz
  TH1F *hTpcVz = (TH1F*)f->Get("hVertexZ");
  hTpcVz->SetMaximum(1.2*hTpcVz->GetMaximum());
  c = draw1D(hTpcVz,Form("v_{z} of primary vertex reconstructed using TPC"),kFALSE,kFALSE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_TpcVz.pdf",run_type));

  // TPC vy vs vx
  TH2F *hTpcVyVx = (TH2F*)f->Get("hVertexXY");
  c = draw2D(hTpcVyVx,Form("v_{y} vs v_{x} of primary vertex reconstructed using TPC"));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_TpcVyVx.pdf",run_type));

  TH2F *hTpcVxVz = (TH2F*)f->Get("hVertexXZ");
  c = draw2D(hTpcVxVz,Form("v_{x} vs v_{z} of primary vertex reconstructed using TPC"));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_TpcVxVz.pdf",run_type));

  TH2F *hTpcVyVz = (TH2F*)f->Get("hVertexYZ");
  c = draw2D(hTpcVyVz,Form("v_{y} vs v_{z} of primary vertex reconstructed using TPC"));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_TpcVyVz.pdf",run_type));

  // TPC-VPD
  TH2F *hVtxZvsVpdVz[2];
  TH1F *hVtxZDiff[2];
  hVtxZvsVpdVz[0] = (TH2F*)f->Get("hVtxZvsVpdVzDefault");
  hVtxZvsVpdVz[1] = (TH2F*)f->Get("hVtxZvsVpdVzClosest");
  hVtxZDiff[0] = (TH1F*)f->Get("hVtxZDiffDefault");
  hVtxZDiff[1] = (TH1F*)f->Get("hVtxZDiffClosest");

  const char *name[2] = {"default","closest"};
  for(int i=0; i<2; i++)
    {
      c = draw2D(hVtxZvsVpdVz[i],"");
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_TpcVpdVz_%s.pdf",run_type,name[i]));

      hVtxZDiff[i]->Rebin(4);
      hVtxZDiff[i]->SetMarkerStyle(21);
      hVtxZDiff[i]->GetXaxis()->SetRangeUser(-15,20);
      TF1 *func = new TF1(Form("func_VzDiff_%d",i),"gaus",-4,3);
      hVtxZDiff[i]->Fit(func,"IR0");
      func->SetLineColor(2);
      c = draw1D(hVtxZDiff[i],Form("TPC vz - VPD vz (%s)",name[i]));
      func->Draw("sames");
      TPaveText *t1 = GetPaveText(0.2,0.3,0.7,0.8,0.04);
      t1->AddText(Form("|#Deltaz| < 5 cm: %4.2f%%",hVtxZDiff[i]->Integral(hVtxZDiff[i]->FindBin(-5),hVtxZDiff[i]->FindBin(5))/hVtxZDiff[i]->GetEntries()*100));
      t1->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_FitDiffVz_%s.pdf",run_type,name[i]));
    }

  // vertex selection
  TH2F *hVtxIndex = (TH2F*)f->Get("hVtxIndClosestVsRank");
  hVtxIndex->GetXaxis()->SetRangeUser(0,12);
  hVtxIndex->GetYaxis()->SetRangeUser(0,12);
  c = draw2D(hVtxIndex,"");
  double good = 0, all = 0;
  for(int bin=1; bin<=hVtxIndex->GetNbinsX(); bin++)
    {
      good += hVtxIndex->GetBinContent(bin,bin);
      for(int bin2=1; bin2<=hVtxIndex->GetNbinsY(); bin2++)
	{
	  all += hVtxIndex->GetBinContent(bin,bin2);
	}
    }
  TPaveText *t1 = GetPaveText(0.25,0.35,0.7,0.8,0.04);
  t1->AddText(Form("Diagonal: %4.2f%%",good/all*100));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_VtxIndex.pdf",run_type));

  // cloest
  TH1F *hVtxIndexClosest = (TH1F*)hVtxIndex->ProjectionX();
  c = draw1D(hVtxIndexClosest,"Index of TPC vertex closest to VPD",true,false);
  TPaveText *t1 = GetPaveText(0.5,0.55,0.7,0.8,0.04);
  t1->AddText(Form("Index=0: %4.2f%%",hVtxIndexClosest->GetBinContent(1)/hVtxIndexClosest->GetEntries()*100));
  t1->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/vtxSel_VtxIndexClosest.pdf",run_type));
}


//================================================
void vtxCuts(const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  THnSparseF *hn = (THnSparseF*)f->Get("mhEventMult_qa_di_mu");
  TH1F *hTpcVpdVz[4];
  hTpcVpdVz[0] = (TH1F*)hn->Projection(5);
  hTpcVpdVz[0]->SetName("hTpcVpdVz_All");

  // ranking > 0 cut
  hn->GetAxis(4)->SetRange(3,3);
  hTpcVpdVz[1] = (TH1F*)hn->Projection(5);
  hTpcVpdVz[1]->SetName("hTpcVpdVz_PosRanking");

  // |TpcVz| < 50 cm
  hn->GetAxis(3)->SetRangeUser(-50,50);
  hTpcVpdVz[2] = (TH1F*)hn->Projection(5);
  hTpcVpdVz[2]->SetName("hTpcVpdVz_TPcVz");

  // |TpcVz-VpdVz| < 5 cm
  hn->GetAxis(5)->SetRangeUser(-5,5);
  hTpcVpdVz[3] = (TH1F*)hn->Projection(5);
  hTpcVpdVz[3]->SetName("hTpcVpdVz_TPcVpdDz");

  TString legName[4] = {"All","Ranking > 0","|TpcVz| < 50 cm","|TpcVz-VpdVz| < 5 cm"};
  TList *list = new TList;
  for(int i=0; i<4; i++)
    {
      printf("[i] %s: %4.2f%%\n",legName[i].Data(),hTpcVpdVz[i]->GetEntries()/hTpcVpdVz[0]->GetEntries()*100);
      if(i<3) list->Add(hTpcVpdVz[i]);
    }
  c = drawHistos(list,"TpcVpdVzDiff",Form("Distribution of TpcVz-VpdVz;TpcVz-VpdVz (cm)"),kFALSE,0,100,kFALSE,-3,8,kTRUE,kTRUE,legName,kTRUE,"",0.6,0.75,0.6,0.8,kFALSE);
  double vzdiff_cut = 5;
  TLine *line = GetLine(-1*vzdiff_cut,0,-1*vzdiff_cut,hTpcVpdVz[0]->GetMaximum()*0.8);
  line->Draw();
  TLine *line = GetLine(vzdiff_cut,0,vzdiff_cut,hTpcVpdVz[0]->GetMaximum()*0.8);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/TpcVpdDz_VtxCut.pdf",run_type));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_vertex/TpcVpdDz_VtxCut.png",run_type));
    }
}


//================================================
void cuts(const Int_t save = 0)
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
