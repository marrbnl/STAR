TFile *f;
Int_t hlt_index = 0;
Int_t trk_index = 0;
const char *run_config = "PrimTrk.ClosePrimVtx.40cm";

//================================================
void qa_track()
{
  gStyle->SetOptStat(1);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);  

  TString cut_name = run_config;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  //f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",run_config),"read");
  //f = TFile::Open("./output/Run13.pp500.jpsi.EventQA.root","read");
  f = TFile::Open("./output/Run13.pp500.jpsi.PID.root","read");

  //qa();
  //qualityCuts();
  //distribution();
  cutCorrelation();
}

//================================================
void distribution(const Int_t save = 0)
{
  gStyle->SetOptStat(0);
  const char *hName[] = {"mhTrkN","mhTrkPt","mhTrkEtaPhi","mhTrkDedx",
			 "mhMthTrkN","mhMthTrkPt","mhMthTrkEtaPhi"};

  for(Int_t i=0; i<7; i++)
    {
      TH1 *h = (TH1*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      if(!h->IsA()->InheritsFrom("TH2")) 
	{
	  TPaveText *t1 = 0;
	  if(hName[i]=="mhTrkN" || hName[i]=="mhMthTrkN")
	    {
	      if(run_type == "Run13_pp500") h->GetXaxis()->SetRangeUser(0,20);
	      t1 = GetPaveText(0.6,0.7,0.7,0.8);
	      t1->AddText(Form("<#> = %2.2f",h->GetMean()));
	    }
	  else if(hName[i]=="mhTrkPt" || hName[i]=="mhMthTrkPt")
	    {
	      h->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");
	      scaleHisto(h,1,1,kTRUE);
	      t1 = GetPaveText(0.6,0.7,0.7,0.8);
	      t1->AddText(Form("<p_{T}> = %2.2f GeV/c",h->GetMean()));
	    }
	  c = draw1D(h,"",kTRUE,kFALSE);
	  if(t1) t1->Draw();
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
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
}

//================================================
void qa(const Int_t save = 1)
{
  const char *hName[] = {"mhTrkEta_qa","mhTrkPhi_qa","mhTrkDca_qa","mhTrkDcaxy_qa",
			 "mhTrkDcaz_qa","mhTrkNhitsFit_qa","mhTrkNhitsDedx_qa","mhTrkFitHitFrac_qa"};

  for(Int_t i=0; i<8; i++)
    {
      TH1 *h = (TH1*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      if(h->IsA()->InheritsFrom("TH2")) 
	{
	  TH2F *h2 = (TH2F*)h;
	  c = draw2D(h2,"");
	}

      TString outname = hName[i];
      outname.ReplaceAll("mh","");
      outname.ReplaceAll("_qa","");
      outname = "qa_" + outname;
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
}

//================================================
void qualityCuts(const Int_t save = 1)
{
  gStyle->SetOptStat(0);
  TH1F *hCuts = (TH1F*)f->Get("hAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;
  cout << "scale = " << scale << endl;

  const char *hName[] = {"mhTrkDca_qa","mhTrkNhitsFit_qa","mhTrkNhitsDedx_qa","mhTrkFitHitFrac_qa"};
  const Int_t cut_index[] = {7,5,6,10};

  for(Int_t i=0; i<4; i++)
    {
      TH2 *h2 = (TH2*)f->Get(Form("%s_%s",hName[i],trigName[kTrigType]));
      TH1 *h1 = (TH1*)h2->ProjectionY(Form("%s_pro_%s",hName[i],trigName[kTrigType]));
      Bool_t logy = kFALSE;
      if(hName[i]=="mhTrkDca_qa") logy = kTRUE;
      c = draw1D(h1,"",logy,kFALSE);

      Double_t cut = hCuts->GetBinContent(cut_index[i])/scale;
      Double_t low_x = h1->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());
      TLine *line = GetLine(cut,0,cut,0.8*h1->GetMaximum(),2,3);
      line->Draw();

      TString outname = hName[i];
      outname.ReplaceAll("mh","");
      outname.ReplaceAll("_qa","");
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/Cut%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/Cut%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
  

  // Track pt for various cuts
  TH2F *hTrkPt_cuts = (TH2F*)f->Get(Form("mhTrkPt_cuts_qa_%s",trigName[kTrigType]));
  draw2D(hTrkPt_cuts);
  const Int_t nCuts = 12;
  TH1F *hTrkPt[nCuts];
  TH1F *hTrkRatio[nCuts];
  const char* option[2] = {"No","Yes"};
  char *cuts[nCuts] = {"All tracks","p_{T} > 1 GeV/c", "|#eta| < 0.8","0 < #varphi < 2#pi",Form("NHitsFit > %1.0f",hCuts->GetBinContent(5)/scale), Form("NHitsDedx > %1.0f",hCuts->GetBinContent(6)/scale), Form("dca < %1.1f cm",hCuts->GetBinContent(7)/scale), Form("%1.1f < n#sigma_{#pi} < %1.1f",hCuts->GetBinContent(8)/scale,hCuts->GetBinContent(9)/scale),Form("HitFraction > %0.2f",hCuts->GetBinContent(10)/scale),Form("Require to match TOF: %s",option[(Int_t)(hCuts->GetBinContent(15)/scale)]),"Matched to MTD hits",Form("|#Deltaz| < %1.0f cm",hCuts->GetBinContent(14)/scale)};
  if(trk_index==1) cuts[6] = "No DCA cut";

  TLegend *leg = new TLegend(0.3,0.63,0.45,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);

  TLegend *leg1 = new TLegend(0.55,0.63,0.7,0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);

  for(Int_t i=0; i<nCuts; i++)
    {
      hTrkPt[i] = (TH1F*)hTrkPt_cuts->ProjectionY(Form("TrkPt_Cut%d",i),i+1,i+1);
      scaleHisto( hTrkPt[i], 1, 1, kTRUE);
      hTrkPt[i]->SetMarkerStyle(21);
      hTrkPt[i]->SetMarkerColor(color[nCuts-1-i]);
      hTrkPt[i]->SetLineColor(color[nCuts-i-1]);
      hTrkPt[i]->GetYaxis()->SetRangeUser(1,1e12);
      if(i==0)
	c = draw1D(hTrkPt[i],Form("Au+Au %s: p_{T} distribution of %s tracks%s;p_{T} (GeV/c)",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kTRUE);
      else
	hTrkPt[i]->Draw("sames");
      if(i<nCuts/2) leg->AddEntry(hTrkPt[i],cuts[i],"PL");
      else          leg1->AddEntry(hTrkPt[i],cuts[i],"PL");
    }
  leg->Draw();
  leg1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkPt_cuts_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkPt_cuts_%s.png",run_type,trigName[kTrigType]));
    }


  for(Int_t i=0; i<nCuts; i++)
    {
      hTrkRatio[i] = (TH1F*)hTrkPt[i]->Clone(Form("hTrkRatio_Cut%d",i));
      if(i==0) hTrkRatio[i]->Divide(hTrkPt[i]);
      else     hTrkRatio[i]->Divide(hTrkPt[i-1]);
      hTrkRatio[i]->GetYaxis()->SetRangeUser(1e-3,2);

      if(i==0)
	{
	  c  = draw1D(hTrkRatio[i],Form("Au+Au %s: fraction of %s tracks surviving a cut%s;p_{T} (GeV/c)",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kFALSE);
	  hTrkRatio[i]->SetName(Form("%s_log",hTrkRatio[i]->GetName()));
	  hTrkRatio[i]->GetYaxis()->SetRangeUser(1e-3,1e2);
	  c1 = draw1D(hTrkRatio[i],Form("Au+Au %s: fraction of %s tracks surviving a cut%s;p_{T} (GeV/c)",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kTRUE);
	}
      else
	{
	  c->cd();
	  hTrkRatio[i]->Draw("sames");
	  c1->cd();
	  hTrkRatio[i]->Draw("sames");
	}
    }
  c->cd();
  leg->Draw();
  leg1->Draw();
  c1->cd();
  leg->Draw();
  leg1->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkFraction_cuts_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkFraction_cuts_%s.png",run_type,trigName[kTrigType]));

      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkFraction_cuts_semilog_%s.pdf",run_type,trigName[kTrigType]));
      c1->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/TrkFraction_cuts_semilog_%s.png",run_type,trigName[kTrigType]));
    }
  
}


//================================================
void cutCorrelation(const Int_t save = 0)
{
  gStyle->SetOptStat(0);
  THnSparseF *hTrkInfo = (THnSparseF*)f->Get(Form("mhTrkInfo_qa_%s",trigName[kTrigType]));
  const Double_t pt_cut = 1;
  hTrkInfo->GetAxis(0)->SetRangeUser(pt_cut,100);
  TH2F *hNHitVsHitFrac = (TH2F*)hTrkInfo->Projection(2,1);
  hNHitVsHitFrac->SetName("hNHitVsHitFrac_all");
  hNHitVsHitFrac->SetTitle(Form("%s: FitHitFrac vs NHitsFit of primary tracks (p_{T} > %1.0f GeV/c);NHitsFit;NHitsFit/NHitsPoss",trigName[kTrigType],pt_cut));
  c = draw2D(hNHitVsHitFrac);
  TLine *line = GetLine(15,0,15,1,1);
  line->Draw();
  line = GetLine(0,0.52,45,0.52,1);
  line->Draw();
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/qa_NHitsFit_vs_NHitsFrac_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/qa_NHitsFit_vs_NHitsFrac_%s.png",run_type,trigName[kTrigType]));
    }

  TList *list = new TList;
  TH1F *hNHitsFit[3];
  for(Int_t i=0; i<3; i++)
    {
      if(i==1) hTrkInfo->GetAxis(2)->SetRangeUser(0,0.52-0.01);
      else if (i==2) hTrkInfo->GetAxis(2)->SetRangeUser(0.52+0.01,1);

      hNHitsFit[i] = (TH1F*)hTrkInfo->Projection(1);
      hNHitsFit[i]->SetName(Form("hNHitsFit_%d",i));
      hNHitsFit[i]->SetMaximum(1.2*hNHitsFit[i]->GetMaximum());
      list->Add(hNHitsFit[i]);
      hTrkInfo->GetAxis(2)->SetRange(0,-1);
    }
  TString legName[3] = {"All","NHitsFit/NHitsPoss < 0.52", "NHitsFit/NHitsPoss > 0.52"};
  c = drawHistos(list,"NHitsFit",Form("%s: # of TPC hits for fitting of primary tracks;NHitsFit",trigName[kTrigType]),kFALSE,0,100,kFALSE,1e-6,0.018,kFALSE,kTRUE,legName,kTRUE,Form("p_{T} > %1.0f GeV/c",pt_cut),0.15,0.25,0.62,0.85,kFALSE);
  if(save) 
    {
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/qa_NHitsFit_CutNHitsFrac_%s.pdf",run_type,trigName[kTrigType]));
      c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/qa_NHitsFit_CutNHitsFrac_%s.png",run_type,trigName[kTrigType]));
    }
}
