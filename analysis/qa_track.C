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
  f = TFile::Open("./output/Run13.pp500.jpsi.EventQA.root","read");

  //qa();
  qualityCuts();
}

//================================================
void qa(const Int_t save = 0)
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
      if(save) 
	{
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.pdf",run_type,outname.Data(),trigName[kTrigType]));
	  c->SaveAs(Form("~/Work/STAR/analysis/Plots/%s/qa_track/%s_%s.png",run_type,outname.Data(),trigName[kTrigType]));
	}
    }
}

//================================================
void qualityCuts(const Int_t save = 0)
{
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
  char *cuts[nCuts] = {"All tracks","p_{T} > 1 GeV/c", "|#eta| < 0.8","0 < #varphi < 2#pi",Form("NHitsFit > %1.0f",nHitCut), Form("NHitsDedx > %1.0f",nDedxHitCut), Form("dca > %1.1f cm",dcaCut), Form("%1.1f < n#sigma_{#pi} < %1.1f",minNsigmaCut,maxNsigmaCut),Form("HitFraction > %0.2f",minFracCut),"Matched to MTD hits","Matched to good hits","|#Deltaz| < 20 cm"};
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
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkPt_cuts_%s.png",run_config,trigName[kTrigType]));


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
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkFraction_cuts_%s.png",run_config,trigName[kTrigType]));
  if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkFraction_cuts_semilog_%s.png",run_config,trigName[kTrigType]));
  
}
