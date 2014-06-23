TFile *f;
const char *cut = "GlobalTrk.HLT";
Int_t hlt_index = 0;
Int_t trk_index = 0;

//================================================
void qa_track()
{
  gStyle->SetOptStat(0);

  TString cut_name = cut;
  if(cut_name.Contains("HLT"))
    hlt_index = 1;
  if(cut_name.Contains("Global"))
    trk_index = 1;

  f = TFile::Open(Form("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.%s.root",cut),"read");
  //trackDistribution();
  qualityCuts();
}

//================================================
void trackDistribution(const Int_t save = 0)
{
  TH1F *hStat = (TH1F*)f->Get("hEventStat");
  // eta vs pt
  TH2F *hTrkEta = (TH2F*)f->Get(Form("hTrkEta_%s",trigName[kTrigType]));
  c = draw2D(hTrkEta,Form("Au+Au %s: #eta vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkEtaVsPt_%s.png",cut,trigName[kTrigType]));
 
  // phi vs pt
  TH2F *hTrkPhi = (TH2F*)f->Get(Form("hTrkPhi_%s",trigName[kTrigType]));
  c = draw2D(hTrkPhi,Form("Au+Au %s: #phi vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkPhiVsPt_%s.png",cut,trigName[kTrigType]));

  // dE/dx vs pt
  TH2F *hTrkDedx = (TH2F*)f->Get(Form("hTrkDedx_%s",trigName[kTrigType]));
  c = draw2D(hTrkDedx,Form("Au+Au %s: dE/dx vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkDedxVsPt_%s.png",cut,trigName[kTrigType]));

  // pt distribution
  TH1F *hTrkPt = (TH1F*)f->Get(Form("hTrkPt_%s",trigName[kTrigType]));
  scaleHisto( hTrkPt, hStat->GetBinContent(kTrigType+7), 1, kTRUE);
  hTrkPt->SetMarkerStyle(21);
  c = draw1D(hTrkPt,Form("Au+Au %s: p_{T} distribution of %s tracks%s;p_{T} (GeV/c);1/N dN/dp_{T}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]),kTRUE,kTRUE);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkPt_%s.png",cut,trigName[kTrigType]));
}

//================================================
void qualityCuts(const Int_t save = 0)
{
  TH1F *hCuts = (TH1F*)f->Get("hAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;
  cout << "scale = " << scale << endl;

  // dca cut
  Double_t dcaCut = hCuts->GetBinContent(7)/scale;
  TH2F *hTrkDca = (TH2F*)f->Get(Form("hTrkDca_qa_%s",trigName[kTrigType]));
  hTrkDca->GetYaxis()->SetRangeUser(0,0.05);
  c = draw2D(hTrkDca,Form("Au+Au %s: dca vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  Double_t low_x = hTrkDca->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hTrkDca->GetXaxis()->GetBinUpEdge(hTrkDca->GetNbinsX());
  TLine *line = GetLine(low_x,dcaCut,up_x,dcaCut,1,3);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.QualityCut_TrkDca_%s.png",cut,trigName[kTrigType]));

  // # of fit hits cut
  Double_t nHitCut = hCuts->GetBinContent(5)/scale;
  TH2F *hTrkNhitsFit = (TH2F*)f->Get(Form("hTrkNhitsFit_qa_%s",trigName[kTrigType]));
  c = draw2D(hTrkNhitsFit,Form("Au+Au %s: NHitsFit vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  Double_t low_x = hTrkNhitsFit->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hTrkNhitsFit->GetXaxis()->GetBinUpEdge(hTrkNhitsFit->GetNbinsX());
  TLine *line = GetLine(low_x,nHitCut,up_x,nHitCut,1,3);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.QualityCut_TrkNhitsFit_%s.png",cut,trigName[kTrigType]));

  // # of de/dx hits cut
  Double_t nDedxHitCut = hCuts->GetBinContent(6)/scale;
  TH2F *hTrkNhitsDedx  = (TH2F*)f->Get(Form("hTrkNhitsDedx_qa_%s",trigName[kTrigType]));
  c = draw2D(hTrkNhitsDedx,Form("Au+Au %s: NHitsDedx vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  Double_t low_x = hTrkNhitsDedx->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hTrkNhitsDedx->GetXaxis()->GetBinUpEdge(hTrkNhitsDedx->GetNbinsX());
  TLine *line = GetLine(low_x,nDedxHitCut,up_x,nDedxHitCut,1,3);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.QualityCut_TrkNhitsDedx_%s.png",cut,trigName[kTrigType]));

  // fraction of fit hits cut
  Double_t minFracCut = hCuts->GetBinContent(10)/scale;
  TH2F *hTrkFitHitFrac = (TH2F*)f->Get(Form("hTrkFitHitFrac_qa_%s",trigName[kTrigType]));
  hTrkFitHitFrac->SetYTitle("fraction");
  c = draw2D(hTrkFitHitFrac,Form("Au+Au %s: hit fraction vs p_{T} of %s tracks%s",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  Double_t low_x = hTrkFitHitFrac->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hTrkFitHitFrac->GetXaxis()->GetBinUpEdge(hTrkFitHitFrac->GetNbinsX());
  TLine *line = GetLine(low_x,minFracCut,up_x,minFracCut,1,3);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.QualityCut_TrkFitHitFrac_%s.png",cut,trigName[kTrigType]));

  // nSigmaPi cut
  TH2F *hTrkNsigmaPi = (TH2F*)f->Get(Form("hTrkNsigmaPi_qa_%s",trigName[kTrigType]));
  Double_t minNsigmaCut = hCuts->GetBinContent(8)/scale;
  Double_t maxNsigmaCut = hCuts->GetBinContent(9)/scale;
  c = draw2D(hTrkNsigmaPi,Form("Au+Au %s: n#sigma_{#pi} vs p_{T} of %s tracks%s;p_{T} (GeV/c);n#sigma_{#pi}",trigName[kTrigType],trk_name[trk_index],hlt_name[hlt_index]));
  Double_t low_x = hTrkNsigmaPi->GetXaxis()->GetBinLowEdge(1);
  Double_t up_x  = hTrkNsigmaPi->GetXaxis()->GetBinUpEdge(hTrkNsigmaPi->GetNbinsX());
  TLine *line = GetLine(low_x,minNsigmaCut,up_x,minNsigmaCut,1,3);
  line->Draw();
  line = GetLine(low_x,maxNsigmaCut,up_x,maxNsigmaCut,1,3);
  line->Draw();
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.QualityCut_TrkNsigmaPi_%s.png",cut,trigName[kTrigType]));

  // Track pt for various cuts
  TH2F *hTrkPt_cuts = (TH2F*)f->Get(Form("hTrkPt_cuts_qa_%s",trigName[kTrigType]));
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
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkPt_cuts_%s.png",cut,trigName[kTrigType]));


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
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkFraction_cuts_%s.png",cut,trigName[kTrigType]));
  if(save) c1->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/%s.TrkFraction_cuts_semilog_%s.png",cut,trigName[kTrigType]));
  
}
