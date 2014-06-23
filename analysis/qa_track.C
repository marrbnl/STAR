TFile *f;
//================================================
void qa_track()
{
  gStyle->SetOptStat(0);

  f = TFile::Open("~/Work/STAR/analysis/Output/jpsi.AuAu200.Run14.DzCut.root","read");
  //trackDistribution();
  qualityCuts();
}

//================================================
void trackDistribution(const Int_t save = 1)
{
  TH1F *hStat = (TH1F*)f->Get("mhEventStat");
  // eta vs pt
  TH2F *hPrimTrkEta[3];
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkEta[i] = (TH2F*)f->Get(Form("hPrimTrkEta_%s",trigName[i]));
      c = draw2D(hPrimTrkEta[i]);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/PrimTrkEtaVsPt_%s.png",trigName[i]));
    }

  // phi vs pt
  TH2F *hPrimTrkPhi[3];
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkPhi[i] = (TH2F*)f->Get(Form("hPrimTrkPhi_%s",trigName[i]));
      c = draw2D(hPrimTrkPhi[i]);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/PrimTrkPhiVsPt_%s.png",trigName[i]));
    }

  // dE/dx vs pt
  TH2F *hPrimTrkDedx[3];
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkDedx[i] = (TH2F*)f->Get(Form("hPrimTrkDedx_%s",trigName[i]));
      c = draw2D(hPrimTrkDedx[i]);
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/PrimTrkDedxVsP_%s.png",trigName[i]));
    }

  // pt distribution
  TH1F *hPrimTrkPt[3];
  TList *list = new TList;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkPt[i] = (TH1F*)f->Get(Form("hPrimTrkPt_%s",trigName[i]));
      scaleHisto( hPrimTrkPt[i], hStat->GetBinContent(i+4), 1, kTRUE);
      list->Add(hPrimTrkPt[i]);
    }
  c = drawHistos(list,"hPrimTrkPt","Au+Au: p_{T} distribution of primary tracks;p_{T} (GeV/c);1/N dN/dp_{T}",kFALSE,0,pi/2,kFALSE,-0.005,0.02,kTRUE,kTRUE,legName,kTRUE,"",0.5,0.68,0.6,0.85);
  if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/PrimTrkPt.png"));
}

//================================================
void qualityCuts(const Int_t save = 1)
{
  TH1F *hCuts = (TH1F*)f->Get("mhAnalysisCuts");
  Double_t scale = hCuts->GetBinContent(3)/1e4;
  cout << "scale = " << scale << endl;

  // dca cut
  TH2F *hPrimTrkDca[3];
  Double_t dcaCut = hCuts->GetBinContent(7)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkDca[i] = (TH2F*)f->Get(Form("hPrimTrkDca_qa_%s",trigName[i]));
      hPrimTrkDca[i]->GetYaxis()->SetRangeUser(0,0.05);
      c = draw2D(hPrimTrkDca[i]);
      Double_t low_x = hPrimTrkDca[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hPrimTrkDca[i]->GetXaxis()->GetBinUpEdge(hPrimTrkDca[i]->GetNbinsX());
      TLine *line = GetLine(low_x,dcaCut,up_x,dcaCut,1,3);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/QualityCut_PrimTrkDca_%s.png",trigName[i]));
    }

  // # of fit hits cut
  TH2F *hPrimTrkNhitsFit[3];
  Double_t nHitCut = hCuts->GetBinContent(5)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkNhitsFit[i] = (TH2F*)f->Get(Form("hPrimTrkNhitsFit_qa_%s",trigName[i]));
      c = draw2D(hPrimTrkNhitsFit[i]);
      Double_t low_x = hPrimTrkNhitsFit[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hPrimTrkNhitsFit[i]->GetXaxis()->GetBinUpEdge(hPrimTrkNhitsFit[i]->GetNbinsX());
      TLine *line = GetLine(low_x,nHitCut,up_x,nHitCut,1,3);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/QualityCut_PrimTrkNhitsFit_%s.png",trigName[i]));
    }

  // # of de/dx hits cut
  TH2F *hPrimTrkNhitsDedx[3];
  Double_t nDedxHitCut = hCuts->GetBinContent(6)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkNhitsDedx[i] = (TH2F*)f->Get(Form("hPrimTrkNhitsDedx_qa_%s",trigName[i]));
      c = draw2D(hPrimTrkNhitsDedx[i]);
      Double_t low_x = hPrimTrkNhitsDedx[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hPrimTrkNhitsDedx[i]->GetXaxis()->GetBinUpEdge(hPrimTrkNhitsDedx[i]->GetNbinsX());
      TLine *line = GetLine(low_x,nDedxHitCut,up_x,nDedxHitCut,1,3);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/QualityCut_PrimTrkNhitsDedx_%s.png",trigName[i]));
    }

  // fraction of fit hits cut
  TH2F *hPrimTrkFitHitFrac[3];
  Double_t minFracCut = hCuts->GetBinContent(10)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkFitHitFrac[i] = (TH2F*)f->Get(Form("hPrimTrkFitHitFrac_qa_%s",trigName[i]));
      hPrimTrkFitHitFrac[i]->SetYTitle("fraction");
      c = draw2D(hPrimTrkFitHitFrac[i]);
      Double_t low_x = hPrimTrkFitHitFrac[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hPrimTrkFitHitFrac[i]->GetXaxis()->GetBinUpEdge(hPrimTrkFitHitFrac[i]->GetNbinsX());
      TLine *line = GetLine(low_x,minFracCut,up_x,minFracCut,1,3);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/QualityCut_PrimTrkFitHitFrac_%s.png",trigName[i]));
    }

  // nSigmaPi cut
  TH2F *hPrimTrkNsigmaPi[3];
  Double_t minNsigmaCut = hCuts->GetBinContent(8)/scale;
  Double_t maxNsigmaCut = hCuts->GetBinContent(9)/scale;
  for(Int_t i=0; i<3; i++)
    {
      hPrimTrkNsigmaPi[i] = (TH2F*)f->Get(Form("hPrimTrkNsigmaPi_qa_%s",trigName[i]));
      c = draw2D(hPrimTrkNsigmaPi[i],Form("%s: n#sigma_{#pi} of primary tracks;p_{T} (GeV/c);n#sigma_{#pi}",trigName[i]));
      Double_t low_x = hPrimTrkNsigmaPi[i]->GetXaxis()->GetBinLowEdge(1);
      Double_t up_x  = hPrimTrkNsigmaPi[i]->GetXaxis()->GetBinUpEdge(hPrimTrkNsigmaPi[i]->GetNbinsX());
      TLine *line = GetLine(low_x,minNsigmaCut,up_x,minNsigmaCut,1,3);
      line->Draw();
      line = GetLine(low_x,maxNsigmaCut,up_x,maxNsigmaCut,1,3);
      line->Draw();
      if(save) c->SaveAs(Form("~/Work/STAR/analysis/Plots/qa_track/QualityCut_PrimTrkNsigmaPi_%s.png",trigName[i]));
    }
}
